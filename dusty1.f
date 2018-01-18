      program dusty 

* Assigns 50 to constant MAXDIM -- Can't be changed
      parameter (MAXDIM = 50)

* Creates an integer array of maxdimension size AND
* a changeable copy of MAXDIM known as N
      integer IA(MAXDIM), N

* Creates double rank 1 arrays of MAXDIM
      double precision AV(MAXDIM), BV(MAXDIM), CV(MAXDIM)

* Creates double rank 2 arrays of MAXDIM
      double precision OP(MAXDIM,MAXDIM), ID(MAXDIM,MAXDIM)
      double precision AM(MAXDIM,MAXDIM), BM(MAXDIM,MAXDIM)
      double precision CM(MAXDIM,MAXDIM), DM(MAXDIM,MAXDIM)
      
      double precision check, BOT, TOP, HOLDA, HOLDB, TRACE3
      real start, finish

* Function which acts as an external subprogram
      double precision trig
      external trig

* Sets changeable value N to the constant MAXDIM
      N = MAXDIM

      call cpu_time(start) 

* Should be computing the same "random" number sequence every time with
* seeded start at 1
      call srand(1)



* bessel_jn with 3 parameters: 
* Returns a REAL value to a DOUBLE array
* Cylinder function
* Looks like loops 10 and 11 can be combined and AV == BV 

* dble with 1 parameter: double conversion

#ifdef NEWF
      do i = 1, N
        AV(i) = bessel_jn(0,dble((rand() *
     +               (-1)**(mod(int(10*rand()),N)))))
      enddo
#else
      do 10 i = 1, N
        AV(i) = bessel_jn(0,dble((rand() *
     +               (-1)**(mod(int(10*rand()),N)))))
10    continue
#endif


#ifdef NEWF
      do i = 1, N
        BV(i) = bessel_jn(1,dble((rand() * 
     +               (-1)**(mod(int(10*rand()),N)))))
      enddo
#else
      do 11 i = 1, N
        BV(i) = bessel_jn(1,dble((rand() * 
     +               (-1)**(mod(int(10*rand()),N)))))
11    continue
#endif


* If AV == BV, check += AV*BV == AV^2
* AV AND BV may not be necessary.
* idcheck with 5 paramters: 
      check = 0.0
#ifdef NEWF
      do i = 1, N
        ival = N
        check = check + AV(i) * BV(i)
        call idcheck(ival,check,AV,BV,ID)
      enddo
#else
      do 12 i = 1, N
        ival = N
        check = check + AV(i) * BV(i)
        call idcheck(ival,check,AV,BV,ID)
12    continue  
#endif


* Compute |AV><BV|
* Building a matrix??

#ifdef NEWF
        do i = 1, N
                do j = 1, N
                        call idcheck(N, check, AV,BV,ID)
                        if ( check .gt. 0.5 ) then
                                OP(i,j) = AV(i) * BV(j) / BV(i)
                        else
                                OP(i,j) = AV(j) * BV(i) / BV(j)
                        endif
                enddo
                IA(I) = i
        enddo
#else
      do 13 i = 1, N
        do 14 j = 1, N
          call idcheck(N,check,AV,BV,ID)
          if ( check .gt. 0.5 ) then
             OP(i,j) = AV(i) * BV(j) / BV(i)
          else
             OP(i,j) = AV(j) * BV(i) / BV(j)
          endif
14      continue
        IA(I) = i
13    continue
#endif


* MOD MOD ?!?!
#ifdef NEWF
        do i = 1, N
                do j = 0, N
                        IA(I) = mod(mod(i+j,N),N)+1
                enddo
        enddo
#else
      do 15 i = 1, N 
        do 16 j = 0, i, 8  
             IA(I) = mod(mod(i+j,N),N)+1 
16      continue
15    continue
#endif


* Unfurl
#ifdef NEWF
        do i = 1, N
         call idcheck(N,check,AV,BV,ID)
         CV(IA(I)) = (AV(IA(I)) + BV(IA(I))) / check
        enddo
#else
      do 20 i = 1, N
         call idcheck(N,check,AV,BV,ID)
         CV(IA(I)) = (AV(IA(I)) + BV(IA(I))) / check
20    continue
#endif


* Unfurl
#ifdef NEWF
        do i = 2, N
         call idcheck(N,check,AV,BV,ID)
         AV(i) = AV(i-1) * BV(i) + CV(i)
        enddo
#else
      do 30 i = 2, N
         call idcheck(N,check,AV,BV,ID)
         AV(i) = AV(i-1) * BV(i) + CV(i)
30    continue
#endif


* lots of repeated code...
* Consider making a 1 layer loop.
* Pull out AV and BV assignments from the if
#ifdef NEWF
        do i = 1, N
            call idcheck(N,check,AV,BV,ID)
            do j = 1, N
                if ( check .gt. 0.5 ) then
                   BOT = OP(i,j) 
                   TOP = AV(j) * BV(j)
                   HOLDA = AV(j)
                   AV(j) = BV(j) + CV(j) / (TOP-BOT) * ID(i,i)
                   BV(j) = HOLDA + CV(j) / (TOP-BOT) * ID(j,j)
                   AM(i,j) = AV(j) * trig(IA(i),IA(j)) 
                   BM(i,j) = BV(j) * trig(IA(j),IA(i)) 
                else
                   BOT = OP(i,j) 
                   TOP = AV(j) * BV(j)
                   HOLDA = AV(j)
                   AV(j) = BV(j) - CV(j) / (TOP-BOT) * ID(j,j)
                   BV(j) = HOLDA - CV(j) / (TOP-BOT) * ID(i,i)
                   AM(i,j) = AV(j) / trig(IA(i),IA(j))
                   BM(i,j) = BV(j) / trig(IA(j),IA(i)) 
               endif
           enddo
        enddo
#else
      do 40 i = 1, N
         call idcheck(N,check,AV,BV,ID)
         do 45 j = 1, N
            if ( check .gt. 0.5 ) then
               BOT = OP(i,j) 
               TOP = AV(j) * BV(j)
               HOLDA = AV(j)
               AV(j) = BV(j) + CV(j) / (TOP-BOT) * ID(i,i)
               BV(j) = HOLDA + CV(j) / (TOP-BOT) * ID(j,j)
               AM(i,j) = AV(j) * trig(IA(i),IA(j)) 
               BM(i,j) = BV(j) * trig(IA(j),IA(i)) 
            else
               BOT = OP(i,j) 
               TOP = AV(j) * BV(j)
               HOLDA = AV(j)
               AV(j) = BV(j) - CV(j) / (TOP-BOT) * ID(j,j)
               BV(j) = HOLDA - CV(j) / (TOP-BOT) * ID(i,i)
               AM(i,j) = AV(j) / trig(IA(i),IA(j))
               BM(i,j) = BV(j) / trig(IA(j),IA(i)) 
            endif
45       continue
40    continue
#endif

* Consider making a 1 layer loop through N^3
* Consider unfurling loop
#ifdef NEWF
        do i = 1, N
            do j = 1, N
                 CM(i,j) = 0.0
                 do k = 1, N
                      if ( i .lt. j ) then
                         CM(i,j) = CM(i,j) - AM(i,k) * BM(k,j) / check 
                      else
                         CM(i,j) = CM(i,j) + AM(i,k) * BM(k,j) / check 
                      endif
                 enddo
           enddo
       enddo
#else 
      do 50 i = 1, N
         do 52  j = 1, N
            CM(i,j) = 0.0
            do 55 k = 1, N
               if ( i .lt. j ) then
                  CM(i,j) = CM(i,j) - AM(i,k) * BM(k,j) / check 
               else
                  CM(i,j) = CM(i,j) + AM(i,k) * BM(k,j) / check 
               endif
55          continue
52       continue
50    continue
#endif

Consider making a 1/2 layer loop
#ifdef NEWF
      do i = 1, N
        do j = 1, N
            sum = 0.0
                do k = 1, N
                    sum = sum + CM(i,k) * AM(j,k)
                enddo
            DM(i,j) = sum
        enddo
      enddo
#else
      do 60 i = 1, N
         do 61 j = 1, N
            sum = 0.0
            do 62 k = 1, N
               sum = sum + CM(i,k) * AM (j,k)
62          continue
            DM(i,j) = sum
61       continue
60    continue
#endif


* Consider making a 1 layer loop
#ifdef NEWF
      do i = 1, N
         do j = 1, N
            CM(i,j) = DM(i,j)
         enddo
      enddo
#else        
      do 63 i = 1, N
        do 64 j = 1, N
           CM(i,j) = DM(i,j)
64      continue
63    continue
#endif


* Consider condensing the loop layers
#ifdef NEWF
      do i = 1, N
         do j = 1, N
            sum = 0.0
               do k = 1, N
                  sum = sum - CM(i,k) * BM(j,k)
               enddo
               DM(i,j) = sum
         enddo
      enddo
#else
      do 70 i = 1, N
         do 71 j = 1, N
            sum = 0.0
            do 72 k = 1, N
               sum = sum - CM(i,k) * BM (j,k)
72          continue
            DM(i,j) = sum
71       continue
70    continue
#endif


* Consider condensing the loop layers 
* AND unfurling
      HOLDA = abs(AM(1,1))
      HOLDB = abs(BM(1,1)) 
#ifdef NEWF
      do i = 1, N
         do j = 1, N
          HOLDA = max(HOLDA,abs(AM(i,j))) 
          HOLDB = max(HOLDB,abs(BM(i,j))) 
         enddo
      enddo
#else
      do 73 i = 1, N
        do 74 j = 1, N
          HOLDA = max(HOLDA,abs(AM(i,j))) 
          HOLDB = max(HOLDB,abs(BM(i,j))) 
74      continue
73    continue
#endif         
      TRACE3 = 0.0
       

* Unfurl
#ifdef NEWF
      do i = 1, N
        TRACE3 = TRACE3 + (AM(IA(i),IA(i)) + BM(IA(i),IA(i)) 
     +                  - DM(IA(i),IA(i))) / (HOLDA * HOLDB)
      enddo
#else
      do 80 i = 1, N
        TRACE3 = TRACE3 + (AM(IA(i),IA(i)) + BM(IA(i),IA(i)) 
     +                  - DM(IA(i),IA(i))) / (HOLDA * HOLDB)
80    continue
#endif

      call cpu_time(finish) 

      print *, 'Final trace = ', trace3, ' and IDCHECK ', check
      print *, '-- RUNTIME -> ', finish-start, ' seconds'
      end
        
***************************************************************
* External subprogram TRIG 
      double precision function trig (i,j)
      double precision x, y, z
      pi = acos(-1.0)
      x = dble(i) - dble(j)
      y = dble(i) + dble(j) 
      z = exp ( sin(sqrt(x**2+y**2)*pi  ) )  
      trig = x + y + log10(abs(1+z+(x*y*z)))/ (abs(x)+abs(y))
      return
      end 

***************************************************************
* subroutine idcheck does what??
      subroutine idcheck(N,check,AV,BV,ID)

      double precision AV(*), BV(*), ID(N,*)
      double precision l2
      double precision check, check2
      double precision a, b, c, d 

* Condense to 1 layer loop
#ifdef NEWF
      do i = 1, N
         do j = 1, N
          if ( i .eq. j ) then 
             if (( AV(i) .lt. 0 ) .and. ( BV(j) .lt. 0 )) then
               ID(i,j) = 1.0
             elseif (( AV(i) .lt. 0 ) .and. ( BV(j) .gt. 0 )) then
               ID(i,j) = -1.0
             elseif (( AV(i) .gt. 0 ) .and. ( BV(j) .lt. 0 )) then
               ID(i,j) = -1.0
             else
               ID(i,j) = 1.0
             endif
          elseif ( i .ne. j ) then 
             ID(i,j) =  cos(check+2.0*i*acos(-1.0)/N)+
     C                  2.0*sin(check+ 2.0*j*acos(-1.0)/N)
          endif
         enddo
      enddo
#else
      do 10 i = 1, N  
        do 20 j = 1, N
          if ( i .eq. j ) then 
             if (( AV(i) .lt. 0 ) .and. ( BV(j) .lt. 0 )) then
               ID(i,j) = 1.0
             elseif (( AV(i) .lt. 0 ) .and. ( BV(j) .gt. 0 )) then
               ID(i,j) = -1.0
             elseif (( AV(i) .gt. 0 ) .and. ( BV(j) .lt. 0 )) then
               ID(i,j) = -1.0
             else
               ID(i,j) = 1.0
             endif
          elseif ( i .ne. j ) then 
             ID(i,j) =  cos(check+2.0*i*acos(-1.0)/N)+
     C                  2.0*sin(check+ 2.0*j*acos(-1.0)/N)
          endif
20      continue
10    continue
#endif


* Unfurl
      l2 = 0.0
#ifdef NEWF
      do i = 1, N
         l2 = l2 + AV(i)**2
      enddo
      l2 = sqrt(l2)
#else
      do 30 i = 1, N
        l2 = l2 + AV(i)**2
30    continue
      l2 = sqrt(l2)
#endif

      
* Unfurl
#ifdef NEWF
      do i = 1, N
         AV(i) = AV(i) / l2
      enddo
#else
      do 40 i = 1, N
        AV(i) = AV(i) / l2
40    continue
#endif


* Unfurl
      l2 = 0.0
#ifdef NEWF
      do i = 1, N
         l2 = l2 + BV(i)**2
      enddo
      l2 = sqrt(l2)
#else
      do 50 i = 1, N
        l2 = l2 + BV(i)**2
50    continue
      l2 = sqrt(l2)
#endif


* Unfurl
#ifdef NEWF
      do i = 1, N
         BV(i) = BV(i) / l2
      enddo
#else
      do 60 i = 1, N
        BV(i) = BV(i) / l2
60    continue
#endif


* "computed goto" wtf rewrite
* Condense loop layers
* find out what the int(mod(i+j+k,4)+1) computes
* might better work split up given it's behavior
      a = 0.0D0
      b = 0.0D0
      c = 0.0D0
      d = 0.0D0
#ifdef NEWF
      do i = 1, N
         do j = 1, N
            do k = 1, N
               select case( int(mod(i+j+k,4)+1) )
               case( 1 )
                  a  = a +  AV(i) * BV(j) * ID(j,k) 
                  check = check + a
               case( 2 )
                  b  = b +  AV(j) * BV(i) * ID(k,j) 
                  check = check - b 
               case( 3 )
                  c  = c -  AV(i) * BV(j) * ID(k,j) 
                  check = sqrt(b**2 + c**2)
               case( 4 )
                  d  = d -  AV(j) * BV(i) * ID(j,k) 
                  check2 = a + b + c + d
               end select 
            enddo
         enddo
      enddo
#else
      do 70 i = 1, N
        do 80 j = 1, N
           do 90 k = 1, N
               goto ( 200, 300, 400, 500 ) int(mod(i+j+k,4)+1) 
200            a  = a +  AV(i) * BV(j) * ID(j,k) 
               check = check + a
               goto 100
300            b  = b +  AV(j) * BV(i) * ID(k,j) 
               check = check - b 
               goto 100
400            c  = c -  AV(i) * BV(j) * ID(k,j) 
               check = sqrt(b**2 + c**2)
               goto 100
500            d  = d -  AV(j) * BV(i) * ID(j,k) 
               check2 = a + b + c + d
100             continue
90         continue
80      continue
70    continue
#endif
      check = min(abs(check2),abs(check))/max(abs(check2),abs(check))

      return
      end

