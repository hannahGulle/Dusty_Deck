      program dusty 

* Assigns 50 to constant MAXDIM -- Can't be changed
      parameter (MAXDIM = 50)
      parameter (pi = acos(-1.0))

* Creates an integer array of maxdimension size AND
* a changeable copy of MAXDIM known as N
      integer IA(MAXDIM), N, i, j, k

* Creates double rank 1 arrays of MAXDIM
      double precision AV(MAXDIM), BV(MAXDIM), CV(MAXDIM)

* Creates double rank 2 arrays of MAXDIM
      double precision OP(MAXDIM,MAXDIM), ID(MAXDIM,MAXDIM)
      double precision AM(MAXDIM,MAXDIM), BM(MAXDIM,MAXDIM)
      double precision CM(MAXDIM,MAXDIM), DM(MAXDIM,MAXDIM)
      
      double precision check, BOT, TOP, HOLDA, HOLDB, TRACE3
      real start, finish

* Function which acts as an external subprogram
      double precision xij, xji, yij, yji, zij, zji, trigij, trigji

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

      do i = 1, N
        AV(i) = bessel_jn(0,dble((rand() *
     +               (-1)**(mod(int(10*rand()),N)))))
      enddo

      do i = 1, N
        BV(i) = bessel_jn(1,dble((rand() * 
     +               (-1)**(mod(int(10*rand()),N)))))
      enddo


* If AV == BV, check += AV*BV == AV^2
* AV AND BV may not be necessary.
* idcheck with 5 paramters: 
      check = 0.0
      do i = 1, N
        ival = N
        check = check + ( AV(i) * BV(i) )
        call idcheck(ival,check,AV,BV,ID)
      enddo


* Compute |AV><BV|
* Building a matrix??

        do i = 1, N
          do j = 1, N
                call idcheck(N, check, AV,BV,ID)
                if ( check .lt. 0.5 ) then
                        OP(i,j) = ( AV(j) * BV(i) ) / BV(j)
                else
                        OP(i,j) = ( AV(i) * BV(j) ) / BV(i)
                endif
                IA(i) = i
           enddo
        enddo

* Mod Mod??
        do i = 1, N
                do j = 0, i, 8
                        IA(i) = mod(mod(i+j,N),N)+1
                enddo
        enddo


* Unfurl
        do i = 1, N
         call idcheck(N,check,AV,BV,ID)
         CV(IA(i)) = (AV(IA(i)) + BV(IA(i))) / check
        enddo


* Unfurl
        do i = 2, N
         call idcheck(N,check,AV,BV,ID)
         AV(i) = AV(i-1) * BV(i) + CV(i)
        enddo


* lots of repeated code...
* Consider making a 1 layer loop.
* Pull out AV and BV assignments from the if
        do i = 1, N
            call idcheck(N,check,AV,BV,ID)
            if ( .not.(check .gt. 0.5) ) then
                    do j = 1, N    
                        BOT = OP(i,j)
                        TOP = AV(j) * BV(j)
                        HOLDA = AV(j)

                        xij = dble( IA(i) ) - dble( IA(j) )
                        xji = dble( IA(j) ) - dble( IA(i) )

                        yij = dble( IA(i) ) + dble( IA(j) )
                        yji = dble( IA(j) ) + dble( IA(i) )

                        zij = exp( sin( sqrt( xij**2 + yij**2 ) * pi ))
                        zji = exp( sin( sqrt( xji**2 + yji**2 ) * pi ))

                        trigij = xij+yij+log10(abs(1+zij+(xij*yij*zij)))
     +                          / (abs(xij)+abs(yij))
                        trigji = xji+yji+log10(abs(1+zji+(xji*yji*zji)))
     +                          / (abs(xji)+abs(yji))
                   
                        AV(j) = BV(j) - CV(j) / (TOP-BOT) * ID(j,j)
                        BV(j) = HOLDA - CV(j) / (TOP-BOT) * ID(i,i)
                        AM(i,j) = AV(j) / trigij
                        BM(i,j) = BV(j) / trigji
                     enddo
                else
                    do j = 1, N    
                        BOT = OP(i,j)
                        TOP = AV(j) * BV(j)
                        HOLDA = AV(j)

                        xij = dble( IA(i) ) - dble( IA(j) )
                        xji = dble( IA(j) ) - dble( IA(i) )

                        yij = dble( IA(i) ) + dble( IA(j) )
                        yji = dble( IA(j) ) + dble( IA(i) )

                        zij = exp( sin( sqrt( xij**2 + yij**2 ) * pi ))
                        zji = exp( sin( sqrt( xji**2 + yji**2 ) * pi ))

                        trigij = xij+yij+log10(abs(1+zij+(xij*yij*zij)))
     +                          / (abs(xij)+abs(yij))
                        trigji = xji+yji+log10(abs(1+zji+(xji*yji*zji)))
     +                          / (abs(xji)+abs(yji))
                   
                        AV(j) = BV(j) + CV(j) / (TOP-BOT) * ID(i,i)
                        BV(j) = HOLDA + CV(j) / (TOP-BOT) * ID(j,j)
                        AM(i,j) = AV(j) * trigij
                        BM(i,j) = BV(j) * trigji
                enddo
           endif
        enddo

* Consider making a 1 layer loop through N^3
* Consider unfurling loop
        do i = 1, N
            do j = 1, N
                 CM(i,j) = 0.0

                 if ( .not.( i .lt. j )) then
                    do k = 1, N
                         CM(i,j) = CM(i,j) + AM(i,k) * BM(k,j) / check 
                    enddo
                 else
                    do k = 1, N
                         CM(i,j) = CM(i,j) - AM(i,k) * BM(k,j) / check
                    enddo
                 endif
           enddo
       enddo

Consider making a 1/2 layer loop
      do i = 1, N
        do j = 1, N
            sum = 0.0
                do k = 1, N
                    sum = sum + CM(i,k) * AM(j,k)
                enddo
            DM(i,j) = sum
        enddo
      enddo


* Consider making a 1 layer loop
      do i = 1, N
         do j = 1, N
            CM(i,j) = DM(i,j)
         enddo
      enddo


* Consider condensing the loop layers
      do i = 1, N
         do j = 1, N
            sum = 0.0
               do k = 1, N
                  sum = sum - CM(i,k) * BM(j,k)
               enddo
               DM(i,j) = sum
         enddo
      enddo


* Consider condensing the loop layers 
* AND unfurling
      HOLDA = abs(AM(1,1))
      HOLDB = abs(BM(1,1)) 
      do i = 1, N
         do j = 1, N
          HOLDA = max(HOLDA,abs(AM(i,j))) 
          HOLDB = max(HOLDB,abs(BM(i,j))) 
         enddo
      enddo
      TRACE3 = 0.0
       

* Unfurl
      do i = 1, N
        TRACE3 = TRACE3 + (AM(IA(i),IA(i)) + BM(IA(i),IA(i)) 
     +                  - DM(IA(i),IA(i))) / (HOLDA * HOLDB)
      enddo

      call cpu_time(finish) 

      print *, 'Final trace = ', trace3, ' and IDCHECK ', check
      print *, '-- RUNTIME -> ', finish-start, ' seconds'
      end
        
***************************************************************
* External subprogram TRIG 
* Deleted and integrated without a function call in the program

***************************************************************
* subroutine idcheck does what??
      subroutine idcheck(N,check,AV,BV,ID)

      double precision AV(*), BV(*), ID(N,*)
      double precision l2
      double precision check, check2
      double precision a, b, c, d 

* Condense to 1 layer loop
      do i = 1, N
         do j = 1, N
          if ( i .ne. j ) then
             ID(i,j) =  cos(check+2.0*i*acos(-1.0)/N)+
     C                  2.0*sin(check+ 2.0*j*acos(-1.0)/N)
          elseif ( i .eq. j ) then 
             if ( (( AV(i) .eq. 0) .and. (BV(j) .eq. 0)) .or.
     C            (( AV(i) .gt. 0) .and. (BV(j) .gt. 0)) ) then
                ID(i,j) = 1.0
             elseif (( AV(i) .gt. 0 ) .and. ( BV(j) .lt. 0 )) then
               ID(i,j) = -1.0
             elseif (( AV(i) .lt. 0 ) .and. ( BV(j) .lt. 0 )) then
               ID(i,j) = 1.0
             elseif (( AV(i) .lt. 0 ) .and. ( BV(j) .gt. 0 )) then
               ID(i,j) = -1.0
             endif
          endif
         enddo
      enddo


* Unfurl
      l2 = 0.0
      do i = 1, N
         l2 = l2 + AV(i)**2
      enddo
      l2 = sqrt(l2)

      
* Unfurl
      do i = 1, N
         AV(i) = AV(i) / l2
      enddo


* Unfurl
      l2 = 0.0
      do i = 1, N
         l2 = l2 + BV(i)**2
      enddo
      l2 = sqrt(l2)


* Unfurl
      do i = 1, N
         BV(i) = BV(i) / l2
      enddo


* "computed goto" wtf rewrite
* Condense loop layers
* find out what the int(mod(i+j+k,4)+1) computes
* might better work split up given it's behavior
      a = 0.0D0
      b = 0.0D0
      c = 0.0D0
      d = 0.0D0
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
      check = min(abs(check2),abs(check))/max(abs(check2),abs(check))

      return
      end

