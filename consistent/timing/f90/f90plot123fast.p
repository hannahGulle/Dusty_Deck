set title "F90 Dusty Deck using\nO1, O2, O3, and Ofast Optimization Flags"
#plot [20:92] [0:150] "f90timing_O0" u 1:2 with linespoints

replot "f90timing_O1" with linespoints
replot "f90timing_O2" with linespoints
replot "f90timing_O3" with linespoints
replot "f90timing_Ofast" with linespoints

set xlabel "Maxsize"
show xlabel
set ylabel "Time (seconds)"
show ylabel
set key inside left
show key

set output "f90time_0-fast.ps"
set terminal postscript enhanced color landscape
replot
