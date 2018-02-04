set title "C Dusty Deck using\nO1, O2, O3, and Ofast Optimization Flags"
plot [20:95] [0:150] "ctiming_O0" u 1:2 with points

replot "ctiming_O1" with linespoints
replot "ctiming_O2" with lines
replot "ctiming_O3" with linespoints
replot "ctiming_Ofast" with lines

set output "ctime_0-fast.ps"
set terminal postscript enhanced color landscape
replot
