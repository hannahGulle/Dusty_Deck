set title "C Dusty Deck using\nO1, O2, O3, and Ofast Optimization Flags"
plot [20:95] [0:150] "ctiming_O1" u 1:2 with points

replot "ctiming_O2" with lines
replot "ctiming_O3" with linespoints
replot "ctiming_Ofast" with lines

set output "ctime_1.ps"
set terminal postscript enhanced color landscape
replot
