set title "C Dusty Deck Comparing O0 with\nO1, O2, O3, and Ofast Optimization Flags"
plot [65:90] [50:120] "ctiming_O0" u 1:2 with linespoints

set xlabel "Maxsize"
show xlabel
set ylabel "Time (seconds)"
show ylabel
set key inside left
show key

replot "ctiming_O1" with linespoints
replot "ctiming_O2" with linespoints
replot "ctiming_O3" with linespoints
replot "ctiming_Ofast" with linespoints

set output "ctime_0-fast_zoomed.ps"
set terminal postscript enhanced color landscape
replot
