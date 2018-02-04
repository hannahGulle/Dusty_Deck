set title "F77 Dusty Deck Comparing O0 with\nO1, O2, O3, and Ofast Optimization Flags"
plot [70:90.5] [50:145] "f77timing_O0" u 1:2 with linespoints

set xlabel "Maxsize"
show xlabel
set ylabel "Time (seconds)"
show ylabel
set key inside left
show key

replot "f77timing_O1" with linespoints
replot "f77timing_O2" with linespoints
replot "f77timing_O3" with linespoints
replot "f77timing_Ofast" with linespoints

set output "f77time_0-fast_zoomed.ps"
set terminal postscript enhanced color landscape
replot
