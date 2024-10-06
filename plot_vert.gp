# gnuplot script to plot vert.dat

# Set title and labels
set title "Vertically Averaged Distribution"
set xlabel "x"
set ylabel "Average u(x)"

# Set output file type and name
set terminal png
set output "vert.png"

# Plot vert.dat
plot "vert.dat" with lines