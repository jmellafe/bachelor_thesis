
set term png
set output "Tiempo consenso.png"

set xlabel "N"
set ylabel "Tiempo consenso"
set logscale y
set logscale x
p "/home/alex/CLionProjects/tfg/results/tiempo_consenso.dat" u 1:2 t ""