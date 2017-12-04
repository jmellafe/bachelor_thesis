
set term png
set output "voter_model_gillespie.png"

set xlabel "Tiempo"
set ylabel "Porcentaje en estado 1"
set yrange [0:1]
p "/home/alex/CLionProjects/tfg/results/voter_model_gillespie.dat" u 1:2 t ""