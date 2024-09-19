#!/bin/bash

awk 'BEGIN {Tsum = 0; Ns = 0};
    /T/{print $3, $9; Tsum += $3; Ns += 1};
    END {print "#", Ns, Tsum/Ns}' OSZICAR > kaka.dat

gnuplot <<EOF
set terminal png
set output 'En_Temp_t.png'
set y2tics
set ylabel "Temperature"
set y2label "Energy"
plot 'kaka.dat' u 1 w l lw .5 axes x1y1 t "Temperature",\
     'kaka.dat' u 2 w l lw .5 axes x1y2 t "Energy"
pause mouse
EOF

tail -n 1 en_temp_t.dat
# feh -xdF En_Temp_t.png

