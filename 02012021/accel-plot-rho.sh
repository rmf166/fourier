#!/bin/sh
rm accel*.pdf
for sc in 2 3 4 5
do
  for sw in 1 2 3 4
  do
    for prb in 1 2 3 4 5 6
    do
      rm plot.p
      echo 'set autoscale' >> plot.p
      echo 'unset logscale; set logscale x' >> plot.p
      echo 'unset label' >> plot.p
      echo 'set xtic auto' >> plot.p
      echo 'set ytic auto' >> plot.p
      if [ ${sw} == 1 ]
      then
        if [ ${sc} == 2 ]; then
          echo 'set title "p = '${prb}', c = 0.8"' >> plot.p
        elif [ ${sc} == 3 ]; then
          echo 'set title "p = '${prb}', c = 0.9"' >> plot.p
        elif [ ${sc} == 4 ]; then
          echo 'set title "p = '${prb}', c = 0.99"' >> plot.p
        elif [ ${sc} == 5 ]; then
          echo 'set title "p = '${prb}', c = 0.9999"' >> plot.p
        fi
      else
        if [ ${sc} == 2 ]; then
          echo 'set title "p = '${prb}', c = 0.8, s = '${sw}'"' >> plot.p
        elif [ ${sc} == 3 ]; then
          echo 'set title "p = '${prb}', c = 0.9, s = '${sw}'"' >> plot.p
        elif [ ${sc} == 4 ]; then
          echo 'set title "p = '${prb}', c = 0.99, s = '${sw}'"' >> plot.p
        elif [ ${sc} == 5 ]; then
          echo 'set title "p = '${prb}', c = 0.9999, s = '${sw}'"' >> plot.p
        fi
      fi
      echo 'set xlabel "{/Symbol t} (mfp)" enhanced' >> plot.p
      echo 'set ylabel "{/Symbol r}" enhanced' >> plot.p
      echo 'set yr [0:1]' >> plot.p
      echo 'set xr [0.01:100]' >> plot.p
      echo 'plot "result-p'${prb}'-s'${sw}'-LC.dat" using 1:'${sc}' title "LC(CMFD)"  with lines linetype 1 lc rgb "blue", \' >> plot.p
      echo '     "result-p'${prb}'-s'${sw}'-LD.dat" using 1:'${sc}' title "LD(CMFD)"  with lines lw 3 linetype 2 lc rgb "blue", \' >> plot.p
      echo '     "result-lp'${prb}'-s'${sw}'-LC.dat" using 1:'${sc}' title "LC(lpCMFD)" with lines linetype 1 lc rgb "green", \' >> plot.p
      echo '     "result-lp'${prb}'-s'${sw}'-LD.dat" using 1:'${sc}' title "LD(lpCMFD)" with lines lw 3 linetype 2 lc rgb "green", \' >> plot.p
      echo '     "../../pfa/results/result-p'${prb}'-s'${sw}'-LC.dat" using 1:'${sc}' title "LC(pCMFD)" with lines linetype 1 lc rgb "red", \' >> plot.p
      echo '     "../../pfa/results/result-p'${prb}'-s'${sw}'-LD.dat" using 1:'${sc}' title "LD(pCMFD)" with lines lw 3 linetype 2 lc rgb "red"' >> plot.p
      echo 'set key left top' >> plot.p
      echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
      echo 'set output "accel-plot-'${sc}'-'${prb}'-'${sw}'.pdf"' >> plot.p
      echo 'replot' >> plot.p
      echo 'set terminal x11' >> plot.p
      gnuplot plot.p
    done
  done
done
pdfunite accel-plot*.pdf accel-final.pdf
rm accel-plot*.pdf
