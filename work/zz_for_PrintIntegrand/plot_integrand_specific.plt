set terminal png
set output '9-8-10-8_integrand.png'

set title "M0nu(RBME) integrand for nr-lr-npr-lpr = (start of file name)"
set key right top box
set xlabel "q [MeV]"
set ylabel "integrand"

set autoscale
xmin=0
xmax=1000
set xrange [xmin:xmax]
set arrow from xmin,0 to xmax,0 nohead linetype 9

plot "9-8-10_M0nu_integrand_F.txt" using 1:2 title "F-part" with lines, \
     "9-8-10_M0nu_integrand_GT.txt" using 1:2 title "GT-part" with lines, \
     "9-8-10-8_M0nu_integrand_T.txt" using 1:2 title "T-part" with lines