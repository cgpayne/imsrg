swapS=${1}-${2}-${3}
swapL=${1}-${2}-${3}-${4}

gnuplot plot_integrand.plt
make trap_integrate
./trap_integrate

mv	integrand.png \
		${swapL}_integrand.png
mv	trap_out.txt \
		${swapL}_trap_out.txt
mv	M0nu_integrand_GT.txt \
		${swapS}_M0nu_integrand_GT.txt
mv	M0nu_integrand_F.txt \
		${swapS}_M0nu_integrand_F.txt
mv	M0nu_integrand_T.txt \
		${swapL}_M0nu_integrand_T.txt
