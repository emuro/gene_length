c-- The vector "mean" stores the mean rho value per mean gene length bin
c-- The vector "var" stores the associated variance per mean gene length bin
c-- The vector "num" stores the number of organisms per mean gene length bin
c-- We define 31 logarithmically distributed bins of mean gene length:

	real mean(31),var(31),num(31)
	character*100 dummy

	bin1=6.11  !bin size (you can choose another)
	bin2=6.11  !bin start (" " " ")

	do i=1,31
	  mean(i)=0.  !Initially vectors are reset to zero
	  var(i)=0.
	  num(i)=0.
	end do

	open(11,file='rho_vs_gene.dat')  !Input data file

	do i=1,3
	  read(11,*) dummy   !We skip the header of rho_vs_gene.dat
	end do

	do i=1,6519
	  read(11,*) r1,r2   !reads mean gene lenght (r1) and rho (r2)

	  Lg=1+int(bin1*(log(r1)-bin2))  !converts a given mean length range in an
				     !integer bin, logarithmically distributed
	  mean(Lg)=mean(Lg)+r2       !sums rho values in the same Lg bin

	  num(Lg)=num(Lg)+1.  !We count the number of organisms in the same Lg bin
	end do
	close(11)

	do i=1,31
	  mean(i)=mean(i)/num(i)  !calculates the mean rho per Lg bin
	end do


	open(11,file='rho_vs_gene.dat')  !Open again data file to calculate variance

	do i=1,3
	  read(11,*) dummy   !We skip the header of rho_vs_gene.dat
	end do

	do i=1,6519
	  read(11,*) r1,r2   !reads mean gene lenght (r1) and rho (r2)

	  Lg=1+int(bin1*(log(r1)-bin2))  !converts a given mean length range in an
				     !integer bin, logarithmically distributed
	  var(Lg)=var(Lg)+((r2-mean(Lg))**2)       !sums (rho-mean)^2 values in the same Lg bin

	end do
	close(11)

	do i=1,31
	  var(i)=var(i)/num(i)  !calculates the variance per Lg bin
	end do

	open(12,file='variance.dat') !Output data file

	do i=1,31

	rLg=exp(bin2+((i-0.5)/bin1))   !recover real values of mean gene length for center of Lg-th bin

	write(12,*) rLg,var(i)
	end do

	close(12)

	end

