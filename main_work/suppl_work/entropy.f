c-- The array "estados" stores the allowed states in rho vs mean gene length data
c-- We define 31 logarithmically distributed bins of mean gene length and
c-- 100 linearly distributed values of rho between (0,1] (101 including 0):

	integer estados(31,101)
	character*100 dummy

	bin1=6.0  !bin size (you can choose another)
	bin2=6.0  !bin start (" " " ")

	do i=1,31
	do j=1,101
	  estados(i,j)=0  !Initially "estados" are reset to zero
	end do
	end do

	open(11,file='rho_vs_gene.dat')  !Input data file

	do i=1,3
	  read(11,*) dummy   !We skip the header of rho_vs_gene.dat
	end do

	do i=1,6519
	  read(11,*) r1,r2   !reads mean gene lenght (r1) and rho (r2)

	  Lg=1+int(bin1*(log(r1)-bin2))  !converts a given mean length range in an
				     !integer bin, logarithmically distributed

	  ist=int(r2*100.)+1  !converts a given rho in a value between 0 and 100
			      !(in order to add it to "estados")

	  estados(Lg,ist)=estados(Lg,ist)+1  !We count the number of organisms
					     !in the same rho and Lg bins
	end do
	close(11)

	open(12,file='entropy.dat') !Output data file

	do i=1,31
	ista=0
	do j=1,101
	  ista=ista+estados(i,j)   !we sum the number of organisms per Lg-th bin
	end do

	entr=0.
	do j=1,101
	 if(estados(i,j).gt.0) then
	  entr=entr-((estados(i,j)/(1.*ista))*log(estados(i,j)/(1.*ista)))  !We calculate the entropy per Lg bin
	 end if
	end do

	rLg=exp(bin2+((i-0.5)/bin1))   !recover real values of mean gene length for center of Lg-th bin

	write(12,*) rLg,entr
	end do

	close(12)

	end

