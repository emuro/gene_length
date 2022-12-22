	real mu,var,rshap,rscal

	rscmin=6.  !min and max for scale param. random values
	rscmax=10.
	rshmin=0.001  !min and max for shape param. random values
	rshmax=1.6

	isort=1   !if =0 generates random numbers arithmetically distributed
		  !if =1, geometrically distributed.
	inum=1000 !How many numbers you want to generate

	do i=1,inum
	  if(isort.eq.0) then

	    !Generates a random scale parameter arithmetically distributed
	    !between rscmin and rscmax:
	    rscal=rscmin+(rscmax-rscmin)*rand()	  

	    !Generates a random shape parameter arithmetically distributed
	    !between rshmin and rshmax:
	    rshap=rshmin+(rshmax-rshmin)*rand()

	  else

	    !Generates a random scale parameter geometrically distributed
	    !between rscmin and rscmax:
	    rscal=exp(log(rscmin)+(log(rscmax)-log(rscmin))*rand())

	    !Generates a random shape parameter geometrically distributed
	    !between rshmin and rshmax:
	    rshap=exp(log(rshmin)+(log(rshmax)-log(rshmin))*rand())

	  endif

c--Calculates the mean and variance for the lognormal distribution
c--from shape and scale parameters:

	   mu=exp(rscal+((rshap**2)/2.))  
	   var=(exp(rshap**2)-1.)*(exp(2.*rscal+rshap**2))
	   write(*,*) mu, var

	end do

	end 