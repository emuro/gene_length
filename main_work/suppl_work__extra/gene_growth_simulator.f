C-- SIMULATOR OF GENE GROWTH BY A MULTIPLICATIVE STOCHASTIC FACTOR

C-- "gen" stores the lengths of 5000 genes for their coding and noncoding parts
C-- Thus, gen(i,1) is the length of the CDS part for the i-th gen, 
C-- and gen(i,2) the length of the nCDS. 
C-- The total lenght of i-th gen is therefore gen(i,1)+gen(i,2)

	integer*8 gen(5000,2)
	integer length,dummy
	character*10 dumm
	real media,varianza,Lc

C--Defining valules for some variables:

	isim=10  !Number of different simulations you want to perform

	iter=100000  !Number of iterations per simulation

	idist=1  !0 if the initial distribution is a uniform distribution
		 !otherwise it is a lognormal distribution

	Lc=1500. !Here we define the critical mean gene length value

	rmin=1.   !rmin and rmax are the limits for the multiplicative random stochastic fact
	rmax=2.   


C--Defining the initial distribution of genes:

	do iu=11,11+isim !Loop for the different simulations

	  if(idist.eq.0) then  !If 1 we generate a uniform distribution with mean 480
	    do i=1,5000
	      gen(i,1)=int(rand()*480.*2.)
	      gen(i,2)=0                   !Initially there are not nCDS
	    end do
	  else		
	    open (10,file='geneslognormal_480.dat')  !File storing an initial lognormal dist.
						     !with mean 480
	    read(10,*) dumm  !We skip the header
	    do i=1,5000
	      read(10,*) length
	      gen(i,1)=length
	      gen(i,2)=0                   !Initially there are not nCDS
	    end do

	    close(10)
	  endif


C--Starting the simulation

	  open(iu)   		!We store the result in the file fort.iu (iu=11, 12...)
				!(each file is an independent simulation)
	  do j=1,iter  !Iterations

	    rmedia=0.
	    do i=1,5000
	      rmedia=rmedia+1.*gen(i,1)+1.*gen(i,2)  !The sum of CDS and nCDS parts is total length of gene
	    end do
	    rmedia=rmedia/5000.			!We calculate the current mean gene length

	    ii=abs(rand()*5000.)+1		!Randomly we choose a gene
	    length=gen(ii,1)+gen(ii,2)		!Current total length of that gene

	    ranfac=(rmin+(rand()*(rmax-rmin))) 	!random stochastic factor between rmin and rmax

	    newlength=nint((1.*length)*ranfac)  !The new gene length is the old length 
						!multiplied by this random factor

	    if(rmedia.le.Lc) then		!If the mean gene length is under Lc
	      gen(ii,1)=newlength		!the growth is in the CDS part
	    else
	      gen(ii,2)=newlength-gen(ii,1)	!otherwise, in the nCDS part
	    endif

	    media=0.
	    rho=0.
	    do i=1,5000
	      media=media+gen(i,1)+gen(i,2)
	      rho=rho+gen(i,2)
	    end do
	    media=media/5000.  !We compute the new mean gene length
	    rho=rho/5000.
	    rho=rho/media      !And the new mean rho value

	    varianza=0.
	    do i=1,5000
	      varianza=varianza+((gen(i,1)+gen(i,2))-media)**2
	    end do
	    varianza=varianza/5000.  !We compute the new variance

	    write(iu,*) media,varianza,rho  !Store in fort.iu the evolution of mean, variance and rho

	    if(rmedia.gt.68287.) then   !Stop iterations when reaching Homo sapiens mean gene length
	      if(iu.eq.11) then
	        open(10,file="final_distribution.dat")
	        do i=1,5000
	          write(10,*) gen(i,1),gen(i,2)  !In one of the simulations, as an example, we store
	        end do                           !the CDS and nCDS lengths of the simulated genes
	        close(10)
	      endif
	      goto 1
	    endif

	  end do  !Ends loop of iterations

1	  continue
	  close(iu)
	end do     !Ends loop of simulations

	end



	