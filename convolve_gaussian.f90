program convolve_gaussian
  !! Convolves a gaussian profile of given fwhm
  !! Requires a data input file [ww(cm-1),xs] and outputs
  !! on a new grid spacing with range the same as original.
  !!
  !! This can not be used as a direct substitute for
  !! convolve_doppler.f90 because there is an energy dependence on
  !! doppler width that is explicitly considered in
  !! convolve_doppler.f90.
  !! 
  !! To compile: gfortran -o convolve_doppler  convolve_doppler.f90
  !! for some reason -O2 of -O3 causes file reading error.
  !! 
  !! To run: SEE USAGE BELOW
  !!
  !! The convolution is done recursively where neighbouring points in the input
  !! data file are used to define a linear interpolation that can be
  !! analytically convolved with a gaussian using the following formulas:
  !! 
  !! gaussian form: f(x)=a*exp(-b*x^2),
  !! linear form: g(x)=c+d*x over limited range [x1,x2] otherwise zero
  !!
  !! Then using the identities:
  !! \int x*exp(-x^2) dx = 1/2*exp(-x^2)
  !! \int_0^x exp*-x'^2 dx' = sqrt(pi/4)*erf(x)
  !!
  !! The final convolved integral is:
  !! [f(x)*g(x)] = sqrt(pi/4/b)*a*(c+d*x)*(erf(sqrt(b)*(x-x1))-erf(sqrt(b)*(x-x2)))
  !!             + a*d/2/b*(exp(-b*(x-x2)^2)-exp(-b*(x-x1)^2))
  !!
  !! Derivation given in Journal IV pg 43
  !!
  !! This is based on an earlier version of convolve_doppler that used a stepwise
  !! rather than linear approximation
  !!
  !! 2008-03-19T14:10:45+1100, Alan Heays
  !! Based on convolve_doppler.f90 - just changed fwhm definition and inputs


  

  implicit none
  integer :: i,j,k,l,m,n=0       
  real*8 :: p,q,r,s,t,u,v,a,b,c,d
  integer,parameter :: WWMAX=100000 !max size, of array
  integer,parameter :: stderr=10,stdout=11
  character(500) :: datafile,str,output_format
  real*8, dimension(WWMAX) :: wwout=0,xsout=0; !read for output
  real*8, dimension(2) :: ww=0,xs=0 !read from data file
  integer*8 :: i1=1,i2=0        !range of output arrays currently in use
  integer :: datafileid,datafile_stat=0
  real*8 :: dww,maxdx,temperature,reduced_mass,fwhm_except_energy,fwhm

  ! prepare to write to stderr,stdout
  ! unit=0 should do this automatically but doesn't seem to work 
  ! for some reaon. 
  open(stderr,file='/dev/stderr',action='write')
  open(stdout,file='/dev/stdout',action='write')

  ! if no arguments given, print usage
  if (iargc().eq.0) then        
     write(stderr,*) "Usage: convolve_doppler fwhm dww [infile]"
     write(stderr,*) 
     write(stderr,*) "   Convolve [ww,xs] data in infile by a doppler profile of width given"
     write(stderr,*) "   in same units as input file."
     write(stderr,*) "   dww is wavenumber grid spacing of final convolved function."     
     write(stderr,*) "   Input data file must be correctly ordered and unique in energy scale."
     write(stderr,*) "   Infile can be specified as an input or provided as standard input."
     stop 1
  end if
  
  ! get some input arguments
  call getarg(1,str);read(str,*) fwhm
  call getarg(2,str);read(str,*) dww !output grid spacing


  ! open data file
  if (iargc().lt.3) then !no input file specified - read from stdin/pipe
     datafile='/dev/stdin' 
  else
     call getarg(3,datafile)
  endif
  open(datafileid,file=datafile,status='old',iostat=datafile_stat,action='read')
  if (datafile_stat.ne.0) then ;
     write(stderr,*) "error: Cannot open input file:",datafile ; 
     stop 1 ; 
  endif

  ! final output format
  output_format = '(f0.4,e18.10)'

  ! read first two input records
  read(datafileid,*,iostat=datafile_stat) ww(1),xs(1)
  read(datafileid,*,iostat=datafile_stat) ww(2),xs(2)

  ! initialise output wwout scale
  i1 = 1 ; i2 = 1 ; wwout(i1) = ww(1)

  ! main loop
  do while (datafile_stat.eq.0) !until EOF or error
     
     ! maximum range of gaussian tail considered
     maxdx = fwhm*6.0

     ! print distant records - ie already beyond range of maxdx
     ! assumes reasonable ordering of input file
     do while ((ww(1)-wwout(i1).gt.maxdx).and.(i1.lt.i2))
        write(stdout,output_format) wwout(i1),xsout(i1)
        i1 = i1 + 1
     end do

     ! add new output records - possibly need to move all output records back to beginning
     ! of storage array if approaching the end
     do while(wwout(i2)-ww(2).lt.maxdx)
        if (i2.ge.WWMAX) then   ! shift saved data back to beginning of array
           if(i1.eq.1) then
              write(stderr,*)"error: arraysize WWMAX not big enough"
              stop 1
           end if
           j=0
           do i=i1,i2
              j=j+1
              wwout(j) = wwout(i) ; xsout(j) = xsout(i) ;
           end do
           i1 = 1; i2 = j;
        end if
        i2 = i2 + 1
        wwout(i2) = wwout(i2-1) + dww
        xsout(i2) = 0.0d0
     end do
 
     ! calculate coefficients of gaussian and linear functions to be convolved together
     ! uses mean energy of input ww(1) ww(2) for doppler with, more accurately should
     ! probably be output energy wwout(i) but then would have to be recalculated at every
     ! step in loop below.
     a = 0.939437278/fwhm
     b = 2.772588722/(fwhm**2.0)
     d = (xs(2)-xs(1))/(ww(2)-ww(1))
     c = xs(1) - d*ww(1)

     ! to avoid too much duplication of calculation in loop
     p = a*dsqrt(3.141592654/4.0/b) ;  q = a*d/2.0/b

     ! loop through output points doing recursive convolution
     do i=i1,i2
        u = dsqrt(b)*(wwout(i)-ww(1))
        v = dsqrt(b)*(wwout(i)-ww(2))
        r = p*(c+d*wwout(i))*(derf(u) - derf(v))
        s = q*(dexp(-u**2.0) - dexp(-v**2.0))        
        xsout(i) = xsout(i) + r + s
     end do

     ! read next input record
     ww(1) = ww(2); xs(1) = xs(2)
     read(datafileid,*,iostat=datafile_stat) ww(2),xs(2)

  end do
     
  ! output whatever is left in output arrays - up to last input ww
  do i=i1,i2
     if (wwout(i).le.ww(1)) write(stdout,output_format) wwout(i),xsout(i)
  end do

  ! close data file
  close(datafileid) 

end program convolve_gaussian


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!An octave version of this program:
!!!!!
!!!!!
!!!!!ww=[1:0.1:50]';
!!!!!xs=ones(size(ww))+gaussian(ww,1,25);
!!!!!
!!!!!wwo = ww;
!!!!!xso = zeros(size(wwo));
!!!!!
!!!!!fwhm = 3;
!!!!!a = 0.939437278/fwhm;
!!!!!b = 2.772588722/(fwhm^2);
!!!!!maxdx = 6*fwhm;
!!!!!
!!!!!for i=1:length(xs)-1
!!!!!  xs1=xs(i);xs2=xs(i+1);
!!!!!  ww1=ww(i);ww2=ww(i+1);
!!!!!
!!!!!  d = (xs2-xs1)/(ww2-ww1);
!!!!!  c = xs1 - d*ww1;
!!!!!  p = sqrt(3.141592654/4.0/b)*a;  
!!!!!  q = a*d/2.0/b;
!!!!!
!!!!!  for j=1:length(xso)
!!!!!    if ww1-wwo(j)<maxdx&&wwo(j)-ww2<maxdx
!!!!!      u = sqrt(b)*(wwo(j)-ww1);
!!!!!      v = sqrt(b)*(wwo(j)-ww2);
!!!!!      r = p*(c+d*wwo(j))*(erf(u) - erf(v));
!!!!!      s = q*(exp(-u^2.0) - exp(-v^2.0));    
!!!!!      xso(j) = xso(j) + r + s;
!!!!!    endif
!!!!!  endfor
!!!!!endfor
!!!!!  
!!!!!      
!!!!!fig(1)
!!!!!plot(ww,xs,"1*-")
!!!!!plot(wwo,xso,"3*-")
!!!!!
!!!!!

