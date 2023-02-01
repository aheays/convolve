!! Doppler broaden [x,y] data.  Outputs the same x grid. Works for any
!! energy/wavenumber/wavelength scaling even if it is nonuniformly
!! sampled.
program convolve_doppler
  
  implicit none
  integer*8 :: i,j,k,l       
  real*8 :: p,q,r,s
  integer*8,parameter :: WWMAX=2e7 !allocated size of data arrays
  integer*8 :: ndata                !actual length of data
  integer*8,parameter :: STDERR=0,STDOUT=6,STDIN=7,DATAFILEID=12
  character(500) :: datafile,str
  real*8, dimension(WWMAX) :: ww=0,xsin=0 !read from data file
  real*8, dimension(WWMAX) :: doppler_profile=0; ! doppler gaussian
  real*8, dimension(WWMAX) :: xsout=0; !read for output after broadening
  integer :: datafile_stat=0
  real*8 :: doppler_fwhm_except_energy,doppler_fwhm
  real*8,parameter :: FWHMS_TO_INCLUDE=5.0d0
  character(500),parameter :: OUTPUT_FORMAT='(f0.8,e18.10)' !final output format

  ! if no arguments given, print usage
  if (iargc().eq.0) then        
     write(STDERR,*) "Usage: convolve_doppler ww_indep_fwhm input_xy_file"
     write(STDERR,*) 
     write(STDERR,*) "   Broaden the data in input_xy_file according to the Doppler"
     write(STDERR,*) "   width ww_indep_fwhm*x. So the actual Doppler width is "
     write(STDERR,*) "   calculated at every point."
     write(STDERR,*) "   Output to stdout."
     write(STDERR,*) 
     stop 1
  end if
  
  ! get some input arguments
  call getarg(1,str);read(str,*) doppler_fwhm_except_energy !width apart from energy factor
  if (iargc().lt.2) then !no input file specified - read from stdin/pipe
     datafile='/dev/stdin' 
  else
     call getarg(2,datafile)
  endif

  ! open data file and read all data
  open(DATAFILEID,file=datafile,status='old',iostat=datafile_stat,action='read')
  if (datafile_stat.ne.0) then
     write(STDERR,*) "error: Cannot open input file:",datafile
     stop 1
  endif
  ndata = 0
  do while (datafile_stat.eq.0) !undatatil EOF or error
     ndata = ndata+1
     if (ndata.gt.WWMAX) then
        write(STDERR,*) "File exceeds maximum array size. Increase WWMAX and recompile."
        stop 1
     end if
     read(DATAFILEID,*,iostat=datafile_stat) ww(ndata),xsin(ndata)
  end do
  ndata = ndata-1
  close(DATAFILEID)
  
  !! loop through input data:
  xsout = 0.0d0
  do i=1,ndata
     !! Get Doppler width - at xsin so approx
     doppler_fwhm = doppler_fwhm_except_energy*ww(i)
     !! Calc Doppler gaussian inside FWHMS_TO_INCLUDE
     !! contribution to below
     do j = i-1,1,-1
        if (dabs(ww(j)-ww(i)).gt.(doppler_fwhm*FWHMS_TO_INCLUDE)) exit
        call calculate_gaussian(ww(j)-ww(i),doppler_profile(j),doppler_fwhm)
     end do
     j = j + 1
     !! contribution to above
     do k = i,ndata
        if (dabs(ww(k)-ww(i)).gt.(doppler_fwhm*FWHMS_TO_INCLUDE)) exit
        call calculate_gaussian(ww(k)-ww(i),doppler_profile(k),doppler_fwhm)
     end do
     k = k - 1
     !! normalise doppler
     doppler_profile(j:k) = doppler_profile(j:k)/sum(doppler_profile(j:k))
     !! broaden input data and add to output data
     xsout(j:k) = xsout(j:k) + xsin(i)*doppler_profile(j:k)
  end do

  !! output data file
  do i = 1,ndata
     write(STDOUT,OUTPUT_FORMAT) ww(i),xsout(i)
  end do

end program convolve_doppler

subroutine calculate_gaussian(x,y,fwhm)
  !! Calculate gaussian. Normalised to peak.
  implicit none
  real*8,intent(in)  :: x
  real*8,intent(in)  :: fwhm
  real*8,intent(out) :: y
  y = dexp(-x**2*2.772588722/fwhm**2) !2.772588722=4log(2)
end subroutine calculate_gaussian
