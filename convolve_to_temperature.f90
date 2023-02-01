program convolve_to_temperature
  !! arguments: T iso res datafile outputfile 
  !!
  !! Input datafile in format [E,xs,Jdd,Jd]. Each [Jdd,Jd] rotational
  !! pair must be as a group in the file, be energy sorted and cover
  !! exactly the same energy range as every other pair.
  !!
  !! Each rotational transition is convolved with an energy dependent
  !! doppler profile and scaled by a boltzmann factor. Doppler width
  !! and boltzmann factor are dependent on temperature(T), in K, and
  !! isotope(iso), which may be 14,1415 or 15 meaning 14N2, 14N15N or
  !! 15N2.
  !!
  !! Output data is in format [E,xs] and is evenly spaced with the
  !! given resolution res and limits defined by the input data.
  !!
  !! Is quite fast but with large input files may run out of
  !! preallocated array space. In which case either cut the input data
  !! into chunks or increase size of parameter WWMAXFINALOUT and
  !! recompile.

  !! To compile: 
  !! gfortran -o convolve_to_temperature  convolve_to_temperature.f90
  !!
  !! To run: 
  !! ./convolve_to_temperature T iso res datafile outputfile
  !!
  !! Formula for doppler (gaussian) FWHM:
  !! fwhm = 2*6.331e-8*sqrt(T(K)*32/RM(atomicunits))*center_energye(cm-1)
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
  !! Derivation given in Journal IV pg 43.

  !! 2008-03-19T14:13:50+1100 Alan Heays
  !! This is based on an earlier version of convolve_doppler that used a stepwise
  !! rather than linear approximation.
  !!
  !! 2008-04-08T16:16:42+1000 Alan Heays
  !! Found that too close together ww points (<1.0e-4cm-1) caused a problem
  !! so I added an if statement so that in this the intervening
  !! cross-section is treated as a constant - the mean value of the
  !! boundary points.
  !! 
  !! 2008-11-20T22:03:56+1100 Alan Heays
  !! Heavily modified to process an input datafile containing multiple
  !! rotational transitions.

  implicit none

  integer :: i,j,k,l,m,n=0       
  real*8 :: p,q,r,s,t,u,v,a,b,c,d
  integer,parameter :: WWMAX=1e7 !max size, of array
  integer,parameter :: WWMAXFINALOUT=1e8 !max size, of output array
  integer,parameter :: DATAFILEID=10
  integer, parameter:: error_unit=0, output_unit=6 !could be replaced by intrinsic module ISO_FORTRAN_ENV
  character(500) :: datafile,str,output_format,outputfile
  real*8, dimension(WWMAX) :: wwout=0,xsout=0; !read for output
  real*8, dimension(WWMAXFINALOUT) :: wwfinalout=0,xsfinalout=0; !read for output
  real*8 :: wwfinalstart                         !beginning of output range
  real*8, dimension(2) :: ww=0,xs=0 !read from data file
  integer*8 :: i1=1,i2=0        !range of output arrays currently in use
  integer :: datafilestat=0
  real*8 :: dww,maxdx,temperature,mass,fwhm_except_energy,fwhm
  real*8 :: boltzmann_factor
  integer :: Jpp,Jp,Jppold,Jpold,isotope,vpp,vppold
  character(500) :: line


  !! open standard error

  !! if incorrect number arguments given, print usage
  if (iargc().ne.4) then        
     write(error_unit,*) "Usage: convolve_to_temperature T iso dww datafile"
     write(error_unit,*) 
     write(error_unit,*) "   Convolve CSE cross section in datafile (format [E,xs,vpp,Jdd,Jd] with"
     write(error_unit,*) "   a doppler profile for temperature T (K) and isotope iso (14/1415/15)."
     write(error_unit,*) "   Also multiplies each rotational transition defined by Jdd and Jd"
     write(error_unit,*) "   by the ground state boltzmann distribution."     
     write(error_unit,*) "   Outputs to stdout in format [E,xs] where energy E gridpoints are"
     write(error_unit,*) "   are separated by dww."
     write(error_unit,*) 
     write(error_unit,*) "   It is required that input data be ordered in Jdd,Jd first and then" 
     write(error_unit,*) "   energy and that each Jdd, Jd pair starts at the same energy."
     stop 1
  end if


  !! TEMPORARY WARNING
  write(*,*) 'WARNING: last time I used this program'
  write(*,*) '2014-11-19T17:29:03+0100, I noticed there was some disagreement'
  write(*,*) 'between the endoced N2 Boltzmann factors and those calculated in'
  write(*,*) 'spectra.py'
  

  !! Get inputs, open files. Assumes input data file [ww,xs,Jdd,Jd],
  !! nonstandard grid. Must be in correct ww order for each Jdd, Jd
  !! pair and each pair has the same ww range.
  call getarg(1,str);read(str,*) temperature !temperature - Kelvin
  call getarg(2,str);read(str,*) isotope !isotope 14,1415,15
  call getarg(3,str);read(str,*) dww !output grid spacing
  call getarg(4,str);read(str,*) datafile ![ww,xs,Jdd,Jd]

  !! determine mass (atomic mass units) from isotope
  if (isotope.eq.14) then ; mass = 28 !preset 14N14N
  else if (isotope.eq.1415) then ;  mass = 29 !preset 14N15N
  else if (isotope.eq.15) then ; mass = 30 !preset 15N15N
  else ; write(error_unit,*) "Bad isotope specification, please use 14,1415, or 15"; stop 1; end if
  
  !! fwhm without energy dependence
  fwhm_except_energy = 2*6.331e-8*dsqrt(temperature*32/mass) 
  
  !! open output file

  !! write header
  write(output_unit,*) "## Generated by convolve_to_temperature."
  write(output_unit,*) "## Input arguments: temperature =",temperature
  write(output_unit,*) "## isotope =",isotope
  write(output_unit,*) "## dww =",dww
  write(output_unit,*) "## input data file =",trim(datafile)

  !! open data file
  open(DATAFILEID,file=datafile,status='old',iostat=datafilestat,action='read')
  if (datafilestat.ne.0) then 
     write(error_unit,*) "error: Cannot open input file:",datafile ; stop 1 ; 
  endif

  !! final output format
  output_format = '(f0.4xe12.4e3)'

  !! read first two input records
  call read_data_skip_comments(datafileid,datafilestat,ww(1),xs(1),vpp,Jpp,Jp)
  call read_data_skip_comments(datafileid,datafilestat,ww(2),xs(2),vpp,Jpp,Jp)

  !! boltzmann_factor for first rotational transition
  call calculate_boltzmann_factor(temperature,isotope,vpp,Jpp,boltzmann_factor)
  
  !! record which rotational transition this is so know when a new one
  !! starts in the input data file
  Jpold=Jp;Jppold=Jpp;vppold=vpp

  !! initialise output wwout scale  for this rotational transition
  wwout=0;xsout=0;
  i1 = 1 ; i2 = 1 ; 
  wwout(i1) = ww(1)

  !! start of energy scale in output file
  wwfinalstart = wwout(1)       

  !! main loop - good luck on figuring out what goes on in here it is!!!!!!!!!!!
  !! ridiculously complicated.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do while (datafilestat.eq.0) !until EOF or error

     !! energy dependent fwhm and maximum range of gaussian tail considered
     !! technically fwhm should be at energies wwout(i) - but not much difference 
     fwhm = fwhm_except_energy*(ww(1)+ww(2))/2.0 
     maxdx = fwhm*6.0

     !! print distant records - ie already beyond range of maxdx
     do while ((ww(1)-wwout(i1).gt.maxdx).and.(i1.lt.i2))
        j=idnint((wwout(i1)-wwfinalstart)/dww)+1
        if (j.gt.WWMAXFINALOUT) then
           write(error_unit,*) "error: too much data, reduce energy range or",&
                & " increase WWMAXFINALOUT in source code and recompile"
           stop 1
        end if

        if (j.lt.1) then
           write(error_unit,*) "error: input data at wavenumber below requested range.",&
                & "This is not allowed. Or modify source to check for this."
           stop 1
        end if

        !! add this rotation to final output data
        wwfinalout(j) = wwout(i1)
        xsfinalout(j) = xsfinalout(j)+xsout(i1)*boltzmann_factor
        i1 = i1 + 1
     end do

     !! add new output records - possibly need to move all output records back to beginning
     !! of storage array if approaching the end
     do while(wwout(i2)-ww(2).lt.maxdx)
        if (i2.ge.WWMAX) then   ! shift saved data back to beginning of array
           if(i1.eq.1) then
              write(error_unit,*) "error: arraysize WWMAX not big enough"
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

     !! calculate coefficients of gaussian and linear functions to be convolved together
     !! uses mean energy of input ww(1) ww(2) for doppler width, more accurately should
     !! probably be output energy wwout(i) but then would have to be recalculated at every
     !! step in loop below.
     a = 0.939437278/fwhm
     b = 2.772588722/(fwhm**2.0)
     
     if (dabs(ww(2)-ww(1)).gt.1.0d-4) then !this 'if' is to combat some kind of machine precision
        d = (xs(2)-xs(1))/(ww(2)-ww(1))    !problem when adjacent ww are too close together
        c = xs(1) - d*ww(1)                !If this occurs assume constant valued xs,
     else                                  !ie. no slope to linear function 
        d = 0.0d0                          
        c = (xs(1)+xs(2))/2
     endif

     !! to avoid too much duplication of calculation in loop
     p = a*dsqrt(3.141592654/4.0/b) ;  q = a*d/2.0/b ; t = dsqrt(b)

     !! loop through output points doing recursive convolution
     do i=i1,i2
        u = t*(wwout(i)-ww(1))
        v = t*(wwout(i)-ww(2))
        r = p*(c+d*wwout(i))*(derf(u) - derf(v))
        s = q*(dexp(-u**2.0) - dexp(-v**2.0))        
        xsout(i) = xsout(i) + r + s
     end do

     !! read next input record
     ww(1) = ww(2); xs(1) = xs(2)
     call read_data_skip_comments(datafileid,datafilestat,ww(2),xs(2),vpp,Jpp,Jp)

     !! end of rotation or end of file,  up to last input ww
     if((datafilestat.ne.0).or.(Jpp.ne.Jppold).or.(Jp.ne.Jpold).or.(vpp.ne.vppold)) then

        !output whatever is left in output arrays 
        do i=i1,i2
           if (wwout(i).lt.ww(1))  then
              j=idnint((wwout(i)-wwfinalstart)/dww)+1
              wwfinalout(j) = wwout(i)
              xsfinalout(j) = xsfinalout(j)+xsout(i)*boltzmann_factor
           endif
        end do

        !! calculate a new boltzmann factor
        call calculate_boltzmann_factor(temperature,isotope,vpp,Jpp,boltzmann_factor)
   
        !! clean up for next loop
        ww(1)=ww(2);xs(1)=xs(2)
        call read_data_skip_comments(datafileid,datafilestat,ww(2),xs(2),vpp,Jpp,Jp) 
        Jpold=Jp;Jppold=Jpp;vppold=vpp
        i1 = 1 ; i2 = 1 ; 
        wwout=0;xsout=0;
        wwout(i1) = ww(1)
 
    end if
  end do
  !! END OF MAIN LOOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! write output file
  do i=1,WWMAXFINALOUT
     !! stop if no cross-section ever added this far into the array
     if (wwfinalout(i).le.0.0d0) exit
     write(output_unit,output_format) wwfinalout(i),xsfinalout(i)
  end do
  
end program convolve_to_temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_boltzmann_factor(temperature,isotope,vwant,Jwant,boltzmann_factor)
  implicit none
  integer,parameter :: LENGTH_J_DATA=51,LENGTH_V_DATA=2
  real*8 :: temperature,partition,k,h,c,constant
  real*8 :: energy(LENGTH_J_DATA,LENGTH_V_DATA)
  real*8 :: nuclear_spin_weight(LENGTH_J_DATA),boltzmann_factor
  integer :: isotope,J,v,i,jtmp,vwant,Jwant

  k = 1.3806505e-23; !J/K.
  h = 6.6260693e-34; !Js
  c = 2.99792458e8; !m/s
  constant = 100.0d0*c*h/k/temperature;

  !! Hardcoded energy levels and nuclear spin weightings.
  if (isotope.eq.14) then
     !! v=0 energy levels
     energy(:,1)=(/ 1175.73663136,1179.71575319,1187.67385905,1199.61067333,1215.52578263,&
          1235.41863576,1259.28854377,1287.13467995,1318.95607986,1354.75164134,&
          1394.52012455,1438.26015202,1485.97020865,1537.64864174,1593.29366108,&
          1652.90333896,1716.47561021,1784.00827227,1855.49898523,1930.94527192,&
          2010.34451791,2093.69397162,2180.9907444,2272.23181052,2367.41400736,&
          2466.53403537,2569.58845824,2676.57370293,2787.48605975,2902.3216825,&
          3021.07658852,3143.74665877,3270.32763796,3400.81513465,3535.20462133,&
          3673.49143452,3815.67077491,3961.73770744,4111.68716142,4265.51393066,&
          4423.21267355,4584.77791323,4750.20403767,4919.48529981,5092.61581768,&
          5269.58957456,5450.40041904,5635.04206524,5823.50809286,6015.79194738,&
          6211.88694017 /)
     !! v=1 energy levels
     energy(:,2)=(/ 3505.62111, 3509.56549, 3517.45412, 3529.28673,&
          3545.06289, 3564.78205, 3588.44353, 3616.04650, 3647.58999,&
          3683.07289, 3722.49397, 3765.85185, 3813.14501, 3864.37178,&
          3919.53039, 3978.61890, 4041.63524, 4108.57720, 4179.44244,&
          4254.22849, 4332.93271, 4415.55236, 4502.08455, 4592.52623,&
          4686.87425, 4785.12529, 4887.27592, 4993.32256, 5103.26149,&
          5217.08886, 5334.80068, 5456.39282, 5581.86101, 5711.20087,&
          5844.40786, 5981.47730, 6122.40438, 6267.18417, 6415.81157,&
          6568.28138, 6724.58824, 6884.72667, 7048.69105, 7216.47560,&
          7388.07445, 7563.48157, 7742.69078, 7925.69579, 8112.49017,&
          8303.06734, 8497.42062/)
     nuclear_spin_weight(1:LENGTH_J_DATA:2)=2
     nuclear_spin_weight(2:LENGTH_J_DATA:2)=1
  elseif (isotope.eq.1415) then
     energy(:,1)=(/1156.091, 1159.9381606, 1167.6323531, 1179.1733198,&
          1194.5606744, 1213.7939017, 1236.8723579,  1263.7952702,&
          1294.5617372, 1329.1707288, 1367.6210863,  1409.9115222,&
          1456.0406204, 1506.0068362, 1559.8084965,  1617.4437994,&
          1678.9108148, 1744.2074841, 1813.3316202,  1886.2809078,&
          1963.0529032, 2043.6450347, 2128.0546022,  2216.2787777,&
          2308.314605,  2404.1590002, 2503.8087513,  2607.2605185,&
          2714.5108343, 2825.5561036, 2940.3926037,  3059.0164844,&
          3181.4237681, 3307.6103499, 3437.5719977,  3571.3043521,  &
          3708.802927,  3850.0631091, 3995.0801585,  4143.8492083,&
          4296.3652652, 4452.6232093, 4612.6177944,  4776.343648,&
          4943.7952713, 5114.9670395, 5289.853202,   5468.4478822,&
          5650.7450779, 5836.7386615, 6026.4223796/)
     energy(:,2)=(/3447.403371, 3451.217515, 3458.845672, 3470.287586, 3485.542868,&
          3504.611005, 3527.491349, 3554.183128, 3584.685440, 3618.997251,&
          3657.117402, 3699.044603, 3744.777436, 3794.314354, 3847.653681,&
          3904.793612, 3965.732215, 4030.467428, 4098.997059, 4171.318791,&
          4247.430175, 4327.328637, 4411.011470, 4498.475843, 4589.718795,&
          4684.737237, 4783.527951, 4886.087593, 4992.412689, 5102.499638,&
          5216.344710, 5333.944050, 5455.293672, 5580.389465, 5709.227190,&
          5841.802478, 5978.110837, 6118.147645, 6261.908152, 6409.387484,&
          6560.580638, 6715.482484, 6874.087767, 7036.391103, 7202.386984,&
          7372.069773, 7545.433710, 7722.472905, 7903.181344, 8087.552888,&
          8275.581270 /)
     nuclear_spin_weight(1:LENGTH_J_DATA)=1
  elseif (isotope.eq.15) then
     energy(:,1)=(/1136.1035572,1139.8187372,1147.2489872,1158.3940672,1173.2536072,&
            1191.8271272,1214.1140272,1240.1135972,1269.8249872,1303.2472472,&
            1340.3792772,1381.2198972,1425.7677672,1474.0214672,1525.9794272,&
            1581.6399572,1641.0012772,1704.0614572,1770.8184572,1841.2701172,&
            1915.4141672,1993.2481972,2074.7696972,2159.9760272,2248.8644372,&
            2341.4320372,2437.6758472,2537.5927472,2641.1794972,2748.4327572,&
            2859.3490472,2973.9247772,3092.1562372,3214.0396072,3339.5709272,&
            3468.7461472,3601.5610672,3738.0113972,3878.0927172,4021.8004772,&
            4169.1300372,4320.0766072,4474.6352972,4632.8010972,4794.5688872,&
            4959.9334172,5128.8893172,5301.4311072,5477.5532072,5657.2498772,&
            5840.5152972/)
     energy(:,2)=(/ 3388.157917, 3391.841777, 3399.209367, 3410.260457, 3424.994687,&
          3443.411567, 3465.510497, 3491.290767, 3520.751517, 3553.891807,&
          3590.710537, 3631.206507, 3675.378397, 3723.224777, 3774.744057,&
          3829.934587, 3888.794537, 3951.321997, 4017.514927, 4087.371167,&
          4160.888417, 4238.064297, 4318.896267, 4403.381697, 4491.517827,&
          4583.301767, 4678.730527, 4777.800977, 4880.509887, 4986.853907,&
          5096.829537, 5210.433197, 5327.661177, 5448.509617, 5572.974587,&
          5701.052017, 5832.737697, 5968.027327, 6106.916477, 6249.400607,&
          6395.475047, 6545.135007, 6698.375597, 6855.191797, 7015.578457,&
          7179.530327, 7347.042037, 7518.108097, 7692.722887, 7870.880697,&
          8052.575667 /)
     nuclear_spin_weight(1:LENGTH_J_DATA:2)=1
     nuclear_spin_weight(2:LENGTH_J_DATA:2)=3
  end if
  !! Calculates partition function.
  partition = 0.0d0
  do i=1,LENGTH_J_DATA
     do jtmp=1,LENGTH_V_DATA
        J = dble(i-1)
        v = dble(jtmp-1)
        partition = partition + nuclear_spin_weight(i)*(2*J+1)*dexp(-energy(i,jtmp)*constant)
     end do
  end do
  !! Return population of desired level.
  boltzmann_factor = nuclear_spin_weight(Jwant+1)*(2*dble(Jwant)+1) &
       *dexp(-energy(Jwant+1,vwant+1)*constant)/partition

end subroutine calculate_boltzmann_factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine read_data_skip_comments(fileid,readstat,ww,xs,vpp,Jpp,Jp)
  implicit none
  real*8 :: ww,xs
  integer :: vpp,Jpp,Jp,fileid,readstat
  character*500 :: raw,line
  read(fileid,'(a)',iostat=readstat) raw
  !! return if error or end of file
  if (readstat.ne.0) return
  !! get rid of leading space
  line=adjustl(raw)
  !! skip is first character # or if line is blank
  if ((len(trim(line)).eq.0).or.(line(1:1).eq.'#')) then
     call read_data_skip_comments(fileid,readstat,ww,xs,vpp,Jpp,Jp)
  else
     ! read(line,*,iostat=readstat) ww,xs,vpp,Jpp,Jp
     read(line,*,iostat=readstat) ww,xs,Jpp,Jp
     vpp = 0
  end if
end subroutine read_data_skip_comments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

