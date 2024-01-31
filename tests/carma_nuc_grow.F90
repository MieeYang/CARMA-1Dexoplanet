!! This code is to test the aerosol freezing routin !
!! Upon execution, a text file (carma_nuctest.txt) is generated.
!! The text file can be read with the IDL procedure read_nuctest.pro.
!!
!! @author  Chuck Bardeen
!! @version July-2009

program carma_nuc_grow_test
  implicit none

  write(*,*) "1-D simulation:"

  call test_nuc_ttl()  
  
  write(*,*) "Done"

end program

!! An initial concentration of sulfate drops at the smallest size, then allow
!! that to nucleate ice and then the ice can grow using a gas. The total mass
!! of ice + gas should be conserved.

subroutine test_nuc_ttl()
  use carma_precision_mod 
  use carma_constants_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmaelement_mod
  use carmagroup_mod
  use carmagas_mod
  use carmasolute_mod
  use carmastate_mod
  use carma_mod
  use clima_mod
  implicit none
  integer                   :: if_res = 0       ! 1 for restart, 0 for initial run

!  integer, parameter        :: NZ           = 1
  integer, parameter        :: NX           = 1
  integer, parameter        :: NY           = 1
  integer, parameter        :: NZ           = 100
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 6
  integer, parameter        :: NBIN         = 50
  integer, parameter        :: NGROUP       = 4
  integer, parameter        :: NSOLUTE      = 1
  integer, parameter        :: NGAS         = 1
  integer, parameter        :: NWAVE        = 50
  integer, parameter        :: LUNOPRT      = 6
  integer, parameter        :: nstep        = 10

  integer, parameter        :: nskip  = 5

  real(kind=f), parameter   :: dtime_old  = 600._f

  real(kind=f), parameter   :: dtime  = 300._f
!  real(kind=f), parameter   :: dtime  = 1._f
!  real(kind=f), parameter   :: dtime  = 5._f
!  real(kind=f), parameter   :: dtime  = 10._f
!  real(kind=f), parameter   :: dtime  = 20._f
!  real(kind=f), parameter   :: dtime  = 100._f
!  real(kind=f), parameter   :: dtime  = 1000._f
!  real(kind=f), parameter   :: dtime  = 1800._f
!  real(kind=f), parameter   :: dtime  = 5000._f
!  real(kind=f), parameter   :: dtime  = 10000._f
!  real(kind=f), parameter   :: dtime  = 50000._f
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 200._f 
!!  real(kind=f), parameter   :: zmin   = 3000._f
  real(kind=f), parameter   :: zmin   = 0._f

!! Used for defining the sulfur concentration
  real(kind=f), parameter             :: n_IN    = 1.e-2_f     !! concentration (cm-3)
  real(kind=f), parameter             :: n_CCN    = 100._f     !! concentration (cm-3)
 
  real(kind=f), parameter             :: r0_IN   = 1e-4_f   !! mean radius (cm)
  real(kind=f), parameter             :: r0_CCN   = 1e-5_f   !! mean radius (cm)

  real(kind=f), parameter             :: fluxbot = 0 !1e4_f
  real(kind=f), parameter             :: concbot = 0 !1e4_f
!!  real(kind=f), parameter             :: ntop    = 1.5e1_f     !! concentration (cm-3) 
!!  real(kind=f), parameter             :: r0   = 1e-5_f   !! mean radius (cm)
!!  real(kind=f), parameter             :: rsig = 1.0_f      !! distribution width
  real(kind=f)                        :: RHO_CN = 2.00_f  

  real(kind=f)                        :: rhop
  real(kind=f)                        :: rhoa
  real(kind=f)                        :: noUse
!! This part for sources

  integer, parameter                  :: I_H2SO4   = 1               !! sulfate aerosol composition
  integer, parameter                  :: I_ICE     = 2               !! ice
  integer, parameter                  :: I_WATER   = 3               !! water
  integer, parameter                  :: I_DUST    = 4               !! dust

  type(carma_type), target            :: carma
  type(carma_type), pointer           :: carma_ptr
  type(carmastate_type)               :: cstate
  integer                             :: rc = 0
  
  real(kind=f), allocatable   :: xc(:,:,:)
  real(kind=f), allocatable   :: dx(:,:,:)
  real(kind=f), allocatable   :: yc(:,:,:)
  real(kind=f), allocatable   :: dy(:,:,:)
  real(kind=f), allocatable   :: zc(:,:,:)
  real(kind=f), allocatable   :: zl(:,:,:)
  real(kind=f), allocatable   :: p(:,:,:)
  real(kind=f), allocatable   :: pl(:,:,:)
  real(kind=f), allocatable   :: t(:,:,:)
  real(kind=f), allocatable   :: wind(:,:,:)
  real(kind=f), allocatable   :: ekz(:,:,:)
  real(kind=f), allocatable   :: relhum(:,:,:)
  real(kind=f), allocatable   :: wml(:)
  real(kind=f), allocatable   :: gr(:)
  real(kind=f), allocatable   :: tinput(:)

  real(kind=f), allocatable   :: mmr(:,:,:,:,:)
  real(kind=f), allocatable   :: mmr_gas(:,:,:,:)
  real(kind=f), allocatable   :: mmr_const(:,:,:,:,:)
  real(kind=f), allocatable   :: satliq(:,:,:,:)
  real(kind=f), allocatable   :: satice(:,:,:,:)

  real(kind=f), allocatable   :: ftopp(:,:)
  real(kind=f), allocatable   :: pcbot(:,:)
  real(kind=f), allocatable   :: fbotp(:,:)
  real(kind=f), allocatable   :: ftopg(:)
  real(kind=f), allocatable   :: fbotg(:)
  real(kind=f), allocatable   :: gcbot(:)
  real(kind=f), allocatable   :: gctop(:)  
  real(kind=f), allocatable   :: pctop(:,:)

  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: dr(:)
  real(kind=f), allocatable   :: rmass(:)

  real(kind=f), allocatable   :: lat(:,:)
  real(kind=f), allocatable   :: lon(:,:)
  
  integer               :: i
  integer               :: ix
  integer               :: iy
  integer               :: iz
  integer               :: istep
  integer               :: istep_old = 1
  integer               :: igas
  integer               :: igroup
  integer               :: ielem
  integer               :: ibin
  integer               :: aeroBin_IN
  integer               :: aeroBin_CCN
  integer               :: isource
  integer               :: ithread
  integer               :: iwave
  integer, parameter    :: lun = 42
  integer, parameter    :: lun2 = 43
  integer, parameter    :: aer = 45
  integer, parameter    :: gasini = 46
  integer, parameter    :: lundiag = 47

  integer               :: ITPP = -1

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat
  real(kind=f)          :: drh
  real(kind=f)          :: prod35
  logical               :: do_explised = .true. 
  real(kind=f)          :: ttt,psati, psatl
 
  real(kind=f)          :: wave(NWAVE)
  real(kind=f)          :: dwave(NWAVE)
  complex(kind=f)       :: refidx_ice(NWAVE)
  complex(kind=f)       :: refidx_water(NWAVE)
  real(kind=f)          :: r_refidx_ice(NWAVE)
  real(kind=f)          :: i_refidx_ice(NWAVE)
  real(kind=f)          :: r_refidx_water(NWAVE)
  real(kind=f)          :: i_refidx_water(NWAVE)
  real(kind=f), allocatable     :: qext(:,:)
  real(kind=f), allocatable     :: ssa(:,:)
  real(kind=f), allocatable     :: asym(:,:)

  ! Set the wavelengths
  data wave &
  /  0.100_f,  0.110_f,  0.130_f,  0.150_f,  0.170_f,  0.190_f, &
     0.340_f,  0.380_f,  0.412_f,  0.440_f,  0.443_f,  0.490_f, &
     0.500_f,  0.531_f,  0.532_f,  0.551_f,  0.555_f,  0.667_f, &
     0.675_f,  0.870_f,  1.020_f,  1.111_f,  1.333_f,  1.562_f, &
     1.640_f,  1.770_f,  2.051_f,  2.210_f,  2.584_f,  3.284_f, &
     3.809_f,  4.292_f,  4.546_f,  4.878_f,  5.128_f,  5.405_f, &
     5.714_f,  6.061_f,  6.452_f,  6.897_f,  7.407_f,  8.333_f, &
     9.009_f, 10.309_f, 12.500_f, 13.889_f, 16.667_f, 20.000_f, &
    26.316_f, 35.714_f /  
  
  data r_refidx_ice &
  / 1.3980_f, 1.3886_f, 1.3692_f, 1.6177_f, 1.4807_f, 1.4101_f, &
    1.3267_f, 1.3215_f, 1.3183_f, 1.3163_f, 1.3161_f, 1.3135_f, &
    1.3130_f, 1.3117_f, 1.3116_f, 1.3110_f, 1.3108_f, 1.3077_f, &
    1.3075_f, 1.3037_f, 1.3012_f, 1.2997_f, 1.2954_f, 1.2903_f, &
    1.2882_f, 1.2839_f, 1.2723_f, 1.2617_f, 1.2087_f, 1.5860_f, &
    1.3803_f, 1.3427_f, 1.3472_f, 1.3392_f, 1.3252_f, 1.3097_f, &
    1.2933_f, 1.3015_f, 1.3201_f, 1.3242_f, 1.3209_f, 1.2964_f, &
    1.2695_f, 1.1502_f, 1.3822_f, 1.5595_f, 1.5294_f, 1.4986_f, &
    1.3799_f, 1.2441_f /    
      
  data i_refidx_ice &
  /   3.1273e-1_f,    2.44e-1_f,    1.73e-1_f,  2.5867e-1_f,  3.3495e-4_f,  1.4156e-8_f, &
        2.0e-11_f,    2.0e-11_f, 2.7622e-11_f, 6.2680e-11_f, 7.1593e-11_f,  4.172e-10_f, &
     5.8890e-10_f,  1.4494e-9_f,  1.4898e-9_f,   2.344e-9_f,   2.564e-9_f,   1.821e-8_f, &
        1.99e-8_f,    2.65e-7_f,    2.25e-6_f,   1.766e-6_f,   1.326e-5_f,  3.8874e-4_f, &
      2.4415e-4_f,  1.4943e-4_f,   1.367e-3_f,  2.3877e-4_f,    7.41e-4_f,    0.11671_f, &
      7.5714e-3_f,  2.1315e-2_f,  3.0921e-2_f,  1.4303e-2_f,   1.243e-2_f,  1.7006e-2_f, &
       3.256e-2_f,   6.439e-2_f,   5.964e-2_f,   5.428e-2_f,  4.5215e-2_f,      0.037_f, &
      3.6847e-2_f,  7.4906e-2_f,      0.422_f,    0.29411_f,    0.12302_f,      0.067_f, &
      3.6117e-2_f,    0.19131_f /

  data r_refidx_water &
  / 1.4766_f, 1.5292_f, 1.6338_f, 1.5219_f, 1.6028_f, 1.4755_f, &
    1.3604_f, 1.3528_f, 1.3482_f, 1.3449_f, 1.3446_f, 1.3402_f, &
    1.3394_f, 1.3372_f, 1.3371_f, 1.3359_f, 1.3356_f, 1.3300_f, &
    1.3297_f, 1.3243_f, 1.3213_f, 1.3196_f, 1.3154_f, 1.3106_f, &
    1.3086_f, 1.3044_f, 1.2945_f, 1.2849_f, 1.2349_f, 1.4389_f, &
    1.3468_f, 1.3201_f, 1.3126_f, 1.3092_f, 1.3010_f, 1.2850_f, &
    1.2547_f, 1.2733_f, 1.3208_f, 1.2986_f, 1.2836_f, 1.2589_f, &
    1.2374_f, 1.1738_f, 1.0979_f, 1.1778_f, 1.3366_f, 1.4677_f, &
    1.5449_f, 1.4945_f /

  data i_refidx_water &
  / 5.4299e-1_f,  4.2956e-1_f,  2.3961e-1_f,   2.772e-1_f,  3.8493e-3_f,   3.807e-7_f, &
     2.758e-9_f,  1.9433e-9_f,  1.3904e-9_f, 9.3290e-10_f, 8.9121e-10_f, 7.3614e-10_f, &
    9.232e-10_f,  1.7925e-9_f,  1.8191e-9_f,  2.5046e-9_f,  2.6748e-9_f,  2.0563e-8_f, &
     2.186e-8_f,  3.7155e-7_f,  2.3793e-6_f,  2.1433e-6_f,  2.4157e-5_f,  1.2171e-4_f, &
    7.9131e-5_f,  1.1193e-4_f,  6.7532e-4_f,  3.3914e-4_f,   2.398e-3_f,  4.5897e-2_f, &
     2.423e-3_f,  8.3906e-3_f,  1.3517e-2_f,  1.4016e-2_f,  1.0869e-2_f,  1.0395e-2_f, &
    2.1369e-2_f,    0.12603_f,     4.2e-2_f,   3.222e-2_f,  3.2434e-2_f,  3.5727e-2_f, &
    3.9987e-2_f,  5.8613e-2_f,      0.259_f,    0.36416_f,    0.42876_f,    0.39334_f, &
       0.3187_f,    0.29709_f /

  ! Construct a complex refractive index.
  refidx_ice   = cmplx(r_refidx_ice, i_refidx_ice)
  refidx_water = cmplx(r_refidx_water, i_refidx_water)

  ! Open the output text file
  open(unit=lun,file="carma_nuc_grow.txt",status="unknown",position="append")
  
  open(unit=lun2,file="carma_nuc_grow_mie.txt",status="unknown",position="append")
  open(unit=lundiag,file="diagnosis.txt",status="unknown",position="append")
 
  ! Convert wavelength um -> cm
  wave = wave * 1e-4_f
  ! Calcualte dwave according to wave
  do iwave = 2,NWAVE-1
    dwave(iwave) = (wave(iwave+1)-wave(iwave-1))/2.
  end do
  dwave(1) = wave (2)-wave(1)
  dwave(NWAVE) = wave(NWAVE)-wave(NWAVE-1)
  
  ! Allocate the arrays that we need for the model
  allocate(xc(NZ,NY,NX), dx(NZ,NY,NX), yc(NZ,NY,NX), dy(NZ,NY,NX), &
           zc(NZ,NY,NX), zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), &
           wind(NZ,NY,NX) , t(NZ,NY,NX), ekz(NZP1,NY,NX), wml(NZ), gr(NZ), &
           tinput(NZ)) 
  allocate(mmr(NZ,NY,NX,NELEM,NBIN))
  allocate(mmr_gas(NZ,NY,NX,NGAS))
  allocate(mmr_const(NZ,NY,NX,NELEM,NBIN))
  allocate(lat(NY,NX), lon(NY,NX))
  allocate(satliq(NZ,NY,NX,NGAS))
  allocate(satice(NZ,NY,NX,NGAS))
  allocate(r(NBIN))
  allocate(dr(NBIN))
  allocate(rmass(NBIN))
  allocate(qext(NWAVE,NBIN))
  allocate(ssa(NWAVE,NBIN)) 
  allocate(asym(NWAVE,NBIN))
  allocate(ftopg(NGAS), fbotg(NGAS), ftopp(NBIN,NELEM), pcbot(NBIN,NELEM), pctop(NBIN,NELEM),&
           fbotp(NBIN,NELEM),gctop(NGAS),gcbot(NGAS))

  ! Vertical profile

  do i = 1, NZ
    zc(i,:,:) = zmin + (deltaz * (i - 0.5_f))
  end do
  call GetClimaAtmosphere(zc(:,1,1), p=p(:,1,1), t=t(:,1,1))

  do i = 1, NZP1
    zl(i,:,:) = zmin + ((i - 1) * deltaz)
  end do
  call GetClimaAtmosphere(zl(:,1,1), p=pl(:,1,1))
  do iz = 1, 10
    ekz(iz,:,:) = 4.e6_f
    wind(iz,:,:) = 0._f
  end do

  do iz = 11, NZ
    ekz(iz,:,:) = 1.e5_f
    wind(iz,:,:) = 0._f    
  end do

  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, LUNOPRT=LUNOPRT, wave=wave, dwave=dwave)
  if (rc /=0) stop "    *** FAILED ***"
	carma_ptr => carma

  ! Define the groups
  call CARMAGROUP_Create(carma, 1, "IN", 1.e-6_f, 2.0_f, I_SPHERE, 1._f, .false., &
                         rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                         scavcoef=0.1_f, do_mie=.false., shortname="IN", &
                         ifallrtn=I_FALLRTN_STD_SHAPE)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAGROUP_Create(carma, 2, "CCN", 1.e-6_f, 2.0_f, I_SPHERE, 1._f, .false., &
                         rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                         scavcoef=0.1_f, do_mie=.false., shortname="CCN", &
                         ifallrtn=I_FALLRTN_STD_SHAPE)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAGROUP_Create(carma, 3, "Ice Crystal", ((RHO_CN/RHO_I)**(1._f/3._f))*2.e-6_f, 2.0_f, I_SPHERE, 3._f, .true., &
                         rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                         scavcoef=0.1_f, refidx=refidx_ice, do_mie=.true., shortname="CRICE", &
                         ifallrtn=I_FALLRTN_STD_SHAPE)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAGROUP_Create(carma, 4, "Water Drop", ((RHO_CN/RHO_W)**(1._f/3._f))*2.e-6_f, 2.0_f, I_SPHERE, 1._f, .false., &
                         rc, do_wetdep=.true., do_drydep=.true., solfac=0.3_f, &
                         scavcoef=0.1_f, refidx=refidx_water, do_mie=.true., shortname="CRDUST", &
                         ifallrtn=I_FALLRTN_STD_SHAPE)
  if (rc /=0) stop "    *** FAILED ***"

  ! Define the elements
  call CARMAELEMENT_Create(carma, 1, 1, "IN", RHO_CN, I_INVOLATILE, I_H2SO4, rc, shortname="CRIN", isolute=1)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAELEMENT_Create(carma, 2, 2, "CCN", RHO_CN, I_INVOLATILE, I_H2SO4, rc, shortname="CRCCN", isolute=1)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAELEMENT_Create(carma, 3, 3, "Ice Crystal", RHO_I, I_VOLATILE, I_ICE, rc, shortname="CRICE")
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAELEMENT_Create(carma, 4, 3, "Core Mass", RHO_CN, I_COREMASS, I_H2SO4, rc, shortname="CRCORE", isolute=1)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAELEMENT_Create(carma, 5, 4, "Liquid Water", RHO_W, I_VOLATILE, I_WATER, rc, shortname="CRLW")
  if (rc /=0) stop "    *** FAILED ***"

  call CARMAELEMENT_Create(carma, 6, 4, "Liquid Core Mass", RHO_CN, I_COREMASS, I_H2SO4, rc, shortname="LCRCORES", isolute=1)
  if (rc /=0) stop "    *** FAILED ***"


  ! Define the Solutes
  call CARMASOLUTE_Create(carma, 1, "Dust", 2, 98._f, RHO_CN, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Define the gases
  call CARMAGAS_Create(carma, 1, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, I_GCOMP_H2O, rc, shortname='Q',ds_threshold=-0.0001_f)
  if (rc /=0) stop "    *** FAILED ***"

  ! Setup the CARMA processes to exercise growth and nucleation.
  call CARMA_AddGrowth(carma, 3, 1, rc)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMA_AddGrowth(carma, 5, 1, rc)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMA_AddNucleation(carma, 1, 4, I_HETNUC, 0._f, rc, igas=1, ievp2elem=1, mucos = 0.945_f)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMA_AddNucleation(carma, 2, 6, I_DROPACT, 0._f, rc, igas=1, ievp2elem=1, mucos = 0.945_f)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMA_AddNucleation(carma, 3, 5, I_ICEMELT, 0._f, rc, igas=1, ievp2elem=1)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMA_AddNucleation(carma, 5, 3, I_DROPFREEZE, 0._f, rc, igas=1, ievp2elem=1)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMA_AddCoagulation(carma, 3, 3, 3, I_COLLEC_DATA, rc)
  if (rc /=0) stop "    *** FAILED ***"

  call CARMA_AddCoagulation(carma, 4, 4, 4, I_COLLEC_DATA, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Setup processes of simulation
!  call CARMA_Initialize(carma, rc, do_coag=.true., do_detrain=.true., do_grow=.true., do_explised=.false., do_vdiff=.true., do_vtran=.true., do_drydep=.false., do_substep=.true.,maxsubsteps = 1, minsubsteps=1,maxretries = 20, itbnd_pc =I_FLUX_SPEC, ibbnd_pc=I_FLUX_SPEC, itbnd_gc=I_FLUX_SPEC, ibbnd_gc = I_FIXED_CONC)

  call CARMA_Initialize(carma, rc, do_coag=.true., do_detrain=.true., &
do_grow=.true., do_explised=do_explised, do_vdiff=.true., do_vtran=.true., &
do_drydep=.false., do_substep=.true.,maxsubsteps = 1000, minsubsteps=1, maxretries = 20, itbnd_pc =I_FIXED_CONC, itbnd_gc=I_FLUX_SPEC, ibbnd_gc = I_FIXED_CONC, ibbnd_pc = I_FIXED_CONC)
  if (rc /=0) stop "    *** FAILED ***"

  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom 
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.
  lat(:,:) = -40.0_f
  lon(:,:) = -105.0_f
  wml(:) = WTMOl_AIR
  gr(:) = GRAV
  ! Horizonal centers
  do ix = 1, NX
    do iy = 1, NY
      dx(:,iy,ix) = deltax
      xc(:,iy,ix) = ix*dx(:,iy,ix) / 2._f
      dy(:,iy,ix) = deltay
      yc(:,iy,ix) = iy*dy(:,iy,ix) / 2._f
    end do
  end do
  

  if (if_res .eq. 0) then

    ! Write output for the test
    write(lun,*) NGROUP, NELEM, NBIN, NGAS, NZ

    do iz = 1, NZ
      write(lun,'(i4,3e10.3)') iz, zc(iz,1,1), p(iz,1,1), t(iz,1,1)
    end do

    do igroup = 1, NGROUP
      call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, rmass=rmass,qext=qext, ssa=ssa, asym=asym)
      if (rc /=0) stop "    *** FAILED ***"
      do ibin = 1, NBIN
        write(lun,'(2i4,2e10.3)') igroup, ibin, r(ibin) * 1e4_f, rmass(ibin)
      end do
    end do

    !!!!!!!!!!!!! Write initial conditions
    call CARMAGROUP_Get(carma, 1, rc, r=r, dr=dr, rmass=rmass,qext=qext,ssa=ssa, asym=asym)
    mmr(:,:,:,:,:) = 0.
    mmr_gas(:,:,:,:) = 0.

    aeroBin_IN = int(log10(r0_IN / (1e-6_f)) / log10(2._f) * 3. + 1.)
    aeroBin_CCN = int(log10(r0_CCN / (1e-6_f)) / log10(2._f) * 3. + 1.)

    do iz = 1,NZ
      ttt = t(iz,1,1)

      psati = 10.0_f * exp(9.550426_f - (5723.265_f / ttt) + (3.53068_f * log(ttt)) - (0.00728332_f * ttt))
      psatl = 10.0_f * exp(54.842763_f - (6763.22_f / ttt) - (4.210_f * log(ttt))+ (0.000367_f * ttt) + &
             (tanh(0.0415_f * (ttt - 218.8_f)) * (53.878_f - (1331.22_f / ttt)-(9.44523_f * log(ttt)) + 0.014025_f* ttt)))

      if (ITPP.eq.0) then
        mmr_gas(iz,:,:,1) = 0.  
      else
        mmr_gas(iz,:,:,1) = 0.95 * min(1.0*(WTMOL_H2O/WTMOL_AIR)*10.*psati/(R_AIR)/t(iz,1,1)/(abs(pl(iz+1,1,1)-pl(iz,1,1))/GRAV/deltaz), &
                      1.0*(WTMOL_H2O/WTMOL_AIR)*10.*psatl/(R_AIR)/t(iz,1,1)/(abs(pl(iz+1,1,1)-pl(iz,1,1))/GRAV/deltaz))

        if (iz.gt.1) then
          if (mmr_gas(iz,1,1,1) .gt. mmr_gas(iz-1,1,1,1)) then   
            mmr_gas(iz,:,:,1) = 0
            ITPP = 0
          end if
        end if
      end if
!      mmr_const(iz,1,1,1,aeroBin_IN) = n_IN * rmass(aeroBin_IN) /  ((abs(pl(iz+1,1,1)-pl(iz,1,1))/GRAV/deltaz)) * 10.
!      mmr_const(iz,1,1,2,aeroBin_CCN) = n_CCN * rmass(aeroBin_CCN) /  ((abs(pl(iz+1,1,1)-pl(iz,1,1))/GRAV/deltaz)) * 10.
    end do

!    mmr(:,1,1,1,:) =  mmr_const(:,1,1,1,:)
!    mmr(:,1,1,2,:) =  mmr_const(:,1,1,2,:)

    ttt = t(1,1,1)
!    psati = 10.0_f * exp(9.550426_f - (5723.265_f / ttt) + (3.53068_f * log(ttt)) - (0.00728332_f * ttt))
    psatl = 10.0_f * exp(54.842763_f - (6763.22_f / ttt) - (4.210_f * log(ttt))+ (0.000367_f * ttt) + &
             (tanh(0.0415_f * (ttt - 218.8_f)) * (53.878_f - (1331.22_f / ttt)-(9.44523_f * log(ttt)) + 0.014025_f* ttt)))
    gcbot(1) = 0.95*(WTMOL_H2O/WTMOL_AIR)*10.*psatl/(R_AIR)/t(1,1,1)/(abs(pl(1+1,1,1)-pl(1,1,1))/GRAV/deltaz)

    write(lun,*) 0

    do iz = 1, NZ
      do ielem = 1, NELEM
        do ibin = 1, NBIN
          write(lun,'(3i4,e10.3)') iz, ielem, ibin, real(mmr(iz,1,1,ielem,ibin))
        end do
      end do
    end do

    do iz = 1, NZ
      do igas = 1, NGAS
        write(lun,'(2i4,3e10.3)') iz, igas, real(mmr_gas(iz,1,1,igas)),0.,0.
      end do
    end do

    write(lun2,*) NGROUP, NWAVE, NBIN
  
    do iwave = 1, NWAVE
      write(lun2,'(i3,e10.3)') iwave, wave(iwave)
    end do
    do igroup = 1, NGROUP

      call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, qext=qext, ssa=ssa, asym=asym)
      do ibin = 1, NBIN
        write(lun2,'(i3,e10.3)') ibin, r(ibin)
      end do

      do iwave = 1, NWAVE
        write(lun2,'(i3,2e10.3,2e10.3)') iwave, refidx_ice(iwave), refidx_water(iwave)
      end do
    end do

    do ibin = 1, NBIN
      pcbot(ibin,:) = 0._f
      pctop(ibin,:) = 0._f
    end do
    pcbot(aeroBin_IN,1) = n_IN
    pctop(aeroBin_IN,1) = n_IN
    pcbot(aeroBin_CCN,2) = n_CCN
    pctop(aeroBin_CCN,2) = 0

    ftopg(1) = 0._f

  else
    aeroBin_IN = int(log10(r0_IN / (1e-6_f)) / log10(2._f) * 3. + 1.)
    aeroBin_CCN = int(log10(r0_CCN / (1e-6_f)) / log10(2._f) * 3. + 1.)
    ttt = t(1,1,1)
!    psat = 10.0_f * exp(9.550426_f - (5723.265_f / ttt) + (3.53068_f * log(ttt)) - (0.00728332_f * ttt))
    psatl = 10.0_f * exp(54.842763_f - (6763.22_f / ttt) - (4.210_f * log(ttt))+ (0.000367_f * ttt) + &
             (tanh(0.0415_f * (ttt - 218.8_f)) * (53.878_f - (1331.22_f / ttt)-(9.44523_f * log(ttt)) + 0.014025_f* ttt)))
    gcbot(1) = 0.95*(WTMOL_H2O/WTMOL_AIR)*10.*psatl/(R_AIR)/t(1,1,1)/(abs(pl(1+1,1,1)-pl(1,1,1))/GRAV/deltaz)

    do ibin = 1, NBIN
      pcbot(ibin,:) = 0._f
      pctop(ibin,:) = 0._f
    end do
    pcbot(aeroBin_IN,1) = n_IN
    pctop(aeroBin_IN,1) = n_IN
    pcbot(aeroBin_CCN,2) = n_CCN
    pctop(aeroBin_CCN,2) = 0

    ftopg(1) = 0._f

    call read_restart(mmr, mmr_gas,satice,satliq,NELEM, NGAS, NBIN, NZ, istep_old)
  end if

  write(*,*) istep_old

  ! Iterate the model over a few time steps.
  do istep = istep_old, istep_old+nstep-1
    write(*,*) istep

    ! Calculate the model time.
    time = (istep - 1) * dtime

    if (mod(istep, nskip) .eq. 0) then
      write(lundiag,'(f12.0)') istep*dtime
      carma%f_lundiag = lundiag
      carma%f_do_printdiag = .TRUE.
    else
      carma%f_do_printdiag = .FALSE.
    end if

!    mmr(:,1,1,1,:) =  mmr_const(:,1,1,1,:)
    !1D
    ix = 1
    iy = 1
    ! Create a CARMASTATE for this column.
 !   tinput = max(t(:,iy,ix)-2., t(:,iy,ix)-2.*(8.*float(istep)/float(nstep)))
    tinput = t(:,iy,ix)
    call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
                        I_CART, I_CART, lat(iy,ix), lon(iy,ix), &
                        xc(:,iy,ix), dx(:,iy,ix), &
                        yc(:,iy,ix), dy(:,iy,ix), &
                        zc(:,iy,ix), zl(:,iy,ix), p(:,iy,ix), &
                        pl(:,iy,ix), tinput(:), wml(:), gr(:), RPLANET, rc, told = tinput, &
                        winds = wind(:,iy,ix), ekz = ekz(:,iy,ix), &
                        pctop=pctop, ftopg=ftopg, pcbot=pcbot, gcbot = gcbot)
    if (rc /=0) stop "    *** FAILED ***"
  
    ! Send the bin mmrs to CARMA
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_SetBin(cstate, ielem, ibin, &
                               mmr(:,iy,ix,ielem,ibin), rc)
        if (rc /=0) stop "    *** FAILED ***"
      end do
    end do
    
    ! Send the gas mmrs to CARMA
    do igas = 1, NGAS
      call CARMASTATE_SetGas(cstate, igas, mmr_gas(:,iy,ix,igas), rc, &
                             mmr_old = mmr_gas, satice_old = satice, &
                             satliq_old = satliq)
      if (rc /=0) stop "    *** FAILED ***"
    end do

    ! Execute the step
    call CARMASTATE_Step(cstate, rc, lndfv=2._f, ocnfv=2._f, icefv=2._f, lndram=100._f, ocnram=100._f, iceram=100._f, lndfrac=0._f, ocnfrac=1._f, icefrac=0._f)
    if (rc /=0) stop "    *** FAILED ***"

    ! Get the updated bin mmr.
    do ielem = 1, NELEM
      do ibin = 1, NBIN
        call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,iy,ix,ielem,ibin), rc)
        if (rc /=0) stop "    *** FAILED ***"
      end do
    end do

    ! Get the updated gas mmr.
    do igas = 1, NGAS
      call CARMASTATE_GetGas(cstate, igas, &
                              mmr_gas(:,iy,ix,igas), rc, &
                              satliq=satliq(:,iy,ix,igas), &
                              satice=satice(:,iy,ix,igas))
      if (rc /=0) stop "    *** FAILED ***"
    end do
    ! Write output for test.
    if (mod(istep, nskip) .eq. 0) then
      write(lun,'(f12.0)') istep*dtime
      write(lun2,'(f12.0)') istep*dtime

      do iz = 1, NZ
        do ielem = 1, NELEM
          do ibin = 1, NBIN
            write(lun,'(3i4,e10.3,e10.3)') iz, ielem, ibin, real(mmr(iz,1,1,ielem,ibin)), real(mmr(iz,1,1,ielem,ibin)*p(iz,1,1) / 287._f / t(iz,1,1))
          end do
        end do
      end do

      do iz = 1, NZ
        do igas = 1, NGAS
          write(lun,'(2i4,3e12.3)') iz, igas, real(mmr_gas(iz,1,1,igas)), satliq(iz,1,1,igas), satice(iz,1,1,igas)
        end do
      end do
    !write mie output  
    do igroup = 1, NGROUP
      call CARMAGROUP_Get(carma, igroup, rc, r=r, dr=dr, qext=qext, ssa=ssa, asym=asym)
      do iwave = 1, NWAVE
        do ibin = 1, NBIN
          write(lun2,'(2i3,3(x,e10.3))') iwave, ibin, qext(iwave,ibin), ssa(iwave,ibin), asym(iwave,ibin)
        end do
      end do
    end do
  
    end if
  end do   ! time loop

  call write_restart(cstate, NELEM, NGAS, NBIN, NZ, istep)
  call write_smart(cstate, carma, wave, NGROUP, NELEM, NGAS, NBIN,NWAVE,NZ,0)

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** FAILED ***"

  ! Close the output file
  close(unit=lun)	
  close(unit=lun2)	

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** FAILED ***"
end subroutine

subroutine read_restart(mmr, mmr_gas,satice,satliq, NELEM, NGAS, NBIN, NZ, istep)
  use carmastate_mod
  use carma_types_mod
  integer, intent(out)                  :: istep
  integer, parameter                    :: res = 44
  real(kind=f), intent(out)                     :: mmr(NZ,1,1,NELEM,NBIN)
  real(kind=f), intent(out)                     :: mmr_gas(NZ,1,1,NGAS)
  real(kind=f), intent(out)                     :: satice(NZ,1,1,NGAS)
  real(kind=f), intent(out)                     :: satliq(NZ,1,1,NGAS)

  open(unit=res, file="carma_restart.txt",status="old")
  read(res,*) istep
  ix=1
  iy=1
  do igas = 1, NGAS
    do iz = 1, NZ
      read(res,*) mmr_gas(iz,iy,ix,igas),satice(iz,iy,ix,igas),satliq(iz,iy,ix,igas)
    end do
  end do

  do ielem = 1, NELEM
    do ibin = 1, NBIN
      do iz = 1, NZ
        read(res,*) mmr(iz,iy,ix,ielem,ibin)
      end do
    end do
  end do

  close(unit=res)
  return
end subroutine

subroutine write_restart(cstate, NELEM, NGAS, NBIN, NZ, istep)
  use carmastate_mod
  use carma_types_mod
  integer                     :: iz
  integer                     :: igas
  integer                     :: igroup
  integer                     :: ielem
  integer                     :: ibin
  integer, parameter          :: res = 44
  real(kind=f), allocatable   :: mmr_gas(:,:,:,:)
  real(kind=f), allocatable   :: satliq(:,:,:,:)
  real(kind=f), allocatable   :: satice(:,:,:,:)
  real(kind=f), allocatable   :: mmr(:,:,:,:,:)
  integer                     :: rc = 0
  type(carmastate_type), intent(in)    :: cstate      !! the carma state

  allocate(mmr(NZ,1,1,NELEM,NBIN))
  allocate(mmr_gas(NZ,1,1,NGAS))
  allocate(satliq(NZ,1,1,NGAS))
  allocate(satice(NZ,1,1,NGAS))

  open(unit=res, file="carma_restart.txt",status="unknown")

  write(res,*) istep

  ix = 1
  iy = 1
  do igas = 1, NGAS
    call CARMASTATE_GetGas(cstate, igas, mmr_gas(:,iy,ix,igas), rc, &
                        satliq=satliq(:,iy,ix,igas), &
                        satice=satice(:,iy,ix,igas))
    if (rc /=0) stop "    *** FAILED ***"
    do iz = 1, NZ
      write(res,*) mmr_gas(iz,iy,ix,igas),satice(iz,iy,ix,igas),satliq(iz,iy,ix,igas)
    end do
  end do

  do ielem = 1, NELEM
    do ibin = 1, NBIN
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,iy,ix,ielem,ibin), rc)
      if (rc /=0) stop "    *** FAILED ***"
      do iz = 1, NZ
        write(res,*) mmr(iz,iy,ix,ielem,ibin)
      end do
    end do
  end do

  close(unit=res)
  return
end subroutine

subroutine write_smart(cstate, carma, wave, NGROUP, NELEM, NGAS, NBIN, NWAVE,NZ,dz)
  use carmastate_mod
  use carma_mod
  use carma_types_mod
  use carmagroup_mod
  integer                     :: iz
  integer                     :: igas
  integer                     :: igroup
  integer                     :: ielem
  integer                     :: ibin
  integer                     :: ifile

  integer, parameter          :: smart1 = 45
  integer, parameter          :: opfile = 46

  character( len = 3 )        :: cfile

  real(kind=f)                :: qsca
  real(kind=f), allocatable   :: temp(:)
  real(kind=f), allocatable   :: pres(:)
  real(kind=f), allocatable   :: rhoa(:)
  real(kind=f), allocatable   :: qext(:,:,:)
  real(kind=f), allocatable   :: ssa(:,:,:)
  real(kind=f), allocatable   :: asym(:,:,:)
  real(kind=f), allocatable   :: r(:,:)
  real(kind=f), allocatable   :: dr(:,:)
  real(kind=f), allocatable   :: rmass(:,:)

  real(kind=f), allocatable   :: mmr_gas(:,:,:,:)
  real(kind=f), allocatable   :: satliq(:,:,:,:)
  real(kind=f), allocatable   :: satice(:,:,:,:)
  real(kind=f), allocatable   :: mmr(:,:,:,:,:)
  integer                     :: rc = 0
  type(carmastate_type), intent(in)    :: cstate      !! the carma state
  type(carma_type), target, intent(in)            :: carma
  real(kind=f), intent(in)    :: wave(NWAVE)

  allocate(mmr(NZ,1,1,NELEM,NBIN))
  allocate(mmr_gas(NZ,1,1,NGAS))
  allocate(temp(NZ), pres(NZ),rhoa(NZ))
  allocate(qext(NGROUP,NWAVE,NBIN))
  allocate(ssa(NGROUP,NWAVE,NBIN))
  allocate(asym(NGROUP,NWAVE,NBIN))

  allocate(r(NGROUP,NBIN))
  allocate(dr(NGROUP,NBIN))
  allocate(rmass(NGROUP,NBIN))

  open(unit=smart1, file="smartColumn.txt",status="unknown")
  ix = 1
  iy = 1

! Read in temp, pressure
  call CARMASTATE_GetState(cstate, rc, temp, pres)
  if (rc /=0) stop "    *** FAILED ***"
! Read in Gas info 
  do igas = 1, NGAS
    call CARMASTATE_GetGas(cstate, igas,mmr_gas(:,iy,ix,igas), rc, &
                        satliq=satliq(:,iy,ix,igas), &
                        satice=satice(:,iy,ix,igas))
    if (rc /=0) stop "    *** FAILED ***"
  end do

! Read in aerosol info
  do ielem = 1, NELEM
    do ibin = 1, NBIN
      call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,iy,ix,ielem,ibin), rc)
      if (rc /=0) stop "    *** FAILED ***"
!      do iz = 1, NZ
!        write(res,*) mmr(iz,iy,ix,ielem,ibin)
!      end do
    end do
  end do

! Get mie info
  do igroup = 1,NGROUP
    call CARMAGROUP_Get(carma, igroup, rc, r=r(igroup,:), dr=dr(igroup,:), &
                        rmass=rmass(igroup,:), qext=qext(igroup,:,:), &
                        ssa=ssa(igroup,:,:), asym=asym(igroup,:,:))
    if (rc /=0) stop "    *** FAILED ***"
  end do

  do IZ = 1,NZ
    if (mod(iz, 1) .eq. 0) then
      write(smart1,'(4e10.3)') pres(NZ-IZ+1), temp(NZ-IZ+1), mmr_gas(NZ-iz+1,1,1,1),1-mmr_gas(NZ-iz+1,1,1,1)
    end if
  end do
  close(unit=smart1)

  ifile = 1
  do igroup=2,3
    do IBIN=1,NBIN
      write( cfile,'(i3)' ) ifile
      open(unit=OPFILE, file="carma_aer"//trim(adjustl( cfile ))//".mie",status="unknown")
      do IWAVE =1,NWAVE
        qsca = qext(igroup,IWAVE,ibin)*ssa(igroup,IWAVE,ibin)
        write(opfile,'(e10.3,3(x,e10.3))') wave(IWAVE)*1e4,qext(igroup,IWAVE,ibin), qsca, asym(igroup,IWAVE,ibin)
      end do
      close(unit=opfile)
      ifile = ifile+1
    end do
  end do

  do iz =1,NZ
    rhoa(iz) = pres(iz) *10. / (R_AIR) / temp(iz)
  end do
  ifile = 1
  do ielem=1,2
    do IBIN=1,NBIN
      write( cfile,'(i3)' ) ifile
      open(unit=opfile,file='cld_tau'//trim(adjustl( cfile))//'.txt',status="unknown")
      do IZ =1,NZ
        taucld = 200.*pi*((r(ielem+1,ibin)/100.)**2)*mmr(NZ-iz+1,iy,ix,ielem*2,ibin)*qext(ielem+1,26,ibin)*rhoa(NZ+1-iz)/(rmass(ielem+1,ibin)/1000.)
!        taucld = 250.*pi*((r(ielem+1,ibin)/100.)**2)*0.000001*qext(ielem+1,26,ibin)*(0.8+0.1*ielem)/(rmass(ielem+1,ibin)/1000.)
        if (mod(iz, 1) .eq. 0) then
          write(opfile,*) pres(NZ-IZ+1), 1.*taucld
        end if
      end do

      close(unit=opfile)
      ifile = ifile+1
    end do
  end do

  return
end subroutine
