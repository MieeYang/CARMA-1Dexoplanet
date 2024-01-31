!! Public domain code to calculate the US Standard Atmosphere pressure and temperature.
!! This underlying routine, Atmosphere was downloaded from the internet at:
!!
!!  http://www.pdas.com/programs/atmos.f90
!!
!! This module wraps the routine in an interface similar to the one used by the CARMA
!! module, which allows columns, vectors of columns, and arrays of columns via the
!! method GetStandardAtmosphere().
module clima_mod
  ! types
  use carma_precision_mod
  use carma_constants_mod
  use carma_planet_mod
  implicit none
  
  private

  interface GetClimaAtmosphere
    module procedure GetClimaAtmosphere_1D
  end interface
  
  public GetClimaAtmosphere

  contains

    subroutine GetClimaAtmosphere_1D(z, p, t)
      real(kind=f), intent(in)              :: z(:)     !! Geometric Altitude (m)
      real(kind=f), optional, intent(out)   :: p(:)     !! pressure (Pa)
      real(kind=f), optional, intent(out)   :: t(:)     !! temperature (K)

      ! Local variables
      integer                               :: NZ
      integer                               :: MNZ 
      integer                               :: iz
      integer                               :: il
      integer                               :: j
      real(kind=f)                          :: botfluxR, pz, tz, qz, wz, frac
      real(kind=f)                          :: pz0
      real(kind=f)                          :: gcmgrav = 9.8
      real(kind=f), allocatable             :: zdense(:)
      real(kind=f), allocatable             :: pdense(:)
      real(kind=f), allocatable             :: rhodense(:)
      real(kind=f), allocatable             :: dsrho(:)
      real(kind=f), allocatable             :: dsz(:)
      real(kind=f), allocatable             :: dsp(:)
      real(kind=f), allocatable             :: dst(:)
      real(kind=f), allocatable             :: dsw(:)
      real(kind=f), allocatable             :: Jin(:)
      real(kind=f), allocatable             :: Zin(:)
      real(kind=f), allocatable             :: Pin(:)
      real(kind=f), allocatable             :: Tin(:)
      real(kind=f), allocatable             :: Qin(:)
      real(kind=f), allocatable             :: Oin(:)
      real(kind=f), allocatable             :: Win(:)
      real(kind=f), allocatable             :: RHOin(:)     !! Density (kg/m3)
      integer, parameter                    :: clm = 45

      NZ = size(z, 1)
      allocate(Jin(101))
      allocate(Zin(101))
      allocate(Pin(101))
      allocate(Tin(101))
      allocate(RHOin(101))
      allocate(dsz(NZ))
      allocate(dst(NZ))
      allocate(dsp(20000))
      allocate(dsw(NZ))
      allocate(pdense(20000))
      allocate(dsrho(20000))
      allocate(zdense(20000))
      allocate(rhodense(20000))


      open(unit=clm, file="clima_allout.tab",status="old")
      do il = 1, 1135
        read(clm,*)
      end do

      do il = 1,101
        read(clm,*)   Jin(102-il), Pin(102-il), Zin(102-il), Tin(102-il)
        Zin(102-il) = Zin(102-il) * 1000.
        Pin(102-il) = Pin(102-il) * 1e5
        RHOin(102-il) = 10000.*Pin(102-il)/(R_AIR)/Tin(102-il)
      end do
      zdense = [(j*2.,j=1,20000,1)]
!      gcmgrav = 9.8
      call spline (Zin, Pin, 101, zdense, pdense, dsp, 20000, 100, 99)
      call spline (Zin, RHOin, 101, zdense, rhodense, dsrho, 20000, 100, 99)
      call spline (Zin, Pin, 101, 0., pz0, dsp(1), 1, 100, 99)
      !! cubic spline interpolation
      do iz = 1, NZ
        pz = pz0
        do j = 1, 20000
          if (z(iz) .lt. zdense(j)) exit
!          pz=pz-2.8*rhodense(j) * gcmgrav
           pz = pdense(j)
        end do
        call spline (Zin, Tin, 101, z(iz), tz, dst(iz), 1, 100, 99)
!        call spline (Zin, Pin, 40, z(iz), pz, dsp(iz), 1, 39, 38)
!        call spline (Zin, Win, 101, z(iz), wz, dsw(iz), 1, 100, 99)
!        call spline (Zin, Qin, 101, z(iz), qz, dsw(iz), 1, 100, 99)
!        do il = 2,39
!          if(z(iz) .lt. Zin(il)) exit
!        end do
!        frac = (z(iz) - Zin(il-1)) / (Zin(il)- Zin(il-1))
!        pz   = Pin(il-1) - frac * (Pin(il-1)-Pin(il))
!        tz   = Tin(il-1) - frac * (Tin(il-1)-Tin(il))
!        wz   = Win(il-1) - frac * (Win(il-1)-Win(il))
        if (present(p)) p(iz) = pz
        if (present(t)) t(iz) = tz
      end do

      close(unit=clm)
  
      return
    end subroutine
  
end module

subroutine spline( x, y, n, sx, f0, f1, m, n1, n2 ) 
    use carma_precision_mod

    implicit none
    integer                :: i, j, k, m, n, n1, n2
    real(kind=f)           :: x(n), y(n), sx(m), f0(m), f1(m)
    real(kind=f)           :: s2(n), h(n1), dy(n1), s(n1), e(n2)
    real(kind=f)           :: z, h1, h2, h3, h4
  
    do i = 1, n1
        h(i)  = x(i+1) - x(i)
!        write(*,*) y(i+1), y(i)
        dy(i) = ( y(i+1) - y(i) ) / h(i)
    end do
    
    s2(1) = 0.d0; s2(n) = 0.d0
    do i = 2, n1
        s2(i) = 6.d0 * ( dy(i) - dy(i-1) )
    end do
    
    z = 0.5d0 / ( h(1) + h(2) )
    s(1) = -h(2) * z
    e(1) = s2(2) * z
    do i = 2, n2
        k    = i - 1
        j    = i + 1
        z    = 1.d0 / ( 2.d0*( h(i)+h(j) ) + h(i)*s(k) )
        s(i) = -h(j) * z
        e(i) = ( s2(j)-h(i)*e(k) ) * z
    end do
    
    s2(n1) = e(n2)
    do i = n2, 2, -1
        k     = i - 1
        s2(i) = s(k)*s2(i+1) + e(k)
    end do
    
    do i = 1, n1
        s(i) = ( s2(i+1) - s2(i) ) / h(i)
    end do
    
    i = 2
    k = 1
    do j = 1, m
        do 
            if ( sx(j) > x(i) ) then
                k = i
                i = i + 1
            else
                exit
            end if
        end do
        h1    = sx(j) - x(k)
        h2    = sx(j) - x(i)
        h3    = h1 * h2
        h4    = s2(k) + h1*s(k)
        z     = ( s2(i) + s2(k) + h4 ) / 6.d0
        f0(j)  = y(k) + h1*dy(k) + h3*Z
        f1(j) = dy(k) + z*( h1+h2 ) + h3 * s(k) / 6.d0
    end do
    
end subroutine spline
