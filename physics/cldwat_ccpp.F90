#undef DEBUG
module cldwat_ccpp
  use shr_kind_mod,  only: r8 => shr_kind_r8
  use wv_saturation, only: estblf, hlatv, tmin, hlatf, rgasv, pcf, &
                           cp, epsqs, ttrice
  
  implicit none

!-----------------------------------------------------------------------
! PUBLIC: Make default data and interfaces private
!-----------------------------------------------------------------------
  private
  
  save
  
  public inimc 
  public cldwat_fice
  
  integer, public::  ktop      ! Level above 10 hPa

  real(r8),public ::  icritc               ! threshold for autoconversion of cold ice
  real(r8),public ::  icritw               ! threshold for autoconversion of warm ice
!!$   real(r8),public,parameter::  conke  = 1.e-6    ! tunable constant for evaporation of precip
!!$   real(r8),public,parameter::  conke  =  2.e-6    ! tunable constant for evaporation of precip
  real(r8),public ::  conke                ! tunable constant for evaporation of precip
  real(r8),public ::  r3lcrit              ! critical radius where liq conversion begins

!-----------------------------------------------------------------------
! PRIVATE: Everything else is private to this module
!-----------------------------------------------------------------------
  real(r8), private:: tmax_fice! max temperature for cloud ice formation
  real(r8), private:: tmin_fice! min temperature for cloud ice formation
  real(r8), private:: tmax_fsnow ! max temperature for transition to convective snow
  real(r8), private:: tmin_fsnow ! min temperature for transition to convective snow
  real(r8), private:: rhonot   ! air density at surface
  real(r8), private:: t0       ! Freezing temperature
  real(r8), private:: cldmin   ! assumed minimum cloud amount
  real(r8), private:: small    ! small number compared to unity
  real(r8), private:: c        ! constant for graupel like snow cm**(1-d)/s
  real(r8), private:: d        ! constant for graupel like snow
  real(r8), private:: esi      ! collection efficient for ice by snow
  real(r8), private:: esw      ! collection efficient for water by snow
  real(r8), private:: nos      ! particles snow / cm**4
  real(r8), private:: pi       ! Mathematical constant
  real(r8), private:: gravit   ! Gravitational acceleration at surface
  real(r8), private:: rh2o
  real(r8), private:: prhonos
  real(r8), private:: thrpd    ! numerical three added to d
  real(r8), private:: gam3pd   ! gamma function on (3+d)
  real(r8), private:: gam4pd   ! gamma function on (4+d)
  real(r8), private:: rhoi     ! ice density
  real(r8), private:: rhos     ! snow density
  real(r8), private:: rhow     ! water density
  real(r8), private:: mcon01   ! constants used in cloud microphysics
  real(r8), private:: mcon02   ! constants used in cloud microphysics
  real(r8), private:: mcon03   ! constants used in cloud microphysics
  real(r8), private:: mcon04   ! constants used in cloud microphysics
  real(r8), private:: mcon05   ! constants used in cloud microphysics
  real(r8), private:: mcon06   ! constants used in cloud microphysics
  real(r8), private:: mcon07   ! constants used in cloud microphysics
  real(r8), private:: mcon08   ! constants used in cloud microphysics

  integer, private ::  k1mb    ! index of the eta level near 1 mb

! Parameters used in findmcnew
  real(r8) :: capnsi               ! sea ice cloud particles / cm3
  real(r8) :: capnc                ! cold and oceanic cloud particles / cm3
  real(r8) :: capnw                ! warm continental cloud particles / cm3
  real(r8) :: kconst               ! const for terminal velocity (stokes regime)
  real(r8) :: effc                 ! collection efficiency
  real(r8) :: alpha                ! ratio of 3rd moment radius to 2nd
  real(r8) :: capc                 ! constant for autoconversion
  real(r8) :: convfw               ! constant used for fall velocity calculation
  real(r8) :: cracw                ! constant used for rain accreting water
  real(r8) :: critpr               ! critical precip rate collection efficiency changes
  real(r8) :: ciautb               ! coefficient of autoconversion of ice (1/s)

#ifdef DEBUG
  integer, private,parameter ::  nlook = 1  ! Number of points to examine
  integer, private ::  ilook(nlook)         ! Longitude index to examine
  integer, private ::  latlook(nlook)       ! Latitude index to examine
  integer, private ::  lchnklook(nlook)     ! Chunk index to examine
  integer, private ::  icollook(nlook)      ! Column index to examine
#endif
 ! Private data
  real(r8), parameter :: unset_r8 = huge(1.0_r8)
  
  contains

!================================================================================================
  subroutine cldwat_fice(ncol, pcols, pver, t, fice, fsnow)
!
! Compute the fraction of the total cloud water which is in ice phase.
! The fraction depends on temperature only. 
! This is the form that was used for radiation, the code came from cldefr originally
! 
! Author: B. A. Boville Sept 10, 2002
!  modified: PJR 3/13/03 (added fsnow to ascribe snow production for convection )
!-----------------------------------------------------------------------
    implicit none

! Arguments
    integer,  intent(in)  :: ncol                 ! number of active columns
    integer,  intent(in)  :: pcols, pver
    real(r8), intent(in)  :: t(pcols,pver)        ! temperature

    real(r8), intent(out) :: fice(pcols,pver)     ! Fractional ice content within cloud
    real(r8), intent(out) :: fsnow(pcols,pver)    ! Fractional snow content for convection

! Local variables
    integer :: i,k                                   ! loop indexes

!-----------------------------------------------------------------------

! Define fractional amount of cloud that is ice
    do k=1,pver
       do i=1,ncol

! If warmer than tmax then water phase
          if (t(i,k) > tmax_fice) then
             fice(i,k) = 0.0_r8

! If colder than tmin then ice phase
          else if (t(i,k) < tmin_fice) then
             fice(i,k) = 1.0_r8

! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else 
             fice(i,k) =(tmax_fice - t(i,k)) / (tmax_fice - tmin_fice)
          end if

! snow fraction partitioning

! If warmer than tmax then water phase
          if (t(i,k) > tmax_fsnow) then
             fsnow(i,k) = 0.0_r8

! If colder than tmin then ice phase
          else if (t(i,k) < tmin_fsnow) then
             fsnow(i,k) = 1.0_r8

! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else 
             fsnow(i,k) =(tmax_fsnow - t(i,k)) / (tmax_fsnow - tmin_fsnow)
          end if

       end do
    end do

    return
  end subroutine cldwat_fice
  
  subroutine inimc( tmeltx, rhonotx, gravitx, rh2ox, hypm, microp_scheme, iulog, pver, masterproc)
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: 
  ! initialize constants for the prognostic condensate
  ! 
  ! Author: P. Rasch, April 1997
  ! 
  !-----------------------------------------------------------------------

     integer k
     real(r8), intent(in) :: tmeltx
     real(r8), intent(in) :: rhonotx
     real(r8), intent(in) :: gravitx
     real(r8), intent(in) :: rh2ox
     real(r8), intent(in), dimension(:) :: hypm
     integer,  intent(in) :: iulog, pver
     logical,  intent(in) :: masterproc

#ifdef UNICOSMP
     real(r8) signgam              ! variable required by cray gamma function
     external gamma
#endif
     character(len=16), intent(in)    :: microp_scheme 

  ! Set following for all physics packages

     tmax_fice = tmeltx    - 10._r8
  !! tmax_fice = tmeltx
  !! tmin_fice = tmax_fice - 20.
     tmin_fice = tmax_fice - 30._r8
     tmax_fsnow = tmeltx
     tmin_fsnow = tmeltx   - 5._r8

  ! Set remaining for RK microphysics

     if( microp_scheme .eq. 'RK' ) then
        rhonot = rhonotx             ! air density at surface (gm/cm3)
        gravit = gravitx
        rh2o   = rh2ox
        rhos = .1_r8                 ! assumed snow density (gm/cm3)
        rhow = 1._r8                 ! water density
        rhoi = 1._r8                 ! ice density
        esi = 1.0_r8                 ! collection efficient for ice by snow
        esw = 0.1_r8                 ! collection efficient for water by snow
        t0 = tmeltx                  ! approximate freezing temp
        cldmin = 0.02_r8             ! assumed minimum cloud amount
        small = 1.e-22_r8            ! a small number compared to unity
        c = 152.93_r8                ! constant for graupel like snow cm**(1-d)/s
        d = 0.25_r8                  ! constant for graupel like snow
        nos = 3.e-2_r8               ! particles snow / cm**4
        pi = 4._r8*atan(1.0_r8)
        prhonos = pi*rhos*nos
        thrpd = 3._r8 + d
        if (d==0.25_r8) then
           gam3pd = 2.549256966718531_r8 ! only right for d = 0.25
           gam4pd = 8.285085141835282_r8
        else
#ifdef UNICOSMP
           call gamma(3._r8+d, signgam, gam3pd)
           gam3pd = sign(exp(gam3pd),signgam)
           call gamma(4._r8+d, signgam, gam4pd)
           gam4pd = sign(exp(gam4pd),signgam)
           write(iulog,*) ' d, gamma(3+d), gamma(4+d) =', gam3pd, gam4pd
#else
           write(iulog,*) ' can only use d ne 0.25 on a cray '
           stop
#endif
        endif
        mcon01 = pi*nos*c*gam3pd/4._r8
        mcon02 = 1._r8/(c*gam4pd*sqrt(rhonot)/(6*prhonos**(d/4._r8)))
        mcon03 = -(0.5_r8+d/4._r8)
        mcon04 = 4._r8/(4._r8+d)
        mcon05 = (3+d)/(4+d)
        mcon06 = (3+d)/4._r8
        mcon07 = mcon01*sqrt(rhonot)*mcon02**mcon05/prhonos**mcon06
        mcon08 = -0.5_r8/(4._r8+d)

  !  find the level about 1mb, we wont do the microphysics above this level
        k1mb = 1
        do k=1,pver-1
           if (hypm(k) < 1.e2_r8 .and. hypm(k+1) >= 1.e2_r8) then
              if (1.e2_r8-hypm(k) < hypm(k+1)-1.e2_r8) then
                 k1mb = k
              else
                 k1mb = k + 1
              end if
              goto 20
           end if
        end do
        if (masterproc) then
           write(iulog,*)'inimc: model levels bracketing 1 mb not found'
        end if
  !     call endrun
        k1mb = 1
  20    if( masterproc ) write(iulog,*)'inimc: model level nearest 1 mb is',k1mb,'which is',hypm(k1mb),'pascals'

        if( masterproc ) write(iulog,*) 'cloud water initialization by inimc complete '

  ! Initialize parameters used by findmcnew
        capnw = 400._r8              ! warm continental cloud particles / cm3
        capnc = 150._r8              ! cold and oceanic cloud particles / cm3
  !     capnsi = 40._r8              ! sea ice cloud particles density  / cm3
        capnsi = 75._r8              ! sea ice cloud particles density  / cm3

        kconst = 1.18e6_r8           ! const for terminal velocity

  !     effc = 1._r8                 ! autoconv collection efficiency following boucher 96
  !     effc = .55*0.05_r8           ! autoconv collection efficiency following baker 93
        effc = 0.55_r8               ! autoconv collection efficiency following tripoli and cotton
  !     effc = 0._r8    ! turn off warm-cloud autoconv
        alpha = 1.1_r8**4
        capc = pi**(-.333_r8)*kconst*effc *(0.75_r8)**(1.333_r8)*alpha  ! constant for autoconversion

  ! critical precip rate at which we assume the collector drops can change the
  ! drop size enough to enhance the auto-conversion process (mm/day)
        critpr = 0.5_r8
        convfw = 1.94_r8*2.13_r8*sqrt(rhow*1000._r8*9.81_r8*2.7e-4_r8)

  ! liquid microphysics
  !     cracw = 6_r8                 ! beheng
        cracw = .884_r8*sqrt(9.81_r8/(rhow*1000._r8*2.7e-4_r8)) ! tripoli and cotton

  ! ice microphysics
        ciautb = 5.e-4_r8

        if ( masterproc ) then
           write(iulog,*)'tuning parameters cldwat: icritw',icritw,'icritc',icritc,'conke',conke,'r3lcrit',r3lcrit
           write(iulog,*)'tuning parameters cldwat: capnw',capnw,'capnc',capnc,'capnsi',capnsi,'kconst',kconst
           write(iulog,*)'tuning parameters cldwat: effc',effc,'alpha',alpha,'capc',capc
           write(iulog,*)'tuning parameters cldwat: critpr',critpr,'convfw',convfw,'cracw',cracw,'ciautb',ciautb
        endif
     endif

     return
  end subroutine inimc
  
  !GJF: model-agnostic subroutines pcond, findmcnew, and findsp from the cldwat module need to be added here to support MP schemes
  
end module cldwat_ccpp
