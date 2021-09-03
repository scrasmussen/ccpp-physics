module zm_conv_common
  
  use shr_kind_mod,    only: r8 => shr_kind_r8
  
  implicit none
  
  public
  
  save
  
  real(r8), parameter :: unset_r8 = huge(1.0_r8)
  real(r8) :: zmconv_c0_lnd = unset_r8    
  real(r8) :: zmconv_c0_ocn = unset_r8    
  real(r8) :: zmconv_ke     = unset_r8    
!------------- zhh added 2013-04-24 -----------------
  real(r8) :: zmconv_dmpdz  = unset_r8    
!------------- zhh added 2013-04-24 -----------------
  
  real(r8) rl         ! wg latent heat of vaporization.
  real(r8) cpres      ! specific heat at constant pressure in j/kg-degk.
  real(r8), parameter :: capelmt = 70._r8  ! threshold value for cape for deep convection.
  real(r8) :: ke           ! Tunable evaporation efficiency set from namelist input zmconv_ke
  real(r8) :: c0_lnd       ! set from namelist input zmconv_c0_lnd
  real(r8) :: c0_ocn       ! set from namelist input zmconv_c0_ocn
  real(r8) tau   ! convective time scale
  real(r8),parameter :: a = 21.656_r8
  real(r8),parameter :: b = 5418._r8
  real(r8),parameter :: c1 = 6.112_r8
  real(r8),parameter :: c2 = 17.67_r8
  real(r8),parameter :: c3 = 243.5_r8
  real(r8) :: tfreez
  real(r8) :: eps1
  real(r8) :: latice
  real(r8) :: cpwv
  real(r8) :: cpliq
  real(r8) :: rh2o
     
  logical :: no_deep_pbl ! default = .false.
                         ! no_deep_pbl = .true. eliminates deep convection entirely within PBL 

!moved from moistconvection.F90
  real(r8) :: rgrav       ! reciprocal of grav
  real(r8) :: rgas        ! gas constant for dry air
  real(r8) :: grav        ! = gravit
  
  integer  limcnv       ! top interface level limit for convection

  real(r8),parameter ::  tiedke_add = 0.5_r8   
  
end module zm_conv_common
