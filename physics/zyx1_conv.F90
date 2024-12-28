
module zyx1_conv

!---------------------------------------------------------------------------------
! Purpose:
!
! ZYX1 convection scheme, code structure modified from zm_conv.F90 in CAM
!
! Zhang Minghua, 2018-07-24 
!
!---------------------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use spmd_utils,      only: masterproc
  use ppgrid,          only: pcols, pver, pverp
  use cldwat,          only: cldwat_fice
  use physconst,       only: cpair, epsilo, gravit, latice, latvap, tmelt, rair, &
                             cpwv, cpliq, rh2o
  use abortutils,      only: endrun
  use cam_logfile,     only: iulog
!zmh
  use mzfunctions_mod, only: zmh_ramp, subcol_kissvec

  implicit none

  save
  private                         ! Make default type private to the module
!
! PUBLIC: interfaces
!
  public zyx1_conv_readnl            ! read zmconv_nl namelist
  public zyx1_conv_init            ! ZYX1 schemea
  public zyx1_conv_tend                ! ZYX1
  public zyx1_conv_evap             ! evaporation of precip from ZM schemea
  public convtran                 ! convective transport
  public momtran                  ! convective momentum transport

!
! Private data
!
   real(r8), parameter :: unset_r8 = huge(1.0_r8)

! == zyxconv_nl namelist default values =================================

   integer  :: niter = 4      ! iteration number of CAPE and vertical velocity
   character(len=16) :: deep_trigger_scheme = 'zmh_triggers' !'cam'
   character(len=16) :: deep_closure_scheme = 'prognostic' !'cam'
   logical  :: zmh_deep_stochastic          = .true. 
   logical  :: zmh_deep_nolimit = .false. 
   logical  :: do_dcape         = .true. 
   real(r8) :: zmh_cape_lmt     = 250._r8
   real(r8) :: zmh_mb           = 1.0_r8
   real(r8) :: zmh_alfa         = 0.3_r8  
   real(r8) :: zmh_tau          = 1800._r8 ! convective time scale
   real(r8) :: zmh_dmpdz        = -1.5E-4_r8 ! initial guess
   real(r8) :: zmh_ke           = 20.0E-6_r8

! for precip microphysics in convective plume   
   real(r8) :: zmh_c0_lnd       = 0.020_r8 !0.0059 
   real(r8) :: zmh_c0_ocn       = 0.030_r8 !0.045

!the parameters below are not used for now 
   real(r8) :: zmh_rain_z0      = 100.0_r8  !meter
   real(r8) :: zmh_rain_zp      = 500._r8 !100.0 !meter over ocean
   real(r8) :: zmh_rain_cr      = 0.03_r8
   real(r8) :: zmh_rain_lnd_z0  = 100.0_r8
   real(r8) :: zmh_rain_lnd_zp  = 500._r8
   real(r8) :: zmh_rain_lnd_cr  = 0.03_r8 !
   real(r8) :: zmh_deep_tau1    = 1200._r8  ! for prognostic cape , time to mb
   real(r8) :: zmh_deep_tau2    = 1200._r8  ! for memory of 1000 J/kg CAPE

! for deep cld amt
   real(r8) :: zmh_dfrc0        = 2.0_r8 !  core cloud fraction coef
   real(r8) :: zmh_dfrcd        = 0.0_r8 !  volume fraction of detained air
   real(r8) :: zmh_dfrcz        = 0.0_r8 !

! for srf inhomogeity
   logical  :: zmh_deep_tpertt  = .false.
   logical  :: zmh_deep_qpertt  = .false.
   logical  :: zmh_deep_qpertf  = .false.
   logical  :: zmh_deep_tpertf  = .false.
   logical  :: zmh_deep_tpertb  = .false.
   logical  :: zmh_deep_qperth  = .false.
   real(r8) :: zmh_deep_xyt     = 0.e3_r8    
   real(r8) :: zmh_deep_xyq     = 0.e3_r8     
   real(r8) :: zmh_deep_tpertmax= 3.0_r8
   real(r8) :: zmh_deep_qpertmax= 2.0e-3_r8

   real(r8) :: zmh_ent_alpha0   = 0.01_r8
   real(r8) :: zmh_ent_turb0    = 5.0e-4_r8
   real(r8) :: zmh_det_turb0    = 0.5e-4_r8
   real(r8) :: zmh_ent_ze       = 1.0_r8
   real(r8) :: zmh_dzt          = 0.2_r8 !0.5 !0.1
   real(r8) :: zmh_dzb          = 0.2_r8 !0.5 !0.9

!  ============================================
   real(r8) rl         ! wg latent heat of vaporization.
   real(r8) cpres      ! specific heat at constant pressure in j/kg-degk.
   real(r8) :: ke     ! Tunable evaporation efficiency set from namelist 
   real(r8),parameter :: a = 21.656_r8
   real(r8),parameter :: b = 5418._r8
   real(r8),parameter :: c1 = 6.112_r8
   real(r8),parameter :: c2 = 17.67_r8
   real(r8),parameter :: c3 = 243.5_r8
   real(r8) :: tfreez
   real(r8) :: eps1

!moved from moistconvection.F90
   real(r8) :: rgrav       ! reciprocal of grav
   real(r8) :: rgas        ! gas constant for dry air
   real(r8) :: grav        ! = gravit
   real(r8) :: cp          ! = cpres = cpair
   
   integer  limcnv       ! top interface level limit for convection

   real(r8),parameter ::  tiedke_add = 0.5_r8   

contains

subroutine zyx1_conv_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'zyxconv_readnl'

   namelist /zyxconv_nl/   deep_trigger_scheme  , deep_closure_scheme    , &
       zmh_deep_nolimit  , do_dcape             , zmh_cape_lmt           , &
       zmh_mb            , zmh_alfa             , zmh_tau                , & 
       zmh_dmpdz         , zmh_ke               , zmh_deep_stochastic    , &
       zmh_c0_lnd        , zmh_c0_ocn           , zmh_rain_z0       , &
       zmh_rain_zp       , zmh_rain_cr          , zmh_rain_lnd_z0        , &
       zmh_rain_lnd_zp   , zmh_rain_lnd_cr      , zmh_deep_tau1          , &
       zmh_deep_tau2     , zmh_dfrc0            , zmh_dfrcd    ,zmh_dfrcz, & 
       zmh_deep_qpertf   , zmh_deep_qperth      , zmh_deep_tpertf        , &
       zmh_deep_tpertb   , zmh_deep_xyt         , zmh_deep_xyq           , &
       zmh_deep_tpertmax , zmh_deep_qpertmax    , zmh_ent_alpha0         , &
       zmh_ent_ze        , zmh_dzt              , zmh_dzb                , &
       zmh_deep_qpertt   , zmh_deep_tpertt      , zmh_ent_turb0, zmh_det_turb0
   !-----------------------------------------------------------------------------

   if(masterproc)then 
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'zyxconv_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, zyxconv_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist2')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

   ke = zmh_ke
   print*,''
   print*,'zyxconv_nl namelist values in atm_in are:'
   print*,'/zyxconv_nl/   deep_trigger_scheme  , deep_closure_scheme    ', &
       'zmh_deep_nolimit  , do_dcape             , zmh_cape_lmt          ' , &
       'zmh_mb            , zmh_alfa             , zmh_tau               ' , & 
       'zmh_dmpdz         , zmh_c0_lnd        , zmh_c0_ocn           ,zmh_rain_z0  ', &
       'zmh_rain_zp       , zmh_rain_cr          , zmh_rain_lnd_z0       ' , &
       'zmh_rain_lnd_zp   , zmh_rain_lnd_cr      , zmh_deep_tau1         ' , &
       'zmh_deep_tau2     , zmh_dfrc0            , zmh_dfrcd    ,zmh_dfrcz', & 
       'zmh_deep_qpertf   , zmh_deep_qperth      , zmh_deep_tpertf        ', &
       'zmh_deep_tpertb   , zmh_deep_xyt         , zmh_deep_xyq           ', &
       'zmh_deep_tpertmax , zmh_deep_qpertmax    , ke , zmh_deep_stochastic',&
       'zmh_ent_alpha0    , zmh_ent_ze        , zmh_dzt              , zmh_dzb',&
       'zmh_deep_qpertt   , zmh_deep_tpertt   , zmh_ent_turb0, , zmh_det_turb0'

   print*, deep_trigger_scheme  , deep_closure_scheme    , &
       zmh_deep_nolimit  , do_dcape             , zmh_cape_lmt           , &
       zmh_mb            , zmh_alfa             , zmh_tau                , & 
       zmh_dmpdz         , zmh_c0_lnd        , zmh_c0_ocn           ,zmh_rain_z0 , &
       zmh_rain_zp       , zmh_rain_cr          , zmh_rain_lnd_z0        , &
       zmh_rain_lnd_zp   , zmh_rain_lnd_cr      , zmh_deep_tau1          , &
       zmh_deep_tau2     , zmh_dfrc0            , zmh_dfrcd    ,zmh_dfrcz, & 
       zmh_deep_qpertf   , zmh_deep_qperth      , zmh_deep_tpertf        , &
       zmh_deep_tpertb   , zmh_deep_xyt         , zmh_deep_xyq           , &
       zmh_deep_tpertmax , zmh_deep_qpertmax    ,ke, zmh_deep_stochastic , &
       zmh_ent_alpha0    , zmh_ent_ze        , zmh_dzt              , zmh_dzb,&
       zmh_deep_qpertt   , zmh_deep_tpertt   , zmh_ent_turb0,  zmh_det_turb0

   end if

#ifdef SPMD
   ! Broadcast namelist variables

   call mpibcast(zmh_deep_nolimit,    1, mpilog,  0, mpicom)
   call mpibcast(deep_trigger_scheme, 1, mpichar,  0, mpicom)
   call mpibcast(deep_closure_scheme, 1, mpichar,  0, mpicom)
   call mpibcast(do_dcape,            1, mpilog,  0, mpicom)
   call mpibcast(zmh_cape_lmt,        1, mpir8,  0, mpicom)
   call mpibcast(zmh_mb,              1, mpir8,  0, mpicom)
   call mpibcast(zmh_alfa,            1, mpir8,  0, mpicom)
   call mpibcast(zmh_tau,             1, mpir8,  0, mpicom)   
   call mpibcast(zmh_dmpdz,           1, mpir8,  0, mpicom)  
   call mpibcast(ke,                  1, mpir8,  0, mpicom)  
   call mpibcast(zmh_deep_stochastic, 1, mpilog,  0, mpicom)  
   call mpibcast(zmh_c0_lnd ,         1, mpir8,  0, mpicom)
   call mpibcast(zmh_c0_ocn ,         1, mpir8,  0, mpicom)
   call mpibcast(zmh_rain_z0,         1, mpir8,  0, mpicom)
   call mpibcast(zmh_rain_zp,         1, mpir8,  0, mpicom)
   call mpibcast(zmh_rain_cr,         1, mpir8,  0, mpicom)
   call mpibcast(zmh_rain_lnd_z0,     1, mpir8,  0, mpicom)
   call mpibcast(zmh_rain_lnd_zp,     1, mpir8,  0, mpicom)
   call mpibcast(zmh_rain_lnd_cr,     1, mpir8,  0, mpicom)
   call mpibcast(zmh_deep_tau1,       1, mpir8,  0, mpicom)
   call mpibcast(zmh_deep_tau2,       1, mpir8,  0, mpicom)

   call mpibcast(zmh_dfrc0,           1, mpir8,  0, mpicom)
   call mpibcast(zmh_dfrcd,           1, mpir8,  0, mpicom)
   call mpibcast(zmh_dfrcz,           1, mpir8,  0, mpicom)

   call mpibcast(zmh_deep_qpertt,     1, mpilog,  0, mpicom)
   call mpibcast(zmh_deep_tpertt,     1, mpilog,  0, mpicom)

   call mpibcast(zmh_deep_qpertf,     1, mpilog,  0, mpicom)
   call mpibcast(zmh_deep_qperth,     1, mpilog,  0, mpicom)
   call mpibcast(zmh_deep_tpertf,     1, mpilog,  0, mpicom)
   call mpibcast(zmh_deep_tpertb,     1, mpilog,  0, mpicom)
   call mpibcast(zmh_deep_xyt,        1, mpir8,  0, mpicom)
   call mpibcast(zmh_deep_xyq,        1, mpir8,  0, mpicom)
   call mpibcast(zmh_deep_tpertmax,   1, mpir8,  0, mpicom)
   call mpibcast(zmh_deep_qpertmax,   1, mpir8,  0, mpicom)

   call mpibcast(zmh_ent_alpha0 ,   1, mpir8,  0, mpicom)
   call mpibcast(zmh_ent_turb0  ,   1, mpir8,  0, mpicom)
   call mpibcast(zmh_det_turb0  ,   1, mpir8,  0, mpicom)
   call mpibcast(zmh_ent_ze     ,   1, mpir8,  0, mpicom)
   call mpibcast(zmh_dzt        ,   1, mpir8,  0, mpicom)
   call mpibcast(zmh_dzb        ,   1, mpir8,  0, mpicom)
#endif

end subroutine zyx1_conv_readnl


subroutine zyx1_conv_init(limcnv_in)

   use dycore,       only: dycore_is, get_resolution

   integer, intent(in)           :: limcnv_in       ! top interface level limit for convection

   ! local variables
   character(len=32)   :: hgrid           ! horizontal grid specifier

   ! Initialization of constants
   limcnv = limcnv_in
   tfreez = tmelt
   eps1   = epsilo
   rl     = latvap
   cpres  = cpair
   rgrav  = 1.0_r8/gravit
   rgas   = rair
   grav   = gravit
   cp     = cpres

end subroutine zyx1_conv_init



subroutine zyx1_conv_tend(lchnk      ,ncol    ,nstep            , &
                    t       ,qh      ,prec    ,jctop   ,jcbot   , &
                    pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                    heat    ,pap     ,paph    ,dpp     ,delt    , &
                    lat     ,lon     ,landfrac,lhflx   ,shflx   , & 
                    psrf    ,bfls_t  ,bfls_q                    , &
                    omega   ,tpert2  ,qpert2  ,gradt   ,gradq   , &  
                    vort3   ,tke     ,massflxbase_p             , &
                    mcon    ,cme     ,cape                      , &
                    dlf     ,pflx    ,zdu     ,rprd             , &
                    mu      ,md      ,du      ,eu      ,ed      , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                    lengath ,ql      ,rliq    ,slflx   ,qtflx   , &
                    deepfrc )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Main driver for zhang-mcfarlane convection scheme 
! 
! Method: 
! performs deep convective adjustment based on mass-flux closure
! algorithm.
! 
! Author:guang jun zhang, m.lazare, n.mcfarlane. CAM Contact: P. Rasch
!
! This is contributed code not fully standardized by the CAM core group.
! All variables have been typed, where most are identified in comments
! The current procedure will be reimplemented in a subsequent version
! of the CAM where it will include a more straightforward formulation
! and will make use of the standard CAM nomenclature
! 
!-----------------------------------------------------------------------
   use constituents, only: pcnst
   use phys_control, only: cam_physpkg_is
   use wv_saturation,only: aqsat
!zmh
   use comsrf,       only: sgh30
   use time_manager, only: get_curr_calday  


!
! ************************ index of variables **********************
!
!  wg * alpha    array of vertical differencing used (=1. for upstream).
!  w  * cape     convective available potential energy.
!  wg * capeg    gathered convective available potential energy.
!  c  * capelmt  threshold value for cape for deep convection.
!  ic  * cpres    specific heat at constant pressure in j/kg-degk.
!  i  * dpp      
!  ic  * delt     length of model time-step in seconds.
!  wg * dp       layer thickness in mbs (between upper/lower interface).
!  wg * dqdt     mixing ratio tendency at gathered points.
!  wg * dsdt     dry static energy ("temp") tendency at gathered points.
!  wg * dudt     u-wind tendency at gathered points.
!  wg * dvdt     v-wind tendency at gathered points.
!  wg * dsubcld  layer thickness in mbs between lcl and maxi.
!  ic  * grav     acceleration due to gravity in m/sec2.
!  wg * du       detrainment in updraft. specified in mid-layer
!  wg * ed       entrainment in downdraft.
!  wg * eu       entrainment in updraft.
!  wg * hmn      moist static energy.
!  wg * hsat     saturated moist static energy.
!  w  * ideep    holds position of gathered points vs longitude index.
!  ic  * pver     number of model levels.
!  wg * j0       detrainment initiation level index.
!  wg * jd       downdraft   initiation level index.
!  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
!  wg * jt       top  level index of deep cumulus convection.
!  w  * lcl      base level index of deep cumulus convection.
!  wg * lclg     gathered values of lcl.
!  w  * lel      index of highest theoretical convective plume.
!  wg * lelg     gathered values of lel.
!  w  * lon      index of onset level for deep convection.
!  w  * maxi     index of level with largest moist static energy.
!  wg * maxg     gathered values of maxi.
!  wg * mb       cloud base mass flux.
!  wg * mc       net upward (scaled by mb) cloud mass flux.
!  wg * md       downward cloud mass flux (positive up).
!  wg * mu       upward   cloud mass flux (positive up). specified
!                at interface
!  ic  * msg      number of missing moisture levels at the top of model.
!  w  * p        grid slice of ambient mid-layer pressure in mbs.
!  i  * pblt     row of pbl top indices.
!  w  * pcpdh    scaled surface pressure.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  wg * pg       grid slice of gathered values of p.
!  w  * q        grid slice of mixing ratio.
!  wg * qd       grid slice of mixing ratio in downdraft.
!  wg * qg       grid slice of gathered values of q.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  wg * qhat     grid slice of upper interface mixing ratio.
!  wg * ql       grid slice of cloud liquid water.
!  wg * qs       grid slice of saturation mixing ratio.
!  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
!  wg * qstpg    grid slice of gathered values of qstp.
!  wg * qu       grid slice of mixing ratio in updraft.
!  ic  * rgas     dry air gas constant.
!  wg * rl       latent heat of vaporization.
!  w  * s        grid slice of scaled dry static energy (t+gz/cp).
!  wg * sd       grid slice of dry static energy in downdraft.
!  wg * sg       grid slice of gathered values of s.
!  wg * shat     grid slice of upper interface dry static energy.
!  wg * su       grid slice of dry static energy in updraft.
!  i/o * t       
!  o  * jctop    row of top-of-deep-convection indices passed out.
!  O  * jcbot    row of base of cloud indices passed out.
!  wg * tg       grid slice of gathered values of t.
!  w  * tl       row of parcel temperature at lcl.
!  wg * tlg      grid slice of gathered values of tl.
!  w  * tp       grid slice of parcel temperatures.
!  wg * tpg      grid slice of gathered values of tp.
!  i/o * u        grid slice of u-wind (real).
!  wg * ug       grid slice of gathered values of u.
!  i/o * utg      grid slice of u-wind tendency (real).
!  i/o * v        grid slice of v-wind (real).
!  w  * va       work array re-used by called subroutines.
!  wg * vg       grid slice of gathered values of v.
!  i/o * vtg      grid slice of v-wind tendency (real).
!  i  * w        grid slice of diagnosed large-scale vertical velocity.
!  w  * z        grid slice of ambient mid-layer height in metres.
!  w  * zf       grid slice of ambient interface height in metres.
!  wg * zfg      grid slice of gathered values of zf.
!  wg * zg       grid slice of gathered values of z.
!
!-----------------------------------------------------------------------
!
! multi-level i/o fields:
!  i      => input arrays.
!  i/o    => input/output arrays.
!  w      => work arrays.
!  wg     => work arrays operating only on gathered points.
!  ic     => input data constants.
!  c      => data constants pertaining to subroutine itself.
!
! input arguments
!
   integer, intent(in) :: lchnk                   ! chunk identifier
   integer, intent(in) :: nstep                   ! step identifier
   integer, intent(in) :: ncol                    ! number of atmospheric columns

   real(r8), intent(in) :: t(pcols,pver)          ! grid slice of temperature at mid-layer.
   real(r8), intent(in) :: qh(pcols,pver,pcnst)   ! grid slice of specific humidity.
   real(r8), intent(in) :: pap(pcols,pver)     
   real(r8), intent(in) :: paph(pcols,pver+1)
   real(r8), intent(in) :: dpp(pcols,pver)        ! local sigma half-level thickness (i.e. dshj).
   real(r8), intent(in) :: zm(pcols,pver)
   real(r8), intent(in) :: geos(pcols)
   real(r8), intent(in) :: zi(pcols,pver+1)
   real(r8), intent(in) :: pblh(pcols)
   real(r8), intent(in) :: tpert2(pcols)
   real(r8), intent(in) :: landfrac(pcols) ! RBN Landfrac
! zmh
   real(r8), intent(in) :: lat(pcols)
   real(r8), intent(in) :: lon(pcols)
   real(r8), intent(in) :: lhflx(pcols)
   real(r8), intent(in) :: shflx(pcols)
   real(r8), intent(in) :: tke(pcols)
   real(r8), intent(in) :: psrf(pcols)
   real(r8), intent(in) :: qpert2(pcols)
   real(r8), intent(in) :: gradt(pcols)
   real(r8), intent(in) :: gradq(pcols)

   real(r8), intent(in) :: bfls_t(pcols,pver)          ! grid slice of temperature at mid-layer.
   real(r8), intent(in) :: bfls_q(pcols,pver)          ! grid slice of temperature at mid-layer.
   real(r8), intent(in) :: omega(pcols,pver)          ! grid slice of temperature at mid-layer.
   real(r8), intent(in) :: vort3(pcols)          ! grid slice of temperature at mid-layer.
   real(r8), intent(in) ::  delt                     ! length of model time-step in seconds.
   real(r8), intent(inout) :: massflxbase_p(pcols,pver)       !

   real(r8) :: tpert(pcols)
   real(r8) :: qpert(pcols)
   real(r8) :: coszrs(pcols)
   real(r8) :: capelmt(pcols)

   real(r8) :: ezmh_deep_tau1, ezmh_deep_tau2, w11, w12,w13,w14
   real(r8) :: calday
   
!
! output arguments
!
   real(r8), intent(out) :: qtnd(pcols,pver)           ! specific humidity tendency (kg/kg/s)
   real(r8), intent(out) :: heat(pcols,pver)           ! heating rate (dry static energy tendency, W/kg)
   real(r8), intent(out) :: mcon(pcols,pverp)
   real(r8), intent(out) :: dlf(pcols,pver)    ! scattrd version of the detraining cld h2o tend
   real(r8), intent(out) :: pflx(pcols,pverp)  ! scattered precip flux at each level
   real(r8), intent(out) :: cme(pcols,pver)
   real(r8), intent(out) :: cape(pcols)        ! w  convective available potential energy.
   real(r8), intent(out) :: zdu(pcols,pver)
   real(r8), intent(out) :: rprd(pcols,pver)     ! rain production rate
! move these vars from local storage to output so that convective
! transports can be done in outside of conv_cam.
   real(r8), intent(out) :: mu(pcols,pver)
   real(r8), intent(out) :: eu(pcols,pver)
   real(r8), intent(out) :: du(pcols,pver)
   real(r8), intent(out) :: md(pcols,pver)
   real(r8), intent(out) :: ed(pcols,pver)
   real(r8), intent(out) :: dp(pcols,pver)       ! 
   real(r8), intent(out) :: dsubcld(pcols)       ! wg layer thickness in mbs between lcl and maxi.
   real(r8), intent(out) :: jctop(pcols)  ! o row of top-of-deep-convection indices passed out.
   real(r8), intent(out) :: jcbot(pcols)  ! o row of base of cloud indices passed out.
   real(r8), intent(out) :: prec(pcols)
   real(r8), intent(out) :: rliq(pcols) ! reserved liquid (not yet in cldliq) for energy integrals
!wxc zmh
   real(r8), intent(out) :: slflx(pcols,pverp)
   real(r8), intent(out) :: qtflx(pcols,pverp)

   real(r8), intent(out) :: deepfrc(pcols,pver)


   real(r8) zs(pcols)
   real(r8) dlg(pcols,pver)    ! gathrd version of the detraining cld h2o tend
   real(r8) pflxg(pcols,pverp) ! gather precip flux at each level
   real(r8) cug(pcols,pver)    ! gathered condensation rate
   real(r8) evpg(pcols,pver)   ! gathered evap rate of rain in downdraft
   real(r8) mumax(pcols)
   integer jt(pcols)                          ! wg top  level index of deep cumulus convection.
   integer maxg(pcols)                        ! wg gathered values of maxi.
   integer ideep(pcols)                       ! w holds position of gathered points vs longitude index.
   integer lengath
!     diagnostic field used by chem/wetdep codes
   real(r8) ql(pcols,pver)                    ! wg grid slice of cloud liquid water.
!
   real(r8) pblt(pcols)           ! i row of pbl top indices.
!-----------------------------------------------------------------------
!
! general work fields (local variables):
!
   real(r8) q(pcols,pver)              ! w  grid slice of mixing ratio.
   real(r8) p(pcols,pver)              ! w  grid slice of ambient mid-layer pressure in mbs.
   real(r8) z(pcols,pver)              ! w  grid slice of ambient mid-layer height in metres.
   real(r8) s(pcols,pver)              ! w  grid slice of scaled dry static energy (t+gz/cp).
   real(r8) tp(pcols,pver)             ! w  grid slice of parcel temperatures.
   real(r8) zf(pcols,pver+1)           ! w  grid slice of ambient interface height in metres.
!zmh
   real(r8) dz(pcols,pver)           ! w  grid slice of ambient interface height in metres.
   real(r8) pf(pcols,pver+1)           ! w  grid slice of ambient interface pressure in mbs.
   real(r8) qstp(pcols,pver)           ! w  grid slice of parcel temp. saturation mixing ratio.

   real(r8) tl(pcols)                  ! w  row of parcel temperature at lcl.

   integer lcl(pcols)                  ! w  base level index of deep cumulus convection.
   integer lel(pcols)                  ! w  index of highest theoretical convective plume.
   integer maxi(pcols)                 ! w  index of level with largest moist static energy.

   integer index(pcols)
   real(r8) precip
!
! gathered work fields:
!
   real(r8) qg(pcols,pver)             ! wg grid slice of gathered values of q.
   real(r8) tg(pcols,pver)             ! w  grid slice of temperature at interface.
   real(r8) pg(pcols,pver)             ! wg grid slice of gathered values of p.
   real(r8) zg(pcols,pver)             ! wg grid slice of gathered values of z.
   real(r8) sg(pcols,pver)             ! wg grid slice of gathered values of s.
   real(r8) tpg(pcols,pver)            ! wg grid slice of gathered values of tp.
   real(r8) zfg(pcols,pver+1)          ! wg grid slice of gathered values of zf.
   real(r8) qstpg(pcols,pver)          ! wg grid slice of gathered values of qstp.
   real(r8) ug(pcols,pver)             ! wg grid slice of gathered values of u.
   real(r8) vg(pcols,pver)             ! wg grid slice of gathered values of v.
   real(r8) cmeg(pcols,pver)

   real(r8) rprdg(pcols,pver)           ! wg gathered rain production rate
   real(r8) capeg(pcols)               ! wg gathered convective available potential energy.
   real(r8) tlg(pcols)                 ! wg grid slice of gathered values of tl.
   real(r8) landfracg(pcols)            ! wg grid slice of landfrac  

   integer lclg(pcols)       ! wg gathered values of lcl.
   integer lelg(pcols)
!
! work fields arising from gathered calculations.
!
   real(r8) dqdt(pcols,pver)           ! wg mixing ratio tendency at gathered points.
   real(r8) dsdt(pcols,pver)           ! wg dry static energy ("temp") tendency at gathered points.
!      real(r8) alpha(pcols,pver)      ! array of vertical differencing used (=1. for upstream).
   real(r8) sd(pcols,pver)             ! wg grid slice of dry static energy in downdraft.
   real(r8) qd(pcols,pver)             ! wg grid slice of mixing ratio in downdraft.
   real(r8) mc(pcols,pver)             ! wg net upward (scaled by mb) cloud mass flux.
   real(r8) qhat(pcols,pver)           ! wg grid slice of upper interface mixing ratio.
   real(r8) qu(pcols,pver)             ! wg grid slice of mixing ratio in updraft.
   real(r8) su(pcols,pver)             ! wg grid slice of dry static energy in updraft.
   real(r8) qs(pcols,pver)             ! wg grid slice of saturation mixing ratio.
   real(r8) shat(pcols,pver)           ! wg grid slice of upper interface dry static energy.
   real(r8) hmn(pcols,pver)            ! wg moist static energy.
   real(r8) hsat(pcols,pver)           ! wg saturated moist static energy.
   real(r8) qlg(pcols,pver)
   real(r8) dudt(pcols,pver)           ! wg u-wind tendency at gathered points.
   real(r8) dvdt(pcols,pver)           ! wg v-wind tendency at gathered points.
!      real(r8) ud(pcols,pver)
!      real(r8) vd(pcols,pver)

   real(r8) mb(pcols)                  ! wg cloud base mass flux.

   integer jlcl(pcols)
   integer j0(pcols)                 ! wg detrainment initiation level index.
   integer jd(pcols)                 ! wg downdraft initiation level index.

   integer i
   integer ii
   integer k
   integer msg                      !  ic number of missing moisture levels at the top of model.
   integer k1,k2,k3,iter
   real (r8) qdifr
   real (r8) sdifr

!zmh local -----------
   logical,  dimension(pcols,10)    :: deep_triggers
   real(r8), dimension(pcols, pver) :: buoy_mid , bfls_buoy_mid
   real(r8), dimension(pcols, pver) :: est      , qsat    , rh
   real(r8), dimension(pcols)       :: bfls_cape, cin     , bfls_cin
   real(r8), dimension(pcols)       :: precw, precws
   real(r8), dimension(pcols)       :: ramp,rampg
   integer,  dimension(pcols)       :: nt_stepping, k750, k950, k100, kpbl
   logical  :: flag(pcols),flagg(pcols)
   integer  :: iflag   !loop index


   real (r8) ent_org(pcols, pver)
   real (r8) det_org(pcols, pver)
   real (r8) ent_turb(pcols, pver)
   real (r8) det_turb(pcols, pver)
   real (r8) ent_tot(pcols, pver)
   real (r8) det_tot(pcols, pver)
   real (r8) w_init(pcols)
   real (r8) rad_init(pcols)
   real (r8) w_up(pcols, pver)
   real (r8) buoy_up(pcols, pver)

   real (r8) dum1 ,dum0, dum2, dum3, dumb
   real (r8) dzb(pcols), dzt(pcols),eps00(pcols)
   integer   jdg(pcols)
   real (r8) eps00g(pcols)
   real (r8) ent_totg(pcols, pver), det_totg(pcols, pver)
   real (r8) rho(pcols, pver)
   real (r8) dsubcldg(pcols)

   integer maxi1(pcols),maxi1g(pcols)
   integer maxi2(pcols),maxi2g(pcols)
   integer lel2(pcols),lel2g(pcols)
   integer lelb(pcols),lelbg(pcols)
   integer jd1(pcols),jd1g(pcols)
   integer jd2(pcols),jd2g(pcols)
   integer jtg(pcols)
   integer k950g(pcols)
   real (r8) latg(pcols),long(pcols)

   integer seed1(pcols),seed2(pcols),seed3(pcols),seed4(pcols)
   real (r8) rand_num(pcols), rand_numw(pcols)

!
!--------------------------Data statements------------------------------
!
! Set internal variable "msg" (convection limit) to "limcnv-1"
!
   msg = limcnv - 1

!zmh
  rand_num(:) = 0.0  
  rand_numw(:) = 0.0  
  if(zmh_deep_stochastic)then 
   do i = 1,ncol
    seed1(i) = abs(t(i,1) - t(i,pver))  * 1000000
    seed2(i) = abs(t(i,3) - t(i,pver))  * 1000000
    seed3(i) = abs(t(i,5) - t(i,pver))  * 1000000
    seed4(i) = abs(t(i,7) - t(i,pver))  * 1000000
   enddo
   call subcol_kissvec(seed1, seed2, seed3, seed4, rand_num)
   call subcol_kissvec(seed4, seed3, seed2, seed1, rand_numw)
  endif

   capelmt(:) = zmh_cape_lmt*(1._r8+rand_num(:)) 

 !adjust for tropical lands
 do i=1,ncol
      if(cos(lat(i)*2.) .ge. 0.) then  !45NS
         if(landfrac(i) .gt. 0.99)then  
              capelmt(i) = capelmt(i)/2.0
!322
!              capelmt(i) = min(capelmt(i), 100.)
         endif
     endif   
 enddo
! initialize necessary arrays.
! zero out variables not used in cam
!
   qtnd(:,:) = 0._r8
   heat(:,:) = 0._r8
   mcon(:,:) = 0._r8
   rliq(:ncol)   = 0._r8
!
! initialize convective tendencies
!
   prec(:ncol) = 0._r8
   do k = 1,pver
      do i = 1,ncol
         dqdt(i,k)  = 0._r8
         dsdt(i,k)  = 0._r8
         dudt(i,k)  = 0._r8
         dvdt(i,k)  = 0._r8
         pflx(i,k)  = 0._r8
         pflxg(i,k) = 0._r8
         cme(i,k)   = 0._r8
         rprd(i,k)  = 0._r8
         zdu(i,k)   = 0._r8
         ql(i,k)    = 0._r8
         qlg(i,k)   = 0._r8
         dlf(i,k)   = 0._r8
         dlg(i,k)   = 0._r8
         deepfrc(i,k)   = 0._r8  !zmh
      end do
   end do
   do i = 1,ncol
      pflx(i,pverp) = 0
      pflxg(i,pverp) = 0
   end do

   slflx(:pcols,:pverp) = 0._r8
   qtflx(:pcols,:pverp) = 0._r8
!
   do i = 1,ncol
      kpbl(i) = pver
      pblt(i) = pver
      dsubcld(i) = 0._r8

      jctop(i) = pver
      jcbot(i) = 1

   end do
!
! calculate local pressure (mbs) and height (m) for both interface
! and mid-layer locations.
!
   do i = 1,ncol
      zs(i) = geos(i)*rgrav
      pf(i,pver+1) = paph(i,pver+1)*0.01_r8
      zf(i,pver+1) = zi(i,pver+1) + zs(i)
   end do
   do k = 1,pver
      do i = 1,ncol
         p(i,k) = pap(i,k)*0.01_r8
         pf(i,k) = paph(i,k)*0.01_r8
         z(i,k) = zm(i,k) + zs(i)
         zf(i,k) = zi(i,k) + zs(i)
      end do
   end do
!
  ! calculate surface inhomogeity
  ! =================================

  qpert(:)   = qpert2(:) *0.0
  tpert(:)   = tpert2(:) *0.0
  calday = get_curr_calday()
  call zenith (calday, lat, lon, coszrs, ncol)

  ! q frontal land only
  if(zmh_deep_qpertf)then
    qpert(:) = qpert(:) + gradq(:)*zmh_deep_xyq*cos(lat(:))*max(coszrs(:),0.0)* &
             (200.0/max(sgh30(:,lchnk),200.)) * min(landfrac(:)*100.,1.)
  endif  
  ! q terrain
  if(zmh_deep_qperth)then
    do i= 1,ncol
       dum1     = z(i,pver-2) - z(i,pver) 
       qpert(i) = qpert(i) + abs(qh(i,pver,1) - qh(i,pver-2,1))/dum1 &
           * min(sgh30(i,lchnk),100._r8) *0.0
     enddo
  endif  

  ! t frontal
  if(zmh_deep_tpertf)then
    tpert(:) = tpert(:) + gradt(:)*zmh_deep_xyt*cos(lat(:))*max(coszrs(:),0.0)* &
                     (200.0/max(sgh30(:,lchnk),200.)) * min(landfrac(:)*100.,1.)
  endif

  ! t breeze
  if(zmh_deep_tpertb)then
    do i=1,ncol
      if(abs(landfrac(i)-0.5)  .le.  0.45 .or. sgh30(i,lchnk) .gt. 50._r8)then
        tpert(i) = tpert(i) + gradt(i)*zmh_deep_xyt*max(coszrs(i),0.)*cos(lat(i))* &
             (200.0/max(sgh30(i,lchnk),200.))
      endif
    enddo
  endif 

  do i=1,ncol
    tpert(i) = min(tpert(i), zmh_deep_tpertmax)
    qpert(i) = min(qpert(i), zmh_deep_qpertmax)
  enddo

   do k = 1,pver
      do i = 1,ncol
         q(i,k) = qh(i,k,1)
         s(i,k) = t(i,k) + (grav/cpres)*z(i,k)
         tp(i,k)=0.0_r8
         shat(i,k) = s(i,k)
         qhat(i,k) = q(i,k)
!zmh
         dz(i,k) = zf(i,k) - zf(i,k+1)
      end do
   end do

!===========================
  if(zmh_deep_tpertt)then 
    do i=1,ncol
     do k=pver,2,-1
      if(zf(i,k+1) .gt. pblh(i)) exit !cycle
       s(i,k)       = s(i,k)  + tpert(i) 
      enddo
    enddo
   tpert(:) = tpert2(:)
  endif
!
  if(zmh_deep_qpertt)then 
    do i=1,ncol
     do k=pver,2,-1
      if(zf(i,k+1) .gt. pblh(i)) exit !cycle
       q(i,k)       = q(i,k)  + qpert(i) 
      enddo
    enddo
   qpert(:) = qpert2(:)
  endif
!
!===========================

   do i = 1,ncol
      capeg(i) = 0._r8
      lclg(i) = 1
      lelg(i) = pver
      maxg(i) = 1
      tlg(i) = 400._r8
      dsubcldg(i) = 0._r8
   end do

! ----------------------   
!  zmh relative humidity
   call aqsat (t    ,pap  ,est    ,qsat    ,pcols   ,  &   !check Pa or mb
              ncol ,pver  ,1       ,pver    )

   rh(:,:) = q(:,:)/qsat(:,:)
   hmn(:,:) = cp*t(:,:) + grav*z(:,:) + rl*q(:,:)
   rho(:,:) = p(:,:)/2.876_r8/t(:,:)

  do i = 1, ncol
    k750(i) = pver
    k950(i) = pver
    k100(i) = pver
    kpbl(i) = pver
    maxi(i)   = pver
    jd(i)   = pver

    precw(i) = 0._r8
    precws(i) = 0._r8

    ! pbl top level 
    do k= pver,msg,-1
     if ((zf(i,k)-pblh(i)) .ge. 0.0_r8) then
         kpbl(i) = k
         exit 
     endif
    end do
    pblt(i) = dble(kpbl(i))  !old inheritage

    ! level below which precipitable water is used in trigger and at which omega 
    ! is used in trigger
    do k= pver,msg,-1
        if(p(i,k)/pf(i,pver+1) .le. 0.8 )  then   
          k750(i) = k
          exit 
        endif
    enddo    
    do k= k750(i), pver
        precw(i)  = precw(i) + q(i,k)*dpp(i,k)/gravit  !  check dpp in kg
        precws(i) = precws(i) + qsat(i,k)*dpp(i,k)/gravit  !  check dpp in kg
    enddo    

   ! level of subcloud layer 
    do k= pver,msg,-1
        if(p(i,k)/pf(i,pver+1) .le. 0.990 )  then   
          k950(i) = k         
          exit 
        endif
    enddo    

    ! highest level parcel can be launched
    do k= pver,msg,-1
       if(p(i,k)<(pf(i,pver+1)-200.))then   ! in mb unit!!
          k100(i) = k         
          exit 
       endif 
    enddo
    k1 = min(kpbl(i),k100(i))

    ! level of maximum moist static energy below level k1
    k  = get_maxk( hmn(i,:),pver, k1, pver)
    !maxi(i) = k
    maxi(i) = k-1

    ! find the top (maxi1) and bottom (maxi2) levels of the launching layer,
    ! set here to 10 mb within the maxi level. 
    ! This is due to high vertical resolution near the surface
    k1 = pver
    do k= maxi(i)+1,pver
     if(p(i,k) - pf(i,maxi(i)+1) .ge. 10.)then
      k1 = k
      exit
     endif
    enddo 
    if(k1==pver)then
      do k = pver,limcnv,-1
        if(pf(i,pver+1) - p(i,k) .ge. 10.)then
          k1 = k
          exit
         endif
       enddo
     endif 
    maxi1(i) = min(k1,maxi(i)) 
    maxi2(i) = max(k1,maxi(i)+1)  !190606 



    ! find level of downdraft
    k1   = get_mink(hmn(i,:),pver, msg, maxi(i)-2)
!zmh added bounds 2019-04-24
    k1   = min(pver-5,k1)
    k1   = max(msg+2,k1)
! test 190602
   k1  = k1-1


    jd(i)    = k1
    do k= k1-1,msg,-1
        if(p(i,k1)-p(i,k) .gt. 150. )  then   
          jd(i) = k-1         
          exit 
        endif
    enddo    

    ! get bottom level (jd1) of downdraft launching layer, set to 50 mb below jd
    jd1(i)    = k1
    do k= k1+1,pver
        if(p(i,k)-p(i,k1) .gt. 50. )  then   
          jd1(i) = k+1         
          exit 
        endif
    enddo    
! test 190602
  jd1(i) = k1-1

   ! level (jd2) for downdraft organized detrainment, set to 25mb above surface
    do k= pver,msg,-1
        if(pf(i,pver+1)-p(i,k) .ge. 100. )  then   
          jd2(i) = k         
          exit 
        endif
    enddo    
    jd2(i) = max(jd(i)+1,jd2(i))

! calculate sub-cloud layer pressure "thickness" for use in
! closure and tendency routines.
!
   do k = k950(i),pver
      dsubcld(i) = dsubcld(i) + dpp(i,k)*0.01
   end do

  enddo  !i


! ----------------------   

! initial guess of entrainment
  ! =================================

     ent_turb(:,:) = zmh_det_turb0 *(1.+0.2*min(rh(:,:),1._r8))  ! p dependence?
     det_turb(:,:) = zmh_det_turb0 *(1.+0.2*min(rh(:,:),1._r8))

     ent_tot(:,:) =  -zmh_dmpdz 
     det_tot(:,:) =  0.0_r8
     w_init(:)    = 1.0_r8

     dzt(:)        = zmh_dzt  !0.10 !normalized thickness of updraft detrainment layer
     dzb(:)        = zmh_dzb  !0.90 !normalized thickness of updraft entraiment layer

  DO ITER = 1, NITER      !iteration for CAPE and vertical velocity

!print*,'maxi=',iter,maxi    ,maxi1    ,maxi2
      ! for initial bouyancy
      call zyx1_buoyan_dilute(lchnk   ,ncol    , limcnv,&
                  q       ,t       ,p       ,z       ,pf       , &
                  tp      ,qstp    ,tl      ,rl      ,cape     , &
                  pblt    ,lcl     ,lel     ,maxi    ,maxi1    ,maxi2  , &
                  rgas    ,grav    ,cpres   ,msg     , tpert   , &
                  qpert   ,ent_tot ,buoy_mid,         cin  )

     ! for vertical velocity and convection top          
     call w_dynamics (ncol ,limcnv ,z ,dz ,buoy_mid ,w_init ,ent_tot,  &
                  maxi ,maxi1, lel ,lel2, lelb, w_up, eps00 )

     do i = 1,ncol
     do k = lel(i),lelb(i)
       w11 = (z(i,k)-z(i,lelb(i)))/max(1.0,z(i,lel(i))-z(i,lelb(i)))
       ent_turb(i,k) = zmh_ent_turb0 * max(w11,0.0) +  zmh_det_turb0 * (1.+0.2*min(rh(i,k),1._r8)) 
       det_turb(i,k) = zmh_ent_turb0 * max(w11,0.0) +  zmh_det_turb0 * (1.+0.2*min(rh(i,k),1._r8)) 
     enddo
     enddo

    if(ITER .ne. NITER)then

     !for total entrainment and detainment rates  , 1 for updraft 
     call ent_dynamics (ncol ,zf  ,dz  ,maxi, maxi1, maxi2,lel, lel, &!this lel here is on purpose
                  dzb  ,dzt , ent_turb  ,det_turb, eps00, ent_tot ,det_tot , 1)

     !update detrainment layer thickness to be between max buoyancy and zero w
!     do i = 1,ncol
!         dzt(i) = (z(i,lel(i)) - z(i,lelb(i)) ) /(z(i,lel(i)) - z(i,maxi2(i)) ) 
!         dzt(i) = max(dzt(i),0.05_r8)
!     enddo

    else  

       call ent_dynamics (ncol ,zf  ,dz  ,maxi, maxi1, maxi2,lel, lel2,&
                 dzb  ,dzt  ,  ent_turb  ,det_turb,eps00, ent_tot ,det_tot , 1)

    endif  


!k=1
!if(k<0)then
!write(*,*)''
!write(*,*)'ITER =',ITER
!write(*,*)'cape, cin,=',cape, cin
!write(*,*)'eps00',eps00
!write(*,*)'mypblh ',pblh
!write(*,*)'lel,lel2, maxi,maxi1,maxi2',lel,lel2, maxi,maxi1,maxi2
!write(*,*)'buoy',buoy_mid*t/9.8
!write(*,*)'ent', ent_tot
!write(*,*)'det', det_tot
!write(*,*)'w_up', w_up
!endif

  ENDDO !iter

  !limit deep convection to be thicker that 200 mb
  do i=1,ncol
   jt(i) = lel(i)
   if( (p(i,maxi1(i)) - p(i,lel(i))) < 200. .or. (maxi1(i) - lel2(i) < 3))then
       cape(i) = 0.0
   endif
  enddo

  lengath = 0

! =================================
! for trigger ---------------------

 flag(:) = .false.
 ramp(:) = 1.0

 select case (deep_trigger_scheme)


case('zmh_triggers') !additional limit

     deep_triggers(:,:) = .true.
     bfls_cape(:) = 0.0_r8


     if(do_dcape)then 
!    -------------------------------------------------------
      call zyx1_buoyan_dilute(lchnk   ,ncol    , limcnv, &
                  bfls_q       ,bfls_t       ,p       ,z       ,pf       , &
                  tp      ,qstp    ,tl      ,rl      ,bfls_cape     , &
                  pblt    ,lclg     ,lelg     ,maxi    ,maxi1    ,maxi2  , &  !the g is for temporary storage
                  rgas    ,grav    ,cpres   ,msg     , tpert,    &
                  qpert   ,ent_tot ,bfls_buoy_mid   , bfls_cin   )

!    -------------------------------------------------------
     end if


    do i=1,ncol
        w11 = tpert2(i) + gradt(i)*zmh_deep_xyt*cos(lat(i))* &
                         (500.0/max(sgh30(i,lchnk),500.))
        w12 = qpert2(i) + gradq(i)*zmh_deep_xyq*cos(lat(i))* &
                         (500.0/max(sgh30(i,lchnk),500.))


       deep_triggers(i,1) = (cape(i)  > capelmt(i) )

       deep_triggers(i,2) = (cape(i)  > bfls_cape(i)*0.5_r8)  

      w13 = 2.         !2 mb/hr  !range expanded 
      w13 = w13/36.    ! Pa/s
      if(landfrac(i) .gt. 0.0) then
        w13 = 0.5/36. 
      endif

      deep_triggers(i,3) = (omega(i,k750(i)) < -w13*(1.+rand_numw(i))*0.5 )             

       ! gradually ramp up with omega
       !if(deep_triggers(i,3))then
       ! w11 = zmh_ramp(-omega(i,k750(i)),0._r8, w13)
       ! ramp(i) = ramp(i)*w11
       !endif


      w14 = 0.70 !5 !0.6
!      if(landfrac(i) .eq. 0.0) then
!        w14 = 0.75
!      endif

      deep_triggers(i,4) = (precw(i) > w14*precws(i)) !precw below k750
      !if(deep_triggers(i,4))then
      !  ramp(i) = ramp(i)*zmh_ramp(precw(i)/precws(i),w14,0.10_r8)
      !endif


!ramp(i) = 1.0
!deep_triggers(i,3) = .true.
!deep_triggers(i,4) = .true.

       ! other optional triggers
!       deep_triggers(i,5) = (pblh(i)-zi(i,pverp)          > 500._r8 )
!       deep_triggers(i,6) = (rh(i,kpbl(i)-2)  > 0.75_r8 )       !rh ,kpbl(i)
!       deep_triggers(i,7) = (tke(i)   > cin(i)  )

!323
!       if(landfrac(i) .le. 0.01)then
!        deep_triggers(i,5) = (vort3(i)*sin(lat(i)) .ge. 0.0)
!       endif

!      if(deep_triggers(i,8))then 
!         w11 = zmh_ramp(vort3(i)*sin(lat(i)),-1.e-7_r8,0.2e-5_r8)
!         w12 = zmh_ramp(-omega(i,k750(i)),0._r8,1.e-5_r8)
!         ramp(i) = ramp(i)*0.5*(w11 + w12)
!       endif
        
      flag(i) = deep_triggers(i,1)
      do iflag = 2, 10
         flag(i) = flag(i) .and. deep_triggers(i,iflag)
      enddo    

      if(flag(i)) then
         lengath = lengath + 1
         index(lengath) = i
       end if

!k=1
!if(k<0)then
!write(*,*)'flag', flag
!write(*,*)'trigg', deep_triggers
!write(*,*)'cape(i)  > capelmt',cape(i),  capelmt
!write(*,*)'cape(i)  > bfls_cape(i)',cape(i)  , bfls_cape(i)
!write(*,*)'cin(i)   < bfls_cin',cin(i)   , bfls_cin(i)
!write(*,*)'omega(i,kpbl(i))',omega(i,kpbl(i)),kpbl(i)
!write(*,*)'pblh(i)',pblh(i)
!write(*,*)'rh(i,kpbl(i)-2)',rh(i,kpbl(i)-2)
!write(*,*)'tke(i)   > cin(i)',tke(i)   , cin(i)
!write(*,*)'precw(i)',precw(i)
!write(*,*)'tpert(i)',tpert(i)
!write(*,*)'qpert(i)',qpert(i)
!write(*,*)'t       ',t
!write(*,*)'bfls_t  ',bfls_t
!rite(*,*)'q       ',q
!rite(*,*)'bfls_q  ',bfls_q
!    write(*,*)'bfls_t',bfls_t
!    write(*,*)'bfls_q',bfls_q
!    write(*,*)'pbl',pblh
!    write(*,*)'tpert',tpert
!    write(*,*)'qpert',qpert
!    write(*,*)'tke',tke

!endif
    end do

   case default !cam  trigger
     do i=1,ncol
      if ( cape(i) > capelmt(i)) then 
         flag(i) = .true.
         lengath = lengath + 1
         index(lengath) = i
       endif
      enddo
  end select  !trigger done
! -------------------------------------------------------

  if(deep_closure_scheme=='prognostic')then
   lengath = ncol
   do i=1,ncol
         index(i) = i
   enddo
  endif

   if (lengath.eq.0) return
   do ii=1,lengath
      i=index(ii)
      ideep(ii)=i
   end do
!
! obtain gathered arrays necessary for ensuing calculations.
!
   do k = 1,pver
      do i = 1,lengath
         dp(i,k) = 0.01_r8*dpp(ideep(i),k)
         qg(i,k) = q(ideep(i),k)
         tg(i,k) = t(ideep(i),k)
         pg(i,k) = p(ideep(i),k)
         zg(i,k) = z(ideep(i),k)
         sg(i,k) = s(ideep(i),k)
         tpg(i,k) = tp(ideep(i),k)
         zfg(i,k) = zf(ideep(i),k)
         qstpg(i,k) = qstp(ideep(i),k)
         ug(i,k) = 0._r8
         vg(i,k) = 0._r8
!zmh
         ent_totg(i,k) = ent_tot(ideep(i),k)
         det_totg(i,k) = det_tot(ideep(i),k)
      end do
   end do
!
   do i = 1,lengath
      zfg(i,pver+1) = zf(ideep(i),pver+1)
   end do
   do i = 1,lengath
      !flagg(i) = flag(ideep(i))
      flagg(i) = deep_triggers(ideep(i),1)
      capeg(i) = cape(ideep(i))
      lclg(i) = lcl(ideep(i))
      lelg(i) = lel(ideep(i))
      maxg(i) = maxi(ideep(i))
      tlg(i) = tl(ideep(i))
      landfracg(i) = landfrac(ideep(i))
!zmh
      rampg(i) = ramp(ideep(i))
      maxi1g(i) = maxi1(ideep(i))
      maxi2g(i) = maxi2(ideep(i))
      jdg(i) = jd(ideep(i))
      jd1g(i) = jd1(ideep(i))
      jd2g(i) = jd2(ideep(i))
      jtg(i) = jt(ideep(i))
      lel2g(i) = lel2(ideep(i))
      lelbg(i) = lelb(ideep(i))
      eps00g(i) = eps00(ideep(i))
      dsubcldg(i) = dsubcld(ideep(i))
      k950g(i) = k950(ideep(i))
      latg(i ) =lat(ideep(i))
      long(i ) =lon(ideep(i))
   end do
!
! calculate sub-cloud layer pressure "thickness" for use in
! closure and tendency routines.
!
!   do k = msg + 1,pver
!      do i = 1,lengath
!         if (k >= maxg(i)) then
!            dsubcld(i) = dsubcld(i) + dp(i,k)
!         end if
!      end do
!   end do

!
! define array of factors (alpha) which defines interfacial
! values, as well as interfacial values for (q,s) used in
! subsequent routines.
!
   do k = msg + 2,pver
      do i = 1,lengath
!            alpha(i,k) = 0.5
         sdifr = 0._r8
         qdifr = 0._r8
         if (sg(i,k) > 0._r8 .or. sg(i,k-1) > 0._r8) &
            sdifr = abs((sg(i,k)-sg(i,k-1))/max(sg(i,k-1),sg(i,k)))
         if (qg(i,k) > 0._r8 .or. qg(i,k-1) > 0._r8) &
            qdifr = abs((qg(i,k)-qg(i,k-1))/max(qg(i,k-1),qg(i,k)))
         if (sdifr > 1.E-6_r8) then
            shat(i,k) = log(sg(i,k-1)/sg(i,k))*sg(i,k-1)*sg(i,k)/(sg(i,k-1)-sg(i,k))
         else
            shat(i,k) = 0.5_r8* (sg(i,k)+sg(i,k-1))
         end if
         if (qdifr > 1.E-6_r8) then
            qhat(i,k) = log(qg(i,k-1)/qg(i,k))*qg(i,k-1)*qg(i,k)/(qg(i,k-1)-qg(i,k))
         else
            qhat(i,k) = 0.5_r8* (qg(i,k)+qg(i,k-1))
         end if
      end do
   end do
!
! obtain cloud properties.
!
   ! indices of the updraft  : maxi2g >= maxg >=maxi1g > lel2g >=lelg
   ! indcies of the downdraft: jdg  <= jdg1   < jd2g  <= pver
   call zmh_plume_prp(lchnk   , &
               qg      ,tg      ,ug      ,vg      ,pg      , &
               zg      ,sg      ,mu      ,eu      ,du      , &
               md      ,ed      ,sd      ,qd      ,mc      , &
               qu      ,su      ,zfg     ,qs      ,hmn     , &
               hsat    ,shat    ,qlg                       , &
               cmeg    ,maxg    ,maxi1g  ,maxi2g  ,lelg    , &
               lel2g   ,jdg     ,jd1g    ,jd2g    ,jlcl    , &
               rl      ,lengath                            , &
               rgas    ,grav    ,cpres   ,msg              , &
               pflxg   ,evpg    ,cug     ,rprdg   ,limcnv  , &
               landfracg        ,tpert   ,qpert   ,buoy_mid, &
               ent_totg,det_totg,eps00g  ,capeg   ,cin, flagg, latg, long) 
              
!
! convert detrainment from units of "1/m" to "1/mb".
!
   do k = msg + 1,pver
      do i = 1,lengath
         du   (i,k) = du   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         eu   (i,k) = eu   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         ed   (i,k) = ed   (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cug  (i,k) = cug  (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         cmeg (i,k) = cmeg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         rprdg(i,k) = rprdg(i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
         evpg (i,k) = evpg (i,k)* (zfg(i,k)-zfg(i,k+1))/dp(i,k)
      end do
   end do

   call closure(lchnk   , &
                qg      ,tg      ,pg      ,zg      ,sg      , &
                tpg     ,qs      ,qu      ,su      ,mc      , &
                du      ,mu      ,md      ,qd      ,sd      , &
                qhat    ,shat    ,dp      ,qstpg   ,zfg     , &
                qlg     ,dsubcldg ,mb      ,capeg   ,tlg     , &
                lclg    ,lelg    ,jt      ,maxg    ,1       , &
                lengath ,rgas    ,grav    ,cpres   ,rl      , &
                msg     ,capelmt ,flagg   )
!

! -------------------------
! -------------------------------------------------------

! closure

  if(deep_trigger_scheme .eq. 'zmh_triggers')then
    do i=1,lengath
       mb(i) = mb(i)*rampg(i)   !gradual transition in trigger
    enddo
  endif

  if(deep_closure_scheme .eq. 'prognostic')then           
    do i=1,lengath
       ezmh_deep_tau1 =  exp(- delt/zmh_deep_tau1)  !time toward mb
       ezmh_deep_tau2 =  exp(- delt/zmh_deep_tau2)

      if( .not. deep_triggers(i,1)) then
       mb(i) = 0.
      else if(flag(i))then
       mb(i) = massflxbase_p(i,1)*(ezmh_deep_tau1+ezmh_deep_tau2)*0.5   &
               + (1.-ezmh_deep_tau1)*mb(i) 
      else
       mb(i) = massflxbase_p(i,1)*(ezmh_deep_tau1+ezmh_deep_tau2)*0.5  
      endif
    enddo
  endif 

  do i=1,lengath
      mb(i) = mb(i)*zmh_mb ! a tunable parameter
      if(mb(i) .le. 1.0e-4_r8)then  !suppress very weak deep convection
          mb(i) = 0._r8
      endif
  enddo
! -------------------------------------------------------
!
! limit cloud base mass flux to theoretical upper bound.
!
   do i=1,lengath
      mumax(i) = 1.0e-10
   end do
   do k=msg + 2,pver
      do i=1,lengath
        mumax(i) = max(mumax(i), mu(i,k)/dp(i,k))
      end do
   end do

! -------------------------
!zmh relaxed limit

 if(zmh_deep_nolimit) then  !substepping to allow large mass flux
   do i=1,lengath
     nt_stepping(i) = int(mumax(i)*delt/0.5_r8)+1
     nt_stepping(i) = min(nt_stepping(i),2)                     !temporily set 2
     mb(i) = min(mb(i),nt_stepping(i)*0.5_r8/(delt*mumax(i)))   
   end do
 else
   do i=1,lengath
      if (mumax(i) > 0._r8) then
         mb(i) = min(mb(i),0.5_r8/(delt*mumax(i)))
      else
         mb(i) = 0._r8
      endif
      nt_stepping(i) = 1
   end do
  endif

! ------------------------- 
   do i = 1, lengath
    massflxbase_p(ideep(i),1) = mb(i)
   enddo

! -------------------------

   do k=msg+1,pver
      do i=1,lengath
         mu   (i,k)  = mu   (i,k)*mb(i)
         md   (i,k)  = md   (i,k)*mb(i)
         mc   (i,k)  = mc   (i,k)*mb(i)
         du   (i,k)  = du   (i,k)*mb(i)
         eu   (i,k)  = eu   (i,k)*mb(i)
         ed   (i,k)  = ed   (i,k)*mb(i)
         cmeg (i,k)  = cmeg (i,k)*mb(i)
         rprdg(i,k)  = rprdg(i,k)*mb(i)
         cug  (i,k)  = cug  (i,k)*mb(i)
         evpg (i,k)  = evpg (i,k)*mb(i)
         pflxg(i,k+1)= pflxg(i,k+1)*mb(i)*100._r8/grav
      end do
   end do
!
! compute temperature and moisture changes due to convection.

   do k=msg+1,pver
      do i=1,lengath
         slflx(ideep(i),k) =   ( mu(i,k)* (su(i,k)-shat(i,k)) &
                         + md(i,k)* (sd(i,k)-shat(i,k)) )*cpair*100._r8/grav
         qtflx(ideep(i),k) =   ( mu(i,k)* (qu(i,k)-qhat(i,k)) &
                         + md(i,k)* (qd(i,k)-qhat(i,k)) )*latvap*100._r8/grav
!   print*,'myqtflx',k,i,ideep(i),qtflx(ideep(i),k)
      end do
   end do

! cloud amount 

   do i=1,lengath
    ii = ideep(i)   
    dum0 = max(zm(ii,lelg(i)) - zm(ii,maxg(i)), 1.0_r8)  !depth

    do k= lelg(i),maxg(i)

         dum1 = (zm(ii,k)-zm(ii,maxg(i)) ) / dum0 !normalized height
         dum1 = 1.0 - zmh_ramp(dum1,0._r8,0.3_r8)  !enhancement near lower plume
         dum1 = 1.0 + zmh_dfrcz*dum1 

         dum2 = dum1*zmh_dfrc0*mu(i,k)*10.24_r8/max(w_up(ii,k),0.1_r8)/rho(ii,k) & 
                + du(i,k)*dz(ii,k)*delt / dp(i,k) * zmh_dfrcd 
                                       !10.24 = 100./g
         deepfrc(ii,k) = min(dum2, 1.0_r8)
      end do !k

         k    = min(lelg(i)+1,pver)
         dum1 = mu(i,k)*delt / dp(i,k)
         dum2 = max(deepfrc(ii,k), zmh_dfrcd*dum1)  
         deepfrc(ii,k) = min(dum2, 1.0_r8)

   end do

!
   call q1q2_zmh(lchnk   , &
                 dqdt    ,dsdt    ,qg      ,qs      ,qu      , &
                 su      ,du      ,qhat    ,shat    ,dp      , &
                 mu      ,md      ,sd      ,qd      ,qlg     , &
                 dsubcldg,jt      ,maxg    ,1       ,lengath , &
                 cpres   ,rl      ,msg     ,          &
                 dlg     ,evpg    ,cug     ,          &
                 delt, nt_stepping,k950g ,flagg)
!
! gather back temperature and mixing ratio.
!
   do k = msg + 1,pver
!DIR$ CONCURRENT
      do i = 1,lengath
!
! q is updated to compute net precip.
!
         q(ideep(i),k) = qh(ideep(i),k,1) + 2._r8*delt*dqdt(i,k)
         qtnd(ideep(i),k) = dqdt (i,k)
         cme (ideep(i),k) = cmeg (i,k)
         rprd(ideep(i),k) = rprdg(i,k)
         zdu (ideep(i),k) = du   (i,k)
         mcon(ideep(i),k) = mc   (i,k)
         heat(ideep(i),k) = dsdt (i,k)*cpres
         dlf (ideep(i),k) = dlg  (i,k)
         pflx(ideep(i),k) = pflxg(i,k)
         ql  (ideep(i),k) = qlg  (i,k)
      end do
   end do

!do i=1,lengath
!if(abs(latg(i)*180./3.1416-0.0) .le. 0.75 .and. abs(long(i)*180./3.1416-70.0) .le. 0.75)then
!write(*,*)'mydlf',i,dlg(i,:)*8.64e7
!write(*,*)'mydqt',i,dqdt(i,:)*8.64e7
!endif
!enddo


!
!DIR$ CONCURRENT
   do i = 1,lengath
      cape(ideep(i)) = capeg(i)   !zmh
      jctop(ideep(i)) = jt(i)
!++bee
      jcbot(ideep(i)) = maxg(i)
!--bee
      pflx(ideep(i),pverp) = pflxg(i,pverp)
   end do

! Compute precip by integrating change in water vapor minus detrained cloud water
   do k = pver,msg + 1,-1
      do i = 1,ncol
         prec(i) = prec(i) - dpp(i,k)* (q(i,k)-qh(i,k,1)) - dpp(i,k)*dlf(i,k)*2*delt
      end do
   end do

! obtain final precipitation rate in m/s.
   do i = 1,ncol
      prec(i) = rgrav*max(prec(i),0._r8)/ (2._r8*delt)/1000._r8
   end do

! Compute reserved liquid (not yet in cldliq) for energy integrals.
! Treat rliq as flux out bottom, to be added back later.
   do k = 1, pver
      do i = 1, ncol
         rliq(i) = rliq(i) + dlf(i,k)*dpp(i,k)/gravit
      end do
   end do
   rliq(:ncol) = rliq(:ncol) /1000._r8


   return
end subroutine zyx1_conv_tend

! ===============================================================
! ===============================================================
subroutine zmh_plume_prp(lchnk   , &
                  q       ,t       ,u       ,v       ,p       , &
                  z       ,s       ,mu      ,eu      ,du      , &
                  md      ,ed      ,sd      ,qd      ,mc      , &
                  qu      ,su      ,zf      ,qst     ,hmn     , &
                  hsat    ,shat    ,ql      ,cmeg             , &
                  jb      ,jb1     ,jb2     ,jt      ,jt2     , &
                  jd      ,jd1     ,jd2     ,jlcl             , &
                  rl      ,il2g    , &
                  rd      ,grav    ,cp      ,msg     , &
                  pflx    ,evp     ,cu      ,rprd    ,limcnv  , &
                  landfrac,tpert   ,qpert   ,buoy_mid         , &
                  ent_tot, det_tot, eps00   ,capeout, cinout,flag, lat, lon) 

!----------------------------------------------------------------------- 
! 
! Plume model
! Author:by Minghua Zhang 2017-02-01
! original code based on CAM
! 
! updraft:  jb2 >= jb >=jb1 launching layer; jt2>jt detrainment trans layer
! downdraft:jd <= jd1       launching layer; jd2<=pver detainment trans layer              
! documentation is yet to be written
!-----------------------------------------------------------------------
!!   use buoysort, only : cal_buoysort

   implicit none

!------------------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                  ! chunk identifier

   real(r8), intent(in) :: q(pcols,pver)         ! spec. humidity of env
   real(r8), intent(in) :: t(pcols,pver)         ! temp of env
   real(r8), intent(in) :: p(pcols,pver)         ! pressure of env
   real(r8), intent(in) :: z(pcols,pver)         ! height of env
   real(r8), intent(in) :: s(pcols,pver)         ! normalized dry static energy of env
   real(r8), intent(in) :: zf(pcols,pverp)       ! height of interfaces
   real(r8), intent(in) :: u(pcols,pver)         ! zonal velocity of env
   real(r8), intent(in) :: v(pcols,pver)         ! merid. velocity of env

   logical,  intent(in) :: flag(pcols)
   real(r8), intent(in) :: landfrac(pcols) ! RBN Landfrac
   real(r8), intent(in) :: tpert(pcols)   !zmh
   real(r8), intent(in) :: qpert(pcols)   !zmh

   integer, intent(inout) :: jt2(pcols)              ! updraft base level
   integer, intent(in) :: jb(pcols)              ! updraft base level
   integer, intent(in) :: jb1(pcols)              ! updraft base level
   integer, intent(in) :: jb2(pcols)              ! updraft base level
   integer, intent(in) :: jd1(pcols)              ! updraft base level
   integer, intent(in) :: jd2(pcols)              ! updraft base level
   integer, intent(inout) :: jt(pcols)              ! updraft plume top
   integer, intent(out) :: jlcl(pcols)            ! updraft lifting cond level
   integer, intent(inout) :: jd(pcols)              ! level of downdraft
   integer, intent(in) :: limcnv                 ! convection limiting level
   integer, intent(in) :: il2g                   !CORE GROUP REMOVE
   integer, intent(in) :: msg                    ! missing moisture vals (always 0)
   real(r8), intent(in) :: rl                    ! latent heat of vap
   real(r8), intent(in) :: shat(pcols,pver)      ! interface values of dry stat energy
   real(r8), intent(in) ::  eps00(pcols)  !organized entrainment coef
   real(r8), intent(in) ::  det_tot(pcols,pver)
   real(r8), intent(in) ::  ent_tot(pcols,pver)

   real(r8), intent(in) :: lat(pcols)
   real(r8), intent(in) :: lon(pcols)

!output
!
   real(r8), intent(out) :: rprd(pcols,pver)     ! rate of production of precip at that layer
   real(r8), intent(out) :: du(pcols,pver)       ! detrainement rate of updraft
   real(r8), intent(out) :: ed(pcols,pver)       ! entrainment rate of downdraft
   real(r8), intent(out) :: eu(pcols,pver)       ! entrainment rate of updraft
   real(r8), intent(out) :: hmn(pcols,pver)      ! moist stat energy of env
   real(r8), intent(out) :: hsat(pcols,pver)     ! sat moist stat energy of env
   real(r8), intent(out) :: mc(pcols,pver)       ! net mass flux
   real(r8), intent(out) :: md(pcols,pver)       ! downdraft mass flux
   real(r8), intent(out) :: mu(pcols,pver)       ! updraft mass flux
   real(r8), intent(out) :: pflx(pcols,pverp)    ! precipitation flux thru layer
   real(r8), intent(out) :: qd(pcols,pver)       ! spec humidity of downdraft
   real(r8), intent(out) :: ql(pcols,pver)       ! liq water of updraft
   real(r8), intent(out) :: qst(pcols,pver)      ! saturation spec humidity of env.
   real(r8), intent(out) :: qu(pcols,pver)       ! spec hum of updraft
   real(r8), intent(out) :: sd(pcols,pver)       ! normalized dry stat energy of downdraft
   real(r8), intent(out) :: su(pcols,pver)       ! normalized dry stat energy of updraft
   real(r8),intent(out)  :: capeout(pcols)
   real(r8),intent(out)  :: cinout(pcols)

   real(r8) :: dd(pcols,pver)       !detrainment of downdraft, not output for now

   real(r8) rd                   ! gas constant for dry air
   real(r8) grav                 ! gravity
   real(r8) cp                   ! heat capacity of dry air
   real(r8) latice

!
! Local workspace

   integer NFREEZE     !number of itegration for freezing calc

   real(r8) gamma(pcols,pver)
   real(r8) dz(pcols,pver)
   real(r8) iprm(pcols,pver)
   real(r8) hu(pcols,pver)
   real(r8) hd(pcols,pver)
   real(r8) eps(pcols,pver)
   real(r8) f(pcols,pver)

   real(r8) qsthat(pcols,pver)
   real(r8) hsthat(pcols,pver)
   real(r8) gamhat(pcols,pver)
   real(r8) cu(pcols,pver)
   real(r8) evp(pcols,pver)
   real(r8) cmeg(pcols,pver)
   real(r8) qds(pcols,pver)


   real(r8) tv(pcols,pver)
   real(r8) tv_up(pcols,pver)
   real(r8) tv_dn(pcols,pver)
   real(r8) buoy_up(pcols,pver)
   real(r8) buoy_mid(pcols,pver)
   real(r8) rho(pcols,pver)
   real(r8) ent_turb_dn(pcols,pver)
   real(r8) det_turb_dn(pcols,pver)
   real(r8) det_tot_dn(pcols,pver)
   real(r8) ent_tot_dn(pcols,pver)
   real(r8) ramp_t1(pcols,pver)
   real(r8) ramp_t2(pcols,pver)

   real(r8) freeze_rate(pcols,pver)

   real(r8) cldradinit(pcols)
   real(r8) cldh(pcols)
   real(r8) c0mask(pcols)

   real(r8) dmpdz, dum1 , dum0, dume, dum2, dum3, dumb
   real(r8) dzb(pcols),dzt(pcols)
   real(r8) hmnf(pcols,pver)
   real(r8) qf(pcols,pver)
   real(r8) sf(pcols,pver)
   integer jd3(pcols) 

!-----------------------

   real(r8) rain_z0(pcols)
   real(r8) rain_zp(pcols)
   integer kount,iter

   real(r8) hmin(pcols)
   real(r8) ftemp(pcols)
   real(r8) eps0(pcols)
   real(r8) rmue(pcols)
   real(r8) zuef(pcols)
   real(r8) zdef(pcols)
   real(r8) epsm(pcols)
   real(r8) ratmjb(pcols)
   real(r8) est(pcols)
   real(r8) alfa(pcols)
   real(r8) ql1,ql2
   real(r8) tu
   real(r8) estu
   real(r8) qstu

   real(r8) small
   real(r8) mdt

   integer i,k

   logical doit(pcols)
   logical done(pcols)
!
!  zmh local 

!------------------------------------------------------------------------------
!
   NFREEZE = 1 !1 is no consideration of freeze 

   latice = 3.34e5_r8

   do i = 1,il2g
!!
! if(abs(lat(i)*180./3.1416-0.0) .le. 0.75 .and. abs(lon(i)*180./3.1416-70.0) .le. 0.75)then  
!  print*,'jtjt2 ',i,jt(i),jt2(i)
! endif
      c0mask(i)  = zmh_c0_ocn * (1._r8-landfrac(i)) +   zmh_c0_lnd * landfrac(i)

      ftemp(i) = 0._r8

      rain_z0(i) = zmh_rain_z0*(1.0_r8 + 0.5*landfrac(i))  !0.5 tunable 
      rain_zp(i) = zmh_rain_zp*(1.0_r8 + 0.5*landfrac(i))
   end do
!
   do k = 1,pver
      do i = 1,il2g
       dz(i,k) = zf(i,k) - zf(i,k+1)
       ramp_t1(i,k) = 1.0_r8 - zmh_ramp(t(i,k)-243.15_r8, 0._r8,25._r8) !freeze
       !ramp_t2(i,k) = 1.0_r8 - zmh_ramp(t(i,k)-243.15_r8, 0._r8,25._r8) !conversion
       ramp_t2(i,k) = 1.0_r8 - zmh_ramp(t(i,k)-253.15_r8, 0._r8,5._r8) !conversion
       ramp_t2(i,k) = sqrt(ramp_t2(i,k))
      end do
   end do

!
! initialize many output and work variables to zero
!
   pflx(:il2g,1) = 0

   do k = 1,pver
      do i = 1,il2g
         mu(i,k) = 0._r8
         f(i,k) = 0._r8
         eps(i,k) = 0._r8
         eu(i,k) = 0._r8
         du(i,k) = 0._r8
         ql(i,k) = 0._r8
         cu(i,k) = 0._r8
         evp(i,k) = 0._r8
         cmeg(i,k) = 0._r8
         qds(i,k) = q(i,k)
         md(i,k) = 0._r8
         ed(i,k) = 0._r8
         sd(i,k) = s(i,k)
         qd(i,k) = q(i,k)
         mc(i,k) = 0._r8
         qu(i,k) = q(i,k)
         su(i,k) = s(i,k)
         est(i) = c1*exp((c2* (t(i,k)-tfreez))/((t(i,k)-tfreez)+c3))
         if ( p(i,k)-est(i) > 0._r8 ) then
            qst(i,k) = eps1*est(i)/ (p(i,k)-est(i))
         else
            qst(i,k) = 1.0_r8
         end if
         gamma(i,k) = qst(i,k)*(1._r8 + qst(i,k)/eps1)*eps1*rl/(rd*t(i,k)**2)*rl/cp
         hmn(i,k) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
         hsat(i,k) = cp*t(i,k) + grav*z(i,k) + rl*qst(i,k)
         hu(i,k) = hmn(i,k)
         hd(i,k) = hmn(i,k)
         rprd(i,k) = 0._r8
      end do
   end do
!

  eps0(:) = 1.e-4_r8 
  mu(:,:) = 0.0_r8
  md(:,:) = 0.0_r8
  eu(:,:) = 0.0_r8
  du(:,:) = 0.0_r8
  ed(:,:) = 0.0_r8
  dd(:,:) = 0.0_r8
  hu = hmn
  su = s
  qu = q
  hd = hmn
  sd = s
  qd = q
  cu = 0.0_r8
  cmeg = 0.0_r8
  rprd = 0.0_r8
  ql = 0.0_r8
  evp = 0.0_r8
  pflx = 0.0_r8

!
   do k=1,msg
      do i=1,il2g
         rprd(i,k) = 0._r8
      end do
   end do
!
! interpolate the layer values of qst, hsat and gamma to
! layer interfaces
!
   do i = 1,il2g
      hsthat(i,msg+1) = hsat(i,msg+1)
      qsthat(i,msg+1) = qst(i,msg+1)
      gamhat(i,msg+1) = gamma(i,msg+1)
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (abs(qst(i,k-1)-qst(i,k)) > 1.E-6_r8) then
            qsthat(i,k) = log(qst(i,k-1)/qst(i,k))*qst(i,k-1)*qst(i,k)/ (qst(i,k-1)-qst(i,k))
         else
            qsthat(i,k) = qst(i,k)
         end if
         hsthat(i,k) = cp*shat(i,k) + rl*qsthat(i,k)
         if (abs(gamma(i,k-1)-gamma(i,k)) > 1.E-6_r8) then
            gamhat(i,k) = log(gamma(i,k-1)/gamma(i,k))*gamma(i,k-1)*gamma(i,k)/ &
                                (gamma(i,k-1)-gamma(i,k))
         else
            gamhat(i,k) = gamma(i,k)
         end if
      end do
   end do
!


!   call cal_buoysort(flagbspdf, bs_cridis, z(i,k), p(i,k)*100._r8, rho(i,k), &
!                     bs_thetal_e, q(i,k), bs_thetal_up, qu(i,k+1)+ql(i,k+1),  &
!                     !q_up(i,k+1)+qliq_up(i,k+1)+qice_up(i,k+1) ), &
!                     bs_wue, bs_xc(i,k), ent_turb(i,k), det_turb(i,k) )

 
!updraft
! jt < jt2 < jb1 < jb <jb2

   do i = 1,il2g
      mu(i,:) = 0.0

      mu(i,jb1(i))=1._r8  
      ! use maximum value , this can be modified

      hu(i,jb1(i)) = hmn(i,jb(i)) + cp*(tiedke_add + tpert(i)) + rl*qpert(i)
      su(i,jb1(i)) = s(i,jb(i))   + tiedke_add +tpert(i)
      qu(i,jb1(i)) = q(i,jb(i))   + qpert(i)

!      hu(i,jb1(i)) = hmn(i,jb1(i)) + cp*(tiedke_add + tpert(i)) + rl*qpert(i)
!      su(i,jb1(i)) = s(i,jb1(i))   + tiedke_add +tpert(i)
!      qu(i,jb1(i)) = q(i,jb1(i))   + qpert(i)

!      hu(i,jb1(i)) = (hmn(i,jb(i)) + hmn(i,jb1(i)) + hmn(i,jb2(i)))/3.0 
!      su(i,jb1(i)) = (s(i,jb(i)) + s(i,jb1(i)) + s(i,jb2(i)))/3.0 
!      qu(i,jb1(i)) = (q(i,jb(i)) + q(i,jb1(i)) + q(i,jb2(i)))/3.0 

      ! below jb1
      do k = jb1(i)+1,jb2(i) 
         dum0 = (zf(i,k) - zf(i,jb2(i)+1) )/(zf(i,jb1(i)) - zf(i,jb2(i)+1) ) 
         hu(i,k) = hu(i,jb1(i))  
         su(i,k) = su(i,jb1(i))
         mu(i,k) = mu(i,jb1(i))*dum0
      enddo
      do k = jb1(i),jb2(i)-1
         eu(i,k) = (mu(i,k) - mu(i,k+1) )/dz(i,k)
         du(i,k) = 0.0_r8
      enddo
 !        eu(i,jb2(i)) = mu(i,k)/dz(i,jb2(i))
         eu(i,jb2(i)) = mu(i,jb2(i))/dz(i,jb2(i))  ! xwang
         du(i,jb2(i)) = 0.0_r8

      ! above jb1
      do k=jb1(i)-1,jt2(i),-1

        dum0 = dz(i,k)
        dum1 = 1._r8/dum0 - (ent_tot(i,k) - det_tot(i,k))/2._r8  
        dum2 = 1._r8/dum0 + (ent_tot(i,k) - det_tot(i,k))/2._r8  
        if(dum1<1.e-6 .or. dum2<1.e-6_r8)then
          mu(i,k) = mu(i,k+1)*exp((ent_tot(i,k) - det_tot(i,k))*dum0  )
          eu(i,k) = mu(i,k+1)*(exp(ent_tot(i,k)*dum0) -1._r8)
          du(i,k) = mu(i,k+1) + eu(i,k) - mu(i,k)
        else
          mu(i,k) = mu(i,k+1)*dum2/dum1
          dumb = (mu(i,k)+mu(i,k+1))/2.0
          eu(i,k) = ent_tot(i,k) * dumb
          du(i,k) = det_tot(i,k) * dumb
        endif
!!
       if(mu(i,k) .le. 1.0e-6_r8)then
         jt(i) = k
         jt2(i) = jt(i)+1
         mu(i,k) = 0.0
         du(i,k) = mu(i,k+1)
         eu(i,k) = 0.0
         exit
       endif
!!!!!!!!!!!!!!!!!!!!!!
       enddo !k

! if(abs(lat(i)*180./3.1416-0.0) .le. 0.75 .and. abs(lon(i)*180./3.1416-70.0) .le. 0.75)then
!  print*,'meme',i, jt(i),jt2(i),jb1(i),jb2(i)
!  print*,'memu',mu(i,:)
!  print*,'meeu',ent_tot(i,:)
!  print*,'medu',mu(i,:)
! endif

       ! for now the following replaces ent and det above jt2
       do k=jt2(i)-1,jt(i),-1
         dum0 = (zf(i,jt(i)) - zf(i,k) )/(zf(i,jt(i)) - zf(i,jt2(i)) ) 
         mu(i,k) = mu(i,jt2(i))*dum0
         du(i,k) = (mu(i,k+1) - mu(i,k))/dz(i,k) 
         eu(i,k) = 0.0_r8
       enddo
     enddo

     freeze_rate(:,:) = 0.0_r8

  DO ITER = 1, NFREEZE
  !=======================
   do i = 1,il2g
       do k=jb1(i)-1,jt(i),-1
        dum0 = dz(i,k)
        if(mu(i,k)>0.02_r8)then
         hu(i,k) = (mu(i,k+1)*hu(i,k+1) +  dum0 * &
         (eu(i,k)*hmn(i,k)-du(i,k)*hu(i,k+1) + latice*freeze_rate(i,k)) )/mu(i,k)
        else
         hu(i,k) = hu(i,k+1) + ent_tot(i,k)*(hmn(i,k) - hu(i,k+1) )*dum0  
        endif

      enddo !k

       mu(i,jt(i)) = 0.0   !
       du(i,jt(i)) = mu(i,jt(i)+1)/dz(i,jt(i))

   enddo

   do i = 1,il2g
      done(i) = .false.
   end do
   kount = 0
   do k = pver,msg + 2,-1
      do i = 1,il2g
         if (( .not. done(i) .and. k > jt(i) .and. k < jb1(i)) ) then   
          if(mu(i,k) > 0.01_r8)then
            su(i,k) = mu(i,k+1)/mu(i,k)*su(i,k+1) + &
                      dz(i,k)/mu(i,k)* (eu(i,k)*s(i,k)-du(i,k)*su(i,k+1) &
                      + latice*freeze_rate(i,k))
            qu(i,k) = mu(i,k+1)/mu(i,k)*qu(i,k+1) + dz(i,k)/mu(i,k)* (eu(i,k)*q(i,k)- &
                            !du(i,k)*qst(i,k))
                            du(i,k)*qu(i,k+1))
          else
            su(i,k) = su(i,k+1) + ent_tot(i,k)*(s(i,k) - su(i,k+1))
            qu(i,k) = qu(i,k+1) + ent_tot(i,k)*(q(i,k) - qu(i,k+1))
          endif

            tu = su(i,k) - grav/cp*zf(i,k)
            estu = c1*exp((c2* (tu-tfreez))/ ((tu-tfreez)+c3))
            qstu = eps1*estu/ ((p(i,k)+p(i,k-1))/2._r8-estu)
            if (qu(i,k) >= qstu) then
               jlcl(i) = k
               kount = kount + 1
               done(i) = .true.
            end if

         end if
      end do
      if (kount >= il2g) goto 691
   end do
691 continue
   do k = msg + 2,pver
      do i = 1,il2g
         !if (k == jb1(i) ) then
         !   qu(i,k) = q(i,jb(i)) + qpert(i)
         !   su(i,k) = (hu(i,k)-rl*qu(i,k))/cp
         !end if
         if ((k > jt(i) .and. k <= jlcl(i)) .and. eps0(i) > 0._r8) then
            su(i,k) = shat(i,k) + (hu(i,k)-hsthat(i,k))/(cp* (1._r8+gamhat(i,k)))
            qu(i,k) = qsthat(i,k) + gamhat(i,k)*(hu(i,k)-hsthat(i,k))/ &
                     (rl* (1._r8+gamhat(i,k)))
         end if
      end do
   end do

! compute condensation in updraft
   do k = pver,msg + 2,-1
      do i = 1,il2g
         if (k >= jt(i) .and. k < jb2(i) ) then
            cu(i,k) = ((mu(i,k)*su(i,k)-mu(i,k+1)*su(i,k+1))/ &
                      dz(i,k)- (eu(i,k)*s(i,k) -du(i,k)*su(i,k+1)))/(rl/cp)
            if (k == jt(i)) cu(i,k) = 0._r8
            cu(i,k) = max(0._r8,cu(i,k))

            freeze_rate(i,k) = cu(i,k) * ramp_t1(i,k)
         end if
      end do
   end do

!k=1
!if(k<0)then
!   write(*,*)
!   write(*,*)'freeze_iter',iter
!   write(*,*)'cu        ',cu*8.64e7
!   write(*,*)'t         ',t
!   write(*,*)'ramp_t    ',ramp_t1
!   write(*,*)'freeze_rate',freeze_rate*8.64e7
!   write(*,*)'hu - h    ',(hu-hmn)/1004.
!   write(*,*)'su - s    ',su-s
!   write(*,*)'qu - q    ',(qu-q)*1000.
!endif

  ENDDO ! freeze iteration
  !=======================

 
!write(*,*)'cu',cu

! compute condensed water, rain-snow  production rate
! accumulate total precipitation (condensation - detrainment of liquid)
! Note ql1 = ql(k) + rprd(k)*dz(k)/mu(k)
! The differencing is somewhat strange (e.g. du(i,k)*ql(i,k+1)) but is
! consistently applied.
!    mu, ql are interface quantities
!    cu, du, eu, rprd are midpoint quantites

   do i = 1,il2g

         dum1 = 1.0_r8
         dum2 = 1.0_r8

    do k = jb(i),jt(i),-1

         if(k .le.jlcl(i) )then
          dum1 = dum1 + eu(i,k)*dz(i,k) 
          dum2 = 1._r8 - (eu(i,k)*dz(i,k))/dum1
         endif

!         dum0 = max(tanh((z(i,k) - zf(i,jlcl(i)) - rain_z0(i))/rain_zp(i)) &
!                  *(1.0-zmh_rain_cr*ramp_t2(i,k)), 0.0_r8)
    
         dum0 = max(tanh((z(i,k) - zf(i,jlcl(i)) - rain_z0(i))/rain_zp(i)), 0.0_r8)
         dum0 = min(1.0-zmh_rain_cr*ramp_t2(i,k) , sqrt(dum0))

         !!rprd(i,k) = dum0 * dum2 * cu(i,k) 
          rprd(i,k) = 0.0

         if (k >= jt(i) .and. k < jb2(i) ) then
            if (mu(i,k) > 0._r8) then
               ql1 = 1._r8/mu(i,k)* (mu(i,k+1)*ql(i,k+1)- &
                     dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*cu(i,k))

               ql1 = max(ql1,0.0)
              !!! ql(i,k) = ql1/ (1._r8+dz(i,k)*c0mask(i))
              !!! rprd(i,k) = c0mask(i)*mu(i,k)*ql(i,k)

               dum0 = c0mask(i)  !*rho(i,k)
               ql(i,k) = ql1/ (1._r8+dz(i,k)*dum0)
               rprd(i,k) = dum0*mu(i,k)*ql(i,k)

               qf(i,k) = 0.0  ! freezing not considered here
                 
           !!    ql1 = 1._r8/mu(i,k)* (mu(i,k+1)*ql(i,k+1)- &
           !!          dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)*(1._r8 - dum0)*cu(i,k))
           !!    ql(i,k) = ql1

           !!    ql2 = 1._r8/mu(i,k)* (mu(i,k+1)*ql(i,k+1)- &
           !!          dz(i,k)*du(i,k)*ql(i,k+1)+dz(i,k)* dum0*cu(i,k))
           !!    qf(i,k) = ql2

            else
               ql(i,k) = 0._r8
               qf(i,k) = 0._r8
            end if
            ql(i,k) = max( ql(i,k),0._r8)
            qf(i,k) = max( qf(i,k),0._r8)

            !rprd(i,k) = c0mask(i)*mu(i,k)*ql(i,k)
         end if
      end do
   end do
!

! =====================
! downdraft

  !dzb(:) = 0.8_r8   !thickness profile assumed, this is the entraining layer
  !dzt(:) = 0.2_r8   !detraining layer

! 190602
  !dzb(:) = zmh_dzb
  !dzt(:) = zmh_dzt

  dzb(:) = 0.5_r8   !thickness profile assumed, this is the entraining layer
  dzt(:) = 0.5_r8   !detraining layer


  jd3(:) = pver

  ent_turb_dn(:,:) = 0.0
  det_turb_dn(:,:) = 0.0

  ent_tot_dn(:,:) = 0.0
  det_tot_dn(:,:) = 0.0

  call ent_dynamics(il2g ,zf ,dz ,jd ,jd1, jd, jd3, jd2, dzb, dzt,&
      ent_turb_dn,det_turb_dn,eps00, ent_tot_dn,det_tot_dn, -1)

 ! jd <= jd1 < jd2 <= jd3 

   do i = 1,il2g
!
      
    if(jt2(i) > jd1(i)) cycle ! low interface of updraft top lower than lower interface of down

      alfa(i) = zmh_alfa
     
      md(i,jd1(i)) = -alfa(i)

    !above jd1
      do k = jd1(i)-1,jd(i),-1
       hd(i,k) = hmn(i,k)
       sd(i,k) = s(i,k)
       qd(i,k) = q(i,k)

       dum0 = (zf(i,jd(i)) - zf(i,k) )/(zf(i,jd(i)) - zf(i,jd1(i)) ) 
       md(i,k) = md(i,jd1(i))*dum0 
       ed(i,k) = (md(i,k)-md(i,k+1))/dz(i,k)
       dd(i,k) = 0.0_r8
      enddo

    !between jd1 and jd2
    do k = jd1(i),jd2(i)-1

        dum0 = dz(i,k)
        dum1 = 1._r8/dum0 - (ent_tot_dn(i,k) - det_tot_dn(i,k))/2._r8  
        dum2 = 1._r8/dum0 + (ent_tot_dn(i,k) - det_tot_dn(i,k))/2._r8  

        if(dum1<1.e-6 .or. dum2<1.e-6_r8)then
          md(i,k+1) = md(i,k)*exp((ent_tot_dn(i,k) - det_tot_dn(i,k))*dum0  )
          ed(i,k) = - md(i,k)*(exp(ent_tot_dn(i,k)*dum0) -1._r8)
          dd(i,k) = - md(i,k) + ed(i,k) + md(i,k+1)
        else
          md(i,k+1) = md(i,k)*dum2/dum1
          dumb = (md(i,k)+md(i,k+1))/2.0
          ed(i,k) = - ent_tot_dn(i,k) * dumb
          dd(i,k) = - det_tot_dn(i,k) * dumb
        endif
     enddo

     !from jd2 to jd3
    do k = jd2(i),jd3(i)-1
          dum0 = (zf(i,k+1) - zf(i,jd3(i)+1) )/(zf(i,jd2(i)) - zf(i,jd3(i)+1) ) 
          md(i,k+1) = dum0*md(i,jd2(i))
          dd(i,k)   = (md(i,k+1) - md(i,k))/dz(i,k)
          ed(i,k)   = 0.0_r8
    enddo
          dd(i,jd3(i)) = - md(i,jd3(i))/dz(i,jd3(i))
          ed(i,jd3(i)) = 0.0 

    do k = jd1(i),jd3(i)-1    
        dum0 = dz(i,k)

        if(md(i,k+1) < -0.01_r8)then 
          hd(i,k+1) = (md(i,k)*hd(i,k) -  dum0*(ed(i,k)*hmn(i,k) - &
                    dd(i,k)*hd(i,k)) )/md(i,k+1)
        else
          hd(i,k+1) = hd(i,k) + ent_tot_dn(i,k)*(hmn(i,k) - hd(i,k)) *dum0
        endif

    enddo !k

   enddo !i

   !!mc = mu + md
   

   do k = msg + 1,pver
      do i = 1,il2g
    if(jt2(i) > jd1(i)) cycle ! low interface of updraft top lower than lower interface of down
         if ((jt2(i)<=jd1(i)) .and. (k >= jd(i) .and. k <= jd3(i)) .and.   md(i,jb(i))< -0.01_r8 ) then
            ratmjb(i) = min(abs(mu(i,jb1(i))/min(md(i,jd2(i)),-0.01)),1._r8)
            md(i,k) = md(i,k)*ratmjb(i)
            ed(i,k) = ed(i,k)*ratmjb(i)
            dd(i,k) = dd(i,k)*ratmjb(i)
         end if
      end do
   end do

!properties of downdraft

   cmeg(:,:) = 0.0_r8
   pflx(:,:) = 0.0_r8
   evp(:,:)  = 0.0_r8
   qd(:,:) = q(:,:)
   sd(:,:) = s(:,:)
   

   do i=1,il2g
!     do k=jt(i),jd(i)+1
     do k=jt(i),jd(i)  ! xwang
         pflx(i,k) = pflx(i,k-1) + rprd(i,k-1)*dz(i,k-1)
     enddo
   enddo

   do i = 1,il2g
      qd(i,jd(i)) = q(i,jd(i)-1)
      sd(i,jd(i)) = (hd(i,jd(i)) - rl*qd(i,jd(i)))/cp
      
    do k = jd(i),jd3(i)-1

      dum1    = zf(i,jd(i))-zf(i,k+1)

      sd(i,k+1) = (s(i,k)+s(i,k+1))*0.5 - 1.5_r8 * zmh_ramp(dum1, 0._r8, 1.e3_r8)  

      qd(i,k+1) = (hd(i,k+1) - cp*sd(i,k+1))/rl

      evp(i,k) = (md(i,k)*qd(i,k) - md(i,k+1)*qd(i,k+1))/dz(i,k) - (ed(i,k)*q(i,k) - &
                    dd(i,k)*qd(i,k))

      evp(i,k) = max(evp(i,k),0._r8) 
      evp(i,k) = min(evp(i,k), pflx(i,k)/dz(i,k))

      rprd(i,k) = rprd(i,k)-evp(i,k)
      cmeg(i,k) = cu(i,k) - evp(i,k)
      pflx(i,k+1) = pflx(i,k) + rprd(i,k)*dz(i,k)

!re-do sd and qd with new evp

      if(md(i,k+1)< -0.01)then
       sd(i,k+1) = ( md(i,k)*sd(i,k) - (-rl/cp*evp(i,k)+ed(i,k)*s(i,k) - &
                    dd(i,k)*sd(i,k))*dz(i,k) )/md(i,k+1)
       qd(i,k+1) = ( md(i,k)*qd(i,k) - (evp(i,k)+ed(i,k)*q(i,k) - &
                    dd(i,k)*qd(i,k))*dz(i,k) )/md(i,k+1)
      else
       dum1 = min(md(i,k+1),-0.01)   
       sd(i,k+1) = sd(i,k) + (ent_tot_dn(i,k)*(s(i,k)-sd(i,k)) &
           + rl/cp*evp(i,k)/dum1 )*dz(i,k)
       qd(i,k+1) = qd(i,k) + (ent_tot_dn(i,k)*(q(i,k)-qd(i,k)) &
           - evp(i,k)/dum1       )*dz(i,k)
      endif

    enddo !k

   enddo !i
!----------------------------------------

   do k = msg + 1,pver
      do i = 1,il2g
         mc(i,k) = mu(i,k) + md(i,k)
!!
!!         du(i,k) = du(i,k) + dd(i,k)
      end do
   end do

 do i = 1,il2g
 capeout(i)  = 0.0_r8
 cinout(i)   = 0.0_r8
  do k = jt(i),jb(i)-1
   tv(i,k)    = t(i,k)* (1._r8+1.608_r8*q(i,k))/ (1._r8+q(i,k))   
   dum1       = (su(i,k)*(z(i,k)-zf(i,k+1))+su(i,k+1)*(zf(i,k)-z(i,k)) )/dz(i,k) 
   dum2       = (qu(i,k)*(z(i,k)-zf(i,k+1))+qu(i,k+1)*(zf(i,k)-z(i,k)) )/dz(i,k) 
   tv_up(i,k) = (dum1 - grav/cp*z(i,k))* (1._r8+1.608_r8*dum2)/ (1._r8+dum2)   

   rho(i,k)   = p(i,k)/rd/tv(i,k)*100.

   dum0 = ql(i,k) + qf(i,k)
!   buoy_up(i,k) = (tv_up(i,k) - tv(i,k))/tv(i,k)*grav - dum0
!   buoy_up(i,k) = (tv_up(i,k) - tv(i,k))/tv(i,k)*grav - dum0*(dum1 - grav/cp*z(i,k))/tv(i,k)*grav    ! xwang
   buoy_up(i,k) = (tv_up(i,k) - tv(i,k))/tv(i,k)*grav - dum0*grav  !new

   if(buoy_up(i,k) > 0. )then 
    capeout(i) = capeout(i) + buoy_up(i,k)*rho(i,k)*dz(i,k)
   else
    if(k > jd(i))then
     cinout(i)  = cinout(i) - buoy_up(i,k)*rho(i,k)*dz(i,k)
    endif
   endif

 enddo
enddo
#define SCM_OUT
#undef SCM_OUT
#ifdef SCM_OUT

if(k<0)then
 do i=1,il2g
 if(abs(lat(i)*180./3.1416-0.0) .le. 0.75 .and. abs(lon(i)*180./3.1416-70.0) .le. 0.75)then
  print*, ''
  print*,'i lon lat ',i, lon(i), lat(i)
  write(*,*)'mu      ',p(i,:)
  write(*,*)'ent_tot ',ent_tot(i,:)
  write(*,*)'det_tot ',det_tot(i,:)
  write(*,*)'mu      ',mu(i,:)
  write(*,*)'eu      ',eu(i,:)
  write(*,*)'du      ',du(i,:)
  write(*,*)'md      ',md(i,:)
  write(*,*)'ent_tot_dn',ent_tot_dn(i,:)
  write(*,*)'det_tot_dn',det_tot_dn(i,:)
  write(*,*)'mc      ',mu(i,:) + md(i,:)
  print*, ''
 endif
 enddo
endif

! for SCM diagnostic printout
if(k<0)then
  hmnf = hmn
  qf = q
  sf = s
 do i = 1, il2g
  do k = jt(i),jb(i)
  hmnf(:,k) = (hmn(:,k) + hmn(:,k-1))*0.5
  qf(:,k)   = (q(:,k)*(z(i,k-1)-zf(i,k)) + q(:,k-1)*(zf(i,k)-z(i,k)) )/(z(i,k-1)-z(i,k))
  sf(:,k)   = (s(:,k)*(z(i,k-1)-zf(i,k)) + s(:,k-1)*(zf(i,k)-z(i,k)) )/(z(i,k-1)-z(i,k))
  enddo
 enddo

  write(*,*)
  write(*,*)'new_xxxxxx'
  write(*,*)'jt,jd,jb    ',jt,jd,jb
  write(*,*)'jb1,jb2,jd2 ',jb1,jb2,jd2
  write(*,*)'eps00',eps00
  write(*,*)'capeout=',capeout
  write(*,*)'cinout=',cinout
  write(*,*)'ent_tot ',ent_tot
  write(*,*)'det_tot ',det_tot
  write(*,*)'mu',mu
  write(*,*)'eu',eu
  write(*,*)'du',du
  write(*,*)'hu    ',hu/1004.
  write(*,*)'hu-hmn',(hu-hmnf)/1004.
  write(*,*)'su    ',su
  write(*,*)'su-s  ',(su-sf)
  write(*,*)'qu    ',qu*1000.
  write(*,*)'qu-q  ',(qu-qf)*1000.
  write(*,*)
  write(*,*)'ent_turb',ent_turb
  write(*,*)'det_turb', det_turb
  write(*,*)'md',md
  write(*,*)'ed',ed
  write(*,*)'dd',dd
  write(*,*)'hd    ',hd/1004.
  write(*,*)'hd-hmn',(hd-hmnf)/1004.
  write(*,*)'sd    ',sd
  write(*,*)'sd-s  ',(sd-sf)
  write(*,*)'qd    ',qd*1000.
  write(*,*)'qd-q  ',(qd-qf)*1000.

  write(*,*)
  write(*,*)'cu    ',cu*86400.*1000.
  write(*,*)'ql    ',ql*1000.
  write(*,*)'evp   ',evp*86400.*1000.
  write(*,*)'pflx  ',pflx*86400.*1000.
  write(*,*)'rprd  ',rprd*86400.*1000.
endif
#ENDIF

   return

end subroutine zmh_plume_prp

! ====================================================
subroutine q1q2_zmh(lchnk   , &
                    dqdt    ,dsdt    ,q       ,qs      ,qu      , &
                    su      ,du      ,qhat0    ,shat0    ,dp      , &
                    mu      ,md      ,sd      ,qd      ,ql      , &
                    dsubcld ,jt      ,mx      ,il1g    ,il2g    , &
                    cp      ,rl      ,msg     ,          &
                    dl      ,evp     ,cu,                &
                    dtime, nt_stepping, k950, flag)


   implicit none

!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: phil rasch dec 19 1995 as q1q2_pjr
!
! Modified by Zhang Minghua on July 15 2018
! 
!-----------------------------------------------------------------------


   real(r8), intent(in) :: cp

   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: il1g
   integer, intent(in) :: il2g
   integer, intent(in) :: msg

   logical,  intent(in) :: flag(pcols)

   real(r8), intent(in) :: q(pcols,pver)
   real(r8), intent(in) :: qs(pcols,pver)
   real(r8), intent(in) :: qu(pcols,pver)
   real(r8), intent(in) :: su(pcols,pver)
   real(r8), intent(in) :: du(pcols,pver)
   real(r8), intent(in) :: qhat0(pcols,pver)
   real(r8), intent(in) :: shat0(pcols,pver)
   real(r8), intent(in) :: dp(pcols,pver)
   real(r8), intent(in) :: mu(pcols,pver)
   real(r8), intent(in) :: md(pcols,pver)
   real(r8), intent(in) :: sd(pcols,pver)
   real(r8), intent(in) :: qd(pcols,pver)
   real(r8), intent(in) :: ql(pcols,pver)
   real(r8), intent(in) :: evp(pcols,pver)
   real(r8), intent(in) :: cu(pcols,pver)
   real(r8), intent(in) :: dsubcld(pcols)
! zmh 
   integer k950(pcols)

   real(r8), intent(in) :: dtime
   integer, intent(in) :: nt_stepping(pcols)

   real(r8),intent(out) :: dqdt(pcols,pver),dsdt(pcols,pver)
   real(r8),intent(out) :: dl(pcols,pver)
   integer kbm
   integer ktm
   integer jt(pcols)
   integer mx(pcols)
!
! work fields:
!
   real(r8) :: qhat(pcols,pver)
   real(r8) :: shat(pcols,pver)
   integer i,m,nt, n
   integer k

   real(r8) emc
   real(r8) rl
   real(r8) dtime2

!-------------------------------------------------------------------
   do k = msg + 1,pver
      do i = il1g,il2g
         dsdt(i,k) = 0._r8
         dqdt(i,k) = 0._r8
         dl(i,k) = 0._r8
      end do
   end do
!
! find the highest level top and bottom levels of convection
!
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

! --------------------------------


 do i = il1g,il2g

 if( .not. flag(i)) cycle

  nt = nt_stepping(i)
  dtime2 = dtime/dble(nt)

  shat(i,:) = shat0(i,:)
  qhat(i,:) = qhat0(i,:)


!write(*,*)' in q1q2'
  do n = 1,nt  
! --------------------------

!   write(*,*)'n=',n
!   write(*,*)'nt_stepping2',nt_stepping(i)
!   write(*,*)'shat(i,:)',shat(i,:)
!   write(*,*)'qhat(i,:)',qhat(i,:)
!   write(*,*)'dsdt(i,:)',dsdt(i,:)
!   write(*,*)'dqdt(i,:)',dqdt(i,:)

   do k = jt(i),pver-1
         emc = -cu (i,k)               &         ! condensation in updraft
               +evp(i,k)                         ! evaporating rain in downdraft

         dsdt(i,k) = -rl/cp*emc &
                     + (+mu(i,k+1)* (su(i,k+1)-shat(i,k+1)) &
                        -mu(i,k)*   (su(i,k)-shat(i,k)) &
                        +md(i,k+1)* (sd(i,k+1)-shat(i,k+1)) &
                        -md(i,k)*   (sd(i,k)-shat(i,k)) &
                       )/dp(i,k)

         dqdt(i,k) = emc + &
                    (+mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)) &
                     -mu(i,k)*   (qu(i,k)-qhat(i,k)) &
                     +md(i,k+1)* (qd(i,k+1)-qhat(i,k+1)) &
                     -md(i,k)*   (qd(i,k)-qhat(i,k)) &
                    )/dp(i,k)

         dl(i,k) = du(i,k)*ql(i,k+1)

      end do !k

   shat(i,:) = shat(i,:) + dsdt(i,:)*dtime2
   qhat(i,:) = qhat(i,:) + dqdt(i,:)*dtime2

   end do  !m for nt
!  -------------------

   dsdt(i,:) = (shat(i,:) - shat0(i,:) )/dtime
   dqdt(i,:) = (qhat(i,:) - qhat0(i,:) )/dtime

!
!DIR$ NOINTERCHANGE!

   k = k950(i)
   dsdt(i,k) = (1._r8/dsubcld(i))* &
!                        (-mu(i,k)* (su(i,k)-shat(i,k)) &
!                         -md(i,k)* (sd(i,k)-shat(i,k)) &
                        (-mu(i,k)* (su(i,k)-shat0(i,k)) &
                         -md(i,k)* (sd(i,k)-shat0(i,k)) &
                        )
   dqdt(i,k) = (1._r8/dsubcld(i))* &
!                        (-mu(i,k)*(qu(i,k)-qhat(i,k)) &
!                         -md(i,k)*(qd(i,k)-qhat(i,k)) &
                        (-mu(i,k)*(qu(i,k)-qhat0(i,k)) &
                         -md(i,k)*(qd(i,k)-qhat0(i,k)) &
                        )

    do k = k950(i)+1,pver
            dsdt(i,k) = dsdt(i,k-1)
            dqdt(i,k) = dqdt(i,k-1)
     end do


   end do  !i
!
   return
end subroutine q1q2_zmh

subroutine zyx1_buoyan_dilute(lchnk   ,ncol    , limcnv, &
                  q       ,t       ,p       ,z       ,pf      , &
                  tp      ,qstp    ,tl      ,rl      ,cape    , &
                  pblt    ,lcl     ,lel     ,mx      ,mx1     ,mx2, &
                  rd      ,grav    ,cp      ,msg     , &
                  tpert   , qpert, ent_tot, buoy, cin)         !zmh
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculates CAPE the lifting condensation level and the convective top
! where buoyancy is first -ve.
! 
! Method: Calculates the parcel temperature based on a simple constant
! entraining plume model. CAPE is integrated from buoyancy.
! 09/09/04 - Simplest approach using an assumed entrainment rate for 
!            testing (dmpdp). 
! 08/04/05 - Swap to convert dmpdz to dmpdp  
!
! SCAM Logical Switches - DILUTE:RBN - Now Disabled 
! ---------------------
! switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
! switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
! switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
! 
! References:
! Raymond and Blythe (1992) JAS 
! 
! Author:
! Richard Neale - September 2004
! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: limcnv

   real(r8), intent(in) :: q(pcols,pver)        ! spec. humidity
   real(r8), intent(in) :: t(pcols,pver)        ! temperature
   real(r8), intent(in) :: p(pcols,pver)        ! pressure
   real(r8), intent(in) :: z(pcols,pver)        ! height
   real(r8), intent(in) :: pf(pcols,pver+1)     ! pressure at interfaces
   real(r8), intent(in) :: pblt(pcols)          ! index of pbl depth
   real(r8), intent(in) :: tpert(pcols)         ! perturbation temperature by pbl processes
!zmh   
   integer, intent(in) :: mx(pcols)         ! perturbation humidity by pbl processes
   integer, intent(in) :: mx1(pcols) ! not used yet 
   integer, intent(in) :: mx2(pcols)   
   real(r8), intent(in) :: qpert(pcols)         ! perturbation humidity by pbl processes
   real(r8), intent(in) :: ent_tot(pcols,pver) 
!
! output arguments
!
   real(r8), intent(out) :: tp(pcols,pver)       ! parcel temperature
   real(r8), intent(out) :: qstp(pcols,pver)     ! saturation mixing ratio of parcel (only above lcl, just q below).
   real(r8), intent(out) :: tl(pcols)            ! parcel temperature at lcl
   real(r8), intent(out) :: cape(pcols)          ! convective aval. pot. energy.
!zmh
   real(r8), intent(out) :: buoy(pcols,pver) 
   real(r8), intent(out) :: cin(pcols)        


   integer lcl(pcols)        !
   integer lel(pcols)        !
   integer lon(pcols)        ! level of onset of deep convection
!   integer mx(pcols)         ! level of max moist static energy
!--------------------------Local Variables------------------------------
!
   real(r8) capeten(pcols,5)     ! provisional value of cape
   real(r8) tv(pcols,pver)       !
   real(r8) tpv(pcols,pver)      !
   !real(r8) buoy(pcols,pver)

   real(r8) a1(pcols)
   real(r8) a2(pcols)
   real(r8) estp(pcols)
   real(r8) pl(pcols)
   real(r8) plexp(pcols)
   real(r8) hmax(pcols)
   real(r8) hmn(pcols)
   real(r8) y(pcols)

   logical plge600(pcols)
   integer knt(pcols)
   integer lelten(pcols,5)

   real(r8) cp
   real(r8) e
   real(r8) grav

   integer i
   integer k
   integer k100(pcols)  !zmh
   integer msg
   integer n

   real(r8) rd
   real(r8) rl
#ifdef PERGRO
   real(r8) rhd
#endif
!
!-----------------------------------------------------------------------
!
   do n = 1,5
      do i = 1,ncol
         lelten(i,n) = pver
         capeten(i,n) = 0._r8
      end do
   end do
!
   do i = 1,ncol
      lon(i) = mx(i)
      knt(i) = 0
      lel(i) = limcnv !zmh
      cape(i) = 0._r8
      hmax(i) = 0._r8
   end do

   tp(:ncol,:) = t(:ncol,:)
   qstp(:ncol,:) = q(:ncol,:)

!!! RBN - Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.
   tv(:ncol,:) = t(:ncol,:) *(1._r8+1.608_r8*q(:ncol,:))/ (1._r8+q(:ncol,:))
   tpv(:ncol,:) = tv(:ncol,:)
   buoy(:ncol,:) = 0._r8

!
!============================================================
   

   do i = 1,ncol ! Initialise LCL variables.
      lcl(i) = mx(i)
      tl(i)  = t(i,mx(i))
      pl(i)  = p(i,mx(i))
      k      = mx(i)
      hmax(i)= cp*t(i,k) + grav*z(i,k) + rl*q(i,k) 
   end do

!
! main buoyancy calculation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DILUTE PLUME CALCULATION USING ENTRAINING PLUME !!!
!!!   RBN 9/9/04   !!!

   call zyx1_plume_dilute(lchnk, ncol, msg, mx, p, t, q, tpert, qpert, ent_tot, &
        tp, tpv, qstp, pl, tl, lcl)


! If lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
   do i = 1,ncol
      plge600(i) = pl(i).ge.600._r8 ! Just change to always allow buoy calculation.
   end do

!
! Main buoyancy calculation.
!
   do k = pver,msg + 1,-1
      do i=1,ncol
         if (k <= mx(i) .and. plge600(i)) then   ! Define buoy from launch level to cloud top.
            tv(i,k) = t(i,k)* (1._r8+1.608_r8*q(i,k))/ (1._r8+q(i,k))
            buoy(i,k) = tpv(i,k) - tv(i,k) !+ tiedke_add  ! +0.5K or not?
         else
            qstp(i,k) = q(i,k)
            tp(i,k)   = t(i,k)            
            tpv(i,k)  = tv(i,k)
         endif
      end do
   end do

!-------------------------------------------------------------------------------

!
   do k = msg + 2,pver
      do i = 1,ncol
         if (k < lcl(i) .and. plge600(i)) then
            if (buoy(i,k+1) > 0. .and. buoy(i,k) <= 0._r8) then
               knt(i) = min(5,knt(i) + 1)
               lelten(i,knt(i)) = k
            end if
         end if
      end do
   end do
!
! calculate convective available potential energy (cape).
!
   do n = 1,5
      do k = msg + 1,pver
         do i = 1,ncol
            if (plge600(i) .and. k <= mx(i) .and. k > lelten(i,n)) then
               capeten(i,n) = capeten(i,n) + rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
         end do
      end do
   end do
!
! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
!
   do n = 1,5
      do i = 1,ncol
         if (capeten(i,n) > cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
         end if
      end do
   end do

! put lower bound on cape for diagnostic purposes.
!
   do i = 1,ncol
      cape(i) = max(cape(i), 0._r8)
   end do
!
!zmh
   cin (:) = 0.0
   do i = 1,ncol
     do k= pver-1,msg + 1,-1
       if(buoy(i,k) .le. 0._r8)then
          cin(i) = cin(i) - rd*min(buoy(i,k),0._r8)*log(pf(i,k+1)/pf(i,k))
       else if (k < (nint(pblt(i))-5) )then
         exit 
       endif 

     enddo  ! k

     if(lel(i)>mx(i)-3)then  !too shallow
       cape(i) = 0._r8
     endif

     do k= pver-1,msg + 1,-1
      buoy(i,k) =  buoy(i,k)*gravit/t(i,k)  !buoyancy force
     enddo

   enddo
!zmh end
!
   return
end subroutine zyx1_buoyan_dilute

! ==========================================================

subroutine zyx1_plume_dilute(lchnk, ncol, msg, klaunch, p, t, q, &
         tpert, qpert,ent_tot, tp, tpv, qstp, pl, tl, lcl)

! Routine  to determine 
!   1. Tp   - Parcel temperature
!   2. qstp - Saturated mixing ratio at the parcel temperature.

!--------------------
implicit none
!--------------------

integer, intent(in) :: lchnk
integer, intent(in) :: ncol
integer, intent(in) :: msg

integer, intent(in), dimension(pcols) :: klaunch(pcols)

real(r8), intent(in), dimension(pcols,pver) :: p
real(r8), intent(in), dimension(pcols,pver) :: t
real(r8), intent(in), dimension(pcols,pver) :: q
real(r8), intent(in), dimension(pcols) :: tpert ! PBL temperature perturbation.
!zmh
real(r8), intent(in), dimension(pcols) :: qpert ! PBL temperature perturbation.
real(r8), intent(in) :: ent_tot(pcols,pver) 

real(r8), intent(inout), dimension(pcols,pver) :: tp    ! Parcel temp.
real(r8), intent(inout), dimension(pcols,pver) :: qstp  ! Parcel water vapour (sat value above lcl).
real(r8), intent(inout), dimension(pcols) :: tl         ! Actual temp of LCL.
real(r8), intent(inout), dimension(pcols) :: pl          ! Actual pressure of LCL. 

integer, intent(inout), dimension(pcols) :: lcl ! Lifting condesation level (first model level with saturation).

real(r8), intent(out), dimension(pcols,pver) :: tpv   ! Define tpv within this routine.

!--------------------

! Have to be careful as s is also dry static energy.


! If we are to retain the fact that CAM loops over grid-points in the internal
! loop then we need to dimension sp,atp,mp,xsh2o with ncol.


real(r8) tmix(pcols,pver)        ! Tempertaure of the entraining parcel.
real(r8) qtmix(pcols,pver)       ! Total water of the entraining parcel.
real(r8) qsmix(pcols,pver)       ! Saturated mixing ratio at the tmix.
real(r8) smix(pcols,pver)        ! Entropy of the entraining parcel.
real(r8) xsh2o(pcols,pver)       ! Precipitate lost from parcel.
real(r8) ds_xsh2o(pcols,pver)    ! Entropy change due to loss of condensate.
real(r8) ds_freeze(pcols,pver)   ! Entropy change sue to freezing of precip.

real(r8) mp(pcols)    ! Parcel mass flux.
real(r8) qtp(pcols)   ! Parcel total water.
real(r8) sp(pcols)    ! Parcel entropy.

real(r8) sp0(pcols)    ! Parcel launch entropy.
real(r8) qtp0(pcols)   ! Parcel launch total water.
real(r8) mp0(pcols)    ! Parcel launch relative mass flux.

real(r8) lwmax      ! Maximum condesate that can be held in cloud before rainout.
real(r8) dmpdp      ! Parcel fractional mass entrainment rate (/mb).
!real(r8) dmpdpc     ! In cloud parcel mass entrainment rate (/mb).
real(r8) dmpdz      ! Parcel fractional mass entrainment rate (/m)
real(r8) dpdz,dzdp  ! Hydrstatic relation and inverse of.
real(r8) senv       ! Environmental entropy at each grid point.
real(r8) qtenv      ! Environmental total water "   "   ".
real(r8) penv       ! Environmental total pressure "   "   ".
real(r8) tenv       ! Environmental total temperature "   "   ".
real(r8) new_s      ! Hold value for entropy after condensation/freezing adjustments.
real(r8) new_q      ! Hold value for total water after condensation/freezing adjustments.
real(r8) dp         ! Layer thickness (center to center)
real(r8) tfguess    ! First guess for entropy inversion - crucial for efficiency!
real(r8) tscool     ! Super cooled temperature offset (in degC) (eg -35).

real(r8) qxsk, qxskp1        ! LCL excess water (k, k+1)
real(r8) dsdp, dqtdp, dqxsdp ! LCL s, qt, p gradients (k, k+1)
real(r8) slcl,qtlcl,qslcl    ! LCL s, qt, qs values.

integer rcall       ! Number of ientropy call for errors recording
integer nit_lheat     ! Number of iterations for condensation/freezing loop.
integer i,k,ii   ! Loop counters.
real(r8) tmp   !zhh debug

!======================================================================
!    SUMMARY
!
!  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
!           and entrains at each level with a specified entrainment rate.
!
! 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.          
!
!======================================================================
!
! Set some values that may be changed frequently.
!

nit_lheat = 2 ! iterations for ds,dq changes from condensation freezing.

lwmax = 1.e-3_r8    ! Need to put formula in for this.
tscool = 0.0_r8   ! Temp at which water loading freezes in the cloud.

qtmix=0._r8
smix=0._r8

qtenv = 0._r8
senv = 0._r8
tenv = 0._r8
penv = 0._r8

qtp0 = 0._r8
sp0  = 0._r8
mp0 = 0._r8

qtp = 0._r8
sp = 0._r8
mp = 0._r8

new_q = 0._r8
new_s = 0._r8

! **** Begin loops ****


do k = pver, msg+1, -1
   do i=1,ncol 

! Initialize parcel values at launch level.

      if (k == klaunch(i)) then 
!zmh 1
         qtp0(i) = q(i,k) + qpert(i)   !
         sp0(i)  = entropy(t(i,k) +tpert(i),p(i,k),qtp0(i))  ! Parcel launch entropy.

!         qtp0(i) = q(i,k-1) + qpert(i)   !
!         sp0(i)  = entropy(t(i,k-1) +tpert(i),p(i,k),qtp0(i))  ! Parcel launch entropy.

         qtp(i) = qtp0(i)
         sp(i)  = sp0(i) 
         
         mp0(i)  = 1._r8       ! Parcel launch relative mass (i.e. 1 parcel stays 1 parcel for dmpdp=0, undilute). 
         smix(i,k)  = sp0(i)
         qtmix(i,k) = qtp0(i)

         tfguess = t(i,k)   + tpert(i)

         rcall = 1
         call ientropy (rcall,i,lchnk,smix(i,k),p(i,k),qtmix(i,k),tmix(i,k),qsmix(i,k),tfguess)
      end if

! Entraining levels
      
      if (k < klaunch(i)) then 

! Set environmental values for this level.                 
         
         dp = (p(i,k)-p(i,k+1)) ! In -ve mb as p decreasing with height - difference between center of layers.

         dmpdz =  - ent_tot(i,k)

         qtenv = 0.5_r8*(q(i,k)+q(i,k+1))         ! Total water of environment.
         tenv  = 0.5_r8*(t(i,k)+t(i,k+1)) 
         penv  = 0.5_r8*(p(i,k)+p(i,k+1))

!         qtenv = 0.5_r8*(q(i,k)+q(i,k))         ! Total water of environment.
!         tenv  = 0.5_r8*(t(i,k)+t(i,k)) 
!         penv  = 0.5_r8*(p(i,k)+p(i,k))

         senv  = entropy(tenv,penv,qtenv)  ! Entropy of environment.   

         dpdz = -(penv*grav)/(rgas*tenv) ! in mb/m since  p in mb.
         dzdp = 1._r8/dpdz                  ! in m/mb
         dmpdp = dmpdz*dzdp              ! /mb Fractional entrainment

         sp(i)  = (sp(i)  - dmpdp*dp*senv)/(1._r8 - dmpdp*dp)
         qtp(i) = (qtp(i) - dmpdp*dp*qtenv)/(1._r8 - dmpdp*dp) 
         smix(i,k) = sp(i)
         qtmix(i,k) = qtp(i)

!write(*,*)'k,mp,smix(i,k),senv'
!write(*,*)k,mp,smix(i,k),senv
!write(*,*)'k,mp,mp0,smix,qtmix,qtenv',k,mp0(i),mp(i),smix(i,k),qtmix(i,k),qtenv 

! Invert entropy from s and q to determine T and saturation-capped q of mixture.
! t(i,k) used as a first guess so that it converges faster.

!!zmh         tfguess = tmix(i,k+1)
         tfguess = tmix(i,k+1) + min(t(i,k)-t(i,k+1),-dzdp*dp*0.005)  !!!!!!!!!!

         rcall = 2
!-------------------------- zhh debug 2013.03.06 -------------------------
!!         if (lchnk==145 .or. lchnk==913 .or. lchnk==1041 .or. lchnk==1937) then
!!         if (lchnk==913) then
!!            print*, 'sp0(i) =', sp0(i), ' sp(i) =', sp(i)            
!!            print*, 'mp0(i) =', mp0(i), ' mp(i) =', mp(i)            
!!            print*, 'lchnk=', lchnk, ' i =',i, ' k =',k
!!            print*, 'sp(i) =', sp(i), ' dp =', dp
!!            print*, 'dzdp =', dzdp, ' senv =', senv
!!            print*, 'q(i,k) =', q(i,k), ' q(i,k+1) =', q(i,k+1)
!!            print*, 't(i,k) =', t(i,k), ' t(i,k+1) =', t(i,k+1)
!!            print*, 'p(i,k) =', p(i,k), ' p(i,k+1) =', p(i,k+1)
!!            print*, 'senv =', senv, ' tenv =', tenv
!!            print*, 'penv =', penv, ' qtenv =', qtenv
!!         end if
!         if ( abs(sp(i)) < 1E6 ) then
!            tmp = 0.
!         else
!            print*, 'lchnk=', lchnk, ' i =',i, ' k =',k
!            print*, 'sp(i) =', sp(i), ' dp =', dp
!            print*, 'dzdp =', dzdp, ' senv =', senv
!            print*, 't(i,k) =', t(i,k), ' t(i,k+1) =', t(i,k+1)
!         end if
!-------------------------- zhh debug 2013.03.06 -------------------------
         call ientropy(rcall,i,lchnk,smix(i,k),p(i,k),qtmix(i,k),tmix(i,k),qsmix(i,k),tfguess)   

!
! Determine if this is lcl of this column if qsmix <= qtmix.
! FIRST LEVEL where this happens on ascending.

         if (qsmix(i,k) <= qtmix(i,k) .and. qsmix(i,k+1) > qtmix(i,k+1)) then
            lcl(i) = k
            qxsk   = qtmix(i,k) - qsmix(i,k)
            qxskp1 = qtmix(i,k+1) - qsmix(i,k+1)
            dqxsdp = (qxsk - qxskp1)/dp
            pl(i)  = p(i,k+1) - qxskp1/dqxsdp    ! pressure level of actual lcl.
            dsdp   = (smix(i,k)  - smix(i,k+1))/dp
            dqtdp  = (qtmix(i,k) - qtmix(i,k+1))/dp
            slcl   = smix(i,k+1)  +  dsdp* (pl(i)-p(i,k+1))  
            qtlcl  = qtmix(i,k+1) +  dqtdp*(pl(i)-p(i,k+1))

!zmh            tfguess = tmix(i,k)
            tfguess = tmix(i,k) + (tmix(i,k+1)-tmix(i,k))&
                     *(pl(i)-p(i,k))/(p(i,k+1)-p(i,k))

            rcall = 3
            call ientropy (rcall,i,lchnk,slcl,pl(i),qtlcl,tl(i),qslcl,tfguess)

!            write(iulog,*)' '
!            write(iulog,*)' p',p(i,k+1),pl(i),p(i,lcl(i))
!            write(iulog,*)' t',tmix(i,k+1),tl(i),tmix(i,lcl(i))
!            write(iulog,*)' s',smix(i,k+1),slcl,smix(i,lcl(i))
!            write(iulog,*)'qt',qtmix(i,k+1),qtlcl,qtmix(i,lcl(i))
!            write(iulog,*)'qs',qsmix(i,k+1),qslcl,qsmix(i,lcl(i))

         endif
!         
      end if !  k < klaunch

 
   end do ! Levels loop
end do ! Columns loop

!!!!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Could stop now and test with this as it will provide some estimate of buoyancy
!! without the effects of freezing/condensation taken into account for tmix.

!! So we now have a profile of entropy and total water of the entraining parcel
!! Varying with height from the launch level klaunch parcel=environment. To the 
!! top allowed level for the existence of convection.

!! Now we have to adjust these values such that the water held in vaopor is < or 
!! = to qsmix. Therefore, we assume that the cloud holds a certain amount of
!! condensate (lwmax) and the rest is rained out (xsh2o). This, obviously 
!! provides latent heating to the mixed parcel and so this has to be added back 
!! to it. But does this also increase qsmix as well? Also freezing processes
 

xsh2o = 0._r8
ds_xsh2o = 0._r8
ds_freeze = 0._r8

!!!!!!!!!!!!!!!!!!!!!!!!!PRECIPITATION/FREEZING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Iterate solution twice for accuracy



do k = pver, msg+1, -1
   do i=1,ncol    
      
! Initialize variables at k=klaunch
      
      if (k == klaunch(i)) then

! Set parcel values at launch level assume no liquid water.            

         tp(i,k)    = tmix(i,k)
! zmh 3        
         !qstp(i,k)  = q(i,k) 
         qstp(i,k)  = q(i,k)  + qpert(i)
         tpv(i,k)   =  (tp(i,k) + tpert(i)) * (1._r8+1.608_r8*qstp(i,k)) / (1._r8+qstp(i,k))
         
      end if

      if (k < klaunch(i)) then
            
! Initiaite loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.

! Iterate nit_lheat times for s,qt changes.

         do ii=0,nit_lheat-1            

! Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).

            xsh2o(i,k) = max (0._r8, qtmix(i,k) - qsmix(i,k) - lwmax)

! Contribution to ds from precip loss of condensate (Accumulated change from smix).(-ve)                     
                     
            ds_xsh2o(i,k) = ds_xsh2o(i,k+1) - cpliq * log (tmix(i,k)/tfreez) * max(0._r8,(xsh2o(i,k)-xsh2o(i,k+1)))
!
! Entropy of freezing: latice times amount of water involved divided by T.
!
 
            if (tmix(i,k) <= tfreez+tscool .and. ds_freeze(i,k+1) == 0._r8) then ! One off freezing of condensate. 
               ds_freeze(i,k) = (latice/tmix(i,k)) * max(0._r8,qtmix(i,k)-qsmix(i,k)-xsh2o(i,k)) ! Gain of LH
            end if
            
            if (tmix(i,k) <= tfreez+tscool .and. ds_freeze(i,k+1) /= 0._r8) then ! Continual freezing of additional condensate.
               ds_freeze(i,k) = ds_freeze(i,k+1)+(latice/tmix(i,k)) * max(0._r8,(qsmix(i,k+1)-qsmix(i,k)))
            end if
            
! Adjust entropy and accordingly to sum of ds (be careful of signs).

            new_s = smix(i,k) + ds_xsh2o(i,k) + ds_freeze(i,k) 

! Adjust liquid water and accordingly to xsh2o.

            new_q = qtmix(i,k) - xsh2o(i,k)

! Invert entropy to get updated Tmix and qsmix of parcel.

            tfguess = tmix(i,k)
            rcall =4
            call ientropy (rcall,i,lchnk,new_s, p(i,k), new_q, tmix(i,k), qsmix(i,k), tfguess)
            
         end do  ! Iteration loop for freezing processes.

! tp  - Parcel temp is temp of mixture.
! tpv - Parcel v. temp should be density temp with new_q total water. 

         tp(i,k)    = tmix(i,k)

! tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)

         if (new_q > qsmix(i,k)) then  ! Super-saturated so condensate present - reduces buoyancy.
            qstp(i,k) = qsmix(i,k)
         else                          ! Just saturated/sub-saturated - no condensate virtual effects.
            qstp(i,k) = new_q
         end if

!zmh 4
         !tpv(i,k) = (tp(i,k)+tpert(i))* (1._r8+1.608_r8*qstp(i,k)) / (1._r8+ new_q) 
         tpv(i,k) = (tp(i,k)+0.*tpert(i))* (1._r8+1.608_r8*qstp(i,k)) / (1._r8+ new_q) 

      end if ! k < klaunch
      
   end do ! Loop for columns
   
end do  ! Loop for vertical levels.


return
end subroutine zyx1_plume_dilute

!-----------------------------------------------------------------------------------------
real(r8) function entropy(TK,p,qtot)
!-----------------------------------------------------------------------------------------
!
! TK(K),p(mb),qtot(kg/kg)
! from Raymond and Blyth 1992
!
     real(r8), intent(in) :: p,qtot,TK
     real(r8) :: qv,qsat,e,esat,L,eref,pref

pref = 1000.0_r8           ! mb
eref = 6.106_r8            ! sat p at tfreez (mb)

L = rl - (cpliq - cpwv)*(TK-tfreez)         ! T IN CENTIGRADE

! Replace call to satmixutils.

esat = c1*exp(c2*(TK-tfreez)/(c3+TK-tfreez))       ! esat(T) in mb
qsat=eps1*esat/(p-esat)                      ! Sat. mixing ratio (in kg/kg).

qv = min(qtot,qsat)                         ! Partition qtot into vapor part only.
e = qv*p / (eps1 +qv)

entropy = (cpres + qtot*cpliq)*log( TK/tfreez) - rgas*log( (p-e)/pref ) + &
        L*qv/TK - qv*rh2o*log(qv/qsat)
!
if (abs(entropy)<1E6) then
!
else
!!   print*, 'entropy = ', entropy
!!   print*, 'esat =', esat, ' qsat =', qsat
!!   print*, '(cpres + qtot*cpliq)*log( TK/tfreez) = ', (cpres + qtot*cpliq)*log( TK/tfreez)
!!   print*, 'rgas*log( (p-e)/pref ) = ', rgas*log( (p-e)/pref )
!!   print*, 'L*qv/TK = ', L*qv/TK
!!   print*, 'qv*rh2o*log(qv/qsat) = ', qv*rh2o*log(qv/qsat)
!!   print*, 'TK =', TK, ' qtot =', qtot
!!   print*, 'rh2o =', rh2o, ' qsat =', qsat
!!   print*, 'qv =', qv, ' qsat =', qsat
!!   print*, 'eps1 =', eps1, ' p =', p
!!   print*, 'esat =', esat
!!   print*, 'c1 =', c1, ' c2 =', c2
!!   print*, 'c3 =', c3, ' tfreez =', tfreez
end if
 
return
end FUNCTION entropy

!
!-----------------------------------------------------------------------------------------
   SUBROUTINE ientropy (rcall,icol,lchnk,s,p,qt,T,qsat,Tfg)
!-----------------------------------------------------------------------------------------
!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg). 
! Inverts entropy, pressure and total water qt 
! for T and saturated vapor mixing ratio
! 

     use phys_grid, only: get_rlon_p, get_rlat_p

     integer, intent(in) :: icol, lchnk, rcall
     real(r8), intent(in)  :: s, p, Tfg, qt
     real(r8), intent(out) :: qsat, T
     real(r8) :: qv,Ts,dTs,fs1,fs2,esat     
     real(r8) :: pref,eref,L,e
     real(r8) :: this_lat,this_lon
     integer :: LOOPMAX,i

LOOPMAX = 100                   !* max number of iteration loops 

! Values for entropy
pref = 1000.0_r8           ! mb ref pressure.
eref = 6.106_r8           ! sat p at tfreez (mb)

! Invert the entropy equation -- use Newton's method

Ts = Tfg                  ! Better first guess based on Tprofile from conv.

converge: do i=0, LOOPMAX

   L = rl - (cpliq - cpwv)*(Ts-tfreez) 

   esat = c1*exp(c2*(Ts-tfreez)/(c3+Ts-tfreez)) ! Bolton (eq. 10)
   qsat = eps1*esat/(p-esat)     
   qv = min(qt,qsat) 
   e = qv*p / (eps1 +qv)  ! Bolton (eq. 16)
   fs1 = (cpres + qt*cpliq)*log( Ts/tfreez ) - rgas*log( (p-e)/pref ) + &
        L*qv/Ts - qv*rh2o*log(qv/qsat) - s
   
   L = rl - (cpliq - cpwv)*(Ts-1._r8-tfreez)         

   esat = c1*exp(c2*(Ts-1._r8-tfreez)/(c3+Ts-1._r8-tfreez))
   qsat = eps1*esat/(p-esat)  
   qv = min(qt,qsat) 
   e = qv*p / (eps1 +qv)
   fs2 = (cpres + qt*cpliq)*log( (Ts-1._r8)/tfreez ) - rgas*log( (p-e)/pref ) + &
        L*qv/(Ts-1._r8) - qv*rh2o*log(qv/qsat) - s 
   
   dTs = fs1/(fs2 - fs1)
   Ts  = Ts+dTs
!-------------------------- zhh debug 2013.03.06 -------------------------
!   if (lchnk==145 .or. lchnk==913 .or. lchnk==1041 .or. lchnk==1937) then
!      if (rcall==2) then
!         
!      end if
!   end if
!-------------------------- zhh debug 2013.03.06 -------------------------
   if (abs(dTs).lt.0.001_r8) exit converge
!!   if (abs(dTs).lt.0.01_r8) exit converge
   if (i .eq. LOOPMAX - 1) then

!zmh     
     Ts = Tfg
     exit converge
!zmh


      this_lat = get_rlat_p(lchnk, icol)*57.296_r8
      this_lon = get_rlon_p(lchnk, icol)*57.296_r8
      write(iulog,*) '*** ZM_CONV: IENTROPY: Failed and about to exit, info follows ****'
      write(iulog,100) 'ZM_CONV: IENTROPY. Details: call#,lchnk,icol= ',rcall,lchnk,icol, &
       ' lat: ',this_lat,' lon: ',this_lon, &
       ' P(mb)= ', p, ' Tfg(K)= ', Tfg, ' qt(g/kg) = ', 1000._r8*qt, &
       ' qsat(g/kg) = ', 1000._r8*qsat,', s(J/kg) = ',s
      call endrun('**** ZM_CONV IENTROPY: Tmix did not converge ****')
   end if
enddo converge

! Replace call to satmixutils.

esat = c1*exp(c2*(Ts-tfreez)/(c3+Ts-tfreez))
qsat=eps1*esat/(p-esat)

qv = min(qt,qsat)                             !       /* check for saturation */
T = Ts 

 100    format (A,I1,I4,I4,7(A,F6.2))

return
end SUBROUTINE ientropy

!============================================i===============================
subroutine ent_dynamics (ncol,zf,dz,jb,jb1,jb2,jt,jt2,dzb,dzt,&
        ent_turb, det_turb,eps00, ent,det,idir)
!    
! ---------------------------------------------------------------------------
! Author: Minghua Zhang 2017-06-01
!
! This subroutine calculates the entrainment and detainment rates of plumes 
! by using geometrc variables of the plume and a closure parametrer eps00    
! The geometry of end layers: 
!     from jb2 to jb1 is the transition from 0 to 1 of mass flux
!     from jt2 to jt  is the transition from 1 to 0 of mass flux
!     jb2 >= jb1 > jt2 >=jt1
! The geometry of the plume trunk is from jb1 to jt2 in which dzb fraction is
!     entraining, dzb fraction is detraing
!
! At the moment, the organized entrainment and detaint profiles are simply 
! specified as function of normalized z with fractional thickness dzb, dzt 
!           
! eps00 is the closure rate from w_dynamics
!
! idir = 1 is for updraft, -1 is for downdraft
! ---------------------------------------------------------------------------

 integer,  intent(in) :: ncol
 integer,  intent(in) :: idir  !1: upward plume, -1 downward plume
 real(r8), intent(in) :: zf(pcols,pver+1)  
 real(r8), intent(in) :: dz(pcols,pver)  
 integer,  intent(in) :: jt(pcols)
 integer,  intent(in) :: jt2(pcols)
 integer,  intent(in) :: jb(pcols)
 integer,  intent(in) :: jb1(pcols)
 integer,  intent(in) :: jb2(pcols)
 real(r8), intent(in) :: dzt(pcols)
 real(r8), intent(in) :: dzb(pcols)
 real(r8), intent(in) :: ent_turb(pcols,pver)  
 real(r8), intent(in) :: det_turb(pcols,pver)  
 real(r8), intent(in) :: eps00(pcols)
 !
 real(r8), intent(out) :: ent(pcols,pver)  
 real(r8), intent(out) :: det(pcols,pver)  

 !local 
 integer i,k
 real(r8) zz(pver),zz2(pver)
 real zt,zb,zt2,zb2 ,dzj,ents,dets,ratio,w11

  ent(:,:) = 0.0_r8
  det(:,:) = 0.0_r8

  do i = 1, ncol

   zz(1:pver) = zf(i,1:pver)-dz(i,1:pver)/2.0_r8

   zt = zz(jt2(i))   
   zb = zz(jb1(i))
    
   zt2 = zt - dzt(i)*(zt-zb)   !altitude where detrain occurs
   zb2 = zb + dzb(i)*(zt-zb)   !altitude where entrainment vanishes

   dzj = max(abs((zb2 - zb)),10._r8) * dble(idir)
   zz2 = (zz-zb)/dzj           ! normalized height relative to launching bot level 
   do k = jt2(i),jb1(i),idir
    if( (zz2(k)+1.e-5 .ge. 0._r8) .and. (zz2(k) .le. 1._r8))then

      w11  =  zz2(k)/zmh_ent_ze  !new

      ent(i,k) = exp(-w11*w11)
    endif
   enddo 

   dzj = max(abs((zt - zt2)),10._r8) * dble(idir)
   zz2 = (zt-zz)/dzj           !normalized height relative to detaining top level
   do k = jt2(i),jb1(i),idir
    if( (zz2(k)+1.e-5 .ge. 0._r8) .and. (zz2(k) .le. 1._r8))then
      det(i,k) = 1._r8 - zz2(k)
    endif
   enddo 

   ! force total organized detainment equals to total entrainment
   ents = sum(  ent(i,:) * dz(i,:)  )
   dets = sum(  det(i,:) * dz(i,:)  )

   ratio = ents/max(dets, 1.e-5_r8)
   det(i,:) = det(i,:)*ratio            

   ! total entraiment
   do k=jt2(i),jb1(i),idir
    ent(i,k) = ent(i,k)*eps00(i) + ent_turb(i,k)
    det(i,k) = det(i,k)*eps00(i) + det_turb(i,k)
   enddo 

   ! entrainment and detrainment in the lauching layers (jb1 to jb2) and 
   ! detraining layer (jt2 to jt)

   dzj      = zf(i,jb1(i)) - zf(i,jb2(i)+1) !new
   if(idir .eq. -1)then
     dzj      = abs( zf(i,jb1(i)+1) - zf(i,jb2(i)) )!new
   endif

   do k=jb1(i)+idir,jb2(i),idir
    ent(i,k) = abs(1._r8/dzj)
    det(i,k) = 0.0_r8
   enddo

   dzj      = zf(i,jt(i)) -  zf(i,jt2(i)+1) !new
   if(idir .eq. -1)then
    dzj      = abs( zf(i,jt(i)+1) -  zf(i,jt2(i)) )!new
   endif
 
   do k=jt(i),jt2(i)-idir,idir
    ent(i,k) = 0.0_r8
    det(i,k) = abs(1._r8/dzj)
   enddo

  enddo

!  write(*,*)
!  write(*,*)'jb,jb1,jb2,jt,jt2',jb,jb1,jb2,jt,jt2
!  write(*,*)'ent',ent
!  write(*,*)'det',det

return
end subroutine ent_dynamics

!================================================================================

subroutine w_dynamics (ncol,limcnv,z,dz,buoy,w_init,ent,jb,k100,&
        jt,jt2,jtb,w_up, eps00)
!
! ---------------------------------------------------------------------------
! Author: Minghua Zhang 2017-06-01
! This subroutine calculates plume vertical velocity, plume top, level of maximun
! buoyance by using bouyance, initial vertical velocity, and entrainment closure 
!

 integer,  intent(in) :: ncol,limcnv
 !integer,  intent(in) :: idir  !1: upward plume, -1 downward plume
 real(r8), intent(in) :: z(pcols,pver)  
 real(r8), intent(in) :: dz(pcols,pver)  
 real(r8), intent(in) :: buoy(pcols,pver)
 real(r8), intent(in) :: w_init(pcols)
 real(r8), intent(in) :: ent(pcols,pver)
 integer,  intent(in) :: jb(pcols)     ! launch level
 integer,  intent(in) :: k100(pcols)   ! launch level bound
 !
 real(r8), intent(out) :: w_up(pcols,pver)  
 integer,  intent(inout):: jt(pcols) !zero vertical w
 integer,  intent(out) :: jt2(pcols) !zero buoyance
 integer,  intent(out) :: jtb(pcols) !maximum B
 real(r8), intent(out) :: eps00(pcols) !entrainment closure parameter 

 !local 
 integer i,k,kmaxb

 real(r8) alpha
 real(r8) dp_radius

 real(r8) work(pver)
 real bmax,wmax,dumb,eps,w2

  alpha = 0.05_r8
  dp_radius = 10.0e3_r8

  w_up(:,:)  = 0.0_r8
  eps00(:)   = 0.0_r8
  jt2(:)     = jt(:)
  jtb(:)     = jt(:)
  
  do i = 1, ncol


   kmaxb = get_maxk(buoy(i,:),pver,jt(i),k100(i)) 
   kmaxb = min(kmaxb,jb(i)-1)

   bmax  = buoy(i,kmaxb) 

   if(bmax < 0._r8) cycle
    
   w_up(i,jb(i)) = w_init(i)

   jt2(i) = jb(i)-1
   do k= jb(i)-1, limcnv, -1 
     dumb = buoy(i,k)
     if((k>kmaxb) .and. (dumb< 0._r8))then
        dumb = 0._r8
     endif
     eps = min(ent(i,k)*dz(i,k),1.0_r8) 

     w2 =  (1._r8 - eps)*w_up(i,k+1)*w_up(i,k+1) &
        + alpha*dumb*dz(i,k)

    if(w2 <  0.04_r8)then 
      exit
    endif
     w_up(i,k) = sqrt(w2)

   enddo  !k

   jt(i) = k+1

   do k=jt(i),pver-1
      if(buoy(i,k)>0)then
         jt2(i) = k
        exit
      endif
   enddo  
   jtb(i) = kmaxb

   jt2(i) = max(jt(i)+1,jt2(i))
   jtb(i) = max(kmaxb, jt2(i))

   ! entrainment closure

   wmax  = w_up(i,kmaxb)

   dumb = z(i,kmaxb) - z(i,jb(i))
   wmax = sum(w_up(i,kmaxb:jb(i)) * dz(i,kmaxb:jb(i)) )/dumb !mean w below Bmax

   wmax  = max(wmax,  1.0_r8)
   wmax  = min(wmax,  15._r8)

   do k=kmaxb,jb(i)
      work(k) = max(buoy(i,k),0._r8)
   enddo
   bmax = sum(work(kmaxb:jb(i))* dz(i,kmaxb:jb(i)) )/dumb   !mean B below Bmax

   w2 = min(bmax, 0.4_r8)
   w2 = max(bmax, 0.01_r8)
     
   ! closure alpha*sqrt(B*H)/(R*W)
   eps00(i) = sqrt(w2*(z(i,jt2(i))-z(i,jb(i))))/dp_radius/wmax * zmh_ent_alpha0
!   print*,''
!  print*,'eps00-1',eps00(i),zmh_ent_alpha0
   eps00(i) = max(eps00(i), 1.0e-5_r8)
   eps00(i) = min(eps00(i), 1.0e-3_r8)
!  print*,'eps00-2',eps00(i)


 enddo   !i
   
return
end subroutine w_dynamics

! ==================================
subroutine closure(lchnk   , &
                   q       ,t       ,p       ,z       ,s       , &
                   tp      ,qs      ,qu      ,su      ,mc      , &
                   du      ,mu      ,md      ,qd      ,sd      , &
                   qhat    ,shat    ,dp      ,qstp    ,zf      , &
                   ql      ,dsubcld ,mb      ,cape    ,tl      , &
                   lcl     ,lel     ,jt      ,mx      ,il1g    , &
                   il2g    ,rd      ,grav    ,cp      ,rl      , &
                   msg     ,capelmt ,flag)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! <Say what the routine does> 
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: G. Zhang and collaborators. CCM contact:P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! We expect to release cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
! 
!-----------------------------------------------------------------------
   use dycore,    only: dycore_is, get_resolution

   implicit none

!
!-----------------------------Arguments---------------------------------
!
   integer, intent(in) :: lchnk                 ! chunk identifier

   real(r8), intent(inout) :: q(pcols,pver)        ! spec humidity
   real(r8), intent(inout) :: t(pcols,pver)        ! temperature
   real(r8), intent(inout) :: p(pcols,pver)        ! pressure (mb)
   real(r8), intent(inout) :: mb(pcols)            ! cloud base mass flux
   real(r8), intent(in) :: z(pcols,pver)        ! height (m)
   real(r8), intent(in) :: s(pcols,pver)        ! normalized dry static energy
   real(r8), intent(in) :: tp(pcols,pver)       ! parcel temp
   real(r8), intent(in) :: qs(pcols,pver)       ! sat spec humidity
   real(r8), intent(in) :: qu(pcols,pver)       ! updraft spec. humidity
   real(r8), intent(in) :: su(pcols,pver)       ! normalized dry stat energy of updraft
   real(r8), intent(in) :: mc(pcols,pver)       ! net convective mass flux
   real(r8), intent(in) :: du(pcols,pver)       ! detrainment from updraft
   real(r8), intent(in) :: mu(pcols,pver)       ! mass flux of updraft
   real(r8), intent(in) :: md(pcols,pver)       ! mass flux of downdraft
   real(r8), intent(in) :: qd(pcols,pver)       ! spec. humidity of downdraft
   real(r8), intent(in) :: sd(pcols,pver)       ! dry static energy of downdraft
   real(r8), intent(in) :: qhat(pcols,pver)     ! environment spec humidity at interfaces
   real(r8), intent(in) :: shat(pcols,pver)     ! env. normalized dry static energy at intrfcs
   real(r8), intent(in) :: dp(pcols,pver)       ! pressure thickness of layers
   real(r8), intent(in) :: qstp(pcols,pver)     ! spec humidity of parcel
   real(r8), intent(in) :: zf(pcols,pver+1)     ! height of interface levels
   real(r8), intent(in) :: ql(pcols,pver)       ! liquid water mixing ratio

   logical,  intent(in) :: flag(pcols)
   real(r8), intent(in) :: cape(pcols)          ! available pot. energy of column
   real(r8), intent(in) :: tl(pcols)
   real(r8), intent(in) :: dsubcld(pcols)       ! thickness of subcloud layer

   integer, intent(in) :: lcl(pcols)        ! index of lcl
   integer, intent(in) :: lel(pcols)        ! index of launch leve
   integer, intent(in) :: jt(pcols)         ! top of updraft
   integer, intent(in) :: mx(pcols)         ! base of updraft
!
!--------------------------Local variables------------------------------
!
   real(r8) dtpdt(pcols,pver)
   real(r8) dqsdtp(pcols,pver)
   real(r8) dtmdt(pcols,pver)
   real(r8) dqmdt(pcols,pver)
   real(r8) dboydt(pcols,pver)
   real(r8) thetavp(pcols,pver)
   real(r8) thetavm(pcols,pver)

   real(r8) dtbdt(pcols),dqbdt(pcols),dtldt(pcols)
   real(r8) beta
   real(r8) capelmt(pcols)
   real(r8) cp
   real(r8) dadt(pcols)
   real(r8) debdt
   real(r8) dltaa
   real(r8) eb
   real(r8) grav

   integer i
   integer il1g
   integer il2g
   integer k, kmin, kmax
   integer msg

   real(r8) rd
   real(r8) rl
! change of subcloud layer properties due to convection is
! related to cumulus updrafts and downdrafts.
! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
! to define betau, betad and f(z).
! note that this implies all time derivatives are in effect
! time derivatives per unit cloud-base mass flux, i.e. they
! have units of 1/mb instead of 1/sec.
!
   do i = il1g,il2g
      mb(i) = 0._r8
      eb = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
      dtbdt(i) = (1._r8/dsubcld(i))* (mu(i,mx(i))*(shat(i,mx(i))-su(i,mx(i)))+ &
                  md(i,mx(i))* (shat(i,mx(i))-sd(i,mx(i))))
      dqbdt(i) = (1._r8/dsubcld(i))* (mu(i,mx(i))*(qhat(i,mx(i))-qu(i,mx(i)))+ &
                 md(i,mx(i))* (qhat(i,mx(i))-qd(i,mx(i))))
      debdt = eps1*p(i,mx(i))/ (eps1+q(i,mx(i)))**2*dqbdt(i)
      dtldt(i) = -2840._r8* (3.5_r8/t(i,mx(i))*dtbdt(i)-debdt/eb)/ &
                 (3.5_r8*log(t(i,mx(i)))-log(eb)-4.805_r8)**2
   end do
!
!   dtmdt and dqmdt are cumulus heating and drying.
!
   do k = msg + 1,pver
      do i = il1g,il2g
         dtmdt(i,k) = 0._r8
         dqmdt(i,k) = 0._r8
      end do
   end do
!
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k == jt(i)) then
            dtmdt(i,k) = (1._r8/dp(i,k))*(mu(i,k+1)* (su(i,k+1)-shat(i,k+1)- &
                          rl/cp*ql(i,k+1))+md(i,k+1)* (sd(i,k+1)-shat(i,k+1)))
            dqmdt(i,k) = (1._r8/dp(i,k))*(mu(i,k+1)* (qu(i,k+1)- &
                         qhat(i,k+1)+ql(i,k+1))+md(i,k+1)*(qd(i,k+1)-qhat(i,k+1)))
         end if
      end do
   end do
!
   beta = 0._r8
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k > jt(i) .and. k < mx(i)) then
            dtmdt(i,k) = (mc(i,k)* (shat(i,k)-s(i,k))+mc(i,k+1)* (s(i,k)-shat(i,k+1)))/ &
                         dp(i,k) - rl/cp*du(i,k)*(beta*ql(i,k)+ (1-beta)*ql(i,k+1))
!          dqmdt(i,k)=(mc(i,k)*(qhat(i,k)-q(i,k))
!     1                +mc(i,k+1)*(q(i,k)-qhat(i,k+1)))/dp(i,k)
!     2                +du(i,k)*(qs(i,k)-q(i,k))
!     3                +du(i,k)*(beta*ql(i,k)+(1-beta)*ql(i,k+1))

            dqmdt(i,k) = (mu(i,k+1)* (qu(i,k+1)-qhat(i,k+1)+cp/rl* (su(i,k+1)-s(i,k)))- &
                          mu(i,k)* (qu(i,k)-qhat(i,k)+cp/rl*(su(i,k)-s(i,k)))+md(i,k+1)* &
                         (qd(i,k+1)-qhat(i,k+1)+cp/rl*(sd(i,k+1)-s(i,k)))-md(i,k)* &
                         (qd(i,k)-qhat(i,k)+cp/rl*(sd(i,k)-s(i,k))))/dp(i,k) + &
                          du(i,k)* (beta*ql(i,k)+(1-beta)*ql(i,k+1))
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k >= lel(i) .and. k <= lcl(i)) then
            thetavp(i,k) = tp(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+1.608_r8*qstp(i,k)-q(i,mx(i)))
            thetavm(i,k) = t(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+0.608_r8*q(i,k))
            dqsdtp(i,k) = qstp(i,k)* (1._r8+qstp(i,k)/eps1)*eps1*rl/(rd*tp(i,k)**2)
!
! dtpdt is the parcel temperature change due to change of
! subcloud layer properties during convection.
!
            dtpdt(i,k) = tp(i,k)/ (1._r8+rl/cp* (dqsdtp(i,k)-qstp(i,k)/tp(i,k)))* &
                        (dtbdt(i)/t(i,mx(i))+rl/cp* (dqbdt(i)/tl(i)-q(i,mx(i))/ &
                         tl(i)**2*dtldt(i)))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = ((dtpdt(i,k)/tp(i,k)+1._r8/(1._r8+1.608_r8*qstp(i,k)-q(i,mx(i)))* &
                          (1.608_r8 * dqsdtp(i,k) * dtpdt(i,k) -dqbdt(i))) - (dtmdt(i,k)/t(i,k)+0.608_r8/ &
                          (1._r8+0.608_r8*q(i,k))*dqmdt(i,k)))*grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do
!
   do k = msg + 1,pver
      do i = il1g,il2g
         if (k > lcl(i) .and. k < mx(i)) then
            thetavp(i,k) = tp(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+0.608_r8*q(i,mx(i)))
            thetavm(i,k) = t(i,k)* (1000._r8/p(i,k))** (rd/cp)*(1._r8+0.608_r8*q(i,k))
!
! dboydt is the integrand of cape change.
!
            dboydt(i,k) = (dtbdt(i)/t(i,mx(i))+0.608_r8/ (1._r8+0.608_r8*q(i,mx(i)))*dqbdt(i)- &
                          dtmdt(i,k)/t(i,k)-0.608_r8/ (1._r8+0.608_r8*q(i,k))*dqmdt(i,k))* &
                          grav*thetavp(i,k)/thetavm(i,k)
         end if
      end do
   end do

!
! buoyant energy change is set to 2/3*excess cape per 3 hours
!
   dadt(il1g:il2g)  = 0._r8
   kmin = minval(lel(il1g:il2g))
   kmax = maxval(mx(il1g:il2g)) - 1
   do k = kmin, kmax
      do i = il1g,il2g
         if ( k >= lel(i) .and. k <= mx(i) - 1) then
            dadt(i) = dadt(i) + dboydt(i,k)* (zf(i,k)-zf(i,k+1))
         endif
      end do
   end do
   do i = il1g,il2g
      !dltaa = -1._r8* (cape(i)-capelmt)
      !dltaa = -1._r8* max((cape(i)-capelmt(i)),0._r8)
      dltaa = -1._r8* max(cape(i)-50._r8,0._r8)
      if (dadt(i) /= 0._r8) mb(i) = max(dltaa/zmh_tau/dadt(i),0._r8)
   end do
!
   return
end subroutine closure


!====================================================
  real(r8) function get_mink(d,np,k1,k2)
!
!====================================================
! to get the index of minimun value of d betwwn k1, k2
!
! Author: Minghua Zhang 2017-06-03
! get index of minimum value of d(k1:k2) 
! ----------------------------------------------------

 real(8)   d(np)
 integer   np,k1,k2
 integer   k3, k
 real(8)   dmn

   dmn = 1.0e10_r8
   k3   = k1
   do k = k1,k2
    if(dmn > d(k))then
      dmn = d(k)
      k3 = k
    endif
   enddo 

   get_mink = k3 

return

end function get_mink

!====================================================
real(r8) function get_maxk(d,np,k1,k2)
! 
!====================================================
! to get the index of maximum value of d(k1:k2)
!
! Author: Minghua Zhang 2017-06-03
!
! ----------------------------------------------------

 real(8)   d(np)
 integer   np,k1,k2
 integer   k3, k
 real(8)   dmn
   
   dmn = -1.0e10_r8
   k3   = k1
   do k = k1,k2
    if(dmn < d(k))then
      dmn = d(k)
      k3 = k
    endif
   enddo 
   get_maxk = k3
return

end function get_maxk

!===============================================================================
subroutine zyx1_conv_evap(ncol,lchnk, &
!zmh     t,pmid,pdel,q, &
     t,pmid,pdel,q, lat, lon,&
     tend_s, tend_s_snwprd, tend_s_snwevmlt, tend_q, &
     prdprec, cldfrc, deltat,  &
     prec, snow, ntprprd, ntsnprd, flxprec, flxsnow )

!-----------------------------------------------------------------------
! Compute tendencies due to evaporation of rain from ZM scheme
!--
! Compute the total precipitation and snow fluxes at the surface.
! Add in the latent heat of fusion for snow formation and melt, since it not dealt with
! in the Zhang-MacFarlane parameterization.
! Evaporate some of the precip directly into the environment using a Sundqvist type algorithm
!-----------------------------------------------------------------------

    use wv_saturation,  only: aqsat
    use phys_grid, only: get_rlat_all_p

!------------------------------Arguments--------------------------------
    integer,intent(in) :: ncol, lchnk             ! number of columns and chunk index
    real(r8),intent(in), dimension(pcols,pver) :: t          ! temperature (K)
    real(r8),intent(in), dimension(pcols,pver) :: pmid       ! midpoint pressure (Pa) 
    real(r8),intent(in), dimension(pcols,pver) :: pdel       ! layer thickness (Pa)
    real(r8),intent(in), dimension(pcols,pver) :: q          ! water vapor (kg/kg)
    real(r8),intent(inout), dimension(pcols,pver) :: tend_s     ! heating rate (J/kg/s)
    real(r8),intent(inout), dimension(pcols,pver) :: tend_q     ! water vapor tendency (kg/kg/s)
    real(r8),intent(out  ), dimension(pcols,pver) :: tend_s_snwprd ! Heating rate of snow production
    real(r8),intent(out  ), dimension(pcols,pver) :: tend_s_snwevmlt ! Heating rate of evap/melting of snow
    


    real(r8), intent(in   ) :: prdprec(pcols,pver)! precipitation production (kg/ks/s)
    real(r8), intent(in   ) :: cldfrc(pcols,pver) ! cloud fraction
    real(r8), intent(in   ) :: deltat             ! time step
    real(r8), intent(in   ) :: lat(pcols)        ! 
    real(r8), intent(in   ) :: lon(pcols)        ! 

    real(r8), intent(inout) :: prec(pcols)        ! Convective-scale preciptn rate
    real(r8), intent(out)   :: snow(pcols)        ! Convective-scale snowfall rate
!
!---------------------------Local storage-------------------------------

    real(r8) :: est    (pcols,pver)    ! Saturation vapor pressure
    real(r8) :: fice   (pcols,pver)    ! ice fraction in precip production
    real(r8) :: fsnow_conv(pcols,pver) ! snow fraction in precip production
    real(r8) :: qsat   (pcols,pver)    ! saturation specific humidity
    real(r8),intent(out) :: flxprec(pcols,pverp)   ! Convective-scale flux of precip at interfaces (kg/m2/s)
    real(r8),intent(out) :: flxsnow(pcols,pverp)   ! Convective-scale flux of snow   at interfaces (kg/m2/s)
    real(r8),intent(out) :: ntprprd(pcols,pver)    ! net precip production in layer
    real(r8),intent(out) :: ntsnprd(pcols,pver)    ! net snow production in layer
    real(r8) :: work1                  ! temp variable (pjr)
    real(r8) :: work2                  ! temp variable (pjr)

    real(r8) :: evpvint(pcols)         ! vertical integral of evaporation
    real(r8) :: evpprec(pcols)         ! evaporation of precipitation (kg/kg/s)
    real(r8) :: evpsnow(pcols)         ! evaporation of snowfall (kg/kg/s)
    real(r8) :: snowmlt(pcols)         ! snow melt tendency in layer
    real(r8) :: flxsntm(pcols)         ! flux of snow into layer, after melting

    real(r8) :: evplimit               ! temp variable for evaporation limits
    real(r8) :: rlat(pcols)

    integer :: i,k                     ! longitude,level indices
    real(r8):: tmp1

!-----------------------------------------------------------------------

! convert input precip to kg/m2/s
    prec(:ncol) = prec(:ncol)*1000._r8

! determine saturation vapor pressure
    call aqsat (t    ,pmid  ,est    ,qsat    ,pcols   , &
         ncol ,pver  ,1       ,pver    )

! determine ice fraction in rain production (use cloud water parameterization fraction at present)
    call cldwat_fice(ncol, t, fice, fsnow_conv)

! zero the flux integrals on the top boundary
    flxprec(:ncol,1) = 0._r8
    flxsnow(:ncol,1) = 0._r8
    evpvint(:ncol)   = 0._r8

    do k = 1, pver
       do i = 1, ncol

! Melt snow falling into layer, if necessary. 
          if (t(i,k) > tmelt) then
             flxsntm(i) = 0._r8
             snowmlt(i) = flxsnow(i,k) * gravit/ pdel(i,k)
          else
             flxsntm(i) = flxsnow(i,k)
             snowmlt(i) = 0._r8
          end if

! relative humidity depression must be > 0 for evaporation
          evplimit = max(1._r8 - q(i,k)/qsat(i,k), 0._r8)

! total evaporation depends on flux in the top of the layer
! flux prec is the net production above layer minus evaporation into environmet

! zmh ke test
! 335
!          evpprec(i) = ke * (1._r8 - cldfrc(i,k)) * evplimit * sqrt(flxprec(i,k)) & 
!               *(1.+10.*(pmid(i,k)/pmid(i,pver))**4 )/11.
!337
          evpprec(i) = ke * (1._r8 - cldfrc(i,k)) * evplimit * sqrt(flxprec(i,k))  !ke=2.e-6

! zmh !!!
!!          evpprec(i) = evpprec(i)/(pmid(i,k)/287.6/t(i,k))
! zmh !!!
!**********************************************************
!!          evpprec(i) = 0.    ! turn off evaporation for now
!**********************************************************

! Don't let evaporation supersaturate layer (approx). Layer may already be saturated.
! Currently does not include heating/cooling change to qsat
          evplimit   = max(0._r8, (qsat(i,k)-q(i,k)) / deltat)

! Don't evaporate more than is falling into the layer - do not evaporate rain formed
! in this layer but if precip production is negative, remove from the available precip
! Negative precip production occurs because of evaporation in downdrafts.
!!$          evplimit   = flxprec(i,k) * gravit / pdel(i,k) + min(prdprec(i,k), 0.)
          evplimit   = min(evplimit, flxprec(i,k) * gravit / pdel(i,k))

! Total evaporation cannot exceed input precipitation
          evplimit   = min(evplimit, (prec(i) - evpvint(i)) * gravit / pdel(i,k))

          evpprec(i) = min(evplimit, evpprec(i))

! evaporation of snow depends on snow fraction of total precipitation in the top after melting
          if (flxprec(i,k) > 0._r8) then
!            evpsnow(i) = evpprec(i) * flxsntm(i) / flxprec(i,k)
!            prevent roundoff problems
             work1 = min(max(0._r8,flxsntm(i)/flxprec(i,k)),1._r8)
             evpsnow(i) = evpprec(i) * work1
          else
             evpsnow(i) = 0._r8
          end if

! vertically integrated evaporation
          evpvint(i) = evpvint(i) + evpprec(i) * pdel(i,k)/gravit

! net precip production is production - evaporation
          ntprprd(i,k) = prdprec(i,k) - evpprec(i)
! net snow production is precip production * ice fraction - evaporation - melting
!pjrworks ntsnprd(i,k) = prdprec(i,k)*fice(i,k) - evpsnow(i) - snowmlt(i)
!pjrwrks2 ntsnprd(i,k) = prdprec(i,k)*fsnow_conv(i,k) - evpsnow(i) - snowmlt(i)
! the small amount added to flxprec in the work1 expression has been increased from 
! 1e-36 to 8.64e-11 (1e-5 mm/day).  This causes the temperature based partitioning
! scheme to be used for small flxprec amounts.  This is to address error growth problems.
#ifdef PERGRO
          work1 = min(max(0._r8,flxsnow(i,k)/(flxprec(i,k)+8.64e-11_r8)),1._r8)
#else
          if (flxprec(i,k).gt.0._r8) then
             work1 = min(max(0._r8,flxsnow(i,k)/flxprec(i,k)),1._r8)
          else
             work1 = 0._r8
          endif
#endif
          work2 = max(fsnow_conv(i,k), work1)
          if (snowmlt(i).gt.0._r8) work2 = 0._r8
!         work2 = fsnow_conv(i,k)
          ntsnprd(i,k) = prdprec(i,k)*work2 - evpsnow(i) - snowmlt(i)
          tend_s_snwprd  (i,k) = prdprec(i,k)*work2*latice
          tend_s_snwevmlt(i,k) = - ( evpsnow(i) + snowmlt(i) )*latice

! precipitation fluxes
          flxprec(i,k+1) = flxprec(i,k) + ntprprd(i,k) * pdel(i,k)/gravit
          flxsnow(i,k+1) = flxsnow(i,k) + ntsnprd(i,k) * pdel(i,k)/gravit

! protect against rounding error
          flxprec(i,k+1) = max(flxprec(i,k+1), 0._r8)
          flxsnow(i,k+1) = max(flxsnow(i,k+1), 0._r8)
! more protection (pjr)
!         flxsnow(i,k+1) = min(flxsnow(i,k+1), flxprec(i,k+1))

! heating (cooling) and moistening due to evaporation 
! - latent heat of vaporization for precip production has already been accounted for
! - snow is contained in prec
          tend_s(i,k)   =-evpprec(i)*latvap + ntsnprd(i,k)*latice
          tend_q(i,k) = evpprec(i)
       end do
    end do

! set output precipitation rates (m/s)
    prec(:ncol) = flxprec(:ncol,pver+1) / 1000._r8
    snow(:ncol) = flxsnow(:ncol,pver+1) / 1000._r8

!**********************************************************
!!$    tend_s(:ncol,:)   = 0.      ! turn heating off
!**********************************************************


  end subroutine zyx1_conv_evap

! ================================================

subroutine convtran(lchnk   , &
                    doconvtran,q       ,ncnst   ,mu      ,md      , &
                    du      ,eu      ,ed      ,dp      ,dsubcld , &
                    jt      ,mx      ,ideep   ,il1g    ,il2g    , &
                    nstep   ,fracis  ,dqdt    ,dpdry   )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective transport of trace species
!
! Mixing ratios may be with respect to either dry or moist air
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use constituents,    only: cnst_get_type_byind
   use ppgrid
   use abortutils, only: endrun

   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncnst                 ! number of tracers to transport
   logical, intent(in) :: doconvtran(ncnst)     ! flag for doing convective transport
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Tracer array including moisture
   real(r8), intent(in) :: mu(pcols,pver)       ! Mass flux up
   real(r8), intent(in) :: md(pcols,pver)       ! Mass flux down
   real(r8), intent(in) :: du(pcols,pver)       ! Mass detraining from updraft
   real(r8), intent(in) :: eu(pcols,pver)       ! Mass entraining from updraft
   real(r8), intent(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
   real(r8), intent(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
   real(r8), intent(in) :: dsubcld(pcols)       ! Delta pressure from cloud base to sfc
   real(r8), intent(in) :: fracis(pcols,pver,ncnst) ! fraction of tracer that is insoluble

   integer, intent(in) :: jt(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: mx(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: ideep(pcols)      ! Gathering array
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index

   real(r8), intent(in) :: dpdry(pcols,pver)       ! Delta pressure between interfaces


! input/output

   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array

!--------------------------Local Variables------------------------------

   integer i                 ! Work index
   integer k                 ! Work index
   integer kbm               ! Highest altitude index of cloud base
   integer kk                ! Work index
   integer kkp1              ! Work index
   integer km1               ! Work index
   integer kp1               ! Work index
   integer ktm               ! Highest altitude index of cloud top
   integer m                 ! Work index

   real(r8) cabv                 ! Mix ratio of constituent above
   real(r8) cbel                 ! Mix ratio of constituent below
   real(r8) cdifr                ! Normalized diff between cabv and cbel
   real(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
   real(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
   real(r8) const(pcols,pver)    ! Gathered tracer array
   real(r8) fisg(pcols,pver)     ! gathered insoluble fraction of tracer
   real(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
   real(r8) dcondt(pcols,pver)   ! Gathered tend array
   real(r8) small                ! A small number
   real(r8) mbsth                ! Threshold for mass fluxes
   real(r8) mupdudp              ! A work variable
   real(r8) minc                 ! A work variable
   real(r8) maxc                 ! A work variable
   real(r8) fluxin               ! A work variable
   real(r8) fluxout              ! A work variable
   real(r8) netflux              ! A work variable

   real(r8) dutmp(pcols,pver)       ! Mass detraining from updraft
   real(r8) eutmp(pcols,pver)       ! Mass entraining from updraft
   real(r8) edtmp(pcols,pver)       ! Mass entraining from downdraft
   real(r8) dptmp(pcols,pver)    ! Delta pressure between interfaces
!-----------------------------------------------------------------------
!
   small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

! Loop ever each constituent
   do m = 2, ncnst
      if (doconvtran(m)) then

         if (cnst_get_type_byind(m).eq.'dry') then
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dpdry(i,k)
                  dutmp(i,k) = du(i,k)*dp(i,k)/dpdry(i,k)
                  eutmp(i,k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                  edtmp(i,k) = ed(i,k)*dp(i,k)/dpdry(i,k)
               end do
            end do
         else
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dp(i,k)
                  dutmp(i,k) = du(i,k)
                  eutmp(i,k) = eu(i,k)
                  edtmp(i,k) = ed(i,k)
               end do
            end do
         endif
!        dptmp = dp

! Gather up the constituent and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
               fisg(i,k) = fracis(ideep(i),k,m)
            end do
         end do

! From now on work only with gathered data

! Interpolate environment tracer values to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc < 0) then
                  cdifr = 0._r8
               else
                  cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
               endif

! If the two layers differ significantly use a geometric averaging
! procedure
               if (cdifr > 1.E-6_r8) then
                  cabv = max(const(i,km1),maxc*1.e-12_r8)
                  cbel = max(const(i,k),maxc*1.e-12_r8)
                  chat(i,k) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))
               end if

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._r8

            end do
         end do

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         do i = il1g,il2g
            mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
            if (mupdudp > mbsth) then
               conu(i,kk) = (+eutmp(i,kk)*fisg(i,kk)*const(i,kk)*dptmp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1))/md(i,k)
            endif
         end do

! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
               if (mupdudp > mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eutmp(i,kk)*fisg(i,kk)* &
                                  const(i,kk)*dptmp(i,kk) )/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-edtmp(i,km1)*fisg(i,km1)*const(i,km1) &
                                  *dptmp(i,km1) )/md(i,k)
               endif
            end do
         end do


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g

! version 1 hard to check for roundoff errors
!               dcondt(i,k) =
!     $                  +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))
!     $                    -mu(i,k)*   (conu(i,k)-chat(i,k))
!     $                    +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))
!     $                    -md(i,k)*   (cond(i,k)-chat(i,k))
!     $                   )/dp(i,k)

! version 2 hard to limit fluxes
!               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)
!     $                 -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
!               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)
!     $                 -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

! version 3 limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
               fluxin =  mu(i,kp1)*conu(i,kp1)+ mu(i,k)*min(chat(i,k),const(i,km1)) &
                         -(md(i,k)  *cond(i,k) + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
               fluxout = mu(i,k)*conu(i,k) + mu(i,kp1)*min(chat(i,kp1),const(i,k)) &
                         -(md(i,kp1)*cond(i,kp1) + md(i,k)*min(chat(i,k),const(i,k)))

               netflux = fluxin - fluxout
               if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                  netflux = 0._r8
               endif
               dcondt(i,k) = netflux/dptmp(i,k)
            end do
         end do
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!DIR$ NOINTERCHANGE
         do k = kbm,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k == mx(i)) then

! version 1
!                  dcondt(i,k) = (1./dsubcld(i))*
!     $              (-mu(i,k)*(conu(i,k)-chat(i,k))
!     $               -md(i,k)*(cond(i,k)-chat(i,k))
!     $              )

! version 2
!                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
!                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
! version 3
                  fluxin =  mu(i,k)*min(chat(i,k),const(i,km1)) - md(i,k)*cond(i,k)
                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*min(chat(i,k),const(i,k))

                  netflux = fluxin - fluxout
                  if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                     netflux = 0._r8
                  endif
!                  dcondt(i,k) = netflux/dsubcld(i)
                  dcondt(i,k) = netflux/dptmp(i,k)
               else if (k > mx(i)) then
!                  dcondt(i,k) = dcondt(i,k-1)
                  dcondt(i,k) = 0._r8
               end if
            end do
         end do

! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8
         do k = 1,pver
            kp1 = min(pver,k+1)
!DIR$ CONCURRENT
            do i = il1g,il2g
               dqdt(ideep(i),k,m) = dcondt(i,k)
            end do
         end do

      end if      ! for doconvtran

   end do

   return
end subroutine convtran

!=========================================================================================

subroutine momtran(lchnk, ncol, &
                    domomtran,q       ,ncnst   ,mu      ,md    , &
                    du      ,eu      ,ed      ,dp      ,dsubcld , &
                    jt      ,mx      ,ideep   ,il1g    ,il2g    , &
                    nstep   ,dqdt    ,pguall     ,pgdall, icwu, icwd, dt, seten    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Convective transport of momentum
!
! Mixing ratios may be with respect to either dry or moist air
! 
! Method: 
! Based on the convtran subroutine by P. Rasch
! <Also include any applicable external references.> 
! 
! Author: J. Richter and P. Rasch
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use constituents,    only: cnst_get_type_byind
   use ppgrid
   use abortutils, only: endrun

   implicit none
!-----------------------------------------------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: ncnst                 ! number of tracers to transport
   logical, intent(in) :: domomtran(ncnst)      ! flag for doing convective transport
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Wind array
   real(r8), intent(in) :: mu(pcols,pver)       ! Mass flux up
   real(r8), intent(in) :: md(pcols,pver)       ! Mass flux down
   real(r8), intent(in) :: du(pcols,pver)       ! Mass detraining from updraft
   real(r8), intent(in) :: eu(pcols,pver)       ! Mass entraining from updraft
   real(r8), intent(in) :: ed(pcols,pver)       ! Mass entraining from downdraft
   real(r8), intent(in) :: dp(pcols,pver)       ! Delta pressure between interfaces
   real(r8), intent(in) :: dsubcld(pcols)       ! Delta pressure from cloud base to sfc
   real(r8), intent(in)    :: dt                       !  time step in seconds : 2*delta_t

   integer, intent(in) :: jt(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: mx(pcols)         ! Index of cloud top for each column
   integer, intent(in) :: ideep(pcols)      ! Gathering array
   integer, intent(in) :: il1g              ! Gathered min lon indices over which to operate
   integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
   integer, intent(in) :: nstep             ! Time step index



! input/output

   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array

!--------------------------Local Variables------------------------------

   integer i                 ! Work index
   integer k                 ! Work index
   integer kbm               ! Highest altitude index of cloud base
   integer kk                ! Work index
   integer kkp1              ! Work index
   integer kkm1              ! Work index
   integer km1               ! Work index
   integer kp1               ! Work index
   integer ktm               ! Highest altitude index of cloud top
   integer m                 ! Work index
   integer ii                 ! Work index

   real(r8) cabv                 ! Mix ratio of constituent above
   real(r8) cbel                 ! Mix ratio of constituent below
   real(r8) cdifr                ! Normalized diff between cabv and cbel
   real(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
   real(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
   real(r8) const(pcols,pver)    ! Gathered wind array
   real(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
   real(r8) dcondt(pcols,pver)   ! Gathered tend array
   real(r8) small                ! A small number
   real(r8) mbsth                ! Threshold for mass fluxes
   real(r8) mupdudp              ! A work variable
   real(r8) minc                 ! A work variable
   real(r8) maxc                 ! A work variable
   real(r8) fluxin               ! A work variable
   real(r8) fluxout              ! A work variable
   real(r8) netflux              ! A work variable

   real(r8) momcu                ! constant for updraft pressure gradient term
   real(r8) momcd                ! constant for downdraft pressure gradient term
   real(r8) sum                  ! sum
   real(r8) sum2                  ! sum2
 
   real(r8) mududp(pcols,pver) ! working variable
   real(r8) mddudp(pcols,pver)     ! working variable

   real(r8) pgu(pcols,pver)      ! Pressure gradient term for updraft
   real(r8) pgd(pcols,pver)      ! Pressure gradient term for downdraft

   real(r8),intent(out) ::  pguall(pcols,pver,ncnst)      ! Apparent force from  updraft PG
   real(r8),intent(out) ::  pgdall(pcols,pver,ncnst)      ! Apparent force from  downdraft PG

   real(r8),intent(out) ::  icwu(pcols,pver,ncnst)      ! In-cloud winds in updraft
   real(r8),intent(out) ::  icwd(pcols,pver,ncnst)      ! In-cloud winds in downdraft

   real(r8),intent(out) ::  seten(pcols,pver) ! Dry static energy tendency
   real(r8)                 gseten(pcols,pver) ! Gathered dry static energy tendency

   real(r8)  mflux(pcols,pverp,ncnst)   ! Gathered momentum flux

   real(r8)  wind0(pcols,pver,ncnst)       !  gathered  wind before time step
   real(r8)  windf(pcols,pver,ncnst)       !  gathered  wind after time step
   real(r8) fkeb, fket, ketend_cons, ketend, utop, ubot, vtop, vbot, gset2
   

!-----------------------------------------------------------------------
!

! Initialize outgoing fields
   pguall(:,:,:)     = 0.0_r8
   pgdall(:,:,:)     = 0.0_r8
! Initialize in-cloud winds to environmental wind
   icwu(:ncol,:,:)       = q(:ncol,:,:)
   icwd(:ncol,:,:)       = q(:ncol,:,:)

! Initialize momentum flux and  final winds
   mflux(:,:,:)       = 0.0_r8
   wind0(:,:,:)         = 0.0_r8
   windf(:,:,:)         = 0.0_r8

! Initialize dry static energy

   seten(:,:)         = 0.0_r8
   gseten(:,:)         = 0.0_r8

! Define constants for parameterization

   momcu = 0.4_r8
   momcd = 0.4_r8

   small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

! Loop ever each wind component
   do m = 1, ncnst                    !start at m = 1 to transport momentum
      if (domomtran(m)) then

! Gather up the winds and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
                wind0(i,k,m) = const(i,k)
            end do
         end do


! From now on work only with gathered data

! Interpolate winds to interfaces

         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g

               ! use arithmetic mean
               chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._r8

            end do
         end do


!
! Pressure Perturbation Term
! 

      !Top boundary:  assume mu is zero 

         k=1
         pgu(:il2g,k) = 0.0_r8
         pgd(:il2g,k) = 0.0_r8

         do k=2,pver-1
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g
            
               !interior points

               mududp(i,k) =  ( mu(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) &
                           +  mu(i,kp1) * (const(i,kp1) - const(i,k))/dp(i,k))

               pgu(i,k) = - momcu * 0.5_r8 * mududp(i,k)
                           

               mddudp(i,k) =  ( md(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) &
                           +  md(i,kp1) * (const(i,kp1) - const(i,k))/dp(i,k))

               pgd(i,k) = - momcd * 0.5_r8 * mddudp(i,k)


            end do
         end do

       ! bottom boundary 
       k = pver
       km1 = max(1,k-1)
       do i=il1g,il2g

          mududp(i,k) =   mu(i,k) * (const(i,k)- const(i,km1))/dp(i,km1)
          pgu(i,k) = - momcu *  mududp(i,k)
          
          mddudp(i,k) =   md(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) 

          pgd(i,k) = - momcd * mddudp(i,k)
          
       end do
       

!
! In-cloud velocity calculations
!

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         kkm1 = max(1,kk-1)
         do i = il1g,il2g
            mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
            if (mupdudp > mbsth) then
                 
               conu(i,kk) = (+eu(i,kk)*const(i,kk)*dp(i,kk)+pgu(i,kk)*dp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-ed(i,km1)*const(i,km1)*dp(i,km1))-pgd(i,km1)*dp(i,km1)/md(i,k)
            endif

                        
         end do



! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkm1 = max(1,kk-1)
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
               if (mupdudp > mbsth) then
            
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eu(i,kk)* &
                                  const(i,kk)*dp(i,kk)+pgu(i,kk)*dp(i,kk))/mupdudp
               endif
            end do

         end do


! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                            
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-ed(i,km1)*const(i,km1) &
                                  *dp(i,km1)-pgd(i,km1)*dp(i,km1) )/md(i,k)

               endif
            end do
         end do


         sum = 0._r8
         sum2 = 0._r8


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g
               ii = ideep(i)
	
! version 1 hard to check for roundoff errors
               dcondt(i,k) =  &
                           +(mu(i,kp1)* (conu(i,kp1)-chat(i,kp1)) &
                           -mu(i,k)*   (conu(i,k)-chat(i,k))      &
                           +md(i,kp1)* (cond(i,kp1)-chat(i,kp1)) &
                           -md(i,k)*   (cond(i,k)-chat(i,k)) &
                          )/dp(i,k)

            end do
         end do

  ! dcont for bottom layer
          !
          !DIR$ NOINTERCHANGE
          do k = kbm,pver
             km1 = max(1,k-1)
             do i = il1g,il2g
                if (k == mx(i)) then

                   ! version 1
                   dcondt(i,k) = (1./dp(i,k))*   &  
                        (-mu(i,k)*(conu(i,k)-chat(i,k)) &
                        -md(i,k)*(cond(i,k)-chat(i,k)) &
                        )
                end if
             end do
          end do

! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8

         do k = 1,pver
            do i = il1g,il2g
               ii = ideep(i)
               dqdt(ii,k,m) = dcondt(i,k)
    ! Output apparent force on the mean flow from pressure gradient
               pguall(ii,k,m) = -pgu(i,k)
               pgdall(ii,k,m) = -pgd(i,k)
               icwu(ii,k,m)   =  conu(i,k)
               icwd(ii,k,m)   =  cond(i,k)
            end do
         end do

          ! Calculate momentum flux in units of mb*m/s2 

          do k = ktm,pver
             do i = il1g,il2g
                ii = ideep(i)
                mflux(i,k,m) = &
                     -mu(i,k)*   (conu(i,k)-chat(i,k))      &
                     -md(i,k)*   (cond(i,k)-chat(i,k))
             end do
          end do


          ! Calculate winds at the end of the time step 

          do k = ktm,pver
             do i = il1g,il2g
                ii = ideep(i)
                km1 = max(1,k-1)
                kp1 = k+1
                windf(i,k,m) = const(i,k)    -   (mflux(i,kp1,m) - mflux(i,k,m)) * dt /dp(i,k)

             end do
          end do

       end if      ! for domomtran
   end do

 ! Need to add an energy fix to account for the dissipation of kinetic energy
    ! Formulation follows from Boville and Bretherton (2003)
    ! formulation by PJR

    do k = ktm,pver
       km1 = max(1,k-1)
       kp1 = min(pver,k+1)
       do i = il1g,il2g

          ii = ideep(i)

          ! calculate the KE fluxes at top and bot of layer 
          ! based on a discrete approximation to b&b eq(35) F_KE = u*F_u + v*F_v at interface
          utop = (wind0(i,k,1)+wind0(i,km1,1))/2.
          vtop = (wind0(i,k,2)+wind0(i,km1,2))/2.
          ubot = (wind0(i,kp1,1)+wind0(i,k,1))/2.
          vbot = (wind0(i,kp1,2)+wind0(i,k,2))/2.
          fket = utop*mflux(i,k,1)   + vtop*mflux(i,k,2)    ! top of layer
          fkeb = ubot*mflux(i,k+1,1) + vbot*mflux(i,k+1,2)  ! bot of layer

          ! divergence of these fluxes should give a conservative redistribution of KE
          ketend_cons = (fket-fkeb)/dp(i,k)

          ! tendency in kinetic energy resulting from the momentum transport
          ketend = ((windf(i,k,1)**2 + windf(i,k,2)**2) - (wind0(i,k,1)**2 + wind0(i,k,2)**2))*0.5/dt

          ! the difference should be the dissipation
          gset2 = ketend_cons - ketend
          gseten(i,k) = gset2

       end do

    end do

    ! Scatter dry static energy to full array
    do k = 1,pver
       do i = il1g,il2g
          ii = ideep(i)
          seten(ii,k) = gseten(i,k)

       end do
    end do

   return
end subroutine momtran

!=========================================================================================


end module zyx1_conv



