! ###########################################################################################
!
! ###########################################################################################
module gfs_mmm_sfclayrev_pre
  use machine, only: kind_phys

  implicit none

  public gfs_mmm_sfclayrev_pre_run

contains

!> \section arg_table_gfs_mmm_sfclayrev_pre_run
!! \htmlinclude gfs_mmm_sfclayrev_pre_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine gfs_mmm_sfclayrev_pre_run(do_mmm_sfclayrev, nCol, zorl, znt, wetness, mavail, pbl_regime, fm, fh, cpm, slimask, xland, &
             wet, dry, icy, tsfc_wat, tsfc_lnd, tsfc_ice, qsfc_wat, qsfc_lnd, qsfc_ice, tsfc, qsfc, &
             ep1, t, q, ps, con_rd, con_cp, hfx_wat, hfx_lnd, hfx_ice, hfx, qfx_wat, qfx_lnd, qfx_ice,qfx, scm_force_flux, &
             qgh, gz1oz0, water_depth, shalwater_depth, svp1, svp2, svp3, svpt0, p1000mb, its, ite, errmsg, errflg)

    ! Input
    logical,                       intent(in)           :: &
        do_mmm_sfclayrev  ! Flag for MMM SFCLAYREV scheme
    logical,                       intent(out)          :: &
        scm_force_flux    ! Flag for not computing surface fluxes but prescribing
    logical,         dimension(:), intent(in)           :: &
        wet,            & ! Flag indicating presence of some ocean or lake surface area fraction
        dry,            & ! Flag indicating presence of some land surface area fraction
        icy               ! Flag indicating presence of some sea ice surface area fraction 

    integer,                       intent(in)           :: &
        nCol              ! Number of horizontal gridpoints
    integer,         dimension(:), intent(in)           :: &
        slimask           ! landmask: sea/land/ice=0/1/2

    real(kind_phys),               intent(in)           :: &
        ep1,            & ! (rv/rd) - 1
        con_rd,         & ! rd
        con_cp            ! cp
    real(kind_phys), dimension(:), intent(in)           :: &
        zorl,           & ! Surface roughness-length (cm)
        tsfc_wat,       & ! Surface skin temperature over water (K)
        tsfc_lnd,       & ! Surface skin temperature over land (K)
        tsfc_ice,       & ! Surface skin temperature over ice (K)
        qsfc_wat,       & ! Surface specific humidity over water (kg kg-1)
        qsfc_lnd,       & ! Surface specific humidity over land (kg kg-1)
        qsfc_ice,       & ! Surface specific humidity over ice(kg kg-1)
        hfx_wat,        & ! Kinematic surface upward sensible heat flux over water (K m s-1)
        hfx_lnd,        & ! Kinematic surface upward sensible heat flux over land (K m s-1)
        hfx_ice,        & ! Kinematic surface upward sensible heat flux over ice (K m s-1)
        qfx_wat,        & ! Kinematic surface upward latent heat flux over water (kg kg-1 m s-1)
        qfx_lnd,        & ! Kinematic surface upward latent heat flux over land (kg kg-1 m s-1)
        qfx_ice,        & ! Kinematic surface upward latent heat flux over ice (kg kg-1 m s-1)
        ps                ! Surface pressure (Pa)
    real(kind_phys), dimension(:), intent(in), optional :: &
        wetness           ! Normalized soil wetness (frac)  
    real(kind_phys), dimension(:,:), intent(in)         :: &
        t                 ! Model layer mean temperature (K)
    real(kind_phys), dimension(:,:,:), intent(in)       :: &
        q                 ! Model layer mean tracer concentration (kg kg-1)

    ! Outputs
    real(kind=kind_phys),           intent(out) :: shalwater_depth, svp1, svp2, svp3, svpt0, p1000mb
    real(kind_phys), dimension(:), intent(out) :: &
        tsfc,           & ! Surface skin temperature (K)
        qsfc,           & ! Surface specific humidity (kg kg-1)
        znt,            & ! Surface roughness-length (m)
        mavail,         & ! Surface moisture availability (frac; [0-1])
        pbl_regime,     & ! PBL regime categories
        fm,             & ! OM parameter for momentum
        fh,             & ! OM parameter for heat
        cpm,            & ! Heat capacity at constant pressure for moist air (J kg-1 K-1)
        xland,          & ! Sea/land/ice mask for MMM physics
        hfx,            & ! Surface kinematic upward sensible heat flux (W m-2)
        qfx,            & ! Surface kinematic upward latent heat flux (W m-2)
        qgh,            & ! Saturation water vapor mixing ratio at lowest model level (kg kg-1)
        gz1oz0,         & ! log(z/z0), where z0 is roughness length
        water_depth       ! Water depth for shallow water roughness scheme

    character(len=*),               intent(out) :: &
        errmsg            ! CCPP error message
    integer,                        intent(out) :: &
        its,            & ! Begin index of horizontal grid
        ite,            & ! End index of horizontal grid
        errflg            ! CCPP error code

    ! Locals
    integer         :: iCol
    real(kind_phys) :: tvcon, rho  

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! Horizontal dimension
    its = 1
    ite = nCol

    ! Consistency checks
    if (.not. do_mmm_sfclayrev) then
        write(errmsg,fmt='(*(a))') 'Logic error: do_mmm_sfclayrev = .false.'
        errflg = 1
        return
    end if

    ! Pass in flag for computing or prescribing surface fluxes
    if (do_mmm_sfclayrev) then
        scm_force_flux = .false.
    endif

    ! Constants from module_sf_mynn.F90
    shalwater_depth=1.0
    svp1 = 0.6112
    svp2 = 17.67
    svp3 = 29.65
    svpt0 = 273.15
    p1000mb = 100000.

    ! Convert roughness length (cm -> m)
    do iCol = 1,nCol
        znt(iCol) = zorl(iCol) * 0.01
    end do
    ! Calculate surface available moisture (fraction [0-1])
    mavail(:) = 1.0     ! WL2023 - should mavail be separated for dry, icy and wet?
    if (present(wetness)) then
        mavail(:) = wetness(:)
    end if

    ! Compute land/sea mask convention from (0-sea/1-land/2-ice) ---> (1-land/2-sea)
    do iCol=1,nCol
        if(slimask(iCol).eq.0)then
            xland(iCol) = 2
        else
            xland(iCol) = 1
        end if
    enddo

    do iCol = 1,nCol
       tvcon = (1. + ep1 * q(iCol,1,1))
       rho   = ps(iCol) / (con_rd * t(iCol,1) * tvcon)
        if (dry(iCol)) then
            hfx(iCol)  = rho * con_cp * hfx_lnd(iCol)
            qfx(iCol)  = rho * qfx_lnd(iCol)
            tsfc(iCol) = tsfc_lnd(iCol)
            qsfc(iCol) = qsfc_lnd(iCol)
        end if

        if (wet(iCol)) then
            hfx(iCol)  = rho * con_cp * hfx_wat(iCol)
            qfx(iCol)  = rho * qfx_wat(iCol)
            tsfc(iCol) = tsfc_wat(iCol)
            qsfc(iCol) = qsfc_wat(iCol)
        end if

        if (icy(iCol)) then
            hfx(iCol)  = rho * con_cp * hfx_ice(iCol)
            qfx(iCol)  = rho * qfx_ice(iCol)
            tsfc(iCol) = tsfc_ice(iCol)
            qsfc(iCol) = qsfc_ice(iCol)
        end if
    end do


  end subroutine gfs_mmm_sfclayrev_pre_run

end module gfs_mmm_sfclayrev_pre
