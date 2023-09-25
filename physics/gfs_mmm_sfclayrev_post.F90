! ###########################################################################################
!
! ###########################################################################################
module gfs_mmm_sfclayrev_post
  use machine, only: kind_phys

  implicit none

  public gfs_mmm_sfclayrev_post_run

contains

!> \section arg_table_gfs_mmm_sfclayrev_post_run
!! \htmlinclude gfs_mmm_sfclayrev_post_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine gfs_mmm_sfclayrev_post_run(nCol, wet, dry, icy, qsfc, qsfc_wat, qsfc_lnd, qsfc_ice,&
                                        znt, zorl, zorl_wat, zorl_lnd, zorl_ice,                &
                                        ust, ust_wat, ust_lnd, ust_ice,                         &
                                        cm, cm_wat, cm_lnd, cm_ice,                             &
                                        cka, cka_wat, cka_lnd, cka_ice,                         &
                                        br, br_wat, br_lnd, br_ice,                             &
                                        stress, stress_wat, stress_lnd, stress_ice,             &
                                        psim, psim_wat, psim_lnd, psim_ice,                     &
                                        psih, psih_wat, psih_lnd, psih_ice,                     &
                                        psim10, psim10_wat, psim10_lnd, psim10_ice,             &
                                        psih2, psih2_wat, psih2_lnd, psih2_ice,                 &
                                        ep1, t, q, ps, con_rd, con_cp, hfx, qfx,      &
                                        hflx, hflx_wat, hflx_lnd, hflx_ice,                     &
                                        qflx, qflx_wat, qflx_lnd, qflx_ice,                     &
                                        errmsg, errflg)

    ! Input
    integer,                       intent(in)           :: &
        nCol              ! Number of horizontal gridpoints

    logical,         dimension(:), intent(in)           :: &
        wet,            & ! Flag indicating presence of some ocean or lake surface area fraction
        dry,            & ! Flag indicating presence of some land surface area fraction
        icy               ! Flag indicating presence of some sea ice surface area fraction 

    real(kind_phys),               intent(in)           :: &
        ep1,            & ! (rv/rd) - 1
        con_rd,         & ! rd
        con_cp            ! cp
    real(kind_phys), dimension(:), intent(in)           :: &
        qsfc,           & ! Surface specific humidity (kg kg-1)
        znt,            & ! Surface roughness length in meter (m)
        ust,            & ! Surface frictional velocity (m s-1)
        cm,             & ! Momentum exchange coefficient at the lowest model level
        cka,            & ! Enthalpy exchange coefficient at the lowest model level
        br,             & ! Bulk Richardson number at the surface
        psim,           & ! Monin-Obukhov similarity function for momentum
        psih,           & ! Monin-Obukhov similarity function for heat
        psim10,         & ! Monin-Obukhov similarity function for momentum at 10m
        psih2,          & ! Monin-Obukhov similarity function for heat at 2m
        ps                ! Surface pressure (Pa)
    real(kind_phys), dimension(:,:), intent(in)         :: &
        t                 ! Model layer mean temperature (K)
    real(kind_phys), dimension(:,:,:), intent(in)       :: &
        q                 ! Model layer mean tracer concentration (kg kg-1)

    ! In/out
    real(kind_phys), dimension(:), intent(inout) :: &
        qsfc_wat,       & ! Surface specific humidity over water (kg kg-1)
        qsfc_lnd,       & ! Surface specific humidity over land (kg kg-1)
        qsfc_ice,       & ! Surface specific humidity over ice (kg kg-1)
        hfx,            & ! Kinematic surface upward sensible heat flux (W m-2)
        qfx,            & ! Kinematic surface upward latent heat flux (W m-2)
        hflx,           & ! Kinematic surface upward sensible heat flux (K m s-1)
        hflx_wat,       & ! Kinematic surface upward sensible heat flux over water (K m s-1)
        hflx_lnd,       & ! Kinematic surface upward sensible heat flux over land (K m s-1)
        hflx_ice,       & ! Kinematic surface upward sensible heat flux over ice (K m s-1)
        qflx,           & ! Kinematic surface upward latent heat flux (kg kg-1 m s-1)
        qflx_wat,       & ! Kinematic surface upward latent heat flux over water (kg kg-1 m s-1)
        qflx_lnd,       & ! Kinematic surface upward latent heat flux over land (kg kg-1 m s-1) 
        qflx_ice,       & ! Kinematic surface upward latent heat flux over ice (kg kg-1 m s-1)
        zorl,           & ! Surface roughness-length in centimeter (cm)
        zorl_wat,       & ! Surface roughness-length in centimeter over water (cm)
        zorl_lnd,       & ! Surface roughness-length in centimeter over land (cm)
        zorl_ice,       & ! Surface roughness-length in centimeter over ice (cm)
        ust_wat,        & ! Surface friction velocity over water (m s-1)
        ust_lnd,        & ! Surface friction velocity over land (m s-1)
        ust_ice,        & ! Surface friction velocity over ice (m s-1)
        cm_wat,         & ! Surface exchange coeff for momentum over water
        cm_lnd,         & ! Surface exchange coeff for momentum over land
        cm_ice,         & ! Surface exchange coeff for momentum over ice
        cka_wat,        & ! Surface exchange coeff heat & moisture over water
        cka_lnd,        & ! Surface exchange coeff heat & moisture over land
        cka_ice,        & ! Surface exchange coeff heat & moisture over ice
        br_wat,         & ! Bulk Richardson number at the surface over water
        br_lnd,         & ! Bulk Richardson number at the surface over land
        br_ice,         & ! Bulk Richardson number at the surface over ice
        stress,         & ! Surface wind stress (m2 s-2)
        stress_wat,     & ! Surface wind stress over water (m2 s-2)
        stress_lnd,     & ! Surface wind stress over land (m2 s-2)
        stress_ice,     & ! Surface wind stress over ice (m2 s-2)
        psim_wat,       & ! Monin-Obukhov similarity function for momentum over water
        psim_lnd,       & ! Monin-Obukhov similarity function for momentum over land
        psim_ice,       & ! Monin-Obukhov similarity function for momentum over ice
        psih_wat,       & ! Monin-Obukhov similarity function for heat over water
        psih_lnd,       & ! Monin-Obukhov similarity function for heat over land
        psih_ice,       & ! Monin-Obukhov similarity function for heat over ice
        psim10_wat,     & ! Monin-Obukhov similarity function for momentum at 10m over water
        psim10_lnd,     & ! Monin-Obukhov similarity function for momentum at 10m over land
        psim10_ice,     & ! Monin-Obukhov similarity function for momentum at 10m over ice
        psih2_wat,      & ! Monin-Obukhov similarity function for heat at 10m over water
        psih2_lnd,      & ! Monin-Obukhov similarity function for heat at 10m over land
        psih2_ice         ! Monin-Obukhov similarity function for heat at 10m over ice

    ! Output
    character(len=*), intent(out) :: &
         errmsg         !
    integer,          intent(out) :: &
         errflg         !

    ! Locals
    integer         :: iCol
    real(kind_phys) :: tvcon, rho

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! Convert roughness length (m -> cm)
    zorl(:) = znt(:)*100.

    ! Surface wind stress (m2 s-2)
    stress(:) = ust(:)*ust(:)
    
    ! Calculate surface kinematic upward sensible and latent heat fluxes (W m-2)
    do iCol = 1,nCol
       tvcon = (1. + ep1 * q(iCol,1,1))
       rho   = ps(iCol) / (con_rd * t(iCol,1) * tvcon)
       hflx(iCol) = hfx(iCol)/rho/con_cp
       qflx(iCol) = qfx(iCol)/rho

        if (dry(iCol)) then
            qsfc_lnd(iCol) = qsfc(iCol)
            hflx_lnd(iCol) = hflx(iCol)
            qflx_lnd(iCol) = qflx(iCol)
            zorl_lnd(iCol) = zorl(iCol)
            ust_lnd(iCol)  = ust(iCol)
            cm_lnd(iCol)  = cm(iCol)
            cka_lnd(iCol)  = cka(iCol)
            br_lnd(iCol)  = br(iCol)
            stress_lnd(iCol) = stress(iCol)
            psim_lnd(iCol)  = psim(iCol)
            psih_lnd(iCol)  = psih(iCol)
            psim10_lnd(iCol)  = psim10(iCol)
            psih2_lnd(iCol)  = psih2(iCol)
        end if

        if (wet(iCol)) then
            qsfc_wat(iCol) = qsfc(iCol)
            hflx_wat(iCol) = hflx(iCol)
            qflx_wat(iCol) = qflx(iCol)
            zorl_wat(iCol) = zorl(iCol)
            ust_wat(iCol)  = ust(iCol)
            cm_wat(iCol)  = cm(iCol)
            cka_wat(iCol)  = cka(iCol)
            br_wat(iCol)  = br(iCol)
            stress_wat(iCol) = stress(iCol)
            psim_wat(iCol)  = psim(iCol)
            psih_wat(iCol)  = psih(iCol)
            psim10_wat(iCol)  = psim10(iCol)
            psih2_wat(iCol)  = psih2(iCol)
        end if

        if (icy(iCol)) then
            qsfc_ice(iCol) = qsfc(iCol)
            hflx_ice(iCol) = hflx(iCol)
            qflx_ice(iCol) = qflx(iCol)
            zorl_ice(iCol) = zorl(iCol)
            ust_ice(iCol)  = ust(iCol)
            cm_ice(iCol)  = cm(iCol)
            cka_ice(iCol)  = cka(iCol)
            br_ice(iCol)  = br(iCol)
            stress_ice(iCol) = stress(iCol)
            psim_ice(iCol)  = psim(iCol)
            psih_ice(iCol)  = psih(iCol)
            psim10_ice(iCol)  = psim10(iCol)
            psih2_ice(iCol)  = psih2(iCol)
        end if

    enddo

  end subroutine gfs_mmm_sfclayrev_post_run

end module gfs_mmm_sfclayrev_post
