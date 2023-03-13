! ###########################################################################################
!
! This is a ccpp-compliant wrapper to call the orographic gravity wave drag scheme implemented
! by NCAR MMM.
!
! ###########################################################################################
module mmm_sfclay
  use machine,      only: kind_phys
  use sf_sfclayrev, only: sf_sfclayrev_run, sf_sfclayrev_init
  implicit none

  public mmm_sfclay_init, mmm_sfclay_run

contains

!> \section arg_table_mmm_sfclay_init
!! \htmlinclude mmm_sfclay_init.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine mmm_sfclay_init(errmsg, errflg)

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    !
    call sf_sfclayrev_init(errmsg, errflg)

  end subroutine mmm_sfclay_init

!> \section arg_table_mmm_sfclay_run
!! \htmlinclude mmm_sfclay_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine mmm_sfclay_run(nCol, nLay, mmm_ogwd_sfcflx, dx, slimask, ps, tsfc, qsfc,  &
        lh, br1, chs2, cqs2, u, v, t, q, p, pi, phii, con_1ograv, con_cp, con_g,       &
       con_rocp, con_rd, con_hvap, ep1, ep2, karman, prsl1, hpbl, stress_wat, stress_lnd,   &
       stress_ice,rmol, zorl, mol, u10, v10, t2m, q2m, hflx, qflx, flhc, flqc, zol, wspd, psim, psih,   &
       isftcflx, iz0tlnd, ustm, ck, cka, cd, cda, errmsg, errflg)

    ! Inputs
    integer, intent(in) :: &
         nCol,         & ! Number of horizontal gridpoints
         nLay,         & ! Number of vertical layers
         mmm_ogwd_sfcflx ! If true, use supplied surface heat/moisture fluxes in scheme.
    real(kind_phys), intent(in) :: &
         con_1ograv,   & ! Physical constant
         con_cp,       & ! Physical constant
         con_g,        & ! Physical constant
         con_rocp,     & ! Physical constant
         con_rd,       & ! Physical constant
         con_hvap,     & ! Physical constant
         ep1,          & ! Physical constant
         ep2,          & ! Physical constant
         karman          ! Physical constant
    integer, dimension(:), intent(in) :: &
         slimask         ! landmask: sea/land/ice=0/1/2          (none)
    real(kind_phys), dimension(:), intent(in) :: &
         prsl1,        & ! Air pressure @ lowest model-layer     (Pa)
         hpbl,         & ! PBL height                            (m)
         stress_wat,   & ! Wind-stress @ surface over water      (m2 s-2)
         stress_lnd,   & ! Wind-stress @ surface over land       (m2 s-2)
         stress_ice,   & ! Wind-stress @ surface over ice        (m2 s-2)
         ps,           & ! Air pressure @ surface                (Pa)
         tsfc,         & ! Air temperature @ surface             (K)
         dx              ! Horizontal grid size                  (m)
    real(kind_phys), dimension(:,:), intent(in) :: &
         u,            & ! u-wind @ layer-center                 (m s-1)
         v,            & ! v-wind @ layer-center                 (m s-1)
         t,            & ! Air temperature @ layer-center        (K)
         q,            & ! Specific humidity @ layer-center      (kg kg-1)
         p,            & ! Air pressure @ layer-center           (Pa)
         pi,           & ! Air pressure @ layer-interface        (Pa)
         phii            ! Geopotential  @ layer-interface       (m2 s-2)

    ! Input/Output
    real(kind_phys), dimension(:), intent(inout) :: &
         qsfc,         & ! Surface specific-humidity             (kg kg-1)
         lh,           & ! Latent heat @ surface                 (W m-2)
         rmol,         & ! Reciprocal of Monin-Obukhov length    (m-1)
         zorl,         & ! Surface roughness-length              (cm)
         mol,          & ! T* (similarity theory)                (K)
         flhc,         & ! Exchange coefficient for heat         (W m-2 K-1)
         flqc,         & ! Exchange coefficient for moisture     (kg m-2 s-1)
         wspd,         & ! Wind speed at lowest model level      (m s-1) 
         br1,          & ! Bulk Richardson number @ surface      (1)
         chs2,         & ! heat exchange coefficient for LSM     (m s-1)
         cqs2,         & ! moisture exchange coefficient for LSM (m s-1)
         zol,          & ! z/L height over Monin-Obukhov length  (m)
         u10,          & ! U-wind @ 10m                          (m s-1)
         v10,          & ! V-wind @ 10m                          (m s-1)
         t2m,          & ! Air-temperature @ 2m                  (K)
         q2m,          & ! Specific-humidity @ 2m                (kg kg-1)
         hflx,         & ! Upward sensible heat flux @ surface   (K m s-1)
         qflx            ! Upward latent heat flux @ surface     (kg kg-1 m s-1)
    ! Outputs
    real(kind_phys), dimension(:), intent(out) :: &
         psim,         & ! MO similarity parameter for momentum  (none)
         psih            ! MO similarity parameter for heat      (none)
    character(len=*), intent(out) :: &
         errmsg          ! CCPP error message
    integer, intent(out) :: &
         errflg          ! CCPP error code

    ! Optional
    integer, intent(in), optional :: &
         isftcflx,     & ! Flag to control thermal roughness-length over water.
         iz0tlnd         ! Flag to control thermal roughness-length over land.
    real(kind_phys), intent(out),   dimension(:), optional :: &
         ck,           & ! Enthalpy exchange coeff @ 10m
         cka,          & ! Enthalpy exchange coeff @ lowest model level
         cd,           & ! Momentum exchange coeff @ 10m
         cda             ! Momentum exchange coeff @ lowest model level 
    real(kind_phys), intent(inout), dimension(:), optional :: &
         ustm            ! u* in similarity theory (m/s) w* added to wspd

    ! Locals
    integer :: iCol, iLay, isfc, isfflx, shalwater_z0
    real(kind_phys) :: tvcon, p1000mb=100000.
    real(kind_phys), dimension(nCol) :: ust, sfc_mavail, xland, rho, chs, cpm, pbl_regime, &
         fm, fh, gz1oz0, water_depth, znt, qgh
    real(kind_phys), dimension(nCol,nLay) :: dz
    real(kind_phys), dimension(nCol,nLay+1) :: zi
    real(kind_phys), parameter :: svp1    = 0.6112
    real(kind_phys), parameter :: svp2    = 17.67
    real(kind_phys), parameter :: svp3    = 29.65
    real(kind_phys), parameter :: svpt0   = 273.15

    real(kind_phys) :: eomeg=1, stbolt=1, shalwater_depth=1 ! NOT USED!!!

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! Logics (Move to typedefs)
    isfflx = 1
    shalwater_z0 = 0

    ! #######################################################################################
    !
    ! Prepare inputs for NCAR MMM surface-layer scheme
    !
    ! #######################################################################################

    ! Compute land/sea mask convention for scheme, from (0-sea/1-land/2-ice) ---> (1-land/2-sea)
    do iCol=1,nCol
       if(slimask(iCol).eq.0) then
          xland(iCol) = 2
       else
          xland(iCol) = 1
       end if
    end do

    ! Compute layer-height
    do iLay = 1,nLay+1
       do iCol = 1,nCol
          zi(iCol,iLay) = phii(iCol,iLay)*con_1ograv
       enddo
    enddo

    ! Compute layer-thickness
    do iLay = 1,nLay
       do iCol = 1,nCol
          dz(iCol,iLay) = abs(zi(iCol,iLay)-zi(iCol,iLay+1))
       enddo
    enddo

    ! Surface friction velocity (from surface wind-stress)
    do iCol = 1,nCol
       if (slimask(iCol) == 0) ust(icol) = sqrt(stress_wat(iCol))
       if (slimask(iCol) == 1) ust(iCol) = sqrt(stress_lnd(iCol))
       if (slimask(iCol) == 2) ust(iCol) = sqrt(stress_ice(iCol))
    enddo

    ! Convert roughness length (cm -> m)
    znt(:) = zorl(:)*0.01

    ! Compute surface moisture availability
    sfc_mavail(:) = 1

    isfc = 1

    ! #######################################################################################
    !
    ! Call NCAR MMM surface-layer scheme
    !
    ! #######################################################################################
    write(*,'(a12,12f10.6)') 'MMM_SFCLAY:',hflx, qflx, tsfc, u10, v10, t2m, t2m, q2m, flhc, flqc, qgh, qsfc
    call sf_sfclayrev_run(u(:,isfc), v(:,isfc), t(:,isfc), q(:,isfc), p(:,isfc), dz(:,isfc),&
         con_cp, con_g, con_rocp, con_rd, con_hvap, prsl1, chs, chs2, cqs2, cpm, hpbl, rmol,&
         znt, ust, sfc_mavail, zol, mol, pbl_regime, psim, psih, fm, fh, xland, hflx, qflx, &
         tsfc, u10, v10, t2m, t2m, q2m, flhc, flqc, qgh, qsfc, lh, gz1oz0, wspd, br1,       &
         isfflx, dx, svp1, svp2, svp3, svpt0, ep1, ep2, karman, eomeg, stbolt, p1000mb,     &
         shalwater_z0, water_depth, shalwater_depth, isftcflx, iz0tlnd, mmm_ogwd_sfcflx,    &
         ustm, ck, cka, cd, cda, 1, nCol, errmsg, errflg)
    write(*,'(a12,12f10.6)') 'MMM_SFCLAY:',hflx, qflx, tsfc, u10, v10, t2m, t2m, q2m, flhc, flqc, qgh, qsfc

    ! #######################################################################################
    !
    ! #######################################################################################

    ! Convert roughness length (m -> cm)
    zorl(:) = znt(:)*100.

  end subroutine mmm_sfclay_run
  !
end module mmm_sfclay
