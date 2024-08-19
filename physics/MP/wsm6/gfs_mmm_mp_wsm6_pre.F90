! ###########################################################################################
!
! ###########################################################################################
module gfs_mmm_mp_wsm6_pre
  use machine, only: kind_phys

  implicit none

  public gfs_mmm_mp_wsm6_pre_run

contains
!> \section arg_table_gfs_mmm_mp_wsm6_pre_run
!! \htmlinclude gfs_mmm_mp_wsm6_pre_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine gfs_mmm_mp_wsm6_pre_run(nCol, nLev, phii, con_1ograv, dz, con_eps, prsl, con_rd, tgrs, spechum, rho, &
                                     re_qc_bg, re_qi_bg, re_qs_bg, re_qc_max, re_qi_max, re_qs_max, & 
                                     convert_dry_rho, qv, qc, qr, qi, qs, qg, &
                                     rainmp_mm, snowmp_mm, graupelmp_mm, &
                                     rain_nonphy_mm, snow_nonphy_mm, graupel_nonphy_mm, errmsg, errflg)

    ! Input variables
    integer,                         intent(in) :: nCol
    integer,                         intent(in) :: nLev
    real(kind_phys),                 intent(in) :: con_eps, con_rd, con_1ograv
    real(kind_phys), dimension(:,:), intent(in) :: &
        tgrs,               &  ! air temperature (K)
        prsl,               &  ! air pressure (Pa) 
        phii                   ! geopotential at model-interfaces (m2 s-2)

    logical,                         intent(in   ) :: convert_dry_rho
    real(kind_phys), dimension(:,:), intent(inout) :: spechum, qv, qc, qr, qi, qs, qg

    ! Output variables
    real(kind_phys),                 intent(out) :: re_qc_bg, re_qi_bg, re_qs_bg
    real(kind_phys),                 intent(out) :: re_qc_max, re_qi_max, re_qs_max
    real(kind_phys), dimension(:),   intent(out) :: rainmp_mm, snowmp_mm, graupelmp_mm, rain_nonphy_mm, snow_nonphy_mm, graupel_nonphy_mm
    real(kind_phys), dimension(:,:), intent(out) :: &
        rho,                &  ! air density (kg/m3)
        dz                     ! layer thickness (m)

    character(len=*),                intent(out) :: &
        errmsg                 ! CCPP error message
    integer,                         intent(out) :: &
        errflg                 ! CCPP error code

    ! Local variables
    integer                               :: i, k
    real(kind_phys),dimension(nCol, nLev+1) :: zi

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! Height (m) at interface from geopotential (m2/s2)
    do k = 1,nLev+1
       do i = 1,nCol
          zi(i,k) = phii(i,k)*con_1ograv
       enddo
    enddo

    ! Compute layer thickness (m)
    do k = 1,nLev
       do i = 1,nCol
          dz(i,k) = abs(zi(i,k)-zi(i,k+1))
       enddo
    enddo

    !> - Density of moist air in kg m-3
    do k=1,nLev
        do i=1,nCol
            ! convert specific humidity to water vapor mixing ratio
            qv(i,k) = spechum(i,k)/(1.0_kind_phys-spechum(i,k))
            ! calculate air density
            rho(i,k)= con_eps * prsl(i,k) / (con_rd * tgrs(i,k) * (qv(i,k) + con_eps))
        end do
    end do

    !> - WSM6 uses "dry" hydrometeor fields, e.g., water vapor mixing ratio (the mass of water vapor per unit mass of dry air)
    !> - CCPP only passes "wet" hydrometeor fields, e.g., specific humidity (the mass of water vapor per unit mass of air incl. vapor)
    !> - Need to do the following conversion
    if (convert_dry_rho) then
        qc = qc/(1.0_kind_phys-spechum)
        qr = qr/(1.0_kind_phys-spechum)
        qi = qi/(1.0_kind_phys-spechum)
        qs = qs/(1.0_kind_phys-spechum)
        qg = qg/(1.0_kind_phys-spechum)
    end if

    ! initialize background effective radius in meter (values from WRF)
    re_qc_bg = 2.49E-6
    re_qi_bg = 4.99E-6
    re_qs_bg = 9.99E-6

    ! initialize maximum effective radius in um (LF2024: have not checked if they are updated in WRF-v4.0)
    re_qc_max = 50.
    re_qi_max = 125.
    re_qs_max = 999.0

    ! intialize precipitation in mm
    rainmp_mm = 0.0
    snowmp_mm  = 0.0
    graupelmp_mm = 0.0
    rain_nonphy_mm = 0.0
    snow_nonphy_mm = 0.0
    graupel_nonphy_mm  = 0.0

  end subroutine gfs_mmm_mp_wsm6_pre_run

end module gfs_mmm_mp_wsm6_pre
