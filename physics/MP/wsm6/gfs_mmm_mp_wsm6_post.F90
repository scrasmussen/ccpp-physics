! ###########################################################################################
!
! ###########################################################################################
module gfs_mmm_mp_wsm6_post
  use machine, only: kind_phys

  implicit none

  public gfs_mmm_mp_wsm6_post_run

contains
!> \section arg_table_gfs_mmm_mp_wsm6_post_run
!! \htmlinclude gfs_mmm_mp_wsm6_post_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine gfs_mmm_mp_wsm6_post_run(nCol, nLev, re_qc, re_qi, re_qs, re_qc_um, re_qi_um, re_qs_um, &
                                      convert_dry_rho, qv, spechum, qc, qr, qi, qs, qg, &
                                      rainmp_mm, snowmp_mm, graupelmp_mm, rainmp, snowmp, graupelmp, &
                                      rain_nonphy_mm, snow_nonphy_mm, graupel_nonphy_mm, &
                                      rain_nonphy, snow_nonphy, graupel_nonphy, prcpmp, errmsg, errflg)

    ! input variables
    integer,                         intent(in) :: nCol
    integer,                         intent(in) :: nLev
    real(kind_phys), dimension(:,:), intent(in) :: re_qc, re_qi, re_qs
    real(kind_phys), dimension(:,:), intent(in) :: qv
    real(kind_phys), dimension(:),   intent(in) :: rainmp_mm, snowmp_mm, graupelmp_mm, rain_nonphy_mm, snow_nonphy_mm, graupel_nonphy_mm
    logical,                         intent(in) :: convert_dry_rho

    ! output variables
    real(kind_phys), dimension(:,:), intent(out) :: re_qc_um, re_qi_um, re_qs_um
    real(kind_phys), dimension(:,:), intent(inout) :: spechum, qc, qr, qi, qs, qg
    real(kind_phys), dimension(:),   intent(out) :: rainmp, snowmp, graupelmp, rain_nonphy, snow_nonphy, graupel_nonphy, prcpmp

    character(len=*),                intent(out) :: &
        errmsg            ! CCPP error message
    integer,                         intent(out) :: &
        errflg            ! CCPP error code
    ! local variables
    integer                               :: i, k

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! Convert unit of meter to micrometer
    do k=1,nLev
        do i=1,nCol
            re_qc_um(i,k) = re_qc(i,k) * 1e6
            re_qi_um(i,k) = re_qi(i,k) * 1e6
            re_qs_um(i,k) = re_qs(i,k) * 1e6
        end do
    end do

    !> - Convert water vapor mixing ratio back to specific humidity
    spechum = qv/(1.0_kind_phys+qv)

    if (convert_dry_rho) then
        qc = qc/(1.0_kind_phys+qv)
        qr = qr/(1.0_kind_phys+qv)
        qi = qi/(1.0_kind_phys+qv)
        qs = qs/(1.0_kind_phys+qv)
        qg = qg/(1.0_kind_phys+qv)
    end if

    !> - Convert unit of precipitation from mm to m
    rainmp = rainmp_mm / 1e3
    snowmp = snowmp_mm / 1e3
    graupelmp = graupelmp_mm / 1e3

    rain_nonphy = rain_nonphy_mm / 1e3
    snow_nonphy = snow_nonphy_mm / 1e3
    graupel_nonphy = graupel_nonphy_mm / 1e3

    !> - Output total amount of precip (rain, snow, graupel) on physics timestep
    prcpmp = rainmp + snowmp + graupelmp

  end subroutine gfs_mmm_mp_wsm6_post_run

end module gfs_mmm_mp_wsm6_post
