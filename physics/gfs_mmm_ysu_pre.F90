! ###########################################################################################
!
! ###########################################################################################
module gfs_mmm_ysu_pre
  use machine, only: kind_phys

  implicit none

  public gfs_mmm_ysu_pre_run

contains

!> \section arg_table_gfs_mmm_ysu_pre_run
!! \htmlinclude gfs_mmm_ysu_pre_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine gfs_mmm_ysu_pre_run(nCol, nLay, ntcw, ntiw, con_1ograv, ep1, con_rd, con_cp, slimask, ps,  &
       xmu, heat, evap, sfc_tau, t, phii, swh, lwh, q, xland, hfx, qfx, ust, rthraten, dz,  &
       znt, qvx, qcx, qix, errmsg, errflg)

    ! Inputs
    integer, intent(in) :: &
         nCol,       & ! Number of horizontal points
         nLay,       & ! Number of vertical grid points.
         ntcw,       & !
         ntiw
    real(kind_phys),intent(in) :: &
         con_1ograv, & ! Physical constant:
         ep1,        & ! Physical constant:
         con_rd,     & ! Physical constant:
         con_cp        ! Physical constant:
    integer, intent(in),dimension(:) :: &
         slimask       ! sea/land/ice mask for UFS applications
    real(kind_phys),intent(in),dimension(:) :: &
         ps,         & ! Surface pressure (Pa)
         xmu,        & ! zenith angle adjustment for shortwave (none)
         heat,       & ! Surface upward sensible heat flux (K m s-1)
         evap,       & ! Surface upward latent heat flux (kg kg-1 m s-1)
         sfc_tau,    & ! Surface wind stress (m2 s-2)
         znt
    real(kind_phys),intent(in),dimension(:,:) :: &
         t,          & ! Temperature @ model layer-centers (K)
         phii,       & ! Geopotential at model-interfaces (m2 s-2)
         swh,        & ! Total sky shortwave heating rate (K s-1)
         lwh           ! Total sky longwave heating rate (K s-1)
    real(kind_phys),intent(in),dimension(:,:,:) :: &
         q             ! tracer concentration (kg kg-1)

    ! Outputs
    real(kind_phys), intent(out), dimension(:) :: &
         xland         ! sea/land/ice mask for MMM physics
    real(kind_phys), intent(out), dimension(:) :: &
         hfx,        & ! Surface upward sensible heat flux (W m-2)
         qfx,        & ! Surface upward latent heat flux (W m-2)
         ust           !
    real(kind_phys), intent(out), dimension(:,:) :: &
         rthraten,   & ! Total radiative heating rate, in potential temp
         dz,         & !
         qvx,        & !
         qcx,        & !
         qix           !
    character(len=*), intent(out) :: &
         errmsg        ! CCPP error message
    integer,          intent(out) :: &
         errflg        ! CCPP error flag

    ! Locals
    integer :: iCol, iLay
    real(kind_phys) :: tvcon, rho
    real(kind_phys),dimension(nCol, nLay+1) :: zi

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! 
    qvx = q(:,:,1)
    qcx = q(:,:,ntcw)
    qix = q(:,:,ntiw)

    ! Compute land/sea mask convention for YSU from (0-sea/1-land/2-ice) ---> (1-land/2-sea)
    do iCol=1,nCol
       if(slimask(iCol).eq.0) then
          xland(iCol) = 2
       else
          xland(iCol) = 1
       end if
    end do

    ! Height (m) at interface from geopotential (m2/s2)
    do iLay = 1,nLay+1
       do iCol = 1,nCol
          zi(iCol,iLay) = phii(iCol,iLay)*con_1ograv
       enddo
    enddo

    ! Compute layer thickness (m)
    do iLay = 1,nLay
       do iCol = 1,nCol
          dz(iCol,iLay) = abs(zi(iCol,iLay)-zi(iCol,iLay+1))
       enddo
    enddo

    ! Total (longwave+shortwave) radiative heating rate (K/s)
    do iLay = 1,nLay
       do iCol = 1,nCol
          rthraten(iCol,iLay) = (swh(iCol,iLay)*xmu(iCol)+lwh(iCol,iLay))
       enddo
    enddo

    ! Kinematic surface fluxes
    do iCol = 1,nCol
       tvcon     = (1.+ep1*q(iCol,1,1))
       rho       = ps(iCol)/(con_rd*t(iCol,1)*tvcon)
       hfx(iCol) = heat(iCol)*rho*con_cp
       qfx(iCol) = evap(iCol)*rho
    enddo

    ! Surface friction velocity (from surface wind-stress)
    ust(:) = sqrt(sfc_tau(:))

  end subroutine gfs_mmm_ysu_pre_run

end module gfs_mmm_ysu_pre
