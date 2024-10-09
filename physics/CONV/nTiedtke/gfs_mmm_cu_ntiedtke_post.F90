! ###########################################################################################
!
! ###########################################################################################
module gfs_mmm_cu_ntiedtke_post
  use machine, only: kind_phys

  implicit none

  public gfs_mmm_cu_ntiedtke_post_run

contains
!> \section arg_table_gfs_mmm_cu_ntiedtke_post_run
!! \htmlinclude gfs_mmm_cu_ntiedtke_post_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine gfs_mmm_cu_ntiedtke_post_run(nCol, nLev, qv, temp, spechum, prevst, prevsq, qc, qi, &
                                          ugrs, vgrs, pres, presi, geoph, geophi, pomg, forcet, forceq, &
                                          raincd_mm, raincd, &
                                          errmsg, errflg)

    ! Input variables
    integer,                         intent(in   ) :: &
        nCol,   & ! Number of horizontal gridpoints
        nLev      ! Number of vertical levels
    real(kind_phys), dimension(:),   intent(in   ) :: raincd_mm
    real(kind_phys), dimension(:,:), intent(in   ), optional :: qv
    
    ! In/out variables
    real(kind_phys), dimension(:,:),   intent(inout) :: temp, spechum, qc, qi, ugrs, vgrs, pres, presi, pomg
    real(kind_phys), dimension(:,:),   intent(inout), optional:: forcet, forceq, geoph, geophi

    ! Output variables
    real(kind_phys), dimension(:),     intent(out  ) :: raincd
    real(kind_phys), dimension(:,:),   intent(out  ), optional :: prevst, prevsq

    character(len=*),                intent(out) :: &
        errmsg            ! CCPP error message
    integer,                         intent(out) :: &
        errflg            ! CCPP error code

    ! Local variables
    integer                                :: i, k, kk
    real(kind_phys), dimension(nCol, nLev) :: pt, pqv, pqc, pqi, pu, pv, prsl, zl, omega, tendt, tendq
    real(kind_phys), dimension(nCol,nLev+1):: prsli, zi

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! Convert water vapor mixing ratio back to specific humidity
    spechum = qv/(1.0_kind_phys+qv)

    ! Variables with vertical layer being reversed back
    do k=1,nLev
        kk = nLev-k+1
        do i=1,nCol
            pt(i,k)    = temp(i,kk)
            pqv(i,k)   = spechum(i,kk)
            pqc(i,k)   = qc(i,kk)
            pqi(i,k)   = qi(i,kk)
            pu(i,k)    = ugrs(i,kk)
            pv(i,k)    = vgrs(i,kk)
            prsl(i,k)  = pres(i,kk)
            prsli(i,k) = presi(i,kk+1)
            zl(i,k)    = geoph(i,kk)
            zi(i,k)    = geophi(i,kk+1)
            omega(i,k) = pomg(i,kk)
            tendt(i,k) = forcet(i,kk)
            tendq(i,k) = forceq(i,kk)
        enddo
    enddo
    prsli(:,nLev+1)=presi(:,1)
       zi(:,nLev+1)=geophi(:,1)

    !output
    temp    = pt
    spechum = pqv
    qc      = pqc
    qi      = pqi
    ugrs    = pu
    vgrs    = pv
    pres    = prsl
    presi   = prsli
    geoph   = zl
    geophi  = zi
    pomg    = omega
    forcet  = tendt
    forceq  = tendq

    ! To calculate tendencies in next iteration
    prevst = temp
    prevsq = spechum

    ! Convert unit of deep convection induced precipitation from mm to m
    raincd = raincd_mm * 1e-3

  end subroutine gfs_mmm_cu_ntiedtke_post_run

end module gfs_mmm_cu_ntiedtke_post
