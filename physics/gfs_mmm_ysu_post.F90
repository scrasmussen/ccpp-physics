! ###########################################################################################
!
! ###########################################################################################
module gfs_mmm_ysu_post
  use machine, only: kind_phys

  implicit none

  public gfs_mmm_ysu_post_run

contains

!> \section arg_table_gfs_mmm_ysu_post_run
!! \htmlinclude gfs_mmm_ysu_post_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine gfs_mmm_ysu_post_run(nCol, nLay, ntcw, ntiw, ntqv, index_of_temperature,       &
       index_of_x_wind, index_of_y_wind, index_of_process_pbl, top_at_1, ysu_timesplit,     &
       flag_for_pbl_generic_tend, lssav, ldiag3d, dtp, dtidx, u, v, t, q, exch_hx, exch_mx, &
       prslk, dudt_pbl, dvdt_pbl, dtdt_pbl, dqvdt_pbl, dqcdt_pbl, dqidt_pbl, dqtdt_pbl,     &
       dtsfc, utnp, vtnp, ttnp, qtnp, u1, v1, t1, q1, qc1, qi1, dtend, sfc_exch_hx,         &
       sfc_exch_mx, errmsg, errflg)

    ! Inputs
    integer, intent(in) :: &
         nCol,                      & !
         nLay,                      & !
         ntcw,                      & !
         ntiw,                      & !
         ntqv,                      & !
         index_of_temperature,      & !
         index_of_x_wind,           & !
         index_of_y_wind,           & !
         index_of_process_pbl         !
    logical, intent(in) :: &
         top_at_1,                  & !
         ysu_timesplit,             & !
         flag_for_pbl_generic_tend, & !
         lssav,                     & !
         ldiag3d                      !
    real(kind_phys), intent(in) :: &
         dtp                          !
    integer,         intent(in), dimension(:,:) :: &
         dtidx                        !
    real(kind_phys), intent(in), dimension(:,:) :: &
         u,                         & !
         v,                         & !
         t,                         & !
         exch_hx,                   & !
         exch_mx,                   & !
         prslk                        !
    real(kind_phys), intent(in), dimension(:,:,:) :: &
         q                            !

    ! In/out
    real(kind_phys), intent(inout), dimension(:) :: &
         dtsfc        !
    real(kind_phys), intent(inout), dimension(:,:) :: &
         dudt_pbl,  & !
         dvdt_pbl,  & !
         dtdt_pbl,  & !
         dqvdt_pbl, & !
         dqcdt_pbl, & !
         dqidt_pbl, & !
         utnp,      & !
         vtnp,      & !
         ttnp,      & !
         u1,        & !
         v1,        & !
         t1,        & !
         q1,        & !
         qc1,       & !
         qi1          !
    real(kind_phys), intent(inout), dimension(:,:,:) :: &
         dqtdt_pbl, & !
         qtnp,      & !
         dtend        !

    ! Outputs
    real(kind_phys), intent(out), dimension(:) :: &
         sfc_exch_hx, & ! DJS2023: These are stored in GFS_diag_type, maybe should be in GFS_interstital_type?
         sfc_exch_mx    ! DJS2023: These are stored in GFS_diag_type, maybe should be in GFS_interstital_type?
    character(len=*), intent(out) :: &
         errmsg         !
    integer,          intent(out) :: &
         errflg         !

    ! Locals
    integer :: iCol, iLay, idtend

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! Save exchange coefficients at the surface.
    sfc_exch_hx(:) = exch_hx(:,1)
    sfc_exch_mx(:) = exch_mx(:,1)

    ! Convert tendencies from potential-temperature -> temperature
    do iCol = 1,nCol
       if (top_at_1) then
          dtsfc(iCol) = dtsfc(iCol) * prslk(iCol,nLay)
       else
          dtsfc(iCol) = dtsfc(iCol) * prslk(iCol,1)
       endif
       do iLay = 1,nLay
          dtdt_pbl(iCol,iLay) = dtdt_pbl(iCol,iLay)*prslk(iCol,iLay)
       enddo
    enddo

    ! Procces-split scheme (default), accumulate tendencies...
    if (.not. ysu_timesplit) then
       do iLay = 1,nLay
          do iCol = 1,nCol
             utnp(iCol,iLay)      = utnp(iCol,iLay)      + dudt_pbl(iCol,iLay)
             vtnp(iCol,iLay)      = vtnp(iCol,iLay)      + dvdt_pbl(iCol,iLay)
             ttnp(iCol,iLay)      = ttnp(iCol,iLay)      + dtdt_pbl(iCol,iLay)
             qtnp(iCol,iLay,1)    = qtnp(iCol,iLay,1)    + dqvdt_pbl(iCol,iLay)
             qtnp(iCol,iLay,ntcw) = qtnp(iCol,iLay,ntcw) + dqcdt_pbl(iCol,iLay)
             qtnp(iCol,iLay,ntiw) = qtnp(iCol,iLay,ntiw) + dqidt_pbl(iCol,iLay)
          enddo
       enddo
    ! Time-splt scheme, advance internal physics state...
    else
       do iLay = 1,nLay
          do iCol = 1,nCol
             u1(iCol,iLay)  = u(iCol,iLay)      + dtp*dudt_pbl(iCol,iLay)
             v1(iCol,iLay)  = v(iCol,iLay)      + dtp*dvdt_pbl(iCol,iLay)
             t1(iCol,iLay)  = t(iCol,iLay)      + dtp*dtdt_pbl(iCol,iLay)
             q1(iCol,iLay)  = q(iCol,iLay,1)    + dtp*dqvdt_pbl(iCol,iLay)
             qc1(iCol,iLay) = q(iCol,iLay,ntcw) + dtp*dqcdt_pbl(iCol,iLay)
             qi1(iCol,iLay) = q(iCol,iLay,ntiw) + dtp*dqidt_pbl(iCol,iLay)
          enddo
       enddo
    endif

    ! Save physics tendencies.
    if(lssav .and. ldiag3d .and. .not. flag_for_pbl_generic_tend) then
       ! Temperature
       idtend = dtidx(index_of_temperature,index_of_process_pbl)
       if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + dtp*dtdt_pbl
       endif
       ! Specific-humidity
       idtend = dtidx(ntqv+100,index_of_process_pbl)
       if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + dtp*dqvdt_pbl
       endif
       ! Zonal-wind
       idtend = dtidx(index_of_x_wind,index_of_process_pbl)
       if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + dtp*dudt_pbl
       endif
       ! Meridional-wind
       idtend = dtidx(index_of_y_wind,index_of_process_pbl)
       if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) + dtp*dvdt_pbl
       endif
       ! Cloud liquid
       ! Cloud ice
    endif

  end subroutine gfs_mmm_ysu_post_run

end module gfs_mmm_ysu_post
