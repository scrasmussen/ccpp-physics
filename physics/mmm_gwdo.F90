! ###########################################################################################
!
! This is a ccpp-compliant wrapper to call the orographic gravity wave drag scheme implemented
! by NCAR MMM.
!
! ###########################################################################################
module mmm_gwdo
  use machine, only: kind_phys
  use bl_gwdo,  only: bl_gwdo_run
  implicit none

  public mmm_gwdo_init, mmm_gwdo_run

contains

!> \section arg_table_mmm_gwdo_init
!! \htmlinclude mmm_gwdo_init.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine mmm_gwdo_init(do_mmm_ogwd, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: do_mmm_ogwd

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    if (do_mmm_ogwd) then
       print*,'Using NCAR MMM OGWD scheme'
    endif

  end subroutine mmm_gwdo_init

!> \section arg_table_mmm_gwdo_run
!! \htmlinclude mmm_gwdo_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine mmm_gwdo_run(nCol, nLay, p, pi, prslk, phii, u, v, t, qv, dtp, con_1ograv,     &
       con_rd, con_cp, con_g, con_rv, con_fvirt, con_pi, dx, var, oc1, oa4, ol4,            &
       dusfc_og, dvsfc_og, dtaux2d_og, dtauy2d_og, utnp, vtnp, u1, v1, mmm_ogwd_timesplit,  &
       do_mmm_ogwd, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: mmm_ogwd_timesplit, do_mmm_ogwd
    integer, intent(in) :: nCol, nLay
    real(kind_phys), intent(in) :: dtp, con_1ograv, con_rd, con_cp, con_g, con_rv, con_fvirt,&
         con_pi
    real(kind_phys),intent(in),dimension(:)   :: dx, var, oc1
    real(kind_phys),intent(in),dimension(:,:) :: u, v, t, qv, p, pi, prslk, phii, oa4, ol4

    ! Input/Output
    real(kind_phys),  intent(inout), dimension(:,:) :: utnp, vtnp, u1, v1

    ! Outputs
    real(kind_phys),  intent(out), dimension(:)   :: dusfc_og, dvsfc_og
    real(kind_phys),  intent(out), dimension(:,:) :: dtaux2d_og, dtauy2d_og
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Local
    integer :: iCol, iLay
    real(kind_phys), dimension(nCol,nLay+1) :: zi
    real(kind_phys), dimension(nCol) :: sina, cosa
    real(kind_phys), dimension(nCol,nLay) :: dudt_gwd, dvdt_gwd

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    if (.not. do_mmm_ogwd) return

    ! #######################################################################################
    !
    ! GFS MMM-OGWD pre (compute inputs for MMM's OGWD scheme...)
    !
    ! #######################################################################################

    ! Height (m) at interface from geopotential (m2/s2)
    do iLay = 1,nLay+1
       do iCol = 1,nCol
          zi(iCol,iLay) = phii(iCol,iLay)*con_1ograv
       enddo
    enddo

    ! #######################################################################################
    !
    ! call bl_gwdo_run()
    !
    ! #######################################################################################
    !
    dudt_gwd(:,:) = 0.
    dvdt_gwd(:,:) = 0.

    call bl_gwdo_run(sina, cosa, dudt_gwd, dvdt_gwd, dtaux2d_og, dtauy2d_og, dusfc_og,      &
         dvsfc_og, u, v, t, qv, pi, p, prslk, zi, var, oc1, oa4(:,1), oa4(:,2), oa4(:,3),   &
         oa4(:,4), ol4(:,1), ol4(:,2), ol4(:,3), ol4(:,4), con_g, con_cp, con_rd, con_rv,   &
         con_fvirt, con_pi, dx, dtp, 1, nCol, nLay, nLay+1, .true., errmsg, errflg)

    ! #######################################################################################
    !
    ! GFS MMM OGWD post (couple OGWD scheme...) 
    !
    ! #######################################################################################    

    ! Procces-split scheme (default), accumulate tendencies...
    if (.not. mmm_ogwd_timesplit) then
       do iLay = 1,nLay
          do iCol = 1,nCol
             utnp(iCol,iLay) = utnp(iCol,iLay) + dudt_gwd(iCol,iLay)
             vtnp(iCol,iLay) = vtnp(iCol,iLay) + dvdt_gwd(iCol,iLay)
          enddo
       enddo
    ! Time-splt scheme, advance internal physics state...
    else
       do iLay = 1,nLay
          do iCol = 1,nCol
             u1(iCol,iLay)  = u(iCol,iLay) + dtp*dudt_gwd(iCol,iLay)
             v1(iCol,iLay)  = v(iCol,iLay) + dtp*dvdt_gwd(iCol,iLay)
          enddo
       enddo
    endif

  end subroutine mmm_gwdo_run
  !
end module mmm_gwdo
