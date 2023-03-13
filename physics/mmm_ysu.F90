! ###########################################################################################
!
! This is a ccpp-compliant wrapper to call the YSU PBL scheme implemented @ NCAR MMM.
!
! ###########################################################################################
module mmm_ysu
  use machine, only: kind_phys
  use bl_ysu,  only: bl_ysu_run
  implicit none

  ! Module constants
  real(kind_phys) :: rovcp, rovg, conw

  public mmm_ysu_init, mmm_ysu_run

contains

!> \section arg_table_mmm_ysu_init
!! \htmlinclude mmm_ysu_init.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine mmm_ysu_init(con_rd, con_cp, con_g, errmsg, errflg)
    ! Inputs
    real(kind_phys),intent(in) :: con_cp, con_g, con_rd

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! Compute module constants.
    ! DJS2023: These may be defined by the host, need to check, if so, pass into _run phase()
    rovcp = con_rd/con_cp
    rovg  = con_rd/con_g
    conw  = 1._kind_phys/con_g

  end subroutine mmm_ysu_init

!> \section arg_table_mmm_ysu_run
!! \htmlinclude mmm_ysu_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine mmm_ysu_run(nCol, nLay, ntrac, slimask, u, v, t, q, ntcw, ntiw, p, pi, phii,   &
       prslk, ps, zorl, psim, psih, heat, evap, sfc_tau, wspd, u10, v10, con_cp, con_g,     &
       con_rd, ep1, ep2, br, karman, dtp, xlv, con_rv, xmu, lwh, swh, do_ysu_cldliq,        &
       do_ysu_cldice, ysu_add_bep, ysu_topdown_pblmix, ysu_timesplit, uo_sfc, vo_sfc, ctopo,&
       ctopo2, frcurb, a_u, a_v, a_t, a_q, a_e, b_u, b_v, b_t, b_q, b_e, dlu, dlg, sfk, vlk,&
       dudt_pbl, dvdt_pbl, dtdt_pbl, dqvdt_pbl, dqcdt_pbl, dqidt_pbl, dqtdt_pbl, exch_hx,   &
       exch_mx, hpbl, kpbl1d, wstar, delta, utnp, vtnp, ttnp,   &
       qtnp, u1, v1, t1, q1, qc1, qi1, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: do_ysu_cldliq, do_ysu_cldice, ysu_topdown_pblmix, ysu_timesplit
    integer, intent(in) :: nCol, nLay, ntrac, ntcw, ntiw
    integer, intent(in),dimension(:) :: slimask
    real(kind_phys),intent(in) :: con_cp, con_g, con_rd, ep1, ep2, karman, con_rv, xlv, dtp
    real(kind_phys),intent(in),dimension(:) :: ps, zorl, psim, psih, heat, evap, sfc_tau,   &
         wspd, br, xmu
    real(kind_phys),intent(in),dimension(:,:) :: u, v, t, p, pi, prslk, swh, lwh, phii
    real(kind_phys),intent(in),dimension(:,:,:) :: q
    real(kind_phys),intent(in),dimension(:), optional :: uo_sfc, vo_sfc, ctopo, ctopo2

    ! From Building Environment Parameterization (BEP) urbran-canopy model (optional)
    logical, intent(in) :: ysu_add_bep
    real(kind_phys),intent(in),dimension(:),   optional :: frcurb
    real(kind_phys),intent(in),dimension(:,:), optional :: a_u, a_v, a_t, a_q, b_u, b_v,    &
         b_t, b_q, sfk, vlk
    real(kind_phys),intent(in),dimension(:,:), optional :: a_e, b_e, dlu, dlg !*NOTE* These are not currently used in bl_ysu_run.

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    integer,          intent(out), dimension(:) :: kpbl1d
    real(kind_phys),  intent(out), dimension(:) :: hpbl, wstar, delta, exch_hx, exch_mx
    real(kind_phys),  intent(out), dimension(:,:) :: dudt_pbl, dvdt_pbl, dtdt_pbl,        &
         dqvdt_pbl, dqcdt_pbl, dqidt_pbl
    real(kind_phys),  intent(out), dimension(:,:,:) :: dqtdt_pbl
    
    ! In/Out
    real(kind_phys),  intent(inout), dimension(:)     :: u10, v10 !*NOTE* These are changed if topographical corrections are provided.
    real(kind_phys),  intent(inout), dimension(:,:)   :: utnp, vtnp, ttnp
    real(kind_phys),  intent(inout), dimension(:,:,:) :: qtnp
    real(kind_phys),  intent(inout), dimension(:,:)   :: u1, v1, t1, q1, qc1, qi1

    ! Locals
    integer :: iCol, iLay, ysu_topdown_pblmix_int !*NOTE* This will go away if/when bl_ysu_run accepts this switch directly as a logical.
    real(kind_phys) :: tvcon
    real(kind_phys),dimension(nCol) :: xland, hfx, qfx, rho, ust, znt
    real(kind_phys),dimension(nCol, nLay) :: rthraten, dz, exch_hx2D, exch_mx2D
    real(kind_phys),dimension(nCol, nLay+1) :: zi
    integer :: i,k

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! #######################################################################################
    !
    ! GFS MMM-YSU-PBL pre (compute inputs for YSU scheme...)
    !
    ! #######################################################################################

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
          zi(iCol,iLay) = phii(iCol,iLay)*conw
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
       rho(iCol) = ps(iCol)/(con_rd*t(iCol,1)*tvcon)
       hfx(iCol) = heat(iCol)*rho(iCol)*con_cp
       qfx(iCol) = evap(iCol)*rho(iCol)
    enddo

    ! Surface roughness length (cm) -> (meters)
    znt(:) = 0.01*zorl(:)

    ! Surface friction velocity (from surface wind-stress)
    ust(:) = sqrt(sfc_tau(:))

    ! YSU scheme needs 1/0 flag, we have .true./.false.
    ysu_topdown_pblmix_int = 0
    if (ysu_topdown_pblmix) ysu_topdown_pblmix_int = 1

    ! #######################################################################################
    !
    ! Call MMM-YSU-PBL scheme
    !
    ! #######################################################################################

    call bl_ysu_run(u, v, t, q(:,:,1), q(:,:,ntcw), q(:,:,ntiw), ntrac, q, p, pi, prslk,    &
         do_ysu_cldliq, do_ysu_cldice, dudt_pbl, dvdt_pbl, dtdt_pbl, dqvdt_pbl, dqcdt_pbl,  &
         dqidt_pbl, dqtdt_pbl, con_cp, con_g, rovcp, con_rd, rovg, ep1, ep2, karman, xlv,   &
         con_rv, dz, ps, znt, ust, hpbl, psim, psih, xland, hfx, qfx, wspd, br, dtp, kpbl1d,&
         exch_hx2D, exch_mx2D, wstar, delta, u10, v10, uo_sfc, vo_sfc, rthraten,            &
         ysu_topdown_pblmix_int, ctopo, ctopo2, a_u, a_v, a_t, a_q, a_e, b_u, b_v, b_t, b_q,&
         b_e, sfk, vlk, dlu, dlg, frcurb, ysu_add_bep, 1, nCol, nLay, nLay+1, errmsg, errflg)

    ! Only need exchange coefficients at the surface.
    exch_hx(:) = exch_hx2D(:,1)
    exch_mx(:) = exch_mx2D(:,1)

    ! #######################################################################################
    !
    ! GFS MMM-YSU-PBL post (couple YSU scheme...) 
    !
    ! #######################################################################################

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

  end subroutine mmm_ysu_run

end module mmm_ysu
