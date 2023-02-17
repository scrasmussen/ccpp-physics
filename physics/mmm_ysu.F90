module mmm_ysu
  use machine, only: kind_phys
  use bl_ysu,  only: bl_ysu_run
  implicit none

  ! Module constants
  real(kind_phys) :: rovcp, rovg

  public mmm_ysu_init, mmm_ysu_run

contains

!> \section arg_table_mmm_ysu_init
!! \htmlinclude mmm_ysu_init.html
!!
  ! #########################################################################################
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

    ! Compute some constants.
    rovcp = con_rd/con_cp
    rovg  = con_rd/con_g

  end subroutine mmm_ysu_init

!> \section arg_table_mmm_ysu_run
!! \htmlinclude mmm_ysu_run.html
!!
  ! #########################################################################################
  ! #########################################################################################
  subroutine mmm_ysu_run(nCol, nLay, ntrac, slimask, u, v, t, q, ntcw, ntiw, p, pi, phii,   &
       prslk, ps, zorl, psim, psih, heat, evap, sfc_tau, wspd, u10, v10, con_cp, con_g,     &
       con_rd, ep1, ep2, br, karman, dtp, xlv, con_rv, xmu, lwh, swh, do_ysu_cldliq,        &
       do_ysu_cldice, ysu_add_bep, ysu_topdown_pblmix, uo_sfc, vo_sfc, ctopo, ctopo2,       &
       frcurb, a_u, a_v, a_t, a_q, a_e, b_u, b_v, b_t, b_q, b_e, dlu, dlg, sfk, vlk,        &
       dudt_pbl, dvdt_pbl, dtdt_pbl, dqvdt_pbl, dqcdt_pbl, dqidt_pbl, dqtdt_pbl, exch_hx,   &
       exch_mx, hpbl, kpbl1d, wstar, delta, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: do_ysu_cldliq, do_ysu_cldice, ysu_add_bep, ysu_topdown_pblmix
    integer, intent(in) :: nCol, nLay, ntrac, ntcw, ntiw
    integer, intent(in),dimension(:) :: slimask
    real(kind_phys),intent(in) :: con_cp, con_g, con_rd, ep1, ep2, karman, con_rv, xlv, dtp
    real(kind_phys),intent(in),dimension(:) :: ps, zorl, psim, psih, heat, evap, sfc_tau,   &
         wspd, br, xmu
    real(kind_phys),intent(in),dimension(:,:) :: u, v, t, p, pi, prslk, swh, lwh, phii
    real(kind_phys),intent(in),dimension(:,:,:) :: q
    real(kind_phys),intent(in),dimension(:), optional :: uo_sfc, vo_sfc, ctopo, ctopo2,     &
         frcurb
    real(kind_phys),intent(in),dimension(:,:), optional :: a_u, a_v, a_t, a_q, a_e, b_u,    &
         b_v, b_t, b_q, b_e, dlu, dlg, sfk, vlk ! DJS: These still need metadata
    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    integer,          intent(inout), dimension(:) :: kpbl1d
    real(kind_phys),  intent(inout), dimension(:) :: u10, v10, hpbl, wstar, delta, exch_hx, &
         exch_mx
    real(kind_phys),  intent(inout), dimension(:,:) :: dudt_pbl, dvdt_pbl, dtdt_pbl,        &
         dqvdt_pbl, dqcdt_pbl, dqidt_pbl
    real(kind_phys),  intent(inout), dimension(:,:,:) :: dqtdt_pbl

    ! Locals
    integer :: iCol, iLay, ysu_topdown_pblmix_int
    real(kind_phys) :: tvcon
    real(kind_phys),dimension(nCol) :: xland, hfx, qfx, rho, govrth, ust, znt
    real(kind_phys),dimension(nCol, nLay) :: theta, rthraten
    real(kind_phys),dimension(nCol, nLay+1) :: zi

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    ! #######################################################################################
    ! GFS MMM-YSU-PBL pre
    ! Compute inputs for YSU scheme...
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
          zi(iCol,iLay) = phii(iCol,iLay)/con_g
       enddo
    enddo

    ! Potential temperature (K)
    do iLay = 1,nLay
       do iCol = 1,nCol
          theta(iCol,iLay)   = t(iCol,iLay)/prslk(iCol,iLay)
       enddo
    enddo
    govrth(:) = con_g/theta(:,1)   ! gravity divided by theta at the surface
    
    ! Total (longwave+shortwave) radiative heating rate from temp/s -> theta/s
    do iLay = 1,nLay
       do iCol = 1,nCol
          rthraten(iCol,iLay) = (swh(iCol,iLay)*xmu(iCol)+lwh(iCol,iLay))/prslk(iCol,iLay)
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

    ! YSU scheme needs 0/1 flag, we have .true./.false.
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
         con_rv, zi, ps, znt, ust, hpbl, psim, psih, xland, hfx, qfx, wspd, br, dtp,        &
         kpbl1d, exch_hx, exch_mx, wstar, delta, u10, v10, uo_sfc, vo_sfc, rthraten,        &
         ysu_topdown_pblmix_int, ctopo, ctopo2, a_u, a_v, a_t, a_q, a_e, b_u, b_v, b_t, b_q,&
         b_e, sfk, vlk, dlu, dlg, frcurb, ysu_add_bep, 1, nCol, nLay, nLay+1, errmsg, errflg)

  end subroutine mmm_ysu_run
end module mmm_ysu
