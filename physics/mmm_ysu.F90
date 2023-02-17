module mmm_ysu
  use machine, only: kind_phys
  use bl_ysu,  only: bl_ysu_run
  implicit none

  public mmm_ysu_run

contains

!> \section arg_table_mmm_ysu_run
!! \htmlinclude mmm_ysu_run.html
!!
  ! #########################################################################################
  ! #########################################################################################
  subroutine mmm_ysu_run(nCol, nLay, nmix, slimask, u, v, t, qx, ntcw, ntiw, p2d, p2di, phii, pi2d, &
       ps, zorl, psim, psih, heat, evap, sfc_tau, wspd, u10, v10, cp, g, rd, ep1, ep2, br, karman, dtp, xlv, rv, &
       xmu, lwh, swh, do_ysu_cldliq, do_ysu_cldice, dudt_pbl, dvdt_pbl, dtdt_pbl, dqvdt_pbl, dqcdt_pbl, dqidt_pbl, dqtdt_pbl, exch_hx, exch_mx, hpbl, kpbl1d, wstar, delta, errmsg, errflg)
    ! Inputs
    logical, intent(in) :: do_ysu_cldliq, do_ysu_cldice
    integer, intent(in) :: nCol, nLay, nmix, ntcw, ntiw
    integer, intent(in),dimension(:) :: slimask
    real(kind_phys),intent(in) :: cp, g, rd, ep1, ep2, karman, rv, xlv, dtp
    real(kind_phys),intent(in),dimension(:) :: ps, zorl, psim, psih, heat, evap, sfc_tau, wspd, br, xmu
    real(kind_phys),intent(in),dimension(:,:) :: u, v, t, p2d, p2di, pi2d, swh, lwh, phii
    real(kind_phys),intent(in),dimension(:,:,:) :: qx
    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    integer,          intent(inout), dimension(:) :: kpbl1d
    
    real(kind_phys),  intent(inout), dimension(:) :: u10, v10, hpbl, wstar, delta, exch_hx, exch_mx
    real(kind_phys),  intent(inout), dimension(:,:) :: dudt_pbl, dvdt_pbl, dtdt_pbl, dqvdt_pbl, dqcdt_pbl, dqidt_pbl
    real(kind_phys),  intent(inout), dimension(:,:,:) :: dqtdt_pbl

    ! Locals
    integer :: iCol, iLay
    real :: tvcon, rovcp, rovg
    real(kind_phys),dimension(nCol) :: xland, hfx, qfx, rhox, govrth, ust, znt, uox, vox, ctopo, ctopo2, frcurb
    real(kind_phys),dimension(nCol, nLay) :: thx, thlix, thvx, rthraten, dz8w2d, a_u, a_v, a_t, a_q, a_e, b_u, b_v, b_t, b_q, b_e, dlu, dlg, sfk, vlk
    integer :: ysu_topdown_pblmix = 1
    logical :: flag_bep = .true.

    ! Compute constants
    rovcp = rd/cp
    rovg  = rd/g

    ! Compute land/sea mask for ysu
    do iCol=1,nCol
       if(slimask(iCol).eq.0) then
          xland(iCol) = 2
       else
          xland(iCol) = 1
       end if
    end do

    do iLay = 1,nLay
       do iCol = 1,nCol
          thx(iCol,iLay)   = t(iCol,iLay)/pi2d(iCol,iLay)
          thlix(iCol,iLay) = (t(iCol,iLay)-xlv*qx(iCol,iLay,ntcw)/cp-2.834E6*qx(iCol,iLay,ntiw)/cp)/pi2d(iCol,iLay)
       enddo
    enddo
    
    do iLay = 1,nLay
       do iCol = 1,nCol
          tvcon = (1.+ep1*qx(iCol,iLay,1))
          thvx(iCol,iLay) = thx(iCol,iLay)*tvcon
       enddo
    enddo

    ! Convert total radiative heating rate from temp/s -> theta/s
    do iLay = 1,nLay
       do iCol = 1,nCol
          rthraten(iCol,iLay) = (swh(iCol,iLay)*xmu(iCol)+lwh(iCol,iLay))/pi2d(iCol,iLay)
       enddo
    enddo

    do iCol = 1,nCol
       tvcon = (1.+ep1*qx(iCol,1,1))
       rhox(iCol) = ps(iCol)/(rd*t(iCol,1)*tvcon)
       govrth(iCol) = g/thx(iCol,1)
       hfx(iCol) = heat(iCol)*rhox(iCol)*cp
       qfx(iCol) = evap(iCol)*rhox(iCol)
       ust(iCol) = sqrt(sfc_tau(iCol))
       znt(iCol) = 0.01*zorl(iCol)
       uox(iCol) = 0.0
       vox(iCol) = 0.0
       ctopo(iCol) = 0.0
       ctopo2(iCol) = 0.0
       a_u(iCol,:) = 0.0
       a_v(iCol,:) = 0.0
       a_t(iCol,:) = 0.0
       a_q(iCol,:) = 0.0
       a_e(iCol,:) = 0.0
       b_u(iCol,:) = 0.0
       b_v(iCol,:) = 0.0
       b_t(iCol,:) = 0.0
       b_q(iCol,:) = 0.0
       b_e(iCol,:) = 0.0
       dlu(iCol,:) = 0.0
       dlg(iCol,:) = 0.0
       sfk(iCol,:) = 0.0
       vlk(iCol,:) = 0.0
       frcurb(iCol) = 0.0
    enddo

    ! Convert height from geopotential
    do iLay = 1,nLay
       do iCol = 1,nCol
          dz8w2d(iCol,iLay) = phii(iCol,iLay)/g
       enddo
    enddo

  ! #########################################################################################
    call bl_ysu_run(u, v, t, qx(:,:,1), qx(:,:,ntcw), qx(:,:,ntiw), nmix, qx, p2d, p2di, &
         pi2d, do_ysu_cldliq, do_ysu_cldice, dudt_pbl, dvdt_pbl, dtdt_pbl, dqvdt_pbl,       &
         dqcdt_pbl, dqidt_pbl, dqtdt_pbl, cp, g, rovcp, rd, rovg, ep1, ep2, karman, xlv, rv,&
         dz8w2d, psfcpa, znt, ust, hpbl, psim, psih, xland, hfx, qfx, wspd, br, dtp, &
         kpbl1d, exch_hx, exch_mx, wstar, delta, u10, v10, &
         uox, &! Optional
         vox, & ! Optional
         rthraten, &
         ysu_topdown_pblmix, &
         ctopo, & ! Optional 
         ctopo2, & ! Optional 
         a_u, & ! Optional
         a_v, & ! Optional 
         a_t, & ! Optional 
         a_q, & ! Optional 
         a_e, & ! Optional 
         b_u, & ! Optional 
         b_v, & ! Optional 
         b_t, & ! Optional 
         b_q, & ! Optional 
         b_e, & ! Optional 
         sfk, & ! Optional 
         vlk, & ! Optional 
         dlu, & ! Optional 
         dlg, & ! Optional 
         frcurb, & ! Optional 
         flag_bep, &
         1, nCol, nLay, nLay+1, errmsg, errflg)
    

  end subroutine mmm_ysu_run
end module mmm_ysu
