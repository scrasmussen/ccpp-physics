!>  \file cires_ugwp.F90
!! This file contains the Unified Gravity Wave Physics (UGWP) scheme by Valery Yudin (University of Colorado, CIRES)

!> This module contains the UGWP v0 scheme by Valery Yudin (University of Colorado, CIRES)
!!
!! See Valery Yudin's presentation at 2017 NGGPS PI meeting:
!! Gravity waves (GWs): Mesoscale GWs transport momentum, energy (heat) , and create eddy mixing in the whole atmosphere domain; Breaking and dissipating GWs deposit: (a) momentum; (b) heat (energy); and create (c) turbulent mixing of momentum, heat, and tracers
!! To properly incorporate GW effects (a-c) unresolved by DYCOREs we need GW physics
!! "Unified": a) all GW effects due to both dissipation/breaking; b) identical GW solvers for all GW sources; c) ability to replace solvers.
!! Unified Formalism:
!! - GW Sources: Stochastic and physics based mechanisms for GW-excitations in the lower atmosphere, calibrated by the high-res analyses/forecasts, and observations (3 types of GW sources: orography, convection, fronts/jets).
!! - GW Propagation: Unified solver for "propagation, dissipation and breaking" excited from all type of GW sources.
!! - GW Effects: Unified representation of GW impacts on the "resolved" flow for all sources (energy-balanced schemes for momentum, heat and mixing).
!! https://www.weather.gov/media/sti/nggps/Presentations%202017/02%20NGGPS_VYUDIN_2017_.pdf
module cires_ugwp

    use machine, only: kind_phys

    use cires_ugwpv0_module, only: knob_ugwp_version, cires_ugwpv0_mod_init, cires_ugwpv0_mod_finalize
    use ugwp_driver_v0
    use gwdps, only: gwdps_run
    use cires_ugwp_triggers

    implicit none

    private

    public cires_ugwp_init, cires_ugwp_run, cires_ugwp_finalize

    logical :: is_initialized = .False.

contains

! ------------------------------------------------------------------------
! CCPP entry points for CIRES Unified Gravity Wave Physics (UGWP) scheme v0
! ------------------------------------------------------------------------
!> The subroutine initializes the CIRES UGWP V0.
!> \section arg_table_cires_ugwp_init Argument Table
!! \htmlinclude cires_ugwp_init.html
!!
    subroutine cires_ugwp_init (me, master, nlunit, input_nml_file, logunit, &
                fn_nml2, lonr, latr, levs, ak, bk, dtp, cdmbgwd, cgwf,       &
                pa_rf_in, tau_rf_in, con_p0, gwd_opt,do_ugwp, errmsg, errflg)

!----  initialization of cires_ugwp
    implicit none

    integer,              intent (in) :: me
    integer,              intent (in) :: master
    integer,              intent (in) :: nlunit
    character(len=*),     intent (in) :: input_nml_file(:)
    integer,              intent (in) :: logunit
    integer,              intent (in) :: lonr
    integer,              intent (in) :: levs
    integer,              intent (in) :: latr
    real(kind=kind_phys), intent (in) :: ak(:), bk(:)
    real(kind=kind_phys), intent (in) :: dtp
    real(kind=kind_phys), intent (in) :: cdmbgwd(:), cgwf(:) ! "scaling" controls for "old" GFS-GW schemes
    real(kind=kind_phys), intent (in) :: pa_rf_in, tau_rf_in
    real(kind=kind_phys), intent (in) :: con_p0
    integer,              intent(in)  :: gwd_opt
    logical,              intent (in) :: do_ugwp

    character(len=*), intent (in) :: fn_nml2
    !character(len=*), parameter   :: fn_nml='input.nml'

    integer :: ios
    logical :: exists
    real    :: dxsg
    integer :: k

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (is_initialized) return
    
    ! Consistency checks
    if (gwd_opt/=1) then
      write(errmsg,'(*(a))') "Logic error: namelist choice of gravity wave &
      & drag is different from cires_ugwp scheme"
      errflg = 1
      return
    end if  

    if (do_ugwp .or. cdmbgwd(3) > 0.0) then
      call cires_ugwpv0_mod_init (me, master, nlunit, input_nml_file, logunit, &
                                fn_nml2, lonr, latr, levs, ak, bk, con_p0, dtp, &
                                cdmbgwd(1:2), cgwf, pa_rf_in, tau_rf_in)
    else
      write(errmsg,'(*(a))') "Logic error: cires_ugwp_init called but do_ugwp is false and cdmbgwd(3) <= 0"
      errflg = 1
      return
    end if

    if (.not.knob_ugwp_version==0) then
      write(errmsg,'(*(a))') 'Logic error: CCPP only supports version zero of UGWP'
      errflg = 1
      return
    end if

    is_initialized = .true.

    end subroutine cires_ugwp_init


! -----------------------------------------------------------------------
! finalize of cires_ugwp   (_finalize)
! -----------------------------------------------------------------------

!> The subroutine finalizes the CIRES UGWP
#if 0
!> \section arg_table_cires_ugwp_finalize Argument Table
!! \htmlinclude cires_ugwp_finalize.html
!!
#endif
    subroutine cires_ugwp_finalize(errmsg, errflg)

    implicit none
!
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not.is_initialized) return

    call cires_ugwpv0_mod_finalize()

    is_initialized = .false.

    end subroutine cires_ugwp_finalize


! -----------------------------------------------------------------------
!    originally from ugwp_driver_v0.f
!    driver of cires_ugwp   (_driver)
! -----------------------------------------------------------------------
!   driver is called after pbl & before chem-parameterizations
! -----------------------------------------------------------------------
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
! -----------------------------------------------------------------------
!>@brief These subroutines and modules execute the CIRES UGWP Version 0.
!> \section gen_cires_ugwp CIRES UGWP V0 Scheme General Algorithm
!! The physics of Non-Orographic Gravity Waves (NGWs) in the UGWP framework 
!!(Yudin et al. 2018 \cite yudin_et_al_2018) is represented by four GW-solvers, introduced in 
!!Lindzen (1981) \cite lindzen_1981, Hines (1997) \cite hines_1997, Alexander 
!!and Dunkerton (1999) \cite alexander_and_dunkerton_1999, and Scinocca (2003) \cite scinocca_2003. 
!!A major modification of these GW solvers was introduced with the addition of the 
!!background dissipation of temperature and winds to the saturation criteria for wave breaking. 
!!This feature is important in the mesosphere and thermosphere for WAM applications and it 
!!considers appropriate scale-dependent dissipation of waves near the model top lid providing 
!!the momentum and energy conservation in the vertical column physics (Shaw and 
!!Shepherd (2009) \cite shaw_and_shepherd_2009). In the UGWP-v0 scheme, a modification of the 
!!Scinocca (2003) \cite scinocca_2003 algorithm for NGWs with non-hydrostatic and rotational 
!!effects for GW propagations and background dissipation is contained in the subroutine 
!!fv3_ugwp_solv2_v0. Future development plans for the UGWP scheme include additional GW-solvers 
!!to be implemented along with physics-based triggering of waves and stochastic approaches 
!!for selection of GW modes characterized by horizontal phase velocities, azimuthal 
!!directions and magnitude of the vertical momentum flux (VMF).
!!
!! In UGWP-v0, the specification for the VMF function is adopted from the 
!! GEOS-5 global atmosphere model of GMAO NASA/GSFC, as described in 
!! Molod et al. (2015) \cite molod_et_al_2015 and employed in the MERRRA-2 
!! reanalysis (Gelaro et al., 2017 \cite gelaro_et_al_2017). The Fortran 
!! subroutine slat_geos5_tamp_v0() describes the latitudinal shape of 
!! VMF-function as displayed in Figure 3 of Molod et al. (2015) 
!! \cite molod_et_al_2015. It shows that the enhanced values of 
!! VMF in the equatorial region gives opportunity to simulate the 
!! QBO-like oscillations in the equatorial zonal winds and lead to more 
!! realistic simulations of the equatorial dynamics in GEOS-5 operational 
!! and MERRA-2 reanalysis products. For the first vertically extended 
!! version of FV3GFS in the stratosphere and mesosphere, this simplified 
!! function of VMF allows us to tune the model climate and to evaluate 
!! multi-year simulations of FV3GFS with the MERRA-2 and ERA-5 reanalysis 
!! products, along with temperature, ozone, and water vapor observations 
!! of current satellite missions. After delivery of the UGWP-code, the 
!! EMC group developed and tested approach to modulate the zonal mean 
!! NGW forcing by 3D-distributions of the total precipitation as a proxy 
!! for the excitation of NGWs by convection and the vertically-integrated  
!! (surface - tropopause) Turbulent Kinetic Energy (TKE). The verification 
!! scores with updated NGW forcing, as reported elsewhere by EMC researchers, 
!! display noticeable improvements in the forecast scores produced by 
!! FV3GFS configuration extended into the mesosphere.
!!
!> \section arg_table_cires_ugwp_run Argument Table
!! \htmlinclude cires_ugwp_run.html
!!
! \section det_cires_ugwp CIRES UGWP V0 Scheme Detailed Algorithm
     subroutine cires_ugwp_run(do_ugwp, me,  master, im,  levs, ntrac, dtp, kdt, lonr, &
         oro, oro_uf, hprime, nmtvr, oc, theta, sigma, gamma, elvmax, clx, oa4,        &
         do_tofd, ldiag_ugwp, cdmbgwd, xlat, xlat_d, sinlat, coslat,                   &
         area, ugrs, vgrs, tgrs, qgrs, prsi, prsl, prslk, phii, phil,                  &
         del, kpbl, dusfcg, dvsfcg, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,                &
         tau_tofd, tau_mtb, tau_ogw, tau_ngw, zmtb, zlwb, zogw,                        &
         dusfc_ms, dvsfc_ms, dusfc_bl, dvsfc_bl,                                       &
         dudt_ogw, dtauy2d_ms, dtaux2d_bl, dtauy2d_bl,                                 &
         dudt_mtb, dudt_tms, du3dt_mtb, du3dt_ogw, du3dt_tms,                          &
         dudt, dvdt, dtdt, rdxzb, con_g, con_pi, con_cp, con_rd, con_rv, con_fvirt,    &
         con_omega, rain, ntke, q_tke, dqdt_tke, lprnt, ipr,                           &
         dtend, dtidx, index_of_x_wind, index_of_y_wind, index_of_temperature,         &
         index_of_process_orographic_gwd, index_of_process_nonorographic_gwd,          &
         ldiag3d, lssav, flag_for_gwd_generic_tend, errmsg, errflg)

    implicit none

    ! interface variables
    integer,                 intent(in) :: me, master, im, levs, ntrac, kdt, lonr, nmtvr
    integer,                 intent(in), dimension(:)       :: kpbl
    real(kind=kind_phys),    intent(in), dimension(:)       :: oro, oro_uf, hprime, oc, theta, sigma, gamma
    logical,                 intent(in)                      :: flag_for_gwd_generic_tend
    ! elvmax is intent(in) for CIRES UGWP, but intent(inout) for GFS GWDPS
    real(kind=kind_phys),    intent(inout), dimension(:)    :: elvmax
    real(kind=kind_phys),    intent(in), dimension(:, :)    :: clx, oa4
    real(kind=kind_phys),    intent(in), dimension(:)       :: xlat, xlat_d, sinlat, coslat, area
    real(kind=kind_phys),    intent(in), dimension(:, :) :: del, ugrs, vgrs, tgrs, prsl, prslk, phil
    real(kind=kind_phys),    intent(in), dimension(:, :) :: prsi, phii
    real(kind=kind_phys),    intent(in), dimension(:,:,:):: qgrs
    real(kind=kind_phys),    intent(in) :: dtp, cdmbgwd(:)
    logical,                 intent(in) :: do_ugwp, do_tofd, ldiag_ugwp

    real(kind=kind_phys),    intent(out), dimension(:)      :: dusfcg, dvsfcg
    real(kind=kind_phys),    intent(out), dimension(:)      :: zmtb, zlwb, zogw, rdxzb
    real(kind=kind_phys),    intent(out), dimension(:)      :: tau_mtb, tau_ogw, tau_tofd, tau_ngw
    real(kind=kind_phys),    intent(out), dimension(:, :):: gw_dudt, gw_dvdt, gw_dtdt, gw_kdis
    real(kind=kind_phys),    intent(out), dimension(:, :):: dudt_mtb, dudt_tms
    real(kind=kind_phys),    intent(out), dimension(:, :), optional :: dudt_ogw
    real(kind=kind_phys),    intent(out), dimension(:), optional :: dusfc_ms, dvsfc_ms, dusfc_bl, dvsfc_bl
    real(kind=kind_phys),    intent(out), dimension(:, :), optional :: dtauy2d_ms
    real(kind=kind_phys),    intent(out), dimension(:, :), optional :: dtaux2d_bl, dtauy2d_bl

    ! dtend is only allocated if ldiag=.true.
    real(kind=kind_phys), optional, intent(inout) :: dtend(:,:,:)
    integer, intent(in)                                      :: dtidx(:,:), &
         index_of_x_wind, index_of_y_wind, index_of_temperature,         &
         index_of_process_orographic_gwd, index_of_process_nonorographic_gwd

    logical,                 intent(in)                         :: ldiag3d, lssav

    ! These arrays only allocated if ldiag_ugwp = .true.
    real(kind=kind_phys),    intent(inout), dimension(:,:), optional :: du3dt_mtb, du3dt_ogw, du3dt_tms

    real(kind=kind_phys),    intent(inout), dimension(:, :):: dudt, dvdt, dtdt

    real(kind=kind_phys),    intent(in) :: con_g, con_pi, con_cp, con_rd, con_rv, con_fvirt, con_omega

    real(kind=kind_phys),    intent(in), dimension(:) :: rain

    integer,                 intent(in) :: ntke
    real(kind=kind_phys),    intent(in), dimension(:,:) :: q_tke, dqdt_tke

    logical, intent(in) :: lprnt
    integer, intent(in) :: ipr

    character(len=*),        intent(out) :: errmsg
    integer,                 intent(out) :: errflg

    ! local variables
    integer :: i, k, idtend
    real(kind=kind_phys), dimension(im)       :: sgh30
    real(kind=kind_phys), dimension(im, levs) :: Pdvdt, Pdudt
    real(kind=kind_phys), dimension(im, levs) :: Pdtdt, Pkdis
    real(kind=kind_phys), dimension(im, levs) :: ed_dudt, ed_dvdt, ed_dtdt
    ! from ugwp_driver_v0.f -> cires_ugwp_initialize.F90 -> module ugwp_wmsdis_init
    real(kind=kind_phys), parameter :: tamp_mpa=30.e-3
    ! switches that activate impact of OGWs and NGWs (WL* how to deal with them? *WL)
    real(kind=kind_phys), parameter :: pogw=1., pngw=1., pked=1.

    real(kind=kind_phys), dimension(:,:), allocatable :: tke
    real(kind=kind_phys), dimension(:),   allocatable :: turb_fac, tem
    real(kind=kind_phys) :: rfac, tx1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! 1) ORO stationary GWs
    !    ------------------
    ! wrap everything in a do_ugwp 'if test' in order not to break the namelist functionality
    if (do_ugwp) then                       ! calling revised old GFS gravity wave drag

      ! topo paras
      ! w/ orographic effects
      if(nmtvr == 14)then
        ! calculate sgh30 for TOFD
        sgh30 = abs(oro - oro_uf)
       ! w/o orographic effects
      else
        sgh30   = 0.
      endif

      zlwb(:)   = 0.

      call GWDPS_V0(im, levs, lonr, do_tofd, Pdvdt, Pdudt, Pdtdt, Pkdis,          &
           ugrs, vgrs, tgrs, qgrs(:,:,1), kpbl, prsi,del,prsl, prslk, phii, phil, &
           dtp, kdt, sgh30, hprime, oc, oa4, clx, theta, sigma, gamma, elvmax,    &
           dusfcg, dvsfcg, xlat_d, sinlat, coslat, area, cdmbgwd(1:2),            &
           me, master, rdxzb, con_g, con_omega, zmtb, zogw, tau_mtb, tau_ogw,     &
           tau_tofd, dudt_mtb, dudt_ogw, dudt_tms)

     else                                    ! calling old GFS gravity wave drag as is

       do k=1,levs
         do i=1,im
           Pdvdt(i,k) = 0.0
           Pdudt(i,k) = 0.0
           Pdtdt(i,k) = 0.0
           Pkdis(i,k) = 0.0
         enddo
       enddo

       if (cdmbgwd(1) > 0.0 .or. cdmbgwd(2) > 0.0) then
         call gwdps_run(im, levs, Pdvdt, Pdudt, Pdtdt,                  &
                    ugrs, vgrs, tgrs, qgrs(:,:,1),                      &
                    kpbl, prsi, del, prsl, prslk, phii, phil, dtp, kdt, &
                    hprime, oc, oa4, clx, theta, sigma, gamma,          &
                    elvmax, dusfcg, dvsfcg, dudt_ogw, dtauy2d_ms,       &
                    dtaux2d_bl, dtauy2d_bl, dusfc_ms, dvsfc_ms,         &
                    dusfc_bl, dvsfc_bl,                                 &
                    con_g,  con_cp, con_rd, con_rv, lonr,               &
                    nmtvr, cdmbgwd, me, lprnt, ipr, rdxzb, ldiag_ugwp,  &
                    errmsg, errflg)
         if (errflg/=0) return
       endif

       tau_mtb   = 0.0  ; tau_ogw   = 0.0 ;  tau_tofd = 0.0
       if (ldiag_ugwp) then
         du3dt_mtb = 0.0  ; du3dt_ogw = 0.0 ;  du3dt_tms= 0.0
       endif

     endif ! do_ugwp


    if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
      idtend = dtidx(index_of_x_wind,index_of_process_orographic_gwd)
      if(idtend>=1) then
         dtend(:,:,idtend) = dtend(:,:,idtend) + Pdudt*dtp
      endif
      idtend = dtidx(index_of_y_wind,index_of_process_orographic_gwd)
      if(idtend>=1) then
         dtend(:,:,idtend) = dtend(:,:,idtend) + Pdvdt*dtp
      endif
      idtend = dtidx(index_of_temperature,index_of_process_orographic_gwd)
      if(idtend>=1) then
         dtend(:,:,idtend) = dtend(:,:,idtend) + Pdtdt*dtp
      endif
    endif
    

    if (cdmbgwd(3) > 0.0) then

      ! 2) non-stationary GW-scheme with GMAO/MERRA GW-forcing
      call slat_geos5_tamp_v0(im, tamp_mpa, xlat_d, tau_ngw)

      if (abs(1.0-cdmbgwd(3)) > 1.0e-6) then
        if (cdmbgwd(4) > 0.0) then
          allocate(turb_fac(im))
          do i=1,im
            turb_fac(i) = 0.0
          enddo
          if (ntke > 0) then
            allocate(tke(im,levs))
            allocate(tem(im))
            tke(:,:) = q_tke(:,:) + dqdt_tke(:,:) * dtp
            tem(:)   = 0.0
            do k=1,(levs+levs)/3
              do i=1,im
                turb_fac(i) = turb_fac(i) + del(i,k) * tke(i,k)
                tem(i)      = tem(i)      + del(i,k)
              enddo
            enddo
            do i=1,im
              turb_fac(i) = turb_fac(i) / tem(i)
            enddo
            deallocate(tke)
            deallocate(tem)
          endif
          rfac = 86400000 / dtp
          do i=1,im
            tx1 = cdmbgwd(4)*min(10.0, max(turb_fac(i),rain(i)*rfac))
            tau_ngw(i) = tau_ngw(i) * max(0.1, min(5.0, tx1))
          enddo
          deallocate(turb_fac)
        endif
        do i=1,im
          tau_ngw(i) = tau_ngw(i) * cdmbgwd(3)
        enddo
      endif

      call fv3_ugwp_solv2_v0(im, levs, dtp, tgrs, ugrs, vgrs,qgrs(:,:,1), &
                             prsl, prsi, phil, xlat_d, sinlat, coslat,    &
                             gw_dudt, gw_dvdt, gw_dtdt, gw_kdis, tau_ngw, &
                             me, master, kdt)

      do k=1,levs
        do i=1,im
          gw_dtdt(i,k) = pngw*gw_dtdt(i,k) + pogw*Pdtdt(i,k)
          gw_dudt(i,k) = pngw*gw_dudt(i,k) + pogw*Pdudt(i,k)
          gw_dvdt(i,k) = pngw*gw_dvdt(i,k) + pogw*Pdvdt(i,k)
          gw_kdis(i,k) = pngw*gw_kdis(i,k) + pogw*Pkdis(i,k)
          ! accumulation of tendencies for CCPP to replicate EMC-physics updates (!! removed in latest code commit to VLAB)
          !dudt(i,k) = dudt(i,k) + gw_dudt(i,k)
          !dvdt(i,k) = dvdt(i,k) + gw_dvdt(i,k)
          !dtdt(i,k) = dtdt(i,k) + gw_dtdt(i,k)
        enddo
      enddo

    else

      do k=1,levs
        do i=1,im
          gw_dtdt(i,k) = Pdtdt(i,k)
          gw_dudt(i,k) = Pdudt(i,k)
          gw_dvdt(i,k) = Pdvdt(i,k)
          gw_kdis(i,k) = Pkdis(i,k)
        enddo
      enddo

    endif

    if (pogw == 0.0) then
      tau_mtb  = 0. ; tau_ogw  = 0. ; tau_tofd = 0.
      dudt_mtb = 0. ; dudt_ogw = 0. ; dudt_tms = 0.
    endif

    if(ldiag3d .and. lssav .and. .not. flag_for_gwd_generic_tend) then
      idtend = dtidx(index_of_x_wind,index_of_process_nonorographic_gwd)
      if(idtend>=1) then
         dtend(:,:,idtend) = dtend(:,:,idtend) + (gw_dudt - Pdudt)*dtp
      endif
      idtend = dtidx(index_of_y_wind,index_of_process_nonorographic_gwd)
      if(idtend>=1) then
         dtend(:,:,idtend) = dtend(:,:,idtend) + (gw_dvdt - Pdvdt)*dtp
      endif
      idtend = dtidx(index_of_temperature,index_of_process_nonorographic_gwd)
      if(idtend>=1) then
         dtend(:,:,idtend) = dtend(:,:,idtend) + (gw_dtdt - Pdtdt)*dtp
      endif
    endif

    end subroutine cires_ugwp_run
end module cires_ugwp
