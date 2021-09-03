module zm_conv_all_post
  
  use shr_kind_mod,    only: r8 => shr_kind_r8
  
  implicit none
  
  contains
!! \section arg_table_zm_conv_all_post_run
!! \htmlinclude zm_conv_all_post_run.html
!!
  subroutine zm_conv_all_post_run(ncol, pcols, pver, pcnst, cam_physpkg, cam_physpkg_cam3, ixcldice, ixcldliq, &
    tend_qv_deep_conv, tend_s_deep_conv, &
    tend_qv_deep_conv_evap, tend_s_deep_conv_evap, tend_s_momtran, tend_u_momtran, tend_v_momtran, tend_q_convtran, &
    tend_u_all, tend_v_all, tend_s_all, tend_q_all, errmsg, errflg)
    
    use iap_ptend_sum, only : physics_ptend_sum
    
    integer,           intent(in   )              :: ncol, pcols, pver, pcnst
    integer, intent(in) :: ixcldice, ixcldliq
    character(len=16), intent(in) :: cam_physpkg, cam_physpkg_cam3
    real(kind=r8), intent(in), dimension(:,:) :: tend_qv_deep_conv, tend_s_deep_conv
    real(kind=r8), intent(in), dimension(:,:) :: tend_qv_deep_conv_evap, tend_s_deep_conv_evap
    real(kind=r8), intent(in), dimension(:,:) :: tend_s_momtran, tend_u_momtran, tend_v_momtran
    real(kind=r8), intent(in), dimension(:,:,:) :: tend_q_convtran
    
    real(kind=r8), intent(inout), dimension(:,:) :: tend_u_all, tend_v_all, tend_s_all
    real(kind=r8), intent(inout), dimension(:,:,:) :: tend_q_all
    
    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg
    
    integer :: i, top, bot
    logical :: lu, lv, ls
    logical, dimension(pcnst) :: lq
    logical :: lu_all, lv_all, ls_all
    logical, dimension(pcnst) :: lq_all
    real(kind=r8), dimension(pcols,pver) :: temp_tend_u, temp_tend_v, temp_tend_s
    real(kind=r8), dimension(pcols,pver,pcnst) :: temp_tend_q
    real(kind=r8), dimension(pcols) :: temp_taux_srf, temp_taux_top, temp_tauy_srf, temp_tauy_top, temp_hflx_srf, temp_hflx_top
    real(kind=r8), dimension(pcols,pcnst) :: temp_cflx_srf, temp_cflx_top
    
    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    top = 1
    bot = pver
    
    temp_taux_srf = 0._r8
    temp_taux_top = 0._r8
    temp_tauy_srf = 0._r8
    temp_tauy_top = 0._r8
    temp_hflx_srf = 0._r8
    temp_hflx_top = 0._r8
    temp_cflx_srf = 0._r8
    temp_cflx_top = 0._r8
    
    ! sum tendencies due to main zm_convr call
    lu = .false.
    lv = .false.
    ls = .true.
    lq(1) = .true.
    lq(2:) = .false.
    temp_tend_u(:,:) = 0._r8
    temp_tend_v(:,:) = 0._r8
    temp_tend_q(:,:,:) = 0._r8
    temp_tend_q(:,:,1) = tend_qv_deep_conv(:,:)
    
    call physics_ptend_sum(ncol, pcnst, top, bot, lu, lv, ls, lq, temp_tend_u, temp_tend_v, &
      tend_s_deep_conv, temp_tend_q, temp_taux_srf, temp_taux_top, temp_tauy_srf, temp_tauy_top, &
      temp_hflx_srf, temp_hflx_top, temp_cflx_srf, temp_cflx_top, &
      lu_all, lv_all, ls_all, lq_all, temp_tend_u, temp_tend_v, tend_s_all, tend_q_all, &
      temp_taux_srf, temp_taux_top, temp_tauy_srf, temp_tauy_top, temp_hflx_srf, temp_hflx_top, &
      temp_cflx_srf, temp_cflx_top)
      
    ! sum tendencies due to zm_conv_evap call
    lu = .false.
    lv = .false.
    ls = .true.
    lq(1) = .true.
    lq(2:) = .false.
    temp_tend_q(:,:,:) = 0._r8
    temp_tend_q(:,:,1) = tend_qv_deep_conv_evap(:,:)
    
    call physics_ptend_sum(ncol, pcnst, top, bot, lu, lv, ls, lq, temp_tend_u, temp_tend_v, &
      tend_s_deep_conv_evap, temp_tend_q, temp_taux_srf, temp_taux_top, temp_tauy_srf, temp_tauy_top, &
      temp_hflx_srf, temp_hflx_top, temp_cflx_srf, temp_cflx_top, &
      lu_all, lv_all, ls_all, lq_all, temp_tend_u, temp_tend_v, tend_s_all, tend_q_all, &
      temp_taux_srf, temp_taux_top, temp_tauy_srf, temp_tauy_top, temp_hflx_srf, temp_hflx_top, &
      temp_cflx_srf, temp_cflx_top)
    
    ! sum tendencies due to zm_conv_momtran call  
    if (cam_physpkg /= cam_physpkg_cam3) then
       lu = .true.
       lv = .true.
       ls = .true.
       lq(:) = .false.
       temp_tend_q(:,:,:) = 0._r8
       
       call physics_ptend_sum(ncol, pcnst, top, bot, lu, lv, ls, lq, tend_u_momtran, tend_v_momtran, &
         tend_s_momtran, temp_tend_q, temp_taux_srf, temp_taux_top, temp_tauy_srf, temp_tauy_top, &
         temp_hflx_srf, temp_hflx_top, temp_cflx_srf, temp_cflx_top, &
         lu_all, lv_all, ls_all, lq_all, tend_u_all, tend_v_all, tend_s_all, temp_tend_q, &
         temp_taux_srf, temp_taux_top, temp_tauy_srf, temp_tauy_top, temp_hflx_srf, temp_hflx_top, &
         temp_cflx_srf, temp_cflx_top)
    end if
    
    ! sum tendencies due to zm_conv_convtran call
    lu = .false.
    lv = .false.
    ls = .false.
    lq(:) = .false.
    lq(ixcldice) = .true.
    lq(ixcldliq) = .true.
    
    call physics_ptend_sum(ncol, pcnst, top, bot, lu, lv, ls, lq, temp_tend_u, temp_tend_v, &
      temp_tend_s, tend_q_convtran, temp_taux_srf, temp_taux_top, temp_tauy_srf, temp_tauy_top, &
      temp_hflx_srf, temp_hflx_top, temp_cflx_srf, temp_cflx_top, &
      lu_all, lv_all, ls_all, lq_all, temp_tend_u, temp_tend_v, temp_tend_s, tend_q_all, &
      temp_taux_srf, temp_taux_top, temp_tauy_srf, temp_tauy_top, temp_hflx_srf, temp_hflx_top, &
      temp_cflx_srf, temp_cflx_top)
    
  end subroutine zm_conv_all_post_run
end module zm_conv_all_post
