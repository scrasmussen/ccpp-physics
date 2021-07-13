module iap_ptend_sum
  
  use shr_kind_mod, only: r8=>shr_kind_r8
  
  public
  
  implicit none
  
  contains
  
  subroutine physics_ptend_sum(ncol, pcnst, top_level, bot_level, lu_in, lv_in, ls_in, lq_in, &
    u_tend_in, v_tend_in, s_tend_in, q_tend_in, taux_srf_in, taux_top_in, tauy_srf_in, tauy_top_in, &
    hflx_srf_in, hflx_top_in, cflx_srf_in, cflx_top_in, &
    lu_sum, lv_sum, ls_sum, lq_sum, u_tend_sum, v_tend_sum, s_tend_sum, q_tend_sum, &
    taux_srf_sum, taux_top_sum, tauy_srf_sum, tauy_top_sum, hflx_sfc_sum, hflx_top_sum, cflx_srf_sum, cflx_top_sum)
!-----------------------------------------------------------------------
! Add ptend fields to ptend_sum for ptend logical flags = .true.
! Where ptend logical flags = .false, don't change ptend_sum
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
    integer, intent(in) :: ncol, pcnst, top_level, bot_level
    logical, intent(in) :: lu_in, lv_in, ls_in, lq_in(pcnst)
    real(kind=r8), intent(in) :: u_tend_in(:,:), v_tend_in(:,:), s_tend_in(:,:), q_tend_in(:,:,:)
    real(kind=r8), intent(in) :: taux_srf_in(:), taux_top_in(:), tauy_srf_in(:), tauy_top_in(:)
    real(kind=r8), intent(in) :: hflx_srf_in(:), hflx_top_in(:), cflx_srf_in(:,:), clfx_top_in(:,:)
    
    logical, intent(inout) :: lu_sum, lv_sum, ls_sum, lq_sum(pcnst)
    real(kind=r8), intent(inout) :: u_tend_sum(:,:), v_tend_sum(:,:), s_tend_sum(:,:), q_tend_sum(:,:,:)
    real(kind=r8), intent(inout) :: taux_srf_sum(:), taux_top_sum(:), tauy_srf_sum(:), tauy_top_sum(:)
    real(kind=r8), intent(inout) :: hflx_srf_sum(:), hflx_top_sum(:), cflx_srf_sum(:,:), clfx_top_sum(:,:)

!---------------------------Local storage-------------------------------
    integer :: i,k,m                               ! column,level,constituent indices

!-----------------------------------------------------------------------    

! Update u,v fields
    if(lu_in) then
       lu_sum = .true.
       do i = 1, ncol
          do k = top_level, bot_level
             u_tend_sum(i,k) = u_tend_sum(i,k) + u_tend_in(i,k)
          end do
          taux_srf_sum(i) = taux_srf_sum(i) + taux_srf_in(i)
          taux_top_sum(i) = taux_top_sum(i) + taux_top_in(i)
       end do
    end if

    if(lv_in) then
       lv_sum = .true.
       do i = 1, ncol
          do k = top_level, bot_level
             v_tend_sum(i,k) = v_tend_sum(i,k) + v_tend_in(i,k)
          end do
          tauy_srf_sum(i) = tauy_srf_sum(i) + tauy_srf_in(i)
          tauy_top_sum(i) = tauy_top_sum(i) + tauy_top_in(i)
       end do
    end if


    if(ls_in) then
       ls_sum = .true.
       do i = 1, ncol
          do k = top_level, bot_level
             s_tend_sum(i,k) = s_tend_sum(i,k) + s_tend_in(i,k)
          end do
          hflux_srf_sum(i) = hflux_srf_sum(i) + hflux_srf_in(i)
          hflux_top_sum(i) = hflux_top_sum(i) + hflux_top_in(i)
       end do
    end if

! Update constituents
    do m = 1, pcnst
       if(lq_in(m)) then
          lq_sum(m) = .true.
          do i = 1,ncol
             do k = top_level, bot_level
                q_tend_sum(i,k,m) = q_tend_sum(i,k,m) + q_tend_in(i,k,m)
             end do
             cflx_srf_sum(i,m) = cflx_srf_sum(i,m) + cflx_srf_in(i,m)
             cflx_top_sum(i,m) = cflx_top_sum(i,m) + cflx_top_in(i,m)
          end do
       end if
    end do

  end subroutine physics_ptend_sum
  
end module iap_ptend_sum