module zm_conv_momtran_post
  
  use shr_kind_mod,    only: r8 => shr_kind_r8
  
  implicit none
  
  public :: zm_conv_momtran_post_init, zm_conv_momtran_post_run

  private
  
  logical :: is_initialized = .False.
  
  contains
!> \section arg_table_zm_conv_momtran_post_init Argument Table
!! \htmlinclude zm_conv_momtran_post_init.html
!!
  subroutine zm_conv_momtran_post_init(cam_physpkg, cam_physpkg_cam3, errmsg, errflg)
    character(len=16), intent(in) :: cam_physpkg, cam_physpkg_cam3
    
    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg
    
    if (is_initialized) return
    
    ! Consistency checks
    if (cam_physpkg == cam_physpkg_cam3) then
       write(errmsg,'(*(a))') "Logic error: choice of physics package (cam3) is incompatible with Zhang-McFarlane convective momentum transport."
       errflg = 1
       return
    end if
    
  end subroutine zm_conv_momtran_post_init

!> \section arg_table_zm_conv_momtran_post_run Argument Table
!! \htmlinclude zm_conv_momtran_post_run.html
!!
  subroutine zm_conv_momtran_post_run(ncol, pcols, pver, pverp, pcnst, ixcldice, ixcldliq, ixnumice, ixnumliq, &
    gravit, rair, cpair, zvir, dt, fv_dycore, microp_scheme, qmin, lnpint, lnpmid, pint, pmid, pdel, rpdel, phis, &
    tend_u, tend_v, tend_s, temp_state_u, temp_state_v, temp_state_s, temp_state_q, temp_state_t, temp_state_zm, &
    temp_state_zi, errmsg, errflg)
    
    use iap_state_update, only : physics_update
    
    integer,           intent(in   )                   :: ncol, pcols, pver, pverp, pcnst, &
                                                          ixcldice, ixcldliq, ixnumice, ixnumliq
    logical,           intent(in   )                   :: fv_dycore
    character(len=16), intent(in   )                   :: microp_scheme
    real(kind=r8),     intent(in   )                   :: gravit, rair, cpair, zvir, dt
    real(kind=r8),     intent(in   ), dimension(:)     :: qmin
    real(kind=r8),     intent(in   ), dimension(:)     :: phis
    real(kind=r8),     intent(in   ), dimension(:,:)   :: lnpint, pint
    real(kind=r8),     intent(in   ), dimension(:,:)   :: lnpmid, pmid, pdel, rpdel
    real(kind=r8),     intent(in   ), dimension(:,:)   :: tend_u, tend_v, tend_s
    real(kind=r8),     intent(inout), dimension(:,:)   :: temp_state_u, temp_state_v, temp_state_s
    real(kind=r8),     intent(inout), dimension(:,:,:) :: temp_state_q
    real(kind=r8),     intent(inout), dimension(:,:)   :: temp_state_t, temp_state_zm, temp_state_zi
    
    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg
    
    integer :: top, bot
    logical :: lu, lv, ls
    logical, dimension(pcnst) :: lq
    character*24 :: name
    real(kind=r8), dimension(pcols,pver,pcnst) :: tend_q
    
    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

    lu = .true.
    lv = .true.
    ls = .true.
    lq(:) = .false.
    name = 'zm_conv_momtran'
    top = 1
    bot = pver
    tend_q(:,:,:) = 0.0_r8
    
    call physics_update(lu, lv, ls, lq, fv_dycore, name, microp_scheme, top, bot, ncol, &
                        pcols, pver, pverp, pcnst, ixcldice, ixcldliq, ixnumice, ixnumliq,     &
                        rair, gravit, cpair, zvir, dt, qmin, tend_u, tend_v, tend_s, tend_q, &
                        lnpint, lnpmid, pint, pmid, pdel, rpdel, phis, temp_state_u, &
                        temp_state_v, temp_state_s, temp_state_q, temp_state_t, temp_state_zm, &
                        temp_state_zi)
    
  end subroutine zm_conv_momtran_post_run
end module zm_conv_momtran_post
