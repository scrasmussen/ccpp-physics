module zm_convr_post
  
  use shr_kind_mod,    only: r8 => shr_kind_r8
  
  implicit none
  
  contains

!> \section arg_table_zm_convr_post_run Argument Table
!! \htmlinclude zm_convr_post_run.html
!!
  subroutine zm_convr_post_run(ncol, pcols, pver, pverp, pcnst, lengath, ixcldice, ixcldliq, ixnumice, ixnumliq, gravit, rair, cpair, zvir, dt, jt, maxg, ideep,   &
    fv_dycore, microp_scheme, qmin, ps, lnpint, lnpmid, pint, pmid, pdel, rpdel, phis, tend_q, tend_s, state_u, state_v, state_s, state_q, mcon, pcont, pconb, dp_cldliq, dp_cldice, &
    temp_state_u, temp_state_v, temp_state_s, temp_state_q, temp_state_t, temp_state_zm, temp_state_zi, errmsg, errflg)
    
    use iap_state_update, only : physics_update
    
    integer,           intent(in   )                   :: ncol, pcols, pver, pverp, pcnst, lengath, &
                                                          ixcldice, ixcldliq, ixnumice, ixnumliq
    integer,           intent(in   ), dimension(:)     :: jt, maxg, ideep
    logical,           intent(in   )                   :: fv_dycore
    character(len=16), intent(in   )                   :: microp_scheme
    real(kind=r8),     intent(in   )                   :: gravit, rair, cpair, zvir, dt
    real(kind=r8),     intent(in   ), dimension(:)     :: qmin
    real(kind=r8),     intent(in   ), dimension(:)     :: ps, phis
    real(kind=r8),     intent(in   ), dimension(:,:)   :: lnpint, pint
    real(kind=r8),     intent(in   ), dimension(:,:)   :: lnpmid, pmid, pdel, rpdel, tend_q, tend_s
    real(kind=r8),     intent(in   ), dimension(:,:)   :: state_u, state_v, state_s
    real(kind=r8),     intent(in   ), dimension(:,:,:) :: state_q
    real(kind=r8),     intent(inout), dimension(:,:)   :: mcon
    real(kind=r8),     intent(  out), dimension(:)     :: pcont, pconb
    real(kind=r8),     intent(  out), dimension(:,:)   :: dp_cldliq, dp_cldice
    real(kind=r8),     intent(  out), dimension(:,:)   :: temp_state_u, temp_state_v, temp_state_s
    real(kind=r8),     intent(  out), dimension(:,:,:) :: temp_state_q
    real(kind=r8),     intent(  out), dimension(:,:)   :: temp_state_t, temp_state_zm, temp_state_zi
    
    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg
    
    integer :: i, top, bot
    logical :: lu, lv, ls
    logical, dimension(pcnst) :: lq
    character*24 :: name
    real(kind=r8), dimension(pcols,pver) :: tend_u, tend_v
    real(kind=r8), dimension(pcols,pver,pcnst) :: tend_q_all    

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    !
    ! Convert mass flux from reported mb/s to kg/m^2/s
    !
    mcon(:ncol,:pver) = mcon(:ncol,:pver) * 100._r8/gravit
    
    pcont(:ncol) = ps(:ncol)
    pconb(:ncol) = ps(:ncol)
    do i = 1,lengath
      if (maxg(i).gt.jt(i)) then
        pcont(ideep(i)) = pmid(ideep(i),jt(i))  ! gathered array (or jctop ungathered)
        pconb(ideep(i)) = pmid(ideep(i),maxg(i))! gathered array
      endif
    end do
    
    ! update temporary physics state variables (for use in subsequent zm_conv routines) with tendencies coming from zm_convr 
    !t, zm, zi - from state1 var; need new vars in phys_int_ephem to hold
    lu = .false.
    lv = .false.
    ls = .true.
    lq(1) = .true.
    lq(2:) = .false.
    name = 'zm_convr'
    top = 1
    bot = pver
    tend_u = 0._r8
    tend_v = 0._r8
    tend_q_all(:,:,:) = 0._r8
    tend_q_all(:,:,1) = tend_q    

    !GJF: replaces the physics_state_copy(state,state1) call in zm_conv_intr.F90/zm_conv_tend
    temp_state_u(:,:) = state_u(:,:)
    temp_state_v(:,:) = state_v(:,:)
    temp_state_s(:,:) = state_s(:,:)
    temp_state_q(:,:,:) = state_q(:,:,:)
    
    call physics_update(lu, lv, ls, lq, fv_dycore, name, microp_scheme, top, bot, ncol, &
                        pcols, pver, pverp, pcnst, ixcldice, ixcldliq, ixnumice, ixnumliq,     &
                        rair, gravit, cpair, zvir, dt, qmin, tend_u, tend_v, tend_s, tend_q_all, &
                        lnpint, lnpmid, pint, pmid, pdel, rpdel, phis, temp_state_u, &
                        temp_state_v, temp_state_s, temp_state_q, temp_state_t, temp_state_zm, &
                        temp_state_zi)
    
    dp_cldliq(:ncol,:) = 0._r8
    dp_cldice(:ncol,:) = 0._r8
    
  end subroutine zm_convr_post_run
end module zm_convr_post
