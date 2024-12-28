module zm_conv_momtran
  
  use shr_kind_mod,    only: r8 => shr_kind_r8
  
  implicit none
  
  public :: zm_conv_momtran_init, zm_conv_momtran_run

  private
  
  logical :: is_initialized = .False.
  integer, parameter :: ncnst = 2
  
  contains
!> \section arg_table_zm_conv_momtran_init Argument Table
!! \htmlinclude zm_conv_momtran_init.html
!!
  subroutine zm_conv_momtran_init(cam_physpkg, cam_physpkg_cam3, errmsg, errflg)
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
    
  end subroutine zm_conv_momtran_init

!> \section arg_table_zm_conv_momtran_run Argument Table
!! \htmlinclude zm_conv_momtran_run.html
!!
  subroutine zm_conv_momtran_run(ncol, pcols, pver, pverp, domomtran, u, v, mu, md, du, eu, ed, dp, dsubcld, jt, mx, ideep, il2g, dudt, dvdt, dsdt, pguall, pgdall, icwu, icwd, dt, errmsg, errflg)
    
    integer, intent(in) :: ncol, pcols, pver, pverp         ! number of atmospheric columns
    logical, intent(in) :: domomtran(ncnst)      ! flag for doing convective transport
    real(kind=r8), intent(in) :: dt
    real(kind=r8), intent(in) :: u(:,:)
    real(kind=r8), intent(in) :: v(:,:)
    real(kind=r8), intent(in) :: mu(:,:)       ! Mass flux up
    real(kind=r8), intent(in) :: md(:,:)       ! Mass flux down
    real(kind=r8), intent(in) :: du(:,:)       ! Mass detraining from updraft
    real(kind=r8), intent(in) :: eu(:,:)       ! Mass entraining from updraft
    real(kind=r8), intent(in) :: ed(:,:)       ! Mass entraining from downdraft
    real(kind=r8), intent(in) :: dp(:,:)       ! Delta pressure between interfaces
    real(kind=r8), intent(in) :: dsubcld(:)    ! Delta pressure from cloud base to sfc
    
    integer, intent(in) :: jt(:)         ! Index of cloud top for each column
    integer, intent(in) :: mx(:)         ! Index of cloud top for each column
    integer, intent(in) :: ideep(:)      ! Gathering array
    integer, intent(in) :: il2g              ! Gathered max lon indices over which to operate
    
    real(kind=r8), intent(inout) :: dudt(:,:), dvdt(:,:), dsdt(:,:)
    real(kind=r8), intent(out) ::  pguall(:,:,:)      ! Apparent force from  updraft PG
    real(kind=r8), intent(out) ::  pgdall(:,:,:)      ! Apparent force from  downdraft PG
    real(kind=r8), intent(out) ::  icwu(:,:,:)      ! In-cloud winds in updraft
    real(kind=r8), intent(out) ::  icwd(:,:,:)      ! In-cloud winds in downdraft
    
    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg
    
    integer, parameter :: il1g = 1   ! Gathered min lon indices over which to operate
    real(kind=r8) :: q(ncol,pver,ncnst) ! Wind array
    real(kind=r8) :: dqdt(ncol,pver,ncnst)  ! Tracer tendency array
    
    !--------------------------Local Variables------------------------------

    integer i                 ! Work index
    integer k                 ! Work index
    integer kbm               ! Highest altitude index of cloud base
    integer kk                ! Work index
    integer kkp1              ! Work index
    integer kkm1              ! Work index
    integer km1               ! Work index
    integer kp1               ! Work index
    integer ktm               ! Highest altitude index of cloud top
    integer m                 ! Work index
    integer ii                 ! Work index

    real(kind=r8) cabv                 ! Mix ratio of constituent above
    real(kind=r8) cbel                 ! Mix ratio of constituent below
    real(kind=r8) cdifr                ! Normalized diff between cabv and cbel
    real(kind=r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
    real(kind=r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
    real(kind=r8) const(pcols,pver)    ! Gathered wind array
    real(kind=r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
    real(kind=r8) dcondt(pcols,pver)   ! Gathered tend array
    real(kind=r8) small                ! A small number
    real(kind=r8) mbsth                ! Threshold for mass fluxes
    real(kind=r8) mupdudp              ! A work variable
    real(kind=r8) minc                 ! A work variable
    real(kind=r8) maxc                 ! A work variable
    real(kind=r8) fluxin               ! A work variable
    real(kind=r8) fluxout              ! A work variable
    real(kind=r8) netflux              ! A work variable

    real(kind=r8) momcu                ! constant for updraft pressure gradient term
    real(kind=r8) momcd                ! constant for downdraft pressure gradient term
    real(kind=r8) sum                  ! sum
    real(kind=r8) sum2                 ! sum2
     
    real(kind=r8) mududp(pcols,pver) ! working variable
    real(kind=r8) mddudp(pcols,pver)     ! working variable

    real(kind=r8) pgu(pcols,pver)      ! Pressure gradient term for updraft
    real(kind=r8) pgd(pcols,pver)      ! Pressure gradient term for downdraft
    
    real(kind=r8)                 gdsdt(pcols,pver) ! Gathered dry static energy tendency
    
    real(kind=r8)  mflux(pcols,pverp,ncnst)   ! Gathered momentum flux

    real(kind=r8)  wind0(pcols,pver,ncnst)       !  gathered  wind before time step
    real(kind=r8)  windf(pcols,pver,ncnst)       !  gathered  wind after time step
    real(kind=r8) fkeb, fket, ketend_cons, ketend, utop, ubot, vtop, vbot, gset2
    
    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    q(:ncol,:pver,1) = u(:ncol,:pver)
    q(:ncol,:pver,2) = v(:ncol,:pver)
    
! Initialize outgoing fields
    pguall(:,:,:)     = 0.0_r8
    pgdall(:,:,:)     = 0.0_r8
! Initialize in-cloud winds to environmental wind
    icwu(:ncol,:,:)   = q(:ncol,:,:)
    icwd(:ncol,:,:)   = q(:ncol,:,:)
! Initialize momentum flux and  final winds
    mflux(:,:,:)      = 0.0_r8
    wind0(:,:,:)      = 0.0_r8
    windf(:,:,:)      = 0.0_r8
! Initialize dry static energy
    dsdt(:,:)         = 0.0_r8
    gdsdt(:,:)        = 0.0_r8

! Define constants for parameterization
    momcu = 0.4_r8
    momcd = 0.4_r8
    small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
    mbsth = 1.e-15_r8

! Find the highest level top and bottom levels of convection
    ktm = pver
    kbm = pver
    do i = il1g, il2g
       ktm = min(ktm,jt(i))
       kbm = min(kbm,mx(i))
    end do
    
! Loop ever each wind component
    do m = 1, ncnst                    !start at m = 1 to transport momentum
      if (domomtran(m)) then
! Gather up the winds and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k)   = q(ideep(i),k,m)
               wind0(i,k,m) = const(i,k)
            end do
         end do

! From now on work only with gathered data

! Interpolate winds to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g
               ! use arithmetic mean
               chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))
! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)
!              provisional tends
               dcondt(i,k) = 0._r8
            end do
         end do
!
! Pressure Perturbation Term
! 
      !Top boundary:  assume mu is zero 
         k=1
         pgu(:il2g,k) = 0.0_r8
         pgd(:il2g,k) = 0.0_r8

         do k=2,pver-1
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g
               !interior points
               mududp(i,k) =  ( mu(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) &
                           +  mu(i,kp1) * (const(i,kp1) - const(i,k))/dp(i,k))
               pgu(i,k)    = - momcu * 0.5_r8 * mududp(i,k)
               mddudp(i,k) =  ( md(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) &
                           +  md(i,kp1) * (const(i,kp1) - const(i,k))/dp(i,k))
               pgd(i,k)    = - momcd * 0.5_r8 * mddudp(i,k)
            end do
         end do

         ! bottom boundary 
         k = pver
         km1 = max(1,k-1)
         do i=il1g,il2g
            mududp(i,k) =   mu(i,k) * (const(i,k)- const(i,km1))/dp(i,km1)
            pgu(i,k)    = - momcu *  mududp(i,k)
            mddudp(i,k) =   md(i,k) * (const(i,k)- const(i,km1))/dp(i,km1) 
            pgd(i,k)    = - momcd * mddudp(i,k)
         end do
!
! In-cloud velocity calculations
!

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         kkm1 = max(1,kk-1)
         do i = il1g,il2g
            mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
            if (mupdudp > mbsth) then 
               conu(i,kk) = (+eu(i,kk)*const(i,kk)*dp(i,kk)+pgu(i,kk)*dp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-ed(i,km1)*const(i,km1)*dp(i,km1))-pgd(i,km1)*dp(i,km1)/md(i,k)
            endif
         end do

! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkm1 = max(1,kk-1)
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + du(i,kk)*dp(i,kk)
               if (mupdudp > mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eu(i,kk)* &
                                  const(i,kk)*dp(i,kk)+pgu(i,kk)*dp(i,kk))/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-ed(i,km1)*const(i,km1) &
                                  *dp(i,km1)-pgd(i,km1)*dp(i,km1) )/md(i,k)
               endif
            end do
         end do

         sum = 0._r8
         sum2 = 0._r8

         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g
               ii = ideep(i)
! version 1 hard to check for roundoff errors
               dcondt(i,k) =  &
                           +(mu(i,kp1)* (conu(i,kp1)-chat(i,kp1)) &
                           -mu(i,k)*   (conu(i,k)-chat(i,k))      &
                           +md(i,kp1)* (cond(i,kp1)-chat(i,kp1)) &
                           -md(i,k)*   (cond(i,k)-chat(i,k)) &
                          )/dp(i,k)
            end do
         end do

  ! dcont for bottom layer
          !
!DIR$ NOINTERCHANGE
         do k = kbm,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k == mx(i)) then
                  ! version 1
                  dcondt(i,k) = (1./dp(i,k))*   &  
                       (-mu(i,k)*(conu(i,k)-chat(i,k)) &
                       -md(i,k)*(cond(i,k)-chat(i,k)) &
                       )
               end if
            end do
         end do

! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8

         do k = 1,pver
            do i = il1g,il2g
               ii = ideep(i)
               dqdt(ii,k,m) = dcondt(i,k)
    ! Output apparent force on the mean flow from pressure gradient
               pguall(ii,k,m) = -pgu(i,k)
               pgdall(ii,k,m) = -pgd(i,k)
               icwu(ii,k,m)   =  conu(i,k)
               icwd(ii,k,m)   =  cond(i,k)
            end do
         end do

          ! Calculate momentum flux in units of mb*m/s2 
         do k = ktm,pver
            do i = il1g,il2g
               ii = ideep(i)
               mflux(i,k,m) = &
                    -mu(i,k)*   (conu(i,k)-chat(i,k))      &
                    -md(i,k)*   (cond(i,k)-chat(i,k))
            end do
         end do

          ! Calculate winds at the end of the time step 
         do k = ktm,pver
            do i = il1g,il2g
               ii = ideep(i)
               km1 = max(1,k-1)
               kp1 = k+1
               windf(i,k,m) = const(i,k)    -   (mflux(i,kp1,m) - mflux(i,k,m)) * dt /dp(i,k)
            end do
         end do
      end if      ! for domomtran
    end do
   
! Need to add an energy fix to account for the dissipation of kinetic energy
  ! Formulation follows from Boville and Bretherton (2003)
  ! formulation by PJR

    do k = ktm,pver
      km1 = max(1,k-1)
      kp1 = min(pver,k+1)
      do i = il1g,il2g
        ii = ideep(i)

        ! calculate the KE fluxes at top and bot of layer 
        ! based on a discrete approximation to b&b eq(35) F_KE = u*F_u + v*F_v at interface
        utop = (wind0(i,k,1)+wind0(i,km1,1))/2.
        vtop = (wind0(i,k,2)+wind0(i,km1,2))/2.
        ubot = (wind0(i,kp1,1)+wind0(i,k,1))/2.
        vbot = (wind0(i,kp1,2)+wind0(i,k,2))/2.
        fket = utop*mflux(i,k,1)   + vtop*mflux(i,k,2)    ! top of layer
        fkeb = ubot*mflux(i,k+1,1) + vbot*mflux(i,k+1,2)  ! bot of layer

        ! divergence of these fluxes should give a conservative redistribution of KE
        ketend_cons = (fket-fkeb)/dp(i,k)

        ! tendency in kinetic energy resulting from the momentum transport
        ketend = ((windf(i,k,1)**2 + windf(i,k,2)**2) - (wind0(i,k,1)**2 + wind0(i,k,2)**2))*0.5/dt

        ! the difference should be the dissipation
        gset2 = ketend_cons - ketend
        gdsdt(i,k) = gset2
      end do
    end do

    ! Scatter dry static energy to full array
    do k = 1,pver
      do i = il1g,il2g
        ii = ideep(i)
        dsdt(ii,k) = gdsdt(i,k)
      end do
    end do

    dudt(:,:) = dqdt(:,:,1)
    dvdt(:,:) = dqdt(:,:,2)

    return

  end subroutine zm_conv_momtran_run
end module zm_conv_momtran