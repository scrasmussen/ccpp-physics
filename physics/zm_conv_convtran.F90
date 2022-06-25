module zm_conv_convtran_1_pre

  use shr_kind_mod,    only: r8 => shr_kind_r8

  implicit none

  public :: zm_conv_convtran_1_pre_run

  private

  contains

!! \section arg_table_zm_conv_convtran_1_pre_run
!! \htmlinclude zm_conv_convtran_1_pre_run.html
!!
  subroutine zm_conv_convtran_1_pre_run(ixcldice, ixcldliq, doconvtran_name, doconvtran_suv, doconvtran_q, dpdry, errmsg, errflg)

    integer,          intent(in)  :: ixcldice
    integer,          intent(in)  :: ixcldliq
    character(len=*), intent(out) :: doconvtran_name
    logical,          intent(out) :: doconvtran_suv(:) ! flag for doing convective transport of suv
    logical,          intent(out) :: doconvtran_q(:)   ! flag for doing convective transport of tracers
    real(kind=r8),    intent(out) :: dpdry(:,:)        ! used in convtran call

    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

    !GJF: these statements were in zm_conv_intr.F90/zm_conv_tend
    doconvtran_name        = 'convtran1'
    doconvtran_suv(1:3)    = .false.
    doconvtran_q(:)        = .false.
    doconvtran_q(ixcldice) = .true.
    doconvtran_q(ixcldliq) = .true.

    ! dpdry is not used in this call to convtran since the
    ! cloud liquid and ice mixing ratios are moist
    dpdry(:,:) = 0._r8

  end subroutine zm_conv_convtran_1_pre_run

end module zm_conv_convtran_1_pre


module zm_conv_convtran_2_pre

  use shr_kind_mod,    only: r8 => shr_kind_r8

  implicit none

  public :: zm_conv_convtran_2_pre_run

  private

  contains

!! \section arg_table_zm_conv_convtran_2_pre_run
!! \htmlinclude zm_conv_convtran_2_pre_run.html
!!
  subroutine zm_conv_convtran_2_pre_run(ideep, il2g, ixcldice, ixcldliq, pdeldry, doconvtran_name, doconvtran_suv, doconvtran_q, dpdry, state_q, temp_state_q, errmsg, errflg)

    integer,          intent(in)  :: ideep(:)          ! Gathering array
    integer,          intent(in)  :: il2g              ! Gathered max lon indices over which to operate
    integer,          intent(in)  :: ixcldice
    integer,          intent(in)  :: ixcldliq
    real(kind=r8),    intent(in)  :: pdeldry(:,:)
    character(len=*), intent(out) :: doconvtran_name
    logical,          intent(out) :: doconvtran_suv(:) ! flag for doing convective transport of suv
    logical,          intent(out) :: doconvtran_q(:)   ! flag for doing convective transport of tracers
    real(kind=r8),    intent(out) :: dpdry(:,:)        ! used in convtran call
    real(kind=r8),    intent(in)  :: state_q(:,:,:)
    real(kind=r8),    intent(out) :: temp_state_q(:,:,:)

    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg

    ! Local variables
    integer, parameter :: il1g = 1
    integer :: i

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

    !DH: these statements were in zm_conv_intr.F90/zm_conv_tend_2
    doconvtran_name        = 'convtran2'
    doconvtran_suv(1:3)    = .false.
    doconvtran_q(:)        = .true.
    doconvtran_q(1)        = .false.
    doconvtran_q(ixcldice) = .false.
    doconvtran_q(ixcldliq) = .false.

    ! initialize dpdry for call to convtran
    ! it is used for tracers of dry mixing ratio type
    dpdry = 0._r8
    do i = il1g,il2g
       !dpdry(i,:) = state%pdeldry(ideep(i,lchnk),:)/100._r8
       dpdry(i,:) = pdeldry(ideep(i),:)/100._r8
    end do

    ! Assign state tracer concentration to temporary tracer concentration
    temp_state_q = state_q

  end subroutine zm_conv_convtran_2_pre_run

end module zm_conv_convtran_2_pre


module zm_conv_convtran

  use shr_kind_mod,    only: r8 => shr_kind_r8

  implicit none

  public :: zm_conv_convtran_run

  private

  contains
!! \section arg_table_zm_conv_convtran_run
!! \htmlinclude zm_conv_convtran_run.html
!!
  subroutine zm_conv_convtran_run(pcols, pver, ixcldice, ixcldliq, cnst_type, q, ncnst, mu, md, du, eu, ed, dp, dsubcld, jt, mx, ideep, il2g, fracis, dqdt, doconvtran, dpdry, errmsg, errflg)

    integer, intent(in) :: pcols, pver
    integer, intent(in) :: ixcldice, ixcldliq
    character(len=3), intent(in) :: cnst_type(:)
    integer, intent(in) :: ncnst                 ! number of tracers to transport
    real(kind=r8), intent(in) :: q(:,:,:)  ! Tracer array including moisture
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
    integer, intent(in) :: il2g          ! Gathered max lon indices over which to operate

    real(kind=r8), intent(in) :: fracis(:,:,:)  ! fraction of tracer that is insoluble
    real(kind=r8), intent(inout) :: dqdt(:,:,:) ! Tracer tendency array
    logical, intent(in) :: doconvtran(:)        ! flag for doing convective transport
    real(kind=r8), intent(in) :: dpdry(:,:)     ! used in convtran call

    ! CCPP error handling
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg

!--------------------------Local Variables------------------------------

    integer i                 ! Work index
    integer k                 ! Work index
    integer kbm               ! Highest altitude index of cloud base
    integer kk                ! Work index
    integer kkp1              ! Work index
    integer km1               ! Work index
    integer kp1               ! Work index
    integer ktm               ! Highest altitude index of cloud top
    integer m                 ! Work index

    real(r8) cabv                 ! Mix ratio of constituent above
    real(r8) cbel                 ! Mix ratio of constituent below
    real(r8) cdifr                ! Normalized diff between cabv and cbel
    real(r8) chat(pcols,pver)     ! Mix ratio in env at interfaces
    real(r8) cond(pcols,pver)     ! Mix ratio in downdraft at interfaces
    real(r8) const(pcols,pver)    ! Gathered tracer array
    real(r8) fisg(pcols,pver)     ! gathered insoluble fraction of tracer
    real(r8) conu(pcols,pver)     ! Mix ratio in updraft at interfaces
    real(r8) dcondt(pcols,pver)   ! Gathered tend array
    real(r8) small                ! A small number
    real(r8) mbsth                ! Threshold for mass fluxes
    real(r8) mupdudp              ! A work variable
    real(r8) minc                 ! A work variable
    real(r8) maxc                 ! A work variable
    real(r8) fluxin               ! A work variable
    real(r8) fluxout              ! A work variable
    real(r8) netflux              ! A work variable

    real(r8) dutmp(pcols,pver)       ! Mass detraining from updraft
    real(r8) eutmp(pcols,pver)       ! Mass entraining from updraft
    real(r8) edtmp(pcols,pver)       ! Mass entraining from downdraft
    real(r8) dptmp(pcols,pver)    ! Delta pressure between interfaces
    !-----------------------------------------------------------------------
    !
    integer, parameter :: il1g = 1   ! Gathered min lon indices over which to operate

    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

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

! Loop ever each constituent
    do m = 2, ncnst
      if (doconvtran(m)) then

         if (cnst_type(m) .eq. 'dry') then !GJF: this should never execute because cloud ice and liquid mixing ratios are 'wet'
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dpdry(i,k)
                  dutmp(i,k) = du(i,k)*dp(i,k)/dpdry(i,k)
                  eutmp(i,k) = eu(i,k)*dp(i,k)/dpdry(i,k)
                  edtmp(i,k) = ed(i,k)*dp(i,k)/dpdry(i,k)
               end do
            end do
         else
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(i,k) = dp(i,k)
                  dutmp(i,k) = du(i,k)
                  eutmp(i,k) = eu(i,k)
                  edtmp(i,k) = ed(i,k)
               end do
            end do
         endif
!        dptmp = dp

! Gather up the constituent and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(i,k) = q(ideep(i),k,m)
               fisg(i,k) = fracis(ideep(i),k,m)
            end do
         end do

! From now on work only with gathered data

! Interpolate environment tracer values to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(i,km1),const(i,k))
               maxc = max(const(i,km1),const(i,k))
               if (minc < 0) then
                  cdifr = 0._r8
               else
                  cdifr = abs(const(i,k)-const(i,km1))/max(maxc,small)
               endif

! If the two layers differ significantly use a geometric averaging
! procedure
               if (cdifr > 1.E-6_r8) then
                  cabv = max(const(i,km1),maxc*1.e-12_r8)
                  cbel = max(const(i,k),maxc*1.e-12_r8)
                  chat(i,k) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(i,k) = 0.5_r8* (const(i,k)+const(i,km1))
               end if

! Provisional up and down draft values
               conu(i,k) = chat(i,k)
               cond(i,k) = chat(i,k)

!              provisional tends
               dcondt(i,k) = 0._r8

            end do
         end do

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         do i = il1g,il2g
            mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
            if (mupdudp > mbsth) then
               conu(i,kk) = (+eutmp(i,kk)*fisg(i,kk)*const(i,kk)*dptmp(i,kk))/mupdudp
            endif
            if (md(i,k) < -mbsth) then
               cond(i,k) =  (-edtmp(i,km1)*fisg(i,km1)*const(i,km1)*dptmp(i,km1))/md(i,k)
            endif
         end do

! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(i,kk) + dutmp(i,kk)*dptmp(i,kk)
               if (mupdudp > mbsth) then
                  conu(i,kk) = (  mu(i,kkp1)*conu(i,kkp1)+eutmp(i,kk)*fisg(i,kk)* &
                                  const(i,kk)*dptmp(i,kk) )/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(i,k) < -mbsth) then
                  cond(i,k) =  (  md(i,km1)*cond(i,km1)-edtmp(i,km1)*fisg(i,km1)*const(i,km1) &
                                  *dptmp(i,km1) )/md(i,k)
               endif
            end do
         end do


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g

! version 1 hard to check for roundoff errors
!               dcondt(i,k) =
!     $                  +(+mu(i,kp1)* (conu(i,kp1)-chat(i,kp1))
!     $                    -mu(i,k)*   (conu(i,k)-chat(i,k))
!     $                    +md(i,kp1)* (cond(i,kp1)-chat(i,kp1))
!     $                    -md(i,k)*   (cond(i,k)-chat(i,k))
!     $                   )/dp(i,k)

! version 2 hard to limit fluxes
!               fluxin =  mu(i,kp1)*conu(i,kp1) + mu(i,k)*chat(i,k)
!     $                 -(md(i,k)  *cond(i,k)   + md(i,kp1)*chat(i,kp1))
!               fluxout = mu(i,k)*conu(i,k)     + mu(i,kp1)*chat(i,kp1)
!     $                 -(md(i,kp1)*cond(i,kp1) + md(i,k)*chat(i,k))

! version 3 limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
               fluxin =  mu(i,kp1)*conu(i,kp1)+ mu(i,k)*min(chat(i,k),const(i,km1)) &
                         -(md(i,k)  *cond(i,k) + md(i,kp1)*min(chat(i,kp1),const(i,kp1)))
               fluxout = mu(i,k)*conu(i,k) + mu(i,kp1)*min(chat(i,kp1),const(i,k)) &
                         -(md(i,kp1)*cond(i,kp1) + md(i,k)*min(chat(i,k),const(i,k)))

               netflux = fluxin - fluxout
               if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                  netflux = 0._r8
               endif
               dcondt(i,k) = netflux/dptmp(i,k)
            end do
         end do
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!DIR$ NOINTERCHANGE
         do k = kbm,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k == mx(i)) then

! version 1
!                  dcondt(i,k) = (1./dsubcld(i))*
!     $              (-mu(i,k)*(conu(i,k)-chat(i,k))
!     $               -md(i,k)*(cond(i,k)-chat(i,k))
!     $              )

! version 2
!                  fluxin =  mu(i,k)*chat(i,k) - md(i,k)*cond(i,k)
!                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*chat(i,k)
! version 3
                  fluxin =  mu(i,k)*min(chat(i,k),const(i,km1)) - md(i,k)*cond(i,k)
                  fluxout = mu(i,k)*conu(i,k) - md(i,k)*min(chat(i,k),const(i,k))

                  netflux = fluxin - fluxout
                  if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                     netflux = 0._r8
                  endif
!                  dcondt(i,k) = netflux/dsubcld(i)
                  dcondt(i,k) = netflux/dptmp(i,k)
               else if (k > mx(i)) then
!                  dcondt(i,k) = dcondt(i,k-1)
                  dcondt(i,k) = 0._r8
               end if
            end do
         end do

! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(:,:,m) = 0._r8
         do k = 1,pver
            kp1 = min(pver,k+1)
!DIR$ CONCURRENT
            do i = il1g,il2g
               dqdt(ideep(i),k,m) = dcondt(i,k)
            end do
         end do

      end if      ! for doconvtran

    end do

    return

  end subroutine zm_conv_convtran_run
end module zm_conv_convtran
