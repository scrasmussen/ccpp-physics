module iap_state_update
  
  use shr_kind_mod, only: r8=>shr_kind_r8
 
  implicit none

  public
  
  contains
  
  subroutine physics_update(lu, lv, ls, lq, fv_dycore, name, microp_scheme, top, bot, ncol, pcols, pver, pverp, pcnst, ixcldice, ixcldliq, ixnumice, ixnumliq, rair, gravit, cpair, zvir, dt, qmin, tend_u, tend_v, tend_s, tend_q, lnpint, lnpmid, pint, pmid, pdel, rpdel, phis, state_u, state_v, state_s, state_q, t, zm, zi)
    logical, intent(in) :: lu, lv, ls, fv_dycore
    logical, intent(in), dimension(:) :: lq
    character*24, intent(in) :: name    ! name of parameterization which produced tendencies.
    character(len=16), intent(in) :: microp_scheme
    integer, intent(in) :: top, bot, ncol, pcols, pver, pverp, pcnst
    integer, intent(in) :: ixcldice, ixcldliq                  ! indices for CLDICE and CLDLIQ
    integer, intent(in) :: ixnumice, ixnumliq                  ! indices for ice and liquid number concentrations
    real(kind=r8), intent(in   )                 :: rair, gravit, cpair, zvir
    
    real(kind=r8), intent(in   )                 :: dt
    real(kind=r8), intent(in   ), dimension(:)   :: qmin
    real(kind=r8), intent(in   ), dimension(:,:) :: tend_u, tend_v, tend_s
    real(kind=r8), intent(in   ), dimension(:,:,:) :: tend_q
    
    real(kind=r8), intent(in   ), dimension(:,:) :: lnpint, pint
    real(kind=r8), intent(in   ), dimension(:,:) :: lnpmid, pmid, pdel, rpdel
    real(kind=r8), intent(in   ), dimension(:)   :: phis
    
    real(kind=r8), intent(inout), dimension(:,:) :: state_u, state_v, state_s
    real(kind=r8), intent(inout), dimension(:,:,:) :: state_q
    real(kind=r8), intent(inout), dimension(:,:) :: t, zm
    real(kind=r8), intent(inout), dimension(:,:) :: zi
    
    integer :: i,k,m
    
    ! Update u,v fields
    if(lu) then
       do k = top, bot
          do i = 1, ncol
             state_u(i,k) = state_u(i,k) + tend_u(i,k) * dt
          end do
       end do
    end if

    if(lv) then
       do k = top, bot
          do i = 1, ncol
             state_v(i,k) = state_v(i,k) + tend_v(i,k) * dt
          end do
       end do
    end if
    
    ! Update dry static energy
    if(ls) then
       do k = top, bot
          do i = 1, ncol
             state_s(i,k) = state_s(i,k) + tend_s(i,k) * dt
          end do
       end do
    end if
    
    do m = 1, pcnst
       if(lq(m)) then
          do k = top, bot
             do i = 1,ncol
                state_q(i,k,m) = state_q(i,k,m) + tend_q(i,k,m) * dt
             end do
          end do

          ! now test for mixing ratios which are too small
          ! don't call qneg3 for number concentration variables
          if (m .ne. ixnumice  .and.  m .ne. ixnumliq) then
            call qneg3(ncol, pcols, pver, m, m, qmin, state_q)
          else
             do k = top, bot
                do i = 1,ncol
                   ! checks for number concentration
                   state_q(i,k,m) = max(1.e-12_r8,state_q(i,k,m))
                   state_q(i,k,m) = min(1.e10_r8,state_q(i,k,m))
                end do
             end do
          end if
       end if
    end do
    
    ! special tests for cloud liquid
    if (ixcldliq > 1) then
       if(lq(ixcldliq)) then
          if (name == 'stratiform' .or. name == 'cldwat'  ) then
#ifdef PERGRO
             where (state_q(:ncol,:pver,ixcldliq) < 1.e-12_r8)
                state_q(:ncol,:pver,ixcldliq) = 0._r8
             end where

             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' .or. microp_scheme .eq. 'M2M') then
                where (state_q(:ncol,:pver,ixcldliq) < 1.e-12_r8)
                   state_q(:ncol,:pver,ixnumliq) = 0._r8
                end where
             end if
#endif
          else if (name == 'convect_deep') then
             where (state_q(:ncol,:pver,ixcldliq) < 1.e-36_r8)
                state_q(:ncol,:pver,ixcldliq) = 0._r8
             end where
             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' .or. microp_scheme .eq. 'M2M' ) then
                where (state_q(:ncol,:pver,ixcldliq) < 1.e-36_r8)
                   state_q(:ncol,:pver,ixnumliq) = 0._r8
                end where
             end if
          end if
       end if
    end if

    ! special tests for cloud ice
    if (ixcldice > 1) then
       if(lq(ixcldice)) then
          if (name == 'stratiform' .or. name == 'cldwat'  ) then
#ifdef PERGRO
             where (state_q(:ncol,:pver,ixcldice) < 1.e-12_r8)
                state_q(:ncol,:pver,ixcldice) = 0._r8
             end where
             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' .or. microp_scheme .eq. 'M2M') then
                where (state_q(:ncol,:pver,ixcldice) < 1.e-12_r8)
                   state_q(:ncol,:pver,ixnumice) = 0._r8
                end where
             end if
#endif
          else if (name == 'convect_deep') then
             where (state_q(:ncol,:pver,ixcldice) < 1.e-36_r8)
                state_q(:ncol,:pver,ixcldice) = 0._r8
             end where
             ! also zero out number concentration
             if ( microp_scheme .eq. 'MG' .or. microp_scheme .eq. 'M2M' ) then
                where (state_q(:ncol,:pver,ixcldice) < 1.e-36_r8)
                   state_q(:ncol,:pver,ixnumice) = 0._r8
                end where
             end if
          end if
       end if
    end if
    
    if (ls .or. lq(1)) then
       call geopotential_dse(pcols, pver, pverp,                 &
            lnpint, lnpmid, pint, pmid, pdel, rpdel,            &
            state_s, state_q(:ncol,:pver,1), phis, rair, gravit, cpair, &
            zvir, t, zi, zm, ncol, fv_dycore)
    end if
    
  end subroutine physics_update
  
  subroutine qneg3(ncol    ,ncold   ,lver    ,lconst_beg  , &
                    lconst_end       ,qmin    ,q       )
    
    integer, intent(in) :: ncol         ! number of atmospheric columns
    integer, intent(in) :: ncold        ! declared number of atmospheric columns
    integer, intent(in) :: lver         ! number of vertical levels in column
    integer, intent(in) :: lconst_beg   ! beginning constituent
    integer, intent(in) :: lconst_end   ! ending    constituent

    real(r8), intent(in) :: qmin(lconst_beg:lconst_end)      ! Global minimum constituent concentration

 !
 ! Input/Output arguments
 !
    real(r8), intent(inout) :: q(ncold,lver,lconst_beg:lconst_end) ! moisture/tracer field
 !
 !---------------------------Local workspace-----------------------------
 !
    integer indx(ncol,lver)  ! array of indices of points < qmin
    integer nval(lver)       ! number of points < qmin for 1 level
    integer nvals            ! number of values found < qmin
    integer nn
    integer iwtmp
    integer i,ii,k           ! longitude, level indices
    integer m                ! constituent index
    integer iw,kw            ! i,k indices of worst violator

    logical found            ! true => at least 1 minimum violator found

    real(r8) worst           ! biggest violator
 !
 !-----------------------------------------------------------------------
 !

    do m=lconst_beg,lconst_end
       nvals = 0
       found = .false.
       worst = 1.e35_r8
       iw = -1
 !
 ! Test all field values for being less than minimum value. Set q = qmin
 ! for all such points. Trace offenders and identify worst one.
 !
!DIR$ preferstream
       do k=1,lver
          nval(k) = 0
!DIR$ prefervector
          nn = 0
          do i=1,ncol
             if (q(i,k,m) < qmin(m)) then
                nn = nn + 1
                indx(nn,k) = i
             end if
          end do
          nval(k) = nn
       end do

       do k=1,lver
          if (nval(k) > 0) then
             found = .true.
             nvals = nvals + nval(k)
             iwtmp = -1
 !cdir nodep,altcode=loopcnt
             do ii=1,nval(k)
                i = indx(ii,k)
                if (q(i,k,m) < worst) then
                   worst = q(i,k,m)
                   iwtmp = ii
                end if
             end do
             if (iwtmp /= -1 ) kw = k
             if (iwtmp /= -1 ) iw = indx(iwtmp,k)
 !cdir nodep,altcode=loopcnt
             do ii=1,nval(k)
                i = indx(ii,k)
                q(i,k,m) = qmin(m)
             end do
          end if
       end do
    end do
 !
    return
  end subroutine qneg3
  
  subroutine geopotential_dse(pcols, pver, pverp,             &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  ,  &
       dse    , q      , phis   , rair   , gravit , cpair  ,  &
       zvir   , t      , zi     , zm     , ncol   , fvdyn)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the temperature  and geopotential height (above the surface) at the
! midpoints and interfaces from the input dry static energy and pressures.
!
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!
! Input arguments
    integer, intent(in) :: pcols, pver, pverp
    integer, intent(in) :: ncol                  ! Number of longitudes
    logical, intent(in) :: fvdyn                 ! True if using finite volume dynamics

    real(r8), intent(in) :: piln (pcols,pverp)   ! Log interface pressures
    real(r8), intent(in) :: pmln (pcols,pver)    ! Log midpoint pressures
    real(r8), intent(in) :: pint (pcols,pverp)   ! Interface pressures
    real(r8), intent(in) :: pmid (pcols,pver)    ! Midpoint pressures
    real(r8), intent(in) :: pdel (pcols,pver)    ! layer thickness
    real(r8), intent(in) :: rpdel(pcols,pver)    ! inverse of layer thickness
    real(r8), intent(in) :: dse  (pcols,pver)    ! dry static energy
    real(r8), intent(in) :: q    (pcols,pver)    ! specific humidity
    real(r8), intent(in) :: phis (pcols)         ! surface geopotential
    real(r8), intent(in) :: rair                 ! Gas constant for dry air
    real(r8), intent(in) :: gravit               ! Acceleration of gravity
    real(r8), intent(in) :: cpair                ! specific heat at constant p for dry air
    real(r8), intent(in) :: zvir                 ! rh2o/rair - 1

! Output arguments

    real(r8), intent(out) :: t(pcols,pver)       ! temperature
    real(r8), intent(out) :: zi(pcols,pverp)     ! Height above surface at interfaces
    real(r8), intent(out) :: zm(pcols,pver)      ! Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
!
    integer  :: i,k                ! Lon, level, level indices
    real(r8) :: hkk(pcols)         ! diagonal element of hydrostatic matrix
    real(r8) :: hkl(pcols)         ! off-diagonal element
    real(r8) :: rog                ! Rair / gravit
    real(r8) :: tv                 ! virtual temperature
    real(r8) :: tvfac              ! Tv/T
!
!-----------------------------------------------------------------------
    rog = rair/gravit

! The surface height is zero by definition.
    do i = 1,ncol
       zi(i,pverp) = 0.0_r8
    end do

! Compute the virtual temperature, zi, zm from bottom up
! Note, zi(i,k) is the interface above zm(i,k)
    do k = pver, 1, -1

! First set hydrostatic elements consistent with dynamics
       if (fvdyn) then
          do i = 1,ncol
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1._r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          end do
       else
          do i = 1,ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          end do
       end if

! Now compute tv, t, zm, zi
       do i = 1,ncol
          tvfac   = 1._r8 + zvir * q(i,k)
          tv      = (dse(i,k) - phis(i) - gravit*zi(i,k+1)) / ((cpair / tvfac) + rair*hkk(i))

          t (i,k) = tv / tvfac

          zm(i,k) = zi(i,k+1) + rog * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog * tv * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_dse
  
end module iap_state_update
