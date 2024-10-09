! ###########################################################################################
!
! ###########################################################################################
module gfs_mmm_cu_ntiedtke_pre
  use machine, only: kind_phys

  implicit none

  public gfs_mmm_cu_ntiedtke_pre_run

contains
!> \section arg_table_gfs_mmm_cu_ntiedtke_pre_run
!! \htmlinclude gfs_mmm_cu_ntiedtke_pre_run.html
!!
  ! #########################################################################################
  !
  ! #########################################################################################
  subroutine gfs_mmm_cu_ntiedtke_pre_run(nCol, nLev, temp, spechum, qv, temp_new, spechum_new,&
                                         flag_init, flag_restart, kdt, fhour, dtp, prevst, prevsq, forcet, forceq, &
                                         qc, qi, ugrs, vgrs, pres, presi, &
                                         geop, geopi, geoph, geophi, con_1ovg,& 
                                         pomg, wet, dry, icy, ep1, con_rd, con_cp, ps,  &
                                         hfx_wat, hfx_lnd, hfx_ice, hfx, qfx_wat, qfx_lnd, qfx_ice, qfx, &
                                         slimask, xland, errmsg, errflg)

    ! Input variables
    logical,                         intent(in   ) :: flag_init, flag_restart

    integer,                         intent(in   ) :: &
        nCol,   & ! Number of horizontal gridpoints
        nLev,   & ! Number of vertical levels
        kdt       ! Index of time step for current iteration

    integer,           dimension(:), intent(in   ) :: slimask  ! landmask: sea/land/ice=0/1/2
                                                                                                            
    real(kind_phys),                 intent(in   ) :: &
        fhour,   & ! Current forecast time (hour)
        dtp,     & ! Physics timestep (second)
        con_1ovg,& ! 1/g
        ep1,     & ! (rv/rd) - 1
        con_rd,  & ! rd
        con_cp     ! cp

    real(kind_phys),   dimension(:), intent(in   ) :: &
        hfx_wat,        & ! Kinematic surface upward sensible heat flux over water (K m s-1)
        hfx_lnd,        & ! Kinematic surface upward sensible heat flux over land (K m s-1)
        hfx_ice,        & ! Kinematic surface upward sensible heat flux over ice (K m s-1)
        qfx_wat,        & ! Kinematic surface upward latent heat flux over water (kg kg-1 m s-1)
        qfx_lnd,        & ! Kinematic surface upward latent heat flux over land (kg kg-1 m s-1)
        qfx_ice,        & ! Kinematic surface upward latent heat flux over ice (kg kg-1 m s-1)
        ps                ! Surface pressure (Pa)

    ! Variables to be vertically reversed    
    real(kind_phys), dimension(:,:), intent(in   ) :: &
        temp,           & ! Air temperature (K)
        spechum,        & ! Specific humidity (kg/kg) 
        geop,           & ! Geopotential (m2 s-2)
        geopi             ! Geopotential at model layer interface (m2 s-2)
    real(kind_phys), dimension(:,:), intent(inout) :: &
        temp_new,       & ! Air temperature of new state (K)
        spechum_new,    & ! Specific humidity of new state (kg/kg)
        qc,             & ! Cloud liquid water mixing ratio of new state (kg/kg)
        qi,             & ! Cloud ice mixing ratio of new state (kg/kg)
        ugrs,           & ! x-wind of new state (m/s)
        vgrs,           & ! y-wind of new state (m/s)
        pres,           & ! Air pressure (Pa)
        presi,          & ! Air Pressure at model layer interface (Pa)
        pomg              ! Layer mean vertical velocity (Pa s-1)

    real(kind_phys),   dimension(:,:), intent(in), optional :: prevst, prevsq ! t and q from previous time step

    logical,         dimension(:),  intent(in)   :: &
        wet,            & ! Flag indicating presence of some ocean or lake surface area fraction
        dry,            & ! Flag indicating presence of some land surface area fraction
        icy               ! Flag indicating presence of some sea ice surface area fraction 

    ! Output variables
    real(kind_phys), dimension(:),  intent(out  ), optional:: &
        hfx,            & ! Surface kinematic upward sensible heat flux (W m-2)
        qfx,            & ! Surface kinematic upward latent heat flux (W m-2)
        xland             ! Sea/land/ice mask for MMM physics
    real(kind_phys), dimension(:,:),intent(out  ), optional:: &
        qv,             & ! Water vapor mixing ratio (kg/kg)
        !pqc,            & ! Cloud liquid water mixing ratio of new state transported by convection (kg/kg)
        !pqi,            & ! Cloud ice mixing ratio of new state transported by convection (kg/kg)
        forcet,         & ! Temperature tendency due to dynamics only (K s-1)
        forceq,         & ! Moisture tendency due to dynamics only (kg kg-1 s-1)
        geoph,          & ! Geopotential height (m)
        geophi            ! Geopoteitial height at model interface (m)
    character(len=*),                intent(out) :: &
        errmsg                 ! CCPP error message
    integer,                         intent(out) :: &
        errflg                 ! CCPP error code

    ! Local variables
    integer                                :: i, k, kk
    real(kind_phys)                        :: tvcon, rho, dtdyn
    real(kind_phys), dimension(nCol, nLev) :: prevsqv, pt, qv_new, pqv, pqc, pqi, pu, pv, prsl, zl, omega, tendt, tendq
    real(kind_phys), dimension(nCol,nLev+1):: prsli, zi

    ! Initialize CCPP error handling
    errmsg = ''
    errflg = 0

    !> convert specific humidity to water vapor mixing ratio
    do k=1,nLev
        do i=1,nCol
            qv_new(i,k) = spechum_new(i,k)/ (1.0_kind_phys - spechum_new(i,k))! to be vertically reversed
            qv(i,k)     = spechum(i,k)    / (1.0_kind_phys - spechum(i,k))    ! to calculate dynamics tendencies
            prevsqv(i,k)= prevsq(i,k)     / (1.0_kind_phys - prevsq(i,k))     ! to calculate dynamics tendencies
        end do
    end do

    ! Calculate tendencies due to nonphys (incl. dynamics+PBL for moisture and dynamics+PBL+radiation for temperature)
      if(flag_init .and. .not.flag_restart) then
        forcet=0.0
        forceq=0.0
      else
        dtdyn=3600.0*(fhour)/kdt
        if(dtp > dtdyn) then
            forcet=(temp - temp_new)/dtp!prevst)/dtp
            forceq=(qv - qv_new)/dtp!prevsqv)/dtp
        else
            forcet=(temp - temp_new)/dtdyn!prevst)/dtdyn
            forceq=(qv - qv_new)/dtdyn!prevsqv)/dtdyn
        endif
      endif

    !> convert heat fluxes
    do i = 1,nCol
       tvcon = (1. + ep1 * spechum(i,1))
       rho   = ps(i) / (con_rd * temp(i,1) * tvcon)
        if (dry(i)) then
            hfx(i)  = rho * con_cp * hfx_lnd(i)
            qfx(i)  = rho * qfx_lnd(i)
        end if
        if (wet(i)) then
            hfx(i)  = rho * con_cp * hfx_wat(i)
            qfx(i)  = rho * qfx_wat(i)
        end if

        if (icy(i)) then
            hfx(i)  = rho * con_cp * hfx_ice(i)
            qfx(i)  = rho * qfx_ice(i)
        end if
    end do

    ! Compute land/sea mask convention from (0-sea/1-land/2-ice) ---> (1-land/2-sea)
    do i=1,nCol
        if(slimask(i).eq.0)then
            xland(i) = 2
        else
            xland(i) = 1
        end if
    enddo

    ! Reverse vertical layers for mass flux calculation
    do k=1,nLev
        kk=nLev-k+1
        do i=1,nCol
            pt(i,k)    = temp_new(i,kk)
            pqv(i,k)   = qv_new(i,kk)
            pqc(i,k)   = qc(i,kk)
            pqi(i,k)   = qi(i,kk)
            pu(i,k)    = ugrs(i,kk)
            pv(i,k)    = vgrs(i,kk)
            prsl(i,k)  = pres(i,kk)
            prsli(i,k) = presi(i,kk+1)
            zl(i,k)    = geop(i,kk) * con_1ovg   ! calculate geopotential height
            zi(i,k)    = geopi(i,kk+1) * con_1ovg! calculate geopotnetial height
            omega(i,k) = pomg(i,kk)
            tendt(i,k) = forcet(i,kk)
            tendq(i,k) = forceq(i,kk)
        enddo
    enddo
    prsli(:,nLev+1)=presi(:,1)
       zi(:,nLev+1)=geopi(:,1) * con_1ovg

    ! Output
    temp_new= pt
    qv      = pqv
    qc      = pqc
    qi      = pqi
    ugrs    = pu
    vgrs    = pv 
    pres    = prsl
    presi   = prsli
    geoph   = zl
    geophi  = zi
    pomg    = omega
    forcet  = tendt
    forceq  = tendq

  end subroutine gfs_mmm_cu_ntiedtke_pre_run

end module gfs_mmm_cu_ntiedtke_pre
