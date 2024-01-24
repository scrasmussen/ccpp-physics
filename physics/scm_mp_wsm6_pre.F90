!>\file wsm6.F90
!! This file runs the Single-Moment 6-class Microphysics scheme (WSM6)


!>\defgroup aathompson Aerosol-Aware Thompson MP Module
!! This module runs the Single-Moment 6-class Microphysics scheme (WSM6)
module scm_mp_wsm6_pre

      use ccpp_kinds, only : kind_phys

      use mp_wsm6, only : mp_wsm6_init
      ! use module_mp_thompson, only : naIN0, naIN1, naCCN0, naCCN1, eps, Nt_c_l, Nt_c_o
      ! use module_mp_thompson, only : re_qc_min, re_qc_max, re_qi_min, re_qi_max, re_qs_min, re_qs_max

      ! use module_mp_thompson_make_number_concentrations, only: make_IceNumber, make_DropletNumber, make_RainNumber

      implicit none

      public :: scm_mp_wsm6_pre_init

      private

      logical :: is_initialized = .False.

      ! integer, parameter :: ext_ndiag3d = 37

   contains

!> This subroutine is a wrapper around the actual mp_wsm6_init().
!! \section arg_table_scm_mp_wsm6_init Argument Table
!! \htmlinclude scm_mp_wsm6_init.html
!!
     subroutine scm_mp_wsm6_pre_init(den0, denr, dens, cl, &
                                 cpv, hail_opt, errmsg, errflg)

         implicit none
         ! FOO TODO: ADD COMMENTS FOR ALL VARIABLES

         ! Input arguments:
         real(kind=kind_phys), intent(in   ) :: den0
         real(kind=kind_phys), intent(in   ) :: denr
         real(kind=kind_phys), intent(in   ) :: dens
         real(kind=kind_phys), intent(in   ) :: cl
         real(kind=kind_phys), intent(in   ) :: cpv
         integer,              intent(in   ) :: hail_opt
         
         ! CCPP error handling
         character(len=*),     intent(  out) :: errmsg
         integer,              intent(  out) :: errflg

         ! ! Interface variables
         ! integer,                   intent(in   ) :: ncol
         ! integer,                   intent(in   ) :: nlev
         ! ! Hydrometeors
         ! logical,                   intent(in   ) :: convert_dry_rho
         ! real(kind_phys),           intent(inout) :: spechum(:,:)
         ! real(kind_phys),           intent(inout) :: qc(:,:)
         ! ! Aerosols
         ! logical,                   intent(in   ) :: is_aerosol_aware
         ! logical,                   intent(in   ) :: merra2_aerosol_aware
         ! real(kind_phys),           intent(inout) :: nc(:,:)

         ! real(kind_phys),           intent(in)    :: aerfld(:,:,:)
         ! ! State variables
         ! real(kind_phys),           intent(in   ) :: phil(:,:)
         ! real(kind_phys),           intent(in   ) :: area(:)
         ! ! MPI information
         ! integer,                   intent(in   ) :: mpicomm
         ! integer,                   intent(in   ) :: mpiroot
         ! ! Threading/blocking information
         ! integer,                   intent(in   ) :: threads
         ! ! Extended diagnostics
         ! logical,                   intent(in   ) :: ext_diag
         ! real(kind_phys),           intent(in   ) :: diag3d(:,:,:)
         ! ! CCPP error handling
         ! character(len=*),          intent(  out) :: errmsg
         ! integer,                   intent(  out) :: errflg

         !

         integer :: i, k

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return

         ! ! Consistency checks
         ! if (imp_physics/=imp_physics_thompson) then
         !    write(errmsg,'(*(a))') "Logic error: namelist choice of microphysics is different from Thompson MP"
         !    errflg = 1
         !    return
         ! end if

         ! if (ext_diag) then
         !    if (size(diag3d,dim=3) /= ext_ndiag3d) then
         !       write(errmsg,'(*(a))') "Logic error: number of diagnostic 3d arrays from model does not match requirements"
         !       errflg = 1
         !       return
         !    end if
         ! end if

         ! if (is_aerosol_aware .and. merra2_aerosol_aware) then
         !    write(errmsg,'(*(a))') "Logic error: Only one Thompson aerosol option can be true, either is_aerosol_aware or merra2_aerosol_aware)"
         !    errflg = 1
         !    return
         ! end if

         ! Call Thompson init
         call mp_wsm6_init(den0, denr, dens, cl, &
                           cpv, hail_opt, errmsg, errflg)

         if (errflg /= 0) return

         ! ! For restart runs, the init is done here
         ! if (restart) then
         !   is_initialized = .true.
         !   return
         ! end if

         ! ! Geopotential height in m2 s-2 to height in m
         ! hgt = phil/con_g

        !  ! Ensure non-negative mass mixing ratios of all water variables
        !  where(spechum<0) spechum = 1.0E-10     ! COMMENT, gthompsn, spechum should *never* be identically zero.
        !  where(qc<0)      qc = 0.0
        !  where(qr<0)      qr = 0.0
        !  where(qi<0)      qi = 0.0
        !  where(qs<0)      qs = 0.0
        !  where(qg<0)      qg = 0.0

        !  !> - Convert specific humidity to water vapor mixing ratio.
        !  !> - Also, hydrometeor variables are mass or number mixing ratio
        !  !> - either kg of species per kg of dry air, or per kg of (dry + vapor).
        !  if (merra2_aerosol_aware) then
        !    call get_niwfa(aerfld, nifa, nwfa, ncol, nlev)
        !  end if


        !  qv = spechum/(1.0_kind_phys-spechum)

        !  if (convert_dry_rho) then
        !    qc = qc/(1.0_kind_phys-spechum)
        !    qr = qr/(1.0_kind_phys-spechum)
        !    qi = qi/(1.0_kind_phys-spechum)
        !    qs = qs/(1.0_kind_phys-spechum)
        !    qg = qg/(1.0_kind_phys-spechum)

        !    ni = ni/(1.0_kind_phys-spechum)
        !    nr = nr/(1.0_kind_phys-spechum)
        !    if (is_aerosol_aware .or. merra2_aerosol_aware) then
        !       nc = nc/(1.0_kind_phys-spechum)
        !       nwfa = nwfa/(1.0_kind_phys-spechum)
        !       nifa = nifa/(1.0_kind_phys-spechum)
        !    end if
        !  end if

        !  ! Density of moist air in kg m-3 and inverse density of air
        !  rho = con_eps*prsl/(con_rd*tgrs*(qv+con_eps))
        !  orho = 1.0/rho

        !  ! Ensure we have 1st guess ice number where mass non-zero but no number.
        !  where(qi .LE. 0.0) ni=0.0
        !  where(qi .GT. 0 .and. ni .LE. 0.0) ni = make_IceNumber(qi*rho, tgrs) * orho
        !  where(qi .EQ. 0.0 .and. ni .GT. 0.0) ni=0.0

        !  ! Ensure we have 1st guess rain number where mass non-zero but no number.
        !  where(qr .LE. 0.0) nr=0.0
        !  where(qr .GT. 0 .and. nr .LE. 0.0) nr = make_RainNumber(qr*rho, tgrs) * orho
        !  where(qr .EQ. 0.0 .and. nr .GT. 0.0) nr=0.0


        !  !..Check for existing aerosol data, both CCN and IN aerosols.  If missing
        !  !.. fill in just a basic vertical profile, somewhat boundary-layer following.
        !  if (is_aerosol_aware) then

        !    ! Potential cloud condensation nuclei (CCN)
        !    if (MAXVAL(nwfa) .lt. eps) then
        !      if (mpirank==mpiroot) write(*,*) ' Apparently there are no initial CCN aerosols.'
        !      do i = 1, ncol
        !        if (hgt(i,1).le.1000.0) then
        !          h_01 = 0.8
        !        elseif (hgt(i,1).ge.2500.0) then
        !          h_01 = 0.01
        !        else
        !          h_01 = 0.8*cos(hgt(i,1)*0.001 - 1.0)
        !        endif
        !        niCCN3 = -1.0*ALOG(naCCN1/naCCN0)/h_01
        !        nwfa(i,1) = naCCN1+naCCN0*exp(-((hgt(i,2)-hgt(i,1))/1000.)*niCCN3)
        !        z1 = hgt(i,2)-hgt(i,1)
        !        nwfa2d(i) = nwfa(i,1) * 0.000196 * (50./z1)
        !        do k = 2, nlev
        !          nwfa(i,k) = naCCN1+naCCN0*exp(-((hgt(i,k)-hgt(i,1))/1000.)*niCCN3)
        !        enddo
        !      enddo
        !    else
        !      if (mpirank==mpiroot) write(*,*) ' Apparently initial CCN aerosols are present.'
        !      if (MAXVAL(nwfa2d) .lt. eps) then
        !        !+---+-----------------------------------------------------------------+
        !        !..Scale the lowest level aerosol data into an emissions rate.  This is
        !        !.. very far from ideal, but need higher emissions where larger amount
        !        !.. of (climo) existing and lesser emissions where there exists fewer to
        !        !.. begin as a first-order simplistic approach.  Later, proper connection to
        !        !.. emission inventory would be better, but, for now, scale like this:
        !        !.. where: Nwfa=50 per cc, emit 0.875E4 aerosols per second per grid box unit
        !        !..        that was tested as ~(20kmx20kmx50m = 2.E10 m**-3)
        !        !+---+-----------------------------------------------------------------+
        !        if (mpirank==mpiroot) write(*,*) ' Apparently there are no initial CCN aerosol surface emission rates.'
        !        do i = 1, ncol
        !           z1 = hgt(i,2)-hgt(i,1)
        !           nwfa2d(i) = nwfa(i,1) * 0.000196 * (50./z1)
        !        enddo
        !      else
        !         if (mpirank==mpiroot) write(*,*) ' Apparently initial CCN aerosol surface emission rates are present.'
        !      endif
        !    endif

        !    ! Potential ice nuclei (IN)
        !    if (MAXVAL(nifa) .lt. eps) then
        !      if (mpirank==mpiroot) write(*,*) ' Apparently there are no initial IN aerosols.'
        !      do i = 1, ncol
        !        if (hgt(i,1).le.1000.0) then
        !           h_01 = 0.8
        !        elseif (hgt(i,1).ge.2500.0) then
        !           h_01 = 0.01
        !        else
        !           h_01 = 0.8*cos(hgt(i,1)*0.001 - 1.0)
        !        endif
        !        niIN3 = -1.0*ALOG(naIN1/naIN0)/h_01
        !        nifa(i,1) = naIN1+naIN0*exp(-((hgt(i,2)-hgt(i,1))/1000.)*niIN3)
        !        nifa2d(i) = 0.
        !        do k = 2, nlev
        !           nifa(i,k) = naIN1+naIN0*exp(-((hgt(i,k)-hgt(i,1))/1000.)*niIN3)
        !        enddo
        !      enddo
        !    else
        !      if (mpirank==mpiroot) write(*,*) ' Apparently initial IN aerosols are present.'
        !      if (MAXVAL(nifa2d) .lt. eps) then
        !        if (mpirank==mpiroot) write(*,*) ' Apparently there are no initial IN aerosol surface emission rates, set to zero.'
        !        ! calculate IN surface flux here, right now just set to zero
        !        nifa2d = 0.
        !      else
        !        if (mpirank==mpiroot) write(*,*) ' Apparently initial IN aerosol surface emission rates are present.'
        !      endif
        !    endif

        !    ! Ensure we have 1st guess cloud droplet number where mass non-zero but no number.
        !    where(qc .LE. 0.0) nc=0.0
        !    where(qc .GT. 0 .and. nc .LE. 0.0) nc = make_DropletNumber(qc*rho, nwfa*rho) * orho
        !    where(qc .EQ. 0.0 .and. nc .GT. 0.0) nc = 0.0

        !    ! Ensure non-negative aerosol number concentrations.
        !    where(nwfa .LE. 0.0) nwfa = 1.1E6
        !    where(nifa .LE. 0.0) nifa = naIN1*0.01

        !    ! Copy to local array for calculating cloud effective radii below
        !    nc_local = nc

        ! else if (merra2_aerosol_aware) then

        !    ! Ensure we have 1st guess cloud droplet number where mass non-zero but no number.
        !    where(qc .LE. 0.0) nc=0.0
        !    where(qc .GT. 0 .and. nc .LE. 0.0) nc = make_DropletNumber(qc*rho, nwfa*rho) * orho
        !    where(qc .EQ. 0.0 .and. nc .GT. 0.0) nc = 0.0

        !  else

        !    ! Constant droplet concentration for single moment cloud water as in
        !    ! module_mp_thompson.F90, only needed for effective radii calculation
        !    nc_local = Nt_c_l/rho

        !  end if

        !  if (convert_dry_rho) then
        !    !qc = qc/(1.0_kind_phys+qv)
        !    !qr = qr/(1.0_kind_phys+qv)
        !    !qi = qi/(1.0_kind_phys+qv)
        !    !qs = qs/(1.0_kind_phys+qv)
        !    !qg = qg/(1.0_kind_phys+qv)

        !    ni = ni/(1.0_kind_phys+qv)
        !    nr = nr/(1.0_kind_phys+qv)
        !    if (is_aerosol_aware .or. merra2_aerosol_aware) then
        !       nc = nc/(1.0_kind_phys+qv)
        !       nwfa = nwfa/(1.0_kind_phys+qv)
        !       nifa = nifa/(1.0_kind_phys+qv)
        !    end if
        !  end if

         is_initialized = .true.

      end subroutine scm_mp_wsm6_pre_init
end module scm_mp_wsm6_pre
