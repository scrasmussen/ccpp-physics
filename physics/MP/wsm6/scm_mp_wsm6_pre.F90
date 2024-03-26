!>\file scm_mp_wsm6_pre.F90
!! This file runs the Single-Moment 6-class Microphysics scheme (WSM6)


!>\defgroup aathompson Aerosol-Aware Thompson MP Module
!! This module runs the Single-Moment 6-class Microphysics scheme (WSM6)
module scm_mp_wsm6_pre

  use ccpp_kinds, only : kind_phys
      use mp_wsm6, only : mp_wsm6_init

      implicit none

      public :: scm_mp_wsm6_pre_init

      private

      logical :: is_initialized = .False.

   contains

!> This subroutine is a wrapper around the actual mp_wsm6_init().
!! \section arg_table_scm_mp_wsm6_init Argument Table
!! \htmlinclude scm_mp_wsm6_init.html
!!
     subroutine scm_mp_wsm6_pre_init(den0, denr, dens, cl, &
                                 cpv, hail_opt, errmsg, errflg)

         implicit none
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

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return

         ! Call wsm6 init
         call mp_wsm6_init(den0, denr, dens, cl, &
                           cpv, hail_opt, errmsg, errflg)

         if (errflg /= 0) return

         is_initialized = .true.

      end subroutine scm_mp_wsm6_pre_init
end module scm_mp_wsm6_pre
