!>\file wsm6.F90
!! This file runs the Single-Moment 6-class Microphysics scheme (WSM6)


!>\defgroup aathompson Aerosol-Aware Thompson MP Module
!! This module runs the Single-Moment 6-class Microphysics scheme (WSM6)
module scm_mp_wsm6_post

      use ccpp_kinds, only : kind_phys

      use mp_wsm6, only : mp_wsm6_finalize
      ! use module_mp_thompson, only : naIN0, naIN1, naCCN0, naCCN1, eps, Nt_c_l, Nt_c_o
      ! use module_mp_thompson, only : re_qc_min, re_qc_max, re_qi_min, re_qi_max, re_qs_min, re_qs_max

      ! use module_mp_thompson_make_number_concentrations, only: make_IceNumber, make_DropletNumber, make_RainNumber

      implicit none

      public :: scm_mp_wsm6_post_finalize

      private

      logical :: is_initialized = .False.

      ! integer, parameter :: ext_ndiag3d = 37

   contains

!> \section arg_table_scm_mp_wsm6_finalize Argument Table
!! \htmlinclude scm_mp_wsm6_finalize.html
!!
      subroutine scm_mp_wsm6_post_finalize(errmsg, errflg)

         implicit none

         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0



         if (.not. is_initialized) return

         call mp_wsm6_finalize(errmsg, errflg)

         is_initialized = .false.

      end subroutine scm_mp_wsm6_post_finalize

end module scm_mp_wsm6_post
