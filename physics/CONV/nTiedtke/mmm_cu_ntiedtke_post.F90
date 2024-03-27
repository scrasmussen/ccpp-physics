!>\file wsm6.F90
!! This file runs the nTietke scheme


!>\defgroup TODO: ADD THIS LINE
!! This module runs the nTietke scheme
module scm_cu_ntiedtke_post

      use ccpp_kinds, only : kind_phys

      use cu_ntiedtke, only : cu_ntiedtke_finalize

      implicit none

      public :: scm_cu_ntiedtke_post_finalize

      private

      logical :: is_initialized = .False.

   contains

!> \section arg_table_scm_cu_ntiedtke_finalize Argument Table
!! \htmlinclude scm_cu_ntiedtke_finalize.html
!!
      subroutine scm_cu_ntiedtke_post_finalize(errmsg, errflg)

         implicit none

         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0



         if (.not. is_initialized) return

         call cu_ntiedtke_finalize(errmsg, errflg)

         is_initialized = .false.

      end subroutine scm_cu_ntiedtke_post_finalize

end module scm_cu_ntiedtke_post
