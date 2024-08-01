!>\file wsm6.F90
!! This file runs the nTiedtke scheme


!>\defgroup TODO: ADD THIS LINE
!! This module runs the nTiedtke scheme
module mmm_cu_ntiedtke_post

      use ccpp_kinds, only : kind_phys

      use mmm_cu_ntiedtke, only : mmm_cu_ntiedtke_finalize

      implicit none

      public :: mmm_cu_ntiedtke_post_finalize

      private

      logical :: is_initialized = .False.

   contains

!> \section arg_table_mmm_cu_ntiedtke_finalize Argument Table
!! \htmlinclude mmm_cu_ntiedtke_finalize.html
!!
      subroutine mmm_cu_ntiedtke_post_finalize(errmsg, errflg)

         implicit none

         character(len=*),          intent(  out) :: errmsg
         integer,                   intent(  out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0



         if (.not. is_initialized) return

         call mmm_cu_ntiedtke_finalize(errmsg, errflg)

         is_initialized = .false.

      end subroutine mmm_cu_ntiedtke_post_finalize

end module mmm_cu_ntiedtke_post
