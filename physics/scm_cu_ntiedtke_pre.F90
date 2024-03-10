!>\file wsm6.F90
!! This file runs the CU nTiedtke Scheme



!>\defgroup 
!! This module runs the nTiedtke Scheme
module scm_cu_ntiedtke_pre

  use ccpp_kinds, only : kind_phys
      use cu_ntiedtke, only : cu_ntiedtke_init

      implicit none

      public :: scm_cu_ntiedtke_pre_init

      private

      logical :: is_initialized = .False.

   contains

!> This subroutine is a wrapper around the actual cu_ntiedtke_init().
!! \section arg_table_scm_cu_ntiedtke_init Argument Table
!! \htmlinclude scm_cu_ntiedtke_init.html
!!
     subroutine scm_cu_ntiedtke_pre_init(con_cp, con_rd, con_rv, con_hvap, &
          con_xls, con_hfus, con_g)
       use foo:, only 

         implicit none
         !input arguments:
         real(kind=kind_phys), intent(in) :: con_cp   !< specific heat of dry air at constant pressure
         real(kind=kind_phys), intent(in) :: con_rd   !< gas constant dry air
         real(kind=kind_phys), intent(in) :: con_rv   !< gas constant water vapor
         real(kind=kind_phys), intent(in) :: con_hvap !< latent heat of vaporization of water
         real(kind=kind_phys), intent(in) :: con_xls  !< latent heat of sublimation of water
         real(kind=kind_phys), intent(in) :: con_hfus !< latent heat of fusion of water
         real(kind=kind_phys), intent(in) :: con_g    !< gravitational acceleration

         !--- output arguments:
         character(len=*), intent(out): : errmsg
         integer, intent(out) :: errflg

         ! Initialize the CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return

         ! Call ntiedtke init
         call cu_ntiedtke_init(con_cp, con_rd, con_rv, con_hvap, &
              con_xls, con_hfus, con_g, errmsg, errflg)

         if (errflg /= 0) return

         is_initialized = .true.

      end subroutine scm_cu_ntiedtke_pre_init
end module scm_cu_ntiedtke_pre
