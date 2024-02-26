! ########################################################################################### 
!
! This controls the working precision within the MMM-physics scheme(s) in the CCPP.
!
! Comments:
! It's common (best) practice for external schemes to parameterize their working precision(s),
! to facilitate portability/interoperability. Then to incorporate these scheme(s) within a 
! host requires a small interface layer that maps the working precision needed by the scheme
! to that provided from the CCPP (in machine.F90).
!
! The name of this interface layer should reflect what is native to the scheme.
!
! *NOTE* This module name (ccpp_kinds) is what is referenced in the MMM-physics submodule,
!        and is a bit too general.
!
!        Two suggestions:
!          - rename ccpp_kinds to mmm_kinds within the mmm_physics repository, and then this
!            module.
!          - use the native types from the MMM MPAS physics. This assumes that the CCPP is
!            using the MPAS authoratative repository. Currently we are using a separate MMM
!            physics repository and not the authoratative.
!        There is also the option of changing the mmm_physics code to import the CCPP working
!        types directly. This is fine for code that is part of the authoratative repository,
!        but not ideal for code that lives in submodules, where the wokring precision names
!        are native to that code. Also, adding "ccpp hooks" within submodules should be
!        avoided and is often met with resistance (rightfully so) from the developer(s).
!
! For example, if using the native MPAS code at it's type names:
!
!!!!MODULE mmm_kinds
!!!!  use machine,  only: kind_phys, kind_dbl_prec, kind_sngl_prec
!!!!  implicit none
!!!!  integer,parameter :: single_precision_for_mpas  = kind_sngl_prec
!!!!  integer,parameter :: double_precision_for_mpas  = kind_dbl_prec
!!!!  integer,parameter :: working_precision_for_mpas = kind_phys
!!!!  ....
!!!!END MODULE mmm_kinds
!
!      
! ###########################################################################################
MODULE ccpp_kinds
  use machine,  only: kind_phys, kind_dbl_prec, kind_sngl_prec
  implicit none
END MODULE ccpp_kinds
