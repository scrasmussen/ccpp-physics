!> \file mmm_kinds.F90
!! File contains ccpp_kinds which is used by the MMM Physics scheme.
! The filename and module are mmm_kinds to clearly indicate which scheme
! uses this file.
module ccpp_kinds
  use machine,  only: kind_phys, kind_dbl_prec, kind_sngl_prec
  implicit none
end module ccpp_kinds

module mmm_kinds
  use ccpp_kinds
end module mmm_kinds
