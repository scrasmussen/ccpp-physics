!> \file scm_kinds.f90
!!  Wrapper with generic interfaces for w3emc library
module w3emc_wrapper
  use iso_fortran_env, only: real32, real64
  implicit none

  interface w3difdat_w
     module procedure :: w3difdat32
     module procedure :: w3difdat64
  end interface w3difdat_w

  interface w3movdat_w
     module procedure :: w3movdat32
     module procedure :: w3movdat64
  end interface w3movdat_w

contains
  subroutine w3difdat32(jdat, idat, it, rinc)
    integer :: jdat(8), idat(8), it
    real(real32) :: rinc(5)
    call w3difdat(jdat, idat, it, rinc)
  end subroutine w3difdat32

  subroutine w3difdat64(jdat, idat, it, rinc)
    integer :: jdat(8), idat(8), it
    real(real64) :: rinc(5)
    call w3difdat(jdat, idat, it, rinc)
  end subroutine w3difdat64

  subroutine w3movdat32(rinc, idat, jdat)
    real(real32) :: rinc(5)
    integer :: idat(8), jdat(8)
    call w3movdat(rinc, idat, jdat)
  end subroutine w3movdat32

  subroutine w3movdat64(rinc, idat, jdat)
    real(real64) :: rinc(5)
    integer :: idat(8), jdat(8)
    call w3movdat(rinc, idat, jdat)
  end subroutine w3movdat64


end module w3emc_wrapper
