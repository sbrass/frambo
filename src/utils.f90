!> \file utils.f90
!! \brief Supply numeric utils.

!> Module: numerics
!> \author: Simon Braß
!> \date: 9.11.2018
module numerics
  use constants

  implicit none

  private

  public faculty
contains
  !> \brief Compute faculty \f$n!\f$ of positive integer \f$n\f$.
  !! \return f
  function faculty (n) result (f)
    integer, intent(in) :: n
    real(default) :: f
    integer :: i
    f = 1
    do i = n, 1, -1
       f = f * i
    end do
  end function faculty
end module numerics
