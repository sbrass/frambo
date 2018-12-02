!> \file c_rambo.f90
!! \brief Provide a C-binding of the RAMBO type and type-bound procedures.
module c_rambo
  use constants
  use iso_c_binding
  use phs_rambo

  implicit none

contains
  function declare_phs_rambo (n_particles, total_momentum, masses) result (ptr) bind(C)
    integer(c_int), value, intent(in) :: n_particles
    real(c_double), dimension(4), intent(in) :: total_momentum
    real(c_double), dimension(*), intent(in) :: masses
    type(c_ptr) :: ptr
    ! type(phs_rambo_t), target :: phs
    type(phs_rambo_t), pointer :: phs
    allocate (phs, source=phs_rambo_t (n_particles, total_momentum, masses(:n_particles)))
    ptr = c_loc (phs)
  end function declare_phs_rambo

  subroutine free_phs_rambo (ptr) bind (C)
    type(c_ptr), intent(in), value :: ptr
    type(phs_rambo_t), pointer :: phs
    call c_f_pointer (ptr, phs)
    if (associated (phs)) then
       deallocate (phs)
    end if
  end subroutine free_phs_rambo

  subroutine write_phs_rambo (ptr) bind (C)
    type(c_ptr), intent(in), value :: ptr
    type(phs_rambo_t), pointer :: phs
    call c_f_pointer (ptr, phs)
    call phs%write ()
  end subroutine write_phs_rambo

  !> \warning We do not check that r_in has the correct size.
  subroutine generate_phs_rambo (ptr, n_r_in, r_in) bind (C)
    type(c_ptr), intent(in), value :: ptr
    integer(c_int), intent(in), value :: n_r_in
    real(c_double), dimension(*), intent(in) :: r_in
    type(phs_rambo_t), pointer :: phs
    call c_f_pointer (ptr, phs)
    call phs%generate (r_in(:n_r_in))
  end subroutine generate_phs_rambo

  function get_weight_phs_rambo (ptr) result (w) bind (C)
    type(c_ptr), intent(in), value :: ptr
    real(c_double) :: w
    type(phs_rambo_t), pointer :: phs
    call c_f_pointer (ptr, phs)
    w = phs%weight ()
  end function get_weight_phs_rambo


end module c_rambo
