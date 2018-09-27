module phs
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use constants
  use momentum

  implicit none

  type, abstract :: phs_t
     integer :: n_particles = 0
     type(momentum4_t) :: total_p = zero_momentum
     type(momentum4_t), dimension(:), allocatable :: p
     real(default) :: jacobian = 0
   contains
     procedure :: init => phs_init
     procedure(phs_write), deferred :: write
     procedure :: base_write => phs_base_write
     procedure :: get_p => phs_get_p
     procedure(phs_generate), deferred :: generate
     procedure :: get_jacobian => phs_get_jacobian
     procedure :: is_valid => phs_is_valid
  end type phs_t

  abstract interface
     subroutine phs_generate (phs, r)
       import
       class(phs_t), intent(inout) :: phs
       real(default), dimension(3 * phs%n_particles - 4), intent(in) :: r
     end subroutine phs_generate

     subroutine phs_write (phs, unit)
       import
       class(phs_t), intent(in) :: phs
       integer, intent(in), optional :: unit
     end subroutine phs_write
  end interface

contains
    subroutine phs_init (phs, n_particles, total_momentum, mass)
    class(phs_t), intent(inout) :: phs
    integer, intent(in) :: n_particles
    real(default), dimension(4) :: total_momentum
    real(default), dimension(n_particles) :: mass
    integer :: i
    phs%n_particles = n_particles
    phs%total_p = momentum4_t (total_momentum)
    allocate (phs%p(n_particles), source=zero_momentum)
    do i = 1, n_particles
       phs%p(i) = momentum4_t (mass(i))
    end do
  end subroutine phs_init

  subroutine phs_base_write (phs, unit)
    class(phs_t), intent(in) :: phs
    integer, intent(in), optional :: unit
    integer :: u, i
    type(momentum4_t) :: total_p
    u = OUTPUT_UNIT; if (present (unit)) u = unit
    write (u, "(A)") "Phasespace basic type:"
    write (u, "(A,I0)") "Number of particles: ", phs%n_particles
    write (u, "(A,F12.4)") "Q**2: ", phs%total_p%p_sq
    call phs%total_p%write (unit)
    total_p = phs%total_p
    do i = 1, phs%n_particles
       call phs%p(i)%write (unit)
       total_p = total_p - phs%p(i)
    end do
    call total_p%write ()
  end subroutine phs_base_write

  function phs_get_p (phs, i) result(p)
    class(phs_t), intent(in) :: phs
    integer, intent(in) :: i
    real(default), dimension(4) :: p
    p = phs%p(i)%get_p ()
  end function phs_get_p

  real(default) function phs_get_jacobian (phs) result (j)
    class(phs_t), intent(in) :: phs
    j = phs%jacobian
  end function phs_get_jacobian

  logical function phs_is_valid (phs) result (flag)
    class(phs_t), intent(in) :: phs
    flag = phs%jacobian /= 0
  end function phs_is_valid
end module phs