module phs
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use constants
  use momentum

  implicit none

  type, abstract :: phs_t
     integer :: n_particles = 0
     type(momentum4_t) :: total_p = zero_momentum
     type(momentum4_t), dimension(:), allocatable :: p
     real(default), dimension(:), allocatable :: p_m
     real(default) :: flux = 0
     real(default) :: volume = 0
     real(default) :: jacobian = 0
   contains
     procedure :: init => phs_init
     procedure(phs_write), deferred :: write
     procedure :: base_write => phs_base_write
     procedure :: get_p => phs_get_p
     procedure(phs_generate), deferred :: generate
     procedure :: get_weight => phs_get_weight
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
    allocate (phs%p_m(n_particles), source=mass)
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
    write (u, "(80(A))") "================================================================================"
    write (u, "(A)") "Phasespace basic type:"
    write (u, "(A,I0)") "Number of particles: ", phs%n_particles
    write (u, "(A,ES12.5)") "Q**2: ", squared (phs%total_p)
    write (u, "(A,ES12.5)") "Volume: ", phs%volume
    write (u, "(A,ES12.5)") "Jacobian: ", phs%jacobian
    call phs%total_p%write (unit)
    total_p = phs%total_p
    do i = 1, phs%n_particles
       ! write (unit, "(A,1X,F12.4)") "On-shell mass:", phs%p_m (i)
       call phs%p(i)%write (unit)
       total_p = total_p - phs%p(i)
    end do
    write (u, "(A,1X,74(A))") "TOTAL", "=========================================================================="
    call total_p%write (unit)
  end subroutine phs_base_write

  function phs_get_p (phs, i) result(p)
    class(phs_t), intent(in) :: phs
    integer, intent(in) :: i
    real(default), dimension(4) :: p
    p = phs%p(i)%get_p ()
  end function phs_get_p

  real(default) function phs_get_weight (phs) result (w)
    class(phs_t), intent(in) :: phs
    w = phs%jacobian * phs%volume
  end function phs_get_weight

  logical function phs_is_valid (phs) result (flag)
    class(phs_t), intent(in) :: phs
    type(momentum4_t) :: total_p
    integer :: i
    total_p = phs%total_p
    do i = 1, phs%n_particles
       total_p = total_p - phs%p(i)
       call total_p%update_mass ()
    end do
    flag = all (phs%p%is_on_shell ())
    flag = flag .and. total_p%get_mass () < 1e5
  end function phs_is_valid
end module phs
