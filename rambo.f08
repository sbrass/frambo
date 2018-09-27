module phs_rambo
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use constants
  use momentum
  use phs

  implicit none

  type, extends(phs_t) :: phs_rambo_t
     private
     real(default) :: weight = 0._default
     real(default), dimension(:), allocatable :: K, M
   contains
     procedure :: write => phs_rambo_write
     procedure :: generate => phs_rambo_generate
     procedure, private :: generate_intermediates => phs_rambo_generate_intermediates
  end type phs_rambo_t

  interface phs_rambo_t
     module procedure :: phs_rambo_init
  end interface phs_rambo_t
contains
  type(phs_rambo_t) function phs_rambo_init (n_particles, total_momentum, masses) result (phs)
    integer, intent(in) :: n_particles
    real(default), dimension(4) :: total_momentum
    real(default), dimension(:), intent(in) :: masses
    call phs%init (n_particles, total_momentum, masses)
    allocate (phs%K(phs%n_particles), source=0._default)
    allocate (phs%M(phs%n_particles), source=0._default)
  end function phs_rambo_init

  subroutine phs_rambo_generate (phs, r)
    class(phs_rambo_t), intent(inout) :: phs
    real(default), dimension(3 * phs%n_particles - 4), intent(in) :: r
    real(default) :: cos_theta, phi, q_abs
    real(default), dimension(3) :: p
    type(momentum4_t) :: Q, QNext
    integer :: i
    associate (n => phs%n_particles, M => phs%M)
      call phs%generate_intermediates (r(1:n - 2))
      Q = phs%total_p
      do i = 2, n
         cos_theta = 2. * r(n - 5 + 2 * i) - 1.
         phi = 2. * PI * r(n - 4 + 2 * i)
         if (phi > PI) phi = phi - PI
         q_abs = 4. * M(i - 1) * rho (M(i - 1), M(i), phs%p(i - 1)%get_mass ())
         p = [cos(phi) * sqrt(1 - cos_theta**2), &
              sin(phi) * sqrt(1 - cos_theta**2), &
              cos_theta]
         p = p * q_abs
         call phs%p(i - 1)%set_p (p)
         call phs%p(i - 1)%boost(Q)
         QNext = Q - phs%p(i - 1)
         call QNext%set_mass (M(i))
         print *, mass(QNext), QNext%get_mass ()
         call QNext%rescale_energy ()
         Q = QNext
      end do
      phs%p(n) = Q
    end associate
  end subroutine phs_rambo_generate

  subroutine phs_rambo_generate_intermediates (phs, r)
    class(phs_rambo_t), intent(inout) :: phs
    real(default), dimension(phs%n_particles - 2), intent(in) :: r
    integer :: i, j
    associate (n => phs%n_particles, K => phs%K, M => phs%M)
      M(1) = phs%total_p%get_mass ()
      M(n) = phs%p(n)%get_mass ()
      ! Prepare K(1) from which the M's are inferred.
      call calculate_k (r)
      do i = 1, n - 1
         M(i) = K(i)
         do j = i, n
            M(i) = M(i) + phs%p(j)%get_mass ()
         end do
      end do
      ! TODO Overall weight factor missing
      phs%weight = 8. * rho(M(n - 1), phs%p(n)%get_mass (), phs%p(n - 1)%get_mass ())
      do i = 2, n - 1
         phs%weight = phs%weight * &
              rho(M(i - 1), M(i), phs%p(i - 1)%get_mass ()) / &
              rho(K(i - 1), K(i), 0._default) * &
              M(i) / K(i)
      end do
      phs%weight = phs%weight * ((1) / M(1))**(2 * n - 4)
    end associate
  contains
    subroutine calculate_k (r)
      real(default), dimension(phs%n_particles - 2), intent(in) :: r
      real(default), dimension(:), allocatable :: u
      integer :: i
      associate (n => phs%n_particles, K => phs%K, M => phs%M)
        K(1) = M(1)
        do i = 1, n - 1
           K(1) = K(1) - phs%p(i)%get_mass()
        end do
        allocate (u(2:n - 1), source=0._default)
        call solve_for_u (r, u)
        do i = 2, n - 1
           K(i) = u(i) * K(i - 1)
        end do
      end associate
    end subroutine calculate_k

    subroutine solve_for_u (r, u)
      real(default), dimension(phs%n_particles - 2), intent(in) :: r
      real(default), dimension(2:phs%n_particles - 1), intent(out) :: u
      integer :: i, j
      real(default) :: f, f_deriv
      ! TODO start value != zero, choose a better one
      u = 0.5_default
      associate (n => phs%n_particles)
        do i = 2, n - 1
           ! Newton's method
           ! TODO parameter for loop number
           do j = 1, 1000
              f = (n + 1 - i) * u(i)**(n - i) &
                   - (n - i) * u(i)**(n + 1 - i) - r(i - 1)
              f_deriv = (n + 1 - i) * (n - i) * u(i)**(n - 1 - i) &
                   - (n - i) * (n + 1 - i) * u(i)**(n - i)
              ! TODO parameter for tolerance
              ! Guard against zero division
              if (f < 1e-23 * f_deriv) exit
              u(i) = u(i) - f / f_deriv
           end do
        end do
      end associate
    end subroutine solve_for_u


  end subroutine phs_rambo_generate_intermediates

  subroutine phs_rambo_write (phs, unit)
    class(phs_rambo_t), intent(in) :: phs
    integer, intent(in), optional :: unit
    integer :: u, i
    u = OUTPUT_UNIT; if (present (unit)) u = unit
    write (u, "(A)") "Phasespace RAMBO type:"
    write (u, "(A,F12.4)") "Weight: ", phs%weight
    write (u, "(11X,A,11X,A,1X,A)") "K", "M"
    do i = 1, phs%n_particles
       write (u, "(2(F12.4))") phs%K(i), phs%M(i)
    end do
    call phs%base_write (unit)
  end subroutine phs_rambo_write
end module phs_rambo
