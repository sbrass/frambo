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
     procedure :: invert => phs_rambo_invert
     procedure, private :: generate_intermediates => phs_rambo_generate_intermediates
     procedure, private :: invert_intermediates => phs_rambo_invert_intermediates
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
         if (phi > PI) phi = phi - 2. * PI
         q_abs = 4. * M(i - 1) * rho (M(i - 1), M(i), phs%p_m(i - 1))
         p = q_abs * [cos(phi) * sqrt(1. - cos_theta**2), &
              sin(phi) * sqrt(1. - cos_theta**2), &
              cos_theta]
         phs%p(i - 1) = momentum4_t ([sqrt (q_abs**2 + phs%p_m(i - 1)**2), p])
         QNext = momentum4_t ([sqrt (q_abs**2 + M(i)**2), -p])
         call phs%p(i - 1)%boost (Q)
         call QNext%boost (Q)
         ! print *, "================================================================================"
         ! call Q%write ()
         ! call phs%p(i - 1)%write ()
         ! call QNext%write()
         ! print *, "================================================================================"
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
      M(n) = phs%p_m(n)
      ! Prepare K(1) from which the M's are inferred.
      call calculate_k (r)
      do i = 1, n - 1
         M(i) = K(i)
         do j = i, n
            M(i) = M(i) + phs%p(j)%get_mass ()
         end do
      end do
      ! TODO Overall weight factor missing
      phs%weight = (1. / (2. * PI)**(3 * n - 4)) / 8.**(n - 1)
      phs%weight = rho(M(n - 1), phs%p(n)%get_mass (), phs%p(n - 1)%get_mass ())
      do i = 2, n - 1
         phs%weight = phs%weight * &
              rho(M(i - 1), M(i), phs%p(i - 1)%get_mass ()) / &
              rho(K(i - 1), K(i), 0._default) * &
              M(i) / K(i)
      end do
      phs%weight = phs%weight * (K(1) / M(1))**(2 * n - 4)
    end associate
  contains
    ! Mi**2 = u_2 * ... * u_i M1**2
    subroutine calculate_k (r)
      real(default), dimension(phs%n_particles - 2), intent(in) :: r
      real(default), dimension(:), allocatable :: u
      integer :: i
      associate (n => phs%n_particles, K => phs%K, M => phs%M)
        K = 0; K(1) = M(1)
        do i = 1, n - 1
           K(1) = K(1) - phs%p(i)%get_mass()
        end do
        allocate (u(2:n - 1), source=0._default)
        call solve_for_u (r, u)
        do i = 2, n - 1
           K(i) = sqrt (u(i) * K(i - 1)**2)
        end do
      end associate
    end subroutine calculate_k

    subroutine solve_for_u (r, u)
      real(default), dimension(phs%n_particles - 2), intent(in) :: r
      real(default), dimension(2:phs%n_particles - 1), intent(out) :: u
      integer :: i, j
      real(default) :: f, f_deriv
      associate (n => phs%n_particles)
        do i = 2, n - 1
           ! Newton's method
           u(i) = r(i - 1)
           do j = 1, 100
              f = (n + 1 - i) * u(i)**(n - i) &
                   - (n - i) * u(i)**(n + 1 - i) - r(i - 1)
              f_deriv = (n + 1 - i) * (n - i) * u(i)**(n - 1 - i) &
                   - (n - i) * (n + 1 - i) * u(i)**(n - i)
              ! TODO epsilon parameter
              u(i) = u(i) - f / f_deriv
              if (abs (f) < 1e-15) exit
           end do
        end do
      end associate
    end subroutine solve_for_u
  end subroutine phs_rambo_generate_intermediates

  subroutine phs_rambo_invert (phs, r)
    class(phs_rambo_t), intent(inout) :: phs
    real(default), dimension(3 * phs%n_particles - 4), intent(out) :: r
    type(momentum4_t), dimension(:), allocatable :: Q
    type(momentum4_t) :: p
    real(default) :: cos_theta, phi
    integer :: i, j
    associate (n => phs%n_particles, M => phs%M)
      allocate (Q(n), source=zero_momentum)
      M(1) = phs%total_p%get_mass ()
      Q(1) = momentum4_t (M(1))
      Q(n) = phs%p(n)
      do i = 2, n - 1
         do j = i, n
            Q(i) = Q(i) + phs%p(j)
         end do
         call Q(i)%update_mass ()
         M(i) = Q(i)%get_mass ()
      end do
      call phs%invert_intermediates (r(1:n - 2))
      do i = 2, n
         p = phs%p(i - 1)
         call p%boost (Q(i - 1), invert = .true.)
         cos_theta = p%cos_theta ()
         phi = p%phi (); if (phi < 0.) phi = phi + 2. * PI
         ! TODO understand atan2 phi
         r(n - 5 + 2 * i) = (cos_theta + 1.) / 2.
         r(n - 4 + 2 * i) = phi / (2. * PI)
      end do
    end associate
  end subroutine phs_rambo_invert

  subroutine phs_rambo_invert_intermediates (phs, r)
    class(phs_rambo_t), intent(inout) :: phs
    real(default), dimension(phs%n_particles - 2), intent(out) :: r
    integer :: i, j
    associate (n => phs%n_particles, K => phs%K, M => phs%M)
      K = M
      do i = 1, n - 1
         do j = i, n
            K(i) = K(i) - phs%p(j)%get_mass ()
         end do
      end do
      call solve_for_r (r)
      ! Calculate weight (for check)
    end associate
  contains
    subroutine solve_for_r (r)
      real(default), dimension(phs%n_particles - 2), intent(out) :: r
      integer :: i
      real(default) :: u
      associate (n => phs%n_particles)
        do i = 2, n - 1
           u = (phs%K(i) / phs%K(i - 1))**2
           r(i - 1) = (n + 1 - i) * u**(n - i) &
                - (n - i) * u**(n + 1 - i)
        end do
      end associate
    end subroutine solve_for_r
  end subroutine phs_rambo_invert_intermediates

  subroutine phs_rambo_write (phs, unit)
    class(phs_rambo_t), intent(in) :: phs
    integer, intent(in), optional :: unit
    integer :: u, i
    u = OUTPUT_UNIT; if (present (unit)) u = unit
    write (u, "(80(A))") "================================================================================"
    write (u, "(A)") "Phasespace RAMBO type:"
    write (u, "(A,F12.4)") "Weight: ", phs%weight
    write (u, "(11X,A,11X,A,1X,A)") "K", "M"
    do i = 1, phs%n_particles
       write (u, "(2(F12.4))") phs%K(i), phs%M(i)
    end do
    call phs%base_write (unit)
  end subroutine phs_rambo_write
end module phs_rambo
