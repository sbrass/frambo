module momentum
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use constants, only: default

  implicit none

  private

  type :: momentum4_t
     private
     ! Mass field
     real(default) :: m = 0
     real(default), dimension(4) :: p = &
          [0._default, 0._default, 0._default, 0._default]
   contains
     procedure :: write => momentum4_write
     procedure :: get_p => momentum4_get_p
     procedure :: get_mass => momentum4_get_mass
     procedure :: get_energy => momentum4_get_energy
     procedure :: get_momentum => momentum4_get_momentum
     procedure :: cos_theta => momentum4_get_cos_theta
     procedure :: phi => momentum4_get_phi
     procedure :: set_mass => momentum4_set_mass
     procedure :: update_mass => momentum4_update_mass
     procedure :: set_p_on_shell => momentum4_set_p_on_shell
     procedure :: is_on_shell => momentum4_is_on_shell
     generic :: operator(+) => add
     generic :: operator(-) => subtract
     procedure, private :: add => momentum4_add
     procedure, private :: subtract => momentum4_subtract
     procedure :: boost => momentum4_boost
  end type momentum4_t

  type(momentum4_t), parameter :: zero_momentum = momentum4_t (0._default, [0._default, 0._default, 0._default, 0._default])

  interface mass
     module procedure :: mass_momentum4, mass_array4
  end interface mass

  interface squared
     module procedure :: squared_momentum4, squared_array4
  end interface squared

  interface momentum4_t
     module procedure :: momentum4_init, momentum4_init_at_rest
  end interface momentum4_t

  public :: momentum4_t, zero_momentum, mass, squared, rho
contains
  real(default) function rho (M1, M2, m)
    real(default), intent(in) :: M1, M2, m
    real(default) :: MP, MM
    ! rho = sqrt ((M1**2 - (M2 + m)**2) * (M1**2 - (M2 - m)**2))
    MP = (M1 - (M2 + m)) * (M1 + (M2 + m))
    MM = (M1 - (M2 - m)) * (M1 + (M2 - m))
    rho = sqrt (MP) * sqrt (MM)
    rho = rho / (8._default * M1**2)
  end function rho

  elemental real(default) function squared_momentum4 (p) result (p_sq)
    type(momentum4_t), intent(in) :: p
    p_sq = squared_array4 (p%p)
  end function squared_momentum4

  pure real(default) function squared_array4 (p) result (p_sq)
    real(default), dimension(4), intent(in) :: p
    p_sq = p(1) * p(1) - dot_product (p(2:), p(2:))
  end function squared_array4

  elemental real(default) function mass_momentum4 (p) result (m)
    type(momentum4_t), intent(in) :: p
    m = mass (p%p)
  end function mass_momentum4

  pure real(default) function mass_array4 (p) result (m)
    real(default), dimension(4), intent(in) :: p
    real(default) :: p_sq
    p_sq = squared (p)
    if (p_sq >= 0) then
       m = sqrt(p_sq)
    else
       m = sqrt (-p_sq)
    end if
  end function mass_array4

  type(momentum4_t) function momentum4_init (v) result (p)
    real(default), dimension(4), intent(in) :: v
    p%p = v
    p%m = mass (v)
  end function momentum4_init

  type(momentum4_t) function momentum4_init_at_rest (m) result (p)
    real(default), intent(in) :: m
    p%p = 0
    p%p(1) = m
    p%m = m
  end function momentum4_init_at_rest

  subroutine momentum4_write (p, unit)
    class(momentum4_t), intent(in) :: p
    integer, intent(in), optional :: unit
    integer :: u
    u = OUTPUT_UNIT; if (present (unit)) u = unit
    write (u, "(12X,A,51X,A)") "M", "P"
    write (u, "(5(1X,F12.4),1X,A,F12.4,A)") p%m, p%p, "(", squared (p) ,")"
  end subroutine momentum4_write

  function momentum4_get_p (p) result (v)
    class(momentum4_t), intent(in) :: p
    real(default), dimension(4) :: v
    v = p%p
  end function momentum4_get_p

  real(default) function momentum4_get_mass (p) result (m)
    class(momentum4_t), intent(in) :: p
    m = p%m
  end function momentum4_get_mass

  real(default) function momentum4_get_energy (p) result (energy)
    class(momentum4_t), intent(in) :: p
    energy = p%p(1)
  end function momentum4_get_energy

  function momentum4_get_momentum (p) result (momentum)
    class(momentum4_t), intent(in) :: p
    real(default), dimension(3) :: momentum
    momentum = p%p(2:4)
  end function momentum4_get_momentum

  function momentum4_get_cos_theta (p) result (ctheta)
    class(momentum4_t), intent(in) :: p
    real(default) :: ctheta
    ctheta = p%p(4) / sqrt (dot_product (p%p(2:4), p%p(2:4)))
  end function momentum4_get_cos_theta

  function momentum4_get_phi (p) result (phi)
    class(momentum4_t), intent(in) :: p
    real(default) :: phi
    phi = atan2 (p%p(3), p%p(2))
  end function momentum4_get_phi

  ! elemental type(momentum4_t) function momentum4_add (p1, p2) result (p)
  type(momentum4_t) function momentum4_add (p1, p2) result (p)
    class(momentum4_t), intent(in) :: p1
    type(momentum4_t), intent(in) :: p2
    p%p = p1%p + p2%p
  end function momentum4_add

  ! elemental type(momentum4_t) function momentum4_subtract (p1, p2) result (p)
  type(momentum4_t) function momentum4_subtract (p1, p2) result (p)
    class(momentum4_t), intent(in) :: p1
    type(momentum4_t), intent(in) :: p2
    p%p = p1%p - p2%p
  end function momentum4_subtract

  subroutine momentum4_set_mass (p, m)
    class(momentum4_t), intent(inout) :: p
    real(default), intent(in) :: m
    p%m = m
  end subroutine momentum4_set_mass

  subroutine momentum4_update_mass (p)
    class(momentum4_t), intent(inout) :: p
    p%m = mass (p)
  end subroutine momentum4_update_mass

  subroutine momentum4_set_p_on_shell (p, v)
    class(momentum4_t), intent(inout) :: p
    real(default), dimension(3), intent(in) :: v
    p%p(2:4) = v
    p%p(1) = sqrt (p%m**2 + dot_product (v, v))
  end subroutine momentum4_set_p_on_shell

  elemental function momentum4_is_on_shell (p) result (flag)
    class(momentum4_t), intent(in) :: p
    logical :: flag
    flag = (mass (p) - p%m) < 1e-5
  end function momentum4_is_on_shell

  subroutine momentum4_boost (p, p_boost, invert)
    class(momentum4_t), intent(inout) :: p
    type(momentum4_t), intent(in) :: p_boost
    logical, intent(in), optional :: invert
    real(default) :: g, b_abs, Estar, ppstar, pp_abs
    real(default), dimension(3) :: b, pt, pp
    if (all (p_boost%p(2:4) == 0)) return
    b = p_boost%p(2:4) / p_boost%p(1)
    b_abs = sqrt (dot_product (b, b))
    if (present(invert)) then
       if (invert) b = - b
    end if
    g = 1. / sqrt (1. - dot_product (b, b))
    ! Split in transversal and parallel components to boost vector
    pp = dot_product (p%p(2:4), b) / dot_product (b, b) * b
    pt = p%p(2:4) - pp
    pp_abs = sqrt (dot_product (pp, pp))
    ! print *, "BOOST: "
    ! print *, "β:            ", b
    ! print *, "p_parallel:   ", pp
    ! print *, "p_transverse: ", pt
    ! print *, "Δ(p, pp+pt):  ", p%p(2:4) - (pp + pt)
    ! Boost parallel to boost vector
    Estar = g * p%p(1) - g * b_abs * pp_abs
    ppstar = -g * b_abs * p%p(1) + g * pp_abs
    p%p(1) = Estar
    p%p(2:4) = ppstar * pp / pp_abs + pt
    ! print *, "BOOST ======================================================================"
    ! print *, p%p(1)
    ! p%p(1) = sqrt (dot_product (p%p(2:4), p%p(2:4)) + p%m**2)
    ! print *, p%p(1)
  end subroutine momentum4_boost
end module momentum
