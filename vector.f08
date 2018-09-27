module momentum
  use ISO_FORTRAN_ENV, only: OUTPUT_UNIT
  use constants, only: default

  implicit none

  private

  type :: momentum4_t
     ! TODO private fields for momentum4_t
     ! private
     real(default) :: p_sq = 0
     real(default), dimension(4) :: p = &
          [0._default, 0._default, 0._default, 0._default]
   contains
     procedure :: write => momentum4_write
     procedure :: get_p => momentum4_get_p
     procedure :: get_mass => momentum4_get_mass
     procedure :: get_energy => momentum4_get_energy
     procedure :: get_momentum => momentum4_get_momentum
     procedure :: rescale_energy => momentum4_rescale_energy
     procedure :: set_p => momentum4_set_p
     procedure :: set_mass => momentum4_set_mass
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
    rho = sqrt ((M1**2 - (M2 + m)**2) * (M1**2 - (M2 - m)**2))
    rho = rho / (8._default * M1**2)
  end function rho

  elemental real(default) function squared_momentum4 (p) result (p_sq)
    type(momentum4_t), intent(in) :: p
    p_sq = p%p(1)**2 - (dot_product (p%p(2:4), p%p(2:4)))
  end function squared_momentum4

  real(default) function squared_array4 (p) result (p_sq)
    real(default), dimension(4), intent(in) :: p
    p_sq = p(1)**2 - (dot_product (p(2:4), p(2:4)))
  end function squared_array4

  real(default) function mass_momentum4 (p) result (m)
    type(momentum4_t), intent(in) :: p
    real(default) :: p_sq
    p_sq = p%p(1)**2 - (dot_product (p%p(2:4), p%p(2:4)))
    if (p_sq >= 0) then
       m = sqrt (p_sq)
    else
       print *, "TAG"
       m = sqrt (-p_sq)
    end if
  end function mass_momentum4

  real(default) function mass_array4 (p) result (m)
    real(default), dimension(4), intent(in) :: p
    real(default) :: p_sq
    p_sq = p(1)**2 - (dot_product (p(2:4), p(2:4)))
    if (p_sq >= 0) then
       m = sqrt(p_sq)
    else
       print *, "TAG"
       m = sqrt (-p_sq)
    end if
  end function mass_array4

  type(momentum4_t) function momentum4_init (v) result (p)
    real(default), dimension(4), intent(in) :: v
    p%p = v
    p%p_sq = squared (v)
  end function momentum4_init

  type(momentum4_t) function momentum4_init_at_rest (m) result (p)
    real(default), intent(in) :: m
    p%p = 0
    p%p_sq = (m**2)
    p%p(1) = m
  end function momentum4_init_at_rest

  subroutine momentum4_write (p, unit)
    class(momentum4_t), intent(in) :: p
    integer, intent(in), optional :: unit
    integer :: u
    u = OUTPUT_UNIT; if (present (unit)) u = unit
    write (u, "(12X,A,12X,A,51X,A)") "M", "~", "P"
    write (u, "(6(1X,F12.4))") p%p_sq, mass(p)**2, p%p
  end subroutine momentum4_write

  function momentum4_get_p (p) result (v)
    class(momentum4_t), intent(in) :: p
    real(default), dimension(4) :: v
    v = p%p
  end function momentum4_get_p

  real(default) function momentum4_get_mass (p) result (m)
    class(momentum4_t), intent(in) :: p
    if (p%p_sq >= 0) then
       m = sqrt(p%p_sq)
    else
       print *, "TAG"
       m = sqrt(-p%p_sq)
    end if
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

  elemental type(momentum4_t) function momentum4_add (p1, p2) result (p)
    class(momentum4_t), intent(in) :: p1
    type(momentum4_t), intent(in) :: p2
    p%p = p1%p + p2%p
    p%p_sq = squared(p)
  end function momentum4_add

  elemental type(momentum4_t) function momentum4_subtract (p1, p2) result (p)
    class(momentum4_t), intent(in) :: p1
    type(momentum4_t), intent(in) :: p2
    p%p = p1%p - p2%p
    p%p_sq = squared(p)
  end function momentum4_subtract

  subroutine momentum4_set_p (p, v)
    class(momentum4_t), intent(inout) :: p
    real(default), dimension(3), intent(in) :: v
    p%p(2:4) = v
    ! TODO Rescale energy to ensure on-shell condition
    call p%rescale_energy ()
  end subroutine momentum4_set_p

  subroutine momentum4_set_mass (p, m)
    class(momentum4_t), intent(inout) :: p
    real(default), intent(in) :: m
    p%p_sq = m**2
  end subroutine momentum4_set_mass

  subroutine momentum4_rescale_energy (p)
    class(momentum4_t), intent(inout) :: p
    real(default) :: e_sq
    e_sq = dot_product(p%p(2:4), p%p(2:4)) + p%p_sq
    if (e_sq >= 0) then
       p%p(1) = sqrt(e_sq)
    else
       e_sq = 0
    end if
  end subroutine momentum4_rescale_energy

  subroutine momentum4_boost (p, p_boost)
    class(momentum4_t), intent(inout) :: p
    type(momentum4_t), intent(in) :: p_boost
    real(default) :: g, b_abs, Estar, ppstar, pp_abs
    real(default), dimension(3) :: b, pt, pp
    if (all (p_boost%p(2:4) == 0)) return
    b = p_boost%p(2:4) / p_boost%p(1)
    b_abs = sqrt (dot_product (b, b))
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
    ! call p%rescale_energy ()
  end subroutine momentum4_boost
end module momentum

! pt.b = 0
! pt = p.b b / |b|
! pp = p - pt!
! E = g*E - g * |b| pp = ... - g * |b| * (p - p.b b)/|b|
