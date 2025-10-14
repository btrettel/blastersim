! Module for forward-mode automatic differentiation.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

module fmad

use prec, only: WP
use checks, only: assert
implicit none
private

public :: sqrt, tanh, log, exp, merge, max, min, abs, sin, cos, tan, atanh
public :: var, stdev
public :: f

! Both the dependent and independent variables need to be of type `ad`.
type, public :: ad
    real(WP)              :: v    ! value
    real(WP), allocatable :: d(:) ! derivative values
contains
    procedure :: init
    procedure :: init_const
    procedure, private :: ad_ad_add, ad_real_add, ad_add_unary
    procedure, private, pass(ad_right) :: real_ad_add
    generic, public :: operator(+) => ad_ad_add, ad_real_add, real_ad_add, ad_add_unary
    procedure, private :: ad_ad_subtract, ad_real_subtract, ad_subtract_unary
    procedure, private, pass(ad_right) :: real_ad_subtract
    generic, public :: operator(-) => ad_ad_subtract, ad_real_subtract, real_ad_subtract, ad_subtract_unary
    procedure, private :: ad_ad_multiply, ad_real_multiply
    procedure, private, pass(ad_right) :: real_ad_multiply
    generic, public :: operator(*) => ad_ad_multiply, ad_real_multiply, real_ad_multiply
    procedure, private :: ad_ad_divide, ad_real_divide
    procedure, private, pass(ad_right) :: real_ad_divide
    generic, public :: operator(/) => ad_ad_divide, ad_real_divide, real_ad_divide
    procedure, private :: ad_real_exponentiate, ad_integer_exponentiate
    generic, public :: operator(**) => ad_real_exponentiate, ad_integer_exponentiate
    procedure, private :: lt_ad
    generic, public :: operator(<) => lt_ad
    procedure, private :: le_ad
    generic, public :: operator(<=) => le_ad
    procedure, private :: gt_ad
    generic, public :: operator(>) => gt_ad
    procedure, private :: ge_ad
    generic, public :: operator(>=) => ge_ad
end type ad

interface sqrt
    module procedure :: ad_sqrt
end interface sqrt

interface tanh
    module procedure :: ad_tanh
end interface tanh

interface log
    module procedure :: ad_log
end interface log

interface exp
    module procedure :: ad_exp
end interface exp

interface merge
    module procedure :: ad_ad_merge
end interface merge

interface max
    module procedure :: ad_ad_max_2
    module procedure :: real_ad_max_2
end interface max

interface min
    module procedure :: ad_ad_min_2
    module procedure :: real_ad_min_2
end interface min

interface abs
    module procedure :: ad_abs
end interface abs

interface sin
    module procedure :: ad_sin
end interface sin

interface cos
    module procedure :: ad_cos
end interface cos

interface tan
    module procedure :: ad_tan
end interface tan

interface atanh
    module procedure :: ad_atanh
end interface atanh

contains

! Constructors for `ad`
! ---------------------

elemental subroutine init(x, v, n, n_d)
    use checks, only: is_close
    
    class(ad), intent(in out) :: x    ! `class` can't be `intent(out)` and `pure`?!?
    real(WP), intent(in)      :: v    ! value of variable to set
    integer, intent(in)       :: n, & ! variable number represented (sets the appropriate derivative)
                                 n_d ! total number of differentiable variables

    integer :: i_dv ! loop index
    
    call assert(n_d >= 0, "fmad (init): n_d must be zero or more")
    
    x%v = v
    
    allocate(x%d(n_d))
    
    if (n_d > 0) then
        call assert(n >= 1, "fmad (init): n must be 1 or more")
        call assert(n <= n_d, "fmad (init): n must be n_d or less")
        do concurrent (i_dv = 1:n_d)
            x%d(i_dv) = merge(1.0_WP, 0.0_WP, i_dv == n)
        end do
        call assert(any(is_close(x%d, 1.0_WP)), "fmad (init): at least one derivative set")
        call assert(any(is_close(x%d, 0.0_WP)) .or. (n_d <= 1), "fmad (init): at least one derivative not set")
    end if
end subroutine init

elemental subroutine init_const(x, v, n_d)
    class(ad), intent(in out) :: x    ! `class` can't be `intent(out)` and `pure`?!?
    real(WP), intent(in)      :: v    ! value of constant to set
    integer, intent(in)       :: n_d ! total number of differentiable variables
    
    call assert(n_d >= 0, "fmad (init): n_d must be zero or more")
    
    allocate(x%d(n_d))

    x%v = v
    if (n_d > 0) then
        x%d = 0.0_WP
    end if
end subroutine init_const

! Operator procedures
! -------------------

elemental function ad_ad_add(ad_left, ad_right)
    ! Adds two `ad`s.
    
    class(ad), intent(in) :: ad_left, ad_right
    
    type(ad) :: ad_ad_add
    
    call assert(allocated(ad_left%d), "fmad (ad_ad_add): ad_left%d must be allocated")
    call assert(allocated(ad_right%d), "fmad (ad_ad_add): ad_right%d must be allocated")
    
    ad_ad_add%v = ad_left%v + ad_right%v
    ad_ad_add%d = ad_left%d + ad_right%d
end function ad_ad_add

elemental function ad_real_add(ad_in, real_in)
    ! Adds an `ad` and a `real`.
    
    class(ad), intent(in) :: ad_in
    real(WP), intent(in)  :: real_in
    
    type(ad) :: ad_real_add
    
    call assert(allocated(ad_in%d), "fmad (ad_real_add): ad_in%d must be allocated")

    ad_real_add%v = ad_in%v + real_in
    ad_real_add%d = ad_in%d
end function ad_real_add

elemental function real_ad_add(real_left, ad_right)
    ! Adds a `real` and an `ad`.
    
    real(WP), intent(in)  :: real_left
    class(ad), intent(in) :: ad_right
    
    type(ad) :: real_ad_add
    
    call assert(allocated(ad_right%d), "fmad (real_ad_add): ad_right%d must be allocated")

    real_ad_add%v = real_left + ad_right%v
    real_ad_add%d = ad_right%d
end function real_ad_add

elemental function ad_ad_subtract(ad_left, ad_right)
    ! Subtracts two `ad`s.
    
    class(ad), intent(in) :: ad_left, ad_right
    
    type(ad) :: ad_ad_subtract
    
    call assert(allocated(ad_left%d), "fmad (ad_ad_subtract): ad_left%d must be allocated")
    call assert(allocated(ad_right%d), "fmad (ad_ad_subtract): ad_right%d must be allocated")

    ad_ad_subtract%v = ad_left%v - ad_right%v
    ad_ad_subtract%d = ad_left%d - ad_right%d
end function ad_ad_subtract

elemental function ad_real_subtract(ad_left, real_right)
    ! Subtracts a `real` from an `ad`.
    
    class(ad), intent(in) :: ad_left
    real(WP), intent(in)  :: real_right
    
    type(ad) :: ad_real_subtract
    
    call assert(allocated(ad_left%d), "fmad (ad_real_subtract): ad_left%d must be allocated")

    ad_real_subtract%v = ad_left%v - real_right
    ad_real_subtract%d = ad_left%d
end function ad_real_subtract

elemental function real_ad_subtract(real_left, ad_right)
    ! Subtracts a `real` from an `ad`.

    real(WP), intent(in)  :: real_left
    class(ad), intent(in) :: ad_right
    
    type(ad) :: real_ad_subtract
    
    call assert(allocated(ad_right%d), "fmad (real_ad_subtract): ad_right%d must be allocated")
    
    real_ad_subtract%v = real_left - ad_right%v
    real_ad_subtract%d = -ad_right%d
end function real_ad_subtract

elemental function ad_subtract_unary(ad_in)
    ! Returns `-rd`.
    
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_subtract_unary
    
    call assert(allocated(ad_in%d), "fmad (ad_subtract_unary): ad_in%d must be allocated")

    ad_subtract_unary%v = -ad_in%v
    ad_subtract_unary%d = -ad_in%d
end function ad_subtract_unary

elemental function ad_add_unary(ad_in)
    ! Returns `+rd`.
    
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_add_unary
    
    call assert(allocated(ad_in%d), "fmad (ad_add_unary): ad_in%d must be allocated")

    ad_add_unary%v = ad_in%v
    ad_add_unary%d = ad_in%d
end function ad_add_unary

elemental function ad_ad_multiply(ad_left, ad_right)
    ! Multiplies two `ad`s.
    
    class(ad), intent(in) :: ad_left, ad_right
    
    type(ad) :: ad_ad_multiply
    
    call assert(allocated(ad_left%d), "fmad (ad_ad_multiply): ad_left%d must be allocated")
    call assert(allocated(ad_right%d), "fmad (ad_ad_multiply): ad_right%d must be allocated")

    ad_ad_multiply%v = ad_left%v * ad_right%v
    ad_ad_multiply%d = ad_left%d * ad_right%v + ad_left%v * ad_right%d
end function ad_ad_multiply

elemental function ad_real_multiply(ad_left, real_right)
    ! Multiplies an `ad` by a `real`.
    
    class(ad), intent(in) :: ad_left
    real(WP), intent(in)  :: real_right
    
    type(ad) :: ad_real_multiply
    
    call assert(allocated(ad_left%d), "fmad (ad_real_multiply): ad_left%d must be allocated")
    
    ad_real_multiply%v = ad_left%v * real_right
    ad_real_multiply%d = ad_left%d * real_right
end function ad_real_multiply

elemental function real_ad_multiply(real_left, ad_right)
    ! Multiplies a `real` by an `ad`.
    
    real(WP), intent(in)  :: real_left
    class(ad), intent(in) :: ad_right
    
    type(ad) :: real_ad_multiply
    
    call assert(allocated(ad_right%d), "fmad (real_ad_multiply): ad_right%d must be allocated")
    
    real_ad_multiply%v = real_left * ad_right%v
    real_ad_multiply%d = real_left * ad_right%d
end function real_ad_multiply

elemental function ad_ad_divide(ad_left, ad_right)
    ! Divides two `ad`.
    
    class(ad), intent(in) :: ad_left, ad_right
    
    type(ad) :: ad_ad_divide
    
    call assert(allocated(ad_left%d), "fmad (ad_ad_divide): ad_left%d must be allocated")
    call assert(allocated(ad_right%d), "fmad (ad_ad_divide): ad_right%d must be allocated")
    
    ad_ad_divide%v = ad_left%v / ad_right%v
    ad_ad_divide%d = (ad_left%d * ad_right%v - ad_left%v * ad_right%d) / (ad_right%v**2)
end function ad_ad_divide

elemental function ad_real_divide(ad_left, real_right)
    ! Divides an `ad` by a `real`.
    
    class(ad), intent(in) :: ad_left
    real(WP), intent(in)  :: real_right
    
    type(ad) :: ad_real_divide
    
    call assert(allocated(ad_left%d), "fmad (ad_real_divide): ad_left%d must be allocated")
    
    ad_real_divide%v = ad_left%v / real_right
    ad_real_divide%d = ad_left%d / real_right
end function ad_real_divide

elemental function real_ad_divide(real_left, ad_right)
    ! Divides a `real` by an `ad`.
    
    real(WP), intent(in)  :: real_left
    class(ad), intent(in) :: ad_right
    
    type(ad) :: real_ad_divide
    
    call assert(allocated(ad_right%d), "fmad (real_ad_divide): ad_right%d must be allocated")
    
    real_ad_divide%v = real_left / ad_right%v
    real_ad_divide%d = -real_left * ad_right%d / (ad_right%v**2)
end function real_ad_divide

elemental function ad_real_exponentiate(ad_in, real_in)
    use checks, only: is_close
    
    ! Exponentiates an `ad` by a `real`.
    
    class(ad), intent(in) :: ad_in
    real(WP), intent(in)  :: real_in
    
    type(ad) :: ad_real_exponentiate
    
    call assert(.not. (is_close(ad_in%v, 0.0_WP) .and. (real_in <= 0.0_WP)), &
                    "fmad (ad_real_exponentiate): exponent is negative or zero and argument is zero")
    call assert(allocated(ad_in%d), "fmad (ad_real_exponentiate): ad_in%d must be allocated")
    
    ad_real_exponentiate%v = ad_in%v**real_in
    ad_real_exponentiate%d = real_in*(ad_in%v**(real_in - 1.0_WP))*ad_in%d
end function ad_real_exponentiate

elemental function ad_integer_exponentiate(ad_in, integer_in)
    use checks, only: is_close
    
    ! Exponentiates an `ad` by an `integer`.
    
    class(ad), intent(in) :: ad_in
    integer, intent(in)   :: integer_in
    
    type(ad) :: ad_integer_exponentiate
    
    call assert(.not. (is_close(ad_in%v, 0.0_WP) .and. (integer_in <= 0)), &
                    "fmad (ad_integer_exponentiate): exponent is negative or zero and argument is zero")
    call assert(allocated(ad_in%d), "fmad (ad_integer_exponentiate): ad_in%d must be allocated")

    ad_integer_exponentiate%v = ad_in%v**integer_in
    ad_integer_exponentiate%d = real(integer_in, WP)*(ad_in%v**(integer_in - 1))*ad_in%d
end function ad_integer_exponentiate

! No `rd**rd` as that's not likely to happen in CFD.

elemental function lt_ad(ad_left, ad_right)
    class(ad), intent(in) :: ad_left
    type(ad), intent(in)  :: ad_right
    
    logical :: lt_ad
    
    lt_ad = ad_left%v < ad_right%v
    
    call assert(allocated(ad_left%d), "fmad (lt_ad): ad_left%d must be allocated")
    call assert(allocated(ad_right%d), "fmad (lt_ad): ad_right%d must be allocated")
end function lt_ad

elemental function le_ad(ad_left, ad_right)
    class(ad), intent(in) :: ad_left
    type(ad), intent(in)  :: ad_right
    
    logical :: le_ad
    
    le_ad = ad_left%v <= ad_right%v
    
    call assert(allocated(ad_left%d), "fmad (le_ad): ad_left%d must be allocated")
    call assert(allocated(ad_right%d), "fmad (le_ad): ad_right%d must be allocated")
end function le_ad

elemental function gt_ad(ad_left, ad_right)
    class(ad), intent(in) :: ad_left
    type(ad), intent(in)  :: ad_right
    
    logical :: gt_ad
    
    gt_ad = ad_left%v > ad_right%v
    
    call assert(allocated(ad_left%d), "fmad (gt_ad): ad_left%d must be allocated")
    call assert(allocated(ad_right%d), "fmad (gt_ad): ad_right%d must be allocated")
end function gt_ad

elemental function ge_ad(ad_left, ad_right)
    class(ad), intent(in) :: ad_left
    type(ad), intent(in)  :: ad_right
    
    logical :: ge_ad
    
    ge_ad = ad_left%v >= ad_right%v
    
    call assert(allocated(ad_left%d), "fmad (ge_ad): ad_left%d must be allocated")
    call assert(allocated(ad_right%d), "fmad (ge_ad): ad_right%d must be allocated")
end function ge_ad

elemental function ad_sqrt(ad_in)
    ! Takes the square root of an `ad`.
    
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_sqrt
    
    call assert(ad_in%v > 0.0_WP, "fmad (ad_sqrt): argument is zero or negative")
    call assert(allocated(ad_in%d), "fmad (ad_sqrt): ad_in%d must be allocated")
    
    ad_sqrt%v = sqrt(ad_in%v)
    ad_sqrt%d = ad_in%d/(2.0_WP * sqrt(ad_in%v))
end function ad_sqrt

elemental function ad_tanh(ad_in)
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_tanh
    
    call assert(allocated(ad_in%d), "fmad (ad_tanh): ad_in%d must be allocated")
    
    ad_tanh%v = tanh(ad_in%v)
    ad_tanh%d = ad_in%d*(1.0_WP - tanh(ad_in%v)**2)
end function ad_tanh

elemental function ad_log(ad_in)
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_log
    
    call assert(ad_in%v > 0.0_WP, "fmad (ad_log): argument is zero or negative")
    call assert(allocated(ad_in%d), "fmad (ad_log): ad_in%d must be allocated")
    
    ad_log%v = log(ad_in%v)
    ad_log%d = ad_in%d/ad_in%v
end function ad_log

elemental function ad_exp(ad_in)
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_exp
    
    call assert(allocated(ad_in%d), "fmad (ad_exp): ad_in%d must be allocated")
    
    ad_exp%v = exp(ad_in%v)
    ad_exp%d = ad_in%d*exp(ad_in%v)
end function ad_exp

pure function ad_ad_merge(ad_left, ad_right, mask)
    class(ad), intent(in) :: ad_left, ad_right
    logical, intent(in)   :: mask
    
    type(ad) :: ad_ad_merge
    
    call assert(allocated(ad_left%d), "fmad (ad_ad_merge): ad_left%d must be allocated")
    call assert(allocated(ad_right%d), "fmad (ad_ad_merge): ad_right%d must be allocated")
    
    if (mask) then
        ad_ad_merge = ad_left
    else
        ad_ad_merge = ad_right
    end if
end function ad_ad_merge

pure function ad_ad_max_2(ad_left, ad_right)
    class(ad), intent(in) :: ad_left, ad_right
    
    type(ad) :: ad_ad_max_2
    
    ad_ad_max_2 = ad_ad_merge(ad_left, ad_right, ad_left%v >= ad_right%v)
end function ad_ad_max_2

pure function real_ad_max_2(real_left, ad_right)
    ! Related: <https://en.wikipedia.org/wiki/Rectifier_(neural_networks)>
    
    real(WP), intent(in)  :: real_left
    class(ad), intent(in) :: ad_right
    
    type(ad) :: real_ad_max_2
    
    type(ad) :: ad_left
    
    call ad_left%init_const(real_left, size(ad_right%d))
    real_ad_max_2 = ad_ad_max_2(ad_left, ad_right)
end function real_ad_max_2

pure function ad_ad_min_2(ad_left, ad_right)
    class(ad), intent(in) :: ad_left, ad_right
    
    type(ad) :: ad_ad_min_2
    
    ad_ad_min_2 = ad_ad_merge(ad_left, ad_right, ad_left%v <= ad_right%v)
end function ad_ad_min_2

pure function real_ad_min_2(real_left, ad_right)
    real(WP), intent(in)  :: real_left
    class(ad), intent(in) :: ad_right
    
    type(ad) :: real_ad_min_2
    
    type(ad) :: ad_left
    
    call ad_left%init_const(real_left, size(ad_right%d))
    real_ad_min_2 = ad_ad_min_2(ad_left, ad_right)
end function real_ad_min_2

elemental function ad_abs(ad_in)
    use checks, only: is_close
    
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_abs
    
    call assert(allocated(ad_in%d), "fmad (ad_abs): ad_in%d must be allocated")
    
    if (.not. is_close(ad_in%v, 0.0_WP)) then
        ad_abs%v = abs(ad_in%v)
        ad_abs%d = ad_in%d*(ad_in%v/abs(ad_in%v))
    else
        call ad_abs%init_const(0.0_WP, size(ad_in%d))
    end if
    
    call assert(ad_abs%v >= 0.0_WP, "fmad (ad_abs): value is negative")
end function ad_abs

elemental function ad_sin(ad_in)
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_sin
    
    call assert(allocated(ad_in%d), "fmad (ad_sin): ad_in%d must be allocated")
    
    ad_sin%v = sin(ad_in%v)
    ad_sin%d = ad_in%d*cos(ad_in%v)
end function ad_sin

elemental function ad_cos(ad_in)
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_cos
    
    call assert(allocated(ad_in%d), "fmad (ad_cos): ad_in%d must be allocated")
    
    ad_cos%v = cos(ad_in%v)
    ad_cos%d = -ad_in%d*sin(ad_in%v)
end function ad_cos

elemental function ad_tan(ad_in)
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_tan
    
    call assert(allocated(ad_in%d), "fmad (ad_tan): ad_in%d must be allocated")
    
    ad_tan%v = tan(ad_in%v)
    ad_tan%d = ad_in%d/(cos(ad_in%v)**2)
end function ad_tan

elemental function ad_atanh(ad_in)
    class(ad), intent(in) :: ad_in
    
    type(ad) :: ad_atanh
    
    call assert(allocated(ad_in%d), "fmad (ad_atanh): ad_in%d must be allocated")
    call assert(abs(ad_in%v) < 1.0_WP, "fmad (ad_atanh): abs(ad_in%v) < 1 violated")
    
    ad_atanh%v = atanh(ad_in%v)
    ad_atanh%d = ad_in%d/(1.0_WP - ad_in%v**2)
end function ad_atanh

pure function var(x, sigmas)
    ! Computes the variance with the FOSM method.
    ! Assumes each random variable is independent.
    ! See putko_approach_2001 section 2.2.
    
    use checks, only: assert_dimension
    
    type(ad), intent(in) :: x(:)
    real(WP), intent(in) :: sigmas(:)
    
    real(WP), allocatable :: var(:)
    
    integer :: i_x, i_dv
    
    call assert(all(sigmas >= 0.0_WP), "ad (var): input sigmas must be >= 0")
    
    allocate(var(size(x)))
    var = 0.0_WP
    do concurrent (i_x = 1:size(x))
        call assert_dimension(x(i_x)%d, sigmas)
        do i_dv = 1, size(x(1)%d) ! SERIAL
            var(i_x) = var(i_x) + (x(i_x)%d(i_dv) * sigmas(i_dv))**2
        end do
    end do
    
    call assert(all(var >= 0.0_WP), "ad (var): output var must be >= 0")
end function var

pure function stdev(x, sigmas)
    type(ad), intent(in) :: x(:)
    real(WP), intent(in) :: sigmas(:)
    
    real(WP), allocatable :: stdev(:)
    
    stdev = sqrt(var(x, sigmas))
    
    call assert(all(stdev >= 0.0_WP), "ad (stdev): output stdev must be >= 0")
end function stdev

pure function f(x, y)
    ! Test function. It's here because nvfortran has a bug if it's an internal procedure in the tests.
    
    type(ad), intent(in) :: x, y
    
    type(ad) :: f
    
    call assert(y%v > 0.0_WP, "fmad (f): y is zero")
    
    f = (2.0_WP * x * y - x**2) / y + y
end function f

end module fmad
