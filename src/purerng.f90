! Module for pure random number generators.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module purerng

use prec, only: I10
implicit none
private

integer, public, parameter :: RNG_LECUYER = 1
integer, public, parameter :: RNG_DETERM  = 2

integer, parameter :: SIZE_LECUYER = 2

integer, parameter :: DETERM_DENOM = 1000

! from lecuyer_efficient_1988, table III
integer, parameter :: L = 2
integer(I10), parameter :: LECUYER_M(L) = [2147483563_I10, 2147483399_I10]
integer(I10), parameter :: LECUYER_A(L) = [40014_I10, 40692_I10]
integer(I10), parameter :: LECUYER_Q(L) = [53668_I10, 52774_I10]
integer(I10), parameter :: LECUYER_R(L) = [12211_I10, 3791_I10]

type, public :: rng_type
    integer(I10), allocatable, private :: seed(:)
    integer, private :: rng_num = RNG_LECUYER
contains
    procedure :: random_number => purerng_random_number
    procedure :: random_seed   => purerng_random_seed
    procedure :: set_rng_num   => set_rng_num
    procedure :: get_rng_num   => get_rng_num
    procedure :: set_determ    => set_determ
    procedure :: int           => rand_int
    procedure :: uniform       => rand_uniform
    procedure :: cauchy        => rand_cauchy
end type rng_type

contains

elemental subroutine purerng_random_number(rng, harvest)
    ! <https://fortranwiki.org/fortran/show/random_number>
    
    use prec, only: WP
    
    class(rng_type), intent(in out) :: rng
    real(WP), intent(out)           :: harvest
    
    select case (rng%rng_num)
        case (RNG_LECUYER)
            call lecuyer(rng, harvest)
        case (RNG_DETERM)
            call determ(rng, harvest)
        case default
            error stop "purerng_random_number: type of random number generator not selected."
    end select
end subroutine purerng_random_number

subroutine purerng_random_seed(rng, seed_size, put, get)
    ! <https://fortranwiki.org/fortran/show/random_seed>
    
    use prec, only: WP
    use checks, only: assert
    
    class(rng_type), intent(in out) :: rng
    
    integer, intent(out), optional                   :: seed_size
    integer(I10), intent(in), optional               :: put(:)
    integer(I10), allocatable, intent(out), optional :: get(:)
    
    integer :: values(8), i
    
    if ((.not. present(seed_size)) .and. (.not. present(put)) .and. (.not. present(get))) then
        if (allocated(rng%seed)) then
            deallocate(rng%seed)
        end if
        
        call date_and_time(values=values)
        
        select case (rng%rng_num)
            case (RNG_LECUYER)
                allocate(rng%seed(SIZE_LECUYER))
                rng%seed(1) = 1_I10 + floor(real(LECUYER_M(1) - 1_I10, WP) * real(values(7), WP) / 60.0_WP, I10)
                rng%seed(2) = 1_I10 + floor(real(LECUYER_M(2) - 1_I10, WP) * real(values(8), WP) / 1000.0_WP, I10)
            case (RNG_DETERM)
                error stop "purerng_random_seed(no args): does not apply to RNG_DETERM"
            case default
                error stop "purerng_random_seed(no args): type of random number generator not selected."
        end select
    end if
    
    if (present(seed_size)) then
        select case (rng%rng_num)
            case (RNG_LECUYER)
                seed_size = SIZE_LECUYER
            case (RNG_DETERM)
                error stop "purerng_random_seed(size): does not apply to RNG_DETERM"
            case default
                error stop "purerng_random_seed(size): type of random number generator not selected."
        end select
    end if
    
    if (present(put)) then
        if (allocated(rng%seed)) then
            deallocate(rng%seed)
        end if
        
        select case (rng%rng_num)
            case (RNG_LECUYER)
                allocate(rng%seed(SIZE_LECUYER))
            case (RNG_DETERM)
                call assert(put(1) >= 2_I10, "purerng (purerng_random_seed): For RNG_DETERM, the first put index is the " // &
                                                "index of put for the next deterministic PRNG draw. Consequently, the first " // &
                                                "index must be 2 or greater as indices 2 to the end are for the return " // &
                                                "values of the deterministic PRNG draws.")
                call assert(put(1) <= size(put, kind=I10), "For RNG_DETERM, the first put index is the " // &
                                                "index of put for the next deterministic PRNG draw. Consequently, the first " // &
                                                "index must be less than or equal to the length of the put array.")
                do concurrent (i = 2:size(put))
                    call assert(put(i) >= 0_I10, "purerng (purerng_random_seed): For RNG_DETERM, return values are " // &
                                                 "proportional to the returned PRNG, and consequently must be >= 0.")
                    call assert(put(i) <= DETERM_DENOM, "purerng (purerng_random_seed): For RNG_DETERM, return values are " // &
                                                 "proportional to the returned PRNG, and consequently must be <= DETERM_DENOM.")
                end do
            case default
                error stop "purerng_random_seed(put): type of random number generator not selected."
        end select
        
        rng%seed = put
    end if
    
    if (present(get)) then
        if (.not. allocated(rng%seed)) then
            error stop "purerng_random_seed(get): rng%seed not set."
        end if
        
        get = rng%seed
    end if
end subroutine purerng_random_seed

elemental subroutine set_rng_num(rng, rng_num)
    class(rng_type), intent(in out) :: rng
    integer, intent(in)             :: rng_num
    
    rng%rng_num = rng_num
end subroutine set_rng_num

elemental subroutine get_rng_num(rng, rng_num)
    class(rng_type), intent(in out) :: rng
    integer, intent(out)            :: rng_num
    
    rng_num = rng%rng_num
end subroutine get_rng_num

subroutine set_determ(rng, harvest)
    use prec, only: WP
    
    class(rng_type), intent(in out) :: rng
    real(WP), intent(in)            :: harvest(:)
    
    integer(I10) :: seed(size(harvest) + 1)
    integer      :: i_harvest
    
    call rng%set_rng_num(RNG_DETERM)
    seed(1) = 2_I10
    do concurrent (i_harvest = 1:size(harvest))
        seed(1 + i_harvest) = nint(real(DETERM_DENOM, WP) * harvest(i_harvest), I10)
    end do
    call rng%random_seed(put=seed)
end subroutine set_determ

elemental subroutine lecuyer(rng, harvest)
    ! Random number generator from lecuyer_efficient_1988, fig. 3.
    
    use prec, only: WP
    use checks, only: assert
    
    type(rng_type), intent(in out) :: rng
    real(WP), intent(out)          :: harvest
    
    integer(I10) :: k(L), z
    
    call assert(allocated(rng%seed), "purerng (lecuyer): seed array must be allocated")
    call assert(size(rng%seed) == 2, "purerng (lecuyer): seed array is the wrong size")
    
    ! lecuyer_efficient_1988 p. 747R
    call assert(rng%seed(1) >= 1_I10, "purerng (lecuyer): seed(1) >= 1, see lecuyer_efficient_1988 p. 747R")
    call assert(rng%seed(1) <= LECUYER_M(1) - 1_I10, "purerng (lecuyer): seed(1) <= LECUYER_M(1) - 1, " // &
                                                        "see lecuyer_efficient_1988 p. 747R")
    call assert(rng%seed(2) >= 1_I10, "purerng (lecuyer): seed(2) >= 1, lecuyer_efficient_1988 p. 747R")
    call assert(rng%seed(2) <= LECUYER_M(2) - 1_I10, "purerng (lecuyer): seed(2) <= LECUYER_M(2) - 1, " // &
                                                        "see lecuyer_efficient_1988 p. 747R")
    
    k        = rng%seed / LECUYER_Q
    rng%seed = LECUYER_A * (rng%seed - k * LECUYER_Q) - k * LECUYER_R
    where (rng%seed < 0)
        rng%seed = rng%seed + LECUYER_M
    end where
    
    z = rng%seed(1) - rng%seed(2)
    ! TODO: Consider rewriting as a `min` to improve branch prediction? Though that would require computing both.
    if (z < 1) then
        z = z + LECUYER_M(1) - 1
    end if
    
    harvest = real(z, WP) / real(LECUYER_M(1), WP)
    
    ! lecuyer_efficient_1988 p. 747R:
    ! > Notice that the function will never return 0.0 or 1.0, as long as `REAL` variables have at least 23-bit mantissa (this is
    ! > the case for most 32-bit machines).
    call assert(harvest > 0.0_WP, "purerng (lecuyer): harvest <= 0")
    call assert(harvest < 1.0_WP, "purerng (lecuyer): harvest >= 1")
    call assert(z > 0_I10, "purerng (lecuyer): z <= 0")
    call assert(z < LECUYER_M(1), "purerng (lecuyer): z >= LECUYER_M(1)")
    
    ! Same `seed` bounds as before.
    call assert(rng%seed(1) >= 1_I10, "purerng (lecuyer): seed(1) >= 1, see lecuyer_efficient_1988 p. 747R")
    call assert(rng%seed(1) <= LECUYER_M(1) - 1_I10, "purerng (lecuyer): seed(1) <= LECUYER_M(1) - 1, " // &
                                                        "see lecuyer_efficient_1988 p. 747R")
    call assert(rng%seed(2) >= 1_I10, "purerng (lecuyer): seed(2) >= 1, lecuyer_efficient_1988 p. 747R")
    call assert(rng%seed(2) <= LECUYER_M(2) - 1_I10, "purerng (lecuyer): seed(2) <= LECUYER_M(2) - 1, " // &
                                                        "see lecuyer_efficient_1988 p. 747R")
end subroutine lecuyer

elemental subroutine determ(rng, harvest)
    use prec, only: WP
    use checks, only: assert
    
    ! deterministic "random" number
    
    type(rng_type), intent(in out) :: rng
    real(WP), intent(out)          :: harvest
    
    harvest = real(rng%seed(rng%seed(1)), WP) / real(DETERM_DENOM, WP)
    
    ! First index of seed is which value is read next.
    rng%seed(1) = rng%seed(1) + 1_I10
    if (rng%seed(1) > size(rng%seed, kind=I10)) then
        rng%seed(1) = 2_I10
    end if
    
    call assert(harvest >= 0.0_WP, "purerng (determ): harvest <= 0")
    call assert(harvest <= 1.0_WP, "purerng (determ): harvest >= 1")
    call assert(rng%seed(1) >= 2_I10, "purerng (determ): z <= 0")
    call assert(rng%seed(1) <= size(rng%seed, kind=I10), "purerng (determ): z >= LECUYER_M(1)")
end subroutine determ

elemental subroutine rand_int(rng, lower_bound, upper_bound, nu)
    use prec, only: WP
    use checks, only: assert
    
    class(rng_type), intent(in out) :: rng
    integer, intent(in)             :: lower_bound, upper_bound
    integer, intent(out)            :: nu
    
    real(WP) :: harvest
    
    call rng%random_number(harvest)
    
    ! The `min` function makes this not return `upper_bound + 1` when `r = 1.0_WP`.
    nu = min(lower_bound + floor(real(upper_bound + 1 - lower_bound, WP) * harvest), upper_bound)
    
    call assert(nu >= lower_bound, "purerng (rand_int): nu >= lower_bound violated")
    call assert(nu <= upper_bound, "purerng (rand_int): nu >= upper_bound violated")
end subroutine rand_int

elemental subroutine rand_uniform(rng, lower_bound, upper_bound, nu)
    ! Returns a uniform random variable.
    
    use prec, only: WP
    use checks, only: assert
    
    class(rng_type), intent(in out) :: rng
    real(WP), intent(in)            :: lower_bound, upper_bound
    real(WP), intent(out)           :: nu
    
    real(WP) :: harvest
    
    call rng%random_number(harvest)
    nu = lower_bound + (upper_bound - lower_bound) * harvest
    
    call assert(nu >= lower_bound, "purerng (rand_uniform): nu >= lower_bound violated")
    call assert(nu <= upper_bound, "purerng (rand_uniform): nu >= upper_bound violated")
end subroutine rand_uniform

elemental subroutine rand_cauchy(rng, m, b, nu)
    ! Returns a Cauchy random variable with inverse transform sampling.
    ! Notation follows <https://mathworld.wolfram.com/CauchyDistribution.html>.
    
    use prec, only: WP, PI
    
    class(rng_type), intent(in out) :: rng
    real(WP), intent(in)            :: m ! median
    real(WP), intent(in)            :: b ! half width
    real(WP), intent(out)           :: nu
    
    real(WP) :: harvest
    
    call rng%random_number(harvest)
    nu = m + b * tan(PI * (harvest - 0.5_WP))
end subroutine rand_cauchy

end module purerng
