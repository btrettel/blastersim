! Module for procedures used for run-time checks.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module checks

use prec, only: WP
implicit none
private

! <https://fortran-lang.discourse.group/t/suggestion-findloc-tolerance/5131/5>
! I started `TOL_FACTOR` at 5, but increased it to 25 based on experience.
real(WP), public, parameter :: TOL_FACTOR = 25.0_WP

public :: abs_tolerance
public :: is_close, all_close
public :: assert
public :: assert_dimension
public :: check
public :: print, write, flush

interface all_close
    module procedure all_close_rank_1
    module procedure all_close_rank_1_rank_0
end interface all_close

interface assert_dimension
    ! LATER: nvfortran doesn't seem to properly support "assumed-rank" arrays, so I made 3 different subroutines.
    ! nvfortran will compile, but the run-time behavior will be wrong.
    
    module procedure assert_dimension_rank_1
    module procedure assert_dimension_rank_2
    module procedure assert_dimension_rank_3
end interface assert_dimension

interface check
    module procedure impure_check
    module procedure pure_check
end interface check

type, public :: logger_type
    integer                        :: unit
    character(:), allocatable      :: message
    type(logger_type), allocatable :: history
end type logger_type

contains

pure function abs_tolerance(input_real_1, input_real_2)
    ! <https://community.intel.com/t5/Intel-Fortran-Compiler/Real-number-tolerance-and-if/m-p/824459#M49399>
    ! <https://fortran-lang.discourse.group/t/suggestion-findloc-tolerance/5131/5>
    ! <https://community.intel.com/t5/Intel-Fortran-Compiler/Floating-point-equivalence-check/m-p/1087149>
    
    real(WP), intent(in) :: input_real_1, input_real_2
    
    real(WP) :: abs_tolerance
    
    abs_tolerance = TOL_FACTOR * spacing(max(abs(input_real_1), abs(input_real_2)))
end function abs_tolerance

elemental function is_close(input_real_1, input_real_2, rel_tol, abs_tol)
    ! Determine whether two reals are close.
    
    ! Interface based on stdlib, though implementation is different.
    ! <https://stdlib.fortran-lang.org/page/specs/stdlib_math.html#is_close-function>
    
    ! Also see:
    ! <https://numpy.org/doc/stable/reference/generated/numpy.isclose.html>
    ! <https://docs.python.org/3/library/math.html#math.isclose>
    
    real(WP), intent(in)           :: input_real_1, input_real_2
    real(WP), intent(in), optional :: rel_tol, abs_tol
    
    logical :: is_close
    
    real(WP) :: rel_tol_, abs_tol_, tol
    
    if (present(rel_tol)) then
        rel_tol_ = rel_tol
    else
        rel_tol_ = 0.0_WP
    end if
    
    if (present(abs_tol)) then
        abs_tol_ = abs_tol
    else
        abs_tol_ = abs_tolerance(input_real_1, input_real_2)
    end if
    
    tol = max(rel_tol_ * abs(input_real_1), rel_tol_ * abs(input_real_2), abs_tol_)
    
    if (abs(input_real_1 - input_real_2) < tol) then
        is_close = .true.
    else
        is_close = .false.
    end if
end function is_close

pure function all_close_rank_1(input_real_1, input_real_2, rel_tol, abs_tol)
    ! like <https://numpy.org/doc/stable/reference/generated/numpy.allclose.html>
    
    real(WP), intent(in)           :: input_real_1(:), input_real_2(:)
    real(WP), intent(in), optional :: rel_tol, abs_tol
    
    logical :: all_close_rank_1
    
    integer :: i, n_match
    
    real(WP) :: rel_tol_, abs_tol_
    
    if (present(rel_tol)) then
        rel_tol_ = rel_tol
    else
        rel_tol_ = 0.0_WP
    end if
    
    if (present(abs_tol)) then
        abs_tol_ = abs_tol
    else
        abs_tol_ = abs_tolerance(maxval(abs(input_real_1)), maxval(abs(input_real_2)))
    end if
    
    call assert_dimension(input_real_1, input_real_2)
    
    all_close_rank_1 = .false.
    n_match = 0
    ! TODO: Convert this to use a reduction when that's possible in all compilers.
    do i = lbound(input_real_1, dim=1), ubound(input_real_1, dim=1) ! SERIAL
        if (is_close(input_real_1(i), input_real_2(i), rel_tol=rel_tol_, abs_tol=abs_tol_)) then
            n_match = n_match + 1
        end if
    end do
    
    if (n_match == size(input_real_1)) then
        all_close_rank_1 = .true.
    else
        all_close_rank_1 = .false.
    end if
end function all_close_rank_1

pure function all_close_rank_1_rank_0(input_real_1, input_real_2, rel_tol, abs_tol)
    real(WP), intent(in)           :: input_real_1(:), input_real_2
    real(WP), intent(in), optional :: rel_tol, abs_tol
    
    logical :: all_close_rank_1_rank_0
    
    integer :: i, n_match
    
    real(WP) :: rel_tol_, abs_tol_
    
    if (present(rel_tol)) then
        rel_tol_ = rel_tol
    else
        rel_tol_ = 0.0_WP
    end if
    
    if (present(abs_tol)) then
        abs_tol_ = abs_tol
    else
        abs_tol_ = abs_tolerance(maxval(abs(input_real_1)), abs(input_real_2))
    end if
    
    all_close_rank_1_rank_0 = .false.
    n_match = 0
    ! TODO: Convert this to use a reduction when that's possible in all compilers.
    do i = lbound(input_real_1, dim=1), ubound(input_real_1, dim=1) ! SERIAL
        if (is_close(input_real_1(i), input_real_2, rel_tol=rel_tol_, abs_tol=abs_tol_)) then
            n_match = n_match + 1
        end if
    end do
    
    if (n_match == size(input_real_1)) then
        all_close_rank_1_rank_0 = .true.
    else
        all_close_rank_1_rank_0 = .false.
    end if
end function all_close_rank_1_rank_0

elemental subroutine assert(condition, message)
    ! To get line number and specific values of variables, use a debugger.
    
    use build, only: DEBUG
    
    logical, intent(in)          :: condition
    character(len=*), intent(in) :: message
    
    character(len=:), allocatable :: full_message
    
    if (DEBUG) then
        if (.not. condition) then
            ! Why not concatenate the strings on the `error stop` line?
            ! That leads to ifx garbling the error message as of version `ifx (IFX) 2024.0.2 20231213`.
            full_message = "***" // new_line("a") // "ASSERTION FAILED. " // message
            
            ! Why is the message in all caps? To make it more noticeable.
            
            error stop full_message
        end if
    end if
end subroutine assert

pure subroutine assert_dimension_rank_1(a, b)
    class(*), intent(in) :: a(:), b(:)
    
    call assert(all(lbound(a) == lbound(b)), "checks (assert_dimension_rank_1): lbound")
    call assert(all(ubound(a) == ubound(b)), "checks (assert_dimension_rank_1): ubound")
end subroutine assert_dimension_rank_1

pure subroutine assert_dimension_rank_2(a, b)
    class(*), intent(in) :: a(:, :), b(:, :)
    
    call assert(all(lbound(a) == lbound(b)), "checks (assert_dimension_rank_2): lbound")
    call assert(all(ubound(a) == ubound(b)), "checks (assert_dimension_rank_2): ubound")
end subroutine assert_dimension_rank_2

pure subroutine assert_dimension_rank_3(a, b)
    class(*), intent(in) :: a(:, :, :), b(:, :, :)
    
    call assert(all(lbound(a) == lbound(b)), "checks (assert_dimension_rank_3): lbound")
    call assert(all(ubound(a) == ubound(b)), "checks (assert_dimension_rank_3): ubound")
end subroutine assert_dimension_rank_3

! check
! -----

! If `condition` is `.false.`, then print a message and increment `rc`.
! If `(rc /= RC_SUCCESS)` later, computation will stop.
! This is used for input validation and other non-assertion checks.
! Making `rc` increment is useful to string up multiple `check`s
! without adding too much logic. The details of the `check` are
! printed or logged, so `rc` does not need to be meaningful beyond
! pass/fail.

subroutine impure_check(condition, message, rc)
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    logical, intent(in)          :: condition ! condition to check
    character(len=*), intent(in) :: message   ! error message to print if `condition` is `.false.`
    integer, intent(in out)      :: rc        ! number of errors encountered
    
    if (.not. condition) then
        write(unit=ERROR_UNIT, fmt="(a)") message
        
        rc = rc + 1
    end if
    
    call assert(rc >= 0, "checks (impure_check): negative rc should be impossible")
end subroutine impure_check

pure subroutine pure_check(condition, message, rc, logger)
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    logical, intent(in)                            :: condition ! condition to check
    character(len=*), intent(in)                   :: message   ! error message to print if `condition` is `.false.`
    integer, intent(in out)                        :: rc        ! number of errors encountered
    type(logger_type), intent(in out), allocatable :: logger
    
    if (.not. condition) then
        call write(ERROR_UNIT, message, logger)
        
        rc = rc + 1
    end if
    
    call assert(rc >= 0, "checks (pure_check): negative rc should be impossible")
end subroutine pure_check

! Pure logger
! -----------

! I would make these type-bound procedures, but it seems you can't make
! type-bound procedures for `allocatable` variables.

pure subroutine print(message, logger)
    use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
    
    character(len=*), intent(in)                   :: message
    type(logger_type), intent(in out), allocatable :: logger
    
    call write(OUTPUT_UNIT, message, logger)
end subroutine print

pure subroutine write(unit, message, logger)
    integer, intent(in)                            :: unit
    character(len=*), intent(in)                   :: message
    type(logger_type), intent(in out), allocatable :: logger
    
    type(logger_type), allocatable :: temp
    
    allocate(temp)
    
    temp%unit    = unit
    temp%message = message
    
    call move_alloc(logger, temp%history)
    call move_alloc(temp,   logger)
end subroutine write

subroutine flush(logger)
    type(logger_type), intent(in out), allocatable :: logger
    
    if (allocated(logger)) then
        call logger_print(logger)
        deallocate(logger)
    end if
end subroutine flush

recursive subroutine logger_print(logger)
    type(logger_type), intent(in), allocatable :: logger
    
    if (allocated(logger%history)) call logger_print(logger%history)
    write(unit=logger%unit, fmt="(a)") logger%message
end subroutine logger_print

end module checks
