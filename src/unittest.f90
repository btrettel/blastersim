! Module for unit testing procedures.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [flt](https://github.com/btrettel/flt)
! License: [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

module unittest

use prec, only: WP
use nmllog, only: log_type, now, TIMESTAMP_LEN, DEBUG_LEVEL
use timer, only: timer_type
use checks, only: assert
implicit none
private

public :: validate_timestamp

character(len=70), parameter :: LONG_LINE = "----------------------------------------------------------------------"

type, public :: test_results_type
    integer                 :: n_failures = 0
    integer                 :: n_tests    = 0
    type(timer_type)        :: wtime
    type(log_type), pointer :: logger => null()
contains
    procedure :: logical_true => logical_true
    procedure :: logical_false => logical_false
    procedure :: real_eq => real_eq
    procedure :: real_ne => real_ne
    procedure :: real_gt => real_gt
    procedure :: real_lt => real_lt
    procedure :: real_ge => real_ge
    procedure :: real_le => real_le
    generic :: integer_eq => integer5_eq, integer10_eq
    procedure, private :: integer5_eq, integer10_eq
    procedure :: integer_ne => integer_ne
    generic :: integer_ge => integer5_ge, integer10_ge
    procedure, private :: integer5_ge, integer10_ge
    generic :: integer_le => integer5_le, integer10_le
    procedure, private :: integer5_le, integer10_le
    ! TODO: integer_gt, integer_lt
    procedure :: character_eq => character_eq
    ! TODO: character_ne
    procedure :: exit_code_eq => exit_code_eq
    procedure :: exit_code_ne => exit_code_ne
    procedure :: start_tests => start_tests
    procedure :: end_tests => end_tests
end type test_results_type

contains

subroutine logical_true(tests, condition, message_in)
    ! Check whether test `condition` is `.true.`, increase `number_of_failures` if `.false.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    logical, intent(in)          :: condition
    character(len=*), intent(in) :: message_in
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=7)              :: variable_type
    logical                       :: test_passes, returned_logical, compared_logical
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_passes, message, returned_logical, compared_logical
    
    timestamp        = now()
    variable_type    = "logical"
    message          = message_in
    returned_logical = condition
    compared_logical = .true.
    
    test_passes = condition
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (logical_true): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (logical_true): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (logical_true): number of failures exceeds number of tests")
end subroutine logical_true

subroutine logical_false(tests, condition, message_in)
    ! Check whether test `condition` is `.false.`, increase `number_of_failures` if `.true.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    logical, intent(in)          :: condition
    character(len=*), intent(in) :: message_in
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=7)              :: variable_type
    logical                       :: test_passes, returned_logical, compared_logical
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_passes, message, returned_logical, compared_logical
    
    timestamp        = now()
    variable_type    = "logical"
    message          = message_in
    returned_logical = condition
    compared_logical = .false.
    
    test_passes = .not. condition
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (logical_false): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (logical_false): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (logical_false): number of failures exceeds number of tests")
end subroutine logical_false

subroutine real_eq(tests, returned_real, compared_real, message_in, abs_tol)
    ! Check whether two reals are close, increase `number_of_failures` if `.false.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    use checks, only: abs_tolerance, is_close
    
    class(test_results_type), intent(in out) :: tests
    
    real(WP), intent(in)         :: returned_real, compared_real
    character(len=*), intent(in) :: message_in
    
    real(WP), intent(in), optional :: abs_tol
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=4)              :: variable_type
    character(len=2)              :: test_operator
    logical                       :: test_passes
    character(len=:), allocatable :: message
    real(WP)                      :: tolerance, difference
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, &
                            returned_real, compared_real, tolerance, difference
    
    difference = abs(returned_real - compared_real)
    call assert(difference >= 0.0_WP, "unittest (real_eq): difference is negative")
    
    test_operator = "=="
    timestamp     = now()
    variable_type = "real"
    message       = message_in
    
    if (present(abs_tol)) then
        test_passes = is_close(returned_real, compared_real, abs_tol=abs_tol)
        tolerance = abs_tol
    else
        test_passes = is_close(returned_real, compared_real)
        tolerance = abs_tolerance(returned_real, compared_real)
    end if
    
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "real returned = ", returned_real
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "real expected = ", compared_real
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "    tolerance = ", tolerance
            if (is_close(abs(compared_real), 0.0_WP)) then
                write(unit=ERROR_UNIT, fmt="(a, g0.8)") "   difference = ", difference
            else
                write(unit=ERROR_UNIT, fmt="(a, g0.8, a, f0.3, a)") "   difference = ", difference, &
                                                            " (", 100.0_WP * difference / abs(compared_real), "%)"
            end if
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (real_eq): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (real_eq): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (real_eq): number of failures exceeds number of tests")
end subroutine real_eq

subroutine real_ne(tests, returned_real, compared_real, message_in, abs_tol)
    ! Check whether two reals are not close, increase `number_of_failures` if `.true.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    use checks, only: abs_tolerance, is_close
    
    class(test_results_type), intent(in out) :: tests
    
    real(WP), intent(in)         :: returned_real, compared_real
    character(len=*), intent(in) :: message_in
    
    real(WP), intent(in), optional :: abs_tol
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=4)              :: variable_type
    character(len=2)              :: test_operator
    logical                       :: test_passes
    character(len=:), allocatable :: message
    real(WP)                      :: tolerance, difference
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, &
                            returned_real, compared_real, tolerance, difference
    
    difference = abs(returned_real - compared_real)
    call assert(difference >= 0.0_WP, "unittest (real_ne): difference is negative")
    
    test_operator = "/="
    timestamp     = now()
    variable_type = "real"
    message       = message_in
    
    if (present(abs_tol)) then
        test_passes = .not. is_close(returned_real, compared_real, abs_tol=abs_tol)
        tolerance = abs_tol
    else
        test_passes = .not. is_close(returned_real, compared_real)
        tolerance = abs_tolerance(returned_real, compared_real)
    end if
    
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "real returned = ", returned_real
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "      /= real = ", compared_real
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "    tolerance = ", tolerance
            if (is_close(abs(compared_real), 0.0_WP)) then
                write(unit=ERROR_UNIT, fmt="(a, g0.8)") "   difference = ", difference
            else
                write(unit=ERROR_UNIT, fmt="(a, g0.8, a, f0.3, a)") "   difference = ", difference, &
                                                            " (", 100.0_WP * difference / abs(compared_real), "%)"
            end if
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (real_ne): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (real_ne): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (real_ne): number of failures exceeds number of tests")
end subroutine real_ne

subroutine real_gt(tests, returned_real, compared_real, message_in)
    ! Check whether `returned_real > compared_real`, increase `number_of_failures` if `.true.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    real(WP), intent(in)         :: returned_real, compared_real
    character(len=*), intent(in) :: message_in
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=4)              :: variable_type
    character(len=1)              :: test_operator
    logical                       :: test_passes
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, &
                            returned_real, compared_real
    
    test_operator = ">"
    timestamp     = now()
    variable_type = "real"
    message       = message_in
    
    test_passes = returned_real > compared_real
    
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "real returned = ", returned_real
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "       > real = ", compared_real
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (real_gt): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (real_gt): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (real_gt): number of failures exceeds number of tests")
end subroutine real_gt

subroutine real_lt(tests, returned_real, compared_real, message_in)
    ! Check whether `returned_real < compared_real`, increase `number_of_failures` if `.true.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    real(WP), intent(in)         :: returned_real, compared_real
    character(len=*), intent(in) :: message_in
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=4)              :: variable_type
    character(len=1)              :: test_operator
    logical                       :: test_passes
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, &
                            returned_real, compared_real
    
    test_operator = "<"
    timestamp     = now()
    variable_type = "real"
    message       = message_in
    
    test_passes = returned_real < compared_real
    
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "real returned = ", returned_real
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "       < real = ", compared_real
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (real_lt): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (real_lt): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (real_lt): number of failures exceeds number of tests")
end subroutine real_lt

subroutine real_ge(tests, returned_real, compared_real, message_in)
    ! Check whether `returned_real >= compared_real`, increase `number_of_failures` if `.true.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    real(WP), intent(in)         :: returned_real, compared_real
    character(len=*), intent(in) :: message_in
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=4)              :: variable_type
    character(len=2)              :: test_operator
    logical                       :: test_passes
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, &
                            returned_real, compared_real
    
    test_operator = ">="
    timestamp     = now()
    variable_type = "real"
    message       = message_in
    
    test_passes = returned_real >= compared_real
    
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "real returned = ", returned_real
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "      >= real = ", compared_real
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (real_ge): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (real_ge): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (real_ge): number of failures exceeds number of tests")
end subroutine real_ge

subroutine real_le(tests, returned_real, compared_real, message_in)
    ! Check whether `returned_real <= compared_real`, increase `number_of_failures` if `.true.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    real(WP), intent(in)         :: returned_real, compared_real
    character(len=*), intent(in) :: message_in
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=4)              :: variable_type
    character(len=2)              :: test_operator
    logical                       :: test_passes
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, &
                            returned_real, compared_real
    
    test_operator = "<="
    timestamp     = now()
    variable_type = "real"
    message       = message_in
    
    test_passes = returned_real <= compared_real
    
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "real returned = ", returned_real
            write(unit=ERROR_UNIT, fmt="(a, g0.8)") "      >= real = ", compared_real
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (real_le): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (real_le): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (real_le): number of failures exceeds number of tests")
end subroutine real_le

subroutine integer5_eq(tests, returned_integer, compared_integer, message_in)
    ! Check whether two integers are identical, increase `number_of_failures` if `.false.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    integer, intent(in)          :: returned_integer, compared_integer
    character(len=*), intent(in) :: message_in
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=7)              :: variable_type
    character(len=2)              :: test_operator
    logical                       :: test_passes
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, returned_integer, compared_integer
    
    timestamp     = now()
    variable_type = "integer"
    message       = message_in
    test_operator = "=="
    
    test_passes = (returned_integer == compared_integer)
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, i7)") "integer returned = ", returned_integer
            write(unit=ERROR_UNIT, fmt="(a, i7)") "integer expected = ", compared_integer
            write(unit=ERROR_UNIT, fmt="(a, i7)") "      difference = ", abs(returned_integer - compared_integer)
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (integer5_eq): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (integer5_eq): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (integer5_eq): number of failures exceeds number of tests")
end subroutine integer5_eq

subroutine integer10_eq(tests, returned_integer, compared_integer, message_in)
    ! Check whether two integers are identical, increase `number_of_failures` if `.false.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    use prec, only: I10
    
    class(test_results_type), intent(in out) :: tests
    
    integer(I10), intent(in)     :: returned_integer, compared_integer
    character(len=*), intent(in) :: message_in
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=7)              :: variable_type
    character(len=2)              :: test_operator
    logical                       :: test_passes
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, returned_integer, compared_integer
    
    timestamp     = now()
    variable_type = "integer"
    message       = message_in
    test_operator = "=="
    
    test_passes = (returned_integer == compared_integer)
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, i11)") "integer returned = ", returned_integer
            write(unit=ERROR_UNIT, fmt="(a, i11)") "integer expected = ", compared_integer
            write(unit=ERROR_UNIT, fmt="(a, i11)") "      difference = ", abs(returned_integer - compared_integer)
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (integer10_eq): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (integer10_eq): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (integer10_eq): number of failures exceeds number of tests")
end subroutine integer10_eq

subroutine integer_ne(tests, returned_integer, compared_integer, message_in)
    ! Check whether `returned_integer` equal than `compared_integer`, increase `number_of_failures` if `.true.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    integer, intent(in)          :: returned_integer, compared_integer
    character(len=*), intent(in) :: message_in
    
    logical :: test_passes
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=7)              :: variable_type
    character(len=2)              :: test_operator
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, returned_integer, compared_integer
    
    timestamp     = now()
    variable_type = "integer"
    test_operator = "/="
    message       = message_in
    
    test_passes = (returned_integer /= compared_integer)
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, i7)") "integer returned = ", returned_integer
            write(unit=ERROR_UNIT, fmt="(a, i7)") "      /= integer = ", compared_integer
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (integer_ne): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (integer_ne): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (integer_ne): number of failures exceeds number of tests")
end subroutine integer_ne

subroutine integer5_ge(tests, returned_integer, compared_integer, message_in)
    ! Check whether `returned_integer >= compared_integer`, increase `number_of_failures` if `.false.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    integer, intent(in)          :: returned_integer, compared_integer
    character(len=*), intent(in) :: message_in
    
    logical :: test_passes
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=7)              :: variable_type
    character(len=2)              :: test_operator
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, returned_integer, compared_integer
    
    timestamp     = now()
    variable_type = "integer"
    test_operator = ">="
    message       = message_in
    
    test_passes = (returned_integer >= compared_integer)
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, i7)") "integer returned = ", returned_integer
            write(unit=ERROR_UNIT, fmt="(a, i7)") "      >= integer = ", compared_integer
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (integer5_ge): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (integer5_ge): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (integer5_ge): number of failures exceeds number of tests")
end subroutine integer5_ge

subroutine integer10_ge(tests, returned_integer, compared_integer, message_in)
    ! Check whether `returned_integer >= compared_integer`, increase `number_of_failures` if `.false.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    use prec, only: I10
    
    class(test_results_type), intent(in out) :: tests
    
    integer(I10), intent(in)     :: returned_integer, compared_integer
    character(len=*), intent(in) :: message_in
    
    logical :: test_passes
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=7)              :: variable_type
    character(len=2)              :: test_operator
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, returned_integer, compared_integer
    
    timestamp     = now()
    variable_type = "integer"
    test_operator = ">="
    message       = message_in
    
    test_passes = (returned_integer >= compared_integer)
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, i11)") "integer returned = ", returned_integer
            write(unit=ERROR_UNIT, fmt="(a, i11)") "      >= integer = ", compared_integer
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (integer10_ge): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (integer10_ge): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (integer10_ge): number of failures exceeds number of tests")
end subroutine integer10_ge

subroutine integer5_le(tests, returned_integer, compared_integer, message_in)
    ! Check whether `returned_integer <= compared_integer`, increase `number_of_failures` if `.false.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    integer, intent(in)          :: returned_integer, compared_integer
    character(len=*), intent(in) :: message_in
    
    logical :: test_passes
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=7)              :: variable_type
    character(len=2)              :: test_operator
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, returned_integer, compared_integer
    
    timestamp     = now()
    variable_type = "integer"
    test_operator = "<="
    message       = message_in
    
    test_passes = (returned_integer <= compared_integer)
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, i7)") "integer returned = ", returned_integer
            write(unit=ERROR_UNIT, fmt="(a, i7)") "      <= integer = ", compared_integer
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (integer5_le): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (integer5_le): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (integer5_le): number of failures exceeds number of tests")
end subroutine integer5_le

subroutine integer10_le(tests, returned_integer, compared_integer, message_in)
    ! Check whether `returned_integer <= compared_integer`, increase `number_of_failures` if `.false.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    use prec, only: I10
    
    class(test_results_type), intent(in out) :: tests
    
    integer(I10), intent(in)     :: returned_integer, compared_integer
    character(len=*), intent(in) :: message_in
    
    logical :: test_passes
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=7)              :: variable_type
    character(len=2)              :: test_operator
    character(len=:), allocatable :: message
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, returned_integer, compared_integer
    
    timestamp     = now()
    variable_type = "integer"
    test_operator = "<="
    message       = message_in
    
    test_passes = (returned_integer <= compared_integer)
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a, i7)") "integer returned = ", returned_integer
            write(unit=ERROR_UNIT, fmt="(a, i7)") "      <= integer = ", compared_integer
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (integer10_le): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (integer10_le): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (integer10_le): number of failures exceeds number of tests")
end subroutine integer10_le

subroutine character_eq(tests, returned_character_in, compared_character_in, message_in)
    ! Check whether two character variables are identical, increase `number_of_failures` if `.false.`.
    
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    
    class(test_results_type), intent(in out) :: tests
    
    character(len=*), intent(in) :: returned_character_in, compared_character_in
    character(len=*), intent(in) :: message_in
    
    logical :: test_passes
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=9)              :: variable_type
    character(len=2)              :: test_operator
    character(len=:), allocatable :: message, returned_character, compared_character
    
    namelist /test_result/ timestamp, variable_type, test_operator, test_passes, message, returned_character, compared_character
    
    timestamp          = now()
    variable_type      = "character"
    test_operator      = "=="
    message            = message_in
    returned_character = trim(returned_character_in)
    compared_character = trim(compared_character_in)
    
    ! It's not clear to me why `returned_character` needs to be trimmed.
    ! The namelist file will have a lot of spaces otherwise.
    
    test_passes = (trim(returned_character) == trim(compared_character))
    write(unit=tests%logger%unit, nml=test_result)
    
    if (.not. test_passes) then
        tests%n_failures = tests%n_failures + 1
        
        if (DEBUG_LEVEL >= tests%logger%stdout_level) then
            write(unit=ERROR_UNIT, fmt="(a)") "character returned = " // returned_character
            write(unit=ERROR_UNIT, fmt="(a)") "character expected = " // compared_character
            write(unit=ERROR_UNIT, fmt="(a, a, a)") "fail: ", message, new_line("a")
        end if
    end if
    
    tests%n_tests = tests%n_tests + 1
    
    call assert(tests%n_tests >= 0, "unittest (character_eq): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (character_eq): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (character_eq): number of failures exceeds number of tests")
end subroutine character_eq

subroutine start_tests(tests, logger)
    use checks, only: assert
    use nmllog, only: DEBUG_LEVEL
    
    class(test_results_type), intent(out) :: tests
    
    type(log_type), intent(in out), target :: logger
    
    logical :: unit_open
    
    inquire(unit=logger%unit, opened=unit_open)
    call assert(unit_open, "unittest (start_tests): logger unit must be open to start tests")
    
    call tests%wtime%start()
    
    tests%logger              => logger
    tests%logger%stdout_level =  DEBUG_LEVEL
    tests%logger%file_level   =  DEBUG_LEVEL
end subroutine start_tests

subroutine end_tests(tests)
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    use checks, only: assert
    
    class(test_results_type), intent(in out) :: tests
    
    integer  :: n_tests, n_failures
    real(WP) :: duration ! in seconds
    logical  :: unit_open
    
    namelist /tests_summary/ n_tests, n_failures, duration
    
    call assert(tests%n_tests >= 0, "unittest (end_tests): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (end_tests): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (end_tests): number of failures exceeds number of tests")
    call assert(associated(tests%logger), "unittest (end_tests): logger is not associated")
    inquire(unit=tests%logger%unit, opened=unit_open)
    call assert(unit_open, "unittest (end_tests): logger unit must be open to end tests")
    
    call tests%wtime%stop()
    duration      = tests%wtime%read()
    n_tests       = tests%n_tests
    n_failures    = tests%n_failures
    write(unit=tests%logger%unit, nml=tests_summary)
    
    write(unit=*, fmt="(a, i0, a, f0.3, a)") "Ran ", tests%n_tests, " tests in ", duration, "s"
    
    if (tests%n_failures /= 0) then
        write(unit=ERROR_UNIT, fmt="(a, i0, a)") "FAILED (failures=", tests%n_failures, ")"
        write(unit=ERROR_UNIT, fmt="(a)") LONG_LINE
        call tests%logger%close()
        stop 1
    else
        write(unit=*, fmt="(a, a)") "All tests passed.", new_line("a")
        write(unit=*, fmt="(a)") "OK"
        write(unit=*, fmt="(a)") LONG_LINE
    end if
    
    nullify(tests%logger)
end subroutine end_tests

subroutine validate_timestamp(tests, timestamp, message)
    type(test_results_type), intent(in out) :: tests
    
    character(len=*), intent(in) :: timestamp, message
    integer :: year, month, day, hour, minutes, seconds, milliseconds, timezone_hour, timezone_minutes, rc_read
    
    ! TODO: If any of these fail, write the timestamp to the log.
    
    call tests%integer_eq(len(trim(timestamp)), TIMESTAMP_LEN, trim(message) // ", timestamp is correct length")
    
    ! There will be a run-time error if the timestamp is too short, so don't do these tests in that case.
    ! The "timestamp is correct length" test will fail, alerting you to that problem if it occurs.
    if (len(trim(timestamp)) >= TIMESTAMP_LEN) then
        read(timestamp(1:4), "(i4)", iostat=rc_read) year
        call tests%integer_eq(rc_read, 0, trim(message) // ", year is numeric")
        if (rc_read == 0) then
            call tests%integer_ge(year, 2000, trim(message) // ", year lower bound")
            call tests%integer_le(year, 2050, trim(message) // ", year upper bound")
        end if
        
        call tests%character_eq(timestamp(5:5), "-", trim(message) // ", timestamp(5)")
        
        read(timestamp(6:7), "(i2)", iostat=rc_read) month
        call tests%integer_eq(rc_read, 0, trim(message) // ", month is numeric")
        if (rc_read == 0) then
            call tests%integer_ge(month, 1, trim(message) // ", month lower bound")
            call tests%integer_le(month, 12, trim(message) // ", month upper bound")
        end if
        
        call tests%character_eq(timestamp(8:8), "-", trim(message) // ", timestamp(8)")
        
        read(timestamp(9:10), "(i2)", iostat=rc_read) day
        call tests%integer_eq(rc_read, 0, trim(message) // ", day is numeric")
        if (rc_read == 0) then
            call tests%integer_ge(day, 1, trim(message) // ", day lower bound")
            call tests%integer_le(day, 31, trim(message) // ", day upper bound")
        end if
        
        call tests%character_eq(timestamp(11:11), "T", trim(message) // ", timestamp(11)")
        
        read(timestamp(12:13), "(i2)", iostat=rc_read) hour
        call tests%integer_eq(rc_read, 0, trim(message) // ", hour is numeric")
        if (rc_read == 0) then
            call tests%integer_ge(hour, 0, trim(message) // ", hour lower bound")
            call tests%integer_le(hour, 23, trim(message) // ", hour upper bound")
        end if
        
        call tests%character_eq(timestamp(14:14), ":", trim(message) // ", timestamp(14)")
        
        read(timestamp(15:16), "(i2)", iostat=rc_read) minutes
        call tests%integer_eq(rc_read, 0, trim(message) // ", minutes are numeric")
        if (rc_read == 0) then
            call tests%integer_ge(minutes, 0, trim(message) // ", minutes lower bound")
            call tests%integer_le(minutes, 59, trim(message) // ", minutes upper bound")
        end if
        
        call tests%character_eq(timestamp(17:17), ":", trim(message) // ", timestamp(17)")
        
        read(timestamp(18:19), "(i2)", iostat=rc_read) seconds
        call tests%integer_eq(rc_read, 0, trim(message) // ", seconds are numeric")
        if (rc_read == 0) then
            call tests%integer_ge(seconds, 0, trim(message) // ", seconds lower bound")
            call tests%integer_le(seconds, 59, trim(message) // ", seconds upper bound")
        end if
        
        call tests%character_eq(timestamp(20:20), ".", trim(message) // ", timestamp(20)")
        
        read(timestamp(21:23), "(i3)", iostat=rc_read) milliseconds
        call tests%integer_eq(rc_read, 0, trim(message) // ", milliseconds are numeric")
        if (rc_read == 0) then
            call tests%integer_ge(milliseconds, 0, trim(message) // ", milliseconds lower bound")
            
            ! tests is commented out as I guess I can't make milliseconds too high.
!            call tests%integer_le(milliseconds, 999, trim(message) // ", milliseconds upper bound")
        end if
        
        !call tests%character_eq(timestamp(24:24), "-", trim(message) // ", timestamp(24)")
        call tests%logical_true((timestamp(24:24) == "-") .or. (timestamp(24:24) == "+"), trim(message) // ", timestamp(24)")
        
        read(timestamp(25:26), "(i2)", iostat=rc_read) timezone_hour
        call tests%integer_eq(rc_read, 0, trim(message) // ", timezone_hour is numeric")
        if (rc_read == 0) then
            call tests%integer_ge(timezone_hour, 0, trim(message) // ", timezone_hour lower bound")
            call tests%integer_le(timezone_hour, 23, trim(message) // ", timezone_hour upper bound")
        end if
        
        call tests%character_eq(timestamp(27:27), ":", trim(message) // ", timestamp(27)")
        
        read(timestamp(28:29), "(i2)", iostat=rc_read) timezone_minutes
        call tests%integer_eq(rc_read, 0, trim(message) // ", timezone_minutes is numeric")
        if (rc_read == 0) then
            call tests%integer_ge(timezone_minutes, 0, trim(message) // ", timezone_minutes lower bound")
            call tests%integer_le(timezone_minutes, 59, trim(message) // ", timezone_minutes upper bound")
        end if
    end if
    
    call assert(tests%n_tests >= 0, "unittest (validate_timestamp): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (validate_timestamp): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (validate_timestamp): number of failures exceeds number of tests")
end subroutine validate_timestamp

subroutine exit_code_eq(tests, command, compared_exit_code, message, output_file, keep_file)
    class(test_results_type), intent(in out) :: tests
    integer, intent(in)                      :: compared_exit_code
    character(len=*), intent(in)             :: command, message, output_file
    logical, intent(in), optional            :: keep_file
    
    integer :: rc_fail, failure_unit
    logical :: keep_file_
    
    if (present(keep_file)) then
        keep_file_ = keep_file
    else
        keep_file_ = .false.
    end if
    
    ! TODO: this would need to be updated for Windows as I'd add .exe to the executable filename there.
    ! TODO: tests also assumes Bash.
    call execute_command_line(command // " 2> " // output_file, exitstat=rc_fail)
    call tests%integer_eq(rc_fail, compared_exit_code, message)

    if (.not. keep_file_) then
        ! Delete output file if test succeeded, otherwise keep it for examination.
        if (rc_fail == compared_exit_code) then
            open(newunit=failure_unit, file=output_file, status="old", action="read")
            close(unit=failure_unit, status="delete")
        end if
    end if
    
    call assert(tests%n_tests >= 0, "unittest (exit_code_eq): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (exit_code_eq): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (exit_code_eq): number of failures exceeds number of tests")
end subroutine exit_code_eq

subroutine exit_code_ne(tests, command, compared_exit_code, message, output_file, keep_file)
    class(test_results_type), intent(in out) :: tests
    integer, intent(in)                      :: compared_exit_code
    character(len=*), intent(in)             :: command, message, output_file
    logical, intent(in), optional            :: keep_file
    
    integer :: rc_fail, failure_unit
    logical :: keep_file_
    
    if (present(keep_file)) then
        keep_file_ = keep_file
    else
        keep_file_ = .false.
    end if
    
    ! TODO: this would need to be updated for Windows as I'd add .exe to the executable filename there.
    ! TODO: tests also assumes Bash.
    call execute_command_line(command // " 2> " // output_file, exitstat=rc_fail)
    call tests%integer_ne(rc_fail, compared_exit_code, message)

    if (.not. keep_file_) then
        ! Delete output file if test succeeded, otherwise keep it for examination.
        if (rc_fail /= compared_exit_code) then
            open(newunit=failure_unit, file=output_file, status="old", action="read")
            close(unit=failure_unit, status="delete")
        end if
    end if
    
    call assert(tests%n_tests >= 0, "unittest (exit_code_ne): negative number of tests")
    call assert(tests%n_failures >= 0, "unittest (exit_code_ne): negative number of failures")
    call assert(tests%n_failures <= tests%n_tests, "unittest (exit_code_ne): number of failures exceeds number of tests")
end subroutine exit_code_ne

end module unittest
