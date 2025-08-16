! Module for structured logging in namelist files.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [flt](https://github.com/btrettel/flt)
! License: [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

! Based roughly on the Python logging module.
! <https://docs.python.org/3/howto/logging.html>

module nmllog

implicit none
private

public :: now

integer, public, parameter :: TIMESTAMP_LEN = 29
integer, public, parameter :: UNIT_CLOSED   = -1

! Needed or else ifort truncates namelist file lines and can't read them back correctly.
! 10000 is arbitrary. Just pick a number larger than anything you expect to use.
! gfortran had a strange run-time error with 1000.
integer, public, parameter :: NML_RECL = 10000

integer, public, parameter :: NOT_SET_LEVEL  = 0
integer, public, parameter :: DEBUG_LEVEL    = 10
integer, public, parameter :: INFO_LEVEL     = 20
integer, public, parameter :: WARNING_LEVEL  = 30
integer, public, parameter :: ERROR_LEVEL    = 40
integer, public, parameter :: CRITICAL_LEVEL = 50

integer, public, parameter :: TIME_AND_LEVEL = 0
integer, public, parameter :: LEVEL_ONLY     = 1
integer, public, parameter :: NO_PREFIX      = 2

character(len=*), public, parameter :: NOT_SET_STRING  = "not set"
character(len=*), public, parameter :: DEBUG_STRING    = "debug"
character(len=*), public, parameter :: INFO_STRING     = "info"
character(len=*), public, parameter :: WARNING_STRING  = "warning"
character(len=*), public, parameter :: ERROR_STRING    = "error"
character(len=*), public, parameter :: CRITICAL_STRING = "critical"

type, public :: log_type
    character(len=:), allocatable :: filename
    integer :: unit          = UNIT_CLOSED
    integer :: file_level    = WARNING_LEVEL
    integer :: stdout_level  = WARNING_LEVEL ! also used for stderr
    integer :: stdout_prefix = LEVEL_ONLY
contains
    procedure :: open => log_open
    procedure :: close => log_close
    procedure :: debug => log_debug
    procedure :: info => log_info
    procedure :: warning => log_warning
    procedure :: error => log_error
    procedure :: critical => log_critical
    procedure :: debug_info => log_debug_info
    procedure :: check => impure_check
end type log_type

type, public :: pure_log_type
    type(pure_log_data_type), pointer :: head => null()
    type(pure_log_data_type), pointer :: tail => null()
    type(log_type)                    :: logger
contains
    procedure :: open => pure_log_open
    procedure :: close => pure_log_close
    procedure :: debug => pure_log_debug
    procedure :: info => pure_log_info
    procedure :: warning => pure_log_warning
    procedure :: error => pure_log_error
    procedure :: critical => pure_log_critical
    procedure :: check => pure_check
end type pure_log_type

type, public :: pure_log_data_type
    type(pure_log_data_type), pointer :: next => null()
    integer                           :: level
    character(len=:), allocatable     :: message
    ! `timestamp` not included as no time procedures can be `pure`
end type pure_log_data_type

contains

function now()
    use checks, only: assert
    
    character(len=5)             :: zone
    integer                      :: values(8)
    character(len=4)             :: year
    character(len=2)             :: month, day, hour, minutes, seconds
    character(len=3)             :: milliseconds
    character(len=TIMESTAMP_LEN) :: now
    
    ! ISO 8601 date-time format.
    ! <https://en.wikipedia.org/wiki/ISO_8601>
    
    call date_and_time(zone=zone, values=values)
    
    write(unit=year,         fmt="(i4.4)") values(1)
    write(unit=month,        fmt="(i2.2)") values(2)
    write(unit=day,          fmt="(i2.2)") values(3)
    write(unit=hour,         fmt="(i2.2)") values(5)
    write(unit=minutes,      fmt="(i2.2)") values(6)
    write(unit=seconds,      fmt="(i2.2)") values(7)
    write(unit=milliseconds, fmt="(i3.3)") values(8)
    
    now = year // "-" // month // "-" // day // "T" // hour // ":" // minutes // ":" // seconds // &
                "." // milliseconds // zone(1:3) // ":" // zone(4:5)
    
    call assert(values(1) > 2000, "nmllog (now): year is in the past")
    call assert(values(1) < 2100, "nmllog (now): year is in the future")
    call assert(values(2) >= 1, "nmllog (now): month is below 1")
    call assert(values(2) <= 12, "nmllog (now): month is above 12")
    call assert(values(3) >= 1, "nmllog (now): day is below 1")
    call assert(values(3) <= 31, "nmllog (now): day is above 31")
    call assert(values(5) >= 0, "nmllog (now): hour is below 0")
    call assert(values(5) <= 23, "nmllog (now): hour is above 23")
    call assert(values(6) >= 0, "nmllog (now): minutes are below 0")
    call assert(values(6) <= 59, "nmllog (now): minutes are above 23")
    call assert(values(7) >= 0, "nmllog (now): seconds are below 0")
    call assert(values(7) <= 59, "nmllog (now): seconds are above 23")
    call assert(values(8) >= 0, "nmllog (now): milliseconds are below 0")
    call assert(values(8) <= 999, "nmllog (now): milliseconds are above 999")
end function now

subroutine log_open(this, filename, level, file_level, stdout_level)
    use checks, only: assert
    
    class(log_type), intent(out) :: this
    
    character(len=*), intent(in)  :: filename
    integer, optional, intent(in) :: level, file_level, stdout_level
    
    if (present(level)) then
        this%file_level   = level
        this%stdout_level = level
    end if
    
    if (present(file_level)) then
        this%file_level = level
    end if
    
    if (present(stdout_level)) then
        this%stdout_level = level
    end if
    
    this%filename = filename
    open(newunit=this%unit, &
            action="write", &
            status="replace", &
            position="rewind", &
            file=trim(this%filename), &
            delim="quote", &
            recl=NML_RECL)
    
    call assert(this%file_level >= 0, "nmllog (log_open): file_level is negative")
    call assert(this%stdout_level >= 0, "nmllog (log_open): stdout_level is negative")
end subroutine log_open

subroutine log_close(this)
    use checks, only: assert
    
    class(log_type), intent(in out) :: this
    
    logical :: unit_open
    
    inquire(unit=this%unit, opened=unit_open)
    call assert(unit_open, "nmllog (log_close): unit must be open to close")
    
    close(unit=this%unit)
    this%unit = UNIT_CLOSED
    
    deallocate(this%filename)
end subroutine log_close

subroutine log_writer(this, message_in, level_code)
    use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT, ERROR_UNIT
    use checks, only: assert
    
    type(log_type), intent(in) :: this
    
    character(len=*), intent(in) :: message_in
    integer, intent(in)          :: level_code
    
    integer :: print_unit
    logical :: unit_open
    
    character(len=TIMESTAMP_LEN)  :: timestamp
    character(len=:), allocatable :: level, message
    
    namelist /log/ timestamp, level, message
    
    inquire(unit=this%unit, opened=unit_open)
    call assert(unit_open, "nmllog (log_writer): unit must be open to write")
    
    timestamp = now()
    
    select case (level_code)
        case (DEBUG_LEVEL)
            level = DEBUG_STRING
        case (INFO_LEVEL)
            level = INFO_STRING
        case (WARNING_LEVEL)
            level = WARNING_STRING
        case (ERROR_LEVEL)
            level = ERROR_STRING
        case (CRITICAL_LEVEL)
            level = CRITICAL_STRING
        case default
            level = NOT_SET_STRING
    end select
    
    message = message_in
    if (level_code >= this%file_level) then
        write(unit=this%unit, nml=log)
    end if
    
    if (level_code >= WARNING_LEVEL) then
        print_unit = ERROR_UNIT
    else
        print_unit = OUTPUT_UNIT
    end if
    
    if (level_code >= this%stdout_level) then
        select case (this%stdout_prefix)
            case (TIME_AND_LEVEL)
                write(unit=print_unit, fmt="(5a)") timestamp, " [", level, "] ", message
            case (LEVEL_ONLY)
                write(unit=print_unit, fmt="(4a)") "[", level, "] ", message
            case (NO_PREFIX)
                write(unit=print_unit, fmt="(a)") message
            case default
                error stop "nmllog (log_writer): invalid stdout_prefix"
        end select
    end if
    
    call assert(this%file_level >= 0, "nmllog (log_writer): file_level is negative")
    call assert(this%stdout_level >= 0, "nmllog (log_writer): stdout_level is negative")
end subroutine log_writer

subroutine log_debug(this, message)
    class(log_type), intent(in) :: this
    
    character(len=*), intent(in) :: message
    
    call log_writer(this, message, DEBUG_LEVEL)
end subroutine log_debug

subroutine log_info(this, message)
    class(log_type), intent(in) :: this
    
    character(len=*), intent(in) :: message
    
    call log_writer(this, message, INFO_LEVEL)
end subroutine log_info

subroutine log_warning(this, message)
    class(log_type), intent(in) :: this
    
    character(len=*), intent(in) :: message
    
    call log_writer(this, message, WARNING_LEVEL)
end subroutine log_warning

subroutine log_error(this, message)
    class(log_type), intent(in) :: this
    
    character(len=*), intent(in) :: message
    
    call log_writer(this, message, ERROR_LEVEL)
end subroutine log_error

subroutine log_critical(this, message)
    class(log_type), intent(in) :: this
    
    character(len=*), intent(in) :: message
    
    call log_writer(this, message, CRITICAL_LEVEL)
end subroutine log_critical

subroutine log_debug_info(this)
    use, intrinsic :: iso_fortran_env, only: compiler_opt => compiler_options, compiler_ver => compiler_version
    use, intrinsic :: ieee_arithmetic, only: ieee_support_datatype, ieee_support_denormal, ieee_support_divide, &
                                                ieee_support_inf, ieee_support_nan, ieee_support_sqrt, ieee_support_standard
    use prec, only: WP
    use checks, only: assert
    
    ! Maybe:
    ! - `epsilon`
    ! - `spacing`
    ! - kind code, range, and huge for `I5`
    
    ! TODO: If `DEBUG_LEVEL >= this%stdout_level` then print this information to stdout.
    ! TODO: Print version number, Git revision, revision date, compile date, and whether there are local code modifications.
    ! Similar to <https://github.com/firemodels/fds/blob/master/Source/cons.f90#L278>.
    
    class(log_type), intent(in) :: this
    
    character(len=:), allocatable :: compiler_options, compiler_version, level
    character(len=TIMESTAMP_LEN)  :: timestamp
    logical                       :: unit_open
    
    !real(WP) :: real_huge
    integer  :: real_kind_code, real_precision, real_range, real_radix, &
                    real_min_exponent, real_max_exponent, &
                    integer_kind_code, integer_range, integer_huge
    logical  :: real_support_datatype, real_support_denormal, real_support_divide, &
                    real_support_inf, real_support_nan, real_support_sqrt, real_support_standard
    
    namelist /debug_info/ timestamp, level, compiler_options, compiler_version, &
                            real_kind_code, real_precision, real_range, real_radix, & !real_huge, &
                            real_min_exponent, real_max_exponent, &
                            real_support_datatype, real_support_denormal, real_support_divide, &
                            real_support_inf, real_support_nan, real_support_sqrt, real_support_standard, &
                            integer_kind_code, integer_range, integer_huge
    
    inquire(unit=this%unit, opened=unit_open)
    call assert(unit_open, "nmllog (log_debug_info): unit must be open to write")
    
    timestamp        = now()
    level            = DEBUG_STRING
    compiler_options = compiler_opt()
    compiler_version = compiler_ver()
    
    real_kind_code        = WP
    real_precision        = precision(1.0_WP)
    real_range            = range(1.0_WP)
    real_radix            = radix(1.0_WP)
    !real_huge             = huge(1.0_WP)
    real_min_exponent     = minexponent(1.0_WP)
    real_max_exponent     = maxexponent(1.0_WP)
    real_support_datatype = ieee_support_datatype(1.0_WP)
    real_support_denormal = ieee_support_denormal(1.0_WP)
    real_support_divide   = ieee_support_divide(1.0_WP)
    real_support_inf      = ieee_support_inf(1.0_WP)
    real_support_nan      = ieee_support_nan(1.0_WP)
    real_support_sqrt     = ieee_support_sqrt(1.0_WP)
    real_support_standard = ieee_support_standard(1.0_WP)
    
    integer_kind_code = kind(1)
    integer_range     = range(1)
    integer_huge      = huge(1)
    
    ! `round="down"` is needed to prevent ifort and ifx from getting irritated when reading in `real_huge`
    ! <https://fortran-lang.discourse.group/t/ifort-namelist-bug/5399>
    ! But: `round="down"` is not supported by flang-7.
    ! Later: When flang-7 or nvfortran or whatever supports `round="down"`, uncomment the `real_huge` lines here and in
    ! test_nmllog.f90.
    write(unit=this%unit, nml=debug_info) !, round="down")
end subroutine log_debug_info

! Pure logging
! ------------

subroutine pure_log_open(pure_logger, logger)
    use checks, only: assert
    
    class(pure_log_type), intent(out) :: pure_logger
    type(log_type), intent(in)        :: logger
    
    logical :: unit_open
    
    call logger%debug("Pure logger started.")
    pure_logger%logger = logger
    
    inquire(unit=pure_logger%logger%unit, opened=unit_open)
    call assert(unit_open, "nmllog (pure_log_open): unit must be open to create pure log")
end  subroutine pure_log_open

subroutine pure_log_close(pure_logger)
    class(pure_log_type), intent(in out) :: pure_logger
    
    type(pure_log_data_type), pointer :: current
    
    ! Write and deallocate all messages.
    do ! SERIAL
        if (.not. associated(pure_logger%head)) then
            exit
        end if
        
        call log_writer(pure_logger%logger, pure_logger%head%message // " [pure]", pure_logger%head%level)
        
        current => pure_logger%head
        pure_logger%head => pure_logger%head%next
        deallocate(current)
    end do
    
    nullify(pure_logger%tail)
    
    call pure_logger%logger%debug("Pure logger ended.")
end subroutine pure_log_close

pure subroutine pure_log_writer(pure_logger, message, level)
    type(pure_log_type), intent(in out) :: pure_logger
    
    character(len=*), intent(in) :: message
    integer, intent(in)          :: level
    
    type(pure_log_data_type), pointer :: pure_log_data
    
    ! TODO: Assert something here. Try something about the association status of a pointer?
    
    allocate(pure_log_data)
    pure_log_data%message = message
    pure_log_data%level   = level
    
    ! TODO: Document this better. It's not clear to me what the logic here is.
    if (.not. associated(pure_logger%head)) then
        pure_logger%head => pure_log_data
    else
        pure_logger%tail%next => pure_log_data
    end if
    
    pure_logger%tail => pure_log_data
end subroutine pure_log_writer

pure subroutine pure_log_debug(this, message)
    class(pure_log_type), intent(in out) :: this
    
    character(len=*), intent(in) :: message
    
    call pure_log_writer(this, message, DEBUG_LEVEL)
end subroutine pure_log_debug

pure subroutine pure_log_info(this, message)
    class(pure_log_type), intent(in out) :: this
    
    character(len=*), intent(in) :: message
    
    call pure_log_writer(this, message, INFO_LEVEL)
end subroutine pure_log_info

pure subroutine pure_log_warning(this, message)
    class(pure_log_type), intent(in out) :: this
    
    character(len=*), intent(in) :: message
    
    call pure_log_writer(this, message, WARNING_LEVEL)
end subroutine pure_log_warning

pure subroutine pure_log_error(this, message)
    class(pure_log_type), intent(in out) :: this
    
    character(len=*), intent(in) :: message
    
    call pure_log_writer(this, message, ERROR_LEVEL)
end subroutine pure_log_error

pure subroutine pure_log_critical(this, message)
    class(pure_log_type), intent(in out) :: this
    
    character(len=*), intent(in) :: message
    
    call pure_log_writer(this, message, CRITICAL_LEVEL)
end subroutine pure_log_critical

! check
! -----

! If `condition` is `.false.`, then print a message and increment `rc`.
! If `(rc /= RC_SUCCESS)` later, computation will stop.
! This is used for input validation and other non-assertion checks.
! Making `rc` increment is useful to string up multiple `check`s
! without adding too much logic. The details of the `check` are
! logged, so there's `rc` does not need to be meaningful beyond
! pass/fail.

subroutine impure_check(logger, condition, message, rc)
    use checks, only: assert
    
    class(log_type), intent(in)  :: logger
    logical, intent(in)          :: condition ! condition to check
    character(len=*), intent(in) :: message   ! error message to print if `condition` is `.false.`
    integer, intent(in out)      :: rc        ! number of errors encountered
    
    if (.not. condition) then
        call logger%error(message)
        
        rc = rc + 1
    end if
    
    call assert(rc >= 0, "nmllog (impure_check): negative rc should be impossible")
end subroutine impure_check

pure subroutine pure_check(logger, condition, message, rc)
    use checks, only: assert
    
    class(pure_log_type), intent(in out) :: logger
    logical, intent(in)                  :: condition ! condition to check
    character(len=*), intent(in)         :: message   ! error message to print if `condition` is `.false.`
    integer, intent(in out)              :: rc        ! number of errors encountered
    
    if (.not. condition) then
        call logger%error(message)
        
        rc = rc + 1
    end if
    
    call assert(rc >= 0, "nmllog (pure_check): negative rc should be impossible")
end subroutine pure_check

end module nmllog
