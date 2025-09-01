! Module for wall-clock timers.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module timer

use prec, only: WP
use checks, only: assert
implicit none
private

! <https://j3-fortran.org/doc/year/21/21-117r3.txt>:
! > Many systems have a clock that is too fast to be effectively supported with default integers. Therefore, it is recommended that
! > all invocations of SYSTEM_CLOCK use integer arguments with a decimal exponent range of at least 18.
! > 
! > Using SYSTEM_CLOCK in this manner lets the processor select the most accurate clock available and minimize how often the COUNT
! > value resets to zero.
! <https://fortran-lang.discourse.group/t/proper-usage-of-system-clock/3245/3>
integer, public, parameter :: TIMER_KIND = selected_int_kind(18)

public :: timeit, sleep, now

integer, public, parameter :: TIMESTAMP_LEN = 29

type, public :: timer_type
    logical             :: active = .false. ! whether the timer is currently timing
    integer(TIMER_KIND) :: sum_count = 0_TIMER_KIND, & ! time before current timer start
                           start_count, count_max
    real(WP)            :: count_rate
    
    ! nagfor complains if `count_rate` is different type than `count` or `count_max`:
    ! <https://fortran-lang.discourse.group/t/proper-usage-of-system-clock/3245>
    
    ! However, the timer is more accurate with a `real` `count_rate`:
    ! <https://stackoverflow.com/a/54275654/1124489>
contains
    procedure :: start => start_timer
    procedure :: stop  => stop_timer
    procedure :: reset => reset_timer
    procedure :: read  => read_timer
end type timer_type

abstract interface
    subroutine timeit_subroutine()
    end subroutine timeit_subroutine
end interface

contains

subroutine start_timer(timer)
    use checks, only: assert, is_close
    
    class(timer_type), intent(in out) :: timer
    
    integer(TIMER_KIND) :: start_count, count_max
    real(WP)            :: count_rate
    
    call assert(.not. timer%active, "timer (start_timer): An active timer can not be started.")
    
    timer%active = .true.
    
    call system_clock(start_count, count_rate, count_max)
    
    call assert(count_rate > 0.0_WP, "timer (start_timer): negative count rate")
    call assert(count_max > 0, "timer (start_timer): negative count max")
    
    ! If the timer was previously used, make sure that the `count_rate` and `count_max` match that previously used.
    if (timer%sum_count /= 0_TIMER_KIND) then
        call assert(timer%count_max == count_max, "timer (start_timer): inconsistent count_max")
        call assert(is_close(timer%count_rate, count_rate), "timer (start_timer): inconsistent count rate")
    end if
    
    timer%start_count = start_count
    timer%count_rate  = count_rate
    timer%count_max   = count_max
end subroutine start_timer

subroutine stop_timer(timer)
    use checks, only: assert, is_close
    
    class(timer_type), intent(in out) :: timer
    
    integer(TIMER_KIND) :: stop_count, count_max, before_count_max
    real(WP)            :: count_rate
    
    call assert(timer%active, "timer (stop_timer): An inactive timer has already stopped.")
    
    timer%active = .false.
    
    call system_clock(stop_count, count_rate, count_max)
    
    call assert(count_rate > 0.0_WP, "timer (stop_timer): negative count rate")
    call assert(count_max > 0, "timer (stop_timer): negative count max")
    
    call assert(is_close(timer%count_rate, count_rate), "timer (stop_timer): inconsistent count_max")
    call assert(timer%count_max == count_max, "timer (stop_timer): inconsistent count rate")
    
    ! Detect if the timer wraps.
    if (stop_count < timer%start_count) then
        ! Then the count wrapped, so adjust count appropriately. This is similar to the following Intel Fortran forum post.
        ! <https://community.intel.com/t5/Intel-Fortran-Compiler/how-to-stop-resetting-the-system-clock/m-p/756868#M12355>
        ! It's not clear to me if this is off by one, but I guess that doesn't matter much.
        
        before_count_max = count_max - timer%start_count
        timer%sum_count = timer%sum_count + stop_count + before_count_max
    else
        ! Then the count did not wrap.
        
        timer%sum_count = timer%sum_count + stop_count - timer%start_count
    end if
end subroutine stop_timer

pure subroutine reset_timer(timer)
    class(timer_type), intent(in out) :: timer
    
    call assert(.not. timer%active, "timer (reset_timer): Timer must be inactive to be reset.")
    
    timer%sum_count = 0_TIMER_KIND
end subroutine reset_timer

pure function read_timer(timer)
    class(timer_type), intent(in) :: timer
    
    real(WP) :: read_timer
    
    call assert(.not. timer%active, "timer (read_timer): Timer must be inactive to be read.")
    
    read_timer = real(timer%sum_count, WP) / timer%count_rate
end function read_timer

function timeit(f, number)
    ! Simplified version of <https://docs.python.org/3/library/timeit.html>.
    
    use prec, only: I9
    
    procedure(timeit_subroutine) :: f
    
    integer(I9), optional :: number
    
    real(WP) :: timeit
    
    type(timer_type) :: wtime
    integer(I9)      :: i, number_
    
    if (present(number)) then
        number_ = number
    else
        number_ = 1000000_I9
    end if
    
    call wtime%start()
    do i = 1, number_ ! SERIAL
        ! Since the subroutine doesn't return anything, won't it be optimized out?
        call f()
    end do
    call wtime%stop()
    
    timeit = wtime%read()
end function timeit

subroutine sleep(duration)
    use checks, only: assert
    
    real(WP), intent(in) :: duration
    
    integer(TIMER_KIND) :: start_count, count_max, stop_count, running_count
    real(WP)            :: count_rate
    
    call assert(duration > 0.0_WP, "duration must be greater than zero.")
    
    call system_clock(start_count, count_rate, count_max)
    stop_count = start_count + ceiling(duration * count_rate, TIMER_KIND)
    
    ! Detect if the timer will wrap.
    if (stop_count > count_max) then
        ! It's not clear to me if this is off by one, but I guess that doesn't matter much.
        stop_count = stop_count - count_max
    end if
    
    do ! SERIAL
        call system_clock(running_count, count_rate, count_max)
        
        if (start_count < stop_count) then
            ! Normal case where the clock doesn't wrap.
            if (running_count > stop_count) then
                exit
            end if
        else
            ! If the clock wraps...
            if ((running_count > stop_count) .and. (running_count < start_count)) then
                exit
            end if
        end if
    end do
end subroutine sleep

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
    
    call assert(values(1) > 2000, "timer (now): year is in the past")
    call assert(values(1) < 2100, "timer (now): year is in the future")
    call assert(values(2) >= 1, "timer (now): month is below 1")
    call assert(values(2) <= 12, "timer (now): month is above 12")
    call assert(values(3) >= 1, "timer (now): day is below 1")
    call assert(values(3) <= 31, "timer (now): day is above 31")
    call assert(values(5) >= 0, "timer (now): hour is below 0")
    call assert(values(5) <= 23, "timer (now): hour is above 23")
    call assert(values(6) >= 0, "timer (now): minutes are below 0")
    call assert(values(6) <= 59, "timer (now): minutes are above 23")
    call assert(values(7) >= 0, "timer (now): seconds are below 0")
    call assert(values(7) <= 59, "timer (now): seconds are above 23")
    call assert(values(8) >= 0, "timer (now): milliseconds are below 0")
    call assert(values(8) <= 999, "timer (now): milliseconds are above 999")
end function now

end module timer
