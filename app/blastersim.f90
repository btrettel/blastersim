! BlasterSim: a spring and pneumatic blaster simulator
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program blastersim

use, intrinsic :: iso_fortran_env, only: IOSTAT_END, ERROR_UNIT, OUTPUT_UNIT
use prec, only: CL
use cli, only: get_input_file_name_from_cli
use io, only: I_BARREL, read_pneumatic_namelist, read_springer_namelist
use cva, only: run_config_type, cv_system_type, run_status_type, T_STOP_DEFAULT, SUCCESS_RUN_RC, TIMEOUT_RUN_RC, run
use stopcodes, only: EX_OK, EX_USAGE
use rev, only: TAG, REVISION_DATE, MODIFIED
implicit none

character(len=CL)                 :: input_file, extra, modified_string
type(run_config_type)             :: config
type(cv_system_type), allocatable :: sys_start, sys_end
integer                           :: rc
type(run_status_type)             :: status

extra = "<http://github.com/btrettel/blastersim/>" // new_line("a") // "Written by Ben Trettel."

call get_input_file_name_from_cli("blastersim", input_file, extra=extra)

if (MODIFIED) then
   modified_string = ", modified"
else
   modified_string = "" 
end if

write(unit=*, fmt="(a)") "BlasterSim " // TAG // " (" // REVISION_DATE // trim(modified_string) // ")"
write(unit=*, fmt="(a)") "Running simulation..."

call read_pneumatic_namelist(trim(input_file), sys_start, config, rc)
if (rc /= 0) then
    if (rc /= IOSTAT_END) stop EX_USAGE, quiet=.true.

    call read_springer_namelist(trim(input_file), sys_start, config, rc)
    if (rc /= 0) then
        if (rc == IOSTAT_END) then
            write(unit=ERROR_UNIT, fmt="(a)") "ERROR: Empty input file? No pneumatic or springer namelists detected."
        end if
        stop EX_USAGE, quiet=.true.
    end if
end if

call run(config, sys_start, sys_end, status)

!tripwire$ begin 640DA585 Update `\secref{run-time-checks}` of verval.tex.
if (status%rc == SUCCESS_RUN_RC) then
    write(unit=OUTPUT_UNIT, fmt="(a)") "SUCCESS!"
    write(unit=OUTPUT_UNIT, fmt="(a, f0.2, a)") "muzzle velocity: ", sys_end%cv(I_BARREL)%x_dot%v%v, " m/s"
    stop EX_OK, quiet=.true.
else
    write(unit=ERROR_UNIT, fmt="(a, i0)") "ERROR: return code ", status%rc
    select case (status%rc)
        case (TIMEOUT_RUN_RC)
            write(unit=ERROR_UNIT, fmt="(a, f3.1, a)") "Projectile did not leave barrel after ", &
                                                            T_STOP_DEFAULT, " seconds."
        case default
            write(unit=ERROR_UNIT, fmt="(a)") "Refer to BlasterSim User's Guide for meaning of return code."
            write(unit=ERROR_UNIT, fmt="(a)") "<http://trettel.us/blastersim/docs/verification.html#run-time-checks>"
    end select
    
    stop EX_USAGE, quiet=.true.
end if
!tripwire$ end

end program blastersim
