! BlasterSim: a spring and pneumatic blaster simulator
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program blastersim

use, intrinsic :: iso_fortran_env, only: ERROR_UNIT, OUTPUT_UNIT
use prec, only: CL
use cli, only: get_input_file_name_from_cli
use io, only: read_springer_namelist
use cva, only: run_config_type, cv_system_type, run_status_type, SUCCESS_RUN_RC, run
use stopcodes, only: EX_OK, EX_USAGE
implicit none

character(len=CL)                 :: input_file, extra
type(run_config_type)             :: config
type(cv_system_type), allocatable :: sys_start, sys_end
integer                           :: rc
type(run_status_type)             :: status

extra = "<http://github.com/btrettel/blastersim/>" // new_line("a") // "Written by Ben Trettel."

call get_input_file_name_from_cli("blastersim", input_file, extra=extra)

call read_springer_namelist(trim(input_file), sys_start, config, rc)
if (rc /= 0) error stop

call run(config, sys_start, sys_end, status)

if (status%rc == SUCCESS_RUN_RC) then
    write(unit=OUTPUT_UNIT, fmt="(a)") "SUCCESS!"
    write(unit=OUTPUT_UNIT, fmt="(a, f0.2, a)") "muzzle velocity: ", sys_end%cv(4)%x_dot%v%v, " m/s"
    stop EX_OK, quiet=.true.
else
    write(unit=ERROR_UNIT, fmt="(a, i0)") "ERROR: return code ", status%rc
    stop EX_USAGE, quiet=.true.
end if

end program blastersim
