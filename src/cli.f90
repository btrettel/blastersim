! Module for a simple command-line interface.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module cli

implicit none
private

public :: get_input_file_name_from_cli

contains

subroutine get_input_file_name_from_cli(prog, input_file_name, extra)
    use, intrinsic :: iso_fortran_env, only: compiler_options, compiler_version, error_unit
    use prec, only: CL
    use build, only: DEBUG
    use rev, only: TAG, REVISION_DATE, MODIFIED
    use checks, only: assert
    use stopcodes, only: EX_OK, EX_NOINPUT

    character(len=*), intent(in)   :: prog
    character(len=CL), intent(out) :: input_file_name
    
    character(len=CL), intent(in), optional :: extra
    
    character(len=CL) :: modified_string
    logical :: input_file_exists
    
    call get_command_argument(1, value=input_file_name)

    if ((len(trim(input_file_name)) == 0) .or. &
            (trim(input_file_name) == "--help") .or. &
            (trim(input_file_name) == "-help") .or. &
            (trim(input_file_name) == "-h") .or. &
            (trim(input_file_name) == "--version") .or. &
            (trim(input_file_name) == "-version") .or. &
            (trim(input_file_name) == "-v")) then
        write(unit=*, fmt="(3a)") "Usage: ", prog, " FILENAME"
        
        if (MODIFIED) then
           modified_string = ", modified"
        else
           modified_string = "" 
        end if
        
        write(unit=*, fmt="(a)") "Version: " // TAG // " (" // REVISION_DATE // trim(modified_string) // ")"
        
        if (DEBUG) then
            write(unit=*, fmt="(a)") "Build: debug"
        else
            write(unit=*, fmt="(a)") "Build: release"
        end if
        
        write(unit=*, fmt="(a, a)") "Compiler: ", compiler_version()
        write(unit=*, fmt="(a, a)") "Compiler flags: ", compiler_options()
        
        if (present(extra)) then
            write(unit=*, fmt="(a)") trim(extra)
        end if
        
        stop EX_OK, quiet=.true.
    end if
    
    inquire(file=input_file_name, exist=input_file_exists)
    
    if (.not. input_file_exists) then
        write(unit=error_unit, fmt="(3a)") 'ERROR: Input file "', trim(input_file_name), '" does not exist.'
        stop EX_NOINPUT, quiet=.true.
    end if
end subroutine get_input_file_name_from_cli

end module cli
