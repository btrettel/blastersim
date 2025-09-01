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

subroutine get_input_file_name_from_cli(prog, input_file_name)
    use, intrinsic :: iso_fortran_env, only: compiler_options, compiler_version, error_unit
    use build, only: DEBUG
    use rev, only: TAG, REVISION_DATE, MODIFIED
    use checks, only: assert
    use stopcodes, only: EX_OK, EX_NOINPUT

    character(len=*), intent(in)               :: prog
    character(len=:), intent(out), allocatable :: input_file_name
    
    character(len=:), allocatable :: modified_string
    integer :: arg_len
    logical :: input_file_exists
    
    call assert(.not. allocated(input_file_name), "cli (get_input_file_name_from_cli): input_file_name already allocated")
    
    call get_command_argument(1, length=arg_len)
    allocate(character(len=arg_len) :: input_file_name)
    call assert(arg_len >= 0, "cli (get_input_file_name_from_cli): arg_len is negative")
    call get_command_argument(1, value=input_file_name)

    if ((len(input_file_name) == 0) .or. &
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
        
        write(unit=*, fmt="(a)") "Version: " // TAG // " (" // REVISION_DATE // modified_string // ")"
        
        if (DEBUG) then
            write(unit=*, fmt="(a)") "Build: debug"
        else
            write(unit=*, fmt="(a)") "Build: release"
        end if
        
        write(unit=*, fmt="(a, a)") "Compiler: ", compiler_version()
        write(unit=*, fmt="(a, a)") "Compiler flags: ", compiler_options()
        
        stop EX_OK
    end if
    call assert(allocated(input_file_name), "cli (get_input_file_name_from_cli): input_file_name should be allocated")
    
    inquire(file=input_file_name, exist=input_file_exists)
    
    if (.not. input_file_exists) then
        write(unit=error_unit, fmt="(3a)") 'ERROR: Input file "', input_file_name, '" does not exist.'
        stop EX_NOINPUT
    end if
end subroutine get_input_file_name_from_cli

end module cli
