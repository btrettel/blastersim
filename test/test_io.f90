! tests for I/O procedures
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_io

use prec, only: WP
use units
use unittest, only: test_results_type
implicit none

type(test_results_type) :: tests

call tests%start_tests("io.nml")

call test_read_springer_namelist(tests)

call tests%end_tests()

contains

subroutine test_read_springer_namelist(tests)
    use prec, only: CL, PI
    use cva, only: cv_system_type, run_config_type
    use port, only: path_join
    use io, only: I_BARREL, read_springer_namelist
    
    type(test_results_type), intent(in out) :: tests
    
    character(len=CL)                 :: input_file, path_array(2)
    type(cv_system_type), allocatable :: sys
    type(run_config_type)             :: config
    integer                           :: rc
    
    integer, parameter :: I_PLUNGER = 2, I_BARREL_ATM  = 3, I_PLUNGER_ATM = 4
    
    type(si_length) :: d_barrel, d_e, d_plunger, l_travel
    type(si_volume) :: vol_dead
    
    path_array(1) = "examples"
    path_array(2) = "springer-example.nml"
    input_file = path_join(path_array)
    
    call read_springer_namelist(input_file, sys, config, rc)
    
    call tests%integer_eq(rc, 0, "test_read_springer_namelist, rc")
    
    d_barrel = sqrt((4.0_WP/PI)*sys%cv(I_BARREL)%csa)
    call tests%real_eq(d_barrel%v%v, 13.0e-3_WP, "test_read_springer_namelist, d_barrel")
    
    d_e = sqrt((4.0_WP/PI)*sys%con(I_PLUNGER, I_BARREL)%a_e)
    call tests%real_eq(d_e%v%v, 5.0e-3_WP, "test_read_springer_namelist, d_e")
    
    d_plunger = sqrt((4.0_WP/PI)*sys%cv(I_PLUNGER)%csa)
    call tests%real_eq(d_plunger%v%v, 30.0e-3_WP, "test_read_springer_namelist, d_plunger")
    
    call tests%real_eq(sys%cv(I_PLUNGER)%k%v%v, 500.0_WP, "test_read_springer_namelist, k")
    call tests%real_eq(sys%cv(I_PLUNGER)%x%v%v, 20.0e-2_WP, "test_read_springer_namelist, l_draw")
    call tests%real_eq(sys%cv(I_PLUNGER)%delta_pre%v%v, 5.0e-2_WP, "test_read_springer_namelist, delta_pre")
    
    l_travel = sys%cv(I_BARREL)%x_stop - sys%cv(I_BARREL)%x
    call tests%real_eq(l_travel%v%v, 50.0e-2_WP, "test_read_springer_namelist, d_plunger")
    
    call tests%real_eq(1.0_WP/sys%cv(I_PLUNGER)%rm_p%v%v, 30.0e-3_WP, "test_read_springer_namelist, m_plunger")
    call tests%real_eq(1.0_WP/sys%cv(I_BARREL)%rm_p%v%v, 1.1e-3_WP, "test_read_springer_namelist, m_proj")
    call tests%real_eq(sys%cv(I_PLUNGER)%m_spring%v%v, 10.0e-3_WP, "test_read_springer_namelist, m_spring")
    call tests%real_eq(sys%cv(I_PLUNGER)%p_fd%v%v, 0.0_WP, "test_read_springer_namelist, p_fd_plunger")
    call tests%real_eq(sys%cv(I_BARREL)%p_fd%v%v, 0.0_WP, "test_read_springer_namelist, p_fd_proj")
    call tests%real_eq(sys%cv(I_PLUNGER)%p_fs%v%v, 0.0_WP, "test_read_springer_namelist, p_fs_plunger")
    call tests%real_eq(sys%cv(I_BARREL)%p_fs%v%v, 0.0_WP, "test_read_springer_namelist, p_fs_proj")
    
    vol_dead = sys%cv(I_BARREL)%x*sys%cv(I_BARREL)%csa
    call tests%real_eq(vol_dead%v%v, 20.0e-6_WP, "test_read_springer_namelist, vol_dead")
    
    call tests%real_eq(sys%con(I_PLUNGER, I_BARREL)%b%v%v, 0.528_WP, "test_read_springer_namelist, b")
    call tests%character_eq(config%id, "springer-example", "test_read_springer_namelist, id")
    call tests%real_eq(sys%cv(I_PLUNGER_ATM)%p_const%v%v, 100000.0_WP, "test_read_springer_namelist, p_atm (1)")
    call tests%real_eq(sys%cv(I_BARREL_ATM)%p_const%v%v, 100000.0_WP, "test_read_springer_namelist, p_atm (2)")
    call tests%real_eq(sys%cv(I_PLUNGER_ATM)%temp_const%v%v, 310.0_WP, "test_read_springer_namelist, temp_atm (1)")
    call tests%real_eq(sys%cv(I_BARREL_ATM)%temp_const%v%v, 310.0_WP, "test_read_springer_namelist, temp_atm (2)")
end subroutine test_read_springer_namelist

end program test_io
