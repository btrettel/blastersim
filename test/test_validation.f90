! validation tests (experimental data)
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_validation

use prec, only: WP
use units
use unittest, only: test_results_type
implicit none

type(test_results_type) :: tests

call tests%start_tests("validation.nml")

call test_2010_08_07(tests)

call tests%end_tests()

contains

subroutine test_2010_08_07(tests)
    use convert
    use gasdata, only: DRY_AIR
    use cva, only: cv_system_type, run_status_type, run
    use prec, only: PI
    use checks, only: assert
    
    type(test_results_type), intent(in out) :: tests

    type(cv_system_type), allocatable :: sys_start, sys_end
    
    type(run_status_type) :: status
    
    type(si_length)          :: d_e, x_1, x_2, d_1, d_2, x_stop_2
    type(si_velocity)        :: x_dot
    type(unitless)           :: y(1)
    type(si_pressure)        :: p_atm, p_1, p_2, p_fs_1, p_fd_1, p_fs_2, p_fd_2
    type(si_temperature)     :: temp_atm
    type(si_area)            :: csa_1, csa_2
    type(si_inverse_mass)    :: rm_p_1, rm_p_2
    type(si_stiffness)       :: k
    type(si_length)          :: x_z
    type(si_volume)          :: vol_1, vol_d
    
    allocate(sys_start)
    allocate(sys_start%cv(2))
    allocate(sys_start%con(2, 2))
    
    sys_start%con(1, 1)%active = .false.
    sys_start%con(2, 1)%active = .false.
    sys_start%con(2, 2)%active = .false.
    
    sys_start%con(1, 2)%active = .true.
    d_e = inch_const(0.125_WP, 0)
    sys_start%con(1, 2)%a_e = (PI/4.0_WP)*square(d_e)
    call sys_start%con(1, 2)%b%v%init_const(0.5_WP, 0)
    
    ! The same for every control volume.
    call x_dot%v%init_const(0.0_WP, 0)
    call y(1)%v%init_const(1.0_WP, 0)
    call p_atm%v%init_const(101350.0_WP, 0)
    call temp_atm%v%init_const(302.6_WP, 0)
    call k%v%init_const(0.0_WP, 0)
    call x_z%v%init_const(0.0_WP, 0)
    
    ! 1: chamber
    ! Appears to be constructed about 6 inches of 3/8" NPT threaded steel nipple from the photo I have.
    ! I guess I bought the pipe from Home Depot or Lowes as it doesn't appear on my old McMaster-Carr orders.
    ! Seems to be similar to <https://www.mcmaster.com/4830K158>.
    
    d_1   = inch_const(0.493_WP, 0)
    csa_1 = (PI/4.0_WP)*square(d_1)
    vol_1 = cubic_inches_const(1.1_WP, 0)
    x_1   = vol_1/csa_1
    call assert(x_1 > inch_const(5.0_WP, 0), "x_1 should be at least 5 inches")
    call assert(x_1 < inch_const(8.0_WP, 0), "x_1 should be less than 8 inches")
    p_1 = psi_const(70.0_WP, 0)
    p_1 = p_1 + p_atm
    call rm_p_1%v%init_const(0.0_WP, 0) ! immobile
    call p_fs_1%v%init_const(0.0_WP, 0)
    call p_fd_1%v%init_const(0.0_WP, 0)
    
    call sys_start%cv(1)%set(x_1, x_dot, y, p_1, temp_atm, csa_1, rm_p_1, p_fs_1, p_fd_1, p_atm, k, x_z, [DRY_AIR], &
                                isentropic_filling=.true.)
    
    ! 2: barrel
    
    d_2   = inch_const(0.527_WP, 0)
    csa_2 = (PI/4.0_WP)*square(d_2)
    vol_d = cubic_inches_const(1.1_WP, 0)
    x_2   = vol_d/csa_2
    p_2   = p_atm
    call rm_p_2%v%init_const(1.0_WP/0.98e-3_WP, 0)
    p_fs_2 = psi_const(0.5_WP, 0) ! estimate
    call p_fd_2%v%init_const(0.0_WP, 0)
    x_stop_2 = x_2 + inch_const(12.0_WP, 0)
    
    call sys_start%cv(2)%set(x_2, x_dot, y, p_2, temp_atm, csa_2, rm_p_2, p_fs_2, p_fd_2, p_atm, k, x_z, [DRY_AIR], x_stop_2)
    
    call run(sys_start, sys_end, status)
    
    call tests%integer_eq(status%rc, 0, "test_2010_08_07, status%rc")
    
    call tests%real_eq(sys_end%cv(1)%x%v%v, x_1%v%v, "test_2010_08_07, chamber end stays still")
    call tests%real_eq(sys_end%cv(1)%x_dot%v%v, 0.0_WP, "test_2010_08_07, chamber end velocity stays zero")
    
    ! characterization test
    ! This is the same number as before the fix in ed8d53d833c364f8f8b8c6788f96bb8d99dacbee.
    ! So that commit didn't break anything worse, at least.
    call tests%real_eq(sys_end%cv(2)%x_dot%v%v, 69.302132379795268_WP, "test_2010_08_07, muzzle velocity")
end subroutine test_2010_08_07

end program test_validation
