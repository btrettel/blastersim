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

call test_2010_08_07_20_psi(tests)
call test_2010_08_07_70_psi(tests)

call tests%end_tests()

contains

!!!!!!!!!!!!!!!!!!!!
! 2010-08-07 tests !
!!!!!!!!!!!!!!!!!!!!

! Why start here? I have validation data I collected back in the day. I wasn't able to model this set of experiments properly in the past (with my theoretical model, not computational), so I should first figure out how to model this properly before making too many assumptions about what the remainder of the model will loo like.

! /home/ben/svn/old/ballistics/data/2010-08-07 muzzle velocity test.gnumeric
! Possibly the test gun: /home/ben/photos/14/DSCN0070.JPG

! I believe the QEV was a Clippard JEV-F2M2.
! I recall doing the tests in the garage in Jefferson. It was hot outside.

! <https://web.archive.org/web/20220927075508/http://btrettel.nerfers.com/forums/topic20.html>:
! > In order to validate my simplified model of the internal ballistics of a pneumatic gun, I did some tests in 2010. I first determined what pressure was necessary for the projectile to leave the barrel. I then fired 15 shots at 30 psi, 40 psi, 50 psi, 60 psi, and 70 psi. I measured the muzzle velocity in ft/s with my chronograph.
! > These tests forced me to conclude that my model was too simple as it did not capture the dynamics of the gun as well as I had hoped. At the moment I'm planning a more complicated model suited better for spud guns, but also applicable for Nerf guns.

! TODO: Check handwritten notes for more information. I can't find anything in the scans I have.

subroutine create_2010_08_07_sys(p_psi, sys_start, x_1_)
    use convert
    use gasdata, only: DRY_AIR
    use prec, only: PI
    use checks, only: assert
    use cva, only: cv_system_type
    
    real(WP), intent(in)                           :: p_psi
    type(cv_system_type), allocatable, intent(out) :: sys_start
    type(si_length), intent(out), optional         :: x_1_
    
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
    d_e = inch_const(0.105_WP, 0)
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
    p_1 = p_atm + psi_const(p_psi, 0)
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
    
    if (present(x_1_)) x_1_ = x_1
end subroutine create_2010_08_07_sys

subroutine test_2010_08_07_20_psi(tests)
    use convert
    use cva, only: cv_system_type, run_status_type, run
    
    type(test_results_type), intent(in out) :: tests

    type(cv_system_type), allocatable :: sys_start, sys_end
    type(run_status_type)             :: status
    
    call create_2010_08_07_sys(20.0_WP, sys_start)
    
    call run(sys_start, sys_end, status)
    
    call tests%integer_eq(status%rc, 1, "test_2010_08_07_20_psi, status%rc (gun did not fire?)")
end subroutine test_2010_08_07_20_psi

subroutine test_2010_08_07_70_psi(tests)
    use convert
    use cva, only: cv_system_type, run_status_type, run
    
    type(test_results_type), intent(in out) :: tests

    type(cv_system_type), allocatable :: sys_start, sys_end
    type(run_status_type)             :: status
    type(si_length)                   :: x_1
    type(si_velocity)                 :: v_exp
    
    call create_2010_08_07_sys(70.0_WP, sys_start, x_1)
    
    call run(sys_start, sys_end, status)
    
    call tests%integer_eq(status%rc, 0, "test_2010_08_07_70_psi, status%rc")
    
    call tests%real_eq(sys_end%cv(1)%x%v%v, x_1%v%v, "test_2010_08_07_70_psi, chamber end stays still")
    call tests%real_eq(sys_end%cv(1)%x_dot%v%v, 0.0_WP, "test_2010_08_07_70_psi, chamber end velocity stays zero")
    
    call tests%real_eq(sys_end%cv(2)%x_dot%v%v, 58.291508555416044_WP, "test_2010_08_07_70_psi, muzzle velocity (characterization)")
    v_exp = fps_const(190.987_WP, 0)
    call tests%real_eq(sys_end%cv(2)%x_dot%v%v, v_exp%v%v, "test_2010_08_07_70_psi, muzzle velocity (validation)", abs_tol=1.0_WP)
end subroutine test_2010_08_07_70_psi

end program test_validation
