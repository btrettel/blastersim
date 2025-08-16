! Module to set precisions, lengths, and mathematical constants.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [flt](https://github.com/btrettel/flt)
! License: [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

module prec

implicit none
private

! Integer precision
! -----------------

integer, public, parameter :: I5 = selected_int_kind(5)
integer, public, parameter :: I9 = selected_int_kind(9) ! This is the most ELF90 allows.
integer, public, parameter :: I10 = selected_int_kind(10)

! Real precision
! --------------

! <https://fortranwiki.org/fortran/show/Real+precision>
! `WP` stands for *working precision* in case I want to change the precision later.
integer, public, parameter :: SP = selected_real_kind(6, 37)    ! single
integer, public, parameter :: DP = selected_real_kind(15, 307)  ! double
!integer, public, parameter :: QP = selected_real_kind(33, 4931) ! quad
integer, public, parameter :: WP = DP

! Number of digits of acceptable precision loss in catastrophic cancellation.
integer, public, parameter :: ACCEPTABLE_LOG10_SPACING_JUMP = floor(real(precision(1.0_WP), WP)/4.0_WP)

! String parameters
! -----------------

integer, public, parameter :: CL = 1024 ! `CL` stands for character length

! Mathematical constants
! ----------------------

! PI: 3.1415926535897932384626433...
! <https://stackoverflow.com/a/49416279/1124489>
real(WP), public, parameter :: PI = 4.0_WP * atan(1.0_WP)

end module prec
