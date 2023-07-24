module mod_constants

	implicit none

	integer(4), parameter :: rp = 4 !(4/8)
	integer(4), parameter :: rp_vtk = 4 !(4/8)

	!
	! Dimensions
	!
	integer(4), parameter :: ndime=3

	!
	! Element characteristics
	!
	integer(4), parameter :: porder=4
	integer(4), parameter :: nnode=(porder+1)**3
	integer(4), parameter :: ngaus=nnode
	integer(4), parameter :: npbou=(porder+1)**2

	!
	! Other constants
	!
	real(rp), parameter :: v_pi = 2.0_rp*asin(1.0_rp) ! Value of Pi

end module mod_constants
