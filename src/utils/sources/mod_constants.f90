module mod_constants

	implicit none

	integer(4), parameter :: rp = 4 !(4/8)
	integer(4), parameter :: ip = 4 !(4/8)
	integer(4), parameter :: rp_vtk = 4 !(4/8)

	!
	! Dimensions
	!
	integer(ip), parameter :: ndime=3

	!
	! Element characteristics
	!
	integer(ip), parameter :: porder=4
	integer(ip), parameter :: nnode=(porder+1)**3
	integer(ip), parameter :: ngaus=nnode
	integer(ip), parameter :: npbou=(porder+1)**2

	!
	! Other constants
	!
	real(rp), parameter :: v_pi = 2.0_rp*asin(1.0_rp) ! Value of Pi

end module mod_constants
