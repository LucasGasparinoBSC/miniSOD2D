module elem_hex

use mod_constants
use mod_maths

implicit none

!integer(ip), parameter :: hex_order_edges(12,2) = transpose(reshape([1,2,1,4,1,5,2,3,2,6,3,4,3,7,4,8,5,6,5,8,6,7,7,8],(/2,12/)))
!integer(ip), parameter :: hex_order_faces(6,4)  = transpose(reshape([1,4,3,2,1,2,6,5,1,5,8,4,2,3,7,6,3,4,8,7,5,6,7,8],(/4,6/)))
integer(ip), allocatable :: hex_order_edges(:,:), hex_order_faces(:,:)

contains

	subroutine init_basic_hex()
		implicit none
		allocate(hex_order_edges(12,2))
		allocate(hex_order_faces(6,4))
		hex_order_edges = transpose(reshape([1,2,1,4,1,5,2,3,2,6,3,4,3,7,4,8,5,6,5,8,6,7,7,8],(/2,12/)))
		hex_order_faces = transpose(reshape([1,4,3,2,1,2,6,5,1,5,8,4,2,3,7,6,3,4,8,7,5,6,7,8],(/4,6/)))
		!$acc enter data copyin(hex_order_edges,hex_order_faces)
	end subroutine init_basic_hex

	subroutine hex_highorder(mporder,mnnode,xi,eta,zeta,atoIJK,N,dN,N_lagrange,dN_lagrange,dlxigp_ip)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Lagrangian HEX64 element model. Built using    !
		! equispaced nodes between [-1,1] on             !
		! (xi,eta,zeta). Ordering follows that of GMSH.  !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		implicit none
		integer(ip),intent(in) :: mporder,mnnode
		real(rp),intent(in)   :: xi, eta, zeta
		integer(ip),intent(in) :: atoIJK(mnnode)
		real(rp),intent(out)  :: N(mnnode), dN(ndime,mnnode),dlxigp_ip(ndime,mporder+1)
		real(rp),intent(out)  :: N_lagrange(mnnode), dN_lagrange(ndime,mnnode)
		real(rp)              :: xi_grid(mporder+1)

		call getEquispaced_roots(mporder,xi_grid)
		call TripleTensorProduct(mporder,mnnode,xi_grid,xi,eta,zeta,atoIJK,N,dN)
		N_lagrange(:) = N(:)
		dN_lagrange(:,:) = dN(:,:)
		call getGaussLobattoLegendre_roots(mporder,xi_grid)
		call TripleTensorProduct(mporder,mnnode,xi_grid,xi,eta,zeta,atoIJK,N,dN,dlxigp_ip)
	end subroutine hex_highorder

	subroutine hexa_edges(mnnode,ielem,nelem,npoin,connec,coord,ncorner,nedge,dist)
		implicit none
		
		integer(ip), intent(in)   :: mnnode,iElem, nelem, npoin
		integer(ip), intent(in)   :: connec(nelem,mnnode)
		real(rp),   intent(in)   :: coord(npoin,ndime)
		integer(ip), intent(out)  :: ncorner,nedge
		real(rp),   intent(out)  :: dist(12,ndime)
		integer(ip)               :: ind(mnnode)
		real(rp)                 :: xp(12,ndime)
		
		ind = connec(ielem,:)
		ncorner = 8
		nedge = 12
		
		xp(1:8,1:ndime) = coord(ind(1:8),1:ndime) ! Corner coordinates
		dist(1,:) = xp(2,:)-xp(1,:)
		dist(2,:) = xp(3,:)-xp(2,:)
		dist(3,:) = xp(4,:)-xp(3,:)
		dist(4,:) = xp(1,:)-xp(4,:)
		
		dist(5,:) = xp(6,:)-xp(5,:)
		dist(6,:) = xp(7,:)-xp(6,:)
		dist(7,:) = xp(8,:)-xp(7,:)
		dist(8,:) = xp(5,:)-xp(8,:)
		
		dist(9,:) = xp(5,:)-xp(1,:)
		dist(10,:) = xp(6,:)-xp(2,:)
		dist(11,:) = xp(7,:)-xp(3,:)
		dist(12,:) = xp(8,:)-xp(4,:)
	
	end subroutine hexa_edges

end module
