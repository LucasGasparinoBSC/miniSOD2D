module mod_mesh

    use mod_constants
    use mod_nvtx
	use mod_maths

    implicit none

	contains

		subroutine create_mesh(nelem,npoin,xyzBase,connec,xyz)

			implicit none
			integer(ip), intent(in)  :: nelem, npoin
			real(rp),   intent(in)  :: xyzBase(nnode,ndime)
			integer(ip), intent(out) :: connec(nelem,nnode)
			real(rp),   intent(out) :: xyz(npoin,ndime)

			call nvtxStartRange("create_mesh")
			call gen_connectivity(nelem,connec)
			call gen_coordinates(nelem,npoin,connec,xyzBase,xyz)
			call nvtxEndRange

		end subroutine create_mesh

		subroutine gen_connectivity(nelem,connec)

			implicit none
			integer(ip), intent(in)  :: nelem
			integer(ip), intent(out) :: connec(nelem,nnode)
			integer(ip)              :: ielem, inode

			call nvtxStartRange("gen_connectivity")
			do ielem = 1, nelem
				do inode = 1,nnode
					connec(ielem,inode) = (ielem-1)*nnode + inode
				end do
			end do
			call nvtxEndRange

		end subroutine gen_connectivity

		subroutine gen_coordinates(nelem,npoin,connec,xyzBase,xyz)

			implicit none
			integer(ip), intent(in)    :: nelem, npoin, connec(nelem,nnode)
			real(rp),   intent(in)    :: xyzBase(nnode,ndime)
			real(rp),   intent(out)   :: xyz(npoin,ndime)
			integer(ip)                :: ielem, inode, idime, ipoin
			real(rp)                  :: xyz0

			xyz0 = 0.0_rp

			call nvtxStartRange("gen_coordinates")
			do ielem = 1, nelem
				do inode = 1,nnode
					xyz(connec(ielem,inode),1) = xyzBase(inode,1)+xyz0
					xyz(connec(ielem,inode),2) = xyzBase(inode,2)+xyz0
					xyz(connec(ielem,inode),3) = xyzBase(inode,3)+xyz0
				end do
				xyz0 = xyz0+3.0_rp
			end do
			call nvtxEndRange

		end subroutine gen_coordinates

		!----------------------------------------------------------------!
		! Tetra 04 mesh                                                  !
		!----------------------------------------------------------------!

		subroutine gen_tet_mesh(nelem,npoin,xyzBase,connec,xyz)

			implicit none
			integer(ip), intent(in)  :: nelem, npoin
			real(rp),   intent(in)  :: xyzBase(4,ndime)
			integer(ip), intent(out) :: connec(nelem,4)
			real(rp),   intent(out) :: xyz(npoin,ndime)

			call nvtxStartRange("create_tet_mesh")
			call gen_tet_connectivity(nelem,connec)
			call gen_tet_coordinates(nelem,npoin,connec,xyzBase,xyz)
			call nvtxEndRange

		end subroutine gen_tet_mesh

		subroutine gen_tet_connectivity(nelem,connec)

			implicit none
			integer(ip), intent(in)  :: nelem
			integer(ip), intent(out) :: connec(nelem,4)
			integer(ip)              :: ielem, inode

			call nvtxStartRange("gen_tet_connectivity")
			do ielem = 1, nelem
				do inode = 1,4
					connec(ielem,inode) = (ielem-1)*4 + inode
				end do
			end do
			call nvtxEndRange

		end subroutine gen_tet_connectivity

		subroutine gen_tet_coordinates(nelem,npoin,connec,xyzBase,xyz)

			implicit none
			integer(ip), intent(in)    :: nelem, npoin, connec(nelem,4)
			real(rp),   intent(in)    :: xyzBase(4,ndime)
			real(rp),   intent(out)   :: xyz(npoin,ndime)
			integer(ip)                :: ielem, inode, idime, ipoin
			real(rp)                  :: xyz0

			xyz0 = 0.0_rp

			call nvtxStartRange("gen_tet_coordinates")
			do ielem = 1, nelem
				do inode = 1,4
					xyz(connec(ielem,inode),1) = xyzBase(inode,1)+xyz0
					xyz(connec(ielem,inode),2) = xyzBase(inode,2)+xyz0
					xyz(connec(ielem,inode),3) = xyzBase(inode,3)+xyz0
				end do
				xyz0 = xyz0+3.0_rp
			end do
			call nvtxEndRange

		end subroutine gen_tet_coordinates

end module mod_mesh