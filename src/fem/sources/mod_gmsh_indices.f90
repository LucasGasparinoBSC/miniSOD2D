module mod_gmsh_indices

use mod_constants
use elem_qua ! Using only the ordering table for edges
use elem_hex ! Using only the ordering tables for edges and faces

	implicit none

	contains

		recursive subroutine genHighOrderHex(p,indexTable)
			implicit none
			integer(ip), intent(in)  :: p
			integer(ip), intent(out) :: indexTable((p+1)**3,3) ! Set as I, J, K
			integer(ip)              :: tableFace((p-1)**2,2), tableVolume((p-1)**3,3)
			integer(ip)              :: inode, iedge, iface, i0, i1, i3, u(3), v(3), i

			! Initialize corner node to 0th position, or generate element of order 0
			indexTable(1,:) = [0,0,0]

			! Generate element of order 1 (corner nodes)
			if (p .gt. 0) then
				!                  i, j, k
				indexTable(2,:) = [p, 0, 0]
				indexTable(3,:) = [p, p, 0]
				indexTable(4,:) = [0, p, 0]
				indexTable(5,:) = [0, 0, p]
				indexTable(6,:) = [p, 0, p]
				indexTable(7,:) = [p, p, p]
				indexTable(8,:) = [0, p, p]
				if (p .gt. 1) then
				! Generate high-order edges
				inode = 9
				do iedge = 1,12
					i0 = hex_order_edges(iedge,1)
					i1 = hex_order_edges(iedge,2)
					u(:) = (indexTable(i1,:) - indexTable(i0,:))/p
					do i = 1,p-1
						indexTable(inode,:) = indexTable(i0,:) + i*u(:)
						inode = inode + 1
					end do
				end do
				! Generate a generic high-order face with p` = p-2
				call genHighOrderQuad(p-2,tableFace)
				tableFace = tableFace + 1
				! Generate faces interior nodes
				do iface = 1,6
					i0 = hex_order_faces(iface,1)
					i1 = hex_order_faces(iface,2)
					i3 = hex_order_faces(iface,4)
					u(:) = (indexTable(i1,:) - indexTable(i0,:))/p
					v(:) = (indexTable(i3,:) - indexTable(i0,:))/p
					do i = 1,((p-1)**2)
						indexTable(inode,:) = indexTable(i0,:) + u(:)*tableFace(i,1) + v(:)*tableFace(i,2)
						inode = inode + 1
					end do
				end do
				! Generate volume nodes
				call genHighOrderHex(p-2,tableVolume)
				tableVolume = tableVolume + 1
				call joinTables([(p-1)**3,3],tableVolume,inode,[(p+1)**3,3],indexTable)
				end if
			end if
		end subroutine genHighOrderHex

		recursive subroutine genHighOrderQuad(p,indexTable)
			implicit none
			integer(ip), intent(in)  :: p
			integer(ip), intent(out) :: indexTable((p+1)**2,2) ! Set as I, J
			integer(ip)              :: tableFace((p-1)**2,2)
			integer(ip)              :: inode, iedge, iface, i0, i1, u(2), i

			indexTable(1,:) = [0,0]
			if (p .gt. 0) then
				indexTable(2,:) = [p,0]
				indexTable(3,:) = [p,p]
				indexTable(4,:) = [0,p]
				if (p .gt. 1) then
				inode = 5
				do iedge = 1,4
					i0 = quad_order_edges(iedge,1)
					i1 = quad_order_edges(iedge,2)
					u(:) = (indexTable(i1,:) - indexTable(i0,:))/p
					do i = 1,p-1
						indexTable(inode,:) = indexTable(i0,:) + i*u(:)
						inode = inode + 1
					end do
				end do
				call genHighOrderQuad(p-2,tableFace)
				tableFace = tableFace + 1
				call joinTables([(p-1)**2,2],tableFace,inode,[(p+1)**2,2],indexTable)
				end if
			end if
		end subroutine genHighOrderQuad

		subroutine joinTables(size1,table1,indexDesti,size2,table2)
			implicit none
			integer(ip), intent(in)    :: indexDesti, size1(2), size2(2)
			integer(ip), intent(in)    :: table1(size1(1),size1(2))
			integer(ip), intent(inout) :: table2(size2(1),size2(2))
			integer(ip)                :: i, j

			j = indexDesti
			do i = 1,size1(1)
				table2(j,:) = table1(i,:)
				j = j + 1
			end do
		end subroutine joinTables

end module mod_gmsh_indices