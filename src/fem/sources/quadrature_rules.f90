module quadrature_rules

use mod_constants
use mod_maths

contains

		!> @brief Computes the integral of a function over a QUA_XX
		!> @details Closed rule quadrature for a function f(x) over a
		!> QUA_X element type. The integral is computed using the
		!> given grid and weights associated with a particular element order.
		!> Notice that the order of the nodes follows that of the element itself.
		!> For boundary elements only!
		!> @param[in] atoIJ Node a to IJ rellationship
		!> @param[out] xgp Quadrature points
		!> @param[out] wgp Quadrature weights
		subroutine GaussLobattoLegendre_qua(mporder,mnpbou,atoIJ,xgp,wgp)
			implicit none
			integer(ip),intent(in) :: mporder,mnpbou,atoIJ(mnpbou)
			real(rp),intent(out)  :: xgp(mnpbou,ndime-1),wgp(mnpbou)
			integer(ip)            :: inode,i,j,lorder(mporder+1)
			real(rp)              :: xi(mporder+1),w1d(mporder+1)

			call getGaussLobattoLegendre_weights_and_roots(mporder,w1d,xi)

			lorder(1) = 1
			lorder(2) = mporder+1
			do i = 3,mporder+1
				lorder(i) = i-1
			end do

			inode = 0
			do i = 1,mporder+1
				do j = 1,mporder+1
					inode = inode + 1
					xgp(atoIJ(inode),1:2) = [xi(lorder(i)), xi(lorder(j))]
					wgp(atoIJ(inode)) = w1d(lorder(i))*w1d(lorder(j))
				end do
			end do
		end subroutine GaussLobattoLegendre_qua

		!> @brief Computes the integral of a function over a HEX_XX
		!> @details Closed rule quadrature for a function f(x) over a
		!> HEX_XX element type. The integral is computed using the
		!> given grid and weights associated with a particular element order.
		!> Notice that the order of the nodes follows that of the element itself.
		!> @param[in] atoIJK Node a to IJK rellationship
		!> @param[out] xgp Quadrature points
		!> @param[out] wgp Quadrature weights
		subroutine GaussLobattoLegendre_hex(mporder,mngaus,atoIJK,xgp,wgp)
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Closed rule quadrature that uses the     !
			! Chebyshev grid+enddpoints as abcissas.   !
			! Weights are obtained by evaluating       !
			! w_j = int(l^n_i(xi_j),-1,1).             !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			implicit none
			integer(ip),intent(in) :: mporder,mngaus,atoIJK(mngaus)
			real(rp),intent(out)  :: xgp(mngaus,ndime),wgp(mngaus)
			integer(ip)            :: inode,i,j,k,lorder(mporder+1)
			real(rp)              :: xi(mporder+1),w1d(mporder+1),w0,w1,w2,w3

			call getGaussLobattoLegendre_weights_and_roots(mporder,w1d,xi)

			lorder(1) = 1
			lorder(2) = mporder+1
			do i = 3,mporder+1
				lorder(i) = i-1
			end do

			inode = 0
			do k = 1,mporder+1
				do i = 1,mporder+1
					do j = 1,mporder+1
						inode = inode + 1
						xgp(atoIJK(inode),1:3) = [xi(lorder(i)), xi(lorder(j)), xi(lorder(k))]
						wgp(atoIJK(inode)) = w1d(lorder(i))*w1d(lorder(j))*w1d(lorder(k))
					end do
				end do
			end do

		end subroutine GaussLobattoLegendre_hex

		subroutine tet_4points(xgp,wgp)
		
			implicit none
			real(rp),intent(out) :: xgp(4,ndime),wgp(4)

			xgp(1,1:3) = [0.5854101966249680_rp, 0.1381966011250110_rp, 0.1381966011250110_rp]
			xgp(2,1:3) = [0.1381966011250110_rp, 0.5854101966249680_rp, 0.1381966011250110_rp]
			xgp(3,1:3) = [0.1381966011250110_rp, 0.1381966011250110_rp, 0.5854101966249680_rp]
			xgp(4,1:3) = [0.1381966011250110_rp, 0.1381966011250110_rp, 0.1381966011250110_rp]

			wgp(1) = 0.2500000000000000_rp
			wgp(2) = 0.2500000000000000_rp
			wgp(3) = 0.2500000000000000_rp
			wgp(4) = 0.2500000000000000_rp

		end subroutine tet_4points

end module quadrature_rules
