!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for introducing mathematical operations that can be          !
! reutilized, such as interpolation. Included routines:               !
!                                                                     !
!   - Interpolation                                                   !
!   - GLL weights & roots                                             !
!   - Polynomial evaluation                                           !
!                                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_maths

use mod_constants

	contains

		!----------------------------------------------------------------------------------
		subroutine getGaussLobattoLegendre_roots(mporder,roots)
			integer(ip),intent(in) :: mporder
			real(rp),dimension(mporder+1),intent(out) :: roots
			real(8),dimension(mporder+1) :: roots_d,legendre_aux

			call ZELEGL(mporder,roots_d,legendre_aux)

			roots = real(roots_d)

		end subroutine getGaussLobattoLegendre_roots
		!----------------------------------------------------------------------------------

		!----------------------------------------------------------------------------------
		subroutine getGaussLobattoLegendre_weights_and_roots(mporder,weights,roots)
			integer(ip),intent(in) :: mporder
			real(rp),dimension(mporder+1),intent(out) :: weights,roots
			real(8),dimension(mporder+1) :: roots_d,weights_d,legendre_aux

			call ZELEGL(mporder,roots_d,legendre_aux)
			call WELEGL(mporder,roots_d,legendre_aux,weights_d)

			roots = real(roots_d)
			weights = real(weights_d)

		end subroutine getGaussLobattoLegendre_weights_and_roots
		!----------------------------------------------------------------------------------

		!----------------------------------------------------------------------------------
		subroutine ZELEGL(N,ET,VN)
		!********************************************************************   
		!   COMPUTES THE NODES RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA   
		!   N  = ORDER OF THE FORMULA                                           
		!   ET = VECTOR OF THE NODES, ET(I), I=0,N                              
		!   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
		!********************************************************************   
			!implicit double precision (A-H,O-Z)                               
			implicit real(8) (A-H,O-Z)
			dimension ET(0:*), VN(0:*)                                        
			integer(ip) :: N, N2, I, IT
			if (N .EQ. 0) return                                              
				N2 = (N-1)/2                                                   
				SN = DFLOAT(2*N-4*N2-3)                                        
				ET(0) = -1.D0                                                  
				ET(N) = 1.D0                                                   
				VN(0) = SN                                                     
				VN(N) = 1.D0                                                   
			if (N .EQ. 1) RETURN                                              
				ET(N2+1) = 0.D0                                                
				X = 0.D0                                                       
				call VALEPO(N,X,Y,DY,D2Y)                                         
				VN(N2+1) = Y                                                   
			if(N .EQ. 2) RETURN                                               
				C  = v_pi/DFLOAT(N)                                              
			do I=1,N2                                                       
				ETX = DCOS(C*DFLOAT(I))                                        
				do IT=1,8                                                       
				call VALEPO(N,ETX,Y,DY,D2Y)                                       
				ETX = ETX-DY/D2Y                                               
				end do                                                          
				ET(I) = -ETX                                                   
				ET(N-I) = ETX                                                  
				VN(I) = Y*SN                                                   
				VN(N-I) = Y                                                    
			end do                                                             
		end subroutine ZELEGL                                             

		subroutine WELEGL(N,ET,VN,WT)
		!********************************************************************** 
		!   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA 
		!   N  = ORDER OF THE FORMULA                                           
		!   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N                       
		!   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N   
		!   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N                            
		!********************************************************************** 
			!IMPLICIT DOUBLE PRECISION (A-H,O-Z)
			implicit real(8) (A-H,O-Z)
			dimension ET(0:*), VN(0:*), WT(0:*)
			integer(ip) :: N, N2, I
			if(N .EQ. 0) return
				N2 = (N-1)/2                                                  
				DN = DFLOAT(N)                                                
				C  = 2.D0/(DN*(DN+1.D0))                                      
			do I=0,N2                                                       
				X = ET(I)                                                     
				Y = VN(I)                                                     
				WTX = C/(Y*Y)                                                 
				WT(I) = WTX                                                   
				WT(N-I) = WTX                                                 
			end do                                                          
			if(N-1 .EQ. 2*N2) return
				X = 0.D0                                                      
				Y = VN(N2+1)                                                  
				WT(N2+1) = C/(Y*Y)                                            
		end subroutine WELEGL

		subroutine VALEPO(N,X,Y,DY,D2Y)                                   
		!*************************************************************          
		!   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N           
		!   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT               
		!   N  = DEGREE OF THE POLYNOMIAL                                       
		!   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED                    
		!   Y  = VALUE OF THE POLYNOMIAL IN X                                   
		!   DY = VALUE OF THE FIRST DERIVATIVE IN X                             
		!   D2Y= VALUE OF THE SECOND DERIVATIVE IN X                            
		!*************************************************************
			!implicit double precision (A-H,O-Z)
			implicit real(8) (A-H,O-Z)
			integer(ip) :: N, i

			Y   = 1.D0
			DY  = 0.D0
			D2Y = 0.D0
			if (N .EQ. 0) return
			Y   = X
			DY  = 1.D0
			D2Y = 0.D0
			if(N .EQ. 1) return
			YP   = 1.D0
			DYP  = 0.D0
			D2YP = 0.D0
			do I=2,N             
				C1 = DFLOAT(I)    
				C2 = 2.D0*C1-1.D0 
				C4 = C1-1.D0      
				YM = Y            
				Y  = (C2*X*Y-C4*YP)/C1                  
				YP = YM                                 
				DYM  = DY                               
				DY   = (C2*X*DY-C4*DYP+C2*YP)/C1        
				DYP  = DYM                              
				D2YM = D2Y                              
				D2Y  = (C2*X*D2Y-C4*D2YP+2.D0*C2*DYP)/C1
				D2YP = D2YM                                                    
			end do                                                                    
		end subroutine VALEPO

		pure subroutine getEquispaced_roots(mporder,xi_equi)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Computes the equispaced loci for the Lagrange    !
			! poly in the interval [-1,1], and includes the    !
			! endpoints. Roots are given in the following      !
			! order:                                           !
			!  - xi(1) = -1                                    !
			!  - xi(2) =  1                                    !
			!    xi(3:porder+1) = xi(1)+[1:N-1]*(2/N)          !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!
			implicit none
			integer(ip),intent(in) :: mporder
			real(rp),intent(out) :: xi_equi(mporder+1)
			integer(ip) :: i

			do i = 1,mporder+1
				xi_equi(i) = -1.0_rp+(2.0_rp*real(i-1,rp)/real(mporder,rp))
			end do

		end subroutine getEquispaced_roots

		pure subroutine eval_lagrangePoly(mporder,xi,xi_p,l_ip)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Evaluates the Lagrange poly of order N at a      !
			! location xi(p), where xi E [-1,1]. Returns an    !
			! array with all possible values l_ip can assume.  !
			! Remark that l_ip = 1 if p==i and 0 if            !
			! p==j~=i.                                         !
			! Expects a grid of the form:                      !
			!                                                  !
			! -1                        1 xi                   !
			!  o----v----v----v----v----o                      !
			!  0    2    3    4    5    1 i                    !
			!                                                  !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			implicit none
			integer(ip),intent(in) :: mporder
			real(rp),intent(in)   :: xi(mporder+1), xi_p
			real(rp),intent(out)  :: l_ip(mporder+1)
			integer(ip)            :: i, j, lorder(mporder+1)

			lorder(1) = 1
			lorder(2) = mporder+1
			do i = 3,mporder+1
				lorder(i) = i-1
			end do
			do i = 0,mporder ! Nodal loop
				l_ip(i+1) = 1.0_rp
				do j = 0,mporder ! Product series
				if (j .ne. (lorder(i+1)-1)) then
					l_ip(i+1) = l_ip(i+1)*((xi_p-xi(j+1))/(xi(lorder(i+1))-xi(j+1)))
				end if
				end do
			end do

		end subroutine eval_lagrangePoly

		pure subroutine eval_lagrangePolyDeriv(mporder,xi,xi_p,dl_ip)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Evaluates the derivative of the Lagrange poly    !
			! of order N at a location xi(p), where            !
			! xi E [-1,1]. Returns an array with all possible  !
			! values dl_ip can assume.                         !
			! Expects a grid of the form:                      !
			!                                                  !
			! -1                        1 xi                   !
			!  o----v----v----v----v----o                      !
			!  0    2    3    4    5    1 i                    !
			!                                                  !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			implicit none
			integer(ip),intent(in) :: mporder
			real(rp),intent(in)   :: xi(mporder+1), xi_p
			real(rp),intent(out)  :: dl_ip(mporder+1)
			integer(ip)            :: i, j, m, lorder(mporder+1)
			real(rp)              :: aux

			lorder(1) = 1
			lorder(2) = mporder+1
			do i = 3,mporder+1
				lorder(i) = i-1
			end do
			do i = 0,mporder ! Nodal loop
				dl_ip(i+1) = 0.0_rp
				do j = 0,mporder ! Product series
				aux = 1.0_rp
				if (j .ne. (lorder(i+1)-1)) then
					do m = 0,mporder
						if (m .ne. j .and. m .ne. lorder(i+1)-1) then
							aux = aux * (xi_p-xi(m+1))/(xi(lorder(i+1))-xi(m+1))
						end if
					end do
					dl_ip(i+1) = dl_ip(i+1) + (1.0_rp/(xi(lorder(i+1))-xi(j+1)))*aux
				end if
				end do
			end do

		end subroutine eval_lagrangePolyDeriv

		pure subroutine TripleTensorProduct(mporder,mnnode,xi_grid,s,t,z,atoIJK,N,dN,dlxigp_ip)

			implicit none

			integer(ip),intent(in)         :: mporder,mnnode,atoIJK(mnnode)
			real(rp),intent(in)           :: s,t,z,xi_grid(mporder+1)
			real(rp),intent(out)          :: N(mnnode),dN(ndime,mnnode)
			real(rp),optional,intent(out) :: dlxigp_ip(ndime,mporder+1)
			integer(ip)                    :: i,j,k,c
			real(rp),dimension(mporder+1) :: lxi_ip, leta_ip, lzeta_ip
			real(rp),dimension(mporder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

			call eval_lagrangePoly(mporder,xi_grid,s,lxi_ip)
			call eval_lagrangePoly(mporder,xi_grid,t,leta_ip)
			call eval_lagrangePoly(mporder,xi_grid,z,lzeta_ip)
			call eval_lagrangePolyDeriv(mporder,xi_grid,s,dlxi_ip)
			call eval_lagrangePolyDeriv(mporder,xi_grid,t,dleta_ip)
			call eval_lagrangePolyDeriv(mporder,xi_grid,z,dlzeta_ip)

			if(present(dlxigp_ip)) then
				dlxigp_ip(1,:) = dlxi_ip(:)
				dlxigp_ip(2,:) = dleta_ip(:)
				dlxigp_ip(3,:) = dlzeta_ip(:)
			end if

			c = 0
			do k = 1,mporder+1
				do i = 1,mporder+1
				do j = 1,mporder+1
					c = c+1
					N(atoIJK(c)) = lxi_ip(i)*leta_ip(j)*lzeta_ip(k)
					dN(1,atoIJK(c)) = dlxi_ip(i)*leta_ip(j)*lzeta_ip(k)
					dN(2,atoIJK(c)) = lxi_ip(i)*dleta_ip(j)*lzeta_ip(k)
					dN(3,atoIJK(c)) = lxi_ip(i)*leta_ip(j)*dlzeta_ip(k)
				end do
				end do
			end do

		end subroutine TripleTensorProduct

		!> @brief Computes the shape functions and their derivatives for QUA boundary elements
		!> @details Given a QUA_XX boundary element of order p, the routine computes
		!> a 1D Lagrangian polynomial over a given grid for x and y coordinate arguments,
		!> then sets N_a = phi_i*phi_j, where atoIJ gives the node a to IJ relationship.
		!> @param [in] xi_grid Either a Lagrange of GLL or Chebyshev 1D grid
		!> @param [in] s The xi isoparametric coordinate of the point where the shape functions are evaluated
		!> @param [in] t The eta isoparametric coordinate of the point where the shape functions are evaluated
		!> @param [in] atoIJ The node a to IJ relationship
		!> @param [out] N The shape functions
		!> @param [out] dN The derivatives of the shape functions
		pure subroutine DoubleTensorProduct(mporder,mnpbou,xi_grid,s,t,atoIJ,N,dN)

			implicit none

			integer(ip), intent(in)       :: mporder,mnpbou,atoIJ(mnpbou)
			real(rp), intent(in)         :: s, t, xi_grid(mporder+1)
			real(rp), intent(out)        :: N(mnpbou), dN(ndime-1,mnpbou)
			integer(ip)                   :: i, j, c
			real(rp), dimension(mporder+1) :: lxi_ip, leta_ip
			real(rp), dimension(mporder+1) :: dlxi_ip, dleta_ip

			call eval_lagrangePoly(mporder,xi_grid,s,lxi_ip)
			call eval_lagrangePoly(mporder,xi_grid,t,leta_ip)
			call eval_lagrangePolyDeriv(mporder,xi_grid,s,dlxi_ip)
			call eval_lagrangePolyDeriv(mporder,xi_grid,t,dleta_ip)

			c = 0
			do i = 1,mporder+1
				do j = 1,mporder+1
				c = c+1
				N(atoIJ(c)) = lxi_ip(i)*leta_ip(j)
				dN(1,atoIJ(c)) = dlxi_ip(i)*leta_ip(j)
				dN(2,atoIJ(c)) = lxi_ip(i)*dleta_ip(j)
				end do
			end do

		end subroutine DoubleTensorProduct

		pure subroutine var_interpolate(mnnode,var,Neval,var_a)
	
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			! Interpolates a variable using element shape      !
			! functions N(xi,eta,zeta). Given a coordinate     !
			! set X = (x,y,z), v(X) = N_a*v_a, where v_a are   !
			! element nodal values.
			! Expects a grid of the form:                      !
			!                                                  !
			! -1                        1 xi                   !
			!  o----v----v----v----v----o                      !
			!  0    2    3    4    5    1 i                    !
			!                                                  !
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
			implicit none
	
			integer(ip),intent(in) :: mnnode
			real(rp),intent(in)   :: var(mnnode),Neval(mnnode)
			real(rp),intent(out)  :: var_a
			integer(ip)            :: inode
	
			var_a = 0.0_rp
			do inode = 1,mnnode
				var_a = var_a+Neval(inode)*var(inode)
			end do
	
		end subroutine var_interpolate
	
end module mod_maths
