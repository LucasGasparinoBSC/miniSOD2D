	module elem_convec

	use mod_nvtx
	use mod_constants

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Computes convective term for Euler/NS equation system, as well             !
		! as for any generic scalar transport that might occur. Based                !
		! on Ljunkvist matrix-free implementation (assembles only rhs vector).       !
		! This module can be passed to CUDA in order to do fine-grained parallelism. !
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		contains

				subroutine full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,u,q,rho,pr,E,Rmass,Rmom,Rener)

					implicit none

					integer(4), intent(in)  :: nelem, npoin
					integer(4), intent(in)  :: connec(nelem,nnode)
					real(rp),    intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
					real(rp),    intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
					real(rp),    intent(in)  :: gpvol(1,ngaus,nelem)
					integer(4), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
					real(rp),    intent(in)  :: q(npoin,ndime), u(npoin,ndime), rho(npoin),pr(npoin), E(npoin)
					real(rp),    intent(out) :: Rmass(npoin)
					real(rp),    intent(out) :: Rmom(npoin,ndime)
					real(rp),    intent(out) :: Rener(npoin)
					integer(4)              :: ielem, igaus, idime, jdime, inode, isoI, isoJ, isoK,kdime,ii
					real(rp)                 :: Re_mom(nnode,ndime)
					real(rp)                 :: Re_mass(nnode), Re_ener(nnode)
					real(rp)                 :: gradIsoRho(ndime),gradIsoP(ndime), gradIsoE(ndime),gradIsoU(ndime,ndime), gradIsoF(ndime,ndime,ndime), gradIsoQ(ndime,ndime), gradIsoFe(ndime,ndime)
					real(rp)                 :: gradRho(ndime),gradP(ndime),gradE(ndime),gradU(ndime,ndime),divF(ndime),divU,divFe,divQ
					real(rp)                 :: ul(nnode,ndime), ql(nnode,ndime), rhol(nnode), prl(nnode),El(nnode),fel(nnode,ndime),fl(nnode,ndime,ndime)
					real(rp), dimension(porder+1) :: dlxi_ip, dleta_ip, dlzeta_ip

					call nvtxStartRange("Full convection")
					!$acc kernels
					Rmom(:,:) = 0.0_rp
					Rmass(:) = 0.0_rp
					Rener(:) = 0.0_rp
					!$acc end kernels

					!$acc parallel loop gang private(Re_ener,Re_mass,Re_mom,ul,ql,rhol,prl,El,fl,fel) !!vector_length(vecLength)
					do ielem = 1,nelem
						!$acc loop vector collapse(2)
						do idime = 1,ndime
						do inode = 1,nnode
							ul(inode,idime) = u(connec(ielem,inode),idime)
							ql(inode,idime) = q(connec(ielem,inode),idime)
							fel(inode,idime) = (E(connec(ielem,inode))+pr(connec(ielem,inode)))*u(connec(ielem,inode),idime)
						end do
						end do
						!$acc loop vector collapse(3)
						do idime = 1,ndime
						do jdime = 1,ndime
							do inode = 1,nnode
								fl(inode,idime,jdime)  = q(connec(ielem,inode),idime)*u(connec(ielem,inode),jdime)
							end do
						end do
						end do
						!$acc loop vector
						do inode = 1,nnode
						rhol(inode) = rho(connec(ielem,inode))
						El(inode) = E(connec(ielem,inode))
						prl(inode) = pr(connec(ielem,inode))
						end do
						!$acc loop vector private(dlxi_ip,dleta_ip,dlzeta_ip, gradIsoRho,gradIsoP, gradIsoE,gradIsoU, gradIsoF, gradIsoQ, gradIsoFe,gradRho,gradP,gradE,gradU,divF,divU,divQ,divFe)
						do igaus = 1,ngaus
						!$acc loop seq
						do ii=1,porder+1
							dlxi_ip(ii) = dlxigp_ip(igaus,1,ii)
							dleta_ip(ii) = dlxigp_ip(igaus,2,ii)
							dlzeta_ip(ii) = dlxigp_ip(igaus,3,ii)
						end do
						isoI = gmshAtoI(igaus) 
						isoJ = gmshAtoJ(igaus) 
						isoK = gmshAtoK(igaus) 

						gradIsoRho(:) = 0.0_rp
						gradIsoP(:) = 0.0_rp
						gradIsoE(:) = 0.0_rp
						gradIsoU(:,:) = 0.0_rp
						gradIsoF(:,:,:) = 0._rp
						gradIsoQ(:,:) = 0._rp
						gradIsoFe(:,:) = 0._rp
						!$acc loop seq
						do ii=1,porder+1
							gradIsoRho(1) = gradIsoRho(1) + dlxi_ip(ii)*rhol(invAtoIJK(ii,isoJ,isoK))
							gradIsoRho(2) = gradIsoRho(2) + dleta_ip(ii)*rhol(invAtoIJK(isoI,ii,isoK))
							gradIsoRho(3) = gradIsoRho(3) + dlzeta_ip(ii)*rhol(invAtoIJK(isoI,isoJ,ii))

							gradIsoP(1) = gradIsoP(1) + dlxi_ip(ii)*prl(invAtoIJK(ii,isoJ,isoK))
							gradIsoP(2) = gradIsoP(2) + dleta_ip(ii)*prl(invAtoIJK(isoI,ii,isoK))
							gradIsoP(3) = gradIsoP(3) + dlzeta_ip(ii)*prl(invAtoIJK(isoI,isoJ,ii))

							gradIsoE(1) = gradIsoE(1) + dlxi_ip(ii)*(prl(invAtoIJK(ii,isoJ,isoK))  + El(invAtoIJK(ii,isoJ,isoK)))
							gradIsoE(2) = gradIsoE(2) + dleta_ip(ii)*(prl(invAtoIJK(isoI,ii,isoK)) + El(invAtoIJK(isoI,ii,isoK)))
							gradIsoE(3) = gradIsoE(3) + dlzeta_ip(ii)*(prl(invAtoIJK(isoI,isoJ,ii))+ El(invAtoIJK(isoI,isoJ,ii)))
							
							!$acc loop seq
							do idime=1,ndime
								gradIsoU(idime,1) = gradIsoU(idime,1) + dlxi_ip(ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
								gradIsoU(idime,2) = gradIsoU(idime,2) + dleta_ip(ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
								gradIsoU(idime,3) = gradIsoU(idime,3) + dlzeta_ip(ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)

								gradIsoQ(idime,1) = gradIsoQ(idime,1) + dlxi_ip(ii)*ql(invAtoIJK(ii,isoJ,isoK),idime)
								gradIsoQ(idime,2) = gradIsoQ(idime,2) + dleta_ip(ii)*ql(invAtoIJK(isoI,ii,isoK),idime)
								gradIsoQ(idime,3) = gradIsoQ(idime,3) + dlzeta_ip(ii)*ql(invAtoIJK(isoI,isoJ,ii),idime)

								gradIsoFe(idime,1) = gradIsoFe(idime,1) + dlxi_ip(ii)*fel(invAtoIJK(ii,isoJ,isoK),idime)
								gradIsoFe(idime,2) = gradIsoFe(idime,2) + dleta_ip(ii)*fel(invAtoIJK(isoI,ii,isoK),idime)
								gradIsoFe(idime,3) = gradIsoFe(idime,3) + dlzeta_ip(ii)*fel(invAtoIJK(isoI,isoJ,ii),idime)
								
								!$acc loop seq
								do jdime=1, ndime
									gradIsoF(idime,jdime,1) = gradIsoF(idime,jdime,1) + dlxi_ip(ii)*fl(invAtoIJK(ii,isoJ,isoK),idime,jdime)
									gradIsoF(idime,jdime,2) = gradIsoF(idime,jdime,2) + dleta_ip(ii)*fl(invAtoIJK(isoI,ii,isoK),idime,jdime)
									gradIsoF(idime,jdime,3) = gradIsoF(idime,jdime,3) + dlzeta_ip(ii)*fl(invAtoIJK(isoI,isoJ,ii),idime,jdime)
								end do
							end do
						end do

						gradRho(:) = 0.0_rp
						gradP(:) = 0.0_rp
						gradE(:) = 0.0_rp
						gradU(:,:) = 0.0_rp
						divF(:) = 0.0_rp
						divQ = 0.0_rp
						divFe = 0.0_rp
						!$acc loop seq
						do idime=1, ndime
							!$acc loop seq
							do jdime=1, ndime
								gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
								gradP(idime)   = gradP(idime) + He(idime,jdime,igaus,ielem) * gradIsoP(jdime)
								gradE(idime)   = gradE(idime) + He(idime,jdime,igaus,ielem) * gradIsoE(jdime)
								divQ = divQ + He(idime,jdime,igaus,ielem) * gradIsoQ(idime,jdime)
								divFe = divFe + He(idime,jdime,igaus,ielem) * gradIsoFe(idime,jdime)
								!$acc loop seq
								do kdime=1,ndime
									gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
									divF(idime) = divF(idime) + He(jdime,kdime,igaus,ielem)*gradIsoF(idime,jdime,kdime)
								end do
							end do
						end do
						divU  = gradU(1,1)  + gradU(2,2)  + gradU(3,3) 
						Re_mass(igaus) = 0.5_rp*(divQ+rhol(igaus)*divU) 
						Re_ener(igaus) = 0.5_rp*(divFe+(El(igaus)+prl(igaus))*divU)
						!$acc loop seq
						do idime=1, ndime
							Re_mom(igaus,idime) = 0.5_rp*(divF(idime)+ql(igaus,idime)*divU)+ gradP(idime)
							Re_mass(igaus) = Re_mass(igaus) + 0.5_rp*(ul(igaus,idime)*gradRho(idime))
							Re_ener(igaus) = Re_ener(igaus) + 0.5_rp*(ul(igaus,idime)*gradE(idime))
							!$acc loop seq
							do jdime=1, ndime
								Re_mom(igaus,idime) = Re_mom(igaus,idime) + 0.5_rp*(ul(igaus,idime)*ul(igaus,jdime)*gradRho(jdime) &
																		+ ql(igaus,jdime)*gradU(idime,jdime))
							end do
							Re_mom(igaus,idime) = gpvol(1,igaus,ielem)*Re_mom(igaus,idime)
						end do
						Re_mass(igaus) = gpvol(1,igaus,ielem)*Re_mass(igaus)
						Re_ener(igaus) = gpvol(1,igaus,ielem)*Re_ener(igaus)
						end do
						!
						! Final assembly
						!
						!$acc loop vector collapse(2)
						do idime = 1,ndime
						do inode = 1,nnode
							!$acc atomic update
							Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime)+Re_mom(inode,idime)
							!$acc end atomic
						end do
						end do
						!$acc loop vector
						do inode = 1,nnode
						!$acc atomic update
						Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode))+Re_mass(inode)
						!$acc end atomic
						!$acc atomic update
						Rener(connec(ielem,inode)) = Rener(connec(ielem,inode))+Re_ener(inode)
						!$acc end atomic
						end do
					end do
					!$acc end parallel loop
					call nvtxEndRange

				end subroutine full_convec_ijk

	end module elem_convec
