module elem_diffu

   use mod_nvtx
   use mod_numerical_params

      ! TODO: Create unit tests for all subroutines

      contains

              subroutine full_diffusion_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,Cp,Pr,rho,u,Tem,mu_fluid,mu_e,mu_sgs,Rmass,Rmom,Rener)
                      implicit none

                      integer(ip), intent(in)  :: nelem, npoin
                      integer(ip), intent(in)  :: connec(nelem,nnode)
                      real(rp),   intent(in)  :: Ngp(ngaus,nnode), dNgp(ndime,nnode,ngaus)
                      real(rp),   intent(in)  :: He(ndime,ndime,ngaus,nelem),xgp(ngaus,ndime),dlxigp_ip(ngaus,ndime,porder+1)
                      real(rp),   intent(in)  :: gpvol(1,ngaus,nelem)
                      integer(ip), intent(in)  :: invAtoIJK(porder+1,porder+1,porder+1),gmshAtoI(nnode), gmshAtoJ(nnode), gmshAtoK(nnode)
                      real(rp),   intent(in)  :: Cp,Pr,rho(npoin),u(npoin,ndime), Tem(npoin), mu_e(nelem,ngaus), mu_sgs(nelem,ngaus)
                      real(rp),   intent(in)  :: mu_fluid(npoin)
                      real(rp),   intent(out) :: Rmass(npoin)
                      real(rp),   intent(out) :: Rmom(npoin,ndime)
                      real(rp),   intent(out) :: Rener(npoin)
                      integer(ip)              :: ielem, igaus, inode, idime, jdime, isoI, isoJ, isoK,kdime,ii
                      integer(ip)              :: ipoin(nnode)
                      real(rp)                :: kappa_e, mu_fgp, mu_egp,divU, tauU(ndime), twoThirds,nu_e,tau(ndime,ndime)
                      real(rp)                :: gradU(ndime,ndime), gradT(ndime),tmp1,vol,arho
                      real(rp)                :: gradIsoRho(ndime),gradIsoT(ndime),gradIsoU(ndime,ndime)
                      real(rp)                :: gradRho(ndime),divDm(ndime),divDr,divDe
                      real(rp)                :: ul(nnode,ndime), rhol(nnode),Teml(nnode),mufluidl(nnode)
                      real(rp)                :: tauXl(nnode,ndime), tauYl(nnode,ndime), tauZl(nnode,ndime)
                      real(rp)                :: gradTl(nnode,ndime),gradRhol(nnode,ndime),tauUl(nnode,ndime)

                      call nvtxStartRange("Full diffusion")
                      twoThirds = 2.0_rp/3.0_rp
                      !$acc kernels present(Rmass,Rmom,Rener)
                      Rmass(:) = 0.0_rp
                      Rmom(:,:) = 0.0_rp
                      Rener(:) = 0.0_rp
                      !$acc end kernels

                      !$acc parallel loop gang private(ul,Teml,rhol,mufluidl,gradRhol,gradTl,tauUl,tauXl,tauYl,tauZl,ipoin) present(connec,rho,u,Tem,mu_fluid,mu_e,mu_sgs,invAtoIJK,gmshAtoI,gmshAtoJ,gmshAtoK,dlxigp_ip,He,gpvol,Rmass,Rmom,Rener)
                      do ielem = 1,nelem
                         !$acc loop vector
                         do inode = 1,nnode
                            ipoin(inode) = connec(ielem,inode)
                         end do
                         !$acc loop vector
                         do inode = 1,nnode
                            rhol(inode) = rho(ipoin(inode))
                            Teml(inode) = Tem(ipoin(inode))
                            mufluidl(inode) = mu_fluid(ipoin(inode))
                         end do
                         !$acc loop vector collapse(2)
                         do inode = 1,nnode
                            do idime = 1,ndime
                               ul(inode,idime) = u(ipoin(inode),idime)
                            end do
                         end do
                         tauXl(:,:) = 0.0_rp
                         tauYl(:,:) = 0.0_rp
                         tauZl(:,:) = 0.0_rp
                         gradTl(:,:) = 0.0_rp
                         gradRhol(:,:) = 0.0_rp
                         tauUl(:,:) = 0.0_rp

                         !$acc loop vector private(tau,gradU,gradT,tauU,gradIsoRho,gradIsoT,gradIsoU,gradRho,divU)
                         do igaus = 1,ngaus
                            nu_e = c_rho*mu_e(ielem,igaus)/rhol(igaus)
                            mu_fgp = mufluidl(igaus)+rhol(igaus)*mu_sgs(ielem,igaus)
                            mu_egp = mu_e(ielem,igaus)
                            kappa_e =mufluidl(igaus)*Cp/Pr+c_ener*mu_e(ielem,igaus)/0.4_rp + rhol(igaus)*mu_sgs(ielem,igaus)/0.9_rp

                            isoI = gmshAtoI(igaus) 
                            isoJ = gmshAtoJ(igaus) 
                            isoK = gmshAtoK(igaus) 

                            gradIsoRho(:) = 0.0_rp
                            gradIsoT(:) = 0.0_rp
                            gradIsoU(:,:) = 0.0_rp
                            !$acc loop seq
                            do ii=1,porder+1
                               gradIsoRho(1) = gradIsoRho(1) + dlxigp_ip(igaus,1,ii)*rhol(invAtoIJK(ii,isoJ,isoK))
                               gradIsoRho(2) = gradIsoRho(2) + dlxigp_ip(igaus,2,ii)*rhol(invAtoIJK(isoI,ii,isoK))
                               gradIsoRho(3) = gradIsoRho(3) + dlxigp_ip(igaus,3,ii)*rhol(invAtoIJK(isoI,isoJ,ii))

                               gradIsoT(1) = gradIsoT(1) + dlxigp_ip(igaus,1,ii)*Teml(invAtoIJK(ii,isoJ,isoK))
                               gradIsoT(2) = gradIsoT(2) + dlxigp_ip(igaus,2,ii)*Teml(invAtoIJK(isoI,ii,isoK))
                               gradIsoT(3) = gradIsoT(3) + dlxigp_ip(igaus,3,ii)*Teml(invAtoIJK(isoI,isoJ,ii))

                               !$acc loop seq
                               do idime=1,ndime
                                  gradIsoU(idime,1) = gradIsoU(idime,1) + dlxigp_ip(igaus,1,ii)*ul(invAtoIJK(ii,isoJ,isoK),idime)
                                  gradIsoU(idime,2) = gradIsoU(idime,2) + dlxigp_ip(igaus,2,ii)*ul(invAtoIJK(isoI,ii,isoK),idime)
                                  gradIsoU(idime,3) = gradIsoU(idime,3) + dlxigp_ip(igaus,3,ii)*ul(invAtoIJK(isoI,isoJ,ii),idime)
                               end do
                            end do

                            gradRho(:) = 0.0_rp
                            gradT(:) = 0.0_rp
                            gradU(:,:) = 0.0_rp
                            !$acc loop seq
                            do idime=1, ndime
                               !$acc loop seq
                               do jdime=1, ndime
                                  gradRho(idime) = gradRho(idime) + He(idime,jdime,igaus,ielem) * gradIsoRho(jdime)
                                  gradT(idime)   = gradT(idime)   + He(idime,jdime,igaus,ielem) * gradIsoT(jdime)
                                  !$acc loop seq
                                  do kdime=1,ndime
                                     gradU(idime,jdime) = gradU(idime,jdime) + He(jdime,kdime,igaus,ielem) * gradIsoU(idime,kdime)
                                  end do
                               end do
                            end do

                            divU = gradU(1,1)+gradU(2,2)+gradU(3,3)

                            tauU(:) = 0.0_rp
                            !$acc loop seq
                            do idime = 1,ndime
                               !$acc loop seq
                               do jdime = 1,ndime
                                  tauU(idime) = tauU(idime) + &
                                     (mu_fgp+mu_egp)*(gradU(idime,jdime)+ gradU(jdime,idime))*ul(igaus,jdime)
                                  tau(idime,jdime) = (mu_fgp+mu_egp)*(gradU(idime,jdime)+gradU(jdime,idime))
                               end do
                               tauU(idime) = tauU(idime)-(mu_fgp)*twoThirds*divU*ul(igaus,idime)
                               tau(idime,idime) = tau(idime,idime)-(mu_fgp)*twoThirds*divU
                            end do

                            !$acc loop seq
                            do idime = 1,ndime
                               tauXl(igaus,idime) =  tau(1,idime)
                               tauYl(igaus,idime) =  tau(2,idime)
                               tauZl(igaus,idime) =  tau(3,idime)
                               gradTl(igaus,idime) =  gradT(idime)
                               gradRhol(igaus,idime) =  gradRho(idime)
                               tauUl(igaus,idime) =  tauU(idime)
                            end do
                         end do

                         !$acc loop vector private(divDm,divDr,divDe) 
                         do igaus = 1,ngaus
                            nu_e = c_rho*mu_e(ielem,igaus)/rhol(igaus)
                            mu_fgp = mufluidl(igaus)+rhol(igaus)*mu_sgs(ielem,igaus)
                            mu_egp = mu_e(ielem,igaus)
                            kappa_e =mufluidl(igaus)*Cp/Pr+c_ener*mu_e(ielem,igaus)*Cp/Pr + Cp*rhol(igaus)*mu_sgs(ielem,igaus)/0.9_rp

                            isoI = gmshAtoI(igaus) 
                            isoJ = gmshAtoJ(igaus) 
                            isoK = gmshAtoK(igaus) 

                            divDe = 0.0_rp
                            divDr = 0.0_rp
                            divDm(:) = 0.0_rp
                            
                            !$acc loop seq
                            do ii=1,porder+1
                               !$acc loop seq
                               do idime=1,ndime
                                  divDr = divDr + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*gradRhol(invAtoIJK(ii,isoJ,isoK),idime)
                                  divDr = divDr + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*gradRhol(invAtoIJK(isoI,ii,isoK),idime)
                                  divDr = divDr + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*gradRhol(invAtoIJK(isoI,isoJ,ii),idime)

                                  divDe = divDe + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*(tauUl(invAtoIJK(ii,isoJ,isoK),idime)+kappa_e*gradTl(invAtoIJK(ii,isoJ,isoK),idime))
                                  divDe = divDe + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*(tauUl(invAtoIJK(isoI,ii,isoK),idime)+kappa_e*gradTl(invAtoIJK(isoI,ii,isoK),idime))
                                  divDe = divDe + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*(tauUl(invAtoIJK(isoI,isoJ,ii),idime)+kappa_e*gradTl(invAtoIJK(isoI,isoJ,ii),idime))
                                 
                                  divDm(1) = divDm(1) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauXl(invAtoIJK(ii,isoJ,isoK),idime)
                                  divDm(1) = divDm(1) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauXl(invAtoIJK(isoI,ii,isoK),idime)
                                  divDm(1) = divDm(1) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauXl(invAtoIJK(isoI,isoJ,ii),idime)

                                  divDm(2) = divDm(2) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauYl(invAtoIJK(ii,isoJ,isoK),idime)
                                  divDm(2) = divDm(2) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauYl(invAtoIJK(isoI,ii,isoK),idime)
                                  divDm(2) = divDm(2) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauYl(invAtoIJK(isoI,isoJ,ii),idime)

                                  divDm(3) = divDm(3) + He(idime,1,invAtoIJK(ii,isoJ,isoK),ielem)*gpvol(1,invAtoIJK(ii,isoJ,isoK),ielem)*dlxigp_ip(invAtoIJK(ii,isoJ,isoK),1,isoI)*tauZl(invAtoIJK(ii,isoJ,isoK),idime)
                                  divDm(3) = divDm(3) + He(idime,2,invAtoIJK(isoI,ii,isoK),ielem)*gpvol(1,invAtoIJK(isoI,ii,isoK),ielem)*dlxigp_ip(invAtoIJK(isoI,ii,isoK),2,isoJ)*tauZl(invAtoIJK(isoI,ii,isoK),idime)
                                  divDm(3) = divDm(3) + He(idime,3,invAtoIJK(isoI,isoJ,ii),ielem)*gpvol(1,invAtoIJK(isoI,isoJ,ii),ielem)*dlxigp_ip(invAtoIJK(isoI,isoJ,ii),3,isoK)*tauZl(invAtoIJK(isoI,isoJ,ii),idime)
                               end do
                            end do

                            !$acc atomic update
                            Rmass(ipoin(igaus)) = Rmass(ipoin(igaus))+nu_e*divDr
                            !$acc end atomic
                            !$acc atomic update
                            Rener(ipoin(igaus)) = Rener(ipoin(igaus))+divDe
                            !$acc end atomic
                            do idime = 1,ndime
                               !$acc atomic update
                               Rmom(ipoin(igaus),idime) = Rmom(ipoin(igaus),idime)+divDm(idime)
                               !$acc end atomic
                            end do
                         end do
                      end do
                      !$acc end parallel loop
                     call nvtxEndRange
              end subroutine full_diffusion_ijk

              subroutine fem_diffu(nelem,npoin,connec,Ngp,dNgp,He,gpvol,Cp,Pr,rho,u,Tem,mu_fluid,mu_e,mu_sgs,Rmass,Rmom,Rener)
                  implicit none
                  integer(ip), intent(in)  :: nelem, npoin
                  integer(ip), intent(in)  :: connec(nelem,4)
                  real(rp),   intent(in)  :: Ngp(4,4), dNgp(ndime,4,4)
                  real(rp),   intent(in)  :: He(ndime,ndime,4,nelem)
                  real(rp),   intent(in)  :: gpvol(1,4,nelem)
                  real(rp),   intent(in)  :: Cp,Pr,rho(npoin),u(npoin,ndime), Tem(npoin), mu_e(nelem,4), mu_sgs(nelem,4)
                  real(rp),   intent(in)  :: mu_fluid(npoin)
                  real(rp),   intent(out) :: Rmass(npoin)
                  real(rp),   intent(out) :: Rmom(npoin,ndime)
                  real(rp),   intent(out) :: Rener(npoin)
                  integer(ip)              :: ielem, igaus, inode, idime, jdime, jnode
                  real(rp)                :: divU, tauU(ndime), tau(ndime,ndime)
                  real(rp)                :: gradU(ndime,ndime), gradT(ndime), gradRho(ndime), gpcar(4,ndime)
                  real(rp)                :: Re_mass(4), Re_ener(4), Re_mom(4,ndime), auxRho, auxU, auxMu, auxNu, auxKappa

                  call nvtxStartRange("FEM diffu")
                  !$acc kernels
                  Rmass(:) = 0.0_rp
                  Rener(:) = 0.0_rp
                  Rmom(:,:) = 0.0_rp
                  !$acc end kernels

                  !$acc parallel loop gang private(Re_mass, Re_ener, Re_mom, gpcar) present(connec,rho,u,Tem,mu_fluid,mu_e,mu_sgs,He,gpvol,Rmass,Rmom,Rener)
                  do ielem = 1,nelem
                     Re_mass(:) = 0.0_rp
                     Re_ener(:) = 0.0_rp
                     Re_mom(:,:) = 0.0_rp
                     !$acc loop seq private(gradRho, gradU, gradT, tau, tauU, divU)
                     do igaus = 1,4
                        !
                        ! Fluid props
                        !
                        auxRho = 0.0_rp
                        auxNu = 0.0_rp
                        auxMu = 0.0_rp
                        auxKappa = 0.0_rp
                        do inode = 1,4
                           auxRho = auxRho + Ngp(inode,igaus) * rho(connec(ielem,inode))
                        end do
                        do inode = 1,4
                           auxNu = auxNu + Ngp(inode,igaus)*mu_e(ielem,inode)*(c_rho / auxRho)
                           auxMu = auxMu + Ngp(inode,igaus) * ( mu_fluid(connec(ielem,inode)) + auxRho*mu_sgs(ielem,inode) + mu_e(ielem,inode) )
                           auxKappa = auxKappa + Ngp(inode,igaus) * ( mu_fluid(connec(ielem,inode))*(Cp/Pr) + c_ener*mu_e(ielem,inode)/0.4_rp + auxRho*mu_sgs(ielem,inode)/0.9_rp )
                        end do
                        !
                        ! GPCAR
                        !
                        do idime = 1,ndime
                           do inode = 1,4
                              gpcar(inode,idime) = DOT_PRODUCT(He(idime,:,igaus,ielem),dNgp(:,inode,igaus))
                           end do
                        end do
                        !
                        ! Gradients and divU
                        !
                        divU = 0.0_rp
                        gradRho(:) = 0.0_rp
                        gradT(:) = 0.0_rp
                        gradU(:,:) = 0.0_rp
                        do idime = 1,ndime
                           do inode = 1,4
                              gradRho(idime) = gradRho(idime) + gpcar(inode,idime)*rho(connec(ielem,inode))
                              gradT(idime)   = gradT(idime)   + gpcar(inode,idime)*Tem(connec(ielem,inode))
                              do jdime = 1,ndime
                                 gradU(idime,jdime) = gradU(idime,jdime) + gpcar(inode,idime)*u(connec(ielem,inode),jdime)
                              end do
                           end do
                        end do
                        divU = gradU(1,1) + gradU(2,2) + gradU(3,3)
                        !
                        ! Tau and TauU (no mu applied)
                        !
                        tau(:,:) = 0.0_rp
                        tauU(:) = 0.0_rp
                        do idime = 1,ndime
                           tau(idime,idime) = -(2.0_rp/3.0_rp)*divU
                           do jdime = 1,ndime
                              tau(idime,jdime) = tau(idime,jdime) + gradU(idime,jdime) + gradU(jdime,idime)
                           end do
                        end do
                        do jdime = 1,ndime
                           auxU = 0.0_rp
                           do inode = 1,4
                              auxU = auxU + Ngp(inode,igaus)*u(connec(ielem,inode),jdime)
                           end do
                           do idime = 1,ndime
                              tauU(idime) = tauU(idime) + tau(idime,jdime)*auxU
                           end do
                        end do
                        !
                        ! Local residuals
                        !
                        do idime = 1,ndime
                           do inode = 1,4
                              Re_mass(inode) = Re_mass(inode) + gpvol(1,igaus,ielem)*auxNu*gpcar(inode,idime)*gradRho(idime)
                              Re_ener(inode) = Re_ener(inode) + gpvol(1,igaus,ielem)*gpcar(inode,idime)*(auxMu*tauU(idime)+auxKappa*gradT(idime))
                              do jdime = 1,ndime
                                 Re_mom(inode,idime) = Re_mom(inode,idime) + gpvol(1,igaus,ielem)*gpcar(inode,jdime)*auxMu*tau(idime,jdime)
                              end do
                           end do
                        end do
                     end do
                     !
                     ! Assembly
                     !
                     do inode = 1,4
                        Rmass(connec(ielem,inode)) = Rmass(connec(ielem,inode)) + Re_mass(inode)
                        Rener(connec(ielem,inode)) = Rener(connec(ielem,inode)) + Re_ener(inode)
                        do idime = 1,ndime
                           Rmom(connec(ielem,inode),idime) = Rmom(connec(ielem,inode),idime) + Re_mom(inode,idime)
                        end do
                     end do
                  end do
                  !$acc end parallel loop
                  call nvtxEndRange

              end subroutine fem_diffu

end module elem_diffu
