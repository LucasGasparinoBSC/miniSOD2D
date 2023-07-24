program main

	use mod_constants
	use mod_nvtx
	use mod_gmsh_indices
	use quadrature_rules
	use elem_qua
	use elem_hex
	use mod_mesh
	use jacobian_oper
	use elem_convec
	use elem_diffu
	use mpi

    implicit none

	!
	! Mesh vars not in mod_constants
	!
	integer(4), parameter   :: nelem = 10000
	integer(4), parameter   :: npoin = nelem * nnode
	integer(4), allocatable :: connec(:,:)
	real(rp)  , allocatable :: coord(:,:), He(:,:,:,:), gpvol(:,:,:)

	!
	! Element characteristics
	!
	integer(4), allocatable :: gmshIJK(:,:), invAtoIJK(:,:,:), AtoIJK(:), AtoI(:), AtoJ(:), AtoK(:), gmsh2sodIJK(:,:)
	real(rp)  , allocatable :: Ngp(:,:), Ngp_l(:,:), dNgp(:,:,:), dNgp_l(:,:,:), dlxigp_ip(:,:,:)
	real(rp)  , allocatable :: wgp(:), xgp(:,:)

	!
	! Loop variables
	!
	integer(4), parameter   :: nruns = 500
	integer(4)              :: idime, ielem, inode, igaus, ipoin, iorder, jnode, i, j, k

	!
	! Case variables and residuals
	!
	real(rp),   allocatable :: u(:,:), q(:,:), rho(:), pr(:), E(:), Tem(:)
	real(rp),   allocatable :: Rmass(:), Rmom(:,:), Rener(:)
	real(rp),   allocatable :: Dmass(:), Dmom(:,:), Dener(:)

	!
	! Fluid properties
	!
	real(rp)                :: Cp, Pra
	real(rp),   allocatable :: mu_fluid(:), mu_e(:,:), mu_sgs(:,:)

	!
	! Timing
	!
	real(8) :: tstart, tend, tconvec, tdiffu, tavg_convec, tavg_diffu, tmax_convec, tmax_diffu, tmin_convec, tmin_diffu
	real(8) :: tconvec_r, tdiffu_r, tavg_convec_r, tavg_diffu_r, tmax_convec_r, tmax_diffu_r, tmin_convec_r, tmin_diffu_r

	!
	! MPI vars
	!
	integer(4) :: ierr, myrank, nprocs

	!
	! Reversed copy of arrays for Guis subs
	!
	integer(4), allocatable :: connec_r(:,:)
	real(rp)  , allocatable :: u_r(:,:), q_r(:,:), Rmom_r(:,:), Dmom_r(:,:), dlxigp_ip_r(:,:,:)

	!
	! Initialize the MPI environment
	!
	call MPI_Init(ierr)

	!
	! Get number of processes and process rank
	!
	call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
	if (nprocs > 1) then
		write(*,*) 'This program only works with 1 MPI process'
		stop
	end if

	!
	! Print run config to screen
	!
	write(*,*) 'Number of elements          = ', nelem
	write(*,*) 'Element order               = ', porder
	write(*,*) 'Number of nodes per element = ', nnode
	write(*,*) 'Nodes on mesh               = ', npoin
	write(*,*) 'Number of runs              = ', nruns
	write(*,*)

	!
	! Generate isopar. element
	!

		! Initialize quad and hex info
		call init_basic_qua()
		call init_basic_hex()

		! Allocate memory for element characteristics vars.
		allocate(gmshIJK(nnode,3), invAtoIJK(porder+1,porder+1,porder+1))
		allocate(AtoIJK(nnode), AtoI(nnode), AtoJ(nnode), AtoK(nnode), gmsh2sodIJK(porder+1,2))
		allocate(Ngp(ngaus,nnode), Ngp_l(ngaus,nnode), dNgp(ndime,nnode,ngaus), dNgp_l(ndime,nnode,ngaus))
		allocate(wgp(ngaus), xgp(ngaus,ndime), dlxigp_ip(ngaus,ndime,porder+1))

		! Generate element gmshIJK
		call nvtxStartRange("Generate gmshIJK Gmsh")
		call genHighOrderHex(porder, gmshIJK)
		call nvtxEndRange
#ifdef __DEBUG__
		do inode = 1, nnode
			write(*,*) 'gmshIJK(',inode,') = ', gmshIJK(inode,:)
		end do
		write(*,*)
#endif

		! Generate gmsh2sodIJK conversion
		call nvtxStartRange("Generate gmsh2sodIJK")
		gmsh2sodIJK(1,:) = [0,1]
		gmsh2sodIJK(2,:) = [porder,2]
		do inode = 3,porder+1
			gmsh2sodIJK(inode,:) = [inode-2,inode]
		end do
		call nvtxEndRange
#ifdef __DEBUG__
		do inode = 1, porder+1
			write(*,*) 'gmsh2sodIJK(',inode,') = ', gmsh2sodIJK(inode,:)
		end do
		write(*,*)
#endif

		! Convert gmshIJK to SOD_IJK
		call nvtxStartRange("Convert gmshIJK to SOD_IJK")
		do inode = 1,nnode
			do i = 1,3
				do iorder = 1,porder+1
					if ( gmshIJK(inode,i) .eq. gmsh2sodIJK(iorder,1) ) then
						gmshIJK(inode,i) = gmsh2sodIJK(iorder,2)
						exit
					end if
				end do
			end do
		end do
		call nvtxEndRange
#ifdef __DEBUG__
		do inode = 1, nnode
			write(*,*) 'gmshIJK(',inode,') = ', gmshIJK(inode,:)
		end do
		write(*,*)
#endif

		! Generate AtoIJK and invAtoIJK
		call nvtxStartRange("Generate invAtoIJK")
		jnode = 1
		do k = 1,porder+1
			do i = 1,porder+1
				do j = 1,porder+1
					do inode = 1,nnode
						if ( gmshIJK(inode,1) .eq. i .and. gmshIJK(inode,2) .eq. j .and. gmshIJK(inode,3) .eq. k ) then
							invAtoIJK(i,j,k) = inode
							exit
						end if
					end do
#ifdef __DEBUG__
					write(*,*) i, j, k, invAtoIJK(i,j,k)
#endif
				AtoIJK(jnode) = invAtoIJK(i,j,k)
				jnode = jnode + 1
				end do
			end do
		end do
		call nvtxEndRange
#ifdef __DEBUG__
		write(*,*)
		do inode = 1, nnode
			write(*,*) 'AtoIJK(',inode,') = ', AtoIJK(inode)
		end do
		write(*,*)
#endif

		! Copy data to GPU
		call nvtxStartRange("Copyin AtoIJK and invAtoIJK")
		!$acc enter data copyin(AtoIJK,invAtoIJK)
		call nvtxEndRange

		! Get separate I, J, K arrays per node
		call nvtxStartRange("Get separate I, J, K arrays per node")
		do inode = 1,nnode
			AtoI(inode) = gmshIJK(inode,1)
			AtoJ(inode) = gmshIJK(inode,2)
			AtoK(inode) = gmshIJK(inode,3)
		end do
		call nvtxEndRange

		! Copy data to GPU
		call nvtxStartRange("Copyin AtoI, AtoJ, AtoK")
		!$acc enter data copyin(AtoI,AtoJ,AtoK)
		call nvtxEndRange

		! Generate isopar. coordinates and GLL quadrature info
		call nvtxStartRange("Generate isopar. coordinates and GLL quadrature info")
		call GaussLobattoLegendre_hex(porder,ngaus,AtoIJK,xgp,wgp)
		call nvtxEndRange
#ifdef __DEBUG__
		write(*,*)
		do igaus = 1,ngaus
			write(*,*) 'xgp(',igaus,') = ', xgp(igaus,:)
		end do
		write(*,*)
		do igaus = 1,ngaus
			write(*,*) 'wgp(',igaus,') = ', wgp(igaus)
		end do
		write(*,*)
#endif

		! Copy data to GPU
		call nvtxStartRange("Copyin xgp, wgp")
		!$acc enter data copyin(xgp,wgp)
		call nvtxEndRange

		! Generate shape functions and derivatives
		call nvtxStartRange("Generate shape functions and derivatives")
		do igaus = 1,ngaus
			call hex_highorder(porder,nnode,xgp(igaus,1),xgp(igaus,2),xgp(igaus,3), &
								AtoIJK,Ngp(igaus,:),dNgp(:,:,igaus), &
								Ngp_l(igaus,:),dNgp_l(:,:,igaus),dlxigp_ip(igaus,:,:))
		end do
		call nvtxEndRange
#ifdef __DEBUG__
		do igaus = 1,ngaus
			do inode = 1,nnode
				write(*,*) 'Ngp(',igaus,',',inode,') = ', Ngp(igaus,inode)
			end do
		end do
		write(*,*)
#endif

		! Copy data to GPU
		call nvtxStartRange("Copyin Ngp, dNgp, Ngp_l, dNgp_l, dlxigp_ip")
		!$acc enter data copyin(Ngp,dNgp,Ngp_l,dNgp_l,dlxigp_ip)
		call nvtxEndRange

	!
	! Generate mesh
	!
	allocate(connec(nelem,nnode))
	allocate(coord(npoin,ndime))
	call nvtxStartRange("Generate mesh")
	call create_mesh(nelem, npoin, xgp, connec, coord)
	call nvtxEndRange
#ifdef __DEBUG__
	do ielem = 1, nelem
		do inode = 1, nnode
			write(*,*) 'connec(',ielem,',',inode,') = ', connec(ielem,inode)
		end do
	end do
	write(*,*)
	do ipoin = 1, npoin
		write(*,*) 'coord(',ipoin,') = ', coord(ipoin,:)
	end do
	write(*,*)
#endif

	! Copy data to GPU
	call nvtxStartRange("Copyin connec, coord")
	!$acc enter data copyin(connec,coord)
	call nvtxEndRange

	!
	! Compute Jcobian info
	!
	allocate(He(ndime,ndime,ngaus,nelem))
	allocate(gpvol(1,ngaus,nelem))
		! Generate GPU memory
		call nvtxStartRange("Generate He and gpvol")
		He = 0.0_rp
		gpvol = 0.0_rp
		!$acc enter data copyin(He,gpvol)
		call nvtxEndRange
	call nvtxStartRange("Compute Jcobian info")
	call elem_jacobian(nelem,npoin,connec,coord,dNgp,wgp,gpvol,He)
	call nvtxEndRange
#ifdef __DEBUG__
	!$acc update host(He,gpvol)
	do ielem = 1,nelem
		do igaus = 1,ngaus
			do i = 1,ndime
				do j = 1,ndime
					write(*,*) 'He(',i,',',j,',',igaus,',',ielem,') = ', He(i,j,igaus,ielem)
				end do
			end do
		end do
	end do
	write(*,*)
	do ielem = 1,nelem
		do igaus = 1,ngaus
			write(*,*) 'gpvol(',igaus,',',ielem,') = ', gpvol(1,igaus,ielem)
		end do
	end do
	write(*,*)
#endif

	!
	! Generate initial conditions
	!
	allocate(u(npoin,ndime), q(npoin,ndime), rho(npoin), pr(npoin), E(npoin), Tem(npoin))
	allocate(mu_fluid(npoin), mu_e(nelem,nnode), mu_sgs(nelem,nnode))
	call nvtxStartRange("Alloc GPU initial conditions")
	!$acc enter data create(u,q,rho,pr,E,Tem,mu_fluid,mu_e,mu_sgs)
	call nvtxEndRange
	call nvtxStartRange("Generate initial conditions")
	!$acc kernels present(u,q,rho,pr,E,Tem,mu_fluid,mu_e,mu_sgs)
	u(:,:) = 1.0_rp
	q(:,:) = 1.0_rp
	rho(:) = 1.0_rp
	pr(:)  = 1.0_rp
	E(:)   = 1.0_rp
	Tem(:)   = 1.0_rp
	mu_fluid(:) = 1.0_rp
	mu_e(:,:) = 1.0_rp
	mu_sgs(:,:) = 1.0_rp
	!$acc end kernels
	call nvtxEndRange

		! change 1st node of every element
		!$acc kernels present(connec,u,q,rho,pr,E,Tem,mu_fluid,mu_e,mu_sgs)
		do ielem = 1,nelem
			u(connec(ielem,1),:) = 10.0_rp
			q(connec(ielem,1),:) = 10.0_rp
			rho(connec(ielem,1)) = 10.0_rp
			pr(connec(ielem,1))  = 10.0_rp
			E(connec(ielem,1))   = 10.0_rp
			Tem(connec(ielem,1))   = 10.0_rp
		end do
		!$acc end kernels

	!
	! Fluid properties
	!
	Cp = 1.0_rp
	Pra = 1.0_rp

	!
	! Generate the reversed arrayys
	!
	allocate(connec_r(nnode,nelem))
	allocate(u_r(ndime,npoin), q_r(ndime,npoin), Rmom_r(ndime,npoin))
	!$acc enter data create(connec_r,u_r,q_r,Rmom_r)

	! Copy connec to connec_r
	!$acc parallel loop collapse(2) present(connec,connec_r)
	do ielem = 1,nelem
		do inode = 1,nnode
			connec_r(inode,ielem) = connec(ielem,inode)
		end do
	end do
	!$acc end parallel loop

	! Copy u to u_r and q to q_r
	!$acc parallel loop collapse(2) present(u,q,u_r,q_r)
	do ipoin = 1,npoin
		do idime = 1,ndime
			u_r(idime,ipoin) = u(ipoin,idime)
			q_r(idime,ipoin) = q(ipoin,idime)
		end do
	end do
	!$acc end parallel loop

	!
	! Call the convective term multiple times
	!
	allocate(Rmass(npoin), Rmom(npoin,ndime), Rener(npoin))
	allocate(Dmass(npoin), Dmom(npoin,ndime), Dener(npoin))
	!$acc enter data create(Rmass,Rmom,Rener,Dmass,Dmom,Dener)
	open(unit=1,file='timers.dat',status='replace')
	tavg_convec = 0.0d0
	tavg_diffu = 0.0d0
	tmax_convec = 0.0d0
	tmax_diffu = 0.0d0
	tmin_convec = 1000000.0d0
	tmin_diffu = 1000000.0d0
	tavg_convec_r = 0.0d0
	tavg_diffu_r = 0.0d0
	tmax_convec_r = 0.0d0
	tmax_diffu_r = 0.0d0
	tmin_convec_r = 1000000.0d0
	tmin_diffu_r = 1000000.0d0
	call nvtxStartRange("Loop kernels")
	do i = 1,nruns
		call nvtxStartRange("Call convective term")
		tstart = MPI_Wtime()
		call full_convec_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,AtoI,AtoJ,AtoK,u,q,rho,pr,E,Rmass,Rmom,Rener)
		tend = MPI_Wtime()
		tconvec = tend-tstart
		tmax_convec = max(tmax_convec,tconvec)
		tmin_convec = min(tmin_convec,tconvec)
		tavg_convec = tavg_convec + tconvec
		call nvtxEndRange
		call nvtxStartRange("Call convec Guillaume")
		tstart = MPI_Wtime()
		call full_convec_ijk_guillaume(nelem,npoin,connec_r,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,AtoI,AtoJ,AtoK,u_r,q_r,rho,pr,E,Rmass,Rmom_r,Rener)
		tend = MPI_Wtime()
		tconvec_r = tend-tstart
		tmax_convec_r = max(tmax_convec_r,tconvec_r)
		tmin_convec_r = min(tmin_convec_r,tconvec_r)
		tavg_convec_r = tavg_convec_r + tconvec_r
		call nvtxEndRange
		call nvtxStartRange("Call diffusive term")
		tstart = MPI_Wtime()
		call full_diffusion_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,AtoI,AtoJ,AtoK,Cp,Pra,rho,u,Tem,mu_fluid,mu_e,mu_sgs,Dmass,Dmom,Dener)
		tend = MPI_Wtime()
		tdiffu = tend-tstart
		tmax_diffu = max(tmax_diffu,tdiffu)
		tmin_diffu = min(tmin_diffu,tdiffu)
		tavg_diffu = tavg_diffu + tdiffu
		write(1,*) i, tconvec, tdiffu
		call nvtxEndRange
	end do
	call nvtxEndRange
	close(1)

	!
	! Write avg times to screen
	!
	tavg_convec = tavg_convec / real(nruns,8)
	tavg_convec_r = tavg_convec_r / real(nruns,8)
	tavg_diffu = tavg_diffu / real(nruns,8)

	write(*,*)
	write(*,*) 'Timings:'
	write(*,*) '----------------------------------------'
	write(*,*) 'Avg. convective time   = ', tavg_convec
	write(*,*) 'Avg. convective_r time = ', tavg_convec_r
	write(*,*) 'Avg. diffusive time    = ', tavg_diffu
	write(*,*) 'Max. convective time   = ', tmax_convec
	write(*,*) 'Max. convective_r time = ', tmax_convec_r
	write(*,*) 'Max. diffusive time    = ', tmax_diffu
	write(*,*) 'Min. convective time   = ', tmin_convec
	write(*,*) 'Min. convective_r time = ', tmin_convec_r
	write(*,*) 'Min. diffusive time    = ', tmin_diffu
	write(*,*) '----------------------------------------'
	write(*,*) 'Variation convec.      = ', (tmax_convec-tmin_convec)/tavg_convec
	write(*,*) 'Variation convec_r.    = ', (tmax_convec_r-tmin_convec_r)/tavg_convec_r
	write(*,*) 'Variation diffu.       = ', (tmax_diffu-tmin_diffu)/tavg_diffu
	write(*,*) '----------------------------------------'

	!
	! Print minimal results set
	!
	call nvtxStartRange("Update host results")
	!$acc update host(Rmass,Rmom,Rmom_r,Rener,Dmass,Dmom,Dener)
	call nvtxEndRange
	write(*,*)
	write(*,*) 'Basic results:'
	write(*,*) '----------------------------------------'
	write(*,*) 'Max Rmass     = ', maxval(Rmass)    , 'Min Rmass     = ', minval(Rmass)
	write(*,*) 'Max Rmom(:,1) = ', maxval(Rmom(:,1)), 'Min Rmom(:,1) = ', minval(Rmom(:,1))
	write(*,*) 'Max Rmom(:,2) = ', maxval(Rmom(:,2)), 'Min Rmom(:,2) = ', minval(Rmom(:,2))
	write(*,*) 'Max Rmom(:,3) = ', maxval(Rmom(:,3)), 'Min Rmom(:,3) = ', minval(Rmom(:,3))
	write(*,*) 'Max Rmom_r(1,:) = ', maxval(Rmom_r(1,:)), 'Min Rmom_r(1,:) = ', minval(Rmom_r(1,:))
	write(*,*) 'Max Rmom_r(2,:) = ', maxval(Rmom_r(2,:)), 'Min Rmom_r(2,:) = ', minval(Rmom_r(2,:))
	write(*,*) 'Max Rmom_r(3,:) = ', maxval(Rmom_r(3,:)), 'Min Rmom_r(3,:) = ', minval(Rmom_r(3,:))
	write(*,*) 'Max Rener     = ', maxval(Rener)    , 'Min Rener     = ', minval(Rener)
	write(*,*) '----------------------------------------'
	write(*,*) 'Max Dmass     = ', maxval(Dmass)    , 'Min Dmass     = ', minval(Dmass)
	write(*,*) 'Max Dmom(:,1) = ', maxval(Dmom(:,1)), 'Min Dmom(:,1) = ', minval(Dmom(:,1))
	write(*,*) 'Max Dmom(:,2) = ', maxval(Dmom(:,2)), 'Min Dmom(:,2) = ', minval(Dmom(:,2))
	write(*,*) 'Max Dmom(:,3) = ', maxval(Dmom(:,3)), 'Min Dmom(:,3) = ', minval(Dmom(:,3))
	write(*,*) 'Max Dener     = ', maxval(Dener)    , 'Min Dener     = ', minval(Dener)

	!
	! Write results to file
	!
#ifdef __DEBUG__
	call nvtxStartRange("Write results to file")
	open(unit=10,file='results_convec.dat',status='replace')
	do ipoin = 1,npoin
		write(10,*) ipoin, Rmass(ipoin), Rmom(ipoin,1), Rmom(ipoin,2), Rmom(ipoin,3), Rener(ipoin)
	end do
	close(10)
	open(unit=11,file='results_diffu.dat',status='replace')
	do ipoin = 1,npoin
		write(11,*) ipoin, Dmass(ipoin), Dmom(ipoin,1), Dmom(ipoin,2), Dmom(ipoin,3), Dener(ipoin)
	end do
	close(11)
	call nvtxEndRange
#endif

	!
	! Finalize MPI environment
	!
	call MPI_Finalize(ierr)

end program main