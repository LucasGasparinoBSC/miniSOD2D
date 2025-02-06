program main

	use mod_constants
	use mod_nvtx
	use mod_gmsh_indices
	use quadrature_rules
	use elem_qua
	use elem_hex
	use elem_tet
	use mod_mesh
	use jacobian_oper
	use elem_convec
	use elem_diffu
	use mpi

    implicit none

	!
	! Mesh vars not in mod_constants
	!
	integer(4), parameter   :: nelem = 4*24000
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
	integer(4), parameter   :: nruns = 100
	integer(4)              :: ielem, inode, igaus, ipoin, iorder, jnode, i, j, k, idime

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
	real(8) :: tconvec_split, tdiffu_split, tavg_convec_split, tavg_diffu_split, tmax_convec_split, tmax_diffu_split, tmin_convec_split, tmin_diffu_split
	real(8) :: tconv_tet, tdiff_tet, tmax_convec_tet, tmax_diffu_tet, tmin_convec_tet, tmin_diffu_tet, tavg_convec_tet, tavg_diffu_tet

	!
	! MPI vars
	!
	integer(4) :: ierr, myrank, nprocs

	!
	! Tetra variables
	!
	integer(4), parameter   :: nnode_t = 4
	integer(4), parameter   :: ngaus_t = 4
	integer(4)              :: npoin_t, nelem_t
	integer(4), allocatable :: connec_t(:,:)
	real(rp)  , allocatable :: coord_t(:,:), He_t(:,:,:,:), gpvol_t(:,:,:), xyzTET(:,:)
	real(rp)  , allocatable :: Ngp_t(:,:), dNgp_t(:,:,:)
	real(rp)  , allocatable :: wgp_t(:), xgp_t(:,:)
	real(rp),   allocatable :: Rmass_s(:), Rmom_s(:,:), Rener_s(:)
	real(rp),   allocatable :: Rmass_t(:), Rmom_t(:,:), Rener_t(:)
	real(rp),   allocatable :: Dmass_t(:), Dmom_t(:,:), Dener_t(:)
	real(rp),   allocatable :: mu_e_t(:,:), mu_sgs_t(:,:)

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
			rho(connec(ielem,1)) = 10.0_rp
			pr(connec(ielem,1))  = 10.0_rp
			E(connec(ielem,1))   = 10.0_rp
			Tem(connec(ielem,1))   = 10.0_rp
		end do
		!$acc end kernels

		!$acc parallel loop collapse(2)
		do ipoin = 1,npoin
			do idime = 1,ndime
				q(ipoin,idime) = rho(ipoin)*u(ipoin,idime)
			end do
		end do
		!$acc end parallel loop

	!
	! Generate TET04 data
	!
		npoin_t = npoin
		nelem_t = npoin_t/nnode_t
		if (MODULO(npoin_t,nnode_t) .ne. 0) then
			write(*,*) 'Error: npoin_t is not divisible by nnode_t'
			stop
		end if

		allocate(Ngp_t(ngaus_t,nnode_t), dNgp_t(ndime,nnode_t,ngaus_t))
		allocate(wgp_t(ngaus_t), xgp_t(ngaus_t,ndime))

		call tet_4points(xgp_t,wgp_t)
		!$acc enter data copyin(xgp_t,wgp_t)

		do igaus = 1,ngaus_t
			call tet_04(xgp_t(igaus,1),xgp_t(igaus,2),xgp_t(igaus,3), &
						Ngp_t(igaus,:),dNgp_t(:,:,igaus))
		end do
		!$acc enter data copyin(Ngp_t,dNgp_t)

		allocate(connec_t(nelem_t,nnode_t))
		allocate(coord_t(npoin_t,ndime))
		allocate(xyzTET(4,3))

		xyzTET(1,1:3) = [0.0_rp, 0.0_rp, 0.0_rp]
		xyzTET(2,1:3) = [1.0_rp, 0.0_rp, 0.0_rp]
		xyzTET(3,1:3) = [0.0_rp, 1.0_rp, 0.0_rp]
		xyzTET(4,1:3) = [0.0_rp, 0.0_rp, 1.0_rp]
		call gen_tet_mesh(nelem_t,npoin_t,xyzTET,connec_t,coord_t)
		!$acc enter data copyin(connec_t,coord_t)

		allocate(He_t(ndime,ndime,ngaus_t,nelem_t))
		allocate(gpvol_t(1,ngaus_t,nelem_t))
		!$acc enter data create(He_t,gpvol_t)
		call tet_jacobian(nelem_t,npoin_t,connec_t,coord_t,dNgp_t,wgp_t,gpvol_t,He_t)

		allocate(mu_e_t(nelem_t,nnode_t), mu_sgs_t(nelem_t,nnode_t))
		!$acc enter data create(mu_e_t,mu_sgs_t)
		!$acc kernels present(mu_e_t,mu_sgs_t)
		mu_e_t(:,:) = 1.0_rp
		mu_sgs_t(:,:) = 1.0_rp
		!$acc end kernels

	!
	! Fluid properties
	!
	Cp = 1.0_rp
	Pra = 1.0_rp

	!
	! Call the convective term multiple times
	!
	allocate(Rmass(npoin), Rmom(npoin,ndime), Rener(npoin))
	allocate(Dmass(npoin), Dmom(npoin,ndime), Dener(npoin))
	allocate(Rmass_s(npoin_t), Rmom_s(npoin_t,ndime), Rener_s(npoin_t))
	allocate(Rmass_t(npoin_t), Rmom_t(npoin_t,ndime), Rener_t(npoin_t))
	allocate(Dmass_t(npoin_t), Dmom_t(npoin_t,ndime), Dener_t(npoin_t))
	!$acc enter data create(Rmass,Rmom,Rener,Dmass,Dmom,Dener)
	!$acc enter data create(Rmass_s,Rmom_s,Rener_s)
	!$acc enter data create(Rmass_t,Rmom_t,Rener_t,Dmass_t,Dmom_t,Dener_t)

	open(unit=1,file='timers.dat',status='replace')
	tavg_convec = 0.0d0
	tavg_convec_split = 0.0d0
	tavg_convec_tet = 0.0d0
	tavg_diffu = 0.0d0
	tavg_diffu_tet = 0.0d0
	tmax_convec = 0.0d0
	tmax_convec_split = 0.0d0
	tmax_convec_tet = 0.0d0
	tmax_diffu = 0.0d0
	tmax_diffu_tet = 0.0d0
	tmin_convec = 1000000.0d0
	tmin_convec_split = 1000000.0d0
	tmin_convec_tet = 1000000.0d0
	tmin_diffu = 1000000.0d0
	tmin_diffu_tet = 1000000.0d0
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

		call nvtxStartRange("Call splited convective term")
		tstart = MPI_Wtime()
		call full_convec_ijk_split(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,AtoI,AtoJ,AtoK,u,q,rho,pr,E,Rmass_s,Rmom_s,Rener_s)
		tend = MPI_Wtime()
		tconvec_split = tend-tstart
		tmax_convec_split = max(tmax_convec_split,tconvec_split)
		tmin_convec_split = min(tmin_convec_split,tconvec_split)
		tavg_convec_split = tavg_convec_split + tconvec_split
		call nvtxEndRange

		call nvtxStartRange("Call convective term TET")
		tstart = MPI_Wtime()
		call fem_convec(nelem_t,npoin_t,connec_t,Ngp_t,dNgp_t,He_t,gpvol_t,u,q,rho,pr,E,Rmass_t,Rmom_t,Rener_t)
		tend = MPI_Wtime()
		tconv_tet = tend-tstart
		tmax_convec_tet = max(tmax_convec_tet,tconv_tet)
		tmin_convec_tet = min(tmin_convec_tet,tconv_tet)
		tavg_convec_tet = tavg_convec_tet + tconv_tet
		call nvtxEndRange

		call nvtxStartRange("Call diffusive term")
		tstart = MPI_Wtime()
		call full_diffusion_ijk(nelem,npoin,connec,Ngp,dNgp,He,gpvol,dlxigp_ip,xgp,invAtoIJK,AtoI,AtoJ,AtoK,Cp,Pra,rho,u,Tem,mu_fluid,mu_e,mu_sgs,Dmass,Dmom,Dener)
		tend = MPI_Wtime()
		tdiffu = tend-tstart
		tmax_diffu = max(tmax_diffu,tdiffu)
		tmin_diffu = min(tmin_diffu,tdiffu)
		tavg_diffu = tavg_diffu + tdiffu
		call nvtxEndRange

		call nvtxStartRange("Call diffusive term TET")
		tstart = MPI_Wtime()
		call fem_diffu(nelem_t,npoin_t,connec_t,Ngp_t,dNgp_t,He_t,gpvol_t,Cp,Pra,rho,u,Tem,mu_fluid,mu_e_t,mu_sgs_t,Dmass_t,Dmom_t,Dener_t)
		tend = MPI_Wtime()
		tdiff_tet = tend-tstart
		tmax_diffu_tet = max(tmax_diffu_tet,tdiff_tet)
		tmin_diffu_tet = min(tmin_diffu_tet,tdiff_tet)
		tavg_diffu_tet = tavg_diffu_tet + tdiff_tet
		call nvtxEndRange

		write(1,10) i, tconvec, tconv_tet, tdiffu, tdiff_tet
	end do
	call nvtxEndRange
	close(1)
	10 format(i3,1x,4(f16.8,1x))

	!
	! Write avg times to screen
	!
	tavg_convec = tavg_convec / real(nruns,8)
	tavg_convec_split = tavg_convec_split / real(nruns,8)
	tavg_convec_tet = tavg_convec_tet / real(nruns,8)
	tavg_diffu = tavg_diffu / real(nruns,8)
	tavg_diffu_tet = tavg_diffu_tet / real(nruns,8)

	write(*,*)
	write(*,*) 'Timings:'
	write(*,*) '----------------------------------------'
	write(*,*) 'Avg. convective time     = ', tavg_convec
	write(*,*) 'Avg. split convective time = ', tavg_convec_split
	write(*,*) 'Avg. TET convective time = ', tavg_convec_tet
	write(*,*) 'Avg. diffusive time      = ', tavg_diffu
	write(*,*) 'Avg. TET diffusive time  = ', tavg_diffu_tet
	write(*,*) 'Max. convective time     = ', tmax_convec
	write(*,*) 'Max. split convective time = ', tmax_convec_split
	write(*,*) 'Max. TET convective time = ', tmax_convec_tet
	write(*,*) 'Max. diffusive time      = ', tmax_diffu
	write(*,*) 'Max. TET diffusive time  = ', tmax_diffu_tet
	write(*,*) 'Min. convective time     = ', tmin_convec
	write(*,*) 'Min. split convective time     = ', tmin_convec_split
	write(*,*) 'Min. TET convective time = ', tmin_convec_tet
	write(*,*) 'Min. diffusive time      = ', tmin_diffu
	write(*,*) 'Min. TET diffusive time  = ', tmin_diffu_tet
	write(*,*) '----------------------------------------'
	write(*,*) 'Variation convec.        = ', (tmax_convec-tmin_convec)/tavg_convec
	write(*,*) 'Variation split convec.  = ', (tmax_convec_split-tmin_convec_split)/tavg_convec_split
	write(*,*) 'Variation TET convec.    = ', (tmax_convec_tet-tmin_convec_tet)/tavg_convec_tet
	write(*,*) 'Variation diffu.         = ', (tmax_diffu-tmin_diffu)/tavg_diffu
	write(*,*) 'Variation TET diffu.     = ', (tmax_diffu_tet-tmin_diffu_tet)/tavg_diffu_tet
	write(*,*) '----------------------------------------'

	!
	! Print minimal results set
	!
	call nvtxStartRange("Update host results")
	!$acc update host(Rmass,Rmom,Rener,Dmass,Dmom,Dener)
	!$acc update host(Rmass_s,Rmom_s,Rener_s)
	!$acc update host(Rmass_t,Rmom_t,Rener_t,Dmass_t,Dmom_t,Dener_t)
	call nvtxEndRange
	write(*,*)
	write(*,*) 'Basic results:'
	write(*,*) '----------------------------------------'
	write(*,*) '--| Convec'
	write(*,*) 'Max Rmass     = ', maxval(Rmass)    , 'Min Rmass     = ', minval(Rmass)
	write(*,*) 'Max Rmom(:,1) = ', maxval(Rmom(:,1)), 'Min Rmom(:,1) = ', minval(Rmom(:,1))
	write(*,*) 'Max Rmom(:,2) = ', maxval(Rmom(:,2)), 'Min Rmom(:,2) = ', minval(Rmom(:,2))
	write(*,*) 'Max Rmom(:,3) = ', maxval(Rmom(:,3)), 'Min Rmom(:,3) = ', minval(Rmom(:,3))
	write(*,*) 'Max Rener     = ', maxval(Rener)    , 'Min Rener     = ', minval(Rener)
	write(*,*) '----------------------------------------'
	write(*,*) '--| Convec_split'
	write(*,*) 'Max Rmass     = ', maxval(Rmass_s)    , 'Min Rmass     = ', minval(Rmass_s)
	write(*,*) 'Max Rmom(:,1) = ', maxval(Rmom_s(:,1)), 'Min Rmom(:,1) = ', minval(Rmom_s(:,1))
	write(*,*) 'Max Rmom(:,2) = ', maxval(Rmom_s(:,2)), 'Min Rmom(:,2) = ', minval(Rmom_s(:,2))
	write(*,*) 'Max Rmom(:,3) = ', maxval(Rmom_s(:,3)), 'Min Rmom(:,3) = ', minval(Rmom_s(:,3))
	write(*,*) 'Max Rener     = ', maxval(Rener_s)    , 'Min Rener     = ', minval(Rener_s)
	write(*,*) '----------------------------------------'
	write(*,*) '--| Diffu'
	write(*,*) 'Max Dmass     = ', maxval(Dmass)    , 'Min Dmass     = ', minval(Dmass)
	write(*,*) 'Max Dmom(:,1) = ', maxval(Dmom(:,1)), 'Min Dmom(:,1) = ', minval(Dmom(:,1))
	write(*,*) 'Max Dmom(:,2) = ', maxval(Dmom(:,2)), 'Min Dmom(:,2) = ', minval(Dmom(:,2))
	write(*,*) 'Max Dmom(:,3) = ', maxval(Dmom(:,3)), 'Min Dmom(:,3) = ', minval(Dmom(:,3))
	write(*,*) 'Max Dener     = ', maxval(Dener)    , 'Min Dener     = ', minval(Dener)
	write(*,*) 'Basic TET results:'
	write(*,*) '----------------------------------------'
	write(*,*) '--| Convec'
	write(*,*) 'Max Rmass     = ', maxval(Rmass_t)    , 'Min Rmass     = ', minval(Rmass_t)
	write(*,*) 'Max Rmom(:,1) = ', maxval(Rmom_t(:,1)), 'Min Rmom(:,1) = ', minval(Rmom_t(:,1))
	write(*,*) 'Max Rmom(:,2) = ', maxval(Rmom_t(:,2)), 'Min Rmom(:,2) = ', minval(Rmom_t(:,2))
	write(*,*) 'Max Rmom(:,3) = ', maxval(Rmom_t(:,3)), 'Min Rmom(:,3) = ', minval(Rmom_t(:,3))
	write(*,*) 'Max Rener     = ', maxval(Rener_t)    , 'Min Rener     = ', minval(Rener_t)
	write(*,*) '----------------------------------------'
	write(*,*) '--| Diffu'
	write(*,*) 'Max Dmass     = ', maxval(Dmass_t)    , 'Min Dmass     = ', minval(Dmass_t)
	write(*,*) 'Max Dmom(:,1) = ', maxval(Dmom_t(:,1)), 'Min Dmom(:,1) = ', minval(Dmom_t(:,1))
	write(*,*) 'Max Dmom(:,2) = ', maxval(Dmom_t(:,2)), 'Min Dmom(:,2) = ', minval(Dmom_t(:,2))
	write(*,*) 'Max Dmom(:,3) = ', maxval(Dmom_t(:,3)), 'Min Dmom(:,3) = ', minval(Dmom_t(:,3))
	write(*,*) 'Max Dener     = ', maxval(Dener_t)    , 'Min Dener     = ', minval(Dener_t)

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
