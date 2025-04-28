module mod_numerical_params

    use mod_constants
    
    implicit none

        ! Time integration algorithm
        integer(ip)  :: flag_rk_order=4 
        integer(ip)  :: flag_implicit=0 !Explicit: RK, Implicit: BDF2


        ! LES 
        integer(ip)  :: flag_les=0
        integer(ip)  :: flag_les_ilsa=0
        real(rp) :: c_sgs = 0.025_rp
        real(rp) :: stau   = 0.022_rp
        real(rp) :: T_ilsa = 1.0_rp

        ! Discretization parameters
        integer(ip)  :: flag_solver_type=1    ! 1 = Lumped, 2 = APINV, 3 = CG
        integer(ip)  :: flag_spectralElem=1  ! 0 for Lagrange type, 1 for Chebyshev type
        integer(ip)  :: flag_normalise_entropy=1
        real(rp) :: ce = 0.1_rp   
        real(rp) :: cmax = 0.5_rp 
        real(rp) :: cglob =1.0_rp
        real(rp) :: c_rho =1.0_rp
        real(rp) :: c_ener = 1.0_rp
        real(rp) :: flag_mu_factor=1.0_rp

        ! Implicit solver
        real(rp)    :: pseudo_max_dt = 1e20
        real(rp)    :: pseudo_cfl = 1.0_rp
        real(rp)    :: pseudo_ftau= 8.0_rp
        integer(ip) ::  maxIter=20
        integer(ip) ::  maxIterNonLineal=20
        integer(ip) ::  pseudo_steps=10
        real(8)   ::  tol=1e-4
        integer(ip) :: flag_use_constant_dt = 0
        integer(ip) :: flag_implicit_repeat_dt_if_not_converged = 1

        !
        ! Reference conditions
        !
        real(rp) :: nscbc_u_inf   = 1.0_rp
        real(rp) :: nscbc_p_inf = 1.0_rp
        real(rp) :: nscbc_gamma_inf = 1.0_rp
        real(rp) :: nscbc_c_inf = 1.0_rp
        real(rp) :: nscbc_rho_inf   = 1.0_rp
        real(rp) :: nscbc_Rgas_inf   = 1.0_rp
        real(rp) :: nscbc_T_H   = 293.0_rp
        real(rp) :: nscbc_T_C   = 293.0_rp

        !
        ! Penalisation buffer zone
        !

        logical :: flag_buffer_on = .false.
        logical :: flag_buffer_on_east = .false.
        logical :: flag_buffer_on_west = .false.
        logical :: flag_buffer_on_north = .false.
        logical :: flag_buffer_on_south = .false.
        logical :: flag_buffer_on_top = .false.
        logical :: flag_buffer_on_bottom = .false.

        real(4) :: flag_buffer_e_min = 0.0_rp
        real(4) :: flag_buffer_e_size= 0.0_rp
        real(4) :: flag_buffer_w_min = 0.0_rp
        real(4) :: flag_buffer_w_size = 0.0_rp

        real(4) :: flag_buffer_n_min = 0.0_rp
        real(4) :: flag_buffer_n_size = 0.0_rp
        real(4) :: flag_buffer_s_min = 0.0_rp
        real(4) :: flag_buffer_s_size = 0.0_rp

        real(4) :: flag_buffer_t_min = 0.0_rp
        real(4) :: flag_buffer_t_size = 0.0_rp
        real(4) :: flag_buffer_b_min = 0.0_rp
        real(4) :: flag_buffer_b_size = 0.0_rp

        !
        ! Wall model averaging
        !
        real(rp)    :: period_walave   = 1.0_rp
        integer(ip)  :: flag_walave     = 0
        integer(ip)  :: flag_walex      = 3

end module mod_numerical_params
