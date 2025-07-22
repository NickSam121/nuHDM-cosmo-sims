module poisson_commons 
  use amr_commons
  use poisson_parameters

  !~~~~~~~~~ begin ~~~~~~~~~
  ! THIS IS IMPORTANT:
  !~~~~~~~~~~~~~~~~~~~~~~~~~
  ! Replaced "allocatable" property by "pointer"
  ! These pointer will pointer to either
  ! (i)  phi_newton, etc., or
  ! (ii) phi_mond, etc.
  !
  real(dp),pointer,dimension(:)  ::phi,phi_old       ! Potential
  real(dp),pointer,dimension(:)  ::rho               ! Density
  real(dp),pointer,dimension(:,:)::f                 ! 3-force
  !
  !~~~~~~~~~~ end ~~~~~~~~~~

  real(dp),allocatable,dimension(:)  ::rho_top   ! Density at last CIC level                                 

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Communicator arrays for multigrid levels
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg

  ! Minimum MG level
  integer :: levelmin_mg

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

end module poisson_commons


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!      ____  _    _ __  __  ____  _   _ _____                        _   _
!     / __ \| |  | |  \/  |/ __ \| \ | |  __ \                      | | (_)
!    | |  | | |  | | \  / | |  | |  \| | |  | |      _ __ ___  _   _| |_ _ _ __   ___  ___
!    | |  | | |  | | |\/| | |  | | . ` | |  | |     | '__/ _ \| | | | __| | '_ \ / _ \/ __|
!    | |__| | |__| | |  | | |__| | |\  | |__| |  _  | | | (_) | |_| | |_| | | | |  __/\__ \
!     \___\_\\____/|_|  |_|\____/|_| \_|_____/  (_) |_|  \___/ \__,_|\__|_|_| |_|\___||___/
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!-----------------------------------------------
!    Parameters
!-----------------------------------------------
module mond_parameters
  use amr_parameters
  implicit none
  logical :: mond = .true., Activate_g_ext = .true., Uniform_DM = .false. !To disable dark energy, set Omega_Lambda_0 = 0.0d0 below.
  real(dp) :: a0 = -1              ! a0 in user time units
  real(dp) :: a0_ms2 = 1.2d-10    ! a0 in m/s^2, 1.2d-10.
  !Activate_g_ext and variables below added by IB.
  real(dp) :: H_0_kms_Mpc = 67.81d0 !67.81d0.
  real(dp) :: Omega_Lambda_0 = 0.685d0 !0.685d0
  !real(dp) :: Omega_Lambda_0 = 6.85d0!test
  real(dp) :: g_ext_dir_L_degrees = 276.0d0 !Galactic longitude of the external field direction, 276.0d0.
  real(dp) :: g_ext_dir_b_degrees = 30.0d0 !Galactic latitude of the external field direction, 30.0d0.
  real(dp) :: g_ext_ms2 = 3.6d-12!actual EF
  !real(dp) :: g_ext_ms2 = 2.4d-11!test EF
  
end module mond_parameters


!-----------------------------------------------
!    Common variables and constants
!-----------------------------------------------
module mond_commons 
  use amr_commons
  use mond_parameters
  implicit none

  !-----------------------------------------------
  !    Constants
  !-----------------------------------------------
  real(dp) :: PI = ACOS(-1.0)

  !-----------------------------------------------
  !    Constants that will be initiliazed
  !    in the subroutine init_mond
  !-----------------------------------------------
  real(dp) :: a0_i   		      ! inverse of a0
  real(dp) :: g_ext_dir(3), g_ext_dir_L, g_ext_dir_b, Degree !Added by IB to handle g_ext direction.
  real(dp) :: K_0_ext = 0.0d0, nu_ext = 0.0d0, g_ext = 0.0d0, g_N_ext = 0.0d0, g_N_ext_x, g_N_ext_y, g_N_ext_z !Added by IB to handle g_ext magnitude.
  real(dp) :: H_0, rho_Lambda_eff = 0.0d0, u_Lambda_eff = 0.0d0 !Added by IB to handle dark energy.
  real(dp) :: FOUR_PI_G

  ! Used in compute_pdm_density_at_levelmin
  integer,dimension(1:2,1:3)    ::ref_nbpg	 ! Projects neighbor parent grid (nbpg): (orientation,dim) |--> 1..6
  integer,dimension(1:2,1:2,1:3)::ref_nbc	 ! Projects neighbor cell (nbc): (row,orientation,dim) |--> 1..12
  integer,dimension(0:12,1:8)   ::ggg,hhh

  ! Used in compute_pdm_density_at_fine_levels
  integer,dimension(1:8, -2:2, -2:2, -2:2) :: grid_ijk   !  (i,j,k) with i,j,k=-2..+2  --> neighbor grid 1..27
  integer,dimension(1:8, -2:2, -2:2, -2:2) :: cell_ijk   !  (i,j,k) with i,j,k=-2..+2  --> neighbor grid's cell 1..8

  !-----------------------------------------------
  !    Variables
  !-----------------------------------------------
  real(dp),pointer,dimension(:)   :: rho_mond  ! Baryonic + phantom dark matter density
  real(dp),pointer,dimension(:)   :: phi_mond  ! MONDian potential
  real(dp),pointer,dimension(:)   :: phi_old_mond
  real(dp),pointer,dimension(:,:) :: f_mond    ! MONDian acceleration

  real(dp),pointer,dimension(:)   :: rho_newton  ! Baryonic matter density
  real(dp),pointer,dimension(:)   :: phi_newton  ! Newtonian potential
  real(dp),pointer,dimension(:)   :: phi_old_newton
  real(dp),pointer,dimension(:,:) :: f_newton    ! Newtonian acceleration

  real(dp),target::rho_mond_tot=0.0D0          ! Mean PDM density in the box, needed by the Poisson solver  
  real(dp),target::rho_newton_tot=0.0D0        ! Mean baryonic density in the box, needed by the Poisson solver  

  logical::connected_Mond = .false.

end module mond_commons

!-----------------------------------------------
!    The interpolation function
!-----------------------------------------------
!
!   Routine that defines nu(x)
!   IMPORTANT: Make sure nu(x)->0 if x>>1
!
!-----------------------------------------------
subroutine get_nu(x, nu)
   use amr_commons
   implicit none
   real(dp),intent(in) :: x
   real(dp),intent(out) :: nu

!  Notice again: it is IMPORTANT that nu(x)->0 if x>>1

   nu = sqrt(0.25d0 + 1.0d0/x) - 0.5d0          ! Simple nu function, made more efficient by IB.
!  nu = sqrt( 0.5 + sqrt( 0.25 + 1.0/x**2 ) ) - 1.0 ! Standard nu function
!  nu = 1.0/(1.0 - exp(-sqrt(x))) - 1.0             ! ...
end subroutine get_nu

subroutine get_nu_K_0(x, nu, K_0) !Added by IB
    use amr_commons, only : dp
    implicit none
    real(dp), intent(in) :: x
    real(dp), intent(out) :: nu, K_0

    !nu is calculated properly this time, so g = nu*g_N.
    nu = sqrt(0.25d0 + 1.0d0/x) + 0.5d0    
    K_0 = -0.5d0/(x*(nu - 0.5d0)*nu)
end subroutine get_nu_K_0

subroutine Sph_xyz(L, b, v) !Added by IB
    use amr_commons, only : dp
    implicit none
    real(dp), intent(in) :: L, b !Must be in radians.
    real(dp), intent(out) :: v(3)
    real(dp) :: cos_b, nu
    cos_b = cos(b)
    v(1) = cos_b*cos(L)
    v(2) = cos_b*sin(L)
    v(3) = sin(b)
end subroutine Sph_xyz

subroutine xyz_Sph(x, y, z, L, b) !Added by IB
    use amr_commons, only : dp
    use mond_commons, only : pi
    implicit none
    real(dp), intent(in) :: x, y, z
    real(dp), intent(out) :: L, b
    real(dp) :: Radian, r_sq
    Radian = 180.0d0/pi;
    r_sq = x*x + y*y;
    b = Radian*asin(z/sqrt(r_sq + z*z));
    if (abs(b) > 89.999999) then
        L = 0.0d0;
    else
        L = Radian*asin(y/sqrt(r_sq)); !Correct only in quadrant 1.
        if (x < 0.0d0) then
            L = 180.0d0 - L; !Quadrants 2 or 3.
        elseif (y < 0.0d0) then
            L = L + 360.0d0; !Quadrant 4.
        end if
    end if
end subroutine xyz_Sph

subroutine g_N_determination(g, g_N, nu, K_0) !Added by IB
    use amr_commons, only : dp, myid
    implicit none
    real(dp), intent(in) :: g
    real(dp), intent(out) :: g_N, nu, K_0
    integer*1 :: i
    real(dp) :: g_N_1, g_N_2, Gradient_inv
    real(dp) :: Error, Error_1, Error_2, Tolerance

    if (g > 1.0d0) then
        g_N_1 = g
    elseif (g > 1.0d-9) then
        g_N_1 = g*g
    else
        g_N = 0
        nu = 0
        K_0 = 0
        write(6, *) 'g_ext assumed disabled, proceeding with g_N = nu = K_0 = 0.'
        return
    end if
    g_N_2 = g_N_1*1.1d0
    Tolerance = g*1.0d-9

    do i = 1, 97
        if (i .eq. 1) then
            g_N = g_N_1
        else if (i .eq. 2) then
            g_N = g_N_2
        else
            Gradient_inv = (g_N_2 - g_N_1)/(Error_2 - Error_1)
            g_N_1 = g_N_2
            Error_1 = Error_2
            g_N_2 = g_N_2 - Gradient_inv*Error_2
            g_N = g_N_2
        end if

        call get_nu_K_0(g_N, nu, K_0)
        Error = nu*g_N - g
        !write(6, *) 'Step i = ', i, ', error in g/a_0 = ', Error, '.'
        if (i .eq. 1) then
            Error_1 = Error
        else
            Error_2 = Error
        end if

        if (abs(Error) < Tolerance) then
            if (i .ge. 2) then
                g_N = g_N - Error*Gradient_inv
            end if
            call get_nu_K_0(g_N, nu, K_0)
            if (myid==1) write(6, *) 'Algorithm converged on step i = ', i, ' out of 97 available, resulting g_N = ', g_N, ' a_0.'
            exit
        end if
     end do
     if ((i .eq. 97) .and. abs(Error) > Tolerance) then
        write(*,*) 'Newton-Raphson failed.'
     endif
end subroutine g_N_determination

!The subroutine output_centres was here.

!https://stackoverflow.com/questions/7876075/getting-free-unit-number-in-fortran
!      integer*4 function get_file_unit (lu_max)
!!
!!   get_file_unit returns a unit number that is not in use
!      integer*4 lu_max,  lu, m, iostat
!      logical   opened
!!
!      m = lu_max  ;  if (m < 1) m = 97
!      do lu = m,1,-1
!         inquire (unit=lu, opened=opened, iostat=iostat)
!         if (iostat.ne.0) cycle
!         if (.not.opened) exit
!      end do
!!
!      get_file_unit = lu
!      return
!      end function get_file_unit

!Can probably also use unit_output_centres, defined in pm_commons (currently commented out). This is a random number chosen via calculator of IT. Likelihood of being in use is very low.

!how to append data to a file on channel 61.
!OPEN(61,file=ofile,action='write',position='append')


!-----------------------------------------------
!    Initialization of static variables
!    This routine is called from adaptive_loop.f90
!-----------------------------------------------
subroutine init_mond
    use amr_commons
    use mond_commons
    use poisson_commons
    implicit none
    integer :: i,j,k,ncell

    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

    !Reset to aexp = 1.0d0.
    scale_l = scale_l/aexp
    scale_t = scale_t/aexp**2
    scale_d = scale_d*aexp**3
    scale_nH = scale_nH*aexp**3
    scale_v = scale_v*aexp
    scale_T2 = scale_T2*aexp**2

    if (ndim < 3) then
        write(*,*) ' ERROR: the "phantom" patch requires NDIM=3'
        call clean_stop
    endif

    ! Allocate additional arrays
    ! In analogy by init_amr.f90
    ncell=ncoarse+twotondim*ngridmax
    allocate(rho_mond    (1:ncell))
    allocate(phi_mond    (1:ncell))
    allocate(phi_old_mond(1:ncell))
    allocate( f_mond (1:ncell,1:3))
    rho_mond=0.0D0
    phi_mond=0.0D0
    f_mond=0.0D0

    ! "Connect" also the Newtonian arrays,
    ! assuming they have already been initialized
    ! in the standard Ramses routines (see init_poisson.f90)
    rho_newton     => rho
    phi_newton     => phi
    phi_old_newton => phi_old
    f_newton       => f

    ! Scale a0 if it is provided in SI units
    if (a0<0 .and. a0_ms2>0) then
        a0 = a0_ms2*100.d0 / scale_l * scale_t**2  ! converts from [cm s^-2]  into  [kpc / [user time unit]^2]
    elseif (a0_ms2<0 .and. a0<0 .and. mond) then
        a0_ms2 = a0/100.d0 * scale_l / scale_t**2
    elseif (myid==1) then 
	write(*,*) ' ERROR: a0 parameter missing'
        call clean_stop
    endif

    !a0_i = 1.d0 / a0 !Disabled by IB as done below.
    if (myid==1) then
      write(*,'(" Initializing MOND extension: a0 = ", E11.4, " m/s^2")') a0_ms2
      write(*,'("                                 = ", E11.4, " [user length unit]/[user time unit]^2")') a0
    endif


    FOUR_PI_G  = 4.0d0*PI*1.0d0  ! G=1 in the Ramses Poisson solver

    if (cosmo) then !Block added by IB for cosmo runs.
        FOUR_PI_G = 1.5D0*omega_m*aexp
        !a0_i = 1.0d0/(aexp*a0), unused in this subroutine and recalculated later.
    else
        FOUR_PI_G  = 4.0d0*PI*1.0d0
        a0_i = 1.0d0/a0 !Not recalculated later.
    end if

    !Block added by IB to handle g_ext.
    Degree = pi/180.0d0
    g_ext_dir_L = g_ext_dir_L_degrees*Degree
    g_ext_dir_b = g_ext_dir_b_degrees*Degree
    call Sph_xyz(g_ext_dir_L, g_ext_dir_b, g_ext_dir)
    !write(6, *) 'g_ext_dir = ',g_ext_dir
    if (Activate_g_ext) then
        g_ext = g_ext_ms2/a0_ms2
        call g_N_determination(g_ext, g_N_ext, nu_ext, K_0_ext)
        g_ext = g_ext*a0
        g_N_ext = g_N_ext*a0
    end if
    g_N_ext_x = g_N_ext*g_ext_dir(1)
    g_N_ext_y = g_N_ext*g_ext_dir(2)
    g_N_ext_z = g_N_ext*g_ext_dir(3)

    if (cosmo) then
      rho_Lambda_eff = 0.0d0
      u_Lambda_eff = 0.0d0
    else !Block added by IB to handle dark energy.
      H_0 = H_0_kms_Mpc*1.0d5/3.08568d24*scale_t
      rho_Lambda_eff = -3.0d0/Four_pi_G*H_0*H_0*Omega_Lambda_0
      u_Lambda_eff = -0.5d0*H_0*H_0*Omega_Lambda_0
    end if

    !  The follwing arrays (ref_nbc, ref_nbpg, ggg, hhh) are used in
    !  compute_pdm_density_at_levelmin

    ! Reference pointer of neighbor cells
    ref_nbc(1,1,1) = 1
    ref_nbc(1,1,2) = 2
    ref_nbc(1,1,3) = 3
    ref_nbc(1,2,1) = 4
    ref_nbc(1,2,2) = 5
    ref_nbc(1,2,3) = 6
    ref_nbc(2,1,1) = 7
    ref_nbc(2,1,2) = 8
    ref_nbc(2,1,3) = 9
    ref_nbc(2,2,1) = 10
    ref_nbc(2,2,2) = 11
    ref_nbc(2,2,3) = 12
    ! Reference pointer of neighbor parent grids
    ref_nbpg(1,1) = 1   ! left, x
    ref_nbpg(1,2) = 2   ! right, x
    ref_nbpg(1,3) = 3   ! left, y
    ref_nbpg(2,1) = 4   ! right, y
    ref_nbpg(2,2) = 5   ! left, z
    ref_nbpg(2,3) = 6   ! right, z

    ggg(      0,        1:8) = 0
    ggg(ref_nbc(1,1,1), 1:8) = (/ ref_nbpg(1,1), 0, ref_nbpg(1,1), 0, ref_nbpg(1,1), 0, ref_nbpg(1,1), 0  /)
    ggg(ref_nbc(1,1,2), 1:8) = (/ ref_nbpg(1,2), ref_nbpg(1,2), 0, 0, ref_nbpg(1,2), ref_nbpg(1,2), 0, 0  /)
    ggg(ref_nbc(1,1,3), 1:8) = (/ ref_nbpg(1,3), ref_nbpg(1,3), ref_nbpg(1,3), ref_nbpg(1,3), 0, 0, 0, 0  /)
    ggg(ref_nbc(1,2,1), 1:8) = (/ 0, ref_nbpg(2,1), 0, ref_nbpg(2,1), 0, ref_nbpg(2,1), 0, ref_nbpg(2,1)  /)
    ggg(ref_nbc(1,2,2), 1:8) = (/ 0, 0, ref_nbpg(2,2), ref_nbpg(2,2), 0, 0, ref_nbpg(2,2), ref_nbpg(2,2)  /)
    ggg(ref_nbc(1,2,3), 1:8) = (/ 0, 0, 0, 0, ref_nbpg(2,3), ref_nbpg(2,3), ref_nbpg(2,3), ref_nbpg(2,3)  /)
    do i=1,2
        do j=1,3
            ggg(ref_nbc(2,i,j), 1:8) = ref_nbpg(i,j) 
        enddo
    enddo

    hhh(      0,        1:8) = (/ 1,2,3,4, 5,6,7,8  /)
    hhh(ref_nbc(1,1,1), 1:8) = (/ 2,1,4,3, 6,5,8,7  /)
    hhh(ref_nbc(1,1,2), 1:8) = (/ 3,4,1,2, 7,8,5,6  /)
    hhh(ref_nbc(1,1,3), 1:8) = (/ 5,6,7,8, 1,2,3,4  /)
    hhh(ref_nbc(1,2,1), 1:8) = (/ 2,1,4,3, 6,5,8,7  /)
    hhh(ref_nbc(1,2,2), 1:8) = (/ 3,4,1,2, 7,8,5,6  /)
    hhh(ref_nbc(1,2,3), 1:8) = (/ 5,6,7,8, 1,2,3,4  /)
    hhh(ref_nbc(2,1,1), 1:8) = (/ 1,2,3,4, 5,6,7,8  /)
    hhh(ref_nbc(2,1,2), 1:8) = (/ 1,2,3,4, 5,6,7,8  /)
    hhh(ref_nbc(2,1,3), 1:8) = (/ 1,2,3,4, 5,6,7,8  /)
    hhh(ref_nbc(2,2,1), 1:8) = (/ 1,2,3,4, 5,6,7,8  /)
    hhh(ref_nbc(2,2,2), 1:8) = (/ 1,2,3,4, 5,6,7,8  /)
    hhh(ref_nbc(2,2,3), 1:8) = (/ 1,2,3,4, 5,6,7,8  /)


    !  The follwing arrays (cell_ijk and grid_ijk) are used in
    !  compute_pdm_density_at_fine_levels

    !!!!! Indices to the 27x8 neighbor cells !!!!!
    cell_ijk(1:8,-2,-2,-2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8,-2,-2,-1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8,-2,-2, 0) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8,-2,-2, 1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8,-2,-2, 2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8,-2,-1,-2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8,-2,-1,-1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8,-2,-1, 0) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8,-2,-1, 1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8,-2,-1, 2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8,-2, 0,-2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8,-2, 0,-1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8,-2, 0, 0) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8,-2, 0, 1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8,-2, 0, 2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8,-2, 1,-2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8,-2, 1,-1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8,-2, 1, 0) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8,-2, 1, 1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8,-2, 1, 2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8,-2, 2,-2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8,-2, 2,-1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8,-2, 2, 0) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8,-2, 2, 1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8,-2, 2, 2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8,-1,-2,-2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8,-1,-2,-1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8,-1,-2, 0) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8,-1,-2, 1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8,-1,-2, 2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8,-1,-1,-2) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8,-1,-1,-1) = (/8,7,6,5,4,3,2,1/)
    cell_ijk(1:8,-1,-1, 0) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8,-1,-1, 1) = (/8,7,6,5,4,3,2,1/)
    cell_ijk(1:8,-1,-1, 2) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8,-1, 0,-2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8,-1, 0,-1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8,-1, 0, 0) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8,-1, 0, 1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8,-1, 0, 2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8,-1, 1,-2) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8,-1, 1,-1) = (/8,7,6,5,4,3,2,1/)
    cell_ijk(1:8,-1, 1, 0) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8,-1, 1, 1) = (/8,7,6,5,4,3,2,1/)
    cell_ijk(1:8,-1, 1, 2) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8,-1, 2,-2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8,-1, 2,-1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8,-1, 2, 0) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8,-1, 2, 1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8,-1, 2, 2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8, 0,-2,-2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 0,-2,-1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 0,-2, 0) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 0,-2, 1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 0,-2, 2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 0,-1,-2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 0,-1,-1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8, 0,-1, 0) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 0,-1, 1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8, 0,-1, 2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 0, 0,-2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 0, 0,-1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 0, 0, 0) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 0, 0, 1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 0, 0, 2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 0, 1,-2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 0, 1,-1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8, 0, 1, 0) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 0, 1, 1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8, 0, 1, 2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 0, 2,-2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 0, 2,-1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 0, 2, 0) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 0, 2, 1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 0, 2, 2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 1,-2,-2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8, 1,-2,-1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8, 1,-2, 0) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8, 1,-2, 1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8, 1,-2, 2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8, 1,-1,-2) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8, 1,-1,-1) = (/8,7,6,5,4,3,2,1/)
    cell_ijk(1:8, 1,-1, 0) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8, 1,-1, 1) = (/8,7,6,5,4,3,2,1/)
    cell_ijk(1:8, 1,-1, 2) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8, 1, 0,-2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8, 1, 0,-1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8, 1, 0, 0) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8, 1, 0, 1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8, 1, 0, 2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8, 1, 1,-2) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8, 1, 1,-1) = (/8,7,6,5,4,3,2,1/)
    cell_ijk(1:8, 1, 1, 0) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8, 1, 1, 1) = (/8,7,6,5,4,3,2,1/)
    cell_ijk(1:8, 1, 1, 2) = (/4,3,2,1,8,7,6,5/)
    cell_ijk(1:8, 1, 2,-2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8, 1, 2,-1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8, 1, 2, 0) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8, 1, 2, 1) = (/6,5,8,7,2,1,4,3/)
    cell_ijk(1:8, 1, 2, 2) = (/2,1,4,3,6,5,8,7/)
    cell_ijk(1:8, 2,-2,-2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 2,-2,-1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 2,-2, 0) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 2,-2, 1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 2,-2, 2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 2,-1,-2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 2,-1,-1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8, 2,-1, 0) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 2,-1, 1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8, 2,-1, 2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 2, 0,-2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 2, 0,-1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 2, 0, 0) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 2, 0, 1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 2, 0, 2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 2, 1,-2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 2, 1,-1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8, 2, 1, 0) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 2, 1, 1) = (/7,8,5,6,3,4,1,2/)
    cell_ijk(1:8, 2, 1, 2) = (/3,4,1,2,7,8,5,6/)
    cell_ijk(1:8, 2, 2,-2) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 2, 2,-1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 2, 2, 0) = (/1,2,3,4,5,6,7,8/)
    cell_ijk(1:8, 2, 2, 1) = (/5,6,7,8,1,2,3,4/)
    cell_ijk(1:8, 2, 2, 2) = (/1,2,3,4,5,6,7,8/)


    grid_ijk(1:8,-2,-2,-2) = (/ 1, 1, 1, 1, 1, 1, 1, 1/)
    grid_ijk(1:8,-2,-2,-1) = (/ 1, 1, 1, 1,10,10,10,10/)
    grid_ijk(1:8,-2,-2, 0) = (/10,10,10,10,10,10,10,10/)
    grid_ijk(1:8,-2,-2, 1) = (/10,10,10,10,19,19,19,19/)
    grid_ijk(1:8,-2,-2, 2) = (/19,19,19,19,19,19,19,19/)
    grid_ijk(1:8,-2,-1,-2) = (/ 1, 1, 4, 4, 1, 1, 4, 4/)
    grid_ijk(1:8,-2,-1,-1) = (/ 1, 1, 4, 4,10,10,13,13/)
    grid_ijk(1:8,-2,-1, 0) = (/10,10,13,13,10,10,13,13/)
    grid_ijk(1:8,-2,-1, 1) = (/10,10,13,13,19,19,22,22/)
    grid_ijk(1:8,-2,-1, 2) = (/19,19,22,22,19,19,22,22/)
    grid_ijk(1:8,-2, 0,-2) = (/ 4, 4, 4, 4, 4, 4, 4, 4/)
    grid_ijk(1:8,-2, 0,-1) = (/ 4, 4, 4, 4,13,13,13,13/)
    grid_ijk(1:8,-2, 0, 0) = (/13,13,13,13,13,13,13,13/)
    grid_ijk(1:8,-2, 0, 1) = (/13,13,13,13,22,22,22,22/)
    grid_ijk(1:8,-2, 0, 2) = (/22,22,22,22,22,22,22,22/)
    grid_ijk(1:8,-2, 1,-2) = (/ 4, 4, 7, 7, 4, 4, 7, 7/)
    grid_ijk(1:8,-2, 1,-1) = (/ 4, 4, 7, 7,13,13,16,16/)
    grid_ijk(1:8,-2, 1, 0) = (/13,13,16,16,13,13,16,16/)
    grid_ijk(1:8,-2, 1, 1) = (/13,13,16,16,22,22,25,25/)
    grid_ijk(1:8,-2, 1, 2) = (/22,22,25,25,22,22,25,25/)
    grid_ijk(1:8,-2, 2,-2) = (/ 7, 7, 7, 7, 7, 7, 7, 7/)
    grid_ijk(1:8,-2, 2,-1) = (/ 7, 7, 7, 7,16,16,16,16/)
    grid_ijk(1:8,-2, 2, 0) = (/16,16,16,16,16,16,16,16/)
    grid_ijk(1:8,-2, 2, 1) = (/16,16,16,16,25,25,25,25/)
    grid_ijk(1:8,-2, 2, 2) = (/25,25,25,25,25,25,25,25/)
    grid_ijk(1:8,-1,-2,-2) = (/ 1, 2, 1, 2, 1, 2, 1, 2/)
    grid_ijk(1:8,-1,-2,-1) = (/ 1, 2, 1, 2,10,11,10,11/)
    grid_ijk(1:8,-1,-2, 0) = (/10,11,10,11,10,11,10,11/)
    grid_ijk(1:8,-1,-2, 1) = (/10,11,10,11,19,20,19,20/)
    grid_ijk(1:8,-1,-2, 2) = (/19,20,19,20,19,20,19,20/)
    grid_ijk(1:8,-1,-1,-2) = (/ 1, 2, 4, 5, 1, 2, 4, 5/)
    grid_ijk(1:8,-1,-1,-1) = (/ 1, 2, 4, 5,10,11,13,14/)
    grid_ijk(1:8,-1,-1, 0) = (/10,11,13,14,10,11,13,14/)
    grid_ijk(1:8,-1,-1, 1) = (/10,11,13,14,19,20,22,23/)
    grid_ijk(1:8,-1,-1, 2) = (/19,20,22,23,19,20,22,23/)
    grid_ijk(1:8,-1, 0,-2) = (/ 4, 5, 4, 5, 4, 5, 4, 5/)
    grid_ijk(1:8,-1, 0,-1) = (/ 4, 5, 4, 5,13,14,13,14/)
    grid_ijk(1:8,-1, 0, 0) = (/13,14,13,14,13,14,13,14/)
    grid_ijk(1:8,-1, 0, 1) = (/13,14,13,14,22,23,22,23/)
    grid_ijk(1:8,-1, 0, 2) = (/22,23,22,23,22,23,22,23/)
    grid_ijk(1:8,-1, 1,-2) = (/ 4, 5, 7, 8, 4, 5, 7, 8/)
    grid_ijk(1:8,-1, 1,-1) = (/ 4, 5, 7, 8,13,14,16,17/)
    grid_ijk(1:8,-1, 1, 0) = (/13,14,16,17,13,14,16,17/)
    grid_ijk(1:8,-1, 1, 1) = (/13,14,16,17,22,23,25,26/)
    grid_ijk(1:8,-1, 1, 2) = (/22,23,25,26,22,23,25,26/)
    grid_ijk(1:8,-1, 2,-2) = (/ 7, 8, 7, 8, 7, 8, 7, 8/)
    grid_ijk(1:8,-1, 2,-1) = (/ 7, 8, 7, 8,16,17,16,17/)
    grid_ijk(1:8,-1, 2, 0) = (/16,17,16,17,16,17,16,17/)
    grid_ijk(1:8,-1, 2, 1) = (/16,17,16,17,25,26,25,26/)
    grid_ijk(1:8,-1, 2, 2) = (/25,26,25,26,25,26,25,26/)
    grid_ijk(1:8, 0,-2,-2) = (/ 2, 2, 2, 2, 2, 2, 2, 2/)
    grid_ijk(1:8, 0,-2,-1) = (/ 2, 2, 2, 2,11,11,11,11/)
    grid_ijk(1:8, 0,-2, 0) = (/11,11,11,11,11,11,11,11/)
    grid_ijk(1:8, 0,-2, 1) = (/11,11,11,11,20,20,20,20/)
    grid_ijk(1:8, 0,-2, 2) = (/20,20,20,20,20,20,20,20/)
    grid_ijk(1:8, 0,-1,-2) = (/ 2, 2, 5, 5, 2, 2, 5, 5/)
    grid_ijk(1:8, 0,-1,-1) = (/ 2, 2, 5, 5,11,11,14,14/)
    grid_ijk(1:8, 0,-1, 0) = (/11,11,14,14,11,11,14,14/)
    grid_ijk(1:8, 0,-1, 1) = (/11,11,14,14,20,20,23,23/)
    grid_ijk(1:8, 0,-1, 2) = (/20,20,23,23,20,20,23,23/)
    grid_ijk(1:8, 0, 0,-2) = (/ 5, 5, 5, 5, 5, 5, 5, 5/)
    grid_ijk(1:8, 0, 0,-1) = (/ 5, 5, 5, 5,14,14,14,14/)
    grid_ijk(1:8, 0, 0, 0) = (/14,14,14,14,14,14,14,14/)
    grid_ijk(1:8, 0, 0, 1) = (/14,14,14,14,23,23,23,23/)
    grid_ijk(1:8, 0, 0, 2) = (/23,23,23,23,23,23,23,23/)
    grid_ijk(1:8, 0, 1,-2) = (/ 5, 5, 8, 8, 5, 5, 8, 8/)
    grid_ijk(1:8, 0, 1,-1) = (/ 5, 5, 8, 8,14,14,17,17/)
    grid_ijk(1:8, 0, 1, 0) = (/14,14,17,17,14,14,17,17/)
    grid_ijk(1:8, 0, 1, 1) = (/14,14,17,17,23,23,26,26/)
    grid_ijk(1:8, 0, 1, 2) = (/23,23,26,26,23,23,26,26/)
    grid_ijk(1:8, 0, 2,-2) = (/ 8, 8, 8, 8, 8, 8, 8, 8/)
    grid_ijk(1:8, 0, 2,-1) = (/ 8, 8, 8, 8,17,17,17,17/)
    grid_ijk(1:8, 0, 2, 0) = (/17,17,17,17,17,17,17,17/)
    grid_ijk(1:8, 0, 2, 1) = (/17,17,17,17,26,26,26,26/)
    grid_ijk(1:8, 0, 2, 2) = (/26,26,26,26,26,26,26,26/)
    grid_ijk(1:8, 1,-2,-2) = (/ 2, 3, 2, 3, 2, 3, 2, 3/)
    grid_ijk(1:8, 1,-2,-1) = (/ 2, 3, 2, 3,11,12,11,12/)
    grid_ijk(1:8, 1,-2, 0) = (/11,12,11,12,11,12,11,12/)
    grid_ijk(1:8, 1,-2, 1) = (/11,12,11,12,20,21,20,21/)
    grid_ijk(1:8, 1,-2, 2) = (/20,21,20,21,20,21,20,21/)
    grid_ijk(1:8, 1,-1,-2) = (/ 2, 3, 5, 6, 2, 3, 5, 6/)
    grid_ijk(1:8, 1,-1,-1) = (/ 2, 3, 5, 6,11,12,14,15/)
    grid_ijk(1:8, 1,-1, 0) = (/11,12,14,15,11,12,14,15/)
    grid_ijk(1:8, 1,-1, 1) = (/11,12,14,15,20,21,23,24/)
    grid_ijk(1:8, 1,-1, 2) = (/20,21,23,24,20,21,23,24/)
    grid_ijk(1:8, 1, 0,-2) = (/ 5, 6, 5, 6, 5, 6, 5, 6/)
    grid_ijk(1:8, 1, 0,-1) = (/ 5, 6, 5, 6,14,15,14,15/)
    grid_ijk(1:8, 1, 0, 0) = (/14,15,14,15,14,15,14,15/)
    grid_ijk(1:8, 1, 0, 1) = (/14,15,14,15,23,24,23,24/)
    grid_ijk(1:8, 1, 0, 2) = (/23,24,23,24,23,24,23,24/)
    grid_ijk(1:8, 1, 1,-2) = (/ 5, 6, 8, 9, 5, 6, 8, 9/)
    grid_ijk(1:8, 1, 1,-1) = (/ 5, 6, 8, 9,14,15,17,18/)
    grid_ijk(1:8, 1, 1, 0) = (/14,15,17,18,14,15,17,18/)
    grid_ijk(1:8, 1, 1, 1) = (/14,15,17,18,23,24,26,27/)
    grid_ijk(1:8, 1, 1, 2) = (/23,24,26,27,23,24,26,27/)
    grid_ijk(1:8, 1, 2,-2) = (/ 8, 9, 8, 9, 8, 9, 8, 9/)
    grid_ijk(1:8, 1, 2,-1) = (/ 8, 9, 8, 9,17,18,17,18/)
    grid_ijk(1:8, 1, 2, 0) = (/17,18,17,18,17,18,17,18/)
    grid_ijk(1:8, 1, 2, 1) = (/17,18,17,18,26,27,26,27/)
    grid_ijk(1:8, 1, 2, 2) = (/26,27,26,27,26,27,26,27/)
    grid_ijk(1:8, 2,-2,-2) = (/ 3, 3, 3, 3, 3, 3, 3, 3/)
    grid_ijk(1:8, 2,-2,-1) = (/ 3, 3, 3, 3,12,12,12,12/)
    grid_ijk(1:8, 2,-2, 0) = (/12,12,12,12,12,12,12,12/)
    grid_ijk(1:8, 2,-2, 1) = (/12,12,12,12,21,21,21,21/)
    grid_ijk(1:8, 2,-2, 2) = (/21,21,21,21,21,21,21,21/)
    grid_ijk(1:8, 2,-1,-2) = (/ 3, 3, 6, 6, 3, 3, 6, 6/)
    grid_ijk(1:8, 2,-1,-1) = (/ 3, 3, 6, 6,12,12,15,15/)
    grid_ijk(1:8, 2,-1, 0) = (/12,12,15,15,12,12,15,15/)
    grid_ijk(1:8, 2,-1, 1) = (/12,12,15,15,21,21,24,24/)
    grid_ijk(1:8, 2,-1, 2) = (/21,21,24,24,21,21,24,24/)
    grid_ijk(1:8, 2, 0,-2) = (/ 6, 6, 6, 6, 6, 6, 6, 6/)
    grid_ijk(1:8, 2, 0,-1) = (/ 6, 6, 6, 6,15,15,15,15/)
    grid_ijk(1:8, 2, 0, 0) = (/15,15,15,15,15,15,15,15/)
    grid_ijk(1:8, 2, 0, 1) = (/15,15,15,15,24,24,24,24/)
    grid_ijk(1:8, 2, 0, 2) = (/24,24,24,24,24,24,24,24/)
    grid_ijk(1:8, 2, 1,-2) = (/ 6, 6, 9, 9, 6, 6, 9, 9/)
    grid_ijk(1:8, 2, 1,-1) = (/ 6, 6, 9, 9,15,15,18,18/)
    grid_ijk(1:8, 2, 1, 0) = (/15,15,18,18,15,15,18,18/)
    grid_ijk(1:8, 2, 1, 1) = (/15,15,18,18,24,24,27,27/)
    grid_ijk(1:8, 2, 1, 2) = (/24,24,27,27,24,24,27,27/)
    grid_ijk(1:8, 2, 2,-2) = (/ 9, 9, 9, 9, 9, 9, 9, 9/)
    grid_ijk(1:8, 2, 2,-1) = (/ 9, 9, 9, 9,18,18,18,18/)
    grid_ijk(1:8, 2, 2, 0) = (/18,18,18,18,18,18,18,18/)
    grid_ijk(1:8, 2, 2, 1) = (/18,18,18,18,27,27,27,27/)
    grid_ijk(1:8, 2, 2, 2) = (/27,27,27,27,27,27,27,27/)

end subroutine init_mond


!-----------------------------------------------
!  Helper routines to switch between 
!  Newtonian/MONDian arrays
!-----------------------------------------------
subroutine connect_Newton()
   use poisson_commons
   use mond_commons
   implicit none
   
   rho     => rho_newton
   phi     => phi_newton
   phi_old => phi_old_newton
   f       => f_newton
   rho_tot => rho_newton_tot
   
   connected_Mond = .false.
end subroutine connect_Newton

subroutine connect_Mond()
   use poisson_commons
   use mond_commons
   implicit none
   
   rho     => rho_mond
   phi     => phi_mond
   phi_old => phi_old_mond
   f       => f_mond
   rho_tot => rho_mond_tot

   connected_Mond = .true.
end subroutine connect_Mond


!-----------------------------------------------
!  Computes rho_mond_tot at level levelmin
!  This is the mean 3-dim. density of PDM 
!  in the intire simulation box.
!  rho_mond_tot is needed by the Poisson solver.
!-----------------------------------------------
subroutine compute_rho_mond_tot(ilevel)
  use amr_commons
  use poisson_commons
  use mond_commons
  implicit none

  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp)::dx,dx3
  integer::ix,iy,iz
  integer::i,istart,ncache,ind,ilevel,igrid,icell,idim,info,count,count_all
  integer,allocatable,dimension(:)::ind_grid
  real(dp)::rmt,rmt_all
  integer,dimension(1:twotondim) :: iskip

#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  if (ilevel .ne. levelmin) then
     return
  endif
  
  if(nboundary>0)then
     rho_tot=0d0
     return
  endif
  
  do ind=1,twotondim
     iskip(ind)=ncoarse+(ind-1)*ngridmax
  enddo
   
  rmt = 0.0
  count = 0

  ! Set position of cell centers relative to grid center
  ! at each level ilevel
  dx=0.5D0**ilevel
  do ind=1,twotondim
      iz=(ind-1)/4
      iy=(ind-1-4*iz)/2
      ix=(ind-1-2*iy-4*iz)
      if(ndim>0) xc(ind,1)=(dble(ix)-0.5D0)*dx
      if(ndim>1) xc(ind,2)=(dble(iy)-0.5D0)*dx
      if(ndim>2) xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do  

  istart=headl(myid,ilevel)
  ncache=numbl(myid,ilevel)

  if(ncache>0)then
     allocate(ind_grid(1:ncache))
     ! Loop over level grids
     igrid=istart
     do i=1,ncache
        ind_grid(i)=igrid
        igrid=next(igrid)
     end do

     ! Loop over grids
     !$____omp parallel do private(ind) reduction(+:rmt) reduction(+:count)
     do i=1,ncache
        ! Loop over cells
        do ind=1,twotondim
           rmt = rmt + rho_mond(iskip(ind) + ind_grid(i))
           count = count + 1
        enddo
     end do
     !$____omp end parallel do
     
     deallocate(ind_grid)
  end if
  
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rmt,   rmt_all,   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(count, count_all, 1, MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
  rho_mond_tot = rmt_all/count_all
#else
  rho_mond_tot = rmt/count
#endif

  if (verbose .and. myid==1) write(*,*) 'rho_mond_tot =',rho_mond_tot

end subroutine


!-----------------------------------------------
!  Compute the PDM density at level ilevel
!  Switches between different routines
!  depending on the grid level
!-----------------------------------------------
subroutine compute_pdm_density(ilevel,icount)
   use amr_commons
   implicit none

   integer::ilevel,icount

   if (ilevel==levelmin) then
      !   compute_pdm_density_at_levelmin:  
      !     Makes use of the already computed acceleration
      !     Used at ilevel==levelmin
      call compute_pdm_density_at_levelmin(ilevel,icount)
   else
      !   compute_pdm_density_at_fine_levels:  
      !     Uses a five-point finite difference approximation
      !     Interpolates the potential at the level boundaries
      !     Used at ilevel>levelmin
      call compute_pdm_density_at_fine_levels(ilevel,icount)
   endif

end subroutine compute_pdm_density


!#########################################################
!##
!##   Compute the PDM density at level levelmin
!##
!#########################################################
!##
!##     Makes use of the already compute acceleration
!##     Used at ilevel==levelmin
!##
!##     The routines compute_pdm_density_at_levelmin and 
!##     compute_pdm_density_at_fine_levels are effectively
!##     equivalent. The problem is that the Newtonian 
!##     potential is not computed for the diagonal boundary 
!##     cells at the coarsest level. Easiest way to avoid
!##     inconveniences is to use compute_pdm_density_at_levelmin at levelmin
!##     and compute_pdm_density_at_levelmin at ilevel>levelmin.
!##
!#########################################################
subroutine compute_pdm_density_at_levelmin(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  use mond_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel,icount

  integer::igrid,ngrid,ncache,i,n,child,iskip,idim

  integer ,save,dimension(1:nvector)      ::ind_grid,ind_cell,ind_cell_father,ig,ih,nbor_grid,nbor_cell
  integer ,save,dimension(1:nvector,0:6)  ::nbor_father_cell,nbor_father_grid
  integer ,save,dimension(1:nvector,0:12) ::cell
  real(dp),save,dimension(1:nvector,0:12) ::p

  real(dp),save,dimension(1:nvector) :: &
         cx1, cy1, cz1, cx2, cy2, cz2, &
         nux1, nuy1, nuz1, nux2, nuy2, nuz2, &
         nux1_2, nuy1_2, nuz1_2, nux2_2, nuy2_2, nuz2_2, &
         gx1, gx2, gy1, gy2, gz1, gz2
  real(dp) :: dx    ! grid step size
  real(dp) :: h_i   ! inverse grid step size
  real(dp) :: factor_a, factor_b  ! weights for the 5-p fda
  
  if(verbose)write(*,111) 'compute_pdm_density_at_levelmin at level',ilevel

  if(ndim .ne. 3) then
     write(*,*) "Error: the MOND module can only be used with ndim=3 !"
     return
  endif

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  h_i = 1.0D0/dx
  factor_a = 27.0d0/24.0d0/dx
  factor_b = 1.0d0/24.0d0/dx

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid

    if (cosmo) then !Block added by IB for cosmo runs.
        FOUR_PI_G = 1.5D0*omega_m*aexp
        a0_i = 1.0d0/(aexp**4*a0)
    else
        FOUR_PI_G  = 4.0d0*PI*1.0d0
        !a0_i = 1.0d0/a0, calculated earlier.
    end if

  ! Compute the Newtonian acceleration, because it is needed below
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call gradient_phi(ind_grid,ngrid,ilevel,icount)
  enddo


  ! Update virtual boundaries
  do idim=1,ndim
     call make_virtual_fine_dp(f_newton(1,idim),ilevel)
  end do
  
  ! Compute the Newtonian grad phi at the physical boundary regions
  call connect_Newton
  call make_boundary_force(ilevel)
  call connect_Mond
  

  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring cells (later, find their grids via "igrid = son(igrid)")
     ! Every grid with index ind_grid(i) has 26 direct neighbours
     ! The neighboring cells can be addressed using the ggg and hhh coefficients (if they exist).
     do i=1,ngrid
        nbor_father_cell(i,0)=father(ind_grid(i)) ! The cell itself
        nbor_father_cell(i,ref_nbpg(1,1))=nbor(ind_grid(i),1) ! left
        nbor_father_cell(i,ref_nbpg(2,1))=nbor(ind_grid(i),2) ! right
        nbor_father_cell(i,ref_nbpg(1,2))=nbor(ind_grid(i),3) ! unten
        nbor_father_cell(i,ref_nbpg(2,2))=nbor(ind_grid(i),4) ! oben
        nbor_father_cell(i,ref_nbpg(1,3))=nbor(ind_grid(i),5) ! vorne
        nbor_father_cell(i,ref_nbpg(2,3))=nbor(ind_grid(i),6) ! hinten
     end do ! do i=1,ngrid

     ! Gather neighboring cells (later, find their grids via "igrid = son(igrid)")
     ! Every grid with index ind_grid(i) has 6 direct neighbours.
     do i=1,ngrid
        do n=0,6
           nbor_father_grid(i,n) = son(nbor_father_cell(i,n))
           if (nbor_father_grid(i,n) == 0) then
              write(*,*) 'Error in subroutine compute_pdm_density_at_levelmin: neighbor grid not available'
              call exit(1)
           endif
        enddo
     enddo

     ! Loop over grid cells:  ind = 1...8
     do child=1,8

        ! Find cell's indices
        iskip=ncoarse+(child-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        enddo

        ! Figure out the potential of the 26 direct neighbor cells (or interpolate)
        do n=0,12
           do i=1,ngrid
              ! ig: selects the correct neighboring grid at level ilevel
              ! ih: selects the correct child cell of the (with ig) selected neighboring grid
              ig(i) = ggg(n,child)  ! 0...6
              ih(i) = hhh(n,child)  ! 1...8
              nbor_grid(i) = nbor_father_grid(i,ig(i))
              nbor_cell(i) = ncoarse+(ih(i)-1)*ngridmax + nbor_father_grid(i,ig(i))
              ! The direct neighbor cell exists, take its value of the potential
              p(i,n) = phi_newton(nbor_cell(i))
              cell(i,n) = nbor_cell(i)
           enddo ! i=1,ngrid
        enddo ! n=1,12

        do i=1,ngrid
             !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
             !!                        Compute \rho_{PDM}                           !!
             !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

             ! Use a five-point FDA to approximate the gradient of phi at points (A/B)_(x/y/z)
             ! - grad phi(x) = - [27/24*(phi(x+h/2) - phi(x-h/2)) - 1/24*(phi(x+1.5h) - phi(x-1.5h))]/h
             ! factor_a/b include already the factor 1/h
             gx2(i) = -aexp**4*g_N_ext_x + ( factor_a*(p(i,ref_nbc(1,2,1)) - p(i,0)) - factor_b*(p(i,ref_nbc(2,2,1)) - p(i,ref_nbc(1,1,1))))
             gx1(i) = -aexp**4*g_N_ext_x + ( factor_a*(p(i,0) - p(i,ref_nbc(1,1,1))) - factor_b*(p(i,ref_nbc(1,2,1)) - p(i,ref_nbc(2,1,1))))
             gy2(i) = -aexp**4*g_N_ext_y + (factor_a*(p(i,ref_nbc(1,2,2)) - p(i,0)) - factor_b*(p(i,ref_nbc(2,2,2)) - p(i,ref_nbc(1,1,2))))
             gy1(i) = -aexp**4*g_N_ext_y + (factor_a*(p(i,0) - p(i,ref_nbc(1,1,2))) - factor_b*(p(i,ref_nbc(1,2,2)) - p(i,ref_nbc(2,1,2))))
             gz2(i) = -aexp**4*g_N_ext_z + (factor_a*(p(i,ref_nbc(1,2,3)) - p(i,0)) - factor_b*(p(i,ref_nbc(2,2,3)) - p(i,ref_nbc(1,1,3))))
             gz1(i) = -aexp**4*g_N_ext_z + (factor_a*(p(i,0) - p(i,ref_nbc(1,1,3))) - factor_b*(p(i,ref_nbc(1,2,3)) - p(i,ref_nbc(2,1,3))))

             !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

             ! B_x
             cy2(i) = -aexp**4*g_N_ext_y + 0.5d0 * ( f_newton(cell(i,0),2) + f_newton(cell(i,ref_nbc(1,2,1)),2) )
             cz2(i) = -aexp**4*g_N_ext_z + 0.5d0 * ( f_newton(cell(i,0),3) + f_newton(cell(i,ref_nbc(1,2,1)),3) )
             cx2(i) = gx2(i)

             ! A_x
             cy1(i) = -aexp**4*g_N_ext_y + 0.5d0 * ( f_newton(cell(i,0),2) + f_newton(cell(i,ref_nbc(1,1,1)),2) )
             cz1(i) = -aexp**4*g_N_ext_z + 0.5d0 * ( f_newton(cell(i,0),3) + f_newton(cell(i,ref_nbc(1,1,1)),3) )
             cx1(i) = gx1(i)

             ! grad(phi)/a0 at point B_x
             nux2(i) = sqrt(cx2(i)*cx2(i) + cy2(i)*cy2(i) + cz2(i)*cz2(i))*a0_i
             ! grad(phi)/a0 at point A_x
             nux1(i) = sqrt(cx1(i)*cx1(i) + cy1(i)*cy1(i) + cz1(i)*cz1(i))*a0_i

             !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

             ! B_y
             cy2(i) = gy2(i)
             cx2(i) = -aexp**4*g_N_ext_x + 0.5d0 * ( f_newton(cell(i,0),1) + f_newton(cell(i,ref_nbc(1,2,2)),1) )
             cz2(i) = -aexp**4*g_N_ext_z + 0.5d0 * ( f_newton(cell(i,0),3) + f_newton(cell(i,ref_nbc(1,2,2)),3) )
   
             ! A_y
             cy1(i) = gy1(i)
             cx1(i) = -aexp**4*g_N_ext_x + 0.5d0 * ( f_newton(cell(i,0),1) + f_newton(cell(i,ref_nbc(1,1,2)),1) )
             cz1(i) = -aexp**4*g_N_ext_z + 0.5d0 * ( f_newton(cell(i,0),3) + f_newton(cell(i,ref_nbc(1,1,2)),3) )

             ! grad(phi)/a0 at point B_y
             nuy2(i) = sqrt(cx2(i)*cx2(i) + cy2(i)*cy2(i) + cz2(i)*cz2(i))*a0_i
             ! grad(phi)/a0 at point A_y
             nuy1(i) = sqrt(cx1(i)*cx1(i) + cy1(i)*cy1(i) + cz1(i)*cz1(i))*a0_i

             !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
   
             ! B_z
             cz2(i) = gz2(i)
             cx2(i) = -aexp**4*g_N_ext_x + 0.5d0 * ( f_newton(cell(i,0),1) + f_newton(cell(i,ref_nbc(1,2,3)),1) )
             cy2(i) = -aexp**4*g_N_ext_y + 0.5d0 * ( f_newton(cell(i,0),2) + f_newton(cell(i,ref_nbc(1,2,3)),2) )
   
             ! A_z
             cz1(i) = gz1(i)
             cx1(i) = -aexp**4*g_N_ext_x + 0.5d0 * ( f_newton(cell(i,0),1) + f_newton(cell(i,ref_nbc(1,1,3)),1) )
             cy1(i) = -aexp**4*g_N_ext_y + 0.5d0 * ( f_newton(cell(i,0),2) + f_newton(cell(i,ref_nbc(1,1,3)),2) )

             ! grad(phi)/a0 at point B_z
             nuz2(i) = sqrt(cx2(i)*cx2(i) + cy2(i)*cy2(i) + cz2(i)*cz2(i))*a0_i
             ! grad(phi)/a0 at point A_z
             nuz1(i) = sqrt(cx1(i)*cx1(i) + cy1(i)*cy1(i) + cz1(i)*cz1(i))*a0_i

             !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

             ! Compute nu(x)
             call get_nu(nux1(i), nux1(i))  ! nu(x) at point A_x
             call get_nu(nuy1(i), nuy1(i))  ! nu(x) at point A_y
             call get_nu(nuz1(i), nuz1(i))  ! nu(x) at point A_z
             call get_nu(nux2(i), nux2(i))  ! nu(x) at point B_x
             call get_nu(nuy2(i), nuy2(i))  ! nu(x) at point B_y
             call get_nu(nuz2(i), nuz2(i))  ! nu(x) at point B_z

             ! Finally: rho_mond = rho + rho_ph
             if (cosmo) then
                 rho_mond(ind_cell(i)) = rho_newton(ind_cell(i)) &
                  + (  nux2(i)*gx2(i) - nux1(i)*gx1(i) + &
                       nuy2(i)*gy2(i) - nuy1(i)*gy1(i) + &
                       nuz2(i)*gz2(i) - nuz1(i)*gz1(i) ) *(h_i/boxlen)/FOUR_PI_G
             else
                 rho_mond(ind_cell(i)) = rho_Lambda_eff + rho_newton(ind_cell(i)) &
                  + (  nux2(i)*gx2(i) - nux1(i)*gx1(i) + &
                       nuy2(i)*gy2(i) - nuy1(i)*gy1(i) + &
                       nuz2(i)*gz2(i) - nuz1(i)*gz1(i) ) *(h_i/boxlen)/FOUR_PI_G
               end if

        enddo ! i=1,ngrid

     enddo ! child=1,8

  end do ! loop over myid grids by vector sweeps

  ! Update boundaries
!   call make_virtual_fine_dp(rho_mond(1),ilevel) This is actually not necessary
111 format('   Entering find_pdm_density_1 for level ',I2)

end subroutine compute_pdm_density_at_levelmin


!#########################################################
!##
!##   Computes the PDM density at level ilevel>levelmin
!##
!#########################################################
!##
!##     Uses a five-point finite difference approximation
!##     Interpolates the potential at the level boundaries
!##     Used at ilevel>levelmin
!##
!#########################################################
subroutine compute_pdm_density_at_fine_levels(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  use mond_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::ilevel,icount
  integer::igrid,ngrid,ncache,i,n,child
  integer::ix,iy,iz
  integer,dimension(1:twotondim)::iskip
  integer ,save,dimension(1:nvector)::ind_grid,ind_cell,ind_father_cell
  integer ,save,dimension(1:nvector,1:twotondim)::nbors_father_grids_tmp ! temporary working space for a subroutine
  integer ,save,dimension(1:nvector,1:threetondim)::nbors_father_cells,nbors_father_grids
  real(dp),save,dimension(1:nvector,1:twotondim,1:threetondim)::pp
  real(dp),save,dimension(1:twotondim,1:3)::xc
  real(dp)::c1,c2,c3,c4

  real(dp),save,dimension(1:nvector) :: &
         cx1, cy1, cz1, cx2, cy2, cz2, &
         nux1, nuy1, nuz1, nux2, nuy2, nuz2, &
         gx1, gx2, gy1, gy2, gz1, gz2
  real(dp) :: dx    ! grid step size
  real(dp) :: h_i   ! inverse grid step size
  real(dp) :: h4_i  ! the inverse of four times the grid step size
  real(dp) :: h24_i ! the inverse of four times the grid step size

  if(verbose)write(*,111) 'compute_pdm_density_at_fine_levels at level',ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  h_i = 1.0D0/dx
  h4_i = 1.0D0/(4.D0 * dx)
  h24_i = 1.0D0/(24.D0 * dx)
  c1 = 27.0d0/24.0d0/dx
  c2 =  1.0d0/24.0d0/dx
  c3 =  1.0d0/24.0d0/dx
  c4 =  8.0d0/24.0d0/dx

    if (cosmo) then !Block added by IB for cosmo runs.
        FOUR_PI_G = 1.5D0*omega_m*aexp
        a0_i = 1.0d0/(aexp**4*a0)
    else
        FOUR_PI_G  = 4.0d0*PI*1.0d0
        !a0_i = 1.0d0/a0, calculated earlier.
    end if

  ! Set position of cell centers relative to grid center
  do child=1,twotondim
      iz=(child-1)/4
      iy=(child-1-4*iz)/2
      ix=(child-1-2*iy-4*iz)
      if(ndim>0) xc(child,1)=(dble(ix)-0.5D0)*dx
      if(ndim>1) xc(child,2)=(dble(iy)-0.5D0)*dx
      if(ndim>2) xc(child,3)=(dble(iz)-0.5D0)*dx
      iskip(child) = ncoarse+(child-1)*ngridmax
  end do
  
  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  
  ! Loop over grids
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     
     ! Gather grids
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! For each grid, get its father cell
     do i=1,ngrid
        ind_father_cell(i)=father(ind_grid(i))
     end do
     
     ! For each grid, get all neighboring father cells
     call get3cubefather(ind_father_cell,nbors_father_cells,nbors_father_grids_tmp,ngrid,ilevel)
     
     ! For each grid, get all neighboring father grids (if exist, else 0)
     do n=1,threetondim
        do i=1,ngrid
           nbors_father_grids(i,n) = son(nbors_father_cells(i,n))
        end do
     end do

     ! For each grid, gather the potential of the all child-cells of the neighboring father grids
     do n=1,threetondim
        do i=1,ngrid
           if (nbors_father_grids(i,n) > 0) then
              ! If the neighboring father grid exists, copy the potential
              do child=1,twotondim
                 pp(i,child,n) = phi_newton(iskip(child)+nbors_father_grids(i,n))
              enddo
           else
              ! If the neighboring father grid does not exist,  
              ! interpolate the potential from the neighboring father cell
              call interpol_phi(nbors_father_cells(i,n), pp(i,1,n), 1, ilevel, icount)
           endif
        enddo
     end do

     ! Loop over cells
     do child=1,twotondim

        ! Get cell index
        do i=1,ngrid
           ind_cell(i)=iskip(child)+ind_grid(i)
        enddo

        do i=1,ngrid

           !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

           gx2(i) =  -aexp**4*g_N_ext_x + &
                    (c1*(pp(i,cell_ijk(child, 1, 0, 0),grid_ijk(child, 1, 0, 0)) - pp(i,cell_ijk(child, 0, 0, 0),grid_ijk(child, 0, 0, 0))) & 
                     - c2*(pp(i,cell_ijk(child, 2, 0, 0),grid_ijk(child, 2, 0, 0)) - pp(i,cell_ijk(child,-1, 0, 0),grid_ijk(child,-1, 0, 0))) )
           gy2(i) =  -aexp**4*g_N_ext_y + &
                    ( c1*(pp(i,cell_ijk(child, 0, 1, 0),grid_ijk(child, 0, 1, 0)) - pp(i,cell_ijk(child, 0, 0, 0),grid_ijk(child, 0, 0, 0))) & 
                     - c2*(pp(i,cell_ijk(child, 0, 2, 0),grid_ijk(child, 0, 2, 0)) - pp(i,cell_ijk(child, 0,-1, 0),grid_ijk(child, 0,-1, 0))) )
           gz2(i) =  -aexp**4*g_N_ext_z + &
                    ( c1*(pp(i,cell_ijk(child, 0, 0, 1),grid_ijk(child, 0, 0, 1)) - pp(i,cell_ijk(child, 0, 0, 0),grid_ijk(child, 0, 0, 0))) & 
                     - c2*(pp(i,cell_ijk(child, 0, 0, 2),grid_ijk(child, 0, 0, 2)) - pp(i,cell_ijk(child, 0, 0,-1),grid_ijk(child, 0, 0,-1))) )

           gx1(i) =  -aexp**4*g_N_ext_x + &
                    ( c1*(pp(i,cell_ijk(child, 0, 0, 0),grid_ijk(child, 0, 0, 0)) - pp(i,cell_ijk(child,-1, 0, 0),grid_ijk(child,-1, 0, 0))) & 
                     - c2*(pp(i,cell_ijk(child, 1, 0, 0),grid_ijk(child, 1, 0, 0)) - pp(i,cell_ijk(child,-2, 0, 0),grid_ijk(child,-2, 0, 0))) )
           gy1(i) =  -aexp**4*g_N_ext_y + &
                    ( c1*(pp(i,cell_ijk(child, 0, 0, 0),grid_ijk(child, 0, 0, 0)) - pp(i,cell_ijk(child, 0,-1, 0),grid_ijk(child, 0,-1, 0))) & 
                     - c2*(pp(i,cell_ijk(child, 0, 1, 0),grid_ijk(child, 0, 1, 0)) - pp(i,cell_ijk(child, 0,-2, 0),grid_ijk(child, 0,-2, 0))) )
           gz1(i) =  -aexp**4*g_N_ext_z + &
                    ( c1*(pp(i,cell_ijk(child, 0, 0, 0),grid_ijk(child, 0, 0, 0)) - pp(i,cell_ijk(child, 0, 0,-1),grid_ijk(child, 0, 0,-1))) & 
                     - c2*(pp(i,cell_ijk(child, 0, 0, 1),grid_ijk(child, 0, 0, 1)) - pp(i,cell_ijk(child, 0, 0,-2),grid_ijk(child, 0, 0,-2))) )
           
           !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

           ! Ax
           cx1(i) = gx1(i)
           cy1(i) = -aexp**4*g_N_ext_y + &
                    c3*(pp(i,cell_ijk(child, 0,-2, 0),grid_ijk(child, 0,-2, 0)) &
                    +pp(i,cell_ijk(child,-1,-2, 0),grid_ijk(child,-1,-2, 0)) &
                    -pp(i,cell_ijk(child, 0, 2, 0),grid_ijk(child, 0, 2, 0)) &
                    -pp(i,cell_ijk(child,-1, 2, 0),grid_ijk(child,-1, 2, 0))) &
                +c4*(pp(i,cell_ijk(child, 0, 1, 0),grid_ijk(child, 0, 1, 0)) &
                    +pp(i,cell_ijk(child,-1, 1, 0),grid_ijk(child,-1, 1, 0)) &
                    -pp(i,cell_ijk(child, 0,-1, 0),grid_ijk(child, 0,-1, 0)) &
                    -pp(i,cell_ijk(child,-1,-1, 0),grid_ijk(child,-1,-1, 0)))
           cz1(i) = -aexp**4*g_N_ext_z + &
                    c3*(pp(i,cell_ijk(child, 0, 0,-2),grid_ijk(child, 0, 0,-2)) &
                    +pp(i,cell_ijk(child,-1, 0,-2),grid_ijk(child,-1, 0,-2)) &
                    -pp(i,cell_ijk(child, 0, 0, 2),grid_ijk(child, 0, 0, 2)) &
                    -pp(i,cell_ijk(child,-1, 0, 2),grid_ijk(child,-1, 0, 2))) &
                +c4*(pp(i,cell_ijk(child, 0, 0, 1),grid_ijk(child, 0, 0, 1)) &
                    +pp(i,cell_ijk(child,-1, 0, 1),grid_ijk(child,-1, 0, 1)) &
                    -pp(i,cell_ijk(child, 0, 0,-1),grid_ijk(child, 0, 0,-1)) &
                    -pp(i,cell_ijk(child,-1, 0,-1),grid_ijk(child,-1, 0,-1)))
           
           ! Bx 
           cx2(i) = gx2(i)
           cy2(i) = -aexp**4*g_N_ext_y + &
                    c3*(pp(i,cell_ijk(child, 0,-2, 0),grid_ijk(child, 0,-2, 0)) &
                    +pp(i,cell_ijk(child, 1,-2, 0),grid_ijk(child, 1,-2, 0)) &
                    -pp(i,cell_ijk(child, 0, 2, 0),grid_ijk(child, 0, 2, 0)) &
                    -pp(i,cell_ijk(child, 1, 2, 0),grid_ijk(child, 1, 2, 0))) &
                +c4*(pp(i,cell_ijk(child, 0, 1, 0),grid_ijk(child, 0, 1, 0)) &
                    +pp(i,cell_ijk(child, 1, 1, 0),grid_ijk(child, 1, 1, 0)) &
                    -pp(i,cell_ijk(child, 0,-1, 0),grid_ijk(child, 0,-1, 0)) &
                    -pp(i,cell_ijk(child, 1,-1, 0),grid_ijk(child, 1,-1, 0)))
           cz2(i) = -aexp**4*g_N_ext_z + &
                    c3*(pp(i,cell_ijk(child, 0, 0,-2),grid_ijk(child, 0, 0,-2)) &
                    +pp(i,cell_ijk(child, 1, 0,-2),grid_ijk(child, 1, 0,-2)) &
                    -pp(i,cell_ijk(child, 0, 0, 2),grid_ijk(child, 0, 0, 2)) &
                    -pp(i,cell_ijk(child, 1, 0, 2),grid_ijk(child, 1, 0, 2))) &
                +c4*(pp(i,cell_ijk(child, 0, 0, 1),grid_ijk(child, 0, 0, 1)) &
                    +pp(i,cell_ijk(child, 1, 0, 1),grid_ijk(child, 1, 0, 1)) &
                    -pp(i,cell_ijk(child, 0, 0,-1),grid_ijk(child, 0, 0,-1)) &
                    -pp(i,cell_ijk(child, 1, 0,-1),grid_ijk(child, 1, 0,-1)))
                    
           ! grad(phi)/a0 at point A_x
           nux1(i) = sqrt(cx1(i)*cx1(i) + cy1(i)*cy1(i) + cz1(i)*cz1(i))*a0_i;
           ! grad(phi)/a0 at point B_x
           nux2(i) = sqrt(cx2(i)*cx2(i) + cy2(i)*cy2(i) + cz2(i)*cz2(i))*a0_i;
           
           !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
                      
           ! Ay 
           cx1(i) = -aexp**4*g_N_ext_x + &
                 c3*(pp(i,cell_ijk(child,-2, 0, 0),grid_ijk(child,-2, 0, 0)) &
                    +pp(i,cell_ijk(child,-2,-1, 0),grid_ijk(child,-2,-1, 0)) &
                    -pp(i,cell_ijk(child, 2, 0, 0),grid_ijk(child, 2, 0, 0)) &
                    -pp(i,cell_ijk(child, 2,-1, 0),grid_ijk(child, 2,-1, 0))) &
                +c4*(pp(i,cell_ijk(child, 1, 0, 0),grid_ijk(child, 1, 0, 0)) &
                    +pp(i,cell_ijk(child, 1,-1, 0),grid_ijk(child, 1,-1, 0)) &
                    -pp(i,cell_ijk(child,-1, 0, 0),grid_ijk(child,-1, 0, 0)) &
                    -pp(i,cell_ijk(child,-1,-1, 0),grid_ijk(child,-1,-1, 0)))
           cy1(i) = gy1(i)
           cz1(i) = -aexp**4*g_N_ext_z + &
                 c3*(pp(i,cell_ijk(child, 0, 0,-2),grid_ijk(child, 0, 0,-2)) &
                    +pp(i,cell_ijk(child, 0,-1,-2),grid_ijk(child, 0,-1,-2)) &
                    -pp(i,cell_ijk(child, 0, 0, 2),grid_ijk(child, 0, 0, 2)) &
                    -pp(i,cell_ijk(child, 0,-1, 2),grid_ijk(child, 0,-1, 2))) &
                +c4*(pp(i,cell_ijk(child, 0, 0, 1),grid_ijk(child, 0, 0, 1)) &
                    +pp(i,cell_ijk(child, 0,-1, 1),grid_ijk(child, 0,-1, 1)) &
                    -pp(i,cell_ijk(child, 0, 0,-1),grid_ijk(child, 0, 0,-1)) &
                    -pp(i,cell_ijk(child, 0,-1,-1),grid_ijk(child, 0,-1,-1)))
           
           ! By 
           cx2(i) = -aexp**4*g_N_ext_x + &
                 c3*(pp(i,cell_ijk(child,-2, 0, 0),grid_ijk(child,-2, 0, 0)) &
                    +pp(i,cell_ijk(child,-2, 1, 0),grid_ijk(child,-2, 1, 0)) &
                    -pp(i,cell_ijk(child, 2, 0, 0),grid_ijk(child, 2, 0, 0)) &
                    -pp(i,cell_ijk(child, 2, 1, 0),grid_ijk(child, 2, 1, 0))) &
                +c4*(pp(i,cell_ijk(child, 1, 0, 0),grid_ijk(child, 1, 0, 0)) &
                    +pp(i,cell_ijk(child, 1, 1, 0),grid_ijk(child, 1, 1, 0)) &
                    -pp(i,cell_ijk(child,-1, 0, 0),grid_ijk(child,-1, 0, 0)) &
                    -pp(i,cell_ijk(child,-1, 1, 0),grid_ijk(child,-1, 1, 0)))
           cy2(i) = gy2(i)
           cz2(i) = -aexp**4*g_N_ext_z + &
                 c3*(pp(i,cell_ijk(child, 0, 0,-2),grid_ijk(child, 0, 0,-2)) &
                    +pp(i,cell_ijk(child, 0, 1,-2),grid_ijk(child, 0, 1,-2)) &
                    -pp(i,cell_ijk(child, 0, 0, 2),grid_ijk(child, 0, 0, 2)) &
                    -pp(i,cell_ijk(child, 0, 1, 2),grid_ijk(child, 0, 1, 2))) &
                +c4*(pp(i,cell_ijk(child, 0, 0, 1),grid_ijk(child, 0, 0, 1)) &
                    +pp(i,cell_ijk(child, 0, 1, 1),grid_ijk(child, 0, 1, 1)) &
                    -pp(i,cell_ijk(child, 0, 0,-1),grid_ijk(child, 0, 0,-1)) &
                    -pp(i,cell_ijk(child, 0, 1,-1),grid_ijk(child, 0, 1,-1)))

           ! grad(phi)/a0 at point A_y
           nuy1(i) = sqrt(cx1(i)*cx1(i) + cy1(i)*cy1(i) + cz1(i)*cz1(i))*a0_i;
           ! grad(phi)/a0 at point B_y
           nuy2(i) = sqrt(cx2(i)*cx2(i) + cy2(i)*cy2(i) + cz2(i)*cz2(i))*a0_i;
                    
           !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!
           
           ! Az 
           cx1(i) = -aexp**4*g_N_ext_x + &
                 c3*(pp(i,cell_ijk(child,-2, 0, 0),grid_ijk(child,-2, 0, 0)) &
                    +pp(i,cell_ijk(child,-2, 0,-1),grid_ijk(child,-2, 0,-1)) &
                    -pp(i,cell_ijk(child, 2, 0, 0),grid_ijk(child, 2, 0, 0)) &
                    -pp(i,cell_ijk(child, 2, 0,-1),grid_ijk(child, 2, 0,-1))) &
                +c4*(pp(i,cell_ijk(child, 1, 0, 0),grid_ijk(child, 1, 0, 0)) &
                    +pp(i,cell_ijk(child, 1, 0,-1),grid_ijk(child, 1, 0,-1)) &
                    -pp(i,cell_ijk(child,-1, 0, 0),grid_ijk(child,-1, 0, 0)) &
                    -pp(i,cell_ijk(child,-1, 0,-1),grid_ijk(child,-1, 0,-1)))
           cy1(i) = -aexp**4*g_N_ext_y + &
                 c3*(pp(i,cell_ijk(child, 0,-2, 0),grid_ijk(child, 0,-2, 0)) &
                    +pp(i,cell_ijk(child, 0,-2,-1),grid_ijk(child, 0,-2,-1)) &
                    -pp(i,cell_ijk(child, 0, 2, 0),grid_ijk(child, 0, 2, 0)) &
                    -pp(i,cell_ijk(child, 0, 2,-1),grid_ijk(child, 0, 2,-1))) &
                +c4*(pp(i,cell_ijk(child, 0, 1, 0),grid_ijk(child, 0, 1, 0)) &
                    +pp(i,cell_ijk(child, 0, 1,-1),grid_ijk(child, 0, 1,-1)) &
                    -pp(i,cell_ijk(child, 0,-1, 0),grid_ijk(child, 0,-1, 0)) &
                    -pp(i,cell_ijk(child, 0,-1,-1),grid_ijk(child, 0,-1,-1)))
           cz1(i) = gz1(i)
           
           ! Bz 
           cx2(i) = -aexp**4*g_N_ext_x + &
                 c3*(pp(i,cell_ijk(child,-2, 0, 0),grid_ijk(child,-2, 0, 0)) &
                    +pp(i,cell_ijk(child,-2, 0, 1),grid_ijk(child,-2, 0, 1)) &
                    -pp(i,cell_ijk(child, 2, 0, 0),grid_ijk(child, 2, 0, 0)) &
                    -pp(i,cell_ijk(child, 2, 0, 1),grid_ijk(child, 2, 0, 1))) &
                +c4*(pp(i,cell_ijk(child, 1, 0, 0),grid_ijk(child, 1, 0, 0)) &
                    +pp(i,cell_ijk(child, 1, 0, 1),grid_ijk(child, 1, 0, 1)) &
                    -pp(i,cell_ijk(child,-1, 0, 0),grid_ijk(child,-1, 0, 0)) &
                    -pp(i,cell_ijk(child,-1, 0, 1),grid_ijk(child,-1, 0, 1)))
           cy2(i) = -aexp**4*g_N_ext_y + &
                 c3*(pp(i,cell_ijk(child, 0,-2, 0),grid_ijk(child, 0,-2, 0)) &
                    +pp(i,cell_ijk(child, 0,-2, 1),grid_ijk(child, 0,-2, 1)) &
                    -pp(i,cell_ijk(child, 0, 2, 0),grid_ijk(child, 0, 2, 0)) &
                    -pp(i,cell_ijk(child, 0, 2, 1),grid_ijk(child, 0, 2, 1))) &
                +c4*(pp(i,cell_ijk(child, 0, 1, 0),grid_ijk(child, 0, 1, 0)) &
                    +pp(i,cell_ijk(child, 0, 1, 1),grid_ijk(child, 0, 1, 1)) &
                    -pp(i,cell_ijk(child, 0,-1, 0),grid_ijk(child, 0,-1, 0)) &
                    -pp(i,cell_ijk(child, 0,-1, 1),grid_ijk(child, 0,-1, 1)))
           cz2(i) = gz2(i)

           ! grad(phi)/a0 at point A_z
           nuz1(i) = sqrt(cx1(i)*cx1(i) + cy1(i)*cy1(i) + cz1(i)*cz1(i))*a0_i;
           ! grad(phi)/a0 at point B_z
           nuz2(i) = sqrt(cx2(i)*cx2(i) + cy2(i)*cy2(i) + cz2(i)*cz2(i))*a0_i;
           
           !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ !!

           ! Compute nu(x)
           call get_nu(nux1(i), nux1(i))  ! nu(x) at point A_x
           call get_nu(nuy1(i), nuy1(i))  ! nu(x) at point A_y
           call get_nu(nuz1(i), nuz1(i))  ! nu(x) at point A_z
           call get_nu(nux2(i), nux2(i))  ! nu(x) at point B_x
           call get_nu(nuy2(i), nuy2(i))  ! nu(x) at point B_y
           call get_nu(nuz2(i), nuz2(i))  ! nu(x) at point B_z

           ! Finally: rho_mond = rho + rho_ph
           if (cosmo) then
               rho_mond(ind_cell(i)) = rho_newton(ind_cell(i)) + &
              ( nux2(i)*gx2(i)-nux1(i)*gx1(i) +       &
                nuy2(i)*gy2(i)-nuy1(i)*gy1(i) +       &
                nuz2(i)*gz2(i)-nuz1(i)*gz1(i)  ) * (h_i/boxlen) / FOUR_PI_G
           else
               rho_mond(ind_cell(i)) = rho_Lambda_eff + rho_newton(ind_cell(i)) + &
              ( nux2(i)*gx2(i)-nux1(i)*gx1(i) +       &
                nuy2(i)*gy2(i)-nuy1(i)*gy1(i) +       &
                nuz2(i)*gz2(i)-nuz1(i)*gz1(i)  ) * (h_i/boxlen) / FOUR_PI_G
            end if
        enddo ! i=1,ngrid
     enddo ! child=1,8
  end do ! loop over myid grids by vector sweeps

  ! Update boundaries
!   call make_virtual_fine_dp(rho_mond(1),ilevel) This is actually not necessary
111 format('   Entering find_pdm_density for level ',I2)
end subroutine compute_pdm_density_at_fine_levels
