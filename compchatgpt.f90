program compressible_ns_solver
    implicit none
    integer, parameter :: nx = 100, ny = 100
    integer, parameter :: nt = 1000
    integer, parameter :: i_max = nx+2, j_max = ny+2 ! ghost cells
    real, parameter :: gamma = 1.4
    real, parameter :: dx = 0.01, dy = 0.01, dt = 1e-4
    real, dimension(i_max, j_max) :: rho, u, v, p, e
    real, dimension(i_max, j_max) :: rho_new, u_new, v_new, p_new, e_new
    integer :: i, j, n

    call initialize(rho, u, v, p)

    do n = 1, nt
        call compute_energy(rho, u, v, p, e)
        call update(rho, u, v, p, e, rho_new, u_new, v_new, p_new, e_new)
        call apply_bc(rho_new, u_new, v_new, p_new)

        rho = rho_new
        u = u_new
        v = v_new
        p = p_new
    end do

    call write_output(rho, u, v, p)
end program compressible_ns_solver

subroutine initialize(rho, u, v, p)
    implicit none
    real, dimension(:,:), intent(out) :: rho, u, v, p
    integer :: i, j

    do j = 1, size(rho,2)
        do i = 1, size(rho,1)
            rho(i,j) = 1.0
            u(i,j) = 0.0
            v(i,j) = 0.0
            p(i,j) = 1.0
        end do
    end do
end subroutine initialize

subroutine compute_energy(rho, u, v, p, e)
    implicit none
    real, dimension(:,:), intent(in) :: rho, u, v, p
    real, dimension(:,:), intent(out) :: e
    integer :: i, j

    do j = 1, size(rho,2)
        do i = 1, size(rho,1)
            e(i,j) = p(i,j)/(gamma - 1.0) + 0.5*rho(i,j)*(u(i,j)**2 + v(i,j)**2)
        end do
    end do
end subroutine compute_energy

subroutine update(rho, u, v, p, e, rho_new, u_new, v_new, p_new, e_new)
    implicit none
    real, dimension(:,:), intent(in) :: rho, u, v, p, e
    real, dimension(:,:), intent(out) :: rho_new, u_new, v_new, p_new, e_new
    integer :: i, j

    ! This is a very simplified Euler update. Real solvers use fluxes (F, G) for each conserved quantity
    do j = 2, size(rho,2)-1
        do i = 2, size(rho,1)-1
            rho_new(i,j) = rho(i,j) - dt/dx * (rho(i+1,j)*u(i+1,j) - rho(i-1,j)*u(i-1,j)) / 2.0
            u_new(i,j)   = u(i,j) - dt/dx * (u(i+1,j)**2 - u(i-1,j)**2) / 2.0
            v_new(i,j)   = v(i,j) - dt/dy * (v(i,j+1)**2 - v(i,j-1)**2) / 2.0
            p_new(i,j)   = p(i,j) - dt/dx * (p(i+1,j)*u(i+1,j) - p(i-1,j)*u(i-1,j)) / 2.0
        end do
    end do
end subroutine update

subroutine apply_bc(rho, u, v, p)
    implicit none
    real, dimension(:,:), intent(inout) :: rho, u, v, p
    integer :: i, j

    ! Simple reflective boundaries
    do j = 1, size(rho,2)
        rho(1,j) = rho(2,j)
        u(1,j) = -u(2,j)
        v(1,j) = v(2,j)
        p(1,j) = p(2,j)

        rho(size(rho,1),j) = rho(size(rho,1)-1,j)
        u(size(rho,1),j) = -u(size(rho,1)-1,j)
        v(size(rho,1),j) = v(size(rho,1)-1,j)
        p(size(rho,1),j) = p(size(rho,1)-1,j)
    end do

    do i = 1, size(rho,1)
        rho(i,1) = rho(i,2)
        u(i,1) = u(i,2)
        v(i,1) = -v(i,2)
        p(i,1) = p(i,2)

        rho(i,size(rho,2)) = rho(i,size(rho,2)-1)
        u(i,size(rho,2)) = u(i,size(rho,2)-1)
        v(i,size(rho,2)) = -v(i,size(rho,2)-1)
        p(i,size(rho,2)) = p(i,size(rho,2)-1)
    end do
end subroutine apply_bc

subroutine write_output(rho, u, v, p)
    implicit none
    real, dimension(:,:), intent(in) :: rho, u, v, p
    integer :: i, j
    open(unit=10, file="output.dat", status="unknown")
    do j = 2, size(rho,2)-1
        do i = 2, size(rho,1)-1
            write(10,*) i, j, rho(i,j), u(i,j), v(i,j), p(i,j)
        end do
    end do
    close(10)
end subroutine write_output
/888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
program navier_stokes_fd
    implicit none
    ! Parameters
    integer, parameter :: nx = 100, ny = 100
    integer, parameter :: i_max = nx+2, j_max = ny+2
    integer, parameter :: nt = 1000
    real, parameter :: dx = 0.01, dy = 0.01, dt = 1e-5
    real, parameter :: gamma = 1.4
    real, parameter :: R = 287.0
    real, parameter :: mu = 1.8e-5
    real, parameter :: kappa = 0.025

    ! State variables
    real :: rho(i_max,j_max), u(i_max,j_max), v(i_max,j_max), T(i_max,j_max)
    real :: p(i_max,j_max), E(i_max,j_max)
    real :: rho_new(i_max,j_max), u_new(i_max,j_max), v_new(i_max,j_max), T_new(i_max,j_max)

    integer :: n

    ! Initialization
    call initialize(rho, u, v, T)
    call compute_pressure_energy(rho, u, v, T, p, E)

    ! Time-stepping
    do n = 1, nt
        call update(rho, u, v, T, p, E, rho_new, u_new, v_new, T_new)
        call apply_bc(rho_new, u_new, v_new, T_new)
        call compute_pressure_energy(rho_new, u_new, v_new, T_new, p, E)

        rho = rho_new
        u = u_new
        v = v_new
        T = T_new
    end do

    call write_output(rho, u, v, T, p)
end program navier_stokes_fd
!---------------------
subroutine initialize(rho, u, v, T)
    implicit none
    real, intent(out) :: rho(:,:), u(:,:), v(:,:), T(:,:)
    integer :: i, j

    do j = 1, size(rho,2)
        do i = 1, size(rho,1)
            rho(i,j) = 1.0
            u(i,j) = 0.0
            v(i,j) = 0.0
            T(i,j) = 300.0
        end do
    end do
end subroutine initialize

subroutine apply_bc(rho, u, v, T)
    implicit none
    real, intent(inout) :: rho(:,:), u(:,:), v(:,:), T(:,:)
    integer :: i, j

    do j = 1, size(rho,2)
        rho(1,j) = rho(2,j)
        u(1,j) = -u(2,j)
        v(1,j) = v(2,j)
        T(1,j) = T(2,j)

        rho(size(rho,1),j) = rho(size(rho,1)-1,j)
        u(size(rho,1),j) = -u(size(rho,1)-1,j)
        v(size(rho,1),j) = v(size(rho,1)-1,j)
        T(size(rho,1),j) = T(size(rho,1)-1,j)
    end do

    do i = 1, size(rho,1)
        rho(i,1) = rho(i,2)
        u(i,1) = u(i,2)
        v(i,1) = -v(i,2)
        T(i,1) = T(i,2)

        rho(i,size(rho,2)) = rho(i,size(rho,2)-1)
        u(i,size(rho,2)) = u(i,size(rho,2)-1)
        v(i,size(rho,2)) = -v(i,size(rho,2)-1)
        T(i,size(rho,2)) = T(i,size(rho,2)-1)
    end do
end subroutine apply_bc
!----------------------------
subroutine compute_pressure_energy(rho, u, v, T, p, E)
    implicit none
    real, intent(in) :: rho(:,:), u(:,:), v(:,:), T(:,:)
    real, intent(out) :: p(:,:), E(:,:)
    integer :: i, j

    do j = 1, size(rho,2)
        do i = 1, size(rho,1)
            p(i,j) = rho(i,j) * R * T(i,j)
            E(i,j) = p(i,j)/(gamma - 1.0) + 0.5 * rho(i,j) * (u(i,j)**2 + v(i,j)**2)
        end do
    end do
end subroutine compute_pressure_energy
!---------------------------
subroutine update(rho, u, v, T, p, E, rho_new, u_new, v_new, T_new)
    implicit none
    real, intent(in) :: rho(:,:), u(:,:), v(:,:), T(:,:), p(:,:), E(:,:)
    real, intent(out) :: rho_new(:,:), u_new(:,:), v_new(:,:), T_new(:,:)
    integer :: i, j
    real :: du_dx, du_dy, dv_dx, dv_dy, dT_dx, dT_dy
    real :: viscous_term_u, viscous_term_v, conduction_T
    real, parameter :: mu = 1.8e-5, kappa = 0.025

    do j = 2, size(rho,2)-1
        do i = 2, size(rho,1)-1
            ! Central difference for gradients
            du_dx = (u(i+1,j) - u(i-1,j)) / (2.0*dx)
            du_dy = (u(i,j+1) - u(i,j-1)) / (2.0*dy)
            dv_dx = (v(i+1,j) - v(i-1,j)) / (2.0*dx)
            dv_dy = (v(i,j+1) - v(i,j-1)) / (2.0*dy)
            dT_dx = (T(i+1,j) - T(i-1,j)) / (2.0*dx)
            dT_dy = (T(i,j+1) - T(i,j-1)) / (2.0*dy)

            ! Viscous stress contributions (simplified isotropic)
            viscous_term_u = mu * (du_dx + du_dy)
            viscous_term_v = mu * (dv_dx + dv_dy)

            ! Heat conduction
            conduction_T = kappa * (dT_dx + dT_dy)

            ! Update equations
            rho_new(i,j) = rho(i,j) - dt * (rho(i,j)*du_dx + rho(i,j)*dv_dy)
            u_new(i,j) = u(i,j) + dt * (-1.0/rho(i,j) * (p(i+1,j)-p(i-1,j))/(2.0*dx) + viscous_term_u/rho(i,j))
            v_new(i,j) = v(i,j) + dt * (-1.0/rho(i,j) * (p(i,j+1)-p(i,j-1))/(2.0*dy) + viscous_term_v/rho(i,j))
            T_new(i,j) = T(i,j) + dt * (conduction_T / (rho(i,j)*R))
        end do
    end do
end subroutine update
!-------------------------------------
subroutine write_output(rho, u, v, T, p)
    implicit none
    real, intent(in) :: rho(:,:), u(:,:), v(:,:), T(:,:), p(:,:)
    integer :: i, j
    open(10, file="solution.dat")
    do j = 2, size(rho,2)-1
        do i = 2, size(rho,1)-1
            write(10,*) i, j, rho(i,j), u(i,j), v(i,j), T(i,j), p(i,j)
        end do
    end do
    close(10)
end subroutine write_output
!-------------------------------------
subroutine compute_vorticity(u, v, omega, nx, ny, dx, dy)
    implicit none
    integer, intent(in) :: nx, ny
    real, intent(in) :: dx, dy
    real, intent(in) :: u(nx, ny), v(nx, ny)
    real, intent(out) :: omega(nx, ny)
    integer :: i, j
    omega=0.0
    ! Compute vorticity using central differences for interior points
    do j = 2, ny - 1
        do i = 2, nx - 1
            omega(i, j) = (v(i + 1, j) - v(i - 1, j)) / (2.0 * dx) - &
                          (u(i, j + 1) - u(i, j - 1)) / (2.0 * dy)
        end do
            i=1
	    omega(i, j) = (v(i + 1, j) - v(i - 0, j)) / (1.0 * dx) - &
                          (u(i, j + 1) - u(i, j - 1)) / (2.0 * dy)
            i=nx
            omega(i, j) = (v(i + 0, j) - v(i - 1, j)) / (1.0 * dx) - &
                          (u(i, j + 1) - u(i, j - 1)) / (2.0 * dy)
    end do
    do i = 2, nx-1
        j=1
        omega(i, j) = (v(i + 1, j) - v(i - 1, j)) / (2.0 * dx) - &
                          (u(i, j + 1) - u(i, j - 0)) / (1.0 * dy)
        j=ny
        omega(i, j) = (v(i + 1, j) - v(i - 1, j)) / (2.0 * dx) - &
                          (u(i, j + 0) - u(i, j - 1)) / (1.0 * dy)
    end do

    ! Optional: Set boundary vorticity to zero or apply appropriate boundary conditions
    do i = 1, nx
        !omega(i, 1) = 0.0
        !omega(i, ny) = 0.0
    end do

    do j = 1, ny
        !omega(1, j) = 0.0
        !omega(nx, j) = 0.0
    end do
end subroutine compute_vorticity
!--------------------------
subroutine compute_stream_function(u, v, psi, nx, ny, dx, dy)
    implicit none
    integer, intent(in) :: nx, ny
    real, intent(in) :: dx, dy
    real, intent(in) :: u(nx, ny), v(nx, ny)
    real, intent(out) :: psi(nx, ny)
    integer :: i, j

    ! Initialize the stream function at the origin
    psi(1,1) = 0.0

    ! Integrate along the first column (j = 1) using u
    do i = 2, nx
        psi(i,1) = psi(i-1,1) + u(i-1,1) * dy
    end do

    ! Integrate along the rows using v
    do i = 1, nx
        do j = 2, ny
            psi(i,j) = psi(i,j-1) - v(i,j-1) * dx
        end do
    end do
end subroutine compute_stream_function
#########################################################
module compressible_ns_rhs_mod
  implicit none
  private
  public :: compute_rhs

contains

  subroutine compute_rhs(q, rhs, dx, dy)
    implicit none
    real, intent(in)  :: q(:, :, :)
    real, intent(out) :: rhs(:, :, :)
    real, intent(in)  :: dx, dy

    integer, parameter :: nvars = 4
    integer :: i, j, n
    integer :: nx, ny
    real, parameter :: gamma = 1.4, R = 287.05
    real, parameter :: mu = 1e-3, kappa = 1e-2

    real, allocatable :: rho(:,:), u(:,:), v(:,:), p(:,:), T(:,:)
    real, allocatable :: tau_xx(:,:), tau_yy(:,:), tau_xy(:,:)
    real, allocatable :: qx(:,:), qy(:,:)
    real :: fx(nvars), fy(nvars)

    nx = size(q, 1)
    ny = size(q, 2)

    allocate(rho(nx,ny), u(nx,ny), v(nx,ny), p(nx,ny), T(nx,ny))
    allocate(tau_xx(nx,ny), tau_yy(nx,ny), tau_xy(nx,ny))
    allocate(qx(nx,ny), qy(nx,ny))

    call apply_boundary_conditions(q)
    call primitive_from_conserved(q, rho, u, v, p, T)
    call compute_viscous_flux(u, v, T, mu, kappa, dx, dy, tau_xx, tau_yy, tau_xy, qx, qy)

    rhs = 0.0

    do j = 2, ny-1
      do i = 2, nx-1
        fx(1) = (q(i+1,j,2) - q(i-1,j,2)) / (2*dx)
        fx(2) = ((q(i+1,j,2)**2/q(i+1,j,1) + p(i+1,j)) - (q(i-1,j,2)**2/q(i-1,j,1) + p(i-1,j))) / (2*dx) - &
                (tau_xx(i+1,j) - tau_xx(i-1,j)) / (2*dx)
        fx(3) = ((q(i+1,j,2)*q(i+1,j,3)/q(i+1,j,1)) - (q(i-1,j,2)*q(i-1,j,3)/q(i-1,j,1))) / (2*dx) - &
                (tau_xy(i+1,j) - tau_xy(i-1,j)) / (2*dx)
        fx(4) = (((q(i+1,j,4) + p(i+1,j)) * u(i+1,j)) - ((q(i-1,j,4) + p(i-1,j)) * u(i-1,j))) / (2*dx) - &
                ((u(i+1,j)*tau_xx(i+1,j) + v(i+1,j)*tau_xy(i+1,j) + qx(i+1,j)) - &
                 (u(i-1,j)*tau_xx(i-1,j) + v(i-1,j)*tau_xy(i-1,j) + qx(i-1,j))) / (2*dx)

        fy(1) = (q(i,j+1,3) - q(i,j-1,3)) / (2*dy)
        fy(2) = ((q(i,j+1,2)*q(i,j+1,3)/q(i,j+1,1)) - (q(i,j-1,2)*q(i,j-1,3)/q(i,j-1,1))) / (2*dy) - &
                (tau_xy(i,j+1) - tau_xy(i,j-1)) / (2*dy)
        fy(3) = ((q(i,j+1,3)**2/q(i,j+1,1) + p(i,j+1)) - (q(i,j-1,3)**2/q(i,j-1,1) + p(i,j-1))) / (2*dy) - &
                (tau_yy(i,j+1) - tau_yy(i,j-1)) / (2*dy)
        fy(4) = (((q(i,j+1,4) + p(i,j+1)) * v(i,j+1)) - ((q(i,j-1,4) + p(i,j-1)) * v(i,j-1))) / (2*dy) - &
                ((u(i,j+1)*tau_xy(i,j+1) + v(i,j+1)*tau_yy(i,j+1) + qy(i,j+1)) - &
                 (u(i,j-1)*tau_xy(i,j-1) + v(i,j-1)*tau_yy(i,j-1) + qy(i,j-1))) / (2*dy)

        do n = 1, nvars
          rhs(i,j,n) = - (fx(n) + fy(n))
        end do
      end do
    end do

    deallocate(rho, u, v, p, T, tau_xx, tau_yy, tau_xy, qx, qy)
  end subroutine compute_rhs

  subroutine primitive_from_conserved(q, rho, u, v, p, T)
    implicit none
    real, intent(in) :: q(:,:,:)
    real, intent(out) :: rho(:,:), u(:,:), v(:,:), p(:,:), T(:,:)
    integer :: i, j
    integer :: nx, ny
    real, parameter :: gamma = 1.4, R = 287.05
    real :: ke

    nx = size(q, 1)
    ny = size(q, 2)

    do j = 1, ny
      do i = 1, nx
        rho(i,j) = q(i,j,1)
        u(i,j) = q(i,j,2) / rho(i,j)
        v(i,j) = q(i,j,3) / rho(i,j)
        ke = 0.5 * rho(i,j) * (u(i,j)**2 + v(i,j)**2)
        p(i,j) = (gamma - 1.0) * (q(i,j,4) - ke)
        T(i,j) = p(i,j) / (rho(i,j) * R)
      end do
    end do
  end subroutine primitive_from_conserved

  subroutine compute_viscous_flux(u, v, T, mu, kappa, dx, dy, tau_xx, tau_yy, tau_xy, qx, qy)
    implicit none
    real, intent(in) :: u(:,:), v(:,:), T(:,:)
    real, intent(in) :: mu, kappa, dx, dy
    real, intent(out) :: tau_xx(:,:), tau_yy(:,:), tau_xy(:,:), qx(:,:), qy(:,:)
    integer :: i, j
    integer :: nx, ny
    real :: dudx, dudy, dvdx, dvdy, dTdx, dTdy

    nx = size(u, 1)
    ny = size(u, 2)

    do j = 2, ny-1
      do i = 2, nx-1
        dudx = (u(i+1,j) - u(i-1,j)) / (2*dx)
        dudy = (u(i,j+1) - u(i,j-1)) / (2*dy)
        dvdx = (v(i+1,j) - v(i-1,j)) / (2*dx)
        dvdy = (v(i,j+1) - v(i,j-1)) / (2*dy)
        dTdx = (T(i+1,j) - T(i-1,j)) / (2*dx)
        dTdy = (T(i,j+1) - T(i,j-1)) / (2*dy)

        tau_xx(i,j) = mu * (2*dudx - 2.0/3.0*(dudx + dvdy))
        tau_yy(i,j) = mu * (2*dvdy - 2.0/3.0*(dudx + dvdy))
        tau_xy(i,j) = mu * (dudy + dvdx)

        qx(i,j) = -kappa * dTdx
        qy(i,j) = -kappa * dTdy
      end do
    end do
  end subroutine compute_viscous_flux

  subroutine apply_boundary_conditions(q)
    implicit none
    real, intent(inout) :: q(:,:,:)
    integer :: i, j
    integer :: nx, ny

    nx = size(q, 1)
    ny = size(q, 2)

    do j = 1, ny
      q(1,j,:)  = q(2,j,:);   q(1,j,2)  = -q(2,j,2)
      q(nx,j,:) = q(nx-1,j,:); q(nx,j,2) = -q(nx-1,j,2)
    end do

    do i = 1, nx
      q(i,1,:)  = q(i,2,:);   q(i,1,3)  = -q(i,2,3)
      q(i,ny,:) = q(i,ny-1,:); q(i,ny,3) = -q(i,ny-1,3)
    end do
  end subroutine apply_boundary_conditions

end module compressible_ns_rhs_mod
############################RK3
  subroutine rk3_step(u, dt, rhs_func)
    implicit none
    real, intent(in) :: dt
    real, dimension(:), intent(inout) :: u
    real, dimension(size(u)) :: u1, u2, dudt
    ! Stage 1
    call rhs_func(u, dudt)
    u1 = u + dt * dudt
    ! Stage 2
    call rhs_func(u1, dudt)
    u2 = 0.75 * u + 0.25 * (u1 + dt * dudt)
    ! Stage 3
    call rhs_func(u2, dudt)
    u  = (1.0/3.0) * u + (2.0/3.0) * (u2 + dt * dudt)
  end subroutine rk3_step
#######################





