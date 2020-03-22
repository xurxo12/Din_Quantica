subroutine Euler_B(psi, psi_old, n, N_steps, V, dt, m, dx)
  complex(8), dimension(n) :: psi, psi_old, psi_new
  real(8)                  :: V(n), dt, m, dx

  psi_new(1)=cmplx(0)
  psi_new(n)=cmplx(0)
  do i=1,N_steps

    psi_new(2:n-1) = psi_old(2:n-1)+cmplx(0,2.d0*dt)*((psi(:n-2)+psi(3:)-2*psi(2:n-1))/(2*m*dx**2)-V(2:n-1)*psi(2:n-1))

    psi_old = psi
    psi = psi_new

  end do
end subroutine

subroutine Euler_RK4(psi, n, N_steps, V, dt, m, dx)
  complex(8), dimension(n) :: psi, Psi_temp
  complex(8), dimension(n-2) :: k1, k2, k3, k4
  real(8) :: dt, m, dx, V(n)

  Psi_temp(1) = cmplx(0)
  Psi_temp(n) = cmplx(0)
  do i=1,N_steps

    k1 = cmplx(0,dt)*((Psi(:n-2)+Psi(3:)-2.*Psi(2:n-1))/(2.*m*dx**2)-V(2:n-1)*Psi(2:n-1))

    Psi_temp(2:n-1) = Psi(2:n-1) + 0.5*k1

    k2 = cmplx(0,dt)*((Psi_temp(:n-2)+Psi_temp(3:)-2.*Psi_temp(2:n-1))/(2.*m*dx**2)-V(2:n-1)*Psi_temp(2:n-1))

    Psi_temp(2:n-1) = Psi(2:n-1) + 0.5*k2

    k3 = cmplx(0,dt)*((Psi_temp(:n-2)+Psi_temp(3:)-2.*Psi_temp(2:n-1))/(2.*m*dx**2)-V(2:n-1)*Psi_temp(2:n-1))

    Psi_temp(2:n-1) = Psi(2:n-1) + k3

    k4 = cmplx(0,dt)*((Psi_temp(:n-2)+Psi_temp(3:)-2.*Psi_temp(2:n-1))/(2.d0*m*dx**2)-V(2:n-1)*Psi_temp(2:n-1))

    psi(2:n-1) = psi(2:n-1)+(k1+2.*k2+2.*k3+k4)/6.
  end do
end subroutine

subroutine Euler_RK4_3D(psi, n, steps, V, dt, m, dx)
  implicit none
  integer :: n, steps, i
  complex(8), dimension(n, n) :: psi , Psi_temp
  complex(8), dimension(n-2, n-2) :: k1 , k2, k3, k4
  complex(8) :: ii
  real(8), dimension(n, n) :: V
  real(8) :: dt, m, dx
  ii = cmplx(0,1)

  Psi_temp(1,:) = cmplx(0)
  Psi_temp(n,:) = cmplx(0)
  Psi_temp(:,1) = cmplx(0)
  Psi_temp(:,n) = cmplx(0)

  do i=1,steps

    k1 = ii*dt*((Psi(2:n-1,:n-2)+Psi(2:n-1,3:)-2.d0*Psi(2:n-1,2:n-1) + &
                Psi(:n-2,2:n-1)+Psi(3:,2:n-1)-2.d0*Psi(2:n-1,2:n-1))/(2.d0*m*dx**2) - &
                V(2:n-1,2:n-1)*Psi(2:n-1,2:n-1))

    Psi_temp(2:n-1, 2:n-1) = Psi(2:n-1, 2:n-1) + 0.5d0*k1

    k2 = ii*dt*((Psi_temp(2:n-1,:n-2)+Psi_temp(2:n-1,3:)-2.d0*Psi_temp(2:n-1,2:n-1)+&
                Psi_temp(:n-2,2:n-1)+Psi_temp(3:,2:n-1)-2.d0*Psi_temp(2:n-1,2:n-1))&
                /(2.d0*m*dx**2) - V(2:n-1,2:n-1)*Psi(2:n-1,2:n-1))

    Psi_temp(2:n-1, 2:n-1) = Psi(2:n-1, 2:n-1) + 0.5d0*k2

    k3 = ii*dt*((Psi_temp(2:n-1,:n-2)+Psi_temp(2:n-1,3:)-2.d0*Psi_temp(2:n-1,2:n-1)+&
                Psi_temp(:n-2,2:n-1)+Psi_temp(3:,2:n-1)-2.d0*Psi_temp(2:n-1,2:n-1))&
                /(2.d0*m*dx**2) - V(2:n-1,2:n-1)*Psi(2:n-1,2:n-1))

    Psi_temp(2:n-1, 2:n-1) = Psi(2:n-1, 2:n-1) + k3

    k4 = ii*dt*((Psi_temp(2:n-1,:n-2)+Psi_temp(2:n-1,3:)-2.d0*Psi_temp(2:n-1,2:n-1)+&
                Psi_temp(:n-2,2:n-1)+Psi_temp(3:,2:n-1)-2.d0*Psi_temp(2:n-1,2:n-1))&
                /(2.d0*m*dx**2) - V(2:n-1,2:n-1)*Psi(2:n-1,2:n-1))

    psi(2:n-1, 2:n-1) = psi(2:n-1, 2:n-1)+(k1+2.d0*k2+2.d0*k3+k4)/6.d0
  end do
end subroutine

subroutine Euler(psi, n, N_steps, V, dt, m, dx)
  complex(8) :: psi(n)
  real(8)    :: dt, m, dx, V(n)

  do i=1,N_steps
    psi(2:n-1) = psi(2:n-1)+cmplx(0,dt)*((psi(:n-2)+psi(3:)-2*psi(2:n-1))/(2*m*dx**2)-V(2:n-1)*psi(2:n-1))
  end do
end subroutine
