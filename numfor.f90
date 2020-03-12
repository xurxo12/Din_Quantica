subroutine Euler_B(psi, psi_old, n, steps, V, dt, m, dx, BC)
  integer :: n, steps, BC
  complex(8), dimension(n) :: psi, psi_old, psi_new
  complex(8) :: ii, constant
  real(8), dimension(n) :: V
  real(8) :: dt, m, dx
  ii = cmplx(0,1)
  constant = ii*2.d0*dt/(2*m*dx**2)

  do i=1,steps

    if (BC==0) then
      psi_new(1)=psi_old(1)+constant*(psi(n)+psi(2)-2*psi(1))-ii*dt*V(1)*psi(1)
      psi_new(n)=psi_old(n)+constant*(psi(n-1)+psi(1)-2*psi(n))-ii*dt*V(n)*psi(n)
    else if (BC==1) then
      psi_new(1)=cmplx(0)
      psi_new(n)=cmplx(0)
    else if (BC==2) then
      psi_new(1)=psi_new(2)
      psi_new(n)=psi_new(n-1)
    end if

    psi_new(2:n-1) = psi_old(2:n-1)+constant*(psi(:n-2)+psi(3:)-2*psi(2:n-1))-ii*2.d0*dt*V(2:n-1)*psi(2:n-1)

    psi_old = psi
    psi = psi_new

  end do

end subroutine

subroutine Euler_RK4(psi, n, steps, V, dt, m, dx)
  implicit none
  integer :: n, steps, i
  complex(8), dimension(n) :: psi, Psi_temp
  complex(8), dimension(n-2) :: k1, k2, k3, k4
  complex(8) :: ii
  real(8), dimension(n) :: V
  real(8) :: dt, m, dx
  ii = cmplx(0,1)

  Psi_temp(1) = cmplx(0)
  Psi_temp(n) = cmplx(0)
  do i=1,steps

    k1 = ii*dt*((Psi(:n-2)+Psi(3:)-2.d0*Psi(2:n-1))/(2.d0*m*dx**2)-V(2:n-1)*Psi(2:n-1))

    Psi_temp(2:n-1) = Psi(2:n-1) + 0.5d0*k1

    k2 = ii*dt*((Psi_temp(:n-2)+Psi_temp(3:)-2.d0*Psi_temp(2:n-1))/(2.d0*m*dx**2)-V(2:n-1)*Psi_temp(2:n-1))

    Psi_temp(2:n-1) = Psi(2:n-1) + 0.5d0*k2

    k3 = ii*dt*((Psi_temp(:n-2)+Psi_temp(3:)-2.d0*Psi_temp(2:n-1))/(2.d0*m*dx**2)-V(2:n-1)*Psi_temp(2:n-1))

    Psi_temp(2:n-1) = Psi(2:n-1) + k3

    k4 = ii*dt*((Psi_temp(:n-2)+Psi_temp(3:)-2.d0*Psi_temp(2:n-1))/(2.d0*m*dx**2)-V(2:n-1)*Psi_temp(2:n-1))

    psi(2:n-1) = psi(2:n-1)+(k1+2.d0*k2+2.d0*k3+k4)/6.d0
    ! Psi(1) = cmplx(0)
    ! Psi(n) = cmplx(0)
  end do

end subroutine

subroutine Euler_RK4_3D(psi, n, steps, dt, m, dx)
  implicit none
  integer :: n, steps, i
  complex(8), dimension(n, n) :: psi , Psi_temp
  complex(8), dimension(n-2, n-2) :: k1 , k2, k3, k4
  complex(8) :: ii
  ! real(8), dimension(n, v) :: V
  !f2py intent(inout) :: psi
  real(8) :: dt, m, dx
  ii = cmplx(0,1)

  Psi_temp(1,:) = cmplx(0)
  Psi_temp(n,:) = cmplx(0)
  Psi_temp(:,1) = cmplx(0)
  Psi_temp(:,n) = cmplx(0)

  do i=1,steps

    k1 = ii*dt*(Psi(2:n-1,:n-2)+Psi(2:n-1,3:)-2.d0*Psi(2:n-1,2:n-1) +&
                Psi(:n-2,2:n-1)+Psi(3:,2:n-1)-2.d0*Psi(2:n-1,2:n-1))/(2.d0*m*dx**2)

    Psi_temp(2:n-1, 2:n-1) = Psi(2:n-1, 2:n-1) + 0.5d0*k1

    k2 = ii*dt*(Psi_temp(2:n-1,:n-2)+Psi_temp(2:n-1,3:)-2.d0*Psi_temp(2:n-1,2:n-1)+&
                Psi_temp(:n-2,2:n-1)+Psi_temp(3:,2:n-1)-2.d0*Psi_temp(2:n-1,2:n-1))&
                /(2.d0*m*dx**2)

    Psi_temp(2:n-1, 2:n-1) = Psi(2:n-1, 2:n-1) + 0.5d0*k2

    k3 = ii*dt*(Psi_temp(2:n-1,:n-2)+Psi_temp(2:n-1,3:)-2.d0*Psi_temp(2:n-1,2:n-1)+&
                Psi_temp(:n-2,2:n-1)+Psi_temp(3:,2:n-1)-2.d0*Psi_temp(2:n-1,2:n-1))&
                /(2.d0*m*dx**2)

    Psi_temp(2:n-1, 2:n-1) = Psi(2:n-1, 2:n-1) + k3

    k4 = ii*dt*(Psi_temp(2:n-1,:n-2)+Psi_temp(2:n-1,3:)-2.d0*Psi_temp(2:n-1,2:n-1)+&
                Psi_temp(:n-2,2:n-1)+Psi_temp(3:,2:n-1)-2.d0*Psi_temp(2:n-1,2:n-1))&
                /(2.d0*m*dx**2)

    psi(2:n-1, 2:n-1) = psi(2:n-1, 2:n-1)+(k1+2.d0*k2+2.d0*k3+k4)/6.d0
  end do
end subroutine

subroutine Euler(psi, n, steps, V, dt, m, dx)
  integer :: n, steps
  double complex, dimension(n) :: psi
  double complex :: ii
  real(8), dimension(n) :: V
  real(8) :: dt, m, dx
  ii = cmplx(0,1)

  do i=1,steps
    psi(2:n-1) = psi(2:n-1)+ii*dt/(2*m*dx**2)*(psi(:n-2)+psi(3:)-2*psi(2:n-1))-ii*dt*V(2:n-1)*psi(2:n-1)
  end do

end subroutine
