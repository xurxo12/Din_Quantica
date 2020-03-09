subroutine Euler_B(psi, psi_old, n, steps, V, dt, m, dx, BC)
  integer :: n, steps, BC
  double complex, dimension(n) :: psi, psi_old, psi_new
  double complex :: constant
  real(8), dimension(n) :: V
  real(8) :: dt, m, dx
  constant = cmplx(0,2)*dt/(2*m*dx**2)

  do i=1,steps

    if (BC==0) then
      psi_new(1)=psi_old(1)+constant*(psi(n)+psi(2)-2*psi(1))-cmplx(0,dt)*V(1)*psi(1)
      psi_new(n)=psi_old(n)+constant*(psi(n-1)+psi(1)-2*psi(n))-cmplx(0,dt)*V(n)*psi(n)
    else if (BC==1) then
      psi_new(1)=cmplx(0.d0,0.d0)
      psi_new(n)=cmplx(0.d0,0.d0)
    else if (BC==2) then
      psi_new(1)=psi_new(2)
      psi_new(n)=psi_new(n-1)
    end if

    psi_new(2:n-1) = psi_old(2:n-1)+constant*(psi(:n-2)+psi(3:)-2*psi(2:n-1))-cmplx(0,2.d0*dt)*V(2:n-1)*psi(2:n-1)

    psi_old = psi
    psi = psi_new

  end do

end subroutine

subroutine Euler_RK4(psi, n, steps, V, dt, m, dx)
  implicit none
  integer :: n, steps, i
  double complex, dimension(n) :: psi
  double complex, dimension(n-2) :: Psi_temp, k1, k2, k3, k4
  real(8), dimension(n) :: V
  real(8) :: dt, m, dx

  Psi_temp(1) = cmplx(0,0.d0)
  Psi_temp(n) = cmplx(0,0.d0)
  write(*,*) 'hello world!'
  write(*,*) dt
  write(*,*) V(1), 'avere'
  write(*,*) Psi_temp(1)
  do i=1,steps

    write(*,*) cmplx(3,dt)
    write(*,*) psi(n/4)

    k1 = cmplx(0,dt)*((psi(:n-2)+psi(3:)-2*psi(2:n-1))/(2*m*dx**2)-V(2:n-1)*psi(2:n-1))
    !
    ! Psi_temp(2:n-1) = Psi(2:n-1) + 0.5d0*k1
    !
    ! k2 = cmplx(0,dt)*((Psi_temp(:n-2)+Psi_temp(3:)-2*Psi_temp(2:n-1))/(2*m*dx**2)-V(2:n-1)*Psi_temp(2:n-1))
    !
    ! Psi_temp(2:n-1) = Psi(2:n-1) + 0.5d0*k2
    !
    ! k3 = cmplx(0,dt)*((Psi_temp(:n-2)+Psi_temp(3:)-2*Psi_temp(2:n-1))/(2*m*dx**2)-V(2:n-1)*Psi_temp(2:n-1))
    !
    ! Psi_temp(2:n-1) = Psi(2:n-1) + k3
    !
    ! k4 = cmplx(0,dt)*((Psi_temp(:n-2)+Psi_temp(3:)-2*Psi_temp(2:n-1))/(2*m*dx**2)-V(2:n-1)*Psi_temp(2:n-1))

    k2 = 0.d0
    k3 = 0.d0
    k4 = 0.d0

    psi(2:n-1) = psi(2:n-1)+(k1+2.d0*k2+3.d0*k3+k4)/6.d0

  end do

end subroutine


subroutine Euler(psi, n, steps, V, dt, m, dx)
  integer :: n, steps
  double complex, dimension(n) :: psi
  real(8), dimension(n) :: V
  real(8) :: dt, m, dx

  do i=1,steps
    psi(2:n-1) = psi(2:n-1)+cmplx(0,dt)/(2*m*dx**2)*(psi(:n-2)+psi(3:)-2*psi(2:n-1))-cmplx(0,dt)*V(2:n-1)*psi(2:n-1)
  end do

end subroutine
