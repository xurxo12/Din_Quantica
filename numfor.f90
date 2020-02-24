subroutine Euler_B(psi, psi_old, n, steps, V, dt, m, dx)
  integer :: n, steps
  double complex, dimension(n) :: psi, psi_old, psi_new
  double complex :: constant
  real(8), dimension(n) :: V
  real(8) :: dt, m, dx
  constant = cmplx(0,1)*dt/(2*m*dx**2)
  ! constant2 = dt

  do i=1,steps

    psi_new(1)=0
    psi_new(n)=0
    psi_new(2:n-1) = psi_old(2:n-1)+constant*(psi(:n-2)+psi(3:)-2*psi(2:n-1))-cmplx(0,dt)*V(2:n-1)*psi(2:n-1)

    psi_old = psi
    psi = psi_new

  end do

end subroutine


subroutine Euler(psi, n, steps, V, dt, m, dx)
  integer :: n, steps
  double complex, dimension(n) :: psi
  real(8), dimension(n) :: V
  real(8) :: dt, m, dx

  do i=1,steps
    psi(2:n-1) = psi(2:n-1)+cmplx(0,1)*dt/(2*m*dx**2)*(psi(:n-2)+psi(3:)-2*psi(2:n-1))-cmplx(0,dt)*V(2:n-1)*psi(2:n-1)
  end do

end subroutine
