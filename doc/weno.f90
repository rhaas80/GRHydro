program weno

  integer, parameter :: nx = 7
  integer, parameter :: ny = 8
  integer, parameter :: nz = 9

  real*8, parameter :: dx = 0.1
  real*8, parameter :: dy = 0.2
  real*8, parameter :: dz = 0.3

  real*8 :: data(nx,ny,nz)
  real*8 :: outdata(nx,ny,nz)
  real*8 :: tmpdata1(nx,ny,nz)
  real*8 :: tmpdata2(nx,ny,nz)

  integer :: i,j,k
  real*8 :: x,y,z

  real*8 :: aFun, fun

  do k = 1+2, nz-2
    z = 4d0 + k * dz
    do j = 1+2,ny-2
      y = 2d0 + j * dy
      do i = 1+2,nx-2
        x = 1d0 + i * dx
        data(i,j,k) = (aFun(x+dx/2)-aFun(x-dx/2)) * &
                    & (aFun(y+dy/2)-aFun(y-dy/2)) * &
                    & (aFun(z+dz/2)-aFun(z-dz/2)) / (dx*dy*dz)
      end do
    end do
  end do

  call apply(data, nx, ny, nz, 0, tmpdata1)
  call apply(tmpdata1, nx, ny, nz, 1, tmpdata2)
  call apply(tmpdata2, nx, ny, nz, 2, outdata)


  do k = 1+2, nz-2
    z = 4d0 + k * dz
    do j = 1+2,ny-2
      y = 2d0 + j * dy
      do i = 1+2,nx-2
        x = 1d0 + i * dx
        if (abs(outdata(i,j,k) - fun(x)*fun(y)*fun(z)) .gt. 1d-6) then
          write(*,*) "Difference at location ",x," is: ", &
                     outdata(i,j,k) - fun(x)*fun(y)*fun(z), " between ", &
                     outdata(i,j,k), " and ", fun(x)*fun(y)*fun(z)
        end if
      end do
    end do
  end do

end program weno

function fun(x)
  implicit none

  real*8 :: fun
  real*8, intent(in) :: x

  fun = x**3
end function

function aFun(x)
  implicit none

  real*8 :: aFun
  real*8, intent(in) :: x

  aFun = (1d0/4.d0)*x**4
end function

subroutine CCTK_ERROR(str)
  implicit none

  character(*) :: str

  write(*,*) str
  stop
end subroutine

!Apply function as in weno.c++
subroutine apply(data, nx, ny, nz, dirn, a_center_xyz)
  use, intrinsic :: ieee_arithmetic

  implicit none

  ! Input/output parameters
  integer, intent(in) :: nx, ny, nz, dirn
  real*8, dimension(nx, ny, nz), intent(in) :: data
  real*8, dimension(nx, ny, nz), intent(out) :: a_center_xyz

  ! Local variables
  integer :: i, j, k, p
  integer :: imin, imax, jmin, jmax, kmin, kmax
  integer, parameter :: stencil_width = 2
  integer :: i_offset(5), j_offset(5), k_offset(5)
  real*8 :: A(5)
  real*8 :: beta1, beta2, beta3
  real*8 :: wbarplus1, wbarplus2, wbarplus3
  real*8 :: iwbarplussum, wplus1, wplus2, wplus3
  real*8, parameter :: weno_eps = 1.0d-10

  ! Beta coefficients for smoothness indicators
  real*8, parameter :: beta_shu(3,6) = reshape([ &
    4.0d0/3.0d0,  4.0d0/3.0d0,  10.0d0/3.0d0, &
    -19.0d0/3.0d0, -13.0d0/3.0d0, -31.0d0/3.0d0, &
    25.0d0/3.0d0,  13.0d0/3.0d0,  25.0d0/3.0d0, &
    11.0d0/3.0d0,  5.0d0/3.0d0,  11.0d0/3.0d0, &
    -31.0d0/3.0d0, -13.0d0/3.0d0, -19.0d0/3.0d0, &
    10.0d0/3.0d0,  4.0d0/3.0d0,  4.0d0/3.0d0 &
  ], [3, 6])

  ! WENO coefficients for de-averaging
  real*8, parameter :: weno_coeffs_de_avg(3,5) = reshape([ &
    -1.0d0/24.0d0,  0.0d0,         0.0d0, &
     1.0d0/12.0d0, -1.0d0/24.0d0,  0.0d0, &
    23.0d0/24.0d0,  13.0d0/12.0d0, 23.0d0/24.0d0, &
     0.0d0,        -1.0d0/24.0d0,  1.0d0/12.0d0, &
     0.0d0,         0.0d0,        -1.0d0/24.0d0 &
  ], [3, 5])

  ! Initialize parameters
  a_center_xyz = ieee_value(a_center_xyz(1,1,1), ieee_signaling_nan)

  ! Check input validity
  if (nx < 2*stencil_width + 1 .or. &
      ny < 2*stencil_width + 1 .or. &
      nz < 2*stencil_width + 1) then
    call CCTK_ERROR("Grid dimensions too small for stencil")
  endif

  if (dirn < 0 .or. dirn > 2) then
    call CCTK_ERROR("Invalid direction")
  endif

  ! Set up stencil offsets based on direction
  ! must be called in order 0,1,2 to correctly compute interior data without computing any NaN
  select case (dirn)
    case (0)  ! x-direction
      i_offset = [-2, -1, 0, 1, 2]
      j_offset = [0, 0, 0, 0, 0]
      k_offset = [0, 0, 0, 0, 0]
      imin = 1 + stencil_width
      imax = nx - stencil_width
      jmin = 1
      jmax = ny
      kmin = 1
      kmax = nz
    case (1)  ! y-direction
      i_offset = [0, 0, 0, 0, 0]
      j_offset = [-2, -1, 0, 1, 2]
      k_offset = [0, 0, 0, 0, 0]
      imin = 1 + stencil_width
      imax = nx - stencil_width
      jmin = 1 + stencil_width
      jmax = ny - stencil_width
      kmin = 1
      kmax = nz
    case (2)  ! z-direction
      i_offset = [0, 0, 0, 0, 0]
      j_offset = [0, 0, 0, 0, 0]
      k_offset = [-2, -1, 0, 1, 2]
      imin = 1 + stencil_width
      imax = nx - stencil_width
      jmin = 1 + stencil_width
      jmax = ny - stencil_width
      kmin = 1 + stencil_width
      kmax = nz - stencil_width
  end select

  ! Main computation loops
  do k = kmin, kmax
    do j = jmin, jmax
      do i = imin, imax
        ! Gather stencil values
        do p = 1, 5
          A(p) = data(i + i_offset(p), j + j_offset(p), k + k_offset(p))
        end do

        ! Calculate smoothness indicators
        beta1 = beta_shu(1,1)*A(1)*A(1) + &
                beta_shu(1,2)*A(1)*A(2) + &
                beta_shu(1,3)*A(2)*A(2) + &
                beta_shu(1,4)*A(1)*A(3) + &
                beta_shu(1,5)*A(2)*A(3) + &
                beta_shu(1,6)*A(3)*A(3)

        beta2 = beta_shu(2,1)*A(2)*A(2) + &
                beta_shu(2,2)*A(2)*A(3) + &
                beta_shu(2,3)*A(3)*A(3) + &
                beta_shu(2,4)*A(2)*A(4) + &
                beta_shu(2,5)*A(3)*A(4) + &
                beta_shu(2,6)*A(4)*A(4)

        beta3 = beta_shu(3,1)*A(3)*A(3) + &
                beta_shu(3,2)*A(3)*A(4) + &
                beta_shu(3,3)*A(4)*A(4) + &
                beta_shu(3,4)*A(3)*A(5) + &
                beta_shu(3,5)*A(4)*A(5) + &
                beta_shu(3,6)*A(5)*A(5)

        ! Calculate WENO weights
        wbarplus1 = -9.0d0/80.0d0 / ((weno_eps + beta1)*(weno_eps + beta1))
        wbarplus2 = 49.0d0/40.0d0 / ((weno_eps + beta2)*(weno_eps + beta2))
        wbarplus3 = -9.0d0/80.0d0 / ((weno_eps + beta3)*(weno_eps + beta3))

        iwbarplussum = 1.0d0 / (wbarplus1 + wbarplus2 + wbarplus3)

        wplus1 = wbarplus1 * iwbarplussum
        wplus2 = wbarplus2 * iwbarplussum
        wplus3 = wbarplus3 * iwbarplussum

        ! Calculate reconstruction
        a_center_xyz(i,j,k) = 0.0d0
        do p = 1, 5
          a_center_xyz(i,j,k) = a_center_xyz(i,j,k) + &
            (wplus1 * weno_coeffs_de_avg(1,p) + &
             wplus2 * weno_coeffs_de_avg(2,p) + &
             wplus3 * weno_coeffs_de_avg(3,p)) * A(p)
        end do

      end do
    end do
  end do

end subroutine apply
