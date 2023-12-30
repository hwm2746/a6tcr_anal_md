module matvec
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function cross(a, b)
  double precision, dimension(3) :: cross
  double precision, dimension(3), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
end function cross

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function matdet3d(a)
  ! 3D matrix determinant
  double precision, dimension(3,3), intent(in) :: a
  double precision matdet3d
  integer i,j
  
  matdet3d= a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
           -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
           +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  !write (*,*) 'Determinant ', matdet3d
  
end function matdet3d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function matinv3d(a) ! result(ainv)
  ! 3D matrix inverse
  double precision, dimension(3,3), intent(in) :: a
  double precision, dimension(3,3):: matinv3d
  double precision rdum

  rdum=1./matdet3d(a)
  if (abs(rdum)<1.0e-32) then
     print *,'ERROR: Determinant too close to zero.'
     call exit(-1)
  endif

  
  matinv3d(1,1)=  (a(2,2)*a(3,3)-a(2,3)*a(3,2))*rdum
  matinv3d(2,1)= -(a(2,1)*a(3,3)-a(2,3)*a(3,1))*rdum
  matinv3d(3,1)=  (a(2,1)*a(3,2)-a(2,2)*a(3,1))*rdum
  matinv3d(1,2)= -(a(1,2)*a(3,3)-a(1,3)*a(3,2))*rdum
  matinv3d(2,2)=  (a(1,1)*a(3,3)-a(1,3)*a(3,1))*rdum
  matinv3d(3,2)= -(a(1,1)*a(3,2)-a(1,2)*a(3,1))*rdum
  matinv3d(1,3)=  (a(1,2)*a(2,3)-a(1,3)*a(2,2))*rdum
  matinv3d(2,3)= -(a(1,1)*a(2,3)-a(1,3)*a(2,1))*rdum
  matinv3d(3,3)=  (a(1,1)*a(2,2)-a(1,2)*a(2,1))*rdum

  return

end function matinv3d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine normalize(a)
  ! normalize vector a
  double precision, dimension(:), intent(inout):: a
  double precision rdum

  rdum=sqrt(dot_product(a,a))
  a=a/rdum
end subroutine normalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function rotate_perp(a,b,theta) result(c)
  ! For orthogonal vectors a and b, rotate b about a by angle theta.
  ! Let normalized versions of a,b,c as a0,b0,c0. Then
  !
  !   (a\times b)\cdot c=sin(theta)
  !   b\cdot c=cos(theta)
  !   a\cdot c=0
  !
  ! The matrix version of the above is inverted to get c.
  double precision, dimension(3), intent(in) :: a, b
  double precision, intent(in) :: theta
  double precision, dimension(3) :: c,a0,b0,c0
  double precision, dimension(3,3) :: m,minv
  double precision norm_b

  a0=a; call normalize(a0)
  norm_b=sqrt(dot_product(b,b));  b0=b/norm_b
  
  m(1,1)=a0(2)*b0(3)-a0(3)*b0(2) ! a0\times b0
  m(1,2)=a0(3)*b0(1)-a0(1)*b0(3)
  m(1,3)=a0(1)*b0(2)-a0(2)*b0(1)

  m(2,1)=b0(1);  m(2,2)=b0(2);  m(2,3)=b0(3);
  m(3,1)=a0(1);  m(3,2)=a0(2);  m(3,3)=a0(3);
  
  c(1)=sin(theta); c(2)=cos(theta); c(3)=0.
  
  minv=matinv3d(m)
  c0=matmul(minv,c); call normalize(c0)
  c=c0*norm_b ! make c and b have the same norm.
  return
  
end function rotate_perp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function rotate_vec(a,b,theta) result(c)
  ! For vectors a and b (not necessarily orthogonal), rotate b about a by
  ! angle theta.  Let normalized versions of a as a0. Then
  !
  !   ba = a0(a0\cdot b) : projection of b on a
  !   bp = b-ba: perp component of b about a
  !   c=rotate_perp(a,bp,theta)+ba
  !
  ! The matrix version of the above is inverted to get c.
  double precision, dimension(3), intent(in) :: a, b
  double precision, intent(in) :: theta
  double precision, dimension(3) :: c,a0,ba,bp
  double precision rdum
  
  a0=a; call normalize(a0)
  ba = dot_product(a0,b) * a0
  bp = b - ba
  
  c = rotate_perp(a,bp,theta) + ba
  return
  
end function rotate_vec


end module matvec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
