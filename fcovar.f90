! Cholesky decomposition of covariance matrix R, where
! R_{ij} = r(|j-i|) = \sum_{m=1}^k a(m) b(m)^|j-i| is real for i,j = 1,2,...
! and |b(m)|<1 for m = 1,...,k.
!
! The factorization is saved in cholsav, which must have size at least
! (2*k+2)*(n+1) real (8).
!
subroutine covchol(k,a,b,n,cholsav)
  implicit none
  integer :: k,n
  real (8) :: cholsav((2*k+2)*(n+1))
  complex (8) :: a(k),b(k)
  cholsav(1)=k
  cholsav(2)=n
  call covcholx(k,a,b,n,cholsav(3),cholsav(3+2*k), &
                cholsav(3+2*k*(1+n)),cholsav(2+(2*k+1)*(1+n)))
end subroutine covchol

subroutine covcholx(k,a,b,n,bsav,d,dsum,dsuminv)
  integer :: k,n, i,j,m
  real (8) :: dsum(n),dsuminv(n),s
  complex (8) :: a(k),b(k),bsav(k),d(k,n),bb(k,k),accum(k,k),v(k)
  bsav=b
  s=sum(a)
  d(:,1)=a/sqrt(s)
  dsum(1)=sum(d(:,1))
  dsuminv(1)=1/dsum(1)
  do i=1, k
     do j=1, k
        bb(j,i)=conjg(b(j))*b(i)
        accum(j,i)=conjg(d(j,1))*d(i,1)*bb(j,i)
     end do
  end do
  !
  ! Proceed row by row to factor
  !
  do m=2, n
     do i=1, k
        v(i)=a(i)-sum(accum(:,i))
     end do
     s=sum(v)
     d(:,m)=v/sqrt(s)
     do i=1, k
        do j=1, k
           accum(j,i)=(accum(j,i)+conjg(d(j,m))*d(i,m))*bb(j,i)
        end do
     end do
     dsum(m)=sum(d(:,m))
     dsuminv(m)=1/dsum(m)
  end do
end subroutine covcholx


! Matrix-vector product U' x = y, where R = U' U is the covariance matrix
! described above and cholsav was produced by covchol above.  The number
! of unknowns n can be less than or equal to the size for which covchol
! was created.
!
subroutine covprodut(n,x,y,cholsav)
  implicit none
  integer :: n,k, nsav
  real (8) :: x(n),y(n),cholsav(*)
  k=cholsav(1)
  nsav=cholsav(2)
  call covprodux(n,x,y,k,cholsav(3),cholsav(3+2*k), &
                 cholsav(3+2*k*(nsav+1)),cholsav(2+(2*k+1)*(nsav+1)))
end subroutine covprodut

subroutine covprodux(n,x,y,k,b,d,dsum,dsuminv)
  implicit none
  integer :: n,k, i
  real (8) :: x(n),y(n),dsum(n),dsuminv(n)
  complex (8) :: b(k),d(k,n),conjb(k),ss(k)
  conjb=conjg(b)
  ss=0
  do i=1, n
     y(i)=x(i)*dsum(i)+sum(ss)
     ss=(ss+x(i)*conjg(d(:,i)))*conjb
  end do
end subroutine covprodux


! Solve linear system U' x = y, where R = U' U is the covariance matrix
! described above and cholsav was produced by covchol above.  The number
! of unknowns n can be less than or equal to the size for which covchol
! was created.  The solution vector x can coincide in memory with the
! right hand side y.
!
subroutine covsolut(n,x,y,cholsav)
  implicit none
  integer :: n,k, nsav
  real (8) :: x(n),y(n),cholsav(*)
  k=cholsav(1)
  nsav=cholsav(2)
  call covsolux(n,x,y,k,cholsav(3),cholsav(3+2*k), &
                cholsav(3+2*k*(nsav+1)),cholsav(2+(2*k+1)*(nsav+1)))
end subroutine covsolut

subroutine covsolux(n,x,y,k,b,d,dsum,dsuminv)
  implicit none
  integer :: n,k, i
  real (8) :: x(n),y(n),dsum(n),dsuminv(n)
  complex (8) :: b(k),d(k,n),conjb(k),ss(k)
  conjb=conjg(b)
  ss=0
  do i=1, n
     x(i)=(y(i)-sum(ss))*dsuminv(i)
     ss=(ss+x(i)*conjg(d(:,i)))*conjb
  end do
end subroutine covsolux


! Error of y relative to x.
!
function relerr(n,x,y) result(r)
  implicit none
  integer :: n,i
  real (8) :: x(n),y(n),r, s1,s2
  s1=0
  s2=0
  do i=1, n
     s1=s1+(x(i)-y(i))**2
     s2=s2+x(i)**2
  end do
  r=sqrt(s1/s2)
end function relerr


! Test performance of Cholesky factorization and matrix-vector
! multiplications y = U' x and z = (U')^{-1} y.
!
program fcovar
  integer :: i
  integer, parameter :: k=11,n=1000000
  real (8) :: x(n),y(n),z(n),cholsav(30000000),relerr,t1,t2
  complex (8) :: a(k),b(k),ima=(0d0,1d0)
  a(1)=10.0+0.0*ima
  a(2)= 0.194487+0.405512*ima
  a(3)= 0.194487-0.405512*ima
  a(4)=-0.4358-0.0374477*ima
  a(5)=-0.4358+0.0374477*ima
  a(6)= 0.4986+0.31128*ima
  a(7)= 0.4986-0.31128*ima
  a(8)= 0.385488-0.00129318*ima
  a(9)= 0.385488+0.00129318*ima
  a(10)=-0.283494-0.291219*ima
  a(11)=-0.283494+0.291219*ima

  b(1)=0.12+0.0*ima
  b(2)=-0.320372+0.797491*ima
  b(3)=-0.320372-0.797491*ima
  b(4)=0.720776+0.102379*ima
  b(5)=0.720776-0.102379*ima
  b(6)=0.370054-0.0357288*ima
  b(7)=0.370054+0.0357288*ima
  b(8)=-0.652465+0.506429*ima
  b(9)=-0.652465-0.506429*ima
  b(10)=0.696761+0.622623*ima
  b(11)=0.696761-0.622623*ima

  call random_number(x)
  call cpu_time(t1)
  call covchol(k,a,b,n,cholsav)
  call cpu_time(t2)
  print *, 'Elapsed: ',t2-t1
  call cpu_time(t1)
  call covprodut(n,x,y,cholsav)
  call cpu_time(t2)
  print *, 'Elapsed: ',t2-t1
  call cpu_time(t1)
  call covsolut(n,z,y,cholsav)
  print *, 'relerr: ', relerr(n,x,z)
  call cpu_time(t2)
  print *, 'Elapsed: ',t2-t1
end program fcovar
