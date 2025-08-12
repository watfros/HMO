program HMO

implicit none  
integer Natom,Nelec
open(1001,file='1.dat') 
!write (*,*) 'Input the number of the atoms'
read(1001,*) Natom , Nelec
write (*,*) 'This molecule has ' , Natom , ' atoms, and ', Nelec , ' electrons.'
close(1001)
call realHMO(Natom,Nelec)

end program

subroutine realHMO(Natom,Nelec)
implicit none   
integer Natom,Nelec,Nocc,i,j,k
real*8  Matrix(Natom,Natom),C(Natom,Natom),Energy(Natom),Pij(Natom,Natom),F(Natom),Norbtelec(Natom),rho(Natom)
open(1001,file='1.dat') 
read(1001,*) 
read(1001,*) Matrix(:,:)
close(1001)
C(:,:)=Matrix(:,:)
CALL DMR(Natom,C,Energy)
C=transpose(C)
write (*,*) 'The Eigen Equation has been solved'
Nocc=(Nelec+1)/2
write (*,*)  Nocc, ' occupied orbitals.'
do i=1,Natom
if (i<=Nocc) then
Norbtelec(i)=2.0
else
Norbtelec(i)=0.0
endif
end do

rho(:)=0.0
do i=1,Natom
do k=1,Natom
rho(i)=rho(i)+Norbtelec(k)*C(k,i)*C(k,i)
end do
end do

Pij(:,:)=0.0
do i=1,Natom
do j=1,Natom
do k=1,Natom
if ( abs (Matrix(i,j) - 0.0D0) > 1.0D-6) then
!write(*,*) i,j
Pij(i,j)=Pij(i,j)+Norbtelec(k)*C(k,i)*C(k,j)
endif
end do
end do
end do

F(:)=0.0
do i=1,Natom
do j=1,Natom
F(i)=F(i)+Pij(i,j)
end do
end do

F(:)=1.7320508-abs(F(:))

open(1002,file='2.txt') 
write (*,*) 'Orbital Energy'
write (*,"(100f12.6)") Energy(:)
write (1002,*) 'Orbital Energy'
write (1002,"(100f12.6)") Energy(:)
write (*,*) 'Orbital coefficient'
write (1002,*) 'Orbital coefficient'
do i=1,Natom
write (*,"(100f12.6)") C(i,:)
write (1002,"(100f12.6)") C(i,:)
end do

Pij=abs(Pij)
write (*,*) 'bond order'
write (1002,*) 'bond order'
do i=1,Natom
write (*,"(100f12.6)") Pij(i,:)
write (1002,"(100f12.6)") Pij(i,:)
end do

write (*,*) 'rho'
write (*,"(100f12.6)") rho(:)
write (1002,*) 'rho'
write (1002,"(100f12.6)") rho(:)

write (*,*) 'free valence'
write (*,"(100f12.6)") F(:)
write (1002,*) 'free valence'
write (1002,"(100f12.6)") F(:)
close(1002)  
pause
end





!对角化
!n:矩阵维数
!as: 双精度 输入时为对称矩阵 输出为特征向量矩阵 列向量
!en： 本征值 从小到大排序



	SUBROUTINE DMR(N,AS,EN)
	INTEGER N
	DOUBLE PRECISION AS(N,N),B(N),EN(N)
	  
	call tred2(as,n,en,b)
      call tqli(en,b,n,as)

        do i=1,n-1
          do j=i+1,n
	      if (en(i).gt.en(j)) then
	        tem=en(i)
		    en(i)=en(j)
		    en(j)=tem
		    do l=1,n
		      tem=as(l,i)
		      as(l,i)=as(l,j)
		      as(l,j)=tem
		    end do
	      end if
	    end do
      end do

      return

      end 

! !!!!!求三对角矩阵的特征值
!
! n  整型变量，输入参数，矩阵阶数 
! d  n个元素的一维实型数组，输入、输出参数，输入时存放对称三对角矩阵的对角元，
!    输出时存放矩阵的特征值
! e  n个元素的一维实型数组，输入参数，存放三对角矩阵的非对角元素，其中e(1)任意
! z  n*n个元素的二维实型数组，输入、输出参数，如果是求三对角矩阵的特征向量，则输入为
!    单位矩阵

      SUBROUTINE tqli(d,e,n,z)
      INTEGER n
      double precision d(n),e(n),z(n,n)
      INTEGER i,iter,k,l,m
      LOGICAL done
      double precision b,c,dd,f,g,p,r,s
      do i=2,n
        e(i-1)=e(i)
      end do
      e(n)=0.d0
      do l=1,n
        iter=0
          do
            done=0
      do m=l,n-1
      dd=dabs(d(m))+dabs(d(m+1))
      if (dabs(e(m))+dd==dd) exit
      end do
	if((dabs(e(m))+dd)/=dd)  m=n
      if(m/=l) then
      if(iter==30) pause 'too many iterations in tqli'
      iter=iter+1
      g=(d(l+1)-d(l))/(2.*e(l))
      r=dsqrt(g**2+1.)
      g=d(m)-d(l)+e(l)/(g+sign(r,g))
      s=1.d0
      c=1.d0
      p=0.d0
      do i=m-1,l,-1
        f=s*e(i)
        b=c*e(i)
        if(dabs(f)>=dabs(g)) then
          c=g/f
		  r=dsqrt(c**2+1.)
		  e(i+1)=f*r
		  s=1/r
		  c=c*s
		else
          s=f/g
		  r=dsqrt(s**2+1.)
		  e(i+1)=g*r
		  c=1/r
		  s=c*s
        endif
        g=d(i+1)-p
        r=(d(i)-g)*s+2.*c*b
        p=s*r
        d(i+1)=g+p
        g=c*r-b
!Omit lines from here ...
        do k=1,n
          f=z(k,i+1)
          z(k,i+1)=s*z(k,i)+c*f
          z(k,i)=c*z(k,i)-s*f
        end do
!... to here when finding only eigenvalues.
      end do
      d(l)=d(l)-p
      e(l)=g
      e(m)=0.d0
      done=-1
      endif
	if (.not.done) exit
       end do
      end do
      END SUBROUTINE tqli

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       SUBROUTINE tred2(a,n,d,e)
         INTEGER n
        double precision a(n,n),d(n),e(n)
       INTEGER i,j,k,l
       double precision f,g,h,hh,scale
       do i=n,2,-1
      l=i-1
       h=0.d0
       scale=0.d0
       if(l>1)then
       do k=1,l
       scale=scale+dabs(a(i,k))
       end do
       if(scale==0.) then
      e(i)=a(i,l)
       else
      do k=1,l
        a(i,k)=a(i,k)/scale
        h=h+a(i,k)**2
      end do
      f=a(i,l)
      g=-sign(dsqrt(h),f)
      e(i)=scale*g
      h=h-f*g
      a(i,l)=f-g
      f=0.d0
      do j=1,l
!Omit following line if finding only eigenvalues
        a(j,i)=a(i,j)/h
        g=0.d0
        do k=1,j
          g=g+a(j,k)*a(i,k)
        end do
        do k=j+1,l
          g=g+a(k,j)*a(i,k)
        end do
        e(j)=g/h
        f=f+e(j)*a(i,j)
      end do
      hh=f/(h+h)
      do j=1,l
        f=a(i,j)
        g=e(j)-hh*f
        e(j)=g
        do k=1,j
          a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
        end do
      end do
       endif
      else
       e(i)=a(i,l)
       endif
       d(i)=h
	end do
!Omit following line if finding only eigenvalues.
	d(1)=0.d0
	e(1)=0.d0
	do i=1,n
!Delete lines from here ...
	l=i-1
	  if(d(i)/=0.)then
	    do j=1,l
	      g=0.d0
	      do k=1,l
        g=g+a(i,k)*a(k,j)
      end do
      do k=1,l
        a(k,j)=a(k,j)-g*a(k,i)
      end do
	    end do
	  endif
!... to here when finding only eigenvalues.
	d(i)=a(i,i)
!Also delete lines from here ...
	a(i,i)=1.d0
	do j=1,l
	a(i,j)=0.d0
	a(j,i)=0.d0
	end do
!... to here when finding only eigenvalues.
	end do
	END SUBROUTINE tred2 
