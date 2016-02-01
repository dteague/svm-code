
implicit none
integer           :: i,j,k,l,m,n,ll
double precision  :: nu,x
double precision, parameter :: space=0.75d0

!added changes in the y and z axis to basis

n=20 !<- possibly want to parameterize this and ll (max l value)?
nu=2.5d0
x=0.d0
ll=0
open(1,file='basis.dat')
write(1,*)n*(ll+1)**2
do l=0,ll
  do m=-l,l
    do i=1,n
 !     do j=1,n
!	do k=1,n
	  write(1,*)l,m
	  write(1,*)nu
	  write(1,*)(i-1)*space,x,x !<- it only gives l,m,nu values to points along the x axis
!        end do
!      end do  
    end do
  end do
end do
close(10)

end
