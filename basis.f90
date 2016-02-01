
implicit none
integer           :: i,j,k,l,m,n,ll
double precision  :: nu,x
double precision, parameter :: space=1.0d0

!added changes in the y and z axis to basis

n=2 !<- possibly want to parameterize this and ll (max l value)?
nu=1.75
x=0.d0
ll=1
open(1,file='basis.dat')
write(1,*)n**3*(ll+1)**2
do l=0,ll
  do m=-l,l
    do i=1,n
      do j=1,n
	do k=1,n
	  write(1,*)l,m
	  write(1,*)nu
	  write(1,*)(i-(1.0*n+1)/2)*space,(j-(1.0*n+1)/2)*space,(k-(1.0*n+1)/2)*space
        end do
      end do  
    end do
  end do
end do
close(10)

end
