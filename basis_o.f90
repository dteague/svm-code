implicit none
integer           :: i,j,k,l,m,n,ll
double precision  :: nu,x
double precision  :: space

!added changes in the y and z axis to basis

open(37, file="basis.inp")

read(37,*) n
read(37,*) ll
read(37,*) x
read(37,*) nu
read(37,*) space


open(1,file='basis.dat')
write(1,*)n**3*(ll+1)**2
do l=0,ll
  do m=-l,l
    do i=-n/2,n/2
      do j=-n/2,n/2
	do k=-n/2,n/2
	  write(1,*)l,m
	  write(1,*)nu
	  write(1,*)space*i,space*j,space*k !<- it only gives l,m,nu values to points along the x axis
        end do
      end do  
    end do
  end do
end do
close(10)

end
