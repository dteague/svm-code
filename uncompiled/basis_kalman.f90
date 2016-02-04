
implicit none
integer           :: i,l,m,n,ll
double precision  :: nu,x

n=15 !<- possibly want to parameterize this and ll (max l value)?
nu=2.5d0
x=0.d0
ll=1
open(1,file='basis.dat')
write(1,*)n*(ll+1)**2
do l=0,ll
  do m=-l,l
    do i=1,n
      write(1,*)l,m
      write(1,*)nu
      write(1,*)(i-1)*0.5d0,x,x !<- it only gives l,m,nu values to points along the x axis
    end do
  end do
end do
close(10)

end
