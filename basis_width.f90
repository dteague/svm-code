implicit none
integer           :: i,l,m,n,ll
double precision  :: nu,x,x0,a0


  open(37, file="basis.inp")

  read(37,*) n
  read(37,*) ll
  read(37,*) x
  read(37,*) x0
  read(37,*) a0



open(1,file='basis.dat')
write(1,*)n*(ll+1)**2
do l=0,ll
  do m=-l,l
    do i=1,n
      nu=1.d0/(a0*x0**(i-1))**2
      
      write(1,*)l,m
      write(1,*)nu
      write(1,*)x,x,x
    end do
  end do
end do
close(1)
close(37)

end
