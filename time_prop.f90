program time_dep_H_t

use input_mod
implicit none

  double precision, allocatable :: h(:,:),o(:,:),ht(:,:)
  complex*16, allocatable       :: hc(:,:),am(:,:),zo(:,:),bm(:,:),cm(:,:),c(:), xi(:),H_xi(:),c_gs(:)
  integer                       :: i,j,k,n,nt,m,l,ij,kl, p,q
  double precision              :: dt
  complex*16                    :: ef, prob, norm
  complex*16, parameter		:: zi = (0.,1.)
  integer, parameter 		:: N_taylor = 4
  complex*16                    :: tc(0:N_taylor)
  integer, parameter		:: eigenvaluesf=16,hamiltonianf=17, overlapf=19, ham_timef=23
  integer, parameter		:: c_groundstatef=29,wavef=19,densityf=50,normf=900
  character*255			:: hamn,overn,cn, htn
 
  call var_input

  if(ntsteps_inp==0) then
    call exit(1)
  end if
  
  
  open(hamiltonianf,file="mat_h.dat")
  open(overlapf,file="mat_o.dat")
  open(c_groundstatef,file="eigenvec.dat")
  open(ham_timef,file="mat_e.dat")
  open(55, file="groundstate.dat")
  open(57, file="ener.dat")
  
  read(57,*) hamn,n
  
  allocate(h(n,n),o(n,n),ht(n,n),c_gs(n))
  allocate(hc(n,n),zo(n,n),am(n,n),cm(n,n),bm(n,n),c(n))
  allocate(xi(n),H_xi(n))

  do i=1, n
    do j=1, n
      read(hamiltonianf,*) h(i,j)
      read(overlapf,*) o(i,j)
      read(ham_timef,*) ht(i,j)
    end do
  end do
  read(c_groundstatef,*) (c_gs(i),i=1,n)
  
  c=c_gs

  nt=ntsteps_inp
  dt=dt_inp
  
  ef=0.05d0  !!!!function for laser!!!!!!!!!!!!!!!!!!!!!!!!!

  
     

    do i=1,nt
      call laser_packet(ef,i*dt)
      hc=h+ef*ht
!   Crank-Nicholson matrix calculated
      am=o-0.5d0*zi*hc*dt
      zo=o+0.5d0*zi*hc*dt
  
      call inv_c(zo,n,bm)
      cm=matmul(bm,am)

      c=matmul(cm,c)  
    
      prob=0
      do q=1,n
        do p=1,n
          prob=prob+conjg(c(q))*c_gs(p)*o(q,p)
        end do
      end do
      write(55,*) real(conjg(prob)*prob)
      
      do q=1,n
        do p=1,n
          norm=norm+conjg(c(q))*c(p)*o(q,p)
        end do
      end do
      write(**,*) real(conjg(norm)*norm)
!  frequency of output
      if(mod(i,T_of_output)==0) then
        write(*,*) ef

      end if
    end do
    

  
  
end program time_dep_H_t

subroutine laser_packet(ef, t)
use input_mod
implicit none

double precision	:: ef, t, a, E_0, freq

  E_0 = sqrt(I_0)*5.329d-9
  freq = 2*3.14159*137*5.29177d-11/lambda

  ef =  E_0 * DSIN(freq*t) * exp(-2*(t-1.5*width)**2d0/(width**2))
  
end subroutine laser_packet

SUBROUTINE lubksb_c(a,n,np,indx,b)
  INTEGER	:: n,np,indx(n)
  complex*16 	:: a(np,np),b(n)
  INTEGER 	:: i,ii,j,ll
  complex*16 	:: sum
  
  ii=0
  do i=1,n
    ll=indx(i)
    sum=b(ll)
    b(ll)=b(i)
    if (ii.ne.0) then
    
      do j=ii,i-1
	sum=sum-a(i,j)*b(j)
      end do
      
    else if (sum.ne.(0.d0,0.d0)) then
      ii=i
    end if
    
    b(i)=sum
  end do
  
  do i=n,1,-1
    sum=b(i)
    do j=i+1,n
      sum=sum-a(i,j)*b(j)
    end do
    b(i)=sum/a(i,i)
  end do
  
  return
       
END SUBROUTINE lubksb_c



SUBROUTINE ludcmp_c(a,n,np,indx,d)
  INTEGER 		:: n,np,indx(n)
  REAL*8 		:: d
  COMPLEX*16		:: a(np,np),sum,du
  integer,PARAMETER	:: Ndim=5000
  Real*8,parameter	:: TINY=1.0d-20
  INTEGER		:: i,imax,j,k
  REAL*8		:: aamax,vv(ndim),dum
  
  d=1.d0
  do i=1,n
    aamax=0.d0
    do j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
    end do
      if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
    vv(i)=1.d0/aamax
  end do
  
  do j=1,n
    do i=1,j-1
      sum=a(i,j)
           
      do k=1,i-1
	sum=sum-a(i,k)*a(k,j)
      end do
           
      a(i,j)=sum
    end do
         
    aamax=0.d0
         
    do i=j,n
      sum=a(i,j)
           
      do k=1,j-1
	sum=sum-a(i,k)*a(k,j)
      end do
           
      a(i,j)=sum
      dum=vv(i)*abs(sum)
           
      if (dum.ge.aamax) then
	imax=i
        aamax=dum
      endif
    end do
 
    if (j.ne.imax)then
      do k=1,n
	du=a(imax,k)
        a(imax,k)=a(j,k)
        a(j,k)=du
      end do
      d=-d
      vv(imax)=vv(j)
    endif
         
    indx(j)=imax
    if(a(j,j).eq.(0.d0,0.d0)) a(j,j)=TINY
         
    if(j.ne.n) then
      du=1.d0/a(j,j)
      do i=j+1,n
	a(i,j)=a(i,j)*du
      end do
    endif
  end do
  
  return
END SUBROUTINE ludcmp_c

! 
subroutine inv_c(a,n,ai)
implicit none
  integer      :: n,i
  complex*16   :: a(n,n),ai(n,n)
  integer      :: indx(n)
  real*8       :: d

  ai=0.d0
  do i=1,n
    ai(i,i)=1.d0
  end do
  call ludcmp_c(a,n,n,indx,d)
  do i=1,n
    call lubksb_c(a,n,n,indx,ai(1,i))
  end do

end subroutine inv_c
