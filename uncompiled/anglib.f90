MODULE ANGLIB
  implicit none
  integer,parameter              :: nmax=100,nnmax=301,n2d=100,indm=10000
  real*8                         :: f(0:nmax),df(-1:nnmax),x2n(0:n2d)
  real*8                         :: dfa(0:2*nmax+1)
  real*8,parameter               :: pi=3.141592653589793d0
  double precision, allocatable  :: gaunt_coeff(:,:,:)
  double precision, allocatable  :: gaunt3_coeff(:,:,:,:,:)
  double precision, allocatable  :: gauntc3(:)
  double precision, allocatable  :: glm_c(:,:)



CONTAINS

! power(x,n)
! Description: finds x**n for real numbers
! parameters: 	x, n
! post: returns x**n taking into acount special cases of 0**n and 0**0 (=1)
function power(x,n)
implicit none
  real*8         :: x,power
  integer        :: n

  if(x.eq.0.d0) then
    if(n.eq.0) then
      power=1.d0
    else
      power=0.d0
    endif
  else
    power=x**n
  endif

end function power


subroutine init_angular
  implicit none
  real*8               :: x
  integer              :: i,j

  dfa(0)=1.d0
  do j=1,nmax
    dfa(j)=(2.d0*j+1.d0)*dfa(j-1)
  end do


  f(0)=0.d0
  df(-1)=0.d0
  df(0)=0.d0
  x2n(0)=0.d0
  f(1)=0.d0
  df(1)=0.d0
  df(2)=dlog(2.d0)
  do i=2,nmax
    x=dble(i)
    f(i)=f(i-1)+dlog(x)
  end do
  do i=3,nnmax
    x=dble(i)
    df(i)=df(i-2)+dlog(x)
  end do
  do i=1,n2d
    x=2.d0
    x2n(i)=x2n(i-1)+dlog(x)
  end do
  write(6,*)'init_angular   :  OK'
end subroutine  init_angular




function glm(l1,m1,l2,m2)
! J. Phys. B 30 2529, Eq. 15b
  implicit none
  double precision :: xf,glm,w
  integer          :: l1,m1,l2,m2

  glm=0.d0
  w=gaunt(l1,l2,m2,l1-l2,m1-m2)
  if(w/=0.d0) then
    xf=df(2*l1+1)-df(2*l2+1)-df(2*(l1-l2)+1)
    glm=w*dexp(xf)
  endif
  
end function glm






!###############################################################!
!			Spherical harmonics                     !
!###############################################################!

! Legendre_poly
! Description: the associated legendre polynomials calculated for l,m at x
! Parameters:	fact -
! 		pll
! 		pmm
! 		pmmpl
! 		somx2
FUNCTION Legendre_poly(l,m,x)
  INTEGER l,m
  REAL*8 Legendre_poly,x
  INTEGER i,ll
  REAL*8 fact,pll,pmm,pmmp1,somx2

  if((m.lt.0).or.(m.gt.l).or.(abs(x).gt.1.)) then         
    write(*,*)  'bad arguments in Legendre_poly'
    stop
  end if

  pmm=1.d0
  if(m.gt.0) then
    somx2=sqrt((1.d0-x)*(1.d0+x))
    fact=1.d0
    do i=1,m
      pmm=-pmm*fact*somx2
      fact=fact+2.d0
    end do
  endif
  if(l.eq.m) then
    Legendre_poly=pmm
  else
    pmmp1=x*(2*m+1)*pmm
    if(l.eq.m+1) then
      Legendre_poly=pmmp1
    else
      do ll=m+2,l
        pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
        pmm=pmmp1
        pmmp1=pll
      end do
      Legendre_poly=pll
    endif
  endif
  return
END function Legendre_poly



! powerc(x,n)
! Description: finds x**n for complex numbers
! parameters: 	x, n
! post: returns x**n taking into acount special cases of 0**n and 0**0 (=1)
function powerc(x,n)
  complex*16     :: x,powerc
  integer        :: n

  if(x.eq.(0.d0,0.d0)) then
    if(n.eq.0) then
      powerc=(1.d0,0.d0)
    else
      powerc=(0.d0,0.d0)
    endif
  else
    powerc=x**n
  endif

end function powerc


! spherical_harmonics
! Description: Calculates the spherical harmonic with l and m=mm
! Parameters:	cost - cos(\theta) but in cartesian coordinates
! 		phi - exp(i*\phi) expressed in cartesian coordinates
! 		zi - imaginary number i
! 		fc - constant for the spherical harmonic
! post: returns the spherical harmonic at r with l and m=mm as a complex number
function spherical_harmonics(r,l,mm)

  implicit none
  real*8                       :: x,y,z,cost,fc,r(3),r2
  complex*16,parameter         :: zi=(0.d0,1.d0)
  integer                      :: i,l,m,mm
  complex*16                   :: spherical_harmonics,phi
  spherical_harmonics=(0.d0,0.d0)

  m=abs(mm)
  if(abs(m).gt.l) then
    spherical_harmonics=(0.d0,0.d0)
    return
  endif

  x=r(1); y=r(2); z=r(3)
  r2=x**2+y**2+z**2

  fc=(2*l+1)/(4.d0*pi)
  do i=1,l-m
    fc=fc*i
  end do
  do i=1,l+m
    fc=fc/i
  end do

  phi=(x+zi*y)/sqrt(x**2+y**2+1.d-20) ! 1.d-20 added to prevent a singular answer
  cost=z/sqrt(r2+1.d-20) ! 1.d-20 added to prevent a singular answer which would crash the program

  spherical_harmonics=sqrt(fc)*powerc(phi,m)*Legendre_poly(l,m,cost)
  
  if(mm.lt.0) spherical_harmonics=(-1)**m*Conjg(spherical_harmonics)

end function spherical_harmonics


!########################################################!
!			Gaunt				 !
!########################################################!


function gaunt(l,l1,m1,l2,m2)
! J. Phys. B 30 2529, Eq. 16
!    gaunt(l,l2,m2,l1,m1)=<l,m|l2,m2|l1,m2> 
!
  implicit none
  real*8                       :: fc,gaunt,w
  integer                      :: l1,l2,l,m1,m2

  gaunt=0.d0
  w=clebgn(1.d0*l1,0.d0,1.d0*l2,0.d0,1.d0*l,0.d0)
  if(w/=0.d0) then
    fc=sqrt((2.d0*l1+1.d0)*(2.d0*l2+1.d0)/(4.d0*pi*(2.d0*l+1.d0)))
    gaunt=fc*clebgn(1.d0*l1,1.d0*m1,1.d0*l2,1.d0*m2,1.d0*l,1.d0*(m1+m2))*w
  endif


end function gaunt


function clebgn(a,az,b,bz,c,cz)
!
!        determines clebsch-gordan coefficients
!
  implicit none
  real*8              :: az,bz,cz,a,b,c,d(6),r,s,clebgn,p,t
  integer             :: k1,k2,k3,k4,k5,k6,i,j,k,mn,mx,m


  if(dabs(az+bz-cz)-0.1d0)1,4,4
1 r=a+b+c
  s=r-dint(r+1.1d0)+1.0d0
  if(s-0.1)9,9,4
9 if(dabs(az)-0.1)2,5,5
2 if(dabs(bz)-0.1)3,5,5
3 if(dsign(1.d0,2.d0*dint(.5d0*r+.2d0)+.5d0-r))4,5,5
4 clebgn=0.0d0
  return
5 d(1)=0.0d0
  d(2)=b-c-az
  d(3)=a-c+bz
  d(4)=a+b-c
  d(5)=a-az
  d(6)=b+bz
  r=dmax1(d(1),d(2),d(3))
  s=dmin1(d(4),d(5),d(6))
  if(s-r+0.1d0)4,6,6
6 mn=idint(r+0.1d0)
  mx=idint(s+0.1d0)
  k1=idint(a+az+1.1d0)
  k2=idint(d(5)+1.1d0)
  k3=idint(d(6)+1.1d0)
  k4=idint(b-bz+1.1d0)
  k5=idint(c+cz+1.1d0)
  k6=idint(c-cz+1.1d0)
  s=0.5*(f(k1-1)+f(k2-1)+f(k3-1)+f(k4-1)+f(k5-1)+f(k6-1))+delt(a,b,c)
  r=0.0d0
  p=dsign(1.d0,2.d0*dint(.5d0*dble(mn)+.2d0)+.5d0-dble(mn))
  do m=mn,mx
    t=-s
    do i=1,3
      j=idint(dble(m)-d(i)+1.1d0)
      k=idint(d(i+3)-dble(m)+1.1d0)
      t=t+f(j-1)+f(k-1)
    end do
    r=r+p*dexp(-t)
    p=-p
  end do
  clebgn=r*dsqrt(2.d0*c+1.d0)
  return
end function clebgn


function delt(a,b,c)
!
!         computes items used by subroutine clebgn
!
implicit none
  real*8               :: delt,a,b,c
  integer              :: i1,i2,i3,i4

  i1=idint(a+b+c+2.1d0)
  i2=idint(a+b-c+1.1d0)
  i3=idint(a+c-b+1.1d0)
  i4=idint(b+c-a+1.1d0)
  delt=0.5*(f(i2-1)+f(i3-1)+f(i4-1)-f(i1-1))
end function delt


subroutine precal_gaunt_3(lao)
implicit none
  double precision, allocatable :: g3(:)
  integer                       :: l1,m1,l2,m2,l3,m3,l4,m4,lm1,lm2,lm3,lm4,lm,lmmax,l,lao

  lmmax=4*lao
  lm=lao**2+(lao+1)+lao

  allocate(gaunt3_coeff(0:lmmax,0:lm,0:lm,0:lm,0:lm))
  allocate(gauntc3(0:lmmax))

  gaunt3_coeff=0.d0
  do l1=0,lao
    do m1=-l1,l1
      lm1=l1**2+(l1+1)+m1
      do l2=0,lao
	do m2=-l2,l2
	  lm2=l2**2+(l2+1)+m2
	  do l3=0,lao
	    do m3=-l3,l3
	      lm3=l3**2+(l3+1)+m3
	      do l4=0,lao
		do m4=-l4,l4
		  lm4=l4**2+(l4+1)+m4
		  call gaunt3(l2,m2,l4,m4,l1,m1,l3,m3)
		  do l=0,l1+l2+l3+l4
		    gaunt3_coeff(l,lm2,lm4,lm1,lm3)=gauntc3(l)
		  end do                
		end do
	      end do
	    end do
	  end do
	end do
      end do
    end do
  end do

end subroutine precal_gaunt_3



subroutine gaunt3(l2,m2,l4,m4,l1,m1,l3,m3)
! J. Phys. B 30 2529, Eq. 38
! J. Phys. B 30 2529, Eq. 16
!    gaunt(l,l2,m2,l1,m1)=<l,m|l2,m2|l1,m2>

  implicit none
  integer           :: l1,l2,l3,l4,m1,m2,m3,m4
  integer           :: l,l21,l43
  double precision  :: g21,g43

  gauntc3=0.d0
  do l21=iabs(l2-l1),l2+l1 
    g21=gaunt(l2,l21,m2-m1,l1,m1)
    do l43=iabs(l4-l3),l4+l3
      g43=gaunt(l4,l43,m4-m3,l3,m3)
      do l=iabs(l43-l21),l43+l21
        gauntc3(l)=gauntc3(l)+g21*g43*gaunt(l,l21,m2-m1,l43,m4-m3)
      end do
    end do
  end do
  
end subroutine gaunt3
  
function gauntc2(l2,m2,k,mu,l1,m1,l)
! J. Phys. B 30 2529, Eq. 25
!    gaunt(l,l2,m2,l1,m1)=<l,m|l2,m2|l1,m1>
!    gauntc2(l2,m2,k,mu,l1,m1,l) = <l2,m2,k,mu|l1,m1,l,m>
!    			 	 = sumofl'(<l2,m2|l1,m1|l',m'><l,m|k,mu|l',m'>
!				 = sumofl'(gaunt(l2,l1,m1,l',m')gaunt(l,k,mu,l',m')
! 
  implicit none
  integer           :: l1,l2,l,k,m1,m2,m,mu
  double precision  :: gauntc2
  integer           :: lp
  double precision  :: g1,g2

  gauntc2=0.d0
  do lp=iabs(l2-l1),l2+l1 
    g1=gaunt(l2,l1,m1,lp,m2-m1)
    g2=gaunt(l,k,mu,lp,m2-m1)
    gauntc2=gauntc2+g1*g2
  end do
  
end function gauntc2
  






!############################################################!
!			Unused code			     !
!############################################################!


function clb(l,lp,ll)
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  c     C(l,lp,L) coefficient (see (A2))
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
  real*8      :: clb
  integer     :: l,lp,ll
  real*8      :: xl,xlp,xll,xnul,hl,hlp,hll
  xl=dfloat(l)
  xlp=dfloat(lp)
  xll=dfloat(ll)
  xnul=0.d0
  hl=dsqrt(2.d0*l+1.d0)
  hlp=dsqrt(2.d0*lp+1.d0)
  hll=dsqrt(2.d0*ll+1.d0)
  clb=0.5d0*clebgn(xl,xnul,xlp,xnul,xll,xnul)*hl*hlp/(dsqrt(pi)*hll)

end function clb



function nine_j(j11,j12,j13,j21,j22,j23,j31,j32,j33)
! 9j symbol (arguments are 2*j)
!     | j11 j12 j13 |
!     | j21 j22 j23 |
!     | j31 j32 j33 |
  implicit none
  real*8    :: nine_j
  integer   :: j11,j12,j13,j21,j22,j23,j31,j32,j33
!
  integer   :: i1,i2,i3,jmin,jmax,ifas,j
  real*8    :: a,b,c

      nine_j=0.d0
      i1=j21-j32
      if(i1.lt.0)i1=-i1
      i2=j33-j11
      if(i2.lt.0)i2=-i2
      i3=j12-j23
      if(i3.lt.0)i3=-i3
      jmin=i1
      if(jmin.lt.i2)jmin=i2
      if(jmin.lt.i3)jmin=i3
      i1=j21+j32
      i2=j33+j11
      i3=j12+j23
      jmax=i1
      if(jmax.gt.i2)jmax=i2
      if(jmax.gt.i3)jmax=i3
      if(jmin.gt.jmax)return
      do j=jmin,jmax,2
        ifas=1
        if(j.ne.2*(j/2))ifas=-1
        a=six_j(j11,j21,j31,j32,j33,j)
        b=six_j(j12,j22,j32,j21,j,j23)
        c=six_j(j13,j23,j33,j,j11,j12)
        nine_j=nine_j+ifas*(dfloat(j)+1.d0)*a*b*c
      end do
end function nine_j


function delr(j1,j2,j3)
  implicit none
  real*8   :: delr
  integer  :: j1,j2,j3
  integer  :: jz1,jz2,jz3,jz4
      jz1=(j1+j2-j3)/2
      if(jz1.lt.0) go to 130
      jz2=(j1-j2+j3)/2
      if(jz2.lt.0) go to 130
      jz3=(j2+j3-j1)/2
      if(jz3.lt.0) go to 130
      jz4=(j1+j2+j3)/2+1
      if(jz4.gt.1000)go to 200
      delr=f(jz1)+f(jz2)+f(jz3)-f(jz4)
      return
130   delr=1.d+20
      return
200   write(6,*)'  too high spins'
      stop
end function delr


function six_j(j1,j2,j5,j4,j3,j6)
! 6j symbol (arguments are 2*j)
!     | j1 j2 j5 |
!     | j4 j3 j6 |
  implicit none
  real*8      :: six_j
  integer     :: j1,j2,j3,j4,j5,j6
!
  real*8      :: z1,z2,phase,xf,fctor
  integer     :: i1,i2,i3,i4,i5,i6,i7,numin,numax,jy1,jy2,jy3,jy4,jy5,jy6,nu
      six_j=0.d0
      z1=delr(j1,j2,j5)
      if(z1.gt.1.d+19) goto 90
      z1=delr(j3,j4,j5)+z1
      if(z1.gt.1.d+19) goto 90
      z2=delr(j1,j3,j6)
      if(z2.gt.1.d+19) goto 90
      z2=delr(j2,j4,j6)+z2
      if(z2.gt.1.d+19) goto 90
      z1=0.5d0*(z1+z2)
      i1=(j1+j2+j4+j3)/2+1
      i2=(j1+j2-j5)/2
      i3=(j4+j3-j5)/2
      i4=(j1+j3-j6)/2
      i5=(j4+j2-j6)/2
      i6=(j1+j4-j5-j6)/2
      i7=(j2+j3-j5-j6)/2
      numin=max0(i6,i7,0)
      numax=min0(i2,i3,i4,i5)
      if(numax.lt.numin) go to 90
      phase=phasef(numin+(j1+j2+j3+j4)/2)
      do nu=numin,numax
        jy1=i2-nu
        jy2=i3-nu
        jy3=i4-nu
        jy4=i5-nu
        jy5=nu-i6
        jy6=nu-i7
        xf=f(jy1)+f(jy2)+f(jy3)+f(jy4)+f(jy5)+f(jy6)
        fctor=xf-z1+f(nu)-f(i1-nu)
        fctor=dexp(fctor)
        six_j=six_j+phase/fctor
        phase=-phase
      end do
90    return
end function six_j


function racah(i1,i2,i3,i4,i5,i6)
  implicit none
  real*8    :: racah
  integer   :: i1,i2,i3,i4,i5,i6
  integer   :: ipf
  ipf=(i1+i2+i5+i4)/2
!
! call six_j and calculate racah
!
  racah=six_j(i1,i2,i3,i4,i5,i6)/phasef(ipf)
  return
end function racah


function three_j(j1,j2,j3,m1,m2)
!
!  three_j is the 3-j symbol. The arguments are obvious, e.g.
!  m3=-m1-m2. j's and m's are doubled.
!  Note that max(j1+j2+j3+1)=2000 (i.e. twice the sum of j's),
!  otherwise change f(0:1000))
!
  implicit none
  integer :: j1,j2,j3,m1,m2
  real*8  :: three_j
!
  integer :: m3,i1,i2,i3,i4,i5,i6,i7,kmin,kmax,isj,k,j1m1,j2m2,j3m3m,j3m3p
  real*8  :: sum,fk,fjm,fsj,den,fj,comfac

      three_j=0.d0
      i1=(j1+j2-j3)/2
      if(i1.lt.0)return
      i2=(j1-j2+j3)/2
      if(i2.lt.0)return
      i3=(-j1+j2+j3)/2
      if(i3.lt.0)return
      fj=f(i1)+f(i2)+f(i3)
      m3=-m1-m2
      i4=(j1-m1)/2
      i5=(j2+m2)/2
      j1m1=(j1+m1)/2
      j2m2=(j2-m2)/2
      j3m3p=(j3+m3)/2
      j3m3m=(j3-m3)/2
      fjm=f(j1m1)+f(i4)+f(i5)+f(j2m2)+f(j3m3p)+f(j3m3m)
      isj=(j1+j2+j3+2)/2
      fsj=f(isj)
      comfac=0.5d0*(fj+fjm-fsj)
      sum=0.d0
      kmax=i1
      if(kmax.gt.i4)kmax=i4
      if(kmax.gt.i5)kmax=i5
      kmin=0
      i6=-(j3-j2+m1)/2
      if(kmin.lt.i6)kmin=i6
      i7=(-j3+j1+m2)/2
      if(kmin.lt.i7)kmin=i7
      fk=phasef(kmin-1)
      do k=kmin,kmax
        fk=-fk
        den=-f(k)-f(i1-k)-f(i4-k)-f(i5-k)-f(-i6+k)-f(-i7+k)
        sum=sum+fk*dexp(den+comfac)
      end do
      three_j=sum*phasef((j1-j2-m3)/2)
      return
end function three_j


function phasef(n)
  implicit none
  real*8  :: phasef
  integer :: n
  phasef=1.d0
  if(n.ne.2*(n/2))phasef=-1.d0
end function phasef



END  MODULE ANGLIB

