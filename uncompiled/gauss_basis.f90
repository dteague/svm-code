MODULE GAUSS_BASIS
  implicit none

  integer,parameter                     :: nmax=100,nnmax=301,n2d=100
  real*8                                :: f(0:nmax),df(-1:nnmax),x2n(0:n2d)
  real*8                                :: dfa(0:2*nmax+1)
  
  real*8,parameter                      :: pi_=3.141592653589793d0
  integer,parameter                     :: lexp=2  !,ngd=(lexp+1)**2*ng
  integer                               :: select
  real*8,dimension(:),allocatable       :: xgr,wgr
  integer,parameter                     :: npow=10,npow_max=100
  real*8,dimension(:),allocatable       :: fa

  integer,parameter                     :: nge=50
!  real*8,parameter                      :: xe0=1.16d0,ae0=0.001d0
   real*8,parameter                      :: xe0=1.15d0,ae0=0.005d0
   real*8                                :: gae(nge)

  real*8,dimension(:,:),allocatable     :: pot_ge
  real*8,dimension(:,:,:),allocatable   :: nl_pot_u,nl_pot_uv
  complex*16                            :: stm(-lexp:lexp,-lexp:lexp),som(-lexp:lexp,-lexp:lexp)
  real*8                                :: spm(-lexp:lexp,-lexp:lexp)
  complex*16                            :: nam(3,-lexp:lexp,-lexp:lexp)
  real*8                                :: cpm(-lexp:lexp,-lexp:lexp)
  complex*16                            :: hom(-lexp:lexp,-lexp:lexp)
  complex*16				:: lim(-lexp:lexp,-lexp:lexp)
  real*8                                :: gis(0:lexp+2,0:lexp+2)
  complex*16                            :: clm(npow,(lexp+1)**2)
  integer                               :: px(npow,(lexp+1)**2),py(npow,(lexp+1)**2),pz(npow,(lexp+1)**2), &
&                                          nlm((lexp+1)**2)
  integer,parameter                     :: n_c=89
  double precision                      :: nu_c(n_c),w_c(n_c)
  double precision,parameter            :: factor=0.1d0

  complex*16,parameter                  :: zi=(0.d0,1.d0)


CONTAINS


! power 
! Description: finds x**n
! parameters: 	x, n
! post: returns x**n taking into acount special cases of 0**n and 0**0 (=1)
function power(x,n)
  implicit none
  real*8       :: power,x
  integer      :: n
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


! gauss_init
! Description: 
! Parameters:	fa - factorial or fa(x) = x!
!		in - number that starts at 1 and goes up, assigning a different number to 
!			to each lm set (ie in(00)=1, in(1-1)=2, in(10)=3... etc)
!		lexp - the maximum l value for the gauss basis
!		sphar -
!		dfa - function makes a product of odd numbers (ie dfa(x) = (2x+1).(2x-1)...5.3.1) or x!!
!		nlm
!		clm
!		clma
!		pxa,pya,pza
!		px,py,pz
subroutine gauss_init
  implicit none
  integer                                  :: i,l,m,na,in,j
  complex*16                               :: clma(npow)
  real*8                                   :: x
  integer                                  :: pxa(npow),pya(npow),pza(npow)
  allocate(fa(0:npow_max))

 !	Create factorial function or x!
  fa(0)=1
  do i=1,npow_max
    fa(i)=fa(i-1)*i
  end do

  do l=0,lexp
    do m=-l,l
      in=l**2+l+m+1
      call sphar(l,m,na,clma,pxa,pya,pza)
      nlm(in)=na
      clm(:,in)=clma(:)
      px(:,in)=pxa(:)
      py(:,in)=pya(:)
      pz(:,in)=pza(:)
    end do
  end do

  call  gequad

! 	create the function of a product of odd numbers or x!!
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

end subroutine  gauss_init




subroutine gauss_potential_m(a,ra,la,b,rb,lb,c,rc)
!
!  exp(-mu*(r-Rc)**2)
!
implicit none
  integer             :: la,ma,lb,mb,ina,inb
  real*8              :: a,b,c,ra(3),rb(3),rc(3)

  integer             :: ia,ib,na,nb
  complex*16          :: clma(npow),clmb(npow),s
  integer             :: pxa(npow),pya(npow),pza(npow)
  integer             :: pxb(npow),pyb(npow),pzb(npow)
  real*8              :: w,cx,cy,cz,sx,sy,sz,xa,xb,xc
  real*8              :: sxm(0:lexp+2,0:lexp+2),sym(0:lexp+2,0:lexp+2),szm(0:lexp+2,0:lexp+2)

  spm=0.d0

  xc=a*b/(a+b+c)*((ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2)
  xb=a*c/(a+b+c)*((ra(1)-rc(1))**2+(ra(2)-rc(2))**2+(ra(3)-rc(3))**2)
  xa=c*b/(a+b+c)*((rc(1)-rb(1))**2+(rc(2)-rb(2))**2+(rc(3)-rb(3))**2)

  w=norma(la,a)*norma(lb,b)*(pi_/(a+b+c))**1.5d0*exp(-xa-xb-xc)
  call gauss_int3(a,b,c,ra(1),rb(1),rc(1),la,lb)
  sxm=gis
  call gauss_int3(a,b,c,ra(2),rb(2),rc(2),la,lb)
  sym=gis
  call gauss_int3(a,b,c,ra(3),rb(3),rc(3),la,lb)
  szm=gis



  do ma=-la,la
    ina=la**2+la+ma+1
    na=nlm(ina)
    clma(:)=Conjg(clm(:,ina))
    pxa(:)=px(:,ina)
    pya(:)=py(:,ina)
    pza(:)=pz(:,ina)

    do mb=-lb,lb
      inb=lb**2+lb+mb+1
      nb=nlm(inb)
      clmb(:)=clm(:,inb)
      pxb(:)=px(:,inb)
      pyb(:)=py(:,inb)
      pzb(:)=pz(:,inb)

      s=(0.d0,0.d0)
      do ia=1,na
	do ib=1,nb
	  if(clma(ia)*clmb(ib).ne.0.d0) then
	    sx=sxm(pxa(ia),pxb(ib))
	    sy=sym(pya(ia),pyb(ib))
	    sz=szm(pza(ia),pzb(ib))
	    s=s+clma(ia)*clmb(ib)*w*sx*sy*sz
	  endif
	end do
      end do
      spm(ma,mb)=s
    end do
  end do

end subroutine gauss_potential_m


subroutine Coulomb_potential_m(a,ra,la,b,rb,lb,rc)
!
!  exp(-mu*(r-Rc)**2)
!
implicit none
  integer             :: la,ma,lb,mb,ina,inb,i_c
  real*8              :: a,b,c,ra(3),rb(3),rc(3)

  integer             :: ia,ib,na,nb
  complex*16          :: clma(npow),clmb(npow),s
  integer             :: pxa(npow),pya(npow),pza(npow)
  integer             :: pxb(npow),pyb(npow),pzb(npow)
  real*8              :: w,cx,cy,cz,sx,sy,sz,xa,xb,xc
  real*8              :: sxm(0:lexp+2,0:lexp+2),sym(0:lexp+2,0:lexp+2),szm(0:lexp+2,0:lexp+2)

  cpm=0.d0
  do i_c=1,n_c
    c=nu_c(i_c)*factor**2

    xc=a*b/(a+b+c)*((ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2)
    xb=a*c/(a+b+c)*((ra(1)-rc(1))**2+(ra(2)-rc(2))**2+(ra(3)-rc(3))**2)
    xa=c*b/(a+b+c)*((rc(1)-rb(1))**2+(rc(2)-rb(2))**2+(rc(3)-rb(3))**2)
  
    w=norma(la,a)*norma(lb,b)*(pi_/(a+b+c))**1.5d0*exp(-xa-xb-xc)
    call gauss_int3(a,b,c,ra(1),rb(1),rc(1),la,lb)
    sxm=gis
    call gauss_int3(a,b,c,ra(2),rb(2),rc(2),la,lb)
    sym=gis
    call gauss_int3(a,b,c,ra(3),rb(3),rc(3),la,lb)
    szm=gis

    do ma=-la,la
      ina=la**2+la+ma+1
      na=nlm(ina)
      clma(:)=Conjg(clm(:,ina))
      pxa(:)=px(:,ina)
      pya(:)=py(:,ina)
      pza(:)=pz(:,ina)

      do mb=-lb,lb
        inb=lb**2+lb+mb+1
        nb=nlm(inb)
        clmb(:)=clm(:,inb)
        pxb(:)=px(:,inb)
        pyb(:)=py(:,inb)
        pzb(:)=pz(:,inb)

        s=(0.d0,0.d0)
        do ia=1,na
          do ib=1,nb
            if(clma(ia)*clmb(ib).ne.(0.d0,0.d0)) then
              sx=sxm(pxa(ia),pxb(ib))
              sy=sym(pya(ia),pyb(ib))
              sz=szm(pza(ia),pzb(ib))
              s=s+clma(ia)*clmb(ib)*w*sx*sy*sz
            endif
          end do
        end do
        cpm(ma,mb)=cpm(ma,mb)+s*factor*w_c(i_c)
      end do
    end do
  end do
  

end subroutine Coulomb_potential_m


subroutine gauss_kinetic_m(a,ra,la,b,rb,lb)
implicit none
  integer             :: la,ma,lb,mb,ina,inb
  real*8              :: a,b,ra(3),rb(3)

  integer             :: ia,ib,na,nb
  complex*16          :: clma(npow),clmb(npow)
  integer             :: pxa(npow),pya(npow),pza(npow)
  integer             :: pxb(npow),pyb(npow),pzb(npow)
  real*8              :: w,cx,cy,cz,sx,sy,sz,xxx,ho,li 
  real*8              :: tx1,tx2,ty1,ty2,tz1,tz2,tx,ty,tz
  real*8              :: qx1,qx2,qy1,qy2,qz1,qz2,qx,qy,qz
  real*8              :: sxm(0:lexp+2,0:lexp+2),sym(0:lexp+2,0:lexp+2),szm(0:lexp+2,0:lexp+2)

  xxx=a*b/(a+b)*((ra(1)-rb(1))**2+(ra(2)-rb(2))**2+(ra(3)-rb(3))**2)
  w=norma(la,a)*norma(lb,b)*(pi_/(a+b))**1.5d0*exp(-xxx)
  som=(0.d0,0.d0)
  nam=(0.d0,0.d0)
  hom=(0.d0,0.d0)
  lim=(0.d0,0.d0)
  stm=(0.d0,0.d0)
  
  tx1=0.d0; tx2=0.d0
  ty1=0.d0; ty2=0.d0
  tz1=0.d0; tz2=0.d0
  qx1=0.d0; qx2=0.d0
  qy1=0.d0; qy2=0.d0
  qz1=0.d0; qz2=0.d0

  call gauss_int(a,b,ra(1),rb(1),la,lb+2)
  sxm=gis
  call gauss_int(a,b,ra(2),rb(2),la,lb+2)
  sym=gis
  call gauss_int(a,b,ra(3),rb(3),la,lb+2)
  szm=gis



  do ma=-la,la
    ina=la**2+la+ma+1
    na=nlm(ina)
    clma(:)=Conjg(clm(:,ina))
    pxa(:)=px(:,ina)
    pya(:)=py(:,ina)
    pza(:)=pz(:,ina)

    do mb=-lb,lb
      inb=lb**2+lb+mb+1
      nb=nlm(inb)
      clmb(:)=clm(:,inb)
      pxb(:)=px(:,inb)
      pyb(:)=py(:,inb)
      pzb(:)=pz(:,inb)


      do ia=1,na
	do ib=1,nb
	  if(clma(ia)*clmb(ib).ne.(0.d0,0.d0)) then
  
	    sx=sxm(pxa(ia),pxb(ib))
	    sy=sym(pya(ia),pyb(ib))
	    sz=szm(pza(ia),pzb(ib))
! overlap
	    som(ma,mb)=som(ma,mb)+clma(ia)*clmb(ib)*w*sx*sy*sz
      
! Kinetic
	    if(pxb(ib)*(pxb(ib)-1).ne.0) tx1=sxm(pxa(ia),pxb(ib)-2)
	    tx2=sxm(pxa(ia),pxb(ib)+2)
	    tx=(pxb(ib)*(pxb(ib)-1)*tx1-2.d0*b*(2*pxb(ib)+1)*sx+4.d0*b**2*tx2)*sy*sz

	    if(pyb(ib)*(pyb(ib)-1).ne.0) ty1=sym(pya(ia),pyb(ib)-2)
	    ty2=sym(pya(ia),pyb(ib)+2)
	    ty=(pyb(ib)*(pyb(ib)-1)*ty1-2.d0*b*(2*pyb(ib)+1)*sy+4.d0*b**2*ty2)*sx*sz

	    if(pzb(ib)*(pzb(ib)-1).ne.0) tz1=szm(pza(ia),pzb(ib)-2)
	    tz2=szm(pza(ia),pzb(ib)+2)
	    tz=(pzb(ib)*(pzb(ib)-1)*tz1-2.d0*b*(2*pzb(ib)+1)*sz+4.d0*b**2*tz2)*sx*sy
	    stm(ma,mb)=stm(ma,mb)-clma(ia)*clmb(ib)*w*(tx+ty+tz)

!  nabla

	    if(pxb(ib).ne.0) qx1=sxm(pxa(ia),pxb(ib)-1)
	    qx2=sxm(pxa(ia),pxb(ib)+1)
	    qx=(pxb(ib)*qx1-2.d0*b*qx2)*sy*sz

	    if(pyb(ib).ne.0) qy1=sym(pya(ia),pyb(ib)-1)
	    qy2=sym(pya(ia),pyb(ib)+1)
	    qy=(pyb(ib)*qy1-2.d0*b*qy2)*sx*sz

	    if(pzb(ib).ne.0) qz1=szm(pza(ia),pzb(ib)-1)
	    qz2=szm(pza(ia),pzb(ib)+1)
	    qz=(pzb(ib)*qz1-2.d0*b*qz2)*sy*sx

	    nam(1,ma,mb)=nam(1,ma,mb)+clma(ia)*clmb(ib)*w*qx
	    nam(2,ma,mb)=nam(2,ma,mb)+clma(ia)*clmb(ib)*w*qy
	    nam(3,ma,mb)=nam(3,ma,mb)+clma(ia)*clmb(ib)*w*qz
	    
!  H.O.
	    ho=tx2*sy*sz+ty2*sx*sz+tz2*sx*sy

	    hom(ma,mb)=hom(ma,mb)+clma(ia)*clmb(ib)*w*ho
	    
! linear x only.  only works for sa=sb=0
	    li=qx2*sy*sz
! 	    li=qy2*sx*sz
! 	    li=qz2*sx*sy

            lim(ma,mb)=lim(ma,mb)+clma(ia)*clmb(ib)*w*li

	  endif
	end do
      end do
    end do
  end do
end subroutine gauss_kinetic_m


subroutine gauss_int(a,b,xa,xb,na,nb)
!
!   integrate((x-xa)**na*(x-xb)**nb*exp(-a*(x-xa)**2-b*(x-xb)**2)
!
  implicit none
  integer            :: na,nb,ka,kb,i,k
  real*8             :: a,b,xa,xb,s,wa,wb,c,pc(0:nb),x,xx
  real*8             :: h(0:na+nb)

  c=xb-xa
  x=b*c/sqrt(a+b)
  xx=0.5d0/sqrt(a+b)
  call c_hermite_pol(x,na+nb,h)
  do i=0,nb
    pc(i)=power(-c,i)
  end do

  gis=0.d0
  do kb=0,nb
    do ka=0,na
      s=0.d0
      do i=0,kb
        s=s+fa(kb)/(fa(i)*fa(kb-i))*pc(i)*power(xx,(ka+kb-i))*h(ka+kb-i)
      end do
      gis(ka,kb)=s
    end do
  end do


end subroutine gauss_int

subroutine gauss_int3(a,b,c,xa,xb,xc,na,nb)
!
!   integrate((x-xa)**na*(x-xb)**nb*exp(-a*(x-xa)**2-b*(x-xb)**2-c*(x-xc)**2)
!
  implicit none
  integer            :: na,nb,ia,ib,ka,kb,i
  real*8             :: a,b,c,xa,xb,xc,s,wa,wb,pc(0:nb),ba,ca,xx,x,h(0:na+nb)


  ba=xb-xa
  ca=xc-xa
  x=(b*ba+c*ca)/sqrt(a+b+c)
  xx=0.5d0/sqrt(a+b+c)
  call c_hermite_pol(x,na+nb,h)
  do i=0,nb
    pc(i)=power(-ba,i)
  end do

  gis=0.d0
  do kb=0,nb
    do ka=0,na
      s=0.d0
      do i=0,kb
        s=s+fa(kb)/(fa(i)*fa(kb-i))*pc(i)*power(xx,(ka+kb-i))*h(ka+kb-i)
      end do
      gis(ka,kb)=s
    end do
  end do

end subroutine gauss_int3


function norma(l,a)
  implicit none
  real*8                     :: a,norma,f
  integer                    :: i,l

  f=1.d0
  do i=1,l
    f=f*(2.d0*i+1)
  end do
  f=f*sqrt(pi_)/2.d0**(l+1)
  norma=sqrt(2.d0*(2.d0*a)**(l+1.5d0)/f)/sqrt(4.d0*pi_)
end function norma


subroutine basint(a,b,n,c)
  implicit none
  real*8                                :: a,b,x,ww,aa
  integer                               :: n,i
  real*8                                :: c(0:n),h(0:n)

  aa=0.5d0/sqrt(a)
  x=b*aa
  call c_hermite_pol(x,n,h)
  ww=1.d0/aa
  do i=0,n
    ww=ww*aa
    c(i)=ww*h(i)
  end do

end subroutine basint


! c_hermite_pol finds the nth hermite polynomial 
! parameters;	x - variable for hermite
!		n - which hermite polynomial
!		h(x) - hermite polynomial, takes 1 variable
! post: returns the nth hermite polynomial at value x
subroutine c_hermite_pol(x,n,h)
  implicit none
  integer                               :: n,i
  real*8                                :: x
  real*8                                :: h(0:n)

  h=0.d0
  h(0)=1.d0
  if(n.ge.1) h(1)=2.d0*x
  if(n.ge.2) then
    do i=1,n-1
      h(i+1)=2.d0*(x*h(i)+i*h(i-1))
    end do
  endif

end subroutine c_hermite_pol

! sphar
! Description:
! Parameters:	clm -
! 		w1
! 		w2
! 		w3
! 		w4
! 		
subroutine sphar(l,m1,n,clm,px,py,pz)
  implicit none
  complex*16            :: clm(npow)
  real*8                :: w1,w2,w3,w4
  integer               :: px(npow),py(npow),pz(npow)
  integer               :: l,m,n,p,q,i,j,k,m1


  clm=(0.d0,0.d0)
  m=m1
  w1=sqrt((2.d0*l+1.d0)*fa(l+m)*fa(l-m))
  n=0
  do q=0,l
    do p=0,l
      if(p+q.le.l.and.p-q.eq.m) then
        w2=(-1)**p/(fa(p)*fa(q)*fa(l-p-q))/2**(p+q)
        if(m.eq.0) then
          do i=0,p
            n=n+1
            w3=fa(p)/(fa(i)*fa(p-i))
            clm(n)=w1*w2*w3
            px(n)=2*i
            py(n)=2*(p-i)
            pz(n)=l-2*p
          end do
        else
          do i=0,p
            do j=0,q
              n=n+1
              w3=fa(p)/(fa(i)*fa(p-i))*fa(q)/(fa(j)*fa(q-j))*(-1)**(q-j)
              k=p+q-i-j
! for real sph
!              if(m1.gt.0) w4=(zi**k+(-zi)**k)/sqrt(2.d0)
!              if(m1.lt.0) w4=(zi**k-(-zi)**k)/sqrt(2.d0)/zi
!              clm(n)=w1*w2*w3*w4
! for complex sph
!              w4=zi**k
              clm(n)=w1*w2*w3*zi**k
              px(n)=i+j
              py(n)=p+q-i-j
              pz(n)=l-p-q
            end do
          end do
        end if
      end if
    end do
  end do

end subroutine sphar


subroutine potential_fit(rp,pot,np,nu,ng,cg)
! fit the local potentials with gaussians
  implicit none
  integer              :: ierr,i,j,k,np,ng
  real*8               :: pot(np),rp(np),cg(ng)
  real*8               :: am(ng,ng),bm(ng),bf(ng)
  real*8               :: r2,nu(ng)




  am=0.d0
  bm=0.d0
  do k=1,np
    r2=rp(k)
    do i=1,ng
      bf(i)=exp(-nu(i)*r2**2)
    end do
    do i=1,ng
      do j=1,ng
        am(j,i)=am(j,i)+bf(i)*bf(j)
      end do
      bm(i)=bm(i)+bf(i)*pot(k)
    end do
  end do

  call gaussj(am,ng,ng,bm,1,1,ierr)

  cg=bm

end subroutine potential_fit










  subroutine gequad
! 
  implicit none        
  double precision      :: urange,drange,acc
!
!       range [10^(-9),1] and accuracy ~10^(-8);
!

         nu_c(1)=4.96142640560223544d19
         nu_c(2)=1.37454269147978052d19
         nu_c(3)=7.58610013441204679d18
         nu_c(4)=4.42040691347806996d18
         nu_c(5)=2.61986077948367892d18
         nu_c(6)=1.56320138155496681d18
         nu_c(7)=9.35645215863028402d17
         nu_c(8)=5.60962910452691703d17
         nu_c(9)=3.3666225119686761d17
         nu_c(10)=2.0218253197947866d17
         nu_c(11)=1.21477756091902017d17
         nu_c(12)=7.3012982513608503d16
         nu_c(13)=4.38951893556421099d16
         nu_c(14)=2.63949482512262325d16
         nu_c(15)=1.58742054072786174d16
         nu_c(16)=9.54806587737665531d15
         nu_c(17)=5.74353712364571709d15
         nu_c(18)=3.455214877389445d15
         nu_c(19)=2.07871658520326804d15
         nu_c(20)=1.25064667315629928d15
         nu_c(21)=7.52469429541933745d14
         nu_c(22)=4.5274603337253175d14
         nu_c(23)=2.72414006900059548d14
         nu_c(24)=1.63912168349216752d14
         nu_c(25)=9.86275802590865738d13
         nu_c(26)=5.93457701624974985d13
         nu_c(27)=3.5709554322296296d13
         nu_c(28)=2.14872890367310454d13
         nu_c(29)=1.29294719957726902d13
         nu_c(30)=7.78003375426361016d12
         nu_c(31)=4.68148199759876704d12
         nu_c(32)=2.8169955024829868d12
         nu_c(33)=1.69507790481958464d12
         nu_c(34)=1.01998486064607581d12
         nu_c(35)=6.13759486539856459d11
         nu_c(36)=3.69320183828682544d11
         nu_c(37)=2.22232783898905102d11
         nu_c(38)=1.33725247623668682d11
         nu_c(39)=8.0467192739036288d10
         nu_c(40)=4.84199582415144143d10
         nu_c(41)=2.91360091170559564d10
         nu_c(42)=1.75321747475309216d10
         nu_c(43)=1.0549735552210995d10
         nu_c(44)=6.34815321079006586d9
         nu_c(45)=3.81991113733594231d9
         nu_c(46)=2.29857747533101109d9
         nu_c(47)=1.38313653595483694d9
         nu_c(48)=8.32282908580025358d8
         nu_c(49)=5.00814519374587467d8
         nu_c(50)=3.01358090773319025d8
         nu_c(51)=1.81337994217503535d8
         nu_c(52)=1.09117589961086823d8
         nu_c(53)=6.56599771718640323d7
         nu_c(54)=3.95099693638497164d7
         nu_c(55)=2.37745694710665991d7
         nu_c(56)=1.43060135285912813d7
         nu_c(57)=8.60844290313506695d6
         nu_c(58)=5.18000974075383424d6
         nu_c(59)=3.116998193057466d6
         nu_c(60)=1.87560993870024029d6
         nu_c(61)=1.12862197183979562d6
         nu_c(62)=679132.441326077231d0
         nu_c(63)=408658.421279877969d0
         nu_c(64)=245904.473450669789d0
         nu_c(65)=147969.568088321005d0
         nu_c(66)=89038.612357311147d0
         nu_c(67)=53577.7362552358895d0
         nu_c(68)=32239.6513926914668d0
         nu_c(69)=19399.7580852362791d0
         nu_c(70)=11673.5323603058634d0
         nu_c(71)=7024.38438577707758d0
         nu_c(72)=4226.82479307685999d0
         nu_c(73)=2543.43254175354295d0
         nu_c(74)=1530.47486269122675d0
         nu_c(75)=920.941785160749482d0
         nu_c(76)=554.163803906291646d0
         nu_c(77)=333.46029740785694d0
         nu_c(78)=200.6550575335041d0
         nu_c(79)=120.741366914147284d0
         nu_c(80)=72.6544243200329916d0
         nu_c(81)=43.7187810415471025d0
         nu_c(82)=26.3071631447061043d0
         nu_c(83)=15.8299486353816329d0
         nu_c(84)=9.52493152341244004d0
         nu_c(85)=5.72200417067776041d0
         nu_c(86)=3.36242234070940928d0
         nu_c(87)=1.75371394604499472d0
         nu_c(88)=0.64705932650658966d0
         nu_c(89)=0.072765905943708247d0
!
         w_c(1)=47.67445484528304247d10
         w_c(2)=11.37485774750442175d9
         w_c(3)=78.64340976880190239d8
         w_c(4)=46.27335788759590498d8
         w_c(5)=24.7380464827152951d8
         w_c(6)=13.62904116438987719d8
         w_c(7)=92.79560029045882433d8
         w_c(8)=52.15931216254660251d8
         w_c(9)=31.67018011061666244d8
         w_c(10)=1.29291036801493046d8
         w_c(11)=1.00139319988015862d8
         w_c(12)=7.75892350510188341d7
         w_c(13)=6.01333567950731271d7
         w_c(14)=4.66141178654796875d7
         w_c(15)=3.61398903394911448d7
         w_c(16)=2.80225846672956389d7
         w_c(17)=2.1730509180930247d7
         w_c(18)=1.68524482625876965d7
         w_c(19)=1.30701489345870338d7
         w_c(20)=1.01371784832269282d7
         w_c(21)=7.86264116300379329d6
         w_c(22)=6.09861667912273717d6
         w_c(23)=4.73045784039455683d6
         w_c(24)=3.66928949951594161d6
         w_c(25)=2.8462050836230259d6
         w_c(26)=2.20777394798527011d6
         w_c(27)=1.71256191589205524d6
         w_c(28)=1.32843556197737076d6
         w_c(29)=1.0304731275955989d6
         w_c(30)=799345.206572271448d0
         w_c(31)=620059.354143595343d0
         w_c(32)=480986.704107449333d0
         w_c(33)=373107.167700228515d0
         w_c(34)=289424.08337412132d0
         w_c(35)=224510.248231581788d0
         w_c(36)=174155.825690028966d0
         w_c(37)=135095.256919654065d0
         w_c(38)=104795.442776800312d0
         w_c(39)=81291.4458222430418d0
         w_c(40)=63059.0493649328682d0
         w_c(41)=48915.9040455329689d0
         w_c(42)=37944.8484018048756d0
         w_c(43)=29434.4290473253969d0
         w_c(44)=22832.7622054490044d0
         w_c(45)=17711.743950151233d0
         w_c(46)=13739.287867104177d0
         w_c(47)=10657.7895710752585d0
         w_c(48)=8267.42141053961834d0
         w_c(49)=6413.17397520136448d0
         w_c(50)=4974.80402838654277d0
         w_c(51)=3859.03698188553047d0
         w_c(52)=2993.51824493299154d0
         w_c(53)=2322.1211966811754d0
         w_c(54)=1801.30750964719641d0
         w_c(55)=1397.30379659817038d0
         w_c(56)=1083.91149143250697d0
         w_c(57)=840.807939169209188d0
         w_c(58)=652.228524366749422d0
         w_c(59)=505.944376983506128d0
         w_c(60)=392.469362317941064d0
         w_c(61)=304.444930257324312d0
         w_c(62)=236.162932842453601d0
         w_c(63)=183.195466078603525d0
         w_c(64)=142.107732186551471d0
         w_c(65)=110.23530215723992d0
         w_c(66)=85.5113346705382257d0
         w_c(67)=66.3325469806696621d0
         w_c(68)=51.4552463353841373d0
         w_c(69)=39.9146798429449273d0
         w_c(70)=30.9624728409162095d0
         w_c(71)=24.018098812215013d0
         w_c(72)=18.6312338024296588d0
         w_c(73)=14.4525541233150501d0
         w_c(74)=11.2110836519105938d0
         w_c(75)=8.69662175848497178d0
         w_c(76)=6.74611236165731961d0
         w_c(77)=5.23307018057529994d0
         w_c(78)=4.05937850501539556d0
         w_c(79)=3.14892659076635714d0
         w_c(80)=2.44267408211071604d0
         w_c(81)=1.89482240522855261d0
         w_c(82)=1.46984505907050079d0
         w_c(83)=1.14019261330527007d0
         w_c(84)=0.884791217422925293d0
         w_c(85)=0.692686387080616483d0
         w_c(86)=0.585244576897023282d0
         w_c(87)=0.576182522545327589d0
         w_c(88)=0.596688817388997178d0
         w_c(89)=0.607879901151108771d0
!
          urange = 1.d0
          drange=1d-08
          acc   =1d-08
!
  end subroutine gequad

  
  

END MODULE GAUSS_BASIS


      SUBROUTINE gaussj(a,n,np,b,m,mp,ierr)
      implicit real*8(a-h,o-z)
      dimension a(np,np),b(np,mp)
      PARAMETER (NMAX=500)
      dimension indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                ierr=1
                return
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) then
        ierr=1
        return
        endif
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END

