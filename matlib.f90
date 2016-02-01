MODULE MATLIB 
USE ANGLIB
USE FUNLIB
implicit none
!
! Based on J. Phys. B 30 2529
!
double precision,allocatable  :: f_fun(:),fl_fun(:,:),il_fun(:,:)
complex*16,allocatable        :: sphabcd(:,:)

CONTAINS

subroutine init_matlib

call init_angular

end subroutine init_matlib





function cnl(n,l,a)
! Eq. (9)  
implicit none
  integer           :: n,l
  double precision  :: a,cnl

  cnl=dexp(x2n(n)+f(n))/(2.d0*a)**(n+l+1.5d0)

end function cnl

function norm(l,a)
! Eq. (10)
implicit none
  integer           :: l
  double precision  :: a,norm
  double precision  :: w

  w=sqrt(pi)*dexp(df(2*l+1)-x2n(l+1))

  norm=sqrt(2*(2.d0*a)**(l+1.5d0)/w)

end function norm


Subroutine precal_ILM(lm,xi,r)
! Eq. (19)
implicit none
  integer           :: lm
  double precision  :: xi,r
  integer           :: l,n,ll
  double precision  :: a,w1


  a=xi*r**2
  il_fun=0.d0
  do ll=0,lm
  do l=0,ll
    if((-1)**(ll-l)>0) then
      n=(ll-l)/2
      il_fun(l,ll)=(-1)**n*cnl(n,l,0.25d0/xi)*exp(-a)*Laguerre(n,l+0.5d0,a)
    endif
  end do
  end do
  
end subroutine precal_ilm  



Subroutine ILM(l1,m1,l2,m2,xi,r,flm)
! Eq. (19)
implicit none
! flm(l)*r**l*Y_lm(r)
  integer           :: l1,m1,l2,m2
  double precision  :: xi,r,flm(0:l1+l2)
  integer           :: l,n,lm1,lm2
  double precision  :: a,w1
!  double precision, external :: Laguerre  

    lm1=l1**2+(l1+1)+m1
    lm2=l2**2+(l2+1)+m2



!  a=xi*r**2
  flm=0.d0
  do l=iabs(l1-l2),l1+l2
    if((-1)**(l1+l2-l)>0) then
!      w1=gaunt(l2,l1,m1,l,m2-m1)
      w1=gaunt_coeff(l,lm1,lm2)
!      n=(l1+l2-l)/2
!      flm(l)=(-1)**n*w1*cnl(n,l,0.25d0/xi)*exp(-a)*Laguerre(n,l+0.5d0,a)
      flm(l)=w1*il_fun(l,l1+l2)
    endif
  end do
  
end subroutine ilm  


! ILM1
! Description: 
! Parameters:	
! Post: 
Subroutine ILM1(l1,m1,l2,m2,xi,r,flm)
! Eq. (19)
implicit none
! flm(l)*r**l*Y_lm(r)
  integer           :: l1,m1,l2,m2
  double precision  :: xi,r,flm(0:l1+l2)
  integer           :: l,n,lm1,lm2
  double precision  :: a,w1
!  double precision, external :: Laguerre  

    lm1=l1**2+(l1+1)+m1
    lm2=l2**2+(l2+1)+m2



  a=xi*r**2
  flm=0.d0
  do l=iabs(l1-l2),l1+l2
    if((-1)**(l1+l2-l)>0) then
!      w1=gaunt(l2,l1,m1,l,m2-m1)
      w1=gaunt_coeff(l,lm1,lm2)
      n=(l1+l2-l)/2
      flm(l)=(-1)**n*w1*cnl(n,l,0.25d0/xi)*exp(-a)*Laguerre(n,l+0.5d0,a)
    endif
  end do
  
end subroutine ilm1  


!!!!!!!!!!NOT USED!!!!!!!!!!
! TLM1
! Description: Tensor operator integral of the tensor operator Y_{k\mu}(\nabla)
!	it is described in equation 23 and 24 of the paper
! Parameters	: w1 - Gaunt coefficient that takes 4 terms (ie <l2m2,k\mu|l1m1,lm>)
!		; 
Subroutine TLM1(l1,m1,l2,m2,k,mu,xi,r,flm)
! Eq. (24)
implicit none
! flm(l)*r**l*Y_lm(r)
  integer           :: l1,m1,l2,m2,k,mu
  double precision  :: xi,r,flm(0:l1+l2+k)
  integer           :: l,n,lm1,lm2
  double precision  :: a,w1
!  double precision, external :: Laguerre  

    lm1=l1**2+(l1+1)+m1
    lm2=l2**2+(l2+1)+m2


  a=xi*r**2
  flm=0.d0
  do l=0,l1+l2+k
    if((-1)**(l1+l2+k-l)>0) then
      w1=gauntc2(l2,m2,k,mu,l1,m1,l)
      n=(l1+l2+k-l)/2
      flm(l)=(-1)**n*w1*cnl(n,l,0.25d0/xi)*exp(-a)*Laguerre(n,l+0.5d0,a)
    endif
  end do
  
end subroutine Tlm1  




Subroutine KLM(l1,m1,l2,m2,xi,r,flm)
! Eq. (22)
implicit none
! flm(l)*r**l*Y_lm(r)
  integer           :: l1,m1,l2,m2
  double precision  :: xi,r,flm(0:l1+l2)
  integer           :: l,n
  double precision  :: a,w1
!  double precision, external :: Laguerre  

  a=xi*r**2
  flm=0.d0
  do l=iabs(l1-l2),l1+l2
    if((-1)**(l1+l2-l)>0) then
      w1=gaunt(l2,l1,m1,l,m2-m1)
      n=(l1+l2-l)/2
      flm(l)=(-1)**n*w1*cnl(n+1,l,0.25d0/xi)*exp(-a)*Laguerre(n+1,l+0.5d0,a)
    endif
  end do
  
end subroutine klm  

Subroutine JLM(l1,m1,l2,m2,xi,r,flm)
! Eq. (33)
implicit none
! flm(l)*r**l*Y_lm(r)
  integer           :: l1,m1,l2,m2
  double precision  :: xi,r,flm(0:l1+l2)
  integer           :: l,n
  double precision  :: a,w1,al,cn,w2,ww
!  double precision, external :: Laguerre  

  a=xi*r**2
  al=0.25d0/xi
  flm=0.d0
  do l=iabs(l1-l2),l1+l2-2
    if((-1)**(l1+l2-l)>0) then
      w1=gaunt(l2,l1,m1,l,m2-m1)
      n=(l1+l2-l)/2
      cn=2**(n-1)*dexp(f(n-1))/(2.d0*al)**(n-1+l+1.5d0)
      w2=exp(-a)*Laguerre(n-1,l+0.5d0,a)*cn
      flm(l)=(-1)**n*w1*w2
    endif
  end do
  l=l1+l2
  w1=gaunt(l2,l1,m1,l,m2-m1)
  if(a==0.d0) then
!  smallgamma(s,x)/x^s=1/s x----> 0
    ww=1.d0/(l+0.5d0)
  else
    ww=1.d0/a**(l+0.5d0)*gammp(l+0.5d0,a)*sqrt(pi)*dexp(f(2*l)-x2n(2*l)-f(l))
  endif  
  flm(l)=w1*(2.d0*xi)**(l+0.5d0)/2.d0*ww
end subroutine jlm


Subroutine ULM(l1,m1,l2,m2,l3,m3,l4,m4,flm)
! Eq. (37)
implicit none
! flm(l)*r**l*Y_lm(r)
  integer           :: l1,m1,l2,m2,l3,m3,l4,m4
  integer           :: lm1,lm2,lm3,lm4
  double precision  :: flm(0:l1+l2+l3+l4)
  integer           :: l,n

 
  lm1=l1**2+(l1+1)+m1
  lm2=l2**2+(l2+1)+m2
  lm3=l3**2+(l3+1)+m3
  lm4=l4**2+(l4+1)+m4

  do l=0,l1+l2+l3+l4
    gauntc3(l)=gaunt3_coeff(l,lm2,lm4,lm1,lm3)
  end do

  flm=0.d0
  do l=0,l1+l2+l3+l4-2
    flm(l)=gauntc3(l)*fl_fun(l,l1+l2+l3+l4)    
  end do
  l=l1+l2+l3+l4
  flm(l)=gauntc3(l)*f_fun(l)    
end subroutine ulm  



Subroutine precal_ULM(lm,xi,r)
implicit none
  integer           :: lm
  double precision  :: xi,r
  integer           :: l,n,ll
  double precision  :: a,w1,al,w2,cn,ww

  a=xi*r**2
  al=0.25d0/xi
  f_fun=0.d0
  fl_fun=0.d0
  do ll=0,lm
    do l=0,ll-2
      if((-1)**(ll-l)>0) then
        n=(ll-l)/2
          cn=2**(n-1)*dexp(f(n-1))/(2.d0*al)**(n-1+l+1.5d0)
          w2=exp(-a)*Laguerre(n-1,l+0.5d0,a)*cn
          fl_fun(l,ll)=(-1)**n*w2
      endif
    end do
  end do
  a=xi*r**2
  do l=0,lm
    if(a==0.d0) then
!    smallgamma(s,x)/x^s=1/s x----> 0
      ww=1.d0/(l+0.5d0)
    else
      ww=1.d0/a**(l+0.5d0)*gammp(l+0.5d0,a)*sqrt(pi)*dexp(f(2*l)-x2n(2*l)-f(l))
    endif  
    f_fun(l)=(2.d0*xi)**(l+0.5d0)/2.d0*ww
  end do  
end subroutine precal_ulm  


Subroutine ULM1(l1,m1,l2,m2,l3,m3,l4,m4,xi,r,flm)
! Eq. (37)
implicit none
! flm(l)*r**l*Y_lm(r)
  integer           :: l1,m1,l2,m2,l3,m3,l4,m4
  double precision  :: xi,r,flm(0:l1+l2+l3+l4)
  integer           :: l,n
  double precision  :: a,w1,al,w2,cn,ww
!  double precision, external :: Laguerre  

  call gaunt3(l2,m2,l4,m4,l1,m1,l3,m3)

  a=xi*r**2
  al=0.25d0/xi
  flm=0.d0
  do l=0,l1+l2+l3+l4-2
    if(gauntc3(l)/=0.d0) then
      n=(l1+l2+l3+l4-l)/2
      cn=2**(n-1)*dexp(f(n-1))/(2.d0*al)**(n-1+l+1.5d0)
      w2=exp(-a)*Laguerre(n-1,l+0.5d0,a)*cn
      flm(l)=(-1)**n*gauntc3(l)*w2
    endif
  end do
  l=l1+l2+l3+l4
  w1=gauntc3(l)
  if(a==0.d0) then
!  smallgamma(s,x)/x^s=1/s x----> 0
    ww=1.d0/(l+0.5d0)
  else
    ww=1.d0/a**(l+0.5d0)*gammp(l+0.5d0,a)*sqrt(pi)*dexp(f(2*l)-x2n(2*l)-f(l))
  endif  
  flm(l)=w1*(2.d0*xi)**(l+0.5d0)/2.d0*ww
    
end subroutine ulm1  

! okmat
! Description:
! Parameters:	w - constants that go in front of Ilm in eq 18 for Slm
! 		xi - constant Xi in eq 18 
! 		r - constant for the difference in the center of the gaussians
! 		rr - norm of the position vector
! 		sph -
! 		kilm
! 		ovlm
! 		ovm
! 		kim
Subroutine okmat(la,ma,lb,mb,alpha,beta,a,b,ovm,kim)
! Eqs. (18) and (21) 
  implicit none
  integer            :: la,ma,lb,mb
  double precision   :: alpha,beta,a(3),b(3)
  complex*16         :: ovm,kim,sph
  double precision   :: xi,r(3),ovlm(0:la+lb),kilm(0:la+lb),rr,w
  integer            :: l

  w=(-1)**lb*(2.d0*pi)**1.5d0*norm(la,0.25d0/alpha)*norm(lb,0.25d0/beta)
  xi=alpha*beta/(alpha+beta)
  r=b-a
  rr=sqrt(r(1)**2+r(2)**2+r(3)**2)
  call ILM1(la,ma,lb,mb,xi,rr,ovlm)
  call KLM(la,ma,lb,mb,xi,rr,kilm)
  ovm=(0.d0,0.d0)
  kim=(0.d0,0.d0)
  do l=iabs(la-lb),la+lb
    sph=spherical_harmonics(r,l,mb-ma)
    ovm=ovm+ovlm(l)*sph*power(rr,l)*w
    kim=kim+kilm(l)*sph*power(rr,l)*w
  end do
end subroutine okmat

Subroutine tenmat(la,ma,lb,mb,alpha,beta,a,b,k,ten)
! Eqs. (23) 
! k=1, x,y,z components
  implicit none
  integer            :: la,ma,lb,mb
  double precision   :: alpha,beta,a(3),b(3)
  complex*16         :: sph,ten(2*k+1)
  double precision   :: xi,r(3),flm(0:la+lb+k),rr,w
  integer            :: l,k,mu

  w=(-1)**lb*(2.d0*pi)**1.5d0*norm(la,0.25d0/alpha)*norm(lb,0.25d0/beta)
  xi=alpha*beta/(alpha+beta)
  r=b-a
  rr=sqrt(r(1)**2+r(2)**2+r(3)**2)
  ten=(0.d0,0.d0)
  do mu=-1,1
    call TLM1(la,ma,lb,mb,k,mu,xi,rr,flm)
    do l=0,la+lb+k
      sph=spherical_harmonics(r,l,mb-ma+mu)
      ten(k+1+mu)=ten(k+1+mu)+flm(l)*sph*power(rr,l)*w
    end do
  end do
  

end subroutine tenmat




subroutine realsp(m,n,w,mm,bk)
  implicit none
  integer              :: m,n,mm(0:1),bk
  complex*16           :: w(0:1)
  complex*16,parameter :: zi=(0.d0,1.d0)

  if(m==0) then
    n=0
    w(0)=(1.d0,0.d0)
    mm(0)=0
  endif
  if(bk==1) then
    if(m<0) then
      n=1
      w(0)=-zi/sqrt(2.d0)
      w(1)=+(-1)**m*zi/sqrt(2.d0)
      mm(0)=abs(m)
      mm(1)=-abs(m)
    endif
    if(m>0) then
      n=1
      w(0)=1/sqrt(2.d0)
      w(1)=(-1)**m/sqrt(2.d0)
      mm(0)=abs(m)
      mm(1)=-abs(m)
    endif
  endif


  if(bk==2) then
    if(m<0) then
      n=1
      w(0)=zi/sqrt(2.d0)
      w(1)=-(-1)**m*zi/sqrt(2.d0)
      mm(0)=abs(m)
      mm(1)=-abs(m)
    endif
    if(m>0) then
      n=1
      w(0)=1/sqrt(2.d0)
      w(1)=(-1)**m/sqrt(2.d0)
      mm(0)=abs(m)
      mm(1)=-abs(m)
    endif
  endif
  
end subroutine realsp





Subroutine erepmat_real(la,ma,lb,mb,lc,mc,ld,md,alpha,beta,gamma,delta,a,b,c,d,numr)
! Eqs. (34-35) 
! modified for real spherical harmonics
! needs optimization
  implicit none
  integer              :: la,ma,lb,mb,lc,mc,ld,md,ico
  double precision     :: alpha,beta,gamma,delta,a(3),b(3),c(3),d(3)
  double precision     :: numr
  complex*16           :: num,wa(0:1),wb(0:1),wc(0:1),wd(0:1)
  integer              :: ia,na,ib,nb,ic,nc,id,nd
  integer              :: mma(0:1),mmb(0:1),mmc(0:1),mmd(0:1)
  real     :: t1,t2
  integer              :: k

  call realsp(ma,na,wa,mma,2)
  call realsp(mb,nb,wb,mmb,1)
  call realsp(mc,nc,wc,mmc,2)
  call realsp(md,nd,wd,mmd,1)

  allocate(sphabcd(0:la+lb+lc+ld,-(la+lb+lc+ld):la+lb+lc+ld))


  numr=0.d0
  do ia=0,na
    do ib=0,nb
      do ic=0,nc
        do id=0,nd
          ico=ia+ib+ic+id
!            call erepmat(ico,la,mma(ia),lb,mmb(ib),lc,mmc(ic),ld,mmd(id),alpha,beta,gamma,delta,a,b,c,d,num)
            numr=numr+wa(ia)*wb(ib)*wc(ic)*wd(id)*num
        end do
      end do
    end do
  end do

  deallocate(sphabcd)
  
end subroutine erepmat_real

Subroutine erepmat(la,ma,lb,mb,lc,mc,ld,md,alpha,beta,gamma,delta,a,b,c,d,num)
! Eqs. (34-35) 
  implicit none
  integer                 :: la,ma,lb,mb,lc,mc,ld,md,ico
  double precision        :: alpha,beta,gamma,delta,a(3),b(3),c(3),d(3)
  complex*16              :: num,sph,sui1,sui2,suu
  double precision        :: xi1,xi2,zeta1,zeta2,rba(3),rdc(3),flm(0:la+lb+lc+ld),iilm1(0:la+lb),iilm2(0:lc+ld)
  double precision        :: rrba,rrdc,rabcd(3),rrabcd
  double precision        :: w,wa,wb,wc,wd,zeta,fc
  double precision        :: fwa,fwb,fwc,fwd
  integer                 :: lap,map,lapp,mapp,lbp,mbp,lbpp,mbpp,mamin,mamax,mbmin,mbmax,l,m
  integer                 :: lcp,mcp,lcpp,mcpp,ldp,mdp,ldpp,mdpp,mcmin,mcmax,mdmin,mdmax
  integer                 :: lmma,lmmap,lmmb,lmmbp,lmmc,lmmcp,lmmd,lmmdp
  complex*16,allocatable  :: part1(:,:,:,:),part2(:,:,:,:)
  complex*16              :: sphab(0:la+lb,-(la+lb):la+lb),sphcd(0:lc+ld,-(lc+ld):lc+ld)

  allocate(part1(-2*lb:2*lb,0:lb,-2*la:2*la,0:la),part2(-2*ld:2*ld,0:ld,-2*lc:2*lc,0:lc))

  allocate(sphabcd(0:la+lb+lc+ld,-(la+lb+lc+ld):la+lb+lc+ld))



  xi1=alpha*beta/(alpha+beta)
  zeta1=alpha+beta
  xi2=gamma*delta/(gamma+delta)
  zeta2=gamma+delta
  rba=b-a
  rdc=d-c
  rrba=sqrt(rba(1)**2+rba(2)**2+rba(3)**2)
  rrdc=sqrt(rdc(1)**2+rdc(2)**2+rdc(3)**2)
  zeta=zeta1*zeta2/(zeta1+zeta2)



!  if(ico==0) then
  w=32.d0*(-1)**(lb+ld)*(2.d0*pi)**6.5d0* &
&   norm(la,0.25d0/alpha)*norm(lb,0.25d0/beta)*norm(lc,0.25d0/gamma)*norm(ld,0.25d0/delta)
  rabcd=(alpha*a+beta*b)/(alpha+beta)-(gamma*c+delta*d)/(gamma+delta)
  rrabcd=sqrt(rabcd(1)**2+rabcd(2)**2+rabcd(3)**2)




  do l=0,la+lb+lc+ld
    fc=power(rrabcd,l)
    do m=-l,l
      sphabcd(l,m)=spherical_harmonics(rabcd,l,m)*fc*w
    end do
  end do

  call precal_ULM(la+lb+lc+ld,zeta,rrabcd)
!  endif



  do l=0,la+lb
    fc=power(rrba,l)
    do m=-l,l
      sphab(l,m)=spherical_harmonics(rba,l,m)*fc
    end do
  end do

  do l=0,lc+ld
    fc=power(rrdc,l)
    do m=-l,l
      sphcd(l,m)=spherical_harmonics(rdc,l,m)*fc
    end do
  end do




  part1=(0.d0,0.d0)
  part2=(0.d0,0.d0)



  call precal_ILM(la+lb,xi1,rrba)
  do lap=0,la
    fwa=(-alpha/zeta1)**lap
    lapp=la-lap
    mamin=max(-lap,ma-la+lap)
    mamax=min(lap,ma+la-lap)
    do map=mamin,mamax
      mapp=ma-map
!      wa=glm(la,ma,lap,map)*(-alpha/zeta1)**lap
      lmma=la**2+(la+1)+ma
      lmmap=lap**2+(lap+1)+map
      wa=glm_c(lmma,lmmap)*fwa
      do lbp=0,lb
        fwb=(beta/zeta1)**lbp
        lbpp=lb-lbp
        mbmin=max(-lbp,mb-lb+lbp)
        mbmax=min(lbp,mb+lb-lbp)
        do mbp=mbmin,mbmax
          mbpp=mb-mbp
          lmmb=lb**2+(lb+1)+mb
          lmmbp=lbp**2+(lbp+1)+mbp
          wb=glm_c(lmmb,lmmbp)*fwb
          call ILM(lapp,mapp,lbpp,mbpp,xi1,rrba,iilm1)
          sui1=(0.d0,0.d0)   
          do l=iabs(lapp-lbpp),lapp+lbpp
!            sph=spherical_harmonics(rba,l,mbpp-mapp)
!            sui1=sui1+iilm1(l)*sph*power(rrba,l)
            sui1=sui1+iilm1(l)*sphab(l,mbpp-mapp)
          end do
          part1(mbp,lbp,map,lap)=sui1*wa*wb
        end do
      end do
    end do
  end do

  call precal_ILM(lc+ld,xi2,rrdc)
  do lcp=0,lc
    fwc=(gamma/zeta2)**lcp
    lcpp=lc-lcp
    mcmin=max(-lcp,mc-lc+lcp)
    mcmax=min(lcp,mc+lc-lcp)
    do mcp=mcmin,mcmax
      mcpp=mc-mcp
      lmmc=lc**2+(lc+1)+mc
      lmmcp=lcp**2+(lcp+1)+mcp
      wc=glm_c(lmmc,lmmcp)*fwc
      do ldp=0,ld
        fwd=(-delta/zeta2)**ldp
        ldpp=ld-ldp
        mdmin=max(-ldp,md-ld+ldp)
        mdmax=min(ldp,md+ld-ldp)
        do mdp=mdmin,mdmax
          mdpp=md-mdp
          lmmd=ld**2+(ld+1)+md
          lmmdp=ldp**2+(ldp+1)+mdp
          wd=glm_c(lmmd,lmmdp)*fwd
          call ILM(lcpp,mcpp,ldpp,mdpp,xi2,rrdc,iilm2)
          sui2=(0.d0,0.d0)   
          do l=iabs(lcpp-ldpp),lcpp+ldpp
!            sph=spherical_harmonics(rdc,l,mdpp-mcpp)
!            sui2=sui2+iilm2(l)*sph*power(rrdc,l)
            sui2=sui2+iilm2(l)*sphcd(l,mdpp-mcpp)
          end do
          part2(mdp,ldp,mcp,lcp)=sui2*wc*wd
        end do
      end do
    end do
  end do

  
  num=(0.d0,0.d0)
  do lap=0,la
    lapp=la-lap
    mamin=max(-lap,ma-la+lap)
    mamax=min(lap,ma+la-lap)
    do map=mamin,mamax
      mapp=ma-map
      do lbp=0,lb
        lbpp=lb-lbp
        mbmin=max(-lbp,mb-lb+lbp)
        mbmax=min(lbp,mb+lb-lbp)
        do mbp=mbmin,mbmax
          mbpp=mb-mbp
          
          do lcp=0,lc
            lcpp=lc-lcp
            mcmin=max(-lcp,mc-lc+lcp)
            mcmax=min(lcp,mc+lc-lcp)
            do mcp=mcmin,mcmax
              mcpp=mc-mcp
              do ldp=0,ld
                ldpp=ld-ldp
                mdmin=max(-ldp,md-ld+ldp)
                mdmax=min(ldp,md+ld-ldp)
                do mdp=mdmin,mdmax
                  mdpp=md-mdp
                  call ULM(lap,map,lbp,mbp,lcp,mcp,ldp,mdp,flm)
                  suu=(0.d0,0.d0)   
                  do l=0,lap+lbp+lcp+ldp
                    suu=suu+flm(l)*sphabcd(l,mbp+mdp-map-mcp)
                  end do
                  sui1=part1(mbp,lbp,map,lap)
                  sui2=part2(mdp,ldp,mcp,lcp)
                  num=num+sui1*sui2*suu
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do



  deallocate(part1,part2)
  deallocate(sphabcd)

  
end subroutine erepmat



Subroutine erepmat1(la,ma,lb,mb,lc,mc,ld,md,alpha,beta,gamma,delta,a,b,c,d,num)
! Eqs. (34-35) 
  implicit none
  integer            :: la,ma,lb,mb,lc,mc,ld,md
  double precision   :: alpha,beta,gamma,delta,a(3),b(3),c(3),d(3)
  complex*16         :: num,sph,sui1,sui2,suu
  double precision   :: xi1,xi2,zeta1,zeta2,rba(3),rdc(3),flm(0:la+lb+lc+ld),iilm1(0:la+lb),iilm2(0:lc+ld)
  double precision   :: rrba,rrdc,rabcd(3),rrabcd
  double precision   :: w,wa,wb,wc,wd,zeta
  integer            :: lap,map,lapp,mapp,lbp,mbp,lbpp,mbpp,mamin,mamax,mbmin,mbmax,l
  integer            :: lcp,mcp,lcpp,mcpp,ldp,mdp,ldpp,mdpp,mcmin,mcmax,mdmin,mdmax
  

  w=32.d0*(-1)**(lb+ld)*(2.d0*pi)**6.5d0* &
&   norm(la,0.25d0/alpha)*norm(lb,0.25d0/beta)*norm(lc,0.25d0/gamma)*norm(ld,0.25d0/delta)
  xi1=alpha*beta/(alpha+beta)
  zeta1=alpha+beta
  xi2=gamma*delta/(gamma+delta)
  zeta2=gamma+delta
  rba=b-a
  rdc=d-c
  rrba=sqrt(rba(1)**2+rba(2)**2+rba(3)**2)
  rrdc=sqrt(rdc(1)**2+rdc(2)**2+rdc(3)**2)
  zeta=zeta1*zeta2/(zeta1+zeta2)
  rabcd=(alpha*a+beta*b)/(alpha+beta)-(gamma*c+delta*d)/(gamma+delta)
  rrabcd=sqrt(rabcd(1)**2+rabcd(2)**2+rabcd(3)**2)
  num=(0.d0,0.d0)
  do lap=0,la
    lapp=la-lap
    mamin=max(-lap,ma-la+lap)
    mamax=min(lap,ma+la-lap)
    do map=mamin,mamax
      mapp=ma-map
      wa=glm(la,ma,lap,map)*(-alpha/zeta1)**lap
      do lbp=0,lb
        lbpp=lb-lbp
        mbmin=max(-lbp,mb-lb+lbp)
        mbmax=min(lbp,mb+lb-lbp)
        do mbp=mbmin,mbmax
          mbpp=mb-mbp
          wb=glm(lb,mb,lbp,mbp)*(beta/zeta1)**lbp
          call ILM(lapp,mapp,lbpp,mbpp,xi1,rrba,iilm1)
          sui1=(0.d0,0.d0)   
          do l=iabs(lapp-lbpp),lapp+lbpp
            sph=spherical_harmonics(rba,l,mbpp-mapp)
            sui1=sui1+iilm1(l)*sph*power(rrba,l)
          end do

          do lcp=0,lc
            lcpp=lc-lcp
            mcmin=max(-lcp,mc-lc+lcp)
            mcmax=min(lcp,mc+lc-lcp)
            do mcp=mcmin,mcmax
              mcpp=mc-mcp
              wc=glm(lc,mc,lcp,mcp)*(gamma/zeta2)**lcp
              do ldp=0,ld
                ldpp=ld-ldp
                mdmin=max(-ldp,md-ld+ldp)
                mdmax=min(ldp,md+ld-ldp)
                do mdp=mdmin,mdmax
                  mdpp=md-mdp
                  wd=glm(ld,md,ldp,mdp)*(-delta/zeta2)**ldp
                  call ILM(lcpp,mcpp,ldpp,mdpp,xi2,rrdc,iilm2)
                  sui2=(0.d0,0.d0)   
                  do l=iabs(lcpp-ldpp),lcpp+ldpp
                    sph=spherical_harmonics(rdc,l,mdpp-mcpp)
                    sui2=sui2+iilm2(l)*sph*power(rrdc,l)
                  end do
                  call ULM1(lap,map,lbp,mbp,lcp,mcp,ldp,mdp,zeta,rrabcd,flm)
                  suu=(0.d0,0.d0)   
                  do l=0,lap+lbp+lcp+ldp
                    sph=spherical_harmonics(rabcd,l,mbp+mdp-map-mcp)
                    suu=suu+flm(l)*sph*power(rrabcd,l)
                  end do
                  num=num+w*wa*wb*wc*wd*sui1*sui2*suu
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
  
end subroutine erepmat1






Subroutine nucmat(la,ma,lb,mb,alpha,beta,a,b,c,num)
! Eqs. (30-33) 
  implicit none
  integer            :: la,ma,lb,mb
  double precision   :: alpha,beta,a(3),b(3),c(3)
  complex*16         :: num,sph,sui,suj
  double precision   :: xi,zeta,r(3),jjlm(0:la+lb),iilm(0:la+lb),rr,w,wa,wb,rabc(3),rrabc
  integer            :: lap,map,lapp,mapp,lbp,mbp,lbpp,mbpp,mamin,mamax,mbmin,mbmax,l
  

  w=8.d0*(-1)**lb*(2.d0*pi)**3.d0*norm(la,0.25d0/alpha)*norm(lb,0.25d0/beta)
  xi=alpha*beta/(alpha+beta)
  zeta=alpha+beta
  r=b-a
  rr=sqrt(r(1)**2+r(2)**2+r(3)**2)
  rabc=(alpha*a+beta*b)/(alpha+beta)-c
  rrabc=sqrt(rabc(1)**2+rabc(2)**2+rabc(3)**2)
  num=(0.d0,0.d0)
  do lap=0,la
    lapp=la-lap
    mamin=max(-lap,ma-la+lap)
    mamax=min(lap,ma+la-lap)
    do map=mamin,mamax
      mapp=ma-map
      wa=glm(la,ma,lap,map)*(-alpha/zeta)**lap
      do lbp=0,lb
        lbpp=lb-lbp
        mbmin=max(-lbp,mb-lb+lbp)
        mbmax=min(lbp,mb+lb-lbp)
        do mbp=mbmin,mbmax
          mbpp=mb-mbp
          wb=glm(lb,mb,lbp,mbp)*(beta/zeta)**lbp
          call ILM1(lapp,mapp,lbpp,mbpp,xi,rr,iilm)
          call JLM(lap,map,lbp,mbp,zeta,rrabc,jjlm)
          sui=(0.d0,0.d0)   
          do l=iabs(lapp-lbpp),lapp+lbpp
            sph=spherical_harmonics(r,l,mbpp-mapp)
            sui=sui+iilm(l)*sph*power(rr,l)
          end do
          suj=(0.d0,0.d0)   
          do l=iabs(lap-lbp),lap+lbp
            sph=spherical_harmonics(rabc,l,mbp-map)
            suj=suj+jjlm(l)*sph*power(rrabc,l)
          end do
          num=num+w*wa*wb*sui*suj
        end do
      end do
    end do
  end do
end subroutine nucmat



Subroutine okmat_real(la,ma,lb,mb,alpha,beta,ra,rb,ovmr,kimr)
! Eqs. (18) and (21) 
! modified for real spherical harmonics
! needs optimization
  implicit none
  integer              :: la,ma,lb,mb
  double precision     :: alpha,beta,ra(3),rb(3),ovmr,kimr
  integer              :: mma,mmb
  complex*16           :: of(4),kf(4)
  complex*16           :: ovm,kim
  complex*16,parameter :: zi=(0.d0,1.d0)
  integer              :: ia,na,ib,nb
  integer              :: m_ma(0:1),m_mb(0:1)
  complex*16           :: wa(0:1),wb(0:1)



!  mma=abs(ma)
!  mmb=abs(mb)
!!   for real spherical harmonics
!  call okmat(la,mma,lb,mmb,alpha,beta,ra,rb,of(1),kf(1))
!  call okmat(la,mma,lb,-mmb,alpha,beta,ra,rb,of(2),kf(2))
!  of(2)=of(2)*(-1)**mb
!  kf(2)=kf(2)*(-1)**mb
!  call okmat(la,-mma,lb,mmb,alpha,beta,ra,rb,of(3),kf(3))
!  of(3)=of(3)*(-1)**ma
!  kf(3)=kf(3)*(-1)**ma
!  call okmat(la,-mma,lb,-mmb,alpha,beta,ra,rb,of(4),kf(4))
!  of(4)=of(4)*(-1)**(ma+mb)
!  kf(4)=kf(4)*(-1)**(ma+mb)
!  call tr_sp(ma,mb,of,ovmr)
!  call tr_sp(ma,mb,kf,kimr)
!  write(200,*)ovmr,kimr

  call realsp(ma,na,wa,m_ma,2)
  call realsp(mb,nb,wb,m_mb,1)

  ovmr=0.d0
  kimr=0.d0
  do ia=0,na
    do ib=0,nb
      call okmat(la,m_ma(ia),lb,m_mb(ib),alpha,beta,ra,rb,ovm,kim)
      ovmr=ovmr+wa(ia)*wb(ib)*ovm
      kimr=kimr+wa(ia)*wb(ib)*kim
    end do
  end do
!  write(200,*)ovmr,kimr
!  write(200,*)'------'





    
end subroutine okmat_real

subroutine tr_sp(ma,mb,f,num)
! transform complex m.e. to real m.e.
implicit none
  integer              :: ma,mb
  double precision     :: num
  complex*16           :: f(4)
  complex*16,parameter :: zi=(0.d0,1.d0)
    
  if(ma.eq.0.and.mb.eq.0) then
    num=f(1)
  endif
  if(ma.gt.0.and.mb.gt.0) then
    num=0.5d0*(f(1)+f(2)+f(3)+f(4))
  endif
  if(ma.eq.0.and.mb.gt.0) then
    num=1.d0/sqrt(2.d0)*(f(1)+f(2))
  endif
  if(ma.eq.0.and.mb.lt.0) then
    num=zi/sqrt(2.d0)*(-f(1)+f(2))
  endif
  if(ma.gt.0.and.mb.eq.0) then
    num=1.d0/sqrt(2.d0)*(f(1)+f(3))
  endif
  if(ma.lt.0.and.mb.eq.0) then
    num=zi/sqrt(2.d0)*(f(1)-f(3))
  endif
  if(ma.lt.0.and.mb.lt.0) then
    num=-0.5d0*(-f(1)+f(2)+f(3)-f(4))
  endif
  if(ma.gt.0.and.mb.lt.0) then
    num=zi*0.5d0*(-f(1)+f(2)-f(3)+f(4))
  endif
  if(ma.lt.0.and.mb.gt.0) then
    num=zi*0.5d0*(f(1)+f(2)-f(3)-f(4))
  endif
end subroutine tr_sp





Subroutine nucmat_real(la,ma,lb,mb,alpha,beta,ra,rb,rc,numr)
! Eqs. (30-33) 
! modified for real spherical harmonics
! needs optimization
  implicit none
  integer              :: la,ma,lb,mb
  double precision     :: alpha,beta,ra(3),rb(3),rc(3),numr
  integer              :: mma,mmb
  complex*16           :: angf(4),num
  integer              :: ia,na,ib,nb
  integer              :: m_ma(0:1),m_mb(0:1)
  complex*16           :: wa(0:1),wb(0:1)

!
!!
!!
!
!  mma=abs(ma)
!  mmb=abs(mb)
!   for real spherical harmonics
!  call nucmat(la,mma,lb,mmb,alpha,beta,ra,rb,rc,angf(1))
!  call nucmat(la,mma,lb,-mmb,alpha,beta,ra,rb,rc,angf(2))
!  angf(2)=angf(2)*(-1)**mb
!  call nucmat(la,-mma,lb,mmb,alpha,beta,ra,rb,rc,angf(3))
!  angf(3)=angf(3)*(-1)**ma
!  call nucmat(la,-mma,lb,-mmb,alpha,beta,ra,rb,rc,angf(4))
!  angf(4)=angf(4)*(-1)**(ma+mb)

!  call tr_sp(ma,mb,angf,numr)

!  write(100,*)la,ma,lb,mb
!  write(100,*)numr

  call realsp(ma,na,wa,m_ma,2)
  call realsp(mb,nb,wb,m_mb,1)

  numr=0.d0
  do ia=0,na
    do ib=0,nb
      call nucmat(la,m_ma(ia),lb,m_mb(ib),alpha,beta,ra,rb,rc,num)
      numr=numr+wa(ia)*wb(ib)*num
    end do
  end do
!  write(100,*)numr 
!  write(100,*)'------'
 
return
    
end subroutine nucmat_real

! precal_gaunt
! Description: 
subroutine precal_gaunt(lao)
implicit none
integer                       :: l1,m1,l2,m2,l3,m3,l4,m4,lm1,lm2,lm3,lm4,lm,lmmax,l,lao


lmmax=2*lao
lm=lao**2+(lao+1)+lao
allocate(gaunt_coeff(0:lmmax,0:lm,0:lm))
allocate(glm_c(0:lm,0:lm))
gaunt_coeff=0.d0
do l1=0,lao
  do m1=-l1,l1
    lm1=l1**2+(l1+1)+m1
    do l2=0,lao
      do m2=-l2,l2
        lm2=l2**2+(l2+1)+m2
        glm_c(lm1,lm2)=glm(l1,m1,l2,m2)
        do l=iabs(l1-l2),l1+l2
          gaunt_coeff(l,lm1,lm2)=gaunt(l2,l1,m1,l,m2-m1)
        end do
      end do
    end do
  end do
end do


allocate(f_fun(0:lmmax),fl_fun(0:lmmax,0:lmmax))
allocate(il_fun(0:lmmax,0:lmmax))


end subroutine precal_gaunt



END MODULE MATLIB 
