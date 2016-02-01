MODULE FUNLIB 
implicit none

CONTAINS



function Laguerre(n,alpha,x)
  implicit none
  real*8   :: Laguerre,alpha,x,bin
  integer  :: i,n
  Laguerre=1.d0; bin=1.d0
  do i=n,1,-1
    bin= bin*(alpha+i)/dfloat(n+1-i)
    Laguerre= bin-x*Laguerre/dfloat(i)
  end do
end function Laguerre




!***********************************************************************
      function bessj(l,x)
!
!     Calculates bessj(l,x) = x*j(l,x), where j(l,x) is the Spherical
!     Bessel Function of integer order l.
!     l=0,2 : The analytic expressions in sin(x) & cos(x) are used.
!     l>2   : The functional value is obtained by downward recursion,
!             if x<l, by upward recursion if x>l, calling 'besrec()'.
!     x<xmin,l>0: bessj()=x**(l+1)/(2*l-1)!!
!     If the argument is excessively small, 'bestiny' sets the output
!     equal to the lowest order term in the powerseries expansion.
!***********************************************************************
      implicit  none
      real*8,parameter    :: xmin=1.0d-10
      real*8              :: bessj,x
      integer             :: l
!      real*8,external     :: besrec,bestiny

      if(l.eq.0)then
         bessj=sin(x)
      else if(x.lt.xmin)then
         bessj=x*bestiny(l,x)
      else if(l.eq.1)then
         bessj=sin(x)/x-cos(x)
      else if(l.eq.2)then
         bessj=(3.0d0/x/x-1.0d0)*sin(x)-3.0d0*cos(x)/x
      else
         bessj=x*besrec(l,x)
      endif

      return
      end function bessj


!***********************************************************************
      function besrec(l,x)
!
!     Recursion for Bessel functions. 'acclog' is about the negative
!     log of the desired accuracy.
!***********************************************************************
      implicit none
      real*8, parameter       :: acclog=1.0d1
      real*8                  :: besrec,x
      integer                 :: l
      real*8                  :: upj,downj,tempj
      integer                 :: m,mstart

      if(x.gt.dble(l))then
!
! upward recursion
!
         downj=sin(x)/x
         tempj=(downj-cos(x))/x
         do 100 m=2,l
            upj=(2.0*m-1.0)/x*tempj-downj
                 downj=tempj
                 tempj=upj
100      continue
         besrec=upj
      else
!
! downward recursion & renormalization
!
         mstart=2*l+int(acclog*sqrt(dble(l)))
200      upj=0.0d0
         tempj=1.0d0
         do 300 m=mstart-1,0,-1
            downj=(2.0*m+3.0)/x*tempj-upj
            if(m.eq.l)besrec=downj
            upj=tempj
            tempj=downj
300      continue
         besrec=besrec/downj*sin(x)/x
      endif

400   return
      end function besrec


!***********************************************************************
      function bestiny(l,x)
!
!     Lowest order powerseries term for spherical Bessel function.
!***********************************************************************
      implicit none
      real*8                   :: bestiny,x
      integer                  :: l,k
      real*8                   :: kdiv
      kdiv=1.0d0
      do 100 k=1,2*l-1,2
         kdiv=kdiv*k
100   continue
      bestiny=x**l/kdiv

      return
      end function bestiny



      FUNCTION gammp(a,x)
      real*8 :: a,gammp,x
      real*8 :: gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END function gammp


      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      real*8 :: a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
      INTEGER i
      real*8 :: an,b,c,d,del,h
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END subroutine gcf

      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      real*8 :: a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
      INTEGER n
      real*8 :: ap,del,sum
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END subroutine gser
      

      FUNCTION gammln(xx)
      real*8 :: gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0, &
      -1.231739572450155d0,.1208650973866179d-2, &
&       -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END function gammln

END MODULE FUNLIB 
