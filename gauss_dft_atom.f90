MODULE GAUSS_DFT_ATOM
  USE GAUSS_BASIS
  USE input_mod
  
  implicit none
! parameters: atomic units
  real*8, parameter            :: E2=1.0,H2M=0.5d0,a_B=1.d0,Ry=0.5d0
! Number of lattice points 
  integer                      :: N_L(3),N_L_points
   
! Lattice index
  integer,allocatable          :: Lattice(:,:),Lattice_inv(:,:,:)
! grid points  
  real*8,allocatable           :: grid_point(:,:)  
! grid spacing
  real*8                       :: grid_step(3),grid_volume
 
! boundary conditions 
  real*8,allocatable           :: V_X0(:,:),V_XN(:,:),V_Y0(:,:),V_YN(:,:),V_Z0(:,:),V_ZN(:,:)
    
  real*8,allocatable           :: rho(:),V_POT(:),wf(:,:,:),V_exchange_up(:)
  real*8,allocatable           :: phi(:),L_phi(:),VH(:),V_ext(:),V_exchange_dw(:)
  real*8,allocatable           :: density(:),density_old(:),density_up(:),H_Phi(:)
  real*8,allocatable           :: density_up_old(:),density_dw(:),density_dw_old(:)
  real*8,allocatable           :: density_new(:)


  
! order of finite  difference
  integer,parameter            :: N_d=4
! max L in multipole expansion
  integer,parameter            :: L_max=4
  real*8,parameter             :: small=1.d-50
! Number of orbitals
  integer                      :: N_orbitals
! single_particle energy
  real*8,allocatable           :: sp_energy(:),Psi(:,:,:)

  integer                      :: N_iteration=4
! Energies 
  real*8                       :: E_hartree,E_exchange,Energy
  integer,parameter            :: N_scf_iter=100
! Gauss basis
  integer                      :: basis_dim
  double precision,allocatable :: basis_nu(:),basis_coef(:,:),basis_r(:,:)
  double precision,allocatable :: dm(:,:),dm_o(:,:),dm_n(:,:),dm_a(:,:),el_rep(:,:,:,:),h_m(:,:),o_m(:,:),l_m(:,:)
  complex*16,allocatable       :: k_m(:,:,:)
  integer,allocatable          :: basis_l(:),basis_m(:)  
! matrix elements
  double precision,allocatable :: mat_vex_up(:,:),mat_vex_dw(:,:),mat_hartree(:,:)
  double precision,allocatable :: mat_hartree_num(:,:)
!
  double precision             :: Z_charge=0.0       !initialized to zero for security  
  double precision,allocatable :: e_sp(:)

 double precision,allocatable  :: c_gs(:)


integer,parameter                    :: N_elements=100
integer                              :: N_Species,N_total_atom
! charge and    index of atoms
integer                              :: z_atom(N_elements),elem(0:N_elements)
! max orbital momentum in basis states
integer,parameter                    :: l_ao_max=2
! number of gaussians in ao expansion
integer                              :: N_gauss
! parameters of gaussians
double precision                     :: a0_gauss,x0_gauss
! expansion coefficient of atomic orbitals
double precision, allocatable        :: ao_basis_c(:,:,:,:)
! number of atomic orbital basis function per atom
integer,allocatable                  :: N_ao(:)
! number of atomic orbital basis function per atom with a given orbital momentum
integer,allocatable                  :: N_ao_l(:,:)
! position of atoms
double precision, allocatable        :: Position(:,:)
! max orbital momentum of a given atom
integer, allocatable                 :: l_max_basis(:)
integer, allocatable                 :: atom_index(:)
! orbital and atom index and inverse table
integer, allocatable                 :: basis_atom(:,:),basis_atom_inv(:,:)
! 
integer, allocatable                 :: basis_ind(:,:,:,:)


! number of grid points for a given atom
  integer, allocatable         :: n_int_points(:)
! integration points for a given atom: defines the index of the integration points
! on the grid
  integer, allocatable         :: int_points_atom(:,:)
! radius of the atomc orbital
  double precision,allocatable :: wf_rad(:)
! basis functions
  complex*16,allocatable       :: Psi_basis(:,:)
! atom pairs with overlapping wave finctions
  integer, allocatable         :: atom_pairs(:,:)
    





CONTAINS


!################################################!
!		Not Used			 !
!################################################!
SUBROUTINE FD_P
  integer :: k1,k2,k3

  rho=-4.d0*pi_*rho
  
  ! Boundary condition in the x direction
  do k2=1,N_L(2)
    do k3=1,N_L(3)
      rho(Lattice_inv(1,k2,k3))=rho(Lattice_inv(1,k2,k3))-V_X0(k2,k3)/grid_step(1)**2/E2
      rho(Lattice_inv(N_L(1),k2,k3))=rho(Lattice_inv(N_L(1),k2,k3))-V_XN(k2,k3)/grid_step(1)**2/E2
    enddo
  enddo
  
  ! Boundary condition in the y direction
  do k1=1,N_L(1)
    do k3=1,N_L(3)
      rho(Lattice_inv(k1,1,k3))=rho(Lattice_inv(k1,1,k3))-V_Y0(k1,k3)/grid_step(2)**2/E2
      rho(Lattice_inv(k1,N_L(2),k3))=rho(Lattice_inv(k1,N_L(2),k3))-V_YN(k1,k3)/grid_step(2)**2/E2
    enddo
  enddo
  
  ! Boundary condition in the z direction
  do k1=1,N_L(1)
    do k2=1,N_L(2)
      rho(Lattice_inv(k1,k2,1))=rho(Lattice_inv(k1,k2,1))-V_Z0(k1,k2)/grid_step(3)**2/E2
      rho(Lattice_inv(k1,k2,N_L(3)))=rho(Lattice_inv(k1,k2,N_L(3)))-V_ZN(k1,k2)/grid_step(3)**2/E2
    enddo
  enddo

  call CoGr(rho,N_L_points)
  v_pot=E2*rho
END SUBROUTINE FD_P


!################################################!
!		Not Used			 !
!################################################!
subroutine calculate_density
  implicit none 
  integer               :: i,j,k,ii,i2,j2
  double precision      :: r(3),x
  complex*16            :: b1,b2,s
  real*4                :: t1,t2
  complex*16,allocatable :: H(:)
  allocate(H(N_L_points))




  dm_a=0.d0
  do k=1,N_orbitals
    do j=1,basis_dim
      do i=1,basis_dim
        dm_a(i,j)=dm_a(i,j)+2.d0*basis_coef(i,k)*basis_coef(j,k)
      end do
    end do
  end do

  density_new=0.d0  
  do i=1,basis_dim
!
    i2=basis_atom(2,i)
    H=0.d0
    do k=1,n_int_points(i2)
      ii=int_points_atom(k,i2)
      H(ii)=psi_basis(k,i)
    end do
!
!    write(6,*)'???',i,sum(Conjg(H(:))*H(:))*grid_volume
    do j=1,basis_dim
!      call cpu_time(t1)
!      s=(0.d0,0.d0)
!      do k=1,N_L_points
!        r(:)=grid_point(:,k)
!        b1=Conjg(basis_func(r,basis_nu(i),basis_l(i),basis_m(i)))
!        b2=basis_func(r,basis_nu(j),basis_l(j),basis_m(j))
!        density_new(k)=density_new(k)+dm_a(i,j)*b1*b2
!        s=s+b1*b2        
!      end do
!      call cpu_time(t2)
!      write(22,*)i,j
!      write(22,*)s*grid_volume
!      write(22,*)t2-t1
!      call cpu_time(t1)
      s=(0.d0,0.d0)
      j2=basis_atom(2,j)
      do k=1,n_int_points(j2)
        ii=int_points_atom(k,j2)
        s=s+Conjg(H(ii))*psi_basis(k,j)
        density_new(ii)=density_new(ii)+dm_a(i,j)*Conjg(H(ii))*psi_basis(k,j)
      end do
!      call cpu_time(t2)
      write(22,*)s*grid_volume
!      write(22,*)t2-t1
    end do
  end do
  density_up=density_new/2.d0
  density_dw=density_new/2.d0

  write(6,*)'???',sum(density_new)*grid_volume,N_orbitals

end subroutine calculate_density



!################################################!
!		Not Used			 !
!################################################!
subroutine matrix_elem
  implicit none 
  integer               :: i,j,k,i2,j2,ii
  double precision      :: r(3),su,sd,sh
  complex*16            :: b1,b2
  complex*16,allocatable :: H(:)
  allocate(H(N_L_points))

  


  do i=1,basis_dim
    i2=basis_atom(2,i)
    H=0.d0
    do k=1,n_int_points(i2)
      ii=int_points_atom(k,i2)
      H(ii)=psi_basis(k,i)
    end do
    do j=1,basis_dim
      su=0.d0
      sd=0.d0
      sh=0.d0
      j2=basis_atom(2,j)
      do k=1,n_int_points(j2)
        ii=int_points_atom(k,j2)
        b1=Conjg(H(ii))
        b2=psi_basis(k,j)
        su=su+b1*b2*v_exchange_up(ii)
        sd=sd+b1*b2*v_exchange_dw(ii)
        sh=sh+b1*b2*vh(ii)
      end do
      mat_vex_dw(j,i)=sd*grid_volume
      mat_vex_up(j,i)=su*grid_volume
      mat_hartree_num(j,i)=sh*grid_volume
    end do
  end do
end subroutine matrix_elem



SUBROUTINE init_lattice
  integer :: k1,k2,k3,num,i,k

!  read(1,*)(N_L(i),i=1,3)
!  read(1,*)(grid_step(i),i=1,3)
  N_L_points=N_L(1)*N_L(2)*N_L(3)
  grid_volume=grid_step(1)*grid_step(2)*grid_step(3)

  allocate(Lattice(3,N_L_points),grid_point(3,N_L_points),Lattice_inv(N_L(1),N_L(2),N_L(3)))


  ! Setup the lattice bookkeeping
  num=0
  do k1=1,N_L(1)
    do k2=1,N_L(2)
      do k3=1,N_L(3)
        num=num+1 
        Lattice(1,num)=k1
        Lattice(2,num)=k2
        Lattice(3,num)=k3
        Lattice_inv(k1,k2,k3)=num
        grid_point(1,num)=-0.5d0*(N_L(1)-1)*grid_step(1)+(k1-1)*grid_step(1)
        grid_point(2,num)=-0.5d0*(N_L(2)-1)*grid_step(2)+(k2-1)*grid_step(2)
        grid_point(3,num)=-0.5d0*(N_L(3)-1)*grid_step(3)+(k3-1)*grid_step(3)
      enddo
    enddo
  enddo
  allocate(density_old(N_L_points),density_new(N_L_points),density(N_L_points),density_up(N_L_points), &
&   density_dw(N_L_points),rho(N_L_points),vh(N_L_points),v_pot(N_L_points), &
&   v_exchange_up(N_L_points),v_exchange_dw(N_L_points))
  allocate(V_X0(N_L(2),N_L(3)),V_XN(N_L(2),N_L(3)),V_Y0(N_L(1),N_L(3)),V_YN(N_L(1),N_L(3)), &
&   V_Z0(N_L(1),N_L(2)),V_ZN(N_L(1),N_L(2)))
  allocate(phi(N_L_points),L_Phi(N_L_Points))
  allocate(wf(-N_d:N_L(1)+N_d,-N_d:N_L(2)+N_d,-N_d:N_L(3)+N_d))
   

END SUBROUTINE init_lattice




  
SUBROUTINE laplace_operator4
  ! 4th order finite difference representation of the Laplacian
  real*8,parameter :: cN0=-205.d0/72.d0,cN1=8.d0/5.d0
  real*8,parameter :: cN2=-1.d0/5.d0,cN3=8.d0/315.d0, cN4=-1.d0/560.d0
  integer          :: i,i1,i2,i3
  real*8           :: K_x,K_y,K_z

  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i); i3=Lattice(3,i)
    wf(i1,i2,i3)=Phi(i)
  enddo

  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i); i3=Lattice(3,i)
    K_x=(cN0* wf(i1,i2,i3)+& 
             cN1*(wf(i1+1,i2,i3)+wf(i1-1,i2,i3))+& 
             cN2*(wf(i1+2,i2,i3)+wf(i1-2,i2,i3))+& 
             cN3*(wf(i1+3,i2,i3)+wf(i1-3,i2,i3))+& 
             cN4*(wf(i1+4,i2,i3)+wf(i1-4,i2,i3)))/Grid_Step(1)**2
    K_y=(cN0* wf(i1,i2,i3)+& 
             cN1*(wf(i1,i2+1,i3)+wf(i1,i2-1,i3))+& 
             cN2*(wf(i1,i2+2,i3)+wf(i1,i2-2,i3))+& 
             cN3*(wf(i1,i2+3,i3)+wf(i1,i2-3,i3))+& 
             cN4*(wf(i1,i2+4,i3)+wf(i1,i2-4,i3)))/Grid_Step(2)**2
    K_z=(cN0* wf(i1,i2,i3)+& 
             cN1*(wf(i1,i2,i3+1)+wf(i1,i2,i3-1))+& 
             cN2*(wf(i1,i2,i3+2)+wf(i1,i2,i3-2))+& 
             cN3*(wf(i1,i2,i3+3)+wf(i1,i2,i3-3))+& 
             cN4*(wf(i1,i2,i3+4)+wf(i1,i2,i3-4)))/Grid_Step(3)**2
    L_Phi(i)=K_x+K_y+K_z
  enddo
END SUBROUTINE laplace_operator4



SUBROUTINE laplace_operator
  ! 1st order finite difference representation of the Laplacian
  integer :: i1,i2,i3,kpoint,i
  real*8  :: k_x,K_y,K_z
  wf=0.d0    
  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i); i3=Lattice(3,i)
    wf(i1,i2,i3)=phi(i)
  end do
  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i); i3=Lattice(3,i)
    K_x=(-2.d0*wf(i1,i2,i3)+wf(i1+1,i2,i3)+wf(i1-1,i2,i3))/grid_step(1)**2
    K_y=(-2.d0*wf(i1,i2,i3)+wf(i1,i2+1,i3)+wf(i1,i2-1,i3))/grid_step(2)**2
    K_z=(-2.d0*wf(i1,i2,i3)+wf(i1,i2,i3+1)+wf(i1,i2,i3-1))/grid_step(3)**2
    L_phi(i)=K_x+K_y+K_z
  enddo
END SUBROUTINE laplace_operator



!################################################!
!		Not Used			 !
!################################################!
SUBROUTINE CoGr(b,NL)
  ! Solution of the Poisson equation using the conjugate gradient method
	integer                          :: NL,iteration,i,j
	real*8                           :: c0,c1,alfa,beta,rr,bn,con,b(NL)
	integer,parameter                :: N_iter=1500
	real*8,parameter                 :: eps=1.d-10
	real*8,dimension(:),allocatable  :: g,d,h,x

	allocate(g(NL),d(NL),h(NL),x(NL))

	bn=sqrt(dot_product(b,b))
	x=0.d0
	phi=x
	if(N_d==4) then
	  call laplace_operator4
	else
	  call laplace_operator
	endif
	g=b+L_phi
	d=g
	c0=dot_product(g,g)

	do iteration=1,N_iter
		con=abs(sqrt(dot_product(g,g))/bn)
		if(con>eps) then
		  phi=d
		  	if(N_d==4) then
					call laplace_operator4
				else
					call laplace_operator
				endif
		  h=L_phi
		  alfa=-c0/dot_product(d,h)
		  x=x+alfa*d
		  g=g+alfa*h
		  c1=dot_product(g,g)
		  beta=c1/c0; c0=c1
		  d=g+beta*d
		endif
	enddo
	b=x
	if(con>eps) then
		write(6,*)'Poisson is not converged!'
	endif

	deallocate(g,d,x,h)
END SUBROUTINE CoGr


!################################################!
!		Not Used			 !
!################################################!
SUBROUTINE bc_multipole
  ! Boundary conditions determined by multipole expansion 
  integer :: k1,k2,k3,i,l,lm
  real*8  :: xx,yy,zz,x,y,z,r,rho_lm((L_max+1)**2)

  rho_lm=0.d0
  do L=0,L_max
    do lm=L**2+1,(L+1)**2
      do i=1,N_L_points
         x=grid_point(1,i); y=grid_point(2,i); z=grid_point(3,i)
         r=sqrt(x*x+y*y+z*z)+small ; xx=x/r ; yy=y/r ; zz=z/r
         rho_lm(lm)=rho_lm(lm)+r**L*Ylm(xx,yy,zz,lm)*rho(i)*grid_volume
       enddo
    enddo
  enddo
  
  ! Boundary condition in the x direction
  do k2=1,N_L(2)
    do k3=1,N_L(3)
      y=grid_point(2,Lattice_inv(1,k2,k3))
      z=grid_point(3,Lattice_inv(1,k2,k3))
      x=grid_point(1,Lattice_inv(1,k2,k3))-grid_step(1)
      V_X0(k2,k3)=mp_pot(x,y,z,rho_lm)
      x=grid_point(1,Lattice_inv(N_L(1),k2,k3))+grid_step(1)
      V_XN(k2,k3)=mp_pot(x,y,z,rho_lm)
    enddo
  enddo
  
  ! Boundary condition in the y direction
  do k1=1,N_L(1)
    do k3=1,N_L(3)
      x=grid_point(1,Lattice_inv(k1,1,k3))
      z=grid_point(3,Lattice_inv(k1,1,k3))
      y=grid_point(2,Lattice_inv(k1,1,k3))-grid_step(2)
      V_Y0(k1,k3)=mp_pot(x,y,z,rho_lm)
      y=grid_point(2,Lattice_inv(k1,N_L(2),k3))+grid_step(2)
      V_YN(k1,k3)=mp_pot(x,y,z,rho_lm)
    enddo
  enddo
  
  ! Boundary condition in the z direction
  do k1=1,N_L(1)
    do k2=1,N_L(2)
      x=grid_point(1,Lattice_inv(k1,k2,1))
      y=grid_point(2,Lattice_inv(k1,k2,1))
      z=grid_point(3,Lattice_inv(k1,k2,1))-grid_step(3)
      V_Z0(k1,k2)=mp_pot(x,y,z,rho_lm)
      z=grid_point(3,Lattice_inv(k1,k2,N_L(3)))+grid_step(3)
      V_ZN(k1,k2)=mp_pot(x,y,z,rho_lm)
    enddo
  enddo
END SUBROUTINE BC_multipole


!################################################!
!		Not Used			 !
!################################################!
FUNCTION mp_pot(x,y,z,rho_lm)
	integer :: L,lm
	real*8  :: mp_pot,xx,yy,zz,r,sump,x,y,z,rho_lm((L_max+1)**2)
  r=sqrt(x*x+y*y+z*z)+small ; xx=x/r ; yy=y/r ; zz=z/r
  sump=0.d0
  do L=0,L_max
    do lm=L**2+1,(L+1)**2
      sump=sump+Ylm(xx,yy,zz,lm)/r**(L+1)*rho_lm(lm)
    enddo
  enddo
  mp_pot= E2*sump
END FUNCTION mp_pot

!################################################!
!		Not Used			 !
!################################################!
! Ylm
! Description: explicity gives the spherical harmonic functions (real only)
! Parameters:	lm - unique number given to each lm pair by the formula l(l+1)+(m+1)
!			the pairs are show next to each Ylm for reference
! post: returns the appropriate Ylm function, but it is only the real part and without the
! 	constant
FUNCTION Ylm(x,y,z,lm)
	real*8,intent(IN) :: x,y,z
	integer,intent(IN)          :: lm
	real*8            :: Ylm,r
	! Spherical harmonics
	! It should be multiplied by sqrt((2*l+1)/4*pi)/r**l
	r=x**2+y**2+z**2
	select case( lm )
		case(1)  ; Ylm=1.d0                                     ! lm=1  (0  0)
		case(2)  ; Ylm=y                                        ! lm=2  (1 -1)
		case(3)  ; Ylm=z                                        ! lm=3  (1  0)
		case(4)  ; Ylm=x                                        ! lm=4  (1  1)
		case(5)  ; Ylm=sqrt(3.d0)*x*y                           ! lm=5  (2 -2)
		case(6)  ; Ylm=sqrt(3.d0)*y*z                           ! lm=6  (2 -1)
		case(7)  ; Ylm=(2*z*z-x*x-y*y)/2.d0                     ! lm=7  (2  0)
		case(8)  ; Ylm=sqrt(3.d0)*x*z                           ! lm=8  (2  1)
		case(9)  ; Ylm=sqrt(3.d0/4.d0)*(x*x-y*y)                ! lm=9  (2  2)
		case(10) ; Ylm=sqrt(5.d0/8.d0)*y*(3*x*x-y*y)            ! lm=10 (3 -3)
		case(11) ; Ylm=sqrt(15.d0)*x*y*z                        ! lm=11 (3 -2)
		case(12) ; Ylm=sqrt(3.d0/8.d0)*y*(4*z*z-x*x-y*y)        ! lm=12 (3 -1)
		case(13) ; Ylm=z*(2*z*z-3*x*x-3*y*y)/2.d0               ! lm=13 (3  0)
		case(14) ; Ylm=sqrt(3.d0/8.d0)*x*(4*z*z-x*x-y*y)        ! lm=14 (3  1)
		case(15) ; Ylm=sqrt(15.d0/4.d0)*z*(x*x-y*y)             ! lm=15 (3  2)
		case(16) ; Ylm=sqrt(5.d0/8.d0)*x*(x*x-3*y*y)            ! lm=16 (3  3)
		case(17) ; Ylm=sqrt(35.d0)/2.d0*x*y*(x**2-y**2)         ! lm=17 (4 -4)
		case(18) ; Ylm=sqrt(35.d0/8.d0)*y*z*(3*x**2-y**2)       ! lm=18 (4 -3)
		case(19) ; Ylm=sqrt(5.d0)/2.d0*x*y*(7*z**2-r**2)        ! lm=19 (4 -2)
		case(20) ; Ylm=sqrt(5.d0/8.d0)*y*(7*z**3-3*z*r**2)      ! lm=20 (4 -1)
		case(21) ; Ylm=(35*z**4-30*z**2*r**2+3.d0*r**2)/8.d0    ! lm=21 (4  0)
		case(22) ; Ylm=sqrt(5.d0/8.d0)*x*(7*z**3-3*z*r**2)      ! lm=22 (4  1)
		case(23) ; Ylm=sqrt(5.d0)/4.d0*(7*z**2-r**2)*(x**2-y**2)! lm=23 (4  2)
		case(24) ; Ylm=sqrt(35.d0/8.d0)*z*x*(x**2-3*y**2)       ! lm=24 (4  3)
		case(25) ; Ylm=sqrt(35.d0)/8.d0*(x**4+y**4-6*x**2*y**2) ! lm=25 (4  4)
	end select
END FUNCTION Ylm


!###############################################!
!		Not Used			!
!###############################################!
SUBROUTINE Hamiltonian_density
  ! Calculate the density dependent part of the Hamiltonian
  rho=density
  
  call bc_multipole
  call fd_p
  vh=-v_pot
  call Exchange_Correlation
END SUBROUTINE Hamiltonian_density


!################################################!
!		Not Used			 !
!################################################!
SUBROUTINE Exchange_Correlation
  ! Calculate the XC potential and energy
  USE XC
  integer :: i
  real*8  :: su,eps_c,eps_x,rho_up,rho_down,v_x_up,v_x_down,v_c_up,v_c_down

  su=0.d0
  do i=1,N_L_points
    rho_up=density_up(i)+1.d-20
    rho_down=density_dw(i)+1.d-20
    call xc_pot(rho_up,rho_down,v_x_up,v_x_down,v_c_up,v_c_down,eps_x,eps_c)
    V_Exchange_up(i)=2*Ry*(v_x_up+v_c_up)
    V_Exchange_dw(i)=2*Ry*(v_x_down+v_c_down)
    su=su+Density(i)*(eps_x+eps_c)-Density_up(i)*V_exchange_up(i)-Density_dw(i)*V_exchange_dw(i)
  enddo
  E_exchange=su*grid_volume
END SUBROUTINE Exchange_correlation


!################################################!
!		Not Used			 !
!################################################!
SUBROUTINE Hartree_Energy
  !     calculate the Hartree-Energy
  E_Hartree=-0.5d0*dot_product(density,VH)*grid_volume
end subroutine Hartree_Energy



  function basis_func(r,nu,l,m)
USE ANGLIB
    integer            :: l,lm,m
    double precision   :: r(3),nu
    double precision   :: x,w1,w2,sh
    complex*16         :: spc,basis_func
    lm=l**2+(l+1)+m
    x=sqrt(r(1)**2+r(2)**2+r(3)**2)
    w1=dexp(0.5d0*(x2n(l)-df(2*l+1)))
    w2=(2.d0*nu/pi_)**0.75d0*w1*(sqrt(2.d0*nu))**l*sqrt(4.d0*pi_)

    spc=spherical_harmonics(r,l,m)
    
!    ww=basis_func_radial(r,nu,l)
!    basis_func=ww*spc

    basis_func=exp(-nu*x**2)*w2*x**l*spc
  end function basis_func


  
  function basis_func_radial(r,nu,l)
    integer            :: l
    double precision   :: r(3),nu
    double precision   :: x,w1,w2,sh
    double precision   :: basis_func_radial

    x=sqrt(r(1)**2+r(2)**2+r(3)**2)
    w1=dexp(0.5d0*(x2n(l)-df(2*l+1)))
    w2=(2.d0*nu/pi_)**0.75d0*w1*(sqrt(2.d0*nu))**l*sqrt(4.d0*pi_)

    basis_func_radial=exp(-nu*x**2)*w2*x**l

  end function basis_func_radial



subroutine hartree_matrix_elements_init
USE MATLIB
  implicit none 
  integer             :: i,j,k,l
  double precision    :: ra(3),rb(3),rc(3),rd(3),a,b,c,d,su,suh
  integer             :: la,ma,lb,mb,lc,mc,ld,md
  complex*16          :: w

  do i=1,basis_dim
          write(6,*)i
    la=basis_l(i)
    ma=basis_m(i)
    a=basis_nu(i)
    ra(:)=basis_r(:,i)    
    do j=1,basis_dim
      lb=basis_l(j)
      mb=basis_m(j)
      b=basis_nu(j)
      rb(:)=basis_r(:,j)    
      do k=i,basis_dim
        lc=basis_l(k)
        mc=basis_m(k)
        c=basis_nu(k)
        rc(:)=basis_r(:,k)    
        do l=1,basis_dim
          ld=basis_l(l)
          md=basis_m(l)
          d=basis_nu(l)
          rd(:)=basis_r(:,l)    
          call erepmat(la,ma,lb,mb,lc,mc,ld,md,a,b,c,d,ra,rb,rc,rd,w)
          el_rep(i,j,k,l)=w
        end do
      end do
    end do
  end do

end  subroutine hartree_matrix_elements_init

!################################################!
!		Not Used			 !
!################################################!
subroutine hartree_matrix_elements
  implicit none 
  integer             :: i,j,k,l
  double precision    :: ra(3),rb(3),rc(3),rd(3),a,b,c,d,w,su,suh
  integer             :: la,ma,lb,mb,lc,mc,ld,md

  suh=0.d0
  do i=1,basis_dim
    do j=1,basis_dim
      su=0.d0
      do k=1,basis_dim
        do l=1,basis_dim
          w=el_rep(l,k,j,i)
          su=su+dm(k,l)*w
        end do
      end do
      mat_hartree(i,j)=su
      suh=suh+dm(i,j)*su
    end do
  end do
  E_hartree=-0.5d0*suh

end  subroutine hartree_matrix_elements

!################################################!
!		Not Used			 !
!################################################!
subroutine wavefn
  implicit none
  double precision   :: r(3),s
  integer            :: i,j,k
  complex*16         :: b
  
  r=0.d0
  do k=1,N_orbitals
    do i=1,1000
      r(1)=i*0.01          !watch out   another parameter with no comment
      s=0.d0
      do j=1,basis_dim
        b=basis_func(r,basis_nu(j),basis_l(j),basis_m(j))
        s=s+basis_coef(j,k)*b
      end do
      write(50,*)r(1),s
    end do
    write(50,*)
  end do
  
end  subroutine wavefn



! matrix_elements_init
! Description: 
! Parameters:	
subroutine matrix_elements_init
USE MATLIB
  implicit none 
  integer             :: i,j,n
  double precision    :: ra(3),rb(3),rc(3),a,b
  complex*16          :: ovmr,kimr,numr,ten(3)
  integer             :: la,ma,lb,mb

  Z_charge=Z_charge_inp    !
!!!!?????
  rc=0.d0
  n=basis_dim
  do i=1,n
    do j=1,n
      a=basis_nu(i)
      b=basis_nu(j)
      la=basis_l(i)
      lb=basis_l(j)
      ma=basis_m(i)
      mb=basis_m(j)
      ra(:)=basis_r(:,i)
      rb(:)=basis_r(:,j)



!      call Coulomb_potential_m(a,ra,la,b,rb,lb,rc)
      call gauss_kinetic_m(a,ra,la,b,rb,lb)
!      h(i,j)=h2m*stm(0,0)-cpm(0,0)*Z_charge
!      o(i,j)=som(0,0)
      call okmat(la,ma,lb,mb,a,b,ra,rb,ovmr,kimr)
      call nucmat(la,ma,lb,mb,a,b,ra,rb,rc,numr)
      call tenmat(la,ma,lb,mb,a,b,ra,rb,1,ten)
      k_m(1,i,j)=sqrt(2.d0*pi/3.d0)*(ten(1)-ten(3))
      k_m(2,i,j)=sqrt(2.d0*pi/3.d0)*(ten(1)+ten(3))*zi
      k_m(3,i,j)=sqrt(4.d0*pi/3.d0)*ten(2)
!      write(78,*)i,j,ma,mb
!      write(78,*)ten
!      write(77,*)i,j,ma,mb
!      write(77,*)k_m(1,i,j),nam(1,ma,mb)
!      write(77,*)k_m(2,i,j),nam(2,ma,mb)
!      write(77,*)k_m(3,i,j),nam(3,ma,mb)
!      write(77,*)som(ma,mb)
!      k_m(1,i,j)=nam(1,ma,mb)
!      k_m(2,i,j)=nam(2,ma,mb)
!      k_m(3,i,j)=nam(3,ma,mb)
!      write(177,*)i,j,ma,mb
!      write(177,*)kimr,ovmr
!      write(177,*)stm(ma,mb),som(ma,mb)

      h_m(i,j)=h2m*kimr-numr*Z_charge
      o_m(i,j)=ovmr
      l_m(i,j)=lim(ma,mb)
    end do
  end do
end subroutine matrix_elements_init


subroutine init_basis
  USE MATLIB
  USE ANGLIB
  implicit none
  integer                              :: k,k1,k2,i,j,l,kk,it,atom,m


  allocate(basis_l(basis_dim),basis_m(basis_dim),basis_nu(basis_dim),basis_coef(basis_dim,N_orbitals),e_sp(N_orbitals))
  allocate(basis_r(3,basis_dim))
  allocate(mat_vex_up(basis_dim,basis_dim),mat_vex_dw(basis_dim,basis_dim),mat_hartree(basis_dim,basis_dim), &
&   dm_a(basis_dim,basis_dim),dm_o(basis_dim,basis_dim),dm_n(basis_dim,basis_dim), &
&   dm(basis_dim,basis_dim),mat_hartree_num(basis_dim,basis_dim),el_rep(basis_dim,basis_dim,basis_dim,basis_dim), &
&   h_m(basis_dim,basis_dim),o_m(basis_dim,basis_dim),k_m(3,basis_dim,basis_dim),l_m(basis_dim,basis_dim))





end subroutine init_basis



!################################################!
!		Not Used			 !
!################################################!
  subroutine set_atom_pairs
    integer :: i,j,it,jt
    real(8), dimension(3) :: p0

    allocate(atom_pairs(N_total_atom,N_total_atom))

    atom_pairs = 0
    do i=1,N_total_atom
       it=atom_index(i)
       p0 = Position(:,i)
       do j=1,N_total_atom
          jt=atom_index(j)
          if (sqrt(sum(p0 - position(:,j))**2) < wf_rad(it)+wf_rad(jt)) then
             atom_pairs(i,j) = 1
          end if
       end do
    end do

  end subroutine set_atom_pairs

  
!################################################!
!		Not Used			 !
!################################################!
  subroutine rad_bas
    USE ANGLIB
    implicit none
    integer                       :: max_p,atom,k,m,i,j,l,it,n
    double precision              :: s,r(3),r2
    complex*16                    :: sph
    double precision, allocatable :: b(:)



    allocate(b(N_gauss))
    allocate(n_int_points(N_total_atom),wf_rad(N_species))

    wf_rad=6.d0



    max_p=0
    n_int_points = 0
    do atom=1,N_total_atom
       k=0
       do i=1,N_L_points
          r(:)=grid_point(:,i)-Position(:,atom)+1.d-12
          r2=sqrt(r(1)**2+r(2)**2+r(3)**2)
          if(r2.lt.wf_rad(atom_index(atom))) k=k+1
       end do
       if(k.gt.max_p) max_p=k
       n_int_points(atom)=k
    end do

    allocate(psi_basis(max_p,basis_dim),int_points_atom(max_p,N_total_atom))

    psi_basis = (0.0d0,0.d0)
    int_points_atom = 0

    do atom=1,N_total_atom
       k=0
       do i=1,N_L_points
          r(:)=grid_point(:,i)-Position(:,atom)+1.d-12
          r2=sqrt(r(1)**2+r(2)**2+r(3)**2)
          if(r2.lt.wf_rad(atom_index(atom))) then
             k=k+1
             int_points_atom(k,atom)=i
          endif
       end do
    end do


   do atom=1,N_total_atom
     it=atom_index(atom)
     do l=0,l_max_basis(it)
       write(6,*)atom,l,n_int_points(atom)
       do k=1,n_int_points(atom)
         i=int_points_atom(k,atom)
         r(:)=grid_point(:,i)-Position(:,atom)+1.d-12
         r2=sqrt(r(1)**2+r(2)**2+r(3)**2)
         do j=1,N_gauss
           b(j)=basis_func_radial(r,basis_nu(j),l)
         end do
         do j=1,N_ao_l(l,it)
           s=0.d0
           do n=1,N_Gauss
             s=s+ao_basis_c(n,j,l,atom)*b(n)
           end do
           do m=-l,l
             sph=spherical_harmonics(r,l,m)             
             Psi_basis(k,basis_ind(m,j,l,atom))=s*sph
           enddo
         enddo
       enddo
     enddo
   end do

   do i=1,basis_dim
     write(6,*)i,sum(Conjg(Psi_basis(:,i))*Psi_basis(:,i))*grid_volume
   end do
   
end subroutine rad_bas


!################################################!
!		Not Used			 !
!################################################!
subroutine num_po
  integer                                      :: i,j,k,ii,i1,i2,j1,j2
  double precision                             :: su
  complex*16, allocatable                      :: H(:)

  allocate(H(N_L_points))

  do i=1,basis_dim
    i1=basis_atom(1,i)
    i2=basis_atom(2,i)
    H=0.d0
    do k=1,n_int_points(i2)
      ii=int_points_atom(k,i2)
      H(ii)=psi_basis(k,i)
    end do
    do j=1,i
      j1=basis_atom(1,j)
      j2=basis_atom(2,j)
!      if(atom_pairs(i2,j2).eq.1) then
        su=0.d0
        do k=1,n_int_points(j2)
          ii=int_points_atom(k,j2)
          su=su+H(ii)*psi_basis(k,j)
        end do
!      endif
    end do
  end do

  deallocate(H)
  
end subroutine num_po



END MODULE GAUSS_DFT_ATOM










