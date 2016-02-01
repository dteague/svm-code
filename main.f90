  
  ! Two electrons and a nucleus using Shifted Correlated Gaussian basis
  ! Code version 1.1   
  ! Authors  Hanif Ahmed, Samuel Lane, Dylan Teague, Jorge Salas and Kalman Varga
  
! Contains the subroutines and functions:

! init_main_1
! calculate_me_He
! calculate_me_Hydrogen
! calculate_wf
! file_name
! 
! time_dep_H_0
! time_dep_H_t
  
  
use input_mod
implicit none
   integer, parameter	:: eigenvaluesf=16,hamiltonianf=17, overlapf=19, ham_timef=23
   integer, parameter	:: c_groundstatef=29,wavef=19,densityf=50,normf=900

  
  call var_input
  call init_main_1
  call time_dep_H_t
!   if(is_helium) then 
!    
!   end if
!   
!   if (.not. is_helium) then
!       
!       call calculate_me_Hydrogen
!   end if
  

end


! init_main_1
! Description; sets up the needed elements for the program to run 
! Parameters:	lao - linear atomic orbitals?
!		basis_dim - dimension of the wavefunction bases
! post: creates the basis, lattice, mattrices etc for the program by reading off bases
! 	values from the file "basis.dat".  Need to run basis.f90 to make this file 

subroutine init_main_1
USE GAUSS_DFT_ATOM
USE LINALG
USE MATLIB
USE GAUSS_BASIS
implicit none

  integer	:: i,j,lao
  logical	:: basis_exist
  
  call gauss_init

  INQUIRE(FILE="basis.dat", EXIST=basis_exist)
  
  if(.not. basis_exist) call exit(1) !exit program is basis.dat doesn't exist
  
  open(1,file='basis.dat')
  read(1,*)basis_dim

  call init_basis
! call init_lattice
  call init_matlib
  lao=l_ao_max
  call precal_gaunt(lao)
  call precal_gaunt_3(lao)

  do i=1,basis_dim
    read(1,*)basis_l(i),basis_m(i)
    read(1,*)basis_nu(i)
    read(1,*)(basis_r(j,i),j=1,3)
  end do

  call hartree_matrix_elements_init
  call matrix_elements_init

end subroutine init_main_1


subroutine file_name(prefix, file, type_file)
use input_mod
use gauss_dft_atom
implicit none
  
  character (len=*)	:: file, prefix, type_file
  character*20 		:: counter, space_string, nu_string, dim_string
  integer		:: ending
  logical		:: is_data

  ending=len_trim(prefix)
  
  select case(type_file)
    case("time_step_data")
      write(counter,'(i5.5)') time_step
      write(file,'(a,a,a)') prefix(1:ending),counter(1:5),'.dat'
      
    case("matrices")
      if(basisf_exist) then
	write(dim_string,"(I3)") basis_dim
	write(nu_string,"(F5.2)") nu_val
	write(space_string,"(F5.2)") space
      	write(file,"(4a,I1,5a)") prefix(1:ending),"_n", trim(adjustl(dim_string)), "_ll", ll_max, &
&	   "_nu",trim(adjustl(nu_string)), "_x",trim(adjustl(space_string)), ".dat"
      else
	write(file,'(a,a,a)') prefix(1:ending),"_n", basis_dim,'.dat'
	write(*,*) "warning, generic name used"
      end if
  end select
    
end subroutine file_name


!#################################################################!
!			    ATOM SETUP				  !
!#################################################################!



! calculate_me_He 
! Description:	calculates the hamiltonian, overlap, as well as eigenvectors for the basis used
! 	and then diagonalizes the matrices to output the eigenvalues
! parameters:	h(:,:) - matrix of the hamiltonian
!		o(:,:) - matrix of the overlap
!		v(:,:) - matrix of eigenvectors, 1st being the eigenvector, 2nd being the ith value
!		e(:) - eigenvalues
!		i,j,k,l, ij, kl - dummy variables for loops in matrix multiplication
! post: gives the energy eigenvalues for the hamiltonian

subroutine calculate_me_He
USE GAUSS_DFT_ATOM
USE LINALG
USE MATLIB
USE GAUSS_BASIS
implicit none

  double precision, allocatable :: h(:,:),o(:,:),v(:,:),e(:)
  integer                       :: i,j,n,k,l,ij,kl
  double precision              :: a,b,ra(3),rb(3),rc(3),numr,ovmr,kimr
  integer                       :: la,ma,lb,mb
  integer, parameter		:: eigenvaluesf=16,hamiltonianf=17, overlapf=19, ham_timef=23
  integer, parameter		:: c_groundstatef=29,wavef=19,densityf=50,normf=900
  double precision, parameter	:: helium_gs = -2.9025
  character*255			:: hamn, overn, cn




  n=(basis_dim*(basis_dim+1))/2	!n is equal to this becuase we assume two fermions so need to include antisymmetric terms
  write(6,*)'dim',n		! ie psi = a(1)b(2)-a(2)b(1)
  
  allocate(h(n,n),o(n,n),e(n),v(n,n))
  
  call file_name("h",hamn,"matrices") 
  call file_name("o",overn,"matrices")
  call file_name("c",cn,"matrices")
  
  open(eigenvaluesf,file="eigenenergies.dat")
  open(hamiltonianf,file=trim(hamn))
  open(overlapf,file=trim(overn))
  open(c_groundstatef,file=trim(cn))


  ij=0
  do i=1,basis_dim
    do j=i,basis_dim
      ij=ij+1
      kl=0
      
      do k=1,basis_dim
        do l=k,basis_dim
          kl=kl+1
          
  ! calculation of hamiltonian by <ij|H(|kl>+|lk>) = <ij|H1+H2|kl>+<ij|H1+H2|lk>
          h(ij,kl)=h_m(i,k)*o_m(j,l)+o_m(i,k)*h_m(j,l)+h_m(i,l)*o_m(j,k)+ &
&           o_m(i,l)*h_m(j,k)+el_rep(i,k,j,l)+el_rep(i,l,j,k)

  ! overlap matrix found by <ij|(|kl>+|lk>) = <ij|kl>+<ij|lk>
          o(ij,kl)=o_m(i,k)*o_m(j,l)+o_m(i,l)*o_m(j,k) 
        end do
      end do
    end do
  end do

  write(hamiltonianf,"(ES20.12)") h
  write(overlapf,"(ES20.12)") o
   
  ! Diagonalize
  call diag1(h,o,n,e,v)
  
  call in_range(helium_gs, e(1))
  
  
  ! Output the eigenvalues (only first 10 since bad approx after first few)
  write(eigenvaluesf,"(ES20.12)") e

!   ground state wave functions 
  write(c_groundstatef,"(ES20.12)") (v(i,1),i=1,n) 
  
  close(hamiltonianf)
  close(eigenvaluesf)
  close(overlapf)
  close(c_groundstatef)
  
deallocate(h,o,e,v)

end subroutine calculate_me_He



! calculate_me_Hydrogen
! Description:	calculates the hamiltonian, overlap, as well as eigenvectors for the basis used
! 	and then diagonalizes the matrices to output the eigenvalues for hydrogen. 
! parameters:	h(:,:) - matrix of the hamiltonian
!		o(:,:) - matrix of the overlap
!		v(:,:) - matrix of eigenvectors, 1st being the eigenvector, 2nd being the ith value
!		e(:) - eigenvalues
! post: gives the energy eigenvalues for the hamiltonian and writes the hamiltonian, overlap, and ground
!	state into dat files

subroutine calculate_me_Hydrogen
USE GAUSS_DFT_ATOM
USE LINALG
USE MATLIB
USE GAUSS_BASIS
use input_mod
implicit none

  double precision, allocatable :: h(:,:),o(:,:),v(:,:),e(:), k(:,:)
  integer                       :: i,j,n
  double precision              :: a,b,ra(3),rb(3),rc(3),numr,ovmr,kimr, k_vec(3), k2
  integer                       :: la,ma,lb,mb
  integer, parameter		:: eigenvaluesf=16,hamiltonianf=17, overlapf=19, ham_timef=23
  integer, parameter		:: c_groundstatef=29,wavef=19,densityf=50,normf=900
  double precision, parameter	:: hydrogen_gs = -0.5d0
  character*255			:: hamn, overn, cn


  n=basis_dim	
  write(6,*)'dim',n		! ie psi = a(1)b(2)-a(2)b(1)
  allocate(h(n,n),o(n,n),e(n),v(n,n),k(n,n))
  
  call file_name("h",hamn,"matrices") 
  call file_name("o",overn,"matrices")
  call file_name("c",cn,"matrices")
  
  open(eigenvaluesf,file="eigenenergies.dat")
  open(hamiltonianf,file=trim(hamn))
  open(overlapf,file=trim(overn))
  open(c_groundstatef,file=trim(cn))

 
  h=h_m
  o=o_m
  
  
  
!   k_vec=k_vec_inp
!   k2=k_vec(1)**2+k_vec(2)**2+k_vec(3)**2
! 
!   do i=1,basis_dim
!     do j=1,basis_dim
!       k(i,j)=k_vec(1)*k_m(1,i,j)+k_vec(2)*k_m(2,i,j)+k_vec(3)*k_m(3,i,j) !adjusted matrix that uses k
!     end do
!   end do
!   
!   h=h+h2m*(-2.d0*zi*k+k2*o)
  
  write(hamiltonianf,"(ES20.12)") h
  write(overlapf,"(ES20.12)") o
  
  ! Diagonalize
  call diag1(h,o,n,e,v)

  call in_range(hydrogen_gs, e(1))
  
  
 ! Output the eigenvalues (only first 10 since bad approx after first few)
  write(eigenvaluesf,"(ES20.12)") e

!   ground state wave functions 
  write(c_groundstatef,"(ES20.12)") (v(i,1),i=1,n) 
  
  close(hamiltonianf)
  close(eigenvaluesf)
  close(overlapf)
  close(c_groundstatef)
  
  deallocate(h,o,e,v)

end subroutine calculate_me_Hydrogen


! in_range
! Description: test to see if eigenenergies calculated are within a correct bond or not
! Parameters:	gs_energy - ground state energy for the atom
!				calc - calculated energy from the basis
!				error_allowed - amount of leway the program gives for the calculated energy
! post: sends error and stops program is energies weren't calculated, were lower than gs or higher
!	than the gs + error.  Otherwise, the program outputs the groundstate energy and then returns

subroutine in_range(gs_energy, calc)
use input_mod
implicit none
  double precision		:: gs_energy, calc
  double precision,parameter	:: error_allowed=0.5d0
  
  INQUIRE(FILE="fort.2", EXIST=file_exists)
  
  if(file_exists) then
    write(*,*) "!!!!!ERROR: MATRIX UNDIAGONALIZABLE, POSSIBLE LINEAR DEPENDENCE EXISTS IN"
    write(*,*) "BASIS, REMOVE FORT.2!!!!!"
    call system('rm fort.2')
    call system('rm fort.9')
    call exit(1)
  end if
  
  if(calc<gs_energy) then
    write(*,*) "!!!!!PROGRAM ERROR: ENERGY IS LOWER THAN ATOM'S GROUND STATE, VARIATIONAL"
    write(*,*) "METHOD FAILED, CHECK CODE!!!!!"
    write(*,*)
    write(*,*) "gs_energy = ", calc
    call exit(1)
  else if(calc>(gs_energy+error_allowed)) then
    write(*,*) "!!!!!WARNING: GROUND STATE HIGHER THAN ALLOWED AMOUNT, CONSIDER A BETTER"
    write(*,*) "BASIS!!!!!"
    write(*,*)
    write(*,*) "gs_energy = ", calc
    return
  else 
    write(*,*) "gs_energy = ", calc
    return
  end if
  
end subroutine in_range

!###############################################################!
!			POTENTIAL SETUP				!
!###############################################################!



subroutine laser_me_setup
USE GAUSS_DFT_ATOM
USE LINALG
USE MATLIB
USE GAUSS_BASIS
use input_mod


  double precision, allocatable :: hef(:,:)
  integer                       :: i,j,k,n,nt,m,l,ij,kl
  integer, parameter	:: eigenvaluesf=16,hamiltonianf=17, overlapf=19, ham_timef=23
  integer, parameter	:: c_groundstatef=29,wavef=19,densityf=50,normf=900
  character*255		:: htn
  
  call file_name("ht",htn, "matrices")
  
  open(ham_timef,file=trim(htn))
  
  if(is_helium) then
    n=(basis_dim*(basis_dim+1))/2 
    allocate(hef(n,n))
!
! setup matrices
!

    ij=0
    do i=1,basis_dim
      do j=i,basis_dim
	    ij=ij+1
	    kl=0
	    do k=1,basis_dim
	      do l=k,basis_dim
	        kl=kl+1
	        hef(ij,kl)=l_m(i,k)*o_m(j,l)+o_m(i,k)*l_m(j,l)+l_m(i,l)*o_m(j,k)+ &
&              o_m(i,l)*l_m(j,k)
	      end do
	    end do
      end do
    end do

    write(ham_timef,"(ES20.12)") hef  
  
  else if(.not. is_helium) then
    n=basis_dim
    allocate(hef(n,n))
    do i=1,n
      do j=1,n
	hef(i,j)=l_m(i,j)
      end do
    end do
    
    write(ham_timef,"(ES20.12)") hef
    
  end if
    
  close(ham_timef)

end subroutine laser_me_setup



! subroutine singleparticle !!!!!!!!!!!!!!!!!!!NAME!!!!!!!!!!!!!!!!!!!!!!!!!
! use input_mod 

! double precision, allocatable :: h(:,:),o(:,:),k(:,:),hc(:,:)
!  integer                       :: i,j,n,k,l,ij,kl
!  double precision              :: a,b,ra(3),rb(3),rc(3),numr,ovmr,kimr, k(3), k_vec, temp
!  integer                       :: la,ma,lb,mb
!  integer, parameter		:: eigenvaluesf=16,hamiltonianf=17, overlapf=19, ham_timef=23
!  integer, parameter		:: c_groundstatef=29,wavef=19,densityf=50,normf=900
!  double precision, parameter	:: hydrogen_gs =0.5d0
!
!   open(hamiltonianf, file="hamiltonain.dat")
!   allocate(h(n,n),o(n,n),k(n,n),hc(n,n))
! 
!   k_vec=k_vec_inp
!   k2=k_vec(1)**2+k_vec(2)**2+k_vec(3)**2
! 
!   do i=1,basis_dim
!     do j=1,basis_dim
!       h(i,j)=h_m(i,j)
!       k(i,j)=k_vec(1)*k_m(1,i,j)+k_vec(2)*k_m(2,i,j)+k_vec(3)*k_m(3,i,j) !adjusted matrix that uses k
!       o(i,j)=o_m(i,j)
!     end do
!   end do
! 
!   hc=h+h2m*(-2.d0*zi*k+k2*o)
!
!
! 
! 
! end subroutine singleparticle



!###############################################################!
!			TIME PROPAGATION			!
!###############################################################!


subroutine time_dep_H_t

USE GAUSS_DFT_ATOM
USE LINALG
USE MATLIB
USE GAUSS_BASIS
use input_mod
implicit none

  double precision, allocatable :: h(:,:),o(:,:),ht(:,:)
  complex*16, allocatable       :: hc(:,:),am(:,:),zo(:,:),bm(:,:),cm(:,:),c(:), xi(:),H_xi(:)
  integer                       :: i,j,k,n,nt,m,l,ij,kl, p,q
  double precision              :: dt
  complex*16                    :: ef, prob, norm1
  integer, parameter 		:: N_taylor = 4
  complex*16                    :: tc(0:N_taylor)
  integer, parameter		:: eigenvaluesf=16,hamiltonianf=17, overlapf=19, ham_timef=23
  integer, parameter		:: c_groundstatef=29,wavef=19,densityf=50,normf=900
  character*255			:: hamn,overn,cn, htn
 
  call need_files(.true.)

  if(ntsteps_inp==0) then
    call exit(1)
  end if
  
  call file_name("h",hamn,"matrices") 
  call file_name("o",overn,"matrices")
  call file_name("c",cn,"matrices")
  call file_name("ht",htn, "matrices")
  
  
  open(hamiltonianf,file=trim(hamn))
  open(overlapf,file=trim(overn))
  open(c_groundstatef,file=trim(cn))
  open(ham_timef,file=trim(htn))
  open(55, file="groundstate.dat")
  
   if(is_helium) then
    n=(basis_dim*(basis_dim+1))/2
  else if(.not. is_helium) then
    n=basis_dim
  else
    write(*,*) "!!!!!ERROR: PROGRAM DOESN'T KNOW WHICH ATOM IS NEEDED, EDIT INPUT.INP!!!!!"
    call exit(1)
  end if
  
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

  !call calculate_wf(c,n,o)
  
  if(taylor) then 
  
    do i=1,N_taylor
      tc(i)=(-zi*dt)**i
      do m=1,i
	tc(i)=tc(i)/m
      end do
    end do

    zo=o
  
    call inv(zo,n,bm)!!!!!!!!!!!!!!!


    do i=1,nt
      call laser_packet(ef, i*dt)
      hc=h+ef*ht
      cm=matmul(bm,hc)
      xi=c
      do m=1,N_taylor
	H_xi=matmul(cm,xi)
	xi=H_xi
	c=c+tc(m)*xi
      end do
      if(mod(i,T_of_output)==0) then
        call calculate_wf(c,n,o)
        write(*,*) ef
      end if
    end do
    
    
!     do i=1,n    
!       write(6,*)c(i)
!     end do

! Crank-Nicholson propagator
  else if(.NOT. taylor) then
    do i=1,nt
      call laser_packet(ef,i*dt)
      hc=h+ef*ht
!   Crank-Nicholson matrix calculated
      am=o-0.5d0*zi*hc*dt
      zo=o+0.5d0*zi*hc*dt
  
      call inv_c(zo,n,bm)
      cm=matmul(bm,am)

      c=matmul(cm,c)  
    
      norm1=0
      prob=0
      do q=1,n
        do p=1,n
          prob=prob+conjg(c(q))*c_gs(p)*o(q,p)
          norm1=norm1+conjg(c(q))*c(p)*o(q,p)
        end do
      end do
      write(55,*) real(conjg(prob)*(prob))
!  frequency of output
      if(mod(i,T_of_output)==0) then
        write(*,*) ef
        !call calculate_wf(c,n,o) 

      end if
    end do
  
!     do i=1,n
!       write(6,*)c(i)
!     end do
    
  endif
  
  
end subroutine time_dep_H_t

subroutine laser_packet(ef, t)
use input_mod
implicit none

double precision	:: ef, t, a, E_0, freq

  E_0 = sqrt(I_0)*5.329d-9
  freq = 2*3.14159*137*5.29177d-11/lambda

  ef =  E_0 * DSIN(freq*t) * exp(-2*(t-1.5*width)**2d0/(width**2))
  
end subroutine laser_packet


! time_dep
! Description: applies the time propagator (either Taylor or Crank-Nicholson) to the initial wavefunction
!	and prints out the result
! Parameters:	h(:,:) - hamiltonian matrix
!		o(:,:) - overlap matrix
!		hc(:,:) - optimized for high enegery hamiltonian that is used with the plane wave ie 
!		(\/ is the del)	    H_k = h2m(\/^2 - 2*i*k.\/ + k^2) = h +h2m*(-2*i*k.\/ + k^2)
!		am(:,:) - second matrix of the Crank-Nicholson Method (ie (1-i*H*dt/(2*hbar))
!		zo(:,:) - non-inverted first matrix of Crank-Nicholson Method (ie (1+i*H*dt/(2*hbar))
!		bm(:,:) - inverted first matrix of Crank-Nicholson or zo^-1 
!		cm(:,:) - [bm]*[am] or the Crank-Nicholson propagator
!		c(:) - vector of the wavefunction
!		k(:,:) - k vector multiplied with the gradient ???????
!		dt - incremental time step
!		nt - number of time steps
! post: prints out the wavefunction the number of time steps after applying the time propagator to the wave
! 	function each time.

subroutine time_dep_H_0 !!!!!!!!!!!!!!!!!NAME!!!!!!!!!!!!!!
!
!  time propagation using gaussians

USE GAUSS_DFT_ATOM
USE LINALG
USE MATLIB
USE GAUSS_BASIS
USE input_mod
implicit none

  double precision, allocatable :: h(:,:),o(:,:)
  complex*16, allocatable       :: hc(:,:),am(:,:),zo(:,:),bm(:,:),cm(:,:),c(:),k(:,:),xi(:),h_xi(:)
  integer                       :: i,j,n,nt,m,l, q,p
  double precision              :: k_vec(3),k2,dt,r(3)
  complex*16                    :: su,s,b, prob
  integer,parameter             :: N_taylor=4 
  complex*16                    :: tc(0:N_taylor)
  integer, parameter		:: eigenvaluesf=16,hamiltonianf=17, overlapf=19, ham_timef=23
  integer, parameter		:: c_groundstatef=29,wavef=19,densityf=50,normf=900
  character*255			:: hamn,overn,cn
 
  
  call need_files(.false.)

!   if one wants to only get eigenenergy, can put ntstep==0
  if(ntsteps_inp==0) then
    call exit(1)
  end if
  
  call file_name("h",hamn,"matrices") 
  call file_name("o",overn,"matrices")
  call file_name("c",cn,"matrices")
  
  
  open(hamiltonianf,file=trim(hamn))
  open(overlapf,file=trim(overn))
  open(c_groundstatef,file=trim(cn))
  
  if(is_helium) then
    n=(basis_dim*(basis_dim+1))/2
  else if(.not. is_helium) then
    n=basis_dim
  else
    write(*,*) "!!!!!ERROR: PROGRAM DOESN'T KNOW WHICH ATOM IS NEEDED, EDIT INPUT.INP!!!!!"
    call exit(1)
  end if
  
  
  write(6,*)'dim',n
  allocate(h(n,n),o(n,n),c_gs(n))
  allocate(hc(n,n),zo(n,n),am(n,n),cm(n,n),bm(n,n),c(n),xi(n),h_xi(n))
!
! setup matrices
!
   do i=1, n
    do j=1, n
      read(hamiltonianf,*) hc(i,j)
      read(overlapf,*) o(i,j)
    end do
  end do
  
  read(c_groundstatef,*) (c_gs(i),i=1,n)
  
  c=c_gs

  call calculate_wf(c,n,o)
  
  nt=ntsteps_inp
  dt=dt_inp


! Taylor propagator - commented since Crank-Nicholson method used instead  

  if(taylor) then 
  
    do i=1,N_taylor
      tc(i)=(-zi*dt)**i
      do m=1,i
	tc(i)=tc(i)/m
      end do
    end do

    zo=o!!!!!!!!!!!!!!!!!!!!!
  
    call inv(zo,n,bm)
    cm=matmul(bm,hc)


    do i=1,nt
      xi=c
      do m=1,N_taylor
	H_xi=matmul(cm,xi)
	xi=H_xi
	c=c+tc(m)*xi
      end do
      if(mod(i,1000)==0) call calculate_wf(c,n,o)
    end do
    
!     do i=1,n
!       write(6,*)c(i)
!     end do
   


! Crank-Nicholson propagator
  else if(.NOT. taylor) then
  
    am=o-0.5d0*zi*hc*dt
    zo=o+0.5d0*zi*hc*dt
  
    call inv_c(zo,n,bm)
  
!   Crank-Nicholson matrix calculated
    cm=matmul(bm,am)

    do i=1,nt
      c=matmul(cm,c)  
!	frequency of output
      if(mod(i,T_of_output)==0) then
        call calculate_wf(c,n,o) 
        prob=0
        do q=1,n
          do p=1,n
            prob=prob+conjg(c(q))*c_gs(p)*o(q,p)
          end do
        end do
        write(*,*) prob
      end if
    end do
  
!     do i=1,n
!       write(6,*)c(i)
!     end do
      
  endif

end subroutine time_dep_H_0




! calculate_wf
! Description: calculates the norm, wavefunction, and density of the given wavefunction along the x axis
! Parameters:	n - dimension size of the basis
! 		c(n) -  coefficients of the wavefunction (ith vector)
! 		su - normalization constant
! 		s - wavefunction in total (ie sum(c_i*psi_i))
! 		b - basis function (note, this is the ith basis function.  The index is hidden in
!			the function basis_func in which the ith cooridinates are entered
! 		r(3) - Distance from the origin
! 		x(3) - total distance or r - R_initial
! post: outputs 3 files called "wavefunction.dat", "density.dat", and "norm.dat" which have the respective 
!	information in each

subroutine calculate_wf(c,n,o)
USE GAUSS_DFT_ATOM
USE MATLIB
USE input_mod
implicit none
  integer                       :: n
  double precision		:: o(n,n)
  complex*16                    :: c(n)
  integer                       :: i,j,l,p,q
  complex*16                    :: su,s,b_j
  double precision              :: r(3),x(3)
  integer, parameter		:: eigenvaluesf=16,hamiltonianf=17, overlapf=19, ham_timef=23
  integer, parameter		:: c_groundstatef=29,wavef=19,densityf=50,normf=900
  character*255			:: waven, densityn, normn, overn
  
  
!  call file_name("wavefunction", waven, "time_step_data")
  call file_name("density", densityn, "time_step_data")
  call file_name("norm", normn, "time_step_data")
  
!  open(wavef,file=trim(waven)) 
  open(densityf,file=trim(densityn))
  open(normf,file=trim(normn))
  
  
  time_step=time_step+1
 
  
!   calculate norm
  su=(0.d0,0.d0)
  do j=1,n
    do l=1,n
      su=su+c(l)*conjg(c(j))*o(l,j)
    end do
  end do
  write(normf,*) real(su)            !this is the norm
    
  r=0.d0
  x=0.d0
!   do l=-half_lat_num,half_lat_num !!!!!!!!!!!!!!!!!!!!!!!!!LIN DEP PROBLEM!!!!!!!!!!!!!!!!!!!!!
!     r(3)=l*space_res(3)
!     do p=-half_lat_num,half_lat_num
!       r(2)=p*space_res(2)
      do q=-half_lat_num,half_lat_num
	r(1)=q*space_res(1)
	call helium_density(c,n,o_m,r,su,s)
	write(densityf,*) r(1), real(s)
	
	
	
!   calculates wavefuncion
! 	s=0.d0
! 	do j=1,basis_dim
! 	    x(:)=r(:)-basis_r(:,j)
! 	    b_j=basis_func(x,basis_nu(j),basis_l(j),basis_m(j))  
! 	  s=s+c(j)*b_j/su**(0.5d0)   ! eq (14.1) in Kalman's book
! 	end do
! 	write(wavef,*) r(1), s
! 	write(densityf,"(F6.2,F6.2,F6.2,ES21.12E4)") r(1),r(2),r(3), real(Conjg(s)*s)
! 	
!       end do
!     end do
  end do


end subroutine calculate_wf


subroutine helium_wave(s, c, r1, r2, su)
USE GAUSS_DFT_ATOM
USE MATLIB
USE input_mod
implicit none

  integer		:: i, j, k
  double precision	:: a(3), b(3), r1(3), r2(3), su
  complex*16		:: a_i, a_j, b_j, b_i, c, s
  
  k=1
  s=0
  do j=1, basis_dim
    a(:)=r1(:)-basis_r(:,j)
    a_j=basis_func(a,basis_nu(j),basis_l(j),basis_m(j)) 
    
    do i=j, basis_dim
      b(:)=r2(:)-basis_r(:,i)
      a_i=basis_func(a,basis_nu(i),basis_l(i),basis_m(i))  
      b_i=basis_func(b,basis_nu(i),basis_l(i),basis_m(i))
      b_j=basis_func(b,basis_nu(j),basis_l(j),basis_m(j))
      
      s = c(k)*(a_i*b_j + a_j*b_i)/su**(0.5d0)
      k=k+1
    end do
  end do

end subroutine

subroutine helium_density(c,n,o,r,su,den)
USE GAUSS_DFT_ATOM
USE MATLIB
USE input_mod
implicit none
  integer		:: i, j, k, l, n
  double precision	:: x1(3),x2(3),x3(3),x4(3), r(3), o(basis_dim, basis_dim)
  complex*16		:: b_i,b_j,b_k,b_l, c(n), su, const,den
  
  den = 0
  
  do i=1, basis_dim
    x1(:)=r(:)-basis_r(:,i)
    b_i=basis_func(x1,basis_nu(i),basis_l(i),basis_m(i))
    
    do j=i, basis_dim
      x2(:)=r(:)-basis_r(:,j)
      b_j=basis_func(x2,basis_nu(j),basis_l(j),basis_m(j))
      
      do k=1, basis_dim
        x3(:)=r(:)-basis_r(:,k)
        b_k=basis_func(x3,basis_nu(k),basis_l(k),basis_m(k))
        
        do l=k, basis_dim
          x4(:)=r(:)-basis_r(:,l)
          b_l=basis_func(x4,basis_nu(l),basis_l(l),basis_m(l))
          const = conjg(c((i-1)*(2*basis_dim-i)/2+j))*c((k-1)*(2*basis_dim-k)/2+l)
          den = den + const*(conjg(b_i)*(b_k*o(j,l)+b_l*o(j,k)))
        end do
      end do
    end do
  end do
  
  den = den/su


end subroutine helium_density


subroutine need_files(is_H_t)
use input_mod
implicit none

  character*255 :: hamn, overn, cn, htn
  logical 	:: ham_exist, over_exist, c_exist, time_exist, is_H_t
  
  call file_name("h",hamn,"matrices") 
  call file_name("o",overn,"matrices")
  call file_name("c",cn,"matrices")
  call file_name("ht",htn, "matrices")
  
  INQUIRE(FILE=trim(hamn), EXIST=ham_exist)
  INQUIRE(FILE=trim(overn), EXIST=over_exist)
  INQUIRE(FILE=trim(cn), EXIST=c_exist)
  INQUIRE(FILE=trim(htn), EXIST=time_exist)
  
  if(.not. (ham_exist .and. over_exist .and. c_exist)) then
    if(is_helium) then
      call calculate_me_He
    else if(.not. is_helium) then
      call calculate_me_Hydrogen
    else
      write(*,*) "!!!!!ERROR: PROGRAM DOESN'T KNOW WHICH ATOM IS NEEDED, EDIT INPUT.INP!!!!!"
      call exit(1)
    end if
  end if
    
  if(is_H_t) then
    if(.not. time_exist) call laser_me_setup
  end if
 
end subroutine need_files
      
  
  


!####################################################################!
!			       NOT USED				     !
!####################################################################!


subroutine calculate_KS_orbitals
USE GAUSS_DFT_ATOM
USE LINALG
USE MATLIB
USE GAUSS_BASIS
implicit none

  double precision, allocatable :: h(:,:),o(:,:),v(:,:),e(:)
  integer                       :: i,j,n,k
  double precision              :: a,b,ra(3),rb(3),rc(3),numr,ovmr,kimr
  integer                       :: la,ma,lb,mb



n=basis_dim
allocate(h(n,n),o(n,n),e(n),v(n,n))
  rc=0.d0
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
      call Coulomb_potential_m(a,ra,la,b,rb,lb,rc)
      call gauss_kinetic_m(a,ra,la,b,rb,lb)
      call okmat_real(la,ma,lb,mb,a,b,ra,rb,ovmr,kimr)
      call nucmat_real(la,ma,lb,mb,a,b,ra,rb,rc,numr)
      write(66,*)i,j
      write(66,*)ovmr,som(ma,mb)
      write(66,*)kimr,stm(ma,mb)
      write(66,*)numr,cpm(ma,mb)
    end do
  end do
  
deallocate(h,o,e,v)

end subroutine calculate_KS_orbitals


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
