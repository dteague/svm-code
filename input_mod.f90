module input_mod

! input file must be called input.inp and there must be a ':','=', or '-'
! between the name of the value and the value itself.  There can be tabs in 
! the input file.
  
  implicit none
  integer		:: ntsteps_inp, T_of_output, half_lat_num, log_num, time_step=0, ll_max
  double precision 	:: dt_inp, k_vec_inp(3), space_res(3), space, nu_val, dummy
  character (len=255)	:: file_line
  logical 		:: taylor, is_helium, is_hydrogen, basisf_exist
  double precision      :: Z_charge_inp, width, lambda, I_0
  LOGICAL 		:: file_exists
  
  integer, parameter	:: logfile = 402
  
  
  CONTAINS

! var_input
! Description: reads in the file "input.inp' and gets the data from the file
! post: Since this is a module, this file stores all of needed variable data
! 	in it to be accessed by other programs that use this module
subroutine var_input
  integer,parameter	:: input = 401, basisf = 403
  
  INQUIRE(FILE="input.inp", EXIST=file_exists)   ! file_exists will be TRUE if the file
						  ! exists and FALSE otherwise
  if(file_exists) then
    continue
  else
    write(*,"(/,5x,A,/)") "!!!!!ERROR: input.inp not found!!!!!"
    call exit(1)
  end if

  open(input, file='input.inp')
 
 
  INQUIRE(FILE="basis.inp", EXIST=file_exists)
  
  
  !used for giving overlap and hamiltonian matrix files unique names
  if(file_exists) then
    basisf_exist = .true.
    
    open(basisf, file="basis.inp")
    read(basisf,*)
    read(basisf,*)ll_max
    read(basisf,*)
    read(basisf,*)nu_val
    read(basisf,*)space
  else
    basisf_exist = .false.
  end if
  
  
  
  do
    read(input,'(a)',end=110) file_line
    call parse(file_line)
  end do
  
  110 close(input)
  !dis_psi and dis_density are useless right now since the program outputs
  !both the wavefunction and the density to files.

   
  !only for now, space_res will be the same in all directions
  space_res(2)=space_res(1)
  space_res(3)=space_res(1)
  
end subroutine

! parse
! Description: parses a line for a tag that matches the data needed and gets
! 	that value
! Parameters:	the_key - tag word to be tested
! 		value_string - string of the value in the line for the tag
! 		file_line - line of text sent to the subroutine to be parsed
! post: Either assigns data to the appropiate tag in the module if the key matches
! 	a tag or it outputs an error line that the key was not recognized
SUBROUTINE parse(file_line)
  implicit none
  integer		:: i
  character*255    	:: the_key,value_string
  character(len=*) 	:: file_line

   open(logfile, file="log.dat")
  
  ! Ignore anything on the line that appears after a !
  i=index(file_line,'!')
  if(i>0) then
    file_line=file_line(1:i-1)
  endif

!	Checks for :,=,- which seperates the name and value
  i=index(file_line,':')
  if(i==0) then
    i=index(file_line,'=')
    if(i==0) then
      i=index(file_line,'-')
      if(i==0) then
	write(logfile,"(a)") "ERROR: Did not read the line"
	write(logfile,*) trim(adjustl(file_line))
	return
      end if
    end if
  end if
  
  read(file_line(1:i-1),'(a)')the_key
  read(file_line(i+1:),'(a)')value_string
    
  call convert_to_lowercase(the_key)

!	select case structure used to test if key has a corresponding tag
!	Allows for multiple names for a key
  select case(trim(adjustl(remove_tabs(the_key))))
    case("nt", "#_timesteps")
      read(value_string,*) ntsteps_inp
      if(ntsteps_inp<0) then
	write(logfile,"(5x,a,5x)") "!!!!!ERROR: number of timesteps is a negative integer!!!!!"
	write(logfile,*) "nt = ", ntsteps_inp
	write(logfile,*) "Defaulting to nt = 100"  !!!!!!must come up with default or just end program!!!!!
	ntsteps_inp = 100
      else
	write(logfile,*) "nt = ", ntsteps_inp
      end if
    
    case("is_helium?", "is_helium", "helium")
      call isBoolean(value_string,is_helium)
      if(is_helium) then
	write(logfile,*) "Used Helium"     
	
      end if

    case("is_hydrogen?", "is_hydrogen", "hydrogen")
      call isBoolean(value_string,is_hydrogen)
      if(is_hydrogen) then
	write(logfile,*) "Used Hydrogen"     
      end if
      
      
    case("taylor", "use_taylor")
      call isBoolean(value_string,taylor)
      write(logfile,*) "Use_taylor? ", taylor
      
    case("dt")
      read(value_string,*) dt_inp
      write(logfile,*) "dt = ", dt_inp
      
    case("z_charge")
      read(value_string,*) Z_charge_inp
      write(logfile,*) "Z_Charge = ", Z_charge_inp
    
    case("t_of_output", "period_output")
      read(value_string,*) T_of_output
      if(T_of_output<1) then
	write(logfile,"(5x,a,5x)") "!!!!!ERROR: Period of output is not a positive integer!!!!!"
	write(logfile,*) "Period_output = ", T_of_output
	write(logfile,*) "Defaulting to T_of_output = 10"  !!!!!!must come up with default or just end program!!!!!
	T_of_output = 10
      else
	write(logfile,*) "Period_output = ", T_of_output
      end if
    
    case("lat_width", "space_res", "space_res(1)")
      read(value_string,*) space_res(1)
      write(logfile,*) "Space_res(1) = ", space_res(1)
    
    case("half_lat_num")	!!!!!!!!!work on variable
      read(value_string,*) half_lat_num
      write(logfile,*) "Half_lattice_num = ", half_lat_num
    
    case("k_x", "k(1)", "k_vec(1)")
      read(value_string,*) k_vec_inp(1)
      write(logfile,*) "k_x = ", k_vec_inp(1)
      
    case("k_y", "k(2)", "k_vec(2)")
      read(value_string,*) k_vec_inp(2)
      write(logfile,*) "k_y = ", k_vec_inp(2)
      
    case("k_z", "k(3)", "k_vec(3)")
      read(value_string,*) k_vec_inp(3)
      write(logfile,*) "k_z = ", k_vec_inp(3)
      
    case("width")
      read(value_string,*) width
      write(logfile,*) "width = ", width
      
    case("lambda")
      read(value_string,*) lambda
      write(logfile,*) "lambda = ", lambda
      
    case("i_0", "i")
      read(value_string,*) I_0
      write(logfile,*) "I_0 = ", I_0
 
    case("disp_wavefunc", "disp_density")	
      continue
      
    case default
       write(logfile,'(2a)')"ERROR: Invalid variable name: ",trim(adjustl(the_key))
       write(*,'(2a)')"!!!!!ERROR: Invalid variable name: ",trim(adjustl(the_key))
       stop
       
  end select
 
END SUBROUTINE parse


! isBoolean
! Desciption: converts a string to a boolean value. 
! post: Tests if the line begins with .t or t.  If so, the string is returned true
! 	In all other cases, the string returns false
subroutine isBoolean(string, bool_val)
  character(len=*)	:: string
  logical		:: bool_val
  
  call convert_to_lowercase(string)
    if((trim(adjustl(string))=="t").OR.(trim(adjustl(string))=="true") &
          .OR.(trim(adjustl(string))==".true.")) then
      bool_val=.TRUE.
    else
      bool_val=.FALSE.
    endif

end subroutine isBoolean

! convert_to_lowercase
! Description: Converts a string to lowercase.  Does not affect non-letter
! 	characters
! post: String of same length that has only lowercase letters
SUBROUTINE convert_to_lowercase(string)
  integer          :: i
  character(len=*) :: string
  character        :: the_char

  do i=1,len(string)
    the_char=string(i:i)
    if((iachar(the_char)>=iachar("A")).AND.(iachar(the_char)<=iachar("Z"))) then
      the_char=achar(iachar(the_char)+(iachar('a')-iachar('A')))
    endif
    string(i:i)=the_char
  end do
END SUBROUTINE convert_to_lowercase


! remove_tabs
! Description: Removes the text after a tab in a string.  Allows for the data
! 	to contain tabs in it to aid in readability
! Post: returns string of length 255 without the text after the first tab.  Just
! 	returns the original string if there are no tabs 
function remove_tabs(string)
  character (len=*)	:: string
  character (len=255) 	:: remove_tabs
  integer		:: i
  
  remove_tabs=string
  i=index(string,'	')
  if(i>0) then
    remove_tabs=string(1:i-1)
  end if
end function remove_tabs

end module




