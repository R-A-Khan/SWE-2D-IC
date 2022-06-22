
!----------------------------------------------------------

module procedures
  implicit none

contains

 !______________________________________________________ 

  function m2e_x(U, shape_U)

    double precision, dimension(:,:), intent(in) :: U
    integer, dimension(:) :: shape_U
    double precision :: m2e_x(shape_U(1),shape_U(2))
    
   
    m2e_x= 5.D-1*(U + CSHIFT(U, SHIFT = -1, DIM = 1))
    
  end function m2e_x
  
 !______________________________________________________ 
  
 function m2e_y(U, shape_U)

    integer, dimension(2) :: shape_U
    double precision, dimension(shape_U(1),shape_U(2)), intent(in) :: U
    double precision :: m2e_y(shape_U(1),shape_U(2))
    
   
    m2e_y=5.D-1*(U + CSHIFT(U, SHIFT = -1, DIM = 2))
    
  end function m2e_y
  
 !______________________________________________________ 


function e2m_x(U, shape_U)

    integer, dimension(2) :: shape_U
    double precision, dimension(shape_U(1),shape_U(2)), intent(in) :: U
    double precision :: e2m_x(shape_U(1),shape_U(2))
    
   
   e2m_x= 5.D-1*(U + CSHIFT(U, SHIFT = 1, DIM = 1))
    
  end function e2m_x
  
 !______________________________________________________ 
 
   function e2m_y(U, shape_U)

    integer, dimension(2) :: shape_U
    double precision, dimension(shape_U(1),shape_U(2)), intent(in) :: U
    double precision :: e2m_y(shape_U(1),shape_U(2))
    
   
   e2m_y= 5.D-1*(U + CSHIFT(U, SHIFT = 1, DIM = 2))
    
  end function e2m_y
  
 !______________________________________________________  

  function cent_diff(U, shape_U, delta, dim)

    integer, dimension(2) :: shape_U
    integer:: dim
    double precision, dimension(shape_U(1),shape_U(2)), intent(in) :: U
    double precision :: cent_diff(shape_U(1),shape_U(2))
    double precision :: delta
    
   if (dim == 1) then
  cent_diff = (CSHIFT(U, SHIFT = 1, DIM = 1)- CSHIFT(U, SHIFT = -1, DIM = 1) )/(2.D0*delta)
   else if (dim ==2) then
  cent_diff = (CSHIFT(U, SHIFT = 1, DIM = 2)- CSHIFT(U, SHIFT = -1, DIM = 2) )/(2.D0*delta)
   else 
   print*, "You have entered an invalid dimension"
   end if
    
  end function cent_diff
  
 !______________________________________________________  
  
  function forw_diff(U, shape_U, delta, dim)

    integer, dimension(2) :: shape_U
    integer:: dim
    double precision, dimension(shape_U(1),shape_U(2)), intent(in) :: U
    double precision :: forw_diff(shape_U(1),shape_U(2))
    double precision :: delta
    
   if (dim == 1) then
   forw_diff = (CSHIFT(U, SHIFT = 1, DIM = 1) - U)/delta
   else if (dim ==2) then
   forw_diff = (CSHIFT(U, SHIFT = 1, DIM = 2) - U)/delta
   else 
   print*, "You have entered an invalid dimension"
   end if
    
  end function forw_diff
  
 !______________________________________________________  
 
   function back_diff(U, shape_U, delta, dim)

    integer, dimension(2) :: shape_U
    integer:: dim
    double precision, dimension(shape_U(1),shape_U(2)), intent(in) :: U
    double precision :: back_diff(shape_U(1),shape_U(2))
    double precision :: delta
    
   if (dim == 1) then
   back_diff = (U - CSHIFT(U, SHIFT = -1, DIM = 1) )/delta
   else if (dim ==2) then
   back_diff = (U - CSHIFT(U, SHIFT = -1, DIM = 2) )/delta
   else 
   print*, "You have entered an invalid dimension"
   end if
    
  end function back_diff
  
 !______________________________________________________    

 function linspace_rk(a,b,n)
  ! Usage: creates an equally spaced array of size n between a and b

  double precision, intent(in) :: a,b
  integer, intent(in) :: n
  double precision :: linspace_rk(n)
  double precision :: h
  integer :: i
  h=(b-a)/(n-1)
  linspace_rk=a+h*(/(i,i=0,n-1)/)
end function linspace_rk

 !_______________________________________________________
 

 
end module