module c_diff
  implicit none

contains

  function cent_diff(U, shape_U, delta, dim)

    integer, dimension(2) :: shape_U
    integer:: dim
    real, dimension(shape_U(1),shape_U(2)), intent(in) :: U
    real :: cent_diff(shape_U(1),shape_U(2))
    real :: delta
    
   if (dim == 1) then
  cent_diff = (CSHIFT(U, SHIFT = 1, DIM = 1)- CSHIFT(U, SHIFT = -1, DIM = 1) )/(2*delta)
   else if (dim ==2) then
  cent_diff = (CSHIFT(U, SHIFT = 1, DIM = 2)- CSHIFT(U, SHIFT = -1, DIM = 2) )/(2*delta)
   else 
   print*, "You have entered an invalid dimension"
   end if
    
  end function cent_diff
  
end module