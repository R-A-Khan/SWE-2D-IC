module bw_diff
  implicit none

contains

  function back_diff(U, shape_U, delta, dim)

    integer, dimension(2) :: shape_U
    integer:: dim
    real, dimension(shape_U(1),shape_U(2)), intent(in) :: U
    real :: back_diff(shape_U(1),shape_U(2))
    real :: delta
    
   if (dim == 1) then
   back_diff = (U - CSHIFT(U, SHIFT = -1, DIM = 1) )/delta
   else if (dim ==2) then
   back_diff = (U - CSHIFT(U, SHIFT = -1, DIM = 2) )/delta
   else 
   print*, "You have entered an invalid dimension"
   end if
    
  end function back_diff
end module