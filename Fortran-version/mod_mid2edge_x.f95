module mid2edge_x
  implicit none

contains

  function m2e_x(U, shape_U)

    real, dimension(:,:), intent(in) :: U
    integer, dimension(:) :: shape_U
    real :: m2e_x(shape_U(1),shape_U(2))
    
   
    m2e_x=0.5*(U + CSHIFT(U, SHIFT = -1, DIM = 1))
    
  end function m2e_x
end module




