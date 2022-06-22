module mid2edge_y
  implicit none

contains

  function m2e_y(U, shape_U)

    integer, dimension(2) :: shape_U
    real, dimension(shape_U(1),shape_U(2)), intent(in) :: U
    real :: m2e_y(shape_U(1),shape_U(2))
    
   
    m2e_y=0.5*(U + CSHIFT(U, SHIFT = -1, DIM = 2))
    
  end function m2e_y
end module

