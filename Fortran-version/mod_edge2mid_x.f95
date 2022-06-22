module edge2mid_x
  implicit none

contains

  function e2m_x(U, shape_U)

    integer, dimension(2) :: shape_U
    real, dimension(shape_U(1),shape_U(2)), intent(in) :: U
    real :: e2m_x(shape_U(1),shape_U(2))
    
   
   e2m_x=0.5*(U + CSHIFT(U, SHIFT = 1, DIM = 1))
    
  end function e2m_x
end module