module procedures
  implicit none

contains

  function m2e_x(U, shape_U)

    real, dimension(:,:), intent(in) :: U
    integer, dimension(:) :: shape_U
    real :: m2e_x(shape_U(1),shape_U(2))
    
   
    m2e_x=0.5*(U + CSHIFT(U, SHIFT = -1, DIM = 1))
    
  end function m2e_x
end module









program test
    use procedures

    implicit none
     real, dimension(4,3) :: a, edge_a
     !integer, dimension(2):: shape_a
     
    
    a = reshape((/ 1,5,4,7,2,8,1,4,3,9,1,11/),(/4,3/))
   ! shape_a = shape(a)
    edge_a =m2e_x(a,shape(a))
    
    print*,  edge_a !reshape(edge_a, (/4,3/))
   ! print*, CSHIFT(a, SHIFT = -1, DIM = 1)
    
end program