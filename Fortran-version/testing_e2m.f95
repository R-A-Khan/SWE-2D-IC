program test
    use edge2mid_y

    implicit none
     real, dimension(4,3) :: a, edge_a
     !integer, dimension(2):: shape_a
     
    
    a = reshape((/ 1,5,4,7,2,8,1,4,3,9,1,11/),(/4,3/))
   ! shape_a = shape(a)
    edge_a =e2m_y(a,shape(a))
    
    print*,  edge_a !reshape(edge_a, (/4,3/))
   ! print*, CSHIFT(a, SHIFT = -1, DIM = 1)
    
end program