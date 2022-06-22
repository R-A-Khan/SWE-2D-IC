program test
    use c_diff

    implicit none
     real, dimension(4,3) :: a, edge_a
     real:: stepsize
     integer:: dimension
     
    
    a = reshape((/ 1,5,4,7,2,8,1,4,3,9,1,11/),(/4,3/))
    stepsize = 3
    dimension = 2
   ! shape_a = shape(a)
    edge_a =cent_diff(a,shape(a), stepsize, dimension)
    
    print*,  edge_a !reshape(edge_a, (/4,3/))
 
    
end program