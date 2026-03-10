program jacobi_eigen
    implicit none

    integer, parameter :: n=5
    integer :: i,j

    real(8) :: A(n,n)   ! The input matrix
    real(8) :: V(n,n)   ! The matrix of eigenvectors
    real(8) :: d(n)     ! The eigenvalues(diagonal elements)

    integer :: iter, p, q
    real(8) :: S
    integer, parameter :: max_iter = 50


    real(8) :: t, c, s, tau
    real(8) :: theta, h

    ! Reading the matrix from the file 'matrix.dat'
    open(unit=10, file='matrix.dat', status='old')
    do i=1,n
        read(10,*) (A(i,j), j=1,n)
    enddo
    close(10)
    
    ! print the input matrix to verify
    !print *, "Input matrix A:"
    !do i=1,n
    !    print *, (A(i,j), j=1,n)
    !enddo


    ! Intialize V to the identity matrix
    do i = 1,n
        do j=1,n
            if (i==j) then
                V(i,j) = 1.0
            else
                V(i,j) = 0.0
            endif
        enddo
    enddo
    
    ! Initialize d to the diagonal elements of A
    do i = 1,n
        d(i) = A(i,i)  ! These values will change during rotation and eventually becone the eigenvalues
    enddo

    ! checking if the matrix is already diagonal
    ! Start the Jacobi iteration loop

    do iter = 1, max_iter
        S = 0.0
        do i = 1, n-1 ! off-diagonal summation
            do j = i+1,n
                S = S+abs(A(i,j))
            enddo
        enddo
        
        ! Check for convergence
        if (S==0.0) then
            exit
        endif
    enddo

    ! Jacobi rotation

    do p = 1, n-1
        do q = p+1, n
            if (abs(A(p,q)) > 0.0) then
                h = d(q) - d(p)

                if (abs(h) + abs(A(p,q)) == abs(h)) then
                    t = A(p,q) / h
                else
                    theta = 0.5 * h / A(p,q)
                    t = 1.0 / (abs(theta) + sqrt(1.0 + theta*theta))
                    if (theta < 0.0) t = -t
                    end if
                


                c = 1.0 / sqrt(1.0 + t*t)
                s = t * c
                tau = s / (1.0 + c)
            endif
        
            if (abs(A(p,q))>0.0) then
                h = t * A(p,q)
                z(p) = z(p) - h
                z(q) = z(q) + h
                d(p) = d(p) - h
                d(q) = d(q) + h
                A(p,q) = 0.0
                

               



            

    

    


        


