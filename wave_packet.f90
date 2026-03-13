program wavepacket
    implicit none

    ! Parameters

    integer, parameter :: Nx = 400
    integer, parameter :: Nk = 600
    integer, parameter :: Nt = 50

    real(8), parameter :: a = 1.0
    real(8), parameter :: xmin = -10.0
    real(8), parameter :: xmax = 10.0
    real(8), parameter :: kmin = -20.0
    real(8), parameter :: kmax = 20.0

    real(8) :: dx,dk,pi
    real(8) :: x(Nx), k(Nk)

    complex(16) :: psi(Nx)
    real(8) :: phi(Nk)
    integer :: i, j
    real(8) :: prob

     real(8) :: tmax, t, dt
    real(8) :: omega
    complex(16) :: sum,term
    integer :: n

    ! Initialize the spatial and momentum grids
    dx = (xmax - xmin) / (Nx - 1)
    dk = (kmax - kmin) / (Nk - 1)

    ! Initialize the spatial grid
    do i = 1, Nx
        x(i) = xmin + (i - 1) * dx
    end do

    ! Initialize the momentum grid
    do j = 1, Nk
        k(j) = kmin + (j - 1) * dk
    end do

   
    pi = acos(-1.0)

    ! compute phi(k)
    do j = 1,Nk
        if (abs(k(j)) .lt. 0.00001) then
            phi(j) = sqrt(2.0/(pi*a)) * (a/2)
        else
            phi(j) = sqrt(2.0/(pi*a)) * sin(k(j)*a/2)/k(j)
        endif
    end do

    ! time parameters

   

    tmax = 5.0
    dt = tmax / Nt

    open(10, file='wavepacket.dat', status='unknown')

    ! Time evolution loops

    do n= 1,Nt
        t = (n-1) * dt
        do i = 1, Nx
            sum = (0.0, 0.0)
            do j = 1,Nk
                omega = k(j)**2 / 2.0
                term = phi(j) * exp((0.0, 1.0) * (k(j)*x(i) - omega*t))
                if ( j ==1 .or. j == Nk) then
                    sum = sum + term
                else if (mod(j,2).eq.0) then
                    sum = sum + 4.0*term
                else
                    sum = sum + 2.0*term
                end if
                
            end do
            psi(i) = (dk/3.0)*sum/sqrt(2.0*pi)

                ! compute probability density
            prob = real(psi(i))**2 + aimag(psi(i))**2
            write(10,*) x(i), prob
        end do
        write(10,*)
        write(10,*)
    end do
    close(10)
end program wavepacket

    