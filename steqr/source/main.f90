program main
    use utility
    implicit none

    integer::N
    real*8, allocatable, dimension(:)::alpha, beta

    real*8, allocatable, dimension(:, :)::eigvec

    integer::i

    N = count_lines("alpha.txt")
    open(unit=99, file="alpha.txt", status="old")
        allocate(alpha(N))
        read(99, *)alpha
    close(99)
    open(unit=99, file= "beta.txt", status="old")
        allocate(beta(N - 1))
        read(99, *)beta
    close(99)

    allocate(eigvec(N, N))
    i = My_dsteqr(alpha, beta, eigvec, N)
    alpha = alpha / 4.556335830019422d-6

    write(*,*)"Energy / cm^-1"//char(9)//"amplitude"//char(9)//"strength"
    do i = 1, N
        write(*,*)alpha(i), char(9), eigvec(1, i), char(9), eigvec(1, i) * eigvec(1, i)
    end do
    
end program main