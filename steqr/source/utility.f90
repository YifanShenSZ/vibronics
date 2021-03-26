module utility
    implicit none

contains
integer function count_lines(file)
    character(*), intent(in)::file
    integer::i
    character*128::chartemp
    count_lines = 0
    open(unit=99, file=file, status="old", iostat=i)
        do
            read(99, *, iostat=i)chartemp
            if (i /= 0) exit
            count_lines = count_lines + 1
        end do
    close(99)
end function count_lines

! diag harvests the eigenvalues
integer function My_dsteqr(diag, subdiag, eigvec, N) result(info)
    integer, intent(in)::N
    real*8, dimension(N), intent(inout)::diag
    real*8, dimension(N - 1), intent(inout)::subdiag
    real*8, dimension(N, N), intent(out)::eigvec
    real*8, dimension(max(1, 2 * N - 2))::work
    call dsteqr('I', N, diag, subdiag, eigvec, N, work, info)
end function My_dsteqr

end module utility