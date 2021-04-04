! diag harvests the eigenvalues, eigvec harvest the eigenvectors
integer*4 function My_dsteqr(diag, subdiag, eigvec, N) result(info)
    integer*4, intent(in)::N
    real*8, dimension(N), intent(inout)::diag
    real*8, dimension(N - 1), intent(inout)::subdiag
    real*8, dimension(N, N), intent(out)::eigvec
    real*8, dimension(max(1, 2 * N - 2))::work
    call dsteqr('I', N, diag, subdiag, eigvec, N, work, info)
end function My_dsteqr