!=========================================================================!
! Subroutine that builds up the Hamiltonian depending on the initial      !
! conditions for the chain:                                               !
!  - CUSTOM                                                               !
!-------------------------------------------------------------------------!
! Parameters                                                              !
!   N            : number of sites                                        !
!   num_basis    : total number of basis vectors                          !
!   Js2D         : array with the coupling pattern                        !
!   basis_vecs   : the array of basis vectors                             !
!   H            : hamiltonian matrix to be returned                      !
!-------------------------------------------------------------------------!
! Written by Luke Mortimer, October 2019                                  !
!=========================================================================!

subroutine build_hamiltonian_custom(N, num_basis, Js2D, basis_vecs, H, eeScale)

    use constants

    ! Inputs
    integer, intent(in)                                  :: N
    integer, intent(in)                                  :: num_basis
    real(kind=dbl), dimension(N, N), intent(in)          :: Js2D
    integer, dimension(num_basis, num_basis), intent(in) :: basis_vecs
    real(kind=dbl), intent(in)                           :: eeScale

    ! Outputs
    real(kind=dbl), dimension(num_basis, num_basis), intent(inout) :: H

    ! Loop counters
    integer :: i, j, k, l

    ! For counting the number of 1s in either basis vector
    integer :: count_i, count_j

    ! For counting the number of differences between two basis vectors
    integer :: count_diff

    ! For storing the combination of two basis vectors
    integer, dimension(N) :: combined

    H = 0

    ! Loop through the basic vectors (squared)
    do i = 1, num_basis
        do j = 1, num_basis

            ! For the diagonal terms, set the relative energies
            if (i == j) then
                
                ! Add each energy for that site
                do k = 1, N
                    if (basis_vecs(i, k) == 1) then
                        H(j, i) = H(j,i) + Js2D(k, k)
                    end if
                end do

                ! Consider the interaction between adjacent excitations TODO 1
                do k = 1, N
                    do l = 1, N
                        if (basis_vecs(i, k) == 1 .and. basis_vecs(i, l) == 1) then
                            H(j, i) = H(j,i) + eeScale*Js2D(k, l)
                        end if
                    end do
                end do

                cycle

            end if

            ! Count the number of 1s in each vector
            count_i = sum(basis_vecs(i, 1:N))
            count_j = sum(basis_vecs(j, 1:N))

            ! If we are not in the same hilbert subspace, cycle
            if (count_i /= count_j) cycle 

            ! Count how many differences from i to j
            count_diff = 0
            do k = 1, N
                if (basis_vecs(i, k) == 1 .and. basis_vecs(j, k) == 0) then
                        count_diff = count_diff + 1
                end if
            end do

            ! Only look at differences of 1
            if (count_diff == 1) then

                ! Combine the two vectors
                combined(1:N) = basis_vecs(i, 1:N) - basis_vecs(j, 1:N)

                ! Search for a -1 and 1 pair, then see if they are coupled
                do k = 1, N

                    if (combined(k) == -1) then

                        do l = 1, N

                            if (combined(l) == 1) then
                                H(j, i) = Js2D(k, l)
                            end if

                        end do

                    endif

                enddo

            endif

        enddo
    enddo

end subroutine







