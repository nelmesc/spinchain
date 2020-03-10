!=========================================================================!
! Subroutine that computes the dynamics of the system given a initial     !
! vector.                                                                 !
!                                                                         !
! REMEMBER: each column of the returned matrix corresponds to the         !
! fidelity of the initial state against the vector with the same number   !
! than that column. First column is time                                  !
!-------------------------------------------------------------------------!
! Return value:                                                           !
!   Fidelity => matrix with the dynamics of the system                    !
!-------------------------------------------------------------------------!
! Written by Marta Estarellas, v0.1, 31/03/2017                           !
!=========================================================================!

subroutine injection_dynamics(HT,hami,eigvals,vectorstotal,c_i,initialtime,modeNum,c_vs_t)

!modules
use constants
use parameters

!inputs
integer, intent(in) :: vectorstotal
integer, dimension(vectorstotal,vectorstotal), intent(in) :: HT
integer, intent(in) :: modeNum

! Optional array which if present returns the coefficients as a function of time
complex(kind=dbl), dimension(:, :), intent(inout), allocatable, optional :: c_vs_t
integer :: step_counter

real(kind=dbl), dimension(vectorstotal), intent(in) :: eigvals

complex(kind=dbl), dimension(vectorstotal,vectorstotal), intent(in) :: hami

complex(kind=dbl), dimension(vectorstotal), intent(inout) :: c_i
real(kind=dbl), intent(out) :: initialtime

!local variables
integer :: i,j,k

real(kind=dbl) :: norm
real(kind=dbl) :: step_size
real(kind=dbl) :: time
real(kind=dbl) :: eof_rho

complex(kind=dbl) :: sum_vec

!arrays
real(kind=dbl), allocatable, dimension(:) :: fidelity
real(kind=dbl), dimension(N) :: siteProb

complex(kind=dbl), allocatable, dimension(:) :: a_m
complex(kind=dbl), allocatable, dimension(:,:) :: Y_o, Y_t
complex(kind=dbl), allocatable, dimension(:,:) :: red_rho
complex(kind=dbl), allocatable, dimension(:) :: initial_state

!allocate arrays
allocate(red_rho(4,4))
allocate(a_m(vectorstotal))
allocate(fidelity(vectorstotal))
allocate(Y_o(vectorstotal,vectorstotal))
allocate(Y_t(vectorstotal,vectorstotal))
allocate(initial_state(vectorstotal))

! If told to return all coefficients, allocate the array and init the counter
if (present(c_vs_t)) then
    if (single) then
        if (.not. allocated(c_vs_t)) allocate(c_vs_t(1, vectorstotal))
    else
        if (.not. allocated(c_vs_t)) allocate(c_vs_t(steps, vectorstotal))
    end if
    c_vs_t = cmplx(0.0_dbl, 0.0_dbl, dbl)
    step_counter = 1
end if

!initialize fidelity
fidelity=0._dbl

!open files
if (.not. use_genetic) then

    open (unit=44,file='dynamics.data',status='unknown')
    open (unit=45,file='exmap.data',status='unknown')
    open (unit=46,file='eof.data',status='unknown')

    !writte comments on the oytputs
    write(44,*) '#DYNAMICS. TIME (1st COL) AND EVOLVED COEFFICIENTS ALL THE BASIS VECTORS ORDERED BY INDEX NUMBER'
    write(45,*) '#DYNAMICS. TIME (1st COL) AND SITE OCCUPATION PROBABILITES ORDERED BY SITE'
    write(46,*) '#EOF. TIME (1st COL) AND EOF BETWEEN Q1 AND Q2'

end if

!normalisation factor
norm=(1._dbl/sqrt(float(numI)))

if (.not. custom) then
    if (.not.read_state) then
        !get initial vector from PARAMETERS file inputs
        call initialState(initialVec)
    else
        call readFromFile(vectorstotal,initial_state,initialtime)
    endif
end if

!Define |¥(0)> = \sum{a_m|m>} and a_m = <m|Y(0)>
!being <¥(0)|=norm*(<initialVec1|+<initialVec2|+...)
if (.not.read_state) then
do i=1,vectorstotal
    a_m(i)=0._dbl
    do j=1,numI
        if (initialVec(modeNum, j) == 0) cycle
        a_m(i) = a_m(i)+initialCoeff(modeNum,j)*(dconjg(hami(initialVec(modeNum,j),i)))
    enddo
enddo
else
!Initial state from the file
do i=1,vectorstotal
    a_m(i)=0._dbl
    do j=1,vectorstotal
        a_m(i) = a_m(i)+initial_state(j)*dconjg(hami(j,i))
    enddo
enddo
endif

step_size = totalTime/steps

if (.not.read_state) then
    initialtime = 0._dbl
endif

if (.not. use_genetic) then
    open (unit=48,file='initial_state.data',status='unknown')
    if (read_state) then
    write(48,*) '#STATE AT TIME T. FIRST ROW IS TIME, FOLLOWING ROWS ARE THE STATE ORDERED IN THE BASIS VECTORS.'
    write(48,*) time
    do i=1,vectorstotal
        write(48,*) initial_state(i)
    enddo
    endif
end if

time = 0._dbl
!main loop for the dynamics
if (.not.single) then

!****************************************************************
!!CALCULATE AT EVERY TIME FIDELITY, EXMAP, RHO, EOF AND ENTROPY *
!****************************************************************
do while (time<=totalTime)

    do i=1,vectorstotal
        do j=1,vectorstotal
            Y_o(i,j) = hami(i,j)*a_m(j)
            Y_t(i,j) = Y_o(i,j)*exp(-1._dbl*im*time*eigvals(j))
        enddo
    enddo

    !|¥(t)> = \sum{c_i|i>}
    c_i=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)
    do i=1,vectorstotal
        sum_vec = cmplx(0._dbl,0._dbl,kind=dbl)
        do j=1,vectorstotal
            sum_vec = sum_vec + Y_t(i,j)
        enddo
        c_i(i) = sum_vec
    enddo

    !Fidelities
    do i=1,vectorstotal
        fidelity(i) = real(dconjg(c_i(i))*c_i(i), dbl)
    enddo

    siteProb=0
    do i=1,vectorstotal
        do k=1,N
            if (HT(i,k)==1) then
                siteProb(k) = siteProb(k) + fidelity(i)
            endif
        enddo
    enddo

    ! Output to file if not optmising
    if (.not. use_genetic) then
        if (pst) then
            write(44,*) (initialtime+time)*J_max, c_i
            write(45,*) (initialtime+time)*J_max, siteProb
        else
            write(44,*) (initialtime+time)*J_max, c_i
            write(45,*) (initialtime+time)*J_max, siteProb
        endif
    end if

    ! If the c_vs_t array is there, copy the info over for this step
    if (present(c_vs_t) .and. step_counter <= steps) then
        c_vs_t(step_counter, 1:vectorstotal) = c_i(1:vectorstotal)
        step_counter = step_counter + 1
    end if

    !calculate eof for all times
    if (eof) then
        red_rho=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)
        call reduced_density_matrix(HT,vectorstotal,red_rho,c_i)
        call entanglement_of_formation(red_rho,eof_rho)
        if (.not. use_genetic) write(46,*) (initialtime+time)*J_max, eof_rho
    endif

    !increase value of time
    time = time + step_size

enddo
else if (single) then

!****************************************************************
!!CALCULATE AT JUST T_A FIDELITY, EXMAP, RHO, EOF AND ENTROPY ***
!****************************************************************

    time = t_A
    do i=1,vectorstotal
        do j=1,vectorstotal
            Y_o(i,j) = hami(i,j)*a_m(j)
            Y_t(i,j) = Y_o(i,j)*exp(-1._dbl*im*time*eigvals(j))
        enddo
    enddo

    !|¥(t)> = \sum{c_i|i>}
    c_i=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)
    do i=1,vectorstotal
        sum_vec = cmplx(0._dbl,0._dbl,kind=dbl)
        do j=1,vectorstotal
            sum_vec = sum_vec + Y_t(i,j)
        enddo
        c_i(i) = sum_vec
    enddo

    siteProb=0
    do i=1,vectorstotal
        do k=1,N
            if (HT(i,k)==1) then
                siteProb(k) = siteProb(k) + fidelity(i)
            endif
        enddo
    enddo

    !Fidelities
    do i=1,vectorstotal
        fidelity(i) = (abs(c_i(i)))**2
    enddo

    if (.not. use_genetic) then
        if (pst) then
            write(44,*) time*J_max, c_i
            write(45,*) time*J_max, siteProb
        else
            write(44,*) time*J_max, c_i
            write(45,*) time*J_max, siteProb
        endif
    end if

    ! If the c_vs_t array is there, copy the info over for the only step
    if (present(c_vs_t)) then
        c_vs_t(1, 1:vectorstotal) = c_i(1:vectorstotal)
    end if

    !calculate eof for ONE time
    if (eof) then
        red_rho=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)
        call reduced_density_matrix(HT,vectorstotal,red_rho,c_i)
        call entanglement_of_formation(red_rho,eof_rho)
        if (.not. use_genetic) write(46,*) time*J_max, eof_rho
    endif
endif

!******************************************
!!OUTPUT STATE AT THE END OF SIMULATION ***
!******************************************
if (.not. use_genetic) then
    open (unit=47,file='final_state.data',status='unknown')
    write(47,*) '#STATE AT TIME T. FIRST ROW IS TIME, FOLLOWING ROWS ARE THE STATE ORDERED IN THE BASIS VECTORS.'
    write(47,*) (time - step_size)
    do i=1,vectorstotal
        write(47,*) c_i(i)
    enddo

    !close files
    close(44)
    close(45)
    close(46)
    close(47)
    close(48)
end if

!deallocate
deallocate(a_m)
deallocate(Y_o)
deallocate(Y_t)
deallocate(red_rho)
deallocate(fidelity)
deallocate(initial_state)

end subroutine
