module noisemod

integer, parameter :: dubp=kind(1.0d0)

! For random numbers
integer, parameter :: bit_32=kind(4)
integer(bit_32) :: ix,iy

real(dubp) :: Noise_SDev=0.5_dubp ! For Normal distribution

end module


module main

contains 

subroutine solve(modeNum, c_vs_t)

!!load subroutines
use dependencies
!!load constants
use constants
!!load initial parameters
use parameters

use noisemod

implicit none

!!This program solves an XY spin chain problem. It diagonalizes the hamiltonian obtaining
!!for chains or networks of different geometries and distributions, and obtains its
!!eigenvectors and eigenvalues as well as evolves it dynamically with a defined set of
!!initial conditions. Properties such Fidelity, Entropy and Entanglement of Formation
!!can be computed by using the relevant subroutines
!a change

! Optional array which if present returns the coefficients as a function of time
complex(kind=dbl), dimension(:, :), allocatable, intent(inout), optional :: c_vs_t

! Optional variable which if present tells which mode to use (defaults to 1)
integer, intent(in) :: modeNum

!***************************************************!
!******************  VARIABLES *********************!
!***************************************************!

!integers
integer :: i,j,k,v,w   !loop dummies
integer :: nit,Ninit,ex      !subroutine Permutations variables
integer :: seed

integer :: vectors1ex       !Initially set to N, reallocate later if needed (E.I.)
integer :: vectors2ex       !Initially set to N, reallocate later if needed (E.I.)
integer :: vectors3ex       !Initially set to N, reallocate later if needed (E.I.)
integer, dimension(:), allocatable :: vectors_ex
integer :: vectorstotal         !Sum of all the vectors
integer :: hub
integer :: len_branch

integer :: info, liwork  !Info in lapack subroutines

integer, dimension(8) :: values !array with date
integer, dimension(:), allocatable :: vec
integer, dimension(branches-1) :: limits
integer, allocatable, dimension (:) :: iwork

integer, allocatable, dimension(:,:) :: HS ! Hilbert subspace matrices
integer, allocatable, dimension(:,:) :: HT ! Total subspace matrix (containing basis vecs)
integer, allocatable, dimension(:,:) :: H1,H2,H3
integer :: vectors_so_far


!floats
real(kind=dbl) :: normal,orto !normaliztion constant
real(kind=dbl) :: r !random number
real(kind=dbl) :: initialtime

real(kind=dbl), dimension(:),   allocatable :: Js
real(kind=dbl), dimension(:,:), allocatable :: Js2D

real(kind=dbl), allocatable, dimension(:) :: eigvals
real(kind=dbl), allocatable, dimension(:) :: rwork

real(kind=dbl), allocatable, dimension(:) :: Noise

!complex
complex(kind=dbl), allocatable, dimension(:) :: work
complex(kind=dbl), allocatable, dimension(:) :: c_i

complex(kind=dbl), allocatable, dimension(:,:) :: hamiD !diagonalized Hamiltonian

!matrices
real(kind=dbl), allocatable, dimension(:,:) :: hami !Hamiltonian

!strings
character :: a
character(len=32) :: tmp, tmp2
character(len=500) :: fmt1,fmt2 !format descriptors

!random number generator
real(dubp), external :: algor_uniform_random

!logicals
logical, parameter, dimension(4) :: topology = (/linear,star,lattice,squared/)
logical, parameter, dimension(7) :: coupling = (/uniform, pst, ssh_a, ssh_b, abc, kitaev, ballistic/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! START PROGRAM AND WRITE OUTPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (output .and. .not. use_genetic) then
!retrieve date
call date_and_time(VALUES=values)
101 FORMAT (1X,59("*"))
102 FORMAT (1X,16("*")," SPIN CHAIN PROGRAM OUTPUT ",16("*"))
103 FORMAT (20X,I2,"/",I2.1,"/",I4,2X,I2,":",I2)
104 FORMAT (1X,59("-"))
open (unit=40,file='spinchain.out',status='replace')
write(40,101)
write(40,102)
write(40,101)
write(40,*) '           © Marta P. Estarellas, 27/07/2016              '
write(40,*) '                   University of York                     '
write(40,103) values(3),values(2),values(1),values(5),values(6)
write(40,104)
write(40,*)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DEFINING THE DESIRED TYPE OF CHAIN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (allocated(vec)) deallocate(vec)
allocate(vec(N))

! Set up variables with N dependence
vectors1ex = N
vectors2ex = N
vectors3ex = N
if(allocated(Js)) deallocate(Js)
if(allocated(Js2D)) deallocate(Js)
allocate(Js(N-1))
allocate(Js2D(N,N))
Js = 0.0_dbl
Js2D = 0.0_dbl

!**********************************************
!this is done in the module called PARAMETERS
!you should ONLY modify that module to set the
!conditions and structure of the chain.
!**********************************************

!if (output .and. .not. use_genetic) write(*,*) '>> Defining System'

!Write initial conditions
if (output .and. .not. use_genetic) then

    201 FORMAT (/"NUMBER OF SITES = ",A)
    write(tmp,'(i4.1)') N
    write(40,201) adjustl(trim(tmp))

    202 FORMAT ("NUMBER OF EXCITATIONS = ",A)
    write(tmp,'(i4.1)') exno
    write(40,202) adjustl(trim(tmp))

    203 FORMAT ("TOPOLOGY = ",A)
    if (custom) then
        tmp = 'CUSTOM'
    else if (linear) then
        tmp = 'LINEAR'
    else if (star) then
        tmp = 'STAR'
    else if (lattice) then
        tmp = 'LATTICE'
    else if (squared) then
        tmp = 'SQUARED'
    end if
    write(40,203) adjustl(trim(tmp))

    204 FORMAT ("COUPLING CONFIGURATION = ",A)
    if (custom) then
        tmp = "CUSTOM"
    else if (uniform) then
        tmp = 'UNIFORM'
    else if (pst) then
        tmp = 'PST'
    else if (ssh_a) then
        tmp = 'SSH - A'
    else if (ssh_b) then
        tmp = 'SSG - B'
    else if (abc) then
        tmp = 'ABC'
    else if (ballistic) then
        tmp = 'BALLISTIC'
    end if
    write(40,204) adjustl(trim(tmp))

    if (custom) then
        write(40, "(A,A)") "CUSTOM NETWORK STRING = ", trim(custom_string)
    end if

    206 FORMAT ("DIAGONAL DISORDER = ",A)
    tmp = 'NO'
    if (random_D) then
        write(tmp,'(f8.2)') E_D
    endif
    write(40,206) adjustl(trim(tmp))

    207 FORMAT ("OFF-DIAGONAL DISORDER = ",A)
    tmp = 'NO'
    if (random_J) then
        write(tmp,'(f8.2)') E_J
    endif
    write(40,207) adjustl(trim(tmp))

    208 FORMAT ("TOTAL TIME = ",A)
    write(tmp,'(f8.2)') totaltime
    write(40,208) adjustl(trim(tmp))

    209 FORMAT ("TIME STEP = ",A)
    write(tmp,'(f8.4)') totaltime/real(steps)
    write(40,209) adjustl(trim(tmp))


    if (eof) then
        210 FORMAT ("EOF = Qubit ",A,'- Qubit ',A)
        write(tmp,'(i3.1)') Q1
        write(tmp2,'(i3.1)') Q2
        write(40,210) adjustl(trim(tmp)), adjustl(trim(tmp2))
    else
        211 FORMAT ("EOF = ",A)
        tmp = 'NO'
        write(40,211) adjustl(trim(tmp))
    endif

    212 FORMAT ("METHOD = ",A)
    if (integration) then
        tmp = 'INT'
    else if (diagonalisation) then
        tmp = 'DIAG'
    end if
    write(40,212) adjustl(trim(tmp))

endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! INITIAL CHECKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (count(topology).ne.1) then
STOP 'ERROR: only one topology can be chosen. Check your PARAMETERS file.'
endif

if (count(coupling).ne.1) then
STOP 'ERROR: only one coupling configurations can be chosen. Check your PARAMETERS file.'
endif

if (integration.and.diagonalisation) then
STOP 'ERROR: only one method can be chosen. Check your PARAMETERS file.'
endif

if (linear) then
    if (ssh_a.or.ssh_b) then
    if (MOD(N-1,4)/=0) then
    STOP 'ERROR: for type (a) ssh chain N needs to be odd and N-1 needs to be divisible by 4.'
    endif
    endif

    if (abc) then
    if (MOD(N-3,4)/=0) then
    STOP 'ERROR: for type ABC linear chain N needs to be odd and N-3 needs to be divisible by 4.'
    endif
    endif

    if (kitaev) then
    if (MOD(N,2)/=0) then
    STOP 'ERROR: for a kitaev chain N needs to be even.'
    endif
    endif
endif

if (star .and. .not. custom) then

    if (branches==0) then
    STOP 'ERROR: For star topology you need to specify how many branches you want.'
    endif

    if (MOD((N-1),branches)/=0) then
    write(tmp,'(i2.1)') branches
    write(*,*) 'ERROR: ',adjustl(trim(tmp)),'-branched networks need (N-1) sites divisible by ',adjustl(trim(tmp)),'.'
    STOP
    endif

    if (abc) then
    if (MOD(N-5,8)/=0) then
    STOP 'ERROR: for type ABC star chain N needs to be odd and N-5 needs to be divisible by 8.'
    endif
    endif

endif

if ((max_eof).and.(single)) then
STOP 'ERROR: If you want to calculate the maximum EOF over a full window you'&
&'cannot do a single point calculation of the dynamics. Set *single* to False.'
endif

!if (output .and. .not. use_genetic) write(*,*) '>> Initial checks'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DEFINING BASIS VECTORS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate number of vectors for each excitation N!/exno!(N-exno)! subspace and the total number
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Allocate the size array
allocate(vectors_ex(exno))

! Determine the sizes
vectorstotal = 0
do i = 1, exno
    vectors_ex(i) = 1
    do j = N+1-i, N
        vectors_ex(i) = vectors_ex(i) * j
    end do
    vectors_ex(i) = vectors_ex(i) / factorial(i)
    vectorstotal = vectorstotal + vectors_ex(i)
end do

! Allocate matrices that will contain all the vectors
allocate(HT(vectorstotal,N))
HT = 0

! Create the subsectors matrices through a recursive call to Permutations
vectors_so_far = 1
do i = 1, exno

    ! Allocate the array to hold the permutations
    allocate(HS(vectors_ex(i),N))

    ! Get the permutations
    nit=1
    Ninit=1
    vec=0
    k=1
    ex=i
    call permutations(ex,nit,vec,N,Ninit,Hs,vectors_ex(i),k)

    ! Update the array of all basis vectors
    HT(vectors_so_far:vectors_so_far+vectors_ex(i)-1,:) = HS
    vectors_so_far = vectors_so_far + vectors_ex(i)

    ! Prevent memory leaks
    deallocate(HS)

end do

! Now we know the vectors, figure out which one AB+C (etc.) corresponds to
if (custom) then
    call get_vector_indices(HT)
end if
if (.not. custom) call initialState(initialVec)
if (.not. use_genetic) then
    do j=1,numModes
        write(40,"(A,I0,A)") "FOR MODE ", j, ":"
        do i=1,numI
            write(40,'(A,i4.1)') "  INITIAL INJECTED VECTOR INDEX = ", initialVec(j,i)
            write(40,'(A,f6.2,A,f6.2,A)') "  WITH COEFFICIENT = (", real(initialCoeff(j,i)), ", ", aimag(initialCoeff(j,i)), ")"
        enddo
        do i=1,numF
            write(40,"(A,i4.1)") "  FINAL VECTOR INDEX = ", targetVec(j,i)
            write(40,'(A,f6.2,A,f6.2,A)') "  WITH COEFFICIENT = (", real(targetCoeff(j,i)), ", ", aimag(targetCoeff(j,i)), ")"
        enddo
    end do

end if
!Stdout vectors matrix
if (output .and. .not. use_genetic) then
    200 FORMAT (/A)
    write(40,FMT=200) 'BASIS VECTORS:'
    do i=1,vectorstotal
        write(40,"(I0,A,500I3)") i,' -->',(HT(i,j),j=1,N)
    enddo
endif

!if (output .and. .not. use_genetic) write(*,*) '>> Basis vectors defined'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! INITIAL STATE NORMALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!if (output .and. .not. use_genetic) write(*,*) '>> Defining initial injection'

!normalization factor depending
!on the number of initial injections

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DEFINE CONNECTIVITY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!initialize
limits = 0
hub = 0
len_branch = 0

if (star .and. .not. custom) then

len_branch = ((2*(N - 1)/branches)+1)
hub = ((len_branch-1)/2) + 1

limits(1) = len_branch
limits(2) = len_branch + 1

do i=3,branches-1
    limits(i) = limits(i-1) + (len_branch-2)
enddo

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DEFINE COUPLING PATTERN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (custom) then

    ! Generate the 2D coupling matrix based on the custom string
    call couplings_custom(Js2D, N, custom_string)
    if (XY_tony) then
        Js2D = Js2D*0.5
    end if
    ! Output the coupling pattern
    if (output .and. .not. use_genetic) then
        write(40,"()") 
        write(40,"(A)") "COUPLING PATTERN:"
        do i = 1, N
            do j = 1, N
                write(40, "(f7.3,A)", advance = "no") Js2D(i,j), " "
            enddo
            write(40, "()") 
        enddo
    endif

else

    call couplings(Js,len_branch,hub,limits)
    if (XY_tony) then
        Js = Js*0.5
    end if
    !Stdout coupling pattern
    301 FORMAT ("(spin",I4,")-(spin",I4,") -->",F6.2)
    if (output .and. .not. use_genetic) then
        write(40,FMT=200) 'COUPLING PATTERN:'
        do i=1,N-1
            write(40,FMT=301) i, i+1, Js(i)
        enddo
    endif

end if

!if (output .and. .not. use_genetic) write(*,*) '>> Coupling pattern defined'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ADD PERTURBATION FACTORS TO THE COUPLINGS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

seed = 0
call noise_sub1 (seed)

if (random_J) then
    do i=1,N-1
        r=algor_uniform_random()
        Js(i)=Js(i)+(r*E_J*J_max)
    enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! BUILD HAMILTONIAN IN THE SPIN BASIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(hami(vectorstotal,vectorstotal))

hami=0.0_dbl
if (custom) then
    call build_hamiltonian_custom(N,vectorstotal,Js2D,HT,hami,eeScale)
else if (linear) then
    call build_hamiltonian_linear(HT,Js,N,vectorstotal,hami)
else if (star) then
    call build_hamiltonian_star(HT,Js,N,vectorstotal,hami,branches,limits,hub,len_branch)
endif

!Stdout Hamiltonian
if (output .and. .not. use_genetic) then
    write(40,FMT=200) 'HAMILTONIAN MATRIX:'
    do i=1,vectorstotal
        write(40,"(30(I0,A))") (ceiling(hami(i,j))," ",j=1,vectorstotal)
    enddo
endif

!if (output .and. .not. use_genetic) write(*,*) '>> Hamiltonian Built'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ADD PERTURBATION FACTORS TO THE DIAGONAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(Noise(vectorstotal))

Noise=0
if (random_D) then
    do i=2,vectorstotal
    w = 0
        do j=1,N
            if (HT(i,j).eq.1) then
                w=w+1 !w counts the number of excitations in a vector
            endif
        enddo

    if (w.eq.1) then
        r=algor_uniform_random()
        Noise(i)=(E_D*r*J_max)
    else
        do j=2,N+1
            do k=1,N
                if ((HT(i,k).eq.1).and.(HT(i,k).eq.HT(j,k))) then
                    Noise(i)=Noise(i)+Noise(j)
                endif
            enddo
        enddo
    endif
    hami(i,i) = hami(i,i) + Noise(i)
    enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! TRANSLATE THE HAMILTONIAN IN THE MJ BASIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DIAGONALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(hamiD(vectorstotal,vectorstotal))

!make copy of Hamiltonian and parse it into complex matrix
hamiD=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)
do i=1,vectorstotal
do j=1,vectorstotal
hamiD(i,j)=(hami(i,j))
enddo
enddo

if (diagonalisation) then

allocate(eigvals(vectorstotal))
allocate(rwork((2*(vectorstotal**2))+5*vectorstotal+1))
allocate(work((vectorstotal**2)+2*vectorstotal))

liwork=5*vectorstotal+3

allocate(iwork(liwork))

!save Hamiltonian matrix
if (files .and. .not. use_genetic) then
    open(unit=89,file='hami.data',status='unknown')
    write(89,*) '#HAMILTONIAN'
    do i=1,vectorstotal
        write(89,*) (hami(i,j),j=1,vectorstotal)
    enddo
    close(89)
endif


!***********************************************************************************************************************
! LAPACK SUBROUTINE ©                                                                                                 !*
! ZHEEV computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix                          !*
!call zheev('V','U',vectorstotal,hamiD,size(hamiD,1),eigvals,work,size(work,1),rwork,info)                            !*
                                                                                                                      !*
call zheevd('V','U',vectorstotal,hamiD,size(hamiD,1),eigvals,work,size(work,1),rwork,size(rwork,1),iwork,liwork,info) !*
if(info/=0) stop 'ERROR in ZHEEV diagonalization'                                                                     !*
                                                                                                                      !*
!***********************************************************************************************************************

!check normalisation eigenvectors
do i=1,vectorstotal
    normal=0.
    do j=1,vectorstotal
        normal=normal+abs(hamiD(i,j))**2
    enddo
    if (abs(1.-normal)>=error) then
        print*, 'ERROR: your eigenvectors are not well normalized'
        STOP
    endif
enddo

!check eigenvectors orthogonality
do v=1,vectorstotal
    do i=1,vectorstotal
        orto=0.
        do j=1,vectorstotal
            orto=orto+real(hamiD(i,j)*dconjg(hamiD(v,j)))
        enddo
        if ((orto>error).and.(v/=i)) then
            print*, 'ERROR: your eigenvectors are not orthogonal'
            STOP
        endif
    enddo
enddo

!!Stdout Eigenvalues
if (output .and. .not. use_genetic) then

    !set formats
    write(tmp,'(i4.1)') vectorstotal
    fmt1='(1X,i4.1,1X,'//tmp//'("(",f7.3,f7.3,")"))'

    !print the headins
    !if N is small,
    !otherwise the
    !file gets too
    !messy
    if ((N.le.20).and.(exno.eq.1)) then
    fmt2='(6X,'
        do i=1,vectorstotal
            write(tmp,'(i4.1)') i
            fmt2=trim(fmt2)//'"Eigenvector'//trim(adjustl(tmp))//':",3X,'
        enddo

        fmt2=trim(fmt2)//")"

        !Eigenvalues
        write(40,FMT=200) 'EIGENVALUES:'
        do i=1,vectorstotal
            write(40,*) eigvals(i)
        enddo

        !Eigenvectors
        write(40,FMT=200) 'EIGENVECTORS'

        if (N.le.20) then
        write(40,fmt2)
        endif
        do i=1,vectorstotal
            write(40,fmt1) i ,(hamiD(i,:))
        enddo
    endif

endif

!!Save data in files
if (files .and. .not. use_genetic) then

    open (unit=41,file='coefficients.data',status='unknown')
    open (unit=42,file='probabilities.data',status='unknown')
    open (unit=43,file='eigenvalues.data',status='unknown')
    open (unit=44,file='maximumprob.data',status='unknown')

    write(41,*) '#COEFFICIENTS. EACH EIGENVECTORS IS A COLUMN. COLUMNS ARE ORDERED BY INCREASING VALUE OF ENERGY.'&
&'ROWS ARE ORDERED BY SITE BASIS VECTORS. ENERGY OF |00..0> IS INCLUDED'
    write(42,*) '#MODULUS SQUARE OF THE COEFFICIENTS. EACH EIGENVECTORS IS A COLUMN. COLUMNS ARE ORDERED'&
&'BY INCREASING VALUE OF ENERGY. ROWS ARE ORDERED BY SITE BASIS VECTORS. VECTOR |00..0> IS INCLUDED'
    write(43,*) '#EIGENVALUES. ENERGY OF |00..0>, E=0, IS INCLUDED'
    write(44,*) '#MAXIMUM PROBABILITIES PER ORDERED SITE BASIS VECTOR FOR THE EIGENVECTORS. FIRST ONE IS THE |00..0> STATE.'


    do i=1,vectorstotal
        write(41,*) real(hamiD(i,:))
        write(42,*) (dble(dconjg(hamiD(i,:))*(hamiD(i,:))))
        write(44,*) MAXVAL(dble((dconjg(hamiD(i,:))*(hamiD(i,:)))))
        write(43,*) eigvals(i)
    enddo

    close(41)
    close(42)
    close(43)
    

endif


if (output .and. .not. use_genetic) then
    !write(*,*) '>> Hamiltonian Diagonalization'
    flush(40)
end if

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! INTEGRATION METHOD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(c_i(vectorstotal))

c_i=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)

if (integration) then

    call time_integration(HT,hamiD,vectorstotal,c_i,modeNum)
    !if (output .and. .not. use_genetic) write(*,*) '>> Dynamics (direct integration method)'

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! DYNAMICS: FIDELITY, ENTROPY AND EOF !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (diagonalisation) then

    if (dynamics) then

        if (present(c_vs_t)) then        
            call injection_dynamics(HT,hamiD,eigvals,vectorstotal,c_i,initialtime,modeNum,c_vs_t)
        else
            call injection_dynamics(HT,hamiD,eigvals,vectorstotal,c_i,initialtime,modeNum)
        end if

        !if (output .and. .not. use_genetic) write(*,*) '>> Dynamics'

    end if

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! CLACULATE EOF AT CERTAIN TIME !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!red_rho=cmplx(0.0_dbl, 0.0_dbl, kind=dbl)
!
!call reduced_density_matrix(HT,vectorstotal,red_rho,c_i)
!
!!!Save reduced density matrix
!if (files) then
!    open(unit=89,file='reduced_rho.data',status='unknown')
!    write(89,*) '#REDUCED DENSITY MATRIX.'
!    do i=1,4
!        write(89,*) (red_rho(i,j),j=1,4)
!    enddo
!    close(89)
!endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! ENTROPY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! PLOTTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Writes in a file data needed for plots
if (.not. use_genetic) then
    open(unit=46,file='info.data',status='unknown')


    write(46, "(A,L)") "GENETIC=", use_genetic

    401 FORMAT ("GRAPHICAL=",L)
    write(46,401) graphical

    402 FORMAT ("REALISATIONS=",A)
    tmp = '1'
    if ((random_D).or.(random_J)) then
    write(tmp,'(i5.4)') num_realisations
    endif
    write(46,402) adjustl(trim(tmp))

    403 FORMAT ("N=",A)
    write(tmp,'(i5.4)') N
    write(46,403) adjustl(trim(tmp))

    404 FORMAT ("VECTORS=",A)
    write(tmp,'(i5.4)') vectorstotal
    write(46,404) adjustl(trim(tmp))

    405 FORMAT ("TOTALTIME=",A)
    write(tmp,'(f8.2)') initialtime+totaltime
    write(46,405) adjustl(trim(tmp))

    406 FORMAT ("TA=",A)
    write(tmp,'(f8.2)') t_A
    write(46,406) adjustl(trim(tmp))

    do j=1, numModes
    407 FORMAT ("INITIALVEC=",A)
    do i=1,numI
    write(tmp,'(i5.2)') initialVec(j,i)
    write(46,407) adjustl(trim(tmp))
    enddo
    enddo

    408 FORMAT ("EOF=",L)
    write(46,408) eof

    409 FORMAT ("SINGLE=",L)
    write(46,409) single

    410 FORMAT ("OFFNOISE=",A)
    write(tmp,'(f8.2)') E_J
    write(46,410) adjustl(trim(tmp))

    411 FORMAT ("DIAGNOISE=",A)
    write(tmp,'(f8.2)') E_D
    write(46,411) adjustl(trim(tmp))

    412 FORMAT ("MAX_EOF=",L)
    write(46,412) max_eof

    413 FORMAT ("METHOD=",A)
    if (integration) then
    tmp = 'INT'
    else if (diagonalisation) then
    tmp = 'DIAG'
    end if
    write(46,413) adjustl(trim(tmp))
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! FREE SPACE AND CLOSE FILES AND CLEAN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (output .and. .not. use_genetic) then
    close(unit=40)
endif

if (graphical) then
    close(unit=46)
endif

!!deallocation **VERY IMPORTANT**
if (allocated(H1))    deallocate(H1)
if (allocated(H3))    deallocate(H2)
if (allocated(H3))    deallocate(H3)
if (allocated(HT))    deallocate(HT)
if (allocated(hami))  deallocate(hami)
if (allocated(Noise)) deallocate(Noise)
if (allocated(c_i))   deallocate(c_i)
if (allocated(hamiD)) deallocate(hamiD)
if (allocated(Js))    deallocate(Js)
if (allocated(Js2D))  deallocate(Js2D)
if (allocated(vec))   deallocate(vec)

if (diagonalisation) then
    if (allocated(eigvals)) deallocate(eigvals)
    if (allocated(rwork))   deallocate(rwork)
    if (allocated(work))    deallocate(work)
    if (allocated(iwork))   deallocate(iwork)
endif

end subroutine

!*********************TO BE REPROGRAMED**********************!

subroutine noise_sub1 (seed)

use noisemod
implicit none
integer :: seed             !may or may not be 32 bit...
!32 bit integer version of seed
integer(kind=bit_32) :: iseed              !random number seed

!f90 intrinsic time call
character(len=10) :: system_time           !length is crucial ...
real (kind=dubp)    :: rtime

if (bit_size(iseed)/=32) call io_abort('Error in algor_set_random_seed - 32-bit integer kind type problem')

if (seed == 0 ) then
call date_and_time(time=system_time)  !character string hhmmss.xxx
read (system_time,*) rtime            !convert to real
rtime = rtime * 1000.0_dubp             !0<rtime<235959999.0 which fits within huge(1)
else
rtime = real(abs(seed),kind=dubp)          !convert seed to real
end if

!and then convert to bit_32 size integer
iseed = int(rtime,kind=bit_32)              !must fit within huge(1)

ix=ieor(777755555_bit_32,iseed)                   !Marsaglia generator
iy=ior(ieor(888889999_bit_32,iseed),1_bit_32)     !Parks-Miller generator

seed = int(rtime)                                 !return the seed that was used in default integer


return
end subroutine noise_sub1

function factorial(a)

    integer, intent(in) :: a
    integer :: factorial, i

    factorial = 1
    do i = 2, a
        factorial = factorial * i
    end do

end function

function algor_uniform_random()
!=========================================================================!
! Return a single random deviate ~ uniform [0,1] or [-1,1] or [-0.5,0.5]  !
! Based on Park-Miller "minimal standard" generator with Schrage's method !
!  to do 32-bit multiplication without requiring higher precision, plus   !
!  additional Marsaglia shift to suppress any weaknesses & correlations.  !
! Using two independent methods greatly increases the period of the       !
!  generator, s.t. resulting period ~2*10^18                              !
! NB Routine is only set to work with 32 bit integers!                    !
!-------------------------------------------------------------------------!
! References:                                                             !
!   S.K. Park and K.W. Miller, Commun. ACM, 31, p1192-1201 (1988)         !
!   L. Schrage, ACM Trans. Math. Soft., 5, p132-138 (1979)                !
!   G. Marsaglia, Linear Algebra and its Applications, 67, p147-156 (1985)!
!-------------------------------------------------------------------------!
! Return value:                                                           !
!   algor_uniform_random => required random deviate                       !
!-------------------------------------------------------------------------!
! Parent module variables used:                                           !
!   ix as next 32-bit integer in Marsaglia generator (updated)            !
!   iy as next 32-bit integer in Park-Miller generator (updated)          !
!-------------------------------------------------------------------------!
! Written by Matt Probert, v0.1, 01/07/2000                               !
!=========================================================================!
use noisemod
implicit none
real(kind=dubp)                 :: algor_uniform_random

!NB We use 3 logical XOR operations to do Marsaglia shift
!=> hard-wire 3 shift values (not all triplets any good)
!=> entire routine preset to only work with 32 bit integers.

!local variables ...
integer(kind=bit_32)            :: iy_tmp       !working value to force integer division
integer(kind=bit_32), parameter :: iy_max=2147483647 !2^31-1
real(kind=dubp), parameter        :: inv_iy_max=1.0_dubp/2147483647.0_dubp
real(kind=dubp)                   :: notright

!do Marsaglia shift sequence, period 2^32-1, to get new ix
ix=ieor(ix,ishft(ix, 13_bit_32))
ix=ieor(ix,ishft(ix,-17_bit_32))
ix=ieor(ix,ishft(ix,  5_bit_32))

!Park-Miller sequence, period iy_max-1, to get new iy
iy_tmp=iy/127773_bit_32                         !NB integer division
iy=16807_bit_32*(iy-iy_tmp*127773_bit_32)-iy_tmp*2836_bit_32  !next value of iy
if (iy < 0_bit_32) iy=iy+iy_max                 !integer modulo iy_max

!Combine ix and iy to get new random number, rescale onto range [0,1]
notright=inv_iy_max*ior(iand(iy_max,ieor(ix,iy)),1_bit_32)       !with masking to ensure non-zero

!Between -1 to 1
!algor_uniform_random=(notright*2.0_dubp)-1_dubp when {0-1}: !notright

!Between -0.5 and 0.5
algor_uniform_random=(notright*1.0_dubp)-0.5_dubp!notright

!Between 0 and 1
!algor_uniform_random=notright

return
end function algor_uniform_random

!**********************************************************************************!

end module 
