module genetic

    include 'mpif.h'

contains 

!=========================================================================!
! Subroutine that solves the system given by the parameters using         !
! a genetic algorithm                                                     !
!-------------------------------------------------------------------------!
! Written by Luke Mortimer, October 2019                                  !
!=========================================================================!

subroutine solve_genetic()

    use main
    use constants
    use parameters
    use dependencies

    implicit none

    ! For storing the initial network its corresponding fitness
    real(kind=dbl)             :: init_fit

    ! Save any bits removed from the given genome, to be added on after
    character(max_string_size) :: pos_direct, init_direct

    ! The array of different genomes that are evolving
    character(max_string_size), dimension(genomes_per_generation) :: population

    ! The array of fitnesses for the corresponding population
    real (kind=dbl),            dimension(genomes_per_generation) :: fitnesses

    ! The temporary array generated using population, before becoming population
    character(max_string_size), dimension(genomes_per_generation) :: population_new

    ! Temp vars used to return the dynamics when creating an animation
    integer :: bestIndex
    real (kind=dbl) :: bestFit
    real (kind=dbl), dimension(steps) :: bestDynams

    ! For keeping track of the best genomes throughout in case an earlier is the best
    real (kind=dbl) :: bestFitThisGen, bestFitOverall
    character(max_string_size) :: bestGenomeThisGen, bestGenomeOverall

    ! Unit and file to write to
    integer, parameter :: geneticFile = 35, animationFile = 36

    ! Loop counters
    integer         :: i, j

    ! Which genomes should be used to generate the next child
    integer         :: parent1, parent2

    ! Store the sum of the fitnesses to save recalculating it
    real (kind=dbl) :: fitness_sum

    ! For getting the final best fidelity and transfer time
    real (kind=dbl), dimension(:), allocatable :: best_fid, best_time, best_qual

    ! For MPI
    integer         :: proc_id, num_proc, mpi_error

    ! True if this node is the root node, saved to reduce proc_id == 0 checks
    logical         :: on_root_node

    ! Keep track of timings
    real            :: start_time, current_time
    integer         :: timeRequired

    ! Where this node should start and stop the processing of population/fitness
    integer         :: array_start, array_end, array_diff

    ! Keep track of whether the mutation amount has been reduced
    logical         :: has_changed_mutate

    ! The seed used for the random number generator
    integer         :: seed

    ! Hold the change in mutate rate as a real
    real(kind=dbl)  :: mutate_amount_real, mutate_delta

    if (.not. stop_after_time) then

        ! Open the output file
        open(unit=geneticFile, file="genetic.out")

        ! Open the file for gif production
        if (create_animation) then
            open(unit=animationFile, file="anim.data")
            write(animationFile, "(f7.3)") totalTime
            write(animationFile, "(i0)") steps
        end if

    end if

    ! Get rid of any <> or # in the genome string
    call process_directives(custom_string, init_direct, pos_direct)

    ! If using multiple nodes and not told to force
    if (.not. force_fitness .and. numModes > 1) then
        linear_fitness = .true.
    else if (.not. force_fitness) then
        linear_fitness = .false.
    end if

    ! Allocate the arrays returning the fidelity for each mode
    allocate(best_fid(numModes))
    allocate(best_time(numModes))
    allocate(best_qual(numModes))

    ! Start the parallelisation
    call MPI_INIT(mpi_error)

    ! Get the rank of the current processor
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_proc, mpi_error)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, mpi_error)

    ! Stop if the amount per genome isn't a multiple of the number of cores
    if (modulo(genomes_per_generation, num_proc) /= 0) then

        ! Stop the parallelisation
        call MPI_FINALIZE(mpi_error)

        ! Stop the program
        stop "Number of genomes per generation must be multiple of number of cores"

    end if

    ! Is this the root (i.e. first, main) node?
    on_root_node = (proc_id == 0)

    ! Where this node should operate in the population/fitness arrays
    array_start = nint((real(genomes_per_generation) / real(num_proc)) * real(proc_id) + 1.0)
    array_end = nint((real(genomes_per_generation) / real(num_proc)) * (real(proc_id) + 1.0))
    array_diff = array_end - array_start + 1

    ! Set the initial mutation amount and the change per iteration
    mutate_amount_real = real(mutate_amount_initial, dbl)
    mutate_amount = int(mutate_amount_real)
    mutate_delta = real(mutate_amount_initial - mutate_amount_final, dbl) / real(max_generations, dbl)

    ! Get the initial time
    call cpu_time(start_time)

    ! Make sure the random number generator is seeded differently per processor
    seed = int(real(start_time, dbl) + 1000.0_dbl * real(proc_id, dbl))
    call noise_sub1(seed)

    ! Only do outputs if on the root node
    if (on_root_node) then

        if (.not. stop_after_time) then

            ! Initial output
            if (use_genetic) write(6, "(A)") "Running first generation for a time estimate"

            write(geneticFile, "(A)") "|------------------------------------------------------|"
            write(geneticFile, "(A)") "|                     PARAMETERS                       |"
            write(geneticFile, "(A)") "|------------------------------------------------------|"
            write(geneticFile, "(A)") ""
            write(geneticFile, "(A,I0)")   "cores                  = ", num_proc
            write(geneticFile, "(A,A)")    "pos directive          = ", trim(pos_direct)
            write(geneticFile, "(A,A)")    "<i|f> directive        = ", trim(init_direct)
            write(geneticFile, "(A,A)")    "stripped genome        = ", trim(custom_string)
            write(geneticFile, "(A,I0)")   "generations            = ", max_generations
            write(geneticFile, "(A,I0)")   "genomes per generation = ", genomes_per_generation
            write(geneticFile, "(A,f5.3)") "mutate chance          = ", mutate_chance
            write(geneticFile, "(A,I0)")   "intial mutate amount   = ", mutate_amount_initial
            write(geneticFile, "(A,I0)")   "final mutate amount    = ", mutate_amount_final
            write(geneticFile, "(A,L)")    "use linear scaling     = ", linear_fitness
            write(geneticFile, "(A,L)")    "only change energies   = ", energies_only
            write(geneticFile, "(A)") ""

            ! Solve the initial system as a reference point and to load initialVec
            init_fit = genetic_fitness(custom_string, best_fid, best_time, best_qual, numModes)

            write(geneticFile, "(A)") ""
            write(geneticFile, "(A)") "|------------------------------------------------------|"
            write(geneticFile, "(A)") "|                 INITIAL SYSTEM                       |"
            write(geneticFile, "(A)") "|------------------------------------------------------|"
            write(geneticFile, "(A)") ""

            write(geneticFile, "(A,A)") "initial genome: ", '"' // trim(init_direct) // trim(custom_string) // &
                                    & trim(pos_direct) // '"'
            write(geneticFile, "(A,f0.2,A)", advance="no") "with fitness: ", init_fit, " ("
            do j = 1, numModes
                write(geneticFile, "(f0.2,A,f0.2,A,f0.2)", advance="no") best_fid(j)*100.0_dbl, "% fidelity ", &
                    & best_qual(j)*100.0_dbl, "% quality at time ", best_time(j)
                if (j /= numModes) write(geneticFile, "(A)",advance="no") ", "
            end do
            write(geneticFile, "(A)") ")"
            write(geneticFile, "(A)") ""

            ! If just running once for fitness, stop here
            if (.not. use_genetic) then

                ! Ensure the file is written as it goes
                flush(geneticFile)

                ! If told to do a single calc, do a dynamics anyway for the sake of graphs
                single = .false.
                call solve(1)

                ! Stop the parallelisation
                call MPI_FINALIZE(mpi_error)

                ! Stop the program
                stop 

            end if

            write(geneticFile, "(A)") "|------------------------------------------------------|"
            write(geneticFile, "(A)") "|                   SYSTEM EVOLUTION                   |"
            write(geneticFile, "(A)") "|------------------------------------------------------|"
            write(geneticFile, "(A)") ""

            write(geneticFile, "(A)") "|------------|--------|---------|--------|-------------|--------|"
            write(geneticFile, "(A)") "| Generation | Worst  | Average |  Best  | Elapsed / s | Mutate |"
            write(geneticFile, "(A)") "|------------|--------|---------|--------|-------------|--------|"

            ! Ensure the file is written as it goes
            flush(geneticFile)

        end if

    end if

    ! Create the initial pool of genomes 
    do i = array_start, array_end

        population(i) = custom_string 

        ! Mutate a certain number of times
        do j = 1, num_initial_mutations
            call genetic_mutate(population(i))
        end do

    end do

    ! Sync the population between nodes
    if (num_proc > 1) then
        call sync_pop(population, array_start, array_end, array_diff)
    end if

    ! Mutatation amount hasn't been changed yet
    has_changed_mutate = .false.

    ! Keep looping until max generators reached
    do i = 1, max_generations-1

        ! Evaluate all of the genomes to get their fitnesses
        do j = array_start, array_end

            fitnesses(j) = genetic_fitness(population(j), best_fid, best_time, best_qual, numModes) 
        end do

        ! Sync the fitnesses between nodes
        if (num_proc > 1) then
            call sync_fit(fitnesses, array_start, array_end, array_diff)
        end if

        ! Sum the fitnesses
        fitness_sum = sum(fitnesses)

        ! Only output info once
        if (on_root_node) then

            ! Get the change in time
            call cpu_time(current_time)

            ! Output info about the current genome
            if (.not. stop_after_time) then

                ! Get the best fitness and update if the best overall 
                bestIndex = maxloc(fitnesses, 1)
                bestFitThisGen = fitnesses(bestIndex)
                bestGenomeThisGen = population(bestIndex)
                if (bestFitThisGen > bestFitOverall) then
                    bestFitOverall = bestFitThisGen
                    bestGenomeOverall = bestGenomeThisGen
                end if

                ! Genetic output
                write(geneticFile, "(A,I4,A,f6.2,A,f6.2,A,f6.2,A,f8.1,A,I6,A)") "|    ", &
                        & i, "    | ", minval(fitnesses), " | ", &
                        & fitness_sum / genomes_per_generation, "  | ", maxval(fitnesses), " |    ", &
                        & current_time-start_time, " | ", mutate_amount, " | "

                ! Output for animation 
                if (create_animation) then
                    
                    bestFit = genetic_fitness(population(bestIndex), best_fid, best_time, best_qual, numModes, bestDynams)
                    write(animationFile, "(A,A,100(f10.5))") trim(population(bestIndex)), trim(pos_direct), bestDynams

                end if

            end if

            ! After the first generation, give a time estimate in a nice human readable form
            if (i == 1) then

                timeRequired = int(real(max_generations)*(current_time-start_time))

                ! If requested, output the time it would take in seconds and then stop
                if (stop_after_time) then

                    write(6, "(I0)") timeRequired

                    ! Stop the parallelisation
                    call MPI_FINALIZE(mpi_error)

                    ! Stop the program
                    stop 

                else
                    write(6, "(A,A,A,A,A)") "Should be finished in ", &
                        & trim(seconds_to_human(timeRequired)), " (at ", trim(finish_time(timeRequired)), ")"
                end if

            end if

            ! Ensure the file is written as it goes
            flush(geneticFile)

        end if

        ! If told to stop early
        if (i == 1 .and. stop_after_time) then

            ! Stop the parallelisation
            call MPI_FINALIZE(mpi_error)

            ! Stop the program
            stop 

        end if

        ! If the maximum fitness goes above a certain value, stop (can be set > 100 to never)
        if (maxval(fitnesses) > stop_after_fit) then
            if (on_root_node) then
                write(6, "(A)") ""
                write(6, "(A,f0.1,A)") "Stopping early since required fitness reached (", stop_after_fit, ")"
            end if
            exit
        end if

        ! Keep creating new genomes until the generation size is reached
        do j = array_start, array_end

            ! Choose the parents
            parent1 = genetic_select(fitnesses, fitness_sum)
            parent2 = genetic_select(fitnesses, fitness_sum)

            ! Crossover the genomes
            population_new(j) = genetic_crossover(population(parent1), population(parent2))

            ! With a certain chance, apply a mutation
            if (algor_uniform_random() < mutate_chance) then
                call genetic_mutate(population_new(j))
            end if

        end do

        ! The new generation becomes the main generation
        do j = 1, genomes_per_generation
            population(j) = population_new(j)
        end do

        ! Sync the population between nodes
        if (num_proc > 1) then
            call sync_pop(population, array_start, array_end, array_diff)
        end if

        ! Lower the mutation rate
        mutate_amount_real = mutate_amount_real - mutate_delta
        mutate_amount = int(mutate_amount_real)

    end do

    ! Evaluation the final generation
    do j = array_start, array_end
        fitnesses(j) = genetic_fitness(population(j), best_fid, best_time, best_qual, numModes) 
    end do

    ! Sync the population between nodes
    if (num_proc > 1) then
        call sync_fit(fitnesses, array_start, array_end, array_diff)
    end if

    ! Only do the final output once
    if (on_root_node) then

        ! Get the change in time
        call cpu_time(current_time)

        ! Sum the fitnesses
        fitness_sum = sum(fitnesses)

        i = max_generations

        ! Output for the final generation
        write(geneticFile, "(A,I4,A,f6.2,A,f6.2,A,f6.2,A,f8.1,A,I6,A)") "|    ", i, "    | ", minval(fitnesses), " | ", &
                    & fitness_sum / genomes_per_generation, "  | ", maxval(fitnesses), " |    ", &
                    & current_time-start_time, " | ", 0, " | "
        write(geneticFile, "(A)") "|------------|--------|---------|--------|-------------|--------|"

        write(geneticFile, "(A)") ""
        write(geneticFile, "(A)") "|------------------------------------------------------|"
        write(geneticFile, "(A)") "|                     FINAL SYSTEM                     |"
        write(geneticFile, "(A)") "|------------------------------------------------------|"
        write(geneticFile, "(A)") ""

        ! Get the best fitness and update if the best overall 
        bestIndex = maxloc(fitnesses, 1)
        bestFitThisGen = fitnesses(bestIndex)
        bestGenomeThisGen = population(bestIndex)
        if (bestFitThisGen > bestFitOverall) then
            bestFitOverall = bestFitThisGen
            bestGenomeOverall = bestGenomeThisGen
        end if

        ! Output the final result
        bestFitOverall = genetic_fitness(bestGenomeOverall, best_fid, best_time, best_qual, numModes)
        write(geneticFile, "(A,A)") "best genome: ", &
                 & '"' // trim(init_direct) // trim(bestGenomeOverall) // trim(pos_direct) // '"'
        write(geneticFile, "(A,f0.2,A)", advance="no") "with fitness: ", bestFitOverall, " ("
        do j = 1, numModes
            write(geneticFile, "(f0.2,A,f0.2,A,f0.2)", advance="no") best_fid(j)*100.0_dbl, "% fidelity ", &
                & best_qual(j)*100.0_dbl, "% quality at time ", best_time(j)
            if (j /= numModes) write(geneticFile, "(A)",advance="no") ", "
        end do
        write(geneticFile, "(A)") ")"
        write(geneticFile, "(A)") ""

        ! Ensure the file is written as it goes
        flush(geneticFile)

        ! Run the final system as normal with full output
        use_genetic = .false.
        single = .false.
        call solve(1)

    end if

    ! Stop the parallelisation
    call MPI_FINALIZE(mpi_error)

    ! Close files
    close(geneticFile)
    if (create_animation) then
        close(animationFile)
    end if

end subroutine

!=========================================================================!
! Subroutine that solves a given system and returns a fitness score       !
!-------------------------------------------------------------------------!
! Parameters                                                              !
!   string      : the genome to use as the network                        !
!-------------------------------------------------------------------------!
! Written by Luke Mortimer, October 2019                                  !
!=========================================================================!

function genetic_fitness(string, fid, time, qual, modes, dynams)

    use constants
    use parameters
    use main

    ! Inputs
    character(*), intent(in) :: string
    integer, intent(in) :: modes

    ! Outputs
    real(kind=dbl) :: genetic_fitness
    real(kind=dbl), dimension(:), intent(out) :: fid, time, qual
    real(kind=dbl), dimension(steps), intent(out), optional :: dynams

    ! Array to hold the various fitnesses for each mode
    real(kind=dbl), dimension(modes) :: fitnesses

    ! Get the coefficients from dynamics to search for the maximum fidelity
    complex(kind=dbl), dimension(:, :), allocatable :: c_vs_t
    
    ! The array combining the various target vecs
    real(kind=dbl), dimension(steps) :: searchArray

    ! Vars for formula simplification
    complex(kind=dbl) :: a, b

    ! Save the fidelity and time until max fidelity for fitness calc
    real(kind=dbl) :: fidelity, transferTime, quality
    complex(kind=dbl) :: complexSum

    ! Determine the index at which the maximum fidelity occurs
    integer :: maxIndex, startIndex, endIndex

    ! Loop counters
    integer :: i, j, l

    ! Only used if animating
    logical :: wasSingle

    ! Set the start and end indices
    startIndex = int((min_time / totalTime) * steps) + 1
    endIndex = int((max_time / totalTime) * steps)

    ! Set the current network to this new string
    custom_string = string

    ! Init the fitness array
    fitnesses = 0.0_dbl
    fid = 0.0_dbl
    time = 0.0_dbl
    qual = 0.0_dbl

    ! If told to give full dynamics, do it
    if (present(dynams)) then
        wasSingle = single
        single = .false.
    end if

    ! Repeat for each mode of operation
    do j = 1, numModes

        ! Solve the system, returning the coefficients for each step
        call solve(j, c_vs_t)

        if (.not. single) then

            ! Combine the target vector coefficients to get the quality factor 
            searchArray = 0.0_dbl
            do l = 1, steps

                complexSum = cmplx(0.0, 0.0)

                ! Sum <b|a>
                do i = 1, numF
                
                    if (targetVec(j, i) == 0) cycle
                    a = c_vs_t(l, targetVec(j, i))
                    b = targetCoeff(j, i)
                    complexSum = complexSum + (a * conjg(b))

                end do

                ! Turn the complex sum into a fidelity |<b|a>|^2
                searchArray(l) = abs(complexSum) ** 2

            end do

            ! If told to, return the fidelity vs time array
            if (present(dynams)) then
                dynams = searchArray
            end if

            ! Find the location where the system is closest to the target state
            maxIndex = startIndex-1+maxloc(searchArray(startIndex:endIndex), 1)
            fidelity = searchArray(maxIndex)
            transferTime = (real(maxIndex, dbl) / real(steps, dbl)) * totalTime

            ! Get the quality at this point
            quality = 1.0_dbl
            do i = 1, numF

                if (targetVec(j, i) == 0) cycle

                a = c_vs_t(maxIndex, targetVec(j, i))
                b = targetCoeff(j, i)
                quality = quality * exp(-abs(atan2(aimag(a), real(a))&
                                                  & -atan2(aimag(b), real(b))))

            end do

        else

            quality = 1.0_dbl
            fidelity = 0.0_dbl
            complexSum = cmplx(0.0, 0.0)

            ! For each component of the target vector
            do i = 1, numF

                if (targetVec(j, i) == 0) cycle
                a = c_vs_t(1, targetVec(j, i))
                b = targetCoeff(j, i)

                ! The quality and fidelity at that specific time
                quality = quality * exp(-abs(atan2(aimag(a), real(a))&
                                                  & -atan2(aimag(b), real(b))))
                complexSum = complexSum + (a * conjg(b))

            end do

            ! Turn the complex sum into a fidelity |<b|a>|^2
            fidelity = abs(complexSum) ** 2

            ! Only have one value if done single point calc
            maxIndex = 1
            transferTime = t_A

        end if

        ! If given variables, put various quantities in them
        qual(j) = quality
        fid(j) = fidelity
        time(j) = transferTime

        ! The fitness from the quality
        if (linear_fitness) then
            fitnesses(j) = 100.0_dbl * fidelity
        else
            fitnesses(j) = 100.0_dbl * exp(fidelity_scale * ((fidelity) - 1.0_dbl))
        end if

        ! Determine the fitness contribution from the time
        if (minimise_time .and. .not. single) then
            fitnesses(j) = fitnesses(j) * exp(time_scale * transferTime)
        end if

    end do

    ! Revert
    if (present(dynams)) then
        single = wasSingle
    end if

    ! Average the fitnesses to get the overall fitness
    genetic_fitness = sum(fitnesses) / modes

end function

!=========================================================================!
! Choose a genome from the given array based on the fitnesses             !
!-------------------------------------------------------------------------!
! Parameters                                                              !
!   genomes      : the array of genomes                                   !
!   fitnesses    : the array of respective fitnesses                      !
!-------------------------------------------------------------------------!
! Written by Luke Mortimer, October 2019                                  !
!=========================================================================!

function genetic_select(fitnesses, fitness_sum)

    use constants
    use parameters
    use main

    ! Inputs
    real(kind=dbl), dimension(:), intent(in) :: fitnesses
    real(kind=dbl), intent(in)               :: fitness_sum

    ! Output
    integer :: genetic_select

    ! The sum which will stop the loop
    real(kind=dbl) :: target_sum

    ! The partial sum of the normalised fitnesses
    real(kind=dbl) :: partial_sum

    ! Loop counters
    integer :: i

    ! Generate the random position between 0 and 1 which will select the parent
    target_sum = (algor_uniform_random() + 0.5_dbl) * fitness_sum

    ! Sum the normalised fitnesses until the target is met, then exit
    partial_sum = 0.0_dbl
    do i = 1, genomes_per_generation
        partial_sum = partial_sum + fitnesses(i)
        if (partial_sum >= target_sum) exit
    end do

    ! Return the index which stopped the loop
    genetic_select = min(i, genomes_per_generation)

end function

!=========================================================================!
! Combine two genomes, returning the result                               !
!-------------------------------------------------------------------------!
! Parameters                                                              !
!   string1    : the first genome to use                                  !
!   string2    : the second genome to use                                 !
!-------------------------------------------------------------------------!
! Written by Luke Mortimer, October 2019                                  !
!=========================================================================!

function genetic_crossover(string1, string2)

    use parameters
    use constants
    use main

    ! Inputs
    character(*), intent(in) :: string1, string2

    ! Output
    character(len(string1)) :: genetic_crossover

    ! Where in the strings they should crossover
    integer :: crossoverIndex

    ! The length of the second (ideally both) string
    integer :: length

    ! Loop counters
    integer :: i

    ! Get the string's length
    length = len_trim(string2)

    ! Ensure the return value starts blank
    genetic_crossover = ""

    if (use_alternating) then

        ! For each character in the strings
        do i = 1, length

            ! With a 50% chance, choose the char from string 1
            if (algor_uniform_random() > 0.0_dbl) then
                genetic_crossover(i:i) = string1(i:i)
            else
                genetic_crossover(i:i) = string2(i:i)
            end if

        end do

    else

        ! Choose the crossover point linearly
        crossoverIndex = int((algor_uniform_random() + 0.5_dbl)*(length-1.0_dbl))+1

        ! Combine the strings at this point
        genetic_crossover(:) = string1(1:crossoverIndex) // string2(crossoverIndex+1:length)

    end if

end function

!=========================================================================!
! Mutate a genome in place                                                !
!-------------------------------------------------------------------------!
! Parameters                                                              !
!   string  : genome to mutate                                            !
!-------------------------------------------------------------------------!
! Written by Luke Mortimer, October 2019                                  !
!=========================================================================!

subroutine genetic_mutate(string)

    use main
    use parameters
    use constants

    ! Inputs (modified)
    character(*), intent(inout) :: string

    ! The index of the 2 digit integer to change
    integer :: numLoc

    ! The change to be applied to the chosen coupling
    integer :: numDelta

    ! Used to store the read coupling strength, which is then modified and rewritten
    integer :: couplingVal

    ! Stores a random number for formula simplification
    integer :: r

    ! Length of the string
    integer :: length

    ! Where the coupling started the mutation negative or not
    logical :: startedNegative

    ! The two coupling characters to write after updating the value
    character(10) :: format_string

    ! Save the length of the string to prevent further len_trim calls
    length = len_trim(string)

    ! How much should the coupling change by? -2, -1, 0, +1, +2
    numDelta = nint(algor_uniform_random() * real(mutate_amount, dbl) * 2.0_dbl)

    ! Pick a random coupling number from the string
    r = nint((algor_uniform_random() + 0.5_dbl) * ((length)/(2+coupling_digits)-1)+1)

    ! Determine where in the string the digits start
    numLoc = (2 + coupling_digits)*(r-1) + 3

    ! If told to only optimise energies, stop if about to change a coupling 
    if (energies_only .and. string(numLoc-2:numLoc-2) /= string(numLoc-1:numLoc-1)) then
        return
    end if

    ! Create a format statement
    write(format_string, "(A,I0,A,I0,A)") "(I", coupling_digits, ".", coupling_digits, ")"

    ! Read the val (can't do I0 internal reads)
    read(string(numLoc:numLoc+coupling_digits-1), format_string) couplingVal

    ! Make the value negative if the letters are the wrong way round 
    if (ichar(string(numLoc-2:numLoc-2)) > ichar(string(numLoc-1:numLoc-1))) then
        couplingVal = -couplingVal
        startedNegative = .true.
    else
        startedNegative = .false.
    end if

    ! Update the val
    if (allow_negative) then
        couplingVal = max(min(couplingVal + numDelta, max_val), -max_val)
    else
        couplingVal = max(min(couplingVal + numDelta, max_val), min_val)
    end if

    ! Don't ever fully remove a node, to prevent later errors
    if (couplingVal == 0) couplingVal = 1

    ! Switch letters if there is a change of sign
    if ((couplingVal < 0 .and. .not. startedNegative) .or. &
      & (couplingVal > 0 .and. startedNegative)) then
        string(numLoc-2:numLoc-1) = string(numLoc-1:numLoc-1) // string(numLoc-2:numLoc-2)
    end if

    ! Ensure the value written to the string is positive
    if (couplingVal < 0) then
        couplingVal = -couplingVal
    end if

    ! Write out the val
    write(string(numLoc:numLoc+coupling_digits-1), format_string) couplingVal

end subroutine

!=========================================================================!
! Sync the fitnesses amongst nodes                                        !
!-------------------------------------------------------------------------!
! Parameters                                                              !
!   fitnesses    : array to sync                                          !
!   array_start  : where in the fitness array to start syncing from       !
!   array_end    : where in the fitness array to stop syncing             !
!   array_diff   : the number of elements to sync                         !
!   on_root_node : whether this is the special node or not                !
!   num_proc     : the total number of nodes                              !
!-------------------------------------------------------------------------!
! Written by Luke Mortimer, October 2019                                  !
!=========================================================================!

subroutine sync_fit(fitnesses, array_start, array_end, array_diff)

    use constants
    use parameters

    implicit none

    real(kind=dbl), dimension(genomes_per_generation), intent(inout) :: fitnesses
    integer, intent(in) :: array_start, array_end, array_diff
    real(kind=dbl), dimension(array_diff) :: fitnessesLocal
    integer :: mpi_error

    fitnessesLocal = fitnesses(array_start:array_end)

    ! Gather the population from every node
    call MPI_ALLGATHER(fitnessesLocal(1), array_diff, MPI_DOUBLE_PRECISION, &
        & fitnesses(1), array_diff, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpi_error)

end subroutine

!=========================================================================!
! Sync the population amongst nodes                                       !
!-------------------------------------------------------------------------!
! Parameters                                                              !
!   population    : array to sync                                         !
!   array_start  : where in the fitness array to start syncing from       !
!   array_end    : where in the fitness array to stop syncing             !
!   array_diff   : the number of elements to sync                         !
!   on_root_node : whether this is the special node or not                !
!   num_proc     : the total number of nodes                              !
!-------------------------------------------------------------------------!
! Written by Luke Mortimer, October 2019                                  !
!=========================================================================!

subroutine sync_pop(population, array_start, array_end, array_diff)

    use parameters

    implicit none

    character(max_string_size), dimension(genomes_per_generation), intent(inout) :: population
    integer, intent(in) :: array_start, array_end, array_diff
    character(max_string_size), dimension(array_diff) :: popLocal
    integer :: mpi_error

    popLocal = population(array_start:array_end)

    ! Make sure all processors are synced here
    call MPI_BARRIER(MPI_COMM_WORLD, mpi_error)

    ! Gather the population from every node
    call MPI_ALLGATHER(popLocal(1), array_diff*max_string_size, MPI_CHARACTER, &
        & population(1), array_diff*max_string_size, &
        & MPI_CHARACTER, MPI_COMM_WORLD, mpi_error)

end subroutine

!=========================================================================!
! Convert some amount of seconds into hours, minutes, seconds etc.        !
!-------------------------------------------------------------------------!
! Parameters                                                              !
!   seconds    : number of seconds to convert                             !
!-------------------------------------------------------------------------!
! Written by Luke Mortimer, October 2019                                  !
!=========================================================================!

function seconds_to_human(seconds)

    integer, intent(in) :: seconds
    character(50) :: seconds_to_human
    integer :: num_hours, num_minutes, num_seconds

    ! Determine the various components
    num_hours = seconds / 3600
    num_minutes = modulo(seconds, 3600) / 60
    num_seconds = modulo(seconds, 60)

    ! Write it to the output string
    if (num_hours > 0) then
        write(seconds_to_human, "(I0,A,I0,A,I0,A)") num_hours, " hours ", num_minutes, " minutes ", num_seconds, " seconds"
    else if (num_minutes > 0) then
        write(seconds_to_human, "(I0,A,I0,A)") num_minutes, " minutes ", num_seconds, " seconds"
    else
        write(seconds_to_human, "(I0,A)") num_seconds, " seconds"
    end if

end function

function finish_time(seconds)

    integer, intent(in) :: seconds
    integer :: seconds_left
    character(50) :: finish_time, time_temp

    integer, dimension(8) :: values
    integer, dimension(8) :: delta_values

    finish_time = ""

    call date_and_time(values=values)
    seconds_left = seconds

    ! Year
    delta_values(1) = seconds_left / int(3.154e+7)
    seconds_left = modulo(seconds_left, int(3.154e+7))

    ! Month
    delta_values(2) = seconds_left / int(2.628e+6)
    seconds_left = modulo(seconds_left, int(2.628e+6))

    ! Day
    delta_values(3) = seconds_left / 86400
    seconds_left = modulo(seconds_left, 86400)

    ! Hour
    delta_values(5) = seconds_left / 3600
    seconds_left = modulo(seconds_left, 3600)

    ! Minute
    delta_values(6) = seconds_left / 60
    seconds_left = modulo(seconds_left, 60)

    values = values + delta_values

    ! Overflow minutes to hours
    if (values(6) >= 60) then
        values(5) = values(5) + 1
        values(6) = 0
    end if

    ! Overflow hours to days
    if (values(5) >= 24) then
        values(3) = values(3) + 1
        values(5) = 0
    end if

    ! Overflow days to months
    if (values(3) > get_days_in_month(values(2), values(1))) then
        values(2) = values(2) + 1
        values(3) = 0
    end if

    ! Overflow months to years
    if (values(2) > 12) then
        values(1) = values(1) + 1
        values(2) = 0
    end if

    time_temp = ""
    write(time_temp, "(I4,A,I2,A,I2,A,I2,A,I2)") values(1),"/",values(2),"/",values(3)," ",values(5),":",values(6)

    ! Remove any spaces apart from where there should be one
    do i = 1, len_trim(time_temp)
        if (i /= 11 .and. time_temp(i:i) == " ") time_temp(i:i) = "0"
    end do

    finish_time = time_temp

end function

! Returns the number of days in a month
function get_days_in_month(month, year)

    integer, intent(in) :: month, year
    integer :: get_days_in_month

    select case (month)

        case (1)
            get_days_in_month = 31
        case (2)
            if (modulo(year, 4) == 0) then
                get_days_in_month = 29
            else
                get_days_in_month = 28
            end if
        case (3)
            get_days_in_month = 31
        case (4)
            get_days_in_month = 30
        case (5)
            get_days_in_month = 31
        case (6)
            get_days_in_month = 30
        case (7)
            get_days_in_month = 31
        case (8)
            get_days_in_month = 31
        case (9)
            get_days_in_month = 30
        case (10)
            get_days_in_month = 31
        case (11)
            get_days_in_month = 30
        case (12)
            get_days_in_month = 31
        case default
            get_days_in_month = 0

    end select

end function

end module
