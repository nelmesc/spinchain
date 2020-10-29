program wrapper

    use constants
    use parameters
    use main
    use genetic
    use dependencies

    ! Vars for processing command args
    character(max_string_size), dimension(:), allocatable :: args
    integer                    :: num_args
    integer                    :: i, j, k
    integer                    :: convertTo = -1

    ! Get the amount of arguments passed to the program
    num_args = command_argument_count()
    allocate(args(num_args))

    ! Save the arguments into an array
    do i = 1, num_args

        ! Get the argument
        call get_command_argument(i, args(i))

    end do
    
    mutate_amount_initial = -1

    ! Process each arg
    i = 1
    do while (i <= num_args)

        ! Process the argument
        if (args(i) == "-o" .or. trim(args(i)) == "--optimise") then

            use_genetic = .true.

        else if (args(i) == "-f" .or. args(i) == "--fast") then

            time_scale = time_scale * 10000.0_dbl

        else if (args(i) == "-w" .or. args(i) == "--window") then

            if (.not. is_real(args(i+1))) then
                print *, "ERROR - the argument following " //trim(args(i))// " should be a float"
                stop
            end if

            totalTime = chars_to_real(args(i+1))

        else if (args(i) == "-e" .or. args(i) == "--eescale") then

            if (.not. is_real(args(i+1))) then
                print *, "ERROR - the argument following " //trim(args(i))// " should be a float"
                stop
            end if

            eeScale = chars_to_real(args(i+1))

        else if (args(i) == "-n" .or. args(i) == "--negative") then

            allow_negative = .true.

        else if (args(i) == "-F" .or. args(i) == "--fidelity") then

            if (.not. is_int(args(i+1))) then
                print *, "ERROR - the argument following " //trim(args(i))// " should be an integer"
                stop
            end if

            stop_after_fid = chars_to_int(args(i+1)) / 100.0

        else if (args(i) == "-g" .or. args(i) == "--genomes") then

            if (.not. is_int(args(i+1))) then
                print *, "ERROR - the argument following " //trim(args(i))// " should be an integer"
                stop
            end if

            genomes_per_generation = chars_to_int(args(i+1))

        else if (args(i) == "-G" .or. args(i) == "--generations") then

            if (.not. is_int(args(i+1))) then
                print *, "ERROR - the argument following " //trim(args(i))// " should be an integer"
                stop
            end if

            max_generations = chars_to_int(args(i+1))

        else if (args(i) == "-T" .or. args(i) == "--timefull") then

            stop_after_time_full = .true.
            use_genetic = .true.

        else if (args(i) == "-t" .or. args(i) == "--time") then

            stop_after_time = .true.
            use_genetic = .true.

        else if (args(i) == "-c" .or. args(i) == "--convert") then

            if (.not. is_int(args(i+1))) then
                print *, "ERROR - the argument following " //trim(args(i))// " should be an integer"
                stop
            end if

            convertTo = chars_to_int(args(i+1))

        else if (args(i) == "-M" .or. args(i) == "--mutate-end") then

            if (.not. is_int(args(i+1))) then
                print *, "ERROR - the argument following " //trim(args(i))// " should be an integer"
                stop
            end if

            mutate_amount_final = chars_to_int(args(i+1))

        else if (args(i) == "-m" .or. args(i) == "--mutate") then

            if (.not. is_int(args(i+1))) then
                print *, "ERROR - the argument following " //trim(args(i))// " should be an integer"
                stop
            end if

            mutate_amount_initial = chars_to_int(args(i+1))

        else if (args(i) == "-e" .or. args(i) == "--energy") then

            energies_only = .true.

        else if (args(i) == "-l" .or. args(i) == "--linear") then

            linear_fitness = .true.
            force_fitness = .true.

        else if (args(i) == "-x" .or. args(i) == "--exponential") then

            linear_fitness = .false.
            force_fitness = .true.

        else if (len_trim(args(i)) >= 5 .and. args(i)(1:1) /= "-" .and. .not. is_int(args(i))) then

            ! Clear the string first, just in case
            custom_string = ""

            ! Strip any spaces or newlines
            k = 1
            do j = 1, len_trim(args(i))
                if (args(i)(j:j) /= " " .and. args(i)(j:j) /= new_line("a")) then
                    custom_string(k:k) = args(i)(j:j)
                    k = k + 1
                end if
            end do

            custom = .true.

        else if (args(i) == "-C" .or. args(i) == "--cores") then

        else if (args(i) == "-s" .or. args(i) == "--scaling") then

        else if (args(i) == "-S" .or. args(i) == "--sites") then

        else if (args(i) == "-E" .or. args(i) == "--eetest") then

        else if (args(i) == "-a" .or. args(i) == "--anim") then

            create_animation = .true.

        else if (.not. is_int(args(i)) .and. .not. is_real(args(i))) then

            print *, "ERROR - unknown argument: ", trim(args(i))
            stop

        end if

        i = i + 1

    end do

    ! If told to search for longer, run for longer too
    if (max_time > totalTime) then
        totalTime = max_time
    end if

    ! Figure out how many digits per coupling there are 
    coupling_digits = get_num_digits(custom_string)

    ! If told to convert 
    if (convertTo > 0) then

        print *, "converting to ", convertTo, " digits"

        ! Convert 
        call convert_digits(custom_string, coupling_digits, convertTo)

        print "(A,A,A)", '"', trim(custom_string), '"'
        stop

    end if

    ! Set things based on the digits per coupling
    if (mutate_amount_initial < 0) then
        mutate_amount_initial = max(10**coupling_digits / 5, 1)
    end if
    mutate_amount = mutate_amount_initial
    max_val = 10**coupling_digits - 1

    ! Either solve the system using a GA or normally, depending
    if (custom) then
        call solve_genetic()
    else
        call solve(1)
    end if

contains

    function get_num_digits(genome_string)

        character(*), intent(in) :: genome_string
        integer :: get_num_digits

        integer :: i
        integer :: start_of_main

        start_of_main = 1

        ! Determine where the main section starts
        do i = 1, len_trim(genome_string)

            if (genome_string(i:i) == ">" .and. genome_string(i+1:i+1) /= "<" .and. genome_string(i+2:i+2) /= "@") then

                start_of_main = i+1
                exit

            end if

        end do

        ! Count the digits until reaching another char
        get_num_digits = 0
        do i = 1, len_trim(genome_string)

            ind = start_of_main+1+i
            asciiVal = ichar(genome_string(ind:ind))

            if (asciiVal >= 48 .and. asciiVal <= 57) then
                get_num_digits = get_num_digits + 1
            else
                exit
            end if

        end do

    end function

    subroutine convert_digits(genome_string, from, to)

        character(*), intent(inout) :: genome_string
        integer, intent(in) :: from, to
        integer :: i, delta, start_of_main

        delta = to-from
        start_of_main = 1

        if (delta == 0) return

        ! Determine where the main section starts
        do i = 1, len_trim(genome_string)

            if (genome_string(i:i) == ">" .and. genome_string(i+1:i+1) /= "<" .and. genome_string(i+2:i+2) /= "@") then

                start_of_main = i+1
                exit

            end if

        end do

        ! Go through the string
        i = start_of_main
        do while (i < len_trim(genome_string))

            ! Break if reaching vis info
            if (genome_string(i:i) == "#") exit

            if (is_int(genome_string(i:i)) .and. .not. is_int(genome_string(i+1:i+1))) then

                if (delta > 0) then
                    do j = 1, delta
                        genome_string(:) = trim(genome_string(:i)) // "0" // trim(genome_string(i+1:))
                        i = i + 1
                    end do
                else if (delta < 0) then
                    do j = 1, -delta
                        genome_string(:) = trim(genome_string(:i-1)) // trim(genome_string(i+1:))
                        i = i - 1
                    end do
                end if

            end if

            i = i + 1

        end do

    end subroutine

    ! See if a character string could be read as an integer
    function is_int(a)

        logical :: is_int
        character(*), intent(in) :: a
        integer :: i, ascii_val

        is_int = .true.
        do i = 1, len_trim(a)

            ! Make sure all of the ascii chars are digits
            ascii_val = ichar(a(i:i))
            if (ascii_val < 48 .or. ascii_val > 57) then
                is_int = .false.
                return
            end if

        end do

    end function

    ! See if a character string could be read as a real
    function is_real(a)

        logical :: is_real
        character(*), intent(in) :: a
        character(len(a)) :: a_copy
        integer :: i, ascii_val

        a_copy = ""

        ! Copy most chars over to the copy of a 
        do i = 1, len(a)

            ! Don't copy ".", "-" or "+"
            if (a(i:i) /= "." .and. a(i:i) /= "-" .and. a(i:i) /= "+") then
                a_copy = trim(a_copy) // a(i:i)
            end if

        end do

        is_real = .true.
        do i = 1, len_trim(a_copy)

            ! Make sure all of the ascii chars are digits
            ascii_val = ichar(a_copy(i:i))
            if (ascii_val < 48 .or. ascii_val > 57) then
                is_real = .false.
                return
            end if

        end do

    end function

end program
