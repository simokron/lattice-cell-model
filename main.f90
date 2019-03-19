!Author: simokron
!
!This is the lattice cell model implementation as part of my master's thesis.

!This module contains the parameters governing the simulation.
module constants
    implicit none

    !L and lambda are the lattice and cell dimensions, respectively; numIters is the number of Monte Carlo steps per frame.
    !beta is the reciprocal temperature; p0 is the concentration of zeroes at t = 0; p1 is the concentration of +1 (and -1 at the moment); phi is to volatility; cutoffConc is the final residual solvent concentration - set to negative number for infinite run-time.
    !To boolean constSeed uses a constant seed for the RNG (for debugging); FBC enables the free boundary conditions.
    !sigma is the spin matrix; numSpins is a tensor of rank 3 which stores the number of spins of each spices per cell.
    integer,parameter :: L = 256, lambda = 4, numIters = 2**22
!    real,parameter :: beta = 0.6, p0 = 0.6, p1 = (1 - p0)/2, phi = 0, cutoffConc = 0.1
    real,parameter :: beta = 0.6, p0 = 0.6, p1 = 0.3, phi = 0, cutoffConc = 0.1
    logical,parameter :: constSeed = .false., FBC = .true., debug = .false.
    integer :: sigma(L,L), numSpins(L/lambda,L/lambda,1:3)
!    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str))) !The result is a 3 x 3 row matrix, i.e. the first three values correspond to the elements in the first row, etc.
!    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !The result is a 3 x 3 row matrix, i.e. the first three values correspond to the elements in the first row, etc.

!    real,dimension(3, 3) :: J_str = transpose(reshape(1.0*[0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str))) !1
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.3*[0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str))) !2
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.08*[0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str))) !4
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.01*[0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str))) !8
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.0015*[0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str))) !16
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.0003*[0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str))) !32

!    real,dimension(3, 3) :: J_str = transpose(reshape(1.3*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !1
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.2*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !2
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.05*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !4
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.01*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !8
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.003*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !16
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.001*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !32

!    real,dimension(3, 3) :: J_str = transpose(reshape((14.1045*exp(-2.5148*lambda) + &
!        0.2353*exp(-0.3906*lambda))*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !scaled

!    real,dimension(3, 3) :: J_str = transpose(reshape((1.2936*lambda**(-2.7085)+0.0061-&
!        0.003)*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !scaled

!    real,dimension(3, 3) :: J_str = transpose(reshape(real(lambda)**(-2)*[0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str))) !J_ORIGINAL SCALED
!    real,dimension(3, 3) :: J_str = transpose(reshape(real(lambda)**(-2)*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !2
!    real,dimension(3, 3) :: J_str = transpose(reshape(real(lambda)**(-2)*0.01* &
!        [0, 80, 160, 80, 0, 80, 160, 80, 0], shape(J_str))) !Based on manual tests. Here be dragons.
    real,dimension(3, 3) :: J_str = transpose(reshape(real(lambda)**(-2)*0.01* &
        [0, 75, 125, 75, 0, 75, 125, 75, 0], shape(J_str))) !These are the best values, IMO.
 
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.01* &
!        [0, 75, 125, 75, 0, 75, 125, 75, 0], shape(J_str))) !UNSCALED

!    real,dimension(3, 3) :: J_str = transpose(reshape(real(lambda)**(-2.7)*0.1*[0, 6, 60, 6, 0, 6, 60, 6, 0], shape(J_str))) !TEST
!    real,dimension(3, 3) :: J_str = transpose(reshape(real(lambda)**(-2.65)*0.1*[0, 6, 60, 6, 0, 6, 60, 6, 0], shape(J_str))) !TEST
!    real,dimension(3, 3) :: J_str = transpose(reshape(real(lambda)**(-2.65)*0.1*[0, 10, 20, 10, 0, 10, 20, 10, 0], shape(J_str))) !TEST

!    real,dimension(3, 3) :: J_str = transpose(reshape(0.005*[0, 10, 12, 10, 0, 10, 12, 10, 0], shape(J_str))) !SCALE C
!    real,dimension(3, 3) :: J_str = transpose(reshape(0.1*[0, 10, 35, 10, 0, 10, 35, 10, 0], shape(J_str))) !1

!    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 1, 6, 1, 0, 3, 6, 3, 0], shape(J_str))) !This will form a 'cap' of +1, cf. fig. 4 in Andrea's paper.
!    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 35, 15, 35, 0, 35, 15, 35, 0], shape(J_str))) !This is the 'strong repulsion' in Andrea's paper. It kinda works but I need much more energy..
!    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 1, 6, 1, -2, 1, 6, 1, -2], shape(J_str))) !This is the 'strong repulsion' in Andrea's paper. It kinda works but I need much more energy..
end module

!This module contains the subroutines responsible for initialising the simulation.
module initialisation
    use constants

contains
    !This subroutine sets up the RNG for the session. If constSeed .eqv. .true., it will use a predefined seed (for debugging).
    subroutine setupRNG()
        integer :: seedConst(33) = (/(i, i=1,33, 1)/)
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        if(constSeed .eqv. .true.) then
            call random_seed(put=seedConst)
        else
            call random_seed(size = n)
            allocate(seed(n))

            call system_clock(count = clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            call random_seed(put = seed)

            deallocate(seed)
        endif
    end subroutine setupRNG

    !This subroutine generates the pseudo-random ternary spin matrix with appropriate concentrations.
    subroutine genSpins(L, sigma, p0, p1)
        integer :: L, sigma(L,L)
        real :: p0, p1

        integer :: i, j
        real :: u

        do i = 1,L
            do j = 1,L
                call random_number(u)

                if(u < p0) then
                    sigma(i,j) = 0
                elseif(u < p0 + p1) then
                    sigma(i,j) = 1
                else
                    sigma(i,j) = -1
                endif
            enddo
        enddo

        return

    end subroutine genSpins
end module initialisation

!This module contains the functions utilised by the dynamics module.
module functions
    use constants

contains
    !This function takes a bond number 'b' and returns the "ij"-coordinates of the originating cell (the bonds point up and to the right) as well as the type of bond (t = 0 (1) for horisontal (vertical) type) as an array.
    function cellCord(b)result(coordResult)
        integer :: t, coordResult(1:3)

        !Here we round down to the current integer (e.g. 0.8 would be 0) to find i_c. This first conditional is the edge-cases where e.g. on a 3 x 3 matrix, the last bond at the first row would have value 6, but then we'd get row 2 from the second conditional since 6/3 = 2.
        if(mod(int(b),(2*(L/lambda))) == 0) then
            i_c = int(aint(b/(2*(L/lambda))))
        else
            i_c = int(aint(b/(2*(L/lambda))) + 1)
        endif

        !Now that we have i_c, we can easily find j_c by determining if it is a horisontal or vertical bond and adjusting the formula correspondingly. This formula was found by considering what consecutive cell-number is associated with a specific set of i_c and j_c, then multiplying by two because you have two bonds per cell.
        !In other words; e.g. given i_c = 2 and j_c = 3, use the fact that L/lambda = 6 (say) to determine the number associated with the cell. The implementation is of course really the reverse; you are given the bond number and want to find which cell it is associated with (and later its pair in the bond!).
        !The two conditionals reflect the two bonds per cell (the formula changes slightly depending on if it's a vertical or horisontal bond).
        if(mod(int(b),2) == 0) then
            t = 0 !Horisontal type.
            j_c = int(b)/2 - (L/lambda)*(i_c-1)
        else
            t = 1 !Vertical type - NOTE THAT THIS IS POINTING DOWNWARDS!! (i.e. LARGER i <=> FURTHER DOWN) - throughout "up" will refer to the geometrical "up", i.e. a lower value of i.
            j_c = (int(b)+1)/2 - (L/lambda)*(i_c-1)
        endif

        !Note that the spin matrix is formatted according to row-column so spin(1,2) would be the spin at row 1, column 2. I.e. i_c represents the 'y' coordinate and j_c represents the 'x' coordinate.
        coordResult = [i_c, j_c, t]

        return
    end

    !This function takes a spin number together with the coordinates of the cell in question and finds the indices of the associated spin within the appropriate cell as an array.
    !The spin indices work with sigma(i_s, j_s).
    function spinCord(s, i_c, j_c)result(coordResult)
        integer :: coordResult(1:2)

        !Now we re-use the earlier ideas to get the 'coordinates' of the spin. Note that I no longer multiply by two, as there is only one "spin per spin" (i.e. we're not considering bonds). In the future I could make a "findI()" function, in which case it would probably make more sense to actually consider bonds in both cases to make the code more general. It shouldn't affect the outcome as the RNG works the same.
        if(mod(int(s),lambda) == 0) then
            i_s = int(aint(s/(lambda)))
        else
            i_s = int(aint(s/(lambda)) + 1)
        endif

        !Now we have i_s and can use a similar expression to find the j_s as j_c earlier. Note that since we don't have bonds, there is only one formulation.
        j_s = int(s) - lambda*(i_s - 1)

        !Now we "transform" the i_s, j_s coordinates to reflect which cell they inhibit. Thus far the coordinates have all been relative the cell, not the whole lattice!
        i_s = i_s + lambda*(i_c - 1)
        j_s = j_s + lambda*(j_c - 1)

        coordResult = [i_s, j_s]

        return
    end

    !This function takes a set of coordinates together with the maximum dimension (e.g. the length of a box) and returns the coordinates of the nearest-neighbours with PBC as an array.
    function PBC(i, j, length)result(PBCResult)
        integer, dimension(1:4) :: PBCResult
        integer :: x_left, x_right, y_up, y_down

        !Next we define some stuff to simplify the PBC (for the cells!). Note again that the 'up' cell is at a row with a lesser value of i_c!
        y_up = i - 1; y_down = i + 1; x_left = j - 1; x_right = j + 1

        !The first one: "if you are in the first column, the neighbour to your left will be the rightmost cell on the same row", and analogously for the rest
        if(j == 1) then
            x_left = length
        elseif(j == length) then
            x_right = 1
        endif
        if(i == 1) then
            y_up = length
        elseif(i == length) then
            y_down = 1
        endif

        PBCResult = [y_up, y_down, x_left, x_right] !This is a 1 x 4 vector. Originally I transposed them and reshaped it, but maaaaan does that take a lot of operations.

        return
    end

    !This function takes the coordinates of a cell and returns the indices of the first spin within that cell (i.e. the one in the top left-hand corner of the cell) as an array.
    function findFirst(i_c, j_c)result(firstResult)
        integer, dimension(1:2) :: firstResult

        !Finding the first coordinates is trivial. For i, simply take the number of rows to a cell and multiply that by the cell index minus one (the minus one is due to the cell index starting at 1 rather than 0) and add 1 because otherwise you get the last row in the cell above, and analogous for j.
        i_first = 1 + lambda*(i_c - 1)
        j_first = 1 + lambda*(j_c - 1)

        firstResult = [i_first, j_first]

        return
    end

    !This function takes the value of L and lambda as well as the initial spin configuration and returns the number of spins of a given species in a rank 3 tensor, where the first the indices relate to the row-column of the cell and the third to the spin species (-1,0,1)+2.
    function initialNum(L, lambda, sigma)result(numResult)
        integer :: L, lambda, sigma(L,L)
        integer, dimension(1:3) :: tempSpins
        integer, dimension(L/lambda,L/lambda,1:3) :: numResult

        integer :: i, j, first(1:2), x1, x2

        !The outer loop runs over the total number of cells.
        do i = 1, L/lambda
            do j = 1, L/lambda
                first = findFirst(i, j)

                !Reset for each run.
                tempSpins = [0,0,0]

                !Here we simply loop over each spin in the cell and counts the number of spins of species \alpha \in {-1, 0, 1}.
                do x1 = first(1), first(1) + lambda - 1
                    do x2 = first(2), first(2) + lambda - 1
                        spin_loop = sigma(x1, x2)
                        if(spin_loop == -1) then
                            tempSpins(1) = tempSpins(1) + 1
                        elseif(spin_loop == 0) then
                            tempSpins(2) = tempSpins(2) + 1
                        else
                            tempSpins(3) = tempSpins(3) + 1
                        endif
                    enddo
                enddo

                !Transform the result to an array.
                do k = 1, 3
                    numResult(i,j,k) = tempSpins(k)
                enddo

            enddo
        enddo

        return
    end

    !This function takes the value of the currently considered spin together with the value of the proposed spin switch, as well as the coordinates of the first spin in the relevant cell and the indices of the spin in question and returns the current energy, along with the proposed energy should the spins be switched in an array.
    function energySelected(spin, spin_p, t, numSpins, i_c, j_c, i_p_c, j_p_c, cellNeighbours)result(energyResult)
        integer :: spin, spin_p, t, numSpins(L/lambda,L/lambda,1:3), i_c, j_c, i_p_c, j_p_c, cellNeighbours(1:12) !Needed to specify the type - Fortran gets real mad if you try and use non-integers for indices (and it seems to always assume real for non-specified types).
        integer :: k, b, h
        integer, dimension(1:3) :: numSpinsA, numSpinsB, numSpinsAprop, numSpinsBprop
        integer, dimension(1:3) :: numSpinsIntA, numSpinsIntB
        real :: E_current, E_proposed, energyResult(1:2) 
        real :: c = 1.0*10**(0) !Allows for tuning of the interaction between cells.
        real :: c2 = 1.0

        !Reset the values before each run.
        E_current = 0
        E_proposed = 0

        !The number of spins in cell A and B can be trivially determined from,together with the number of spins of species k the currently considered spin interacts with (just the number of spins in that cell of the species, minus one to avoid self-interaction).
        do k = 1,3
            numSpinsA(k) = numSpins(i_c,j_c,k)
            numSpinsB(k) = numSpins(i_p_c,j_p_c,k)
            numSpinsAprop(k) = numSpins(i_c,j_c,k)
            numSpinsBprop(k) = numSpins(i_p_c,j_p_c,k)
            numSpinsIntA(k) = numSpins(i_c,j_c,k)
            numSpinsIntB(k) = numSpins(i_p_c,j_p_c,k)
        enddo
        do k = 1,3
            if(k == spin + 2) numSpinsIntA(k) = numSpinsIntA(k) - 1
            if(k == spin_p +2) numSpinsIntB(k) = numSpinsIntB(k) - 1
        enddo

        !In the proposed switch, the number of spins in each cell is altered thus
        do k = 1,3
            if(k == spin + 2) then
                numSpinsAprop(k) = numSpinsAprop(k) - 1
                numSpinsBprop(k) = numSpinsBprop(k) + 1
            endif
            if(k == spin_p + 2) then
                numSpinsAprop(k) = numSpinsAprop(k) + 1
                numSpinsBprop(k) = numSpinsBprop(k) - 1
            endif
        enddo

        !The energy is determined using the number of spins per cell, since their location doesn't matter in mean-field.
        do k = 1, 3
            E_current = E_current + c2*numSpinsIntA(k)*J_str(spin+2,k) + c2*numSpinsIntB(k)*J_str(spin_p+2,k)
            E_proposed = E_proposed + c2*numSpinsIntA(k)*J_str(spin_p+2,k) + c2*numSpinsIntB(k)*J_str(spin+2,k)
        enddo

        !And the nearest-neighbouring cell-cell interaction (note that the FBC changes the calculations with the neighbours slightly).
        !This section is hard to follow without the aid of the figure showing the nearest-neighbour cluster - see the thesis.
        do k = 1, 3
            do b = 1, 3
                E_current = E_current + c*numSpinsA(k)*numSpinsB(b)*J_str(k,b)
                E_proposed = E_proposed + c*numSpinsAprop(k)*numSpinsBprop(b)*J_str(k,b)
            enddo
        enddo
        do h = 1, 5, 2
            if(FBC .eqv. .true.) then
                if(i_c == L/lambda .and. t == 0) then
                    if(h == 5) then
                        go to 20
                    endif
                elseif(i_c == L/lambda .and. t == 1) then
                    if(h == 3) then
                        go to 20
                    endif
                elseif(i_c == 1 .and. t == 0) then
                    if(h == 1) then
                        go to 20
                    endif
                endif
            endif
            do k = 1, 3
                do b = 1, 3
                    E_current = E_current + c*numSpinsA(k)*numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                    E_proposed = E_proposed + c*numSpinsAprop(k)*numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                enddo
            enddo
20          continue
        enddo
        do h = 7, 11, 2
            if(FBC .eqv. .true.) then
                if(i_c == L/lambda .and. t == 0) then
                    if(h == 7) then
                        go to 30
                    endif
                elseif(i_c == 1 .and. t == 0) then
                    if(h == 11) then
                        go to 30
                    endif
                elseif(i_c == 1+1 .and. t == 1) then
                    if(h == 9) then
                        go to 30
                    endif
                endif
            endif
            do k = 1, 3
                do b = 1, 3
                    E_current = E_current + c*numSpinsB(k)*numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                    E_proposed = E_proposed + c*numSpinsBprop(k)*numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                enddo
            enddo
30          continue
        enddo

        energyResult = [E_current, E_proposed]

        return
    end

end module functions

!This module contains the dynamics in the simulation.
module dynamics
    use constants
    use functions

contains
    !This subroutine handles the evaporation of the top row by replacing the spin at (i,j) by +1 or -1 according to the appropriate concentration and returning the spin matrix.
    subroutine evap(i_s, j_s, i_c, j_c, sigma, numSpins)
        integer :: i_s, j_s, i_c, j_c, sigma(L,L), numSpins(L/lambda, L/lambda, 1:3)
        !Simply get a random number between 0 and 1 and compare it to the normalised concentrations, preserving the initial ratios of -1 and +1 species.
        call random_number(P)
        numSpins(i_c,j_c,0+2) = numSpins(i_c,j_c,0+2) - 1
        if(P < p1/(1 - p0)) then
            sigma(i_s,j_s) = 1
            numSpins(i_c,j_c,1+2) = numSpins(i_c,j_c,1+2) + 1
        else
            sigma(i_s,j_s) = -1
            numSpins(i_c,j_c,-1+2) = numSpins(i_c,j_c,-1+2) + 1
        endif
    end subroutine evap
    
    !Pseudo-random evaporation test (probably shite).
    subroutine evapRand(i_s, j_s, i_c, j_c, sigma, numSpins)
        integer :: i_s, j_s, i_c, j_c, sigma(L,L), numSpins(L/lambda, L/lambda, 1:3), spinT
        integer :: i_c_evap, j_c_evap, i_s_evap, j_s_evap, spinTemp_evap(1:2), spin_evap, k
        real :: u, conc

        k = 0
        conc = real(count(sigma == 0))/real(L**2) !Conc after evap event.
        
60      call random_number(u); i_c_evap = 1 + floor((L/lambda)*u) !This yields an integer between 1 and L/lambda
        call random_number(u); j_c_evap = 1 + floor((L/lambda)*u)

        !Now we randomly select a spin inside of the first cell of the bond.
        call random_number(u); s_evap = 1 + floor((lambda**2)*u) !This yields an integer between 1 and lambda^2 (i.e. the total number of spins in a cell).

        !And determine the "xy-coordinates" of the spin.
        spinTemp_evap = spinCord(s_evap, i_c_evap, j_c_evap)
        i_s_evap = spinTemp_evap(1); j_s_evap = spinTemp_evap(2)

!        print *, 'i_s_evap = ', i_s_evap
!        print *, 'j_s_evap = ', j_s_evap

        !Finally we have the actual spin at i_s and j_s, which we store in spin.
        spin_evap = sigma(i_s_evap, j_s_evap)
        k = k + 1

        if(spin_evap /= 0 .and. k < (1/conc) .and. k < 1000) goto 60

        spinT = sigma(i_s, j_s)

        call updateSigma(spinT, spin_evap, sigma, i_s, j_s, i_s_evap, j_s_evap)
        call recalcSpins(spinT, spin_evap, numSpins, i_c, j_c, i_c_evap, j_c_evap)
    end subroutine evapRand

    !This subroutine updates sigma, should the proposed move be approved (trivial, but convenient to have as a subroutine).
    subroutine updateSigma(spin, spin_p, sigma, i_s, j_s, i_p_s, j_p_s)
        integer :: spin, spin_p, sigma(L,L), i_s, j_s, i_p_s, j_p_s

        !Simply switch the spins at i_s, j_s and i_p_s, j_p_s via the spin and spin_p variables.
        sigma(i_s, j_s) = spin_p
        sigma(i_p_s, j_p_s) = spin

        return
    end subroutine updateSigma

    !This subroutine 'recalculates' the numSpins tensor (really it just updates it - there is no need for a bunch of calculations).
    subroutine recalcSpins(spin, spin_p, numSpins, i_c, j_c, i_p_c, j_p_c)
        integer :: spin, spin_p, numSpins(L/lambda,L/lambda,1:3), i_c, j_c, i_p_c, j_p_c
        integer :: k

        !Simply loop over the number of spin species and adjust accordingly.
        do k = 1,3
            if(k == spin + 2) then
                numSpins(i_c,j_c,k) = numSpins(i_c,j_c,k) - 1
                numSpins(i_p_c,j_p_c,k) = numSpins(i_p_c,j_p_c,k) + 1
            endif
            if(k == spin_p + 2) then
                numSpins(i_c,j_c,k) = numSpins(i_c,j_c,k) + 1
                numSpins(i_p_c,j_p_c,k) = numSpins(i_p_c,j_p_c,k) - 1
            endif
        enddo

        return
    end subroutine recalcSpins

    !This subroutine is really the main part of the program. It essentially takes a spin matrix (state) and returns an updated version of the spin matrix after numIters iterations.
    !Note that I currently do not keep track of any observables like energy etc for the whole system, but it would be trivial to implement
    subroutine metropolis(L, lambda, sigma, numSpins, numIters, beta)
        !These were defined in the constants module.
        integer :: L, lambda, sigma(L,L), numSpins(L/lambda,L/lambda,1:3), numIters
        real :: beta

        !This section is a bit of a mess still. The general rules are;
        !index "c" := related to cells, so i_c, j_c will hold the "xy-coordinates" of the currently considered cell.
        !index "s" := related to spins, so i_s, j_s will hold the "xy-coordinates" of the currently considered spin (note however that these are really two different sets of coordinates as there exist fewer cells than spins except in the edge-case of lambda = 1). For the spins, these are literally the coordinates of the spins, and thus sigma(i_s,j_s) will return the value of the spin currently being considered. This of course makes no sense for spin(i_c,j_c).
        !index "p" indicates a ''paired'' quantity, i.e. i_p_c, j_p_c will be the "xy-coordinates" of the bond in the *paired* cell in the bond.
        !The odd ones out then are y_down, y_up, etc. which relates to the nearest-neighbouring cells.
        !t stands for ''type'' and it is used to keep track of if the bond is a vertical or horisontal one (used for the volatility which is currently disabled).
        !b is the currently considered bond. It is really an integer but because of the way I generate it, it needs to be a real for now.
        !u is used in the RNG.
        !w and P are used to check for acceptance criteria (Boltzmaan factor).
        integer :: i_c, j_c, k, spin, y_down, y_up, x_left, x_right, t
        integer :: i_s, j_s, i_p_s, j_p_s, y_p_down, y_p_up, x_p_left, x_p_right
        integer :: i_p_c, j_p_c, spin_p
        integer, dimension(1:4) :: PBCTemp !Temporary array for storing the coordinates from the PBC function.
        integer, dimension(1:2) :: spinTemp !Temporary array for storing the coordinates of the spin.
        integer, dimension(1:3) :: cell1 !Temporary array for storing the coordinates and type of bond of the cells.
        integer, dimension(1:12) :: cellNeighbours !Array to store the coordinates of the 6 neighbours (see image in thesis).
        real, dimension(1:2) :: energy !Array for storing the energy before and after proposed spin switch.
        real :: b, s, u, dE, w, P

        !This is the actual Metropolis algorithm. Again, somewhat messy.
        do k = 1,numIters
            !-Cell and spin selection-------------------------------------------
            !First we randomly select a bond.
    10      call random_number(u); b = 1 + floor(2*((L/lambda)**2)*u) !This yields an integer between 1 and 2*(L/lambda)^2 (i.e. the total number of bonds between cells). Whence b will be the currently considered bond number.

            !Then we determine the "ij-coordinates" of the origin cell of the bond.
            cell1 = cellCord(b)
            i_c = cell1(1); j_c = cell1(2); t = cell1(3)

            !Now we randomly select a spin inside of the first cell of the bond.
            call random_number(u); s = 1 + floor((lambda**2)*u) !This yields an integer between 1 and lambda^2 (i.e. the total number of spins in a cell).

            !And determine the "x1x2-coordinates" of the spin.
            spinTemp = spinCord(s, i_c, j_c)
            i_s = spinTemp(1); j_s = spinTemp(2)

            !Finally we have the actual spin at i_s and j_s, which we store in spin.
            spin = sigma(i_s, j_s)

            !Now we call the PBC() function which returns a 2 x 2 matrix with the "xy-cordinates" of the nearest-neighbours.
            PBCtemp = PBC(i_c, j_c, L/lambda)
            y_up = PBCTemp(1); y_down = PBCTemp(2); x_left = PBCTemp(3); x_right = PBCTemp(4)

            !Next we find the position of the other cell associated with the bond, which we have denoted _p for 'pair'.
            !This is simple: "if we have a horisontal bond, the pair must be the cell to the right; otherwise it is the cell directly above".
            !That is, each cell has a vertical bond pointing upwards and a horisontal bond pointing to the right (it could equally well have been the left; physics is invariant to arbitrary coordinate definitions).
            !Then we just call PBC like before and store the nearest neighbouring cells. Note that this is technically overkill, as one of the nearest neighbours will be the originating cell, but it's so quick that it doesn't matter - there are WAAAAAAAY worse parts of this program.
            if(t == 0) then
                i_p_c = i_c; j_p_c = x_right
            else
                i_p_c = y_up; j_p_c = j_c
            endif
            PBCtemp = PBC(i_p_c, j_p_c, L/lambda)
            y_p_up = PBCTemp(1); y_p_down = PBCTemp(2); x_p_left = PBCTemp(3); x_p_right = PBCTemp(4)

            !Now we randomly select a spin inside of the paired cell of the bond.
            call random_number(u); s = 1 + floor((lambda**2)*u) !This *should* be an integer between 1 and lambda^2 (i.e. the total number of spins in a cell).

            !And determine the "x1x2-coordinates" of the spin.
            spinTemp = spinCord(s, i_p_c, j_p_c)
            i_p_s = spinTemp(1); j_p_s = spinTemp(2)

            spin_p = sigma(i_p_s, j_p_s) !Now we have the actual spin at i_p and j_p, which we store in spin_p.

            !-Evaporation of top row--------------------------------------------
            if(i_c == 1 .and. spin == 0) then
                call evap(i_s, j_s, i_c, j_c, sigma, numSpins)
                GO TO 10
            endif

            !-Solvent volatility------------------------------------------------
            !This is the upwards drift due to the volatility.
            if(phi == 0) GO TO 30
            if(t == 1 .and. spin_p == 0 .and. spin /= 0) then
                call random_number(P) !Compare to a pseudo-random number between 0 and 1.
                if(P > phi) then
                    call updateSigma(spin, spin_p, sigma, i_s, j_s, i_p_s, j_p_s)
                    call recalcSpins(spin, spin_p, numSpins, i_c, j_c, i_p_c, j_p_c)
                endif
                GO TO 10
            endif
30          continue

            if(i_c == 1 .and. t == 1) GO TO 10 !Avoid moves across top/bottom boundary.
            
            !-Dynamics----------------------------------------------------------
            !Before we compute the energy, we must find the nearest-neighbouring cells.
            if(t == 0) cellNeighbours = [[y_up,j_c], [i_c,x_left], [y_down,j_c], &
                [y_p_down,j_p_c], [i_p_c,x_p_right], [y_p_up,j_p_c]]
            if(t == 1) cellNeighbours = [[i_c,x_left], [y_down,j_c], [i_c,x_right], &
                [i_p_c,x_p_right], [y_p_up,j_p_c], [i_p_c,x_p_left]]

            !Now we simply call the energySelected function, with the needed information and BAM! Just like that, we're done. Thanks, objective oriented programming.
            energy = energySelected(spin, spin_p, t, numSpins, i_c, j_c, i_p_c, j_p_c, cellNeighbours)
            dE = energy(2) - energy(1)

            !Now we simply insert the acceptance criteria from the model.
            w = exp(-beta*dE)

            !Compare to a pseudo-random number between 0 and 1.
            call random_number(P)

            !And switch spins if it is energetically favourable or uphill if P <= w.
            if(dE <= 0 .or. P <= w) then
                call updateSigma(spin, spin_p, sigma, i_s, j_s, i_p_s, j_p_s)
                call recalcSpins(spin, spin_p, numSpins, i_c, j_c, i_p_c, j_p_c)
            endif

        enddo
        
        return

    end subroutine metropolis

end module dynamics

!This module takes care of number to string conversion for creation of directories etc.
module strings
    use constants

    ! GLOBAL FUNCTIONS
    public :: num2str

    ! Everything else is private
    private

    interface num2str
      module procedure num2str_int
      module procedure num2str_real
    end interface 

contains

    function num2str_int(number)
        implicit none
        integer,intent(in) :: number
        character(len=6)   :: num2str_int
        character(len=6)   :: tmp

        write(tmp,'(I6)')number
        num2str_int = tmp
    end function

    function num2str_real(number)
        implicit none
        real,intent(in)    :: number
        character(len=6)   :: num2str_real
        character(len=6)   :: tmp
        
        write(tmp,'(F6.4)')number
        num2str_real = tmp
    end function
end module

!The main program controls the simulation in the sense that it initialises the system, writes .dat files, informs the user about approximate time remaining etc.
program main
    use strings
    use initialisation
    use functions
    use dynamics
    implicit none

    !-Input parameters----------------------------------------------------------
    !n, i and j are used in loops; stat is used to remove old frames.
    !file_id and file_name hold info for saving/erasing the frames.
    !start, finish, etc. are only used to keep track of the time and to provide some nice prompts to the user; conc is the current concentration of zeros; MCS is the current number of MCS (for prompt).
    integer :: i, j, n, stat, nI
    character(128) :: file_id, file_name, folderName
    character :: d
    logical :: file_exists = .true., folder_exists = .true.
    real :: start, finish, numHour, numMin, numSec, conc, MCS

    !The debugging mode simply uses /frames and deletes all of the files in there without warning!
    !In non-debugging mode, this creates a directory with the correct name (if it does not currently exist).
    if(debug .eqv. .true.) then
        folderName = 'frames'
        GO TO 70
    else
        nI = aint(log(real(numIters))/log(2.))
        folderName = 'lambda_' // trim(adjustl(num2str(lambda))) // '-L_' // trim(adjustl(num2str(L))) // '-J_' // &
            trim(adjustl(num2str(J_str(1,1)*lambda**2))) // '_' // trim(adjustl(num2str(J_str(1,2)*lambda**2))) // '_' // &
            trim(adjustl(num2str(J_str(1,3)*lambda**2))) // '-numIters_2-' // trim(adjustl(num2str(nI)))
        if(FBC .eqv. .true.) folderName = trim(folderName) // '-FBC'
        if(FBC .eqv. .false.) folderName = trim(folderName) // '-PBC'
        inquire(file = './' // trim(folderName), exist = folder_exists)
        if(folder_exists .eqv. .true.) GO TO 70
        call execute_command_line ('mkdir ' // folderName)
    endif

    !This clears the directory using magic.
70  n = 0
    do while (file_exists .eqv. .true.)
        n = n + 1
        write(file_id, '(i0)') n
        file_name = trim(adjustl(folderName)) // '/frame-' // trim(adjustl(file_id)) // '.dat'
        inquire(file = trim(file_name), exist = file_exists)
        if(debug .eqv. .false.) then
            if(n == 1 .and. file_exists .eqv. .true.) then
                print '(//)'; print *, 'The directroy ' // trim(folderName) // ' is not empty!'
                print '(">Delete old files? (y/n)")'; read *, d
                if(d == 'y') then
                    print '(//,"Deleting files...")'
                    GO TO 80
                else
                    print '(//,"Aborting!")'
                    GO TO 90
                endif
            endif
        endif
80      open(10, iostat=stat, file = trim(file_name), form = 'formatted', status = 'old')
        if(stat == 0) close(10, status='delete')
    enddo
    if(d == 'y') print '("Files deleted.")'

!    !Debugging stuff - basically gives one the ability to continue an aborted simulation from a specified state (REMEMBER TO COMMENT OUT THE STUFF ABOVE - IT *WILL* REMOVE ALL OF THE FILES OTHERWISE! AND ALSO THE GENSPINS CALL BELOW).
!    write(file_id, '(i0)') 470
!    file_name = 'frames/frame-' // trim(adjustl(file_id)) // '.dat'
!    open(10, file = trim(file_name), form = 'formatted')
!    do i = 1,L
!        read(10,*) (sigma(i,j), j = 1,L)
!    enddo
!    close(10)

    call setupRNG() !Sets up the RNG (either with a constant seed or using UNIX time for m-m-m-aximum pseudo-randomness!).
    call genSpins(L, sigma, p0, p1) !Generates a pseudo-random L x L "tenary spin matrix".

    numSpins = initialNum(L, lambda, sigma) !Determines the initial number of spins.
    conc = real(count(sigma == 0))/real(L**2) !Calculates the intial concentration.

    print '(//,"Starting dynamics...")'
    call cpu_time(start) !Keeps track of the time.
    !This is the main loop; it writes the .dat files and calls the metropolis() subroutine to alter the spin matrix and it also informs the user about the number of frames, time remaining etc.
    n = 1
    do while (conc > cutoffConc)
        write(file_id, '(i0)') n
        file_name = trim(adjustl(folderName)) // '/frame-' // trim(adjustl(file_id)) // '.dat'
        open(10, file = trim(file_name), form = 'formatted')
        do i = 1,L
            write(10,*) (sigma(i,j), j = 1,L)
        enddo
        close(10)

        call metropolis(L, lambda, sigma, numSpins, numIters, beta) !Returns the spin matrix after numIters.
        conc = real(count(sigma == 0))/real(L**2) !Recalculate the concentration after numIters and compare to cutoffConc.

        !Just some nice feedback to the user.
        if(mod(n, 10) == 0 .or. n == 2) then
!            print '(/,"Frame ", 18x, i6)',n
            MCS = numIters*n/L**2
            print '(/,"MCS: ", 17x, F10.0)', MCS
            print '("Concentration of zeros: ", F10.2)', conc
            call cpu_time(finish)
            numHour = aint(finish/3600)
            numMin = aint((finish - (numHour)*3600)/60)
            numSec = aint((finish - (numHour)*3600 - numMin*60))
            print '("Total elapsed time: ", 9x,i2," hour(s), ",i2," minute(s) and ",i2," second(s).")'&
            ,int(numHour),int(floor(numMin)),int(floor(numSec))
        endif
        n = n + 1
    enddo

    !Final feedback after completion.
    call cpu_time(finish)
    numHour = aint(finish/3600)
    numMin = aint((finish - (numHour)*3600)/60)
    numSec = aint((finish - (numHour)*3600 - numMin*60))
    print '(//,"Completed after: ", 12x,i2," hour(s), ",i2," minute(s) and ",i2," second(s).")'&
    ,int(numHour),int(floor(numMin)),int(floor(numSec))

90 continue
end program main
