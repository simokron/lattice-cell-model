!Author: simokron
!
!This is the lattice cell model implementation as part of my master's thesis.

!This module contains the parameters governing the simulation.
module constants
    implicit none

    integer,parameter :: L = 128, lambda = 1, numIters = 2**22
    real,parameter :: beta = 0.6, p0 = 0.4, p1 = (1 - p0)/2, phi = 0, cutoffConc = -0.1
!    real,parameter :: beta = 0.6, p0 = 0.4, p1 = 0.5, phi = 0, cutoffConc = 0.1
    integer :: sigma(L,L), n
    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str))) !The result is a 3 x 3 row matrix, i.e. the first three values correspond to the elements in the first row, etc.
!    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 1, 6, 1, 0, 3, 6, 3, 0], shape(J_str))) !This will form a 'cap' of +1, cf. fig. 4 in Andrea's paper.
!    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 35, 15, 35, 0, 35, 15, 35, 0], shape(J_str))) !This is the 'strong repulsion' in Andrea's paper. It kinda works but I need much more energy..
end module

!This module contains the subroutine responsible for initialising the simulation.
module initialisation
    use constants

contains
    !Thus subroutine generates the pseudo-random ternary spin matrix with appropriate concentrations.
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
    !This function takes a bond number 'b' and returns the "xy"-coordinates of the originating cell (the bonds point up and to the right) as well as the type of bond (t = 0 (1) for horisontal (vertical) type) as an array.
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

    !This function takes the value of the currently considered spin together with the value of the proposed spin switch, as well as the coordinates of the first spin in the relevant cell and the indices of the spin in question and returns the current energy, along with the proposed energy should the spins be switched in an array.
    function energySelected(spin, spin_p, i_first, j_first, i_s, j_s)result(energyResult)
        integer :: spin, spin_p, i_first, j_first, i_s, j_s !Needed to specify the type - Fortran gets real mad if you try and use non-integers for indices (and it seems to always assume real for non-specified types).
        real, dimension (1:2) :: energyResult
        integer :: i_loop, j_loop, spin_loop
        real :: E_current, E_proposed!, d_ik

        !Reset the value before each run.
        E_current = 0
        E_proposed = 0
        
        !Now we just loop for all of the spins in the cell (excluding the selected one) and calculate the energy from the selected spin interacting with the "loop spin".
        !Here I reset the energy for each numIters loop.
        do i_loop = i_first, i_first + lambda - 1
            do j_loop = j_first, j_first + lambda - 1
                if(i_loop == i_s .and. j_loop == j_s) GO TO 20
!                d_ik = sqrt(real(i_s - i_loop)**2 + real(j_s - j_loop)**2)
                spin_loop = sigma(i_loop, j_loop)
                E_current = E_current + J_str(spin + 2, spin_loop + 2)/2
                E_proposed = E_proposed + J_str(spin_p + 2, spin_loop + 2)/2
20              CONTINUE
            enddo
        enddo  
        
        energyResult = [E_current, E_proposed]
        
        return
    end

    !This function takes the index of the first spin in a cell and returns the 'spin' associated with that cell as an integer with value -1, 0 or +1.
    !The 'spin' associated with a cell is determined from the concentrations - e.g. a 2 x 2 cell with 75% species +1 and 25% species 0 will in 75% of cases be considered a 'spin +1 cell' in the cell-cell interaction.
    function findSpins(i_first, j_first)result(spinResult)
        integer :: spinResult
        integer :: i_loop, j_loop, spin_loop
        real, dimension(1:3) :: numSpins
        real, dimension(1:3) :: conc

        !Reset for each run.
        numSpins = [0,0,0]

        !Now we loop for all of the spins in the cell and count the number of spins with each value.
        do i_loop = i_first, i_first + lambda - 1
            do j_loop = j_first, j_first + lambda - 1
                spin_loop = sigma(i_loop, j_loop)
                k = k + 1
                if(spin_loop == -1) then
                    numSpins(1) = numSpins(1) + 1
                elseif(spin_loop == 0) then
                    numSpins(2) = numSpins(2) + 1
                else
                    numSpins(3) = numSpins(3) + 1
                endif
                    
            enddo
        enddo
        
        !And convert the values to concentrations by dividing by lambda^2 (total number of spins in the cell).
        conc = [numSpins(1), numSpins(2), numSpins(3)]/lambda**2

        !And compare the concentrations to a pseudo-random number between 0 and 1 to determine the 'spin of the cell'.
        call random_number(u)
        if(u < conc(1)) then
            spinResult = -1
        elseif(u < conc(1) + conc(2)) then
            spinResult = 0
        else
            spinResult = 1
        endif
        
        return
    end

    !This function takes the spin of the eight cells involved in the cell-cell nearest-neighbour interaction and returns the energy of the current configuration, as well as that after the proposed spin-switch as an array.
    !The cellSpins array has 10 entries due to the stochastic nature of selecting the 'cell spins' - the first (last) two entries are the current (proposed) 'cell spins' according to findSpins(), whilst the remaining nearest-neighbour cells are invariant to the proposed spin-switch.
    function energyCells(cellSpins)result(energyResult)
        integer, dimension (1:10) :: cellSpins 
        real, dimension (1:2) :: energyResult
        integer :: h
        real :: F_current, F_proposed
    
        !Reset before every run.
        F_current = 0
        F_proposed = 0
        
        !Now we loop for all of the spins in the nearest-neighbour cell-cell interaction and determine the current and proposed energies.
        F_current = F_current + J_str(cellSpins(1) + 2, cellSpins(2) + 2)/2 &
        + J_str(cellSpins(2) + 2, cellSpins(1) + 2)/2
        F_proposed = F_proposed + J_str(cellSpins(9) + 2, cellSpins(10) + 2)/2 &
        + J_str(cellSpins(10) + 2, cellSpins(9) + 2)/2
        do h = 3, 5
            F_current = F_current + J_str(cellSpins(1) + 2, cellSpins(h) + 2)/2
            F_proposed = F_proposed + J_str(cellSpins(9) + 2, cellSpins(h) + 2)/2
        enddo  
        do h = 6, 8
            F_current = F_current + J_str(cellSpins(2) + 2, cellSpins(h) + 2)/2
            F_proposed = F_proposed + J_str(cellSpins(10) + 2, cellSpins(h) + 2)/2
        enddo
        
        energyResult = [F_current, F_proposed]
        
        return
    end

end module functions

!This module contains the dynamics in the simulation.
module dynamics
    use constants
    use functions
    
contains
    !This subroutine handles the evaporation of the top row by replacing the spin at (i,j) by +1 or -1 according to the appropriate concentration and returning the spin matrix.
    subroutine evap(i, j)
        !Simply get a random number between 0 and 1 and compare it to the normalised concentrations, preserving the initial ratios of -1 and +1 species.
        call random_number(P)
        if(P < p1/(1 - p0)) then
            sigma(i,j) = 1
        else
            sigma(i,j) = -1
        endif
    end subroutine evap

    !This subroutine is really the main part of the program. It essentially takes a spin matrix (state) and returns an updated version of the spin matrix after numIters iterations.
    !Note that I currently do not keep track of any observables like energy etc for the whole system, but it would be trivial to implement
    subroutine metropolis(L, lambda, sigma, numIters, beta)
        !These were defined in the constants module.
        integer :: L, lambda, sigma(L,L), numIters
        real :: beta
    
        !This section is a bit of a mess still. The general rules are;
        !index "c" := related to cells, so i_c, j_c will hold the "xy-coordinates" of the currently considered cell.
        !index "s" := related to spins, so i_s, j_s will hold the "xy-coordinates" of the currently considered spin (note however that these are really two different sets of coordinates as there exist fewer cells than spins except in the edge-case of lambda = 1). For the spins, these are literally the coordinates of the spins, and thus sigma(i_s,j_s) will return the value of the spin currently being considered. This of course makes no sense for spin(i_c,j_c).
        !index "p" indicates a ''paired'' quantity, i.e. i_p_c, j_p_c will be the "xy-coordinates" of the bond in the *paired* cell in the bond.
        !The odd ones out then are y_down, y_up, etc. which relates to the nearest-neighbouring cells (and thus also the PBC). In the future I will make a function/subroutine for determining the nearest-neighbours of the cells, as well as the PBC.
        !t stands for ''type'' and it is used to keep track of if the bond is a vertical or horisontal one (used for the volatility which is currently disabled).
        !b is the currently considered bond. It is really an integer but because of the way I generate it, it needs to be a real for now.
        !u is used in the RNG.
        !E_current1 etc. relate to the current and proposed energies of the system (due to a spin move! I am yet to implement cell-cell interaction).
        !w and P are used to check for acceptance criteria (Boltzmaan factor).
        integer :: i_c, j_c, k, h, spin, y_down, y_up, x_left, x_right, t
        integer :: i_s, j_s, i_p_s, j_p_s, y_p_down, y_p_up, x_p_left, x_p_right
        integer :: i_p_c, j_p_c, spin_p
        integer, dimension(1:4) :: PBCTemp
        integer, dimension(1:2) :: first, first_p, spinTemp
        integer, dimension(1:3) :: cell1
        integer, dimension(1:12) :: cellNeighbours
        integer, dimension(1:16) :: cellFirsts
        integer, dimension(1:10) :: cellSpins
        real, dimension(1:2) :: Energy_c, Energy_p_c, cellEnergy
        real :: b, s, u, dE, w, P, P_c, P_comb, E_current, E_proposed
    
        !This is the actual Metropolis algorithm. Again, somewhat messy.
        do k = 1,numIters
            !First we randomly select a bond.
    10      call random_number(u); b = 1 + floor(2*((L/lambda)**2)*u) !This yields an integer between 1 and 2*(L/lambda)^2 (i.e. the total number of bonds between cells). Whence b will be the currently considered bond number.
    
            !Then we determine the "xy-coordinates" of the origin cell of the bond.
            cell1 = cellCord(b)
            i_c = cell1(1); j_c = cell1(2); t = cell1(3)

            !Now we randomly select a spin inside of the first cell of the bond.
            call random_number(u); s = 1 + floor((lambda**2)*u) !This yields an integer between 1 and lambda^2 (i.e. the total number of spins in a cell).

            !And determine the "xy-coordinates" of the spin.
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
            !With this knowledge, we can also trivially determine the nearest-neighbours of the paired cell without using PBC() again since it adds lots of overhead.
            if(t == 0) then
                i_p_c = i_c; j_p_c = x_right
            else
                i_p_c = y_up; j_p_c = j_c
            endif
            PBCtemp = PBC(i_p_c, j_p_c, L/lambda)
            y_p_up = PBCTemp(1); y_p_down = PBCTemp(2); x_p_left = PBCTemp(3); x_p_right = PBCTemp(4)
    
            !Now we randomly select a spin inside of the paired cell of the bond.
            call random_number(u); s = 1 + floor((lambda**2)*u) !This *should* be an integer between 1 and lambda^2 (i.e. the total number of spins in a cell).
            
            !And determine the "xy-coordinates" of the spin.
            spinTemp = spinCord(s, i_p_c, j_p_c)
            i_p_s = spinTemp(1); j_p_s = spinTemp(2)
    
            spin_p = sigma(i_p_s, j_p_s) !Now we have the actual spin at i_p and j_p, which we store in spin_p.

            !-Evaporation of top row--------------------------------------------
            if(i_c == 1 .and. spin == 0) then
                call evap(i_s,j_s)
                GO TO 10
            endif
    
            !-Solvent volatility------------------------------------------------
            !This is the upwards drift. It's pretty janky at the moment since it just switches spin and spin_p with some probability phi when spin = 0 and spin_p /= 0.
            if(t == 1 .and. spin == 0 .and. spin_p /= 0 .and. i_c /= 1) then
                if(phi == 0) GO TO 30
                call random_number(P) !Compare to a pseudo-random number between 0 and 1.
                if(P < phi) then
                    sigma(i_s, j_s) = spin_p; sigma(i_p_s, j_p_s) = spin;
                endif
                GO TO 10
            endif
30          continue

            if(i_c == 1 .and. t == 1) GO TO 10 !Avoid moves across top/bottom boundary.
    
            !-Dynamics----------------------------------------------------------
            !First we find the "first" spin in the current cell (i.e. the "xy-coordinates" of that spin).
            first = findFirst(i_c, j_c)
            
            !Then we determine the energy in that cell (both the current and the proposed energy).
            energy_c = energySelected(spin, spin_p, first(1), first(2), i_s, j_s)
    
            !Now we do the exact same thing for the spins in the paired cell (note that spin and spin_p have the opposite meaning from the paired spin's point of view).
            first_p = findFirst(i_p_c, j_p_c)
            energy_p_c = energySelected(spin_p, spin, first_p(1), first_p(2), i_p_s, j_p_s)
            
            !This gives us the current energy, as well as the proposed energy. Note that these energies are due to the spin move and do not correspond to cell-cell energies!
            E_current = energy_c(1) + energy_p_c(1)
            E_proposed = energy_c(2) + energy_p_c(2)
    
            !Finally we compute the energy difference should the spin be switched.
            dE = E_proposed - E_current

            !For the second conditional, we must consider the inter-cell energy. Luckily we can re-use some of the logic from before.
            !To start with, we must find the "xy-coordinates" of the nearest-neighbouring cells. 
            if(t == 0) cellNeighbours = [[y_up,j_c], [y_down,j_c], [i_c,x_left], &
                [y_p_up,j_p_c], [y_p_down,j_p_c], [i_p_c,x_p_right]]
            if(t == 1) cellNeighbours = [[y_down,j_c], [i_c,x_left], [i_c,x_right], &
                [y_p_up,j_p_c], [i_p_c,x_p_left], [i_p_c,x_p_right]]

            !Now we need to find the ''dominating spin'' in each cell
            cellFirsts([1,2]) = [first]; cellFirsts([3,4]) = [first_p]
            do h = 5, 15, 2
                cellFirsts([h,h+1]) = findFirst(cellNeighbours(h-4), cellneighbours(h-4+1))
            enddo

!            print *, "t = ", t
!            print *, "x_p_right = ", x_p_right
!            print *, "test1 = ", [first, first_p, cellNeighbours]
!            print *, "test2 = ", cellFirsts

            do h = 1, 8
                cellSpins(h) = findSpins(cellFirsts(2*h-1),cellFirsts(2*h))
            enddo
            
            !And the ''dominating spin'' in the relevant cells should they switch (i.e.\ cell 1 and 2)
            sigma(i_s, j_s) = spin_p; sigma(i_p_s, j_p_s) = spin
            do h = 9, 10
                cellspins(h) = findSpins(cellFirsts(2*h-17),cellFirsts(2*h-16))
            enddo
            sigma(i_s, j_s) = spin; sigma(i_p_s, j_p_s) = spin_p

            !And the cell-cell energy if they change places
            cellEnergy = energyCells(cellSpins)
            dF = cellEnergy(2) - cellEnergy(1)

            !Now we simply insert the acceptance criteria from the model.
            w = exp(-beta*dE)
            w_c = exp(-beta*dF)

            !Compare to a pseudo-random number between 0 and 1.
            call random_number(P)
            call random_number(P_c)
            call random_number(P_comb)

            if(dE <= 0 .and. dF <= 0) then
                sigma(i_s, j_s) = spin_p
                sigma(i_p_s, j_p_s) = spin
            elseif(dE <= 0 .and. dF > 0) then
                if(w_c >= P_c) then
                    sigma(i_s, j_s) = spin_p
                    sigma(i_p_s, j_p_s) = spin
                endif
            elseif(dE > 0 .and. dF <= 0) then
                if(w >= P) then
                    sigma(i_s, j_s) = spin_p
                    sigma(i_p_s, j_p_s) = spin
                endif
            elseif(dE > 0 .and. dF > 0) then
                if(w*w_c >= P_comb) then
                    sigma(i_s, j_s) = spin_p
                    sigma(i_p_s, j_p_s) = spin
                endif
            endif

        enddo
    
        return

    end subroutine metropolis 

end module dynamics

!The main program controls the simulation in the sense that it initialises the system, writes .dat files, informs the user about approximate time remaining etc.
program main
    use initialisation
    use dynamics
    implicit none

    !-Input parameters----------------------------------------------------------
    !Here L is the number of sites along one line of the 2D plane. So L = 128 yields a 128 x 128 spin matrix. In the future, non-square lattices should be implemented (super easy fix).
    !Lambda is the number of spins along one line in each 2D square cell (in the non-square case later, I think the cells should still be squares, but that's debatable I guess).
    !numIters represents the number of iterations between each exported frame and numFrames is the total number of frames.
    !beta is 1/kB*T, p0 and p1 are the initial concentrations of 0 and +1 spins, respectively. Currently I only set p0 and let p1 be half what's left (always 50/50 of -1 and +1 species of remainder after solvent has been introduced).
    !J_str is the matrix representation of the interaction parameters (J_00 = J_11 = J_-1-1 = 0; J_01 = J_10 = J_0-1 = J_-10 = 1; J_-11 = J_1-1 = 6 in this case).
    !sigma will be an L x L spin matrix (values -1, 0 or 1 at all places).
    !n, i and j are used in loops.
    !stat is used to remove old frames.
    !file_id and file_name hold info for saving/erasing the frames.
    !start, finish, timeLeft etc. are only used to keep track of the time and to provide some nice prompts to the user.
    integer :: i, j, stat
    character(32) :: file_id, file_name
    logical :: file_exists = .true.
    real :: start, finish, numHour, numMin, numSec, conc, MCS

    call random_seed()
    call genSpins(L, sigma, p0, p1) !Generates a pseudo-random L x L "tenary spin matrix".

    conc = real(count(sigma == 0))/real(L**2)

    !This clears the frames/ directory!
    n = 0
    do while (file_exists .eqv. .true.)
        n = n + 1
        write(file_id, '(i0)') n
        file_name = 'frames/frame-' // trim(adjustl(file_id)) // '.dat'
        inquire(file = trim(file_name), exist = file_exists)
        open(10, iostat=stat, file = trim(file_name), form = 'formatted', status = 'old')
        if(stat == 0) close(10, status='delete') 
    enddo

    !Debugging stuff - basically gives one the ability to continue an aborted simulation from a specified state (REMEMBER TO COMMENT OUT THE STUFF ABOVE - IT *WILL* REMOVE ALL OF THE FILES OTHERWISE!).
!    write(file_id, '(i0)') 170
!    file_name = 'frames/frame-' // trim(adjustl(file_id)) // '.dat'
!    open(10, file = trim(file_name), form = 'formatted')
!    do i = 1,L
!        read(10,*) (sigma(i,j), j = 1,L)
!    enddo
!    close(10)

    call cpu_time(start) !Keeps track of the time.
    !This is the main loop; it writes the .dat files and calls the metropolis() subroutine to alter the spin matrix and it also informs the user about the number of frames, time remaining etc.
    n = 1
    do while (conc >= cutoffConc)
        write(file_id, '(i0)') n
        file_name = 'frames/frame-' // trim(adjustl(file_id)) // '.dat'
        open(10, file = trim(file_name), form = 'formatted')
        do i = 1,L
            write(10,*) (sigma(i,j), j = 1,L)
        enddo
        close(10)
        
        call metropolis(L, lambda, sigma, numIters, beta) !Returns the spin matrix after numIters.

        conc = real(count(sigma == 0))/real(L**2)

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

    call cpu_time(finish)
    numHour = aint(finish/3600)
    numMin = aint((finish - (numHour)*3600)/60)
    numSec = aint((finish - (numHour)*3600 - numMin*60))
    print '(//,"Completed after: ", 12x,i2," hour(s), ",i2," minute(s) and ",i2," second(s).")'&
    ,int(numHour),int(floor(numMin)),int(floor(numSec))
    
end program main
