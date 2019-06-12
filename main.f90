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
    integer,parameter :: L = 512, lambda = 1
!    character(128) :: prefix = 'automatedRun/1024/'
!    character(128) :: prefix = 'debug/'
!    character(128) :: prefix = 'recreation/'
!    character(128) :: prefix = 'J_str/'
!    character(128) :: prefix = 'PBCvsFBC/'
!    character(128) :: prefix = 'solventDistribution/'
!    character(128) :: prefix = 'topView/'
    character(128) :: prefix = 'topView-Emilio/'

!    real,parameter :: beta = 0.6, p0 = 0.6, p1 = (1 - p0)/2, phi = 0, cutoffConc = 0.1, U = 2.0
    real,parameter :: beta = 0.6, p0 = 0.6, p1 = 0.30, phi = 0.0, cutoffConc = 0.1, U = 2.0
    logical,parameter :: constSeed = .false., FBC = .false., topView = .true., noEvap = .true.
    integer :: sigma(L,L), numSpins(L/lambda,L/lambda,1:3)
    integer, allocatable :: numIters

!    real,dimension(3, 3) :: J_str = transpose(reshape([0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str))) !J_ORIGINAL SCALED
!    real,dimension(3, 3) :: J_str = transpose(reshape([0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !2

    real,dimension(3, 3) :: J_str = transpose(reshape(0.5*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !2

!    real,dimension(3, 3) :: J_str = transpose(reshape(real(lambda)**(-2)*[0, 1, 0, 1, 0, 1, 0, 1, 0], shape(J_str))) !+1 and -1 are functionally the same.

!    real,dimension(3, 3) :: J_str = transpose(reshape(real(lambda)**(-2)*[0, 1, 2, 1, -1, 1, 2, 1, 0], shape(J_str))) !Solvent likes each other - should behave more like liquid.
!    real,dimension(3, 3) :: J_str = transpose(reshape((1/10.)*real(lambda)**(-2)*[0, 10, 20, 10, 5, 10, 20, 10, 0], shape(J_str))) !Solvent repels each other.
!    real,dimension(3, 3) :: J_str = transpose(reshape((1/10.)*real(lambda)**(-2)*[-5, 10, 20, 10, 0, 10, 20, 10, -5], shape(J_str))) !The other phases like themselves.

!    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 1, 6, 1, 0, 3, 6, 3, 0], shape(J_str))) !This will form a 'cap' of +1, cf. fig. 4 in Andrea's paper.
!    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 35, 15, 35, 0, 35, 15, 35, 0], shape(J_str))) !This is the 'strong repulsion' in Andrea's paper. It kinda works but I need much more energy..
!    integer,dimension(3, 3) :: J_str = transpose(reshape([0, 1, 6, 1, -2, 1, 6, 1, -2], shape(J_str))) !This is the 'strong repulsion' in Andrea's paper. It kinda works but I need much more energy.

end module

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

        integer :: j, i
        real :: u

        do j = 1,L
            do i = 1,L
                call random_number(u)

                if(u < p0) then
                    sigma(j,i) = 0
                elseif(u < p0 + p1) then
                    sigma(j,i) = 1
                else
                    sigma(j,i) = -1
                endif
            enddo
        enddo

        return

    end subroutine genSpins
end module initialisation

!This module contains the functions utilised by the dynamics module.
module functions
    use constants
    use strings

contains
    !This function takes a cell (or spin) number 'c' and returns the "ij"-coordinates of the originating cell as an array.
    function findCoord(c, length)result(coordResult)
        real :: c
        integer :: length, coordResult(1:2)
        integer :: j, i

        !Here we round down to the current integer (e.g. 0.8 would be 0) to find j. This first conditional is the edge-cases where e.g. on a 3 x 3 matrix, the last cell at the first row would have value 6, but then we'd get row 2 from the second conditional since 6/3 = 2.
        if(mod(int(c),length) == 0) then
            j = int(aint(c/length))
        else
            j = int(aint(c/length) + 1)
        endif

        !Now we simply determine the column i_c from
        i = int(c) - length*(j-1)

!        print *, "c = ", c
!        print *, "j = ", j
!        print *, "i = ", i

        !Note that the spin matrix is formatted according to row-column so spin(1,2) would be the spin at row 1, column 2.
        !This means that we follow the convention from the thesis - i.e. j describes rows and i describes columns, increasing from the northwestern corner.
        coordResult = [j, i]

        return
    end

    !This function takes a spin number together with the coordinates of the cell in question and finds the indices of the associated spin within the appropriate cell as an array.
    !The spin indices work with sigma(i_s, j_s).
    function spinCoord(s, j_c, i_c)result(coordResult)
        real :: s
        integer :: j_c, i_c
        integer :: coordResult(1:2)
        integer :: j_s, i_s
        integer, dimension(1:2) :: spinTemp !Temporary array for storing the coordinates and type of bond of the cells.

        !Now we re-use the earlier function
        spinTemp = findCoord(s, lambda)
        j_s = spinTemp(1); i_s = spinTemp(2)

        !Now we "transform" the i_s, j_s coordinates to reflect which cell they inhibit. Thus far the coordinates have all been relative the cell, not the whole lattice!
        j_s = j_s + lambda*(j_c - 1)
        i_s = i_s + lambda*(i_c - 1)

!        print *, "s = ", s
!        print *, "j_s = ", j_s
!        print *, "i_s = ", i_s

        coordResult = [j_s, i_s]

        return
    end

    !This function takes a set of coordinates together with the maximum dimension (e.g. the length of a box) and returns the coordinates of the nearest-neighbours with PBC as an array.
    function PBC(j, i, length)result(PBCResult)
        integer :: j, i, length
        integer, dimension(1:4) :: PBCResult
        integer :: j_u, j_d, i_l, i_r

        !Next we define some stuff to simplify the PBC (for the cells!).
        j_u = j - 1; j_d = j + 1; i_l = i - 1; i_r = i + 1

        !The first one: "if you are in the first column, the neighbour to your left will be the rightmost cell on the same row", and analogously for the rest.
        if(i == 1) then
            i_l = length
        elseif(i == length) then
            i_r = 1
        endif
        if(j == 1) then
            j_u = length
        elseif(j == length) then
            j_d = 1
        endif

        PBCResult = [j_u, i_l, j_d, i_r]

        return
    end

    !This function takes the coordinates of a cell and returns the indices of the first spin within that cell (i.e. the one in the top left-hand corner of the cell) as an array.
    function findFirst(j_c, i_c)result(firstResult)
        integer :: j_c, i_c
        integer, dimension(1:2) :: firstResult

        !Finding the first coordinates is trivial. For j, simply take the number of rows to a cell and multiply that by the cell index minus one (the minus one is due to the cell index starting at 1 rather than 0) and add 1 because otherwise you get the last row in the cell above, and analogous for i.
        j_first = 1 + lambda*(j_c - 1)
        i_first = 1 + lambda*(i_c - 1)

        firstResult = [j_first, i_first]

        return
    end

    !This function takes the value of L and lambda as well as the initial spin configuration and returns the number of spins of a given species in a rank 3 tensor, where the first the indices relate to the row-column of the cell and the third to the spin species (-1,0,1)+2.
    function initialNum(L, lambda, sigma)result(numResult)
        integer :: L, lambda, sigma(L,L)
        integer, dimension(1:3) :: tempSpins
        integer, dimension(L/lambda,L/lambda,1:3) :: numResult

        integer :: j, i, first(1:2), x2, x1

        !The outer loop runs over the total number of cells.
        do j = 1, L/lambda
            do i = 1, L/lambda
                first = findFirst(j, i)

                !Reset for each run.
                tempSpins = [0,0,0]

                !Here we simply loop over each spin in the cell and count the number of spins of species \alpha \in {-1, 0, 1}.
                do x2 = first(1), first(1) + lambda - 1
                    do x1 = first(2), first(2) + lambda - 1
                        spin_loop = sigma(x2, x1)
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
                    numResult(j,i,k) = tempSpins(k)
                enddo

            enddo
        enddo

        return
    end

    !This function takes the value of the currently considered spin together with the value of the proposed spin switch, as well as the coordinates of the first spin in the relevant cell and the indices of the spin in question and returns the current energy, along with the proposed energy should the spins be switched in an array.
    function energySelected(spin, spin_p, t, numSpins, j_c, i_c, j_c_p, i_c_p, cellNeighbours)result(energyResult)
        integer :: spin, spin_p, t, numSpins(L/lambda,L/lambda,1:3), j_c, i_c, j_c_p, i_c_p, cellNeighbours(1:12) !Needed to specify the type - Fortran gets real mad if you try and use non-integers for indices (and it seems to always assume real for non-specified types).
        integer :: k, b, h
        integer, dimension(1:2) :: coordNum
        integer, dimension(1:3) :: numSpinsA, numSpinsB, numSpinsAprop, numSpinsBprop
        integer, dimension(1:3) :: numSpinsIntA, numSpinsIntB
        real :: E_current, E_proposed, energyResult(1:2)
        real, dimension(1:2) :: c

        !Reset the values before each run.
        E_current = 0; E_proposed = 0; coordNum = [4, 4]

        !The number of spins in cell A and B can be trivially determined using numSpins, together with the number of spins of species k (the currently considered spin) interacts with (just the number of spins in that cell of the species, minus one to avoid self-interaction).
        do k = 1,3
            numSpinsA(k) = numSpins(j_c,i_c,k)
            numSpinsB(k) = numSpins(j_c_p,i_c_p,k)
            numSpinsAprop(k) = numSpins(j_c,i_c,k)
            numSpinsBprop(k) = numSpins(j_c_p,i_c_p,k)
            numSpinsIntA(k) = numSpins(j_c,i_c,k)
            numSpinsIntB(k) = numSpins(j_c_p,i_c_p,k)
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

        !The internal energy is determined using the number of spins per cell, since their location doesn't matter in mean-field.
        do k = 1, 3
            E_current = E_current + numSpinsIntA(k)*J_str(spin+2,k)*U + numSpinsIntB(k)*J_str(spin_p+2,k)*U
            E_proposed = E_proposed + numSpinsIntA(k)*J_str(spin_p+2,k)*U + numSpinsIntB(k)*J_str(spin+2,k)*U
        enddo

        !And the nearest-neighbouring cell-cell interaction (note that the FBC changes the calculations with the neighbours slightly).
        !This section is hard to follow without the aid of the figure showing the nearest-neighbour cluster - see the thesis.

        !First we adjust the coordNum.
        do h = 1, 11, 2
            if(FBC .eqv. .true.) then
                if(topView .eqv. .false.) then
                    if(j_c == L/lambda) then
                        if(t == 1) then
                            if(h == 3) coordNum(1) = coordNum(1) - 1

                        elseif(t == 2) then
                            if(h == 5) coordNum(2) = coordNum(2) - 1
                            if(h == 7) coordNum(1) = coordNum(1) - 1

                        elseif(t == 3) then

                        elseif(t == 4) then
                            if(h == 5) coordNum(1) = coordNum(1) - 1
                            if(h == 7) coordNum(2) = coordNum(2) - 1

                        endif

                    elseif(j_c == 1) then
                        if(t == 1) then

                        elseif(t == 2) then
                            if(h == 1) coordNum(2) = coordNum(2) - 1
                            if(h == 11) coordNum(1) = coordNum(1) - 1

                        elseif(t == 3) then
                            if(h == 9) coordNum(1) = coordNum(1) - 1

                        elseif(t == 4) then
                            if(h == 1) coordNum(1) = coordNum(1) - 1
                            if(h == 11) coordNum(2) = coordNum(2) - 1

                        endif

                    elseif(j_c == L/lambda - 1 .and. t == 3) then
                        if(h == 3) coordNum(2) = coordNum(2) - 1

                    elseif(j_c == 1+1 .and. t == 1) then
                        if(h == 9) coordNum(2) = coordNum(2) - 1

                    endif
                endif
            endif
        enddo

        !Adjust the scaling factor w.r.t. the coordination number.
        do k = 1,2
            c(k) = 4./coordNum(k)
        enddo

        !The energy between '1' and '2'.
        do k = 1, 3
            do b = 1, 3
                E_current = E_current + (c(1)/c(2))*numSpinsA(k)*numSpinsB(b)*J_str(k,b)
                E_proposed = E_proposed + (c(2)/c(1))*numSpinsAprop(k)*numSpinsBprop(b)*J_str(k,b)
            enddo
        enddo

        !The energy between '1' and '3', '4', '5'.
        do h = 1, 5, 2
            if(FBC .eqv. .true.) then
                if(topView .eqv. .false.) then
                    if(j_c == L/lambda) then
                        if(t == 1) then
                            if(h == 3) GO TO 15

                        elseif(t == 2) then
                            if(h == 5) GO TO 15

                        elseif(t == 3) then

                        elseif(t == 4) then
                            if(h == 5) GO TO 15

                        endif

                    elseif(j_c == 1) then
                        if(t == 1) then

                        elseif(t == 2) then
                            if(h == 1) GO TO 15

                        elseif(t == 3) then

                        elseif(t == 4) then
                            if(h == 1) GO TO 15

                        endif

                    elseif(j_c == L/lambda - 1 .and. t == 3) then
                        if(h == 3) GO TO 15

                    elseif(j_c == 1+1 .and. t == 1) then

                    endif
                endif
            endif
            do k = 1, 3
                do b = 1, 3
                    if(t == 1 .or. t == 4) then
                        E_current = E_current + c(1)*numSpinsA(k)*numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                        E_proposed = E_proposed + c(1)*numSpinsAprop(k)* &
                        numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                    else
                        E_current = E_current + c(2)*numSpinsB(k)*numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                        E_proposed = E_proposed + c(2)*numSpinsBprop(k)* &
                        numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                    endif
                enddo
            enddo
15          continue
        enddo


        !The energy between '2' and '6', '7', '8'.
        do h = 7, 11, 2
            if(FBC .eqv. .true.) then
                if(topView .eqv. .false.) then
                    if(j_c == L/lambda) then
                        if(t == 1) then

                        elseif(t == 2) then
                            if(h == 7) GO TO 25

                        elseif(t == 3) then

                        elseif(t == 4) then
                            if(h == 7) GO TO 25

                        endif

                    elseif(j_c == 1) then
                        if(t == 1) then

                        elseif(t == 2) then
                            if(h == 11) GO TO 25

                        elseif(t == 3) then
                            if(h == 9) GO TO 25

                        elseif(t == 4) then
                            if(h == 11) GO TO 25

                        endif

                    elseif(j_c == L/lambda - 1 .and. t == 3) then
                    elseif(j_c == 1+1 .and. t == 1) then
                        if(h == 9) GO TO 25
                    endif
                endif
            endif
            do k = 1, 3
                do b = 1, 3
                    if(t == 1 .or. t == 4) then
                        E_current = E_current + c(2)*numSpinsB(k)*numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                        E_proposed = E_proposed + c(2)*numSpinsBprop(k)* &
                        numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                    else
                        E_current = E_current + c(1)*numSpinsA(k)*numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                        E_proposed = E_proposed + c(1)*numSpinsAprop(k)* &
                        numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                    endif
                enddo
            enddo
25          continue
        enddo

        energyResult = [E_current, E_proposed] !Return the energy.

        return
    end

end module functions

!This module contains the dynamics in the simulation.
module dynamics
    use constants
    use functions

contains
    !This subroutine handles the evaporation of the top row by replacing the spin at (i,j) by +1 or -1 according to the appropriate concentration and returning the spin matrix.
    subroutine evap(j_s, i_s, j_c, i_c, sigma, numSpins)
        integer :: j_s, i_s, j_c, i_c, sigma(L,L), numSpins(L/lambda, L/lambda, 1:3)
        !Simply get a random number between 0 and 1 and compare it to the normalised concentrations, preserving the initial ratios of -1 and +1 species.
        call random_number(P)
        numSpins(j_c,i_c,0+2) = numSpins(j_c,i_c,0+2) - 1
        if(P < p1/(1 - p0)) then
            sigma(j_s,i_s) = 1
            numSpins(j_c,i_c,1+2) = numSpins(j_c,i_c,1+2) + 1
        else
            sigma(j_s,i_s) = -1
            numSpins(j_c,i_c,-1+2) = numSpins(j_c,i_c,-1+2) + 1
        endif
    end subroutine evap

    !This subroutine updates sigma, should the proposed move be approved (trivial, but convenient to have as a subroutine).
    subroutine updateSigma(spin, spin_p, sigma, j_s, i_s, j_s_p, i_s_p)
        integer :: spin, spin_p, sigma(L,L), j_s, i_s, j_s_p, i_s_p

        !Simply switch the spins at i_s, j_s and i_s_p, j_s_p via the spin and spin_p variables.
        sigma(j_s, i_s) = spin_p
        sigma(j_s_p, i_s_p) = spin

        return
    end subroutine updateSigma

    !This subroutine 'recalculates' the numSpins tensor (really it just updates it - there is no need for a bunch of calculations).
    subroutine recalcSpins(spin, spin_p, numSpins, j_c, i_c, j_c_p, i_c_p)
        integer :: spin, spin_p, numSpins(L/lambda,L/lambda,1:3), j_c, i_c, j_c_p, i_c_p
        integer :: k

        !Simply loop over the number of spin species and adjust accordingly.
        do k = 1,3
            if(k == spin + 2) then
                numSpins(j_c,i_c,k) = numSpins(j_c,i_c,k) - 1
                numSpins(j_c_p,i_c_p,k) = numSpins(j_c_p,i_c_p,k) + 1
            endif
            if(k == spin_p + 2) then
                numSpins(j_c,i_c,k) = numSpins(j_c,i_c,k) + 1
                numSpins(j_c_p,i_c_p,k) = numSpins(j_c_p,i_c_p,k) - 1
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
        integer :: j_c, i_c, k, spin, j_u, j_d, i_l, i_r, t
        integer :: j_c_p, i_c_p, spin_p
        integer :: j_s, i_s, j_s_p, i_s_p, j_u_p, j_d_p, i_l_p, i_r_p
        integer, dimension(1:4) :: PBCTemp !Temporary array for storing the coordinates from the PBC function.
        integer, dimension(1:2) :: spinTemp !Temporary array for storing the coordinates of the spin.
        integer, dimension(1:2) :: cellTemp !Temporary array for storing the coordinates and type of bond of the cells.
        integer, dimension(1:12) :: cellNeighbours !Array to store the coordinates of the 6 neighbours (see image in thesis).
        real, dimension(1:2) :: energy !Array for storing the energy before and after proposed spin switch.
        real :: c, s, u, dE, w, P

        !This is the actual Metropolis algorithm. Again, somewhat messy.
        do k = 1,numIters
            !-Cell and spin selection-------------------------------------------
!            !First we randomly select a bond.
!10          call random_number(u); b = 1 + floor(2*((L/lambda)**2)*u) !This yields an integer between 1 and 2*(L/lambda)^2 (i.e. the total number of bonds between cells). Whence b will be the currently considered bond number.

!            !Then we determine the "ij-coordinates" of the origin cell of the bond.
!            cell1 = cellCord(b)
!            i_c = cell1(1); j_c = cell1(2); t = cell1(3)

            !First we randomly select a cell.
10          call random_number(u); c = 1 + floor(((L/lambda)**2)*u) !This yields an integer between 1 and (L/lambda)^2 (i.e. the total number of cells). Whence c will be the currently considered cell number.

            !Then we determine the "ij-coordinates" of the origin cell of the bond (note that Fortran uses row-major order, i.e. sigma(j,i), where j is the row).
            cellTemp = findCoord(c, L/lambda)
            j_c = cellTemp(1); i_c = cellTemp(2)

            !Now we randomly select a spin inside of the first cell of the bond.
            call random_number(u); s = 1 + floor((lambda**2)*u) !This yields an integer between 1 and lambda^2 (i.e. the total number of spins in a cell).

            !And determine the "x1x2-coordinates" of the spin (again, note that we have row-column, so x2 is the row in this case, i.e. sigma(x1,x2)).
            spinTemp = spinCoord(s, j_c, i_c)
            j_s = spinTemp(1); i_s = spinTemp(2)

            !Finally we have the actual spin at (i_s,j_s), which we store in spin.
            spin = sigma(j_s, i_s)

            !Now we call the PBC() function which returns a 2 x 2 matrix with the "ij-cordinates" of the nearest-neighbours.
            PBCtemp = PBC(j_c, i_c, L/lambda)
            j_u = PBCTemp(1); i_l = PBCTemp(2); j_d = PBCTemp(3); i_r = PBCTemp(4)

            !Now we decide on a nearest-neighbouring cell.
            call random_number(u); t = 1 + floor(4*u) !This yields an integer between 1 and 4 (i.e. the number of nearest-neighbours). Whence t will be the direction where t goes from t = 1 at 12 o-clock to t = 4 to the right.

            !Next we find the position of the other cell associated with the bond, which we have denoted _p for 'pair'.
            if(t == 1) then
                j_c_p = j_u; i_c_p = i_c
            elseif(t == 2) then
                j_c_p = j_c; i_c_p = i_l
            elseif(t == 3) then
                j_c_p = j_d; i_c_p = i_c
            else
                j_c_p = j_c; i_c_p = i_r
            endif

            !And then we re-run PBC() to find the nearest-neighbours. In principle, this could be hard-coded, but it's so quick that it doesn't really matter.
            PBCtemp = PBC(j_c_p, i_c_p, L/lambda)
            j_u_p = PBCTemp(1); i_l_p = PBCTemp(2); j_d_p = PBCTemp(3); i_r_p = PBCTemp(4)

            !Now we randomly select a spin inside of the paired cell of the bond.
            call random_number(u); s = 1 + floor((lambda**2)*u) !This yields an integer between 1 and lambda^2 (i.e. the total number of spins in a cell).

            !And determine the "x1x2-coordinates" of the spin.
            spinTemp = spinCoord(s, j_c_p, i_c_p)
            j_s_p = spinTemp(1); i_s_p = spinTemp(2)

            spin_p = sigma(j_s_p, i_s_p) !Now we have the actual spin at (i_s_p,j_s_p), which we store in spin_p.

            !-Evaporation-------------------------------------------------------
            if(topView .eqv. .false.) then
                if(noEvap .eqv. .false.) then
                    if(j_c == 1 .and. spin == 0 .and. t == 1) then
                        call evap(j_s, i_s, j_c, i_c, sigma, numSpins)
                        GO TO 10
                    endif
                endif
            elseif(topView .eqv. .true.) then
                if(noEvap .eqv. .false.) then
                    if(spin == 0) then
                        call random_number(P) !Compare to a pseudo-random number between 0 and 1.
                        if(P < 0.0002) call evap(j_s, i_s, j_c, i_c, sigma, numSpins)
!                        if(P < 0.002) call evap(j_s, i_s, j_c, i_c, sigma, numSpins) !DEBUG
                        GO TO 10
                    endif
                endif
            endif

            !-Avoid moves across top/bottom boundary----------------------------
            if(topView .eqv. .false.) then
                if(j_c == 1 .and. t == 1) then
                    GO TO 10
                elseif(j_c == L/lambda .and. t == 3) then
                    GO TO 10
                endif
            endif

            !-Dynamics----------------------------------------------------------
            !Before we compute the energy, we must find the nearest-neighbouring cells.
            if(t == 1) then
                cellNeighbours = [[j_c,i_l], [j_d,i_c], [j_c,i_r], [j_c_p,i_r_p], [j_u_p,i_c_p], [j_c_p,i_l_p]]
            elseif(t == 2) then
                cellNeighbours = [[j_u_p,i_c_p], [j_c_p,i_l_p], [j_d_p,i_c_p], [j_d,i_c], [j_c,i_r], [j_u,i_c]]
            elseif(t == 3) then
                cellNeighbours = [[j_c_p,i_l_p], [j_d_p,i_c_p], [j_c_p,i_r_p], [j_c,i_r], [j_u,i_c], [j_c,i_l]]
            elseif(t == 4) then
                cellNeighbours = [[j_u,i_c], [j_c,i_l], [j_d,i_c], [j_d_p,i_c_p], [j_c_p,i_r_p], [j_u_p,i_c_p]]
            endif

            !Now we simply call the energySelected function, with the needed information and BAM! Just like that, we're done. Thanks, objective oriented programming.
            energy = energySelected(spin, spin_p, t, numSpins, j_c, i_c, j_c_p, i_c_p, cellNeighbours)
            dE = energy(2) - energy(1)

!            print *, dE

            !Now we simply insert the acceptance criteria from the model.
            w = exp(-beta*dE*real(lambda)**(-2))

            !Compare to a pseudo-random number between 0 and 1.
            call random_number(P)

            !And switch spins if it is energetically favourable or uphill if P <= w.
            if(dE < 0 .or. P <= w) then
                call updateSigma(spin, spin_p, sigma, j_s, i_s, j_s_p, i_s_p)
                call recalcSpins(spin, spin_p, numSpins, j_c, i_c, j_c_p, i_c_p)
            endif

            !-Solvent volatility------------------------------------------------
            !This is the upwards drift due to the volatility.
            if(phi == 0) GO TO 30
            spin = sigma(j_s, i_s)
            spin_p = sigma(j_s_p, i_s_p)
            if(t == 1 .and. spin == 0 .and. spin_p /= 0 .or. t == 3 .and. spin_p == 0 .and. spin /= 0) then
                call random_number(P) !Compare to a pseudo-random number between 0 and 1.
                if(P < phi) then
                    call updateSigma(spin, spin_p, sigma, j_s, i_s, j_s_p, i_s_p)
                    call recalcSpins(spin, spin_p, numSpins, j_c, i_c, j_c_p, i_c_p)
                endif
            endif
30          continue

        enddo

        return

    end subroutine metropolis

end module dynamics

!The main program controls the simulation in the sense that it initialises the system, writes .dat files, informs the user about approximate time remaining etc.
program main
    use strings
    use initialisation
    use functions
    use dynamics
    use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
    implicit none

    !-Input parameters----------------------------------------------------------
    !n, i and j are used in loops; stat is used to remove old frames.
    !file_id and file_name hold info for saving/erasing the frames.
    !start, finish, etc. are only used to keep track of the time and to provide some nice prompts to the user; conc is the current concentration of zeros; MCS is the current number of MCS (for prompt).
    integer :: j, i, n, stat, nI
    character(128) :: file_id, file_name, folderName
    character :: d, d2
    logical :: file_exists = .true., folder_exists = .true.
    real :: start, finish, numHour, numMin, numSec, conc, MCS
    INTEGER (KIND=2) :: expon

    if(L == 1024) then
        expon = nint(16.9955*(L/lambda)**(0.1102))
    elseif(L == 512) then
        expon = nint(13.8647*(L/lambda)**(0.1309))
        expon = 26 !HARDCODED
    elseif(L == 256) then
        expon = nint(11.3520*(L/lambda)**(0.1621))
    elseif(L == 128) then
        expon = nint(7.1576*(L/lambda)**(0.2355))
    else
        print '("Please fit numIters data. Using L = 512 values...",/,"Press ENTER to continue.")'
        read(stdin,*)
        expon = nint(13.8647*(L/lambda)**(0.1309))
    endif
    if(expon > 30) then
        print '("Integer overflow - using maximum numIters exponent of 30...",/,"Press ENTER to continue.")'
        read(stdin,*)
        expon = 30
    endif
    numIters = 2**(expon)

    !This creates a directory with the correct name (if it does not currently exist).
    nI = nint(log(real(numIters))/log(2.))
    folderName = 'lambda_' // trim(adjustl(num2str(lambda))) // '-L_' // trim(adjustl(num2str(L))) // '-J_' // &
        trim(adjustl(num2str(J_str(1,1)))) // '_' // trim(adjustl(num2str(J_str(1,2)))) // '_' // &
        trim(adjustl(num2str(J_str(1,3)))) // '-numIters_2-' // trim(adjustl(num2str(nI))) // '-initialDist_' // &
        trim(adjustl(num2str(nint(p0*100.)))) // '_' // trim(adjustl(num2str(nint(p1*100)))) // '_' // &
        trim(adjustl(num2str(nint((1 - p0 - p1)*100)))) // '-U_' // trim(adjustl(num2str(U)))
    if(topView .eqv. .false.) then
        if(FBC .eqv. .true.) folderName = trim(folderName) // '-FBC'
        if(FBC .eqv. .false.) folderName = trim(folderName) // '-PBC'
    else
        folderName = trim(folderName) // '-topView'
    endif
    folderName = trim(prefix) // trim(folderName)
    inquire(file = './' // trim(folderName), exist = folder_exists)
    if(folder_exists .eqv. .true.) GO TO 70
    print '(//,"Creating directory...")'
    call execute_command_line ('mkdir ' // prefix)
    call execute_command_line ('mkdir ' // folderName)

    !This clears the directory using magic.
70  n = 0
    do while (file_exists .eqv. .true.)
        n = n + 1
        write(file_id, '(i0)') n
        file_name = trim(adjustl(folderName)) // '/frame-' // trim(adjustl(file_id)) // '.dat'
        inquire(file = trim(file_name), exist = file_exists)
        if(n == 1 .and. file_exists .eqv. .true.) then
            print '(//)'; print *, 'The directroy ' // trim(folderName) // ' is not empty!'
            print '(">Delete old files? (y/n)")'; read *, d
            if(d == 'y') then
                print '(//,"Deleting files...")'
                GO TO 80
            else
                print '(//,"Continue simulation from last frame? (y/n)")'; read *, d2
                if(d2 == 'y') then
                    EXIT
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

    !This starts the simulation from the last frame in the directory.
    if(d2 == 'y') then
        !First we find the correct 'n' for the latest frame.
        n = 0
        do while (file_exists .eqv. .true.)
            n = n + 1
            write(file_id, '(i0)') n
            file_name = trim(adjustl(folderName)) // '/frame-' // trim(adjustl(file_id)) // '.dat'
            inquire(file = trim(file_name), exist = file_exists)
        enddo

        !Then we read the frame into memory.
        n = n-1
        write(file_id, '(i0)') n
        file_name = trim(adjustl(folderName)) // '/frame-' // trim(adjustl(file_id)) // '.dat'
        open(10, file = trim(file_name), form = 'formatted')
        do j = 1,L
            read(10,*) (sigma(j,i), i = 1,L)
        enddo
        close(10)

        !Then we setup the RNG and start the dynamics.
        call setupRNG()
        GO TO 85
    endif

    call setupRNG() !Sets up the RNG (either with a constant seed or using UNIX time for m-m-m-aximum pseudo-randomness!).
    call genSpins(L, sigma, p0, p1) !Generates a pseudo-random L x L "tenary spin matrix".

85  numSpins = initialNum(L, lambda, sigma) !Determines the initial number of spins.
    conc = real(count(sigma == 0))/real(L**2) !Calculates the intial concentration.

    print '(//,"Starting dynamics...")'
    call cpu_time(start) !Keeps track of the time.
    !This is the main loop; it writes the .dat files and calls the metropolis() subroutine to alter the spin matrix and it also informs the user about the number of frames, time remaining etc.
    if(d2 /= 'y') n = 1
    do while (conc > cutoffConc)
        write(file_id, '(i0)') n
        file_name = trim(adjustl(folderName)) // '/frame-' // trim(adjustl(file_id)) // '.dat'
        open(10, file = trim(file_name), form = 'formatted')
        do j = 1,L
            write(10,*) (sigma(j,i), i = 1,L)
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
