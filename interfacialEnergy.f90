!Author: simokron
!
!This is a complementary program to the main one as part of my master's thesis.

!This module contains the constants used in the calculations.
module constants
    implicit none

    integer :: L, lambda
    real :: J_str(3,3)
    real, parameter :: beta = 0.6, cutoffConc = 0.1
    logical, parameter :: FBC = .true., topView = .false.
    integer, dimension (:,:), allocatable :: sigma
    real, dimension (:,:), allocatable :: energyCells
    integer, dimension (:,:,:), allocatable :: numSpins
    character(256), allocatable :: dir
!    character(256) :: prefix = 'automatedRun/128/'
!    character(256) :: folder = 'lambda_2-L_128-J_0.0000_1.0000_2.0000-numIters_2-19-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_4-L_128-J_0.0000_1.0000_2.0000-numIters_2-16-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_8-L_128-J_0.0000_1.0000_2.0000-numIters_2-14-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_16-L_128-J_0.0000_1.0000_2.0000-numIters_2-12-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_32-L_128-J_0.0000_1.0000_2.0000-numIters_2-10-initialDist_60_20_20-FBC'

    character(256) :: prefix = 'automatedRun/256/'
!    character(256) :: folder = 'lambda_2-L_256-J_0.0000_1.0000_2.0000-numIters_2-25-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_4-L_256-J_0.0000_1.0000_2.0000-numIters_2-22-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_8-L_256-J_0.0000_1.0000_2.0000-numIters_2-20-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_16-L_256-J_0.0000_1.0000_2.0000-numIters_2-18-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_32-L_256-J_0.0000_1.0000_2.0000-numIters_2-16-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_64-L_256-J_0.0000_1.0000_2.0000-numIters_2-14-initialDist_60_20_20-FBC'

!    character(256) :: folder = 'lambda_2-L_256-J_0.0000_1.0000_2.0000-numIters_2-25-initialDist_60_20_20-PBC'
!    character(256) :: folder = 'lambda_4-L_256-J_0.0000_1.0000_2.0000-numIters_2-22-initialDist_60_20_20-PBC'
!    character(256) :: folder = 'lambda_8-L_256-J_0.0000_1.0000_2.0000-numIters_2-20-initialDist_60_20_20-PBC'
!    character(256) :: folder = 'lambda_16-L_256-J_0.0000_1.0000_2.0000-numIters_2-18-initialDist_60_20_20-PBC'
!    character(256) :: folder = 'lambda_32-L_256-J_0.0000_1.0000_2.0000-numIters_2-16-initialDist_60_20_20-PBC'
    character(256) :: folder = 'lambda_64-L_256-J_0.0000_1.0000_2.0000-numIters_2-14-initialDist_60_20_20-PBC'

!    character(256) :: prefix = 'automatedRun/512/'
!    character(256) :: folder = 'lambda_2-L_512-J_0.0000_1.0000_2.0000-numIters_2-29-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_4-L_512-J_0.0000_1.0000_2.0000-numIters_2-26-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_8-L_512-J_0.0000_1.0000_2.0000-numIters_2-24-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_16-L_512-J_0.0000_1.0000_2.0000-numIters_2-22-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_32-L_512-J_0.0000_1.0000_2.0000-numIters_2-20-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_64-L_512-J_0.0000_1.0000_2.0000-numIters_2-18-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_128-L_512-J_0.0000_1.0000_2.0000-numIters_2-17-initialDist_60_20_20-FBC'

!    character(256) :: prefix = 'automatedRun/1024/'
!    character(256) :: folder = 'lambda_2-L_1024-J_0.0000_1.0000_2.0000-numIters_2-30-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_4-L_1024-J_0.0000_1.0000_2.0000-numIters_2-30-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_8-L_1024-J_0.0000_1.0000_2.0000-numIters_2-29-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_16-L_1024-J_0.0000_1.0000_2.0000-numIters_2-27-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_32-L_1024-J_0.0000_1.0000_2.0000-numIters_2-25-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_64-L_1024-J_0.0000_1.0000_2.0000-numIters_2-23-initialDist_60_20_20-FBC'
!    character(256) :: folder = 'lambda_128-L_1024-J_0.0000_1.0000_2.0000-numIters_2-21-initialDist_60_20_20-FBC'


!    character(256) :: prefix = 'PBCvsFBC/'
!    character(256) :: folder = 'lambda_4-L_256-J_0.0000_1.0000_2.0000-numIters_2-22-initialDist_60_20_20-PBC'
!    character(256) :: folder = 'lambda_4-L_256-J_0.0000_1.0000_2.0000-numIters_2-22-initialDist_60_20_20-FBC'

end module

!This module sets up the parameters.
module initialisation
    use constants

contains
    !This subroutine finds the parameters from the dir name.
    subroutine findParamaters(dir, lambda, L)
        character(256) :: dir
        integer :: lambda, L

        character(256) :: tStr
        integer :: ipos, iposLim, t(1:16)

        ipos = index(dir,"lambda_")
        iposLim = index(trim(dir(ipos +len("lambda_"):)),"-L")
        tStr = trim(dir(ipos + len("lambda_"):ipos + len("lambda_") + iposLim - 2))
        read(tStr(1:),'(i6)') t(1)
        lambda = t(1)

        ipos = index(dir,"L_")
        iposLim = index(trim(dir(ipos +len("L_"):)),"-J")
        tStr = trim(dir(ipos + len("L_"):ipos + len("L_") + iposLim - 2))
        read(tStr(1:),'(i6)') t(1)
        L = t(1)

!        ipos = index(dir,"numIters_")
!        tStr = trim(dir(ipos +len("numIters_"):ipos + len("numIters_")+3))
!        read(tStr(1:1),'(i6)') t(1)
!        read(tStr(3:),'(i6)') t(2)
!        numIters = t(1)**t(2)
    end subroutine findParamaters
end module

!This module contains the functions utilised throughout.
module functions
    use constants

contains
    !This function takes a set of coordinates together with the maximum dimension (e.g. the length of a box) and returns the coordinates of the nearest-neighbours with PBC as an array.
    function PBC(j, i, length)result(PBCResult)
        integer :: j, i, length
        integer, dimension(1:4) :: PBCResult
        integer :: j_u, j_d, i_l, i_r

        !Next we define some stuff to simplify the PBC (for the cells!).
        j_u = j - 1; j_d = j + 1; i_l = i - 1; i_r = i + 1

        !The first one: "if you are in the first column, the neighbour to your left will be the rightmost cell on the same row", and analogously for the rest
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

    function calcEnergy(sigma)result(energyResult)
        integer :: sigma(L,L)

        integer :: numSpins(L/lambda,L/lambda,1:3), coordNum
        integer :: j_c, i_c, j_d, j_u, i_l, i_r, k, b, h
        real :: E_current
        integer, dimension(1:4) :: PBCTemp
        integer, dimension(1:8) :: cellNeighbours
        real :: energyResult(L/lambda,L/lambda)

        numSpins = initialNum(L, lambda, sigma) !Determines the initial number of spins.

        do j_c = 1, L/lambda
            do i_c = 1, L/lambda
                PBCtemp = PBC(j_c, i_c, L/lambda)
                j_u = PBCTemp(1); i_l = PBCTemp(2); j_d = PBCTemp(3); i_r = PBCTemp(4)

                cellNeighbours = [[j_u, i_c], [j_c, i_l], [j_d, i_c], [j_c, i_r]]

                !Reset for each cell.
                E_current = 0; coordNum = 4

                do h = 1, 7, 2
                    if(FBC .eqv. .true. .and. topView .eqv. .false.) then
                        if(j_c == L/lambda .and. h == 5) then
                            coordNum = coordNum - 1
                        elseif(j_c == 1 .and. h == 1) then
                            coordNum = coordNum - 1
                        endif
                    endif

                    c = 4./coordNum
                enddo

                do h = 1, 7, 2
                    if(FBC .eqv. .true. .and. topView .eqv. .false.) then
                        if(j_c == L/lambda .and. h == 5) then
                            GO TO 20
                        elseif(j_c == 1 .and. h == 1) then
                            GO TO 20
                        endif
                    endif
                    do k = 1, 3
                        do b = 1, 3
                            E_current = E_current + c*numSpins(j_c, i_c, k)*&
                                numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)/2
                        enddo
                    enddo
        20          continue
                enddo

                energyResult(j_c, i_c) = E_current
            enddo
        enddo

        return
    end

end module functions

!This module contains the subroutines utilised throughout.
module subroutines
    use constants
    use functions

contains
    subroutine nFinder(dir, n, L)
        character(256) :: dir
        integer :: n, L

        character(256) :: file_id, file_name
        logical :: file_exists = .true.

        n = 0
        do while (file_exists .eqv. .true.)
            n = n + 1
            write(file_id, '(i0)') n
            file_name = trim(dir) // '/frame-' // trim(adjustl(file_id)) // '.dat'
            inquire(file = trim(file_name), exist = file_exists)
        enddo
        n = n - 1

        do while (conc < cutoffConc)
            write(file_id, '(i0)') n
            file_name = trim(dir) // '/frame-' // trim(adjustl(file_id)) // '.dat'
            open(10, file = trim(file_name), form = 'formatted')
            do j = 1,L
                read(10,*) (sigma(j,i), i = 1,L)
            enddo
            close(10)

            conc = real(count(sigma == 0))/real(L**2)
            n = n - 1
        enddo
        n = n + 1

    end subroutine nFinder

end module

!This is the main program.
program interfacialEnergy
    use constants
    use initialisation
    use subroutines
    implicit none

    integer :: n, j, i
    character(256) :: file_id, file_name

    dir = trim(prefix) // trim(folder)

    call findParamaters(dir, lambda, L)
    allocate(sigma(L,L))
    allocate(numSpins(L/lambda,L/lambda,1:3))
    allocate(energyCells(L/lambda,L/lambda))

    J_str = transpose(reshape(real(lambda)**(-2)*[0, 1, 2, 1, 0, 1, 2, 1, 0], shape(J_str))) !2
!    J_str = transpose(reshape([0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str)))

    call nFinder(dir, n, L)

    write(file_id, '(i0)') n
    file_name = trim(dir) // '/frame-' // trim(adjustl(file_id)) // '.dat'
    open(10, file = trim(file_name), form = 'formatted')
    do j = 1,L
        read(10,*) (sigma(j,i), i = 1,L)
    enddo
    close(10)

    energyCells = calcEnergy(sigma)
    print *, 'H_int = ', sum(energyCells)/L**2

end program interfacialEnergy
