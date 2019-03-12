!Author: simokron
!
!This is a complementary program to the main one as part of my master's thesis.

!This module contains the constants used in the calculations.
module constants
    implicit none
    
    integer :: L, lambda
    real :: J_str(3,3)
    real, parameter :: beta = 0.6, cutoffConc = 0.1
    logical,parameter :: FBC = .true.
    integer, dimension (:,:), allocatable :: sigma
    real, dimension (:,:), allocatable :: energyCells
    integer, dimension (:,:,:), allocatable :: numSpins
!    character(256) :: dir = 'evapRand/se_lambda_4-L_1024_frames_v3-c_1-J_squareScaled-numIters_2-27-FBC_evapRand'
!    character(256) :: dir = 'evapRand/se_lambda_8-L_1024_frames_v3-c_1-J_squareScaled-numIters_2-26-FBC_evapRand'
!    character(256) :: dir = 'evapRand/se_lambda_16-L_1024_frames_v3-c_1-J_squareScaled-numIters_2-24-FBC_evapRand'
!    character(256) :: dir = 'evapRand/se_lambda_32-L_1024_frames_v3-c_1-J_squareScaled-numIters_2-22-FBC_evapRand'
!    character(256) :: dir = 'evapRand/se_lambda_64-L_1024_frames_v3-c_1-J_squareScaled-numIters_2-20-FBC_evapRand'
    character(256) :: dir = 'evapRand/se_lambda_128-L_1024_frames_v3-c_1-J_squareScaled-numIters_2-18-FBC_evapRand'

!    character(256) :: dir = 's_lambda_4-L_512_frames_v3-c_1-J_squareScaled-numIters_2-24-FBC'
!    character(256) :: dir = 's_lambda_8-L_512_frames_v3-c_1-J_squareScaled-numIters_2-24-FBC'
!    character(256) :: dir = 's_lambda_16-L_512_frames_v3-c_1-J_squareScaled-numIters_2-23-FBC'
!    character(256) :: dir = 's_lambda_32-L_512_frames_v3-c_1-J_squareScaled-numIters_2-20-FBC'
!    character(256) :: dir = 's_lambda_64-L_512_frames_v3-c_1-J_squareScaled-numIters_2-18-FBC'

!    character(256) :: dir = 'p0_04/s_lambda_8-L_512_frames_v3-c_1-J_squareScaled-numIters_2-24-FBC'
!    character(256) :: dir = 'p0_06/s_lambda_8-L_512_frames_v3-c_1-J_squareScaled-numIters_2-24-FBC'

!    character(256) :: dir = ''
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
        tStr = trim(dir(ipos +len("lambda_"):ipos + len("lambda_") + iposLim - 2))
        read(tStr(1:),'(i6)') t(1)
        lambda = t(1)
        
        ipos = index(dir,"L_")
        iposLim = index(trim(dir(ipos +len("L_"):)),"_frames")
        tStr = trim(dir(ipos +len("L_"):ipos + len("L_") + iposLim - 2))
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

    !This function takes the value of L and lambda as well as the initial spin configuration and returns the number of spins of a given species in a 3D array, where the first the indices relate to the row-column of the cell and the third to the spin species (-1,0,1)+2.
    function initialNum(L, lambda, sigma)result(numResult)
        integer :: L, lambda, sigma(L,L)
        integer, dimension(1:3) :: tempSpins
        integer, dimension(L/lambda,L/lambda,1:3) :: numResult

        integer :: i, j, first(1:2), x1, x2

        do i = 1, L/lambda
            do j = 1, L/lambda
                first = findFirst(i, j)

                !Reset for each run.
                tempSpins = [0,0,0]

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

                do k = 1, 3
                    numResult(i,j,k) = tempSpins(k)
                enddo

            enddo
        enddo

        return
    end

    function calcEnergy(sigma)result(energyResult)
        integer :: sigma(L,L)
    
        integer :: numSpins(L/lambda,L/lambda,1:3)
        integer :: i_c, j_c, y_down, y_up, x_left, x_right, k, b, h
        real :: E_current
        integer, dimension(1:4) :: PBCTemp
        integer, dimension(1:8) :: cellNeighbours
        real :: energyResult(L/lambda,L/lambda)

        numSpins = initialNum(L, lambda, sigma) !Determines the initial number of spins.
        
        do i_c = 1, L/lambda
            do j_c = 1, L/lambda
                PBCtemp = PBC(i_c, j_c, L/lambda)
                y_up = PBCTemp(1); y_down = PBCTemp(2); x_left = PBCTemp(3); x_right = PBCTemp(4)
                
                cellNeighbours = [[i_c, x_right], [y_up, j_c], [i_c, x_left], [y_down, j_c]]
                
                !Reset for each cell.
                E_current = 0

                do h = 1, 7, 2
                    if(FBC .eqv. .true.) then
                        if(i_c == L/lambda) then
                            if(h == 7) then
                                go to 20
                            endif
                        elseif(i_c == 1) then
                            if(h == 3) then
                                go to 20
                            endif
                        endif
                    endif
                    do k = 1, 3
                        do b = 1, 3
                            E_current = E_current + numSpins(i_c, j_c, k)*&
                                numSpins(cellNeighbours(h), cellNeighbours(h+1), b)*J_str(k,b)
                        enddo
                    enddo
        20          continue
                enddo
                
                energyResult(i_c, j_c) = E_current
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
            do i = 1,L
                read(10,*) (sigma(i,j), j = 1,L)
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

    integer :: n, i, j
    character(256) :: file_id, file_name
    
    call findParamaters(dir, lambda, L)
    allocate(sigma(L,L))
    allocate(numSpins(L/lambda,L/lambda,1:3))
    allocate(energyCells(L/lambda,L/lambda))
    J_str = transpose(reshape(real(lambda)**(-2)*0.01* &
        [0, 75, 125, 75, 0, 75, 125, 75, 0], shape(J_str))) !1
!    J_str = transpose(reshape([0, 1, 6, 1, 0, 1, 6, 1, 0], shape(J_str)))

    call nFinder(dir, n, L)

    write(file_id, '(i0)') n
    file_name = trim(dir) // '/frame-' // trim(adjustl(file_id)) // '.dat'
    open(10, file = trim(file_name), form = 'formatted')
    do i = 1,L
        read(10,*) (sigma(i,j), j = 1,L)
    enddo
    close(10)

    energyCells = calcEnergy(sigma)
    print *, 'H_int = ', sum(energyCells)/L**2

end program interfacialEnergy
