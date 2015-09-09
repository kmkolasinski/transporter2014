MODULE moddevice

    use modjed
    use modutils
    use modmixers
    USE, INTRINSIC :: ISO_C_BINDING
!DEC$ IF DEFINED  (USE_UMF_PACK)
    use mUMFPACK
!DEC$ ENDIF
    implicit none
    private
!DEC$ IF DEFINED  (USE_PARDISO)
    EXTERNAL PARDISO_GETENV
    INTEGER  PARDISO_GETENV

    EXTERNAL PARDISO_SETENV
    INTEGER  PARDISO_SETENV

    INTEGER PARDISO_OOC_FILE_NAME
    PARAMETER ( PARDISO_OOC_FILE_NAME = 1 )
!DEC$ ENDIF
    integer :: nx,ny,nz
    doubleprecision :: dx
    integer,dimension(:,:,:),allocatable         :: DEVICE_FLAGS
    doubleprecision,dimension(:,:,:),allocatable :: DEVICE_DENS
    doubleprecision,dimension(:,:,:),allocatable :: DEVICE_POT


    integer,dimension(:,:,:), allocatable        :: GINDEX ! INDEKSUJE GEOMETRIE UKLADU (-1) OZNACZA POZA UKLADEM
    doubleprecision,dimension(:),allocatable     :: MATA   ! GLOWNA MACIERZ PROGRAMU W FORMACIE (ROW,COL,VALS), TUTAJ TYLKO VALS
    integer,dimension(:,:),allocatable           :: IDXA   ! INDEKSY MACIERZY (ROW,COL)
    integer,allocatable :: HBROWS(:)
    integer :: matAsize
    integer :: DEVICE_MESH_SIZE
    doubleprecision,dimension(:),allocatable     :: VECPOT

    !type(simplemixer) :: mixer
    type(andersonmixer) :: mixer

    type rectancular_material
        logical :: bActive
        doubleprecision :: xstart,ystart,zstart,xstop,ystop,zstop,rho0,bparam
    endtype rectancular_material

    integer :: no_tf_dens,no_cnst_dens,no_cnst_pot
    integer,parameter :: MAX_NO_TF_DENS = 100
    type(rectancular_material) tfdens(MAX_NO_TF_DENS)
    type(rectancular_material) cnstdens(MAX_NO_TF_DENS)
    type(rectancular_material) cnstpot(MAX_NO_TF_DENS)

    public :: device_create, device_free
    public :: device_solve,device_solve_poisson
    public :: device_add_const_dens,device_add_const_pot,device_add_tf_dens
    public :: DEVICE_FLAGS , DEVICE_DENS , DEVICE_POT
    contains


    subroutine device_create(pnx,pny,pnz)
        integer :: pnx,pny,pnz
        integer :: i,j,k
        print*,"Tworzenie ukladu:",pnx,pny,pnz

        nx = pnx
        ny = pny
        nz = pnz
        DEVICE_MESH_SIZE = nx*ny*nz
        call device_free()
        allocate(DEVICE_FLAGS(nx,ny,nz))
        allocate(DEVICE_DENS(nx,ny,nz))
        allocate(DEVICE_POT(nx,ny,nz))
        allocate(GINDEX    (nx,ny,nz))
        allocate(IDXA    (DEVICE_MESH_SIZE*7,2))
        allocate(HBROWS  (DEVICE_MESH_SIZE*7))
        allocate(MATA    (DEVICE_MESH_SIZE*7))
        allocate(VECPOT  (DEVICE_MESH_SIZE))
        DEVICE_POT   = 0
        DEVICE_DENS  = 0

        DEVICE_FLAGS = DFLAG_NEUMAN
        do i = 2 , nx-1
        do j = 2 , ny-1
        do k = 2 , nz-1
            DEVICE_FLAGS(i,j,k) = DFLAG_NORMAL
        enddo
        enddo
        enddo
        dx = atomic_DX * L2LR

        call mixer%init(DEVICE_MESH_SIZE,0.8D0,4)

        no_tf_dens    = 0
        no_cnst_dens  = 0
        no_cnst_pot   = 0
    end subroutine device_create

    subroutine device_free()
        call mixer%free_mixer()
        if(allocated(DEVICE_FLAGS)) deallocate(DEVICE_FLAGS)
        if(allocated(DEVICE_DENS))  deallocate(DEVICE_DENS)
        if(allocated(DEVICE_POT))   deallocate(DEVICE_POT)
        if(allocated(GINDEX))       deallocate(GINDEX)
        if(allocated(IDXA))         deallocate(IDXA)
        if(allocated(MATA))         deallocate(MATA)
        if(allocated(VECPOT))       deallocate(VECPOT)
        if(allocated(HBROWS))       deallocate(HBROWS)
    end subroutine device_free

    subroutine device_add_const_dens(xstart,xstop,ystart,ystop,zstart,zstop,dvalue)
        doubleprecision :: xstart,ystart,zstart,xstop,ystop,zstop,dvalue
        integer :: i,j,k
        doubleprecision :: volume , nvalue


        nvalue = dvalue * (1 / (L2LR**3))


        no_cnst_dens = no_cnst_dens + 1
        cnstdens(no_cnst_dens)%rho0   =  nvalue
        cnstdens(no_cnst_dens)%bparam =  0.0 ! not used here
        cnstdens(no_cnst_dens)%bActive= .true.
        cnstdens(no_cnst_dens)%xstart = xstart
        cnstdens(no_cnst_dens)%ystart = ystart
        cnstdens(no_cnst_dens)%zstart = zstart
        cnstdens(no_cnst_dens)%xstop  = xstop
        cnstdens(no_cnst_dens)%ystop  = ystop
        cnstdens(no_cnst_dens)%zstop  = zstop


        volume = 0


        do i =  xstart / atomic_DX , xstop / atomic_DX
            do j =  ystart / atomic_DX , ystop / atomic_DX
                do k =  zstart / atomic_DX , zstop / atomic_DX
                    ! zeby nie wyjsc poza uklad
                    if(i > 0 .and. i <= nx .and. j > 0 .and. j <= ny .and. k > 0 .and. k <= nz) then
                        DEVICE_DENS(i,j,k) = nvalue
                        volume = volume + dx**3
                    endif
                enddo
            enddo
        enddo
        print"(A,e12.2,A,f8.2)"," device:: Added constant density:",dvalue," with total charge:",nvalue*volume

    end subroutine device_add_const_dens

    subroutine device_add_tf_dens(xstart,xstop,ystart,ystop,zstart,zstop,dvalue)
        doubleprecision :: xstart,ystart,zstart,xstop,ystop,zstop,dvalue
        integer :: i,j,k
        doubleprecision :: volume , nvalue

        volume = 0
        nvalue = dvalue * (1 / (L2LR**3))
        no_tf_dens = no_tf_dens + 1

        tfdens(no_tf_dens)%rho0   =  nvalue
        tfdens(no_tf_dens)%bparam =  2.0 * M_EFF / (3.0 * M_PI) **(2.0/3.0)


        tfdens(no_tf_dens)%bActive= .true.

        tfdens(no_tf_dens)%xstart = xstart
        tfdens(no_tf_dens)%ystart = ystart
        tfdens(no_tf_dens)%zstart = zstart

        tfdens(no_tf_dens)%xstop  = xstop
        tfdens(no_tf_dens)%ystop  = ystop
        tfdens(no_tf_dens)%zstop  = zstop

        do i =  xstart / atomic_DX , xstop / atomic_DX
            do j =  ystart / atomic_DX , ystop / atomic_DX
                do k =  zstart / atomic_DX , zstop / atomic_DX
                    ! zeby nie wyjsc poza uklad
                    if(i > 0 .and. i <= nx .and. j > 0 .and. j <= ny .and. k > 0 .and. k <= nz) then
                    DEVICE_DENS(i,j,k) = -nvalue
                    volume = volume + dx**3
                    endif
                enddo
            enddo
        enddo
        print"(A,e12.2,A,f8.2)"," device:: Added Thomas-Fermi density:",dvalue," with total charge:",nvalue*volume
        print*,"        Udepl:",tfdens(no_tf_dens)%rho0**(2.0/3.0)/tfdens(no_tf_dens)%bparam,"[meV]"

    end subroutine device_add_tf_dens


    subroutine device_add_const_pot(xstart,xstop,ystart,ystop,zstart,zstop,dvalue)
        doubleprecision :: xstart,ystart,zstart,xstop,ystop,zstop,dvalue
        integer :: i,j,k

        no_cnst_pot = no_cnst_pot + 1
        cnstpot(no_cnst_pot)%rho0   =  dvalue/Rd
        cnstpot(no_cnst_pot)%bparam =  0.0 ! not used here
        cnstpot(no_cnst_pot)%bActive= .true.
        cnstpot(no_cnst_pot)%xstart = xstart
        cnstpot(no_cnst_pot)%ystart = ystart
        cnstpot(no_cnst_pot)%zstart = zstart
        cnstpot(no_cnst_pot)%xstop  = xstop
        cnstpot(no_cnst_pot)%ystop  = ystop
        cnstpot(no_cnst_pot)%zstop  = zstop

        do i =  xstart / atomic_DX , xstop / atomic_DX
            do j =  ystart / atomic_DX , ystop / atomic_DX
                do k =  zstart / atomic_DX , zstop / atomic_DX
                    ! zeby nie wyjsc poza uklad
                    if(i > 0 .and. i <= nx .and. j > 0 .and. j <= ny .and. k > 0 .and. k <= nz) then
                        DEVICE_POT  (i,j,k) = dvalue/Rd
                        DEVICE_FLAGS(i,j,k) = DFLAG_DIRICHLET
                    endif
                enddo
            enddo
        enddo
        print"(A,e12.2,A,f8.2)"," device:: Added constant potential:",dvalue

    end subroutine device_add_const_pot


    subroutine device_set_geometry(grid_scale)
        integer :: grid_scale
        integer :: i,j,k,d,g
        doubleprecision :: xstart,ystart,zstart,xstop,ystop,zstop,aux


        DEVICE_DENS = 0
        do i =  2 , nx-1
            do j =  2 , ny-1
                do k =  2 , nz -1
                     DEVICE_FLAGS(i,j,k) = DFLAG_NORMAL
                enddo
            enddo
        enddo


        g = grid_scale
        do d = 1 , no_cnst_dens

            xstart = cnstdens(d)%xstart
            ystart = cnstdens(d)%ystart
            zstart = cnstdens(d)%zstart

            xstop = cnstdens(d)%xstop
            ystop = cnstdens(d)%ystop
            zstop = cnstdens(d)%zstop
            aux   = cnstdens(d)%rho0


            do i =  xstart / atomic_DX / g , xstop / atomic_DX / g
                do j =  ystart / atomic_DX / g , ystop / atomic_DX / g
                    do k =  zstart / atomic_DX / g , zstop / atomic_DX / g
                        ! zeby nie wyjsc poza uklad
                        if(i > 0 .and. i <= nx/g .and. j > 0 .and. j <= ny/g .and. k > 0 .and. k <= nz/g) then
                            DEVICE_DENS(i*g,j*g,k*g) = aux
                        endif
                    enddo
                enddo
            enddo
        enddo ! end of const dens

        do d = 1 , no_cnst_pot

            xstart = cnstpot(d)%xstart
            ystart = cnstpot(d)%ystart
            zstart = cnstpot(d)%zstart

            xstop = cnstpot(d)%xstop
            ystop = cnstpot(d)%ystop
            zstop = cnstpot(d)%zstop
            aux   = cnstpot(d)%rho0


            do i =  xstart / atomic_DX / g , xstop / atomic_DX / g
                do j =  ystart / atomic_DX / g , ystop / atomic_DX / g
                    do k =  zstart / atomic_DX / g , zstop / atomic_DX / g
                        ! zeby nie wyjsc poza uklad
                        if(i > 0 .and. i <= nx/g .and. j > 0 .and. j <= ny/g .and. k > 0 .and. k <= nz/g) then
                            DEVICE_POT(i*g,j*g,k*g)   = aux
                            DEVICE_FLAGS(i*g,j*g,k*g) = DFLAG_DIRICHLET
                        endif
                    enddo
                enddo
            enddo
        enddo ! end of const dens


    end subroutine device_set_geometry


    subroutine device_solve()

        integer :: i,j,k,iter,d

        doubleprecision :: xstart,ystart,zstart,xstop,ystop,zstop,dvalue,aux,U
        MIXERS_SHOW_DEBUG = .true.

        !call device_solve_poisson_directly(DIRECT_POISSON_PREP_MAT)
        do iter = 0 , 500
        print*,"Iteracja:",iter," Last residuum:",mixer%get_last_residuum()
        call mixer%set_input_vec(pack(DEVICE_DENS, .true. ))

        call device_set_geometry(2)
        call modutils_3darray2VTK(DEVICE_DENS/cnstdens(1)%rho0,atomic_DX,"rho")
        call modutils_3darray2VTK(DEVICE_POT,dx,"pot")
        !stop

        do d = 1 , no_tf_dens

            xstart = tfdens(d)%xstart
            ystart = tfdens(d)%ystart
            zstart = tfdens(d)%zstart

            xstop = tfdens(d)%xstop
            ystop = tfdens(d)%ystop
            zstop = tfdens(d)%zstop
            aux =   -abs(tfdens(d)%rho0)**(2.0/3.0)/ tfdens(d)%bparam
            do i =  xstart / atomic_DX , xstop / atomic_DX
                do j =  ystart / atomic_DX , ystop / atomic_DX
                    do k =  zstart / atomic_DX , zstop / atomic_DX
                        ! zeby nie wyjsc poza uklad
                        if(i > 0 .and. i <= nx .and. j > 0 .and. j <= ny .and. k > 0 .and. k <= nz) then
                            U = -DEVICE_POT(i,j,k) ! energia potencjalna
                            if( U > aux ) then
                                DEVICE_DENS(i,j,k) = tfdens(d)%rho0 &
                                - ( tfdens(d)%rho0**(2.0/3.0) + tfdens(d)%bparam * U )**(3.0/2.0)

                                !DEVICE_DENS(i,j,k) = DEVICE_DENS(i,j,k) / ( exp( (U - aux)/0.001 )+1)
                            else
                                DEVICE_DENS(i,j,k) =  tfdens(d)%rho0
                            endif

                            DEVICE_DENS(i,j,k) = tfdens(d)%rho0 / ( exp( (U - aux)/0.01 )+1)
                            DEVICE_DENS(i,j,k) = -DEVICE_DENS(i,j,k)
                            !print*,i,j,k,DEVICE_DENS(i,j,k),"r0=",tfdens(d)%rho0
                        endif
                    enddo
                enddo
            enddo
        enddo ! end of

        call mixer%set_output_vec(pack(DEVICE_DENS, .true. ))
        call mixer%mix()
        DEVICE_DENS = reshape(mixer%mixedVec,(/nx,ny,nz/))
        write(54232,*),iter ,mixer%get_last_residuum()

        !if(mixer%get_last_residuum() < 1.0D-3) then
        call device_solve_poisson()
        !else
        !    call device_solve_poisson_directly(DIRECT_POISSON_SOLVE)
        !endif
        if(mixer%get_last_residuum() < 1.0D-10 ) exit


        !if(mod(iter,5) == 0) then
            call write_to_file(321,"rho.txt",DEVICE_DENS(:,:,nz-15)/tfdens(1)%rho0,nx,ny)
            call write_to_file(321,"pot.txt",DEVICE_POT(:,:,nz-15),nx,ny)
        !endif

        enddo
        !call device_solve_poisson_directly(DIRECT_POISSON_FREE_SOLVER)

    end subroutine device_solve


    subroutine device_solve_poisson()

        integer :: i,j,k,iter

        doubleprecision :: delta2pi
        doubleprecision :: dx2
        doubleprecision :: aux1,aux2,delta1,delta2,eps,omega


        delta2pi = 4*dx*dx*M_PI
        dx2      = dx*dx

        delta1 = 5.0
        delta2 = 1.0
        eps    = 0.001
        iter = 0

        print*,"device:: solver start"
        do while( abs(delta2-delta1) > eps .or. iter < 10)
            delta2 = delta1

            call device_iterate_solution(delta1,2)

            iter = iter + 1
            if(mod(iter,500) == 0) then
                print*, iter , "eps:" ,  abs(delta2-delta1)
            endif
        enddo ! end of while

    end subroutine device_solve_poisson


    subroutine device_solve_poisson_grid(grid_scale)
        integer :: grid_scale
        integer :: i,j,k,iter

        doubleprecision :: delta1,delta2,eps


        delta1 = 5.0
        delta2 = 1.0
        eps    = 0.001
        iter = 0

        print*,"device:: solver start with grid:",grid_scale



        do while( abs(delta2-delta1) > eps .or. iter < 10)
            delta2 = delta1
            call device_iterate_solution(delta1,grid_scale)
            iter = iter + 1
            if(mod(iter,500) == 0) then
                print*, iter , "eps:" ,  abs(delta2-delta1)
            endif
        enddo ! end of while

    end subroutine device_solve_poisson_grid



    subroutine device_iterate_solution(delta1,grid_scale)
        doubleprecision ,intent(inout) :: delta1
        integer :: grid_scale
        integer :: i,j,k,g

        doubleprecision :: delta2pi
        doubleprecision :: dx2
        doubleprecision :: aux1,aux2,omega


        omega    = 1.8
        delta2pi = 4*(dx*grid_scale)*(dx*grid_scale)*M_PI/E_MAT
        dx2      = (dx*grid_scale)*(dx*grid_scale)

        delta1  = 0
        g = grid_scale
        do i = grid_scale , nx - grid_scale+1 , grid_scale
        do j = grid_scale , ny - grid_scale+1 , grid_scale
        do k = grid_scale , nz - grid_scale+1 , grid_scale
            if(DEVICE_FLAGS(i,j,k) == DFLAG_NORMAL) then
                aux1 = (DEVICE_POT(i+g,j,k)+DEVICE_POT(i-g,j,k) + &
                        DEVICE_POT(i,j+g,k)+DEVICE_POT(i,j-g,k) + &
                        DEVICE_POT(i,j,k+g)+DEVICE_POT(i,j,k-g))
                DEVICE_POT(i,j,k) = (1-omega)*DEVICE_POT(i,j,k)  &
                            + (omega)*(aux1 + delta2pi * DEVICE_DENS(i,j,k))/6.0
            else if(DEVICE_FLAGS(i,j,k) == DFLAG_NEUMAN) then
                aux1 = 0
                aux2 = 0
                if(i == 1 ) then      ; aux1 = aux1 + DEVICE_POT(i+g,j,k) ; aux2 = aux2 + 1;
                else if(i == nx) then ; aux1 = aux1 + DEVICE_POT(i-g,j,k) ; aux2 = aux2 + 1;
                else if(j == 1 ) then ; aux1 = aux1 + DEVICE_POT(i,j+g,k) ; aux2 = aux2 + 1;
                else if(j == ny) then ; aux1 = aux1 + DEVICE_POT(i,j-g,k) ; aux2 = aux2 + 1;
                else if(k == 1 ) then ; aux1 = aux1 + DEVICE_POT(i,j,k+g) ; aux2 = aux2 + 1;
                else if(k == nz) then ; aux1 = aux1 + DEVICE_POT(i,j,k-g) ; aux2 = aux2 + 1;
                endif
                DEVICE_POT(i,j,k) = aux1
            else if(DEVICE_FLAGS(i,j,k) == DFLAG_DIRICHLET) then
                print*,"dirichlet:",DEVICE_POT(i,j,k)
            endif
            delta1 = delta1 +  DEVICE_POT(i,j,k)**2
        enddo
        enddo
        enddo
    end subroutine device_iterate_solution

    subroutine device_solve_poisson_directly(opt)
        integer :: opt
        integer :: i,j,k,iter,itmp
        doubleprecision :: dx2

        dx2      = dx*dx
        if(opt == DIRECT_POISSON_FREE_SOLVER) then
            call solve_system(DEVICE_MESH_SIZE,MATASIZE,IDXA(:,2),HBROWS,MATA(:),VECPOT,3)
            return
        endif
        if(opt == DIRECT_POISSON_PREP_MAT) then
        GINDEX = 0
        iter = 0
        do i = 1 , nx
        do j = 1 , ny
        do k = 1 , nz
            if(DEVICE_FLAGS(i,j,k) >= DFLAG_NORMAL) then
                iter = iter + 1
                GINDEX(i,j,k) = iter
            endif
        enddo
        enddo
        enddo

        VECPOT =  0
        itmp   = 0
        do i = 1 , nx
        do j = 1 , ny
        do k = 1 , nz
            if(DEVICE_FLAGS(i,j,k) == DFLAG_NORMAL) then

                itmp = itmp + 1
                matA(itmp)   = (1.0) / dx2
                idxA(itmp,1) = GINDEX(i  , j  , k)
                idxA(itmp,2) = GINDEX(i  , j  , k-1)

                itmp = itmp + 1
                matA(itmp)   = (1.0) / dx2
                idxA(itmp,1) = GINDEX(i  , j  , k)
                idxA(itmp,2) = GINDEX(i  , j-1, k)

                itmp = itmp + 1
                matA(itmp)   = (1.0) / dx2
                idxA(itmp,1) = GINDEX(i  , j, k)
                idxA(itmp,2) = GINDEX(i-1, j, k)

                itmp = itmp + 1
                matA(itmp)   = -(6.0) / dx2
                idxA(itmp,1) = GINDEX(i  , j  , k)
                idxA(itmp,2) = GINDEX(i  , j, k)

                itmp = itmp + 1
                matA(itmp)   = (1.0) / dx2
                idxA(itmp,1) = GINDEX(i  , j, k)
                idxA(itmp,2) = GINDEX(i+1, j, k)

                itmp = itmp + 1
                matA(itmp)   = (1.0) / dx2
                idxA(itmp,1) = GINDEX(i  , j  , k)
                idxA(itmp,2) = GINDEX(i  , j+1, k)


                itmp = itmp + 1
                matA(itmp)   = (1.0) / dx2
                idxA(itmp,1) = GINDEX(i  , j  , k)
                idxA(itmp,2) = GINDEX(i  , j  , k+1)


            else if(DEVICE_FLAGS(i,j,k) == DFLAG_NEUMAN) then
                if(i == 1 ) then ;
                    itmp = itmp + 1
                    matA(itmp)   = 1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i  , j, k)

                    itmp = itmp + 1
                    matA(itmp)   = -1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i+1, j, k)

                else if(i == nx) then ;

                    itmp = itmp + 1
                    matA(itmp)   = 1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i  , j, k)

                    itmp = itmp + 1
                    matA(itmp)   = -1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i-1, j, k)

                else if(j == 1 ) then
                    itmp = itmp + 1
                    matA(itmp)   = 1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i  , j, k)

                    itmp = itmp + 1
                    matA(itmp)   = -1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i  , j+1, k)
                else if(j == ny) then ;
                    itmp = itmp + 1
                    matA(itmp)   = 1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i  , j, k)

                    itmp = itmp + 1
                    matA(itmp)   = -1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i, j-1, k)
                else if(k == 1 ) then ;
                    itmp = itmp + 1
                    matA(itmp)   = 1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i  , j, k)

                    itmp = itmp + 1
                    matA(itmp)   = -1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i  , j, k+1)
                else if(k == nz) then ;
                    itmp = itmp + 1
                    matA(itmp)   = 1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i  , j, k)

                    itmp = itmp + 1
                    matA(itmp)   = -1.0
                    idxA(itmp,1) = GINDEX(i  , j, k)
                    idxA(itmp,2) = GINDEX(i  , j, k-1)
                endif
            else if(DEVICE_FLAGS(i,j,k) == DFLAG_DIRICHLET) then

                    itmp = itmp + 1
                    matA(itmp)   = (1.0)
                    idxA(itmp,1) = GINDEX(i, j, k)
                    idxA(itmp,2) = GINDEX(i, j, k)

            endif
        enddo
        enddo
        enddo
        matAsize = itmp
        print*,"Liczba niezerowych elementow macierzy:",matAsize," szacowana:",DEVICE_MESH_SIZE*7
        call convert_to_HB(MATASIZE,IDXA,HBROWS)
        call solve_system(DEVICE_MESH_SIZE,MATASIZE,IDXA(:,2),HBROWS,MATA(:),VECPOT,1)
        endif

        ! przeliczanie gestosci
        VECPOT = 0
        do i = 1 , nx
        do j = 1 , ny
        do k = 1 , nz
            if(DEVICE_FLAGS(i,j,k) == DFLAG_NORMAL) then
                VECPOT(GINDEX(i,j,k)) = - 4 * M_PI * DEVICE_DENS(i,j,k) / E_MAT
            else if(DEVICE_FLAGS(i,j,k) == DFLAG_DIRICHLET) then
                VECPOT(GINDEX(i,j,k)) = DEVICE_POT(i,j,k)
            endif
        enddo
        enddo
        enddo

        print*,"device:: direct poisson solve"
        call solve_system(DEVICE_MESH_SIZE,MATASIZE,IDXA(:,2),HBROWS,MATA(:),VECPOT,2)

        do i = 1 , nx
        do j = 1 , ny
        do k = 1 , nz
            DEVICE_POT(i,j,k) = VECPOT(GINDEX(i,j,k))
        enddo
        enddo
        enddo

    end subroutine device_solve_poisson_directly







    ! ==========================================================================
    !
    !
    !                          SUPER LU
    !
    !
    !
    ! ==========================================================================

    subroutine convert_to_HB(no_vals,rows_cols,out_rows)
          integer,intent(in)                  :: no_vals
          integer,intent(inout),dimension(:,:)  :: rows_cols
          integer,intent(inout),dimension(:) :: out_rows
          integer :: iterator, irow ,  from , to
          integer :: i, n

          n        = no_vals
          iterator = 0
          irow     = 0
          do i = 1 , n
              if( rows_cols(i,1) /= irow ) then
                iterator = iterator + 1
                out_rows(iterator) = i
                irow = rows_cols(i,1)
              endif
          enddo
          out_rows(iterator+1) = n + 1


!DEC$ IF DEFINED  (USE_UMF_PACK)
      irow = size(out_rows)-1
        ! sortowanie  kolumn
      do i = 1 , irow-1
      from = out_rows(i)
      to   = out_rows(i+1)-1
          call sort_col_vals(IDXA(from:to,2),cmatA(from:to))
      enddo

      ! przesuwanie indeksow do zera
      out_rows       = out_rows -1
      rows_cols(:,2) = rows_cols(:,2) -1
!DEC$ ENDIF

!DEC$ IF DEFINED  (USE_PARDISO)

      irow = size(out_rows)-1
        ! sortowanie  kolumn
      do i = 1 , irow-1
      from = out_rows(i)
      to   = out_rows(i+1)-1
          call sort_col_vals(IDXA(from:to,2),matA(from:to))
      enddo

!DEC$ ENDIF



      end subroutine convert_to_HB

    subroutine sort_col_vals(cols,vals)
            integer,intent(inout),dimension(:)    :: cols
            doubleprecision,intent(inout),dimension(:) :: vals
            integer :: tmp_col
            doubleprecision :: tmp_val
            integer :: i  , j , n
            logical :: test
            n = size(cols)

            test = .true.

            ! sortowanie bombelkowe
            do while(test)
              test = .false.
              do i = 1 , n-1
                if( cols(i) > cols(i+1)  ) then
                tmp_col   = cols(i)
                cols(i)   = cols(i+1)
                cols(i+1) = tmp_col

                tmp_val   = vals(i)
                vals(i)   = vals(i+1)
                vals(i+1) = tmp_val

                test = .true.
                exit
                endif
              enddo
            enddo
    end subroutine sort_col_vals



    subroutine solve_system(no_rows,no_vals,colptr,rowind,values,b,iopt)
        integer,intent(in)                 :: no_rows
        integer,intent(in)                 :: no_vals
        integer,intent(in),dimension(:)    :: colptr,rowind
        doubleprecision,intent(in),dimension(:) :: values
        doubleprecision,intent(inout),dimension(:) :: b
        integer :: iopt
        integer n, nnz, nrhs, ldb

        integer, save    ::  info = 0
        integer*8 , save :: factors = 0

        doubleprecision,save :: total_time
!DEC$ IF DEFINED  (USE_UMF_PACK)
        ! UMFPACK constants
        type(c_ptr),save :: symbolic,numeric
        ! zero-based arrays
        real(8),save :: control(0:UMFPACK_CONTROL-1),umf_info(0:UMFPACK_INFO-1)
        doubleprecision,allocatable,dimension(:),save :: b_sol

!DEC$ ENDIF

!DEC$ IF DEFINED  (USE_PARDISO)

        INTEGER*8,save  :: pt(64)
        INTEGER,save    :: phase
        INTEGER,save    :: maxfct, mnum, mtype, error, msglvl
        INTEGER,save    :: iparm(64)
        doubleprecision,allocatable,dimension(:),save :: b_sol

        INTEGER    ,save::  idum(1)
        doubleprecision ,save::  ddum(1)

!DEC$ ENDIF

        n    = no_rows
        nnz  = no_vals
        ldb  = n
        nrhs = 1


!DEC$ IF DEFINED  (USE_UMF_PACK)
      selectcase (iopt)
      case (1)
            total_time = get_clock();
            allocate(b_sol(size(b)))

            call umf4cdef (control)
            call umf4csym (n,n, rowind, colptr, values, symbolic, control, umf_info)
            call umf4cnum (rowind, colptr, values, symbolic, numeric, control, umf_info)
            call umf4cfsym (symbolic)
            !total_time =  umf_info(UMFPACK_NUMERIC_TIME)+umf_info(UMFPACK_SYMBOLIC_TIME)
            if (umf_info(UMFPACK_STATUS) .eq. 0) then
                if(TRANS_DEBUG) then
                     write (*,*) 'Factorization succeeded. Mem needed:', umf_info(UMFPACK_PEAK_MEMORY)/8.0/1024/1024 , "[MB]"
                endif
            else
                 write(*,*) 'UMFERROR: INFO from factorization = ', umf_info(UMFPACK_STATUS)
            endif

      case(2)
            b_sol = 0
            call umf4csolr (UMFPACK_Aat, rowind, colptr, values, b_sol, b, numeric, control, umf_info)
            b  = b_sol;

            if (umf_info(UMFPACK_STATUS) .eq. 0) then
                if(TRANS_DEBUG) then
                    write (*,*) 'Solve succeeded. Time needed:',umf_info(UMFPACK_SOLVE_WALLTIME)
                endif
            else
                 write(*,*) 'UMF ERROR: INFO from solve = ', umf_info(UMFPACK_STATUS)
            endif

      case(3)
            print*,"UMFPACK Solved:"
            print*,"Solve time needed:",get_clock()-total_time,"[s]"
            call umf4cfnum (numeric)
            deallocate(b_sol)
      endselect

!DEC$ ELSE IF DEFINED  (USE_PARDISO)


 selectcase (iopt)
      case (1)
      allocate(b_sol(size(b)))
          total_time = get_clock();
          maxfct = 1 ! in many application this is 1
          mnum   = 1 ! same here
          iparm = 0
          iparm(1) = 1 ! no solver default
          iparm(2) = 2 ! fill-in reordering from METIS
          iparm(3) = 1 ! numbers of processors, value of OMP_NUM_THREADS
          iparm(4) = 0 ! 0 - no iterative-direct algorithm, if 1 multirecursive iterative algorithm 61, 31 para me
          iparm(5) = 0 ! no user fill-in reducing permutation
          iparm(6) = 0 ! =0 solution on the first n compoments of x
          iparm(7) = 0 ! not in use
          iparm(8) = 2 ! numbers of iterative refinement steps
          iparm(9) = 0 ! not in use
          iparm(10) = 10 ! perturbe the pivot elements with 1E-13
          iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
          iparm(12) = 0 ! not in use
          iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric).
          iparm(14) = 0 ! Output: number of perturbed pivots
          iparm(15) = 0 ! not in use
          iparm(16) = 0 ! not in use
          iparm(17) = 0 ! not in use
          iparm(18) = -1 ! Output: number of nonzeros in the factor LU
          iparm(19) = -1 ! Output: Mflops for LU factorization
          iparm(20) = 0 ! Output: Numbers of CG Iterations
          iparm(32) = 0 ! if 1 use multirecursive iterative algorithm
           error = 0 ! initialize error flag
          msglvl = 0 ! print statistical information
          mtype     = 11      ! real unsymmetric matrix

          phase     = 11      ! only reordering and symbolic factorization
          CALL pardiso (pt, maxfct, mnum, mtype, phase, n, values, rowind, colptr,&
                       idum, nrhs, iparm, msglvl, ddum, ddum, error)

          !WRITE(*,*) 'Reordering completed ... '

          IF (error .NE. 0) THEN
            WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
          END IF


          !WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

    !C.. Factorization.
          phase     = 22  ! only factorization
          CALL pardiso (pt, maxfct, mnum, mtype, phase, n, values, rowind, colptr,&
                       idum, nrhs, iparm, msglvl, ddum, ddum, error)

          !WRITE(*,*) 'Factorization completed ...  '
          IF (error .NE. 0) THEN
             WRITE(*,*) 'The following ERROR was detected: ', error
            STOP 1
          ENDIF

          WRITE(*,*) 'Peak memory usage   = ',max (IPARM(15), IPARM(16)+IPARM(17))/1024.0,"[MB]"

      case(2)
          b_sol = 0
    !C.. Back substitution and iterative refinement
          phase     = 33  ! only factorization
          iparm(8)  = 3   ! max numbers of iterative refinement steps
          CALL pardiso (pt, maxfct, mnum, mtype, phase, n, values, rowind, colptr,&
                       idum, nrhs, iparm, msglvl, b, b_sol, error)

          b  = b_sol;
          !WRITE(*,*) 'Solve completed ... '
          IF (error .NE. 0) THEN
             WRITE(*,*) 'The following ERROR was detected: ', error

          ENDIF

      case(3)
    !C.. Termination and release of memory
            phase     = -1           ! release internal memory
            CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum,&
                       idum, nrhs, iparm, msglvl, ddum, ddum, error)

            print*,"PARDISO Solved:"
            print*,"Solve time needed:",get_clock()-total_time,"[s]"
            deallocate(b_sol)
      endselect


!DEC$ ELSE
      selectcase (iopt)
      case (1)
      total_time = get_clock();
! First, factorize the matrix. The factors are stored in *factors* handle.
      !iopt = 1
      call c_fortran_zgssv( iopt, n, nnz, nrhs, values, colptr , rowind , b, ldb,factors, info )
!
      if (info .eq. 0) then
         write (*,*) 'Factorization succeeded'
      else
         write(*,*) 'INFO from factorization = ', info
      endif
      case(2)
! Second, solve the system using the existing factors.
!      iopt = 2
      call c_fortran_zgssv( iopt, n, nnz, nrhs, values, colptr,rowind ,  b, ldb,factors, info )
!
      if (info .eq. 0) then
!         write (*,*) 'Solve succeeded'
!         write (*,*) (b(i), i=1, n)
      else
         write(*,*) 'INFO from triangular solve = ', info
      endif
      case(3)
! Last, free the storage allocated inside SuperLU
!      iopt = 3
      call c_fortran_zgssv( iopt, n, nnz, nrhs, values, colptr,rowind, b, ldb,factors, info )
      print*,"SuperLU Solved:"
      print*,"Solve time needed:",get_clock()-total_time,"[s]"
      endselect
!DEC$ ENDIF

      endsubroutine solve_system


end MODULE  moddevice
