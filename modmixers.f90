MODULE modmixers
    implicit none
    private
    logical :: MIXERS_SHOW_DEBUG = .true.

    type simplemixer

        doubleprecision, allocatable, dimension(:,:) :: vecIn,vecOut
        doubleprecision, allocatable, dimension(:)   :: mixedVec
        double precision :: omega
        integer          :: nx
        integer          :: M ! liczba pamietanych krokow
        integer          :: iter

        contains
        ! -------------------------------------------------
        ! Procedury klasy
        ! -------------------------------------------------
        procedure, public, pass(mixer) :: init !()
        procedure, public, pass(mixer) :: free_mixer !()
        procedure, public, pass(mixer) :: set_input_vec! (mixer,iVec)
        procedure, public, pass(mixer) :: set_output_vec! (mixer,oVec)
        procedure, public, pass(mixer) :: mix! (mixer)
        procedure, public, pass(mixer) :: get_last_residuum! (mixer)



    endtype simplemixer
    type, extends(simplemixer) :: andersonmixer



        contains
        ! -------------------------------------------------
        ! Procedury klasy
        ! -------------------------------------------------
!        procedure, public, pass(zrodlo) :: zrodlo_zwolnij_pamiec!() wiadomo

    endtype andersonmixer


    type, extends(simplemixer) :: broydenmixer



        contains
        ! -------------------------------------------------
        ! Procedury klasy
        ! -------------------------------------------------
!        procedure, public, pass(zrodlo) :: zrodlo_zwolnij_pamiec!() wiadomo

    endtype broydenmixer


    public :: simplemixer
    public :: andersonmixer
    public :: broydenmixer
    public :: MIXERS_SHOW_DEBUG

    contains


    subroutine init(mixer,nx,omega,M)
        class(simplemixer) :: mixer
        integer :: nx
        doubleprecision :: omega
        integer,optional :: M
        call mixer%free_mixer()

        mixer%M     = 1
        select type (mixer)
            type is (simplemixer)
            ! if anderson store M elements
            class is (andersonmixer)
                print*,"mixer:: initializing anderson mixer"
                mixer%M = 2
                if(present(M)) mixer%M     = M


        end select
        mixer%iter  = 0
        mixer%nx    = nx
        mixer%omega = omega
        allocate(mixer%vecIn (nx,mixer%M))
        allocate(mixer%vecOut(nx,mixer%M))
        allocate(mixer%mixedVec(nx))

        mixer%vecIn  = 0
        mixer%vecOut = 0
    endsubroutine init


    doubleprecision function get_last_residuum(mixer) result(rval)
        class(simplemixer) :: mixer
        rval = sum((mixer%vecIn(:,1)-mixer%vecOut(:,1))**2)/sum((mixer%vecIn(:,1)+mixer%vecOut(:,1))**2)
    endfunction get_last_residuum


    subroutine free_mixer(mixer)
        class(simplemixer) :: mixer
        if(allocated(mixer%vecIn ))  deallocate(mixer%vecIn)
        if(allocated(mixer%vecOut))  deallocate(mixer%vecOut)
        if(allocated(mixer%mixedVec))deallocate(mixer%mixedVec)
    endsubroutine free_mixer

    subroutine set_input_vec(mixer,inVec)
        class(simplemixer) :: mixer
        doubleprecision, dimension(:) :: inVec
        integer :: m

        select type (mixer)

            type is (simplemixer)
                if(MIXERS_SHOW_DEBUG) then
                print*,"mixer:: simple:: pushing input vectors"
                endif
                mixer%vecIn(:,1) = inVec
            ! if anderson store M elements
            class is (andersonmixer)
                if(MIXERS_SHOW_DEBUG) then
                print*,"mixer:: anderson:: pushing input vectors"
                endif
                do m = mixer%M , 2 , -1
                mixer%vecIn(:,m) = mixer%vecIn(:,m-1)
                enddo
                mixer%vecIn(:,1) = inVec
        end select

    endsubroutine set_input_vec

    subroutine set_output_vec(mixer,outVec)
        class(simplemixer) :: mixer
        doubleprecision, dimension(:) :: outVec
        integer :: m

        select type (mixer)


            type is (simplemixer)
                if(MIXERS_SHOW_DEBUG) then
                print*,"mixer:: simple:: pushing output vectors"
                endif
                mixer%vecOut(:,1) = outVec
            ! if anderson store M elements
            class is (andersonmixer)
                if(MIXERS_SHOW_DEBUG) then
                print*,"mixer:: anderson:: pushing output vectors"
                endif
                do m = mixer%M , 2 , -1
                mixer%vecOut(:,m) = mixer%vecIn(:,m-1)
                enddo
                mixer%vecOut(:,1) = outVec
        end select
    endsubroutine set_output_vec

    subroutine mix(mixer)
        class(simplemixer) :: mixer
        double precision ::  abeta , anorm , alpha , aux1 ,aux2
        doubleprecision :: beta(mixer%M) , deltaFm , deltaFm1 , F
        integer :: i , k , niter

        mixer%iter = mixer%iter + 1

        select type (mixer)
            type is (simplemixer)
                if(MIXERS_SHOW_DEBUG) then
                    print*,"mixer:: simple:: calculating output vec"
                endif
                mixer%mixedVec = (1-mixer%omega)*mixer%vecIn(:,1) + mixer%omega * mixer%vecOut(:,1)

            ! if anderson store M elements
            class is (andersonmixer)
                if(MIXERS_SHOW_DEBUG) then
                    print*,"mixer:: anderson:: calculating output vec"
                endif
                if(mixer%iter < 3) then
                    mixer%mixedVec = (1-mixer%omega)*mixer%vecIn(:,1) + &
                                        mixer%omega * mixer%vecOut(:,1)
                else if(mixer%M == 1) then
                    mixer%mixedVec = (1-mixer%omega)*mixer%vecIn(:,1) + &
                    mixer%omega * mixer%vecOut(:,1)
                else
                    mixer%mixedVec = 0
                    beta = 0
                    niter = min(mixer%M-1,mixer%iter)
                    do k = 1 , niter
                        abeta = sum( (mixer%vecOut(:,k) - mixer%vecIn(:,k))*( &
                                      mixer%vecOut(:,k) - mixer%vecIn(:,k) -  &
                                     (mixer%vecOut(:,k+1) - mixer%vecIn(:,k+1)) ) &
                                     )
                        anorm = sum(abs((mixer%vecOut(:,k)-mixer%vecIn(:,k)) -  &
                                       (mixer%vecOut(:,k+1)-mixer%vecIn(:,k+1)) )**2)

                        beta(k+1)  = - abeta/anorm * 0.05
                    enddo

                    beta(1) = 1 - sum(beta(2:mixer%M))
                    print"(A,100f16.4)","mixer:: anderson betas:",beta
                    alpha = mixer%omega
                    do i = 1 , mixer%nx
                        aux1 = 0
                        aux2 = 0
                        do k = 1 , mixer%M
                            aux1 = aux1 + ( beta(k) ) * mixer%vecIn(i,k)
                            aux2 = aux2 + ( beta(k) ) * mixer%vecOut(i,k)
                        enddo
                        mixer%mixedVec(i) = (1-alpha)*aux1 + alpha * aux2
                    enddo


!                    anorm = sum( (mixer%vecOut(:,1) - mixer%vecIn(:,1))*( &
!                                  mixer%vecOut(:,1) - mixer%vecIn(:,1) -  &
!                                 (mixer%vecOut(:,2) - mixer%vecIn(:,2)) ) &
!                                 )
!                    abeta = sum(abs((mixer%vecOut(:,1)-mixer%vecIn(:,1)) -  &
!                                   (mixer%vecOut(:,2)-mixer%vecIn(:,2)) )**2)
!                    beta(1)  = anorm/abeta * 0.1
!                    alpha = mixer%omega
!                    print*,"beta=",beta
!                    do i = 1 , mixer%nx
!                        aux1 = ( 1 - beta(1) ) * mixer%vecIn(i,1)  + beta(1) * mixer%vecIn ( i , 2)
!                        aux2 = ( 1 - beta(1) ) * mixer%vecOut(i,1) + beta(1) * mixer%vecOut( i , 2)
!                        mixer%mixedVec(i) = (1-alpha)*aux1 + alpha * aux2
!                    enddo

!                else
!                    print*,"mixer:: anderson for M > 2 is not implemented"
!                    stop
                endif

        end select


    end subroutine mix



end module modmixers
