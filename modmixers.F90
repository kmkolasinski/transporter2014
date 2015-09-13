MODULE modmixers
    implicit none
    private
    logical :: SCF_MIXER_SHOW_DEBUG = .false.

    ENUM,BIND(C)
        ENUMERATOR :: SCF_MIXER_SIMPLE            = 1
        ENUMERATOR :: SCF_MIXER_LOCAL_SIMPLE      = 2
        ENUMERATOR :: SCF_MIXER_ANDERSON          = 3
        ENUMERATOR :: SCF_MIXER_EXTENDED_ANDERSON = 4
        ENUMERATOR :: SCF_MIXER_BROYDEN           = 5
    END ENUM
    PUBLIC :: SCF_MIXER_SIMPLE , SCF_MIXER_LOCAL_SIMPLE , SCF_MIXER_ANDERSON , SCF_MIXER_EXTENDED_ANDERSON , SCF_MIXER_BROYDEN


    type scfmixer

        ! GENERAL VARIABLES USED BY ALL MIXERS
        doubleprecision, allocatable, dimension(:)   :: mixedVec ! mixed vector

        doubleprecision, allocatable, dimension(:,:) :: vecIn    ! SCF black box input vector  X
        doubleprecision, allocatable, dimension(:,:) :: vecOut   ! SCF black boc output vector F(X)
        doubleprecision, allocatable, dimension(:,:) :: vecRes   ! residual vector

        double precision :: omega             ! mixing parameter
        integer          :: NO_SITES          ! the lenght of the vector
        integer          :: NO_MEM_STEPS      ! number of backward interation steps
        integer          :: CURRENT_ITERATION ! iterator
        integer          :: MIXER_TYPE        ! same as variable name


        ! SIMPLE LOCAL MIXER VARIABLES
        ! vector of size NO_SITES containg the local mixing parameters
        doubleprecision, allocatable, dimension(:)   :: vecLocalOmega
        ! maximum value of "omega" in vector vecLocalOmega
        doubleprecision                              :: max_omega

        ! EXTENDED ANDERSON MIXER VARIABLES
        ! non linear mixing parameter
        doubleprecision :: omega0

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
        procedure, public, pass(mixer) :: get_current_residuum! (mixer)

    endtype scfmixer


    public :: scfmixer
    public :: SCF_MIXER_SHOW_DEBUG

    contains


    character(32) function mixer_get_name(mixer_type)
        integer :: mixer_type
        mixer_get_name = "SCF_MIXER_NONE"
        select case(mixer_type)
            case(SCF_MIXER_SIMPLE)
                mixer_get_name = "SCF_MIXER_SIMPLE"
            case(SCF_MIXER_LOCAL_SIMPLE)
                mixer_get_name = "SCF_MIXER_LOCAL_SIMPLE"
            case(SCF_MIXER_ANDERSON)
                mixer_get_name = "SCF_MIXER_ANDERSON"
            case(SCF_MIXER_EXTENDED_ANDERSON)
                mixer_get_name = "SCF_MIXER_EXTENDED_ANDERSON"
            case(SCF_MIXER_BROYDEN)
                mixer_get_name = "SCF_MIXER_BROYDEN"

        end select
    end function mixer_get_name


    ! ------------------------------------------------------------------------- !
    ! init - initializes all necesarry variables
    ! ------------------------------------------------------------------------- !
    subroutine init(mixer,mixer_type,nx,omega,max_omega,M,omega0)
        class(scfmixer) :: mixer
        integer :: mixer_type,nx
        doubleprecision          :: omega
        doubleprecision,optional :: max_omega
        integer,optional         :: M
        doubleprecision,optional :: omega0



        mixer%MIXER_TYPE         = mixer_type
        mixer%NO_MEM_STEPS       = 1 ! default value can be change according to method
        mixer%CURRENT_ITERATION  = 0
        mixer%NO_SITES           = nx
        mixer%omega              = omega
        mixer%max_omega          = omega
        mixer%omega0             = 0.01 ! smal value by default

        call mixer%free_mixer() ! just in case free memory

        select case(mixer%MIXER_TYPE)
            case(SCF_MIXER_SIMPLE)

                ! nothing needs to be done here ...

            case(SCF_MIXER_LOCAL_SIMPLE)
                mixer%NO_MEM_STEPS  = 2 ! forcing number of memory steps
                if(present(max_omega)) mixer%max_omega  = max_omega

!                allocate(mixer%vecLocalOmega(nx))
!                mixer%vecLocalOmega = omega


            case(SCF_MIXER_ANDERSON)

                mixer%NO_MEM_STEPS                 = 2
                if(present(M)) mixer%NO_MEM_STEPS  = M


            case(SCF_MIXER_EXTENDED_ANDERSON)

                mixer%NO_MEM_STEPS                 = 5
                if(present(M)) mixer%NO_MEM_STEPS  = M
                if(present(omega0)) mixer%omega0   = omega0

            case(SCF_MIXER_BROYDEN)

                mixer%NO_MEM_STEPS                 = 5
                if(present(M)) mixer%NO_MEM_STEPS  = M
                if(present(omega0)) mixer%omega0   = omega0
        end select

        !if(SCF_MIXER_SHOW_DEBUG) then
            print*,"mixer:: initializing ",mixer_get_name(mixer%MIXER_TYPE)
            print*,"    Memory size:", mixer%NO_MEM_STEPS
            print*,"    Vector size:", mixer%NO_SITES
            print*,"    Omega min  :", mixer%omega

            if(mixer%MIXER_TYPE == SCF_MIXER_LOCAL_SIMPLE) then
                print*,"    Omega max  :", mixer%max_omega
            endif
            if(mixer%MIXER_TYPE == SCF_MIXER_EXTENDED_ANDERSON) then
                print*,"    Omega0     :", mixer%omega0
            endif
            if(mixer%MIXER_TYPE == SCF_MIXER_BROYDEN) then
                print*,"    Omega0     :", mixer%omega0
            endif
        !endif


        allocate(mixer%vecIn (nx,mixer%NO_MEM_STEPS))
        allocate(mixer%vecOut(nx,mixer%NO_MEM_STEPS))
        allocate(mixer%vecRes(nx,mixer%NO_MEM_STEPS))
        allocate(mixer%mixedVec(nx))
        allocate(mixer%vecLocalOmega(nx))
        mixer%vecLocalOmega = omega
        mixer%vecIn         = 0
        mixer%vecOut        = 0
        mixer%vecRes        = 0
        mixer%mixedVec      = 0


    endsubroutine init


    ! ------------------------------------------------------------------------- !
    ! get_last_residuum - returns the last normalized value of the residual vector
    !   Can be used for the simulation stopping criterium.
    ! ------------------------------------------------------------------------- !
    doubleprecision function get_last_residuum(mixer) result(rval)
        class(scfmixer) :: mixer
        rval = sum((mixer%vecRes(:,1))**2)/sum((mixer%vecRes(:,1))**2)
    endfunction get_last_residuum

    ! ------------------------------------------------------------------------- !
    ! get_current_residuum - returns the last normalized value of the residual vector
    !   Can be used for the simulation stopping criterium.
    ! ------------------------------------------------------------------------- !
    doubleprecision function get_current_residuum(mixer) result(rval)
        class(scfmixer) :: mixer
        rval = sum((mixer%mixedVec(:)-mixer%vecIn(:,1))**2)/sum((mixer%mixedVec(:)+mixer%vecIn(:,1))**2)
    endfunction get_current_residuum

    ! ------------------------------------------------------------------------- !
    ! free_mixer - frees memory. Use it wisely!
    ! ------------------------------------------------------------------------- !
    subroutine free_mixer(mixer)
        class(scfmixer) :: mixer

        if(SCF_MIXER_SHOW_DEBUG) then
            print*,"mixer:: cleaning memory ",mixer_get_name(mixer%MIXER_TYPE)
        endif
        if(allocated(mixer%vecIn ))         deallocate(mixer%vecIn)
        if(allocated(mixer%vecOut))         deallocate(mixer%vecOut)
        if(allocated(mixer%mixedVec))       deallocate(mixer%mixedVec)
        if(allocated(mixer%vecRes))         deallocate(mixer%vecRes)
        if(allocated(mixer%vecLocalOmega))  deallocate(mixer%vecLocalOmega)
    endsubroutine free_mixer

    ! ------------------------------------------------------------------------- !
    ! set_input_vec - set the current vector which is used for SCF calculation
    !   This the input vector for the SCF black box. This can be an electrostatic
    !   potential or e.g. electron density.
    ! ------------------------------------------------------------------------- !
    subroutine set_input_vec(mixer,inVec)
        class(scfmixer) :: mixer
        doubleprecision, dimension(:) :: inVec
        integer :: m

        if(SCF_MIXER_SHOW_DEBUG) then
            print*,"mixer:: ",TRIM(mixer_get_name(mixer%MIXER_TYPE)),":: pushing input vector"
        endif

        select case(mixer%MIXER_TYPE)
            case(SCF_MIXER_SIMPLE)

                mixer%vecIn(:,1) = inVec

            case(SCF_MIXER_LOCAL_SIMPLE,SCF_MIXER_ANDERSON,SCF_MIXER_EXTENDED_ANDERSON,SCF_MIXER_BROYDEN)

                do m = mixer%NO_MEM_STEPS , 2 , -1
                    mixer%vecIn(:,m) = mixer%vecIn(:,m-1)
                enddo
                mixer%vecIn(:,1) = inVec

        end select
    endsubroutine set_input_vec

    ! ------------------------------------------------------------------------- !
    ! set_output_vec - same as above, but for the vector which the ouput of the
    !   SCF black box.
    ! ------------------------------------------------------------------------- !
    subroutine set_output_vec(mixer,outVec)
        class(scfmixer) :: mixer
        doubleprecision, dimension(:) :: outVec
        integer :: m

        if(SCF_MIXER_SHOW_DEBUG) then
            print*,"mixer:: ",TRIM(mixer_get_name(mixer%MIXER_TYPE)),":: pushing output vector"
        endif

        select case(mixer%MIXER_TYPE)
            case(SCF_MIXER_SIMPLE)

                mixer%vecOut(:,1) = outVec

            case(SCF_MIXER_LOCAL_SIMPLE,SCF_MIXER_ANDERSON,SCF_MIXER_EXTENDED_ANDERSON,SCF_MIXER_BROYDEN)

                do m = mixer%NO_MEM_STEPS , 2 , -1
                    mixer%vecOut(:,m) = mixer%vecOut(:,m-1)
                enddo
                mixer%vecOut(:,1) = outVec

        end select

    endsubroutine set_output_vec


    ! ------------------------------------------------------------------------- !
    ! mix() - updates the residual vectors and calculates the form of the mixed
    !   vector using propper algorithm. Calulated vector is saved to mixedVec.

    ! References: for more details of method implemented here see:
    ! [1] D. D. Johnson, Phys. Rev. B 38, 12807 (1988)
    ! [2] V Eyert, JOURNAL OF COMPUTATIONAL PHYSICS 124, 271–285 (1996)
    ! [3] O. Certik, Master Thesis: http://www.ondrejcertik.com/media/cookbook/master.pdf
    ! ------------------------------------------------------------------------- !
    ! SOME NOTES:
    ! In order to have comparision with the [2] we use the following notation:
    ! |X_l > = vecIn(:,l)  - where increasing values of l denote older iterations.
    !                       e.g. l=1 keeps values from current iteration.
    ! |Y_l > = vecOut(:,l) - same as above but it keeps black box output vectors
    ! |F_l > = vecRes(:,l) - residual vector
    ! ------------------------------------------------------------------------- !

    subroutine mix(mixer)
        class(scfmixer)  :: mixer
        double precision :: aux1, aux2
        integer          :: i,k

        mixer%CURRENT_ITERATION = mixer%CURRENT_ITERATION + 1

        ! Update the residual vectors
        do k = 1 , mixer%NO_MEM_STEPS
            mixer%vecRes(:,k) = mixer%vecOut(:,k) - mixer%vecIn(:,k)
        enddo


        if(SCF_MIXER_SHOW_DEBUG) then
            print*,"mixer:: ",TRIM(mixer_get_name(mixer%MIXER_TYPE)),":: calculating output vec"
        endif

        ! ---------------------------------------------
        ! Choose the proper method
        ! ---------------------------------------------
        select case(mixer%MIXER_TYPE)


            ! ---------------------------------------------
            ! 1. The simplest possible way to mix vectors
            ! ---------------------------------------------
            case(SCF_MIXER_SIMPLE)
                ! Compare it to Eq. (3.1) from Ref. [2].
                if(mixer%CURRENT_ITERATION == 1) then
                    mixer%mixedVec = mixer%vecOut(:,1)
                else
                    mixer%mixedVec = mixer%vecIn(:,1) + mixer%omega * mixer%vecRes(:,1)
                endif

            ! ---------------------------------------------
            ! 2. For the description of the simple local
            !    mixing scheme see [3] Eq (3.28).
            ! ---------------------------------------------
            case(SCF_MIXER_LOCAL_SIMPLE)
                if(mixer%CURRENT_ITERATION == 1) then
                    mixer%mixedVec = mixer%vecOut(:,1)
!                else if(mixer%CURRENT_ITERATION < 2) then ! for first two iters use simple scheme
!                    mixer%mixedVec = mixer%vecIn(:,1) + mixer%omega * mixer%vecRes(:,1)
                else ! otherwise do the local mixing according the algorithm described in [3]
                    mixer%mixedVec = 0
                    do i = 1 , mixer%NO_SITES
                        aux1 = mixer%vecRes(i,2)*mixer%vecRes(i,1)
                        if(aux1 > 0) then
                            mixer%vecLocalOmega(i) = mixer%vecLocalOmega(i) + mixer%omega
                            if(mixer%vecLocalOmega(i) > mixer%max_omega) mixer%vecLocalOmega(i) = mixer%max_omega
                        else
                            mixer%vecLocalOmega(i) = mixer%omega
                        endif
                        aux2 = mixer%vecLocalOmega(i)
                        mixer%mixedVec(i) = mixer%vecIn(i,1) + aux2 * mixer%vecRes(i,1)
                    enddo
                endif ! end of if (iter > 2)


            ! ---------------------------------------------
            ! 3. We use the implementaion proposed in [2]
            ! ---------------------------------------------
            case(SCF_MIXER_ANDERSON)
                call anderson_method(mixer)

            ! ---------------------------------------------
            ! 4. We use the implementaion proposed in [2]
            ! ---------------------------------------------
            case(SCF_MIXER_EXTENDED_ANDERSON)
                call extended_anderson_method(mixer)

            ! ---------------------------------------------
            ! 5. We use the implementaion proposed in [3]
            ! ---------------------------------------------
            case(SCF_MIXER_BROYDEN)
                call broyden_method(mixer)

        end select


    end subroutine mix


    ! ------------------------------------------------------------------------- !
    ! anderson_method() - implementation of the anderson mixing scheme taken from
    !       Reference [1]
    ! [1] V Eyert, JOURNAL OF COMPUTATIONAL PHYSICS 124, 271–285 (1996)
    ! ------------------------------------------------------------------------- !
    subroutine anderson_method(mixer)
        class(scfmixer)  :: mixer
        integer :: niter,i,j,k
        doubleprecision :: aux1, aux2
        doubleprecision, allocatable :: matTheta(:,:) , vecTheta(:,:)

        if(mixer%CURRENT_ITERATION == 1) then
            mixer%mixedVec = mixer%vecOut(:,1)
        else if(mixer%NO_MEM_STEPS == 1) then
            mixer%mixedVec = mixer%vecIn(:,1) + mixer%omega * mixer%vecRes(:,1)
        else


            ! If iteration number is smaller than NO_MEM_STEPS creat
            ! system of size niter.
            niter = min(mixer%NO_MEM_STEPS-1,mixer%CURRENT_ITERATION)

            allocate(matTheta(niter,niter))
            allocate(vecTheta(niter,1))

            ! Calculate the theta vector from (4.3) Ref [1].
            ! We enumerate the vectors from 1 not from 0 like in C++
            matTheta = 0
            do i = 1 , niter
            do j = 1 , niter
                matTheta(i,j) = sum( (mixer%vecRes(:,1)-mixer%vecRes(:,1+i)) * &
                                     (mixer%vecRes(:,1)-mixer%vecRes(:,1+j)))
            enddo
                vecTheta(i,1) = sum( ( mixer%vecRes(:,1)-mixer%vecRes(:,1+i) )*mixer%vecRes(:,1) )
            enddo


            ! solve system of linear equations (4.3)
            call dgaussj(matTheta,niter,niter,vecTheta,1,1)


            do i = 1 , mixer%NO_SITES
                aux1 = mixer%vecIn (i,1)
                aux2 = mixer%vecRes(i,1)
                do k = 1 , niter
                    ! Eq. (4.1)
                    aux1 = aux1 + ( vecTheta(k,1) ) * (mixer%vecIn (i,k+1)-mixer%vecIn (i,1))
                    ! Eq. (4.2)
                    aux2 = aux2 + ( vecTheta(k,1) ) * (mixer%vecRes(i,k+1)-mixer%vecRes(i,1))
                enddo
                ! Eq. (4.4)
                mixer%mixedVec(i) = aux1 + mixer%omega * aux2
            enddo


            deallocate(matTheta)
            deallocate(vecTheta)


        endif

    end subroutine anderson_method



    ! ------------------------------------------------------------------------- !
    ! extended_anderson_method() - implementation of the anderson mixing scheme
    !   taken from Reference [1]
    !
    ! [1] V Eyert, JOURNAL OF COMPUTATIONAL PHYSICS 124, 271–285 (1996)
    ! ------------------------------------------------------------------------- !
    subroutine extended_anderson_method(mixer)
        class(scfmixer)  :: mixer
        integer :: niter,n,m,k,i
        doubleprecision :: aux1, aux2
        doubleprecision, allocatable :: matGamma(:,:) , vecGamma(:,:)


        if(mixer%CURRENT_ITERATION == 1) then
            mixer%mixedVec = mixer%vecOut(:,1)
        else if(mixer%CURRENT_ITERATION < 3 .or. mixer%NO_MEM_STEPS == 1) then
            mixer%mixedVec = mixer%vecIn(:,1) + mixer%omega * mixer%vecRes(:,1)
        else


            ! If iteration number is smaller than NO_MEM_STEPS creat
            ! system of size niter.
            niter = min(mixer%NO_MEM_STEPS-1,mixer%CURRENT_ITERATION)

            allocate(matGamma(niter,niter))
            allocate(vecGamma(niter,1))

            ! Calculate the theta vector from (7.6) Ref [1].
            ! We enumerate the vectors from 1 not from 0 like in C++
            matGamma = 0
            do n = niter + 1 , 2 , -1 ! from l-M to l-1 = 2, because last iteration has l=1
            do m = niter + 1 , 2 , -1
                matGamma(n-1,m-1) = sum( (mixer%vecRes(:,n-1)-mixer%vecRes(:,n))  * &
                                         (mixer%vecRes(:,m-1)-mixer%vecRes(:,m))) * &
                                         ( 1.0 + delta(n,m)*(mixer%omega0)**2 )
            enddo
                vecGamma(n-1,1)    = sum( ( mixer%vecRes(:,n-1)-mixer%vecRes(:,n) )*mixer%vecRes(:,1) )
            enddo


            ! solve system of linear equations (7.6)
            call dgaussj(matGamma,niter,niter,vecGamma,1,1)

            ! calculation of new vector from Eq. (7.7)
            do i = 1 , mixer%NO_SITES

                mixer%mixedVec(i) = mixer%vecIn (i,1) + mixer%omega * mixer%vecRes(i,1)
                aux1 = 0
                aux2 = 0

                do m = niter + 1 , 2 , -1
                    aux1 = aux1 + ( vecGamma(m-1,1) ) * (mixer%vecIn (i,m-1)-mixer%vecIn (i,m))
                    aux2 = aux2 + ( mixer%omega * vecGamma(m-1,1) ) * (mixer%vecRes(i,m-1)-mixer%vecRes(i,m))
                enddo

                mixer%mixedVec(i) = mixer%mixedVec(i) - (aux1 + aux2)
            enddo


            deallocate(matGamma)
            deallocate(vecGamma)


        endif

    end subroutine extended_anderson_method


    ! ------------------------------------------------------------------------- !
    ! broyden method() - based on the [3]
    ! ------------------------------------------------------------------------- !
    subroutine broyden_method(mixer)
        class(scfmixer)  :: mixer
        integer :: niter,i,j
        doubleprecision fnorm
        doubleprecision, allocatable :: vecu(:,:) , vecv(:,:)


        niter = min(mixer%NO_MEM_STEPS-1,mixer%CURRENT_ITERATION-1)
        allocate(vecu(mixer%NO_SITES,niter),vecv(mixer%NO_SITES,niter))


        do i=1,niter
          fnorm=sum((mixer%vecRes(:,i)-mixer%vecRes(:,i+1))**2)
          vecv(:,i)=(mixer%vecRes(:,i)-mixer%vecRes(:,i+1))/(fnorm)
        enddo


        do i=niter,1,-1
          vecu(:,i)=(mixer%vecIn(:,i)-mixer%vecIn(:,i+1))+mixer%omega*(mixer%vecRes(:,i)-mixer%vecRes(:,i+1))


          do j=i+1,niter
            vecu(:,i)=vecu(:,i)-&
            vecu(:,j)*sum(vecv(:,j)*(mixer%vecRes(:,i)-mixer%vecRes(:,i+1)))
          enddo

        enddo

        mixer%mixedVec = mixer%vecIn(:,1) + mixer%omega * mixer%vecRes(:,1)

        do i=1,niter
          mixer%mixedVec=mixer%mixedVec-&
          vecu(:,i)*sum(vecv(:,i)*mixer%vecRes(:,1))
        enddo

        deallocate(vecu,vecv)


    end subroutine broyden_method

! ----------------------------------------------------------------------------------
!                                       TOOLS
! ----------------------------------------------------------------------------------
! From the numerical recipes
SUBROUTINE dgaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np
      doubleprecision a(np,np),b(np,mp)
      integer,PARAMETER :: NMAX=500
      INTEGER :: i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      doubleprecision :: big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge. abs(big))then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END subroutine dgaussj



   ! Kronecker delta
    integer function delta(i,j)
        integer :: i,j

        delta = 0
        if( i == j ) delta = 1
    end function delta

end module modmixers
