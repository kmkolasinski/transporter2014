MODULE spinmodpop
    USE, INTRINSIC :: ISO_C_BINDING
    use            :: modjed
    use            :: modutils
    implicit none
    private

    double precision :: DX  , Ef , t0 , tc0 , BZ , BX , BY , hny , Emax , KIN
    integer          :: N , wypisz = 0



    complex*16,dimension(:,:),allocatable        :: Hamiltonian
    integer,dimension(:,:),allocatable           :: GINDEX ! (wezek,spin)
    double precision,dimension(:),allocatable    :: UVEC,Ey ! potencjal i pole elektryczne
    double precision,dimension(:,:),allocatable  :: SUVEC ! z potencjalem do spinu
    complex*16,dimension(:,:),allocatable        :: Chi_m_in ! mod wchodzacy do ukladu
    complex*16,dimension(:),allocatable          :: K_m_in
    complex*16,dimension(:,:),allocatable        :: Chi_m_out ! wychodzacy
    complex*16,dimension(:),allocatable          :: K_m_out

    complex*16,dimension(:,:,:),allocatable    :: ChiMod    ! (x,mod,kierunek) - spinor
    complex*16,dimension(:,:),allocatable      :: ChiKvec   ! (mod,kierunek)
    complex*16,dimension(:,:),allocatable      :: ChiLambda ! (mod,kierunek)
    integer                                    :: LICZBA_MODOW,L_M,LICZBA_MODOW_EVANESCENTYCH


    ! -------------------------------------------------
    !                    LAPACK
    ! -------------------------------------------------
    !     .. Local Scalars ..
      INTEGER          INFO, LWORK, LRWORK, LIWORK, IL, IU, M , LWMAX
      DOUBLE PRECISION ABSTOL, VL, VU

!     .. Local Arrays ..
      INTEGER,allocatable,dimension(:)          :: ISUPPZ, IWORK
      DOUBLE PRECISION,allocatable,dimension(:) :: MODPOP_VALS( : ), RWORK( : )
      COMPLEX*16,allocatable,dimension(:)       :: Z( :, : ), WORK( : )


      public :: spinmodpop_inicjalizacja
      public :: spinmodpop_relacja_dyspersji
      public :: spinmodpop_zwalnienie_pamieci
      public :: spinmodpop_calc_modes_from_wfm!(pDx,pN,pEf,pB,pUvec,pbHorizontal)
      public :: spinmodpop_calc_TB_current
      public :: LICZBA_MODOW,L_M,LICZBA_MODOW_EVANESCENTYCH,ChiMod,ChiKvec,ChiLambda,MODPOP_VALS
      public :: SUVEC
    contains

! ------------------------------------------------------------------------------------ !
! Inicjalizcja problemu wlasnego dla relacji dyspersji ze spinami (rashba + lateral)
! podajemy skok siatki w nm
! liczba oczek siatki w kierunku Y
! energie Fermiego w mevach
! pole magnetyczne w teslach
! potencjal w mev
! ------------------------------------------------------------------------------------ !
subroutine spinmodpop_inicjalizacja(pDx,pN,pEf,pB,pUvec)
    double precision,intent(in)               :: pDX
    integer,intent(in)                        :: pN
    double precision,intent(in)               :: pEf
    double precision,intent(in)               :: pB
    double precision,intent(in),dimension(pN) :: pUVEC
    integer    :: num , iter


    DX    = pDX*L2LR
    N     = pN-2
    Ef    = pEf/Rd/1000 ! bo Ef w meV
    LWMAX = 100*N*2
    t0    = 0.5/DX/DX
    BZ    = BtoDonorB(pB)
    BX    = BtoDonorB(atomic_Bx)
    BY    = BtoDonorB(atomic_By)
    hny   = (N+1)/2.0!

    if(TRANS_DEBUG)then
        print*,"Relacja dyspersji:"
        print*,"    N  :",N
        print*,"    hny:",hny
        print*,"    t0 :",t0
        print*,"    Rd :",Rd
        print*,"    N  :",N
        print*,"    Ef :",Ef
        print*,"    DX :",DX
        print*,"    hny:",hny
        print*,"    B  :",pB
        print*,"    Kmax:",3.14159/pDX
    endif

    allocate(Uvec(N))
    allocate(Ey(N))
    allocate(GINDEX(N,2))
    allocate(Hamiltonian(2*N,2*N))

    iter = 1
    do num = 1 , N
        Uvec(num)        = pUVEC(num+1)/Rd/1000 ! bo Ef w meV
        GINDEX(num,1) = iter
        GINDEX(num,2) = iter + N
        iter = iter + 1
    enddo

    ! obliczanie natezenia pola elektrycznego

    do num = 2 , N-1
        Ey(num) = (uvec(num+1) - uvec(num-1))/2/dx
    enddo
    Ey(1) = Ey(2)
    Ey(N) = Ey(N-1)


    Hamiltonian = 0

    allocate(ISUPPZ( 2*N ))
    allocate(IWORK( LWMAX ))
    allocate(MODPOP_VALS( 2*N ))
    allocate(RWORK( LWMAX ))
    allocate(Z( 2*N, 2*N ))
    allocate(WORK( LWMAX ))
    !
    !     Query the optimal workspace.
    !
    LWORK  = -1
    LRWORK = -1
    LIWORK = -1
!
    ABSTOL = -1.0
    VL     = -2*Ef
    VU     = +2*Ef

    call spinmodpop_utworz_hamiltonian((0.0D0,0.0D0))

    open(unit = 33345, file= "Hamiltonian.txt" )
    do iter = 1 , 2*N
    write(33345,"(5000f10.4)"),dble(Hamiltonian(iter,:))
    enddo
    close(33345)

    CALL ZHEEVR( "N", 'Values', 'Lower', 2*N, Hamiltonian, 2*N, VL, VU, IL,&
    &             IU, ABSTOL, M, MODPOP_VALS, Z, 2*N, ISUPPZ, WORK, LWORK, RWORK,&
    &             LRWORK, IWORK, LIWORK, INFO )

    if( INFO /= 0 ) then
        print*,"  SPINMODPOP: Problem zle zdefiniowany INFO:",INFO
        stop
    endif
    LWORK  = MIN( LWMAX, INT( WORK( 1 ) ) )
    LRWORK = MIN( LWMAX, INT( RWORK( 1 ) ) )
    LIWORK = MIN( LWMAX, IWORK( 1 ) )



end subroutine spinmodpop_inicjalizacja


! --------------------------------------------
! Tworzenie hamiltonianu
! --------------------------------------------
subroutine spinmodpop_utworz_hamiltonian(pk)
    complex*16 :: pk
    integer :: i,j
    complex*16 ::  Dy , Cx , Kx , Kmx
    complex*16 :: tc0

    tc0 = so_rashba*0.5/dx
    Hamiltonian = 0
    do i = 1 , N
        ! spin up
        Kx  = exp(+pk*DX)
        Kmx = exp(-pk*DX)
        Cx  = exp(II*DX*DX*BZ*(i-hNY))
        Dy  = exp(II*Ey(i)*so_loc*DX)

        Hamiltonian(GINDEX(i,1),GINDEX(i,1)) = UVEC(i) + SUVEC(i,+1) + 4*t0 + 0.5 * 0.5 * G_LAN * M_EFF * BZ - t0*Kmx*Cx*conjg(Dy) - t0*Kx*conjg(Cx)*Dy
        Hamiltonian(GINDEX(i,1),GINDEX(i,2)) = tc0*Kx*conjg(Cx) - tc0*Kmx*Cx + 0.5 * 0.5 * G_LAN * M_EFF * (Bx  - II*By)

        if ( i > 1 ) then
            Hamiltonian(GINDEX(i,1),GINDEX(i-1,1)) = -t0
            Hamiltonian(GINDEX(i,1),GINDEX(i-1,2)) = II*tc0
        endif
        if ( i < N ) then
            Hamiltonian(GINDEX(i,1),GINDEX(i+1,1)) = -t0
            Hamiltonian(GINDEX(i,1),GINDEX(i+1,2)) = -II*tc0
        endif

        ! spin down
        Hamiltonian(GINDEX(i,2),GINDEX(i,2)) = UVEC(i) + SUVEC(i,+1) + 4*t0 - 0.5* 0.5 * G_LAN * M_EFF * BZ - t0*Kmx*Cx*(Dy) - t0*Kx*conjg(Cx)*conjg(Dy)
        Hamiltonian(GINDEX(i,2),GINDEX(i,1)) = -tc0*Kx*conjg(Cx) + tc0*Kmx*Cx + 0.5* 0.5 * G_LAN * M_EFF * (Bx  + II*By)

        if ( i > 1 ) then
            Hamiltonian(GINDEX(i,2),GINDEX(i-1,2)) = -t0
            Hamiltonian(GINDEX(i,2),GINDEX(i-1,1)) = II*tc0
        endif
        if ( i < N ) then
            Hamiltonian(GINDEX(i,2),GINDEX(i+1,2)) = -t0
            Hamiltonian(GINDEX(i,2),GINDEX(i+1,1)) = -II*tc0
        endif
!        write(222,"(i,A,10f16.6)"),i,"up=", Kmx*Cx*conjg(Dy) + Kx*conjg(Cx)*Dy - (Kmx*Cx*(Dy) + Kx*conjg(Cx)*conjg(Dy)) , Kmx*conjg(Dy) + Kx*Dy , Kmx*(Dy) + Kx*conjg(Dy)

    enddo




!    Kmx = 0
!    do i = 1 , 2 * N
!    do j = i , 2 * N
!
!        Kmx = Kmx + abs( Hamiltonian(i,j) - conjg(Hamiltonian(j,i)) )
!
!    enddo
!    enddo
!    print*,abs(pk),abs(Kmx)

endsubroutine spinmodpop_utworz_hamiltonian


subroutine spinmodpop_relacja_dyspersji(pkmin,pkmax,pdk,pEmax,nazwa_pliku)
    double precision :: pkmin,pkmax,pdk,pEmax

    character(*) :: nazwa_pliku
    double precision :: kmin,kmax,dk,Emax
    complex*16 :: kvec
    doubleprecision :: skank, YA, YB, polaryzacje(1000,3) , energie(1000)
    integer :: i, max_mods
    complex*16 , dimension(:) , allocatable :: XYup , XYdwn , Zup , Zdwn




    open(unit = 12321, file= nazwa_pliku )
    open(unit = 12322, file= nazwa_pliku(1:len(nazwa_pliku)-4)//"_polar_z.txt")
    open(unit = 12323, file= nazwa_pliku(1:len(nazwa_pliku)-4)//"_polar_x.txt")
    open(unit = 12324, file= nazwa_pliku(1:len(nazwa_pliku)-4)//"_polar_y.txt")

    kmax = pkmax*LR2L
    kmin = pkmin*LR2L
    dk   = pdk*LR2L
    Emax = pEmax/Rd/1000.0
    allocate(XYup(N))
    allocate(XYdwn(N))
    allocate(Zup(N))
    allocate(Zdwn(N))

    kvec = II*0.0
    call spinmodpop_utworz_hamiltonian(kvec)
    ABSTOL = -1.0
    !     Set VL, VU to compute eigenvalues in half-open (VL,VU] interval
    VL = -Emax
    VU = Emax
    !
    !     Solve eigenproblem.
    !
    CALL ZHEEVR( "V", 'Values', 'Lower', 2*N, Hamiltonian, 2*N, VL, VU, IL,&
    &             IU, ABSTOL, M, MODPOP_VALS, Z, 2*N, ISUPPZ, WORK, LWORK, RWORK,&
    &             LRWORK, IWORK, LIWORK, INFO )
    max_mods = M

    do skank = kmin , kmax , dk

        kvec = II*skank
        call spinmodpop_utworz_hamiltonian(kvec)
        ABSTOL = -1.0
        !     Set VL, VU to compute eigenvalues in half-open (VL,VU] interval
        VL = -2*Emax
        VU = Emax
        !
        !     Solve eigenproblem.
        !
        MODPOP_VALS = 0
        CALL ZHEEVR( "V", 'Values', 'Lower', 2*N, Hamiltonian, 2*N, VL, VU, IL,&
        &             IU, ABSTOL, M, MODPOP_VALS, Z, 2*N, ISUPPZ, WORK, LWORK, RWORK,&
        &             LRWORK, IWORK, LIWORK, INFO )
!
        polaryzacje = -100
        energie     = 0
        do i = 1 , M
            Zup  = 0
            Zdwn = 0
            Zup (1:N) = Z(1:N    ,i)
            Zdwn(1:N) = Z(N+1:2*N,i)

            ! Kierunek Z
            YA = sum(abs(Zup)**2)
            YB = sum(abs(Zdwn)**2)
            polaryzacje(i,1) = (YA-YB)/(YA+YB) ! kierunek z

            ! Kierunek X
            YA = sum(abs(Zup+Zdwn)**2) ! up
            YB = sum(abs(Zup-Zdwn)**2) ! down
            polaryzacje(i,2) = (YA-YB)/(YA+YB) ! kierunek x

            ! Kierunek Y
            YA = sum(abs(Zup+II*Zdwn)**2)
            YB = sum(abs(Zup-II*Zdwn)**2)
            polaryzacje(i,3) = (YA-YB)/(YA+YB) ! kierunek y

            energie(i) = MODPOP_VALS(i)*Rd*1000.0
        enddo


        write(12321,"(1000e20.10)"),skank*L2LR,atomic_Ef,M+0.0,energie(1:max_mods)
        write(12322,"(1000e20.10)"),skank*L2LR,atomic_Ef,M+0.0,polaryzacje(1:max_mods,1)
        write(12323,"(1000e20.10)"),skank*L2LR,atomic_Ef,M+0.0,polaryzacje(1:max_mods,2)
        write(12324,"(1000e20.10)"),skank*L2LR,atomic_Ef,M+0.0,polaryzacje(1:max_mods,3)

    enddo

    if(allocated(GINDEX))      print*,GINDEX

    close(12321)
    close(12322)
    close(12323)
    close(12324)

    deallocate(XYup)
    deallocate(XYdwn)
    deallocate(Zup)
    deallocate(Zdwn)

end subroutine spinmodpop_relacja_dyspersji


subroutine spinmodpop_zwalnienie_pamieci()

    if(TRANS_DEBUG == .true. ) print*,"    SPINMODPOP:Zwalnianie pamieci."
    if(allocated(UVEC))        deallocate(Uvec)



    if(allocated(Ey))          deallocate(Ey)



    if(allocated(Hamiltonian)) deallocate(Hamiltonian)




    if(allocated(ISUPPZ))deallocate(ISUPPZ)
    if(allocated(IWORK ))deallocate(IWORK )
    if(allocated(MODPOP_VALS))deallocate(MODPOP_VALS)
    if(allocated(RWORK ))deallocate(RWORK )
    if(allocated(Z     ))deallocate(Z     )
    if(allocated(WORK  ))deallocate(WORK  )

    if(allocated(Chi_m_IN))  deallocate(Chi_m_IN)
    if(allocated(K_m_IN))    deallocate(K_m_IN)
    if(allocated(Chi_m_OUT)) deallocate(Chi_m_OUT)
    if(allocated(K_m_OUT))   deallocate(K_m_OUT)

    if(allocated(ChiMod ))    deallocate(ChiMod)
    if(allocated(ChiKvec))    deallocate(ChiKvec)
    if(allocated(ChiLambda))  deallocate(ChiLambda)

    if(allocated(GINDEX))      print*,GINDEX
    if(allocated(GINDEX))      deallocate(GINDEX)
    if(TRANS_DEBUG == .true. ) print*,"    SPINMODPOP:Zwalnianie pamieci zakonczone."


end subroutine spinmodpop_zwalnienie_pamieci



! -----------------------------------------------------------------------------------------------
!                      Wyznaczanie modow metoda Wave-Function-Matching
! Podajemy krok siatki w nm, rozmiar wejscia w oczkach siatki (wliczajac punkty gdzie
! f.f. przyjmuje wartosc zero - warunek Dirichleta).
! pEf   - energia Fermiego w meV
! pB    - pole magnetyczne w Teslach
! pUvec - wektor z potencjalem w danym wejsciu
!
! Funkcja wyznacza wszystkie mody lezace ponizej energii Fermiego plus wszystkie
! mody evanescentne. Funkcja zostala napisana na podstawie artykulu:
! "Calculating Scattering Matrices by Wave Function Matching" sekcja:
! 1.1.3 Wave function matching in three dimensions, wykorzystuje przede wszystkim wzor (52)
! -----------------------------------------------------------------------------------------------
subroutine spinmodpop_calc_modes_from_wfm(pDx,pN,pEf,pB,pUvec,pbHorizontal)
    double precision,intent(in)               :: pDX
    integer,intent(in)                        :: pN
    double precision,intent(in)               :: pEf
    double precision,intent(in)               :: pB
    double precision,intent(in),dimension(pN) :: pUVEC
    logical, intent(in) :: pbHorizontal

    ! zmienne pomocnicze
    complex*16, allocatable , dimension(:,:)   :: MB,MA
    complex*16, allocatable , dimension(:,:,:) :: Mdiag,Mham,Mtau,MatS,MatBs
    integer :: i,num_in,num_out,iter

    INTEGER                                      :: LDVL, LDVR , LWMAX , LWORK , INFO
    COMPLEX*16 , dimension(:) ,     allocatable  :: ALPHA , BETA , WORK
    double precision, dimension(:), allocatable  :: RWORK
    double precision, dimension(:,:), allocatable  :: spinuvec
    COMPLEX*16 , dimension(:) ,     allocatable  :: tmp_vec , tmp_vec2
    COMPLEX*16 , dimension(:,:)   , allocatable  :: Z
    COMPLEX*16      :: DUMMY(1,1),lambda , YcY
    complex*16      :: kvec , Cx , Sy
    doubleprecision :: dkvec , polarA,polarB


    ! konwersja jednostek do jednostek donorowych
    DX    = pDX*L2LR
    N     = pN-2
    Ef    = pEf/Rd/1000 ! bo Ef w meV
    t0    = 0.5/DX/DX
    tc0   = 0.5/DX
    BZ    = BtoDonorB(pB)
    BX    = BtoDonorB(atomic_Bx)
    BY    = BtoDonorB(atomic_By)
    hny   = (N+1)/2.0! dobor odpowieniego cechowania



    ! alokacja tablic

    if(allocated(Uvec))       deallocate(Uvec)
    if(allocated(spinuvec))       deallocate(spinuvec)
    if(allocated(Ey))  deallocate(Ey)
    if(allocated(GINDEX))  deallocate(GINDEX)
    if(allocated(Hamiltonian))  deallocate(Hamiltonian)
    if(allocated(Mtau))  deallocate(Mtau)
    if(allocated(Mdiag))  deallocate(Mdiag)
    if(allocated(Mham))  deallocate(Mham)
    if(allocated(MatS))  deallocate(MatS)
    if(allocated(MatBs))  deallocate(MatBs)
    if(allocated(tmp_vec))  deallocate(tmp_vec)
    if(allocated(tmp_vec2))  deallocate(tmp_vec2)
    if(allocated(MA))  deallocate(MA)
    if(allocated(MB))  deallocate(MB)
    if(allocated(Z))  deallocate(Z)
    if(allocated(ALPHA))  deallocate(ALPHA)
    if(allocated(BETA))  deallocate(BETA)
    if(allocated(RWORK))  deallocate(RWORK)
    if(allocated(WORK))  deallocate(WORK)


    allocate(Uvec(N))
    allocate(spinuvec(N,-1:1))
    allocate(Ey(N))
    allocate(GINDEX(N,2))
    allocate(Hamiltonian(2*N,2*N))
    allocate(Mtau (N,N,2))
    allocate(Mdiag(N,N,2))
    allocate(Mham (N,N,2))
    allocate(MatS (N,N,2))
    allocate(MatBs(N,N,2))
    allocate(tmp_vec (2*N))
    allocate(tmp_vec2(2*N))
    allocate(MA (4*N,4*N))
    allocate(MB (4*N,4*N))
    allocate(Z  (4*N,4*N))

    ! Ustalenie parametrow LAPACKA
    LWMAX = 50 * N * 2
    LDVL  = 4  * N
    LDVR  = 4  * N
    ! Alokacja macierzy LAPACKA
    allocate(ALPHA(4*N))
    allocate(BETA(4*N))
    allocate(RWORK(16*N))
    allocate(WORK(LWMAX))


    do i = 1 , N
        Uvec(i)        = pUVEC(i+1)/Rd/1000 ! bo Ef w meV
        spinuvec(i,:)  = SUVEC(i+1,:) ! bo Ef w donorowych
    enddo



    iter = 1
    do i = 1 , N
        GINDEX(i,1) = iter
        GINDEX(i,2) = iter + N
        iter = iter + 1
    enddo

    do i = 2 , N-1
        Ey(i) = (uvec(i+1) - uvec(i-1))/2/dx
    enddo
    Ey(1) = Ey(2)
    Ey(N) = Ey(N-1)






    ! Tworzenie rownania wlasnego  - wzor (52)
    Mdiag = 0
    Mham  = 0
    Mtau  = 0
    MatBs = 0
    MatS  = 0

! Przygotowanie podmacierzy diagonalnej 1 i macierzy B
!        Cx  = exp(II*DX*DX*BZ*(i-hNY))
!        Dy  = exp(II*Ey(i)*so_loc*DX)

    if(pbHorizontal) then
        do i = 1 , N
            Mdiag(i,i,1) =     1 ! up
            Mdiag(i,i,2) =     1 ! down
            Mtau (i,i,1) =  - t0 * exp(+II*(DX*DX*BZ*(i-hNY))) * exp(-II*DX*Ey(i)*so_loc)
            Mtau (i,i,2) =  - t0 * exp(+II*(DX*DX*BZ*(i-hNY))) * exp(+II*DX*Ey(i)*so_loc)
            MatBs(i,i,1) =  so_rashba*tc0 * exp(+II*(DX*DX*BZ*(i-hNY)))
            MatBs(i,i,2) = -so_rashba*tc0 * exp(+II*(DX*DX*BZ*(i-hNY)))
        enddo
    else
        print*,"SpinModPop::Wersja nie obslugiwana!"
        stop
        do i = 1 , N
            Mdiag(i,i,1) =     1
            Mtau (i,i,1) =  - t0 * exp(-II*(DX*DX*BZ*(i-hNY)))
        enddo
    endif


! Przygotowanei podmacierzy (E-H) i macierzy od raszby
    do i = 1 , N
        Mham(i,i,1)   = UVEC(i) + spinuvec(i,+1) + 4*t0 + 0.5 * 0.5 * G_LAN * M_EFF * BZ - Ef
        Mham(i,i,2)   = UVEC(i) + spinuvec(i,-1) + 4*t0 - 0.5 * 0.5 * G_LAN * M_EFF * BZ - Ef
        if(i < N) Mham(i,i+1,:) = -t0
        if(i > 1) Mham(i,i-1,:) = -t0

        if(i < N) MatS(i,i+1,:) = -II*so_rashba*tc0
        if(i > 1) MatS(i,i-1,:) = +II*so_rashba*tc0
        ! zeeman spin down
        MatS(i,i,2) = 0.5 * 0.5 * G_LAN * M_EFF * ( BX - II * BY)
        MatS(i,i,1) = 0.5 * 0.5 * G_LAN * M_EFF * ( BX + II * BY)
    enddo

    ! Wypelnienie macierzy:
    MA  = 0
    MB  = 0

    ! spin up
    MA(1:N     ,1:N)     = Mham(:,:,1)
    MA(N+1:2*N ,1:N)     = Mdiag(:,:,1)
    MA(1:N,N+1:2*N)      = Mtau(:,:,1)

    MB(1:N,1:N)          = - conjg(Mtau(:,:,1))
    MB(N+1:2*N,N+1:2*N)  = Mdiag(:,:,1)

    ! rashba down
    MA(1:N,2*N+1  :2*N+N)   = MatS (:,:,2)
    MA(1:N,2*N+1+N:2*N+2*N) = MatBs(:,:,2)
    MB(1:N,2*N+1  :2*N+N)   = conjg(MatBS(:,:,2))

    ! spin down
    MA(2*N+1:2*N+N     ,2*N+1:2*N+N)     = Mham(:,:,2)
    MA(2*N+N+1:2*N+2*N ,2*N+1:2*N+N)     = Mdiag(:,:,2)
    MA(2*N+1:2*N+N,2*N+N+1:2*N+2*N)      = Mtau(:,:,2)

    MB(2*N+1:2*N+N,2*N+1:2*N+N)          = - conjg(Mtau(:,:,2))
    MB(2*N+N+1:2*N+2*N,2*N+N+1:2*N+2*N)  = Mdiag(:,:,2)

    ! rashba up
    MA(2*N+1:2*N+N,1:N)     = MatS (:,:,1)
    MA(2*N+1:2*N+N,N+1:2*N) = MatBs(:,:,1)
    MB(2*N+1:2*N+N,1:N)     = conjg(MatBS(:,:,1))

!    print*,"Tworzenie macierzy=",get_clock()


!    open(unit = 33345, file= "MA.txt" )
!    do i = 1 , 4*N
!    write(33345,"(5000f10.4)"),abs(MA(i,:))
!    enddo
!    close(33345)
!    open(unit = 33345, file= "MB.txt" )
!    do i = 1 , 4*N
!    write(33345,"(5000f10.4)"),abs(MB(i,:))
!    enddo
!    close(33345)
    LWORK = -1

!        przygotowanie
    CALL ZGGEV("N","N", 4*N, MA, 4*N, MB,4*N, ALPHA,BETA, &
&   DUMMY, 1, Z, LDVR, WORK, LWORK, RWORK, INFO )

    LWORK = MIN(LWMAX, INT( WORK(1)))
!    print*,"LMAX Query=",get_clock()

    ! rozwiazanie problemu wlasnego
    CALL ZGGEV("N","V", 4*N, MA, 4*N, MB , 4*N, ALPHA,BETA, &
&   DUMMY, 1, Z, LDVR, WORK, LWORK, RWORK, INFO )

    ! sprawdzamy czy uklad zostal skontruowany poprawnie
    if( INFO /= 0 ) then
        print*,"Modpop::WFM::Nie udalo sie poprawnie znalezc modow"
        stop
    endif

!    print*,"Rozwiazanie rownania wlasnego=",get_clock()
!do i = 1 , 4*N
!if(abs(Beta(i))>1e-16) then
!lambda = (ALPHA(i)/BETA(i))
!    print*,i,abs(lambda),(log(lambda)/II/DX)
!endif
!enddo

    LICZBA_MODOW = 0
    do i = 1 , 4*N
        if(abs(Beta(i))>1e-16) then ! zgodnie z przykladem LApacka
            lambda= (ALPHA(i)/BETA(i))
            YcY = DX*(sum(abs(Z(N+1:2*N,i))**2) + sum(abs(Z(2*N+N+1:2*N+2*N,i))**2))
            Z(:,i) = Z(:,i)/sqrt(YcY)
            if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then ! jesli nie mod evanescentny
                !kvec = (log(lambda)/II/DX)   ! spradzamy znak Kvec zeby policzyc ile jest wektorow
                                             ! w jedna strone


                dkvec = (log(lambda)/II/DX)

!                kvec  = spinmodpop_calc_TB_current(dkvec,Z(N+1:2*N,i),Z(2*N+N+1:2*N+2*N,i),N)
                kvec  = spinmodpop_calc_current(dkvec,Z(N+1:2*N,i),Z(2*N+N+1:2*N+2*N,i),N)

                if(dble(kvec) > 0) LICZBA_MODOW = LICZBA_MODOW + 1
            endif
        endif
    enddo




   ! if(LICZBA_MODOW == 0) return;

    if(allocated(ChiMod ))    deallocate(ChiMod)
    if(allocated(ChiKvec))    deallocate(ChiKvec)
    if(allocated(ChiLambda))  deallocate(ChiLambda)

    allocate(ChiMod (2*(N+2),2*N,2))
    allocate(ChiKvec(2*N,2))
    allocate(ChiLambda(2*N,2))


    ChiMod      = 0
    ChiKvec     = 0
    ChiLambda   = 0
    if( TRANS_DEBUG ) then
        print*,""
        print*,"WEKTORY FALOWE (IN/OUT):"
    endif
    num_in  = 0
    num_out = 0

    ! szukamy spinorow
    do i = 1 , 4*N
        if(abs(Beta(i))>1e-16) then
            lambda = (ALPHA(i)/BETA(i))

            if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then
                kvec         = (log(lambda)/DX)


                dkvec = (log(lambda)/II/DX)

                dkvec = spinmodpop_calc_current(dkvec,Z(N+1:2*N,i),Z(2*N+N+1:2*N+2*N,i),N)
!                dkvec = spinmodpop_calc_TB_current(dkvec,Z(N+1:2*N,i),Z(2*N+N+1:2*N+2*N,i),N)

                if( dkvec > 0) then
                    num_in                   = num_in + 1

                    ChiKvec(num_in,1)          = kvec
                    ChiLambda(num_in,1)        = lambda
                    ChiMod(2  :N+1,num_in,1)   = Z(N+1:2*N,i)
                    ChiMod(N+4:2*N+3,num_in,1) = Z(2*N+N+1:2*N+2*N,i)

                    YcY                    = sum(abs(ChiMod(:,num_in,1))**2)*DX
                    ChiMod(:,num_in,1)     = ChiMod(:,num_in,1)/sqrt(YcY)

!                    print"(A,i4,A,2f12.6,A)","  K_IN (",num_in,")=",kvec*L2LR,"[nm]"
                else
                    num_out                    = num_out + 1
                    ChiKvec(num_out,2)         = kvec
                    ChiLambda(num_out,2)       = lambda

                    ChiMod(2  :N+1,num_out,2)   = Z(N+1:2*N,i)
                    ChiMod(N+4:2*N+3,num_out,2) = Z(2*N+N+1:2*N+2*N,i)

                    YcY                        = sum(abs(ChiMod(:,num_out,2))**2)*DX
                    ChiMod(:,num_out,2)         = ChiMod(:,num_out,2)/sqrt(YcY)

!                    print"(A,i4,A,2f12.6,A)","  K_OUT(",num_out,")=",kvec*L2LR,"[nm]"
                endif
            endif
        endif ! end of if Beta
    enddo



    call spinmodpop_sort_vectors_by_values(ChiMod(:,1:LICZBA_MODOW,1),ChiKvec(1:LICZBA_MODOW,1),ChiLambda(1:LICZBA_MODOW,1),+1,0)
    call spinmodpop_sort_vectors_by_values(ChiMod(:,1:LICZBA_MODOW,2),ChiKvec(1:LICZBA_MODOW,2),ChiLambda(1:LICZBA_MODOW,2),-1,0)

    if( TRANS_DEBUG ) then
    print*,"-------------------------------------------------------------"
    print*,"                        SPINORY"
    print*," K wave.  :         Input mod.      |         Output mod."
    print*,"-------------------------------------------------------------"

    do i = 1 , LICZBA_MODOW
         iter = 2*(N+2)
         polarA = sum(conjg(ChiMod(1:iter/2,i,1))*ChiMod(1:iter/2,i,1)) - sum(conjg(ChiMod(iter/2+1:iter,i,1))*ChiMod(iter/2+1:iter,i,1))
         polarB = sum(conjg(ChiMod(1:iter/2,i,2))*ChiMod(1:iter/2,i,2)) - sum(conjg(ChiMod(iter/2+1:iter,i,2))*ChiMod(iter/2+1:iter,i,2))
         print"(A,i4,A,3f10.4,A,3f10.4)","K[",i,"][nm]:",ChiKvec(i,1)*L2LR,dble(polarA)*DX," | ",ChiKvec(i,2)*L2LR,dble(polarB)*DX
    enddo

    open(unit = 333, file= "modup+.txt" )
    do i = 1 , 2*(N+2)
    write(333,"(1000f20.6)"),i*DX*LR2L,abs(ChiMod(i,1:LICZBA_MODOW,1))**2
    enddo
    open(unit = 333, file= "modup-.txt" )
    do i = 1 , 2*(N+2)
    write(333,"(1000f20.6)"),i*DX*LR2L,abs(ChiMod(i,1:LICZBA_MODOW,2))**2
    enddo

    print*,""
    print*,"WEKTORY EVANESCENTNE (IN/OUT):"

    endif ! end of if(TRANS_DEBUG)


        do i = 4*N , 1 , -1
            if(abs(Beta(i))>1e-16) then
                lambda= (ALPHA(i)/BETA(i))
                kvec  = (log(lambda)/DX)

                if( abs(lambda) > 1 + 1E-6 ) then
                    num_out                  = num_out + 1
!                    K_m_out(num_out)         = kvec
!                    Chi_m_out(2:N+1,num_out) = Z(N+1:2*N,i)
!                    YcY                      = sum(abs(Chi_m_out(:,num_out))**2)*DX
!                    Chi_m_out(:,num_out)     = Chi_m_out(:,num_out)/sqrt(YcY)

                    ChiKvec(num_out,2)         = kvec
                    ChiLambda(num_out,2)       = lambda

                    ChiMod(2  :N+1,num_out,2)   = Z(N+1:2*N,i)
                    ChiMod(N+4:2*N+3,num_out,2) = Z(2*N+N+1:2*N+2*N,i)

                    YcY                         = sum(abs(ChiMod(:,num_out,2))**2)*DX
                    ChiMod(:,num_out,2)         = ChiMod(:,num_out,2)/sqrt(YcY)

                endif
            endif ! end of if beta
        enddo

        do i = 1 , 4*N
            if(abs(Beta(i))>1e-16) then
                lambda= (ALPHA(i)/BETA(i))
                kvec  = (log(lambda)/DX)
!                print*,i,kvec*L2LR
                if(  abs(lambda) < 1 - 1E-6 ) then

                    num_in                   = num_in + 1
!                    K_M_IN(num_in)           = kvec
!                    Chi_M_IN(2:N+1,num_in)   = Z(N+1:2*N,i)
!                    YcY                      = sum(abs(Chi_M_IN(:,num_in))**2)*DX
!                    Chi_M_IN(:,num_in)       = Chi_M_IN(:,num_in)/sqrt(YcY)


                    ChiKvec(num_in,1)          = kvec
                    ChiLambda(num_in,1)        = lambda
                    ChiMod(2  :N+1,num_in,1)   = Z(N+1:2*N,i)
                    ChiMod(N+4:2*N+3,num_in,1) = Z(2*N+N+1:2*N+2*N,i)

                    YcY                    = sum(abs(ChiMod(:,num_in,1))**2)*DX
                    ChiMod(:,num_in,1)     = ChiMod(:,num_in,1)/sqrt(YcY)


                endif
            endif ! end of if beta
        enddo

    call spinmodpop_sort_vectors_by_values(ChiMod(:,LICZBA_MODOW+1:2*N,1),ChiKvec(LICZBA_MODOW+1:2*N,1),ChiLambda(LICZBA_MODOW+1:2*N,1),+1,1)
    call spinmodpop_sort_vectors_by_values(ChiMod(:,LICZBA_MODOW+1:2*N,2),ChiKvec(LICZBA_MODOW+1:2*N,2),ChiLambda(LICZBA_MODOW+1:2*N,2),-1,1)


!
!        call modpop_sort_vectors_by_values(Chi_M_IN (:,LICZBA_MODOW+1:N),K_M_IN (LICZBA_MODOW+1:N),-1,1)
!        call modpop_sort_vectors_by_values(Chi_M_OUT(:,LICZBA_MODOW+1:N),K_M_OUT(LICZBA_MODOW+1:N),-1,1)
!
        LICZBA_MODOW_EVANESCENTYCH = 0
        ! bierzemy tylko te mody ktore nie maja czesci falowej (tj tylko czyste evanescentne exp(+/-kx))
        do i = LICZBA_MODOW + 1 , 2*N
            !if( abs(imag(K_M_OUT(i))) < 1.0E-6 .and. abs(imag(K_M_IN(i))) < 1.0E-6 .and. &
            ! &  abs(dble(K_M_OUT(i))) > 1.0E-6 .and. abs(dble(K_M_IN(i))) > 1.0E-6  ) then
                LICZBA_MODOW_EVANESCENTYCH = LICZBA_MODOW_EVANESCENTYCH + 1
!                Chi_M_IN (:,LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH) = Chi_M_IN(:,i)
!                K_M_IN   (LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH)   = K_M_IN(i)
!!
!                Chi_M_OUT(:,LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH) = Chi_M_OUT(:,i)
!                K_M_OUT  (LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH)   = K_M_OUT(i)

                !Chi_M_OUT(:,LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH) = Chi_M_IN(:,i)
                !K_M_OUT  (LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH)   = -K_M_IN(i)
                !call modpop_calc_mode_from_k((K_M_OUT(i)),Chi_m_out(:,i),DX,N+2,Ef,BZ,Uvec,pbHorizontal);
           !endif

        enddo
        LICZBA_MODOW_EVANESCENTYCH = LICZBA_MODOW_EVANESCENTYCH/8
        !print*,"NOWAN=",LICZBA_MODOW_EVANESCENTYCH
        if ( TRANS_DEBUG  ) then
        print*,"-------------------------------------------------------------"
        print*," K evan.  :         Input mod.      |         Output mod."
        print*,"-------------------------------------------------------------"

        do i =  LICZBA_MODOW + 1 , LICZBA_MODOW_EVANESCENTYCH + LICZBA_MODOW

             iter = 2*(N+2)
             polarA = sum(conjg(ChiMod(1:iter/2,i,1))*ChiMod(1:iter/2,i,1)) - sum(conjg(ChiMod(iter/2+1:iter,i,1))*ChiMod(iter/2+1:iter,i,1))
             polarB = sum(conjg(ChiMod(1:iter/2,i,2))*ChiMod(1:iter/2,i,2)) - sum(conjg(ChiMod(iter/2+1:iter,i,2))*ChiMod(iter/2+1:iter,i,2))
             print"(A,i4,A,3f10.4,A,3f10.4)","K[",i,"][nm]:",ChiKvec(i,1)*L2LR,dble(polarA)*DX," | ",ChiKvec(i,2)*L2LR,dble(polarB)*DX
        enddo
        print*,"-------------------------------------------------------------"
        endif




    deallocate(Uvec       )
    deallocate(spinuvec   )
    deallocate(Ey         )
    deallocate(GINDEX     )
    deallocate(Hamiltonian)
    deallocate(Mtau )
    deallocate(Mdiag)
    deallocate(Mham )
    deallocate(MatS )
    deallocate(MatBs)
    deallocate(tmp_vec )
    deallocate(tmp_vec2)
    deallocate(MA )
    deallocate(MB )
    deallocate(Z  )

    ! Alokacja macierzy LAPACKA
    deallocate(ALPHA)
    deallocate(BETA)
    deallocate(RWORK)
    deallocate(WORK)

end subroutine spinmodpop_calc_modes_from_wfm

doubleprecision function spinmodpop_calc_current(kvec,chiup,chidown,N) result(rval)
    double precision :: kvec
    complex*16,dimension(:) :: chiup,chidown
    integer :: N ! rozmiar modu dla spinu
    doubleprecision :: rhoup, rhodwn
    complex*16 :: tmpval ,soval  , Xup,Xdwn
    integer :: i
    tmpval = 0
    soval  = 0
    do i = 1 , N
        Xup    = chiup  (i)
        Xdwn   = chidown(i)

        rhoup  = abs(chiup  (i))**2
        rhodwn = abs(chidown(i))**2
        tmpval = tmpval -  rhoup * sin(dx*dx*(i-hny)*Bz-kvec*dx)
        tmpval = tmpval -  rhodwn* sin(dx*dx*(i-hny)*Bz-kvec*dx)
        soval  = soval + ( (rhoup-rhodwn)*so_loc*Ey(i) + II*so_rashba*( Xdwn*conjg(Xup) - conjg(Xdwn)*Xup )  )*dx
    enddo
!    print*,"TD:",kvec,dble(tmpval) ,dble(soval)
    rval = dble(tmpval + soval)

end function spinmodpop_calc_current


doubleprecision function spinmodpop_calc_TB_current(kvec,chiup,chidown,N) result(rval)
    double precision :: kvec
    complex*16     ,dimension(:) :: chiup,chidown

    integer :: N ! rozmiar modu dla spinu
    doubleprecision :: rhoup, rhodwn
    complex*16 :: tmpval ,soval  , Xup,Xdwn , tlj(2,2,2)
    integer :: i
    tmpval = 0
    soval  = 0


    do i = 1 , N
        ! spin up-up
        tlj(1,1,1) = - t0 *  exp(-II*(DX*DX*BZ*(i-hNY))) * exp(+II*DX*Ey(i)*so_loc) ! l+1
        tlj(1,1,2) = - t0 *  exp(+II*(DX*DX*BZ*(i-hNY))) * exp(-II*DX*Ey(i)*so_loc) ! l-1

        ! spin down-down
        tlj(2,2,1) = - t0 *  exp(-II*(DX*DX*BZ*(i-hNY))) * exp(-II*DX*Ey(i)*so_loc) ! l+1
        tlj(2,2,2) = - t0 *  exp(+II*(DX*DX*BZ*(i-hNY))) * exp(+II*DX*Ey(i)*so_loc) ! l-1

        ! spin up-down
        tlj(1,2,1) = +so_rashba*tc0 * exp(-II*(DX*DX*BZ*(i-hNY)))
        tlj(1,2,2) = -so_rashba*tc0 * exp(+II*(DX*DX*BZ*(i-hNY)))! l-1

        ! spin up-down
        tlj(2,1,1) = -so_rashba*tc0 * exp(-II*(DX*DX*BZ*(i-hNY)))
        tlj(2,1,2) = +so_rashba*tc0 * exp(+II*(DX*DX*BZ*(i-hNY)))! l-1
        ! wezly do przodu
        tmpval = tmpval + dx*II*( tlj(1,1,1)*conjg(chiup(i))*YY(kvec,chiup(i),1)     - conjg(tlj(1,1,1)*conjg(chiup(i))*YY(kvec,chiup(i),1)) )
        tmpval = tmpval + dx*II*( tlj(2,2,1)*conjg(chidown(i))*YY(kvec,chidown(i),1) - conjg(tlj(2,2,1)*conjg(chidown(i))*YY(kvec,chidown(i),1)) )

        soval = soval + dx*II*( tlj(1,2,1)*conjg(chiup(i))*YY(kvec,chidown(i),1)   - conjg(tlj(1,2,1)*conjg(chiup(i))*YY(kvec,chidown(i),1)) )
        soval = soval + dx*II*( tlj(2,1,1)*conjg(chidown(i))*YY(kvec,chiup(i),1)   - conjg(tlj(2,1,1)*conjg(chidown(i))*YY(kvec,chiup(i),1)) )

        ! wezly do tylu
!        tmpval = tmpval - dx*II*( tlj(1,1,2)*conjg(chiup(i))*YY(kvec,chiup(i),-1)     - conjg(tlj(1,1,2)*conjg(chiup(i))*YY(kvec,chiup(i),-1)) )
!        tmpval = tmpval - dx*II*( tlj(2,2,2)*conjg(chidown(i))*YY(kvec,chidown(i),-1) - conjg(tlj(2,2,2)*conjg(chidown(i))*YY(kvec,chidown(i),-1)) )
!
!        soval = soval - dx*II*( tlj(1,2,2)*conjg(chiup(i))*YY(kvec,chidown(i),-1)   - conjg(tlj(1,2,2)*conjg(chiup(i))*YY(kvec,chidown(i),-1)) )
!        soval = soval - dx*II*( tlj(2,1,2)*conjg(chidown(i))*YY(kvec,chiup(i),-1)   - conjg(tlj(2,1,2)*conjg(chidown(i))*YY(kvec,chiup(i),-1)) )
    enddo
!    print*,"TB:",kvec,dble(tmpval) ,dble(soval)
    rval = dble(tmpval + soval)
end function spinmodpop_calc_TB_current

complex*16 function YY(kvec,chi,i) result(rval)
    double precision :: kvec
    complex*16 :: chi
    integer :: i

    rval = exp(II*kvec*(i)*dx)*chi

end function YY


subroutine spinmodpop_sort_vectors_by_values(vectors,vals,lambdas,order,takeReal)
    complex*16,dimension(:,:) :: vectors
    complex*16,dimension(:)   :: vals
    complex*16,dimension(:)   :: lambdas
    integer :: order,takeReal
    complex*16,dimension(:),allocatable  :: tmpvec
    complex*16 :: tmpval,tmplambda

    integer :: i , j , n , nvec , imin


    n    = size(vals)
    nvec = size(vectors(:,1))


    allocate(tmpvec(nvec))
!    if(TRANS_DEBUG) then
!     print*,"Sortowanie wektorow:"
!     print*,"  Liczba wektor  :"  ,n
!     print*,"  Rozmiar wektora:",nvec
!    endif


    do i = 1 , n
        imin   = i
        do j = i+1 , n
           if(abs(vals(j)) > 1.0e-10) then
           if(order > 0 ) then
                if(takeReal == 0) then
                    if( imag((vals(imin))) < imag((vals(j))) ) imin = j
                else
                    if( (real(vals(imin))) < (real(vals(j))) ) imin = j
                endif
           endif
           if(order < 0 ) then
                if(takeReal == 0) then
                    if( imag((vals(imin))) > imag((vals(j))) ) imin = j
                else
                    if( (real(vals(imin))) > (real(vals(j))) ) imin = j
                endif
           endif
           endif
        enddo

        tmpval     = vals(i)
        vals(i)    = vals(imin)
        vals(imin) = tmpval

        tmpval        = lambdas(i)
        lambdas(i)    = lambdas(imin)
        lambdas(imin) = tmpval

        tmpvec          = vectors(:,i)
        vectors(:,i)    = vectors(:,imin)
        vectors(:,imin) = tmpvec
    enddo

    deallocate(tmpvec)
end subroutine spinmodpop_sort_vectors_by_values

endmodule spinmodpop
