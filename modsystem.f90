module modsystem
    use modzrodlo
    implicit none
    private

     ! -----------------------------------------------------------------
     ! Zmienne systemowe
     ! -----------------------------------------------------------------
    integer :: no_zrodel      ! liczba zrodel w ukladzie
    integer :: TRANS_MAXN     ! liczba oczek brane do rozwiazania ukladu rownan
    integer :: Nx,Ny          ! wymiar ukladu
    double precision :: DX,Ef,Bz

    type(czrodlo),dimension(:),allocatable        :: zrodla ! tablica z obiektami zrodla
    integer,dimension(:,:), allocatable           :: GFLAGS ! ZAWIERA FLAGI UKLADU (DIRICHLET,WEJSCIA,POZA)
    integer,dimension(:,:), allocatable           :: GINDEX ! INDEKSUJE GEOMETRIE UKLADU (-1) OZNACZA POZA UKLADEM
    double precision,dimension(:,:), allocatable  :: UTOTAL ! MACIERZ POTENCJALU EFEKTYWNEGO


    ! rodzaje zapisu do pliku
    ENUM,BIND(C)
        ENUMERATOR :: ZAPISZ_FLAGI      = 0
        ENUMERATOR :: ZAPISZ_POTENCJAL  = 1
        ENUMERATOR :: ZAPISZ_KONTUR     = 2
        ENUMERATOR :: ZAPISZ_INDEKSY    = 3
    END ENUM
    public :: ZAPISZ_FLAGI , ZAPISZ_POTENCJAL , ZAPISZ_KONTUR , ZAPISZ_INDEKSY
    public :: zrodla,UTOTAL,GFLAGS,GINDEX
    public :: system_inicjalizacja , system_zwalnienie_pamieci , system_inicjalizacja_ukladu
    public :: system_zapisz_do_pliku
    contains


    ! =========================================================
    !               INICJALIZACJA SYSTEMU
    ! =========================================================
    ! PROCEDURA INICJALIZACYJNA WYWOLYWANA JAKO PIERWSZA
    ! PNX,PNY       - ROZMIAR PROBLEMU, WYMIAR MACIERZY GINDEX, GFLAGS
    ! pLiczbaZrodel - LICZBA WSZYTKICH ZRODEL W UKLADZIE
    ! pDX           - Podajemy krok siatki w [nm]
    ! ----------------------------------------------------------
    subroutine system_inicjalizacja(pnx,pny,pLiczbaZrodel,pdx)
        integer,intent(in)             :: pnx,pny
        integer,intent(in)             :: pLiczbaZrodel
        double precision, intent(in)   :: pdx

        if(TRANS_DEBUG==.true.) then
            print*,"System: Inicjalizacja:";
            print*,"    nx          =",pnx
            print*,"    ny          =",pny
            print*,"    liczbaZrodel=",pLiczbaZrodel
        endif
        nx = pnx
        ny = pny
        dx = pdx*L2LR

        ! alokacja pamieci
        allocate( GFLAGS(nx,ny))
        allocate( GINDEX(nx,ny))
        allocate( UTOTAL(nx,ny))
        GFLAGS = B_NORMAL
        GINDEX = 0
        UTOTAL = 0
        ! Tworzymy uklad...
        no_zrodel = pLiczbaZrodel
        allocate(zrodla(no_zrodel))
        ! okreslanie rzeczywistego rozmiaru problemu, punkty o indeksach -1 leza poza ukladem, nie beda
        ! potrzebne wiec zatem do obliczen.
        TRANS_MAXN = 0;
    end subroutine system_inicjalizacja


    ! =========================================================
    !
    !           CZYSZCZENIE ZAALOKOWANEJ PAMIECI
    !       funkcja wywolywana jest na samym koncu po wykonaniu
    !       obliczen
    ! =========================================================

    subroutine system_zwalnienie_pamieci()
        integer :: i
        if(TRANS_DEBUG==.true.) print*,"System: Zwalnianie pamieci"
        if(allocated(GFLAGS))   deallocate(GFLAGS)
        if(allocated(GINDEX))   deallocate(GINDEX)
        if(allocated(UTOTAL))   deallocate(UTOTAL)
        do i = 1 , no_zrodel
            call zrodla(i)%zrodlo_zwolnij_pamiec()
        enddo
        if(allocated(zrodla)) deallocate(zrodla)
    end subroutine system_zwalnienie_pamieci

    ! ------------------------------------------------------------
    ! Funkcja tworzy tablice flag oraz gindex na podstawie
    ! zmodyfikowanej w mainie tablicy gflags
    ! Podajemy:
    ! in_len - dlugosc na jaka od wejsc zostanie utworzona flaga B_NORMAL
    ! smooth_rad - promien wygladzania ostrych katow w ukladzie
    ! smooth_iter - liczba iteracji wygladzania
    ! ------------------------------------------------------------
    subroutine system_inicjalizacja_ukladu(in_len,smooth_rad,smooth_iter)
        integer ,intent(in) :: in_len,smooth_rad,smooth_iter
        integer :: soft_iter , iii , jjj , norm_iter , empty_iter
        integer :: nrz , i , j , ni, nj , mi , mj , iter

        print * , "System: Tworzenie tablic GINDEX oraz GFLAGS"
        ! dodawanie zrodej wejsciowych do warunkow brzegowych
        do nrz = 1 , no_zrodel
            do i = 2 , zrodla(nrz)%N - 1
                ni = zrodla(nrz)%polozenia(i,1)
                nj = zrodla(nrz)%polozenia(i,2)
                GFLAGS(ni,nj) = B_WEJSCIE
            enddo
            if(zrodla(nrz)%bKierunek)then
                ni = zrodla(nrz)%polozenia(1,1) + 1
                nj = zrodla(nrz)%polozenia(1,2)
                mi = ni + in_len - 1
                mj = zrodla(nrz)%polozenia(zrodla(nrz)%N,2)
                GFLAGS(ni:mi,nj:mj) = B_NORMAL
                print*,ni,mi,nj,mj
            else
                ni = zrodla(nrz)%polozenia(1,1) - 1
                nj = zrodla(nrz)%polozenia(1,2)
                mi = ni - in_len + 1
                mj = zrodla(nrz)%polozenia(zrodla(nrz)%N,2)
                GFLAGS(mi:ni,nj:mj) = B_NORMAL
                print*,ni,mi,nj,mj
            endif
        enddo

        ! -------------------------------------------------------------
		! wygladzanie krawedzi
        ! -------------------------------------------------------------
		do soft_iter = 1 , smooth_iter
        do i =  1+smooth_rad*2 , nx-smooth_rad*2-1
        do j =  1+smooth_rad , ny-smooth_rad-1
         norm_iter  = 0
		 empty_iter = 0
			do iii =  -smooth_rad , smooth_rad , 1
			do jjj =  -smooth_rad , smooth_rad , 1
				if(GFLAGS(i+iii,j+jjj) == B_EMPTY )  empty_iter = empty_iter + 1
				if(GFLAGS(i+iii,j+jjj) == B_NORMAL ) norm_iter = norm_iter   + 1
			enddo
			enddo
			if(empty_iter < norm_iter .and. GFLAGS(i,j) == B_EMPTY ) GFLAGS(i,j) = B_NORMAL
        enddo
        enddo
		enddo ! end of iter

		! tworzenie odpowiednich warunkow brzegowych (tak zeby nie bylo poloczenia
		! B_normal z B_Empty)
        do i = 2 , nx-1
        do j = 2 , ny-1
        		if (GFLAGS(i,j) == B_NORMAL) then
        			iter = (GFLAGS(i-1,j)-B_EMPTY) * (GFLAGS(i+1,j)-B_EMPTY) * (GFLAGS(i,j-1)-B_EMPTY) * (GFLAGS(i,j+1)-B_EMPTY)
        			if ( iter == 0 ) GFLAGS(i,j) = B_DIRICHLET
                endif
        enddo
        enddo

        ! usuwanie zbednych flag dirichleta
        do i = 1 , nx
        do j = 1 , ny
        		if (GFLAGS(i,j) == B_DIRICHLET) then
                    iter = 1
                    if(i>1 ) iter = iter*GFLAGS(i-1,j)*(GFLAGS(i-1,j)-B_WEJSCIE)
                    if(j>1 ) iter = iter*GFLAGS(i,j-1)*(GFLAGS(i,j-1)-B_WEJSCIE)
                    if(i<NX) iter = iter*GFLAGS(i+1,j)*(GFLAGS(i+1,j)-B_WEJSCIE)
                    if(j<NY) iter = iter*GFLAGS(i,j+1)*(GFLAGS(i,j+1)-B_WEJSCIE)
        			if ( iter /= 0 ) GFLAGS(i,j) = B_EMPTY
                endif
        enddo
        enddo

        ! nadawanie indeksow
        GINDEX = 0
        iter   = 1
        do i = 1 , nx
        do j = 1 , ny
           if( GFLAGS(i,j) == B_WEJSCIE .or. GFLAGS(i,j) == B_NORMAL .or. GFLAGS(i,j) == B_DIRICHLET  ) then
           GINDEX(i,j) = iter
           iter = iter + 1
           endif
        enddo
        enddo

    end subroutine system_inicjalizacja_ukladu

    ! ------------------------------------------------------------ -------
    ! Funkcja zapisuje do pliku nazwa dane wskazane przez typ (patrz enum)
    ! ------------------------------------------------------------ -------
    subroutine system_zapisz_do_pliku(nazwa,typ)
    character(*) , intent(in) :: nazwa
    integer                   :: typ

    double precision :: fval
    integer          :: i,j

    open(unit=86554,file=nazwa)
    select case(typ)
    ! ---------------------------------------
    case(ZAPISZ_KONTUR)
        do i = 1 , NX
        do j = 1 , NY
                fval = 0;
                if(GINDEX(i,j) > 0) fval = 1
                write(86554,"(3e20.6)"),i*DX*Lr2L,j*DX*Lr2L,fval
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case(ZAPISZ_FLAGI)
        do i = 1 , NX
        do j = 1 , NY
                fval = GFLAGS(i,j)
                write(86554,"(3e20.6)"),i*DX*Lr2L,j*DX*Lr2L,fval
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case(ZAPISZ_POTENCJAL)
        do i = 1 , NX
        do j = 1 , NY
                fval = UTOTAL(i,j)
                write(86554,"(3e20.6)"),i*DX*Lr2L,j*DX*Lr2L,fval
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case(ZAPISZ_INDEKSY)
        do i = 1 , NX
        do j = 1 , NY
                fval = GINDEX(i,j)
                write(86554,"(3e20.6)"),i*DX*Lr2L,j*DX*Lr2L,fval
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case default
            print*,"System: Zapisz do pliku - podano zly argument."
    end select

    close(86554)
    endsubroutine system_zapisz_do_pliku


end module
