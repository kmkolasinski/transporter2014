module modsystem
    use modzrodlo
    implicit none
    private

     ! -----------------------------------------------------------------
     ! Zmienne systemowe
     ! -----------------------------------------------------------------
    integer :: no_zrodel      ! liczba zrodel w ukladzie
    integer :: no_abs_zrodel  ! libcza wirtualnych zrodel (w.b. Nowaka)
    integer :: TRANS_MAXN     ! liczba oczek brane do rozwiazania ukladu rownan
    integer :: MATASIZE
    integer :: Nx,Ny          ! wymiar ukladu
    double precision :: DX,Ef,Bz
    doubleprecision :: TRANS_R ! prawdopodobienstwo odbicia
    doubleprecision :: TRANS_T ! prawdopodobienstwo przejscia


    type(czrodlo),dimension(:),allocatable        :: zrodla    ! TABLICA Z OBIEKTAMI ZRODLA
    type(cabs_zrodlo),dimension(:),allocatable     :: abs_zrodla ! TABLICA Z OBIEKTAMI WIRTUALNYCH W.B.
    integer,dimension(:,:), allocatable           :: GFLAGS ! ZAWIERA FLAGI UKLADU (DIRICHLET,WEJSCIA,POZA)
    integer,dimension(:,:), allocatable           :: ZFLAGS ! ZAWIERA FLAGI WEJSC WARTOSC FLAGI OZNACZA NUMER WEJSCIA - 1,2, ...
    integer,dimension(:,:), allocatable           :: GINDEX ! INDEKSUJE GEOMETRIE UKLADU (-1) OZNACZA POZA UKLADEM
    double precision,dimension(:,:), allocatable  :: UTOTAL ! MACIERZ POTENCJALU EFEKTYWNEGO
    complex*16,dimension(:),allocatable           :: CMATA  ! GLOWNA MACIERZ PROGRAMU W FORMACIE (ROW,COL,VALS), TUTAJ TYLKO VALS
    integer,dimension(:,:),allocatable            :: IDXA   ! INDEKSY MACIERZY (ROW,COL)
    complex*16,dimension(:),allocatable           :: VPHI   ! SZUKANA FUNKCJA FALOWA, PRZELICZANA POTEM NA LDOS
    complex*16,dimension(:,:),allocatable         ::  PHI   ! FUNKCJA FALOWA ZAPISANA NA DWUWYMIAROWEJ MACIERZY
    double precision,dimension(:,:), allocatable  :: CURRENT! FUNKCJA GESTOSCI PRADU PRAWDOPODOBIENSTWA

    ! rodzaje zapisu do pliku
    ENUM,BIND(C)
        ENUMERATOR :: ZAPISZ_FLAGI      = 0
        ENUMERATOR :: ZAPISZ_POTENCJAL  = 1
        ENUMERATOR :: ZAPISZ_KONTUR     = 2
        ENUMERATOR :: ZAPISZ_INDEKSY    = 3
        ENUMERATOR :: ZAPISZ_PHI        = 4
        ENUMERATOR :: ZAPISZ_J          = 5
    END ENUM

    public :: ZAPISZ_FLAGI , ZAPISZ_POTENCJAL , ZAPISZ_KONTUR , ZAPISZ_INDEKSY , ZAPISZ_PHI , ZAPISZ_J
    public :: zrodla,UTOTAL,GFLAGS,GINDEX
    public :: system_inicjalizacja , system_zwalnienie_pamieci , system_inicjalizacja_ukladu
    public :: system_dodaj_abs_zrodlo !(pY1,pYN,pX1,pEf,pKierunek)
    public :: system_zapisz_do_pliku , system_rozwiaz_problem
    public :: system_gauss , system_dodaj_lorentza , system_dodaj_pionowy_slupek_potencjalu
    public :: system_fermi, system_fermiT , system_dfermidE
    public :: TRANS_T,TRANS_R
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
        allocate( ZFLAGS(nx,ny))
        allocate( GINDEX(nx,ny))
        allocate( UTOTAL(nx,ny))
        allocate(    PHI(nx,ny))
        allocate(CURRENT(nx,ny))
        GFLAGS = B_NORMAL
        ZFLAGS = 0
        GINDEX = 0
        UTOTAL = 0
        ! Tworzymy uklad...
        no_zrodel = pLiczbaZrodel
        allocate(zrodla(no_zrodel))
        ! okreslanie rzeczywistego rozmiaru problemu, punkty o indeksach -1 leza poza ukladem, nie beda
        ! potrzebne wiec zatem do obliczen.
        TRANS_MAXN = 0;
    end subroutine system_inicjalizacja

    ! ========================================================================================
    !               FUNKCJA ALOKUJE NOWE ZRODLO Z
    !   TRANSPARENTNYMI WARUNKAMI BRZEGOWYMI. ALOKACJA PRZEPRO-
    ! WADZANA JEST DYNAMICZNIE, DLATEGO NIE TRZEBA PODAWAC INDEKSOW.
    ! ========================================================================================
    ! Ustawiamy wczesniej zaalokowane zrodlo przez polecenie systemowe.
    ! Parametry:
    ! pY1,pYN - polozenia y na siatce numerycznej liczone od 1 do NY, w przypadku zrodla poziomego
    !           oznaczaja polozenia X1 oraz XN
    ! pX1 - polozenie X zrodla, dla zrodla poziomego polozenie Y1
    ! pEf -Ef [meV]
    ! pKierunek - enum ZRODLO_KIERUNEK_PRAWO/LEWO/GORA/DOL - ustala w ktora skierowane jest zrodlo
    ! --------------------------------------------------------------------
    subroutine system_dodaj_abs_zrodlo(pY1,pYN,pX1,pEf,pKierunek)
        integer,intent(in)         ::  pY1,pYN,pX1
        doubleprecision,intent(in) ::  pEf
        integer,intent(in)         ::  pKierunek ! enum
        type(cabs_zrodlo),dimension(:),allocatable     :: tmp_abs_zrodla
        integer :: i
        no_abs_zrodel = size(abs_zrodla)


        allocate(tmp_abs_zrodla(no_abs_zrodel+1))
        do i = 1 , no_abs_zrodel
            call tmp_abs_zrodla(i)%abs_zrodlo_skopiuj(abs_zrodla(i))
            call abs_zrodla(i)%abs_zrodlo_zwolnij_pamiec()
        enddo

        if(allocated(abs_zrodla)) deallocate(abs_zrodla)
        allocate(abs_zrodla(no_abs_zrodel+1))
        do i = 1 , no_abs_zrodel
            call abs_zrodla(i)%abs_zrodlo_skopiuj(tmp_abs_zrodla(i))
        enddo

        call abs_zrodla(no_abs_zrodel+1)%abs_zrodlo_ustaw(pY1,pYN,pX1,pEf,pKierunek)
        no_abs_zrodel = no_abs_zrodel + 1
    end subroutine system_dodaj_abs_zrodlo



    ! =========================================================
    !
    !           CZYSZCZENIE ZAALOKOWANEJ PAMIECI
    !       funkcja wywolywana jest na samym koncu po wykonaniu
    !       obliczen
    ! =========================================================

    subroutine system_zwalnienie_pamieci()
        integer :: i
        if(TRANS_DEBUG==.true.) print*,"System: Zwalnianie pamieci"
        print*,"    Czyszczenie tablic..."
        if(allocated(GFLAGS))   deallocate(GFLAGS)
        if(allocated(ZFLAGS))   deallocate(ZFLAGS)
        if(allocated(GINDEX))   deallocate(GINDEX)
        if(allocated(UTOTAL))   deallocate(UTOTAL)
        if(allocated(CURRENT))  deallocate(CURRENT)
        if(allocated(PHI))      deallocate(PHI)
        print*,"    Czyszczenie zrodel..."
        do i = 1 , no_zrodel
            call zrodla(i)%zrodlo_zwolnij_pamiec()
        enddo
        print*,"    Czyszczenie abs_zrodel..."
        do i = 1 , no_abs_zrodel
            call abs_zrodla(i)%abs_zrodlo_zwolnij_pamiec()
        enddo
        if(allocated(zrodla))    deallocate(zrodla)
        if(allocated(abs_zrodla)) deallocate(abs_zrodla)
    end subroutine system_zwalnienie_pamieci

    ! ------------------------------------------------------------
    ! Funkcja tworzy tablice flag oraz gindex na podstawie
    ! zmodyfikowanej w mainie tablicy gflags
    ! Podajemy:
    ! in_len - dlugosc na jaka od wejsc zostanie utworzona flaga B_NORMAL
    ! smooth_rad - promien wygladzania ostrych katow w ukladzie
    ! smooth_iter - liczba iteracji wygladzania
    ! ------------------------------------------------------------caliente
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
                ZFLAGS(ni,nj) = nrz ! przypisujemy fladze numer wejscia
            enddo


            select case (zrodla(nrz)%bKierunek)
            case (ZRODLO_KIERUNEK_PRAWO)
                ni = zrodla(nrz)%polozenia(1,1) + 1
                nj = zrodla(nrz)%polozenia(1,2)
                mi = ni + in_len - 1
                mj = zrodla(nrz)%polozenia(zrodla(nrz)%N,2)

                GFLAGS(ni:mi,nj:mj) = B_NORMAL
            case (ZRODLO_KIERUNEK_LEWO)
                ni = zrodla(nrz)%polozenia(1,1) - 1
                nj = zrodla(nrz)%polozenia(1,2)
                mi = ni - in_len + 1
                mj = zrodla(nrz)%polozenia(zrodla(nrz)%N,2)

                GFLAGS(mi:ni,nj:mj) = B_NORMAL
            case (ZRODLO_KIERUNEK_GORA)
                ni = zrodla(nrz)%polozenia(1,1)
                nj = zrodla(nrz)%polozenia(1,2) + 1
                mi = zrodla(nrz)%polozenia(zrodla(nrz)%N,1)
                mj = nj + in_len - 1

                GFLAGS(ni:mi,nj:mj) = B_NORMAL
            case (ZRODLO_KIERUNEK_DOL)
                ni = zrodla(nrz)%polozenia(1,1)
                nj = zrodla(nrz)%polozenia(1,2) - 1
                mi = zrodla(nrz)%polozenia(zrodla(nrz)%N,1)
                mj = nj - in_len + 1

                GFLAGS(ni:mi,mj:nj) = B_NORMAL
            endselect


        enddo
        ! przypisywanie  transparentnych w.b.
        do nrz = 1 , no_abs_zrodel
            do i = 2 , abs_zrodla(nrz)%N - 1
                ni = abs_zrodla(nrz)%polozenia(i,1)
                nj = abs_zrodla(nrz)%polozenia(i,2)
                GFLAGS(ni,nj) = B_TRANSPARENT
                ZFLAGS(ni,nj) = nrz ! przypisujemy fladze numer wejscia
            enddo
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
           if( GFLAGS(i,j)  == B_WEJSCIE      .or. &
              & GFLAGS(i,j) == B_NORMAL       .or. &
              & GFLAGS(i,j) == B_TRANSPARENT  .or. &
              & GFLAGS(i,j) == B_DIRICHLET  ) then
           GINDEX(i,j) = iter
           iter = iter + 1
           endif
        enddo
        enddo

    end subroutine system_inicjalizacja_ukladu


    ! --------------------------------------------------------------------
    ! Rozwiazywanie problemu dla zrodla o podanym numerze - nrz
    ! --------------------------------------------------------------------
    subroutine system_rozwiaz_problem(nrz,TR_MAT)
        integer,intent(in)  :: nrz
        double precision,dimension(:,:), allocatable  :: TR_MAT
        integer,allocatable :: HBROWS(:)
        integer :: i,j,ni,nj,mod_in

        Ef  = atomic_Ef/1000.0/Rd
        BZ  = BtoDonorB(atomic_Bz)
        print*,"! ----------------------------------------------- !"
        print*,"! Rozpoczecie obliczen dla zrodla=", nrz
        print*,"!   Ef=", Ef*Rd*1000    , "[meV]"
        print*,"!   Bz=", DonorBtoB(Bz) , "[T]"
        print*,"! ----------------------------------------------- !"



        call oblicz_rozmiar_macierzy()

        allocate(CMATA(MATASIZE))
        allocate(IDXA (MATASIZE,2))
        allocate(VPHI(TRANS_MAXN))
        allocate(HBROWS(TRANS_MAXN+1))

        ! tablice z wspolczynnikami T i R
        if(allocated(TR_MAT))deallocate(TR_MAT)
        allocate(TR_MAT(no_zrodel,zrodla(nrz)%liczba_modow))

        TR_MAT = 0

        CURRENT = 0;
        PHI     = 0

        TRANS_R = 0
        TRANS_T = 0
        ! --------------------------------------------------------------------
        !
        ! --------------------------------------------------------------------
        call wypelnij_macierz()
        call convert_to_HB(MATASIZE,IDXA,HBROWS)

        do mod_in = 1 , zrodla(nrz)%liczba_modow

            zrodla(nrz)%ck(:)      = 0
            zrodla(nrz)%ck(mod_in) = 1

            call zrodla(nrz)%zrodlo_oblicz_Fj(dx)
            VPHI = 0
            do i = 2 , zrodla(nrz)%N - 1
                ni = zrodla(nrz)%polozenia(i,1)
                nj = zrodla(nrz)%polozenia(i,2)
                VPHI(GINDEX(ni,nj)) = zrodla(nrz)%Fj(i)
            enddo


            call solve_system(TRANS_MAXN,MATASIZE,IDXA(:,2),HBROWS,CMATA(:),VPHI)
            call oblicz_TR(nrz,mod_in)
            call system_oblicz_J()

            do i = 1 ,no_zrodel
                TR_MAT(i,mod_in) = sum(abs(zrodla(i)%Jout(:)) )
            enddo


            do i = 1 , Nx
            do j = 1 , Ny
                if(GINDEX(i,j) > 0) PHI(i,j) = PHI(i,j) + abs( VPHI(GINDEX(i,j)) )**2
            enddo
            enddo

        enddo ! end of petla po modach

        deallocate(CMATA)
        deallocate(IDXA)
        deallocate(VPHI)
        deallocate(HBROWS)

    end subroutine system_rozwiaz_problem

    ! ---------------------------------------------------------------------
    !   Wyliczanie rozmiaru macierzy: zakladamy, ze na :
    ! 1. warunek Dirichleta potrzeba jedna wartosc
    ! 2. normalny punkt wewnatrz ukladu potrzeba 5 wartosci (4 siasiadow + 1 wezel centralny)
    ! 3. dla wejsc potrzeba N-2 + 1 wezlow: gdzie N-2 oznacza ze bierzemy
    !    wezly wejsciowe od 2 .. N -1, a +1 wynika z faktu, ze mamy zawsze
    !    polaczenie z wezlem (i+1,j), lub (i-1,j)
    ! ---------------------------------------------------------------------

    subroutine oblicz_rozmiar_macierzy()
        integer :: i , j
        MATASIZE   = 0
        TRANS_MAXN = 0
        do i = 1 , NX
        do j = 1 , Ny
            select case(GFLAGS(i,j))
            case(B_DIRICHLET)
                    MATASIZE   = MATASIZE + 1
                    TRANS_MAXN = TRANS_MAXN + 1
            case(B_WEJSCIE)
                    MATASIZE   = MATASIZE + zrodla(ZFLAGS(i,j))%N + 1 - 2 ! +1 na relacje prawo/lewo/gora/dol -2 na obciecie indeksow (1 i N)
                    TRANS_MAXN = TRANS_MAXN + 1
            case(B_TRANSPARENT)
                    MATASIZE   = MATASIZE + 2 ! na transparentne w.b. potrzebujemy dwie komurki
                    TRANS_MAXN = TRANS_MAXN + 1
            case(B_NORMAL)
                    MATASIZE = MATASIZE + 5
                    TRANS_MAXN = TRANS_MAXN + 1
            end select
        enddo
        enddo
        print*,"Rozmiar macierzy - H:",MATASIZE
        print*,"Rozmiar problemu - Y:",TRANS_MAXN
    end subroutine oblicz_rozmiar_macierzy


    ! ---------------------------------------------------------------------
    !   Funkcja wypelnia macierz hamiltonianu w sposob rzadki, zapisujac
    ! kolejne niezerowe elementy H do dwoch wektorow:
    ! cmatA(i) - i-ta wartosc niezerowa w macierzy
    ! idxA(i,1) - w jakim wierszu macierzy sie znajduje i-ty nie zerowy element
    ! idxA(i,2) - oraz w jakiej kolumnie
    ! W zaleznosci od flagi wezla elementy macierzy H sa wyznaczane w
    ! odpowiedni sposob.
    ! ---------------------------------------------------------------------

    subroutine wypelnij_macierz()

        ! zmienne pomocniczne
        integer          :: i,j,itmp,ni,nj,pnj,nn,ln,nzrd,pni
        complex*16       :: post
        doubleprecision  :: kvec
        itmp  = 1
        do i = 1, nx
        do j = 1, ny
            if(GINDEX(i,j) > 0) then
                if( GFLAGS(i,j) == B_DIRICHLET ) then

                    cmatA(itmp)  = CMPLX(1.0)
                    idxA(itmp,1) = GINDEX(i, j)
                    idxA(itmp,2) = GINDEX(i, j)
                    itmp = itmp + 1

                elseif( GFLAGS(i,j) == B_NORMAL) then
                    cmatA(itmp)   = CMPLX( 2.0/DX/DX + UTOTAL(i,j) - Ef )
                    idxA (itmp,1) = GINDEX(i,j)
                    idxA (itmp,2) = GINDEX(i,j)
                    itmp = itmp + 1

                    cmatA(itmp)   = CMPLX(-0.5/DX/DX*EXP(+II*DX*DX*(j)*BZ)  )
                    idxA (itmp,1) = GINDEX(i  ,j)
                    idxA (itmp,2) = GINDEX(i-1,j)
                    itmp = itmp + 1

                    cmatA(itmp)   = CMPLX(-0.5/DX/DX*EXP(-II*DX*DX*(j)*BZ)  )
                    idxA (itmp,1) = GINDEX(i  ,j)
                    idxA (itmp,2) = GINDEX(i+1,j)
                    itmp = itmp + 1

                    cmatA(itmp) = CMPLX(-0.5/DX/DX)
                    idxA(itmp,1) = GINDEX(i,j)
                    idxA(itmp,2) = GINDEX(i,j+1)
                    itmp = itmp + 1

                    cmatA(itmp) = CMPLX(-0.5/DX/DX)
                    idxA(itmp,1) = GINDEX(i,j)
                    idxA(itmp,2) = GINDEX(i,j-1)
                    itmp = itmp + 1

                ! ----------------------------------------------------------------------
                ! Obsluga wejsc
                ! ----------------------------------------------------------------------
                else if( GFLAGS(i,j) == B_WEJSCIE) then
                    nzrd = ZFLAGS(i,j)



                    ni   = i
                    nj   = j ! globalne polozenia na siatce

                    select case (zrodla(nzrd)%bKierunek)
                    ! ------------------------------------------------------------------
                    ! Zrodla prawe:
                    !
                    !                      ------------------------>
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_PRAWO)
                        ln   = j - zrodla(nzrd)%polozenia(1,2) + 1 ! lokalny indeks


                        cmatA(itmp)   = CMPLX(-0.5/DX/DX*2*cos(DX*DX*(nj)*BZ))
                        idxA (itmp,1) = GINDEX(i  ,j)
                        idxA (itmp,2) = GINDEX(i+1,j)
                        itmp = itmp + 1
                        post = 0.5/DX/DX*(EXP(+II*DX*DX*(nj)*BZ)) ! zmienna pomocnicza
                    ! ------------------------------------------------------------------
                    ! Zrodla lewe:
                    !
                    !                      <------------------------
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_LEWO)

                        ln   = j - zrodla(nzrd)%polozenia(1,2) + 1

                        cmatA(itmp)   = CMPLX(-0.5/DX/DX*2*cos(DX*DX*(nj)*BZ))
                        idxA (itmp,1) = GINDEX(i  ,j)
                        idxA (itmp,2) = GINDEX(i-1,j)
                        itmp = itmp + 1
                        post =-0.5/DX/DX*(EXP(-II*DX*DX*(nj)*BZ))
                    ! ------------------------------------------------------------------
                    ! Zrodla dolne:
                    !
                    !                      ------------------------>
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_GORA)

                        ln   = i - zrodla(nzrd)%polozenia(1,1) + 1

                        cmatA(itmp)   = CMPLX(-0.5/DX/DX*2)
                        idxA (itmp,1) = GINDEX(i  ,j)
                        idxA (itmp,2) = GINDEX(i,j+1)
                        itmp = itmp + 1

                        post = 0.5/DX/DX
                    ! ------------------------------------------------------------------
                    ! Zrodla gorne:
                    !
                    !                      <------------------------
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_DOL)

                        ln   = i - zrodla(nzrd)%polozenia(1,1) + 1

                        cmatA(itmp)   = CMPLX(-0.5/DX/DX*2)
                        idxA (itmp,1) = GINDEX(i  ,j)
                        idxA (itmp,2) = GINDEX(i,j-1)
                        itmp = itmp + 1

                        post = -0.5/DX/DX
                    endselect

                    select case (zrodla(nzrd)%bKierunek)
                    ! ------------------------------------------------------------------
                    !
                    !                      <------------------------
                    !                      ------------------------>
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_PRAWO,ZRODLO_KIERUNEK_LEWO)
                    do nn = 2 , zrodla(nzrd)%N - 1
                    pnj   = zrodla(nzrd)%polozenia(nn,2)
                    if( ln == nn  ) then

                        cmatA(itmp)   = CMPLX( 2.0/DX/DX + UTOTAL(i,j) - Ef ) &
                           & + post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
                        idxA (itmp,1) = GINDEX(ni,nj)
                        idxA (itmp,2) = GINDEX(ni,pnj)
                        itmp = itmp + 1

                    else if( nn == ln+1  ) then
                        cmatA(itmp) = CMPLX(-0.5/DX/DX) &
                        &   + post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
                        idxA(itmp,1) = GINDEX(ni,nj)
                        idxA(itmp,2) = GINDEX(ni,pnj)
                        itmp = itmp + 1
                    else if( nn == ln-1  ) then

                        cmatA(itmp) = CMPLX(-0.5/DX/DX) &
                        &   + post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
                        idxA(itmp,1) = GINDEX(ni,nj)
                        idxA(itmp,2) = GINDEX(ni,pnj)
                        itmp = itmp + 1
                    else
                        cmatA(itmp)  = post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
                        idxA(itmp,1) = GINDEX(ni,nj)
                        idxA(itmp,2) = GINDEX(ni,pnj)
                        itmp = itmp + 1
                    endif
                    enddo
                    ! ------------------------------------------------------------------
                    !
                    !                      <------------------------
                    !                      ------------------------>
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_GORA,ZRODLO_KIERUNEK_DOL)

                    do nn = 2 , zrodla(nzrd)%N - 1
                    pni   = zrodla(nzrd)%polozenia(nn,1)

                    if( ln == nn  ) then

                        cmatA(itmp)   = CMPLX( 2.0/DX/DX + UTOTAL(i,j) - Ef ) &
                           & + post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
                        idxA (itmp,1) = GINDEX(ni,nj)
                        idxA (itmp,2) = GINDEX(pni,nj)
                        itmp = itmp + 1

                    else if( nn == ln+1  ) then
                        cmatA(itmp) = CMPLX(-0.5/DX/DX)*(EXP(-II*DX*DX*(nj)*BZ)) &
                        &   + post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
                        idxA(itmp,1) = GINDEX(ni,nj)
                        idxA(itmp,2) = GINDEX(pni,nj)
                        itmp = itmp + 1
                    else if( nn == ln-1  ) then

                        cmatA(itmp) = CMPLX(-0.5/DX/DX)*(EXP(+II*DX*DX*(nj)*BZ)) &
                        &   + post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
                        idxA(itmp,1) = GINDEX(ni,nj)
                        idxA(itmp,2) = GINDEX(pni,nj)
                        itmp = itmp + 1
                    else
                        cmatA(itmp)  = post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
                        idxA(itmp,1) = GINDEX(ni,nj)
                        idxA(itmp,2) = GINDEX(pni,nj)
                        itmp = itmp + 1
                    endif

                    enddo
                    endselect

                ! ----------------------------------------------------------------------- !
                ! endif ! end of if WEJSCIE, transparentne warunki brzegowe
                ! ----------------------------------------------------------------------- !
                else if( GFLAGS(i,j) == B_TRANSPARENT) then

                    nzrd = ZFLAGS(i,j)
                    kvec = abs_zrodla(nzrd)%kvec
                    ni   = i
                    nj   = j ! globalne polozenia na siatce

                    cmatA(itmp)   = 1
                    idxA (itmp,1) = GINDEX(i, j)
                    idxA (itmp,2) = GINDEX(i, j)
                    itmp = itmp + 1

                    cmatA(itmp)   =-EXP(II*KVEC*DX)
                    idxA (itmp,1) = GINDEX(i, j)
                    select case (abs_zrodla(nzrd)%bKierunek)
                    ! ------------------------>
                    case (ZRODLO_KIERUNEK_PRAWO)
                        idxA (itmp,2) = GINDEX(i+1,j)
                    !<------------------------
                    case (ZRODLO_KIERUNEK_LEWO)
                        idxA (itmp,2) = GINDEX(i-1,j)
                    ! ------------------------>
                    case (ZRODLO_KIERUNEK_GORA)
                        idxA (itmp,2) = GINDEX(i, j+1)
                    ! <------------------------
                    case (ZRODLO_KIERUNEK_DOL)
                        idxA (itmp,2) = GINDEX(i, j-1)
                    endselect
                    itmp = itmp + 1

                    endif ! rodzaj komorki - ostatni else B_TRANSPARENT
            endif ! end if GINDEX > 0
        enddo ! end of j
        enddo ! end of i

        if ( itmp -1 - MATASIZE /= 0 ) then
            print *, "Blad rozmiar maicierzy niezgodny z przewidywanym:" , itmp - 1
            stop
        else
            print *, "Macierz skonstruowana poprawnie..."
        endif


    end subroutine wypelnij_macierz

    ! -------------------------------------------------------------------
    ! Procedura oblicza prawdopodobienstwo przejscia T oraz odbicia R
    ! dla podanego zrodla wejsciowego nrz i danego modu wchodzacego mod_in
    ! -------------------------------------------------------------------
    subroutine oblicz_TR(nrz,mod_in)
        integer,intent(in) :: nrz,mod_in
        integer :: i
        doubleprecision :: Jin
        ! obliczamy apliduty prawdopodobienstwa
        do i = 1 , no_zrodel
            call zrodla(i)%zrodlo_oblicz_dk(VPHI,GINDEX,dx)
            !call zrodla(i)%zrodlo_wypisz_ckdk()
        enddo
        ! na podstawie amplitud prawdopodobienstwa obliczamy strumienie
        do i = 1 , no_zrodel
            call zrodla(i)%zrodlo_oblicz_JinJout(dx)
        enddo

        ! na podstawie strumieni obliczamy transmisje oraz prawdopodobienstwo odbicia
        Jin = zrodla(nrz)%Jin(mod_in)
        do i = 1 , no_zrodel
            zrodla(i)%Jin(:)  = zrodla(i)%Jin(:)   / Jin
            zrodla(i)%Jout(:) = zrodla(i)%Jout(:)  / Jin
            call zrodla(i)%zrodlo_wypisz_JinJout()
        enddo

        TRANS_R = TRANS_R + sum( abs(zrodla(nrz)%Jout(:)) )
        !TRANS_T = 0
        do i = 1 ,no_zrodel
            if(i /= nrz) TRANS_T = TRANS_T + sum( abs(zrodla(i)%Jout(:)) )
        enddo
        print*,"T = ", TRANS_T
        print*,"R = ", TRANS_R
        print*,"W = ", TRANS_R + TRANS_T

    end subroutine oblicz_TR
    ! -------------------------------------------------------------------
    ! Zwraca wartosc rozkladu gaussa w punkcie (x,y) dla gaussa o srodku
    ! w (xpos,ypos) oraz sigmie = sigma i amplidudzie = amplitude.
    ! x,y,xpos,ypos - podajemy w indeksach siatki (i,j)
    ! -------------------------------------------------------------------
    doubleprecision function system_gauss(x,y,xpos,ypos,sigma,amplitude)
        doubleprecision :: x,y,xpos,ypos,sigma,amplitude
            system_gauss = amplitude*exp(-sigma*(( x - xpos )**2+( y - ypos )**2)*DX*DX)
    end function system_gauss

    ! rozklad femiego
    double precision function system_fermi(E,Ef) result( rval )
    double precision :: E,Ef
        rval = 1.0/(exp( (E-Ef)/kbT ) + 1)
    end function  system_fermi
    double precision function system_fermiT(E,Ef,T) result( rval )
    double precision :: E,Ef,T
        rval = 1.0/(exp( (E-Ef)/T ) + 1)
    end function  system_fermiT
    ! pochodna rozkladu rozklad femiego
    double precision function system_dfermidE(E,Ef,dE) result( rval )
    double precision :: E,Ef,dE
        rval = (system_fermi(E+dE,Ef)-system_fermi(E-dE,Ef))/dE
    end function  system_dfermidE

    ! -------------------------------------------------------------------
    ! Podajemy    w jednostkach SI U w meV dx,dy w nm.
    ! -------------------------------------------------------------------
    subroutine system_dodaj_lorentza(U,ldx,ldy,x0,y0)
        double precision, intent(in) :: U,ldx,ldy,x0,y0
        integer :: i,j,rozbieg,ni,nj,mi,mj,nrz
        double precision :: x,y,ddx,ddy,dU,dx0,dy0



        dU = U/1000.0/Rd
        ddx= ldx*L2LR
        ddy= ldy*L2LR
        dx0= x0*L2LR
        dy0= y0*L2LR

        if(TRANS_DEBUG == 1 ) then
            print*," Dodawanie Lorentza:",U,ldx,ldy,x0,y0
            print*," Dodawanie Lorentza[Donor]:",dU,ddx,ddy,dx0,dy0
        endif

        do i = 1 , NX
        do j = 1 , NY
            x = i*DX
            y = j*DX
            UTOTAL(i,j) = UTOTAL(i,j) + dU/( 1 + ((x-dx0)/ddx)**2 + ((y-dy0)/ddy)**2 )
        enddo
        enddo

        ! zerowanie potencjalu przy wejsciach
        rozbieg = 20
        do nrz = 1 , no_zrodel
            select case (zrodla(nrz)%bKierunek)
            case (ZRODLO_KIERUNEK_PRAWO)
                ni = zrodla(nrz)%polozenia(1,1) + 1
                nj = zrodla(nrz)%polozenia(1,2)
                mi = ni + rozbieg - 1
                mj = zrodla(nrz)%polozenia(zrodla(nrz)%N,2)
                UTOTAL(ni:mi,nj:mj) = 0
            case (ZRODLO_KIERUNEK_LEWO)
                ni = zrodla(nrz)%polozenia(1,1) - 1
                nj = zrodla(nrz)%polozenia(1,2)
                mi = ni - rozbieg + 1
                mj = zrodla(nrz)%polozenia(zrodla(nrz)%N,2)
                UTOTAL(ni:mi,nj:mj) = 0
            case (ZRODLO_KIERUNEK_GORA)
                ni = zrodla(nrz)%polozenia(1,1)
                nj = zrodla(nrz)%polozenia(1,2) + 1
                mi = zrodla(nrz)%polozenia(zrodla(nrz)%N,1)
                mj = nj + rozbieg - 1
                UTOTAL(ni:mi,nj:mj) = 0
            case (ZRODLO_KIERUNEK_DOL)
                ni = zrodla(nrz)%polozenia(1,1)
                nj = zrodla(nrz)%polozenia(1,2) - 1
                mi = zrodla(nrz)%polozenia(zrodla(nrz)%N,1)
                mj = nj - rozbieg + 1
                UTOTAL(ni:mi,nj:mj) = 0
            endselect
        enddo



    end subroutine system_dodaj_lorentza

    ! -------------------------------------------------------------------
    ! Dodaje prostokatny obszar potencjalu z wygladzonymi za pomoca funkcji Diraca krawedziami.
    ! Podajemy
    ! xpos  - polozenie x [nm] paska potencjalu
    ! ypos1 - poczatek polozenia od dolu ([nm])
    ! ypos2 - koniec polozenia paska - gora [nm]
    ! amp   - amplituda w [meV]
    ! width - szerokosc paska w [nm]
    ! temp  - stopien rozmycia (posiada interpretacje temperatury z rozkladu Diraca)
    !         wyzsza temperatura wieksze rozmycie
    ! -------------------------------------------------------------------
    subroutine system_dodaj_pionowy_slupek_potencjalu(xpos,ypos1,ypos2,amp,width,temp)
        doubleprecision , intent(in) :: xpos,ypos1,ypos2,amp,width,temp
        double precision :: xp,yp,pvalue
        integer :: i,j

        do i = 1 , nx
        do j = 1 , ny
            xp = i * dx * LR2L
            yp = j * dx * LR2L

            pvalue = system_fermiT(xp,xpos-width/2,temp) + (1 - system_fermiT(xp,xpos+width/2,temp))
            pvalue = 1 - pvalue
            pvalue = pvalue * ( system_fermiT(-yp,-ypos1,temp) * system_fermiT(yp,ypos2,temp)   )

            UTOTAL(i,j) = UTOTAL(i,j) + pvalue * amp / Rd / 1000.0;
        enddo
        enddo
    end subroutine system_dodaj_pionowy_slupek_potencjalu


    ! -------------------------------------------------------------------
    ! Oblicza gestosc pradu prawdopodobienstwa na podstawie obecnej f.f.
    ! -------------------------------------------------------------------
    subroutine system_oblicz_J()
          integer          :: i,j
          complex*16       :: ycy,cx,jx,jy


          ! -----------------------------------------------------------------
          do i = 2 , nx-1
          do j = 2 , ny-1
              if(GFLAGS(i,j) == 0) then

                cx  = EXP(+II*DX*DX*(j)*BZ)

                ycy = conjg(cx)*VPHI(GINDEX(i+1,j))-VPHI(GINDEX(i-1,j))*cx
                jx  = (-1.0/4.0/DX*II)*(conjg(VPHI(GINDEX(i,j)))*ycy - VPHI(GINDEX(i,j))*conjg(ycy))
                ycy = VPHI(GINDEX(i,j+1))-VPHI(GINDEX(i,j-1))
                jy  = (-1.0/4.0/DX*II)*(conjg(VPHI(GINDEX(i,j)))*ycy - VPHI(GINDEX(i,j))*conjg(ycy))
                CURRENT(i,j) = CURRENT(i,j) + sqrt(DBLE(jx*jx + jy*jy))
              endif
          enddo
          enddo

    endsubroutine system_oblicz_J
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
                write(86554,"(3f20.6)"),i*DX*Lr2L,j*DX*Lr2L,fval
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
    case(ZAPISZ_PHI)
        do i = 1 , NX
        do j = 1 , NY
                fval = abs(PHI(i,j))
                write(86554,"(3e20.6)"),i*DX*Lr2L,j*DX*Lr2L,fval
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case(ZAPISZ_J)
        do i = 1 , NX
        do j = 1 , NY
                fval = CURRENT(i,j)
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
          integer,intent(in),dimension(:,:)  :: rows_cols
          integer,intent(inout),dimension(:) :: out_rows
          integer :: iterator, irow
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
      end subroutine convert_to_HB


    subroutine solve_system(no_rows,no_vals,colptr,rowind,values,b)
        integer,intent(in)                 :: no_rows
        integer,intent(in)                 :: no_vals
        integer,intent(in),dimension(:)    :: colptr,rowind
        complex*16,intent(in),dimension(:) :: values
        complex*16,intent(inout),dimension(:) :: b

        integer n, nnz, nrhs, ldb, info, iopt
        integer*8 factors
!
!      call zhbcode1(n, n, nnz, values, rowind, colptr)
!

        n    = no_rows
        nnz  = no_vals
        ldb  = n
        nrhs = 1


! First, factorize the matrix. The factors are stored in *factors* handle.
      iopt = 1
      call c_fortran_zgssv( iopt, n, nnz, nrhs, values, colptr , rowind , b, ldb,factors, info )
!
      if (info .eq. 0) then
         write (*,*) 'Factorization succeeded'
      else
         write(*,*) 'INFO from factorization = ', info
      endif
!
! Second, solve the system using the existing factors.
      iopt = 2
      call c_fortran_zgssv( iopt, n, nnz, nrhs, values, colptr,rowind ,  b, ldb,factors, info )
!
      if (info .eq. 0) then
         write (*,*) 'Solve succeeded'
!         write (*,*) (b(i), i=1, n)
      else
         write(*,*) 'INFO from triangular solve = ', info
      endif

! Last, free the storage allocated inside SuperLU
      iopt = 3
      call c_fortran_zgssv( iopt, n, nnz, nrhs, values, colptr,rowind, b, ldb,factors, info )

      endsubroutine solve_system

end module
