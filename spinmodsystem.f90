module modspinsystem
    use modspinzrodlo
    use modjed
    use modutils
!DEC$ IF DEFINED  (USE_UMF_PACK)
    use mUMFPACK
!DEC$ ENDIF

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
    doubleprecision  :: TRANS_R ! prawdopodobienstwo odbicia
    doubleprecision  :: TRANS_T ! prawdopodobienstwo przejscia
    doubleprecision :: Bz , Bx , By , Ef , Dx

    type(cspinzrodlo),dimension(:),allocatable    :: zrodla    ! TABLICA Z OBIEKTAMI ZRODLA
!    type(cabs_zrodlo),dimension(:),allocatable     :: abs_zrodla ! TABLICA Z OBIEKTAMI WIRTUALNYCH W.B.
    integer,dimension(:,:), allocatable           :: GFLAGS ! ZAWIERA FLAGI UKLADU (DIRICHLET,WEJSCIA,POZA)
    integer,dimension(:,:), allocatable           :: ZFLAGS ! ZAWIERA FLAGI WEJSC WARTOSC FLAGI OZNACZA NUMER WEJSCIA - 1,2, ...
    integer,dimension(:,:,:), allocatable         :: GINDEX ! INDEKSUJE GEOMETRIE UKLADU (-1) OZNACZA POZA UKLADEM, jako trzeci argument dochodzi spin -1,+1
    double precision,dimension(:,:), allocatable  :: UTOTAL ! MACIERZ POTENCJALU EFEKTYWNEGO
    double precision,dimension(:,:), allocatable  :: DUTOTAL ! MACIERZ POTENCJALU EFEKTYWNEGO w jednostkach donorowych
    double precision,dimension(:,:), allocatable  :: DEX,DEY ! MACIERZ POLA ELEKTRYCZNEGO w jednostkach donorowych
    complex*16,dimension(:),allocatable           :: CMATA  ! GLOWNA MACIERZ PROGRAMU W FORMACIE (ROW,COL,VALS), TUTAJ TYLKO VALS
    integer,dimension(:,:),allocatable            :: IDXA   ! INDEKSY MACIERZY (ROW,COL)
    complex*16,dimension(:),allocatable           :: VPHI   ! SZUKANA FUNKCJA FALOWA, PRZELICZANA POTEM NA LDOS
    complex*16,dimension(:,:,:),allocatable       :: PHI   ! FUNKCJA FALOWA ZAPISANA NA DWUWYMIAROWEJ MACIERZY
    double precision,dimension(:,:,:), allocatable:: CURRENT! FUNKCJA GESTOSCI PRADU PRAWDOPODOBIENSTWA
    double precision,dimension(:,:,:), allocatable:: CURRENTXY ! WETOR GESTOSCI PRADU PRAWDO
    double precision,dimension(:,:)  , allocatable:: DIVJ ! DYWERGENCJA WEKTORA J
    double precision,dimension(:,:,:), allocatable:: POLARYZACJE ! ZAWIERA LOKALNE POLARYZACJE SPINOWE (x,y,z)
!    complex*16,dimension(:,:,:),allocatable       :: WAVEFUNC  ! FUNKCJA FALOWA ZAPISANA NA DWUWYMIAROWEJ MACIERZY, TRZECIA KOLUMNA TO N-TY MOD

    ! Problem wlasny
!    integer                                       :: Widmo_NoStates
!    double precision,dimension(:), allocatable          :: Widmo_Evals
!    complex*16,dimension(:,:) ,allocatable        :: Widmo_Vecs

    ! rodzaje zapisu do pliku
    ENUM,BIND(C)
        ENUMERATOR :: ZAPISZ_FLAGI       = 0
        ENUMERATOR :: ZAPISZ_POTENCJAL   = 1
        ENUMERATOR :: ZAPISZ_KONTUR      = 2
        ENUMERATOR :: ZAPISZ_INDEKSY     = 3
        ENUMERATOR :: ZAPISZ_PHI         = 4
        ENUMERATOR :: ZAPISZ_J_ALL       = 5
        ENUMERATOR :: ZAPISZ_J_TOTAL     = 10
        ENUMERATOR :: ZAPISZ_DIVJ        = 11
        ENUMERATOR :: ZAPISZ_POLARYZACJE = 12
        ENUMERATOR :: ZAPISZ_WAVEFUNC    = 6
        ENUMERATOR :: ZAPISZ_STANY_WLASNE= 7
        ENUMERATOR :: ZAPISZ_WIDMO_VRTCAL= 8
        ENUMERATOR :: ZAPISZ_WIDMO_HRZNTL= 9
    END ENUM

    public :: ZAPISZ_FLAGI , ZAPISZ_POTENCJAL , ZAPISZ_KONTUR , ZAPISZ_INDEKSY , ZAPISZ_PHI
    public :: ZAPISZ_J_ALL, ZAPISZ_DIVJ, ZAPISZ_J_TOTAL, ZAPISZ_WAVEFUNC , ZAPISZ_POLARYZACJE
    public :: zrodla, UTOTAL, GFLAGS, GINDEX !,PHI
    public :: spinsystem_inicjalizacja , spinsystem_zwalnienie_pamieci , spinsystem_inicjalizacja_ukladu
!    public :: system_dodaj_abs_zrodlo !(pY1,pYN,pX1,pEf,pKierunek)
    public :: spinsystem_zapisz_do_pliku , spinsystem_rozwiaz_problem
    public :: spinsystem_gauss, spinsystem_dodaj_lorentza, spinsystem_dodaj_pionowy_slupek_potencjalu
    public :: spinsystem_fermi, spinsystem_fermiT, spinsystem_dfermidE
!    public :: system_widmo,system_zapisz_widmo_do_pliku,Widmo_NoStates,Widmo_Evals,Widmo_Vecs
!    public :: ZAPISZ_STANY_WLASNE , ZAPISZ_WIDMO_VRTCAL , ZAPISZ_WIDMO_HRZNTL
!    public :: TRANS_T,TRANS_R
    contains


    ! =========================================================
    !               INICJALIZACJA SYSTEMU
    ! =========================================================
    ! PROCEDURA INICJALIZACYJNA WYWOLYWANA JAKO PIERWSZA
    ! PNX,PNY       - ROZMIAR PROBLEMU, WYMIAR MACIERZY GINDEX, GFLAGS
    ! pLiczbaZrodel - LICZBA WSZYTKICH ZRODEL W UKLADZIE
    ! ----------------------------------------------------------
    subroutine spinsystem_inicjalizacja(pnx,pny,pLiczbaZrodel)
        integer,intent(in)             :: pnx,pny
        integer,intent(in)             :: pLiczbaZrodel


        if(TRANS_DEBUG==.true.) then
            print*,"System: Inicjalizacja:";
            print*,"    nx          =",pnx
            print*,"    ny          =",pny
            print*,"    liczbaZrodel=",pLiczbaZrodel
        endif
        nx = pnx
        ny = pny


        ! alokacja pamieci
        allocate( GFLAGS(nx,ny))
        allocate( ZFLAGS(nx,ny))
        allocate( GINDEX(nx,ny,-1:1))
        allocate( UTOTAL(nx,ny))
        allocate(DUTOTAL(nx,ny))
        allocate(    PHI(nx,ny,-1:1))
        allocate(    DEX(nx,ny))
        allocate(    DEY(nx,ny))
        allocate(CURRENT(nx,ny,5))
        allocate(DIVJ(nx,ny))
        allocate(CURRENTXY(nx,ny,2))
        allocate(POLARYZACJE(nx,ny,3))


        GFLAGS = B_NORMAL
        ZFLAGS = 0
        GINDEX = 0
        UTOTAL = 0
        DUTOTAL= 0
        PHI    = 0
        ! Tworzymy uklad...
        no_zrodel = pLiczbaZrodel
        allocate(zrodla(no_zrodel))
        ! okreslanie rzeczywistego rozmiaru problemu, punkty o indeksach -1 leza poza ukladem, nie beda
        ! potrzebne wiec zatem do obliczen.
        TRANS_MAXN = 0;
    end subroutine spinsystem_inicjalizacja
!
!    ! ========================================================================================
!    !               FUNKCJA ALOKUJE NOWE ZRODLO Z
!    !   TRANSPARENTNYMI WARUNKAMI BRZEGOWYMI. ALOKACJA PRZEPRO-
!    ! WADZANA JEST DYNAMICZNIE, DLATEGO NIE TRZEBA PODAWAC INDEKSOW.
!    ! ========================================================================================
!    ! Ustawiamy wczesniej zaalokowane zrodlo przez polecenie systemowe.
!    ! Parametry:
!    ! pY1,pYN - polozenia y na siatce numerycznej liczone od 1 do NY, w przypadku zrodla poziomego
!    !           oznaczaja polozenia X1 oraz XN
!    ! pX1 - polozenie X zrodla, dla zrodla poziomego polozenie Y1
!    ! pEf -Ef [meV]
!    ! pKierunek - enum ZRODLO_KIERUNEK_PRAWO/LEWO/GORA/DOL - ustala w ktora skierowane jest zrodlo
!    ! --------------------------------------------------------------------
!    subroutine system_dodaj_abs_zrodlo(pY1,pYN,pX1,pEf,pKierunek)
!        integer,intent(in)         ::  pY1,pYN,pX1
!        doubleprecision,intent(in) ::  pEf
!        integer,intent(in)         ::  pKierunek ! enum
!        type(cabs_zrodlo),dimension(:),allocatable     :: tmp_abs_zrodla
!        integer :: i
!        no_abs_zrodel = size(abs_zrodla)
!
!
!        allocate(tmp_abs_zrodla(no_abs_zrodel+1))
!        do i = 1 , no_abs_zrodel
!            call tmp_abs_zrodla(i)%abs_zrodlo_skopiuj(abs_zrodla(i))
!            call abs_zrodla(i)%abs_zrodlo_zwolnij_pamiec()
!        enddo
!
!        if(allocated(abs_zrodla)) deallocate(abs_zrodla)
!        allocate(abs_zrodla(no_abs_zrodel+1))
!        do i = 1 , no_abs_zrodel
!            call abs_zrodla(i)%abs_zrodlo_skopiuj(tmp_abs_zrodla(i))
!        enddo
!
!        call abs_zrodla(no_abs_zrodel+1)%abs_zrodlo_ustaw(pY1,pYN,pX1,pEf,pKierunek)
!        no_abs_zrodel = no_abs_zrodel + 1
!    end subroutine system_dodaj_abs_zrodlo

    ! =========================================================
    !
    !           CZYSZCZENIE ZAALOKOWANEJ PAMIECI
    !       funkcja wywolywana jest na samym koncu po wykonaniu
    !       obliczen
    ! =========================================================

    subroutine spinsystem_zwalnienie_pamieci()
        integer :: i
        if(TRANS_DEBUG==.true.) print*,"System: Zwalnianie pamieci"
        print*,"    Czyszczenie tablic..."
        if(allocated(GFLAGS))   deallocate(GFLAGS)
        if(allocated(ZFLAGS))   deallocate(ZFLAGS)
        if(allocated(GINDEX))   deallocate(GINDEX)
        if(allocated(UTOTAL))   deallocate(UTOTAL)
        if(allocated(DUTOTAL))  deallocate(DUTOTAL)
        if(allocated(CURRENT))  deallocate(CURRENT)
        if(allocated(CURRENTXY))deallocate(CURRENTXY)
        if(allocated(POLARYZACJE))deallocate(POLARYZACJE)
        if(allocated(DIVJ))     deallocate(DIVJ)
        if(allocated(PHI))      deallocate(PHI)
        if(allocated(DEX))      deallocate(DEX)
        if(allocated(DEY))      deallocate(DEY)
!        if( allocated(WAVEFUNC) ) deallocate(WAVEFUNC)
!        if(allocated(Widmo_Evals)) deallocate(Widmo_Evals)
!        if(allocated(Widmo_Vecs)) deallocate(Widmo_Vecs)
        print*,"    Czyszczenie zrodel..."
        do i = 1 , no_zrodel
            call zrodla(i)%spinzrodlo_zwolnij_pamiec()
        enddo
        print*,"    Czyszczenie abs_zrodel..."
        print*,"Brak transparentnych"
!        do i = 1 , no_abs_zrodel
!            call abs_zrodla(i)%abs_zrodlo_zwolnij_pamiec()
!        enddo
        if(allocated(zrodla))    deallocate(zrodla)
!        if(allocated(abs_zrodla)) deallocate(abs_zrodla)
    end subroutine spinsystem_zwalnienie_pamieci

    ! ------------------------------------------------------------
    ! Funkcja tworzy tablice flag oraz gindex na podstawie
    ! zmodyfikowanej w mainie tablicy gflags
    ! Podajemy:
    ! in_len - dlugosc na jaka od wejsc zostanie utworzona flaga B_NORMAL
    ! smooth_rad - promien wygladzania ostrych katow w ukladzie
    ! smooth_iter - liczba iteracji wygladzania
    ! ------------------------------------------------------------caliente
    subroutine spinsystem_inicjalizacja_ukladu(in_len,smooth_rad,smooth_iter)
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
                print*,"tutaj moze byc zle liczone"
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
        print*,"Brak transparentnych"
!        ! przypisywanie  transparentnych w.b.
!        do nrz = 1 , no_abs_zrodel
!            do i = 2 , abs_zrodla(nrz)%N - 1
!                ni = abs_zrodla(nrz)%polozenia(i,1)
!                nj = abs_zrodla(nrz)%polozenia(i,2)
!                GFLAGS(ni,nj) = B_TRANSPARENT
!                ZFLAGS(ni,nj) = nrz ! przypisujemy fladze numer wejscia
!            enddo
!        enddo

        ! -------------------------------------------------------------
        ! wygladzanie krawedzi
        ! -------------------------------------------------------------
		do soft_iter = 1 , smooth_iter
        do i =  1+smooth_rad , nx-smooth_rad-1
        do j =  1+smooth_rad , ny-smooth_rad-1
         norm_iter  = 0
		 empty_iter = 0
			do iii =  -smooth_rad , smooth_rad , 1
			do jjj =  -smooth_rad , smooth_rad , 1
				if(GFLAGS(i+iii,j+jjj) == B_EMPTY )  empty_iter = empty_iter + 1
				if(GFLAGS(i+iii,j+jjj) == B_NORMAL ) norm_iter = norm_iter   + 1
			enddo
			enddo
			if(empty_iter < norm_iter .and. GFLAGS(i,j) == B_EMPTY  ) GFLAGS(i,j) = B_NORMAL
            if(empty_iter > norm_iter .and. GFLAGS(i,j) == B_NORMAL ) GFLAGS(i,j) = B_EMPTY
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
        do iii = -1 , 1 , 2

        do i = 1 , nx
        do j = 1 , ny
           if( GFLAGS(i,j)  == B_WEJSCIE      .or. &
              & GFLAGS(i,j) == B_NORMAL       .or. &
              & GFLAGS(i,j) == B_TRANSPARENT  .or. &
              & GFLAGS(i,j) == B_DIRICHLET  ) then
                GINDEX(i,j,-iii) = iter
                iter = iter + 1
           endif
        enddo
        enddo

        enddo ! end of iii - numerowanie po spinach

    end subroutine spinsystem_inicjalizacja_ukladu


    ! --------------------------------------------------------------------
    ! Rozwiazywanie problemu dla zrodla o podanym numerze - nrz
    ! --------------------------------------------------------------------
    subroutine spinsystem_rozwiaz_problem(nrz,TR_MAT)
        integer,intent(in)  :: nrz
        double precision,dimension(:,:), allocatable  :: TR_MAT
        integer,allocatable :: HBROWS(:)
        doubleprecision :: YA, YB
        complex*16 ::  Zup , Zdwn
        integer :: i,j,ni,nj,mod_in,dir

        Ef  = atomic_Ef/1000.0/Rd
        BZ  = BtoDonorB(atomic_Bz)
        Bx  = BtoDonorB(atomic_Bx)
        By  = BtoDonorB(atomic_By)
        DX  = atomic_DX*L2LR

        print*,"! ----------------------------------------------- !"
        print*,"! Rozpoczecie obliczen dla zrodla=", nrz
        print*,"!   Ef=", Ef*Rd*1000    , "[meV]"
        print*,"!   Bz=", DonorBtoB(Bz) , "[T]"
        print*,"! ----------------------------------------------- !"


        call reset_clock()
        call oblicz_rozmiar_macierzy()
!        print*,"Oblicz rozmiar macierzy:", get_clock()
        call reset_clock()


        allocate(CMATA(MATASIZE))
        allocate(IDXA (MATASIZE,2))
        allocate(VPHI(TRANS_MAXN))
        allocate(HBROWS(TRANS_MAXN+1))
!        if( allocated(WAVEFUNC) ) deallocate(WAVEFUNC)
!        allocate(WAVEFUNC(nx,ny,zrodla(nrz)%liczba_modow))

        ! tablice z wspolczynnikami T i R
        if(allocated(TR_MAT))deallocate(TR_MAT)
        allocate(TR_MAT(no_zrodel,zrodla(nrz)%liczba_modow))
        TR_MAT = 0
        CURRENT = 0;
        DIVJ    = 0
        CURRENTXY = 0
        POLARYZACJE = 0
        PHI     = 0
        TRANS_R = 0
        TRANS_T = 0

        DUTOTAL = UTOTAL/1000.0/Rd

        ! pole elekrtyczne jest trzymane w jednostkach donorowych
        do i = 1 , NX
        do j = 2 , NY-1
            DEY(i,j) = (DUTOTAL(i,j+1) - DUTOTAL(i,j-1))/2/DX
        enddo
        enddo
        DEY(:,1)  = DEY(:,2)
        DEY(:,NY) = DEY(:,NY-1)

        do i = 2 , NX-1
        do j = 1 , NY
            DEX(i,j) = (DUTOTAL(i+1,j) - DUTOTAL(i-1,j))/2/DX
        enddo
        enddo
        DEX(1,:)  = DEX(2   ,:)
        DEX(NX,:) = DEX(NX-1,:)


        ! --------------------------------------------------------------------
        !
        ! --------------------------------------------------------------------
        call wypelnij_macierz()
        call convert_to_HB(MATASIZE,IDXA,HBROWS)
!        print*,"Tworzenie macierzy:", get_clock()
        call reset_clock()
        call solve_system(TRANS_MAXN,MATASIZE,IDXA(:,2),HBROWS,CMATA(:),VPHI,1)
!        print*,"Pierwszy solve:", get_clock()
        call reset_clock()

        mod_in = 1
        do mod_in = 1 , zrodla(nrz)%liczba_modow
        VPHI = 0
        zrodla(nrz)%ck(:)      = 0
        zrodla(nrz)%ck(mod_in) = 1
        call zrodla(nrz)%spinzrodlo_oblicz_Fj()
        do i = 2 , zrodla(nrz)%N - 1
            ni = zrodla(nrz)%polozenia(i,1)
            nj = zrodla(nrz)%polozenia(i,2)
            dir = zrodla(nrz)%dir

            VPHI(GINDEX(ni,nj,+1)) = -( zrodla(nrz)%Fj(i,+1)*Tu(-dir,+1,ni,nj) + zrodla(nrz)%Fj(i,-1)*Su(-dir,-1,ni,nj))
            VPHI(GINDEX(ni,nj,-1)) = -( zrodla(nrz)%Fj(i,+1)*Su(-dir,+1,ni,nj) + zrodla(nrz)%Fj(i,-1)*Tu(-dir,-1,ni,nj))
        enddo



        call solve_system(TRANS_MAXN,MATASIZE,IDXA(:,2),HBROWS,CMATA(:),VPHI,2)
        call system_oblicz_J()
        call oblicz_TR(nrz,mod_in)

        do i = 1 ,no_zrodel
            TR_MAT(i,mod_in) = sum( abs(zrodla(i)%TR(1:zrodla(i)%liczba_modow,-zrodla(i)%dir)) )
        enddo

        do i = 1 , Nx
        do j = 1 , Ny
            if(GINDEX(i,j,+1) > 0) PHI(i,j,+1) = PHI(i,j,+1) + abs( VPHI(GINDEX(i,j,+1)) )**2
            if(GINDEX(i,j,-1) > 0) PHI(i,j,-1) = PHI(i,j,-1) + abs( VPHI(GINDEX(i,j,-1)) )**2

            if(GINDEX(i,j,+1) > 0) then
            Zup  = VPHI(GINDEX(i,j,+1))
            Zdwn = VPHI(GINDEX(i,j,-1))
            ! Kierunek Z
            YA = (abs(Zup)**2)
            YB = (abs(Zdwn)**2)
            POLARYZACJE(i,j,1) = POLARYZACJE(i,j,1) + (YA-YB)!/(YA+YB) ! kierunek z

            ! Kierunek X
            YA = (abs(Zup+Zdwn)**2) ! up
            YB = (abs(Zup-Zdwn)**2) ! down
            POLARYZACJE(i,j,2) = POLARYZACJE(i,j,2) + (YA-YB)!/(YA+YB) ! kierunek x

            ! Kierunek Y
            YA = (abs(Zup+II*Zdwn)**2)
            YB = (abs(Zup-II*Zdwn)**2)
            POLARYZACJE(i,j,3) = POLARYZACJE(i,j,3) + (YA-YB)!/(YA+YB) ! kierunek y
            endif
!            if(GINDEX(i,j) > 0) WAVEF  UNC(i,j,mod_in) = VPHI(GINDEX(i,j))
        enddo
        enddo


        enddo ! end of petla po modach


        print*,"T = ", TRANS_T
        print*,"R = ", TRANS_R
        print*,"W = ", TRANS_R + TRANS_T
        call solve_system(TRANS_MAXN,MATASIZE,IDXA(:,2),HBROWS,CMATA(:),VPHI,3)
	    zrodla(nrz)%ck(:)      = 0

        deallocate(CMATA)
        deallocate(IDXA)
        deallocate(VPHI)
        deallocate(HBROWS)

    end subroutine spinsystem_rozwiaz_problem
!
    ! ---------------------------------------------------------------------
    !   Wyliczanie rozmiaru macierzy: zakladamy, ze na :
    ! 1. warunek Dirichleta potrzeba jedna wartosc
    ! 2. normalny punkt wewnatrz ukladu potrzeba 9 wartosci (4 siasiadow + 1 wezel centralny + 4 sasiadow od spinu)
    ! 3. dla wejsc potrzeba N-2 + 1 wezlow: gdzie N-2 oznacza ze bierzemy
    !    wezly wejsciowe od 2 .. N -1, a +2 wynika z faktu, ze mamy zawsze
    !    polaczenie z wezlem (i+1,j), lub (i-1,j) razy spin
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
!                    print*,"dirichlet:",1
            case(B_WEJSCIE)
                    MATASIZE   = MATASIZE + 2*(zrodla(ZFLAGS(i,j))%N  - 2) + 2 ! +2 na relacje prawo/lewo/gora/dol -2 na obciecie indeksow (1 i N)
                    TRANS_MAXN = TRANS_MAXN + 1
!                    print*,"wejscie:",2*(zrodla(ZFLAGS(i,j))%N  - 2) + 2,zrodla(ZFLAGS(i,j))%N
            case(B_TRANSPARENT)
                    MATASIZE   = MATASIZE + 2 ! na transparentne w.b. potrzebujemy dwie komurki
                    TRANS_MAXN = TRANS_MAXN + 1
!                    print*,"transp:",2
            case(B_NORMAL)
                    MATASIZE = MATASIZE + 9 + 1 !+1 od oddzialywania z polem (Zemann)
                    TRANS_MAXN = TRANS_MAXN + 1
!                    print*,"normal:",9
            end select
        enddo
        enddo
        MATASIZE   = MATASIZE   * 2 ! *2 bo spin
        TRANS_MAXN = TRANS_MAXN * 2
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
    complex*16 function T0(d,s,u,v) result(rval)
        integer :: d,s,u,v
        rval = 4 / 2.0 / DX / DX + DUTOTAL(u,v) + 0.5*s*G_LAN*M_EFF*BZ - Ef
    end function T0

    complex*16 function Tu(d,s,u,v) result(rval)
        integer :: d,s,u,v
        rval = - 1.0 / 2.0 / DX / DX  * EXP(-d*II*DX*DX*(v)*BZ) * EXP(+d*s*II*DX*DEY(u,v)*so_loc)
    end function Tu

    complex*16 function Tv(d,s,u,v) result(rval)
        integer :: d,s,u,v
        rval = - 1.0 / 2.0 / DX / DX * EXP(-d*s*II*DX*DEX(u,v)*so_loc)
    end function Tv

    complex*16 function Su(d,s,u,v) result(rval)
        integer :: d,s,u,v
        rval = - s*d*so_rashba / 2.0 / DX  * EXP(-d*II*DX*DX*(v)*BZ)
    end function Su

    complex*16 function Sv(d,s,u,v) result(rval)
        integer :: d,s,u,v
        rval = - d*II*so_rashba / 2.0 / DX
    end function Sv

    complex*16 function Sb(d,s,u,v) result(rval)
        integer :: d,s,u,v
        rval = 0.5*G_LAN*M_EFF*(Bx + s * II * By )
    end function Sb

    subroutine wypelnij_macierz()

        ! zmienne pomocniczne
        integer          :: itmp,ni,nj,pnj,nn,ln,nzrd,pni,s,u,v,dir
        complex*16       :: post
        doubleprecision  :: kvec
        itmp  = 1
        do s = +1 , -1 , -2 ! petla po spinach
        do u = 1, nx
        do v = 1, ny
            if(GINDEX(u,v,s) > 0) then
                if( GFLAGS(u,v) == B_DIRICHLET ) then

                    cmatA(itmp)  = CMPLX(1.0)
                    idxA(itmp,1) = GINDEX(u, v , s)
                    idxA(itmp,2) = GINDEX(u, v , s)
                    itmp = itmp + 1

                elseif( GFLAGS(u,v) == B_NORMAL) then

                    cmatA(itmp)   = T0(0,s,u,v)
                    idxA (itmp,1) = GINDEX(u,v , s)
                    idxA (itmp,2) = GINDEX(u,v , s)
                    itmp = itmp + 1

                    cmatA(itmp)   = Tu(-1,s,u,v)
                    idxA (itmp,1) = GINDEX(u  ,v,s)
                    idxA (itmp,2) = GINDEX(u-1,v,s)
                    itmp = itmp + 1

                    cmatA(itmp)   = Tu(+1,s,u,v)
                    idxA (itmp,1) = GINDEX(u  ,v,s)
                    idxA (itmp,2) = GINDEX(u+1,v,s)
                    itmp = itmp + 1

                    cmatA(itmp)  = Tv(+1,s,u,v)
                    idxA(itmp,1) = GINDEX(u,v  ,s)
                    idxA(itmp,2) = GINDEX(u,v+1,s)
                    itmp = itmp + 1

                    cmatA(itmp)  = Tv(-1,s,u,v)
                    idxA(itmp,1) = GINDEX(u,v  ,s)
                    idxA(itmp,2) = GINDEX(u,v-1,s)
                    itmp = itmp + 1

                    ! oddzialywanie typu rashba
                    cmatA(itmp)   = Su(-1,-s,u,v)
                    idxA (itmp,1) = GINDEX(u  ,v,+s)
                    idxA (itmp,2) = GINDEX(u-1,v,-s)
                    itmp = itmp + 1

                    cmatA(itmp)   = Su(+1,-s,u,v)
                    idxA (itmp,1) = GINDEX(u  ,v,+s)
                    idxA (itmp,2) = GINDEX(u+1,v,-s)
                    itmp = itmp + 1

                    cmatA(itmp)  = Sv(+1,-s,u,v)
                    idxA(itmp,1) = GINDEX(u,v  ,+s)
                    idxA(itmp,2) = GINDEX(u,v+1,-s)
                    itmp = itmp + 1

                    cmatA(itmp)  = Sv(-1,-s,u,v)
                    idxA(itmp,1) = GINDEX(u,v  ,+s)
                    idxA(itmp,2) = GINDEX(u,v-1,-s)
                    itmp = itmp + 1

                    ! odzialywanie z polem w kierunku x i y (zeemann)
                    cmatA(itmp)  = Sb(0,-s,0,0)
                    idxA(itmp,1) = GINDEX(u,v  ,+s)
                    idxA(itmp,2) = GINDEX(u,v  ,-s)
                    itmp = itmp + 1




                ! ----------------------------------------------------------------------
                ! Obsluga wejsc
                ! ----------------------------------------------------------------------
                else if( GFLAGS(u,v) == B_WEJSCIE) then
                    nzrd = ZFLAGS(u,v)
                    dir  = zrodla(nzrd)%dir


                    ni   = u
                    nj   = v ! globalne polozenia na siatce

                    select case (zrodla(nzrd)%bKierunek)
                    ! ------------------------------------------------------------------
                    ! Zrodla prawe:
                    !
                    !                      ------------------------>
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_PRAWO,ZRODLO_KIERUNEK_LEWO)

                        ln   = v - zrodla(nzrd)%polozenia(1,2) + 1 ! lokalny indeks

                        cmatA(itmp)   = (Tu(+dir,+s,u,v) + Tu(-dir,+s,u,v))
!                        if( ZRODLO_KIERUNEK_PRAWO == zrodla(nzrd)%bKierunek ) cmatA(itmp)  = 0.0
                        idxA (itmp,1) = GINDEX(u    ,v,+s)
                        idxA (itmp,2) = GINDEX(u+dir,v,+s)
                        itmp = itmp + 1

                        cmatA(itmp)   = (Su(+dir,-s,u,v) + Su(-dir,-s,u,v))
!                        if( ZRODLO_KIERUNEK_PRAWO == zrodla(nzrd)%bKierunek ) cmatA(itmp)  = 0.0
                        idxA (itmp,1) = GINDEX(u    ,v,+s)
                        idxA (itmp,2) = GINDEX(u+dir,v,-s)
                        itmp = itmp + 1

                    ! ------------------------------------------------------------------
                    ! Zrodla dolne:
                    !
                    !                      ------------------------>
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_GORA)

!                        ln   = i - zrodla(nzrd)%polozenia(1,1) + 1
!
!                        cmatA(itmp)   = CMPLX(-0.5/DX/DX*2)
!                        idxA (itmp,1) = GINDEX(i  ,j)
!                        idxA (itmp,2) = GINDEX(i,j+1)
!                        itmp = itmp + 1
!
!                        post = 0.5/DX/DX
                    ! ------------------------------------------------------------------
                    ! Zrodla gorne:
                    !
                    !                      <------------------------
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_DOL)
!
!                        ln   = i - zrodla(nzrd)%polozenia(1,1) + 1
!
!                        cmatA(itmp)   = CMPLX(-0.5/DX/DX*2)
!                        idxA (itmp,1) = GINDEX(i  ,j)
!                        idxA (itmp,2) = GINDEX(i,j-1)
!                        itmp = itmp + 1
!
!                        post = -0.5/DX/DX
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

                        cmatA(itmp)   = T0(0,s,u,v) +&
                                zrodla(nzrd)%Sigma(ln,nn,+s,+s)*Tu(-dir,+s,u,v) + &
                                zrodla(nzrd)%Sigma(ln,nn,-s,+s)*Su(-dir,-s,u,v)



!                        if( ZRODLO_KIERUNEK_PRAWO == zrodla(nzrd)%bKierunek ) cmatA(itmp)  = 1.0

                        idxA (itmp,1) = GINDEX(ni,nj ,+s)
                        idxA (itmp,2) = GINDEX(ni,pnj,+s)
                        itmp = itmp + 1


                        cmatA(itmp)  = Sb(0,-s,0,0) + &
                            zrodla(nzrd)%Sigma(ln,nn,+s,-s)*Tu(-dir,+s,u,v) + &
                            zrodla(nzrd)%Sigma(ln,nn,-s,-s)*Su(-dir,-s,u,v)


!                        if( ZRODLO_KIERUNEK_PRAWO == zrodla(nzrd)%bKierunek ) cmatA(itmp)  = 0.0
                        idxA(itmp,1) = GINDEX(ni,nj ,+s)
                        idxA(itmp,2) = GINDEX(ni,pnj,-s)
                        itmp = itmp + 1

                    else if( nn == ln+1  ) then

                        cmatA(itmp)   = Tv(+1,s,u,v) +&
                                zrodla(nzrd)%Sigma(ln,nn,+s,+s)*Tu(-dir,+s,u,v) + &
                                zrodla(nzrd)%Sigma(ln,nn,-s,+s)*Su(-dir,-s,u,v)


!                        if( ZRODLO_KIERUNEK_PRAWO == zrodla(nzrd)%bKierunek ) cmatA(itmp)  = 0.0

                        idxA(itmp,1) = GINDEX(ni,nj ,+s)
                        idxA(itmp,2) = GINDEX(ni,pnj,+s)
                        itmp = itmp + 1



                        cmatA(itmp)  = Sv(+1,-s,u,v) + &
                            zrodla(nzrd)%Sigma(ln,nn,+s,-s)*Tu(-dir,+s,u,v) + &
                            zrodla(nzrd)%Sigma(ln,nn,-s,-s)*Su(-dir,-s,u,v)


!                        if( ZRODLO_KIERUNEK_PRAWO == zrodla(nzrd)%bKierunek ) cmatA(itmp)  = 0.0

                        idxA(itmp,1) = GINDEX(ni,nj ,+s)
                        idxA(itmp,2) = GINDEX(ni,pnj,-s)
                        itmp = itmp + 1


                    else if( nn == ln-1  ) then

                        cmatA(itmp)   = Tv(-1,s,u,v) +&
                                zrodla(nzrd)%Sigma(ln,nn,+s,+s)*Tu(-dir,+s,u,v) + &
                                zrodla(nzrd)%Sigma(ln,nn,-s,+s)*Su(-dir,-s,u,v)


!                        if( ZRODLO_KIERUNEK_PRAWO == zrodla(nzrd)%bKierunek ) cmatA(itmp)  = 0.0

                        idxA(itmp,1) = GINDEX(ni,nj ,+s)
                        idxA(itmp,2) = GINDEX(ni,pnj,+s)
                        itmp = itmp + 1

                        cmatA(itmp)  = Sv(-1,-s,u,v) + &
                            zrodla(nzrd)%Sigma(ln,nn,+s,-s)*Tu(-dir,+s,u,v) + &
                            zrodla(nzrd)%Sigma(ln,nn,-s,-s)*Su(-dir,-s,u,v)


!                        if( ZRODLO_KIERUNEK_PRAWO == zrodla(nzrd)%bKierunek ) cmatA(itmp)  = 0.0

                        idxA(itmp,1) = GINDEX(ni,nj ,+s)
                        idxA(itmp,2) = GINDEX(ni,pnj,-s)
                        itmp = itmp + 1

                    else
                        cmatA(itmp)   = &
                                zrodla(nzrd)%Sigma(ln,nn,+s,+s)*Tu(-dir,+s,u,v) + &
                                zrodla(nzrd)%Sigma(ln,nn,-s,+s)*Su(-dir,-s,u,v)


!                        if( ZRODLO_KIERUNEK_PRAWO == zrodla(nzrd)%bKierunek ) cmatA(itmp)  = 0.0

                        idxA(itmp,1) = GINDEX(ni,nj ,+s)
                        idxA(itmp,2) = GINDEX(ni,pnj,+s)
                        itmp = itmp + 1

                        cmatA(itmp)  = &
                            zrodla(nzrd)%Sigma(ln,nn,+s,-s)*Tu(-dir,+s,u,v) + &
                            zrodla(nzrd)%Sigma(ln,nn,-s,-s)*Su(-dir,-s,u,v)


!                        if( ZRODLO_KIERUNEK_PRAWO == zrodla(nzrd)%bKierunek ) cmatA(itmp)  = 0.0

                        idxA(itmp,1) = GINDEX(ni,nj ,+s)
                        idxA(itmp,2) = GINDEX(ni,pnj,-s)
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
!
!                    do nn = 2 , zrodla(nzrd)%N - 1
!                    pni   = zrodla(nzrd)%polozenia(nn,1)
!
!                    if( ln == nn  ) then
!
!                        cmatA(itmp)   = CMPLX( 2.0/DX/DX + UTOTAL(i,j) - Ef ) &
!                           & + post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
!                        idxA (itmp,1) = GINDEX(ni,nj)
!                        idxA (itmp,2) = GINDEX(pni,nj)
!                        itmp = itmp + 1
!
!                    else if( nn == ln+1  ) then
!                        cmatA(itmp) = CMPLX(-0.5/DX/DX)*(EXP(-II*DX*DX*(nj)*BZ)) &
!                        &   + post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
!                        idxA(itmp,1) = GINDEX(ni,nj)
!                        idxA(itmp,2) = GINDEX(pni,nj)
!                        itmp = itmp + 1
!                    else if( nn == ln-1  ) then
!
!                        cmatA(itmp) = CMPLX(-0.5/DX/DX)*(EXP(+II*DX*DX*(nj)*BZ)) &
!                        &   + post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
!                        idxA(itmp,1) = GINDEX(ni,nj)
!                        idxA(itmp,2) = GINDEX(pni,nj)
!                        itmp = itmp + 1
!                    else
!                        cmatA(itmp)  = post*zrodla(nzrd)%zrodlo_alfa_v_i(dx,ln,nn)
!                        idxA(itmp,1) = GINDEX(ni,nj)
!                        idxA(itmp,2) = GINDEX(pni,nj)
!                        itmp = itmp + 1
!                    endif
!
!                    enddo
                    endselect

                ! ----------------------------------------------------------------------- !
                ! endif ! end of if WEJSCIE, transparentne warunki brzegowe
                ! ----------------------------------------------------------------------- !
                else if( GFLAGS(u,v) == B_TRANSPARENT) then

!                    nzrd = ZFLAGS(i,j)
!                    kvec = abs_zrodla(nzrd)%kvec
!                    ni   = i
!                    nj   = j ! globalne polozenia na siatce
!
!                    cmatA(itmp)   = 1
!                    idxA (itmp,1) = GINDEX(i, j)
!                    idxA (itmp,2) = GINDEX(i, j)
!                    itmp = itmp + 1
!
!                    cmatA(itmp)   =-EXP(II*KVEC*DX)
!                    idxA (itmp,1) = GINDEX(i, j)
!                    select case (abs_zrodla(nzrd)%bKierunek)
!                    ! ------------------------>
!                    case (ZRODLO_KIERUNEK_PRAWO)
!                        idxA (itmp,2) = GINDEX(i+1,j)
!                    !<------------------------
!                    case (ZRODLO_KIERUNEK_LEWO)
!                        idxA (itmp,2) = GINDEX(i-1,j)
!                    ! ------------------------>
!                    case (ZRODLO_KIERUNEK_GORA)
!                        idxA (itmp,2) = GINDEX(i, j+1)
!                    ! <------------------------
!                    case (ZRODLO_KIERUNEK_DOL)
!                        idxA (itmp,2) = GINDEX(i, j-1)
!                    endselect
!                    itmp = itmp + 1

                    endif ! rodzaj komorki - ostatni else B_TRANSPARENT
            endif ! end if GINDEX > 0
        enddo ! end of j
        enddo ! end of i
        enddo ! end of s - spiny

        if ( itmp -1 - MATASIZE /= 0 ) then
            print *, "Blad rozmiar maicierzy niezgodny z przewidywanym:" , itmp - 1
!            stop
             MATASIZE = itmp - 1
        else
            print *, "Macierz skonstruowana poprawnie..."
        endif

!    do u = 1 , MATASIZE
!        write(222,"(i10,A,2i5,A,2f10.4)"),u,"(row,col)=",IDXA(u,:)," value=",CMATA(u)
!    enddo


    end subroutine wypelnij_macierz

    ! -------------------------------------------------------------------
    ! Procedura oblicza prawdopodobienstwo przejscia T oraz odbicia R
    ! dla podanego zrodla wejsciowego nrz i danego modu wchodzacego mod_in
    ! -------------------------------------------------------------------
    subroutine oblicz_TR(nrz,mod_in)
        integer,intent(in) :: nrz,mod_in
        integer :: i,no_modes,dir
        doubleprecision :: Jin


        ! obliczamy apliduty prawdopodobienstwa
        do i = 1 , no_zrodel
            call zrodla(i)%spinzrodlo_oblicz_dk(VPHI,GINDEX,nx,ny)
        enddo


        ! na podstawie strumieni obliczamy transmisje oraz prawdopodobienstwo odbicia
        Jin = zrodla(nrz)%ChiCurr(mod_in,zrodla(nrz)%dir)
        call zrodla(nrz)%spinzrodlo_wypisz_JinJout()


        do i = 1 , no_zrodel
             zrodla(i)%TR = 0
             no_modes = zrodla(i)%liczba_modow
             dir = zrodla(i)%dir
             zrodla(i)%TR(1:no_modes,+dir) = abs(zrodla(i)%ck(1:no_modes))**2 * zrodla(i)%ChiCurr(1:no_modes,+dir) / Jin
             zrodla(i)%TR(1:no_modes,-dir) = abs(zrodla(i)%dk(1:no_modes))**2 * zrodla(i)%ChiCurr(1:no_modes,-dir) / Jin

        enddo

        TRANS_R = TRANS_R + sum( abs(zrodla(nrz)%TR(1:zrodla(nrz)%liczba_modow,-zrodla(nrz)%dir)) )

        do i = 1 ,no_zrodel
            dir = zrodla(i)%dir
            no_modes = zrodla(i)%liczba_modow
            if(i /= nrz) TRANS_T = TRANS_T + sum( abs(zrodla(i)%TR(1:no_modes,-dir)) )
        enddo


    end subroutine oblicz_TR
    ! -------------------------------------------------------------------
    ! Zwraca wartosc rozkladu gaussa w punkcie (x,y) dla gaussa o srodku
    ! w (xpos,ypos) oraz sigmie = sigma i amplidudzie = amplitude.
    ! x,y,xpos,ypos - podajemy w indeksach siatki (i,j)
    ! -------------------------------------------------------------------
    doubleprecision function spinsystem_gauss(x,y,xpos,ypos,sigma,amplitude)
        doubleprecision :: x,y,xpos,ypos,sigma,amplitude
            spinsystem_gauss = amplitude*exp(-sigma*(( x - xpos )**2+( y - ypos )**2)*DX*DX)
    end function spinsystem_gauss

    ! rozklad femiego
    double precision function spinsystem_fermi(E,Ef) result( rval )
    double precision :: E,Ef
        rval = 1.0/(exp( (E-Ef)/kbT ) + 1)
    end function  spinsystem_fermi
    double precision function spinsystem_fermiT(E,Ef,T) result( rval )
    double precision :: E,Ef,T
        rval = 1.0/(exp( (E-Ef)/T ) + 1)
    end function  spinsystem_fermiT
    ! pochodna rozkladu rozklad femiego
    double precision function spinsystem_dfermidE(E,Ef,dE) result( rval )
    double precision :: E,Ef,dE
        rval = (spinsystem_fermi(E+dE,Ef)-spinsystem_fermi(E-dE,Ef))/dE
    end function  spinsystem_dfermidE

    ! -------------------------------------------------------------------
    ! Podajemy    w jednostkach SI U w meV dx,dy w nm.
    ! -------------------------------------------------------------------
    subroutine spinsystem_dodaj_lorentza(U,ldx,ldy,x0,y0)
        double precision, intent(in) :: U,ldx,ldy,x0,y0
        integer :: i,j,rozbieg,ni,nj,mi,mj,nrz
        double precision :: x,y,ddx,ddy,dU,dx0,dy0



        dU = U!/1000.0/Rd
        ddx= ldx!*L2LR
        ddy= ldy!*L2LR
        dx0= x0!*L2LR
        dy0= y0!*L2LR

        if(TRANS_DEBUG == 1 ) then
            print*," Dodawanie Lorentza:",U,ldx,ldy,x0,y0
            print*," Dodawanie Lorentza[Donor]:",dU,ddx,ddy,dx0,dy0
        endif

        do i = 1 , NX
        do j = 1 , NY
            x = i*atomic_DX
            y = j*atomic_DX
            DUTOTAL(i,j) =  dU/( 1 + ((x-dx0)/ddx)**2 + ((y-dy0)/ddy)**2 )
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
                DUTOTAL(ni:mi,nj:mj) = 0
            case (ZRODLO_KIERUNEK_LEWO)
                ni = zrodla(nrz)%polozenia(1,1) - 1
                nj = zrodla(nrz)%polozenia(1,2)
                mi = ni - rozbieg + 1
                mj = zrodla(nrz)%polozenia(zrodla(nrz)%N,2)
                DUTOTAL(mi:ni,nj:mj) = 0
            case (ZRODLO_KIERUNEK_GORA)
                ni = zrodla(nrz)%polozenia(1,1)
                nj = zrodla(nrz)%polozenia(1,2) + 1
                mi = zrodla(nrz)%polozenia(zrodla(nrz)%N,1)
                mj = nj + rozbieg - 1
                DUTOTAL(ni:mi,nj:mj) = 0
            case (ZRODLO_KIERUNEK_DOL)
                ni = zrodla(nrz)%polozenia(1,1)
                nj = zrodla(nrz)%polozenia(1,2) - 1
                mi = zrodla(nrz)%polozenia(zrodla(nrz)%N,1)
                mj = nj - rozbieg + 1
                DUTOTAL(ni:mi,nj:mj) = 0
            endselect
        enddo

        UTOTAL = UTOTAL + DUTOTAL

    end subroutine spinsystem_dodaj_lorentza

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
    subroutine spinsystem_dodaj_pionowy_slupek_potencjalu(xpos,ypos1,ypos2,amp,width,temp)
        doubleprecision , intent(in) :: xpos,ypos1,ypos2,amp,width,temp
        double precision :: xp,yp,pvalue
        integer :: i,j

        do i = 1 , nx
        do j = 1 , ny
            xp = i * atomic_DX
            yp = j * atomic_DX

            pvalue = spinsystem_fermiT(xp,xpos-width/2,temp) + (1 - spinsystem_fermiT(xp,xpos+width/2,temp))
            pvalue = 1 - pvalue
            pvalue = pvalue * ( spinsystem_fermiT(-yp,-ypos1,temp) * spinsystem_fermiT(yp,ypos2,temp)   )

            UTOTAL(i,j) = UTOTAL(i,j) + pvalue * amp
        enddo
        enddo
    end subroutine spinsystem_dodaj_pionowy_slupek_potencjalu


    ! -------------------------------------------------------------------
    ! Oblicza gestosc pradu prawdopodobienstwa na podstawie obecnej f.f.
    ! -------------------------------------------------------------------
    subroutine system_oblicz_J()
          integer          :: i,j,s
          complex*16       :: ycy,cx,jx(-1:1),jy(-1:1),rhoup , rhodwn,Xup,Xdwn
          complex*16       :: jrashba_x , jrashba_y , jlateral_x , jlateral_y
          complex*16       :: tj0


          ! -----------------------------------------------------------------
          tj0 = (-1.0/4.0/DX*II)
          CURRENTXY = 0
          do i = 2 , nx-1
          do j = 2 , ny-1

              jx = 0
              jy = 0
              if(GFLAGS(i,j) == 0) then
                cx  = EXP(+II*DX*DX*(j)*BZ)
                do s = +1 , -1 , -2


                ycy = conjg(cx)*VPHI(GINDEX(i+1,j,s))-VPHI(GINDEX(i-1,j,s))*cx
                jx(s)  =  tj0*(conjg(VPHI(GINDEX(i,j,s)))*ycy - VPHI(GINDEX(i,j,s))*conjg(ycy))

                ycy = VPHI(GINDEX(i,j+1,s))-VPHI(GINDEX(i,j-1,s))
                jy(s)  =  tj0*(conjg(VPHI(GINDEX(i,j,s)))*ycy - VPHI(GINDEX(i,j,s))*conjg(ycy))

                enddo ! end of s

                Xup  = VPHI(GINDEX(i,j,+1))
                Xdwn = VPHI(GINDEX(i,j,-1))

                rhoup  = abs(Xup)**2
                rhodwn = abs(Xdwn)**2

                jrashba_x = II*so_rashba*( Xdwn*conjg(Xup) - conjg(Xdwn)*Xup )
                jrashba_y =    so_rashba*( Xdwn*conjg(Xup) + conjg(Xdwn)*Xup )

                jlateral_x= +(rhoup-rhodwn)*so_loc*DEY(i,j)
                jlateral_y= -(rhoup-rhodwn)*so_loc*DEX(i,j)


                CURRENT(i,j,1)  = CURRENT(i,j,1) + sqrt(jx(+1)**2 + jy(+1)**2)
                CURRENT(i,j,2)  = CURRENT(i,j,2) + sqrt(jx(-1)**2 + jy(-1)**2)
                CURRENT(i,j,3)  = CURRENT(i,j,3) + sqrt(jrashba_x**2  + jrashba_y**2 )
                CURRENT(i,j,4)  = CURRENT(i,j,4) + sqrt(jlateral_x**2 + jlateral_y**2)
                CURRENT(i,j,5)  = CURRENT(i,j,5) + sqrt( (jx(+1)+jx(-1)+jrashba_x+jlateral_x)**2 + &
                                                         (jy(+1)+jy(-1)+jrashba_y+jlateral_y)**2 )


                CURRENTXY(i,j,1) = jx(+1) + jx(-1) + jrashba_x + jlateral_x
                CURRENTXY(i,j,2) = jy(+1) + jy(-1) + jrashba_y + jlateral_y
              endif
          enddo
          enddo
          call system_oblicz_divj()
    endsubroutine system_oblicz_J


    subroutine system_oblicz_divj()
          integer :: i,j
          do i = 3 , nx-2
          do j = 3 , ny-2
              if(GFLAGS(i,j) == 0) then
                DIVJ(i,j) = DIVJ(i,j) +  (CURRENTXY(i+1,j,1)-CURRENTXY(i-1,j,1))/2/DX +&
                                         (CURRENTXY(i,j+1,2)-CURRENTXY(i,j-1,2))/2/DX
              endif
          enddo
          enddo
    end subroutine system_oblicz_divj

    ! ------------------------------------------------------------ -------
    ! Funkcja zapisuje do pliku nazwa dane wskazane przez typ (patrz enum)
    ! ------------------------------------------------------------ -------
    subroutine spinsystem_zapisz_do_pliku(nazwa,typ,xstart,xstop,ystart,ystop)
    character(*) , intent(in) :: nazwa
    integer                   :: typ
    integer,optional,intent(in)::xstart,xstop,ystart,ystop

    double precision :: fval,fval1,fval2,fval_up,fval_dwn
    integer          :: i,j,x1,x2,y1,y2

    x1 = 1
    x2 = Nx
    y1 = 1
    y2 = Ny
    if(present(xstart)) x1 = xstart
    if(present(xstop))  x2 = xstop
    if(present(ystart)) y1 = ystart
    if(present(ystop))  y2 = ystop


    open(unit=86554,file=nazwa)
    select case(typ)
    ! ---------------------------------------
    case(ZAPISZ_KONTUR)
        do i = x1 , x2
        do j = y1 , y2
                fval = 0;
                if(GFLAGS(i,j) >= B_NORMAL) fval = 1
                write(86554,"(3f20.6)"),i*atomic_DX,j*atomic_DX,fval
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case(ZAPISZ_FLAGI)
        do i = x1 , x2
        do j = y1 , y2
                fval = GFLAGS(i,j)
                write(86554,"(3f20.6)"),i*atomic_DX,j*atomic_DX,fval
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case(ZAPISZ_POTENCJAL)
        do i = x1 , x2
        do j = y1 , y2
                fval = UTOTAL(i,j)
                write(86554,"(3f20.6)"),i*atomic_DX,j*atomic_DX,fval
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case(ZAPISZ_INDEKSY)
        do i = x1 , x2
        do j = y1 , y2

                write(86554,"(2f20.6,2i10)"),i*atomic_DX,j*atomic_DX,GINDEX(i,j,+1),GINDEX(i,j,-1)
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case(ZAPISZ_PHI)

        do i = x1 , x2
        do j = y1 , y2
!                fval = abs(PHI(i,j))
                write(86554,"(5f20.6)"),i*atomic_DX,j*atomic_DX,abs(PHI(i,j,+1))+abs(PHI(i,j,-1)),abs(PHI(i,j,+1)),abs(PHI(i,j,-1))
        enddo
            write(86554,*),""
        enddo


    ! ---------------------------------------
    case(ZAPISZ_POLARYZACJE)

        do i = x1 , x2
        do j = y1 , y2
!                fval = abs(PHI(i,j))
                write(86554,"(5f20.6)"),i*atomic_DX,j*atomic_DX,POLARYZACJE(i,j,1),POLARYZACJE(i,j,2),POLARYZACJE(i,j,3)
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case(ZAPISZ_J_TOTAL)

        do i = x1 , x2
        do j = y1 , y2
                fval = CURRENT(i,j,5)
                write(86554,"(5f20.6)"),i*atomic_DX,j*atomic_DX,fval
        enddo
            write(86554,*),""
        enddo
    case(ZAPISZ_DIVJ)

        do i = x1 , x2
        do j = y1 , y2
                fval = DIVJ(i,j)
                write(86554,"(5f20.6)"),i*atomic_DX,j*atomic_DX,fval
        enddo
            write(86554,*),""
        enddo
    case(ZAPISZ_J_ALL)

        write(86554,*),"#x[nm]  y[nm]   j(x,y)  jrashba(x,y)    jlateral(x,y)"
        do i = x1 , x2
        do j = y1 , y2
                fval     =  CURRENT(i,j,5)
                fval_up  =  CURRENT(i,j,1)
                fval_dwn =  CURRENT(i,j,2)
                fval1 =  CURRENT(i,j,3)
                fval2 =  CURRENT(i,j,4)
                write(86554,"(10f20.6)"),i*atomic_DX,j*atomic_DX,fval,fval_up,fval_dwn,fval1,fval2
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case(ZAPISZ_WAVEFUNC)
        print*,"brak wavefunc"
        do i = x1 , x2
        do j = y1 , y2
!                write(86554,"(300e20.6)"),i*atomic_DX,j*atomic_DX,WAVEFUNC(i,j,:)
        enddo
            write(86554,*),""
        enddo
    case(ZAPISZ_STANY_WLASNE)
        print*,"brak stanow wlasnych"
        do i = x1 , x2
        do j = y1 , y2
!                if(GFLAGS(i,j) == B_NORMAL) then
!                    write(86554,"(1000e20.6)"),i*atomic_DX,j*atomic_DX,abs(Widmo_Vecs(GINDEX(i,j),1:Widmo_NoStates))**2
!                else
!                    write(86554,"(1000e20.6)"),i*atomic_DX,j*atomic_DX,abs(Widmo_Vecs(1,1:Widmo_NoStates))*0
!                endif
        enddo
            write(86554,*),""
        enddo
    ! ---------------------------------------
    case default
            print*,"System: Zapisz do pliku - podano zly argument."
    end select

    close(86554)
    endsubroutine spinsystem_zapisz_do_pliku
!
!    ! ------------------------------------------------------------ -------
!    ! Funkcja zapisuje do pliku nazwa dane wskazane przez typ (patrz enum)
!    ! przeznaczona jest dla wynikow otrzymanych za pomoca solvera stanow
!    ! wlasnych.
!    ! stan_start - numer stanu od ktorego bedzie szlo zapisywanie
!    ! stan_end   - koncowy stan
!    ! xstart,xstop,ystart,ystop - zakres rysowania
!    ! ------------------------------------------------------------ -------
!    subroutine system_zapisz_widmo_do_pliku(nazwa,typ,stan_start,stan_end,xstart,xstop,ystart,ystop)
!    character(*) , intent(in) :: nazwa
!    integer                   :: typ
!    integer,optional,intent(in)::stan_start,stan_end,xstart,xstop,ystart,ystop
!    integer          :: i,j,x1,x2,y1,y2,s1,s2
!    doubleprecision  :: zeros(Widmo_NoStates)
!    x1 = 1
!    x2 = Nx
!    y1 = 1
!    y2 = Ny
!    s1 = 1
!    s2 = Widmo_NoStates
!    if(Widmo_NoStates <= 0) then
!        print*,"Error: Nie mozna zapisac stanow wlasnych poniewaz ich liczba wynosi 0.."
!        stop
!    endif
!
!    if(s1 > Widmo_NoStates) s1 = Widmo_NoStates
!    if(s2 > Widmo_NoStates) s2 = Widmo_NoStates
!
!    if(present(xstart)) x1 = xstart
!    if(present(xstop))  x2 = xstop
!    if(present(ystart)) y1 = ystart
!    if(present(ystop))  y2 = ystop
!    if(present(stan_start))  s1 = stan_start
!    if(present(stan_end))    s2 = stan_end
!
!    zeros = 0
!
!
!    open(unit=86554,file=nazwa)
!    select case(typ)
!    ! ---------------------------------------
!    case(ZAPISZ_WIDMO_VRTCAL)
!        do i = s1 , s2
!            write(86554,*),Widmo_Evals(i)*1000.0*Rd
!        enddo
!    case(ZAPISZ_WIDMO_HRZNTL)
!        write(86554,"(3000e20.6)"),DBLE(Widmo_Evals(s1:s2))*1000.0*Rd
!    case(ZAPISZ_STANY_WLASNE)
!        do i = x1 , x2
!        do j = y1 , y2
!                if(GFLAGS(i,j) == B_NORMAL) then
!                    write(86554,"(1000e20.6)"),i*DX*Lr2L,j*DX*Lr2L,abs(Widmo_Vecs(GINDEX(i,j),s1:s2))**2
!                else
!                    write(86554,"(1000e20.6)"),i*DX*Lr2L,j*DX*Lr2L,zeros(s1:s2)
!                endif
!        enddo
!            write(86554,*),""
!        enddo
!    ! ---------------------------------------
!    case default
!            print*,"System: Zapisz widmo do pliku - podano zly argument."
!    end select
!
!    close(86554)
!    endsubroutine system_zapisz_widmo_do_pliku
!
!
!    ! =======================================================================================
!    !
!    !
!    !                          FEAST - STANY WLASNE
!    !
!    ! =======================================================================================
!    ! Znajduje stany wlasne ukladu zamknietego, takiego same dla ktorego odbywa sie
!    ! transport.
!    ! Emin      - minimalna energia od ktorej beda szukane wartosci wlasne [meV]
!    ! Emax      - maksymalna energia [meV]
!    ! NoStates  - spodziewana liczba stanow w tym zakresie
!    ! liczba_konturow[def- 8] - opcjonalna - liczba "controu points" wykonywanych przez feasta.
!    !                    dopuszczalne wartosci to: {3,4,5,6,8,10,12,16,20,24,32,40,48}
!    ! wypisz_informacje[def - 0] - opcjonalna - informuje czy beda wypisywane informacje do konsoli (0 lub 1)
!    ! maks_iter[def - 10] - opcjonalna - maksymalna liczba iteracji wykonywanych przez feasta
!    ! =======================================================================================
!    subroutine system_widmo(pEmin,pEmax,NoStates,pliczba_konturow,pwypisz_informacje,pmaks_iter)
!        doubleprecision   :: pEmin, pEmax
!        integer           :: NoStates
!        integer,optional  :: pliczba_konturow,pwypisz_informacje,pmaks_iter
!
!
!        integer :: fpm(128)
!        integer :: i,j,info,itmp,nw,M0,loop,no_evals,iter
!        integer :: liczba_konturow,wypisz_informacje,maks_iter
!        doubleprecision :: epsout
!        doubleprecision :: Emin, Emax , Ecurr
!
!
!        integer,allocatable                    :: HBROWS(:)
!        complex*16,dimension(:,:), allocatable :: EVectors
!        double precision,dimension(:), allocatable   :: Evalues,Rerrors
!        integer,dimension(:,:),allocatable     :: WINDEX
!        call reset_clock()
!        ! Przejscie do jednostek donorowych
!        Emin = pEmin / 1000.0 / Rd
!        Emax = pEmax / 1000.0 / Rd
!
!        BZ  = BtoDonorB(atomic_Bz)
!        ! Na problem wlasny sa osobne indeksy
!        allocate(WINDEX(nx,ny))
!        WINDEX = 0
!        iter   = 1
!        do i = 1 , nx
!        do j = 1 , ny
!           if( GFLAGS(i,j) == B_NORMAL  ) then
!                WINDEX(i,j) = iter
!                iter = iter + 1
!           endif
!        enddo
!        enddo
!        TRANS_MAXN = iter-1
!
!        ! ustalanie domyslnej liczby konturow i wypisywania
!        if(.not. present(pliczba_konturow)) then
!            liczba_konturow = 8
!        else
!            liczba_konturow = pliczba_konturow
!        endif
!        if(.not. present(pwypisz_informacje)) then
!            wypisz_informacje = 0
!        else
!            wypisz_informacje = pwypisz_informacje
!        endif
!        if(.not. present(pmaks_iter)) then
!            maks_iter = 20
!        else
!            maks_iter = pmaks_iter
!        endif
!
!        call feastinit(fpm)
!
!        fpm(1)=wypisz_informacje ! nie wypisuj informacji
!        fpm(2)=liczba_konturow   ! liczba konturow
!        fpm(3)=12                ! wykladnik bledu ponizej ktorego procedura sie zatrzymuje: e=10^(fpm(3))
!        fpm(4)=maks_iter         ! maksymalna liczba iteracji po ktorej jak sie nie zbiegnie to proced. sie zatrzyma
!        fpm(5)=0                 ! startujemy z domyslnymi wektorami (jak 1 to z dostarczonymi)
!        fpm(6)=0                 ! kryterium zbieznosci poprzez residuum (0 albo 1)
!
!        allocate(CMATA(TRANS_MAXN*5))
!        allocate(IDXA (TRANS_MAXN*5,2))
!        itmp  = 1
!        do i = 2, nx-1
!        do j = 2, ny-1
!        if(WINDEX(i,j) > 0) then
!        if(WINDEX(i,j-1) > 0) then
!                cmatA(itmp) = CMPLX(-0.5/DX/DX)
!                idxA(itmp,1) = WINDEX(i,j)
!                idxA(itmp,2) = WINDEX(i,j-1)
!                itmp = itmp + 1
!        endif
!        if(WINDEX(i-1,j) > 0) then
!
!                cmatA(itmp)   = CMPLX(-0.5/DX/DX*EXP(+II*DX*DX*(j)*BZ)  )
!                idxA (itmp,1) = WINDEX(i  ,j)
!                idxA (itmp,2) = WINDEX(i-1,j)
!                itmp = itmp + 1
!        endif
!        if(WINDEX(i,j) > 0) then
!
!                cmatA(itmp)   = CMPLX( 2.0/DX/DX + UTOTAL(i,j) )
!                idxA (itmp,1) = WINDEX(i,j)
!                idxA (itmp,2) = WINDEX(i,j)
!                itmp = itmp + 1
!        endif
!        if(WINDEX(i+1,j) > 0) then
!
!                cmatA(itmp)   = CMPLX(-0.5/DX/DX*EXP(-II*DX*DX*(j)*BZ)  )
!                idxA (itmp,1) = WINDEX(i  ,j)
!                idxA (itmp,2) = WINDEX(i+1,j)
!                itmp = itmp + 1
!        endif
!        if(WINDEX(i,j+1) > 0) then
!
!                cmatA(itmp) = CMPLX(-0.5/DX/DX)
!                idxA(itmp,1) = WINDEX(i,j)
!                idxA(itmp,2) = WINDEX(i,j+1)
!                itmp = itmp + 1
!        endif
!        endif ! end of jesli jestesmy na komurce B_NORMAl
!        enddo
!        enddo
!
!
!        itmp        = itmp - 1
!        MATASIZE    = itmp
!        nw          = itmp
!        if(wypisz_informacje==1) then
!            print*,"--------------------------------------------------"
!            print*,"Widmo:"
!            print*,"--------------------------------------------------"
!            print*,"Rozmiar problemu N:     ",TRANS_MAXN
!            print*,"Zakres energii:         ",Emin*1000.0*Rd," do ",Emax*1000.0*Rd,"w meV"
!            print*,"Liczba konturow:        ",liczba_konturow
!            print*,"Maksymalna liczba iter. :",maks_iter
!            print*,"Zalozona liczba stanow: ",NoStates
!        endif
!
!        ! --------------------------------------------system_inicjalizacja_ukladu---------------
!        !
!        ! -----------------------------------------------------------
!        allocate(HBROWS(TRANS_MAXN+1))
!        call convert_to_HB(MATASIZE,IDXA,HBROWS)
!
!
!        ! zgadujemy liczbe stanow
!        M0  = NoStates
!
!        allocate(EVectors(TRANS_MAXN,M0))
!        allocate(Evalues(M0))
!        allocate(Rerrors(M0))
!
!
!
!        call zfeast_hcsrev('F',&                ! - 'F' oznacza ze podawana jest pelna macierz
!                              TRANS_MAXN,&    ! - rozmiar problemu (ile wezlow z flaga B_NORMAL)
!                              CMATA(1:nw),&    ! - kolejne nie zerowe wartosci w macierzy H
!                              HBROWS,&         ! - numeracja wierszy (rodzaj zapisu macierzy rzakidch)
!                              idxA(1:nw,2),&   ! - indeksy kolumn odpowiadaja tablicy wartosci CMATA
!                              fpm,&            ! - wektor z konfiguracja procedury
!                              epsout,&         ! - Residuum wyjsciowe
!                              loop, &          ! - Koncowa liczba iteracji
!                              Emin,&           ! - Minimalna energia przeszukiwania
!                              Emax,&           ! - Maksymalna energia
!                              M0,&             ! - Spodziewana liczba modow w zakresie (Emin,Emax)
!                              Evalues,&        ! - Wektor z otrzymanymi wartosciami wlasnymi
!                              EVectors,&       ! - Macierz z wektorami (kolejne kolumny odpowiadaja kolejnym wartoscia z tablicy Evalues)
!                              no_evals,&       ! - Liczba otrzymanych wartosci z przedziale (Emin,Emax)
!                              Rerrors,&        ! - Wektor z bledami dla kolejnych wartosci wlasnych
!                              info)            ! - Ewentualne informacje o bledach
!
!
!
!        if(wypisz_informacje==1) then
!            print*,"Eps wyjsciowy           :",  epsout
!            print*,"Liczba iteracji         :",  loop
!            print*,"Znaleziona l. stanow    :",  no_evals
!            print*,"Info                    :",  info
!            print*,"Czas obliczen [s]       :",  get_clock()
!            print*,"--------------------------------------------------"
!        endif
!
!        Widmo_NoStates = no_evals
!
!        ! ----------------------------------------------------------------------------------
!        ! Obsluga bledow:
!        ! ----------------------------------------------------------------------------------
!        selectcase(info)
!        case( 202 )
!            print*," Error : Problem with size of the system n (n≤0) "
!            stop
!        case( 201 )
!            print*," Error : Problem with size of initial subspace m0 (m0≤0 or m0>n) "
!            stop
!        case( 200 )
!            print*," Error : Problem with emin,emax (emin≥emax) "
!            stop
!        case(100:199)
!            print"(A,I4,A)"," Error : Problem with ",info-100,"-th value of the input Extended Eigensolver parameter (fpm(i)). Only the parameters in use are checked. "
!            Widmo_NoStates = 0
!            stop
!        case( 4 )
!            print*," Warning : Successful return of only the computed subspace after call withfpm(14) = 1 "
!            Widmo_NoStates = 0
!
!        case( 3 )
!            print*," Warning : Size of the subspace m0 is too small (m0<m) "
!            Widmo_NoStates = 0
!
!        case( 2 )
!            print*," Warning : No Convergence (number of iteration loops >fpm(4))"
!            Widmo_NoStates = 0
!        case( 1 )
!            print*," Warning : No eigenvalue found in the search interval. See remark below for further details. "
!            Widmo_NoStates = 0
!        case( 0 )
!            print*,               "---------------------------------------------"
!            print"(A,i12)",       "Widmo: Znaleziono stanow :",Widmo_NoStates
!            print"(A,f12.3)",     "       W czasie T [s]    :",get_clock()
!            print"(A,e12.4)",     "       Z bledem epsout   :",epsout
!            print"(A,i12)",       "       W po liczbie iter.:",loop
!            print*,               "---------------------------------------------"
!        case( -1 )
!            print*," Error : Internal error for allocation memory. "
!            stop
!        case( -2 )
!            print*," Error : Internal error of the inner system solver. Possible reasons: not enough memory for inner linear system solver or inconsistent input. "
!            stop
!        case( -3 )
!            print*," Error : Internal error of the reduced eigenvalue solver Possible cause: matrix B may not be positive definite. It can be checked with LAPACK routines, if necessary."
!            stop
!        case(-199:-100)
!            print"(A,I4,A)"," Error : Problem with the ",-info-100,"-th argument of the Extended Eigensolver interface. "
!            stop
!        endselect
!
!        ! -----------------------------------------------------------------
!        ! Kopiowanie wynikow do odpowiednich tablic
!        ! -----------------------------------------------------------------
!
!
!        if(allocated(Widmo_Evals)) deallocate(Widmo_Evals)
!        if(allocated(Widmo_Vecs))  deallocate(Widmo_Vecs)
!
!
!
!!        Widmo_NoStates = no_evals
!        if(Widmo_NoStates > 0 ) then
!            call oblicz_rozmiar_macierzy();
!
!            allocate(Widmo_Vecs(TRANS_MAXN,Widmo_NoStates))
!            do i = 1 , nx
!            do j = 1 , ny
!                if(WINDEX(i,j) > 0) then
!                    Widmo_Vecs(GINDEX(i,j),1:Widmo_NoStates) = EVectors(WINDEX(i,j),1:Widmo_NoStates)
!                endif
!            enddo
!            enddo
!
!            allocate(Widmo_Evals(Widmo_NoStates))
!            Widmo_Evals(1:Widmo_NoStates) = Evalues(1:Widmo_NoStates)
!        endif
!
!        deallocate(CMATA)
!        deallocate(IDXA)
!        deallocate(HBROWS)
!        deallocate(EVectors)
!        deallocate(Evalues)
!        deallocate(Rerrors)
!        deallocate(WINDEX)
!
!    end subroutine system_widmo
!
!
!
!
    ! ==========================================================================
    !
    !
    !                          SUPER LU
    !
    !
    !
    ! ==========================================================================

    subroutine convert_to_HB(no_vals,rows_cols,out_rows)
          integer,intent(in)                      :: no_vals
          integer,intent(inout),dimension(:,:)  :: rows_cols
          integer,intent(inout),dimension(:)    :: out_rows
          integer :: iterator, irow , from , to
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
          if(TRANS_SOLVER == USE_UMFPACK) then
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

          endif
!DEC$ ENDIF



      end subroutine convert_to_HB

    subroutine sort_col_vals(cols,vals)
            integer,intent(inout),dimension(:)    :: cols
            complex*16,intent(inout),dimension(:) :: vals
            integer :: tmp_col
            complex*16 :: tmp_val
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
        complex*16,intent(in),dimension(:) :: values
        complex*16,intent(inout),dimension(:) :: b
        integer :: iopt
        integer n, nnz, nrhs, ldb

        integer, save    ::  info = 0
        integer*8 , save :: factors = 0
        complex*16,allocatable,dimension(:),save :: b_sol
        doubleprecision,save :: total_time

!DEC$ IF DEFINED  (USE_UMF_PACK)
        ! UMFPACK constants
        type(c_ptr),save :: symbolic,numeric
        ! zero-based arrays
        real(8),save :: control(0:UMFPACK_CONTROL-1),umf_info(0:UMFPACK_INFO-1)
!DEC$ ENDIF


        n    = no_rows
        nnz  = no_vals
        ldb  = n
        nrhs = 1


!DEC$ IF DEFINED  (USE_UMF_PACK)
    selectcase (iopt)
      case (1)
            allocate(b_sol(size(b)))

            call umf4cdef (control)
            call umf4csym (n,n, rowind, colptr, values, symbolic, control, umf_info)
            call umf4cnum (rowind, colptr, values, symbolic, numeric, control, umf_info)
            call umf4cfsym (symbolic)
            total_time =  umf_info(UMFPACK_NUMERIC_TIME)+umf_info(UMFPACK_SYMBOLIC_TIME)
            if (umf_info(UMFPACK_STATUS) .eq. 0) then
                if(TRANS_DEBUG) then
                     write (*,*) 'Factorization succeeded. Time needed:',&
                    total_time, " Mem needed:", umf_info(UMFPACK_PEAK_MEMORY)/8.0/1024/1024 , "[MB]"
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
                total_time  = total_time + umf_info(UMFPACK_SOLVE_WALLTIME)
            else
                 write(*,*) 'UMF ERROR: INFO from solve = ', umf_info(UMFPACK_STATUS)
            endif

      case(3)
            print*,"UMFPACK Solved:"
            print*,"      Total time needed:",total_time,"[s]"
            !print*,"      S&N memory needed:",umf_info(UMFPACK_PEAK_MEMORY)/1024/1024,"[MB]"
            call umf4cfnum (numeric)
            deallocate(b_sol)
      endselect
!DEC$ ELSE
    selectcase (iopt)
      case (1)
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
      endselect
!DEC$ ENDIF


      endsubroutine solve_system

end module modspinsystem

