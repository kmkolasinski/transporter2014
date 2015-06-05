module modspinzrodlo
    use modutils
    use modjed
    use spinmodpop
    use ifport
    use modinip
implicit none




    ! =========================================================
    !                   Klasa zrodla
    ! =========================================================
    ! przechowuje informacje o zrodle
    type cspinzrodlo

        integer :: N ! liczba oczek siatki
        integer :: bKierunek      ! ustala zwrot zrodla, 1 - w kierunku dodatnich wartosci x lub y, 0 - przeciwnie
        double precision :: hnY , hnX   ! wysokosc zrodla w oczkach siatki, moze byc polowkowa
        integer :: dir ! orientacja zrodla (-1,+1)
        double precision,dimension(2) :: r1  ! nie uzywane w programie polozenie pierwszego punktu w nm
        double precision,dimension(2) :: r2  ! nie uzywane w programie polozenie ostatniego punktu w nm

        integer,dimension(:,:),allocatable :: polozenia; ! wekor polozen (n,1) ,(n,2) na siatce numerycznej mapowanie na indeks
                                                         ! bedzie odbywac sie przez tablice GINDEX(i,j)
        complex*16,dimension(:,:,:),allocatable   :: ChiMod      ! macierz modow wchodzacych, (N,M,K) - M liczba modow , K - kierunek (1 - in , -1 - out)
        complex*16,dimension(:,:,:),allocatable   :: ChiModUp    ! macierz modow wchodzacych spin up, (N,M,K) - M liczba modow , K - kierunek (1 - in , -1 - out)
        complex*16,dimension(:,:,:),allocatable   :: ChiModDown  ! macierz modow wchodzacych spin down, (N,M,K) - M liczba modow , K - kierunek (1 - in , -1 - out)
        complex*16,dimension(:,:),allocatable     :: ChiKvec     ! wektory falowe zwiazane z modami wchodz (M,K) - M - mod, K- kierunek (+1,-1).
        complex*16,dimension(:,:),allocatable     :: ChiLambda   ! lamda zwiazana z danym modem = exp(k*dx) , gdzie k moze byc zespolone
        doubleprecision,dimension(:,:),allocatable:: ChiCurr     ! prady wejsciowe (M-mod kierunek, K - kierunek (+1,-1))
        doubleprecision,dimension(:,:),allocatable:: TR          ! (mod,-1-R,+1-T) - Prawdopodobienstwo przejscia/odbicia

        doubleprecision,dimension(:),allocatable  :: Uvec
        doubleprecision,dimension(:,:),allocatable  :: SpinUVec
        doubleprecision,dimension(:),allocatable  :: Ey
!        complex*16,dimension(:),allocatable       :: deltamk  ! tablica pomocnicza przyspieszajaca liczenie elementow macierzowych
!        complex*16,dimension(:),allocatable       :: deltaink  ! tablica pomocnicza przyspieszajaca liczenie elementow macierzowych
!        complex*16,dimension(:),allocatable       :: deltaoutk  ! tablica pomocnicza przyspieszajaca liczenie elementow macierzowych

        integer                                   :: liczba_modow ! liczba dostepnym w zrodle modow
        integer                                   :: liczba_evans ! liczba modow evanescentnych

        complex*16,dimension(:),allocatable       :: ck       ! amplitudy wchodzace
        complex*16,dimension(:),allocatable       :: dk       ! amplitudy odbite


!        complex*16,dimension(200,200)             :: m_r,m_t  ! macierz wspolczynnikow odbicia i transmisji
        complex*16,dimension(:,:),allocatable     :: Aij
        complex*16,dimension(:,:),allocatable     :: Sij
        complex*16,dimension(:,:,:,:),allocatable :: Sigma ! (v,j,spin1,spin2) , spiny (-1,+1)
        complex*16,dimension(:,:,:),allocatable   :: SijChiAuxMat ! macierz pomocnicza zavierajaca iloczyny macierzy Sij i Wektorow Chi
        complex*16,dimension(:),allocatable       :: SijAijCkAuxVec ! pomocniczy wektor z obliczonymi iloczynami
        complex*16,dimension(:,:),allocatable     :: Fj ! wektor wyrazow wolnych (i,spin), gdzie spin -1,+1

        logical :: bZaalokowane = .false. ! flaga ktora mowi o tym czy zrodlo ma juz zaalokowana pamiec
        contains
        ! -------------------------------------------------
        ! Procedury zrodla
        ! -------------------------------------------------
        procedure, public, pass(zrodlo) :: spinzrodlo_zwolnij_pamiec!() wiadomo
        procedure, public, pass(zrodlo) :: spinzrodlo_alokuj_pamiec !(pN,lM) podajemy liczbe punktow oraz liczbe modow
        procedure, public, pass(zrodlo) :: spinzrodlo_wypisz_info   !() wypisuje parametry zrodla
        procedure, public, pass(zrodlo) :: spinzrodlo_wypisz_ckdk   !() wypisyje amplitudy ck, dk
        procedure, public, pass(zrodlo) :: spinzrodlo_wypisz_JinJout!() wypisyje strumien wejsciowy oraz wyjsciowy (ich stosunki)
        procedure, public, pass(zrodlo) :: spinzrodlo_zapisz_mody!(zrodlo,up_filename,down_filename,dx,bSaveEvan)
!        procedure, public, pass(zrodlo) :: zrodlo_zapisz_m_rt   !(mod) zapisuje biezace amplitudy od macierzy m_r i m_t dla zadanego modu
        procedure, public, pass(zrodlo) :: spinzrodlo_ustaw   !(mod) zapisuje biezace amplitudy od macierzy m_r i m_t dla zadanego modu
        procedure, public, pass(zrodlo) :: spinzrodlo_relacja_dyspersji!(zrodlo,pdx,pEf,pBz,pkmin,pkmax,pdk,pEmax,nazwa_pliku)
!        procedure, public, pass(zrodlo) :: zrodlo_alfa_v_i   !
        procedure, public, pass(zrodlo) :: spinzrodlo_oblicz_Fj   !
        procedure, public, pass(zrodlo) :: spinzrodlo_oblicz_dk!(zrodlo,VPHI,GINDEX,nx,ny)   !
        procedure, public, pass(zrodlo) :: spinzrodlo_oblicz_JinJout!(zrodlo)

    end type cspinzrodlo

!    ! =========================================================
!    !                   Klasa abszrodlo
!    ! wykorzystywana jest to symulacji transparentnych w.b zaproponowanych
!    ! przez nowaka
!    ! =========================================================
!    type cabs_zrodlo
!        integer :: N ! liczba oczek siatki
!        integer :: bKierunek      ! ustala zwrot zrodla, 1 - w kierunku dodatnich wartosci x lub y, 0 - przeciwnie
!        integer,dimension(:,:),allocatable :: polozenia; ! wekor polozen (n,1) ,(n,2) na siatce numerycznej mapowanie na indeks
!        doubleprecision :: kvec ! wirtualny wektor falowy
!        contains
!        procedure, public, pass(zrodlo) :: abs_zrodlo_zwolnij_pamiec
!        procedure, public, pass(zrodlo) :: abs_zrodlo_ustaw
!        procedure, public, pass(zrodlo) :: abs_zrodlo_skopiuj!(zrodlo,od_zrodla)
!    end type cabs_zrodlo


contains

!------------------------------------------------------------
!                 Funkcje klasy zrodla
!------------------------------------------------------------
    subroutine spinzrodlo_zwolnij_pamiec(zrodlo)
        class(cspinzrodlo) :: zrodlo

        if(TRANS_DEBUG==.true.) print*,"Spinzrodlo: Zwalanianie pamieci"
        if(allocated(zrodlo%ChiKvec))    deallocate(zrodlo%ChiKvec);
        if(allocated(zrodlo%Uvec))       deallocate(zrodlo%Uvec);
        if(allocated(zrodlo%SpinUVec))   deallocate(zrodlo%SpinUVec);
        if(allocated(zrodlo%Ey))         deallocate(zrodlo%Ey);
        if(allocated(zrodlo%ChiLambda))  deallocate(zrodlo%ChiLambda);
        if(allocated(zrodlo%ChiMod))     deallocate(zrodlo%ChiMod);
        if(allocated(zrodlo%ChiModUp))   deallocate(zrodlo%ChiModUp);
        if(allocated(zrodlo%ChiModDown)) deallocate(zrodlo%ChiModDown);
        if(allocated(zrodlo%ChiCurr))    deallocate(zrodlo%ChiCurr);
        if(allocated(zrodlo%polozenia))  deallocate(zrodlo%polozenia);
        if(allocated(zrodlo%ck))         deallocate(zrodlo%ck);
        if(allocated(zrodlo%dk))         deallocate(zrodlo%dk);
        if(allocated(zrodlo%Aij))        deallocate(zrodlo%Aij);
        if(allocated(zrodlo%Sij))        deallocate(zrodlo%Sij);
        if(allocated(zrodlo%Sigma))      deallocate(zrodlo%Sigma);
        if(allocated(zrodlo%Fj))         deallocate(zrodlo%Fj);
        if(allocated(zrodlo%SijAijCkAuxVec))deallocate(zrodlo%SijAijCkAuxVec);
        if(allocated(zrodlo%SijChiAuxMat))  deallocate(zrodlo%SijChiAuxMat);
        if(allocated(zrodlo%TR))         deallocate(zrodlo%TR);

        zrodlo%bZaalokowane = .false.

    end subroutine spinzrodlo_zwolnij_pamiec


! ----------------------------------------------------------
! Alokuje zrodlo dla zadanej liczby oczek siatki - pN
! oraz lM - liczby modow oraz liczby modow evanescentnych - lEvanMods
! ----------------------------------------------------------
    subroutine spinzrodlo_alokuj_pamiec(zrodlo,pN,lM,lEvanMods)
        class(cspinzrodlo) :: zrodlo
        integer,intent(in) :: pN,lM,lEvanMods

        zrodlo%N            = pN;
        zrodlo%liczba_modow = lM;
        zrodlo%liczba_evans = lEvanMods;

        if(TRANS_DEBUG==.true.) print*,"Spinzrodlo: Alokowanie pamieci zrodla spinorowego"

        if(zrodlo%bZaalokowane == .true.) then
            call zrodlo%spinzrodlo_zwolnij_pamiec()
        endif

        allocate(zrodlo%polozenia (pN  ,2) );
        allocate(zrodlo%Uvec      (pN  ) );
        allocate(zrodlo%SpinUVec  (pN ,-1:1 ) );
        allocate(zrodlo%Ey        (pN  ) );
        allocate(zrodlo%ChiMod    (2*pN,lM+lEvanMods,-1:1));   ! 2 * pN bo spinor
        allocate(zrodlo%ChiModUp  (pN  ,lM+lEvanMods,-1:1));   ! 1 * pN bo up
        allocate(zrodlo%ChiModDown(pN  ,lM+lEvanMods,-1:1));   ! 1 * pN bo down

        allocate(zrodlo%ChiCurr   (lM+lEvanMods,-1:1));
        allocate(zrodlo%ChiKvec   (lM+lEvanMods,-1:1));
        allocate(zrodlo%ChiLambda (lM+lEvanMods,-1:1));
        allocate(zrodlo%TR        (lM,-1:1));
        allocate(zrodlo%ck        (lM+lEvanMods));
        allocate(zrodlo%dk        (lM+lEvanMods));


        allocate(zrodlo%Sij       (lM+lEvanMods,lM+lEvanMods));
        allocate(zrodlo%Aij       (lM+lEvanMods,lM));
        allocate(zrodlo%Sigma     (pN,pN,-1:+1,-1:+1));
        allocate(zrodlo%Fj        (pN,-1:+1));
        allocate(zrodlo%SijAijCkAuxVec (lM+lEvanMods));
        allocate(zrodlo%SijChiAuxMat (lM+lEvanMods,pN,-1:+1));

        zrodlo%ChiCurr    = 0
        zrodlo%ChiKvec    = 0
        zrodlo%ChiLambda  = 0
        zrodlo%ChiMod     = 0
        zrodlo%ChiModUp   = 0
        zrodlo%ChiModDown = 0
        zrodlo%ck         = 0
        zrodlo%dk         = 0
        zrodlo%Uvec       = 0
        zrodlo%SpinUVec   = 0
        zrodlo%Sij        = 0
        zrodlo%Aij        = 0
        zrodlo%Sigma      = 0
        zrodlo%Fj         = 0
        zrodlo%bZaalokowane = .true.
    end subroutine spinzrodlo_alokuj_pamiec

    subroutine spinzrodlo_wypisz_info(zrodlo)
        class(cspinzrodlo) :: zrodlo
        integer :: i
        double precision :: polarA, polarB
        if(TRANS_DEBUG==.true.) then
        print*,"Spinzrodlo: Wypisywanie informacji o zrodle:"
        print*,"    N       = ",zrodlo%N
        print*,"    L. modow= ",zrodlo%liczba_modow
        print*,"    L. evanm= ",zrodlo%liczba_evans

        call modjed_jaki_kierunek(zrodlo%bKierunek)
        print*,"    HNX     = ",zrodlo%hnX
        print*,"    HNY     = ",zrodlo%hnY
        write(*,"(A,2f10.4,A)"),"    R1(x,y) = (",zrodlo%r1*Lr2L,")"
        write(*,"(A,2f10.4,A)"),"    R2(x,y) = (",zrodlo%r2*Lr2L,")"
        print*,"    Mody wejsciowe/wyjsciowe:",zrodlo%liczba_modow," (polaryzacja) "


        do i = 1 , zrodlo%liczba_modow

             polarA = sum(abs(zrodlo%ChiModUp(:,i,+1))**2) - sum(abs(zrodlo%ChiModDown(:,i,1))**2)
             polarA = polarA/(sum(abs(zrodlo%ChiModUp(:,i,+1))**2) + sum(abs(zrodlo%ChiModDown(:,i,+1))**2))

             polarB = sum(abs(zrodlo%ChiModUp(:,i,-1))**2) - sum(abs(zrodlo%ChiModDown(:,i,-1))**2)
             polarB = polarB/(sum(abs(zrodlo%ChiModUp(:,i,-1))**2) + sum(abs(zrodlo%ChiModDown(:,i,-1))**2))

             print"(A,i4,A,2f8.4,A,f8.4,A,2f8.4,A,f8.4,A)","K[",i,"][nm]:",&
                    zrodlo%ChiKvec(i,+1)*L2LR," p(",dble(polarA),") | ",&
                    zrodlo%ChiKvec(i,-1)*L2LR," p(",dble(polarA),")"
        enddo

        print*,"    Mody wejsciowe/wyjsciowe evanescentne:",zrodlo%liczba_evans!, "BRAK OBSLUGI"
!        do i = zrodlo%liczba_modow + 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
!            print"(A,i4,A,2f12.6,A,2f12.6)","K[",i,"][nm]:",zrodlo%ChiKvec(i,+1)*L2LR," | ",zrodlo%ChiKvec(i,-1)*L2LR
!        enddo
        endif
    endsubroutine spinzrodlo_wypisz_info



    subroutine spinzrodlo_wypisz_ckdk(zrodlo)
        class(cspinzrodlo) :: zrodlo
        integer :: i
        if(TRANS_DEBUG==.true.) then
        print*, "! ----------------------------------------- !"
        if(zrodlo%bKierunek == ZRODLO_KIERUNEK_PRAWO .or. zrodlo%bKierunek == ZRODLO_KIERUNEK_GORA ) then
            print*, "       ------>       |        <------"
        else
            print*, "       <------       |        ------>"
        endif
        print*, "! ----------------------------------------- !"
        do i = 1 , zrodlo%liczba_modow
            write(*,"(A,i4,A,f8.4,A,i4,A,f8.4)"),"   Ck(",i,")=",abs(zrodlo%ck(i))**2,"  |  Dk(",i,")=",abs(zrodlo%dk(i))**2
        enddo
        endif
    endsubroutine spinzrodlo_wypisz_ckdk


    subroutine spinzrodlo_wypisz_JinJout(zrodlo)
        class(cspinzrodlo) :: zrodlo
        integer :: i
        if(TRANS_DEBUG==.true.) then
        print*, "! ----------------------------------------- !"
        if(zrodlo%bKierunek == ZRODLO_KIERUNEK_PRAWO .or. zrodlo%bKierunek == ZRODLO_KIERUNEK_GORA  ) then
            print*, "        ------>       |        <------"
        else
            print*, "        <------       |        ------>"
        endif
        print*, "! ----------------------------------------- !"
        do i = 1 , zrodlo%liczba_modow
            write(*,"(A,i4,A,f8.4,A,i4,A,f8.4)"),"   Jin(",i,")=",zrodlo%ChiCurr(i,+1),"  |  Jout(",i,")=",zrodlo%ChiCurr(i,-1)
        enddo
        print*, "! ----------------------------------------- !"
        write(*,"(A,f8.4,A,f8.4)"),"    SUMA(in)=",sum(zrodlo%ChiCurr(:,+1)),"  |   SUMA(out)=",sum(zrodlo%ChiCurr(:,-1))
        endif
    endsubroutine spinzrodlo_wypisz_JinJout
!
!    subroutine spinzrodlo_zapisz_m_rt(zrodlo,mod)
!        class(cspinzrodlo) :: zrodlo
!        integer,intent(in) :: mod ! numer modu wchodzacego
!        integer            :: i
!        do i = 1 , zrodlo%liczba_modow
!            zrodlo%m_r(mod,i) = zrodlo%dk(i)
!            zrodlo%m_t(mod,i) = zrodlo%ck(i)
!        enddo
!    endsubroutine spinzrodlo_zapisz_m_rt
!
    ! -------------------------------------------------------------------------
    ! Zapisuje mody do pliku o nazwie filename, podajemy rozwniez krok siatki w nm
    ! bSaveEvan - optionalnie mozemy zdecydowac czy chcemy zapisywac evanescente czy tylko
    ! fale plaskie, domyslnie mody evanescentne nie sa zapisywane
    ! -------------------------------------------------------------------------
    subroutine spinzrodlo_zapisz_mody(zrodlo,up_filename,down_filename,bSaveEvan)
        class(cspinzrodlo) :: zrodlo
        character(*) :: up_filename,down_filename
        logical,optional,intent(in) :: bSaveEvan
        integer :: i,no_modow
        logical :: bSaveE

        print*,"Spinzrodlo: Zapisywanie modow do pliku:",up_filename," oraz ",down_filename
        open(4193,file=up_filename)
        open(4194,file=down_filename)
        if(present(bSaveEvan)) then
            bSaveE = bSaveEvan ! domyslnie nie bedziemy zapisywac modow evanescentym
        else
            bSaveE = .false.
        endif
        no_modow = zrodlo%liczba_modow
        ! spin up
        do i = 1 , zrodlo%N
        if(bSaveE)           write(4193,"(500f20.8)"),zrodlo%polozenia(i,:)*atomic_DX,zrodlo%Uvec(i),zrodlo%SpinUVec(i,1),zrodlo%SpinUVec(i,-1),abs(zrodlo%ChiModUp(i,:,+1))**2,abs(zrodlo%ChiModUp(i,:,-1))**2
        if(.not. bSaveE)     write(4193,"(500f20.8)"),zrodlo%polozenia(i,:)*atomic_DX,zrodlo%Uvec(i),zrodlo%SpinUVec(i,1),zrodlo%SpinUVec(i,-1),abs(zrodlo%ChiModUp(i,1:no_modow,+1))**2,abs(zrodlo%ChiModUp(i,1:no_modow,-1))**2
        enddo
        ! spin down
        do i = 1 , zrodlo%N
        if(bSaveE)           write(4194,"(500f20.8)"),zrodlo%polozenia(i,:)*atomic_DX,zrodlo%Uvec(i),zrodlo%SpinUVec(i,1),zrodlo%SpinUVec(i,-1),abs(zrodlo%ChiModDown(i,:,+1))**2,abs(zrodlo%ChiModDown(i,:,-1))**2
        if(.not. bSaveE)     write(4194,"(500f20.8)"),zrodlo%polozenia(i,:)*atomic_DX,zrodlo%Uvec(i),zrodlo%SpinUVec(i,1),zrodlo%SpinUVec(i,-1),abs(zrodlo%ChiModDown(i,1:no_modow,+1))**2,abs(zrodlo%ChiModDown(i,1:no_modow,-1))**2
        enddo

        close(4193)
        close(4194)
    end subroutine spinzrodlo_zapisz_mody

    ! --------------------------------------------------------------------
    ! Obliczanie relacji dyspersji dla podanego zrodla:
    ! Parametry:
    ! pkmin,pkmax - podajemy minimalna wartosc wektora k w jednostkach 1/nm
    ! pdk         - rozdzielczosc skanu
    ! pEmax       - maksymalna energia do ktorej liczone sa stany. Zwykle > Ef.
    ! nazwa_pliku - nazwa pliku do jakiego zostanie zapisana relacja dyspersji
    ! --------------------------------------------------------------------
    subroutine spinzrodlo_relacja_dyspersji(zrodlo,pkmin,pkmax,pdk,pEmax,nazwa_pliku)
        class(cspinzrodlo)   ::  zrodlo
        double precision     :: pdx,pEf,pBz,pkmin,pkmax,pdk,pEmax
        character(*)         :: nazwa_pliku
        integer :: i
        print*,"Spinzrodlo: Zapisywanie relacji dyspresji do pliku:",nazwa_pliku

        if(allocated(SUVEC)) deallocate(SUVEC)
        allocate(SUVEC(zrodlo%N,-1:1))
        do i = 1 , zrodlo%N
            SUVEC(i,+1) = SUTOTAL(zrodlo%polozenia(i,1),zrodlo%polozenia(i,2),+1) ! w jednostkach donorowych
            SUVEC(i,-1) = SUTOTAL(zrodlo%polozenia(i,1),zrodlo%polozenia(i,2),-1) ! w jednostkach donorowych
        enddo


        call spinmodpop_inicjalizacja(atomic_DX,zrodlo%N,atomic_Ef,atomic_Bz,zrodlo%Uvec)
        call spinmodpop_relacja_dyspersji(pkmin,pkmax,pdk,pEmax,nazwa_pliku)
        call spinmodpop_zwalnienie_pamieci()
        deallocate(SUVEC)

    end subroutine spinzrodlo_relacja_dyspersji

    ! --------------------------------------------------------------------
    ! Ustawiamy wczesniej zaalokowane zrodlo przez polecenie systemowe.
    ! Parametry:
    ! pY1,pYN - polozenia y na siatce numerycznej liczone od 1 do NY, w przypadku zrodla poziomego
    !           oznaczaja polozenia X1 oraz XN
    ! pX1 - polozenie X zrodla, dla zrodla poziomego polozenie Y1
    ! pKierunek - enum ZRODLO_KIERUNEK_PRAWO/LEWO/GORA/DOL - ustala w ktora skierowane jest zrodlo
    ! pUTOTAL   - referencja do potencjalu ukladu
    ! --------------------------------------------------------------------
    subroutine spinzrodlo_ustaw(zrodlo,pY1,pYN,pX1,pKierunek,pUTOTAL)
        class(cspinzrodlo)         ::  zrodlo
        integer,intent(in)         ::  pY1,pYN,pX1
        integer,intent(in)         ::  pKierunek ! enum
        double precision,dimension(:,:) :: pUTOTAL ! calkowity potencjal w [meV]
        double precision :: pUvec(pYN - pY1 + 1)

        double precision :: dx, Ef, Bz
        complex*16,dimension(:,:),allocatable :: tempB
        complex*16 :: ctmp , lambda , lambda0
        integer :: N,N2N ! liczba oczek dla zrodla
        integer :: lModow,lModowEvan,dir
        ! zmienne pomocnicze
        integer :: i,j,k,ntmp,v,p

        dx  = atomic_DX*L2LR    ! konwertujemy do jednostek donorowych
        Ef  = atomic_Ef/1000.0/Rd
        BZ  = BtoDonorB(atomic_Bz)
        call reset_clock()
        print*,"Spinzrodlo: Tworzenie nowego zrodla."
        print*,"    Maksymalna wartosc Bz",DonorBtoB(2*3.14159/dx/dx) , "[T]"

        N = pYN - pY1 + 1; ! to jest tak samo liczone ale zmienne pYN , pY1 maja inna interpretacje



        if( pKierunek == ZRODLO_KIERUNEK_PRAWO .or. pKierunek == ZRODLO_KIERUNEK_LEWO ) then

            if(allocated(SUVEC)) deallocate(SUVEC)
            allocate(SUVEC(N,-1:1))
            do i = 1 , N
                pUVEC(i) = pUTOTAL(pX1,pY1 + i - 1)
                SUVEC(i,+1) = SUTOTAL(pX1,pY1 + i - 1,+1) ! w jednostkach donorowych
                SUVEC(i,-1) = SUTOTAL(pX1,pY1 + i - 1,-1) ! w jednostkach donorowych
            enddo


            call spinmodpop_calc_modes_from_wfm(atomic_DX,N,atomic_Ef,atomic_Bz,pUVEC,.true.)
        else ! dla zrodel gora dol
            do i = 1 , N
                pUVEC(i) = pUTOTAL(pY1 + i - 1,pX1)
            enddo
            call spinmodpop_calc_modes_from_wfm(atomic_DX,N,atomic_Ef,atomic_Bz,pUVEC,.false.)
        endif
        deallocate(SUVEC)


        lModow     = LICZBA_MODOW
        lModowEvan = LICZBA_MODOW_EVANESCENTYCH

        print*,"    Liczba modow     =",lModow
        print*,"    Liczba modow evan=",lModowEvan

!        print*,"Szukanie modow A = ", get_clock()
        call zrodlo%spinzrodlo_alokuj_pamiec(N,lModow,lModowEvan)

        zrodlo%Uvec = pUVEC ! w jednostkach atomowych


        ! pole elekrtyczne jest trzymane w jednostkach donorowych
        do i = 2 , N-1
            zrodlo%Ey(i) = (zrodlo%Uvec(i+1) - zrodlo%Uvec(i-1))/2/DX/1000.0/Rd
        enddo
        zrodlo%Ey(1) = zrodlo%Ey(2)
        zrodlo%Ey(N) = zrodlo%Ey(N-1)


        if(lModow > 0) then
            ! POBIERAMY MODY POPRZECZNE ORAZ ICH WEKTORY FALOWE
            ! pobieramy od razu z evanescentnymi wektorami
            zrodlo%ChiKvec  (1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,+1) = ChiKvec(1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,1)
            zrodlo%ChiKvec  (1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,-1) = ChiKvec(1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,2)

            zrodlo%ChiLambda(1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,+1) = ChiLambda(1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,1)
            zrodlo%ChiLambda(1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,-1) = ChiLambda(1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,2)

            zrodlo%ChiMod(:,1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,+1) = ChiMod(:,1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,1)
            zrodlo%ChiMod(:,1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,-1) = ChiMod(:,1:LICZBA_MODOW+LICZBA_MODOW_EVANESCENTYCH,2)

            zrodlo%ChiModUp  (1:N,:,:) = zrodlo%ChiMod(1:N,:,:)
            zrodlo%ChiModDown(1:N,:,:) = zrodlo%ChiMod(N+1:2*N,:,:)
        else
            print*,"    Liczba modow wynosi zero - transportu brak"

        endif

        ! MODUL RELACJI DYSPERSJI JUZ NIE POTRZEBNY
        call spinmodpop_zwalnienie_pamieci()
!        print*,"Szukanie modow=", get_clock()
        call reset_clock()

        ! ustawianie parametrow zrodel
        zrodlo%bKierunek = pKierunek

        if( pKierunek == ZRODLO_KIERUNEK_PRAWO .or. pKierunek == ZRODLO_KIERUNEK_LEWO ) then
            zrodlo%polozenia(:,1) = pX1
            do i = 1 , N
                zrodlo%polozenia(i,2) = pY1 + i - 1
            enddo
            zrodlo%r1             = (/pX1,pY1/)*DX
            zrodlo%r2             = (/pX1,pYN/)*DX
            zrodlo%hny            = (pYN + pY1)/2.0*DX

            do i = 1 , N
                zrodlo%SpinUVec(i,:)  = SUTOTAL(zrodlo%polozenia(i,1),zrodlo%polozenia(i,2),:)
            enddo


        else ! dla zrodel gora dol
            zrodlo%polozenia(:,2) = pX1
            do i = 1 , N
                zrodlo%polozenia(i,1) = pY1 + i - 1
            enddo
            zrodlo%r1             = (/pY1,pX1/)*DX
            zrodlo%r2             = (/pYN,pX1/)*DX
            zrodlo%hnx            = (pYN + pY1)/2.0*DX
        endif


         if(lModow == 0)  return

        call zrodlo%spinzrodlo_wypisz_info()

        call zrodlo%spinzrodlo_oblicz_JinJout()

        select case (zrodlo%bKierunek)
        case (ZRODLO_KIERUNEK_PRAWO)
            dir = +1
            zrodlo%dir = +1
        case (ZRODLO_KIERUNEK_LEWO)
            dir = -1
            zrodlo%dir = -1
        case (ZRODLO_KIERUNEK_GORA)
        case (ZRODLO_KIERUNEK_DOL)
        endselect

        ntmp =  zrodlo%liczba_modow + zrodlo%liczba_evans
        select case (zrodlo%bKierunek)
        case (ZRODLO_KIERUNEK_PRAWO,ZRODLO_KIERUNEK_LEWO)
            ! ----------------------------------------------------------------------
            ! Macierze warunkow brzegowych dla zrodel skierowanych w prawo ------>
            ! Aij    = < X(-i) | X(+j) >
            ! Sij^-1 = < X(-i) | X(-j) >
            ! ----------------------------------------------------------------------
            ! ----------------------------------------------------------------------
            ! Macierze warunkow brzegowych dla zrodel skierowanych w lewo  <-------
            ! Aij    = < X(+i) | X(-j) >
            ! Sij^-1 = < X(+i) | X(+j) >
            ! ----------------------------------------------------------------------
            do i = 1 , ntmp
            do j = 1 , lModow
                zrodlo%Aij(i,j) = sum( conjg(zrodlo%ChiMod(:,i,-dir))*zrodlo%ChiMod(:,j,+dir) )*DX
            end do
            end do

            allocate(tempB(ntmp,1))
            do i = 1 , ntmp
            do j = 1 , ntmp
                zrodlo%Sij(i,j) = sum( conjg(zrodlo%ChiMod(:,i,-dir))*zrodlo%ChiMod(:,j,-dir) )*DX
            end do
            end do
            tempB = 1
            call zgaussj(zrodlo%Sij(1:ntmp,1:ntmp),ntmp,ntmp,tempB(1:ntmp,1),1,1)
            deallocate(tempB)


            ! Wyznaczanie macierzy pomocniczej zawierajacej iloczyny macierzy Sij i wektora Chi
            zrodlo%SijChiAuxMat = 0
            do k = 1 , ntmp
            do i = 1 , zrodlo%N

                lambda0 = zrodlo%ChiLambda(k,-dir) * exp(II*BZ*zrodlo%hnY*DX)
                lambda = ( lambda0)**( +1)  - ( lambda0)**( -1)
                do p = 1 , ntmp
                    zrodlo%SijChiAuxMat(k,i,+1) = zrodlo%SijChiAuxMat(k,i,+1) +&
                        lambda*DX*( zrodlo%Sij(k,p)*conjg(zrodlo%ChiModUp  (i,p,-dir)) )

                    zrodlo%SijChiAuxMat(k,i,-1) = zrodlo%SijChiAuxMat(k,i,-1) +&
                        lambda*DX*( zrodlo%Sij(k,p)*conjg(zrodlo%ChiModDown(i,p,-dir)) )
                end do
            end do
            end do

            zrodlo%Sigma = 0
            do v = 1 , zrodlo%N
            do i = 1 , zrodlo%N

                do k = 1 , ntmp
                  zrodlo%Sigma(v,i,+1,+1) = zrodlo%Sigma(v,i,+1,+1) +&
                            zrodlo%ChiModUp(v,k,-dir) * zrodlo%SijChiAuxMat(k,i,+1)

                  zrodlo%Sigma(v,i,+1,-1) = zrodlo%Sigma(v,i,+1,-1) +&
                            zrodlo%ChiModUp(v,k,-dir) * zrodlo%SijChiAuxMat(k,i,-1)

                  zrodlo%Sigma(v,i,-1,-1) = zrodlo%Sigma(v,i,-1,-1) +&
                            zrodlo%ChiModDown(v,k,-dir) * zrodlo%SijChiAuxMat(k,i,-1)

                  zrodlo%Sigma(v,i,-1,+1) = zrodlo%Sigma(v,i,-1,+1) +&
                            zrodlo%ChiModDown(v,k,-dir) * zrodlo%SijChiAuxMat(k,i,+1)


!                lambda0 = zrodlo%ChiLambda(k,-dir) * exp(II*BZ*zrodlo%hnY*DX)
!                lambda = ( lambda0)**( +1)  - ( lambda0)**( -1)
!                do p = 1 , ntmp
!                zrodlo%Sigma(v,i,+1,+1) = zrodlo%Sigma(v,i,+1,+1) +&
!                                  zrodlo%ChiModUp(v,k,-dir) * lambda*DX*( zrodlo%Sij(k,p)* &
!                            conjg(zrodlo%ChiModUp(i,p,-dir)) )
!
!                zrodlo%Sigma(v,i,+1,-1) = zrodlo%Sigma(v,i,+1,-1) +&
!                                  zrodlo%ChiModUp(v,k,-dir) * lambda*DX*( zrodlo%Sij(k,p)* &
!                            conjg(zrodlo%ChiModDown(i,p,-dir)) )
!
!                zrodlo%Sigma(v,i,-1,-1) = zrodlo%Sigma(v,i,-1,-1) +&
!                                  zrodlo%ChiModDown(v,k,-dir) * lambda*DX*( zrodlo%Sij(k,p)* &
!                            conjg(zrodlo%ChiModDown(i,p,-dir)) )
!
!                zrodlo%Sigma(v,i,-1,+1) = zrodlo%Sigma(v,i,-1,+1) +&
!                                  zrodlo%ChiModDown(v,k,-dir) * lambda*DX*( zrodlo%Sij(k,p)* &
!                            conjg(zrodlo%ChiModUp(i,p,-dir)) )
!                enddo

                enddo
            enddo
            enddo
            zrodlo%Sigma = - dir*zrodlo%Sigma

        case (ZRODLO_KIERUNEK_GORA)
            print*,"MODZRODLO: Zrodlo wejsciowe nie obslugiwane!"
            stop

!
!            ! Wzynaczanie wektora pomocniczego
!            zrodlo%deltamk(:) = 0 ! dajemy zero bo nie bedzie uzywane, bo tutaj cechowanie zalezy jeszcze od x
!
!            ! ----------------------------------------------------------------------
!            ! Macierze warunkow brzegowych dla zrodel skierowanych w gore ------>
!            ! Aij    = < X(-i) | X(+j) >
!            ! Sij^-1 = < X(-i) | X(-j) >
!            ! ----------------------------------------------------------------------
!            do i = 1 , ntmp
!            do j = 1 , lModow
!                zrodlo%Aij(i,j) = 0
!                do k = 1 , zrodlo%N
!
!
!                    ctmp = exp( (zrodlo%polozenia(k,2))*DX*(&
!                            &II*zrodlo%k_m_in(j) +&
!                            &II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
!
!
!                    zrodlo%Aij(i,j) = zrodlo%Aij(i,j) + ctmp*( conjg(zrodlo%Chi_m_out(k,i))*zrodlo%Chi_m_in (k,j) )*DX
!                enddo
!
!                !zrodlo%Aij(i,j) = sum( conjg(zrodlo%Chi_m_out(:,i))*zrodlo%Chi_m_in (:,j) )*DX
!            end do
!            end do
!
!            allocate(tempB(ntmp,1))
!            do i = 1 , ntmp
!            do j = 1 , ntmp
!                zrodlo%Sij(i,j) = 0
!                do k = 1 , zrodlo%N
!                    if( j <= lModow) then
!                    ctmp = exp( (zrodlo%polozenia(k,2))*DX*(-II*zrodlo%k_m_in(j) +&
!                            &II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
!
!                    else
!                    ctmp = exp( zrodlo%polozenia(k,2)*DX*(&
!                            &-zrodlo%k_m_in(j) +&
!                            &II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
!                    endif
!
!!                    ctmp = 1
!
!                    zrodlo%Sij(i,j) = zrodlo%Sij(i,j) + ctmp*( conjg(zrodlo%Chi_m_out(k,i))*zrodlo%Chi_m_out(k,j) )*DX
!                enddo
!
!            end do
!            end do
!            tempB = 1
!            call zgaussj(zrodlo%Sij(1:ntmp,1:ntmp),ntmp,ntmp,tempB(1:ntmp,1),1,1)
!            deallocate(tempB)
!
!            ! Wyznaczanie macierzy pomocniczej zawierajacej iloczyny macierzy Sij i wektora Chi
!            do i = 1 , ntmp
!            do j = 1 , N
!                zrodlo%SijChiAuxMat(i,j) = sum( zrodlo%Sij(i,:) * conjg(zrodlo%Chi_m_out(j,:)) )
!            end do
!            end do


            ! ---------------------------------------------------------------- !
            case (ZRODLO_KIERUNEK_DOL)
            print*,"MODZRODLO: Zrodlo wejsciowe nie obslugiwane!"
            stop
!            ! Wzynaczanie wektora pomocniczego
!            zrodlo%deltamk(:) = 0 ! dajemy zero bo nie bedzie uzywane, bo tutaj cechowanie zalezy jeszcze od x
!
!            ! ----------------------------------------------------------------------
!            ! Macierze warunkow brzegowych dla zrodel skierowanych w dol  <-------
!            ! Aij    = < X(+i) | X(-j) >
!            ! Sij^-1 = < X(+i) | X(+j) >
!            ! ----------------------------------------------------------------------
!            do i = 1 , ntmp
!            do j = 1 , lModow
!                zrodlo%Aij(i,j) = 0
!                do k = 1 , zrodlo%N
!!                    ctmp = exp( zrodlo%polozenia(k,2)*DX*(&
!!                            &-II*zrodlo%k_m_in(j) +&
!!                            & II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
!
!                    ctmp = exp( zrodlo%polozenia(k,2)*DX*(-2*II*zrodlo%k_m_in(j) ) )
!
!                    zrodlo%Aij(i,j) = zrodlo%Aij(i,j) + ctmp*( conjg(zrodlo%Chi_m_in(k,i))*zrodlo%Chi_m_out (k,j) )*DX
!                enddo
!
!
!            end do
!            end do
!
!            allocate(tempB(ntmp,1))
!            do i = 1 , ntmp
!            do j = 1 , ntmp
!                zrodlo%Sij(i,j) = 0
!                do k = 1 , zrodlo%N
!!                    if( j <= lModow) then
!!                    ctmp = exp( (zrodlo%polozenia(k,2))*DX*( &
!!                            &+II*zrodlo%k_m_in(j) + &
!!                            & II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
!!                    else
!!                    ctmp = exp( (zrodlo%polozenia(k,2))*DX*( &
!!                            & zrodlo%k_m_out(j) + &
!!                            & II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
!!                    endif
!
!                    if( j <= lModow) then
!                    ctmp = exp( (zrodlo%polozenia(k,2))*DX*(II*zrodlo%k_m_in(j) +&
!                            &II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
!
!                    else
!                    ctmp = exp( zrodlo%polozenia(k,2)*DX*(&
!                            &-zrodlo%k_m_in(j) +&
!                            &II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
!                    endif
!
!
!                    zrodlo%Sij(i,j) = zrodlo%Sij(i,j) + ctmp*( conjg(zrodlo%Chi_m_in(k,i))*zrodlo%Chi_m_in(k,j) )*DX
!
!                enddo
!
!            end do
!            end do
!            tempB = 1
!            call zgaussj(zrodlo%Sij(1:ntmp,1:ntmp),ntmp,ntmp,tempB(1:ntmp,1),1,1)
!            deallocate(tempB)
!
!            ! Wyznaczanie macierzy pomocniczej zawierajacej iloczyny macierzy Sij i wektora Chi
!            do i = 1 , ntmp
!            do j = 1 , N
!                zrodlo%SijChiAuxMat(i,j) = sum( zrodlo%Sij(i,:) * conjg(zrodlo%Chi_m_in(j,:)) )
!            end do
!            end do


        case default
            print*,"Modzrodlo:: Nie ma takiego typu zrodla jak",zrodlo%bKierunek
            stop
        endselect



    end subroutine spinzrodlo_ustaw



    ! --------------------------------------------------------------------
    ! Oblicza wartosci dla wektora wyrazow wolnych jak rozwiazywany jest
    ! uklad rownan na Psi(u,v). Jako argument podajemy jedynie dx w
    ! jednostkach donorowych. Utworzony wektor przechowywany jest w tablicy
    ! Fj(v). Procedura uzupelnia wektor wyrazow wolnych w zaleznosci
    ! od tego jaki jest kierunek zrodla. Przed wywolaniem tej procedury
    ! nalezy ustawic odpowiednie wartosci amplitud wejsciowych ck(i).
    ! --------------------------------------------------------------------
    subroutine spinzrodlo_oblicz_Fj(zrodlo)
        class(cspinzrodlo)  ::  zrodlo
        doubleprecision :: pdx
        integer         :: i , pi, pj, k , p, q , dir
        complex*16      :: post , Xkn , deltapk , deltamk , kvec , lambda , lambda0
        doubleprecision :: dx, Ef, Bz , ypos , xpos , bpart

        dx  = atomic_DX*L2LR
        Ef  = atomic_Ef/1000.0/Rd
        BZ  = BtoDonorB(atomic_Bz)


        select case (zrodlo%bKierunek)
        case (ZRODLO_KIERUNEK_PRAWO)
            dir = +1
        case (ZRODLO_KIERUNEK_LEWO)
            dir = -1
        case (ZRODLO_KIERUNEK_GORA)
        case (ZRODLO_KIERUNEK_DOL)
        endselect

        ! obliczanie wektora pomocniczego przy liczeniu wyrazu wolnego
        ! nie zalezu on od kierunku zrodla
        do k = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
            zrodlo%SijAijCkAuxVec(k) = 0
            do p = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
            do q = 1 , zrodlo%liczba_modow
                zrodlo%SijAijCkAuxVec(k) = zrodlo%SijAijCkAuxVec(k) + zrodlo%Sij(k,p)*zrodlo%Aij(p,q)*zrodlo%ck(q)
            enddo
            enddo
        enddo


        do i = 1 , zrodlo%N
                zrodlo%Fj(i,:) = 0
                pi         =  zrodlo%polozenia(i,1)
                pj         =  zrodlo%polozenia(i,2)
                xpos       =  (pi) * dx
                ypos       =  (pj) * dx

                select case (zrodlo%bKierunek)
                ! ---------------------------------------------------------------- !
                ! Wypelniamy wektor wyrazow wolnych dla zrodej skierowanych w prawo
                !
                !                        --------------->
                !
                case (ZRODLO_KIERUNEK_PRAWO,ZRODLO_KIERUNEK_LEWO)

                    do k = 1 , zrodlo%liczba_modow
                        lambda0         = zrodlo%ChiLambda(k,+dir) * exp(II*BZ*zrodlo%hnY*DX)
                        lambda          = ( lambda0 )**( +1) - ( lambda0 )**( -1)

                        zrodlo%Fj(i,+1) = zrodlo%Fj(i,+1) + zrodlo%ck(k)*lambda*zrodlo%ChiModUp(i,k,+dir)
                        zrodlo%Fj(i,-1) = zrodlo%Fj(i,-1) + zrodlo%ck(k)*lambda*zrodlo%ChiModDown(i,k,+dir)
                    enddo

                    do k = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
                        Xkn             = zrodlo%SijAijCkAuxVec(k)
                        lambda0         = zrodlo%ChiLambda(k,-dir) * exp(II*BZ*zrodlo%hnY*DX)
                        lambda          = ( lambda0 )**( +1) - ( lambda0 )**( -1)

                        zrodlo%Fj(i,+1) = zrodlo%Fj(i,+1) - Xkn*lambda*zrodlo%ChiModUp(i,k,-dir)
                        zrodlo%Fj(i,-1) = zrodlo%Fj(i,-1) - Xkn*lambda*zrodlo%ChiModDown(i,k,-dir)
                    enddo

                    zrodlo%Fj(i,:) = -dir*zrodlo%Fj(i,:)

                case (ZRODLO_KIERUNEK_GORA)
                    print*,"SPINMODZRODLO: Zrodlo wejsciowe nie obslugiwane!"
                    stop
!                    post             = -(0.5/DX/DX)
!                    bpart            = Bz*( xpos  - zrodlo%hnx )
!                    do k = 1 , zrodlo%liczba_modow
!
!                        kvec         = zrodlo%k_m_in(k)
!                        deltamk      = exp( (ypos +  DX)*(+II*kvec + II*bpart) ) &
!                                    &- exp( (ypos -  DX)*(+II*kvec + II*bpart) )
!                        zrodlo%Fj(i) = zrodlo%Fj(i) + zrodlo%ck(k)*deltamk*zrodlo%Chi_m_in (i,k)
!
!                    enddo
!                    do k = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
!
!                        if( k <= zrodlo%liczba_modow) then
!                            deltamk = exp( (ypos +  DX)*(-II*zrodlo%k_m_in(k) + II*bpart) ) - exp(  (ypos-DX)*(-II*zrodlo%k_m_in(k) + II*bpart) )
!                        else
!                            deltamk = exp( (ypos + DX)*(-zrodlo%k_m_in(k) + II*bpart)  ) - exp( (ypos-DX)*( -zrodlo%k_m_in(k) + II*bpart) )
!                        endif
!
!                        Xkn          = zrodlo%SijAijCkAuxVec(k)
!                        zrodlo%Fj(i) = zrodlo%Fj(i) - Xkn*deltamk*zrodlo%Chi_m_out(i,k)
!                    enddo
                case (ZRODLO_KIERUNEK_DOL)
                    print*,"SPINMODZRODLO: Zrodlo wejsciowe nie obslugiwane!"
                    stop
                endselect
!                zrodlo%Fj(i) = zrodlo%Fj(i)*post

        enddo ! end of do i=1, zrodlo%N

    end subroutine spinzrodlo_oblicz_Fj

   subroutine spinzrodlo_oblicz_dk(zrodlo,VPHI,GINDEX,nx,ny)
        class(cspinzrodlo)          :: zrodlo
        complex*16,dimension(:) :: VPHI
        integer,dimension(nx,ny,-1:1):: GINDEX
        integer :: nx,ny
        doubleprecision :: dx

        integer :: p,q,k,i,pi,pj,dir
        complex*16 :: dk , alpha_p, beta_p


        dx  = atomic_DX * L2LR
        dir = zrodlo%dir

        do k = 1 , zrodlo%liczba_modow
            dk = 0
            do p = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
                alpha_p = 0
                do i = 2 , zrodlo%N-1
                    pi   = zrodlo%polozenia(i,1)
                    pj   = zrodlo%polozenia(i,2)

                    alpha_p = alpha_p + dx*conjg(zrodlo%ChiModUp  (i,p,-dir))*VPHI(GINDEX(pi,pj,+1))
                    alpha_p = alpha_p + dx*conjg(zrodlo%ChiModDown(i,p,-dir))*VPHI(GINDEX(pi,pj,-1))


                enddo
                beta_p = 0
                do q = 1 , zrodlo%liczba_modow
                    beta_p = beta_p + zrodlo%Aij(p,q)*zrodlo%ck(q)
                enddo
                dk = dk + zrodlo%Sij(k,p)*( alpha_p - beta_p )
            enddo ! end do p
            zrodlo%dk(k) = dk

        enddo ! end do k

        call zrodlo%spinzrodlo_wypisz_ckdk()

    end subroutine spinzrodlo_oblicz_dk
!
    ! wzory brane z mojej pracy inzynierskiej (strona 12 - tabela)
    subroutine spinzrodlo_oblicz_JinJout(zrodlo)
        class(cspinzrodlo)    :: zrodlo
        double precision  :: dx
        integer          :: k , pj , i
        double precision :: Jin , Jout , Bz ,kvecin , kvecout , rhoup , rhodwn, t0 , tc0 , kvec , hny
        complex*16 ::  Xup,Xdwn,tlj(2,2,2)
        BZ  = BtoDonorB(atomic_Bz)
        DX  = atomic_DX * L2LR
        t0  = 0.5/DX/DX
        tc0 = 0.5/DX

        do k = 1 , zrodlo%liczba_modow
            Jin  = 0
            Jout = 0
            kvecin = imag(zrodlo%ChiKvec(k,+1))
            kvecout= imag(zrodlo%ChiKvec(k,-1))

            select case (zrodlo%bKierunek)
            ! ------------------------------------------------------------------
            !
            !                          ---------------->
            !                          <----------------
            !
            ! ------------------------------------------------------------------
            case (ZRODLO_KIERUNEK_PRAWO,ZRODLO_KIERUNEK_LEWO)
                hny    = zrodlo%hnY/dx

                do i = 1 , zrodlo%N
                    pj    = zrodlo%polozenia(i,2)

!                    Xup    = zrodlo%ChiModUp  (i,k,1)
!                    Xdwn   = zrodlo%ChiModDown(i,k,1)
!
!                    rhoup  = abs(Xup)**2
!                    rhodwn = abs(Xdwn)**2
!                    Jin    = Jin -  rhoup * sin(dx*dx*(pj-zrodlo%hnY/dx)*Bz-kvecin*dx)
!                    Jin    = Jin -  rhodwn* sin(dx*dx*(pj-zrodlo%hnY/dx)*Bz-kvecin*dx)
!                    Jin    = Jin + dx*( (rhoup-rhodwn)*so_loc*zrodlo%Ey(i) + II*so_rashba*( Xdwn*conjg(Xup) - conjg(Xdwn)*Xup ) )
!
!
!                    Xup    = zrodlo%ChiModUp  (i,k,-1)
!                    Xdwn   = zrodlo%ChiModDown(i,k,-1)
!
!                    rhoup  = abs(Xup)**2
!                    rhodwn = abs(Xdwn)**2
!                    Jout    = Jout -  rhoup * sin(dx*dx*(pj-zrodlo%hnY/dx)*Bz-kvecout*dx)
!                    Jout    = Jout -  rhodwn* sin(dx*dx*(pj-zrodlo%hnY/dx)*Bz-kvecout*dx)
!                    Jout    = Jout + dx*( (rhoup-rhodwn)*so_loc*zrodlo%Ey(i) + II*so_rashba*( Xdwn*conjg(Xup) - conjg(Xdwn)*Xup ) )


        Xup    = zrodlo%ChiModUp  (i,k,1)
        Xdwn   = zrodlo%ChiModDown(i,k,1)
        kvec   = imag(zrodlo%ChiKvec(k,+1))

        tlj(1,1,1) = - t0 *  exp(-II*(DX*DX*BZ*(pj-hNY))) * exp(+II*DX*zrodlo%Ey(i)*so_loc) ! l+1
        tlj(1,1,2) = - t0 *  exp(+II*(DX*DX*BZ*(pj-hNY))) * exp(-II*DX*zrodlo%Ey(i)*so_loc) ! l-1

        ! spin down-down
        tlj(2,2,1) = - t0 *  exp(-II*(DX*DX*BZ*(pj-hNY))) * exp(-II*DX*zrodlo%Ey(i)*so_loc) ! l+1
        tlj(2,2,2) = - t0 *  exp(+II*(DX*DX*BZ*(pj-hNY))) * exp(+II*DX*zrodlo%Ey(i)*so_loc) ! l-1

        ! spin up-down
        tlj(1,2,1) = +so_rashba*tc0 * exp(-II*(DX*DX*BZ*(pj-hNY)))
        tlj(1,2,2) = -so_rashba*tc0 * exp(+II*(DX*DX*BZ*(pj-hNY)))! l-1

        ! spin up-down
        tlj(2,1,1) = -so_rashba*tc0 * exp(-II*(DX*DX*BZ*(pj-hNY)))
        tlj(2,1,2) = +so_rashba*tc0 * exp(+II*(DX*DX*BZ*(pj-hNY)))! l-1
        ! wezly do przodu
        Jin = Jin + dx*II*( tlj(1,1,1)*conjg(Xup)*YY(kvec,Xup,1)     - conjg(tlj(1,1,1)*conjg(Xup)*YY(kvec,Xup,1)) )
        Jin = Jin + dx*II*( tlj(2,2,1)*conjg(Xdwn)*YY(kvec,Xdwn,1) - conjg(tlj(2,2,1)*conjg(Xdwn)*YY(kvec,Xdwn,1)) )

        Jin = Jin + dx*II*( tlj(1,2,1)*conjg(Xup)*YY(kvec,Xdwn,1)   - conjg(tlj(1,2,1)*conjg(Xup)*YY(kvec,Xdwn,1)) )
        Jin = Jin + dx*II*( tlj(2,1,1)*conjg(Xdwn)*YY(kvec,Xup,1)   - conjg(tlj(2,1,1)*conjg(Xdwn)*YY(kvec,Xup,1)) )


        Xup    = zrodlo%ChiModUp  (i,k,-1)
        Xdwn   = zrodlo%ChiModDown(i,k,-1)
        kvec   = imag(zrodlo%ChiKvec(k,-1))

        ! wezly do przodu
        Jout = Jout + dx*II*( tlj(1,1,1)*conjg(Xup)*YY(kvec,Xup,1)     - conjg(tlj(1,1,1)*conjg(Xup)*YY(kvec,Xup,1)) )
        Jout = Jout + dx*II*( tlj(2,2,1)*conjg(Xdwn)*YY(kvec,Xdwn,1) - conjg(tlj(2,2,1)*conjg(Xdwn)*YY(kvec,Xdwn,1)) )

        Jout = Jout + dx*II*( tlj(1,2,1)*conjg(Xup)*YY(kvec,Xdwn,1)   - conjg(tlj(1,2,1)*conjg(Xup)*YY(kvec,Xdwn,1)) )
        Jout = Jout + dx*II*( tlj(2,1,1)*conjg(Xdwn)*YY(kvec,Xup,1)   - conjg(tlj(2,1,1)*conjg(Xdwn)*YY(kvec,Xup,1)) )



                enddo ! end of i
            case (ZRODLO_KIERUNEK_GORA,ZRODLO_KIERUNEK_DOL)

                print*,"Nie obslugiwane liczenie pradow w kierunku gora dol"
                stop
!                do i = 1 , zrodlo%N
!                    pj    = zrodlo%polozenia(i,1)
!                    Jin   = Jin  + abs(zrodlo%Chi_m_in (i,k))**2*sin( +dx*dx*pj*Bz + dx*(+kvec - zrodlo%hnX*Bz) )
!                    Jout  = Jout + abs(zrodlo%Chi_m_out(i,k))**2*sin( +dx*dx*pj*Bz + dx*(-kvec - zrodlo%hnX*Bz) )
!                enddo ! end of i
            endselect

        zrodlo%ChiCurr(k,+1) = Jin
        zrodlo%ChiCurr(k,-1) = Jout

        enddo ! end of k

!        zrodlo%ChiCurr(:,-1) = zrodlo%ChiCurr(:,-1) / zrodlo%ChiCurr(1,+1)
!        zrodlo%ChiCurr(:,+1) = zrodlo%ChiCurr(:,+1) / zrodlo%ChiCurr(1,+1)

    end subroutine spinzrodlo_oblicz_JinJout


complex*16 function YY(kvec,chi,i) result(rval)
    double precision :: kvec
    complex*16 :: chi
    integer :: i
    doubleprecision :: dx
    dx = atomic_DX * L2LR

    rval = exp(II*kvec*(i)*dx)*chi

end function YY

!! ---------------------------------------------------------------------------------------------
!!
!!   Funkcje klasy abszrodla
!!
!! ---------------------------------------------------------------------------------------------
!    subroutine abs_zrodlo_zwolnij_pamiec(zrodlo)
!        class(cabs_zrodlo) :: zrodlo
!
!        if(TRANS_DEBUG==.true.) print*,"Abszrodlo: Zwalanianie pamieci"
!        if(allocated(zrodlo%polozenia)) deallocate(zrodlo%polozenia);
!
!
!    end subroutine abs_zrodlo_zwolnij_pamiec
!
!
!! --------------------------------------------------------------------------------- !
!! Podobnie jak to jest dla normalnej klasy zrodlo. Potrzebna jest tylko informacja
!! o energii fermiego
!! --------------------------------------------------------------------------------- !
!
!    subroutine abs_zrodlo_ustaw(zrodlo,pY1,pYN,pX1,pEf,pKierunek)
!        class(cabs_zrodlo)             ::  zrodlo
!        integer,intent(in)         ::  pY1,pYN,pX1
!        doubleprecision,intent(in) ::  pEf
!        integer,intent(in)         ::  pKierunek ! enum
!
!        double precision :: Ef
!        integer :: N,i
!
!        print*,"Dodawanie zrodla z transparentnymi w.b"
!
!        Ef  = pEf/1000.0/Rd
!
!        N = pYN - pY1 + 1; ! to jest tak samo liczone ale zmienne pYN , pY1 maja inna interpretacje
!        zrodlo%kvec = sqrt(2*Ef)
!        zrodlo%N    = N
!
!        if(allocated(zrodlo%polozenia)) deallocate(zrodlo%polozenia)
!
!        allocate(zrodlo%polozenia(N,2))
!        zrodlo%polozenia = 0
!        ! ustawianie parametrow zrodel
!        zrodlo%bKierunek = pKierunek
!
!        if( pKierunek == ZRODLO_KIERUNEK_PRAWO .or. pKierunek == ZRODLO_KIERUNEK_LEWO ) then
!            zrodlo%polozenia(:,1) = pX1
!            do i = 1 , N
!                zrodlo%polozenia(i,2) = pY1 + i - 1
!            enddo
!        else ! dla zrodel gora dol
!            zrodlo%polozenia(:,2) = pX1
!            do i = 1 , N
!                zrodlo%polozenia(i,1) = pY1 + i - 1
!            enddo
!        endif
!
!        print*,"kvec      = " , zrodlo%kvec * L2LR , "[nm]"
!        print*,"kierunek  ="  , zrodlo%bKierunek
!        print*,"N         ="  , zrodlo%N
!
!
!    end subroutine abs_zrodlo_ustaw
!
!! --------------------------------------------------------------------------------- !
!! Kopuje dane z jednego zrodla do drugiego
!! --------------------------------------------------------------------------------- !
!
!    subroutine abs_zrodlo_skopiuj(zrodlo,od_zrodla)
!        class(cabs_zrodlo)       ::  zrodlo,od_zrodla
!
!        zrodlo%bKierunek = od_zrodla%bKierunek
!        zrodlo%kvec      = od_zrodla%kvec
!        zrodlo%N         = od_zrodla%N
!        if(allocated(zrodlo%polozenia)) deallocate(zrodlo%polozenia)
!        allocate(zrodlo%polozenia(zrodlo%N,2))
!
!        zrodlo%polozenia = od_zrodla%polozenia
!    end subroutine abs_zrodlo_skopiuj

endmodule modspinzrodlo
