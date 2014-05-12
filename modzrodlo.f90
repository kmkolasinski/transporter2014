module modzrodlo
    use modutils
    use modjed
    use modpop
    use ifport
    use modinip
implicit none


    ! =========================================================
    !                   Klasa zrodla
    ! =========================================================
    ! przechowuje informacje o zrodle
    type czrodlo
        integer :: N ! liczba oczek siatki
        logical :: bWejscie       ! jesli wejscie to zrodlo bedzie traktowane jako zrodlo wejsciowe
        logical :: bKierunek      ! ustala zwrot zrodla, 1 - w kierunku dodatnich wartosci x lub y, 0 - przeciwnie
        double precision :: hnY   ! wysokosc zrodla w oczkach siatki, moze byc polowkowa
        double precision,dimension(2) :: r1  ! nie uzywane w programie polozenie pierwszego punktu w nm
        double precision,dimension(2) :: r2  ! nie uzywane w programie polozenie ostatniego punktu w nm
        integer,dimension(:,:),allocatable :: polozenia; ! wekor polozen (n,1) ,(n,2) na siatce numerycznej mapowanie na indeks
                                                         ! bedzie odbywac sie przez tablice GINDEX(i,j)
        complex*16,dimension(:,:),allocatable     :: Chi_m_in ! macierz modow wchodzacych, (N,M) - M liczba modow
        complex*16,dimension(:,:),allocatable     :: Chi_m_out! mody wychodzace
        double precision,dimension(:),allocatable :: k_m_in   ! wektory falowe zwiazane z modami wchodz.
        double precision,dimension(:),allocatable :: k_m_out  ! z wychodzacymi
        integer                                   :: liczba_modow ! liczba dostepnym w zrodle modow
        integer                                   :: liczba_evans ! liczba modow evanescentnych
        complex*16,dimension(:),allocatable       :: ck       ! amplitudy wchodzace
        complex*16,dimension(:),allocatable       :: dk       ! amplitudy odbite
        integer                                   :: rozbieg  ! zakres rzutowania w procesie wyznaczania ck, dk
        complex*16,dimension(200,200)             :: m_r,m_t  ! macierz wspolczynnikow odbicia i transmisji
        complex*16,dimension(:,:),allocatable     :: Smat
        contains
        ! -------------------------------------------------
        ! Procedury zrodla
        ! -------------------------------------------------
        procedure, public, pass(modn) :: zrodlo_zwolnij_pamiec!() wiadomo
        procedure, public, pass(modn) :: zrodlo_alokuj_pamiec !(pN,lM) podajemy liczbe punktow oraz liczbe modow
        procedure, public, pass(modn) :: zrodlo_wypisz_info   !() wypisuje parametry zrodla
        procedure, public, pass(modn) :: zrodlo_wypisz_ckdk   !() wypisyje amplitudy ck, dk
        procedure, public, pass(modn) :: zrodlo_zapisz_m_rt   !(mod) zapisuje biezace amplitudy od macierzy m_r i m_t dla zadanego modu
        procedure, public, pass(modn) :: zrodlo_ustaw   !(mod) zapisuje biezace amplitudy od macierzy m_r i m_t dla zadanego modu
        procedure, public, pass(modn) :: zrodlo_zapisz_mody   !(mod,filename,dx) zapisuje  mody danego zrodla do pliku, podajemy dx [nm]

    end type czrodlo

contains

!------------------------------------------------------------
!                 Funkcje klasy zrodla
!------------------------------------------------------------
    subroutine zrodlo_zwolnij_pamiec(modn)
        class(czrodlo) :: modn

        if(TRANS_DEBUG==1) print*,"Zrodlo: Zwalanianie pamieci"
        if(allocated(modn%polozenia)) deallocate(modn%polozenia);
        if(allocated(modn%Chi_m_in))  deallocate(modn%Chi_m_in);
        if(allocated(modn%Chi_m_out)) deallocate(modn%Chi_m_out);
        if(allocated(modn%k_m_in))    deallocate(modn%k_m_in);
        if(allocated(modn%k_m_out))   deallocate(modn%k_m_out);
        if(allocated(modn%ck))        deallocate(modn%ck);
        if(allocated(modn%dk))        deallocate(modn%dk);
        if(allocated(modn%Smat))      deallocate(modn%Smat);
    end subroutine zrodlo_zwolnij_pamiec

! ----------------------------------------------------------
! Alokuje zrodlo dla zadanej liczby oczek siatki - pN
! oraz lM - liczby modow oraz liczby modow evanescentnych - lEvanMods
! ----------------------------------------------------------
    subroutine zrodlo_alokuj_pamiec(modn,pN,lM,lEvanMods)
        class(czrodlo) :: modn
        integer,intent(in) :: pN,lM,lEvanMods

        modn%N            = pN;
        modn%liczba_modow = lM;
        modn%liczba_evans = lEvanMods;

        if(TRANS_DEBUG==1) print*,"Zrodlo: Alokowanie pamieci"
        allocate(modn%polozenia(pN,2) );
        allocate(modn%Chi_m_in (pN,lM+lEvanMods));
        allocate(modn%Chi_m_out(pN,lM+lEvanMods));
        allocate(modn%k_m_in   (lM+lEvanMods));
        allocate(modn%k_m_out  (lM+lEvanMods));
        allocate(modn%ck       (lM+lEvanMods));
        allocate(modn%dk       (lM+lEvanMods));
        allocate(modn%Smat(pN,pN))
        modn%ck = 0
        modn%dk = 0

    end subroutine zrodlo_alokuj_pamiec

    subroutine zrodlo_wypisz_info(modn)
        class(czrodlo) :: modn
        integer :: i
        if(TRANS_DEBUG==.true.) then
        print*,"Zrodlo: Wypisywanie informacji o zrodle:"
        print*,"    N       = ",modn%N
        print*,"    L. modow= ",modn%liczba_modow
        print*,"    Rozbieg = ",modn%rozbieg
        print*,"    Wejscie = ",modn%bWejscie
        print*,"    Kierunek= ",modn%bKierunek
        print*,"    HNY     = ",modn%hnY
        write(*,"(A,2f10.4,A)"),"    R1(x,y) = (",modn%r1,")"
        write(*,"(A,2f10.4,A)"),"    R2(x,y) = (",modn%r2,")"
        print*,"    Mody wejsciowe:"
        do i = 1 , modn%liczba_modow
            print*,"        Kin [",i,"]=",modn%k_m_in(i)*L2LR
        enddo
        print*,"    Mody wyjsciowe:"
        do i = 1 , modn%liczba_modow
            print*,"        Kout[",i,"]=",modn%k_m_out(i)*L2LR
        enddo
        endif
    endsubroutine zrodlo_wypisz_info

    subroutine zrodlo_wypisz_ckdk(modn)
        class(czrodlo) :: modn
        integer :: i
        if(TRANS_DEBUG==1) then
        do i = 1 , modn%liczba_modow
            write(*,"(A,i4,A,f12.4,A,i4,A,f12.4)"),"Ck(",i,")=",abs(modn%ck(i))**2,"  Dk(",i,")=",abs(modn%dk(i))**2
        enddo
        endif
    endsubroutine zrodlo_wypisz_ckdk

    subroutine zrodlo_zapisz_m_rt(modn,mod)
        class(czrodlo) :: modn
        integer,intent(in) :: mod ! numer modu wchodzacego
        integer            :: i
        do i = 1 , modn%liczba_modow
            modn%m_r(mod,i) = modn%dk(i)
            modn%m_t(mod,i) = modn%ck(i)
        enddo
    endsubroutine zrodlo_zapisz_m_rt

    ! -------------------------------------------------------------------------
    ! Zapisuje mody do pliku o nazwie filename, podajemy rozwniez krok siatki w nm
    ! -------------------------------------------------------------------------
    subroutine zrodlo_zapisz_mody(modn,filename,dx)
        class(czrodlo) :: modn
        character(*) :: filename
        double precision :: dx
        integer :: i
        print*,"Zrodlo: Zapisywanie modow do pliku:",filename
        open(4193,file=filename)
        do i = 1 , modn%N
            write(4193,"(300e20.8)"),modn%polozenia(i,2)*dx,abs(modn%Chi_m_in(i,:))**2,abs(modn%Chi_m_out(i,:))**2
        enddo
        close(4193)
    end subroutine zrodlo_zapisz_mody


    ! --------------------------------------------------------------------
    ! Ustawiamy wczesniej zaalokowane zrodlo przez polecenie systemowe.
    ! Parametry:
    ! pY1,pYN - polozenia y na siatce numerycznej liczone od 1 do NY
    ! pX1 - polozenie X zrodla
    ! pDX,pEf,pBz - DX [nm] , Ef [meV] , BZ [T]
    ! pKierunek - true = prawo , false = lewo : okresla skierowanie zrodla
    ! pWejscie  - czy zrodlo jest wejsciowe, czy wyjsciowe (true,false)
    ! pRozbieg  - liczba oczek siatki brana do liczenia amplitud rozpraszania
    ! pUTOTAL   - referencja do potencjalu ukladu
    ! --------------------------------------------------------------------
    subroutine zrodlo_ustaw(modn,pY1,pYN,pX1,pDX,pEf,pBz,pKierunek,pWejscie,pRozbieg,pUTOTAL)
       class(czrodlo)             ::  modn
       integer,intent(in)         ::  pY1,pYN,pX1
       doubleprecision,intent(in) ::  pDx,pEf,pBz
       logical,intent(in)         ::  pKierunek,pWejscie
       integer,intent(in)         ::  pRozbieg
       double precision,dimension(:,:) :: pUTOTAL ! calkowity potencjal w [meV]
        double precision :: pUvec(pYN - pY1 + 1)
        double precision :: dx, Ef, Bz
        integer :: N ! liczba oczek dla zrodla
        integer :: lModow
        ! zmienne pomocnicze
        integer :: i

        dx  = pdx*L2LR    ! konwertujemy do jednostek donorowych
        Ef  = pEf/1000.0/Rd
        BZ  = BtoDonorB(pBz)

        N = pYN - pY1 + 1;

        do i = 1 , N
            pUVEC(i) = pUTOTAL(pX1,pY1 + i - 1) !*1000*Rd ! wracamy do jed atomowych bo modpop takie potrzebuje
        enddo
        call modpop_inicjalizacja(pDX,N,pEf,pBz,pUVEC)
        call modpop_liczba_podpasm(lModow)
        !call modpop_relacja_dyspersji(6,"rel.txt")
        call modn%zrodlo_alokuj_pamiec(N,lModow,0)

        modn%polozenia(:,1) = pX1
        do i = 1 , N
            modn%polozenia(i,2) = pY1 + i - 1
        enddo
        ! POBIERAMY MODY POPRZECZNE ORAZ ICH WEKTORY FALOWE
        call modpop_get_km(lModow,modn%k_m_in,modn%k_m_out)
        call modpop_get_chi(lModow,N,modn%Chi_m_in,modn%Chi_m_out)

        ! ustawianie parametrow zrodel
        modn%bWejscie       = pWejscie
        modn%bKierunek      = pKierunek
        modn%r1             = (/pX1,pY1/)*DX
        modn%r2             = (/pX1,pYN/)*DX
        modn%hny            = (pYN + pY1)/2.0
        modn%rozbieg        = pRozbieg
        modn%m_r  = 0
        modn%m_t  = 0
        modn%Smat = 0
        call modn%zrodlo_wypisz_info()

        ! MODUL RELACJI DYSPERSJI JUZ NIE POTRZEBNY
        call modpop_zwalnienie_pamieci()
    end subroutine zrodlo_ustaw

endmodule modzrodlo
