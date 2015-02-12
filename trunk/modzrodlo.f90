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
        integer :: bKierunek      ! ustala zwrot zrodla, 1 - w kierunku dodatnich wartosci x lub y, 0 - przeciwnie
        double precision :: hnY , hnX   ! wysokosc zrodla w oczkach siatki, moze byc polowkowa
        double precision,dimension(2) :: r1  ! nie uzywane w programie polozenie pierwszego punktu w nm
        double precision,dimension(2) :: r2  ! nie uzywane w programie polozenie ostatniego punktu w nm
        integer,dimension(:,:),allocatable :: polozenia; ! wekor polozen (n,1) ,(n,2) na siatce numerycznej mapowanie na indeks
                                                         ! bedzie odbywac sie przez tablice GINDEX(i,j)
        complex*16,dimension(:,:),allocatable     :: Chi_m_in ! macierz modow wchodzacych, (N,M) - M liczba modow
        complex*16,dimension(:,:),allocatable     :: Chi_m_out! mody wychodzace
        complex*16,dimension(:),allocatable       :: k_m_in   ! wektory falowe zwiazane z modami wchodz.
        complex*16,dimension(:),allocatable       :: k_m_out  ! z wychodzacymi
        complex*16,dimension(:),allocatable       :: deltamk  ! tablica pomocnicza przyspieszajaca liczenie elementow macierzowych
        complex*16,dimension(:),allocatable       :: deltaink  ! tablica pomocnicza przyspieszajaca liczenie elementow macierzowych
        complex*16,dimension(:),allocatable       :: deltaoutk  ! tablica pomocnicza przyspieszajaca liczenie elementow macierzowych
        integer                                   :: liczba_modow ! liczba dostepnym w zrodle modow
        integer                                   :: liczba_evans ! liczba modow evanescentnych
        complex*16,dimension(:),allocatable       :: ck       ! amplitudy wchodzace
        complex*16,dimension(:),allocatable       :: dk       ! amplitudy odbite
        doubleprecision,dimension(:),allocatable  :: Jin,Jout ! prady wejsciowe zwiazane z ck i wyjsciowe - dk

        complex*16,dimension(200,200)             :: m_r,m_t  ! macierz wspolczynnikow odbicia i transmisji
        complex*16,dimension(:,:),allocatable     :: Aij
        complex*16,dimension(:,:),allocatable     :: Sij
        complex*16,dimension(:,:),allocatable     :: SijChiAuxMat ! macierz pomocnicza zavierajaca iloczyny macierzy Sij i Wektorow Chi
        complex*16,dimension(:),allocatable       :: SijAijCkAuxVec ! pomocniczy wektor z obliczonymi iloczynami
        complex*16,dimension(:),allocatable       :: Fj ! wektor wyrazow wolnych
        logical :: bZaalokowane = .false. ! flaga ktora mowi o tym czy zrodlo ma juz zaalokowana pamiec
        contains
        ! -------------------------------------------------
        ! Procedury zrodla
        ! -------------------------------------------------
        procedure, public, pass(zrodlo) :: zrodlo_zwolnij_pamiec!() wiadomo
        procedure, public, pass(zrodlo) :: zrodlo_alokuj_pamiec !(pN,lM) podajemy liczbe punktow oraz liczbe modow
        procedure, public, pass(zrodlo) :: zrodlo_wypisz_info   !() wypisuje parametry zrodla
        procedure, public, pass(zrodlo) :: zrodlo_wypisz_ckdk   !() wypisyje amplitudy ck, dk
        procedure, public, pass(zrodlo) :: zrodlo_wypisz_JinJout!() wypisyje strumien wejsciowy oraz wyjsciowy
        procedure, public, pass(zrodlo) :: zrodlo_zapisz_m_rt   !(mod) zapisuje biezace amplitudy od macierzy m_r i m_t dla zadanego modu
        procedure, public, pass(zrodlo) :: zrodlo_ustaw   !(mod) zapisuje biezace amplitudy od macierzy m_r i m_t dla zadanego modu
        procedure, public, pass(zrodlo) :: zrodlo_zapisz_mody   !(mod,filename,dx) zapisuje  mody danego zrodla do pliku, podajemy dx [nm]
        procedure, public, pass(zrodlo) :: zrodlo_alfa_v_i   !
        procedure, public, pass(zrodlo) :: zrodlo_oblicz_Fj   !
        procedure, public, pass(zrodlo) :: zrodlo_oblicz_dk!(zrodlo,VPHI,GINDEX,dx)   !
        procedure, public, pass(zrodlo) :: zrodlo_oblicz_JinJout!(zrodlo,dx)

    end type czrodlo

    ! =========================================================
    !                   Klasa abszrodlo
    ! wykorzystywana jest to symulacji transparentnych w.b zaproponowanych
    ! przez nowaka
    ! =========================================================
    type cabs_zrodlo
        integer :: N ! liczba oczek siatki
        integer :: bKierunek      ! ustala zwrot zrodla, 1 - w kierunku dodatnich wartosci x lub y, 0 - przeciwnie
        integer,dimension(:,:),allocatable :: polozenia; ! wekor polozen (n,1) ,(n,2) na siatce numerycznej mapowanie na indeks
        doubleprecision :: kvec ! wirtualny wektor falowy
        contains
        procedure, public, pass(zrodlo) :: abs_zrodlo_zwolnij_pamiec
        procedure, public, pass(zrodlo) :: abs_zrodlo_ustaw
        procedure, public, pass(zrodlo) :: abs_zrodlo_skopiuj!(zrodlo,od_zrodla)
    end type cabs_zrodlo


contains

!------------------------------------------------------------
!                 Funkcje klasy zrodla
!------------------------------------------------------------
    subroutine zrodlo_zwolnij_pamiec(zrodlo)
        class(czrodlo) :: zrodlo

        if(TRANS_DEBUG==.true.) print*,"Zrodlo: Zwalanianie pamieci"
        if(allocated(zrodlo%polozenia)) deallocate(zrodlo%polozenia);
        if(allocated(zrodlo%Chi_m_in))  deallocate(zrodlo%Chi_m_in);
        if(allocated(zrodlo%Chi_m_out)) deallocate(zrodlo%Chi_m_out);
        if(allocated(zrodlo%k_m_in))    deallocate(zrodlo%k_m_in);
        if(allocated(zrodlo%k_m_out))   deallocate(zrodlo%k_m_out);
        if(allocated(zrodlo%deltamk))   deallocate(zrodlo%deltamk);
        if(allocated(zrodlo%deltaink))   deallocate(zrodlo%deltaink);
        if(allocated(zrodlo%deltaoutk))   deallocate(zrodlo%deltaoutk);
        if(allocated(zrodlo%ck))        deallocate(zrodlo%ck);
        if(allocated(zrodlo%dk))        deallocate(zrodlo%dk);
        if(allocated(zrodlo%Jin))       deallocate(zrodlo%Jin);
        if(allocated(zrodlo%Jout))      deallocate(zrodlo%Jout);
        if(allocated(zrodlo%Aij))       deallocate(zrodlo%Aij);
        if(allocated(zrodlo%Sij))       deallocate(zrodlo%Sij);
        if(allocated(zrodlo%SijChiAuxMat))       deallocate(zrodlo%SijChiAuxMat);
        if(allocated(zrodlo%Fj))        deallocate(zrodlo%Fj);
        if(allocated(zrodlo%SijAijCkAuxVec))     deallocate(zrodlo%SijAijCkAuxVec);

        zrodlo%bZaalokowane = .false.

    end subroutine zrodlo_zwolnij_pamiec

! ----------------------------------------------------------
! Alokuje zrodlo dla zadanej liczby oczek siatki - pN
! oraz lM - liczby modow oraz liczby modow evanescentnych - lEvanMods
! ----------------------------------------------------------
    subroutine zrodlo_alokuj_pamiec(zrodlo,pN,lM,lEvanMods)
        class(czrodlo) :: zrodlo
        integer,intent(in) :: pN,lM,lEvanMods

        zrodlo%N            = pN;
        zrodlo%liczba_modow = lM;
        zrodlo%liczba_evans = lEvanMods;

        if(TRANS_DEBUG==.true.) print*,"Zrodlo: Alokowanie pamieci"

        if(zrodlo%bZaalokowane == .true.) then
            call zrodlo%zrodlo_zwolnij_pamiec()
        endif

        allocate(zrodlo%polozenia(pN,2) );
        allocate(zrodlo%Chi_m_in (pN,lM+lEvanMods));
        allocate(zrodlo%Chi_m_out(pN,lM+lEvanMods));
        allocate(zrodlo%k_m_in   (lM+lEvanMods));
        allocate(zrodlo%k_m_out  (lM+lEvanMods));
        allocate(zrodlo%deltamk  (lM+lEvanMods));
        allocate(zrodlo%deltaink (lM+lEvanMods));
        allocate(zrodlo%deltaoutk(lM+lEvanMods));
        allocate(zrodlo%ck       (lM+lEvanMods));
        allocate(zrodlo%dk       (lM+lEvanMods));
        allocate(zrodlo%Jin      (lM+lEvanMods));
        allocate(zrodlo%Jout     (lM+lEvanMods));
        allocate(zrodlo%Aij(lM+lEvanMods,lM))
        allocate(zrodlo%Sij(lM+lEvanMods,lM+lEvanMods))
        allocate(zrodlo%SijChiAuxMat(lM+lEvanMods,pN))
        allocate(zrodlo%Fj(pN))
        allocate(zrodlo%SijAijCkAuxVec(lM+lEvanMods))


        zrodlo%k_m_in  = 0
        zrodlo%k_m_out = 0
        zrodlo%Chi_m_in  = 0
        zrodlo%Chi_m_out = 0
        zrodlo%m_r  = 0
        zrodlo%m_t  = 0
        zrodlo%ck   = 0
        zrodlo%dk   = 0
        zrodlo%Jin  = 0
        zrodlo%Jout = 0
        zrodlo%bZaalokowane = .true.
    end subroutine zrodlo_alokuj_pamiec

    subroutine zrodlo_wypisz_info(zrodlo)
        class(czrodlo) :: zrodlo
        integer :: i
        if(TRANS_DEBUG==.true.) then
        print*,"Zrodlo: Wypisywanie informacji o zrodle:"
        print*,"    N       = ",zrodlo%N
        print*,"    L. modow= ",zrodlo%liczba_modow
        print*,"    L. evanm= ",zrodlo%liczba_evans
        print*,"    Kierunek= ",zrodlo%bKierunek
        print*,"    HNX     = ",zrodlo%hnX
        print*,"    HNY     = ",zrodlo%hnY
        write(*,"(A,2f10.4,A)"),"    R1(x,y) = (",zrodlo%r1,")"
        write(*,"(A,2f10.4,A)"),"    R2(x,y) = (",zrodlo%r2,")"
        print*,"    Mody wejsciowe/wyjsciowe:",zrodlo%liczba_modow
        do i = 1 , zrodlo%liczba_modow
            print"(A,i4,A,2f12.6,A,2f12.6)","K[",i,"][nm]:",zrodlo%k_m_in(i)*L2LR," | ",zrodlo%k_m_out(i)*L2LR
        enddo
        print*,"    Mody wejsciowe/wyjsciowe evanescentne:",zrodlo%liczba_evans
        do i = zrodlo%liczba_modow + 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
            print"(A,i4,A,2f12.6,A,2f12.6)","K[",i,"][nm]:",zrodlo%k_m_in(i)*L2LR," | ",zrodlo%k_m_out(i)*L2LR
        enddo
        endif
    endsubroutine zrodlo_wypisz_info



    subroutine zrodlo_wypisz_ckdk(zrodlo)
        class(czrodlo) :: zrodlo
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
    endsubroutine zrodlo_wypisz_ckdk


    subroutine zrodlo_wypisz_JinJout(zrodlo)
        class(czrodlo) :: zrodlo
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
            write(*,"(A,i4,A,f8.4,A,i4,A,f8.4)"),"   Jin(",i,")=",zrodlo%Jin(i),"  |  Jout(",i,")=",zrodlo%Jout(i)
        enddo
        print*, "! ----------------------------------------- !"
        write(*,"(A,f8.4,A,f8.4)"),"    SUMA(in)=",sum(zrodlo%Jin(:)),"  |   SUMA(out)=",sum(zrodlo%Jout(:))
        endif
    endsubroutine zrodlo_wypisz_JinJout

    subroutine zrodlo_zapisz_m_rt(zrodlo,mod)
        class(czrodlo) :: zrodlo
        integer,intent(in) :: mod ! numer modu wchodzacego
        integer            :: i
        do i = 1 , zrodlo%liczba_modow
            zrodlo%m_r(mod,i) = zrodlo%dk(i)
            zrodlo%m_t(mod,i) = zrodlo%ck(i)
        enddo
    endsubroutine zrodlo_zapisz_m_rt

    ! -------------------------------------------------------------------------
    ! Zapisuje mody do pliku o nazwie filename, podajemy rozwniez krok siatki w nm
    ! bSaveEvan - optionalnie mozemy zdecydowac czy chcemy zapisywac evanescente czy tylko
    ! fale plaskie, domyslnie mody evanescentne nie sa zapisywane
    ! -------------------------------------------------------------------------
    subroutine zrodlo_zapisz_mody(zrodlo,filename,dx,bSaveEvan)
        class(czrodlo) :: zrodlo
        character(*) :: filename
        double precision :: dx
        logical,optional,intent(in) :: bSaveEvan
        integer :: i,no_modow
        logical :: bSaveE

        print*,"Zrodlo: Zapisywanie modow do pliku:",filename
        open(4193,file=filename)
        if(present(bSaveEvan)) then
            bSaveE = bSaveEvan ! domyslnie nie bedziemy zapisywac modow evanescentym
        else
            bSaveE = .false.
        endif
        no_modow = zrodlo%liczba_modow

        do i = 1 , zrodlo%N
        if(bSaveE)           write(4193,"(300f20.8)"),zrodlo%polozenia(i,:)*dx,abs(zrodlo%Chi_m_in(i,:))**2,abs(zrodlo%Chi_m_out(i,:))**2
        if(.not. bSaveE)     write(4193,"(300f20.8)"),zrodlo%polozenia(i,:)*dx,abs(zrodlo%Chi_m_in(i,1:no_modow ))**2,abs(zrodlo%Chi_m_out(i,1:no_modow))**2
        enddo
        close(4193)
    end subroutine zrodlo_zapisz_mody


    ! --------------------------------------------------------------------
    ! Ustawiamy wczesniej zaalokowane zrodlo przez polecenie systemowe.
    ! Parametry:
    ! pY1,pYN - polozenia y na siatce numerycznej liczone od 1 do NY, w przypadku zrodla poziomego
    !           oznaczaja polozenia X1 oraz XN
    ! pX1 - polozenie X zrodla, dla zrodla poziomego polozenie Y1
    ! pDX,pEf,pBz - DX [nm] , Ef [meV] , BZ [T]
    ! pKierunek - enum ZRODLO_KIERUNEK_PRAWO/LEWO/GORA/DOL - ustala w ktora skierowane jest zrodlo
    ! pUTOTAL   - referencja do potencjalu ukladu
    ! --------------------------------------------------------------------
    subroutine zrodlo_ustaw(zrodlo,pY1,pYN,pX1,pDX,pEf,pBz,pKierunek,pUTOTAL)
        class(czrodlo)             ::  zrodlo
        integer,intent(in)         ::  pY1,pYN,pX1
        doubleprecision,intent(in) ::  pDx,pEf,pBz
        integer,intent(in)         ::  pKierunek ! enum
        double precision,dimension(:,:) :: pUTOTAL ! calkowity potencjal w [meV]
        double precision :: pUvec(pYN - pY1 + 1)

        double precision :: dx, Ef, Bz
        complex*16,dimension(:,:),allocatable :: tempB
        complex*16 :: ctmp
        integer :: N ! liczba oczek dla zrodla
        integer :: lModow,lModowEvan
        ! zmienne pomocnicze
        integer :: i,j,k,ntmp

        dx  = pdx*L2LR    ! konwertujemy do jednostek donorowych
        Ef  = pEf/1000.0/Rd
        BZ  = BtoDonorB(pBz)

        print*,"Okres B",DonorBtoB(2*3.14159/dx/dx) , "[T]"

        N = pYN - pY1 + 1; ! to jest tak samo liczone ale zmienne pYN , pY1 maja inna interpretacje



        if( pKierunek == ZRODLO_KIERUNEK_PRAWO .or. pKierunek == ZRODLO_KIERUNEK_LEWO ) then
            do i = 1 , N
                pUVEC(i) = pUTOTAL(pX1,pY1 + i - 1)
            enddo
            call modpop_calc_modes_from_wfm(pDX,N,pEf,pBz,pUVEC,.true.)
        else ! dla zrodel gora dol
            do i = 1 , N
                pUVEC(i) = pUTOTAL(pY1 + i - 1,pX1)
            enddo
            call modpop_calc_modes_from_wfm(pDX,N,pEf,pBz,pUVEC,.false.)
        endif

        call modpop_liczba_podpasm(lModow,lModowEvan)
        print*,"Liczba modow evan=",lModowEvan
!        lModowEvan = 0

        call zrodlo%zrodlo_alokuj_pamiec(N,lModow,lModowEvan)
        if(lModow > 0) then
        ! POBIERAMY MODY POPRZECZNE ORAZ ICH WEKTORY FALOWE
        ! pobieramy od razu z evanescentnymi wektorami
        call modpop_get_km (lModow+lModowEvan,zrodlo%k_m_in,zrodlo%k_m_out)
        call modpop_get_chi(lModow+lModowEvan,N,zrodlo%Chi_m_in,zrodlo%Chi_m_out)

        endif
        ! MODUL RELACJI DYSPERSJI JUZ NIE POTRZEBNY
        call modpop_zwalnienie_pamieci()


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
        else ! dla zrodel gora dol
            zrodlo%polozenia(:,2) = pX1
            do i = 1 , N
                zrodlo%polozenia(i,1) = pY1 + i - 1
            enddo
            zrodlo%r1             = (/pY1,pX1/)*DX
            zrodlo%r2             = (/pYN,pX1/)*DX
            zrodlo%hnx            = (pYN + pY1)/2.0*DX
        endif


        call zrodlo%zrodlo_wypisz_info()




        ntmp =  zrodlo%liczba_modow + zrodlo%liczba_evans
        select case (zrodlo%bKierunek)
        ! ---------------------------------------------------------------- !
        case (ZRODLO_KIERUNEK_PRAWO)

            ! Wzynaczanie wektora pomocniczego
            do k = 1 , ntmp
                zrodlo%deltaink(k)  = exp( +DX*(zrodlo%k_m_in(k)  + II*Bz*zrodlo%hny)) - exp( -DX*(zrodlo%k_m_in(k)  + II*Bz*zrodlo%hny))
                zrodlo%deltaoutk(k) = exp( +DX*(zrodlo%k_m_out(k) + II*Bz*zrodlo%hny)) - exp( -DX*(zrodlo%k_m_out(k) + II*Bz*zrodlo%hny))
            enddo
            ! ----------------------------------------------------------------------
            ! Macierze warunkow brzegowych dla zrodel skierowanych w prawo ------>
            ! Aij    = < X(-i) | X(+j) >
            ! Sij^-1 = < X(-i) | X(-j) >
            ! ----------------------------------------------------------------------
            do i = 1 , ntmp
            do j = 1 , lModow
                zrodlo%Aij(i,j) = sum( conjg(zrodlo%Chi_m_out(:,i))*zrodlo%Chi_m_in (:,j) )*DX
            end do
            end do

            allocate(tempB(ntmp,1))
            do i = 1 , ntmp
            do j = 1 , ntmp
                zrodlo%Sij(i,j) = sum( conjg(zrodlo%Chi_m_out(:,i))*zrodlo%Chi_m_out(:,j) )*DX
            end do
            end do
            tempB = 1
            call zgaussj(zrodlo%Sij(1:ntmp,1:ntmp),ntmp,ntmp,tempB(1:ntmp,1),1,1)
            deallocate(tempB)

            ! Wyznaczanie macierzy pomocniczej zawierajacej iloczyny macierzy Sij i wektora Chi
            do i = 1 , ntmp
            do j = 1 , N
                zrodlo%SijChiAuxMat(i,j) = sum( zrodlo%Sij(i,:)*conjg(zrodlo%Chi_m_out(j,:)) )
            end do
            end do


        ! ---------------------------------------------------------------- !
        case (ZRODLO_KIERUNEK_LEWO)
            ! Wzynaczanie wektora pomocniczego
            do k = 1 , ntmp
                zrodlo%deltaink(k)  = exp( +DX*(zrodlo%k_m_in(k)  + II*Bz*zrodlo%hny)) - exp( -DX*(zrodlo%k_m_in(k)  + II*Bz*zrodlo%hny))
                zrodlo%deltaoutk(k) = exp( +DX*(zrodlo%k_m_out(k) + II*Bz*zrodlo%hny)) - exp( -DX*(zrodlo%k_m_out(k) + II*Bz*zrodlo%hny))
            enddo

            ! ----------------------------------------------------------------------
            ! Macierze warunkow brzegowych dla zrodel skierowanych w lewo  <-------
            ! Aij    = < X(+i) | X(-j) >
            ! Sij^-1 = < X(+i) | X(+j) >
            ! ----------------------------------------------------------------------

            do i = 1 , ntmp
            do j = 1 , lModow
                zrodlo%Aij(i,j) = sum( conjg(zrodlo%Chi_m_in(:,i))*zrodlo%Chi_m_out(:,j) )*DX
            end do
            end do

            allocate(tempB(ntmp,1))
            do i = 1 , ntmp
            do j = 1 , ntmp
                zrodlo%Sij(i,j) = sum( conjg(zrodlo%Chi_m_in(:,i))*zrodlo%Chi_m_in (:,j) )*DX
            end do
            end do
            tempB = 1
            call zgaussj(zrodlo%Sij(1:ntmp,1:ntmp),ntmp,ntmp,tempB(1:ntmp,1),1,1)
            deallocate(tempB)

            ! Wyznaczanie macierzy pomocniczej zawierajacej iloczyny macierzy Sij i wektora Chi
            do i = 1 , ntmp
            do j = 1 , N
                zrodlo%SijChiAuxMat(i,j) = sum( zrodlo%Sij(i,:)*conjg(zrodlo%Chi_m_in(j,:)) )
            end do
            end do
            ! ---------------------------------------------------------------- !
        case (ZRODLO_KIERUNEK_GORA)
            print*,"MODZRODLO: Zrodlo wejsciowe nie obslugiwane!"
            stop


            ! Wzynaczanie wektora pomocniczego
            zrodlo%deltamk(:) = 0 ! dajemy zero bo nie bedzie uzywane, bo tutaj cechowanie zalezy jeszcze od x

            ! ----------------------------------------------------------------------
            ! Macierze warunkow brzegowych dla zrodel skierowanych w gore ------>
            ! Aij    = < X(-i) | X(+j) >
            ! Sij^-1 = < X(-i) | X(-j) >
            ! ----------------------------------------------------------------------
            do i = 1 , ntmp
            do j = 1 , lModow
                zrodlo%Aij(i,j) = 0
                do k = 1 , zrodlo%N


                    ctmp = exp( (zrodlo%polozenia(k,2))*DX*(&
                            &II*zrodlo%k_m_in(j) +&
                            &II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )


                    zrodlo%Aij(i,j) = zrodlo%Aij(i,j) + ctmp*( conjg(zrodlo%Chi_m_out(k,i))*zrodlo%Chi_m_in (k,j) )*DX
                enddo

                !zrodlo%Aij(i,j) = sum( conjg(zrodlo%Chi_m_out(:,i))*zrodlo%Chi_m_in (:,j) )*DX
            end do
            end do

            allocate(tempB(ntmp,1))
            do i = 1 , ntmp
            do j = 1 , ntmp
                zrodlo%Sij(i,j) = 0
                do k = 1 , zrodlo%N
                    if( j <= lModow) then
                    ctmp = exp( (zrodlo%polozenia(k,2))*DX*(-II*zrodlo%k_m_in(j) +&
                            &II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )

                    else
                    ctmp = exp( zrodlo%polozenia(k,2)*DX*(&
                            &-zrodlo%k_m_in(j) +&
                            &II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
                    endif

!                    ctmp = 1

                    zrodlo%Sij(i,j) = zrodlo%Sij(i,j) + ctmp*( conjg(zrodlo%Chi_m_out(k,i))*zrodlo%Chi_m_out(k,j) )*DX
                enddo

            end do
            end do
            tempB = 1
            call zgaussj(zrodlo%Sij(1:ntmp,1:ntmp),ntmp,ntmp,tempB(1:ntmp,1),1,1)
            deallocate(tempB)

            ! Wyznaczanie macierzy pomocniczej zawierajacej iloczyny macierzy Sij i wektora Chi
            do i = 1 , ntmp
            do j = 1 , N
                zrodlo%SijChiAuxMat(i,j) = sum( zrodlo%Sij(i,:) * conjg(zrodlo%Chi_m_out(j,:)) )
            end do
            end do


            ! ---------------------------------------------------------------- !
            case (ZRODLO_KIERUNEK_DOL)
            print*,"MODZRODLO: Zrodlo wejsciowe nie obslugiwane!"
            stop
            ! Wzynaczanie wektora pomocniczego
            zrodlo%deltamk(:) = 0 ! dajemy zero bo nie bedzie uzywane, bo tutaj cechowanie zalezy jeszcze od x

            ! ----------------------------------------------------------------------
            ! Macierze warunkow brzegowych dla zrodel skierowanych w dol  <-------
            ! Aij    = < X(+i) | X(-j) >
            ! Sij^-1 = < X(+i) | X(+j) >
            ! ----------------------------------------------------------------------
            do i = 1 , ntmp
            do j = 1 , lModow
                zrodlo%Aij(i,j) = 0
                do k = 1 , zrodlo%N
!                    ctmp = exp( zrodlo%polozenia(k,2)*DX*(&
!                            &-II*zrodlo%k_m_in(j) +&
!                            & II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )

                    ctmp = exp( zrodlo%polozenia(k,2)*DX*(-2*II*zrodlo%k_m_in(j) ) )

                    zrodlo%Aij(i,j) = zrodlo%Aij(i,j) + ctmp*( conjg(zrodlo%Chi_m_in(k,i))*zrodlo%Chi_m_out (k,j) )*DX
                enddo


            end do
            end do

            allocate(tempB(ntmp,1))
            do i = 1 , ntmp
            do j = 1 , ntmp
                zrodlo%Sij(i,j) = 0
                do k = 1 , zrodlo%N
!                    if( j <= lModow) then
!                    ctmp = exp( (zrodlo%polozenia(k,2))*DX*( &
!                            &+II*zrodlo%k_m_in(j) + &
!                            & II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
!                    else
!                    ctmp = exp( (zrodlo%polozenia(k,2))*DX*( &
!                            & zrodlo%k_m_out(j) + &
!                            & II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
!                    endif

                    if( j <= lModow) then
                    ctmp = exp( (zrodlo%polozenia(k,2))*DX*(II*zrodlo%k_m_in(j) +&
                            &II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )

                    else
                    ctmp = exp( zrodlo%polozenia(k,2)*DX*(&
                            &-zrodlo%k_m_in(j) +&
                            &II*Bz*( zrodlo%polozenia(k,1)*DX - zrodlo%hnX ) ) )
                    endif


                    zrodlo%Sij(i,j) = zrodlo%Sij(i,j) + ctmp*( conjg(zrodlo%Chi_m_in(k,i))*zrodlo%Chi_m_in(k,j) )*DX

                enddo

            end do
            end do
            tempB = 1
            call zgaussj(zrodlo%Sij(1:ntmp,1:ntmp),ntmp,ntmp,tempB(1:ntmp,1),1,1)
            deallocate(tempB)

            ! Wyznaczanie macierzy pomocniczej zawierajacej iloczyny macierzy Sij i wektora Chi
            do i = 1 , ntmp
            do j = 1 , N
                zrodlo%SijChiAuxMat(i,j) = sum( zrodlo%Sij(i,:) * conjg(zrodlo%Chi_m_in(j,:)) )
            end do
            end do


        case default
            print*,"Modzrodlo:: Nie ma takiego typu zrodla jak",zrodlo%bKierunek
            stop
        endselect



    end subroutine zrodlo_ustaw



    ! --------------------------------------------------------------------
    ! Oblicza wartosci dla wektora wyrazow wolnych jak rozwiazywany jest
    ! uklad rownan na Psi(u,v). Jako argument podajemy jedynie dx w
    ! jednostkach donorowych. Utworzony wektor przechowywany jest w tablicy
    ! Fj(v). Procedura uzupelnia wektor wyrazow wolnych w zaleznosci
    ! od tego jaki jest kierunek zrodla. Przed wywolaniem tej procedury
    ! nalezy ustawic odpowiednie wartosci amplitud wejsciowych ck(i).
    ! --------------------------------------------------------------------
    subroutine zrodlo_oblicz_Fj(zrodlo,pdx)
        class(czrodlo)  ::  zrodlo
        doubleprecision :: pdx
        integer         :: i , pi, pj, k , p, q
        complex*16      :: post , Xkn , deltapk , deltamk , kvec
        doubleprecision :: dx, Ef, Bz , ypos , xpos , bpart
        dx  = pdx
        Ef  = atomic_Ef/1000.0/Rd
        BZ  = BtoDonorB(atomic_Bz)

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
                zrodlo%Fj(i) = 0
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
                case (ZRODLO_KIERUNEK_PRAWO)
                    post             = -(0.5/DX/DX)*EXP(II*DX*DX*pj*BZ)
                    do k = 1 , zrodlo%liczba_modow
                        zrodlo%Fj(i) = zrodlo%Fj(i) + zrodlo%ck(k)*zrodlo%deltaink(k)*zrodlo%Chi_m_in (i,k)
                    enddo
                    do k = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
                        Xkn          = zrodlo%SijAijCkAuxVec(k)
                        deltamk      = zrodlo%deltaoutk(k)
                        zrodlo%Fj(i) = zrodlo%Fj(i) - Xkn*deltamk*zrodlo%Chi_m_out(i,k)
                    enddo
                ! ---------------------------------------------------------------- !
                ! Wypelniamy wektor wyrazow wolnych dla zrodej skierowanych w lewo
                !
                !                        <---------------
                !
                case (ZRODLO_KIERUNEK_LEWO)
                        post         = (0.5/DX/DX)*EXP(-II*DX*DX*pj*BZ)
                    do k = 1 , zrodlo%liczba_modow
                        zrodlo%Fj(i) = zrodlo%Fj(i) + zrodlo%ck(k)*zrodlo%deltaoutk(k)*zrodlo%Chi_m_out(i,k)
                    enddo
                    do k = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
                        Xkn          = zrodlo%SijAijCkAuxVec(k)
                        deltamk      = zrodlo%deltaink(k)
                        zrodlo%Fj(i) = zrodlo%Fj(i) - Xkn*deltamk*zrodlo%Chi_m_in(i,k)
                    enddo

                ! ---------------------------------------------------------------- !
                ! Wypelniamy wektor wyrazow wolnych dla zrodej skierowanych w gore
                !
                !                        --------------->
                !
                case (ZRODLO_KIERUNEK_GORA)
                    print*,"MODZRODLO: Zrodlo wejsciowe nie obslugiwane!"
                    stop
                    post             = -(0.5/DX/DX)
                    bpart            = Bz*( xpos  - zrodlo%hnx )
                    do k = 1 , zrodlo%liczba_modow

                        kvec         = zrodlo%k_m_in(k)
                        deltamk      = exp( (ypos +  DX)*(+II*kvec + II*bpart) ) &
                                    &- exp( (ypos -  DX)*(+II*kvec + II*bpart) )
                        zrodlo%Fj(i) = zrodlo%Fj(i) + zrodlo%ck(k)*deltamk*zrodlo%Chi_m_in (i,k)

                    enddo
                    do k = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans

                        if( k <= zrodlo%liczba_modow) then
                            deltamk = exp( (ypos +  DX)*(-II*zrodlo%k_m_in(k) + II*bpart) ) - exp(  (ypos-DX)*(-II*zrodlo%k_m_in(k) + II*bpart) )
                        else
                            deltamk = exp( (ypos + DX)*(-zrodlo%k_m_in(k) + II*bpart)  ) - exp( (ypos-DX)*( -zrodlo%k_m_in(k) + II*bpart) )
                        endif

                        Xkn          = zrodlo%SijAijCkAuxVec(k)
                        zrodlo%Fj(i) = zrodlo%Fj(i) - Xkn*deltamk*zrodlo%Chi_m_out(i,k)
                    enddo
                case (ZRODLO_KIERUNEK_DOL)
                    print*,"MODZRODLO: Zrodlo wejsciowe nie obslugiwane!"
                    stop
                endselect
                zrodlo%Fj(i) = zrodlo%Fj(i)*post

        enddo ! end of do i=1, zrodlo%N

    end subroutine zrodlo_oblicz_Fj



    ! --------------------------------------------------------------------
    ! Oblicza wyraz alpha(v,i) dany wzrorem:
    !
    !   alpha(v,i) = DX * sum_k,p {  delta_k * Chi(k,v) * S(k,p) * Chi^*(p,i)  }
    !
    !   gdzie: delta_k = 2i*sin(kDX + DX*B*H).
    ! Funkcja ta potrzebna jest do poprawnego wypelnienia macierzy H podczas
    ! tworzenia ukladu rownan. Wiecej szczegolow w notatkach.
    ! ----------------------------------------------------------------------
    complex*16 function zrodlo_alfa_v_i(zrodlo,pdx,v,i)
        class(czrodlo)  :: zrodlo
        doubleprecision :: pdx
        integer         :: v,i
        doubleprecision :: dx ,Ef ,Bz , bpart , ypos
        integer         :: k
        complex*16      :: Xkn , deltamk

        dx  = pdx
        Ef  = atomic_Ef/1000.0/Rd
        BZ  = BtoDonorB(atomic_Bz)

        Xkn  = 0

        select case (zrodlo%bKierunek)
        ! -----------------------------------------------------------------
        !       Wyznaczamy alpha(v,i) dla wejsc skierowanych w prawo
        !
        !                       --------------------->
        !
        ! -----------------------------------------------------------------
        case (ZRODLO_KIERUNEK_PRAWO)

            do k = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
                Xkn = Xkn + zrodlo%deltaoutk(k)*zrodlo%Chi_m_out(v,k)*zrodlo%SijChiAuxMat(k,i);
            enddo
        ! -----------------------------------------------------------------
        !       Wyznaczamy alpha(v,i) dla wejsc skierowanych w lewo
        !
        !                       <---------------------
        !
        ! -----------------------------------------------------------------
        case (ZRODLO_KIERUNEK_LEWO)
            do k = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
                Xkn = Xkn + zrodlo%deltaink(k)*zrodlo%Chi_m_in(v,k)*zrodlo%SijChiAuxMat(k,i);
            enddo

        ! -----------------------------------------------------------------
        !       Wyznaczamy alpha(v,i) dla wejsc skierowanych w gore
        !
        !                       --------------------->
        !
        ! -----------------------------------------------------------------

        case (ZRODLO_KIERUNEK_GORA)
            print*,"MODZRODLO: Zrodlo wejsciowe nie obslugiwane!"
            stop
            bpart = Bz*( zrodlo%polozenia(v,1)*dx -  zrodlo%hnx)
            ypos  = (zrodlo%polozenia(v,2))*dx
            do k = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
                deltamk = 0
                if( k <= zrodlo%liczba_modow) then
                    deltamk = exp( (ypos + DX)*(-II*zrodlo%k_m_in(k) + II*bpart) ) - exp(  (ypos-DX)*(-II*zrodlo%k_m_in(k) + II*bpart) )
                else
                    deltamk = exp( (ypos + DX)*(-zrodlo%k_m_in(k) + II*bpart)  )   - exp( (ypos-DX)*( -zrodlo%k_m_in(k) + II*bpart) )
                endif

!                deltamk = exp( (DX)*(-II*zrodlo%k_m_in(k)) ) - exp(  (-DX)*(-II*zrodlo%k_m_in(k) ) )

                Xkn = Xkn + deltamk*zrodlo%Chi_m_out(v,k) * zrodlo%SijChiAuxMat(k,i);

            enddo

        ! -----------------------------------------------------------------
        !       Wyznaczamy alpha(v,i) dla wejsc skierowanych w dol
        !
        !                       <---------------------
        !
        ! -----------------------------------------------------------------
        case (ZRODLO_KIERUNEK_DOL)
            print*,"MODZRODLO: Zrodlo wejsciowe nie obslugiwane!"
            stop
            bpart = Bz*( zrodlo%polozenia(v,1)*dx -  zrodlo%hnx)
            ypos  = (zrodlo%polozenia(v,2))*dx
            do k = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
                deltamk = 0
                if( k <= zrodlo%liczba_modow) then
                    deltamk = exp( (ypos + DX)*(II*zrodlo%k_m_in(k) + II*bpart) ) - exp(  (ypos-DX)*(II*zrodlo%k_m_in(k) + II*bpart) )
                else
                    deltamk = exp( (ypos + DX)*(zrodlo%k_m_in(k) + II*bpart)  )   - exp( (ypos-DX)*( zrodlo%k_m_in(k) + II*bpart) )
                endif

                Xkn = Xkn + deltamk*zrodlo%Chi_m_in(v,k) * zrodlo%SijChiAuxMat(k,i);
            enddo

        endselect

        zrodlo_alfa_v_i = Xkn*dx

    end function zrodlo_alfa_v_i


    subroutine zrodlo_oblicz_dk(zrodlo,VPHI,GINDEX,dx)
        class(czrodlo)          :: zrodlo
        complex*16,dimension(:) :: VPHI
        integer,dimension(:,:)  :: GINDEX
        double precision        :: dx

        integer :: p,q,k,i,pi,pj
        complex*16 :: dk , alpha_p, beta_p


        do k = 1 , zrodlo%liczba_modow
            dk = 0
            do p = 1 , zrodlo%liczba_modow + zrodlo%liczba_evans
                alpha_p = 0
                do i = 2 , zrodlo%N-1
                    pi   = zrodlo%polozenia(i,1)
                    pj   = zrodlo%polozenia(i,2)

                    select case (zrodlo%bKierunek)
                    ! ------------------------------------------------------------------
                    !
                    !                          <----------------
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_PRAWO,ZRODLO_KIERUNEK_GORA)
                        alpha_p = alpha_p + dx*conjg(zrodlo%Chi_m_out(i,p))*VPHI(GINDEX(pi,pj))
                    ! ------------------------------------------------------------------
                    !
                    !                          ---------------->
                    !
                    ! ------------------------------------------------------------------
                    case (ZRODLO_KIERUNEK_LEWO,ZRODLO_KIERUNEK_DOL)
                        alpha_p = alpha_p + dx*conjg(zrodlo%Chi_m_in(i,p))*VPHI(GINDEX(pi,pj))
                    endselect

                enddo
                beta_p = 0
                do q = 1 , zrodlo%liczba_modow
                    beta_p = beta_p + zrodlo%Aij(p,q)*zrodlo%ck(q)
                enddo
                dk = dk + zrodlo%Sij(k,p)*( alpha_p - beta_p )
            enddo ! end do p
            zrodlo%dk(k) = dk

        enddo ! end do k

        call zrodlo%zrodlo_wypisz_ckdk()

    end subroutine zrodlo_oblicz_dk


    ! wzory brane z mojej pracy inzynierskiej (strona 12 - tabela)
    subroutine zrodlo_oblicz_JinJout(zrodlo,dx)
        class(czrodlo)    :: zrodlo
        double precision  :: dx

        integer          :: k , pj , i
        double precision :: Jin , Jout , Bz ,kvec

        BZ  = BtoDonorB(atomic_Bz)



        do k = 1 , zrodlo%liczba_modow
            Jin  = 0
            Jout = 0
            kvec = abs(zrodlo%k_m_in(k))

            select case (zrodlo%bKierunek)
            ! ------------------------------------------------------------------
            !
            !                          ---------------->
            !                          <----------------
            !
            ! ------------------------------------------------------------------
            case (ZRODLO_KIERUNEK_PRAWO,ZRODLO_KIERUNEK_LEWO)
                do i = 1 , zrodlo%N
                    pj    = zrodlo%polozenia(i,2)
                    Jin   = Jin  + abs(zrodlo%Chi_m_in (i,k))**2*sin( -dx*dx*pj*Bz + dx*(+kvec + zrodlo%hnY*Bz) )
                    Jout  = Jout + abs(zrodlo%Chi_m_out(i,k))**2*sin( -dx*dx*pj*Bz + dx*(-kvec + zrodlo%hnY*Bz) )
                enddo ! end of i
            case (ZRODLO_KIERUNEK_GORA,ZRODLO_KIERUNEK_DOL)
                do i = 1 , zrodlo%N
                    pj    = zrodlo%polozenia(i,1)
                    Jin   = Jin  + abs(zrodlo%Chi_m_in (i,k))**2*sin( +dx*dx*pj*Bz + dx*(+kvec - zrodlo%hnX*Bz) )
                    Jout  = Jout + abs(zrodlo%Chi_m_out(i,k))**2*sin( +dx*dx*pj*Bz + dx*(-kvec - zrodlo%hnX*Bz) )
                enddo ! end of i
            endselect

            Jin   = - Jin  * dx * abs(zrodlo%ck(k))**2
            Jout  = - Jout * dx * abs(zrodlo%dk(k))**2
            zrodlo%Jin(k)  = Jin
            zrodlo%Jout(k) = Jout
        enddo ! end of k


    end subroutine

! ---------------------------------------------------------------------------------------------
!
!   Funkcje klasy abszrodla
!
! ---------------------------------------------------------------------------------------------
    subroutine abs_zrodlo_zwolnij_pamiec(zrodlo)
        class(cabs_zrodlo) :: zrodlo

        if(TRANS_DEBUG==.true.) print*,"Abszrodlo: Zwalanianie pamieci"
        if(allocated(zrodlo%polozenia)) deallocate(zrodlo%polozenia);


    end subroutine abs_zrodlo_zwolnij_pamiec


! --------------------------------------------------------------------------------- !
! Podobnie jak to jest dla normalnej klasy zrodlo. Potrzebna jest tylko informacja
! o energii fermiego
! --------------------------------------------------------------------------------- !

    subroutine abs_zrodlo_ustaw(zrodlo,pY1,pYN,pX1,pEf,pKierunek)
        class(cabs_zrodlo)             ::  zrodlo
        integer,intent(in)         ::  pY1,pYN,pX1
        doubleprecision,intent(in) ::  pEf
        integer,intent(in)         ::  pKierunek ! enum

        double precision :: Ef
        integer :: N,i

        print*,"Dodawanie zrodla z transparentnymi w.b"

        Ef  = pEf/1000.0/Rd

        N = pYN - pY1 + 1; ! to jest tak samo liczone ale zmienne pYN , pY1 maja inna interpretacje
        zrodlo%kvec = sqrt(2*Ef)
        zrodlo%N    = N

        if(allocated(zrodlo%polozenia)) deallocate(zrodlo%polozenia)

        allocate(zrodlo%polozenia(N,2))
        zrodlo%polozenia = 0
        ! ustawianie parametrow zrodel
        zrodlo%bKierunek = pKierunek

        if( pKierunek == ZRODLO_KIERUNEK_PRAWO .or. pKierunek == ZRODLO_KIERUNEK_LEWO ) then
            zrodlo%polozenia(:,1) = pX1
            do i = 1 , N
                zrodlo%polozenia(i,2) = pY1 + i - 1
            enddo
        else ! dla zrodel gora dol
            zrodlo%polozenia(:,2) = pX1
            do i = 1 , N
                zrodlo%polozenia(i,1) = pY1 + i - 1
            enddo
        endif

        print*,"kvec      = " , zrodlo%kvec * L2LR , "[nm]"
        print*,"kierunek  ="  , zrodlo%bKierunek
        print*,"N         ="  , zrodlo%N


    end subroutine abs_zrodlo_ustaw

! --------------------------------------------------------------------------------- !
! Kopuje dane z jednego zrodla do drugiego
! --------------------------------------------------------------------------------- !

    subroutine abs_zrodlo_skopiuj(zrodlo,od_zrodla)
        class(cabs_zrodlo)       ::  zrodlo,od_zrodla

        zrodlo%bKierunek = od_zrodla%bKierunek
        zrodlo%kvec      = od_zrodla%kvec
        zrodlo%N         = od_zrodla%N
        if(allocated(zrodlo%polozenia)) deallocate(zrodlo%polozenia)
        allocate(zrodlo%polozenia(zrodlo%N,2))

        zrodlo%polozenia = od_zrodla%polozenia
    end subroutine abs_zrodlo_skopiuj

endmodule modzrodlo
