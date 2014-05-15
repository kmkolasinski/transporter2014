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
        doubleprecision,dimension(:),allocatable  :: Jin,Jout ! prady wejsciowe zwiazane z ck i wyjsciowe - dk
        integer                                   :: rozbieg  ! zakres rzutowania w procesie wyznaczania ck, dk
        complex*16,dimension(200,200)             :: m_r,m_t  ! macierz wspolczynnikow odbicia i transmisji
        complex*16,dimension(:,:),allocatable     :: Aij
        complex*16,dimension(:,:),allocatable     :: Sij
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
        if(allocated(zrodlo%ck))        deallocate(zrodlo%ck);
        if(allocated(zrodlo%dk))        deallocate(zrodlo%dk);
        if(allocated(zrodlo%Jin))       deallocate(zrodlo%Jin);
        if(allocated(zrodlo%Jout))      deallocate(zrodlo%Jout);
        if(allocated(zrodlo%Aij))       deallocate(zrodlo%Aij);
        if(allocated(zrodlo%Sij))       deallocate(zrodlo%Sij);
        if(allocated(zrodlo%Fj))        deallocate(zrodlo%Fj);

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
        allocate(zrodlo%ck       (lM+lEvanMods));
        allocate(zrodlo%dk       (lM+lEvanMods));
        allocate(zrodlo%Jin      (lM+lEvanMods));
        allocate(zrodlo%Jout     (lM+lEvanMods));
        allocate(zrodlo%Aij(lM,lM))
        allocate(zrodlo%Sij(lM,lM))
        allocate(zrodlo%Fj(pN))

        zrodlo%bZaalokowane = .true.
    end subroutine zrodlo_alokuj_pamiec

    subroutine zrodlo_wypisz_info(zrodlo)
        class(czrodlo) :: zrodlo
        integer :: i
        if(TRANS_DEBUG==.true.) then
        print*,"Zrodlo: Wypisywanie informacji o zrodle:"
        print*,"    N       = ",zrodlo%N
        print*,"    L. modow= ",zrodlo%liczba_modow
        print*,"    Rozbieg = ",zrodlo%rozbieg
        print*,"    Kierunek= ",zrodlo%bKierunek
        print*,"    HNY     = ",zrodlo%hnY
        write(*,"(A,2f10.4,A)"),"    R1(x,y) = (",zrodlo%r1,")"
        write(*,"(A,2f10.4,A)"),"    R2(x,y) = (",zrodlo%r2,")"
        print*,"    Mody wejsciowe:"
        do i = 1 , zrodlo%liczba_modow
            print*,"        Kin [",i,"]=",zrodlo%k_m_in(i)*L2LR
        enddo
        print*,"    Mody wyjsciowe:"
        do i = 1 , zrodlo%liczba_modow
            print*,"        Kout[",i,"]=",zrodlo%k_m_out(i)*L2LR
        enddo
        endif
    endsubroutine zrodlo_wypisz_info



    subroutine zrodlo_wypisz_ckdk(zrodlo)
        class(czrodlo) :: zrodlo
        integer :: i
        if(TRANS_DEBUG==.true.) then
        print*, "! ----------------------------------------- !"
        if(zrodlo%bKierunek == .true. ) then
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
        if(zrodlo%bKierunek == .true. ) then
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
    ! -------------------------------------------------------------------------
    subroutine zrodlo_zapisz_mody(zrodlo,filename,dx)
        class(czrodlo) :: zrodlo
        character(*) :: filename
        double precision :: dx
        integer :: i
        print*,"Zrodlo: Zapisywanie modow do pliku:",filename
        open(4193,file=filename)
        do i = 1 , zrodlo%N
            write(4193,"(300e20.8)"),zrodlo%polozenia(i,2)*dx,abs(zrodlo%Chi_m_in(i,:))**2,abs(zrodlo%Chi_m_out(i,:))**2
        enddo
        close(4193)
    end subroutine zrodlo_zapisz_mody


    ! --------------------------------------------------------------------
    ! Ustawiamy wczesniej zaalokowane zrodlo przez polecenie systemowe.
    ! Parametry:
    ! pY1,pYN - polozenia y na siatce numerycznej liczone od 1 do NY
    ! pX1 - polozenie X zrodla
    ! pDX,pEf,pBz - DX [nm] , Ef [meV] , BZ [T]
    ! pKierunek - enum ZRODLO_KIERUNEK_PRAWO/LEWO - ustala w ktora skierowane jest zrodlo
    ! pWejscie  - czy zrodlo jest wejsciowe, czy wyjsciowe (true,false)
    ! pRozbieg  - liczba oczek siatki brana do liczenia amplitud rozpraszania
    ! pUTOTAL   - referencja do potencjalu ukladu
    ! --------------------------------------------------------------------
    subroutine zrodlo_ustaw(zrodlo,pY1,pYN,pX1,pDX,pEf,pBz,pKierunek,pRozbieg,pUTOTAL)
        class(czrodlo)             ::  zrodlo
        integer,intent(in)         ::  pY1,pYN,pX1
        doubleprecision,intent(in) ::  pDx,pEf,pBz
        integer,intent(in)         ::  pKierunek ! enum
        integer,intent(in)         ::  pRozbieg
        double precision,dimension(:,:) :: pUTOTAL ! calkowity potencjal w [meV]
        double precision :: pUvec(pYN - pY1 + 1)
        double precision :: dx, Ef, Bz
        complex*16,dimension(:,:),allocatable :: tempB
        integer :: N ! liczba oczek dla zrodla
        integer :: lModow
        ! zmienne pomocnicze
        integer :: i,j

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
        call zrodlo%zrodlo_alokuj_pamiec(N,lModow,0)

        zrodlo%polozenia(:,1) = pX1
        do i = 1 , N
            zrodlo%polozenia(i,2) = pY1 + i - 1
        enddo
        ! POBIERAMY MODY POPRZECZNE ORAZ ICH WEKTORY FALOWE
        call modpop_get_km(lModow,zrodlo%k_m_in,zrodlo%k_m_out)
        call modpop_get_chi(lModow,N,zrodlo%Chi_m_in,zrodlo%Chi_m_out)

        ! ustawianie parametrow zrodel
        if(pKierunek == ZRODLO_KIERUNEK_PRAWO) then
            zrodlo%bKierunek      = .true.
        else if(pKierunek == ZRODLO_KIERUNEK_LEWO) then
            zrodlo%bKierunek      = .false.
        else
            print*,"Zrodlo: Error nie znany kierunek zrodla. Do wyboru: PRAWO/LEWO"
            stop
        endif

        zrodlo%r1             = (/pX1,pY1/)*DX
        zrodlo%r2             = (/pX1,pYN/)*DX
        zrodlo%hny            = (pYN + pY1)/2.0*DX
        zrodlo%rozbieg        = pRozbieg
        zrodlo%m_r  = 0
        zrodlo%m_t  = 0
        zrodlo%ck   = 0
        zrodlo%dk   = 0
        zrodlo%Jin  = 0
        zrodlo%Jout = 0
        call zrodlo%zrodlo_wypisz_info()

        ! MODUL RELACJI DYSPERSJI JUZ NIE POTRZEBNY
        call modpop_zwalnienie_pamieci()


        if(zrodlo%bKierunek == .true.) then
        ! ----------------------------------------------------------------------
        ! Macierze warunkow brzegowych dla zrodel skierowanych w praweo ------>
        ! Aij    = < X(-i) | X(+j) >
        ! Sij^-1 = < X(-i) | X(-j) >
        ! ----------------------------------------------------------------------
        allocate(tempB(lModow,1))
        do i = 1 , lModow
        do j = 1 , lModow
            zrodlo%Aij(i,j) = sum( conjg(zrodlo%Chi_m_out(:,i))*zrodlo%Chi_m_in (:,j) )*DX
            zrodlo%Sij(i,j) = sum( conjg(zrodlo%Chi_m_out(:,i))*zrodlo%Chi_m_out(:,j) )*DX
        end do
        end do
        tempB = 1
        call zgaussj(zrodlo%Sij,lModow,lModow,tempB,1,1)
        deallocate(tempB)

        else
        ! ----------------------------------------------------------------------
        ! Macierze warunkow brzegowych dla zrodel skierowanych w lewo  <-------
        ! Aij    = < X(+i) | X(-j) >
        ! Sij^-1 = < X(+i) | X(+j) >
        ! ----------------------------------------------------------------------
        allocate(tempB(lModow,1))
        do i = 1 , lModow
        do j = 1 , lModow
            zrodlo%Aij(i,j) = sum( conjg(zrodlo%Chi_m_in(:,i))*zrodlo%Chi_m_out(:,j) )*DX
            zrodlo%Sij(i,j) = sum( conjg(zrodlo%Chi_m_in(:,i))*zrodlo%Chi_m_in (:,j) )*DX
        end do
        end do
        tempB = 1
        call zgaussj(zrodlo%Sij,lModow,lModow,tempB,1,1)
        deallocate(tempB)
        endif ! end of if jaki kierunek zrodla



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
        complex*16      :: post , Xkn , deltapk , deltamk
        doubleprecision :: dx, Ef, Bz , kvec
        dx  = pdx
        Ef  = atomic_Ef/1000.0/Rd
        BZ  = BtoDonorB(atomic_Bz)


        if(zrodlo%bKierunek == .true.) then
        ! -------------------------------------------------------------------
        ! Wypelniamy wektor wyrazow wolnych dla zrodej skierowanych w prawo
        !
        !                        --------------->
        !
        do i = 1 , zrodlo%N
                zrodlo%Fj(i) = 0
                pi         =  zrodlo%polozenia(i,1)
                pj         =  zrodlo%polozenia(i,2)
                post       = -(0.5/DX/DX)*EXP(II*DX*DX*pj*BZ)

                do k = 1 , zrodlo%liczba_modow
                    Xkn = 0
                    do p = 1 , zrodlo%liczba_modow
                    do q = 1 , zrodlo%liczba_modow
                        Xkn = Xkn + zrodlo%Sij(k,p)*zrodlo%Aij(p,q)*zrodlo%ck(q)
                    enddo
                    enddo
                    kvec   = zrodlo%k_m_in(k)
                    deltapk = 2*II*sin(+kvec*DX + DX*Bz*zrodlo%hny )
                    deltamk = 2*II*sin(-kvec*DX + DX*Bz*zrodlo%hny )
                    zrodlo%Fj(i) = zrodlo%Fj(i)  &
                         &   + zrodlo%ck(k)*deltapk*zrodlo%Chi_m_in (i,k)  &
                         &   - Xkn*       deltamk*zrodlo%Chi_m_out(i,k)
                     !print*,i,deltapk/2/DX
                enddo
                zrodlo%Fj(i) = zrodlo%Fj(i)*post
        enddo ! end of do i=1

        else
        ! -------------------------------------------------------------------
        ! Wypelniamy wektor wyrazow wolnych dla zrodej skierowanych w lewo
        !
        !                        <---------------
        !
        do i = 1 , zrodlo%N
                zrodlo%Fj(i) = 0
                pi         =  zrodlo%polozenia(i,1)
                pj         =  zrodlo%polozenia(i,2)
                post       = (0.5/DX/DX)*EXP(-II*DX*DX*pj*BZ)

                do k = 1 , zrodlo%liczba_modow
                    Xkn = 0
                    do p = 1 , zrodlo%liczba_modow
                    do q = 1 , zrodlo%liczba_modow
                        Xkn = Xkn + zrodlo%Sij(k,p)*zrodlo%Aij(p,q)*zrodlo%ck(q)
                    enddo
                    enddo
                    kvec   = zrodlo%k_m_in(k)
                    deltapk = 2*II*sin(+kvec*DX + DX*Bz*zrodlo%hny )
                    deltamk = 2*II*sin(-kvec*DX + DX*Bz*zrodlo%hny )
                    zrodlo%Fj(i) = zrodlo%Fj(i)  &
                         &   + zrodlo%ck(k)*deltamk*zrodlo%Chi_m_out(i,k)  &
                         &   - Xkn*       deltapk*zrodlo%Chi_m_in (i,k)
                     !print*,i,deltapk/2/DX
                enddo
                zrodlo%Fj(i) = zrodlo%Fj(i)*post
        enddo ! end of do i=1
        endif

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
        doubleprecision :: dx ,Ef ,Bz ,kvec
        integer         :: k,p
        complex*16      :: Xkn , deltamk , deltapk

        dx  = pdx
        Ef  = atomic_Ef/1000.0/Rd
        BZ  = BtoDonorB(atomic_Bz)

        if(zrodlo%bKierunek == .true.) then
        ! -----------------------------------------------------------------
        !       Wyznaczamy alpha(v,i) dla wejsc skierowanych w prawo
        !
        !                       --------------------->
        !
        ! -----------------------------------------------------------------
            Xkn  = 0
            do k = 1 , zrodlo%liczba_modow
                kvec = zrodlo%k_m_in(k)
                deltamk = 2*II*sin(-kvec*DX + DX*Bz*zrodlo%hny )
            do p = 1 , zrodlo%liczba_modow
                Xkn = Xkn + deltamk*zrodlo%Chi_m_out(v,k)*zrodlo%Sij(k,p)*conjg(zrodlo%Chi_m_out(i,p))
            enddo
            enddo
            zrodlo_alfa_v_i = Xkn*dx

        else
        ! -----------------------------------------------------------------
        !       Wyznaczamy alpha(v,i) dla wejsc skierowanych w lewo
        !
        !                       <---------------------
        !
        ! -----------------------------------------------------------------
            Xkn  = 0
            do k = 1 , zrodlo%liczba_modow
                kvec = zrodlo%k_m_in(k)
                deltapk = 2*II*sin(kvec*DX + DX*Bz*zrodlo%hny )
            do p = 1 , zrodlo%liczba_modow
                Xkn = Xkn + deltapk*zrodlo%Chi_m_in(v,k)*zrodlo%Sij(k,p)*conjg(zrodlo%Chi_m_in(i,p))
            enddo
            enddo
            zrodlo_alfa_v_i = Xkn*dx

        endif
    end function zrodlo_alfa_v_i


    subroutine zrodlo_oblicz_dk(zrodlo,VPHI,GINDEX,dx)
        class(czrodlo)          :: zrodlo
        complex*16,dimension(:) :: VPHI
        integer,dimension(:,:)  :: GINDEX
        double precision        :: dx

        integer :: p,q,k,i,pi,pj
        complex*16 :: dk , alpha_p, beta_p

        if(zrodlo%bKierunek == .true.) then
        ! ------------------------------------------------------------------
        !
        !                          <----------------
        !
        ! ------------------------------------------------------------------
        do k = 1 , zrodlo%liczba_modow
            dk = 0
            do p = 1 , zrodlo%liczba_modow
                alpha_p = 0
                do i = 1 , zrodlo%N
                    pi   = zrodlo%polozenia(i,1)
                    pj   = zrodlo%polozenia(i,2)
                    alpha_p = alpha_p + dx*conjg(zrodlo%Chi_m_out(i,p))*VPHI(GINDEX(pi,pj))
                enddo
                beta_p = 0
                do q = 1 , zrodlo%liczba_modow
                    beta_p = beta_p + zrodlo%Aij(p,q)*zrodlo%ck(q)
                enddo
                dk = dk + zrodlo%Sij(k,p)*(  alpha_p - beta_p )
            enddo ! end do p
            zrodlo%dk(k) = dk
        enddo ! end do k
        else ! else of if (kierunek)
        ! ------------------------------------------------------------------
        !
        !                          ---------------->
        !
        ! ------------------------------------------------------------------
            do k = 1 , zrodlo%liczba_modow
                dk = 0
                do p = 1 , zrodlo%liczba_modow
                    alpha_p = 0
                    do i = 1 , zrodlo%N
                        pi   = zrodlo%polozenia(i,1)
                        pj   = zrodlo%polozenia(i,2)
                        alpha_p = alpha_p + dx*conjg(zrodlo%Chi_m_in(i,p))*VPHI(GINDEX(pi,pj))
                    enddo
                    beta_p = 0
                    do q = 1 , zrodlo%liczba_modow
                        beta_p = beta_p + zrodlo%Aij(p,q)*zrodlo%ck(q)
                    enddo
                    dk = dk + zrodlo%Sij(k,p)*(  alpha_p - beta_p )
                enddo ! end do p
                zrodlo%dk(k) = dk
            enddo ! end do k

        endif ! end of if(kierunek)



    end subroutine zrodlo_oblicz_dk


    ! wzory brane z mojej pracy inzynierskiej (strona 12 - tabela)
    subroutine zrodlo_oblicz_JinJout(zrodlo,dx)
        class(czrodlo)    :: zrodlo
        double precision  :: dx

        integer          :: k , pj , i
        double precision :: Jin , Jout , Bz ,kvec

        BZ  = BtoDonorB(atomic_Bz)

        if(zrodlo%bKierunek == .true. ) then
        ! ------------------------------------------------------------------
        !
        !                          ---------------->
        !
        ! ------------------------------------------------------------------
            do k = 1 , zrodlo%liczba_modow
                Jin  = 0
                Jout = 0
                kvec = zrodlo%k_m_in(k)
                do i = 1 , zrodlo%N
                    pj   = zrodlo%polozenia(i,2)
                    Jin   = Jin  + abs(zrodlo%Chi_m_in(i,k))**2*sin( -dx*dx*pj*Bz + dx*(+kvec + zrodlo%hnY*Bz) )
                    Jout  = Jout + abs(zrodlo%Chi_m_out (i,k))**2*sin( -dx*dx*pj*Bz + dx*(-kvec + zrodlo%hnY*Bz) )

                enddo ! end of i

                Jin   = - Jin  * dx * abs(zrodlo%ck(k))**2
                Jout  = - Jout * dx * abs(zrodlo%dk(k))**2
                zrodlo%Jin(k)  = Jin
                zrodlo%Jout(k) = Jout
            enddo ! end of k

        else ! else of if kierunek
        ! ------------------------------------------------------------------
        !
        !                          <----------------
        !
        ! ------------------------------------------------------------------
            do k = 1 , zrodlo%liczba_modow
                Jin  = 0
                Jout = 0
                kvec = zrodlo%k_m_in(k)
                do i = 1 , zrodlo%N
                    pj   = zrodlo%polozenia(i,2)
                    Jout = Jout + abs(zrodlo%Chi_m_in(i,k))**2*sin( -dx*dx*pj*Bz + dx*(+kvec + zrodlo%hnY*Bz) )
                    Jin  = Jin  + abs(zrodlo%Chi_m_out (i,k))**2*sin( -dx*dx*pj*Bz + dx*(-kvec + zrodlo%hnY*Bz) )
                enddo ! end of i

                Jin   = - Jin  * dx * abs(zrodlo%ck(k))**2
                Jout  = - Jout * dx * abs(zrodlo%dk(k))**2
                zrodlo%Jin(k)  = Jin
                zrodlo%Jout(k) = Jout
            enddo ! end of k
        endif ! end of if kierunek

    end subroutine



endmodule modzrodlo
