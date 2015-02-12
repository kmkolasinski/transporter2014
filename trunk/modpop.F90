!
! File:   modulModPoprzeczny.F90
! Author: mkk
!
! Created on 2 lipiec 2013, 15:17
!

MODULE modpop
    USE, INTRINSIC :: ISO_C_BINDING
    use            :: modjed
    implicit none
    private



    double precision :: DX  , Ef , t0 , BZ , hny , Emax , KIN
    integer          :: N , wypisz = 0
    logical          :: bPoziome, bEmax

    complex*16,dimension(:,:),allocatable        :: Hamiltonian
    double precision,dimension(:),allocatable    :: UVEC
    complex*16,dimension(:,:),allocatable        :: Chi_m_in ! mod wchodzacy do ukladu
    complex*16,dimension(:),allocatable          :: K_m_in
    complex*16,dimension(:,:),allocatable        :: Chi_m_out ! wychodzacy
    complex*16,dimension(:),allocatable          :: K_m_out
    integer                                      :: LICZBA_MODOW,L_M,LICZBA_MODOW_EVANESCENTYCH


    ! -------------------------------------------------
    !                    LAPACK
    ! -------------------------------------------------
    !     .. Local Scalars ..
      INTEGER          INFO, LWORK, LRWORK, LIWORK, IL, IU, M , LWMAX
      DOUBLE PRECISION ABSTOL, VL, VU
!
!     .. Local Arrays ..
      INTEGER,allocatable,dimension(:)          :: ISUPPZ, IWORK
      DOUBLE PRECISION,allocatable,dimension(:) :: W( : ), RWORK( : )
      COMPLEX*16,allocatable,dimension(:)       :: Z( :, : ), WORK( : )


      public :: modpop_inicjalizacja
      public :: modpop_znajdz_wartosci_wlasne
      public :: modpop_relacja_dyspersji
      public :: modpop_zwalnienie_pamieci
      public :: modpop_liczba_podpasm
      public :: modpop_zapisz_wektory
      public :: modpop_get_km
      public :: modpop_get_chi
      public :: modpop_calc_modes_from_wfm



    contains


    subroutine modpop_inicjalizacja(pDx,pN,pEf,pB,pUvec)
        double precision,intent(in)               :: pDX
        integer,intent(in)                        :: pN
        double precision,intent(in)               :: pEf
        double precision,intent(in)               :: pB
        double precision,intent(in),dimension(pN) :: pUVEC
        complex*16 :: compk
        integer    :: num
        bPoziome = .false. ! ustalamy na sztywno zrodla prawo - lewo

        DX    = pDX*L2LR
        N     = pN-2
        Ef    = pEf/Rd/1000 ! bo Ef w meV
        LWMAX = 50*N
        t0    = 0.5/DX/DX
        BZ    = BtoDonorB(pB)
        hny   = (N+1)/2.0!
        wypisz = 0
        if(TRANS_DEBUG)then
            print*,"Relacja dyspersji:"
            print*,"    N  :",N
            print*,"    hny:",hny
            print*,"t0 :",t0
            print*,"Rd :",Rd
            print*,"N  :",N
            print*,"Ef :",Ef
            print*,"DX :",DX
            print*,"hny:",hny
            print*,"B  :",pB
            print*,"Kmax:",3.14159/pDX
        endif

        allocate(Uvec(N))
        allocate(Hamiltonian(N,N))


        do num = 1 , N
            Uvec(num)        = pUVEC(num+1)/Rd/1000 ! bo Ef w meV
        enddo

        Hamiltonian = 0

        allocate(ISUPPZ( N ))
        allocate(IWORK( LWMAX ))
        allocate(W( N ))
        allocate(RWORK( LWMAX ))
        allocate(Z( N, N ))
        allocate(WORK( LWMAX ))

  !
  !     Query the optimal workspace.
  !


        LWORK  = -1
        LRWORK = -1
        LIWORK = -1

        ABSTOL = -1.0
        VL     =  0.0
        VU     =  Ef

        call modpop_utworz_hamiltonian((0.0D0,0.0D0))



        CALL ZHEEVR( "N", 'Values', 'Lower', N, Hamiltonian, N, VL, VU, IL,&
       &             IU, ABSTOL, M, W, Z, N, ISUPPZ, WORK, LWORK, RWORK,&
       &             LRWORK, IWORK, LIWORK, INFO )


        LWORK  = MIN( LWMAX, INT( WORK( 1 ) ) )
        LRWORK = MIN( LWMAX, INT( RWORK( 1 ) ) )
        LIWORK = MIN( LWMAX, IWORK( 1 ) )



        call modpop_znajdz_wartosci_wlasne("N",(0.0D0,0.0D0),Ef,.true.)
        bEmax = .true.
        compk = 3.14159/DX
        call modpop_znajdz_wartosci_wlasne("N",compk,5000*Ef,.false.)
        bEmax = .false.

        if(LICZBA_MODOW == 0) return;



!        allocate(Chi_m_in (N+2,LICZBA_MODOW))
!        allocate(K_m_in   (LICZBA_MODOW)  )
!        allocate(Chi_m_out(N+2,LICZBA_MODOW))
!        allocate(K_m_out  (LICZBA_MODOW)  )


!		call modpop_relacja_dyspersji(8,"rel.txt")

!        print*,"Wektory (IN):"
!        do num = 1, LICZBA_MODOW
!            call modpop_znajdz_k_in(num)
!        enddo
!        print*,"EWektory (IN):"
!        do num = LICZBA_MODOW + 1, LICZBA_MODOW
!            call modpop_znajdz_evan_k_in(num)
!        enddo
!
!        print*,"Wektory (OUT):"
!        do num = 1, LICZBA_MODOW
!            call modpop_znajdz_k_out(num)
!        enddo
!        print*,"EWektory (OUT):"
!        do num = LICZBA_MODOW + 1, LICZBA_MODOW
!            call modpop_znajdz_evan_k_out(num)
!        enddo
!        wypisz = 1



    end subroutine modpop_inicjalizacja


   subroutine modpop_relacja_dyspersji(dlugosc,nazwa_pliku)
       integer,intent(in) :: dlugosc
        character(LEN=dlugosc),intent(in) :: nazwa_pliku
        character(LEN=dlugosc+5):: evan_plik
        double precision :: K_VEC , GORNA_GRANICA
        complex*16        :: CKVEC
        integer           :: i

        wypisz = 0
        if(LICZBA_MODOW == 0) return;

        K_VEC         = -M_PI/DX !(K_m_OUT(1)+0.3*K_m_OUT(1))
        GORNA_GRANICA =  M_PI/DX !K_m_IN(1)+0.3*K_m_IN(1)

        print*,"Dlugosc:",    dlugosc
        print*,"Plik   :",    nazwa_pliku
        evan_plik = "evan_"//nazwa_pliku

        open(unit=1111,file=nazwa_pliku(1:dlugosc))
        open(unit=1112,file=evan_plik(1:dlugosc+5))

        do while( K_VEC <= GORNA_GRANICA )
        CKVEC = K_VEC
        call modpop_znajdz_wartosci_wlasne("N",CKVEC,Emax,.false.)

        if(L_M > 0 )then
            write(1111,"(2e16.8)",advance="no"),K_VEC*L2LR,Ef*Rd*1000
            do i = 1,L_M
                write(1111,"(e16.8)",advance="no"),W(i)*Rd*1000
            enddo
            write(1111,"(A)")," "
        endif

        CKVEC = K_VEC*II
        call modpop_znajdz_wartosci_wlasne("N",CKVEC,Emax,.false.)

        if(L_M > 0 )then
            write(1112,"(2e16.8)",advance="no"),K_VEC*L2LR,Ef*Rd*1000

            do i = 1,L_M
                write(1112,"(e16.8)",advance="no"),W(i)*Rd*1000
            enddo
            write(1112,"(A)")," "
        endif

        K_VEC = K_VEC + GORNA_GRANICA/200.0


        enddo
        close(1111)
        close(1112)


    end subroutine modpop_relacja_dyspersji


    subroutine modpop_zapisz_wektory(dlugosc,nazwa_pliku,p1,p2)
       integer,intent(in) :: dlugosc
        character(LEN=dlugosc),intent(in) :: nazwa_pliku
        double precision, intent(in) :: p1,p2

        character(LEN=(dlugosc+1)) :: wlasc_n_pliku
        integer :: i

        if(LICZBA_MODOW == 0) return;
        wlasc_n_pliku(1:dlugosc) = nazwa_pliku

        wlasc_n_pliku((dlugosc-3):dlugosc+1) = "+.txt"
        print*,"Zapis wektorow (+) do:",wlasc_n_pliku

        open(unit=1111,file=wlasc_n_pliku)
        !do i = 1 , N + 2
        do i = 1 , N + 2
            write(1111,"(2e20.6)",advance='no'),p1 + (i-1)*DX*LR2L,p2
            !write(1111,"(1000f20.10)"),abs(Chi_m_IN(i,:))**2!,DBLE(Chi_m_IN(i,m)),IMAG(Chi_m_IN(i,m))
            write(1111,"(1000f20.10)"),abs(Chi_m_IN(i,:))**2,DBLE(Chi_m_IN(i,:)),IMAG(Chi_m_IN(i,:))
        enddo
        close(1111)

        wlasc_n_pliku((dlugosc-3):dlugosc+1) = "-.txt"

        print*,"Zapis wektorow (-) do:",wlasc_n_pliku

        open(unit=1111,file=wlasc_n_pliku)
        !do i = 1 , N + 2
        do i = 1 , N + 2
            write(1111,"(2e20.6)",advance='no'),p1 + (i-1)*DX*LR2L,p2
            write(1111,"(1000f20.10)"),abs(Chi_m_OUT(i,:))**2,DBLE(Chi_m_OUT(i,:)),IMAG(Chi_m_OUT(i,:))
        enddo
        close(1111)


    end subroutine modpop_zapisz_wektory

    ! -----------------------------------------------------
    ! TYP: V - wszytko , N tylko wartosci wlasne
    ! -----------------------------------------------------
    subroutine modpop_znajdz_wartosci_wlasne(typ,pk,pEf,bZnajdzLiczbeModow)
        character(*),intent(in)      :: typ
        complex*16, intent(in) :: pk
        double precision, intent(in) :: pEf
        logical,optional,intent(in)  :: bZnajdzLiczbeModow


        call modpop_utworz_hamiltonian( pk )

        ABSTOL = -1.0
  !     Set VL, VU to compute eigenvalues in half-open (VL,VU] interval
        VL = -10000.0
        VU = pEf


  !
  !     Solve eigenproblem.
  !
        CALL ZHEEVR( typ, 'Values', 'Lower', N, Hamiltonian, N, VL, VU, IL,&
       &             IU, ABSTOL, M, W, Z, N, ISUPPZ, WORK, LWORK, RWORK,&
       &             LRWORK, IWORK, LIWORK, INFO )

  !
  !     Check for convergence.
  !
        IF( INFO.GT.0 ) THEN
           WRITE(*,*)'Obliczanie modu poprzecznego zawiodlo na ',pk, " przy Ef=", Ef
           STOP
        END IF
  !
  !     Print the number of eigenvalues found.
  !

        if(present(bZnajdzLiczbeModow)) then
            if(bZnajdzLiczbeModow .eqv. .true. ) then
                LICZBA_MODOW = M
                WRITE(*,'(A,I2)')' Liczba modow:', LICZBA_MODOW
            endif
        endif
        L_M = M

        if(bEmax .eqv. .true.) then
	        Emax = W(1)*1000*Rd
        endif

    end subroutine modpop_znajdz_wartosci_wlasne

    ! --------------------------------------------
    ! Tworzenie hamiltonianu
    ! --------------------------------------------
    subroutine modpop_utworz_hamiltonian(pk)
        complex*16, intent(in) :: pk

            integer :: i
            complex*16 :: alpha



        if(bPoziome) then
            alpha              = exp( II*(pk*DX + DX*DX*BZ*(1-hNY)) ) + exp(-II*(pk*DX + DX*DX*BZ*(1-hNY)) )
            Hamiltonian(1,2)   = -t0
            Hamiltonian(1,1)   =  UVEC(1) + 4*t0 - t0*alpha
        do i = 2 , N - 1
            alpha              = exp( II*(pk*DX + DX*DX*BZ*(i-hNY)) ) + exp(-II*(pk*DX + DX*DX*BZ*(i-hNY)) )
            Hamiltonian(i,i-1) = -t0
            Hamiltonian(i,i+1) = -t0
            Hamiltonian(i,i)   =  UVEC(i) + 4*t0 - t0*alpha
        enddo
            alpha              = exp( II*(pk*DX + DX*DX*BZ*(N-hNY)) ) + exp(-II*(pk*DX + DX*DX*BZ*(N-hNY)) )
            Hamiltonian(N,N-1) = -t0
            Hamiltonian(N,N)   =  UVEC(N) + 4*t0 - t0*alpha
        else
            alpha              = exp( II*(pk*DX - DX*DX*BZ*(1-hNY)) ) + exp(-II*(pk*DX - DX*DX*BZ*(1-hNY)) )
            Hamiltonian(1,2)   = -t0
            Hamiltonian(1,1)   =  UVEC(1) + 4*t0 - t0*alpha
        do i = 2 , N - 1
            alpha              = exp( II*(pk*DX - DX*DX*BZ*(i-hNY)) ) + exp(-II*(pk*DX - DX*DX*BZ*(i-hNY)) )
            Hamiltonian(i,i-1) = -t0
            Hamiltonian(i,i+1) = -t0
            Hamiltonian(i,i)   =  UVEC(i) + 4*t0 - t0*alpha
        enddo
            alpha              = exp( II*(pk*DX - DX*DX*BZ*(N-hNY)) ) + exp(-II*(pk*DX - DX*DX*BZ*(N-hNY)) )
            Hamiltonian(N,N-1) = -t0
            Hamiltonian(N,N)   =  UVEC(N) + 4*t0 - t0*alpha
        endif

       if(wypisz ==  1)      then
            print*,1,1,Hamiltonian(1,1)
            print*,1,2,Hamiltonian(1,2)
        do i = 2 , N - 1
            print*,i,i-1,Hamiltonian(i,i-1)
            print*,i,i,Hamiltonian(i,i)
            print*,i,i+1,Hamiltonian(i,i+1)
        enddo
            print*,N,N-1,Hamiltonian(N,N-1)
            print*,N,N,Hamiltonian(N,N)
       endif

    endsubroutine modpop_utworz_hamiltonian


    subroutine modpop_zwalnienie_pamieci()


        print*,"MODPOP:Zwalnianie pamieci."

        if(allocated(UVEC))        deallocate(Uvec)
        if(allocated(Hamiltonian)) deallocate(Hamiltonian)


        if(allocated(ISUPPZ))deallocate(ISUPPZ)
        if(allocated(IWORK ))deallocate(IWORK )
        if(allocated(W     ))deallocate(W     )
        if(allocated(RWORK ))deallocate(RWORK )
        if(allocated(Z     ))deallocate(Z     )
        if(allocated(WORK  ))deallocate(WORK  )

        if(allocated(Chi_m_IN))  deallocate(Chi_m_IN)
        if(allocated(K_m_IN))    deallocate(K_m_IN)
        if(allocated(Chi_m_OUT)) deallocate(Chi_m_OUT)
        if(allocated(K_m_OUT))   deallocate(K_m_OUT)

    end subroutine modpop_zwalnienie_pamieci



    subroutine modpop_liczba_podpasm(no_podpasm,no_evan)
        integer , intent(inout) :: no_podpasm,no_evan
        no_podpasm = LICZBA_MODOW
        no_evan    = LICZBA_MODOW_EVANESCENTYCH

    end subroutine modpop_liczba_podpasm

    subroutine modpop_get_km(pLiczbaModow,pKm_in,pKm_out)
        integer, intent(in)                              :: pLiczbaModow
        complex*16,dimension(pLiczbaModow),intent(inout) :: pKm_in,pKm_out

        pKm_in (1:pLiczbaModow) = K_m_in (1:pLiczbaModow)
        pKm_out(1:pLiczbaModow) = K_m_out(1:pLiczbaModow)

    end subroutine modpop_get_km


    subroutine modpop_get_chi(pLiczbaModow,pN,pchi_m_in,pchi_m_out)
        integer, intent(in)                                 :: pLiczbaModow,pN
        complex*16,dimension(pN,pLiczbaModow),intent(inout) :: pchi_m_in,pchi_m_out

        pchi_m_in (1:pN,1:pLiczbaModow) = Chi_m_in (1:pN,1:pLiczbaModow)
        pchi_m_out(1:pN,1:pLiczbaModow) = Chi_m_out(1:pN,1:pLiczbaModow)

    end subroutine modpop_get_chi


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
    subroutine modpop_calc_modes_from_wfm(pDx,pN,pEf,pB,pUvec,pbHorizontal)
        double precision,intent(in)               :: pDX
        integer,intent(in)                        :: pN
        double precision,intent(in)               :: pEf
        double precision,intent(in)               :: pB
        double precision,intent(in),dimension(pN) :: pUVEC
        logical, intent(in) :: pbHorizontal

        ! zmienne pomocnicze
        complex*16, allocatable , dimension(:,:) :: Mdiag,Mham,MB,MA , Mtau
        integer :: i,num_in,num_out

        INTEGER                                      :: LDVL, LDVR , LWMAX , LWORK , INFO
        COMPLEX*16 , dimension(:) ,     allocatable  :: ALPHA , BETA , WORK
        double precision, dimension(:), allocatable  :: RWORK
        COMPLEX*16 , dimension(:,:)   , allocatable  :: Z
        COMPLEX*16 :: DUMMY(1,1),lambda , YcY
        complex*16 :: kvec
        doubleprecision :: dkvec
        bPoziome = .false. ! ustalamy na sztywno zrodla prawo - lewo

        ! konwersja jednostek do jednostek donorowych
        DX    = pDX*L2LR
        N     = pN-2
        Ef    = pEf/Rd/1000 ! bo Ef w meV
        t0    = 0.5/DX/DX
        BZ    = BtoDonorB(pB)
        hny   = (N+1)/2.0! dobor odpowieniego cechowania



        ! alokacja tablic
        if(.not. allocated(Uvec)) allocate(Uvec(N))

        do i = 1 , N
            Uvec(i)        = pUVEC(i+1)/Rd/1000 ! bo Ef w meV
        enddo


        allocate(Mtau(N,N))
        allocate(Mdiag(N,N))
        allocate(Mham (N,N))
        allocate(MA (2*N,2*N))
        allocate(MB (2*N,2*N))
        allocate(Z  (2*N,2*N))


        ! Tworzenie rownania wlasnego  - wzor (52)
        Mdiag = 0
        Mham  = 0
        Mtau  = 0

        ! Przygotowanie podmacierzy diagonalnej 1 i macierzy B

        if(pbHorizontal) then
            do i = 1 , N
                Mdiag(i,i) =     1
                Mtau (i,i) =  - t0 * exp(+II*(DX*DX*BZ*(i-hNY)))
            enddo
        else
            do i = 1 , N
                Mdiag(i,i) =     1
                Mtau (i,i) =  - t0 * exp(-II*(DX*DX*BZ*(i-hNY)))
            enddo
        endif
        ! Przygotowanei podmacierzy (E-H)
        do i = 1 , N
            Mham(i,i)   = 4*t0 + Uvec(i) - Ef
            if(i < N) Mham(i,i+1) = -t0
            if(i > 1) Mham(i,i-1) = -t0
        enddo

        ! Wypelnienie macierzy:
        MA  = 0
        MB  = 0

        MA(N+1:2*N,1:N)     =  Mdiag
        MA(1:N,N+1:2*N)     =  Mtau
        MA(N+1:2*N,N+1:2*N) =  Mham

        MB(1:N,1:N)         =  Mdiag
        MB(N+1:2*N,N+1:2*N) = -conjg(Mtau)


        ! Ustalenie parametrow LAPACKA
        LWMAX = 50 * N
        LDVL  = 2  * N
        LDVR  = 2  * N
        ! Alokacja macierzy LAPACKA
        allocate(ALPHA(2*N))
        allocate(BETA(2*N))
        allocate(RWORK(8*N))
        allocate(WORK(LWMAX))



        LWORK = -1

!        przygotowanie
        CALL ZGGEV("N","N", 2*N, MA, 2*N, MB,2*N, ALPHA,BETA, &
    &   DUMMY, 1, Z, LDVR, WORK, LWORK, RWORK, INFO )

        LWORK = MIN(LWMAX, INT( WORK(1)))


        ! rozwiazanie problemu wlasnego
        CALL ZGGEV("N","V", 2*N, MA, 2*N, MB , 2*N, ALPHA,BETA, &
    &   DUMMY, 1, Z, LDVR, WORK, LWORK, RWORK, INFO )

        ! sprawdzamy czy uklad zostal skontruowany poprawnie
        if( INFO /= 0 ) then
            print*,"Modpop::WFM::Nie udalo sie poprawnie znalezc modow"
            stop
        endif


        ! Obliczanie liczby modow ,
        ! pamietamy ze: lambda = exp(i*k*DX)
        LICZBA_MODOW = 0
        do i = 1 , 2*N
            if(abs(Beta(i))>1e-16) then ! zgodnie z przykladem LApacka
                lambda= (ALPHA(i)/BETA(i))

                if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then ! jesli nie mod evanescentny
                    kvec  = log(lambda)/II/DX   ! spradzamy znak Kvec zeby policzyc ile jest wektorow
                                                ! w jedna strone

                    if(dble(kvec) > 0) LICZBA_MODOW = LICZBA_MODOW + 1
                endif
            endif
        enddo



        print*,"Liczba modow:",LICZBA_MODOW
        if(LICZBA_MODOW == 0) return;

        if(allocated(Chi_m_in))  deallocate(Chi_m_in)
        if(allocated(K_m_in))    deallocate(K_m_in)
        if(allocated(Chi_m_out)) deallocate(Chi_m_out)
        if(allocated(K_m_out))   deallocate(K_m_out)


        allocate(Chi_m_in (N+2,N))
        allocate(K_m_in   (N)   )
        allocate(Chi_m_out(N+2,N))
        allocate(K_m_out  (N)   )

        Chi_m_in  = 0
        Chi_m_out = 0
        print*,""
        print*,"WEKTORY FALOWE (IN/OUT):"


        num_in  = 0
        num_out = 0
        do i = 1 , 2*N
            if(abs(Beta(i))>1e-16) then
                lambda = (ALPHA(i)/BETA(i))

                if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then
                    kvec  = (log(lambda)/DX)
                    !kvec  = (log(lambda)/II/DX) ! stare
                    if(imag(kvec) > 0) then
                    !if(dble(kvec) > 0) then
                        num_in                   = num_in + 1
                        K_M_IN(num_in)           = kvec
                        Chi_M_IN(2:N+1,num_in)   = Z(N+1:2*N,i)
                        YcY                      = sum(abs(Chi_M_IN(:,num_in))**2)*DX
                        Chi_M_IN(:,num_in)       = Chi_M_IN(:,num_in)/sqrt(YcY)
                        !print"(A,i4,A,2f12.6,A)","  K_IN (",num_in,")=",K_M_IN(num_in)*L2LR,"[nm]"
                    else
                        num_out                  = num_out + 1
                        K_m_out(num_out)         = kvec
                        Chi_m_out(2:N+1,num_out) = Z(N+1:2*N,i)
                        YcY                      = sum(abs(Chi_m_out(:,num_out))**2)*DX
                        Chi_m_out(:,num_out)     = Chi_m_out(:,num_out)/sqrt(YcY)
                        !print"(A,i4,A,2f12.6,A)","  K_OUT(",num_out,")=",K_m_out(num_out)*L2LR,"[nm]"
                    endif
                endif
            endif
        enddo

        call modpop_sort_vectors_by_values(Chi_M_IN (:,1:LICZBA_MODOW),K_M_IN (1:LICZBA_MODOW),+1,0)
        call modpop_sort_vectors_by_values(Chi_M_OUT(:,1:LICZBA_MODOW),K_M_OUT(1:LICZBA_MODOW),+1,0)

!        if( TRANS_DEBUG ) then
!        print*,"-------------------------------------------------------------"
!        print*," K wave.  :         Input mod.      |         Output mod."
!        print*,"-------------------------------------------------------------"
!        do i = 1 , LICZBA_MODOW
!             print"(A,i4,A,2f12.6,A,2f12.6)","K[",i,"][nm]:",K_M_IN(i)*L2LR," | ",K_M_OUT(i)*L2LR
!        enddo
!        print*,"-------------------------------------------------------------"
!
!
!        print*,""
!        print*,"WEKTORY EVANESCENTNE (IN/OUT):"
!
!        endif ! end of if(TRANS_DEBUG)

        do i = 2*N , 1 , -1
            if(abs(Beta(i))>1e-16) then
                lambda= (ALPHA(i)/BETA(i))
                kvec  = (log(lambda)/DX)

                if( abs(lambda) > 1 + 1E-6 ) then
                    num_out                  = num_out + 1
                    K_m_out(num_out)         = kvec
                    Chi_m_out(2:N+1,num_out) = Z(N+1:2*N,i)
                    YcY                      = sum(abs(Chi_m_out(:,num_out))**2)*DX
                    Chi_m_out(:,num_out)     = Chi_m_out(:,num_out)/sqrt(YcY)
                endif
            endif ! end of if beta
        enddo

        do i = 1 , 2*N
            if(abs(Beta(i))>1e-16) then
                lambda= (ALPHA(i)/BETA(i))
                kvec  = (log(lambda)/DX)
!                print*,i,kvec*L2LR
                if(  abs(lambda) < 1 - 1E-6 ) then

                    num_in                   = num_in + 1
                    K_M_IN(num_in)           = kvec
                    Chi_M_IN(2:N+1,num_in)   = Z(N+1:2*N,i)
                    YcY                      = sum(abs(Chi_M_IN(:,num_in))**2)*DX
                    Chi_M_IN(:,num_in)       = Chi_M_IN(:,num_in)/sqrt(YcY)
                endif
            endif ! end of if beta
        enddo

        call modpop_sort_vectors_by_values(Chi_M_IN (:,LICZBA_MODOW+1:N),K_M_IN (LICZBA_MODOW+1:N),-1,1)
        call modpop_sort_vectors_by_values(Chi_M_OUT(:,LICZBA_MODOW+1:N),K_M_OUT(LICZBA_MODOW+1:N),-1,1)

        LICZBA_MODOW_EVANESCENTYCH = 0
        ! bierzemy tylko te mody ktore nie maja czesci falowej (tj tylko czyste evanescentne exp(+/-kx))
        do i = LICZBA_MODOW + 1 , N
            !if( abs(imag(K_M_OUT(i))) < 1.0E-6 .and. abs(imag(K_M_IN(i))) < 1.0E-6 .and. &
            ! &  abs(dble(K_M_OUT(i))) > 1.0E-6 .and. abs(dble(K_M_IN(i))) > 1.0E-6  ) then
                LICZBA_MODOW_EVANESCENTYCH = LICZBA_MODOW_EVANESCENTYCH + 1
!                Chi_M_IN (:,LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH) = Chi_M_IN(:,i)
!                K_M_IN   (LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH)   = K_M_IN(i)
!!
!                Chi_M_OUT(:,LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH) = Chi_M_OUT(:,i)
!                K_M_OUT  (LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH)   = K_M_OUT(i)

                Chi_M_OUT(:,LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH) = Chi_M_IN(:,i)
                K_M_OUT  (LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH)   = -K_M_IN(i)
                call modpop_calc_mode_from_k((K_M_OUT(i)),Chi_m_out(:,i),DX,N+2,Ef,BZ,Uvec,pbHorizontal);
           !endif

        enddo
        LICZBA_MODOW_EVANESCENTYCH = 1*LICZBA_MODOW_EVANESCENTYCH/8

        !i = 3
        !call modpop_calc_mode_from_k((K_M_IN(i)),Chi_m_out(:,i),DX,N+2,Ef,BZ,Uvec,pbHorizontal);

!        call zrodlo_mode_from_kvec(K_M_IN(i),N+2,UVEC,Chi_M_IN(:,i))
!        open(unit = 333, file= "evan.txt" )
!        do dkvec = -3.14159/DX/10 , 3.14159/DX/10 , 0.01
!            K_M_IN(i) = - II * CMPLX(0.0D0,dkvec)
!            call zrodlo_mode_from_kvec(K_M_IN(i),N+2,UVEC,Chi_M_IN(:,i))
!        enddo


!        if ( TRANS_DEBUG ) then
!        print*,"-------------------------------------------------------------"
!        print*," K evan.  :         Input mod.      |         Output mod."
!        print*,"-------------------------------------------------------------"
!        do i = LICZBA_MODOW + 1 , LICZBA_MODOW_EVANESCENTYCH + LICZBA_MODOW
!             print"(A,i4,A,2f12.6,A,2f12.6)","K[",i,"][nm]:",K_M_IN(i)*L2LR," | ",K_M_OUT(i)*L2LR
!        enddo
!        print*,"-------------------------------------------------------------"
!        endif


        deallocate(ALPHA)
        deallocate(BETA)
        deallocate(RWORK)
        deallocate(WORK)
        deallocate(Z)


        deallocate(Mdiag)
        deallocate(Mham )
        deallocate(MA   )
        deallocate(MB   )
        deallocate(Uvec )

    end subroutine modpop_calc_modes_from_wfm


    subroutine modpop_sort_vectors_by_values(vectors,vals,order,takeReal)
        complex*16,dimension(:,:) :: vectors
        complex*16,dimension(:)   :: vals
        integer :: order,takeReal
        complex*16,dimension(:),allocatable  :: tmpvec
        complex*16 :: tmpval

        integer :: i , j , n , nvec , imin


        n    = size(vals)
        nvec = size(vectors(:,1))


        allocate(tmpvec(nvec))
        if(TRANS_DEBUG) then
         print*,"Sortowanie wektorow:"
         print*,"  Liczba wektor  :"  ,n
         print*,"  Rozmiar wektora:",nvec
        endif


        do i = 1 , n
            imin   = i
            do j = i+1 , n
               if(abs(vals(j)) > 1.0e-10) then
               if(order > 0 ) then
                    if(takeReal == 0) then
                        if( abs((vals(imin))) < abs((vals(j))) ) imin = j
                    else
                        if( abs(real(vals(imin))) < abs(real(vals(j))) ) imin = j
                    endif
               endif
               if(order < 0 ) then
                    if(takeReal == 0) then
                        if( abs((vals(imin))) > abs((vals(j))) ) imin = j
                    else
                        if( abs(real(vals(imin))) > abs(real(vals(j))) ) imin = j
                    endif
               endif
               endif
            enddo

            tmpval     = vals(i)
            vals(i)    = vals(imin)
            vals(imin) = tmpval

            tmpvec          = vectors(:,i)
            vectors(:,i)    = vectors(:,imin)
            vectors(:,imin) = tmpvec
        enddo

        deallocate(tmpvec)
    end subroutine modpop_sort_vectors_by_values

subroutine modpop_calc_mode_from_k(kvec,modek,pDx,pN,pEf,pB,pUvec,pbHorizontal)
complex*16,intent(in)               :: kvec
complex*16,intent(inout),dimension(:)     :: modek
double precision,intent(in)               :: pDX
integer,intent(in)                        :: pN
double precision,intent(in)               :: pEf
double precision,intent(in)               :: pB
double precision,intent(in),dimension(pN) :: pUVEC
logical, intent(in) :: pbHorizontal
doubleprecision :: t0 , YcY ,delta
complex*16 :: Cstar , Cnostar
integer :: i,j

complex*16,dimension(pN) :: a , c , b , r , Upet , u

Upet = 0
Upet(2:pN-1) = pUVEC(:)

if(pbHorizontal == .true.) then
!    do i = 1 , pN
!    write(111,"(4f20.8)"),i*1.0D0,modek(i),abs(modek(i))**2
!    enddo

    t0 = 0.5/pdx/pdx

!    b= 0
!
!    do i = 2 , pN-1
!        b(i) = sin(i*3.1)
!    enddo
!    r            =  0
!    !r(3) = 1
!    YcY = sum(abs(b(:))**2)*pDX
!    b(:) = b(:)/sqrt(YcY)
!    u = 0
!    delta = 1
!    do while( delta > 1.0E-16 )
!
!    delta = abs(sum(abs(u-b)**2)/sum(abs(u+b)**2))
!    u = b
!
!    do i = 2 , pN-1
!        Cstar = t0*exp(kvec*pdx)!*exp(-II*pDX*pDX*pB*(i-(pN+1)/2.0))
!        b(i) = t0*(b(i-1)+b(i+1))/( 4*t0 + Upet(i)*0 - pEf - Cstar - conjg(Cstar) )
!    enddo
!    YcY = sum(abs(b(:))**2)*pDX
!    b(:) = b(:)/sqrt(YcY)
!    enddo
!    modek = u



    !call  tridag(a,b,c,r,u,pN)
    !modek  = u


!    modek = 0
!    modek(1) = 0
!    modek(2) = 1.0
!    do i = 3 , pN-1
!        Cstar   = exp(+kvec*pdx)*exp(-II*pDX*pDX*pB*(i-1-(pN+1)/2.0))
!        Cnostar = exp(-kvec*pdx)*exp(+II*pDX*pDX*pB*(i-1-(pN+1)/2.0))
!        modek(i) = -modek(i-2) + modek(i-1)*( 4 + UVEC(i-1)/t0 - pEf/t0 - Cstar - Cnostar )
!    enddo
!    modek(pN)  = 0
!    modek(pN-1) = 1.0
!    do i = pN-2 , 2 , -1
!        Cstar = exp(+kvec*pdx)*exp(-II*pDX*pDX*pB*(i+1-(pN+1)/2.0))
!        modek(i) = -modek(i+2) + modek(i+1)*( 4 + UVEC(i+1)/t0 - pEf/t0 - Cstar - conjg(Cstar) )
!    enddo

    a = modek
    do i = 1 , pN
        modek(i) = a(pN+1-i)
    enddo

    YcY = sum(abs(modek(:))**2)*pDX
    modek(:) = modek(:)/sqrt(YcY)

!    do i = 1 , pN
!    write(112,"(4f20.8)"),i*1.0D0,modek(i),abs(modek(i))**2
!    enddo

else
    print*,"ERROR::modpop_calc_mode_from_k nie ma jeszcze wejsc gora/dol"
    stop
endif
!
!
end subroutine modpop_calc_mode_from_k

!
!      SUBROUTINE tridag(a,b,c,r,u,n)
!      INTEGER n,NMAX
!      complex*16 a(n),b(n),c(n),r(n),u(n)
!      PARAMETER (NMAX=50000)
!      INTEGER j
!      complex*16 bet,gam(NMAX)
!      if(b(1).eq.0.)pause 'tridag: rewrite equations'
!      bet=b(1)
!      u(1)=r(1)/bet
!      do 11 j=2,n
!        gam(j)=c(j-1)/bet
!        bet=b(j)-a(j)*gam(j)
!        if(bet.eq.0.)pause 'tridag failed'
!        u(j)=(r(j)-a(j)*u(j-1))/bet
!11    continue
!      do 12 j=n-1,1,-1
!        u(j)=u(j)-gam(j+1)*u(j+1)
!12    continue
!      return
!      ENDSUBROUTINE tridag

! dopuszczalne wartosci to: {3,4,5,6,8,10,12,16,20,24,32,40,48}
  subroutine zrodlo_mode_from_kvec(kvec,pN,pUVEC,modek,pliczba_konturow,pwypisz_informacje,pmaks_iter)
        complex*16 :: kvec
        integer :: pN
        double precision,intent(in),dimension(pN) :: pUVEC
        complex*16,intent(inout),dimension(:)     :: modek
        integer           :: NoStates
        integer,optional  :: pliczba_konturow,pwypisz_informacje,pmaks_iter


        integer :: fpm(128),MATASIZE,TRANS_MAXN
        integer :: i,j,info,itmp,nw,M0,loop,no_evals,iter,Widmo_NoStates
        integer :: liczba_konturow,wypisz_informacje,maks_iter
        doubleprecision :: epsout
        doubleprecision :: Emin, Emax , Ecurr , YcY


        integer,allocatable                    :: HBROWS(:),IDXA(:,:)
        complex*16,dimension(:,:), allocatable :: EVectors
        complex*16,dimension(:), allocatable   :: CMatA
        double precision,dimension(:), allocatable   :: Evalues,Rerrors

        complex*16 :: Cstar , Cnostar


!        call reset_clock()
!        ! Przejscie do jednostek donorowych

        Emin = Ef - Ef / 2.0
        Emax = Ef + Ef / 2.0
        NoStates   = pN/2
        TRANS_MAXN = pN-2
        print*,kvec*L2LR

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
        ! ustalanie domyslnej liczby konturow i wypisywania
        if(.not. present(pliczba_konturow)) then
            liczba_konturow = 16
        else
            liczba_konturow = pliczba_konturow
        endif
        if(.not. present(pwypisz_informacje)) then
            wypisz_informacje = 0
        else
            wypisz_informacje = pwypisz_informacje
        endif
        if(.not. present(pmaks_iter)) then
            maks_iter = 200
        else
            maks_iter = pmaks_iter
        endif

        call feastinit(fpm)

        fpm(1)=wypisz_informacje ! nie wypisuj informacji
        fpm(2)=liczba_konturow   ! liczba konturow
        fpm(3)=12                ! wykladnik bledu ponizej ktorego procedura sie zatrzymuje: e=10^(fpm(3))
        fpm(4)=maks_iter         ! maksymalna liczba iteracji po ktorej jak sie nie zbiegnie to proced. sie zatrzyma
        fpm(5)=0                 ! startujemy z domyslnymi wektorami (jak 1 to z dostarczonymi)
        fpm(6)=0                 ! kryterium zbieznosci poprzez residuum (0 albo 1)
!
        allocate(CMATA(TRANS_MAXN*3))
        allocate(IDXA (TRANS_MAXN*3,2))

        itmp = 1
        do i = 2 , pN - 1
            iter = i - 1
            if(i > 2) then
                cmatA(itmp) = -t0
                idxA(itmp,1) = iter
                idxA(itmp,2) = iter-1
                itmp = itmp + 1
            endif
                Cstar   = exp(+kvec*DX)*exp(-II*DX*DX*BZ*(iter-hny))
                Cnostar = exp(-kvec*DX)*exp(+II*DX*DX*BZ*(iter-hny))
                cmatA(itmp) = 4*t0 + pUVEC(iter) - t0*(Cstar + Cnostar)
                idxA(itmp,1) = iter
                idxA(itmp,2) = iter
                itmp = itmp + 1

            if(i < pN-1) then
                cmatA(itmp) = -t0
                idxA(itmp,1) = iter
                idxA(itmp,2) = iter+1
                itmp = itmp + 1
            endif
        enddo


        itmp        = itmp - 1
        MATASIZE    = itmp
        nw          = itmp
        if(wypisz_informacje==1) then
            print*,"--------------------------------------------------"
            print*,"Widmo:"
            print*,"--------------------------------------------------"
            print*,"Rozmiar problemu N:     ",TRANS_MAXN
            print*,"Zakres energii:         ",Emin*1000.0*Rd," do ",Emax*1000.0*Rd,"w meV"
            print*,"Liczba konturow:        ",liczba_konturow
            print*,"Maksymalna liczba iter. :",maks_iter
            print*,"Zalozona liczba stanow: ",NoStates
        endif

        ! -----------  ----system_inicjalizacja_ukladu---------------
        !
        ! -----------------------------------------------------------
        allocate(HBROWS(TRANS_MAXN+1))
        call modpop_convert_to_HB(MATASIZE,IDXA,HBROWS)
!

        ! zgadujemy liczbe stanow
        M0  = NoStates

        allocate(EVectors(TRANS_MAXN,M0))
        allocate(Evalues(M0))
        allocate(Rerrors(M0))



        call zfeast_hcsrev('F',&                ! - 'F' oznacza ze podawana jest pelna macierz
                              TRANS_MAXN,&    ! - rozmiar problemu (ile wezlow z flaga B_NORMAL)
                              CMATA(1:nw),&    ! - kolejne nie zerowe wartosci w macierzy H
                              HBROWS,&         ! - numeracja wierszy (rodzaj zapisu macierzy rzakidch)
                              idxA(1:nw,2),&   ! - indeksy kolumn odpowiadaja tablicy wartosci CMATA
                              fpm,&            ! - wektor z konfiguracja procedury
                              epsout,&         ! - Residuum wyjsciowe
                              loop, &          ! - Koncowa liczba iteracji
                              Emin,&           ! - Minimalna energia przeszukiwania
                              Emax,&           ! - Maksymalna energia
                              M0,&             ! - Spodziewana liczba modow w zakresie (Emin,Emax)
                              Evalues,&        ! - Wektor z otrzymanymi wartosciami wlasnymi
                              EVectors,&       ! - Macierz z wektorami (kolejne kolumny odpowiadaja kolejnym wartoscia z tablicy Evalues)
                              no_evals,&       ! - Liczba otrzymanych wartosci z przedziale (Emin,Emax)
                              Rerrors,&        ! - Wektor z bledami dla kolejnych wartosci wlasnych
                              info)            ! - Ewentualne informacje o bledach



        if(wypisz_informacje==1) then
            print*,"Eps wyjsciowy           :",  epsout
            print*,"Liczba iteracji         :",  loop
            print*,"Znaleziona l. stanow    :",  no_evals
            print*,"Info                    :",  info
            print*,"--------------------------------------------------"
        endif

        Widmo_NoStates = no_evals

        ! ----------------------------------------------------------------------------------
        ! Obsluga bledow:
        ! ----------------------------------------------------------------------------------
        selectcase(info)
        case( 202 )
            print*," Error : Problem with size of the system n (n≤0) "
            stop
        case( 201 )
            print*," Error : Problem with size of initial subspace m0 (m0≤0 or m0>n) "
            stop
        case( 200 )
            print*," Error : Problem with emin,emax (emin≥emax) "
            stop
        case(100:199)
            print"(A,I4,A)"," Error : Problem with ",info-100,"-th value of the input Extended Eigensolver parameter (fpm(i)). Only the parameters in use are checked. "
            Widmo_NoStates = 0
            stop
        case( 4 )
            print*," Warning : Successful return of only the computed subspace after call withfpm(14) = 1 "
            Widmo_NoStates = 0

        case( 3 )
            print*," Warning : Size of the subspace m0 is too small (m0<m) "
            Widmo_NoStates = 0

        case( 2 )
            print*," Warning : No Convergence (number of iteration loops >fpm(4))"
            Widmo_NoStates = 0
        case( 1 )
            print*," Warning : No eigenvalue found in the search interval. See remark below for further details. "
            Widmo_NoStates = 0
        case( 0 )
            print*,               "---------------------------------------------"
            print"(A,i12)",       "Widmo: Znaleziono stanow :",Widmo_NoStates
            print"(A,e12.4)",     "       Z bledem epsout   :",epsout
            print"(A,i12)",       "       W po liczbie iter.:",loop
            print*,               "---------------------------------------------"
        case( -1 )
            print*," Error : Internal error for allocation memory. "
            stop
        case( -2 )
            print*," Error : Internal error of the inner system solver. Possible reasons: not enough memory for inner linear system solver or inconsistent input. "
            stop
        case( -3 )
            print*," Error : Internal error of the reduced eigenvalue solver Possible cause: matrix B may not be positive definite. It can be checked with LAPACK routines, if necessary."
            stop
        case(-199:-100)
            print"(A,I4,A)"," Error : Problem with the ",-info-100,"-th argument of the Extended Eigensolver interface. "
            stop
        endselect

        do i = 1 , Widmo_NoStates
            print*,i, Evalues(i)*1000.0*Rd
            YcY = sum(abs(EVectors(:,i))**2)*DX
            EVectors(:,i) = EVectors(:,i)/sqrt(YcY)

        enddo

!        write(333,"(40f20.8)"),kvec*L2LR,Ef*1000.0*Rd,Evalues(1:Widmo_NoStates)*1000.0*Rd
        do i = 1 , pN-2
        write(112,"(40f20.8)"),i*2.0D0,abs(modek(i+1))**2,abs(EVectors(i,1:Widmo_NoStates))**2
        enddo



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
        deallocate(CMATA)
        deallocate(IDXA)
        deallocate(HBROWS)
        deallocate(EVectors)
        deallocate(Evalues)
        deallocate(Rerrors)


    end subroutine zrodlo_mode_from_kvec

    subroutine modpop_convert_to_HB(no_vals,rows_cols,out_rows)
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
      end subroutine modpop_convert_to_HB






END MODULE modpop



















