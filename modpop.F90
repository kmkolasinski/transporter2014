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

        !Uvec        = pUVEC/Rd/1000
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
        integer :: i,j,num_in,num_out

        INTEGER                                      :: LDVL, LDVR , LWMAX , LWORK , INFO
        COMPLEX*16 , dimension(:) ,     allocatable  :: ALPHA , BETA , WORK
        double precision, dimension(:), allocatable  :: RWORK
        COMPLEX*16 , dimension(:,:)   , allocatable  :: Z
        COMPLEX*16 :: DUMMY(1,1),lambda , YcY
        complex*16 :: kvec
        bPoziome = .false. ! ustalamy na sztywno zrodla prawo - lewo

        ! konwersja jednostek do jednostek donorowych
        DX    = pDX*L2LR
        N     = pN-2
        Ef    = pEf/Rd/1000 ! bo Ef w meV
        t0    = 0.5/DX/DX
        BZ    = BtoDonorB(pB)
        hny   = (N+1)/2.0!


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
        do i = 1 , N
            Mdiag(i,i) =     1
            Mtau (i,i) =  - t0 * exp(II*(DX*DX*BZ*(i-hNY)))
        enddo
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
                !print"(i4,12f16.10)",i,abs(lambda),log(lambda)/II/DX*L2LR
                if(  abs(abs(lambda) - 1.0) < 1.0E-6 ) then
                    kvec  = log(lambda)/II/DX*L2LR
                    !print"(i4,12f16.10)",i,abs(lambda),log(lambda)/II/DX*L2LR
                    if(dble(kvec) > 0)LICZBA_MODOW = LICZBA_MODOW + 1
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
                    kvec  = (log(lambda)/II/DX)
                    if(dble(kvec) > 0) then
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

        call modpop_sort_vectors_by_values(Chi_M_IN (:,1:LICZBA_MODOW),K_M_IN (1:LICZBA_MODOW),+1)
        call modpop_sort_vectors_by_values(Chi_M_OUT(:,1:LICZBA_MODOW),K_M_OUT(1:LICZBA_MODOW),+1)


        print*,"-------------------------------------------------------------"
        print*," K wave.  :         Input mod.      |         Output mod."
        print*,"-------------------------------------------------------------"
        do i = 1 , LICZBA_MODOW
             print"(A,i4,A,2f12.6,A,2f12.6)","K[",i,"][nm]:",K_M_IN(i)*L2LR," | ",K_M_OUT(i)*L2LR
        enddo
        print*,"-------------------------------------------------------------"


        print*,""
        print*,"WEKTORY EVANESCENTNE (IN/OUT):"

        do i = N , 1 , -1
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

        do i = N , 2*N
            if(abs(Beta(i))>1e-16) then
                lambda= (ALPHA(i)/BETA(i))
                kvec  = (log(lambda)/DX)
                if(  abs(lambda) < 1 - 1E-6 ) then

                    num_in                   = num_in + 1
                    K_M_IN(num_in)           = kvec
                    Chi_M_IN(2:N+1,num_in)   = Z(N+1:2*N,i)
                    YcY                      = sum(abs(Chi_M_IN(:,num_in))**2)*DX
                    Chi_M_IN(:,num_in)       = Chi_M_IN(:,num_in)/sqrt(YcY)


                endif
            endif ! end of if beta
        enddo

        call modpop_sort_vectors_by_values(Chi_M_IN (:,LICZBA_MODOW+1:N),K_M_IN (LICZBA_MODOW+1:N),-1)
        call modpop_sort_vectors_by_values(Chi_M_OUT(:,LICZBA_MODOW+1:N),K_M_OUT(LICZBA_MODOW+1:N),-1)

        LICZBA_MODOW_EVANESCENTYCH = 0
        ! bierzemy tylko te mody ktore nie maja czesci falowej (tj tylko czyste evanescentne exp(+/-kx))
        do i = LICZBA_MODOW + 1 , N
            if( abs(imag(K_M_OUT(i))) < 1.0E-6 ) then
                LICZBA_MODOW_EVANESCENTYCH = LICZBA_MODOW_EVANESCENTYCH + 1
                Chi_M_IN (:,LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH) = Chi_M_IN(:,i)
                K_M_IN   (LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH)   = K_M_IN(i)
                Chi_M_OUT(:,LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH) = Chi_M_OUT(:,i)
                K_M_OUT  (LICZBA_MODOW + LICZBA_MODOW_EVANESCENTYCH)   = K_M_OUT(i)
            endif
        enddo

        print*,"-------------------------------------------------------------"
        print*," K evan.  :         Input mod.      |         Output mod."
        print*,"-------------------------------------------------------------"
        do i = LICZBA_MODOW + 1 , LICZBA_MODOW_EVANESCENTYCH + LICZBA_MODOW
             print"(A,i4,A,2f12.6,A,2f12.6)","K[",i,"][nm]:",K_M_IN(i)*L2LR," | ",K_M_OUT(i)*L2LR
        enddo
        print*,"-------------------------------------------------------------"

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


    subroutine modpop_sort_vectors_by_values(vectors,vals,order)
        complex*16,dimension(:,:) :: vectors
        complex*16,dimension(:)   :: vals
        integer :: order
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
                    if( abs(vals(imin)) < abs(vals(j)) ) imin = j
               endif
               if(order < 0 ) then
                    if( abs(vals(imin)) > abs(vals(j)) ) imin = j
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


END MODULE modpop



















