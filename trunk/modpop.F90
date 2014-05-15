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
    integer           :: NO_EVAN_MODES
    complex*16,dimension(:,:),allocatable        :: Hamiltonian
    double precision,dimension(:),allocatable    :: UVEC
    complex*16,dimension(:,:),allocatable        :: Chi_m_in ! mod wchodzacy do ukladu
    double precision,dimension(:),allocatable    :: K_m_in
    complex*16,dimension(:,:),allocatable        :: Chi_m_out ! wychodzacy
    double precision,dimension(:),allocatable    :: K_m_out
    integer                                      :: LICZBA_MODOW,L_M


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
      public :: NO_EVAN_MODES


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

        NO_EVAN_MODES = 0
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
        if(LICZBA_MODOW+NO_EVAN_MODES > N) NO_EVAN_MODES = N  - LICZBA_MODOW
        print*,"Liczba modow evanescentnych:",NO_EVAN_MODES

        allocate(Chi_m_in (N+2,LICZBA_MODOW+NO_EVAN_MODES))
        allocate(K_m_in   (LICZBA_MODOW+NO_EVAN_MODES)  )
        allocate(Chi_m_out(N+2,LICZBA_MODOW+NO_EVAN_MODES))
        allocate(K_m_out  (LICZBA_MODOW+NO_EVAN_MODES)  )


!		call modpop_relacja_dyspersji(8,"rel.txt")

        print*,"Wektory (IN):"
        do num = 1, LICZBA_MODOW
            call modpop_znajdz_k_in(num)
        enddo
        print*,"EWektory (IN):"
        do num = LICZBA_MODOW + 1, LICZBA_MODOW + NO_EVAN_MODES
            call modpop_znajdz_evan_k_in(num)
        enddo

        print*,"Wektory (OUT):"
        do num = 1, LICZBA_MODOW
            call modpop_znajdz_k_out(num)
        enddo
        print*,"EWektory (OUT):"
        do num = LICZBA_MODOW + 1, LICZBA_MODOW + NO_EVAN_MODES
            call modpop_znajdz_evan_k_out(num)
        enddo
        wypisz = 1



!        print*,1.0/dx/dx/2*Rd*1000.0
!        print*,dx*dx*BZ*LR2L
        !call modpop_zapisz_wektory(7,"mod.txt",0.0D0,1.0D0)


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
            write(1111,"(1000f20.10)"),abs(Chi_m_IN(i,:))**2!,DBLE(Chi_m_IN(i,m)),IMAG(Chi_m_IN(i,m))
        enddo
        close(1111)

        wlasc_n_pliku((dlugosc-3):dlugosc+1) = "-.txt"

        print*,"Zapis wektorow (-) do:",wlasc_n_pliku

        open(unit=1111,file=wlasc_n_pliku)
        !do i = 1 , N + 2
        do i = 1 , N + 2
            write(1111,"(2e20.6)",advance='no'),p1 + (i-1)*DX*LR2L,p2
            write(1111,"(1000f20.10)"),abs(Chi_m_OUT(i,:))**2!,DBLE(Chi_m_OUT(i,m)),IMAG(Chi_m_OUT(i,m))
        enddo
        close(1111)


    end subroutine modpop_zapisz_wektory

   subroutine modpop_znajdz_k_in(num)
        integer, intent(in) :: num
        double precision    ::  k1 , k2 , kdiff , NK


        double precision :: YcY
        complex*16        :: cnk

        k1    = 0
        k2    = 3.14159/DX
        kdiff = 1.0
        do while( kdiff > 1e-6  )
            nk    = (k1+k2)/2
            Emax = Ef*2
            cnk  = nk
            call modpop_znajdz_wartosci_wlasne("N",cnk,Emax,.false.)
            if(L_M < num ) then
                k2 = nk
            else
                if ( W(num) > Ef  ) then
                    k2 = nk
                else
                    k1 = nk
                endif
            endif
            kdiff = abs( (k2-k1)/k2)
        enddo


        cnk = nk
        call modpop_znajdz_wartosci_wlasne("V",CNK,Emax,.false.)
        Chi_M_IN(:,num)     = 0
        K_M_IN(num)         = NK



!        !Chi_M_IN(:,num) = Z(:,num)
        Chi_M_IN(2:N+1,num) = Z(1:N,num)

        YcY          = sum(abs(Chi_M_IN(:,num))**2)*DX
        Chi_M_IN(:,num) = Chi_M_IN(:,num)/sqrt(YcY)!/SQRT(2*3.14159*NK)     ! dzielimy do rachunku zaburzen

        print*,"K_IN(",num,")=",K_M_IN(num)*L2LR,"[nm]"



    end subroutine modpop_znajdz_k_in

    subroutine modpop_znajdz_evan_k_in(num)
        integer, intent(in) :: num
        double precision    ::  k1 , k2 , kdiff , NK


        double precision :: YcY
        complex*16        :: cnk

        k1    = 0
        k2    = 3.14159/DX
        kdiff = 1.0

        do while( kdiff > 1e-6  )

            nk   = (k1+k2)/2
            Emax =  Ef*2*NO_EVAN_MODES
            cnk  = -nk*II

            call modpop_znajdz_wartosci_wlasne("N",cnk,Emax,.false.)
            if ( W(num) < Ef  ) then
                k2 = nk
            else
                k1 = nk
            endif

            kdiff = abs( (k2-k1)/k2)
        enddo


        cnk = nk
        call modpop_znajdz_wartosci_wlasne("V",CNK,Emax,.false.)


!        Chi_M_OUT(:,num)     = 0
!        K_M_OUT(num)         = NK
!        Chi_M_OUT(2:N+1,num) = Z(1:N,num)
!        YcY              = sum(abs(Chi_M_OUT(:,num))**2)*DX
!        Chi_M_OUT(:,num) = Chi_M_OUT(:,num)/sqrt(YcY)!/SQRT(2*3.14159*K_M_IN(num))
!        print*,"EK_OUT(",num,")=",K_M_OUT(num)*L2LR,"[nm]"

        Chi_M_IN(:,num)     = 0
        K_M_IN(num)         = NK
        Chi_M_IN(2:N+1,num) = Z(1:N,num)
        YcY                 = sum(abs(Chi_M_IN(:,num))**2)*DX
        Chi_M_IN(:,num)     = Chi_M_IN(:,num)/sqrt(YcY)!/SQRT(2*3.14159*NK)     ! dzielimy do rachunku zaburzen
        print*,"EK_IN(",num,")=",K_M_IN(num)*L2LR,"[nm]"



    end subroutine modpop_znajdz_evan_k_in


    subroutine modpop_znajdz_k_out(num)
        integer, intent(in) :: num
        !double precision    :: K_VEC , k1 , k2 , kdiff , NK
        double precision    ::  NK

        complex*16        :: cnk
        double precision :: YcY

!        k1    = 0
!        k2    =-3.14159/DX
!        kdiff = 1.0
!
!        do while( kdiff > 1.0e-6  )
!            nk    = (k1+k2)/2
!            K_VEC = nk
!            cnk   = nk
!            call modpop_znajdz_wartosci_wlasne("N",cnk,Emax,.false.)
!            if(L_M < num ) then
!                k2 = nk
!            else
!                if ( W(num) > Ef  ) then
!                    k2 = nk
!                else
!                    k1 = nk
!                endif
!            endif
!            kdiff = abs( (k2-k1)/k2)
!        enddo

    	nk  = -K_M_IN(num)
        cnk = nk
        call modpop_znajdz_wartosci_wlasne("V",CNK,Emax,.false.)

        Chi_M_OUT(:,num)     = 0
        K_M_OUT(num)         = NK
        Chi_M_OUT(2:N+1,num) = Z(1:N,num)

        YcY              = sum(abs(Chi_M_OUT(:,num))**2)*DX
        Chi_M_OUT(:,num) = Chi_M_OUT(:,num)/sqrt(YcY)!/SQRT(2*3.14159*K_M_IN(num))

        print*,"K_OUT(",num,")=",K_M_OUT(num)*L2LR,"[nm]"



    end subroutine modpop_znajdz_k_out

 subroutine modpop_znajdz_evan_k_out(num)
        integer, intent(in) :: num
        !double precision    :: K_VEC , k1 , k2 , kdiff , NK
        double precision    ::  NK

        complex*16        :: cnk
        double precision :: YcY
!
!        k1    = 0
!        k2    =-3.14159/DX
!        kdiff = 1.0
!
!        do while( kdiff > 1.0e-6  )
!            nk    = (k1+k2)/2
!            K_VEC = nk
!            cnk   = nk*II
!            call modpop_znajdz_wartosci_wlasne("N",cnk,Emax,.false.)
!            if ( W(num) < Ef  ) then
!                k2 = nk
!            else
!                k1 = nk
!            endif
!            kdiff = abs( (k2-k1)/k2)
!        enddo


		nk  = -K_M_IN(num)
        cnk = II*nk
        call modpop_znajdz_wartosci_wlasne("V",CNK,Emax,.false.)

        Chi_M_OUT(:,num)     = 0
        K_M_OUT(num)         = NK
        Chi_M_OUT(2:N+1,num) = Z(1:N,num)
        YcY              = sum(abs(Chi_M_OUT(:,num))**2)*DX
        Chi_M_OUT(:,num) = Chi_M_OUT(:,num)/sqrt(YcY)!/SQRT(2*3.14159*K_M_IN(num))
        print*,"EK_OUT(",num,")=",K_M_OUT(num)*L2LR,"[nm]"

!        Chi_M_IN(:,num)     = 0
!        K_M_IN(num)         = NK
!        Chi_M_IN(2:N+1,num) = Z(1:N,num)
!        YcY                 = sum(abs(Chi_M_IN(:,num))**2)*DX
!        Chi_M_IN(:,num)     = Chi_M_IN(:,num)/sqrt(YcY)!/SQRT(2*3.14159*NK)     ! dzielimy do rachunku zaburzen
!        print*,"EK_IN(",num,")=",K_M_IN(num)*L2LR,"[nm]"




    end subroutine modpop_znajdz_evan_k_out

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

!        if(bPoziome) then
!            Hamiltonian(1,2)   = -t0
!            Hamiltonian(1,1)   =  UVEC(1) + 4*t0 - 2*t0*cos(pk*DX + DX*DX*BZ*(1-hNY))
!        do i = 2 , N - 1
!            Hamiltonian(i,i-1) = -t0
!            Hamiltonian(i,i+1) = -t0
!            Hamiltonian(i,i)   =  UVEC(i) + 4*t0 - 2*t0*cos(pk*DX + DX*DX*BZ*(i-hNY) )
!        enddo
!            Hamiltonian(N,N-1) = -t0
!            Hamiltonian(N,N)   =  UVEC(N) + 4*t0 - 2*t0*cos(pk*DX + DX*DX*BZ*(N-hNY))
!        else
!            Hamiltonian(1,2)   = -t0
!            Hamiltonian(1,1)   =  UVEC(1) + 4*t0 - 2*t0*cos(pk*DX - DX*DX*BZ*(1-hNY))
!        do i = 2 , N - 1
!            Hamiltonian(i,i-1) = -t0
!            Hamiltonian(i,i+1) = -t0
!            Hamiltonian(i,i)   =  UVEC(i) + 4*t0 - 2*t0*cos(pk*DX - DX*DX*BZ*(i-hNY) )
!        enddo
!            Hamiltonian(N,N-1) = -t0
!            Hamiltonian(N,N)   =  UVEC(N) + 4*t0 - 2*t0*cos(pk*DX - DX*DX*BZ*(N-hNY))
!        endif


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



    subroutine modpop_liczba_podpasm(no_podpasm)
        integer , intent(inout) :: no_podpasm
        no_podpasm = LICZBA_MODOW

    end subroutine modpop_liczba_podpasm

    subroutine modpop_get_km(pLiczbaModow,pKm_in,pKm_out)
        integer, intent(in)                                    :: pLiczbaModow
        double precision,dimension(pLiczbaModow+NO_EVAN_MODES),intent(inout) :: pKm_in,pKm_out

        pKm_in (1:pLiczbaModow+NO_EVAN_MODES) = K_m_in (1:pLiczbaModow+NO_EVAN_MODES)
        pKm_out(1:pLiczbaModow+NO_EVAN_MODES) = K_m_out(1:pLiczbaModow+NO_EVAN_MODES)

    end subroutine modpop_get_km


    subroutine modpop_get_chi(pLiczbaModow,pN,pchi_m_in,pchi_m_out)
        integer, intent(in)                                 :: pLiczbaModow,pN
        complex*16,dimension(pN,pLiczbaModow+NO_EVAN_MODES),intent(inout) :: pchi_m_in,pchi_m_out

        pchi_m_in (1:pN,1:pLiczbaModow+NO_EVAN_MODES) = Chi_m_in (1:pN,1:pLiczbaModow+NO_EVAN_MODES)
        pchi_m_out(1:pN,1:pLiczbaModow+NO_EVAN_MODES) = Chi_m_out(1:pN,1:pLiczbaModow+NO_EVAN_MODES)

    end subroutine modpop_get_chi


END MODULE modpop
