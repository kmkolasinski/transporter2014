MODULE modjed
    USE, INTRINSIC :: ISO_C_BINDING
    implicit none


    double precision,parameter  :: M_PI      = 3.1415926
    ! promien bohra w nm
    double precision,parameter  :: A0        = 0.0529177249
    ! jednostka energi (atomowych) eV
    double precision,parameter  :: E0        = 27.211384523
    complex*16 ,parameter       :: II        = CMPLX(0,1)
    double precision  :: KbT
    double precision  :: Temp
    double precision  :: Rd
    double precision  :: L2LR
    double precision  :: LR2L
    double precision  :: M_EFF ! masa efektywna
    double precision  :: E_MAT ! przenikalnosc
    double precision  :: G_LAN ! czynnik landego
    double precision :: atomic_Ef   ! energia elektronu [meV]
    double precision :: atomic_Bz   ! polemagnetyczne w [T]
    double precision :: atomic_Bx   ! polemagnetyczne w [T]
    double precision :: atomic_By   ! polemagnetyczne w [T]
    double precision :: atomic_Rashba   ! sprzezenie spin-orbita typu rashba [nm^2]
    double precision :: atomic_LOC      ! sprzezenie spin-orbita typu lateral [nm^2]
    double precision :: atomic_DX
    double precision :: so_rashba,so_loc ! w jednostkach donorowych
    double precision :: so_Fz
    double precision :: so_alpha3D

    logical,parameter :: TRANS_DEBUG = .false.

    ENUM,BIND(C)
        ENUMERATOR :: B_NORMAL          = 0 ! FLAGA OZNACZA ZE TO NIE JEST BRZEG
        ENUMERATOR :: B_DIRICHLET       = 2
        ENUMERATOR :: B_WEJSCIE         = 16
        ENUMERATOR :: B_TRANSPARENT     = 32
        ENUMERATOR :: B_EMPTY           = -1
    END ENUM

    ENUM,BIND(C)
        ENUMERATOR :: ZRODLO_KIERUNEK_PRAWO= 0
        ENUMERATOR :: ZRODLO_KIERUNEK_LEWO = 1
        ENUMERATOR :: ZRODLO_KIERUNEK_GORA = 2
        ENUMERATOR :: ZRODLO_KIERUNEK_DOL  = 3
    END ENUM

    contains
    subroutine modjed_jaki_kierunek(kierunek)
        integer :: kierunek
        selectcase(kierunek)
        case(ZRODLO_KIERUNEK_PRAWO)
            print*,"    Kierunek= PRAWO -->"
        case(ZRODLO_KIERUNEK_LEWO)
            print*,"    Kierunek= LEWO <--"
        case(ZRODLO_KIERUNEK_GORA)
            print*,"    Kierunek= GORA -->"
        case(ZRODLO_KIERUNEK_DOL)
            print*,"    Kierunek= DOL <--"
        endselect
    endsubroutine modjed_jaki_kierunek

    subroutine modjed_ustaw_InGaAs()
        ! na podstawie -- http://www.ioffe.ru/SVA/NSM/Semicond/GaInAs/
        ! lande factor -- http://arxiv.org/pdf/1406.2848v1.pdf
        ! http://journals.aps.org/prb/abstract/10.1103/PhysRevB.35.7729

        G_LAN      = -8.97
        so_Fz      = 20.0 ! w meV/nm
        so_alpha3D = 0.572 ! w nm^2

        call modjed_ustaw_konwersje_jednostek(0.0465D0,12.0D0)
    end subroutine modjed_ustaw_InGaAs

    subroutine modjed_ustaw_InSb()
        ! na podstawie pracy doktorskiej (colwiz zakladka SO)
        G_LAN = -50.0
!        G_LAN = -10.0
        call modjed_ustaw_konwersje_jednostek(0.014D0,16.0D0)
    end subroutine modjed_ustaw_InSb

    subroutine modjed_ustaw_konwersje_jednostek(pM_EFF,pE_MAT)
        double precision,intent(in) :: pM_EFF
        double precision,intent(in) :: pE_MAT

        print*,"Zmiana jednostek..."

        E_MAT     = pE_MAT
        M_EFF     = pM_EFF
        Temp      = 0.0
        Rd        = E0*M_EFF/E_MAT/E_MAT
        L2LR      = M_EFF/E_MAT/A0
        LR2L      = 1.0/L2LR
        kbT       = 8.617D-5*Temp/Rd

        atomic_Rashba = so_alpha3D*so_Fz
        atomic_LOC    = so_alpha3D
        so_rashba     = atomic_Rashba * L2LR  / 1000.0 / Rd
        so_loc        = atomic_LOC * L2LR * L2LR
        print*,"Rashba  =",so_rashba," [donorowe]"
        print*,"Lateral =",so_loc   ," [donorowe]"

    end subroutine modjed_ustaw_konwersje_jednostek



    ! Funkcja konwertujaca pole magnetyczne w teslach do
    ! pola w jednostkach donorowych.
    double precision function BtoDonorB(inBZ) result(rval)
      double precision, intent(in) :: inBZ
      rval = 2*inBZ/m_eff*5.78838e-5/Rd
    endfunction BtoDonorB

    double precision function DonorBtoB(inBZ) result(rval)
      double precision, intent(in) :: inBZ
      rval = inBZ*m_eff/5.78838e-5*Rd/2
    endfunction DonorBtoB

END MODULE modjed
