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
    double precision  :: M_EFF
    double precision  :: E_MAT
    double precision  :: G_LAN

    logical,parameter :: TRANS_DEBUG = .true.

    ENUM,BIND(C)
        ENUMERATOR :: B_NORMAL          = 0 ! FLAGA OZNACZA ZE TO NIE JEST BRZEG
        ENUMERATOR :: B_DIRICHLET       = 2
        ENUMERATOR :: B_WEJSCIE         = 16
        ENUMERATOR :: B_EMPTY           = -1
    END ENUM

    contains


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
