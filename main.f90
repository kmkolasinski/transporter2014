
program transporter
 use modutils
 use modpop
 use modjed
 use modinip
 use modsystem
 implicit none


! -----------------------------------------------
! Deklaracje zmiennych i tablic
! -----------------------------------------------
 integer,parameter :: nx            = 400;
 integer,parameter :: ny            = 50;
 integer,parameter :: liczba_zrodel = 2;
 double precision,parameter :: dx   = 2;

! -----------------------------------------------
! Zmienne wczytywane z config.ini
! -----------------------------------------------
 double precision :: meff ! masa efektywna
 double precision :: eps  ! przenikalnosc
 double precision :: glan ! czynnik landego
 double precision :: Ef   ! energia elektronu [meV]
 double precision :: Bz   ! polemagnetyczne w [T]


! -----------------------------------------------
! Wczytywanie ustawien
! -----------------------------------------------
call setIniFilename("config.ini")
call getDoubleValue("Dane","m_eff",meff)
call getDoubleValue("Dane","eps",eps)
call getDoubleValue("Dane","g_lan",glan)
call getDoubleValue("Dane","Ef",Ef)
call getDoubleValue("Dane","Bz",Bz)


call modjed_ustaw_konwersje_jednostek(meff,eps);
call system_inicjalizacja(NX,NY,liczba_zrodel,DX);

! ---------------------------------------
! Ustawiamy zrodla
! ---------------------------------------
call zrodla(1)%zrodlo_ustaw(1,NY,1,dx,ef,bz,.true.,.true.,2,UTOTAL)
call zrodla(1)%zrodlo_zapisz_mody("mody1.txt",dx)

call zrodla(2)%zrodlo_ustaw(1,NY,NX,dx,ef,bz,.false.,.false.,2,UTOTAL)
call zrodla(2)%zrodlo_zapisz_mody("mody2.txt",dx)

call utworz_system()

call system_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI);
call system_zapisz_do_pliku("indeksy.txt",ZAPISZ_INDEKSY);



call system_zwalnienie_pamieci()
contains

subroutine utworz_system()
    integer :: i,j
    ! prosty test

    !GFLAGS(NX/2-20:NX/2+20,1:NY/2) = B_EMPTY
    !GFLAGS(1:80,:)                 = B_EMPTY
    !GFLAGS(NX-80:NX,:)             = B_EMPTY

    ! dajemy na brzegach dirichleta
    GFLAGS(1,:)  = B_DIRICHLET
    GFLAGS(NX,:) = B_DIRICHLET
    GFLAGS(:,1)  = B_DIRICHLET
    GFLAGS(:,NY) = B_DIRICHLET
    call system_inicjalizacja_ukladu(80,2,10)
end subroutine utworz_system


end program transporter

