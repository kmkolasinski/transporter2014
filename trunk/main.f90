
program transporter
 use modutils
 use modpop
 use modjed
 use modinip
 use modsystem
 use ifport
 implicit none


! -----------------------------------------------
! Deklaracje zmiennych i tablic
! -----------------------------------------------
 integer,parameter :: nx            = 201;
 integer,parameter :: ny            = 201;
 integer,parameter :: liczba_zrodel = 3;
 double precision,parameter :: dx   = 2;
 integer :: zwidth
 integer :: i,j
! -----------------------------------------------
! Zmienne wczytywane z config.ini
! -----------------------------------------------
call setIniFilename("config.ini")
call getDoubleValue("Dane","m_eff",M_EFF)
call getDoubleValue("Dane","eps",E_MAT)
call getDoubleValue("Dane","g_lan",G_LAN)
call getDoubleValue("Dane","Ef",atomic_Ef)
call getDoubleValue("Dane","Bz",atomic_Bz)


call modjed_ustaw_konwersje_jednostek(M_EFF,E_MAT);
call system_inicjalizacja(NX,NY,liczba_zrodel,DX);


zwidth = 10

call zrodla(1)%zrodlo_ustaw(NY/2-zwidth+50,NY/2+zwidth+50,1 ,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
call zrodla(3)%zrodlo_ustaw(50,NY-50,nx ,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!call zrodla(4)%zrodlo_ustaw(NY/2-zwidth-50,NY/2+zwidth-50,1 ,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%zrodlo_ustaw(NY/2+NY/4,NY-10,NX,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!call zrodla(3)%zrodlo_ustaw(NY/2-zwidth-50,NY/2+zwidth-50,NX,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,2,UTOTAL)

!call zrodla(2)%zrodlo_ustaw(NX/2-zwidth,NX/2+zwidth,10,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_GORA,UTOTAL)
call zrodla(2)%zrodlo_ustaw(50,150,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_GORA,UTOTAL)


call utworz_system()

call system_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI);
call system_zapisz_do_pliku("indeksy.txt",ZAPISZ_INDEKSY);

call system_rozwiaz_problem(1)
call system_zapisz_do_pliku("phi1.txt",ZAPISZ_PHI);

call zrodla(1)%zrodlo_zapisz_mody("mod1.txt",dx,.true.)
call zrodla(2)%zrodlo_zapisz_mody("mod2.txt",dx,.true.)


contains

subroutine utworz_system()
    integer :: i,j , wjazd , promien
    ! prosty test
    wjazd = 50
    !GFLAGS(200+NX/2-20:200+NX/2,1:NY) = B_EMPTY

!    do i = 1 , nx
!    do j = 1 , ny
!        promien = sqrt((i - nx/2.0)**2 + (j - ny/2.0)**2)
!        if( promien > ny/2   ) then
!!        if( promien > ny/2 .or. promien < ny/2-40  ) then
!             GFLAGS(i,j) = B_EMPTY
!        endif
!    enddo
!    enddo

    !UTOTAL(NX/2-10:NX/2+20,NY/2-20+50:NY/2+50) = atomic_Ef/1000.0/Rd

    GFLAGS(1:wjazd,:)                 = B_EMPTY
    !GFLAGS(NX-wjazd:NX,:)             = B_EMPTY

    GFLAGS(:,1:wjazd)                 = B_EMPTY

    ! dajemy na brzegach dirichleta
    GFLAGS(1,:)  = B_DIRICHLET
    GFLAGS(NX,:) = B_DIRICHLET
    GFLAGS(:,1)  = B_DIRICHLET
    GFLAGS(:,NY) = B_DIRICHLET
    call system_inicjalizacja_ukladu(wjazd,2,10)
end subroutine utworz_system


end program transporter

