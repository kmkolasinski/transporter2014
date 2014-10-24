
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
 integer :: nx            = 250;
 integer :: ny            = 50;
 integer :: liczba_zrodel = 2;
 double precision :: dx   = 4;
 double precision,dimension(:,:), allocatable  :: TR_MAT
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

open(unit = 4356, file= "stabilizacja.txt" )
do i = 60 , 150

nx = 380/4 + 5
call modjed_ustaw_konwersje_jednostek(M_EFF,E_MAT);
call system_inicjalizacja(NX,NY,liczba_zrodel,DX);

call zrodla(1)%zrodlo_ustaw(15,ny-15,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
call zrodla(2)%zrodlo_ustaw(15,ny-15,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
UTOTAL = 0
call utworz_system()
call system_dodaj_pionowy_slupek_potencjalu((nx/2.0-20.)*dx,0.0D0,ny*dx,2.0D0,25.0D0,3.0d0)
call system_dodaj_pionowy_slupek_potencjalu((nx/2.0+20.)*dx,0.0D0,ny*dx,2.0D0,25.0D0,3.0d0)

call system_widmo(3.0D0,4.0D0,100,8,1,25)
call system_zapisz_widmo_do_pliku("widmo.txt",ZAPISZ_STANY_WLASNE)


write(4356,"(100e20.6)"),i*dx,DBLE(Widmo_Evals(1:Widmo_NoStates))*1000.0*Rd
call system_zwalnienie_pamieci()
enddo
close(4356)
stop

open(unit=222,file="T.txt")
!do atomic_Bz = 0 , 0.3 , 0.002

!call system_rozwiaz_problem(1,TR_MAT)
!TRANS_T = sum(TR_MAT(2,:))
!TRANS_R = sum(TR_MAT(1,:))
!
!!write(222,"(20e20.8)"),atomic_Ef,atomic_Bz,TRANS_T,TRANS_R
!!close(222)
!
!call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI);
!call system_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI);
call system_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
!
!!!call system_zapisz_do_pliku("J.txt",ZAPISZ_J);


call system_zwalnienie_pamieci()
if(allocated(TR_MAT))deallocate(TR_MAT)

contains

subroutine utworz_system()
    integer :: i,j , wjazd , promien
    ! prosty test
    wjazd = 15

    GFLAGS(1:15,:)       = B_EMPTY
    GFLAGS(NX-15:NX,:)   = B_EMPTY

    GFLAGS(:,1:14)       = B_EMPTY
    GFLAGS(:,ny-14:ny)   = B_EMPTY

!    GFLAGS(:,1:wjazd)       = B_EMPTY
!    GFLAGS(:,ny-wjazd:ny)   = B_EMPTY

    ! dajemy na brzegach dirichleta
    GFLAGS(1,:)  = B_DIRICHLET
    GFLAGS(NX,:) = B_DIRICHLET
    GFLAGS(:,1)  = B_DIRICHLET
    GFLAGS(:,NY) = B_DIRICHLET
    call system_inicjalizacja_ukladu(wjazd,0,0)
end subroutine utworz_system



end program transporter

