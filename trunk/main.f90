
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
 integer,parameter :: nx            = 250;
 integer,parameter :: ny            = 250;
 integer,parameter :: liczba_zrodel = 4;
 double precision,parameter :: dx   = 4;
 double precision,dimension(:,:), allocatable  :: TR_MAT
 integer :: zwidth
 integer :: i,j
 doubleprecision :: qpc_w
 doubleprecision :: qpc2_w,qpc2_ypos
! -----------------------------------------------
! Zmienne wczytywane z config.ini
! -----------------------------------------------
call setIniFilename("config.ini")
call getDoubleValue("Dane","m_eff",M_EFF)
call getDoubleValue("Dane","eps",E_MAT)
call getDoubleValue("Dane","g_lan",G_LAN)
call getDoubleValue("Dane","Ef",atomic_Ef)
call getDoubleValue("Dane","Bz",atomic_Bz)
call getDoubleValue("Dane","qpc_w",qpc_w)

call modjed_ustaw_konwersje_jednostek(M_EFF,E_MAT);
call system_inicjalizacja(NX,NY,liczba_zrodel,DX);



call zrodla(1)%zrodlo_ustaw(ny/2-40,ny/2+40,1 ,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
call zrodla(2)%zrodlo_ustaw(50,ny-50,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
call zrodla(3)%zrodlo_ustaw(50,nx-80,1 ,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_GORA,UTOTAL)
call zrodla(4)%zrodlo_ustaw(50,nx-80,ny,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_DOL,UTOTAL)



!call system_dodaj_lorentza(15.0D0,10.0D0,10.0D0,400.0D0,550.0D0)


call system_dodaj_pionowy_slupek_potencjalu(47*dx,40*dx,ny/2*dx-qpc_w/2,50.0D0,30.0D0,3.0D0)
call system_dodaj_pionowy_slupek_potencjalu(47*dx,ny/2*dx+qpc_w/2,(ny-40)*dx,50.0D0,30.0D0,3.0D0)

qpc2_ypos = ny/2*dx - 0
qpc2_w = 50 ! nm
call system_dodaj_pionowy_slupek_potencjalu((nx-77)*dx,40*dx,qpc2_ypos-qpc2_w/2,50.0D0,30.0D0,3.0D0)
call system_dodaj_pionowy_slupek_potencjalu((nx-77)*dx,qpc2_ypos+qpc2_w/2,(ny-40)*dx,50.0D0,30.0D0,3.0D0)

call utworz_system()


call system_rozwiaz_problem(1,TR_MAT)

print*,sum(TR_MAT(1,:))
print*,sum(TR_MAT(2,:))
print*,sum(TR_MAT(3,:))
print*,sum(TR_MAT(4,:))
print*,sum(TR_MAT(:,:))
print*,TRANS_R,TRANS_T

call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI);
call system_zapisz_do_pliku("J.txt",ZAPISZ_J);


!call system_zapisz_do_pliku("flagi1.txt",ZAPISZ_FLAGI);
!call system_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)

deallocate(TR_MAT)

contains

subroutine utworz_system()
    integer :: i,j , wjazd , promien
    ! prosty test
    wjazd = 80

    GFLAGS(1:50,:)          = B_EMPTY
    GFLAGS(NX-wjazd:NX,:)   = B_EMPTY


    ! dajemy na brzegach dirichleta
    GFLAGS(1,:)  = B_DIRICHLET
    GFLAGS(NX,:) = B_DIRICHLET
    GFLAGS(:,1)  = B_DIRICHLET
    GFLAGS(:,NY) = B_DIRICHLET
    call system_inicjalizacja_ukladu(wjazd,0,10)
end subroutine utworz_system



end program transporter

