
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
 integer,parameter :: liczba_zrodel = 6;
 double precision,parameter :: dx   = 4;
 double precision,dimension(:,:), allocatable  :: TR_MAT
 integer :: zwidth
 integer :: i,j
 doubleprecision :: qpc_w
 doubleprecision :: qpc2_w,qpc2_ypos
 integer :: input_w , in_w
 integer :: output_w , out_w
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





input_w = 60
in_w = (ny - input_w - 20)/2

output_w = 80
out_w = (ny - output_w - 20)/2

print*, "input_w=", input_w
print*, "in_w   =", in_w

call system_dodaj_abs_zrodlo(1,nx,1 ,atomic_Ef,ZRODLO_KIERUNEK_GORA)
call system_dodaj_abs_zrodlo(1,nx,ny,atomic_Ef,ZRODLO_KIERUNEK_DOL )


open(unit=222,file="T.txt")
!do atomic_Bz = 0 , 0.3 , 0.002

call zrodla(1)%zrodlo_ustaw(ny/2-input_w/2,ny/2+input_w/2,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
call zrodla(4)%zrodlo_ustaw(2,in_w,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)



call zrodla(3)%zrodlo_ustaw(ny-in_w,ny-1,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)

call zrodla(2)%zrodlo_ustaw(ny/2-output_w/2,ny/2+output_w/2,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
call zrodla(5)%zrodlo_ustaw(2,out_w,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
call zrodla(6)%zrodlo_ustaw(ny-out_w,ny-1,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)






!call system_dodaj_lorentza(15.0D0,10.0D0,10.0D0,400.0D0,550.0D0)

UTOTAL = 0
call system_dodaj_pionowy_slupek_potencjalu(47*dx,(ny/2-input_w/2)*dx,ny/2*dx-qpc_w/2,50.0D0,30.0D0,3.0D0)
call system_dodaj_pionowy_slupek_potencjalu(47*dx,(ny/2+1)*dx+qpc_w/2,(ny/2+input_w/2)*dx,50.0D0,30.0D0,3.0D0)

qpc2_ypos = ny/2*dx - 0
qpc2_w = 80 ! nm
call system_dodaj_pionowy_slupek_potencjalu((nx-77)*dx,(ny/2-output_w/2)*dx,qpc2_ypos-qpc2_w/2,50.0D0,30.0D0,3.0D0)
call system_dodaj_pionowy_slupek_potencjalu((nx-77)*dx,qpc2_ypos+qpc2_w/2,(ny/2+output_w/2)*dx,50.0D0,30.0D0,3.0D0)





call utworz_system()
call system_rozwiaz_problem(1,TR_MAT)


TRANS_T = sum(TR_MAT(2,:))
TRANS_R = sum(TR_MAT(1,:))


write(222,"(20e20.8)"),atomic_Ef,atomic_Bz,qpc_w,TRANS_T,TRANS_R,zrodla(1)%liczba_modow - TRANS_T - TRANS_R , zrodla(1)%liczba_modow + 0.0



!enddo
close(222)


call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI);
!call system_zapisz_do_pliku("J.txt",ZAPISZ_J);
!call system_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI);
!call system_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)

deallocate(TR_MAT)

contains

subroutine utworz_system()
    integer :: i,j , wjazd , promien
    ! prosty test
    wjazd = 100

    GFLAGS(1:50,:)       = B_EMPTY
    GFLAGS(NX-80:NX,:)   = B_EMPTY
!    GFLAGS(:,1:wjazd)       = B_EMPTY
!    GFLAGS(:,ny-wjazd:ny)   = B_EMPTY

    ! dajemy na brzegach dirichleta
    GFLAGS(1,:)  = B_DIRICHLET
    GFLAGS(NX,:) = B_DIRICHLET
    GFLAGS(:,1)  = B_DIRICHLET
    GFLAGS(:,NY) = B_DIRICHLET
    call system_inicjalizacja_ukladu(wjazd,0,10)
end subroutine utworz_system



end program transporter

