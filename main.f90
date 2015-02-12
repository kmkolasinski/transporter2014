
program transporter
 use modutils
! use modpop
! use spinmodpop
 use modspinzrodlo
 use modspinsystem
 use modjed
 use modinip
! use modsystem
 use ifport
 implicit none


! -----------------------------------------------
! Deklaracje zmiennych i tablic
! -----------------------------------------------
 integer :: nx            = 50;
 integer :: ny            = 50;
 integer :: liczba_zrodel = 2;
 double precision :: dx   = 2 , pdx
 double precision,dimension(:,:), allocatable  :: TR_MAT
 integer :: i,j


! -----------------------------------------------
! Zmienne wczytywane z config.ini
! -----------------------------------------------
call setIniFilename("config.ini")
call getIntValue("Dane","nx",nx)
call getIntValue("Dane","ny",ny)
call getDoubleValue("Dane","dx",atomic_DX)
call getDoubleValue("Dane","m_eff",M_EFF)
call getDoubleValue("Dane","eps",E_MAT)
call getDoubleValue("Dane","g_lan",G_LAN)
call getDoubleValue("Dane","Ef",atomic_Ef)
call getDoubleValue("Dane","Bz",atomic_Bz)
call getDoubleValue("Dane","rashba",atomic_Rashba)
call getDoubleValue("Dane","loc",atomic_LOC)


call modjed_ustaw_konwersje_jednostek(M_EFF,E_MAT);


dx = atomic_DX
call spinsystem_inicjalizacja(NX,NY,liczba_zrodel);


UTOTAL= 0
do i = 1 , ny
!    UTOTAL(:,i) = (i-0.0*ny/2.0-0.5)**2*0.00-0.0
!    write(444,*),i,UTOTAL(1,i)
enddo

print*,"asd1",get_clock()
call reset_clock()
!call zrodla(1)%spinzrodlo_ustaw(3,30,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%spinzrodlo_ustaw(NY-30,NY-3,nx-2,ZRODLO_KIERUNEK_LEWO,UTOTAL)

call zrodla(1)%spinzrodlo_ustaw(3,NY-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
call zrodla(2)%spinzrodlo_ustaw(3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)

!call zrodla(1)%spinzrodlo_relacja_dyspersji(-0.5D0,0.5D0,0.01D0,2*atomic_Ef,"rel.txt")
!call zrodla(1)%spinzrodlo_zapisz_mody("modup.txt","moddown.txt",.true.)
!call zrodlo%spinzrodlo_relacja_dyspersji(-0.2D0,+0.2D0,0.001D0,4*atomic_Ef,"relacja.txt")

print*,"asd2",get_clock()
call reset_clock()
!call zrodla(1)%spinzrodlo_wypisz_JinJout()
!call zrodlo%spinzrodlo_zapisz_mody("modsup.txt","modsdown.txt",.true.);



print*,"asd3",get_clock()
call utworz_system()
call spinsystem_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI)
call spinsystem_zapisz_do_pliku("indeksy.txt",ZAPISZ_INDEKSY)

call spinsystem_rozwiaz_problem(1,TR_MAT)
call spinsystem_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
call spinsystem_zwalnienie_pamieci()

stop
!! -------------------------------- Skan energii --------------------------------
!open(unit=222,file="T.txt")
!do atomic_Ef = 1.0 , 5.0 , 0.01
!
!call zrodla(1)%zrodlo_ustaw(15,ny-15,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%zrodlo_ustaw(15,ny-15,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!UTOTAL = 0
!call system_dodaj_lorentza(2.0D0,10.0D0,10.0D0,nx/2*dx,ny/2*dx)
!
!call utworz_system()
!call system_rozwiaz_problem(1,TR_MAT)
!
!TRANS_T = sum(TR_MAT(2,:))
!TRANS_R = sum(TR_MAT(1,:))
!
!write(222,"(20e20.8)"),atomic_Ef,atomic_Bz,TRANS_T,TRANS_R,TRANS_T+TRANS_R
!
!enddo
!close(222)
!call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
!stop

!
!
!!do atomic_Ef = 0.0 , 0.5 , 0.5/100.0
!call zrodla(1)%zrodlo_ustaw(3,ny-2,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!
!call zrodla(1)%zrodlo_zapisz_mody("mody.txt",dx,.true.)
!call zrodla(2)%zrodlo_ustaw(3,ny-2,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!pdx = dx*L2LR
!UTOTAL(nx/2:nx/2,1:ny/2-10) = 15 * 0.5 / pdx / pdx
!UTOTAL(nx/2:nx/2,ny/2+30:ny) = 15 * 0.5 / pdx / pdx
!call utworz_system()
!!call system_dodaj_lorentza(10.0D0,30.0D0,30.0D0,nx/2*dx,ny/4*dx)
!call system_rozwiaz_problem(1,TR_MAT)
!TRANS_T = sum(TR_MAT(2,:))
!TRANS_R = sum(TR_MAT(1,:))
!write(333,*),atomic_Ef,TRANS_T
!!enddo
!call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI);
!call system_zapisz_do_pliku("j.txt",ZAPISZ_J);
!call system_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL);
!call system_zapisz_do_pliku("kontur.txt",ZAPISZ_FLAGI);
!
!
!
!call system_zwalnienie_pamieci()
!if(allocated(TR_MAT))deallocate(TR_MAT)

contains


subroutine utworz_system()
    integer :: i,j , wjazd , promien
    ! prosty test
    wjazd = nx/2-15
    GFLAGS = B_EMPTY
    GFLAGS(nx/4:nx-nx/4,3:ny-3)  = B_NORMAL
    !GFLAGS(1:55,:)       = B_EMPTY
    GFLAGS(NX/2-15:NX/2+15,1:NY/2)   = B_EMPTY
!
!    GFLAGS(:,1:14)       = B_EMPTY
!    GFLAGS(:,ny-14:ny)   = B_EMPTY



    ! dajemy na brzegach dirichleta
    GFLAGS(1,:)  = B_DIRICHLET
    GFLAGS(NX,:) = B_DIRICHLET
    GFLAGS(:,1)  = B_DIRICHLET
    GFLAGS(:,NY) = B_DIRICHLET
    call spinsystem_inicjalizacja_ukladu(wjazd,0,0)
end subroutine utworz_system



end program transporter

