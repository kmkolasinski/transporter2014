
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
 integer,parameter :: nx            = 401;
 integer,parameter :: ny            = 201;
 integer,parameter :: liczba_zrodel = 2;
 double precision,parameter :: dx   = 2;
 integer :: zwidth
 double precision, dimension(:,:),allocatable :: Utmp
 integer :: i,j,riter,niter
 double precision :: T_SUM
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


zwidth = 15
call zrodla(1)%zrodlo_ustaw(NY/2-zwidth,NY/2+zwidth,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,2,UTOTAL)
call zrodla(2)%zrodlo_ustaw(NY/2-zwidth,NY/2+zwidth,NX,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,2,UTOTAL)
call utworz_system()
call system_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI);

! -------------------------------------------------------------------
!                       Ustawiamy zrodla
! -------------------------------------------------------------------
allocate(Utmp(nx,ny))
Utmp = UTOTAL


open(unit = 111, file= "TodIter.txt" )
!do atomic_Bz = 0.0,0.3,0.002
do atomic_Ef = 2.0,6.0,0.05
niter = 1
T_SUM = 0
riter = 1
!do riter = 1 , niter

UTOTAL = Utmp
do i = 100 , nx-100
do j = 1 , ny
    UTOTAL(i,j) = UTOTAL(i,j) + 2*(rand() - 0.5)*atomic_Ef/1000.0/Rd/3000.0
end do
end do


call zrodla(1)%zrodlo_ustaw(NY/2-zwidth,NY/2+zwidth,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,2,UTOTAL)
!call zrodla(1)%zrodlo_zapisz_mody("mody1.txt",dx)
call zrodla(2)%zrodlo_ustaw(NY/2-zwidth,NY/2+zwidth,NX,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,2,UTOTAL)
!call zrodla(2)%zrodlo_zapisz_mody("mody2.txt",dx)
call system_rozwiaz_problem(1)
T_SUM = T_SUM + TRANS_T
write(111,"(i10,4e20.8)"),riter,atomic_Ef,T_SUM/riter,TRANS_T+TRANS_R
enddo

close(111)
open(unit = 111, file= "Tave.txt" )
write(111,"(3e20.8)"),atomic_Ef,T_SUM/niter
close(111)

call system_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI);
call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI);
call system_zwalnienie_pamieci()
deallocate(Utmp)

contains

subroutine utworz_system()
    integer :: i,j , wjazd , promien
    ! prosty test
    wjazd = 210
    !GFLAGS(200+NX/2-20:200+NX/2,1:NY) = B_EMPTY

!    do i = 1 , nx
!    do j = 1 , ny
!        promien = sqrt((i - nx/2.0)**2 + (j - ny/2.0)**2)
!        if( promien > ny/2 .or. promien < ny/2-20  ) then
!             GFLAGS(i,j) = B_EMPTY
!        endif
!    enddo
!    enddo
    UTOTAL(NX/2-10-50:NX/2+10-50,:) = atomic_Ef/1000.0/Rd/5
    UTOTAL(NX/2-10+50:NX/2+10+50,:) = atomic_Ef/1000.0/Rd/5

    GFLAGS(1:wjazd,:)                 = B_EMPTY
    GFLAGS(NX-wjazd:NX,:)             = B_EMPTY

    ! dajemy na brzegach dirichleta
    GFLAGS(1,:)  = B_DIRICHLET
    GFLAGS(NX,:) = B_DIRICHLET
    GFLAGS(:,1)  = B_DIRICHLET
    GFLAGS(:,NY) = B_DIRICHLET
    call system_inicjalizacja_ukladu(wjazd,2,10)
end subroutine utworz_system


end program transporter

