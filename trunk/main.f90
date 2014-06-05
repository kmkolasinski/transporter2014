
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
 integer :: i,j,riter,niter,giter,iter
 double precision :: T_SUM
 doubleprecision :: xp,yp,x0,y0,amp,sigma
 character(len=11) :: plik
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
iter = 1
do atomic_Bz = 0.0,0.10,0.01
call zrodla(1)%zrodlo_ustaw(1,NY,1,dx,atomic_Ef,-atomic_Bz,ZRODLO_KIERUNEK_PRAWO,2,UTOTAL)

if(iter < 10) then
    write(plik(1:1),"(i1)"),iter
    plik(2:9) = "-ave.txt"
else if(iter < 100) then
    write(plik(1:2),"(i2)"),iter
    plik(3:10) = "-ave.txt"
else if(iter < 1000) then
    write(plik(1:3),"(i3)"),iter
    plik(4:11) = "-ave.txt"
endif
print*,"dane/"//plik
call zrodla(1)%zrodlo_zapisz_mody("dane/"//plik,dx)

iter = iter +  1
enddo
stop
call zrodla(1)%zrodlo_ustaw(NY/2-zwidth,NY/2+zwidth,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,2,UTOTAL)
call zrodla(2)%zrodlo_ustaw(NY/2-zwidth,NY/2+zwidth,NX,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,2,UTOTAL)
call utworz_system()
call system_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI);

! -------------------------------------------------------------------
!                       Ustawiamy zrodla
! -------------------------------------------------------------------
allocate(Utmp(nx,ny))
Utmp = UTOTAL


open(unit = 111, file= "TodIter4.txt" )
!do atomic_Bz = 0.0,0.3,0.002
do atomic_Bz = 0.2,0.4,0.002
niter = 50
T_SUM = 0
riter = 1
do riter = 1 , niter



UTOTAL = 0
do giter  = 1, 5
    x0  = rand()*(nx-100)+50
    y0  = rand()*(ny)
    amp   = -0.5*atomic_Ef/1000.0/Rd
    sigma = 0.1

    do i = 50 , nx-50
    do j = 1 , ny
        xp = i
        yp = j
       UTOTAL(i,j) = UTOTAL(i,j) + system_gauss(xp,yp,x0,y0,sigma,amp)
    end do
    end do
enddo

call zrodla(1)%zrodlo_ustaw(NY/2-zwidth,NY/2+zwidth,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,2,UTOTAL)
!call zrodla(1)%zrodlo_zapisz_mody("mody1.txt",dx)
call zrodla(2)%zrodlo_ustaw(NY/2-zwidth,NY/2+zwidth,NX,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,2,UTOTAL)
!call zrodla(2)%zrodlo_zapisz_mody("mody2.txt",dx)
call system_rozwiaz_problem(1)
T_SUM = T_SUM + TRANS_T
enddo
write(111,"(i10,4e20.8)"),riter,atomic_Bz,T_SUM/niter,TRANS_T+TRANS_R
enddo

close(111)
open(unit = 111, file= "Tave.txt" )
write(111,"(3e20.8)"),atomic_Ef,T_SUM/niter
close(111)

call system_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI);
call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI);
call system_zapisz_do_pliku("utotal.txt",ZAPISZ_POTENCJAL);
call system_zwalnienie_pamieci()
deallocate(Utmp)

contains

subroutine utworz_system()
    integer :: i,j , wjazd , promien
    ! prosty test
    wjazd = 110
    !GFLAGS(200+NX/2-20:200+NX/2,1:NY) = B_EMPTY

    do i = 1 , nx
    do j = 1 , ny
        promien = sqrt((i - nx/2.0)**2 + (j - ny/2.0)**2)
        if( promien > ny/2 .or. promien < ny/2-40  ) then
             GFLAGS(i,j) = B_EMPTY
        endif
    enddo
    enddo
    !UTOTAL(NX/2-150:NX/2+150,:) = -atomic_Ef/1000.0/Rd*5
    !UTOTAL(NX/2-10+50:NX/2+10+50,:) = atomic_Ef/1000.0/Rd/5

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

