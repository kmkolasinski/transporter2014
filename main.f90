
program transporter
 use modutils
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
 double precision :: dx   = 2 , pdx , omega ,omega2 , x , y , gamma , xpos
 double precision,dimension(:,:), allocatable  :: TR_MAT
 integer :: i,j
 double precision :: width , G21(-1:1) , G23(-1:1)


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
call getDoubleValue("Dane","omega",omega)
call getDoubleValue("Dane","omega2",omega2)
call getDoubleValue("Dane","gamma",gamma)
call getDoubleValue("Dane","xpos",xpos)


!call modjed_ustaw_konwersje_jednostek(M_EFF,E_MAT);
!call modjed_ustaw_InGaAs()
call modjed_ustaw_InSb()
dx = atomic_DX
call spinsystem_inicjalizacja(NX,NY,liczba_zrodel);


UTOTAL = 0
do i = 1 , nx
do j = 1 , ny
    x = i * dx
    y = j * dx
    UTOTAL(i,j) = ( 0.5*(omega/1000.0/Rd)**2*((y-ny*dx/2)**2) + (omega2/1000.0/Rd)*(y-ny*dx/2) )* &
                   exp( -0.5*( x - xpos )**2/(gamma)**2 )
enddo
enddo

!call zrodla(1)%spinzrodlo_ustaw(3,NY/2-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%spinzrodlo_ustaw(NY/2-3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)
call reset_clock()

call zrodla(1)%spinzrodlo_ustaw(3,NY-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
call zrodla(2)%spinzrodlo_ustaw(3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)

!call zrodla(3)%spinzrodlo_ustaw(NY/2-20+40,NY/2+20+40,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)

!call zrodla(1)%spinzrodlo_relacja_dyspersji(-0.5D0,0.5D0,0.001D0,atomic_Ef*2,"rel1.txt")
!call spinsystem_dodaj_lorentza(5.0D0,50.0D0,50.0D0,nx/2*atomic_DX,ny/2*atomic_DX)

call utworz_system()



    open(unit = 222, file= "T.txt" )

    call spinsystem_rozwiaz_problem(1,TR_MAT)
!    G21(+1) = TR_MAT(1,1)
!    G21(-1) = TR_MAT(1,2)
!    G23(+1) = TR_MAT(3,1)
!    G23(-1) = TR_MAT(3,2)
!    print*,"G21(+1)=",G21(+1)
!    print*,"G21(-1)=",G21(-1)
!
!    print*,"G23(+1)=",G23(+1)
!    print*,"G23(-1)=",G23(-1)
!    print*,"W=",width,"R=",sum(TR_MAT(1,:)),"T=",sum(TR_MAT(2,:))
    write(222,*),omega,sum(TR_MAT(2,:)),sum(TR_MAT(1,:))
!    write(222,"(20f16.8)"),omega,G21(+1),G21(-1),G23(+1),G23(-1)

    close(222)
!call spinsystem_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
!call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
call spinsystem_zapisz_do_pliku("j.txt",ZAPISZ_J_ALL)
!call spinsystem_zapisz_do_pliku("kon.txt",ZAPISZ_KONTUR)
!call spinsystem_zapisz_do_pliku("divj.txt",ZAPISZ_DIVJ)
!call spinsystem_zapisz_do_pliku("pol.txt",ZAPISZ_POLARYZACJE)


call spinsystem_zwalnienie_pamieci()
if(allocated(TR_MAT))deallocate(TR_MAT)


!call system_inicjalizacja(NX,NY,liczba_zrodel,atomic_DX);
!
!call system_dodaj_abs_zrodlo(nx/2-15,nx-1,1,atomic_Ef,ZRODLO_KIERUNEK_GORA)
!call system_dodaj_abs_zrodlo(nx/2-15,nx-1,ny-2,atomic_Ef,ZRODLO_KIERUNEK_DOL)
!
!UTOTAL= 0
!
!call zrodla(1)%zrodlo_ustaw(NY/2-5,NY/2+5,1 ,atomic_DX,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%zrodlo_ustaw(2      ,NY-1   ,nx,atomic_DX,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO ,UTOTAL)
!
!
!!call system_dodaj_lorentza(10.0D0,50.0D0,50.0D0,250.0D0,400.0D0)
!
!
!
!call utworz_system()
!!call spinsystem_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI)
!!call spinsystem_zapisz_do_pliku("indeksy.txt",ZAPISZ_INDEKSY)
!
!width = 100
!!do width = 0.0 , 100.0 , 1.0D0
!!    UTOTAL = 0
!!    call spinsystem_dodaj_pionowy_slupek_potencjalu( 100.0D0,-100.0D0,300.0D0-width/2,14.0D0,50.0D0,4.3D0)
!!    call spinsystem_dodaj_pionowy_slupek_potencjalu( 100.0D0,300.0D0+width/2,700.0D0,14.0D0,50.0D0,4.3D0)
!
!    call system_rozwiaz_problem(1,TR_MAT)
!    print*,"W=",width,"R=",sum(TR_MAT(1,:)),"T=",sum(TR_MAT(2,:))
!    write(222,*),width,sum(TR_MAT(1,:)),sum(TR_MAT(2,:))
!!enddo
!
!call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
!call system_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
!call system_zapisz_do_pliku("j.txt",ZAPISZ_J)
!call system_zapisz_do_pliku("kon.txt",ZAPISZ_KONTUR)
!call system_zwalnienie_pamieci()
!if(allocated(TR_MAT))deallocate(TR_MAT)



contains


subroutine utworz_system()
    integer :: i,j , wjazd , promien
    ! prosty test
    wjazd  = nx/4+20
    GFLAGS = B_EMPTY
!    GFLAGS(nx/8:nx-nx/4,3:NY-3)  = B_NORMAL


    wjazd  = nx/2+20
    GFLAGS = B_EMPTY

!    GFLAGS(nx/8:nx-nx/4,3:NY-3)  = B_NORMAL


!    GFLAGS(nx/8:nx/2,ny/2-35:ny/2+35)  = B_NORMAL
    !GFLAGS(1:55,:)       = B_EMPTY
!    GFLAGS(NX/2-15:NX/2+15,1:NY/2-10)   = B_EMPTY
    !GFLAGS(NX/2-15:NX/2+15,NY/2+10:NY)   = B_EMPTY
!
!    GFLAGS(:,1:14)       = B_EMPTY
!    GFLAGS(:,ny-14:ny)   = B_EMPTY



    ! dajemy na brzegach dirichleta
    GFLAGS(1,:)  = B_DIRICHLET
    GFLAGS(NX,:) = B_DIRICHLET
    GFLAGS(:,1)  = B_DIRICHLET
    GFLAGS(:,NY) = B_DIRICHLET
!    call system_inicjalizacja_ukladu(wjazd,2,4)
    call spinsystem_inicjalizacja_ukladu(wjazd,4,4)
end subroutine utworz_system



end program transporter

