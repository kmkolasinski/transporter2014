
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
doubleprecision :: xpos , deltax , x1,x2
doubleprecision :: la1,lx1,ly1,la2,lx2,ly2
doubleprecision :: ypos

    la1 = 5
    lx1 = 15
    ly1 = 15
    la2 = 5
    lx2 = 15
    ly2 = 15
! -----------------------------------------------
! Zmienne wczytywane z config.ini
! -----------------------------------------------
call setIniFilename("config.ini")
call getDoubleValue("Dane","m_eff",M_EFF)
call getDoubleValue("Dane","eps",E_MAT)
call getDoubleValue("Dane","g_lan",G_LAN)
call getDoubleValue("Dane","Ef",atomic_Ef)
call getDoubleValue("Dane","Bz",atomic_Bz)
call getDoubleValue("Dane","ypos",ypos)



call modjed_ustaw_konwersje_jednostek(M_EFF,E_MAT);
call system_inicjalizacja(NX,NY,liczba_zrodel,DX);


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


call zrodla(1)%zrodlo_ustaw(15,ny-15,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
call zrodla(2)%zrodlo_ustaw(15,ny-15,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
call utworz_system()

! -----------------------------------------------------------------------
! Rysowanie rezonansow dla ukladu pojedycznego tipu
! -----------------------------------------------------------------------
!!!
!open(unit=222,file="deltaClear.txt")
!xpos = nx/2*dx
!do deltax = 0 , 100.0 , 0.5
!    UTOTAL = 0
!
!    call system_dodaj_lorentza(la1,lx1,ly1,xpos+deltax,ypos)
!    call system_dodaj_lorentza(la2,lx2,ly2,xpos-deltax,ypos)
!    call system_rozwiaz_problem(1,TR_MAT)
!
!    TRANS_T = sum(TR_MAT(2,:))
!    TRANS_R = sum(TR_MAT(1,:))
!    write(222,"(20e20.8)"),deltax,TRANS_T,TRANS_R,TRANS_T+TRANS_R
!
!enddo
!
!close(222)
!stop
! -----------------------------------------------------------------
! Skany
! -----------------------------------------------------------------

open(unit=222,file="T.txt")
deltax = 0
xpos   = nx/2 * dx
x1     = 10.0
x2     = 50.0
call znajdz_maksimum(x1,x2,xpos,ypos)

UTOTAL = 0
call system_dodaj_lorentza(2.0D0,10.0D0,10.0D0,nx/2*dx-80,ny/2.*dx)
call system_dodaj_lorentza(-2.0D0,10.0D0,10.0D0,nx/2*dx+20,ny/3.*dx)
call system_dodaj_lorentza(3.0D0,10.0D0,10.0D0,nx/2*dx+100,3.*ny/4.*dx)

!call system_dodaj_lorentza(2.0D0,10.0D0,10.0D0,nx/2*dx-80,ny/2*dx)
!call system_dodaj_lorentza(2.0D0,10.0D0,10.0D0,nx/2*dx+80,ny/3*dx)

do xpos = 200 , nx*dx-200 , dx

    x1 = x1 - 5
    x2 = x2 + 5
    call znajdz_maksimum(x1,x2,xpos,ypos)
    !call system_rozwiaz_problem(1,TR_MAT)

    TRANS_T = sum(TR_MAT(2,:))
    TRANS_R = sum(TR_MAT(1,:))
    write(222,"(20e20.8)"),xpos,ypos,x1,TRANS_T,TRANS_R,TRANS_T+TRANS_R , UTOTAL(xpos/dx,ny/2)

enddo

close(222)


!UTOTAL = 0
!call system_dodaj_lorentza(2.0D0,10.0D0,10.0D0,nx/2*dx,ny/2*dx)
!call system_rozwiaz_problem(1,TR_MAT)
!call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
!call system_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
stop
!open(unit = 4356, file= "stabilizacja.txt" )
!do nx = 60 , 150
!
!!nx = 380/4 + 5
!call modjed_ustaw_konwersje_jednostek(M_EFF,E_MAT);
!call system_inicjalizacja(NX,NY,liczba_zrodel,DX);
!
!call zrodla(1)%zrodlo_ustaw(15,ny-15,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%zrodlo_ustaw(15,ny-15,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!
!call utworz_system()
!UTOTAL = 0
!call system_dodaj_pionowy_slupek_potencjalu((nx/2.0-20.)*dx,0.0D0,ny*dx,2.0D0,30.0D0,2.0d0)
!call system_dodaj_pionowy_slupek_potencjalu((nx/2.0+20.)*dx,0.0D0,ny*dx,2.0D0,30.0D0,2.0d0)
!
!call system_widmo(0.0D0,10.0D0,200,8,1,20)
!!call system_zapisz_widmo_do_pliku("widmo.txt",ZAPISZ_STANY_WLASNE)
!
!
!write(4356,"(200e20.6)"),nx*dx,Widmo_Evals(1:Widmo_NoStates)*1000.0*Rd
!call system_zwalnienie_pamieci()
!enddo
!close(4356)
!stop


!do atomic_Bz = 0 , 0.3 , 0.002

!UTOTAL = 0
!call modjed_ustaw_konwersje_jednostek(M_EFF,E_MAT);
!call system_inicjalizacja(NX,NY,liczba_zrodel,DX);
!call system_dodaj_pionowy_slupek_potencjalu((nx/2.0-20.)*dx,0.0D0,ny*dx,2.0D0,30.0D0,2.0d0)
!call system_dodaj_pionowy_slupek_potencjalu((nx/2.0+20.)*dx,0.0D0,ny*dx,2.0D0,30.0D0,2.0d0)
!
!
!open(unit=222,file="T.txt")
!do atomic_Ef = 1.0 , 10.0 , 0.005
!
!call zrodla(1)%zrodlo_ustaw(15,ny-15,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%zrodlo_ustaw(15,ny-15,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!call utworz_system()
!call system_rozwiaz_problem(1,TR_MAT)
!
!TRANS_T = sum(TR_MAT(2,:))
!TRANS_R = sum(TR_MAT(1,:))
!
!write(222,"(20e20.8)"),atomic_Ef,atomic_Bz,TRANS_T,TRANS_R
!
!enddo
!close(222)
!stop
!atomic_Ef = 2.66
!call zrodla(1)%zrodlo_ustaw(15,ny-15,1,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%zrodlo_ustaw(15,ny-15,nx,dx,atomic_Ef,atomic_Bz,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!call utworz_system()
!call system_rozwiaz_problem(1,TR_MAT)
!
!call system_widmo(0.0D0,10.0D0,200,8,1,20)
!call system_zapisz_widmo_do_pliku("widmo.txt",ZAPISZ_STANY_WLASNE)
!call system_zapisz_widmo_do_pliku("E.txt",ZAPISZ_WIDMO_VRTCAL)
!
!call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI);
!!call system_zapisz_do_pliku("flagi.txt",ZAPISZ_FLAGI);
!call system_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
!!call system_zapisz_do_pliku("J.txt",ZAPISZ_J);


call system_zwalnienie_pamieci()
if(allocated(TR_MAT))deallocate(TR_MAT)

contains

subroutine znajdz_maksimum(d1,d2,xpos,ypos)
    doubleprecision :: d1,d2,xpos,ypos
    doubleprecision :: dlambda,T1,Tpl,Tml,T2,d_x,deltax1,deltax2,Tmax,deltamax

    integer :: i
    dlambda = 0.01
    T1 = 0
    T2 = 0
    d_x = 0
    deltax1 = d1
    deltax2 = d2

    Tmax = 0
    deltamax = 0
    do d_x = deltax1 , deltax2 , abs(deltax1 - deltax2)/10

        call system_dodaj_lorentza(la1,lx1,ly1,xpos+d_x,ypos)
        call system_dodaj_lorentza(la2,lx2,ly2,xpos-d_x,ypos)
        call system_rozwiaz_problem(1,TR_MAT)
        call system_dodaj_lorentza(-la1,lx1,ly1,xpos+d_x,ypos)
        call system_dodaj_lorentza(-la2,lx2,ly2,xpos-d_x,ypos)
        Tml = sum(TR_MAT(2,:))
        if(Tmax < Tml) then
                deltamax = d_x
                Tmax = Tml
        endif
    enddo

    deltax1 = deltamax - abs(d1 - d2)/10
    deltax2 = deltamax + abs(d1 - d2)/10


    d_X = deltax1
    call system_dodaj_lorentza(la1,lx1,ly1,xpos+d_x,ypos)
    call system_dodaj_lorentza(la2,lx2,ly2,xpos-d_x,ypos)
    call system_rozwiaz_problem(1,TR_MAT)
    call system_dodaj_lorentza(-la1,lx1,ly1,xpos+d_x,ypos)
    call system_dodaj_lorentza(-la2,lx2,ly2,xpos-d_x,ypos)
    Tml = sum(TR_MAT(2,:))

    d_X = deltax2
    call system_dodaj_lorentza(la1,lx1,ly1,xpos+d_x,ypos)
    call system_dodaj_lorentza(la2,lx2,ly2,xpos-d_x,ypos)
    call system_rozwiaz_problem(1,TR_MAT)
    call system_dodaj_lorentza(-la1,lx1,ly1,xpos+d_x,ypos)
    call system_dodaj_lorentza(-la2,lx2,ly2,xpos-d_x,ypos)
    Tpl = sum(TR_MAT(2,:))


    do i = 1 , 50

    d_X = (deltax1 + deltax2)/2 - dlambda
    call system_dodaj_lorentza(la1,lx1,ly1,xpos+d_x,ypos)
    call system_dodaj_lorentza(la2,lx2,ly2,xpos-d_x,ypos)
    call system_rozwiaz_problem(1,TR_MAT)
    call system_dodaj_lorentza(-la1,lx1,ly1,xpos+d_x,ypos)
    call system_dodaj_lorentza(-la2,lx2,ly2,xpos-d_x,ypos)
    T1 = sum(TR_MAT(2,:))

    d_x = (deltax1 + deltax2)/2 + dlambda
    call system_dodaj_lorentza(la1,lx1,ly1,xpos+d_x,ypos)
    call system_dodaj_lorentza(la2,lx2,ly2,xpos-d_x,ypos)
    call system_rozwiaz_problem(1,TR_MAT)
    call system_dodaj_lorentza(-la1,lx1,ly1,xpos+d_x,ypos)
    call system_dodaj_lorentza(-la2,lx2,ly2,xpos-d_x,ypos)
    T2 = sum(TR_MAT(2,:))



    if( (T2-T1) > 0 ) then
        deltax1  = d_x - dlambda
        Tml      = T1
    else
        deltax2  = d_x - dlambda
        Tpl      = T1
    endif

    if( abs(deltax1-deltax2) < dlambda ) dlambda = abs(deltax1-deltax2)/10

    if( abs(deltax1-deltax2)/abs(deltax1) < 1.0e-4 )  exit

    !write(333,"(i,20e20.8)"),i,d_x,T1,abs(deltax1-deltax2)/abs(deltax1)
    enddo

    d1 = deltax1
    d2 = deltax1

end subroutine znajdz_maksimum


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

