
program transporter
 use modutils
 use modspinsystem
 use modjed
 use modinip
 use modspindft
! use modspinzrodlo
! use modsystem
 use ifport
 implicit none


! -----------------------------------------------
! Deklaracje zmiennych i tablic
! -----------------------------------------------
 integer :: nx            = 50;
 integer :: ny            = 50;
 integer :: liczba_zrodel = 2;
 double precision :: dx   = 2 , pdx , omega  , x , y , xpos
 double precision,dimension(:,:), allocatable  :: TR_MAT
 integer :: i,j, sum_mod
 double precision :: width  , sigmax , sigmay , kvec , poleB , dEy , widthT
 double precision :: fread(4) , lx , ly, utip , lxmin ,lxmax , lymax , lymin
 character(len=16) :: file_name

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
call getDoubleValue("Dane","Bx",atomic_Bx)
call getDoubleValue("Dane","By",atomic_By)


call getDoubleValue("Dane","so_alpha3D",so_alpha3D)
call getDoubleValue("Dane","so_Fz",so_Fz)
call getDoubleValue("Dane","omega",omega)
call getDoubleValue("Dane","sigmax",sigmax)
call getDoubleValue("Dane","sigmay",sigmay)
call getDoubleValue("Dane","xpos",xpos)
call getDoubleValue("Dane","dEy",dEy)


call modjed_ustaw_konwersje_jednostek(0.0465D0,12.4D0);

!call modjed_ustaw_InGaAs()
!call modjed_ustaw_InSb()

open(unit = 22211, file= "T.txt" )

call spinsystem_inicjalizacja(NX,NY,liczba_zrodel);


dx      = atomic_DX
Utip    = 30.0
poleB	= atomic_Bz
Utip    = 5.0
do so_alpha3D = 0.0 , 1.0 , 0.02

    atomic_Rashba = so_alpha3D*so_Fz
    atomic_LOC    = so_alpha3D
    so_rashba     = atomic_Rashba * L2LR  / 1000.0 / Rd
    so_loc        = atomic_LOC    * L2LR * L2LR

    call zrodla(1)%spinzrodlo_ustaw(3,NY-3,1 ,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
    call zrodla(2)%spinzrodlo_ustaw(3,NY-3,nx,ZRODLO_KIERUNEK_LEWO ,UTOTAL)
    call utworz_system(nx)
    UTOTAL  =    0.0
    lymin   =   30.0
    lymax   =  250.0
    widthT  =   40.0
    lx      = nx/2.0*atomic_DX
    ly      = ny/2.0*atomic_DX

    do i = nx/10 , nx - nx/10
    do j = 1 , ny
        x = i * dx
        y = j * dx
        UTOTAL(i,j) = gauss_gate(omega,xpos,0.0D0,sigmax,sigmay,x,y) + &
                      gauss_gate(omega,xpos,ny*dx,sigmax,sigmay,x,y) + j*dx * dEy


        UTOTAL(i,j) = UTOTAL(i,j)  + Utip/( 1 + ((x-lx)/widthT)**2 + ((y-ly)/widthT)**2 )

    enddo
    enddo

    call spinsystem_rozwiaz_problem(1,TR_MAT)
    !call spinsystem_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
    write(22211,"(20e20.6)"),poleB,so_alpha3D,sum(TR_MAT(2,:)),sum(TR_MAT(1,:)),zrodla(1)%liczba_modow+0.0

enddo
close(22211)

!open(unit=222,file="T.txt")
!
!Utip  =   5.0
!lymin = -80.0
!lymax = 200.0
!widthT=  80.0
!lx    = 200.0
!
!!do lx = lxmin , lxmax , dx
!do ly = lymin , lymax , dx*2
!
!
!    call spinsystem_dodaj_lorentza(+Utip,widthT,widthT,lx,ly)
!    call utworz_system(nx)
!    call spinsystem_rozwiaz_problem(1,TR_MAT)
!    TRANS_T = sum(TR_MAT(2,:))
!    TRANS_R = sum(TR_MAT(1,:))
!    write(222,"(20e20.8)"),lx,ly,TRANS_T,TRANS_R,zrodla(1)%liczba_modow - TRANS_T - TRANS_R
!    call spinsystem_dodaj_lorentza(-Utip,widthT,widthT,lx,ly)
!    call spinsystem_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
!
!enddo
!!write(222,*),""
!!enddo
!close(222)




!write(22211,"(20e20.6)"),omega,atomic_Ef,poleB,sum(TR_MAT(2,:)),sum(TR_MAT(1,:))


!close(22211)
!call spinsystem_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
!call spinsystem_zapisz_do_pliku("j.txt",ZAPISZ_J_ALL)
!call spinsystem_zapisz_do_pliku("tot.txt",ZAPISZ_POTENCJAL)
if(allocated(TR_MAT))deallocate(TR_MAT)
call spindft_free()


contains




subroutine utworz_system(nx)
    integer :: nx
    integer :: i,j , wjazd , promien
    ! prosty test


    wjazd  = nx/2+3
    GFLAGS = B_EMPTY

    ! dajemy na brzegach dirichleta
    GFLAGS(1,:)  = B_DIRICHLET
    GFLAGS(NX,:) = B_DIRICHLET
    GFLAGS(:,1)  = B_DIRICHLET
    GFLAGS(:,NY) = B_DIRICHLET
!    call system_inicjalizacja_ukladu(wjazd,2,4)
    call spinsystem_inicjalizacja_ukladu(wjazd,4,0)
end subroutine utworz_system

double precision function gauss_gate(U0,xp,yp,sigmax,sigmay,x,y) result(rval)
        doubleprecision, intent(in):: U0,xp,yp,sigmax,sigmay,x,y

        rval =  U0 * exp( -(( x - xp)/(2*sigmax))**2 ) * exp( -(( y - yp)/(2*sigmay))**2 )

endfunction gauss_gate




end program transporter
