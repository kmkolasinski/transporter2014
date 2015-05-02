
program transporter
 use modutils
 use modspinsystem
 use modjed
 use modinip
 use modspinzrodlo
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
 integer :: i,j, sum_mod
 double precision :: width , G21(-1:1) , G23(-1:1) , sigmax , sigmay , kvec
 type(cspinzrodlo) :: qpc_zrodlo

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
call getDoubleValue("Dane","omega2",omega2)
call getDoubleValue("Dane","gamma",gamma)
call getDoubleValue("Dane","sigmax",sigmax)
call getDoubleValue("Dane","sigmay",sigmay)
call getDoubleValue("Dane","xpos",xpos)


call modjed_ustaw_konwersje_jednostek(0.0465D0,12.0D0);

!call modjed_ustaw_InGaAs()
!call modjed_ustaw_InSb()

print*,"solver=",TRANS_SOLVER
!TRANS_SOLVER = USE_UMFPACK
dx = atomic_dx
call spinsystem_inicjalizacja(NX,NY,liczba_zrodel);



!
!open(unit = 222, file= "T.txt" )
!open(unit = 333, file= "K.txt" )
!!do  atomic_Bz = -0.0 , 3.0 , 0.05
!
UTOTAL = 0
do i = 1 , nx
do j = 1 , ny
    x = i * dx
    y = j * dx
    !UTOTAL(i,j) = gauss_gate(omega,xpos,0.0D0,sigmax,sigmay,x,y) + &
    !              gauss_gate(omega,xpos,ny*dx,sigmax,sigmay,x,y)
    !UTOTAL(i,j) = ( 0.5*(omega/1000.0/Rd)**2*((y-ny*dx/2)**2)  )* &
    !               exp( -0.5*( x - xpos )**2/(gamma)**2 ) + (omega2/1000.0/Rd)*(y-0*ny*dx/2)

enddo
enddo
!
!call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
!call sledzenie_kf()
!
!stop

!call zrodla(1)%spinzrodlo_ustaw(3,NY/2-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%spinzrodlo_ustaw(NY/2-3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)
call reset_clock()

call zrodla(1)%spinzrodlo_ustaw(3,NY-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
call zrodla(2)%spinzrodlo_ustaw(3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)

!call qpc_zrodlo%spinzrodlo_ustaw(3,NY-3,nx/2,ZRODLO_KIERUNEK_PRAWO,UTOTAL)


!call zrodla(1) %spinzrodlo_relacja_dyspersji(-0.5D0,0.5D0,0.001D0,atomic_Ef*2,"rel1.txt")
!call qpc_zrodlo%spinzrodlo_relacja_dyspersji(-0.5D0,0.5D0,0.001D0,atomic_Ef*2,"rel1.txt")


    call utworz_system()
    call spinsystem_rozwiaz_problem(1,TR_MAT)
    G21(+1) = 0
    G21(-1) = 0

    do sum_mod  = 1  , zrodla(1)%liczba_modow -1 , 2
        G21(+1) = G21(+1) + TR_MAT(2,sum_mod+0)
        G21(-1) = G21(-1) + TR_MAT(2,sum_mod+1)
    enddo
!    print*,"W=",width,"R=",sum(TR_MAT(1,:)),"T=",sum(TR_MAT(2,:))," inne = " , G21(+1)+G21(-1)
    write(222,"(20e20.6)"),omega,atomic_Bz,sum(TR_MAT(2,:)),sum(TR_MAT(1,:)) !,G21(+1)-G21(-1),G21(+1),G21(-1)
!    write(333,"(30e20.6)"),omega,atomic_Bz,imag(qpc_zrodlo%ChiKvec(1:qpc_zrodlo%liczba_modow,+1))*L2LR



!call qpc_zrodlo%spinzrodlo_zwolnij_pamiec()
!enddo ! end of B loop

    close(222)
    close(333)
call spinsystem_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
call spinsystem_zapisz_do_pliku("j.txt",ZAPISZ_J_ALL)
call spinsystem_zapisz_do_pliku("kon.txt",ZAPISZ_KONTUR)
call spinsystem_zapisz_do_pliku("divj.txt",ZAPISZ_DIVJ)
call spinsystem_zapisz_do_pliku("pol.txt",ZAPISZ_POLARYZACJE)


call spinsystem_zwalnienie_pamieci()
if(allocated(TR_MAT))deallocate(TR_MAT)

print*,"TOTAL_TIME:" , get_clock()

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

double precision function gauss_gate(U0,xp,yp,sigmax,sigmay,x,y) result(rval)
        doubleprecision, intent(in):: U0,xp,yp,sigmax,sigmay,x,y

        rval =  U0 * exp( -(( x - xp)/(2*sigmax))**2 ) * exp( -(( y - yp)/(2*sigmay))**2 )

endfunction gauss_gate


subroutine sledzenie_kf()
double precision :: aux , ub , W1(50) , W2(50) , kvecs(50)
integer :: no_kvecs , iter
! przeliczamy jednostki
call modjed_ustaw_konwersje_jednostek(0.0465D0,12.0D0);

!so_rashba     = 0
!so_loc        = 0
atomic_Bz     = 0

call qpc_zrodlo%spinzrodlo_ustaw(3,NY-3,nx/2,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
print*,"g_land=",G_LAN


no_kvecs = qpc_zrodlo%liczba_modow/2
do iter = 1 , qpc_zrodlo%liczba_modow
    kvecs(iter) = imag(qpc_zrodlo%ChiKvec(iter,+1)) * L2LR
    print*,"kvec1=",kvecs(iter)
enddo

print*,"liczba sledzonych modow=",no_kvecs

!kvec = 0
!call qpc_zrodlo%spinzrodlo_relacja_dyspersji(-0.5D0,0.5D0,0.001D0,atomic_Ef*2,"rel1.txt")

open(unit = 333, file= "ZeemanEodB.txt" )
do atomic_Bz = 0.01 , 10.0 , 0.05

    call spinmodpop_inicjalizacja(atomic_DX,qpc_zrodlo%N,atomic_Ef,atomic_Bz,qpc_zrodlo%Uvec)

    W1 = 0
    W2 = 0
    do iter = 1 , no_kvecs
        call spinmodpop_relacja_dyspersji(kvecs(2*iter-1),kvecs(2*iter-1) ,0.001D0,atomic_Ef*2,"rel.txt")
        W1(iter) = (MODPOP_VALS(2*iter-1))*Rd*1000.0
        call spinmodpop_relacja_dyspersji(kvecs(2*iter-0),kvecs(2*iter-0),0.001D0,atomic_Ef*2,"rel.txt")
        W2(iter) = (MODPOP_VALS(2*iter-0))*Rd*1000.0
    enddo


    aux =  ((M_EFF)*BtoDonorB(atomic_Bz))*Rd*1000.0

    write(333,"(100e20.6)"),atomic_Bz,omega,G_LAN,&
                    sum((W1(1:no_kvecs)-W2(1:no_kvecs))/aux)/no_kvecs,&
                        (W1(1:no_kvecs)-W2(1:no_kvecs))/aux,&
                        W1(1:no_kvecs),W2(1:no_kvecs)



    call spinmodpop_zwalnienie_pamieci()
enddo

close(333)

call qpc_zrodlo%spinzrodlo_zwolnij_pamiec()
end subroutine sledzenie_kf


end program transporter
