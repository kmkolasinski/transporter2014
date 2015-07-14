
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
 double precision :: dx   = 2 , pdx , omega ,omega2 , x , y , gamma , xpos
 double precision,dimension(:,:), allocatable  :: TR_MAT
 integer :: i,j, sum_mod
 double precision :: width , G21(-1:1) , G23(-1:1) , sigmax , sigmay , kvec , poleB
 double precision :: fread(4)
 character(len=16) :: file_name
 doubleprecision :: lx0 , la0, ld0 , xtip , ytip
 doubleprecision, allocatable ,dimension(:,:) :: PPOT
 doubleprecision :: time1 , time2
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


call modjed_ustaw_konwersje_jednostek(0.0465D0,12.4D0);

!call modjed_ustaw_InGaAs()
!call modjed_ustaw_InSb()



UTOTAL = 0
allocate(PPOT(nx,ny))
PPOT = 0


!call system_inicjalizacja(NX,NY,liczba_zrodel);
!
!call zrodla(1)%zrodlo_ustaw(3,NY-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%zrodlo_ustaw(3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!call system_dodaj_lorentza(1.0D0,20.0D0,20.0D0,400.0D0,50.0D0);
!call utworz_system()
!call system_rozwiaz_problem(1,TR_MAT)
!
!call system_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
!call system_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
!call system_zapisz_do_pliku("kon.txt",ZAPISZ_KONTUR)
!UTOTAL = 0
!call system_dodaj_lorentza(1.0D0,20.0D0,20.0D0,400.0D0,60.0D0);
!!call system_rozwiaz_problem(1,TR_MAT)
!!call system_zapisz_do_pliku("d_phi.txt",ZAPISZ_PHI)
!call system_rozwiaz_problem_iteracyjnie(1,1,TR_MAT)
!!call system_rozwiaz_problem(1,TR_MAT)
!call system_zapisz_do_pliku("i_phi.txt",ZAPISZ_PHI)
!call system_zapisz_do_pliku("i_pot.txt",ZAPISZ_POTENCJAL)
!call system_zapisz_do_pliku("i_kon.txt",ZAPISZ_KONTUR)
!
!call system_zwalnienie_pamieci()
!if(allocated(TR_MAT))deallocate(TR_MAT)



call spinsystem_inicjalizacja(NX,NY,liczba_zrodel);

dx = atomic_DX
UTOTAL = 0
do i = 1 , nx
do j = 1 , ny
    x = i * dx
    y = j * dx
    UTOTAL(i,j) = gauss_gate(omega,xpos,0.0D0,sigmax,sigmay,x,y) + &
                  gauss_gate(omega,xpos,ny*dx,sigmax,sigmay,x,y)

enddo
enddo

!call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)

!UTOTAL = 0
!SUTOTAL = 0
!open(unit=432,file="dft_uint.txt")
!do i = 1 , nx
!do j = 1 , ny
!    read(432,"(4e20.6)"),fread
!    !print*,fread
!    SUTOTAL(i,j,+1) = fread(3) / 1000.0 / Rd
!    SUTOTAL(i,j,-1) = fread(4) / 1000.0 / Rd
!enddo
!    read(432,"(A)")
!enddo
!call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
!
!atomic_Ef = 57.8624513068335

call zrodla(1)%spinzrodlo_ustaw(3,NY-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
call zrodla(2)%spinzrodlo_ustaw(3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)
call utworz_system(nx)

lx0  = 5.0
la0  = 0.5
ld0  = 5.0 ! nm
!open(unit = 543, file= "T0Ef15.txt" )
!do xtip = 150 , 250 , 10
!do ytip = 0   , 160 , 10
!
!
!    do i = 1 , nx
!    do j = 1 , ny
!        PPOT(i,j) = gauss_gate(la0,xtip,ytip,ld0,ld0,i*dx,j*dx)
!    enddo
!    enddo
!
!    UTOTAL = UTOTAL + PPOT
!    call spinsystem_rozwiaz_problem(1,TR_MAT)
!    UTOTAL = UTOTAL - PPOT
!    write(543,"(20e20.6)"),xtip,ytip,sum(TR_MAT(2,:)),sum(TR_MAT(1,:))
!    print *, xtip,ytip,sum(TR_MAT(2,:)),sum(TR_MAT(1,:))
!
!    !call spinsystem_zapisz_do_pliku("phi0.txt",ZAPISZ_PHI)
!
!enddo
!    write(543,*),""
!enddo
!close(543)
time1 = get_clock()
call spinsystem_rozwiaz_problem(1,TR_MAT,.true.)
!call spinsystem_rozwiaz_problem(1,TR_MAT)
time1 = get_clock() - time1

time2 = get_clock()
PPOT = 0
call spinsystem_rozwiaz_problem_z_poprawki(1,TR_MAT,PPOT)
!call spinsystem_rozwiaz_problem(1,TR_MAT)
time2 = get_clock() - time2
print*,"========================================"
print*,"Czas wykonania pelnej symulacji:",time1
print*,"Czas wykonania sym. z  poprawki:",time2
print*,"========================================"
call spinsystem_zapisz_do_pliku("phi0.txt",ZAPISZ_PHI)

!open(unit = 543, file= "TpEf15II.txt" )
!do xtip = 150 , 250 , 10
!do ytip = 0   , 160 , 10
!
!    do i = 1 , nx
!    do j = 1 , ny
!        PPOT(i,j) = gauss_gate(la0,xtip,ytip,ld0,ld0,i*dx,j*dx)/Rd/1000.0
!    enddo
!    enddo
!    call spinsystem_rozwiaz_problem_z_poprawki(1,TR_MAT,PPOT)
!
!    write(543,"(20e20.6)"),xtip,ytip,sum(TR_MAT(2,:)),sum(TR_MAT(1,:))
!    print *, xtip,ytip,sum(TR_MAT(2,:)),sum(TR_MAT(1,:))
!
!    !call spinsystem_zapisz_do_pliku("phi0.txt",ZAPISZ_PHI)
!
!enddo
!    write(543,*),""
!enddo
!close(543)

!call spinsystem_rozwiaz_problem(1,TR_MAT,.true.)
!call spinsystem_rozwiaz_problem_z_poprawki(1,TR_MAT,PPOT)


call spinsystem_finalizuj_poprawki()




!call spinsystem_zapisz_do_pliku("phiP.txt",ZAPISZ_PHI)
!call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)

deallocate(PPOT)
call spinsystem_zwalnienie_pamieci()

!
!call solve_trans_system()

stop


call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
!call spinsystem_zapisz_do_pliku("kon.txt",ZAPISZ_KONTUR)
!call spinsystem_zapisz_do_pliku("flag.txt",ZAPISZ_FLAGI)

!TRANS_EIGPROBLEM_PERIODIC_X = .true.
!call spindft_initialize()
!call spindft_solve_temp_annealing()
!atomic_Ef = DFT_FINDED_EF
!call solve_trans_system()

!call zrodla(1)%spinzrodlo_ustaw(3,NY-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%spinzrodlo_ustaw(3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)



call spinsystem_rozwiaz_problem(1,TR_MAT)
call spinsystem_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)


!call spinsystem_zapisz_do_pliku("pot_dft.txt",ZAPISZ_POTENCJAL)
open(unit = 222, file= "T.txt" )
write(222,"(20e20.6)"),omega,atomic_Ef,atomic_Bz,sum(TR_MAT(2,:)),sum(TR_MAT(1,:)) !,G21(+1)-G21(-1),G21(+1),G21(-1)
close(222)

if(allocated(TR_MAT))deallocate(TR_MAT)


stop

!atomic_Ef = DFT_FINDED_EF
!call zrodla(1)%spinzrodlo_ustaw(3,NY-3,5,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!write(111,*),atomic_Ef,DFT_CURR_RESIDUUM,zrodla(1)%liczba_modow
!
!
!UTOTAL = 0
!do i = 1 , nx
!do j = 1 , ny
!    x = i * dx
!    y = j * dx
!    UTOTAL(i,j) = gauss_gate(omega,xpos,0.0D0,sigmax,sigmay,x,y) + &
!                  gauss_gate(omega,xpos,ny*dx,sigmax,sigmay,x,y)
!enddo
!enddo
!call spindft_free()
!call spindft_initialize()
!call spindft_fix_ef(atomic_Ef)
!call spindft_solve_temp_annealing()
!
!
!
!call spindft_free()
!atomic_Ef = DFT_FINDED_EF
!call zrodla(1)%spinzrodlo_ustaw(3,NY-3,5,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!write(111,*),atomic_Ef,DFT_CURR_RESIDUUM,zrodla(1)%liczba_modow
!!call zrodla(2)%spinzrodlo_ustaw(3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!!
!!call spinsystem_rozwiaz_problem(1,TR_MAT)
!!call spinsystem_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
!call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
!!call spinsystem_zapisz_do_pliku("kon.txt",ZAPISZ_KONTUR)
!
!!call spinsystem_widmo(0.0D0,2*atomic_Ef,150,0,8)
!!call spinsystem_widmo(0.0D0,2*atomic_Ef,150,2,8)
!!call spinsystem_zapisz_widmo_do_pliku("stany",ZAPISZ_STANY_WLASNE)
!!call spinsystem_zapisz_widmo_do_pliku("widmo.txt",ZAPISZ_WIDMO_VRTCAL)


call spinsystem_zwalnienie_pamieci()
if(allocated(TR_MAT))deallocate(TR_MAT)


!
!
!!open(unit = 222, file= "T.txt" )
!!open(unit = 333, file= "K.txt" )
!!!do  atomic_Bz = -0.0 , 3.0 , 0.05
!!
!UTOTAL = 0
!do i = 1 , nx
!do j = 1 , ny
!    x = i * dx
!    y = j * dx
!
!
!    UTOTAL(i,j) = gauss_gate(omega,xpos,0.0D0,sigmax,sigmay,x,y) + &
!                  gauss_gate(omega,xpos,ny*dx,sigmax,sigmay,x,y)
!    !UTOTAL(i,j) = ( 0.5*(omega/1000.0/Rd)**2*((y-ny*dx/2)**2)  )* &
!    !               exp( -0.5*( x - xpos )**2/(gamma)**2 ) + (omega2/1000.0/Rd)*(y-0*ny*dx/2)
!
!enddo
!enddo
!!
!!call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
!!call sledzenie_kf()
!!
!!stop
!
!!call zrodla(1)%spinzrodlo_ustaw(3,NY/2-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!!call zrodla(2)%spinzrodlo_ustaw(NY/2-3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!call reset_clock()
!
!call zrodla(1)%spinzrodlo_ustaw(3,NY-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!call zrodla(2)%spinzrodlo_ustaw(3,NY-3,nx,ZRODLO_KIERUNEK_LEWO,UTOTAL)
!
!
!
!
!    call utworz_system()
!    call spinsystem_rozwiaz_problem(1,TR_MAT)
!    G21(+1) = 0
!    G21(-1) = 0
!
!    do sum_mod  = 1  , zrodla(1)%liczba_modow -1 , 2
!        G21(+1) = G21(+1) + TR_MAT(2,sum_mod+0)
!        G21(-1) = G21(-1) + TR_MAT(2,sum_mod+1)
!    enddo
!!    print*,"W=",width,"R=",sum(TR_MAT(1,:)),"T=",sum(TR_MAT(2,:))," inne = " , G21(+1)+G21(-1)
!    write(222,"(20e20.6)"),omega,atomic_Bz,sum(TR_MAT(2,:)),sum(TR_MAT(1,:)) !,G21(+1)-G21(-1),G21(+1),G21(-1)
!!    write(333,"(30e20.6)"),omega,atomic_Bz,imag(qpc_zrodlo%ChiKvec(1:qpc_zrodlo%liczba_modow,+1))*L2LR
!
!
!    close(222)
!    close(333)
!call spinsystem_zapisz_do_pliku("phi.txt",ZAPISZ_PHI)
!call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL)
!call spinsystem_zapisz_do_pliku("j.txt",ZAPISZ_J_ALL)
!call spinsystem_zapisz_do_pliku("kon.txt",ZAPISZ_KONTUR)
!call spinsystem_zapisz_do_pliku("divj.txt",ZAPISZ_DIVJ)
!call spinsystem_zapisz_do_pliku("pol.txt",ZAPISZ_POLARYZACJE)
!
!
!call spinsystem_zwalnienie_pamieci()
!if(allocated(TR_MAT))deallocate(TR_MAT)



contains

subroutine solve_trans_system()

    doubleprecision, allocatable :: TRANS_SUTOTAL(:,:,:)
    integer :: X_REAPEAT , i ,j , r , nnx
    doubleprecision :: w_omega , ave_pot(-1:1)

    X_REAPEAT = 1
    nnx =nx+2*nx*X_REAPEAT
    allocate(TRANS_SUTOTAL(nnx,ny,-1:1))
    if(allocated(TR_MAT))deallocate(TR_MAT)

    TRANS_SUTOTAL(:,:,1) = SUTOTAL(1,ny/2,1)
    TRANS_SUTOTAL(:,:,-1) = SUTOTAL(1,ny/2,-1)

    ! licze sredni potencjal na wejsciu
    ave_pot = 0
    j = 0
    do i = ny/2  - 5 , ny/2 + 5
        j = j + 1
        ave_pot(+1) = ave_pot(+1) + SUTOTAL(1,i,+1)
        ave_pot(-1) = ave_pot(-1) + SUTOTAL(1,i,-1)
    enddo
    ave_pot = ave_pot / j



    do i = 1 , nnx
        w_omega = (i-1.0)/(nx-1.0)
        if( i > nnx - nx * X_REAPEAT) w_omega = -(i-nnx)/(nx-1.0)

        w_omega = 2*(w_omega - 0.5)

        if(w_omega > 1)  w_omega = 1.0
        if(w_omega < 0)  w_omega = 0.0

        w_omega = 1.0

        TRANS_SUTOTAL(i,:,+1) = w_omega*SUTOTAL(1,:,+1) + (1-w_omega) * ave_pot(+1)
        TRANS_SUTOTAL(i,:,-1) = w_omega*SUTOTAL(1,:,-1) + (1-w_omega) * ave_pot(-1)
    enddo

    do i = 1 , nx
    do j = 1 , ny
        TRANS_SUTOTAL(i+(X_REAPEAT)*nx,j,:) = SUTOTAL(i,j,:)
    enddo
    enddo

    call spinsystem_zwalnienie_pamieci()
    call spinsystem_inicjalizacja(nnx,NY,liczba_zrodel);


    UTOTAL  = 0
    SUTOTAL = TRANS_SUTOTAL
    do i = 1 , nnx
    do j = 1 , ny
        x = i * dx
        y = j * dx
        UTOTAL(i,j) = gauss_gate(omega,xpos + nx*atomic_DX*X_REAPEAT/2,0.0D0,sigmax,sigmay,x,y) + &
                      gauss_gate(omega,xpos + nx*atomic_DX*X_REAPEAT/2,ny*dx,sigmax,sigmay,x,y)
    enddo
    enddo


    call spinsystem_zapisz_do_pliku("pot.txt",ZAPISZ_POTENCJAL,nx,nnx-nx)

    call zrodla(1)%spinzrodlo_ustaw(3,NY-3,1  ,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
    call zrodla(2)%spinzrodlo_ustaw(3,NY-3,nnx,ZRODLO_KIERUNEK_LEWO,UTOTAL)


    call utworz_system(nnx)

    call spinsystem_rozwiaz_problem(1,TR_MAT)
    call spinsystem_zapisz_do_pliku("phi.txt",ZAPISZ_PHI,nx,nnx-nx)

    !call zrodla(1)%spinzrodlo_relacja_dyspersji(-0.5D0,0.5D0,0.01D0,atomic_Ef*3,"rel.txt")
!    call zrodla(1)%spinzrodlo_zapisz_mody("mup.txt","mdwn.txt")

    deallocate(TRANS_SUTOTAL)
    call spinsystem_zwalnienie_pamieci()

end subroutine solve_trans_system



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
    call spinsystem_inicjalizacja_ukladu(wjazd,4,4)
end subroutine utworz_system

double precision function gauss_gate(U0,xp,yp,sigmax,sigmay,x,y) result(rval)
        doubleprecision, intent(in):: U0,xp,yp,sigmax,sigmay,x,y

        rval =  U0 * exp( -(( x - xp)/(2*sigmax))**2 ) * exp( -(( y - yp)/(2*sigmay))**2 )

endfunction gauss_gate

!
!subroutine sledzenie_kf()
!double precision :: aux , ub , W1(50) , W2(50) , kvecs(50)
!integer :: no_kvecs , iter
!! przeliczamy jednostki
!call modjed_ustaw_konwersje_jednostek(0.0465D0,12.0D0);
!
!!so_rashba     = 0
!!so_loc        = 0
!atomic_Bz     = 0
!
!call qpc_zrodlo%spinzrodlo_ustaw(3,NY-3,nx/2,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
!print*,"g_land=",G_LAN
!
!
!no_kvecs = qpc_zrodlo%liczba_modow/2
!do iter = 1 , qpc_zrodlo%liczba_modow
!    kvecs(iter) = imag(qpc_zrodlo%ChiKvec(iter,+1)) * L2LR
!    print*,"kvec1=",kvecs(iter)
!enddo
!
!print*,"liczba sledzonych modow=",no_kvecs
!
!!kvec = 0
!!call qpc_zrodlo%spinzrodlo_relacja_dyspersji(-0.5D0,0.5D0,0.001D0,atomic_Ef*2,"rel1.txt")
!
!open(unit = 333, file= "ZeemanEodB.txt" )
!do atomic_Bz = 0.01 , 10.0 , 0.05
!
!    call spinmodpop_inicjalizacja(atomic_DX,qpc_zrodlo%N,atomic_Ef,atomic_Bz,qpc_zrodlo%Uvec)
!
!    W1 = 0
!    W2 = 0
!    do iter = 1 , no_kvecs
!        call spinmodpop_relacja_dyspersji(kvecs(2*iter-1),kvecs(2*iter-1) ,0.001D0,atomic_Ef*2,"rel.txt")
!        W1(iter) = (MODPOP_VALS(2*iter-1))*Rd*1000.0
!        call spinmodpop_relacja_dyspersji(kvecs(2*iter-0),kvecs(2*iter-0),0.001D0,atomic_Ef*2,"rel.txt")
!        W2(iter) = (MODPOP_VALS(2*iter-0))*Rd*1000.0
!    enddo
!
!
!    aux =  ((M_EFF)*BtoDonorB(atomic_Bz))*Rd*1000.0
!
!    write(333,"(100e20.6)"),atomic_Bz,omega,G_LAN,&
!                    sum((W1(1:no_kvecs)-W2(1:no_kvecs))/aux)/no_kvecs,&
!                        (W1(1:no_kvecs)-W2(1:no_kvecs))/aux,&
!                        W1(1:no_kvecs),W2(1:no_kvecs)
!
!
!
!    call spinmodpop_zwalnienie_pamieci()
!enddo
!
!close(333)
!
!call qpc_zrodlo%spinzrodlo_zwolnij_pamiec()
!end subroutine sledzenie_kf


end program transporter
