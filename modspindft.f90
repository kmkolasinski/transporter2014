module modspindft
    use modspinsystem
    use modjed
    use modinip
    use modutils
    use modmixers
    use xc
    use ifport
    implicit none
    private
    doubleprecision, allocatable :: u_exchange(:,:)         ! tablica przechowujaca potencjal wymiany (:,-1) - spin down , (:,+1) - spin up
    doubleprecision, allocatable :: u_correlation(:,:)      ! jak wyzej dla potencjalu korelacyjnego
    doubleprecision, allocatable :: u_hartree(:)            ! potencjal Hartree of gazu elektronowego
    doubleprecision, allocatable :: u_donor_hartree(:)      ! potencjal od donorow - ten jest staly liczony jest raz przed symulacja
    doubleprecision, allocatable :: rho(:,:)                ! gestosc w nastepnej iteracji
    doubleprecision, allocatable :: new_rho(:,:)            ! gestosc w obecnej iteracji
    doubleprecision, allocatable :: rho_tot(:)              ! calkowita gestosc
    doubleprecision, allocatable :: xc_res(:,:)             ! tablica pomocnicza wykorzystywana przez modul XC
    doubleprecision, allocatable :: VSQRTL(:,:,:)           ! tablica z przeliczonymi odleglosciami jak w calce kulombowskiej
    integer,dimension(:,:),allocatable   :: DFTINDEX        ! tablica zawierajaca unikalne indeksy na siatce. Odwolania postawi rho(GINDEX(i,j),-1)

    integer          :: DFT_NO_INTEGRAL_REPEAT_X ! w przypadku periodycznych warunkow brzegowych
    integer          :: DFT_TRANSMAX             ! rozmiar macierzy - liczba wezlow w ukladzie razy spin
    integer          :: NX , NY                  ! podzial ukladu szerokosc wysokosc w oczkach
    integer          :: DFT_NO_STATES            ! poczatkowa oraz pozniej obecna liczba stanow szukanych przez feasta
    integer          :: GLOBAL_ITER              ! liczba iteracji dla stalej temperatury
    integer          :: GLOBAL_TEMP_ITER         ! calkowita liczba iteracji
    integer          :: DFT_MAX_ITER            ! maksymalna liczba iteracji
    integer          :: DFT_IMPROVE_FEAST_STEPS ! co ile beda zmieniane ustawienia feasta tak aby liczona liczba stanow byla zawsze jak najmniejsza
    integer          :: DFT_TEMP_NO_STEPS       ! liczba iteracji schladzania
    integer          :: DFT_USE_CORE_FREEZING   ! jesli 1 to gdy symulacja osiagnie pewien stopien uzbieznienia czesc stanow zostaje zamrozonych
    integer          :: DFT_USE_INITIAL_GUESS   ! jesli 1 to jako start zostanie uzyta zgadywana gestosc elektronowa, zwykle dziala i zbiega sie szybciej
    logical          :: DFT_FIX_CORE_STATES     ! jesli residuum bedzie male do zamrozimy
    integer          :: DFT_MEM_STEPS           ! liczba krokow wstecz pamietane przez metode do mieszania gestosci


    double precision :: DFT_W_PARAM , DFT_MAX_W_PARAM
    double precision :: DFT_AVERAGE_DELTA_ENERGY ! srednia odleglosc miedzy stanami
    double precision :: DFT_CURR_DELTA_ENERGY    ! obecna zmiana w energii Fermiego
    double precision :: DFT_NO_DONORS            ! liczba donorow w ukladzie
    double precision :: DFT_Z_SPACER             ! odleglosc warstwy donorow od gazu elektronowego w [nm]
    double precision :: DFT_TEMP                 ! maksymalna wartosc temperatury w [K]
    double precision :: DFT_TEMP_MIN             ! minimalna wartosc temperatury w [K]
    double precision :: DFT_CURR_TEMP            ! obecna temperatura w [K]
    double precision :: DFT_RESIDUUM             ! maksymalne residuum ponizej ktorego symulacja zostaje uznana za uzbiezniona
    double precision :: DFT_CURR_RESIDUUM        ! obecne residuum
    doubleprecision  :: DFT_FINDED_EF            ! znaleziona wartosc energii Fermiego po uzbieznieniu w [meV]

    ! inne zmienne pomocnicze
    doubleprecision  :: Ef , old_Ef , Bz , Bx , By , DX , DFT_ATOMIC_EF
    doubleprecision  :: max_stanow , max_Ef , dft_donor_dens

    ! klasa odpowiedzialna za mieszanie gestosci z iteracji na iteracje
    type(scfmixer) :: mixer

    ! robienie zmiennych publiczymi
    public :: spindft_initialize , spindft_free , spindft_solve
    public :: spindft_solve_temp_annealing
    public :: spindft_zapisz_rho , spindft_wczytaj_rho_poczatkowe
    public :: DFT_TEMP ,  DFT_TEMP_MIN  , DFT_MAX_ITER , DFT_CURR_RESIDUUM , DFT_FINDED_EF
    public :: DFT_NO_DONORS , DFTINDEX , new_rho
    contains


! ======================================================================== !
! =====================     Inicjalizacja problemu  ====================== !
!
! Funkcja inicjalizuje tablice oraz przelicza wstepne dane na podstawie
! informacji o systemie. Aby ta funkcja wywolala sie poprawnie w pierwszej
! kolejnosci musza zostac wywolane funkcje inicjalizujace system do
! obliczen transportu elektronowego. Np.
! call spinsystem_inicjalizacja(NX,NY,liczba_zrodel);
! call zrodla(1)%spinzrodlo_ustaw(4,NY-3,1,ZRODLO_KIERUNEK_PRAWO,UTOTAL)
! call utworz_system(nx)
! call spindft_initialize()
! ...
! ======================================================================== !
    subroutine spindft_initialize()
        integer ::  i , j , iter , mx


        DFT_FIX_CORE_STATES   = .false.
        DFT_FIX_CORE_STATES   = .false.
        DFT_USE_INITIAL_GUESS = .false.
        DFT_USE_CORE_FREEZING = 0

        BZ  = BtoDonorB(atomic_Bz)
        Bx  = BtoDonorB(atomic_Bx)
        By  = BtoDonorB(atomic_By)
        DX  = atomic_DX*L2LR
        nx = size(GINDEX,1)
        ny = size(GINDEX,2)


        print*,"Inicjalizacja problemu SPIN-DFT:"

        call getDoubleValue("DFT","temp_start",DFT_TEMP)
        call getDoubleValue("DFT","temp_min",DFT_TEMP_MIN)

        call getDoubleValue("DFT","z_spacer",DFT_Z_SPACER)
        call getDoubleValue("DFT","no_donors",DFT_NO_DONORS)
        call getDoubleValue("DFT","w_min",DFT_W_PARAM)
        call getDoubleValue("DFT","w_max",DFT_MAX_W_PARAM)
        call getDoubleValue("DFT","init_max_energy",DFT_ATOMIC_EF)
        call getDoubleValue("DFT","residuum",DFT_RESIDUUM)


        call getIntValue   ("DFT","temp_no_steps",DFT_TEMP_NO_STEPS)
        call getIntValue   ("DFT","max_iter",DFT_MAX_ITER)
        call getIntValue   ("DFT","mem_steps",DFT_MEM_STEPS)
        call getIntValue   ("DFT","init_states",DFT_NO_STATES)
        call getIntValue   ("DFT","improve_feast",DFT_IMPROVE_FEAST_STEPS)
        call getIntValue   ("DFT","use_core_freezing",DFT_USE_CORE_FREEZING)
        call getIntValue   ("DFT","use_initial_guess",DFT_USE_INITIAL_GUESS)

        Ef = DFT_ATOMIC_EF/1000.0/Rd

        ! czyszczenie tablic przed ewentualna ponowna inicjalizacja
        call spindft_free()

        ! przeliczanie tablicy indeksow dla problemu wlasnego
        allocate(DFTINDEX(nx,ny))

        DFTINDEX = 0
        iter   = 1
        do i = 1 , nx
        do j = 1 , ny
            ! w problemie wlasnym iteresujace sa tylko punkty wewnatrz ukladu flaga B_NORMAL
            if( GFLAGS(i,j) == B_NORMAL  ) then
                DFTINDEX(i,j) = iter
                iter = iter + 1
           endif
        enddo
        enddo

        ! gdy jestesmy w trybie periodycznosci to dodatkowo dochodzi mapowanie periodycznosci
        ! w kierunku X
        if(TRANS_EIGPROBLEM_PERIODIC_X) then
            do j = 1 , ny
                DFTINDEX(nx,j) = DFTINDEX(2   ,j)
                DFTINDEX(1 ,j) = DFTINDEX(nx-1,j)
            enddo
        endif

        DFT_TRANSMAX = iter - 1 ! ilosc oczek na pojedynczy spin

        ! obliczanie gestosci donorow na podstwie zadanej liczby elektronow i rozmiaru ukladu
        ! Zalozenie: donory aktywuja sie tam gdzie gaz elektronowy
        dft_donor_dens = DFT_NO_DONORS / DFT_TRANSMAX / DX / DX


        allocate(u_hartree      (DFT_TRANSMAX))
        allocate(u_donor_hartree(DFT_TRANSMAX))
        allocate(u_exchange     (DFT_TRANSMAX,-1:1))
        allocate(u_correlation  (DFT_TRANSMAX,-1:1))
        allocate(rho            (DFT_TRANSMAX,-1:1))
        allocate(new_rho        (DFT_TRANSMAX,-1:1))
        allocate(rho_tot        (DFT_TRANSMAX))
        allocate(xc_res         (DFT_TRANSMAX,3))



        ! jesli jestemy w trybie periodycznosci to calkowanie po gestosci
        ! bedzie odbywac sie po pudlach periodycznosci w kierynku X.
        ! Zalozenie: na sztywno ustawione 8 pudel periodycznosci
        if(TRANS_EIGPROBLEM_PERIODIC_X) then
            DFT_NO_INTEGRAL_REPEAT_X = 8
            allocate(VSQRTL(0:nx*DFT_NO_INTEGRAL_REPEAT_X,0:ny,0:1)) ! dwie warstwy 0 (na gaz) i 1 (na donory)
        else
            allocate(VSQRTL(0:nx,0:ny,0:1)) ! dwie warstwy 0 (na gaz) i 1 (na donory)
            DFT_NO_INTEGRAL_REPEAT_X = 1
        endif


        ! obliczanie pierwiastkow do calku kulomboskiej
        do i = 0 , nx * DFT_NO_INTEGRAL_REPEAT_X
        do j = 0 , ny
            ! w plaszczyznie gazu
            if( i == 0 .and. j == 0 ) then
                VSQRTL(0,0,0) = (4*DX*log(sqrt(2.0)+1))/DX/DX
            else
                VSQRTL(i,j,0) = 1.0/sqrt( (i*DX)**2 + (j*DX)**2  )
            endif
            ! do donorow
            VSQRTL(i,j,1) = 1.0/sqrt( (i*DX)**2 + (j*DX)**2 + (DFT_Z_SPACER*L2LR)**2 )

        enddo
        enddo
        call calcHartreePotential(bCalcHartreeDonor=.true.)


        GLOBAL_ITER = 1
        rho_tot     = 0
        rho         = 0
        new_rho     = 0

        DFT_FINDED_EF = 0


        ! Inicjalizacja procedury miksujacej
        call mixer%init(SCF_MIXER_EXTENDED_ANDERSON,DFT_TRANSMAX*2,DFT_W_PARAM,DFT_MAX_W_PARAM,DFT_MEM_STEPS)


    end subroutine spindft_initialize

! ======================================================================== !
! Pseudonaukowa metoda szacowania gestosci poczatkowej, oparta na
! wlasnych obserwacjach oraz na teorii LDA i promienu Seitza.
! ======================================================================== !
    subroutine spindft_calculate_approximated_initial_dens()

        integer :: i,j,radius_step,xi,yi,ave_counter,xii,k,smoothstep
        doubleprecision :: ave_dens,ave_radius,ave_rho,dval,dvalx,dvaly,min_dist
        double precision, allocatable :: SEITZ_DIST(:,:,:)

        print*," <c> Obliczanie szacowanej gestosci elektronowej na start."
        ave_dens    = DFT_NO_DONORS / DFT_TRANSMAX / (DX)**2 ! Srednia gestosc ukladu
        ave_radius  = 1.0/sqrt( M_PI * ave_dens)             ! Promien seitza dla tej gestosci
        radius_step = ave_radius/DX                          ! Przliczenie promienia na liczbe oczek siatki


        allocate(SEITZ_DIST(NX,NY,2))
        SEITZ_DIST = 0
        ! robienie mapy odleglosci od krawedzi i przeliczanie jej
        ! jednostkacj r/rs - gdzie rs to primien Seizta dla 2D
        do i = 1 , NX
        do j = 1 , NY
            min_dist = max(nx*dx,ny*dx)

            do xi = 1 , NX
            do yi = 1 , NY
                if(TRANS_EIGPROBLEM_PERIODIC_X) then
                if(GFLAGS(xi,yi) == B_EMPTY) then
                    dvalx = sqrt(abs(xi-i+0.0)**2 + abs(yi-j+0.0)**2)*dx-2*dx
                    if(dvalx < min_dist) min_dist = dvalx
                endif
                else
                if(GFLAGS(xi,yi) /= B_NORMAL) then
                    dvalx = sqrt(abs(xi-i+0.0)**2 + abs(yi-j+0.0)**2)*dx-2*dx
                    if(dvalx < min_dist) min_dist = dvalx
                endif
                endif
            enddo
            enddo
            SEITZ_DIST(i,j,1) = max(min_dist/ave_radius,0.0)
        enddo
        enddo

        ! male wygladzenie otrzymanej mapy w celu otrzymania lepszych
        ! przyjemniejszych dla oka rezultatow
        smoothstep = max(1,radius_step/4)
        do k = 1 , 1 ! liczba wygladzen - tutaj tylko jeden ale moze byc dowolna
        do i = 1 , NX
        do j = 1 , NY
            ave_rho     = 0
            ave_counter = 0
            do xi = i-smoothstep , i+smoothstep
            do yi = j-smoothstep , j+smoothstep
                xii = xi
                if(TRANS_EIGPROBLEM_PERIODIC_X) then
                    xii = abs(mod(xi-1,NX))+1
                endif
                if(xii>0 .and. xii<=NX .and. yi>0 .and. yi<=NY) then
                        ave_rho     = ave_rho     + SEITZ_DIST(xii,yi,1)
                        ave_counter = ave_counter + 1
                endif
            enddo
            enddo
            SEITZ_DIST(i,j,2) = ave_rho/ave_counter
        enddo
        enddo
            SEITZ_DIST(:,:,1) = SEITZ_DIST(:,:,2)
        enddo

        ! przeliczanie gestosci elektronowej na podstawie wyssanego z palca
        ! wzoru.
        rho_tot = 0
        do i = 1 , NX
        do j = 1 , NY
            if(DFTINDEX(i,j) > 0) then
                dvaly                   = SEITZ_DIST(i,j,1)
                ! pik w okolicy r/rs=1 o szerokosci 1 potem sie wyplaszcza do wartosci arbitralnie dobranej
                dval                    = 1.4*dvaly**2  * exp(-((dvaly)**2)) + 1.0/(1+exp(-(dvaly)**2)) - 0.5
                rho_tot(DFTINDEX(i,j))  =  dval
            endif
        enddo
        enddo

        ! normalizacja gestosci do liczby donorow
        rho_tot = rho_tot / sum(rho_tot) / dx**2 * DFT_NO_DONORS
        ! dzielenie gestosci na dwa
        rho(:,-1)   = rho_tot/2
        rho(:,+1)   = rho_tot/2
        new_rho     = rho


        open(unit=6789,file="init_dens.txt")
        do i = 1 , NX
        do j = 1 , NY
            if( DFTINDEX(i,j) > 0) then
                write(6789,*),i*atomic_DX,j*atomic_DX,rho_tot(DFTINDEX(i,j))
            else
                write(6789,*),i*atomic_DX,j*atomic_DX,0.0
            endif
        enddo
            write(6789,*),""
        enddo
        close(6789)

        deallocate(SEITZ_DIST)


    end subroutine spindft_calculate_approximated_initial_dens

! ======================================================================== !
! =================     Start problemu z symulowanym      =================!
! =================             wyzarzaniem               =================!
!
!
!
! ======================================================================== !
    subroutine spindft_solve_temp_annealing()
        integer :: itemp , icopy
        doubleprecision :: tomega,dtmp

        ! jesli tak wybrano to gestosc poczatkowa jest szacowana
        if(DFT_USE_INITIAL_GUESS) then
            call spindft_calculate_approximated_initial_dens()
        endif

        GLOBAL_TEMP_ITER = 0
        open(unit = 8998, file= "dft_iter.txt" )


        do itemp = 1 , DFT_TEMP_NO_STEPS
            tomega = (itemp-1.0)/(DFT_TEMP_NO_STEPS-1.0)

            if(DFT_TEMP_NO_STEPS == 1) tomega = 0.0

            tomega        = 1  -  (1-tomega)**2
            DFT_CURR_TEMP = DFT_TEMP*(1 - tomega) + tomega*DFT_TEMP_MIN

            call spindft_solve()

            print*,"Symulacja uzbiezniona, po krokach:",GLOBAL_ITER
            print*,"Zakonczony krok czasowy          :",DFT_TEMP_NO_STEPS
            print*,"Z residuum rownym                :",DFT_CURR_RESIDUUM
            print*,"Z temperatura                    :",DFT_CURR_TEMP
            DFT_FINDED_EF = Ef*1000*Rd
            print*,"Znaleziona energia Fermiego      :",DFT_FINDED_EF
            print*,"Liczba elektronow up             :",sum(new_rho(:,+1))*DX*DX
            print*,"Liczba elektronow up             :",sum(new_rho(:,-1))*DX*DX

        enddo ! end of do no temp steps
        print*,"-----------------------------------------------------"
        print*,"Symulacja uzbiezniona w sumie po krokrach:",GLOBAL_TEMP_ITER
        print*,"Z residuum rownym                        :",DFT_CURR_RESIDUUM
        DFT_FINDED_EF = Ef*1000*Rd
        print*,"Znaleziona energia Fermiego              :",DFT_FINDED_EF

    end subroutine spindft_solve_temp_annealing


! ======================================================================== !
! =============  Start problemu dla ustalonej wartosci T    ===============!
!
! ======================================================================== !
    subroutine spindft_solve()
        integer          :: i,j,m
        double precision :: bz_copy
        ! obliczenia startowe dla gestosci elektronowej takiej samej
        ! jak gestosc donorow
        GLOBAL_ITER = 1
        kbT       = 8.617D-5/Rd*DFT_CURR_TEMP
        old_Ef    = 100.0
        bz_copy   = atomic_Bz


        ! Pierwsza iteracja:
        call spindft_ks_iteration(0)
        call calcComplexEleDensity(reset=.true.)


        DFT_FINDED_EF = 0
        do GLOBAL_ITER = 2 , DFT_MAX_ITER
            GLOBAL_TEMP_ITER  = GLOBAL_TEMP_ITER + 1
            print*,"**************************************************"
            print*,"    Iteracja DFT:", GLOBAL_ITER
            print*,"    Temperatura :", kbT * Rd / 8.617D-5
            print*,"**************************************************"
            ! polaryzacja stanow w pierwszych iteracjach
            if(GLOBAL_TEMP_ITER < 3) then
               atomic_Bz = 0.5
            else
               atomic_Bz = bz_copy
            endif

            ! W pewnych warunkach przeliczane sa liczba stanow i maksymalna
            ! energia dla feasta.
            if( ( MOD(GLOBAL_ITER,DFT_IMPROVE_FEAST_STEPS) == 0 &
                  .and. DFT_CURR_DELTA_ENERGY < DFT_AVERAGE_DELTA_ENERGY ) &
                  .and. GLOBAL_ITER       > 5 &
                  .and. DFT_FIX_CORE_STATES == .false. ) then

                print*,"Zmiana parametrow feast-a:"
                print*,"    Ef z     :",DFT_ATOMIC_EF," na:",max_Ef * Rd * 1500.0
                print*,"    L. stanow:",DFT_NO_STATES," na:",max_stanow  + 10
                DFT_ATOMIC_EF = max_Ef * Rd * 1000.0
                DFT_NO_STATES = max_stanow  + 10
                call spindft_ks_iteration(0)

            else
                call spindft_ks_iteration(1)
            endif


            ! Obliczanie gestosci na podstawie znalezionych stanow K-S
            call calcComplexEleDensity()


            ! Debugowanie:
            write(8998,"(20e20.8)"),GLOBAL_TEMP_ITER+0.0D0,DFT_CURR_TEMP,GLOBAL_ITER+0.0,Ef*1000*Rd,DFT_CURR_RESIDUUM

            open(unit=222,file="dft_rho.txt")
            do i = 1 , nx
            do j = 1 , ny
                if(DFTINDEX(i,j) > 0 ) &
                write(222,"(20f20.6)"),i*atomic_DX,j*atomic_DX,new_rho(DFTINDEX(i,j),1),new_rho(DFTINDEX(i,j),-1)
            enddo
                write(222,*),""
            enddo
            close(222)

!            open(unit=555,file="dft_uint.txt")
!            do i = 1 , nx
!            do j = 1 , ny
!                if(DFTINDEX(i,j) > 0 ) then
!                    write(555,"(5f20.6)"),i*atomic_DX,j*atomic_DX,SUTOTAL(i,j,1)*Rd*1000.0,SUTOTAL(i,j,-1)*Rd*1000.0
!                else
!                    write(555,"(5f20.6)"),i*atomic_DX,j*atomic_DX,0.0D0,0.0D0
!                endif
!            enddo
!                write(555,*),""
!            enddo
!            close(555)

            if(  DFT_CURR_RESIDUUM < DFT_RESIDUUM .and. GLOBAL_ITER > 5 ) then
                print*,"Finished..."
                exit
            endif

        enddo

        ! zwalnianie pieniedzy
        call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,2,8)

    end subroutine spindft_solve


! ======================================================================== !
! =============         Zwalnianie pamieci                  ===============!
!
! ======================================================================== !
    subroutine spindft_free()
        print*,"SPIN-DFT free memory"
        if(allocated(DFTINDEX))   deallocate(DFTINDEX)
        if(allocated(u_hartree))  deallocate(u_hartree)
        if(allocated(u_donor_hartree))  deallocate(u_donor_hartree)
        if(allocated(u_exchange))  deallocate(u_exchange)
        if(allocated(u_correlation)) deallocate(u_correlation)
        if(allocated(rho))         deallocate(rho)
        if(allocated(new_rho))     deallocate(new_rho)
        if(allocated(rho_tot))     deallocate(rho_tot)
        if(allocated(xc_res))      deallocate(xc_res)
        if(allocated(VSQRTL))      deallocate(VSQRTL)

        call mixer%free_mixer()
        ! Zwalnianie pamieci
        call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,2,8)
    end subroutine spindft_free



! ======================================================================== !
! =============                Iteracja KS                  ===============!
! 1. Obliczanie potencjalow: wymiany, korelacyjny i kulombowski
! ======================================================================== !
    subroutine spindft_ks_iteration(opt)
        integer :: opt
        integer :: i,j,s
        logical :: CZY_OSTATNIO_ZAMROZONO_CORE_STATES

        ! Zmienne wykorzystywywane przez algorytm do mrozenia stanow rdzeniowych
        integer                     :: no_core_states , Core_Iter , Core_Overlap_Test , propper_num_of_states
        complex*16,allocatable      :: Core_Widmo_Vecs(:,:) , Core_Widmo_Vecs_Tmp(:,:)
        doubleprecision,allocatable :: Core_Widmo_Evals(:)  , Core_Widmo_Evals_Tmp(:)
        doubleprecision             :: Core_Last_Enery


        call XAttaccalite_1(new_rho(:,1),new_rho(:,-1),xc_res)
        u_exchange(:,+1) = xc_res(:,1)
        u_exchange(:,-1) = xc_res(:,2)

        call CAttaccalite_1(new_rho(:,1),new_rho(:,-1),xc_res)
        u_correlation(:,+1) = xc_res(:,1)
        u_correlation(:,-1) = xc_res(:,2)

        call calcHartreePotential()

        ! sumowanie wszystkich przyczynkow
        do s = 1 , -1  , -2
        do i = 1 , nx
        do j = 1 , ny
            if(DFTINDEX(i,j) > 0 ) then
                SUTOTAL(i,j,s) =   u_exchange(DFTINDEX(i,j),s) &
                                 + u_correlation(DFTINDEX(i,j),s) &
                                 - u_hartree(DFTINDEX(i,j))
            endif
        enddo
        enddo
        enddo

        ! uzywaj zamrazania stanow rdzeniowych/walencyjnych
        if(DFT_USE_CORE_FREEZING == 1) then
        print*," <c> Zamrazanie stanow walencyjnych wlaczone!"

        ! Zamrazanie stanow rdzeniowych
        CZY_OSTATNIO_ZAMROZONO_CORE_STATES = DFT_FIX_CORE_STATES
        DFT_FIX_CORE_STATES = .false.
        Core_Last_Enery     = 0.0
        no_core_states      = 0
        ! jesli obecna zmiana energii jest mniejsza niz srednia odleglosc miedzy stanami
        ! to mrozenie jest mozliwe. Warunek konieczny!
        if(  DFT_CURR_DELTA_ENERGY < DFT_AVERAGE_DELTA_ENERGY/100 .and. GLOBAL_ITER > 5 ) then
            print*,"Warning freezing of core states is enabled!: Curr Res:",DFT_CURR_RESIDUUM
            DFT_FIX_CORE_STATES = .true.
        endif

        ! Jesli jest spelniony warunek konieczny
        if(DFT_FIX_CORE_STATES) then
            no_core_states = 1
            do i = 1 , Widmo_NoStates
                if( Widmo_Evals(i) < Ef - 3*KbT ) then
                    no_core_states = i
                endif
            enddo
            ! szacowanie polozenia ostatniego stanu rdzeniowego
            ! ponizej tego stanu funkcje KS nie beda ponownie liczone
            no_core_states = max(no_core_states - 5,0)

            ! if there is no core states do nothing
            if(no_core_states == 0) then
                DFT_FIX_CORE_STATES = .false.
                print*," <warning> No core states equal zero!"
            else

                ! dokladne szukanie odpowiedniej energii i ostatniego
                ! stanu rdzeniowego. Stan "k" ten jest wybierany jako pierszy
                ! stan napotkany przy przeszukiwaniu od no_core_states do 1
                ! ktory jest oddalony od staniu  k+1 o energie znacznie wieksza
                ! niz srednia energia miedzy wszystkimi stanami w poprzedniej iteracji
                Core_Last_Enery = Widmo_Evals(no_core_states) + 0.1
                DFT_FIX_CORE_STATES = .false.
                do i = no_core_states , 2 , -1
                    if( abs(Widmo_Evals(i) - Widmo_Evals(i-1)) >  DFT_AVERAGE_DELTA_ENERGY ) then
                        Core_Last_Enery = (Widmo_Evals(i) + Widmo_Evals(i-1))/2
                        no_core_states  = i-1
                        DFT_FIX_CORE_STATES = .true.
                        exit
                    endif
                enddo

                print*,"Last core state:",no_core_states
                print*,"Last core energ:",Core_Last_Enery*Rd*1000
                print*,"Curr Fermi ener:",Ef*Rd*1000
                print*,"DeltaFermi ener:",abs(Ef-old_Ef)*Rd*1000
                print*,"Zmiana parametrow feast-a:"
                print*,"    Ef z     :",DFT_ATOMIC_EF," na:",max_Ef * Rd * 1000.0
                print*,"    L. stanow:",DFT_NO_STATES," na:",max_stanow  + 10
                DFT_ATOMIC_EF = max_Ef * Rd * 1000.0
                DFT_NO_STATES = max_stanow  + 10
            endif

        endif

        ! W razie powodzenia tworzona jest kopia stanow rdzeniowych
        if(DFT_FIX_CORE_STATES) then

            allocate(Core_Widmo_Evals(Widmo_NoStates))
            allocate(Core_Widmo_Vecs(size(Widmo_Vecs,1),Widmo_NoStates))

            Core_Widmo_Vecs  = Widmo_Vecs
            Core_Widmo_Evals = Widmo_Evals
        endif

        ! Liczba stanow do policzenenia minus wszystkie stany rdzeniowe
        propper_num_of_states = DFT_NO_STATES - no_core_states
        ! Jesli w poporzedniej iteracji mrozenie stanow bylo wlaczone to trzeba wyczyscic stany
        ! i zresetowac obliczenia
        if(CZY_OSTATNIO_ZAMROZONO_CORE_STATES == .true. .and. DFT_FIX_CORE_STATES==.false.) then
            print*," <c> Resetowanie problemu wlasnego"
            call spinsystem_widmo(Core_Last_Enery*Rd*1000,DFT_ATOMIC_EF,propper_num_of_states,2,8)
            call spinsystem_widmo(Core_Last_Enery*Rd*1000,DFT_ATOMIC_EF,propper_num_of_states,0,8)
        else
            print*," <c> Rozwiazywanie problemu wlasnego."
            call spinsystem_widmo(Core_Last_Enery*Rd*1000,DFT_ATOMIC_EF,propper_num_of_states,opt,8)
        endif
        ! Funkcja spinsystem_widmo w razie gdy podana przewidywana liczba stanow jest za mala
        ! moze zmienic wartosc propper_num_of_states, co mozna wykorzystac w nastepnej iteracji
        DFT_NO_STATES = propper_num_of_states + no_core_states


        ! Tu moze sie okazac ze po rozwiazaniu problemu wlasnego znaleziona liczba stanow
        ! moze byc mniejsza niz liczba donorow. Rachunki sa wtedy powtarzane dla nowych
        ! wartosci DFT_ATOMIC_EF i DFT_NO_STATES co powinno pozwolic na znalezienie odpowiedniej
        ! liczby stanow.
667     if(Widmo_NoStates + no_core_states < DFT_NO_DONORS) then

            DFT_ATOMIC_EF = DFT_ATOMIC_EF * (1 + 0.5)
            DFT_NO_STATES = DFT_NO_STATES * (1 + 0.5)
            print*,"Za mala energia. Nowa Ef:", DFT_ATOMIC_EF
            propper_num_of_states = DFT_NO_STATES-no_core_states
            call spinsystem_widmo(Core_Last_Enery*Rd*1000,DFT_ATOMIC_EF,propper_num_of_states,2,8)
            call spinsystem_widmo(Core_Last_Enery*Rd*1000,DFT_ATOMIC_EF,propper_num_of_states,0,8)

            goto 667
        endif

        ! Tutaj obliczone stany z okolicy energii Fermiego sa laczone
        ! ze stanami rdzeniowymi
        if(DFT_FIX_CORE_STATES) then

            print*," <c> Laczenie stanow. Obliczono: ",Widmo_NoStates," nowych stanow."

            ! alokowanie tablic pomocniczych
            allocate(Core_Widmo_Evals_Tmp(Widmo_NoStates + no_core_states))
            allocate(Core_Widmo_Vecs_Tmp(size(Widmo_Vecs,1),Widmo_NoStates + no_core_states))
            ! kopiowanie tylko stanow rdzeniowych
            do i = 1 , no_core_states
                Core_Widmo_Evals_Tmp(i) = Core_Widmo_Evals(i)
                Core_Widmo_Vecs_Tmp(:,i)= Core_Widmo_Vecs(:,i)
            enddo

            ! uzupenianie bazy o nowo obliczone stany - powyzej stanow rdzeniowych
            Core_Iter = 0
            do i = 1 , Widmo_NoStates
                Core_Iter = Core_Iter + 1
                Core_Widmo_Vecs_Tmp (:,no_core_states+Core_Iter) = Widmo_Vecs(:,i)
                Core_Widmo_Evals_Tmp(no_core_states+Core_Iter)   = Widmo_Evals(i)
            enddo

            Widmo_NoStates = no_core_states + Core_Iter ! Nowa liczba stanow
            ! Zwalnianie starych tablic i alokacja
            deallocate(Widmo_Evals)
            deallocate(Widmo_Vecs)
            allocate(Widmo_Evals(Widmo_NoStates))
            allocate(Widmo_Vecs(size(Core_Widmo_Vecs_Tmp,1),Widmo_NoStates))

            ! kopiowanie obliczonych i sklejonych stanow do wlasciwych tablic
            do i = 1 , Widmo_NoStates
                Widmo_Evals(i)  = Core_Widmo_Evals_Tmp(i)
                Widmo_Vecs(:,i) = Core_Widmo_Vecs_Tmp(:,i)
            enddo


            ! dealokacja
            deallocate(Core_Widmo_Evals)
            deallocate(Core_Widmo_Vecs)
            deallocate(Core_Widmo_Evals_Tmp)
            deallocate(Core_Widmo_Vecs_Tmp)
        endif ! end of if sklejanie funkcji rdzeniowych z obliczonymi


        else ! uzywaj stanow rdzeniowych
        print*," <c> Liczenie bez zamrazania stanow rdzeniowych!"
        print*," <c> Rozwiazywanie problemu wlasnego."
        call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,opt,8)

668     if(Widmo_NoStates < DFT_NO_DONORS) then
            print*,"Za mala energia."
            DFT_ATOMIC_EF = DFT_ATOMIC_EF * (1 + 0.5)
            DFT_NO_STATES = DFT_NO_STATES * (1 + 0.5)

            call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,2,8)
            call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,0,8)
            goto 668
        endif

        endif ! end of liczenie stanow rdzeniowych


    end subroutine spindft_ks_iteration


! ======================================================================
! OBLICZANIE WARTOSCI POTENCJALU ODDZIALYWANIA - HARTREE'EGO (KOHN-SHAMA)
! ======================================================================
    subroutine calcHartreePotential(bCalcHartreeDonor)
        logical,optional :: bCalcHartreeDonor


        integer :: m , k , ss , ib , jb , i ,j
        doubleprecision :: qhart
        print*," <c> Obliczanie potencjalu oddzialywania."

        m     = DFT_NO_INTEGRAL_REPEAT_X-1
        ss    = -m*(NX-1)


        ! oblicz potencjal od donorow
        if(present(bCalcHartreeDonor)) then
        u_donor_hartree = 0
        do ib  = 1 , NX
        do jb  = 1 , NY
            if(DFTINDEX(ib,jb) > 0) then
            qhart = 0.0
            do i = 1-m*NX , (m+1)*NX
            do j = 1 , NY
                 if(DFTINDEX(mod(abs(1-NX*m-i),NX)+1,j) > 0) then
                    qhart = qhart + dft_donor_dens*VSQRTL(abs( ib - i ), abs(j - jb), 1 )
                 endif
            enddo
            enddo
            u_donor_hartree(DFTINDEX(ib,jb)) = qhart*DX*DX
            endif
        enddo
        enddo
        ! po obliczeniu nie licz nic wiecej
        return
        endif


        ! Wykonanie dwuwymiarowej calki kulombowskiej
        do ib  = 1 , NX
        do jb  = 1 , NY
            if(DFTINDEX(ib,jb) > 0) then
            qhart = 0.0
            do i = 1-m*NX , (m+1)*NX
            do j = 1 , NY
                 if(DFTINDEX(mod(abs(1-NX*m-i),NX)+1,j) > 0) then
                    qhart = qhart - rho_tot(DFTINDEX(mod(abs(1-NX*m-i),NX)+1,j))*VSQRTL(abs(ib - i) ,abs(jb - j), 0 )
                 endif
            enddo
            enddo
            u_hartree(DFTINDEX(ib,jb)) = qhart*DX*DX + u_donor_hartree(DFTINDEX(ib,jb))
            endif
        enddo
        enddo


    end subroutine calcHartreePotential


! ======================================================================
! OBLICZANIE CALKOWITEJ GESTOSCI ELEKTRONOWEJ I GESTOSCI DLA POSZCZEGOLNYCH
!         SPINOW (UP/DOWN). Rozklad gestosci obliczany jest w oparciu o rozklad
!         fermiego diraca:
!
!                   Nele = SUM( 1/(exp(Ei-E_f/Kb*T)+1)*|Yi|^2 )
!
!         Energia fermiego i temperatura ustalane sa na poczatku programu.
! ======================================================================

    subroutine calcComplexEleDensity(reset)

    logical,intent(in),optional :: reset

    double precision :: fermi,Na,Efa,Efb,dens_ele,ndonor,dEF
    integer          :: val , s , i , j  , max_iter , ef_iter
    logical          :: fixTest
    double precision :: polaryzacje(-1:1)


    print*," <c> Obliczanie calkowitej gestosci elektronowej."

    ! Obliczanie E_f wyruwnojacej ladunek, normalizacja
    ndonor   = DFT_NO_DONORS
    old_Ef   = Ef
    max_iter = 100
    ef_iter  = 1
    Ef       = Widmo_Evals(ndonor)

    Efa      = Widmo_Evals(max(ndonor-10,1.0))
    Efb      = Widmo_Evals(min(ndonor+10.0,0.0+Widmo_NoStates))
    Na       = (calcfermiSum(Efa)-ndonor)

    ! Szukanie odpowiedniej Energi Fermiego metoda podzialu
!    do while( abs(calcfermiSum(Ef)-ndonor) > 1.0E-8 )
!
!      Ef = (Efa+Efb)/2
!
!      if ( Na*(calcfermiSum(Ef)-ndonor) < 0.0 ) then
!          Efb = Ef
!      else
!          Efa = Ef
!      endif
!      Na  = (calcfermiSum(Efa)-ndonor)
!      ef_iter = ef_iter + 1
!      if(ef_iter > max_iter) then
!        print*,"Error: nie udalo sie znalezc odpowiedniej energii."
!        Ef  = Widmo_Evals(ndonor)
!
!        exit
!      endif
!    enddo

    dEF = 0.0001/Rd/1000.0
    do while( abs(calcfermiSum(Ef)-ndonor) > 1.0E-10 )
      Efa = calcfermiSum(Ef+dEF)-ndonor
      Efb = calcfermiSum(Ef-dEF)-ndonor
      Ef  = Ef - (calcfermiSum(Ef)-ndonor) / ( ( Efa - Efb  ) / (2 * dEF) )
      ef_iter = ef_iter + 1
      print*,ef_iter,Ef*1000.0*Rd
      if(ef_iter > max_iter) then
        print*,"Error: nie udalo sie znalezc odpowiedniej energii."
        Ef  = Widmo_Evals(ndonor)
        exit
      endif
    enddo

    ! Ta funkcja nadpisuje tablice rho z nowa gestoscia otrzymana z rozkladu
    ! Fermiego-Diraca. Ta gestosc 'rho' stanowi output z czarnej skrzynki DFT
    ! Tymczasem 'rho_new' stanowilo input do tej czarnej skrzynki, w tej
    ! iteracji. Funkcja mix_densities dokona mieszania gestosci do nastepnej
    ! itercji zapisujac wynik w tablicy 'rho_new'
    dens_ele = calcfermiNele(Ef)


    ! SIMPLE MIXINIG
    if(present(reset)) then
       !new_rho = rho
       call mix_densities()
    else
       !new_rho = rho * DFT_W_PARAM + (1-DFT_W_PARAM) * new_rho
       call mix_densities()
    endif



    rho_tot  = new_rho(:,-1) + new_rho(:,+1)
    ! Sumowanie po gestosci elektronowej->liczba elektronow w pudle
    dens_ele = sum(rho_tot)*dx*dx


    if(TRANS_DEBUG == .true.) then
    ! Energie znajduja sie w pierwszej kolumnie eigvals
    print*," <i> Stany energetyczne i obsadzenia:"
    DO val = 1,Widmo_NoStates,Widmo_NoStates/15 + 1
        fermi = 1.0/(exp(((Widmo_Evals(val))-Ef)/KbT)+1.0)
        polaryzacje = 0
        do s = +1 , -1 , -2
        do i = 1 , nx
        do j = 1 , ny
            if( DFTINDEX(i,j) > 0 ) then
                 polaryzacje(s) = polaryzacje(s) + abs(Widmo_Vecs(GINDEX(i,j,s),val))**2
            endif
        enddo
        enddo
        enddo ! end of iii - numerowanie po spinach
        polaryzacje(1) = (polaryzacje(1)-polaryzacje(-1))/(sum(polaryzacje))
        print"(A,i4,A,e16.6,A,e16.6,A,e16.6)","     E[",val,"][meV]=", &
     &            1000*(DBLE(Widmo_Evals(val))  )*Rd,"  f[%]=",fermi," p[%]", &
     &           polaryzacje(1)
    end do ! end of do(wartosci wlasne)

    endif


    print*,"Obliczona  Energia Fermiego:",Ef*1000*Rd    ,"[meV]"
    print*,"Maksymalna Energia Fermiego:",max_Ef*1000*Rd,"[meV]"
    print*,"Maksymalna Liczba Stanow   :",max_stanow
    print*,"Srednia delta energi stanow:",DFT_AVERAGE_DELTA_ENERGY*1000*Rd  ,"[meV]"
    print*,"Obecna delta energii       :",DFT_CURR_DELTA_ENERGY*1000*Rd     ,"[meV]"

    if( abs(dens_ele - ndonor) > 0.01) then

        print*,"ERROR:Zle wyliczana warunek rownowagi ladunkowej"
        print*,"    ====================================="
        print("(A,f16.6)"),"     = L. don. w ukadzie:",ndonor
        print("(A,f16.6)"),"     = Roznica ele - don:",dens_ele - ndonor
        print("(A,f16.6)"),"     = L. rho up        :",sum(new_rho(:,+1))*DX*DX
        print("(A,f16.6)"),"     = L. rho down      :",sum(new_rho(:,-1))*DX*DX
        print*,"    ====================================="
    endif

    end subroutine calcComplexEleDensity



! ====================================================================== !
! =====================     FUNKCJA POMOCNICZA    ====================== !
!
! Funkcja obliczajaca liczbe eleketronow z uwzglednieniem rozkladu
!        Fermiego-Diraca, jako parametr wejsciowy podajemy energie
!        fermiego "Efermi", obliczamy calkowita gestosc dla danej
!        energii i zwaraca liczbe elektronow, z otrzymanej gestosci.
!
! ====================================================================== !
    double precision function calcFermiNele(Efermi) result(r_val)
        double precision,intent(in) :: Efermi
        integer          :: val,i,s,j
        double precision :: Ei,fermi
        rho  = 0      ! zerowanie gestosci
        max_stanow = 0
        max_Ef     = 0
        DFT_AVERAGE_DELTA_ENERGY = 0
        DO val = 1,Widmo_NoStates
            Ei = Widmo_Evals(val)   ! Energie znajduja sie w pierwszej kolumnie eigvals

            ! Obsiadanie stanow zgodnie ze statystyka Fermiego
            if( (Ei-Efermi)/KbT > 50 ) then
                fermi = 0.0
            else
                fermi = 1.0/(exp((Ei-Efermi)/KbT)+1.0)
            endif

            do s = +1 , -1 , -2
            do i = 1 , nx
            do j = 1 , ny
                if( DFTINDEX(i,j) > 0 ) then
                     rho(DFTINDEX(i,j),s) = rho(DFTINDEX(i,j),s) + &
                            fermi * abs(Widmo_Vecs(GINDEX(i,j,s),val))**2
                endif

            enddo
            enddo
            enddo ! end of iii - numerowanie po spinach

            ! Jesli obsadzenie jest male to mozna oszacowac maksymalna
            ! wartosc energii oraz liczbe stanow do FEASTA aby ten
            ! liczyl szybciej
            if(fermi > 0.000001) then
                max_stanow = val
                max_Ef     = Widmo_Evals(max_stanow)*1.0
            endif
            ! zliczanie sredniej odleglosci miedzy stanami
            if(val < Widmo_NoStates) then
            DFT_AVERAGE_DELTA_ENERGY = DFT_AVERAGE_DELTA_ENERGY + (Widmo_Evals(val+1) - Widmo_Evals(val))
            endif
        end do ! end of do(wartosci wlasne)

        ! sredniowanie
        DFT_AVERAGE_DELTA_ENERGY = abs(DFT_AVERAGE_DELTA_ENERGY) / (Widmo_NoStates-1)
        DFT_CURR_DELTA_ENERGY    = abs(Ef - old_Ef)
        r_val = sum(rho)*dx*dx
    end function calcFermiNele



! ====================================================================== !
! =====================     FUNKCJA POMOCNICZA 2  ====================== !
!
! Funkcja obliczajaca liczbe stanow obsadzonych z rozkladu F-D dla danej
!        Energii fermiego Efermi
!
! ====================================================================== !
    double precision function calcFermiSum(Efermi) result(rval)
        double precision,intent(in) :: Efermi

        integer          :: val
        double precision :: Nstates
        double precision :: Ei,fermi

        Nstates    = 0
        DO val = 1, Widmo_NoStates
          Ei = Widmo_Evals(val)
          ! Obsadzanie stanow zgodnie ze statystyka Fermiego
          if( (Ei-Efermi)/KbT > 50 ) then
              fermi = 0.0
          else
              fermi   = 1.0/(exp((Ei-Efermi)/KbT)+1.0)
          endif
          Nstates   = fermi  + Nstates
        end do ! end of do(wartosci wlasne)

        rval = Nstates

    end function calcFermiSum



! ======================================================================
!
!           Mieszanie gestosci wedlug zastosowanego schematu
!
! ======================================================================
    subroutine mix_densities()
        doubleprecision,allocatable :: tmpArray(:)

        allocate(tmpArray(2*DFT_TRANSMAX))


        tmpArray(1:DFT_TRANSMAX)                = new_rho(:,+1)
        tmpArray(DFT_TRANSMAX+1:2*DFT_TRANSMAX) = new_rho(:,-1)

        call mixer%set_input_vec (tmpArray)


        tmpArray(1:DFT_TRANSMAX)                = rho(:,+1)
        tmpArray(DFT_TRANSMAX+1:2*DFT_TRANSMAX) = rho(:,-1)

        call mixer%set_output_vec(tmpArray)
        deallocate(tmpArray)
        call mixer%mix()

        new_rho(:,+1) = mixer%mixedVec(1:DFT_TRANSMAX)
        new_rho(:,-1) = mixer%mixedVec(DFT_TRANSMAX+1:2*DFT_TRANSMAX)


        DFT_CURR_RESIDUUM =  abs((old_Ef - Ef)/(old_Ef + Ef))


    end subroutine mix_densities




    subroutine spindft_zapisz_rho()
        integer :: i
        open(unit = 432, file= "r.dat" )
        do i = 1 , DFT_TRANSMAX
            write(432,"(2f20.6)"),new_rho(i,+1),new_rho(i,-1)
        enddo
        close(432)
    end subroutine spindft_zapisz_rho

    subroutine spindft_wczytaj_rho_poczatkowe()
        integer :: i
        doubleprecision :: rvals(2)
        open(unit = 432, file= "r.dat" )
        do i = 1 , DFT_TRANSMAX
            read(432,"(2f20.6)"),rvals
            new_rho(i,+1) = rvals(1)
            new_rho(i,-1) = rvals(2)
        enddo
        close(432)
        rho_tot =   new_rho(i,+1) + new_rho(i,-1)
        call spindft_ks_iteration(0)


    end subroutine spindft_wczytaj_rho_poczatkowe



!    subroutine spindft_ks_from_basis_iteration()
!
!        integer :: i,j,s,m,n
!        doubleprecision,allocatable :: ks_potential(:)
!
!        print*," <c> Rozwiazywanie problemu wlasnego dla rozwiniecia w bazie"
!
!        allocate(ks_potential(size(Basis_EVecs,1)))
!
!        call XAttaccalite_1(new_rho(:,1),new_rho(:,-1),xc_res)
!        u_exchange(:,+1) = xc_res(:,1)
!        u_exchange(:,-1) = xc_res(:,2)
!
!        call CAttaccalite_1(new_rho(:,1),new_rho(:,-1),xc_res)
!        u_correlation(:,+1) = xc_res(:,1)
!        u_correlation(:,-1) = xc_res(:,2)
!
!        call calcHartreePotential()
!
!        ! sumowanie wszystkich przyczynkow
!        do s = 1 , -1  , -2
!        do i = 1 , nx
!        do j = 1 , ny
!            if(DFTINDEX(i,j) > 0 ) then
!
!                SUTOTAL(i,j,s) =   u_exchange(DFTINDEX(i,j),s) &
!                                 + u_correlation(DFTINDEX(i,j),s) &
!                                 - u_hartree(DFTINDEX(i,j))
!
!                ks_potential(GINDEX(i,j,s)) = SUTOTAL(i,j,s) - Basis_InitialPotential(GINDEX(i,j,s))
!
!            endif
!        enddo
!        enddo
!        enddo
!
!
!        do m = 1 , Basis_No_States
!        do n = 1 , Basis_No_States
!            Basis_KSHamiltonian(m,n) = sum(ks_potential*conjg(Basis_EVecs(:,m))*Basis_EVecs(:,n))*dx*dx
!            if( m == n ) Basis_KSHamiltonian(m,n) = Basis_KSHamiltonian(m,n) + Basis_EVals(m)
!        enddo
!        enddo
!
!
!        call spindft_solve_dense_evproblem(Basis_KSHamiltonian,Basis_KSEVectors,Basis_KSEnergies)
!
!
!        Widmo_Evals = Basis_KSEnergies
!
!        Widmo_Vecs  = 0
!        do n = 1 , Basis_No_States
!        do m = 1 , Basis_No_States
!            Widmo_Vecs(:,n) = Widmo_Vecs(:,n) + Basis_KSEVectors(n,m)*Basis_EVecs(:,m)
!        enddo
!        enddo
!
!        deallocate(ks_potential)
!    end subroutine spindft_ks_from_basis_iteration
!
!! =============================================================================
!!   Obliczanie problemu wlasnego dla zagadnienia z gestymi macierzami
!! =============================================================================
!
!    subroutine spindft_solve_dense_evproblem(Hamiltonian,Evecs,Evals)
!        complex*16,dimension(:,:)    :: Hamiltonian,Evecs
!        doubleprecision,dimension(:) :: Evals
!
!        INTEGER          N
!        INTEGER          LDA
!        INTEGER          LWMAX
!
!        !     .. Local Scalars ..
!        INTEGER          INFO, LWORK
!        !     .. Local Arrays ..
!        !*     RWORK dimension should be at least MAX(1,3*N-2)
!        DOUBLE PRECISION,allocatable,dimension(:) ::  W, RWORK
!        COMPLEX*16,allocatable ::      A( :, : ), WORK( : )
!
!
!        print*," <c> Rozwiazywanie problemu wlasnego dla funkcji bazowych."
!
!        N     = size(Hamiltonian,1)
!        LDA   = N
!        LWORK = 3*N - 2
!
!
!        allocate(W(N))
!        allocate(A(LDA,N))
!
!        A = Hamiltonian
!
!        allocate(WORK(LWORK))
!        allocate(RWORK(LWORK))
!!*
!!*     Query the optimal workspace.
!!*
!        LWORK = -1
!        CALL ZHEEV( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK, INFO )
!
!        LWORK = WORK( 1 )
!
!        deallocate(WORK)
!        allocate(WORK(LWORK))
!
!
!
!        !print*,"optimal work place:",LWORK,N
!
!!*
!!*     Solve eigenproblem.
!!*
!        CALL ZHEEV( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,INFO )
!!*
!!*     Check for convergence.
!!*
!        IF( INFO.GT.0 ) THEN
!         WRITE(*,*)'Error: The algorithm failed to compute eigenvalues. Info:',INFO
!         STOP
!        else
!         print*," <c> Rozwiazanie OK..."
!        END IF
!
!        deallocate(W)
!        deallocate(A)
!        deallocate(WORK)
!        deallocate(RWORK)
!
!    end subroutine spindft_solve_dense_evproblem
!!*  =============================================================================
!!*

endmodule
