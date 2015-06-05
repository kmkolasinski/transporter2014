module modspindft
    use modspinsystem
    use modjed
    use modinip
    use modutils
    use xc
    use ifport
    implicit none
    private
    doubleprecision, allocatable :: u_exchange(:,:)
    doubleprecision, allocatable :: u_correlation(:,:)
    doubleprecision, allocatable :: u_hartree(:) , u_donor_hartree(:)
    doubleprecision, allocatable :: rho(:,:),new_rho(:,:),rho_tot(:) , xc_res(:,:)
    doubleprecision, allocatable :: rho_mem(:,:,:) , rho_mem_in(:,:,:)  , broyden_R(:,:,:) , diagonal_Jm(:,:) , broyden_U(:,:,:) , broyden_V(:,:,:)

    doubleprecision, allocatable :: VSQRTL(:,:,:)


    ! Broyden Eli
    doubleprecision, allocatable :: rho_pp(:,:) , rho_p(:,:) , rho_or_p(:,:) , rhoU_or(:,:)
    doubleprecision, allocatable :: dR_u_p(:,:) , u_p(:,:,:) , v_p(:,:,:)  ,u_n(:,:) , v_n(:,:)


    integer,dimension(:,:),allocatable   :: DFTINDEX
    integer          :: DFT_NO_INTEGRAL_REPEAT_X ! w przypadku periodycznych warunkow brzegowych
    integer          :: DFT_TRANSMAX , NX , NY , DFT_NO_STATES , GLOBAL_ITER , GLOBAL_TEMP_ITER
    integer          :: DFT_DENS_NO_MEM , DFT_MAX_ITER
    integer          :: DFT_IMPROVE_FEAST_STEPS  , DFT_TEMP_NO_STEPS
    logical          :: DFT_FIX_EF ! jesli true to wtedy gestosc nie jest wyliczana na podstawie bilansu ladonkow
    double precision :: DFT_W_PARAM , DFT_MAX_W_PARAM
    double precision :: DFT_NO_DONORS , DFT_Z_SPACER
    double precision :: DFT_TEMP ,  DFT_TEMP_MIN  , DFT_CURR_TEMP  , DFT_RESIDUUM , DFT_CURR_RESIDUUM
    doubleprecision  :: DFT_FINDED_EF , DFT_FIXED_ENERGY
    doubleprecision  :: Ef , Bz , Bx , By , DX , DFT_ATOMIC_EF
    doubleprecision  :: max_stanow , max_Ef , dft_donor_dens

    public :: spindft_initialize , spindft_free , spindft_solve , spindft_fix_ef
    public :: spindft_solve_temp_annealing
    public :: spindft_zapisz_rho , spindft_wczytaj_rho_poczatkowe
    public :: DFT_TEMP ,  DFT_TEMP_MIN  , DFT_MAX_ITER , DFT_CURR_RESIDUUM , DFT_FINDED_EF ,DFT_FIX_EF
    contains



    subroutine spindft_initialize()
        integer ::  i , j , iter , mx

        DFT_FIX_EF = .false.

        BZ  = BtoDonorB(atomic_Bz)
        Bx  = BtoDonorB(atomic_Bx)
        By  = BtoDonorB(atomic_By)
        DX  = atomic_DX*L2LR
        nx = size(GINDEX,1)
        ny = size(GINDEX,2)


        print*,"Inicjalizacja problemu SPIN-DFT"

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
        call getIntValue   ("DFT","init_states",DFT_NO_STATES)
        call getIntValue   ("DFT","improve_feast",DFT_IMPROVE_FEAST_STEPS)
        print*,"Problem uzbiezniany na podstawie liczby donorow: DFT_FIX_EF = false"

        DFT_DENS_NO_MEM = 2!DFT_MAX_ITER

        Ef = DFT_ATOMIC_EF/1000.0/Rd

        call spindft_free()

        allocate(DFTINDEX(nx,ny))

        DFTINDEX = 0
        iter   = 1
        do i = 1 , nx
        do j = 1 , ny
            if( GFLAGS(i,j) == B_NORMAL  ) then
                DFTINDEX(i,j) = iter
                iter = iter + 1
           endif
        enddo
        enddo

        if(TRANS_EIGPROBLEM_PERIODIC_X) then
            do j = 1 , ny
                DFTINDEX(nx,j) = DFTINDEX(2   ,j)
                DFTINDEX(1 ,j) = DFTINDEX(nx-1,j)
            enddo
        endif

        DFT_TRANSMAX = iter - 1 ! ilosc oczek na pojedynczy spin

        dft_donor_dens = DFT_NO_DONORS / DFT_TRANSMAX / DX / DX

        allocate(u_hartree(DFT_TRANSMAX))
        allocate(u_donor_hartree(DFT_TRANSMAX))
        allocate(u_exchange(DFT_TRANSMAX,-1:1))
        allocate(u_correlation(DFT_TRANSMAX,-1:1))
        allocate(rho(DFT_TRANSMAX,-1:1))

        allocate(new_rho(DFT_TRANSMAX,-1:1))
        allocate(rho_tot(DFT_TRANSMAX))
        allocate(xc_res(DFT_TRANSMAX,3))

        allocate(rho_mem(DFT_TRANSMAX,-1:1,DFT_DENS_NO_MEM))
        allocate(rho_mem_in(DFT_TRANSMAX,-1:1,DFT_DENS_NO_MEM))
        allocate(broyden_R(DFT_TRANSMAX,-1:1,DFT_DENS_NO_MEM))
        allocate(broyden_U(DFT_TRANSMAX,-1:1,DFT_DENS_NO_MEM))
        allocate(broyden_V(DFT_TRANSMAX,-1:1,DFT_DENS_NO_MEM))
        allocate(diagonal_Jm(DFT_TRANSMAX,-1:1))

        allocate(rho_pp(DFT_TRANSMAX,-1:1))
        allocate(rho_p(DFT_TRANSMAX,-1:1))
        allocate(rho_or_p(DFT_TRANSMAX,-1:1))
        allocate(rhoU_or(DFT_TRANSMAX,-1:1))
        allocate(dR_u_p(DFT_TRANSMAX,-1:1))
        allocate(v_n(DFT_TRANSMAX,-1:1))
        allocate(u_n(DFT_TRANSMAX,-1:1))
        allocate(u_p(DFT_TRANSMAX,-1:1,DFT_DENS_NO_MEM))
        allocate(v_p(DFT_TRANSMAX,-1:1,DFT_DENS_NO_MEM))
        v_n = 0
        u_n = 0
        v_p = 0
        u_p = 0
        dR_u_p = 0


        if(TRANS_EIGPROBLEM_PERIODIC_X) then
            DFT_NO_INTEGRAL_REPEAT_X = 8
            allocate(VSQRTL(0:nx*DFT_NO_INTEGRAL_REPEAT_X,0:ny,0:1)) ! dwie warstwy 0 (na gaz) i 1 (na donory)
        else
            allocate(VSQRTL(0:nx,0:ny,0:1)) ! dwie warstwy 0 (na gaz) i 1 (na donory)
            DFT_NO_INTEGRAL_REPEAT_X = 1
        endif


        ! obliczanie pierwiastkow
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
        diagonal_Jm = 0
        diagonal_Jm(:,-1) = DFT_W_PARAM
        diagonal_Jm(:,+1) = DFT_W_PARAM


        rho_tot   = 0 !      dft_donor_dens
        rho       = 0 ! 0.5*dft_donor_dens
        new_rho   = 0 !      dft_donor_dens
!        do i = 1 , DFT_TRANSMAX
!            new_rho(i,-1) = rand()
!            new_rho(i,+1) = rand()
!        enddo
!        new_rho =new_rho / sum(new_rho)*dx*dx * DFT_NO_DONORS
!        rho_tot  = new_rho(:,-1) + new_rho(:,+1)

        rho_mem   = 0
        rho_mem_in= 0
        broyden_R = 0
        broyden_U = 0
        broyden_V = 0


        DFT_FINDED_EF = 0


    end subroutine spindft_initialize

    ! podajemy energie w jednostach [meV]
    subroutine spindft_fix_ef(fixed_ef)
        doubleprecision ::fixed_ef
        DFT_FIX_EF = .true.
        Ef = fixed_ef / 1000.0 / Rd
        DFT_FIXED_ENERGY = Ef
    end subroutine


    subroutine spindft_solve_temp_annealing()
        integer :: itemp
        doubleprecision :: tomega
        GLOBAL_TEMP_ITER = 0
        open(unit = 8998, file= "dft_iter.txt" )
        do itemp = 1 , DFT_TEMP_NO_STEPS
            tomega = (itemp-1.0)/(DFT_TEMP_NO_STEPS-1.0)

            if(DFT_TEMP_NO_STEPS == 1) tomega = 0.0

            tomega = 1-(1-tomega)**2


            DFT_CURR_TEMP = DFT_TEMP*(1 - tomega) + tomega*DFT_TEMP_MIN
            call spindft_solve()


            print*,"Symulacja uzbiezniona, po krokach:",GLOBAL_ITER
            print*,"Zakonczony krok czasowy          :",DFT_TEMP_NO_STEPS
            print*,"Z residuum rownym                :",DFT_CURR_RESIDUUM
            print*,"Z temperatura                    :",DFT_CURR_TEMP
            DFT_FINDED_EF = Ef*1000*Rd
            print*,"Znaleziona energia Fermiego      :",DFT_FINDED_EF

        enddo
        print*,"-----------------------------------------------------"
        print*,"Symulacja uzbiezniona w sumie po krokrach:",GLOBAL_TEMP_ITER
        print*,"Z residuum rownym:",DFT_CURR_RESIDUUM
        DFT_FINDED_EF = Ef*1000*Rd
        print*,"Znaleziona energia Fermiego:",DFT_FINDED_EF

    end subroutine spindft_solve_temp_annealing



    subroutine spindft_solve()
        integer :: i,j,m
        ! obliczenia startowe dla gestosci elektronowej takiej samej
        ! jak gestosc donorow
        GLOBAL_ITER = 1
        call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,0,8)
668     if(Widmo_NoStates < DFT_NO_DONORS) then
            print*,"Za mala energia."
            DFT_ATOMIC_EF = DFT_ATOMIC_EF * (1 + 0.5)
            DFT_NO_STATES = DFT_NO_STATES * (1 + 0.5)

            call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,2,8)
            call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,0,8)
            goto 668
        endif

        call calcComplexEleDensity(reset =.true.)

        DFT_FINDED_EF = 0
        do GLOBAL_ITER = 2 , DFT_MAX_ITER
            GLOBAL_TEMP_ITER  = GLOBAL_TEMP_ITER + 1

            print*,"Iteracja DFT:",GLOBAL_ITER
            kbT       = 8.617D-5/Rd*DFT_CURR_TEMP
            print*,"Temperatura", kbT * Rd / 8.617D-5


            if(MOD(GLOBAL_ITER,DFT_IMPROVE_FEAST_STEPS) == 0 .or. GLOBAL_ITER == 5) then
                print*,"Zmiana parametrow feast-a:"

                DFT_ATOMIC_EF = max_Ef * Rd * 1000.0
                DFT_NO_STATES = max_stanow*1.5
                call spindft_ks_iteration(0)
            else
                call spindft_ks_iteration(1)
            endif


            call calcComplexEleDensity()
            rho = new_rho
            do i = 1 , nx
            do j = 1 , ny
                if(DFTINDEX(i,j) > 0 .and. DFTINDEX(i,ny+1-j) > 0 )  then
                    new_rho(DFTINDEX(i,j),:)  = (rho(DFTINDEX(i,j),:) + rho(DFTINDEX(i,ny+1-j),:))/2
                endif
            enddo
            enddo

            write(8998,"(20f20.6)"),GLOBAL_TEMP_ITER+0.0D0,DFT_CURR_TEMP,GLOBAL_ITER+0.0,Ef*1000*Rd,DFT_CURR_RESIDUUM


            open(unit=222,file="dft_rho.txt")
            do i = 1 , nx
            do j = 1 , ny
                if(DFTINDEX(i,j) > 0 ) &
                write(222,"(20f20.6)"),i*atomic_DX,j*atomic_DX,new_rho(DFTINDEX(i,j),1),new_rho(DFTINDEX(i,j),-1)
            enddo
                write(222,*),""
            enddo
            close(222)

            open(unit=555,file="dft_uint.txt")
            do i = 1 , nx
            do j = 1 , ny
                if(DFTINDEX(i,j) > 0 ) then
                    write(555,"(5f20.6)"),i*atomic_DX,j*atomic_DX,SUTOTAL(i,j,1)*Rd*1000.0,SUTOTAL(i,j,-1)*Rd*1000.0
                else
                    write(555,"(5f20.6)"),i*atomic_DX,j*atomic_DX,0.0D0,0.0D0
                endif
            enddo
                write(555,*),""
            enddo
            close(555)

            if(  DFT_CURR_RESIDUUM < DFT_RESIDUUM  ) then
                print*,"Finished"
                exit
            endif

        enddo


        call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,2,8)

    end subroutine spindft_solve



    subroutine spindft_free()
        if(allocated(DFTINDEX))   deallocate(DFTINDEX)
        if(allocated(u_hartree))  deallocate(u_hartree)
        if(allocated(u_donor_hartree))  deallocate(u_donor_hartree)
        if(allocated(u_exchange))  deallocate(u_exchange)
        if(allocated(u_correlation)) deallocate(u_correlation)
        if(allocated(rho))         deallocate(rho)
        if(allocated(rho_mem))     deallocate(rho_mem)
        if(allocated(rho_mem_in))  deallocate(rho_mem_in)
        if(allocated(broyden_R))   deallocate(broyden_R)
        if(allocated(broyden_U))   deallocate(broyden_U)
        if(allocated(broyden_V))   deallocate(broyden_V)
        if(allocated(diagonal_Jm)) deallocate(diagonal_Jm)
        if(allocated(new_rho))     deallocate(new_rho)
        if(allocated(rho_tot))     deallocate(rho_tot)
        if(allocated(xc_res))      deallocate(xc_res)
        if(allocated(VSQRTL))      deallocate(VSQRTL)


        if(allocated(rho_pp)) deallocate(rho_pp)
        if(allocated(rho_p)) deallocate(rho_p)
        if(allocated(rho_or_p)) deallocate(rho_or_p)
        if(allocated(rhoU_or)) deallocate(rhoU_or)
        if(allocated(dR_u_p)) deallocate(dR_u_p)
        if(allocated(v_n)) deallocate(v_n)
        if(allocated(u_n)) deallocate(u_n)
        if(allocated(u_p)) deallocate(u_p)
        if(allocated(v_p)) deallocate(v_p)


        ! Zwalnianie pamieci
        call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,2,8)
    end subroutine spindft_free


    subroutine spindft_ks_iteration(opt)
        integer :: opt
        integer :: i,j,s

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
                SUTOTAL(i,j,s) = + u_exchange(DFTINDEX(i,j),s) &
                                 + u_correlation(DFTINDEX(i,j),s) &
                                 - u_hartree(DFTINDEX(i,j))
            endif
        enddo
        enddo
        enddo

        print*," <c> Rozwiazywanie problemu wlasnego."
        call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,opt,8)

667     if(Widmo_NoStates < DFT_NO_DONORS) then
            print*,"Za mala energia."
            DFT_ATOMIC_EF = DFT_ATOMIC_EF * (1 + 0.5)
            DFT_NO_STATES = DFT_NO_STATES * (1 + 0.5)

            call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,2,8)
            call spinsystem_widmo(0.0D0,DFT_ATOMIC_EF,DFT_NO_STATES,0,8)
            goto 667
        endif

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



        do ib  = 1 , NX
        do jb  = 1 , NY
            if(DFTINDEX(ib,jb) > 0) then
            qhart = 0.0
            do i = 1-m*NX , (m+1)*NX
            do j = 1 , NY
                 if(DFTINDEX(mod(abs(1-NX*m-i),NX)+1,j) > 0) then
                    qhart = qhart - rho_tot(DFTINDEX(mod(abs(1-NX*m-i),NX)+1,j))*VSQRTL(abs(ib - i) ,abs(jb - j), 0 )
                    !qhart = qhart + dft_donor_dens*VSQRTL(abs( ib - i ), abs(j - jb), 1 )
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

    double precision :: fermi,Na,Efa,Efb,dens_ele,ndonor
    integer          :: val , s , i , j
    logical          :: fixTest
    double precision :: polaryzacje(-1:1)


    print*," <c> Obliczanie calkowitej gestosci elektronowej."


    ! Obliczanie E_f wyruwnojacej ladunek, normalizacja

    ndonor = DFT_NO_DONORS

    ! jesli energia fermiego jest ustalona to wtedy nie liczymy gestosci
    ! z zachowania ladunku
    fixTest = GLOBAL_TEMP_ITER < 10

    if(DFT_FIX_EF == .false. .or. fixTest) then
        Efa = -100
        Efb =  100.0
        Na  = (calcfermiSum(Efa)-ndonor)
        do while( abs(calcfermiSum(Ef)-ndonor) > 1.0E-8 )
          Ef = (Efa+Efb)/2
          if ( Na*(calcfermiSum(Ef)-ndonor) < 0.0 ) then
              Efb = Ef
          else
              Efa = Ef
          endif
          Na  = (calcfermiSum(Efa)-ndonor)
        enddo
    else
        Ef = DFT_FIXED_ENERGY
    endif



    dens_ele = calcfermiNele(Ef)


    ! SIMPLE MIXINIG
    if(present(reset)) then
       call mix_densities(rho,new_rho)
    else
       !new_rho = rho*(DFT_W_PARAM ) + new_rho*(1.0 - DFT_W_PARAM )
       call mix_densities(rho,new_rho)
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

    print*,"Obliczona  Energia Fermiego:",Ef*1000*Rd,"[meV]"
    print*,"Maksymalna Energia Fermiego:",max_Ef*1000*Rd,"[meV]"
    print*,"Maksymalna Liczba Stanow   :",max_stanow

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
        DO val = 1,Widmo_NoStates
            Ei = Widmo_Evals(val)   ! Energie znajduja sie w pierwszej kolumnie eigvals

            ! Obsiadanie stanow zgodnie ze statystyka Fermiego
            if( (Ei-Efermi)/KbT > 500 ) then
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

            if(fermi > 0.00001) then
                max_stanow = val
                max_Ef     = Widmo_Evals(max_stanow)*1.1
            endif
        end do ! end of do(wartosci wlasne)

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
          if( (Ei-Efermi)/KbT > 500 ) then
              fermi = 0.0
          else
              fermi   = 1.0/(exp((Ei-Efermi)/KbT)+1.0)
          endif
          Nstates   = fermi  + Nstates
        end do ! end of do(wartosci wlasne)

        rval = Nstates

    end function calcFermiSum


!    subroutine mix_densities(rho_or,rho)
!        doubleprecision, allocatable :: rho_or(:,:),rho(:,:)
!        integer :: mem , s , i , it_u
!        doubleprecision ::  calka_dRU , calka_dRD
!
!! ! 	it=0 - w. poczatkowy, nie ma mieszania
!! ! 	it=1 - rho_1 = rho_0 normalnie zmieszane z rho_or_1
!! ! 	it=2 - rho_2 = rho_1 normalnie zmieszane z rho_or_2
!! ! ! 	  mieszanie:
!	    if(GLOBAL_ITER .ne. 1 ) then
!        if(GLOBAL_ITER  <= 3  ) then
!            write(*,*) GLOBAL_ITER," jeszcze zwykle mieszanie"
!            rho  =rho_or*DFT_W_PARAM + rho_p*(1.D0-DFT_W_PARAM)
!        else
!            write(*,*) GLOBAL_ITER," mieszanie Broydena"
!            it_u = GLOBAL_ITER - 2
!
!            calka_dRU=0.D0; calka_dRD=0.D0
!            dR_u_p = rho_or - rho_p - rho_or_p + rho_pp
!
!            do i=1,DFT_TRANSMAX
!              calka_dRU=calka_dRU+abs(dR_u_p(i,+1))**2
!              calka_dRD=calka_dRD+abs(dR_u_p(i,-1))**2
!            enddo
!
!            call multiply_Jac(it_u-1,u_p,v_p,dR_u_p,u_n)
!            u_p(:,:,it_u)=  rho_p - rho_pp - u_n
!
!            v_p(:,+1,it_u)=dR_u_p(:,+1)/calka_dRU
!            v_p(:,-1,it_u)=dR_u_p(:,-1)/calka_dRD
!
!            u_n = rho_or - rho_p
!            call multiply_Jac(it_u,u_p,v_p,u_n,v_n)
!            rho= -v_n + rho_p
!
!!! ! ! 		zaalokowac tablice do it_u+1
!!		call expand_uv(it_u,uU_p,uD_p,vU_p,vD_p)
!	      endif
!	    else
!          write(*,*) GLOBAL_ITER," przepisywanie gestosci"
!          rho = rho_or
!	    endif
!! ! 	  przepisanie tablic:
!	    rho_pp     =rho_p
!	    rho_p      =rho
!	    rho_or_p   =rho_or
!
!        DFT_CURR_RESIDUUM = sum(abs(rho_or-rho))*dx*dx


    subroutine mix_densities(in_rho,out_rho)
        doubleprecision, allocatable :: in_rho(:,:),out_rho(:,:)
        integer :: mem , s , i , it_u

        print*,"Residuum:"
        ! przesuwanie tablic
        do mem = DFT_DENS_NO_MEM , 2 , -1
            rho_mem_in(:,:,mem) = rho_mem_in(:,:,mem-1)
            rho_mem  (:,:,mem)  = rho_mem  (:,:,mem-1)
            broyden_R(:,:,mem)  = broyden_R(:,:,mem-1)
            broyden_U(:,:,mem)  = broyden_U(:,:,mem-1)
            broyden_V(:,:,mem)  = broyden_V(:,:,mem-1)
            print*,"Past:",mem," residuum:",sum(abs(broyden_V(:,1,mem)))*dx*dx
        enddo

        ! w tablicy
        broyden_R(:,:,1)  = in_rho - rho_mem(:,:,1) ! bierzemy roznice Fx - x




        DFT_CURR_RESIDUUM = sum(abs(broyden_R(:,1,1)))*dx*dx
        print*,"Past:",1," residuum:",DFT_CURR_RESIDUUM
        print*,""

        if(GLOBAL_ITER == 1) then
            out_rho = in_rho
        else  if(GLOBAL_ITER == 2) then! proste mieszanie
            out_rho = rho_mem(:,:,1) + DFT_W_PARAM * broyden_R(:,:,1)
        else
            ! mieszanie na podstawie: http://www.ondrejcertik.com/media/cookbook/master.pdf
            ! wzor (3.28)
            do i = 1 , DFT_TRANSMAX
                do s = 1 , -1 , -2
                if( broyden_R(i,s,1)*broyden_R(i,s,2) > 0 ) then
                    diagonal_Jm(i,s) = diagonal_Jm(i,s) + DFT_W_PARAM
                    if(diagonal_Jm(i,s) > DFT_MAX_W_PARAM) diagonal_Jm(i,s) = DFT_MAX_W_PARAM
                else
                    diagonal_Jm(i,s) =  DFT_W_PARAM
                endif
                enddo ! of s
            enddo
            out_rho = rho_mem(:,:,1) + diagonal_Jm * broyden_R(:,:,1)

!            do s = 1 , -1 , -2
!                broyden_U(:,s,1) = (rho_mem(:,s,2) - rho_mem(:,s,3)) + DFT_W_PARAM * ( broyden_R(:,s,1) - broyden_R(:,s,2) )
!                do mem = 2 , DFT_DENS_NO_MEM
!                    broyden_U(:,s,1) = broyden_U(:,s,1) - broyden_U(:,s,mem) * sum(broyden_V(:,s,mem)* ( broyden_R(:,s,1) - broyden_R(:,s,2) ) )
!                enddo
!                broyden_V(:,s,1) = ( broyden_R(:,s,1) - broyden_R(:,s,2) )/sum(abs(broyden_R(:,s,1) - broyden_R(:,s,2))**2)
!            enddo
!
!            do s = 1 , -1 , -2
!                out_rho(:,s) =  rho_mem(:,s,1) + DFT_W_PARAM * ( broyden_R(:,s,1) )
!                do mem = 1 , DFT_DENS_NO_MEM
!                    out_rho(:,s) = out_rho(:,s) - broyden_U(:,s,mem) * sum(broyden_V(:,s,mem)* ( broyden_R(:,s,1) ) )
!                enddo
!            enddo


        endif

        rho_mem(:,:,1)   = out_rho(:,:) ! aktualizujemy gestosc na pierwszym miejscu

    end subroutine mix_densities


	subroutine multiply_Jac(size_u,tab_u,tab_v,tab_i,tab_o)

	  integer, intent(in)::size_u
	  double precision,dimension(:,:),allocatable,intent(in)::tab_i
	  double precision,dimension(:,:,:),allocatable,intent(in)::tab_u,tab_v
	  double precision,dimension(:,:),allocatable::tab_o
	  integer ::  is,ij1,ij2

	  tab_o= -DFT_W_PARAM * tab_i

	  do is=1,size_u
	    do ij1=1,DFT_TRANSMAX
	      do ij2=1,DFT_TRANSMAX
            tab_o(ij1,+1)=tab_o(ij1,+1)+tab_u(ij1,+1,is)*tab_v(ij2,+1,is)*tab_i(ij2,+1)
            tab_o(ij1,-1)=tab_o(ij1,-1)+tab_u(ij1,-1,is)*tab_v(ij2,-1,is)*tab_i(ij2,-1)
	      enddo
	    enddo
	  enddo

	end subroutine multiply_Jac


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


endmodule
