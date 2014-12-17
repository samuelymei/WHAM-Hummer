!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for the calculation of the free energy from biased simulation like umbrella sampling !
! and temperature replica exchange molecular dynamics simualtion                              !
! Written by                                                                                  !
!                                        Ye Mei                                               !
!                                      11/17/2014                                             !
!                           East China Normal University                                      !
! Reference:                                                                                  !
!      J. Comput. Chem. 33, 453-465 (2012)                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
module WHAM
  use precision_m
  use control
  use simulation
  use react_coord_bin
  implicit none
  
  private
  integer(kind=4), public :: NumW  ! number of simulations
  integer(kind=4), public :: NumJ  ! number of reaction coordinate bins
  integer(kind=4), public :: NumK  ! number of energy bins
  integer(kind=4), public :: NumB  ! number of bins
  real(kind=fp_kind) :: tolerance = 1.0D-4  ! convergence criterion for free energy 
  integer(kind=4), parameter :: MaxITS = 1000          ! max number of iterations
  real(kind=fp_kind), public :: T_target = 300.0      ! target temperature
  integer(kind=4), public :: idOutputFile = 11
  integer(kind=4) :: NumBootstrap = 100

  real(kind=fp_kind) :: beta ! inverse temperature

  integer(kind=4) :: indexW ! index of windows/simulations
  integer(kind=4) :: indexS ! index of snapshot
  integer(kind=4) :: indexK ! index of energy bin
  integer(kind=4) :: indexJ ! index of reaction coordinate bin

  integer(kind=4) :: indexB ! index of combined Bin

  real(kind=fp_kind), allocatable :: logUnbiasedDensity(:)

  public :: startWHAM, finalizeWHAM
contains
  subroutine startWHAM(fdataid, foutid, temperature, tol, iBootstrap, nBootstrap)
    implicit none
    integer(kind=4), intent(in) :: fdataid
    integer(kind=4), intent(in) :: foutid
    real(kind=fp_kind), intent(in) :: temperature, tol
    integer(kind=4), intent(in) :: iBootstrap, nBootstrap
    idOutputFile = foutid
    T_target = temperature
    tolerance = tol
    NumBootstrap = nBootstrap
    read(fdataid,*)NumW, NumJ, NumK
    write(6,'(A,I6)')'Number of simulations:', NumW
    write(6,'(A,I6)')'Number of reaction coordinate bins:', NumJ
    write(6,'(A,I6)')'Number of energy bins:', NumK
    write(6,'(A,F10.3)')'Target temperature:', T_target
    write(6,'(A,F10.3)')'Convergence criterion:', tolerance
    if(debugLevel > 1) write(6,'(A)') 'call readReactCoordBinInfo'
    call readReactCoordBinInfo(fdataid, NumW, NumJ)
    nSimulation = NumW
    if(debugLevel > 1) write(6,'(A)') 'call readSimulationInfo'
    call readSimulationInfo(fdataid)
    call calculatePMF()
    if( iBootstrap == 1 ) then
      call bootstrap()
    else
      write(idOutputFile,'(2F10.4)')(reactCoordBin(indexB, 1)%binRC, logUnbiasedDensity(indexB), indexB = 1, NumB)
    end if
  end subroutine startWHAM

  subroutine calculatePMF()
    if(debugLevel > 1) write(6,'(A)') 'call initBins'
    call initBins(NumJ, NumK, NumB, T_target)

    if(allocated(logUnbiasedDensity))deallocate(logUnbiasedDensity)
    allocate(logUnbiasedDensity(NumB))

    if(debugLevel > 1) write(6,'(A)') 'call iteration'
    call iteration
  end subroutine calculatePMF

  subroutine bootstrap()
    use precision_m
    use simulation
    implicit none
    integer(kind=4) :: indexBootstrap
    integer(kind=4) :: iPointerSnapshot
    type(simulation_info), allocatable :: simulationsOrigin(:)
    real(kind=fp_kind) :: myrand
    real(kind=fp_kind), allocatable :: logUnbiasedDensityBootstrap(:,:)
    real(kind=fp_kind), allocatable :: expectValueBootstrap(:), GStandardErrBootstrap(:)
    real(kind=fp_kind), allocatable :: iLogUnbiasedDensityBootstrap(:)
    real(kind=fp_kind) :: expectValue, standardErr
    real(kind=fp_kind), allocatable :: logUnbiasedDensityOrigin(:)
    if(debugLevel>1)write(*,'(A)')'Begin bootstrapping'
    if(debugLevel == 1) debugLevel = 0
! copy simulations to simulationsOrigin
    allocate(simulationsOrigin(nSimulation))
    allocate(logUnbiasedDensityBootstrap(NumB,NumBootstrap))
    allocate(expectValueBootstrap(NumB))
    allocate(GStandardErrBootstrap(NumB))
    allocate(iLogUnbiasedDensityBootstrap(NumBootstrap))
    allocate(logUnbiasedDensityOrigin(NumB))
    logUnbiasedDensityOrigin = logUnbiasedDensity
    do indexW = 1, nSimulation
      simulationsOrigin(indexW)%nSnapshots = simulations(indexW)%nSnapshots
      simulationsOrigin(indexW)%beta = simulations(indexW)%beta
      allocate(simulationsOrigin(indexW)%snapshots(simulations(indexW)%nSnapshots))
      simulationsOrigin(indexW)%snapshots(:)%jReactCoordBin = &
            & simulations(indexW)%snapshots(:)%jReactCoordBin
      simulationsOrigin(indexW)%snapshots(:)%energyUnbiased = &
            & simulations(indexW)%snapshots(:)%energyUnbiased
    end do

    do indexBootstrap = 1, NumBootstrap

      do indexW = 1, nSimulation
        nullify(simulations(indexW)%snapshots)
        nullify(simulations(indexW)%bins)
      end do
      deallocate(simulations)
      allocate(simulations(nSimulation))

      do indexW = 1, nSimulation
        simulations(indexW)%nSnapshots = simulationsOrigin(indexW)%nSnapshots
        simulations(indexW)%beta = simulationsOrigin(indexW)%beta
        allocate(simulations(indexW)%snapshots(simulations(indexW)%nSnapshots))
        do indexS = 1, simulations(indexW)%nSnapshots
          iPointerSnapshot = int(myrand()*simulations(indexW)%nSnapshots) + 1
          simulations(indexW)%snapshots(indexS)%jReactCoordBin = &
              & simulationsOrigin(indexW)%snapshots(iPointerSnapshot)%jReactCoordBin
          simulations(indexW)%snapshots(indexS)%energyUnbiased = &
              & simulationsOrigin(indexW)%snapshots(iPointerSnapshot)%energyUnbiased
        end do
      end do   
      if(debugLevel > 1) write(*,'(A)') 'call calculatePMF'
      call calculatePMF()

      logUnbiasedDensityBootstrap(:,indexBootstrap) = logUnbiasedDensity(:)
    end do
    do indexB = 1, NumB
      iLogUnbiasedDensityBootstrap(:) = logUnbiasedDensityBootstrap(indexB,:)
      call Mean_and_StandardErr(NumBootstrap,iLogUnbiasedDensityBootstrap,expectValue,standardErr)
      expectValueBootstrap(indexB) = expectValue
      GStandardErrBootstrap(indexB) = standardErr
    end do
    write(idOutputFile, '(4F10.4)')(reactCoordBin(indexB, 1)%binRC, logUnbiasedDensityOrigin(indexB), &
            & expectValueBootstrap(indexB), GStandardErrBootstrap(indexB), indexB = 1, NumB)
    deallocate(simulationsOrigin)
    deallocate(logUnbiasedDensityBootstrap)
    deallocate(expectValueBootstrap)
    deallocate(GStandardErrBootstrap)
    deallocate(iLogUnbiasedDensityBootstrap)
    deallocate(logUnbiasedDensityOrigin)
  end subroutine bootstrap

  subroutine finalizeWHAM
    implicit none
    if(allocated(logUnbiasedDensity))deallocate(logUnbiasedDensity)
    call deleteReactCoordBinInfo
    call deleteSimulationInfo
  end subroutine finalizeWHAM

  subroutine iteration
    use constant, only : kB
    implicit none
    integer(kind=4) :: iIteration
    real(kind=fp_kind) :: deltaG(NumB)
    real(kind=fp_kind) :: aCap
    real(kind=fp_kind) :: dAcapdDeltaG(NumB)
    real(kind=fp_kind) :: eps, xtol
    integer(kind=4) :: iprint(2)
    integer(kind=4) :: iflag
    logical :: diagco
    real(kind=fp_kind) :: DIAG(NumW)
    integer(kind=4), parameter :: Mcorrection = 5
    real(kind=fp_kind) :: w((NumW+1)*(2*Mcorrection+1))
    external LB2
    diagco = .false.
    iprint(1) = 0
    iprint(2) = 0
    eps = 1.0D-4
    xtoL = 1.D-14
    iflag = 0
    deltaG = 0.d0
    do iIteration = 1, MaxITS
      if(debugLevel > 0) write(*,'(A, I)') 'Iteration', iIteration
      call AofDeltaG( deltaG, aCap, dAcapdDeltaG )
      if(debugLevel > 1) write(*,'(A, F10.4)') 'aCap=', aCap
      if(debugLevel > 1) write(*,'(8F10.3)') dAcapdDeltaG
      if(debugLevel > 1) then
        iprint(1) = 1
        iprint(2) = 0
      end if
      call LBFGS(NumW, Mcorrection, deltaG, aCap, dAcapdDeltaG, diagco, diag, iprint, eps, xtoL, w, iflag)
      if(iflag<=0)exit
      if(iIteration == MaxITS)write(*,'(A)') 'Max number of iterations has been reached'
    end do
    call deltaG2G(deltaG)
  end subroutine iteration

  subroutine AofDeltaG( deltaG, aCap, dAcapdDeltaG )
    implicit none
    real(kind=fp_kind), intent(in) :: deltaG(NumW)
    real(kind=fp_kind), intent(out) :: aCap, dAcapdDeltaG(NumW)
    integer(kind=4) :: indexB, indexW
    integer(kind=4) :: indexW2
    integer(kind=4) :: totalHistogram(NumB)
    real(kind=fp_kind) :: g(NumW)
    real(kind=fp_kind) :: denominator(NumB)
    real(kind=fp_kind) :: sumOfFraction

    g(1) = 0.d0
    do indexB = 2, NumB
      g(indexB) = g(indexB-1) + deltaG(indexB-1)
    end do 

    totalHistogram = 0 
    do indexB = 1, NumB
      do indexW = 1, NumW
        totalHistogram(indexB) = totalHistogram(indexB) & 
                             & + simulations(indexW)%bins(indexB)%histogram
      end do
    end do

    denominator = 0.d0
    do indexB = 1, NumB
      do indexW = 1, NumW
        denominator(indexB) = denominator(indexB) + &
                            & simulations(indexW)%nEffectivenSnapshots * &
                            & simulations(indexW)%bins(indexB)%biasingFactor * exp(g(indexW))
      end do
    end do

    aCap = 0.d0
    do indexW = 2, NumW
      aCap = aCap - simulations(indexW)%nEffectivenSnapshots * g(indexW)
    end do

    do indexB = 1, NumB
      aCap = aCap - totalHistogram(indexB) * &
                & log( totalHistogram(indexB)/denominator(indexB) )
    end do
    
    dAcapdDeltaG = 0.d0
    do indexW = 1, NumW
      do indexW2 = indexW + 1, NumW
        sumOfFraction = 0.d0
        do indexB = 1, NumB
          sumOfFraction = sumOfFraction + &
              & totalHistogram(indexB)*simulations(indexW2)%bins(indexB)%biasingFactor/denominator(indexB)
        end do
        dAcapdDeltaG(indexW) = dAcapdDeltaG(indexW) + &
          & simulations(indexW2)%nEffectivenSnapshots * ( exp(g(indexW2)) * sumOfFraction - 1.d0 )
      end do
    end do
  end subroutine AofDeltaG

  subroutine deltaG2G(deltaG)
    use constant, only : kB
    implicit none
    real(kind=fp_kind) :: deltaG(NumW)
    real(kind=fp_kind) :: g(NumW)
    real(kind=fp_kind) :: f(NumW)
    real(kind=fp_kind) :: unbiasedDensity(NumB)
    integer(kind=4) :: totalHistogram(NumB)
    real(kind=fp_kind) :: denominator(NumB)
    integer(kind=4) :: indexB, indexW
    

    g(1) = 0.d0
    do indexB = 2, NumB
      g(indexB) = g(indexB-1) + deltaG(indexB-1)
    end do 

    f = exp(g)

    totalHistogram = 0 
    do indexB = 1, NumB
      do indexW = 1, NumW
        totalHistogram(indexB) = totalHistogram(indexB) & 
                             & + simulations(indexW)%bins(indexB)%histogram
      end do
    end do

    do indexB = 1, NumB
      denominator(indexB) = 0.d0
      do indexW = 1, nSimulation
        denominator(indexB) = denominator(indexB) + &
           & simulations(indexW)%nEffectivenSnapshots * f(indexW) * &
           & simulations(indexW)%bins(indexB)%biasingFactor
      end do
    end do
    unbiasedDensity = totalHistogram / denominator

    logUnbiasedDensity = - kB * T_target * log(unbiasedDensity)
    logUnbiasedDensity = logUnbiasedDensity - minval(logUnbiasedDensity)
  end subroutine deltaG2G
end module WHAM
