!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for the calculation of the free energy from biased simulation like umbrella sampling !
! and temperature replica exchange molecular dynamics simualtion                              !
! Written by                                                                                  !
!                                        Ye Mei                                               !
!                           East China Normal University                                      !
!                                      11/17/2014                                             !
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
 
  integer(kind=4), allocatable :: totalHistogram(:)

  public :: startWHAM, finalizeWHAM
  public :: getDaDg, aCap
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
    write(6,'(A,E10.3)')'Convergence criterion:', tolerance
    if(debugLevel > 1) write(6,'(A)') 'call readReactCoordBinInfo'
    call readReactCoordBinInfo(fdataid, NumW, NumJ)
    nSimulation = NumW
    if(debugLevel > 1) write(6,'(A)') 'call readSimulationInfo'
    call readSimulationInfo(fdataid)
    call calculatePMF()
    print*,'iBootstrap=',iBootstrap
    if( iBootstrap == 1 ) then
      call bootstrap()
    else
      write(idOutputFile,'(2F10.4)')(reactCoordBin(indexJ, 1)%binRC, logUnbiasedDensity(indexJ), indexJ = 1, NumJ)
    end if
  end subroutine startWHAM

  subroutine calculatePMF()
    if(debugLevel > 1) write(6,'(A)') 'call initBins'
    call initBins(NumJ, NumK, NumB, T_target)

    if(allocated(logUnbiasedDensity))deallocate(logUnbiasedDensity)
    allocate(logUnbiasedDensity(NumJ))

    if(debugLevel > 1) write(6,'(A)') 'call iteration'
    call calculateTotalHistogram
    call iteration
!    call iteration2
  end subroutine calculatePMF

  subroutine bootstrap()
    implicit none
    integer(kind=4) :: indexBootstrap
    integer(kind=4) :: iPointerSnapshot
    type(simulation_info), allocatable :: simulationsOrigin(:)
    real(kind=fp_kind) :: myrand
    real(kind=fp_kind), allocatable :: logUnbiasedDensityBootstrap(:,:)
    real(kind=fp_kind), allocatable :: expectValueBootstrap(:), GStandardDevBootstrap(:)
    real(kind=fp_kind), allocatable :: iLogUnbiasedDensityBootstrap(:)
    real(kind=fp_kind) :: expectValue, standardErr, standardDev
    real(kind=fp_kind), allocatable :: logUnbiasedDensityOrigin(:)
    if(debugLevel > 1)write(*,'(A)')'Begin bootstrapping'
    if(debugLevel == 1) debugLevel = 0
! copy simulations to simulationsOrigin
    allocate(simulationsOrigin(nSimulation))
    allocate(logUnbiasedDensityBootstrap(NumJ,NumBootstrap))
    allocate(expectValueBootstrap(NumJ))
    allocate(GStandardDevBootstrap(NumJ))
    allocate(iLogUnbiasedDensityBootstrap(NumBootstrap))
    allocate(logUnbiasedDensityOrigin(NumJ))
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
      write(*,'(A,I10)')'Bootstrap step:', indexBootstrap

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
    do indexJ = 1, NumJ
      iLogUnbiasedDensityBootstrap(:) = logUnbiasedDensityBootstrap(indexJ,:)
!      call Mean_and_StandardErr(NumBootstrap,iLogUnbiasedDensityBootstrap,expectValue,standardErr)
      call Mean_and_StandardDev(NumBootstrap,iLogUnbiasedDensityBootstrap,expectValue,standardDev)
      expectValueBootstrap(indexJ) = expectValue
      GStandardDevBootstrap(indexJ) = standardDev
    end do
    write(idOutputFile, '(4F10.4)')(reactCoordBin(indexJ, 1)%binRC, logUnbiasedDensityOrigin(indexJ), &
            & expectValueBootstrap(indexJ), GStandardDevBootstrap(indexJ), indexJ = 1, NumJ)
    deallocate(simulationsOrigin)
    deallocate(logUnbiasedDensityBootstrap)
    deallocate(expectValueBootstrap)
    deallocate(GStandardDevBootstrap)
    deallocate(iLogUnbiasedDensityBootstrap)
    deallocate(logUnbiasedDensityOrigin)
  end subroutine bootstrap

  subroutine finalizeWHAM
    implicit none
    if(allocated(logUnbiasedDensity))deallocate(logUnbiasedDensity)
    if(allocated(totalHistogram))deallocate(totalHistogram)
    call deleteReactCoordBinInfo
    call deleteSimulationInfo
    write(*,'(A)')'WHAM Done'
  end subroutine finalizeWHAM

  subroutine iteration2
    use constant, only : kB
    implicit none
    integer(kind=4) :: iIteration
    real(kind=fp_kind) :: deltaG(NumW-1)
    real(kind=fp_kind) :: dAcapdDeltaG(NumW-1)
    integer(kind=4) :: iter
    real(kind=fp_kind) :: fret
    real(kind=fp_kind) :: gtol
    deltaG = -1.0d0
    gtol=1.D-8
    fret = 0.d0
    if(debugLevel > 1) write(*,'(A, E20.8)') 'aCap=', aCap(deltaG)
    call getDadG(deltaG,dAcapdDeltaG)
    if(debugLevel > 1) write(*,'(8E14.6)') dAcapdDeltaG
    call dfpmin(deltaG,NumW-1,gtol,iter,fret,aCap,getDaDg)
    if(debugLevel > 1) write(*,'(A,I4,A)') 'dfpmin converged in ',iter,' cycles'
    call deltaG2G(deltaG)
  end subroutine iteration2

  subroutine iteration
    use constant, only : kB
    implicit none
    integer(kind=4) :: iIteration
    real(kind=fp_kind) :: deltaG(NumW-1)
    real(kind=fp_kind) :: dAcapdDeltaG(NumW-1)
    real(kind=fp_kind) :: eps, xtol
    integer(kind=4) :: iprint(2)
    integer(kind=4) :: iflag
    logical :: diagco
    real(kind=fp_kind) :: DIAG(NumW-1)
    integer(kind=4), parameter :: Mcorrection = 7
    real(kind=fp_kind) :: w((NumW+1)*(2*Mcorrection+1))
    external LB2
    diagco = .false.
    iprint(1) = 0
    iprint(2) = 0
    eps = 1.0D-8
    xtoL = 1.D-13
    iflag = 0
    deltaG = -1.0d0
    do iIteration = 1, MaxITS
      if(debugLevel > 0) write(*,'(A, I5)') 'Iteration', iIteration
      if(debugLevel > 1) write(*,'(A, E20.8)') 'aCap=', aCap(deltaG)
      call getDadG(deltaG,dAcapdDeltaG)
      if(debugLevel > 1) write(*,'(8E14.6)') dAcapdDeltaG
      if(debugLevel > 1) then
        iprint(1) = 1
        iprint(2) = 0
      end if
      call LBFGS(NumW-1, Mcorrection, deltaG, aCap, dAcapdDeltaG, diagco, diag, iprint, eps, xtoL, w, iflag)
      if(iflag<=0)exit
      if(iIteration == MaxITS)write(*,'(A)') 'Max number of iterations has been reached'
    end do
    call deltaG2G(deltaG)
  end subroutine iteration

  function aCap(deltaG)
    implicit none
    real(kind=fp_kind), intent(in) :: deltaG(NumW-1)
    real(kind=fp_kind) :: aCap
    integer(kind=4) :: indexB, indexW
    integer(kind=4) :: indexW2
    real(kind=fp_kind) :: g(NumW)
    real(kind=fp_kind) :: denominator(NumB)
    g(1) = 0.d0
    do indexW = 2, NumW
      g(indexW) = g(indexW-1) + deltaG(indexW-1)
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

  end function aCap

  subroutine getDaDg(deltaG,dAdG)
    implicit none
    real(kind=fp_kind), intent(in) :: deltaG(NumW-1)
    real(kind=fp_kind), intent(out) :: dAdG(NumW-1)
    integer(kind=4) :: indexB, indexW
    integer(kind=4) :: indexW2
    real(kind=fp_kind) :: g(NumW)
    real(kind=fp_kind) :: denominator(NumB)
    real(kind=fp_kind) :: sumOfFraction

    g(1) = 0.d0
    do indexW = 2, NumW
      g(indexW) = g(indexW-1) + deltaG(indexW-1)
    end do 

    denominator = 0.d0
    do indexB = 1, NumB
      do indexW = 1, NumW
        denominator(indexB) = denominator(indexB) + &
                            & simulations(indexW)%nEffectivenSnapshots * &
                            & simulations(indexW)%bins(indexB)%biasingFactor * exp(g(indexW))
      end do
    end do

    dAdG = 0.d0
    do indexW = 1, NumW
      do indexW2 = indexW + 1, NumW
        sumOfFraction = 0.d0
        do indexB = 1, NumB
          sumOfFraction = sumOfFraction + &
              & totalHistogram(indexB)*simulations(indexW2)%bins(indexB)%biasingFactor/denominator(indexB)
        end do
        dAdG(indexW) = dAdG(indexW) + &
          & simulations(indexW2)%nEffectivenSnapshots * ( exp(g(indexW2)) * sumOfFraction - 1.d0 )
      end do
    end do
  end subroutine getDaDg

  subroutine deltaG2G(deltaG)
    use constant, only : kB
    implicit none
    real(kind=fp_kind) :: deltaG(NumW-1)
    real(kind=fp_kind) :: g(NumW)
    real(kind=fp_kind) :: f(NumW)
    real(kind=fp_kind) :: unbiasedDensity(NumB)
    real(kind=fp_kind) :: denominator(NumB)
    integer(kind=4) :: indexB, indexW
    if(debugLevel > 1) write(*,'(A)')'call deltaG2G'

    g(1) = 0.d0
    do indexW = 2, NumW
      g(indexW) = g(indexW-1) + deltaG(indexW-1)
    end do 

    f = exp(g)

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

  subroutine calculateTotalHistogram
    implicit none
    if(allocated(totalHistogram))deallocate(totalHistogram)
    allocate(totalHistogram(NumB))
    totalHistogram = 1
    do indexB = 1, NumB
      do indexW = 1, NumW
        totalHistogram(indexB) = totalHistogram(indexB) &
                             & + simulations(indexW)%bins(indexB)%histogram
      end do
    end do
  end subroutine calculateTotalHistogram
end module WHAM
