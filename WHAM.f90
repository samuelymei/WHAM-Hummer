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
  use simulation
  use react_coord_bin
  implicit none
  
  private
  integer(kind=4), public :: NumW  ! number of simulations
  integer(kind=4), public :: NumJ  ! number of reaction coordinate bins
  integer(kind=4), public :: NumK  ! number of energy bins
  integer(kind=4), public :: NumB  ! number of bins
  real(kind=fp_kind), parameter :: TOLERANCE = 1.0D-4  ! convergence criterion for free energy 
  integer(kind=4), parameter :: MaxITS = 1000          ! max number of iterations
  real(kind=fp_kind), public :: T_target = 3.0D2      ! target temperature

  real(kind=fp_kind) :: beta ! inverse temperature

  integer(kind=4) :: indexW ! index of windows/simulations
  integer(kind=4) :: indexS ! index of snapshot
  integer(kind=4) :: indexK ! index of energy bin
  integer(kind=4) :: indexJ ! index of reaction coordinate bin

  integer(kind=4) :: indexB ! index of combined Bin

  real(kind=fp_kind), allocatable :: unbiasedDensity(:)
  
  public :: startWHAM, finalizeWHAM
contains
  subroutine startWHAM(fid)
    implicit none
    integer(kind=4), intent(in) :: fid
    read(fid,*)NumW, NumJ, NumK, T_target
    write(6,'(A,I6)')'Number of simulations:', NumW
    write(6,'(A,I6)')'Number of reaction coordinate bins:', NumJ
    write(6,'(A,I6)')'Number of energy bins:', NumK
    write(6,'(A,F10.3)')'Target temperature:', T_target
    print*,'call readReactCoordBinInfo'
    call readReactCoordBinInfo(fid, NumW, NumJ)
    nSimulation = NumW
    print*,'call readSimulationInfo'
    call readSimulationInfo(fid)
    print*,'call initBins'
    call initBins(NumJ, NumK, NumB, T_target)
    print*,'iteration'
    call iteration
  end subroutine startWHAM

  subroutine finalizeWHAM
    implicit none
    if(allocated(unbiasedDensity))deallocate(unbiasedDensity)
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
    iprint(1) = 1
    iprint(2) = 0
    eps = 1.0D-4
    xtoL = 1.D-14
    iflag = 0
    deltaG = 0.d0
    do iIteration = 1, MaxITS
      print*,'Iteration', iIteration
      call AofDeltaG( deltaG, aCap, dAcapdDeltaG )
      print*, 'aCap=', aCap
      print*, dAcapdDeltaG
      call LBFGS(NumW, Mcorrection, deltaG, aCap, dAcapdDeltaG, diagco, diag, iprint, eps, xtoL, w, iflag)
      if(iflag<=0)exit
      if(iIteration == MaxITS)print*,'Max number of iterations has been reached'
    end do
    call writedensity(deltaG)
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

  subroutine writedensity(deltaG)
    use constant, only : kB
    implicit none
    real(kind=fp_kind) :: deltaG(NumW)
    real(kind=fp_kind) :: g(NumW)
    real(kind=fp_kind) :: f(NumW)
    real(kind=fp_kind) :: unbiasedDensity(NumB)
    real(kind=fp_kind) :: logUnbiasedDensity(NumB)
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
    write(99,'(2F10.4)')(reactCoordBin(indexB, 1)%binRC, logUnbiasedDensity(indexB), indexB = 1, NumB)

  end subroutine writedensity
end module WHAM
