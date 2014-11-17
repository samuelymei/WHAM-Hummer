!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for the calculation of the free energy from biased simulation like umbrella sampling !
! and temperature replica exchange molecular dynamics simualtion                              !
! Written by                                                                                  !
!                                        Ye Mei                                               !
!                                      10/17/2014                                             !
!                           East China Normal University                                      !
! Reference:                                                                                  !
!      J. Phys. Chem. B, 109, 6722-6731 (2005)                                                !
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
    real(kind=fp_kind) :: numerator(NumB), denominator(NumB)
    real(kind=fp_kind) :: freeenergyMin
    real(kind=fp_kind) :: unbiasedDensityRMSD
    real(kind=fp_kind) :: unbiasedDensityOld(NumB)
    real(kind=fp_kind) :: sumOfUnbiasedDensity
    real(kind=fp_kind) :: freeenergy(NumW)
    integer(kind=4) :: iIteration
    logical :: converged
    real(kind=fp_kind) :: rmsd
    real(kind=fp_kind) :: logUnbiasedDensity(NumB)
    allocate(unbiasedDensity(NumB))
    freeenergy = 1.d0   ! assign an initial guess of the free energy
    unbiasedDensity = 1.d0/NumB
    converged = .false.

    do indexB = 1, NumB
      numerator(indexB) = 0.d0
      do indexW = 1, nSimulation
        numerator(indexB) = numerator(indexB) + & 
           & simulations(indexW)%bins(indexB)%histogram
      end do
    end do

    do iIteration = 1, MAXITS
       unbiasedDensityOld = unbiasedDensity
      do indexB = 1, NumB
        denominator(indexB) = 0.d0
        do indexW = 1, nSimulation
          denominator(indexB) = denominator(indexB) + &
             & simulations(indexW)%nEffectivenSnapshots * freeenergy(indexW) * &
             & simulations(indexW)%bins(indexB)%biasingFactor
        end do
      end do
      unbiasedDensity = numerator / denominator
 
      sumOfUnbiasedDensity = sum(unbiasedDensity)
      unbiasedDensity = unbiasedDensity / sumOfUnbiasedDensity

      freeenergy = 0.d0
      do indexW = 1, nSimulation
        do indexB = 1, NumB
          freeenergy(indexW) = freeenergy(indexW) + &
            & simulations(indexW)%bins(indexB)%biasingFactor * unbiasedDensity(indexB)
        end do    
      end do
      freeenergy = 1.d0 / freeenergy

      unbiasedDensityRMSD = rmsd(NumB, unbiasedDensity, unbiasedDensityOld)
      write(6,'(A,1X,I4,A,E12.4)')'Iteration ', iIteration, ': RMSD of unbiased density:', unbiasedDensityRMSD

      if( unbiasedDensityRMSD < TOLERANCE ) converged = .true.
      if( converged ) then
        logUnbiasedDensity = - kB * T_target * log(unbiasedDensity)
        logUnbiasedDensity = logUnbiasedDensity - minval(logUnbiasedDensity)
        write(99,'(2F10.4)')(reactCoordBin(indexB, 1)%binRC, logUnbiasedDensity(indexB), indexB = 1, NumB)
        exit
      end if
      if ( iIteration == MAXITS ) then
         write(6,*)'Convergence failure'
      end if
    end do
  end subroutine iteration
end module WHAM
