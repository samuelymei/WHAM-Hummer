module simulation
  use precision_m
  use snapshot
  use react_coord_bin
  use bin
  private
  integer(kind=4), public :: nSimulation

  type :: simulation_info
    real(kind=fp_kind) :: beta  ! inverse temperature of this simulation
    integer(kind=4) :: nSnapshots ! Number of snapshots in this simulation
    integer(kind=4) :: nEffectivenSnapshots ! snapshots outside the RC bins are removed
    type ( snapshot_info ), pointer :: snapshots(:) ! point to the snapshots
    type ( bin_info ), pointer :: bins(:)
  end type simulation_info
  type ( simulation_info ), allocatable, public :: simulations(:) ! to save the information of all the simulations
  real(kind=fp_kind) :: energyMin, energyMax
  public :: readSimulationInfo, deleteSimulationInfo, initBins

contains

  subroutine readSimulationInfo(fid)
    use constant, only : kB
    implicit none
    integer(kind=4), intent(in) :: fid
    real(kind=fp_kind) :: simulationTemperature
    integer(kind=4) :: indexW, indexS
    allocate(simulations(nSimulation))
    do indexW = 1, nSimulation
      read(fid, *) simulationTemperature, simulations(indexW)%nSnapshots
      simulations(indexW)%beta = 1.d0 / ( kB * simulationTemperature )
      write(6,'(A,I3,A,F10.6,A,I10)') ' Simulation ', indexW, ':    beta = ', simulations(indexW)%beta, &
         & ', Number of snapshots:', simulations(indexW)%nSnapshots
      allocate(simulations(indexW)%snapshots(simulations(indexW)%nSnapshots))
      do indexS = 1, simulations(indexW)%nSnapshots
        read(fid,*) simulations(indexW)%snapshots(indexS)%jReactCoordBin, & 
                  & simulations(indexW)%snapshots(indexS)%energyUnbiased
      end do
    end do
    call energyMinMax
  end subroutine readSimulationInfo

  subroutine deleteSimulationInfo
    implicit none
    integer(kind=4) :: indexW
    do indexW = 1, nSimulation
      deallocate(simulations(indexW)%snapshots)
    end do
    if(allocated(simulations))deallocate(simulations)
  end subroutine deleteSimulationInfo

  subroutine energyMinMax
    implicit none
    integer(kind=4) :: indexW, indexS
    energyMin = simulations(1)%snapshots(1)%energyUnbiased
    energyMax = simulations(1)%snapshots(1)%energyUnbiased
    do indexW = 1, nSimulation
      do indexS = 1, simulations(indexW)%nSnapshots
        if(simulations(indexW)%snapshots(indexS)%energyUnbiased < energyMin) & 
           & energyMin = simulations(indexW)%snapshots(indexS)%energyUnbiased
        if(simulations(indexW)%snapshots(indexS)%energyUnbiased > energyMax) & 
           & energyMax = simulations(indexW)%snapshots(indexS)%energyUnbiased
      end do
    end do
    write(6,'(A,G12.5,A,G12.5)')'Energy range: ', energyMin, '~', energyMax
  end subroutine energyMinMax

  subroutine initBins(NumRCbin, NumEbin, NumBin, T_target)
    use constant
    implicit none
    integer(kind=4), intent(in) :: NumRCbin, NumEbin
    integer(kind=4), intent(out) :: NumBin
    real(kind=fp_kind), intent(in) :: T_target

    real(kind=fp_kind) :: energy(NumEbin)
    real(kind=fp_kind) :: deltaE
    integer(kind=4) :: indexRCbin, indexEbin
    real(kind=fp_kind) :: delta 
    integer(kind=4) :: indexW, indexB
    real(kind=fp_kind) :: beta_target
    real(kind=fp_kind) :: biasingPotential

    delta = (energyMax - energyMin)/(100*NumEbin)
    deltaE = (energyMax - energyMin + delta)/NumEbin
    do indexEbin = 1, NumEbin
      energy(indexEbin) = energyMin + deltaE * (indexEbin-0.5)
    end do

    NumBin = NumRCbin * NumEbin
    write(6,*)'     NumRCbin      NumEbin      NumBin     nSimulation '
    write(6,*)NumRCbin, NumEbin, NumBin, nSimulation
    do indexW = 1, nSimulation
      allocate(simulations(indexW)%bins(NumBin))
    end do

    beta_target = 1.d0 / (kB*T_target)
    write(6,'(A,F10.6)') 'Beta of the target temperature:', beta_target

    indexB = 0
    do indexRCbin = 1, NumRCbin
      do indexEbin = 1, NumEbin
        indexB = indexB + 1 
        do indexW = 1, nSimulation
          simulations(indexW)%bins(indexB)%energyBiasing = reactCoordBin(indexRCbin,indexW)%energyBiasing
          simulations(indexW)%bins(indexB)%energy = energy(indexEbin)
          simulations(indexW)%bins(indexB)%energyBinWidth = deltaE
          simulations(indexW)%bins(indexB)%biasingFactor = &
            & exp(-(simulations(indexW)%beta-beta_target)*simulations(indexW)%bins(indexB)%energy) * &
            & exp(-simulations(indexW)%beta*simulations(indexW)%bins(indexB)%energyBiasing)
          simulations(indexW)%bins(indexB)%jReactCoord = indexRCbin
          simulations(indexW)%bins(indexB)%kEnergy = indexEbin
        end do
      end do
    end do
    print*,'call buildHistogram'
    call buildHistogram(NumRCbin, NumEbin, NumBin)
  end subroutine initBins

  subroutine buildHistogram(NumRCbin, NumEbin, NumBin)
    implicit none
    integer(kind=4), intent(in) :: NumRCbin, NumEbin, NumBin
    integer(kind=4) :: indexW, indexB, indexS
    integer(kind=4) :: indexRCbin, indexEbin
    integer(kind=4) :: iTemp
    logical :: isInRange
    do indexW = 1, nSimulation
      simulations(indexW)%bins(:)%histogram = 0
      do indexS = 1, simulations(indexW)%nSnapshots
        if(simulations(indexW)%snapshots(indexS)%jReactCoordBin < 0 .or. &
         & simulations(indexW)%snapshots(indexS)%jReactCoordBin > NumRCbin ) cycle

        if(NumEbin == 1) then  ! NOT T-WHAM
          indexEbin = 1
        else
          indexEbin = -999
          do iTemp = 1, NumEbin
            if ( abs( simulations(indexW)%snapshots(indexS)%energyUnbiased - &
                    & simulations(indexW)%bins(iTemp)%energy ) <= &
                      simulations(indexW)%bins(iTemp)%energyBinWidth / 2.d0 ) then
              indexEbin = iTemp
              cycle
            end if
          end do
        end if

        indexRCbin = simulations(indexW)%snapshots(indexS)%jReactCoordBin
        if(isInRange(indexEbin,1,NumEbin) .and. isInRange(indexRCbin,1,NumRCbin)) then
          indexB = ( indexRCbin - 1 ) * NumEbin + indexEbin
          simulations(indexW)%bins(indexB)%histogram = &
                     & simulations(indexW)%bins(indexB)%histogram + 1
        end if
      end do
      simulations(indexW)%nEffectivenSnapshots = sum(simulations(indexW)%bins(:)%histogram)
      write(6,'(A,I4,A,I10)')'Effective number of snapshots in simulation', indexW, ':', &
        & simulations(indexW)%nEffectivenSnapshots
    end do
  end subroutine buildHistogram

end module simulation

