module snapshot
  use precision_m
  implicit none
  public 
  integer(kind=4) :: nSnapshots
  type :: snapshot_info
    real(kind=fp_kind) :: energyUnbiased ! Energy
    integer(kind=4) :: jReactCoordBin    ! index of RC bin
    integer(kind=4) :: kEnergyBin        ! index of Energy bin
  end type snapshot_info
end module snapshot


