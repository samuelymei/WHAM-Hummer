module bin
  use precision_m
  implicit none
  type :: bin_info
    real(kind=fp_kind) :: energy
    real(kind=fp_kind) :: energyBinWidth
    real(kind=fp_kind) :: energyBiasing
    real(kind=fp_kind) :: biasingFactor
    integer(kind=4) :: histogram
    integer(kind=4) :: jReactCoord
    integer(kind=4) :: kEnergy
  end type bin_info
end module bin
