module react_coord_bin
  use precision_m
  implicit none
  private
  type :: reactCoordBin_info
    real(kind=fp_kind) :: binRC
    real(kind=fp_kind) :: energyBiasing
  end type reactCoordBin_info
  type ( reactCoordBin_info ), allocatable, public :: reactCoordBin(:, :)
  public :: readReactCoordBinInfo, deleteReactCoordBinInfo
contains
  subroutine readReactCoordBinInfo( fid, NumW, NumJ )
    implicit none
    integer(kind=4), intent(in) :: fid, NumW, NumJ
    integer(kind=4) :: indexW
    allocate(reactCoordBin(NumJ, NumW))
    read(fid,*)reactCoordBin(:, 1)%binRC
    do indexW = 2, NumW
      reactCoordBin(:, indexW)%binRC = reactCoordBin(:, 1)%binRC
    end do
    do indexW = 1, NumW
      read(fid,*)reactCoordBin(:, indexW)%energyBiasing
    end do
  end subroutine readReactCoordBinInfo

  subroutine deleteReactCoordBinInfo
    implicit none
    deallocate(reactCoordBin)
  end subroutine deleteReactCoordBinInfo
end module react_coord_bin
