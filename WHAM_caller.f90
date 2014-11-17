program WHAM_caller
  use WHAM
  implicit none
  integer(kind=4) :: fid
  character(len=60) :: datafile
  fid = 10
  read*,datafile
  open(fid, file = datafile)
  call startWHAM(fid)
  call finalizeWHAM
  close(fid)
end program WHAM_caller
