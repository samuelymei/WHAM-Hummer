program WHAM_caller
  use precision_m
  use control
  use WHAM
  implicit none
  integer(kind=4) :: narg, iarg
  character(len=20) :: arg
  integer(kind=4) :: fmetaid, foutputid
  character(len=60) :: metafile, outputfile
  real(kind=fp_kind) :: tolerance
  real(kind=fp_kind) :: temperature
  character(len=200) :: command_help
  integer(kind=4) :: iBootstrap
  character(len=200):: thisProgram
  character(len=255) :: cmd
  logical :: enoughArg
  call initialControlOpt(temperature,tolerance,iBootstrap,metafile,outputfile)
  call get_command(cmd)
  call get_command_argument(0, thisProgram)
  command_help = 'Usage: '//trim(thisProgram)//" [-s] [-d debugLevel] -T temperature -t tolerance -f metafile -o outputfile -h"
  narg = command_argument_count()
  iarg = 1
  iBootstrap = 0
  debugLevel = 1
  do while (iarg <= narg)
    call get_command_argument(iarg, arg)
    if ( arg == '-s' .or. arg == '-S' ) then
      iBootstrap = 1
    else if ( arg == '-d' ) then
      if( .not. enoughArg(narg, iarg+1) ) then
        write(*,'(A)') command_help
        stop
      end if
      iarg = iarg + 1
      call get_command_argument(iarg, arg)
      read(arg, *) debugLevel
    else if( arg == '-T' ) then
      if( .not. enoughArg(narg, iarg+1) ) then
        write(*,'(A)') command_help
        stop
      end if
      iarg = iarg + 1
      call get_command_argument(iarg, arg)
      read(arg,*) temperature
    else if( arg == '-t' ) then
      if( .not. enoughArg(narg, iarg+1) ) then
        write(*,'(A)') command_help
        stop
      end if
      iarg = iarg + 1
      call get_command_argument(iarg, arg)
      read(arg,*) tolerance
    else if( arg == '-f' ) then
      if( .not. enoughArg(narg, iarg+1) ) then
        write(*,'(A)') command_help
        stop
      end if
      iarg = iarg + 1
      call get_command_argument(iarg, arg)
      metafile = arg
    else if( arg == '-o' ) then
      if( .not. enoughArg(narg, iarg+1) ) then
        write(*,'(A)') command_help
        stop
      end if
      iarg = iarg + 1
      call get_command_argument(iarg, arg)
      outputfile = arg
    else if( arg == '-h' .or. arg == '-H' ) then
      write(*,'(A)') command_help
      stop
    else
      write(*,'(A)') command_help
      stop
    end if
    iarg = iarg + 1
  end do
  write(*,*)trim(cmd)
  write(*,'(A,F10.3)')'Target temperature:', temperature
  write(*,'(A,E10.3)')'Tolerance:', tolerance
  write(*,'(A,A)') 'Metafile: ', metafile
  write(*,*)
  fmetaid = 10
  foutputid = 11
  open(fmetaid, file = metafile, status = 'OLD')
  open(foutputid, file = outputfile)
  call startWHAM( fmetaid, foutputid, tolerance, iBootstrap)
  call finalizeWHAM
  close(fmetaid)
  close(foutputid)
end program WHAM_caller

function enoughArg( narg, nargInNeed )
  implicit none
  logical :: enoughArg
  integer(kind=4), intent(in) :: narg, nargInNeed
  enoughArg = .true.
  if( narg < nargInNeed ) enoughArg = .false.
  return
end function enoughArg

subroutine initialControlOpt( temperature, tolerance, iBootstrap, metafile, outputfile)
  use precision_m
  implicit none
  integer(kind=4), intent(out) :: iBootstrap
  real(kind=fp_kind), intent(out) :: temperature, tolerance
  character(len=*), intent(out) :: metafile, outputfile
  temperature = 300.d0
  tolerance = 1.D-4
  iBootstrap = 0
  metafile = 'metafile'
  outputfile = 'wham.out'
end subroutine initialControlOpt
