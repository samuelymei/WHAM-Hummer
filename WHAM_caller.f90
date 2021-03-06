program WHAM_caller
  use precision_m
  use control
  use WHAM
  implicit none
  integer(kind=4) :: narg, iarg
  character(len=20) :: arg
  integer(kind=4) :: fmetaid, foutputid
  character(len=60) :: metafile, outputfile
  integer(kind=4) :: nBootstrap
  real(kind=fp_kind) :: tolerance
  real(kind=fp_kind) :: temperature
  character(len=400) :: command_help
  integer(kind=4) :: iBootstrap
  character(len=200):: thisProgram
  character(len=255) :: cmd
  integer(kind=4) :: istatus
  logical :: enoughArg
  call initialControlOpt(temperature, tolerance, iBootstrap, nBootstrap, metafile, outputfile)
  call get_command(cmd)
  call get_command_argument(0, thisProgram)
  command_help = 'Usage: '//trim(thisProgram)//' [-s nBootstrap] [-d debugLevel] [-T temperature] [-t tolerance] [-f metafile] [-o outputfile] [-h]'&
              &  //NEW_LINE('A')//'  Default:'&
              &  //NEW_LINE('A')//'    nBootstrap = 0' &
              &  //NEW_LINE('A')//'    debugLevel = 1' &
              &  //NEW_LINE('A')//'    temperature = 300.0' &
              &  //NEW_LINE('A')//'    tolerance = 1.D-3' &
              &  //NEW_LINE('A')//'    metafile = metafile' &
              &  //NEW_LINE('A')//'    outputfile = wham.out'


  narg = command_argument_count()
  iarg = 1
  iBootstrap = 0
  debugLevel = 1
  do while (iarg <= narg)
    call get_command_argument(iarg, arg)
    if ( arg == '-s' .or. arg == '-S' ) then
      if( .not. enoughArg(narg, iarg+1) ) call quitme(command_help)
      call getnextarg( iarg, arg, istatus )
      if( istatus /= 0 ) call quitme(command_help)
      read(arg, *) nBootstrap
      if(nBootstrap > 0) iBootstrap = 1

    else if ( arg == '-d' ) then
      if( .not. enoughArg(narg, iarg+1) ) call quitme(command_help)
      call getnextarg( iarg, arg, istatus )
      if( istatus /= 0 ) call quitme(command_help)
      read(arg, *) debugLevel

    else if( arg == '-T' ) then
      if( .not. enoughArg(narg, iarg+1) ) call quitme(command_help)
      call getnextarg( iarg, arg, istatus )
      if( istatus /= 0 ) call quitme(command_help)
      read(arg,*) temperature

    else if( arg == '-t' ) then
      if( .not. enoughArg(narg, iarg+1) ) call quitme(command_help)
      call getnextarg( iarg, arg, istatus )
      if( istatus /= 0 ) call quitme(command_help)
      read(arg,*) tolerance

    else if( arg == '-f' ) then
      if( .not. enoughArg(narg, iarg+1) ) call quitme(command_help)
      call getnextarg( iarg, arg, istatus )
      if( istatus /= 0 ) call quitme(command_help)
      metafile = trim(arg)

    else if( arg == '-o' ) then
      if( .not. enoughArg(narg, iarg+1) ) call quitme(command_help)
      call getnextarg( iarg, arg, istatus )
      if( istatus /= 0 ) call quitme(command_help)
      outputfile = trim(arg)

    else if( arg == '-h' .or. arg == '-H' ) then
      call quitme(command_help)

    else
      call quitme(command_help)
    end if
    iarg = iarg + 1
  end do
  write(*,*)trim(cmd)
  write(*,'(A,F10.3)')'Target temperature:', temperature
  write(*,'(A,E10.3)')'Tolerance:', tolerance
  write(*,'(A,A)') 'Metafile: ', metafile
  if(iBootstrap == 1) then
    write(*,'(A,1X,I5,A)') 'Do bootstrap:', nBootstrap, ' Times'
  else
    write(*,'(A)') 'Do bootstrap: No'
  end if
  write(*,*)
  fmetaid = 10
  foutputid = 11
  open(fmetaid, file = metafile, status = 'OLD')
  open(foutputid, file = outputfile)
  call startWHAM( fmetaid, foutputid, temperature, tolerance, iBootstrap, nBootstrap)
  call finalizeWHAM
  close(fmetaid)
  close(foutputid)
end program WHAM_caller

subroutine getnextarg( iarg, arg, istatus )
  implicit none
  integer(kind=4), intent(in out) :: iarg
  character(len=*), intent(out) :: arg
  integer(kind=4), intent(out) :: istatus
  iarg = iarg + 1
  call get_command_argument( iarg, arg, status=istatus )
end subroutine getnextarg

function enoughArg( narg, nargInNeed )
  implicit none
  logical :: enoughArg
  integer(kind=4), intent(in) :: narg, nargInNeed
  enoughArg = .true.
  if( narg < nargInNeed ) enoughArg = .false.
  return
end function enoughArg

subroutine quitme(commandline)
  implicit none
  character(len=*) :: commandline
  write(*,'(A)') commandline
  stop
end subroutine quitme

subroutine initialControlOpt( temperature, tolerance, iBootstrap, nBootstrap, metafile, outputfile)
  use precision_m
  implicit none
  integer(kind=4), intent(out) :: iBootstrap, nBootstrap
  real(kind=fp_kind), intent(out) :: temperature, tolerance
  character(len=*), intent(out) :: metafile, outputfile
  temperature = 300.d0
  tolerance = 1.D-3
  iBootstrap = 0
  nBootstrap = 0
  metafile = 'metafile'
  outputfile = 'wham.out'
end subroutine initialControlOpt
