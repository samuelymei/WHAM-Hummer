function isInRange(i, imin, imax)
  implicit none
  logical :: isInRange
  integer(kind=4), intent(in) :: i, imin, imax
  isInRange = .false.
  if(i<=imax .and. i>=imin) isInRange = .true.
end function isInRange

function rmsd(n,a1,a2)
  use precision_m
  implicit none
  integer(kind=4), intent(in) :: n
  real(kind=fp_kind), intent(in) :: a1(n), a2(n)
  real(kind=fp_kind) :: rmsd
  integer(kind=4) :: i
  if(n==0)then
    rmsd = 0.d0
    return
  end if
  rmsd = 0.d0
  rmsd = sqrt(sum((a1-a2)**2)/n)
end function rmsd

