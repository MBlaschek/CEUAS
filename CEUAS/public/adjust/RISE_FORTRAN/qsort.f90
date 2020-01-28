
module sort
contains
! ***********************************
! *
  Subroutine Qsort(X, Ipt)
! *
! ***********************************
! * Sort Array X(:) in ascendent order 
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt.
! ***********************************
 
    Type Limits
       Integer :: Ileft, Iright
    End Type Limits
 
    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer, Parameter :: Isw = 10
 
    Real (kind=8), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
 
    Integer :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Integer, Allocatable :: IIpt(:)
    Type (Limits), Allocatable :: Stack(:)
 
 
    Allocate(Stack(Size(X)))
 
    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I
 
       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
 
       Do While (Stack(ISpos)%Ileft /= 0)
 
          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn, Ipt)
 
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do
 
    Else
 
       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
 
       Allocate(IIpt(10))
       Do While (Stack(ISpos)%Ileft /= 0)
!          Write(*,*)Ispos, ISmax
 
          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, IIpt, Ileft, Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn)
 
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do
       Deallocate(IIpt)
 
    End If
 
    Deallocate(Stack)
 
    Return
 
  CONTAINS
 
    ! ***********************************
    Integer Function ChoosePiv(XX, IIleft, IIright) Result (IIpv)
    ! ***********************************
    ! * Choose a Pivot element from XX(Ileft:Iright)
    ! * for Qsort. This routine chooses the median
    ! * of the first, last and mid element of the 
    ! * list.
    ! ***********************************
 
      Real (kind=8), Intent (in) :: XX(:)
      Integer, Intent (in) :: IIleft, IIright
 
      Real (kind=8) :: XXcp(3)
      Integer :: IIpt(3), IImd
 
      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IImd)
      XXcp(3) = XX(IIright)
      IIpt = (/1,2,3/)
 
      CALL InsrtLC(XXcp, IIpt, 1, 3)
 
      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select
 
      Return
    End Function ChoosePiv
 
    ! ***********************************
    Subroutine InsrtLC(XX, IIpt, IIl, IIr)
    ! ***********************************
    ! * Perform an insertion sort of the list 
    ! * XX(:) between index values IIl and IIr.
    ! * IIpt(:) returns the permutations
    ! * made to sort.
    ! ***********************************
 
      Real (kind=8), Intent (inout) :: XX(:)
      Integer, Intent (inout) :: IIpt(:)
      Integer, Intent (in) :: IIl, IIr
 
      Real (kind=8) :: RRtmp
      Integer :: II, JJ, IItmp
 
 
      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL Swap_IN(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do
 
      Return
    End Subroutine InsrtLC
 
  End Subroutine Qsort
 
! ***********************************
! *
  Integer Function Partition(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)
! *
! ***********************************
! * This routine arranges the array X
! * between the index values Ileft and Iright
! * positioning elements smallers than
! * X(Ipv) at the left and the others 
! * at the right.
! * Internal routine used by Qsort.
! ***********************************
 
    Real (kind=8), Intent (inout) :: X(:)
    Integer, Intent (in) :: Ileft, Iright, Ipv
    Integer, Intent (inout), Optional :: Ipt(:)
 
    Real (kind=8) :: Rpv
    Integer :: I
 
    Rpv = X(Ipv)
    CALL Swap(X, Ipv, Iright)
    If (Present(Ipt)) CALL Swap_IN(Ipt, Ipv, Iright)
    Ipvfn = Ileft
 
    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             CALL Swap_IN(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If
 
    CALL Swap(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL Swap_IN(Ipt, Ipvfn, Iright)
 
    Return
  End Function Partition
 
! ***********************************
! *
  Subroutine Swap(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************
 
    Real (kind=8), Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J
 
    Real (kind=8) :: Itmp
 
    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp
 
    Return
  End Subroutine Swap
 
! ***********************************
! *
  Subroutine Swap_IN(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************
 
    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J
 
    Integer :: Itmp
 
    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp
 
    Return
  End Subroutine Swap_IN
end module sort
