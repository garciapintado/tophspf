MODULE moduleHspfPwater


  USE nrtype
  USE topo_module, ONLY : xlowtopo,ylowtopo,xhitopo,yhitopo, &
                          dxtopo, dytopo, mxtopo, mytopo         ! to get geometadata and DTM
  IMPLICIT NONE

  PRIVATE


  PUBLIC :: ConstructSURS
  PUBLIC :: DestructSURS

 ! REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: surs       !actual surface detention storage       [mm] storage

  CONTAINS

    SUBROUTINE ConstructSURS

    INTEGER :: i

    ALLOCATE(surs(mytopo(1),mxtopo(1))) ! as topoHSPF storage (not good)
    surs = -0.02d0
    surs(1:150,1:50) = RESHAPE( (/(i, i=1,7500)/)*1.0d-05,[150,50], ORDER=[2,1]) ! fill in byrow=TRUE

    !surs(1:INT(mytopo(1)/2),INT(mxtopo(1)/2):mxtopo(1)) = 0.02d0    !positive state in the upper-right cuadrant
  END SUBROUTINE ConstructSURS

  SUBROUTINE DestructSURS
    DEALLOCATE(surs)
  END SUBROUTINE DestructSURS


END MODULE modulehspfpwater
