MODULE selfcons
  USE nrtype
  USE modsim
  USE vars
  USE intldc
  USE intpole
  IMPLICIT NONE
    
CONTAINS

  FUNCTION scenergy(k)
    USE intpole
    IMPLICIT NONE
    REAL(QP), INTENT(IN) :: k
    REAL(QP) :: scenergy

    scenergy=k*1.0_qp

  END FUNCTION scenergy 
  FUNCTION continuum() 
    IMPLICIT NONE
    
    ! Read info file
    call rdInfo()

    

  END FUNCTION continuum
END MODULE selfcons