#include "constants.h"

MODULE forcesAndPotential
	real(kind=PREC) :: minU

CONTAINS

SUBROUTINE totalPotential(x, y, z, totalU, t)
	IMPLICIT NONE
	REAL(KIND=PREC), INTENT(IN) :: x, y, z
	REAL(KIND=PREC), INTENT(OUT) :: totalU
	REAL(KIND=PREC), OPTIONAL, INTENT(IN) :: t
	
	IF(PRESENT(t)) THEN
		CALL potential(x, y, z, totalU, t)
	ELSE
		CALL potential(x, y, z, totalU)
	END IF
	totalU = totalU - minU
END SUBROUTINE totalPotential

SUBROUTINE totalForce(x, y, z, fx, fy, fz, totalU, t)
	IMPLICIT NONE
	REAL(KIND=PREC), INTENT(IN) :: x, y, z
	REAL(KIND=PREC), INTENT(OUT) :: fx, fy, fz
	REAL(KIND=PREC), INTENT(OUT) :: totalU
	REAL(KIND=PREC), OPTIONAL, INTENT(IN) :: t
	
	IF(PRESENT(t)) THEN
		CALL force(x, y, z, fx, fy, fz, totalU, t)
	ELSE
		CALL force(x, y, z, fx, fy, fz, totalU, 0, 0)
	END IF
	totalU = totalU - minU
END SUBROUTINE totalForce

!SUBROUTINE totalForceDan(x, y, z, fx, fy, fz, totalU, t)
!	IMPLICIT NONE
!	real(kind=PREC), intent(in) :: x, y, z
!	real(kind=PREC), intent(out) :: fx, fy, fz
!	real(kind=PREC), intent(out) :: totalU
!	real(kind=PREC), optional, intent(in) :: t
!	
!	IF(PRESENT(t)) THEN
!		CALL force_dan(x, y, z, fx, fy, fz, t)
!	ELSE
!		CALL force_dan(x, y, z, fx, fy, fz)
!	END IF
!
!END SUBROUTINE totalForceDan


SUBROUTINE calcEnergy(state, energy, t)
	IMPLICIT NONE
	REAL(KIND=PREC), DIMENSION(6), INTENT(IN) :: state
	REAL(KIND=PREC), INTENT(OUT) :: energy
	REAL(KIND=PREC), OPTIONAL, INTENT(IN) :: t
	REAL(KIND=PREC) :: totalU
	
	IF(PRESENT(t)) THEN
		CALL totalPotential(state(1), state(2), state(3), totalU, t)
	ELSE
		CALL totalPotential(state(1), state(2), state(3), totalU)
	END IF
	
	energy = totalU + SUM(state(4:6)*state(4:6))/(2.0_8*MASS_N)
END SUBROUTINE calcEnergy

END MODULE
