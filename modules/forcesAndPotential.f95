#include "constants.h"

MODULE forcesAndPotential
	
	USE convertTrace
	REAL(KIND=PREC) :: minU

CONTAINS

SUBROUTINE totalPotential(x, y, z, totalU, t, trX, trY, trZ, n)
	IMPLICIT NONE
	REAL(KIND=PREC), INTENT(IN) :: x, y, z
	REAL(KIND=PREC), INTENT(OUT) :: totalU
	INTEGER, OPTIONAL, INTENT(IN) :: n
	REAL(KIND=PREC), OPTIONAL, INTENT(IN) :: t
	REAL(KIND=PREC), OPTIONAL, INTENT(IN), ALLOCATABLE :: trX(:), trY(:), trZ(:)
	REAL(KIND=PREC) ::  xf, yf, zf
	
	IF(PRESENT(t)) THEN
		IF(PRESENT(n)) THEN
		!	ALLOCATE(trX(n))
		!	ALLOCATE(trY(n))
		!	ALLOCATE(trZ(n))
			CALL shiftTrace(t,trX,trY,trZ, n, x,y,z, xf,yf,zf)
		ELSE ! If we can't properly shift the traces, just don't shift!
			xf = x
			yf = y
			zf = z
		END IF
		CALL potential(xf, yf, zf, totalU, t)
	ELSE
		CALL potential(x, y, z, totalU, 0.0_8)
	END IF
	totalU = totalU - minU
END SUBROUTINE totalPotential

SUBROUTINE totalForce(x, y, z, fx, fy, fz, totalU, t, trX, trY, trZ, n)
	IMPLICIT NONE
	REAL(KIND=PREC), INTENT(IN) :: x, y, z
	REAL(KIND=PREC), INTENT(OUT) :: fx, fy, fz
	REAL(KIND=PREC), INTENT(OUT) :: totalU
	INTEGER, OPTIONAL, INTENT(IN) :: n
	REAL(KIND=PREC), OPTIONAL, INTENT(IN) :: t
	REAL(KIND=PREC), OPTIONAL, INTENT(IN), ALLOCATABLE :: trX(:), trY(:), trZ(:)
	REAL(KIND=PREC) ::  xf, yf, zf
	
	IF(PRESENT(t)) THEN
		IF(PRESENT(n)) THEN
		!	ALLOCATE(trX(n))
		!	ALLOCATE(trY(n))
		!	ALLOCATE(trZ(n))
			CALL shiftTrace(t,trX,trY,trZ, n, x,y,z, xf,yf,zf)
		ELSE ! If we can't properly shift the traces, just don't shift!
			xf = x
			yf = y
			zf = z
		END IF
		
		CALL force(xf, yf, zf, fx, fy, fz, totalU, t)
	ELSE
		CALL force(x, y, z, fx, fy, fz, totalU, 0.0_8)
	END IF
	
	totalU = totalU - minU
END SUBROUTINE totalForce

SUBROUTINE calcEnergy(state, energy, t, trX, trY,trZ, n)
	IMPLICIT NONE
	REAL(KIND=PREC), DIMENSION(6), INTENT(IN) :: state
	REAL(KIND=PREC), INTENT(OUT) :: energy
	INTEGER, OPTIONAL, INTENT(IN) :: n
	REAL(KIND=PREC), OPTIONAL, INTENT(IN) :: t
	REAL(KIND=PREC), OPTIONAL, INTENT(IN), ALLOCATABLE :: trX(:), trY(:), trZ(:)
	REAL(KIND=PREC) ::  xf, yf, zf, totalU
		
	IF(PRESENT(t)) THEN
		IF(PRESENT(n)) THEN
		!	ALLOCATE(trX(n))
		!	ALLOCATE(trY(n))
		!	ALLOCATE(trZ(n))
			CALL shiftTrace(t,trX,trY,trZ, n, state(1),state(2),state(3), xf,yf,zf)
		ELSE
			xf = state(1)
			yf = state(2)
			zf = state(3)
		END IF
		CALL totalPotential(xf, yf, zf, totalU, t)
	ELSE
		CALL totalPotential(state(1), state(2), state(3), totalU, 0.0_8)
	END IF
	
	energy = totalU + SUM(state(4:6)*state(4:6))/(2.0_8*MASS_N)
END SUBROUTINE calcEnergy

END MODULE
