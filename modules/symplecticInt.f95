#include "constants.h"

MODULE symplecticInt

CONTAINS

SUBROUTINE symplecticStep(state, deltaT, energy, t, trX,trY,trZ,lenTr)
	!Take 1 step of length deltaT using symplectic integration scheme.
	!See Candy & Rozmus 1991 and McLachlan & Atela 1992
	USE constants
	USE forcesAndPotential
	USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT
	IMPLICIT NONE
	REAL(KIND=PREC), INTENT(INOUT), DIMENSION(6) :: state
	REAL(KIND=PREC), INTENT(IN) :: deltaT
	REAL(KIND=PREC), INTENT(OUT) :: energy
	REAL(KIND=PREC), OPTIONAL, INTENT(INOUT) :: t
	REAL(KIND=PREC), INTENT(IN), ALLOCATABLE, &
						DIMENSION(:), OPTIONAL :: trX(:), trY(:), trZ(:)
	INTEGER, OPTIONAL :: lenTr

	REAL(KIND=PREC) :: fx, fy, fz, totalU
	INTEGER :: n
	
	!IF(PRESENT(lenTr)) THEN
	!	ALLOCATE trX(lenTr)
	!	ALLOCATE trY(lenTr)
	!	ALLOCATE trZ(lenTr)
	!END IF
					
	!REAL(KIND=PREC), OPTIONAL, INTENT(IN) :: freq

	
	!Only gets compiled if VARIABLE isn't defined at the top
	!For each step of the integration,
		!Calc force
	!Unrolling the loop - first iteration gives energy but subsequent won't
	n = 1
	IF(PRESENT(t)) THEN
		IF(PRESENT(lenTr)) THEN
			CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU, t, &
							trX, trY, trZ, lenTr)
		ELSE 
			CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU, t)
		END IF
		energy = totalU + SUM(state(4:6)*state(4:6))/(2.0_8*MASS_N)! &
!			+ state(1)*GRAV*MASS_N
		t = t + a(n)*deltaT
	ELSE
		CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU)
		energy = totalU + SUM(state(4:6)*state(4:6))/(2.0_8*MASS_N)! &
!			+ state(1)*GRAV*MASS_N
	END IF
	!CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU)
	!energy = totalU + SUM(state(4:6)*state(4:6))/(2.0_8*MASS_N) &
	!	+ state(1)*GRAV*MASS_N - uZero
	!PRINT *, "Force: ", fx, fy, fx
	!Use integration constants (a's and b's)
	!Reference:
	state(4) = state(4) + b(n)*fx*deltaT
	state(5) = state(5) + b(n)*fy*deltaT
	state(6) = state(6) + b(n)*fz*deltaT
	state(1) = state(1) + a(n)*state(4)*deltaT/MASS_N
	state(2) = state(2) + a(n)*state(5)*deltaT/MASS_N
	state(3) = state(3) + a(n)*state(6)*deltaT/MASS_N

	DO n=2,4,1
		!Calc force
		IF(PRESENT(t)) THEN
			IF(PRESENT(lenTr)) THEN
				CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU, t, &
								trX, trY, trZ, lenTr)
			ELSE 
				CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU, t)
			END IF
			t = t + a(n)*deltaT
		ELSE
			CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU)
		END IF
		!CALL totalForce(state(1), state(2), state(3), fx, fy, fz, totalU)
		!PRINT *, "Force: ", fx, fy, fx
		!Use integration constants (a's and b's)
		!Reference:
		!px=px+b(n)*fx*deltaT
		!py=py+b(n)*fy*deltaT
		!pz=pz+b(n)*fz*deltaT
		!x=x+a(n)*px*deltaT/MASS_N
		!y=y+a(n)*py*deltaT/MASS_N
		!z=z+a(n)*pz*deltaT/MASS_N
		state(4) = state(4) + b(n)*fx*deltaT
		state(5) = state(5) + b(n)*fy*deltaT
		state(6) = state(6) + b(n)*fz*deltaT
		state(1) = state(1) + a(n)*state(4)*deltaT/MASS_N
		state(2) = state(2) + a(n)*state(5)*deltaT/MASS_N
		state(3) = state(3) + a(n)*state(6)*deltaT/MASS_N
	END DO
END SUBROUTINE symplecticStep

END MODULE
