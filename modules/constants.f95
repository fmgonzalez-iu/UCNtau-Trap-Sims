#include "constants.h"

MODULE constants
	IMPLICIT NONE
	REAL(KIND=PREC) :: PI
	REAL(KIND=PREC), DIMENSION(4) :: a
	REAL(KIND=PREC), DIMENSION(4) :: b
	REAL(KIND=PREC) :: dt, liptime
	INTEGER :: nsep
	SAVE
END MODULE
