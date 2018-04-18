#include "constants.h"

MODULE convertTrace

CONTAINS

SUBROUTINE loadTrace(x, y, z, n)
	IMPLICIT NONE
	INTEGER, INTENT(OUT) :: n
	REAL(KIND=PREC), INTENT(OUT), DIMENSION(:), ALLOCATABLE :: x(:), y(:), z(:)
	
	INTEGER :: i, error, nX, nY, nZ
	REAL(KIND=PREC) :: dummy
	
	OPEN(UNIT=10, FILE="./Chris_Acc_Trap_Displacements/xvals.bin", &
		FORM="UNFORMATTED",ACCESS="STREAM", STATUS="OLD", &
		ACTION = "READ")
	!Figure how long the data set is:
	nX = 0
	DO
		READ(10,IOSTAT=error), dummy
		IF (error < 0) EXIT
		nX = nX + 1
	END DO
		
	REWIND(10)
	!Allocate and read data from file
	ALLOCATE(x(nX))
	DO i=1,nX,1
		READ(10,IOSTAT=error) x(i)
		IF (error .NE. 0) EXIT
		
	END DO
	CLOSE(10)
	
	! Same thing in Y
	OPEN(UNIT=11, FILE="./Chris_Acc_Trap_Displacements/yvals.bin", &
		FORM="UNFORMATTED", ACCESS="STREAM", STATUS="OLD", &
		ACTION = "READ")
	nY = 0
	DO
		READ(11,IOSTAT = error), dummy
		IF (error < 0) EXIT
		nY = nY + 1
	END DO
	
	IF (nX .NE. nY) THEN
		PRINT *, "Length of Y file does not match length of X file!"
		STOP
	END IF
	
	REWIND(11)
	ALLOCATE(y(nY))
	DO i=1,nY
		READ(11,IOSTAT=error) y(i)
		IF (error .NE. 0) EXIT
	END DO
	CLOSE(11)
	
	! Same thing in Z
	OPEN(UNIT=12, FILE="./Chris_Acc_Trap_Displacements/zvals.bin", &
		FORM="UNFORMATTED", ACCESS="STREAM", STATUS="OLD", &
		ACTION = "READ")
	nZ = 0
	DO
		READ(12, IOSTAT = error), dummy
		IF (error < 0) EXIT
		nZ = nZ + 1
	END DO
	
	IF (nX .NE. nZ) THEN
		PRINT *, "Length of Z file does not match length of X file!"
		STOP
	END IF
	
	REWIND(12)
	ALLOCATE(z(nZ))
	DO i=1,nZ
		READ(12, IOSTAT=error) z(i)
		IF (error .NE. 0) EXIT
	END DO
	CLOSE(12)
	
	n = nX
    PRINT *, "The number of events in heating file is: "
	PRINT *, n	
END SUBROUTINE loadTrace

SUBROUTINE shiftTrace(t,trX,trY,trZ,n, xi,yi,zi, xf,yf,zf)
	IMPLICIT NONE
	INTEGER, INTENT(IN) ::  n
	REAL(KIND=PREC), INTENT(IN) :: t
	REAL(KIND=PREC), INTENT(IN), ALLOCATABLE :: trX(:), trY(:), trZ(:)
	REAL(KIND=PREC), INTENT(IN) :: xi, yi, zi
	REAL(KIND=PREC), INTENT(OUT) :: xf, yf, zf
	
	INTEGER :: iLow, iHigh
	REAL(KIND=PREC) :: sampDT, frac
	
	sampDT = 0.0004 !Constant sampling rate
		
	iLow  = FLOOR(t/sampDT)
	iHigh = iLow + 1
	frac = (t - iLow*sampDT)/(sampDT) ! Load fractions as function of low/high indices	
	iLow  = MOD(iLow,n) + 1 !convert to the actual amount -- add one for FORTRAN
	iHigh = MOD(iHigh,n)+ 1
	
	xf = xi + trX(iLow) + frac*(trX(iHigh) - trX(iLow)) !shift our variables
	yf = yi + trY(iLow) + frac*(trY(iHigh) - trY(iLow))
	zf = zi + trZ(iLow) + frac*(trZ(iHigh) - trZ(iLow))

END SUBROUTINE shiftTrace

END MODULE
