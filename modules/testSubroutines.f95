#include "constants.h"

MODULE testSubroutines

CONTAINS

!-----------------Dagger absorption property functions------------------
FUNCTION WAVENUM(ePerp, u) RESULT(ret)
    COMPLEX, INTENT(IN) :: ePerp, u
    COMPLEX :: ret
    ret = SQRT((2*(MASS_N/HBAR)/HBAR)*(ePerp - u))
END FUNCTION WAVENUM

FUNCTION GAMMAK(kn, knm1) RESULT(ret)
    COMPLEX, INTENT(IN) :: kn, knm1
    COMPLEX :: ret
    ret = knm1/kn
END FUNCTION GAMMAK

FUNCTION M(kn, knm1, z) RESULT(ret)
    COMPLEX, INTENT(IN) :: kn, knm1, z
    COMPLEX, DIMENSION(2,2) :: ret
    ret(1,1) = (1.0/2.0)*(1.0 + GAMMAK(kn,knm1))*EXP((0.0,1.0)*(knm1-kn)*z);
    ret(1,2) = (1.0/2.0)*(1.0 - GAMMAK(kn,knm1))*EXP(-(0.0,1.0)*(knm1+kn)*z);
    ret(2,1) = (1.0/2.0)*(1.0 - GAMMAK(kn,knm1))*EXP((0.0,1.0)*(knm1+kn)*z);
    ret(2,2) = (1.0/2.0)*(1.0 + GAMMAK(kn,knm1))*EXP(-(0.0,1.0)*(knm1-kn)*z);
END FUNCTION M

!-----------------Subroutines called many times-------------------------
SUBROUTINE zOffDipCalc(t, z, holdIn)
	REAL(KIND=PREC), INTENT(IN) :: t
	REAL(KIND=PREC), INTENT(OUT) :: z
	REAL(KIND=PREC), INTENT(IN), OPTIONAL :: holdIn
	
	INTEGER :: nDips, i
	REAL(KIND=PREC) :: speed, holdT
	REAL(KIND=PREC), DIMENSION(4) :: dipHeights, dipEnds
			
	IF(PRESENT(holdIn)) THEN
		holdT=holdIn
	ELSE
		holdT=900
	END IF

	!FUN PROJECT: Make this variable so we don't have to recompile to switch # of dips
	nDips = 4 
	
	dipHeights = (/ 0.49_8, 0.380_8, 0.250_8, 0.010_8 /) !3 dip vectors -- from Octet_3dip.txt
	dipEnds = (/ 0.0_8, 40.0_8, 60.0_8, 210.0_8 /) 
	
!    dipHeights = (/0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010/) !9 dip vectors
!    dipEnds = (/0.0,  40.0,  80.0,  100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 300.0/) 

!    dipHeights = (/0.49, 0.250, 0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010/) 
!    dipEnds =     (/0.0_8,  200.0_8,  200.0+holdT, 200.0+holdT+20.0, 200.0+holdT+40.0, 200.0+holdT+50.0, &
!                    200.0+holdT+60.0, 200.0+holdT+70.0, 200.0+holdT+80.0, 200.0+holdT+90.0, &
!                    200.0+holdT+100.0, 200.0+holdT+120.0/) !9 dip PSE vectors
    	
	IF (t .GT. dipEnds(nDips)) THEN
		z = 0.01
		RETURN
	END IF
	
	DO i=1,nDips,1
		IF (dipEnds(i) .GT. t) THEN
			EXIT
		END IF
	END DO
	
	speed = SIGN(1.0_8, dipHeights(i-1) - dipHeights(i))*0.49_8/13.0_8
		
	z = dipHeights(i-1) - speed*(t-dipEnds(i-1))
	
    IF ((speed .GT. 0 .AND. z .LT. dipHeights(i)) .OR. (speed .LT. 0 .AND. z .GT. dipHeights(i))) THEN
        z = dipHeights(i)
	END IF
END SUBROUTINE zOffDipCalc

SUBROUTINE reflect(state, norm, tang)
	USE trackGeometry
	USE constants
	
	REAL(KIND=PREC), DIMENSION(6), INTENT(INOUT) :: state
	REAL(KIND=PREC), DIMENSION(3), INTENT(IN) :: norm, tang
	REAL(KIND=PREC) :: u1, u2, theta, phi, pN, pT, pTprime, pLen, pTarget
	REAL(KIND=PREC), DIMENSION(3) :: tangPrime, newPdir
	
	pTarget = SQRT(state(4)**2 + state(5)**2 + state(6)**2)
	
	CALL cross(norm, tang, tangPrime)
	
	CALL RANDOM_NUMBER(u1)
	CALL RANDOM_NUMBER(u2)
	theta = ASIN(SQRT(u1))
	phi = 2.0_8 * PI * u2

	pN = COS(theta)
	pT = SIN(theta)*COS(phi)
	pTprime = SIN(theta)*SIN(phi)
	
	newPdir = pN*norm + pT*tang + pTprime*tangPrime
	state(4) = newPdir(1)
	state(5) = newPdir(2)
	state(6) = newPdir(3)
	
	pLen = SQRT(state(4)**2 + state(5)**2 + state(6)**2)
	state(4) = state(4) * pTarget/pLen
	state(5) = state(5) * pTarget/pLen
	state(6) = state(6) * pTarget/pLen
	
END SUBROUTINE reflect

SUBROUTINE absorb(ePerp, prob)
    USE constants
    REAL(KIND=PREC), INTENT(IN) :: ePerp
    REAL(KIND=PREC), INTENT(OUT) :: prob
    REAL(KIND=8) :: voxide, woxide, vboron, wboron, vznd, wzns
    COMPLEX, DIMENSION(4) :: pots, zs
    COMPLEX, DIMENSION(2,2) :: mbar
    COMPLEX :: ePerp_c
    INTEGER :: i
    
    ePerp_c = CMPLX(ePerp, 0.0)
    
    voxide = (2*PI*((HBAR/MASS_N)*HBAR))*ABORON*NBORONB2O3
    voxide = voxide + (2*PI*((HBAR/MASS_N)*HBAR))*AOXYGEN*NOXYGENB2O3
    woxide = (HBAR/2)*NBORONB2O3*2200*SIGMABORON + (HBAR/2)*NOXYGENB2O3*SIGMAOXYGEN
    vboron = (2*PI*((HBAR/MASS_N)*HBAR))*ABORON*NBORON
    wboron = (HBAR/2)*NBORON*2200*SIGMABORON
    vzns = (2*PI*((HBAR/MASS_N)*HBAR))*AZINC*NZINC
    vzns = vzns + (2*PI*((HBAR/MASS_N)*HBAR))*ASULFUR*NSULFUR
    wzns = (HBAR/2)*NZINC*SIGMAZINC + (HBAR/2)*NSULFUR*SIGMASULFUR
    
    pots(1) = (0,0)
    pots(2) = CMPLX(voxide, -woxide)
    pots(3) = CMPLX(vboron, -wboron)
    pots(4) = CMPLX(vzns, -wzns)
    zs(1) = (0.0, 0.0)
    zs(2) = (0.0, 0.0)
    zs(3) = (20e-9)
    zs(4) = (10000e-9)
    mbar(1,1) = (1, 0)
    mbar(1,2) = (0, 0)
    mbar(2,1) = (0, 0)
    mbar(2,2) = (1, 0)
    
    DO i = 4,2,-1
        mbar = MATMUL(mbar, M(WAVENUM(ePerp_c, pots(i)), WAVENUM(ePerp_c, pots(i-1)), zs(i-1)))
    END DO
    
    prob = 1.0_8 - REALPART(CONJG(-mbar(2,1)/mbar(2,2))*(-mbar(2,1)/mbar(2,2)))
END SUBROUTINE absorb

SUBROUTINE check_upscatter(state, blockHit, prob)

	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	REAL(KIND=PREC), INTENT(INOUT), DIMENSION(6) :: state
	LOGICAL, INTENT(INOUT) :: blockHit
	REAL(KIND=PREC), INTENT(IN), OPTIONAL :: prob	
	
	REAL(KIND=PREC) :: hitU, upscatterProb
	REAL(KIND=PREC), DIMENSION(6) :: shiftState
	REAL(KIND=PREC), DIMENSION(3,3) :: rotation, invRotation
	REAL(KIND=PREC), DIMENSION(3) :: blockSize, blockPos

	IF(PRESENT(prob)) THEN
		upscatterProb = prob
	ELSE
		upscatterProb = 1.0_8
	END IF
	
	rotation = RESHAPE((/ 0.305229_8, 0.944162_8, 0.124068_8, -0.939275_8,&
				 0.319953_8, -0.124068_8, -0.156836_8, -0.0786645_8, 0.984487_8 /),&
				 SHAPE(rotation))
	invRotation = RESHAPE((/ 0.305229_8, -0.939275_8, -0.156836_8, 0.944162_8,&
				 0.319953_8, -0.0786645_8, 0.124068_8, -0.124068_8, 0.984487_8 /),&
				 SHAPE(invRotation))
	blockSize = (/ 0.0125_8, 0.0125_8, 0.00625_8 /)
	blockPos = (/ 0.2207_8, 0.127_8, -1.4431_8/)
	
	!Shift the neutron position to "block coordinates"
	shiftState = (/ MATMUL(rotation, state(1:3) - blockPos), state(4:6) /)
	
	!Check if the neutron is inside the block
	IF (ABS(shiftState(1)) .GT. blockSize(1) .OR. &
			ABS(shiftState(2)) .GT. blockSize(2) .OR. &
			ABS(shiftState(3)) .GT. blockSize(3)) THEN
		state = state
	ELSE
		CALL RANDOM_NUMBER(hitU)
		IF (hitU < upscatterProb) THEN 
			blockHit = .TRUE.
		!If we don't upscatter, reflect off the block.
		ELSE IF (shiftState(1) .GT. 0 .AND. shiftState(4) .LT. 0) THEN 
			CALL reflect(shiftState, (/ 1.0_8, 0.0_8, 0.0_8 /), (/ 0.0_8, 0.0_8, 1.0_8 /))
			state = (/ MATMUL(invRotation, shiftState(1:3)) + blockPos, shiftState(4:6) /)
		ELSE IF (shiftState(1) .LT. 0 .AND. shiftState(4) .GT. 0) THEN
			CALL reflect(shiftState, (/ -1.0_8, 0.0_8, 0.0_8 /), (/ 0.0_8, 0.0_8, 1.0_8 /))
			state = (/ MATMUL(invRotation, shiftState(1:3)) + blockPos, shiftState(4:6) /)
		ELSE IF (shiftState(2) .GT. 0 .AND. shiftState(5) .LT. 0) THEN 
			CALL reflect(shiftState, (/ 0.0_8, 1.0_8, 0.0_8 /), (/ 0.0_8, 0.0_8, 1.0_8 /))
			state = (/ MATMUL(invRotation, shiftState(1:3)) + blockPos, shiftState(4:6) /)
		ELSE IF (shiftState(2) .LT. 0 .AND. shiftState(5) .GT. 0) THEN
			CALL reflect(shiftState, (/ 0.0_8, -1.0_8, 0.0_8 /), (/ 0.0_8, 0.0_8, 1.0_8 /))
			state = (/ MATMUL(invRotation, shiftState(1:3)) + blockPos, shiftState(4:6) /)
		ELSE IF (shiftState(3) .GT. 0 .AND. shiftState(6) .LT. 0) THEN 
			CALL reflect(shiftState, (/ 0.0_8, 0.0_8, 1.0_8 /), (/ 1.0_8, 0.0_8, 0.0_8 /))
			state = (/ MATMUL(invRotation, shiftState(1:3)) + blockPos, shiftState(4:6) /)
		ELSE 
			PRINT *, "UHOH"
		END IF
	END IF
END SUBROUTINE check_upscatter

!-----------------Subroutines that track neutrons-----------------------
SUBROUTINE trackDaggerHitTime(state, holdT)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	REAL(KIND=PREC), DIMENSION(6), INTENT(INOUT) :: state
	REAL(KIND=PREC), INTENT(IN), OPTIONAL :: holdT
	
	REAL(KIND=PREC) :: t, fracTravel, predX, predZ, energy, zOff, zeta
	REAL(KIND=PREC) :: cleaningTime, settlingTime
	REAL(KIND=PREC), DIMENSION(6) :: prevState
	REAL(KIND=4), DIMENSION(50) :: hitT, hitE
		
	INTEGER :: i, numSteps, nHit
	
	nHit = 0
	hitT = 0.0_8
	hitE = 0.0_8
	t = 0.0_8
	
	!Hard coded in cleaning time, defaults to a 20s hold.
	cleaningTime = 50.0_8
	IF (PRESENT(holdT)) THEN
		settlingTime = cleaningTime + holdT
	ELSE
		settlingTime = cleaningTime + 20.0_8
	END IF
	
	numSteps = settlingTime/dt
	DO i=1,numSteps,1
		CALL symplecticStep(state, dt, energy)
		t = t + dt
	END DO
		
	DO
		prevState = state
		CALL symplecticStep(state, dt, energy)
		t = t + dt
		IF (SIGN(1.0_8, state(2)) .NE. SIGN(1.0_8, prevState(2))) THEN
			fracTravel = ABS(prevState(2))/(ABS(state(2)) + ABS(prevState(2)))
			predX = prevState(1) + fracTravel * (state(1) - prevState(1))
			predZ = prevState(3) + fracTravel * (state(3) - prevState(3))		
			
			CALL zOffDipCalc(t - settlingTime, zOff)
			IF (predX > 0.0_8) THEN
				zeta = 0.5_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 1.0_8)**2)
			ELSE
				zeta = 1.0_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 0.5_8)**2)
			END IF
			!TD offset from central axis: 6" ~0.1524m
		IF (predX > -0.3524_8 .AND. predX < 0.0476_8 .AND. zeta > 0.0_8 .AND. predZ < (-1.5_8 + zOff + 0.2_8)) THEN
			nHit = nHit + 1	
			hitT(nHit) = t - settlingTime
			hitE(nHit) = state(5)*state(5)/(2.0_8*MASS_N)
			IF (nHit .EQ. 50) THEN
				EXIT
			END IF
			IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
				CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
				state = prevState
			ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
				CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
				state = prevState
			ELSE
				PRINT *, "UHOH"
			END IF
!			WRITE(1) t - (20.0_8 + 50.0_8), predX, predZ - zOff
!			EXIT
			ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8) .AND. &
				predZ < (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
				ABS(predX + 0.1524_8) < (0.40_8 + 2.0179_8*(predZ + 1.5_8 - zOff - 0.2_8))/2.0_8) THEN
				IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
					CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
					CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				END IF
!				PRINT *, "BOUNCE LOWER"
!				PRINT *, predX, predZ, zOff
			ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
				predZ < (-1.5_8 + zOff + 0.2_8 + 0.2667_8) .AND. &
				ABS(predX + 0.1524_8) < 0.69215_8/2.0_8) THEN
				IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
					CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
					CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				END IF
!				PRINT *, "BOUNCE UPPER"
!				PRINT *, predX, predZ, zOff
			END IF			
			
			IF (t > 2000) THEN
				EXIT
			END IF
		END IF
	END DO
	WRITE(1) energy, hitT, hitE
END SUBROUTINE trackDaggerHitTime

SUBROUTINE fixedEffDaggerHitTime(state, holdT)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	REAL(KIND=PREC), DIMENSION(6), INTENT(INOUT) :: state
	REAL(KIND=PREC), INTENT(IN), OPTIONAL :: holdT
	
	INTEGER :: i, numSteps, nHit, nHitHouseLow, nHitHouseHigh
	REAL(KIND=PREC) :: t, fracTravel, predX, predZ, energy, zOff, zeta
	REAL(KIND=PREC) :: cleaningTime, settlingTime, absProb, absU, deathTime
	REAL(KIND=PREC), DIMENSION(6) :: prevState
		
	t = 0.0_8
	nHit = 0
	nHitHouseLow = 0
	nHitHouseHigh = 0
	
	!Hardcoded in cleaning time, holdT defaults to 20s
	cleaningTime = 50.0_8
	IF (PRESENT(holdT)) THEN
		settlingTime = cleaningTime + holdT
	ELSE
		settlingTime = cleaningTime + 20.0_8
	END IF
	
	!FixedEff will also have neutrons decay	
	CALL RANDOM_NUMBER(deathTime)
	deathTime = -877.7*LOG(deathTime)

	numSteps = settlingTime/dt
	DO i=1,numSteps,1
		CALL symplecticStep(state, dt, energy)
		t = t + dt
	END DO
	
	DO
		prevState = state
		CALL symplecticStep(state, dt, energy)
		t = t + dt
		IF (SIGN(1.0_8, state(2)) .NE. SIGN(1.0_8, prevState(2))) THEN
			fracTravel = ABS(prevState(2))/(ABS(state(2)) + ABS(prevState(2)))
			predX = prevState(1) + fracTravel * (state(1) - prevState(1))
			predZ = prevState(3) + fracTravel * (state(3) - prevState(3))
			
			CALL zOffDipCalc(t - settlingTime, zOff)
			IF (predX > 0.0_8) THEN
				zeta = 0.5_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 1.0_8)**2)
			ELSE
				zeta = 1.0_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 0.5_8)**2)
			END IF
			!TD offset from central axis: 6" ~0.1524m
			IF (predX > -0.3524_8 .AND. predX < 0.0476_8 .AND. zeta > 0.0_8 .AND. predZ < (-1.5_8 + zOff + 0.2_8)) THEN
				nHit = nHit + 1
				CALL absorb(state(5)*state(5)/(2.0_8*MASS_N), absProb)
				CALL RANDOM_NUMBER(absU)
				IF (absU .LT. absProb) THEN
					WRITE(1) t - settlingTime, energy, state(5)*state(5)/(2.0_8*MASS_N), &
					predX, 0.0_8, predZ, zOff, nHit, nHitHouseLow, nHitHouseHigh
					EXIT
				END IF
                
				IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
					CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
					CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE
					PRINT *, "UHOH"
				END IF
            ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8) .AND. &
					predZ < (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
					ABS(predX + 0.1524_8) < (0.40_8 + 2.0179_8*(predZ + 1.5_8 - zOff - 0.2_8))/2.0_8) THEN
				IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
					CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
				CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				END IF
				nHitHouseLow = nHitHouseLow + 1
			ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
					predZ < (-1.5_8 + zOff + 0.2_8 + 0.2667_8) .AND. &
					ABS(predX + 0.1524_8) < 0.69215_8/2.0_8) THEN
				IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
					CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
					CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				END IF
				nHitHouseHigh = nHitHouseHigh + 1
			END IF
			IF (t-settlingTime > deathTime) THEN
				EXIT
			END IF
		END IF
	END DO
END SUBROUTINE fixedEffDaggerHitTime

SUBROUTINE trackDaggerAndBlock(state, holdT)

	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	REAL(KIND=PREC), DIMENSION(6), INTENT(INOUT) :: state
	REAL(KIND=PREC), INTENT(IN), OPTIONAL :: holdT
	
	INTEGER :: i, numSteps, nHit, nHitHouseLow, nHitHouseHigh
	LOGICAL :: blockHit, dagHit
	REAL(KIND=PREC) :: t, fracTravel, predX, predZ, energy, zOff, zeta, blockProb
	REAL(KIND=PREC) :: cleaningTime, settlingTime, absProb, absU, deathTime
	REAL(KIND=PREC), DIMENSION(6) :: prevState
	
	nHit = 0
	nHitHouseLow = 0
	nHitHouseHigh = 0
	blockHit = .FALSE.	
	dagHit = .FALSE. 
	t = 0.0_8
	
	!Hard coded in cleaning time, defaults to 20s hold
	cleaningTime = 50.0_8 
	IF (PRESENT(holdT)) THEN
		settlingTime = cleaningTime + holdT
	ELSE
		settlingTime = cleaningTime + 20.0_8
	END IF 
	
	!Block upscatter probability (100% upscatter hardcoded in for now)
	blockProb = 1.0_8
	
	!Functionality for neutron decay. Presently commented out, set to kill after 2000s
	!CALL RANDOM_NUMBER(deathTime)
	!deathTime = -877.7*LOG(deathTime)
	deathTime = 2000.0

	numSteps = settlingTime/dt
	DO i=1,numSteps,1
		CALL symplecticStep(state, dt, energy)
		t = t + dt
		CALL check_upscatter(state, blockHit, blockProb)
		IF (blockHit) THEN
		!	WRITE(2) t, energy, nHit, nHitHouseLow, nHitHouseHigh, &
		!	state(1), state(2), state(3)
			EXIT
		END IF
	END DO
	
	DO
		prevState = state
		CALL symplecticStep(state, dt, energy)
		t = t + dt
		
		IF ((.NOT. blockHit) .AND. (.NOT. dagHit)) THEN
			CALL check_upscatter(state, blockHit, blockProb) 
			IF (blockHit) THEN
			!	WRITE(2) t, energy, nHit, nHitHouseLow, nHitHouseHigh, &
			!	state(1), state(2), state(3)
			END IF 
		END IF
		IF (blockHit) THEN 
			PRINT *, "hit a block!"
			EXIT
		END IF
		
		IF (SIGN(1.0_8, state(2)) .NE. SIGN(1.0_8, prevState(2))) THEN
			fracTravel = ABS(prevState(2))/(ABS(state(2)) + ABS(prevState(2)))
			predX = prevState(1) + fracTravel * (state(1) - prevState(1))
			predZ = prevState(3) + fracTravel * (state(3) - prevState(3))
			
			CALL zOffDipCalc(t - settlingTime, zOff)
			IF (predX > 0.0_8) THEN
				zeta = 0.5_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 1.0_8)**2)
			ELSE
				zeta = 1.0_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 0.5_8)**2)
			END IF
			!TD offset from central axis: 6" ~0.1524m
			IF (predX > -0.3524_8 .AND. predX < 0.0476_8 .AND. &
					zeta > 0.0_8 .AND. predZ < (-1.5_8 + zOff + 0.2_8)) THEN
				nHit = nHit + 1
				CALL absorb(state(5)*state(5)/(2.0_8*MASS_N), absProb)
				CALL RANDOM_NUMBER(absU)
				IF (absU < absProb) THEN
					WRITE(1) t - settlingTime, energy, state(5)*state(5)/(2.0_8*MASS_N), &
					predX, 0.0_8, predZ, zOff, nHit, nHitHouseLow, nHitHouseHigh
					dagHit = .TRUE.
					!PRINT *, "hit dagger"
					EXIT
				END IF
                
				IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
					CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
					CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE
					PRINT *, "UHOH"
				END IF
            ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8) .AND. &
					predZ < (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
					ABS(predX + 0.1524_8) < (0.40_8 + 2.0179_8*(predZ + 1.5_8 - zOff - 0.2_8))/2.0_8) THEN
				IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
					CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
				CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				END IF
				nHitHouseLow = nHitHouseLow + 1
			ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
					predZ < (-1.5_8 + zOff + 0.2_8 + 0.2667_8) .AND. &
					ABS(predX + 0.1524_8) < 0.69215_8/2.0_8) THEN
				IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
					CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
					CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), (/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				END IF
				nHitHouseHigh = nHitHouseHigh + 1
			END IF
			IF (t-settlingTime > deathTime) THEN
				EXIT
			END IF
		END IF
	END DO
	
END SUBROUTINE trackDaggerAndBlock

!-----------------Subroutines for debug purposes------------------------
SUBROUTINE compPots()
	USE forcesAndPotential
	IMPLICIT NONE
	REAL(KIND=PREC) :: x, y, z, fx, fy, fz, totalU
	INTEGER :: i
	
	z = -1.49_8
	x = 0.1_8
	y = 0.1_8
	DO i=0,20,1
		CALL totalForce(x, y, z + i*(0.01_8/20.0_8), fx, fy, fz, totalU)
		PRINT *, totalU, fx, fy, fz
	END DO
END SUBROUTINE compPots

SUBROUTINE trackEnergyGain(state, energy_start, energy_end, sympT, freq)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	REAL(KIND=PREC), DIMENSION(6), INTENT(INOUT) :: state
	REAL(KIND=PREC), INTENT(OUT) :: energy_start, energy_end
	REAL(KIND=PREC), OPTIONAL, INTENT(INOUT) :: sympT
	REAL(KIND=PREC), OPTIONAL, INTENT(IN) :: freq
    REAL(KIND=PREC) :: eta, prevEta, t, energy, t_end
	INTEGER :: i, numSteps, triggered, rising
    
    t = 0.0_8
    
    triggered = 0
    rising = 0
    prevEta = 10.0_8
    eta = 10.0_8

	numSteps = 1000e0/dt
!	totalU = 0.0_8
	
!	CALL calcEnergy(state, energy_start)
    energy_start = 0.0_8
!	energy = 0.0
	
	DO i=1,numSteps,1
		IF(PRESENT(sympT)) THEN
			CALL symplecticStep(state, dt, energy, sympT, freq)
		ELSE
			CALL symplecticStep(state, dt, energy)
			t = t+dt
		END IF
!        IF (triggered .EQ. 0) THEN
        prevEta = eta
        IF (state(1) > 0) THEN
            eta = 1.0_8 - SQRT(state(1)**2 + (SQRT(state(2)**2 + state(3)**2) - 0.5_8)**2)
        END IF
        IF (state(1) <= 0) THEN
            eta = 0.5_8 - SQRT(state(1)**2 + (SQRT(state(2)**2 + state(3)**2) - 1.0_8)**2)
        END IF
        IF ((rising .EQ. 0) .AND. (prevEta < eta)) THEN
            rising = 1
        END IF
        IF ((rising .EQ. 1) .AND. (prevEta > eta) .AND. (triggered .EQ. 0)) THEN
!            PRINT *, t, energy
            triggered = 1
            rising = 0
            energy_start = energy
        END IF
        IF ((rising .EQ. 1) .AND. (prevEta > eta) .AND. (triggered .EQ. 1)) THEN
!            PRINT *, t, energy
            t_end = t
            rising = 0
            energy_end = energy
        END IF
!        END IF
		IF (100.0_8*energy/(MASS_N*GRAV) > 38.0_8 + 5.0_8) THEN
!			PRINT *, "DEAD"
			EXIT
		END IF
	END DO
!	PRINT *, energy_start, energy_end
!	PRINT *, t_end, energy_end
END SUBROUTINE trackEnergyGain

SUBROUTINE testEnergyGain(freq, height, sympT, eStart, eEnd)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	REAL(KIND=PREC), INTENT(IN) :: freq, height
	REAL(KIND=PREC), INTENT(INOUT) :: sympT
	REAL(KIND=PREC), INTENT(OUT) :: eStart, eEnd

	REAL(KIND=PREC), DIMENSION(6) :: state
	INTEGER :: i, numSteps
	
	state = (/0.05_8, 0.0_8, height, 0.0_8, 0.0_8, 0.0_8/)
	CALL calcEnergy(state, eStart)

	numSteps = 10e0/dt
	
	DO i=1,numSteps,1
		CALL symplecticStep(state, dt, eEnd, sympT, freq)
		IF (state(3) >= -1.5_8 + ((height + 1.5) / 4.0_8) .AND. state(6) > 0) THEN
			EXIT
		END IF
	END DO
!	PRINT *, energy_start, energy_end
END SUBROUTINE testEnergyGain

SUBROUTINE trackAndPrint(state, sympT)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	REAL(KIND=PREC), DIMENSION(6), INTENT(INOUT) :: state
	REAL(KIND=PREC), OPTIONAL, INTENT(IN) :: sympT
	
	REAL(KIND=PREC) :: pr, rdot, r, pphi, phidot, phi, ptheta, thetadot, theta,&
		ldot, l, phiOffset, prevPhi, totalU, totalKE, t, energy, &
		fx_dan, fy_dan, fz_dan, e_dan, fx_nate, fy_nate, fz_nate, e_nate
	INTEGER :: i, numSteps, numPoints, modulus
		
!	numSteps = 250e0_8/dt
	numSteps = 1000e0/dt
	t = 0.0_8
	totalU = 0.0_8
	
!	CALL calcEnergy(state, energy)
	energy = 0.0
	
	DO i=1,numSteps,1
!		IF(t > 475.0_8 .AND. t < 480.0_8) THEN
		IF(INT(dt*10_8*i)-INT(dt*10_8*(i-1)) .NE. 0) THEN
!		IF(1 .EQ. 1) THEN
!			energy = totalU + SUM(state(1,4:6)**2)/(2.0_8*MASS_N)

			PRINT *, dt*i, state(1), state(2), state(3),&
			state(4)/MASS_N, state(5)/MASS_N, state(6)/MASS_N, energy,&
			totalU!, fx, fy, fz

			!PRINT *, dt*i, energy
		END IF
		IF(PRESENT(sympT)) THEN
			CALL symplecticStep(state, dt, energy, t, 60.0_8)
		ELSE
			CALL symplecticStep(state, dt, energy)
!			CALL totalForce(state(1), state(2), state(3), fx_nate, fy_nate, fz_nate, e_nate)
!			CALL totalForceDan(state(1), state(2), state(3), fx_dan, fy_dan, fz_dan, e_dan)
!			PRINT *, t, (fx_nate-fx_dan)/fx_nate, (fy_nate-fy_dan)/fy_nate,&
!			(fz_nate-fz_dan)/fz_nate, (e_nate-e_dan)/e_nate
			t = t+dt
		END IF
	END DO
	!PRINT *, t, state(1), state(2), state(3), state(4), state(5), state(6)
END SUBROUTINE trackAndPrint

SUBROUTINE calcx0Mesh()
	USE forcesAndPotential
	USE constants
	!Calculate a field map on one plane for diagnostics
	IMPLICIT NONE
	REAL(KIND=PREC) :: x, y, z, fx, fy, fz, totalUplus
	REAL(KIND=PREC) :: x0, y0, z0
	INTEGER :: xIt, yIt, zIt, aIt
	
	x0 = -2
	y0 = -2
	z0 = -2
	
	DO xIt=0,100,1
		DO yIt=0,100,1
			DO zIt=0,100,1
				x = x0 + 4 * xIt/100.0_8
				y = y0 + 4 * yIt/100.0_8
				z = z0 + 4 * zIt/100.0_8
				CALL totalPotential(x, y, z, totalUplus)
				PRINT *, x, y, z, totalUplus/MASS_N
			END DO
		END DO
	END DO
END SUBROUTINE calcx0Mesh

END MODULE
