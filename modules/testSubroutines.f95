#include "constants.h"

MODULE testSubroutines

CONTAINS

!-----------------------------------------------------------------------
!-----------------Dagger absorption property functions------------------
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
!-----------------Subroutines called many times-------------------------
!-----------------------------------------------------------------------
! Check if we hit the cleaner and clean our neutrons
SUBROUTINE cleaning(state, prevState, cleanZ, cleanHit)
	
	REAL(KIND=PREC), DIMENSION(6), INTENT(IN) :: state
	REAL(KIND=PREC), DIMENSION(6), INTENT(IN) :: prevState
	REAL(KIND=PREC), INTENT(IN) :: cleanZ
	LOGICAL, INTENT(OUT) :: cleanHit
			
	! Initialize variables
	cleanHit = .FALSE.
			
	! Stop trajectory if we pass through the cleaner, which is only on +x
	! Cleaner has 100% efficiency. Could probably work in a reflection prob too?
	IF (state(2) > 0) THEN
		IF((prevState(3) < -1.5 + cleanZ) .AND. (state(3) > -1.5 + cleanZ)) THEN
			cleanHit = .TRUE.
		ELSE IF((prevState(3) > -1.5 + cleanZ) .AND. (state(3) < -1.5 + cleanZ)) THEN
			cleanHit = .TRUE.
		END IF
	END IF
		
END SUBROUTINE cleaning
! Check if we hit the dagger
SUBROUTINE checkDagHit(state, prevState, zOff, dagHit)

	REAL(KIND=PREC), DIMENSION(6), INTENT(INOUT) :: state, prevState
	REAL(KIND=PREC), INTENT(IN) :: zOff
	INTEGER, INTENT(OUT) :: dagHit
	
	REAL(KIND=PREC) :: fracTravel, predX, predZ, zeta, absU, absProb
	
	! Initialize variables:
	! dagHit is a switch variable for position the neutron hits. Truth table:
	!  -1: ERROR!!!
	! 	0: No hits
	! 	1: Hits the dagger, absorbed
	!	2: Hits the dagger, reflected
	! 	3: Hits the low housing
	!	4: Hits the high housing
	! I've copied this truth table into the main program for posterity.
	dagHit = 0
	
	! Check if we've crossed the plane. This initialized dagger hit
	IF (SIGN(1.0_8, state(2)) .NE. SIGN(1.0_8, prevState(2))) THEN
		fracTravel = ABS(prevState(2))/(ABS(state(2)) + ABS(prevState(2)))
		predX = prevState(1) + fracTravel * (state(1) - prevState(1))
		predZ = prevState(3) + fracTravel * (state(3) - prevState(3))
		
		! calculate zeta
		IF (predX > 0.0_8) THEN
			zeta = 0.5_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 1.0_8)**2)
		ELSE
			zeta = 1.0_8 - SQRT(predX**2 + (ABS(predZ - zOff) - 0.5_8)**2)
		END IF
		
		! Check if we hit the dagger first.
		!TD offset from central axis: 6" ~0.1524m
		IF (predX > -0.3524_8 .AND. predX < 0.0476_8 .AND. &
				zeta > 0.0_8 .AND. predZ < (-1.5_8 + zOff + 0.2_8)) THEN
					
			! Need to check if we are absorbed or reflected
			CALL absorb(state(5)*state(5)/(2.0_8*MASS_N), absProb)
			CALL RANDOM_NUMBER(absU)
			IF (absU < absProb) THEN
				! Hit the dagger and absorbed!
				dagHit = 1
			ELSE
				! Hit the dagger and reflected!
				dagHit = 2
				IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
					CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), &
											(/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
					CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), &
											(/0.0_8, 0.0_8, 1.0_8/))
					state = prevState
				ELSE
					! ERROR FLAG!!!
					dagHit = -1
				END IF
			END IF
		! Now, check if we hit the housing low			
		ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8) .AND. &
				predZ < (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
				ABS(predX + 0.1524_8) < (0.40_8 + 2.0179_8* &
										(predZ+1.5_8- zOff-0.2_8))/2.0_8) THEN
			
			dagHit = 3
			IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
				CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), &
										(/0.0_8, 0.0_8, 1.0_8/))
				state = prevState
			ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
				CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), &
										(/0.0_8, 0.0_8, 1.0_8/))
				state = prevState
			ELSE
				! ERROR FLAG!!!
				dagHit = -1
			END IF
		! Finally, check if we hit the housing high
		ELSE IF (predZ >= (-1.5_8 + zOff + 0.2_8 + 0.14478_8) .AND. &
				predZ < (-1.5_8 + zOff + 0.2_8 + 0.2667_8) .AND. &
				ABS(predX + 0.1524_8) < 0.69215_8/2.0_8) THEN
			
			dagHit = 4	
			IF (prevState(2) > 0 .AND. prevState(5) < 0) THEN
				CALL reflect(prevState, (/0.0_8, 1.0_8, 0.0_8/), &
										(/0.0_8, 0.0_8, 1.0_8/))
				state = prevState
			ELSE IF (prevState(2) < 0 .AND. prevState(5) > 0) THEN
				CALL reflect(prevState, (/0.0_8, -1.0_8, 0.0_8/), &
										(/0.0_8, 0.0_8, 1.0_8/))
				state = prevState
			ELSE
				! ERROR FLAG!!!
				dagHit = -1
			END IF
			
		END IF
	END IF

END SUBROUTINE checkDagHit
! Calculate dagger position
SUBROUTINE zOffDipCalc(t, z, nDips, holdIn, pseType)

	REAL(KIND=PREC), INTENT(IN) :: t
	REAL(KIND=PREC), INTENT(OUT) :: z
	INTEGER, OPTIONAL :: nDips, pseType
	REAL(KIND=PREC), INTENT(IN), OPTIONAL :: holdIn

	INTEGER :: i, pseDips
	REAL(KIND=PREC) :: speed, holdT
	REAL(KIND=PREC), DIMENSION(:), ALLOCATABLE :: dipHeights, dipEnds
	LOGICAL :: normRun
	
	!I've introduced the mapping integer 'pseType' as a flag. Truth table:
	!	CASE 1: Normal running
	!	CASE 2: Deep cleaning (clean at 25 cm) [for phase space evolution]
	!	CASE 3: Height dependent time constant data
	!NB: if pseType is present, then so is holdIn
	IF(PRESENT(pseType)) THEN
		SELECT CASE (pseType)
			CASE (1)
				normRun = .TRUE.
				holdT = holdIn
				pseDips = 0
			CASE (2)
				normRun = .TRUE.
				holdT = holdIn
				pseDips = 2
			CASE (3)
				normRun = .FALSE.
				holdT = holdIn
				pseDips = 0
		END SELECT
	ELSE
		normRun = .TRUE.
		holdT = 900.0
		pseDips = 0
	END IF
		
	IF(.NOT. PRESENT(nDips)) THEN
		nDips = 4 !Default to 4-Dip (which is really 3-Dip!)
	END IF
	
	! Previous switches allow us to determine the number of dips we need. 	
	ALLOCATE(dipHeights(nDips+pseDips))
	ALLOCATE(dipEnds(nDips+pseDips))
			
	! Our heights are hard-coded in, but we can use the previous cases to choose between them.
	IF (normRun) THEN
		! Select number of dips from normal running
		SELECT CASE (nDips)
			CASE (4) !3-Dips from Octet_3dip.txt
				dipHeights((1+pseDips):(nDips+pseDips)) = &
									(/ 0.49_8, 0.380_8, 0.250_8, 0.010_8 /) 
				dipEnds((1+pseDips):(nDips+pseDips)) = &
									(/ 0.0_8, 40.0_8, 60.0_8, 210.0_8 /) 
			CASE (10)!9-Dips from Octet_9dip.txt
				dipHeights((1+pseDips):(nDips+pseDips)) = &
									(/0.49, 0.380, 0.250, 0.180, 0.140, &
										0.110, 0.080, 0.060, 0.040, 0.010/)
				dipEnds((1+pseDips):(nDips+pseDips)) = &
									(/0.0, 40.0, 80.0, 100.0, 120.0, 140.0, &
										160.0, 180.0, 200.0, 300.0/) 
			CASE DEFAULT
				PRINT *, "ERROR: Unknown number of dips in dagger movement! Ending run!"
				STOP
		END SELECT
	ELSE 
		! Select number of dips for time constant runs
		! Automate this??? -- No good easy way, it's actually a weird distribution of heights
		SELECT CASE (nDips)
			CASE (4) 
				dipHeights = (/0.490, 0.380, 0.250, 0.010/)
				dipEnds = (/0.0, 40.0, 400.0, 500.0/)
			CASE (5)
				dipHeights = (/0.490, 0.380, 0.250, 0.180, 0.010/)
				dipEnds = (/0.0, 40.0, 80.0, 400.0, 500.0/)
			CASE (6) 
				dipHeights = (/0.490, 0.380, 0.250, 0.180, 0.140, &
								0.010/)
				dipEnds = (/0.0, 40.0, 80.0, 100.0, 400.0, &
								500.0 /)
			CASE (7)
				dipHeights = (/0.490, 0.380, 0.250, 0.180, 0.140, &
								0.110, 0.010 /)
				dipEnds = (/0.0, 40.0, 80.0, 100.0, 120.0, &
							400.0, 500.0/)
			CASE (8)
				dipHeights = (/0.490, 0.380, 0.250, 0.180, 0.140, &
								0.110, 0.080, 0.010 /)
				dipEnds = (/0.0, 40.0, 80.0, 100.0, 120.0, &
							140.0, 400.0, 500.0/)
			CASE (9)
				dipHeights = (/0.490, 0.380, 0.250, 0.180, 0.140, &
								0.110, 0.080, 0.060, 0.010/)
				dipEnds = (/0.0, 40.0, 80.0, 100.0, 120.0, &
							140.0, 160.0, 400.0, 500.0/)
			CASE (10)
				dipHeights = (/0.490, 0.380, 0.250, 0.180, 0.140, &
								0.110, 0.080, 0.060, 0.040, 0.010/)
				dipEnds = (/0.0, 40.0, 80.0, 100.0, 120.0, &
							140.0, 160.0, 180.0, 400.0, 500.0/)
			CASE (11)
				! Case 11 is weird, since it double drops 0.01 
				dipHeights = (/0.490, 0.380, 0.250, 0.180, 0.140, &
								0.110, 0.080, 0.060, 0.040, 0.010, 0.010/)
				dipEnds = (/0.0, 40.0, 80.0, 100.0, 120.0, &
							140.0, 160.0, 180.0, 200.0, 400.0, 500.0/)
			CASE DEFAULT
				PRINT *, "ERROR: Unknown number of dips in dagger movement! Ending run!"
				STOP
		END SELECT
	END IF
	
	!Check if we're doing deep cleaning
	IF(pseDips .EQ. 2) THEN
		dipHeights(1:2) = (/0.49, 0.250/)
		dipEnds(1:2) = (/0.0_8, 200.0_8/)
		DO i=3,(nDips+pseDips),1
			dipEnds(i) = dipEnds(i)/2 + (200.0+holdT)
		END DO
	END IF
	
	IF (t .GT. dipEnds(nDips+pseDips)) THEN
		z = 0.01
		RETURN
	END IF
	
	DO i=1,(nDips+pseDips),1
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
! General reflection subroutine
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
! General absorption subroutine
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
    zs(3) = (5.76556e-9) !This is the Boron layer thickness
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
! Check block scattering. Set blockScale to 0.0 to turn this off.
SUBROUTINE check_upscatter(prevState, state, blockScale, blockHit, ePerp)

	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	REAL(KIND=PREC), INTENT(IN), DIMENSION(6) :: prevState
	REAL(KIND=PREC), INTENT(INOUT), DIMENSION(6) :: state
	REAL(KIND=PREC), INTENT(IN) :: blockScale
	LOGICAL, INTENT(INOUT) :: blockHit
	REAL(KIND=PREC), INTENT(OUT) :: ePerp
		
	INTEGER :: ii
	REAL(KIND=PREC) :: hitU, upscatterProb
	REAL(KIND=PREC), DIMENSION(6) :: shiftPrevState, shiftState
	REAL(KIND=PREC), DIMENSION(3,3) :: rotation, invRotation
	REAL(KIND=PREC), DIMENSION(3) :: blockSize, blockPos, fracTravel
		
	! These numbers are calculated in a Mathematica notebook.
	! I'm keeping the matrices in here because I don't want to lose them, even if they're slightly unphysical.
	
	!rotation = RESHAPE((/ 0.3052294878125326_8, 0.9441621677380962_8, 0.1240675653899834_8, &
	!				-0.9392750072476546_8, 0.3199526527130549_8, -0.1240675653899834_8, &
	!				-0.1568356481467703_8, -0.07866448394274296_8, 0.9844869112570286_8 /),&
	!			SHAPE(rotation))
	!invRotation = RESHAPE((/ 0.4981846821401102_8, 0.8070369212756123_8, 0.3170227597175609_8, &
	!				-0.7463198129200769_8, 0.5852378206028750_8, -0.3170227597175609_8, &
	!				-0.4413827809553729_8, -0.07866448394274296_8, 0.8938641617394241_8 /),&
	!			SHAPE(rotation))
	invRotation = RESHAPE((/ 0.4981846187446881_8, 0.8070370832294234_8, 0.3170224470581768_8, &
					-0.7463198129200769_8, 0.5852378206028750_8, -0.3170224470581768_8, &
					-0.4413823244195138_8, -0.07866452547073647_8, 0.8938643835182667_8 /),&
				SHAPE(invRotation))
	!invRotation = RESHAPE((/ 0.3052294878125326_8, -0.9392750072476545_8, -0.1568356481467703_8, &
	!				0.9441621677380965_8, 0.3199526527130550_8, -0.07866448394274296_8, &
	!				0.1240675653899834_8, -0.1240675653899834_8, 0.9844869112570285_8 /), &
	!			SHAPE(invRotation))
	!rotation = RESHAPE((/ 0.4981846821401102_8, -0.7463198129200770_8, -0.4413827809553729_8, &
	!				0.8070369212756124_8, 0.5852378206028751_8, -0.07866448394274295_8, &
	!				0.3170227597175610_8, -0.3170227597175610_8, 0.8938641617394240_8 /), &
	!			SHAPE(rotation))
	rotation = RESHAPE((/ 0.4981846187446881_8, -0.7463198129200769_8, -0.4413823244195138_8, &
					0.8070370832294234_8, 0.5852378206028751_8, -0.07866452547073647_8, &
					0.3170227597175610_8, -0.3170224470581768_8, 0.8938643835182667_8 /), &
				SHAPE(rotation))
	!blockPos = (/ 0.2207_8, 0.127_8, -1.4431_8/)
	blockPos = (/ 0.2175882546128424, 0.1264454150954313, -1.436798256096196 /) !block raised by 0.7 mm (most physical)
	!blockPos = (/ 0.2172351487533068, 0.1263824834750547, -1.436083164589382 /) !block raised by 1.5 mm 
	!blockPos = (/ 0.2165730752666775, 0.1262644866848486, -1.434742368014104 /) ! block raised by 3.0 mm
	
	! Size of block -- set with a scaling factor
	blockSize = (/ 0.0127_8, 0.0127_8, 0.00635_8 /)
	DO ii=1,3,1
		blockSize(ii) = blockScale*blockSize(ii)
	END DO
	
	!Shift the neutron position to "block coordinates"
	fracTravel = 0.0_8
	shiftPrevState = (/ MATMUL(rotation, prevState(1:3) - blockPos), MATMUL(rotation,prevState(4:6)) /)
	shiftState = (/ MATMUL(rotation, state(1:3) - blockPos), MATMUL(rotation,state(4:6)) /)
	
	!Check if the neutron is inside the block
	IF (ABS(shiftState(1)) .GT. blockSize(1) .OR. &
			ABS(shiftState(2)) .GT. blockSize(2) .OR. &
			ABS(shiftState(3)) .GT. blockSize(3)) THEN
		state = state
	ELSE
		blockHit = .TRUE.
		! See how far we've traveled through the block. If we're purely inside or outside, this is zero
		fracTravel(1) = (SIGN(blockSize(1), shiftState(1)) - shiftPrevState(1)) / &
						ABS(shiftState(1) - shiftPrevState(1))
		IF (ABS(fracTravel(1)) > 1) THEN
			fracTravel(1) = 0.0_8
		END IF
		fracTravel(2) = (SIGN(blockSize(2), shiftState(2)) - shiftPrevState(2)) / &
						ABS(shiftState(2) - shiftPrevState(2))
		IF (ABS(fracTravel(2)) > 1) THEN
			fracTravel(2) = 0.0_8
		END IF
		fracTravel(3) = (SIGN(blockSize(3), shiftState(3)) - shiftPrevState(3)) / &
						ABS(shiftState(3) - shiftPrevState(3))
		IF (ABS(fracTravel(3)) > 1) THEN
			fracTravel(3) = 0.0_8
		END IF		
		
		!If we don't upscatter, reflect off the block.
		IF ((ABS(fracTravel(3)) > ABS(fracTravel(1))) .AND. &
				(ABS(fracTravel(3)) > ABS(fracTravel(2)))) THEN
			ePerp = shiftPrevState(6)*shiftPrevState(6)/(2.0_8*MASS_N)
			CALL reflect(shiftPrevState, (/ 0.0_8, 0.0_8, 1.0_8 /), (/ 1.0_8, 0.0_8, 0.0_8 /))
			state = (/ MATMUL(invRotation, shiftPrevState(1:3)) + blockPos, &
					MATMUL(invRotation,shiftPrevState(4:6)) /)
			!PRINT *," Reflection in +z "
		ELSE IF ((ABS(fracTravel(1)) > ABS(fracTravel(2))) .AND. &
				(ABS(fracTravel(1)) > ABS(fracTravel(3))) .AND. &
				(fracTravel(1) < 0)) THEN
			ePerp = shiftPrevState(4)*shiftPrevState(4)/(2.0_8*MASS_N)
			CALL reflect(shiftPrevState, (/ 1.0_8, 0.0_8, 0.0_8 /), (/ 0.0_8, 0.0_8, 1.0_8 /))
			state = (/ MATMUL(invRotation, shiftPrevState(1:3)) + blockPos, &
					MATMUL(invRotation,shiftPrevState(4:6)) /)
			!PRINT *," Reflection in +x "
		ELSE IF ((ABS(fracTravel(1)) > ABS(fracTravel(2))) .AND. &
				(ABS(fracTravel(1)) > ABS(fracTravel(3))) .AND. &
				(fracTravel(1) > 0)) THEN
			ePerp = shiftPrevState(4)*shiftPrevState(4)/(2.0_8*MASS_N)
			CALL reflect(shiftPrevState, (/ -1.0_8, 0.0_8, 0.0_8 /), (/ 0.0_8, 0.0_8, 1.0_8 /))
			state = (/ MATMUL(invRotation, shiftPrevState(1:3)) + blockPos, &
					MATMUL(invRotation,shiftPrevState(4:6)) /)
			!PRINT *," Reflection in -x "
		ELSE IF ((ABS(fracTravel(2)) > ABS(fracTravel(1))) .AND. &
				(ABS(fracTravel(2)) > ABS(fracTravel(3))) .AND. &
				(fracTravel(2) < 0)) THEN
			ePerp = shiftPrevState(5)*shiftPrevState(5)/(2.0_8*MASS_N)
			CALL reflect(shiftPrevState, (/ 0.0_8, 1.0_8, 0.0_8 /), (/ 0.0_8, 0.0_8, 1.0_8 /))
			state = (/ MATMUL(invRotation, shiftPrevState(1:3)) + blockPos, &
					MATMUL(invRotation,shiftPrevState(4:6)) /)
			!PRINT *," Reflection in +y "
		ELSE IF ((ABS(fracTravel(2)) > ABS(fracTravel(1))) .AND. &
				(ABS(fracTravel(2)) > ABS(fracTravel(3))) .AND. &
				(fracTravel(2) > 0)) THEN
			ePerp = shiftPrevState(5)*shiftPrevState(5)/(2.0_8*MASS_N)
			CALL reflect(shiftPrevState, (/ 0.0_8, -1.0_8, 0.0_8 /), (/ 0.0_8, 0.0_8, 1.0_8 /))
			state = (/ MATMUL(invRotation, shiftPrevState(1:3)) + blockPos, &
					MATMUL(invRotation,shiftPrevState(4:6)) /)
			!PRINT *," Reflection in -y "
		ELSE 
			CALL reflect(shiftPrevState, (/ 0.0_8, 0.0_8, -1.0_8 /), (/ 1.0_8, 0.0_8, 0.0_8 /))
			state = (/ MATMUL(invRotation, shiftPrevState(1:3)) + blockPos, &
					MATMUL(invRotation,shiftPrevState(4:6)) /) 
			blockHit = .FALSE. !Managed to sneak under the block, don't record it...
			!PRINT *, "UHOH"
		END IF		
	END IF
END SUBROUTINE check_upscatter

!-----------------------------------------------------------------------
!------- MAIN TRACKING PROGRAM! (PUT ALL FUNCTIONALITY IN HERE) --------
!-----------------------------------------------------------------------
SUBROUTINE trackDaggerFull(state, holdT, blockScale, nDips,pse, trX,trY,trZ,nTr)

	USE symplecticInt
	USE constants
	USE forcesAndPotential
	USE convertTrace
	IMPLICIT NONE
	
	! The only OPTIONAL input is the heating. Everything else should be auto-set.
	REAL(KIND=PREC), DIMENSION(6), INTENT(INOUT) :: state
	REAL(KIND=PREC), INTENT(IN) :: holdT, blockScale
	INTEGER, INTENT(IN) :: nDips, pse
	REAL(KIND=PREC), INTENT(IN), ALLOCATABLE, &
						DIMENSION(:), OPTIONAL :: trX(:), trY(:), trZ(:)
	INTEGER, INTENT(IN), OPTIONAL :: nTr
	
	INTEGER :: i, numSteps, nHit, nHitHouseLow, nHitHouseHigh, nBlockHit, dagHit
	LOGICAL :: exitFlag
	REAL(KIND=PREC) :: t, fracTravel, predX, predY, predZ, energy, zOff, zeta
	REAL(KIND=PREC) :: trapFill, beamFill, cleaningTime, settlingTime, deathTime
	REAL(KIND=PREC) :: cleanZ, countingTime, deepClean, absProb, rNum, ePerp
	REAL(KIND=PREC), DIMENSION(6) :: prevState, iniState
		
	!Initialize our variables, including our initial state (which is saved)
	nHit = 0
	nHitHouseLow = 0
	nHitHouseHigh = 0
	nBlockHit = 0
	t = 0.0_8
	exitFlag = .FALSE.
	dagHit = 0
	iniState = state
	
	! Constants for timing and heights and stuff
	trapFill = 70.0_8 !Trap filling time constant -- will change w/ roundhouse
	beamFill = 150.0_8 !Amount of time beam is on/trap door is open
	cleaningTime = 50.0_8 ! Cleaning time. Might make this (another!) variable
	cleanZ = 0.380_8
	deathTime = 3000.0_8 ! Max time of calculation. 
	deepClean = 200.0_8 ! Time that deep cleaning takes
			
	! Choose our required settling time -- time neutron is in trap during filling
	DO 
		CALL RANDOM_NUMBER(rNum)
		settlingTime = -trapFill*LOG(rNum)
		IF (settlingTime < beamFill) THEN
			EXIT
		END IF
	END DO
			
	! Start with filling time:
	numSteps = settlingTime/dt
	t = (beamFill - settlingTime)
	
	DO i=1,numsteps,1
		! Catch to end run
		IF (exitFlag) THEN
			t = t + (numSteps- i)*dt 
			EXIT
		END IF
		
		prevState = state
		! Step either with or without heating
		IF (PRESENT(nTr)) THEN
			CALL symplecticStep(state, dt, energy, t, trX, trY, trZ,nTr)
		ELSE
			CALL symplecticStep(state, dt, energy, t)
		END IF
		! Check for block scatter
		IF (blockScale .NE. 0) THEN
			CALL check_upscatter(prevState, state, blockScale, exitFlag, ePerp)
			! If we hit the block, use the bool to write and continue
			IF (exitFlag) THEN
				nBlockHit = nBlockHit + 1
				WRITE(2) t, energy, nHit, nHitHouseLow, nHitHouseHigh, &
					nBlockHit, state(1), state(2), state(3), ePerp, &
					iniState(4), iniState(5), iniState(6), settlingTime
				exitFlag = .FALSE.
			END IF
		END IF
		
		! Check cleaning -- both dagger (first) and cleaner (second).
		! No dagger hit will be recorded, just flagged 
		CALL checkDagHit(state, prevState, cleanZ, dagHit)
		IF (dagHit .EQ. 1) THEN
			dagHit = 1
			exitFlag = .TRUE.
		ELSE
			! Reset: we don't care if the UCN was reflected or anything.
			dagHit = 0
			CALL cleaning(state, prevState, cleanZ, exitFlag)
		END IF
		
	END DO
	! Now do cleaning time -- same as settling time but without any dagger hits
	numSteps = cleaningTime/dt
	DO i=1,numSteps,1
		! Catch to end run
		IF (exitFlag) THEN
			t = t + (cleaningTime - i*dt) 
			EXIT
		END IF
		
		prevState = state
		! Step either with or without heating
		IF (PRESENT(nTr)) THEN
			CALL symplecticStep(state, dt, energy, t, trX, trY, trZ,nTr)
		ELSE
			CALL symplecticStep(state, dt, energy, t)
		END IF
		! Check for block scatter
		IF (blockScale .NE. 0) THEN
			CALL check_upscatter(prevState, state, blockScale, exitFlag, ePerp)
			! If we hit the block, use the bool to write and continue
			IF (exitFlag) THEN
				nBlockHit = nBlockHit + 1
				WRITE(2) t, energy, nHit, nHitHouseLow, nHitHouseHigh, &
					nBlockHit, state(1), state(2), state(3), ePerp, &
					iniState(4), iniState(5), iniState(6), settlingTime
				exitFlag = .FALSE.
			END IF
		END IF
		! Check if we hit the cleaner. Dagger is out of the way now.
		CALL cleaning(state, prevState, cleanZ, exitFlag)
		
	END DO
	
	! This is the holding time -- no dagger and no cleaner. 
	! Only happens in normal running (no deep cleaning)

	! I've introduced the mapping integer 'pse' as a flag. Truth table:
	!	CASE 1: Normal running
	!	CASE 2: Deep cleaning (clean at 25 cm) [for phase space evolution]
	!	CASE 3: Height dependent time constant data
	! We'll do a normal holding time if we DON'T have a Deep Cleaning pse run.
	IF ((pse .EQ. 1) .OR. (pse .EQ. 3)) THEN
		numSteps = holdT/dt
		DO i=1,numSteps,1
			! Catch to end run -- basically if we got cleaned already
			IF (exitFlag) THEN
				t = t + (holdT - i*dt) 
				EXIT
			END IF
			
			prevState = state
			! Step either with or without heating
			IF (PRESENT(nTr)) THEN
				CALL symplecticStep(state, dt, energy, t, trX, trY, trZ,nTr)
			ELSE
				CALL symplecticStep(state, dt, energy, t)
			END IF
			
			IF (blockScale .NE. 0) THEN
				CALL check_upscatter(prevState, state, blockScale, exitFlag, ePerp)
				! If we hit the block, use the bool to write and continue
				IF (exitFlag) THEN
					nBlockHit = nBlockHit + 1
					WRITE(2) t, energy, nHit, nHitHouseLow, nHitHouseHigh, &
						nBlockHit, state(1), state(2), state(3), ePerp, &
						iniState(4), iniState(5), iniState(6), settlingTime
					exitFlag = .FALSE.
				END IF
			END IF
			
		END DO
		! Offset for our counting 
		countingTime = cleaningTime + beamFill + holdT
		
	ELSE IF (pse .EQ. 2) THEN
		! Deep cleaning has dagger movement before hold
		countingTime = cleaningTime + beamFill
	END IF
			
	! Now we do our unload, with the dagger moving.
	DO
		! Check (before running) if we're already done.
		IF (exitFlag) THEN
			EXIT
		END IF
		
		prevState = state
		! Step either with or without heating
		IF (PRESENT(nTr)) THEN
			CALL symplecticStep(state, dt, energy, t, trX, trY, trZ,nTr)
		ELSE
			CALL symplecticStep(state, dt, energy, t)
		END IF
		
		IF (blockScale .NE. 0) THEN
			!Only way we can make it this far is if we didn't tag exitFlag
			CALL check_upscatter(prevState, state, blockScale, exitFlag, ePerp)
			! If we hit the block, use the bool to write and continue
			IF (exitFlag) THEN
				nBlockHit = nBlockHit + 1
				WRITE(2) t, energy, nHit, nHitHouseLow, nHitHouseHigh, &
					nBlockHit, state(1), state(2), state(3), ePerp, &
					iniState(4), iniState(5), iniState(6), settlingTime
				exitFlag = .FALSE.
			END IF
		END IF
	
		! Put the cleaner in for deep cleaning
		IF ((pse .EQ. 2) .AND. ((t-countingTime) .LT. deepClean)) THEN
			CALL cleaning(state, prevState, cleanZ, exitFlag)
		END IF
		
		! Calculate the dagger height
		CALL zOffDipCalc(t - countingTime, zOff, nDips, holdT, pse)
		! Check if we hit the dagger		
		CALL checkDagHit(state, prevState, zOff, dagHit)
		
		! dagHit is a switch variable for position the neutron hits. Truth table:
		!  -1: ERROR!!!
		! 	0: No hits
		! 	1: Hits the dagger, absorbed
		!	2: Hits the dagger, reflected
		! 	3: Hits the low housing
		!	4: Hits the high housing
		! I've copied this truth table into the checkDagHit subroutine for posterity.
		SELECT CASE (dagHit)
			CASE (-1)
				PRINT *, "ERROR: Something weird happened with the dagger"	
				EXIT
			CASE (1)
				nHit = nHit + 1
				! Calculate dagger hit components
				fracTravel = ABS(prevState(2))/(ABS(state(2)) + ABS(prevState(2)))
				predX = prevState(1) + fracTravel * (state(1) - prevState(1))
				predY = prevState(2) + fracTravel * (state(2) - prevState(2))
				predZ = prevState(3) + fracTravel * (state(3) - prevState(3))
				ePerp = state(5)*state(5)/(2.0_8*MASS_N)
				WRITE(1) t - countingTime, energy, ePerp, predX, predY, predZ, &
					zOff, nHit, nHitHouseLow, nHitHouseHigh, nBlockHit, &
					iniState(4), iniState(5), iniState(6), settlingTime
				exitFlag = .TRUE.
			CASE (2)
				nHit = nHit + 1
			CASE (3)
				nHitHouseLow = nHitHouseLow + 1
			CASE (4)
				nHitHouseHigh = nHitHouseHigh + 1
			CASE DEFAULT
				CONTINUE
		END SELECT
				
		IF (t-settlingTime > deathTime) THEN
			exitFlag = .TRUE.
		END IF
		
	END DO
	
END SUBROUTINE trackDaggerFull

!-----------------------------------------------------------------------
!-----------------Subroutines for debug purposes------------------------
!-----------------------------------------------------------------------
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
			CALL symplecticStep(state, dt, energy, sympT)
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
		CALL symplecticStep(state, dt, eEnd, sympT)
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
	
	CALL calcEnergy(state, energy)
	energy = 0.0

	PRINT *, 0.0_8, state(1), state(2), state(3),&
			state(4)/MASS_N, state(5)/MASS_N, state(6)/MASS_N, energy,&
			totalU
	
	DO i=1,numSteps,1
!		IF(t > 475.0_8 .AND. t < 480.0_8) THEN
		IF(INT(dt*100_8*i)-INT(dt*100_8*(i-1)) .NE. 0) THEN
!		IF(1 .EQ. 1) THEN
!			energy = totalU + SUM(state(1,4:6)**2)/(2.0_8*MASS_N)

			PRINT *, dt*i, state(1), state(2), state(3),&
			state(4)/MASS_N, state(5)/MASS_N, state(6)/MASS_N, energy,&
			totalU!, fx, fy, fz

			!PRINT *, dt*i, energy
		END IF
		IF(PRESENT(sympT)) THEN
			CALL symplecticStep(state, dt, energy, t)
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
