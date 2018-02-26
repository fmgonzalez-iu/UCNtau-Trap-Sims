#include "constants.h"

MODULE testSubroutines

CONTAINS

SUBROUTINE compPots()
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC) :: x, y, z, fx, fy, fz, totalU
	integer :: i
	
	z = -1.49_8
	x = 0.1_8
	y = 0.1_8
	DO i=0,20,1
		CALL totalForce(x, y, z + i*(0.01_8/20.0_8), fx, fy, fz, totalU)
		PRINT *, totalU, fx, fy, fz
	END DO
END SUBROUTINE compPots

SUBROUTINE zOffDipCalc(t, z)
	real(kind=PREC), intent(in) :: t
	real(kind=PREC), intent(out) :: z
	
	integer :: nDips = 4
	real(kind=PREC) :: speed
	real(kind=PREC), dimension(4) :: dipHeights
	real(kind=PREC), dimension(4) :: dipEnds
	
	integer :: i
	
!	dipHeights = (/0.49, 0.380, 0.250, 0.180, 0.140, 0.110, 0.080, 0.060, 0.040, 0.010/)
	dipHeights = (/0.49_8, 0.380_8, 0.250_8, 0.010_8 /)
!	dipEnds =     (/0.0,  40.0,  80.0,  100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 300.0/)
	dipEnds =     (/0.0_8,  40.0_8, 60.0_8,  210.0_8 /)
	
	IF (t > dipEnds(nDips)) THEN
		z = 0.01
		RETURN
	END IF
	
	speed = 0.49_8/13.0_8
	
	DO i=1,nDips,1
		IF (dipEnds(i) > t) THEN
			EXIT
		END IF
	END DO
	
	z = dipHeights(i-1) - speed*(t-dipEnds(i-1))
	
	IF (z < dipHeights(i)) THEN
		z = dipHeights(i)
	END IF
END SUBROUTINE zOffDipCalc

SUBROUTINE reflect(state, norm, tang)
	USE trackGeometry
	USE constants
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC), dimension(3), intent(in) :: norm, tang
	real(kind=PREC) :: u1, u2, theta, phi, pN, pT, pTprime, pLen, pTarget
	real(kind=PREC), dimension(3) :: tangPrime, newPdir
	
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

SUBROUTINE trackDaggerHitTime(state)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC), dimension(6) :: prevState

	real(kind=PREC) :: t, fracTravel, predX, predZ, energy, zOff, zeta
	real(kind=PREC) :: settlingTime
	real(kind=4), dimension(50) :: hitT
	real(kind=4), dimension(50) :: hitE
	
	integer :: i, numSteps, nHit
	
	nHit = 0
	!PRINT *, "Started dagger hit time function"
	hitT = 0.0_8
	hitE = 0.0_8
	
	t = 0.0_8
	
	settlingTime = 20.0_8 + 50.0_8
	
	numSteps = settlingTime/dt
	DO i=1,numSteps,1
		CALL symplecticStep(state, dt, energy)
		t = t + dt
	END DO
	
	!PRINT *, "got to end of settling time!"
	
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
			IF (ABS(predX) < .2 .AND. zeta > 0.0_8 .AND. predZ < (-1.5_8 + zOff + 0.2_8)) THEN
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
!				WRITE(1) t - (20.0_8 + 50.0_8), predX, predZ - zOff
!				EXIT
			END IF
			
			IF (t > 2000) THEN
				EXIT
			END IF
		END IF
	END DO
	!PRINT *, "Got to end of dagger detection!"
	
	WRITE(1) energy, hitT, hitE
END SUBROUTINE trackDaggerHitTime

SUBROUTINE trackDaggerHitTimeFixedEff(state)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC), dimension(6) :: prevState

	real(kind=PREC) :: t, fracTravel, predX, predZ, energy, zOff, zeta, hitU
	real(kind=PREC) :: settlingTime

	integer :: i, numSteps, nHit
	
	nHit = 0
	
	t = 0.0_8
	
	settlingTime = 20.0_8 + 50.0_8
	
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
			IF (ABS(predX) < .2 .AND. zeta > 0.0_8 .AND. predZ < (-1.5_8 + zOff + 0.2_8)) THEN
                CALL RANDOM_NUMBER(hitU)
                IF (hitU < 0.155) THEN
                    WRITE(1) t - settlingTime, energy, state(5)*state(5)/(2.0_8*MASS_N)
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
!				WRITE(1) t - (20.0_8 + 50.0_8), predX, predZ - zOff
!				EXIT
			END IF
			
			IF (t > 2000) THEN
				EXIT
			END IF
		END IF
	END DO
END SUBROUTINE trackDaggerHitTimeFixedEff

SUBROUTINE trackDaggerAndBlock(state)

	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC), dimension(6) :: prevState

	real(kind=PREC) :: t, fracTravel, predX, predZ, energy, zOff, zeta, hitU
	real(kind=PREC) :: settlingTime
	real(kind=4), dimension(50) :: hitT
	real(kind=4), dimension(50) :: hitE
	
	integer :: i, numSteps, nHit
	LOGICAL :: blockHit, dagHit

	nHit = 0
		
	hitT = 0.0_8
	hitE = 0.0_8
	
	t = 0.0_8
	blockHit = .FALSE.	
	dagHit = .FALSE.

	settlingTime = 20.0_8 + 50.0_8
!	PRINT *, "loaded dagger and block!"
!	PRINT *, state
	numSteps = settlingTime/dt
	
	DO i=1,numSteps,1
		IF (blockHit) THEN
			EXIT
		END IF
		CALL symplecticStep(state, dt, energy)
		t = t + dt
		CALL check_upscatter(state, t, energy, blockHit)
	END DO
	
!	PRINT *, "Got to end of settling period detection!"
	DO
		IF (blockHit) THEN
			EXIT
		END IF
		prevState = state
		CALL symplecticStep(state, dt, energy)
		t = t + dt
		CALL check_upscatter(state, t, energy, blockHit) 
		IF (blockHit) THEN 
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
			IF (ABS(predX) < .2 .AND. zeta > 0.0_8 .AND. predZ < (-1.5_8 + zOff + 0.2_8)) THEN
				nHit = nHit + 1
				hitT(nHit) = t - settlingTime
				hitE(nHit) = state(5)*state(5)/(2.0_8*MASS_N)
				IF (nHit .EQ. 50) THEN
					dagHit = .TRUE.
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
			END IF
			
			IF (t > 2000) THEN
				dagHit = .TRUE.
				EXIT
			END IF
		END IF
		
	END DO
	IF (dagHit) THEN	
		PRINT *, "Got to end of dagger detection!"
		WRITE(1) energy, hitT, hitE
	END IF
END SUBROUTINE trackDaggerAndBlock

SUBROUTINE trackEnergyGain(state, energy_start, energy_end, sympT, freq)
	USE symplecticInt
	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC), intent(out) :: energy_start, energy_end
	real(kind=PREC), optional, intent(inout) :: sympT
	real(kind=PREC), optional, intent(in) :: freq
    real(kind=PREC) :: eta, prevEta, t, energy, t_end
	integer :: i, numSteps, triggered, rising
    
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
		IF(present(sympT)) THEN
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
	real(kind=PREC), intent(in) :: freq, height
	real(kind=PREC), intent(inout) :: sympT
	real(kind=PREC), intent(out) :: eStart, eEnd

	real(kind=PREC), dimension(6) :: state
	integer :: i, numSteps
	
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
	real(kind=PREC), dimension(6), intent(inout) :: state
	real(kind=PREC), optional, intent(in) :: sympT
	
	real(kind=PREC) :: pr, rdot, r, pphi, phidot, phi, ptheta, thetadot, theta,&
		ldot, l, phiOffset, prevPhi, totalU, totalKE, t, energy, &
		fx_dan, fy_dan, fz_dan, e_dan, fx_nate, fy_nate, fz_nate, e_nate
	integer :: i, numSteps, numPoints, modulus
		
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
		IF(present(sympT)) THEN
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
	real(kind=PREC) :: x, y, z, fx, fy, fz, totalUplus
	real(kind=PREC) :: x0, y0, z0
	integer :: xIt, yIt, zIt, aIt
	
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

SUBROUTINE check_upscatter(state, t, energy, blockHit)

	USE constants
	USE forcesAndPotential
	IMPLICIT NONE
	
	real(kind=PREC), intent(inout), dimension(6) :: state
	real(kind=PREC), dimension(6) :: shiftState
		
	real(kind=PREC) :: t, energy
	real(kind=PREC) :: hitU, upscatterProb
		
	real(kind=PREC), dimension(3,3) :: rotation, invRotation
	real(kind=PREC), dimension(3) :: blockSize
	real(kind=PREC), dimension(3) :: blockPos
	LOGICAL :: blockHit
	
	upscatterProb = 1.0_8
	rotation = RESHAPE((/ 0.305229_8, 0.944162_8, 0.124068_8, -0.939275_8, &
				0.319953_8, -0.124068_8, -0.156836_8, -0.0786645_8, 0.984487_8 /), &
				SHAPE(rotation))
	invRotation = RESHAPE((/ 0.305229_8, -0.939275_8, -0.156836_8, 0.944162_8, &
				0.319953_8, -0.0786645_8, 0.124068_8, -0.124068_8, 0.984487_8 /), &
				SHAPE(invRotation))
	blockSize = (/ 0.0125_8, 0.0125_8, 0.00625_8 /)
	blockPos = (/ 0.2207_8, 0.127_8, -1.4431_8/)
	
	!Shift the neutron position to "block coordinates"
	shiftState = (/ MATMUL(rotation, state(1:3) - blockPos), state(4:6) /)
	!Check if the neutron is inside the block
	IF (ABS(shiftState(1)) .GT. blockSize(1)) THEN
		state = state
	ELSE IF (ABS(shiftState(2)) .GT. blockSize(2)) THEN
		state = state
	ELSE IF (ABS(shiftState(3)) .GT. blockSize(3)) THEN
		state = state
	ELSE
		PRINT *, "Block scattering!"
		CALL RANDOM_NUMBER(hitU)
		IF (hitU < upscatterProb) THEN 
			WRITE(2) t, energy, state(5)*state(5)/(2.0_8*MASS_N)
			WRITE(2) state(1), state(2), state(3)
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
	!Same but in y
!	ELSE IF (ABS(shiftState(2)) .LT. blockSize(2)) THEN
!		!PRINT *, "Block scattering in y!"
!		CALL RANDOM_NUMBER(hitU)
!		IF (hitU < upscatterProb) THEN 
!			WRITE(2) t, energy, state(5)*state(5)/(2.0_8*MASS_N)	
!			STOP
!		ELSE IF (shiftState(2) .GT. 0 .AND. shiftState(5) .LT. 0) THEN 
!			CALL reflect(shiftState, (/ 0.0_8, 1.0_8, 0.0_8 /), (/ 0.0_8, 0.0_8, 1.0_8 /))
!			state = (/ MATMUL(invRotation, shiftState(1:3)) + blockPos, shiftState(4:6) /)
!		ELSE IF (shiftState(2) .LT. 0 .AND. shiftState(5) .GT. 0) THEN
!			CALL reflect(shiftState, (/ 0.0_8, -1.0_8, 0.0_8 /), (/ 0.0_8, 0.0_8, 1.0_8 /))
!			state = (/ MATMUL(invRotation, shiftState(1:3)) + blockPos, shiftState(4:6) /)
!		ELSE 
!			PRINT *, "UHOH"
!		END IF
!	!Also z, but there's only one possible z coord
!	ELSE IF (ABS(shiftState(3)) .LT. blockSize(3)) THEN
!		!PRINT *, "Block scattering in z!"
!		CALL RANDOM_NUMBER(hitU)
!		IF (hitU < upscatterProb) THEN 
!			WRITE(2), t, energy, state(5)*state(5)/(2.0_8*MASS_N)
!			STOP
!		ELSE IF (shiftState(3) .GT. 0 .AND. shiftState(6) .LT. 0) THEN 
!			CALL reflect(shiftState, (/ 0.0_8, 0.0_8, 1.0_8 /), (/ 1.0_8, 0.0_8, 0.0_8 /))
!			state = (/ MATMUL(invRotation, shiftState(1:3)) + blockPos, shiftState(4:6) /)
!		ELSE 
!			PRINT *, "UHOH"
!		END IF
!	END IF
	
END SUBROUTINE check_upscatter

END MODULE
