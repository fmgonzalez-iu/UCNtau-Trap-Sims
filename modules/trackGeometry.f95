#include "constants.h"

MODULE trackGeometry

CONTAINS

SUBROUTINE randomPointTrap(x,y,z,px,py,pz, minE,maxE,dist)

	USE forcesAndPotential
	USE constants
	IMPLICIT NONE
	
	REAL(KIND=PREC), INTENT(OUT) :: x, y, z, px, py, pz
	REAL(KIND=PREC), INTENT(IN) :: minE, maxE
	INTEGER, INTENT(IN) :: dist
	
	REAL(KIND=PREC) :: minEnergy, maxEnergy, energy, totalU, & !Used for energy
					target_p, p_len, max_p !Used for momentum scaling
	REAL(KIND=PREC) :: u1, u2, phi, theta !Used for lambertian generation of track directions
	
	! Scale our input max/min to get our energies
	minEnergy = GRAV*MASS_N*minE
	maxEnergy = GRAV*MASS_N*maxE
	max_p = SQRT(2.0_8*MASS_N*maxEnergy)
	
	! Generate particle above energy threshold
	DO
		CALL RANDOM_NUMBER(energy)
		energy = maxEnergy * energy
		IF (energy > minEnergy) THEN
			EXIT
		END IF
	END DO

	! Generate particle positon on trap door
	DO
		! This is the minimum potential point
		z = -1.464413669130002_8
		CALL RANDOM_NUMBER(x)
		x = x*0.15_8 - 0.075_8
		CALL RANDOM_NUMBER(y)
		y = y*0.15_8 - 0.075_8
		! Make sure said particle can reach this value
		CALL totalPotential(x, y, z, totalU)
		IF (totalU < energy) THEN
			EXIT
		END IF
	END DO
		
	target_p = SQRT(2.0_8*MASS_N*(energy - totalU))
	
	! Set particle's direction
	CALL RANDOM_NUMBER(u1)
	CALL RANDOM_NUMBER(u2)
	SELECT CASE (dist)
		CASE (1)
			!cos(theta)*sin(theta) distribution
			theta = ASIN(SQRT(u1))
			phi = 2.0_8 * PI * u2
			
			px = SIN(theta)*COS(phi)
			py = SIN(theta)*SIN(phi)
			pz = COS(theta)
			
		CASE (2)
			!Isotropic emission
			u1 = u1 * 2 - 1.0_8
			phi = u2 * 2.0_8 * PI
			
			px = SQRT(1-u1*u1)*COS(phi)
			py = SQRT(1-u1*u1)*SIN(phi)
			pz = u1
		
		CASE (3) 
			!Flat in theta emission -- DEFAULT from test_C_eval (for now)
			theta = u1 * PI / 2.0_8
			phi = 2.0_8 * PI * u2

		    px = SIN(theta)*COS(phi)
		    py = SIN(theta)*SIN(phi)
		    pz = COS(theta)
		
		CASE DEFAULT
			PRINT *, "ERROR: Something happened with generating distribution tracks!"
			STOP
	
	END SELECT
    
    ! Scale momentum by energy    
    p_len = SQRT(px*px + py*py + pz*pz)

    px = (target_p/p_len)*px
    py = (target_p/p_len)*py
    pz = (target_p/p_len)*pz
!    PRINT *, x, y, z, px, py, pz, (px*px + py*py + pz*pz)/(2.0_8*MASS_N), totalU, theta
END SUBROUTINE randomPointTrap

SUBROUTINE randomPointTrapOpt(x,y,z,px,py,pz, minE,maxE,dist)

	USE forcesAndPotential
	USE constants
	IMPLICIT NONE
	
	REAL(KIND=PREC), INTENT(OUT) :: x, y, z, px, py, pz
	REAL(KIND=PREC), INTENT(IN) :: minE, maxE
	INTEGER, INTENT(IN) :: dist
	
	REAL(KIND=PREC) :: minEnergy, maxEnergy, energy, totalU, & !Used for energy
					target_p, p_len, max_p !Used for momentum scaling
	REAL(KIND=PREC) :: u1, u2, phi, theta !Used for lambertian generation of track directions
	REAL(KIND=PREC) :: eScale, thScale, eInt
	
	! Scaling factors hard coded in here to reduce errors
	eScale = 1.2
	thScale = 1.28
	
	! Scale our input max/min to get our energies
	minEnergy = GRAV*MASS_N*minE
	maxEnergy = GRAV*MASS_N*maxE
	max_p = SQRT(2.0_8*MASS_N*maxEnergy)
	
	! Generate particle above energy threshold
	DO
		CALL RANDOM_NUMBER(eInt)
		energy = maxEnergy * (eInt ** eScale)
		IF (energy > minEnergy) THEN
			EXIT
		END IF
	END DO

	! Generate particle positon on trap door
	DO
		! This is the minimum potential point
		z = -1.464413669130002_8
		CALL RANDOM_NUMBER(x)
		x = x*0.15_8 - 0.075_8
		CALL RANDOM_NUMBER(y)
		y = y*0.15_8 - 0.075_8
		! Make sure said particle can reach this value
		CALL totalPotential(x, y, z, totalU)
		IF (totalU < energy) THEN
			EXIT
		END IF
	END DO
		
	target_p = SQRT(2.0_8*MASS_N*(energy - totalU))
	
	! Set particle's direction
	CALL RANDOM_NUMBER(u1)
	CALL RANDOM_NUMBER(u2)
	SELECT CASE (dist)
		CASE (1)
			!cos(theta)*sin(theta) distribution
			theta = ASIN(SQRT(u1))
			phi = 2.0_8 * PI * u2
			
			px = SIN(theta)*COS(phi)
			py = SIN(theta)*SIN(phi)
			pz = COS(theta)**thScale
			
		CASE (2)
			!Isotropic emission
			u1 = u1 * 2 - 1.0_8
			phi = u2 * 2.0_8 * PI
			
			px = SQRT(1-u1*u1)*COS(phi)
			py = SQRT(1-u1*u1)*SIN(phi)
			pz = u1
		
		CASE (3) 
			!Flat in theta emission -- DEFAULT from test_C_eval (for now)
			theta = u1 * PI / 2.0_8
			phi = 2.0_8 * PI * u2

		    px = SIN(theta)*COS(phi)
		    py = SIN(theta)*SIN(phi)
		    pz = COS(theta)
		
		CASE DEFAULT
			PRINT *, "ERROR: Something happened with generating distribution tracks!"
			STOP
	
	END SELECT
    
    ! Scale momentum by energy    
    p_len = SQRT(px*px + py*py + pz*pz)

    px = (target_p/p_len)*px
    py = (target_p/p_len)*py
    pz = (target_p/p_len)*pz
!    PRINT *, x, y, z, px, py, pz, (px*px + py*py + pz*pz)/(2.0_8*MASS_N), totalU, theta
END SUBROUTINE randomPointTrapOpt

!
SUBROUTINE randomPointTrapOptimum(x,y,z,px,py,pz)
    USE forcesAndPotential
    USE constants
    IMPLICIT NONE
    REAL(KIND=PREC), INTENT(OUT) :: x, y, z, px, py, pz
    REAL(KIND=PREC) :: zeta, maj_r, min_r, totalU, energy, max_p, maxEnergy, &
                    target_p, p_reject_val, p_len, en_reject_val
    REAL(KIND=PREC) :: u1, u2, phi, theta !Used for lambertian generation of track directions
    ! These maxenergy, minenergy things are kept from previous runs. With post-processing this
    ! is less necessary, so I'm just dumping a bunch of commented numbers here
    !	maxEnergy = GRAV*MASS_N*0.38_8 !Assume cleaning height of 38cm
!	maxEnergy = GRAV*MASS_N*0.345_8 !Assume cleaning height of 38cm - 3.5cm to the zero energy point
!	maxEnergy = GRAV*MASS_N*0.3444136691300193_8 !Assume Nathan's calc is 100% right
!	maxEnergy = GRAV*MASS_N*0.45_8 !ASSUME NOTHING! CLEAN IT UP YOURSELF!
    !A UCN of energy 34.5cm can reach 38cm w.r.t. the bottom of the trap.
    !maxEnergy = GRAV*MASS_N*0.345_8 !Assume cleaning height of 38cm - 3.5cm to the zero energy point
    maxEnergy = GRAV*MASS_N*0.3444136691300193_8 !More precise maxEnergy using more precise zmin.
    !A UCN of energy 34.5cm can reach 38cm w.r.t. the bottom of the trap.
    max_p = SQRT(2.0_8*MASS_N*maxEnergy)
    
    DO
        CALL RANDOM_NUMBER(energy)
        energy = maxEnergy * energy
!        IF (energy < 11.7/JTONEV) THEN
		IF (energy < 11.641026_8/JTONEV) THEN
            CYCLE
        END IF
        
        CALL RANDOM_NUMBER(en_reject_val)
        IF (en_reject_val < energy/maxEnergy) THEN
            EXIT
        END IF
    END DO

    DO
        z = -1.464413669130002_8
        CALL RANDOM_NUMBER(x)
        x = x*0.15_8 - 0.075_8
        CALL RANDOM_NUMBER(y)
        y = y*0.15_8 - 0.075_8

        CALL totalPotential(x, y, z, totalU)
        IF (totalU < energy) THEN
            EXIT
        END IF
    END DO
        
    target_p = SQRT(2.0_8*MASS_N*(energy - totalU))

    !cos(theta)*sin(theta) distribution
    CALL RANDOM_NUMBER(u1)
    CALL RANDOM_NUMBER(u2)
    theta = ASIN(SQRT(u1))
    phi = 2.0_8 * PI * u2

    px = SIN(theta)*COS(phi)
    py = SIN(theta)*SIN(phi)
    pz = COS(theta)

    p_len = SQRT(px*px + py*py + pz*pz)

    px = (target_p/p_len)*px
    py = (target_p/p_len)*py
    pz = (target_p/p_len)*pz
END SUBROUTINE randomPointTrapOptimum

SUBROUTINE cross(a, b, res)
	IMPLICIT NONE
	REAL(KIND=PREC), DIMENSION(3), INTENT(IN) :: a, b
	REAL(KIND=PREC), DIMENSION(3), INTENT(OUT) :: res
	res(1) = a(2)*b(3)-a(3)*b(2)
	res(2) = a(1)*b(3)-a(3)*b(1)
	res(3) = a(1)*b(2)-a(2)*b(1)
END SUBROUTINE cross

SUBROUTINE resetState(state)
	!Reset the 2 tracks back along their separation vectors.
	IMPLICIT NONE
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	real(kind=PREC), dimension(2,6) :: delta
	real(kind=PREC), dimension(1,6) :: rejection
	real(kind=PREC) :: lenA, lenB, lenrejAB
	
	!Calculate difference vectors
	delta(1,:) = (state(2,:)-state(1,:))!&
	!/METRIC_SCALE
	delta(2,:) = (state(3,:)-state(1,:))!&
	!/METRIC_SCALE
	
	!Calculate the rejection, or perp. part of A from B
	rejection(1,:) = delta(2,:)&
					-delta(1,:)&
						*(SUM(delta(1,:)*delta(2,:))/SUM(delta(1,:)**2))
	
	! Length is just the sum of the square of the difference vector components
	lenA = SQRT(SUM(delta(1,:)**2))
	lenrejAB = SQRT(SUM(rejection(1,:)**2))
	
	!Normalize and re-scale the difference to the initial delta
	delta(1,:) = delta(1,:) / lenA
	rejection(1,:) = rejection(1,:) / lenrejAB
	delta(1,:) = (EPSILON&
		!*METRIC_SCALE)&
		* delta(1,:))
	rejection(1,:) = (EPSILON&
		!*METRIC_SCALE)&
		* rejection(1,:))
	
	!Set the secondary trajectories to be ref+delta
	state(2,:) = state(1,:) + delta(1,:)
	state(3,:) = state(1,:) + rejection(1,:)
	
!	PRINT *, "Reset Separations:"
!	PRINT *, SQRT(SUM(((state(1,:)-state(2,:))/METRIC_SCALE)**2))
!	PRINT *, SQRT(SUM(((state(1,:)-state(3,:))/METRIC_SCALE)**2))
END SUBROUTINE resetState

SUBROUTINE createInitialTrajectories(state)
	IMPLICIT NONE
	real(kind=PREC), dimension(3,6), intent(inout) :: state
	real(kind=PREC), dimension(3) :: randomP, perp1, perp2, a, b
	real(kind=PREC) :: pNorm, pLen
	
	! Generate a random vector
!	randomP(1) = RAND(0)
!	randomP(2) = RAND(0)
!	randomP(3) = RAND(0)
	CALL RANDOM_NUMBER(randomP)
	
	!Normalize and scale it so it's similar to the reference trajectory
	pLen = SQRT(SUM(state(1,4:6)**2))
	pNorm = SQRT(SUM(randomP**2))

	randomP = randomP/pNorm
	randomP = pLen*randomP
	
	!Cross Product, generates 2 perpendicular vectors
	CALL cross(state(1,4:6), randomP, perp1)
	CALL cross(state(1,4:6), perp1, perp2)
	
	!Now that we have 2 perpindicular vectors, start with the initial p vector
	!For both trajectories. Then subtract delta from the vector and add
	!Delta in the perpindicular direction.
	
	!Now we have 2 perpindicular vectors to the trajectory in P space.
	!Want to subtract a from P and add b in perp direction
	!Constrains distance to be epsilon, size to still be P
	!N.B. Only approximate separation in metric space.
	a=EPSILON*EPSILON*4.0e-27_8*4.0e-27_8/(2.0e0_8*pNorm)
	b=EPSILON*4.0e-27_8*SQRT(1.0e0_8-EPSILON*EPSILON*4.0e-27_8*4.0e-27_8/(4.0e0_8*pNorm*pNorm))
    !a=(EPSILON**2)*(PSCALE**2)/(2.0e0_8*pNorm)
	!b=EPSILON*PSCALE*SQRT(1.0e0_8-(EPSILON**2)*(PSCALE**2)/(4.0e0_8*pNorm*pNorm))
	
	!Normalize and scale by b
	pNorm = SQRT(SUM(perp1**2))
	perp1 = perp1 / pNorm
	perp1 = perp1 * b
	
	pNorm = SQRT(SUM(perp2**2))
	perp2 = perp2 / pNorm
	perp2 = perp2 * b
	
	!Subtract off a and add in b in perp direction
	state(2,4:6) = (1.0e0_8-a/pLen)*state(1,4:6)+perp1
	state(3,4:6) = (1.0e0_8-a/pLen)*state(1,4:6)+perp2
	
	!keep real space coord.
	state(2,1:3) = state(1,1:3)
	state(3,1:3) = state(1,1:3)
!	PRINT *, "INIT Separations:"
!	PRINT *, SQRT(SUM(((state(1,:)-state(2,:))/METRIC_SCALE)**2))
!	PRINT *, SQRT(SUM(((state(1,:)-state(3,:))/METRIC_SCALE)**2))
END SUBROUTINE createInitialTrajectories

END MODULE
