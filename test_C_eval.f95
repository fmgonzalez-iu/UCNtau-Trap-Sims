#include "modules/constants.h"

PROGRAM track
	! COMPILER ISSUE: My personal laptop can't do mpi multithreading. All
	!				  MPI functions/subroutines are commented out for now.
	! 				  Additionally, there are some PRINT warnings that should
	!				  be commented for BigRed2, or they'll write default values.
	! FILE CHANGES: The module convertTrace has the heating files hard-coded in.
	!				This means that for running on a different machine, need to
	!				change which file path we're using. (Eventually will have 
	!				to do this with magnetic fields too.)
	! STILL TO DO: Lyapunov exponent and other cleaner studies don't compile.
	!			   Again, they're commented out for for now, but will need
	!			   to be implemented in the future, possibly with more flags.
	
	!USE mpi
	USE constants
	USE testSubroutines
	USE convertTrace
	USE forcesAndPotential
	USE trackGeometry
	!USE lyapunov

	!-------------------------------------------------------------------
	!------------------ Initialize Variables ---------------------------
	!-------------------------------------------------------------------
	IMPLICIT NONE
!	REAL(KIND=PREC) :: x, y, z, fx, fy, fz, totalU, energy, sympT
!	REAL(KIND=PREC) :: energy_start, energy_end, maxEgain
!	REAL(KIND=PREC) :: freq, height, holdTime, bScale
!	REAL(KIND=PREC), DIMENSION(:), ALLOCATABLE :: trX(:), trY(:), trZ(:)
!	REAL(KIND=PREC), ALLOCATABLE :: states(:,:)
!	REAL(KIND=PREC) :: res_lyap(17)
	
!	CHARACTER(LEN=256) :: arg, fName, fName2, rankString
!	INTEGER :: i, j, k, seedLen, seedOff, lengthTr
!	INTEGER, DIMENSION(32) :: rngSeed
!	INTEGER :: rankMPI, sizeMPI, tag, next, from, ierr, workerIt, trajPerWorker
!	INTEGER :: ntraj, nDips, pse

	! New variable lists because I'm trying to optimize -- others are for other code
	REAL(KIND=PREC) :: holdTime, bScale, minE, maxE, heatF
	REAL(KIND=PREC) :: res_lyap(17)
	REAL(KIND=PREC) :: res_clean(10)
	REAL(KIND=PREC), DIMENSION(:), ALLOCATABLE :: trX(:), trY(:), trZ(:)
	REAL(KIND=PREC), ALLOCATABLE :: states(:,:)
	LOGICAL :: nDecay
	INTEGER :: ntraj, nDips, pse, dist, &
			   rankMPI, sizeMPI, ierr, trajPerWorker, &
			   seedOff, seedLen, &
			   lengthTr, &
			   i
	INTEGER, DIMENSION(32) :: rngSeed
	CHARACTER(LEN=256) :: arg, fName, fName2, rankString
	
	rankMPI = 1
	sizeMPI = 1
	
!	CALL MPI_INIT(ierr)
!	CALL MPI_COMM_RANK(MPI_COMM_WORLD, rankMPI, ierr)
!	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, sizeMPI, ierr)

	!-------------------------------------------------------------------
	!------------------ Load Command Line Args -------------------------
	!-------------------------------------------------------------------
	IF (IARGC() .GT. 12 .OR. IARGC() .LT. 3) THEN
		PRINT *, "ERROR: Not enough or too many arguments!"
		PRINT *, "REQUIRED: timestep n_traj OUTFILE"
		PRINT *, "OPTIONAL: (holdTime nDips PSE) (minE maxE dist heatF) (bScale OUTFILE2)"
		PRINT *, "Also note that this order has changed recently!"
!		CALL MPI_FINALIZE(ierr)
		CALL EXIT(0)
	END IF
	
	! I'm hardcoding this in because I'm too lazy to modify command line
	nDecay = .TRUE.
	
	CALL GETARG(1, arg)
	READ(arg,*) dt
	
	CALL GETARG(2, arg)
	READ(arg,*) ntraj
	
	CALL GETARG(3, fName)
	WRITE (rankString, "(I0)") rankMPI
	fName = TRIM(fName) // TRIM(rankString)
		
	IF (IARGC() .GE. 4) THEN
		CALL GETARG(4, arg)
		READ(arg,*) holdTime
	ELSE
		PRINT *, "DEFAULT: holdTime defaulting to 20.0s"
		holdTime = 20.0_8
	END IF
	
	IF (IARGC() .GE. 5) THEN
		CALL GETARG(5, arg)
		READ(arg,*) nDips
		nDips = nDips + 1 !First "dip" is actually second element of vector
	ELSE 
		PRINT *, "DEFAULT: nDips defaulting to 3"
		nDips = 4
	END IF
	
	IF (IARGC() .GE. 6) THEN
		CALL GETARG(6, arg)
		READ(arg,*) pse
		IF (pse .LT. 1 .OR. pse .GT. 3) THEN
			!	CASE 1: Normal running
			!	CASE 2: Deep cleaning (clean at 25 cm) [for phase space evolution]
			!	CASE 3: Height dependent time constant data
			PRINT *, "WARNING: chosen pse not available! Defaulting to normal running"
			pse = 1
		END IF
	ELSE
		PRINT *, "DEFAULT: pse defaulting to 1 (normal running)"
		pse = 1
	END IF
	
	IF (IARGC() .GE. 7) THEN
		CALL GETARG(7, arg)
		READ(arg,*) minE
	ELSE
		PRINT *, "DEFAULT: minE defaulting to 1 cm"
		minE = 0.01_8
	END IF
	
	IF (IARGC() .GE. 8) THEN
		CALL GETARG(8, arg)
		READ(arg,*) maxE
	ELSE
		PRINT *, "DEFAULT: maxE defaulting to 45 cm"
		maxE = 0.45_8
	END IF
	
	IF (IARGC() .GE. 9) THEN
		CALL GETARG(9, arg)
		READ(arg,*) dist
		IF (dist .LT. 0 .OR. dist .GT. 3) THEN
			PRINT *, "WARNING: Unknown spectrum distribution! Defaulting to flat in theta!"
			dist = 3
		END IF
	ELSE
		PRINT *, "DEFAULT: dist defaulting to 'flat in theta'!"
		dist = 3
	END IF 
	
	IF (IARGC() .GE. 10) THEN
		CALL GETARG(10, arg)
		READ(arg,*) heatF
	ELSE
		PRINT *, "DEFAULT: heatF defaulting to 0 (no heating)!"
		heatF = 0.0_8
	END IF
	
	IF (IARGC() .GE. 11) THEN	
		CALL GETARG(11, arg)
		READ(arg,*) bScale
		IF (bScale .LT. 0.0_8) THEN
			PRINT *, "WARNING: block scaling cannot be negative! Removing block!"
			bScale = 0.0_8
		END IF
	ELSE
		PRINT *, "DEFAULT: bScale defaulting to 0.0 (no block)"
		bScale = 0.0_8
	END IF
	
	IF (IARGC() .GE. 12) THEN
		CALL GETARG(12, fName2)
		WRITE (rankString, "(I0)") rankMPI
		fName2 = TRIM(fName2) // TRIM(rankString)		
	END IF
		
	! Write the filenames, and open their files
	PRINT *, fName
	OPEN(UNIT=1,FILE=fName, FORM='UNFORMATTED')
	! There's a flag here for the block.
	IF (IARGC() .GE. 12) THEN
		PRINT *, fName2
		OPEN(UNIT=2, FILE=fName2, FORM='UNFORMATTED')	
	END IF
	
	!-------------------------------------------------------------------
	!---------------- Start creating trajectories ----------------------
	!-------------------------------------------------------------------
	! Set values of global variables -- don't need to declare these (or dt)
	PI=4.0e0_8*ATAN(1.0e0_8)
	a(1)=.5153528374311229364e0_8
	a(2)=-.085782019412973646e0_8
	a(3)=.4415830236164665242e0_8
	a(4)=.1288461583653841854e0_8
	b(1)=.1344961992774310892e0_8
	b(2)=-.2248198030794208058e0_8
	b(3)=.7563200005156682911e0_8
	b(4)=.3340036032863214255e0_8
	! I don't know what minU is?
	minU = -2.390352484438862e-26_8 !For lambda = 0.05114, brem=1.35
	
	! Old but OK random number generator. Might want to improve later?
	CALL RANDOM_SEED(SIZE=seedLen)
	seedOff = 0
	IF (seedLen + seedOff > 32) THEN
		PRINT *, "Error! The requested length of seed is too long"
		CALL EXIT(0)
	END IF
	rngSeed(1) = 4434
	DO i=2,seedLen,1
		rngSeed(i) = MOD((48271*rngSeed(i-1)), 2147483647)
	END DO
	CALL RANDOM_SEED(put=rngSeed(1+seedOff:seedLen+seedOff))

	! Figure out trajectories
	trajPerWorker = ntraj/sizeMPI	
	ALLOCATE(states(ntraj,6))
	DO i=1,ntraj,1
		CALL randomPointTrap(states(i,1), states(i,2), states(i,3), &
								states(i,4), states(i,5), states(i,6), &
								minE, maxE, dist)
	END DO
	
	! For Lyapunov:
!	ALLOCATE(states(ntraj,3,6))
!	DO i=1,ntraj,1
!		CALL randomPointTrap(states(i,1), states(i,2), states(i,3), &
!								states(i,4), states(i,5), states(i,6), &
!								minE, maxE, dist)
!		CALL createInitialTrajectories(states(i,:,:))
!	END DO

	! Load heating (If heating is on!)
	IF (heatF .NE. 0.0_8) THEN
		CALL loadTrace(trX,trY,trZ,lengthTr,heatF)
	ELSE
		! If no heating, allocate our traces as 1D to free memory
		ALLOCATE(trX(1))
		ALLOCATE(trY(1))
		ALLOCATE(trZ(1))
	END IF
	
	!-------------------------------------------------------------------
	!------------------- Begin running tracking code -------------------
	!-------------------------------------------------------------------
	
	DO i=1,ntraj,1
	!DO i=trajPerWorker*rankMPI+1,trajPerWorker*(rankMPI+1),1
!		PRINT *, "Beginning full runs!"
		! It's faster to not carry traces (a big matrix) in the calc if there's no heating!
		IF (heatF .EQ. 0.0_8) THEN
!			PRINT *, "Heating is presently OFF"
			CALL trackDaggerFull(states(i, :),holdTime, bScale, nDips, pse, nDecay)
		ELSE
!			PRINT *, "Heating is presently ON"
			CALL trackDaggerFull(states(i, :),holdTime, bScale, nDips, pse, &
								nDecay, trX, trY, trZ, lengthTr)
		END IF 
	END DO
	
	!-------------------------------------------------------------------
	!------------------- Optional, unimplemented code ------------------
	!-------------------------------------------------------------------

	! I don't know what minU actually does, but these look like fancy numbers
	! so I don't want to delete them... 
!	minU = -2.4283243003838247e-26_8
!	minU = -2.390245661413933e-26_8 !For lambda = 0.0508, brem=1.4


	

	
!    CALL trackAndPrint(states(ntraj,1,:))

	!DO i=1,ntraj,1
	!DO i=trajPerWorker*rankMPI+1,trajPerWorker*(rankMPI+1),1


		
!		PRINT *, "Beginning other types of runs! (unimplemented)"
!		CALL calcLyapunov(states(i,:,:), res_lyap)
!		WRITE(1), res_lyap
!		CALL calcCleanTime(states(i,:,:), res_clean)
!		WRITE(1), res_clean
!		PRINT *, i
!	END DO
	
!	DO i=1,81,1 !Freq
!		DO j=0,39,1 !Height
!			maxEgain = 0.0_8
!			DO k=0,40,1 !Phase
!				freq = i*2000_8/(80_8)
!				height = -1.4 + 0.4 * (j/40.0)
!				sympT = 0.0_8 + (1.0_8/freq) * k / 20.0
!				CALL testEnergyGain(freq, height, sympT, energy_start, energy_end)
!				IF ((energy_end - energy_start) > maxEgain) THEN
!					maxEgain = (energy_end - energy_start)
!				END IF
!!				PRINT *, freq, height, j / 20.0 * 2.0 * PI, energy_start*JTONEV, (energy_end - energy_start)*JTONEV
!			END DO
!			PRINT *, freq, height, j / 20.0 * 2.0 * PI, energy_start*JTONEV, maxEgain*JTONEV
!!			PRINT *, i, j, maxEgain*JTONEV
!		END DO
!	END DO

!	CALL trackAndPrint(states(ntraj, :), 0.0_8)

	!-------------------------------------------------------------------
	!--------------------- END OF PROGRAM ------------------------------
	!-------------------------------------------------------------------
	
!	CALL MPI_FINALIZE(ierr)
	
END PROGRAM track
