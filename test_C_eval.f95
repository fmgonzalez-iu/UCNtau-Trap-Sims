#include "modules/constants.h"

PROGRAM track
!	USE mpi
	USE constants
	USE testSubroutines
	USE forcesAndPotential
	USE trackGeometry

	IMPLICIT NONE
	REAL(KIND=PREC) :: x, y, z, fx, fy, fz, totalU, energy, sympT
	REAL(KIND=PREC) :: energy_start, energy_end, maxEgain
	REAL(KIND=PREC) :: freq, height
	REAL(KIND=PREC), ALLOCATABLE :: states(:,:)
	CHARACTER(LEN=256) :: arg
	CHARACTER(LEN=256) :: fName, fName2
	CHARACTER(LEN=256) :: rankString
	INTEGER :: i, j, k
	INTEGER :: seedLen
	INTEGER, DIMENSION(32) :: rngSeed
	INTEGER :: rank, size, tag, next, from, ierr, workerIt, trajPerWorker
	INTEGER :: ntraj

!	CALL MPI_INIT(ierr)
!	CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
!	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
	!rank = 2
	IF (IARGC() .NE. 4 .AND. IARGC() .NE. 3) THEN
		PRINT *, "Error! Not enough or too many arguments!"
		PRINT *, "timestep n_traj OUTFILE (OUTFILE2)"
!		CALL MPI_FINALIZE(ierr)
		CALL EXIT(0)
	END IF
	
	CALL GETARG(1, arg)
	READ(arg,*) dt
	CALL GETARG(2, arg)
	READ(arg,*) ntraj
	CALL GETARG(3, fName)
!	READ(arg, FORMAT=A) fName
	
	WRITE (rankString, "(I0)") rank
	
	fName = TRIM(fName) // TRIM(rankString)
		
	PRINT *, fName
	OPEN(UNIT=1,FILE=fName, FORM='UNFORMATTED')
	
	IF (IARGC() .EQ. 4) THEN
		CALL GETARG(4, fName2)
		WRITE (rankString, "(I0)") rank
		fName2 = TRIM(fName2) // TRIM(rankString)		
		PRINT *, fName2
		OPEN(UNIT=2, FILE=fName2, FORM='UNFORMATTED')
	END IF
	
	trajPerWorker = ntraj/size
	
!	minU = -2.4283243003838247e-26_8
	minU = -2.390245661413933e-26_8
	minU = -2.390352484438862e-26_8 !For lambda = 0.05114, brem=1.35
	
	ALLOCATE(states(ntraj,6))
		
	PI=4.0e0_8*ATAN(1.0e0_8)
	a(1)=.5153528374311229364e0_8
	a(2)=-.085782019412973646e0_8
	a(3)=.4415830236164665242e0_8
	a(4)=.1288461583653841854e0_8
	b(1)=.1344961992774310892e0_8
	b(2)=-.2248198030794208058e0_8
	b(3)=.7563200005156682911e0_8
	b(4)=.3340036032863214255e0_8
!	PRINT *, SUM(a), SUM(b)
!	a(1)=(1.0_8/6.0_8)&
!		*(2.0_8+2.0_8**(1.0_8/3.0_8)&
!			+2.0_8**(-1.0_8/3.0_8))
!	a(2)=(1.0_8/6.0_8)&
!		*(1.0_8-2.0_8**(1.0_8/3.0_8)&
!			-2.0_8**(-1.0_8/3.0_8))
!	a(3)=a(2)
!	a(4)=a(1)
!	b(1)=0.0_8
!	b(2)=1.0_8/(2.0_8-2.0_8**(1.0_8/3.0_8))
!	b(3)=1.0_8/(1.0_8-2.0_8**(2.0_8/3.0_8))
!	b(4)=b(2)

	CALL RANDOM_SEED(size=seedLen)
	IF (seedLen > 32) THEN
		PRINT *, "Error! The requested length of seed is too long"
		CALL EXIT(0)
	END IF
!	PRINT *, "called a seed!"
	!I'm not going to care about proper types since it's just for seed values
	rngSeed(1) = 4434
	DO i=2,seedLen,1
		rngSeed(i) = MOD((48271*rngSeed(i-1)), 2147483647)
	END DO
	CALL RANDOM_SEED(put=rngSeed(1:seedLen))
	
!	DO i=1,100,1
!		CALL zOffDipCalc(i*(250.0_8/100.0_8), z)
!		PRINT *, z, i*(250.0_8/100.0_8)
!	END DO
			
	DO i=1,ntraj,1
		!CALL randomPointTrap(states(i,1), states(i,2), states(i,3), states(i,4), states(i,5), states(i,6))
		CALL randomPointTrapOptimum(states(i,1), states(i,2), states(i,3),&
		 states(i,4), states(i,5), states(i,6))
	END DO
	!PRINT *, states(1,:)
	!PRINT *, "called a random point trap"
	
!	PRINT *, trajPerWorker*rank+1, trajPerWorker*(rank+1)
	DO i=trajPerWorker*rank+1,trajPerWorker*(rank+1),1
!		sympT = 0.0_8
!       CALL trackAndPrint(states(i,:))
!		CALL trackEnergyGain(states(i,:), energy_start, energy_end, sympT, 30.0_8)
!		CALL trackEnergyGain(states(i,:), energy_start, energy_end)
!		PRINT *, rank, i, energy_start, energy_end, (energy_end - energy_start)/energy_start
!		PRINT *, rank, i, (energy_end - energy_start)/energy_start
		IF (IARGC() .EQ. 3) THEN
!			PRINT *, "Starting run without block!"
			CALL trackDaggerHitTime(states(i, :))
		ELSE IF (IARGC() .EQ. 4) THEN
!			PRINT *, "Starting run with Aluminum block!"
			CALL trackDaggerAndBlock(states(i, :))
		ELSE
!			PRINT *, "Something wrong with iargc!"
!			CALL MPI_FINALIZE(ierr)
			CALL EXIT(0)
		END IF 
!		CALL trackDaggerHitTimeFixedEff(states(i, :))
	END DO
    
!    CALL trackDaggerHitTime(states(ntraj, :))
    
!    CALL trackAndPrint(states(ntraj,:))

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
!	CALL trackDaggerHitTime(states(ntraj, :))
	
	
!	CALL MPI_FINALIZE(ierr)
	
END PROGRAM track
