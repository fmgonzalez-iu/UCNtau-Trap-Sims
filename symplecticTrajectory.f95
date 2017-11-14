#include "modules/constants.h"

!SUBROUTINE totalForce(x, y, z, fx, fy, fz, loops)
!	!Calculate force from B-field and add in gravity
!	IMPLICIT NONE
!	real(kind=PREC), intent(in) :: x, y, z
!	real(kind=PREC), intent(out) :: fx, fy, fz
!	real(kind=PREC), dimension(3,NWIRES+1), intent(in) :: loops
!	real(kind=PREC) :: k
!	
!	k=MASS_N
!	fy=0
!	fz=0
!	fx=-k*x
!END SUBROUTINE totalForce

PROGRAM track
	USE mpi
	USE constants
	USE ellipticInt
	USE forcesAndPotential
	USE wires
	USE lyapunov
	USE symplecticInt
	USE testSubroutines
	USE trackGeometry
	USE loopSetup
	IMPLICIT NONE
	integer :: n
	character(len=64) :: arg
	real(kind=PREC) :: x, y, z, px, py, pz
	real(kind=PREC) :: t, ke, v, e, eStart, maxDifference, minDifference
	real(kind=PREC) :: diff, movingAvg
	real(kind=PREC) :: kTest, eTest
	!real(kind=PREC), dimension(1000,3,6) :: states
	real(kind=PREC), allocatable :: states(:,:,:)
	real(kind=PREC), dimension(10) :: results
	integer, dimension(32) :: rngSeed
	integer(kind=8) :: i
	integer(kind=8) :: numStep
	integer :: numPoints
	integer :: xIt, zIt
	integer :: seedLen
	integer :: ntraj
	integer :: rank, size, tag, next, from, ierr, workerIt, trajPerWorker, trajIndex

	CALL MPI_INIT(ierr)
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
	
	IF (IARGC() .NE. 10) THEN
		PRINT *, "Error! Not enough or too many arguments!"
		PRINT *, "nwires major_r minor_r xmax dt ntraj liptime nsep pertStrength trajectory_index"
		CALL MPI_FINALIZE(ierr)
		CALL EXIT(0)
	END IF
	
	CALL GETARG(1, arg)
	READ(arg,*) nwires
	!PRINT *, "nwires:",nwires
	CALL GETARG(2, arg)
	READ(arg,*) torus_major_r
	!PRINT *, "major_r:",torus_major_r
	CALL GETARG(3, arg)
	READ(arg,*) torus_minor_r
	!PRINT *, "minor_r:",torus_minor_r
!	CALL GETARG(4, arg)
!	READ(arg,*) xOff
!	!PRINT *,
!	CALL GETARG(5, arg)
!	READ(arg,*) yOff
!	!PRINT *,
	CALL GETARG(4, arg)
	!PRINT *, "xmax:",xmax
	READ(arg,*) xmax
	CALL GETARG(5, arg)
	READ(arg,*) dt
	CALL GETARG(6, arg)
	READ(arg,*) ntraj
	CALL GETARG(7, arg)
	READ(arg,*) liptime
	CALL GETARG(8, arg)
	READ(arg,*) nsep
	CALL GETARG(9, arg)
	READ(arg,*) pertStrength
	CALL GETARG(10, arg)
	READ(arg,*) trajIndex
	
	!Inside constants subroutine
	PI=4.0e0_8*ATAN(1.0e0_8)
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
	a(1)=.5153528374311229364e0_8
	a(2)=-.085782019412973646e0_8
	a(3)=.4415830236164665242e0_8
	a(4)=.1288461583653841854e0_8
	b(1)=.1344961992774310892e0_8
	b(2)=-.2248198030794208058e0_8
	b(3)=.7563200005156682911e0_8
	b(4)=.3340036032863214255e0_8
	
	trajPerWorker = ntraj/size
	
	allocate(loops(5,nwires+1))
	CALL setUpLoops()
	
	nPert = 0
	pertStrength = 0.0_8

!	CALL setUpLoopsExplicit(N8WIRESUI/PI, N8WIRESUMIN*MASS_N,&
!		N8WIRESULIP*MASS_N)
!	CALL createFinitePertWires(0.5_8*N8WIRESUI/(2.0_8*PI))
!	pertStrength = 0.5_8*N8WIRESUI/(2.0_8*PI)
	
!	PRINT *, SUM(a)
!	PRINT *, SUM(b)
	
	!Since, in principle, the # generator may need multiple seeds
	!Use a quick & easy RNG (park & miller).
	CALL RANDOM_SEED(size=seedLen)
	IF (seedLen > 32) THEN
		PRINT *, "Error! The requested length of seed is too long"
		CALL EXIT(0)
	END IF
	!I'm not going to care about proper types since it's just for seed values
	rngSeed(1) = 4434
	DO i=2,seedLen,1
		rngSeed(i) = MOD((48271*rngSeed(i-1)), 2147483647)
	END DO
	CALL RANDOM_SEED(put=rngSeed(1:seedLen))

	allocate(states(ntraj,3,6))

	DO i=1,ntraj,1
		CALL randomPoint(states(i,1,1), states(i,1,2), states(i,1,3),&
			states(i,1,4), states(i,1,5), states(i,1,6), uTrapMax)
		!PRINT *, states(i,1,1)
		!CALL createInitialTrajectories(states(i,:,:))
	END DO

!	DO i=trajPerWorker*rank+1,trajPerWorker*(rank+1),1
!		!WRITE(*,"(I5,I5)") rank, i
!		CALL calcLyapunov(states(i,:,:), results)
!		WRITE(*,"(2I6,10ES28.19E3)") rank, i, results
!		!PRINT *, i
!	END DO

	IF(rank .ne. 0) THEN
		CALL MPI_FINALIZE(ierr)
		CALL EXIT(0)
	END IF
	
	!CALL calcx0Mesh()
	CALL trackAndPrint(states(trajIndex,:,:))
	!CALL testVecExp()
	
	
	
	!CALL trackLyapunov(states(trajIndex,:,:))
	!CALL trackSumLyapunov(states(trajIndex,:,:))
	!CALL calcx0Mesh()
	!CALL lowFieldSearch()
	
	CALL MPI_FINALIZE(ierr)
	
END PROGRAM track