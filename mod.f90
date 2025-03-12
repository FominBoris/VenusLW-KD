! MODULE MOD_IN_GAS_CLOUD ! Information from atmospheric profile
module mod_in_gas_cloud   ! MOD_IN_GAS_CLOUD
       	CHARACTER*36 PASS_QUADR
! Attention NANGLE2 <= N_LEGENDRE (line 173)
	CHARACTER ATM_PATH*30
	CHARACTER*80 TITLE
   INTEGER*4 :: NGAS,JMAX ,ICL_MIN,ICL_MAX, ISO2, IH2O, ICO2
     REAL*4, ALLOCATABLE :: Z(:),T1(:),RO1(:,:),RABMA(:)
  REAL*8, ALLOCATABLE ::   P1(:),P_LN(:)
  PARAMETER (NANGLE=30,NANGLE2=NANGLE*2) ! 20 Quadrature !!!!! ######
  PARAMETER (NTH=1, NCOMP = 3)  ! for KD  H20,CO2 and SO2 only NTH=20481 H=...
	CHARACTER*5  MOLECULE(NCOMP)
   REAL*8 COSINE(NANGLE),WEIGHT(NANGLE)
   DIMENSION RK(NTH)
!   END MODULE MOD_IN_GAS_CLOUD
   end module mod_in_gas_cloud
   SUBROUTINE ATM_PROF_READING
   USE MOD_IN_GAS_CLOUD
   CHARACTER*1 A_ZO
   CHARACTER ATPA*60  !###
   DIMENSION Z_CL(99) ! MAXimal namber of layers in cloud
    OPEN(110,FILE='ATM_AEROSOL.MODEL') ! input File - AEROSOL (CLOUD) model
    READ(110,*) NLAYER  ! Number of Layers with different Phase Functions
NLAY=NLAYER ; N_CL=0
  IF(NLAY>0)THEN
! ###   NLAY=ABS(NLAYER)
   READ(110,*)WN5 ! here not essential
!--- 2023 ----!
   IPER=(NLAYER-0.1)/10+1
     DO I_R=1,IPER ; READ(110,555)E55 ; END DO ! here E55 not essential
555 FORMAT(10E12.5)
!-------------!
	  READ(110,*)N_TEC,ZDOW,ZU      ! Boundaries in km
    N_CL=2 ; Z_CL(1)=ZDOW ; Z_CL(2)=ZU
 IF(NLAY > 1)THEN
      READ(110,*)NFRAC     ! Number of fractions
      READ(110,*)WFRA ! here not essential       ! their Weights
        DO J=1,NFRAC
        READ(110,*)A_ZO ! here not essential ! path to aerosol model in database
        END DO
        DO I=2,NLAY
	  READ(110,*)N_TEC,ZDOW,ZU      ! Boundaries in km
      READ(110,*)NFRAC  ! Number of fractions
      READ(110,*)WFRA ! here not essential       ! their Weights
        DO J=1,NFRAC
        READ(110,*)A_ZO ! here not essential ! path to aerosol model in database
        END DO
   IF(ZDOW /= Z_CL(N_CL)) THEN
           N_CL=N_CL+1 ; Z_CL(N_CL)=ZDOW
   END IF
           N_CL=N_CL+1 ; Z_CL(N_CL)=ZU
       END DO
 END IF
END IF
    CLOSE (110)
 IF(NLAY>0)THEN
       IF(N_CL.LT.2.OR.N_CL.GT.99) THEN
       WRITE(*,*)' CHECK CLOUD/AEROSOL MODEL:',NLAY,ZDOW,ZU
	   STOP
       END IF
 END IF
! -------
       OPEN(66,FILE='k_coef.chk')
77 FORMAT(A20)
OPEN(55,FILE='./Atmospheres/'//ATM_PATH)
!*
	READ(55,'(A)')TITLE
	WRITE(66,'(A)')TITLE
	READ(55,*)NGAS,JMAX
JMAXC=JMAX+N_CL
		DO I=1,NGAS
		READ(55,'(A)')MOLECULE(I)
		END DO
        ALLOCATE (Z(JMAXC),P1(JMAXC),T1(JMAXC),P_LN(JMAXC))
        ALLOCATE (RO1(NGAS,JMAXC),STAT=IERR)
ALLOCATE (RABMA(JMAXC),STAT=IERR)
         IF(IERR/=0)THEN
          WRITE(*,*)' Allocaion is wrong !!!'
          STOP
         END IF

	DO J=1,JMAX
	READ(55,*)Z(J),P1(J),T1(J),( RO1(I,J), I = 1, NGAS )
	END DO
		CLOSE(55)
! Inserting CLOUD into atmospheric model
! *************************************************
IF(NLAY>0)THEN
             GAS_CLOUD=0.00001 ! [KM] (If /Z_CL-Z/ < GAS_CLOUD levels are the same)
			 ICL_MAX=-1 ; ICL_MIN=100000
        DO ICL=1,N_CL
		ZC=Z_CL(ICL)
		 DO JGA=2,JMAX
		 JGAM=JGA-1
		 IF(Z(JGAM)<=ZC.AND.ZC<Z(JGA)) EXIT
         END DO
! Z(J-1)=ZC
          IF(ABS(Z(JGAM)-ZC)<GAS_CLOUD) THEN
		  Z(JGAM)=Z_CL(ICL)
		   IF(JGA<ICL_MIN)ICL_MIN=JGAM
           IF(JGA>ICL_MAX)ICL_MAX=JGAM
		  ELSE
! Z(J)=ZC
           IF(ABS(Z(JGA)-ZC)<GAS_CLOUD) THEN
		   Z(JGA)=Z_CL(ICL)
		   IF(JGA<ICL_MIN)ICL_MIN=JGA
           IF(JGA>ICL_MAX)ICL_MAX=JGA
		   ELSE
! Z(J-1)< ZC < Z(J) - RAZDVIDGKA
     		   DO JJ=JMAX,JGA,-1
           Z(JJ+1)=Z(JJ)
           P1(JJ+1)=P1(JJ)
		   T1(JJ+1)=T1(JJ)
		       DO III=1,NGAS
		   RO1(III,JJ+1)=RO1(III,JJ)
		       END DO
		   END DO
		   Z(JGA)=ZC
           C1=(Z(JGA+1)-ZC)/(Z(JGA+1)-Z(JGAM))
		   C2=1.-C1
		   P1(JGA)=P1(JGAM)*C1+P1(JGA+1)*C2
		   T1(JGA)=T1(JGAM)*C1+T1(JGA+1)*C2
		       DO III=1,NGAS
		   RO1(III,JGA)=RO1(III,JGAM)*C1+RO1(III,JGA+1)*C2
		       END DO
			   continue
              JMAX=JMAX+1
		   IF(JGA<ICL_MIN)ICL_MIN=JGA
           IF(JGA>ICL_MAX)ICL_MAX=JGA
		   END IF
          END IF
		END DO
 !### WRITE(66,*)'CLOUD between: ',ICL_MIN,ICL_MAX,Z(ICL_MIN),Z(ICL_MAX)
 END IF

	WRITE(*,*)'The number of gases (NGAS) : ',NGAS
	WRITE(*,*)'The number of levels (JMAX) : ',JMAX
WRITE(66,*)NGAS,JMAX
		DO I=1,NGAS
WRITE(66,566)MOLECULE(I)
566 FORMAT(A5)
		END DO
	DO J=1,JMAX
	P1(J)=P1(J)*1013.
    P_LN(J)=DLOG(P1(J))

!###	WRITE(66,*)J,Z(J),P1(J),T1(J), (RO1(I,J), I = 1, NGAS )
	WRITE(66,676)Z(J),P1(J),T1(J),( RO1(I,J), I = 1, NGAS )
676 FORMAT(F10.3,E12.5,F8.2,10E12.4)
	END DO
                CLOSE(66)
 !  ************************************** For KD ***************************** !
 ISO2=0
 IH2O=0
 ICO2=0
                  DO I=1,NGAS
                  IF(MOLECULE(I).EQ.'H2O')IH2O=I
                  IF(MOLECULE(I).EQ.'CO2') ICO2=I
                  IF(MOLECULE(I).EQ.'SO2') ISO2=I
                   END DO
  END SUBROUTINE ATM_PROF_READING

 ! ----------------------------------- !
 MODULE M_C
 REAL*4 X0_,Y0_,Z0_,P1_0,P2_0,P3_0,P4_0 !Coordinates and Stokes parameters of initial ray.
 REAL*8 A0_,B0_,C0_ &  !  COS of initial ray.
        ,X_RAND        ! For the random data
 INTEGER*4 N_POINTS ! The real number of scattering points in the trajectory.
 PARAMETER (IMAX_MC=998,IMAX_M_=IMAX_MC+2) !The maximal number of scattering points in the trajectory (cut off).
 REAL*4 X_(IMAX_M_),Y_(IMAX_M_),Z_(IMAX_M_),COS1_(IMAX_M_),COS2_(IMAX_M_),COS3_(IMAX_M_) &
 ,P1_(IMAX_M_),P2_(IMAX_M_),P3_(IMAX_M_),P4_(IMAX_M_) ! Trajectory points.
  ALLOCATABLE LMINV(:),CL0V(:),CL1V(:) ! for LW
  END MODULE M_C

   MODULE A_MOD ! 27 Nov.,2003-24 August,2004 (+ Creating database of Legendre series)
   REAL*8 Pi,ST_ANG
   PARAMETER  (IHOL=80,Pi=3.14159265359 &
   ,N_LEGENDRE=201 &  ! ### the LAST NUMBER ( >= NANGLE2)  of the term in the Legendre series
   ,IM4=1000 & ! DIMENSION of the Arrays for ANGLE scattering simulation.
   ,I_FCTR=4) ! I_FCTR=1 for other FORTRAN.
  ALLOCATABLE :: E55(:),ZDOWN(:),ZUP(:) &! Extinction and boundaries in each layer
 ,PH_F_I(:),ST_ANG(:),WFRAC(:,:) & ! Array over STandard ANGle grid, weights of the fractions.
 , PROBFUN(:,:),SCATT(:),ABSS(:),EXTT(:),SERLEG(:,:)
  INTEGER*4 NLAY, IMEN_PHF, IDIF
 ! Number of Layers with different Phase Functions, DIM for STANDARD grid, Number of Aerosols
  REAL*8 WN55 ! Given wavenumber cm^-1 (usually 18181.818... = 0.55 mkm)
  END MODULE A_MOD
  MODULE INITIALFORKD ! Information from atmospheric profile
      IMPLICIT REAL*8 (A-H,O-Z)
       PARAMETER (NCOMP = 3)
	CHARACTER*23 LINE_PATH, ATM_PATH
	CHARACTER*80 TITLE
	CHARACTER MOLECULE(NCOMP)*5,MLC*7
  INTEGER*4 NGAS,JMAX,ISO2,IH2O,ICO2
    REAL*8, ALLOCATABLE :: Z(:),P(:),T(:),RO(:,:),P_LN(:)
REAL*8,ALLOCATABLE :: TAUMA(:),D_TAU(:),FLUXUP(:),FLUXDO(:),QT(:),PL(:)
  CONTAINS
   SUBROUTINE ATM_PROF_READING___
      IMPLICIT REAL*8 (A-H,O-Z)
        OPEN(66,FILE='k_coef.chk')
	OPEN(99,FILE='k_coef.in')
	READ(99,'(A)')ATM_PATH
		WRITE(66,*)' The atmospheric conditions from file = ',ATM_PATH
	CLOSE(99)
	OPEN(55,FILE='./Atmospheres/'//ATM_PATH)
!*
		READ(55,5580)TITLE
5580 FORMAT(A80)
	WRITE(*,*)TITLE
pause 7
	WRITE(66,*)TITLE
	READ(55,*)NGAS,JMAX
	WRITE(*,*)'The number of gases (NGAS) : ',NGAS
	WRITE(*,*)'The number of levels (JMAX) : ',JMAX
	WRITE(66,*)'The number of gases (NGAS) : ',NGAS
	WRITE(66,*)'The number of levels (JMAX) : ',JMAX
		DO I=1,NGAS
		READ(55,'(A)')MOLECULE(I)
		END DO
	WRITE(66,*)'Account atmosphere gases		: ',	&
				(MOLECULE(I),I=1,NGAS)

        ALLOCATE (Z(JMAX),P(JMAX),T(JMAX),TAUMA(JMAX),D_TAU(JMAX) &
		,P_LN(JMAX),FLUXUP(JMAX),FLUXDO(JMAX),QT(JMAX),PL(JMAX))
        ALLOCATE (RO(NCOMP,JMAX),STAT=IERR)
         IF(IERR/=0)THEN
          WRITE(*,*)' Allocaion is wrong !!!'
          STOP
         END IF
open(155,file='z-p-T-pln')
	DO J=1,JMAX
	READ(55,*)Z(J),P(J),T(J),(RO(I,J), I = 1, NGAS )
	WRITE(66,*)Z(J),P(J),T(J),(RO(I,J), I = 1, NGAS )
	P(J)=P(J)*1013.
    P_LN(J)=DLOG(P(J))
!*	WRITE(*,*)Z(J),P(J),T(J),( RO(I,J), I = 1, NGAS )
write(155,*)Z(J),P(J),T(J),P_LN(J)
	END DO
	stop
		CLOSE(55)
                CLOSE(66)

 ISO2=0
 IH2O=0
 ICO2=0
                  DO I=1,NGAS
                  IF(MOLECULE(I).EQ.'H2O')IH2O=I
                  IF(MOLECULE(I).EQ.'CO2') ICO2=I
                  IF(MOLECULE(I).EQ.'SO2') ISO2=I
                   END DO

! ----------------------------------------------- !
  END SUBROUTINE ATM_PROF_READING___
  END MODULE INITIALFORKD


