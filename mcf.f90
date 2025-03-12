! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
! Principal Features (3-version):
! 1. No random simulations for wavenumber points (warying weight is used)
! 2. Homogeneouzly distributed Z0_ points  (warying weight is used)
! For DELTA-test see Line 362 (Use commentaries to swich on/off test,
! and  DELEETE files in ./M-C_OPT_BASE.)
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
 ! ***###   Lines for replacing RAND() by CALL RANDOM_NUMBER()
    SUBROUTINE QUAD_DEF ! Quadrature inserting
  USE MOD_IN_GAS_CLOUD
   PASS_QUADR='./LEGENDRE_SERIES/G_L_QUADRATURE.  '
   IF(NANGLE < 10)THEN
  WRITE(PASS_QUADR(34:34),'(I1)')NANGLE
  ELSE
      IF(NANGLE < 100) THEN
  WRITE(PASS_QUADR(34:35),'(I2)')NANGLE
      ELSE
  WRITE(PASS_QUADR(34:36),'(I3)')NANGLE
	  END IF
  END IF
		OPEN(500,FILE=PASS_QUADR)
 DO IA=1,NANGLE ; READ(500,*)COSINE(IA),WEIGHT(IA)
 END DO
        CLOSE(500)
	END SUBROUTINE QUAD_DEF

! ----------------------------------------------- !
   SUBROUTINE FILES_CREATING   ! 10 Nov.,2003
   USE MOD_IN_GAS_CLOUD
   USE  A_MOD
 ! Creating and reading FILES with Atmospheric Optical Properties for Monte-Carlo !
! ***** Legendre series ******** ! 24 August, 2004.
    REAL*8 COP,ANG,C,A_LEGENDRE
  CHARACTER*80 AZL,NAME_W,NAME_W2,NAME_W3
  CHARACTER*2 HVOST(99)
  CHARACTER*80 A_ZOL ! Aerosol Name-Path
  ALLOCATABLE  COP(:),PHF(:)  ! Aerosol Angle grid and Phase function in permanent DATABASE.
  ALLOCATABLE A_ZOL(:,:)
  DIMENSION A_LEGENDRE(0:N_LEGENDRE),A4(0:N_LEGENDRE)
! 1 --------------------------------------- !

      IM_PHF=1 ; ALLOCATE(COP(0:IM_PHF),PHF(0:IM_PHF))
    OPEN(10,FILE='ATM_AEROSOL.MODEL') ! input File - AEROSOL (CLOUD) in Atmosphere
    OPEN(9,FILE='AER_MOD.PROTOCOL')   ! output File for control

! --- Standard ANGLE grid --- !
       OPEN(8,FILE='./M-C_OPT_BASE/STAND_ANGLE.GRID')
      DO IMEN_PHF=0,1000000
       READ(8,*)ANG
	   IF(ANG > 179.99) EXIT
      END DO
   ALLOCATE (PH_F_I(0:IMEN_PHF),ST_ANG(0:IMEN_PHF))
        REWIND(8)
        DO I=0,IMEN_PHF
        READ(8,*)ANG
       ST_ANG(I)=DCOS(Pi*ANG/180.D0)
        END DO
 ST_ANG(IMEN_PHF)=-1.D0
        CLOSE(8)
! *****  Reading information about  AEROSOL (CLOUD) in Atmosphere ***** !
    READ(10,*) NLAYER ; WRITE(9,*)NLAYER ! Number of Layers with different Phase Functions
   NLAY=ABS(NLAYER)
   DO NL=1,NLAY
   HVOST(NL)='  '
   IF(NL < 10 )THEN
   WRITE(HVOST(NL)(1:1),'(I1)')NL
                  ELSE
   WRITE(HVOST(NL)(1:2),'(I2)')NL
   END IF
   END DO
 ALLOCATE(E55(NLAY),ZDOWN(NLAY),ZUP(NLAY)) !&&&,NUFRAC(NLAY)) ! Ext.(0.55 mkm), Boundaries (KM), Frac.Numb.

ALLOCATE (PROBFUN(NLAY,0:IM4),SCATT(NLAY),ABSS(NLAY),EXTT(NLAY),SERLEG(NLAY,0:NANGLE2))
! This part to DEFINE MAXimal number of fractions in complex aerosol (for dimension) !
   READ(10,*)WN55 ; WN55=-WN55
   NSTR=(NLAY-0.1)/10 +1
   DO N__R=1,NSTR ;  READ(10,*)E555  ; END DO  ! ### 5
      ILM=0
      DO I=1,NLAY
	  READ(10,*)N_TEC
      READ(10,*)IL
      IF(IL > ILM) ILM=IL
      READ(10,*) WWW
        DO J=1,IL
        READ(10,8)AZL
        END DO
 8     FORMAT(A80)
      END DO
    ALLOCATE (A_ZOL(NLAY,ILM),WFRAC(NLAY,ILM)) ! Names of aerosol fractions and weights
    REWIND(10)
    READ(10,*) NLAYER             ! Number of Layers
   READ(10,*)WN55 ; WRITE(9,*)WN55, ' WaveNumber (cm^-1)' ; WN55=-WN55
   READ(10,6555)E55
   WRITE(9,6555)E55 ! Extinctions at 0.55 mkm in each layer
6555 FORMAT(10E12.5)
      DO I=1,NLAY
	  READ(10,*)N_TEC,ZDOWN(I),ZUP(I)      ! Boundaries in km
	  WRITE(9,*)N_TEC,ZDOWN(I),ZUP(I)
      READ(10,*)NFRAC ; WRITE(9,*)NFRAC        ! Number of fractions
      READ(10,*) (WFRAC(I,J),J=1,NFRAC)        ! their Weights
      WRITE(9,*) (WFRAC(I,J),J=1,NFRAC)
                SW=0. ! Sum of weights control
	        DO J=1,NFRAC
			SW=SW+WFRAC(I,J)
	        END DO
			WRITE(9,*)' SUM Weights = ',SW
			IF(SW < 0.9999 .OR. SW > 1.0001) THEN
			WRITE(*,*)I,'-th fraction! ATTENTION: SUM of WEIGHTS = ', SW
			PAUSE
			END IF
        DO J=1,NFRAC
        READ(10,8)A_ZOL(I,J) ; WRITE(9,8)A_ZOL(I,J) ! path to aerosol model in database
        END DO
       END DO
    CLOSE (10)
    CLOSE(9)
! 2 ***** Information about AEROSOL (CLOUD) in Atmosphere has been written ***** !
  END SUBROUTINE FILES_CREATING


  SUBROUTINE LAYER_PROPERTY
! **** Creating Aerosol Optical Properties for the given WaveNumber. **** !
  USE MOD_IN_GAS_CLOUD
  USE A_MOD
  CHARACTER FLNCL*2
  REAL*4,SAVE :: WNOLD,WN_KD
  PARAMETER (N_WN=16) ! &&& Number K-terms.
  DIMENSION WN_KD(N_WN),FLNCL(N_WN) ! &&& Clouds WaveNumbers in K-terms.
  COMMON/W_N/WN

!%%%%%%%%%%%%%%%%%%%%%%
DATA WN_KD/100.,300.,450.,550.,680.,830.,1000.,1150.,1250.,1350.,1600., &
           2000.,2800.,3500.,4300.,5000./
DATA WN_OLD/0.0/
DATA FLNCL/'_1','_2','_3','_4','_5','_6','_7','_8','_9','10','11','12','13','14','15','16'/
!%%%%%%%%% READING !!! %%%%%%%%%%%%%!
IF(WN/=WNOLD)THEN
WNOLD=WN
DO JFL=1,N_WN
IF(WN==WN_KD(JFL))EXIT
END DO
   IRESAE=85*3*I_FCTR
   IREPROB=1001*I_FCTR ! FORTRAN
   IRESRG=85*I_FCTR ! FORTRAN   !  (can be used  < 85 as in QUADRA)
OPEN(2053, ACCESS='DIRECT', FORM='UNFORMATTED',RECL = IRESAE, FILE='./M-C_OPT_BASE/S&A&E.'//FLNCL(JFL),ERR=251)
OPEN(2051, ACCESS='DIRECT', FORM='UNFORMATTED',RECL = IREPROB, FILE='./M-C_OPT_BASE/PROBF.'//FLNCL(JFL),ERR=251)
OPEN(2052, ACCESS='DIRECT', FORM='UNFORMATTED',RECL = IRESRG, FILE='./M-C_OPT_BASE/SRG.'//FLNCL(JFL),ERR=251)
DO JJLL=1,NLAY
      READ(2051,REC=JJLL,ERR=251)PROBFUN(JJLL,:)
      READ(2052,REC=JJLL,ERR=251)SERLEG(JJLL,:)
END DO
      READ(2053,REC=1,ERR=251)SCATT(:),ABSS(:),EXTT(:)
END IF
!%%%%%%%%%%%%%%%%%%%%%%
RETURN
251 WRITE(*,*)' ************ ERROR *************** '
  END

 SUBROUTINE SC_AB_A(Z_I,SC_A,AB_A,L_NUMB)
! Subroutine to calculate aerosol scatt. absorp. coeffitients and layer number.
 USE A_MOD
 DATA ZI_OLD,L_OLD/-1000.0,1/
  SC_A=0.0 ; AB_A=0.0
 IF(Z_I < ZDOWN(1).OR.Z_I > ZUP(NLAY))RETURN
  IF(Z_I <= ZI_OLD)THEN
        DO L=L_OLD,1,-1
        IF(Z_I <= ZUP(L) .AND. ZDOWN(L) <= Z_I ) EXIT
        END DO
  IF(L < 1) RETURN
  ELSE
        DO L=L_OLD,NLAY
        IF(Z_I <= ZUP(L) .AND. ZDOWN(L) <= Z_I) EXIT
        END DO
  IF(L > NLAY) RETURN
 END IF
  SC_A=SCATT(L) ; AB_A=ABSS(L)
  ZI_OLD=Z_I ; L_OLD=L ; L_NUMB=L
 END

	FUNCTION TEMPER(ZZZ)
 USE MOD_IN_GAS_CLOUD
	REAL*8 temper
!*
		DO J=2,JMAX
		IF(ZZZ >= Z(J-1).AND.ZZZ <= Z(J)) GOTO 1
		END DO
		WRITE(*,*)' TEMPER '
		STOP
 1		C1=(Z(J)-ZZZ)/(Z(J)-Z(J-1))
		C2=1.-C1
		TEMPER=C1*T1(J-1)+C2*T1(J)
		END
! ***************************************************** !

 SUBROUTINE I_INTO_LEG(RINT,FN) ! Intensity into Legendre series
 USE MOD_IN_GAS_CLOUD
 DIMENSION FN(0:NANGLE2),RINT(NANGLE2)
     FN=0.
   DO I=1,NANGLE
	  I_=NANGLE2-I+1
      XX=COSINE(I)
      WW=WEIGHT(I)
      FP=RINT(I)
      FM=RINT(I_)
      P_0=1.
      P_1=XX
      FN(0)=FN(0)+WW*P_0*(FP+FM)
      FN(1)=FN(1)+WW*P_1*(FP-FM)
      SGN=-1.    !D0
      A1=1.      !D0
      A2=1.      !D0
      DO J=2,NANGLE2
      A3=A1
      A1=A1+1.   !D0
      A2=A2+2.   !D0
      SGN=-SGN
      PP=(A2*XX*P_1-A3*P_0)/A1
      P_0=P_1
      P_1=PP
      FN(J)=FN(J)+WW*PP*(FP+FM*SGN)
      END DO
      END DO
      A1=0.5     !D0
      DO J=0,NANGLE2
      FN(J)=FN(J)*A1
      A1=A1+1.   !  D0
      END DO
 END
   FUNCTION COS_INIT(FN_,RAND_VAL) ! Random initial COS(Zenith Angle)
USE MOD_IN_GAS_CLOUD
INTEGER,SAVE :: IBG=0
REAL*4,SAVE :: F_PROB_L
REAL*4,SAVE :: X_PROB
PARAMETER (L=50) ! 50 L-number of points in Probability Function table
DIMENSION FN_(0:NANGLE2),F_N(0:NANGLE2),X_PROB(0:L),F_PROB_L(0:NANGLE2,0:L),F_PROB(0:L)
   IF(IBG == 0)THEN ! "Probability Function" for Legendre Polinoms
   IBG=1
  STEP=2./L; A=-1.
  X_PROB(0)=A
  DO I=1,L! Inteval (A=-1,B)
  B=A+STEP*I
  X_PROB(I)=B
  F_N=0.
   DO JQ=1,NANGLE ! Quadrature
! (-1,B) -> (-1,1)
   WES=WEIGHT(JQ)*(B+1.)/2.
!***********
      P_0=1.
   F_N(0)=F_N(0)+(WES+WES)
! (+)POINT
      XX=((B+1.)*COSINE(JQ)+(B-1.))/2.
      P_1=XX
      F_N(1)=F_N(1)+WES*P_1
      A1=1.      !D0
      A2=1.      !D0
      DO J=2,NANGLE2
      A3=A1
      A1=A1+1.   !D0
      A2=A2+2.   !D0
      PP=(A2*XX*P_1-A3*P_0)/A1
      P_0=P_1
      P_1=PP
      F_N(J)=F_N(J)+WES*P_1
      END DO
! (-)POINT
      P_0=1.
      XX=(-(B+1.)*COSINE(JQ)+(B-1.))/2.
      P_1=XX
      F_N(1)=F_N(1)+WES*P_1
      A1=1.      !D0
      A2=1.      !D0
      DO J=2,NANGLE2
      A3=A1
      A1=A1+1.   !D0
      A2=A2+2.   !D0
      PP=(A2*XX*P_1-A3*P_0)/A1
      P_0=P_1
      P_1=PP
      F_N(J)=F_N(J)+WES*P_1
      END DO
!***********
   END DO   ! Quadrature
F_PROB_L(:,I)=F_N
  END DO ! Inteval (A=-1,B)
  X_PROB(L)=1.
  DO J=0,NANGLE2 ! Factor (J+1/2)
  F_PROB_L(J,:)=F_PROB_L(J,:)/(J+0.5)
  END DO
  END IF
 ! Carrent Probability Function
      F_PROB=F_PROB_L(0,:)*FN_(0)
	  DO J=1,NANGLE2
      F_PROB=F_PROB+F_PROB_L(J,:)*FN_(J)
	  END DO
	  ORM=F_PROB(L)
	  F_PROB=F_PROB/ORM
! Simple smoothing
      DO J=1,L
	  IF(F_PROB(J) < F_PROB(J-1))F_PROB(J)=F_PROB(J-1)
	  END DO
! --------------------- !
 IA=0 ; IB=L
 DO WHILE(IB-IA>1)
 IC=(IA+IB)/2
 FC= F_PROB(IC)
 IF(RAND_VAL <= FC) THEN
 IB=IC
 ELSE
 IA=IC
 END IF
 END DO
 FA=F_PROB(IA) ; FB=F_PROB(IB)
 C1=(FB-RAND_VAL)/(FB-FA) ; C2=1.-C1
 COS_INIT=C1*X_PROB(IA)+C2*X_PROB(IB)
  END

! Some modification for the LONGWAVE (no Raylegh etc.)
! 11 Dec., 2003 - 3 Sept., 2004.
      SUBROUTINE MONTE_CARLO(NEWDIR,REFLECT)
!*------------------------------------------------------------------------*
!* Subroutine to obtain photon's trajectory by Monte-Carlo technique      *
!* See in MODULE M_C                                                      *
!* X_(i) is x- coordinate of i-th scattering points, (i=1,2,...,N_POINTS),*
!* Y_(i),Z_(i) are y,z- coordinates, respectively.                        *
!* COS1_(i),COS2_(i),COS3_(i) and P1_(i),...,P4_(i) are cosines of the    *
!*    photon's direction and its Stokes parameters after i-th scattering. *
!* Now it's considered 4-Stokes param. only (for spherical particl. etc.),*
!*  but their number may be easily increased (see COMMON/RAY_ST/ ...).    *
!* NEWDIR  -  ... angles and Stokes parameters after the scattering.      *
!* REFLECT -  ... angles and St. par. after the reflection at the surface.*
!*------------------------------------------------------------------------*
USE M_C ; USE A_MOD
USE MOD_IN_GAS_CLOUD
  REAL*8 A,B,C,CHECK,RAND

  DATA K/0/
         IF(K == 0)THEN
          K=1
         CHECK=A0_**2+B0_**2+C0_**2
          IF(DABS(1.D0-CHECK).GT.1D-10)THEN
      WRITE(*,*)' Sum of the initial cosines = ',CHECK,A0_,B0_,C0_
      WRITE(*,*)' Pause (only once) ! '
         PAUSE
          END IF
         END IF
 ! --------- Start of a new trajectory ------------ !
          I=1
  X_I=X0_ ;  Y_I=Y0_ ;  Z_I=Z0_ ;  X_(I)=X_I ;  Y_(I)=Y_I ; Z_(I)=Z_I
   A=A0_ ;   B=B0_ ;   C=C0_ ;  COS1_(I)=A ;  COS2_(I)=B ;  COS3_(I)=C
Q1=P1_0; Q2=P2_0; Q3=P3_0; Q4=P4_0; P1_(I)=Q1; P2_(I)=Q2; P3_(I)=Q3; P4_(I)=Q4

! ***** 351 NEW TRAJECTORY *****!
 1    I=I+1
 ! ---------  After i-th scatterring --------------- !
!  New Z-coordinate
CALL RNB(RAND)
      DTAU=-DLOG(RAND) ! ***
      DTAU=DTAU*DABS(C) ! Taking into account of the cosine.
      ISIG=1
      IF(C < 0.D0) ISIG=-1
 CALL NEW_Z(Z_I,DTAU,ISIG,IBOUND) !  Variant with X,Y -see in previous versions
    Z_(I)=Z_I
! <<< IBOUND=-1,0,+1, if a new point is ON a SURFACE, INSIDE and OUTSIDE ATMOSPHERE >>>
      N_POINTS=I
      IF(IBOUND.EQ.1)RETURN ! Photon outside atmosphere - End of this trajjectory.

     IF(IBOUND.EQ.0)THEN
     CALL NEWDIR(X_I,Y_I,Z_I,A,B,C,Q1,Q2,Q3,Q4)! Scattering inside atmosphere - New direction.
                       ELSE
  					  Z_I=Z(1) ; Z_(I)=Z_I
      CALL REFLECT(X_I,Y_I,Z_I,A,B,C,Q1,Q2,Q3,Q4)! Reflectance at the surface - New direction.
       Z_I=ZDOWN(1)
                      END IF
    COS1_(I)=A ;  COS2_(I)=B ;  COS3_(I)=C
    P1_(I)=Q1  ;  P2_(I)=Q2  ;  P3_(I)=Q3 ;  P4_(I)=Q4
      IF(I.GE.IMAX_MC)THEN ! End of trajectory (artificial).
          COS1_(I)=0.D0
          COS2_(I)=0.D0
          COS3_(I)=0.9999999
          COS1_(I+1)=0.D0
          COS2_(I+1)=0.D0
          COS3_(I+1)=0.9999999
          Z_(I+1)=Z(JMAX)
          X_(I+1)=X_(I)
          Y_(I+1)=Y_(I)
          N_POINTS=I+1
          RETURN
       END IF
                   GO TO 1
        END


  SUBROUTINE NEW_DIRECTION(X_I,Y_I,Z_I,A,B,C,Q1,Q2,Q3,Q4)
! New Photon's direction (A,B,C (old) -> A,B,C (new) ) Aerosol+Molecular scattering.
  USE A_MOD
  USE M_C
  REAL*8 A,B,C,W1,W2,W3,W4,W5,CUG,RAND
  COMMON/SCAT_POINT/LNUMB
! ----------------------------- *
!  New direction of the photon (Without polarization) *
  1     CALL RNB(RAND)
    W1=1.D0-2.D0*RAND  ! COS & SIN : see Marchuk p.10 . ! ***###
   CALL RNB(RAND)
          W2=1.D0-2.D0*RAND
          W3=W1**2 +W2**2
          IF(W3.GT.1D0.OR.W3.EQ.0.d0)GOTO 1
          CALL RNB(RAND)
             IND=(IM4-1)*RAND
      	         CUG=PROBFUN(LNUMB,IND) ! Aerosol scattering
          W3=DSQRT((1.D0-CUG**2)/W3)
          W4=W1*W3
          W5=W2*W3
          W1=A*W4-B*W5
!  New cosines *
          W3=(CUG-W1/(1.D0+DABS(C)))
          A=A*W3+W4
          B=B*W3-W5
          C=C*CUG-W1*SIGN(1.D0,C)
  END

      SUBROUTINE LAMBERT(X_I,Y_I,Z_I,A,B,C,P1,P2,P3,P4)
USE M_C
! Photon's reflection at the surface - Lambert's low.
! New direct.and polar.(A,B,C,P1,...,P4 old > A,B,C,P1,...,P4 new).
! X_I,Y_I,Z_I- for future applications.
        REAL*8 A,B,C,W1,W2,W3,RAND
      P1=1.
      P2=1.
      P3=0.
      P4=0.
CALL RNB(RAND)
      W1=RAND
CALL RNB(RAND)
      W2=RAND
      C=MAX(W1,W2)  ! modelling P(X)=2X (see Ermakov p.25)
  1   CALL RNB(RAND)
          W1=1.D0-2.D0*RAND
 CALL RNB(RAND)
          W2=1.D0-2.D0*RAND
          W3=W1**2 +W2**2
          IF(W3.GT.1D0)GOTO 1
          W3=DSQRT((1.D0-C**2)/W3)
          A=W1*W3
          B=W2*W3
       END SUBROUTINE LAMBERT

 SUBROUTINE NEW_Z(Z_Z,DTAU,ISIG,IBOUND) ! 22 Oct.,2004
! ------------------------------------------------------------ *
! To calculate Z-coordinate for the plane-parallel atmosphere. *
! C0 - cosine between Z-axis and a photon's direction.         *
! IBOUND = (-1,0,+1) if Z at the SURFACE and IN or OUT atmosph.*
! DTAU=TAU*ABS(C0)  ; ISIG=C0/ABS(C0) (for UP/DOWN)            *
! ------------------------------------------------------------ *
 USE A_MOD
 USE MOD_IN_GAS_CLOUD
 REAL*4,SAVE :: WOLD
 REAL*4,SAVE :: SC_LOC
 INTEGER,SAVE :: L_NUM
      DIMENSION SC_LOC(200),L_NUM(200)
      COMMON/W_N/WAVE ! Wavenumber
      COMMON/SCAT_POINT/LNUMB
      DATA WOLD/-1./

 IF(WOLD.NE.WAVE)THEN ! New wavenumber
       WOLD=WAVE
       DO L=1,NLAY
        ZZZ=(ZDOWN(L)+ZUP(L))*0.5
     CALL SC_AB_A(ZZZ,S_SCAT,S_ABS,LNB)
        SC_LOC(L)=S_SCAT ; L_NUM(L)=LNB
       END DO
 END IF
! ------
 IF(Z_Z > ZUP(NLAY)) Z_Z=ZUP(NLAY)-0.000001
 IF(Z_Z < ZDOWN(1)) Z_Z=ZDOWN(1)+ 0.000001
 DO IL=1,NLAY ; IF(Z_Z >= ZDOWN(IL) .AND. Z_Z <= ZUP(IL)) EXIT ; END DO

IF(ISIG == 1) THEN

! --------- UP ---------- !
   DO IC=1,1000000
       ZIL=ZUP(IL)
       DELTAZ=ZIL-Z_Z
       DT=DELTAZ*SC_LOC(IL)
       DTAU=DTAU-DT
	   IF(DTAU <= 0. .OR. IL == NLAY) EXIT
        IL=IL+1
        Z_Z=ZIL
   END DO
     IF(DTAU > 0.) THEN
       IBOUND =1 ;  Z_Z=1000. ! Space
     ELSE
       Z_Z=ZIL+DTAU/SC_LOC(IL)
       IBOUND=0
      LNUMB=L_NUM(IL)
     END IF
ELSE

! -------- DOWN --------- !
     DO IC=1,1000000
         ZIL=ZDOWN(IL)
         DELTAZ=Z_Z-ZIL
          DT=SC_LOC(IL)*DELTAZ
         DTAU=DTAU-DT
        IF(DTAU <= 0. .OR. IL == 1)EXIT
        Z_Z=ZIL
        IL=IL-1
       END DO
      IF(DTAU > 0) THEN ! Photon has fallen down on the surface
          IBOUND=-1 ;   Z_Z=Z(1) ;      LNUMB=1
	  ELSE
         Z_Z=ZIL-DTAU/SC_LOC(IL)
         IBOUND=0
      LNUMB=L_NUM(IL)
	  END IF

 END IF

 END

       SUBROUTINE  RNB(RAND_)
 USE M_C
      REAL*8 dm37,RAND_
      DATA dm37 /137438953472.D0/
           X_RAND=X_RAND*3125.D0
          X_RAND=DMOD(X_RAND,dm37)
          RAND_=X_RAND/dm37
      END

! *** 22.08.2023 ***!

!*** ! Calculates initial fluxes-intensities (0-scattering) for Monte-Carlo  ***
! and IMMEDIATELY (at each wavenumber point) scattering of the radiation.    ***
 SUBROUTINE FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NROZ,JMAXX)
 !(10 Sept.,2003- 5 Aug.,2004- 17 Nov.,2004)
!*** The LONGWAVE flux calculations by the "fast" technique                  ***
!*** The trapezoidal rule for integration over altitude                      ***
!*** "UNIVERSAL" (adapting spectral grid H=Const depending on dopler width)  ***
!*** Automatic taking into account of the number of atmosph. levels number   ***
!*** Wiscombe's algorithm (linear in tau). Double Precision. Zi = Zj         ***
!*******************************************************************************
 USE MOD_IN_GAS_CLOUD ; USE A_MOD
	IMPLICIT REAL*8 (A-H,O-Z)
 REAL*4 PL,ZUKA,S_SCAT,SUCA,ABDO,EMISSIVITY, &
 AB_AEROS,SC_AEROS,EXT_AEROS,S_3,FN,Tsurf,CC0,COS_INIT
 INTEGER,SAVE :: IBEG,JZM,JZM_1
 REAL*8,SAVE :: CL0,CL1
	PARAMETER(DELT_CRT=0.0001,TLIM=0.001,	&
	ND=200,NDCLO=120) ! VVV 90 ND,NDCLO - MAX number of Z-points in atmosphere and cloud
!* ND	- maximal number (for DIMENSION) of the J-grid.
!* NQ	- ... of the wavenumber grid (per inverse cm).
!*.... The name of FILE with absorption coefficients at J-levels:
!* ---------------------------------------------------------------*
 PARAMETER (STRONG_ABS=50.) ! <<<S>>> if GASEOUS Absorption > STRONG_ABS scattering will be omitted
	DIMENSION FLUXUP(JMAXX),FLUXDO(JMAXX),FLUXUP_(JMAXX),FLUXDO_(JMAXX) &
     ,EMIS_ON(NANGLE),FN(0:NANGLE2),FNCLO(0:NANGLE2,NDCLO) &
     ,FPROB_I(NDCLO) &! Array for simulation of the starting point
     ,Emis_Up(ND,NANGLE),Emis_Down(ND,NANGLE),EMISSIVITY(ND,NANGLE2),&
  E_TAU(ND,NANGLE),FLU(ND),FLD(ND),TAUMA_E(ND),CL0(ND),CL1(ND),TAUMA_A(ND), & ! <<<S>>>
  TAUMA(ND), &
                  FLU_(ND),FLD_(ND), FLUS(ND),FLDS(ND) ! For scattered radiation
 COMMON/CL/S_3(ND),AB_AEROS(ND),SC_AEROS(ND),EXT_AEROS(ND),ISUM,N_CATAL(ND)! N_CATAL eg 1,2,5,6,..
                !N_CATAL Levels with scattering (Number of such levels = ISUM)
 COMMON/KD/PL(150)
		DATA IBEG/0/
! --------------------------------- !
		ABDO=0.   ! ATTENTION !!!
		IF(IBEG == 0)THEN
		IBEG=1
        CALL QUAD_DEF ! Quadrature inserting
						JZM=JMAX
         IF(JZM>ND)THEN
          WRITE(*,*)'******* Increase ND-parameter !!! *******'
          STOP
         END IF

			JZM_1=JZM-1
!* --------------------------------------------- *
!******************* CL0,CL1 for TAUMA - matrix ****************
	CL0(1) = 0.D0
	CL1(1) = 0.D0
!*
	DO 1 J = 2, JZM
!*
	Z1=Z(J-1)
	Z2=Z(J)
	PSIL=Z1
	PSIL1=Z2
	GA=Z2-Z1
	GB=0.5*(Z2+Z1)*GA
	CC=PSIL-PSIL1
	CL1(J)=(PSIL*GA-GB)/CC
	CL0(J)=(GB-PSIL1*GA)/CC
!*
1	CONTINUE

! COS(1) > COS(2) ... etc.
         DO WHILE (IPR/=NANGLE)
		 IPR=1
		  DO IA=2,NANGLE
          CO1=COSINE(IA-1)
		  CO2=COSINE(IA)
		    IF(CO1<CO2)THEN
            W1=WEIGHT(IA-1)
			WEIGHT(IA-1)=WEIGHT(IA)
			WEIGHT(IA)=W1
			COSINE(IA-1)=CO2
			COSINE(IA)=CO1
			ELSE
			IPR=IPR+1
			END IF
		  END DO
		 END DO
	END IF
!C ***********************************
!C
!C *** STEP BY STEP
 !   ***  TAKING INTO ACCOUNT OF AEROSOL SCATTERING AND ABSORPTION ***
    CALL LAYER_PROPERTY ! Setting of aerosol opt. prop. at given Wave-Number
         TAUMA(JZM)=0.
		 TAUMA_E(JZM)=0.
         TAUMA_A(JZM)=0. ! <<<S>>>
    DO J=JZM_1,1,-1
TAUMA(J)=TAUMA(J+1)+0.5*(RABMA(J+1)+RABMA(J))*(Z(J+1)-Z(J))
         ZUKA=(Z(J)+Z(J+1))*0.5
         CALL SC_AB_A(ZUKA,S_SCAT,SUCA,LNUMB)
	    TAUMA_E(J)=TAUMA_E(J+1) + (S_SCAT+SUCA)*(Z(J+1)-Z(J))
        TAUMA_A(J)=TAUMA_A(J+1) + SUCA*(Z(J+1)-Z(J)) ! <<<S>>>
    AB_AEROS(J)=SUCA
    SC_AEROS(J)=S_SCAT
    EXT_AEROS(J)=S_SCAT+SUCA
     END DO
!  upon the wave number intervals where Planck=const ***
    SUM_SCAT=0.d0
ISUM=0
III=1
IF(EXT_AEROS(III) > 0.0) THEN
ISUM=ISUM+1 ; N_CATAL(ISUM)=III
END IF
	 DO III=2,JZM_1
 IF(EXT_AEROS(III) > 0.0.OR.EXT_AEROS(III-1) > 0.0) THEN
 ISUM=ISUM+1 ; N_CATAL(ISUM)=III
 END IF
     END DO
     S_3=0.
	 DO I=2,ISUM
	 J1=N_CATAL(I-1) ; J2=N_CATAL(I)
	  IF(J2-J1 == 1)THEN
	  S_3(I)=S_3(I-1)+SC_AEROS(J1)*(Z(J2)-Z(J1))
	  END IF
	 END DO
!*
		DO JJ=1,JZM
                FLU(JJ)=0.D0
                FLD(JJ)=0.D0
                FLU_(JJ)=0.D0
                FLD_(JJ)=0.D0
		END DO
  PL_Surf=PL(1) ! <===========================
! ------------------------------
!* Loop over wavenumber points *

! *****Include or not scattering *****
INC_SCA=1 ; J=N_CATAL(ISUM) ! <<<S>>>
!      *************************
IF(INC_SCA==1) THEN ! <<<S>>>
 TAUMA=TAUMA+TAUMA_E   ! Initial clouud/aerosol optical thicknesses (EXTINCTION)
 ELSE
 TAUMA=TAUMA+TAUMA_A
 END IF	! <<<S>>>
    DO J=JZM_1,1,-1 !
     TAT=-(TAUMA(J)-TAUMA(J+1))
     DO IAN=1,NANGLE
     E_TAU(J,IAN)=DEXP(TAT/COSINE(IAN))
     END DO
    END DO
! ------------- Emissivities --------------------
!  Downward
      DO J=1,JZM_1
		BB2=PL(J)
		BB1=PL(J+1)
		TATA2=TAUMA(J)
		TATA1=TAUMA(J+1)
		DELT_TAU=TATA2-TATA1
! ------- !
FCTR=1.
IF(EXT_AEROS(J) > 0.000001.AND.INC_SCA==1) &
                 FCTR=(RABMA(J)+AB_AEROS(J))/(RABMA(J)+EXT_AEROS(J)) ! <<<S>>>
! ------- !
			DO IAN=1,NANGLE
			X_X=COSINE(IAN)
		IF(DELT_TAU < DELT_CRT)THEN
		EMIS=0.5D0*(BB1+BB2)*DELT_TAU      ! X_X/X_X
		ELSE
		ALPHA=(BB1*TATA2-BB2*TATA1)/DELT_TAU
		BETTA=(BB2-BB1)/DELT_TAU
		EMIS=(ALPHA-BETTA*X_X+BETTA*TATA2 +	&
	(BETTA*X_X-BETTA*TATA1-ALPHA)*DEXP((TATA1-TATA2)/X_X))*X_X
		END IF
        Emis_Down(J,IAN)=EMIS*WEIGHT(IAN) & ! All weights
 * FCTR ! Correction of the emissivity (K+A)/(K+A+S)
	        	END DO
       END DO
! Upward
			DO IAN=1,NANGLE ! Surface
            Emis_Up(1,IAN)=PL_Surf*COSINE(IAN)*WEIGHT(IAN)! All weights
                        END DO
	DO J=2,JZM
		BB2=PL(J-1)
		BB1=PL(J)
		TATA2=TAUMA(J-1)
		TATA1=TAUMA(J)
		DELT_TAU=TATA2-TATA1
		FCTR=1.0
! ------- !
FCTR=1.
IF(EXT_AEROS(J-1)>0.000001.AND.INC_SCA == 1) & ! <<<S>>>
 FCTR=(RABMA(J)+AB_AEROS(J-1))/(RABMA(J)+EXT_AEROS(J-1))! <<<S>>>
! ------- !
			DO IAN=1,NANGLE
			X_X=COSINE(IAN)
		IF(DELT_TAU < DELT_CRT)THEN
		EMIS=0.5*(BB1+BB2)*DELT_TAU
		ELSE
		ALPHA=(BB1*TATA2-BB2*TATA1)/DELT_TAU
		BETTA=(BB2-BB1)/DELT_TAU
		EMIS=(ALPHA+BETTA*X_X+BETTA*TATA1 -	&
 	(BETTA*X_X+BETTA*TATA2+ALPHA)*DEXP((TATA1-TATA2)/X_X))*X_X
		END IF
                Emis_Up(J,IAN)=EMIS*WEIGHT(IAN) & ! All weights
 * FCTR ! Correction of the emissivity (K+A)/(K+A+S)
            END DO
        END DO
 ! --------------------------------------------------------
!****** Downward fluxes
        EMIS_ON=0.D0
	DO  I = JZM_1,1,-1  !   Loop over ZI-points ***
             F_=0.D0
                      DO IAN=1,NANGLE
                      EMIS_ON(IAN)=Emis_DOWN(I,IAN)+E_TAU(I,IAN)*EMIS_ON(IAN)
             F_=F_+EMIS_ON(IAN)
 EMISSIVITY(I,NANGLE2+1-IAN)=EMIS_ON(IAN)/(WEIGHT(IAN)*COSINE(IAN))
                      END DO
	              FLD(I)=FLD(I)+F_
     END DO! End of the loop over Zi points *
	 F_DOWN_S=F_*2. ! Downward Flux at the surface
!****** Upward fluxes ******
             F_=0.D0
              DO IAN=1,NANGLE
EMIS_ON(IAN)=(1.-ABDO)*Emis_UP(1,IAN)+ABDO*F_DOWN_S*COSINE(IAN)*WEIGHT(IAN)  ! ALBEDO is taken into account, LAMBERT low
              F_=F_+EMIS_ON(IAN)
 EMISSIVITY(1,IAN)=EMIS_ON(IAN)/(WEIGHT(IAN)*COSINE(IAN))
              END DO
              FLU(1)=FLU(1)+F_
	DO  I = 2,JZM  !   Loop over ZI-points ***
             F_=0.D0
              DO IAN=1,NANGLE
             EMIS_ON(IAN)=Emis_UP(I,IAN)+E_TAU(I-1,IAN)*EMIS_ON(IAN)
              F_=F_+EMIS_ON(IAN)
 EMISSIVITY(I,IAN)=EMIS_ON(IAN)/(WEIGHT(IAN)*COSINE(IAN))
              END DO
              FLU(I)=FLU(I)+F_
     END DO! End of the loop over Zi points *
! **************************** !

 IF(INC_SCA==1) THEN !  <<<S>>> Including absorption or not
! Emissivities into Legendre series and recording !
  SU_I=0.
	 DO III=1,ISUM
 N_C=N_CATAL(III)

 CALL I_INTO_LEG(EMISSIVITY(N_C,:),FN)
 FNCLO(:,III)=FN
 FPROB_I(III)=FN(0)
  IF(III == 1) THEN
        FB=FN(0)
  ELSE
 FA=FB ; FB=FN(0)
 END IF
 IF(III > 1)SU_I=SU_I+0.5*(FA+FB)*(Z(N_C)-Z(N_C-1))*SC_AEROS(N_C-1)
    END DO
 G=SU_I
 G=G*PI*2.
SUM_SCAT=SUM_SCAT+G
 CALL CHAIN2(FLUS,FLDS,ABDO,NROZ,FPROB_I,FNCLO,JMAXX) !  Monte-Carlo
       DO J=1,JZM
       FLU_(J)=FLU_(J)+FLUS(J)
       FLD_(J)=FLD_(J)+FLDS(J)
       END DO
 END IF ! <<<S>>>
! ---------------------- !
        DO J=1,JZM
       FLUXUP(J)=FLUXUP(J)+FLU(J)
       FLUXDO(J)=FLUXDO(J)+FLD(J)
       FLUXUP_(J)=FLUXUP_(J)+FLU_(J)
       FLUXDO_(J)=FLUXDO_(J)+FLD_(J)
       END DO
END
!***************************************************
SUBROUTINE CHAIN2(FLUXUP,FLUXDO,A_BEDO,NROZ,FPROB_I,FNCLO,JMAXX) ! 18.11.2004
!*********************************************************************
! To calculate LW fluxes in the SCATTERING PLANE-PARALLEL atmosphere !                                      *
!*********************************************************************
USE MOD_IN_GAS_CLOUD
USE M_C
USE A_MOD
	IMPLICIT REAL*8 (A-H,O-Z)
REAL*4 PL,A_BEDO, &
AB_AEROS,SC_AEROS,EXT_AEROS,S_3,FN,Tsurf,CC0,COS_INIT
     PARAMETER (ND=200,NDCLO=120,QQ_MIN=1E-4) ! see in FLUX_0_SCAT
	 ! *** if QQ<QQ_MIN = > go to NEW Trajectory *** !

!* ---------------------------------------------------------------*
  SAVE :: LMIN,CL1I,CL0I,CL1,CL0
  DIMENSION FLUXUP(JMAXX),FLUXDO(JMAXX),DIST(JMAX),AERRAB(JMAX) &
  ,LMIN(ND),CL1I(ND),CL0I(ND),CL1(ND),CL0(ND),TAUMA(ND) &
  ,FPROB_I(NDCLO),FNCLO(0:NANGLE2,NDCLO) & ! Array for simulation of the starting point
  ,FN(0:NANGLE2),FN1(0:NANGLE2),FN2(0:NANGLE2)
   DIMENSION AAC(JMAX)
       COMMON/R_A/VSTART,VFINISH
        COMMON/CL/S_3(ND),AB_AEROS(ND),SC_AEROS(ND),EXT_AEROS(ND),ISUM,N_CATAL(ND)
COMMON/KD/PL(150)
       EXTERNAL NEW_Z,NEW_DIRECTION,LAMBERT
	   DATA IBG/0/
    SUM_PH=0.D0
	DO J=1,JMAX ; FLUXUP(J)=0.D0 ; FLUXDO(J)=0.D0 ; END DO
	   IF(IBG.EQ.0)THEN
	   IBG=1
 ALLOCATE (LMINV(IMAX_M_),CL0V(IMAX_M_),CL1V(IMAX_M_),STAT=IERR)
         IF(IERR/=0)THEN
          WRITE(*,*)' Allocaion is wrong (CHAIN) !!!'
          STOP
         END IF
 !******************* CL0,CL1 for TAUMA - matrix  ****************
       CL0(1) = 0.
      CL1(1) = 0.
      DO J = 2, JMAX
      Z1=Z(J-1)
      Z2=Z(J)
      PSIL=Z1
      PSIL1=Z2
      GA=Z2-Z1
      GB=0.5*(Z2+Z1)*GA
      CC=PSIL-PSIL1
      CL1(J)=(PSIL*GA-GB)/CC
      CL0(J)=(GB-PSIL1*GA)/CC
     END DO

!  ***  I- MESH  ***
      DO 12 I=1,JMAX
      ZZZ=Z(I)
      LMIN(I)=JMAX-1
      DO 13 J=2,JMAX
      IF(ZZZ.LT.Z(J))GOTO 14
   13 CONTINUE
      GOTO 12
   14 LMIN(I)=J-1
   12 CONTINUE
      DO 16 I=1,JMAX
      J=LMIN(I)
      Z1=Z(J)
      Z2=Z(J+1)
      PSIL=Z1
      PSIL1=Z2
      ZZ=Z(I)
      GA=ZZ-Z1
      GB=0.5*(ZZ+Z1)*GA
      CC=PSIL-PSIL1
      PSIZ=ZZ
      CL0I(I)=(GB-PSIL1*GA)/CC
      CL1I(I)=(PSIL*GA-GB)/CC
   16 CONTINUE
   JZM_1=JMAX-1
	   END IF

!*  ***  TAU - MATRIX ***
      TAUMA(1)=0.
      DO  J = 2, JMAX
TAUMA(J)=TAUMA(J-1)+(RABMA(J)+RABMA(J-1))/2.*(Z(J)-Z(J-1))
      END DO
!    TAKING INTO ACCOUNT OF AEROSOL ABSORPTION ***
         SUCTAU=0.
         DO J=2,JMAX
        AERRAB(J)=AB_AEROS(J-1)
         SUCTAU=SUCTAU+(Z(J)-Z(J-1))*AB_AEROS(J-1)
           TAUMA(J)=TAUMA(J)+SUCTAU
         END DO
          AERRAB(1)=AERRAB(2)
 ! --------  CHAIN - CYCLE ---------------------------------
!  ***  Z_ - MESH ***
      ROZ=2.0/NROZ !      ! Attention - 2-ka
      DO 128 NMONTE=1,NROZ
!### if(nmonte/10*10.eq.nmonte)write(*,*)nmonte,nroz
! ****************************************** !
! ------- Starting Z0_ point -------- !
   X0_=0. ; Y0_=0.
   CALL RNB(RAND)
   PZ=RAND*S_3(ISUM)
  DO J7=2,ISUM
 IF(S_3(J7) >= PZ)EXIT
 END DO
 IF(J7 > ISUM) I=ISUM
 C1=(S_3(J7)-PZ)/(S_3(J7)-S_3(J7-1)) ; C2=(1.-C1)
 J2=N_CATAL(J7) ; J1=J2-1
 Z0_=C1*Z(J1)+C2*Z(J2)
QQ=ROZ*(FPROB_I(J7-1)*C1+FPROB_I(J7)*C2)*S_3(ISUM)
       DO N_L=1,1000000
        IF(Z0_ <= ZUP(N_L) .AND. ZDOWN(N_L) <= Z0_ ) EXIT
        END DO
 C2=(ZUP(N_L)-Z0_)/(ZUP(N_L)-ZDOWN(N_L)) ; C1=1.-C2
! ----------- Initial Photon's direction C0_ ------------ !
   FN1=FNCLO(:,J7-1)
   FN2=FNCLO(:,J7)
   FN=FN1*C1+FN2*C2
    DO J7=0,NANGLE2
	FN(J7)=FN(J7)*SERLEG(N_L,J7)
    END DO
   CALL RNB(RAND) ; CC0=RAND
 C0_=COS_INIT(FN,CC0)
    B0_=0.
      A0_=DSQRT(1.D0-C0_**2)
! -------------------------------------------------------- !
    CALL MONTE_CARLO(NEW_DIRECTION,LAMBERT)
!------------------ Venus 29 August, 2023  --------! Trajectory cut of if SURFACE.
DO JJJ=1,N_POINTS
IF(Z_(JJJ)==0.0)EXIT
END DO
N_POINTS2=JJJ+1
IF(N_POINTS2<N_POINTS)N_POINTS=N_POINTS2
                                             Z_(N_POINTS)=Z(JMAX)
!* -----------------------------------------------------------
!  *** PARAMETERS FOR IV - MESH ***
!  *** LMINV( ),CL0V( ),CL1V( ) ***
      DO 322 IV=1,N_POINTS
      ZZZ=Z_(IV)
      LMINV(IV)=JMAX-1
      DO 323 J=2,JMAX
      IF(ZZZ.LT.Z(J))GOTO 324
  323 CONTINUE
      GOTO 322
  324 LMINV(IV)=J-1
  322 CONTINUE
      DO 316 I=1,N_POINTS
      J=LMINV(I)
      Z1=Z(J)
      Z2=Z(J+1)
      ZZ=Z_(I)
      CC=Z2-Z1
      CCZ=ZZ-Z1
      CCZZ=(ZZ+Z1)*CCZ*0.5
      CL0V(I)=(Z2*CCZ-CCZZ)/CC
      CL1V(I)=(CCZZ-Z1*CCZ)/CC
  316 CONTINUE
!  ***********************************
!  ***  ABSORPTION ONLY FOR WEITS ***
!
      DO 21 IV=2,N_POINTS
      J1=LMINV(IV-1)
      J2=LMINV(IV)
             TUKA1=(Z_(IV-1)-Z(J1))*AB_AEROS(J1)
             TUKA2=(Z_(IV)-Z(J2))*AB_AEROS(J2)
      CCL00=CL0V(IV-1)
      CCL10=CL1V(IV-1)
      CCL01=CL0V(IV)
      CCL11=CL1V(IV)
             ABSTE=ABS(COS3_(IV-1))
                   IF(ABSTE.LT.1E-9) THEN
                      TETCOS=1E9
                      IF(COS3_(IV-1).LT.0.)TETCOS=-TETCOS
                   ELSE
      TETCOS=1./COS3_(IV-1)
                   END IF

            IF(ABSTE.LT.1.0E-04)THEN
                                     ABSTE = 0.0
                                ELSE
                                     ABSTE = 1.0/ABSTE
                                ENDIF

      ZZZ=Z_(IV-1)
      ZZZZ=Z_(IV)
      IF=0    ;  IB=0
      DIST=0. ;  PIST=0.
      DO IIII=1,JMAX
    FFF=(ZZZ-Z(IIII))*(ZZZZ-Z(IIII))
      IF(FFF.LE.0..AND.IB.EQ.0)IB=IIII
      IF(FFF.LE.0.)IF=IIII
      END DO

      ALB=1.
      IF(Z_(IV).EQ.Z(1))ALB=A_BEDO

!*  ----------------------------------------
      TAU1=TAUMA(J1)+RABMA(J1)*CCL00+RABMA(J1+1)*CCL10
      TAU2=TAUMA(J2)+RABMA(J2)*CCL01+RABMA(J2+1)*CCL11
                  TAU1=TAU1+TUKA1
                  TAU2=TAU2+TUKA2
!   ***  Calculation at Zi - LEVELS ***
      IF(IB.EQ.0)GOTO 1325
      DO 3222 II=IB,IF
      JJ=LMIN(II)
      TAU=TAUMA(JJ)+RABMA(JJ)*CL0I(II)+RABMA(JJ+1)*CL1I(II)
             TAU=TAU+AERRAB(JJ+1)*(Z(II)-Z(JJ))
 ETATATETCOS=0.
 TATATETCOS=(TAU1-TAU)*TETCOS
 IF(TATATETCOS<0.0)ETATATETCOS=DEXP(TATATETCOS)
       WES=QQ*ETATATETCOS
      DIST(II)=DIST(II)+WES
   3222 CONTINUE
 1325 CONTINUE
      YL=(TAU1-TAU2)*TETCOS
QQYL=0.
IF(YL<0.)QQYL=DEXP(YL)
            QQ=QQ*QQYL*ALB
            IF(IB.EQ.0)GOTO 21
      DO II=IB,IF
      IF(TETCOS.GT.0.)THEN
       FLUXUP(II)=FLUXUP(II)+DIST(II)
                      ELSE
       FLUXDO(II)=FLUXDO(II)+DIST(II)
                      END IF
      END DO
IF(WES<1E-6)GOTO 128
   21 CONTINUE   ! *** Loop over trajectories points

  128 CONTINUE   !  *** Loop over photons in the given range

!  WRITE(*,*)' Scattered Energy ',SUM_PH*2.0*PI  ; pause
      END


