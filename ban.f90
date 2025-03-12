 SUBROUTINE BAND_A_1 (NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *
!* using K-distribution technique IN 10 - 200 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL,TAUMA
COMMON/KD/PL(150)
! --------------------------------------------------------- *
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'a',NCH) ! Planck function for region 'a'
	END DO
! * ---------------  H2O ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       H2O1=0.; H2O2=0.
       IF(IH2O > 0)THEN
      H2O1=RO1(IH2O,J) ;  H2O2=RO1(IH2O,J+1)
      END IF
!        H2O
      	A1=0.
      IF(NCH == 1)THEN
VOL1=ChA_1_H2O(P_LN(J))*H2O1
VOL2=ChA_1_H2O(P_LN(J+1))*H2O2
RABMA(J+1)=VOL2
 END IF

! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChA_1_H2O(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChA_1_H2O
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'a.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChA_1_H2O =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChA_1_H2O =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChA_1_H2O=EXP(ABSOR_LOG)
! --------------------------------------------
          END


 SUBROUTINE BAND_A_2(NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *
!* using K-distribution technique IN 10 - 200 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL,TAUMA
COMMON/KD/PL(150)
! --------------------------------------------------------- *
RABMA=0.
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'a',NCH) ! Planck function for region 'a'
	END DO
! * ---------------  H2O ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       H2O1=0.; H2O2=0.
       IF(IH2O > 0)THEN
      H2O1=RO1(IH2O,J) ;  H2O2=RO1(IH2O,J+1)
      END IF
!        H2O
      	A1=0.
      IF(NCH == 2)THEN
VOL1=ChA_2_H2O(P_LN(J))*H2O1
VOL2=ChA_2_H2O(P_LN(J+1))*H2O2
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChA_2_H2O(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChA_2_H2O
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'a.2_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChA_2_H2O =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChA_2_H2O =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChA_2_H2O=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_B_1 (NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *
!* using K-distribution technique IN 10 - 200 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL,TAUMA
COMMON/KD/PL(150)
! --------------------------------------------------------- *
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'b',NCH) ! Planck function for region 'b'
	END DO
! * ---------------  H2O ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       H2O1=0.; H2O2=0.
       IF(IH2O > 0)THEN
      H2O1=RO1(IH2O,J) ;  H2O2=RO1(IH2O,J+1)
      END IF
!        H2O
      	A1=0.
      IF(NCH == 1)THEN
VOL1=ChB_1_H2O(P_LN(J))*H2O1
VOL2=ChB_1_H2O(P_LN(J+1))*H2O2
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChB_1_H2O(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChB_1_H2O
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'b.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChB_1_H2O =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChB_1_H2O =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChB_1_H2O=EXP(ABSOR_LOG)
! --------------------------------------------
          END


 SUBROUTINE BAND_B_2(NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *
!* using K-distribution technique IN 10 - 200 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL,TAUMA
COMMON/KD/PL(150)
! --------------------------------------------------------- *
RABMA=0.
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'b',NCH) ! Planck function for region 'A'
	END DO
! * ---------------  H2O ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       H2O1=0.; H2O2=0.
       IF(IH2O > 0)THEN
      H2O1=RO1(IH2O,J) ;  H2O2=RO1(IH2O,J+1)
      END IF
!        H2O
      	A1=0.
      IF(NCH == 2)THEN
VOL1=ChB_2_H2O(P_LN(J))*H2O1
VOL2=ChB_2_H2O(P_LN(J+1))*H2O2
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChB_2_H2O(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChB_2_H2O
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'b.2_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChB_2_H2O =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChB_2_H2O =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChB_2_H2O=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_C(NCH)
USE MOD_IN_GAS_CLOUD  ! RO [mol/(cm^2*km].
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *
!* using K-distribution technique IN 400 - 500 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'c',NCH) ! Planck function for region 'W'
	END DO
! * ------------------------------------- *
         JMAX_1=JMAX-1
          DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)  ! RO [mol/(cm^2*km].
      END IF
!     ------------- CO2 ------------------ !
      	A1=0.
      IF(NCH == 1)THEN
VOL1=ChC_CO2(P_LN(J))*CO21
VOL2=ChC_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
  	    END IF
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChC_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChC_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1
! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
      OPEN (10,FILE='./F/'//'c.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChC_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChC_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChC_CO2=EXP(ABSOR_LOG)
          END

 SUBROUTINE BAND_D(NCH) ! 24 July,2024.
USE MOD_IN_GAS_CLOUD  ! RO [mol/(cm^2*km].
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *
!* using K-distribution technique IN 400 - 500 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'d',NCH) ! Planck function for region 'W'
	END DO
! * ------------------------------------- *
         JMAX_1=JMAX-1
          DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)  ! RO [mol/(cm^2*km].
      END IF
!     ------------- CO2 ------------------ !
      	A1=0.
      IF(NCH == 1)THEN
VOL1=Chd_CO2(P_LN(J))*CO21
VOL2=Chd_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
  	    END IF
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END
         FUNCTION Chd_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,Chd_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1
! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
      OPEN(10,FILE='./F/'//'d.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Chd_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Chd_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Chd_CO2=EXP(ABSOR_LOG)
          END
 SUBROUTINE BAND_X_1(NCH)  ! 24 Oct.,2024.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 500 - 600 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
 write(*,*)'BAND_X_1'
       DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'e',NCH) ! Planck function for region 'C'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 1)THEN
VOL1=ChX_1_CO2(P_LN(J))*CO21
VOL2=ChX_1_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
 write(*,*)'*** BAND_X_1'
       END

         FUNCTION ChX_1_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChX_1_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
		COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
OPEN(678,FILE='./F/e(T)00._1')
DO I=0,NPD
READ(678,*,END=677)WAPPA(I),APPA(I)
END DO
677 IAP=I-1
OPEN(10,FILE='./F/'//'e.0_1')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
CALL T_CORRECTION(A1,WAPPA,APPA,IAP)  ! <<<*********************** T-Correction ! ************************* !
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChX_1_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChX_1_CO2 =AFL ; RETURN ; END IF
! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChX_1_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_X_2(NCH)
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 10 - 200 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
write(*,*)'BAND_X_2'
RABMA=0.
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'e',NCH) ! Planck function for region 'C'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2> 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        H2O
      	A1=0.
      IF(NCH == 2)THEN
VOL1=ChX_2_CO2(P_LN(J))*CO21
VOL2=ChX_2_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
write(*,*)'***BAND_X_2'
       END

         FUNCTION ChX_2_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChX_2_CO2
DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
!### OPEN(678,FILE='./F/X.1_2')
OPEN(678,FILE='./F/e(T)00._2')
DO I=0,NPD
READ(678,*,END=677)WAPPA(I),APPA(I)
END DO
677 IAP=I-1
 OPEN(10,FILE='./F/'//'e.0_2')  !     HAUS
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
CALL T_CORRECTION(A1,WAPPA,APPA,IAP)  ! <<<***
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChX_2_CO2=ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChX_2_CO2=AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChX_2_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_X_3(NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 10 - 200 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
write(*,*)'BAND_X_3'
! --------------------------------------------------------- *
RABMA=0.
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'e',NCH) ! Planck function for region 'C'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 3)THEN
VOL1=ChX_3_CO2(P_LN(J))*CO21
VOL2=ChX_3_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF

! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
write(*,*)'***BAND_X_3'
       END

         FUNCTION ChX_3_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChX_3_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
OPEN(678,FILE='./F/e(T)00._3')
DO I=0,NPD
READ(678,*,END=677)WAPPA(I),APPA(I)
END DO
677 IAP=I-1
 OPEN(10,FILE='./F/'//'e.0_3')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
  CALL T_CORRECTION(A1,WAPPA,APPA,IAP)  ! <<<***
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChX_3_CO2=ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChX_3_CO2=AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChX_3_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_X_4(NCH)  ! 27Jul.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 500 - 600 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
write(*,*)'BAND_X_4'
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'e',NCH) ! Planck function for region 'C'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 4)THEN
VOL1=ChX_4_CO2(P_LN(J))*CO21
VOL2=ChX_4_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
write(*,*)'*** BAND_X_4'
       END

        FUNCTION ChX_4_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
          REAL*8 P8,ChX_4_CO2
         PARAMETER (NPD=64,JM=64) ! Number of points in internal P-grid (LOG-scale).
         DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
COMMON/Tcor/WH(0:JM),DTdT(0:JM)
         SAVE IBG,NP,A1
            DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
OPEN(678,FILE='./F/e(T)00._4') !APPA11.0c4')
DO I=0,NPD
READ(678,*,END=677)WAPPA(I),APPA(I)
END DO
677 IAP=I-1
OPEN(10,FILE='./F/'//'e.0_4') ! HAUSS
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
 CALL T_CORRECTION(A1,WAPPA,APPA,IAP)  ! <<<*********************** T-Correction ! ************************* !
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChX_4_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChX_4_CO2 =AFL ; RETURN ; END IF
! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChX_4_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_X_5(NCH)  ! 27Jul.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 500 - 600 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
write(*,*)'BAND_X_5'
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'e',NCH) ! Planck function for region 'C'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 5)THEN
VOL1=ChX_5_CO2(P_LN(J))*CO21
VOL2=ChX_5_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
write(*,*)'### BAND_X_5'
       END

         FUNCTION ChX_5_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChX_5_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
OPEN(678,FILE='./F/e(T)00._5')
DO I=0,NPD
READ(678,*,END=677)WAPPA(I),APPA(I)
END DO
677 IAP=I-1
OPEN(10,FILE='./F/'//'e.0_5')    ! HAUS !!!
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
 CALL T_CORRECTION(A1,WAPPA,APPA,IAP)  ! <<<*******
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChX_5_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChX_5_CO2 =AFL ; RETURN ; END IF
! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChX_5_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_X_6(NCH)  ! 27Jul.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 500 - 600 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
 write(*,*)'BAND_X_6'
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'e',NCH) ! Planck function for region 'C'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 6)THEN
VOL1=ChX_6_CO2(P_LN(J))*CO21
VOL2=ChX_6_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChX_6_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChX_6_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1
! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
OPEN(678,FILE='./F/e(T)00._6')
DO I=0,NPD
READ(678,*,END=677)WAPPA(I),APPA(I)
END DO
677 IAP=I-1
 OPEN(10,FILE='./F/'//'e.0_6')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
  CALL T_CORRECTION(A1,WAPPA,APPA,IAP)  ! <<<*******
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChX_6_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChX_6_CO2 =AFL ; RETURN ; END IF
! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChX_6_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_X_7(NCH)  ! 27Jul.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 500 - 600 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
write(*,*)'BAND_X_7'
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'e',NCH) ! Planck function for region 'C'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 7)THEN
VOL1=ChX_7_CO2(P_LN(J))*CO21
VOL2=ChX_7_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
write(*,*)'BAND_X_7 fin'
       END

         FUNCTION ChX_7_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         REAL*8 P8,ChX_7_CO2
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
OPEN(678,FILE='./F/e(T)00._7')
DO I=0,NPD
READ(678,*,END=677)WAPPA(I),APPA(I)
END DO
677 IAP=I-1
OPEN(10,FILE='./F/'//'e.0_7')   !  HAUS !!!
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
  CALL T_CORRECTION(A1,WAPPA,APPA,IAP)  ! <<<*******
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChX_7_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChX_7_CO2 =AFL ; RETURN ; END IF
! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChX_7_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_X_8(NCH)  ! 27Jul.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 500 - 600 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'e',NCH) ! Planck function for region 'C'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 8)THEN
VOL1=ChX_8_CO2(P_LN(J))*CO21
VOL2=ChX_8_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChX_8_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         REAL*8 P8,ChX_8_CO2
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
         COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
OPEN(678,FILE='./F/e(T)11._8')
    DO I=0,NPD
    READ(678,*,END=677)WAPPA(I),APPA(I)
    END DO
    677 IAP=I-1
OPEN(10,FILE='./F/'//'e.1_8')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
 CALL T_CORRECTION(A1,WAPPA,APPA,IAP)  ! <<<*******
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChX_8_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChX_8_CO2 =AFL ; RETURN ; END IF
! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChX_8_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_X_9(NCH)  ! 27Jul.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 500 - 600 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'e',NCH) ! Planck function for region 'C'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 9)THEN
VOL1=ChX_9_CO2(P_LN(J))*CO21
VOL2=ChX_9_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChX_9_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         REAL*8 P8,ChX_9_CO2
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
OPEN(678,FILE='./F/e(T)11._9')
    DO I=0,NPD
 READ(678,*,END=677)WAPPA(I),APPA(I)
    END DO
    677 IAP=I-1
 OPEN(10,FILE='./F/'//'e.1_9')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
 CALL T_CORRECTION(A1,WAPPA,APPA,IAP)  ! <<<*******
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChX_9_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChX_9_CO2 =AFL ; RETURN ; END IF
! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChX_9_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_X10(NCH)  ! 27Jul.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 500 - 600 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'e',NCH) ! Planck function for region 'C'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 10)THEN
VOL1=ChX10_CO2(P_LN(J))*CO21
VOL2=ChX10_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF


! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChX10_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         REAL*8 P8,ChX10_CO2
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
OPEN(10,FILE='./F/'//'e.110')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
!A1=-45.5  ! wwwwwwwwwwwwww -42.5 wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww !
 A1=-44.5  ! wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww !

          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChX10_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChX10_CO2 =AFL ; RETURN ; END IF
! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChX10_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_Y_1(NCH)  !  7 August, 2024.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 760 - 900 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'f',NCH) ! Planck function for region 'f'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 1)THEN
VOL1=ChY_1_CO2(P_LN(J))*CO21
VOL2=ChY_1_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChY_1_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChY_1_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
		COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1
! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
OPEN(10,FILE='./F/'//'f.1_')
          DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChY_1_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChY_1_CO2 =AFL ; RETURN ; END IF
! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChY_1_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_Y_2(NCH)  !  7 August, 2024.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longXave fluxes and cooling rates *
!* using K-distribution technique IN 760 - 900 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
      DO  L=1,JMAX
       PL(L)=PLANCK_oper(T1(L),'f',NCH) ! Planck function for region 'f'
	END DO
! * ---------------  CO2 ---------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
VOL1=ChY_2_CO2(P_LN(J))*CO21
VOL2=ChY_2_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

        FUNCTION ChY_2_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120,JM=64) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChY_2_CO2
DIMENSION PPP(0:NPD),A1(0:NPD),WAPPA(0:NPD),APPA(0:NPD)
COMMON/Tcor/WH(0:JM),DTdT(0:JM)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
 OPEN(10,FILE='./F/'//'f.2_')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChY_2_CO2=ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChY_2_CO2=AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChY_2_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_g_1(NCH)
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 only
!* Using K-distribution technique in 900-1100  cm^-1 REGION *  =>1KD term number 1 ! <===
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
 RABMA=0.
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'g',NCH) ! Planck function for region 'g'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
VOL1=Ch5_1_CO2(P_LN(J))*CO21
VOL2=Ch5_1_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION Ch5_1_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch5_1_CO2,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'g.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch5_1_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch5_1_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch5_1_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

SUBROUTINE BAND_g_2(NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 only
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/ PL(150)
 RABMA=0.
! --------------------------------------------------------- *
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'g',NCH) ! Planck function for region 'g'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
VOL1=Ch5_2_CO2(P_LN(J))*CO21
VOL2=Ch5_2_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1

!* * * * * * * * * * * * * * * *
       END

         FUNCTION Ch5_2_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,Ch5_2_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
                 OPEN (10,FILE='./F/'//'g.2_') ! File with the data for interpolation.

		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch5_2_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch5_2_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch5_2_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_g_3(NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 only
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/ PL(150)
! --------------------------------------------------------- *
 RABMA=0.
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'g',NCH) ! Planck function for region 'g'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
 VOL1=Ch5_3_CO2(P_LN(J))*CO21
VOL2=Ch5_3_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1

!* * * * * * * * * * * * * * * *
       END

         FUNCTION Ch5_3_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,Ch5_3_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
                    OPEN (10,FILE='./F/'//'g.3_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch5_3_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch5_3_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch5_3_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_g_4(NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 only
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/ PL(150)
! --------------------------------------------------------- *
  RABMA=0.
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'g',NCH) ! Planck function for region 'g'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
VOL1=Ch5_4_CO2(P_LN(J))*CO21
VOL2=Ch5_4_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1

!* * * * * * * * * * * * * * * *
       END

         FUNCTION Ch5_4_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,Ch5_4_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
                 OPEN (10,FILE='./F/'//'g.4_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch5_4_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch5_4_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch5_4_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_g_5(NCH)
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 only
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/ PL(150)
! --------------------------------------------------------- *
 RABMA=0.
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'g',NCH) ! Planck function for region 'g'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
VOL1=Ch5_5_CO2(P_LN(J))*CO21
VOL2=Ch5_5_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION Ch5_5_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch5_5_CO2,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'g.5_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch5_5_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch5_5_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch5_5_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_h(NCH)
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 and SO2 only
!* Using K-distribution technique in 1100 -1210  cm^-1 REGION (F) *  =>1KD term number 1 ! <===
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
 RABMA=0.
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'h',NCH) ! Planck function for region 'F'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
VOL1=Ch6_1_CO2(P_LN(J))*CO21
VOL2=Ch6_1_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
           END DO
RABMA(1)=VOL1
! * ------------------------------------------------------------------- *
! * --------------   SO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       SO21=0.; SO22=0.
       IF(ISO2 > 0)THEN
      SO21=RO1(ISO2,J) ;  SO22=RO1(ISO2,J+1)
      END IF
!        SO2
VOL1=Ch6_1_SO2(P_LN(J))*SO21
VOL2=Ch6_1_SO2(P_LN(J+1))*SO22
RABMA(J+1)=RABMA(J+1)+VOL2
           END DO
RABMA(1)=RABMA(1)+VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION Ch6_1_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch6_1_CO2,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'h_c.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch6_1_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch6_1_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch6_1_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

         FUNCTION Ch6_1_SO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch6_1_SO2,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'h_s.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch6_1_SO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch6_1_SO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch6_1_SO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_i(NCH)
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 and SO2 only
!* Using K-distribution technique in 1100 -1210  cm^-1 REGION (F) *  =>1KD term number 1 ! <===
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
 RABMA=0.
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'i',NCH) ! Planck function for region 'F'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
VOL1=Ch7_1_CO2(P_LN(J))*CO21
VOL2=Ch7_1_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
! * --------------   SO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       SO21=0.; SO22=0.
       IF(ISO2 > 0)THEN
      SO21=RO1(ISO2,J) ;  SO22=RO1(ISO2,J+1)
      END IF
!        SO2
VOL1=Ch7_1_SO2(P_LN(J))*SO21
VOL2=Ch7_1_SO2(P_LN(J+1))*SO22
RABMA(J+1)=RABMA(J+1)+VOL2
           END DO
RABMA(1)=RABMA(1)+VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION Ch7_1_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch7_1_CO2,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'i_c.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch7_1_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch7_1_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch7_1_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

         FUNCTION Ch7_1_SO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch7_1_SO2,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'i_s.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch7_1_SO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch7_1_SO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch7_1_SO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_j(NCH)
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 and SO2 only
!* Using K-distribution technique in 1100 -1210  cm^-1 REGION (F) *  =>1KD term number 1 ! <===
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
 RABMA=0.
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'j',NCH) ! Planck function for region 'F'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
VOL1=Ch8_1_CO2(P_LN(J))*CO21
VOL2=Ch8_1_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
! * --------------   SO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       SO21=0.; SO22=0.
       IF(ISO2 > 0)THEN
      SO21=RO1(ISO2,J) ;  SO22=RO1(ISO2,J+1)
      END IF
!        SO2
VOL1=Ch8_1_SO2(P_LN(J))*SO21
VOL2=Ch8_1_SO2(P_LN(J+1))*SO22
RABMA(J+1)=RABMA(J+1)+VOL2
           END DO
RABMA(1)=RABMA(1)+VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION Ch8_1_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch8_1_CO2,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'j_c.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch8_1_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch8_1_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch8_1_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
         FUNCTION Ch8_1_SO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch8_1_SO2,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'j_s.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch8_1_SO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch8_1_SO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch8_1_SO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_k(NCH)
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2, SO2 and H2O only
!* Using K-distribution technique in 1400 -1800  cm^-1 REGION (g) *  =>1KD term number 1 ! <===
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/PL(150)
! --------------------------------------------------------- *
 RABMA=0.
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'k',NCH) ! Planck function for region 'F'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
VOL1=Ch9_1_CO2(P_LN(J))*CO21
VOL2=Ch9_1_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1
! * --------------   SO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       SO21=0.; SO22=0.
       IF(ISO2 > 0)THEN
      SO21=RO1(ISO2,J) ;  SO22=RO1(ISO2,J+1)
      END IF
!        SO2
VOL1=Ch9_1_SO2(P_LN(J))*SO21
VOL2=Ch9_1_SO2(P_LN(J+1))*SO22
RABMA(J+1)=RABMA(J+1)+VOL2
           END DO
RABMA(1)=RABMA(1)+VOL1
! * --------------   H2O ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       H2O1=0.; H2O2=0.
       IF(IH2O > 0)THEN
      H2O1=RO1(IH2O,J) ;  H2O2=RO1(IH2O,J+1)
      END IF
!        H2O
VOL1=Ch9_1_H2O(P_LN(J))*H2O1
VOL2=Ch9_1_H2O(P_LN(J+1))*H2O2
RABMA(J+1)=RABMA(J+1)+VOL2
           END DO
RABMA(1)=RABMA(1)+VOL1
!* * * * * * * * * * * * * * * *
       END

         FUNCTION Ch9_1_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch9_1_CO2,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'k_c.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch9_1_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch9_1_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch9_1_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
         FUNCTION Ch9_1_SO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch9_1_SO2,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'k_s.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch9_1_SO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch9_1_SO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch9_1_SO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
         FUNCTION Ch9_1_H2O(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 Ch9_1_H2O,P8
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'k_h.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; Ch9_1_H2O =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; Ch9_1_H2O =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
Ch9_1_H2O=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_l(NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 only in this version
!* using K-distribution technique IN 1800-2650  cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/ PL(150)
! --------------------------------------------------------- *
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'l',NCH) ! Planck function for region 'C'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 1)THEN
VOL1=ChJ_CO2(P_LN(J))*CO21
VOL2=ChJ_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1

!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChJ_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChJ_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'l.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChJ_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChJ_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChJ_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_m(NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 only in this version
!* using K-distribution technique IN 2650 - 3000 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/ PL(150)
! --------------------------------------------------------- *
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'m',NCH) ! Planck function for region 'C'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 1)THEN
VOL1=ChK_CO2(P_LN(J))*CO21
VOL2=ChK_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1

!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChK_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChK_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'m.1_') ! File with the data for interpolation.
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChK_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChK_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChK_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_n(NCH) ! 9 Feb.,2003.
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 only in this version
!* using K-distribution technique IN 3000 - 4100 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/ PL(150)
! --------------------------------------------------------- *
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'n',NCH) ! Planck function for region 'C'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 1)THEN
VOL1=ChL_CO2(P_LN(J))*CO21
VOL2=ChL_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1

!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChL_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChL_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
          OPEN (10,FILE='./F/'//'n.1_')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChL_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChL_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChL_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

 SUBROUTINE BAND_o(NCH) ! .
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 only in this version
!* using K-distribution technique IN 4100 - 4400 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/ PL(150)
! --------------------------------------------------------- *
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'o',NCH) ! Planck function for region 'C'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 1)THEN
VOL1=ChM_CO2(P_LN(J))*CO21
VOL2=ChM_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1

!* * * * * * * * * * * * * * * *
! * ---------------- H2O --------------------- *
            JMAX_1=JMAX-1
            DO J=JMAX_1,1,-1
          H2O1=0.; H2O2=0.
          IF(IH2O > 0)THEN
         H2O1=RO1(IH2O,J) ;  H2O2=RO1(IH2O,J+1)
         END IF
   !        H2O
         	A1=0.
         IF(NCH == 1)THEN
    VOL1=ChM_H2O(P_LN(J))*H2O1
    VOL2=ChM_H2O(P_LN(J+1))*H2O2
   RABMA(J+1)=RABMA(J+1)+VOL2
    END IF
! * ------------------------------------------------------------------- *
              END DO
    RABMA(1)=RABMA(1)+VOL1

!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChM_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChM_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
       OPEN (10,FILE='./F/'//'o_c.1_')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChM_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChM_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChM_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END

! ***  H2O *** !
           FUNCTION ChM_H2O(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChM_H2O
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
       OPEN (10,FILE='./F/'//'o_h.1_')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChM_H2O =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChM_H2O =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChM_H2O=EXP(ABSOR_LOG)
! --------------------------------------------
          END
 SUBROUTINE BAND_p(NCH) !
USE MOD_IN_GAS_CLOUD
!*************************************************************
!* This program calculates longwave fluxes and cooling rates *  => CO2 only in this version
!* using K-distribution technique IN 4400 - 6000 cm^-1 REGION *
!*************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
	  REAL*4 PL
COMMON/KD/ PL(150)
! --------------------------------------------------------- *
         DO  L=1,JMAX
        PL(L)=PLANCK_oper(T1(L),'p',NCH) ! Planck function for region 'N'
	END DO
! * --------------   CO2 ----------------------- *
         JMAX_1=JMAX-1
         DO J=JMAX_1,1,-1
       CO21=0.; CO22=0.
       IF(ICO2 > 0)THEN
      CO21=RO1(ICO2,J) ;  CO22=RO1(ICO2,J+1)
      END IF
!        CO2
      	A1=0.
      IF(NCH == 1)THEN
VOL1=ChN_CO2(P_LN(J))*CO21
VOL2=ChN_CO2(P_LN(J+1))*CO22
RABMA(J+1)=VOL2
 END IF
! * ------------------------------------------------------------------- *
           END DO
RABMA(1)=VOL1

!* * * * * * * * * * * * * * * *
       END

         FUNCTION ChN_CO2(P8) ! ln(P) [mbar] RO [mol/(cm^2*km].
!*** CROSS SECTIONS ***!
 USE INITIALforKD                          !   (0,1,...,NP) P&T - interpolation.
         PARAMETER (NPD=120) ! Number of points in internal P-grid (LOG-scale).
         REAL*8 P8,ChN_CO2
         DIMENSION PPP(0:NPD),A1(0:NPD)
		 SAVE IBG,NP,A1

! ---------------------------------------------------------------------- !
          DATA IBG/0/
! -----  Initial information reading --------
         IF(IBG.EQ.0)THEN
         IBG=1
         OPEN (10,FILE='./F/'//'p.1_')
		  DO I=0,NPD
          READ (10,*,END=345)PPP(I),A1(I)
          END DO
345 NP=I-1
          CLOSE(10)
          PSL=PPP(0)  ; ASL=EXP(A1(0))
          PFL=PPP(NP) ; AFL=EXP(A1(NP))
  		 END IF
! --------------------------------------------
IF(P8>=PSL)THEN ; ChN_CO2 =ASL ; RETURN ; END IF
IF(P8<=PFL)THEN ; ChN_CO2 =AFL ; RETURN ; END IF

! ------   linear log(P)interpolation ------
DO J=1,NP
IF(PPP(J-1)>=P8.AND.PPP(J)<=P8)EXIT
END DO
CLIN2=(PPP(J-1)-P8)/(PPP(J-1)-PPP(J)) ; CLIN1=1.-CLIN2
ABSOR_LOG=CLIN1*A1(J-1)+CLIN2*A1(J)
ChN_CO2=EXP(ABSOR_LOG)
! --------------------------------------------
          END
























































