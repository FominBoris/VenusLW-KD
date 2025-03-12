SUBROUTINE GRID_66V11
USE MOD_IN_GAS_CLOUD
PARAMETER (JM=64,JATM=66)
DIMENSION ZH(0:JM),TH(0:JM),Wv(0:JM),TVv(0:JM),TVh(0:JM) &
,Z_ATM(JATM),T_ATM(JATM),W_ATM(JATM),P_ATM(JATM)                                               
COMMON/Tcor/WH(0:JM),DTdT(0:JM)

! ------  now VH and TH from Venus.11 ------ !
DATA WH/ 11.36924, 11.24064, 11.10933, 10.97523, 10.83816, 10.69797, 10.55463, 10.40800, 10.25774, 10.10361, &
  9.94512,  9.78212,  9.61436,  9.44131,  9.26277,  9.07825,  8.88741,  8.68988,  8.48502,  8.27205, &
  8.05030,  7.81907,  7.57735,  7.32438,  7.05984,  6.78193,  6.48816,  6.17539,  5.83871,  5.47712, &
  5.09705,  4.70841,  4.32086,  3.93960,  3.56413,  3.19215,  2.81985,  2.44327,  2.05744,  1.65647, &
  1.23760,  0.80021,  0.34311, -0.13299, -0.62411, -1.12676, -1.63861, -2.15840, -2.68551, -3.21838, &
  -3.78580, -4.40976, -5.01552, -5.58775, -6.17393, -6.78529, -7.59069, -8.92231,-10.34946,-11.76139, &
  -13.05004,-14.23103,-15.28519,-16.22346,-17.03439/
DATA TH/726.60,712.40,696.70,681.20,665.90,650.70,635.65,620.70,605.20,588.90,572.50,555.90,539.10,522.25,  &
  505.35,488.40,471.65,454.95,438.00,420.75,403.25,385.90,368.65,351.75,335.40,318.75,300.10,279.60, &
  258.80,242.95,234.30,232.70,236.30,240.15,243.10,243.90,242.05,237.70,229.45,219.35,210.00,200.70, &
  191.95,185.10,180.25,176.60,173.55,170.80,168.50,166.75,164.60,161.55,157.50,152.50,147.50,142.50, &
  134.50,123.00,121.00,130.00,138.50,144.00,146.50,143.50,137.00/
!---------------------- 

DATA ZH/ 1.0,  3.0,  5.0,  7.0,  9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0,&
 31.0, 33.0, 35.0, 37.0, 39.0, 41.0, 43.0, 45.0, 47.0, 49.0, 51.0, 53.0, 55.0, 57.0, 59.0,&
 61.0, 63.0, 65.0, 67.0, 69.0, 71.0, 73.0, 75.0, 77.0, 79.0, 81.0, 83.0, 85.0, 87.0, 89.0,&
 91.0, 93.0, 95.0, 97.0, 99.0,101.0,103.0,105.0,107.0,109.0,111.0,114.0,118.0,122.0,126.0,&
130.0,134.0,138.0,142.0,146.0/,&
Z_ATM/0.0,  2.0,  4.0,  6.0,  8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, &
 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, &
 44.0, 46.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 62.0, 64.0, &
 66.0, 68.0, 70.0, 72.0, 74.0, 76.0, 78.0, 80.0, 82.0, 84.0, 86.0, &
 88.0, 90.0, 92.0, 94.0, 96.0, 98.0,100.0,102.0,104.0,106.0,108.0, &
110.0,112.0,116.0,120.0,124.0,128.0,132.0,136.0,140.0,144.0,148.0/

!*** PROSTOY OTBOR 117 -> 66 in ATMOSPHERIC Profile ***!
DO J=1,JATM
ZZZ=Z_ATM(J)
  DO I=1,JMAX
IF(ZZZ==Z(I))EXIT
  END DO
T_ATM(J)=T1(I) ; W_ATM(J)=P_LN(I); P_ATM(J)=P1(I)
END DO
J=J-1
IF(J/=JATM)THEN
WRITE(*,*)'*** Problems in SIMPLE CHOISE ***'  
STOP
END IF

!*** Simple INTERPOLATION for 'W-S(Z)....'  => Wv(0:JM),TVv(0:JM) (initial Wv - grid!)
DO J=0,JM
TVv(J)=(T_ATM(J+1)+T_ATM(J+2))/2.
PPPPP=(P_ATM(J+1)+P_ATM(J+2))/2.
Wv(J)=ALOG(PPPPP)
!### Wv(J)=(W_ATM(J+1)+W_ATM(J+2))/2. ! *** NOT USED ***
END DO
! --- Itog:  Wv(...) eto LogP na STANDARTNOY setke Z kak v 'W-S(Z)...' dlia DANNOY Atmosphere ---!
!  ---     TVv(...)  todge dlia T ---!    

!***  TVv(Wv-grid)  -> TVh(WH-grid)  ***!
TTTT=TH(0);  DTdT(0)=(TVv(0)-TTTT)/TTTT
TVh(0)=TTTT
DO J=1,JM
WW_=WH(J)
  DO JJ=1,JM 
IF(Wv(JJ)<=WW_)EXIT  
  END DO
  IF(JJ>JM)JJ=JM
C2=(Wv(JJ-1)-WW_)/(Wv(JJ-1)-Wv(JJ))  ; C1=1.-C2
TVh(J)=C1*TVv(JJ-1)+C2*TVv(JJ)
 DTdT(J)=(TVh(J)-TH(J))/TH(J)
END DO
! --- Itog:  DTdT(,,,)   eto (TVh-TH)/TH dlia rscheta T-correction (na WH(...) setke) --- !
END

SUBROUTINE GRID_66
USE MOD_IN_GAS_CLOUD
PARAMETER (JM=64,JATM=66)
DIMENSION ZH(0:JM),TH(0:JM),Wv(0:JM),TVv(0:JM),TVh(0:JM) &
,Z_ATM(JATM),T_ATM(JATM),W_ATM(JATM),P_ATM(JATM)                                               
COMMON/Tcor/WH(0:JM),DTdT(0:JM)
DATA WH/11.36930,11.24071,11.10970,10.97532,10.83688,10.69559,10.55084,10.40299,10.25263,10.09848,&
 9.93814, 9.77306, 9.60600, 9.43373, 9.25221, 9.06383, 8.87161, 8.67273, 8.46615, 8.25217,&
 8.03131, 7.80409, 7.57069, 7.33000, 7.07774, 6.81424, 6.54128, 6.24926, 5.93778, 5.60788,&
 5.26294, 4.90181, 4.53285, 4.15757, 3.77301, 3.38464, 2.98845, 2.58047, 2.16415, 1.73896,&
 1.29634, 0.83315, 0.34738,-0.15816,-0.68196,-1.21884,-1.76584,-2.31796,-2.86974,-3.41556,&
 -3.94922,-4.47346,-5.01553,-5.58776,-6.17394,-6.78529,-7.59071,-8.92231,-10.34946,-11.76139,&
 -13.05004,-14.23103,-15.28519,-16.22346,-17.03439/, &
 ZH/ 1.0,  3.0,  5.0,  7.0,  9.0, 11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 29.0,&
 31.0, 33.0, 35.0, 37.0, 39.0, 41.0, 43.0, 45.0, 47.0, 49.0, 51.0, 53.0, 55.0, 57.0, 59.0,&
 61.0, 63.0, 65.0, 67.0, 69.0, 71.0, 73.0, 75.0, 77.0, 79.0, 81.0, 83.0, 85.0, 87.0, 89.0,&
 91.0, 93.0, 95.0, 97.0, 99.0,101.0,103.0,105.0,107.0,109.0,111.0,114.0,118.0,122.0,126.0,&
130.0,134.0,138.0,142.0,146.0/,&
TH/725.00,709.00,693.00,677.00,661.00,645.00,630.00,616.00,602.00,587.00,570.50,553.50,537.00,520.50,&
503.50,486.50,469.50,453.50,438.00,423.00,409.50,397.00,385.00,372.50,357.50,340.50,322.00,301.50,&
282.50,268.50,257.50,248.00,240.00,236.50,235.00,230.50,225.50,220.00,214.00,207.00,198.50,189.00,&
180.00,173.50,169.00,165.50,163.00,162.00,162.50,166.50,167.50,162.50,157.50,152.50,147.50,142.50,&
134.50,123.00,121.00,130.00,138.50,144.00,146.50,143.50,137.00/,&
Z_ATM/0.0,  2.0,  4.0,  6.0,  8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, &
 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0, 42.0, &
 44.0, 46.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 62.0, 64.0, &
 66.0, 68.0, 70.0, 72.0, 74.0, 76.0, 78.0, 80.0, 82.0, 84.0, 86.0, &
 88.0, 90.0, 92.0, 94.0, 96.0, 98.0,100.0,102.0,104.0,106.0,108.0, &
110.0,112.0,116.0,120.0,124.0,128.0,132.0,136.0,140.0,144.0,148.0/

!*** PROSTOY OTBOR 117 -> 66 in ATMOSPHERIC Profile ***!
DO J=1,JATM
ZZZ=Z_ATM(J)
  DO I=1,JMAX
IF(ZZZ==Z(I))EXIT
  END DO
T_ATM(J)=T1(I) ; W_ATM(J)=P_LN(I); P_ATM(J)=P1(I)
END DO
J=J-1
IF(J/=JATM)THEN
WRITE(*,*)'*** Problems in SIMPLE CHOISE ***'  
STOP
END IF

!*** Simple INTERPOLATION for 'W-S(Z)....'  => Wv(0:JM),TVv(0:JM) (initial Wv - grid!)
DO J=0,JM
TVv(J)=(T_ATM(J+1)+T_ATM(J+2))/2.
PPPPP=(P_ATM(J+1)+P_ATM(J+2))/2.
Wv(J)=ALOG(PPPPP)
!### Wv(J)=(W_ATM(J+1)+W_ATM(J+2))/2. ! *** NOT USED ***
END DO
! --- Itog:  Wv(...) eto LogP na STANDARTNOY setke Z kak v 'W-S(Z)...' dlia DANNOY Atmosphere ---!
!  ---     TVv(...)  todge dlia T ---!    

!***  TVv(Wv-grid)  -> TVh(WH-grid)  ***!
TTTT=TH(0);  DTdT(0)=(TVv(0)-TTTT)/TTTT
TVh(0)=TTTT
DO J=1,JM
WW_=WH(J)
  DO JJ=1,JM 
IF(Wv(JJ)<=WW_)EXIT  
  END DO
  IF(JJ>JM)JJ=JM
C2=(Wv(JJ-1)-WW_)/(Wv(JJ-1)-Wv(JJ))  ; C1=1.-C2
TVh(J)=C1*TVv(JJ-1)+C2*TVv(JJ)
 DTdT(J)=(TVh(J)-TH(J))/TH(J)
END DO
! --- Itog:  DTdT(,,,)   eto (TVh-TH)/TH dlia rscheta T-correction (na WH(...) setke) --- !
END
!--------------------------------------------------------------------------!

SUBROUTINE T_CORRECTION(SHV,WAPPA,APPA,IAPPA)
PARAMETER (JM=64)
DIMENSION SHV(0:JM),WAPPA(0:JM),APPA(0:JM)

COMMON/Tcor/WH(0:JM),DTdT(0:JM)
WBEG=WAPPA(0)

! ---- Before Correction ---- ! 
DO J=0,JM
IF(WH(J)<=WBEG)EXIT
END DO

! ---- Correction ---- ! 
JBEG=J
II=-1
DO J=JBEG,JM
IF(II>IAPPA)EXIT 
II=II+1
SHV(J)=SHV(J)+APPA(II)*DTdT(J)
END DO
END

 PROGRAM Ven_LW_KD_MODEL ! 18 Aug.,2023.
USE A_MOD 
USE MOD_IN_GAS_CLOUD
USE M_C
!*************************************************************************
!* This program calculates longwave fluxes and cooling rates          *
!* using K-distribution technique   and MC for CLOUD  treatment     *
!*************************************************************************
REAL*8 VSTART,VFINISH,VEND,VF,FUP,FDO,S1,S2,FLUXUP,FLUXDO,FLUXUP_,FLUXDO_
CHARACTER BANDS*1,BA*1,FIR*25
    PARAMETER (IOUT=598,NB_D=21) !
DIMENSION BANDS(NB_D)
DIMENSION DP(200)    ! < ***
COMMON/KD/PL(150) 
COMMON/W_N/WN4   
ALLOCATABLE FUP(:),FDO(:),XO(:),FLUXUP(:),FLUXDO(:),FLUXUP_(:),FLUXDO_(:),ZI(:)
  CP_CONST=3.963 ; WRITE(*,*)' *** ',CP_CONST,' ->  Venus ***'
 !* ---------------------------------- *
  X_RAND=130000000001.D0
! -------------- SETTINGS ------------ !
 OPEN(47,FILE='CONTROL_KD_LW.INP') ; OPEN(48,FILE='CHECK.DATA')
READ(47,*)NGAME_1 ; NGAME= NGAME_1 ;  WRITE(48,*)' NGAME = ',NGAME
READ(47,190)ATM_PATH ;  WRITE(48,190)ATM_PATH  
 190 FORMAT(A30)
 ALB=0. ! ### READ(47,*)ALB    ; WRITE(48,*)ALB ! <*** Const in this Modification)
 !#### READ(47,*)Tsurf ;  WRITE(48,*)Tsurf 
READ(47,*)NBD
DO J=1,NBD
READ(47,107)BANDS(J)
END DO
107 FORMAT(A1)
WRITE(48,*)NBD,(BANDS(J),J=1,NBD)
WRITE(*,*)NBD,(BANDS(J),J=1,NBD)
  READ(47,567)FIR  ; WRITE(48,567)FIR
567 FORMAT(A25)
CLOSE(47) ; CLOSE(48)

CALL FILES_CREATING   

CALL ATM_PROF_READING 
ALLOCATE (FUP(JMAX),FDO(JMAX),XO(JMAX), &
FLUXUP(JMAX),FLUXDO(JMAX),FLUXUP_(JMAX),FLUXDO_(JMAX),ZI(JMAX),STAT=IERR)
! ZI(JMAX) seems can be omitted!

Tsurf=T1(1)
 CALL DP_Corr(Z,P1,T1,DP,JMAX)
!--------------------------------------!
CALL  GRID_66 ! ### ATTENTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(*,*)' *** CALL  GRID_66 !!! for Haus (00) ***'
! ------------------------------------ !    
	FUP=0.D0 ; FDO=0.D0 
DO N=1,NBD
BA=BANDS(N)
WRITE(*,*)NBD,N,BA

! intervals a, b, c, d, ...
! ------------------------------ a 10-200 -------------------------------------------------- !
  IF(BA=='a') THEN !  10-200  cm-1
FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
WN4=100. 
DO NCH=1,2
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
IF(NCH==1)  CALL BAND_A_1(NCH)
IF(NCH==2)  CALL BAND_A_2(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX)
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO
END DO
  END IF

! ------------------------------ b  200-400 -------------------------------------------------- !
  IF(BA=='b') THEN ! 200-400  cm-1
FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
WN4=300. ! set optical properties at this value
DO NCH=1,2
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
! there are two k-terms here:
IF(NCH==1)  CALL BAND_B_1(NCH) ! input: effective cross-section and profile of CO2 and yields absorption coeffcients in km-1
IF(NCH==2)  CALL BAND_B_2(NCH) 
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) ! radiative transfer solution
! Flux is split into two terms: flux from unscattered photons (clear sky) and scattering part (Monte-Carlo term)
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO ! sum of fluxes for all k-terms
END DO   
  END IF

! ------------------------------ c 400-500 -------------------------------------------------- !
  IF(BA=='c') THEN ! 400-500  cm-1
FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
WN4=450. 
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
NCH=1 ; CALL BAND_c(NCH) 
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO 
  END IF

! ------------------------------ d  500 -600  (1 Channel) -------------------------------------------------- !
  IF(BA=='d') THEN !  500-600
WN4=550. 
NCH=1 
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
 CALL BAND_d(NCH) 
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO 
  END IF

! !###  ------------------------------ r  500 -600  (7 Channels) -------------------------------------------------- !
 !###  IF(BA=='r') THEN !  500-600
!###  WN4=550. 
!###  FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
!###  DO NCH=1,7
!###           FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
!###  IF(NCH==1)  CALL BAND_r_1(NCH) 
!###  IF(NCH==2)  CALL BAND_r_2(NCH) 
!###  IF(NCH==3)  CALL BAND_r_3(NCH) 
!###  IF(NCH==4)  CALL BAND_r_4(NCH) 
!###  IF(NCH==5)  CALL BAND_r_5(NCH) 
!###  IF(NCH==6)  CALL BAND_r_6(NCH) 
!###  IF(NCH==7)  CALL BAND_r_7(NCH) 
!###  CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) 
!###  FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
!###             FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO 
!###  END DO   
 !###   END IF

! ------------------------------ e - (10 CHANNELs, T correction) -------------------------------------------------- !
  IF(BA=='e') THEN !  600-760
WN4=680. 
IPR=0
DO NCH=1,10 
!#### WRITE(*,*)NCH
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
  IF(NCH==1) CALL BAND_X_1(NCH) 
  IF(NCH==2) CALL BAND_X_2(NCH) 
  IF(NCH==3) CALL BAND_X_3(NCH) 
  IF(NCH==4) CALL BAND_X_4(NCH) 
  IF(NCH==5) CALL BAND_X_5(NCH) 
  IF(NCH==6) CALL BAND_X_6(NCH) 
  IF(NCH==7) CALL BAND_X_7(NCH) 
! =======================================  00 ->11
IF(NCH>=8.AND.IPR==0)THEN
WRITE(*,*)'  CALL  GRID_66 for V11  ***'
IPR=1 ; CALL  GRID_66V11
END IF
  IF(NCH==8) CALL BAND_X_8(NCH) 
  IF(NCH==9) CALL BAND_X_9(NCH) 
  IF(NCH==10) CALL BAND_X10(NCH) 
  CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO 
END DO    ! <*** DO NCH=1,10
CALL  GRID_66
WRITE(*,*)'  CALL  GRID_66 !!! again for Haus (00)  ***'
  END IF

! ------------------------------ f   760-900 cm-1 -------------------------------------------------- !
  IF(BA=='f') THEN !  760-900
WN4=830. 
DO NCH=1,2 
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
 IF(NCH==1) CALL BAND_Y_1(NCH) 
 IF(NCH==2) CALL BAND_Y_2(NCH) 
 CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) 
 FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO 
END DO    ! <*** DO NCH=1,2
  END IF

!  ---------------------   g (900 - 1100 cm^-1) --------- !
IF(BA=='g') THEN
WN4=1000.  
DO NCH=1,5
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
IF(NCH==1) CALL BAND_g_1(NCH)
IF(NCH==2)  CALL BAND_g_2(NCH)
IF(NCH==3) CALL BAND_g_3(NCH)
IF(NCH==4)  CALL BAND_g_4(NCH)
IF(NCH==5) CALL BAND_g_5(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO 
END DO    ! <*** DO NCH=1,5
END IF

!  ---------------------   h (1100 - 1210 cm^-1) --------- !
IF(BA=='h') THEN
WN4=1150.0  
NCH=1
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
 CALL BAND_h(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO 
END IF

! ------------------------- i  (1210 - 1300 cm^-1) ------------------------ !
IF(BA=='i') THEN
WN4=1250.0  
NCH=1
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
CALL BAND_i(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO 
END IF

! ------ j 1300 - 1400  -----------------------!
IF(BA=='j') THEN
WN4=1350.0  
NCH=1
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
CALL BAND_j(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO 
END IF

! ------ k 1400 - 1800 ---------------- !
IF(BA=='k') THEN
WN4=1600.0  
NCH=1
         FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
CALL BAND_k(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
           FUP=FUP+FLUXUP ; FDO=FDO+FLUXDO 
END IF

! -------------------------  l  1800-2650 ------------------------------------------------------- !
  IF(BA=='l') THEN  !    ### letter (L but little)
       FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
 WN4=2000. ; NCH=1 ; CALL BAND_l(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) !###,Tsurf) ! 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
			FUP=FUP+FLUXUP
    		                FDO=FDO+FLUXDO
   END IF

! ------------------------- m 2650-3000 --------------------------------- !
  IF(BA=='m') THEN 
FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
 WN4=2800. ; NCH=1 ; CALL BAND_m(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) !###,Tsurf) ! 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
			FUP=FUP+FLUXUP
    		                FDO=FDO+FLUXDO
  END IF

!------------------------- n 3000-4100 ------------------------------------------------------- !
  IF(BA=='n') THEN
  FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
 WN4=3500. ; NCH=1 ; CALL BAND_n(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX)  
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
			FUP=FUP+FLUXUP
    		                FDO=FDO+FLUXDO
  END IF
! ------------------------- o 4100 - 4400 o ------------------------------------------------------- !
  IF(BA=='o') THEN
 FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
WN4=4300. ; NCH=1 ;  CALL BAND_o(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) !###,Tsurf) ! 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
			FUP=FUP+FLUXUP
    		                FDO=FDO+FLUXDO
  END IF
! ------------------------- p 4400- 6000 ------------------------------------------------------- !
  IF(BA=='p') THEN 
FLUXUP=0.D0 ; FLUXDO=0.D0 ; FLUXUP_=0.D0  ; FLUXDO_=0.D0 
 WN4=5000. ; NCH=1 ; CALL BAND_p(NCH)
CALL FLUX_AB_SCAT(FLUXUP,FLUXDO,FLUXUP_,FLUXDO_,NGAME,JMAX) !###,Tsurf) ! 
FLUXUP=FLUXUP+FLUXUP_ ; FLUXDO=FLUXDO+FLUXDO_
			FUP=FUP+FLUXUP
    		                FDO=FDO+FLUXDO
  END IF

END DO ! *** DO N=1,NBD
! = = = = = = = =!


!*			Writing			*
		OPEN(IOUT,FILE=FIR)
			DO J=JMAX,1,-1 !  100,99,98...km
        QT=0.
       IF(J>1) THEN
        S1=FDO(J-1)-FUP(J-1)
        S2=FDO(J)-FUP(J)
        QT=-(S1-S2)/(P1(J)-P1(J-1))*CP_CONST
     QTcor=(S1-S2)/DP(J)*CP_CONST !### /1013.25
	 if(j==jmax) then ; qt=0. ; qtcor=0. ; end if
		XO(J)=QTcor
         END IF
	WRITE(IOUT,167)Z(J),FUP(J),FDO(J),QTcor,QT ! Here ZI = ZJ is used
 167		FORMAT(F7.2,2E13.5,3X,2E14.4)
			END DO
			XO(1)=XO(2)
			CLOSE(IOUT)

!###		OPEN(IOUT,FILE='__'//FIR)
!###!###			DO J=JMAX,1,-1 !  100,99,98...km
!###DO J=2,JMAX !  100,99,98...km
 !###       QT=0.
!###!###        IF(J>1) THEN
!###        S1=FDO(J-1)-FUP(J-1)
 !###       S2=FDO(J)-FUP(J)
!###        QT=-(S1-S2)/(P1(J)-P1(J-1))*CP_CONST
!###     QTcor=(S1-S2)/DP(J)*CP_CONST !### /1013.25
!###	 if(j==jmax) then ; qt=0. ; qtcor=0. ; end if
!###		XO(J)=QTcor
 !###       END IF
!###	WRITE(IOUT,167)Z(J),FUP(J),FDO(J),QTcor,QT ! Here ZI = ZJ is used
!###			END DO
!###			XO(1)=XO(2)
!###CLOSE(IOUT)

END
        FUNCTION PLANCK_oper(TT,H,N) 
! 3.01.2002 - Integrated Planck function for "separate" intervals *
        IMPLICIT REAL*8 (A-H,O-Z)
REAL*4 TT
        CHARACTER NAME*12,Chan(20)*3,H*1,HO*1
    	PARAMETER (NT=8001)
        DIMENSION T(NT),P(NT)
        DATA NO,HO/0,'0'/
 DATA Chan/'.1_','.2_','.3_','.4_','.5_','.6_','.7_','.8_','.9_' &
 ,'.10','.11','.12','.13','.14','.15','.16','.17','.18','.19','.20'/
        IF(N/=NO.OR.HO/=H)THEN
		NO=N ; HO=H
       NAME='./F/PLANCK_'//H
        OPEN(10,FILE=NAME//Chan(N))
        READ(10,*)N_T
	IF(NT.NE.N_T)THEN
	WRITE(*,*)'N_T NT',N_T,NT
	STOP
	END IF
         DO I=1,N_T
         READ(10,*)T(I),P(I)
         END DO
	          CLOSE(10)
          HH=(T(N_T)-T(1))/(N_T-1)
          TS=T(1)
        END IF 
         I=(TT-TS)/HH+1
          C=(TT-T(I))/HH
         PLANCK_oper= (1.D0-C)*P(I)+C*P(I+1)
         END      
SUBROUTINE DP_Corr(Z,P,T,DP,JMAX)
! *** 18.08.2023.  P smothing in Atm. Prof. ***
CHARACTER N_V*2,H_NAM*40,GAS*3
REAL*8 P
PARAMETER (VEN_M=44,GVEN=8.87,R=8.314462)!  with CLOUD levels. 
DIMENSION Z(JMAX),P(JMAX),T(JMAX),DP(200)
OPEN(111,FILE='DP_cor&bad') 
!------------------------------------------------------------!
! ***********  SMOOTHING ********* !
! *** Delta P-smoothing! ***
DO J=2,JMAX
TEMP=(T(J)+T(J-1))/2.
DZ=Z(J)-Z(J-1)
ZNAM=VEN_M*GVEN*DZ
DEL=(R*TEMP)
RATIO=ZNAM/DEL*0.98227
BOL=EXP(-RATIO)
DPCOR=P(J-1)*(1.-BOL)
DPBAD=(P(J-1)-P(J))
DP(J)=DPCOR
WRITE(111,*)Z(J),DPCOR,DPBAD
END DO
CLOSE(111)
END



