! -------------------------------------------------------------------------
!
SUBROUTINE RUFFSOIL(F,MV,T,KSIGMA,GND_SIG,THETA,EPS_TOP,&
RHO, SAND, SILT, CLAY,R_H_MOD,R_V_MOD,NUM)  !ADDED MAY18,2015
!
! -------------------------------------------------------------------------
!
!     CODE ORIGINALLY OBTAINED IN MATLAB FROM PULLIAINEN
!
!     FUNCTION FOR CALCULATING REFLECTIVITIES OF ROUGH, BARE SOILS
!     ACCORDING TO THEORY BY WEGMULLER & MÄTZLER
!     F = FREQUENCY [HZ]
!     MV = VOLUMETRIC MOISTURE [0..1]
!     T = TEMP [K]
!     KSIGMA = NORMALIZED SURFACE SDEV [m]
!     THETA = NADIR ANGLE [DEG]
!     EPS_TOP = EPSILON OF OVERLYING MEDIUM
!     EPS_SOIL = DEFINE EPSILON EXPLICITLY
!
!     Nov24,2015, New variables:
!     RHO = SOIL BULK DENSITY [kg/m^3]
!     SAND = SAND CONTENT [%]
!     SILT = SILT CONTENT [%]
!     CLAY = CLAY CONTENT [%]
!
!    VERSION HISTORY:
!      1.0    ?? ?.?.?
!      2.0    MD 1 APR 05 TRANSLATED TO FORTRAN FROM MATLAB


IMPLICIT NONE
REAL,INTENT(IN) :: F,MV,T,KSIGMA,THETA
REAL,INTENT(IN) :: RHO,SAND,SILT,CLAY
REAL,INTENT(IN) ::  GND_SIG
COMPLEX,INTENT(IN) :: EPS_TOP
REAL,INTENT(OUT) :: R_H_MOD,R_V_MOD
REAL GND_TEMP,GND_EPS,PI,FRESNEL_H,FRESNEL_V,THETA_R
COMPLEX EPS_SOIL,EPS_EFF
Integer :: NUM



PI=3.14159
GND_TEMP=T-273.15 ! K TO C


CALL EPSSOIL2(MV,GND_TEMP,F,RHO,SAND,SILT,CLAY,EPS_SOIL)

EPS_EFF=EPS_SOIL/EPS_TOP

CALL GAMMAH(EPS_EFF,THETA,FRESNEL_H)
CALL GAMMAV(EPS_EFF,THETA,FRESNEL_V)

THETA_R=THETA/180*PI

R_H_MOD=FRESNEL_H*EXP(-KSIGMA**((0.1*COS(THETA_R))**0.5))


IF (THETA<=60) THEN
    R_V_MOD=R_H_MOD*COS(THETA_R)**0.65
ELSEIF(THETA==70)THEN
    R_V_MOD=R_H_MOD*0.621
ELSE


!Revised by Jinmei, to be revised if necessary
PRINT *, 'ERROR: PROPOGATION ANGLE LARGER THAN THE ONE WHICH COULD BE HANDLED BY SOIL ROUGHNESS MODEL'
    R_V_MOD=R_H_MOD*0.621
END IF

!Jinmei added, roughness smaller than 0.0001 mm will be considered smooth
IF(GND_SIG.LE.0.0001e-3) THEN  !The Avg jump of soil was at least 0.0001 m, i.e. 0.1 mm
    R_V_MOD=FRESNEL_V
ENDIF

IF(NUM.eq.2 .and. .false.)Then
Print *,'EPS_SOIL',EPS_SOIL
Print *,'EPS_TOP',EPS_TOP
Print *,'EPS_EFF',EPS_EFF
Endif

END SUBROUTINE RUFFSOIL






! -------------------------------------------------------------------------
!
SUBROUTINE EPSSOIL2(MV0,T,F,RHOS_KG,SAND,SILT,CLAY,EPSS_VAR)
!
! -------------------------------------------------------------------------
!
!     CODE ORIGINALLY OBTAINED IN MATLAB FROM PULLIAINEN
!
!     FUNCTION FOR CALCULATING EPSILON FOR SOIL USING FREQUENCY,
!     TEMPERATURE AND VOLUMETRIC SOIL MOSITURE. USES EPSW.M FOR
!     DIELECTRICITY OF WATER.
!
!     BY J. PULLIAINEN (MOD. BY K. TIGERSTEDT)
!
!     MV - SOIL VOLUMETRIC MOISTURE[0..1]
!     T - SOIL TEMPERATURE [C]
!     F - FREQUENCY[HZ]
!     RHOS_KG [KG/M3] SOIL BULK DENSITY
!     SAND - SAND CONTENT[%]
!     SILT - SILT CONTENT[%]
!     CLAY - CLAY CONTENT[%]
!
!     NOTE: SOME FINNISH COMMENTS WERE NOT COPIED IN ENTIRETY -MD
!   VERSION HISTORY:
!      1.0    JP ?.?.?
!      1.1    KT ?.?.?
!      2.0    MD 1 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!      3.0    JINMEI, REVISED ACCORDING TO THE DOBSON MODEL

IMPLICIT NONE

REAL, INTENT(IN) :: MV0,T,F,RHOS_KG,SAND,SILT,CLAY
COMPLEX,INTENT(OUT) :: EPSS_VAR
REAL RHOS,MV
REAL BETA,ALFA,REI,A,B,C,M,IEI_S,A_P,B_P,C_P,IEI_P,&
DELTA_IEI,INVT,ALF,B_1,B_2,BB,BET_M,BET_DELTA,BET,IEI,SS,EW_R,EW_I
REAL BETA_R, BETA_I, EPSALF_R, EPSS_R, EPSS_I
REAL SIGMAE, TEMP, TK
INTEGER OPT_CONSIDER_ICE
COMPLEX EW,EPSALF

MV=MV0
! July,2015, revised to prevent Inf values
IF(MV0==0.0) MV=0.00001

ALFA=0.65
!BETA=1.09-0.11*SAND+0.18*CLAY
RHOS=RHOS_KG/1000.0


! OPT_CONSIDER_ICE=0 means, you think the MV is the equivalent unfrozen water content
! OPT_CONSIDER_ICE=1 means, you think the water is totaly frozen
! I have a better code with part of water, part of water if necessary
OPT_CONSIDER_ICE=0
IF(T>0 .OR. OPT_CONSIDER_ICE==0)THEN
    TK=T+273.15
    M=F/1E9
    CALL EPSW(M,TK,EW_R,EW_I)
    !PRINT *,'EW_I',EW_I
    !Jinmei added, because according to Dobson, soilwater in soil in not pure water
    IF (M>=1.4) THEN
        SIGMAE=-1.645 +1.9390*RHOS -0.0225622*SAND +0.01594*CLAY
    ELSE
        SIGMAE=0.0467 +0.2204*RHOS -0.004111*SAND +0.006614*CLAY
    ENDIF
    EW_I=EW_I+SIGMAE*(1-RHOS/2.65)/(0.05563132*M*MV)
ELSE
    !ASSUME ALL WATER HAD BEEN FROZEN TO ICE
    REI=3.1884+9.1E-4*T          ! MATZLER AND WEGMULLER 1987
    A=0.0026                     ! IMPURE ICE -5 ASTETTA (MATZLER)
    B=0.00023
    C=0.87
    M=F/1E9
    IEI_S=A/M+B*M**C
    A_P=6E-4                     ! PURE ICE -5 ASTETTA (MATZLER)
    B_P=6.5E-5
    C_P=0.7
    IEI_P=A_P/M+B_P*M**C_P
    DELTA_IEI=IEI_S-IEI_P

    ! HUFFORD 1991
    INVT=300/(T+273)-1
    ALF=(0.00504+0.0062*INVT)*EXP(-22.1*INVT)

    ! (MISHIMA,MATZLER)
    B_1=0.0207
    B_2=1.16E-11
    BB=335
    BET_M=B_1/(T+273)*EXP(BB/(T+273))/(EXP(BB/(T+273))-1)**2+B_2*M**2
    BET_DELTA=EXP(-10.02+0.0364*T)
    BET=BET_M+BET_DELTA
    IEI=ALF/M+BET*M
    SS=10
    IEI=IEI+DELTA_IEI*SS/13

    EW_R=REI
    EW_I=IEI
END IF

EW=CMPLX(EW_R,(-1*EW_I))

!this was coded wrong!
!EPSALF=1+0.65*RHOS+MV**BETA*(EW**ALFA-1)
!EPSS_VAR=EPSALF**(1/ALFA)

BETA_R=1.2748 - 0.00519 * SAND - 0.00152 * CLAY
BETA_I=1.33797 - 0.00603 * SAND - 0.00166 * CLAY

TEMP=RHOS/2.65*(4.7**ALFA-1.0)
EPSALF_R=1+TEMP+MV**BETA_R*(EW_R**ALFA)-MV
EPSS_R=EPSALF_R**(1/ALFA)
IF(M<1.4)THEN
EPSS_R=EPSS_R*1.15-0.68
ENDIF

EPSS_I=MV**(BETA_I/ALFA)*EW_I

EPSS_VAR=CMPLX(EPSS_R,(-1*EPSS_I))

!PRINT *,'EW_R',EW_R,'EW_I',EW_I
!PRINT *,'SIGMAE',SIGMAE
!PRINT *,'BETA_R',BETA_R,'BETA_I',BETA_I

END SUBROUTINE EPSSOIL2