! -------------------------------------------------------------------------
!
SUBROUTINE GAMMAV(EPSS,THETA,GAMMAV_VAR)
!
! -------------------------------------------------------------------------
!
!     CODE ORIGINALLY OBTAINED IN MATLAB FROM PULLIAINEN
!
!    VERSION HISTORY:
!      1.0    JPK 5.12.89
!      1.1    JP 17.03.93
!      2.0    MD 1 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!
!     NOTE: ADDITIONAL FINNISH COMMENTS NOT COPIED HERE

IMPLICIT NONE

REAL,INTENT(IN) :: THETA
COMPLEX,INTENT(IN) :: EPSS
REAL,INTENT(OUT) :: GAMMAV_VAR
REAL THETA_RAD,PI,COSTHETA
COMPLEX NELIO

PI=3.14159
THETA_RAD=THETA/180*PI
COSTHETA=COS(THETA_RAD)
NELIO=SQRT(EPSS-SIN(THETA_RAD)**2)
GAMMAV_VAR=(ABS((EPSS*COSTHETA-NELIO)/(EPSS*COSTHETA+NELIO)))**2

END SUBROUTINE GAMMAV


