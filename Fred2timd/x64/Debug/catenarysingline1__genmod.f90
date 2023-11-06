        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:21:20 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CATENARYSINGLINE1__genmod
          INTERFACE 
            SUBROUTINE CATENARYSINGLINE1(XF,ZF,L,EA,W,N,CB,HF,VF,HA,VA, &
     &NOTHREAD)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8) :: XF
              REAL(KIND=8) :: ZF
              REAL(KIND=8), INTENT(IN) :: L
              REAL(KIND=8), INTENT(IN) :: EA
              REAL(KIND=8), INTENT(IN) :: W
              REAL(KIND=8), INTENT(IN) :: CB
              REAL(KIND=8), INTENT(OUT) :: HF
              REAL(KIND=8), INTENT(OUT) :: VF
              REAL(KIND=8), INTENT(OUT) :: HA
              REAL(KIND=8), INTENT(OUT) :: VA
              INTEGER(KIND=4), INTENT(IN) :: NOTHREAD
            END SUBROUTINE CATENARYSINGLINE1
          END INTERFACE 
        END MODULE CATENARYSINGLINE1__genmod
