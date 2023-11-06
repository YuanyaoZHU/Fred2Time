        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE IWHC__genmod
          INTERFACE 
            SUBROUTINE IWHC(XLC,TSP,NSTP,H,HS,TP,BETA,ETA)
              INTEGER(KIND=4) :: NSTP
              REAL(KIND=8) :: XLC(3)
              REAL(KIND=8) :: TSP
              REAL(KIND=8) :: H
              REAL(KIND=8) :: HS
              REAL(KIND=8) :: TP
              REAL(KIND=8) :: BETA
              REAL(KIND=8) :: ETA(0:NSTP)
            END SUBROUTINE IWHC
          END INTERFACE 
        END MODULE IWHC__genmod
