        !COMPILER-GENERATED INTERFACE MODULE: Fri May 30 12:00:54 2014
        MODULE CONVRETAR__genmod
          INTERFACE 
            SUBROUTINE CONVRETAR(L,TAG,TSP,NSTP,RFK,DX,FTAR,RK)
              INTEGER(KIND=4) :: NSTP
              INTEGER(KIND=4) :: L
              INTEGER(KIND=4) :: TAG
              REAL(KIND=8) :: TSP
              REAL(KIND=8) :: RFK(6,6,0:NSTP)
              REAL(KIND=8) :: DX(6,0:NSTP)
              REAL(KIND=8) :: FTAR(6)
              REAL(KIND=8) :: RK(6,4,0:NSTP)
            END SUBROUTINE CONVRETAR
          END INTERFACE 
        END MODULE CONVRETAR__genmod
