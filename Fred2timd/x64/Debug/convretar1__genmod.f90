        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CONVRETAR1__genmod
          INTERFACE 
            SUBROUTINE CONVRETAR1(L,TAG,TSP,NSTP,RFK,DX,FTAR,RK)
              INTEGER(KIND=4) :: NSTP
              INTEGER(KIND=4) :: L
              INTEGER(KIND=4) :: TAG
              REAL(KIND=8) :: TSP
              REAL(KIND=8) :: RFK(6,6,0:2*NSTP)
              REAL(KIND=8) :: DX(6,0:NSTP)
              REAL(KIND=8) :: FTAR(6)
              REAL(KIND=8) :: RK(6,4,0:NSTP)
            END SUBROUTINE CONVRETAR1
          END INTERFACE 
        END MODULE CONVRETAR1__genmod
