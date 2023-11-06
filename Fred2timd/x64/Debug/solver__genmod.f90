        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVER__genmod
          INTERFACE 
            SUBROUTINE SOLVER(IP,A,B,NA,NB,N,NRHS,NSYS)
              INTEGER(KIND=4) :: NSYS
              INTEGER(KIND=4) :: NB
              INTEGER(KIND=4) :: NA
              INTEGER(KIND=4) :: IP
              REAL(KIND=8) :: A(NA,NA,NSYS)
              REAL(KIND=8) :: B(NA,NB,NSYS)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NRHS
            END SUBROUTINE SOLVER
          END INTERFACE 
        END MODULE SOLVER__genmod
