        !COMPILER-GENERATED INTERFACE MODULE: Fri May 30 12:00:56 2014
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
