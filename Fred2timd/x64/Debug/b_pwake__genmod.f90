        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE B_PWAKE__genmod
          INTERFACE 
            SUBROUTINE B_PWAKE(IFHEADING,IFWAKE,RR,X01,Y01,Z01,ZHUB,U0, &
     &CT,R0,DV01)
              LOGICAL(KIND=4), INTENT(IN) :: IFHEADING
              LOGICAL(KIND=4), INTENT(IN) :: IFWAKE
              REAL(KIND=8) :: RR
              REAL(KIND=8) :: X01
              REAL(KIND=8) :: Y01
              REAL(KIND=8) :: Z01
              REAL(KIND=8) :: ZHUB
              REAL(KIND=8) :: U0
              REAL(KIND=8) :: CT
              REAL(KIND=8) :: R0
              REAL(KIND=8) :: DV01
            END SUBROUTINE B_PWAKE
          END INTERFACE 
        END MODULE B_PWAKE__genmod
