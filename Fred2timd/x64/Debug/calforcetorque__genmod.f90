        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALFORCETORQUE__genmod
          INTERFACE 
            SUBROUTINE CALFORCETORQUE(WINDDRT,TRSTFORC,AEROTORQ,XRL,XPC,&
     &FORCE)
              REAL(KIND=8), INTENT(IN) :: WINDDRT
              REAL(KIND=8), INTENT(IN) :: TRSTFORC
              REAL(KIND=8), INTENT(IN) :: AEROTORQ
              REAL(KIND=8), INTENT(IN) :: XRL(3)
              REAL(KIND=8), INTENT(IN) :: XPC(3)
              REAL(KIND=8), INTENT(OUT) :: FORCE(6)
            END SUBROUTINE CALFORCETORQUE
          END INTERFACE 
        END MODULE CALFORCETORQUE__genmod
