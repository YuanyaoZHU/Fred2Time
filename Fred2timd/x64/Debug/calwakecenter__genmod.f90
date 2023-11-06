        !COMPILER-GENERATED INTERFACE MODULE: Fri Mar 17 20:19:43 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CALWAKECENTER__genmod
          INTERFACE 
            SUBROUTINE CALWAKECENTER(WINDDRT,X1,Y1,X2,Y2,X,Y)
              REAL(KIND=8), INTENT(IN) :: WINDDRT
              REAL(KIND=8), INTENT(IN) :: X1
              REAL(KIND=8), INTENT(IN) :: Y1
              REAL(KIND=8), INTENT(IN) :: X2
              REAL(KIND=8), INTENT(IN) :: Y2
              REAL(KIND=8), INTENT(OUT) :: X
              REAL(KIND=8), INTENT(OUT) :: Y
            END SUBROUTINE CALWAKECENTER
          END INTERFACE 
        END MODULE CALWAKECENTER__genmod
