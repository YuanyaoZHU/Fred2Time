! ====================================================================   
      SUBROUTINE Trans1to2( Theta1, Theta2, Theta3, TransMat )
!  1to2: matrix
!  2to1: inverse matrix or transpose matrix

      ! This routine computes the 3x3 transformation matrix, TransMat,
      !   to a coordinate system x (with orthogonal axes x1, x2, x3)
      !   resulting from three rotations (Theta1, Theta2, Theta3) about the
      !   orthogonal axes (X1, X2, X3) of coordinate system X.  All angles
      !   are assummed to be small, as such, the order of rotations does
      !   not matter and Euler angles do not need to be used.  This routine
      !   is used to compute the transformation matrix (TransMat) between
      !   undeflected (X) and deflected (x) coordinate systems.  In matrix
      !   form:
      !      {x1}   [TransMat(Theta1, ] {X1-XG1}       (1to2)
      !      {x2} = [         Theta2, ]*{X2-XG2}
      !      {x3}   [         Theta3 )] {X3-XG3}
      
      !   and the inverse transform:
      !      {X1}   [XG1]+[TransMat(Theta1, ]^T {x1}       (2to1)
      !      {X2} = [XG2]+[         Theta2, ]^T*{x2}
      !      {X3}   [XG3]+[         Theta3 )]^T {x3}      

      ! The transformation matrix, TransMat, is the closest orthonormal
      !   matrix to the nonorthonormal, but skew-symmetric, Bernoulli-Euler
      !   matrix:
      !          [   c5c6   c4s6+s4s5c6   s4s6-c4s5c6 ]
      !      A = [  -c5s6   c4c6-s4s5s6   s4c6+c4s5s6 ]
      !          [    s5       -s4c5          c4c5    ]
      ! where ci=cos(Thetai), si=sin(Thetai),ti=tan(Thetai)
       
      ! NOTE: This routine is a MODIFICATION of SUBROUTINE SmllRotTrans() given in
      !       FAST source file FAST.f90.


       IMPLICIT                        NONE

      ! Passed Variables:

       REAL(8), INTENT(IN )      :: Theta1    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: Theta2    ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: Theta3    ! The large rotation about X3, (rad).
       REAL(8), INTENT(OUT)      :: TransMat (3,3)                                  ! The resulting transformation matrix from X to x, (-).


       ! Local Variables:

       REAL(8)                   :: C4,C5,C6,S4,S5,S6                               

!      ! Compute some intermediate results:
!
       C4 = DCOS(Theta1)
       C5 = DCOS(Theta2)
       C6 = DCOS(Theta3)
       
       S4 = DSIN(Theta1)
       S5 = DSIN(Theta2)
       S6 = DSIN(Theta3)   
!
!
       ! Define the transformation matrix:
!
       TransMat(1,1) = C5*C6
       TransMat(2,2) = C4*C6-S4*S5*S6
       TransMat(3,3) = C4*C5
       TransMat(1,2) = C4*S6+S4*S5*C6
       TransMat(2,1) = -C5*S6
       TransMat(1,3) = S4*S6-C4*S5*C6
       TransMat(3,1) = S5
       TransMat(2,3) = S4*C6+C4*S5*S6
       TransMat(3,2) = -S4*C5

!
      RETURN
      END SUBROUTINE Trans1to2
      
! ====================================================================   
      SUBROUTINE Trans3to4_1( YawAng, TransMat )


       IMPLICIT                        NONE

      ! Passed Variables:

       REAL(8), INTENT(IN )      :: YawAng    ! wind turbine yaw angle, (rad).

       REAL(8), INTENT(OUT)      :: TransMat (3,3)                                  ! The resulting transformation matrix from X to x, (-).


       ! Local Variables:

       REAL(8)                   :: C1,S1                             

!      ! Compute some intermediate results:
!
       C1 = DCOS(YawAng)
       S1 = DSIN(YawAng) 
!
!
       ! Define the transformation matrix:
!
       TransMat = 0.D0
 
       TransMat(1,1) = C1
       TransMat(2,2) = C1
       TransMat(3,3) = 1.D0
       TransMat(1,2) = S1
       TransMat(2,1) = -S1

!
      RETURN
      END SUBROUTINE Trans3to4_1
      
      
 ! ====================================================================   
      SUBROUTINE Trans3to4_2( TiltAng, TransMat )


       IMPLICIT                        NONE

      ! Passed Variables:

       REAL(8), INTENT(IN )      :: TiltAng    ! wind turbine tilt angle, (rad).

       REAL(8), INTENT(OUT)      :: TransMat (3,3)                                  ! The resulting transformation matrix from X to x, (-).


       ! Local Variables:

       REAL(8)                   :: C1,S1                             

!      ! Compute some intermediate results:
!
       C1 = DCOS(TiltAng)
       S1 = DSIN(TiltAng) 
!
!
       ! Define the transformation matrix:
!
       TransMat = 0.D0
       
       TransMat(1,1) = C1
       TransMat(2,2) = 1.D0
       TransMat(3,3) = C1
       TransMat(1,3) = -S1
       TransMat(3,1) = S1

!
      RETURN
      END SUBROUTINE Trans3to4_2     
      
! ====================================================================   
      SUBROUTINE Trans3to4( YawAng, TiltAng, TransMat )
!  3to4: matrix
!  4to3: inverse matrix or transpose matrix


       IMPLICIT                        NONE

      ! Passed Variables:

       REAL(8), INTENT(IN )      :: YawAng    ! wind turbine yaw angle, (rad).
       REAL(8), INTENT(IN )      :: TiltAng    ! wind turbine tilt angle, (rad).

       REAL(8), INTENT(OUT)      ::  TransMat(3,3)                                 ! The resulting transformation matrix from X to x, (-).


       ! Local Variables:

       REAL(8)                   :: TransMat1(3,3), TransMat2(3,3)                          

!      ! Compute some intermediate results:
!
       CALL Trans3to4_1( YawAng, TransMat1 )
       CALL Trans3to4_2( TiltAng, TransMat2 )  
       
       CALL BRMUL(TransMat2,TransMat1,3,3,3,TransMat) 

!
      RETURN
      END SUBROUTINE Trans3to4
      
!
!  *************************************************
!  *         PRODUCT OF TWO REAL MATRIX.           *
!  *************************************************  

      SUBROUTINE BRMUL(A,B,M,N,K,C)
	  IMPLICIT  NONE     

!      C(M,K) IS OUTPUT
!
      INTEGER,INTENT(IN)::M,N,K
      REAL*8,INTENT(IN)::A(1:M,1:N),B(1:N,1:K)
	  REAL*8,INTENT(OUT)::C(1:M,1:K)
	  INTEGER I,J,L

	  DO 50 I=1,M
	  DO 50 J=1,K
	     C(I,J)=0.0
	   DO 10 L=1,N
	     C(I,J)=C(I,J)+A(I,L)*B(L,J)
10	   CONTINUE
50	  CONTINUE

	  RETURN
      END SUBROUTINE BRMUL 
      
      
      
 ! ====================================================================   
      SUBROUTINE Trans4to5( WingAng, TransMat )
!  4to5: matrix
!  5to4: inverse matrix or transpose matrix


       IMPLICIT                        NONE

      ! Passed Variables:

       REAL(8), INTENT(IN )      :: WingAng    ! Wing angle

       REAL(8), INTENT(OUT)      :: TransMat (3,3)                                  ! The resulting transformation matrix from X to x, (-).


       ! Local Variables:

       REAL(8)                   :: C1,S1                             

!      ! Compute some intermediate results:
!
       C1 = DCOS(WingAng)
       S1 = DSIN(WingAng) 
!
!
       ! Define the transformation matrix:
!
       TransMat = 0.D0
       
       TransMat(1,1) = 1.D0
       TransMat(2,2) = C1
       TransMat(3,3) = C1
       TransMat(2,3) = S1
       TransMat(3,2) = -S1

!
      RETURN
      END SUBROUTINE Trans4to5
      
! ====================================================================   
      SUBROUTINE Trans5to6( ConeAng, TransMat )
!  5to6: matrix
!  6to5: inverse matrix or transpose matrix


       IMPLICIT                        NONE

      ! Passed Variables:

       REAL(8), INTENT(IN )      :: ConeAng    ! Cone angle

       REAL(8), INTENT(OUT)      :: TransMat (3,3)                                  ! The resulting transformation matrix from X to x, (-).


       ! Local Variables:

       REAL(8)                   :: C1,S1                             

!      ! Compute some intermediate results:
!
       C1 = DCOS(ConeAng)
       S1 = DSIN(ConeAng) 
!
!
       ! Define the transformation matrix:
!
       TransMat = 0.D0
       
       TransMat(1,1) = C1
       TransMat(2,2) = 1.D0
       TransMat(3,3) = C1
       TransMat(1,3) = -S1
       TransMat(3,1) = S1

!
      RETURN
      END SUBROUTINE Trans5to6    
      
      
! ====================================================================   
      SUBROUTINE Trans1to6( YawAng, TiltAng, WingAng, ConeAng, Theta1, Theta2, Theta3, TransMat )
!  1to6: matrix
!  6to1: inverse matrix or transpose matrix


       IMPLICIT                        NONE

      ! Passed Variables:

       REAL(8), INTENT(IN )      :: YawAng    ! wind turbine yaw angle, (rad).
       REAL(8), INTENT(IN )      :: TiltAng    ! wind turbine tilt angle, (rad).
       REAL(8), INTENT(IN )      :: WingAng    ! Wing angle
       REAL(8), INTENT(IN )      :: ConeAng    ! Cone angle
       REAL(8), INTENT(IN )      :: Theta1    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: Theta2    ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: Theta3    ! The large rotation about X3, (rad).
       
       REAL(8), INTENT(OUT)      ::  TransMat(3,3)            ! The resulting transformation matrix from X to x, (-).


       ! Local Variables:

       REAL(8)                   :: TransMat1(3,3), TransMat2(3,3), TransMat3(3,3),  TransMat4(3,3),  TransMat5(3,3),  TransMat6(3,3)                             

!      ! Compute some intermediate results:
!
       CALL Trans1to2( Theta1, Theta2, Theta3, TransMat1 )
       CALL Trans3to4( YawAng, TiltAng, TransMat2 )
       
       CALL Trans4to5( WingAng, TransMat4 )
       CALL Trans5to6( ConeAng, TransMat6 )
       
       CALL BRMUL(TransMat2,TransMat1,3,3,3,TransMat3) 
       CALL BRMUL(TransMat4,TransMat3,3,3,3,TransMat5) 
       CALL BRMUL(TransMat6,TransMat5,3,3,3,TransMat )        
!
      RETURN
      END SUBROUTINE Trans1to6   
      
! 
! ********************************************* 
! * Inverse of a matric  [A] by               * 
! *           LU decomposition                * 
! * From 'NUMERICAL RECIPES'     pp. 35-37    * 
! ********************************************* 
! 
       SUBROUTINE INVERSE(a,y,n,np)
       INTEGER n,np,indx(n)
       REAL*8 a(np,np),b(np),y(np,np),d

	 CALL ludcmp(a,n,np,indx,d)
	  y=0.0d0
       do 10 i=1, n
	  y(i,i)=1.0d0
10	 continue
!
	 do 100 j=1, n
	  do 20 i=1, n
	   b(i)=y(j,i)
20	  continue
!
        call  lubksb(a,n,np,indx,b)
!
	  do 30 i=1, n
	   y(j,i)=b(i)
30	  continue

100	 continue

      END	 SUBROUTINE INVERSE
! 
! ********************************************* 
! * SOLUTION OF LINEAR EQUATIONS  [A]{X}= {B} * 
! *           LU DECOMPOSITION                * 
! * FROM 'NUMERICAL RECIPES'     pp. 35-37    * 
! ********************************************* 
! 
      SUBROUTINE ludcmp(a,n,np,indx,d)

      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (dAbs(a(i,j)).gt.aamax) aamax=dAbs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dAbs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax) then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n) then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END	 SUBROUTINE ludcmp


      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12     continue
!       
	 do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14     continue
       return
      END	 SUBROUTINE lubksb        
! ====================================================================   
      SUBROUTINE Trans6to1( YawAng, TiltAng, WingAng, ConeAng, Theta1, Theta2, Theta3, TransMat )


       IMPLICIT                        NONE

      ! Passed Variables:

       REAL(8), INTENT(IN )      :: YawAng    ! wind turbine yaw angle, (rad).
       REAL(8), INTENT(IN )      :: TiltAng    ! wind turbine tilt angle, (rad).
       REAL(8), INTENT(IN )      :: WingAng    ! Wing angle
       REAL(8), INTENT(IN )      :: ConeAng    ! Cone angle
       REAL(8), INTENT(IN )      :: Theta1    ! The large rotation about X1, (rad).
       REAL(8), INTENT(IN )      :: Theta2    ! The large rotation about X2, (rad).
       REAL(8), INTENT(IN )      :: Theta3    ! The large rotation about X3, (rad).
       
       REAL(8), INTENT(OUT)      ::  TransMat(3,3)         ! The resulting transformation matrix from X to x, (-).


       ! Local Variables:

       REAL(8)                   :: TransMat1(3,3)
       
!      ! Compute some intermediate results:
!
       CALL Trans1to6( YawAng, TiltAng, WingAng, ConeAng, Theta1, Theta2, Theta3, TransMat1 )
       CALL INVERSE(TransMat1,TransMat,3,3)

!
      RETURN
      END SUBROUTINE Trans6to1
      
      
    