
!   N 自变量个数
!   M 单形定点数
SUBROUTINE MATRIX_SIM_NON
    USE GREENMOD
	USE CUMOD
	
    INTEGER:: NSIK
    EXTERNAL:: FS
    
    NL=NTPN
    M=NBPOINT
    ALLOCATE(SP(NL),SPN(NL),SPX(NL),SPY(NL),SPZ(NL),SPXX(NL),SPXY(NL),SPXZ(NL),SPYY(NL),SPYZ(NL))
    ALLOCATE(VERR(NL),MAL(NL,NL))	
   
    K1=0
    DO IP=1,NL

    IF(IP+K1.LE.M) VERR(IP+K1)=U*VECN(IP+K1,1)
    IF(IP+K1.GT.M) VERR(IP+K1)=0.

    END DO

    DO I=1,M
	    DO J=1,NL
            IF(I.EQ.J) THEN
		        MAL(I,J)=4.*PI*SE(I) !+SUMHIJ(I)
                !MAL(I,J)=2.*PI
            ELSE
                MAL(I,J)=SH(I,J)
            END IF
        END DO
    END DO

    DO I=M+1,NL
	    DO J=1,NL
                MAL(I,J)=U**2/G*SXX(I,J)+SZ(I,J)
        END DO
    END DO

    CALL MATRIX_PREGMRES

    SID=0.005
    SIU=1.5
    SIV=0.5
    EPS=1.E-9
    
    WRITE(*,*)"CALL 单纯形算法"

    CALL SIMPLE_NON(NL,NL+1,QG,SID,SIU,SIV,EPS,NSIK,FS)

    WRITE(*,*)"NSIK          "

END SUBROUTINE    

SUBROUTINE SIMPLE_NON(N,M,X0,D,U,V,EPS,K,FS)

	DIMENSION X(N),F(M),XX(N,M),XT(N),XF(N),X0(N),VERR0(N),DVERP(N),DVERQ(N)
	!DOUBLE PRECISION 
    REAL:: X,Z,F,XX,XT,XF,FR,FL,FG,FT,FF,FS
	INTEGER R,G
	K=1

    D=0.
    DO I=1,N-1
        D=D+ABS(X0(I+1)-X0(I))
    END DO
    D=D/(N-1)/1000.

	!FR=SQRT(1.0D0*M)
	FR=SQRT(1.0*M)
	FL=D*(FR-1.0)/(1.414*N)          !Q
	FG=D*(FR+N-1.0)/(1.414*N)        !P

	DO 10 I=2,M
	DO 10 J=1,N
10	XX(J,I)=X0(J)+FL

    WRITE(*,*)"POINT     1"

	DO 20 I=2,M
20	XX(I-1,I)=X0(I-1)+FG

    WRITE(*,*)"POINT     2"

    GOTO 716
	DO 41 I=1,2
	  DO 31 J=1,N
31	  X(J)=XX(J,I)
	  F(I)=FS(N,X,VERR0)
41	CONTINUE

	DO 42 I=3,M
	  F(I)=0.
      CALL DFS(N,DFL,DFG,1,I-1,DVERP,DVERQ)
      DO J=1,N
	    IF(J.EQ.1) THEN
            F(I)=F(I)+ABS(VERR0(J)-DVERP(J)+DVERQ(J))
        ELSE IF(J.EQ.I-1) THEN
            F(I)=F(I)+ABS(VERR0(J)-DVERQ(J)+DVERP(J))
        ELSE
            F(I)=F(I)+ABS(VERR0(J))
        END IF
      END DO
      
    !WRITE(*,*)VERR0(I-1)
42	CONTINUE
    716 CONTINUE

    !GOTO 616
!   ORIGIN
	DO 40 I=1,M
	  DO 30 J=1,N
30	  X(J)=XX(J,I)
	  F(I)=FS(N,X,VERR0)
40	CONTINUE
!   END ORIGIN
    616CONTINUE

    WRITE(*,*)"START LOOP    "

50	FR=F(1)
	FL=F(1)
	R=1
	L=1

    WRITE(*,*)K

	DO 60 I=2,M
	  IF (F(I).GT.FR) THEN
	    R=I
	    FR=F(I)
	  END IF
	  IF (F(I).LT.FL) THEN
	    L=I
	    FL=F(I)
	  END IF
60	CONTINUE

	G=1
	FG=F(1)
	J=1
	IF (R.EQ.1) THEN
	  G=2
	  FG=F(2)
	  J=2
	END IF

	DO 70 I=J+1,M
	  IF ((I.NE.R).AND.(F(I).GT.FG)) THEN
	    G=I
	    FG=F(I)
	  END IF
70	CONTINUE

	DO 90 J=1,N
	  XF(J)=0.0
	  DO 80 I=1,M
	    IF (I.NE.R) XF(J)=XF(J)+XX(J,I)/N
80	  CONTINUE
	  XT(J)=2.0*XF(J)-XX(J,R)
90	CONTINUE
	FT=FS(N,XT,VERR0)

	IF (FT.LT.F(L)) THEN
	  DO 100 J=1,N
100	  XF(J)=(1.0+U)*XT(J)-U*XF(J)
	  FF=FS(N,XF,VERR0)
	  IF (FF.LT.F(L)) THEN
	    DO 110 J=1,N
110	    XX(J,R)=XF(J)
	    F(R)=FF
	  ELSE
	    DO 120 J=1,N
120	    XX(J,R)=XT(J)
	    F(R)=FT
	  END IF
	ELSE IF (FT.LE.F(G)) THEN
	  DO 130 J=1,N
130	  XX(J,R)=XT(J)
	  F(R)=FT
	ELSE
	  IF (FT.LE.F(R)) THEN
	    DO 140 J=1,N
140	    XX(J,R)=XT(J)
	    F(R)=FT
	  END IF
	  DO 150 J=1,N
150	  XF(J)=V*XX(J,R)+(1.0-V)*XF(J)
	  FF=FS(N,XF,VERR0)
	  IF (FF.GT.F(R)) THEN
        WRITE(*,*)"PROBLEM  "
        
        !GOTO 170
!ORIGIN 
	    DO 170 I=1,M
	      DO 160 J=1,N
	        XX(J,I)=(XX(J,I)+XX(J,L))/2.0
	        X(J)=XX(J,I)
160	      CONTINUE
	      F(I)=FS(N,X,VERR0)
170	    CONTINUE
!ENE ORIGIN

	  ELSE
	    DO 180 J=1,N
180	    XX(J,R)=XF(J)
	    F(R)=FF
	  END IF
	END IF
	FF=0.0
	FT=0.0
	DO 190 I=1,M
	  FF=FF+F(I)/M
	  FT=FT+F(I)*F(I)
190	CONTINUE
	FT=(FT-M*FF*FF)/N

    WRITE(*,*)"EPS   ",K,FT

	IF (FT.GE.EPS) THEN
	  K=K+1
	  IF (K.LT.2001) GOTO 50
	END IF
	DO 210 J=1,N
	  X(J)=0.0
	  DO 200 I=1,M
200	  X(J)=X(J)+XX(J,I)/M
210	CONTINUE
	Z=FS(N,X)
    WRITE(*,*)"LAST EPS     ",Z

	RETURN
	END

FUNCTION FS(NL,SSIGM,SVER)     !计算线性问题
!**************************************************************
!CGENERATE THE REAL MATRIX OF THE INTEGRATION EQUATION
!**************************************************************
	USE GREENMOD
	USE CUMOD

    DIMENSION SSIGM(NL),SVER(NL)
    INTEGER:: M,NL


    M=NBPOINT
    !SIGM=0.
    !SIGM0=0.
    
    !GOTO 909

    !DO ITER=1,10
 

    SPN=0.;SP=0.
    DO IP=1,M
        DO J=1,NL
            SP(IP)=SP(IP)+SA(IP,J)*SSIGM(J)    !求解I点速度势沿X向的导数（自由面）
            IF(IP.NE.J) THEN
                SPN(IP)=SPN(IP)+SH(IP,J)*SSIGM(J)    !求解I点速度势沿法向导数（物面）
            ELSE
                SPN(IP)=SPN(IP)+4.*PI*SE(IP)*SSIGM(J)    !求解I点速度势沿法向导数（物面）
            END IF
        END DO
    END DO

    SPX=0.;SPY=0.;SPZ=0.;SPXX=0.;SPXY=0.;SPXZ=0.;SPYY=0.;SPYZ=0.
    DO IP=M+1,NL
        DO J=1,NL
            SP(IP)=SP(IP)+SA(IP,J)*SSIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPX(IP)=SPX(IP)+SX(IP,J)*SSIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPY(IP)=SPY(IP)+SY(IP,J)*SSIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPZ(IP)=SPZ(IP)+SZ(IP,J)*SSIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPXX(IP)=SPXX(IP)+SXX(IP,J)*SSIGM(J)
            SPXY(IP)=SPXY(IP)+SXY(IP,J)*SSIGM(J)
            SPXZ(IP)=SPXZ(IP)+SXZ(IP,J)*SSIGM(J)
            SPYY(IP)=SPYY(IP)+SYY(IP,J)*SSIGM(J)
            SPYZ(IP)=SPYZ(IP)+SYZ(IP,J)*SSIGM(J)
        END DO
    END DO

    DO IP=1,M
        VERR(IP)=(U*VECN(IP,1)-SPN(IP))
        !VERR(IP)=-U*VECN(IP,1)
    END DO
    
    DO IP=M+1,NL
        VERR(IP)=1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U**2*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPY(IP)*SPXY(IP)-2.*SPX(IP)*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U**2*SPY(IP)*SPXY(IP)-2.*SPY(IP)*SPX(IP)*SPXY(IP)-2.*SPY(IP)*SPY(IP)*SPYY(IP)-2.*SPY(IP)*SPZ(IP)*SPYZ(IP))&
        +SPZ(IP)

        !VERR(IP)=(1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        !-1./2./G*(2.*U**2*SPX(IP)*SPXX(IP))&
        !-1./2./G*(2.*U**2*SPY(IP)*SPXY(IP))&
        !+SPZ(IP))     
        !VERR(IP)=0.
    END DO

    FS=0.
    DO IP=1,NL
        FS=FS+ABS(VERR(IP))
    END DO
    !FS=SQRT(FS)
    SVER=VERR

	RETURN
END function

SUBROUTINE DFS(NL,DFL,DFG,MIP,MIQ,DVERP,DVERQ)     !计算线性问题
	!FL=D*(FR-1.0)/(1.414*N)          !Q
	!FG=D*(FR+N-1.0)/(1.414*N)        !P

	USE GREENMOD
	USE CUMOD

    INTEGER:: M,MIP,MIQ
    DIMENSION:: DVERP(NL),DVERQ(NL)

    M=NBPOINT
    NL=NTPN
    !SIGM=0.
    !SIGM0=0.
    
    !GOTO 909

    !DO ITER=1,10
    

    DO IP=1,M
        IF(IP.NE.MIP) THEN
            DVERP(IP)=U*VECN(IP,1)-SH(IP,MIP)*DFG    !求解I点速度势沿法向导数（物面）
        ELSE
            DVERP(IP)=U*VECN(IP,1)-4.*PI*SE(IP)*DFG    !求解I点速度势沿法向导数（物面）
        END IF
    
        IF(IP.NE.MIQ) THEN
            DVERQ(IP)=U*VECN(IP,1)-SH(IP,MIQ)*DFG   !求解I点速度势沿法向导数（物面）
        ELSE
            DVERQ(IP)=U*VECN(IP,1)-4.*PI*SE(IP)*DFL    !求解I点速度势沿法向导数（物面）
        END IF
    END DO

    DO IP=M+1,NL
            SPX(IP)=SX(IP,MIP)*DFG    !求解I点速度势沿X向的导数（自由面）
            SPY(IP)=SY(IP,MIP)*DFG    !求解I点速度势沿X向的导数（自由面）
            SPZ(IP)=SZ(IP,MIP)*DFG    !求解I点速度势沿X向的导数（自由面）
            SPXX(IP)=SXX(IP,MIP)*DFG
            SPXY(IP)=SXY(IP,MIP)*DFG
            SPXZ(IP)=SXZ(IP,MIP)*DFG
            SPYY(IP)=SYY(IP,MIP)*DFG
            SPYZ(IP)=SYZ(IP,MIP)*DFG
    END DO

    DO IP=M+1,NL
        DVERP(IP)=1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U**2*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPY(IP)*SPXY(IP)-2.*SPX(IP)*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U**2*SPY(IP)*SPXY(IP)-2.*SPY(IP)*SPX(IP)*SPXY(IP)-2.*SPY(IP)*SPY(IP)*SPYY(IP)-2.*SPY(IP)*SPZ(IP)*SPYZ(IP))&
        +SPZ(IP)
    END DO

!QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
    DO IP=M+1,NL
            SPX(IP)=SX(IP,MIQ)*DFL    !求解I点速度势沿X向的导数（自由面）
            SPY(IP)=SY(IP,MIQ)*DFL    !求解I点速度势沿X向的导数（自由面）
            SPZ(IP)=SZ(IP,MIQ)*DFL    !求解I点速度势沿X向的导数（自由面）
            SPXX(IP)=SXX(IP,MIQ)*DFL
            SPXY(IP)=SXY(IP,MIQ)*DFL
            SPXZ(IP)=SXZ(IP,MIQ)*DFL
            SPYY(IP)=SYY(IP,MIQ)*DFL
            SPYZ(IP)=SYZ(IP,MIQ)*DFL
    END DO

    DO IP=M+1,NL
        DVERQ(IP)=1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U**2*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPY(IP)*SPXY(IP)-2.*SPX(IP)*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U**2*SPY(IP)*SPXY(IP)-2.*SPY(IP)*SPX(IP)*SPXY(IP)-2.*SPY(IP)*SPY(IP)*SPYY(IP)-2.*SPY(IP)*SPZ(IP)*SPYZ(IP))&
        +SPZ(IP)
    END DO


	RETURN
END SUBROUTINE