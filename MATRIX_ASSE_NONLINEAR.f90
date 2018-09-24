SUBROUTINE MATRIX_NON     !计算线性问题
!**************************************************************
!CGENERATE THE REAL MATRIX OF THE INTEGRATION EQUATION
!**************************************************************
	USE GREENMOD
	USE CUMOD

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2
    REAL,ALLOCATABLE:: DL(:),NYI(:,:),NRI(:,:),NRT(:,:),verr0(:)
    REAL:: MGX_X,MGX_Y,MGX_Z,MGY_Y,MGY_Z,MGX_XX,MGY_XY,MGZ_XZ,MGX_XY,MGY_YY,MGZ_YZ,MGX_XZ,MGY_YZ
    REAL:: NORM1,NORM2

    ALLOCATE(DL(NBPOINT))

    !CALL PADCOEF
    !CALL DIFFCORD_DECOMP   ! 计算坐标对计算平面导数
    !CALL DIFFCORD_MATRIX_DECOMP ! 计算差分算子


!采用bin文件格式可以减小文件大小，加快读写速度 
    open(16,file='coefuse.bin',form='binary',access='sequential',&
          status='old')
    rewind(16)
	WRITE(*,*)
	WRITE(*,*)'SUB MATRIX.........'
	!ALLOCATE(SH(NTPN,NTPN),SX(NTPN,NTPN),SZ(NTPN,NTPN),SA(NTPN,NTPN),SA1(NTPN,NTPN),SXX(NTPN,NTPN))
	
!读入对角线元素        
    !IF(NIT.EQ.1) THEN
    !OPEN(115,FILE='TESTMATRIXSUM.DAT')
    !DO I=1,NBPOINT    
    !    READ(115,*)J,DL(I)
    !END DO

    !CLOSE(115)
    !END IF

	DO 10 I=1,NTPN
		DO 10 J=1,NTPN
		!READ(16)SA(I,J),SX(I,J),SH(I,J),SZ(I,J),SXX(I,J)   !SA1船体表面诱导系数，SA自由表面诱导系数 
!CREAD(16,*)SH(I,J),SA(I,J)
10	CONTINUE

	NL=NTPN
    M=NBPOINT
    L=NFPOINT

!  NTPN+NRPOINT IN SINGLE CPU
   IF(MOD(NL,NTHREAD).NE.0) THEN
     DO I=1,NTHREAD-1
       IF(MOD(NL,NTHREAD).EQ.I) THEN
         DO J=1,NTHREAD-1
            NFS(J)=(NL-I)/NTHREAD
         END DO
         NFS(NTHREAD)=NFS(NTHREAD-1)+I
         EXIT
       END IF
     END DO
   ELSE 
     DO I=1,NTHREAD
       NFS(I)=(NL)/NTHREAD
     END DO
   END IF

	WRITE(*,*)'    NL=     ',NL,'      NR=     ',NR

    GOTO 20

	!ALLOCATE(MAR(L+NL,NL),VER(L+NL),VERR(L+NL))	
    MAR=0.
    VER=0.
    VERR=0.

	ALLOCATE(MAL(NL,NL),MAL1(NL,NL),VERR(NL),MAR(NL,NL))	
    MAL=0.
    MAR=0.

!产生方程右边的向量		
	DO I=1,L              !第一部分B1
		VERR(I)=0.0
    END DO

    DO I=L+1,L+M           !第二部分B2
		VERR(I)=U*VECN(I-L,1)
    END DO

    DO I=L+M+1,L+NL
		VERR(I)=0.0   !第三部分B3
    END DO	

!ASSIGN THE INITIAL VALUE OF THE RIGHT VECTOR
	OPEN(2,FILE='VERR.DAT')
    
!形成左边的矩阵	
	ALLOCATE(MAL(L+NL,L+NL),MAL1(L+NL,L+NL))	

!形成A11	
    DO I=1,L
	    DO J=1,L
		IF (I.EQ.J)THEN
			MAL(I,J)=-G !*4.*PI
		ELSE 
			MAL(I,J)=0.
		END IF
        END DO
    END DO

!形成A12和A13	
    DO I=1,L
        DO J=L+1,L+NL
            DO K=M+1,M+L
                MAL(I,J)=MAL(I,J)+DMX(I+M,K)*SA(K,J-L)*U
            END DO
        END DO
    END DO

    DO I=1,L
	    DO J=L+1,L+NL
		    !MAL(I,J)=SX(I+M,J-L)*U
        END DO
    END DO

!形成A21	
    DO I=L+1,L+M
	    DO J=1,L
		    MAL(I,J)=0.
        END DO
    END DO

!形成A22和A23	
    DO I=L+1,L+M
	    DO J=L+1,L+NL
            IF(I-L.EQ.J-L) THEN
		        MAL(I,J)=SE(I-L)*4*PI !+SUMHIJ(I-L)
            ELSE
                MAL(I,J)=SH(I-L,J-L)
            END IF
        END DO
    
    END DO

!形成A31	
    DO I=L+M+1,NL+L
	    DO J=1,L
		    !MAL(I,J)=AEX(I-L,J+M)*U
            MAL(I,J)=DMX(I-L,J+M)*U
        END DO
    END DO

!形成A32和A33	
    DO I=L+M+1,NL+L
	    DO J=1+L,NL+L
            !IF(I-L.EQ.J-L) THEN
		    !    MAL(I,J)=4*PI*SE(I-L) !+SZ(I-L,J-L)
            !ELSE
		        MAL(I,J)=SZ(I-L,J-L)
            !END IF
        END DO
    END DO

    DO I=1,NL !+L
        !DO J=1,NL !+L
            WRITE(2,*)I,I,MAL(I,I)
            !IF(ABS(DMX(I,J)).GE.1.E-6) WRITE(2,*)I,J,DMX(I,J)
        !END DO
    END DO

    DO I=M+1,NL
        DO J=M+1,NL
            !IF(ABS(AEX(I,J)).GE.1.E-6) WRITE(2,*)I-M,J-M,AEX(I,J)
        END DO
    END DO


    !goto 30
! 测试程序段
        
    DEALLOCATE(MAL,VERR)

	20 ALLOCATE(MAR(NL+500,NL+500),VER(NL+500),MAL1(NL+500,NL+500))
    ALLOCATE(MAL(NL+500,NL+500),VERR(NL+500))	
    ALLOCATE(NYI(NL+500,1),NRI(NL+500,1),NRT(1,NL+500),VERR0(NL+500))
    MAL=0.

!START PARALLEL
    CALL OMP_SET_NUM_THREADS(NTHREAD)

!$omp parallel private(CP,k1,k2,IP,J,K)
    
    CP=OMP_GET_THREAD_NUM()+1

    K1=SUM(NFS(1:CP))-NFS(CP) 
    K2=SUM(NFS(1:CP))

    DO IP=1,NFS(CP)
    
	!DO I=1,M             !第一部分B1
	!	VERR(I)=U*VECN(I,1)
    !END DO
    IF(IP+K1.LE.M) VERR(IP+K1)=U*VECN(IP+K1,1)

    !DO I=M+1,NL           !第二部分B2
	!	VERR(I)=0.
    !END DO
    IF(IP+K1.GT.M) VERR(IP+K1)=0.

    !DO I=1,M
	!    DO J=1,NL
    !        IF(I.EQ.J) THEN
	!	        MAL(I,J)=4.*PI*SE(I) !+SUMHIJ(I)
    !            !MAL(I,J)=2.*PI
    !        ELSE
    !            MAL(I,J)=SH(I,J)
    !        END IF
    !    END DO
    !END DO
    IF(IP+K1.LE.M) THEN
        DO J=1,NL
            IF(IP+K1.EQ.J) THEN
            !IF(IP+K1.EQ.J.AND.IDS(IP+K1).NE.1) THEN
		        MAL(IP+K1,J)=4.*PI*SE(IP+K1) !+SUMHIJ(IP+K1)
                !MAL(IP+K1,J)=2.*PI
            ELSE
                MAL(IP+K1,J)=SH(IP+K1,J)
            END IF
        END DO    
    END IF

    !MAL1=0
    !MAR=0.

    !DO I=M+1,NL
	    !DO J=1,NL
            !MAR(I,J)=0.
            !DO K=M+1,NL
                !MAL1(I,J)=MAL1(I,J)+SX(I,K)*APX(K,J)
                !MAL1(I,J)=MAL1(I,J)+SX(I,K)*DMX(K,J)
                
                !MAR(I,J)=MAR(I,J)+DMX(I,K)*SA(K,J)
            !END DO

            !MAL1(I,J)=0.
            !DO K=M+1,NL
                !MAL1(I,J)=MAL1(I,J)+SX(I,K)*APX(K,J)
                !MAL1(I,J)=MAL1(I,J)+SX(I,K)*DMX(K,J)
                !MAL1(I,J)=MAL1(I,J)+DMX(I,K)*SX(K,J)
                
                !MAL1(I,J)=MAL1(I,J)+DMX(I,K)*SX(K,J)
                
                !MAL1(I,J)=MAL1(I,J)+DMX(I,K)*MAR(K,J)
            !END DO

            !IF(I.EQ.J) THEN
		    !    MAL(I,J)=U**2/G*MAL1(I,J)+4.*PI*SE(I) !+SUMHIJ(I-L)
                !MAL(I,J)=U**2/G*SXX(I,J)+4.*PI*SE(I) !+SUMHIJ(I-L)
            !ELSE
                !MAL(I,J)=U**2/G*MAL1(I,J)+SZ(I,J)
                
                !MAL(I,J)=U**2/G*SXX(I,J)+SZ(I,J)
            !END IF
        
        !END DO
    !END DO
    IF(IP+K1.GT.M) THEN
	    DO J=1,NL
            !MAR(IP+K1,J)=0.
            DO K=M+1,NL
                !MAL1(I,J)=MAL1(I,J)+SX(I,K)*APX(K,J)
                !MAL1(I,J)=MAL1(I,J)+SX(I,K)*DMX(K,J)
                !MAR(IP+K1,J)=MAR(IP+K1,J)+DMX(IP+K1,K)*SA(K,J)
            END DO

            !MAL1(IP+K1,J)=0.
            DO K=M+1,NL
                !IF (ABS(DMX(IP+K1,K)).GE.1.E-6) THEN
                    !MAL1(IP+K1,J)=MAL1(IP+K1,J)+DMX(IP+K1,K)*SX(K,J)
                !END IF
            END DO

            !IF(IP+K1.NE.J) THEN
                !MAL(IP+K1,J)=U**2/G*MAL1(IP+K1,J)+SZ(IP+K1,J)
                !MAL(IP+K1,J)=U**2/G*SY(IP+K1,J)+SZ(IP+K1,J)
                MAL(IP+K1,J)=U**2/G*SXX(IP+K1,J)+SZ(IP+K1,J)
            !ELSE
            !    MAL(IP+K1,J)=U**2/G*MAL1(IP+K1,J)+4.*PI*SE(IP+K1)
            !END IF

        END DO
    END IF
    END DO
!$omp end parallel

    !MAL1=MAL
    !CALL BRINV(MAL1,NL)
    !QG=0
    !DO I=1,NL
    !    DO J=1,NL 
    !        QG(I)=QG(I)+MAL1(I,J)*VERR(J)
    !    END DO
    !END DO
    CALL MATRIX_PREGMRES

    ALLOCATE(SP(NL+500),SPN(NL+500),SPX(NL+500),SPY(NL+500),SPZ(NL+500),SPXX(NL+500),SPXY(NL+500),SPXZ(NL+500),SPYY(NL+500),SPYZ(NL+500))
    !ALLOCATE(MGX_X(NL),MGX_Y(NL),MGX_Z(NL),MGY_Y(NL),MGY_Z(NL),MGX_XX(NL),MGY_XY(NL),MGZ_XZ(NL),MGX_XY(NL),MGY_YY(NL),MGZ_YZ(NL))
    !ALLOCATE(MGX_XZ(NL),MGY_YZ(NL))
    ALLOCATE(SIGM0(NL+500),SIGM(NL+500))

    SIGM0=0.
    SIGM=QG
    VERR0=VERR
    !SIGM=0.
    !SIGM0=0.
    
    !GOTO 909

    !DO ITER=1,10
    
    SPN=0.;SP=0.
    DO IP=1,M
        DO J=1,NL
            SP(IP)=SP(IP)+SA(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            IF(IP.NE.J) THEN
                SPN(IP)=SPN(IP)+SH(IP,J)*SIGM(J)    !求解I点速度势沿法向导数（物面）
            ELSE
                SPN(IP)=SPN(IP)+4.*PI*SE(IP)*SIGM(J)    !求解I点速度势沿法向导数（物面）
            END IF
        END DO
    END DO

    OPEN(2,FILE='PHIS.DAT')
    DO I=1,M
        WRITE(2,*)I,SP(I),SIGM(I)
    END DO
    CLOSE(2)
    
    WRITE(*,*)"VERR..........OK"

    !GOTO 919

    MAL=0.
    DO IP=1,M
        DO J=1,NL
            IF(IP.EQ.J) THEN
            !IF(IP+K1.EQ.J.AND.IDS(IP+K1).NE.1) THEN
		        MAL(IP,J)=-4.*PI*SE(IP) !+SUMHIJ(IP+K1)
                !MAL(IP+K1,J)=2.*PI
            ELSE
                MAL(IP,J)=-SH(IP,J)
            END IF
        END DO  
    END DO

    DO IP=M+1,NL
        DO J=1,NL
            MAL(IP,J)=1./(2.*G)*2.*U**2*SXX(IP,J)+SZ(IP,J)
        END DO
    END DO

    !GOTO 909

    SPX=0.;SPY=0.;SPZ=0.;SPXX=0.;SPXY=0.;SPXZ=0.;SPYY=0.;SPYZ=0.
    DO IP=M+1,NL
        DO J=1,NL
            SP(IP)=SP(IP)+SA(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPX(IP)=SPX(IP)+SX(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPY(IP)=SPY(IP)+SY(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPZ(IP)=SPZ(IP)+SZ(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPXX(IP)=SPXX(IP)+SXX(IP,J)*SIGM(J)
            SPXY(IP)=SPXY(IP)+SXY(IP,J)*SIGM(J)
            SPXZ(IP)=SPXZ(IP)+SXZ(IP,J)*SIGM(J)
            SPYY(IP)=SPYY(IP)+SYY(IP,J)*SIGM(J)
            SPYZ(IP)=SPYZ(IP)+SYZ(IP,J)*SIGM(J)
        END DO
    END DO

    DO IP=1,M
        VERR(IP)=(U*VECN(IP,1)-SPN(IP))
        !VERR(IP)=U*VECN(IP,1)
    END DO
    
    DO IP=M+1,NL
        VERR(IP)=(1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U**2*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPY(IP)*SPXY(IP)-2.*SPX(IP)*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U**2*SPY(IP)*SPXY(IP)-2.*SPY(IP)*SPX(IP)*SPXY(IP)-2.*SPY(IP)*SPY(IP)*SPYY(IP)-2.*SPY(IP)*SPZ(IP)*SPYZ(IP))&
        +SPZ(IP))

        !VERR(IP)=(1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        !-1./2./G*(2.*U**2*SPX(IP)*SPXX(IP))&
        !-1./2./G*(2.*U**2*SPY(IP)*SPXY(IP))&
        !+SPZ(IP))        
        !VERR(IP)=0.
    END DO

    DO IP=M+1,NL
		IF((((IP)/(NL/4))*(NL/4)).EQ.IP)THEN
			WRITE(*,200)IP,NL 
    200	    FORMAT('MATRIX NON',3X,'6HINODE=',I8,5X,I8)
		END IF
        DO J=1,NL
            GOTO 902
            DO K=1,NL
                !GOTO 902
                MGX_X=0.;MGX_Y=0.;MGX_Z=0.;MGY_Y=0.;MGY_Z=0.;MGX_XX=0.;MGY_XY=0.;MGZ_XZ=0.;MGX_XY=0.;MGY_YY=0.;MGZ_YZ=0.;;MGX_XZ=0.;;MGY_YZ=0.
                GOTO 901
                DO K1=1,NL
                  MGX_X=MGX_X+&
                      (SX(IP,K)*SX(IP,K1)*SIGM(K1)+SX(IP,K)*SX(IP,K1)*SIGM(K1))    !PHIX PHIX
                  MGX_Y=MGX_Y+&
                      (SX(IP,K)*SY(IP,K1)*SIGM(K1)+SY(IP,K)*SX(IP,K1)*SIGM(K1))    !PHIX PHIY
                  MGX_Z=MGX_Z+&
                      (SX(IP,K)*SZ(IP,K1)*SIGM(K1)+SZ(IP,K)*SX(IP,K1)*SIGM(K1))    !PHIX PHIX
                  MGY_Y=MGY_Y+&
                      (SY(IP,K)*SY(IP,K1)*SIGM(K1)+SY(IP,K)*SY(IP,K1)*SIGM(K1))    !PHIY PHIY
                  MGY_Z=MGY_Z+&
                      (SY(IP,K)*SZ(IP,K1)*SIGM(K1)+SZ(IP,K)*SY(IP,K1)*SIGM(K1))    !PHIY PHIZ
                  MGX_XX=MGX_XX+&
                       SX(IP,K)*SXX(IP,K1)*SIGM(K1)+SXX(IP,K)*SX(IP,K1)*SIGM(K1)    !PHIX PHIXX
                  MGY_XY=MGY_XY+&
                       SY(IP,K)*SXY(IP,K1)*SIGM(K1)+SXY(IP,K)*SY(IP,K1)*SIGM(K1)    !PHIY PHIXY
                  MGZ_XZ=MGZ_XZ+&
                       SZ(IP,K)*SXZ(IP,K1)*SIGM(K1)+SXZ(IP,K)*SZ(IP,K1)*SIGM(K1)    !PHIZ PHIXZ
                  MGX_XY=MGX_XY+&
                       SX(IP,K)*SXY(IP,K1)*SIGM(K1)+SXY(IP,K)*SX(IP,K1)*SIGM(K1)    !PHIX PHIXY
                  MGY_YY=MGY_YY+&
                       SY(IP,K)*SYY(IP,K1)*SIGM(K1)+SYY(IP,K)*SY(IP,K1)*SIGM(K1)    !PHIY PHIYY
                  MGZ_YZ=MGZ_YZ+&
                       SZ(IP,K)*SYZ(IP,K1)*SIGM(K1)+SYZ(IP,K)*SZ(IP,K1)*SIGM(K1)    !PHIZ PHIYZ
                END DO
                901 CONTINUE
                MAL(IP,J)=MAL(IP,J)+&
               (-U/G*SX(IP,J)*SXX(IP,K)-U/G*SXX(IP,J)*SX(IP,K)&    !-2U/2G PHIX PHIXX
                -U/G*SY(IP,J)*SXY(IP,K)-U/G*SXY(IP,J)*SY(IP,K)&    !-2U/2G PHIY PHIXY
                -U/G*SZ(IP,J)*SXZ(IP,K)-U/G*SXZ(IP,J)*SZ(IP,K)&    !-2U/2G PHIZ PHIXZ
                -U/G*SX(IP,J)*SXX(IP,K)-U/G*SXX(IP,J)*SX(IP,K)&    !-2U/2G PHIX PHIXX
                -U/G*SY(IP,J)*SXY(IP,K)-U/G*SXY(IP,J)*SY(IP,K))*SIGM(K) !&    !-2U/2G PHIY PHIXY
                
                !+1./G*SX(IP,J)*MGX_XX*SIGM(K)+1./G*SX(IP,J)*MGX_XX*SIGM(K)&
                !+1./G*SXX(IP,J)*MGX_X*SIGM(K)&    !2/2G  PHIX PHIX PHIXX

                !+2./G*SX(IP,J)*MGY_XY*SIGM(K)+2./G*SY(IP,J)*MGX_XY*SIGM(K)&
                !+2./G*SXY(IP,J)*MGX_Y*SIGM(K)&    !4/2G  PHIX PHIY PHIXY

                !+1./G*SX(IP,J)*MGZ_XZ*SIGM(K)+1./G*SZ(IP,J)*MGX_XZ*SIGM(K)&
                !+1./G*SXZ(IP,J)*MGX_Z*SIGM(K)&    !2/2G  PHIX PHIZ PHIXZ

                !+1./G*SY(IP,J)*MGY_YY*SIGM(K)+1./G*SY(IP,J)*MGY_YY*SIGM(K)&
                !+1./G*SYY(IP,J)*MGY_Y*SIGM(K)&    !2/2G  PHIY PHIY PHIYY

                !+1./G*SY(IP,J)*MGZ_YZ*SIGM(K)+1./G*SZ(IP,J)*MGY_YZ*SIGM(K)&
                !+1./G*SYZ(IP,J)*MGY_Z*SIGM(K))    !2/2G  PHIY PHIZ PHIYZ
                
                902 CONTINUE
            END DO
        END DO
    END DO

    WRITE(*,*)"MAL..........OK"


    OPEN(601,FILE='DIAG_MATRIX.DAT')
    DO IP=M+1,NL 
        WRITE(601,*)IP,MAL(IP,IP)
    END DO
        
    !909 CONTINUE

    !CALL MATRIX_PREGMRES
    !CALL BLLU
    !MAL1=MAL
    !CALL BRINV(MAL1,NL)
    !QG=0
    !DO I=1,NL
    !    DO J=1,NL 
    !        QG(I)=QG(I)+MAL1(I,J)*VERR(J)
    !    END DO
    !END DO
    CALL MATRIX_PREGMRES

    SIGM0=SIGM
    VERR0=VERR

	OPEN(1,FILE='PHI.DAT')
	DO I=1,NL !+NFPOINT
        !READ(1,*)QG(I)    
    END DO
	CLOSE(1)
    SIGM=SIGM0-QG

    NORM1=0.
    NORM2=0
    DO IP=1,NL
        NORM1=NORM1+(SIGM0(IP)-SIGM(IP))**2
        !NORM2=NORM2+SIGM(IP)**2
        NORM2=NORM2+VERR(IP)**2
    END DO
    NORM1=SQRT(NORM1)/SQRT(NORM2)
    NORM1=SQRT(NORM2)/NL
    NORM1=MAXVAL(VERR(:))
    WRITE(*,*)"NEWTON RESIDUAL",ITER,NORM1
    !IF(NORM1.LE.1.E-5) EXIT

    !END DO

    !STOP
    909 CONTINUE
    !GOTO 919
!****************************************************************************************
    DO ITER=1,40

    SPN=0.;SP=0.
    DO IP=1,M
        DO J=1,NL
            SP(IP)=SP(IP)+SA(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            IF(IP.NE.J) THEN
                SPN(IP)=SPN(IP)+SH(IP,J)*SIGM(J)    !求解I点速度势沿法向导数（物面）
            ELSE
                SPN(IP)=SPN(IP)+4.*PI*SE(IP)*SIGM(J)    !求解I点速度势沿法向导数（物面）
            END IF
        END DO
    END DO

    SPX=0.;SPY=0.;SPZ=0.;SPXX=0.;SPXY=0.;SPXZ=0.;SPYY=0.;SPYZ=0.
    DO IP=M+1,NL
        DO J=1,NL
            SP(IP)=SP(IP)+SA(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPX(IP)=SPX(IP)+SX(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPY(IP)=SPY(IP)+SY(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPZ(IP)=SPZ(IP)+SZ(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            SPXX(IP)=SPXX(IP)+SXX(IP,J)*SIGM(J)
            SPXY(IP)=SPXY(IP)+SXY(IP,J)*SIGM(J)
            SPXZ(IP)=SPXZ(IP)+SXZ(IP,J)*SIGM(J)
            SPYY(IP)=SPYY(IP)+SYY(IP,J)*SIGM(J)
            SPYZ(IP)=SPYZ(IP)+SYZ(IP,J)*SIGM(J)
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
    WRITE(*,*)"VERR..........OK"

    OPEN(601,FILE='DIAG_MATRIX.DAT')
    DO IP=M+1,NL 
        WRITE(601,*)IP,SX(IP,IP),SPX(IP),SXY(IP,IP)
    END DO
        
        NYI(:,1)=VERR(:)-VERR0(:)
        NRT(1,:)=SIGM(:)-SIGM0(:)
        NRI(:,1)=SIGM(:)-SIGM0(:)
        RTR=0.
        DO I=1,NL
            RTR=RTR+(SIGM(I)-SIGM0(I))**2
        END DO
        MAL=MAL+MATMUL(NYI-MATMUL(MAL,NRI),NRT)/RTR

    !909 CONTINUE

    !CALL MATRIX_PREGMRES
    !CALL BLLU
    !MAL1=MAL
    !CALL BRINV(MAL1,NL)
    !QG=0
    !DO I=1,NL
    !    DO J=1,NL 
    !        QG(I)=QG(I)+MAL1(I,J)*VERR(J)
    !    END DO
    !END DO
    CALL MATRIX_PREGMRES

    SIGM0=SIGM
    VERR0=VERR
	!OPEN(1,FILE='PHI.DAT')
	!DO I=1,NL !+NFPOINT
    !    READ(1,*)QG(I)    
    !END DO
	!CLOSE(1)
    SIGM=SIGM0-QG

    NORM1=0.
    NORM2=0.
    DO IP=1,NL
        NORM1=NORM1+(SIGM0(IP)-SIGM(IP))**2
        !NORM2=NORM2+SIGM(IP)**2
        NORM2=NORM2+VERR(IP)**2
    END DO
    NORM1=SQRT(NORM1)/SQRT(NORM2)
    NORM1=SQRT(NORM2)/NL
    NORM1=MAXVAL(VERR(:))
    WRITE(*,*)"NEWTON RESIDUAL",ITER,NORM1
    IF(NORM1.LE.1.E-6) EXIT

    END DO

    !STOP
    
    919 CONTINUE

30  DEALLOCATE(MAR,VER,MAL1,IDS)											
	CLOSE(16)
!释放内存

	!DEALLOCATE(SH,SA,SA1,SX,SZ,SXX)
	WRITE(*,*)'		END MATRIX'

	RETURN
END SUBROUTINE
