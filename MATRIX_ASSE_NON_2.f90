SUBROUTINE MATRIX_NON_2     !计算线性问题
!**************************************************************
!CGENERATE THE REAL MATRIX OF THE INTEGRATION EQUATION
!**************************************************************
	USE GREENMOD
	USE CUMOD

    IMPLICIT NONE

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2
    REAL,ALLOCATABLE:: DL(:),NYI(:,:),NRI(:,:),NRT(:,:),verr0(:)
    REAL:: MGX_X,MGX_Y,MGX_Z,MGY_Y,MGY_Z,MGX_XX,MGY_XY,MGZ_XZ,MGX_XY,MGY_YY,MGZ_YZ,MGX_XZ,MGY_YZ
    REAL:: NORM1,NORM2,RTR
    INTEGER:: NAMRK
    INTEGER:: I,J,K,K1,IP,M,L,NL,CP,K2,NMARK,ITER
    INTEGER:: NPREC,NITER
    REAL:: SORW,ALPHA
    CHARACTER*7 TITLE,TITLE2

    ALLOCATE(DL(NBPOINT))

    M=NBPOINT
    L=NFPOINT
    NL=M+L

    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(99,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\'//TRIM(ADJUSTL(TITLE))//'\WAVE_PRO_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
    IF(ITTE.EQ.1) OPEN(117,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\NEWTON_RESIDU.DAT')


    IF(MDIFF.EQ.2) THEN
!   用差分求解
    
    WRITE(*,*)"用差分求解"
    ALLOCATE(DVALU(M+L,M+L))
    DVALU(1:M+L,1:M+L)=SX(1:M+L,1:M+L)
    CALL DIFF_COEMATRIX_DECOMP
    DO I=M+1,NL
        !WRITE(*,*)SXX(I,I),DVMX(I,I)
	    DO J=1,NL
            SXX(I,J)=DVMX(I,J)
            SXY(I,J)=DVMY(I,J)
        END DO
    END DO
    DEALLOCATE(DVMX,DVMY)

    DVALU(1:M+L,1:M+L)=SY(1:M+L,1:M+L)
    CALL DIFF_COEMATRIX_DECOMP
    DO I=M+1,NL
	    DO J=1,NL
            SYY(I,J)=DVMY(I,J)
        END DO
    END DO
    DEALLOCATE(DVMX,DVMY)
    
    DVALU(1:M+L,1:M+L)=SZ(1:M+L,1:M+L)
    CALL DIFF_COEMATRIX_DECOMP
    DO I=M+1,NL
	    DO J=1,NL
            SXZ(I,J)=DVMX(I,J)
            SYZ(I,J)=DVMY(I,J)
        END DO
    END DO
    DEALLOCATE(DVMX,DVMY)
    DEALLOCATE(DVALU)
    END IF
!   END

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



	20 ALLOCATE(MAR(NL,NL),VER(NL),MAL1(NL,NL))
    IF(NIT.EQ.1.AND.ITTE.EQ.1) THEN
        ALLOCATE(MAL(NL+500,NL+500),VERR(NL+500))	
    END IF
    ALLOCATE(NYI(NL,1),NRI(NL,1),NRT(1,NL),VERR0(NL))
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
    
    SORW=0.6
    SORW=0.8
    !SORW=1.6

    !MAL1=MAL
    !CALL BRINV(MAL1,NL)
    !QG=0
    !DO I=1,NL
    !    DO J=1,NL 
    !        QG(I)=QG(I)+MAL1(I,J)*VERR(J)
    !    END DO
    !END DO
    !CALL MATRIX_PREGMRES

    
    IF(ITTE.EQ.1) THEN
    NMARK=1
    CALL MATRIX_PRE0GMRES(NMARK)
    ELSE
            ALPHA=1.
            !IF(NTPN.LE.4000) THEN
            IF(FR.LE.0.35) THEN
            IF(ITTE.LE.6) THEN
                !ALPHA=1./6.*(ITTE)
            END IF
            ELSE
            IF(ITTE.LE.9) THEN
                !ALPHA=1./9.*(ITTE)
            END IF
            END IF

            SP0(NBPOINT+1:NBPOINT+NFPOINT)=SP(NBPOINT1+1:NBPOINT1+NFPOINT)*ALPHA
            SPX0(NBPOINT+1:NBPOINT+NFPOINT)=SPX(NBPOINT1+1:NBPOINT1+NFPOINT)*ALPHA
            SPY0(NBPOINT+1:NBPOINT+NFPOINT)=SPY(NBPOINT1+1:NBPOINT1+NFPOINT)*ALPHA
            SPZ0(NBPOINT+1:NBPOINT+NFPOINT)=SPZ(NBPOINT1+1:NBPOINT1+NFPOINT)*ALPHA
            SPXX0(NBPOINT+1:NBPOINT+NFPOINT)=SPXX(NBPOINT1+1:NBPOINT1+NFPOINT)*ALPHA
            SPXY0(NBPOINT+1:NBPOINT+NFPOINT)=SPXY(NBPOINT1+1:NBPOINT1+NFPOINT)*ALPHA
            SPXZ0(NBPOINT+1:NBPOINT+NFPOINT)=SPXZ(NBPOINT1+1:NBPOINT1+NFPOINT)*ALPHA
            SPYY0(NBPOINT+1:NBPOINT+NFPOINT)=SPYY(NBPOINT1+1:NBPOINT1+NFPOINT)*ALPHA
            SPYZ0(NBPOINT+1:NBPOINT+NFPOINT)=SPYZ(NBPOINT1+1:NBPOINT1+NFPOINT)*ALPHA  
            !END IF
    END IF


    IF(NIT.EQ.1.AND.ITTE.EQ.1) THEN
        ALLOCATE(SP(NL+500),SPN(NL+500),SPX(NL+500),SPY(NL+500),SPZ(NL+500),SPXX(NL+500),SPXY(NL+500),SPXZ(NL+500),SPYY(NL+500),SPYZ(NL+500))
        ALLOCATE(SP0(NL+500),SPX0(NL+500),SPY0(NL+500),SPZ0(NL+500),SPXX0(NL+500),SPXY0(NL+500),SPXZ0(NL+500),SPYY0(NL+500),SPYZ0(NL+500))
        ALLOCATE(SIGM0(NL+500),SIGM(NL+500))
    END IF
    !ALLOCATE(MGX_X(NL),MGX_Y(NL),MGX_Z(NL),MGY_Y(NL),MGY_Z(NL),MGX_XX(NL),MGY_XY(NL),MGZ_XZ(NL),MGX_XY(NL),MGY_YY(NL),MGZ_YZ(NL))
    !ALLOCATE(MGX_XZ(NL),MGY_YZ(NL))

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

    IF(ITTE.EQ.1) GOTO 919
    !GOTO 919

    !OPEN(2,FILE='PHIS.DAT')
    !DO I=1,M
        !WRITE(2,*)I,SP(I),SIGM(I)
    !END DO
    !CLOSE(2)
    
    WRITE(*,*)"VERR..........OK"

    MAL=0.
    DO IP=1,M
        DO J=1,NL
            IF(IP.EQ.J) THEN
            !IF(IP+K1.EQ.J.AND.IDS(IP+K1).NE.1) THEN
		        MAL(IP,J)=4.*PI*SE(IP) !+SUMHIJ(IP+K1)
                !MAL(IP+K1,J)=2.*PI
            ELSE
                MAL(IP,J)=SH(IP,J)
            END IF
        END DO  
    END DO

    DO IP=M+1,NL
        DO J=1,NL
            MAL(IP,J)=1./(2.*G)*2.*U**2*SXX(IP,J)+SZ(IP,J)
        END DO
    END DO

    DO IP=1,M
        !VERR(IP)=(U*VECN(IP,1)-SPN(IP))
        VERR(IP)=U*VECN(IP,1)
    END DO
    

    OPEN(981,FILE='TEST_VERR.DAT')
    DO IP=M+1,NL
        !VERR(IP)=(1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        !-1./2./G*(2.*U**2*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPY(IP)*SPXY(IP)-2.*SPX(IP)*SPZ(IP)*SPXZ(IP))&
        !-1./2./G*(2.*U**2*SPY(IP)*SPXY(IP)-2.*SPY(IP)*SPX(IP)*SPXY(IP)-2.*SPY(IP)*SPY(IP)*SPYY(IP)-2.*SPY(IP)*SPZ(IP)*SPYZ(IP))&
        !+SPZ(IP))

        VERR(IP)=-1./2./G*(-2.*U*SPX0(IP)*SPXX0(IP)-2.*U*SPY0(IP)*SPXY0(IP)-2.*U*SPZ0(IP)*SPXZ0(IP))&
        +1./2./G*(2.*U*SPX0(IP)*SPXX0(IP)-2.*SPX0(IP)*SPX0(IP)*SPXX0(IP)-2.*SPX0(IP)*SPY0(IP)*SPXY0(IP)-2.*SPX0(IP)*SPZ0(IP)*SPXZ0(IP))&
        +1./2./G*(2.*U*SPY0(IP)*SPXY0(IP)-2.*SPY0(IP)*SPX0(IP)*SPXY0(IP)-2.*SPY0(IP)*SPY0(IP)*SPYY0(IP)-2.*SPY0(IP)*SPZ0(IP)*SPYZ0(IP))
        

        !VERR(IP)=(1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        !-1./2./G*(2.*U**2*SPX(IP)*SPXX(IP))&
        !-1./2./G*(2.*U**2*SPY(IP)*SPXY(IP))&
        !+SPZ(IP))        
        !VERR(IP)=0.
        !WRITE(*,*)VERR(IP)
        !WRITE(981,*)IP,VERR(IP)
    END DO
    !VERR0=VERR

    !GOTO 909

    WRITE(*,*)"MAL..........OK"

        
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
    !CALL MATRIX_PREGMRES
    !NMARK=3
    NMARK=1
    CALL MATRIX_PRE0GMRES(NMARK)
    SIGM0=SIGM
    VERR0=VERR
    DO I=1,NTPN
        !SIGM(I)=SIGM0(I)-SORW*QG(I)
        SIGM(I)=QG(I)
    END DO

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
        !VERR(IP)=(U*VECN(IP,1)-SPN(IP))
        VERR(IP)=(-U*VECN(IP,1)+SPN(IP))
        !VERR(IP)=-U*VECN(IP,1)
        !WRITE(981,*)IP,VERR(IP)
    END DO

    DO IP=M+1,NL
        VERR(IP)=(1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U**2*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPY(IP)*SPXY(IP)-2.*SPX(IP)*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U**2*SPY(IP)*SPXY(IP)-2.*SPY(IP)*SPX(IP)*SPXY(IP)-2.*SPY(IP)*SPY(IP)*SPYY(IP)-2.*SPY(IP)*SPZ(IP)*SPYZ(IP))&
        +SPZ(IP))
        WRITE(981,*)IP,VERR(IP)
    END DO

    NORM1=0.
    NORM2=0
    DO IP=1,NL
        !NORM1=NORM1+(SIGM0(IP)-SIGM(IP))**2
        !NORM2=NORM2+SIGM(IP)**2
        !NORM2=NORM2+VERR(IP)**2
    END DO
    !NORM1=SQRT(NORM1)/SQRT(NORM2)
    !NORM1=SQRT(NORM2)/NL
    NORM1=MAXVAL(ABS(VERR(1:NL)))
    WRITE(*,*)"NEWTON RESIDUAL",ITER,NORM1
    WRITE(117,*)"NEWTON RESIDUAL",ITTE,NORM1
    !IF(NORM1.LE.1.E-5) EXIT

    !END DO


    !STOP
    909 CONTINUE
    !GOTO 919
!****************************************************************************************
    NITER=ITTE*4 
    NITER=60
    DO ITER=1,NITER

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

    IF(MOD(ITTE,50).NE.0) GOTO 919

    DO IP=1,M
        !VERR(IP)=(U*VECN(IP,1)-SPN(IP))
        VERR(IP)=(-U*VECN(IP,1)+SPN(IP))
        !VERR(IP)=-U*VECN(IP,1)
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
    !WRITE(*,*)"VERR..........OK"

        
        !NYI(:,1)=VERR(:)-VERR0(:)
        !NRT(1,:)=SIGM(:)-SIGM0(:)
        !NRI(:,1)=SIGM(:)-SIGM0(:)
        !RTR=0.
        DO I=1,NL
        !    RTR=RTR+(SIGM(I)-SIGM0(I))**2
        END DO
        !MAL=MAL+MATMUL(NYI-MATMUL(MAL,NRI),NRT)/RTR

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
    !CALL MATRIX_PREGMRES
    NMARK=3
    NPREC=80
    !IF(ITER.GE.NPREC.AND.MOD(ITER,NPREC).EQ.0) NMARK=2;WRITE(*,*)ITER
    !IF(ITER.GE.NPREC.AND.MOD(ITER-1,NPREC).EQ.0) NMARK=1;WRITE(*,*)ITER
    CALL MATRIX_PRE0GMRES(NMARK)

    SIGM0=SIGM
    VERR0=VERR
	!OPEN(1,FILE='PHI.DAT')
	!DO I=1,NL !+NFPOINT
    !    READ(1,*)QG(I)    
    !END DO
	!CLOSE(1)
    DO I=1,NTPN
    SIGM(I)=SIGM0(I)-SORW*QG(I)
    END DO

    NORM1=0.
    NORM2=0.
    DO IP=1,NL
        NORM1=NORM1+(SIGM0(IP)-SIGM(IP))**2
        NORM2=NORM2+SIGM(IP)**2
        !NORM2=NORM2+VERR(IP)**2
    END DO
    NORM1=SQRT(NORM1)/SQRT(NORM2)
    !NORM1=SQRT(NORM2)/NL
    NORM1=MAXVAL(ABS(VERR(1:NL)))
    WRITE(*,*)"NEWTON RESIDUAL",ITER,NORM1
    IF(NORM1.LE.5.E-3) THEN
        WRITE(117,*)"NEWTON RESIDUAL",ITTE,NORM1
        EXIT
    END IF

    !IF(ITER.GE.NITER) THEN
    !    WRITE(*,*)"NO CONVERGENCE"; 
        !PAUSE
    !END IF

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

SUBROUTINE MATRIX_PRE0GMRES(NMARK)
!**************************************************************
!	GENERATE AND SOLVE THE INTEGRATION EQUATION
!**************************************************************


  USE GREENMOD
  USE CUMOD

  implicit none 

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF


  REAL TPR(30),HGIJ,HGIJ1,TOL1,RNORM1,BNORM,RNORM,TPVE
  INTEGER K,THR1,THR2,CP,ITR,ITER,BM,METH,i,j,M,L,KSTART,ITEROUT,N
  REAL,PARAMETER ::ITERMAX=1000
  REAL:: TOL  !TOL=1.0E-6
  REAL,ALLOCATABLE:: TPM(:)
  CHARACTER*8 CHAR_TIME1,CHAR_TIME2
  CHARACTER*8 CHAR_DATE1,CHAR_DATE2
  REAL:: TIME1,TIME2
  INTEGER:: HOUR1,MINI1,SECO1,HOUR2,MINI2,SECO2
  INTEGER:: NMARK,NTHREAD1

  NTHREAD1=NTHREAD
  NTHREAD=1
  
  IF(ITTE.EQ.1) THEN
    TOL=1.E-7
  ELSE
    TOL=1.E-6
    IF(NMARK.EQ.1) THEN    
        !TOL=1.E-7
    END IF
  END IF

  OPEN(66,FILE="TEMP.DAT")

  DO I=1,NTPN
    DO J=1,NTPN
        !WRITE(66,*)I,J,MAL(I,J)
    END DO
  END DO

  BM=NTPN/2. !200
  BM=2000
  BM=50
  BM=35
  write(*,*)"bm",bm
  M=NBPOINT
  L=NFPOINT
  N=M+L !+L

    OPEN(1,FILE='PHI.DAT')


  !WRITE(*,*)"M+L",M+L
  ALLOCATE(CI(M+L,M+L),APG(M+L,M+L),WI(M+L,3),BV(M+L,BM+50),BH(BM+1,BM+50))
  ALLOCATE(BB(M+L),RG(M+L),PG(M+L),BNORM1(NTHREAD),GG(BM+1+50),YG(BM+1+50),CG(BM+50),SG(BM+50))
  ALLOCATE(TPM(N))



  !VERR=0
  TOL1=1
  ITEROUT=0
  KSTART=0
  CI=0
  CG=0
  SG=0

  METH=0         
  IF(NMARK.EQ.1) THEN      
  CALL PRE_CONDITION(METH)
  END IF

  DO I=1,N
    !WRITE(66,*)MAL(I,I),DAG(I,I),DA_1(I,I)
  END DO
  WRITE(*,*)"OK"

! 预条件结束
    CALL TIME(CHAR_TIME1)
	CALL DATE(CHAR_DATE1)


!   加预设条件   
    
    DO I=1,NPS(CP)
      DO J=1,N
        !APG(THR1+I,J)=0.
        DO K=1,N
            !APG(THR1+I,J)=APG(THR1+I,J)+CI(THR1+I,K)*MAL(K,J)
        END DO
      END DO
    END DO
    
    
    !AG(THR1+1:THR2,:)=APG(THR1+1:THR2,:)
    
    !AG(1:NTPN,:)=MAL(1:NTPN,:)
    !APG(1:NTPN,:)=MAL(1:NTPN,:)


    !BB(THR1+1:THR2)=APG(1,THR1+1:THR2)
    BB(1:N)=VERR(1:N)

    BNORM1(1)=0
    DO I=1,N !NPS(CP)
      BNORM1(1)=BNORM1(1)+BB(I)*BB(I)
    END DO
!   向量的模 
    !CALL FNORM(BB(K1:K2),BB(K1:K2),N,1,BNORM1(CP))


	!BNORM=SQRT(SUM(BNORM1(1:NTHREAD)))
    BNORM=BNORM1(1)

    IF(NIT.EQ.1) THEN
      !QG=0.
    END IF
    QG=0
    !RG=BB

    RG=0.
    DO I=1,N
      DO J=1,N
        !RG(I)=RG(I)+BB(I)-AG(I,J)*QG(J)
        RG(I)=RG(I)-MAL(I,J)*QG(J)
      END DO
      RG(I)=RG(I)+BB(I)
    END DO


    RNORM1=0.
    DO I=1,N
      RNORM1=RNORM1+RG(I)*RG(I)    !改动前
      !RNORM1=RNORM1+ABS(RG(I))
    END DO

!   预条件
    METH=1               
    BAX=RG
    CALL PRE_CONDITION(METH)

    DO I=1,N
        !WRITE(66,*)X_BAR(I),BAX(I),UA(I,I)
    END DO
    RG=RG0

    !STOP

    IF(SQRT(RNORM1).LE.TOL) GOTO 109    !改动前
    !IF(RNORM1/N.LE.TOL) GOTO 109    

 88   CONTINUE
 
!   从新设置重启动数
    IF(ITEROUT*BM+ITER.GE.5*BM) THEN 
        TOL=TOL+2.E-6
    END IF
!**********************************************

    RNORM1=0
    DO I=1,N
      RNORM1=RNORM1+RG(I)*RG(I)
    END DO
      
!   SET UP V^1 AND G^0


    !WRITE(*,*)'K1',K1,NPS(CP) 

	
    !CALL FNORM(R,R,N,NS,RNORM1)
!C	WRITE(6,*) RNORM1
	RNORM=SQRT(RNORM1)
	
    DO I=1,N
      PG(I)=RG(I)/RNORM
	END DO

    GG=0

    GG(1)=RNORM
!   BEGIN ITERATION
    ITER=1

    !CG=0.;SG=0.

!	BWRITE(66,*) 'BETA',G(1)
    DO WHILE((ITER.LE.BM).AND.(TOL1.GT.TOL).AND.(ITEROUT*BM+ITER.LE.ITERMAX))
    !DO WHILE((ITER.LE.BM).AND.(TOL1.GT.TOL))


        BV(1:NTPN,ITER)=PG(1:NTPN)
        !CALL EQUAL(BV(:,1,ITER),P,N,NS)
 !  WRITE(66,*)'ITER',ITER
!	WRITE(66,*) ((BV(I1,J1,ITER),I1=1,N),J1=1,NS)

!   FORM  AV(ITER)	 
        APG(1,1:NTPN)=0
        DO I=1,NTPN
          DO J=1,N
            APG(1,I)=APG(1,I)+MAL(I,J)*PG(J)
          END DO
        END DO


	  !CALL MATAV(A,P,AP,N,NS)  !A*P
!	IF(ITER.EQ.1) WRITE(*,*) 'AP', AP(3,1)
!   INITIALIZE V^(ITER+1) WITH AV(ITER)
        PG(1:NTPN)=APG(1,1:NTPN)
        !WI(K1+1:K2,1)=APG(1,K1+1:K2)
     !GOTO 11
! MAKE V^(ITER+1) ORTHOGONAL TO V^(I),I<=ITER

!   预条件
    METH=1               
    BAX=PG
    CALL PRE_CONDITION(METH)
    PG=RG0
        DO  J=1,ITER
! 1)正交化
	   !CALL FNORM(BV(1,1,J),AP,N,NS,HIJ)
           
            HGIJ=0.    
            DO I=1,N
              !HGIJ=HGIJ+BV(I,J)*APG(1,I)
              HGIJ=HGIJ+BV(I,J)*PG(I)
              !HGIJ=HGIJ+BV(1,J)*APG(1,I)
              !WRITE(66,*)"ABS(GG(ITER+1))",bV(I,J),PG(I)
            END DO
            BH(J,ITER)=HGIJ

     !WRITE(*,*)"ABS(GG(ITER+1))",V(I,J),PG(I)
          !TPR(CP)=0
          !DO I=1,NPS(CP)
            !TPR(CP)=TPR(CP)+BV(K1+I,J)*PG(1,K1+I)
          !END DO      
         !DO I=1,N
           !WI(I,NS,J1+1)=WI(I,NS,J1)-BV(I,NS,J)*TPR
         !END DO
         !WRITE(*,*)SQRT(RNORM/RNORM1)

          DO  I=1,NTPN
		    !P(I,K)=P(I,K)-TPR*BV(I,K,J)
            PG(I)=PG(I)-HGIJ*BV(I,J)
	      ENDDO 
	      !BH(J,ITER)=TPR

          !BH(J,ITER)=HGIJ


          !END DO 
          !CALL FNORM(WI(:,:,J1),WI(:,:,J1),N,NS,RNORM)
          !CALL FNORM(WI(:,:,J1+1),WI(:,:,J1+1),N,NS,RNORM1)
          !IF(SQRT(RNORM/RNORM1).LE.TH) EXIT         
	    END DO 
        !WRITE(*,*)"BH(J,ITER)",BH(J,ITER)



       !WRITE(*,*)"PG",PG(1:N)

       RNORM1=0
       DO I=1,N
         RNORM1=RNORM1+PG(I)*PG(I)
       END DO     

       !PG(K1+1:K2)=APG(1,K1+1:K2)

!NORMALIZE V^(ITER+1)
! 2)标准化
	 !CALL FNORM(P,P,N,NS,RNORM1)
	 !RNORM =REAL(RNORM1)
     !RNORM=SQRT(RNORM1)

       DO  I=1,N
		  PG(I)=PG(I)/SQRT(RNORM1)
	   END DO

     BH(ITER+1,ITER)=SQRT(RNORM1) 
!	 WRITE(66,*)'H OLD BEFORE ', (BH(I,ITER),I=1,ITER+1)    

! APPLY ROTATIONS TO NEW H COLUMN

! 3)更新V和H
     DO I=1,ITER-1
	   HGIJ=BH(I,ITER)
	   HGIJ1=BH(I+1,ITER)
	   BH(I,ITER)=HGIJ*CG(I)+HGIJ1*SG(I)
	   BH(I+1,ITER)=-HGIJ*SG(I)+HGIJ1*CG(I)
	 ENDDO
!	 WRITE(66,*)'H OLD AFTER ', (BH(I,ITER),I=1,ITER)  
!    COMPUTE NEW ROTATIONS
     HGIJ=BH(ITER,ITER)
	 HGIJ1=BH(ITER+1,ITER) 	   
     
     CALL GIVENS(HGIJ,HGIJ1,CG(ITER),SG(ITER))
!    WRITE(66,*) 'GIVENS',C(ITER),S(ITER)
!    APPLY NEW ROTATIONS
     BH(ITER,ITER)= HGIJ*CG(ITER)+HGIJ1*SG(ITER)
     BH(ITER+1,ITER)=0.
!    WRITE(66,*)'H NEW AFTER ', (BH(I,ITER),I=1,ITER+1) 

     !WRITE(*,*)"ABS(GG(ITER+1))",HGIJ
     HGIJ=GG(ITER)
     !WRITE(*,*)"ABS(GG(ITER+1))",HGIJ
     GG(ITER)=CG(ITER)*HGIJ
	 GG(ITER+1)=-SG(ITER)*HGIJ
     TOL1=ABS(GG(ITER+1))/BNORM
!    WRITE(66,*)'G  AFTER ', (G(I ),I=1,ITER+1) 
!C   TOL1=CABS(G(ITER+1))
!    WRITE(6,*) 'ITER',ITER
!	 WRITE(6,*) 'TOL1',TOL1
!	 WRITE(66,*) 'TOLERANCE',TOL1
!    WRITE(151,*) ITEROUT*M+ITER,LOG10(TOL1)
	 ITER=ITER+1
     
        !WRITE(66,*) 'TOL1',TOL1
        !WRITE(211,*)ITER,TOL1
        
        !WRITE(*,*)ITER,TOL1
        IF(ITER.EQ.BM) WRITE(*,*)ITER,TOL1

    END DO


    ITER=ITER-1
    !WRITE(66,*)
    !WRITE(66,*) 'TOL1',BM,TOL1
 
!   COMPUTE  LEAST SQUARE  SOLUTION
    DO I=1,ITER
	  YG(I)=GG(I)
	ENDDO
!	WRITE(66,*)'G'
!	WRITE(66,*) (G(I),I=1,ITER)
	DO I=ITER ,1,-1 
 	   YG(I)=YG(I)/BH(I,I)   
 	   DO J=I-1,1,-1
          YG(J)=YG(J)-BH(J,I)*YG(I)
       ENDDO
    ENDDO
!	WRITE(66,*) 'Y'
!	WRITE(66,*)(Y(J),J=1,ITER)
!	 WRITE(66,*)'BV'
	
    !DO  K=1,ITER
!	 WRITE(66,*)'ITER',K
	
!       WRITE(66,*)(BV(I,1,K),I=1,N)
	!ENDDO

!    ACQUIRE ITERATION SOLUTION
	 DO I=1,N
	   DO K=1,ITER
         QG(I)=QG(I)+YG(K)*BV(I,K)
	   ENDDO
	 ENDDO
!	 WRITE(66,*) 'Q'
!    WRITE(66,*) (Q(I,1),I=1,N)

     IF((TOL1.LT.TOL).AND.(ITEROUT*BM+ITER.LE.ITERMAX)) THEN
     !IF((TOL1.LT.TOL)) THEN
!      ITERATIONS END
	 GOTO 98
	 ENDIF
     
     KSTART=0
	 IF((TOL1.GT.TOL).AND.(ITEROUT*BM+ITER.LT.ITERMAX)) KSTART=1
     !IF((TOL1.GT.TOL)) KSTART=1

     IF( (ITEROUT*BM+ITER.GT.ITERMAX).AND.(TOL1.GT.TOL) ) THEN
     
     !IF( (TOL1.GT.TOL) ) THEN
!	   WRITE(6,*) 'WARNING NO CONVERGENCE  AFTER ' 
!	   WRITE(6,*) ITEROUT*M+ITER,'ITERATIONS'
!	   WRITE(6,*) 'TOLRENCE  IS',TOL1
!C	   WRITE(*,*) 'ITEROUT',ITEROUT
       WRITE(*,*) 'TOLRENCE  IS',TOL1

       GOTO  109
	 ENDIF

    !WRITE(*,*) 'TOLRENCE  IS',TOL1

	IF(KSTART.EQ.1) THEN

	ITEROUT=ITEROUT+1
!C  　    
!       WRITE(6,*) 'ITEROUT',ITEROUT     

    APG=0
    DO I=1,N
      DO J=1,N
        APG(1,I)=APG(1,I)+AG(I,J)*QG(J)
      END DO
    END DO
                                                                                                                                                                                                                                                          
    !CALL  MATAV(A,Q,AP,N,NS)
    DO  I=1,N
	  RG(I)=BB(I)-APG(1,I)
	END DO

!   预条件
    METH=1               
    BAX=RG
    CALL PRE_CONDITION(METH)
    RG=RG0

!   ITERATIONS BEGIN AGAIN
      
	GOTO 88
    ENDIF
  98  CONTINUE
!      WRITE(6,*) 'CONGRATULATIONS CONVERGENCE ACQUIRE ' 
!	   WRITE(6,*)    'AFTER', ITEROUT*M+ITER,' ITERATIONS'
	  !WRITE(*,*) 'TOLRENCE  IS',TOL1  

  99  CONTINUE 
    BB=QG

	ITR=ITEROUT*BM+ITER
	!CLOSE(66)
    WRITE(111,*) 'ITERATIONS     ',ITR
    WRITE(111,*) 

  109  CONTINUE
  
  WRITE(*,*)"ITR=   ",ITR,TOL1     
        
  DO I=1,M+L !+L
     IF(I.LE.NBPOINT) THEN
        !PHIS(I)=QG(I)
        !WRITE(*,*)X(I)
     ELSE
        !PHIS_N(I)=QG(I)
     END IF
     WRITE(1,*)QG(I)
     !WRITE(*,*)I,QG(I)
  END DO
  
  !STOP

  !WRITE(*,*)QG(30)

  !GMRES 结束
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
	WRITE(111,*)'BEGIN TIME        : ', CHAR_TIME1,'  ON  ',CHAR_DATE1
    WRITE(111,*)'END   TIME  SOLVER: ', CHAR_TIME2,'  ON  ',CHAR_DATE2
    READ(CHAR_TIME1(1:2),*)HOUR1
    READ(CHAR_TIME1(4:5),*)MINI1
    READ(CHAR_TIME1(7:8),*)SECO1
    READ(CHAR_TIME2(1:2),*)HOUR2
    READ(CHAR_TIME2(4:5),*)MINI2
    READ(CHAR_TIME2(7:8),*)SECO2
    TGMR=(HOUR2-HOUR1)*3600+(MINI2-MINI1)*60+SECO2-SECO1

    WRITE(311,*)FR,TILU,TGMR

  !预条件释放内存
    METH=2
    IF(NMARK.EQ.2.OR.ITTE.EQ.1) THEN      
        !CALL PRE_CONDITION(METH)
    END IF

!   END GMRES MUDULE

14	CLOSE(16)   
    CLOSE(1)
!C 释放内存
  DEALLOCATE(CI,APG,WI,BV,BH)
  DEALLOCATE(BB,RG,PG,BNORM1,GG,YG,CG,SG)
  DEALLOCATE(TPM)
  
  WRITE(*,*)"END GMRES"
    NTHREAD=NTHREAD1
  
END SUBROUTINE



