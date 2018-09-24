SUBROUTINE MATRIX_NON_4     !计算线性问题
!**************************************************************
!CGENERATE THE REAL MATRIX OF THE INTEGRATION EQUATION
!**************************************************************
	USE GREENMOD
	USE CUMOD
    USE INPUTDATA
    !IMPLICIT NONE

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2
    CHARACTER*1 creturn
    REAL,ALLOCATABLE:: DL(:),NYI(:,:),NRI(:,:),NRT(:,:),verr0(:),RTH(:)
    REAL:: MGX_X,MGX_Y,MGX_Z,MGY_Y,MGY_Z,MGX_XX,MGY_XY,MGZ_XZ,MGX_XY,MGY_YY,MGZ_YZ,MGX_XZ,MGY_YZ
    REAL:: NORM1,NORM2,RTR,NORMB(2)
    REAL:: SORW
    INTEGER:: NAMRK
    INTEGER:: I,J,K,K1,IP,M,L,NL,CP,K2,NMARK,ITER
    INTEGER:: NPREC,NITER
    REAL:: XTR1(100),ZTR1(100),TAOX1(100),ALPHA
    CHARACTER*7 TITLE,TITLE2
    ALLOCATE(DL(NBPOINT))

    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(99,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\'//TRIM(ADJUSTL(TITLE))//'\WAVE_PRO_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
    !IF(ITTE.EQ.1) OPEN(117,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\NEWTON_RESIDU.DAT')


    !CALL PADCOEF
    !CALL DIFFCORD_DECOMP   ! 计算坐标对计算平面导数
    !CALL DIFFCORD_MATRIX_DECOMP ! 计算差分算子


!采用bin文件格式可以减小文件大小，加快读写速度 
    open(16,file='coefuse.bin',form='binary',access='sequential',&
          status='old')
    rewind(16)
	!WRITE(*,*)
	!WRITE(*,*)'SUB MATRIX.........'
	!ALLOCATE(SH(NTPN,NTPN),SX(NTPN,NTPN),SZ(NTPN,NTPN),SA(NTPN,NTPN),SA1(NTPN,NTPN),SXX(NTPN,NTPN))
	
    M=NBPOINT
    L=NFPOINT
    NL=NTPN
    !NL=M+L+NBTPOINT

    
    IF(MDIFF.EQ.2) THEN
!   用差分求解
    
    !WRITE(*,*)"用差分求解"
    ALLOCATE(DVALU(M+L,NL))
    DVALU(1:M+L,1:NL)=SX(1:M+L,1:NL)
    CALL DIFF_COEMATRIX_DECOMP
    DO I=M+1,M+L
        !WRITE(*,*)SXX(I,I),DVMX(I,I)
	    DO J=1,NL
            SXX(I,J)=DVMX(I,J)
            SXY(I,J)=DVMY(I,J)
        END DO
    END DO
    DEALLOCATE(DVMX,DVMY)

    DVALU(1:M+L,1:NL)=SY(1:M+L,1:NL)
    CALL DIFF_COEMATRIX_DECOMP
    DO I=M+1,M+L
	    DO J=1,NL
            !SXY(I,J)=DVMX(I,J)
            SYY(I,J)=DVMY(I,J)
        END DO
    END DO
    DEALLOCATE(DVMX,DVMY)
    
    DVALU(1:M+L,1:NL)=SZ(1:M+L,1:NL)
    CALL DIFF_COEMATRIX_DECOMP
    DO I=M+1,M+L
	    DO J=1,NL
            SXZ(I,J)=DVMX(I,J)
            SYZ(I,J)=DVMY(I,J)
        END DO
    END DO
    DEALLOCATE(DVMX,DVMY)
    DEALLOCATE(DVALU)
    END IF
!   END
    creturn = achar(13)
    writE(*,'(A\)')creturn
    !write( * , '(a,a)' , advance='no' ) creturn ,"PROGESS  : [##################..] 90%"
    jm=27
    write( * , '(5a,f6.2,1a)',advance="no")'PROGRESS:' , '[' , & !// 如您的编译器不支持，请用上方语句代替
    repeat('#' , jm ) , repeat( '.' , 30-jm ) , '] ' , 90.00 , "%" 

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

	NL=NTPN !+NBTPOINT
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
        ALLOCATE(MAL(NL+500,NL+500),MAL0(NL+500,NL+500),VERR(NL+500))	
    END IF
    ALLOCATE(NYI(NL,1),NRI(NL,1),NRT(1,NL),VERR0(NL),RTH(NL))
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
    IF(IP+K1.LE.M.OR.IP+K1.GT.M+L) THEN
        DO J=1,NL
            IF(IP+K1.EQ.J) THEN
            !IF(IP+K1.EQ.J.AND.IDS(IP+K1).NE.1) THEN
		        MAL(IP+K1,J)=4.*PI*SE(IP+K1)+SUMHIJ(IP+K1)
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
    IF(IP+K1.GT.M.AND.IP+K1.LE.M+L) THEN
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
    DO I=1,NL
        DO J=1,NL 
    !        QG(I)=QG(I)+MAL1(I,J)*VERR(J)
            !IF(ABS(MAL(I,J)).GT.1000) THEN
            !WRITE(*,*)MAL(I,J)
            !END IF
        END DO
        !WRITE(*,*)I,MAL(I,I)
    END DO
    !CALL MATRIX_PREGMRES
    
    !IF(ITTE.LE.20) THEN
    !IF(ITTE.EQ.1.AND.FR.GE.0.15) THEN
    
    !IF(ITTE.EQ.1.AND.FR.GE.0.12) THEN
    
    IF(ITTE.EQ.1) THEN
    !IF(ITTE.LE.2) THEN
    
    !IF(ITTE.EQ.1.OR.ITTE.EQ.2) THEN
    !IF(ITTE.EQ.100) THEN
    IF(MTROM.EQ.1) CALL TRANSOM_MATRIX
    !NMARK=1
    !CALL MATRIX_PRE0GMRES(NMARK)
    
    !WRITE(*,*)"CALL MATLAB"
    !CALL MATRIX_MATLAB
    !CALL FGMRES_MKL
    CALL LU_MKL
    !CALL BLLU

    ELSE
    
            !IF(ITTE.GE.4) THEN
            ALPHA=1.
            !IF(NTPN.LE.4000) THEN
            IF(FR.LE.0.35) THEN
            IF(ITTE.LE.4) THEN
                ALPHA=1./4.*(ITTE)
            END IF
            ELSE
            IF(ITTE.LE.6) THEN
                ALPHA=1./6.*(ITTE)
            END IF
            END IF
            
            SP0N_1(1:NFPOINT)=SP0N(1:NFPOINT)
            SPX0N_1(1:NFPOINT)=SPX0N(1:NFPOINT)
            SPY0N_1(1:NFPOINT)=SPY0N(1:NFPOINT)
            SPZ0N_1(1:NFPOINT)=SPZ0N(1:NFPOINT)
            SPXX0N_1(1:NFPOINT)=SPXX0N(1:NFPOINT)
            SPXY0N_1(1:NFPOINT)=SPXY0N(1:NFPOINT)
            SPXZ0N_1(1:NFPOINT)=SPXZ0N(1:NFPOINT)
            SPYY0N_1(1:NFPOINT)=SPYY0N(1:NFPOINT)
            SPYZ0N_1(1:NFPOINT)=SPYZ0N(1:NFPOINT)

            IF(ITTE.EQ.2) THEN
                SP0N_1=0.;SPX0N_1=0.;SPY0N_1=0.;SPZ0N_1=0.
                SPXX0N_1=0.;SPXY0N_1=0.;SPXZ0N_1=0.;SPYY0N_1=0.;SPYZ0N_1=0.
            END IF

            SP0N(1:NFPOINT)=SP(NBPOINT1+1:NBPOINT1+NFPOINT)
            SPX0N(1:NFPOINT)=SPX(NBPOINT1+1:NBPOINT1+NFPOINT)
            SPY0N(1:NFPOINT)=SPY(NBPOINT1+1:NBPOINT1+NFPOINT)
            SPZ0N(1:NFPOINT)=SPZ(NBPOINT1+1:NBPOINT1+NFPOINT)
            SPXX0N(1:NFPOINT)=SPXX(NBPOINT1+1:NBPOINT1+NFPOINT)
            SPXY0N(1:NFPOINT)=SPXY(NBPOINT1+1:NBPOINT1+NFPOINT)
            SPXZ0N(1:NFPOINT)=SPXZ(NBPOINT1+1:NBPOINT1+NFPOINT)
            SPYY0N(1:NFPOINT)=SPYY(NBPOINT1+1:NBPOINT1+NFPOINT)
            SPYZ0N(1:NFPOINT)=SPYZ(NBPOINT1+1:NBPOINT1+NFPOINT) 

            SP0(1:NFPOINT)=SP0N_1(1:NFPOINT)+(SP0N(1:NFPOINT)-SP0N_1(1:NFPOINT))*NSWIN(16,1)
            SPX0(1:NFPOINT)=SPX0N_1(1:NFPOINT)+(SPX0N(1:NFPOINT)-SPX0N_1(1:NFPOINT))*NSWIN(16,1)
            SPY0(1:NFPOINT)=SPY0N_1(1:NFPOINT)+(SPY0N(1:NFPOINT)-SPY0N_1(1:NFPOINT))*NSWIN(16,1)
            SPZ0(1:NFPOINT)=SPZ0N_1(1:NFPOINT)+(SPZ0N(1:NFPOINT)-SPZ0N_1(1:NFPOINT))*NSWIN(16,1)
            SPXX0(1:NFPOINT)=SPXX0N_1(1:NFPOINT)+(SPXX0N(1:NFPOINT)-SPXX0N_1(1:NFPOINT))*NSWIN(16,1)
            SPXY0(1:NFPOINT)=SPXY0N_1(1:NFPOINT)+(SPXY0N(1:NFPOINT)-SPXY0N_1(1:NFPOINT))*NSWIN(16,1)
            SPXZ0(1:NFPOINT)=SPXZ0N_1(1:NFPOINT)+(SPXZ0N(1:NFPOINT)-SPXZ0N_1(1:NFPOINT))*NSWIN(16,1)
            SPYY0(1:NFPOINT)=SPYY0N_1(1:NFPOINT)+(SPYY0N(1:NFPOINT)-SPYY0N_1(1:NFPOINT))*NSWIN(16,1)
            SPYZ0(1:NFPOINT)=SPYZ0N_1(1:NFPOINT)+(SPYZ0N(1:NFPOINT)-SPYZ0N_1(1:NFPOINT))*NSWIN(16,1)
      
            
            !END IF

    !IF(ITTE.LE.4) THEN
        MAL0=MAL
        !CALL BRINV(MAL0(1:NTPN,1:NTPN),NTPN)
        DO I=1,NL
        !WRITE(*,*)SUM(MAL0(I,1:NL))
        END DO
    !END IF
    SORW=0.6
    !SORW=0.5
    !SORW=1.
    QG=0.
    DO I=1,NL
        DO J=1,NL 
            QG(I)=QG(I)+SORW*MAL0(I,J)*VERR(J)
        END DO
    END DO
    END IF

    IF(NIT.EQ.1.AND.ITTE.EQ.1) THEN
        ALLOCATE(SP(NL+500),SPN(NL+500),SPX(NL+500),SPY(NL+500),SPZ(NL+500),SPXX(NL+500),SPXY(NL+500),SPXZ(NL+500),SPYY(NL+500),SPYZ(NL+500))
        ALLOCATE(SP0(NL+500),SPX0(NL+500),SPY0(NL+500),SPZ0(NL+500),SPXX0(NL+500),SPXY0(NL+500),SPXZ0(NL+500),SPYY0(NL+500),SPYZ0(NL+500))
        ALLOCATE(SP0N(NL+500),SPX0N(NL+500),SPY0N(NL+500),SPZ0N(NL+500),SPXX0N(NL+500),SPXY0N(NL+500),SPXZ0N(NL+500),SPYY0N(NL+500),SPYZ0N(NL+500))
        ALLOCATE(SP0N_1(NL+500),SPX0N_1(NL+500),SPY0N_1(NL+500),SPZ0N_1(NL+500),SPXX0N_1(NL+500),SPXY0N_1(NL+500),SPXZ0N_1(NL+500),SPYY0N_1(NL+500),SPYZ0N_1(NL+500))
        ALLOCATE(SIGM0(NL+500),SIGM(NL+500))
    END IF

    IF(ITTE.LE.1) THEN
    !IF(ITTE.LE.2) THEN
            SP0=0
            SPX0=0
            SPY0=0
            SPZ0=0
            SPXX0=0
            SPXY0=0
            SPXZ0=0
            SPYY0=0
            SPYZ0=0
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
    DO IP=M+1,M+L
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

    DO IP=M+L+1,NL
        DO J=1,NL
            SP(IP)=SP(IP)+SA(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            IF(IP.NE.J) THEN
                SPN(IP)=SPN(IP)+SH(IP,J)*SIGM(J)    !求解I点速度势沿法向导数（物面）
            ELSE
                SPN(IP)=SPN(IP)+4.*PI*SE(IP)*SIGM(J)    !求解I点速度势沿法向导数（物面）
            END IF
        END DO
    END DO    
    
    !GOTO 898
    IF(FR.GE.NSWIN(12,1)) THEN
        ALLOCATE(SMPHI(NTPN))    
        SMPHI(M+1:M+L)=SPX(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPX(M+1:NL)=SMPHI(M+1:M+L)
        SMPHI(M+1:M+L)=SPY(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPY(M+1:NL)=SMPHI(M+1:M+L)
        SMPHI(M+1:M+L)=SPZ(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPZ(M+1:NL)=SMPHI(M+1:M+L)
        SMPHI(M+1:M+L)=SPXX(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPXX(M+1:M+L)=SMPHI(M+1:M+L)
        SMPHI(M+1:M+L)=SPXY(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPXY(M+1:M+L)=SMPHI(M+1:M+L)
        SMPHI(M+1:M+L)=SPXZ(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPXZ(M+1:M+L)=SMPHI(M+1:M+L)
        DEALLOCATE(SMPHI) 
    END IF
    898 CONTINUE    
    
    IF(ITTE.EQ.1) GOTO 919
    !IF(ITTE.LE.2) GOTO 919
    !IF(ITTE.LE.20) GOTO 919
    !GOTO 919

    !OPEN(2,FILE='PHIS.DAT')
    !DO I=1,M
        !WRITE(2,*)I,SP(I),SIGM(I)
    !END DO
    !CLOSE(2)
    
    !WRITE(*,*)"VERR..........OK"

    MAL=0.
    DO IP=1,M
        DO J=1,NL
            IF(IP.EQ.J) THEN
            !IF(IP+K1.EQ.J.AND.IDS(IP+K1).NE.1) THEN
		        MAL(IP,J)=4.*PI*SE(IP)+SUMHIJ(IP)
                !MAL(IP+K1,J)=2.*PI
            ELSE
                MAL(IP,J)=SH(IP,J)
            END IF
        END DO  
    END DO

    DO IP=M+1,M+L
        DO J=1,NL
            !MAL(IP,J)=-1./2./G*(2.*U*SX(IP,J)*SPXX0(IP)-2.*SX(IP,J)*SPX0(IP)*SPXX0(IP)-2.*SX(IP,J)*SPY0(IP)*SPXY0(IP)-2.*SX(IP,J)*SPZ0(IP)*SPXZ0(IP))&
            !          -1./2./G*(2.*U*SY(IP,J)*SPXY0(IP)-2.*SY(IP,J)*SPX0(IP)*SPXY0(IP)-2.*SY(IP,J)*SPY0(IP)*SPYY0(IP)-2.*SY(IP,J)*SPZ0(IP)*SPYZ0(IP))&
            !          +SZ(IP,J)&
            !          +1./2./G*(2.*U**2*SXX(IP,J))&
            !          +1./2./G*(-2.*U*SX(IP,J)*SPXX0(IP)-2.*U*SY(IP,J)*SPXY0(IP)-2.*U*SZ(IP,J)*SPXZ0(IP))
            
            !MAL(IP,J)=-1./2./G*(2.*U*SPX0(IP)*SXX(IP,J))&
            !          -1./2./G*(2.*U*SPY0(IP)*SXY(IP,J))&
            !          +SZ(IP,J)&
            !          +1./2./G*(2.*U**2*SXX(IP,J)) !& 
            
            
            MAL(IP,J)=1./2./G*(2.*U**2*SXX(IP,J))+SZ(IP,J) !&
                     !-1./2./G*(4.*U*SPX0(IP-M)*SXX(IP,J))
            
            !MAL(IP,J)=U**2/G*SXX(IP,J)+SZ(IP,J)
        END DO
    END DO

    DO IP=M+L+1,NL
        DO J=1,NL
            IF(IP.EQ.J) THEN
            !IF(IP+K1.EQ.J.AND.IDS(IP+K1).NE.1) THEN
		        MAL(IP,J)=4.*PI*SE(IP)+SUMHIJ(IP)
                !MAL(IP+K1,J)=2.*PI
            ELSE
                MAL(IP,J)=SH(IP,J)
            END IF
        END DO  
    END DO    
    
    DO IP=1,M
        !VERR(IP)=(U*VECN(IP,1)-SPN(IP))
        VERR(IP)=U*VECN(IP,1)
    END DO
    
    DO IP=M+1,M+L
        !VERR(IP)=-1./2./G*(-2.*U*SPX0(IP)*SPXX0(IP)-2.*U*SPY0(IP)*SPXY0(IP)-2.*U*SPZ0(IP)*SPXZ0(IP))&
        !         +1./2./G*(-2.*SPX0(IP)*SPX0(IP)*SPXX0(IP)-2.*SPX0(IP)*SPY0(IP)*SPXY0(IP)-2.*SPX0(IP)*SPZ0(IP)*SPXZ0(IP))&
        !         +1./2./G*(-2.*SPY0(IP)*SPX0(IP)*SPXY0(IP)-2.*SPY0(IP)*SPY0(IP)*SPYY0(IP)-2.*SPY0(IP)*SPZ0(IP)*SPYZ0(IP))

        !VERR(IP)=-1./2./G*(-2.*U*SPX0(IP)*SPXX0(IP)-2.*U*SPY0(IP)*SPXY0(IP)-2.*U*SPZ0(IP)*SPXZ0(IP))&
        !+1./2./G*(2.*U*SPX0(IP)*SPXX0(IP)-2.*SPX0(IP)*SPX0(IP)*SPXX0(IP)-2.*SPX0(IP)*SPY0(IP)*SPXY0(IP)-2.*SPX0(IP)*SPZ0(IP)*SPXZ0(IP))&
        !+1./2./G*(2.*U*SPY0(IP)*SPXY0(IP)-2.*SPY0(IP)*SPX0(IP)*SPXY0(IP)-2.*SPY0(IP)*SPY0(IP)*SPYY0(IP)-2.*SPY0(IP)*SPZ0(IP)*SPYZ0(IP))
        
        !VERR(IP)=0.
    END DO

    !检查这里
    
    !IF(NSWIN(13,1)/PANIN(1,3).GT.1.8) THEN
    IF((U/SQRT(G*NSWIN(13,1))).LT.0.942) THEN
    !IF(NSWIN(13,1)/PANIN(1,3).GT.1.0) THEN
    WRITE(*,*)"NORMAL"
    DO IP=1,NFPOINT
        
        VERR(IP+M)=-1./2./G*(-2.*U*SPX0(IP)*SPXX0(IP)-2.*U*SPY0(IP)*SPXY0(IP)-2.*U*SPZ0(IP)*SPXZ0(IP))&
        +1./2./G*(2.*U*SPX0(IP)*SPXX0(IP)-2.*SPX0(IP)*SPX0(IP)*SPXX0(IP)-2.*SPX0(IP)*SPY0(IP)*SPXY0(IP)-2.*SPX0(IP)*SPZ0(IP)*SPXZ0(IP))&
        +1./2./G*(2.*U*SPY0(IP)*SPXY0(IP)-2.*SPY0(IP)*SPX0(IP)*SPXY0(IP)-2.*SPY0(IP)*SPY0(IP)*SPYY0(IP)-2.*SPY0(IP)*SPZ0(IP)*SPYZ0(IP))
         
        !VERR(IP+M)=0.
    END DO
    ELSE
    DO IP=1,NFPOINT    
        VERR(IP+M)=-1./2./G*(-2.*U*SPY0(IP)*SPXY0(IP)-2.*U*SPZ0(IP)*SPXZ0(IP))&
        +1./2./G*(-2.*SPX0(IP)*SPX0(IP)*SPXX0(IP)-2.*SPX0(IP)*SPY0(IP)*SPXY0(IP)-2.*SPX0(IP)*SPZ0(IP)*SPXZ0(IP))&
        +1./2./G*(2.*U*SPY0(IP)*SPXY0(IP)-2.*SPY0(IP)*SPX0(IP)*SPXY0(IP)-2.*SPY0(IP)*SPY0(IP)*SPYY0(IP)-2.*SPY0(IP)*SPZ0(IP)*SPYZ0(IP))
        
        !VERR(IP+M)=-1./2./G*(-2.*U*SPX0(IP)*SPXX0(IP))&
        !           +1./2./G*(2.*U*SPY0(IP)*SPXY0(IP))
                
        VERR(IP+M)=0.
    END DO
    END IF
    
    IF(MDFREE.EQ.3) VERR(M+1:M+L)=0.  !LINEAR FREESURFACE 
    
    DO IP=M+L+1,NL
        !VERR(IP)=(U*VECN(IP,1)-SPN(IP))
        VERR(IP)=0.
    END DO    
    
    DO IP=1,NFPOINT
        !VERR(IP+M)=-1./2./G*(-2.*U*(SPX0N(IP)*SPXX0N(IP)+SPX0N_1(IP)*SPXX0N_1(IP))/2.&
        !                     -2.*U*(SPY0N(IP)*SPXY0N(IP)+SPY0N_1(IP)*SPXY0N_1(IP))/2.&
        !                     -2.*U*(SPZ0N(IP)*SPXZ0N(IP)+SPZ0N_1(IP)*SPXZ0N_1(IP))/2.)&
        !            +1./2./G*(2.*U*(SPX0N(IP)*SPXX0N(IP)+SPX0N_1(IP)*SPXX0N_1(IP))/2.&
        !                     -2.*(SPX0N(IP)*SPX0N(IP)*SPXX0N(IP)+SPX0N_1(IP)*SPX0N_1(IP)*SPXX0N_1(IP))/2.&
        !                     -2.*(SPX0N(IP)*SPY0N(IP)*SPXY0N(IP)+SPX0N_1(IP)*SPY0N_1(IP)*SPXY0N_1(IP))/2.&
        !                     -2.*(SPX0N(IP)*SPZ0N(IP)*SPXZ0N(IP)+SPX0N_1(IP)*SPZ0N_1(IP)*SPXZ0N_1(IP))/2.)&
        !            +1./2./G*(2.*U*(SPY0N(IP)*SPXY0N(IP)+SPY0N(IP)*SPXY0N_1(IP))/2.&
        !                     -2.*(SPY0N(IP)*SPX0N(IP)*SPXY0N(IP)+SPY0N_1(IP)*SPX0N_1(IP)*SPXY0N_1(IP))/2.&
        !                     -2.*(SPY0N(IP)*SPY0N(IP)*SPYY0N(IP)+SPY0N_1(IP)*SPY0N_1(IP)*SPYY0N_1(IP))/2.&
        !                     -2.*(SPY0N(IP)*SPZ0N(IP)*SPYZ0N(IP)+SPY0N_1(IP)*SPZ0N_1(IP)*SPYZ0N_1(IP))/2.)
    END DO

    !VERR0=VERR

    !GOTO 909

    !WRITE(*,*)"MAL..........OK"

        
    !909 CONTINUE

    !IF(ITTE.LE.4) THEN
        !MAL0=MAL
        !CALL BRINV(MAL0(1:NL,1:NL),NL)
    !END IF
    !SORW=0.5
    !OPEN(25,FILE='TEST.DAT')
    QG=0.
    DO I=1,NL
        DO J=1,NL 
            !QG(I)=QG(I)+MAL0(I,J)*VERR(J)
            !WRITE(25,*)J,NBPOINT,MAL0(I,J)
        END DO
        !WRITE(25,*)QG(I)
        !STOP
        !WRITE(*,*)MAL0(I,I+1),VERR(I),QG(I)
    END DO
    !CALL MATRIX_PREGMRES
    !NMARK=3
    !CALL MATRIX_PRE0GMRES(NMARK)
    
    
    IF(MTROM.EQ.1) CALL TRANSOM_MATRIX

    !NMARK=1
    !CALL MATRIX_PRE0GMRES(NMARK)
    !WRITE(*,*)"CALL MATLAB"
    !CALL MATRIX_MATLAB
    !CALL FGMRES_MKL
    CALL LU_MKL
    !CALL BLLU
    SIGM0=SIGM
    SIGM=QG


    909 CONTINUE
    !GOTO 919
!****************************************************************************************
    NITER=ITTE*4 

    !END DO

    NITER=300
    MAL=MAL0
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
    DO IP=M+1,M+L
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

    DO IP=M+L+1,NL
        DO J=1,NL
            SP(IP)=SP(IP)+SA(IP,J)*SIGM(J)    !求解I点速度势沿X向的导数（自由面）
            IF(IP.NE.J) THEN
                SPN(IP)=SPN(IP)+SH(IP,J)*SIGM(J)    !求解I点速度势沿法向导数（物面）
            ELSE
                SPN(IP)=SPN(IP)+4.*PI*SE(IP)*SIGM(J)    !求解I点速度势沿法向导数（物面）
            END IF
        END DO
    END DO    
    
    IF(FR.GE.NSWIN(12,1)) THEN
        ALLOCATE(SMPHI(NTPN))    
        SMPHI(M+1:M+L)=SPX(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPX(M+1:M+L)=SMPHI(M+1:M+L)
        SMPHI(M+1:M+L)=SPY(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPY(M+1:M+L)=SMPHI(M+1:M+L)
        SMPHI(M+1:M+L)=SPZ(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPZ(M+1:M+L)=SMPHI(M+1:M+L)
        SMPHI(M+1:M+L)=SPXX(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPXX(M+1:M+L)=SMPHI(M+1:M+L)
        SMPHI(M+1:M+L)=SPXY(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPXY(M+1:M+L)=SMPHI(M+1:M+L)
        SMPHI(M+1:M+L)=SPXZ(M+1:M+L)
        CALL SMOOTHING_RECT_PHI
        SPXZ(M+1:M+L)=SMPHI(M+1:M+L)
        DEALLOCATE(SMPHI) 
    END IF
    
    DO IP=1,M
        VERR(IP)=(U*VECN(IP,1)-SPN(IP))
        !VERR(IP)=-U*VECN(IP,1)
    END DO
    
    DO IP=M+1,M+L
        VERR(IP)=(1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPX(IP)*SPXX(IP)-2.*SPX(IP)*SPY(IP)*SPXY(IP)-2.*SPX(IP)*SPZ(IP)*SPXZ(IP))&
        -1./2./G*(2.*U*SPY(IP)*SPXY(IP)-2.*SPY(IP)*SPX(IP)*SPXY(IP)-2.*SPY(IP)*SPY(IP)*SPYY(IP)-2.*SPY(IP)*SPZ(IP)*SPYZ(IP))&
        +SPZ(IP))
        
        
        !VERR(IP)=(1./2./G*(2.*U**2*SPXX(IP)-2.*U*SPX(IP)*SPXX(IP)-2.*U*SPY(IP)*SPXY(IP)-2.*U*SPZ(IP)*SPXZ(IP))&
        !-1./2./G*(2.*U**2*SPX(IP)*SPXX(IP))&
        !-1./2./G*(2.*U**2*SPY(IP)*SPXY(IP))&
        !+SPZ(IP))      
        !VERR(IP)=0.
    END DO

    DO IP=M+L+1,NL
        VERR(IP)=(0.-SPN(IP))
        !VERR(IP)=-U*VECN(IP,1)
    END DO    
    
    NORM1=0.
    NORM2=0
    DO IP=1,NL
        NORM1=NORM1+(SIGM0(IP)-SIGM(IP))**2
        NORM2=NORM2+SIGM(IP)**2
        !NORM2=NORM2+VERR(IP)**2
    END DO
    NORMB(1)=MAXVAL(ABS(VERR(1:NBPOINTW)))
    !NORMB(2)=MAXVAL(ABS(VERR(NBPOINT+1:NL)))

    IF(MTROM.EQ.1) THEN
        NORMB(2)=MAXVAL(ABS(VERR(NBPOINT+1:NL-NZS*(NFXF+1))))
    ELSE
        !NORMB(2)=MAXVAL(ABS(VERR(NBPOINT+1:NTPN)))
        NORMB(2)=MAXVAL(ABS(VERR(NBPOINT+1:M+L)))
    END IF

    NORM1=MAXVAL(NORMB)
    NEWTON_RESIDU=NORMB(2) !NORM1
    !WRITE(*,*)"NEWTON RESIDUAL",ITER,NORM1
    !WRITE(117,*)"NEWTON RESIDUAL",ITTE,NORM1
    !IF(NORM1.LE.1.E-5) EXIT
    !STOP

    GOTO 919

    WRITE(*,*)"VERR..........OK"

        NYI=0.
        NRT=0.
        NRI=0.

        NYI(1:NL,1)=VERR(1:NL)-VERR0(1:NL)
        NRT(1,1:NL)=SIGM(1:NL)-SIGM0(1:NL)
        NRI(1:NL,1)=SIGM(1:NL)-SIGM0(1:NL)
        RTR=0.
        RTH=0.
        DO I=1,NL
            DO J=1,NL
                RTH(I)=RTH(I)+NRT(1,J)*MAL(J,I)
            END DO
        END DO
        DO I=1,NL
            RTR=RTR+RTH(I)*NYI(I,1)
        END DO

        IF(MOD(ITER,120).EQ.0) THEN
            MAL=MAL+MATMUL(MATMUL(NRI-MATMUL(MAL,NYI),NRT),MAL)/RTR !MATMUL(MATMUL(NRT,MAL),NYI)
        END IF

    !IF(ITER.GE.30) SORW=0.5
    SIGM0=SIGM
    VERR0=VERR
    QG=0
    DO I=1,NL
        DO J=1,NL 
            QG(I)=QG(I)+SORW*MAL(I,J)*VERR(J)
        END DO
    END DO
    SIGM=SIGM0-QG

    NORM1=0.
    NORM2=0.
    DO IP=1,NL
        NORM1=NORM1+(SIGM0(IP)-SIGM(IP))**2
        NORM2=NORM2+SIGM(IP)**2
        !NORM2=NORM2+VERR(IP)**2
    END DO
    NORM1=SQRT(NORM1)/SQRT(NORM2)
    !NORM1=SQRT(NORM2)/NL
    NORMB(1)=MAXVAL(ABS(VERR(1:NBPOINTW)))
    IF(MTROM.EQ.1) THEN
        NORMB(2)=MAXVAL(ABS(VERR(NBPOINT+1:NL-NZS*(NFXF+1))))
    ELSE
        NORMB(2)=MAXVAL(ABS(VERR(NBPOINT+1:NL)))
    END IF
    NORM1=MAXVAL(NORMB)
    !WRITE(*,*)"NEWTON RESIDUAL",ITER,NORM1
    IF(NORM1.LE.1.E-6) EXIT

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
	!DEALLOCATE(SXY,SXZ,SYY,SYZ)
    DEALLOCATE(SXX,SXY,SXZ,SYY,SYZ)
       
    
    
    !WRITE(*,*)'		END MATRIX'

	RETURN
END SUBROUTINE

SUBROUTINE TRANSOM_MATRIX     !计算线性问题
!**************************************************************
!CGENERATE THE REAL MATRIX OF THE INTEGRATION EQUATION
!**************************************************************
	USE GREENMOD
	USE CUMOD

    REAL:: ETATR(20)

    OF=0.55
    !OF=0.5

    CJ=-(1.-OF**2)/(OF+OF**2)/DFX
    CJ_1=-OF**2/(OF+OF**2)/DFX
    CTR=1./(OF+OF**2)/DFX

    CALL ESPL2(XTR(1:NTAO),ZTR(1:NTAO),NTAO,CORD(NTPN-NZS+1:NTPN,2),NZS,ETATR(1:NZS))
    DO I=1,NZS
        IF(CORD(NTPN-NZS+I,1).GE.XTR(NTAO)) ETATR(I)=ZTR(NTAO)
        !ETATR(I)=ZTR(1)
    END DO

    GOTO 10
    DO I=NTPN-NZS+1,NTPN
        DO J=1,NTPN
            MAL(I,J)=U/G*SX(I,J)
        END DO
        !VERR(I)=ETATR(I-NTPN+NZS)
    END DO
    10 CONTINUE

    GOTO 20
    DO I=NTPN-NZS*2+1,NTPN-NZS
        DO J=1,NTPN
            !MAL(I,J)=U/G*(SXX(I,J)-CJ*SX(I,J)-CJ_1*SX(I-NZS,J))
            MAL(I,J)=U/G*(SXX(I,J)-2.*SX(I,J)/DFX)
        END DO
        !VERR(I)=CTR*ETATR(I-NTPN+NZS*2)
        VERR(I)=-(2.*ETATR(I-NTPN+NZS*2)/DFX+TAOX(1))
    END DO
    20 CONTINUE

    !GOTO 30
    DO I=NTPN-NZS+1,NTPN
        DO J=1,NTPN
            MAL(I,J)=U/G*(SXX(I,J)-CJ*SX(I,J)-CJ_1*SX(I-NZS,J))
            !MAL(I,J)=U/G*(SXX(I,J)-2.*SX(I,J)/DFX)
        END DO
        VERR(I)=CTR*ETATR(I-NTPN+NZS)
        !VERR(I)=-(2.*ETATR(I-NTPN+NZS)/DFX+TAOX(1))
    END DO
    30 CONTINUE

    GOTO 40
    I=NBPOINT+(NFXF)*NFY+1
        DO J=1,NTPN
            MAL(I,J)=U/G*(SXX(I,J)-CJ*SX(I,J)-CJ_1*SX(I-NFY,J))
        END DO
        VERR(I)=CTR*ETATR(NZS)
    40 CONTINUE

END SUBROUTINE