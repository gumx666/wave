SUBROUTINE MATRIX_DAWSOM  !计算DAWSON问题
!**************************************************************
!CGENERATE THE REAL MATRIX OF THE INTEGRATION EQUATION
!**************************************************************
	USE GREENMOD
	USE CUMOD

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

    REAL,ALLOCATABLE:: DL(:)
    
    ALLOCATE(DL(NBPOINT))

    !CALL PADCOEF
    !CALL DIFFCORD_DECOMP   ! 计算坐标对计算平面导数
    !CALL DIFFCORD_MATRIX_DECOMP ! 计算差分算子


!采用bin文件格式可以减小文件大小，加快读写速度 
!读入对角线元素        
    !IF(NIT.EQ.1) THEN
    OPEN(115,FILE='TESTMATRIXSUM.DAT')
    DO I=1,NBPOINT    
        READ(115,*)J,DL(I)
    END DO
    CLOSE(115)
    !END IF

    open(16,file='coefuse.bin',form='binary',access='sequential',&
          status='old')
    rewind(16)
	WRITE(*,*)
	WRITE(*,*)'SUB MATRIX.........'
	!ALLOCATE(SH(NTPN,NTPN),SX(NTPN,NTPN),SZ(NTPN,NTPN),SA(NTPN,NTPN),SA1(NTPN,NTPN),SY(NTPN,NTPN))

	DO 10 I=1,NTPN
		DO 10 J=1,NTPN
		!READ(16)SA(I,J),SX(I,J),SH(I,J),SZ(I,J),SY(I,J)   !SA1船体表面诱导系数，SA自由表面诱导系数 
!CREAD(16,*)SH(I,J),SA(I,J)
10	CONTINUE

	NL=NTPN
    M=NBPOINT
    L=NFPOINT

	WRITE(*,*)'    NL=     ',NL,'      NR=     ',NR

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

	ALLOCATE(MAL(NL,NL),MAL1(NL,NL),VERR(NL),MAR(NL,NL))	
    MAL=0.
    MAR=0.

!START PARALLEL
    CALL OMP_SET_NUM_THREADS(NTHREAD)

!$omp parallel private(CP,k1,k2,IP,J,K)
    
    CP=OMP_GET_THREAD_NUM()+1

    K1=SUM(NFS(1:CP))-NFS(CP) 
    K2=SUM(NFS(1:CP))

    DO IP=1,NFS(CP)

	!DO I=1,M             !第一部分物面
	!	VERR(I)=0. !U*VECN(I,1)
    !END DO
    IF(K1+IP.LE.M) VERR(IP+K1)=0.

    !DO I=M+1,NL
    !    DO J=M+1,NL
    !        MAR(I,1)=MAR(I,1)+DMX(I,J)*((DPHIS_X(J)-U)**2+DPHIS_Y(J)**2)
    !        MAR(I,2)=MAR(I,2)+DMY(I,J)*((DPHIS_X(J)-U)**2+DPHIS_Y(J)**2)
    !    END DO
    !END DO

    !DO I=M+1,NL           !第二部分自由面
	!	VERR(I)=-0.5*((DPHIS_X(I)-U)*MAR(I,1)+DPHIS_Y(I)*MAR(I,2))
    !END DO
    IF(K1+IP.GT.M) THEN
        DO J=M+1,M+L
            MAR(IP+K1,1)=MAR(IP+K1,1)+DMX(IP+K1,J)*((DPHIS_X(J)-U)**2+DPHIS_Y(J)**2)
            MAR(IP+K1,2)=MAR(IP+K1,2)+DMY(IP+K1,J)*((DPHIS_X(J)-U)**2+DPHIS_Y(J)**2)
        END DO
        VERR(IP+K1)=-0.5*((DPHIS_X(IP+K1)-U)*MAR(IP+K1,1)+DPHIS_Y(IP+K1)*MAR(IP+K1,2)) 
    END IF

    !DO I=1,M            !第一部分物面
	!    DO J=1,NL
    !        IF(I.EQ.J) THEN
		        !MAL(I,J)=DL(I) !4.*PI*SE(I) !+SUMHIJ(I)
	!	        MAL(I,J)=4.*PI*SE(I) !+SUMHIJ(I)
                !WRITE(*,*)MAL(I,J)
    !        ELSE
    !            MAL(I,J)=SH(I,J)
    !        END IF
    !    END DO
    !END DO
    IF(K1+IP.LE.M) THEN
	    DO J=1,NL
            !IF(IP+K1.EQ.J.AND.ABS(CORD(K1+IP,3)).GE.1.E-4) THEN
            IF(IP+K1.EQ.J) THEN
		        MAL(IP+K1,J)=4.*PI*SE(IP+K1) !+SUMHIJ(I)
                !WRITE(*,*)MAL(I,J)
            ELSE
                MAL(IP+K1,J)=SH(IP+K1,J)
            END IF
        END DO
    END IF

    !MAL1(K1+IP:K2,:)=0
    !DO I=M+1,NL           !第二部分自由面
	!    DO J=1,NL
    !        MAL(I,J)=G*SZ(I,J)
!
    !        MAL1(I,J)=0.
    !        DO K=M+1,NL
    !            MAL1(I,J)=MAL1(I,J)+DMX(I,K)*SX(K,J)
    !        END DO
    !        MAL(I,J)=MAL(I,J)+(DPHIS_X(I)-U)**2*MAL1(I,J)
!
    !        MAL1(I,J)=0.
    !        DO K=M+1,NL
    !            MAL1(I,J)=MAL1(I,J)+DMX(I,K)*SY(K,J)
    !        END DO
    !        MAL(I,J)=MAL(I,J)+2.*(DPHIS_X(I)-U)*DPHIS_Y(I)*MAL1(I,J)
!
    !        MAL1(I,J)=0.
    !        DO K=M+1,NL
    !            MAL1(I,J)=MAL1(I,J)+DMY(I,K)*SY(K,J)
    !        END DO
    !        MAL(I,J)=MAL(I,J)+DPHIS_Y(I)**2*MAL1(I,J)

    !        MAL(I,J)=MAL(I,J)+0.5*SX(I,J)*MAR(I,1)
    !        MAL(I,J)=MAL(I,J)+0.5*SY(I,J)*MAR(I,2)
    !    END DO
    !END DO
    IF(K1+IP.GT.M) THEN
	    DO J=1,NL
            MAL(IP+K1,J)=G*SZ(IP+K1,J)
!
            MAL1(IP+K1,J)=0.
            DO K=M+1,NL
                MAL1(IP+K1,J)=MAL1(IP+K1,J)+DMX(IP+K1,K)*SX(K,J)
            END DO
            MAL(IP+K1,J)=MAL(IP+K1,J)+(DPHIS_X(IP+K1)-U)**2*MAL1(IP+K1,J)
!
            MAL1(IP+K1,J)=0.
            DO K=M+1,NL
                MAL1(IP+K1,J)=MAL1(IP+K1,J)+DMX(IP+K1,K)*SY(K,J)
            END DO
            MAL(IP+K1,J)=MAL(IP+K1,J)+2.*(DPHIS_X(IP+K1)-U)*DPHIS_Y(IP+K1)*MAL1(IP+K1,J)
!
            MAL1(IP+K1,J)=0.
            DO K=M+1,NL
                MAL1(IP+K1,J)=MAL1(IP+K1,J)+DMY(IP+K1,K)*SY(K,J)
            END DO
            MAL(IP+K1,J)=MAL(IP+K1,J)+DPHIS_Y(IP+K1)**2*MAL1(IP+K1,J)
            
            MAL(IP+K1,J)=MAL(IP+K1,J)+0.5*SX(IP+K1,J)*MAR(IP+K1,1)
            MAL(IP+K1,J)=MAL(IP+K1,J)+0.5*SY(IP+K1,J)*MAR(IP+K1,2)

!           新方程
            MAL(IP+K1,J)=MAL(IP+K1,J)+(DPHIS_X(IP+K1)-U)*DPHIS_XY(IP+K1)*SY(IP+K1,J)
!
            MAL(IP+K1,J)=MAL(IP+K1,J)+DPHIS_Y(IP+K1)*DPHIS_XY(IP+K1)*SX(IP+K1,J)
!
            MAL(IP+K1,J)=MAL(IP+K1,J)+DPHIS_Y(IP+K1)*DPHIS_YY(IP+K1)*SY(IP+K1,J)

        END DO
    END IF
    END DO
!$omp end parallel 
    
    !GOTO 600
    OPEN(991,FILE='CORNERNODS.DAT')
    DO I=1,NFX
        IP=nbpoint+NFXF*NFY+(I-1)*NFY+1
        J=(I-1)*NZ+1
        
        !SX(IP,J)=0.
        !SY(IP,J)=0.
        !SZ(IP,J)=0.
        !DO K=M+1,NL
        !    SX(IP,J)=SX(IP,J)+DMX(IP,K)*SA(K,J)
        !    SY(IP,J)=SY(IP,J)+DMY(IP,K)*SA(K,J)
        !END DO
        
        !SX(IP,J)=SX(IP,J)+4.*PI*SE(J)*VECN(IP,1)
        !SZ(IP,J)=SZ(IP,J)+4.*PI*SE(J)*VECN(IP,3)
    END DO

    DO I=1,NFX
        IP=nbpoint+NFXF*NFY+(I-1)*NFY+1
        J=(I-1)*NZ+1
 
            !MAL(IP,J)=G*SZ(IP,J)
!
            !MAL1(IP,J)=0.
            !DO K=M+1,NL
            !    MAL1(IP,J)=MAL1(IP,J)+DMX(IP,K)*SX(K,J)
            !END DO
            !MAL(IP,J)=MAL(IP,J)+(DPHIS_X(IP)-U)**2*MAL1(IP,J)
!
            !MAL1(IP,J)=0.
            !DO K=M+1,NL
            !    MAL1(IP,J)=MAL1(IP,J)+DMY(IP,K)*SX(K,J)
            !END DO
            !MAL(IP,J)=MAL(IP,J)+2.*(DPHIS_X(IP)-U)*DPHIS_Y(IP)*MAL1(IP,J)
!
            !MAL1(IP,J)=0.
            !DO K=M+1,NL
            !    MAL1(IP,J)=MAL1(IP,J)+DMY(IP,K)*SY(K,J)
            !END DO
            !MAL(IP,J)=MAL(IP,J)+DPHIS_Y(IP)**2*MAL1(IP,J)
            
            !MAL(IP,J)=MAL(IP,J)+0.5*SX(IP,J)*MAR(IP,1)
            !MAL(IP,J)=MAL(IP,J)+0.5*SY(IP,J)*MAR(IP,2)

!           新方程
            !MAL(IP,J)=MAL(IP,J)+(DPHIS_X(IP)-U)*DPHIS_XY(IP)*SY(IP,J)
!
            !MAL(IP,J)=MAL(IP,J)+DPHIS_Y(IP)*DPHIS_XY(IP)*SX(IP,J)
!
            !MAL(IP,J)=MAL(IP,J)+DPHIS_Y(IP)*DPHIS_YY(IP)*SY(IP,J)

        WRITE(991,119)IP,SX(IP,J),SY(IP,J),SZ(IP,J),SZ(IP,IP),MAL(IP,J)
        119FORMAT(I8,6F15.6)
        !WRITE(991,*)CORD(IP,:)
        !WRITE(991,*)CORD(J,:)
        
        !WRITE(991,*)IP,SX(IP,J),SE(IP)*4.*PI
        !WRITE(991,*)IP,SX(IP,J),SX(IP,IP)
        !WRITE(991,*)IP,SX(IP,J),SX(IP,IP)
        
        !WRITE(*,*)CORD(IP,:)
        !WRITE(*,*)CORD(JP,:)
    END DO
    CLOSE(991)
    600 CONTINUE

30  DEALLOCATE(MAR,MAL1,DL,IDS)											
	CLOSE(16)
!释放内存

	!DEALLOCATE(SH,SA,SA1,SX,SZ,SY)
	WRITE(*,*)'		END MATRIX'

	RETURN
END SUBROUTINE

SUBROUTINE MATRIX_DBFLOW   !计算叠模势
	USE GREENMOD
	USE CUMOD

!C 采用BIN文件格式可以减小文件大小，加快读写速度 
    OPEN(16,FILE='COEFUSE_SP.BIN',FORM='BINARY',ACCESS='SEQUENTIAL',& 
         STATUS='OLD')
    REWIND(16)

    M=NBPOINT
    L=NFPOINT

	WRITE(*,*)'SUB MATRIX.........'
	ALLOCATE(SH(M,M),SA(M,M))

	DO 10 I=1,M
		DO 10 J=1,M
		READ(16)SA(I,J),SH(I,J)
        !WRITE(*,*)I,J,SH(I,J),SA(I,J)
   		!READ(18,*)SH(I,J),SA(I,J)
10	CONTINUE

	ALLOCATE(MAR(M,M),VER(M),VERR(M))	
! 产生方程右边的向量
! 前M+N+L个为与物体表面法向导数相关的项
	DO 51 I=1,M

	DO 51 J=1,M
		MAR(I,J)=SA(I,J)
51	CONTINUE

	DO 52 I=1,M
		VER(I)=U*VECN(I,1)
52	CONTINUE

!C	ASSIGN THE INITIAL VALUE OF THE RIGHT VECTOR
	DO 60 I=1,M
		VERR(I)=0.0
60	CONTINUE

	DO 61 I=1,M
		DO 62 J=1,M
			VERR(I)=VERR(I)+MAR(I,J)*VER(J)
62		CONTINUE
		!WRITE(*,*)'VERR(',I,')=',VERR(I)
61	CONTINUE

	DEALLOCATE(MAR,VER)

  ! 形成左边的矩阵	
	ALLOCATE(MAL(M,M))	
	OPEN(22,FILE='TESTMATRIXSUM.DAT')
	DO 22 I=1,M
		TESTSUM=4*PI-SUMHIJ(I)
		WRITE(22,*)I,TESTSUM,SUMHIJ(I)
        !WRITE(*,*)'I=',I,'SUMHIJ=',SUMHIJ(I)
        !WRITE(22,101)I,SUMHIJ(I),SH(I,I),SUMHIJ1(I)/4./PI
22	CONTINUE
	CLOSE(22)

    101FORMAT('I=',I5,'   SUMHIJ=',F15.6,'   TESTSUM=',F15.6,F15.6)
    			
	DO 12 I=1,M
	DO 13 J=1,M
		IF (I.EQ.J)THEN
            MAL(I,J)=4*PI-SUMHIJ(I) !OLD
            !MAL(I,J)=SUMHIJ(I) !OLD
        ELSE 
			MAL(I,J)=SH(I,J)
		END IF
13	CONTINUE
12	CONTINUE

	CLOSE(16)
!C 释放内存
	WRITE(*,*)'TEST'
	DEALLOCATE(SUMHIJ)
	DEALLOCATE(SH,SA)

    CALL EQUATION        !求解叠模流动

	RETURN
END SUBROUTINE

SUBROUTINE MATRIX_HOB      !计算线性问题
!**************************************************************
!CGENERATE THE REAL MATRIX OF THE INTEGRATION EQUATION
!**************************************************************
	USE GREENMOD
	USE CUMOD

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

    REAL,ALLOCATABLE:: DL(:)
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

	20 ALLOCATE(MAR(NL,NL),VER(NL),MAL1(NL,NL))
    ALLOCATE(MAL(NL,NL),VERR(NL))	
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

    IF(MDIFF.EQ.2) THEN
!   用差分求解
    
    WRITE(*,*)"用差分求解"

    ALLOCATE(DVALU(M+L,M+L))
    DVALU=SX
    CALL DIFF_COEMATRIX_DECOMP

    DO I=M+1,NL
	    DO J=1,NL
            MAL(I,J)=U**2/G*DVMX(I,J)+SZ(I,J)
        END DO
    END DO

    DEALLOCATE(DVALU)
    END IF
!   END


    !STOP


30  DEALLOCATE(MAR,VER,MAL1,IDS)											
	CLOSE(16)
!释放内存

    !NBPOINT=2
    !NFPOINT=2
    !NTPN=4

    !VERR(1)=4.
    !VERR(2)=7.
    !VERR(3)=-1.
    !VERR(4)=0
    !MAL(1,1)=7;MAL(1,2)=2;MAL(1,3)=1;MAL(1,4)=-2
    !MAL(2,1)=9;MAL(2,2)=15;MAL(2,3)=3;MAL(2,4)=-2
    !MAL(3,1)=-2;MAL(3,2)=-2;MAL(3,3)=11;MAL(3,4)=5
    !MAL(4,1)=1;MAL(4,2)=3;MAL(4,3)=2;MAL(4,4)=13
    !MAL(1,1)=7;MAL(1,2)=2;MAL(1,3)=1;MAL(1,4)=-2
    !MAL(2,1)=0;MAL(2,2)=15;MAL(2,3)=3;MAL(2,4)=-2
    !MAL(3,1)=0;MAL(3,2)=0;MAL(3,3)=11;MAL(3,4)=5
    !MAL(4,1)=0;MAL(4,2)=0;MAL(4,3)=0;MAL(4,4)=13

	!DEALLOCATE(SH,SA,SA1,SX,SZ,SXX)
	WRITE(*,*)'		END MATRIX'

	RETURN
END SUBROUTINE

SUBROUTINE PADCOEF  !空间一阶导系数
    USE GREENMOD
    USE CUMOD
    USE NMLMOD

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

    !IMPLICIT NONE
    !WFA,VOL,XBUO(3)    
    COMMON /GWA/GA(400),GW(400)
	DIMENSION ISM(3,2)
	DATA ISM/1,1,1,1,-1,1/

    REAL:: XI,ETA,AJ
    REAL:: SN(9),SNX(9),SNE(9),SR(3),DN(3)
	REAL:: A(3,3),B(3,3),C(2,2),VV(3)
	REAL:: IS(3),JS(3)
    REAL:: TF1(6),TPF(6)
    REAL:: FCORL(9,3),FCNL(9,3)
	INTEGER I,J,K,L
    INTEGER THR1,THR2,CP
	REAL:: T,D    
    REAL:: TPVX(3,2),TPN,TPE,TPTN,TPTE,TPHE,TPHN
	REAL:: XXE(2,3),TSD(6,6)
    REAL:: DALP,DV
	
    !M=NBPOINT
    !L=NFPOINT

    DNODE=0
    PHIS_X=0
    PHIS_Y=0
    PHIS_Z=0
    PHIS_TX=0
    ETAW_X=0
    ETAW_Y=0

    !GOTO 20

    !START PARALLEL
    CALL OMP_SET_NUM_THREADS(NTHREAD)

!!!$OMP PARALLEL PRIVATE(CP,THR1,THR2,FCORL,FCNL,ISYM,XMN,XI,ETA,AJ,SN,SNX,SNE,SR,DN,TF1,TPF)
    
    !WRITE(*,*)"WF1",WF1

    CP=OMP_GET_THREAD_NUM()+1

    THR1=SUM(NBS(1:CP))-NBS(CP) 
    THR2=SUM(NBS(1:CP))
    
    !WF1(:,CP)=0
    !WRITE(*,*)CP
    THR1=0
    CP=1

    !GOTO 40  !用方法2求导
    !节点编号、单元内节点编号、重单元数目

    ALLOCATE(DPXC(NTPN,9,4),DPYC(NTPN,9,4),DPZC(NTPN,9,4),DPNC(NTPN,9,4),DEXC(NTPN,9,4),DEYC(NTPN,9,4))
    ALLOCATE(APX(NTPN,NTPN),APY(NTPN,NTPN),APZ(NTPN,NTPN),AEX(NTPN,NTPN),AEY(NTPN,NTPN),APEX(NTPN,NTPN),AENX(NTPN,NTPN),APN(NTPN,NTPN))
    ALLOCATE(DPM(NTPN,9,4))

    DPXC=0.
    DPYC=0.
    DPZC=0.
    DPNC=0.

    GOTO 50

    DO IE=1,NTP !NBS(CP)
    !DO IE=1,NTP !NBS(CP)
!		依次将每一面元的四个顶点的坐标赋值到CORL(4,3)中
!		依次将每一面元的四个顶点的法向矢量赋值到CNL(4,3)中
		
		DO I1=1,PAN(IE)
			DO J1=1,3
				FCORL(I1,J1)=CORD(INE(THR1+IE,I1),J1)
				FCNL(I1,J1)=VECN(INE(THR1+IE,I1),J1) 
   	        END DO
            !FCORL(I1,3)=ETAW(INE(THR1+IE,I1))
        END DO
            

        DO IP=1,PAN(IE)

		IF(PAN(IE).EQ.3) THEN
            !GOTO 100
            DNODE(INE(IE,IP))=DNODE(INE(IE,IP))+1.
            IF(IP.EQ.1) THEN
                XI=1.
                ETA=0
            END IF
            IF(IP.EQ.2) THEN
                XI=0
                ETA=1.
            END IF
            IF(IP.EQ.3) THEN
                XI=0.
                ETA=0.
            END IF
            !CALL ISOPAR_3(XI,ETA,FCORL(1:PAN(IE),:),SN(1:PAN(IE)),SNX(1:PAN(IE)),SNE(1:PAN(IE))&
            !                           ,AJ,SR,DN)
        END IF

		IF(PAN(IE).EQ.4) THEN
            DNODE(INE(IE,IP))=DNODE(INE(IE,IP))+1.
            IF(IP.EQ.1) THEN
                XI=-1.
                ETA=-1.
            END IF
            IF(IP.EQ.2) THEN
                XI=1.
                ETA=-1.
            END IF
            IF(IP.EQ.3) THEN
                XI=1.
                ETA=1.
            END IF
            IF(IP.EQ.4) THEN
                XI=-1.
                ETA=1.
            END IF        
            CALL ISOPAR_4(XI,ETA,FCORL(1:PAN(IE),:),SN(1:PAN(IE)),SNX(1:PAN(IE)),SNE(1:PAN(IE))&
                                       ,AJ,SR,DN)
        END IF
        
        IF(PAN(IE).EQ.9) THEN 
            DNODE(INE(IE,IP))=DNODE(INE(IE,IP))+1.
            IF(IP.EQ.1) THEN
                XI=-1.
                ETA=-1.
            END IF
            IF(IP.EQ.2) THEN
                XI=0.
                ETA=-1.
            END IF
            IF(IP.EQ.3) THEN
                XI=1.
                ETA=-1.
            END IF
            IF(IP.EQ.4) THEN
                XI=1.
                ETA=0.
            END IF
            IF(IP.EQ.5) THEN
                XI=1.
                ETA=1.
            END IF
            IF(IP.EQ.6) THEN
                XI=0.
                ETA=1.
            END IF
            IF(IP.EQ.7) THEN
                XI=-1.
                ETA=1.
            END IF
            IF(IP.EQ.8) THEN
                XI=-1.
                ETA=0.
            END IF    
            IF(IP.EQ.9) THEN
                XI=0.
                ETA=0.
            END IF                      
            CALL ISOPAR_9(XI,ETA,FCORL(1:PAN(IE),:),SN(1:PAN(IE)),SNX(1:PAN(IE)),SNE(1:PAN(IE))&
                                       ,AJ,SR,DN)
        END IF     
                  
        DO I=1,3
		    XXE(1,i)=0.0
			XXE(2,i)=0.0
        END DO

		DO I=1,3
	        DO J=1,PAN(IE)
		        xxe(1,i)=xxe(1,i)+Fcorl(j,i)*snx(j)
   				xxe(2,i)=xxe(2,i)+Fcorl(j,i)*sne(j)
            END DO
        END DO
	
    	DO I=1,3
		    DO J=1,2
			    B(I,J)=XXE(J,I)
            END DO
        END DO

    	DO I=1,2
		    DO J=1,2
                C(I,J)=XXE(J,I)
            END DO
        END DO

		DO I=1,3
		    B(I,3)=DN(I)
        END DO

		DO I=1,3
		    DO J=1,3
			    A(I,J)=B(I,J)
            END DO
        END DO
        
        CALL BRINV(A,3)
        IF(IE.GT.NBBLOCK) CALL BRINV(C,2)


		DO I1=1,PAN(IE)
            DPXC(INE(IE,IP),I1,DNODE(INE(IE,IP)))=(A(1,1)*SNX(I1)+A(2,1)*SNE(I1))  !PHIS空间偏导数对应的PHIS前的乘数
            DPYC(INE(IE,IP),I1,DNODE(INE(IE,IP)))=(A(1,2)*SNX(I1)+A(2,2)*SNE(I1)) !节点编号、单元内节点编号、重单元数目
            DPZC(INE(IE,IP),I1,DNODE(INE(IE,IP)))=(A(1,3)*SNX(I1)+A(2,3)*SNE(I1))  
            
            DEXC(INE(IE,IP),I1,DNODE(INE(IE,IP)))=(C(1,1)*SNX(I1)+C(2,1)*SNE(I1))  !PHIS空间偏导数对应的PHIS前的乘数
            DEYC(INE(IE,IP),I1,DNODE(INE(IE,IP)))=(C(1,2)*SNX(I1)+C(2,2)*SNE(I1))
            DPNC(INE(IE,IP),I1,DNODE(INE(IE,IP)))=(A(3,1))*SN(I1)    !PHIS_N前的乘数
            DPM(INE(IE,IP),I1,DNODE(INE(IE,IP)))=INE(IE,I1) !系数对应的编号，节点编号、单元内节点编号、重单元数目 
        END DO
        END DO
    END DO        


    DO IP=NBPOINT+1,NTPN  
        !DO I1=1,PAN(IE)    
           !DPC(IP,I1,1)=DPC(IP,I1,1)/DNODE(I) 
           !DPC(IP,I1,2)=DPC(IP,I1,2)/DNODE(I)
        !END DO
    END DO

50  CALL TRANMATRIX_PX2P
END SUBROUTINE

SUBROUTINE TRANMATRIX_PX2P
    USE GREENMOD
    USE CUMOD
    USE NMLMOD

    GOTO 50

    APX=0.
    APY=0.
    AEX=0.
    AEY=0.
    APZ=0.
    APN=0.
    DO IP=1,NTPN
        DO J=1,INT(DNODE(IP))
            DO K=1,9 !PAN(K,) DNODE(IP) 面元内节点数目
                APX(IP,DPM(IP,K,J))=APX(IP,DPM(IP,K,J))+DPXC(IP,K,J)/DNODE(IP)
                APY(IP,DPM(IP,K,J))=APY(IP,DPM(IP,K,J))+DPYC(IP,K,J)/DNODE(IP)
                APZ(IP,DPM(IP,K,J))=APZ(IP,DPM(IP,K,J))+DPZC(IP,K,J)/DNODE(IP)
                AEX(IP,DPM(IP,K,J))=AEX(IP,DPM(IP,K,J))+DEXC(IP,K,J)/DNODE(IP)
                AEY(IP,DPM(IP,K,J))=AEY(IP,DPM(IP,K,J))+DEYC(IP,K,J)/DNODE(IP)
                APN(IP,DPM(IP,K,J))=APN(IP,DPM(IP,K,J))+DPNC(IP,K,J)/DNODE(IP)

                !APX(IP,DPM(IP,K,J))=DPXC(IP,K,J)
                !APY(IP,DPM(IP,K,J))=DPYC(IP,K,J)
                !APZ(IP,DPM(IP,K,J))=DPZC(IP,K,J)
                !AEX(IP,DPM(IP,K,J))=DEXC(IP,K,J)
                !AEY(IP,DPM(IP,K,J))=DEYC(IP,K,J)
                !APN(IP,DPM(IP,K,J))=DPNC(IP,K,J)

                !ETAW_X=AEX*ETAW
                !PHIS_X=APX*PHIS+APN*PHIS
            END DO
        END DO
    END DO
        
    APEX=0.
    AENX=0.
    DO I=1,NTPN
        DO J=1,NTPN
            DO K=1,NTPN
                APEX(I,J)=APEX(I,J)+AEX(I,K)*APX(K,J)
                AENX(I,J)=AENX(I,J)+AEX(I,K)*APN(K,J)
            END DO
        END DO
        AENX(I,I)=AENX(I,I)+1.
    END DO    
          
    !CALL BRINV(AENX(NBPOINT+1:NTPN,NBPOINT+1:NTPN),NFPOINT)      
            
50  DEALLOCATE(DPXC,DPYC,DPZC,DEXC,DEYC,DPM,DPNC)       
END SUBROUTINE

SUBROUTINE EQUATION

	USE EQUATIONMOD
	USE GREENMOD

	DOUBLE PRECISION T
	CHARACTER*8 CHAR_TIME3,CHAR_TIME4
	CHARACTER*9 CHAR_DATE3,CHAR_DATE4

	OPEN(1,FILE='DB_FLOW.DAT')
	NNN=NBPOINT !NTPN !+NFPOINT
	WRITE(*,*) 
	write(*,102)NNN,NNN
102	format(1x,'SUB equation:'i5,'*',i5)
	!CALL TIME(char_time3)
	!CALL DATE(CHAR_DATE3)
	PRINT *, '    BEGIN TIME: ', CHAR_TIME3,' ON  ',CHAR_DATE3
	
	ALLOCATE(X(NNN),JS(NNN))

!c	CALL AGAUS(A,B,NNN,x,L,JS)
	N=NNN

	L=1
	DO 250 K=1,N-1
	  D=0.0
	  DO 210 I=K,N
	  DO 210 J=K,N
	    IF (ABS(MAL(I,J)).GT.D) THEN
	      D=ABS(MAL(I,J))
	      JS(K)=J
	      IS=I
	    END IF
210	  CONTINUE
	  IF (D+1.0.EQ.1.0) THEN
	    L=0
	  ELSE
	    IF (JS(K).NE.K) THEN
	      DO 220 I=1,N
	        T=MAL(I,K)
	        MAL(I,K)=MAL(I,JS(K))
	        MAL(I,JS(K))=T
220	      CONTINUE
	    END IF
	    IF (IS.NE.K) THEN
	      DO 230 J=K,N
	        T=MAL(K,J)
	        MAL(K,J)=MAL(IS,J)
	        MAL(IS,J)=T
230	      CONTINUE
	      T=VERR(K)
	      VERR(K)=VERR(IS)
	      VERR(IS)=T
	    END IF
	  END IF
	  IF (L.EQ.0) THEN
	    WRITE(*,299)
!C	    RETURN
		GOTO 31
	  END IF
	  DO 240 J=K+1,N
	    MAL(K,J)=MAL(K,J)/MAL(K,K)
240	  CONTINUE
	  VERR(K)=VERR(K)/MAL(K,K)
	  DO 260 I=K+1,N
	    DO 270 J=K+1,N
	      MAL(I,J)=MAL(I,J)-MAL(I,K)*MAL(K,J)
270	    CONTINUE
	    VERR(I)=VERR(I)-MAL(I,K)*VERR(K)
260	  CONTINUE
250	CONTINUE
	IF (ABS(MAL(N,N))+1.0.EQ.1.0) THEN
	  L=0
	  WRITE(*,299)
!C	  RETURN
	  GOTO 31
	END IF
	X(N)=VERR(N)/MAL(N,N)
	DO 280 I=N-1,1,-1
	  T=0.0
	  DO 290 J=I+1,N
	    T=T+MAL(I,J)*X(J)
290	  CONTINUE
	  X(I)=VERR(I)-T
280	CONTINUE
299	FORMAT(1X,' FAIL ')
	JS(N)=N
	DO 291 K=N,1,-1
	  IF (JS(K).NE.K) THEN
	    T=X(K)
	    X(K)=X(JS(K))
	    X(JS(K))=T
	  END IF
291	CONTINUE

    IF(L.NE.0)THEN
		WRITE(*,*)'    SUCCESSFULLY SOLVE THE INTEGRATION EQUATION'
		DO 300 I=1,NNN
		WRITE(1,*)X(I)
!CWRITE(*,*)X(I)
300		CONTINUE
	ELSE 
		WRITE(*,*)'    THERE IS SOME PROBLEM IN SOLVING THE EQUATION'
	END IF
	!CALL TIME(char_time4)
	!CALL DATE(CHAR_DATE4)
	PRINT *, '    END TIME: ', CHAR_TIME4,'  ON  ',CHAR_DATE4
31	CLOSE(1)	
!C 释放内存
	DEALLOCATE(X,JS,MAL,VERR)

	RETURN
END

SUBROUTINE COND_MAL_NPRE
    USE GREENMOD

    REAL,ALLOCATABLE:: TPMAL(:,:),SUMMAL(:)
    REAL:: NORM1A,NORM1A_

    ALLOCATE(SUMMAL(NTPN),TPMAL(NTPN,NTPN))

    SUMMAL=0.
    DO J=1,NTPN
        DO I=1,NTPN
            SUMMAL(J)=SUMMAL(J)+ABS(MAL(J,I))
        END DO
    END DO

    NORM1A=0.
    DO I=1,NTPN
        IF(NORM1A.LE.SUMMAL(I)) NORM1A=SUMMAL(I)
    END DO
    DO J=1,NTPN
        DO I=1,NTPN
            !NORM1A=NORM1A+MAL(I,J)**2
        END DO
    END DO
    !NORM1A=SQRT(NORM1A)

    !WRITE(*,*)"NORM1A MATRIX A",NORM1A

    DO I=1,NTPN
        DO J=1,NTPN
            TPMAL(I,J)=MAL(I,J)
        END DO
    END DO
    CALL BRINV(TPMAL,NTPN)
    SUMMAL=0.
    DO J=1,NTPN
        DO I=1,NTPN
            SUMMAL(J)=SUMMAL(J)+ABS(TPMAL(J,I))
        END DO
    END DO

    NORM1A_=0.
    DO I=1,NTPN
        IF(NORM1A_.LE.SUMMAL(I)) NORM1A_=SUMMAL(I)
    END DO
    DO J=1,NTPN
        DO I=1,NTPN
            !NORM1A_=NORM1A_+TPMAL(I,J)**2
        END DO
    END DO
    !NORM1A_=SQRT(NORM1A_)

    WRITE(*,*)NORM1A,NORM1A_,NORM1A*NORM1A_
    WRITE(111,121)FR,NORM1A,NORM1A_,NORM1A*NORM1A_
    121FORMAT("ORI_MAL   ",5F15.3)

    DEALLOCATE(SUMMAL,TPMAL)

END SUBROUTINE


SUBROUTINE COND_MAL
    USE GREENMOD

    REAL,ALLOCATABLE:: TPMAL(:,:),SUMMAL(:),TPMAM(:,:),TPMAM_(:,:)
    REAL:: NORM1A,NORM1A_

    ALLOCATE(SUMMAL(NTPN),TPMAL(NTPN,NTPN),TPMAM(NTPN,NTPN),TPMAM_(NTPN,NTPN))

    TPMAM=0.
    DO I=1,NTPN
        DO J=1,NTPN
            DO K=1,NTPN
                TPMAM(I,J)=TPMAM(I,J)+LAG(I,K)*UA(K,J)
            END DO
        END DO
    END DO
    TPMAM_=TPMAM
    CALL BRINV(TPMAM_,NTPN)

    TPMAL=0.
    DO I=1,NTPN
        DO J=1,NTPN
            DO K=1,NTPN
                TPMAL(I,J)=TPMAL(I,J)+TPMAM_(I,K)*MAL(K,J)
            END DO
        END DO
    END DO

    OPEN(11,FILE='TEST_DIAG.DAT')
    SUMMAL=0.
    DO J=1,NTPN
        DO I=1,NTPN
            SUMMAL(J)=SUMMAL(J)+ABS(TPMAL(J,I))
        END DO
        WRITE(11,*)J,TPMAL(J,J)
    END DO

    NORM1A=0.
    DO I=1,NTPN
        IF(NORM1A.LE.SUMMAL(I)) NORM1A=SUMMAL(I)
    END DO
    DO J=1,NTPN
        DO I=1,NTPN
            !NORM1A=NORM1A+TPMAL(I,J)**2
        END DO
    END DO
    !NORM1A=SQRT(NORM1A)

    !WRITE(*,*)"NORM1A MATRIX A",NORM1A

!***************************************************************

    DO I=1,NTPN
        DO J=1,NTPN
            !TPMAL(I,J)=MAL(I,J)
        END DO
    END DO
    CALL BRINV(TPMAL,NTPN)

    TPMAM_=0.
    DO I=1,NTPN
        DO J=1,NTPN
            DO K=1,NTPN
                !TPMAM_(I,J)=TPMAM_(I,J)+TPMAM(I,K)*TPMAL(K,J)
            END DO
        END DO
    END DO

    SUMMAL=0.
    DO J=1,NTPN
        DO I=1,NTPN
            !SUMMAL(J)=SUMMAL(J)+ABS(TPMAM_(J,I))
            SUMMAL(J)=SUMMAL(J)+ABS(TPMAL(J,I))
        END DO
        WRITE(11,*)J,TPMAL(J,J)
    END DO
    CLOSE(11)

    NORM1A_=0.
    DO I=1,NTPN
        IF(NORM1A_.LE.SUMMAL(I)) NORM1A_=SUMMAL(I)
    END DO
    DO J=1,NTPN
        DO I=1,NTPN
            !NORM1A_=NORM1A_+TPMAM_(I,J)**2
        END DO
    END DO
    !NORM1A_=SQRT(NORM1A_)

    WRITE(*,*)NORM1A,NORM1A_,NORM1A*NORM1A_
    WRITE(111,121)NORM1A,NORM1A_,NORM1A*NORM1A_
    WRITE(*,*)
    121FORMAT('PRE_CON    ',5F15.3)

    DEALLOCATE(SUMMAL,TPMAL,TPMAM,TPMAM_)

END SUBROUTINE