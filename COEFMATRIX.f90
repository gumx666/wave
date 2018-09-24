SUBROUTINE COEFMATRIX_SP
!   ***************************************************
!   GENERATE THE GLOBAL MATRIX OF 1st INTEGRAL EQUATION  
!   ***************************************************
    


    USE GREENMOD

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

    !IMPLICIT NONE
    CHARACTER MARK
	INTEGER:: ISM(3,4)
	DATA ISM/1,1,1,1,-1,1,1,-1,-1,1,1,-1/
	INTEGER:: NLOOP,K,K1,K2,CP,NBLOOP,NPLOOP,IP,I3,I1,J1,IE,J,ISYM,ISINGULAR,I,I2
    REAL ::FCORL(9,3),FCNL(9,3),FPS(3),T,TL,FAIJ(4,6),BETA,LOFF
    REAL ::HIJX(9),HIJN(9),HIJZ(9),HYX(9),HYN(9),HYZ(9),FHIJ(9) !诱导系数 

    WRITE(*,*) 
    WRITE(*,*)'SUB COEFMATRIX..............'

    open(18,file='coefuse_SP.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(18)
    open(19,file='coefuse_T.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(19)
    
    open(118,file='coefuse.DAT')

    NRPOINT=0.

	ALLOCATE (PS(3),CORL(9,3),CNL(9,3))
	ALLOCATE (SSH1(NTPN+NRPOINT),SSH(NTPN+NRPOINT),SSA(NTPN+NRPOINT),SSA1(NTPN+NRPOINT),HIJ(9),AIJ(9),AIJ1(9),HIJ1(9))
	ALLOCATE (SUMHIJ(NTPN+NRPOINT))
    WRITE(*,*)"NTPN+NRPOINT",NTPN,NRPOINT
    ALLOCATE (RSH(NTPN+NRPOINT,NTPN+NRPOINT,1),RSA(NTPN+NRPOINT,NTPN+NRPOINT),RFSA(NBPOINT,6))

    RSH=0.;RSA=0.;RFSA=0.

    m=nbpoint
    L=NFPOINT
    NRPOINT=0.
    NRBLOCK=0.
!  NTPN+NRPOINT IN SINGLE CPU
   IF(MOD(M+L+NRPOINT,NTHREAD).NE.0) THEN
     DO I=1,NTHREAD-1
       IF(MOD(M+L+NRPOINT,NTHREAD).EQ.I) THEN
         DO J=1,NTHREAD-1
            NPS1(J)=(M+L+NRPOINT-I)/NTHREAD
         END DO
         NPS1(NTHREAD)=NPS1(NTHREAD-1)+I
         EXIT
       END IF
     END DO
   ELSE 
     DO I=1,NTHREAD
       NPS1(I)=(M+L+NRPOINT)/NTHREAD
     END DO
   END IF

    DO IP=1,NTPN
		SUMHIJ(IP)=0.0
    END DO

    IF(STAGE.EQ.2) THEN
        NPLOOP=NTPN
        NBLOOP=NTP
    END IF

    IF(STAGE.EQ.1) THEN
        NPLOOP=NBPOINT
        NBLOOP=NBBLOCK
    END IF

    K1=0
    CP=1

!!START PARALLEL
!    CALL OMP_SET_NUM_THREADS(NTHREAD)

!!$omp parallel private(CP,k1,k2,FCORL,FCNL,ISYM,FPS,T,HIJ,AIJ,IP,I3,IE,I1,I2,J1,ISINGULAR,TL,K,FAIJ,FHIJ)
    
!    CP=OMP_GET_THREAD_NUM()+1

!    K1=SUM(NPS1(1:CP))-NPS1(CP) 
!    K2=SUM(NPS1(1:CP))
    !K1=0
    !WRITE(*,*)"K1   ",K1

    !ALLOCATE (HIJ(9),AIJ(9))

    DO IP=1,NPLOOP
    !DO IP=1,NPS1(CP)

        !LOFF=LDC*DMIP(IP)**BETA
		IF((((IP+K1)/100)*100).EQ.IP+K1)THEN
			WRITE(*,200)IP+K1,NPLOOP !NTPN
		END IF
  200	FORMAT('16H    MOVING BALL',3X,'6HINODE=',I8,5X,I8)
   
!		INITIAL THE DIMENSIONS
!		SSH(K)中存放有奇点部分的偶极子效应
!		SSA(K)中存放有奇点部分的源效应
        DO 5 I3=1,NPLOOP
			SSH(I3)=0.
            SSH1(I3)=0.
			SSA(I3)=0.
			SSA1(I3)=0.0
    5		CONTINUE						  
!		将每一个点的坐标赋值到ps（3）中
        
        BETA=0.5
        !LOFF=DMIP(IE)**BETA  !单元内潜上置距离  
        DO 10 J1=1,3
   10		FPS(J1)=CORD(IP+K1,J1) !+VECN(IP,J1)*LOFF
!		LOOP FOR ELEMENTS
	    !FPS(3)=FPS(3)+0.05
    	DO IE=1,NBLOOP
!			WRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!			依次将每一面元的PAN(IE)个顶点的坐标赋值到corl(PAN(IE),3)中
!			依次将每一面元的PAN(IE)个顶点的法向矢量赋值到cnl(PAN(IE),3)中
			DO ISYM=1,4
				
			DO I1=1,PAN(IE)
				I2=I1
                IF(PAN(IE).EQ.3) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=5-I1	
                END IF                
                IF(PAN(IE).EQ.4) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				    IF(MOD(ISYM,4).EQ.0.AND.(I1.NE.1)) I2=6-I1	
                END IF

                IF(PAN(IE).EQ.8) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                END IF

                IF(PAN(IE).EQ.9) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                END IF

                DO J1=1,3
                    FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)
                    !STAGE=1 MIRROR OF FREE SURFACE    STAGE=2 MIRROR OF BOTTOM
                    IF(STAGE.EQ.2) THEN
                        IF(ISYM.EQ.3.OR.ISYM.EQ.4) THEN
                            FCORL(I1,3)=FCORL(I1,3)-2.*HZ
                        END IF
				    END IF
                    FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM)  
                END DO
            END DO


!			判断是否是奇点以调用不同的程序，判断标志为ISINGULAR
            !WRITE(*,*)"ps",ps
            !WRITE(*,*)"fcorl",fcorl(1:4,:)
			ISINGULAR=0
			DO 50 I1=1,PAN(IE)
				T=0.
				DO 55 J1=1,3
   55				T=T+(FCORL(I1,J1)-FPS(J1))**2
				    T=SQRT(T)
				    IF(T.LT.1.E-7) THEN
					    ISINGULAR=1
				    END IF
   50		CONTINUE

!c			是奇点则调用cofahs，否则调用cofah,aij为源效应，hij为偶极子效应 
			DO I=1,PAN(IE)
				FHIJ(I)=0.0
				AIJ(I)=0.0
            END DO
              					
			IF(ISINGULAR.EQ.1) THEN
                !IF(PAN(IE).EQ.4) CALL COFAHS_4(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),HIJ(1:PAN(IE)))
                IF(PAN(IE).EQ.9) CALL COFAHS_SP9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE)))
                IF(PAN(IE).EQ.8) CALL COFAHS_SP8(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE)))
                !IF(IP.EQ.1.AND.IE.LE.1) WRITE(*,*)IP,IE,HIJ(1:PAN(IE))
			ELSE
                !IF(PAN(IE).EQ.4) CALL COFAH_4(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),HIJ(1:PAN(IE)))
                IF(PAN(IE).EQ.9) CALL COFAH_SP9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE)))
			    IF(PAN(IE).EQ.8) CALL COFAH_SP8(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE)))
			ENDIF

            IF(IP+K1.EQ.1) THEN
                !CALL COFAH_4SN(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),XGRA,FAIJ)   !与PHIS_T的边界积分方程有关
            END IF

! 			将每一面元中各顶点处的积分值找到对应的总序列中的位置，并进行叠加
! 			SSH(k)中存放偶极子效应
! 			SSA(k)中存放界面上的源效应
!			SUMHIJ(NTPN)用于存放对角线上的量

			DO I1=1,PAN(IE)
                !WRITE(*,*)IE
				K=INE(IE,I1)
				I2=I1
                IF(PAN(IE).EQ.3) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=5-I1	
                END IF                
                IF(PAN(IE).EQ.4) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				    IF(MOD(ISYM,4).EQ.0.AND.(I1.NE.1)) I2=6-I1	
                END IF
                IF(PAN(IE).EQ.8) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                END IF
                IF(PAN(IE).EQ.9) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                END IF

                RSA(IP+K1,K)=RSA(IP+K1,K)+AIJ(I2)

                TL=SQRT((CORD(IP+K1,1)-CORD(K,1))**2+(CORD(IP+K1,2)-CORD(K,2))**2+(CORD(IP+K1,2)-CORD(K,2))**2)
                IF(K.NE.IP+K1)THEN
                    !IF(TL.GT.1.E-6) THEN
                        !IF(IP+K1.LE.NTPN) THEN
				            !SSH(K)=SSH(K)+HIJ(I2)
                            RSH(IP+K1,K,1)=RSH(IP+K1,K,1)+FHIJ(I2)   !D(1/R)/D(X)
                            !END IF
				    !END IF
                END IF
				    
                !IF(K.EQ.IP+K1)THEN
                IF(K.NE.IP+K1)THEN
                !IF(TL.GT.1.E-6) THEN
                !IF(TL.LE.1.E-6) THEN
					SUMHIJ(IP+K1)=SUMHIJ(IP+K1)+FHIJ(I2)   !D(1/R)/D(N)  
                END IF

                !END IF
                !IF(IP.EQ.1) WRITE(*,*)SSA(K),SSH(K)
                IF(ISYM.LE.2) THEN
                    !RSH(IP+K1,K)=SSH(K)
                    !RSA(IP+K1,K)=SSA(K)
                END IF  
                
              
                    
            END DO
        END DO
        END DO
    END DO
!!$OMP END PARALLEL

    IF(STAGE.EQ.2) THEN
        NPLOOP=NBPOINT+NFPOINT
    END IF

    DO I=1,NPLOOP
        DO J=1,NPLOOP
        	WRITE(18)RSA(I,J),RSH(I,J,1)  
        END DO
        111FORMAT(10F15.6)
        !WRITE(*,*)RSH(I,3,1),SUMHIJ(I)
    END DO    


	DEALLOCATE(PS,CORL,CNL,SSH1,SSH,SSA,SSA1,HIJ,AIJ,RSH,RSA,RFSA,HIJ1,AIJ1)
	WRITE(*,*)'		END COEFMATRIX'

	CLOSE(18)
    CLOSE(118)
    RETURN
END SUBROUTINE

SUBROUTINE COEFMATRIX                                                                                                                                   
!   ***************************************************
!   GENERATE THE GLOBAL MATRIX OF 1st INTEGRAL EQUATION  
!   ***************************************************
    
    USE GREENMOD
    USE INPUTDATA
    
#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

    !IMPLICIT NONE
    CHARACTER MARK
	INTEGER:: ISM(3,4)
	DATA ISM/1,1,1,1,-1,1,1,-1,-1,1,1,-1/
	INTEGER:: NLOOP,K,K1,K2,CP,NBLOOP,NPLOOP,IP,I3,I1,J1,IE,J,ISYM,ISINGULAR,I,I2,I11,I22,J11,IAF(2),CAF,IBO(4),IUP1(3),IUP2(3)
    REAL ::FCORL(9,3),FCNL(9,3),FPS(3),T,TL,FAIJ(4,6),BETA,LOFF,CORL1(9,3),CORL2(9,3),SR(3),D1,D2,D3,D4,XI(6),ETA(6)
    REAL ::HIJX(9),HIJN(9),HIJZ(9),HYX(9),HYN(9),HYZ(9),FHIJ(9,4),VV(3) !诱导系数 
    REAL ::SN(9),SNX(9),SNE(9),AJ,DN(3)
    REAL,ALLOCATABLE:: CORD11(:,:)
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2

    ALLOCATE(CORD11(NTPN,3))
	ALLOCATE(SH(NTPN,NTPN),SX(NTPN,NTPN),SZ(NTPN,NTPN),SA(NTPN,NTPN),SA1(NTPN,NTPN),SY(NTPN,NTPN))

    CORD11=CORD

    WRITE(*,*)"DFY=     ",DFY,0.6*DFX
    WRITE(*,*) 
    WRITE(*,*)'SUB COEFMATRIX..............'

    open(18,file='coefuse.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(18)
    open(19,file='coefuse_T.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(19)
    
    open(118,file='coefuse.DAT')

    NRPOINT=0.

    ALLOCATE (IDS(NTPN))
	ALLOCATE (PS(3),CORL(9,3),CNL(9,3))
	ALLOCATE (SSH1(NTPN+NRPOINT),SSH(NTPN+NRPOINT),SSA(NTPN+NRPOINT),SSA1(NTPN+NRPOINT),HIJ(9),AIJ(9),AIJ1(9),HIJ1(9))
	ALLOCATE (SUMHIJ(NTPN+NRPOINT))
    WRITE(*,*)"NTPN+NRPOINT",NTPN,NRPOINT
    ALLOCATE (RSH(NTPN+NRPOINT,NTPN+NRPOINT,4),RSA(NTPN+NRPOINT,NTPN+NRPOINT),RFSA(NBPOINT,6))

    OPEN(4,FILE='INDUCED_VELO.DAT')

    RSH=0.;RSA=0.;RFSA=0.

    m=nbpoint
    L=NFPOINT
    NRPOINT=0.
    NRBLOCK=0.
!  NTPN+NRPOINT IN SINGLE CPU
   IF(MOD(M+L+NRPOINT,NTHREAD).NE.0) THEN
     DO I=1,NTHREAD-1
       IF(MOD(M+L+NRPOINT,NTHREAD).EQ.I) THEN
         DO J=1,NTHREAD-1
            NPS1(J)=(M+L+NRPOINT-I)/NTHREAD
         END DO
         NPS1(NTHREAD)=NPS1(NTHREAD-1)+I
         EXIT
       END IF
     END DO
   ELSE 
     DO I=1,NTHREAD
       NPS1(I)=(M+L+NRPOINT)/NTHREAD
     END DO
   END IF

    DO IP=1,NTPN
		SUMHIJ(IP)=0.0
    END DO

    IF(STAGE.EQ.2) THEN
        NPLOOP=NTPN
        NBLOOP=NTP
    END IF

    IF(STAGE.EQ.1) THEN
        NPLOOP=NBPOINT
        NBLOOP=NBBLOCK
    END IF

    K1=0
    CP=1

!!START PARALLEL
!    CALL OMP_SET_NUM_THREADS(NTHREAD)

!!$omp parallel private(CP,k1,k2,FCORL,FCNL,ISYM,FPS,T,HIJ,AIJ,IP,I3,IE,I1,I2,J1,ISINGULAR,TL,K,FAIJ,FHIJ)
    
!    CP=OMP_GET_THREAD_NUM()+1

!    K1=SUM(NPS1(1:CP))-NPS1(CP) 
!    K2=SUM(NPS1(1:CP))
    !K1=0
    !WRITE(*,*)"K1   ",K1

    !ALLOCATE (HIJ(9),AIJ(9))

    IDS=0
    DO IP=1,NPLOOP
    !DO IP=1,NPS1(CP)

        !LOFF=LDC*DMIP(IP)**BETA
		IF((((IP+K1)/500)*500).EQ.IP+K1)THEN
			WRITE(*,200)IP+K1,NPLOOP !NTPN
		END IF
  200	FORMAT('16H    MOVING BALL',3X,'6HINODE=',I8,5X,I8)
   
!		INITIAL THE DIMENSIONS
!		SSH(K)中存放有奇点部分的偶极子效应
!		SSA(K)中存放有奇点部分的源效应
        DO 5 I3=1,NPLOOP
			SSH(I3)=0.
            SSH1(I3)=0.
			SSA(I3)=0.
			SSA1(I3)=0.0
    5		CONTINUE						  
!		将每一个点的坐标赋值到ps（3）中
        
        BETA=0.5
        !LOFF=DMIP(IE)**BETA  !单元内潜上置距离  
        DO 10 J1=1,3
   10		FPS(J1)=CORD(IP+K1,J1) !+VECN(IP,J1)*LOFF
!		LOOP FOR ELEMENTS
	    !FPS(3)=FPS(3)+0.05
    	DO IE=1,NBLOOP
!			WRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!			依次将每一面元的PAN(IE)个顶点的坐标赋值到corl(PAN(IE),3)中
!			依次将每一面元的PAN(IE)个顶点的法向矢量赋值到cnl(PAN(IE),3)中
			
            !IF(MNKDB.EQ.2) THEN
            !    IF(IE.LE.NBBLOCK) THEN
             !       KK=4
            !    ELSE
            !        KK=2
            !    END IF
            !ELSE IF(MNKDB.EQ.1) THEN
            !    KK=2
            !END IF

            
            DO ISYM=1,KK
				
			DO I1=1,PAN(IE)
				I2=I1
                IF(PAN(IE).EQ.3) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=5-I1	
                END IF                
                IF(PAN(IE).EQ.4) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				    IF(MOD(ISYM,4).EQ.0.AND.(I1.NE.1)) I2=6-I1	
                END IF
                IF(PAN(IE).EQ.8) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                END IF
                IF(PAN(IE).EQ.9) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                END IF

                DO J1=1,3
                    IF(IE.GT.NBBLOCK) THEN
                        !LOFF=1.0*DMIP(INE(IE,I1))**BETA  !单元内潜上置距离  
                        !LOFF=1.0*DFX !LINEAR OK
                        !LOFF=0.8*DFX
                        !LOFF=1.2*DFX !S60 OK
                        !LOFF=1.3*DFX
                        !LOFF=1.5*DFX
                        !LOFF=0.2*DFX
                        !LOFF=0.4*DFX
                        LOFF=0.8*DFX
                        LOFF=1.0*DFX
                        !LOFF=0.5*DFX
                        !LOFF=1.5*DFX
                        
                        !LOFF=1.0*DFX
                        !LOFF=0.9*DFX
                        !LOFF=1.2*DFX   !yuzheng
                        !LOFF=0.5*DFX
                    ELSE
                        LOFF=0
                        !LOFF=DFX !S60 OK
                        !LOFF=DFY
                        !LOFF=1.3*DFY
                    END IF
                    !WRITE(*,*)DFX !DMIE(IE),BETA, LOFF
                    
                    IF(IE.GT.NBBLOCK) THEN
                    !IF(ABS(CORD(INE(IE,I2),2)).GE.1.E-6) THEN
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)+VECN(INE(IE,I2),J1)*LOFF*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM) 
                        
                        IF(ISYM.GT.2) THEN
                            FCORL(I1,3)=CORD(INE(IE,I2),3)*ISM(3,ISYM)+VECN(INE(IE,I2),3)*LOFF*ISM(3,ISYM)-2.*NSWIN(13,1)
                        END IF
                    !ELSE
                    !    IF(J1.NE.2) THEN
                    !        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)+&
                    !                    (VECN(INE(IE,I1),J1)+VECN(INE(IE,I1),J1))*LOFF*ISM(J1,ISYM)/2.
                    !    ELSE
                    !        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)+&
                    !                    (VECN(INE(IE,I1),J1)-VECN(INE(IE,I1),J1))*LOFF*ISM(J1,ISYM)/2.                 
                    !    END IF

                    
                    !*********************************
                    goto 600
                    !ELSE IF(IP.GT.NBPOINT.AND.IE.LE.NBBLOCK.AND.&
                    !         ABS(CORD(INE(IE,I2),3)).LE.1.E-5) THEN
                    
                    !IDS(IE)=1

                    !ELSE IF(IP.GT.NBPOINT.AND.IE.LE.NBBLOCK.AND.&
                    !         SQRT((CORD(INE(IE,I2),1)-CORD(IP,1))**2+&
                    !              (CORD(INE(IE,I2),2)-CORD(IP,2))**2+&
                    !              (CORD(INE(IE,I2),3)-CORD(IP,3))**2).LE.1.E-5) THEN

                    !ELSE IF(IE.LE.NBBLOCK.AND.&
                    !         ABS(CORD(INE(IE,I2),3)).LE.1.E-5) THEN

                    !ELSE IF(IP.GT.NBPOINT.AND.IE.LE.NBBLOCK.AND.&
                    !         ABS(CORD(INE(IE,I2),3)).LE.1.E-5.AND.&
                    !         ABS(CORD(INE(IE,I2),2)).LE.1.E-4) THEN

                    !        DO J11=1,3
                    !            FCORL(I2,J11)=CORD(INE(IE,I2),J11)*ISM(J11,ISYM)
                    !            FCORL(I2,J11)=CORD(INE(IE,I2),J11)*ISM(J11,ISYM)+&
                    !                          VECN(INE(IE,I2),J11)*LOFF*ISM(J11,ISYM)*0.3
                    !            FCNL(I2,J11)=VECN(INE(IE,I2),J11)*ISM(J11,ISYM)
                    !        END DO

                    !GOTO 600
                    !*********************************
                        
                        
                        !LOFF1=ABS(CORD(INE(IE,1),3)-CORD(INE(IE,2),3))
                        !DO I11=1,8
                        !    LOFF=ABS(CORD(INE(IE,I11),3)-CORD(INE(IE,I11+1),3))
                        !    IF(LOFF.GE.LOFF1) THEN
                        !        LOFF1=LOFF
                        !    END IF
                        !END DO

                        !GOTO 600

                        DO I11=1,9
                            I22=I11
                            IF(ISYM.EQ.2.AND.I11.NE.1.AND.(I11.NE.9)) I22=10-I11   !保证镜像后的面元节点按逆时针排列
                            IF(ISYM.EQ.4.AND.I11.NE.1.AND.(I11.NE.9)) I22=10-I11   !保证镜像后的面元节点按逆时针排列

                            !IDS(INE(IE,I22))=1

                            VV(1)=-CORD(INE(IE,I22),1)/SQRT(CORD(INE(IE,I22),1)**2+CORD(INE(IE,I22),2)**2)
                            VV(2)=-CORD(INE(IE,I22),2)/SQRT(CORD(INE(IE,I22),1)**2+CORD(INE(IE,I22),2)**2)
                            VV(3)=0.

                            !按面元法向内收
                            DO J11=1,3
                                !FCORL(I11,J1)=CORD(INE(IE,I22),J1)*ISM(J1,ISYM)
                                !FCORL(I11,3)=CORD(INE(IE,I22),3)*ISM(3,ISYM)+ISM(3,ISYM)*2.*dfx !LOFF
                                !write(*,*)CORD(INE(IE,I2),1:3)
                                FCORL(I11,J11)=CORD(INE(IE,I22),J11)*ISM(J11,ISYM)
                                !FCORL(I11,J11)=CORD(INE(IE,I22),J11)*ISM(J11,ISYM)+&
                                !               VECN(INE(IE,I22),J11)*LOFF*ISM(J11,ISYM)*0.4
                                !FCORL(I11,J11)=CORD(INE(IE,I22),J11)*ISM(J11,ISYM)*0.9
                                
                                !FCORL(I11,J11)=CORD(INE(IE,I22),J11)*ISM(J11,ISYM)-&
                                !               VECN(INE(IE,I22),J11)*LOFF*ISM(J11,ISYM)*0.5
                                FCNL(I11,J11)=VECN(INE(IE,I22),J11)*ISM(J11,ISYM)
                            END DO

                            FCORL(I11,2)=CORD(INE(IE,I22),2)*ISM(2,ISYM)+VECN(INE(IE,I22),2)/ABS(VECN(INE(IE,I22),2))&
                                         *LOFF*ISM(2,ISYM)*0.4

                            !IF(ABS(CORD(INE(IE,I22),3)).LE.1.E-5) THEN
                            !    IDS(INE(IE,I22))=1
                            !    FCORL(I11,2)=CORD(INE(IE,I22),2)*ISM(2,ISYM)+&
                            !                   VECN(INE(IE,I22),2)*LOFF*ISM(2,ISYM)*0.5
                            !END IF
                        END DO

                        GOTO 11

                        !计算中点坐标
                        DO J11=1,3
                            FCORL(2,J11)=(FCORL(1,J11)+FCORL(3,J11))/2.
                            FCORL(4,J11)=(FCORL(3,J11)+FCORL(5,J11))/2.
                            FCORL(6,J11)=(FCORL(5,J11)+FCORL(7,J11))/2.
                            FCORL(8,J11)=(FCORL(7,J11)+FCORL(1,J11))/2.
                            FCORL(9,J11)=(FCORL(1,J11)+FCORL(3,J11)+FCORL(5,J11)+FCORL(7,J11))/4.
                        END DO

                        IF(ISYM.EQ.1) THEN
                            DO I11=1,9
                                CORD11(INE(IE,I11),:)=FCORL(I11,:)
                                !WRITE(*,*)CORD11(INE(IE,I11),:)
                            END DO
                            !WRITE(*,*)"OKKKKKKKKKKKKKKKKKKKKKKKK"
                        END IF

                        GOTO 11

                        IBO(1)=1
                        IBO(2)=3
                        IBO(3)=5
                        IBO(4)=7
                        DO I11=3,1,-1
                            DO I22=1,I11
                                IF(FCORL(IBO(I22),3).GT.FCORL(IBO(I22+1),3)) THEN
                                    K=IBO(I22)
                                    IBO(I22)=IBO(I22+1)
                                    IBO(I22+1)=K
                                END IF
                            END DO
                        END DO

                        DO I=1,4
                            !WRITE(*,*)I,IBO(I),FCORL(IBO(I),3)
                        END DO

                        D1=0.4

                        IF(IBO(1).EQ.1) THEN
                            IF(IBO(2).EQ.3) THEN
                                IUP1(1)=7;IUP1(2)=6;IUP1(3)=5
                                IUP2(1)=8;IUP2(2)=9;IUP2(3)=4
                                XI(1)=-1.;ETA(1)=D1
                                XI(2)=0.;ETA(2)=D1
                                XI(3)=1.;ETA(3)=D1
                                XI(4)=-1.;ETA(4)=(-1.+D1)/2.
                                XI(5)=0.;ETA(5)=(-1.+D1)/2.
                                XI(6)=1.;ETA(6)=(-1.+D1)/2.
                            END IF
                            IF(IBO(2).EQ.7) THEN
                                IUP1(1)=3;IUP1(2)=4;IUP1(3)=5
                                IUP2(1)=2;IUP2(2)=9;IUP2(3)=6
                                ETA(1)=-1.;XI(1)=D1
                                ETA(2)=0.;XI(2)=D1
                                ETA(3)=1.;XI(3)=D1
                                ETA(4)=-1.;XI(4)=(-1.+D1)/2.
                                ETA(5)=0.;XI(5)=(-1.+D1)/2.
                                ETA(6)=1.;XI(6)=(-1.+D1)/2.
                            END IF
                        ELSE IF(IBO(1).EQ.3) THEN
                            IF(IBO(2).EQ.5) THEN
                                IUP1(1)=1;IUP1(2)=8;IUP1(3)=7
                                IUP2(1)=2;IUP2(2)=9;IUP2(3)=6
                                ETA(1)=-1.;XI(1)=-D1
                                ETA(2)=0.;XI(2)=-D1
                                ETA(3)=1.;XI(3)=-D1
                                ETA(4)=-1.;XI(4)=(1.-D1)/2.
                                ETA(5)=0.;XI(5)=(1.-D1)/2.
                                ETA(6)=1.;XI(6)=(1.-D1)/2.
                            END IF
                            IF(IBO(2).EQ.1) THEN
                                IUP1(1)=7;IUP1(2)=6;IUP1(3)=5
                                IUP2(1)=8;IUP2(2)=9;IUP2(3)=4
                                XI(1)=-1.;ETA(1)=D1
                                XI(2)=0.;ETA(2)=D1
                                XI(3)=1.;ETA(3)=D1
                                XI(4)=-1.;ETA(4)=(-1.+D1)/2.
                                XI(5)=0.;ETA(5)=(-1.+D1)/2.
                                XI(6)=1.;ETA(6)=(-1.+D1)/2.
                            END IF
                        ELSE IF(IBO(1).EQ.5) THEN
                            IF(IBO(2).EQ.7) THEN
                                IUP1(1)=1;IUP1(2)=2;IUP1(3)=3
                                IUP2(1)=8;IUP2(2)=9;IUP2(3)=4
                                XI(1)=-1.;ETA(1)=-D1
                                XI(2)=0.;ETA(2)=-D1
                                XI(3)=1.;ETA(3)=-D1
                                XI(4)=-1.;ETA(4)=(1.-D1)/2.
                                XI(5)=0.;ETA(5)=(1.-D1)/2.
                                XI(6)=1.;ETA(6)=(1.-D1)/2.
                            END IF
                            IF(IBO(2).EQ.3) THEN
                                IUP1(1)=1;IUP1(2)=8;IUP1(3)=7
                                IUP2(1)=2;IUP2(2)=9;IUP2(3)=6
                                ETA(1)=-1.;XI(1)=-D1
                                ETA(2)=0.;XI(2)=-D1
                                ETA(3)=1.;XI(3)=-D1
                                ETA(4)=-1.;XI(4)=(1.-D1)/2.
                                ETA(5)=0.;XI(5)=(1.-D1)/2.
                                ETA(6)=1.;XI(6)=(1.-D1)/2.
                            END IF
                        ELSE IF(IBO(1).EQ.7) THEN
                            IF(IBO(2).EQ.1) THEN
                                IUP1(1)=3;IUP1(2)=4;IUP1(3)=5
                                IUP2(1)=2;IUP2(2)=9;IUP2(3)=6
                                ETA(1)=-1.;XI(1)=D1
                                ETA(2)=0.;XI(2)=D1
                                ETA(3)=1.;XI(3)=D1
                                ETA(4)=-1.;XI(1)=(-1.+D1)/2.
                                ETA(5)=0.;XI(2)=(-1.+D1)/2.
                                ETA(6)=1.;XI(3)=(-1.+D1)/2.
                            END IF
                            IF(IBO(2).EQ.5) THEN
                                IUP1(1)=1;IUP1(2)=2;IUP1(3)=3
                                IUP2(1)=8;IUP2(2)=9;IUP2(3)=4
                                XI(1)=-1.;ETA(1)=-D1
                                XI(2)=0.;ETA(2)=-D1
                                XI(3)=1.;ETA(3)=-D1
                                XI(4)=-1.;ETA(1)=(1.-D1)/2.
                                XI(5)=0.;ETA(2)=(1.-D1)/2.
                                XI(6)=1.;ETA(3)=(1.-D1)/2.
                            END IF
                        END IF
                               
                        CORL2=FCORL       
                        CORL1=FCORL      
                                
                        DO I11=1,3
                            CALL ISOPAR_9(XI(I11),ETA(I11),CORL2,SN,SNX,SNE,AJ,SR,DN)
                            CORL1(IUP1(I11),:)=SR(:)
                        END DO

                        DO I11=4,6
                            CALL ISOPAR_9(XI(I11),ETA(I11),CORL2,SN,SNX,SNE,AJ,SR,DN)
                            CORL1(IUP2(I11-3),:)=SR(:)
                        END DO

                        DO I11=1,9
                            FCORL(I11,:)=CORL1(I11,:)
                        END DO

                        IF(ISYM.EQ.1) THEN
                            DO I11=1,9
                                CORD11(INE(IE,I11),:)=FCORL(I11,:)
                                !WRITE(*,*)CORD11(INE(IE,I11),:)
                            END DO
                            !WRITE(*,*)"OKKKKKKKKKKKKKKKKKKKKKKKK"
                        END IF

                        !STOP

                        GOTO 11

                        600 CONTINUE
                    ELSE
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM)
                        
                        IF(ISYM.GT.2) THEN
                            FCORL(I1,3)=CORD(INE(IE,I2),3)*ISM(3,ISYM)-2.*NSWIN(13,1)
                        END IF
                    END IF

                    !FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)
                    !FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM)

                    !STAGE=1 MIRROR OF FREE SURFACE    STAGE=2 MIRROR OF BOTTOM
                END DO

                !IF(IP.GT.NBPOINT+(NFXF+NFX+1)*NFY) THEN    !配置点前移
                IF(INE(IE,I2).GT.NBPOINT) THEN
                !IF(INE(IE,I2).GT.NBPOINT.AND.IP.GT.NBPOINT) THEN
                    !LOFF=0.5*DMIP(INE(IE,I1))**BETA
                    !LOFF=0.4*DFX !DAWSON OK
                    LOFF=0.25*DFX !LINEAR OK
                    !LOFF=1.0*DFX !wigley OK
                    !LOFF=0.8*DFX !yuzheng OK
                    !LOFF=0.25*DFX
                    !LOFF=0.1*DFX !LINEAR OK
                    !LOFF=0.
                    !LOFF=0.6*DFX 
                    !LOFF=0.4*DFX
                    FCORL(I1,1)=CORD(INE(IE,I2),1)-LOFF
                    !FPS(1)=CORD(IP,1)+LOFF
                    !LOFF=0.6*DFX
                END IF

                !物面水线节点内收
                IF(ABS(CORD(INE(IE,I2),3)).LE.1.E-4.AND.INE(IE,I2).LE.NBPOINT) THEN
                    !LOFF=SQRT(CORD(INE(IE,I2),1)**2+CORD(INE(IE,I2),2)**2)
                    !FCORL(I1,1)=(LOFF-0.5*DFX)*CORD(INE(IE,I2),1)/LOFF
                    !FCORL(I1,2)=(LOFF-0.5*DFX)*CORD(INE(IE,I2),2)/LOFF

                    !FCORL(I1,1)=CORD(INE(IE,I2),1)*0.85
                    !FCORL(I1,2)=CORD(INE(IE,I2),2)*0.85
                    !FCORL(I1,3)=CORD(INE(IE,I2),3)-RD/20.
                END IF
            END DO

            11 continue


!			判断是否是奇点以调用不同的程序，判断标志为ISINGULAR
            !WRITE(*,*)"ps",ps
            !WRITE(*,*)"fcorl",fcorl(1:4,:)
			ISINGULAR=0
			DO I1=1,PAN(IE)
				T=0.
				DO J1=1,3
    				T=T+(FCORL(I1,J1)-FPS(J1))**2
                END DO
				T=SQRT(T)
                IF(T.LT.1.E-5) THEN
					ISINGULAR=1
				END IF
            END DO

!c			是奇点则调用cofahs，否则调用cofah,aij为源效应，hij为偶极子效应 
			DO I=1,PAN(IE)
				FHIJ(I,:)=0.0
				AIJ(I)=0.0
            END DO
              					
			IF(ISINGULAR.EQ.1) THEN
                !IF(PAN(IE).EQ.4) CALL COFAHS_4(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),HIJ(1:PAN(IE)))
                !IF(PAN(IE).EQ.9) CALL COFAHS_9NN(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                IF(IP.GT.NBPOINT.AND.ABS(CORD(IP,3)).LE.1.E-5) THEN  
                    CALL COFAHS_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                    !CALL COFAHN_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FHIJ(1:PAN(IE),:))
                ELSE
                    CALL COFAHS_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                END IF
                !IF(PAN(IE).EQ.8) THEN
                !    CALL COFAHS_8(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                    !WRITE(*,*)IP,IE,PAN(IE)
                !END IF
                !IF(IP.EQ.1.AND.IE.LE.1) WRITE(*,*)IP,IE,HIJ(1:PAN(IE))
			ELSE IF(ISINGULAR.EQ.0) THEN
                !IF(PAN(IE).EQ.4) CALL COFAH_4(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),HIJ(1:PAN(IE)))
                !CALL COFAH_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                CALL COFAH_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                !FHIJ=0.
            ENDIF

            IF(IP+K1.EQ.1) THEN
                !CALL COFAH_4SN(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),XGRA,FAIJ)   !与PHIS_T的边界积分方程有关
            END IF

! 			将每一面元中各顶点处的积分值找到对应的总序列中的位置，并进行叠加
! 			SSH(k)中存放偶极子效应
! 			SSA(k)中存放界面上的源效应
!			SUMHIJ(NTPN)用于存放对角线上的量

			DO I1=1,PAN(IE)
                !WRITE(*,*)IE
				K=INE(IE,I1)
				I2=I1
                IF(PAN(IE).EQ.3) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=5-I1	
                END IF                
                IF(PAN(IE).EQ.4) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				    IF(MOD(ISYM,4).EQ.0.AND.(I1.NE.1)) I2=6-I1	
                END IF
                IF(PAN(IE).EQ.8) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                END IF
                IF(PAN(IE).EQ.9) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                END IF


                TL=SQRT((CORD(IP+K1,1)-CORD(K,1))**2+(CORD(IP+K1,2)-CORD(K,2))**2+(CORD(IP+K1,2)-CORD(K,2))**2)
                
                !IF(K.NE.IP+K1)THEN
                    !IF(TL.GT.1.E-6) THEN
                        !IF(IP+K1.LE.NTPN) THEN
				            !SSH(K)=SSH(K)+HIJ(I2)
                
                !IF(ISYM.LE.2.OR.IE.LE.NBBLOCK) THEN
                !IF(IP.LE.NBPOINT.OR.IE.GT.NBBLOCK) THEN
                            RSA(IP+K1,K)=RSA(IP+K1,K)+AIJ(I2)
                
                !IF(K.NE.IP+K1)THEN
                            RSH(IP+K1,K,1)=RSH(IP+K1,K,1)+FHIJ(I2,1)   !D(1/R)/D(X)
                            RSH(IP+K1,K,2)=RSH(IP+K1,K,2)&
                                           +FHIJ(I2,1)*VECN(IP+K1,1)&
                                           +FHIJ(I2,2)*VECN(IP+K1,2)&
                                           +FHIJ(I2,3)*VECN(IP+K1,3)   !D(1/R)/D(N)
                            RSH(IP+K1,K,3)=RSH(IP+K1,K,3)+FHIJ(I2,3)   !D(1/R)/D(Z)
                            RSH(IP+K1,K,4)=RSH(IP+K1,K,4)+FHIJ(I2,4)   !D(1/R)/D(Y)
                !END IF
                        !END IF
				    !END IF
                !END IF
				    
                IF(K.EQ.IP+K1)THEN
                !IF(K.NE.IP+K1)THEN
                !IF(TL.GT.1.E-6) THEN
                !IF(TL.LE.1.E-6) THEN
					SUMHIJ(IP+K1)=SUMHIJ(IP+K1)+FHIJ(I2,1)*VECN(IP+K1,1)&
                                               +FHIJ(I2,2)*VECN(IP+K1,2)&
                                               +FHIJ(I2,3)*VECN(IP+K1,3)   !D(1/R)/D(N) 
                    !SUMHIJ(IP+K1)=SUMHIJ(IP+K1)+FHIJ(I2,2)                           
                END IF

                !END IF
                !IF(IP.EQ.1) WRITE(*,*)SSA(K),SSH(K)
                IF(ISYM.LE.2) THEN
                    !RSH(IP+K1,K)=SSH(K)
                    !RSA(IP+K1,K)=SSA(K)
                END IF
                

                IF(IP.EQ.NBPOINT+NFXF*NFY+1.AND.K.EQ.1.AND.ISYM.EQ.1) THEN
                    !WRITE(*,*)IP,K,AIJ(I2)
                    !WRITE(*,*)IP,K,FHIJ(I2,2)
                    !WRITE(*,*)CORD(IP,1:3)
                    !WRITE(*,*)CORD(K,1:3)
                    !STOP
                END IF   
                      
            END DO
        END DO
        END DO
    END DO
!!$OMP END PARALLEL

    IF(STAGE.EQ.2) THEN
        NPLOOP=NBPOINT+NFPOINT
    END IF

    DO I=1,NPLOOP
        DO J=1,NPLOOP
        	!WRITE(18)RSA(I,J),RSH(I,J,1),RSH(I,J,2),RSH(I,J,3),RSH(I,J,4)   
            SA(I,J)=RSA(I,J)
            SX(I,J)=RSH(I,J,1)
            SH(I,J)=RSH(I,J,2)
            SZ(I,J)=RSH(I,J,3)
            SY(I,J)=RSH(I,J,4)
        END DO
        !WRITE(118,111)SUMHIJ(I),RSH(I,1,1),RSH(I,2,1),RSH(I,3,1),RSA(I,2),CORD(I,1:3) !,RSH(I,I,2),RSH(I,I,4)
        111FORMAT(10F15.6)
        !WRITE(*,*)SUMHIJ(I)
    END DO    


!c	释放内存
	DEALLOCATE(PS,CORL,CNL,SSH1,SSH,SSA,SSA1,HIJ,AIJ,RSH,RSA,RFSA,HIJ1,AIJ1,CORD11)
	WRITE(*,*)'		END COEFMATRIX'

    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
    !PRINT *, 'END   TIME  COEFMATRIX: ', CHAR_TIME2,'  ON  ',CHAR_DATE2    

	CLOSE(18)
    CLOSE(118)
    CLOSE(4)
    RETURN
END SUBROUTINE

SUBROUTINE COFAH_4(CORDA,CONM,PS,AIJ,HIJ)
!    ***********************************************      
!     CALCULATE THE AIJ AND HIJ TERMS IN THE MATRIX     
!      CORDA(4,3) THE CORD. OF 4 NODES (CLOCKWISE)      
!      CONM(4,3) THE NORMAL VECTOR OF 4 NODES          4---3
!      PS(3) THE CORD. OF FIELD POINT                  |   |
!      AIJ(4): PANEL SOURCE EFFECT                     1---2
!      HIJ(4): PANEL DIPOLE EFFECT
!   ************************************************
	USE NMLMOD
    REAL:: XI,ETA,AJ
    REAL:: CORDA(4,3),PS(3),CONM(4,3),AIJ(4),HIJ(4,8),PR(3)
    REAL:: R(4),SRR(3),AY(4),HY(4,8)
    REAL:: SN(4),SNX(4),SNE(4),SR(3),DN(3)
    COMMON /GWA/GA(400),GW(400)

!	初始化每一面元的四个顶点的源效应AIJ（4）和偶极子效应HIJ（4）
    DO I=1,4
        AIJ(I)=0.
        HIJ(I,:)=0.
    END DO  
!	计算四个顶点的平均坐标，并计算此点到源点的距离R1
	R1=0
    DO I=1,3
		PR(I)=0.
		DO J=1,4
			PR(I)=PR(I)+CORDA(J,I)
        END DO
		PR(I)=PR(I)/4.-PS(I)
		R1=R1+PR(I)*PR(I)
    END DO
      
	R1=SQRT(R1)      
!	分别计算四条边的长度，记在R（4）中     
    DO I=1,4
        R(I)=0.
    END DO

	DO J=2,4
	    DO I=1,3
			R(J)=R(J)+(CORDA(J,I)-CORDA(J-1,I))**2
        END DO
    END DO
	DO I=1,3
        R(1)=R(1)+(CORDA(4,I)-CORDA(1,I))**2
    END DO

!	求得平均边长RV，并根据R0＝R1/RV来确定选用高斯点的个数
    RV=0.
    DO I=1,4
        RV=RV+SQRT(R(I))
    END DO
    RV=RV/4.
    R0=R1/RV

    NGUS=21
    IF(R0.LT.1.6) GOTO 30
    NGUS=11
    IF(R0.LT.1.9) GOTO 30
    NGUS=9
    IF(R0.LT.2.4) GOTO 30
    NGUS=7
    IF(R0.LT.3.4) GOTO 30
    NGUS=5
    IF(R0.LT.6.4) GOTO 30
    NGUS=3
    30 CONTINUE

!	K=2+3+...+(N-1)+1=N*(N-1)/2
    K=NGUS*(NGUS-1)/2
!   GAUSS-LEGENDRE NUMERICAL INTEGRATION

!	对X积分     
    DO IX=1,NGUS
		XI=GA(K+IX-1)
!	对Y积分
		DO IY=1,NGUS
			ETA=GA(K+IY-1)
!	调用ISOPAR，求得形状函数J及单位法向矢量等     
			CALL ISOPAR_4(XI,ETA,CORDA,SN,SNX,SNE,AJ,SR,DN)

!      
			RH=0.
			DO I=1,3
				SRR(I)=SR(I)-PS(I)   
   				RH=RH+SRR(I)*SRR(I)
            END DO
			RH=SQRT(RH)
			!RH3=0.0-RH**3
			RH3=RH**3
     
!			NML=1 使用ISOPAR的曲线平面上的精确法向矢量，
!               否则使用四节点单元的等参法向矢量
			IF(NML.NE.1) THEN
				DO I1=1,3  
					DN(I1)=0.
					DO I2=1,4	
   						DN(I1)=DN(I1)+CONM(I2,I1)*SN(I2)
    				END DO
	            END DO
    		END IF
!       
			BNR=SRR(1)*DN(1)+SRR(2)*DN(2)+SRR(3)*DN(3)
   
			DO J=1,4
				!AY(J)=GW(K+IY-1)*SN(J)*AJ/RH
   				!HY(J)=GW(K+IY-1)*SN(J)*BNR*AJ/RH3
                
                AY(J)=GW(K+IY-1)*SN(J)*AJ/RH
                HY(J,1)=GW(K+IY-1)*SN(J)*SRR(1)*AJ/RH3   !DX
                HY(J,2)=GW(K+IY-1)*SN(J)*SRR(2)*AJ/RH3   !DY      
                HY(J,3)=GW(K+IY-1)*SN(J)*SRR(3)*AJ/RH3   !DZ
           
                HY(J,4)=GW(K+IY-1)*SN(J)*3.*SRR(1)**2*AJ/RH**5&
                    -GW(K+IY-1)*SN(J)*AJ/RH**3    !DXX

        !HY(J,4)=3.*SRR(1)**2*AJ/(RH**5)
        !HY(J,4)=HY(J,4)*GW(K+IY-1)*SN(J)-GW(K+IY-1)*SN(J)*AJ/(RH3)

        HY(J,5)=3.*SRR(1)*SRR(2)*AJ/(RH**5)    !DXY
        HY(J,5)=HY(J,5)*GW(K+IY-1)*SN(J)

        HY(J,6)=3.*SRR(1)*SRR(3)*AJ/(RH**5)    !DXZ
        HY(J,6)=HY(J,6)*GW(K+IY-1)*SN(J)   

        HY(J,7)=GW(K+IY-1)*SN(J)*3.*SRR(2)**2*AJ/(RH**5)&
                -GW(K+IY-1)*SN(J)*AJ/RH**3    !DYY
        !HY(J,7)=3.*SRR(2)**2*AJ/(RH**5)
        !HY(J,7)=HY(J,7)*GW(K+IY-1)*SN(J)-GW(K+IY-1)*SN(J)*AJ/(RH3)    !DYY

        HY(J,8)=3.*SRR(2)*SRR(3)*AJ/(RH**5)    !DYZ
        HY(J,8)=HY(J,8)*GW(K+IY-1)*SN(J)    !DYZ                
                
            END DO
!
			DO J=1,4
				!AIJ(J)=AIJ(J)+GW(K+IX-1)*AY(J)
				!HIJ(J)=HIJ(J)+GW(K+IX-1)*HY(J)
                

        AIJ(J)=AIJ(J)+GW(K+IX-1)*AY(J)
        HIJ(J,1)=HIJ(J,1)+GW(K+IX-1)*HY(J,1)
        HIJ(J,2)=HIJ(J,2)+GW(K+IX-1)*HY(J,2)
        HIJ(J,3)=HIJ(J,3)+GW(K+IX-1)*HY(J,3)
        HIJ(J,4)=HIJ(J,4)+GW(K+IX-1)*HY(J,4)

        HIJ(J,5)=HIJ(J,5)+GW(K+IX-1)*HY(J,5)
        HIJ(J,6)=HIJ(J,6)+GW(K+IX-1)*HY(J,6)
        HIJ(J,7)=HIJ(J,7)+GW(K+IX-1)*HY(J,7)
        HIJ(J,8)=HIJ(J,8)+GW(K+IX-1)*HY(J,8)
            END DO
        END DO
    END DO
    
    RETURN
END SUBROUTINE

SUBROUTINE COFAHS_4(CORDA,CONM,PS,AIJ,HIJ)
!   ***********************************************      
!   CALCULATE THE AIJ AND HIJ TERMS IN THE MATRIX   
!   THE PS IS ONE OF THE NODES OF ELEMENT          
!   CORD(4,3) THE CORD. OF 8 NODES (CLOCKWISE)     
!   CONM(4,3) THE NORMAL VECTOR OF 8 NODES         4---3
!   PS(3) THE CORD. OF FIELE POINT                 |   |
!   AIJ(4): PANEL SOURCE EFFECT                    1---2
!   HIJ(4): PANEL DIPOLE EFFECT
!   ************************************************
!	使用了蔡瑞英，曾昭景，黄文龙，《边界元法程序涉及及其工程应用》中方法
	USE NMLMOD
!   COMMON NML
    REAL:: CORDA(4,3),CONM(4,3),CORDB(4,3),CONM1(4,3),PS(3)
    REAL:: AIJ(4),HIJ(4,8),AIJ1(4),HIJ1(4,8)
    REAL:: AJ(3)
    REAL:: SNX(4),SNE(4)
    REAL:: SN1(4)
    REAL:: SR1(3)
    REAL:: DN1(3)
    COMMON /GWA/GA(400),GW(400)

	REAL:: XL,FL,WX,WY,FL1,XW,YW,XQP,YQP,ZQP,RQ2,RQ1,FL2,ATT,HTT
	REAL:: FJ(2),XT(2),YT(2)

!	初始化    
    R1=0.
    ISN=0.
	PI=3.1416
!	分别计算源点PS(3)到此面元四个顶点的距离，以判断奇点的位置
    DO I=1,4
	    T=(CORDA(I,1)-PS(1))**2+(CORDA(I,2)-PS(2))**2&
     	  +(CORDA(I,3)-PS(3))**2
		T=SQRT(T)
	    IF(T.LT.1.E-6) ISN=I
	    AIJ1(I)=0.
	    HIJ1(I,:)=0.
	    AIJ(I)=0.0
	    HIJ(I,:)=0.0
    END DO
!	如果没有找到奇点，提示出错并停止程序     
    IF(ISN.EQ.0) THEN
	    WRITE(*,*)"ISN",ISN,T
        WRITE(*,*)CORDA
        WRITE(*,*)
        WRITE(*,*)PS
        PAUSE 'WARNING: ISN=0, THE SOURCE POINT OUT OF THE ELEMENT'
	    STOP
    END IF	
	
	DO I=1,4
!         将奇点转到第一点，如奇点在1，则CORDB的1234对应于CORDA中的1，2，3，4，
!                           若奇点在2，则CORDB的1234对应于CORDA中的2，3，4，1
!                           若奇点在3，则CORDB的1234对应于CORDA中的3，4，1，2
!                           若奇点在4，则CORDB的1234对应于CORDA中的4，1，2，3
!         以便在下面的积分中将正方形区域转化为三角形区域
        K=MOD(ISN+I-2,4)+1
	    DO J=1,3
            CORDB(I,J)=CORDA(K,J)
	        CONM1(I,J)=CONM(K,J)
        END DO
    END DO

!   G-L数值积分，使用11个GAUSS点
    NGUS=11
    KGUS=NGUS*(NGUS-1)/2

!	对KSI积分
	DO IX=1,NGUS
		XL=GA(KGUS+IX-1)
		WX=GW(KGUS+IX-1)
!		对ETA积分
		DO IY=1,NGUS
			YL=GA(KGUS+IY-1)
			WY=GW(KGUS+IY-1)

			FL1=TAN(0.125*PI*(1.0+XL))
			FJ(1)=0.125*(1.0+YL)*PI*(1.0+FL1*FL1) !DE1DE2
			FJ(2)=FJ(1)
			
			XT(1)=YL                  !E1
			YT(1)=(1.0+YL)*FL1-1      !E2
			XT(2)=YT(1)
			YT(2)=XT(1)
				
			DO K=1,2
				XW=XT(K)
				YW=YT(K)
! 				分别就两种情形调用ISOPAR，求得形状函数J和单位法向矢量等
				CALL ISOPAR_4(XW,YW,CORDB,SN1,SNX,SNE,AJ(1),SR1,DN1)
				COSBX=DN1(1)
				COSBY=DN1(2)
				COSBZ=DN1(3)

				XQP=SR1(1)-PS(1)
				YQP=SR1(2)-PS(2)
				ZQP=SR1(3)-PS(3)

				RQ2=XQP*XQP+YQP*YQP+ZQP*ZQP
				RQ1=SQRT(RQ2)
				
!C				FL1=WX*WY*AJ(1)*FJ(K)/(4*PI)
				FL1=WX*WY*AJ(1)*FJ(K)
				FL2=(XQP*COSBX+YQP*COSBY+ZQP*COSBZ)/RQ1

				HTT=-FL1*FL2/RQ2
				ATT=FL1/RQ1
				DO J=1,4
					IF(J.NE.1)THEN
						!HIJ1(J)=HIJ1(J)+SN1(J)*HTT
                        
                        HIJ1(J,1)=HIJ1(J,1)+SN1(J)*(FL1*XQP/RQ1/RQ2)   !DX
                        HIJ1(J,2)=HIJ1(J,2)+SN1(J)*(FL1*YQP/RQ1/RQ2)    !DY      
                        HIJ1(J,3)=HIJ1(J,3)+SN1(J)*(FL1*ZQP/RQ1/RQ2)    !DZ
					ENDIF
					AIJ1(J)=AIJ1(J)+SN1(J)*ATT
                END DO
            END DO
        END DO
    END DO

!   ASSEMBLE THE LOCAL COEF. MATRIX
!	求得最终的AIJ(4)和HIJ(4)
	DO J=1,4
		K=MOD(ISN+J-2,4)+1
		!AIJ(K)=AIJ1(J)
		!HIJ(K)=HIJ1(J)   
        
        AIJ(K)=AIJ1(J)
        HIJ(K,1)=HIJ1(J,1)
        HIJ(K,2)=HIJ1(J,2)
        HIJ(K,3)=HIJ1(J,3)
        HIJ(K,4)=HIJ1(J,4)

        HIJ(J,5)=HIJ(J,5)+GW(K+IX-1)*HIJ1(J,5)
        HIJ(J,6)=HIJ(J,6)+GW(K+IX-1)*HIJ1(J,6)
        HIJ(J,7)=HIJ(J,7)+GW(K+IX-1)*HIJ1(J,7)
        HIJ(J,8)=HIJ(J,8)+GW(K+IX-1)*HIJ1(J,8)        
    END DO
      
	RETURN
END SUBROUTINE

SUBROUTINE COEFMATRIX_PARA_1
!   ***************************************************
!   GENERATE THE GLOBAL MATRIX OF 1st INTEGRAL EQUATION  
!   ***************************************************
    
    USE GREENMOD

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

    !IMPLICIT NONE
    CHARACTER MARK
	INTEGER:: ISM(3,4)
	DATA ISM/1,1,1,1,-1,1,1,-1,-1,1,1,-1/
	INTEGER:: NLOOP,K,K1,K2,CP,NBLOOP,NPLOOP,IP,I3,I1,J1,IE,J,ISYM,ISINGULAR,I,I2,I11,I22,J11,IAF(2),CAF,IBO(4),IUP1(3),IUP2(3)
    REAL ::FCORL(9,3),FCNL(9,3),FPS(3),T,TL,FAIJ(9),BETA,LOFF,CORL1(9,3),CORL2(9,3),SR(3),D1,D2,D3,D4,XI(6),ETA(6)
    REAL ::HIJX(9),HIJN(9),HIJZ(9),HYX(9),HYN(9),HYZ(9),FHIJ(9,8),VV(3) !诱导系数 
    REAL ::SN(9),SNX(9),SNE(9),AJ,DN(3)
    REAL,ALLOCATABLE:: CORD11(:,:)
    
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2
    
    ALLOCATE(CORD11(NTPN+500,3))
	ALLOCATE(SH(NTPN+500,NTPN+500),SX(NTPN+500,NTPN+500),SZ(NTPN+500,NTPN+500),SA(NTPN+500,NTPN+500),SA1(NTPN+500,NTPN+500),SY(NTPN+500,NTPN+500))
	IF(NIT.EQ.1.AND.ITTE.EQ.1) THEN
        ALLOCATE(SXX(NTPN+500,NTPN+500),SXY(NTPN+500,NTPN+500),SXZ(NTPN+500,NTPN+500),SYY(NTPN+500,NTPN+500),SYZ(NTPN+500,NTPN+500))
    END IF

    CORD11=CORD

    !WRITE(*,*)"DFY=     ",DFY,0.6*DFX
    !WRITE(*,*) 
    !WRITE(*,*)'SUB COEFMATRIX..............'

    open(18,file='coefuse.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(18)
    open(19,file='coefuse_T.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(19)
    
    open(118,file='coefuse.DAT')

    NRPOINT=0.

    ALLOCATE (IDS(NTPN))
	ALLOCATE (PS(3),CORL(9,3),CNL(9,3))
	ALLOCATE (SSH1(NTPN+NRPOINT),SSH(NTPN+NRPOINT),SSA(NTPN+NRPOINT),SSA1(NTPN+NRPOINT),HIJ(9),AIJ(9),AIJ1(9),HIJ1(9))
	ALLOCATE (SUMHIJ(NTPN+NRPOINT))
    !WRITE(*,*)"NTPN+NRPOINT",NTPN,NRPOINT
    ALLOCATE (RSH(NTPN+NRPOINT,NTPN+NRPOINT,4),RSA(NTPN+NRPOINT,NTPN+NRPOINT),RFSA(NBPOINT,6))

    OPEN(4,FILE='INDUCED_VELO.DAT')

    RSH=0.;RSA=0.;RFSA=0.
    SX=0.;SY=0.;SZ=0.;SXX=0.;SXY=0.;SXZ=0.;SYY=0.;SYZ=0.;SH=0.

    m=nbpoint
    L=NFPOINT
    NRPOINT=0.
    NRBLOCK=0.
!  NTPN+NRPOINT IN SINGLE CPU
   IF(MOD(M+L+NRPOINT,NTHREAD).NE.0) THEN
     DO I=1,NTHREAD-1
       IF(MOD(M+L+NRPOINT,NTHREAD).EQ.I) THEN
         DO J=1,NTHREAD-1
            NPS1(J)=(M+L+NRPOINT-I)/NTHREAD
         END DO
         NPS1(NTHREAD)=NPS1(NTHREAD-1)+I
         EXIT
       END IF
     END DO
   ELSE 
     DO I=1,NTHREAD
       NPS1(I)=(M+L+NRPOINT)/NTHREAD
     END DO
   END IF

    DO IP=1,NTPN
		SUMHIJ(IP)=0.0
    END DO

    IF(STAGE.EQ.2) THEN
        NPLOOP=NTPN
        NBLOOP=NTP
    END IF

    IF(STAGE.EQ.1) THEN
        NPLOOP=NBPOINT
        NBLOOP=NBBLOCK
    END IF

    !K1=0
    CP=1

    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
    !PRINT *, 'END   TIME  DISTRIBUTE: ', CHAR_TIME2,'  ON  ',CHAR_DATE2

!START PARALLEL
    CALL OMP_SET_NUM_THREADS(NTHREAD)


!!!$omp parallel private(CP,k1,k2,FCORL,FCNL,ISYM,FPS,T,HIJ,AIJ,IP,I3,IE,I1,I2,J1,ISINGULAR,TL,K,FAIJ,FHIJ)

    
    CP=OMP_GET_THREAD_NUM()+1

    K1=SUM(NPS1(1:CP))-NPS1(CP) 
    K2=SUM(NPS1(1:CP))
    K1=0
    !WRITE(*,*)"K1   ",K1
    !write(*,*)NPS1(CP)

    !ALLOCATE (HIJ(9),AIJ(9))

    IDS=0
    DO IP=1,NPLOOP
    !DO IP=1,NPS1(CP)

        !LOFF=LDC*DMIP(IP)**BETA
		!IF((((IP+K1)/(ntpn/4))*(NTPN/4)).EQ.IP+K1)THEN
			!WRITE(*,200)IP+K1,NPLOOP !NTPN
		!END IF
  200	FORMAT('16H    MOVING BALL',3X,'6HINODE=',I8,5X,I8)
   
!		INITIAL THE DIMENSIONS
!		SSH(K)中存放有奇点部分的偶极子效应
!		SSA(K)中存放有奇点部分的源效应
        DO 5 I3=1,NPLOOP
			SSH(I3)=0.
            SSH1(I3)=0.
			SSA(I3)=0.
			SSA1(I3)=0.0
    5		CONTINUE						  
!		将每一个点的坐标赋值到ps（3）中
        
        BETA=0.5
        !LOFF=DMIP(IE)**BETA  !单元内潜上置距离  
        DO 10 J1=1,3
   10		FPS(J1)=CORD(IP+K1,J1) !+VECN(IP,J1)*LOFF
!		LOOP FOR ELEMENTS
	    !FPS(3)=FPS(3)+0.05
    	DO IE=1,NBLOOP !+(NFY-1)/2
!			WRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!			依次将每一面元的PAN(IE)个顶点的坐标赋值到corl(PAN(IE),3)中
!			依次将每一面元的PAN(IE)个顶点的法向矢量赋值到cnl(PAN(IE),3)中
			
            IF(MNKDB.EQ.2) THEN
                IF(IE.LE.NBBLOCK) THEN
                    KK=4
                ELSE
                    KK=2
                END IF
            ELSE IF(MNKDB.EQ.1) THEN
                KK=2
            END IF

            DO ISYM=1,KK
				
			DO I1=1,PAN(IE)
				I2=I1
                IF(PAN(IE).EQ.3) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=5-I1	
                END IF                
                IF(PAN(IE).EQ.4) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				    IF(MOD(ISYM,4).EQ.0.AND.(I1.NE.1)) I2=6-I1	
                END IF
                IF(PAN(IE).EQ.8) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                END IF
                IF(PAN(IE).EQ.9) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                END IF

                DO J1=1,3
                    IF(IE.GT.NBBLOCK) THEN
                        !LOFF=0.8*DFX
                        
                        LOFFZ=1.0*DFX
                        !LOFFZ=1.2*DFX
                        IF(FR.GE.0.3) THEN
                            LOFFZ=1.4*DFX
                        ELSE IF(FR.GE.0.6) THEN
                            LOFFZ=1.2*DFX
                        ELSE IF(FR.GE.0.7) THEN
                            LOFFZ=1.2*DFX
                        ELSE
                            LOFFZ=1.2*DFX
                        END IF
                        !LOFFZ=0.4*DFY
                        !LOFFZ=1.0*DFX
                    ELSE
                        LOFF=DFY
                        !LOFF=1.3*DFY
                    END IF
                    LOFFZ=1.5*DFX !KCS
                    !WRITE(*,*)DFX !DMIE(IE),BETA, LOFF
                    
                    IF(IE.GT.NBBLOCK) THEN
                    !IF(ABS(CORD(INE(IE,I2),2)).GE.1.E-6) THEN
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)+VECN(INE(IE,I2),J1)*LOFFZ*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM) 

                        !FCORL(I1,2)=CORD(INE(IE,I2),2)*ISM(J1,ISYM)-LOFFY*ISM(J1,ISYM)

                        if(fr.lT.0.3) FCORL(I1,3)=LOFFZ
                        
                        FCORL(I1,3)=LOFFZ
                    ELSE
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM)
                    END IF
                    
                    !IF(INE(IE,I2).GT.NBPOINT+(NFXF+NFX)*NFY) THEN
                        !FCORL(I1,1)=CORD(INE(IE,I2),1)-LOFFX
                    !END IF    
                END DO

                IF(INE(IE,I2).GT.NBPOINT) THEN
                !IF(IP+K1.GT.NBPOINT) THEN
                !IF(INE(IE,I2).GT.NBPOINT+(NFXF+NFX)*NFY) THEN
                    !LOFF=0.8*DFX !yuzheng OK
                    !LOFFX=0.25*DFX
                    IF(MDIFF.EQ.1) THEN
                        LOFFX=1.0*DFX !LINEAR OK
                        LOFFX=0.8*DFX !LINEAR OK
                    ELSE
                        !LOFFX=0.4*DFX !LINEAR OK
                        !LOFFX=0.1*DFX !LINEAR OK
                        LOFFX=0.3*DFX !LINEAR OK
                        !LOFFX=0.4*DFX !LINEAR OK
                        !LOFFX=1.25*DFX !LINEAR OK
                        !write(*,*)LOFFX
                    END IF
                    
                    FCORL(I1,1)=CORD(INE(IE,I2),1)-LOFFX
                END IF

            END DO
            11 continue

!			判断是否是奇点以调用不同的程序，判断标志为ISINGULAR
            !WRITE(*,*)"ps",ps
            !WRITE(*,*)"fcorl",fcorl(1:4,:)
			ISINGULAR=0
			DO I1=1,PAN(IE)
				T=0.
				DO J1=1,3
    				T=T+(FCORL(I1,J1)-FPS(J1))**2
                END DO
				T=SQRT(T)
                IF(T.LT.1.E-5) THEN
					ISINGULAR=1
				END IF
            END DO
            
!c			是奇点则调用cofahs，否则调用cofah,aij为源效应，hij为偶极子效应 
            
			DO I=1,PAN(IE)
				!FHIJ(I,:)=0.0
            END DO
			!AIJ(I)=0.0

			IF(ISINGULAR.EQ.1) THEN
                !IF(PAN(IE).EQ.4) CALL COFAHS_4(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),HIJ(1:PAN(IE)))
                !IF(PAN(IE).EQ.9) CALL COFAHS_9NN(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                !IF(IP+K1.GT.NBPOINT) THEN  
                !    CALL COFAHS_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FAIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                !    CALL COFAHN_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FHIJ(1:PAN(IE),:))
                !ELSE
                    CALL COFAHS_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FAIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                    !FHIJ=0.
                !END IF
                !IF(PAN(IE).EQ.8) THEN
                !    CALL COFAHS_8(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                    !WRITE(*,*)IP,IE,PAN(IE)
                !END IF
                !IF(IP.EQ.1.AND.IE.LE.1) WRITE(*,*)IP,IE,HIJ(1:PAN(IE))
			ELSE IF(ISINGULAR.EQ.0) THEN
                !IF(PAN(IE).EQ.4) CALL COFAH_4(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),HIJ(1:PAN(IE)))
                CALL COFAH_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FAIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                !FHIJ=0.
            ENDIF
            !WRITE(*,*)IP


            IF(IP+K1.EQ.1) THEN
                !CALL COFAH_4SN(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),XGRA,FAIJ)   !与PHIS_T的边界积分方程有关
            END IF

! 			将每一面元中各顶点处的积分值找到对应的总序列中的位置，并进行叠加
! 			SSH(k)中存放偶极子效应
! 			SSA(k)中存放界面上的源效应
!			SUMHIJ(NTPN)用于存放对角线上的量

			DO I1=1,PAN(IE)
                !WRITE(*,*)IE
				K=INE(IE,I1)
				I2=I1
                IF(PAN(IE).EQ.3) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=5-I1	
                END IF                
                IF(PAN(IE).EQ.4) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				    IF(MOD(ISYM,4).EQ.0.AND.(I1.NE.1)) I2=6-I1	
                END IF
                IF(PAN(IE).EQ.8) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                END IF
                IF(PAN(IE).EQ.9) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                END IF

                DO KK=1,9
                    TL=SQRT((CORD(IP+K1,1)-CORD(INE(IE,KK),1))**2+(CORD(IP+K1,2)-CORD(INE(IE,KK),2))**2+(CORD(IP+K1,3)-CORD(INE(IE,KK),3))**2)
                    IF(TL.LE.1.E-6) THEN
                        EXIT
                    END IF
                END DO

                !IF(K.NE.IP+K1)THEN
                    !IF(TL.GT.1.E-6) THEN
                        !IF(IP+K1.LE.NTPN) THEN
				            !SSH(K)=SSH(K)+HIJ(I2)
                
                !IF(ISYM.LE.2.OR.IE.LE.NBBLOCK) THEN
                !IF(IP.LE.NBPOINT.OR.IE.GT.NBBLOCK) THEN
                IF(K.LE.NTPN) THEN
                !IF(ISYM.LE.2) THEN
                            RSA(IP+K1,K)=RSA(IP+K1,K)+FAIJ(I2)
                
                !IF(K.NE.IP+K1)THEN
                            RSH(IP+K1,K,1)=RSH(IP+K1,K,1)+FHIJ(I2,1)   !D(1/R)/D(X)
                            RSH(IP+K1,K,2)=RSH(IP+K1,K,2)&
                                           +FHIJ(I2,1)*VECN(IP+K1,1)&
                                           +FHIJ(I2,2)*VECN(IP+K1,2)&
                                           +FHIJ(I2,3)*VECN(IP+K1,3)   !D(1/R)/D(N)
                            RSH(IP+K1,K,3)=RSH(IP+K1,K,3)+FHIJ(I2,3)   !D(1/R)/D(Z)
                            RSH(IP+K1,K,4)=RSH(IP+K1,K,4)+FHIJ(I2,4)   !D(1/R)/D(Y)


                            SX(IP+K1,K)=SX(IP+K1,K)+FHIJ(I2,1)
                            SY(IP+K1,K)=SY(IP+K1,K)+FHIJ(I2,2)
                            SZ(IP+K1,K)=SZ(IP+K1,K)+FHIJ(I2,3)
                            SXX(IP+K1,K)=SXX(IP+K1,K)+FHIJ(I2,4)
                            SXY(IP+K1,K)=SXY(IP+K1,K)+FHIJ(I2,5)
                            SXZ(IP+K1,K)=SXZ(IP+K1,K)+FHIJ(I2,6)
                            SYY(IP+K1,K)=SYY(IP+K1,K)+FHIJ(I2,7)
                            SYZ(IP+K1,K)=SYZ(IP+K1,K)+FHIJ(I2,8)
                            SH(IP+K1,K)=SH(IP+K1,K)&
                                           +FHIJ(I2,1)*VECN(IP+K1,1)&
                                           +FHIJ(I2,2)*VECN(IP+K1,2)&
                                           +FHIJ(I2,3)*VECN(IP+K1,3)

                END IF
                        !END IF
				    !END IF
                !END IF
				    
                !IF(K.EQ.IP+K1)THEN
                IF(K.NE.IP+K1)THEN
                !IF(TL.GT.1.E-6) THEN
                !IF(TL.LE.1.E-6) THEN
					SUMHIJ(IP+K1)=SUMHIJ(IP+K1)+FHIJ(I2,1)*VECN(IP+K1,1)&
                                               +FHIJ(I2,2)*VECN(IP+K1,2)&
                                               +FHIJ(I2,3)*VECN(IP+K1,3)   !D(1/R)/D(N) 
                    !SUMHIJ(IP+K1)=SUMHIJ(IP+K1)+FHIJ(I2,2)                           
                END IF

                !END IF
                !IF(IP.EQ.1) WRITE(*,*)SSA(K),SSH(K)
                IF(ISYM.LE.2) THEN
                    !RSH(IP+K1,K)=SSH(K)
                    !RSA(IP+K1,K)=SSA(K)
                END IF
            END DO
        END DO
        END DO
    END DO
!!$OMP END PARALLEL

    IF(STAGE.EQ.2) THEN
        NPLOOP=NBPOINT+NFPOINT
    END IF

    DO I=1,NPLOOP
        DO J=1,NPLOOP
        	!WRITE(18)RSA(I,J),RSH(I,J,1),RSH(I,J,2),RSH(I,J,3),RSH(I,J,4)   
            SA(I,J)=RSA(I,J)
            !SX(I,J)=RSH(I,J,1)
            !SH(I,J)=RSH(I,J,2)
            !SZ(I,J)=RSH(I,J,3)
            !SXX(I,J)=RSH(I,J,4)
        END DO
        !WRITE(118,*)I,SX(I,I),SX(I,I+NFY)
        !111FORMAT(10F15.6)
        !WRITE(*,*)SUMHIJ(I)
    END DO    

    DO I=NBPOINT+(NFXF+NFX-1)*NFY+1,NPLOOP-100,NFY
        DO J=1,4
        	!WRITE(18)RSA(I,J),RSH(I,J,1),RSH(I,J,2),RSH(I,J,3),RSH(I,J,4)   
            !SA(I,J)=RSA(I,J)
            !SX(I,J)=RSH(I,J,1)
            !SH(I,J)=RSH(I,J,2)
            !SZ(I,J)=RSH(I,J,3)
            !SXX(I,J)=RSH(I,J,4)
            !WRITE(118,111)I,I+(J-1)*NFY,SX(I+(J-1)*NFY,I)
        END DO
        !WRITE(118,111)I,SX(I,I),SX(I+NFY,I),SX(I+2.*NFY,I),SX(I+3.*NFY,I),SX(I+4.*NFY,I)
        111FORMAT(2I8,10F15.6)
        !WRITE(*,*)SUMHIJ(I)
    END DO  
    !STOP

!c	释放内存
	DEALLOCATE(PS,CORL,CNL,SSH1,SSH,SSA,SSA1,HIJ,AIJ,RSH,RSA,RFSA,HIJ1,AIJ1,CORD11)
	!WRITE(*,*)'		END COEFMATRIX'

    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
    !PRINT *, 'END   TIME  COEFMATRIX: ', CHAR_TIME2,'  ON  ',CHAR_DATE2

	CLOSE(18)
    CLOSE(118)
    CLOSE(4)
    RETURN
END SUBROUTINE

SUBROUTINE COEFMATRIX_PARA
!   ***************************************************
!   GENERATE THE GLOBAL MATRIX OF 1st INTEGRAL EQUATION  
!   ***************************************************
    
    USE GREENMOD
    USE INPUTDATA

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

    !IMPLICIT NONE
    CHARACTER MARK
	INTEGER:: ISM(3,4)
	DATA ISM/1,1,1,1,-1,1,1,-1,-1,1,1,-1/
	INTEGER:: NLOOP,K,K1,K2,CP,NBLOOP,NPLOOP,IP,I3,I1,J1,IE,J,ISYM,ISINGULAR,I,I2,I11,I22,J11,IAF(2),CAF,IBO(4),IUP1(3),IUP2(3)
    REAL ::FCORL(9,3),FCNL(9,3),FPS(3),T,TL,FAIJ(9),BETA,LOFF,CORL1(9,3),CORL2(9,3),SR(3),D1,D2,D3,D4,XI(6),ETA(6)
    REAL ::HIJX(9),HIJN(9),HIJZ(9),HYX(9),HYN(9),HYZ(9),FHIJ(9,8),VV(3) !诱导系数 
    REAL ::SN(9),SNX(9),SNE(9),AJ,DN(3)
    REAL,ALLOCATABLE:: CORD11(:,:)
    
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2
    CHARACTER*1 creturn
    
    ALLOCATE(CORD11(NTPN+500,3))
	ALLOCATE(SH(NTPN+500,NTPN+500),SX(NTPN+500,NTPN+500),SZ(NTPN+500,NTPN+500),SA(NTPN+500,NTPN+500),SA1(NTPN+500,NTPN+500),SY(NTPN+500,NTPN+500))
	!IF(NIT.EQ.1.AND.ITTE.EQ.1) THEN
        ALLOCATE(SXX(NTPN+500,NTPN+500),SXY(NTPN+500,NTPN+500),SXZ(NTPN+500,NTPN+500),SYY(NTPN+500,NTPN+500),SYZ(NTPN+500,NTPN+500))
    !END IF
    
    
    CORD11=CORD

    !WRITE(*,*)"DFY=     ",DFY,0.6*DFX
    !WRITE(*,*) 
    !WRITE(*,*)'SUB COEFMATRIX..............'

    open(18,file='coefuse.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(18)
    open(19,file='coefuse_T.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(19)
    
    open(118,file='coefuse.DAT')

    NRPOINT=0.

    ALLOCATE (IDS(NTPN))
	ALLOCATE (PS(3),CORL(9,3),CNL(9,3))
	ALLOCATE (SSH1(NTPN+NRPOINT),SSH(NTPN+NRPOINT),SSA(NTPN+NRPOINT),SSA1(NTPN+NRPOINT),HIJ(9),AIJ(9),AIJ1(9),HIJ1(9))
	ALLOCATE (SUMHIJ(NTPN+NRPOINT))
    !WRITE(*,*)"NTPN+NRPOINT",NTPN,NRPOINT
    ALLOCATE (RSH(NTPN+NRPOINT,NTPN+NRPOINT,4),RSA(NTPN+NRPOINT,NTPN+NRPOINT),RFSA(NBPOINT,6))

    OPEN(4,FILE='INDUCED_VELO.DAT')

    RSH=0.;RSA=0.;RFSA=0.
    SX=0.;SY=0.;SZ=0.;SXX=0.;SXY=0.;SXZ=0.;SYY=0.;SYZ=0.;SH=0.

    m=nbpoint
    L=NFPOINT
    NRPOINT=0.
    NRBLOCK=0.
!  NTPN+NRPOINT IN SINGLE CPU
   IF(MOD(NTPN,NTHREAD).NE.0) THEN
     DO I=1,NTHREAD-1
       IF(MOD(NTPN,NTHREAD).EQ.I) THEN
         DO J=1,NTHREAD-1
            NPS1(J)=(NTPN-I)/NTHREAD
         END DO
         NPS1(NTHREAD)=NPS1(NTHREAD-1)+I
         EXIT
       END IF
     END DO
   ELSE 
     DO I=1,NTHREAD
       NPS1(I)=(NTPN)/NTHREAD
     END DO
   END IF

    DO IP=1,NTPN
		SUMHIJ(IP)=0.0
    END DO

    IF(STAGE.EQ.2) THEN
        NPLOOP=NTPN
        NBLOOP=NTP
    END IF

    IF(STAGE.EQ.1) THEN
        NPLOOP=NBPOINT
        NBLOOP=NBBLOCK
    END IF

    creturn = char(10) ! generate carriage return
    
    !K1=0
    CP=1
    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
    !PRINT *, 'END   TIME  DISTRIBUTE: ', CHAR_TIME2,'  ON  ',CHAR_DATE2

            IF(NSWIN(13,1).LE.500) THEN
                KK=4
            ELSE
                KK=2
            END IF
            !KK=2
            
            !write(*,*)kk,nswin(13,1)
    
!START PARALLEL
    CALL OMP_SET_NUM_THREADS(NTHREAD)


!$omp parallel private(CP,k1,k2,FCORL,FCNL,ISYM,FPS,T,HIJ,AIJ,IP,I3,IE,I1,I2,J1,ISINGULAR,TL,K,FAIJ,FHIJ,jm)

    
    CP=OMP_GET_THREAD_NUM()+1

    K1=SUM(NPS1(1:CP))-NPS1(CP) 
    K2=SUM(NPS1(1:CP))
    !K1=0
    !WRITE(*,*)"K1   ",K1
    !write(*,*)NPS1(CP)

    !ALLOCATE (HIJ(9),AIJ(9))

    IDS=0
    !DO IP=1,NPLOOP
    DO IP=1,NPS1(CP)

        !LOFF=LDC*DMIP(IP)**BETA
		!IF((((IP+K1)/(ntpn/4))*(NTPN/4)).EQ.IP+K1)THEN
			!WRITE(*,200)IP+K1,NPLOOP !NTPN
		!END IF
        
        !进度条
       
        IF(CP.EQ.1) THEN
          IF((IP/(NPS1(CP)/20))*(NPS1(CP)/20).EQ.IP)THEN
            creturn = char(13)

          jm=int(30./real(NPS1(CP))*IP*0.8)
          writE(*,'(A\)')creturn
          write( * , '(5a,f6.2,1a)',advance="no")'PROGRESS:' , '[' , & !// 如您的编译器不支持，请用上方语句代替
          repeat('#' , jm ) , repeat( '.' , 30-jm ) , '] ' , IP*80.0/NPS1(CP) , "%"  
          END IF
       END IF 
          
  200	FORMAT('16H    MOVING BALL',3X,'6HINODE=',I8,5X,I8)
   
!		INITIAL THE DIMENSIONS
!		SSH(K)中存放有奇点部分的偶极子效应
!		SSA(K)中存放有奇点部分的源效应
        DO 5 I3=1,NPLOOP
			SSH(I3)=0.
            SSH1(I3)=0.
			SSA(I3)=0.
			SSA1(I3)=0.0
    5		CONTINUE						  
!		将每一个点的坐标赋值到ps（3）中
        
        BETA=0.5
        !LOFF=DMIP(IE)**BETA  !单元内潜上置距离  
        DO 10 J1=1,3
   10		FPS(J1)=CORD(IP+K1,J1) !+VECN(IP,J1)*LOFF
!		LOOP FOR ELEMENTS
	    !FPS(3)=FPS(3)+0.05
    	DO IE=1,NBLOOP !+(NFY-1)/2
!			WRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!			依次将每一面元的PAN(IE)个顶点的坐标赋值到corl(PAN(IE),3)中
!			依次将每一面元的PAN(IE)个顶点的法向矢量赋值到cnl(PAN(IE),3)中

            DO ISYM=1,KK
				
			DO I1=1,PAN(IE)
				I2=I1
                IF(PAN(IE).EQ.3) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=5-I1	
                END IF                
                IF(PAN(IE).EQ.4) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				    IF(MOD(ISYM,4).EQ.0.AND.(I1.NE.1)) I2=6-I1	
                END IF
                IF(PAN(IE).EQ.8) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                END IF
                
                IF(PAN(IE).EQ.9) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                END IF

                DO J1=1,3
                    IF(IE.GT.NBBLOCK) THEN
                        !LOFF=0.8*DFX
                        
                        LOFFZ=1.0*DFX
                        !LOFFZ=1.2*DFX
                        IF(FR.GE.0.3) THEN
                            LOFFZ=1.4*DFX
                        ELSE IF(FR.GE.0.6) THEN
                            LOFFZ=1.2*DFX
                        ELSE IF(FR.GE.0.7) THEN
                            LOFFZ=1.2*DFX
                        ELSE
                            LOFFZ=1.2*DFX
                        END IF
                        LOFFZ=NSWIN(9,1)*DFX
                        !IF(NSWIN(13,1).LE.500) THEN
                        !    LOFFZ=2.0*DFX
                        !END IF
                    ELSE
                        LOFF=DFY
                        !LOFF=1.3*DFY
                    END IF
                    !LOFFZ=1.5*DFX !KCS
                    !WRITE(*,*)DFX !DMIE(IE),BETA, LOFF
                    
                    IF(IE.GT.NBBLOCK.AND.IE.LE.NBBLOCK+NFBLOCK) THEN
                    !IF(ABS(CORD(INE(IE,I2),2)).GE.1.E-6) THEN
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)+VECN(INE(IE,I2),J1)*LOFFZ*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM) 

                        !FCORL(I1,2)=CORD(INE(IE,I2),2)*ISM(J1,ISYM)-LOFFY*ISM(J1,ISYM)

                        if(fr.lT.0.3) FCORL(I1,3)=LOFFZ
                        
                        IF(NSWIN(13,1).GE.500) THEN
                            FCORL(I1,3)=LOFFZ
                        END IF
                        !FCORL(I1,3)=LOFFZ
                        
                        IF(ISYM.GT.2) THEN
                            FCORL(I1,3)=-2.*NSWIN(13,1)+CORD(INE(IE,I2),3)*ISM(3,ISYM)+LOFFZ*ISM(3,ISYM)
                            !FCORL(I1,3)=-2.*NSWIN(13,1)+LOFFZ*ISM(3,ISYM)
                        END IF
                    ELSE
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM)
                        
                        IF(ISYM.GT.2) THEN
                            FCORL(I1,3)=CORD(INE(IE,I2),3)*ISM(3,ISYM)-2.*NSWIN(13,1)
                        END IF
                    END IF
                    
                    !IF(INE(IE,I2).GT.NBPOINT+(NFXF+NFX)*NFY) THEN
                        !FCORL(I1,1)=CORD(INE(IE,I2),1)-LOFFX
                    !END IF    
                END DO

                IF(INE(IE,I2).GT.NBPOINT.AND.INE(IE,I2).LE.NBPOINT+NFPOINT) THEN
                !IF(IP+K1.GT.NBPOINT) THEN
                !IF(INE(IE,I2).GT.NBPOINT+(NFXF+NFX)*NFY) THEN
                    !LOFF=0.8*DFX !yuzheng OK
                    !LOFFX=0.25*DFX
                    IF(MDIFF.EQ.1) THEN
                        LOFFX=1.0*DFX !LINEAR OK
                        LOFFX=0.8*DFX !LINEAR OK
                    ELSE
                        !LOFFX=0.4*DFX !LINEAR OK
                        !LOFFX=0.1*DFX !LINEAR OK
                        !LOFFX=0.3*DFX !LINEAR OK
                        LOFFX=NSWIN(10,1)*DFX
                        !LOFFX=0.4*DFX !LINEAR OK
                        !LOFFX=1.25*DFX !LINEAR OK
                        !write(*,*)LOFFX
                    END IF
                    
                    FCORL(I1,1)=CORD(INE(IE,I2),1)-LOFFX
                END IF

            END DO
            11 continue

!			判断是否是奇点以调用不同的程序，判断标志为ISINGULAR
            !WRITE(*,*)"ps",ps
            !WRITE(*,*)"fcorl",fcorl(1:4,:)
			ISINGULAR=0
			DO I1=1,PAN(IE)
				T=0.
				DO J1=1,3
    				T=T+(FCORL(I1,J1)-FPS(J1))**2
                END DO
				T=SQRT(T)
                IF(T.LT.1.E-5) THEN
					ISINGULAR=1
				END IF
            END DO
            
!c			是奇点则调用cofahs，否则调用cofah,aij为源效应，hij为偶极子效应 
            
			DO I=1,PAN(IE)
				!FHIJ(I,:)=0.0
            END DO
			!AIJ(I)=0.0

			IF(ISINGULAR.EQ.1) THEN
                !IF(PAN(IE).EQ.4) CALL COFAHS_4(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),HIJ(1:PAN(IE)))
                !IF(PAN(IE).EQ.9) CALL COFAHS_9NN(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                !IF(IP+K1.GT.NBPOINT) THEN  
                !    CALL COFAHS_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FAIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                !    CALL COFAHN_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FHIJ(1:PAN(IE),:))
                !ELSE
                    CALL COFAHS_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FAIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                    !FHIJ=0.
                !END IF
                !IF(PAN(IE).EQ.8) THEN
                !    CALL COFAHS_8(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                    !WRITE(*,*)IP,IE,PAN(IE)
                !END IF
                !IF(IP.EQ.1.AND.IE.LE.1) WRITE(*,*)IP,IE,HIJ(1:PAN(IE))
			ELSE IF(ISINGULAR.EQ.0) THEN
                !IF(PAN(IE).EQ.4) CALL COFAH_4(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),HIJ(1:PAN(IE)))
                CALL COFAH_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FAIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                !FHIJ=0.
            ENDIF
            !WRITE(*,*)IP


            IF(IP+K1.EQ.1) THEN
                !CALL COFAH_4SN(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),XGRA,FAIJ)   !与PHIS_T的边界积分方程有关
            END IF

! 			将每一面元中各顶点处的积分值找到对应的总序列中的位置，并进行叠加
! 			SSH(k)中存放偶极子效应
! 			SSA(k)中存放界面上的源效应
!			SUMHIJ(NTPN)用于存放对角线上的量

			DO I1=1,PAN(IE)
                !WRITE(*,*)IE
				K=INE(IE,I1)
				I2=I1
                IF(PAN(IE).EQ.3) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=5-I1	
                END IF                
                IF(PAN(IE).EQ.4) THEN
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				    IF(MOD(ISYM,4).EQ.0.AND.(I1.NE.1)) I2=6-I1	
                END IF
                IF(PAN(IE).EQ.8) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                END IF
                IF(PAN(IE).EQ.9) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                END IF

                !DO KK=1,9
                !    TL=SQRT((CORD(IP+K1,1)-CORD(INE(IE,KK),1))**2+(CORD(IP+K1,2)-CORD(INE(IE,KK),2))**2+(CORD(IP+K1,3)-CORD(INE(IE,KK),3))**2)
                !    IF(TL.LE.1.E-6) THEN
                !        EXIT
                !    END IF
                !END DO

                !IF(K.NE.IP+K1)THEN
                    !IF(TL.GT.1.E-6) THEN
                        !IF(IP+K1.LE.NTPN) THEN
				            !SSH(K)=SSH(K)+HIJ(I2)
                
                !IF(ISYM.LE.2.OR.IE.LE.NBBLOCK) THEN
                !IF(IP.LE.NBPOINT.OR.IE.GT.NBBLOCK) THEN
                IF(K.LE.NTPN) THEN
                !IF(ISYM.LE.2) THEN
                            RSA(IP+K1,K)=RSA(IP+K1,K)+FAIJ(I2)
                
                !IF(K.NE.IP+K1)THEN
                            RSH(IP+K1,K,1)=RSH(IP+K1,K,1)+FHIJ(I2,1)   !D(1/R)/D(X)
                            RSH(IP+K1,K,2)=RSH(IP+K1,K,2)&
                                           +FHIJ(I2,1)*VECN(IP+K1,1)&
                                           +FHIJ(I2,2)*VECN(IP+K1,2)&
                                           +FHIJ(I2,3)*VECN(IP+K1,3)   !D(1/R)/D(N)
                            RSH(IP+K1,K,3)=RSH(IP+K1,K,3)+FHIJ(I2,3)   !D(1/R)/D(Z)
                            RSH(IP+K1,K,4)=RSH(IP+K1,K,4)+FHIJ(I2,4)   !D(1/R)/D(Y)


                            SX(IP+K1,K)=SX(IP+K1,K)+FHIJ(I2,1)
                            SY(IP+K1,K)=SY(IP+K1,K)+FHIJ(I2,2)
                            SZ(IP+K1,K)=SZ(IP+K1,K)+FHIJ(I2,3)
                            SXX(IP+K1,K)=SXX(IP+K1,K)+FHIJ(I2,4)
                            SXY(IP+K1,K)=SXY(IP+K1,K)+FHIJ(I2,5)
                            SXZ(IP+K1,K)=SXZ(IP+K1,K)+FHIJ(I2,6)
                            SYY(IP+K1,K)=SYY(IP+K1,K)+FHIJ(I2,7)
                            SYZ(IP+K1,K)=SYZ(IP+K1,K)+FHIJ(I2,8)
                            SH(IP+K1,K)=SH(IP+K1,K)&
                                           +FHIJ(I2,1)*VECN(IP+K1,1)&
                                           +FHIJ(I2,2)*VECN(IP+K1,2)&
                                           +FHIJ(I2,3)*VECN(IP+K1,3)          

                END IF
                        !END IF
				    !END IF
                !END IFe:\
				    
                !IF(K.EQ.IP+K1)THEN
                IF(K.EQ.IP+K1.AND.ISYM.GE.2)THEN
                !IF(K.NE.IP+K1)THEN
                !IF(TL.GT.1.E-6) THEN
                IF(ABS(CORD(IP+K1,2)).GE.1.E-5) THEN
					SUMHIJ(IP+K1)=SUMHIJ(IP+K1)+FHIJ(I2,1)*VECN(IP+K1,1)&
                                               +FHIJ(I2,2)*VECN(IP+K1,2)&
                                               +FHIJ(I2,3)*VECN(IP+K1,3)   !D(1/R)/D(N) 
                    !SUMHIJ(IP+K1)=SUMHIJ(IP+K1)+FHIJ(I2,2)  
                    SUMHIJ(IP+K1)=0.
                END IF    
                END IF

                !END IF
                !IF(IP.EQ.1) WRITE(*,*)SSA(K),SSH(K)
                IF(ISYM.LE.2) THEN
                    !RSH(IP+K1,K)=SSH(K)
                    !RSA(IP+K1,K)=SSA(K)
                END IF
            END DO
        END DO
        END DO
    END DO
!$OMP END PARALLEL

    IF(STAGE.EQ.2) THEN
        NPLOOP=NBPOINT+NFPOINT
    END IF

    DO I=1,NPLOOP
        DO J=1,NPLOOP
        	!WRITE(18)RSA(I,J),RSH(I,J,1),RSH(I,J,2),RSH(I,J,3),RSH(I,J,4)   
            SA(I,J)=RSA(I,J)
            !SX(I,J)=RSH(I,J,1)
            !SH(I,J)=RSH(I,J,2)
            !SZ(I,J)=RSH(I,J,3)
            !SXX(I,J)=RSH(I,J,4)
        END DO
        !WRITE(118,*)I,SX(I,I),SX(I,I+NFY)
        !111FORMAT(10F15.6)
        !WRITE(*,*)SUMHIJ(I)
    END DO    

    DO I=NBPOINT+(NFXF+NFX-1)*NFY+1,NPLOOP-100,NFY
        DO J=1,4
        	!WRITE(18)RSA(I,J),RSH(I,J,1),RSH(I,J,2),RSH(I,J,3),RSH(I,J,4)   
            !SA(I,J)=RSA(I,J)
            !SX(I,J)=RSH(I,J,1)
            !SH(I,J)=RSH(I,J,2)
            !SZ(I,J)=RSH(I,J,3)
            !SXX(I,J)=RSH(I,J,4)
            !WRITE(118,111)I,I+(J-1)*NFY,SX(I+(J-1)*NFY,I)
        END DO
        !WRITE(118,111)I,SX(I,I),SX(I+NFY,I),SX(I+2.*NFY,I),SX(I+3.*NFY,I),SX(I+4.*NFY,I)
        111FORMAT(2I8,10F15.6)
        !WRITE(*,*)SUMHIJ(I)
    END DO  
    !STOP

    
    GOTO 500
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)     
    PRINT *, 'BEGIN   TIME  MATRIX_PRODUCT: ', CHAR_TIME2,'  ON  ',CHAR_DATE2    
    ABC=0.

    DO I1=1,2
    CALL SAXPY(NTPN,SX(1:NTPN,1:NTPN),SX(1:NTPN,1:NTPN),1,SA(1:NTPN,1:NTPN),1)
    END DO
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)     
    PRINT *, 'END   TIME  MATRIX_PRODUCT_1: ', CHAR_TIME2,'  ON  ',CHAR_DATE2   
    
    DO I1=1,2
    DO I=1,NTPN
    DO J=1,NTPN
    DO K=1,NTPN
       SA(I,J)=SA(I,J)+SX(I,K)*SX(K,J) 
    END DO
    END DO
    END DO    
    END DO    

    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)     
    PRINT *, 'END   TIME  MATRIX_PRODUCT_2: ', CHAR_TIME2,'  ON  ',CHAR_DATE2   
    STOP    
    500 CONTINUE     
    
!c	释放内存
	DEALLOCATE(PS,CORL,CNL,SSH1,SSH,SSA,SSA1,HIJ,AIJ,RSH,RSA,RFSA,HIJ1,AIJ1,CORD11)
	!WRITE(*,*)'		END COEFMATRIX'
    
    
    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
    !PRINT *, 'END   TIME  COEFMATRIX: ', CHAR_TIME2,'  ON  ',CHAR_DATE2

	CLOSE(18)
    CLOSE(118)
    CLOSE(4)
    RETURN
END SUBROUTINE

SUBROUTINE COEFMATRIX_PARA_4
!   ***************************************************
!   GENERATE THE GLOBAL MATRIX OF 1st INTEGRAL EQUATION  
!   ***************************************************
    
    USE GREENMOD

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

    !IMPLICIT NONE
    CHARACTER MARK
	INTEGER:: ISM(3,4)
	DATA ISM/1,1,1,1,-1,1,1,-1,-1,1,1,-1/
	INTEGER:: NLOOP,K,K1,K2,CP,NBLOOP,NPLOOP,IP,I3,I1,J1,IE,J,ISYM,ISINGULAR,I,I2,I11,I22,J11,IAF(2),CAF,IBO(4),IUP1(3),IUP2(3)
    REAL ::FCORL(9,3),FCNL(9,3),FPS(3),T,TL,FAIJ(9),BETA,LOFF,CORL1(9,3),CORL2(9,3),SR(3),D1,D2,D3,D4,XI(6),ETA(6)
    REAL ::HIJX(9),HIJN(9),HIJZ(9),HYX(9),HYN(9),HYZ(9),FHIJ(9,8),VV(3) !诱导系数 
    REAL ::SN(9),SNX(9),SNE(9),AJ,DN(3)
    REAL,ALLOCATABLE:: CORD11(:,:)
    INTEGER,DIMENSION(:,:),ALLOCATABLE:: INE_4(:,:)
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*8 CHAR_DATE1,CHAR_DATE2

    ALLOCATE(CORD11(NTPN+500,3))
	ALLOCATE(SH(NTPN+500,NTPN+500),SX(NTPN+500,NTPN+500),SZ(NTPN+500,NTPN+500),SA(NTPN+500,NTPN+500),SA1(NTPN+500,NTPN+500),SY(NTPN+500,NTPN+500))
	IF(NIT.EQ.1.AND.ITTE.EQ.1) THEN
        ALLOCATE(SXX(NTPN+500,NTPN+500),SXY(NTPN+500,NTPN+500),SXZ(NTPN+500,NTPN+500),SYY(NTPN+500,NTPN+500),SYZ(NTPN+500,NTPN+500))
    END IF

    CORD11=CORD

    WRITE(*,*)"DFY=     ",DFY,0.6*DFX
    WRITE(*,*) 
    WRITE(*,*)'SUB COEFMATRIX..............'

    open(18,file='coefuse.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(18)
    open(19,file='coefuse_T.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(19)
    
    open(118,file='coefuse.DAT')

    NRPOINT=0.

    ALLOCATE (IDS(NTPN))
	ALLOCATE (PS(3),CORL(9,3),CNL(9,3))
	ALLOCATE (SSH1(NTPN+NRPOINT),SSH(NTPN+NRPOINT),SSA(NTPN+NRPOINT),SSA1(NTPN+NRPOINT),HIJ(9),AIJ(9),AIJ1(9),HIJ1(9))
	ALLOCATE (SUMHIJ(NTPN+NRPOINT))
    WRITE(*,*)"NTPN+NRPOINT",NTPN,NRPOINT
    ALLOCATE (RSH(NTPN+NRPOINT,NTPN+NRPOINT,4),RSA(NTPN+NRPOINT,NTPN+NRPOINT),RFSA(NBPOINT,6))

    OPEN(4,FILE='INDUCED_VELO.DAT')

    RSH=0.;RSA=0.;RFSA=0.
    SX=0.;SY=0.;SZ=0.;SXX=0.;SXY=0.;SXZ=0.;SYY=0.;SYZ=0.;SH=0.

    m=nbpoint
    L=NFPOINT
    NRPOINT=0.
    NRBLOCK=0.
!  NTPN+NRPOINT IN SINGLE CPU
   IF(MOD(M+L+NRPOINT,NTHREAD).NE.0) THEN
     DO I=1,NTHREAD-1
       IF(MOD(M+L+NRPOINT,NTHREAD).EQ.I) THEN
         DO J=1,NTHREAD-1
            NPS1(J)=(M+L+NRPOINT-I)/NTHREAD
         END DO
         NPS1(NTHREAD)=NPS1(NTHREAD-1)+I
         EXIT
       END IF
     END DO
   ELSE 
     DO I=1,NTHREAD
       NPS1(I)=(M+L+NRPOINT)/NTHREAD
     END DO
   END IF

    DO IP=1,NTPN
		SUMHIJ(IP)=0.0
    END DO

    IF(STAGE.EQ.2) THEN
        NPLOOP=NTPN
        NBLOOP=NTP
    END IF

    IF(STAGE.EQ.1) THEN
        NPLOOP=NBPOINT
        NBLOOP=NBBLOCK
    END IF

    ALLOCATE(INE_4(NTP*4,4))
    DO IE=1,NTP
        INE_4(IE,1)=INE(IE,1);INE_4(IE,2)=INE(IE,2);INE_4(IE,3)=INE(IE,9);INE_4(IE,4)=INE(IE,8)
        INE_4(IE+NTP,1)=INE(IE,2);INE_4(IE+NTP,2)=INE(IE,3);INE_4(IE+NTP,3)=INE(IE,4);INE_4(IE+NTP,4)=INE(IE,9)
        INE_4(IE+NTP*2,1)=INE(IE,4);INE_4(IE+NTP*2,2)=INE(IE,5);INE_4(IE+NTP*2,3)=INE(IE,6);INE_4(IE+NTP*2,4)=INE(IE,9)
        INE_4(IE+NTP*3,1)=INE(IE,6);INE_4(IE+NTP*3,2)=INE(IE,7);INE_4(IE+NTP*3,3)=INE(IE,8);INE_4(IE+NTP*3,4)=INE(IE,9)
    END DO
    
    !K1=0
    CP=1


!START PARALLEL
    CALL OMP_SET_NUM_THREADS(NTHREAD)


!$omp parallel private(CP,k1,k2,FCORL,FCNL,ISYM,FPS,T,HIJ,AIJ,IP,I3,IE,I1,I2,J1,ISINGULAR,TL,K,FAIJ,FHIJ)

    
    CP=OMP_GET_THREAD_NUM()+1

    K1=SUM(NPS1(1:CP))-NPS1(CP) 
    K2=SUM(NPS1(1:CP))
    !K1=0
    !WRITE(*,*)"K1   ",K1
    !write(*,*)NPS1(CP)

    !ALLOCATE (HIJ(9),AIJ(9))

    IDS=0
    !DO IP=1,NPLOOP
    DO IP=1,NPS1(CP)

        !LOFF=LDC*DMIP(IP)**BETA
		IF((((IP+K1)/(ntpn/4))*(NTPN/4)).EQ.IP+K1)THEN
			WRITE(*,200)IP+K1,NPLOOP !NTPN
		END IF
  200	FORMAT('16H    MOVING BALL',3X,'6HINODE=',I8,5X,I8)
   
!		INITIAL THE DIMENSIONS
!		SSH(K)中存放有奇点部分的偶极子效应
!		SSA(K)中存放有奇点部分的源效应
        DO 5 I3=1,NPLOOP
			SSH(I3)=0.
            SSH1(I3)=0.
			SSA(I3)=0.
			SSA1(I3)=0.0
    5		CONTINUE						  
!		将每一个点的坐标赋值到ps（3）中
        
        BETA=0.5
        !LOFF=DMIP(IE)**BETA  !单元内潜上置距离  
        DO 10 J1=1,3
   10		FPS(J1)=CORD(IP+K1,J1) !+VECN(IP,J1)*LOFF
!		LOOP FOR ELEMENTS
	    !FPS(3)=FPS(3)+0.05
    	DO IE=1,NBLOOP*4
!			WRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!			依次将每一面元的PAN(IE)个顶点的坐标赋值到corl(PAN(IE),3)中
!			依次将每一面元的PAN(IE)个顶点的法向矢量赋值到cnl(PAN(IE),3)中
			
            IF(MNKDB.EQ.2) THEN
                IF(IE.LE.NBBLOCK) THEN
                    KK=4
                ELSE
                    KK=2
                END IF
            ELSE IF(MNKDB.EQ.1) THEN
                KK=2
            END IF

            DO ISYM=1,KK
				
			DO I1=1,4
				I2=I1
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				    IF(MOD(ISYM,4).EQ.0.AND.(I1.NE.1)) I2=6-I1	

                DO J1=1,3
                    !IF(IE.GT.NBBLOCK) THEN
                    IF(INE_4(IE,I2).GT.NBPOINT) THEN
                        !LOFF=0.8*DFX
                        
                        LOFFZ=1.0*DFX
                        !LOFFZ=1.2*DFX
                        IF(FR.GE.0.4) THEN
                            LOFFZ=1.2*DFX
                        ELSE IF(FR.GE.0.6) THEN
                            LOFFZ=1.5*DFX
                        ELSE IF(FR.GE.0.7) THEN
                            LOFFZ=1.5*DFX
                        ELSE
                            LOFFZ=1.0*DFX
                        END IF
                        !LOFFY=0.4*DFY
                        !LOFFZ=0.1*DFX
                    ELSE
                        LOFFZ=DFY
                        !LOFF=1.3*DFY
                    END IF
                    !WRITE(*,*)DFX !DMIE(IE),BETA, LOFF
                    
                    !IF(IE.GT.NBBLOCK) THEN
                    IF(INE_4(IE,I2).GT.NBPOINT) THEN
                    !IF(ABS(CORD(INE(IE,I2),2)).GE.1.E-6) THEN
                        FCORL(I1,J1)=CORD(INE_4(IE,I2),J1)*ISM(J1,ISYM)+VECN(INE_4(IE,I2),J1)*LOFFZ*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE_4(IE,I2),J1)*ISM(J1,ISYM) 

                        !FCORL(I1,2)=CORD(INE(IE,I2),2)*ISM(J1,ISYM)-LOFFY*ISM(J1,ISYM)

                        if(fr.lT.0.3) FCORL(I1,3)=LOFFZ
                        
                        FCORL(I1,3)=LOFFZ
                    ELSE
                        FCORL(I1,J1)=CORD(INE_4(IE,I2),J1)*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE_4(IE,I2),J1)*ISM(J1,ISYM)
                    END IF
                    
                    !IF(INE(IE,I2).GT.NBPOINT+(NFXF+NFX)*NFY) THEN
                        !FCORL(I1,1)=CORD(INE(IE,I2),1)-LOFFX
                    !END IF    
                END DO

                IF(INE_4(IE,I2).GT.NBPOINT) THEN
                !IF(IP+K1.GT.NBPOINT) THEN
                !IF(INE(IE,I2).GT.NBPOINT+(NFXF+NFX)*NFY) THEN
                    !LOFF=0.8*DFX !yuzheng OK
                    !LOFFX=0.25*DFX
                    IF(MDIFF.EQ.1) THEN
                        LOFFX=1.0*DFX !LINEAR OK
                        !LOFFX=1.2*DFX !LINEAR OK
                    ELSE
                        !LOFFX=0.4*DFX !LINEAR OK
                        LOFFX=0.1*DFX !LINEAR OK
                        !LOFFX=0.3*DFX !LINEAR OK
                        !LOFFX=0.4*DFX !LINEAR OK
                        !LOFFX=1.0*DFX !LINEAR OK
                        !write(*,*)LOFFX
                    END IF
                    
                    FCORL(I1,1)=CORD(INE_4(IE,I2),1)-LOFFX
                END IF

            END DO
            11 continue

!			判断是否是奇点以调用不同的程序，判断标志为ISINGULAR
            !WRITE(*,*)"ps",ps
            !WRITE(*,*)"fcorl",fcorl(1:4,:)
			ISINGULAR=0
			DO I1=1,4
				T=0.
				DO J1=1,3
    				T=T+(FCORL(I1,J1)-FPS(J1))**2
                END DO
				T=SQRT(T)
                IF(T.LT.1.E-5) THEN
					ISINGULAR=1
				END IF
            END DO
            
!c			是奇点则调用cofahs，否则调用cofah,aij为源效应，hij为偶极子效应 
            
			!AIJ(I)=0.0

			IF(ISINGULAR.EQ.1) THEN
                !IF(PAN(IE).EQ.4) CALL COFAHS_4(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),HIJ(1:PAN(IE)))
                !IF(PAN(IE).EQ.9) CALL COFAHS_9NN(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                !IF(IP+K1.GT.NBPOINT) THEN  
                !    CALL COFAHS_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FAIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                !    CALL COFAHN_9(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,FHIJ(1:PAN(IE),:))
                !ELSE
                CALL COFAHS_4(FCORL(1:4,:),FCNL(1:4,:),FPS,FAIJ(1:4),FHIJ(1:4,:))
                    !FHIJ=0.
                !END IF
                !IF(PAN(IE).EQ.8) THEN
                !    CALL COFAHS_8(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),FHIJ(1:PAN(IE),:))
                    !WRITE(*,*)IP,IE,PAN(IE)
                !END IF
                !IF(IP.EQ.1.AND.IE.LE.1) WRITE(*,*)IP,IE,HIJ(1:PAN(IE))
			ELSE IF(ISINGULAR.EQ.0) THEN
                !IF(PAN(IE).EQ.4) CALL COFAH_4(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),FPS,AIJ(1:PAN(IE)),HIJ(1:PAN(IE)))
                CALL COFAH_4(FCORL(1:4,:),FCNL(1:4,:),FPS,FAIJ(1:4),FHIJ(1:4,:))
                !FHIJ=0.
            ENDIF
            !WRITE(*,*)IP


            IF(IP+K1.EQ.1) THEN
                !CALL COFAH_4SN(FCORL(1:PAN(IE),:),FCNL(1:PAN(IE),:),XGRA,FAIJ)   !与PHIS_T的边界积分方程有关
            END IF

! 			将每一面元中各顶点处的积分值找到对应的总序列中的位置，并进行叠加
! 			SSH(k)中存放偶极子效应
! 			SSA(k)中存放界面上的源效应
!			SUMHIJ(NTPN)用于存放对角线上的量

			DO I1=1,4
                !WRITE(*,*)IE
				K=INE_4(IE,I1)
				I2=I1      
				    IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				    IF(MOD(ISYM,4).EQ.0.AND.(I1.NE.1)) I2=6-I1	


                !IF(K.NE.IP+K1)THEN
                    !IF(TL.GT.1.E-6) THEN
                        !IF(IP+K1.LE.NTPN) THEN
				            !SSH(K)=SSH(K)+HIJ(I2)
                
                !IF(ISYM.LE.2.OR.IE.LE.NBBLOCK) THEN
                !IF(IP.LE.NBPOINT.OR.IE.GT.NBBLOCK) THEN
                IF(ISYM.LE.2) THEN
                            RSA(IP+K1,K)=RSA(IP+K1,K)+FAIJ(I2)
                
                !IF(K.NE.IP+K1)THEN
                            RSH(IP+K1,K,1)=RSH(IP+K1,K,1)+FHIJ(I2,1)   !D(1/R)/D(X)
                            RSH(IP+K1,K,2)=RSH(IP+K1,K,2)&
                                           +FHIJ(I2,1)*VECN(IP+K1,1)&
                                           +FHIJ(I2,2)*VECN(IP+K1,2)&
                                           +FHIJ(I2,3)*VECN(IP+K1,3)   !D(1/R)/D(N)
                            RSH(IP+K1,K,3)=RSH(IP+K1,K,3)+FHIJ(I2,3)   !D(1/R)/D(Z)
                            RSH(IP+K1,K,4)=RSH(IP+K1,K,4)+FHIJ(I2,4)   !D(1/R)/D(Y)


                            SX(IP+K1,K)=SX(IP+K1,K)+FHIJ(I2,1)
                            SY(IP+K1,K)=SY(IP+K1,K)+FHIJ(I2,2)
                            SZ(IP+K1,K)=SZ(IP+K1,K)+FHIJ(I2,3)
                            SXX(IP+K1,K)=SXX(IP+K1,K)+FHIJ(I2,4)
                            SXY(IP+K1,K)=SXY(IP+K1,K)+FHIJ(I2,5)
                            SXZ(IP+K1,K)=SXZ(IP+K1,K)+FHIJ(I2,6)
                            SYY(IP+K1,K)=SYY(IP+K1,K)+FHIJ(I2,7)
                            SYZ(IP+K1,K)=SYZ(IP+K1,K)+FHIJ(I2,8)
                            SH(IP+K1,K)=SH(IP+K1,K)&
                                           +FHIJ(I2,1)*VECN(IP+K1,1)&
                                           +FHIJ(I2,2)*VECN(IP+K1,2)&
                                           +FHIJ(I2,3)*VECN(IP+K1,3)

                END IF
                        !END IF
				    !END IF
                !END IF
				    
                !IF(K.EQ.IP+K1)THEN
                IF(K.NE.IP+K1)THEN
                !IF(TL.GT.1.E-6) THEN
                !IF(TL.LE.1.E-6) THEN
					SUMHIJ(IP+K1)=SUMHIJ(IP+K1)+FHIJ(I2,1)*VECN(IP+K1,1)&
                                               +FHIJ(I2,2)*VECN(IP+K1,2)&
                                               +FHIJ(I2,3)*VECN(IP+K1,3)   !D(1/R)/D(N) 
                    !SUMHIJ(IP+K1)=SUMHIJ(IP+K1)+FHIJ(I2,2)                           
                END IF

                !END IF
                !IF(IP.EQ.1) WRITE(*,*)SSA(K),SSH(K)
                IF(ISYM.LE.2) THEN
                    !RSH(IP+K1,K)=SSH(K)
                    !RSA(IP+K1,K)=SSA(K)
                END IF
            END DO
        END DO
        END DO
    END DO
!$OMP END PARALLEL

    IF(STAGE.EQ.2) THEN
        NPLOOP=NBPOINT+NFPOINT
    END IF

    DO I=1,NPLOOP
        DO J=1,NPLOOP
        	!WRITE(18)RSA(I,J),RSH(I,J,1),RSH(I,J,2),RSH(I,J,3),RSH(I,J,4)   
            SA(I,J)=RSA(I,J)
            !SX(I,J)=RSH(I,J,1)
            !SH(I,J)=RSH(I,J,2)
            !SZ(I,J)=RSH(I,J,3)
            !SXX(I,J)=RSH(I,J,4)
        END DO
        !WRITE(118,*)I,SX(I,I),SX(I,I+NFY)
        !111FORMAT(10F15.6)
        !WRITE(*,*)SUMHIJ(I)
    END DO    

    DO I=NBPOINT+(NFXF+NFX-1)*NFY+1,NPLOOP-100,NFY
        DO J=1,4
        	!WRITE(18)RSA(I,J),RSH(I,J,1),RSH(I,J,2),RSH(I,J,3),RSH(I,J,4)   
            !SA(I,J)=RSA(I,J)
            !SX(I,J)=RSH(I,J,1)
            !SH(I,J)=RSH(I,J,2)
            !SZ(I,J)=RSH(I,J,3)
            !SXX(I,J)=RSH(I,J,4)
            !WRITE(118,111)I,I+(J-1)*NFY,SX(I+(J-1)*NFY,I)
        END DO
        !WRITE(118,111)I,SX(I,I),SX(I+NFY,I),SX(I+2.*NFY,I),SX(I+3.*NFY,I),SX(I+4.*NFY,I)
        111FORMAT(2I8,10F15.6)
        !WRITE(*,*)SUMHIJ(I)
    END DO  
    !STOP

    

    
!c	释放内存
	DEALLOCATE(PS,CORL,CNL,SSH1,SSH,SSA,SSA1,HIJ,AIJ,RSH,RSA,RFSA,HIJ1,AIJ1,CORD11)
	WRITE(*,*)'		END COEFMATRIX'


	CLOSE(18)
    CLOSE(118)
    CLOSE(4)
    RETURN
END SUBROUTINE

