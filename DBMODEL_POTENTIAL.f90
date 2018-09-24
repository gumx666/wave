SUBROUTINE DBPOTENTIAL
	USE GREENMOD
	USE CUMOD

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

    REAL:: XI,ETA,AJ
    COMMON /GWA/GA(400),GW(400)
	INTEGER:: ISM(3,4)
	DATA ISM/1,1,1,1,-1,1,1,-1,-1,1,1,-1/

    REAL,DIMENSION(:,:):: XXE(2,3)
    REAL,DIMENSION(:):: SN(9),SNX(9),SNE(9),SR(3),DN(3),VV(3),FPS(3),SRR(3)
	REAL,DIMENSION(:,:):: A(3,3),B(3,3)
    REAL,DIMENSION(:,:):: FCORL(9,3),FCNL(9,3)
	DIMENSION IS(3),JS(3)
	INTEGER:: L,ISYM
	!DOUBLE PRECISION A,T,D
	NSHIP=NX*NZ
	ALLOCATE (CORL(9,3),CNL(9,3),PS(3))
	ALLOCATE (SP(NTPN))
    ALLOCATE (DPHIS(NTPN),DPHIS_X(NTPN),DPHIS_Y(NTPN))
    ALLOCATE (DPHIS_XX(NTPN),DPHIS_XY(NTPN),DPHIS_YY(NTPN))

    M=NBPOINT
    L=NFPOINT

! 读入叠模速度势
	OPEN(1,FILE='DB_FLOW.DAT')
	DO I=1,M
        READ(1,*)SP(I)    
    END DO
	CLOSE(1)
	
! 读入诱导系数
    OPEN(16,FILE='COEFUSE_SP.BIN',FORM='BINARY',ACCESS='SEQUENTIAL',&
          STATUS='OLD')
    REWIND(16)
	WRITE(*,*)
	WRITE(*,*)'SUB MATRIX.........'
	ALLOCATE(SH(M,M),SA(M,M))

    THR1=0.

! 读入影响系数
	DO 10 I=1,M
		DO 10 J=1,M
		READ(16)SA(I,J),SH(I,J)
        !CREAD(16,*)SH(I,J),SA(I,J)
10	CONTINUE
    CLOSE(16)

    DPHIS_X=0.
    DPHIS_Y=0.
    DPHIS=0.

    DO I=1,M
        DPHIS(I)=SP(I)
    END DO

!  NTPN+NRPOINT IN SINGLE CPU
    IF(MOD(L,NTHREAD).NE.0) THEN
        DO I=1,NTHREAD-1
            IF(MOD(L,NTHREAD).EQ.I) THEN
                DO J=1,NTHREAD-1
                    NFS(J)=(L-I)/NTHREAD
                END DO
                NFS(NTHREAD)=NFS(NTHREAD-1)+I
                EXIT
            END IF
        END DO
    ELSE 
        DO I=1,NTHREAD
            NFS(I)=(L)/NTHREAD
        END DO
    END IF

	AIJX=0.0

!$omp parallel private(CP,k1,k2,IP,J,J1,K,IE,ISYM,FPS,AIJY,FCORL,FCNL,ISINGULAR,IX,IY,XI,ETA,PHIKSI,RH,RH3,SRR,&
                       BNR,AY,HY,I1,I2,T,SN,SNX,SNE,AJ,SR,DN)
    
    CP=OMP_GET_THREAD_NUM()+1

    K1=SUM(NFS(1:CP))-NFS(CP)+NBPOINT
    K2=SUM(NFS(1:CP))+NBPOINT

    !write(*,*)k2

    DO IP=1,NFS(CP)
    !DO IP=M+1,M+L
        FPS(:)=CORD(IP+K1,:)
        
        IF((((IP+K1-M)/100)*100).EQ.IP+K1-M)THEN
			WRITE(*,200)IP+K1-M,L !NTPN
		END IF
        200	FORMAT('16H    MOVING BALL',3X,'6HINODE=',I8,5X,I8)
        
        DO ISYM=1,4

	    DO IE=1,NBBLOCK
! WRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
! 依次将每一面元的四个顶点的坐标赋值到CORL(4,3)中
! 依次将每一面元的四个顶点的法向矢量赋值到CNL(4,3)中
            DO I1=1,PAN(IE)
                I2=I1
                
                IF(PAN(IE).EQ.9) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                END IF               
                
                IF(PAN(IE).EQ.8) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                END IF

                DO J1=1,3
                    FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)
		            FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM)  
                END DO
            END DO

            ISINGULAR=0
	        DO I1=1,PAN(IE)
				I2=I1
                IF(PAN(IE).EQ.9) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                END IF               
                
                IF(PAN(IE).EQ.8) THEN
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                END IF

                T=0.
		        DO J1=1,3
   				    T=T+(FCORL(I1,J1)-FPS(J1))**2
		        END DO	
            	T=SQRT(T)
				IF(T.LT.1.E-5.AND.IP+K1.NE.INE(IE,I2)) THEN
					ISINGULAR=IP+K1
                    DPHIS(IP+K1)=SP(INE(IE,I2))
                    !WRITE(*,*)IP,PS(1:2)
                    GOTO 100  
				END IF
            END DO

		    NGUS=11
		    K=NGUS*(NGUS-1)/2
!GAUSS-LEGENDRE NUMERICAL INTEGRATION

!对X积分     
		    DO IX=1,NGUS
		        XI=GA(K+IX-1)
!对Y积分
		        DO IY=1,NGUS
			        ETA=GA(K+IY-1)
!调用ISOPAR，求得形状函数J及单位法向矢量等     
			        IF(PAN(IE).EQ.9) CALL ISOPAR_9(XI,ETA,FCORL,SN,SNX,SNE,AJ,SR,DN)
			        !IF(PAN(IE).EQ.8) CALL ISOPAR_9(XI,ETA,FCORL,SN,SNX,SNE,AJ,SR,DN)
			        IF(PAN(IE).EQ.8) CALL ISOPAR_8(XI,ETA,FCORL(1:PAN(IE),:),SN(1:PAN(IE)),SNX(1:PAN(IE)),SNE(1:PAN(IE)),AJ,SR,DN) 

                    PHIKSI=0.
			        DO I1=1,PAN(IE)
			            I2=I1
                        IF(PAN(IE).EQ.9) THEN
                            IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                            IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                        END IF               
                
                        IF(PAN(IE).EQ.8) THEN
                            IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                            IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                            IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                            IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                        END IF
                        PHIKSI=PHIKSI+SN(I1)*SP(INE(IE,I2))
                    END DO

                    RH=0.
                    DO I=1,3
                        SRR(I)=SR(I)-FPS(I)   
                        RH=RH+SRR(I)*SRR(I)
                    END DO
                    RH=SQRT(RH)
	                RH3=-RH**3

                    bnr=srr(1)*dn(1)+srr(2)*dn(2)+srr(3)*dn(3)

                    AIJY=0.
			        DO I1=1,PAN(IE)
                        I2=I1
                        IF(PAN(IE).EQ.9) THEN
                            IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                            IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                        END IF               
                
                        IF(PAN(IE).EQ.8) THEN
                            IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                            IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                            IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                            IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
                        END IF
                        AY=gw(k+iy-1)*SN(I1)*aj/rh           !源
                        !HY=GW(K+IY-1)*SN(I1)*SRR(1)*AJ/RH3   !偶极
                        HY=GW(K+IY-1)*SN(I1)*BNR*AJ/RH3   !偶极
			            
                        AIJY=AIJY+(AY*(U*VECN(INE(IE,I2),1))-HY*SP(INE(IE,I2)))/PI2/2.
                        !AIJY=AIJY+(AY*(U*VECN(INE(IE,I2),1))-HY*SP(INE(IE,I2)))/PI2

			            !AIJY=AIJY+(AY*(U*DN(1))-HY*SP(INE(IE,I2)))/PI2/2.
                    END DO
                    DPHIS(IP+K1)=DPHIS(IP+K1)+GW(k+ix-1)*AIJY
                    
                END DO
                END DO
            END DO
        END DO
        100 CONTINUE
    END DO
!$omp end parallel    

    !对速度势滤波
    ALLOCATE(ETAW(M+L))
    ETAW(M+1:M+L)=DPHIS(M+1:M+L)
    !CALL SMOOTHING_RECT
    DPHIS(M+1:M+L)=ETAW(M+1:M+L)
    DEALLOCATE(ETAW)
    !完成

	WRITE(*,*)'END MATRIX.........'

    DO I=NBPOINT+1,NTPN
        DO J=NBPOINT+1,NTPN
            DPHIS_X(I)=DPHIS_X(I)+DMX(I,J)*DPHIS(J)
            DPHIS_Y(I)=DPHIS_Y(I)+DMY(I,J)*DPHIS(J)
            !DPHIS_X(I)=DPHIS_X(I)+APX(I,J)*DPHIS(J)
            !DPHIS_Y(I)=DPHIS_Y(I)+APY(I,J)*DPHIS(J)
        END DO
    END DO  

    DPHIS_XX=0.
    DPHIS_XY=0.
    DPHIS_YY=0.
    DO I=NBPOINT+1,NTPN
        DO J=NBPOINT+1,NTPN
            DPHIS_XX(I)=DPHIS_XX(I)+DMX(I,J)*DPHIS_X(J)
            DPHIS_XY(I)=DPHIS_XY(I)+DMX(I,J)*DPHIS_Y(J)
            DPHIS_YY(I)=DPHIS_YY(I)+DMY(I,J)*DPHIS_Y(J)
            !DPHIS_X(I)=DPHIS_X(I)+APX(I,J)*DPHIS(J)
            !DPHIS_Y(I)=DPHIS_Y(I)+APY(I,J)*DPHIS(J)
        END DO
    END DO  


    OPEN(5,FILE='BF_SURFACE.DAT')
101	format(1x,'title="episode solid mesh"'/1x,'variables="x","Y",& 
     "Z","NX","NY","NZ"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')             
         
    WRITE(5,101)NFPOINT,NFBLOCK*4 
    DO I=NBPOINT+1,NTPN
        WRITE(5,102)CORD(I,1:2),DPHIS(I),DPHIS_X(I),DPHIS_Y(I),1
    END DO
    102 FORMAT(6F15.6)	
    
    DO I=NBBLOCK+1,NTP
			WRITE(5,*)INE(I,1)-NBPOINT,INE(I,2)-NBPOINT,INE(I,9)-NBPOINT,INE(I,8)-NBPOINT
            WRITE(5,*)INE(I,8)-NBPOINT,INE(I,9)-NBPOINT,INE(I,6)-NBPOINT,INE(I,7)-NBPOINT
            WRITE(5,*)INE(I,2)-NBPOINT,INE(I,3)-NBPOINT,INE(I,4)-NBPOINT,INE(I,9)-NBPOINT
            WRITE(5,*)INE(I,4)-NBPOINT,INE(I,5)-NBPOINT,INE(I,6)-NBPOINT,INE(I,9)-NBPOINT
    END DO

    CLOSE(5)
    CLOSE(16)

	DEALLOCATE(SP)
	DEALLOCATE(CORL,CNL,PS)
    DEALLOCATE(SA,SH)

    RETURN

END SUBROUTINE
