SUBROUTINE DIFF_COEMATRIX
    USE GREENMOD
	USE CUMOD

    REAL,ALLOCATABLE ::DME(:,:),DMET(:,:),D_XI(:,:),D_ET(:,:)
    REAL:: D1,DA,DB,U30,U31,U32,U33,C30,C31,C32
    INTEGER:: J,J1,J2,J3

    M=NBPOINT
    L=NFPOINT
    NL=M+L

    ALLOCATE(DME(NL,NL),DMET(NL,NL),D_XI(NL,NL),D_ET(NL,NL),DVMX(NL,NL),DVMY(NL,NL))
    ALLOCATE(DCORD(NL,NL),ECORD(NL,NL))

! 初始化 
    DME=0.
    DMET=0.  
    DVMX=0.
    DVMY=0.

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 横向中心差分
    DO IX=1,NFXF+NFXA+NFX 
        DO IY=1,NFY      !奇数行
            IF(IY.EQ.1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                J2=J1+1

                D1=CORD(J1,2)-CORD(J,2)
                DA=CORD(J2,2)-CORD(J1,2)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,4)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                DMET(J,J2)=U32
                DMET(J,J1)=U31
                DMET(J,J)=U30
                DO JP=1,NL
                    D_ET(J,JP)=U32*DVALU(J2,JP)+U31*DVALU(J1,JP)+U30*DVALU(J,JP)
                END DO
            ELSE IF(IY.GE.2.AND.IY.LE.NFY-1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                J_1=J-1
                
                D1=CORD(J,2)-CORD(J_1,2)
                DA=CORD(J1,2)-CORD(J,2)
                CALL DIFF_COE_C3(D1,DA,C30,C31,C32)
                DCORD(J,4)=CORD(J1,2)*C32+CORD(J,2)*C31+CORD(J_1,2)*C30
                DMET(J,J1)=C32
                DMET(J,J)=C31
                DMET(J,J_1)=C30
                DO JP=1,NL
                    D_ET(J,JP)=C32*DVALU(J1,JP)+C31*DVALU(J,JP)+C30*DVALU(J_1,JP)
                END DO
            ELSE IF(IY.EQ.NFY) THEN
                J=M+(IX-1)*NFY+IY
                J_1=J-1
                J_2=J_1-1

                D1=CORD(J_1,2)-CORD(J,2)
                DA=CORD(J_2,2)-CORD(J_1,2)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                !DCORD(J,4)=CORD(J_2,2)*U32+CORD(J_1,2)*U31+CORD(J,2)*U30
                !DMET(J,J_2)=U32
                !DMET(J,J_1)=U31
                !DMET(J,J)=U30
                DCORD(J,4)=(CORD(J,2)-CORD(J_1,2))/(CORD(J,2)-CORD(J_1,2))
                DO JP=1,NL
                    !D_ET(J,JP)=U32*DVALU(J_2,JP)+U31*DVALU(J_1,JP)+U30*DVALU(J,JP)
                    D_ET(J,JP)=(DVALU(J,JP)-DVALU(J_1,JP))/(CORD(J,2)-CORD(J_1,2))
                END DO
            END IF
        END DO
    END DO

 

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分 船尾
    !OPEN(14,FILE='DIFF_OPERA.DAT')
    DO IX=1,NFXF+NFXA+NFX     !沿横向求DX/DE
        IF(IX.EQ.NFXF+NFXA+NFX) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J_1=J-NFY
                J_2=J_1-NFY

                D1=CORD(J_1,1)-CORD(J,1)
                DA=CORD(J_2,1)-CORD(J_1,1)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,1)=CORD(J_2,1)*U32+CORD(J_1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J_2,2)*U32+CORD(J_1,2)*U31+CORD(J,2)*U30
                DME(J,J_2)=U32
                DME(J,J_1)=U31
                DME(J,J)=U30
                DO JP=1,NL
                    D_XI(J,JP)=U32*DVALU(J_2,JP)+U31*DVALU(J_1,JP)+U30*DVALU(J,JP)
                END DO
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-1) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J_1=J-NFY
                J1=J+NFY

                D1=CORD(J,1)-CORD(J_1,1)
                DA=CORD(J1,1)-CORD(J,1)
                CALL DIFF_COE_C3(D1,DA,C30,C31,C32)
                DCORD(J,1)=CORD(J1,1)*C32+CORD(J,1)*C31+CORD(J_1,1)*C30
                DCORD(J,2)=CORD(J1,2)*C32+CORD(J,2)*C31+CORD(J_1,2)*C30
                DME(J,J1)=C32
                DME(J,J)=C31
                DME(J,J_1)=C30
                DO JP=1,NL
                    D_ET(J,JP)=C32*DVALU(J1,JP)+C31*DVALU(J,JP)+C30*DVALU(J_1,JP)
                END DO
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-2) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,1)=CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                DME(J,J2)=U32
                DME(J,J1)=U31
                DME(J,J)=U30
                DO JP=1,NL
                    D_XI(J,JP)=U32*DVALU(J2,JP)+U31*DVALU(J1,JP)+U30*DVALU(J,JP)
                END DO
            END DO
        ELSE IF(IX.LE.NFXF+NFXA+NFX-3) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY
                J3=J2+NFY

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                DB=CORD(J3,1)-CORD(J2,1)     
                !WRITE(*,*)D1,DA,DB
                !CALL DIFF_COE_U4(D1,DA,DB,U30,U31,U32,U33)
                CALL DIFF_COE_D4(J,J1,J2,J3,U30,U31,U32,U33)
                DCORD(J,1)=CORD(J3,1)*U33+CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J3,2)*U33+CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                DME(J,J3)=U33
                DME(J,J2)=U32
                DME(J,J1)=U31
                DME(J,J)=U30
                DO JP=1,NL
                    D_XI(J,JP)=U33*DVALU(J3,JP)+U32*DVALU(J2,JP)+U31*DVALU(J1,JP)+U30*DVALU(J,JP)
                END DO
                !WRITE(14,102)U33*6,U32*6,U31*6,U30*6
                !102format(4f15.6)         
            END DO
        END IF
    END DO
    !CLOSE(14)


    DO I=M+1,NL
        ECORD(I,1)=DCORD(I,4)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))  !EPS_X
        ECORD(I,2)=-DCORD(I,3)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))  !ETA_Y
        ECORD(I,3)=-DCORD(I,2)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2)) 

        ECORD(I,4)=DCORD(I,1)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))
        !WRITE(*,*)DCORD(I,1)
        !其余项为0
    END DO

    DO I=1,NL
        DO J=1,NL
            !DMX(I,J)=DME(I,J)*ECORD(I,1)+DMET(I,J)*ECORD(I,3)
            !DMY(I,J)=DME(I,J)*ECORD(I,2)+DMET(I,J)*ECORD(I,4)
            DVMX(I,J)=D_XI(I,J)*ECORD(I,1)+D_ET(I,J)*ECORD(I,3)
            DVMY(I,J)=D_XI(I,J)*ECORD(I,2)+D_ET(I,J)*ECORD(I,4)
        END DO
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA
! ECORD(I,J)  J=1 DE/DX  J=2 DE/DY  J=3 DETA/DX  J=4 DETA/DY
    DEALLOCATE(DME,DMET,DCORD,ECORD,D_XI,D_ET)
END SUBROUTINE



SUBROUTINE DIFF_COEMATRIX_DECOMP   
    USE GREENMOD
	USE CUMOD

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF    
    
    REAL,ALLOCATABLE ::DME(:,:),DMET(:,:),D_XI(:,:),D_ET(:,:)
    REAL:: D1,DA,DB,U30,U31,U32,U33,C30,C31,C32
    INTEGER:: J,J1,J2,J3

    M=NBPOINT
    L=NFPOINT
    NL=NTPN

    ALLOCATE(DME(NL,NL),DMET(NL,NL),D_XI(NL,NL),D_ET(NL,NL),DVMX(NL,NL),DVMY(NL,NL))
    ALLOCATE(DCORD(NL,NL),ECORD(NL,NL))

! 初始化 
    DME=0.
    DMET=0.  
    DVMX=0.
    DVMY=0.

    
!$omp parallel private(CP,k1,k2,JP,IX,IY,J,J1,J2,j3,D1,DA,J_1,J_2,U30,U31,U32,C30,C31,C32,U33)

    
    CP=OMP_GET_THREAD_NUM()+1

    K1=SUM(NPS1(1:CP))-NPS1(CP) 
    K2=SUM(NPS1(1:CP))
    !K1=0
    !WRITE(*,*)"K1   ",K1
    !write(*,*)cp,NPS1(CP)

    !ALLOCATE (HIJ(9),AIJ(9))

    IDS=0
    !DO IP=1,NPLOOP
    DO  JP=1,NPS1(CP)    
! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 横向中心差分
    !DO JP=1,NL
    DO IX=1,NFXF+NFXA+NFX 
        DO IY=1,NFY      !奇数行
            IF(IY.EQ.1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                J2=J1+1

                !D1=CORD(J1,2)-CORD(J,2)
                !DA=CORD(J2,2)-CORD(J1,2)
                !CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                !DCORD(J,4)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DMET(J,J2)=U32
                !DMET(J,J1)=U31
                !DMET(J,J)=U30
                !DO JP=1,NL
                !    D_ET(J,JP)=U32*DVALU(J2,JP)+U31*DVALU(J1,JP)+U30*DVALU(J,JP)
                !END DO

                D1=CORD(J1,2)-CORD(J,2)
                DA=CORD(J2,2)-CORD(J1,2)
                !CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,4)=(CORD(J1,2)-CORD(J,2))/D1
                !DMET(J,J2)=U32
                !DMET(J,J1)=1./D1
                !DMET(J,J)=-1./D1
                !DO JP=1,NL
                    D_ET(J,K1+JP)=(DVALU(J1,K1+JP)-DVALU(J,K1+JP))/D1
                !END DO

            ELSE IF(IY.GE.2.AND.IY.LE.NFY-1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                J_1=J-1
                
                D1=CORD(J,2)-CORD(J_1,2)
                DA=CORD(J1,2)-CORD(J,2)
                CALL DIFF_COE_C3(D1,DA,C30,C31,C32)
                DCORD(J,4)=CORD(J1,2)*C32+CORD(J,2)*C31+CORD(J_1,2)*C30
                !DMET(J,J1)=C32
                !DMET(J,J)=C31
                !DMET(J,J_1)=C30
                !DO JP=1,NL
                    D_ET(J,K1+JP)=C32*DVALU(J1,K1+JP)+C31*DVALU(J,K1+JP)+C30*DVALU(J_1,K1+JP)
                !END DO
            ELSE IF(IY.EQ.NFY) THEN
                J=M+(IX-1)*NFY+IY
                J_1=J-1
                J_2=J_1-1

                D1=CORD(J_1,2)-CORD(J,2)
                DA=CORD(J_2,2)-CORD(J_1,2)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,4)=CORD(J_2,2)*U32+CORD(J_1,2)*U31+CORD(J,2)*U30
                !DMET(J,J_2)=U32
                !DMET(J,J_1)=U31
                !DMET(J,J)=U30
                !DO JP=1,NL
                    D_ET(J,K1+JP)=U32*DVALU(J_2,K1+JP)+U31*DVALU(J_1,K1+JP)+U30*DVALU(J,K1+JP)
                !END DO
            END IF
        END DO
    END DO

    IF(MTROM.EQ.1) THEN !TRANSOM
    
    DO IX=1,NFXF+1
        DO IY=1,NZS      !奇数行
            IF(IY.EQ.1) THEN
                J=M+(NFXF+NFX+NFXA)*NFY+(IX-1)*NZS+IY
                J1=J+1
                J2=J1+1

                !D1=CORD(J1,2)-CORD(J,2)
                !DA=CORD(J2,2)-CORD(J1,2)
                !CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                !DCORD(J,4)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DMET(J,J2)=U32
                !DMET(J,J1)=U31
                !DMET(J,J)=U30
                !DO JP=1,NL
                    !D_ET(J,JP)=U32*DVALU(J2,JP)+U31*DVALU(J1,JP)+U30*DVALU(J,JP)
                !END DO

                D1=CORD(J1,2)-CORD(J,2)
                DA=CORD(J2,2)-CORD(J1,2)
                !CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,4)=(CORD(J1,2)-CORD(J,2))/D1
                !DMET(J,J2)=U32
                !DMET(J,J1)=1./D1
                !DMET(J,J)=-1./D1
                !DO JP=1,NL
                    D_ET(J,k1+JP)=(DVALU(J1,K1+JP)-DVALU(J,K1+JP))/D1
                !END DO
            ELSE IF(IY.EQ.NZS) THEN
                J=M+(NFXF+NFX+NFXA)*NFY+(IX-1)*NZS+IY
                J1=M+(IX-1)*NFY+1
                J_1=J-1
                
                !D1=CORD(J,2)-CORD(J_1,2)
                !DA=CORD(J1,2)-CORD(J,2)
                !CALL DIFF_COE_C3(D1,DA,C30,C31,C32)
                !DCORD(J,4)=CORD(J1,2)*C32+CORD(J,2)*C31+CORD(J_1,2)*C30
                !DMET(J,J1)=C32
                !DMET(J,J)=C31
                !DMET(J,J_1)=C30
                !DO JP=1,NL
                !    D_ET(J,JP)=C32*DVALU(J1,JP)+C31*DVALU(J,JP)+C30*DVALU(J_1,JP)
                !END DO

                !D1=CORD(J1,2)-CORD(J,2)
                !DA=CORD(J2,2)-CORD(J1,2)
                !CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,4)=(CORD(J,2)-CORD(J_1,2))/D1
                !DMET(J,J2)=U32
                !DMET(J,J)=1./D1
                !DMET(J,J_1)=-1./D1
                !DO JP=1,NL
                    D_ET(J,K1+JP)=(DVALU(J,K1+JP)-DVALU(J_1,K1+JP))/D1
                !END DO
            ELSE 
                J=M+(NFXF+NFX+NFXA)*NFY+(IX-1)*NZS+IY
                J1=J+1
                J_1=J-1
                
                D1=CORD(J,2)-CORD(J_1,2)
                DA=CORD(J1,2)-CORD(J,2)
                CALL DIFF_COE_C3(D1,DA,C30,C31,C32)
                DCORD(J,4)=CORD(J1,2)*C32+CORD(J,2)*C31+CORD(J_1,2)*C30
                !DMET(J,J1)=C32
                !DMET(J,J)=C31
                !DMET(J,J_1)=C30
                !DO JP=1,NL
                    D_ET(J,K1+JP)=C32*DVALU(J1,K1+JP)+C31*DVALU(J,K1+JP)+C30*DVALU(J_1,K1+JP)
                !END DO
            END IF
        END DO
    END DO
    END IF

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分 船尾
    !OPEN(14,FILE='DIFF_OPERA.DAT')
    
    DO IX=1,NFXF     !沿横向求DX/DE
        IF(IX.GE.NFXF+2) THEN 
        !IF(IX.LE.NFXF+2) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY

                !D1=CORD(J1,1)-CORD(J,1)
                !DA=CORD(J2,1)-CORD(J1,1)
                !CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                !DCORD(J,1)=CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                !DCORD(J,2)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DME(J,J2)=U32
                !DME(J,J1)=U31
                !DME(J,J)=U30
                !DO JP=1,NL
                !    D_XI(J,JP)=U32*DVALU(J2,JP)+U31*DVALU(J1,JP)+U30*DVALU(J,JP)
                !END DO

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)

                !DME(J,J1)=1./D1
                !DME(J,J)=-1./D1
                DCORD(J,1)=(CORD(J1,1)-CORD(J,1))/(CORD(J1,1)-CORD(J,1))
                DCORD(J,2)=(CORD(J1,2)-CORD(J,2))/(CORD(J1,1)-CORD(J,1))
                !DO JP=1,NL
                    D_XI(J,K1+JP)=(DVALU(J1,K1+JP)-DVALU(J,K1+JP))/(CORD(J1,1)-CORD(J,1))
                !END DO
            END DO
        ELSE 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY
                J3=J2+NFY

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                DB=CORD(J3,1)-CORD(J2,1)     
                !WRITE(*,*)D1,DA,DB
                !CALL DIFF_COE_U4(D1,DA,DB,U30,U31,U32,U33)
                CALL DIFF_COE_D4(J,J1,J2,J3,U30,U31,U32,U33)
                DCORD(J,1)=CORD(J3,1)*U33+CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J3,2)*U33+CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DME(J,J3)=U33
                !DME(J,J2)=U32
                !DME(J,J1)=U31
                !DME(J,J)=U30
                !DO JP=1,NL
                    D_XI(J,K1+JP)=U33*DVALU(J3,K1+JP)+U32*DVALU(J2,K1+JP)+U31*DVALU(J1,K1+JP)+U30*DVALU(J,K1+JP)
                !END DO
                !WRITE(14,102)U33*6,U32*6,U31*6,U30*6
                !102format(4f15.6)         
            END DO
        END IF         
    END DO

! 纵向迎风差分 船中
    !OPEN(14,FILE='DIFF_OPERA.DAT')
    DO IX=NFXF+1,NFXF+NFX-1     !沿横向求DX/DE
        !IF(IX.EQ.NFXF+NFX-1) THEN 
        IF(IX.EQ.NFXF+NFX) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY


                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)

                !DME(J,J1)=1./D1
                !DME(J,J)=-1./D1
                DCORD(J,1)=(CORD(J1,1)-CORD(J,1))/(CORD(J1,1)-CORD(J,1))
                DCORD(J,2)=(CORD(J1,2)-CORD(J,2))/(CORD(J1,1)-CORD(J,1))
                !DO JP=1,NL
                    D_XI(J,K1+JP)=(DVALU(J1,K1+JP)-DVALU(J,K1+JP))/(CORD(J1,1)-CORD(J,1))
                !END DO
            END DO
        !ELSE IF(IX.EQ.NFXF+NFX-2) THEN 
        ELSE IF(IX.EQ.NFXF+NFX+1) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,1)=CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DME(J,J2)=U32
                !DME(J,J1)=U31
                !DME(J,J)=U30
                !DO JP=1,NL
                    D_XI(J,K1+JP)=U32*DVALU(J2,K1+JP)+U31*DVALU(J1,K1+JP)+U30*DVALU(J,K1+JP)
                !END DO
            END DO
        ELSE 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY
                J3=J2+NFY

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                DB=CORD(J3,1)-CORD(J2,1)     
                !WRITE(*,*)D1,DA,DB
                !CALL DIFF_COE_U4(D1,DA,DB,U30,U31,U32,U33)
                CALL DIFF_COE_D4(J,J1,J2,J3,U30,U31,U32,U33)
                DCORD(J,1)=CORD(J3,1)*U33+CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J3,2)*U33+CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DME(J,J3)=U33
                !DME(J,J2)=U32
                !DME(J,J1)=U31
                !DME(J,J)=U30
                !DO JP=1,NL
                    D_XI(J,K1+JP)=U33*DVALU(J3,K1+JP)+U32*DVALU(J2,K1+JP)+U31*DVALU(J1,K1+JP)+U30*DVALU(J,K1+JP)
                !END DO
                !WRITE(14,102)U33*6,U32*6,U31*6,U30*6
                !102format(4f15.6)         
            END DO
        END IF
    END DO

! 纵向迎风差分 船首
    !OPEN(14,FILE='DIFF_OPERA.DAT')
    DO IX=NFXF+NFX,NFXF+NFXA+NFX     !沿横向求DX/DE
        IF(IX.EQ.NFXF+NFXA+NFX) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J_1=J-NFY
                J2=J1+NFY


                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)

                !DME(J,J)=1./D1
                !DME(J,J_1)=-1./D1
                DCORD(J,1)=(CORD(J,1)-CORD(J_1,1))/(CORD(J,1)-CORD(J_1,1))
                DCORD(J,2)=(CORD(J,2)-CORD(J_1,2))/(CORD(J,1)-CORD(J_1,1))
                !DO JP=1,NL
                    D_XI(J,K1+JP)=(DVALU(J,K1+JP)-DVALU(J_1,K1+JP))/(CORD(J,1)-CORD(J_1,1))
                !END DO
            END DO

        ELSE IF(IX.EQ.NFXF+NFXA+NFX-1) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY


                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)

                !DME(J,J1)=1./D1
                !DME(J,J)=-1./D1
                DCORD(J,1)=(CORD(J1,1)-CORD(J,1))/(CORD(J1,1)-CORD(J,1))
                DCORD(J,2)=(CORD(J1,2)-CORD(J,2))/(CORD(J1,1)-CORD(J,1))
                !DO JP=1,NL
                    D_XI(J,K1+JP)=(DVALU(J1,K1+JP)-DVALU(J,K1+JP))/(CORD(J1,1)-CORD(J,1))
                !END DO
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-2) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,1)=CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DME(J,J2)=U32
                !DME(J,J1)=U31
                !DME(J,J)=U30
                !DO JP=1,NL
                    D_XI(J,K1+JP)=U32*DVALU(J2,K1+JP)+U31*DVALU(J1,K1+JP)+U30*DVALU(J,K1+JP)
                !END DO
            END DO
        ELSE IF(IX.LE.NFXF+NFXA+NFX-3) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY
                J3=J2+NFY

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                DB=CORD(J3,1)-CORD(J2,1)     
                !WRITE(*,*)D1,DA,DB
                !CALL DIFF_COE_U4(D1,DA,DB,U30,U31,U32,U33)
                CALL DIFF_COE_D4(J,J1,J2,J3,U30,U31,U32,U33)
                DCORD(J,1)=CORD(J3,1)*U33+CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J3,2)*U33+CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DME(J,J3)=U33
                !DME(J,J2)=U32
                !DME(J,J1)=U31
                !DME(J,J)=U30
                !DO JP=1,NL
                    D_XI(J,K1+JP)=U33*DVALU(J3,K1+JP)+U32*DVALU(J2,K1+JP)+U31*DVALU(J1,K1+JP)+U30*DVALU(J,K1+JP)
                !END DO
                !WRITE(14,102)U33*6,U32*6,U31*6,U30*6
                !102format(4f15.6)         
            END DO
        END IF
    END DO

    !STOP

    !CLOSE(14)
    !TRANSOM
    IF(MTROM.EQ.1) THEN
    DO IX=1,NFXF+1     !沿横向求DX/DE        IF(IX.EQ.NFXF) THEN 
        IF(IX.EQ.NFXF+1) THEN
        !IF(IX.LE.NFXF+1) THEN
            DO IY=1,NZS    !纵向
                J=M+(NFXF+NFX+NFXA)*NFY+(IX-1)*NZS+IY
                J1=J-NZS
                J2=J1+NZS

                !D1=CORD(J1,1)-CORD(J,1)
                !DA=CORD(J2,1)-CORD(J1,1)
                !CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                !DCORD(J,1)=CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                !DCORD(J,2)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DME(J,J2)=U32
                !DME(J,J1)=U31
                !DME(J,J)=U30
                !DO JP=1,NL
                !    D_XI(J,JP)=U32*DVALU(J2,JP)+U31*DVALU(J1,JP)+U30*DVALU(J,JP)
                !END DO

                !D1=CORD(J1,1)-CORD(J,1)
                !DA=CORD(J2,1)-CORD(J1,1)

                !DME(J,J1)=1./D1
                !DME(J,J)=-1./D1
                DCORD(J,1)=(CORD(J1,1)-CORD(J,1))/(CORD(J1,1)-CORD(J,1))
                DCORD(J,2)=(CORD(J1,2)-CORD(J,2))/(CORD(J1,1)-CORD(J,1))
                !DO JP=1,NL
                    D_XI(J,K1+JP)=(DVALU(J1,K1+JP)-DVALU(J,K1+JP))/(CORD(J1,1)-CORD(J,1))
                !END DO
            END DO
        ELSE IF(IX.EQ.NFXF) THEN 
        !ELSE IF(IX.GT.NFXF+2) THEN 
            DO IY=1,NZS    !纵向
                J=M+(NFXF+NFX+NFXA)*NFY+(IX-1)*NZS+IY
                J1=J+NZS
                J2=J1+NZS

                !D1=CORD(J1,1)-CORD(J,1)
                !DA=CORD(J2,1)-CORD(J1,1)
                !CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                !DCORD(J,1)=CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                !DCORD(J,2)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DME(J,J2)=U32
                !DME(J,J1)=U31
                !DME(J,J)=U30
                !DO JP=1,NL
                !    D_XI(J,JP)=U32*DVALU(J2,JP)+U31*DVALU(J1,JP)+U30*DVALU(J,JP)
                !END DO

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)

                !DME(J,J1)=1./D1
                !DME(J,J)=-1./D1
                DCORD(J,1)=(CORD(J1,1)-CORD(J,1))/(CORD(J1,1)-CORD(J,1))
                DCORD(J,2)=(CORD(J1,2)-CORD(J,2))/(CORD(J1,1)-CORD(J,1))
                !DO JP=1,NL
                    D_XI(J,K1+JP)=(DVALU(J1,K1+JP)-DVALU(J,K1+JP))/(CORD(J1,1)-CORD(J,1))
                !END DO
            END DO
        ELSE IF(IX.EQ.NFXF-1) THEN 
        !ELSE IF(IX.GT.NFXF+1) THEN 
            DO IY=1,NZS    !纵向
                J=M+(NFXF+NFX+NFXA)*NFY+(IX-1)*NZS+IY
                J1=J+NZS
                J2=J1+NZS

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,1)=CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DME(J,J2)=U32
                !DME(J,J1)=U31
                !DME(J,J)=U30
                !DO JP=1,NL
                    D_XI(J,K1+JP)=U32*DVALU(J2,K1+JP)+U31*DVALU(J1,K1+JP)+U30*DVALU(J,K1+JP)
                !END DO
            END DO
        ELSE 
        !ELSE IF(IX.GT.NFXF+1) THEN 
            DO IY=1,NZS    !纵向
                J=M+(NFXF+NFX+NFXA)*NFY+(IX-1)*NZS+IY
                J1=J+NZS
                J2=J1+NZS
                J3=J2+NZS

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                DB=CORD(J3,1)-CORD(J2,1)     
                !WRITE(*,*)D1,DA,DB
                !CALL DIFF_COE_U4(D1,DA,DB,U30,U31,U32,U33)
                CALL DIFF_COE_D4(J,J1,J2,J3,U30,U31,U32,U33)
                DCORD(J,1)=CORD(J3,1)*U33+CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J3,2)*U33+CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                !DME(J,J3)=U33
                !DME(J,J2)=U32
                !DME(J,J1)=U31
                !DME(J,J)=U30
                !DO JP=1,NL
                    D_XI(J,K1+JP)=U33*DVALU(J3,K1+JP)+U32*DVALU(J2,K1+JP)+U31*DVALU(J1,K1+JP)+U30*DVALU(J,K1+JP)
                !END DO
                !WRITE(14,102)U33*6,U32*6,U31*6,U30*6
                !102format(4f15.6)         
            END DO
        END IF         
    END DO
    END IF
    END DO
!$omp end parallel
    
    
99    DO I=M+1,M+L
        ECORD(I,1)=DCORD(I,4)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))  !EPS_X
        ECORD(I,2)=-DCORD(I,3)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))  !ETA_Y
        ECORD(I,3)=-DCORD(I,2)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2)) 
        ECORD(I,4)=DCORD(I,1)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))
        !WRITE(*,*)DCORD(I,1)
        !其余项为0
    END DO

    DO I=1,M+L
        DO J=1,NL
            !DMX(I,J)=DME(I,J)*ECORD(I,1)+DMET(I,J)*ECORD(I,3)
            !DMY(I,J)=DME(I,J)*ECORD(I,2)+DMET(I,J)*ECORD(I,4)
            DVMX(I,J)=D_XI(I,J)*ECORD(I,1)+D_ET(I,J)*ECORD(I,3)
            DVMY(I,J)=D_XI(I,J)*ECORD(I,2)+D_ET(I,J)*ECORD(I,4)
        END DO
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA
! ECORD(I,J)  J=1 DE/DX  J=2 DE/DY  J=3 DETA/DX  J=4 DETA/DY
    DEALLOCATE(DME,DMET,DCORD,ECORD,D_XI,D_ET)
END SUBROUTINE