SUBROUTINE DIFFCORD_MATRIX_NUN_1
    USE GREENMOD
	USE CUMOD

    REAL,ALLOCATABLE ::DME(:,:),DMET(:,:)
    REAL:: D1,DA,DB,U30,U31,U32,U33,C30,C31,C32
    INTEGER:: J,J1,J2,J3

    M=NBPOINT
    L=NFPOINT
    NL=M+L

    ALLOCATE(DME(NL,NL),DMET(NL,NL),DMX(NL,NL),DMY(NL,NL))
    ALLOCATE(DCORD(NL,NL),ECORD(NL,NL))

! 初始化 
    DME=0.
    DMET=0.  
    DMX=0.
    DMY=0.

    !DO K=1,NL !NL
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
            ELSE IF(IY.EQ.NFY) THEN
                J=M+(IX-1)*NFY+IY
                J_1=J-1
                J_2=J_1-1

                D1=CORD(J_1,2)-CORD(J,2)
                DA=CORD(J_2,2)-CORD(J_1,2)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,4)=CORD(J_2,2)*U32+CORD(J_1,2)*U31+CORD(J,2)*U30
                DMET(J,J_2)=U32
                DMET(J,J_1)=U31
                DMET(J,J)=U30
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
                !WRITE(14,102)U33*6,U32*6,U31*6,U30*6
                !102format(4f15.6)         
            END DO
        END IF
    END DO
    !CLOSE(14)
    !END DO

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
            DMX(I,J)=DME(I,J)*ECORD(I,1)+DMET(I,J)*ECORD(I,3)
            DMY(I,J)=DME(I,J)*ECORD(I,2)+DMET(I,J)*ECORD(I,4)
        END DO
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA
! ECORD(I,J)  J=1 DE/DX  J=2 DE/DY  J=3 DETA/DX  J=4 DETA/DY
    DEALLOCATE(DME,DMET,DCORD,ECORD)
END SUBROUTINE

SUBROUTINE DIFFCORD_MATRIX_NUN
    USE GREENMOD
	USE CUMOD

    REAL,ALLOCATABLE ::DME(:,:),DMET(:,:)
    REAL:: D1,DA,DB,U30,U31,U32,U33,C30,C31,C32
    INTEGER:: J,J1,J2,J3

    M=NBPOINT
    L=NFPOINT
    NL=M+L

    ALLOCATE(DME(NL,NL),DMET(NL,NL),DMX(NL,NL),DMY(NL,NL))
    ALLOCATE(DCORD(NL,NL),ECORD(NL,NL))

! 初始化 
    DME=0.
    DMET=0.  
    DMX=0.
    DMY=0.

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
            ELSE IF(IY.EQ.NFY) THEN
                J=M+(IX-1)*NFY+IY
                J_1=J-1
                J_2=J_1-1

                D1=CORD(J_1,2)-CORD(J,2)
                DA=CORD(J_2,2)-CORD(J_1,2)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,4)=CORD(J_2,2)*U32+CORD(J_1,2)*U31+CORD(J,2)*U30
                DMET(J,J_2)=U32
                DMET(J,J_1)=U31
                DMET(J,J)=U30
            END IF
        END DO
    END DO

! 横向中心差分船尾网格
    DO IX=1,NFXF+1 
        DO IY=1,NFYS      !奇数行
            IF(IY.EQ.1) THEN
                J=M+NFPOINT0+(IX-1)*NFYS+IY
                J1=J+1
                J2=J1+1
                
                D1=CORD(J1,2)-CORD(J,2)
                DA=CORD(J2,2)-CORD(J1,2)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,4)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                DMET(J,J2)=U32
                DMET(J,J1)=U31
                DMET(J,J)=U30
            ELSE IF(IY.GE.2.AND.IY.LE.NFYS-1) THEN
                J=M+NFPOINT0+(IX-1)*NFYS+IY
                J1=J+1
                J_1=J-1
                D1=CORD(J,2)-CORD(J_1,2)
                DA=CORD(J1,2)-CORD(J,2)
                CALL DIFF_COE_C3(D1,DA,C30,C31,C32)
                DCORD(J,4)=CORD(J1,2)*C32+CORD(J,2)*C31+CORD(J_1,2)*C30
                DMET(J,J1)=C32
                DMET(J,J)=C31
                DMET(J,J_1)=C30
            ELSE IF(IY.EQ.NFYS) THEN
                J=M+NFPOINT0+(IX-1)*NFYS+IY
                J_1=J-1
                J1=M+(IX-1)*NFY+1

                D1=CORD(J,2)-CORD(J_1,2)
                DA=CORD(J1,2)-CORD(J,2)
                CALL DIFF_COE_C3(D1,DA,C30,C31,C32)
                DCORD(J,4)=CORD(J1,2)*C32+CORD(J,2)*C31+CORD(J_1,2)*C30
                DMET(J,J1)=C32
                DMET(J,J)=C31
                DMET(J,J_1)=C30
            END IF
        END DO
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分 船尾
    !OPEN(14,FILE='DIFF_OPERA.DAT')
    DO IX=1,NFXF     !沿横向求DX/DE
        IF(IX.EQ.NFXF) THEN 
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
            END DO
        ELSE IF(IX.EQ.NFXF-1) THEN 
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
            END DO
        ELSE IF(IX.LE.NFXF-2) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY
                J3=J2+NFY

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                DB=CORD(J3,1)-CORD(J2,1)     
                !WRITE(*,*)D1,DA,DB
                CALL DIFF_COE_U4(D1,DA,DB,U30,U31,U32,U33)
                DCORD(J,1)=CORD(J3,1)*U33+CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J3,2)*U33+CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                DME(J,J3)=U33
                DME(J,J2)=U32
                DME(J,J1)=U31
                DME(J,J)=U30  
                !WRITE(14,102)U33*6,U32*6,U31*6,U30*6
                !102format(4f15.6)         
            END DO
        END IF
    END DO
    !CLOSE(14)

! 纵向迎风差分 船中
    !OPEN(14,FILE='DIFF_OPERA.DAT')
    DO IX=NFXF+1,NFXF+NFX-1     !沿横向求DX/DE
        IF(IX.EQ.NFXF+NFX-1) THEN 
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
            END DO
        ELSE IF(IX.EQ.NFXF+NFX-2) THEN 
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
            END DO
        ELSE IF(IX.LE.NFXF+NFX-3) THEN 
            DO IY=1,NFY    !纵向
                J=M+(IX-1)*NFY+IY   !I节点编号
                J1=J+NFY
                J2=J1+NFY
                J3=J2+NFY

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                DB=CORD(J3,1)-CORD(J2,1)     
                !WRITE(*,*)D1,DA,DB
                CALL DIFF_COE_U4(D1,DA,DB,U30,U31,U32,U33)
                DCORD(J,1)=CORD(J3,1)*U33+CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J3,2)*U33+CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                DME(J,J3)=U33
                DME(J,J2)=U32
                DME(J,J1)=U31
                DME(J,J)=U30  
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
                J_2=J_1-NFY

                D1=CORD(J_1,1)-CORD(J,1)
                DA=CORD(J_2,1)-CORD(J_1,1)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,1)=CORD(J_2,1)*U32+CORD(J_1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J_2,2)*U32+CORD(J_1,2)*U31+CORD(J,2)*U30
                DME(J,J_2)=U32
                DME(J,J_1)=U31
                DME(J,J)=U30
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
                CALL DIFF_COE_U4(D1,DA,DB,U30,U31,U32,U33)
                DCORD(J,1)=CORD(J3,1)*U33+CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J3,2)*U33+CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                DME(J,J3)=U33
                DME(J,J2)=U32
                DME(J,J1)=U31
                DME(J,J)=U30  
                !WRITE(14,102)U33*6,U32*6,U31*6,U30*6
                !102format(4f15.6)         
            END DO
        END IF
    END DO

!   船尾部分
    DO IX=1,NFXF+1     !沿横向求DX/DE
        IF(IX.EQ.NFXF+1) THEN 
            DO IY=1,NFYS    !纵向
                J=M+NFPOINT0+(IX-1)*NFYS+IY   !I节点编号
                J_1=J-NFYS
                J_2=J_1-NFYS

                D1=CORD(J_1,1)-CORD(J,1)
                DA=CORD(J_2,1)-CORD(J_1,1)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,1)=CORD(J_2,1)*U32+CORD(J_1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J_2,2)*U32+CORD(J_1,2)*U31+CORD(J,2)*U30
                DME(J,J_2)=U32
                DME(J,J_1)=U31
                DME(J,J)=U30
            END DO
        ELSE IF(IX.EQ.NFXF) THEN 
            DO IY=1,NFYS    !纵向
                J=M+NFPOINT0+(IX-1)*NFYS+IY   !I节点编号
                J1=J+NFYS
                J_1=J-NFYS

                D1=CORD(J,1)-CORD(J_1,1)
                DA=CORD(J1,1)-CORD(J,1)
                CALL DIFF_COE_C3(D1,DA,C30,C31,C32)
                DCORD(J,1)=CORD(J1,1)*C32+CORD(J,1)*C31+CORD(J_1,1)*C30
                DCORD(J,2)=CORD(J1,2)*C32+CORD(J,2)*C31+CORD(J_1,2)*C30
                DME(J,J1)=C32
                DME(J,J)=C31
                DME(J,J_1)=C30
            END DO
        ELSE IF(IX.EQ.NFXF-1) THEN 
            DO IY=1,NFYS    !纵向
                J=M+NFPOINT0+(IX-1)*NFYS+IY   !I节点编号
                J1=J+NFYS
                J2=J1+NFYS

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                CALL DIFF_COE_U3(D1,DA,U30,U31,U32)
                DCORD(J,1)=CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                DME(J,J2)=U32
                DME(J,J1)=U31
                DME(J,J)=U30
            END DO
        ELSE IF(IX.LE.NFXF-2) THEN 
            DO IY=1,NFYS    !纵向
                J=M+NFPOINT0+(IX-1)*NFYS+IY   !I节点编号
                J1=J+NFYS
                J2=J1+NFYS
                J3=J2+NFYS

                D1=CORD(J1,1)-CORD(J,1)
                DA=CORD(J2,1)-CORD(J1,1)
                DB=CORD(J3,1)-CORD(J2,1)     
                !WRITE(*,*)D1,DA,DB
                CALL DIFF_COE_U4(D1,DA,DB,U30,U31,U32,U33)
                DCORD(J,1)=CORD(J3,1)*U33+CORD(J2,1)*U32+CORD(J1,1)*U31+CORD(J,1)*U30
                DCORD(J,2)=CORD(J3,2)*U33+CORD(J2,2)*U32+CORD(J1,2)*U31+CORD(J,2)*U30
                DME(J,J3)=U33
                DME(J,J2)=U32
                DME(J,J1)=U31
                DME(J,J)=U30  
            END DO
        END IF
    END DO

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
            DMX(I,J)=DME(I,J)*ECORD(I,1)+DMET(I,J)*ECORD(I,3)
            DMY(I,J)=DME(I,J)*ECORD(I,2)+DMET(I,J)*ECORD(I,4)
        END DO
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA
! ECORD(I,J)  J=1 DE/DX  J=2 DE/DY  J=3 DETA/DX  J=4 DETA/DY
    DEALLOCATE(DME,DMET,DCORD,ECORD)
END SUBROUTINE



SUBROUTINE DIFF_COE_C3(D1,DA,C30,C31,C32)
    REAL:: D1,DA,C30,C31,C32

    DA=DA/D1
    D1=D1/D1

    C32=1./(DA+DA**2)*(-1.)
    C31=(DA**2-1.)/(DA+DA**2)*(-1.)
    C30=-DA**2/(DA+DA**2)*(-1.)

END SUBROUTINE

SUBROUTINE DIFF_COE_U3(D1,DA,U30,U31,U32)
    REAL:: D1,DA,U30,U31,U32,S

    DA=DA/D1
    D1=D1/D1

    S=1.+DA

    U32=1./(S-S**2)
    U31=-S**2/(S-S**2)
    U30=(S**2-1)/(S-S**2)

END SUBROUTINE

SUBROUTINE DIFF_COE_U4(D1,DA,DB,U30,U31,U32,U33)
    REAL:: D1,DA,DB,U30,U31,U32,U33,S,P

    DA=DA/D1
    DB=DB/D1
    D1=D1/D1

    S=1.+DA
    P=1.+DA+DB
    
    U33=(S**2-S**3)/((S-S**3)*(P**2-P**3)-(P-P**3)*(S**2-S**3))
    U32=-(P**2-P**3)/((S-S**3)*(P**2-P**3)-(P-P**3)*(S**2-S**3))
    U31=(S**3*(P**2-P**3)-P**3*(S**2-S**3))&
        /((S-S**3)*(P**2-P**3)-(P-P**3)*(S**2-S**3))
    U30=-((S**3-1)*(P**2-P**3)-(P**3-1)*(S**2-S**3))&
        /((S-S**3)*(P**2-P**3)-(P-P**3)*(S**2-S**3))

       

END SUBROUTINE


SUBROUTINE DIFF_COE_D4(J,J1,J2,J3,U30,U31,U32,U33)
    USE GREENMOD
    
    REAL:: DD,U30,U31,U32,U33,LJ,LJ1,LJ2,LJ3
    INTEGER:: J,J1,J2,J3

    LJ=0. !CORD(J,1)
    LJ1=1. !CORD(J1,1)
    LJ2=(CORD(J2,1)-CORD(J1,1))/(CORD(J1,1)-CORD(J,1))+LJ1
    LJ3=(CORD(J3,1)-CORD(J2,1))/(CORD(J1,1)-CORD(J,1))+LJ2

    DD=(LJ1-LJ)*(LJ2-LJ)*(LJ3-LJ)*(LJ3-LJ1)*&
       (LJ2-LJ1)*(LJ3-LJ2)*(LJ3+LJ2+LJ1-3.*LJ)

    U33=(LJ1-LJ)**2*(LJ2-LJ)**2*(LJ2-LJ1)*(LJ2+LJ1-2.*LJ)/DD
    U32=-(LJ1-LJ)**2*(LJ3-LJ)**2*(LJ3-LJ1)*(LJ3+LJ1-2.*LJ)/DD
    U31=(LJ2-LJ)**2*(LJ3-LJ)**2*(LJ3-LJ2)*(LJ3+LJ2-2.*LJ)/DD
    U30=-(U33+U32+U31)
    !WRITE(*,*)U30*6.,U31*6.,U32*6.,U33*6.

    !U30=1.71
    !U31=-2.64
    !U32=1.14
    !U33=-0.21
END SUBROUTINE