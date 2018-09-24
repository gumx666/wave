SUBROUTINE DIFFCORD_MATRIX_DECOMP
    USE GREENMOD
	USE CUMOD

    REAL,ALLOCATABLE ::DME(:,:),DMET(:,:)

    M=NBPOINT
    L=NFPOINT
    NL=M+L

    ALLOCATE(DME(NL,NL),DMET(NL,NL),DMX(NL,NL),DMY(NL,NL))

! 初始化 
    DME=0.
    DMET=0.  
    DMX=0.
    DMY=0.

! CALCULATE DY/DETA

    DET=1. !/(NFY-1)

	C31=1./2./DET
	C32=-1./2./DET

	C21=1./DET
	C22=-1./DET

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 横向中心差分
    DO IX=1,NFXF+NFXA+NFX 
        DO IY=1,NFY      !奇数行
            IF(IY.EQ.1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                J2=J1+1
                
                DMET(J,J)=3./2.
                DMET(J,J1)=-4./2.
                DMET(J,J2)=1./2.


                !DMET(J,J1)=C21
                !DMET(J,J)=C22
            ELSE IF(IY.GE.2.AND.IY.LE.NFY-1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                J_1=J-1
                !DCORD(J,4)=CORD(J1,2)*C31+CORD(J_1,1)*C32
                DMET(J,J1)=C31
                DMET(J,J_1)=C32
            ELSE IF(IY.EQ.NFY) THEN
                J=M+(IX-1)*NFY+IY
                J_1=J-1
                !DCORD(J,4)=CORD(J,2)*C21+CORD(J_1,1)*C22
                DMET(J,J)=C21
                DMET(J,J_1)=C22
            END IF
        END DO
    END DO

! CALCULATE DX/DEP
	C41=10.0/6.
	C42=-15.0/6.
	C43=6.0/6.
	C44=-1.0/6.
     	
	!C41=1.71
	!C42=-2.64
	!C43=1.14
	!C44=-0.21

	C31=3.0/2.
	C32=-4.0/2.
	C33=1.0/2.

	C21=1.
	C22=-1.

    WRITE(*,*)'NFXF+NFXA+NFX',NFXF+NFXA+NFX

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分
    DO IX=1,NFXF     !沿横向求DX/DE
        IF(IX.EQ.NFXF) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I_1=I-NFY
                !DCORD(I,1)=CORD(I,1)*C21+CORD(I1,1)*C22
                !DCORD(I,2)=CORD(I,2)*C21+CORD(I1,2)*C22
                DME(I,I)=C21
                DME(I,I_1)=C22
            END DO
        ELSE IF(IX.EQ.NFXF-1) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                !DCORD(I,1)=CORD(I,1)*C31+CORD(I1,1)*C32+CORD(I2,1)*C33
                !DCORD(I,2)=CORD(I,2)*C31+CORD(I1,2)*C32+CORD(I2,2)*C33
                DME(I,I1)=C21
                DME(I,I)=C22
            END DO
        ELSE IF(IX.LE.NFXF-2) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                !DCORD(I,1)=CORD(I,1)*C31+CORD(I1,1)*C32+CORD(I2,1)*C33
                !DCORD(I,2)=CORD(I,2)*C31+CORD(I1,2)*C32+CORD(I2,2)*C33
                DME(I,I)=C31
                DME(I,I1)=C32
                DME(I,I2)=C33
            END DO
        ELSE IF(IX.LE.NFXF-3) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                I3=I2+NFY
                !DCORD(I,1)=CORD(I,1)*C41+CORD(I1,1)*C42+CORD(I2,1)*C43+CORD(I3,1)*C44
                !DCORD(I,2)=CORD(I,2)*C41+CORD(I1,2)*C42+CORD(I2,2)*C43+CORD(I3,2)*C44
                DME(I,I)=C41
                DME(I,I1)=C42
                DME(I,I2)=C43
                DME(I,I3)=C44
            END DO
        END IF
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分
    DO IX=NFXF+1,NFXF+NFX-1     !沿横向求DX/DE
        IF(IX.EQ.NFXF+NFX-1) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                !DCORD(I,1)=CORD(I,1)*C21+CORD(I1,1)*C22
                !DCORD(I,2)=CORD(I,2)*C21+CORD(I1,2)*C22
                DME(I,I)=C21
                DME(I,I1)=C22
            END DO
        ELSE IF(IX.EQ.NFXF+NFX-2) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                !DCORD(I,1)=CORD(I,1)*C31+CORD(I1,1)*C32+CORD(I2,1)*C33
                !DCORD(I,2)=CORD(I,2)*C31+CORD(I1,2)*C32+CORD(I2,2)*C33
                DME(I,I)=C31
                DME(I,I1)=C32
                DME(I,I2)=C33
            END DO
        ELSE IF(IX.LE.NFXF+NFX-3) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                I3=I2+NFY
                !DCORD(I,1)=CORD(I,1)*C41+CORD(I1,1)*C42+CORD(I2,1)*C43+CORD(I3,1)*C44
                !DCORD(I,2)=CORD(I,2)*C41+CORD(I1,2)*C42+CORD(I2,2)*C43+CORD(I3,2)*C44
                DME(I,I)=C41
                DME(I,I1)=C42
                DME(I,I2)=C43
                DME(I,I3)=C44
            END DO
        END IF
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分
    DO IX=NFXF+NFX,NFXF+NFXA+NFX     !沿横向求DX/DE
        IF(IX.EQ.NFXF+NFXA+NFX) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I_1=I-NFY
                !DCORD(I,1)=CORD(I_1,1)*C21+CORD(I,1)*C22   !DF/DX
                !DCORD(I,2)=CORD(I_1,2)*C21+CORD(I,2)*C22   !DF/DX
                DME(I,I_1)=C21
                DME(I,I)=C22
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-1) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                !DCORD(I,1)=CORD(I,1)*C21+CORD(I1,1)*C22
                !DCORD(I,2)=CORD(I,2)*C21+CORD(I1,2)*C22
                DME(I,I)=C21
                DME(I,I1)=C22
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-2) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                !DCORD(I,1)=CORD(I,1)*C31+CORD(I1,1)*C32+CORD(I2,1)*C33
                !DCORD(I,2)=CORD(I,2)*C31+CORD(I1,2)*C32+CORD(I2,2)*C33
                DME(I,I)=C31
                DME(I,I1)=C32
                DME(I,I2)=C33
            END DO
        ELSE IF(IX.LE.NFXF+NFXA+NFX-3) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                I3=I2+NFY
                !DCORD(I,1)=CORD(I,1)*C41+CORD(I1,1)*C42+CORD(I2,1)*C43+CORD(I3,1)*C44
                !DCORD(I,2)=CORD(I,2)*C41+CORD(I1,2)*C42+CORD(I2,2)*C43+CORD(I3,2)*C44
                DME(I,I)=C41
                DME(I,I1)=C42
                DME(I,I2)=C43
                DME(I,I3)=C44
            END DO
        END IF
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

SUBROUTINE DIFFCORD_DECOMP
    USE GREENMOD
	USE CUMOD

    M=NBPOINT
    L=NFPOINT
    NL=M+L
    ALLOCATE(DCORD(NL,NL),ECORD(NL,NL))


! 初始化    
    DCORD=0
    ECORD=0

! CALCULATE DY/DETA

    DET=1. !/(NFY-1)

	C31=1./2./DET
	C32=-1./2./DET

	C21=1./DET
	C22=-1./DET

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 横向中心差分
    DO IX=1,NFXF+NFXA+NFX 
        DO IY=1,NFY      !奇数行
            IF(IY.EQ.1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                J2=J1+1
                !DCORD(J,4)=CORD(J1,2)*C21+CORD(J,2)*C22
                DCORD(J,4)=(CORD(J,2)*3.-CORD(J1,2)*4.+CORD(J2,2))/2.
            ELSE IF(IY.GE.2.AND.IY.LE.NFY-1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                J_1=J-1
                DCORD(J,4)=CORD(J1,2)*C31+CORD(J_1,2)*C32
            ELSE IF(IY.EQ.NFY) THEN
                J=M+(IX-1)*NFY+IY
                J_1=J-1
                DCORD(J,4)=CORD(J,2)*C21+CORD(J_1,2)*C22
            END IF
        END DO
    END DO

! CALCULATE DX/DEP
	C41=10.0/6.
	C42=-15.0/6.
	C43=6.0/6.
	C44=-1.0/6.
     	
	!C41=1.71
	!C42=-2.64
	!C43=1.14
	!C44=-0.21
    	
	C31=3.0/2.
	C32=-4.0/2.
	C33=1.0/2.

	C21=1.
	C22=-1.


! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分
    DO IX=1,NFXF       !沿横向求DX/DE
        IF(IX.EQ.NFXF) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I_1=I-NFY
                DCORD(I,1)=CORD(I,1)*C21+CORD(I_1,1)*C22   !DF/DX
                DCORD(I,2)=CORD(I,2)*C21+CORD(I_1,2)*C22   !DF/DX
            END DO
        ELSE IF(IX.EQ.NFXF-1) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                DCORD(I,1)=CORD(I,1)*C22+CORD(I1,1)*C21   !DF/DX
                DCORD(I,2)=CORD(I,2)*C22+CORD(I1,2)*C21   !DF/DX
            END DO
        ELSE IF(IX.EQ.NFXF-2) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                DCORD(I,1)=CORD(I,1)*C31+CORD(I1,1)*C32+CORD(I2,1)*C33
                DCORD(I,2)=CORD(I,2)*C31+CORD(I1,2)*C32+CORD(I2,2)*C33
            END DO
        ELSE IF(IX.LE.NFXF-3) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                I3=I2+NFY
                DCORD(I,1)=CORD(I,1)*C41+CORD(I1,1)*C42+CORD(I2,1)*C43+CORD(I3,1)*C44
                DCORD(I,2)=CORD(I,2)*C41+CORD(I1,2)*C42+CORD(I2,2)*C43+CORD(I3,2)*C44
            END DO
        END IF
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分
    DO IX=NFXF+1,NFXF+NFX-1    !沿横向求DX/DE
        IF(IX.EQ.NFXF+NFX-1) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                DCORD(I,1)=CORD(I,1)*C21+CORD(I1,1)*C22   !DF/DX
                DCORD(I,2)=CORD(I,2)*C21+CORD(I1,2)*C22   !DF/DX
            END DO
        ELSE IF(IX.EQ.NFXF+NFX-2) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                DCORD(I,1)=CORD(I,1)*C31+CORD(I1,1)*C32+CORD(I2,1)*C33
                DCORD(I,2)=CORD(I,2)*C31+CORD(I1,2)*C32+CORD(I2,2)*C33
            END DO
        ELSE IF(IX.LE.NFXF+NFX-3) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                I3=I2+NFY
                DCORD(I,1)=CORD(I,1)*C41+CORD(I1,1)*C42+CORD(I2,1)*C43+CORD(I3,1)*C44
                DCORD(I,2)=CORD(I,2)*C41+CORD(I1,2)*C42+CORD(I2,2)*C43+CORD(I3,2)*C44
            END DO
        END IF
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分
    DO IX=NFXF+NFX,NFXF+NFXA+NFX    !沿横向求DX/DE
        IF(IX.EQ.NFXF+NFXA+NFX) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I_1=I-NFY
                DCORD(I,1)=CORD(I_1,1)*C21+CORD(I,1)*C22   !DF/DX
                DCORD(I,2)=CORD(I_1,2)*C21+CORD(I,2)*C22   !DF/DX
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-1) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                DCORD(I,1)=CORD(I,1)*C21+CORD(I1,1)*C22
                DCORD(I,2)=CORD(I,2)*C21+CORD(I1,2)*C22
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-2) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                DCORD(I,1)=CORD(I,1)*C31+CORD(I1,1)*C32+CORD(I2,1)*C33
                DCORD(I,2)=CORD(I,2)*C31+CORD(I1,2)*C32+CORD(I2,2)*C33
            END DO
        ELSE IF(IX.LE.NFXF+NFXA+NFX-3) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                I3=I2+NFY
                DCORD(I,1)=CORD(I,1)*C41+CORD(I1,1)*C42+CORD(I2,1)*C43+CORD(I3,1)*C44
                DCORD(I,2)=CORD(I,2)*C41+CORD(I1,2)*C42+CORD(I2,2)*C43+CORD(I3,2)*C44
            END DO
        END IF
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA
! ECORD(I,J)  J=1 DE/DX  J=2 DE/DY  J=3 DETA/DX  J=4 DETA/DY

    DO I=M+1,NL
        ECORD(I,1)=DCORD(I,4)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))  !EPS_X
        ECORD(I,2)=-DCORD(I,3)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))  !ETA_Y
        ECORD(I,3)=-DCORD(I,2)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2)) 
        ECORD(I,4)=DCORD(I,1)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))
        !其余项为0
    END DO
END SUBROUTINE

SUBROUTINE DIFFCORD_MATRIX
    USE GREENMOD
	USE CUMOD

    REAL,ALLOCATABLE ::DME(:,:),DMET(:,:)

    M=NBPOINT
    L=NFPOINT
    NL=M+L

    ALLOCATE(DME(NL,NL),DMET(NL,NL),DMX(NL,NL),DMY(NL,NL))

! 初始化 
    DME=0.
    DMET=0.  
    DMX=0.
    DMY=0.

! CALCULATE DY/DETA

    DET=1. !/(NFY-1)

	C31=1./2./DET
	C32=-1./2./DET

	C21=1./DET
	C22=-1./DET

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 横向中心差分
    DO IX=1,NFXF+NFXA+NFX 
        DO IY=1,NFY      !奇数行
            IF(IY.EQ.1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                !DCORD(J,4)=CORD(J1,2)*C21+CORD(J,1)*C22
                DMET(J,J1)=C21
                DMET(J,J)=C22
            ELSE IF(IY.GE.2.AND.IY.LE.NFY-1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                J_1=J-1
                !DCORD(J,4)=CORD(J1,2)*C31+CORD(J_1,1)*C32
                DMET(J,J1)=C31
                DMET(J,J_1)=C32
            ELSE IF(IY.EQ.NFY) THEN
                J=M+(IX-1)*NFY+IY
                J_1=J-1
                !DCORD(J,4)=CORD(J,2)*C21+CORD(J_1,1)*C22
                DMET(J,J)=C21
                DMET(J,J_1)=C22
            END IF
        END DO
    END DO

! CALCULATE DX/DEP
	C41=10.0/6.
	C42=-15.0/6.
	C43=6.0/6.
	C44=-1.0/6.
     	
	C41=1.71
	C42=-2.64
	C43=1.14
	C44=-0.21

	C31=3.0/2.
	C32=-4.0/2.
	C33=1.0/2.

	C21=1.
	C22=-1.

    write(*,*)'NFXF+NFXA+NFX',NFXF+NFXA+NFX

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分
    DO IX=1,NFXF+NFXA+NFX     !沿横向求DX/DE
        IF(IX.EQ.NFXF+NFXA+NFX) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I_1=I-NFY
                !DCORD(I,1)=CORD(I_1,1)*C21+CORD(I,1)*C22   !DF/DX
                !DCORD(I,2)=CORD(I_1,2)*C21+CORD(I,2)*C22   !DF/DX
                DME(I,I_1)=C21
                DME(I,I)=C22
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-1) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                !DCORD(I,1)=CORD(I,1)*C21+CORD(I1,1)*C22
                !DCORD(I,2)=CORD(I,2)*C21+CORD(I1,2)*C22
                DME(I,I)=C21
                DME(I,I1)=C22
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-2) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                !DCORD(I,1)=CORD(I,1)*C31+CORD(I1,1)*C32+CORD(I2,1)*C33
                !DCORD(I,2)=CORD(I,2)*C31+CORD(I1,2)*C32+CORD(I2,2)*C33
                DME(I,I)=C31
                DME(I,I1)=C32
                DME(I,I2)=C33
            END DO
        ELSE IF(IX.LE.NFXF+NFXA+NFX-3) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                I3=I2+NFY
                !DCORD(I,1)=CORD(I,1)*C41+CORD(I1,1)*C42+CORD(I2,1)*C43+CORD(I3,1)*C44
                !DCORD(I,2)=CORD(I,2)*C41+CORD(I1,2)*C42+CORD(I2,2)*C43+CORD(I3,2)*C44
                DME(I,I)=C41
                DME(I,I1)=C42
                DME(I,I2)=C43
                DME(I,I3)=C44
            END DO
        END IF
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

SUBROUTINE DIFFCORD
    USE GREENMOD
	USE CUMOD

    M=NBPOINT
    L=NFPOINT
    NL=M+L
    ALLOCATE(DCORD(NL,NL),ECORD(NL,NL))


! 初始化    
    DCORD=0
    ECORD=0

! CALCULATE DY/DETA

    DET=1. !/(NFY-1)

	C31=1./2./DET
	C32=-1./2./DET

	C21=1./DET
	C22=-1./DET

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 横向中心差分
    DO IX=1,NFXF+NFXA+NFX 
        DO IY=1,NFY      !奇数行
            IF(IY.EQ.1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                DCORD(J,4)=CORD(J1,2)*C21+CORD(J,2)*C22
            ELSE IF(IY.GE.2.AND.IY.LE.NFY-1) THEN
                J=M+(IX-1)*NFY+IY
                J1=J+1
                J_1=J-1
                DCORD(J,4)=CORD(J1,2)*C31+CORD(J_1,2)*C32
            ELSE IF(IY.EQ.NFY) THEN
                J=M+(IX-1)*NFY+IY
                J_1=J-1
                DCORD(J,4)=CORD(J,2)*C21+CORD(J_1,2)*C22
            END IF
        END DO
    END DO

! CALCULATE DX/DEP
	C41=10.0/6.
	C42=-15.0/6.
	C43=6.0/6.
	C44=-1.0/6.
     	
	C41=1.71
	C42=-2.64
	C43=1.14
	C44=-0.21
    	
	C31=3.0/2.
	C32=-4.0/2.
	C33=1.0/2.

	C21=1.
	C22=-1.


! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA   
! 纵向迎风差分
    DO IX=1,NFXF+NFXA+NFX    !沿横向求DX/DE
        IF(IX.EQ.NFXF+NFXA+NFX) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I_1=I-NFY
                DCORD(I,1)=CORD(I_1,1)*C21+CORD(I,1)*C22   !DF/DX
                DCORD(I,2)=CORD(I_1,2)*C21+CORD(I,2)*C22   !DF/DX
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-1) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                DCORD(I,1)=CORD(I,1)*C21+CORD(I1,1)*C22
                DCORD(I,2)=CORD(I,2)*C21+CORD(I1,2)*C22
            END DO
        ELSE IF(IX.EQ.NFXF+NFXA+NFX-2) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                DCORD(I,1)=CORD(I,1)*C31+CORD(I1,1)*C32+CORD(I2,1)*C33
                DCORD(I,2)=CORD(I,2)*C31+CORD(I1,2)*C32+CORD(I2,2)*C33
            END DO
        ELSE IF(IX.LE.NFXF+NFXA+NFX-3) THEN 
            DO IY=1,NFY    !纵向
                I=M+(IX-1)*NFY+IY   !I节点编号
                I1=I+NFY
                I2=I1+NFY
                I3=I2+NFY
                DCORD(I,1)=CORD(I,1)*C41+CORD(I1,1)*C42+CORD(I2,1)*C43+CORD(I3,1)*C44
                DCORD(I,2)=CORD(I,2)*C41+CORD(I1,2)*C42+CORD(I2,2)*C43+CORD(I3,2)*C44
            END DO
        END IF
    END DO

! DCORD(I,J)  J=1 DX/DE  J=2 DY/DE  J=3 DX/DETA  J=4 DY/DETA
! ECORD(I,J)  J=1 DE/DX  J=2 DE/DY  J=3 DETA/DX  J=4 DETA/DY

    DO I=M+1,NL
        ECORD(I,1)=DCORD(I,4)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))  !EPS_X
        ECORD(I,2)=-DCORD(I,3)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))  !ETA_Y
        ECORD(I,3)=-DCORD(I,2)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2)) 
        ECORD(I,4)=DCORD(I,1)/(DCORD(I,1)*DCORD(I,4)-DCORD(I,3)*DCORD(I,2))
        !其余项为0
    END DO
END SUBROUTINE
