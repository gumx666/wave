SUBROUTINE SOLID_ANGLES
    USE GREENMOD
    USE CUMOD

    INTEGER:: I1,I2,IE,IP,K(30),L(30),N,DN(15000),MNP(15000,10),DNODE1(15000)
    REAL:: ALP(30),TA1(30,3),TA2(30,3),NJ,NJT,TINT(3)
    REAL,ALLOCATABLE:: CORDA(:,:),VEA(:,:,:),TCORD(:,:)
    integer,ALLOCATABLE:: INEA(:,:),PANK(:)
    integer:: INP(15000,30,2)
    INTEGER:: ISM(3,4)
	DATA ISM/1,1,1,1,-1,1,1,-1,-1,1,1,-1/
    OPEN(11,FILE='SANGLES.DAT')
    
    GOTO 123

    ALLOCATE(TCORD(NTPN,3))
    DO I=1,NTPN
        TCORD(I,:)=CORD(I,:)
        CORD(I,:)=PCORD(I,:)
    END DO



    DNODE=0
    ISYM=2
    NTPN=NTPN*ISYM
    NTP=NTP*ISYM
    ALLOCATE(CORDA(NTPN,3),INEA(NTP,9),VEA(NTP,9,3),PANK(NTP))

    DO I=1,NTPN/ISYM
        DO J=1,ISYM
            DO J1=1,3
                CORDA(NTPN/ISYM*(J-1)+I,J1)=CORD(I,J1)*ISM(J1,J) 
            END DO
        END DO
    END DO

    DO J=1,ISYM	
        DO I=1,NTP/ISYM
            PANK(NTP/ISYM*(J-1)+I)=PAN(I)
            DO J1=1,PAN(I)
                I2=J1
                IF(PAN(I).EQ.4) THEN
				    IF(MOD(J,2).EQ.0.AND.(J1.NE.1)) I2=6-J1	
		            !IF(MOD(J,4).EQ.0.AND.(J1.NE.1)) I2=6-J1
                ELSE IF(PAN(I).EQ.3) THEN
				    IF(MOD(J,2).EQ.0.AND.(J1.NE.1)) I2=5-J1
                ELSE IF(PAN(I).EQ.9) THEN
                    IF(MOD(J,2).EQ.2.AND.(J1.NE.1).AND.(J1.NE.9)) I2=10-J1   !保证镜像后的面元节点按逆时针排列
                END IF
                INEA(NTP/ISYM*(J-1)+I,J1)=INE(I,I2)+NTPN/ISYM*(J-1) 

                VEA(NTP/ISYM*(J-1)+I,J1,1)=VE(I,I2,1)*ISM(1,J) 
                VEA(NTP/ISYM*(J-1)+I,J1,2)=VE(I,I2,2)*ISM(2,J) 
                VEA(NTP/ISYM*(J-1)+I,J1,3)=VE(I,I2,3)*ISM(3,J)
            END DO
            !WRITE(*,*)NTP/ISYM*(J-1)+I 
        END DO
    END DO

    DN=0
    DNODE=0
    DO I1=1,NTPN/ISYM
        DO I2=1,NTP !/ISYM
            DO I3=1,PANK(I2)
                IF(INEA(I2,I3).EQ.I1) THEN
                    DNODE(I1)=DNODE(I1)+1
                    INP(I1,INT(DNODE(I1)),1)=I2
                    INP(I1,INT(DNODE(I1)),2)=I3
                    GOTO 90
                END IF
            END DO
        END DO

        90 DO I2=1,NTP !/ISYM
            DO I3=1,PANK(I2)
            IF(SQRT((CORDA(I1,1)-CORDA(INEA(I2,I3),1))**2+(CORDA(I1,2)-CORDA(INEA(I2,I3),2))**2+(CORDA(I1,3)-CORDA(INEA(I2,I3),3))**2).LE.1.E-5&
               .AND.I2.NE.INP(I1,1,1)) THEN
                DNODE(I1)=DNODE(I1)+1
                !WRITE(*,*)DNODE(I1)
                INP(I1,INT(DNODE(I1)),1)=I2
                INP(I1,INT(DNODE(I1)),2)=I3
                !DN(I1)=DN(I1)+1 
                !MNP(I1,DN(I1))=I2    !与I1重合点的编号为MNP(I1,DN(I1))
                !WRITE(*,*)I1,MNP(I1,DN(I1))     
            END IF
            END DO
            !if(INP(I1,INT(DNODE(I1)),1).eq.0) WRITE(*,*)INP(I1,INT(DNODE(I1)),1),INP(I1,INT(DNODE(I1)),2)
        END DO         
    END DO

    !DNODE=0
    !DO IE=1,NTP !/ISYM
    !    DO J=1,4
            !IP=INEA(IE,J)
            !DNODE(IP)=DNODE(IP)+1
            !INP(IP,DNODE(IP),1)=IE            !判断公用的单元
            !INP(IP,DNODE(IP),2)=J             !判断公用点在单元内的节点编号
    !    END DO
    !END DO

    !DNODE1(1:NTPN)=DNODE(1:NTPN)

    !DO IP=1,NTPN !/ISYM
        !IF(DN(IP).GE.1) THEN
            !DO I=1,DN(IP)
                    !DO J=1,DNODE1(MNP(IP,I))
                        !DNODE(IP)=DNODE(IP)+1
                        !WRITE(*,*)DNODE(IP),DNODE1(MNP(IP,I)),MNP(IP,I)
                        !INP(IP,DNODE(IP),1)=INP(MNP(IP,I),J,1)
                        !INP(IP,DNODE(IP),2)=INP(MNP(IP,I),J,2)
                        !WRITE(*,*)INP(IP,DNODE(IP),2),INP(IP,DNODE(IP),1)
                    !END DO
            !END DO
        !END IF
    !END DO

    !DO IP=1,NTPN/ISYM
    DO IP=1,NBPOINT !NTPN/ISYM-NFPOINT
        N=0
        !WRITE(*,*)DNODE(IP)
        IF(ABS(DNODE(IP)-1.).LE.1.E-4) THEN
            SE(IP)=PI
        ELSE IF(ABS(DNODE(IP)-2.).LE.1.E-4) THEN     !两个公用单元情况
            K(1)=INP(IP,1,1)
            K(2)=INP(IP,2,1)
            L(1)=INP(IP,1,2)
            L(2)=INP(IP,2,2)

            DO IE=1,DNODE(IP)  
                IF(PANK(INP(IP,IE,1)).EQ.9) THEN   !9节点单元               
                
                IF(INP(IP,IE,2).EQ.2) THEN      !判断公共节点的两边
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),1),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),3),J)
                    END DO
                END IF
                !CASE(3)
                IF(INP(IP,IE,2).EQ.4) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),3),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),5),J)
                    END DO
                END IF
                IF(INP(IP,IE,2).EQ.6) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),5),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),7),J)
                    END DO
                END IF
                IF(INP(IP,IE,2).EQ.8) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),1),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),7),J)
                    END DO
                END IF

                TA1(IE,:)=CORDA(INEA(INP(IP,IE,1),9),:)
                END IF
            END DO

            N=0
            ALP=0

            I1=1                    !第I1个公用单元
            I2=2
            K(1)=INP(IP,1,1)       !单元编号
            K(2)=INP(IP,2,1)
            L(1)=INP(IP,1,2)       !共用点在单元内的编号
            L(2)=INP(IP,2,2)

            !write(*,*)k(1),l(1)

            NJT=((TA1(I1,1)-CORDA(IP,1))*VEA(K(2),L(2),1)+&   !单元2法向在单元1切线上的投影
                (TA1(I1,2)-CORDA(IP,2))*VEA(K(2),L(2),2)+&
                (TA1(I1,3)-CORDA(IP,3))*VEA(K(2),L(2),3))/&
                SQRT((TA1(I1,1)-CORDA(IP,1))**2+(TA1(I1,2)-CORDA(IP,2))**2+(TA1(I1,3)-CORDA(IP,3))**2)
                
            NJ=(TA1(I1,1)-CORDA(IP,1))*VEA(K(1),L(1),1)+&    !单元1法向在单元1切线上的投影
               (TA1(I1,2)-CORDA(IP,2))*VEA(K(1),L(1),2)+&
               (TA1(I1,3)-CORDA(IP,3))*VEA(K(1),L(1),3)/&
                SQRT((TA1(I1,1)-CORDA(IP,1))**2+(TA1(I1,2)-CORDA(IP,2))**2+(TA1(I1,3)-CORDA(IP,3))**2)
                
                    
                        
                        N=1 !N+1
                        !WRITE(*,*)'NJT',NJT

                        IF(NJ.GT.NJT) THEN
                            SGN=1.
                        ELSE IF(NJ.LT.NJT) THEN
                            SGN=-1.
                        ELSE
                            SGN=0.
                        END IF
                    
                        NJ=VEA(K(1),L(1),1)*VEA(K(2),L(2),1)+VEA(K(1),L(1),2)*VEA(K(2),L(2),2)+VEA(K(1),L(1),3)*VEA(K(2),L(2),3)
                        
                        IF(ABS(NJ-1.0).LE.1.E-5) THEN
                            ALP(N)=PI !+SGN*ACOS(NJ)
                        ELSE IF(ABS(NJ+1.0).LE.1.E-5) THEN
                            ALP(N)=0. !2.*PI
                        ELSE
                            ALP(N)=PI+SGN*ACOS(NJ)
                            !WRITE(*,*)SGN,ACOS(NJ),NJT
                        END IF

                        GOTO 70


            NJ=0.0
            DO J=1,3
                NJ=NJ+VEA(K(1),L(1),J)*VEA(K(2),L(2),J)
            END DO 
            IF (ABS(NJ-1.).LE.1.E-4) THEN
                ALP(1)=PI
            ELSE
                ALP(1)=PI+ACOS(NJ) !+SGN*ACOS(NJ)
            END IF
            SE(IP)=ALP(1) !+PI
            N=1

            SE(IP)=ALP(1)/2.
            70 SE(IP)=SUM(ALP(1:N))-(N-2)*PI
            SE(IP)=PI2-ALP(1) !SE(IP)/2.
        ELSE
            DO IE=1,DNODE(IP)  
                IF(PANK(INP(IP,IE,1)).EQ.4) THEN   !4节点单元               
                
                IF(INP(IP,IE,2).EQ.1) THEN      !判断公共节点的两边
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),2),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),4),J)
                    END DO
                END IF
                !CASE(3)
                IF(INP(IP,IE,2).EQ.3) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),4),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),2),J)
                    END DO
                END IF
                IF(INP(IP,IE,2).EQ.2) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),3),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),1),J)
                    END DO
                END IF
                IF(INP(IP,IE,2).EQ.4) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),1),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),3),J)
                    END DO
                END IF

                ELSE IF(PANK(INP(IP,IE,1)).EQ.3) THEN     !3节点单元  

                IF(INP(IP,IE,2).EQ.1) THEN      !判断公共节点的两边
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),2),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),3),J)
                    END DO
                END IF
                !CASE(3)
                IF(INP(IP,IE,2).EQ.2) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),1),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),3),J)
                    END DO
                END IF
                IF(INP(IP,IE,2).EQ.3) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),1),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),2),J)
                    END DO
                END IF

                ELSE IF(PANK(INP(IP,IE,1)).EQ.9) THEN   !9节点单元               
                
                IF(INP(IP,IE,2).EQ.1) THEN      !判断公共节点的两边
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),2),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),8),J)
                    END DO
                END IF
                !CASE(3)
                IF(INP(IP,IE,2).EQ.3) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),4),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),2),J)
                    END DO
                END IF
                IF(INP(IP,IE,2).EQ.5) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),4),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),6),J)
                    END DO
                END IF
                IF(INP(IP,IE,2).EQ.7) THEN
                    DO J=1,3
                        TA1(IE,J)=CORDA(INEA(INP(IP,IE,1),6),J)
                        TA2(IE,J)=CORDA(INEA(INP(IP,IE,1),8),J)
                    END DO
                END IF
                END IF
            END DO

            N=0
            ALP=0

            I1=1                    !第I1个公用单元
            K(1)=INP(IP,I1,1)       !单元编号
            L(1)=INP(IP,I1,2)       !共用点在单元内的编号

            !write(*,*)k(1),l(1)

            NJT=(((CORDA(IP,2)-TA1(I1,2))*(CORDA(IP,3)-TA2(I1,3))-(CORDA(IP,3)-TA1(I1,3))*(CORDA(IP,2)-TA2(I1,2)))*VEA(K(1),L(1),1)+&
                ((CORDA(IP,3)-TA1(I1,3))*(CORDA(IP,1)-TA2(I1,1))-(CORDA(IP,1)-TA1(I1,1))*(CORDA(IP,3)-TA2(I1,3)))*VEA(K(1),L(1),2)+&
                ((CORDA(IP,1)-TA1(I1,1))*(CORDA(IP,2)-TA2(I1,2))-(CORDA(IP,2)-TA1(I1,2))*(CORDA(IP,1)-TA2(I1,1)))*VEA(K(1),L(1),3)) !/&
                !SQRT(VE(K(1),L(1),1)**2+VE(K(1),L(1),2)**2+VE(K(1),L(1),3)**2)/&
                !SQRT(((CORDA(IP,2)-TA1(I1,2))*(CORDA(IP,3)-TA2(I1,3))-(CORDA(IP,3)-TA1(I1,3))*(CORDA(IP,2)-TA2(I1,2)))**2+&
                !((CORDA(IP,3)-TA1(I1,3))*(CORDA(IP,1)-TA2(I1,1))-(CORDA(IP,1)-TA1(I1,1))*(CORDA(IP,3)-TA2(I1,3)))**2+&
                !((CORDA(IP,1)-TA1(I1,1))*(CORDA(IP,2)-TA2(I1,2))-(CORDA(IP,2)-TA1(I1,2))*(CORDA(IP,1)-TA2(I1,1)))**2)
                        
            !WRITE(*,*)NJT
            IF(NJT.GE.0.) THEN  !判断顺序
                TINT(:)=TA1(I1,:)
            ELSE
                TINT(:)=TA2(I1,:)
            END IF

            !TINT(:)=TA1(I1,:)

            50 CONTINUE
            !WRITE(*,*)"IP",IP,TINT(:)
            DO I2=1,DNODE(IP)
                    K(1)=INP(IP,I1,1)
                    K(2)=INP(IP,I2,1)
                    L(1)=INP(IP,I1,2)
                    L(2)=INP(IP,I2,2)
                    
                    IF(I1.NE.I2) THEN

                    IF(SQRT((TINT(1)-TA1(I2,1))**2+(TINT(2)-TA1(I2,2))**2+(TINT(3)-TA1(I2,3))**2).LE.1E-5.OR.&
                       SQRT((TINT(1)-TA2(I2,1))**2+(TINT(2)-TA2(I2,2))**2+(TINT(3)-TA2(I2,3))**2).LE.1E-5) THEN
                        
                        NJT=(CORDA(IP,1)-TINT(1))*(VEA(K(1),L(1),2)*VEA(K(2),L(2),3)-VEA(K(1),L(1),3)*VEA(K(2),L(2),2))+&
                            (CORDA(IP,2)-TINT(2))*(VEA(K(1),L(1),3)*VEA(K(2),L(2),1)-VEA(K(1),L(1),1)*VEA(K(2),L(2),3))+&    
                            (CORDA(IP,3)-TINT(3))*(VEA(K(1),L(1),1)*VEA(K(2),L(2),2)-VEA(K(1),L(1),2)*VEA(K(2),L(2),1))
                        
                        N=N+1
                        !WRITE(*,*)'NJT',NJT

                        IF(NJT.LT.0.) THEN
                            SGN=-1.
                        ELSE IF(NJT.GT.0.) THEN
                            SGN=1.
                        ELSE
                            SGN=0.
                        END IF
                    
                        NJ=VEA(K(1),L(1),1)*VEA(K(2),L(2),1)+VEA(K(1),L(1),2)*VEA(K(2),L(2),2)+VEA(K(1),L(1),3)*VEA(K(2),L(2),3)
                        
                        IF(ABS(NJ-1.0).LE.1.E-4) THEN
                            ALP(N)=PI !+SGN*ACOS(NJ)
                        ELSE IF(ABS(NJ+1.0).LE.1.E-4) THEN
                            ALP(N)=2.*PI
                        ELSE
                            ALP(N)=PI+SGN*ACOS(NJ)
                        END IF
                        
                        I1=I2
                        !WRITE(*,*)"I1",I1

                        IF(SQRT((TINT(1)-TA1(I2,1))**2+(TINT(2)-TA1(I2,2))**2+(TINT(3)-TA1(I2,3))**2).LE.1E-5) THEN
                            TINT(:)=TA2(I2,:)
                        ELSE
                            TINT(:)=TA1(I2,:)
                            !WRITE(*,*)"TA2"
                        END IF

                        IF(I1.EQ.1) THEN
                            GOTO 60
                        ELSE
                            GOTO 50
                        END IF

                    END IF
                        
                    END IF
            END DO

            60 CONTINUE
            SE(IP)=SUM(ALP(1:N))-(N-2)*PI
            SE(IP)=SE(IP)/2.
            !SE(IP)=0.5
            !WRITE(*,*)N,ALP(1:N)
        END IF
        !WRITE(*,101)IP,DN(IP),N,DNODE(IP),SE(IP)/PI2,CORD(IP,:)  
        !WRITE(*,*)SE(IP)/PI2

        !IF(ABS(CORD(IP,3)).LE.1.E-5) SE(IP)=SE(IP)*2.  !改动至此
        DO I1=1,NFX
            IP1=MWAP(I1)
            !SE(IP1)=SE(IP1)*2
        END DO
        
        WRITE(11,101)IP,DN(IP),N,DNODE(IP),SE(IP)/PI2,CORD(IP,:)
        SE(IP)=SE(IP)/PI2
!*************************************
        !SE(IP)=0.5
!*************************************        
        101FORMAT(3I4,5F12.4,9I9) 
        !END IF 
        !WRITE(*,*)SE(IP)

        20 CONTINUE                 
    END DO

    !NTPN=NBPOINT+NFPOINT
    !NTP=NBBLOCK+NFBLOCK
    NTPN=NTPN/ISYM
    NTP=NTP/ISYM

    DO I=1,NTPN
        CORD(I,:)=TCORD(I,:)
    END DO
    DEALLOCATE(TCORD)

    DEALLOCATE(CORDA,INEA,VEA,PANK)

    123 SE=0.5

END SUBROUTINE
    