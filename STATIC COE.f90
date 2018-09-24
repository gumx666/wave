SUBROUTINE SIMSTATIC
USE INPUTDATA    
!*****************************
!STA(1)=WFA
!STV(1,1)=V
!STV(1:3,2)=XB,YB,ZB
!*****************************  
    USE GREENMOD
    USE CUMOD
    REAL ::TEMP,STA,STV(3,6),STMA(3),STA1,WP(200,3)
    REAL,ALLOCATABLE:: CORDF(:,:,:),TCORD(:,:)
    INTEGER,ALLOCATABLE:: MARK(:),MARK1(:)
    INTEGER ::I,J,K,COUT,NBBLOCKT
    CHARACTER ::M,C
  
    ALLOCATE (MARK(NBPOINT),MARK1(NBPOINT),TCORD(NTPN,3))

    NBBLOCKT=NBBLOCK
    NBBLOCK=NBBLOCKW

    !OPEN(11,FILE='MASS.DAT')
  
    DO WHILE(.TRUE.)
        !READ(11,*)C
        IF(C=='C') EXIT
    END DO
    
    XGRA(1:3)=NSWIN(6,1:3)
    !READ(11,*)XGRA(1:3)
    !IF(ITTE.GE.2) XGRA=0.
    IF(MDFREE.EQ.1) THEN
!   重心变换
    TPX=XGRA(1)
    TPZ=XGRA(3)
    XGRA(3)=SINK-TPX*SIN(TRI)+TPZ*COS(TRI)
    XGRA(1)=TPX*COS(TRI)+TPZ*SIN(TRI)
    !XGRA(3)=SINK+TPZ
    END IF

    DO I=1,NTPN
        TCORD(I,:)=CORD(I,:)
        CORD(I,:)=PCORD(I,:)
    END DO
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO WHILE(.TRUE.)
        READ(11,*)C
        IF(C=='M') EXIT
    END DO
    READ(11,*)MASS(1,1)
    DO WHILE(.TRUE.)
        READ(11,*)C
        IF(C=='I') EXIT
    END DO
    READ(11,*)IXX,IYY,IZZ,IXY,IYZ,IXZ

    READ(11,*)C
    WRITE(*,*)C

    !READ(11,*)RESTO(1,1:6)
    !READ(11,*)RESTO(2,1:6)
    !READ(11,*)RESTO(3,1:6)
    !READ(11,*)RESTO(4,1:6)
    !READ(11,*)RESTO(5,1:6)
    !READ(11,*)RESTO(6,1:6)

    DO I=2,3
        MASS(I,I)=MASS(1,1)
    END DO
    MASS(4,4)=MASS(1,1)*IXX**2
    MASS(5,5)=MASS(1,1)*IYY**2
    MASS(6,6)=MASS(1,1)*IZZ**2  

    WRITE(*,112)XGRA
    WRITE(*,108)MASS(1,1),MASS(2,2),MASS(3,3)
    WRITE(*,107)MASS(4,4),MASS(5,5),MASS(6,6)
    108 FORMAT(1X,"MASS FROM DISTRUBITION   ",3X,3E15.4)
    107 FORMAT(1X,"INITIAL                  ",3X,3E15.4)
    112 FORMAT(1X,"CENTER OF GRAVATY        ",3X,3E15.4)

    WFA=0
    VOL=0
    XBUO=0
    LXX=0
    LYY=0
    LXY=0
    
    !COMPUTATE WATER PLANE AREA
    M='S'
    COUT=(NFX-1)/2
    DO I=1,COUT
        MARK1(I)=NBBLOCKP+NFXF/2*(NFY-1)/2+(NFY-1)/2*(I-1)+1
    END DO
  
    K=1
    MARK=MARK1
    DO I=2,COUT
        DO J=1,I-1
            IF(MARK1(J).EQ.MARK1(I)) GOTO 12
        END DO
        K=K+1
        MARK(K)=MARK1(I)
    12 CONTINUE
    END DO


    ALLOCATE (CORDF(NBBLOCK,3,9),CORL(9,3))
    DO I=1,COUT 
        !CORDF(I,:,1)=CORD(INE(MARK(I),3),:)
        !CORDF(I,:,2)=CORD(INE(MARK(I),2),:)  
        !CORDF(I,:,3)=CORD(INE(MARK(I),1),:) 
        !CORDF(I,:,4)=CORDF(I,:,3)
        !CORDF(I,:,5)=CORDF(I,:,3)
        !CORDF(I,:,6)=CORDF(I,:,2)
        !CORDF(I,:,7)=CORDF(I,:,1)
        !CORDF(I,:,8)=CORDF(I,:,1)
        !CORDF(I,:,9)=CORDF(I,:,2)
        !CORDF(I,2,5)=0.
        !CORDF(I,2,6)=0.
        !CORDF(I,2,7)=0.
        !CORDF(I,2,4)=CORDF(I,2,3)/2.
        !CORDF(I,2,9)=CORDF(I,2,2)/2.
        !CORDF(I,2,8)=CORDF(I,2,1)/2.
        !CORDF(I,:,1)
    END DO

    DO 30 IE=1,COUT
!C			WRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!c			依次将每一面元的四个顶点的坐标赋值到corl(4,3)中
        DO 40 J1=1,9
            CORL(J1,:)=CORDF(IE,:,J1)
        40 CONTINUE
        
        CALL BVOLUM_9(XGRA,M,CORL,STA,STV(:,:),STMA)
        
        WFA=WFA+STA
        LXX=LXX+STMA(1)
        LYY=LYY+STMA(2)
        LXY=LXY+STMA(3)
    30 CONTINUE
    
    WFA=0.
    LXX=0.
    LYY=0.
    LXY=0.
    DO I=1,NFX-1
        !IF(MDFREE.EQ.2) THEN
        WP(I,1)=0.5*(HULL_WP(I,1)+HULL_WP0(I,1))
        WP(I,2)=0.5*(HULL_WP(I,2)+HULL_WP0(I,2))
        WP(I,3)=0.
        WP(I+1,1)=0.5*(HULL_WP(I+1,1)+HULL_WP0(I+1,1))
        WP(I+1,2)=0.5*(HULL_WP(I+1,2)+HULL_WP0(I+1,2))
        WP(I+1,3)=0.
        !ELSE
        !IF(MDFREE.EQ.1) THEN
        
        WP(I,1)=HULL_WP(I,1)
        WP(I,2)=HULL_WP(I,2)
        WP(I+1,1)=HULL_WP(I+1,1)
        WP(I+1,2)=HULL_WP(I+1,2)
        !END IF

        WFA=WFA+(WP(I+1,2)+WP(I,2))*(WP(I+1,1)-WP(I,1))/2.
        LXX=LXX+((WP(I+1,2)+WP(I,2))/2.)**3*(WP(I+1,1)-WP(I,1))/3.
        LYY=LYY+((WP(I+1,1)+WP(I,1))/2.)**2*((WP(I+1,2)+WP(I,2))/2.)*(WP(I+1,1)-WP(I,1))
        LXY=LXY+((WP(I+1,1)+WP(I,1))/2.)*((WP(I+1,2)+WP(I,2))/2.)*(WP(I+1,1)-WP(I,1))
    END DO


    WFA=WFA*2.
    LXX=2.*LXX;LYY=2.*LYY;LXY=2.*LXY
  
    !LXX=6.0788E5
    !LYY=3.801E7
    !LXY=0
    
    !XFA=-16.456-XGRA(1)
    XFA=0
    YFA=0

    WRITE(*,109)LXX,LYY,LXY
    109 FORMAT(1X,"LXX  ",E15.4,"     LYY  ",E15.4,"     LXY  ",E15.4)

    WRITE(*,*)"ELEMENT NUMBER_WFA ",COUT
    WRITE(*,*)"WATERPLANE_AREA=      ",WFA
    
    !COMPUTATE DISPLACEMENT
    

    M='V'
    GOTO 333 !TEST
    !M='S'
    COUT=0
    DO I=1,NBBLOCK+NFBLOCK
        TEMP=0
        DO J=1,PAN(I) !NODENUM
            TEMP=TEMP+ABS(CORD(INE(I,J),3))
        END DO
        IF(TEMP.LE.1E-5.OR.ABS(TEMP-NODENUM*RD).LE.1E-5) GOTO 21    
        COUT=COUT+1
        MARK1(COUT)=I
        21 CONTINUE
    END DO
  
    K=1
    MARK=MARK1
    DO I=2,COUT
        DO J=1,I-1
            IF(MARK1(J).EQ.MARK1(I)) GOTO 22
        END DO
        K=K+1
        MARK(K)=MARK1(I)
        22 CONTINUE
    END DO

    COUT=K

    DO I=1,COUT
        DO J=1,PAN(MARK(I)) !NODENUM 
            CORDF(I,:,J)=CORD(INE(MARK(I),J),:)
        END DO 
    END DO

    333 CONTINUE

    SAREA=0
    MASS(1,1)=0.
    !WRITE(*,*)"COUT",COUT,NBBLOCK
    DO 33 IE=1,NBBLOCK !COUT
!	    WRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!		依次将每一面元的四个顶点的坐标赋值到corl(4,3)中
        DO 44 J1=1,PAN(IE)
            !CORL(J1,:)=CORDF(IE,:,J1)
            CORL(J1,:)=CORD(INE(IE,J1),:)
        44 CONTINUE

        IF(PAN(IE).EQ.4) CALL BVOLUM_4(M,CORL,STA,STV(:,:),STMA)
        IF(PAN(IE).EQ.9) THEN
            CALL BVOLUM_9V(M,CORL,STA1)
            CALL BVOLUM_9(XGRA,M,CORL,STA,STV(:,:),STMA)
        END IF
        
        IF(PAN(IE).EQ.8) THEN
            !CALL BVOLUM_8V(M,CORL(1:PAN(IE),:),CORDO(ICORD(IE),:,:),STA1)
            CALL BVOLUM_8V_OLD(M,CORL(1:PAN(IE),:),CORDO(ICORD(IE),:,:),STA1)
            !WRITE(*,*)"STA1",STA1
            !STOP
            CALL BVOLUM_8_OLD(M,CORL(1:PAN(IE),:),STA,STV(:,:),STMA)
        END IF
        

        SAREA=SAREA+STA1*2.
        VOL=VOL+STV(1,1)
        XBUO(1)=XBUO(1)+STV(1,2)
        XBUO(3)=XBUO(3)+STV(3,2)
        MASS(1,1)=MASS(1,1)+STV(1,2)

        !WRITE(*,999)IE,STV(1,1),STV(1,2),STV(1,2),STV(1,2)
        999FORMAT(I8,6F15.6)
    33 CONTINUE

    IF(ITTE.EQ.1) SAREA0=SAREA

        WRITE(*,*)"____________________________"

    VOL=ABS(VOL)*2.  
    XBUO(1)=XBUO(1)*2./VOL
    XBUO(3)=XBUO(3)*2./VOL
    MASS(1,1)=MASS(1,1)*2*P1

    !XGRA(1)=XBUO(1)
    !程序计算

    !MASS(1,1)=VOL*P1
    MASS(2,2)=VOL*P1
    MASS(3,3)=MASS(2,2)
    MASS(4,4)=MASS(2,2)*IXX**2
    MASS(5,5)=MASS(2,2)*IYY**2
    MASS(6,6)=MASS(2,2)*IZZ**2      

    IF(ITTE.EQ.1) THEN
        MASS0=MASS
    END IF
    
    !IF(ITTE.EQ.1) THEN
    RESTO(3,3)=P1*G*WFA
    RESTO(5,5)=P1*G*VOL*(XBUO(3)-XGRA(3))+P1*G*LYY !*WFA
    RESTO(3,5)=-P1*G*LXY
    RESTO(5,3)=-P1*G*LXY
    !!!!!!!!!!!!!!
   
    !SOFT SPRING
    RESTO(1,1)=(1.5*MASS(1,1)+MASS(1,2))*(CFR+U*WN)
    RESTO(2,2)=(1.6*MASS(1,1)+MASS(1,2))*(CFR+U*WN)
    RESTO(6,6)=(1.6*MASS(6,6)+MASS(1,2))*(CFR+U*WN)
    !END IF

    WRITE(*,*)"MASS"
    WRITE(*,118)MASS(1,:)
    WRITE(*,118)MASS(2,:)
    WRITE(*,118)MASS(3,:)
    WRITE(*,118)MASS(4,:)
    WRITE(*,118)MASS(5,:)
    WRITE(*,118)MASS(6,:)

    WRITE(*,*)"RESTO"
    WRITE(*,118)RESTO(1,:)
    WRITE(*,118)RESTO(2,:)
    WRITE(*,118)RESTO(3,:)
    WRITE(*,118)RESTO(4,:)
    WRITE(*,118)RESTO(5,:)
    WRITE(*,118)RESTO(6,:)
    118 FORMAT(1X,6E13.4)

    DO I=1,NTPN
        !CORD(I,1)=CORD(I,1)-XGRA(1)
        !CORD(I,2)=CORD(I,2)-XGRA(2)
        !CORD(I,3)=CORD(I,3)-XGRA(3)
    END DO
    !XBUO=XBUO-XGRA

    WRITE(*,*)"SURFACE_AREA=      ",SAREA
    !WRITE(*,*)"ELEMENT NUMBER_VOL=",COUT
    WRITE(*,105),VOL,XBUO
    105 FORMAT(1X,"DISPLACEMENT=  ",E10.4,3F10.5)
    WRITE(*,106),VOL*P1
    106 FORMAT(1X,"DISP    MASS=  ",F16.4)
    WRITE(*,110),LXX/VOL
    110 FORMAT(1X,"BMxx=          ",4F10.4)
    WRITE(*,111),LYY/VOL
    111 FORMAT(1X,"BMyy=          ",4F10.4)
    WRITE(*,122),Lxy/WFA
    122 FORMAT(1X,"XF=          ",4F10.4)
    !CHANGE COORDINATE

    !自由面抬高
    DO I=NBPOINT+1,NTPN
        !CORD(I,3)=CORD(I,3)+RD/40.
    END DO


    DO I=1,NTPN
        CORD(I,:)=TCORD(I,:)
    END DO
    NBBLOCK=NBBLOCKT

    CLOSE(11)
    DEALLOCATE (TCORD)
    DEALLOCATE (CORDF,MARK,MARK1)
    DEALLOCATE (CORL)
END SUBROUTINE


SUBROUTINE SIMSTATIC_OLD
!*****************************
!STA(1)=WFA
!STV(1,1)=V
!STV(1:3,2)=XB,YB,ZB
!*****************************  
    USE GREENMOD
    USE CUMOD
    USE INPUTDATA    
    REAL ::TEMP,STA,STV(3,6),STMA(3),STA1
    REAL,ALLOCATABLE:: CORDF(:,:,:),TCORD(:,:)
    INTEGER,ALLOCATABLE:: MARK(:),MARK1(:)
    INTEGER ::I,J,K,COUT
    CHARACTER ::M,C
  
    ALLOCATE (MARK(NBPOINT),MARK1(NBPOINT),TCORD(NTPN,3))


    !OPEN(11,FILE='MASS.DAT')
  
    !DO WHILE(.TRUE.)
        !READ(11,*)C
        !IF(C=='C') EXIT
    !END DO
    
    !READ(11,*)XGRA(1:3)
    !XGRA=0.
    XGRA(1:3)=NSWIN(6,1:3)
    IF(MDFREE.EQ.1) THEN
!   重心变换
    TPX=XGRA(1)
    TPZ=XGRA(3)
    XGRA(3)=SINK-TPX*SIN(TRI)+TPZ*COS(TRI)
    XGRA(1)=TPX*COS(TRI)+TPZ*SIN(TRI)
    !XGRA(3)=SINK+TPZ
    END IF

    DO I=1,NTPN
        TCORD(I,:)=CORD(I,:)
        CORD(I,:)=PCORD(I,:)
    END DO
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !DO WHILE(.TRUE.)
    !    READ(11,*)C
    !    IF(C=='M') EXIT
    !END DO
    !READ(11,*)MASS(1,1)
    !DO WHILE(.TRUE.)
    !    READ(11,*)C
    !    IF(C=='I') EXIT
    !END DO
    !READ(11,*)IXX,IYY,IZZ,IXY,IYZ,IXZ

    !READ(11,*)C
    !WRITE(*,*)C

    !READ(11,*)RESTO(1,1:6)
    !READ(11,*)RESTO(2,1:6)
    !READ(11,*)RESTO(3,1:6)
    !READ(11,*)RESTO(4,1:6)
    !READ(11,*)RESTO(5,1:6)
    !READ(11,*)RESTO(6,1:6)

    DO I=2,3
        MASS(I,I)=MASS(1,1)
    END DO
    MASS(4,4)=MASS(1,1)*IXX**2
    MASS(5,5)=MASS(1,1)*IYY**2
    MASS(6,6)=MASS(1,1)*IZZ**2  

    !WRITE(*,112)XGRA
    !WRITE(*,108)MASS(1,1),MASS(2,2),MASS(3,3)
    !WRITE(*,107)MASS(4,4),MASS(5,5),MASS(6,6)
    108 FORMAT(1X,"MASS FROM DISTRUBITION   ",3X,3E15.4)
    107 FORMAT(1X,"INITIAL                  ",3X,3E15.4)
    112 FORMAT(1X,"CENTER OF GRAVATY        ",3X,3E15.4)

    WFA=0
    VOL=0
    XBUO=0
    LXX=0
    LYY=0
    LXY=0
    
    !COMPUTATE WATER PLANE AREA
    M='S'
    COUT=(NFX-1)/2
    DO I=1,COUT
        MARK1(I)=NBBLOCKP+NFXF/2*(NFY-1)/2+(NFY-1)/2*(I-1)+1
    END DO
  
    K=1
    MARK=MARK1
    DO I=2,COUT
        DO J=1,I-1
            IF(MARK1(J).EQ.MARK1(I)) GOTO 12
        END DO
        K=K+1
        MARK(K)=MARK1(I)
    12 CONTINUE
    END DO


    ALLOCATE (CORDF(NBBLOCK,3,9),CORL(9,3))
    DO I=1,COUT 
        CORDF(I,:,1)=CORD(INE(MARK(I),3),:)
        CORDF(I,:,2)=CORD(INE(MARK(I),2),:)  
        CORDF(I,:,3)=CORD(INE(MARK(I),1),:) 
        CORDF(I,:,4)=CORDF(I,:,3)
        CORDF(I,:,5)=CORDF(I,:,3)
        CORDF(I,:,6)=CORDF(I,:,2)
        CORDF(I,:,7)=CORDF(I,:,1)
        CORDF(I,:,8)=CORDF(I,:,1)
        CORDF(I,:,9)=CORDF(I,:,2)
        CORDF(I,2,5)=0.
        CORDF(I,2,6)=0.
        CORDF(I,2,7)=0.
        CORDF(I,2,4)=CORDF(I,2,3)/2.
        CORDF(I,2,9)=CORDF(I,2,2)/2.
        CORDF(I,2,8)=CORDF(I,2,1)/2.

        CORDF(I,3,:)=0.
    END DO

    DO 30 IE=1,COUT
!C			WRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!c			依次将每一面元的四个顶点的坐标赋值到corl(4,3)中
        DO 40 J1=1,9
            CORL(J1,:)=CORDF(IE,:,J1)
        40 CONTINUE
        
        CALL BVOLUM_9(XGRA,M,CORL,STA,STV(:,:),STMA)
        
        WFA=WFA+STA
        LXX=LXX+STMA(1)
        LYY=LYY+STMA(2)
        LXY=LXY+STMA(3)
    30 CONTINUE
    
    WFA=WFA*2.
    LXX=2.*LXX;LYY=2.*LYY;LXY=2.*LXY
  
    !LXX=6.0788E5
    !LYY=3.801E7
    !LXY=0
    
    !XFA=-16.456-XGRA(1)
    XFA=0
    YFA=0

    !WRITE(*,109)LXX,LYY,LXY
    109 FORMAT(1X,"LXX  ",E15.4,"     LYY  ",E15.4,"     LXY  ",E15.4)

    !WRITE(*,*)"ELEMENT NUMBER_WFA ",COUT
    !WRITE(*,*)"WATERPLANE_AREA=      ",WFA
    
    !COMPUTATE DISPLACEMENT
    

    M='V'
    GOTO 333 !TEST
    !M='S'
    COUT=0
    DO I=1,NBBLOCK+NFBLOCK
        TEMP=0
        DO J=1,PAN(I) !NODENUM
            TEMP=TEMP+ABS(CORD(INE(I,J),3))
        END DO
        IF(TEMP.LE.1E-5.OR.ABS(TEMP-NODENUM*RD).LE.1E-5) GOTO 21    
        COUT=COUT+1
        MARK1(COUT)=I
        21 CONTINUE
    END DO
  
    K=1
    MARK=MARK1
    DO I=2,COUT
        DO J=1,I-1
            IF(MARK1(J).EQ.MARK1(I)) GOTO 22
        END DO
        K=K+1
        MARK(K)=MARK1(I)
        22 CONTINUE
    END DO

    COUT=K

    DO I=1,COUT
        DO J=1,PAN(MARK(I)) !NODENUM 
            CORDF(I,:,J)=CORD(INE(MARK(I),J),:)
        END DO 
    END DO

    333 CONTINUE

    SAREA=0
    MASS(1,1)=0.
    !WRITE(*,*)"COUT",COUT,NBBLOCK
    DO 33 IE=1,NBBLOCK !COUT
!	    WRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!		依次将每一面元的四个顶点的坐标赋值到corl(4,3)中
        DO 44 J1=1,PAN(IE)
            !CORL(J1,:)=CORDF(IE,:,J1)
            CORL(J1,:)=CORD(INE(IE,J1),:)
        44 CONTINUE

        IF(PAN(IE).EQ.4) CALL BVOLUM_4(M,CORL,STA,STV(:,:),STMA)
        IF(PAN(IE).EQ.9) THEN
            CALL BVOLUM_9V(M,CORL,STA1)
            CALL BVOLUM_9(XGRA,M,CORL,STA,STV(:,:),STMA)
        END IF
        
        IF(PAN(IE).EQ.8) THEN
            !CALL BVOLUM_8V(M,CORL(1:PAN(IE),:),CORDO(ICORD(IE),:,:),STA1)
            CALL BVOLUM_8V_OLD(M,CORL(1:PAN(IE),:),CORDO(ICORD(IE),:,:),STA1)
            !WRITE(*,*)"STA1",STA1
            !STOP
            CALL BVOLUM_8_OLD(M,CORL(1:PAN(IE),:),STA,STV(:,:),STMA)
        END IF
        

        SAREA=SAREA+STA1*2.
        VOL=VOL+STV(1,1)
        XBUO(1)=XBUO(1)+STV(1,2)
        XBUO(3)=XBUO(3)+STV(3,2)
        MASS(1,1)=MASS(1,1)+STV(1,2)

        !WRITE(*,999)IE,STV(1,1),STV(1,2),STV(1,2),STV(1,2)
        999FORMAT(I8,6F15.6)
    33 CONTINUE

    IF(ITTE.EQ.1) SAREA0=SAREA

        !WRITE(*,*)"____________________________"

    VOL=ABS(VOL)*2.  
    XBUO(1)=XBUO(1)*2./VOL
    XBUO(3)=XBUO(3)*2./VOL
    MASS(1,1)=MASS(1,1)*2*P1

    !XGRA(1)=XBUO(1)
    !程序计算

    !MASS(1,1)=VOL*P1
    MASS(2,2)=VOL*P1
    MASS(3,3)=MASS(2,2)
    MASS(4,4)=MASS(2,2)*IXX**2
    MASS(5,5)=MASS(2,2)*IYY**2
    MASS(6,6)=MASS(2,2)*IZZ**2      

    IF(ITTE.EQ.1) THEN
        MASS0=MASS
    END IF
    
    !IF(ITTE.EQ.1) THEN
    RESTO(3,3)=P1*G*WFA
    RESTO(5,5)=P1*G*VOL*(XBUO(3)-XGRA(3))+P1*G*LYY !*WFA
    RESTO(3,5)=-P1*G*LXY
    RESTO(5,3)=-P1*G*LXY
    !!!!!!!!!!!!!!
   
    !SOFT SPRING
    RESTO(1,1)=(1.5*MASS(1,1)+MASS(1,2))*(CFR+U*WN)
    RESTO(2,2)=(1.6*MASS(1,1)+MASS(1,2))*(CFR+U*WN)
    RESTO(6,6)=(1.6*MASS(6,6)+MASS(1,2))*(CFR+U*WN)
    !END IF

    !WRITE(*,*)"MASS"
    !WRITE(*,118)MASS(1,:)
    !WRITE(*,118)MASS(2,:)
    !WRITE(*,118)MASS(3,:)
    !WRITE(*,118)MASS(4,:)
    !WRITE(*,118)MASS(5,:)
    !WRITE(*,118)MASS(6,:)

    !WRITE(*,*)"RESTO"
    !WRITE(*,118)RESTO(1,:)
    !WRITE(*,118)RESTO(2,:)
    !WRITE(*,118)RESTO(3,:)
    !WRITE(*,118)RESTO(4,:)
    !WRITE(*,118)RESTO(5,:)
    !WRITE(*,118)RESTO(6,:)
    118 FORMAT(1X,6E13.4)

    DO I=1,NTPN
        !CORD(I,1)=CORD(I,1)-XGRA(1)
        !CORD(I,2)=CORD(I,2)-XGRA(2)
        !CORD(I,3)=CORD(I,3)-XGRA(3)
    END DO
    !XBUO=XBUO-XGRA

    !WRITE(*,*)"SURFACE_AREA=      ",SAREA
    !WRITE(*,*)"ELEMENT NUMBER_VOL=",COUT
    !WRITE(*,105),VOL,XBUO
    105 FORMAT(1X,"DISPLACEMENT=  ",E10.4,3F10.5)
    !WRITE(*,106),VOL*P1
    106 FORMAT(1X,"DISP    MASS=  ",F16.4)
    !WRITE(*,110),LXX/VOL
    110 FORMAT(1X,"BMxx=          ",4F10.4)
    !WRITE(*,111),LYY/VOL
    111 FORMAT(1X,"BMyy=          ",4F10.4)
    !WRITE(*,122),Lxy/WFA
    122 FORMAT(1X,"XF=          ",4F10.4)
    !CHANGE COORDINATE

    !自由面抬高
    DO I=NBPOINT+1,NTPN
        !CORD(I,3)=CORD(I,3)+RD/40.
    END DO


    DO I=1,NTPN
        CORD(I,:)=TCORD(I,:)
    END DO

    CLOSE(11)
    DEALLOCATE (TCORD)
    DEALLOCATE (CORDF,MARK,MARK1)
    DEALLOCATE (CORL)
END SUBROUTINE

SUBROUTINE BVOLUM_4(M,CORDA,TWFA,TVOL,TAMON)
!c    ***********************************************      
!c      calculate the WFA in the matrix     
!c      corda(4,3) the cord. of 4 nodes (clockwise)      
!c      conm(4,3) the normal vector of 4 nodes          4---3
!c      ps(3) the cord. of field point                  |   |
!c      aij(4): panel source effect                     1---2
!c      hij(4): panel dipole effect
!c   ************************************************
    USE NMLMOD
    REAL:: XI,ETA,AJ,TEMP(8),TEMP1(8),TWFA,TVOL(3,6),TAMON(3)
    REAL:: corda(4,3),tcorda(4,3),ps(3),conm(4,3),pr(3)
    REAL:: r(4),srr(3),ay(4),hy(4)
    REAL:: sn(4),snx(4),sne(4),sr(3),dn(3)
    COMMON /gwa/ga(400),gw(400)
    CHARACTER M

!c	初始化每一面元的四个顶点的源效应aij（4）和偶极子效应hij（4）
    IF(M.EQ.'V') THEN
        tcorda=corda
        corda(:,2)=0
    END IF 

    twfa=0
    tvol=0
    temp=0
    temp1=0
    TAMON=0
    ngus=21

!c	k=2+3+...+(n-1)+1=n*(n-1)/2
    k=ngus*(ngus-1)/2
!c  gauss-legendre numerical integration
    
!c	对x积分     
    do 40 ix=1,ngus
	    xi=ga(k+ix-1)
!c	对y积分
		do 40 iy=1,ngus
            temp=0
            eta=ga(k+iy-1)
!c	调用isopar，求得形状函数J及单位法向矢量等     
			call isopar_4(xi,eta,corda,sn,snx,sne,aj,sr,dn)
      
			do 60 j=1,4
				!temp(1)=temp(1)+sn(j)
                temp(1)=1. !temp(1)+sn(j)
                temp(2)=temp(2)+ABS(tcorda(j,2))*sn(j)
                temp(3)=temp(3)+corda(j,1)*ABS(tcorda(j,2))*sn(j)
                temp(4)=temp(4)+corda(j,3)*ABS(tcorda(j,2))*sn(j)
                
                temp(5)=temp(5)+corda(j,2)**2*sn(j)
                temp(6)=temp(6)+corda(j,1)**2*sn(j)
                temp(7)=temp(7)+corda(j,1)*corda(j,2)*sn(j)

                temp(8)=temp(8)+sn(j)
    60		continue
            temp1(1)=temp(1)*gw(k+iy-1)*aj
            temp1(2)=temp(2)*gw(k+iy-1)*aj
            temp1(3)=temp(3)*gw(k+iy-1)*aj
            temp1(4)=temp(4)*gw(k+iy-1)*aj
!c
            temp1(5)=temp(5)*gw(k+iy-1)*aj
            temp1(6)=temp(6)*gw(k+iy-1)*aj
            temp1(7)=temp(7)*gw(k+iy-1)*aj
			!do 40 j=1,4
			twfa=twfa+temp1(1)*gw(k+ix-1)
			tvol(1,1)=tvol(1,1)+temp1(2)*gw(k+ix-1)
            tvol(1,2)=tvol(1,2)+temp1(3)*gw(k+ix-1)
            tvol(3,2)=tvol(3,2)+temp1(4)*gw(k+ix-1)

			TAMON(1)=TAMON(1)+temp1(5)*gw(k+ix-1)
            TAMON(2)=TAMON(2)+temp1(6)*gw(k+ix-1)
            TAMON(3)=TAMON(3)+temp1(7)*gw(k+ix-1)
    40continue
      return
END SUBROUTINE

SUBROUTINE BVOLUM_8(M,CORDA,TWFA,TVOL,TAMON)
!c    ***********************************************      
!c     calculate the aij and hij terms in the matrix   
!c      the ps is one of the nodes of element          4---7---3
!c      cord(8,3) the cord. of 8 nodes (clockwise)     |       |
!c      conm(8,3) the normal vector of 8 nodes         8       6
!c      ps(3) the cord. of fiele point                 |       |
!c      aij(8): panel source effect                    1---5---2
!c      hij(8): panel dipole effect
!c   ************************************************
!cdp   implicit REAL(a-h,o-z)
    USE NMLMOD
    REAL:: XI,ETA,TEMP(8),TEMP1(8),TWFA,TVOL(3,6),TAMON(3)
    REAL:: corda(8,3),tcorda(8,3),ps(3),conm(8,3),pr(3)
    REAL:: r(8),srr(3),ay(8),hy(8),p1,p2
    COMMON /gwa/ga(400),gw(400)
    CHARACTER M

      REAL:: XI1,XI2,XI3,XI4,ETA1,ETA2,ETA3,ETA4
      REAL:: aij1(8),hij1(8),cord1(8,3),conm1(8,3)
      REAL:: aj(4)
      REAL:: snx(8),sne(8),sm(8,4)
      REAL:: sn1(8),sn2(8),sn3(8),sn(8,4),sn4(8)
      equivalence (sn(1,1),sn1),(sn(1,2),sn2),(sn(1,3),sn3),(sn(1,4),sn4)
      REAL:: sr1(3),sr2(3),sr3(3),sr(3,4),sr4(3)
      equivalence (sr(1,1),sr1),(sr(1,2),sr2),(sr(1,3),sr3),(sr(1,4),sr4)
      REAL:: dn1(3),dn2(3),dn3(3),dn(3,4),dn4(3)
      equivalence (dn(1,1),dn1),(dn(1,2),dn2),(dn(1,3),dn3),(dn(1,4),dn4)
      common nbp,nfp,nbe,nfe,np,mp,nf,radius
      REAL:: CORDT(3,2)
      REAL:: SNT(3,4),SNXT(3),SNET(3),AJT(4),SRT(3,4),DNT(3,4)

      !write(*,*)"test"

      D1=0.9
      D2=D1
      D3=D2
      D4=D3

      r1=0.
      isn=2

!	初始化每一面元的四个顶点的源效应aij（4）和偶极子效应hij（4）
    IF(M.EQ.'V') THEN
        tcorda=corda
        !corda(:,2)=0
        !WRITE(*,*)"corda",corda(:,3)
    END IF 

    twfa=0
    tvol=0
    temp=0
    temp1=0
    TAMON=0

    ngus=21

! 	k=2+3+...+(n-1)+1=n*(n-1)/2
    KGUS=ngus*(ngus-1)/2
!   gauss-legendre numerical integration
    
         
	SELECT CASE(ISN)
    CASE(1)
        DXI=-D1;DET=-D3
    CASE(2)
        DXI=D2;DET=-D3
    CASE(3)
        DXI=D2;DET=D4
    CASE(4)
        DXI=-D1;DET=D4
    CASE(5)
        DXI=0.;DET=-D3
    CASE(6)
        DXI=D2;DET=0.
    CASE(7)
        DXI=0.;DET=D4
    CASE(8) 
        DXI=-D1;DET=0.
    END SELECT       
            
    DO IX=1,NGUS
        XI=GA(KGUS+IX-1)
        P1=0.5*(1.0+GA(KGUS+IX-1))
        DO IY=1,NGUS
            ETA=GA(KGUS+IY-1)
            P2=0.5*(1.0+GA(KGUS+IY-1))
            temp=0.
! 	调用isopar，求得形状函数J及单位法向矢量等    

            CORDT(1,1)=DXI;CORDT(2,1)=-1.;CORDT(3,1)=-1.
            CORDT(1,2)=DET;CORDT(2,2)=1.;CORDT(3,2)=-1.
            !CALL isopar_3(xi,eta,cordT,snT(:,1),snxT,sneT,ajT(1),srT(:,1),dnT(:,1))
            CORDT(1,1)=DXI;CORDT(2,1)=-1.;CORDT(3,1)=1.
            CORDT(1,2)=DET;CORDT(2,2)=-1.;CORDT(3,2)=-1.
            !CALL isopar_3(xi,eta,cordT,snT(:,2),snxT,sneT,ajT(2),srT(:,2),dnT(:,2))
            CORDT(1,1)=DXI;CORDT(2,1)=1.;CORDT(3,1)=1.
            CORDT(1,2)=DET;CORDT(2,2)=-1.;CORDT(3,2)=1.
            !CALL isopar_3(xi,eta,cordT,snT(:,3),snxT,sneT,ajT(3),srT(:,3),dnT(:,3))
            CORDT(1,1)=DXI;CORDT(2,1)=1.;CORDT(3,1)=-1.
            CORDT(1,2)=DET;CORDT(2,2)=1.;CORDT(3,2)=1.
            !CALL isopar_3(xi,eta,cordT,snT(:,4),snxT,sneT,ajT(4),srT(:,4),dnT(:,4))


            XI1=(1.-P1)*(DXI)+P1*(1-P2)*(-1.)+P1*P2*(-1.)
            ETA1=(1.-P1)*(DET)+P1*(1-P2)*(1.)+P1*P2*(-1.)
            XI2=(1.-P1)*(DXI)+P1*(1-P2)*(-1.)+P1*P2*(1.)
            ETA2=(1.-P1)*(DET)+P1*(1-P2)*(-1.)+P1*P2*(-1.)
            XI3=(1.-P1)*(DXI)+P1*(1-P2)*(1.)+P1*P2*(1.)
            ETA3=(1.-P1)*(DET)+P1*(1-P2)*(-1.)+P1*P2*(1.)
            XI4=(1.-P1)*(DXI)+P1*(1-P2)*(1.)+P1*P2*(-1.)
            ETA4=(1.-P1)*(DET)+P1*(1-P2)*(1.)+P1*P2*(1.)

            CALL ISOPAR_8(XI1,ETA1,CORDA,SN1,SNX,SNE,AJ(1),SR1,DN(:,1))
            CALL ISOPAR_8(XI2,ETA2,CORDA,SN2,SNX,SNE,AJ(2),SR2,DN(:,2))
            CALL ISOPAR_8(XI3,ETA3,CORDA,SN3,SNX,SNE,AJ(3),SR3,DN(:,3))
	        CALL ISOPAR_8(XI4,ETA4,CORDA,SN4,SNX,SNE,AJ(4),SR4,DN(:,4)) 
 
            !WRITE(*,209)SR,AJ
            209FORMAT(3F15.6,E15.8) 

			do j=1,8
                DO k=1,4

                    IF(K.EQ.1) AM=ABS(-1.-DXI)
                    IF(K.EQ.2) AM=ABS(-1.-DET)
                    IF(K.EQ.3) AM=ABS(1.-DXI)
                    IF(K.EQ.4) AM=ABS(1.-DET)

                    !WRITE(*,*)AM/2.*P1,AJT(K)

                temp(1)=temp(1)+sn(j,K)*aj(k)*AJT(K)
                !temp(2)=temp(2)+ABS(tcorda(j,2))*sn(j)
                !temp(3)=temp(3)+corda(j,1)*ABS(tcorda(j,2))*sn(j)
                !temp(4)=temp(4)+corda(j,3)*ABS(tcorda(j,2))*sn(j)
                
                temp(2)=temp(2)+tcorda(j,3)*dn(3,K)*sn(j,K)*aj(k)*AJT(K)
                temp(3)=temp(3)-corda(j,1)*tcorda(j,3)*dn(3,K)*sn(j,K)*aj(k)*AJT(K)
                temp(4)=temp(4)-corda(j,3)*tcorda(j,2)*dn(2,K)*sn(j,K)*aj(k)*AJT(K)
                

                temp(5)=temp(5)+corda(j,2)**2*sn(j,K)*aj(k)*AJT(K)
                temp(6)=temp(6)+corda(j,1)**2*sn(j,K)*aj(k)*AJT(K)
                !temp(7)=temp(7)+corda(j,1)*corda(j,2)*sn(j)
                temp(7)=temp(7)+corda(j,1)*sn(j,K)*aj(k)*AJT(K)
                temp(8)=temp(8)+sn(j,K)*aj(k)*AJT(K)
                END DO
            END DO

            temp1(1)=temp(1)*gw(KGUS+iy-1)
            temp1(2)=temp(2)*gw(KGUS+iy-1)
            temp1(3)=temp(3)*gw(KGUS+iy-1)
            temp1(4)=temp(4)*gw(KGUS+iy-1)
!c
            temp1(5)=temp(5)*gw(KGUS+iy-1)
            temp1(6)=temp(6)*gw(KGUS+iy-1)
            temp1(7)=temp(7)*gw(KGUS+iy-1)

			!do 40 j=1,4
			twfa=twfa+temp1(1)*gw(KGUS+ix-1)
			tvol(1,1)=tvol(1,1)+temp1(2)*gw(KGUS+ix-1)
            tvol(1,2)=tvol(1,2)+temp1(3)*gw(KGUS+ix-1)
            tvol(3,2)=tvol(3,2)+temp1(4)*gw(KGUS+ix-1)

			TAMON(1)=TAMON(1)+temp1(5)*gw(KGUS+ix-1)
            TAMON(2)=TAMON(2)+temp1(6)*gw(KGUS+ix-1)
            TAMON(3)=TAMON(3)+temp1(7)*gw(KGUS+ix-1)
        end do
    end do

      return
END SUBROUTINE

SUBROUTINE BVOLUM_8_OLD(M,CORDA,TWFA,TVOL,TAMON)
!c    ***********************************************      
!c      calculate the WFA in the matrix     
!c      corda(9,3) the cord. of 9 nodes (clockwise)      
!c      conm(9,3) the normal vector of 9 nodes         
!c      ps(3) the cord. of field point                 
!c      aij(9): panel source effect                     
!c      hij(9): panel dipole effect
!c   ************************************************
    USE NMLMOD
    REAL:: XI,ETA,AJ,TEMP(8),TEMP1(8),TWFA,TVOL(3,6),TAMON(3)
    REAL:: corda(8,3),tcorda(8,3),ps(3),conm(8,3),pr(3)
    REAL:: r(8),srr(3),ay(8),hy(8)
    REAL:: sn(8),snx(8),sne(8),sr(3),dn(3)
    COMMON /gwa/ga(400),gw(400)
    CHARACTER M

!c	初始化每一面元的四个顶点的源效应aij（4）和偶极子效应hij（4）
    IF(M.EQ.'V') THEN
        tcorda=corda
        !corda(:,2)=0
        !WRITE(*,*)"corda",corda(:,3)
    END IF 

    twfa=0
    tvol=0
    temp=0
    temp1=0
    TAMON=0
    ngus=11

! 	k=2+3+...+(n-1)+1=n*(n-1)/2
    k=ngus*(ngus-1)/2
!   gauss-legendre numerical integration
    
! 	对x积分     
    do 40 ix=1,ngus
	    xi=ga(k+ix-1)
! 	对y积分
		do 40 iy=1,ngus
            temp=0
            eta=ga(k+iy-1)
! 	调用isopar，求得形状函数J及单位法向矢量等     
			call isopar_8(xi,eta,corda,sn,snx,sne,aj,sr,dn)
            
            !WRITE(*,209)SR,AJ
            209FORMAT(3F15.6,E15.8) 

			do 60 j=1,8
				temp(1)=temp(1)+sn(j)
                !temp(2)=temp(2)+ABS(tcorda(j,2))*sn(j)
                !temp(3)=temp(3)+corda(j,1)*ABS(tcorda(j,2))*sn(j)
                !temp(4)=temp(4)+corda(j,3)*ABS(tcorda(j,2))*sn(j)
                
                temp(2)=temp(2)+tcorda(j,3)*dn(3)*sn(j)
                temp(3)=temp(3)-corda(j,1)*tcorda(j,3)*dn(3)*sn(j)
                temp(4)=temp(4)-corda(j,3)*tcorda(j,2)*dn(2)*sn(j)
                

                temp(5)=temp(5)+corda(j,2)**2*sn(j)
                temp(6)=temp(6)+corda(j,1)**2*sn(j)
                !temp(7)=temp(7)+corda(j,1)*corda(j,2)*sn(j)
                temp(7)=temp(7)+corda(j,1)*sn(j)
                temp(8)=temp(8)+sn(j)
    60		continue

            temp1(1)=temp(1)*gw(k+iy-1)*aj
            temp1(2)=temp(2)*gw(k+iy-1)*aj
            temp1(3)=temp(3)*gw(k+iy-1)*aj
            temp1(4)=temp(4)*gw(k+iy-1)*aj
!c
            temp1(5)=temp(5)*gw(k+iy-1)*aj
            temp1(6)=temp(6)*gw(k+iy-1)*aj
            temp1(7)=temp(7)*gw(k+iy-1)*aj

			!do 40 j=1,4
			twfa=twfa+temp1(1)*gw(k+ix-1)
			tvol(1,1)=tvol(1,1)+temp1(2)*gw(k+ix-1)
            tvol(1,2)=tvol(1,2)+temp1(3)*gw(k+ix-1)
            tvol(3,2)=tvol(3,2)+temp1(4)*gw(k+ix-1)

			TAMON(1)=TAMON(1)+temp1(5)*gw(k+ix-1)
            TAMON(2)=TAMON(2)+temp1(6)*gw(k+ix-1)
            TAMON(3)=TAMON(3)+temp1(7)*gw(k+ix-1)
    40continue

      return
END SUBROUTINE

SUBROUTINE BVOLUM_9(XGRA,M,CORDA,TWFA,TVOL,TAMON)
!c    ***********************************************      
!c      calculate the WFA in the matrix     
!c      corda(9,3) the cord. of 9 nodes (clockwise)      
!c      conm(9,3) the normal vector of 9 nodes         
!c      ps(3) the cord. of field point                 
!c      aij(9): panel source effect                     
!c      hij(9): panel dipole effect
!c   ************************************************
    USE NMLMOD
    REAL:: XI,ETA,AJ,TEMP(9),TEMP1(9),TWFA,TVOL(3,6),TAMON(3)
    REAL:: corda(9,3),tcorda(9,3),ps(3),conm(9,3),pr(3)
    REAL:: r(9),srr(3),ay(9),hy(9)
    REAL:: sn(9),snx(9),sne(9),sr(3),dn(3),XGRA(3)
    COMMON /gwa/ga(400),gw(400)
    CHARACTER M

!c	初始化每一面元的四个顶点的源效应aij（4）和偶极子效应hij（4）
    IF(M.EQ.'V') THEN
        tcorda=corda
        !corda(:,2)=0
        !WRITE(*,*)"corda",corda(:,3)
    END IF 

    twfa=0
    tvol=0
    temp=0
    temp1=0
    TAMON=0
    ngus=21

!c	k=2+3+...+(n-1)+1=n*(n-1)/2
    k=ngus*(ngus-1)/2
!c  gauss-legendre numerical integration
    
!c	对x积分     
    do 40 ix=1,ngus
	    xi=ga(k+ix-1)
!c	对y积分
		do 40 iy=1,ngus
            temp=0
            eta=ga(k+iy-1)
!c	调用isopar，求得形状函数J及单位法向矢量等     
			call isopar_9(xi,eta,corda,sn,snx,sne,aj,sr,dn)
      
			!do 60 j=1,9
			!	temp(1)=temp(1)+sn(j)
            !    !temp(2)=temp(2)+ABS(tcorda(j,2))*sn(j)
            !    !temp(3)=temp(3)+corda(j,1)*ABS(tcorda(j,2))*sn(j)
            !    !temp(4)=temp(4)+corda(j,3)*ABS(tcorda(j,2))*sn(j)
            !    
            !    temp(2)=temp(2)+tcorda(j,3)*dn(3)*sn(j)
            !    temp(3)=temp(3)-corda(j,1)*tcorda(j,3)*dn(3)*sn(j)
            !    temp(4)=temp(4)-corda(j,3)*tcorda(j,2)*dn(2)*sn(j)
            !    
            !
            !    temp(5)=temp(5)+corda(j,2)**2*sn(j)
            !    temp(6)=temp(6)+corda(j,1)**2*sn(j)
            !    !temp(7)=temp(7)+corda(j,1)*corda(j,2)*sn(j)
            !    temp(7)=temp(7)+corda(j,1)*sn(j)
            !    temp(8)=temp(8)+sn(j)
!    60		continue

            do 60 j=1,9
		        temp(1)=1.
                
                !temp(2)=temp(2)+tcorda(j,3)*dn(3)*sn(j)
                temp(2)=SR(3)*DN(3)
                !temp(2)=SR(1)*DN(1)
                temp(3)=-(SR(1))*SR(3)*DN(3)
                temp(4)=-SR(3)*SR(2)*DN(2)
                !temp(3)=temp(3)-corda(j,1)*tcorda(j,3)*dn(3)*sn(j)
                !temp(4)=temp(4)-corda(j,3)*tcorda(j,2)*dn(2)*sn(j)
                
                temp(5)=SR(2)**2
                temp(6)=(SR(1))**2
                temp(7)=(SR(1))
                temp(8)=1.

                temp(9)=SR(3)*DN(1)
                !WRITE(*,*)XGRA(1)
                !temp(5)=temp(5)+corda(j,2)**2*sn(j)
                !temp(6)=temp(6)+corda(j,1)**2*sn(j)
                !temp(7)=temp(7)+corda(j,1)*sn(j)
                !temp(8)=temp(8)+sn(j)
    60		continue

            temp1(1)=temp(1)*gw(k+iy-1)*aj
            temp1(2)=temp(2)*gw(k+iy-1)*aj
            temp1(3)=temp(3)*gw(k+iy-1)*aj
            temp1(4)=temp(4)*gw(k+iy-1)*aj
!c
            temp1(5)=temp(5)*gw(k+iy-1)*aj
            temp1(6)=temp(6)*gw(k+iy-1)*aj
            temp1(7)=temp(7)*gw(k+iy-1)*aj
            temp1(9)=temp(9)*gw(k+iy-1)*aj

			!do 40 j=1,4
			twfa=twfa+temp1(1)*gw(k+ix-1)
			tvol(1,1)=tvol(1,1)+temp1(2)*gw(k+ix-1)
            tvol(1,2)=tvol(1,2)+temp1(3)*gw(k+ix-1)
            tvol(3,2)=tvol(3,2)+temp1(4)*gw(k+ix-1)

            tvol(1,3)=tvol(1,3)+temp1(9)*gw(k+ix-1)

			TAMON(1)=TAMON(1)+temp1(5)*gw(k+ix-1)
            TAMON(2)=TAMON(2)+temp1(6)*gw(k+ix-1)
            TAMON(3)=TAMON(3)+temp1(7)*gw(k+ix-1)
    40continue

      return
END SUBROUTINE

SUBROUTINE BVOLUM_8V(M,CORDA,CORDAO,TWFA)
!c    ***********************************************      
!c     calculate the aij and hij terms in the matrix   
!c      the ps is one of the nodes of element          4---7---3
!c      cord(8,3) the cord. of 8 nodes (clockwise)     |       |
!c      conm(8,3) the normal vector of 8 nodes         8       6
!c      ps(3) the cord. of fiele point                 |       |
!c      aij(8): panel source effect                    1---5---2
!c      hij(8): panel dipole effect
!c   ************************************************
!cdp   implicit REAL(a-h,o-z)
    USE NMLMOD
    REAL:: XI,ETA,TEMP(8),TEMP1(8),TWFA,TVOL(3,6),TAMON(3)
    REAL:: corda(8,3),tcorda(8,3),ps(3),conm(8,3),pr(3),cordaO(9,3)
    REAL:: r(8),srr(3),ay(8),hy(8),p1,p2
    COMMON /gwa/ga(400),gw(400)
    CHARACTER M

      REAL:: XI1,XI2,XI3,XI4,ETA1,ETA2,ETA3,ETA4
      REAL:: aij1(8),hij1(8),cord1(8,3),conm1(8,3)
      REAL:: aj(4)
      REAL:: snx(8),sne(8),sm(8,4)
      REAL:: sn1(8),sn2(8),sn3(8),sn(8,4),sn4(8)
      equivalence (sn(1,1),sn1),(sn(1,2),sn2),(sn(1,3),sn3),(sn(1,4),sn4)
      REAL:: sr1(3),sr2(3),sr3(3),sr(3,4),sr4(3)
      equivalence (sr(1,1),sr1),(sr(1,2),sr2),(sr(1,3),sr3),(sr(1,4),sr4)
      REAL:: dn1(3),dn2(3),dn3(3),dn(3,4),dn4(3)
      equivalence (dn(1,1),dn1),(dn(1,2),dn2),(dn(1,3),dn3),(dn(1,4),dn4)
      common nbp,nfp,nbe,nfe,np,mp,nf,radius
      REAL:: CORDT(3,2)
      REAL:: SNT(3,4),SNXT(3),SNET(3),AJT(4),SRT(3,4),DNT(3,4)

      !write(*,*)"test"

      D1=0.9
      D2=D1
      D3=D2
      D4=D3

      r1=0.
      isn=1

!	初始化每一面元的四个顶点的源效应aij（4）和偶极子效应hij（4）

    twfa=0
    tvol=0
    temp=0
    temp1=0
    TAMON=0
    ngus=21


! 	k=2+3+...+(n-1)+1=n*(n-1)/2
    KGUS=NGUS*(NGUS-1)/2
!   gauss-legendre numerical integration
    
         
	SELECT CASE(ISN)
    CASE(1)
        DXI=-D1;DET=-D3
    CASE(2)
        DXI=D2;DET=-D3
    CASE(3)
        DXI=D2;DET=D4
    CASE(4)
        DXI=-D1;DET=D4
    CASE(5)
        DXI=0.;DET=-D3
    CASE(6)
        DXI=D2;DET=0.
    CASE(7)
        DXI=0.;DET=D4
    CASE(8) 
        DXI=-D1;DET=0.
    END SELECT       
           
    DO IX=1,NGUS
        XI=GA(KGUS+IX-1)
        P1=0.5*(1.0+GA(KGUS+IX-1))
        DO IY=1,NGUS
            ETA=GA(KGUS+IY-1)
            P2=0.5*(1.0+GA(KGUS+IY-1))
            temp=0.
! 	调用isopar，求得形状函数J及单位法向矢量等    

            CORDT(1,1)=DXI;CORDT(2,1)=-1.;CORDT(3,1)=-1.
            CORDT(1,2)=DET;CORDT(2,2)=1.;CORDT(3,2)=-1.
            !CALL isopar_3(xi,eta,cordT,snT(:,1),snxT,sneT,ajT(1),srT(:,1),dnT(:,1))
            !WRITE(*,*)AJT(1)
            CORDT(1,1)=DXI;CORDT(2,1)=-1.;CORDT(3,1)=1.
            CORDT(1,2)=DET;CORDT(2,2)=-1.;CORDT(3,2)=-1.
            !CALL isopar_3(xi,eta,cordT,snT(:,2),snxT,sneT,ajT(2),srT(:,2),dnT(:,2))
            CORDT(1,1)=DXI;CORDT(2,1)=1.;CORDT(3,1)=1.
            CORDT(1,2)=DET;CORDT(2,2)=-1.;CORDT(3,2)=1.
            !CALL isopar_3(xi,eta,cordT,snT(:,3),snxT,sneT,ajT(3),srT(:,3),dnT(:,3))
            CORDT(1,1)=DXI;CORDT(2,1)=1.;CORDT(3,1)=-1.
            CORDT(1,2)=DET;CORDT(2,2)=1.;CORDT(3,2)=1.
            !CALL isopar_3(xi,eta,cordT,snT(:,4),snxT,sneT,ajT(4),srT(:,4),dnT(:,4))

            XI1=(1.-P1)*(DXI)+P1*(1-P2)*(-1.)+P1*P2*(-1.)
            ETA1=(1.-P1)*(DET)+P1*(1-P2)*(1.)+P1*P2*(-1.)
            XI2=(1.-P1)*(DXI)+P1*(1-P2)*(-1.)+P1*P2*(1.)
            ETA2=(1.-P1)*(DET)+P1*(1-P2)*(-1.)+P1*P2*(-1.)
            XI3=(1.-P1)*(DXI)+P1*(1-P2)*(1.)+P1*P2*(1.)
            ETA3=(1.-P1)*(DET)+P1*(1-P2)*(-1.)+P1*P2*(1.)
            XI4=(1.-P1)*(DXI)+P1*(1-P2)*(1.)+P1*P2*(-1.)
            ETA4=(1.-P1)*(DET)+P1*(1-P2)*(1.)+P1*P2*(1.)

            CALL ISOPAR_8(XI1,ETA1,CORDA,SN1,SNX,SNE,AJ(1),SR1,DN(:,1))
            CALL ISOPAR_8(XI2,ETA2,CORDA,SN2,SNX,SNE,AJ(2),SR2,DN(:,2))
            CALL ISOPAR_8(XI3,ETA3,CORDA,SN3,SNX,SNE,AJ(3),SR3,DN(:,3))
	        CALL ISOPAR_8(XI4,ETA4,CORDA,SN4,SNX,SNE,AJ(4),SR4,DN(:,4)) 
 
            !WRITE(*,209)SR,AJ
            209FORMAT(3F15.6,E15.8) 

			!do j=1,8
                DO k=1,4

                    IF(K.EQ.1) AM=ABS(-1.-DXI)
                    IF(K.EQ.2) AM=ABS(-1.-DET)
                    IF(K.EQ.3) AM=ABS(1.-DXI)
                    IF(K.EQ.4) AM=ABS(1.-DET)

                    !WRITE(*,*)AM/2.*P1,AJT(K)

                    temp(1)=temp(1)+aj(k)*AJT(K)
                    !temp(1)=aj(1)*AJT(1)
                    !temp(1)=temp(1)+aj(4)**AM/2.*P1 !AJT(4)

                END DO
            !END DO

            temp1(1)=temp(1)*gw(kGUS+iy-1)
			twfa=twfa+temp1(1)*gw(kGUS+ix-1)

        end do
    end do

      return
END SUBROUTINE


SUBROUTINE BVOLUM_8V_OLD(M,CORDA,CORDAO,TWFA)
!c    ***********************************************      
!c      calculate the WFA in the matrix     
!c      corda(8,3) the cord. of 8 nodes (clockwise)      
!c      conm(8,3) the normal vector of 9 nodes         
!c      ps(3) the cord. of field point                 
!c      aij(8): panel source effect                     
!c      hij(8): panel dipole effect
!c   ************************************************
    USE NMLMOD
    REAL:: XI,ETA,AJ,TEMP(7),TEMP1(7),TWFA,TVOL(3,6),TAMON(3),cordaO(9,3)
    REAL:: corda(8,3),tcorda(8,3),ps(3),conm(8,3),pr(3)
    REAL:: r(8),srr(3),ay(8),hy(8)
    REAL:: sn(8),snx(8),sne(8),sr(3),dn(3)
    COMMON /gwa/ga(400),gw(400)
    CHARACTER M

!	初始化每一面元的四个顶点的源效应aij（4）和偶极子效应hij（4）

    twfa=0
    tvol=0
    temp=0
    temp1=0
    TAMON=0
    ngus=11

!c	k=2+3+...+(n-1)+1=n*(n-1)/2
    k=ngus*(ngus-1)/2
!c  gauss-legendre numerical integration
    
!c	对x积分     
    do 40 ix=1,ngus
	    xi=ga(k+ix-1)
!c	对y积分
		do 40 iy=1,ngus
            temp=0
            eta=ga(k+iy-1)
!c	调用isopar，求得形状函数J及单位法向矢量等     
			call isopar_8(xi,eta,corda,sn,snx,sne,aj,sr,dn)
			!call isopar_9(xi,eta,cordaO,sn,snx,sne,aj,sr,dn)
      
			do 60 j=1,8
				temp(1)=aj !temp(1)+sn(j)
    60		continue

            temp1(1)=temp(1)*gw(k+iy-1)
			twfa=twfa+temp1(1)*gw(k+ix-1)
    40continue

      return
END SUBROUTINE


SUBROUTINE BVOLUM_9V(M,CORDA,TWFA)
!c    ***********************************************      
!c      calculate the WFA in the matrix     
!c      corda(9,3) the cord. of 9 nodes (clockwise)      
!c      conm(9,3) the normal vector of 9 nodes         
!c      ps(3) the cord. of field point                 
!c      aij(9): panel source effect                     
!c      hij(9): panel dipole effect
!c   ************************************************
    USE NMLMOD
    REAL:: XI,ETA,AJ,TEMP(7),TEMP1(7),TWFA,TVOL(3,6),TAMON(3)
    REAL:: corda(9,3),tcorda(9,3),ps(3),conm(9,3),pr(3)
    REAL:: r(9),srr(3),ay(9),hy(9)
    REAL:: sn(9),snx(9),sne(9),sr(3),dn(3)
    COMMON /gwa/ga(400),gw(400)
    CHARACTER M

!	初始化每一面元的四个顶点的源效应aij（4）和偶极子效应hij（4）

    twfa=0
    tvol=0
    temp=0
    temp1=0
    TAMON=0
    ngus=11

!c	k=2+3+...+(n-1)+1=n*(n-1)/2
    k=ngus*(ngus-1)/2
!c  gauss-legendre numerical integration
    
!c	对x积分     
    do 40 ix=1,ngus
	    xi=ga(k+ix-1)
!c	对y积分
		do 40 iy=1,ngus
            temp=0
            eta=ga(k+iy-1)
!c	调用isopar，求得形状函数J及单位法向矢量等     
			call isopar_9(xi,eta,corda,sn,snx,sne,aj,sr,dn)
      
			do 60 j=1,9
				temp(1)=1. !temp(1)+sn(j)
    60		continue

            temp1(1)=temp(1)*gw(k+iy-1)*aj
			twfa=twfa+temp1(1)*gw(k+ix-1)
    40continue

      return
END SUBROUTINE

