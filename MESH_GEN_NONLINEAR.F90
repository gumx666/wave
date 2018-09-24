SUBROUTINE OFFSHIP_NOL
    !船体和自由面的水线进行统一编号
    !需要对水线上的点提供物面的法向矢量，如果是对自由面
    !上的面元,在计算时另外提供法向矢量


    !GENERATE THE OFFSET DATA OF THE VERTEXES AND THE ELEMENTS
    !NTPN(NUMBER OF POINTS),NTP(NUMBER OF PANELS)
	USE GREENMOD
	USE CUMOD

    REAL:: FAC,FACF,FACB,RATIO,LB,LF,LA,WX,WY,WZ,WXX
    REAL ::TMX,TMY,TMZ,TVX,TVY,TVZ,THETA,tppn
    INTEGER :: INEB,COUT

	DIMENSION X(5000,5000),Y(5000,5000),Z(5000,5000) 
    DIMENSION XS(25000),YS(25000),ZS(25000),INEB(25000,9)    
    DIMENSION WX(-200:500),WY(-200:500),WZ(5000),WXX(500)
    REAL WX1(200),WY1(200),WZ1(200),WX2(200,3),UI(200)                     
	REAL VX(500,500),VY(500,500),VZ(500,500)
    REAL XTR1(100),ZTR1(100),TAOX1(100)
    INTEGER NWP
	DIMENSION AX(500),AY(500,500)
    CHARACTER*7:: TITLE,TITLE2
    REAL:: WPVI1(200,3),WPVI2(200,3)
    INTEGER:: NPVI1,NPVI2
    INTEGER,DIMENSION(:):: MARKWAP1(100)
    INTEGER:: NWAPL1
	ALLOCATE(FX(5000,5000),FY(5000,5000))
    
    !PANELINDEX.DAT:
    !FREEINDEX.DAT: POINTS AND INDEX OF FREESURFACE
    !ALLINDEX.DAT: GENERATE ALL THE INDEX AND POINTS ON TWO FACES
    !GENERATE THE INTERFACE DATA FOR SECOND LAYER INTEGRATION
    !用于offsetin子程序的输入文件                                                                                                                                                                                                                                                   
	OPEN(11,FILE='PANELINDEX.DAT')
	OPEN(12,FILE='FREEINDEX.DAT')
	OPEN(13,FILE='ALLINDEX.DAT')

    !用于形成界面和物面及总体的网格图象
    !PANELPIC.DAT: DRAW THE PICTURE OF THE PANELS ON THE EPISODE
    !ALLPIC.DAT: DRAW THE OVERALL PICTURE OF THE PANELS 
	OPEN(14,FILE='PANELPIC.DAT')
	OPEN(15,FILE='FREEPIC.DAT')
	OPEN(16,FILE='ALLPIC.DAT')

    !用于形成物面法向矢量的分布图，辨别矢量方向的正确与否
	!OPEN(17,FILE='PANELVECTOR.DAT')
	!OPEN(18,FILE='FREEVECTOR.DAT')
	!OPEN(19,FILE='ALLVECTOR.DAT')
	!WRITE(*,*) 
	!WRITE(*,*)'BEGIN SHIP MESH GENERATE 9NODES'

    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(17,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\PANELVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    !OPEN(18,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREEVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    !OPEN(19,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\ALLVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(17,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\PANELVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(18,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREEVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(19,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\ALLVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    

	IF(ABS(RL/(NX-1)-DFX).GT.1E-6)THEN
		!WRITE(*,*)'ERROR	NX?	
        !STOP
	ENDIF

    !NFX=NX
    THETA=30.
    THETAF=18.

    FAC=2.0
    FACF=1.0
    !FACB=2.0
    FACB=1.0
    
    RATIO=0.5  !linear ok

    LS=RL
    LF=1.0*DARL !1.4*DARL
    LA=0.24*DARL !0.6*DARL
    LB=DARB !RATIO*DARL

    DPX=LA
    DPY=LB

    !界面上的网格选取多大还需要继续考虑
    !此参数作为全局变量，在greenmod中定义，

    !物面总的节点数	
    !NBPOINT=NX*NZ !+NX*(NZ-1)
    !物面总的面元数
	!NBBLOCK=(NX-1)*(NZ-1)/4 !*2
    
    !OPEN(10,FILE='S60.NEU')
    !OPEN(10,FILE='yuzhengchuan.NEU')
    !OPEN(10,FILE='KCS.NEU')
    !OPEN(10,FILE='KCS1.NEU')
    !OPEN(10,FILE='RS9N.neu')
    IF(ITTE.EQ.1) THEN
    !IF(ITTE.EQ.1.AND.NIT.EQ.1) THEN
    NWP=50
    WZ=0
    !WZ=0.01
    WPLZ=0
    DO I=1,NWP
        WX(I)=0.5-(I-1)*1./(NWP-1)
    END DO
    ELSE
    
    NWP=NFX
    DO I=1,NWP
        IP=NBPOINT+(NFXF+I-1)*NFY+1
        WX(I)=PCORD(IP,1)
        WZ(I)=WPLZ(I+NFXF)
        !IF(ITTE.LE.6) WZ(I)=WPLZ(I)*(0.6+0.4/7.*ITTE) 
        !IF(MDFREE.EQ.1) 
        !IF(NSWIN(13,1)/PANIN(1,3).LE.1.8) WZ(I)=WPLAZ(I+NFXF)
        IF(MDFREE.EQ.3) WZ(I)=0.
    END DO

    DO I=1,NWP
        WX1(I)=WX(NWP-I+1);WZ1(I)=WZ(NWP-I+1)
    END DO
    DO I=1,NWP
        WX(I)=WX1(I);WZ(I)=WZ1(I)
        !WRITE(*,*)WZ(I)
    END DO
    END IF

    IF(ITTE.EQ.1) THEN
    SINK=0.
    TRI=0.
    END IF

    IF(MDFREE.EQ.2) THEN
    SINK=0.
    TRI=0.
    END IF

    !WRITE(*,*)"NFX ************",NFX
    IF(XMDEN.NE.2) THEN
    !DFX=LB/40.
    ZFAC=0.75*DFX
    ZFAC1=ZFAC
    !WRITE(*,*)"ZFAC          ",ZFAC
    CALL SHIPHULL(SINK,TRI,NWP,WX(1:NWP),WZ(1:NWP),NFX1,WX1,WY1,XTR1,ZTR1,TAOX1,NTAO,ZFAC1,WPVI1,NPVI1,WPVI2,NPVI2,NZ1,NWAPL1,MARKWAP1)
    NWAPL=NWAPL1
    MARKWAP=MARKWAP1
    XTR=XTR1;ZTR=ZTR1;TAOX=TAOX1
    NBZ=NZ1
    END IF
        
    
    IF(XMDEN.EQ.1) THEN
    OPEN(10,FILE='SHIPHULL_MESH.DAT')
    END IF
    !WRITE(*,*)"NFX.............",NFX1

    IF(XMDEN.EQ.2) THEN
    OPEN(10,FILE='PANEL.NEU')
    DO I=1,6
        READ(10,*)
    END DO
    !读入船体面元参数
    READ(10,*)NBPOINT,NBBLOCK
    READ(10,*)
    READ(10,*)
    END IF
    
    
    !OPEN(1001,FILE='VIRTUALHULL.DAT')
    !WRITE(1001,*)NPVI1
    DO I=1,NPVI1
    !    WRITE(1001,*)WPVI1(I,:)
    END DO
    !WRITE(1001,*)NPVI2
    DO I=1,NPVI2
    !    WRITE(1001,*)WPVI2(I,:)
    END DO
    CLOSE(1001)

    
    !DO I=1,NFX
        !WRITE(*,*)I,WX1(I),WY1(I)
    !END DO
    
    IF(XMDEN.EQ.1) THEN
    !读入船体面元参数
    NBPOINT1=NBPOINT
    READ(10,*)NBPOINT,NBBLOCK
    IF(ITTE.EQ.1) NBPOINT1=NBPOINT
    !自由面总的节点数
	!NFPOINT=NFX*NFY
    !自由面总的面元数
	!NFBLOCK=(NFX-1)*(NFY-1)/4
    !控制面总的节点数
	!NRPOINT=0. !NZ*(2*NFY+NFX)
    !控制面总的面元数
	!NRBLOCK=0. !(NZ-1)*(2*NFY-2+NFX-1)/4
	
    !读入坐标数据
    DO I=1,NBPOINT
        READ(10,*)XS(I),YS(I),ZS(I)
        !WRITE(11,104)XS(I),YS(I),ZS(I),0,0,0
        !WRITE(13,104)XS(I),YS(I),ZS(I),0,0,0
	    !WRITE(17,104)XS(I),YS(I),ZS(I),0,0,0
        !WRITE(18,104)XS(I),YS(I),ZS(I),0,0,0
	    !WRITE(19,104)XS(I),YS(I),ZS(I),0,0,0
    END DO
    ELSE
        DO I=1,NBPOINT
        READ(10,*)TPPN,XS(I),YS(I),ZS(I)
        !WRITE(*,*)XS(I),YS(I),ZS(I)
        XS(I)=XS(I)/SCALE2
        YS(I)=YS(I)/SCALE2
        ZS(I)=ZS(I)/SCALE2
        END DO
    END IF

    WY=0.
    
    IF(XMDEN.EQ.1) THEN
    DO I=1,NBBLOCK
        READ(10,*)INEB(I,1:9)
    END DO
    ELSE
    READ(10,*)
    READ(10,*)
    DO I=1,NBBLOCK
        READ(10,*)TPPN,TPPN,TPPN,INEB(I,1:9)
    END DO
    END IF
    CLOSE(10)

    IF(XMDEN.EQ.2) THEN
        COUT=0
        DO I=1,NBPOINT
            IF(ABS(ZS(I)).LE.1.E-6) THEN
                COUT=COUT+1
                WZ(COUT)=ZS(I);WX1(COUT)=XS(I);WY1(COUT)=YS(I)
                !WRITE(*,*)WX(COUT)
            END IF
        END DO

        DO I=COUT-1,1,-1
            DO  J=1,I     
                IF(WX1(J).GT.WX1(J+1)) THEN
                    TPX=WX1(J);TPY=WY1(J);TPZ=WZ(J)
                    WX1(J)=WX1(J+1);WY1(J)=WY1(J+1);WZ(J)=WZ(J+1)
                    WX1(J+1)=TPX;WY1(J+1)=TPY;WZ(J+1)=TPZ
                END IF
            END DO
        END DO  
    END IF
    
    IF(XMDEN.EQ.1) THEN
        WPLX(1:NFX1)=WX1(1:NFX1);WPLY(1:NFX1)=WY1(1:NFX1)
    ELSE IF(XMDEN.EQ.2) THEN
        NFX1=COUT
        WPLX(1:NFX1)=WX1(1:NFX1);WPLY(1:NFX1)=WY1(1:NFX1)
        DO I=1,NFX1
            WRITE(*,*)I,WX1(I),COUT
        end do
    END IF
    
    IF(ITTE.EQ.1) THEN
    GOTO 199

    !寻找水线节点
    COUT=0
    DO I=1,NBPOINT
        IF(ABS(ZS(I)).LE.1.E-6) THEN
            COUT=COUT+1
            WZ(COUT)=ZS(I);WX(COUT)=XS(I);WY(COUT)=YS(I)
            !WRITE(*,*)WX(COUT)
        END IF
    END DO

    DO I=COUT-1,1,-1
        DO  J=1,I     
            IF(WX(J).GT.WX(J+1)) THEN
                TPX=WX(J);TPY=WY(J);TPZ=WZ(J)
                WX(J)=WX(J+1);WY(J)=WY(J+1);WZ(J)=WZ(J+1)
                WX(J+1)=TPX;WY(J+1)=TPY;WZ(J+1)=TPZ
                !WRITE(*,*) WX(I)
            END IF
        END DO
    END DO
    !GOTO 199

!   光滑水线曲率
    DO I=NFX/2,1,-1
        IF(ABS((WY1(I+1)-WY(I))/(WX1(I+1)-WX1(I))).GE.TAN(THETA/180.*PI)) THEN
            WY1(I)=WY1(I+1)-(WY1(I+2)-WY1(I+1))
            WX1(I)=WX1(I+1)-(WX1(I+2)-WX1(I+1))
        END IF
        WRITE(20,121)I,WX1(I),WY1(I),WZ1(I)
    END DO

    OPEN(20,FILE='WATERLINE.DAT')
    DO I=1,COUT
        WRITE(20,121)I,WX1(I),WY1(I),WZ1(I)
    END DO
    CLOSE(20)
    121FORMAT(I8,4F15.6)

    199 CONTINUE
    END IF

!   结束


    !设置水线节点
    !NFX=COUT

    !DFX=(WX(COUT)-WX(1))/(NFX-1)
    !DO I=1,COUT
    !    WXX(I)=WX(1)+(I-1)*DFX
    !END DO
    !调用样条函数
    !CALL ESPL2(WX(1:COUT),WY(1:COUT),COUT,WXX(1:NFX),NFX,WY(1:NFX))
    !WX(1:NFX)=WXX(1:NFX)

    !NFX=(NFXF/1.8)
    !NFX=(NFXF/2.0)
    
    !NFX=INT(NFY/LB)

    !NFX=NFY
    !NFX=(NFXF/1.6)
    !NFX=(NFXF/2.0)

    !NFX=NFX1
    IF(MOD(NFX,2).EQ.0) NFX=NFX+1
    !WRITE(*,*)"NFX.............",NFX
    !IF(ABS(WY1(1)).GE.0.025) THEN
    IF(ABS(WY1(1)).GE.10) THEN
        MTROM=1    !加方尾
        NTFY=2
        NFPOINT=(NFX+NFXF+NFXA)*NFY+(NFXF+1)*NTFY
        NFBLOCK=(NFY-1)*(NFX+NFXF+NFXA-1)/4+(NFXF)*NTFY/4
    ELSE
        NFPOINT=(NFX+NFXF+NFXA)*NFY
        NFBLOCK=(NFY-1)*(NFX+NFXF+NFXA-1)/4
    END IF
    
    NRPOINT=0. !(NFX+NFXF+NFXA)*NZ+2*NFY*NZ
    NRBLOCK=0. !(NFX+NFXF+NFXA-1)*(NZ-1)/4+2*(NFY-1)*(NZ-1)/4

    !总节点数＝物面总节点数＋自由面总节点数－水线节点数
	NAPOINT=NBPOINT+NFPOINT+NRPOINT
    !总面元数＝物面总面元数＋自由面总面元数
	NABLOCK=NBBLOCK+NFBLOCK+NRBLOCK

102 format(1x,'title="panel on episode solid"'/1x,'variables="x","Y",& 
     "Z"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')
101	format(1x,'title="episode solid mesh"'/1x,'variables="x","Y",& 
     "Z","NX","NY","NZ"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')
   
    !WRITE(*,*)"NNODES  =    ",NFPOINT+NBPOINT
    !WRITE(*,*)"NPANELS =    ",NFBLOCK+NBBLOCK

	WRITE(11,*)NBPOINT,NBBLOCK*4
	WRITE(12,*)NFPOINT,NFBLOCK
	WRITE(13,*)NAPOINT,NABLOCK

	WRITE(14,102)NBPOINT,NBBLOCK*4
	WRITE(15,102)NFPOINT,NFBLOCK*4
	WRITE(16,102)NRPOINT,NRBLOCK*4

	WRITE(17,101)NBPOINT,NBBLOCK*4
	WRITE(18,101)NFPOINT+NBPOINT,(NFBLOCK+NBBLOCK)*4
	WRITE(19,101)NAPOINT,NABLOCK*4

    IF(ITTE.EQ.1..AND.NIT.EQ.1) THEN
    ALLOCATE(ETAW(NBPOINT+2*NFPOINT+500),ETAW0(NBPOINT+2*NFPOINT+500),ETAW1(NBPOINT+2*NFPOINT+500))
    ALLOCATE(SE(NBPOINT+2*NFPOINT+500))
    ALLOCATE(CORD0(NBPOINT+500,3),VECN0(NBPOINT+500,3))
    ALLOCATE(QG(NBPOINT+2*NFPOINT+500))
    ALLOCATE(DCORZ(NBPOINT+500))
    ETAW=0.
    ETAW0=0.
    ETAW1=0.
    ETAW0=0.
    WPLAZ=0.
    DCORZ=0. 
    END IF
    
    DO I=1,NBPOINT
        WRITE(11,104)XS(I),YS(I),ZS(I),0,0,0
        WRITE(13,104)XS(I),YS(I),ZS(I),0,0,0
	    WRITE(17,104)XS(I),YS(I),ZS(I),0,0,0
        WRITE(18,104)XS(I),YS(I),ZS(I),0,0,0
	    WRITE(19,104)XS(I),YS(I),ZS(I),0,0,0
    END DO

103	FORMAT(1X,I8,6F15.6)
104	FORMAT(1X,6F15.6)

    !GOTO 91

    FZ=0
    !DFX=LF/(NFXF)
    !DFY=LB/(NFY-1)
	!DFY=DFX
	K2=0	
	VX2=0.0
	VY2=0.0
	VZ2=1.0

    GOTO 299
    !WRITE(*,*)"DFY",LB,DFY
	DO I=1,NFXF
        WX1(I-NFXF)=WX1(1)-DFX*(NFXF-I+1)
            !WRITE(*,*)I
    END DO
    !光顺船尾后的水线
    DO I=1,-NFXF,-1
        IF(ABS((WY1(I+1)-WY(I))/(WX(I+1)-WX(I))).GE.TAN(THETA/180.*PI)) THEN
            WY(I)=WY(I+1)-(WY(I+2)-WY(I+1))
            !WRITE(*,*)I
        END IF
        !WRITE(20,121)I,WX(I),WY(I),WZ(I)
    END DO

    299 CONTINUE

!   AVERAGE WATERPLANE
    DO J=1,NFX
        UI(J)=1./(NFX-1)*(J-1)
    END DO
    
    WX2=0
    WX2(:,1)=WX1(:)
    WX2(:,2)=WY1(:)
    IF(MTROM.NE.1) THEN
    DO I1=NFX1/2,1,-1
        IF(ABS((WX2(I1+1,2)-WX2(I1,2))/(WX2(I1+1,1)-WX2(I1,1))).GE.TAN(THETA/180.*PI)) THEN
        !IF(ABS((WY1(I+1)-WY(I))/(WX(I+1)-WX(I))).GE.TAN(THETA/180.*PI)) THEN
            DO I=I1,1,-1
            WX2(I,2)=WX2(I+1,2)-SQRT((WX2(I+2,1)-WX2(I+1,1))**2+(WX2(I+2,2)-WX2(I+1,2))**2)*SIN(THETA/180.*PI)
            WX2(I,1)=WX2(I+1,1)-SQRT((WX2(I+2,1)-WX2(I+1,1))**2+(WX2(I+2,2)-WX2(I+1,2))**2)*COS(THETA/180.*PI)
            IF(WX2(I,2).LE.0) THEN
                WX2(I,2)=0.
                !THETA=0.
            END IF
            END DO
            EXIT
        END IF
    END DO
    ELSE
    DO I1=NFX1/2,1,-1
        IF(ABS((WX2(I1+1,2)-WX2(I1,2))/(WX2(I1+1,1)-WX2(I1,1))).GE.TAN(THETA/180.*PI)) THEN
        !IF(ABS((WY1(I+1)-WY(I))/(WX(I+1)-WX(I))).GE.TAN(THETA/180.*PI)) THEN
            DO I=I1,1,-1
            WX2(I,2)=WX2(I+1,2)
            !WX2(I,1)=WX2(I+1,1)-SQRT((WX2(I+2,1)-WX2(I+1,1))**2+(WX2(I+2,2)-WX2(I+1,2))**2)*COS(THETA/180.*PI)
            IF(WX2(I,2).LE.0) THEN
                WX2(I,2)=0.
                !THETA=0.
            END IF
            END DO
            EXIT
        END IF
    END DO
    END IF


    GOTO 807
!   圆形球鼻艏
    !WX2=0
    !WX2(:,1)=WX1(:)
    !WX2(:,2)=WY1(:)
    DO I1=NFX1/2+1,NFX1
        IF(ABS((WX2(I1+1,2)-WX2(I1,2))/(WX2(I1+1,1)-WX2(I1,1))).GE.TAN(THETAF/180.*PI)) THEN
        !IF(ABS((WY1(I+1)-WY(I))/(WX(I+1)-WX(I))).GE.TAN(THETA/180.*PI)) THEN
            DO I=I1,NFX1
            WX2(I,2)=WX2(I-1,2)-SQRT((WX2(I-2,1)-WX2(I-1,1))**2+(WX2(I-2,2)-WX2(I-1,2))**2)*SIN(THETA/180.*PI)
            WX2(I,1)=WX2(I-1,1)+SQRT((WX2(I-2,1)-WX2(I-1,1))**2+(WX2(I-2,2)-WX2(I-1,2))**2)*COS(THETA/180.*PI)
            IF(WX2(I,2).LE.0) THEN
                !WX2(I,2)=0.
                !THETA=0.
            END IF
            END DO
            EXIT
        END IF
    END DO   
    807 CONTINUE
!

    IF(MTROM.NE.1) THEN
    IF(WX2(1,2).GE.0.025) THEN
        WX2(1,1)=WX2(1,1)-0.08
        WX2(1,2)=WX2(1,2)*0.3
    END IF
    
    DO I=1,200
    IF(WX2(I,1).LE.-0.35) THEN
        WX2(I,1)=WX2(I,1) !*1.05
    END IF
    END DO
    ELSE
        WX2(I,1)=WX2(I,1)
    END IF

    !WPVI2,NPVI2
    CALL NURBS_INT(NFX1-1,WX2(1:NFX1,:),NFX-1,UI(1:NFX),WX2(1:NFX,:))
    HULL_WP(1:NFX,:)=WX2(1:NFX,:)
    IF(ITTE.EQ.1) HULL_WP0(1:NFX,:)=WX2(1:NFX,:)

    DFX=ABS(WX2(NFX/2,1)-WX2(NFX/2-1,1)) !*1.1 !*1.4 !*1.2    !kcs
    !DFX=1.0/50.

    IF(MTROM.EQ.1) THEN
        WX2(1,1)=WX2(1,1)-DFX
    END IF

    CALL NURBS_INT(NFX-1,WX2(1:NFX,:),NFX-1,UI(1:NFX),WX2(1:NFX,:))
    HULL_WP(1:NFX,:)=WX2(1:NFX,:)
    IF(ITTE.EQ.1) HULL_WP0(1:NFX,:)=WX2(1:NFX,:)
    DFX=ABS(WX2(NFX/2,1)-WX2(NFX/2-1,1)) !*1.1 !*1.4 !*1.2    !kcs    
    
    !IF(ITTE.GE.2) THEN
    !DO I=1,NPVI2
    !    WX2(I,:)=WPVI2(NPVI2-I+1,:)
    !END DO
    !CALL NURBS_INT(NPVI2-1,WX2(1:NPVI2,:),NFX-1,UI(1:NFX),WX2(1:NFX,:))
    !END IF
! 

    !DFX=ABS(X(1,1)-X(2,1))
    !DFX=ABS(WX(COUT/2)-WX(COUT/2-1))*1.0 !*1.4 !*1.2
    !DFX=LB/(NFY-1) !*1.1 !*1.4 !*1.2    !kcs
    
    !DFX=ABS(WX(COUT/2)-WX(COUT/2-1))*1.2 !*1.4 !*1.2    !YUZHENG


!FIRST PART OF THE FREE SURFACE (UPSTREAM, BEFORE THE BOW)
	DO I=1,NFXF
		DO J=1,NFY  
            
            !FX(I,J)=OUTX(NFX,1)+LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            
            !FX(I,J)=WX1(1)-LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            !FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            !FX(I,J)=X(1,1)-LF*REAL(NFXF-I+1)/(NFXF)

            IF(MTROM.NE.1) THEN
            NSJ=3
            IF(ABS(WX2(1,2)).GE.0.01) THEN
                IF(I.GE.NFXF-NSJ) THEN
                    FX(I,J)=WX2(1,1)-DFX*(NFXF-I+1)
                    FY(I,J)=WX2(1,2)-DFX*(NFXF-I+1)*TAN(THETAF/180.*PI)+(LB-WX2(1,2)+DFX*(NFXF-I+1)*TAN(THETAF/180.*PI))*REAL(J-1)/(NFY-1)
                    !WRITE(*,*)WX2(1,2),DFX*(NFXF-I+1)
                    IF(WX2(1,2).LE.DFX*(NFXF-I+1)*TAN(THETAF/180.*PI)) FY(I,J)=(LB)*REAL(J-1)/(NFY-1)
                ELSE
                    FX(I,J)=WX2(1,1)-DFX*(NFXF-I+1)
                    FY(I,J)=(LB)*REAL(J-1)/(NFY-1)
                    !FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
                END IF                
            ELSE 
                FX(I,J)=WX2(1,1)-DFX*(NFXF-I+1)
                FY(I,J)=(LB)*REAL(J-1)/(NFY-1)
                !FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            END IF
            
            ELSE
                FX(I,J)=WX2(1,1)-DFX*(NFXF-I+1)
                FY(I,J)=WX2(1,2)+(LB-WX2(1,2))*REAL(J-1)/(NFY-1)
                !FY(I,J)=WX2(1,2)+(LB-WX2(1,2))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
                !WRITE(*,*)FY(I,1)
            END IF


            !IF(MOD(I,2).EQ.0) THEN
            !    FX(I,J)=WX1(1)-ABS(WX1(1)-WX1(2))*(NFXF-I+1)-ABS(WX1(1)-WX1(2))*(NFXF-I)*0.2
            !ELSE
            !    FX(I,J)=WX1(1)-ABS(WX1(1)-WX1(2))*(NFXF-I+1)-ABS(WX1(1)-WX1(2))*(NFXF-I-1)*0.2
            !END IF

            !FY(I,J)=WY1(I-NFXF)+(LB-WY1(I-NFXF))*REAL(J-1)/(NFY-1)
            
            IP=(I-1)*NFY+J
            FZ=ETAW(NBPOINT1+IP)
            	
            VX2=0.0
		    VY2=0.0
			VZ2=1.0
      
            K1=K1+1

			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    GOTO 399
    DO I=1,NFXF,2 
        DO J=1,NFY
            IF(I.GT.1) THEN
                FX(I-1,J)=(FX(I-2,J)+FX(I,J))/2.
                FY(I-1,J)=(FY(I-2,J)+FY(I,J))/2.
            END IF
        END DO
    END DO
	DO J=1,NFY,2 
        DO I=1,NFXF 
            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
            END IF
        END DO
    END DO

    DO I=1,NFXF-1
		DO J=1,NFY
            K1=K1+1

            IP=(I-1)*NFY+J
            FZ=ETAW(NBPOINT1+IP)
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    DO J=1,NFY
            FX(1,J)=(FX(NFXF-1,J)+WX2(1,1))/2.
            !FX(1,J)=(FX(NFXF-1,J)+X(1,1))/2.
            
            FY(1,J)=FY(NFXF-1,J)
            !FY(1,J)=(LB)*REAL(J-1)/(NFY-1)

            IP=(NFXF-1)*NFY+J
            FZ=ETAW(NBPOINT1+IP)
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
    END DO
    399 CONTINUE

   

!SECOND PART OF THE FREE SURFACE (ALONG SHIP)
	DO I=1,NFX
		DO J=1,NFY !,2  
            !DFY=(LB-TY(NFX+NFXF-I+1))/(NFY-1)

            !FX(I,J)=OUTX(NFX-I+1,1)
            !FX(I,J)=X(I,1)
            
            !FY(I,J)=WY(I)+(LB-WY(I))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            FY(I,J)=WX2(I,2)+(LB-WX2(I,2))*REAL(J-1)/(NFY-1)
            !FY(I,J)=WX2(I,2)+(LB-WX2(I,2))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            !FY(I,J)=Y(I,1)+(LB-Y(I,2))*REAL(J-1)/(NFY-1)
            !FX(I,J)=X(I,1)
            FX(I,J)=WX2(I,1)
            
            !FX(I,J)=WX1(I)

            !FY(I,J)=Y(I,1)+(LB-Y(I,2))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            VX2=0.0
		    VY2=0.0
			VZ2=1.0

            IF(J.GT.1) THEN
                !FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                !FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
                !WRITE(13,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
			    !WRITE(18,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
			    !WRITE(19,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
            END IF
         
            IP=NFXF*NFY+(I-1)*NFY+J
            FZ=ETAW(NBPOINT1+IP)

            K1=K1+1
			!WRITE(12,103)K2,FX(I),FY(I,J),FZ,VX2,VY2,VZ2
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
            !END IF
        END DO
    END DO

    DFY=ABS(FY(1,1)-FY(1,2))
    DFY=LB/20.
    DFY=1.0/22.
    DFY=1.0/25.
    !WRITE(*,*)"DFY     ",DFY

!THIRD PART OF THE FREE SURFACE (UPSTREAM, AFTER THE STERN)
    !DFX=LA/(NFXA)
    !DFY=LB/(NFY-1)

	DO I=1,NFXA !,2
		DO J=1,NFY !,2  
            !FX(I,J)=TX(1)-(I-NFX)*DFX
			!FY(I,J)=(J-1)*DFY

            !FX(I,J)=OUTX(1,1)-LA*(EXP(REAL(I)/(NFXA)*FAC)-1.)/(EXP(FAC)-1.)
            
            !FX(I,J)=WX(NFX)+LA*(EXP(REAL(I)/(NFXA)*FACF)-1.)/(EXP(FACF)-1.)
            
            !FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*REAL(I)/(NFXA)

            !FX(I,J)=WX1(NFX)+DFX*REAL(I)
            
            FX(I,J)=WX2(NFX,1)+DFX*REAL(I)
            !FX(I,J)=WX2(NFX,1)+ABS(WX2(NFX,1)-WX2(NFX-1,1))*REAL(I)
            
            FY(I,J)=LB*REAL(J-1)/(NFY-1)
            !FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            !WRITE(*,*)WX1(NFX)
            !FY(I,J)=(J-1)*DFY

            VX2=0.0
		    VY2=0.0
			VZ2=1.0

            IP=(NFXF+NFX)*NFY+(I-1)*NFY+J
            !write(*,*)NBPOINT+IP
            FZ=ETAW(NBPOINT1+IP)

            K1=K1+1
			!WRITE(12,103)K2,FX(I),FY(I,J),FZ,VX2,VY2,VZ2
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
            !END IF
        END DO
    END DO

    GOTO 499
    DO I=2,NFXA,2 
	    DO J=1,NFY 
            IF(I.GT.2) THEN
                FX(I-1,J)=(FX(I-2,J)+FX(I,J))/2.
                FY(I-1,J)=(FY(I-2,J)+FY(I,J))/2.
            END IF
        END DO
    END DO
	DO J=1,NFY,2 
        DO I=1,NFXA
            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
            END IF
        END DO
    END DO

    DO J=1,NFY
            FX(1,J)=(FX(2,J)+WX(NFX))/2.
            FY(1,J)=FY(2,J)
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
    END DO

    DO I=2,NFXA
		DO J=1,NFY
            K1=K1+1
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2

			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO
    499 CONTINUE

!TRANSOM
    IF(MTROM.EQ.1) THEN
    DO I=1,NFXF+1
        DO J=1,NTFY
            FX(I,J)=WX1(1)-DFX*(NFXF-I+1)
            FY(I,J)=WY1(1)*REAL(J-1)/(NTFY)
            
            IP=(NFXF+NFX+NFXA)*NFY+(I-1)*NTFY+J
            FZ=ETAW(NBPOINT1+IP)
            	
            VX2=0.0
		    VY2=0.0
			VZ2=1.0
      
            K1=K1+1

			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO
    END IF
!

!THE PANEL INDEX ON THE SHIP SURFACE
 	K2=0
	DO I=1,NBBLOCK !NX-1,2
		!DO J=1,NZ-1,2
			K2=K2+1
		    ND1=J+(I-1)*NZ
			ND2=ND1+1
            ND3=ND2+1
            ND4=ND3+NZ
            ND5=ND3+2*NZ
            ND6=ND5-1
            ND7=ND6-1
            ND8=ND7-NZ
			ND9=ND8+1

            ND1=INEB(I,1)
            ND2=INEB(I,2)
            ND3=INEB(I,3)
            ND4=INEB(I,4)
            ND5=INEB(I,5)
            ND6=INEB(I,6)
            ND7=INEB(I,7)
            ND8=INEB(I,8)
            ND9=INEB(I,9)

			WRITE(13,105)9,ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
            !WRITE(13,105)9,ND1,ND8,ND7,ND6,ND5,ND4,ND3,ND2,ND9

			!WRITE(16,*)ND1,ND2,ND3,ND4
			WRITE(11,*)ND1,ND2,ND9,ND8
            WRITE(11,*)ND8,ND9,ND6,ND7
            WRITE(11,*)ND2,ND3,ND4,ND9
            WRITE(11,*)ND9,ND4,ND5,ND6
			WRITE(17,*)ND1,ND2,ND9,ND8
            WRITE(17,*)ND8,ND9,ND6,ND7
            WRITE(17,*)ND2,ND3,ND4,ND9
            WRITE(17,*)ND9,ND4,ND5,ND6
			WRITE(18,*)ND1,ND2,ND9,ND8
            WRITE(18,*)ND8,ND9,ND6,ND7
            WRITE(18,*)ND2,ND3,ND4,ND9
            WRITE(18,*)ND9,ND4,ND5,ND6
			WRITE(19,*)ND1,ND2,ND9,ND8
            WRITE(19,*)ND8,ND9,ND6,ND7
            WRITE(19,*)ND2,ND3,ND4,ND9
            WRITE(19,*)ND9,ND4,ND5,ND6
        !END DO
    END DO
    105 FORMAT(1X,10I8)

!----	THE PANEL INDEX ON THE FREE SURFACE

	DO I=1,NFX+NFXA+NFXF-1,2
		DO J=1,NFY-1,2
            K3=K3+1
		    ND1=NBPOINT+J+(I-1)*NFY
			ND2=ND1+NFY !NBPOINT+J+I*NFY
            ND3=ND2+NFY
			ND4=ND3+1 !NBPOINT+J+1+I*NFY
			ND5=ND4+1
            ND6=ND5-NFY !NBPOINT+J+1+(I-1)*NFY
            ND7=ND6-NFY
            ND8=ND7-1
            ND9=ND8+NFY

			!WRITE(12,*)K3,NI,NI+1,NI+NFX+1,NI+NFX-1
			!WRITE(15,*)NI,NI+1,NI+NFX+1,NI+NFX-1
			!WRITE(18,*)NI,NI+1,NI+NFX+1,NI+NFX-1

            WRITE(13,105)9,ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
			WRITE(18,*)ND1,ND2,ND9,ND8
            WRITE(18,*)ND8,ND9,ND6,ND7
            WRITE(18,*)ND2,ND3,ND4,ND9
            WRITE(18,*)ND9,ND4,ND5,ND6
			WRITE(19,*)ND1,ND2,ND9,ND8
            WRITE(19,*)ND8,ND9,ND6,ND7
            WRITE(19,*)ND2,ND3,ND4,ND9
            WRITE(19,*)ND9,ND4,ND5,ND6
        END DO
    END DO


!TRANSOM
    IF(MTROM.EQ.1) THEN
	DO I=1,NFXF-1,2
		DO J=1,NTFY-1,2
            IF(J.NE.NTFY-1) THEN
            K3=K3+1
		    ND1=NBPOINT+(NFX+NFXF+NFXA)*NFY+(I-1)*NTFY+J
			ND2=ND1+NTFY !NBPOINT+J+I*NFY
            ND3=ND2+NTFY
			ND4=ND3+1 !NBPOINT+J+1+I*NFY
			ND5=ND4+1
            ND6=ND5-NTFY !NBPOINT+J+1+(I-1)*NFY
            ND7=ND6-NTFY
            ND8=ND7-1
            ND9=ND8+NTFY
            WRITE(13,105)9,ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
			WRITE(18,*)ND1,ND2,ND9,ND8
            WRITE(18,*)ND8,ND9,ND6,ND7
            WRITE(18,*)ND2,ND3,ND4,ND9
            WRITE(18,*)ND9,ND4,ND5,ND6
			WRITE(19,*)ND1,ND2,ND9,ND8
            WRITE(19,*)ND8,ND9,ND6,ND7
            WRITE(19,*)ND2,ND3,ND4,ND9
            WRITE(19,*)ND9,ND4,ND5,ND6

            ELSE
            K3=K3+1
		    ND1=NBPOINT+(NFX+NFXF+NFXA)*NFY+(I-1)*NTFY+J
			ND2=ND1+NTFY !NBPOINT+J+I*NFY
            ND3=ND2+NTFY
			ND4=ND3+1 !NBPOINT+J+1+I*NFY
            ND9=ND4-NTFY
            ND8=ND9-NTFY

			ND7=NBPOINT+(I-1)*NFY+1
            ND6=ND7+NFY !NBPOINT+J+1+(I-1)*NFY
            ND5=ND6+NFY !NBPOINT+J+1+(I-1)*NFY
            WRITE(13,105)9,ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
			WRITE(18,*)ND1,ND2,ND9,ND8
            WRITE(18,*)ND8,ND9,ND6,ND7
            WRITE(18,*)ND2,ND3,ND4,ND9
            WRITE(18,*)ND9,ND4,ND5,ND6
			WRITE(19,*)ND1,ND2,ND9,ND8
            WRITE(19,*)ND8,ND9,ND6,ND7
            WRITE(19,*)ND2,ND3,ND4,ND9
            WRITE(19,*)ND9,ND4,ND5,ND6
            END IF
        END DO
    END DO

    END IF
!
    !NFX=NFX+NFXF+NFXA

    DEALLOCATE(FX,FY)
	CLOSE(19)
	CLOSE(18)
	CLOSE(17)
	CLOSE(16)
	CLOSE(15)
	CLOSE(14)
	CLOSE(13)
	CLOSE(12)
	CLOSE(11)
END SUBROUTINE

SUBROUTINE OFFWIGLEY_9NODES_NON
    !船体和自由面的水线进行统一编号
    !需要对水线上的点提供物面的法向矢量，如果是对自由面
    !上的面元,在计算时另外提供法向矢量


    !GENERATE THE OFFSET DATA OF THE VERTEXES AND THE ELEMENTS
    !NTPN(NUMBER OF POINTS),NTP(NUMBER OF PANELS)
	USE GREENMOD
	USE CUMOD

    REAL:: FAC,FACF,FACB,RATIO,LB,LF,LA
    REAL ::TMX,TMY,TMZ,TVX,TVY,TVZ
    CHARACTER*7:: TITLE,TITLE2

	DIMENSION X(500,500),Y(500,500),Z(500,500)                       
	DIMENSION VX(500,500),VY(500,500),VZ(500,500)
	DIMENSION AX(500),AY(500,500)
    REAL:: TPWPLAZ(500)
	ALLOCATE(FX(500,500),FY(500,500))

    !PANELINDEX.DAT:
    !FREEINDEX.DAT: POINTS AND INDEX OF FREESURFACE
    !ALLINDEX.DAT: GENERATE ALL THE INDEX AND POINTS ON TWO FACES
    !GENERATE THE INTERFACE DATA FOR SECOND LAYER INTEGRATION
    !用于offsetin子程序的输入文件                                                                                                                                                                                                                                                   
	OPEN(11,FILE='PANELINDEX.DAT')
	OPEN(12,FILE='FREEINDEX.DAT')
	OPEN(13,FILE='ALLINDEX.DAT')

    !用于形成界面和物面及总体的网格图象
    !PANELPIC.DAT: DRAW THE PICTURE OF THE PANELS ON THE EPISODE
    !ALLPIC.DAT: DRAW THE OVERALL PICTURE OF THE PANELS 
	OPEN(14,FILE='PANELPIC.DAT')
	OPEN(15,FILE='FREEPIC.DAT')
	OPEN(16,FILE='ALLPIC.DAT')

    !用于形成物面法向矢量的分布图，辨别矢量方向的正确与否
	!OPEN(17,FILE='PANELVECTOR.DAT')
	!OPEN(18,FILE='FREEVECTOR.DAT')
	IF(MMESH.EQ.2) OPEN(19,FILE='ALLVECTOR.DAT')
	!WRITE(*,*) 
	!WRITE(*,*)'BEGIN SHIP MESH GENERATE 9NODES'
    
    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(17,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\PANELVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    !OPEN(18,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREEVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    !IF(MMESH.EQ.1) OPEN(19,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\ALLVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')

    OPEN(17,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\PANELVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(18,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREEVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(19,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\ALLVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
      
    
	IF(ABS(RL/(NX-1)-DFX).GT.1E-6)THEN
		!WRITE(*,*)'ERROR	NX?	
        !STOP
	ENDIF



    IF(ITTE.EQ.1) THEN
    WPLZ=0
    END IF

    NFX=NX

    FAC=2.0
    FACF=1.3
    FACB=1.0
    RATIO=0.4  !linear ok
    
    LS=RL
    LF=1.0*DARL !1.4*DARL
    LA=0.24*DARL !0.6*DARL
    LB=DARB !RATIO*DARL

    DPX=LA
    DPY=LB

    
    IF(MMESH.NE.2) THEN
    NBPOINT1=NBPOINT
    IF(ITTE.EQ.1) NBPOINT1=NBPOINT
    END IF

    !界面上的网格选取多大还需要继续考虑
    !此参数作为全局变量，在greenmod中定义，
    !在interfacelattice和interfacewave子程序中可以直接引用 
    !物面总的节点数	
	NBPOINT=NX*NZ !+NX*(NZ-1)
    !物面总的面元数
	NBBLOCK=(NX-1)*(NZ-1)/4 !*2
    !自由面总的节点数
	NFPOINT=NFX*NFY
    !自由面总的面元数
	NFBLOCK=(NFX-1)*(NFY-1)/4
    !控制面总的节点数
	NRPOINT=0. !NZ*(2*NFY+NFX)
    !控制面总的面元数
	NRBLOCK=0. !(NZ-1)*(2*NFY-2+NFX-1)/4

    NFPOINT=(NFX+NFXF+NFXA)*NFY
    NFBLOCK=(NFY-1)*(NFX+NFXF+NFXA-1)/4
    
    IF(MMESH.EQ.2) GOTO 2002
    IF(ITTE.EQ.1.AND.NIT.EQ.1) THEN
    ALLOCATE(ETAW(NBPOINT+2*NFPOINT+500),ETAW0(NBPOINT+2*NFPOINT+500),ETAW1(NBPOINT+2*NFPOINT+500))
    ALLOCATE(SE(NBPOINT+2*NFPOINT+500))
    ALLOCATE(CORD0(NBPOINT+500,3),VECN0(NBPOINT+500,3))
    ALLOCATE(QG(NBPOINT+2*NFPOINT+500))
    ALLOCATE(DCORZ(NBPOINT+500))
    ETAW=0.
    ETAW0=0.
    ETAW1=0.
    ETAW0=0.
    WPLAZ=0. 
    DCORZ=0.
    END IF
    2002 CONTINUE

    IF(MDFREE.EQ.2) THEN
    SUAT=SUA
    DFFT=DFF
    DFF=RD
    SUA=0
    END IF

    IF(ITTE.LE.2) THEN
    DFF=RD
    SUA=0
    WPLAZ=0.
    WPLZ=0.
    END IF
    
    NRPOINT=0. !(NFX+NFXF+NFXA)*NZ+2*NFY*NZ
    NRBLOCK=0. !(NFX+NFXF+NFXA-1)*(NZ-1)/4+2*(NFY-1)*(NZ-1)/4

    !总节点数＝物面总节点数＋自由面总节点数－水线节点数
	NAPOINT=NBPOINT+NFPOINT+NRPOINT
    !总面元数＝物面总面元数＋自由面总面元数
	NABLOCK=NBBLOCK+NFBLOCK+NRBLOCK

	WRITE(11,*)NBPOINT,NBBLOCK*4
	WRITE(12,*)NFPOINT,NFBLOCK
	WRITE(13,*)NAPOINT,NABLOCK

	WRITE(14,102)NBPOINT,NBBLOCK*4
	WRITE(15,102)NFPOINT,NFBLOCK*4
	WRITE(16,102)NRPOINT,NRBLOCK*4

	WRITE(17,101)NBPOINT,NBBLOCK*4
	WRITE(18,101)NFPOINT+NBPOINT,(NFBLOCK+NBBLOCK)*4
	WRITE(19,101)NAPOINT,NABLOCK*4
102 format(1x,'title="panel on episode solid"'/1x,'variables="x","Y",& 
     "Z"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')
101	format(1x,'title="episode solid mesh"'/1x,'variables="x","Y",& 
     "Z","NX","NY","NZ"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')
   
    !WRITE(*,*)"NNODES  =    ",NFPOINT+NBPOINT
    !WRITE(*,*)"NPANELS =    ",NFBLOCK+NBBLOCK
	
    IF(MDFREE.EQ.3) WPLZ=0
    IF(MMESH.eq.2) THEN
        WPLZ=0.
    end if

    TPWPLAZ=WPLAZ
    !WPLAZ=0. 
    !WPLZ=0.
    
    NWPLAZ=0
	DO I=1,NX
		!TMX=-RL/2.0+RL/(NX-1)*(I-1)
        
        IF(I.LE.(NX-1)/2) THEN
            !TMX=-RL/2.0+RL/2.*(EXP(REAL(I-1)/((NX-1)/2.)*FAC)-1.)/(EXP(FAC)-1.)
            TMX=-RL/2.0+RL/2.*REAL(I-1)/((NX-1)/2.)
        ELSE IF(I.EQ.(NX-1)/2+1) THEN
            TMX=0.
        ELSE
            I1=NX-I+1
            !TMX=RL/2.0-RL/2.*(EXP(REAL(I1-1)/((NX-1)/2.)*FAC)-1.)/(EXP(FAC)-1.)
            TMX=RL/2.0-RL/2.*REAL(I1-1)/((NX-1)/2.)
        END IF


        !RD=DFF/COS(SUA)+TMX*TAN(SUA)

        !FX(NFX-I+1,J)=((DARL-RL)*(EXP(REAL(I-1)/(NFX-1)*FAC)-1.)/(EXP(FAC)-1.)+RL)*COS((-PI/(NFY-1))*(J-1))

		DO J=1,NZ
			TMZ=((-TMX*TAN(SUA)-DFF)+WPLZ(nfxf+I)/COS(SUA))/(NZ-1)*(J-1)-TMX*TAN(SUA) !+WPLAZ(I)/COS(SUA)

			!TMZ=-TMX*TAN(SUA)-RD/(NZ-1)*(J-1)
			
            !TMZ=(-TMX*TAN(SUA)-DFF)/(NZ-1)*(J-1)-TMX*TAN(SUA) !+WPLAZ(I)/COS(SUA)
			
            !TMY=RB/2.0*(1-(2*TMX/RL)**2)*(1-(TMZ/RD)**2) !*(1+0.2*(2*TMX/RL)**2)
			TMZ1=(-TMX*TAN(SUA)+DFF)+WPLZ(nfxf+I)/COS(SUA)
			TMZ=(-TMX*TAN(SUA)+WPLZ(nfxf+I)/COS(SUA))-TMZ1/(NZ-1)*(J-1)
            
            !DFF=RD-SINK/COS(SUA)

            IF(ABS(TMZ1).GE.RD) THEN
                IF(TMZ.LE.RD-DFF) THEN
                    TMY=RB/2.0*(1-(2*TMX/RL)**2)*(1-((TMZ+DFF-RD)/RD)**2) !*(1+0.2*(2*TMX/RL)**2)
			    ELSE
                    TMY=RB/2.0*(1-(2*TMX/RL)**2) !*(1-(TMZ/RD)**2)
                END IF
            ELSE 
                TMY=RB/2.0*(1-(2*TMX/RL)**2)*(1-((TMZ+DFF-RD)/RD)**2) !*(1+0.2*(2*TMX/RL)**2)
            END IF

            X(I,J)=TMX*COS(SUA)+TMZ*SIN(SUA)
            Y(I,J)=TMY
            Z(I,J)=TMZ*COS(SUA)+TMX*SIN(SUA) !+WPLAZ(I) !-RD/30. !-0.1*RD
            !Z(I,J)=TMZ
            !IF(J.EQ.1) THEN
                !Z(I,J)=Z(I,J)-DMESH
            !END IF

			TVX=-4*RB*TMX/(RL**2)*(1-(TMZ/RD)**2) !*(1+0.2*(2*TMX/RL)**2)+&
                !0.8*RB*TMX/(RL**2)*RB/2.0*(1-(2*TMX/RL)**2)*(1-(TMZ/RD)**2)
            TVY=-1.0
			TVZ=-RB*TMZ/(RD**2)*(1-(2*TMX/RL)**2)*(1+0.2*(2*TMX/RL)**2)
			
            VX1=TVX*COS(SUA)+TVZ*SIN(SUA)
            VY1=TVY
            VZ1=TVZ*COS(SUA)+TVX*SIN(SUA)

            VN=SQRT(VX1**2+VY1**2+VZ1**2)
			VX(I,J)=VX1/VN
			VY(I,J)=VY1/VN
			VZ(I,J)=VZ1/VN

            IF(MMESH.NE.2) THEN
            IF(J.EQ.1) THEN
                WPLX(I)=X(I,J)
                WPLY(I)=Y(I,J)
                WPLZ(I)=Z(I,J)
            END IF
            END IF
        END DO
    END DO

    !X(2,4)=X(2,6);Y(2,4)=Y(2,6);Z(2,4)=Z(2,6)
    !X(4,3)=X(4,4);Y(4,3)=Y(4,4);Z(4,3)=Z(4,4)

	K1=0
	DO I=1,NX
        NWPLAZ(I)=(I-1)*NZ+1
		DO J=1,NZ
			K1=K1+1
			!WRITE(13,104)X(I,J),Y(I,J),Z(I,J)-0.05*RD,VX(I,J),VY(I,J),VZ(I,J)
			!WRITE(17,104)X(I,J),Y(I,J),Z(I,J)-0.05*RD,VX(I,J),VY(I,J),VZ(I,J)
            !WRITE(18,104)X(I,J),Y(I,J),Z(I,J)-0.05*RD,VX(I,J),VY(I,J),VZ(I,J)
			!WRITE(19,104)X(I,J),Y(I,J),Z(I,J)-0.05*RD,VX(I,J),VY(I,J),VZ(I,J)

            WRITE(11,104)X(I,J),Y(I,J),Z(I,J),VX(I,J),VY(I,J),VZ(I,J)
            WRITE(13,104)X(I,J),Y(I,J),Z(I,J),VX(I,J),VY(I,J),VZ(I,J)
			WRITE(17,104)X(I,J),Y(I,J),Z(I,J),VX(I,J),VY(I,J),VZ(I,J)
            WRITE(18,104)X(I,J),Y(I,J),Z(I,J),VX(I,J),VY(I,J),VZ(I,J)
			WRITE(19,104)X(I,J),Y(I,J),Z(I,J),VX(I,J),VY(I,J),VZ(I,J)
        END DO
    END DO

    WPLAZ=TPWPLAZ
    
    DFX=ABS(X(1,1)-X(2,1)) !*1.2 !*1.2
    !DFX=1./20.

103	FORMAT(1X,I8,6F15.6)
104	FORMAT(1X,6F15.6)

    !GOTO 91

    FZ=0
    !DFX=LF/(NFXF)
    DFY=LB/(NFY-1)  !  old
    !DFY=LB/30.  !  new

	!DFY=DFX
	K2=0	
	VX2=0.0
	VY2=0.0
	VZ2=1.0

    !WRITE(*,*)"DFY",LB,DFY

!FIRST PART OF THE FREE SURFACE (UPSTREAM, BEFORE THE BOW)
	DO I=1,NFXF,2
		DO J=1,NFY,2  
            
            !FX(I,J)=OUTX(NFX,1)+LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            
            !FX(I,J)=X(1,1)-LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)

            !FX(I,J)=X(1,1)-LF*REAL(NFXF-I+1)/(NFXF)
            FX(I,J)=X(1,1)-DFX*(NFXF-I+1)

            IF(MGRID.EQ.1) THEN
            FY(I,J)=LB*REAL(J-1)/(NFY-1)
            ELSE	
            FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            END IF
            
            VX2=0.0
		    VY2=0.0
			VZ2=1.0
      
            K1=K1+1

			!WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			!WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			!WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    DO I=1,NFXF,2 
        DO J=1,NFY
            IF(I.GT.1) THEN
                FX(I-1,J)=(FX(I-2,J)+FX(I,J))/2.
                FY(I-1,J)=(FY(I-2,J)+FY(I,J))/2.
            END IF
        END DO
    END DO
	DO J=1,NFY,2 
        DO I=1,NFXF 
            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
            END IF
        END DO
    END DO

    DO I=1,NFXF-1
		DO J=1,NFY
            K1=K1+1
            IFE=NBPOINT+(I-1)*NFY+J
            FZ=ETAW(IFE)
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    DO J=1,NFY
            FX(1,J)=(FX(NFXF-1,J)+X(1,1))/2.
            FY(1,J)=FY(NFXF-1,J)

            IFE=NBPOINT+(NFXF-1)*NFY+J
            FZ=ETAW(IFE)
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
    END DO


!SECOND PART OF THE FREE SURFACE (ALONG SHIP)
	DO I=1,NFX
		DO J=1,NFY,2  
            !DFY=(LB-TY(NFX+NFXF-I+1))/(NFY-1)

            !FX(I,J)=OUTX(NFX-I+1,1)
            !FX(I,J)=X(I,1)
            !FY(I,J)=OUTX(NFX-I+1,2)+(LB-OUTX(NFX-I+1,2))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            !FY(I,J)=Y(I,1)+(LB-Y(I,2))*REAL(J-1)/(NFY-1)
            
            FX(I,J)=X(I,1)
            IF(MGRID.EQ.1) THEN
            FY(I,J)=Y(I,1)+(LB-Y(I,2))*REAL(J-1)/(NFY-1)
            ELSE
            FY(I,J)=Y(I,1)+(LB-Y(I,2))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            END IF
            
            VX2=0.0
		    VY2=0.0
			VZ2=1.0

            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.

                IFE=NBPOINT+(NFXF)*NFY+(I-1)*NFY+J
                FZ=ETAW(IFE)
                WRITE(13,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
			    WRITE(18,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
			    WRITE(19,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
            END IF

            IF(J.EQ.1) THEN
                !FY(I,J)=FY(I,J)+DMESH
            END IF

            K1=K1+1
			!WRITE(12,103)K2,FX(I),FY(I,J),FZ,VX2,VY2,VZ2
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
            IFE=NBPOINT+(NFXF)*NFY+(I-1)*NFY+J
            FZ=ETAW(IFE)
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
            !END IF
        END DO
    END DO

!THIRD PART OF THE FREE SURFACE (UPSTREAM, AFTER THE STERN)
    !DFX=LA/(NFXA)
    DFY=LB/(NFY-1)
    !DFY=LB/28
    DFY=LB/20.
    DFY=1.0/23.

	DO I=2,NFXA,2
		DO J=1,NFY,2  
            !FX(I,J)=TX(1)-(I-NFX)*DFX
			!FY(I,J)=(J-1)*DFY

            !FX(I,J)=OUTX(1,1)-LA*(EXP(REAL(I)/(NFXA)*FAC)-1.)/(EXP(FAC)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*(EXP(REAL(I)/(NFXA)*FACF)-1.)/(EXP(FACF)-1.)
            
            !FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*REAL(I)/(NFXA)
            FX(I,J)=X(NFX,1)+DFX*REAL(I)
            
            
            IF(MGRID.EQ.1) THEN
            FY(I,J)=LB*REAL(J-1)/(NFY-1)
            ELSE
            FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            END IF
            !FY(I,J)=(J-1)*DFY

            VX2=0.0
		    VY2=0.0
			VZ2=1.0
            !END IF
        END DO
    END DO

    DO I=2,NFXA,2 
	    DO J=1,NFY 
            IF(I.GT.2) THEN
                FX(I-1,J)=(FX(I-2,J)+FX(I,J))/2.
                FY(I-1,J)=(FY(I-2,J)+FY(I,J))/2.
            END IF
        END DO
    END DO
	DO J=1,NFY,2 
        DO I=1,NFXA
            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
            END IF
        END DO
    END DO

    DO J=1,NFY
            FX(1,J)=(FX(2,J)+X(NX,1))/2.
            FY(1,J)=FY(2,J)
            IFE=NBPOINT+(NFXF+NFX)*NFY+J
            FZ=ETAW(IFE)
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
    END DO

    DO I=2,NFXA
		DO J=1,NFY
            K1=K1+1
            IFE=NBPOINT+(NFXF+NFX)*NFY+(I-1)*NFY+J
            FZ=ETAW(IFE)
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

!THE PANEL INDEX ON THE SHIP SURFACE
 	K2=0
	DO I=1,NX-1,2
		DO J=1,NZ-1,2
			K2=K2+1
		    ND1=J+(I-1)*NZ
			ND2=ND1+1
            ND3=ND2+1
            ND4=ND3+NZ
            ND5=ND3+2*NZ
            ND6=ND5-1
            ND7=ND6-1
            ND8=ND7-NZ
			ND9=ND8+1
			WRITE(13,105)9,ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
			!WRITE(16,*)ND1,ND2,ND3,ND4
			WRITE(11,*)ND1,ND2,ND9,ND8
            WRITE(11,*)ND8,ND9,ND6,ND7
            WRITE(11,*)ND2,ND3,ND4,ND9
            WRITE(11,*)ND9,ND4,ND5,ND6
			WRITE(17,*)ND1,ND2,ND9,ND8
            WRITE(17,*)ND8,ND9,ND6,ND7
            WRITE(17,*)ND2,ND3,ND4,ND9
            WRITE(17,*)ND9,ND4,ND5,ND6
			WRITE(18,*)ND1,ND2,ND9,ND8
            WRITE(18,*)ND8,ND9,ND6,ND7
            WRITE(18,*)ND2,ND3,ND4,ND9
            WRITE(18,*)ND9,ND4,ND5,ND6
			WRITE(19,*)ND1,ND2,ND9,ND8
            WRITE(19,*)ND8,ND9,ND6,ND7
            WRITE(19,*)ND2,ND3,ND4,ND9
            WRITE(19,*)ND9,ND4,ND5,ND6
        END DO
    END DO
    105 FORMAT(1X,10I8)

!----	THE PANEL INDEX ON THE FREE SURFACE

	DO I=1,NFX+NFXA+NFXF-1,2
		DO J=1,NFY-1,2
            K3=K3+1
		    ND1=NBPOINT+J+(I-1)*NFY
			ND2=ND1+NFY !NBPOINT+J+I*NFY
            ND3=ND2+NFY
			ND4=ND3+1 !NBPOINT+J+1+I*NFY
			ND5=ND4+1
            ND6=ND5-NFY !NBPOINT+J+1+(I-1)*NFY
            ND7=ND6-NFY
            ND8=ND7-1
            ND9=ND8+NFY

			!WRITE(12,*)K3,NI,NI+1,NI+NFX+1,NI+NFX-1
			!WRITE(15,*)NI,NI+1,NI+NFX+1,NI+NFX-1
			!WRITE(18,*)NI,NI+1,NI+NFX+1,NI+NFX-1

            WRITE(13,105)9,ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
			WRITE(18,*)ND1,ND2,ND9,ND8
            WRITE(18,*)ND8,ND9,ND6,ND7
            WRITE(18,*)ND2,ND3,ND4,ND9
            WRITE(18,*)ND9,ND4,ND5,ND6
			WRITE(19,*)ND1,ND2,ND9,ND8
            WRITE(19,*)ND8,ND9,ND6,ND7
            WRITE(19,*)ND2,ND3,ND4,ND9
            WRITE(19,*)ND9,ND4,ND5,ND6
        END DO
    END DO

    !NFX=NFX+NFXF+NFXA

    DEALLOCATE(FX,FY)
	CLOSE(19)
	CLOSE(18)
	CLOSE(17)
	CLOSE(16)
	CLOSE(15)
	CLOSE(14)
	CLOSE(13)
	CLOSE(12)
	CLOSE(11)

    
    IF(MDFREE.EQ.2) THEN
    SUA=SUAT
    DFF=DFFT
    END IF
END SUBROUTINE