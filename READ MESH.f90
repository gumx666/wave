SUBROUTINE OFFSETIN
!!C!C!C!C!C!C!C!C!C!C!C!C!C!C!C!C!C**
!*    INPUT NODE CORD. AND ELEMENT'S NODE CODING
!!C!C!C!C!C!C!C!C!C!C!C!C!C!C!C!C!C**
!CNTPN: THE NUMBER OF THE VERTEXES
!CNTP: THE NUMBER OF THE ELEMENTS
!CCORD(NTPN,3):COORDINATES OF THE VERTEXES
!CINE(NTP,4):POINTER OF EACH ELEMENT
      USE GREENMOD
      USE INPUTDATA
    CHARACTER*7:: TITLE,TITLE2
	  !WRITE(*,*) 
	  !WRITE(*,100)
100   FORMAT(1X,'BEGIN OFFSETIN')
      OPEN(11,FILE='ALLINDEX.DAT')
!!CALLINDEX.DAT: GENERATE ALL THE INDEX AND POINTS ON TWO FACES
!!C设定最大面元数和顶点数，读入物体上的顶点数，面元数，界面上的顶点数，面元数
	NTPMAX=60000
	NTPNMAX=60000
	TPNTP=0
    READ(11,*)NTPN,NTP
    IF(NTPN.GT.NTPNMAX.OR.NTP.GT.NTPMAX) STOP 'TOO MANY POINTS OR PANELS'
    
    !SHALLOW WATER
    IF(NSWIN(13,1).LE.500) THEN
        !NBTPOINT=NSWIN(14,2)*NSWIN(14,1)
        !NBTBLOCK=(NSWIN(14,2)-1)*(NSWIN(14,1)-1)/4
        !NTPN=NTPN+NBTPOINT
        !NTP=NTP+NBTBLOCK
    END IF
    !
    
    IF(NIT.EQ.1.AND.ITTE.EQ.1) THEN
	    ALLOCATE(CORD(NTPN+500,3),VECN(NTPN+500,3),INE(NTP+500,9),PAN(NTP+500),VE(NTP+500,9,3),PCORD(NTPN+500,3))
        ALLOCATE(DMIP(NTPN+500),DMIE(NTP+500))
    END IF
!C读入第一层流体中物体上及界面上各点的坐标和面元顶点索引

    IF(ITTE.EQ.1) THEN
        DFYYY=0.3*DFY   !wigley
        DFYYY=(NSWIN(11,1))*DFY
        
        !DFYYY=0.35*DFY !S60
        !DFYYY=0.4*DFY !KCS
        !DFYYY=0.45*DFY
        !DFYYY=1.0*DFY 
        !DFYYY=0.35*DFY   !wigley
        IF(FR.GE.0.7) DFYYY=NSWIN(11,1)*DFY !DFYYY=0.4*DFY
    ELSE
        DFYYY=(NSWIN(11,1))*DFY
        !DFYYY=0.45*DFY !KCS
        !DFYYY=0.5*DFY !WIGLEY
        
        !DFYYY=0.4*DFY
        
        !DFYYY=0.35*DFY !S60
        !DFYYY=0.42*DFY
        !DFYYY=0.6*DFY   !wigley
        !DFYYY=0.3*DFY   !wigley
        DFYYY=NSWIN(11,1)*DFY
        
        IF(FR.GE.0.7) DFYYY=(NSWIN(11,1)+0.05)*DFY
        
        !DFYYY=0.3*DFY   !wigley
    END IF
    LOFFY=DFYYY
    
	DO I=1,NBPOINT+NFPOINT
	    READ(11,*)CORD(I,1),CORD(I,2),CORD(I,3),VECN(I,1),VECN(I,2),VECN(I,3)
        !IF(I.GT.NBPOINT.AND.ITTE.NE.1) CORD(I,3)=ETAW(I-NBPOINT+NBPOINT1)
    END DO
  
    !DO I=1,NWAPL
    !    WRITE(*,*)CORD(MARKWAP(I),:)
    !END DO
    !STOP
    
	DO I=1,NBBLOCK+NFBLOCK
        READ(11,*)PAN(I),INE(I,1:PAN(I))
        !READ(11,*)PAN(I),INE(I,1:4)
        !WRITE(*,*)PAN(I),INE(I,1:4)
        IF(PAN(I).EQ.9) THEN
            TPNTP=TPNTP+4
        ELSE
            TPNTP=TPNTP+1
        END IF
        !READ(11,*)K,(INE(I,J),J=1,PAN(I))
	END DO

    !IF(ITTE.EQ.1) THEN
        NBPOINTW=NBPOINT
        NBBLOCKW=NBBLOCK
    !END IF

    !IF(ITTE.GE.2)  CALL VIRTMESH
    
    PAN=9

    ITTT=3
    !IF(NSWIN(13,1)/RD.LE.2.0) ITTT=30
    !ITTT=30
    IF(ITTE.LE.ITTT) THEN
	DO I=1,NTPN
	    PCORD(I,:)=CORD(I,:)
        NBPOINTP=NBPOINT
        NBBLOCKP=NBBLOCK
    END DO
    END IF

    !IF(ITTE.GE.4) THEN 
    IF(ITTE.GE.ITTT+1) THEN 
    DO I=1,NFPOINT
	    CORD(I+NBPOINT,1:2)=PCORD(I+NBPOINTP,1:2)
        !IF(ITTE.LE.3) THEN
        !    CORD(I+NBPOINT,3)=ETAW1(I) !ETAW(I+NBPOINT1)
        !ELSE
        !    CORD(I+NBPOINT,3)=(0.5*ETAW0(I)+0.5*ETAW1(I))
        !END IF

        IF(MOD(ITTE,2).EQ.0) THEN
            !CORD(I+NBPOINT,3)=(0.4*ETAW0(I)+0.6*ETAW1(I))
        END IF
    END DO 
    END IF

	DO I=1,NTPN
        IF(I.GT.NBPOINT) THEN
            CORD(I,2)=CORD(I,2)+DFYYY
        END IF
        
        !IF(I.GT.NBPOINT+(NFXF+NFX)*NFY) THEN   
    END DO
    
    DO I=1,NFX
        !IP=NBPOINT+(NFXF+I-1)*NFY+1
        !CORD(IP,2)=CORD(IP,2)+DFYYY
        !CORD(IP+1,2)=CORD(IP+1,2)+DFYYY
    END DO
        
    IF(ITTE.GE.2) THEN
	DO I=1,NTPN
        IF(I.LE.NBPOINT) THEN
            !CORD(I,3)=CORD(I,3)+0.0625/10.
        END IF       
    END DO    
    END IF
    
	DO I=1,NFX
        !IP=NBPOINT+NFXF*NFY+(I-1)*NFY+1
        !    CORD(IP,2)=CORD(IP,2)+DFYYY
        !    CORD(IP+1,2)=(CORD(IP,2)+CORD(IP+2,2))/2.
            !CORD(I,1)=CORD(I,1)+LOFFX
    END DO


    IF(ITTE.EQ.1) CALL SMOOTH_MESH_1
    IF(MTROM.EQ.1) CALL TRANSOMMESH(DFYYY)
    !IF(NSWIN(13,1).LE.500) CALL BOTTOM

    !IF(ITTE.GE.3) THEN 
    IF(ITTE.GE.2) THEN 
    DO I=1,NFPOINT
        CORD(I+NBPOINT,3)=ETAW0(I)+NSWIN(15,1)*(ETAW1(I)-ETAW0(I))
        !CORD(I+NBPOINT,3)=ETAW1(I)
        !CORD(I+NBPOINT,3)=0.
    END DO
    
    DO I=1,NBPOINT
        !CORD(I,3)=CORD(I,3)+0.0625/15.
    END DO   
    ELSE
    DO I=1,NFPOINT
        !CORD(I+NBPOINT,3)=ETAW0(I)+0.9*(ETAW1(I)-ETAW0(I))
        CORD(I+NBPOINT,3)=0.
    END DO 
    END IF

    IF(MDFREE.EQ.3) CORD(NBPOINT+1:NBPOINT+NFPOINT,3)=0    !LINEAR FREESURFACE
    
    !WRITE(*,*)"AVERAGE      ",SUM(CORD(NBPOINT+1:NTPN,3))/NFPOINT
    
    VE=3
    DO I=1,NBBLOCK+NFBLOCK
        DO J=1,PAN(I)
            VE(I,J,:)=VECN(INE(I,J),:)
            !WRITE(*,*)SQRT(VE(I,J,1)**2+VE(I,J,2)**2+VE(I,J,3)**2)
        END DO
    END DO

!   变换间断单元
    !CALL DISCOUNT
    
    !STOP
    IF(NIT.EQ.1.AND.ITTE.EQ.1) THEN
        ALLOCATE(DNODE(NTPN+500))
    END IF
    
    TPNTP=0
	DO I=1,NTP
        IF(PAN(I).EQ.9) THEN
            TPNTP=TPNTP+4
        ELSE
            TPNTP=TPNTP+1
        END IF
        !READ(11,*)K,(INE(I,J),J=1,PAN(I))
	END DO


    !CORD(INE(41,1),:)=CORD(INE(41,2),:)
    !CORD(INE(41,3),:)=CORD(INE(41,2),:)
    !CORD(INE(41,4),:)=(CORD(INE(41,3),:)+CORD(INE(41,5),:))/2.
    !CORD(INE(41,8),:)=(CORD(INE(41,1),:)+CORD(INE(41,7),:))/2.

    
    CALL ELEMENTAREA
    CALL VECTOR
    
    
    !SHIP MEAN WETTED SORFACE
    IF(ITTE.EQ.1) THEN
	    CORD0(1:NBPOINT,:)=CORD(1:NBPOINT,:)
	    VECN0(1:NBPOINT,:)=VECN(1:NBPOINT,:)
    END IF

    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(17,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\PANELVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    !OPEN(18,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREEVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    !OPEN(19,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\ALLVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')

    OPEN(17,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\PANELVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(18,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREEVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(19,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\ALLVECTOR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    
    
    
	WRITE(17,101)NBPOINT,NBBLOCK*4
	WRITE(18,101)NFPOINT,(NFBLOCK)*4
	WRITE(19,101)NTPN,NTP*4
101	format(1x,'title="episode solid mesh"'/1x,'variables="x","Y",& 
     "Z","NX","NY","NZ"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')

    DO I=1,NBPOINT
        WRITE(17,111)CORD(I,1),CORD(I,2),CORD(I,3),VECN(I,1:3)
        WRITE(19,111)CORD(I,1),CORD(I,2),CORD(I,3),VECN(I,1:3)
    END DO

    DO I=NBPOINT+1,NTPN
        !WRITE(18,111)CORD(I,1),CORD(I,2)-DFYYY,CORD(I,3),VECN(I,1:3)
        !WRITE(19,111)CORD(I,1),CORD(I,2)-DFYYY,CORD(I,3),VECN(I,1:3)
        WRITE(18,111)CORD(I,1),CORD(I,2),CORD(I,3),VECN(I,1:3)
        WRITE(19,111)CORD(I,1),CORD(I,2),CORD(I,3),VECN(I,1:3)
    END DO
111FORMAT(6F15.6)



    DO I=1,NBBLOCK
        WRITE(17,*)INE(I,1),INE(I,2),INE(I,9),INE(I,8)
        WRITE(17,*)INE(I,2),INE(I,3),INE(I,4),INE(I,9)
        WRITE(17,*)INE(I,9),INE(I,4),INE(I,5),INE(I,6)
        WRITE(17,*)INE(I,8),INE(I,9),INE(I,6),INE(I,7)
        WRITE(19,*)INE(I,1),INE(I,2),INE(I,9),INE(I,8)
        WRITE(19,*)INE(I,2),INE(I,3),INE(I,4),INE(I,9)
        WRITE(19,*)INE(I,4),INE(I,5),INE(I,6),INE(I,9)
        WRITE(19,*)INE(I,6),INE(I,7),INE(I,8),INE(I,9)
    END DO

    DO I=NBBLOCK+1,NTP
        WRITE(18,*)INE(I,1)-NBPOINT,INE(I,2)-NBPOINT,INE(I,9)-NBPOINT,INE(I,8)-NBPOINT
        WRITE(18,*)INE(I,2)-NBPOINT,INE(I,3)-NBPOINT,INE(I,4)-NBPOINT,INE(I,9)-NBPOINT
        WRITE(18,*)INE(I,4)-NBPOINT,INE(I,5)-NBPOINT,INE(I,6)-NBPOINT,INE(I,9)-NBPOINT
        WRITE(18,*)INE(I,6)-NBPOINT,INE(I,7)-NBPOINT,INE(I,8)-NBPOINT,INE(I,9)-NBPOINT
        WRITE(19,*)INE(I,1),INE(I,2),INE(I,9),INE(I,8)
        WRITE(19,*)INE(I,2),INE(I,3),INE(I,4),INE(I,9)
        WRITE(19,*)INE(I,4),INE(I,5),INE(I,6),INE(I,9)
        WRITE(19,*)INE(I,6),INE(I,7),INE(I,8),INE(I,9)
    END DO

    CLOSE(17)
    CLOSE(18)
    CLOSE(19)

!CWRITE(*,*)'SUBROUTINE OFFSETIN FINISHED'	
	CLOSE(11)
END SUBROUTINE OFFSETIN

SUBROUTINE BOTTOM
    USE INPUTDATA
    USE GREENMOD
    
    REAL:: DTFX,DTFY
    !NSWIN(13,1)--DEPTH
    
    !NSWIN(14,1)--NBX; NSWIN(14,2)--NBY 
    
    DTFX=ABS(CORD(NBPOINT+NFPOINT,1)-CORD(NBPOINT+1,1))/(NSWIN(14,1)-1)
    DTFY=ABS(CORD(NBPOINT+NFY,2)-CORD(NBPOINT+1,2))/(NSWIN(14,2)-1)
    DO I=1,NSWIN(14,1)
        DO J=1,NSWIN(14,2)
            CORD(NBPOINT+NFPOINT+(I-1)*NSWIN(14,2)+J,1)=CORD(NBPOINT+1,1)+(I-1)*DTFX
            CORD(NBPOINT+NFPOINT+(I-1)*NSWIN(14,2)+J,2)=(J-1)*DTFY
            CORD(NBPOINT+NFPOINT+(I-1)*NSWIN(14,2)+J,3)=-NSWIN(13,1)
        END DO
    END DO
    
    IE=NBBLOCK+NFBLOCK
    DO I=1,NSWIN(14,1)-1,2
        DO J=1,NSWIN(14,2)-1,2
            IE=IE+1
            INE(IE,1)=NBPOINT+NFPOINT+(I-1)*NSWIN(14,2)+J
            INE(IE,2)=INE(IE,1)+1
            INE(IE,3)=INE(IE,2)+1
            INE(IE,4)=INE(IE,3)+NSWIN(14,2)
            INE(IE,5)=INE(IE,4)+NSWIN(14,2)
            INE(IE,6)=INE(IE,5)-1
            INE(IE,7)=INE(IE,6)-1
            INE(IE,8)=INE(IE,7)-NSWIN(14,2)
            INE(IE,9)=INE(IE,8)+1
        END DO
    END DO
    
END SUBROUTINE

SUBROUTINE TRANSOMMESH(DFYYY)
    USE GREENMOD
    REAL:: WPVI1(200,3),WPVI2(200,3),WPVI3(200,3)
    INTEGER:: NPVI1,NPVI2
    REAL TPCOR(500,3),UI(500)
    !REAL,ALLOCATABLE:: TP(:,:)
    !INTEGER,DIMENSION(:,:),ALLOCATABLE:: INTP 
    REAL:: DFYYY
    
    !CALL SHIPHULL(SINK,TRI,NWP,WX(1:NWP),WZ(1:NWP),NFX1,WX1,WY1,XTR1,ZTR1,TAOX1,NTAO,ZFAC1,WPVI1,NPVI1,WPVI2,NPVI2,NZ1)
    
    NZS=3
    DYS=(CORD(NBPOINT+1,2)-DFYYY)/(NZS-1)
    !DYS=(CORD(NBPOINT+1,2))/(NZS-1)
    DO I=1,NFXF+1
        IR=NBPOINT+(I-1)*NFY+1
        DO J=1,NZS
            IS=NTPN+(I-1)*NZS+J
            CORD(IS,1)=CORD(IR,1)
            !CORD(IS,1)=CORD(NBPOINT+1,1)+(CORD(NBPOINT+NFXF*NFY+1,1)-DFX-CORD(NBPOINT+1,1))/(NFXF)*(I-1)
            CORD(IS,2)=(J-1)*DYS
        END DO
    END DO

    !WRITE(*,*)
    !WRITE(*,*)"INITIAL    NFPOINT",NFPOINT
    !WRITE(*,*)

    NSPOINT=NZS*(NFXF+1)
    !NSBLOCK=(NZS*NFXF)/4
    NSBLOCK=(NZS-1)*NFXF/4

    DO I=1,(NFXF),2
        DO J=1,(NZS-1),2
        IE=(I-1)/2*(NZS-1)/2+(J-1)/2+1+NBBLOCK+NFBLOCK
        INE(IE,1)=NBPOINT+NFPOINT+(I-1)*NZS+J
        INE(IE,2)=INE(IE,1)+NZS
        INE(IE,3)=INE(IE,2)+NZS
        INE(IE,4)=INE(IE,3)+1
        INE(IE,5)=INE(IE,4)+1
        INE(IE,6)=INE(IE,5)-NZS
        INE(IE,7)=INE(IE,6)-NZS
        INE(IE,8)=INE(IE,7)-1
        INE(IE,9)=INE(IE,8)+NZS
        !WRITE(*,*)IE

        IF(J.EQ.NZS-1) THEN
            !INE(IE,5)=NBPOINT+(I+1)*NFY+1
            !INE(IE,6)=INE(IE,5)-NFY
            !INE(IE,7)=INE(IE,6)-NFY
        END IF
        END DO
    END DO

    !WRITE(*,*)"NODES ON TRANSOM:  ",NSPOINT
    
    NFPOINT=NFPOINT+NSPOINT
    NFBLOCK=NFBLOCK+NSBLOCK
    NTPN=NTPN+NSPOINT
    NTP=NTP+NSBLOCK
    !write(*,*)NTP

    !DEALLOCATE(TP,INTP)
END SUBROUTINE

SUBROUTINE SMOOTH_MESH
    USE GREENMOD
    REAL TPCOR(500,3),UI(500)

    DO J=1,NFXF+NFX+NFXA
        UI(J)=1./(NFXF+NFX+NFXA-1)*(J-1)
    END DO
    
    DO J=2,2
        DO I=1,NFXF+NFX+NFXA
            TPCOR(I,:)=CORD(NBPOINT+(I-1)*NFY+J,:)
        END DO
        CALL NURBS_INT(NFXF+NFX+NFXA-1,TPCOR(1:NFXF+NFX+NFXA,:),NFXF+NFX+NFXA-1,UI(1:NFXF+NFX+NFXA),TPCOR(1:NFXF+NFX+NFXA,:))

        DO I=1,NFXF+NFX+NFXA
            CORD(NBPOINT+(I-1)*NFY+J,:)=TPCOR(I,:)
        END DO
    END DO
    DO J=3,3
        DO I=1,NFXF+NFX+NFXA
            CORD(NBPOINT+(I-1)*NFY+J,1)=(2.*CORD(NBPOINT+(I-1)*NFY+J-1,1)-CORD(NBPOINT+(I-1)*NFY+J-2,1))
            CORD(NBPOINT+(I-1)*NFY+J,2)=(2.*CORD(NBPOINT+(I-1)*NFY+J-1,2)-CORD(NBPOINT+(I-1)*NFY+J-2,2))
        END DO
    END DO
    
    DO J=5,NFY,2
        DO I=1,NFXF+NFX+NFXA
            TPCOR(I,:)=CORD(NBPOINT+(I-1)*NFY+J,:)
        END DO
        !WRITE(*,*)NFXF+NFX+NFXA,TPCOR(1:NFXF+NFX+NFXA,1)
        !WRITE(*,*)
        !STOP
        CALL NURBS_INT(NFXF+NFX+NFXA-1,TPCOR(1:NFXF+NFX+NFXA,:),NFXF+NFX+NFXA-1,UI(1:NFXF+NFX+NFXA),TPCOR(1:NFXF+NFX+NFXA,:))
        DO I=1,NFXF+NFX+NFXA
            CORD(NBPOINT+(I-1)*NFY+J,:)=TPCOR(I,:)
        END DO
    END DO
    DO J=4,NFY-1,2
        DO I=1,NFXF+NFX+NFXA
            CORD(NBPOINT+(I-1)*NFY+J,1)=(CORD(NBPOINT+(I-1)*NFY+J+1,1)+CORD(NBPOINT+(I-1)*NFY+J-1,1))/2.
            CORD(NBPOINT+(I-1)*NFY+J,2)=(CORD(NBPOINT+(I-1)*NFY+J+1,2)+CORD(NBPOINT+(I-1)*NFY+J-1,2))/2.
        END DO
    END DO

END SUBROUTINE

SUBROUTINE SMOOTH_MESH_1
    USE GREENMOD
    REAL TPCOR(500,3),UI(500)

    !GOTO 20
    !GOTO 10
    DO J=1,NFXF+NFX+NFXA
        UI(J)=1./(NFXF+NFX+NFXA-1)*(J-1)
    END DO
    DO J=1,NFX+NFXA
        !UI(J)=1./(NFX+NFXA-1)*(J-1)
    END DO
    DO J=1,NFY
        DO I=1,NFXF+NFX-1
            TPCOR(I,:)=CORD(NBPOINT+(I-1)*NFY+J,:)
            TPCOR(I,3)=0.
        END DO
        DO I=NFXF+NFX+1,NFXF+NFX+NFXA
            TPCOR(I-1,:)=CORD(NBPOINT+(I-1)*NFY+J,:)
            TPCOR(I-1,3)=0.
        END DO
        
        DO I=1,NFXF+NFX+NFXA
            !TPCOR(I,:)=CORD(NBPOINT+(I-1)*NFY+J,:)
            !TPCOR(I,3)=0.
        END DO

        !CALL NURBS_INT(NFX+NFXA-2,TPCOR(NFXF+1:NFXF+NFX+NFXA-1,:),NFX+NFXA-1,UI(1:NFX+NFXA),TPCOR(NFXF+1:NFXF+NFX+NFXA,:))
        CALL NURBS_INT(NFXF+NFX+NFXA-2,TPCOR(1:NFXF+NFX+NFXA-1,:),NFXF+NFX+NFXA-1,UI(1:NFXF+NFX+NFXA),TPCOR(1:NFXF+NFX+NFXA,:))
        !CALL NURBS_INT(NFXF+NFX+NFXA-1,TPCOR(1:NFXF+NFX+NFXA,:),NFXF+NFX+NFXA-1,UI(1:NFXF+NFX+NFXA),TPCOR(1:NFXF+NFX+NFXA,:))
        DO I=1,NFXF+NFX+NFXA
            CORD(NBPOINT+(I-1)*NFY+J,1:2)=TPCOR(I,1:2)
        END DO
    END DO
    10 CONTINUE

    GOTO 20
    FACB=2.0
    FACB=1.5
    FACB=1.2
    FACB=0.1
    DO J=1,NFY
        UI(J)=1./(NFY-1)*(J-1)
        UI(J)=(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
    END DO
    DO I=1,NFXF+NFX+NFXA
        DO J=1,NFY
            TPCOR(J,:)=CORD(NBPOINT+(I-1)*NFY+J,:)
            TPCOR(I,3)=0.
        END DO
        CALL NURBS_INT(NFY-1,TPCOR(1:NFY,:),NFY-1,UI(1:NFY),TPCOR(1:NFY,:))
        DO J=1,NFY
            CORD(NBPOINT+(I-1)*NFY+J,:)=TPCOR(J,1:2)
        END DO
    END DO

    DO J=2,NFY-1,2
        DO I=1,NFXF+NFX+NFXA
            !CORD(NBPOINT+(I-1)*NFY+J,1)=(CORD(NBPOINT+(I-1)*NFY+J+1,1)+CORD(NBPOINT+(I-1)*NFY+J-1,1))/2.
            CORD(NBPOINT+(I-1)*NFY+J,2)=(CORD(NBPOINT+(I-1)*NFY+J+1,2)+CORD(NBPOINT+(I-1)*NFY+J-1,2))/2.
        END DO
    END DO
    20 CONTINUE
    
       
END SUBROUTINE

SUBROUTINE DISCOUNT
    USE GREENMOD

    REAL:: CORD1(3)
    INTEGER:: ISID(3),IWAT(3),IVAN(9),COUT,IDCON(100),IDCON1(100)
    REAL:: AJ,XI,ETA,D1,D2,D3,D4
    REAL:: sn(9),snx(9),sne(9),sr(3),dn(3),CORL1(8,3)
    REAL,ALLOCATABLE:: CORDTP(:,:)
    ALLOCATE(CORL(9,3))
    ALLOCATE(CORDTP(NTPN+500,3))

    D1=0.9
    D2=D1
    D3=D2
    D4=D3

!   寻找需要改变的间断单元
    COUT=0
    DO IE=1,NBBLOCK
        DO J=1,9
            IF(ABS(CORD(INE(IE,J),3)).LE.1.E-4) THEN
                PAN(IE)=8
                COUT=COUT+1
                IDCON(COUT)=IE
                IDCON1(COUT)=IE
                ICORD(IE)=COUT
                DO K=1,9
                    CORDO(ICORD(IE),K,1:3)=CORD(INE(IE,K),:)
                END DO
                GOTO 10
            END IF
        END DO
        10 CONTINUE
    END DO


!   对间断单元按X坐标排序
    DO I=COUT-1,1,-1
        DO  J=1,I     
            IF(CORD(INE(IDCON1(J),9),1).GT.CORD(INE(IDCON1(J+1),9),1)) THEN
                TPX=IDCON1(J)
                IDCON1(J)=IDCON1(J+1)
                IDCON1(J+1)=TPX
            END IF
        END DO
    END DO
    !WRITE(*,*)"COUT        ",COUT
    !DO I=1,COUT
        !WRITE(*,*)IDCON(I),IDCON1(I),CORD(INE(IDCON1(I),9),1)
    !END DO

!   改变首尾间断单元
        IE=IDCON1(1)
        DO J=1,7,2
            IF(ABS(CORD(INE(IE,J),3)).LE.1.E-4.AND.ABS(CORD(INE(IE,J),2)).LE.1.E-4) THEN
                ISID(1)=J
            END IF
        END DO

!       将面源节点内收
            DO J=1,9
                CORL(J,:)=CORD(INE(IE,J),:)
            END DO

!           J=1
            XI=-D1
            ETA=-D3
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(1,:)=SR(:)
!           J=2
            XI=0.
            ETA=-D3
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(5,:)=SR(:)
!           J=3
            XI=D2
            ETA=-D3
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(2,:)=SR(:)
!           J=4
            XI=D2
            ETA=0.
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(6,:)=SR(:)
!           J=5
            XI=D2
            ETA=D4
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(3,:)=SR(:)
!           J=6
            XI=0.
            ETA=D4
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(7,:)=SR(:)
!           J=7
            XI=-D1
            ETA=D4
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(4,:)=SR(:)
!           J=8
            XI=-D1
            ETA=0.
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(8,:)=SR(:)

!       查看面源多余的边界点 
        IVAN=0
        IF(ISID(1).EQ.1) THEN
            IVAN(1)=1;IVAN(2)=2;IVAN(3)=8;IVAN(4)=9
        END IF
        IF(ISID(1).EQ.3) THEN
            IVAN(1)=3;IVAN(2)=4;IVAN(3)=2;IVAN(4)=9
        END IF
        IF(ISID(1).EQ.5) THEN
            IVAN(1)=5;IVAN(2)=6;IVAN(3)=4;IVAN(4)=9           
        END IF
        IF(ISID(1).EQ.7) THEN
            IVAN(1)=7;IVAN(2)=8;IVAN(3)=6;IVAN(4)=9
        END IF 

!       去掉面元多余的点
!       对IVAN(I)排序     
        DO I=4-1,1,-1
            DO J=1,I     
                IF(INE(IE,IVAN(J)).GT.INE(IE,IVAN(J+1))) THEN
                    TPX=IVAN(J)
                    IVAN(J)=IVAN(J+1)
                    IVAN(J+1)=TPX
                END IF
            END DO
        END DO


        CORDTP(1:NBPOINT,:)=CORD(1:NBPOINT,:)
        DO I=4,1,-1
            DO J=INE(IE,IVAN(I)),NBPOINT-1
                CORDTP(J,:)=CORDTP(J+1,:)
            END DO

            DO J=1,NBBLOCK
                DO L=1,9
                    IF(INE(J,L).GT.INE(IE,IVAN(I))) THEN
                        INE(J,L)=INE(J,L)-1
                    END IF
                END DO
            END DO
        END DO
            DO J=NBBLOCK+1,NTP
                DO L=1,9
                    INE(J,L)=INE(J,L)-4
                END DO
            END DO

        DO J=NBPOINT+1,NTPN
            CORDTP(J-4,:)=CORD(J,:)
        END DO

        NBPOINT=NBPOINT-4
        NTPN=NBPOINT+NFPOINT

        CORD(1:NTPN,:)=CORDTP(1:NTPN,:)

!       赋值给间断元
        CORDTP(1:NBPOINT,:)=CORD(1:NBPOINT,:)
        DO I=1,8
            CORDTP(NBPOINT+1:NBPOINT+8,:)=CORL1(1:8,:)
            INE(IE,I)=NBPOINT+I
        END DO
        CORDTP(NBPOINT+9:NTPN+8,:)=CORD(NBPOINT+1:NTPN,:)
        
        DO J=NBBLOCK+1,NTP
            DO L=1,9
                INE(J,L)=INE(J,L)+8
            END DO
        END DO

        NBPOINT=NBPOINT+8
        NTPN=NBPOINT+NFPOINT
        
        CORD(1:NTPN,:)=CORDTP(1:NTPN,:)

!   改变中间的间断元
    DO K=2,COUT-1
        IE=IDCON1(K)
        PAN(IE)=8
        COUT1=0
        DO J=1,7,2
            IF(ABS(CORD(INE(IE,J),3)).LE.1.E-4) THEN
                COUT1=COUT1+1
                ISID(COUT1)=J
                !WRITE(*,*)"ISID(COUT1)      ",ISID(COUT1)
            END IF
        END DO

        IF(CORD(INE(IE,ISID(2)),1).LE.CORD(INE(IE,ISID(1)),1)) THEN
            ISID(1)=ISID(2)
        END IF

        IF(K.NE.COUT) THEN
!       将面源节点内收
        IF(PAN(IE).EQ.8) THEN
            DO J=1,9
                CORL(J,:)=CORD(INE(IE,J),:)
            END DO

!           J=1
            XI=-D1
            ETA=-D3
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(1,:)=SR(:)
!           J=2
            XI=0.
            ETA=-D3
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(5,:)=SR(:)
!           J=3
            XI=D2
            ETA=-D3
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(2,:)=SR(:)
!           J=4
            XI=D2
            ETA=0.
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(6,:)=SR(:)
!           J=5
            XI=D2
            ETA=D4
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(3,:)=SR(:)
!           J=6
            XI=0.
            ETA=D4
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(7,:)=SR(:)
!           J=7
            XI=-D1
            ETA=D4
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(4,:)=SR(:)
!           J=8
            XI=-D1
            ETA=0.
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(8,:)=SR(:)
        END IF
        END IF

!       查看面源多余的边界点 
        IVAN=0
        IF(ISID(1).EQ.1) THEN
            IVAN(1)=1;IVAN(2)=2;IVAN(3)=8;IVAN(4)=9
        END IF
        IF(ISID(1).EQ.3) THEN
            IVAN(1)=3;IVAN(2)=4;IVAN(3)=2;IVAN(4)=9
        END IF
        IF(ISID(1).EQ.5) THEN
            IVAN(1)=5;IVAN(2)=6;IVAN(3)=4;IVAN(4)=9           
        END IF
        IF(ISID(1).EQ.7) THEN
            IVAN(1)=7;IVAN(2)=8;IVAN(3)=6;IVAN(4)=9
        END IF 

!       去掉面元多余的点
!       对IVAN(I)排序     
        DO I=4-1,1,-1
            DO J=1,I     
                IF(INE(IE,IVAN(J)).GT.INE(IE,IVAN(J+1))) THEN
                    TPX=IVAN(J)
                    IVAN(J)=IVAN(J+1)
                    IVAN(J+1)=TPX
                END IF
            END DO
        END DO

        CORDTP(1:NBPOINT,:)=CORD(1:NBPOINT,:)
        DO I=4,1,-1
            DO J=INE(IE,IVAN(I)),NBPOINT-1
                CORDTP(J,:)=CORDTP(J+1,:)
            END DO

            DO J=1,NBBLOCK
                DO L=1,9
                    IF(INE(J,L).GT.INE(IE,IVAN(I))) THEN
                        INE(J,L)=INE(J,L)-1
                    END IF
                END DO
            END DO
        END DO
            DO J=NBBLOCK+1,NTP
                DO L=1,9
                    INE(J,L)=INE(J,L)-4
                END DO
            END DO

        DO J=NBPOINT+1,NTPN
            CORDTP(J-4,:)=CORD(J,:)
        END DO

        NBPOINT=NBPOINT-4
        NTPN=NBPOINT+NFPOINT

        CORD(1:NTPN,:)=CORDTP(1:NTPN,:)

!       赋值给间断元
        CORDTP(1:NBPOINT,:)=CORD(1:NBPOINT,:)
        DO I=1,8
            CORDTP(NBPOINT+1:NBPOINT+8,:)=CORL1(1:8,:)
            INE(IE,I)=NBPOINT+I
        END DO
        CORDTP(NBPOINT+9:NTPN+8,:)=CORD(NBPOINT+1:NTPN,:)
        DO J=NBBLOCK+1,NTP
            DO L=1,9
                INE(J,L)=INE(J,L)+8
            END DO
        END DO

        NBPOINT=NBPOINT+8
        NTPN=NBPOINT+NFPOINT
        
        CORD(1:NTPN,:)=CORDTP(1:NTPN,:)
    END DO

!   尾部
        IE=IDCON1(COUT)
        DO J=2,8,2
            IF(ABS(CORD(INE(IE,J),3)).LE.1.E-4) THEN
                ISID(1)=J
            END IF
        END DO

!       将面源节点内收
        IF(PAN(IE).EQ.8) THEN
            DO J=1,9
                CORL(J,:)=CORD(INE(IE,J),:)
            END DO

!           J=1
            XI=-D1
            ETA=-D3
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(1,:)=SR(:)
!           J=2
            XI=0.
            ETA=-D3
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(5,:)=SR(:)
!           J=3
            XI=D2
            ETA=-D3
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(2,:)=SR(:)
!           J=4
            XI=D2
            ETA=0.
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(6,:)=SR(:)
!           J=5
            XI=D2
            ETA=D4
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(3,:)=SR(:)
!           J=6
            XI=0.
            ETA=D4
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(7,:)=SR(:)
!           J=7
            XI=-D1
            ETA=D4
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(4,:)=SR(:)
!           J=8
            XI=-D1
            ETA=0.
            CALL ISOPAR_9(XI,ETA,CORL,SN,SNX,SNE,AJ,SR,DN)
            CORL1(8,:)=SR(:)
        END IF

!       查看面源多余的边界点 
        IVAN=0
        IF(ISID(1).EQ.2) THEN
            IVAN(1)=1;IVAN(2)=2;IVAN(3)=3;IVAN(4)=4;IVAN(5)=8;IVAN(6)=9
        END IF
        IF(ISID(1).EQ.4) THEN
            IVAN(1)=2;IVAN(2)=3;IVAN(3)=4;IVAN(4)=5;IVAN(5)=6;IVAN(6)=9
        END IF
        IF(ISID(1).EQ.6) THEN
            IVAN(1)=4;IVAN(2)=5;IVAN(3)=6;IVAN(4)=7;IVAN(5)=8;IVAN(6)=9         
        END IF
        IF(ISID(1).EQ.8) THEN
            IVAN(1)=1;IVAN(2)=2;IVAN(3)=6;IVAN(4)=7;IVAN(5)=8;IVAN(6)=9
        END IF 

!       去掉面元多余的点
!       对IVAN(I)排序     
        DO I=6-1,1,-1
            DO J=1,I     
                IF(INE(IE,IVAN(J)).GT.INE(IE,IVAN(J+1))) THEN
                    TPX=IVAN(J)
                    IVAN(J)=IVAN(J+1)
                    IVAN(J+1)=TPX
                END IF
            END DO
        END DO

        CORDTP(1:NBPOINT,:)=CORD(1:NBPOINT,:)
        DO I=6,1,-1
            DO J=INE(IE,IVAN(I)),NBPOINT-1
                CORDTP(J,:)=CORDTP(J+1,:)
            END DO

            DO J=1,NBBLOCK
                DO L=1,9
                    IF(INE(J,L).GT.INE(IE,IVAN(I))) THEN
                        INE(J,L)=INE(J,L)-1
                    END IF
                END DO
            END DO
        END DO
            DO J=NBBLOCK+1,NTP
                DO L=1,9
                    INE(J,L)=INE(J,L)-6
                END DO
            END DO

        DO J=NBPOINT+1,NTPN
            CORDTP(J-6,:)=CORD(J,:)
        END DO

        NBPOINT=NBPOINT-6
        NTPN=NBPOINT+NFPOINT

        CORD(1:NTPN,:)=CORDTP(1:NTPN,:)

!       赋值给间断元
        CORDTP(1:NBPOINT,:)=CORD(1:NBPOINT,:)
        DO I=1,8
            CORDTP(NBPOINT+1:NBPOINT+8,:)=CORL1(1:8,:)
            INE(IE,I)=NBPOINT+I
        END DO
        CORDTP(NBPOINT+9:NTPN+8,:)=CORD(NBPOINT+1:NTPN,:)
        DO J=NBBLOCK+1,NTP
            DO L=1,9
                INE(J,L)=INE(J,L)+8
            END DO
        END DO

        NBPOINT=NBPOINT+8
        NTPN=NBPOINT+NFPOINT
        
        CORD(1:NTPN,:)=CORDTP(1:NTPN,:)

    DEALLOCATE(CORL,CORDTP)

    RETURN
END SUBROUTINE

SUBROUTINE VECTOR
    USE GREENMOD
    REAL:: AJ,XI,ETA,D1,D2,D3,D4
    REAL:: sn(9),snx(9),sne(9),sr(3),dn(3)
    ALLOCATE(CORL(9,3))

    D1=0.9
    D2=D1
    D3=D2
    D4=D3

    OPEN(10,FILE='HULLPANEL_NVEC.DAT')
    OPEN(11,FILE='ALLPANEL_NVEC.DAT')

    !CORD(INE(15,5),:)=0.5*(CORD(INE(15,4),:)+CORD(INE(15,6),:))
    VECN=0.
    DNODE=0.
    DMA=0.
    DO I=1,NTP
        
        DO J=1,PAN(I)
            CORL(J,:)=CORD(INE(I,J),:)
            DNODE(INE(I,J))=DNODE(INE(I,J))+1.
        END DO

        IF(PAN(I).EQ.3) THEN
        XI=1;eta=0.
        !CALL ISOPAR_3(xi,eta,corL(1:PAN(I),:),sn(1:PAN(I)),snx(1:PAN(I)),sne(1:PAN(I)),aj,sr,dn(1:PAN(I)))
        !IF(I.LE.5) WRITE(*,*)SN(1)
        VE(I,1,:)=DN(:)
        VECN(INE(I,1),:)=VECN(INE(I,1),:)+DN(:)/SN(1)
        XI=0;eta=1
        !CALL ISOPAR_3(xi,eta,corL(1:PAN(I),:),sn(1:PAN(I)),snx(1:PAN(I)),sne(1:PAN(I)),aj,sr,dn(1:PAN(I)))
        !IF(I.LE.5) WRITE(*,*)SN(1)
        VE(I,2,:)=DN(:)
        VECN(INE(I,2),:)=VECN(INE(I,2),:)+DN(:)/SN(2)
        XI=0;eta=0
        !CALL ISOPAR_3(xi,eta,corL(1:PAN(I),:),sn(1:PAN(I)),snx(1:PAN(I)),sne(1:PAN(I)),aj,sr,dn(1:PAN(I)))
        !IF(I.LE.5) WRITE(*,*)SN(1)
        VE(I,3,:)=DN(:)
        VECN(INE(I,3),:)=VECN(INE(I,3),:)+DN(:)/SN(3)
        END IF

        IF(PAN(I).EQ.4) THEN
        XI=-1;eta=-1
        CALL ISOPAR_4(xi,eta,corL(1:PAN(I),:),sn(1:PAN(I)),snx(1:PAN(I)),sne(1:PAN(I)),aj,sr,dn(1:PAN(I)))
        !IF(I.LE.5) WRITE(*,*)SN(1)
        VE(I,1,:)=DN(:)
        VECN(INE(I,1),:)=VECN(INE(I,1),:)+DN(:)/SN(1)
        XI=1;eta=-1
        CALL ISOPAR_4(xi,eta,corL(1:PAN(I),:),sn(1:PAN(I)),snx(1:PAN(I)),sne(1:PAN(I)),aj,sr,dn(1:PAN(I)))
        !IF(I.LE.5) WRITE(*,*)SN(2)
        VE(I,2,:)=DN(:)
        VECN(INE(I,2),:)=VECN(INE(I,2),:)+DN(:)/SN(2)
        XI=1;eta=1
        CALL ISOPAR_4(xi,eta,corL(1:PAN(I),:),sn(1:PAN(I)),snx(1:PAN(I)),sne(1:PAN(I)),aj,sr,dn(1:PAN(I)))
        !IF(I.LE.5) WRITE(*,*)SN(3)
        VE(I,3,:)=DN(:)
        VECN(INE(I,3),:)=VECN(INE(I,3),:)+DN(:)/SN(3)
        XI=-1;eta=1
        CALL ISOPAR_4(xi,eta,corL(1:PAN(I),:),sn(1:PAN(I)),snx(1:PAN(I)),sne(1:PAN(I)),aj,sr,dn(1:PAN(I)))
        !IF(I.LE.5) WRITE(*,*)SN(4)
        VE(I,4,:)=DN(:)
        VECN(INE(I,4),:)=VECN(INE(I,4),:)+DN(:)/SN(4)
        END IF

        IF(PAN(I).EQ.8) THEN
        XI=-D1;eta=-D3
        !WRITE(*,*)"OR",CORL(1,:)
        CALL ISOPAR_8(xi,eta,corL(1:8,:),sn(1:8),snx(1:8),sne(1:8),aj,sr,dn(1:8))
        !IF(I.LE.5) WRITE(*,*)SN(1)
        VE(I,1,:)=DN(:)
        VECN(INE(I,1),:)=VECN(INE(I,1),:)+DN(:)/SN(1)
        VECN(INE(I,1),:)=DN(:)/SN(1)
        XI=D2;eta=-D3
        !WRITE(*,*)"OR",CORL(2,:)
        CALL ISOPAR_8(xi,eta,corL(1:8,:),sn(1:8),snx(1:8),sne(1:8),aj,sr,dn(1:8))
        !IF(I.LE.5) WRITE(*,*)SN(2)
        VE(I,2,:)=DN(:)
        VECN(INE(I,2),:)=VECN(INE(I,2),:)+DN(:)/SN(2)
        VECN(INE(I,2),:)=DN(:)/SN(2)
        XI=D2;eta=D4
        !WRITE(*,*)"OR",CORL(3,:)
        CALL ISOPAR_8(xi,eta,corL(1:8,:),sn(1:8),snx(1:8),sne(1:8),aj,sr,dn(1:8))
        !IF(I.LE.5) WRITE(*,*)SN(3)
        VE(I,3,:)=DN(:)
        VECN(INE(I,3),:)=VECN(INE(I,3),:)+DN(:)/SN(3)
        VECN(INE(I,3),:)=DN(:)/SN(3)
        XI=-D1;eta=D4
        !WRITE(*,*)"OR",CORL(4,:)
        CALL ISOPAR_8(xi,eta,corL(1:8,:),sn(1:8),snx(1:8),sne(1:8),aj,sr,dn(1:8))
        !IF(I.LE.5) WRITE(*,*)SN(4)
        VE(I,4,:)=DN(:)
        VECN(INE(I,4),:)=VECN(INE(I,4),:)+DN(:)/SN(4)
        VECN(INE(I,4),:)=DN(:)/SN(4)
        XI=0.;eta=-D3
        !WRITE(*,*)"OR",CORL(5,:)
        CALL ISOPAR_8(xi,eta,corL(1:8,:),sn(1:8),snx(1:8),sne(1:8),aj,sr,dn(1:8))
        !IF(I.LE.5) WRITE(*,*)SN(5)
        VE(I,5,:)=DN(:)
        VECN(INE(I,5),:)=VECN(INE(I,5),:)+DN(:)/SN(5)
        VECN(INE(I,5),:)=DN(:)/SN(5)
        XI=D2;eta=0.
        !WRITE(*,*)"OR",CORL(6,:)
        CALL ISOPAR_8(xi,eta,corL(1:8,:),sn(1:8),snx(1:8),sne(1:8),aj,sr,dn(1:8))
        !IF(I.LE.5) WRITE(*,*)SN(6)
        VE(I,6,:)=DN(:)
        VECN(INE(I,6),:)=VECN(INE(I,6),:)+DN(:)/SN(6)
        VECN(INE(I,6),:)=DN(:)/SN(6)
        XI=0.;eta=D4
        CALL ISOPAR_8(xi,eta,corL(1:8,:),sn(1:8),snx(1:8),sne(1:8),aj,sr,dn(1:8))
        !IF(I.LE.5) WRITE(*,*)SN(7)
        VE(I,7,:)=DN(:)
        VECN(INE(I,7),:)=VECN(INE(I,7),:)+DN(:)/SN(7)
        VECN(INE(I,7),:)=DN(:)/SN(7)
        XI=-D1;eta=0.
        CALL ISOPAR_8(xi,eta,corL(1:8,:),sn(1:8),snx(1:8),sne(1:8),aj,sr,dn(1:8))
        !IF(I.LE.5) WRITE(*,*)SN(8)
        VE(I,8,:)=DN(:)
        VECN(INE(I,8),:)=VECN(INE(I,8),:)+DN(:)/SN(8)
        VECN(INE(I,8),:)=DN(:)/SN(8)
        20 CONTINUE
        END IF

        IF(PAN(I).EQ.9) THEN
        XI=-1;eta=-1
        CALL ISOPAR_9(xi,eta,corL,sn,snx,sne,aj,sr,dn)
        !IF(I.LE.5) WRITE(*,*)SN(1)
        VE(I,1,:)=DN(:)
        VECN(INE(I,1),:)=VECN(INE(I,1),:)+DN(:)/SN(1)
        !VECN(INE(I,1),:)=DN(:)/SN(1)
        XI=0;eta=-1
        CALL ISOPAR_9(xi,eta,corL,sn,snx,sne,aj,sr,dn)
        !IF(I.LE.5) WRITE(*,*)SN(2)
        VE(I,2,:)=DN(:)
        VECN(INE(I,2),:)=VECN(INE(I,2),:)+DN(:)/SN(2)
        !VECN(INE(I,2),:)=DN(:)/SN(2)
        XI=1;eta=-1
        CALL ISOPAR_9(xi,eta,corL,sn,snx,sne,aj,sr,dn)
        !IF(I.LE.5) WRITE(*,*)SN(3)
        VE(I,3,:)=DN(:)
        VECN(INE(I,3),:)=VECN(INE(I,3),:)+DN(:)/SN(3)
        !VECN(INE(I,3),:)=DN(:)/SN(3)
        XI=1;eta=0
        CALL ISOPAR_9(xi,eta,corL,sn,snx,sne,aj,sr,dn)
        !IF(I.LE.5) WRITE(*,*)SN(4)
        VE(I,4,:)=DN(:)
        VECN(INE(I,4),:)=VECN(INE(I,4),:)+DN(:)/SN(4)
        !VECN(INE(I,4),:)=DN(:)/SN(4)
        XI=1;eta=1
        CALL ISOPAR_9(xi,eta,corL,sn,snx,sne,aj,sr,dn)
        !IF(I.LE.5) WRITE(*,*)SN(5)
        VE(I,5,:)=DN(:)
        VECN(INE(I,5),:)=VECN(INE(I,5),:)+DN(:)/SN(5)
        !VECN(INE(I,5),:)=DN(:)/SN(5)
        XI=0;eta=1
        CALL ISOPAR_9(xi,eta,corL,sn,snx,sne,aj,sr,dn)
        !IF(I.LE.5) WRITE(*,*)SN(6)
        VE(I,6,:)=DN(:)
        VECN(INE(I,6),:)=VECN(INE(I,6),:)+DN(:)/SN(6)
        !VECN(INE(I,6),:)=DN(:)/SN(6)
        XI=-1;eta=1
        CALL ISOPAR_9(xi,eta,corL,sn,snx,sne,aj,sr,dn)
        !IF(I.LE.5) WRITE(*,*)SN(7)
        VE(I,7,:)=DN(:)
        VECN(INE(I,7),:)=VECN(INE(I,7),:)+DN(:)/SN(7)
        !VECN(INE(I,7),:)=DN(:)/SN(7)
        XI=-1;eta=0
        CALL ISOPAR_9(xi,eta,corL,sn,snx,sne,aj,sr,dn)
        !IF(I.LE.5) WRITE(*,*)SN(8)
        VE(I,8,:)=DN(:)
        VECN(INE(I,8),:)=VECN(INE(I,8),:)+DN(:)/SN(8)
        !VECN(INE(I,8),:)=DN(:)/SN(8)
        XI=0;eta=0
        CALL ISOPAR_9(xi,eta,corL,sn,snx,sne,aj,sr,dn)
        !IF(I.LE.5) WRITE(*,*)SN(9)
        VE(I,9,:)=DN(:)
        VECN(INE(I,9),:)=VECN(INE(I,9),:)+DN(:)/SN(9)
        !VECN(INE(I,9),:)=DN(:)/SN(9)
        END IF
    END DO

    !WRITE(10,101)NTPN,4*NTP
    !WRITE(11,101)NTPN,4*NTP
    WRITE(10,101)NBPOINT,TPNTP-4*NFBLOCK
    WRITE(11,101)NTPN,TPNTP !4*NTP
 101FORMAT(1x,'title="episode solid mesh"'/1x,'variables="x","Y",& 
     "Z","NX","NY","NZ"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')
 102FORMAT(6F15.6)

    DEALLOCATE(CORL)
    !CALL ELEMENTAREA

    DO I=1,NTPN
        !VECN(I,:)=VECN(I,:)/DNODE(I)
        VECN(I,:)=VECN(I,:)/SQRT(VECN(I,1)**2+VECN(I,2)**2+VECN(I,3)**2)
        !VECN(I,:)=VECN(I,:)
        !DMIP(I)=DMIP(I)/DNODE(I)
        !WRITE(*,*)DNODE(I)
    END DO


!   法向量除以模

    DO I=1,NTPN
        !IF(VECN(I,1).GT.VECN(I,2).AND.VECN(I,1).GT.VECN(I,3)
        !WRITE(10,102)CORD(I,:),VECN(I,:)
    END DO

    DO I=1,NBPOINT
        WRITE(10,102)CORD(I,:),VECN(I,:)
    END DO 

    DO I=1,NBBLOCK   
        IF(PAN(I).EQ.9) THEN
	    WRITE(10,*)INE(I,1),INE(I,2),INE(I,9),INE(I,8)
        WRITE(10,*)INE(I,8),INE(I,9),INE(I,6),INE(I,7)
        WRITE(10,*)INE(I,2),INE(I,3),INE(I,4),INE(I,9)
        WRITE(10,*)INE(I,9),INE(I,4),INE(I,5),INE(I,6)
        END IF
        IF(PAN(I).EQ.4.OR.PAN(I).EQ.8) THEN
	    WRITE(10,*)INE(I,1),INE(I,2),INE(I,3),INE(I,4)
        END IF
        IF(PAN(I).EQ.3) THEN
	    WRITE(10,*)INE(I,1),INE(I,2),INE(I,3),INE(I,3)
        END IF
    END DO

    DO I=1,NTPN
        IF(I.LE.NBPOINT) THEN
            WRITE(11,102)CORD(I,:),VECN(I,:)
        ELSE
            WRITE(11,102)CORD(I,:),0,0,1
        END IF
    END DO 
    DO I=1,NTP
        IF(PAN(I).EQ.9) THEN
	    WRITE(11,*)INE(I,1),INE(I,2),INE(I,9),INE(I,8)
        WRITE(11,*)INE(I,8),INE(I,9),INE(I,6),INE(I,7)
        WRITE(11,*)INE(I,2),INE(I,3),INE(I,4),INE(I,9)
        WRITE(11,*)INE(I,9),INE(I,4),INE(I,5),INE(I,6)
        END IF
        IF(PAN(I).EQ.4.OR.PAN(I).EQ.8) THEN
	    WRITE(11,*)INE(I,1),INE(I,2),INE(I,3),INE(I,4)
        END IF
        IF(PAN(I).EQ.3) THEN
	    WRITE(11,*)INE(I,1),INE(I,2),INE(I,3),INE(I,3)
        END IF
    END DO

    CLOSE(10)
END SUBROUTINE
  

SUBROUTINE ELEMENTAREA   !计算单元内潜或上置距离

	USE GREENMOD
	USE CUMOD
	USE NMLMOD

    !IMPLICIT NONE
    !WFA,VOL,XBUO(3)    
    COMMON /GWA/GA(400),GW(400)

    REAL:: XI,ETA,AJ
    REAL:: SN(9),SNX(9),SNE(9),SR(3),DN(3)
	REAL:: A(3,3),B(3,3),VV(3)
	REAL:: IS(3),JS(3)
    REAL:: TF1(6),TPF(6)
    REAL:: FCORL(9,3),FCNL(9,3)
    REAL:: TPA
	INTEGER L,K1
    INTEGER THR1,THR2,CP
	REAL:: T,D

    !WRITE(*,*) 
    !WRITE(*, *)'SUB ELEMENT AREA..............'

    open(18,file='coefuse.bin',form='binary',access='sequential',&
        status='unknown')
    rewind(18)
    
    open(118,file='coefuse.DAT')

	!ALLOCATE (PS(3),CORL(9,3),CNL(9,3),VE(NTP,4,3))
	ALLOCATE (HIJ(9),AIJ(9))

    K1=0

    IF(STAGE.EQ.2) THEN
        NPLOOP=NTPN
        NBLOOP=NTP
    END IF

    IF(STAGE.EQ.1) THEN
        NPLOOP=NBPOINT
        NBLOOP=NBBLOCK
    END IF

    DMIP=0.
    DMIE=0.
    DO IE=1,NTP !NBBLOCK
		NGUS=11
		K=NGUS*(NGUS-1)/2

		DO I1=1,PAN(IE)
			I2=I1
			DO J1=1,3
		       FCORL(I1,J1)=CORD(INE(IE,I1),J1)
			   FCNL(I1,J1)=VECN(INE(IE,I1),J1) 
            END DO
        END DO

!		GAUSS-LEGENDRE NUMERICAL INTEGRATION
!		对X积分     
		DO 42 IX=1,NGUS
		    XI=GA(K+IX-1)
!		    对Y积分
		    DO 42 IY=1,NGUS
			    ETA=GA(K+IY-1)
!			    调用ISOPAR，求得形状函数J及单位法向矢量等     
			    CALL ISOPAR_9(XI,ETA,FCORL,SN,SNX,SNE,AJ,SR,DN)

                TPA=0.		
			    DO 24 I1=1,9
                    TPA=TPA+SN(I1)
                    !WRITE(*,*)"TPF",TPF
   24		    CONTINUE
     		
                TPA=GW(K+IY-1)*AJ
                DMIE(IE)=DMIE(IE)+TPA*GW(K+IX-1)
   42		CONTINUE
        !DMIE(IE)=AJ

        DO I1=1,PAN(IE)
            !WRITE(*,*)DMIE(IE)
			K=INE(IE,I1)
            DMIP(K)=DMIP(K)+DMIE(IE)/4.
        END DO
    END DO

    !DEALLOCATE (PS,CORL,CNL)
	DEALLOCATE (HIJ,AIJ)

    RETURN
END SUBROUTINE