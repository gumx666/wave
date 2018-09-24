MODULE MESHDATA
    INTEGER:: NX,NZ,BX,NBZ,TPBU
    INTEGER:: NSX,STETY
    INTEGER:: NSEH,NSEB,NSES,NSEBS,NBX,NBSX
    INTEGER:: I,J,K
    INTEGER:: NUP,NLOW,NBLOW,NBUP
    INTEGER:: MSTEN,MBULB,MBULS
    INTEGER:: NWP
    INTEGER,DIMENSION(:):: NPOI(500),SYMME(3)
    REAL:: ZSTE,ZBUL,ZBOSS
    REAL:: LS,BS,TS
    REAL:: STENCUT,BULCUT1,BULCUT2
    REAL:: TRANS(3),SCALE
    REAL,ALLOCATABLE:: VERH(:,:,:)
    REAL,ALLOCATABLE:: VERZ(:,:,:)
    REAL,ALLOCATABLE:: VERI(:,:,:)
    REAL:: UI(0:200)
    REAL:: WX(200),WY(200),WZ(200)
    REAL:: WXVI(200),WYVI(200),WZVI(200)
    REAL:: WX1(200),WY1(200),WZ1(200)
    REAL:: SINK,TRI
    CHARACTER*50 TABFILE
END MODULE

SUBROUTINE SHIPHULL(SINK1,TRI1,NTP,WX3,WZ3,N,WX2,WY2,YTR,ZTR,TAOX,NTAO,ZFAC,WPVI1,NPVI1,WPVI2,NPVI2,NZ1,NWAPL,MARKWAP)
    USE MESHDATA
    REAL:: WX2(200),WY2(200)
    INTEGER NTP
    DIMENSION WX3(NTP),WZ3(NTP)
    DIMENSION MWAP(200)
    REAL:: TAO(100,3),YTR(100),ZTR(100),TAOX(100)
    REAL:: ZFAC
    REAL:: WPVI1(200,3),WPVI2(200,3)
    INTEGER,DIMENSION(:):: MARKWAP(100)
    INTEGER:: NWAPL
    INTEGER:: NPVI1,NPVI2

    
    NWP=NTP
    SINK=SINK1;TRI=TRI1

    WX=0;WY=0;WZ=0
    WX2=0.;WY2=0.
    !NWP=51
    !OPEN(20,FILE='WP.DAT')
    DO I=1,NWP
        WX(I)=WX3(I);WZ(I)=WZ3(I)
        WXVI(I)=WX3(I);WZVI(I)=WZ3(I)+ZFAC
        !WRITE(*,*)WX(I),WZ(I)
    !    READ(20,*)WX(I),WZ(I)
    END DO

    DO I=1,NWP
        WZ(I)=WZ(I)-SINK
        WZVI(I)=WZVI(I)-SINK
    END DO
    DO I=1,NWP
        TPX=WX(I);TPZ=WZ(I)
        WX(I)=TPX*COS(TRI)-TPZ*SIN(TRI)
        WZ(I)=TPX*SIN(TRI)+TPZ*COS(TRI)

        TPX=WXVI(I);TPZ=WZVI(I)
        WXVI(I)=TPX*COS(TRI)-TPZ*SIN(TRI)
        WZVI(I)=TPX*SIN(TRI)+TPZ*COS(TRI)
        !WRITE(1111,*)WX(I),WZ(I)
    END DO
    !CLOSE(1111)
!

    !WRITE(*,*)"«Î ‰»Î¥¨√˚£∫"
    !READ(*,*)TABFILE
    !TABFILE='S60'
    TABFILE='KCS'
    !TABFILE='KVLCC2'
    !TABFILE='5415'

    CALL READFILE 
    CALL READOFF 

    CALL CUTHULL_VI(WPVI2,NPVI2)
    
    DO I=1,NWP
        WXVI(I)=WX3(I);WZVI(I)=WZ3(I)+ZFAC/2.
    END DO
    DO I=1,NWP
        WZVI(I)=WZVI(I)-SINK
    END DO
    DO I=1,NWP
        TPX=WXVI(I);TPZ=WZVI(I)
        WXVI(I)=TPX*COS(TRI)-TPZ*SIN(TRI)
        WZVI(I)=TPX*SIN(TRI)+TPZ*COS(TRI)
    END DO
    CALL CUTHULL_VI(WPVI1,NPVI1)

    CALL CUTHULL_INI !CUT WETTED SURFACE
    CALL HULLMESH
    IF(MSTEN.GT.0) CALL STERNMESH
    IF(MBULB.GT.0) CALL BULBMESH
    IF(MBULS.GT.0) CALL BOSSMESH

    CALL MESHASSE
    
    IF(MSTEN.GT.0) THEN
    N=NX+NSX-1 
    DO I=1,NX
        WX2(N-I+1)=VERI(I,NZ,1)
        WY2(N-I+1)=VERI(I,NZ,2)
        !WRITE(*,*)WX2(N-I+1),WY2(N-I+1)
    END DO 
    DO I=NX+1,N
        WX2(N-I+1)=VERI(I+1,NUP,1)
        WY2(N-I+1)=VERI(I+1,NUP,2)
        !WRITE(*,*)WX2(N-I+1),WY2(N-I+1)
    END DO
    
    ELSE
    N=NX 
    DO I=1,N
        WX2(N-I+1)=VERI(I,NZ,1)
        WY2(N-I+1)=VERI(I,NZ,2)
        !WRITE(*,*)WX2(N-I+1),WY2(N-I+1)
    END DO 

    MWAP(I)=I*NZ
    END IF

    IF(MSTEN.GT.0) THEN
        I=NX+NSX
        I1=I-1
        DO J=1,NUP
            TAO(J,1)=(VERI(I,J,1)-VERI(I1,J,1))/&
            SQRT((VERI(I,J,1)-VERI(I1,J,1))**2+(VERI(I,J,2)-VERI(I1,J,2))**2+(VERI(I,J,3)-VERI(I1,J,3))**2)
            TAO(J,2)=(VERI(I,J,2)-VERI(I1,J,2))/&
            SQRT((VERI(I,J,1)-VERI(I1,J,1))**2+(VERI(I,J,2)-VERI(I1,J,2))**2+(VERI(I,J,3)-VERI(I1,J,3))**2)
            TAO(J,3)=(VERI(I,J,3)-VERI(I1,J,3))/&
            SQRT((VERI(I,J,1)-VERI(I1,J,1))**2+(VERI(I,J,2)-VERI(I1,J,2))**2+(VERI(I,J,3)-VERI(I1,J,3))**2)
            ZTR(J)=VERI(I,J,3)
            YTR(J)=VERI(I,J,2)
            TAOX(J)=TAO(J,1)
            TAOX(J)=(VERI(I,J,3)-VERI(I1,J,3))/(VERI(I,J,1)-VERI(I1,J,1))
        END DO
        NTAO=NUP
    ELSE
        I=NX+NSX
        I1=I-1
        DO J=1,NZ
            TAO(J,1)=(VERI(I,J,1)-VERI(I1,J,1))/&
            SQRT((VERI(I,J,1)-VERI(I1,J,1))**2+(VERI(I,J,2)-VERI(I1,J,2))**2+(VERI(I,J,3)-VERI(I1,J,3))**2)
            TAO(J,2)=(VERI(I,J,2)-VERI(I1,J,2))/&
            SQRT((VERI(I,J,1)-VERI(I1,J,1))**2+(VERI(I,J,2)-VERI(I1,J,2))**2+(VERI(I,J,3)-VERI(I1,J,3))**2)
            TAO(J,3)=(VERI(I,J,3)-VERI(I1,J,3))/&
            SQRT((VERI(I,J,1)-VERI(I1,J,1))**2+(VERI(I,J,2)-VERI(I1,J,2))**2+(VERI(I,J,3)-VERI(I1,J,3))**2)
            ZTR(J)=VERI(I,J,3)
            YTR(J)=VERI(I,J,2)
            TAOX(J)=TAO(J,1)
            TAOX(J)=(VERI(I,J,3)-VERI(I1,J,3))/(VERI(I,J,1)-VERI(I1,J,1))
            WRITE(*,*)YTR(J),ZTR(J)
        END DO 
        NTAO=NZ  
    END IF
    !NTAO=NUP
    DEALLOCATE(VERH,VERZ,VERI) 
    NZ1=NZ
    
    IF(MSTERN.EQ.0) THEN
        NWAPL=NX    
        DO I=1,NX
            MARKWAP(I)=I*NZ
        END DO    
    ELSE
        NWAPL=NX+NSX    
        DO I=1,NX
            MARKWAP(I)=I*NZ
        END DO  
        DO I=NX+1,NX+NSX
            MARKWAP(I)=NX*NZ+I*NUP
        END DO
    END IF
    
    !WRITE(*,*)N
    
    RETURN
END SUBROUTINE

SUBROUTINE READFILE    
    USE MESHDATA
    USE INPUTDATA
    CHARACTER*2 CTEM
    
    !OPEN(401,FILE="OFFSET_PARA.DAT")
    !OPEN(401,FILE=TRIM(TABFILE)//'_INPUT.DAT')

    !READ(401,*)CTEM !HULL
    !READ(401,*)LS,BS,TS
    LS=PANIN(1,1)
    BS=PANIN(1,2)
    TS=PANIN(1,3)
    
    !READ(401,*)
    !READ(401,*)CTEM !HULL
    !READ(401,*)NX,NZ
    NX=PANIN(2,1)
    NZ=PANIN(2,2)

    !READ(401,*)
    !READ(401,*)CTEM !STERN
    !READ(401,*)MSTEN,NSX,STETY,STENCUT
    MSTEN=PANIN(3,1)
    NSX=PANIN(3,2)
    STETY=PANIN(3,3)
    STENCUT=0.8
    
    !READ(401,*)    
    !READ(401,*)CTEM !BULB
    !READ(401,*)MBULB,NBX !,BULCUT1,BULCUT2
    BULCUT1=0.4;BULCUT2=0.4
    MBULB=PANIN(4,1)
    NBX=PANIN(4,2)
    
    !READ(401,*)
    !READ(401,*)CTEM !BULBSTERN
    !READ(401,*)MBULS,NBSX !,NBSZ
    MBULS=PANIN(5,1)
    NBSX=PANIN(5,2)

    !READ(401,*)
    !READ(401,*)CTEM !ZSTERN,ZBULB
    !READ(401,*)ZSTE,ZBUL,ZBOSS
    ZSTE=PANIN(6,1)
    ZBUL=PANIN(6,2)
    ZBOSS=PANIN(6,3)

    !READ(401,*)
    !READ(401,*)CTEM !TRANSLET
    !READ(401,*)TRANS(1:3)
    !WRITE(*,*)TRANS
    TRANS(1:3)=PANIN(7,1:3)

    !READ(401,*)
    !READ(401,*)CTEM !SCALE
    !READ(401,*)SCALE
    SCALE=PANIN(8,1)
    
    !READ(401,*)
    !READ(401,*)CTEM !SYMM
    !READ(401,*)SYMME
    SYMME(:)=PANIN(9,1:3)

    TS=TS/SCALE
    ZSTE=(ZSTE+TRANS(3))/SCALE !-TS
    ZBUL=(ZBUL+TRANS(3))/SCALE !-TS
    ZBOSS=(ZBOSS+TRANS(3))/SCALE !-TS
    !WRITE(*,*)ZSTE,ZBUL,ZBOSS
    CLOSE(401)
END SUBROUTINE

SUBROUTINE READOFF 
    USE MESHDATA
    USE FOLDERMOD
    CHARACTER*2 CTEM
    
    !OPEN(411,FILE="OFFSET_INPUT.DAT")
    !OPEN(411,FILE=TRIM(TABFILE)//'_OFFSET.DAT')
    OPEN(411,FILE=''//TRIM(ADJUSTL(ADDRESS1))//'\'//TRIM(ADJUSTL(FOLDERNAME1))//'_OFFIN.DAT')

    
    !READ(411,*)CTEM,NGR
    
    READ(411,*)CTEM,NSEH
    
    ALLOCATE(VERH(NSEH+500,200,3),VERZ(NSEH+500,NZ+50,3),VERI(NX+500,NZ+50,3))    
    
    DO I=1,NSEH
        READ(411,*)NPOI(I)
        !WRITE(*,*)NPOI(I)
        DO J=1,NPOI(I)
            READ(411,*)VERH(I,J,:)
        END DO
        !WRITE(*,*)VERH(I,1,1),NPOI(I)
    END DO   
    
    IF(MSTEN.GT.0) THEN
    READ(411,*)CTEM,NSES
    DO I=1,NSES
        READ(411,*)NPOI(I+NSEH)
        DO J=1,NPOI(I+NSEH)
            READ(411,*)VERH(I+NSEH,J,:)
        END DO
    END DO  
    END IF
         
    IF(MBULB.GT.0) THEN
    READ(411,*)CTEM,NSEB
    DO I=1,NSEB
        READ(411,*)NPOI(I+NSEH+NSES)
        DO J=1,NPOI(I+NSEH+NSES)
            READ(411,*)VERH(I+NSEH+NSES,J,:)
        END DO
    END DO  
    END IF

    IF(MBULS.GT.0) THEN
    READ(411,*)CTEM,NSEBS
    DO I=1,NSEBS
        READ(411,*)NPOI(I+NSEH+NSEB+NSES)
        DO J=1,NPOI(I+NSEH+NSEB+NSES)
            READ(411,*)VERH(I+NSEH+NSEB+NSES,J,:)
        END DO
    END DO  
    END IF

    DO I=1,NSEH+NSEB+NSES+NSEBS
        DO J=1,200
            VERH(I,J,:)=VERH(I,J,:)+TRANS(:)
        END DO
    END DO

    DO I=1,NSEH+NSEB+NSES+NSEBS
        DO J=1,200
            VERH(I,J,:)=VERH(I,J,:)/SCALE
            DO K=1,3
                VERH(I,J,K)=VERH(I,J,K)*SYMME(K)
            END DO
        END DO
    END DO
    
    CLOSE(411)
END SUBROUTINE

SUBROUTINE CUTHULL_VI(WPVI,NPVI)
    USE MESHDATA
 
    REAL:: UIA,UIB,UI0
    REAL:: C(3)
    REAL:: TP(100,3),X(200),Z(200) !,WX1(200),WZ1(200)
    REAL:: WPVI(200,3)
    INTEGER:: NPVI
    INTEGER:: COUT



!   STERN
    COUT=0
    NPVI=0
    IF(MSTEN.GT.0) THEN
    DO I=NSEH+1,NSEH+NSES
        IF(VERH(I,1,3).LT.WZVI(NWP)) THEN
        GOTO 601
        UIA=0.0
        UIB=1.0
        UI0=0.2   
        DO J=1,100
            CALL NURBS_INT_SI(NPOI(I)-1,VERH(I,1:NPOI(I),:),UI0,C)
            !WRITE(*,*)J,UI0
            WPVI(I,1)=C(1)
            WPVI(I,2)=C(2)
            WPVI(I,3)=C(3)

            IF(ABS(C(3)-WZVI(NWP)).LE.1.E-8) THEN
                NP=1
                DO K=1,NPOI(I)
                    IF(VERH(I,K,3).LE.WZVI(NWP)) THEN
                        NP=K
                    END IF
                END DO
                !NPOI(I)=NP+1
                
                !VERH(I,NPOI(I),:)=C(:)
                !VERH(I,NPOI(I),3)=0
                
                !WRITE(*,*)C(3),NPOI(I)
                WX1(NWP)=C(1)
                WY1(NWP)=C(2)
                WZ1(NWP)=C(3)

                EXIT
            ELSE IF(C(3)-WZVI(NWP).GT.0) THEN
                UIB=UI0
                UI0=(UIA+UIB)/2.
            ELSE IF(C(3)-WZVI(NWP).LT.0) THEN
                UIA=UI0
                UI0=(UIB+UIA)/2.
            END IF
        END DO
        601 CONTINUE
        ELSE
            !WRITE(*,*)"OK"
            COUT=COUT+1
        END IF
    END DO
    NPVI=NSES-COUT
    END IF

!   SHIP HULL
    DO I=1,1 !NSEH
        UIA=0.0
        UIB=1.0
        UI0=0.2   
        !WRITE(*,*)NPOI(I)
        DO J=1,100
            CALL NURBS_INT_SI(NPOI(I)-1,VERH(I,1:NPOI(I),:),UI0,C)
                !WRITE(*,*)J,UI0
            IF(ABS(C(3)-WZVI(1)).LE.1.E-8) THEN
                NP=1
                DO K=1,NPOI(I)
                    IF(VERH(I,K,3).LE.WZVI(1)) THEN
                        NP=K
                    END IF
                END DO
                !NPOI(I)=NP+1
                !NPOI(I)=NP !+1
                !VERH(I,NPOI(I),:)=C(:)
                !VERH(I,NPOI(I),3)=0
                
                WX1(1)=C(1)
                WY1(1)=C(2)
                WZ1(1)=C(3)
                EXIT
            ELSE IF(C(3)-WZVI(1).GT.0) THEN
                UIB=UI0
                UI0=(UIA+UIB)/2.
            ELSE IF(C(3)-WZVI(1).LT.0) THEN
                UIA=UI0
                UI0=(UIB+UIA)/2.
            END IF
        END DO
        
    END DO


    DO I=2,NWP-1
        X(I)=WXVI(NWP-I+1)
        Z(I)=WZVI(NWP-I+1)
    END DO
    X(1)=WX1(NWP);Z(1)=WZ1(NWP)
    X(NWP)=WX1(1);Z(NWP)=WZ1(1)

    !WRITE(*,*)Z(1:NWP)
    CALL ESPL2(X(1:NWP),Z(1:NWP),NWP,WX1(1:NSEH+NPVI),NSEH+NPVI,WZ1(1:NSEH+NPVI))
    !WX1(1)=WXVI(NWP)
    !WZ1(1)=WZVI(NWP)
    !WX1(NSEH+NSES)=WXVI(1)
    !WZ1(NSEH+NSES)=WZVI(1)
    !WRITE(*,*)
    !WRITE(*,*)WZ1(1:NSEH+NPVI)


    DO I=1,NSEH+NPVI
        UIA=0.0
        UIB=1.0
        UI0=0.2  
        !WRITE(*,*)NPOI(I) 
        DO J=1,100
            CALL NURBS_INT_SI(NPOI(I)-1,VERH(I,1:NPOI(I),:),UI0,C)
                !WRITE(*,*)J,UI0
            WPVI(I,1)=C(1)
            WPVI(I,2)=C(2)
            WPVI(I,3)=C(3)

            IF(ABS(C(3)-WZ1(I)).LE.1.E-8) THEN
                NP=1
                DO K=1,NPOI(I)
                    IF(VERH(I,K,3).LE.WZ1(I)) THEN
                        NP=K
                    END IF
                END DO
                !NPOI(I)=NP+1
                !NPOI(I)=NP+1
                !VERH(I,NPOI(I),:)=C(:)
                !WRITE(*,*)WPVI(I,3)
                !VERH(I,NPOI(I),3)=0
                EXIT
            ELSE IF(C(3)-WZ1(I).GT.0) THEN
                UIB=UI0
                UI0=(UIA+UIB)/2.
            ELSE IF(C(3)-WZ1(I).LT.0) THEN
                UIA=UI0
                UI0=(UIB+UIA)/2.
            END IF
        END DO
        
        IF(J.GE.100) THEN
            WPVI(I,1)=VERH(I,NPOI(I),1)
            WPVI(I,2)=VERH(I,NPOI(I),2)
            WPVI(I,3)=WZ1(I)
        END IF

        !WRITE(*,*)VERH(I,NPOI(I),3),J,NPOI(I)
        !WRITE(*,*)
    END DO
    NPVI=NPVI+NSEH

    DO I=1,NPVI
        WPVI(I,3)=WPVI(I,3)+SINK
    END DO

    DO I=1,NPVI
            TPX=WPVI(I,1);TPZ=WPVI(I,3)
            WPVI(I,1)=TPX*COS(TRI)+TPZ*SIN(TRI)
            WPVI(I,3)=-TPX*SIN(TRI)+TPZ*COS(TRI)
            !WRITE(*,*)WPVI(I,3)
    END DO
    !WRITE(*,*)
END SUBROUTINE

SUBROUTINE CUTHULL_INI
    USE MESHDATA
 
    REAL:: UIA,UIB,UI0
    REAL:: C(3)
    REAL:: TP(100,3),X(200),Z(200)
    INTEGER:: COUT

!   STERN
    COUT=0
    IF(MSTEN.GT.0) THEN
    DO I=NSEH+1,NSEH+NSES
        IF(VERH(I,1,3).LT.WZ(NWP)) THEN
        GOTO 601
        UIA=0.0
        UIB=1.0
        UI0=0.2   
        DO J=1,100
            CALL NURBS_INT_SI(NPOI(I)-1,VERH(I,1:NPOI(I),:),UI0,C)
            !WRITE(*,*)J,UI0
            IF(ABS(C(3)-WZ(NWP)).LE.1.E-8) THEN
                NP=1
                DO K=1,NPOI(I)
                    IF(VERH(I,K,3).LE.WZ(NWP)) THEN
                        NP=K
                    END IF
                END DO
                NPOI(I)=NP !+1
                VERH(I,NPOI(I),:)=C(:)
                VERH(I,NPOI(I),3)=0
                !WRITE(*,*)C(3),NPOI(I)
                EXIT
            ELSE IF(C(3)-WZ(NWP).GT.0) THEN
                UIB=UI0
                UI0=(UIA+UIB)/2.
            ELSE IF(C(3)-WZ(NWP).LT.0) THEN
                UIA=UI0
                UI0=(UIB+UIA)/2.
            END IF
        END DO
        601 CONTINUE
        ELSE
            !WRITE(*,*)"OK"
            COUT=COUT+1
        END IF
    END DO

    IF(COUT.EQ.0) GOTO 901

    NLS=101
    DO J=1,NLS
        UI(J)=1./(NLS-1)*(J-1)
    END DO
    CALL NURBS_INT(NSES-1,VERH(NSEH+1:NSEH+NSES,1,:),NLS-1,UI(1:NLS),VERH(NSEH+NSES-COUT+1,1:NLS,:))
    JA=0
    JB=0
    DO I=1,NLS
        IF(VERH(NSEH+NSES-COUT+1,I,1).LT.VERH(NSEH+NSES-COUT,1,1)) THEN
            JA=I
            EXIT
        END IF
    END DO
    DO I=1,NLS
        !WRITE(*,*)I,VERH(NSEH+NSES-COUT+1,I,:)
        IF(VERH(NSEH+NSES-COUT+1,I,3).LT.WZ(NWP)) THEN
            !WRITE(*,*)VERH(NSEH+NSES-COUT+1,I,3)
            JB=I
        END IF
    END DO
    
    !WRITE(*,*)VERH(NSEH+NSES-COUT,1,3),VERH(NSEH+NSES-COUT+1,JA,3)
    WRITE(*,*)JA,JB
    IF(JA.EQ.JB) THEN
        COUT=COUT+1
        GOTO 991
    END IF

    IF(ABS(VERH(NSEH+NSES-COUT+1,JB,1)-VERH(NSEH+NSES-COUT,NPOI(NSEH+NSES-COUT),1))/&
       ABS(VERH(NSEH+NSES-COUT+1,JB,2)-VERH(NSEH+NSES-COUT,NPOI(NSEH+NSES-COUT),2)).LE.TAN(20./180.*3.14)) THEN
       COUT=COUT+1
       GOTO 991
    END IF

    DO I=1,JB-JA+1
        TP(I,:)=VERH(NSEH+NSES-COUT+1,I+JA-1,:)
    END DO
    DO I=1,JB-JA+1
        VERH(NSEH+NSES-COUT+1,I,:)=TP(I,:)
    END DO

    VERH(NSEH+NSES-COUT+1,JB-JA+2,3)=WZ(NWP)
    VERH(NSEH+NSES-COUT+1,JB-JA+2,2)=0.
    VERH(NSEH+NSES-COUT+1,JB-JA+2,1)=VERH(NSEH+NSES-COUT+1,JB-JA+1,1)+&
                                    (VERH(NSEH+NSES-COUT+1,JB-JA+1,1)-VERH(NSEH+NSES-COUT+1,JB-JA,1))/&
                                    (VERH(NSEH+NSES-COUT+1,JB-JA+1,3)-VERH(NSEH+NSES-COUT+1,JB-JA,3))*&
                                    (VERH(NSEH+NSES-COUT+1,JB-JA+2,3)-VERH(NSEH+NSES-COUT+1,JB-JA+1,3))
    NPOI(NSEH+NSES-COUT+1)=JB-JA+2

    DO I=1,JB-JA+2
        !WRITE(*,*)VERH(NSEH+NSES-COUT+1,I,:)
    END DO

    991 CONTINUE

    COUT=COUT-1
!   CHANGE SECTIONS
    DO I=NSEH+NSES-COUT+1,NSEH+NSES+NSEB+NSEBS
        VERH(I,1:200,:)=VERH(I+COUT,1:200,:)
        NPOI(I)=NPOI(I+COUT)
    END DO
    NSES=NSES-COUT

    901 CONTINUE

!   LAST POINT ON WAVE PROFIILE
    WX1(NSEH+NSES)=VERH(NSEH+NSES,NPOI(NSEH+NSES),1)
    WY1(NSEH+NSES)=VERH(NSEH+NSES,NPOI(NSEH+NSES),2)
    !WZ1(NSEH+NSES)=VERH(NSEH+NSES,NPOI(NSEH+NSES),3)

    END IF

    !WRITE(*,*)NSES

!   SHIP HULL
    DO I=1,1 !NSEH
        UIA=0.0
        UIB=1.0
        UI0=0.2   
        !WRITE(*,*)NPOI(I)
        DO J=1,100
            CALL NURBS_INT_SI(NPOI(I)-1,VERH(I,1:NPOI(I),:),UI0,C)
                !WRITE(*,*)J,UI0
            IF(ABS(C(3)-WZ(1)).LE.1.E-8) THEN
                NP=1
                DO K=1,NPOI(I)
                    IF(VERH(I,K,3).LE.WZ(1)) THEN
                        NP=K
                    END IF
                END DO
                !NPOI(I)=NP+1
                NPOI(I)=NP !+1
                VERH(I,NPOI(I),:)=C(:)
                !VERH(I,NPOI(I),3)=0
                EXIT
            ELSE IF(C(3)-WZ(1).GT.0) THEN
                UIB=UI0
                UI0=(UIA+UIB)/2.
            ELSE IF(C(3)-WZ(1).LT.0) THEN
                UIA=UI0
                UI0=(UIB+UIA)/2.
            END IF
        END DO
        !WRITE(*,*)VERH(I,1:NPOI(I),3)
        !WRITE(*,*)
        !WRITE(*,*)VERH(I,NPOI(I),3),VERH(I,NPOI(I)-1,3)
        !WRITE(*,*)
        
        WX1(1)=VERH(1,NPOI(1),1)
        WY1(1)=VERH(1,NPOI(1),2)
        WZ1(1)=VERH(1,NPOI(1),3)
    END DO

!   CHANGE WATER LINE
    DO I=2,NWP-1
        WX(I)=(WX1(NSEH+NSES)-WX1(1))/(WX(NWP)-WX(1))*WX(I)
        !WY(I)=(-WX1(NSEH+NSES)+WX1(1))/(WX(NWP)-WX(1))*WY(I)
        !WZ(I)=(-WX1(NSEH+NSES)+WX1(1))/(WX(NWP)-WX(1))*WZ(I)
        !WRITE(*,*)WX(I)
    END DO
    WX(1)=VERH(1,NPOI(1),1) !;WY(1)=VERH(1,NPOI(1),2);WZ(1)=VERH(1,NPOI(1),3)
    WX(NWP)=VERH(NSEH+NSES,NPOI(NSEH+NSES),1) !;WY(NWP)=VERH(1,NPOI(NSEH+NSES),2);WZ(NWP)=VERH(1,NPOI(NSEH+NSES),3)

    DO I=1,NSEH+NSES
        WX1(I)=VERH(I,NPOI(I),1);WY1(I)=VERH(I,NPOI(I),2)
        !WRITE(*,*)WX(1:NWP)
    END DO
    DO I=1,NWP
        X(I)=WX(NWP-I+1)
        Z(I)=WZ(NWP-I+1)
    END DO
    !WRITE(*,*)Z(1:NWP)

    CALL ESPL2(X(1:NWP),Z(1:NWP),NWP,WX1(1:NSEH+NSES),NSEH+NSES,WZ1(1:NSEH+NSES))
    DO I=2,NSEH+NSES-1,2
        WX1(I)=(WX1(I+1)+WX1(I-1))/2.
        !WRITE(*,*)WZ1(I)
    END DO
    !WZ1(NSEH+NSES)=WZ(NWP)
    !WRITE(*,*)
    !WRITE(*,*)VERH(1,NPOI(1),1),VERH(NSEH+NSES,NPOI(NSEH+NSES),1)
    !WRITE(*,*)WZ1(1:NSEH+NSES)

!
    DO I=2,NSEH+NSES
        UIA=0.0
        UIB=1.0
        UI0=0.2  
        !WRITE(*,*)NPOI(I) 
        DO J=1,100
            CALL NURBS_INT_SI(NPOI(I)-1,VERH(I,1:NPOI(I),:),UI0,C)
                !WRITE(*,*)J,UI0
            IF(ABS(C(3)-WZ1(I)).LE.1.E-8) THEN
                NP=1
                DO K=1,NPOI(I)
                    IF(VERH(I,K,3).LE.WZ1(I)) THEN
                        NP=K
                    END IF
                END DO
                !NPOI(I)=NP+1
                NPOI(I)=NP !+1
                VERH(I,NPOI(I),:)=C(:)
                !VERH(I,NPOI(I),3)=0
                EXIT
            ELSE IF(C(3)-WZ1(I).GT.0) THEN
                UIB=UI0
                UI0=(UIA+UIB)/2.
            ELSE IF(C(3)-WZ1(I).LT.0) THEN
                UIA=UI0
                UI0=(UIB+UIA)/2.
            END IF
        END DO
        
        !WRITE(*,*)VERH(I,NPOI(I),3),J,NPOI(I)
        !WRITE(*,*)
    END DO
END SUBROUTINE

SUBROUTINE CUTHULL
    USE MESHDATA
 
    REAL:: UIA,UIB,UI0
    REAL:: C(3)


!   FIRST SECTION
    UIA=0.3
    UIB=1.0
    UI0=0.3   
    DO I=1,100
        CALL NURBS_INT_SI(NPOI(1)-1,VERH(1,1:NPOI(1),:),UI0,C)
        IF(ABS(C(3)-WZ(NWP)).LE.1.E-7) THEN
            NP=1
            DO J=1,NPOI(1)
                IF(VERH(1,J,3).LE.C(3)) THEN
                    NP=J
                END IF
            END DO
            NPOI(1)=NP+1
            VERH(1,NPOI(1),:)=C(:)
            EXIT
        ELSE IF(C(3).GT.WZ(NWP)) THEN
            UI0=(UIA+UI0)/2.
        ELSE IF(C(3).LT.WZ(NWP)) THEN
            UI0=(UIB+UI0)/2.
        END IF
    END DO

!   LAST SECTION
    IF(MSTEN.EQ.0) THEN
    UIA=0.3
    UIB=1.0
    UI0=0.3   
    DO I=1,100
        CALL NURBS_INT_SI(NPOI(NSEH)-1,VERH(NSEH,1:NPOI(NSEH),:),UI0,C)
        IF(ABS(C(3)-WZ(NWP)).LE.1.E-7) THEN
            NP=1
            DO J=1,NPOI(NSEH)
                IF(VERH(NSEH,J,3).LE.C(3)) THEN
                    NP=J
                END IF
            END DO
            NPOI(NSEH)=NP+1
            VERH(NSEH,NPOI(NSEH),:)=C(:)
            EXIT
        ELSE IF(C(3).GT.WZ(NWP)) THEN
            UI0=(UIA+UI0)/2.
        ELSE IF(C(3).LT.WZ(NWP)) THEN
            UI0=(UIB+UI0)/2.
        END IF
    END DO
    END IF
END SUBROUTINE

SUBROUTINE BOSSMESH
    USE MESHDATA

!   Z DIRECTION
    DO J=1,NLOW
        UI(J)=1./(NLOW-1)*(J-1)
    END DO
    
    !WRITE(*,*)NLOW
    DO I=NSEH+NSES+NSEB+1,NSEH+NSES+NSEB+NSEBS
        CALL NURBS_INT(NPOI(I)-1,VERH(I,1:NPOI(I),:),NLOW-1,UI(1:NLOW),VERZ(I,1:NLOW,:))
        !WRITE(*,*)NPOI(I),VERH(I,1:NPOI(I),1)
        !WRITE(*,*)VERZ(I,1:NLOW,3)
        !WRITE(*,*)
    END DO
!
!   X DIRECTION
    DO J=1,NBSX
        UI(J)=1./(NBSX-1)*(J-1)
    END DO
    !WRITE(*,*)NSX,NSES
    
    !CHANGE LAST SECTION OH HULL
    VERZ(NSEH+NSES+NSEB+1,1:NLOW,:)=VERI(NX,1:NLOW,:)
    !

    DO I=1,NLOW
        CALL NURBS_INT(NSEBS-1,VERZ(NSEH+NSES+NSEB+1:NSEH+NSES+NSEB+NSEBS,I,:),NBSX-1,UI(1:NBSX),VERI(NX+NSX+NBX+1:NX+NSX+NBX+NBSX,I,:))
        !WRITE(*,*)VERZ(NSEH+NSES+NSEB+1:NSEH+NSES+NSEB+NSEBS,I,1)
        !WRITE(*,*)
    END DO

    
END SUBROUTINE

SUBROUTINE BULBMESH
    USE MESHDATA

    INTEGER:: NST1,NST2
        NBX0=NBX
        NBX=41
    
!   Z DIRECTION
    DO J=1,NBLOW
        UI(J)=1./(NBLOW-1)*(J-1)
    END DO

    DO I=NSEH+NSES+1,NSEH+NSES+NSEB        
        CALL NURBS_INT(NPOI(I)-1,VERH(I,1:NPOI(I),:),NBLOW-1,UI(1:NUP),VERZ(I,1:NBLOW,:))
        !WRITE(*,*)NPOI(I),VERH(I,1:NPOI(I),1)
        !WRITE(*,*)VERZ(I,1:NBLOW,1)
        !WRITE(*,*)
    END DO
!
!   X DIRECTION
    DO J=1,NBX
        UI(J)=1./(NBX-1)*(J-1)
    END DO
    !WRITE(*,*)NSX,NSES
    
    !CHANGE FIRST SECTION OH HULL
    VERZ(NSEH+NSES+NSEB+1,1:NBLOW,:)=VERI(1,1:NBLOW,:)
    !WRITE(*,*)NBX-1,NSEB
    !
    DO I=1,NBLOW
        CALL NURBS_INT(NSEB,VERZ(NSEH+NSES+1:NSEH+NSES+NSEB+1,I,:),NBX-1,UI(1:NBX),VERI(NX+NSX+1:NX+NSX+NBX,I,:))
        !WRITE(*,*)VERI(NX+NSX+1:NX+NSX+NBX,I,1)
        !WRITE(*,*)
    END DO

    !GOTO 50
!***********************************************************
        !CUT LINE
        NST1=INT(NBX*BULCUT1)
        NST2=INT(NBX*BULCUT2)
        NST0=NST1
        IF(NST2.GE.NST1) NST0=NST2
        !WRITE(*,*)NST1,NST2,NST0

        DO I=2,NST1
            VERZ(NX+NSX+NBX,NST1-I+1,:)=VERI(NX+NSX+I,1,:)
            !WRITE(*,*)VERI(NX+NSX+I,1,1)
        END DO
        !WRITE(*,*)
        DO I=1,NST2
            VERZ(NX+NSX+NBX,I+NST1-1,:)=VERI(NX+NSX+I,NBLOW,:)
            !WRITE(*,*)VERI(NX+NSX+NST1-I+1,NBLOW,1)
        END DO
        VERI(NX+NSX+NST0,1:NST1+NST2-1,:)=VERZ(NX+NSX+NBX,1:NST1+NST2-1,:)
        !WRITE(*,*)VERI(NX+NSX+NST0,1:NST1+NST2-1,1)

        DO J=1,NBLOW
            UI(J)=1./(NBLOW-1)*(J-1)
        END DO
        CALL NURBS_INT(NST1+NST2-2,VERI(NX+NSX+NST0,1:NST1+NST2-1,:),NBLOW-1,UI(1:NBLOW),VERI(NX+NSX+NST0,1:NBLOW,:))
        !WRITE(*,*)VERI(NX+NSX+NST0,1:NBLOW,:)
        
        DO J=1,NBX0
            !WRITE(*,*)NBX0
            UI(J)=1./(NBX0-1)*(J-1)
        END DO
        DO I=1,NBLOW
            CALL NURBS_INT(NBX-NST0,VERI(NX+NSX+NST0:NX+NSX+NBX,I,:),NBX0-1,UI(1:NBX0),VERI(NX+NSX+1:NX+NSX+NBX0,I,:))
            !WRITE(*,*)VERI(NX+NSX+NST0:NX+NSX+NBX,I,:)
        END DO
        NBX=NBX0
    50 CONTINUE

END SUBROUTINE

SUBROUTINE STERNMESH
    USE MESHDATA

    
    IF(STETY.EQ.0.AND.STENCUT.GT.0) THEN
        NSX0=NSX
        !NSX=40
    END IF
    
!   NO STERN
!   Z DIRECTION
    DO J=1,NUP
        UI(J)=1./(NUP-1)*(J-1)
    END DO

    DO I=NSEH+1,NSEH+NSES
        CALL NURBS_INT(NPOI(I)-1,VERH(I,1:NPOI(I),:),NUP-1,UI(1:NUP),VERZ(I,1:NUP,:))
        !WRITE(*,*)I,VERH(I,1:NPOI(I),3)
        !WRITE(*,*)
    END DO

!
!   X DIRECTION
    DO J=1,NSX
        UI(J)=1./(NSX-1)*(J-1)
    END DO
    !WRITE(*,*)NSX,NSES
    
    !CHANGE LAST SECTION OH HULL
    VERZ(1,1:NZ,:)=VERZ(NSEH,1:NZ,:)
    !VERZ(NSEH,1:NUP,:)=VERZ(1,NLOW:NZ,:)
    VERZ(NSEH,1:NUP,:)=VERZ(1,NZ-NUP+1:NZ,:)

    !
    DO I=1,NUP
        CALL NURBS_INT(NSES,VERZ(NSEH:NSEH+NSES,I,:),NSX-1,UI(1:NSX),VERI(NX+1:NX+NSX,I,:))
    END DO

    GOTO 50
    
    IF(STETY.EQ.0.AND.STENCUT.GT.0) THEN
        !CUT LINE
        NST0=INT(NSX*STENCUT)

    
        DO J=1,NUP
            UI(J)=1./(NUP-1)*(J-1)
        END DO
        CALL NURBS_INT(NSX-NST0,VERI(NX+NST0:NX+NSX,1,:),NUP-1,UI(1:NUP),VERI(NX+NST0,1:NUP,:))
        
        
        DO J=1,NSX0
            UI(J)=1./(NSX0-1)*(J-1)
        END DO
        DO I=1,NUP
            CALL NURBS_INT(NST0-1,VERI(NX+1:NX+NST0,I,:),NSX0-1,UI(1:NSX0),VERI(NX+1:NX+NSX0,I,:))
        END DO
        NSX=NSX0
    END IF
    50 CONTINUE
END SUBROUTINE

SUBROUTINE HULLMESH
    USE MESHDATA
    REAL TPVER(200,3)

    IF(MSTEN.EQ.0.AND.MBULB.EQ.0) THEN
!   NO STERN
!   Z DIRECTION
    DO J=1,NZ
        UI(J)=1./(NZ-1)*(J-1)
    END DO

    DO I=1,NSEH
        CALL NURBS_INT(NPOI(I)-1,VERH(I,1:NPOI(I),:),NZ-1,UI(1:NZ),VERZ(I,1:NZ,:))
        DO J=1,NZ
            IF(VERZ(I,J,3).LT.-TS) THEN
                !VERZ(I,J,3)=-TS
            END IF
        END DO
        !WRITE(*,*)
        !WRITE(*,*)VERH(I,1,1),NPOI(I)
        !WRITE(*,*)VERH(I,1:NPOI(I),3)
        !WRITE(*,*)I
        !WRITE(*,*)VERZ(I,1:NZ,1)
        !STOP
    END DO
!
!   X DIRECTION
    DO J=1,NX
        UI(J)=1./(NX-1)*(J-1)
    END DO

    DO I=1,NZ
        CALL NURBS_INT(NSEH-1,VERZ(1:NSEH,I,:),NX-1,UI(1:NX),VERI(1:NX,I,:))
        !WRITE(*,*)VERZ(1:NSEH,I,1)
        !WRITE(*,*)
        !WRITE(*,*)I
        !WRITE(*,*)VERI(1:NX,I,1)
    END DO
    END IF
!**********************************************

    IF(MSTEN.GT.0.AND.MBULB.EQ.0) THEN
!   ONLY STERN
    NST0=1
    DO J=2,NPOI(NSEH)
        IF(ABS(VERH(NSEH,NST0,3)-ZSTE).GE.ABS(VERH(NSEH,J,3)-ZSTE)) THEN
            NST0=J
        END IF
    END DO

    D1=0.
    DO J=1,NST0-1
        D1=D1+SQRT((VERH(NSEH,J,1)-VERH(NSEH,J+1,1))**2+&
                   (VERH(NSEH,J,2)-VERH(NSEH,J+1,2))**2+&
                   (VERH(NSEH,J,3)-VERH(NSEH,J+1,3))**2)
    END DO

    D2=0.
    DO J=NST0,NPOI(NSEH)-1
        D2=D2+SQRT((VERH(NSEH,J,1)-VERH(NSEH,J+1,1))**2+&
                   (VERH(NSEH,J,2)-VERH(NSEH,J+1,2))**2+&
                   (VERH(NSEH,J,3)-VERH(NSEH,J+1,3))**2)
    END DO

    NLOW=INT(D1/(D1+D2)*(NZ-1))+1
    NUP=NZ-NLOW+1
    IF(MOD(NUP,2).EQ.0) THEN
        NUP=NUP-1
        NLOW=NLOW+1
    END IF

    !WRITE(*,*)NST0,NLOW,NUP

!   Z DIRECTION
    DO J=1,NZ
        UI(J)=1./(NZ-1)*(J-1)
    END DO

    DO I=1,NSEH-1
        CALL NURBS_INT(NPOI(I)-1,VERH(I,1:NPOI(I),:),NZ-1,UI(1:NZ),VERZ(I,1:NZ,:))
        DO J=1,NZ
            IF(VERZ(I,J,3).LT.-TS) THEN
                !VERZ(I,J,3)=-TS
            END IF
        END DO
    END DO

    !LOWER
    DO J=1,NLOW
        UI(J)=1./(NLOW-1)*(J-1)
    END DO
    DO I=NSEH,NSEH
        CALL NURBS_INT(NST0-1,VERH(I,1:NST0,:),NLOW-1,UI(1:NLOW),VERZ(I,1:NLOW,:))
    END DO
    
    !UPPER
    DO J=1,NUP
        UI(J)=1./(NUP-1)*(J-1)
    END DO

    DO I=NSEH,NSEH
        CALL NURBS_INT(NPOI(I)-NST0,VERH(I,NST0:NPOI(I),:),NZ-NLOW,UI(1:NUP),VERZ(I,NLOW:NZ,:))
    END DO
!
!   X DIRECTION
    DO J=1,NX
        UI(J)=1./(NX-1)*(J-1)
    END DO

    DO I=1,NZ
        CALL NURBS_INT(NSEH-1,VERZ(1:NSEH,I,:),NX-1,UI(1:NX),VERI(1:NX,I,:))
    END DO
    !VERI(NSEH,NLOW,2)=0.

    END IF


!*************************************************************************************************

!   STERN AND BULB
    IF(MSTEN.GT.0.AND.MBULB.GT.0) THEN
!   CUT STERN
    NST0S=1
    DO J=2,NPOI(NSEH)
        IF(VERH(NSEH,J,3).GE.ZSTE) THEN
            NST0S=J
            EXIT
        END IF
    END DO

    NST0B=1
    DO J=2,NPOI(NSEH)
        IF(VERH(NSEH,J,3).LE.ZBOSS) THEN
            NST0B=J
        END IF
    END DO

    D1=0.
    DO J=1,NST0B-1
        D1=D1+SQRT((VERH(NSEH,J,1)-VERH(NSEH,J+1,1))**2+&
                   (VERH(NSEH,J,2)-VERH(NSEH,J+1,2))**2+&
                   (VERH(NSEH,J,3)-VERH(NSEH,J+1,3))**2)
    END DO
    D2=0.
    DO J=NST0S,NPOI(NSEH)-1
        D2=D2+SQRT((VERH(NSEH,J,1)-VERH(NSEH,J+1,1))**2+&
                   (VERH(NSEH,J,2)-VERH(NSEH,J+1,2))**2+&
                   (VERH(NSEH,J,3)-VERH(NSEH,J+1,3))**2)
    END DO
    D12=0.
    DO J=NST0B,NST0S-1
        D12=D12+SQRT((VERH(NSEH,J,1)-VERH(NSEH,J+1,1))**2+&
                     (VERH(NSEH,J,2)-VERH(NSEH,J+1,2))**2+&
                     (VERH(NSEH,J,3)-VERH(NSEH,J+1,3))**2)
    END DO

    NLOW=INT(D1/(D1+D2+D12)*(NZ-1))+1
    NUP=INT(D2/(D1+D2+D12)*(NZ-1))+1
    IF(MOD(NUP,2).EQ.0) THEN
        NUP=NUP+1
    END IF
    IF(MOD(NLOW,2).EQ.0) THEN
        NLOW=NLOW-1
    END IF
    NIN=NZ-NLOW-NUP+2
    
!   CUT BULB
    NSTB=1
    DO J=2,NPOI(1)
        !WRITE(*,*)ABS(VERH(1,NSTB,3)-ZBUL),ABS(VERH(1,J,3)-ZBUL),NSTB,NPOI(1)
        !IF(ABS(VERH(1,NSTB,3)-ZBUL).GE.ABS(VERH(1,J,3)-ZBUL)) THEN
        IF(VERH(1,J,3).GE.ZBUL) THEN
            NSTB=J-1
            EXIT
        END IF
    END DO
    D1=0.
    DO J=1,NSTB-1
        D1=D1+SQRT((VERH(1,J,1)-VERH(1,J+1,1))**2+&
                   (VERH(1,J,2)-VERH(1,J+1,2))**2+&
                   (VERH(1,J,3)-VERH(1,J+1,3))**2)
    END DO
    D2=0.
    DO J=NSTB,NPOI(1)-1
        D2=D2+SQRT((VERH(1,J,1)-VERH(1,J+1,1))**2+&
                   (VERH(1,J,2)-VERH(1,J+1,2))**2+&
                   (VERH(1,J,3)-VERH(1,J+1,3))**2)
    END DO
    NBLOW=INT(D1/(D1+D2)*(NZ-1))+1
    NBUP=NZ-NBLOW+1
    IF(MOD(NBLOW,2).EQ.0) THEN
        NBUP=NBUP-1
        NBLOW=NBLOW+1
    END IF

    !WRITE(*,*)NLOW,NUP,NZ
    !STOP

!   Z DIRECTION
    DO J=1,NZ
        UI(J)=1./(NZ-1)*(J-1)
    END DO

    DO I=2,NSEH-1
        CALL NURBS_INT(NPOI(I)-1,VERH(I,1:NPOI(I),:),NZ-1,UI(1:NZ),VERZ(I,1:NZ,:))
        DO J=1,NZ
            IF(VERZ(I,J,3).LT.-TS) THEN
                VERZ(I,J,3)=-TS
            END IF
        END DO
        !WRITE(*,*)VERZ(I,1:NZ,3)
    END DO

    
!   LAST SECTION
    !LOWER
    DO J=1,NLOW
        UI(J)=1./(NLOW-1)*(J-1)
    END DO
    TPVER(1:NST0B,:)=VERH(NSEH,1:NST0B,:)
    TPVER(NST0B+1,3)=ZBOSS;TPVER(NST0B+1,1)=VERH(NSEH,NST0B,1);TPVER(NST0B+1,2)=0.
    DO I=NSEH,NSEH
        CALL NURBS_INT(NST0B,TPVER(1:NST0B+1,:),NLOW-1,UI(1:NLOW),VERZ(I,1:NLOW,:))
        !WRITE(*,*)TPVER(1:NST0B+1,1)
    END DO

    !UPPER
    DO J=1,NUP
        UI(J)=1./(NUP-1)*(J-1)
    END DO
    TPVER(2:NPOI(NSEH)-NST0S+2,:)=VERH(NSEH,NST0S:NPOI(NSEH),:)
    TPVER(1,3)=ZSTE;TPVER(1,1)=VERH(NSEH,NST0S,1);TPVER(1,2)=0. !VERH(NSEH,NST0S,2)
    DO I=NSEH,NSEH
        !WRITE(*,*)NZ-NUP+1,NZ
        !WRITE(*,*)TPVER(1:NPOI(NSEH)-NST0S+2,3)
        CALL NURBS_INT(NPOI(I)-NST0S+1,TPVER(1:NPOI(NSEH)-NST0S+2,:),NUP-1,UI(1:NUP),VERZ(I,NZ-NUP+1:NZ,:))
        !WRITE(*,*)VERZ(I,NZ-NUP+1:NZ,2)
    END DO
    
    
    !WRITE(*,*)NST0S,NST0B
    !INTER
    DO J=1,NIN
        UI(J)=1./(NIN-1)*(J-1)
    END DO
    TPVER(NST0B+1:NST0S-1,:)=VERH(NSEH,NST0B+1:NST0S-1,:)
    TPVER(NST0B+1:NST0S-1,2)=0
    TPVER(NST0B,3)=ZBOSS;TPVER(NST0B,1)=VERH(NSEH,NST0B,1);TPVER(NST0B,2)=0.
    TPVER(NST0S,3)=ZSTE;TPVER(NST0S,1)=VERH(NSEH,NST0S,1);TPVER(NST0S,2)=0. !VERH(NSEH,NST0S,2)
    DO I=NSEH,NSEH
        !WRITE(*,*)ZSTE,ZBOSS
        CALL NURBS_INT(NST0S-NST0B,TPVER(NST0B:NST0S,:),NIN-1,UI(1:NIN),VERZ(I,NLOW:NLOW+NIN-1,:))
        !WRITE(*,*)VERZ(I,NLOW:NLOW+NIN,3)
    END DO
    !STOP
    
    !PAUSE
!   FAST SECTION
    !LOWER
    DO J=1,NBLOW
        UI(J)=1./(NBLOW-1)*(J-1)
    END DO
    !WRITE(*,*)NSTB
    TPVER(1:NSTB,:)=VERH(1,1:NSTB,:)
    TPVER(NSTB+1,3)=ZBUL;TPVER(NSTB+1,1)=VERH(1,NSTB,1);TPVER(NSTB+1,2)=0.
    DO I=1,1
        !CALL NURBS_INT(NSTB-1,VERH(I,1:NSTB,:),NBLOW-1,UI(1:NBLOW),VERZ(I,1:NBLOW,:))
        CALL NURBS_INT(NSTB,TPVER(1:NSTB+1,:),NBLOW-1,UI(1:NBLOW),VERZ(I,1:NBLOW,:))
    END DO
    !UPPER
    DO J=1,NBUP
        UI(J)=1./(NBUP-1)*(J-1)
    END DO
    VERH(1,NSTB,3)=ZBUL
    DO I=1,1
        CALL NURBS_INT(NPOI(I)-NSTB,VERH(I,NSTB:NPOI(I),:),NZ-NBLOW,UI(1:NBUP),VERZ(I,NBLOW:NZ,:))
        !WRITE(*,*)VERZ(I,NBLOW:NZ,3)
    END DO

!
!   X DIRECTION
    DO J=1,NX
        UI(J)=1./(NX-1)*(J-1)
    END DO

    DO I=1,NZ
        CALL NURBS_INT(NSEH-1,VERZ(1:NSEH,I,:),NX-1,UI(1:NX),VERI(1:NX,I,:))
    END DO

    VERI(1,NBLOW:NZ,2)=0.
    !VERI(NX,NBLOW,2)=0.

    END IF

END SUBROUTINE

SUBROUTINE MESHASSE
    USE MESHDATA
    

    DO I=1,NX+500
        DO J=1,NZ+50
            VERI(I,J,3)=VERI(I,J,3)+SINK
        END DO
    END DO

    DO I=1,NX+500
        DO J=1,NZ+50
            TPX=VERI(I,J,1);TPZ=VERI(I,J,3)
            VERI(I,J,1)=TPX*COS(TRI)+TPZ*SIN(TRI)
            VERI(I,J,3)=-TPX*SIN(TRI)+TPZ*COS(TRI)
        END DO
    END DO


	OPEN(119,FILE='SHIPHULL_MESH_4.DAT')
	OPEN(120,FILE='SHIPHULL_MESH.DAT')
    
	IF(MSTEN.EQ.0.AND.MBULB.EQ.0.AND.MBULS.EQ.0) THEN
        WRITE(119,101)NX*NZ,(NX-1)*(NZ-1)
        WRITE(120,101)NX*NZ,(NX-1)*(NZ-1)/4
	END IF
   	IF(MSTEN.GT.0.AND.MBULB.EQ.0.AND.MBULS.EQ.0) THEN
        WRITE(119,101)NX*NZ+NSX*NUP-NUP,(NX-1)*(NZ-1)+(NSX-1)*(NUP-1)
        WRITE(120,101)NX*NZ+NSX*NUP-NUP,(NX-1)*(NZ-1)/4+(NSX-1)*(NUP-1)/4
    END IF
	IF(MSTEN.GT.0.AND.MBULB.GT.0.AND.MBULS.EQ.0) THEN
        WRITE(119,101)NX*NZ+(NSX-1)*NUP+(NBX-1)*NBLOW,(NX-1)*(NZ-1)+(NSX-1)*(NUP-1)+(NBX-1)*(NBLOW-1)
        WRITE(120,101)NX*NZ+(NSX-1)*NUP+(NBX-1)*NBLOW,(NX-1)*(NZ-1)/4+(NSX-1)*(NUP-1)/4+(NBX-1)*(NBLOW-1)/4
    END IF
	IF(MSTEN.GT.0.AND.MBULB.GT.0.AND.MBULS.GT.0) THEN
        WRITE(119,101)NX*NZ+(NSX-1)*NUP+(NBX-1)*NBLOW+(NBSX-1)*NLOW,(NX-1)*(NZ-1)+(NSX-1)*(NUP-1)+(NBX-1)*(NBLOW-1)+(NBSX-1)*(NLOW-1)
        WRITE(120,101)NX*NZ+(NSX-1)*NUP+(NBX-1)*NBLOW+(NBSX-1)*NLOW,(NX-1)*(NZ-1)/4+(NSX-1)*(NUP-1)/4+(NBX-1)*(NBLOW-1)/4+(NBSX-1)*(NLOW-1)/4
    END IF

101	format(1x,2i8)

    DO I=1,NX
        DO J=1,NZ
            WRITE(119,104)VERI(I,J,1),VERI(I,J,2),VERI(I,J,3),0,0,0
            WRITE(120,104)VERI(I,J,1),VERI(I,J,2),VERI(I,J,3),0,0,0
        END DO
    END DO

    IF(MSTEN.GT.0) THEN
    DO I=NX+2,NX+NSX
        DO J=1,NUP
            WRITE(119,104)VERI(I,J,1),VERI(I,J,2),VERI(I,J,3),0,0,0
            WRITE(120,104)VERI(I,J,1),VERI(I,J,2),VERI(I,J,3),0,0,0
        END DO
    END DO
    END IF

    IF(MBULB.GT.0) THEN
    DO I=NX+NSX+1,NX+NSX+NBX-1
        DO J=1,NBLOW
            WRITE(119,104)VERI(I,J,1),VERI(I,J,2),VERI(I,J,3),0,0,0
            WRITE(120,104)VERI(I,J,1),VERI(I,J,2),VERI(I,J,3),0,0,0
        END DO
    END DO
    END IF

    IF(MBULS.GT.0) THEN
    DO I=NX+NSX+NBX+2,NX+NSX+NBX+NBSX
        DO J=1,NLOW
            WRITE(119,104)VERI(I,J,1),VERI(I,J,2),VERI(I,J,3),0,0,0
            WRITE(120,104)VERI(I,J,1),VERI(I,J,2),VERI(I,J,3),0,0,0
        END DO
    END DO
    END IF

104	FORMAT(1X,6F15.6)
106 FORMAT(1X,9I8)    

    DO I=1,NX-1,2
        DO J=1,NZ-1,2
		    ND1=J+(I-1)*NZ
            ND2=ND1+1
            ND3=ND2+1
            ND4=ND3+NZ
            ND5=ND4+NZ
            ND6=ND5-1
            ND7=ND6-1
            ND8=ND7-NZ
            ND9=ND8+1

            WRITE(119,*)ND1,ND2,ND3,ND4
            WRITE(120,106)ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
        END DO
    END DO

    IF(MSTEN.GT.0) THEN
    DO I=1,NSX-1,2
        DO J=1,NUP-1,2
		    IF(I.EQ.1) THEN
                !ND1=(NX-1)*NZ+J+NLOW-1
                ND1=(NX-1)*NZ+J+(NZ-NUP+1)-1
            ELSE
                ND1=(NX)*(NZ)+J+(I-2)*NUP
            END IF

            ND2=ND1+1
            ND3=ND2+1
            ND4=ND3+NUP
            ND5=ND4+NUP
            ND6=ND5-1
            ND7=ND6-1
            ND8=ND7-NUP
            ND9=ND8+1

            WRITE(119,*)ND1,ND2,ND3,ND4
            WRITE(120,106)ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
        END DO
    END DO
    END IF
    
    !WRITE(*,*)NX*NZ !+(NSX-2)*NUP

    IF(MBULB.GT.0) THEN
    DO I=1,NBX-1,2
        DO J=1,NBLOW-1,2


		    IF(I.EQ.NBX-2) THEN
                ND1=NX*NZ+(NSX-1)*NUP+J+(I-1)*NBLOW
                ND2=ND1+1
                ND3=ND2+1
                ND4=ND3+NBLOW
                ND9=ND4-1
                ND8=ND9-1
                ND7=J
                ND6=ND7+1
                ND5=ND6+1
            ELSE
                ND1=NX*NZ+(NSX-1)*NUP+J+(I-1)*NBLOW
                ND2=ND1+1
                ND3=ND2+1
                ND4=ND3+NBLOW
                ND5=ND4+NBLOW
                ND6=ND5-1
                ND7=ND6-1
                ND8=ND7-NBLOW
                ND9=ND8+1
            END IF

            WRITE(119,*)ND1,ND2,ND3,ND4
            WRITE(120,106)ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
            !write(*,*)ND1
        END DO
    END DO
    END IF


    IF(MBULS.GT.0) THEN
    DO I=1,NBSX-1,2
        DO J=1,NLOW-1,2
            
		    IF(I.EQ.1) THEN
                ND1=(NX-1)*NZ+J
                ND2=ND1+1
                ND3=ND2+1
                ND4=NX*NZ+(NSX-1)*NUP+(NBX-1)*NBLOW+(I-1)*NLOW+J+2
                ND5=ND4+NLOW
                ND6=ND5-1
                ND7=ND6-1
                ND8=ND7-NLOW
                ND9=ND8+1
            ELSE
                ND1=NX*NZ+(NSX-1)*NUP+(NBX-1)*NBLOW+(I-2)*NLOW+J
                ND2=ND1+1
                ND3=ND2+1
                ND4=ND3+NLOW
                ND5=ND4+NLOW
                ND6=ND5-1
                ND7=ND6-1
                ND8=ND7-NLOW
                ND9=ND8+1
            END IF


            WRITE(119,*)ND1,ND2,ND3,ND4
            WRITE(120,106)ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
        END DO
    END DO
    END IF
    CLOSE(119)
    CLOSE(120)
END SUBROUTINE


SUBROUTINE NURBS_INT_SI(N,VER,UI0,C)
    
    DIMENSION VER(0:N,3),D(N),U0(0:N),U(0:N+6),RIG(0:N+2)
    DIMENSION P(0:N+2,3) !,UI(0:M+6)
    INTEGER N
    REAL NIP,UI0,C(3)
    REAL,ALLOCATABLE ::MA(:,:)
    ALLOCATE(MA(N+3,N+3))


    DO I=1,N
        D(I)=SQRT((VER(I,1)-VER(I-1,1))**2+(VER(I,2)-VER(I-1,2))**2+(VER(I,3)-VER(I-1,3))**2)
    END DO

    DTOL=SUM(D(:))

    U0(0)=0.
    U0(N)=1.
    DO K=1,N-1
        U0(K)=U0(K-1)+D(K)/DTOL
    END DO

    U(0:3)=0.
    U(N+3:N+6)=1.
    DO J=1,N-1
        U(J+3)=U0(J)
    END DO
    
    P=0.
    DO J=1,3
        D0=(VER(1,J)-VER(0,J))/(U0(1)-U0(0))*U(4)/3.
        DN=(VER(N,J)-VER(N-1,J))/(U0(N)-U0(N-1))*(1.-U(N+2))/3.
            
        !d0=0.
        !dn=0.

        MA=0.
        MA(1,1)=1.
        RIG(0)=VER(0,J)
        DO K=1,N-1
            CALL BASE_FUNC(N,U,U0(K),K,MA(K+1,K+1))
            CALL BASE_FUNC(N,U,U0(K),K+1,MA(K+1,K+2))
            CALL BASE_FUNC(N,U,U0(K),K+2,MA(K+1,K+3))
            RIG(K)=VER(K,J)
            !WRITE(*,*)MA(K+1,K+1),MA(K+1,K+2),MA(K+1,K+3)
        END DO
        RIG(N)=VER(N,J)
        RIG(N+1)=D0
        RIG(N+2)=DN

        MA(N+1,N+3)=1.
        MA(N+2,1)=-1.
        MA(N+2,2)=1.
        MA(N+3,N+2)=-1.
        MA(N+3,N+3)=1.
    
        !write(*,*)MA(N+1,:)
        !write(*,*)
        CALL BRINV(MA,N+3)
        DO I=1,N+3
            DO K=1,N+3
                P(I-1,J)=P(I-1,J)+MA(I,K)*RIG(K-1)
            END DO
            !WRITE(*,*)P(I-1,J)
        END DO
        !write(*,*)
    END DO
        !STOP


    !STOP
    C=0.
    DO I=0,N+2 
        !WRITE(*,*)P(I,2)
        CALL BASE_FUNC(N,U,UI0,I,NIP)
        C(:)=C(:)+NIP*P(I,:)
        !C(J,:)=NIP
        !if(UI0(J).eq.1) WRITE(*,*)NIP,UI0(J)
    END DO

    IF(UI0.EQ.1) THEN
        C(:)=VER(N,:)
    END IF


    DEALLOCATE(MA)
    RETURN
END SUBROUTINE

SUBROUTINE NURBS_INT(N,VER,M,UI0,C)
    
    DIMENSION VER(0:N,3),D(N),U0(0:N),U(0:N+6),RIG(0:N+2),UI0(0:M)
    DIMENSION P(0:N+2,3),UI(0:M+6),C(0:M,3)
    INTEGER N
    REAL NIP
    REAL,ALLOCATABLE ::MA(:,:)
    ALLOCATE(MA(N+3,N+3))


    DO I=1,N
        D(I)=SQRT((VER(I,1)-VER(I-1,1))**2+(VER(I,2)-VER(I-1,2))**2+(VER(I,3)-VER(I-1,3))**2)
    END DO

    DTOL=SUM(D(:))

    U0(0)=0.
    U0(N)=1.
    DO K=1,N-1
        U0(K)=U0(K-1)+D(K)/DTOL
    END DO

    U(0:3)=0.
    U(N+3:N+6)=1.
    DO J=1,N-1
        U(J+3)=U0(J)
    END DO
    
    P=0.
    DO J=1,3
        D0=(VER(1,J)-VER(0,J))/(U0(1)-U0(0))*U(4)/3.
        DN=(VER(N,J)-VER(N-1,J))/(U0(N)-U0(N-1))*(1.-U(N+2))/3.
            
        !d0=0.
        !dn=0.

        MA=0.
        MA(1,1)=1.
        RIG(0)=VER(0,J)
        DO K=1,N-1
            CALL BASE_FUNC(N,U,U0(K),K,MA(K+1,K+1))
            CALL BASE_FUNC(N,U,U0(K),K+1,MA(K+1,K+2))
            CALL BASE_FUNC(N,U,U0(K),K+2,MA(K+1,K+3))
            RIG(K)=VER(K,J)
            !WRITE(*,*)MA(K+1,K+1),MA(K+1,K+2),MA(K+1,K+3)
        END DO
        RIG(N)=VER(N,J)
        RIG(N+1)=D0
        RIG(N+2)=DN

        MA(N+1,N+3)=1.
        MA(N+2,1)=-1.
        MA(N+2,2)=1.
        MA(N+3,N+2)=-1.
        MA(N+3,N+3)=1.
    
        !write(*,*)MA(N+1,:)
        !write(*,*)
        CALL BRINV(MA,N+3)
        DO I=1,N+3
            DO K=1,N+3
                P(I-1,J)=P(I-1,J)+MA(I,K)*RIG(K-1)
            END DO
            !WRITE(*,*)P(I-1,J)
        END DO
        !write(*,*)
    END DO
        !STOP


    !STOP
    C=0.
    DO J=0,M
        DO I=0,N+2 
        !WRITE(*,*)P(I,2)
            CALL BASE_FUNC(N,U,UI0(J),I,NIP)
            C(J,:)=C(J,:)+NIP*P(I,:)
            !C(J,:)=NIP
            !if(UI0(J).eq.1) WRITE(*,*)NIP,UI0(J)
        END DO
        !WRITE(*,*)
    END DO
    IF(UI0(M).EQ.1) THEN
        C(M,:)=VER(N,:)
    END IF


    DEALLOCATE(MA)
    RETURN
END SUBROUTINE


SUBROUTINE BASE_FUNC(N,U,UT,I,NIB)
    DIMENSION U(0:N+6)
    INTEGER N,P
    !DOUBLE PRECISION NB
    REAL NIB,UT,NBA,NBB
    REAL,ALLOCATABLE:: NB(:,:)
    ALLOCATE(NB(0:N+6,0:3))

    NB=0.

    IF(UT.GE.U(I).AND.UT.LT.U(I+1)) THEN
        NB(I,0)=1.
    ELSE
        NB(I,0)=0.
    END IF

    IF(UT.GE.U(I+1).AND.UT.LT.U(I+2)) THEN
        NB(I+1,0)=1.
    ELSE
        NB(I+1,0)=0.
    END IF

    IF(UT.GE.U(I+2).AND.UT.LT.U(I+3)) THEN
        NB(I+2,0)=1.
    ELSE
        NB(I+2,0)=0.
    END IF

    IF(UT.GE.U(I+3).AND.UT.LT.U(I+4)) THEN
        NB(I+3,0)=1.
    ELSE
        NB(I+3,0)=0.
    END IF

    DO P=1,3
        DO J=I,I+3-P
            !NB(J,P)=(UT-U(J))/(U(J+P)-U(J))*NB(J,P-1)+(U(J+P+1)-UT)/(U(J+P+1)-U(J+1))*NB(J+1,P-1)
            
            NBA=(UT-U(J))/(U(J+P)-U(J))*NB(J,P-1)
            NBB=(U(J+P+1)-UT)/(U(J+P+1)-U(J+1))*NB(J+1,P-1)

            IF(ABS(U(J+P)-U(J)).LE.1.E-9) THEN
                NBA=0.
            END IF

            IF(ABS(U(J+P+1)-U(J+1)).LE.1.E-9) THEN
                NBB=0.
            END IF

            !IF(U(J+P+1)-UT.EQ.0.OR.NB(J+1,P-1).EQ.0) THEN
            !    NB(J,P)=(UT-U(J))/(U(J+P)-U(J))*NB(J,P-1)
            !    IF(UT-U(J).EQ.0.OR.NB(J,P-1).EQ.0) THEN
            !        NB(J,P)=0.
            !    END IF
            !ELSE IF(UT-U(J).EQ.0.OR.NB(J,P-1).EQ.0) THEN
            !    NB(J,P)=(U(J+P)-UT)/(U(J+P+1)-U(J+1))*NB(J+1,P-1)
            !END IF
            !NB(J,P)=(UT-U(J))/(U(J+P)-U(J))*NB(J,P-1)+(U(J+P+1)-UT)/(U(J+P+1)-U(J+1))*NB(J+1,P-1)
            NB(J,P)=NBA+NBB
            !WRITE(*,*)J,P,NB(J,P)
            !WRITE(*,*)J+P+1
        END DO
    END DO
    !STOP

    NIB=NB(I,3)
    DEALLOCATE(NB)

    RETURN
END SUBROUTINE




