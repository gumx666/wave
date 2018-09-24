SUBROUTINE OFFWIGLEY_3_9NODES
    !船体和自由面的水线进行统一编号
    !需要对水线上的点提供物面的法向矢量，如果是对自由面
    !上的面元,在计算时另外提供法向矢量


    !GENERATE THE OFFSET DATA OF THE VERTEXES AND THE ELEMENTS
    !NTPN(NUMBER OF POINTS),NTP(NUMBER OF PANELS)
	USE GREENMOD
	USE CUMOD

    REAL:: FAC,FACF,FACB,RATIO,LB,LF,LA
    REAL ::TMX,TMY,TMZ,TVX,TVY,TVZ

	DIMENSION X(500,500),Y(500,500),Z(500,500)                       
	DIMENSION VX(500,500),VY(500,500),VZ(500,500)
	DIMENSION AX(500),AY(500,500)
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
	OPEN(17,FILE='PANELVECTOR.DAT')
	OPEN(18,FILE='FREEVECTOR.DAT')
	OPEN(19,FILE='ALLVECTOR.DAT')
	WRITE(*,*) 
	WRITE(*,*)'BEGIN SHIP MESH GENERATE 9NODES'

    !DFF=0.002699+RD
    !SUA=-0.001896 

	IF(ABS(RL/(NX-1)-DFX).GT.1E-6)THEN
		!WRITE(*,*)'ERROR	NX?	
        !STOP
	ENDIF

    NFX=NX

    FAC=2.0
    FACF=1.2
    FACB=1.5
    RATIO=0.4  !linear ok
    
    LS=RL
    LF=1.0*DARL !1.4*DARL
    LA=0.24*DARL !0.6*DARL
    LB=DARB !RATIO*DARL

    DPX=LA
    DPY=LB

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
   
    WRITE(*,*)"NNODES  =    ",NFPOINT+NBPOINT
    WRITE(*,*)"NPANELS =    ",NFBLOCK+NBBLOCK
	
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

        RD=DFF/COS(SUA)+TMX*TAN(SUA)

        !FX(NFX-I+1,J)=((DARL-RL)*(EXP(REAL(I-1)/(NFX-1)*FAC)-1.)/(EXP(FAC)-1.)+RL)*COS((-PI/(NFY-1))*(J-1))

		DO J=1,NZ
			TMZ=-TMX*TAN(SUA)-RD/(NZ-1)*(J-1)
            !TMZ=-TMX*TAN(SUA)-RD*(EXP(REAL(J-1)/((NX-1))*FAC)-1.)/(EXP(FAC)-1.)
			TMY=RB/2.0*(1-(2*TMX/RL)**2)*(1-(TMZ/RD)**2)*(1+0.2*(2*TMX/RL)**2)
			
            X(I,J)=TMX*COS(SUA)+TMZ*SIN(SUA)
            Y(I,J)=TMY
            Z(I,J)=TMZ*COS(SUA)+TMX*SIN(SUA) !-RD/30. !-0.1*RD
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

        END DO
    END DO

    !X(2,4)=X(2,6);Y(2,4)=Y(2,6);Z(2,4)=Z(2,6)
    !X(4,3)=X(4,4);Y(4,3)=Y(4,4);Z(4,3)=Z(4,4)

	K1=0
	DO I=1,NX
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

    DFX=ABS(X(1,1)-X(2,1))*1.0

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

    WRITE(*,*)"DFY",LB,DFY

!FIRST PART OF THE FREE SURFACE (UPSTREAM, BEFORE THE BOW)
	DO I=1,NFXF,2
		DO J=1,NFY,2  
            
            !FX(I,J)=OUTX(NFX,1)+LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            
            FX(I,J)=X(1,1)-LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            !FX(I,J)=X(1,1)-LF*REAL(NFXF-I+1)/(NFXF)
            !FX(I,J)=X(1,1)-DFX*(NFXF-I+1)

            !FY(I,J)=LB*REAL(J-1)/(NFY-1)
            	
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
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    DO J=1,NFY
            FX(1,J)=(FX(NFXF-1,J)+X(1,1))/2.
            FY(1,J)=FY(NFXF-1,J)
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
            FY(I,J)=Y(I,1)+(LB-Y(I,2))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            VX2=0.0
		    VY2=0.0
			VZ2=1.0

            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
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
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
            !END IF
        END DO
    END DO

!THIRD PART OF THE FREE SURFACE (UPSTREAM, AFTER THE STERN)
    !DFX=LA/(NFXA)
    !DFY=LB/(NFY-1)
    DFY=ABS(FY(1,1)-FY(1,2))

	DO I=2,NFXA,2
		DO J=1,NFY,2  
            !FX(I,J)=TX(1)-(I-NFX)*DFX
			!FY(I,J)=(J-1)*DFY

            !FX(I,J)=OUTX(1,1)-LA*(EXP(REAL(I)/(NFXA)*FAC)-1.)/(EXP(FAC)-1.)
            
            FX(I,J)=X(NFX,1)+LA*(EXP(REAL(I)/(NFXA)*FACF)-1.)/(EXP(FACF)-1.)
            
            FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*REAL(I)/(NFXA)
            !FX(I,J)=X(NFX,1)+DFX*REAL(I)
            
            !FY(I,J)=LB*REAL(J-1)/(NFY-1)
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
END SUBROUTINE


SUBROUTINE OFFWIGLEY_9NODES_N
    !船体和自由面的水线进行统一编号
    !需要对水线上的点提供物面的法向矢量，如果是对自由面
    !上的面元,在计算时另外提供法向矢量


    !GENERATE THE OFFSET DATA OF THE VERTEXES AND THE ELEMENTS
    !NTPN(NUMBER OF POINTS),NTP(NUMBER OF PANELS)
	USE GREENMOD
	USE CUMOD

    REAL:: FAC,FACF,FACB,RATIO,LB,LF,LA
    REAL ::TMX,TMY,TMZ,TVX,TVY,TVZ
    REAL:: WX(100),WY(100),WXX(100),WYY(100)

	DIMENSION X(500,500),Y(500,500),Z(500,500)                       
	DIMENSION VX(500,500),VY(500,500),VZ(500,500)
	DIMENSION AX(500),AY(500,500)
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
	OPEN(17,FILE='PANELVECTOR.DAT')
	OPEN(18,FILE='FREEVECTOR.DAT')
	OPEN(19,FILE='ALLVECTOR.DAT')
	WRITE(*,*) 
	WRITE(*,*)'BEGIN SHIP MESH GENERATE 9NODES'

    !DFF=0.002699+RD
    !SUA=-0.001896 

	IF(ABS(RL/(NX-1)-DFX).GT.1E-6)THEN
		!WRITE(*,*)'ERROR	NX?	
        !STOP
	ENDIF

    FAC=2.0
    FACF=1.8
    FACB=1.5
    RATIO=0.4  !linear ok
    
    LS=RL
    LF=1.0*DARL !1.4*DARL
    LA=0.24*DARL !0.6*DARL
    LB=DARB !RATIO*DARL

    DPX=LA
    DPY=LB

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
   
    WRITE(*,*)"NNODES  =    ",NFPOINT+NBPOINT
    WRITE(*,*)"NPANELS =    ",NFBLOCK+NBBLOCK
	
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

        RD=DFF/COS(SUA)+TMX*TAN(SUA)

        !FX(NFX-I+1,J)=((DARL-RL)*(EXP(REAL(I-1)/(NFX-1)*FAC)-1.)/(EXP(FAC)-1.)+RL)*COS((-PI/(NFY-1))*(J-1))

		DO J=1,NZ
			TMZ=-TMX*TAN(SUA)-RD/(NZ-1)*(J-1)
            !TMZ=-TMX*TAN(SUA)-RD*(EXP(REAL(J-1)/((NX-1))*FAC)-1.)/(EXP(FAC)-1.)
			TMY=RB/2.0*(1-(2*TMX/RL)**2)*(1-(TMZ/RD)**2) !*(1+0.2*(2*TMX/RL)**2)
			
            X(I,J)=TMX*COS(SUA)+TMZ*SIN(SUA)
            Y(I,J)=TMY
            Z(I,J)=TMZ*COS(SUA)+TMX*SIN(SUA) !-RD/30. !-0.1*RD
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

        END DO
    END DO

    !X(2,4)=X(2,6);Y(2,4)=Y(2,6);Z(2,4)=Z(2,6)
    !X(4,3)=X(4,4);Y(4,3)=Y(4,4);Z(4,3)=Z(4,4)

	K1=0
	DO I=1,NX
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

    !DFX=ABS(X(1,1)-X(2,1))*1.0
    DFX=RL/(NFX-1)

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

    WRITE(*,*)"DFY",LB,DFY

    DO I=1,NX
        WX(I)=X(I,1)
        WY(I)=Y(I,1)
    END DO

    DO I=1,NFX
        WXX(I)=WX(1)+(I-1)*DFX
    END DO
    !调用样条函数
    CALL ESPL2(WX(1:NX),WY(1:NX),NX,WXX(1:NFX),NFX,WYY(1:NFX))
    !WX(1:NX)=WXX(1:NX)


!FIRST PART OF THE FREE SURFACE (UPSTREAM, BEFORE THE BOW)
	DO I=1,NFXF,2
		DO J=1,NFY,2  
            
            !FX(I,J)=OUTX(NFX,1)+LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            
            !FX(I,J)=X(1,1)-LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            !FX(I,J)=X(1,1)-LF*REAL(NFXF-I+1)/(NFXF)
            FX(I,J)=X(1,1)-DFX*(NFXF-I+1)

            !FY(I,J)=LB*REAL(J-1)/(NFY-1)
            	
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
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    DO J=1,NFY
            FX(1,J)=(FX(NFXF-1,J)+X(1,1))/2.
            FY(1,J)=FY(NFXF-1,J)
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
            
            !FY(I,J)=WYY(I)+(LB-WYY(I))*REAL(J-1)/(NFY-1)

            !FX(I,J)=X(I,1)
            FX(I,J)=X(1,1)+DFX*(I-1)
            FY(I,J)=WYY(I)+(LB-WYY(I))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            VX2=0.0
		    VY2=0.0
			VZ2=1.0

            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
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
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
            !END IF
        END DO
    END DO

!THIRD PART OF THE FREE SURFACE (UPSTREAM, AFTER THE STERN)
    !DFX=LA/(NFXA)
    !DFY=LB/(NFY-1)
    DFY=ABS(FY(1,1)-FY(1,2))

	DO I=2,NFXA,2
		DO J=1,NFY,2  
            !FX(I,J)=TX(1)-(I-NFX)*DFX
			!FY(I,J)=(J-1)*DFY

            !FX(I,J)=OUTX(1,1)-LA*(EXP(REAL(I)/(NFXA)*FAC)-1.)/(EXP(FAC)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*(EXP(REAL(I)/(NFXA)*FACF)-1.)/(EXP(FACF)-1.)
            
            FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*REAL(I)/(NFXA)
            FX(I,J)=X(NX,1)+DFX*I
            
            FY(I,J)=LB*REAL(J-1)/(NFY-1)
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
END SUBROUTINE


SUBROUTINE OFFSHIP
    !船体和自由面的水线进行统一编号
    !需要对水线上的点提供物面的法向矢量，如果是对自由面
    !上的面元,在计算时另外提供法向矢量


    !GENERATE THE OFFSET DATA OF THE VERTEXES AND THE ELEMENTS
    !NTPN(NUMBER OF POINTS),NTP(NUMBER OF PANELS)
	USE GREENMOD
	USE CUMOD

    REAL:: FAC,FACF,FACB,RATIO,LB,LF,LA,WX,WY,WZ,WXX
    REAL ::TMX,TMY,TMZ,TVX,TVY,TVZ,THETA
    INTEGER :: INEB,COUT

	DIMENSION X(5000,5000),Y(5000,5000),Z(5000,5000) 
    DIMENSION XS(25000),YS(25000),ZS(25000),INEB(25000,9)    
    DIMENSION WX(-200:500),WY(-200:500),WZ(5000),WXX(500)                     
	DIMENSION VX(500,500),VY(500,500),VZ(500,500)
	DIMENSION AX(500),AY(500,500)
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
	OPEN(17,FILE='PANELVECTOR.DAT')
	OPEN(18,FILE='FREEVECTOR.DAT')
	OPEN(19,FILE='ALLVECTOR.DAT')
	WRITE(*,*) 
	WRITE(*,*)'BEGIN SHIP MESH GENERATE 9NODES'

	IF(ABS(RL/(NX-1)-DFX).GT.1E-6)THEN
		!WRITE(*,*)'ERROR	NX?	
        !STOP
	ENDIF

    !NFX=NX
    THETA=28.

    FAC=2.0
    FACF=1.8
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
    OPEN(10,FILE='KCS.NEU')
    !OPEN(10,FILE='KCS1.NEU')
    !OPEN(10,FILE='RS9N.neu')
    DO I=1,6
        READ(10,*)
    END DO
    !读入船体面元参数
    READ(10,*)NBPOINT,NBBLOCK
    
    !自由面总的节点数
	!NFPOINT=NFX*NFY
    !自由面总的面元数
	!NFBLOCK=(NFX-1)*(NFY-1)/4
    !控制面总的节点数
	!NRPOINT=0. !NZ*(2*NFY+NFX)
    !控制面总的面元数
	!NRBLOCK=0. !(NZ-1)*(2*NFY-2+NFX-1)/4
	
    READ(10,*)
    READ(10,*)
    !读入坐标数据
    DO I=1,NBPOINT
        READ(10,*)ZS(I),XS(I),YS(I),ZS(I)
        !WRITE(11,104)XS(I),YS(I),ZS(I),0,0,0
        !WRITE(13,104)XS(I),YS(I),ZS(I),0,0,0
	    !WRITE(17,104)XS(I),YS(I),ZS(I),0,0,0
        !WRITE(18,104)XS(I),YS(I),ZS(I),0,0,0
	    !WRITE(19,104)XS(I),YS(I),ZS(I),0,0,0
    END DO

    WY=0.

    !读入节点排列数据
    READ(10,*)
    READ(10,*)
    DO I=1,NBBLOCK
        READ(10,*)A,A,A,INEB(I,1:7)
        READ(10,*)INEB(I,8:9)
    END DO
    CLOSE(10)

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
    DO I=COUT/2,1,-1
        IF(ABS((WY(I+1)-WY(I))/(WX(I+1)-WX(I))).GE.TAN(THETA/180.*PI)) THEN
            WY(I)=WY(I+1)-(WY(I+2)-WY(I+1))
            WX(I)=WX(I+1)-(WX(I+2)-WX(I+1))
        END IF
        WRITE(20,121)I,WX(I),WY(I),WZ(I)
    END DO

    OPEN(20,FILE='WATERLINE.DAT')
    DO I=1,COUT
        WRITE(20,121)I,WX(I),WY(I),WZ(I)
    END DO
    CLOSE(20)
    121FORMAT(I8,4F15.6)

    199 CONTINUE

!   结束

    !DFX=ABS(X(1,1)-X(2,1))
    !DFX=ABS(WX(COUT/2)-WX(COUT/2-1))*1.0 !*1.4 !*1.2
    DFX=ABS(WX(COUT/2)-WX(COUT/2-1))*1.1 !*1.4 !*1.2    !kcs
    
    !DFX=ABS(WX(COUT/2)-WX(COUT/2-1))*1.2 !*1.4 !*1.2    !YUZHENG

    !设置水线节点
    NFX=COUT

    !DFX=(WX(COUT)-WX(1))/(NFX-1)
    !DO I=1,COUT
    !    WXX(I)=WX(1)+(I-1)*DFX
    !END DO
    !调用样条函数
    !CALL ESPL2(WX(1:COUT),WY(1:COUT),COUT,WXX(1:NFX),NFX,WY(1:NFX))
    !WX(1:NFX)=WXX(1:NFX)


    NFPOINT=(NFX+NFXF+NFXA)*NFY
    NFBLOCK=(NFY-1)*(NFX+NFXF+NFXA-1)/4
    
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
   
    WRITE(*,*)"NNODES  =    ",NFPOINT+NBPOINT
    WRITE(*,*)"NPANELS =    ",NFBLOCK+NBBLOCK

	WRITE(11,*)NBPOINT,NBBLOCK*4
	WRITE(12,*)NFPOINT,NFBLOCK
	WRITE(13,*)NAPOINT,NABLOCK

	WRITE(14,102)NBPOINT,NBBLOCK*4
	WRITE(15,102)NFPOINT,NFBLOCK*4
	WRITE(16,102)NRPOINT,NRBLOCK*4

	WRITE(17,101)NBPOINT,NBBLOCK*4
	WRITE(18,101)NFPOINT+NBPOINT,(NFBLOCK+NBBLOCK)*4
	WRITE(19,101)NAPOINT,NABLOCK*4

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

    WRITE(*,*)"DFY",LB,DFY

	DO I=1,NFXF
        WX(I-NFXF)=WX(1)-DFX*(NFXF-I+1)
            !WRITE(*,*)I
    END DO


    !光顺船尾后的水线
    DO I=1,-NFXF,-1
        IF(ABS((WY(I+1)-WY(I))/(WX(I+1)-WX(I))).GE.TAN(THETA/180.*PI)) THEN
            WY(I)=WY(I+1)-(WY(I+2)-WY(I+1))
            WRITE(*,*)I
        END IF
        !WRITE(20,121)I,WX(I),WY(I),WZ(I)
    END DO

!FIRST PART OF THE FREE SURFACE (UPSTREAM, BEFORE THE BOW)
	DO I=1,NFXF
		DO J=1,NFY  
            
            !FX(I,J)=OUTX(NFX,1)+LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            
            !FX(I,J)=WX(1)-LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            !FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            !FX(I,J)=X(1,1)-LF*REAL(NFXF-I+1)/(NFXF)
             
            FX(I,J)=WX(1)-DFX*(NFXF-I+1)
            !FY(I,J)=LB*REAL(J-1)/(NFY-1)
            FY(I,J)=WY(I-NFXF)+(LB-WY(I-NFXF))*REAL(J-1)/(NFY-1)
            	
            VX2=0.0
		    VY2=0.0
			VZ2=1.0
      
            K1=K1+1

			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    GOTO 299
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
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    DO J=1,NFY
            FX(1,J)=(FX(NFXF-1,J)+WX(1))/2.
            !FX(1,J)=(FX(NFXF-1,J)+X(1,1))/2.
            FY(1,J)=FY(NFXF-1,J)
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(1,J),FY(1,J),FZ,VX2,VY2,VZ2
    END DO
    299 CONTINUE

!SECOND PART OF THE FREE SURFACE (ALONG SHIP)
	DO I=1,NFX
		DO J=1,NFY,2  
            !DFY=(LB-TY(NFX+NFXF-I+1))/(NFY-1)

            !FX(I,J)=OUTX(NFX-I+1,1)
            !FX(I,J)=X(I,1)
            
            !FY(I,J)=WY(I)+(LB-WY(I))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            FY(I,J)=WY(I)+(LB-WY(I))*REAL(J-1)/(NFY-1)
            
            !FY(I,J)=Y(I,1)+(LB-Y(I,2))*REAL(J-1)/(NFY-1)
            !FX(I,J)=X(I,1)
            FX(I,J)=WX(I)

            !FY(I,J)=Y(I,1)+(LB-Y(I,2))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            VX2=0.0
		    VY2=0.0
			VZ2=1.0

            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
                WRITE(13,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
			    WRITE(18,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
			    WRITE(19,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
            END IF

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
    DFY=LB/40.
    WRITE(*,*)"DFY     ",DFY

!THIRD PART OF THE FREE SURFACE (UPSTREAM, AFTER THE STERN)
    !DFX=LA/(NFXA)
    !DFY=LB/(NFY-1)

	DO I=2,NFXA,2
		DO J=1,NFY,2  
            !FX(I,J)=TX(1)-(I-NFX)*DFX
			!FY(I,J)=(J-1)*DFY

            !FX(I,J)=OUTX(1,1)-LA*(EXP(REAL(I)/(NFXA)*FAC)-1.)/(EXP(FAC)-1.)
            
            !FX(I,J)=WX(NFX)+LA*(EXP(REAL(I)/(NFXA)*FACF)-1.)/(EXP(FACF)-1.)
            
            !FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*REAL(I)/(NFXA)

            FX(I,J)=WX(NFX)+DFX*REAL(I)
            FY(I,J)=LB*REAL(J-1)/(NFY-1)
            
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

			!WRITE(13,105)9,ND1,ND2,ND3,ND4,ND5,ND6,ND7,ND8,ND9
            WRITE(13,105)9,ND1,ND8,ND7,ND6,ND5,ND4,ND3,ND2,ND9

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

SUBROUTINE OFFWIGLEY_9NODES
    !船体和自由面的水线进行统一编号
    !需要对水线上的点提供物面的法向矢量，如果是对自由面
    !上的面元,在计算时另外提供法向矢量


    !GENERATE THE OFFSET DATA OF THE VERTEXES AND THE ELEMENTS
    !NTPN(NUMBER OF POINTS),NTP(NUMBER OF PANELS)
	USE GREENMOD
	USE CUMOD

    REAL:: FAC,FACF,FACB,RATIO,LB,LF,LA
    REAL ::TMX,TMY,TMZ,TVX,TVY,TVZ

	DIMENSION X(500,500),Y(500,500),Z(500,500)                       
	DIMENSION VX(500,500),VY(500,500),VZ(500,500)
	DIMENSION AX(500),AY(500,500)
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
	OPEN(17,FILE='PANELVECTOR.DAT')
	OPEN(18,FILE='FREEVECTOR.DAT')
	OPEN(19,FILE='ALLVECTOR.DAT')
	WRITE(*,*) 
	WRITE(*,*)'BEGIN SHIP MESH GENERATE 9NODES'

    !DFF=0.002699+RD
    !SUA=-0.001896 

	IF(ABS(RL/(NX-1)-DFX).GT.1E-6)THEN
		!WRITE(*,*)'ERROR	NX?	
        !STOP
	ENDIF

    NFX=NX

    FAC=2.0
    FACF=1.1
    FACB=1.2
    RATIO=0.4  !linear ok
    
    LS=RL
    LF=1.0*DARL !1.4*DARL
    LA=0.24*DARL !0.6*DARL
    LB=DARB !RATIO*DARL

    DPX=LA
    DPY=LB

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
   
    WRITE(*,*)"NNODES  =    ",NFPOINT+NBPOINT
    WRITE(*,*)"NPANELS =    ",NFBLOCK+NBBLOCK
	
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

        RD=DFF/COS(SUA)+TMX*TAN(SUA)

        !FX(NFX-I+1,J)=((DARL-RL)*(EXP(REAL(I-1)/(NFX-1)*FAC)-1.)/(EXP(FAC)-1.)+RL)*COS((-PI/(NFY-1))*(J-1))

		DO J=1,NZ
			TMZ=-TMX*TAN(SUA)-RD/(NZ-1)*(J-1)
            !TMZ=-TMX*TAN(SUA)-RD*(EXP(REAL(J-1)/((NX-1))*FAC)-1.)/(EXP(FAC)-1.)
			TMY=RB/2.0*(1-(2*TMX/RL)**2)*(1-(TMZ/RD)**2) !*(1+0.2*(2*TMX/RL)**2)
			
            X(I,J)=TMX*COS(SUA)+TMZ*SIN(SUA)
            Y(I,J)=TMY
            Z(I,J)=TMZ*COS(SUA)+TMX*SIN(SUA) !-RD/30. !-0.1*RD
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

        END DO
    END DO

    !X(2,4)=X(2,6);Y(2,4)=Y(2,6);Z(2,4)=Z(2,6)
    !X(4,3)=X(4,4);Y(4,3)=Y(4,4);Z(4,3)=Z(4,4)

	K1=0
	DO I=1,NX
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

    DFX=ABS(X(1,1)-X(2,1))*1.0

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

    WRITE(*,*)"DFY",LB,DFY

!FIRST PART OF THE FREE SURFACE (UPSTREAM, BEFORE THE BOW)
	DO I=1,NFXF,2
		DO J=1,NFY,2  
            
            !FX(I,J)=OUTX(NFX,1)+LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            
            !FX(I,J)=X(1,1)-LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            FX(I,J)=X(1,1)-LF*REAL(NFXF-I+1)/(NFXF)
            !FX(I,J)=X(1,1)-DFX*(NFXF-I+1)

            !FY(I,J)=LB*REAL(J-1)/(NFY-1)
            	
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
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    DO J=1,NFY
            FX(1,J)=(FX(NFXF-1,J)+X(1,1))/2.
            FY(1,J)=FY(NFXF-1,J)
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
            FY(I,J)=Y(I,1)+(LB-Y(I,2))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            VX2=0.0
		    VY2=0.0
			VZ2=1.0

            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
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
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
            !END IF
        END DO
    END DO

!THIRD PART OF THE FREE SURFACE (UPSTREAM, AFTER THE STERN)
    !DFX=LA/(NFXA)
    !DFY=LB/(NFY-1)
    DFY=ABS(FY(1,1)-FY(1,2))
    DFY=LB/40.

	DO I=2,NFXA,2
		DO J=1,NFY,2  
            !FX(I,J)=TX(1)-(I-NFX)*DFX
			!FY(I,J)=(J-1)*DFY

            !FX(I,J)=OUTX(1,1)-LA*(EXP(REAL(I)/(NFXA)*FAC)-1.)/(EXP(FAC)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*(EXP(REAL(I)/(NFXA)*FACF)-1.)/(EXP(FACF)-1.)
            
            FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            FX(I,J)=X(NFX,1)+LA*REAL(I)/(NFXA)
            !FX(I,J)=X(NFX,1)+DFX*REAL(I)
            
            !FY(I,J)=LB*REAL(J-1)/(NFY-1)
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
END SUBROUTINE

SUBROUTINE OFFWIGLEY_9NODES_UN
    !船体和自由面的水线进行统一编号
    !需要对水线上的点提供物面的法向矢量，如果是对自由面
    !上的面元,在计算时另外提供法向矢量


    !GENERATE THE OFFSET DATA OF THE VERTEXES AND THE ELEMENTS
    !NTPN(NUMBER OF POINTS),NTP(NUMBER OF PANELS)
	USE GREENMOD
	USE CUMOD

    REAL:: FAC,FACF,FACB,RATIO,LB,LF,LA
    REAL ::TMX,TMY,TMZ,TVX,TVY,TVZ

	DIMENSION X(500,500),Y(500,500),Z(500,500)                       
	DIMENSION VX(500,500),VY(500,500),VZ(500,500)
	DIMENSION AX(500),AY(500,500)
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
	OPEN(17,FILE='PANELVECTOR.DAT')
	OPEN(18,FILE='FREEVECTOR.DAT')
	OPEN(19,FILE='ALLVECTOR.DAT')
	WRITE(*,*) 
	WRITE(*,*)'BEGIN SHIP MESH GENERATE 9NODES'


	IF(ABS(RL/(NX-1)-DFX).GT.1E-6)THEN
		!WRITE(*,*)'ERROR	NX?	
        !STOP
	ENDIF

    NFX=NX

    FAC=2.0
    FACF=1.3
    FACB=1.2
    RATIO=0.4  !linear ok
    
    LS=RL
    LF=1.0*DARL !1.4*DARL
    LA=0.24*DARL !0.6*DARL
    LB=DARB !RATIO*DARL

    DPX=LA
    DPY=LB

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
   
    WRITE(*,*)"NNODES  =    ",NFPOINT+NBPOINT
    WRITE(*,*)"NPANELS =    ",NFBLOCK+NBBLOCK
	
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


        RD=DFF/COS(SUA)+TMX*TAN(SUA)

        !FX(NFX-I+1,J)=((DARL-RL)*(EXP(REAL(I-1)/(NFX-1)*FAC)-1.)/(EXP(FAC)-1.)+RL)*COS((-PI/(NFY-1))*(J-1))

		DO J=1,NZ
			TMZ=-TMX*TAN(SUA)-RD/(NZ-1)*(J-1)
            !TMZ=-TMX*TAN(SUA)-RD*(EXP(REAL(J-1)/((NX-1))*FAC)-1.)/(EXP(FAC)-1.)
			TMY=RB/2.0*(1-(2*TMX/RL)**2)*(1-(TMZ/RD)**2) !*(1+0.2*(2*TMX/RL)**2)
			
            X(I,J)=TMX*COS(SUA)+TMZ*SIN(SUA)
            Y(I,J)=TMY
            Z(I,J)=TMZ*COS(SUA)+TMX*SIN(SUA) !-RD/30. !-0.1*RD
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

        END DO
    END DO

    !X(2,4)=X(2,6);Y(2,4)=Y(2,6);Z(2,4)=Z(2,6)
    !X(4,3)=X(4,4);Y(4,3)=Y(4,4);Z(4,3)=Z(4,4)

	K1=0
	DO I=1,NX
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

    DFX=ABS(X(1,1)-X(2,1)) !*1.2
    !DFX=1./20.

103	FORMAT(1X,I8,6F15.6)
104	FORMAT(1X,6F15.6)

    !GOTO 91

    FZ=0
    !DFX=LF/(NFXF)
    DFY=LB/(NFY-1)  !  old
    DFY=LB/40.  !  new

	!DFY=DFX
	K2=0	
	VX2=0.0
	VY2=0.0
	VZ2=1.0

    WRITE(*,*)"DFY",LB,DFY

!FIRST PART OF THE FREE SURFACE (UPSTREAM, BEFORE THE BOW)
	DO I=1,NFXF,2
		DO J=1,NFY,2  
            
            !FX(I,J)=OUTX(NFX,1)+LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            
            !FX(I,J)=X(1,1)-LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            !FX(I,J)=X(1,1)-LF*REAL(NFXF-I+1)/(NFXF)
            FX(I,J)=X(1,1)-DFX*(NFXF-I+1)

            !FY(I,J)=LB*REAL(J-1)/(NFY-1)
            	
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
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    DO J=1,NFY
            FX(1,J)=(FX(NFXF-1,J)+X(1,1))/2.
            FY(1,J)=FY(NFXF-1,J)
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
            FY(I,J)=Y(I,1)+(LB-Y(I,2))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            VX2=0.0
		    VY2=0.0
			VZ2=1.0

            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
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
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
            !END IF
        END DO
    END DO

!THIRD PART OF THE FREE SURFACE (UPSTREAM, AFTER THE STERN)
    !DFX=LA/(NFXA)
    DFY=LB/(NFY-1)
    !DFY=LB/21

	DO I=2,NFXA,2
		DO J=1,NFY,2  
            !FX(I,J)=TX(1)-(I-NFX)*DFX
			!FY(I,J)=(J-1)*DFY

            !FX(I,J)=OUTX(1,1)-LA*(EXP(REAL(I)/(NFXA)*FAC)-1.)/(EXP(FAC)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*(EXP(REAL(I)/(NFXA)*FACF)-1.)/(EXP(FACF)-1.)
            
            FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*REAL(I)/(NFXA)
            FX(I,J)=X(NFX,1)+DFX*REAL(I)
            
            !FY(I,J)=LB*REAL(J-1)/(NFY-1)
            
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
END SUBROUTINE


SUBROUTINE ESPL2(X,Y,N,XX,M,S)

	DIMENSION X(N),Y(N),XX(M),DY(N),DDY(N)
    DIMENSION S(M),DS(M),DDS(M),H(N)
	!DOUBLE PRECISION X,Y,XX,DY,DDY,S,DS,DDS,H,DY1,DYN,T,H0,H1,BETA,ALPHA
	
    DY1=0.
    DYN=0.

    DY(1)=-0.5
	H0=X(2)-X(1)
	H(1)=3.0*(Y(2)-Y(1))/(2.0*H0)-DY1*H0/4.0
	DO 10 J=2,N-1
	  H1=X(J+1)-X(J)
	  ALPHA=H0/(H0+H1)
	  BETA=(1.0-ALPHA)*(Y(J)-Y(J-1))/H0
	  BETA=3.0*(BETA+ALPHA*(Y(J+1)-Y(J))/H1)
	  DY(J)=-ALPHA/(2.0+(1.0-ALPHA)*DY(J-1))
	  H(J)=(BETA-(1.0-ALPHA)*H(J-1))
	  H(J)=H(J)/(2.0+(1.0-ALPHA)*DY(J-1))
	  H0=H1
10	CONTINUE
	DY(N)=(3.0*(Y(N)-Y(N-1))/H1+DYN*H1/2.0-H(N-1)) /(2.0+DY(N-1))
	DO 20 J=N-1,1,-1
20	DY(J)=DY(J)*DY(J+1)+H(J)
	DO 30 J=1,N-1
30	H(J)=X(J+1)-X(J)
	DO 40 J=1,N-1
	  H1=H(J)*H(J)
	  DDY(J)=6.0*(Y(J+1)-Y(J))/H1-2.0*(2.0*DY(J)+DY(J+1))/H(J)
40	CONTINUE
	H1=H(N-1)*H(N-1)
	DDY(N)=6.0*(Y(N-1)-Y(N))/H1+2.0*(2.0*DY(N)+DY(N-1))/H(N-1)
	T=0.0
	DO 50 I=1,N-1
	  H1=0.5*H(I)*(Y(I)+Y(I+1))
	  H1=H1-H(I)*H(I)*H(I)*(DDY(I)+DDY(I+1))/24.0
	  T=T+H1
50	CONTINUE
	DO 70 J=1,M
	  IF (XX(J).GE.X(N)) THEN
	    I=N-1
	  ELSE
	    I=1
60	    IF (XX(J).GT.X(I+1)) THEN
	      I=I+1
	      GOTO 60
	    END IF
	  END IF
	  H1=(X(I+1)-XX(J))/H(I)
	  S(J)=(3.0*H1*H1-2.0*H1*H1*H1)*Y(I)
	  S(J)=S(J)+H(I)*(H1*H1-H1*H1*H1)*DY(I)
	  DS(J)=6.0*(H1*H1-H1)*Y(I)/H(I)
	  DS(J)=DS(J)+(3.0*H1*H1-2.0*H1)*DY(I)
	  DDS(J)=(6.0-12.0*H1)*Y(I)/(H(I)*H(I))
	  DDS(J)=DDS(J)+(2.0-6.0*H1)*DY(I)/H(I)
	  H1=(XX(J)-X(I))/H(I)
	  S(J)=S(J)+(3.0*H1*H1-2.0*H1*H1*H1)*Y(I+1)
	  S(J)=S(J)-H(I)*(H1*H1-H1*H1*H1)*DY(I+1)
	  DS(J)=DS(J)-6.0*(H1*H1-H1)*Y(I+1)/H(I)
	  DS(J)=DS(J)+(3.0*H1*H1-2.0*H1)*DY(I+1) 
	  DDS(J)=DDS(J)+(6.0-12.0*H1)*Y(I+1)/(H(I)*H(I))
	  DDS(J)=DDS(J)-(2.0-6.0*H1)*DY(I+1)/H(I)
70	CONTINUE
	RETURN
END SUBROUTINE




SUBROUTINE OFFSHIP_NU
    !船体和自由面的水线进行统一编号
    !需要对水线上的点提供物面的法向矢量，如果是对自由面
    !上的面元,在计算时另外提供法向矢量


    !GENERATE THE OFFSET DATA OF THE VERTEXES AND THE ELEMENTS
    !NTPN(NUMBER OF POINTS),NTP(NUMBER OF PANELS)
	USE GREENMOD
	USE CUMOD

    REAL:: FAC,FACF,FACB,RATIO,LB,LF,LA,WX,WY,WZ,WXX
    REAL ::TMX,TMY,TMZ,TVX,TVY,TVZ
    INTEGER :: INEB,COUT

	DIMENSION X(5000,5000),Y(5000,5000),Z(5000,5000) 
    DIMENSION XS(25000),YS(25000),ZS(25000),INEB(25000,9)    
    DIMENSION WX(5000),WY(5000),WZ(5000),WXX(500)                     
	DIMENSION VX(500,500),VY(500,500),VZ(500,500)
	DIMENSION AX(500),AY(500,500)
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
	OPEN(17,FILE='PANELVECTOR.DAT')
	OPEN(18,FILE='FREEVECTOR.DAT')
	OPEN(19,FILE='ALLVECTOR.DAT')
	WRITE(*,*) 
	WRITE(*,*)'BEGIN SHIP MESH GENERATE 9NODES'

	IF(ABS(RL/(NX-1)-DFX).GT.1E-6)THEN
		!WRITE(*,*)'ERROR	NX?	
        !STOP
	ENDIF

    !NFX=NX

    FAC=2.0
    FACF=1.0
    !FACB=2.0
    FACB=1.2
    
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
    OPEN(10,FILE='yuzhengchuan.NEU')
    !OPEN(10,FILE='KCS.NEU')
    !OPEN(10,FILE='RS9N.neu')
    DO I=1,6
        READ(10,*)
    END DO
    !读入船体面元参数
    READ(10,*)NBPOINT,NBBLOCK
    
    !自由面总的节点数
	!NFPOINT=NFX*NFY
    !自由面总的面元数
	!NFBLOCK=(NFX-1)*(NFY-1)/4
    !控制面总的节点数
	!NRPOINT=0. !NZ*(2*NFY+NFX)
    !控制面总的面元数
	!NRBLOCK=0. !(NZ-1)*(2*NFY-2+NFX-1)/4
	
    READ(10,*)
    READ(10,*)
    !读入坐标数据
    DO I=1,NBPOINT
        READ(10,*)ZS(I),XS(I),YS(I),ZS(I)
        !WRITE(11,104)XS(I),YS(I),ZS(I),0,0,0
        !WRITE(13,104)XS(I),YS(I),ZS(I),0,0,0
	    !WRITE(17,104)XS(I),YS(I),ZS(I),0,0,0
        !WRITE(18,104)XS(I),YS(I),ZS(I),0,0,0
	    !WRITE(19,104)XS(I),YS(I),ZS(I),0,0,0
    END DO



    !读入节点排列数据
    READ(10,*)
    READ(10,*)
    DO I=1,NBBLOCK
        READ(10,*)A,A,A,INEB(I,1:7)
        READ(10,*)INEB(I,8:9)
    END DO
    CLOSE(10)

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

    OPEN(20,FILE='WATERLINE.DAT')
    DO I=1,COUT
        WRITE(20,121)I,WX(I),WY(I),WZ(I)
    END DO
    CLOSE(20)
    121FORMAT(I8,4F15.6)

!   光滑水线曲率
    DO I=COUT/2,1,-1
        IF(ABS((WY(I+1)-WY(I))/(WX(I+1)-WX(I))).GE.TAN(30./180.*PI)) THEN
            WY(I)=WY(I+1)-(WY(I+2)-WY(I+1))
            WX(I)=WX(I+1)-(WX(I+2)-WX(I+1))
        END IF

        !WRITE(20,121)I,WX(I),WY(I),WZ(I)
    END DO
    
!   结束

    !DFX=ABS(X(1,1)-X(2,1))
    DFX=ABS(WX(COUT/2)-WX(COUT/2-1))*1.3 !*1.4 !*1.2

    !设置水线节点
    NFX=COUT

    !DFX=(WX(COUT)-WX(1))/(NFX-1)
    !DO I=1,COUT
    !    WXX(I)=WX(1)+(I-1)*DFX
    !END DO
    !调用样条函数
    !CALL ESPL2(WX(1:COUT),WY(1:COUT),COUT,WXX(1:NFX),NFX,WY(1:NFX))
    !WX(1:NFX)=WXX(1:NFX)


    NFPOINT=(NFX+NFXF+NFXA)*NFY
    NFBLOCK=(NFY-1)*(NFX+NFXF+NFXA-1)/4
    
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
   
    WRITE(*,*)"NNODES  =    ",NFPOINT+NBPOINT
    WRITE(*,*)"NPANELS =    ",NFBLOCK+NBBLOCK

	WRITE(11,*)NBPOINT,NBBLOCK*4
	WRITE(12,*)NFPOINT,NFBLOCK
	WRITE(13,*)NAPOINT,NABLOCK

	WRITE(14,102)NBPOINT,NBBLOCK*4
	WRITE(15,102)NFPOINT,NFBLOCK*4
	WRITE(16,102)NRPOINT,NRBLOCK*4

	WRITE(17,101)NBPOINT,NBBLOCK*4
	WRITE(18,101)NFPOINT+NBPOINT,(NFBLOCK+NBBLOCK)*4
	WRITE(19,101)NAPOINT,NABLOCK*4

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

    WRITE(*,*)"DFY",LB,DFY

!FIRST PART OF THE FREE SURFACE (UPSTREAM, BEFORE THE BOW)
	DO I=1,NFXF,2
		DO J=1,NFY,2  
            
            !FX(I,J)=OUTX(NFX,1)+LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            
            FX(I,J)=WX(1)-LF*(EXP(REAL(NFXF-I+1)/(NFXF)*FACF)-1.)/(EXP(FACF)-1.)
            !FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            !FX(I,J)=X(1,1)-LF*REAL(NFXF-I+1)/(NFXF)
             
            !FX(I,J)=WX(1)-DFX*(NFXF-I+1)
            FY(I,J)=LB*REAL(J-1)/(NFY-1)
            	
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
			!WRITE(11,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(13,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(18,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
			WRITE(19,104)FX(I,J),FY(I,J),FZ,VX2,VY2,VZ2
        END DO
    END DO

    DO J=1,NFY
            FX(1,J)=(FX(NFXF-1,J)+WX(1))/2.
            !FX(1,J)=(FX(NFXF-1,J)+X(1,1))/2.
            FY(1,J)=FY(NFXF-1,J)
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
            
            FY(I,J)=WY(I)+(LB-WY(I))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            !FY(I,J)=WY(I)+(LB-WY(I))*REAL(J-1)/(NFY-1)
            
            !FY(I,J)=Y(I,1)+(LB-Y(I,2))*REAL(J-1)/(NFY-1)
            !FX(I,J)=X(I,1)
            FX(I,J)=WX(I)

            !FY(I,J)=Y(I,1)+(LB-Y(I,2))*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)

            VX2=0.0
		    VY2=0.0
			VZ2=1.0

            IF(J.GT.1) THEN
                FX(I,J-1)=(FX(I,J-2)+FX(I,J))/2.
                FY(I,J-1)=(FY(I,J-2)+FY(I,J))/2.
                WRITE(13,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
			    WRITE(18,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
			    WRITE(19,104)FX(I,J-1),FY(I,J-1),FZ,VX2,VY2,VZ2
            END IF

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
    DFY=LB/40.

!THIRD PART OF THE FREE SURFACE (UPSTREAM, AFTER THE STERN)
    !DFX=LA/(NFXA)
    !DFY=LB/(NFY-1)

	DO I=2,NFXA,2
		DO J=1,NFY,2  
            !FX(I,J)=TX(1)-(I-NFX)*DFX
			!FY(I,J)=(J-1)*DFY

            !FX(I,J)=OUTX(1,1)-LA*(EXP(REAL(I)/(NFXA)*FAC)-1.)/(EXP(FAC)-1.)
            
            FX(I,J)=WX(NFX)+LA*(EXP(REAL(I)/(NFXA)*FACF)-1.)/(EXP(FACF)-1.)
            
            !FY(I,J)=LB*(EXP(REAL(J-1)/(NFY-1)*FACB)-1.)/(EXP(FACB)-1.)
            
            !FX(I,J)=X(NFX,1)+LA*REAL(I)/(NFXA)

            !FX(I,J)=WX(NFX)+DFX*REAL(I)
            FY(I,J)=LB*REAL(J-1)/(NFY-1)
            
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