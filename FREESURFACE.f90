SUBROUTINE SURFACEWAVE_HOB
!**************************************************************
!CGENERATE THE FREESURFACE POTENT AND SURFACE WAVES
!**************************************************************
	USE GREENMOD
	USE CUMOD
    USE INPUTDATA
    CHARACTER*7 TITLE,TITLE2
    REAL:: PHIWA(200)

    WRITE(TITLE,'(F7.3)')FR    
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(4,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREESURFACE_LINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    !OPEN(5,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREESURFACE_NONLIN_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    !OPEN(6,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\WAVEPRO_LINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    
    OPEN(4,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREESURFACE_LINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(5,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREESURFACE_NONLIN_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    
    OPEN(6,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\POTENTIAL_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    
    

	!WRITE(*,*) 
	!WRITE(*,*)'SUB SURFACEWAVE......'
	NNF=NFX*NFY
    ALLOCATE(VRIGHT(NBPOINT))
	!ALLOCATE(SP(NTPN))

    M=NBPOINT
    L=NFPOINT

!WRITE(*,*)'FINISHED CALCULATING THE FREE SURFACE POTENTIAL'
!WRITE(*,*)'HAVING STORED IN THE FILE FREEPHI.DAT'

101	format(1x,'title="episode solid mesh"'/1x,'variables="X","Y",& 
     "Z","NX","NY","NZ"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')             
111	format(1x,'title="episode solid mesh"'/1x,'variables="x","Y",& 
     "Z","NX","NY","NZ","SP"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')       
  
     WRITE(6,111)NBPOINT,NBBLOCK*4
     
    !CALL PADCOEF
    
    ETAW=0.
    IF(MTROM.EQ.1) THEN
        WRITE(4,101)NFPOINT,NFBLOCK*4 !+NFXF
    ELSE
        WRITE(4,101)NFPOINT+NFXF+NFX+NFXA,NFBLOCK*4+NFXF+NFX+NFXA-1
    END IF
    
    DO I=NBPOINT+1,M+L
        DO J=NBPOINT+1,NTPN
            !采用差分求导
            ETAW(I)=ETAW(I)+DMX(I,J)*SP(J)*U/G
        END DO
        !根据源分布模型求导

        ETAW(I)=SPX(I)*U/G !-0.5*(SPX(I)**2+SPY(I)**2+SPZ(I)**2)/G
    END DO   

    IF(FR.GE.NSWIN(12,1)) CALL SMOOTHING_RECT   !滤波函数
    
    DO I=1,NBPOINT
        WRITE(6,112)CORD(I,:),VECN(I,:),SP(I)
    END DO
    DO I=1,NBBLOCK
        WRITE(6,*)INE(I,1),INE(I,2),INE(I,9),INE(I,8)
        WRITE(6,*)INE(I,2),INE(I,3),INE(I,4),INE(I,9)
        WRITE(6,*)INE(I,9),INE(I,4),INE(I,5),INE(I,6)
        WRITE(6,*)INE(I,8),INE(I,9),INE(I,6),INE(I,7)
    END DO
    
    
    DO I=NBPOINT+1,M+L
        WRITE(4,102)CORD(I,1:2),ETAW(I),0,0,1
        !WRITE(4,102)CORD(I,1:2),SIGM(I),0,0,1
        !WRITE(4,102)CORD(I,1:2),SP(I),0,0,1	
        !CORD(I,3)=ETAW(I)
    END DO          
    
    
    DO I=4-NFXF,NFX
        IP=NBPOINT+(NFXF+I-1)*NFY+1
        IP1=IP+1
        IP0=IP1+1
        IP_1=IP-NFY+1
        DYL=SQRT((CORD(IP_1,1)-CORD(IP,1))**2+(CORD(IP_1,2)-CORD(IP,2))**2)/(CORD(IP,2)-CORD(IP_1,2))&
            *(PCORD(IP,2)-CORD(IP,2))
        DXL=(CORD(IP,1)-CORD(IP_1,1))/(CORD(IP,2)-CORD(IP_1,2))&
            *(PCORD(IP,2)-CORD(IP,2))

        
        !ETA=ETAW(IP)
        !IF(CORD(IP,1).GE.0.) THEN
        ETA=(PCORD(IP,2)-CORD(IP1,2))*(PCORD(IP,2)-CORD(IP,2))/(CORD(IP0,2)-CORD(IP1,2))/(CORD(IP0,2)-CORD(IP,2))*ETAW(IP0)&
           +(PCORD(IP,2)-CORD(IP0,2))*(PCORD(IP,2)-CORD(IP,2))/(CORD(IP1,2)-CORD(IP0,2))/(CORD(IP1,2)-CORD(IP,2))*ETAW(IP1)&
           +(PCORD(IP,2)-CORD(IP0,2))*(PCORD(IP,2)-CORD(IP1,2))/(CORD(IP,2)-CORD(IP0,2))/(CORD(IP,2)-CORD(IP1,2))*ETAW(IP)
        ETA=ETAW(IP)+(ETAW(IP)-ETAW(IP1))/(CORD(IP,2)-CORD(IP1,2))*(PCORD(IP,2)-CORD(IP,2))
        !END IF

        ETA1=ETAW(IP)+(ETAW(IP)-ETAW(IP_1))/&
             SQRT((CORD(IP_1,1)-CORD(IP,1))**2+(CORD(IP_1,2)-CORD(IP,2))**2)*&
             SQRT(DXL**2+DYL**2)

        IF(I.NE.NFX) THEN
            !ETA=(PHIWA(I+1)-PHIWA(I))*(PCORD(IP+NFY,1)-PCORD(IP,1))/&
            !    ((PCORD(IP+NFY,1)-PCORD(IP,1))**2+(PCORD(IP+NFY,2)-PCORD(IP,2))**2)*U/G
            !ETA=(PHIWA(I+1)-PHIWA(I))/&
            !    SQRT((PCORD(IP+NFY,1)-PCORD(IP,1))**2+(PCORD(IP+NFY,2)-PCORD(IP,2))**2)*U/G
        END IF

        !WRITE(6,*)CORD(IP,1),ETAW(IP),ETA
        !WRITE(6,*)CORD(IP,1)+DXL,ETA1,CORD(IP,1),ETA
    END DO
    
    !CALL SMOOTHING_RECT   !滤波函数
   
    
    DO IP=1,NFX+NFXF+NFXA 
        WRITE(4,102)CORD(NBPOINT+(IP-1)*NFY+1,1),PCORD(NBPOINTP+(IP-1)*NFY+1,2),WPLZ(IP),0,0,1
        !WRITE(5,102)CORD(NBPOINT+(IP-1)*NFY+1,1),PCORD(NBPOINTP+(IP-1)*NFY+1,2),WPLZ(IP),0,0,1
        !WRITE(4,102)CORD(NBPOINT+(IP-1)*NFY+1,1),PCORD(NBPOINTP+(IP-1)*NFY+1,2),WFZ(IP),0,0,1
        !WRITE(4,102)CORD(I,1:2),SP(I-NBPOINT),0,0,1	
        !CORD(I,3)=ETAW(I)
    END DO 
    
    102 FORMAT(6F15.6)	
    112 FORMAT(7F15.6)	
    
    ETAW=0.
    IF(MTROM.EQ.1) THEN
        WRITE(5,101)NFPOINT,NFBLOCK*4 !+NFXF
    ELSE
        WRITE(5,101)NFPOINT+NFXF+NFX+NFXA,NFBLOCK*4+NFXF+NFX+NFXA-1
    END IF
    DO I=NBPOINT+1,M+L
        DO J=NBPOINT+1,NTPN
            !采用差分求导
            ETAW(I)=ETAW(I)+DMX(I,J)*SP(J)*U/G
        END DO
        !根据源分布模型求导

        ETAW(I)=SPX(I)*U/G-0.5*(SPX(I)**2+SPY(I)**2+SPZ(I)**2)/G
        !ETAW(I)=SPX(I)*U/G
    END DO   

    IF(FR.GE.NSWIN(12,1)) CALL SMOOTHING_RECT 
    
    DO J=1,NFY
        !IP=NBPOINT+(NFX+NFXF+NFXA-1)*NFY+J

        !-P1*(0.5*(VV(1)**2+VV(2)**2&    !非线性阻力
     	!		 +VV(3)**2)-U*VV(1))*DN(1)-P1*G*(SR(3))*DN(1)
        !WRITE(*,*)0.5*((-U+SPX(IP))**2+SPY(IP)**2+SPZ(IP)**2)/G+ETAW(IP)-U**2/2./G
    END DO

    DO I=1,NFX+NFXF+NFXA
        IP=NBPOINT+(NFXF+I-1)*NFY+1
        IP1=IP+1
        IP0=IP1+1
        IP=NBPOINT+(I-1)*NFY+1
        
        ETA=ETAW(IP)
        !IF(CORD(IP,1).GE.0.) THEN
            ETA=(PCORD(IP,2)-CORD(IP1,2))*(PCORD(IP,2)-CORD(IP,2))/(CORD(IP0,2)-CORD(IP1,2))/(CORD(IP0,2)-CORD(IP,2))*ETAW(IP0)&
               +(PCORD(IP,2)-CORD(IP0,2))*(PCORD(IP,2)-CORD(IP,2))/(CORD(IP1,2)-CORD(IP0,2))/(CORD(IP1,2)-CORD(IP,2))*ETAW(IP1)&
               +(PCORD(IP,2)-CORD(IP0,2))*(PCORD(IP,2)-CORD(IP1,2))/(CORD(IP,2)-CORD(IP0,2))/(CORD(IP,2)-CORD(IP1,2))*ETAW(IP)
        !END IF
        ETA=ETAW(IP)+(ETAW(IP)-ETAW(IP1))/(CORD(IP,2)-CORD(IP1,2))*(PCORD(IP,2)-CORD(IP,2))
        
        WPLAZ(I)=ETAW(IP)
        WPLAX(I)=CORD(IP,1)
        IF(CORD(IP,1).GE.0.38) WPLAZ(nfxf+I)=ETA
        !WRITE(7,*)CORD(IP,1),ETAW(IP),WPLAZ(I)
        WPLAZ(I)=ETAW(IP)
        !IF(I.GE.NFX-1) WPLAZ(I)=ETA
    END DO

    IF(NSWIN(13,1).LE.1000) THEN   !CRITICAL SHALLOW WATER
        WPLZ(1:NFX+NFXF+NFXA)=WPLAZ(1:NFX+NFXF+NFXA)
    END IF
    
    DO I=1,NFX
	    !WRITE(7,*)CORD(NWPLAZ(I),1),WPLAZ(I)
    END DO

    !if(fr.le.0.2) CALL SMOOTHING_RECT   !滤波函数

    
    
    IF(MTROM.EQ.1) THEN
    DO I=NBPOINT+1,M+L
        IF(I.LE.NBPOINT+(NFX+NFXF+NFXA)*NFY) THEN
        WRITE(5,102)CORD(I,1),CORD(I,2)-LOFFY,ETAW(I),0,0,1
        ELSE
        WRITE(5,102)CORD(I,1),CORD(I,2),ETAW(I),0,0,1
        END IF
        !WRITE(4,102)CORD(I,1:2),SP(I-NBPOINT),0,0,1	
        !CORD(I,3)=ETAW(I)
    END DO
    ELSE
    DO I=NBPOINT+1,M+L
        WRITE(5,102)CORD(I,1:2),ETAW(I),0,0,1
        !WRITE(4,102)CORD(I,1:2),SP(I-NBPOINT),0,0,1	
        !CORD(I,3)=ETAW(I)
    END DO
    DO IP=1,NFX+NFXF+NFXA 
        WRITE(5,102)CORD(NBPOINT+(IP-1)*NFY+1,1),PCORD(NBPOINTP+(IP-1)*NFY+1,2),WPLZ(IP),0,0,1
        !WRITE(5,102)CORD(NBPOINT+(IP-1)*NFY+1,1),PCORD(NBPOINTP+(IP-1)*NFY+1,2),WPLZ(IP),0,0,1
        !WRITE(4,102)CORD(NBPOINT+(IP-1)*NFY+1,1),PCORD(NBPOINTP+(IP-1)*NFY+1,2),WFZ(IP),0,0,1
        !WRITE(4,102)CORD(I,1:2),SP(I-NBPOINT),0,0,1	
    END DO  
    END IF
    
    K3=0
    DO I=NBBLOCK+1,NBBLOCK+NFBLOCK
			WRITE(4,*)INE(I,1)-NBPOINT,INE(I,2)-NBPOINT,INE(I,9)-NBPOINT,INE(I,8)-NBPOINT
            WRITE(4,*)INE(I,8)-NBPOINT,INE(I,9)-NBPOINT,INE(I,6)-NBPOINT,INE(I,7)-NBPOINT
            WRITE(4,*)INE(I,2)-NBPOINT,INE(I,3)-NBPOINT,INE(I,4)-NBPOINT,INE(I,9)-NBPOINT
            WRITE(4,*)INE(I,4)-NBPOINT,INE(I,5)-NBPOINT,INE(I,6)-NBPOINT,INE(I,9)-NBPOINT
    END DO
    
    DO I=NBBLOCK+1,NBBLOCK+NFBLOCK
			WRITE(5,*)INE(I,1)-NBPOINT,INE(I,2)-NBPOINT,INE(I,9)-NBPOINT,INE(I,8)-NBPOINT
            WRITE(5,*)INE(I,8)-NBPOINT,INE(I,9)-NBPOINT,INE(I,6)-NBPOINT,INE(I,7)-NBPOINT
            WRITE(5,*)INE(I,2)-NBPOINT,INE(I,3)-NBPOINT,INE(I,4)-NBPOINT,INE(I,9)-NBPOINT
            WRITE(5,*)INE(I,4)-NBPOINT,INE(I,5)-NBPOINT,INE(I,6)-NBPOINT,INE(I,9)-NBPOINT
    END DO
    
    IF(MTROM.NE.1) THEN
    DO I=1,(NFXF+NFX+NFXA-1)
        ND1=I+NFPOINT
        ND2=I+1+NFPOINT
        ND3=I*NFY+1
        ND4=ND3-NFY
        WRITE(4,*)ND1,ND2,ND3,ND4
        WRITE(5,*)ND1,ND2,ND3,ND4
    END DO
    END IF
    
    IF(MTROM.EQ.1) THEN
        DO I=1,NFXF
                ND1=NFPOINT-NZS*(NFXF+1)+(I)*NZS
                ND2=ND1+NZS
                ND3=(I)*NFY+1
                ND4=ND3-NFY
			    !WRITE(4,*)ND1,ND2,ND3,ND4
                !WRITE(5,*)ND1,ND2,ND3,ND4
        END DO
    END IF

	!WRITE(*,*)'		END SURFACEWAVE'

    
   ! IF(ITTE.EQ.1) THEN 
    !IF(MOD(ITTE,3).EQ.0) THEN 
    ETAW0=ETAW1
    DO I=1,NFPOINT
        ETAW1(I)=ETAW(I+NBPOINT)
    END DO
    !WPLZ=0
   ! END IF
    
    CALL ESPL2(WPLAX(1:NFX),WPLAZ(1:NFX),NFX,CORD(1:NBPOINT,1),NBPOINT,DCORZ(1:NBPOINT))
    !WRITE(*,*)DCORZ(1:NBPOINT)
    !STOP

    !DEALLOCATE(APX,APY,APZ,AEX,AEY,APEX,AENX,APN)
	!DEALLOCATE(SH,SA,SX,SZ,SY)										
	!DEALLOCATE(ETAW)	
    !DEALLOCATE(SPX,SPY,SPZ)
    !WPLAZ=0.

	CLOSE(4)
    CLOSE(5)
	CLOSE(3)
	CLOSE(1)
    CLOSE(6)
    CLOSE(7)
    CLOSE(16)
END	SUBROUTINE 
    
SUBROUTINE SURFACEWAVE_HOB_1
!**************************************************************
!CGENERATE THE FREESURFACE POTENT AND SURFACE WAVES
!**************************************************************
	USE GREENMOD
	USE CUMOD
    CHARACTER*7 TITLE,TITLE2
    REAL:: PHIWA(200)

    WRITE(TITLE,'(F7.3)')FR    
    WRITE(TITLE2,'(I2)')ITTE
    OPEN(4,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREESURFACE_LINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(5,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\FREESURFACE_NONLIN_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(6,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\WAVEPRO_LINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')


	WRITE(*,*) 
	WRITE(*,*)'SUB SURFACEWAVE......'
	NNF=NFX*NFY
    ALLOCATE(VRIGHT(NBPOINT))
	!ALLOCATE(SP(NTPN))

    !ETAW0=ETAW

! 读入影响系数
	DO 10 I=1,NTPN
		DO 10 J=1,NTPN
		!READ(16)SA(I,J),SX(I,J),SH(I,J),SZ(I,J),SY(I,J)   !SA1船体表面诱导系数，SA自由表面诱导系数 
        !CREAD(16,*)SH(I,J),SA(I,J)
10	CONTINUE

    DO I=1,NBPOINT
        DO IP=1,NFX
            J=NBPOINT+(NFXF+IP-1)*NFY+1
            IF(SQRT(PCORD(I,1)-PCORD(J,1))**2+(PCORD(I,2)-PCORD(J,2))**2+(PCORD(I,3)-PCORD(J,3))**2.LE.1.E-5) THEN
            !WRITE(*,*)I,J,SP(I),SP(J)
            !SPX(J)=SPX(I)
            !PHIWA(IP)=SP(I)
            !WRITE(*,*)IP,PHIWA(IP)
            END IF
        END DO
    END DO

!WRITE(*,*)'FINISHED CALCULATING THE FREE SURFACE POTENTIAL'
!WRITE(*,*)'HAVING STORED IN THE FILE FREEPHI.DAT'

101	format(1x,'title="episode solid mesh"'/1x,'variables="x","Y",& 
     "Z","NX","NY","NZ"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')             
  
    !CALL PADCOEF
    
    ETAW=0.
    IF(MTROM.EQ.1) THEN
        WRITE(4,101)NFPOINT,NFBLOCK*4+NFXF
    ELSE
        WRITE(4,101)NFPOINT,NFBLOCK*4
    END IF

    DO I=NBPOINT+1,NTPN
        DO J=NBPOINT+1,NTPN
            !采用差分求导
            ETAW(I)=ETAW(I)+DMX(I,J)*SP(J)*U/G
        END DO
        !根据源分布模型求导

        ETAW(I)=SPX(I)*U/G !-0.5*(SPX(I)**2+SPY(I)**2+SPZ(I)**2)/G
    END DO   

    DO I=4-NFXF,NFX
        IP=NBPOINT+(NFXF+I-1)*NFY+1
        IP1=IP+1
        IP0=IP1+1
        IP_1=IP-NFY+1
        DYL=SQRT((CORD(IP_1,1)-CORD(IP,1))**2+(CORD(IP_1,2)-CORD(IP,2))**2)/(CORD(IP,2)-CORD(IP_1,2))&
            *(PCORD(IP,2)-CORD(IP,2))
        DXL=(CORD(IP,1)-CORD(IP_1,1))/(CORD(IP,2)-CORD(IP_1,2))&
            *(PCORD(IP,2)-CORD(IP,2))

        
        !ETA=ETAW(IP)
        !IF(CORD(IP,1).GE.0.) THEN
        ETA=(PCORD(IP,2)-CORD(IP1,2))*(PCORD(IP,2)-CORD(IP,2))/(CORD(IP0,2)-CORD(IP1,2))/(CORD(IP0,2)-CORD(IP,2))*ETAW(IP0)&
           +(PCORD(IP,2)-CORD(IP0,2))*(PCORD(IP,2)-CORD(IP,2))/(CORD(IP1,2)-CORD(IP0,2))/(CORD(IP1,2)-CORD(IP,2))*ETAW(IP1)&
           +(PCORD(IP,2)-CORD(IP0,2))*(PCORD(IP,2)-CORD(IP1,2))/(CORD(IP,2)-CORD(IP0,2))/(CORD(IP,2)-CORD(IP1,2))*ETAW(IP)
        ETA=ETAW(IP)+(ETAW(IP)-ETAW(IP1))/(CORD(IP,2)-CORD(IP1,2))*(PCORD(IP,2)-CORD(IP,2))
        !END IF

        ETA1=ETAW(IP)+(ETAW(IP)-ETAW(IP_1))/&
             SQRT((CORD(IP_1,1)-CORD(IP,1))**2+(CORD(IP_1,2)-CORD(IP,2))**2)*&
             SQRT(DXL**2+DYL**2)

        IF(I.NE.NFX) THEN
            !ETA=(PHIWA(I+1)-PHIWA(I))*(PCORD(IP+NFY,1)-PCORD(IP,1))/&
            !    ((PCORD(IP+NFY,1)-PCORD(IP,1))**2+(PCORD(IP+NFY,2)-PCORD(IP,2))**2)*U/G
            !ETA=(PHIWA(I+1)-PHIWA(I))/&
            !    SQRT((PCORD(IP+NFY,1)-PCORD(IP,1))**2+(PCORD(IP+NFY,2)-PCORD(IP,2))**2)*U/G
        END IF

        WRITE(6,*)CORD(IP,1),ETAW(IP),ETA
        !WRITE(6,*)CORD(IP,1)+DXL,ETA1,CORD(IP,1),ETA
    END DO
    
    !CALL SMOOTHING_RECT   !滤波函数

    DO I=NBPOINT+1,NTPN
        WRITE(4,102)CORD(I,1:2),ETAW(I),0,0,1
        !WRITE(4,102)CORD(I,1:2),SP(I-NBPOINT),0,0,1	
        !CORD(I,3)=ETAW(I)
    END DO           	
    102 FORMAT(6F15.6)	
    
    ETAW=0.
    IF(MTROM.EQ.1) THEN
        WRITE(5,101)NFPOINT,NFBLOCK*4+NFXF
    ELSE
        WRITE(5,101)NFPOINT,NFBLOCK*4
    END IF
    DO I=NBPOINT+1,NTPN
        DO J=NBPOINT+1,NTPN
            !采用差分求导
            ETAW(I)=ETAW(I)+DMX(I,J)*SP(J)*U/G
        END DO
        !根据源分布模型求导

        ETAW(I)=SPX(I)*U/G-0.5*(SPX(I)**2+SPY(I)**2+SPZ(I)**2)/G
        !ETAW(I)=SPX(I)*U/G
    END DO   
    
    DO J=1,NFY
        !IP=NBPOINT+(NFX+NFXF+NFXA-1)*NFY+J

        !-P1*(0.5*(VV(1)**2+VV(2)**2&    !非线性阻力
     	!		 +VV(3)**2)-U*VV(1))*DN(1)-P1*G*(SR(3))*DN(1)
        !WRITE(*,*)0.5*((-U+SPX(IP))**2+SPY(IP)**2+SPZ(IP)**2)/G+ETAW(IP)-U**2/2./G
    END DO

    DO I=1,NFX
        IP=NBPOINT+(NFXF+I-1)*NFY+1
        IP1=IP+1
        IP0=IP1+1

        ETA=ETAW(IP)
        !IF(CORD(IP,1).GE.0.) THEN
            ETA=(PCORD(IP,2)-CORD(IP1,2))*(PCORD(IP,2)-CORD(IP,2))/(CORD(IP0,2)-CORD(IP1,2))/(CORD(IP0,2)-CORD(IP,2))*ETAW(IP0)&
               +(PCORD(IP,2)-CORD(IP0,2))*(PCORD(IP,2)-CORD(IP,2))/(CORD(IP1,2)-CORD(IP0,2))/(CORD(IP1,2)-CORD(IP,2))*ETAW(IP1)&
               +(PCORD(IP,2)-CORD(IP0,2))*(PCORD(IP,2)-CORD(IP1,2))/(CORD(IP,2)-CORD(IP0,2))/(CORD(IP,2)-CORD(IP1,2))*ETAW(IP)
        !END IF
        ETA=ETAW(IP)+(ETAW(IP)-ETAW(IP1))/(CORD(IP,2)-CORD(IP1,2))*(PCORD(IP,2)-CORD(IP,2))
        
        WPLAZ(I)=ETAW(IP)
        WPLAX(I)=CORD(IP,1)
        IF(CORD(IP,1).GE.0.38) WPLAZ(I)=ETA
        !WRITE(7,*)CORD(IP,1),ETAW(IP),WPLAZ(I)
        WPLAZ(I)=ETAW(IP)
        !IF(I.GE.NFX-1) WPLAZ(I)=ETA
    END DO

    DO I=1,NFX
	    !WRITE(7,*)CORD(NWPLAZ(I),1),WPLAZ(I)
    END DO

    !if(fr.le.0.2) CALL SMOOTHING_RECT   !滤波函数

    DO I=NBPOINT+1,NTPN
        WRITE(5,102)CORD(I,1:2),ETAW(I),0,0,1
        !WRITE(4,102)CORD(I,1:2),SP(I-NBPOINT),0,0,1	
        !CORD(I,3)=ETAW(I)
    END DO           	

    K3=0
    DO I=NBBLOCK+1,NTP
			WRITE(4,*)INE(I,1)-NBPOINT,INE(I,2)-NBPOINT,INE(I,9)-NBPOINT,INE(I,8)-NBPOINT
            WRITE(4,*)INE(I,8)-NBPOINT,INE(I,9)-NBPOINT,INE(I,6)-NBPOINT,INE(I,7)-NBPOINT
            WRITE(4,*)INE(I,2)-NBPOINT,INE(I,3)-NBPOINT,INE(I,4)-NBPOINT,INE(I,9)-NBPOINT
            WRITE(4,*)INE(I,4)-NBPOINT,INE(I,5)-NBPOINT,INE(I,6)-NBPOINT,INE(I,9)-NBPOINT
    END DO

    DO I=NBBLOCK+1,NTP
			WRITE(5,*)INE(I,1)-NBPOINT,INE(I,2)-NBPOINT,INE(I,9)-NBPOINT,INE(I,8)-NBPOINT
            WRITE(5,*)INE(I,8)-NBPOINT,INE(I,9)-NBPOINT,INE(I,6)-NBPOINT,INE(I,7)-NBPOINT
            WRITE(5,*)INE(I,2)-NBPOINT,INE(I,3)-NBPOINT,INE(I,4)-NBPOINT,INE(I,9)-NBPOINT
            WRITE(5,*)INE(I,4)-NBPOINT,INE(I,5)-NBPOINT,INE(I,6)-NBPOINT,INE(I,9)-NBPOINT
    END DO

    IF(MTROM.EQ.1) THEN
        DO I=1,NFXF
                ND1=NFPOINT-NZS*(NFXF+1)+(I)*NZS
                ND2=ND1+NZS
                ND3=(I)*NFY+1
                ND4=ND3-NFY
			    WRITE(4,*)ND1,ND2,ND3,ND4
                WRITE(5,*)ND1,ND2,ND3,ND4
        END DO
    END IF

	WRITE(*,*)'		END SURFACEWAVE'

    ETAW0=ETAW1
    DO I=1,NFPOINT
        ETAW1(I)=ETAW(I+NBPOINT)
    END DO

    CALL ESPL2(WPLAX(1:NFX),WPLAZ(1:NFX),NFX,CORD(1:NBPOINT,1),NBPOINT,DCORZ(1:NBPOINT))
    !WRITE(*,*)DCORZ(1:NBPOINT)
    !STOP

    !DEALLOCATE(APX,APY,APZ,AEX,AEY,APEX,AENX,APN)
	!DEALLOCATE(SH,SA,SX,SZ,SY)										
	!DEALLOCATE(ETAW)	
    !DEALLOCATE(SPX,SPY,SPZ)
    !WPLAZ=0.

	CLOSE(4)
    CLOSE(5)
	CLOSE(3)
	CLOSE(1)
    CLOSE(6)
    CLOSE(7)
    CLOSE(16)
END	SUBROUTINE 

SUBROUTINE TEST1
!**************************************************************
!CGENERATE THE FREESURFACE POTENT AND SURFACE WAVES
!**************************************************************
	USE GREENMOD
	USE CUMOD
	REAL :: PHI(200,200),PKSI(200,200)
	REAL :: PHI2(200,200),PHIDAOSHU(200,200)

	OPEN(4,FILE='TEST.DAT')

	ALLOCATE(SP(NTPN))

!C读取求解结果
!C物面，自由面各节点的速度势
	DO I=1,NTPN
		SP(I)=SIN(3.*(CORD(I,2)+CORD(I,1)))
    END DO


!C从求解结果中获得自由面节点的速度势
!C物面非水线部分有K个节点	
	!K=NX*(NZ-1)
    K=NX*NZ

!CWRITE(*,*)'FINISHED CALCULATING THE FREE SURFACE POTENTIAL'
!CWRITE(*,*)'HAVING STORED IN THE FILE FREEPHI.DAT'
101	format(1x,'title="episode solid mesh"'/1x,'variables="x","Y",& 
     "Z","NX","NY","NZ"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral')             
  
    !CALL PADCOEF

    WRITE(4,101)NFPOINT,NFBLOCK*4

    DO I=NBPOINT+1,NTPN
        PKSI(1,1)=0
        DO J=NBPOINT+1,NTPN
            !PKSI(1,1)=PKSI(1,1)+APEX(I,J)*SP(J)
            PKSI(1,1)=PKSI(1,1)+DMY(I,J)*SP(J)
        END DO
        WRITE(4,102)CORD(I,1:2),PKSI(1,1),3.*COS(3.*(CORD(I,2)+CORD(I,1))),0,1	
    END DO           	
    102 FORMAT(6F15.6)											
		
    K3=0
    DO I=NBBLOCK+1,NTP
			WRITE(4,*)INE(I,1)-NBPOINT,INE(I,2)-NBPOINT,INE(I,9)-NBPOINT,INE(I,8)-NBPOINT
            WRITE(4,*)INE(I,8)-NBPOINT,INE(I,9)-NBPOINT,INE(I,6)-NBPOINT,INE(I,7)-NBPOINT
            WRITE(4,*)INE(I,2)-NBPOINT,INE(I,3)-NBPOINT,INE(I,4)-NBPOINT,INE(I,9)-NBPOINT
            WRITE(4,*)INE(I,4)-NBPOINT,INE(I,5)-NBPOINT,INE(I,6)-NBPOINT,INE(I,9)-NBPOINT
    END DO

	DEALLOCATE(SP)

	WRITE(*,*)'		END SURFACEWAVE'

    DEALLOCATE(APX,APY,APZ,AEX,AEY,APEX,AENX,APN)

	CLOSE(4)
	CLOSE(3)
	CLOSE(1)
END	SUBROUTINE 

SUBROUTINE MATRIX
!**************************************************************
!CGENERATE THE REAL MATRIX OF THE INTEGRATION EQUATION
!**************************************************************
	USE GREENMOD
	USE CUMOD



!C采用bin文件格式可以减小文件大小，加快读写速度 
      open(16,file='coefuse.bin',form='binary',access='sequential',&
          status='old')
      rewind(16)
	WRITE(*,*)
	WRITE(*,*)'SUB MATRIX.........'
	ALLOCATE(SH(NTPN+500,NTPN+500),SA(NTPN+500,NTPN+500),SA1(NTPN+500,NTPN+500))
	WRITE(*,*)'    DFX=',DFX,'  DFY=',DFY

	DO 10 I=1,NTPN
		DO 10 J=1,NTPN
		READ(16)SH(I,J),SA(I,J),SA1(I,J)
!!CREAD(16,*)SH(I,J),SA(I,J)
10	CONTINUE


	WRITE(*,*)'    NBPOINT=',NBPOINT
	WRITE(*,*)'    NFPOINT=',NFPOINT
	WRITE(*,*)'    NX=     ',NX
	NL=NTPN
	NR=NTPN
    L=NFPOINT


    SA=0.0
    DO I=1,NL
        SA(I,I)=0.5
    END DO


	WRITE(*,*)'    NL=     ',NL,'      NR=     ',NR

	ALLOCATE(MAR(NL,NL),VER(NL),VERR(NL))	
    MAR=0.
    VER=0.
    VERR=0.

!产生方程右边的向量

	
	DO 50 I=1,NL
		DO 51 J=1,NL
			MAR(I,J)=SA1(I,J)
51		CONTINUE
50	CONTINUE

	DO 55 I=1,NL
		VER(I)=0.0
55	CONTINUE
	
!CM为物面除水线面上的节点数，位于总节点序列的１－Ｍ
	!M=NX*(NZ-1)
    M=NX*NZ		
	DO 52 I=1,M
		VER(I)=0.0-U*VECN(I,1)
52	CONTINUE

!C物面上位于水线上的节点数，位于总节点序列的M+NFX0*NFY+1+(I-1)*NFY
!C由于在对节点法向矢量赋值时对水线上的点的法向矢量为物面的法向矢量，所以可以直接使用
	DO 53 II=1,NX
		I=M+NFX0*NFY+1+(II-1)*NFY
		!VER(I)=0.0-U*VECN(I,1)
53	CONTINUE

	OPEN(1,FILE='VER.DAT')
	DO 54 I=1,NL
		WRITE(1,*)I,VER(I)
54	CONTINUE
	CLOSE(1)	

!CASSIGN THE INITIAL VALUE OF THE RIGHT VECTOR
	OPEN(2,FILE='VERR.DAT')
	DO 60 I=1,NL
		VERR(I)=0.0
60	CONTINUE
	DO 61 I=1,NL
		DO 62 J=1,NL
			VERR(I)=VERR(I)+MAR(I,J)*VER(J)
62		CONTINUE
		!WRITE(2,*)'VERR(',I,')=',VERR(I)
61	CONTINUE
!101	FORMAT(16F12.6)
	DEALLOCATE(MAR,VER)
	!CLOSE(2)


!形成左边的矩阵	
	ALLOCATE(MAL(NL+500,NL+500),MAL1(NL+500,NL+500))	
	DO 12 I=1,NL
	DO 13 J=1,NL
		IF (I.EQ.J)THEN
			!MAL(I,J)=2*PI-SUMHIJ(I)
            !MAL(I,J)=2*PI
            MAL(I,J)=4.*PI*SE(I)
		ELSE 
			MAL(I,J)=SH(I,J)
		END IF
13	CONTINUE

12	CONTINUE

	S1=UOG*DFX*DFX   !UOG=G/U**2
!A1=10.0/(6.0*S1)
!A2=0.0-25.0/(6.0*S1)
!A3=20.0/(6.0*S1)
!A4=0.0-5.0/(6.0*S1)
		
	A1=2.0/S1         ! 1/S1=U**2/G/DX 
	A2=0.0-5.0/S1
	A3=4.0/S1
	A4=0.0-1.0/S1

	B1=0.0
	B2=0.0-2.0/S1
	B3=3.0/S1
	B4=0.0-1.0/S1


	
!对第三块：自由面速度势系数的修正：
!DO 15 I=1,NL
!DO 15 J=M+1,NL
!JJ=J-M
!!CIF(JJ.LE.NFY)THEN
!!C	MAL(I,J)=MAL(I,J)+B1*SA(I,J)+B2*SA(I,J+NFY)
!!C  &				+B3*SA(I,J+2*NFY)+B4*SA(I,J+3*NFY)
!!CELSE IF	((L-JJ).GE.(3*NFY))THEN
!!CIF	((L-JJ).GE.(3*NFY))THEN
!!C	MAL(I,J)=MAL(I,J)+A1*SA(I,J)+A2*SA(I,J+NFY)
!!C  &				+A3*SA(I,J+2*NFY)+A4*SA(I,J+3*NFY)
!!CELSE IF((L-JJ).GE.(2*NFY))THEN
!!C	MAL(I,J)=MAL(I,J)+A1*SA(I,J)+A2*SA(I,J+NFY)
!!C  &				+A3*SA(I,J+2*NFY)
!!CELSE IF((L-JJ).GE.NFY)THEN
!!C	MAL(I,J)=MAL(I,J)+A1*SA(I,J)+A2*SA(I,J+NFY)
!!CELSE 
!!C	MAL(I,J)=MAL(I,J)+A1*SA(I,J)
!!CEND IF
!C15	CONTINUE
    MAL1=0.

!CTHE SECOND COLUMN ON THE FREE SURFACE
!PHIXX(2,J)=A1*PHI(2,J)+B2*PHI(1,J)
	DO 320 I=2,2                 
	DO 321 J=1,NFY
		IM0=M+(I-1)*NFY+J   !  NFY+1——2*NFY
		IM1=M+(I-2)*NFY+J   !  1——NFY
		DO 322 IX=1,NTPN
			MAL1(IX,IM0)=MAL1(IX,IM0)+SA(IX,IM0)*A1
			MAL1(IX,IM1)=MAL1(IX,IM1)+SA(IX,IM0)*B2
322		CONTINUE
321	CONTINUE
320	CONTINUE
!  THE THIRD COLUMN
!PHIXX(3,J)=A1*PHI(3,J)+A2*PHI(2,J)+B3*PHI(1,J)
	DO 330 I=3,3
	DO 331 J=1,NFY
		IM0=M+(I-1)*NFY+J   !  2*NFY+1——3*NFY  
		IM1=M+(I-2)*NFY+J   !  NFY+1——2*NFY   
		IM2=M+(I-3)*NFY+J   !  1——NFY
		DO 332 IX=1,NTPN
			MAL1(IX,IM0)=MAL1(IX,IM0)+SA(IX,IM0)*A1
			MAL1(IX,IM1)=MAL1(IX,IM1)+SA(IX,IM0)*A2
			MAL1(IX,IM2)=MAL1(IX,IM2)+SA(IX,IM0)*B3
332		CONTINUE
331	CONTINUE
330	CONTINUE

!FROM THE THIRD COLUMN TO THE LAST COLUMN
	DO 20 I=4,NFX
	DO 20 J=1,NFY
		DO 21 JJ=1,NFY-1
			IF(FY(I-1,1).GE.FY(I,J))THEN
				J11=1
				J12=2
			ELSE IF((FY(I-1,JJ).LT.FY(I,J)).AND.&
     			(FY(I-1,JJ+1).GE.FY(I,J)))THEN
				J11=JJ
				J12=JJ+1
			ENDIF
			IF(FY(I-2,1).GE.FY(I,J))THEN
				J21=1
				J22=2
			ELSE IF((FY(I-2,JJ).LT.FY(I,J)).AND.&
     			(FY(I-2,JJ+1).GE.FY(I,J)))THEN
				J21=JJ
				J22=JJ+1
			ENDIF
			IF(FY(I-3,1).GE.FY(I,J))THEN
				J31=1
				J32=2
			ELSE IF((FY(I-3,JJ).LT.FY(I,J)).AND.&
     			(FY(I-3,JJ+1).GE.FY(I,J)))THEN
				J31=JJ
				J32=JJ+1
			ENDIF
21		CONTINUE
		D11=(FY(I-1,J12)-FY(I,J))/(FY(I-1,J12)-FY(I-1,J11))  !DY
		D12=(FY(I,J)-FY(I-1,J11))/(FY(I-1,J12)-FY(I-1,J11))
		D21=(FY(I-2,J22)-FY(I,J))/(FY(I-2,J22)-FY(I-2,J21))
		D22=(FY(I,J)-FY(I-2,J21))/(FY(I-2,J22)-FY(I-2,J21))
		D31=(FY(I-3,J32)-FY(I,J))/(FY(I-3,J32)-FY(I-3,J31))
		D32=(FY(I,J)-FY(I-3,J31))/(FY(I-3,J32)-FY(I-3,J31))
		IM0=M+(I-1)*NFY+J       !  NFY+J——3*NFY 
		IM11=M+(I-2)*NFY+J11
		IM12=M+(I-2)*NFY+J12
		IM21=M+(I-3)*NFY+J21
		IM22=M+(I-3)*NFY+J22
		IM31=M+(I-4)*NFY+J31
		IM32=M+(I-4)*NFY+J32

        !测试程序段
        !D11=1.
        !D12=1.
        !D21=1.
        !D22=1.
        !D31=1.
        !D32=1.
        !结束

		DO 22 IX=1,NTPN
			MAL1(IX,IM0)=MAL1(IX,IM0)+SA(IX,IM0)*A1
			MAL1(IX,IM11)=MAL1(IX,IM11)+SA(IX,IM0)*A2*D11
			MAL1(IX,IM12)=MAL1(IX,IM12)+SA(IX,IM0)*A2*D12
			MAL1(IX,IM21)=MAL1(IX,IM21)+SA(IX,IM0)*A3*D21
			MAL1(IX,IM22)=MAL1(IX,IM22)+SA(IX,IM0)*A3*D22
			MAL1(IX,IM31)=MAL1(IX,IM31)+SA(IX,IM0)*A4*D31
			MAL1(IX,IM32)=MAL1(IX,IM32)+SA(IX,IM0)*A4*D32
22		CONTINUE
20	CONTINUE														

    DO I=M+1,NTPN
        DO J=M+1,M+L
            MAL(I,J)=MAL(I,J)+MAL1(I,J)
            IF(ABS(MAL1(I,J)).GE.1.E-6) WRITE(2,*)I-M,J-M,MAL1(I,J),DNODE(I) 
        END DO
    END DO
     
    DO I=1,NL
        !WRITE(2,*)MAL1(I,M+1)
    END DO

	CLOSE(16)
!C释放内存
	DEALLOCATE(SH,SA,SA1,MAL1)
	WRITE(*,*)'		END MATRIX'

	RETURN
END SUBROUTINE

SUBROUTINE SMOOTHING_RECT
    USE GREENMOD
    REAL:: FIC(3,0:3)
    REAL,ALLOCATABLE:: FILT_X(:),FILT_Y(:)
    INTEGER:: IX,IX_1,IX1,IX_2,IX2,IX_3,IX3,NL

    ALLOCATE(FILT_X((NFX+NFXA+NFXF)*NFY),FILT_Y((NFX+NFXA+NFXF)*NFY))

    NL=NBPOINT
    FIC=0
    
    FIC(1,0)=2.0*1/2.0**2
    FIC(1,1)=1.0*1/2.0**2
    
    FIC(2,0)=10.0*1/2.0**4
    FIC(2,1)=4.0*1/2.0**4
    FIC(2,2)=-1.0*1/2.0**4

    FIC(3,0)=44.0*1/2.0**6
    FIC(3,1)=15.0*1/2.0**6
    FIC(3,2)=-6.0*1/2.0**6
    FIC(3,3)=1.0*1/2.0**6
	
    DO I=1,NFY
        DO J=1,NFX+NFXA+NFXF
            K=J-1
            IX=K*NFY+I
            IX_1=(K-1)*NFY+I
            IX1=(K+1)*NFY+I
            IX_2=(K-2)*NFY+I
            IX2=(K+2)*NFY+I
            IX_3=(K-3)*NFY+I
            IX3=(K+3)*NFY+I
            
            IF(J.GE.4.AND.J.LE.NFX+NFXA+NFXF-3) THEN
                FILT_X(IX)=FIC(3,0)*ETAW(NL+IX)+&
                           FIC(3,1)*ETAW(NL+IX_1)+FIC(3,1)*ETAW(NL+IX1)+&
                           FIC(3,2)*ETAW(NL+IX_2)+FIC(3,2)*ETAW(NL+IX2)+&
                           FIC(3,3)*ETAW(NL+IX_3)+FIC(3,3)*ETAW(NL+IX3)    
            ELSE IF(J.EQ.3.OR.J.EQ.NFX+NFXA+NFXF-2) THEN
                FILT_X(IX)=FIC(2,0)*ETAW(NL+IX)+&
                           FIC(2,1)*ETAW(NL+IX_1)+FIC(2,1)*ETAW(NL+IX1)+&
                           FIC(2,2)*ETAW(NL+IX_2)+FIC(2,2)*ETAW(NL+IX2)
            ELSE IF(J.EQ.2.OR.J.EQ.NFX+NFXA+NFXF-1) THEN
                FILT_X(IX)=FIC(1,0)*ETAW(NL+IX)+&
                           FIC(1,1)*ETAW(NL+IX_1)+FIC(1,1)*ETAW(NL+IX1)    
            ELSE
                FILT_X(IX)=ETAW(NL+IX)
            END IF
        END DO
    END DO

    DO I=1,NFX+NFXA+NFXF
		DO J=1,NFY    
            K=I-1
            IX=K*NFY+J
            IX_1=K*NFY+J-1
            IX1=K*NFY+J+1
            IX_2=K*NFY+J-2
            IX2=K*NFY+J+2
            IX_3=K*NFY+J-3
            IX3=K*NFY+J+3
            
            IF(J.GE.4.AND.J.LE.NFY-3) THEN
                FILT_Y(IX)=FIC(3,0)*FILT_X(IX)+&
                           FIC(3,1)*FILT_X(IX_1)+FIC(3,1)*FILT_X(IX1)+&
                           FIC(3,2)*FILT_X(IX_2)+FIC(3,2)*FILT_X(IX2)+&
                           FIC(3,3)*FILT_X(IX_3)+FIC(3,3)*FILT_X(IX3)    
            ELSE IF(J.EQ.3.OR.J.EQ.NFY-2) THEN
                FILT_Y(IX)=FIC(2,0)*FILT_X(IX)+&
                           FIC(2,1)*FILT_X(IX_1)+FIC(2,1)*FILT_X(IX1)+&
                           FIC(2,2)*FILT_X(IX_2)+FIC(2,2)*FILT_X(IX2)
            ELSE IF(J.EQ.2.OR.J.EQ.NFY-1) THEN
                FILT_Y(IX)=FIC(1,0)*FILT_X(IX)+&
                           FIC(1,1)*FILT_X(IX_1)+FIC(1,1)*FILT_X(IX1)    
            ELSE
                FILT_Y(IX)=FILT_X(IX)
            END IF

        END DO
    END DO

    DO I=1,NFPOINT
        ETAW(I+NL)=FILT_Y(I)
    END DO

    DEALLOCATE(FILT_X,FILT_Y)
END SUBROUTINE

SUBROUTINE SMOOTHING_RECT_PHI
    USE GREENMOD
    REAL:: FIC(3,0:3)
    REAL,ALLOCATABLE:: FILT_X(:),FILT_Y(:)
    INTEGER:: IX,IX_1,IX1,IX_2,IX2,IX_3,IX3,NL

    ALLOCATE(FILT_X((NFX+NFXA+NFXF)*NFY),FILT_Y((NFX+NFXA+NFXF)*NFY))

    NL=NBPOINT
    FIC=0
    
    FIC(1,0)=2.0*1/2.0**2
    FIC(1,1)=1.0*1/2.0**2
    
    FIC(2,0)=10.0*1/2.0**4
    FIC(2,1)=4.0*1/2.0**4
    FIC(2,2)=-1.0*1/2.0**4

    FIC(3,0)=44.0*1/2.0**6
    FIC(3,1)=15.0*1/2.0**6
    FIC(3,2)=-6.0*1/2.0**6
    FIC(3,3)=1.0*1/2.0**6
	
    DO I=1,NFY
        DO J=1,NFX+NFXA+NFXF
            K=J-1
            IX=K*NFY+I
            IX_1=(K-1)*NFY+I
            IX1=(K+1)*NFY+I
            IX_2=(K-2)*NFY+I
            IX2=(K+2)*NFY+I
            IX_3=(K-3)*NFY+I
            IX3=(K+3)*NFY+I
            
            IF(J.GE.4.AND.J.LE.NFX+NFXA+NFXF-3) THEN
                FILT_X(IX)=FIC(3,0)*SMPHI(NL+IX)+&
                           FIC(3,1)*SMPHI(NL+IX_1)+FIC(3,1)*SMPHI(NL+IX1)+&
                           FIC(3,2)*SMPHI(NL+IX_2)+FIC(3,2)*SMPHI(NL+IX2)+&
                           FIC(3,3)*SMPHI(NL+IX_3)+FIC(3,3)*SMPHI(NL+IX3)    
            ELSE IF(J.EQ.3.OR.J.EQ.NFX+NFXA+NFXF-2) THEN
                FILT_X(IX)=FIC(2,0)*SMPHI(NL+IX)+&
                           FIC(2,1)*SMPHI(NL+IX_1)+FIC(2,1)*SMPHI(NL+IX1)+&
                           FIC(2,2)*SMPHI(NL+IX_2)+FIC(2,2)*SMPHI(NL+IX2)
            ELSE IF(J.EQ.2.OR.J.EQ.NFX+NFXA+NFXF-1) THEN
                FILT_X(IX)=FIC(1,0)*SMPHI(NL+IX)+&
                           FIC(1,1)*SMPHI(NL+IX_1)+FIC(1,1)*SMPHI(NL+IX1)    
            ELSE
                FILT_X(IX)=SMPHI(NL+IX)
            END IF
        END DO
    END DO

    DO I=1,NFX+NFXA+NFXF
		DO J=1,NFY    
            K=I-1
            IX=K*NFY+J
            IX_1=K*NFY+J-1
            IX1=K*NFY+J+1
            IX_2=K*NFY+J-2
            IX2=K*NFY+J+2
            IX_3=K*NFY+J-3
            IX3=K*NFY+J+3
            
            IF(J.GE.4.AND.J.LE.NFY-3) THEN
                FILT_Y(IX)=FIC(3,0)*FILT_X(IX)+&
                           FIC(3,1)*FILT_X(IX_1)+FIC(3,1)*FILT_X(IX1)+&
                           FIC(3,2)*FILT_X(IX_2)+FIC(3,2)*FILT_X(IX2)+&
                           FIC(3,3)*FILT_X(IX_3)+FIC(3,3)*FILT_X(IX3)    
            ELSE IF(J.EQ.3.OR.J.EQ.NFY-2) THEN
                FILT_Y(IX)=FIC(2,0)*FILT_X(IX)+&
                           FIC(2,1)*FILT_X(IX_1)+FIC(2,1)*FILT_X(IX1)+&
                           FIC(2,2)*FILT_X(IX_2)+FIC(2,2)*FILT_X(IX2)
            ELSE IF(J.EQ.2.OR.J.EQ.NFY-1) THEN
                FILT_Y(IX)=FIC(1,0)*FILT_X(IX)+&
                           FIC(1,1)*FILT_X(IX_1)+FIC(1,1)*FILT_X(IX1)    
            ELSE
                FILT_Y(IX)=FILT_X(IX)
            END IF

        END DO
    END DO

    DO I=1,NFPOINT
        SMPHI(I+NL)=FILT_Y(I)
    END DO

    DEALLOCATE(FILT_X,FILT_Y)
END SUBROUTINE