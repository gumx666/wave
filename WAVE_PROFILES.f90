SUBROUTINE WAVE_PROFILE
	USE GREENMOD
	USE CUMOD
    REAL:: XI,ETA,AJ,AIWL1,AIWL
    common /gwa/ga(400),gw(400)
	DIMENSION ISM(3,2)
	DATA ISM/1,1,1,1,-1,1/

    REAL,dimension(:,:):: xxe(2,3)
    REAL,dimension(:):: sn(9),snx(9),sne(9),sr(3),dn(3),VV(3)
    REAL,dimension(:):: sn3(3),snx3(3),sne3(3),sr3(2),dn3(2)
	REAL,dimension(:,:):: a(3,3),b(3,3)
	DIMENSION IS(3),JS(3),FRES(2),FRES0(2)
	INTEGER:: L,ISYM
    CHARACTER*7 TITLE,TITLE2

    REAL ::FCORL(9,3),FCNL(9,3),FPS(3),T,TL,FAIJ(9),BETA,LOFF
    REAL ::HIJX(9),HIJN(9),HIJZ(9),HYX(9),HYN(9),HYZ(9),FHIJ(9,8) !诱导系数 
    REAL ::WX(500),WZ(500),WXL(500),WZL(500)
    REAL ::WIX(500),WIZ(500)

    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(99,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\'//TRIM(ADJUSTL(TITLE))//'\WAVE_PRO_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
    !OPEN(106,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\NEW_WPRO_LINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    !OPEN(107,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\NEW_WPRO_NONLINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(106,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\WPRO_LINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(107,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\WPRO_NONLINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
      

	!DOUBLE PRECISION A,T,D
	!WRITE(*,*) 
	!WRITE(*,*)'SUB PRESSURE......'
	NSHIP=NX*NZ
	ALLOCATE (CORL(9,3),CNL(9,3))
	!ALLOCATE (SP(NTPN))

! 读入源强

	
! 读入诱导系数
    open(16,file='coefuse.bin',form='binary',access='sequential',&
          status='old')
    rewind(16)
	!WRITE(*,*)
	!WRITE(*,*)'SUB MATRIX.........'
	!ALLOCATE(SH(NTPN,NTPN),SX(NTPN,NTPN),SZ(NTPN,NTPN),SA(NTPN,NTPN),SY(NTPN,NTPN))
	!ALLOCATE(PHIS(NTPN))
    !ALLOCATE(SPX(NTPN),SPY(NTPN),SPZ(NTPN))

    THR1=0.

	AIJX=0.0
	AIJX1=0.0
    AIJX2=0.0
    AIJX3=0.0
    AIJX4=0.0
    AIJX5=0.0

    IF(ITTE.EQ.1) THEN
        WZ0=0.
        WPLAZ=0.
    END IF

    !WRITE(106,*)WPLX(1),WPLAZ(1)
    !WRITE(107,*)WPLX(1),WPLAZ(1)
    
    !WX(1)=WPLX(1);WZ(1)=WPLAZ(1)
    !WXL(1)=WPLX(1);WZL(1)=WPLAZ(1)

    WX(1)=PCORD(NBPOINT+NFXF*NFY+1,1);WZ(1)=WPLAZ(1)
    WXL(1)=PCORD(NBPOINT+NFXF*NFY+1,1);WZL(1)=WPLAZ(1)

!********************************************************
!计算波剖面
    THR1=0.
    DO IP=1,NFX+NFXF+NFXA !-1
        FPS(1)=WPLX(IP)
        FPS(2)=WPLY(IP)+0.2*DFY !+0.2*DFY
        FPS(3)=WPLZ(IP) 
        
        !FPS(1)=PCORD(NBPOINT+(NFXF+IP-1)*NFY+1,1)
        FPS(1)=CORD(NBPOINT+(IP-1)*NFY+1,1) !-0.25*DFY
        
        IF(ITTE.LE.1) THEN
        IF(IP.GT.NFXF.AND.IP.LE.NFXF+NFX) THEN
            IF(FR.LT.0.4) THEN
                FPS(2)=PCORD(NBPOINTP+(IP-1)*NFY+1,2)+0.1*DFY
            ELSE IF(FR.GE.0.4.AND.FR.LT.0.6) THEN
                FPS(2)=PCORD(NBPOINTP+(IP-1)*NFY+1,2)+0.2*DFY
            ELSE IF(FR.GE.0.6.AND.FR.LT.0.7) THEN
                FPS(2)=PCORD(NBPOINTP+(IP-1)*NFY+1,2)+0.3*DFY
                !FPS(2)=CORD(NBPOINT+(IP-1)*NFY+1,2)
                !WRITE(*,*)CORD(NBPOINT+(IP-1)*NFY+1,2),FPS(2)
            ELSE
                FPS(2)=PCORD(NBPOINTP+(IP-1)*NFY+1,2)+0.4*DFY
            END IF
        ELSE
            FPS(2)=PCORD(NBPOINTP+(IP-1)*NFY+1,2)
        END IF
        
        !FPS(2)=CORD(NBPOINT+(NFXF+IP-1)*NFY+1,2)-0.25*DFY
        !FPS(3)=WPLAZ(IP)
        FPS(3)=WZ0(IP) 
        !FPS(3)=0.
        ELSE
            FPS(1:3)=CORD(NBPOINT+(IP-1)*NFY+1,1:3)
            FPS(2)=PCORD(NBPOINTP+(IP-1)*NFY+1,2)+0.15*DFY
        END IF
        
        
        WPLX(IP)=FPS(1)
        WPLY(IP)=FPS(2)
        WPLZ(IP)=FPS(3)

        WETAL=0.
        WETA=0.
        WETAV1=0.
        WETAV2=0.
        WETAV3=0.      

        DO KK=1,1
        IF(KK.EQ.1) THEN
            STAR=1
            ENDD=NTP !NBBLOCKW
        ELSE
            STAR=NBBLOCK+1
            ENDD=NTP
        END IF

        DO IE=STAR,ENDD
        !DO IE=1,NTP
            DO ISYM=1,2
                DO I1=1,9
                    I2=I1
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    
                    DO J1=1,3
                    IF(IE.GT.NBBLOCK) THEN
                    !IF(ABS(CORD(INE(IE,I2),2)).GE.1.E-6) THEN
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)+VECN(INE(IE,I2),J1)*LOFFZ*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM) 
                        
                        !FCORL(I1,2)=CORD(INE(IE,I2),2)*ISM(J1,ISYM)-LOFFY*ISM(J1,ISYM)
                        
                        if(fr.lT.0.3) FCORL(I1,3)=LOFFZ
                        FCORL(I1,3)=LOFFZ
                    ELSE
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM)
                    END IF
                    END DO

                    IF(INE(IE,I2).GT.NBPOINT) THEN
                        !LOFF=1.0*DFX 
                        FCORL(I1,1)=CORD(INE(IE,I2),1)-LOFFX
                    END IF
                END DO

                    CALL COFAH_9(FCORL(1:9,:),FCNL(1:9,:),FPS,FAIJ(1:9),FHIJ(1:9,:))

                    !WETAV1=0.
                    !WETAV2=0.
                    !WETAV3=0.
                    DO I1=1,9
                        I2=I1
                        IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                        IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                        
                        WETAL=WETAL+U/G*FHIJ(I1,1)*SIGM(INE(IE,I2))
                        WETA=WETA+U/G*FHIJ(I1,1)*SIGM(INE(IE,I2))!&
                          !-0.5*((FHIJ(I2,1)*SIGM(INE(IE,I1)))**2&
                          !+     (FHIJ(I2,2)*SIGM(INE(IE,I1)))**2&
                          !+     (FHIJ(I2,3)*SIGM(INE(IE,I1)))**2)/G
                        WETAV1=WETAV1+FHIJ(I1,1)*SIGM(INE(IE,I2))
                        WETAV2=WETAV2+FHIJ(I1,2)*SIGM(INE(IE,I2))
                        WETAV3=WETAV3+FHIJ(I1,3)*SIGM(INE(IE,I2))
                    END DO
                    !WETA=WETA-0.5/G*(WETAV1**2+WETAV2**2+WETAV3**2)
                    !WRITE(*,*)WETAV1,WETAV2,WETAV3
                END DO
            END DO
        END DO
                    
        WETA=WETA-0.5/G*(WETAV1**2+WETAV2**2+WETAV3**2)

        !WRITE(106,*)WPLX(IP),WETAL
        !WRITE(107,*)WPLX(IP),WETA
        WX(IP)=WPLX(IP);WZ(IP)=WETA
        !IF(ITE.GE.3) WZ(IP)=WETA
        WXL(IP)=WPLX(IP);WZL(IP)=WETAL
        
        WX(IP)=FPS(1);WZ(IP)=WETA
        WXL(IP)=FPS(1);WZL(IP)=WETAL

        !WRITE(*,*)WX(IP),WETA
    END DO       
	DEALLOCATE(CORL,CNL)
    !STOP
    !WRITE(106,*)WPLX(NFX),WPLAZ(NFX)
    !WRITE(107,*)WPLX(NFX),WPLAZ(NFX)
    
    !WX(NFX)=WPLX(NFX);WZ(NFX)=WPLAZ(NFX)
    !WXL(NFX)=WPLX(NFX);WZL(NFX)=WPLAZ(NFX)

    IF(ITTE.NE.1) THEN
    !WX(NFXF+NFX)=PCORD(NBPOINTP+(NFXF+NFX)*NFY+1,1);
    WZ(NFXF+NFX)=WPLAZ(NFXF+NFX)
    !WXL(NFXF+NFX)=PCORD(NBPOINTP+(NFXF+NFX)*NFY+1,1);
    WZL(NFXF+NFX)=WPLAZ(NFXF+NFX)

    !WX(NFXF+1)=PCORD(NBPOINT+NFXF*NFY+1,1);
    WZ(NFXF+1:NFXF+4)=WPLAZ(NFXF+1:NFXF+4)
    !WXL(NFXF+1)=PCORD(NBPOINT+NFXF*NFY+1,1);
    WZL(NFXF+1:NFXF+4)=WPLAZ(NFXF+1:NFXF+4)    
    
    END IF
    
    IF(ITTE.EQ.1) THEN   
    WZ(NFXF+NFX)=0.5*(WZ(NFXF+NFX+1)+WZ(NFXF+NFX-1))
    WZL(NFXF+NFX)=0.5*(WZL(NFXF+NFX+1)+WZL(NFXF+NFX-1))
    WZ(NFXF+1)=0.5*(WZ(NFXF+2)+WZ(NFXF))
    WZL(NFXF+1)=0.5*(WZL(NFXF+2)+WZL(NFXF))
    END IF
    !WZ(1:NFX)=WPLAZ(1:NFX)

    !POTENTIAL ON HULL 
    DO I=2,NFX-1
        !IP=(I-1)*NZ+1
        !WZ(NFXF+I)=U/G*(SP(IP+NZ)-SP(IP-NZ))/(CORD(IP+NZ,1)-CORD(IP-NZ,1))
    END DO
    
    DO IP=1,151
        WIX(IP)=-0.5+(1./150.)*(IP-1)
    END DO
    CALL ESPL2(WXL(NFXF+1:NFXF+NFX),WZL(NFXF+1:NFXF+NFX),NFX,WIX(1:151),151,WIZ(1:151))
    DO IP=1,151
        IF(IP.LE.NFX+NFXF.AND.IP.GT.NFXF) THEN
            WRITE(106,112)WIX(IP),WIZ(IP),WXL(IP),WZL(IP) 
        ELSE
            WRITE(106,112)WIX(IP),WIZ(IP)
        END IF
    END DO
    WPLZ(1:NFXF+NFX+NFXA)=WZ(1:NFXF+NFX+NFXA)
    
    
    112FORMAT(4F18.8)

    CALL ESPL2(WX(NFXF+1:NFXF+NFX),WZ(NFXF+1:NFXF+NFX),NFX,WIX(1:151),151,WIZ(1:151))
    DO IP=1,151
        IF(IP.LE.NFX+NFXF.AND.IP.GT.NFXF) THEN
            WRITE(107,112)WIX(IP),WIZ(IP),WX(IP),WZ(IP)
        ELSE
            WRITE(107,112)WIX(IP),WIZ(IP)
        END IF
    END DO

    WZ0=WZ

	!DEALLOCATE(SH,SA,SX,SZ,SY)
    !DEALLOCATE(SPX,SPY,SPZ)
	!DEALLOCATE(CORL,CNL)

    !OPEN(1111,FILE="WAVE_PROFILE.DAT")
    !DO I=1,NBPOINT
    !    IF(ABS(CORD(I,3)).LE.1.E-5) THEN
            !WRITE(1111,*)CORD(I,1),CORD(I,2),SP(I)
    !    END IF
    !END DO
    !STOP
    
	!WRITE(*,*)'		END MATRIX'
    
    CLOSE(16)
    CLOSE(17)
    CLOSE(99)
    RETURN
END SUBROUTINE
    
SUBROUTINE WAVE_PROFILE_1
	USE GREENMOD
	USE CUMOD
    REAL:: XI,ETA,AJ,AIWL1,AIWL
    common /gwa/ga(400),gw(400)
	DIMENSION ISM(3,2)
	DATA ISM/1,1,1,1,-1,1/

    REAL,dimension(:,:):: xxe(2,3)
    REAL,dimension(:):: sn(9),snx(9),sne(9),sr(3),dn(3),VV(3)
    REAL,dimension(:):: sn3(3),snx3(3),sne3(3),sr3(2),dn3(2)
	REAL,dimension(:,:):: a(3,3),b(3,3)
	DIMENSION IS(3),JS(3),FRES(2),FRES0(2)
	INTEGER:: L,ISYM
    CHARACTER*7 TITLE,TITLE2

    REAL ::FCORL(9,3),FCNL(9,3),FPS(3),T,TL,FAIJ(9),BETA,LOFF
    REAL ::HIJX(9),HIJN(9),HIJZ(9),HYX(9),HYN(9),HYZ(9),FHIJ(9,8) !诱导系数 
    REAL ::WX(500),WZ(500),WXL(500),WZL(500)
    REAL ::WIX(500),WIZ(500)

    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(99,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\'//TRIM(ADJUSTL(TITLE))//'\WAVE_PRO_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
    OPEN(106,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\NEW_WPRO_LINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(107,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\NEW_WPRO_NONLINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')


	!DOUBLE PRECISION A,T,D
	WRITE(*,*) 
	WRITE(*,*)'SUB PRESSURE......'
	NSHIP=NX*NZ
	ALLOCATE (CORL(9,3),CNL(9,3))
	!ALLOCATE (SP(NTPN))

! 读入源强

	
! 读入诱导系数
    open(16,file='coefuse.bin',form='binary',access='sequential',&
          status='old')
    rewind(16)
	WRITE(*,*)
	WRITE(*,*)'SUB MATRIX.........'
	!ALLOCATE(SH(NTPN,NTPN),SX(NTPN,NTPN),SZ(NTPN,NTPN),SA(NTPN,NTPN),SY(NTPN,NTPN))
	!ALLOCATE(PHIS(NTPN))
    !ALLOCATE(SPX(NTPN),SPY(NTPN),SPZ(NTPN))

    THR1=0.

	AIJX=0.0
	AIJX1=0.0
    AIJX2=0.0
    AIJX3=0.0
    AIJX4=0.0
    AIJX5=0.0

    IF(ITTE.EQ.1) THEN
        WZ0=0.
    END IF

    !WRITE(106,*)WPLX(1),WPLAZ(1)
    !WRITE(107,*)WPLX(1),WPLAZ(1)
    
    !WX(1)=WPLX(1);WZ(1)=WPLAZ(1)
    !WXL(1)=WPLX(1);WZL(1)=WPLAZ(1)

    WX(1)=PCORD(NBPOINT+NFXF*NFY+1,1);WZ(1)=WPLAZ(1)
    WXL(1)=PCORD(NBPOINT+NFXF*NFY+1,1);WZL(1)=WPLAZ(1)

!********************************************************
!计算波剖面
    THR1=0.
    DO IP=1,NFX !-1
        FPS(1)=WPLX(IP)
        FPS(2)=WPLY(IP)+0.2*DFY !+0.2*DFY
        FPS(3)=WPLZ(IP) 
        
        !FPS(1)=PCORD(NBPOINT+(NFXF+IP-1)*NFY+1,1)
        FPS(1)=CORD(NBPOINT+(NFXF+IP-1)*NFY+1,1) !-0.25*DFY
        IF(FR.GE.0.35) THEN
        !FPS(2)=PCORD(NBPOINTP+(NFXF+IP-1)*NFY+1,2)+0.2*DFY
        FPS(2)=PCORD(NBPOINTP+(NFXF+IP-1)*NFY+1,2)+0.1*DFY
        !FPS(2)=CORD(NBPOINT+(NFXF+IP-1)*NFY+1,2)-0.1*DFY
        ELSE
        !FPS(2)=PCORD(NBPOINTP+(NFXF+IP-1)*NFY+1,2)+0.2*DFY
        FPS(2)=PCORD(NBPOINTP+(NFXF+IP-1)*NFY+1,2)+0.1*DFY
        !FPS(2)=CORD(NBPOINT+(NFXF+IP-1)*NFY+1,2)-0.1*DFY
        END IF
        !FPS(2)=CORD(NBPOINT+(NFXF+IP-1)*NFY+1,2)-0.25*DFY
        
        !FPS(3)=WPLAZ(IP)
        FPS(3)=WZ0(IP) 
        !FPS(3)=0.

        WPLX(IP)=FPS(1)
        WPLY(IP)=FPS(2)
        WPLZ(IP)=FPS(3)

        WETAL=0.
        WETA=0.

        DO KK=1,1
        IF(KK.EQ.1) THEN
            STAR=1
            ENDD=NTP !NBBLOCKW
        ELSE
            STAR=NBBLOCK+1
            ENDD=NTP
        END IF

        DO IE=STAR,ENDD
        !DO IE=1,NTP
            DO ISYM=1,2
                DO I1=1,9
                    I2=I1
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    
                    DO J1=1,3
                    IF(IE.GT.NBBLOCK) THEN
                    !IF(ABS(CORD(INE(IE,I2),2)).GE.1.E-6) THEN
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)+VECN(INE(IE,I2),J1)*LOFFZ*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM) 
                        
                        !FCORL(I1,2)=PCORD(INE(IE,I2),2)*ISM(J1,ISYM)+LOFFY*ISYM
                        
                        if(fr.lT.0.3) FCORL(I1,3)=LOFFZ
                        FCORL(I1,3)=LOFFZ
                    ELSE
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM)
                    END IF
                    END DO

                    IF(INE(IE,I2).GT.NBPOINT) THEN
                        !LOFF=1.0*DFX 
                        FCORL(I1,1)=CORD(INE(IE,I2),1)-LOFFX
                    END IF
                END DO

                    CALL COFAH_9(FCORL(1:9,:),FCNL(1:9,:),FPS,FAIJ(1:9),FHIJ(1:9,:))

                    WETAV1=0.
                    WETAV2=0.
                    WETAV3=0.
                    DO I1=1,9
                        I2=I1
                        IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                        IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                        
                        WETAL=WETAL+U/G*FHIJ(I1,1)*SIGM(INE(IE,I2))
                        WETA=WETA+U/G*FHIJ(I1,1)*SIGM(INE(IE,I2))!&
                          !-0.5*((FHIJ(I2,1)*SIGM(INE(IE,I1)))**2&
                          !+     (FHIJ(I2,2)*SIGM(INE(IE,I1)))**2&
                          !+     (FHIJ(I2,3)*SIGM(INE(IE,I1)))**2)/G
                        WETAV1=WETAV1+FHIJ(I1,1)*SIGM(INE(IE,I2))
                        WETAV2=WETAV2+FHIJ(I1,2)*SIGM(INE(IE,I2))
                        WETAV3=WETAV3+FHIJ(I1,3)*SIGM(INE(IE,I2))
                    END DO
                    WETA=WETA-0.5/G*(WETAV1**2+WETAV2**2+WETAV3**2)
                    !WRITE(*,*)WETAV1,WETAV2,WETAV3
                END DO
            END DO
        END DO

        !WRITE(106,*)WPLX(IP),WETAL
        !WRITE(107,*)WPLX(IP),WETA
        WX(IP)=WPLX(IP);WZ(IP)=WETA
        !IF(ITE.GE.3) WZ(IP)=WETA
        WXL(IP)=WPLX(IP);WZL(IP)=WETAL
        
        WX(IP)=FPS(1);WZ(IP)=WETA
        WXL(IP)=FPS(1);WZL(IP)=WETAL

        !WRITE(*,*)WX(IP),WETA
    END DO       
	DEALLOCATE(CORL,CNL)
    !STOP
    !WRITE(106,*)WPLX(NFX),WPLAZ(NFX)
    !WRITE(107,*)WPLX(NFX),WPLAZ(NFX)
    
    !WX(NFX)=WPLX(NFX);WZ(NFX)=WPLAZ(NFX)
    !WXL(NFX)=WPLX(NFX);WZL(NFX)=WPLAZ(NFX)

    WX(NFX)=PCORD(NBPOINTP+(NFXF+NFX)*NFY+1,1);WZ(NFX)=WPLAZ(NFX)
    WXL(NFX)=PCORD(NBPOINTP+(NFXF+NFX)*NFY+1,1);WZL(NFX)=WPLAZ(NFX)

    WX(1)=PCORD(NBPOINT+NFXF*NFY+1,1);WZ(1:4)=WPLAZ(1:4)
    WXL(1)=PCORD(NBPOINT+NFXF*NFY+1,1);WZL(1:4)=WPLAZ(1:4)
    !WZ(1:NFX)=WPLAZ(1:NFX)

    DO IP=1,151
        WIX(IP)=-0.5+(1./150.)*(IP-1)
    END DO
    CALL ESPL2(WXL(1:NFX),WZL(1:NFX),NFX,WIX(1:151),151,WIZ(1:151))
    DO IP=1,151
        IF(IP.LE.NFX) THEN
            WRITE(106,112)WIX(IP),WIZ(IP),WXL(IP),WZL(IP)
            WPLZ(IP)=WZ(IP)
        ELSE
            WRITE(106,112)WIX(IP),WIZ(IP)
        END IF
    END DO

    112FORMAT(4F18.8)

    CALL ESPL2(WX(1:NFX),WZ(1:NFX),NFX,WIX(1:151),151,WIZ(1:151))
    DO IP=1,151
        IF(IP.LE.NFX) THEN
            WRITE(107,112)WIX(IP),WIZ(IP),WX(IP),WZ(IP)
        ELSE
            WRITE(107,112)WIX(IP),WIZ(IP)
        END IF
    END DO

    WZ0=WZ

	!DEALLOCATE(SH,SA,SX,SZ,SY)
    !DEALLOCATE(SPX,SPY,SPZ)
	!DEALLOCATE(CORL,CNL)

    !OPEN(1111,FILE="WAVE_PROFILE.DAT")
    !DO I=1,NBPOINT
    !    IF(ABS(CORD(I,3)).LE.1.E-5) THEN
            !WRITE(1111,*)CORD(I,1),CORD(I,2),SP(I)
    !    END IF
    !END DO
    !STOP
    
	WRITE(*,*)'		END MATRIX'
    
    CLOSE(16)
    CLOSE(17)
    CLOSE(99)
    RETURN
END SUBROUTINE

SUBROUTINE WATER_LINE
	USE GREENMOD
	USE CUMOD
    REAL:: XI,ETA,AJ,AIWL1,AIWL
    common /gwa/ga(400),gw(400)
	DIMENSION ISM(3,2)
	DATA ISM/1,1,1,1,-1,1/

    REAL,dimension(:,:):: xxe(2,3)
    REAL,dimension(:):: sn(9),snx(9),sne(9),sr(3),dn(3),VV(3)
    REAL,dimension(:):: sn3(3),snx3(3),sne3(3),sr3(2),dn3(2)
	REAL,dimension(:,:):: a(3,3),b(3,3)
	DIMENSION IS(3),JS(3),FRES(2),FRES0(2)
	INTEGER:: L,ISYM
    CHARACTER*7 TITLE,TITLE2

    REAL ::FCORL(9,3),FCNL(9,3),FPS(3),T,TL,FAIJ(9),BETA,LOFF
    REAL ::HIJX(9),HIJN(9),HIJZ(9),HYX(9),HYN(9),HYZ(9),FHIJ(9,8) !诱导系数 
    REAL ::WX(100),WZ(100),WXL(100),WZL(100)
    REAL ::WIX(200),WIZ(200)

        REAL TPCOR(500,3),UI(500)

    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(99,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\'//TRIM(ADJUSTL(TITLE))//'\WAVE_PRO_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
    OPEN(106,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\NEW_WPRO_LINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(107,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\NEW_WPRO_NONLINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')


	!DOUBLE PRECISION A,T,D
	WRITE(*,*) 
	WRITE(*,*)'SUB PRESSURE......'
	NSHIP=NX*NZ
	ALLOCATE (CORL(9,3),CNL(9,3))
	!ALLOCATE (SP(NTPN))

! 读入源强

	
! 读入诱导系数
    open(16,file='coefuse.bin',form='binary',access='sequential',&
          status='old')
    rewind(16)
	WRITE(*,*)
	WRITE(*,*)'SUB MATRIX.........'
	!ALLOCATE(SH(NTPN,NTPN),SX(NTPN,NTPN),SZ(NTPN,NTPN),SA(NTPN,NTPN),SY(NTPN,NTPN))
	!ALLOCATE(PHIS(NTPN))
    !ALLOCATE(SPX(NTPN),SPY(NTPN),SPZ(NTPN))

    THR1=0.

	AIJX=0.0
	AIJX1=0.0
    AIJX2=0.0
    AIJX3=0.0
    AIJX4=0.0
    AIJX5=0.0

    IF(ITTE.EQ.1) THEN
        WZ0=0.
    END IF

    !WRITE(106,*)WPLX(1),WPLAZ(1)
    !WRITE(107,*)WPLX(1),WPLAZ(1)
    
    !WX(1)=WPLX(1);WZ(1)=WPLAZ(1)
    !WXL(1)=WPLX(1);WZL(1)=WPLAZ(1)

    WX(1)=PCORD(NBPOINT+NFXF*NFY+1,1);WZ(1)=WPLAZ(1)
    WXL(1)=PCORD(NBPOINT+NFXF*NFY+1,1);WZL(1)=WPLAZ(1)

!********************************************************
!计算波剖面
    THR1=0.
    DO IP=1,NFX+NFXF+NFXA !-1
        
        FPS(1)=CORD(NBPOINT+(IP-1)*NFY+1,1)
        FPS(2)=PCORD(NBPOINTP+(IP-1)*NFY+1,2)
        FPS(3)=0. 
        WETA=0.
        WIX(IP-NFXF)=FPS(1)
        
        DO KK=1,1
        IF(KK.EQ.1) THEN
            STAR=1
            ENDD=NTP !NBBLOCKW
        ELSE
            STAR=NBBLOCK+1
            ENDD=NTP
        END IF

        DO IE=STAR,ENDD
        !DO IE=1,NTP
            DO ISYM=1,2
                DO I1=1,9
                    I2=I1
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    
                    DO J1=1,3
                    IF(IE.GT.NBBLOCK) THEN
                    !IF(ABS(CORD(INE(IE,I2),2)).GE.1.E-6) THEN
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)+VECN(INE(IE,I2),J1)*LOFFZ*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM) 
                        
                        !FCORL(I1,2)=PCORD(INE(IE,I2),2)*ISM(J1,ISYM)+LOFFY*ISYM
                        
                        if(fr.lT.0.3) FCORL(I1,3)=LOFFZ
                        FCORL(I1,3)=LOFFZ
                    ELSE
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM)
                    END IF
                    END DO

                    IF(INE(IE,I2).GT.NBPOINT) THEN
                        !LOFF=1.0*DFX 
                        FCORL(I1,1)=CORD(INE(IE,I2),1)-LOFFX
                    END IF
                END DO

                    CALL COFAH_9(FCORL(1:9,:),FCNL(1:9,:),FPS,FAIJ(1:9),FHIJ(1:9,:))

                    WETAV1=0.
                    WETAV2=0.
                    WETAV3=0.
                    DO I1=1,9
                        I2=I1
                        IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                        IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                        
                        WETA=WETA+FAIJ(I1)*SIGM(INE(IE,I2))!&

                    END DO
                    !WRITE(*,*)WETAV1,WETAV2,WETAV3
                END DO
            END DO
        END DO

        WFZ(IP)=WETA
        

        WRITE(*,*)IP,WETA
    END DO       
	DEALLOCATE(CORL,CNL)

    !GOTO 20
    !GOTO 10
    DO J=1,NFX
        UI(J)=1./(NFX-1)*(J-1)
    END DO
    
    DO I=1,NWAPL
        TPCOR(NWAPL-I+1,:)=CORD(MARKWAP(I),:)
        TPCOR(NWAPL-I+1,3)=SP(MARKWAP(I))
        WRITE(*,*)TPCOR(I,1),SP(MARKWAP(I))
    END DO
    WRITE(*,*)
    CALL ESPL2(TPCOR(1:NWAPL,1),TPCOR(1:NWAPL,3),NWAPL,WIX(1:NFX),NFX,WIZ(1:NFX))
    DO IP=1,NFX
        WFZ(NFXF+IP)=WIZ(IP)
        WRITE(*,*)WIX(IP),WFZ(NFXF+IP)
    END DO        

    RETURN
END SUBROUTINE

SUBROUTINE WAVE_PROFILE_OLD
	USE GREENMOD
	USE CUMOD
    REAL:: XI,ETA,AJ,AIWL1,AIWL
    common /gwa/ga(400),gw(400)
	DIMENSION ISM(3,2)
	DATA ISM/1,1,1,1,-1,1/

    REAL,dimension(:,:):: xxe(2,3)
    REAL,dimension(:):: sn(9),snx(9),sne(9),sr(3),dn(3),VV(3)
    REAL,dimension(:):: sn3(3),snx3(3),sne3(3),sr3(2),dn3(2)
	REAL,dimension(:,:):: a(3,3),b(3,3)
	DIMENSION IS(3),JS(3),FRES(2),FRES0(2)
	INTEGER:: L,ISYM
    CHARACTER*7 TITLE,TITLE2

    REAL ::FCORL(9,3),FCNL(9,3),FPS(3),T,TL,FAIJ(9),BETA,LOFF
    REAL ::HIJX(9),HIJN(9),HIJZ(9),HYX(9),HYN(9),HYZ(9),FHIJ(9,8) !诱导系数 
    REAL ::WX(500),WZ(500),WXL(500),WZL(500)
    REAL ::WIX(200),WIZ(200)

    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(99,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\'//TRIM(ADJUSTL(TITLE))//'\WAVE_PRO_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
    OPEN(106,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\NEW_WPRO_LINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(107,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\NEW_WPRO_NONLINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')


	!DOUBLE PRECISION A,T,D
	WRITE(*,*) 
	WRITE(*,*)'SUB PRESSURE......'
	NSHIP=NX*NZ
	ALLOCATE (CORL(9,3),CNL(9,3))
	!ALLOCATE (SP(NTPN))

! 读入源强

	
! 读入诱导系数
    open(16,file='coefuse.bin',form='binary',access='sequential',&
          status='old')
    rewind(16)
	WRITE(*,*)
	WRITE(*,*)'SUB MATRIX.........'
	!ALLOCATE(SH(NTPN,NTPN),SX(NTPN,NTPN),SZ(NTPN,NTPN),SA(NTPN,NTPN),SY(NTPN,NTPN))
	!ALLOCATE(PHIS(NTPN))
    !ALLOCATE(SPX(NTPN),SPY(NTPN),SPZ(NTPN))

    THR1=0.

	AIJX=0.0
	AIJX1=0.0
    AIJX2=0.0
    AIJX3=0.0
    AIJX4=0.0
    AIJX5=0.0

    IF(ITTE.EQ.1) THEN
        WZ0=0.
    END IF

    !WRITE(106,*)WPLX(1),WPLAZ(1)
    !WRITE(107,*)WPLX(1),WPLAZ(1)
    
    !WX(1)=WPLX(1);WZ(1)=WPLAZ(1)
    !WXL(1)=WPLX(1);WZL(1)=WPLAZ(1)

    !WX(1)=PCORD(NBPOINT+NFXF*NFY+1,1);WZ(1)=WPLAZ(1)
    !WXL(1)=PCORD(NBPOINT+NFXF*NFY+1,1);WZL(1)=WPLAZ(1)

!********************************************************
!计算波剖面
    THR1=0.
    DO IP=1,NFX !-1
        FPS(1)=WPLX(IP)
        FPS(2)=WPLY(IP)+0.2*DFY !+0.2*DFY
        FPS(3)=WPLZ(IP) 
        
        !FPS(1)=PCORD(NBPOINT+(NFXF+IP-1)*NFY+1,1)
        FPS(1)=CORD(NBPOINT+(NFXF+IP-1)*NFY+1,1) !-0.25*DFY
        IF(FR.GE.0.35) THEN
        !FPS(2)=PCORD(NBPOINTP+(NFXF+IP-1)*NFY+1,2)+0.2*DFY
        FPS(2)=CORD(NBPOINT+(NFXF+IP-1)*NFY+1,2)-0.2*DFY
        ELSE
        !FPS(2)=PCORD(NBPOINTP+(NFXF+IP-1)*NFY+1,2)+0.2*DFY
        FPS(2)=CORD(NBPOINT+(NFXF+IP-1)*NFY+1,2)-0.2*DFY
        END IF
        !FPS(2)=CORD(NBPOINT+(NFXF+IP-1)*NFY+1,2)-0.25*DFY
        
        !FPS(3)=WPLAZ(IP)
        FPS(3)=WZ0(IP) 
        !FPS(3)=0.

        WETAL=0.
        WETA=0.

        DO KK=1,1
        IF(KK.EQ.1) THEN
            STAR=1
            ENDD=NTP !NBBLOCKW
        ELSE
            STAR=NBBLOCK+1
            ENDD=NTP
        END IF

        DO IE=STAR,ENDD
        !DO IE=1,NTP
            DO ISYM=1,2
                DO I1=1,9
                    I2=I1
                    IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                    
                    DO J1=1,3
                    IF(IE.GT.NBBLOCK) THEN
                    !IF(ABS(CORD(INE(IE,I2),2)).GE.1.E-6) THEN
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)+VECN(INE(IE,I2),J1)*LOFFZ*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM) 
                        
                        !FCORL(I1,2)=PCORD(INE(IE,I2),2)*ISM(J1,ISYM)+LOFFY*ISYM
                        
                        if(fr.lT.0.3) FCORL(I1,3)=LOFFZ
                        FCORL(I1,3)=LOFFZ
                    ELSE
                        FCORL(I1,J1)=CORD(INE(IE,I2),J1)*ISM(J1,ISYM)
                        FCNL(I1,J1)=VECN(INE(IE,I2),J1)*ISM(J1,ISYM)
                    END IF
                    END DO

                    IF(INE(IE,I2).GT.NBPOINT) THEN
                        !LOFF=1.0*DFX 
                        FCORL(I1,1)=CORD(INE(IE,I2),1)-LOFFX
                    END IF
                END DO

                    CALL COFAH_9(FCORL(1:9,:),FCNL(1:9,:),FPS,FAIJ(1:9),FHIJ(1:9,:))

                    WETAV1=0.
                    WETAV2=0.
                    WETAV3=0.
                    DO I1=1,9
                        I2=I1
                        IF(ISYM.EQ.2.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                        IF(ISYM.EQ.4.AND.I1.NE.1.AND.(I1.NE.9)) I2=10-I1   !保证镜像后的面元节点按逆时针排列
                        
                        WETAL=WETAL+U/G*FHIJ(I1,1)*SIGM(INE(IE,I2))
                        WETA=WETA+U/G*FHIJ(I1,1)*SIGM(INE(IE,I2))!&
                          !-0.5*((FHIJ(I2,1)*SIGM(INE(IE,I1)))**2&
                          !+     (FHIJ(I2,2)*SIGM(INE(IE,I1)))**2&
                          !+     (FHIJ(I2,3)*SIGM(INE(IE,I1)))**2)/G
                        WETAV1=WETAV1+FHIJ(I1,1)*SIGM(INE(IE,I2))
                        WETAV2=WETAV2+FHIJ(I1,2)*SIGM(INE(IE,I2))
                        WETAV3=WETAV3+FHIJ(I1,3)*SIGM(INE(IE,I2))
                    END DO
                    WETA=WETA-0.5/G*(WETAV1**2+WETAV2**2+WETAV3**2)
                END DO
            END DO
        END DO

        !WRITE(106,*)WPLX(IP),WETAL
        !WRITE(107,*)WPLX(IP),WETA
        WX(IP)=WPLX(IP);WZ(IP)=WETA
        !IF(ITE.GE.3) WZ(IP)=WETA
        WXL(IP)=WPLX(IP);WZL(IP)=WETAL
        
        WX(IP)=FPS(1);WZ(IP)=WETAL
        WXL(IP)=FPS(1);WZL(IP)=WETAL

        !WRITE(*,*)WX(IP),WETA
    END DO       
	DEALLOCATE(CORL,CNL)
    !STOP
    !WRITE(106,*)WPLX(NFX),WPLAZ(NFX)
    !WRITE(107,*)WPLX(NFX),WPLAZ(NFX)
    
    !WX(NFX)=WPLX(NFX);WZ(NFX)=WPLAZ(NFX)
    !WXL(NFX)=WPLX(NFX);WZL(NFX)=WPLAZ(NFX)

    WX(NFX)=PCORD(NBPOINTP+(NFXF+NFX)*NFY+1,1);WZ(NFX)=WPLAZ(NFX)
    WXL(NFX)=PCORD(NBPOINTP+(NFXF+NFX)*NFY+1,1);WZL(NFX)=WPLAZ(NFX)

    DO IP=1,151
        WIX(IP)=-0.5+(1./150.)*(IP-1)
    END DO
    CALL ESPL2(WXL(1:NFX),WZL(1:NFX),NFX,WIX(1:151),151,WIZ(1:151))
    DO IP=1,151
        IF(IP.LE.NFX) THEN
            WRITE(106,112)WIX(IP),WIZ(IP),WXL(IP),WZL(IP)
            WPLZ(IP)=WZ(IP)
        ELSE
            WRITE(106,112)WIX(IP),WIZ(IP)
        END IF
    END DO

    112FORMAT(4F18.8)

    CALL ESPL2(WX(1:NFX),WZ(1:NFX),NFX,WIX(1:151),151,WIZ(1:151))
    DO IP=1,151
        IF(IP.LE.NFX) THEN
            WRITE(107,112)WIX(IP),WIZ(IP),WX(IP),WZ(IP)
        ELSE
            WRITE(107,112)WIX(IP),WIZ(IP)
        END IF
    END DO

    WZ0=WZ

	!DEALLOCATE(SH,SA,SX,SZ,SY)
    !DEALLOCATE(SPX,SPY,SPZ)

	WRITE(*,*)'		END MATRIX'
    
    CLOSE(16)
    CLOSE(17)
    CLOSE(99)
    RETURN
END SUBROUTINE