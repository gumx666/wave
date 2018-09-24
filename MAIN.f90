PROGRAM NONLINEAR_SHIPWAVE

	USE GREENMOD
	USE CUMOD
	USE NMLMOD
    USE INPUTDATA
    
    REAL ::FFR(30)
    INTEGER ::NFR
    INTEGER ::MCONV
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2
	CHARACTER*8 TOTA_TIME1,TOTA_TIME2
	CHARACTER*9 TOTA_DATE1,TOTA_DATE2
    CHARACTER*7 TITLE
    CHARACTER*1 creturn,br
    REAL:: TIME1,TIME2
    !CALL TIME(CHAR_TIME1)
	!CALL DATE(CHAR_DATE1)
    
    WRITE(*,*)"\\\\S\\\H\\\I\\\P\\\\\\\\\\W\\\A\\\V\\\E\\\\\\"
    WRITE(*,*)"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
    WRITE(*,*)"\\\JIAODA FASTSHIP SOLVER       VERSION 1.0\\\"
    WRITE(*,*)"\\\PANEL METHOD FOR THE NONLINEAR SHIPWAVES\\\"
    WRITE(*,*)"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
    WRITE(*,*)"\\\\S\\\H\\\I\\\P\\\\\\\\\\W\\\A\\\V\\\E\\\\\\"
    WRITE(*,*)
    
!$OMP PARALLEL
!$OMP SINGLE
    NTHREAD = OMP_GET_NUM_THREADS()
!$OMP END SINGLE
!$OMP END PARALLEL

!建立结果文件夹
    !CALL NEWFOLDER

    CALL TIME(TOTA_TIME1)
	CALL DATE(TOTA_DATE1)
    
!!C公共参数	
	NML=1
	G=9.80665
      PI=4.*atan(1.)
      PIH=pi/2.0
      PI2=pi*2.0
	P1=1000.0

!!C从数据文件：input.txt中输入主要参数
!!CU,DX,DY,NX,NY,NFX

    CALL DATAIN

    OPEN(9999,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'_LOG.DAT')
    
    DMESH=0. !RD/110.

    MNKDB=1 ! NK=1    DB=2

    MDIFF=NSWIN(8,1) ! ANA=1   NUM=2
    
    MDFREE=NSWIN(2,1) ! FREE=1  FIXED=2
    
    MGRID=2 ! UNIFORM GRID=1
    
    NTHREAD=NSWIN(1,1)
    
    !PRINT *, 'WE ARE USING',NTHREAD,' THREAD(S)'
	!U=FR*SQRT(G*RL)


!!CGET THE GAUSS POINTS AND WEIGHTS
    CALL GASIN

!!CGENERATE THE POINTS AND LATTICE


    OPEN(201,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\RESISTANCE_NON.DAT')
    OPEN(202,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\SINKAGE&TRIM.DAT')
    OPEN(203,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\RESISTANCE_LIN.DAT')
    OPEN(204,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\SINKAGE&TRIM_LIN.DAT')
    
    !CALL INITIALIZATION   !初始化
    NFR=NSWIN(3,1)
    RL=PANIN(1,1)
    RB=PANIN(1,2)
    RD=PANIN(1,3)
    NITTE=NSWIN(7,1)
    MTROM=NSWIN(5,1)
    NFXF=PANIN(10,1)
    NFX=PANIN(10,2)
    NFXA=PANIN(10,3)
    NFY=PANIN(10,4)
    DARB=PANIN(11,1)
    
    WRITE(*,*)    
    WRITE(*,*)"INITIALIZATION FINISHES     :"
    WRITE(*,*)    
    WRITE(*,*)"COMPUTATION STARTS          :"
    WRITE(*,*)    
    WRITE(*,*)"THREADS USED                :",NTHREAD
        
    WRITE(*,*)'从何种文件生成网格：1——型线；2——网格'
    READ(*,*)XMDEN
    
    
    DO NIT=1,NFR     !航速循环
        
    
    !每个迭代步前初始化
    SUA=0.0 
    DFF=RD
    ETAW=0.
    ETAW0=0.
    WPLAZ=0.  
    MMESH=1
    !WRITE(TITLE,'(I2)')NIT
  
    !OPEN(111,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\MATRIX_CONDITION.DAT')
    !OPEN(211,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\RESIDUAL.DAT')
    !OPEN(211,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\RESIDUAL.DAT')
    !OPEN(311,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\CUP_TIME.DAT')
    
    !FR=FFR(NIT)
    FR=NSWIN(4,NIT)
    WRITE(TITLE,'(F7.3)')FR
    result =systemqq( 'MD '//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'')    
    WRITE(*,*)
    WRITE(*,*)"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
    WRITE(*,*)"FROUDE NUMBER, INDEX    :",FR,NIT
    WRITE(*,*)
       
    WRITE(9999,*)"*******************************************************"
    WRITE(9999,*)"          FROUDE NUMBER :                 ",FR
    WRITE(9999,*)"*******************************************************"
    
    !NITTE=1
    DO ITTE=1,NITTE   !迭代循环
      !MTROM=2 !TRAMSOM=1 
      !FR=0.1+0.01*(NIT-1)
      WRITE(9999,*)
      WRITE(9999,*)
      WRITE(9999,*)"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
      WRITE(9999,*)"         NUMBER OF ITERATION :          ",ITTE
      WRITE(9999,*)"\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\"
      
      CALL TIME(CHAR_TIME1)
	  CALL DATE(CHAR_DATE1)
	  !WRITE(*,*)'BEGIN TIME   ITE: ', CHAR_TIME1,'  ON  ',CHAR_DATE1
      
      U=FR*SQRT(G)
	  UOG=G/(U*U)
	  RLAMDA=2*PI*(FR**2)*RL

      OPEN(997,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\RESISC_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
      OPEN(998,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\SINK&TRIM_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
      OPEN(718,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\RESIDUAL.DAT')
      
      creturn = achar(13) ! generate carriage return
   
      !STOP
       
      !CALL OFFWIGLEY_9NODES
      IF(FOLDERNAME.EQ.'WIGLEY'.OR.FOLDERNAME.EQ.'wigley') THEN
      NX=PANIN(2,1)
      NZ=PANIN(2,2)
      CALL OFFWIGLEY_9NODES_NON
      ELSE
      !CALL OFFWIGLEY_3_9NODES
      !CALL OFFSHIP
      !CALL TIME(TOTA_TIME1)
      !CALL DATE(TOTA_DATE1)
      CALL CPU_TIME(TIME1)
      CALL OFFSHIP_NOL
      CALL CPU_TIME(TIME2)
      !CALL TIME(TOTA_TIME2)
	  !CALL DATE(TOTA_DATE2)
	  !WRITE(*,*)'BEGIN TIME   MATRIX: ', TOTA_TIME1,'  ON  ',TOTA_DATE1
      !WRITE(*,*)'END   TIME   ITE: ', TOTA_TIME2,'  ON  ',TOTA_TIME2
      WRITE(*,*)"TIME",TIME2-TIME1
      END IF

      CALL OFFSETIN
      WRITE(*,*)"NUMBER OF ITERATION, NTPN  :",ITTE,NTPN

      CALL SIMSTATIC_OLD
    
      
      WRITE(9999,*)"TOTAL NUMBER OF NODES AND PANEL   :        ",NTPN,NTP
      WRITE(9999,*)"NUMBER OF NODES AND PANEL ON FREE :        ",NFPOINT,NFBLOCK
      WRITE(9999,*)"NUMBER OF NODES AND PANEL ON HULL :        ",NFBLOCK,NBBLOCK
      WRITE(9999,*)
      WRITE(9999,*)
      WRITE(9999,*)"WETTED SURFACE AREA  :        ",SAREA,"m^2 "
      WRITE(9999,*)"DISPLACEMENT VOLUME  :        ",VOL,"m^3 "
      WRITE(9999,*)"DISPLACEMENT MASS    :        ",VOL*P1,"kg"
      WRITE(9999,*)"CENTER OF BUOYANCY   :        ",XBUO(1),"m ",XBUO(2),"m ",XBUO(3),"m "
      WRITE(9999,*)"CENTER OF GRAVITY    :        ",XGRA(1),"m ",XGRA(2),"m ",XGRA(3),"m "
      WRITE(9999,*)"WATER PLANE AREA     :        ",WFA,"m^2"
      WRITE(9999,*)"WATER PLANE Lyy,Lxy  :        ",LYY,LXY,"m^4"

      CALL SOLID_ANGLES
    
      CALL TIME(CHAR_TIME1)
	  CALL DATE(CHAR_DATE1)

      CALL DIFFCORD_MATRIX_NUN_1

      !CALL PHIMAX   !计算最大相角
      !STOP
      
      !WRITE(*,*)"START COEFMATRIX_PARA"
      STAGE=2
	  CALL COEFMATRIX_PARA
   
      !CALL TIME(CHAR_TIME1)
	  !CALL DATE(CHAR_DATE1)

      CALL MATRIX_NON_4
      
      CALL TIME(CHAR_TIME2)
	  CALL DATE(CHAR_DATE2)
	  WRITE(9999,*)'BEGIN TIME   MATRIX: ', CHAR_TIME1,'  ON  ',CHAR_DATE1
      WRITE(9999,*)'END   TIME   MATRIX: ', CHAR_TIME2,'  ON  ',CHAR_DATE2

      CALL WAVE_PROFILE
      CALL PRESSURE_4
      CALL BVELOCITY
      
      !WRITE(*,*)"CALL SURFACEWAVE_HOB"
      CALL SURFACEWAVE_HOB
    
      WRITE(997,*)ITTE,-CW !,-CW1
      WRITE(998,*)ITTE,-SINK/RL,-TRI !,-SINKL/RL,-TRI,-TRIL

      IF(ITTE.EQ.1) THEN
        WRITE(203,*)FR,-CW !,-CW1
        WRITE(204,*)FR,-SINK/RL,-TRI !-SINKL/RL,-TRI,-TRIL
      END IF

      DEALLOCATE(VRIGHT)
      DEALLOCATE(DMX,DMY)

      !IF(MDIFF.EQ.2) DEALLOCATE(DVMX,DVMY)
    
      DEALLOCATE(SUMHIJ)
      IF(MNKDB.EQ.2) THEN
        DEALLOCATE(DPHIS,DPHIS_X,DPHIS_Y)
        DEALLOCATE(DPHIS_XX,DPHIS_XY,DPHIS_YY)
      END IF
      !DEALLOCATE(SH,SA,SA1,SX,SZ,SY)
      DEALLOCATE(SH,SA,SA1,SX,SZ,SY)
      
        writE(*,'(A\)')creturn
        jm=30
        write( * , '(5a,f7.2,1a)')'PROGRESS:' , '[' , & !// 如您的编译器不支持，请用上方语句代替
        repeat('#' , jm ) , repeat( '.' , 30-jm ) , '] ' , 100.00 , "%"       
      
      CALL CONVER_CRET(MCONV)
      IF(MCONV.EQ.1) GOTO 99 
  
      CALL TIME(TOTA_TIME2)
	  CALL DATE(TOTA_DATE2)
      !WRITE(*,*)'END   TIME   ITE: ', CHAR_TIME2,'  ON  ',CHAR_DATE2
      
      END DO

    99CONTINUE
    WRITE(201,*)FR,-CW !,-CW1
    !WRITE(202,*)FR,-SINK/RL,-TRI !-SINKL/RL,-TRI,-TRIL
    WRITE(202,*)FR,-SINK,-TRI !-SINKL/RL,-TRI,-TRIL
    
    END DO
    
    DEALLOCATE(QG)
!!C记录程序运行时间在Time.out文件中
!	CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
	!PRINT *, 'BEGIN TIME: ', CHAR_TIME1,'  ON  ',CHAR_DATE1
    !PRINT *, 'END   TIME: ', CHAR_TIME2,'  ON  ',CHAR_DATE2
	!OPEN(11,FILE='TIME.out')
    !  WRITE(11,*)'BEGIN TIME:',CHAR_TIME1,'  ON  ',CHAR_DATE1
	!WRITE(11,*)'END   TIME:',CHAR_TIME2,'  ON  ',CHAR_DATE2
	CLOSE(11)
    CLOSE(9999)
	!WRITE(*,*)'		END ALL'
    
    
    CALL TIME(TOTA_TIME2)
	CALL DATE(TOTA_DATE2)
    WRITE(*,*)
    WRITE(*,*)
	PRINT *, 'BEGIN TIME: ', TOTA_TIME1,'  ON  ',TOTA_DATE1
    PRINT *, 'END   TIME: ', TOTA_TIME2,'  ON  ',TOTA_DATE2
  
    pause
END PROGRAM 

SUBROUTINE NEWFOLDER

    USE DFLIB  !for systemqq
    USE DFPORT !for system
    USE GREENMOD
    LOGICAL(4) result
    INTEGER(4) I, errnum


    WRITE(*,*)"请输入名称"
    READ(*,*)NAME
    !!result =systemqq( 'md e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'')
    
    !result =systemqq( 'md '//TRIM(ADJUSTL(NAME))//'')
    !WRITE(*,*)"文件夹已建立"
    !result = SYSTEMQQ('cd e:\newfolder')
    !if(result .eqv. .True.) print *, 'command ok!'

    !I = SYSTEM("dir > e:\newfolder\file.lst.txt")
    !If (I .eq. -1) then
    !errnum = ierrno( )
    !print *, 'Error ', errnum
    !end if

    !OPEN(10,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\TEST.txt')!
    !WRITE(10,*)'TEST'
    !CLOSE(10)
END SUBROUTINE


subroutine gasin()
!C***************************************************
!!!!C input the gauss integration points and weights 
!C***************************************************    
!cdp   implicit REAL (a-h,o-z)   

    USE GREENMOD
      common /gwa/ga(400),gw(400)
!c
      !open(11,file='gauss.in',status='old')  
      OPEN(11,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\gauss.in')
      rewind(11)
!!!!!C
      k=0
      do 10 i=2,21
		do 10 j=1,i
			k=k+1
			read(11,*) ai,ga(k),gw(k)
   10 continue
!!!C
      close(11)
     
      return
 end

SUBROUTINE INITIALIZATION
    USE GREENMOD
    
    CALL OFFWIGLEY_9NODES_UN


    ALLOCATE(ETAW(NBPOINT+NFPOINT+500),ETAW0(NBPOINT+NFPOINT+500),ETAW1(NBPOINT+NFPOINT+500))
    ALLOCATE(SE(NBPOINT+NFPOINT+500))
    ALLOCATE(CORD0(NBPOINT+500,3),VECN0(NBPOINT+500,3))
    ALLOCATE(QG(NBPOINT+NFPOINT+500))
    ETAW=0.
    ETAW0=0.
    ETAW1=0.
    ETAW0=0.
    WPLAZ=0. 
END SUBROUTINE

SUBROUTINE CONVER_CRET(MCONV)
    USE GREENMOD
    
    INTEGER:: MCONV
    REAL:: RESETA,RESSIN,RESTRI
    REAL:: ERRF,NEWTON_ERROR

    MCONV=0
    SUA_I=SUA
    DFF_I=DFF
    IF(FR.LE.0.7) THEN
    !IF(FR.LT.0.3) THEN
    ERRF=0.003*U**2/G
    NEWTON_ERROR=4.E-3
    ELSE
    ERRF=0.004*U**2/G
    !ERRF=0.005*U**2/G
    NEWTON_ERROR=1.
    END IF
    
    SUA=-TRI
    DFF=RD-SINK/COS(SUA)
    
    RESSIN=ABS((DFF-DFF_I)/DFF)
    RESTRI=ABS((SUA-SUA_I)/SUA)

    RESETA=0.
    SETA=0.
    DO I=NBPOINT+1,NBPOINT+(NFXF+NFX+NFXA)*NFY
        RESETA=RESETA+(ETAW0(I)-ETAW1(I))**2
        SETA=SETA+ETAW1(I)**2
    END DO
    !RESETA=MAXVAL(ABS(ETAW0(NBPOINT+1:NBPOINT+(NFXF+NFX+NFXA)*NFY)-ETAW(NBPOINT+1:NBPOINT+(NFXF+NFX+NFXA)*NFY)))
    RESETA=MAXVAL(ABS(ETAW0(1:(NFXF+NFX+NFXA)*NFY)-ETAW1(1:(NFXF+NFX+NFXA)*NFY)))
    !RESETA=MAXVAL(ABS(ETAW0-ETAW1))
    !RESETA=SQRT(RESETA)
    
    !WRITE(9999,*)"*******************************************************"
    WRITE(9999,111)ITTE,RESSIN,RESTRI,RESETA,ERRF,NEWTON_RESIDU
    !WRITE(*,*)"*******************************************************"

    WRITE(718,111)ITTE,RESSIN,RESTRI,RESETA,ERRF,NEWTON_RESIDU
    111FORMAT(' RESIDUAL(SINK,TRIM,ETAW)  ',I6,7E15.6)

    !IF(RESSUN.LE.1.E-3.AND.RESTRI.LE.5.E-2.AND.RESETA.LE.5.E-4) THEN 
    
    !IF(ITTE.GE.6.AND.RESSIN.LE.1.E-3.AND.RESTRI.LE.1.E-1.AND.RESETA.LE.ERRF.AND.NEWTON_RESIDU.LE.NEWTON_ERROR) THEN 
    IF(ITTE.GE.6.AND.RESSIN.LE.1.E-3.AND.RESTRI.LE.1.E-1.AND.RESETA.LE.ERRF) THEN
    !IF(RESSUN.LE.1.E-3.AND.RESTRI.LE.1.E-1.AND.RESETA.LE.ERRF) THEN 
    
    !IF(RESSUN.LE.1.E-3.AND.RESTRI.LE.1.E-1.AND.RESETA.LE.ERRF.AND.NEWTON_RESIDU.LE.1.E-4) THEN 
        MCONV=1
    END IF 
END SUBROUTINE

SUBROUTINE DATAIN
    USE INPUTDATA
    USE GREENMOD
    USE FOLDERMOD
    USE DFLIB  !for systemqq
    USE DFPORT !for system
    LOGICAL(4) result
    INTEGER(4) I, errnum
    
    WRITE(*,*)"PLEASE INPUT FOLDER ADDRESS :"
    READ(*,*)ADDRESS
    result =systemqq( 'md '//TRIM(ADJUSTL(ADDRESS))//'')
    !result =systemqq( 'md e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'')
    WRITE(*,*)"PLEASE INPUT SHIP NAME      :"
    READ(*,*)FOLDERNAME
    result =systemqq( 'md '//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'')
    !result =systemqq( 'md '//TRIM(ADJUSTL(NAME))//'')
    !WRITE(*,*)"文件夹已建立"
    ADDRESS1=ADDRESS
    FOLDERNAME1=FOLDERNAME
    
    
    OPEN(11,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'_PANIN.DAT')
    OPEN(12,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'_NSWIN.DAT')
    READ(11,*)
    READ(11,*)PANIN(1,1:3)           
    DO I=1,2; READ(11,*); END DO
    READ(11,*)PANIN(2,1:2)
    DO I=1,2; READ(11,*); END DO
    READ(11,*)PANIN(3,1:3)
    DO I=1,2; READ(11,*); END DO
    READ(11,*)PANIN(4,1:2)
    DO I=1,2; READ(11,*); END DO
    READ(11,*)PANIN(5,1:2)
    DO I=1,2; READ(11,*); END DO
    READ(11,*)PANIN(6,1:3)
    DO I=1,2; READ(11,*); END DO
    READ(11,*)PANIN(7,1:3)
    DO I=1,2; READ(11,*); END DO
    READ(11,*)PANIN(8,1)
    DO I=1,2; READ(11,*); END DO
    READ(11,*)PANIN(9,1:3)
    DO I=1,2; READ(11,*); END DO
    READ(11,*)PANIN(10,1:4)
    DO I=1,2; READ(11,*); END DO
    READ(11,*)PANIN(11,1:1)
    
    READ(12,*)
    READ(12,*)NSWIN(1,1)            !NUMBER OF THREADS USED
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(2,1)            !FIXED OR FREE MODEL 
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(3,1)            !
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(4,1:NSWIN(3,1)) !
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(13,1)           !WATER DEPTH    
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(5,1)            !TRANSOM 
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(6,1:3)          !XGRA
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(7,1)          
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(8,1)           
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(9,1)            !RAISED PANEL DISTANCE ON FREE SURFACE   
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(10,1)           !COLLECTION POINT SHIFTING
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(11,1)           !Y DISTANCE OF COLLECTION POINT
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(12,1)           !SMOOTHING FN 
    DO I=1,2; READ(12,*); END DO    !RELAXATION FACTOR FREE SURFACE
    READ(12,*)NSWIN(15,1)
    DO I=1,2; READ(12,*); END DO    !RELAXATION FACTOR FOR INFLUENCE COEFFICIENTS
    READ(12,*)NSWIN(16,1)
    IF(NSWIN(13,1).LE.500) THEN     !IF SHALLOW WATER
    DO I=1,2; READ(12,*); END DO
    READ(12,*)NSWIN(14,1),NSWIN(14,2)       
    END IF
    
    SCALE2=PANIN(8,1)
    
    CLOSE(11)
    CLOSE(12)
END SUBROUTINE





	


	
