SUBROUTINE PRE_CONDITION(METH)
    USE GREENMOD
    USE CUMOD

    IMPLICIT NONE

    INTEGER:: METH
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*8 CHAR_DATE1,CHAR_DATE2

    REAL:: TIME1,TIME2,T1,T2
    INTEGER:: I,J,K
    INTEGER:: LK,UJ
    INTEGER:: LEVP
    REAL,DIMENSION(:),POINTER:: PL,PR,PU
    INTEGER:: HOUR1,MINI1,SECO1,HOUR2,MINI2,SECO2
    INTEGER:: M,L,N,COUT

!   原始参数
    !REAL:: ERR=1.0E-5,TAU=1.E-4      !ok1     ！优先级2
    !REAL:: ERR=1.0E-4,TAU=1.E-3       !ok2    ！优先级1
    !!REAL:: ERR=1.0E-4,TAU=1.E-2              ！优先级3
!   原始参数

    REAL:: ERR=1.0E-3,TAU=1.E-3
    
    !IF(FR.LE.0.15) THEN
    !    ERR=1.E-6
    !    TAU=1.E-5
    !END IF

    !NBPOINT=3
    !NFPOINT=3

    M=NBPOINT
    L=NFPOINT
    N=M+L
  
    !MAL(1,1)=9;MAL(1,2)=4;MAL(1,3)=4;MAL(1,4)=2;MAL(1,5)=4;MAL(1,6)=8
    !MAL(2,1)=3;MAL(2,2)=3;MAL(2,3)=12;MAL(2,4)=6;MAL(2,5)=6;MAL(2,6)=4
    !MAL(3,1)=2;MAL(3,2)=4;MAL(3,3)=5;MAL(3,4)=2;MAL(3,5)=2;MAL(3,6)=1
    !MAL(4,1)=4;MAL(4,2)=2;MAL(4,3)=1;MAL(4,4)=7;MAL(4,5)=9;MAL(4,6)=9
    !MAL(5,1)=2;MAL(5,2)=6;MAL(5,3)=4;MAL(5,4)=3;MAL(5,5)=7;MAL(5,6)=5
    !MAL(6,1)=4;MAL(6,2)=7;MAL(6,3)=8;MAL(6,4)=6;MAL(6,5)=2;MAL(6,6)=3


    IF(METH.EQ.0) THEN


    IF(ITTE.EQ.1.AND.NIT.EQ.1) THEN
    ALLOCATE(LAG(M+L+500,M+L+500),UA(M+L+500,M+L+500),&
                !DDA(M+L+500,M+L+500),DA_1(M+L+500,M+L+500),X_BAR(N+500),BAX(N+500),RG0(N+500),S_LN(N+500,N+500),S_UN(N+500,N+500),&
                !NORMJ(N+500),ROW(N+500))
                X_BAR(N+500),BAX(N+500),RG0(N+500),NORMJ(N+500),ROW(N+500))
    
    END IF
    
    DAG=0.
    DDA=0.
    LAG=0.
    UA=0.
    DA_1=0.
    NORMJ=0.

    LEVP=5
    !LEV=9999999


    CALL TIME(CHAR_TIME1)
	CALL DATE(CHAR_DATE1)
    DO I=1,N
        !DO J=1,N
            !NORMJ(I)=NORMJ(I)+MAL(I,J)**2
            !IF(ABS(MAL(I,J)).LE.ERR) LEV(I,J)=0
        !END DO
        !NORMJ(I)=SQRT(NORMJ(I))
        
        NORMJ(I)=MAXVAL(MAL(I,1:N))
    END DO

    !WRITE(*,*)SUM(KL(1:n)),SUM(KU(1:n))
    

  !DO I=1,N
    !LAG(I,I)=1.
  !  DO J=1,I-1
        !LAG(I,J)=DDA(I,J)
  !  END DO
    !LAG(I,1:I-1)=DDA(I,1:I-1)

  !  DO J=I,N
        !UA(I,J)=DDA(I,J)
  !  END DO
  !END DO 

    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
	PRINT *, 'BEGIN TIME        : ', CHAR_TIME1,'  ON  ',CHAR_DATE1
    PRINT *, 'END   TIME  PRE_PRECONDITION  : ', CHAR_TIME2,'  ON  ',CHAR_DATE2

    WRITE(*,*)"START ILUT"

    CALL TIME(CHAR_TIME1)
	CALL DATE(CHAR_DATE1)

    COUT=0


    T1=0.
    T2=0.

    UA(1,:)=MAL(1,:)
    ROW=0.
    DO I=2,N
		IF((((I)/2000)*2000).EQ.I)THEN
			WRITE(*,200)I,N !NTPN
		END IF
  200	FORMAT('ILU',3X,'6HINODE=',I8,5X,I8)
        ROW(:)=MAL(I,:)

        !CALL CPU_TIME(TIME1)

        DO K=1,I-1
            IF(ABS(ROW(K)).GT.ERR) THEN
                ROW(K)=ROW(K)/UA(K,K)
                !WRITE(*,*)"MAL(I,I)  ",UA(K,K),MAL(K,K)
                !ROW(K)=ROW(K)/MAL(K,K)
                
                IF(ABS(ROW(K)).LT.TAU*NORMJ(I)) THEN     !DROPPING RULE OF TAO
                    ROW(K)=0.
                ELSE
                    ROW(K+1:N)=ROW(K+1:N)-ROW(K)*UA(K,K+1:N)
                END IF
                
                !IF(ABS(ROW(K)).GE.TAU*NORMJ(I)) THEN     !DROPPING RULE OF TAO
                !    ROW(K+1:N)=ROW(K+1:N)-ROW(K)*UA(K,K+1:N)
                    !DO J=K+1,N 
                        !IF(UA(K,J).GE.1.E-7) ROW(J)=ROW(J)-ROW(K)*UA(K,J)
                    !    ROW(J)=ROW(J)-ROW(K)*UA(K,J)
                    !END DO
                !END IF

                !DO J=1,N
                    !LEV(I,:)=LEV(I,K)+LEV(K,:)+1
                !END DO
            END IF
        END DO

        !CALL CPU_TIME(TIME2)
        !T1=T1+TIME2-TIME1
        
        !CALL CPU_TIME(TIME1)
        DO K=1,N
            !IF(ABS(ROW(K)).LT.TAU*NORMJ(I).AND.K.NE.I&
            !   .AND.LEV(I,K).GT.LEVP) ROW(K)=0.
            
            !IF(ABS(ROW(K)).LT.TAU*NORMJ(I).AND.K.NE.I.OR.ABS(ROW(K)).LT.ERR) THEN
            IF(ABS(ROW(K)).LE.ERR) THEN
                !ROW(K)=0.
            ELSE IF(K.LE.I-1) THEN
                LAG(I,K)=ROW(K)
            ELSE IF(K.GE.I) THEN
                UA(I,K)=ROW(K)
            END IF   
        END DO
        !CALL CPU_TIME(TIME2)
        !T2=T2+TIME2-TIME1

        !PR=>ROW(1:I-1)
        !PL=>LAG(I,1:I-1)
        !PL=PR
        !PR=>ROW(I:N)
        !PU=>UA(I,I:N)
        !PU=PR
        
        !LAG(I,1:I-1)=ROW(1:I-1)
        !UA(I,I:N)=ROW(I:N)
        ROW(1:N)=0.
    END DO


    !WRITE(*,*)T1,T2

    GOTO 99

  

    99 CONTINUE

    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
	!PRINT *, 'BEGIN TIME        : ', CHAR_TIME1,'  ON  ',CHAR_DATE1
    PRINT *, 'END   TIME  ILU   : ', CHAR_TIME2,'  ON  ',CHAR_DATE2
    !WRITE(111,*) 'BEGIN TIME        : ', CHAR_TIME1,'  ON  ',CHAR_DATE1
    WRITE(111,*) 'END   TIME  ILU   : ', CHAR_TIME2,'  ON  ',CHAR_DATE2
    
    READ(CHAR_TIME1(1:2),*)HOUR1
    READ(CHAR_TIME1(4:5),*)MINI1
    READ(CHAR_TIME1(7:8),*)SECO1
    READ(CHAR_TIME2(1:2),*)HOUR2
    READ(CHAR_TIME2(4:5),*)MINI2
    READ(CHAR_TIME2(7:8),*)SECO2
    TILU=(HOUR2-HOUR1)*3600+(MINI2-MINI1)*60+SECO2-SECO1

  DO I=1,N
    LAG(I,I)=1.
    !UA(I,I:N)=DDA(I,I:N)

    !DO J=I,N
        !UA(I,J)=DDA(I,J)
    !END DO
  END DO        

  !DO I=1,N
    !DO J=1,N
    !    IF(ABS(MAL(I,J)).GT.ERR) COUT=COUT+1 
        !COUT=COUT+1 
    !END DO
  !END DO

  !DO I=1,N
    !WRITE(*,111)LAG(I,2)
  !END DO
  WRITE(*,*)

  !DO I=1,N
    !WRITE(*,111)UA(I,I)
  !END DO
  WRITE(*,*)

  !DO I=1,N
    !do j=1,n
        !MAL(I,J)=0.
        !DO K=1,N
            !MAL(I,J)=MAL(I,J)+LAG(I,K)*UA(K,J)
        !END DO
    !END DO
    !WRITE(*,111)MAL(I,1:N)
  !END DO

  !STOP

  !WRITE(*,*)"NON ZERO NUMBERS",COUT


  111format(6f15.6)  

  ELSE IF(METH.EQ.1) THEN


  X_BAR=0.
  RG0=0.


  X_BAR(1)=BAX(1)/(LAG(1,1))
  DO I=2,N
    DO J=1,I-1
        X_BAR(I)=X_BAR(I)-(LAG(I,J))*X_BAR(J)
    END DO
    X_BAR(I)=(BAX(I)+X_BAR(I))/(LAG(I,I))
  END DO

  RG0(N)=X_BAR(N)/UA(N,N)
  DO I=N-1,1,-1
    DO J=I+1,N
        RG0(I)=RG0(I)-(UA(I,J))*RG0(J)
    END DO
    RG0(I)=(X_BAR(I)+RG0(I))/(UA(I,I))
  END DO


  END IF

  IF(METH.EQ.2) THEN

    !CALL COND_MAL

    DEALLOCATE(LAG,UA,X_BAR,BAX,RG0)
  END IF

END SUBROUTINE


SUBROUTINE PRE_CONDITION_1(METH)
  USE GREENMOD
  USE CUMOD
  INTEGER:: METH

  M=NBPOINT
  L=NFPOINT
  N=M+L

  IF(METH.EQ.0) THEN

  ALLOCATE(DAG(M+L,M+L),LAG(M+L,M+L),UA(M+L,M+L),&
            DDA(M+L,M+L),DA_1(M+L,M+L),X_BAR(N),BAX(N),RG0(N))
  DAG=0.
  DDA=0.
  LAG=0.
  UA=0.
  DA_1=0.


  DO I=1,N
    !DAG(I,I)=MAL(I,I)
    IF(I.GT.NBPOINT) THEN
    DO J=1,N
        IF(ABS(DAG(I,I)).LE.ABS(MAL(I,J))) THEN
            !DAG(I,I)=MAL(I,J)
        END IF
    END DO
    END IF
     !DA_1(I,I)=1./DAG(I,I)
  END DO

  DO I=1,N
    DAG(I,I)=MAL(I,I)
  END DO

  DO I=1,N
    DAG(I,I)=1./DAG(I,I)
    DO J=I+1,N
        IF(ABS(MAL(I,J)).GE.1.E-6.AND.ABS(MAL(J,I)).GE.1.E-6) THEN
            DAG(J,J)=DAG(J,J)-MAL(J,I)*DAG(I,I)*MAL(I,J)
        END IF
    END DO
  END DO    

  DO I=1,N
    DA_1(I,I)=1./DAG(I,I)
  END DO

  DO I=1,N
    IF(I.GT.NBPOINT) THEN
    DO J=1,N
        IF(ABS(DAG(I,I)).LE.ABS(MAL(I,J))) THEN
            DAG(I,I)=MAL(I,J)
        END IF
    END DO
    END IF
  END DO

  DO I=1,N
    IF(I.GT.1) THEN
        DO J=1,I-1
            LAG(I,J)=MAL(I,J)
        END DO
     END IF
     
     IF(I.LE.N) THEN
        DO J=I+1,N
            UA(I,J)=MAL(I,J)
        END DO
     END IF
  END DO
       
  DO I=1,N
    DO J=1,N
        DDA(I,J)=0.
        DO K=1,N
            DDA(I,J)=DDA(I,J)+DA_1(I,K)*UA(K,J)
        END DO
    END DO
    DDA(I,I)=DDA(I,I)+1.
  END DO 
  
  !write(*,*)
  do i=1,n
  !write(*,111)lag(i,:)
  end do
  
  !write(*,*)
  do i=1,n
  !write(*,111)dag(i,:)
  end do
  
  !write(*,*)
  do i=1,n
  !write(*,111)ua(i,:)
  end do

  111format(5f15.6)

  ELSE IF(METH.EQ.1) THEN


  X_BAR=0.
  RG0=0.

  X_BAR(1)=BAX(1)/(DAG(1,1)+LAG(1,1))
  DO I=2,N
    DO J=1,I-1
        X_BAR(I)=X_BAR(I)-(DAG(I,J)+LAG(I,J))*X_BAR(J)
    END DO
    X_BAR(I)=(BAX(I)+X_BAR(I))/(DAG(I,I)+LAG(I,I))
  END DO

  RG0(N)=X_BAR(N)/DDA(N,N)
  DO I=N-1,1,-1
    DO J=I+1,N
        RG0(I)=RG0(I)-(DDA(I,J))*RG0(J)
    END DO
    RG0(I)=(X_BAR(I)+RG0(I))/(DDA(I,I))
  END DO

  END IF

  IF(METH.EQ.2) THEN

    DEALLOCATE(DAG,LAG,UA,DDA,DA_1,X_BAR,BAX,RG0)
  END IF

END SUBROUTINE

SUBROUTINE MATRIX_PREGMRES
!**************************************************************
!	GENERATE AND SOLVE THE INTEGRATION EQUATION
!**************************************************************


  USE GREENMOD
  USE CUMOD

  implicit none 

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF


  REAL TPR(30),HGIJ,HGIJ1,TOL1,RNORM1,BNORM,RNORM,TPVE
  INTEGER K,THR1,THR2,CP,ITR,ITER,BM,METH,i,j,M,L,KSTART,ITEROUT,N
  REAL,PARAMETER ::ITERMAX=4000,TOL=5.0E-6 !TOL=1.0E-6
  REAL,ALLOCATABLE:: TPM(:)
  CHARACTER*8 CHAR_TIME1,CHAR_TIME2
  CHARACTER*8 CHAR_DATE1,CHAR_DATE2
  REAL:: TIME1,TIME2
  INTEGER:: HOUR1,MINI1,SECO1,HOUR2,MINI2,SECO2

  OPEN(66,FILE="TEMP.DAT")

  DO I=1,NTPN
    DO J=1,NTPN
        !WRITE(66,*)I,J,MAL(I,J)
    END DO
  END DO

  BM=NTPN/2. !200
  BM=2000
  BM=35
  write(*,*)"bm",bm
  M=NBPOINT
  L=NFPOINT
  N=M+L !+L

    OPEN(1,FILE='PHI.DAT')

!  NTPN+NFPOINT IN SINGLE CPU
20 IF(MOD(N,NTHREAD).NE.0) THEN
     DO I=1,NTHREAD-1
       IF(MOD(N,NTHREAD).EQ.I) THEN
         DO J=1,NTHREAD-1
            NPS(J)=(N-I)/NTHREAD
         END DO
         NPS(NTHREAD)=NPS(NTHREAD-1)+I
         EXIT
       END IF
     END DO
   ELSE 
     DO I=1,NTHREAD
       NPS(I)=(N)/NTHREAD
     END DO
   END IF

  !WRITE(*,*)"M+L",M+L
  ALLOCATE(CI(M+L,M+L),APG(M+L,M+L),AG(M+L,M+L),WI(M+L,3),BV(M+L,BM),BH(BM+1,BM))
  ALLOCATE(BB(M+L),RG(M+L),PG(M+L),BNORM1(NTHREAD),GG(BM+1),YG(BM+1),CG(BM),SG(BM))
  ALLOCATE(TPM(N))

  DO I=1,N
    DO J=I+1,N
        IF(ABS(MAL(I,I)).LE.ABS(MAL(J,I))) THEN
            !TPM(1:N)=MAL(I,1:N)
            !MAL(I,1:N)=MAL(J,1:N)
            !MAL(J,1:N)=TPM(1:N)

            !TPVE=VERR(I)
            !VERR(I)=VERR(J)
            !VERR(J)=TPVE
        END IF
    END DO
  END DO


  !VERR=0
  TOL1=1
  ITEROUT=0
  KSTART=0
  CI=0
  CG=0
  SG=0

  METH=0               
  CALL PRE_CONDITION(METH)
  DO I=1,N
    !WRITE(66,*)MAL(I,I),DAG(I,I),DA_1(I,I)
  END DO
  WRITE(*,*)"OK"

! 预条件结束
    CALL TIME(CHAR_TIME1)
	CALL DATE(CHAR_DATE1)


  CALL OMP_SET_NUM_THREADS(NTHREAD)

!$OMP PARALLEL PRIVATE(CP,THR1,THR2,I,J,K)
    
    CP=OMP_GET_THREAD_NUM()+1

    !WRITE(*,*)"CP=",CP,NPS(CP)

    THR1=SUM(NPS(1:CP))-NPS(CP) 
    THR2=SUM(NPS(1:CP))
    
    !WRITE(*,*)"K1=",K1

!   GMRES MUDULE

    DO I=1,NPS(CP)
      !CI(THR1+I,THR1+I)=1.0 !/A(I,I)
      !K=SUM(NPS(1:CP))-NPS(CP)+I       !SERIAL NUMBER
      !CI(THR1+I,THR1+I)=1.0/MAL(THR1+I,THR1+I)   !预处理
      CI(THR1+I,THR1+I)=1.0
    END DO

!   加预设条件   
    
    DO I=1,NPS(CP)
      DO J=1,N
        !APG(THR1+I,J)=0.
        DO K=1,N
            !APG(THR1+I,J)=APG(THR1+I,J)+CI(THR1+I,K)*MAL(K,J)
        END DO
      END DO
    END DO
    
    
    !AG(THR1+1:THR2,:)=APG(THR1+1:THR2,:)
    AG(THR1+1:THR2,:)=MAL(THR1+1:THR2,:)
    APG(THR1+1:THR2,:)=MAL(THR1+1:THR2,:)


    DO I=1,NPS(CP)
      !APG(1,THR1+I)=0.
      DO J=1,N
        !APG(1,THR1+I)=APG(1,THR1+I)+CI(THR1+I,J)*VERR(J)
      END DO
    END DO


!$OMP END PARALLEL

    !BB(THR1+1:THR2)=APG(1,THR1+1:THR2)
    BB(1:N)=VERR(1:N)

    BNORM1(1)=0
    DO I=1,N !NPS(CP)
      BNORM1(1)=BNORM1(1)+BB(I)*BB(I)
    END DO
!   向量的模 
    !CALL FNORM(BB(K1:K2),BB(K1:K2),N,1,BNORM1(CP))


	!BNORM=SQRT(SUM(BNORM1(1:NTHREAD)))
    BNORM=BNORM1(1)

    IF(NIT.EQ.1) THEN
      !QG=0.
    END IF
    QG=0
    !RG=BB

    RG=0.
    DO I=1,N
      DO J=1,N
        !RG(I)=RG(I)+BB(I)-AG(I,J)*QG(J)
        RG(I)=RG(I)-AG(I,J)*QG(J)
      END DO
      RG(I)=RG(I)+BB(I)
    END DO


    RNORM1=0.
    DO I=1,N
      RNORM1=RNORM1+RG(I)*RG(I)    !改动前
      !RNORM1=RNORM1+ABS(RG(I))
    END DO

!   预条件
    METH=1               
    BAX=RG
    CALL PRE_CONDITION(METH)
    
    DO I=1,N
        !WRITE(66,*)X_BAR(I),BAX(I),UA(I,I)
    END DO
    RG=RG0

    !STOP

    IF(SQRT(RNORM1).LE.TOL) GOTO 109    !改动前
    !IF(RNORM1/N.LE.TOL) GOTO 109    

 88   CONTINUE


    RNORM1=0
    DO I=1,N
      RNORM1=RNORM1+RG(I)*RG(I)
    END DO
      
!   SET UP V^1 AND G^0


    !WRITE(*,*)'K1',K1,NPS(CP) 

	
    !CALL FNORM(R,R,N,NS,RNORM1)
!C	WRITE(6,*) RNORM1
	RNORM=SQRT(RNORM1)
	
    DO I=1,N
      PG(I)=RG(I)/RNORM
	END DO

    GG=0

    GG(1)=RNORM
!   BEGIN ITERATION
    ITER=1

    !CG=0.;SG=0.

!	BWRITE(66,*) 'BETA',G(1)
    DO WHILE((ITER.LE.BM).AND.(TOL1.GT.TOL).AND.(ITEROUT*BM+ITER.LE.ITERMAX))
    !DO WHILE((ITER.LE.BM).AND.(TOL1.GT.TOL))

!$OMP PARALLEL PRIVATE(CP,THR1,THR2,I,J)
    
    CP=OMP_GET_THREAD_NUM()+1

    THR1=SUM(NPS(1:CP))-NPS(CP) 
    THR2=SUM(NPS(1:CP))

    !WRITE(*,*)"ITER",ITER

      !DO CP=1,NCUP   !CPU NUMBER

!   SAVE  P AS THE V(ITER)
        !K1=SUM(NPS(1:CP))-NPS(CP)
        !K2=SUM(NPS(1:CP))

        BV(THR1+1:THR2,ITER)=PG(THR1+1:THR2)
        !CALL EQUAL(BV(:,1,ITER),P,N,NS)
 !  WRITE(66,*)'ITER',ITER
!	WRITE(66,*) ((BV(I1,J1,ITER),I1=1,N),J1=1,NS)

!   FORM  AV(ITER)	 
        APG(1,THR1+1:THR2)=0
        DO I=1,NPS(CP)
          DO J=1,N
            APG(1,THR1+I)=APG(1,THR1+I)+AG(THR1+I,J)*PG(J)
          END DO
        END DO

!$OMP BARRIER 

	  !CALL MATAV(A,P,AP,N,NS)  !A*P
!	IF(ITER.EQ.1) WRITE(*,*) 'AP', AP(3,1)
!   INITIALIZE V^(ITER+1) WITH AV(ITER)
        PG(THR1+1:THR2)=APG(1,THR1+1:THR2)
        !WI(K1+1:K2,1)=APG(1,K1+1:K2)
     !GOTO 11
! MAKE V^(ITER+1) ORTHOGONAL TO V^(I),I<=ITER
 
!$OMP BARRIER 

!$OMP SINGLE 
!   预条件
    METH=1               
    BAX=PG
    CALL PRE_CONDITION(METH)
    PG=RG0
!$OMP END SINGLE 
     
        DO  J=1,ITER
! 1)正交化
	   !CALL FNORM(BV(1,1,J),AP,N,NS,HIJ)
           
!$OMP SINGLE  
            HGIJ=0.    
            DO I=1,N
              !HGIJ=HGIJ+BV(I,J)*APG(1,I)
              HGIJ=HGIJ+BV(I,J)*PG(I)
              !HGIJ=HGIJ+BV(1,J)*APG(1,I)
              !WRITE(66,*)"ABS(GG(ITER+1))",bV(I,J),PG(I)
            END DO
            BH(J,ITER)=HGIJ

     !WRITE(*,*)"ABS(GG(ITER+1))",V(I,J),PG(I)
!$OMP END SINGLE 
          !TPR(CP)=0
          !DO I=1,NPS(CP)
            !TPR(CP)=TPR(CP)+BV(K1+I,J)*PG(1,K1+I)
          !END DO      
         !DO I=1,N
           !WI(I,NS,J1+1)=WI(I,NS,J1)-BV(I,NS,J)*TPR
         !END DO
         !WRITE(*,*)SQRT(RNORM/RNORM1)

          DO  I=1,NPS(CP)
		    !P(I,K)=P(I,K)-TPR*BV(I,K,J)
            PG(THR1+I)=PG(THR1+I)-HGIJ*BV(THR1+I,J)
	      ENDDO 
	      !BH(J,ITER)=TPR

          !BH(J,ITER)=HGIJ


!$OMP BARRIER

          !END DO 
          !CALL FNORM(WI(:,:,J1),WI(:,:,J1),N,NS,RNORM)
          !CALL FNORM(WI(:,:,J1+1),WI(:,:,J1+1),N,NS,RNORM1)
          !IF(SQRT(RNORM/RNORM1).LE.TH) EXIT         
	    END DO 
        !WRITE(*,*)"BH(J,ITER)",BH(J,ITER)

!$OMP END PARALLEL 


       !WRITE(*,*)"PG",PG(1:N)

       RNORM1=0
       DO I=1,N
         RNORM1=RNORM1+PG(I)*PG(I)
       END DO     

       !PG(K1+1:K2)=APG(1,K1+1:K2)

!NORMALIZE V^(ITER+1)
! 2)标准化
	 !CALL FNORM(P,P,N,NS,RNORM1)
	 !RNORM =REAL(RNORM1)
     !RNORM=SQRT(RNORM1)

       DO  I=1,N
		  PG(I)=PG(I)/SQRT(RNORM1)
	   END DO

     BH(ITER+1,ITER)=SQRT(RNORM1) 
!	 WRITE(66,*)'H OLD BEFORE ', (BH(I,ITER),I=1,ITER+1)    

! APPLY ROTATIONS TO NEW H COLUMN

! 3)更新V和H
     DO I=1,ITER-1
	   HGIJ=BH(I,ITER)
	   HGIJ1=BH(I+1,ITER)
	   BH(I,ITER)=HGIJ*CG(I)+HGIJ1*SG(I)
	   BH(I+1,ITER)=-HGIJ*SG(I)+HGIJ1*CG(I)
	 ENDDO
!	 WRITE(66,*)'H OLD AFTER ', (BH(I,ITER),I=1,ITER)  
!    COMPUTE NEW ROTATIONS
     HGIJ=BH(ITER,ITER)
	 HGIJ1=BH(ITER+1,ITER) 	   
     
     CALL GIVENS(HGIJ,HGIJ1,CG(ITER),SG(ITER))
!    WRITE(66,*) 'GIVENS',C(ITER),S(ITER)
!    APPLY NEW ROTATIONS
     BH(ITER,ITER)= HGIJ*CG(ITER)+HGIJ1*SG(ITER)
     BH(ITER+1,ITER)=0.
!    WRITE(66,*)'H NEW AFTER ', (BH(I,ITER),I=1,ITER+1) 

     !WRITE(*,*)"ABS(GG(ITER+1))",HGIJ
     HGIJ=GG(ITER)
     !WRITE(*,*)"ABS(GG(ITER+1))",HGIJ
     GG(ITER)=CG(ITER)*HGIJ
	 GG(ITER+1)=-SG(ITER)*HGIJ
     TOL1=ABS(GG(ITER+1))/BNORM
!    WRITE(66,*)'G  AFTER ', (G(I ),I=1,ITER+1) 
!C   TOL1=CABS(G(ITER+1))
!    WRITE(6,*) 'ITER',ITER
!	 WRITE(6,*) 'TOL1',TOL1
!	 WRITE(66,*) 'TOLERANCE',TOL1
!    WRITE(151,*) ITEROUT*M+ITER,LOG10(TOL1)
	 ITER=ITER+1
     
        !WRITE(66,*) 'TOL1',TOL1
        WRITE(*,*)ITER,TOL1
        !WRITE(211,*)ITER,TOL1

    END DO


    ITER=ITER-1
    !WRITE(66,*)
    !WRITE(66,*) 'TOL1',BM,TOL1
 
!   COMPUTE  LEAST SQUARE  SOLUTION
    DO I=1,ITER
	  YG(I)=GG(I)
	ENDDO
!	WRITE(66,*)'G'
!	WRITE(66,*) (G(I),I=1,ITER)
	DO I=ITER ,1,-1 
 	   YG(I)=YG(I)/BH(I,I)   
 	   DO J=I-1,1,-1
          YG(J)=YG(J)-BH(J,I)*YG(I)
       ENDDO
    ENDDO
!	WRITE(66,*) 'Y'
!	WRITE(66,*)(Y(J),J=1,ITER)
!	 WRITE(66,*)'BV'
	
    !DO  K=1,ITER
!	 WRITE(66,*)'ITER',K
	
!       WRITE(66,*)(BV(I,1,K),I=1,N)
	!ENDDO

!    ACQUIRE ITERATION SOLUTION
	 DO I=1,N
	   DO K=1,ITER
         QG(I)=QG(I)+YG(K)*BV(I,K)
	   ENDDO
	 ENDDO
!	 WRITE(66,*) 'Q'
!    WRITE(66,*) (Q(I,1),I=1,N)

     IF((TOL1.LT.TOL).AND.(ITEROUT*BM+ITER.LE.ITERMAX)) THEN
     !IF((TOL1.LT.TOL)) THEN
!      ITERATIONS END
	 GOTO 98
	 ENDIF
     
     KSTART=0
	 IF((TOL1.GT.TOL).AND.(ITEROUT*BM+ITER.LT.ITERMAX)) KSTART=1
     !IF((TOL1.GT.TOL)) KSTART=1

     IF( (ITEROUT*BM+ITER.GT.ITERMAX).AND.(TOL1.GT.TOL) ) THEN
     
     !IF( (TOL1.GT.TOL) ) THEN
!	   WRITE(6,*) 'WARNING NO CONVERGENCE  AFTER ' 
!	   WRITE(6,*) ITEROUT*M+ITER,'ITERATIONS'
!	   WRITE(6,*) 'TOLRENCE  IS',TOL1
!C	   WRITE(*,*) 'ITEROUT',ITEROUT
       WRITE(*,*) 'TOLRENCE  IS',TOL1
       GOTO  109
	 ENDIF

    !WRITE(*,*) 'TOLRENCE  IS',TOL1

	IF(KSTART.EQ.1) THEN

	ITEROUT=ITEROUT+1
!C  　    
!       WRITE(6,*) 'ITEROUT',ITEROUT     

    APG=0
    DO I=1,N
      DO J=1,N
        APG(1,I)=APG(1,I)+AG(I,J)*QG(J)
      END DO
    END DO
                                                                                                                                                                                                                                                          
    !CALL  MATAV(A,Q,AP,N,NS)
    DO  I=1,N
	  RG(I)=BB(I)-APG(1,I)
	END DO

!   预条件
    METH=1               
    BAX=RG
    CALL PRE_CONDITION(METH)
    RG=RG0

!   ITERATIONS BEGIN AGAIN
      
	GOTO 88
    ENDIF
  98  CONTINUE
!      WRITE(6,*) 'CONGRATULATIONS CONVERGENCE ACQUIRE ' 
!	   WRITE(6,*)    'AFTER', ITEROUT*M+ITER,' ITERATIONS'
	  !WRITE(*,*) 'TOLRENCE  IS',TOL1  

  99  CONTINUE 
    BB=QG

	ITR=ITEROUT*BM+ITER
	!CLOSE(66)
    WRITE(111,*) 'ITERATIONS     ',ITR
    WRITE(111,*) 

  109  CONTINUE
  
  WRITE(*,*)"ITR=   ",ITR,TOL1     
        
  DO I=1,M+L !+L
     IF(I.LE.NBPOINT) THEN
        !PHIS(I)=QG(I)
        !WRITE(*,*)X(I)
     ELSE
        !PHIS_N(I)=QG(I)
     END IF
     WRITE(1,*)QG(I)
     !WRITE(*,*)I,QG(I)
  END DO
  
  !STOP

  !WRITE(*,*)QG(30)

  !GMRES 结束
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
	WRITE(111,*)'BEGIN TIME        : ', CHAR_TIME1,'  ON  ',CHAR_DATE1
    WRITE(111,*)'END   TIME  SOLVER: ', CHAR_TIME2,'  ON  ',CHAR_DATE2
    READ(CHAR_TIME1(1:2),*)HOUR1
    READ(CHAR_TIME1(4:5),*)MINI1
    READ(CHAR_TIME1(7:8),*)SECO1
    READ(CHAR_TIME2(1:2),*)HOUR2
    READ(CHAR_TIME2(4:5),*)MINI2
    READ(CHAR_TIME2(7:8),*)SECO2
    TGMR=(HOUR2-HOUR1)*3600+(MINI2-MINI1)*60+SECO2-SECO1

    WRITE(311,*)FR,TILU,TGMR

  !预条件释放内存
    METH=2      
    CALL PRE_CONDITION(METH)

!   END GMRES MUDULE

14	CLOSE(16)   
    CLOSE(1)
!C 释放内存
  DEALLOCATE(CI,APG,AG,WI,BV,BH)
  DEALLOCATE(BB,RG,PG,BNORM1,GG,YG,CG,SG)
  DEALLOCATE(TPM)
END SUBROUTINE



SUBROUTINE MATRIX_GMRES
!**************************************************************
!	GENERATE AND SOLVE THE INTEGRATION EQUATION
!**************************************************************
  USE GREENMOD
  USE CUMOD

#IFDEF _OPENMP 
   INCLUDE 'OMP_LIB.H'  !NEEDED FOR OMP_GET_NUM_THREADS()
#ENDIF

   
  REAL TPR(30),HGIJ,HGIJ1,TOL1
  INTEGER K,THR1,THR2,CP,ITR,ITER,BM,KSTART
  PARAMETER (ITERMAX=15000,TOL=1.0E-6)
! 产生入射势
! ALLOCATE(PHIW(NTPN),PHIW_N(NTPN),PHIW_T(NTPN))
  
  BM=2500
  M=NBPOINT
  L=NFPOINT
  N=M+L !+L

    OPEN(1,FILE='PHI.DAT')

!  NTPN+NFPOINT IN SINGLE CPU
20 IF(MOD(N,NTHREAD).NE.0) THEN
     DO I=1,NTHREAD-1
       IF(MOD(N,NTHREAD).EQ.I) THEN
         DO J=1,NTHREAD-1
            NPS(J)=(N-I)/NTHREAD
         END DO
         NPS(NTHREAD)=NPS(NTHREAD-1)+I
         EXIT
       END IF
     END DO
   ELSE 
     DO I=1,NTHREAD
       NPS(I)=(N)/NTHREAD
     END DO
   END IF

  !WRITE(*,*)"M+L",M+L
  ALLOCATE(CI(M+L,M+L),APG(M+L,M+L),AG(M+L,M+L),WI(M+L,3),BV(M+L,BM),BH(BM+1,BM))
  ALLOCATE(BB(M+L),RG(M+L),PG(M+L),BNORM1(NTHREAD),GG(BM+1),YG(BM+1),CG(BM),SG(BM))
  !ALLOCATE(QG(M+L))

  DA=0.
  DDA=0.
  LA=0.
  UA=0.

  !VERR=0
  TOL1=1
  ITEROUT=0
  KSTART=0
  CI=0

  

  WRITE(*,*)"OK"

  CALL OMP_SET_NUM_THREADS(NTHREAD)

!$OMP PARALLEL PRIVATE(CP,THR1,THR2,I,J,K,SB2)
    
    CP=OMP_GET_THREAD_NUM()+1

    !WRITE(*,*)"CP=",CP,NPS(CP)

    THR1=SUM(NPS(1:CP))-NPS(CP) 
    THR2=SUM(NPS(1:CP))
    
    !WRITE(*,*)"K1=",K1

!   GMRES MUDULE

    DO I=1,NPS(CP)
      !CI(THR1+I,THR1+I)=1.0 !/A(I,I)
      !K=SUM(NPS(1:CP))-NPS(CP)+I       !SERIAL NUMBER
      !CI(THR1+I,THR1+I)=1.0/MAL(THR1+I,THR1+I)   !预处理
      CI(THR1+I,THR1+I)=1.0
    END DO

!   加预设条件   
    
    DO I=1,NPS(CP)
      DO J=1,N
        APG(THR1+I,J)=0.
        DO K=1,N
            APG(THR1+I,J)=APG(THR1+I,J)+CI(THR1+I,K)*MAL(K,J)
        END DO
      END DO
    END DO
    
    
    AG(THR1+1:THR2,:)=APG(THR1+1:THR2,:)

!$OMP BARRIER

    DO I=1,NPS(CP)
      APG(1,THR1+I)=0.
      DO J=1,N
        APG(1,THR1+I)=APG(1,THR1+I)+CI(THR1+I,J)*VERR(J)
      END DO
    END DO

    BB(THR1+1:THR2)=APG(1,THR1+1:THR2)

    BNORM1(CP)=0
    DO I=1,NPS(CP)
      BNORM1(CP)=BNORM1(CP)+BB(THR1+I)*BB(THR1+I)
    END DO
!   向量的模 
    !CALL FNORM(BB(K1:K2),BB(K1:K2),N,1,BNORM1(CP))

!$OMP END PARALLEL

	BNORM=SQRT(SUM(BNORM1(1:NTHREAD)))

    IF(NIT.EQ.1) THEN
      QG=0
    END IF
    !QG=0
    !RG=BB

    RG=0

    DO I=1,N
      DO J=1,N
        !RG(I)=RG(I)+BB(I)-AG(I,J)*QG(J)
        RG(I)=RG(I)-AG(I,J)*QG(J)
      END DO
      RG(I)=RG(I)+BB(I)
    END DO
   
    RNORM1=0
    DO I=1,N
      RNORM1=RNORM1+RG(I)*RG(I)    !改动前
      !RNORM1=RNORM1+ABS(RG(I))
    END DO
    
    IF(SQRT(RNORM1).LE.TOL) GOTO 109    !改动前
    !IF(RNORM1/N.LE.TOL) GOTO 109    

 88   CONTINUE
    
    RNORM1=0
    DO I=1,N
      RNORM1=RNORM1+RG(I)*RG(I)
    END DO
      
!   SET UP V^1 AND G^0

    !WRITE(*,*)'K1',K1,NPS(CP) 

	
    !CALL FNORM(R,R,N,NS,RNORM1)
!C	WRITE(6,*) RNORM1
	RNORM=SQRT(RNORM1)
	
    DO I=1,N
      PG(I)=RG(I)/RNORM
	END DO

    GG=0

    GG(1)=RNORM
!   BEGIN ITERATION
    ITER=1


!	BWRITE(66,*) 'BETA',G(1)
    DO WHILE((ITER.LE.BM).AND.(TOL1.GT.TOL).AND.(ITEROUT*BM+ITER.LE.ITERMAX))
    !DO WHILE((ITER.LE.BM).AND.(TOL1.GT.TOL))

!$OMP PARALLEL PRIVATE(CP,THR1,THR2,I,J)
    
    CP=OMP_GET_THREAD_NUM()+1

    THR1=SUM(NPS(1:CP))-NPS(CP) 
    THR2=SUM(NPS(1:CP))

    !WRITE(*,*)"ITER",ITER

      !DO CP=1,NCUP   !CPU NUMBER

!   SAVE  P AS THE V(ITER)
        !K1=SUM(NPS(1:CP))-NPS(CP)
        !K2=SUM(NPS(1:CP))

        BV(THR1+1:THR2,ITER)=PG(THR1+1:THR2)
        !CALL EQUAL(BV(:,1,ITER),P,N,NS)
 !  WRITE(66,*)'ITER',ITER
!	WRITE(66,*) ((BV(I1,J1,ITER),I1=1,N),J1=1,NS)

!   FORM  AV(ITER)	 
        APG(1,THR1+1:THR2)=0
        DO I=1,NPS(CP)
          DO J=1,N
            APG(1,THR1+I)=APG(1,THR1+I)+AG(THR1+I,J)*PG(J)
          END DO
        END DO

!$OMP BARRIER 

	  !CALL MATAV(A,P,AP,N,NS)  !A*P
!	IF(ITER.EQ.1) WRITE(*,*) 'AP', AP(3,1)
!   INITIALIZE V^(ITER+1) WITH AV(ITER)
        PG(THR1+1:THR2)=APG(1,THR1+1:THR2)
        !WI(K1+1:K2,1)=APG(1,K1+1:K2)
     !GOTO 11
! MAKE V^(ITER+1) ORTHOGONAL TO V^(I),I<=ITER
 
!!$OMP BARRIER 
     
        DO  J=1,ITER
! 1)正交化
	   !CALL FNORM(BV(1,1,J),AP,N,NS,HIJ)
           
!$OMP SINGLE  
            HGIJ=0    
            DO I=1,N
              HGIJ=HGIJ+BV(I,J)*APG(1,I)
              !HGIJ=HGIJ+BV(1,J)*APG(1,I)
            END DO
            BH(J,ITER)=HGIJ
!$OMP END SINGLE 
          !TPR(CP)=0
          !DO I=1,NPS(CP)
            !TPR(CP)=TPR(CP)+BV(K1+I,J)*PG(1,K1+I)
          !END DO      
         !DO I=1,N
           !WI(I,NS,J1+1)=WI(I,NS,J1)-BV(I,NS,J)*TPR
         !END DO
         !WRITE(*,*)SQRT(RNORM/RNORM1)

          DO  I=1,NPS(CP)
		    !P(I,K)=P(I,K)-TPR*BV(I,K,J)
            PG(THR1+I)=PG(THR1+I)-HGIJ*BV(THR1+I,J)
	      ENDDO 
	      !BH(J,ITER)=TPR

          !BH(J,ITER)=HGIJ


!$OMP BARRIER

          !END DO 
          !CALL FNORM(WI(:,:,J1),WI(:,:,J1),N,NS,RNORM)
          !CALL FNORM(WI(:,:,J1+1),WI(:,:,J1+1),N,NS,RNORM1)
          !IF(SQRT(RNORM/RNORM1).LE.TH) EXIT         
	    END DO 
        !WRITE(*,*)"BH(J,ITER)",BH(J,ITER)

!$OMP END PARALLEL 


       !WRITE(*,*)"PG",PG(1:N)

       RNORM1=0
       DO I=1,N
         RNORM1=RNORM1+PG(I)*PG(I)
       END DO     

       !PG(K1+1:K2)=APG(1,K1+1:K2)

!NORMALIZE V^(ITER+1)
! 2)标准化
	 !CALL FNORM(P,P,N,NS,RNORM1)
	 !RNORM =REAL(RNORM1)
     !RNORM=SQRT(RNORM1)

       DO  I=1,N
		  PG(I)=PG(I)/SQRT(RNORM1)
	   END DO

     BH(ITER+1,ITER)=SQRT(RNORM1) 
!	 WRITE(66,*)'H OLD BEFORE ', (BH(I,ITER),I=1,ITER+1)    

! APPLY ROTATIONS TO NEW H COLUMN

! 3)更新V和H
     DO I=1,ITER-1
	   HGIJ=BH(I,ITER)
	   HGIJ1=BH(I+1,ITER)
	   BH(I,ITER)=HGIJ*CG(I)+HGIJ1*SG(I)
	   BH(I+1,ITER)=-HGIJ*SG(I)+HGIJ1*CG(I)
	 ENDDO
!	 WRITE(66,*)'H OLD AFTER ', (BH(I,ITER),I=1,ITER)  
!    COMPUTE NEW ROTATIONS
     HGIJ=BH(ITER,ITER)
	 HGIJ1=BH(ITER+1,ITER) 	   
     
     CALL GIVENS(HGIJ,HGIJ1,CG(ITER),SG(ITER))
!    WRITE(66,*) 'GIVENS',C(ITER),S(ITER)
!    APPLY NEW ROTATIONS
     BH(ITER,ITER)= HGIJ*CG(ITER)+HGIJ1*SG(ITER)
     BH(ITER+1,ITER)=0.
!    WRITE(66,*)'H NEW AFTER ', (BH(I,ITER),I=1,ITER+1) 

     HGIJ=GG(ITER)
     GG(ITER)=CG(ITER)*HGIJ
	 GG(ITER+1)=-SG(ITER)*HGIJ
     TOL1=ABS(GG(ITER+1))/BNORM
!    WRITE(66,*)'G  AFTER ', (G(I ),I=1,ITER+1) 
!C   TOL1=CABS(G(ITER+1))
!    WRITE(6,*) 'ITER',ITER
!	 WRITE(6,*) 'TOL1',TOL1
!	 WRITE(66,*) 'TOLERANCE',TOL1
!    WRITE(151,*) ITEROUT*M+ITER,LOG10(TOL1)
	 ITER=ITER+1
     
     IF(MOD(ITER,10).EQ.0) WRITE(*,*)ITER,TOL1
        !WRITE(211,*)ITER,TOL1

    END DO


    ITER=ITER-1
!C  WRITE(6,*) 'TOL1',TOL1
 
!   COMPUTE  LEAST SQUARE  SOLUTION
    DO I=1,ITER
	  YG(I)=GG(I)
	ENDDO
!	WRITE(66,*)'G'
!	WRITE(66,*) (G(I),I=1,ITER)
	DO I=ITER ,1,-1 
 	   YG(I)=YG(I)/BH(I,I)   
 	   DO J=I-1,1,-1
          YG(J)=YG(J)-BH(J,I)*YG(I)
       ENDDO
    ENDDO
!	WRITE(66,*) 'Y'
!	WRITE(66,*)(Y(J),J=1,ITER)
!	 WRITE(66,*)'BV'
	
    !DO  K=1,ITER
!	 WRITE(66,*)'ITER',K
	
!       WRITE(66,*)(BV(I,1,K),I=1,N)
	!ENDDO

!    ACQUIRE ITERATION SOLUTION
	 DO I=1,N
	   DO K=1,ITER
         QG(I)=QG(I)+YG(K)*BV(I,K)
	   ENDDO
	 ENDDO
!	 WRITE(66,*) 'Q'
!    WRITE(66,*) (Q(I,1),I=1,N)

     IF((TOL1.LT.TOL).AND.(ITEROUT*BM+ITER.LE.ITERMAX)) THEN
     !IF((TOL1.LT.TOL)) THEN
!      ITERATIONS END
	 GOTO 98
	 ENDIF
     
     KSTART=0
	 IF((TOL1.GT.TOL).AND.(ITEROUT*BM+ITER.LT.ITERMAX)) KSTART=1
     !IF((TOL1.GT.TOL)) KSTART=1

     IF( (ITEROUT*BM+ITER.GT.ITERMAX).AND.(TOL1.GT.TOL) ) THEN
     
     !IF( (TOL1.GT.TOL) ) THEN
!	   WRITE(6,*) 'WARNING NO CONVERGENCE  AFTER ' 
!	   WRITE(6,*) ITEROUT*M+ITER,'ITERATIONS'
!	   WRITE(6,*) 'TOLRENCE  IS',TOL1
!C	   WRITE(*,*) 'ITEROUT',ITEROUT
       WRITE(*,*) 'TOLRENCE  IS',TOL1
       GOTO  109
	 ENDIF

    !WRITE(*,*) 'TOLRENCE  IS',TOL1

	IF(KSTART.EQ.1) THEN

	ITEROUT=ITEROUT+1
!C  　    
!       WRITE(6,*) 'ITEROUT',ITEROUT     

    APG=0
    DO I=1,N
      DO J=1,N
        APG(1,I)=APG(1,I)+AG(I,J)*QG(J)
      END DO
    END DO
                                                                                                                                                                                                                                                          
    !CALL  MATAV(A,Q,AP,N,NS)
    DO  I=1,N
	  RG(I)=BB(I)-APG(1,I)
	END DO
!     ITERATIONS BEGIN AGAIN
      
	GOTO 88
    ENDIF
  98  CONTINUE
!      WRITE(6,*) 'CONGRATULATIONS CONVERGENCE ACQUIRE ' 
!	   WRITE(6,*)    'AFTER', ITEROUT*M+ITER,' ITERATIONS'
	  !WRITE(*,*) 'TOLRENCE  IS',TOL1  

  99  CONTINUE 
    BB=QG

	ITR=ITEROUT*BM+ITER
	!CLOSE(66)


  109  CONTINUE
  
  WRITE(*,*)"ITR=   ",ITR,TOL1     
        
  DO I=1,M+L !+L
     IF(I.LE.NBPOINT) THEN
        !PHIS(I)=QG(I)
        !WRITE(*,*)X(I)
     ELSE
        !PHIS_N(I)=QG(I)
     END IF
     WRITE(1,*)QG(I)
     !WRITE(*,*)I,QG(I)
  END DO
  !STOP

  !WRITE(*,*)QG(30)

!   END GMRES MUDULE

14	CLOSE(16)   
    CLOSE(1)
!C 释放内存
  DEALLOCATE(CI,APG,AG,WI,BV,BH)
  DEALLOCATE(BB,RG,PG,BNORM1,GG,YG,CG,SG,QG)
  DEALLOCATE(MAL,VERR)
END SUBROUTINE

SUBROUTINE GIVENS(Z1,Z2,C,S)
!
       REAL  Z1,Z2,S
       REAL  C
!
       REAL  NORMZ
       INTRINSIC ABS , SQRT 
!
       REAL  ZERO, ONE
       PARAMETER (ZERO = 0.0E0, ONE = 1.0E0 )
       REAL  DZRO,DONE
       PARAMETER (DZRO = 0.0E0 ,DONE = 1.0E0 ) 
!
       NORMZ =SQRT(Z1**2+Z2**2)
       IF (ABS (Z1).NE.DZRO) THEN
         S=Z2*ABS(Z1)/(Z1*NORMZ)
         C=ABS(Z1)/NORMZ
       ELSE IF (ABS(Z2).NE.DZRO) THEN
         C=DZRO
         S=Z2/ABS(Z2)
       ELSE
         C=DONE
         S=ZERO
       ENDIF
!
       RETURN
END SUBROUTINE   