SUBROUTINE LU_MKL
    USE GREENMOD
    
    INTEGER,DIMENSION(:),ALLOCATABLE:: IPIV
    INTEGER:: INFO1,INFO2
    REAL*8,DIMENSION(:,:),ALLOCATABLE:: MATA
    REAL*8,DIMENSION(:),ALLOCATABLE:: MATB,MATX
    
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2
    CHARACTER TRANS
   
    ALLOCATE(IPIV(NTPN))
    ALLOCATE(MATA(NTPN,NTPN),MATB(NTPN),MATX(NTPN))
    
    DO I=1,NTPN
        MATA(I,1:NTPN)=MAL(I,1:NTPN)
    END DO
    MATB(1:NTPN)=VERR(1:NTPN) 
    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2) 
    !PRINT *, 'BEGIN TIME  LU_MKL: ', CHAR_TIME2,'  ON  ',CHAR_DATE2     
    
    TRANS='N'
    call DGETRF(NTPN,NTPN,MATA,NTPN,IPIV,INFO1)
    call DGETRS(TRANS, NTPN,1, MATA, NTPN, IPIV, MATB, NTPN, INFO2)
    
    QG(1:NTPN)=MATB(1:NTPN)
    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)     
    !PRINT *, 'END   TIME  LU_MKL: ', CHAR_TIME2,'  ON  ',CHAR_DATE2 
    
    DEALLOCATE(MATA,MATB,MATX)
    DEALLOCATE(IPIV)
END SUBROUTINE
    
     
    
SUBROUTINE BLLU
	USE GREENMOD

    INTEGER:: N,K,NNN
    REAL,ALLOCATABLE:: ML(:,:),MU(:,:)
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2

    !OPEN(991,FILE='LU_FMATRIX.DAT')
    OPEN(1,FILE='PHI.DAT')
    
    NNN=NTPN !+NFPOINT
    ALLOCATE(ML(NNN,NNN),MU(NNN,NNN),RG(NNN))
    ML=0.
    MU=0.
    RG=0.


    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)     
    PRINT *, 'BEGIN   TIME  BLLU: ', CHAR_TIME2,'  ON  ',CHAR_DATE2     
    
	DO 100 K=1,NNN-1
	  !IF (ABS(MAL(K,K))+1.0.EQ.1.0) THEN
	  !  WRITE(*,*) '***FAIL***'
	  !  JT=0
	  !END IF
!	  DO 10 I=K+1,NNN
      MAL(I,K)=MAL(I,K)/MAL(K,K)
10	  MAL(I,K)=MAL(I,K)/MAL(K,K)
      

!      DO 20 I=K+1,NNN
!      DO 20 J=K+1,NNN
!20	   MAL(I,J)=MAL(I,J)-MAL(I,K)*MAL(K,J)
      
      CALL SGEMM('N','N',NNN-K,NNN-K,1,-1.,MAL(K+1:NNN,K),NNN-K,MAL(K,K+1:NNN),1,1.,MAL(K+1:NNN,K+1:NNN),NNN-K)
      !CALL SAXPY(NNN-K,MAL(K+1:NNN,K+1:NNN),MAL(K+1:NNN,K),1,MAL(I,K+1:NNN),1)

      !CALL VSMUL(NNN-K,MAL(K,K+1:NNN),ML(I,K),TPM(K+1:NNN))
      !MAL(I,K+1:NNN)=MAL(I,K+1:NNN)-TPM(K+1:NNN)
      
      !CALL SAXPY(1,MAL(I,K+1:NNN),ML(I,K),1,MAL(I,K+1:NNN),1)   
      !CALL SGEMV('N',NNN-K,1,1,MAL(I,K+1:NNN),NNN-K,ML(I,K),1,1,MAL(I,K+1:NNN),1)

!10    CONTINUE
!     DO 20 I=K+1,NNN
!	  DO 20 J=K+1,NNN
!20	  MAL(I,J)=MAL(I,J)-MAL(I,K)*MAL(K,J)
100	CONTINUE


        !IF(((K/1000)*1000).EQ.K)THEN
		!	WRITE(*,101)K,NNN 
		!END IF
        !101	FORMAT('   LU',3X,'I=',I8,5X,I8)
200	CONTINUE

    !DO I=1,NNN
    !    DO J=1,I-1
    !        ML(I,J)=MAL(I,J)
    !        MU(I,J)=0.0
    !    END DO
        
    !    ML(I,I)=1.0
    !    MU(I,I)=MAL(I,I)
    !    DO J=I+1,NNN
    !        ML(I,J)=0.0
    !        MU(I,J)=MAL(I,J)
    !    END DO
        !WRITE(991,*)MAL(I,I),ML(I,I),MU(I,I)
    !END DO     
    
    
    CLOSE(991)
    JT=1

    DO I=2,NNN
        DO J=1,I-1
            VERR(I)=VERR(I)-MAL(I,J)*VERR(J)
        END DO
    END DO

    QG=0
    QG(NNN)=VERR(NNN)/MAL(NNN,NNN)
    DO I=NNN-1,1,-1
        QG(I)=VERR(I) !/U(I,I)
        DO J=I+1,NNN
            QG(I)=QG(I)-MAL(I,J)*QG(J) !/U(I,I)
        END DO
        QG(I)=QG(I)/MAL(I,I)
    END DO

    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)     
    PRINT *, 'END   TIME  BLLU: ', CHAR_TIME2,'  ON  ',CHAR_DATE2     
    
    DO I=1,NNN
		!WRITE(1,*)RG(I)
    END DO
    CLOSE(1)
    DEALLOCATE(ML,MU,RG)
	DEALLOCATE(MAL,VERR)
END SUBROUTINE    

SUBROUTINE AGSDL(A,B,N,EPS)
	DIMENSION A(N,N),B(N),X(N)
	DOUBLE PRECISION A,B,X,T,S,P,Q

    OPEN(1,FILE='PHI.DAT')

	DO 5 I=1,N
	  IF (ABS(A(I,I))+1.0.EQ.1.0) THEN
	    L=0
	    WRITE(*,100)
        !WRITE(*,*)A(I,I)
	    RETURN
	  END IF
5	CONTINUE
100	FORMAT(1X,'  FAIL')
	L=100
	DO 10 I=1,N
10	X(I)=0.0
20	P=0.0
	L=L-1
	DO 50 I=1,N
	  T=X(I)
	  S=0.0
	  DO 30 J=1,N
	    IF (J.NE.I) S=S+A(I,J)*X(J)
30	  CONTINUE
	  X(I)=(B(I)-S)/A(I,I)
	  Q=ABS(X(I)-T)/(1+ABS(X(I)))
	  IF (Q.GT.P) P=Q
50	CONTINUE
	IF ((P.GE.EPS).AND.(L.NE.0)) GOTO 20
	IF (L.EQ.0) WRITE(*,100)

    DO I=1,N
		WRITE(1,*)X(I)
    END DO
    CLOSE(1)

	!DEALLOCATE(MAL,VERR)
	RETURN
END
