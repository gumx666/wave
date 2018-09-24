SUBROUTINE PRE_CONDITION_2
    USE GREENMOD
    USE CUMOD

    !IMPLICIT NONE

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

    REAL:: ERR=1E-2,TAU=1.E-3
    
    !IF(FR.LE.0.15) THEN
    !    ERR=1.E-6
    !    TAU=1.E-5
    !END IF

    !NBPOINT=3
    !NFPOINT=3

    M=NBPOINT
    L=NFPOINT
    N=M+L
  



    IF(ITTE.EQ.1.AND.NIT.EQ.1) THEN
    ALLOCATE(LAGA(NTPN*NTPN/2),UAA(NTPN*NTPN/2),NORMJ(N+500),ROW(N+500),&
    IALA(NTPN+1),IAUA(NTPN+1),JALA(NTPN*NTPN/2),JAUA(NTPN*NTPN/2),UA(NTPN,NTPN))
    END IF
    
    
    !DAG=0.
    !DDA=0.
    !LAG=0.
    !UA=0.
    !DA_1=0.
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

    COUT=0


    T1=0.
    T2=0.

    
    !UA(1,:)=MAL(1,:)
    UA(1,:)=MAL(1,:)
    ROW=0.
    LAGA(1)=1
    UAA(1:N)=MAL(1,1:N)
    IALA(1)=1
    IAUA(1)=1
    JALA(1)=1
    DO J=1,NTPN
    JAUA(J)=J
    END DO
    ILA=1
    IUA=NTPN
    JLA=1
    JUA=NTPN    
    DO I=2,N
		!IF((((I)/2000)*2000).EQ.I)THEN
	    !		WRITE(*,200)I,N !NTPN
		!END IF
    !200	FORMAT('ILU',3X,'6HINODE=',I8,5X,I8)
        ROW(:)=MAL(I,:)

        
        IALA(I)=JLA+1
        IAUA(I)=JUA+1
        DO K=1,I-1
            IF(ABS(ROW(K)).GT.ERR) THEN
                ROW(K)=ROW(K)/UA(K,K)
                !WRITE(*,*)"MAL(I,I)  ",UA(K,K),MAL(K,K)
                !ROW(K)=ROW(K)/MAL(K,K)
                
                IF(ABS(ROW(K)).LT.TAU*NORMJ(I)) THEN     !DROPPING RULE OF TAO
                    ROW(K)=0.
                ELSE
                    !ROW(K+1:N)=ROW(K+1:N)-ROW(K)*UA(K,K+1:N)
                    DO K1=IAUA(K)+1,IAUA(K+1)-1
         
                        ROW(JAUA(K1))=ROW(JAUA(K1))-ROW(K)*UA(K,JAUA(K1))
                        !ROW(K+1:N)
                    END DO
                END IF
                
            END IF
        END DO

        
                JUA=JUA+1
                UAA(JUA)=ROW(I)
                JAUA(JUA)=I
                UA(I,I)=ROW(I)
        DO K=1,N
            !IF(ABS(ROW(K)).LT.TAU*NORMJ(I).AND.K.NE.I&
            !   .AND.LEV(I,K).GT.LEVP) ROW(K)=0.
            
            !IF(ABS(ROW(K)).LT.TAU*NORMJ(I).AND.K.NE.I.OR.ABS(ROW(K)).LT.ERR) THEN
            !IF(ABS(ROW(K)).LT.TAU*NORMJ(I).AND.ABS(ROW(K)).LT.ERR.AND.I.NE.K) THEN
            
            IF(ABS(ROW(K)).LE.ERR) THEN
                !ROW(K)=0.
            ELSE IF(K.LE.I-1) THEN
                !LAG(I,K)=ROW(K)
                JLA=JLA+1
                LAGA(JLA)=ROW(K)
                JALA(JLA)=K
                !LAG(I,K)=ROW(K)
            ELSE IF(K.GT.I) THEN
                JUA=JUA+1
                UAA(JUA)=ROW(K)
                JAUA(JUA)=K
                UA(I,K)=ROW(K)
            END IF   
        END DO
                JLA=JLA+1
                LAGA(JLA)=1
                JALA(JLA)=I
              
        
        !LAG(I,1:I-1)=ROW(1:I-1)
        !UA(I,I:N)=ROW(I:N)
        ROW(1:N)=0.
    END DO
    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)     
    PRINT *, 'END   TIME  ilut: ', CHAR_TIME2,'  ON  ',CHAR_DATE2    
    
        IALA(N+1)=JLA+1
        IAUA(N+1)=JUA+1

    !WRITE(*,*)T1,T2
    DO I=1,NTPN
        DO J=1,NTPN
            !LAGA((I-1)*NTPN+J)=LAG(I,J)
            !JALA((I-1)*NTPN+J)=J
            !UAA((I-1)*NTPN+J)=UA(I,J)
            !JAUA((I-1)*NTPN+J)=J
        END DO
        !IALA(I)=(I-1)*NTPN+1
        !IAUA(I)=(I-1)*NTPN+1
        !WRITE(*,*)IA(I)
    END DO
    
    !IALA(NTPN+1)=NTPN*NTPN+1
    !IAUA(NTPN+1)=NTPN*NTPN+1

 

END SUBROUTINE


SUBROUTINE FGMRES_MKL
    
    USE GREENMOD

INCLUDE "mkl_rci.fi"    



      INTEGER SIZE,IERR
      PARAMETER (SIZE=128)
      INTEGER IPAR(SIZE)
      DOUBLE PRECISION DPAR(SIZE), TMP(NTPN*(2*150+1)+(150*(150+9))/2+1)
      DOUBLE PRECISION RESIDUAL(NTPN)
      DOUBLE PRECISION RHS(NTPN), B(NTPN)
      DOUBLE PRECISION TRVEC(NTPN)
      DOUBLE PRECISION TOL
!---------------------------------------------------------------------------
! Some additional variables to use with the RCI (P)FGMRES solver
!---------------------------------------------------------------------------
      INTEGER ITERCOUNT, EXPECTED_ITERCOUNT
      PARAMETER (EXPECTED_ITERCOUNT=4)
      INTEGER RCI_REQUEST, I
      DOUBLE PRECISION DVAR
      INTEGER MATSIZE, INCX, REF_NIT
      INTEGER MAXFIL    
!---------------------------------------------------------------------------
! An external BLAS function is taken from MKL BLAS to use
! with the RCI (P)FGMRES solver
!---------------------------------------------------------------------------
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2
    
    
    
    INTEGER,DIMENSION(:),ALLOCATABLE:: IPIV,IA,JA
    INTEGER:: INFO1,INFO2
    REAL*8,DIMENSION(:),ALLOCATABLE:: A
    REAL*8,DIMENSION(:),ALLOCATABLE:: MATB,MATX
    
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2
    CHARACTER TRANS
   
    !PAUSE
    ALLOCATE(IPIV(NTPN))
    ALLOCATE(A(NTPN*NTPN),IA(NTPN+1),JA(NTPN*NTPN),MATB(NTPN),MATX(NTPN))
    !PAUSE
    
    DO I=1,NTPN
        DO J=1,NTPN
            A((I-1)*NTPN+J)=MAL(I,J)
            JA((I-1)*NTPN+J)=J
        END DO
        IA(I)=(I-1)*NTPN+1
        !WRITE(*,*)IA(I)
    END DO
    IA(NTPN+1)=NTPN*NTPN+1
    MATB(1:NTPN)=VERR(1:NTPN) 

    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2) 
    PRINT *, 'BEGIN TIME  GMRES_MKL: ', CHAR_TIME2,'  ON  ',CHAR_DATE2     
    
    TRANS='N'
    
!*******************************************************************************
!  Content:
!  Intel MKL RCI (P)FGMRES ((Preconditioned) Flexible Generalized Minimal
!                                                       RESidual method) example
!*******************************************************************************

!---------------------------------------------------------------------------
!  Example program for solving non-symmetric indefinite system of equations
!  Fully advanced case: full functionality of RCI FGMRES solver is exploited
!---------------------------------------------------------------------------


      DO I=1,NTPN
        MATX(I)=0.
      ENDDO

      PRINT *,'--------------------------------------------------------'
      PRINT *,'The FULLY ADVANCED example of usage of RCI FGMRES solver'
      PRINT *,'   to solve a non-symmetric indefinite non-degenerate'
      PRINT *,'          algebraic system of linear equations'
      PRINT *,'--------------------------------------------------------'
!---------------------------------------------------------------------------
! Initialize variables and the right hand side through matrix-vector product
!---------------------------------------------------------------------------
      
    !CALL MKL_DCSRGEMV('N', NTPN, A, IA, JA, MATX, MATB)

!---------------------------------------------------------------------------
! Save the right-hand side in vector B for future use
!---------------------------------------------------------------------------
	 !CALL DCOPY(N, RHS, 1, B, 1)
!---------------------------------------------------------------------------
! Initialize the initial guess
!---------------------------------------------------------------------------

      B(1:NTPN)=VERR(1:NTPN) 
      RHS(1:NTPN)=VERR(1:NTPN)
          

!---------------------------------------------------------------------------
! Initialize the solver
!---------------------------------------------------------------------------
      CALL DFGMRES_INIT(NTPN, MATX, RHS, RCI_REQUEST, IPAR, DPAR, TMP)
      IF (RCI_REQUEST.NE.0) GOTO 999

      IPAR(31)=1
	  DPAR(31)=1.D-5
	  TOL=1.d-6
	  MAXFIL=1     
      MAXFIL=NTPN/5
      !CALL DCSRILUT(NTPN, A, IA, JA, BILUT, IBILUT, JBILUT,TOL, MAXFIL, IPAR, DPAR, IERR)
      !WRITE(*,*)IEER       
      
      CALL PRE_CONDITION_2
      CALL TIME(CHAR_TIME2)
	  CALL DATE(CHAR_DATE2) 
      PRINT *, 'END TIME  ILUT: ', CHAR_TIME2,'  ON  ',CHAR_DATE2    
      
!---------------------------------------------------------------------------
! Set the desired parameters:
! do the restart after 2 iterations
! LOGICAL parameters:
! do not do the stopping test for the maximal number of iterations
! do the Preconditioned iterations of FGMRES method
! DOUBLE PRECISION parameters
! set the relative tolerance to 1.0D-3 instead of default value 1.0D-6
!---------------------------------------------------------------------------
      IPAR(15)=15
      IPAR(8)=0
      IPAR(11)=1
      DPAR(1)=1.0D-7      
!---------------------------------------------------------------------------
! Check the correctness and consistency of the newly set parameters
!---------------------------------------------------------------------------
      CALL DFGMRES_CHECK(NTPN, MATX, RHS, RCI_REQUEST, IPAR, DPAR, TMP)
      !WRITE(*,*)"RCI_REQUEST      ",RCI_REQUEST

      IF (RCI_REQUEST.NE.0) GOTO 999
!---------------------------------------------------------------------------
! Print the info about the RCI FGMRES method
!---------------------------------------------------------------------------
      PRINT *, ''
      PRINT *,'Some info about the current run of RCI FGMRES method:'
      PRINT *, ''
      IF (IPAR(8).NE.0) THEN
         WRITE(*,'(A,I1,A,A)') 'As IPAR(8)=',IPAR(8),', the automatic',' test for the maximal number of iterations will be'
         PRINT *,'performed'
      ELSE
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(8)=',IPAR(8),', the automatic',' test for the maximal number of iterations will be'
      	PRINT *,'skipped'
      ENDIF
      PRINT *,'+++'
      IF (IPAR(9).NE.0) THEN
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(9)=',IPAR(9),', the automatic',' residual test will be performed'
      ELSE
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(9)=',IPAR(9),', the automatic',' residual test will be skipped'
      ENDIF
      PRINT *,'+++'
      IF (IPAR(10).NE.0) THEN
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(10)=',IPAR(10),', the',' user-defined stopping test will be requested via'
      	PRINT *,'RCI_REQUEST=2'
      ELSE
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(10)=',IPAR(10),', the',' user-defined stopping test will not be requested, thus,'
      	PRINT *,'RCI_REQUEST will not take the value 2'
      ENDIF
      PRINT *,'+++'
      IF (IPAR(11).NE.0) THEN
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(11)=',IPAR(11),', the',' Preconditioned FGMRES iterations will be performed, thus,'
      	WRITE(*,'(A,A)') 'the preconditioner action will be requested', ' via RCI_REQUEST=3'
      ELSE
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(11)=',IPAR(11),', the',' Preconditioned FGMRES iterations will not be performed,'
      	WRITE(*,'(A)') 'thus, RCI_REQUEST will not take the value 3'
      ENDIF
      PRINT *,'+++'
      IF (IPAR(12).NE.0) THEN
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(12)=',IPAR(12),', the automatic',' test for the norm of the next generated vector is'
      	WRITE(*,'(A,A)') 'not equal to zero up to rounding and',' computational errors will be performed,'
      	PRINT *,'thus, RCI_REQUEST will not take the value 4'
      ELSE
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(12)=',IPAR(12),', the automatic',' test for the norm of the next generated vector is'
      	WRITE(*,'(A,A)') 'not equal to zero up to rounding and',' computational errors will be skipped,'
      	WRITE(*,'(A,A)') 'thus, the user-defined test will be requested', ' via RCI_REQUEST=4'
      ENDIF
      PRINT *,'+++'
!---------------------------------------------------------------------------
! Compute the solution by RCI (P)FGMRES solver with preconditioning
! Reverse Communication starts here
!---------------------------------------------------------------------------
1     CALL DFGMRES(NTPN, MATX, RHS, RCI_REQUEST, IPAR, DPAR, TMP)
      !WRITE(*,*)"RCI_REQUEST      ",RCI_REQUEST
!---------------------------------------------------------------------------
! If RCI_REQUEST=0, then the solution was found with the required precision
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.0) GOTO 3
!---------------------------------------------------------------------------
! If RCI_REQUEST=1, then compute the vector A*TMP(IPAR(22))
! and put the result in vector TMP(IPAR(23))
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.1) THEN
      	CALL MKL_DCSRGEMV('N',NTPN, A, IA, JA, TMP(IPAR(22)), TMP(IPAR(23)))
      	GOTO 1
      ENDIF      
!---------------------------------------------------------------------------
! If RCI_request=2, then do the user-defined stopping test
! The residual stopping test for the computed solution is performed here
!---------------------------------------------------------------------------
! NOTE: from this point vector B(N) is no longer containing the right-hand
! side of the problem! It contains the current FGMRES approximation to the
! solution. If you need to keep the right-hand side, save it in some other
! vector before the call to DFGMRES routine. Here we saved it in vector
! RHS(N). The vector B is used instead of RHS to preserve the original
! right-hand side of the problem and guarantee the proper restart of FGMRES
! method. Vector B will be altered when computing the residual stopping
! criterion!
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.2) THEN
! Request to the DFGMRES_GET routine to put the solution into B(N) via IPAR(13)
      	IPAR(13)=1
! Get the current FGMRES solution in the vector B(N)
      	CALL DFGMRES_GET(NTPN, MATX, B, RCI_REQUEST, IPAR, DPAR, TMP, ITERCOUNT)
! Compute the current true residual via MKL (Sparse) BLAS routines
      	CALL MKL_DCSRGEMV('N', NTPN, A, IA, JA, B, RESIDUAL)
      	CALL DAXPY(NTPN, -1.0D0, RHS, 1, RESIDUAL, 1)
      	DVAR=DNRM2(NTPN, RESIDUAL, 1)
        !WRITE(*,*)"DVAR  ",DVAR
      	IF (DVAR.LT.1.0E-3) THEN
      	   GOTO 3
      	ELSE
      	   GOTO 1
      	ENDIF
      ENDIF         
!---------------------------------------------------------------------------
! If RCI_REQUEST=3, then apply the preconditioner on the vector
! TMP(IPAR(22)) and put the result in vector TMP(IPAR(23))
! Here is the recommended usage of the result produced by ILUT routine
! via standard MKL Sparse Blas solver routine mkl_dcsrtrsv.
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.3) THEN
       !CALL MKL_DCSRTRSV('L','N','U',NTPN,BILUT,IBILUT,JBILUT,TMP(IPAR(22)),TRVEC)
       !CALL MKL_DCSRTRSV('U','N','N',NTPN,BILUT,IBILUT,JBILUT,TRVEC,TMP(IPAR(23)))
       CALL MKL_DCSRTRSV('L','N','U',NTPN,LAGA,IALA,JALA,TMP(IPAR(22)),TRVEC)
       CALL MKL_DCSRTRSV('U','N','N',NTPN,UAA,IAUA,JAUA,TRVEC,TMP(IPAR(23)))
       GOTO 1
      ENDIF   
!---------------------------------------------------------------------------
! If RCI_REQUEST=4, then check if the norm of the next generated vector is
! not zero up to rounding and computational errors. The norm is contained
! in DPAR(7) parameter
!---------------------------------------------------------------------------  
      IF (RCI_REQUEST.EQ.4) THEN
      	IF (DPAR(7).LT.1.0D-12) THEN
      	   GOTO 3
      	ELSE
      	   GOTO 1
      	ENDIF
!---------------------------------------------------------------------------
! If RCI_REQUEST=anything else, then DFGMRES subroutine failed
! to compute the solution vector: COMPUTED_SOLUTION(N)
!---------------------------------------------------------------------------
      ELSE
      	GOTO 999
      ENDIF      
!---------------------------------------------------------------------------
! Reverse Communication ends here
! Get the current iteration number and the FGMRES solution. (DO NOT FORGET to
! call DFGMRES_GET routine as computed_solution is still containing
! the initial guess!). Request to DFGMRES_GET to put the solution into
! vector COMPUTED_SOLUTION(N) via IPAR(13)
!---------------------------------------------------------------------------
3     IPAR(13)=0
      CALL DFGMRES_GET(NTPN, MATX, RHS, RCI_REQUEST, IPAR, DPAR, TMP, ITERCOUNT)      
      WRITE(*,*)"ITERCOUNT      ",ITERCOUNT

      GOTO 9999
999   WRITE( *,'(A,I2)') 'The solver has returned the ERROR code ',RCI_REQUEST
9999 CONTINUE      
      !WRITE(*,*)"ITERCOUNT      ",ITERCOUNT
    
      DO I=1,NBPOINT
      !WRITE(*,*)RHS(I),B(I),VERR(I)
      END DO
      !STOP 
    
    QG(1:NTPN)=MATX(1:NTPN)
    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)     
    PRINT *, 'END   TIME  GMRES_MKL: ', CHAR_TIME2,'  ON  ',CHAR_DATE2 
    
    DEALLOCATE(A,MATB,MATX)
    DEALLOCATE(IPIV)
END SUBROUTINE
    
    
SUBROUTINE FGMRES_MKL_1
    
    USE GREENMOD

INCLUDE "mkl_rci.fi"    



      INTEGER SIZE,IERR
      PARAMETER (SIZE=128)
      INTEGER IPAR(SIZE)
      DOUBLE PRECISION DPAR(SIZE), TMP(NTPN*(2*NTPN+1)+(NTPN*(NTPN+9))/2+1)
      DOUBLE PRECISION RESIDUAL(NTPN)
      DOUBLE PRECISION RHS(NTPN), B(NTPN)
      INTEGER JBILUT(NTPN*NTPN),IBILUT(NTPN+1)
      DOUBLE PRECISION BILUT(NTPN*NTPN),TRVEC(NTPN)
      DOUBLE PRECISION TOL
!---------------------------------------------------------------------------
! Some additional variables to use with the RCI (P)FGMRES solver
!---------------------------------------------------------------------------
      INTEGER ITERCOUNT, EXPECTED_ITERCOUNT
      PARAMETER (EXPECTED_ITERCOUNT=4)
      INTEGER RCI_REQUEST, I
      DOUBLE PRECISION DVAR
      INTEGER MATSIZE, INCX, REF_NIT
      INTEGER MAXFIL    
!---------------------------------------------------------------------------
! An external BLAS function is taken from MKL BLAS to use
! with the RCI (P)FGMRES solver
!---------------------------------------------------------------------------
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2
    
    
    
    INTEGER,DIMENSION(:),ALLOCATABLE:: IPIV,IA,JA
    INTEGER:: INFO1,INFO2
    REAL*8,DIMENSION(:),ALLOCATABLE:: A
    REAL*8,DIMENSION(:),ALLOCATABLE:: MATB,MATX
    
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2
    CHARACTER TRANS
   
    ALLOCATE(IPIV(NTPN))
    ALLOCATE(A(NTPN*NTPN),IA(NTPN+1),JA(NTPN*NTPN),MATB(NTPN),MATX(NTPN))

    
    DO I=1,NTPN
        DO J=1,NTPN
            A((I-1)*NTPN+J)=MAL(I,J)
            JA((I-1)*NTPN+J)=J
        END DO
        IA(I)=(I-1)*NTPN+1
        !WRITE(*,*)IA(I)
    END DO
    IA(NTPN+1)=NTPN*NTPN+1
    MATB(1:NTPN)=VERR(1:NTPN) 

    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2) 
    PRINT *, 'BEGIN TIME  GMRES_MKL: ', CHAR_TIME2,'  ON  ',CHAR_DATE2     
    
    TRANS='N'
    
!*******************************************************************************
!  Content:
!  Intel MKL RCI (P)FGMRES ((Preconditioned) Flexible Generalized Minimal
!                                                       RESidual method) example
!*******************************************************************************

!---------------------------------------------------------------------------
!  Example program for solving non-symmetric indefinite system of equations
!  Fully advanced case: full functionality of RCI FGMRES solver is exploited
!---------------------------------------------------------------------------


      DO I=1,NTPN
        MATX(I)=0.
      ENDDO

      PRINT *,'--------------------------------------------------------'
      PRINT *,'The FULLY ADVANCED example of usage of RCI FGMRES solver'
      PRINT *,'   to solve a non-symmetric indefinite non-degenerate'
      PRINT *,'          algebraic system of linear equations'
      PRINT *,'--------------------------------------------------------'
!---------------------------------------------------------------------------
! Initialize variables and the right hand side through matrix-vector product
!---------------------------------------------------------------------------
      
    !CALL MKL_DCSRGEMV('N', NTPN, A, IA, JA, MATX, MATB)

!---------------------------------------------------------------------------
! Save the right-hand side in vector B for future use
!---------------------------------------------------------------------------
	 !CALL DCOPY(N, RHS, 1, B, 1)
!---------------------------------------------------------------------------
! Initialize the initial guess
!---------------------------------------------------------------------------

      B(1:NTPN)=VERR(1:NTPN) 
      RHS(1:NTPN)=VERR(1:NTPN)
          

!---------------------------------------------------------------------------
! Initialize the solver
!---------------------------------------------------------------------------
      CALL DFGMRES_INIT(NTPN, MATX, RHS, RCI_REQUEST, IPAR, DPAR, TMP)
      IF (RCI_REQUEST.NE.0) GOTO 999

      IPAR(31)=1
	  DPAR(31)=1.D-5
	  TOL=1.d-6
	  MAXFIL=1     
      MAXFIL=NTPN/5
      CALL DCSRILUT(NTPN, A, IA, JA, BILUT, IBILUT, JBILUT,TOL, MAXFIL, IPAR, DPAR, IERR)
      !WRITE(*,*)IEER
      CALL TIME(CHAR_TIME2)
	  CALL DATE(CHAR_DATE2) 
      PRINT *, 'END TIME  ILUT: ', CHAR_TIME2,'  ON  ',CHAR_DATE2           
      
!---------------------------------------------------------------------------
! Set the desired parameters:
! do the restart after 2 iterations
! LOGICAL parameters:
! do not do the stopping test for the maximal number of iterations
! do the Preconditioned iterations of FGMRES method
! DOUBLE PRECISION parameters
! set the relative tolerance to 1.0D-3 instead of default value 1.0D-6
!---------------------------------------------------------------------------
      IPAR(15)=10
      IPAR(8)=0
      IPAR(11)=1
      DPAR(1)=1.0D-6      
!---------------------------------------------------------------------------
! Check the correctness and consistency of the newly set parameters
!---------------------------------------------------------------------------
      CALL DFGMRES_CHECK(NTPN, MATX, RHS, RCI_REQUEST, IPAR, DPAR, TMP)
      !WRITE(*,*)"RCI_REQUEST      ",RCI_REQUEST

      IF (RCI_REQUEST.NE.0) GOTO 999
!---------------------------------------------------------------------------
! Print the info about the RCI FGMRES method
!---------------------------------------------------------------------------
      PRINT *, ''
      PRINT *,'Some info about the current run of RCI FGMRES method:'
      PRINT *, ''
      IF (IPAR(8).NE.0) THEN
         WRITE(*,'(A,I1,A,A)') 'As IPAR(8)=',IPAR(8),', the automatic',' test for the maximal number of iterations will be'
         PRINT *,'performed'
      ELSE
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(8)=',IPAR(8),', the automatic',' test for the maximal number of iterations will be'
      	PRINT *,'skipped'
      ENDIF
      PRINT *,'+++'
      IF (IPAR(9).NE.0) THEN
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(9)=',IPAR(9),', the automatic',' residual test will be performed'
      ELSE
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(9)=',IPAR(9),', the automatic',' residual test will be skipped'
      ENDIF
      PRINT *,'+++'
      IF (IPAR(10).NE.0) THEN
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(10)=',IPAR(10),', the',' user-defined stopping test will be requested via'
      	PRINT *,'RCI_REQUEST=2'
      ELSE
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(10)=',IPAR(10),', the',' user-defined stopping test will not be requested, thus,'
      	PRINT *,'RCI_REQUEST will not take the value 2'
      ENDIF
      PRINT *,'+++'
      IF (IPAR(11).NE.0) THEN
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(11)=',IPAR(11),', the',' Preconditioned FGMRES iterations will be performed, thus,'
      	WRITE(*,'(A,A)') 'the preconditioner action will be requested', ' via RCI_REQUEST=3'
      ELSE
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(11)=',IPAR(11),', the',' Preconditioned FGMRES iterations will not be performed,'
      	WRITE(*,'(A)') 'thus, RCI_REQUEST will not take the value 3'
      ENDIF
      PRINT *,'+++'
      IF (IPAR(12).NE.0) THEN
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(12)=',IPAR(12),', the automatic',' test for the norm of the next generated vector is'
      	WRITE(*,'(A,A)') 'not equal to zero up to rounding and',' computational errors will be performed,'
      	PRINT *,'thus, RCI_REQUEST will not take the value 4'
      ELSE
      	WRITE(*,'(A,I1,A,A)') 'As IPAR(12)=',IPAR(12),', the automatic',' test for the norm of the next generated vector is'
      	WRITE(*,'(A,A)') 'not equal to zero up to rounding and',' computational errors will be skipped,'
      	WRITE(*,'(A,A)') 'thus, the user-defined test will be requested', ' via RCI_REQUEST=4'
      ENDIF
      PRINT *,'+++'
!---------------------------------------------------------------------------
! Compute the solution by RCI (P)FGMRES solver with preconditioning
! Reverse Communication starts here
!---------------------------------------------------------------------------
1     CALL DFGMRES(NTPN, MATX, RHS, RCI_REQUEST, IPAR, DPAR, TMP)
      WRITE(*,*)"RCI_REQUEST      ",RCI_REQUEST
!---------------------------------------------------------------------------
! If RCI_REQUEST=0, then the solution was found with the required precision
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.0) GOTO 3
!---------------------------------------------------------------------------
! If RCI_REQUEST=1, then compute the vector A*TMP(IPAR(22))
! and put the result in vector TMP(IPAR(23))
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.1) THEN
      	CALL MKL_DCSRGEMV('N',NTPN, A, IA, JA, TMP(IPAR(22)), TMP(IPAR(23)))
      	GOTO 1
      ENDIF      
!---------------------------------------------------------------------------
! If RCI_request=2, then do the user-defined stopping test
! The residual stopping test for the computed solution is performed here
!---------------------------------------------------------------------------
! NOTE: from this point vector B(N) is no longer containing the right-hand
! side of the problem! It contains the current FGMRES approximation to the
! solution. If you need to keep the right-hand side, save it in some other
! vector before the call to DFGMRES routine. Here we saved it in vector
! RHS(N). The vector B is used instead of RHS to preserve the original
! right-hand side of the problem and guarantee the proper restart of FGMRES
! method. Vector B will be altered when computing the residual stopping
! criterion!
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.2) THEN
! Request to the DFGMRES_GET routine to put the solution into B(N) via IPAR(13)
      	IPAR(13)=1
! Get the current FGMRES solution in the vector B(N)
      	CALL DFGMRES_GET(NTPN, MATX, B, RCI_REQUEST, IPAR, DPAR, TMP, ITERCOUNT)
! Compute the current true residual via MKL (Sparse) BLAS routines
      	CALL MKL_DCSRGEMV('N', NTPN, A, IA, JA, B, RESIDUAL)
      	CALL DAXPY(NTPN, -1.0D0, RHS, 1, RESIDUAL, 1)
      	DVAR=DNRM2(NTPN, RESIDUAL, 1)
        WRITE(*,*)"DVAR  ",DVAR
      	IF (DVAR.LT.1.0E-3) THEN
      	   GOTO 3
      	ELSE
      	   GOTO 1
      	ENDIF
      ENDIF         
!---------------------------------------------------------------------------
! If RCI_REQUEST=3, then apply the preconditioner on the vector
! TMP(IPAR(22)) and put the result in vector TMP(IPAR(23))
! Here is the recommended usage of the result produced by ILUT routine
! via standard MKL Sparse Blas solver routine mkl_dcsrtrsv.
!---------------------------------------------------------------------------
      IF (RCI_REQUEST.EQ.3) THEN
       CALL MKL_DCSRTRSV('L','N','U',NTPN,BILUT,IBILUT,JBILUT,TMP(IPAR(22)),TRVEC)
       CALL MKL_DCSRTRSV('U','N','N',NTPN,BILUT,IBILUT,JBILUT,TRVEC,TMP(IPAR(23)))
       GOTO 1
      ENDIF   
!---------------------------------------------------------------------------
! If RCI_REQUEST=4, then check if the norm of the next generated vector is
! not zero up to rounding and computational errors. The norm is contained
! in DPAR(7) parameter
!---------------------------------------------------------------------------  
      IF (RCI_REQUEST.EQ.4) THEN
      	IF (DPAR(7).LT.1.0D-12) THEN
      	   GOTO 3
      	ELSE
      	   GOTO 1
      	ENDIF
!---------------------------------------------------------------------------
! If RCI_REQUEST=anything else, then DFGMRES subroutine failed
! to compute the solution vector: COMPUTED_SOLUTION(N)
!---------------------------------------------------------------------------
      ELSE
      	GOTO 999
      ENDIF      
!---------------------------------------------------------------------------
! Reverse Communication ends here
! Get the current iteration number and the FGMRES solution. (DO NOT FORGET to
! call DFGMRES_GET routine as computed_solution is still containing
! the initial guess!). Request to DFGMRES_GET to put the solution into
! vector COMPUTED_SOLUTION(N) via IPAR(13)
!---------------------------------------------------------------------------
3     IPAR(13)=0
      CALL DFGMRES_GET(NTPN, MATX, RHS, RCI_REQUEST, IPAR, DPAR, TMP, ITERCOUNT)      
      WRITE(*,*)"ITERCOUNT      ",ITERCOUNT

999   WRITE( *,'(A,I2)') 'The solver has returned the ERROR code ',RCI_REQUEST
    
      !WRITE(*,*)"ITERCOUNT      ",ITERCOUNT
    
      DO I=1,NBPOINT
      !WRITE(*,*)RHS(I),B(I),VERR(I)
      END DO
      !STOP 
    
    QG(1:NTPN)=MATX(1:NTPN)
    
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)     
    PRINT *, 'END   TIME  GMRES_MKL: ', CHAR_TIME2,'  ON  ',CHAR_DATE2 
    
    DEALLOCATE(A,MATB,MATX)
    DEALLOCATE(IPIV)
END SUBROUTINE