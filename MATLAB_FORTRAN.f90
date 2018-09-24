SUBROUTINE MATRIX_MATLAB
    USE GREENMOD
    
    REAL*8,DIMENSION(:,:),ALLOCATABLE:: MATA
    REAL*8,DIMENSION(:),ALLOCATABLE:: MATB,MATX
    INTEGER*8 engOpen, engClose,mxCreateDoubleMatrix ! ! 注：在64位机器里，要声明为integer*8
    INTEGER*8 mxGetPr
    INTEGER*8 engPutVariable,engGetVariable,engEvalString,engGetMatrix
    INTEGER*8 ep, ASBM, BSBM,XSBM, reg, status,LC ! 相当于指针的作用
    
	CHARACTER*8 CHAR_TIME1,CHAR_TIME2
	CHARACTER*9 CHAR_DATE1,CHAR_DATE2

    ALLOCATE(MATA(NTPN,NTPN),MATB(NTPN),MATX(NTPN))
     
    DO I=1,NTPN
        MATA(I,1:NTPN)=MAL(I,1:NTPN)
    END DO
    MATB(1:NTPN)=VERR(1:NTPN)
    
   
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
    PRINT *, 'END   TIME  INPUT: ', CHAR_TIME2,'  ON  ',CHAR_DATE2 
    
    ep = engOpen('matlab')  !ep= engOpen('\0')   ! 打开引擎

    !WRITE(*,*)"OKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK"
    
    asbm =mxCreateDoubleMatrix(NTPN,NTPN,0)    ! 创建双精度数组第三个参数0表示实数;否则用mxCOMPLEX表示复数
    Bsbm =mxCreateDoubleMatrix(NTPN,1,0)
    !Xsbm =mxCreateDoubleMatrix(NTPN,1,0)
    
    call mxCopyReal8ToPtr(MATA,mxGetPr(ASBM),NTPN*NTPN) ! 给新创建的数组赋值,把GG中的值转化为matlab能识别的语言存在asbm中  
    call mxCopyReal8ToPtr(MATB,mxGetPr(BSBM),NTPN) 
    
    status =engPutVariable(ep, 'ASBM', ASBM) ! 将新创建的数组植入matlab引擎（工作区）
    !WRITE(*,*)STATUS
    
    status =engPutVariable(ep, 'BSBM', BSBM) ! 将新创建的数组植入matlab引擎（工作区）

    WRITE(*,*)"COPY IN OVER"
    CALL TIME(CHAR_TIME2)
	CALL DATE(CHAR_DATE2)
    PRINT *, 'END   TIME  COPYIN: ', CHAR_TIME2,'  ON  ',CHAR_DATE2 


    mmm=engEvalString(ep, 'xsbmm = ASBM\BSBM')
    lc=engGetVariable(ep,'xsbmm') 

    call mxCopyPtrToReal8(mxGetPr(lc),MATX,NTPN) 
    
    WRITE(*,*)"COPY OUT OVER"
    !stop

    
    call mxDestroyArray(asbm)
    call mxDestroyArray(bsbm)
    status = engClose(ep)
    
    OPEN(1,FILE='PHI.DAT')
    
    DO I=1,NTPN
        WRITE(1,*)MATX(I)
        QG(I)=MATX(I)
    END DO
    CLOSE(1)
    !QG(NBPOINT+1:NBPOINT+NFY)=0.
    
    !DEALLOCATE(MAL,VERR)
    DEALLOCATE(MATA,MATB,MATX)
    
    
    
END SUBROUTINE