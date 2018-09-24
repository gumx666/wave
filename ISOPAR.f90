SUBROUTINE isopar_3(xi,cord,sn,snx,sne,aj,sr,dn)
!c     *************************************************
!c     the isoparametric element with 8 nodes: cord(8,3)
!c     calculate the value in the local point (xi,eta)
!c     (1) the shape function  sn(8)
!c     (2) derivatives of shape function snx(8),sne(8)
!c     (3) jacobian matrix aj
!c     (4) position vector: sr(3)
!c     (5) normal vector: dn(3)
!c     *************************************************
!cdp   implicit REAL (a-h,o-z)
      REAL:: cord(3,2),sn(3),snx(3),sne(3),sr(2),dn(2)
      REAL:: xxe(2)
      REAL:: AJ,XI,TANTH

    SN(1)=0.5*XI*(XI-1.)
    SN(2)=1.-XI**2
    SN(3)=0.5*XI*(XI+1.)
    
    SNX(1)=0.5*(2.*XI-1.)
    SNX(2)=-2.*XI
    SNX(3)=0.5*(1.+2.*XI)
    
      do 20 i=1,2
      sr(i)=0.
      xxe(i)=0.
      do 20 j=1,3
!c     xxe(1,*): gradient of x to xi (vector)
!c     xxe(2,*): gradient of x to eta (vector)
      sr(i)=sr(i)+cord(j,i)*sn(j)
   20 xxe(i)=xxe(i)+cord(j,i)*snx(j)
    
        
!c     jacobian value
  
      aj=xxe(1)**2+xxe(2)**2
      if(aj.le.0.0) then
 !      write(*,*) xi,eta
 !      write(*,*) aj 
 !      do 1000 i=1,9
 !1000  write(*,501) i,cord(i,1),cord(i,2),cord(i,3)
 501   FORMAT(I4,3F10.5)
       !write(*,*) 'aj=',aj
       !write(*,*) sr
       !write(*,*)
       !write(*,*)
!c      read(*,*)
       !pause 'warning: the jacobian value is letter then zero'
      endif
      aj=sqrt(aj)

!     unit normal vector
      IF(ABS(XXE(1)).LE.1.E-5) THEN
        TANTH=0.
      ELSE
        TANTH=XXE(2)/XXE(1)
      END IF
      DN(2)=-SQRT(1./(1.+TANTH**2))
      DN(1)=-TANTH*DN(2)
   
      return
END SUBROUTINE

SUBROUTINE isopar_9(xi,eta,cord,sn,snx,sne,aj,sr,dn)
!c     *************************************************
!c     the isoparametric element with 8 nodes: cord(8,3)
!c     calculate the value in the local point (xi,eta)
!c     (1) the shape function  sn(8)
!c     (2) derivatives of shape function snx(8),sne(8)
!c     (3) jacobian matrix aj
!c     (4) position vector: sr(3)
!c     (5) normal vector: dn(3)
!c     *************************************************
!cdp   implicit REAL (a-h,o-z)
      REAL:: cord(9,3),sn(9),snx(9),sne(9),sr(3),dn(3)
      REAL:: plx(9),ple(9),xxe(2,3),cc(3)
      REAL:: AJ,XI,ETA
      data plx /-1.,0.,1.,1.,1.,0.,-1.,-1.,0/
      data ple /-1.,-1.,-1.,0.,1.,1.,1.,0.,0/

      do 5 i=1,7,2
	    b=plx(i)*xi
	    c=ple(i)*eta
		sn(i)=0.25*b*c*(1.+b)*(1.+c)
		snx(i)=0.25*plx(i)*c*(1+c)*(1+2*b)
		sne(i)=0.25*ple(i)*b*(1+b)*(1+2*c)		
5	continue

	do 10 i=2,6,4
		c=ple(i)*eta
		sn(i)=0.5*(1-xi*xi)*(1+c)*c
		snx(i)=-xi*(1+c)*c
		sne(i)=0.5*(1-xi*xi)*ple(i)*(1+2*c)
10	continue
	do 11 i=4,8,4
		b=plx(i)*xi
		sn(i)=0.5*(1-eta*eta)*(1+b)*b
		snx(i)=0.5*(1-eta*eta)*plx(i)*(1+2*b)
		sne(i)=-eta*(1+b)*b
11	continue
		sn(9)=(1-xi*xi)*(1-eta*eta)
		snx(9)=-2.*xi*(1-eta*eta)
		sne(9)=-2.*eta*(1-xi*xi)
    
      do 20 i=1,3
      cc(i)=0.
      sr(i)=0.
      xxe(1,i)=0.
      xxe(2,i)=0.
      do 20 j=1,9
!c     xxe(1,*): gradient of x to xi (vector)
!c     xxe(2,*): gradient of x to eta (vector)
      sr(i)=sr(i)+cord(j,i)*sn(j)
      xxe(1,i)=xxe(1,i)+cord(j,i)*snx(j)
   20 xxe(2,i)=xxe(2,i)+cord(j,i)*sne(j)
  
      do 30 i=1,3
      cc(1)=cc(1)+xxe(1,i)*xxe(1,i)
      cc(2)=cc(2)+xxe(2,i)*xxe(2,i)
      cc(3)=cc(3)+xxe(1,i)*xxe(2,i)
   30 continue
   
!c     jacobian value
  
      aj=cc(1)*cc(2)-cc(3)*cc(3)
      if(aj.le.0.0) then
 !      write(*,*) xi,eta
 !      write(*,*) aj 
 !      do 1000 i=1,9
 !1000  write(*,501) i,cord(i,1),cord(i,2),cord(i,3)
 501   FORMAT(I4,3F10.5)
       write(*,*) 'aj=',aj
       !write(*,*) sr
       !write(*,*)
       !write(*,*)
!c      read(*,*)
       pause 'warning: the jacobian value is letter then zero'
      endif
      aj=sqrt(aj)
!c
!c     unit normal vector      
!c
      dn(1)=xxe(1,2)*xxe(2,3)-xxe(1,3)*xxe(2,2)
      dn(2)=xxe(1,3)*xxe(2,1)-xxe(1,1)*xxe(2,3)
      dn(3)=xxe(1,1)*xxe(2,2)-xxe(1,2)*xxe(2,1)      
      do 40 i=1,3
   40 dn(i)=dn(i)/aj
   
      return
END SUBROUTINE

SUBROUTINE ISOPAR_4(XI,ETA,CORDA,SN,SNX,SNE,AJ,SR,DN)
!     *************************************************
!     THE ISOPARAMETRIC ELEMENT WITH 4 NODES: CORDA(4,3)
!     CALCULATE THE VALUE IN THE LOCAL POINT (XI,ETA)
!     (1) THE SHAPE FUNCTION  SN(4)
!     (2) DERIVATIVES OF SHAPE FUNCTION SNX(4),SNE(4)
!     (3) JACOBIAN MATRIX AJ
!     (4) POSITION VECTOR: SR(3)
!     (5) NORMAL VECTOR: DN(3)
!     *************************************************
    REAL:: XI,ETA,AJ
    REAL:: CORDA(4,3),SN(4),SNX(4),SNE(4),SR(3),DN(3),WFAJ(4)
    REAL:: PLX(4),PLE(4),XXE(2,3),CC(3)
    DATA PLX /-1.,1.,1.,-1./
    DATA PLE /-1.,-1.,1.,1./

!	SN=N(4),SNX=NXI(4),SNE=NETA(4)
    DO I=1,4
	    B=PLX(I)*XI
		C=PLE(I)*ETA
		SN(I)=.25*(1.+B)*(1.+C)
		SNX(I)=.25*PLX(I)*(1.+C)
    	SNE(I)=.25*PLE(I)*(1.+B)
    END DO
!    
    DO I=1,3
		CC(I)=0.
		SR(I)=0.
		XXE(1,I)=0.
		XXE(2,I)=0.
    END DO

	DO I=1,3
		DO J=1,4
!			XXE(1,*): XXI,YXI,ZXI
!			XXE(2,*): XETA,YETA,ZETA
			SR(I)=SR(I)+CORDA(J,I)*SN(J)
			XXE(1,I)=XXE(1,I)+CORDA(J,I)*SNX(J)
   			XXE(2,I)=XXE(2,I)+CORDA(J,I)*SNE(J)
        END DO
    END DO

!	CC(1)为X,Y,Z对XI的导数的平方和
!	CC(2)为X,Y,Z对ETA的导数的平方和
!	CC(3)为XXI*XETA+YXI*YETA+ZXI*ZETA   
    DO I=1,3
		CC(1)=CC(1)+XXE(1,I)*XXE(1,I)
	    CC(2)=CC(2)+XXE(2,I)*XXE(2,I)
		CC(3)=CC(3)+XXE(1,I)*XXE(2,I)
    END DO
     
!   JACOBIAN VALUE
    AJ=CC(1)*CC(2)-CC(3)*CC(3)

    IF(AJ.LE.0.0) THEN
        !GOTO 50
	    WRITE(*,*) XI,ETA
		WRITE(*,*) 
		DO I=1,4
 		    WRITE(*,*) I,CORDA(I,1),CORDA(I,2),CORDA(I,3)
		END DO
        WRITE(*,*) "4 NODES PANEL"
        WRITE(*,*) 'AJ=',AJ
		WRITE(*,*) SR
		PAUSE 'WARNING: THE JACOBIAN VALUE IS LETTER THEN ZERO'
    ENDIF

    AJ=SQRT(AJ)

!	NORMAL VECTOR      
!	DN(1)=YXI*ZETA-ZXI*YETA
!	DN(2)=ZXI*XETA-XXI*ZETA
!	DN(3)=XXI*YETA-XETA*YXI
    DN(1)=XXE(1,2)*XXE(2,3)-XXE(1,3)*XXE(2,2)
    DN(2)=XXE(1,3)*XXE(2,1)-XXE(1,1)*XXE(2,3)
    DN(3)=XXE(1,1)*XXE(2,2)-XXE(1,2)*XXE(2,1)       
!	UNIT NORMAL VECTOR
    DO I=1,3
        DN(I)=DN(I)/AJ
    END DO  
    50  CONTINUE
    
    RETURN
END SUBROUTINE