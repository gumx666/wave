SUBROUTINE COFAHS_8(CORD,CONM,PS,AIJ,HIJ)
!    ***********************************************      
!     CALCULATE THE AIJ AND HIJ TERMS IN THE MATRIX   
!      THE PS IS ONE OF THE NODES OF ELEMENT          4---7---3
!      CORD(8,3) THE CORD. OF 8 NODES (CLOCKWISE)     |       |
!      CONM(8,3) THE NORMAL VECTOR OF 8 NODES         8       6
!      PS(3) THE CORD. OF FIELE POINT                 |       |
!      AIJ(8): PANEL SOURCE EFFECT                    1---5---2
!      HIJ(8): PANEL DIPOLE EFFECT
!   ************************************************
!   IMPLICIT REAL(A-H,O-Z)
    REAL:: XI1,XI2,XI3,XI4,ETA1,ETA2,ETA3,ETA4
    REAL:: CORD(8,3),CONM(8,3),PS(3),AIJ(8),HIJ(8,4)
    REAL:: AIJ1(8),HIJ1(8,3),CORD1(8,3),CONM1(8,3)
    REAL:: AY(8),HY(8,3),RH(4),RH3(4),BNR(4),BNR0(4),AJ(4)
    REAL:: SNX(8,4),SNE(8,4),SM(8,4)
    REAL:: SN1(8),SN2(8),SN3(8),SN(8,4),SN4(8)
    EQUIVALENCE (SN(1,1),SN1),(SN(1,2),SN2),(SN(1,3),SN3),(SN(1,4),SN4)
    REAL:: SR1(3),SR2(3),SR3(3),SR(3,4),SR4(3)
    EQUIVALENCE (SR(1,1),SR1),(SR(1,2),SR2),(SR(1,3),SR3),(SR(1,4),SR4)
    REAL:: DN1(3),DN2(3),DN3(3),DN(3,4),DN4(3)
    !EQUIVALENCE (DN(1,1),DN1),(DN(1,2),DN2),(DN(1,3),DN3),(DN(1,4),DN4)
    REAL:: FP_12(4,3),FP_02(4,3)
    REAL:: IN1(4),IN2(3)
    REAL:: HN(8,4),PSR(3,4),RP(4)
    REAL:: P0,P1,P2,DXI,DET 
    REAL:: SN0(8,4),DN0(3,4),RP0(4),PSR0(3,4),HN0(8,4),AJ0(4),SR0(3,4)
    REAL:: SNX0(8,4),SNE0(8,4),JP(3,4),AJP(4),JP0(3,4),AJP0(4),XX1(3,2),XX2(3,2)
    REAL:: D1,D2,D3,D4
    COMMON NML,NBP,NFP,NBE,NFE,NP,MP,NF,RADIUS
    COMMON /GWA/GA(400),GW(400)
    INTEGER:: J
    
      D1=0.9
      D2=D1
      D3=D2
      D4=D3

    R1=0.
    ISN=0
    DO I=1,8
        T=(CORD(I,1)-PS(1))**2+(CORD(I,2)-PS(2))**2+(CORD(I,3)-PS(3))**2
        T=SQRT(T)
        IF(T.LT.1.E-5) ISN=I
        AIJ1(I)=0.
        HIJ1(I,:)=0.
    END DO
!    
    IF(ISN.EQ.0) THEN
        PAUSE 'WARNING: ISN=0, THE SOURCE POINT OUT OF THE ELEMENT'
        STOP
    ENDIF	
!      
    NGUS=21
    !NGUS=21
!   
    KGUS=NGUS*(NGUS-1)/2
!   GAUSS-LEGENDRE NUMERICAL INTEGRATION
	SELECT CASE(ISN)
    CASE(1)
        DXI=-D1;DET=-D3
    CASE(2)
        DXI=D2;DET=-D3
    CASE(3)
        DXI=D2;DET=D4
    CASE(4)
        DXI=-D1;DET=D4
    CASE(5)
        DXI=0.;DET=-D3
    CASE(6)
        DXI=D2;DET=0.
    CASE(7)
        DXI=0.;DET=D4
    CASE(8)
        DXI=-D1;DET=0.
    END SELECT       
            
    DO IX=1,NGUS
        P1=0.5*(1.0+GA(KGUS+IX-1))
        DO IY=1,NGUS
            P2=0.5*(1.0+GA(KGUS+IY-1))
     
            XI1=(1.-P1)*(DXI)+P1*(1-P2)*(-1.)+P1*P2*(-1.)
            ETA1=(1.-P1)*(DET)+P1*(1-P2)*(1.)+P1*P2*(-1.)
            XI2=(1.-P1)*(DXI)+P1*(1-P2)*(-1.)+P1*P2*(1.)
            ETA2=(1.-P1)*(DET)+P1*(1-P2)*(-1.)+P1*P2*(-1.)
            XI3=(1.-P1)*(DXI)+P1*(1-P2)*(1.)+P1*P2*(1.)
            ETA3=(1.-P1)*(DET)+P1*(1-P2)*(-1.)+P1*P2*(1.)
            XI4=(1.-P1)*(DXI)+P1*(1-P2)*(1.)+P1*P2*(-1.)
            ETA4=(1.-P1)*(DET)+P1*(1-P2)*(1.)+P1*P2*(1.)

            CALL ISOPAR_8(XI1,ETA1,CORD,SN1,SNX,SNE,AJ(1),SR1,DN(:,1))
            CALL ISOPAR_8(XI2,ETA2,CORD,SN2,SNX,SNE,AJ(2),SR2,DN(:,2))
            CALL ISOPAR_8(XI3,ETA3,CORD,SN3,SNX,SNE,AJ(3),SR3,DN(:,3))
	        CALL ISOPAR_8(XI4,ETA4,CORD,SN4,SNX,SNE,AJ(4),SR4,DN(:,4))         

            DO J=1,4
                RH(J)=0.
                DO I=1,3
                    RH(J)=RH(J)+(SR(I,J)-PS(I))*(SR(I,J)-PS(I))
                END DO
                RH(J)=SQRT(RH(J))
                RH3(J)=0.0-RH(J)**3                                       !SOLVE FOR R AND R^3(XI,ETA)
            END DO
!           WRITE(*,*)P1/RH(1),1./RP(1)
!           NML=1 USING THE NORMAL VECTOR ACCURATELY ON ISOPAR. CURVE SURFACE
!           ELSE USING THE ISOPARAMETER NORMAL VECTOR BY 8 NODE NOMAL
!      
!           WRITE(*,*)RP(2),RH(2)/P1

            IF(NML.NE.1) THEN
                DO 50 J=1,4
                    DO 50 I1=1,3  
                        DN(I1,J)=0.
                        DO 50 I2=1,8	
                        50  DN(I1,J)=DN(I1,J)+CONM(I2,I1)*SN(I2,I1)          !SOLVE FOR DN
            END IF

            DO I=1,8
                AY(I)=0.
                HY(I,:)=0.
            END DO
         
            DO J=1,8
                DO K=1,4
                    IF(K.EQ.1) AM=ABS(-1.-DXI)
                    IF(K.EQ.2) AM=ABS(-1.-DET)
                    IF(K.EQ.3) AM=ABS(1.-DXI)
                    IF(K.EQ.4) AM=ABS(1.-DET)

                    AY(J)=AY(J)+GW(KGUS+IY-1)*SN(J,K)*AJ(K)*(p1/rh(k))*2.*AM/4.
                    !AY(J)=AY(J)+GW(KGUS+IY-1)*SN(J,K)*AJ(K)/RP(k)*2.*AM/4.
                END DO
            END DO

            DO J=1,8
                AIJ1(J)=AIJ1(J)+GW(KGUS+IX-1)*AY(J)
                HIJ1(J,1)=0.
                HIJ1(J,2)=0.
                HIJ1(J,3)=0.
            END DO
        END DO
    END DO 

    DO J=1,8
      	AIJ(J)=AIJ1(J)
        HIJ(J,1)=HIJ1(J,1)
        HIJ(J,2)=HIJ1(J,2)
        HIJ(J,3)=HIJ1(J,3)
    END DO	

    RETURN
END SUBROUTINE

SUBROUTINE COFAH_8(CORD,CONM,PS,AIJ,HIJ) 
!C    ***********************************************      
!C     CALCULATE THE AIJ AND HIJ TERMS IN THE MATRIX    4---8---3  
!C      CORD(8,3) THE CORD. OF 8 NODES (CLOCKWISE)      |       |
!C      CONM(8,3) THE NORMAL VECTOR OF 8 NODES          5       7
!C      PS(3) THE CORD. OF FIELE POINT                  |       |
!C      AIJ(8): PANEL SOURCE EFFECT                     1---6---2
!C      HIJ(8): PANEL DIPOLE EFFECT
!C   ************************************************
!CDP   IMPLICIT REAL(A-H,O-Z)
      REAL:: CORD(8,3),PS(3),CONM(8,3),AIJ(8),HIJ(8,4),PR(3)
      REAL:: R(4),SRR(3),AY(8),HY(8,4)
      REAL:: SN(8),SNX(8),SNE(8),SR(3),DN(3)
      REAL:: AJ,XI,ETA
      REAL:: PI4=4.*3.1415926
      COMMON NML,NBP,NFP,NBE,NFE,NP,MP,NF,RADIUS
      COMMON /GWA/GA(400),GW(400)
!C      
      r1=0.
      do 5 i=1,8
      aij(i)=0.
    5 hij(i,:)=0.
      
      do 10 i=1,3
      pr(i)=0.
      do 12 j=1,8
   12 pr(i)=pr(i)+cord(j,i)
      pr(i)=pr(i)/8.-ps(i)
      r1=r1+pr(i)*pr(i)
   10 continue
      r1=sqrt(r1)      
      
      do 15 i=1,4
   15 r(i)=0.
      do 20 i=1,3
      do 22 j=3,7,2
   22 r(i)=r(i)+(cord(j,i)-cord(j-2,i))**2
      r(4)=r(4)+(cord(7,i)-cord(1,i))**2
   20 continue

      rv=0.
      do 25 i=1,4
   25 rv=rv+sqrt(r(i))
      rv=rv/8.
      r0=r1/rv

      ngus=21
      if(r0.lt.1.6) goto 30
      ngus=11
      if(r0.lt.1.9) goto 30
      ngus=9
      if(r0.lt.2.4) goto 30
      ngus=7
      if(r0.lt.3.4) goto 30
      ngus=5
      if(r0.lt.6.4) goto 30
      ngus=3
   30 continue

      ngus=21
      K=NGUS*(NGUS-1)/2

!C     GAUSS-LEGENDRE NUMERICAL INTEGRATION
     
      DO 40 IX=1,NGUS
        XI=GA(K+IX-1)
        DO 40 IY=1,NGUS
            ETA=GA(K+IY-1)
      
            CALL ISOPAR_8(XI,ETA,CORD,SN,SNX,SNE,AJ,SR,DN)
      
            RH=0.
            DO 45 I=1,3
                SRR(I)=SR(I)-PS(I)   
         45 RH=RH+SRR(I)*SRR(I)
            RH=SQRT(RH)
! LEE ORIGIN
	        RH3=RH**3
!     
!     NML=1 USING THE NORMAL VECTOR ACCURATELY ON ISOPAR. CURVE SURFACE
!     ELSE USING THE ISOPARAMETER NORMAL VECTOR BY 8 NODE NOMAL
!      
      IF(NML.NE.1) THEN
       DO 48 I1=1,3  
       DN(I1)=0.
       DO 48 I2=1,8	
   48  DN(I1)=DN(I1)+CONM(I2,I1)*SN(I2)
      ENDIF
!       
      BNR=SRR(1)*DN(1)+SRR(2)*DN(2)+SRR(3)*DN(3)
!      
      DO J=1,8
        AY(J)=GW(K+IY-1)*SN(J)*AJ/RH
        HY(J,1)=GW(K+IY-1)*SN(J)*SRR(1)*AJ/RH3
        HY(J,2)=GW(K+IY-1)*SN(J)*SRR(2)*AJ/RH3      
        HY(J,3)=GW(K+IY-1)*SN(J)*SRR(3)*AJ/RH3
           
        HY(J,4)=GW(K+IY-1)*SN(J)*3.*SRR(1)**2*AJ/RH**5&
                -GW(K+IY-1)*SN(J)*AJ/RH**3    !DXX
      END DO
!C
      DO 40 J=1,8
        AIJ(J)=AIJ(J)+GW(K+IX-1)*AY(J)
        HIJ(J,1)=HIJ(J,1)+GW(K+IX-1)*HY(J,1)
        HIJ(J,2)=HIJ(J,2)+GW(K+IX-1)*HY(J,2)
        HIJ(J,3)=HIJ(J,3)+GW(K+IX-1)*HY(J,3)
        HIJ(J,4)=HIJ(J,4)+GW(K+IX-1)*HY(J,4)

      40 CONTINUE

      RETURN
END SUBROUTINE

SUBROUTINE isopar_8(xi,eta,cord,sn,snx,sne,aj,sr,dn)
!     *************************************************
!     the isoparametric element with 8 nodes: cord(8,3)
!     calculate the value in the local point (xi,eta)
!     (1) the shape function  sn(8)
!     (2) derivatives of shape function snx(8),sne(8)
!     (3) jacobian matrix aj
!     (4) position vector: sr(3)
!     (5) normal vector: dn(3)
!     *************************************************
!     implicit REAL (a-h,o-z)
      REAL:: cord(8,3),sn(8),snx(8),sne(8),sr(3),dn(3)
      REAL:: plx(8),ple(8),xxe(2,3),cc(3)
      REAL:: AJ,XI,ETA,D1,D2,D3,D4

      D1=0.9
      D2=D1
      D3=D2
      D4=D3

      !SHAPE FUNCTION
      SN(1)=1./(D1+D2)/(D3+D4)*(D1-XI)*(D3-ETA)*&
            (-1.-XI/D2-ETA/D4)
      SN(2)=1./(D1+D2)/(D3+D4)*(D2+XI)*(D3-ETA)*&
            (-1.+XI/D1-ETA/D4)
      SN(3)=1./(D1+D2)/(D3+D4)*(D2+XI)*(D4+ETA)*&
            (-1.+XI/D1+ETA/D3)
      SN(4)=1./(D1+D2)/(D3+D4)*(D1-XI)*(D4+ETA)*&
            (-1.-XI/D2+ETA/D3)

      SN(5)=1./D1/D2/(D3+D4)*(D1-XI)*(D2+XI)*(D3-ETA)
      SN(6)=1./D3/D4/(D1+D2)*(D2+XI)*(D4+ETA)*(D3-ETA)
      SN(7)=1./D1/D2/(D3+D4)*(D1-XI)*(D2+XI)*(D4+ETA)
      SN(8)=1./D3/D4/(D1+D2)*(D1-XI)*(D4+ETA)*(D3-ETA)

      !DERIATIVES OF SHAPE FUNCTION
      SNX(1)=-1./(D1+D2)/(D3+D4)*(D3-ETA)*(-1.-XI/D2-ETA/D4)&
             -1./D2*1./(D1+D2)/(D3+D4)*(D1-XI)*(D3-ETA)
      SNX(2)=1./(D1+D2)/(D3+D4)*(D3-ETA)*(-1.+XI/D1-ETA/D4)&
             +1./D1*1./(D1+D2)/(D3+D4)*(D2+XI)*(D3-ETA)
      SNX(3)=1./(D1+D2)/(D3+D4)*(D4+ETA)*(-1.+XI/D1+ETA/D3)&
             +1./D1*1./(D1+D2)/(D3+D4)*(D2+XI)*(D4+ETA)
      SNX(4)=-1./(D1+D2)/(D3+D4)*(D4+ETA)*(-1.-XI/D2+ETA/D3)&
             -1./D2*1./(D1+D2)/(D3+D4)*(D1-XI)*(D4+ETA)
      SNX(5)=-1./D1/D2/(D3+D4)*(D2+XI)*(D3-ETA)&
             +1./D1/D2/(D3+D4)*(D1-XI)*(D3-ETA)  
      SNX(6)=1./D3/D4/(D1+D2)*(D4+ETA)*(D3-ETA)  
      SNX(7)=-1./D1/D2/(D3+D4)*(D2+XI)*(D4+ETA)&
             +1./D1/D2/(D3+D4)*(D1-XI)*(D4+ETA)  
      SNX(8)=-1./D3/D4/(D1+D2)*(D4+ETA)*(D3-ETA)   

      SNE(1)=-1./(D1+D2)/(D3+D4)*(D1-XI)*(-1.-XI/D2-ETA/D4)&
             -1./D4*1./(D1+D2)/(D3+D4)*(D1-XI)*(D3-ETA)
      SNE(2)=-1./(D1+D2)/(D3+D4)*(D2+XI)*(-1.+XI/D1-ETA/D4)&
             -1./D4*1./(D1+D2)/(D3+D4)*(D2+XI)*(D3-ETA)
      SNE(3)=+1./(D1+D2)/(D3+D4)*(D2+XI)*(-1.+XI/D1+ETA/D3)&
             +1./D3*1./(D1+D2)/(D3+D4)*(D2+XI)*(D4+ETA)
      SNE(4)=+1./(D1+D2)/(D3+D4)*(D1-XI)*(-1.-XI/D2+ETA/D3)&
             +1./D3*1./(D1+D2)/(D3+D4)*(D1-XI)*(D4+ETA)
      SNE(5)=-1./D1/D2/(D3+D4)*(D1-XI)*(D2+XI)
      SNE(6)=1./D3/D4/(D1+D2)*(D2+XI)*(D3-ETA)&
             -1./D3/D4/(D1+D2)*(D2+XI)*(D4+ETA)
      SNE(7)=1./D1/D2/(D3+D4)*(D1-XI)*(D2+XI)
      SNE(8)=1./D3/D4/(D1+D2)*(D1-XI)*(D3-ETA)&
             -1./D3/D4/(D1+D2)*(D1-XI)*(D4+ETA)
    
      do 20 i=1,3
        cc(i)=0.
        sr(i)=0.
        xxe(1,i)=0.
        xxe(2,i)=0.
      do 20 j=1,8
!     xxe(1,*): gradient of x to xi (vector)
!     xxe(2,*): gradient of x to eta (vector)
      sr(i)=sr(i)+cord(j,i)*sn(j)
      xxe(1,i)=xxe(1,i)+cord(j,i)*snx(j)
   20 xxe(2,i)=xxe(2,i)+cord(j,i)*sne(j)
  
      !WRITE(*,*)SR(:)

      do 30 i=1,3
        cc(1)=cc(1)+xxe(1,i)*xxe(1,i)
        cc(2)=cc(2)+xxe(2,i)*xxe(2,i)
        cc(3)=cc(3)+xxe(1,i)*xxe(2,i)
   30 continue
   
!c     jacobian value
  
      aj=cc(1)*cc(2)-cc(3)*cc(3)
      if(aj.le.0.0) then
       !write(*,*) xi,eta
       !write(*,*) 
       do 1000 i=1,8
 1000  write(*,501) i,cord(i,1),cord(i,2),cord(i,3)
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

