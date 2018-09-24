SUBROUTINE COFAH_9(CORD,CONM,PS,AIJ,HIJ) 
!C    ***********************************************      
!C     CALCULATE THE AIJ AND HIJ TERMS IN THE MATRIX    7---6---5  
!C      CORD(8,3) THE CORD. OF 8 NODES (CLOCKWISE)      |       |
!C      CONM(8,3) THE NORMAL VECTOR OF 8 NODES          8   9   4
!C      PS(3) THE CORD. OF FIELE POINT                  |       |
!C      AIJ(8): PANEL SOURCE EFFECT                     1---2---3
!C      HIJ(8): PANEL DIPOLE EFFECT
!C   ************************************************
      !IMPLICIT REAL(A-H,O-Z)
      !IMPLICIT NONE
      REAL:: CORD(9,3),PS(3),CONM(9,3),AIJ(9),HIJ(9,8),PR(3)
      REAL:: R(4),SRR(3),AY(9),HY(9,8)
      REAL:: SN(9),SNX(9),SNE(9),SR(3),DN(3)
      REAL:: AJ,XI,ETA
      REAL:: R1,RV,R0,RH,RH3
      REAL:: BNR
      REAL:: PI4=4.*3.1415926
      INTEGER:: I,J,K,NGUS,IX,IY
      COMMON NML,NBP,NFP,NBE,NFE,NP,MP,NF,RADIUS
      COMMON /GWA/GA(400),GW(400)
!C      
      R1=0.
      DO 5 I=1,9
      AIJ(I)=0.
    5 HIJ(I,:)=0.
!C      
      DO 10 I=1,3
      PR(I)=0.
      DO 12 J=1,9
   12 PR(I)=PR(I)+CORD(J,I)
      PR(I)=PR(I)/9.-PS(I) 
      R1=R1+PR(I)*PR(I)
   10 CONTINUE
      R1=SQRT(R1)      

      DO 15 I=1,4
   15 R(I)=0.
      DO 20 I=1,3

      R(1)=R(1)+(CORD(1,I)-CORD(3,I))**2
      R(2)=R(2)+(CORD(3,I)-CORD(5,I))**2
      R(3)=R(3)+(CORD(5,I)-CORD(7,I))**2
      R(4)=R(4)+(CORD(7,I)-CORD(1,I))**2

   20 CONTINUE

      RV=0.
      DO 25 I=1,4
   25 RV=RV+SQRT(R(I))
      RV=RV/8.
      R0=R1/RV

      NGUS=21
      !NGUS=11
      IF(R0.LT.1.6) GOTO 30
      NGUS=11
      IF(R0.LT.1.9) GOTO 30
      NGUS=9
      IF(R0.LT.2.4) GOTO 30
      NGUS=7
      IF(R0.LT.3.4) GOTO 30
      NGUS=5
      IF(R0.LT.6.4) GOTO 30
      NGUS=3
   30 CONTINUE
   
      !NGUS=11
      K=NGUS*(NGUS-1)/2

!C     GAUSS-LEGENDRE NUMERICAL INTEGRATION
     
      DO 40 IX=1,NGUS
      XI=GA(K+IX-1)
      DO 40 IY=1,NGUS
      ETA=GA(K+IY-1)
      
      CALL ISOPAR_9(XI,ETA,CORD,SN,SNX,SNE,AJ,SR,DN)
      
      RH=0.
      DO 45 I=1,3
      SRR(I)=SR(I)-PS(I)   
   45 RH=RH+SRR(I)**2
      RH=SQRT(RH)
! LEE ORIGIN
	  RH3=RH**3
!       
      BNR=SRR(1)*DN(1)+SRR(2)*DN(2)+SRR(3)*DN(3)
!      

      DO J=1,9
        !IF(IY.NE.11) THEN
        AY(J)=GW(K+IY-1)*SN(J)*AJ/RH
        HY(J,1)=GW(K+IY-1)*SN(J)*SRR(1)*AJ/RH3   !DX
        HY(J,2)=GW(K+IY-1)*SN(J)*SRR(2)*AJ/RH3   !DY      
        HY(J,3)=GW(K+IY-1)*SN(J)*SRR(3)*AJ/RH3   !DZ
           
        HY(J,4)=GW(K+IY-1)*SN(J)*3.*SRR(1)**2*AJ/RH**5&
                -GW(K+IY-1)*SN(J)*AJ/RH**3    !DXX

        !HY(J,4)=3.*SRR(1)**2*AJ/(RH**5)
        !HY(J,4)=HY(J,4)*GW(K+IY-1)*SN(J)-GW(K+IY-1)*SN(J)*AJ/(RH3)

        HY(J,5)=3.*SRR(1)*SRR(2)*AJ/(RH**5)    !DXY
        HY(J,5)=HY(J,5)*GW(K+IY-1)*SN(J)

        HY(J,6)=3.*SRR(1)*SRR(3)*AJ/(RH**5)    !DXZ
        HY(J,6)=HY(J,6)*GW(K+IY-1)*SN(J)   

        HY(J,7)=GW(K+IY-1)*SN(J)*3.*SRR(2)**2*AJ/(RH**5)&
                -GW(K+IY-1)*SN(J)*AJ/RH**3    !DYY
        !HY(J,7)=3.*SRR(2)**2*AJ/(RH**5)
        !HY(J,7)=HY(J,7)*GW(K+IY-1)*SN(J)-GW(K+IY-1)*SN(J)*AJ/(RH3)    !DYY

        HY(J,8)=3.*SRR(2)*SRR(3)*AJ/(RH**5)    !DYZ
        HY(J,8)=HY(J,8)*GW(K+IY-1)*SN(J)    !DYZ
        !END IF
      END DO
!C
      DO 40 J=1,9
        AIJ(J)=AIJ(J)+GW(K+IX-1)*AY(J)
        HIJ(J,1)=HIJ(J,1)+GW(K+IX-1)*HY(J,1)
        HIJ(J,2)=HIJ(J,2)+GW(K+IX-1)*HY(J,2)
        HIJ(J,3)=HIJ(J,3)+GW(K+IX-1)*HY(J,3)
        HIJ(J,4)=HIJ(J,4)+GW(K+IX-1)*HY(J,4)

        HIJ(J,5)=HIJ(J,5)+GW(K+IX-1)*HY(J,5)
        HIJ(J,6)=HIJ(J,6)+GW(K+IX-1)*HY(J,6)
        HIJ(J,7)=HIJ(J,7)+GW(K+IX-1)*HY(J,7)
        HIJ(J,8)=HIJ(J,8)+GW(K+IX-1)*HY(J,8)

      40 CONTINUE

      RETURN
END SUBROUTINE

SUBROUTINE COFAHN_9(CORD,CONM,PS,AIJ,HIJ) 
!C    ***********************************************      
!C     CALCULATE THE AIJ AND HIJ TERMS IN THE MATRIX    7---6---5  
!C      CORD(8,3) THE CORD. OF 8 NODES (CLOCKWISE)      |       |
!C      CONM(8,3) THE NORMAL VECTOR OF 8 NODES          8   9   4
!C      PS(3) THE CORD. OF FIELE POINT                  |       |
!C      AIJ(8): PANEL SOURCE EFFECT                     1---2---3
!C      HIJ(8): PANEL DIPOLE EFFECT
!C   ************************************************
      !IMPLICIT REAL(A-H,O-Z)
      !IMPLICIT NONE
      REAL:: CORD(9,3),PS(3),CONM(9,3),AIJ(9),HIJ(9,8),PR(3)
      REAL:: R(4),SRR(3),AY(9),HY(9,8)
      REAL:: SN(9),SNX(9),SNE(9),SR(3),DN(3)
      REAL:: AJ,XI,ETA
      REAL:: R1,RV,R0,RH,RH3
      REAL:: BNR
      REAL:: PI4=4.*3.1415926
      INTEGER:: I,J,K,NGUS,IX,IY
      COMMON NML,NBP,NFP,NBE,NFE,NP,MP,NF,RADIUS
      COMMON /GWA/GA(400),GW(400)
!C      
      R1=0.
      DO 5 I=1,9
      AIJ(I)=0.
    5 HIJ(I,:)=0.
!C      
      DO 10 I=1,3
      PR(I)=0.
      DO 12 J=1,9
   12 PR(I)=PR(I)+CORD(J,I)
      PR(I)=PR(I)/9.-PS(I) 
      R1=R1+PR(I)*PR(I)
   10 CONTINUE
      R1=SQRT(R1)      

      DO 15 I=1,4
   15 R(I)=0.
      DO 20 I=1,3

      R(1)=R(1)+(CORD(1,I)-CORD(3,I))**2
      R(2)=R(2)+(CORD(3,I)-CORD(5,I))**2
      R(3)=R(3)+(CORD(5,I)-CORD(7,I))**2
      R(4)=R(4)+(CORD(7,I)-CORD(1,I))**2

   20 CONTINUE

      RV=0.
      DO 25 I=1,4
   25 RV=RV+SQRT(R(I))
      RV=RV/8.
      R0=R1/RV

      NGUS=21
      IF(R0.LT.1.6) GOTO 30
      NGUS=11
      IF(R0.LT.1.9) GOTO 30
      NGUS=9
      IF(R0.LT.2.4) GOTO 30
      NGUS=7
      IF(R0.LT.3.4) GOTO 30
      NGUS=5
      IF(R0.LT.6.4) GOTO 30
      NGUS=3
   30 CONTINUE
   
      NGUS=21
      K=NGUS*(NGUS-1)/2

!C     GAUSS-LEGENDRE NUMERICAL INTEGRATION
     
      DO 40 IX=1,NGUS
      XI=GA(K+IX-1)
      DO 40 IY=1,NGUS
      ETA=GA(K+IY-1)
      
      CALL ISOPAR_9(XI,ETA,CORD,SN,SNX,SNE,AJ,SR,DN)
      
      RH=0.
      DO 45 I=1,3
      SRR(I)=SR(I)-PS(I)   
   45 RH=RH+SRR(I)**2
      RH=SQRT(RH)
! LEE ORIGIN
	  RH3=RH**3
!       
      BNR=SRR(1)*DN(1)+SRR(2)*DN(2)+SRR(3)*DN(3)
!      

      DO J=1,9
        !IF(IY.NE.11) THEN
        AY(J)=GW(K+IY-1)*SN(J)*AJ/RH
        HY(J,1)=GW(K+IY-1)*SN(J)*SRR(1)*AJ/RH3   !DX
        HY(J,2)=GW(K+IY-1)*SN(J)*SRR(2)*AJ/RH3   !DY      
        HY(J,3)=GW(K+IY-1)*SN(J)*SRR(3)*AJ/RH3   !DZ
           
        HY(J,4)=GW(K+IY-1)*SN(J)*3.*SRR(1)**2*AJ/RH**5&
                -GW(K+IY-1)*SN(J)*AJ/RH**3    !DXX

        !HY(J,4)=3.*SRR(1)**2*AJ/(RH**5)
        !HY(J,4)=HY(J,4)*GW(K+IY-1)*SN(J)-GW(K+IY-1)*SN(J)*AJ/(RH3)

        HY(J,5)=3.*SRR(1)*SRR(2)*AJ/(RH**5)    !DXY
        HY(J,5)=HY(J,5)*GW(K+IY-1)*SN(J)

        HY(J,6)=3.*SRR(1)*SRR(3)*AJ/(RH**5)    !DXZ
        HY(J,6)=HY(J,6)*GW(K+IY-1)*SN(J)   

        !HY(J,7)=GW(K+IY-1)*SN(J)*3.*SRR(2)**2*AJ/(RH**5)&
        !        -GW(K+IY-1)*SN(J)*AJ/(RH3)    !DYY
        HY(J,7)=3.*SRR(2)**2*AJ/(RH**5)
        HY(J,7)=HY(J,7)*GW(K+IY-1)*SN(J)-GW(K+IY-1)*SN(J)*AJ/(RH3)    !DYY

        HY(J,8)=3.*SRR(2)*SRR(3)*AJ/(RH**5)    !DYZ
        HY(J,8)=HY(J,8)*GW(K+IY-1)*SN(J)    !DYZ
        !END IF
      END DO
!C
      DO 40 J=1,9
        AIJ(J)=AIJ(J)+GW(K+IX-1)*AY(J)
        HIJ(J,1)=HIJ(J,1)+GW(K+IX-1)*HY(J,1)
        HIJ(J,2)=HIJ(J,2)+GW(K+IX-1)*HY(J,2)
        HIJ(J,3)=HIJ(J,3)+GW(K+IX-1)*HY(J,3)
        HIJ(J,4)=HIJ(J,4)+GW(K+IX-1)*HY(J,4)

        HIJ(J,5)=HIJ(J,5)+GW(K+IX-1)*HY(J,5)
        HIJ(J,6)=HIJ(J,6)+GW(K+IX-1)*HY(J,6)
        HIJ(J,7)=HIJ(J,7)+GW(K+IX-1)*HY(J,7)
        HIJ(J,8)=HIJ(J,8)+GW(K+IX-1)*HY(J,8)

      40 CONTINUE

      RETURN
END SUBROUTINE

!SUBROUTINE COFAHS_9(CORD,CONM,PS,AIJ)     
SUBROUTINE COFAHS_9(CORD,CONM,PS,AIJ,HIJ)
!    ***********************************************      
!     CALCULATE THE AIJ AND HIJ TERMS IN THE MATRIX   
!      THE PS IS ONE OF THE NODES OF ELEMENT          7---6---5
!      CORD(8,3) THE CORD. OF 8 NODES (CLOCKWISE)     |       |
!      CONM(8,3) THE NORMAL VECTOR OF 8 NODES         8   9   4
!      PS(3) THE CORD. OF FIELE POINT                 |       |
!      AIJ(8): PANEL SOURCE EFFECT                    1---2---3
!      HIJ(8): PANEL DIPOLE EFFECT
!   ************************************************
!   IMPLICIT REAL(A-H,O-Z)
    REAL:: XI1,XI2,XI3,XI4,ETA1,ETA2,ETA3,ETA4
    REAL:: CORD(9,3),CONM(9,3),PS(3),AIJ(9),HIJ(9,8)
    REAL:: AIJ1(9),HIJ1(9,3),CORD1(9,3),CONM1(9,3)
    REAL:: AY(9),HY(9,3),RH(4),RH3(4),BNR(4),BNR0(4),AJ(4)
    REAL:: SNX(9,4),SNE(9,4),SM(9,4)
    REAL:: SN1(9),SN2(9),SN3(9),SN(9,4),SN4(9)
    EQUIVALENCE (SN(1,1),SN1),(SN(1,2),SN2),(SN(1,3),SN3),(SN(1,4),SN4)
    REAL:: SR1(3),SR2(3),SR3(3),SR(3,4),SR4(3)
    EQUIVALENCE (SR(1,1),SR1),(SR(1,2),SR2),(SR(1,3),SR3),(SR(1,4),SR4)
    REAL:: DN1(3),DN2(3),DN3(3),DN(3,4),DN4(3)
    EQUIVALENCE (DN(1,1),DN1(1)),(DN(1,2),DN2(1)),(DN(1,3),DN3(1)),(DN(1,4),DN4(1))
    EQUIVALENCE (DN(2,1),DN1(2)),(DN(2,2),DN2(2)),(DN(2,3),DN3(2)),(DN(2,4),DN4(2))
    EQUIVALENCE (DN(3,1),DN1(3)),(DN(3,2),DN2(3)),(DN(3,3),DN3(3)),(DN(3,4),DN4(3))
    REAL:: FP_12(4,3),FP_02(4,3)
    REAL:: IN1(4),IN2(3)
    REAL:: HN(9,4),PSR(3,4),RP(4)
    REAL:: P0,P1,P2 
    REAL:: SN0(9,4),DN0(3,4),RP0(4),PSR0(3,4),HN0(9,4),AJ0(4),SR0(3,4)
    REAL:: SNX0(9,4),SNE0(9,4),JP(3,4),AJP(4),JP0(3,4),AJP0(4),XX1(3,2),XX2(3,2)
    COMMON NML,NBP,NFP,NBE,NFE,NP,MP,NF,RADIUS
    COMMON /GWA/GA(400),GW(400)
    INTEGER:: J,ISN
    
    R1=0.
    ISN=0
    DO I=1,9
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
    !NGUS=11
    NGUS=21
!   
    KGUS=NGUS*(NGUS-1)/2
!   GAUSS-LEGENDRE NUMERICAL INTEGRATION
!     
!   CORNER NODE 1, 3, 5, 7
!   ********************************
    IF((MOD(ISN,2).EQ.1).AND.(ISN.NE.9)) THEN
      
    DO I=1,8
        K=MOD(ISN+I-2,8)+1
        DO J=1,3
            CORD1(I,J)=CORD(K,J)
            CONM1(I,J)=CONM(K,J)
        END DO
    END DO

	DO J=1,3
		CORD1(9,J)=CORD(9,J)
		CONM1(9,J)=CONM(9,J)
    END DO
	
    IN1=0
    IN2=0              
    DO IX=1,NGUS
        P1=0.5*(1.0+GA(KGUS+IX-1))
        P0=0.5*(1.0+GA(KGUS))
        P0=0.
        !IN1=0
        DO IY=1,NGUS
            P2=0.5*(1.0+GA(KGUS+IY-1))
     
            XI1=-1.+2.*P1
            ETA1=-1.+2.*P1*P2
            XI2=-1.+2.*P1-2.*P1*P2
            ETA2=-1.+2.*P1 

            CALL ISOPAR_9(XI1,ETA1,CORD1,SN1,SNX(:,1),SNE(:,1),AJ(1),SR1,DN1)
            CALL ISOPAR_9(XI2,ETA2,CORD1,SN2,SNX(:,2),SNE(:,2),AJ(2),SR2,DN2)

            !DN(:,1)=DN1(:)
            !DN(:,2)=DN2(:)
            !DN(:,4)=DN4(:)

            CALL ISOPAR_1_1(P1,P2,CORD1,HN(:,1),PSR(:,1),RP(1))       !RP为距离 ，PSR为坐标点
            CALL ISOPAR_1_2(P1,P2,CORD1,HN(:,2),PSR(:,2),RP(2))

            XI1=-1.+2.*P0
            ETA1=-1.+2.*P0*P2
            XI2=-1.+2.*P0-2.*P0*P2
            ETA2=-1.+2.*P0 

            CALL ISOPAR_9(XI1,ETA1,CORD1,SN0(:,1),SNX0(:,1),SNE0(:,1),AJ0(1),SR0(:,1),DN0(:,1))
            CALL ISOPAR_9(XI2,ETA2,CORD1,SN0(:,2),SNX0(:,2),SNE0(:,2),AJ0(2),SR0(:,2),DN0(:,2))

            !P0=0.
            CALL ISOPAR_1_1(P0,P2,CORD1,HN0(:,1),PSR0(:,1),RP0(1))
            CALL ISOPAR_1_2(P0,P2,CORD1,HN0(:,2),PSR0(:,2),RP0(2))

            DO J=1,2
                DO K=1,9
                    SM(K,J)=SN(K,J)
                    IF(K.EQ.1) SM(K,J)=SM(K,J)-1.                         !SHAPE FUNCTION(N(P1,P2))          
                END DO
            END DO          

            DO J=1,2
                RH(J)=0.
                DO I=1,3
                    RH(J)=RH(J)+(SR(I,J)-PS(I))*(SR(I,J)-PS(I))
                END DO
                RH(J)=SQRT(RH(J))
                RH3(J)=0.0-RH(J)**3                                       !SOLVE FOR R AND R^3(XI,ETA)
            END DO
            !WRITE(*,*)P1/RH(1),1./RP(1)


            DO J=1,2       
                BNR(J)=0.
                BNR0(J)=0.
                DO I=1,3
                    BNR(J)=BNR(J)+PSR(I,J)*DN(I,J)
                    BNR0(J)=BNR0(J)+PSR0(I,J)*DN0(I,J)                            !R(P1,P2)*DN
                    !BNR(J)=BNR(J)+(SR(I,J)-PS(I))*DN(I,J)
                END DO
            END DO

            DO I=1,9
                AY(I)=0.
                HY(I,:)=0.
            END DO
       
            !IN1=0.  
            DO J=1,9
                DO K=1,2
                    AY(J)=AY(J)+GW(KGUS+IY-1)*SN(J,K)*AJ(K)/RP(K)
!   LEE ORIGIN
!   W=(SM(J,K)/RH(K))*(BNR(K)/RH(K))*AJ(K)*(P1/RH(K))
                    IF(J.NE.1) THEN        !w=sm(j,k)*bnr(k)*aj(k)*p1/rh3(k)
                        W=-HN(J,K)*AJ(K)*PSR(1,K)/RP(K)**3*2.   !(X-XI)/R**3)
                        HY(J,1)=HY(J,1)+GW(KGUS+IY-1)*W 
                        
                        W=-HN(J,K)*AJ(K)*PSR(2,K)/RP(K)**3*2.   !(Y-YI)/R**3
                        HY(J,2)=HY(J,2)+GW(KGUS+IY-1)*W 
                        
                        W=-HN(J,K)*AJ(K)*PSR(3,K)/RP(K)**3*2.   !(Z-ZI)/R**3
                        HY(J,3)=HY(J,3)+GW(KGUS+IY-1)*W 
                    ELSE
                        !IF(IX.NE.1) THEN
                        FP_12(K,1)=-SN(J,K)*PSR(1,K)*AJ(K)/(RP(K)**3)*2. 
                        FP_12(K,2)=-SN(J,K)*PSR(2,K)*AJ(K)/(RP(K)**3)*2. 
                        FP_12(K,3)=-SN(J,K)*PSR(3,K)*AJ(K)/(RP(K)**3)*2. 
                        
                        FP_02(K,1)=-SN0(J,K)*AJ0(K)*PSR0(1,K)/(RP0(K)**3)*2.
                        FP_02(K,2)=-SN0(J,K)*AJ0(K)*PSR0(2,K)/(RP0(K)**3)*2.
                        FP_02(K,3)=-SN0(J,K)*AJ0(K)*PSR0(3,K)/(RP0(K)**3)*2.
                        
                        HY(J,1)=HY(J,1)+GW(KGUS+IY-1)*(FP_12(K,1)-FP_02(K,1))/P1                    
                        HY(J,2)=HY(J,2)+GW(KGUS+IY-1)*(FP_12(K,2)-FP_02(K,2))/P1                    
                        HY(J,3)=HY(J,3)+GW(KGUS+IY-1)*(FP_12(K,3)-FP_02(K,3))/P1                    
                        !END IF
                        IF(IX.EQ.1) THEN    
                            IN2(1)=IN2(1)+GW(KGUS+IY-1)*FP_02(K,1)*LOG(RP0(K))
                            IN2(2)=IN2(2)+GW(KGUS+IY-1)*FP_02(K,2)*LOG(RP0(K))
                            IN2(3)=IN2(3)+GW(KGUS+IY-1)*FP_02(K,3)*LOG(RP0(K))
                        END IF
                        !ELSE

                        !FP_12(K,1)=-HN(J,K)*PSR(1,K)*AJ(K)/(RP(K)**3)*2. 
                        !FP_12(K,2)=-HN(J,K)*PSR(2,K)*AJ(K)/(RP(K)**3)*2.
                        !FP_12(K,3)=-HN(J,K)*PSR(3,K)*AJ(K)/(RP(K)**3)*2.

                        !FP_02(K,1)=-HN0(J,K)*AJ0(K)*PSR0(1,K)/(RP0(K)**3)*2.
                        !FP_02(K,2)=-HN0(J,K)*AJ0(K)*PSR0(2,K)/(RP0(K)**3)*2.
                        !FP_02(K,3)=-HN0(J,K)*AJ0(K)*PSR0(3,K)/(RP0(K)**3)*2.
                        
                        !HY(J,1)=HY(J,1)+GW(KGUS+IY-1)*(FP_12(K,1)-FP_02(K,1))                   
                        !HY(J,2)=HY(J,2)+GW(KGUS+IY-1)*(FP_12(K,2)-FP_02(K,2))                   
                        !HY(J,3)=HY(J,3)+GW(KGUS+IY-1)*(FP_12(K,3)-FP_02(K,3))  

                        !IF(IX.EQ.1) THEN    
                        !    IN2(1)=IN2(1)+GW(KGUS+IY-1)*LOG(RP0(K))*(-SN0(J,K)*AJ0(K)*PSR0(1,K)/(RP0(K)**3)*2.)
                        !    IN2(2)=IN2(2)+GW(KGUS+IY-1)*LOG(RP0(K))*(-SN0(J,K)*AJ0(K)*PSR0(2,K)/(RP0(K)**3)*2.)
                        !    IN2(3)=IN2(3)+GW(KGUS+IY-1)*LOG(RP0(K))*(-SN0(J,K)*AJ0(K)*PSR0(3,K)/(RP0(K)**3)*2.)
                        !END IF  
                        !END IF
                    END IF
                END DO
            END DO

            DO J=1,9
                AIJ1(J)=AIJ1(J)+GW(KGUS+IX-1)*AY(J)
                HIJ1(J,1)=HIJ1(J,1)+GW(KGUS+IX-1)*HY(J,1)
                HIJ1(J,2)=HIJ1(J,2)+GW(KGUS+IX-1)*HY(J,2)
                HIJ1(J,3)=HIJ1(J,3)+GW(KGUS+IX-1)*HY(J,3)
            END DO
        END DO
    END DO        
    HIJ1(1,1)=HIJ1(1,1)+IN2(1)
    HIJ1(1,2)=HIJ1(1,2)+IN2(2)
    HIJ1(1,3)=HIJ1(1,3)+IN2(3)

!   ASSEMBLE THE LOCAL COEF. MATRIX
 
    DO J=1,8
        K=MOD(ISN+J-2,8)+1
        AIJ(K)=AIJ1(J)
        HIJ(K,1)=HIJ1(J,1)
        HIJ(K,2)=HIJ1(J,2)
        HIJ(K,3)=HIJ1(J,3)
    END DO      
	AIJ(9)=AIJ1(9)
	HIJ(9,1)=HIJ1(9,1)
	HIJ(9,2)=HIJ1(9,2)  
	HIJ(9,3)=HIJ1(9,3) 

!   MIDSIDE NODE 2, 4, 6, 8    
!   ***************************
    ELSE IF(ISN.NE.9)THEN 
      
    DO I=1,8
        K=MOD(ISN+I-2,8)+1
        K1=MOD(I,8)+1
        DO J=1,3
            CORD1(K1,J)=CORD(K,J)
            CONM1(K1,J)=CONM(K,J)  
        END DO

    END DO

   
	DO J=1,3
		CORD1(9,J)=CORD(9,J)
		CONM1(9,J)=CONM(9,J)
    END DO

    IN2=0.
    IN1=0.                
    DO IX=1,NGUS
        P1=0.5*(1.0+GA(KGUS+IX-1))
        P0=0.5*(1.0+GA(KGUS))
        P0=0.
        DO IY=1,NGUS
            P2=0.5*(1.0+GA(KGUS+IY-1)) 
            XI1=P1
            ETA1=-1.+2.*P1*P2
            XI2=P1-2.*P1*P2
            ETA2=-1.+2.*P1 
            XI3=-P1
            ETA3=-1.+2.*P1-2.*P1*P2
            CALL ISOPAR_9(XI1,ETA1,CORD1,SN1,SNX,SNE,AJ(1),SR1,DN1)
            CALL ISOPAR_9(XI2,ETA2,CORD1,SN2,SNX,SNE,AJ(2),SR2,DN2)
            CALL ISOPAR_9(XI3,ETA3,CORD1,SN3,SNX,SNE,AJ(3),SR3,DN3)

            CALL ISOPAR_2_1(P1,P2,CORD1,HN(:,1),PSR(:,1),RP(1))
            CALL ISOPAR_2_2(P1,P2,CORD1,HN(:,2),PSR(:,2),RP(2))
            CALL ISOPAR_2_3(P1,P2,CORD1,HN(:,3),PSR(:,3),RP(3))

            XI1=P0
            ETA1=-1.+2.*P0*P2
            XI2=P0-2.*P0*P2
            ETA2=-1.+2.*P0 
            XI3=-P0
            ETA3=-1.+2.*P0-2.*P0*P2
            CALL ISOPAR_9(XI1,ETA1,CORD1,SN0(:,1),SNX,SNE,AJ0(1),SR0(:,1),DN0(:,1))
            CALL ISOPAR_9(XI2,ETA2,CORD1,SN0(:,2),SNX,SNE,AJ0(2),SR0(:,2),DN0(:,2))
            CALL ISOPAR_9(XI3,ETA3,CORD1,SN0(:,3),SNX,SNE,AJ0(3),SR0(:,3),DN0(:,3))

            CALL ISOPAR_2_1(P0,P2,CORD1,HN0(:,1),PSR0(:,1),RP0(1))
            CALL ISOPAR_2_2(P0,P2,CORD1,HN0(:,2),PSR0(:,2),RP0(2))
            CALL ISOPAR_2_3(P0,P2,CORD1,HN0(:,3),PSR0(:,3),RP0(3))

            DO J=1,3
                DO K=1,9
                    SM(K,J)=SN(K,J)
                    IF(K.EQ.2) SM(K,J)=SM(K,J)-1.

                END DO
            END DO        

            DO J=1,3
                RH(J)=0.
                DO I=1,3
                    RH(J)=RH(J)+(SR(I,J)-PS(I))*(SR(I,J)-PS(I))
                END DO
                RH(J)=SQRT(RH(J))
	            RH3(J)=0.0-RH(J)**3
            END DO

            DO J=1,3       
                BNR(J)=0.  
                BNR0(J)=0.
                DO I=1,3
                    BNR(J)=BNR(J)+PSR(I,J)*DN(I,J)
                    BNR0(J)=BNR0(J)+PSR0(I,J)*DN0(I,J)
                    !BNR(J)=BNR(J)+(SR(I,J)-PS(I))*DN(I,J)
                END DO
            END DO

            DO J=1,9
                !IF(J.NE.2) WRITE(*,*)J,HN(J,3)*P1,SN(J,3)
                !IF(J.EQ.2) WRITE(*,*)J,HN(J,3)*P1+1,SN(J,3)          
            END DO

            DO I=1,9
                AY(I)=0.
                HY(I,:)=0.
            END DO
         
            DO J=1,9
                !WRITE(*,*)J,HN(J,3)*P1+1.,SN(J,3)
                DO K=1,3
                    !WRITE(*,*)K,PSR(1,K)*P1,SR(1,K)-PS(1)
                    CF=1.
                    IF(K.EQ.2) THEN
                        CF=1.0
                    ELSE
                        CF=0.5
                    END IF
                    AY(J)=AY(J)+GW(KGUS+IY-1)*SN(J,K)*AJ(K)/RP(K)*CF
!   LEE ORIGIN
                    IF(J.NE.2) THEN
                        W=-HN(J,K)*AJ(K)*PSR(1,K)/RP(K)**3*2.*CF   !(X-XI)/R**3
                        HY(J,1)=HY(J,1)+GW(KGUS+IY-1)*W 
                        
                        W=-HN(J,K)*AJ(K)*PSR(2,K)/RP(K)**3*2.*CF   !(Y-YI)/R**3
                        HY(J,2)=HY(J,2)+GW(KGUS+IY-1)*W 
                        
                        W=-HN(J,K)*AJ(K)*PSR(3,K)/RP(K)**3*2.*CF   !(Z-ZI)/R**3
                        HY(J,3)=HY(J,3)+GW(KGUS+IY-1)*W 
                    ELSE
                        !IF(IX.NE.1) THEN
                        FP_12(K,1)=-SN(J,K)*PSR(1,K)*AJ(K)/(RP(K)**3)*2.*CF  
                        FP_12(K,2)=-SN(J,K)*PSR(2,K)*AJ(K)/(RP(K)**3)*2.*CF  
                        FP_12(K,3)=-SN(J,K)*PSR(3,K)*AJ(K)/(RP(K)**3)*2.*CF  
                        
                        FP_02(K,1)=-SN0(J,K)*AJ0(K)*PSR0(1,K)/(RP0(K)**3)*2.*CF 
                        FP_02(K,2)=-SN0(J,K)*AJ0(K)*PSR0(2,K)/(RP0(K)**3)*2.*CF 
                        FP_02(K,3)=-SN0(J,K)*AJ0(K)*PSR0(3,K)/(RP0(K)**3)*2.*CF 
                        
                        HY(J,1)=HY(J,1)+GW(KGUS+IY-1)*(FP_12(K,1)-FP_02(K,1))/P1                    
                        HY(J,2)=HY(J,2)+GW(KGUS+IY-1)*(FP_12(K,2)-FP_02(K,2))/P1                    
                        HY(J,3)=HY(J,3)+GW(KGUS+IY-1)*(FP_12(K,3)-FP_02(K,3))/P1                    
                        !END IF
                        IF(IX.EQ.1) THEN    
                            IN2(1)=IN2(1)+GW(KGUS+IY-1)*FP_02(K,1)*LOG(RP0(K))
                            IN2(2)=IN2(2)+GW(KGUS+IY-1)*FP_02(K,2)*LOG(RP0(K))
                            IN2(3)=IN2(3)+GW(KGUS+IY-1)*FP_02(K,3)*LOG(RP0(K))
                        END IF
                        !ELSE

                        !FP_12(K,1)=-HN(J,K)*PSR(1,K)*AJ(K)/(RP(K)**3)*2.*CF  
                        !FP_12(K,2)=-HN(J,K)*PSR(2,K)*AJ(K)/(RP(K)**3)*2.*CF 
                        !FP_12(K,3)=-HN(J,K)*PSR(3,K)*AJ(K)/(RP(K)**3)*2.*CF 

                        !FP_02(K,1)=-HN0(J,K)*AJ0(K)*PSR0(1,K)/(RP0(K)**3)*2.*CF 
                        !FP_02(K,2)=-HN0(J,K)*AJ0(K)*PSR0(2,K)/(RP0(K)**3)*2.*CF 
                        !FP_02(K,3)=-HN0(J,K)*AJ0(K)*PSR0(3,K)/(RP0(K)**3)*2.*CF 
                        
                        !HY(J,1)=HY(J,1)+GW(KGUS+IY-1)*(FP_12(K,1)-FP_02(K,1))                   
                        !HY(J,2)=HY(J,2)+GW(KGUS+IY-1)*(FP_12(K,2)-FP_02(K,2))                   
                        !HY(J,3)=HY(J,3)+GW(KGUS+IY-1)*(FP_12(K,3)-FP_02(K,3))  

                        !IF(IX.EQ.1) THEN    
                        !    IN2(1)=IN2(1)+GW(KGUS+IY-1)*LOG(RP0(K))*(-SN0(J,K)*AJ0(K)*PSR0(1,K)/(RP0(K)**3)*2.)
                        !    IN2(2)=IN2(2)+GW(KGUS+IY-1)*LOG(RP0(K))*(-SN0(J,K)*AJ0(K)*PSR0(2,K)/(RP0(K)**3)*2.)
                        !    IN2(3)=IN2(3)+GW(KGUS+IY-1)*LOG(RP0(K))*(-SN0(J,K)*AJ0(K)*PSR0(3,K)/(RP0(K)**3)*2.)
                        !END IF  
                        !END IF

                    END IF
                END DO
            END DO

  !170 HY(J)=HY(J)+GW(KGUS+IY-1)*W*CF

            DO J=1,9
                AIJ1(J)=AIJ1(J)+GW(KGUS+IX-1)*AY(J)
                HIJ1(J,1)=HIJ1(J,1)+GW(KGUS+IX-1)*HY(J,1)
                HIJ1(J,2)=HIJ1(J,2)+GW(KGUS+IX-1)*HY(J,2)
                HIJ1(J,3)=HIJ1(J,3)+GW(KGUS+IX-1)*HY(J,3)
            END DO
        END DO
    END DO        
    HIJ1(2,1)=HIJ1(2,1)+IN2(1)
    HIJ1(2,2)=HIJ1(2,2)+IN2(2)
    HIJ1(2,3)=HIJ1(2,3)+IN2(3)

!  ASSEMBLE THE LOCAL COEF. MATRIX
    !WRITE(*,*)HIJ1(1)
    
    DO J=1,8
      	K=MOD(ISN+J-2,8)+1
      	K1=MOD(J,8)+1
      	AIJ(K)=AIJ1(K1)
        HIJ(K,1)=HIJ1(K1,1)
        HIJ(K,2)=HIJ1(K1,2)
        HIJ(K,3)=HIJ1(K1,3)
    END DO		
	AIJ(9)=AIJ1(9)
	HIJ(9,1)=HIJ1(9,1)
	HIJ(9,2)=HIJ1(9,2)  
	HIJ(9,3)=HIJ1(9,3) 

	ELSE
   
    DO 235 I=1,9
      DO 235 J=1,3
		CORD1(I,J)=CORD(I,J)
		CONM1(I,J)=CONM(I,J)    
 235	CONTINUE   
   

    IN2=0.
    IN1=0
    DO IX=1,NGUS
        P1=0.5*(1.0+GA(KGUS+IX-1))
        P0=0.5*(1.0+GA(KGUS))
        P0=0.
        DO IY=1,NGUS
            P2=0.5*(1.0+GA(KGUS+IY-1))
    
            XI1=-P1+2.*P1*P2
            ETA1=-P1
            XI2=P1
            ETA2=-P1+2.*P1*P2 
            XI3=P1-2.*P1*P2
            ETA3=P1
	        XI4=-P1
	        ETA4=P1-2*P1*P2
            CALL ISOPAR_9(XI1,ETA1,CORD1,SN1,SNX,SNE,AJ(1),SR1,DN1)
            CALL ISOPAR_9(XI2,ETA2,CORD1,SN2,SNX,SNE,AJ(2),SR2,DN2)
            CALL ISOPAR_9(XI3,ETA3,CORD1,SN3,SNX,SNE,AJ(3),SR3,DN3)
	        CALL ISOPAR_9(XI4,ETA4,CORD1,SN4,SNX,SNE,AJ(4),SR4,DN4)

            !DN(:,1)=DN1(:)
            !DN(:,2)=DN2(:)
            !DN(:,3)=DN3(:)
            !DN(:,4)=DN4(:)

            CALL ISOPAR_3_1(P1,P2,CORD1,HN(:,1),PSR(:,1),RP(1))
            CALL ISOPAR_3_2(P1,P2,CORD1,HN(:,2),PSR(:,2),RP(2))
            CALL ISOPAR_3_3(P1,P2,CORD1,HN(:,3),PSR(:,3),RP(3)) 
            CALL ISOPAR_3_4(P1,P2,CORD1,HN(:,4),PSR(:,4),RP(4))  

            XI1=-P0+2.*P0*P2
            ETA1=-P0
            XI2=P0
            ETA2=-P0+2.*P0*P2 
            XI3=P0-2.*P0*P2
            ETA3=P0
	        XI4=-P0
	        ETA4=P0-2*P0*P2
            CALL ISOPAR_9(XI1,ETA1,CORD1,SN0(:,1),SNX,SNE,AJ0(1),SR0(:,1),DN0(:,1))
            CALL ISOPAR_9(XI2,ETA2,CORD1,SN0(:,2),SNX,SNE,AJ0(2),SR0(:,2),DN0(:,2))
            CALL ISOPAR_9(XI3,ETA3,CORD1,SN0(:,3),SNX,SNE,AJ0(3),SR0(:,3),DN0(:,3))
            CALL ISOPAR_9(XI4,ETA4,CORD1,SN0(:,4),SNX,SNE,AJ0(4),SR0(:,4),DN0(:,4))

            CALL ISOPAR_3_1(P0,P2,CORD1,HN0(:,1),PSR0(:,1),RP0(1))
            CALL ISOPAR_3_2(P0,P2,CORD1,HN0(:,2),PSR0(:,2),RP0(2))
            CALL ISOPAR_3_3(P0,P2,CORD1,HN0(:,3),PSR0(:,3),RP0(3))
            CALL ISOPAR_3_4(P0,P2,CORD1,HN0(:,4),PSR0(:,4),RP0(4))

            DO J=1,9
                !IF(J.NE.9) WRITE(*,*)J,HN(J,4)*P1,SN(J,4)
                !IF(J.EQ.9) WRITE(*,*)J,HN(J,4)*P1+1,SN(J,4)          
            END DO

            IF(IX.EQ.1) THEN
                !SN0=SN
                !HN0=HN
                !DN0=DN
                !RP0=RP
                !PSR0=PSR
                !AJ0=AJ
            END IF

            DO 242 J=1,4
                DO 242 K=1,9
                    SM(K,J)=SN(K,J)
                    IF(K.EQ.9) SM(K,J)=SM(K,J)-1.
                    !IF(K.EQ.9) TSN(K,J)=TSN(K,J)-1.
        242 CONTINUE      

            DO 245 J=1,4
                RH(J)=0.0
                DO 246 I=1,3
                246 RH(J)=RH(J)+(SR(I,J)-PS(I))*(SR(I,J)-PS(I))
                    RH(J)=SQRT(RH(J))
	                RH3(J)=0.0-RH(J)**3
        245	CONTINUE


            DO J=1,4      
                BNR(J)=0.
                BNR0(J)=0.
                DO I=1,3
            !255    BNR(J)=BNR(J)+(SR(I,J)-PS(I))*DN(I,J)
                    BNR(J)=BNR(J)+PSR(I,J)*DN(I,J)
                    BNR0(J)=BNR0(J)+PSR0(I,J)*DN0(I,J)
                END DO
            END DO
     
            DO I=1,9
                AY(I)=0.
                HY(I,:)=0.
            END DO
      
            IN1=0
            DO J=1,9
                DO K=1,4
                    !WRITE(*,*)K,PSR(3,K)*P1,SR(3,K)-PS(3)
	                CF=0.5
                    AY(J)=AY(J)+GW(KGUS+IY-1)*SN(J,K)*AJ(K)/RP(K)*CF 
!   LEE ORIGIN
                    IF(J.NE.9) THEN
                        W=-HN(J,K)*AJ(K)*PSR(1,K)/RP(K)**3*2.*CF   !(X-XI)/R**3
                        HY(J,1)=HY(J,1)+GW(KGUS+IY-1)*W 
                        
                        W=-HN(J,K)*AJ(K)*PSR(2,K)/RP(K)**3*2.*CF   !(Y-YI)/R**3
                        HY(J,2)=HY(J,2)+GW(KGUS+IY-1)*W 
                        
                        W=-HN(J,K)*AJ(K)*PSR(3,K)/RP(K)**3*2.*CF   !(Z-ZI)/R**3
                        HY(J,3)=HY(J,3)+GW(KGUS+IY-1)*W 
                    ELSE
                        !IF(IX.NE.1) THEN
                        FP_12(K,1)=-SN(J,K)*PSR(1,K)*AJ(K)/(RP(K)**3)*2.*CF  
                        FP_12(K,2)=-SN(J,K)*PSR(2,K)*AJ(K)/(RP(K)**3)*2.*CF  
                        FP_12(K,3)=-SN(J,K)*PSR(3,K)*AJ(K)/(RP(K)**3)*2.*CF  
                        
                        FP_02(K,1)=-SN0(J,K)*AJ0(K)*PSR0(1,K)/(RP0(K)**3)*2.*CF 
                        FP_02(K,2)=-SN0(J,K)*AJ0(K)*PSR0(2,K)/(RP0(K)**3)*2.*CF 
                        FP_02(K,3)=-SN0(J,K)*AJ0(K)*PSR0(3,K)/(RP0(K)**3)*2.*CF 
                        
                        HY(J,1)=HY(J,1)+GW(KGUS+IY-1)*(FP_12(K,1)-FP_02(K,1))/P1                    
                        HY(J,2)=HY(J,2)+GW(KGUS+IY-1)*(FP_12(K,2)-FP_02(K,2))/P1                    
                        HY(J,3)=HY(J,3)+GW(KGUS+IY-1)*(FP_12(K,3)-FP_02(K,3))/P1                    
                        
                        IF(IX.EQ.1) THEN    
                            IN2(1)=IN2(1)+GW(KGUS+IY-1)*FP_02(K,1)*LOG(RP0(K))
                            IN2(2)=IN2(2)+GW(KGUS+IY-1)*FP_02(K,2)*LOG(RP0(K))
                            IN2(3)=IN2(3)+GW(KGUS+IY-1)*FP_02(K,3)*LOG(RP0(K))
                        END IF
                        !ELSE

                        !FP_12(K,1)=-HN(J,K)*PSR(1,K)*AJ(K)/(RP(K)**3)*2.*CF  
                        !FP_12(K,2)=-HN(J,K)*PSR(2,K)*AJ(K)/(RP(K)**3)*2.*CF 
                        !FP_12(K,3)=-HN(J,K)*PSR(3,K)*AJ(K)/(RP(K)**3)*2.*CF 

                        !FP_02(K,1)=-HN0(J,K)*AJ0(K)*PSR0(1,K)/(RP0(K)**3)*2.*CF 
                        !FP_02(K,2)=-HN0(J,K)*AJ0(K)*PSR0(2,K)/(RP0(K)**3)*2.*CF 
                        !FP_02(K,3)=-HN0(J,K)*AJ0(K)*PSR0(3,K)/(RP0(K)**3)*2.*CF 
                        
                        !HY(J,1)=HY(J,1)+GW(KGUS+IY-1)*(FP_12(K,1)-FP_02(K,1))                   
                        !HY(J,2)=HY(J,2)+GW(KGUS+IY-1)*(FP_12(K,2)-FP_02(K,2))                   
                        !HY(J,3)=HY(J,3)+GW(KGUS+IY-1)*(FP_12(K,3)-FP_02(K,3))  

                        !IF(IX.EQ.1) THEN    
                        !    IN2(1)=IN2(1)+GW(KGUS+IY-1)*LOG(RP0(K))*(-SN0(J,K)*AJ0(K)*PSR0(1,K)/(RP0(K)**3)*2.)
                        !    IN2(2)=IN2(2)+GW(KGUS+IY-1)*LOG(RP0(K))*(-SN0(J,K)*AJ0(K)*PSR0(2,K)/(RP0(K)**3)*2.)
                        !    IN2(3)=IN2(3)+GW(KGUS+IY-1)*LOG(RP0(K))*(-SN0(J,K)*AJ0(K)*PSR0(3,K)/(RP0(K)**3)*2.)
                        !END IF  
                        !END IF

                    END IF

                END DO
            END DO
!C
            DO J=1,9
                AIJ1(J)=AIJ1(J)+GW(KGUS+IX-1)*AY(J)
                HIJ1(J,1)=HIJ1(J,1)+GW(KGUS+IX-1)*HY(J,1)
                HIJ1(J,2)=HIJ1(J,2)+GW(KGUS+IX-1)*HY(J,2)
                HIJ1(J,3)=HIJ1(J,3)+GW(KGUS+IX-1)*HY(J,3)
            END DO
        END DO
    END DO        
    HIJ1(9,1)=HIJ1(9,1)+IN2(1)
    HIJ1(9,2)=HIJ1(9,2)+IN2(2)
    HIJ1(9,3)=HIJ1(9,3)+IN2(3)

      DO 280 J=1,9
        AIJ(J)=AIJ1(J)
        HIJ(J,1)=HIJ1(J,1)
        HIJ(J,2)=HIJ1(J,2)
        HIJ(J,3)=HIJ1(J,3)        
280	CONTINUE	
  
      ENDIF
      
      RETURN
END SUBROUTINE


SUBROUTINE ISOPAR_1_1(P1,P2,CORD,SN,SR,RP)
!    *************************************************
!    THE ISOPARAMETRIC ELEMENT WITH 8 NODES: CORD(8,3)
!    CALCULATE THE VALUE IN THE LOCAL POINT (XI,ETA)
!   (1) THE SHAPE FUNCTION  SN(8)
!   (2) DERIVATIVES OF SHAPE FUNCTION SNX(8),SNE(8)
!   (3) JACOBIAN MATRIX AJ
!   (4) POSITION VECTOR: SR(3)
!   (5) NORMAL VECTOR: DN(3)
!   *************************************************
!   IMPLICIT REAL (A-H,O-Z)      
    
    REAL:: CORD(9,3),SN(9),SNX(9),SNE(9),SR(3),DN(3)
    REAL:: PLX(9),PLE(9),XXE(2,3),CC(3)
    REAL:: AJ
    REAL:: A,B,C,D
    DATA PLX /-1.,0.,1.,1.,1.,0.,-1.,-1.,0/
    DATA PLE /-1.,-1.,-1.,0.,1.,1.,1.,0.,0/
    REAL:: P1,P2,RP

    !A=P1;B=2.*P1;C=P1*P2;D=2.*P1*P2
    
    SN(1)=(-1.+2.*P1)*(-1.+P1*P2)*(-1.+2.*P1*P2)&
          -2.*(-1.+P1*P2)*(-1.+2.*P1*P2)&
          +P2*(-1.+2.*P1*P2)&
          -2.*P2
    !SN(1)=(-1.+P1)*(-1.+2.*P1)*(-1.+P1*P2)*(-1.+2.*P1*P2)
    SN(2)=-4*(-1.+P1)*(-1.+P1*P2)*(-1.+2.*P1*P2)
    SN(3)=(-1.+2.*P1)*(-1.+P1*P2)*(-1.+2.*P1*P2)
    SN(4)=-4*P1*(-1.+2*P1)*P2*(-1.+P1*P2)
    SN(5)=P1*(-1.+2.*P1)*P2*(-1.+2.*P1*P2)
    SN(6)=-4.*(-1.+P1)*P1*P2*(-1.+2.*P1*P2)
    SN(7)=(-1.+P1)*(-1.+2.*P1)*P2*(-1.+2.*P1*P2)
    SN(8)=-4*(-1.+P1)*(-1.+2.*P1)*P2*(-1.+P1*P2)
    SN(9)=16*(-1.+P1)*P1*P2*(-1.+P1*P2)

    SR=0.;RP=0.
    DO I=1,9
        DO J=1,3
            SR(J)=SR(J)-SN(I)*CORD(I,J)
        END DO
    END DO
    RP=SQRT(SR(1)**2+SR(2)**2+SR(3)**2)

END SUBROUTINE
 

SUBROUTINE ISOPAR_1_2(P1,P2,CORD,SN,SR,RP)
!    *************************************************
!    THE ISOPARAMETRIC ELEMENT WITH 8 NODES: CORD(8,3)
!    CALCULATE THE VALUE IN THE LOCAL POINT (XI,ETA)
!   (1) THE SHAPE FUNCTION  SN(8)
!   (2) DERIVATIVES OF SHAPE FUNCTION SNX(8),SNE(8)
!   (3) JACOBIAN MATRIX AJ
!   (4) POSITION VECTOR: SR(3)
!   (5) NORMAL VECTOR: DN(3)
!   *************************************************
!   IMPLICIT REAL (A-H,O-Z)      
    
    REAL:: CORD(9,3),SN(9),SNX(9),SNE(9),SR(3),DN(3)
    REAL:: PLX(9),PLE(9),XXE(2,3),CC(3)
    REAL:: AJ
    DATA PLX /-1.,0.,1.,1.,1.,0.,-1.,-1.,0/
    DATA PLE /-1.,-1.,-1.,0.,1.,1.,1.,0.,0/
    REAL:: P1,P2,RP

    SN(1)=(-1.+2.*P1)*(1.+P1*(-1.+P2))*(1.+2.*P1*(-1.+P2))&
          -2.*(1.+P1*(-1.+P2))*(1.+2.*P1*(-1.+P2))&
          +(-1.+P2)*(1.+2.*P1*(-1.+P2))&
          +2.*(-1.+P2)
    !SN(1)=(-1.+P1)*(-1.+2.*P1)*(1.+P1*(-1.+P2))*(1.+2.*P1*(-1.+P2))
    SN(2)=-4.*(-1.+P1)*(-1.+2.*P1)*(1.+P1*(-1.+P2))*(-1.+P2)
    SN(3)=(-1.+P1)*(-1+2*P1)*(1+2*P1*(-1+P2))*(-1+P2)
    SN(4)=-4.*(-1.+P1)*P1*(1+2*P1*(-1.+P2))*(-1.+P2)
    SN(5)=P1*(-1.+2*P1)*(1.+2*P1*(-1.+P2))*(-1.+P2)
    SN(6)=-4*P1*(-1.+2*P1)*(1.+P1*(-1.+P2))*(-1.+P2)
    SN(7)=(-1.+2*P1)*(1+P1*(-1.+P2))*(1.+2.*P1*(-1.+P2))
    SN(8)=-4*(-1.+P1)*(1+P1*(-1.+P2))*(1.+2.*P1*(-1.+P2))
    SN(9)=16*(-1.+P1)*P1*(1+P1*(-1.+P2))*(-1.+P2)

    SR=0.;RP=0.
    DO I=1,9
        DO J=1,3
            SR(J)=SR(J)-SN(I)*CORD(I,J)
        END DO
    END DO
    RP=SQRT(SR(1)**2+SR(2)**2+SR(3)**2)
END SUBROUTINE 


SUBROUTINE ISOPAR_2_1(P1,P2,CORD,SN,SR,RP)
!    *************************************************
!    THE ISOPARAMETRIC ELEMENT WITH 8 NODES: CORD(8,3)
!    CALCULATE THE VALUE IN THE LOCAL POINT (XI,ETA)
!   (1) THE SHAPE FUNCTION  SN(8)
!   (2) DERIVATIVES OF SHAPE FUNCTION SNX(8),SNE(8)
!   (3) JACOBIAN MATRIX AJ
!   (4) POSITION VECTOR: SR(3)
!   (5) NORMAL VECTOR: DN(3)
!   *************************************************
!   IMPLICIT REAL (A-H,O-Z)      
    
    REAL:: CORD(9,3),SN(9),SNX(9),SNE(9),SR(3),DN(3)
    REAL:: PLX(9),PLE(9),XXE(2,3),CC(3)
    REAL:: AJ,P1,P2,RP
    DATA PLX /-1.,0.,1.,1.,1.,0.,-1.,-1.,0/
    DATA PLE /-1.,-1.,-1.,0.,1.,1.,1.,0.,0/ 

    !1/4 (1 - p1) p1 (2 - 2 p1 p2) (-1 + 2 p1 p2)
    SN(1)=1./2.*(-1.+P1)*(-1.+P1*P2)*(-1.+2.*P1*P2)
    !-(1/2) (1 - p1^2) (2 - 2 p1 p2) (-1 + 2 p1 p2)
    SN(2)=-P1*(-1.+P1*P2)*(-1.+2.*P1*P2)&
          +P2*(-1.+2.*P1*P2)&
          -2.*P2
    !-1/4 p1 (1 + p1) (2 - 2 p1 p2) (-1 + 2 p1 p2)
    SN(3)=(1./2.)*(1.+P1)*(-1.+P1*P2)*(-1. + 2. *P1 *P2)
    !-2 p1^2 (1 + p1) p2 (-1 + p1 p2)
    SN(4)=-2.*P1*(1.+P1)*P2*(-1.+P1*P2)
    !1/2 p1^2 (1 + p1) p2 (-1 + 2 p1 p2)
    SN(5)=1./2.*P1*(1.+P1)*P2*(-1.+2*P1*P2)
    !p1 (1 - p1^2) p2 (-1 + 2 p1 p2)
    SN(6)=-(-1.+P1**2)*P2*(-1.+2*P1*P2)
    !-(1/2) (1 - p1) p1^2 p2 (-1 + 2 p1 p2)
    SN(7)=(1./2.)*(-1.+P1)*P1*P2*(-1.+2.*P1*P2)
    !-2 (-1 + p1) p1^2 p2 (-1 + p1 p2)
    SN(8)=-2.*(-1.+P1)*P1*P2*(-1.+P1*P2)
    !4 p1 (-1 + p1^2) p2 (-1 + p1 p2)
    SN(9)=4.*(-1.+P1**2)*P2*(-1.+P1*P2)

    SR=0.;RP=0.
    DO I=1,9
        DO J=1,3
            SR(J)=SR(J)-SN(I)*CORD(I,J)
        END DO
    END DO
    RP=SQRT(SR(1)**2+SR(2)**2+SR(3)**2)
END SUBROUTINE 

SUBROUTINE ISOPAR_2_2(P1,P2,CORD,SN,SR,RP)
!    *************************************************
!    THE ISOPARAMETRIC ELEMENT WITH 8 NODES: CORD(8,3)
!    CALCULATE THE VALUE IN THE LOCAL POINT (XI,ETA)
!   (1) THE SHAPE FUNCTION  SN(8)
!   (2) DERIVATIVES OF SHAPE FUNCTION SNX(8),SNE(8)
!   (3) JACOBIAN MATRIX AJ
!   (4) POSITION VECTOR: SR(3)
!   (5) NORMAL VECTOR: DN(3)
!   *************************************************
!   IMPLICIT REAL (A-H,O-Z)      
    
    REAL:: CORD(9,3),SN(9),SNX(9),SNE(9),SR(3),DN(3)
    REAL:: PLX(9),PLE(9),XXE(2,3),CC(3)
    REAL:: AJ,P1,P2,RP
    DATA PLX /-1.,0.,1.,1.,1.,0.,-1.,-1.,0/
    DATA PLE /-1.,-1.,-1.,0.,1.,1.,1.,0.,0/

    SN(1)=1./2. *(-1. + p1)*(-1.+ 2.* p1)* (-1. + 2.* p2)* (1. + p1* (-1. + 2.* p2))
    !SN(2)=(-1. + 2.* P1) *(1. - (P1 - 2. *P1 *P2)**2)&
    !      - 2. *(1. - (P1 - 2. *P1 *P2)**2)&
    !      - (1. - 2. *P2)*(P1 - 2. *P1 *P2)

    SN(2)=(-1.+2.*P1)*(1.-(P1-2.*P1*P2)**2)&
          -2.*(1.-(P1-2.*P1*P2)**2)&
          -(1.-2.*P2)*(P1-2.*P1*P2)
    !SN(2)=-(1.-P1)*(-1.+2*P1)*(1.-(P1-2.*P1*P2)**2)
    !-(1/4) (2. - 2. P1) (-1. + 2. P1) (P1 - 2. P1 P2) (1 + P1 - 2. P1 P2)
    SN(3)=-(1./4.) *(2. - 2.* p1)* (-1. + 2.* p1)* (1. - 2.*p2) *(1. + p1 - 2.* p1* p2)
    !SN(3)=-(1./2.)*(1.-P1)*(-1+2.*P1)*(1.-2.*P2)*(1+P1-2*P1*P2)
    !1/2 (1 - (-1. + 2. P1)^2) (P1 - 2. P1 P2) (1 + P1 - 2. P1 P2)
    SN(4)=-2.*(-1.+P1)*P1*(-1.+2.*P2)*(-1.+P1*(-1.+2.*P2))
    SN(5)=1./2. *(-1. + 2. *P1)* (P1 - 2. *P1 *P2)* (1. + P1 - 2.* P1* P2)
    !SN(6)=1./2. *(-1. + 2.* P1)*  2.* (1. - (P1 - 2.* P1 *P2)**2)
    SN(6)=(-1.+2.*P1)*(1.-(P1-2.*P1*P2)**2)
    !SN(7)=-(1./4.)* (-1. + 2. *P1)* 2.*(P1 - 2.* P1 *P2) *(1. - P1 + 2. *P1 *P2)
    SN(7)=-(1./2.)*(-1.+2.*P1)*(P1-2.*P1*P2)*(1.+P1*(-1.+2.*P2))
    !SN(8)=-(1./2.) *(1. - (-1. + 2.* P1)**2)* (P1 - 2.* P1 *P2)* (1. - P1+ 2.* P1* P2)
    SN(8)=-2.*(-1.+P1)*P1*(-1.+2.*P2)*(1.+P1*(-1.+2.*P2))
    SN(9)=-4.*(-1.+p1)*(1.-(P1-2.*p1*p2)**2)
    !SN(9)= (4. - 4. *P1) *(1.-1.* (P1 - 2. *P1 *P2)**2)

    SR=0.;RP=0.
    DO I=1,9
        DO J=1,3
            SR(J)=SR(J)-SN(I)*CORD(I,J)
        END DO
    END DO
    RP=SQRT(SR(1)**2+SR(2)**2+SR(3)**2)

END SUBROUTINE 

SUBROUTINE ISOPAR_2_3(P1,P2,CORD,SN,SR,RP)
!    *************************************************
!    THE ISOPARAMETRIC ELEMENT WITH 8 NODES: CORD(8,3)
!    CALCULATE THE VALUE IN THE LOCAL POINT (XI,ETA)
!   (1) THE SHAPE FUNCTION  SN(8)
!   (2) DERIVATIVES OF SHAPE FUNCTION SNX(8),SNE(8)
!   (3) JACOBIAN MATRIX AJ
!   (4) POSITION VECTOR: SR(3)
!   (5) NORMAL VECTOR: DN(3)
!   *************************************************
!   IMPLICIT REAL (A-H,O-Z)      
    
    REAL:: CORD(9,3),SN(9),SNX(9),SNE(9),SR(3),DN(3)
    REAL:: PLX(9),PLE(9),XXE(2,3),CC(3)
    REAL:: AJ,P1,P2,RP
    DATA PLX /-1.,0.,1.,1.,1.,0.,-1.,-1.,0/
    DATA PLE /-1.,-1.,-1.,0.,1.,1.,1.,0.,0/

    SN(1)=1./2.*(1.+P1)*(1.+P1*(-1.+P2))*(1.+2*P1*(-1.+P2))  
    SN(2)=-P1*(1.+P1*(-1+P2))*(1.+2*P1*(-1+P2))&
          +(-1.+P2)*(1.+2.*P1*(-1.+P2))&
          +2.*(-1.+P2)
    
    !SN(2)=(1.-P1**2)*(1.+P1*(-1+P2))*(1.+2*P1*(-1+P2))
    SN(3)=1./2.*(-1.+P1)*(1.+P1*(-1.+P2))*(1.+2.*P1*(-1.+P2))
    SN(4)=-2.*(-1.+P1)*P1*(1.+P1*(-1.+P2))*(-1.+P2)
    SN(5)=1./2.*(-1.+P1)*P1*(1+2.*P1*(-1+P2))*(-1.+P2)
    SN(6)=-(-1.+P1**2)*(1.+2*P1*(-1+P2))*(-1.+P2)
    SN(7)=1./2.*P1*(1.+P1)*(1+2.*P1*(-1.+P2))*(-1.+P2)
    SN(8)=-2.*P1*(1.+P1)*(1.+P1*(-1.+P2))*(-1.+P2)
    SN(9)=4.*(-1.+P1**2)*(1.+P1*(-1.+P2))*(-1.+P2)

    SR=0.;RP=0.
    DO I=1,9
        DO J=1,3
            SR(J)=SR(J)-SN(I)*CORD(I,J)
        END DO
    END DO
    RP=SQRT(SR(1)**2+SR(2)**2+SR(3)**2)

END SUBROUTINE    
      
SUBROUTINE ISOPAR_3_1(P1,P2,CORD,SN,SR,RP)
!    *************************************************
!    THE ISOPARAMETRIC ELEMENT WITH 8 NODES: CORD(8,3)
!    CALCULATE THE VALUE IN THE LOCAL POINT (XI,ETA)
!   (1) THE SHAPE FUNCTION  SN(8)
!   (2) DERIVATIVES OF SHAPE FUNCTION SNX(8),SNE(8)
!   (3) JACOBIAN MATRIX AJ
!   (4) POSITION VECTOR: SR(3)
!   (5) NORMAL VECTOR: DN(3)
!   *************************************************
!   IMPLICIT REAL (A-H,O-Z)      
    
    REAL:: CORD(9,3),SN(9),SNX(9),SNE(9),SR(3),DN(3)
    REAL:: PLX(9),PLE(9),XXE(2,3),CC(3)
    REAL:: AJ,P1,P2,RP
    DATA PLX /-1.,0.,1.,1.,1.,0.,-1.,-1.,0/
    DATA PLE /-1.,-1.,-1.,0.,1.,1.,1.,0.,0/

    SN(1)=-1./4*P1*(1+P1)*(-1+2*P2)*(1+P1-2*P1*P2)
    SN(2)=1./2*(1+P1)*(1-P1**2*(1-2*P2)**2)
    SN(3)=1./4*P1*(1+P1)*(-1+2*P2)*(1+P1*(-1+2*P2))
    SN(4)=-1./2*(-1+P1**2)*(-1+2*P2)*(1+P1*(-1+2*P2))
    SN(5)=1./4*(-1+P1)*P1*(-1+2*P2)*(1+P1*(-1+2*P2))
    SN(6)=-1./2*(-1+P1)*(-1+P1**2*(1-2*P2)**2)
    SN(7)=1./4*(1-P1)*P1*(-1+2*P2)*(1+P1-2*P1*P2)
    SN(8)=-1./2*(-1+P1**2)*(-1+2*P2)*(-1+P1*(-1+2*P2))
    
    SN(9)=-P1*(1-P1**2*(1-2*P2)**2)&
          -P1*(1-2*P2)**2
    
    !SN(9)=(1-P1**2)*(1-P1**2*(1-2*P2)**2)

    SR=0.;RP=0.
    DO I=1,9
        DO J=1,3
            SR(J)=SR(J)-SN(I)*CORD(I,J)
        END DO
    END DO
    RP=SQRT(SR(1)**2+SR(2)**2+SR(3)**2)

END SUBROUTINE       

SUBROUTINE ISOPAR_3_2(P1,P2,CORD,SN,SR,RP)
!    *************************************************
!    THE ISOPARAMETRIC ELEMENT WITH 8 NODES: CORD(8,3)
!    CALCULATE THE VALUE IN THE LOCAL POINT (XI,ETA)
!   (1) THE SHAPE FUNCTION  SN(8)
!   (2) DERIVATIVES OF SHAPE FUNCTION SNX(8),SNE(8)
!   (3) JACOBIAN MATRIX AJ
!   (4) POSITION VECTOR: SR(3)
!   (5) NORMAL VECTOR: DN(3)
!   *************************************************
!   IMPLICIT REAL (A-H,O-Z)      
    
    REAL:: CORD(9,3),SN(9),SNX(9),SNE(9),SR(3),DN(3)
    REAL:: PLX(9),PLE(9),XXE(2,3),CC(3)
    REAL:: AJ,P1,P2,RP
    DATA PLX /-1.,0.,1.,1.,1.,0.,-1.,-1.,0/
    DATA PLE /-1.,-1.,-1.,0.,1.,1.,1.,0.,0/

    SN(1)=1./4*(1-P1)*P1*(-1+2*P2)*(1+P1-2*P1*P2)
    SN(2)=-(1./2)*(-1+P1**2)*(-1+2*P2)*(-1+P1*(-1+2*P2))
    SN(3)=-1./4*P1*(1+P1)*(-1+2*P2)*(1+P1-2*P1*P2)
    SN(4)=1./2*(1+P1)*(1-P1**2*(1-2*P2)**2)
    SN(5)=1./4*P1*(1+P1)*(-1+2*P2)*(1+P1*(-1+2*P2))
    SN(6)=-(1./2)*(-1+P1**2)*(-1+2*P2)*(1+P1*(-1+2*P2))
    SN(7)=(1./4)*(-1+P1)*P1*(-1+2*P2)*(1+P1*(-1+2*P2))
    SN(8)=-(1./2)*(-1+P1)*(-1+P1**2*(1-2*P2)**2)
    
    SN(9)=-P1*(1-P1**2*(1-2*P2)**2)&
          -P1*(1-2*P2)**2
    
    !SN(9)=(1-P1**2)*(1-P1**2*(1-2*P2)**2)

    SR=0.;RP=0.
    DO I=1,9
        DO J=1,3
            SR(J)=SR(J)-SN(I)*CORD(I,J)
        END DO
    END DO
    RP=SQRT(SR(1)**2+SR(2)**2+SR(3)**2)

END SUBROUTINE         

SUBROUTINE ISOPAR_3_3(P1,P2,CORD,SN,SR,RP)
!    *************************************************
!    THE ISOPARAMETRIC ELEMENT WITH 8 NODES: CORD(8,3)
!    CALCULATE THE VALUE IN THE LOCAL POINT (XI,ETA)
!   (1) THE SHAPE FUNCTION  SN(8)
!   (2) DERIVATIVES OF SHAPE FUNCTION SNX(8),SNE(8)
!   (3) JACOBIAN MATRIX AJ
!   (4) POSITION VECTOR: SR(3)
!   (5) NORMAL VECTOR: DN(3)
!   *************************************************
!   IMPLICIT REAL (A-H,O-Z)      
    
    REAL:: CORD(9,3),SN(9),SNX(9),SNE(9),SR(3),DN(3)
    REAL:: PLX(9),PLE(9),XXE(2,3),CC(3)
    REAL:: AJ,P1,P2,RP
    DATA PLX /-1.,0.,1.,1.,1.,0.,-1.,-1.,0/
    DATA PLE /-1.,-1.,-1.,0.,1.,1.,1.,0.,0/

    SN(1)=1./4*(-1+P1)*P1*(-1+2*P2)*(1+P1*(-1+2*P2))
    SN(2)=-(1./2)*(-1+P1)*(-1+P1**2*(1-2*P2)**2)
    SN(3)=-1./4*(1-P1)*(P1-2*P1*P2)*(1+P1-2*P1*P2)
    SN(4)=1./2*(1-P1**2)*(1-2*P2)*(1+P1-2*P1*P2)
    SN(5)=1./4*(1+P1)*(P1-2*P1*P2)*(1+P1-2*P1*P2)
    SN(6)=(1./2)*(1+P1)*(1-(P1-2*P1*P2)**2)
    SN(7)=-(1./4)*(1+P1)*(P1-2*P1*P2)*(1+P1*(-1+2*P2))
    SN(8)=-(1./2)*(-1+P1**2)*(-1+2*P2)*(1+P1*(-1+2*P2))
    
    SN(9)=-P1*(1-(P1-2*P1*P2)**2)&
          -(1-2*P2)*(P1-2*P1*P2)
    
    !SN(9)=(1-P1**2)*(1-(P1-2*P1*P2)**2)

    SR=0.;RP=0.
    DO I=1,9
        DO J=1,3
            SR(J)=SR(J)-SN(I)*CORD(I,J)
        END DO
    END DO
    RP=SQRT(SR(1)**2+SR(2)**2+SR(3)**2)

END SUBROUTINE   
 
SUBROUTINE ISOPAR_3_4(P1,P2,CORD,SN,SR,RP)
!    *************************************************
!    THE ISOPARAMETRIC ELEMENT WITH 8 NODES: CORD(8,3)
!    CALCULATE THE VALUE IN THE LOCAL POINT (XI,ETA)
!   (1) THE SHAPE FUNCTION  SN(8)
!   (2) DERIVATIVES OF SHAPE FUNCTION SNX(8),SNE(8)
!   (3) JACOBIAN MATRIX AJ
!   (4) POSITION VECTOR: SR(3)
!   (5) NORMAL VECTOR: DN(3)
!   *************************************************
!   IMPLICIT REAL (A-H,O-Z)      
    
    REAL:: CORD(9,3),SN(9),SNX(9),SNE(9),SR(3),DN(3)
    REAL:: PLX(9),PLE(9),XXE(2,3),CC(3)
    REAL:: AJ,P1,P2,RP
    DATA PLX /-1.,0.,1.,1.,1.,0.,-1.,-1.,0/
    DATA PLE /-1.,-1.,-1.,0.,1.,1.,1.,0.,0/

    SN(1)=-1./4*(1+P1)*(P1-2*P1*P2)*(1+P1*(-1+2*P2))
    SN(2)=-(1./2)*(-1+P1**2)*(-1+2*P2)*(1+P1*(-1+2*P2))
    SN(3)=(1./4)*P1*(-1+P1)*(-1+2*P2)*(1+P1*(-1+2*P2))
    SN(4)=-1./2*(-1+P1)*(-1+P1**2*(1-2*P2)**2)
    SN(5)=-1./4*(1-P1)*(P1-2*P1*P2)*(1+P1-2*P1*P2)
    SN(6)=1./2*(1-P1**2)*(1-2*P2)*(1+P1-2*P1*P2)
    SN(7)=(1./4)*(1+P1)*(P1-2*P1*P2)*(1+P1-2*P1*P2)
    SN(8)=(1./2)*(1+P1)*(1-(P1-2*P1*P2)**2)
    
    SN(9)=-P1*(1.-(P1-2.*P1*P2)**2)&
          -(1.-2.*P2)*(P1-2*P1*P2)
    
    !SN(9)=(1-P1**2)*(1-(P1-2*P1*P2)**2)

    SR=0.;RP=0.
    DO I=1,9
        DO J=1,3
            SR(J)=SR(J)-SN(I)*CORD(I,J)
        END DO
    END DO
    RP=SQRT(SR(1)**2+SR(2)**2+SR(3)**2)

END SUBROUTINE 
       

