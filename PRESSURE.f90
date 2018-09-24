SUBROUTINE PRESSURE
	USE GREENMOD
	USE CUMOD
    REAL:: XI,ETA,AJ,AIWL1,AIWL
    common /gwa/ga(400),gw(400)
	DIMENSION ISM(3,2)
	DATA ISM/1,1,1,1,-1,1/

    REAL,dimension(:,:):: xxe(2,3)
    REAL,dimension(:):: sn(9),snx(9),sne(9),sr(3),dn(3),VV(3)
    REAL,dimension(:):: sn3(3),snx3(3),sne3(3),sr3(2),dn3(2),FCORL(3,2)
	REAL,dimension(:,:):: a(3,3),b(3,3)
	DIMENSION IS(3),JS(3),FRES(2),FRES0(2)
    REAL:: NX1L
	INTEGER:: L,ISYM
    CHARACTER*7 TITLE,TITLE2


    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(99,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\'//TRIM(ADJUSTL(TITLE))//'\WAVE_PRO_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
    OPEN(7,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\WAVEPRO_NONLINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')


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

	!AIJX2=0.0
    !AIJX3=0.0
    !AIJX4=0.0
    !AIJX5=0.0

    !AIJX2=-MASS0(3,3)*G
    !AIJX3=-MASS0(3,3)*G
    !AIJX4=-MASS0(5,5)*G
    !AIJX5=-MASS0(5,5)*G

    FRES=0.

    IF(ITTE.EQ.1) THEN   !LINEAR CALCULATION
        FRES0=0.
        DO IE=1,NBBLOCK
            DO ISYM=1,2
		        DO I1=1,PAN(IE)
			        I2=I1
			        IF(PAN(THR1+IE).EQ.9) THEN
                        IF(ISYM.EQ.2.AND.(I1.NE.1).AND.(I1.NE.9)) I2=10-I1
                    END IF
            
                    DO J1=1,3
				        CORL(I1,J1)=CORD0(INE(THR1+IE,I2),J1)*ISM(J1,ISYM)
				        CNL(I1,J1)=VECN0(INE(THR1+IE,I2),J1)*ISM(J1,ISYM)  
                    END DO
                END DO

		        ngus=11
		        k=ngus*(ngus-1)/2
                !对x积分     
		        do ix=1,ngus
		        xi=ga(k+ix-1)
                !对y积分
		        AY=0.0
		        do iy=1,ngus
			        eta=ga(k+iy-1)
!调用isopar，求得形状函数J及单位法向矢量等     
			        IF(PAN(IE).EQ.9) call isopar_9(xi,eta,corl,sn,snx,sne,aj,sr,dn)
			
  			        FRES0(1)=FRES0(1)+GW(k+ix-1)*GW(k+iy-1)*AJ*(-SR(3))*DN(3)
  			        FRES0(2)=FRES0(2)+GW(k+ix-1)*GW(k+iy-1)*AJ*(-SR(3))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3))
                END DO
                END DO
            END DO
        END DO
    END IF

!********************************************************
!计算波剖面
    THR1=0.
	ISYM=1
    DO IWP=1,NFX
    DO IE=1,NBBLOCK
        DO J=1,9
        IF(INE(IE,J).EQ.NWPLAZ(IWP)) THEN
            JIE=J
            GOTO 501
        END IF
        END DO
        
        GOTO 701
        501 CONTINUE

        !WRITE(*,*)JIE,IWP,NWPLAZ(IWP)
		DO I1=1,PAN(IE)
			I2=I1
			!IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1
			IF(PAN(THR1+IE).EQ.9) THEN
                IF(ISYM.EQ.2.AND.(I1.NE.1).AND.(I1.NE.9)) I2=10-I1
            END IF

            DO J1=1,3
				CORL(I1,J1)=CORD(INE(THR1+IE,I2),J1)*ISM(J1,ISYM)
				CNL(I1,J1)=VECN(INE(THR1+IE,I2),J1)*ISM(J1,ISYM)  
            END DO
        END DO

        SELECT CASE(JIE)
        CASE(1)
        XI=-1.;ETA=-1.
        CASE(2)
        XI=0.;ETA=-1.
        CASE(3)
        XI=1.;ETA=-1.
        CASE(4)
        XI=1.;ETA=0.
        CASE(5)
        XI=1.;ETA=1.
        CASE(6)
        XI=0.;ETA=1.
        CASE(7)
        XI=-1.;ETA=1.
        CASE(8)
        XI=-1.;ETA=0.
        END SELECT
        CALL isopar_9(xi,eta,corl,sn,snx,sne,aj,sr,dn)

	    do i=1,3
		    xxe(1,i)=0.0
		    xxe(2,i)=0.0
        END DO

	    do i=1,3
			do j=1,PAN(IE)
!!!C	xxe(1,*): xxi,yxi,zxi
!!!C	xxe(2,*): xeta,yeta,zeta
				xxe(1,i)=xxe(1,i)+corl(j,i)*snx(j)
   				xxe(2,i)=xxe(2,i)+corl(j,i)*sne(j)
            END DO
        END DO
!!!C
	    DO I=1,3
			DO J=1,2
				B(I,J)=XXE(J,I)
            END DO
        END DO

		DO I=1,3
				B(I,3)=DN(I)
        END DO

		DO I=1,3
			DO J=1,3
				A(I,J)=B(I,J)
!!!C	write(*,*)a(i,j)
            END DO
        END DO

		CALL BRINV(A,3)
			
		PHIKSI=0.0
		PHIETA=0.0
		DO I1=1,PAN(IE)
		    I2=I1
		    IF((ISYM.EQ.2).AND.(I1.NE.1).AND.(I1.NE.9)) I2=10-I1
			
            PHIKSI=PHIKSI+SNX(I1)*SP(INE(IE,I2))
			PHIETA=PHIETA+SNE(I1)*SP(INE(IE,I2))
        END DO
		DO I=1,3
			VV(I)=A(1,I)*PHIKSI+A(2,I)*PHIETA+A(3,I)*U*DN(1)  !VV速度势的散度
        END DO

		WPLAZ(IWP)=0.0-(0.5*(VV(1)**2+VV(2)**2&   
     			+VV(3)**2)-U*VV(1))/G
        WPLAZ(IWP)=0.0+(U*VV(1))/G	
        
        701 CONTINUE  
    END DO
    END DO

!**************************************************************************


	DO 31 IE=1,NBBLOCK
    DO 31 ISYM=1,2
!!CWRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!!C依次将每一面元的四个顶点的坐标赋值到corl(4,3)中
!!C依次将每一面元的四个顶点的法向矢量赋值到cnl(4,3)中
		
		DO I1=1,PAN(IE)
			I2=I1
			!IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1
			IF(PAN(THR1+IE).EQ.9) THEN
                IF(ISYM.EQ.2.AND.(I1.NE.1).AND.(I1.NE.9)) I2=10-I1
            END IF
            IF(PAN(IE).EQ.8) THEN
                IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
            END IF

            DO J1=1,3
				CORL(I1,J1)=CORD(INE(THR1+IE,I2),J1)*ISM(J1,ISYM)
				CNL(I1,J1)=VECN(INE(THR1+IE,I2),J1)*ISM(J1,ISYM)  
            END DO
        END DO

		ngus=11
		k=ngus*(ngus-1)/2
!gauss-legendre numerical integration

!对x积分     
		do 42 ix=1,ngus
		xi=ga(k+ix-1)
!对y积分
		AY=0.0
		do 42 iy=1,ngus
			eta=ga(k+iy-1)
!调用isopar，求得形状函数J及单位法向矢量等     
			IF(PAN(IE).EQ.9) call isopar_9(xi,eta,corl,sn,snx,sne,aj,sr,dn)
			!IF(PAN(IE).EQ.8) call isopar_8(xi,eta,corl(1:PAN(IE),:),sn(1:PAN(IE)),snx(1:PAN(IE)),sne(1:PAN(IE)),aj,sr,dn)
			do 210 i=1,3
				xxe(1,i)=0.0
				xxe(2,i)=0.0
210			continue

			do 21 i=1,3
			do 21 j=1,PAN(IE)
!!!C	xxe(1,*): xxi,yxi,zxi
!!!C	xxe(2,*): xeta,yeta,zeta
				xxe(1,i)=xxe(1,i)+corl(j,i)*snx(j)
   				xxe(2,i)=xxe(2,i)+corl(j,i)*sne(j)
   21			continue
!!!C
			DO 22 I=1,3
			DO 22 J=1,2
				B(I,J)=XXE(J,I)
22			CONTINUE

			DO 23 I=1,3
				B(I,3)=DN(I)
23			CONTINUE

			DO 	230 I=1,3
			DO  230 J=1,3
				A(I,J)=B(I,J)
!!!C	write(*,*)a(i,j)
230			CONTINUE

			CALL BRINV(A,3)
			
			PHIKSI=0.0
			PHIETA=0.0
            PHIKSX=0.
            PHIKSY=0.
            PHIKSZ=0.
			DO 24 I1=1,PAN(IE)
			    I2=I1
			    IF((ISYM.EQ.2).AND.(I1.NE.1).AND.(I1.NE.9)) I2=10-I1
				
                PHIKSI=PHIKSI+SNX(I1)*SP(INE(IE,I2))
				PHIETA=PHIETA+SNE(I1)*SP(INE(IE,I2))

                PHIKSX=PHIKSX+SN(I1)*SPX(INE(IE,I2))
                PHIKSY=PHIKSY+SN(I1)*SPY(INE(IE,I2))
                PHIKSZ=PHIKSZ+SN(I1)*SPZ(INE(IE,I2))
24			CONTINUE
			DO I=1,3
				VV(I)=A(1,I)*PHIKSI+A(2,I)*PHIETA+A(3,I)*U*DN(1)  !VV速度势的散度
            END DO
            
            !VV(1)=PHIKSX
            !VV(2)=PHIKSY
            !VV(3)=PHIKSZ

            !此处有问题

			BNR=0.0-P1*(0.5*(VV(1)**2+VV(2)**2&    !非线性阻力
     			+VV(3)**2)-U*VV(1))*DN(1)
            BNR1=0.0+P1*(U*VV(1))*DN(1)	   !线性阻力
            
            BNR2=0.0-P1*(0.5*(VV(1)**2+VV(2)**2&        !非线性垂向力
     			    +VV(3)**2)-U*VV(1))*DN(3) !-P1*G*(SR(3))*DN(3)  
            !BNR2=-P1*G*(SR(3))*DN(3)  

            BNR3=0.0+P1*(U*VV(1))*DN(3) !-P1*G*(SR(3))*DN(3) 	            !线性垂向力

            BNR4=0.0-P1*(0.5*(VV(1)*VV(1)+VV(2)*VV(2)&           !非线性纵摇力
     			 +VV(3)*VV(3))-U*VV(1))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3)) !&
                 !-P1*G*(SR(3))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3))

            BNR5=0.0+P1*(U*VV(1))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3)) !&
                 !-P1*G*(SR(3))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3))
            !BNR5=DN(2)*SR(2)

			AIJY=gw(k+iy-1)*aj*BNR
			AIJY1=GW(K+IY-1)*AJ*BNR1
            AIJY2=GW(K+IY-1)*AJ*BNR2
            AIJY3=GW(K+IY-1)*AJ*BNR3
            AIJY4=GW(K+IY-1)*AJ*BNR4
            AIJY5=GW(K+IY-1)*AJ*BNR5

  			AIJX=AIJX+gw(k+ix-1)*AIJY
			AIJX1=AIJX1+GW(K+IX-1)*AIJY1
			AIJX2=AIJX2+GW(K+IX-1)*AIJY2
			AIJX3=AIJX3+GW(K+IX-1)*AIJY3
			AIJX4=AIJX4+GW(K+IX-1)*AIJY4
			AIJX5=AIJX5+GW(K+IX-1)*AIJY5

            IF(ITTE.NE.1) THEN   !NONLINEAR CALCULATION
  			    FRES(1)=FRES(1)+GW(k+ix-1)*GW(k+iy-1)*AJ*(-SR(3))*DN(3)
  			    FRES(2)=FRES(2)+GW(k+ix-1)*GW(k+iy-1)*AJ*(-SR(3))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3))
            END IF

            !AIJX  非线性阻力
            !AIJY1 线性阻力
            !AIJY2 非线性垂向力
            !AIJY3 线性垂向力
            !AIJY4 非线性纵摇力
            !AIJY5 线性纵摇力
            IF(ISYM.EQ.1) THEN 
                IF(MOD(IE-1,(NZ-1)/2).EQ.0) THEN
                    IF(IX.EQ.1) THEN
                        !WRITE(99,*)SR(1),U/G*VV(1),&
                        !                 (-0.5*(VV(1)**2+VV(2)**2&    !非线性阻力
     			        !                  +VV(3)**2)+U*VV(1))/G
                    END IF
                END IF
            END IF

   42		continue
   31	CONTINUE

!计算水线积分
    AIWL=0.
    AIWL1=0.
    DO I=1,(NFX-1)/2
        IE=NBBLOCK+NFXF*(NFY-1)/4+(I-1)*(NFY-1)/2+1
		
        ETW1_1=SPX(INE(IE,1))*U/G
        ETW2_1=SPX(INE(IE,2))*U/G
        ETW3_1=SPX(INE(IE,3))*U/G
        ETW1=SPX(INE(IE,1))*U/G
        ETW2=SPX(INE(IE,2))*U/G
        ETW3=SPX(INE(IE,3))*U/G
        !ETW1=SPX(INE(IE,1))*U/G-0.5*(SPX(INE(IE,1))**2+SPY(INE(IE,1))**2+SPZ(INE(IE,1))**2)/G
        !ETW2=SPX(INE(IE,2))*U/G-0.5*(SPX(INE(IE,2))**2+SPY(INE(IE,2))**2+SPZ(INE(IE,2))**2)/G
        !ETW3=SPX(INE(IE,3))*U/G-0.5*(SPX(INE(IE,3))**2+SPY(INE(IE,3))**2+SPZ(INE(IE,3))**2)/G

        FCORL(1,1:2)=CORD(INE(IE,1),1:2)
        FCORL(2,1:2)=CORD(INE(IE,2),1:2)
        FCORL(3,1:2)=CORD(INE(IE,3),1:2)      

        NGUS=8
		K=NGUS*(NGUS-1)/2
!对x积分     
		DO IX=1,NGUS
		    XI=GA(K+IX-1)

            CALL isopar_3(xi,FCORL,sn3,snx3,sne3,aj,sr3,dn3)
            
            PHIKSI1=ETW1_1**2*SN3(1)+ETW2_1**2*SN3(2)+ETW3_1**2*SN3(3)
            PHIKSI=ETW1**2*SN3(1)+ETW2**2*SN3(2)+ETW3**2*SN3(3)
            
            !WRITE(*,*)DN3(1),DN3(2),sr3

            AIWL=AIWL-P1*G*GW(K+IX-1)*PHIKSI*DN3(1)*AJ
            AIWL1=AIWL1-P1*G*GW(K+IX-1)*PHIKSI1*DN3(1)*AJ
        END DO
    END DO
    AIWL=AIWL/(0.5*P1*U*U*SAREA0)
    AIWL1=AIWL1/(0.5*P1*U*U*SAREA0)

    WRITE(*,*)AIJX2,AIJX5
    !STOP

    IF(ITTE.EQ.1) THEN   !LINEAR CALCULATION
        SINK=(AIJX2*RESTO(5,5)-AIJX4*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        SINKL=(AIJX3*RESTO(5,5)-AIJX5*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))

        TRI=(RESTO(3,3)*AIJX4-RESTO(5,3)*AIJX2)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        TRIL=(RESTO(3,3)*AIJX5-RESTO(5,3)*AIJX3)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))

        !SINK0=SINK;SINKL0=SINKL
        !TRI0=TRI;TRIL0=TRIL
    ELSE
        SINK=(AIJX2*RESTO(5,5)-AIJX4*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        SINKL=(AIJX3*RESTO(5,5)-AIJX5*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))

        TRI=(RESTO(3,3)*AIJX4-RESTO(5,3)*AIJX2)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        TRIL=(RESTO(3,3)*AIJX5-RESTO(5,3)*AIJX3)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
    END IF
            !SINK=AIJX2/RESTO(3,3)
            !SINKL=AIJX3/RESTO(3,3)
            !TRI=AIJX4/RESTO(5,5)
            !TRIL=AIJX5/RESTO(5,5)

	WRITE(*,*)'AIJX=',AIJX/(0.5*P1*U*U*SAREA0)
	WRITE(*,*)'AIJX1=',AIJX1/(0.5*P1*U*U*SAREA0)
    WRITE(*,*)'SINKING=',AIJX2/RESTO(3,3)
    WRITE(*,*)'trim=',AIJX3/RESTO(5,5) !/0.75

!CSAREA=0.5*PI*((RB/RA)**2)*((2*RA)**2)+0.5*PI*(RB/RA)*(2*RA)
!!!C &	*ASIN(SQRT(1-(RB/RA)**2))/SQRT(1-(RB/RA)**2)
	!SAREA=0.0743953*2
	WRITE(*,*)'AREA OF WIGLEY IS',SAREA0

	CW=AIJX/(0.5*P1*U*U*SAREA0)
	CW1=AIJX1/(0.5*P1*U*U*SAREA0)
	!CW1=AIJX1 !/(0.5*P1*U*U*SAREA)

    !CW=AIJX/(0.5*P1*U*U*WFA)
	!CW1=AIJX1/(0.5*P1*U*U*WFA)
		
    WRITE(*,*)FR,'R/(0.5PU^2S)=',CW
	WRITE(*,*)FR,'R/(0.5PU^2S)LINEAR=',CW1
	!OPEN(2,FILE='R.DAT')
	
    !WRITE(201,101)FR,CW,CW1,AIWL,AIWL1,CW-AIWL,CW1-AIWL1
    101format(F15.6,'  R/(0.5PU^2S)=',F15.6,'    R/(0.5PU^2S)LINEAR=  ',5F15.6)
    !WRITE(202,102)FR,AIJX2/RESTO(3,3)/RL,AIJX3/RESTO(3,3)/RL,&
    !                         AIJX4/RESTO(5,5),AIJX5/RESTO(5,5) 
    
    !WRITE(202,102)FR,SINK/RL,SINKL/RL,TRI,TRIL
    
    102FORMAT(5F15.6)
	!WRITE(2,*)NI,'R/(0.5PU^2S)LINEAR=',CW1
	!CLOSE(2)
	!DEALLOCATE(SP)
	DEALLOCATE(CORL,CNL)
    
	!DEALLOCATE(SH,SA,SX,SZ,SY)
    !DEALLOCATE(SPX,SPY,SPZ)

	WRITE(*,*)'		END MATRIX'
    
    CLOSE(16)
    CLOSE(99)
    RETURN
END SUBROUTINE


SUBROUTINE PRESSURE_4
	USE GREENMOD
	USE CUMOD
    USE INPUTDATA
    REAL:: XI,ETA,AJ,AIWL1,AIWL
    common /gwa/ga(400),gw(400)
	DIMENSION ISM(3,2)
	DATA ISM/1,1,1,1,-1,1/

    REAL,dimension(:,:):: xxe(2,3)
    REAL,dimension(:):: sn(9),snx(9),sne(9),sr(3),dn(3),VV(3),SPL4(9),SPL(9)
    REAL,dimension(:):: sn3(3),snx3(3),sne3(3),sr3(2),dn3(2),FCORL(3,2)
	REAL,dimension(:,:):: a(3,3),b(3,3),CORL4(4,3)
	DIMENSION IS(3),JS(3),FRES(2),FRES0(2)
	INTEGER:: L,ISYM,IPAN
    CHARACTER*7 TITLE,TITLE2


    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(99,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\'//TRIM(ADJUSTL(TITLE))//'\WAVE_PRO_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
    !OPEN(7,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\WAVEPRO_NONLINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    
    !OPEN(7,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\WAVEPRO_NONLINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
      

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
    AIJX6=0.0
    AIJX7=0.0
    AIJX8=0.0
    !XGRA=0.
    !AIJX2=-MASS0(3,3)*G
    !AIJX3=-MASS0(3,3)*G
    !AIJX4=-MASS0(5,5)*G
    !AIJX5=-MASS0(5,5)*G

    !WRITE(*,*)'NUMBER OF WETTED SURFACE PANELS: ',NBBLOCKW*4

	DO 31 IE=1,NBBLOCKW
    DO 31 ISYM=1,1
        DO 31 IPAN=1,4 
        !DO 31 IPAN=2,2

!!CWRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!!C依次将每一面元的四个顶点的坐标赋值到corl(4,3)中
!!C依次将每一面元的四个顶点的法向矢量赋值到cnl(4,3)中
		
		DO I1=1,PAN(IE)
			I2=I1
			!IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1
			IF(PAN(THR1+IE).EQ.9) THEN
                IF(ISYM.EQ.2.AND.(I1.NE.1).AND.(I1.NE.9)) I2=10-I1
            END IF
            IF(PAN(IE).EQ.8) THEN
                IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
            END IF

            DO J1=1,3
				CORL(I1,J1)=CORD(INE(THR1+IE,I2),J1)*ISM(J1,ISYM)
				CNL(I1,J1)=VECN(INE(THR1+IE,I2),J1)*ISM(J1,ISYM)  
            END DO
            SPL(I1)=SP(INE(IE,I2))
        END DO
  
        SELECT CASE(IPAN)
        CASE(1)
        CORL4(1,:)=CORL(1,:);CORL4(2,:)=CORL(2,:);CORL4(3,:)=CORL(9,:);CORL4(4,:)=CORL(8,:)
        SPL4(1)=SPL(1);SPL4(2)=SPL(2);SPL4(3)=SPL(9);SPL4(4)=SPL(8)

        CASE(2)
        CORL4(1,:)=CORL(2,:);CORL4(2,:)=CORL(3,:);CORL4(3,:)=CORL(4,:);CORL4(4,:)=CORL(9,:)
        SPL4(1)=SPL(2);SPL4(2)=SPL(3);SPL4(3)=SPL(4);SPL4(4)=SPL(9)
        
        CASE(3)
        CORL4(1,:)=CORL(4,:);CORL4(2,:)=CORL(5,:);CORL4(3,:)=CORL(6,:);CORL4(4,:)=CORL(9,:)
        SPL4(1)=SPL(4);SPL4(2)=SPL(5);SPL4(3)=SPL(6);SPL4(4)=SPL(9)
        
        CASE(4)
        CORL4(1,:)=CORL(6,:);CORL4(2,:)=CORL(7,:);CORL4(3,:)=CORL(8,:);CORL4(4,:)=CORL(9,:)
        SPL4(1)=SPL(6);SPL4(2)=SPL(7);SPL4(3)=SPL(8);SPL4(4)=SPL(9)
        END SELECT

		ngus=21
		k=ngus*(ngus-1)/2
!gauss-legendre numerical integration

!对x积分     
		do 42 ix=1,ngus
		xi=ga(k+ix-1)
!对y积分
        !xi=1.

		AY=0.0
		do 42 iy=1,ngus
			eta=ga(k+iy-1)

            !eta=1.
!调用isopar，求得形状函数J及单位法向矢量等     
			!IF(PAN(IE).EQ.9) call isopar_9(xi,eta,corl,sn,snx,sne,aj,sr,dn)
            
			call isopar_4(xi,eta,corl4(1:4,:),sn(1:4),snx(1:4),sne(1:4),aj,sr,dn)
			!IF(PAN(IE).EQ.8) call isopar_8(xi,eta,corl(1:PAN(IE),:),sn(1:PAN(IE)),snx(1:PAN(IE)),sne(1:PAN(IE)),aj,sr,dn)
			do 210 i=1,3
				xxe(1,i)=0.0
				xxe(2,i)=0.0
210			continue

			do 21 i=1,3
			do 21 j=1,4 !PAN(IE)
!!!C	xxe(1,*): xxi,yxi,zxi
!!!C	xxe(2,*): xeta,yeta,zeta
				xxe(1,i)=xxe(1,i)+corl4(j,i)*snx(j)
   				xxe(2,i)=xxe(2,i)+corl4(j,i)*sne(j)
   21			continue
!!!C
			DO 22 I=1,3
			DO 22 J=1,2
				B(I,J)=XXE(J,I)
22			CONTINUE

			DO 23 I=1,3
				B(I,3)=DN(I)
23			CONTINUE

			DO 	230 I=1,3
			DO  230 J=1,3
				A(I,J)=B(I,J)
!!!C	write(*,*)a(i,j)
230			CONTINUE

			CALL BRINV(A,3)
			
			PHIKSI=0.0
			PHIETA=0.0
			DO 24 I1=1,4
			    !I2=I1
			    !IF((ISYM.EQ.2).AND.(I1.NE.1).AND.(I1.NE.9)) I2=10-I1
				IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1	
				
                !PHIKSI=PHIKSI+SNX(I1)*SP(INE(IE,I2))
				!PHIETA=PHIETA+SNE(I1)*SP(INE(IE,I2))
                PHIKSI=PHIKSI+SNX(I1)*SPL4(I1)
				PHIETA=PHIETA+SNE(I1)*SPL4(I1)
24			CONTINUE
			DO I=1,3
				VV(I)=A(1,I)*PHIKSI+A(2,I)*PHIETA+A(3,I)*U*DN(1)  !VV速度势的散度
            END DO
            
            !VV(1)=PHIKSX
            !VV(2)=PHIKSY
            !VV(3)=PHIKSZ

            !此处有问题
			BNR=0.0-P1*(0.5*(VV(1)**2+VV(2)**2&    !非线性阻力
     			 +VV(3)**2)-U*VV(1))*DN(1)-P1*G*(SR(3))*DN(1) !-P1*G*(SR(3)-DCORZ(INE(IE,9)))*DN(1) !
			!BNR=0.0-P1*(0.5*(VV(1)**2+VV(2)**2&    !非线性阻力
     		!	 +VV(3)**2)-U*VV(1))*DN(1) !-P1*G*(SR(3))*DN(1) !-P1*G*(SR(3)-DCORZ(INE(IE,9)))*DN(1) !
			!BNR=-P1*G*(SR(3))*DN(1) !-P1*G*(SR(3)-DCORZ(INE(IE,9)))*DN(1) !
            
            BNR1=0.0+P1*(U*VV(1))*DN(1)	   !线性阻力
            
            BNR2=0.0-P1*(0.5*(VV(1)**2+VV(2)**2&        !非线性垂向力
     			    +VV(3)**2)-U*VV(1))*DN(3) !-P1*G*(SR(3))*DN(3)  
            !BNR2=-P1*G*(SR(3))*DN(3)  

            BNR3=0.0+P1*(U*VV(1))*DN(3) !-P1*G*(SR(3))*DN(3) 	            !线性垂向力

            BNR4=0.0-P1*(0.5*(VV(1)*VV(1)+VV(2)*VV(2)&           !非线性纵摇力
     			 +VV(3)*VV(3))-U*VV(1))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3)) !&
     			 !+VV(3)*VV(3))-U*VV(1))*((SR(3))*DN(1)-(SR(1))*DN(3)) !&
                 !-P1*G*(SR(3))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3))

            BNR5=0.0+P1*(U*VV(1))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3)) !&
            !BNR5=0.0+P1*(U*VV(1))*((SR(3))*DN(1)-(SR(1))*DN(3))
                 !-P1*G*(SR(3))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3))
            !BNR5=1 !DN(2)*SR(2)

            BNR6=-P1*G*SR(3)*DN(3)
            BNR7=-P1*G*SR(3)*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3))
            !BNR7=-P1*G*SR(3)*((SR(3))*DN(1)-(SR(1))*DN(3))
            BNR8=-P1*G*(SR(3))*DN(1)

            !BNR6=-P1*G*SR(3)*DN(3)
            !BNR7=-P1*G*SR(2)*DN(2)

			AIJY=gw(k+iy-1)*aj*BNR
			AIJY1=GW(K+IY-1)*AJ*BNR1
            AIJY2=GW(K+IY-1)*AJ*BNR2
            AIJY3=GW(K+IY-1)*AJ*BNR3
            AIJY4=GW(K+IY-1)*AJ*BNR4
            AIJY5=GW(K+IY-1)*AJ*BNR5

            AIJY6=GW(K+IY-1)*AJ*BNR6
            AIJY7=GW(K+IY-1)*AJ*BNR7
            AIJY8=GW(K+IY-1)*AJ*BNR8

  			AIJX=AIJX+gw(k+ix-1)*AIJY
			AIJX1=AIJX1+GW(K+IX-1)*AIJY1
			AIJX2=AIJX2+GW(K+IX-1)*AIJY2
			AIJX3=AIJX3+GW(K+IX-1)*AIJY3
			AIJX4=AIJX4+GW(K+IX-1)*AIJY4
			AIJX5=AIJX5+GW(K+IX-1)*AIJY5
			IF(ISYM.EQ.1) THEN !HALF BODY
            AIJX6=AIJX6+GW(K+IX-1)*AIJY6
			AIJX7=AIJX7+GW(K+IX-1)*AIJY7
			AIJX8=AIJX8+GW(K+IX-1)*AIJY8
            END IF

            IF(MOD(IE,5).EQ.0) THEN
                !WRITE(*,*)IE,SR(3),SR(1),(0.0-P1*(0.5*(VV(1)**2+VV(2)**2&    !非线性阻力
     			! +VV(3)**2)-U*VV(1))*DN(1)-P1*G*(SR(3)))/(0.5*P1*U*U*SAREA0) !-1./(2*G)*P1*U**2
            END IF

            !IF(ITTE.NE.1) THEN   !NONLINEAR CALCULATION
  			!    FRES(1)=FRES(1)+GW(k+ix-1)*GW(k+iy-1)*AJ*(-SR(3))*DN(3)
  			!    FRES(2)=FRES(2)+GW(k+ix-1)*GW(k+iy-1)*AJ*(-SR(3))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3))
            !END IF

            !AIJX  非线性阻力
            !AIJY1 线性阻力
            !AIJY2 非线性垂向力
            !AIJY3 线性垂向力
            !AIJY4 非线性纵摇力
            !AIJY5 线性纵摇力
            !AIJY6 浮力，z
            !AIJY7 浮力，yy

   42		continue
   31	CONTINUE

    !WRITE(*,*)AIJX6,AIJX7
!对x积分 
    AIJX=AIJX*2.
    AIJX1=AIJX1*2.  
    AIJX2=AIJX2*2.
    AIJX3=AIJX3*2.  
    AIJX4=AIJX4*2.
    AIJX5=AIJX5*2.    
    AIJX6=AIJX6*2.
    AIJX7=AIJX7*2.
    AIJX8=AIJX8*2.
    !WRITE(*,*)AIJX6*7.2786**3/2.,AIJX7,-MASS0(3,3)*XGRA(1)
    !WRITE(*,*)AIJX6,AIJX,MASS0(1,1)
    !STOP

    IF(ITTE.EQ.1) THEN
        MASS0(1,1)=AIJX8
        MASS0(3,3)=-AIJX6/G
        MASS0(5,5)=-AIJX7/G
    END IF

    !AIJX=AIJX+MASS0(1,1)*G

    IF(ITTE.EQ.1) THEN   !LINEAR CALCULATION

        !AIJX2=AIJX2+(AIJX6-MASS0(3,3)*G)
        !AIJX3=AIJX3+(AIJX6-MASS0(3,3)*G)
        !AIJX4=AIJX4+AIJX7 !+(MASS0(3,3)*XGRA(1)*G)
        !AIJX5=AIJX5+AIJX7 !+(MASS0(3,3)*XGRA(1)*G)
        SINK=(AIJX2*RESTO(5,5)-AIJX4*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        SINKL=(AIJX3*RESTO(5,5)-AIJX5*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))

        TRI=(RESTO(3,3)*AIJX4-RESTO(5,3)*AIJX2)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        TRIL=(RESTO(3,3)*AIJX5-RESTO(5,3)*AIJX3)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))

        !SINK0=SINK;SINKL0=SINKL
        !TRI0=TRI;TRIL0=TRIL
    ELSE
        !AIJX2=AIJX2+(AIJX6+MASS0(3,3)*G)
        !AIJX3=AIJX3+(AIJX6+MASS0(3,3)*G)
        !AIJX4=AIJX4+AIJX7 !+(MASS0(5,5)*G)
        !AIJX5=AIJX5+AIJX7 !+(MASS0(5,5)*G)
        
        !IF(MDFREE.EQ.1) THEN
        IF(MDFREE.EQ.1.OR.MDFREE.EQ.3) THEN
        
        AIJX2T=AIJX2
        AIJX3T=AIJX3
        AIJX4T=AIJX4
        AIJX5T=AIJX5


        AIJX2=AIJX2+AIJX6+MASS0(3,3)*G
        AIJX3=AIJX3+AIJX6+MASS0(3,3)*G
        AIJX4=AIJX4+AIJX7+MASS0(5,5)*G
        AIJX5=AIJX5+AIJX7+MASS0(5,5)*G
        !AIJX2=AIJX2-AIJX20+AIJX6-AIJX60
        !AIJX3=AIJX3-AIJX30+AIJX6-AIJX60
        !AIJX4=AIJX4-AIJX40+AIJX7-AIJX70
        !AIJX5=AIJX5-AIJX50+AIJX7-AIJX70

        !SINK=SINK+0.75*(AIJX2*RESTO(5,5)-AIJX4*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        SINKL=SINKL+(AIJX3*RESTO(5,5)-AIJX5*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))

        !TRI=TRI+0.75*(RESTO(3,3)*AIJX4-RESTO(5,3)*AIJX2)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        TRIL=TRIL+(RESTO(3,3)*AIJX5-RESTO(5,3)*AIJX3)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        
        SINK=SINK+NSWIN(16,1)*(AIJX2*RESTO(5,5)-AIJX4*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        TRI=TRI+NSWIN(16,1)*(RESTO(3,3)*AIJX4-RESTO(5,3)*AIJX2)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))   
        
        !SINK=(AIJX2*RESTO(5,5)-AIJX4*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        !SINKL=(AIJX3*RESTO(5,5)-AIJX5*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        !TRI=(RESTO(3,3)*AIJX4-RESTO(5,3)*AIJX2)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        !TRIL=(RESTO(3,3)*AIJX5-RESTO(5,3)*AIJX3)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
    
        ELSE 
        !AIJX2=AIJX2+AIJX6+MASS0(3,3)*G
        !AIJX4=AIJX4+AIJX7+MASS0(5,5)*G

        SINK=(AIJX2*RESTO(5,5)-AIJX4*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        SINKL=(AIJX3*RESTO(5,5)-AIJX5*RESTO(3,5))/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))

        TRI=(RESTO(3,3)*AIJX4-RESTO(5,3)*AIJX2)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        TRIL=(RESTO(3,3)*AIJX5-RESTO(5,3)*AIJX3)/(RESTO(3,3)*RESTO(5,5)-RESTO(5,3)*RESTO(3,5))
        END IF
        
    END IF
 
        AIJX20=AIJX2T
        AIJX30=AIJX3T
        AIJX40=AIJX4T
        AIJX50=AIJX5T
        AIJX60=AIJX6
        AIJX70=AIJX7
        !AIJX60=-MASS(3,3)*G
        !AIJX70=-MASS(5,5)*G
            !SINK=AIJX2/RESTO(3,3)
            !SINKL=AIJX3/RESTO(3,3)
            !TRI=AIJX4/RESTO(5,5)
            !TRIL=AIJX5/RESTO(5,5)

    ALWX1=0.
    IF(ITTE.GE.3) THEN   !OLD
    !IF(ITTE.GE.1) THEN
    DO I=NFXF+1,NFXF+NFX-1
        DWL=SQRT((WPLX(I)-WPLX(I+1))**2+(WPLY(I)-WPLY(I+1))**2)

        !ALWX1=ALWX1+2.*0.5*P1*G*0.5*(WPLZ(I)**2+WPLZ(I+1)**2)*DWL*(WPLY(I+1)-WPLY(I))/SQRT((WPLX(I)-WPLX(I+1))**2+(WPLY(I)-WPLY(I+1))**2)
        !WRITE(*,*)WPLY(I+1),WPLY(I)
    END DO
    END IF

    
    IF(ITTE.EQ.1) THEN
    DO I=1,NTAO
        DWL=SQRT((XTR(I)-XTR(I+1))**2)
        ALW=ALW+2.*0.5*P1*G*0.5*(ABS(ZTR(I))*ZTR(I)+ABS(ZTR(I+1))*ZTR(I+1))*DWL

        !ALW=ALW+2.*0.5*P1*G*0.5*(WPLZ(I)**2+WPLZ(I+1)**2)*DWL*(WPLY(I+1)-WPLY(I))/SQRT((WPLX(I)-WPLX(I+1))**2+(WPLY(I)-WPLY(I+1))**2)
        
        !WRITE(*,*)WPLY(I+1),WPLY(I)
    END DO
    END IF    
    
    !WRITE(*,*)'ALW=',ALW/(0.5*P1*U*U*SAREA0)
    !WRITE(*,*)
	!WRITE(*,*)'AIJX=',(AIJX+ALW)/(0.5*P1*U*U*SAREA0) !,ALW/(0.5*P1*U*U*SAREA0)
	!WRITE(*,*)'AIJX1=',AIJX1/(0.5*P1*U*U*SAREA0)
    !WRITE(*,*)'SINKING=',AIJX2/RESTO(3,3)
    !WRITE(*,*)'trim=',AIJX3/RESTO(5,5) !/0.75
    !STOP
    
!CSAREA=0.5*PI*((RB/RA)**2)*((2*RA)**2)+0.5*PI*(RB/RA)*(2*RA)
!!!C &	*ASIN(SQRT(1-(RB/RA)**2))/SQRT(1-(RB/RA)**2)
	!SAREA=0.0743953*2
	!WRITE(*,*)'AREA OF WIGLEY IS',SAREA

    
    !IF(ITTTE.GT.1) THEN
	CW=(AIJX+ALW-MASS0(1,1))/(0.5*P1*U*U*SAREA0)
    !ELSE
	!CW=(AIJX+ALW)/(0.5*P1*U*U*SAREA0)
    !END IF
	CW1=AIJX1/(0.5*P1*U*U*SAREA0)
	!CW1=AIJX1 !/(0.5*P1*U*U*SAREA)

    !CW=AIJX/(0.5*P1*U*U*WFA)
	!CW1=AIJX1/(0.5*P1*U*U*WFA)
	write(*,*)CW,CW1

    
    WRITE(9999,*)
    WRITE(9999,*)
    WRITE(9999,*)"CW OF PRESSURE INTEGRAL :     ",CW !AIJX/(0.5*P1*U*U*SAREA0)
    IF(MTROM.EQ.1) THEN
    WRITE(9999,*)"CW OF TRANSOM CORRECTION:     ",ALW/(0.5*P1*U*U*SAREA0)
	END IF
    WRITE(9999,*)"SINKAGE OF Lpp/2        :     ",-SINK
    WRITE(9999,*)"TRIM BY STERN           :     ",-TRI
	
    !WRITE(201,101)FR,CW,CW1,AIWL,AIWL1,CW-AIWL,CW1-AIWL1
    101format(F15.6,'  R/(0.5PU^2S)=',F15.6,'    R/(0.5PU^2S)LINEAR=  ',5F15.6)
    !WRITE(202,102)FR,AIJX2/RESTO(3,3)/RL,AIJX3/RESTO(3,3)/RL,&
    !                         AIJX4/RESTO(5,5),AIJX5/RESTO(5,5) 
    
    !WRITE(202,102)FR,SINK/RL,SINKL/RL,TRI,TRIL
    
    102FORMAT(5F15.6)
	!WRITE(2,*)NI,'R/(0.5PU^2S)LINEAR=',CW1
	!CLOSE(2)
	!DEALLOCATE(SP)
	DEALLOCATE(CORL,CNL)
    
	!DEALLOCATE(SH,SA,SX,SZ,SY)
    !DEALLOCATE(SPX,SPY,SPZ)

	!WRITE(*,*)'		END MATRIX'
    
    CLOSE(16)
    CLOSE(99)
    RETURN
END SUBROUTINE

SUBROUTINE RESTORING

	USE GREENMOD
	USE CUMOD
    REAL:: XI,ETA,AJ,AIWL1,AIWL
    common /gwa/ga(400),gw(400)
	DIMENSION ISM(3,2)
	DATA ISM/1,1,1,1,-1,1/

    REAL,dimension(:,:):: xxe(2,3),TPWPLZ(100)
    REAL,dimension(:):: sn(9),snx(9),sne(9),sr(3),dn(3),VV(3),SPL4(9),SPL(9)
    REAL,dimension(:):: sn3(3),snx3(3),sne3(3),sr3(2),dn3(2),FCORL(3,2)
	REAL,dimension(:,:):: a(3,3),b(3,3),CORL4(4,3)
	DIMENSION IS(3),JS(3),FRES(2),FRES0(2)
	INTEGER:: L,ISYM,IPAN,TPOIN
    CHARACTER*7 TITLE,TITLE2
    REAL,ALLOCATABLE:: CORDR(:,:)
    INTEGER,DIMENSION(:,:),ALLOCATABLE:: INER

    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(99,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\'//TRIM(ADJUSTL(TITLE))//'\WAVE_PRO_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
    !OPEN(7,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\WAVEPRO_NONLINEAR_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')


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

    MMESH=2
    NBPOINTR=NBPOINT
    NBBLOCKR=NBBLOCK        !WPLAZ=0.    TPWPLZ=WPLZ    WPLZ(:)=0.    !CALL OFFSHIP_NOL
    CALL OFFWIGLEY_9NODES_NON
    OPEN(11,FILE='ALLINDEX.DAT')
    MMESH=1    

    READ(11,*)TPOIN,TP

	ALLOCATE(CORDR(TPOIN,3),INER(NBBLOCK,9))
!C读入第一层流体中物体上及界面上各点的坐标和面元顶点索引


	DO I=1,TPOIN
	    READ(11,*)CORDR(I,1),CORDR(I,2),CORDR(I,3) !,VECN(I,1),VECN(I,2),VECN(I,3)
    END DO
  
	DO I=1,NBBLOCK
        READ(11,*)TP,INER(I,1:9)
	END DO

    THR1=0.

	AIJX=0.0
	AIJX1=0.0
    AIJX2=0.0
    AIJX3=0.0
    AIJX4=0.0
    AIJX5=0.0
    AIJX6=0.0
    AIJX7=0.0
    AIJX8=0.0
    !XGRA=0.
    !AIJX2=-MASS0(3,3)*G
    !AIJX3=-MASS0(3,3)*G
    !AIJX4=-MASS0(5,5)*G
    !AIJX5=-MASS0(5,5)*G

    !WRITE(*,*)'NUMBER OF WETTED SURFACE PANELS: ',NBBLOCKW*4

	DO 31 IE=1,NBBLOCK
    DO 31 ISYM=1,1
        DO 31 IPAN=1,4 
        !DO 31 IPAN=2,2

!!CWRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!!C依次将每一面元的四个顶点的坐标赋值到corl(4,3)中
!!C依次将每一面元的四个顶点的法向矢量赋值到cnl(4,3)中
		
		DO I1=1,9
			I2=I1

            DO J1=1,3
				CORL(I1,J1)=CORDR(INER(THR1+IE,I2),J1)*ISM(J1,ISYM)
				CNL(I1,J1)=VECN(INER(THR1+IE,I2),J1)*ISM(J1,ISYM)  
            END DO
        END DO
  
        SELECT CASE(IPAN)
        CASE(1)
        CORL4(1,:)=CORL(1,:);CORL4(2,:)=CORL(2,:);CORL4(3,:)=CORL(9,:);CORL4(4,:)=CORL(8,:)

        CASE(2)
        CORL4(1,:)=CORL(2,:);CORL4(2,:)=CORL(3,:);CORL4(3,:)=CORL(4,:);CORL4(4,:)=CORL(9,:)
        
        CASE(3)
        CORL4(1,:)=CORL(4,:);CORL4(2,:)=CORL(5,:);CORL4(3,:)=CORL(6,:);CORL4(4,:)=CORL(9,:)
        
        CASE(4)
        CORL4(1,:)=CORL(6,:);CORL4(2,:)=CORL(7,:);CORL4(3,:)=CORL(8,:);CORL4(4,:)=CORL(9,:)
        END SELECT

		ngus=21
		k=ngus*(ngus-1)/2
!gauss-legendre numerical integration

!对x积分     
		do 42 ix=1,ngus
		xi=ga(k+ix-1)
!对y积分
        !xi=1.

		AY=0.0
		do 42 iy=1,ngus
			eta=ga(k+iy-1)

            !eta=1.
!调用isopar，求得形状函数J及单位法向矢量等     
            
			call isopar_4(xi,eta,corl4(1:4,:),sn(1:4),snx(1:4),sne(1:4),aj,sr,dn)
			!IF(PAN(IE).EQ.8) call isopar_8(xi,eta,corl(1:PAN(IE),:),sn(1:PAN(IE)),snx(1:PAN(IE)),sne(1:PAN(IE)),aj,sr,dn)
			

            BNR6=-P1*G*SR(3)*DN(3)
            BNR7=-P1*G*SR(3)*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3))
            !BNR7=-P1*G*SR(3)*((SR(3))*DN(1)-(SR(1))*DN(3))


            AIJY6=GW(K+IY-1)*AJ*BNR6
            AIJY7=GW(K+IY-1)*AJ*BNR7

            AIJX6=AIJX6+GW(K+IX-1)*AIJY6
			AIJX7=AIJX7+GW(K+IX-1)*AIJY7

            !IF(ITTE.NE.1) THEN   !NONLINEAR CALCULATION
  			!    FRES(1)=FRES(1)+GW(k+ix-1)*GW(k+iy-1)*AJ*(-SR(3))*DN(3)
  			!    FRES(2)=FRES(2)+GW(k+ix-1)*GW(k+iy-1)*AJ*(-SR(3))*((SR(3)-XGRA(3))*DN(1)-(SR(1)-XGRA(1))*DN(3))
            !END IF

            !AIJX  非线性阻力
            !AIJY1 线性阻力
            !AIJY2 非线性垂向力
            !AIJY3 线性垂向力
            !AIJY4 非线性纵摇力
            !AIJY5 线性纵摇力
            !AIJY6 浮力，z
            !AIJY7 浮力，yy

   42		continue
   31	CONTINUE

    !WRITE(*,*)AIJX6,AIJX7
!对x积分 
    
    AIJX6=AIJX6*2.
    AIJX7=AIJX7*2.
    !WRITE(*,*)AIJX6*7.2786**3/2.,AIJX7,-MASS0(3,3)*XGRA(1)
    !WRITE(*,*)AIJX6,AIJX,MASS0(1,1)
    !STOP

    MASS(3,3)=-AIJX6/G
    MASS(5,5)=-AIJX7/G

    WRITE(*,*)"DISPLACEMENT:         ",MASS(3,3)
    !STOP

    DEALLOCATE(CORDR,INER)
	DEALLOCATE(CORL,CNL)
    
	!DEALLOCATE(SH,SA,SX,SZ,SY)
    !DEALLOCATE(SPX,SPY,SPZ)
    
    !WPLZ=TPWPLZ
    NBPOINT=NBPOINTR
    NBBLOCK=NBBLOCKR
	WRITE(*,*)'		END MATRIX'
    !stop

    CLOSE(11)
    !CLOSE(99)
    RETURN


END SUBROUTINE

SUBROUTINE BRINV(A,N)
	DIMENSION A(N,N),IS(N),JS(N)
	!DOUBLE PRECISION A,T,D
        
    DO I=1,N
        !IF(N.EQ.6) WRITE(*,112)A(I,:)
    END DO
    112FORMAT(6F15.4)
    
    L=1
	DO 100 K=1,N
	  D=0.0
	  DO 10 I=K,N
	  DO 10 J=K,N
	    IF (ABS(A(I,J)).GT.D) THEN
	      D=ABS(A(I,J))
	      IS(K)=I
	      JS(K)=J
	    END IF
10	  CONTINUE
	  IF (D+1.0.EQ.1.0) THEN
	    L=0
	    !WRITE(*,20)
        !WRITE(*,*)A
        
        !STOP
	    RETURN
	  END IF
20	  FORMAT(1X,'ERR**NOT INV')
	  DO 30 J=1,N
	    T=A(K,J)
	    A(K,J)=A(IS(K),J)
	    A(IS(K),J)=T
30	  CONTINUE
	  DO 40 I=1,N
	    T=A(I,K)
	    A(I,K)=A(I,JS(K))
	    A(I,JS(K))=T
40	  CONTINUE
	  A(K,K)=1/A(K,K)
	  DO 50 J=1,N
	    IF (J.NE.K) THEN
	      A(K,J)=A(K,J)*A(K,K)
	    END IF
50	  CONTINUE
	  DO 70 I=1,N
	    IF (I.NE.K) THEN
	      DO 60 J=1,N
	        IF (J.NE.K) THEN
	          A(I,J)=A(I,J)-A(I,K)*A(K,J)
	        END IF
60	      CONTINUE
	    END IF
70	  CONTINUE
	  DO 80 I=1,N
	    IF (I.NE.K) THEN
	      A(I,K)=-A(I,K)*A(K,K)
	    END IF
80	  CONTINUE
100	CONTINUE
	DO 130 K=N,1,-1
	  DO 110 J=1,N
	    T=A(K,J)
	    A(K,J)=A(JS(K),J)
	    A(JS(K),J)=T
110	  CONTINUE
	  DO 120 I=1,N
	    T=A(I,K)
	    A(I,K)=A(I,IS(K))
	    A(I,IS(K))=T
120	  CONTINUE
130	CONTINUE
	RETURN
END SUBROUTINE

SUBROUTINE BVELOCITY
	USE GREENMOD
	USE CUMOD
    REAL:: XI,ETA,AJ,AIWL1,AIWL
    common /gwa/ga(400),gw(400)
	DIMENSION ISM(3,2)
	DATA ISM/1,1,1,1,-1,1/

    REAL,dimension(:,:):: xxe(2,3)
    REAL,dimension(:):: sn(9),snx(9),sne(9),sr(3),dn(3),VV(3)
    REAL,dimension(:):: sn3(3),snx3(3),sne3(3),sr3(2),dn3(2),FCORL(3,2)
	REAL,dimension(:,:):: a(3,3),b(3,3)
	DIMENSION IS(3),JS(3),FRES(2),FRES0(2)
    REAL:: NX1L
	INTEGER:: L,ISYM
    CHARACTER*7 TITLE,TITLE2


    WRITE(TITLE,'(F7.3)')FR
    WRITE(TITLE2,'(I2)')ITTE
    !OPEN(99,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\'//TRIM(ADJUSTL(TITLE))//'\WAVE_PRO_FN='//TRIM(ADJUSTL(TITLE))//'.DAT')
    !OPEN(7,FILE='e:\SHIP_WDRAG_JOURNAL_NON\'//TRIM(ADJUSTL(NAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\VELOCITY_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    OPEN(7,FILE=''//TRIM(ADJUSTL(ADDRESS))//'\'//TRIM(ADJUSTL(FOLDERNAME))//'\FN_'//TRIM(ADJUSTL(TITLE))//'\VELOCITY_ITER='//TRIM(ADJUSTL(TITLE2))//'.DAT')
    
    
    101	format(1x,'title="episode solid mesh"'/1x,'variables="X","Y",& 
     "Z","VX","VY","VZ","NX","NY","NZ"'/1x,'zone t="panel",n=',i8,',e=',i8,',f=fepoint,& 
      et=quadrilateral') 
    WRITE(7,101)NBPOINT,NBBLOCK*4 

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

	!AIJX2=0.0
    !AIJX3=0.0
    !AIJX4=0.0
    !AIJX5=0.0

    !AIJX2=-MASS0(3,3)*G
    !AIJX3=-MASS0(3,3)*G
    !AIJX4=-MASS0(5,5)*G
    !AIJX5=-MASS0(5,5)*G

    FRES=0.


	DO 31 IE=1,NBBLOCK
    DO 31 ISYM=1,1
!!CWRITE(*,*)'		ELEMENT NO. ',IE,'    OF',NTP
!!C依次将每一面元的四个顶点的坐标赋值到corl(4,3)中
!!C依次将每一面元的四个顶点的法向矢量赋值到cnl(4,3)中
		
		DO I1=1,PAN(IE)
			I2=I1
			!IF(MOD(ISYM,2).EQ.0.AND.(I1.NE.1)) I2=6-I1
			IF(PAN(THR1+IE).EQ.9) THEN
                IF(ISYM.EQ.2.AND.(I1.NE.1).AND.(I1.NE.9)) I2=10-I1
            END IF
            IF(PAN(IE).EQ.8) THEN
                IF(ISYM.EQ.2.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                IF(ISYM.EQ.2.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列

                IF(ISYM.EQ.4.AND.I1.NE.1.AND.I1.LE.4) I2=6-I1   !保证镜像后的面元节点按逆时针排列
                IF(ISYM.EQ.4.AND.I1.GE.5) I2=13-I1   !保证镜像后的面元节点按逆时针排列
            END IF

            DO J1=1,3
				CORL(I1,J1)=CORD(INE(THR1+IE,I2),J1)*ISM(J1,ISYM)
				CNL(I1,J1)=VECN(INE(THR1+IE,I2),J1)*ISM(J1,ISYM)  
            END DO
        END DO

        DO 31 IK=1,9
                SELECT CASE(IK)
                CASE(1)
                IX=-1.;ETA=-1.
                CASE(2)
                IX=0.;ETA=-1.
                CASE(3)
                IX=1.;ETA=-1.
                CASE(4)
                IX=1.;ETA=0.
                CASE(5)
                IX=1.;ETA=1.
                CASE(6)
                IX=0.;ETA=1.
                CASE(7)
                IX=-1.;ETA=1.
                CASE(8)
                IX=-1.;ETA=0.
                CASE(9)
                IX=0.;ETA=0.
                END SELECT        
                
!调用isopar，求得形状函数J及单位法向矢量等     
			call isopar_9(xi,eta,corl,sn,snx,sne,aj,sr,dn)
			!IF(PAN(IE).EQ.8) call isopar_8(xi,eta,corl(1:PAN(IE),:),sn(1:PAN(IE)),snx(1:PAN(IE)),sne(1:PAN(IE)),aj,sr,dn)
			do 210 i=1,3
				xxe(1,i)=0.0
				xxe(2,i)=0.0
210			continue

			do 21 i=1,3
			do 21 j=1,PAN(IE)
!!!C	xxe(1,*): xxi,yxi,zxi
!!!C	xxe(2,*): xeta,yeta,zeta
				xxe(1,i)=xxe(1,i)+corl(j,i)*snx(j)
   				xxe(2,i)=xxe(2,i)+corl(j,i)*sne(j)
   21			continue
!!!C
			DO 22 I=1,3
			DO 22 J=1,2
				B(I,J)=XXE(J,I)
22			CONTINUE

			DO 23 I=1,3
				B(I,3)=DN(I)
23			CONTINUE

			DO 	230 I=1,3
			DO  230 J=1,3
				A(I,J)=B(I,J)
!!!C	write(*,*)a(i,j)
230			CONTINUE

			CALL BRINV(A,3)
			
			PHIKSI=0.0
			PHIETA=0.0
            PHIKSX=0.
            PHIKSY=0.
            PHIKSZ=0.
			DO 24 I1=1,PAN(IE)
			    I2=I1
			    IF((ISYM.EQ.2).AND.(I1.NE.1).AND.(I1.NE.9)) I2=10-I1
				
                PHIKSI=PHIKSI+SNX(I1)*SP(INE(IE,I2))
				PHIETA=PHIETA+SNE(I1)*SP(INE(IE,I2))

                PHIKSX=PHIKSX+SN(I1)*SPX(INE(IE,I2))
                PHIKSY=PHIKSY+SN(I1)*SPY(INE(IE,I2))
                PHIKSZ=PHIKSZ+SN(I1)*SPZ(INE(IE,I2))
24			CONTINUE
			DO I=1,3
				VV(I)=A(1,I)*PHIKSI+A(2,I)*PHIETA+A(3,I)*U*DN(1)  !VV速度势的散度
            END DO
            
            BVEL(INE(IE,IK),1)=-U+VV(1)
            BVEL(INE(IE,IK),2)=VV(2)
            BVEL(INE(IE,IK),3)=VV(3)


   31	CONTINUE

    DO I=1,NBPOINT
        WRITE(7,125)CORD(I,1:3),BVEL(I,1:3),VECN(I,1:3)      
    END DO
125 FORMAT(9F16.5)   
    DO I=1,NBBLOCK
        WRITE(7,*)INE(I,1),INE(I,2),INE(I,9),INE(I,8)
        WRITE(7,*)INE(I,2),INE(I,3),INE(I,4),INE(I,9)
        WRITE(7,*)INE(I,9),INE(I,4),INE(I,5),INE(I,6)
        WRITE(7,*)INE(I,8),INE(I,9),INE(I,6),INE(I,7)
    END DO
   
	DEALLOCATE(CORL,CNL)
    
	!DEALLOCATE(SH,SA,SX,SZ,SY)
    !DEALLOCATE(SPX,SPY,SPZ)

	WRITE(*,*)'		END MATRIX'
    
    CLOSE(7)
    CLOSE(16)
    CLOSE(99)
    RETURN
END SUBROUTINE


SUBROUTINE PHIMAX
	USE GREENMOD
	USE CUMOD
    REAL:: XI,ETA,AJ,AIWL1,AIWL
    common /gwa/ga(400),gw(400)
	DIMENSION ISM(3,2)
	DATA ISM/1,1,1,1,-1,1/

    REAL,dimension(:,:):: xxe(2,3)
    REAL,dimension(:):: sn(9),snx(9),sne(9),sr(3),dn(3),VV(3)
    REAL,dimension(:):: sn3(3),snx3(3),sne3(3),sr3(2),dn3(2),FCORL(3,2)
	REAL,dimension(:,:):: a(3,3),b(3,3)

    OPEN(9001,FILE='PHI_MAX.DAT')
    
	ALLOCATE (CORL(9,3),CNL(9,3))

!CSAREA=0.5*PI*((RB/RA)**2)*((2*RA)**2)+0.5*PI*(RB/RA)*(2*RA)
!!!C &	*ASIN(SQRT(1-(RB/RA)**2))/SQRT(1-(RB/RA)**2)
	!SAREA=0.0743953*2
	WRITE(*,*)'AREA OF WIGLEY IS',SAREA

    THR1=0.
    NPHI=3001
    DO IPHI=1,NPHI
        DPHI=(19.5-2.5)/180.*PI/(NPHI-1)
        PSIM=2.5/180.*PI+DPHI*(IPHI-1)
    
        DETPSI=SQRT(1.+4.*(TAN(PSIM))**2+SQRT(1.-8.*(TAN(PSIM))**2))
        
        PHIPP=2.*SQRT(2.)*SQRT(1.-8.*TAN(PSIM)**2)*ABS(SIN(PSIM))/DETPSI
        PHID=-(3.-SQRT(1.-8.*TAN(PSIM)**2))/8./SQRT(2.)*COS(PSIM)/ABS(TAN(PSIM))*DETPSI
        QDD=(1.+SQRT(1.-8.*TAN(PSIM)**2))/(4.*TAN(PSIM))
        PHIPP=(QDD*(3.+2.*QDD**2)*SIN(PSIM)-COS(PSIM))/((1.+QDD**2)**1.5)
        !WRITE(*,*)QDD,PHIPP,PHIPP1
        
        AIJX=0.
        AIJX1=0.
        
        
        DO IE1=1,NBBLOCK
        
		DO I1=1,9
			I2=I1

            DO J1=1,3
				CORL(I1,J1)=CORD(INE(IE1,I2),J1) !*ISM(J1,ISYM)
				CNL(I1,J1)=VECN(INE(IE1,I2),J1) !*ISM(J1,ISYM)  
            END DO
        END DO

		ngus=9
		k=ngus*(ngus-1)/2
!gauss-legendre numerical integration

!对x积分     
		do 42 ix=1,ngus
		xi=ga(k+ix-1)
!对y积分
		do 42 iy=1,ngus
			eta=ga(k+iy-1)
!调用isopar，求得形状函数J及单位法向矢量等     
			call isopar_9(xi,eta,corl,sn,snx,sne,aj,sr,dn)
			!IF(PAN(IE).EQ.8) call isopar_8(xi,eta,corl(1:PAN(IE),:),sn(1:PAN(IE)),snx(1:PAN(IE)),sne(1:PAN(IE)),aj,sr,dn)

            
            !VV(1)=PHIKSX
            !VV(2)=PHIKSY
            !VV(3)=PHIKSZ

            !此处有问题
            
            PHIDXX=SQRT(1.+QDD**2)*SR(1)/(FR**2)
            PHIDYY=SQRT(1.+QDD**2)*SR(2)*QDD/(FR**2)
            
			BNR=2./FR**4*EXP((1.+QDD**2)*(SR(3))/FR**2)*COS(PHIDYY)*COS(PHIDXX)*(-DN(1))
			BNR1=-2./FR**4*EXP((1.+QDD**2)*(SR(3))/FR**2)*COS(PHIDYY)*SIN(PHIDXX)*(-DN(1))
			!BNR=1./FR**4*EXP((1.+QDD**2)*(SR(3))/FR**2)*&
            !    (COS(PHIDYY)*COS(PHIDXX)-SIN(PHIDYY)*SIN(PHIDXX))*(-DN(1))
			!BNR1=-1./FR**4*EXP((1.+QDD**2)*(SR(3))/FR**2)*&
            !    (COS(PHIDYY)*SIN(PHIDXX)+COS(PHIDXX)*SIN(PHIDYY))*(-DN(1))

			AIJY=gw(k+iy-1)*aj*BNR
			AIJY1=gw(k+iy-1)*aj*BNR1
			!AIJY=gw(k+iy-1)*aj*COS(PHIDYY)*COS(PHIDXX)
			!AIJY1=gw(k+iy-1)*aj*COS(PHIDYY)*SIN(PHIDXX)

  			AIJX=AIJX+gw(k+ix-1)*AIJY
  			AIJX1=AIJX1+gw(k+ix-1)*AIJY1
            !WRITE(*,*)AIJX,AIJX1
   42		continue
        END DO
        
        
   
   901  ETAMAX=SQRT(2./PI)*FR**2*SQRT(AIJX**2+AIJX1**2)*SQRT((1.+QDD**2)/PHIPP)
        !ETAMAX=SQRT(2./PI)*FR**2*SQRT(AIJX1**2)*SQRT((1.+QDD**2)/PHIPP)
   
        WRITE(9001,*)PSIM/PI*180.,ETAMAX
        !WRITE(*,*)PSIM/PI*180.,AIJX,AIJX1
        !WRITE(*,*)AIJX,AIJX1
    END DO
    CLOSE(9001)
    
	DEALLOCATE (CORL,CNL)
    RETURN
END SUBROUTINE