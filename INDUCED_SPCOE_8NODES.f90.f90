SUBROUTINE cofah_SP8(cord,conm,ps,aij,hij) 
!c    ***********************************************      
!c     calculate the aij and hij terms in the matrix    4---7---3  
!c      cord(8,3) the cord. of 8 nodes (clockwise)      |       |
!c      conm(8,3) the normal vector of 8 nodes          8       6
!c      ps(3) the cord. of fiele point                  |       |
!c      aij(8): panel source effect                     1---5---2
!c      hij(8): panel dipole effect
!c   ************************************************
!cdp   implicit REAL(a-h,o-z)
      REAL:: cord(8,3),ps(3),conm(8,3),aij(8),hij(8),pr(3)
      REAL:: r(4),srr(3),ay(8),hy(8)
      REAL:: sn(8),snx(8),sne(8),sr(3),dn(3)
      REAL:: AJ,XI,ETA
      common nml,nbp,nfp,nbe,nfe,np,mp,nf,radius
      common /gwa/ga(400),gw(400)
   
      r1=0.
      do 5 i=1,8
      aij(i)=0.
    5 hij(i)=0.
      
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
      k=ngus*(ngus-1)/2

!c     gauss-legendre numerical integration
     
      do 40 ix=1,ngus
      xi=ga(k+ix-1)
      do 40 iy=1,ngus
      eta=ga(k+iy-1)
      
      !call isopar_8(xi,eta,cord,sn,snx,sne,aj,sr,dn)
      
      rh=0.
      do 45 i=1,3
      srr(i)=sr(i)-ps(i)   
   45 rh=rh+srr(i)*srr(i)
      rh=sqrt(rh)
!C LEE ORIGIN
	  RH3=RH**3
      rh3=0.0-rh**3
!c     
!c     nml=1 using the normal vector accurately on isopar. curve surface
!c      else using the isoparameter normal vector by 8 node nomal
!c      
      if(nml.ne.1) then
       do 48 i1=1,3  
       dn(i1)=0.
       do 48 i2=1,8	
   48  dn(i1)=dn(i1)+conm(i2,i1)*sn(i2)
      endif
!c       
      bnr=srr(1)*dn(1)+srr(2)*dn(2)+srr(3)*dn(3)
!c      
      do 60 j=1,8
      ay(j)=gw(k+iy-1)*sn(j)*aj/rh
   60 hy(j)=gw(k+iy-1)*sn(j)*bnr*aj/rh3
!c
      do 40 j=1,8
      aij(j)=aij(j)+gw(k+ix-1)*ay(j)
   40 hij(j)=hij(j)+gw(k+ix-1)*hy(j)
   
      return
END SUBROUTINE

SUBROUTINE cofahs_SP8(cord,conm,ps,aij,hij)
!c    ***********************************************      
!c     calculate the aij and hij terms in the matrix   
!c      the ps is one of the nodes of element          4---7---3
!c      cord(8,3) the cord. of 8 nodes (clockwise)     |       |
!c      conm(8,3) the normal vector of 8 nodes         8       6
!c      ps(3) the cord. of fiele point                 |       |
!c      aij(8): panel source effect                    1---5---2
!c      hij(8): panel dipole effect
!c   ************************************************
!cdp   implicit REAL(a-h,o-z)
      REAL:: XI1,XI2,XI3,XI4,ETA1,ETA2,ETA3,ETA4
      REAL:: cord(8,3),conm(8,3),ps(3),aij(8),hij(8)
      REAL:: aij1(8),hij1(8),cord1(8,3),conm1(8,3)
      REAL:: ay(8),hy(8),rh(4),rh3(4),bnr(4),aj(4)
      REAL:: snx(8),sne(8),sm(8,4)
      REAL:: sn1(8),sn2(8),sn3(8),sn(8,4),sn4(8)
      equivalence (sn(1,1),sn1),(sn(1,2),sn2),(sn(1,3),sn3),(sn(1,4),sn4)
      REAL:: sr1(3),sr2(3),sr3(3),sr(3,4),sr4(3)
      equivalence (sr(1,1),sr1),(sr(1,2),sr2),(sr(1,3),sr3),(sr(1,4),sr4)
      REAL:: dn1(3),dn2(3),dn3(3),dn(3,4),dn4(3)
      equivalence (dn(1,1),dn1),(dn(1,2),dn2),(dn(1,3),dn3),(dn(1,4),dn4)
      common nml,nbp,nfp,nbe,nfe,np,mp,nf,radius
      common /gwa/ga(400),gw(400)

      D1=0.9
      D2=D1
      D3=D2
      D4=D3

      r1=0.
      isn=0
      do 5 i=1,8
      t=(cord(i,1)-ps(1))**2+(cord(i,2)-ps(2))**2+(cord(i,3)-ps(3))**2
      t=sqrt(t)
      if(t.lt.1.e-6) isn=i
      aij1(i)=0.
    5 hij1(i)=0.
!c    
      if(isn.eq.0) then
        pause 'warning: isn=0, the source point out of the element'
        stop
      endif	
!      
      ngus=9
!   
      kgus=ngus*(ngus-1)/2
!
!     gauss-legendre numerical integration
   
    do 235 i=1,8
      do 235 j=1,3
		cord1(I,j)=cord(I,j)
		conm1(I,j)=conm(I,j)    
 235	CONTINUE   
         
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

            !CALL ISOPAR_8(XI1,ETA1,CORD,SN1,SNX,SNE,AJ(1),SR1,DN(:,1))
            !CALL ISOPAR_8(XI2,ETA2,CORD,SN2,SNX,SNE,AJ(2),SR2,DN(:,2))
            !CALL ISOPAR_8(XI3,ETA3,CORD,SN3,SNX,SNE,AJ(3),SR3,DN(:,3))
	        !CALL ISOPAR_8(XI4,ETA4,CORD,SN4,SNX,SNE,AJ(4),SR4,DN(:,4))     

            do 245 j=1,4
                rh(j)=0.
                do 246 i=1,3
                246 rh(j)=rh(j)+(sr(i,j)-ps(i))*(sr(i,j)-ps(i))
                rh(j)=sqrt(rh(j))
	            RH3(J)=0.0-RH(J)**3
        245	CONTINUE
   
!     nml=1 using the normal vector accurately on isopar. curve surface
!     else using the isoparameter normal vector by 8 node nomal
!      
            if(nml.ne.1) then
                do 250 j=1,4
                    do 250 i1=1,3 
                        dn(i1,j)=0.
                        do 250 i2=1,8	
                250  dn(i1,j)=dn(i1,j)+conm1(i2,i1)*sn(i2,i1)
            endif

            do 255 j=1,4      
                bnr(j)=0.
                do 255 i=1,3
                255 bnr(j)=bnr(j)+(sr(i,j)-ps(i))*dn(i,j)
     
                do 260 i=1,8
                    ay(i)=0.
                260 hy(i)=0.
         
                DO j=1,8
                    DO k=1,4
                    IF(K.EQ.1) AM=ABS(-1.-DXI)
                    IF(K.EQ.2) AM=ABS(-1.-DET)
                    IF(K.EQ.3) AM=ABS(1.-DXI)
                    IF(K.EQ.4) AM=ABS(1.-DET)

                    AY(J)=AY(J)+GW(KGUS+IY-1)*SN(J,K)*AJ(K)*(p1/rh(k))*2.*AM/4.

                    END DO
                END DO

            DO J=1,8
                AIJ1(J)=AIJ1(J)+GW(KGUS+IX-1)*AY(J)
                HIJ1(J)=0.
            END DO
        END DO
    END DO 
      
    DO J=1,8
      	AIJ(J)=AIJ1(J)
        HIJ(J)=HIJ1(J)
    END DO	

      return
END SUBROUTINE
       