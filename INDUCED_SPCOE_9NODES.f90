SUBROUTINE cofah_SP9(cord,conm,ps,aij,hij) 
!c    ***********************************************      
!c     calculate the aij and hij terms in the matrix    7---6---5  
!c      cord(8,3) the cord. of 8 nodes (clockwise)      |       |
!c      conm(8,3) the normal vector of 8 nodes          8   9   4
!c      ps(3) the cord. of fiele point                  |       |
!c      aij(8): panel source effect                     1---2---3
!c      hij(8): panel dipole effect
!c   ************************************************
!cdp   implicit REAL(a-h,o-z)
      REAL:: cord(9,3),ps(3),conm(9,3),aij(9),hij(9),pr(3)
      REAL:: r(4),srr(3),ay(9),hy(9)
      REAL:: sn(9),snx(9),sne(9),sr(3),dn(3)
      REAL:: AJ,XI,ETA
      common nml,nbp,nfp,nbe,nfe,np,mp,nf,radius
      common /gwa/ga(400),gw(400)
!c      
      r1=0.
      do 5 i=1,9
      aij(i)=0.
    5 hij(i)=0.
!c      
      do 10 i=1,3
      pr(i)=0.
      do 12 j=1,9
   12 pr(i)=pr(i)+cord(j,i)
      pr(i)=pr(i)/9.-ps(i)
      r1=r1+pr(i)*pr(i)
   10 continue
      r1=sqrt(r1)      

      do 15 i=1,4
   15 r(i)=0.
      do 20 i=1,3
      R(1)=R(1)+(CORD(1,I)-CORD(3,I))**2
      R(2)=R(2)+(CORD(3,I)-CORD(5,I))**2
      R(3)=R(3)+(CORD(5,I)-CORD(7,I))**2
      R(4)=R(4)+(CORD(7,I)-CORD(1,I))**2
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
   
      k=ngus*(ngus-1)/2

!c     gauss-legendre numerical integration
     
      do 40 ix=1,ngus
      xi=ga(k+ix-1)
      do 40 iy=1,ngus
      eta=ga(k+iy-1)
      
      call isopar_9(xi,eta,cord,sn,snx,sne,aj,sr,dn)
      
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
       do 48 i2=1,9	
   48  dn(i1)=dn(i1)+conm(i2,i1)*sn(i2)
      endif
!c       
      bnr=srr(1)*dn(1)+srr(2)*dn(2)+srr(3)*dn(3)
!c      
      do 60 j=1,9
      ay(j)=gw(k+iy-1)*sn(j)*aj/rh
   60 hy(j)=gw(k+iy-1)*sn(j)*bnr*aj/rh3
!c
      do 40 j=1,9
      aij(j)=aij(j)+gw(k+ix-1)*ay(j)
   40 hij(j)=hij(j)+gw(k+ix-1)*hy(j)
   
      return
END SUBROUTINE

SUBROUTINE cofahs_SP9(cord,conm,ps,aij,hij)
!c    ***********************************************      
!c     calculate the aij and hij terms in the matrix   
!c      the ps is one of the nodes of element          7---6---5
!c      cord(8,3) the cord. of 8 nodes (clockwise)     |       |
!c      conm(8,3) the normal vector of 8 nodes         8   9   4
!c      ps(3) the cord. of fiele point                 |       |
!c      aij(8): panel source effect                    1---2---3
!c      hij(8): panel dipole effect
!c   ************************************************
!cdp   implicit REAL(a-h,o-z)
      REAL:: XI1,XI2,XI3,XI4,ETA1,ETA2,ETA3,ETA4
      REAL:: cord(9,3),conm(9,3),ps(3),aij(9),hij(9)
      REAL:: aij1(9),hij1(9),cord1(9,3),conm1(9,3)
      REAL:: ay(9),hy(9),rh(4),rh3(4),bnr(4),aj(4)
      REAL:: snx(9),sne(9),sm(9,4)
      REAL:: sn1(9),sn2(9),sn3(9),sn(9,4),sn4(9)
      equivalence (sn(1,1),sn1),(sn(1,2),sn2),(sn(1,3),sn3),(sn(1,4),sn4)
      REAL:: sr1(3),sr2(3),sr3(3),sr(3,4),sr4(3)
      equivalence (sr(1,1),sr1),(sr(1,2),sr2),(sr(1,3),sr3),(sr(1,4),sr4)
      REAL:: dn1(3),dn2(3),dn3(3),dn(3,4),dn4(3)
      equivalence (dn(1,1),dn1),(dn(1,2),dn2),(dn(1,3),dn3),(dn(1,4),dn4)
      common nml,nbp,nfp,nbe,nfe,np,mp,nf,radius
      common /gwa/ga(400),gw(400)
!c      
      r1=0.
      isn=0
      do 5 i=1,9
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
!c      
      ngus=11
!c   
      kgus=ngus*(ngus-1)/2
!c
!c     gauss-legendre numerical integration
!c     
!c     corner node 1, 3, 5, 7
!c  ********************************
      if((mod(isn,2).eq.1).and.(isn.ne.9)) then
      
      do 30 i=1,8
      k=mod(isn+i-2,8)+1
      do 30 j=1,3
      cord1(i,j)=cord(k,j)
   30 conm1(i,j)=conm(k,j)
	do 31 j=1,3
		cord1(9,j)=cord(9,j)
		conm1(9,j)=conm(9,j)
31	continue
	        
      do 40 ix=1,ngus
      p1=.5*(1.0+ga(kgus+ix-1))
      do 40 iy=1,ngus
      p2=.5*(1.0+ga(kgus+iy-1))
     
      xi1=-1.+2.*p1
      eta1=-1.+2.*p1*p2
      xi2=-1.+2.*p1-2.*p1*p2
      eta2=-1.+2.*p1 
      call isopar_9(xi1,eta1,cord1,sn1,snx,sne,aj(1),sr1,dn1)
      call isopar_9(xi2,eta2,cord1,sn2,snx,sne,aj(2),sr2,dn2)
      do 42 j=1,2
      do 42 k=1,9
      sm(k,j)=sn(k,j)
      if(k.eq.1) sm(k,j)=sm(k,j)-1.
   42 continue      

      do 45 j=1,2
      rh(j)=0.
      do 46 i=1,3
   46 rh(j)=rh(j)+(sr(i,j)-ps(i))*(sr(i,j)-ps(i))
      rh(j)=sqrt(rh(j))
!CLEE ORIGIN
!C	RH3(J)=RH(J)**3
      rh3(j)=0.0-rh(j)**3
45	CONTINUE
    
!c     nml=1 using the normal vector accurately on isopar. curve surface
!c      else using the isoparameter normal vector by 8 node nomal
!c      
      if(nml.ne.1) then
       do 50 j=1,2
       do 50 i1=1,3  
       dn(i1,j)=0.
       do 50 i2=1,9	
   50  dn(i1,j)=dn(i1,j)+conm1(i2,i1)*sn(i2,i1)
      endif

      do 55 j=1,2       
      bnr(j)=0.
      do 55 i=1,3
   55 bnr(j)=bnr(j)+(sr(i,j)-ps(i))*dn(i,j)
      
      do 60 i=1,9
      ay(i)=0.
   60 hy(i)=0.
         
      do 70 j=1,9
      do 70 k=1,2
      ay(j)=ay(j)+gw(kgus+iy-1)*sn(j,k)*aj(k)*(p1/rh(k))
!CLEE ORIGIN
!C    w=(sm(j,k)/rh(k))*(bnr(k)/rh(k))*aj(k)*(p1/rh(k))
      w=sm(j,k)*bnr(k)*aj(k)*p1/rh3(k)

   70 hy(j)=hy(j)+gw(kgus+iy-1)*w

      do 40 j=1,9
      aij1(j)=aij1(j)+gw(kgus+ix-1)*ay(j)
   40 hij1(j)=hij1(j)+gw(kgus+ix-1)*hy(j)

!c     assemble the local coef. matrix
 
      do 80 j=1,8
      k=mod(isn+j-2,8)+1
      aij(k)=aij1(j)
   80 hij(k)=hij1(j)  
	aij(9)=aij1(9)
	hij(9)=hij1(9) 

!c     midside node 2, 4, 6, 8    
!c    ***************************
      else if(isn.ne.9)then 
!c      
      do 135 i=1,8
      k=mod(isn+i-2,8)+1
      k1=mod(i,8)+1
      do 135 j=1,3
      cord1(k1,j)=cord(k,j)
  135 conm1(k1,j)=conm(k,j)    
	do 136 j=1,3
		cord1(9,j)=cord(9,j)
		conm1(9,j)=conm(9,j)
136	continue
            
      do 140 ix=1,ngus
      p1=.5*(1.0+ga(kgus+ix-1))
      do 140 iy=1,ngus
      p2=.5*(1.0+ga(kgus+iy-1))
     
      xi1=p1
      eta1=-1.+2.*p1*p2
      xi2=p1-2.*p1*p2
      eta2=-1.+2.*p1 
      xi3=-p1
      eta3=-1.+2.*p1-2.*p1*p2
      call isopar_9(xi1,eta1,cord1,sn1,snx,sne,aj(1),sr1,dn1)
      call isopar_9(xi2,eta2,cord1,sn2,snx,sne,aj(2),sr2,dn2)
      call isopar_9(xi3,eta3,cord1,sn3,snx,sne,aj(3),sr3,dn3)
      do 142 j=1,3
      do 142 k=1,9
      sm(k,j)=sn(k,j)
      if(k.eq.2) sm(k,j)=sm(k,j)-1.
  142 continue      

      do 145 j=1,3
      rh(j)=0.
      do 146 i=1,3
  146 rh(j)=rh(j)+(sr(i,j)-ps(i))*(sr(i,j)-ps(i))
      rh(j)=sqrt(rh(j))
!CLEE ORIGIN
!C      rh3(j)=rh(j)**3
	RH3(J)=0.0-RH(J)**3
145	CONTINUE
     
!c     nml=1 using the normal vector accurately on isopar. curve surface
!c      else using the isoparameter normal vector by 8 node nomal
!c      
      if(nml.ne.1) then
       do 150 j=1,3
       do 150 i1=1,3  
       dn(i1,j)=0.
       do 150 i2=1,9	
  150  dn(i1,j)=dn(i1,j)+conm1(i2,i1)*sn(i2,i1)
      endif

      do 155 j=1,3       
      bnr(j)=0.
      do 155 i=1,3
  155 bnr(j)=bnr(j)+(sr(i,j)-ps(i))*dn(i,j)
      
      do 160 i=1,9
      ay(i)=0.
  160 hy(i)=0.
         
      do 170 j=1,9
      do 170 k=1,3
      cf=1.
      if(k.ne.2) cf=.5
      ay(j)=ay(j)+gw(kgus+iy-1)*sn(j,k)*aj(k)*(p1/rh(k))*cf
!CLEE ORIGIN 
!C      w=(sm(j,k)/rh(k))*(bnr(k)/rh(k))*aj(k)*(p1/rh(k))
      w=sm(j,k)*bnr(k)*aj(k)*p1/rh3(k)

  170 hy(j)=hy(j)+gw(kgus+iy-1)*w*cf

      do 140 j=1,9
      aij1(j)=aij1(j)+gw(kgus+ix-1)*ay(j)
  140 hij1(j)=hij1(j)+gw(kgus+ix-1)*hy(j)

!c     assemble the local coef. matrix
   
      do 180 j=1,8
      k=mod(isn+j-2,8)+1
      k1=mod(j,8)+1
      aij(k)=aij1(k1)
  180 hij(k)=hij1(k1)
	aij(9)=aij1(9)
	hij(9)=hij1(9)
	
	else
   
    do 235 i=1,9
      do 235 j=1,3
		cord1(I,j)=cord(I,j)
		conm1(I,j)=conm(I,j)    
 235	CONTINUE   
         
      do 240 ix=1,ngus
      p1=.5*(1.0+ga(kgus+ix-1))
      do 240 iy=1,ngus
      p2=.5*(1.0+ga(kgus+iy-1))
    
      xi1=-p1+2.*p1*p2
      eta1=-p1
      xi2=p1
      eta2=-p1+2.*p1*p2 
      xi3=p1-2.*p1*p2
      eta3=p1
	xi4=-p1
	eta4=p1-2*p1*p2
      call isopar_9(xi1,eta1,cord1,sn1,snx,sne,aj(1),sr1,dn1)
      call isopar_9(xi2,eta2,cord1,sn2,snx,sne,aj(2),sr2,dn2)
      call isopar_9(xi3,eta3,cord1,sn3,snx,sne,aj(3),sr3,dn3)
	  call isopar_9(xi4,eta4,cord1,sn4,snx,sne,aj(4),sr4,dn4)
      do 242 j=1,4
      do 242 k=1,9
      sm(k,j)=sn(k,j)
      if(k.eq.9) sm(k,j)=sm(k,j)-1.
  242 continue      

      do 245 j=1,4
      rh(j)=0.
      do 246 i=1,3
  246 rh(j)=rh(j)+(sr(i,j)-ps(i))*(sr(i,j)-ps(i))
      rh(j)=sqrt(rh(j))
!CLEE ORIGIN
!C      rh3(j)=rh(j)**3
	RH3(J)=0.0-RH(J)**3
245	CONTINUE
   
!c     nml=1 using the normal vector accurately on isopar. curve surface
!c      else using the isoparameter normal vector by 8 node nomal
!c      
      if(nml.ne.1) then
       do 250 j=1,4
       do 250 i1=1,3 
       dn(i1,j)=0.
       do 250 i2=1,9	
  250  dn(i1,j)=dn(i1,j)+conm1(i2,i1)*sn(i2,i1)
      endif

      do 255 j=1,4      
      bnr(j)=0.
      do 255 i=1,3
  255 bnr(j)=bnr(j)+(sr(i,j)-ps(i))*dn(i,j)
     
      do 260 i=1,9
      ay(i)=0.
  260 hy(i)=0.
         
      do 270 j=1,9
      do 270 k=1,4
!c      cf=1.
!c      if(k.ne.2) cf=.5
	cf=0.5
      ay(j)=ay(j)+gw(kgus+iy-1)*sn(j,k)*aj(k)*(p1/rh(k))*cf
!CLEE ORIGIN 
!C      w=(sm(j,k)/rh(k))*(bnr(k)/rh(k))*aj(k)*(p1/rh(k))
      w=sm(j,k)*bnr(k)*aj(k)*p1/rh3(k)

  270 hy(j)=hy(j)+gw(kgus+iy-1)*w*cf
!c
      do 240 j=1,9
      aij1(j)=aij1(j)+gw(kgus+ix-1)*ay(j)
  240 hij1(j)=hij1(j)+gw(kgus+ix-1)*hy(j)
!c
!c     assemble the local coef. matrix
!c   

!C      do 280 j=1,9
!c      k=mod(isn+j-2,8)+1
!c      k1=mod(j,8)+1
!C      aij(K)=aij1(k1)
!C  280 hij(k)=hij1(k1)
      do 280 j=1,9
      aij(J)=aij1(J)
      hij(J)=hij1(J)	         
280	CONTINUE	
  
      endif
      
      return
END SUBROUTINE
       