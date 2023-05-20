!    Last change:  PA   05 June 2022    3:00 am
!    PROGRAM FOR COMPRESSIBLE FLOW IN CURVILINEAR CORDINATE SYSTEM
!    This program is used to validate with kakade by changing the orientation of cylinder
!    This program uses BC suitable for high heating 
program compressible
implicit none

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: FC,GC,FD,GD,U_VEC_OLD,S
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: U_VEC,F,G,H,VAR,V_CONT,q_vec
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: P,E,ROW,VIS,COND
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: ZI_X,ZI_Y,ETA_X,ETA_Y,X,Y,SP
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: NS_X,NS_Y,TS_X,TS_Y,TA_X,TA_Y,NA_X,NA_Y,A,TV,p_mean

DOUBLE PRECISION    :: LAMBDA(2),PHI(2),RN(4),VT(4),VN(4),SI(2),G_ETA(2),FX_I(4),GY_I(4)
DOUBLE PRECISION    :: DVT_DN,FV,PHIF_R,SP_EN,SONIC,DROP,etol,SUM1,DN,CORR,V_NOLD
DOUBLE PRECISION    :: Tinf,uinf,vinf,Pinf,deninf,Einf,alpha,TERM,TIME,AAA,BBB
DOUBLE PRECISION    :: R_MAX,R_LOCAL,M_LOCAL,MN,MAXM,MINM,RPLUSU,RMINUSU
DOUBLE PRECISION    :: XG,YG,X_ZIG,X_ETAG,Y_ZIG,Y_ETAG,GRADZI_T
DOUBLE PRECISION    :: GRAD_ETA,COEFF_ZI,COEFF_ETA,COEFF_ZI_ETA,DF_ZI,DF_ETA,DG_ZI,DG_ETA
DOUBLE PRECISION    :: DU_DX,DU_DY,DV_DX,DV_DY,DT_DX,DT_DY,DF_DX,DG_DY
DOUBLE PRECISION    :: DU_DZI,DV_DZI,DVN_DZI,DU_DETA,DV_DETA,DT_DZI,DT_DETA
DOUBLE PRECISION    :: X_ZI,X_ETA,Y_ZI,Y_ETA,ZI,ETA,DZI_X_ZI,DZI_Y_ZI,DETA_X_ETA,DETA_Y_ETA
DOUBLE PRECISION    :: U_INT,V_INT,P_INT,VIS_INT,COND_INT,V_NINF,V_N,V_T,RNP,RNM,MEAN
DOUBLE PRECISION    :: WPPU,WPU,WU,WF,FF,Z1,Z2,U_SN,URU1,URU2,URU,URC,UR1,UR2,URU_L,tstart,tend
DOUBLE PRECISION    :: ZI_XA,ZI_YA,ETA_XA,ETA_YA,FDN,GDN,FDC,GDC,FCN,GCN,FCC,GCC,QQ,HH

DOUBLE PRECISION    :: CX,CY,CXP,CYP,CXV,CYV,CXD,CYD,CL,CD,CDP,CDV,CDD,ANUSS,FUN1,FUN2
DOUBLE PRECISION    :: QFLUX,QFLUXP,JACOB,JACOBP,DILAT,DILATP,GRAETA,GRAETAP,VORT,VORTP,dev,dev_max
DOUBLE PRECISION    :: PX,PY,DILATX,DILATY,VORTX,VORTY,HEATF,SHEAR,Nu_local,ML,p_min,p_max
DOUBLE PRECISION    :: PXY,VORTXY,DILATXY,CXYP,CXYV,CXYD,CMP,CMV,CMD,CM,XX_MAX,UTP,UTPP,VAR1,VAR2,VAR3

!---------------------------PVM+H-----------------------------------------------------------
! DOUBLE PRECISION    :: AW,BW,ML(3),MLP(3),MLM(3),FP,FM,FPMM,FMMP,WH,WHP,WHM,WHA,PT,PTP,PTM
! DOUBLE PRECISION    :: Piph,Pimh,DP_DZI(4),DP_DETA(4),DP_DZ(4)
!-------------------------------------------------------------------------------------------

INTEGER             :: I,J,L,W_L,IP,IPP,IM,IOS
INTEGER             :: nstep,maxstep,st,np(2)
INTEGER             :: IC1,IC2,IC3,IC4,ITE,ITU1,ITU2,ITU3,ITU4,ITL1,ITL2,ITL3,ITL4

!--------------------------------------------------------------------------------------
 CHARACTER(LEN=21) ::  MODEL="VISCOUS"
!CHARACTER(LEN=100) :: OUTPUT
 CHARACTER(LEN=20)  :: TVD_E="no"
 CHARACTER(LEN=20)  :: ENTH_TOT="no"
 CHARACTER(LEN=20)  :: VEL_ETA="mean"
 LOGICAL           ::  LOGIC1,LOGIC2
 LOGICAL           ::  INITIAL = .FALSE. 

!================================PARAMETER=============================================
DOUBLE PRECISION,PARAMETER :: THETA=0.0d0,PHI_R=-7.0d0,S_=110.0d0,GAMA=1.40d0,R=287.0d0,U_inf=1.0d0
DOUBLE PRECISION,PARAMETER :: DT=0.00001d0,gw=0.0d0,EPS1=0.0d0,TS=1.0d0+EPS1,T0=273.0d0,T_FS=300.0d0	
DOUBLE PRECISION,PARAMETER :: A1=0.028109,A2=1.0739,A3=-0.099183,C1=0.012012,Cv=R/(GAMA-1)
DOUBLE PRECISION,PARAMETER :: C2=0.065283,C3=-0.015755,MACH=0.1d0,Cp=GAMA*Cv,Tr=T_FS/T0 
DOUBLE PRECISION,PARAMETER :: RE=1000.0d0,PR=0.71d0,GR_IN=1.0d0,PI=acos(-1.0d0)
DOUBLE PRECISION,PARAMETER :: APH1=C1/2.0,APH2=C2/3.0,APH3=C3/4.0,FR =SQRT(1.0d0/GR_IN)
DOUBLE PRECISION,PARAMETER :: EK=0.5*MACH*MACH*GAMA*(GAMA-1.0D0),PK=1.0D0/(GAMA*MACH*MACH)
DOUBLE PRECISION,PARAMETER :: BK=GAMA*(GAMA-1.0)*MACH*MACH/RE,AA=0.2d0,PP=1.2d0
DOUBLE PRECISION,PARAMETER :: FFACT=PP*(AA**(PP-1.0)),THETAR=THETA*PI/180
DOUBLE PRECISION,PARAMETER :: B1=1.0d0,B2=5.47886E-03,B3=2.33633E-02,B4=7.12231E-03,B5=6.89677E-04
!---------------------------------------------------------------------------------------   
    
    maxstep=40000000  ! MAXIMUM TIME LOOPS

    nstep=0

!=================READ INPUT FILE===============================

 OPEN(2,FILE='INP.DAT',STATUS='UNKNOWN',IOSTAT=ios)
 READ(2,*)np(1),np(2),ZI,ETA
 READ(2,*)IC1,IC2,IC3,IC4
 CLOSE(2)

!----------------ALLOCATION-------------------------
 ALLOCATE(X(np(1),np(2)),Y(np(1),np(2)),ZI_X(np(1),np(2)),ZI_Y(np(1),np(2)))
 ALLOCATE(ETA_X(np(1),np(2)),ETA_Y(np(1),np(2)))
!-----------------------------------------
  
 OPEN(2,FILE='INP.DAT',STATUS='UNKNOWN',IOSTAT=ios)
 READ(2,*)
 READ(2,*)
       
 DO J=1,np(2)
 DO I=1,np(1)
 READ(2,*)AAA,BBB,X(I,J),Y(I,J)
 END DO
 END DO

 DO J=1,np(2)
 DO I=1,np(1)
 READ(2,*)ZI_X(I,J),ZI_Y(I,J),ETA_X(I,J),ETA_Y(I,J)
 END DO
 END DO

 CLOSE(2)

!  Finding the Trailing edge I-Index
! =========================================================================
      XX_MAX=X(1,1)
      ITE=1
      DO I=2,np(1)-1
      IF(X(I,1)>XX_MAX)THEN
      ITE=I
      XX_MAX=X(I,1) 
      ENDIF
      END DO

      IF(ITE==1)THEN
       ITU1=2
       ITU2=3
       ITU3=4
       ITU4=5
       ITL1=np(1)-1
       ITL2=np(1)-2
       ITL3=np(1)-3
       ITL4=np(1)-4
      ELSE
       ITU1=ITE-1
       ITU2=ITE-2
       ITU3=ITE-3
       ITU4=ITE-4
       ITL1=ITE+1
       ITL2=ITE+2
       ITL3=ITE+3
       ITL4=ITE+4
      ENDIF

!      write(*,*)ITE
!      PAUSE
!========================================================================== 
 

!================TRUNCATION================================================

 !  np(2)=309        ! artificial boundary cut at 60 non-dimensional distance
   
!==================ALLOCATION==============================================
   
 ALLOCATE(FC(np(1),np(2),4),GC(np(1),np(2),4))
 ALLOCATE(FD(np(1),np(2),4),GD(np(1),np(2),4))
 ALLOCATE(U_VEC_OLD(np(1),np(2),4),S(np(1),np(2),4))

 ALLOCATE(U_VEC(np(1),np(2),4),VAR(np(1),np(2),3))
 ALLOCATE(V_CONT(np(1),np(2),2))
 ALLOCATE(q_vec(np(1),np(2),4))

 ALLOCATE(P(np(1),np(2)),E(np(1),np(2)),ROW(np(1),np(2)))
 ALLOCATE(VIS(np(1),np(2)),COND(np(1),np(2)))

 ALLOCATE(F(np(1),3,2),G(np(1),3,2))

 ALLOCATE(NS_X(np(1)),NS_Y(np(1)),TS_X(np(1)),TS_Y(np(1)),TA_X(np(1)),TA_Y(np(1)))
 ALLOCATE(NA_X(np(1)),NA_Y(np(1)),A(np(1)))

 ALLOCATE(SP(np(1),2))
 ALLOCATE(p_mean(np(1)),TV(np(1)))

!===========Rotating the cylinder with mesh===========================	

!     READING METRICS FILE ALREADY GENRATED BY ANOTHER PROGRAM
!     ZI_X=DERIVATIVE OF "ZI" W.R.T X
!     ZI_Y==DERIVATIVE OF "ZI" W.R.T  Y
!     ETA_X==DERIVATIVE OF "ETA" W.R.T X
!     ETA_Y==DERIVATIVE OF "ETA" W.R.T Y

    DO J=1,np(2)
    DO I=1,np(1)

    XG=X(I,J)
    YG=Y(I,J)
    X(I,J)=XG*COS(PHI_R*PI/180)-YG*SIN(PHI_R*PI/180)
    Y(I,J)=XG*SIN(PHI_R*PI/180)+YG*COS(PHI_R*PI/180)  

    JACOB=ZI_X(I,J)*ETA_Y(I,J)-ZI_Y(I,J)*ETA_X(I,J)
    X_ZIG=ETA_Y(I,J)/JACOB
    X_ETAG=-ZI_Y(I,J)/JACOB
    Y_ZIG=-ETA_X(I,J)/JACOB
    Y_ETAG=ZI_X(I,J)/JACOB

    X_ZI=X_ZIG*COS(PHI_R*PI/180)-Y_ZIG*SIN(PHI_R*PI/180)
    X_ETA=X_ETAG*COS(PHI_R*PI/180)-Y_ETAG*SIN(PHI_R*PI/180)
    Y_ZI=X_ZIG*SIN(PHI_R*PI/180)+Y_ZIG*COS(PHI_R*PI/180)
    Y_ETA=X_ETAG*SIN(PHI_R*PI/180)+Y_ETAG*COS(PHI_R*PI/180)
        		
    JACOB=X_ZI*Y_ETA-X_ETA*Y_ZI
    ZI_X(I,J)=Y_ETA/JACOB
    ETA_X(I,J)=-Y_ZI/JACOB
    ZI_Y(I,J)=-X_ETA/JACOB
    ETA_Y(I,J)=X_ZI/JACOB
    END DO
    END DO 

!===========DETERMINING R_MAX======================================

    R_MAX=0.d0
    DO J=1,np(2)
    DO I=1,np(1)-1
    R_LOCAL=SQRT(X(I,J)**2+Y(I,J)**2)
    IF(R_LOCAL.GE.R_MAX)R_MAX=R_LOCAL 
    END DO
    END DO
	
!========Calculating the UNIT normal and tangential vectors at the solid boundary

    DO I=1,np(1)
    GRAD_ETA=SQRT((ETA_X(I,1)**2) + (ETA_Y(I,1)**2))
    NS_X(I)=ETA_X(I,1)/GRAD_ETA
    NS_Y(I)=ETA_Y(I,1)/GRAD_ETA
    TS_X(I)=-NS_Y(I)
    TS_Y(I)=NS_X(I)
    END DO

!========CALCULATING UNIT NORMAL (INTO THE FLOW DOMAIN) AT THE ARTIFICIAL BOUNDARY

    DO I=1,np(1)
    GRAD_ETA=SQRT(ETA_X(I,np(2))**2+ETA_Y(I,np(2))**2)
    NA_X(I)=-ETA_X(I,np(2))/GRAD_ETA
    NA_Y(I)=-ETA_Y(I,np(2))/GRAD_ETA

    TA_X(I)=-NA_Y(I)
    TA_Y(I)=NA_X(I)
    END DO

!==========free-stream conditions (dimensionless)=======================

    alpha=theta*PI/180
    uinf=U_inf*COS(alpha)
    vinf=U_inf*SIN(alpha)
    Pinf=0.d0
    Tinf=1.d0
    deninf=1.d0
    Einf=Tinf+APH1*(Tinf-1)**2+APH2*(Tinf-1)**3+APH3*(Tinf-1)**4+EK*(uinf**2+vinf**2)
    
!===============SETTING INITIAL CONDITIONS===================

    IF(INITIAL) THEN
    time=0.d0
    
    DO J=1,np(2)  !bs,es 
    DO I=1,np(1) 

    VAR(I,J,1)=uinf
    VAR(I,J,2)=vinf
    VAR(I,J,3)=Tinf
    ROW(I,J)=1.d0     
    E(I,J)=Einf
    P(I,J)=Pinf 
    
  IF(MODEL=='VISCOUS') THEN   
    VIS(I,J)=1.d0
    COND(I,J)=1.d0
  ELSE
    VIS(I,J)=0.d0
    COND(I,J)=0.d0
  ENDIF
    END DO
    END DO
      
  ELSE

!======================Read OUTPUT File=====================
 OPEN(unit=7,file='OUTPUT',STATUS='unknown')
 
 READ(7,100) Time,nstep
 READ(7,*) 
 READ(7,*)

 DO J=1,np(2)
 READ(7,*)
 DO I=1,np(1)
 READ(7,*,iostat=st) X(I,J),Y(I,J),ROW(I,J),P(I,J),VAR(I,J,3),VAR(I,J,1),VAR(I,J,2),E(I,J)
 IF(st/=0) THEN
 WRITE(6,'(A,4(1X,I4))') "Error in reading combined input file at",I,J,st
 stop
 ENDIF

 END DO	    
 END DO

 CLOSE(7)

 100	format(8X,F20.10,I9,1X)


  GO TO 101
!======================================================================= 
 END IF

!----------------------------------------------------------------------------------------
!       SETTING BOUNDARY CONDITIONS AT SOLID SURFACE
!----------------------------------------------------------------------------------------
 

   IF(MODEL=='VISCOUS') THEN
    J=1
    DO I=1,np(1)-1

    VAR(I,J,1)=0.0d0
    VAR(I,J,2)=0.0d0 
    VAR(I,J,3)=TS                       
    E(I,J)=VAR(I,J,3)+APH1*(VAR(I,J,3)-1)**2+APH2*(VAR(I,J,3)-1)**3 &
             +APH3*(VAR(I,J,3)-1)**4+EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)
    VIS(I,J)=(VAR(I,J,3)**1.5)*(Tr+(S_/T0))/(VAR(I,J,3)*Tr+(S_/T0))
    COND(I,J)=A1+A2*VAR(I,J,3)+A3*VAR(I,J,3)**2   
    END DO  ! I LOOP ENDS
!-------------------------PERIODICITY-----------------------------------
    VAR(np(1),1,1:3)=VAR(1,1,1:3)
    E(np(1),1)=E(1,1)

    VIS(np(1),1)=VIS(1,1)
    COND(np(1),1)=COND(1,1)
               
!================PRESSURE BOUNDARY CONDITION============================

!CALCULATING FLUX VECTORS AND SOURCE VECTOR (WITHOUT PRESSURE TERMS)

    DO J=1,3 
    DO I=1,np(1)-1

       IP=I+1
       IM=I-1
    IF(I.EQ.1)IM=np(1)-1

    DU_DZI=(VAR(IP,J,1)-VAR(IM,J,1))/(2*ZI)
    DV_DZI=(VAR(IP,J,2)-VAR(IM,J,2))/(2*ZI)
   
    IF(J.EQ.1) THEN
    DU_DETA=(VAR(I,J+1,1)-VAR(I,J,1))/ETA
    DV_DETA=(VAR(I,J+1,2)-VAR(I,J,2))/ETA
    ENDIF

    IF(J.GT.1) THEN
    DU_DETA=(VAR(I,J+1,1)-VAR(I,J-1,1))/(2*ETA)
    DV_DETA=(VAR(I,J+1,2)-VAR(I,J-1,2))/(2*ETA)   
    ENDIF
          
    DU_DX=ZI_X(I,J)*DU_DZI+ETA_X(I,J)*DU_DETA
    DU_DY=ZI_Y(I,J)*DU_DZI+ETA_Y(I,J)*DU_DETA

    DV_DX=ZI_X(I,J)*DV_DZI+ETA_X(I,J)*DV_DETA
    DV_DY=ZI_Y(I,J)*DV_DZI+ETA_Y(I,J)*DV_DETA
      

    DILAT=DU_DX+DV_DY
    LAMBDA(1)=DU_DX-DILAT/3
    LAMBDA(2)=DV_DY-DILAT/3
    
    SHEAR=DV_DX+DU_DY

    F(I,J,1)=ROW(I,J)*VAR(I,J,1)*VAR(I,J,1)
    F(I,J,2)=ROW(I,J)*VAR(I,J,2)*VAR(I,J,1)
    G(I,J,1)=F(I,J,2)
    G(I,J,2)=ROW(I,J)*VAR(I,J,2)*VAR(I,J,2)

    F(I,J,1)=F(I,J,1)-2*VIS(I,J)*LAMBDA(1)/RE
    F(I,J,2)=F(I,J,2)-VIS(I,J)*SHEAR/RE
    G(I,J,1)=G(I,J,1)-VIS(I,J)*SHEAR/RE
    G(I,J,2)=G(I,J,2)-2*VIS(I,J)*LAMBDA(2)/RE
    
    END DO ! I LOOP ENDS
!-------------------------PERIODICITY-----------------------------------    
    F(np(1),J,1:2)=F(1,J,1:2)
    G(np(1),J,1:2)=G(1,J,1:2)
    END DO ! J LOOP ENDS

!-----------------------------------------------------------------------
    J=1
    DO I=1,np(1)-1
    
       IP=I+1
       IM=I-1
    IF(I.EQ.1)IM=np(1)-1

    SUM1=0.0d0
    JACOB=ZI_X(I,J)*ETA_Y(I,J)-ZI_Y(I,J)*ETA_X(I,J)

    G_ETA(1)=-ZI_Y(I,J)/JACOB
    G_ETA(2)=+ZI_X(I,J)/JACOB
    DN=ETA
    
    SP(I,1)=0.0d0
    SP(I,2)=gw*(EPS1/(1+EPS1))/(FR*FR)
    IF(I.EQ.1) THEN
    SP(np(1),1)=SP(1,1)
    SP(np(1),2)=SP(1,2)
    ENDIF

    DO L=1,2
    DF_ZI=(F(IP,J,L)-F(IM,J,L))/(2*ZI)
    DF_ETA=(-3*F(I,J,L)+4*F(I,J+1,L)-F(I,J+2,L))/(2*ETA)
    DG_ZI=(G(IP,J,L)-G(IM,J,L))/(2*ZI)
    DG_ETA=(-3*G(I,J,L)+4*G(I,J+1,L)-G(I,J+2,L))/(2*ETA)

    DF_DX=ZI_X(I,J)*DF_ZI+ETA_X(I,J)*DF_ETA
    DG_DY=ZI_Y(I,J)*DG_ZI+ETA_Y(I,J)*DG_ETA

    SUM1=SUM1+(-DF_DX-DG_DY+SP(I,L))*G_ETA(L)
    END DO
    
    CORR=gw*DN*(1.d0/PK)*G_ETA(2)/((1+EPS1)*FR*FR)

    P(I,J)=(P(I,J+1)-DN*SUM1)/(1-CORR)
    ROW(I,J)=((P(I,J)/PK)+1)/VAR(I,J,3)

    END DO  ! I LOOP ENDS
!-------------------------PERIODICITY-----------------------------------
    P(np(1),1)=P(1,1)
    ROW(np(1),1)=ROW(1,1)
!-----------------------------------------------------------------------




    ELSE

   J=1
    DO I=1,np(1)-1
     ROW(I,J)=(4*ROW(I,J+1)-ROW(I,J+2))/3.0
!    IF(I.EQ.ITE)ROW(i,j)=2*ROW(i,j+1)-ROW(i,j+2) 

     UTP=VAR(I,J+1,1)*TS_X(I)+VAR(I,J+1,2)*TS_Y(I)
     UTPP=VAR(I,J+2,1)*TS_X(I)+VAR(I,J+2,2)*TS_Y(I)
           
     UTP=ROW(I,J+1)*UTP
     UTPP=ROW(I,J+2)*UTPP

     TERM=((4*UTP-UTPP)/3.0)/(TS_X(I)*NS_Y(I)-NS_X(I)*TS_Y(I))

     VAR(I,J,1)=NS_Y(I)*TERM/ROW(I,J)
     VAR(I,J,2)=-NS_X(I)*TERM/ROW(I,J)

!          IF(I.EQ.ITE)THEN
!           U(I,J)=0.0
!           V(I,J)=0.0
!          ENDIF

!    T(I,J)=(4*T(I,J+1)-T(I,J+2))/3.0!-T(I,3)

!     E(I,J)=T(I,J)+EK*(U(I,J)**2+V(I,J)**2)

     TERM=(4*U_VEC(I,J+1,4)-U_VEC(I,J+2,4))/3.0
     E(I,J)=TERM/ROW(I,J)
     SP_EN=E(I,J)-EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)
     VAR(I,J,3)=1.0d0+B1*(SP_EN-1)+B2*(SP_EN-1)**2+B3*(SP_EN-1)**3+B4*(SP_EN-1)**4+B5*(SP_EN-1)**5
     P(I,J)=PK*(ROW(I,J)*VAR(I,J,3)-1) 
  
     VIS(I,J)=0.d0
     COND(I,J)=0.d0
       
   END DO
  ENDIF

         VAR(np(1),J,1)=VAR(1,J,1)
         VAR(np(1),J,2)=VAR(1,J,2)
         VAR(np(1),J,3)=VAR(1,J,3)     
         ROW(np(1),J)=ROW(1,J)
         E(np(1),J)=E(1,J)
         P(np(1),J)=P(1,J)
         VIS(np(1),J)=VIS(1,J)
         COND(np(1),J)=COND(1,J)
!-----------------------------------------------------------------------
!       SETTING BOUNDARY CONDITIONS AT ARTIFICIAL SURFACE
!-----------------------------------------------------------------------
    J=np(2)
 DO I=1,np(1)-1

    A(I)=(uinf*NA_X(I)+vinf*NA_Y(I))

 IF (A(I).GE.0.d0) THEN
    VAR(I,J,1)=uinf
    VAR(I,J,2)=vinf
    ROW(I,J)=1.d0 !+A0*SIN((2.D0*PI/LAMBDA_Z)*Z(K)) ! PERTURBATION 
    VAR(I,J,3)=1.d0
    P(I,J)=(ROW(I,J)*VAR(I,J,3)-1)*PK
    E(I,J)=VAR(I,J,3)+APH1*(VAR(I,J,3)-1)**2+APH2*(VAR(I,J,3)-1)**3 &
             +APH3*(VAR(I,J,3)-1)**4+EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)  
 ELSE
    VAR(I,J,1)=2*VAR(I,J-1,1)-VAR(I,J-2,1)
    VAR(I,J,2)=2*VAR(I,J-1,2)-VAR(I,J-2,2)
    ROW(I,J)  =2*ROW(I,J-1)-ROW(I,J-2)
    VAR(I,J,3)=2.0D0*VAR(I,J-1,3)-VAR(I,J-2,3)
    E(I,J)=VAR(I,J,3)+APH1*(VAR(I,J,3)-1)**2+APH2*(VAR(I,J,3)-1)**3 &
               +APH3*(VAR(I,J,3)-1)**4+EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)
    P(I,J)=(ROW(I,J)*VAR(I,J,3)-1)*PK
 END IF

 IF(MODEL=='VISCOUS')THEN
    VIS(I,J)=(VAR(I,J,3)**1.5)*(Tr+(S_/T0))/(VAR(I,J,3)*Tr+(S_/T0))
    COND(I,J)=A1+A2*VAR(I,J,3)+A3*VAR(I,J,3)**2
 ELSE
     VIS(I,J)=0.d0
    COND(I,J)=0.d0  
 ENDIF

 END DO

!-------------------------PERIODICITY-----------------------------------
    J=np(2)    
    VAR(np(1),J,1:3)=VAR(1,J,1:3)
    ROW(np(1),J)=ROW(1,J)
    E(np(1),J)=E(1,J)
    P(np(1),J)=P(1,J)
!IF (MODEL.EQ.'VISCOUS') THEN
    VIS(np(1),J)=VIS(1,J)
    COND(np(1),J)=COND(1,J)
!ENDIF  

!-----------------------------------------------------------------------
!     TIME LOOP BEGINS
!-----------------------------------------------------------------------
101   CONTINUE
    
    etol=1.0E-06
    SI(1)=0.1D0  
    SI(2)=0.1D0

DO WHILE (nstep < maxstep)

	IF (mod(nstep,100)==0) CALL cpu_time(tStart)
	nstep = nstep+1
	time  = time + dt

!IF(id==0) WRITE(*,*)nstep

    DO J=1,np(2)!bs,es
    DO I=1,np(1)-1
    U_VEC_OLD(I,J,1)=ROW(I,J)
    U_VEC_OLD(I,J,2)=ROW(I,J)*VAR(I,J,1)
    U_VEC_OLD(I,J,3)=ROW(I,J)*VAR(I,J,2)
    U_VEC_OLD(I,J,4)=ROW(I,J)*E(I,J)
    
    V_CONT(I,J,1)=VAR(I,J,1)*ZI_X(I,J)+VAR(I,J,2)*ZI_Y(I,J)
    V_CONT(I,J,2)=VAR(I,J,1)*ETA_X(I,J)+VAR(I,J,2)*ETA_Y(I,J)

    DO L=1,4
    U_VEC(I,J,L)=U_VEC_OLD(I,J,L)
    END DO
    
   IF(MODEL=='VISCOUS') THEN
    VIS(I,J)=(VAR(I,J,3)**1.5)*(Tr+(S_/T0))/(VAR(I,J,3)*Tr+(S_/T0))  
    COND(I,J)=A1+A2*VAR(I,J,3)+A3*VAR(I,J,3)**2
   ELSE
    VIS(I,J)=0.d0
    COND(I,J)=0.d0  
   ENDIF

!   ENDIF

    END DO !I LOOP ENDS
!-------------------------PERIODICITY-----------------------------------   
    DO L=1,4
    U_VEC(np(1),J,L)=U_VEC(1,J,L)
    END DO
    V_CONT(np(1),J,1:2)=V_CONT(1,J,1:2)   
!   IF (MODEL.EQ.'VISCOUS') THEN
    VIS(np(1),J)=VIS(1,J)
    COND(np(1),J)=COND(1,J)
!   ENDIF  
    END DO  !J LOOP ENDS 
!-----------------------------------------------------------------------  


    DO W_L=1,2 ! STARTING PREDICTOR-CORRECTOR CYCLE

    DO J=1,np(2)
    DO I=1,np(1)
     q_vec(I,J,1)=ROW(I,J)
     q_vec(I,J,2)=ROW(I,J)*VAR(I,J,1)
     q_vec(I,J,3)=ROW(I,J)*VAR(I,J,2)
     q_vec(I,J,4)=ROW(I,J)*E(I,J)
    IF(ENTH_TOT=='yes')then
     q_vec(I,J,4)=q_vec(I,J,4)+(GAMA-1)*(P(I,J)*(1./PK)+1)
    endif
     
    END DO
    END DO
!-----------------------------------------------------------------------
!     COMPUTING CONVECTIVE FLUXES AT I+1/2 AND J+1/2
!-----------------------------------------------------------------------

!***************** INTERCEL VELOCITY CALCULATION U_INTCEL **************************
     
    DO J=2,np(2)-1!bn,en
    DO I=1,np(1)-1

    MEAN=0.5*(V_CONT(I,J,1)+V_CONT(I+1,J,1))

    IP=I+1
    IPP=I+2
    IM=I-1
    IF(I.EQ.1)IM=np(1)-1
    IF(I.EQ.np(1)-1)IPP=2
     
    RPLUSU=V_CONT(IPP,J,1)-V_CONT(IP,J,1)
    RMINUSU=V_CONT(IP,J,1)-V_CONT(I,J,1)
    DROP=ABS(RPLUSU)+ABS(RMINUSU)

    IF(DROP.LE.etol) THEN
    WPPU=0.2d0
    ELSE
    WPPU=ABS(RPLUSU-RMINUSU)/DROP
    END IF ! K1

    RPLUSU=V_CONT(IP,J,1)-V_CONT(I,J,1)
    RMINUSU=V_CONT(I,J,1)-V_CONT(IM,J,1)
    DROP=ABS(RPLUSU)+ABS(RMINUSU)

    IF(DROP.LE.etol) THEN
    WPU=0.2d0
    ELSE
    WPU=ABS(RPLUSU-RMINUSU)/DROP
    END IF !K2
    WU=WPU
    IF(WPPU.ge.WPU)WU=WPPU

    IF(MEAN.GT.0.d0) ff=1.d0
    IF(MEAN.LT.0.d0) ff=-1.d0
    IF(MEAN.EQ.0.d0) ff=0.d0

    z1=ABS(V_CONT(I,J,1))
    z2=ABS(V_CONT(IP,J,1))
    IF (z1.GT.z2) THEN
    U_SN=z1               
    ELSE
    U_SN=z2                
    END IF
    WF= U_SN/(U_SN+SI(1))
    
    URU1=-V_CONT(IM,J,1)+9*V_CONT(I,J,1)+9*V_CONT(IP,J,1)-V_CONT(IPP,J,1)
    
    URU2=-V_CONT(IM,J,1)+3*V_CONT(I,J,1)-3*V_CONT(IP,J,1)+V_CONT(IPP,J,1)

    URU=(1.d0/16)*(URU1+ff*URU2)

    URC=(1.d0/16)*(-V_CONT(IM,J,1)+9*V_CONT(I,J,1)+9*V_CONT(IP,J,1)-V_CONT(IPP,J,1))
    
    URU=URC+WF*(URU-URC)
    UR1=URU+WU*(URC-URU)
    UR2=URC+WU*(URU-URC)

!*************** RANGE BOUNDEDNESS CRITERIA *********************

    MAXM=V_CONT(I+1,J,1)
    MINM=V_CONT(I,J,1)
    IF(MINM.GT.MAXM) THEN
    MAXM=V_CONT(I,J,1)
    MINM=V_CONT(I+1,J,1)
    ENDIF

    LOGIC1=(UR1.GE.MINM).AND.(UR1.LE.MAXM)
    LOGIC2=(UR2.GE.MINM).AND.(UR2.LE.MAXM)

    IF(LOGIC1.OR.LOGIC2) THEN
       IF(LOGIC1)U_INT=UR1
       IF(LOGIC2)U_INT=UR2
       IF(LOGIC1.AND.LOGIC2)U_INT=0.5*(UR1+UR2)
    ELSE
       U_INT=0.5*(V_CONT(I+1,J,1)+V_CONT(I,J,1))
    ENDIF


!******************* INTERCELL CONVECTIVE FLUX FC ********************
       
    DO L=1,4
 
    RPLUSU=q_vec(IP,J,L)-q_vec(I,J,L)
    RMINUSU=q_vec(I,J,L)-q_vec(IM,J,L)
    DROP=ABS(RPLUSU)+abs(RMINUSU)  

    IF(DROP.LE.etol) THEN	
    WPPU=0.2d0
    ELSE
    WPPU= ABS(RPLUSU-RMINUSU)/DROP
    END IF                                             

    RPLUSU=q_vec(IPP,J,L)-q_vec(IP,J,L)
    RMINUSU=q_vec(IP,J,L)-q_vec(I,J,L)
    DROP=ABS(RPLUSU)+abs(RMINUSU) 
 
    IF(DROP.LE.etol) THEN	
    WPU=0.2d0
    ELSE
    WPU=ABS(RPLUSU-RMINUSU)/DROP
    END IF !f2

    WU=WPPU
    IF(WPU.GT.WPPU) WU=WPU

    IF(MEAN.GT.0.d0) ff=1.d0
    IF(MEAN.LT.0.d0) ff=-1.d0
    IF(MEAN.EQ.0.d0) ff=0.d0

    z1=ABS(V_CONT(I,J,1))
    z2=ABS(V_CONT(IP,J,1))
    IF (z1.GT.z2) THEN
    U_SN=z1               
    ELSE
    U_SN=z2                
    END IF
    WF= U_SN/(U_SN+SI(1))

    URU_L=0.5*(q_vec(I,J,L)+q_vec(I+1,J,L))+ff*0.5*(q_vec(I,J,L)-q_vec(I+1,J,L))

    URU1=-q_vec(IM,J,L)+9*q_vec(I,J,L)+9*q_vec(IP,J,L)-q_vec(IPP,J,L)

    URU2=-q_vec(IM,J,L)+3*q_vec(I,J,L)-3*q_vec(IP,J,L)+q_vec(IPP,J,L)

    URU= (1.d0/16)*(URU1+ff*URU2)

    URC=(1.d0/16)*(-q_vec(IM,J,L)+9*q_vec(I,J,L)+9*q_vec(IP,J,L)-q_vec(IPP,J,L))

    URU=URC+WF*(URU-URC)
    PHIF_R=URU+WU*(URU_L-URU)

!*********************RANGE BOUNDEDNESS*******************************

    MAXM=q_vec(IP,J,L)
    MINM=q_vec(I,J,L)
    IF(MINM.GT.MAXM) THEN
    MAXM=q_vec(I,J,L)
    MINM=q_vec(IP,J,L)
    ENDIF

    LOGIC1=(PHIF_R.GE.MINM).AND.(PHIF_R.LE.MAXM)
  
    IF(LOGIC1) THEN

    FC(I,J,L)=U_INT*PHIF_R
    ELSE
    FC(I,J,L)=U_INT*0.5*(q_vec(IP,J,L)+q_vec(I,J,L))
    ENDIF
    END DO

    END DO !I LOOP ENDS
!-------------------------PERIODICITY-----------------------------------    
    DO L=1,4
    FC(np(1),J,L)=FC(1,J,L) 
    END DO 
!-----------------------------------------------------------------------    
    END DO ! J LOOP ENDS

!****************** INTERCEL VELOCITY CALCULATION V_INTCEL *************
    DO J=1,np(2)-1!bs,en
    DO I=1,np(1)-1

    MEAN=0.5*(V_CONT(I,J,2)+V_CONT(I,J+1,2))

    IF(J.EQ.1.OR.J.EQ.np(2)-1) THEN
    V_INT=(V_CONT(I,J,2)+V_CONT(I,J+1,2))/2.0D0
    ELSE
     
    RPLUSU=V_CONT(I,J+2,2)-V_CONT(I,J+1,2)
    RMINUSU=V_CONT(I,J+1,2)-V_CONT(I,J,2)
    DROP=ABS(RPLUSU)+abs(RMINUSU)

    IF(DROP.LE.etol) THEN  
    WPPU=0.2D0
    ELSE
    WPPU=ABS(RPLUSU-RMINUSU)/DROP
    END IF !K3
    RPLUSU=V_CONT(I,J+1,2)-V_CONT(I,J,2)
    RMINUSU=V_CONT(I,J,2)-V_CONT(I,J-1,2)
    DROP=ABS(RPLUSU)+abs(RMINUSU)

    IF(DROP.LE.etol) THEN 
    WPU=0.2D0
    ELSE
    WPU=ABS(RPLUSU-RMINUSU)/DROP
    END IF !K4
    WU=WPU
    IF(WPPU.ge.WPU)WU=WPPU
            
    IF(MEAN.GT.0) ff=1.0D0
    IF(MEAN.LT.0) ff=-1.0D0
    IF(MEAN.EQ.0) ff=0.0D0

    z1=ABS(V_CONT(I,J,2))
    z2=ABS(V_CONT(I,J+1,2))
    IF (z1.GT.z2) THEN
     U_SN=z1                
    ELSE
     U_SN=z2               
    END IF
    WF=  U_SN/(U_SN+SI(2))
    URU1=-V_CONT(I,J-1,2)+9.0D0*V_CONT(I,J,2) &
                    +9.0D0*V_CONT(I,J+1,2)-V_CONT(I,J+2,2)
    URU2=-V_CONT(I,J-1,2)+3.0D0*V_CONT(I,J,2) &
                    -3.0D0*V_CONT(I,J+1,2)+V_CONT(I,J+2,2)

    URU=(1.0D0/16.0D0)*( URU1+ff* URU2)

    URC=(1.0D0/16.0D0)*(-V_CONT(I,J-1,2)+9.0D0*V_CONT(I,J,2) &
                     +9.0D0*V_CONT(I,J+1,2)-V_CONT(I,J+2,2))
     
    URU= URC+WF*( URU- URC)
    UR1= URU+WU*( URC- URU)
    UR2= URC+WU*( URU- URC)

!*************** RANGE BOUNDEDNESS CRITERIA *********************

    MAXM=V_CONT(I,J+1,2)
    MINM=V_CONT(I,J,2)
    IF(MINM.GT.MAXM) THEN
    MAXM=V_CONT(I,J,2)
    MINM=V_CONT(I,J+1,2)
    ENDIF

    LOGIC1=(UR1.GE.MINM).AND.(UR1.LE.MAXM)
    LOGIC2=(UR2.GE.MINM).AND.(UR2.LE.MAXM)

    IF(LOGIC1.OR.LOGIC2) THEN
       IF(LOGIC1)V_INT=UR1
       IF(LOGIC2)V_INT=UR2
       IF(LOGIC1.AND.LOGIC2)V_INT=(UR1+UR2)/2.0D0
    ELSE
       V_INT=(V_CONT(I,J+1,2)+V_CONT(I,J,2))/2.0D0
    ENDIF

     IF(VEL_ETA=='mean')THEN
     V_INT=(V_CONT(I,J,2)+V_CONT(I,J+1,2))/2.0D0
     ENDIF

    ENDIF

!******************* INTERCELL CONVECTIVE FLUX GC ********************

    IF(J.EQ.1.OR.J.EQ.np(2)-1) THEN

    DO L=1,4
    IF(MEAN.GE.0.d0) THEN
    GC(I,J,L)=V_INT*q_vec(I,J,L)
    ELSE
    GC(I,J,L)=V_INT*q_vec(I,J+1,L)
    ENDIF
    END DO

    ELSE   !  IF(J.EQ.1.OR.J.EQ.np(2)-1) THEN

    DO L=1,4
        
    RPLUSU=q_vec(I,J+1,L)-q_vec(I,J,L)
    RMINUSU=q_vec(I,J,L)-q_vec(I,J-1,L)
    DROP=ABS(RPLUSU)+ABS(RMINUSU)  

    IF(DROP.LE.etol) THEN	
    WPPU=0.2d0
    ELSE
    WPPU= ABS(RPLUSU-RMINUSU)/DROP
    END IF                                        

    RPLUSU=q_vec(I,J+2,L)-q_vec(I,J+1,L)
    RMINUSU=q_vec(I,J+1,L)- q_vec(I,J,L)
    DROP=ABS(RPLUSU)+ABS(RMINUSU)  
 
    IF(DROP.LE.etol) THEN	
    WPU=0.2d0
    ELSE
    WPU= ABS(RPLUSU-RMINUSU)/DROP
    END IF
                                      
    WU=WPPU
    IF(WPU.GT.WPPU)WU=WPU     
         
    IF(MEAN.GT.0.d0) ff=1.d0
    IF(MEAN.LT.0.d0) ff=-1.d0
    IF(MEAN.EQ.0.d0) ff=0.d0

    z1=ABS(V_CONT(I,J,2))
    z2=ABS(V_CONT(I,J+1,2))
    IF (z1.GT.z2) THEN
     U_SN=z1                
    ELSE
     U_SN=z2               
    END IF
    WF=U_SN/(U_SN+SI(2))

    URU_L=0.5*(q_vec(I,J,L)+q_vec(I,J+1,L))+ff*0.5*(q_vec(I,J,L)-q_vec(I,J+1,L))
                      
    URU1=-q_vec(I,J-1,L)+9*q_vec(I,J,L)+9*q_vec(I,J+1,L)-q_vec(I,J+2,L)
    
    URU2=-q_vec(I,J-1,L)+3*q_vec(I,J,L)-3*q_vec(I,J+1,L)+q_vec(I,J+2,L)

    URU=(1.d0/16)*(URU1+ff*URU2)

    URC=(1.d0/16)*(-q_vec(I,J-1,L)+9*q_vec(I,J,L)+9*q_vec(I,J+1,L)-q_vec(I,J+2,L))

    URU=URC+WF*(URU-URC)

    PHIF_R=URU+WU*(URU_L-URU) 

!*********************RANGE BOUNDEDNESS*******************************
          
    MAXM=q_vec(I,J+1,L)
    MINM=q_vec(I,J,L)
    IF(MINM.GT.MAXM) THEN
    MAXM=q_vec(I,J,L)
    MINM=q_vec(I,J+1,L)
    ENDIF

    LOGIC1=(PHIF_R.GE.MINM).AND.(PHIF_R.LE.MAXM)
  
    IF(LOGIC1) THEN
    GC(I,J,L)=V_INT*PHIF_R
    ELSE
    GC(I,J,L)=V_INT*0.5*(q_vec(I,J+1,L)+q_vec(I,J,L))
    ENDIF

    END DO

    ENDIF

    END DO
    END DO
 
!     ----------------------------------------------------------------------------
!      COMPUTING DIFFUSIVE FLUXES AND SOURCE VECTOR
!     ----------------------------------------------------------------------------
!========================FD============================================
    DO J=1,np(2)-1!bs,en     !bs,es
    DO I=1,np(1)-1
    
    IP=I+1
    IM=I-1
    IF(I.EQ.1)IM=np(1)-1
!-----------------------------------------------------------------------   
    ZI_XA=0.5*(ZI_X(I,J)+ZI_X(IP,J))
    ZI_YA=0.5*(ZI_Y(I,J)+ZI_Y(IP,J))       
    ETA_XA=0.5*(ETA_X(I,J)+ETA_X(IP,J))
    ETA_YA=0.5*(ETA_Y(I,J)+ETA_Y(IP,J))  
!-----------------------------------------------------------------------
    DU_DZI=(VAR(IP,J,1)-VAR(I,J,1))/ZI
    DV_DZI=(VAR(IP,J,2)-VAR(I,J,2))/ZI
    DT_DZI=(VAR(IP,J,3)-VAR(I,J,3))/ZI
!-----------------------------------------------------------------------    
    IF(J.EQ.1) THEN
    FUN1=(VAR(I,J+1,1)-VAR(I,J,1))/ETA
    FUN2=(VAR(IP,J+1,1)-VAR(IP,J,1))/ETA
    DU_DETA=0.5*(FUN1+FUN2)
    
    FUN1=(VAR(I,J+1,2)-VAR(I,J,2))/ETA
    FUN2=(VAR(IP,J+1,2)-VAR(IP,J,2))/ETA
    DV_DETA=0.5*(FUN1+FUN2)
    
    FUN1=(VAR(I,J+1,3)-VAR(I,J,3))/ETA
    FUN2=(VAR(IP,J+1,3)-VAR(IP,J,3))/ETA
    DT_DETA=0.5*(FUN1+FUN2)   
       ELSE      
    FUN1=(VAR(I,J+1,1)-VAR(I,J-1,1))/(2*ETA)
    FUN2=(VAR(IP,J+1,1)-VAR(IP,J-1,1))/(2*ETA)
    DU_DETA=0.5*(FUN1+FUN2)
    
    FUN1=(VAR(I,J+1,2)-VAR(I,J-1,2))/(2*ETA)
    FUN2=(VAR(IP,J+1,2)-VAR(IP,J-1,2))/(2*ETA)
    DV_DETA=0.5*(FUN1+FUN2)
    
    FUN1=(VAR(I,J+1,3)-VAR(I,J-1,3))/(2*ETA)
    FUN2=(VAR(IP,J+1,3)-VAR(IP,J-1,3))/(2*ETA)
    DT_DETA=0.5*(FUN1+FUN2)
    ENDIF 
!-----------------------------------------------------------------------         
    DU_DX=DU_DZI*ZI_XA+DU_DETA*ETA_XA
    DV_DX=DV_DZI*ZI_XA+DV_DETA*ETA_XA   
    DT_DX=DT_DZI*ZI_XA+DT_DETA*ETA_XA
          
    DU_DY=DU_DZI*ZI_YA+DU_DETA*ETA_YA  
    DV_DY=DV_DZI*ZI_YA+DV_DETA*ETA_YA
    DT_DY=DT_DZI*ZI_YA+DT_DETA*ETA_YA
!-----------------------------------------------------------------------   
    DILAT=DU_DX+DV_DY
    LAMBDA(1)=DU_DX-DILAT/3
    LAMBDA(2)=DV_DY-DILAT/3
    
    SHEAR=DV_DX+DU_DY
 
    VIS_INT=0.5*(VIS(I,J)+VIS(IP,J))    
    COND_INT=0.5*(COND(I,J)+COND(IP,J)) 
    U_INT=0.5*(VAR(I,J,1)+VAR(IP,J,1))  
    V_INT=0.5*(VAR(I,J,2)+VAR(IP,J,2))   
    P_INT=0.5*(P(I,J)+P(IP,J))         
    
    PHI(1)=-VIS_INT*BK*(2.d0*U_INT*LAMBDA(1)+V_INT*SHEAR)
    PHI(2)=-VIS_INT*BK*(2.d0*V_INT*LAMBDA(2)+U_INT*SHEAR)
!-----------------------------------------------------------------------    
    FD(I,J,1)=0.d0
    FD(I,J,2)=ZI_XA*P_INT-VIS_INT*(2*ZI_XA*LAMBDA(1)+ZI_YA*SHEAR)/RE
    FD(I,J,3)=ZI_YA*P_INT-VIS_INT*(2*ZI_YA*LAMBDA(2)+ZI_XA*SHEAR)/RE
    IF(ENTH_TOT=='yes')THEN
    FD(I,J,4)=-GAMA*COND_INT &
               *(ZI_XA*DT_DX+ZI_YA*DT_DY)/(RE*PR)+ZI_XA*PHI(1)+ZI_YA*PHI(2)
    ELSE
    FD(I,J,4)=(GAMA-1)*(ZI_XA*U_INT+ZI_YA*V_INT)*(P_INT*(1.d0/PK)+1)-GAMA*COND_INT &
               *(ZI_XA*DT_DX+ZI_YA*DT_DY)/(RE*PR)+ZI_XA*PHI(1)+ZI_YA*PHI(2)
    ENDIF       
!-----------------------------------------------------------------------       
    END DO !I LOOP ENDS
    
!----------------PERIODICITY--------------------
    DO L=1,4
    FD(np(1),J,L)=FD(1,J,L)   
    END DO  
!-----------------------------------------------    
    END DO !J LOOP ENDS  

!========================GD============================================
    DO J=1,np(2)-1!bs,en     !bs,es
    DO I=1,np(1)-1
    
    IP=I+1
    IM=I-1
    IF(I.EQ.1)IM=np(1)-1
!-----------------------------------------------------------------------   
    ZI_XA=0.5*(ZI_X(I,J)+ZI_X(I,J+1))
    ZI_YA=0.5*(ZI_Y(I,J)+ZI_Y(I,J+1))       
    ETA_XA=0.5*(ETA_X(I,J)+ETA_X(I,J+1))
    ETA_YA=0.5*(ETA_Y(I,J)+ETA_Y(I,J+1))   
!-----------------------------------------------------------------------   
    FUN1=(VAR(IP,J,1)-VAR(IM,J,1))/(2*ZI)
    FUN2=(VAR(IP,J+1,1)-VAR(IM,J+1,1))/(2*ZI)
    DU_DZI=0.5*(FUN1+FUN2)
    
    FUN1=(VAR(IP,J,2)-VAR(IM,J,2))/(2*ZI)
    FUN2=(VAR(IP,J+1,2)-VAR(IM,J+1,2))/(2*ZI)
    DV_DZI=0.5*(FUN1+FUN2)

    FUN1=(VAR(IP,J,3)-VAR(IM,J,3))/(2*ZI)
    FUN2=(VAR(IP,J+1,3)-VAR(IM,J+1,3))/(2*ZI)
    DT_DZI=0.5*(FUN1+FUN2)        
!-----------------------------------------------------------------------   
    DU_DETA=(VAR(I,J+1,1)-VAR(I,J,1))/ETA   
    DV_DETA=(VAR(I,J+1,2)-VAR(I,J,2))/ETA     
    DT_DETA=(VAR(I,J+1,3)-VAR(I,J,3))/ETA 
!-----------------------------------------------------------------------        
    DU_DX=DU_DZI*ZI_XA+DU_DETA*ETA_XA
    DV_DX=DV_DZI*ZI_XA+DV_DETA*ETA_XA  
    DT_DX=DT_DZI*ZI_XA+DT_DETA*ETA_XA
          
    DU_DY=DU_DZI*ZI_YA+DU_DETA*ETA_YA  
    DV_DY=DV_DZI*ZI_YA+DV_DETA*ETA_YA
    DT_DY=DT_DZI*ZI_YA+DT_DETA*ETA_YA
!-----------------------------------------------------------------------   
    DILAT=DU_DX+DV_DY
    LAMBDA(1)=DU_DX-DILAT/3
    LAMBDA(2)=DV_DY-DILAT/3
    
    SHEAR=DV_DX+DU_DY
    
    VIS_INT=0.5*(VIS(I,J)+VIS(I,J+1)) 
    COND_INT=0.5*(COND(I,J)+COND(I,J+1)) 
    U_INT=0.5*(VAR(I,J,1)+VAR(I,J+1,1))  
    V_INT=0.5*(VAR(I,J,2)+VAR(I,J+1,2))   
    P_INT=0.5*(P(I,J)+P(I,J+1))         
    
    PHI(1)=-VIS_INT*BK*(2.d0*U_INT*LAMBDA(1)+V_INT*SHEAR)
    PHI(2)=-VIS_INT*BK*(2.d0*V_INT*LAMBDA(2)+U_INT*SHEAR)   
!-----------------------------------------------------------------------
    GD(I,J,1)=0.d0
    GD(I,J,2)=ETA_XA*P_INT-VIS_INT*(2.d0*ETA_XA*LAMBDA(1)+ETA_YA*SHEAR)/RE
    GD(I,J,3)=ETA_YA*P_INT-VIS_INT*(2.d0*ETA_YA*LAMBDA(2)+ETA_XA*SHEAR)/RE
    IF(ENTH_TOT=='yes')THEN
    GD(I,J,4)=-GAMA*COND_INT &
               *(ETA_XA*DT_DX+ETA_YA*DT_DY)/(RE*PR)+ETA_XA*PHI(1)+ETA_YA*PHI(2)    
    ELSE
    GD(I,J,4)=(GAMA-1)*(ETA_XA*U_INT+ETA_YA*V_INT)*(P_INT*(1.d0/PK)+1)-GAMA*COND_INT &
               *(ETA_XA*DT_DX+ETA_YA*DT_DY)/(RE*PR)+ETA_XA*PHI(1)+ETA_YA*PHI(2)
    ENDIF
!-----------------------------------------------------------------------    
    END DO
    END DO    

!========================SOURCE=========================================
    DO J=2,np(2)-1!bn,en
    DO I=1,np(1)-1
    
    IP=I+1
    IM=I-1
    IF(I.EQ.1)IM=np(1)-1 
!-----------------------------------------------------------------------  
    DZI_X_ZI=(ZI_X(IP,J)-ZI_X(IM,J))/(2*ZI)
    DZI_Y_ZI=(ZI_Y(IP,J)-ZI_Y(IM,J))/(2*ZI)

    DETA_X_ETA=(ETA_X(I,J+1)-ETA_X(I,J-1))/(2*ETA)
    DETA_Y_ETA=(ETA_Y(I,J+1)-ETA_Y(I,J-1))/(2*ETA)

    HH=DZI_X_ZI+DETA_X_ETA
    QQ=DZI_Y_ZI+DETA_Y_ETA
    
    JACOB=ZI_X(I,J)*ETA_Y(I,J)-ZI_Y(I,J)*ETA_X(I,J)
!-----------------------------------------------------------------------   
    S(I,J,1)=0.d0
    S(I,J,2)=0.d0
    S(I,J,3)=-gw*(ROW(I,J)-1)/(FR*FR)
    S(I,J,4)=gw*(1-ROW(I,J))*VAR(I,J,2)*BK*RE/(FR*FR)
!-----------------------------------------------------------------------    
    DO L=1,4
    FCN=0.5*(FC(IM,J,L)+FC(I,J,L))   !FC at node in body-fitted coordinate
    GCN=0.5*(GC(I,J-1,L)+GC(I,J,L))  !GC at node in body-fitted coordinate
    
    FCC=(ETA_Y(I,J)*FCN-ZI_Y(I,J)*GCN)/JACOB  !FC at node in cartesian coordinate
    GCC=(ZI_X(I,J)*GCN-ETA_X(I,J)*FCN)/JACOB  !GC at node in cartesian coordinate   

    FDN=0.5*(FD(IM,J,L)+FD(I,J,L))   !FD at node in body-fitted coordinate
    GDN=0.5*(GD(I,J-1,L)+GD(I,J,L))  !GD at node in body-fitted coordinate
    
    FDC=(ETA_Y(I,J)*FDN-ZI_Y(I,J)*GDN)/JACOB  !FD at node in cartesian coordinate
    GDC=(ZI_X(I,J)*GDN-ETA_X(I,J)*FDN)/JACOB  !GD at node in cartesian coordinate

    S(I,J,L)=S(I,J,L)+(FCC+FDC)*HH+(GCC+GDC)*QQ
    END DO 

    END DO
    END DO
!-----------------------------------------------------------------------

!========SELECTING PREDICTOR OR COLLECTOR==================
!~     AW=0.4d0  !0.2
!~     BW=0.1d0  !0.01, 0.12, 0.05

    IF(W_L==1) THEN
    
    DO J=2,np(2)-1!bn,en
    DO I=1,np(1)-1    

      IP=I+1
      IM=I-1
    IF(I.EQ.1)IM=np(1)-1

!~ !========================PVM+H==========================================
!~ 
!~ !------------------------ZI-----------------------------------------
!~    ML(1) =MACH*V_CONT(I,J,1)/SQRT(VAR(I,J,3))
!~    MLP(1)=MACH*V_CONT(IP,J,1)/SQRT(VAR(IP,J,3))   
!~    MLM(1)=MACH*V_CONT(IM,J,1)/SQRT(VAR(IM,J,3))
!~    
!~ !------------------------ETA-----------------------------------------
!~    ML(2) =MACH*V_CONT(I,J,2)/SQRT(VAR(I,J,3))
!~    MLP(2)=MACH*V_CONT(I,J+1,2)/SQRT(VAR(I,J+1,3))   
!~    MLM(2)=MACH*V_CONT(I,J-1,2)/SQRT(VAR(I,J-1,3))
!~    
!~ !------------------------Z-----------------------------------------
!~    ML(3) =MACH*/SQRT(VAR(I ,J,3))
!~    MLP(3)=MACH*/SQRT(VAR(I,J,3))   
!~    MLM(3)=MACH*/SQRT(VAR(I,J,3))
!~  
!~ !-----------------weight functions-------------------   
!~   DO L=1,3 
!~   
!~    WH =AW*EXP(-(ML(L)/BW)**2)
!~    WHP=AW*EXP(-(MLP(L)/BW)**2)
!~    WHM=AW*EXP(-(MLM(L)/BW)**2)  
!~    WHA=(WHM+WH+WHP)/3  
!~    
!~    IF(ABS(ML(L)).LE. 1.d0) THEN
!~     FP=0.5*(1.d0+SIN(ML(L)*PI/2))
!~     FM=0.5*(1.d0-SIN(ML(L)*PI/2))
!~    ENDIF 
!~    
!~    IF(ABS(ML(L)).GT. 1.d0) THEN
!~     IF(ML(L).GT.0) ff=1.d0
!~     IF(ML(L).LT.0) ff=-1.d0
!~     IF(ML(L).EQ.0) ff=0.d0
!~     
!~     FP=0.5*(1.d0+ff)
!~     FM=0.5*(1.d0-ff)
!~    ENDIF 
!~    
!~    IF(ABS(MLP(L)).LE. 1.d0) FMMP=0.5*(1.d0-SIN(MLP(L)*PI/2))
!~    
!~    IF(ABS(MLP(L)).GT. 1.d0) THEN
!~     IF(MLP(L).GT.0) ff=1.d0
!~     IF(MLP(L).LT.0) ff=-1.d0
!~     IF(MLP(L).EQ.0) ff=0.d0
!~     
!~     FMMP=0.5*(1.d0-ff)
!~    ENDIF
!~    
!~    IF(ABS(MLM(L)).LE. 1.d0) FPMM=0.5*(1.d0+SIN(MLM(L)*PI/2))
!~    
!~    IF(ABS(MLM(L)).GT. 1.d0) THEN   
!~     IF(MLM(L).GT.0) ff=1.d0
!~     IF(MLM(L).LT.0) ff=-1.d0
!~     IF(MLM(L).EQ.0) ff=0.d0
!~    
!~     FPMM=0.5*(1.d0+ff)
!~    ENDIF
!~ !-----------------------------------------------------     
!~    IF (L==1) THEN
!~    PT =ZI_X(I ,J)*P(I ,J)
!~    PTP=ZI_X(IP,J)*P(IP,J)
!~    PTM=ZI_X(IM,J)*P(IM,J)
!~ 
!~    Piph=0.5*(ZI_X(IP,J)+ZI_X(I ,J))*(FP*P(I,J)+FMMP*P(IP,J))
!~    Pimh=0.5*(ZI_X(I ,J)+ZI_X(IM,J))*(FPMM*P(IM,J)+FM*P(I,J))
!~    DP_DZI(2)=WHA*(Piph-Pimh)/ZI+(1.d0-WHA)*(PTP-PT)/ZI  
!~    
!~    PT =ZI_Y(I ,J)*P(I ,J)
!~    PTP=ZI_Y(IP,J)*P(IP,J)
!~    PTM=ZI_Y(IM,J)*P(IM,J)
!~  
!~    Piph=0.5*(ZI_Y(IP,J)+ZI_Y(I ,J))*(FP*P(I,J)+FMMP*P(IP,J))
!~    Pimh=0.5*(ZI_Y(I ,J)+ZI_Y(IM,J))*(FPMM*P(IM,J)+FM*P(I,J))
!~    DP_DZI(3)=WHA*(Piph-Pimh)/ZI+(1.d0-WHA)*(PTP-PT)/ZI 
!~      
!~    DP_DZI(1)=0.d0  
!~    DP_DZI(4)=0.d0
!~    DP_DZI(5)=0.d0
!~    ENDIF
!~  
!~    IF (L==2) THEN  
!~    PT =ETA_X(I,J  )*P(I,J  )
!~    PTP=ETA_X(I,J+1)*P(I,J+1)
!~    PTM=ETA_X(I,J-1)*P(I,J-1)
!~    
!~    Piph=0.5*(ETA_X(I,J+1)+ETA_X(I,J  ))*(FP*P(I,J)+FMMP*P(I,J+1))
!~    Pimh=0.5*(ETA_X(I,J  )+ETA_X(I,J-1))*(FPMM*P(I,J-1)+FM*P(I,J))
!~    DP_DETA(2)=WHA*(Piph-Pimh)/ETA+(1.d0-WHA)*(PTP-PT)/ETA  
!~  
!~    PT =ETA_Y(I,J  )*P(I,J  )
!~    PTP=ETA_Y(I,J+1)*P(I,J+1)
!~    PTM=ETA_Y(I,J-1)*P(I,J-1)
!~    
!~    Piph=0.5*(ETA_Y(I,J+1)+ETA_Y(I,J  ))*(FP*P(I,J)+FMMP*P(I,J+1))
!~    Pimh=0.5*(ETA_Y(I,J  )+ETA_Y(I,J-1))*(FPMM*P(I,J-1)+FM*P(I,J))
!~    DP_DETA(3)=WHA*(Piph-Pimh)/ETA+(1.d0-WHA)*(PTP-PT)/ETA  
!~     
!~    DP_DETA(1)=0.d0
!~    DP_DETA(4)=0.d0
!~    DP_DETA(5)=0.d0 
!~    ENDIF
!~ 
!~    IF (L==3) THEN    
!~    Piph=FP*P(I,J)+FMMP*P(I,J)
!~    Pimh=FPMM*P(I,J)+FM*P(I,J) 
!~    DP_DZ(4)=WHA*(Piph-Pimh)/DZ+(1.d0-WHA)*(P(I,J)-P(I,J))/DZ
!~    DP_DZ(1)=0.d0
!~    DP_DZ(2)=0.d0
!~    DP_DZ(3)=0.d0 
!~    DP_DZ(5)=0.d0   
!~    ENDIF  
!~ 
!~    END DO
 !---------------------------------------------------------------------
    DO L=1,4  
    U_VEC(I,J,L)=U_VEC(I,J,L)- DT &
     *((FC(I,J,L)-FC(IM,J,L))/ZI  + (FD(I,J,L)-FD(IM,J,L))/ZI  &               
      +(GC(I,J,L)-GC(I,J-1,L))/ETA+ (GD(I,J,L)-GD(I,J-1,L))/ETA - S(I,J,L))
    END DO

    END DO  !I LOOP ENDS
!-------------------------PERIODICITY-----------------------------------
    DO L=1,4
    U_VEC(np(1),J,L)=U_VEC(1,J,L)
    END DO
    END DO  !J LOOP ENDS
!-------------------------------------------------------------------

    DO J=2,np(2)-1!bn,en
    DO I=1,np(1)-1
      
    ROW(I,J)=U_VEC(I,J,1)
    VAR(I,J,1)=U_VEC(I,J,2)/ROW(I,J)
    VAR(I,J,2)=U_VEC(I,J,3)/ROW(I,J)
    E(I,J)=U_VEC(I,J,4)/ROW(I,J)
    SP_EN=E(I,J)-EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)
 !  VAR(I,J,3)=-0.049208+1.082337*SP_EN-0.034044*SP_EN**2-0.00022219*SP_EN**3
    VAR(I,J,3)=1.0d0+B1*(SP_EN-1)+B2*(SP_EN-1)**2+B3*(SP_EN-1)**3+B4*(SP_EN-1)**4+B5*(SP_EN-1)**5
    P(I,J)=(ROW(I,J)*VAR(I,J,3)-1.d0)*PK
    V_CONT(I,J,1)=VAR(I,J,1)*ZI_X(I,J)+VAR(I,J,2)*ZI_Y(I,J)        
    V_CONT(I,J,2)=VAR(I,J,1)*ETA_X(I,J)+VAR(I,J,2)*ETA_Y(I,J)

    IF(MODEL=='VISCOUS') THEN
    VIS(I,J)=(VAR(I,J,3)**1.5)*(Tr+(S_/T0))/(VAR(I,J,3)*Tr+(S_/T0))  
    COND(I,J)=A1+A2*VAR(I,J,3)+A3*VAR(I,J,3)**2
    ELSE
    VIS(I,J)=0.d0  
    COND(I,J)=0.d0
    ENDIF
  
    END DO  !I LOOP ENDS    
    
!-------------------------PERIODICITY-----------------------------------
    VAR(np(1),J,1:3)=VAR(1,J,1:3)
    ROW(np(1),J)=ROW(1,J)
    E(np(1),J)=E(1,J)
    P(np(1),J)=P(1,J)
    V_CONT(np(1),J,1:2)=V_CONT(1,J,1:2)
!   IF (MODEL.EQ.'VISCOUS') THEN
    VIS(np(1),J)=VIS(1,J)
    COND(np(1),J)=COND(1,J)
!   ENDIF
!-----------------------------------------------------------------------
    END DO  !J LOOP ENDS
!-------------------------------------------------------------------

   ELSE    ! Selecting predictor or corrector
    
    DO J=2,np(2)-1!bn,en
    DO I=1,np(1)-1
    
      IP=I+1
      IM=I-1
    IF(I.EQ.1)IM=np(1)-1
    
!~ !========================PVM+H==========================================
!~ 
!~ !------------------------ZI-----------------------------------------
!~    ML(1) =MACH*V_CONT(I,J,1)/SQRT(VAR(I,J,3))
!~    MLP(1)=MACH*V_CONT(IP,J,1)/SQRT(VAR(IP,J,3))   
!~    MLM(1)=MACH*V_CONT(IM,J,1)/SQRT(VAR(IM,J,3))
!~    
!~ !------------------------DETA-----------------------------------------
!~    ML(2) =MACH*V_CONT(I,J,2)/SQRT(VAR(I,J,3))
!~    MLP(2)=MACH*V_CONT(I,J+1,2)/SQRT(VAR(I,J+1,3))   
!~    MLM(2)=MACH*V_CONT(I,J-1,2)/SQRT(VAR(I,J-1,3))
!~    
!~ !------------------------Z-----------------------------------------
!~    ML(3) =MACH*/SQRT(VAR(I ,J,3))
!~    MLP(3)=MACH*/SQRT(VAR(I,J,3))   
!~    MLM(3)=MACH*/SQRT(VAR(I,J,3))
!~  
!~ !-----------------weight functions-------------------   
!~   DO L=1,3 
!~   
!~    WH =AW*EXP(-(ML(L)/BW)**2)
!~    WHP=AW*EXP(-(MLP(L)/BW)**2)
!~    WHM=AW*EXP(-(MLM(L)/BW)**2)  
!~    WHA=(WHM+WH+WHP)/3 
!~     
!~    IF(ABS(ML(L)).LE. 1.d0) THEN
!~     FP=0.5*(1.d0+SIN(ML(L)*PI/2))
!~     FM=0.5*(1.d0-SIN(ML(L)*PI/2))
!~    ENDIF 
!~    
!~    IF(ABS(ML(L)).GT. 1.d0) THEN
!~     IF(ML(L).GT.0) ff=1.d0
!~     IF(ML(L).LT.0) ff=-1.d0
!~     IF(ML(L).EQ.0) ff=0.d0
!~     
!~     FP=0.5*(1.d0+ff)
!~     FM=0.5*(1.d0-ff)
!~    ENDIF 
!~    
!~    IF(ABS(MLP(L)).LE. 1.d0) FMMP=0.5*(1.d0-SIN(MLP(L)*PI/2))
!~    
!~    IF(ABS(MLP(L)).GT. 1.d0) THEN
!~     IF(MLP(L).GT.0) ff=1.d0
!~     IF(MLP(L).LT.0) ff=-1.d0
!~     IF(MLP(L).EQ.0) ff=0.d0
!~     
!~     FMMP=0.5*(1.d0-ff)
!~    ENDIF
!~    
!~    IF(ABS(MLM(L)).LE. 1.d0) FPMM=0.5*(1.d0+SIN(MLM(L)*PI/2))
!~    
!~    IF(ABS(MLM(L)).GT. 1.d0) THEN   
!~     IF(MLM(L).GT.0) ff=1.d0
!~     IF(MLM(L).LT.0) ff=-1.d0
!~     IF(MLM(L).EQ.0) ff=0.d0
!~    
!~     FPMM=0.5*(1.d0+ff)
!~    ENDIF
!~     
!~ !-----------------------------------------------------    
!~    IF (L==1) THEN
!~    PT =ZI_X(I ,J)*P(I ,J)
!~    PTP=ZI_X(IP,J)*P(IP,J)
!~    PTM=ZI_X(IM,J)*P(IM,J)
!~ 
!~    Piph=0.5*(ZI_X(IP,J)+ZI_X(I ,J))*(FP*P(I,J)+FMMP*P(IP,J))
!~    Pimh=0.5*(ZI_X(I ,J)+ZI_X(IM,J))*(FPMM*P(IM,J)+FM*P(I,J))
!~    DP_DZI(2)=WHA*(Piph-Pimh)/ZI+(1.d0-WHA)*(PT-PTM)/ZI  
!~    
!~    PT =ZI_Y(I ,J)*P(I ,J)
!~    PTP=ZI_Y(IP,J)*P(IP,J)
!~    PTM=ZI_Y(IM,J)*P(IM,J)
!~  
!~    Piph=0.5*(ZI_Y(IP,J)+ZI_Y(I ,J))*(FP*P(I,J)+FMMP*P(IP,J))
!~    Pimh=0.5*(ZI_Y(I ,J)+ZI_Y(IM,J))*(FPMM*P(IM,J)+FM*P(I,J))
!~    DP_DZI(3)=WHA*(Piph-Pimh)/ZI+(1.d0-WHA)*(PT-PTM)/ZI 
!~      
!~    DP_DZI(1)=0.d0  
!~    DP_DZI(4)=0.d0
!~    DP_DZI(5)=0.d0
!~    ENDIF
!~  
!~    IF (L==2) THEN  
!~    PT =ETA_X(I,J  )*P(I,J  )
!~    PTP=ETA_X(I,J+1)*P(I,J+1)
!~    PTM=ETA_X(I,J-1)*P(I,J-1)
!~    
!~    Piph=0.5*(ETA_X(I,J+1)+ETA_X(I,J  ))*(FP*P(I,J)+FMMP*P(I,J+1))
!~    Pimh=0.5*(ETA_X(I,J  )+ETA_X(I,J-1))*(FPMM*P(I,J-1)+FM*P(I,J))
!~    DP_DETA(2)=WHA*(Piph-Pimh)/ETA+(1.d0-WHA)*(PT-PTM)/ETA  
!~  
!~    PT =ETA_Y(I,J  )*P(I,J  )
!~    PTP=ETA_Y(I,J+1)*P(I,J+1)
!~    PTM=ETA_Y(I,J-1)*P(I,J-1)
!~    
!~    Piph=0.5*(ETA_Y(I,J+1)+ETA_Y(I,J  ))*(FP*P(I,J)+FMMP*P(I,J+1))
!~    Pimh=0.5*(ETA_Y(I,J  )+ETA_Y(I,J-1))*(FPMM*P(I,J-1)+FM*P(I,J))
!~    DP_DETA(3)=WHA*(Piph-Pimh)/ETA+(1.d0-WHA)*(PT-PTM)/ETA  
!~     
!~    DP_DETA(1)=0.d0
!~    DP_DETA(4)=0.d0
!~    DP_DETA(5)=0.d0 
!~    ENDIF
!~ 
!~    IF (L==3) THEN    
!~    Piph=FP*P(I,J)+FMMP*P(I,J)
!~    Pimh=FPMM*P(I,J)+FM*P(I,J) 
!~    DP_DZ(4)=WHA*(Piph-Pimh)/DZ+(1.d0-WHA)*(P(I,J)-P(I,J))/DZ
!~    DP_DZ(1)=0.d0
!~    DP_DZ(2)=0.d0
!~    DP_DZ(3)=0.d0 
!~    DP_DZ(5)=0.d0   
!~    ENDIF  
!~ 
!~    END DO   
  
 !---------------------------------------------------------------------
    DO L=1,4
    U_VEC(I,J,L)=(U_VEC(I,J,L)+U_VEC_OLD(I,J,L))/2-0.5*DT &
      *((FC(I,J,L)-FC(IM,J,L))/ZI + (FD(I,J,L)-FD(IM,J,L))/ZI  &
      +(GC(I,J,L)-GC(I,J-1,L))/ETA+ (GD(I,J,L)-GD(I,J-1,L))/ETA - S(I,J,L))
    END DO
    
    END DO ! I LOOP
!-------------------------PERIODICITY-----------------------------------
    DO L=1,4
    U_VEC(np(1),J,L)=U_VEC(1,J,L)
    END DO
!-----------------------------------------------------------------------      
    END DO ! J LOOP
!----------------------------------------------------------------------- 

    END IF   ! PREDICTOR-CORRECTOR
            
    END DO   ! predictor-corrector loop ends
!-----------------------------------------------------------------------

    DO J=2,np(2)-2!bn,en
    DO I=1,np(1)-1
     
    ROW(I,J)=U_VEC(I,J,1)
    VAR(I,J,1)=U_VEC(I,J,2)/ROW(I,J)      
    VAR(I,J,2)=U_VEC(I,J,3)/ROW(I,J)
    E(I,J)=U_VEC(I,J,4)/ROW(I,J)
    SP_EN=E(I,J)-EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)
 !  VAR(I,J,3)=-0.049208+1.082337*SP_EN-0.034044*SP_EN**2-0.00022219*SP_EN**3
    VAR(I,J,3)=1.0d0+B1*(SP_EN-1)+B2*(SP_EN-1)**2+B3*(SP_EN-1)**3+B4*(SP_EN-1)**4+B5*(SP_EN-1)**5 
    P(I,J)=(ROW(I,J)*VAR(I,J,3)-1.d0)*PK

    V_CONT(I,J,1)=VAR(I,J,1)*ZI_X(I,J)+VAR(I,J,2)*ZI_Y(I,J)
    V_CONT(I,J,2)=VAR(I,J,1)*ETA_X(I,J)+VAR(I,J,2)*ETA_Y(I,J)

   IF(MODEL=='VISCOUS') THEN
    VIS(I,J)=(VAR(I,J,3)**1.5)*(Tr+(S_/T0))/(VAR(I,J,3)*Tr+(S_/T0))
    COND(I,J)=A1+A2*VAR(I,J,3)+A3*VAR(I,J,3)**2
   ELSE
    VIS(I,J)=0.d0  
    COND(I,J)=0.d0

   ENDIF

    END DO  !I LOOP    
    
!-------------------------PERIODICITY-----------------------------------
    VAR(np(1),J,1:3)=VAR(1,J,1:3)
    ROW(np(1),J)=ROW(1,J)
    E(np(1),J)=E(1,J)
    P(np(1),J)=P(1,J)
    V_CONT(np(1),J,1:2)=V_CONT(1,J,1:2)
!   IF (MODEL.EQ.'VISCOUS') THEN
    VIS(np(1),J)=VIS(1,J)
    COND(np(1),J)=COND(1,J)
!   ENDIF
!-----------------------------------------------------------------------
    END DO  !J LOOP


    IF(MODEL=='VISCOUS')THEN
    ELSE
    if(TVD_E=='yes')then
      DO J=2,30
      DO I=ITU4,ITL4
      ip=i+1
      im=i-1
      if(i.eq.1)im=np(1)-1
      
      VAR1=E(IP,J)
      VAR2=E(I,J)
      VAR3=E(IM,J)
      
      p_max=MAX(VAR1,VAR2,VAR3)
      p_min=MIN(VAR1,VAR2,VAR3)
      p_mean(i)=(VAR1+VAR2+VAR3)/3.0
      dev=dabs(E(I,J)-p_mean(i))
      dev_max=dabs(p_max-p_min)

      if(dev_max==0)Then
       TV(I)=0
      else
       TV(I)=dev/dev_max
      endif

      if(TV(I)>=0.1)THEN
       E(I,J)=p_mean(i)
       SP_EN=E(I,J)-EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)
       VAR(I,J,3)=1.0d0+B1*(SP_EN-1)+B2*(SP_EN-1)**2+B3*(SP_EN-1)**3+B4*(SP_EN-1)**4+B5*(SP_EN-1)**5     
       ROW(I,J)=((P(I,J)/PK)+1)/VAR(I,J,3)
!       P(I,J)=PK*(ROW(I,J)*T(I,J)-1)
       U_VEC(I,J,1)=ROW(I,J)
       U_VEC(I,J,2)=ROW(I,J)*VAR(I,J,1)
       U_VEC(I,J,3)=ROW(I,J)*VAR(I,J,2)
       U_VEC(I,J,4)=ROW(I,J)*E(I,J)
      ENDIF

      END DO
       E(np(1),J)=E(1,J)
       ROW(np(1),J)=ROW(1,J)
       VAR(np(1),J,3)=VAR(1,J,3)

       U_VEC(np(1),J,1)=U_VEC(1,J,1)
       U_VEC(np(1),J,2)=U_VEC(1,J,2)
       U_VEC(np(1),J,3)=U_VEC(1,J,3)
       U_VEC(np(1),J,4)=U_VEC(1,J,4)
      END DO

!     -------------------------------------------------------------------
!         correcting pressures for grid scale osccilations
!     -------------------------------------------------------------------
     
      DO J=2,30
      DO I=ITU4,ITL4
      ip=i+1
      im=i-1
      if(i.eq.1)im=np(1)-1
      
      VAR1=P(IP,J)
      VAR2=P(I,J)
      VAR3=P(IM,J)
      
      p_max=MAX(VAR1,VAR2,VAR3)
      p_min=MIN(VAR1,VAR2,VAR3)
      p_mean(i)=(VAR1+VAR2+VAR3)/3.0
      dev=dabs(P(I,J)-p_mean(i))
      dev_max=dabs(p_max-p_min)

      if(dev_max==0)Then
       TV(I)=0
      else
       TV(I)=dev/dev_max
      endif

      if(TV(I)>=0.1)THEN
       P(I,J)=p_mean(i)
!       T(I,J)=E(I,J)-EK*(U(I,J)**2+V(I,J)**2)
       ROW(I,J)=((P(I,J)/PK)+1)/VAR(I,J,3)
!       P(I,J)=PK*(ROW(I,J)*T(I,J)-1)
       U_VEC(I,J,1)=ROW(I,J)
       U_VEC(I,J,2)=ROW(I,J)*VAR(I,J,1)
       U_VEC(I,J,3)=ROW(I,J)*VAR(I,J,2)
       U_VEC(I,J,4)=ROW(I,J)*E(I,J)
      ENDIF

      END DO
       E(np(1),J)=E(1,J)
       ROW(np(1),J)=ROW(1,J)
       VAR(np(1),J,3)=VAR(1,J,3)
       P(np(1),J)=P(1,J)
       U_VEC(np(1),J,1)=U_VEC(1,J,1)
       U_VEC(np(1),J,2)=U_VEC(1,J,2)
       U_VEC(np(1),J,3)=U_VEC(1,J,3)
       U_VEC(np(1),J,4)=U_VEC(1,J,4)
      END DO
   endif  ! TVD
  ENDIF ! VISCOUS/INVISCID

!=======SETTING BOUNDARY CONDITIONS AT SOLID SURFACE====================
 
    IF(MODEL=='VISCOUS')THEN
    J=1
    DO I=1,np(1)-1

    VAR(I,J,1)=0.0d0
    VAR(I,J,2)=0.0d0 
    VAR(I,J,3)=TS                       
    E(I,J)=VAR(I,J,3)+APH1*(VAR(I,J,3)-1)**2+APH2*(VAR(I,J,3)-1)**3 &
             +APH3*(VAR(I,J,3)-1)**4+EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)
    VIS(I,J)=(VAR(I,J,3)**1.5)*(Tr+(S_/T0))/(VAR(I,J,3)*Tr+(S_/T0))
    COND(I,J)=A1+A2*VAR(I,J,3)+A3*VAR(I,J,3)**2   
    END DO  ! I LOOP ENDS
!-------------------------PERIODICITY-----------------------------------
    VAR(np(1),1,1:3)=VAR(1,1,1:3)
    E(np(1),1)=E(1,1)

    VIS(np(1),1)=VIS(1,1)
    COND(np(1),1)=COND(1,1)
               
!================PRESSURE BOUNDARY CONDITION============================

!CALCULATING FLUX VECTORS AND SOURCE VECTOR (WITHOUT PRESSURE TERMS)

    DO J=1,3 
    DO I=1,np(1)-1

       IP=I+1
       IM=I-1
    IF(I.EQ.1)IM=np(1)-1

    DU_DZI=(VAR(IP,J,1)-VAR(IM,J,1))/(2*ZI)
    DV_DZI=(VAR(IP,J,2)-VAR(IM,J,2))/(2*ZI)
   
    IF(J.EQ.1) THEN
    DU_DETA=(VAR(I,J+1,1)-VAR(I,J,1))/ETA
    DV_DETA=(VAR(I,J+1,2)-VAR(I,J,2))/ETA
    ENDIF

    IF(J.GT.1) THEN
    DU_DETA=(VAR(I,J+1,1)-VAR(I,J-1,1))/(2*ETA)
    DV_DETA=(VAR(I,J+1,2)-VAR(I,J-1,2))/(2*ETA)   
    ENDIF
          
    DU_DX=ZI_X(I,J)*DU_DZI+ETA_X(I,J)*DU_DETA
    DU_DY=ZI_Y(I,J)*DU_DZI+ETA_Y(I,J)*DU_DETA

    DV_DX=ZI_X(I,J)*DV_DZI+ETA_X(I,J)*DV_DETA
    DV_DY=ZI_Y(I,J)*DV_DZI+ETA_Y(I,J)*DV_DETA
      

    DILAT=DU_DX+DV_DY
    LAMBDA(1)=DU_DX-DILAT/3
    LAMBDA(2)=DV_DY-DILAT/3
    
    SHEAR=DV_DX+DU_DY

    F(I,J,1)=ROW(I,J)*VAR(I,J,1)*VAR(I,J,1)
    F(I,J,2)=ROW(I,J)*VAR(I,J,2)*VAR(I,J,1)
    G(I,J,1)=F(I,J,2)
    G(I,J,2)=ROW(I,J)*VAR(I,J,2)*VAR(I,J,2)

    F(I,J,1)=F(I,J,1)-2*VIS(I,J)*LAMBDA(1)/RE
    F(I,J,2)=F(I,J,2)-VIS(I,J)*SHEAR/RE
    G(I,J,1)=G(I,J,1)-VIS(I,J)*SHEAR/RE
    G(I,J,2)=G(I,J,2)-2*VIS(I,J)*LAMBDA(2)/RE
    
    END DO ! I LOOP ENDS
!-------------------------PERIODICITY-----------------------------------    
    F(np(1),J,1:2)=F(1,J,1:2)
    G(np(1),J,1:2)=G(1,J,1:2)
    END DO ! J LOOP ENDS

!-----------------------------------------------------------------------
    J=1
    DO I=1,np(1)-1
    
       IP=I+1
       IM=I-1
    IF(I.EQ.1)IM=np(1)-1

    SUM1=0.0d0
    JACOB=ZI_X(I,J)*ETA_Y(I,J)-ZI_Y(I,J)*ETA_X(I,J)

    G_ETA(1)=-ZI_Y(I,J)/JACOB
    G_ETA(2)=+ZI_X(I,J)/JACOB
    DN=ETA
    
    SP(I,1)=0.0d0
    SP(I,2)=gw*(EPS1/(1+EPS1))/(FR*FR)
    IF(I.EQ.1) THEN
    SP(np(1),1)=SP(1,1)
    SP(np(1),2)=SP(1,2)
    ENDIF

    DO L=1,2
    DF_ZI=(F(IP,J,L)-F(IM,J,L))/(2*ZI)
    DF_ETA=(-3*F(I,J,L)+4*F(I,J+1,L)-F(I,J+2,L))/(2*ETA)
    DG_ZI=(G(IP,J,L)-G(IM,J,L))/(2*ZI)
    DG_ETA=(-3*G(I,J,L)+4*G(I,J+1,L)-G(I,J+2,L))/(2*ETA)

    DF_DX=ZI_X(I,J)*DF_ZI+ETA_X(I,J)*DF_ETA
    DG_DY=ZI_Y(I,J)*DG_ZI+ETA_Y(I,J)*DG_ETA

    SUM1=SUM1+(-DF_DX-DG_DY+SP(I,L))*G_ETA(L)
    END DO
    
    CORR=gw*DN*(1.d0/PK)*G_ETA(2)/((1+EPS1)*FR*FR)

    P(I,J)=(P(I,J+1)-DN*SUM1)/(1-CORR)
    ROW(I,J)=((P(I,J)/PK)+1)/VAR(I,J,3)

    END DO  ! I LOOP ENDS
!-------------------------PERIODICITY-----------------------------------
    P(np(1),1)=P(1,1)
    ROW(np(1),1)=ROW(1,1)
!-----------------------------------------------------------------------


     ELSE

     J=1
     DO I=1,np(1)-1
     ROW(I,J)=(4*ROW(I,J+1)-ROW(I,J+2))/3.0
!    IF(I.EQ.ITE)ROW(i,j)=2*ROW(i,j+1)-ROW(i,j+2) 

     UTP=VAR(I,J+1,1)*TS_X(I)+VAR(I,J+1,2)*TS_Y(I)
     UTPP=VAR(I,J+2,1)*TS_X(I)+VAR(I,J+2,2)*TS_Y(I)
           
     UTP=ROW(I,J+1)*UTP
     UTPP=ROW(I,J+2)*UTPP

     TERM=((4*UTP-UTPP)/3.0)/(TS_X(I)*NS_Y(I)-NS_X(I)*TS_Y(I))

     VAR(I,J,1)=NS_Y(I)*TERM/ROW(I,J)
     VAR(I,J,2)=-NS_X(I)*TERM/ROW(I,J)

!          IF(I.EQ.ITE)THEN
!           U(I,J)=0.0
!           V(I,J)=0.0
!          ENDIF

!    T(I,J)=(4*T(I,J+1)-T(I,J+2))/3.0!-T(I,3)

!     E(I,J)=T(I,J)+EK*(U(I,J)**2+V(I,J)**2)

     TERM=(4*U_VEC(I,J+1,4)-U_VEC(I,J+2,4))/3.0
     E(I,J)=TERM/ROW(I,J)
     SP_EN=E(I,J)-EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)
     VAR(I,J,3)=1.0d0+B1*(SP_EN-1)+B2*(SP_EN-1)**2+B3*(SP_EN-1)**3+B4*(SP_EN-1)**4+B5*(SP_EN-1)**5
     P(I,J)=PK*(ROW(I,J)*VAR(I,J,3)-1) 
  
     VIS(I,J)=0.d0
     COND(I,J)=0.d0
   END DO

   ENDIF
 
         VAR(np(1),J,1)=VAR(1,J,1)
         VAR(np(1),J,2)=VAR(1,J,2)
         VAR(np(1),J,3)=VAR(1,J,3)     
         ROW(np(1),J)=ROW(1,J)
         E(np(1),J)=E(1,J)
         P(np(1),J)=P(1,J)
         VIS(np(1),1)=VIS(1,1)
         COND(np(1),1)=COND(1,1)

      if(TVD_E=='yes')then

      IF(MODEL=='INVISCID')THEN
      J=1
      DO I=ITU4,ITL4
      ip=i+1
      im=i-1
      if(i.eq.1)im=np(1)-1
      
      VAR1=E(IP,J)
      VAR2=E(I,J)
      VAR3=E(IM,J)
      
      p_max=MAX(VAR1,VAR2,VAR3)
      p_min=MIN(VAR1,VAR2,VAR3)
      p_mean(i)=(VAR1+VAR2+VAR3)/3.0
      dev=dabs(E(I,J)-p_mean(i))
      dev_max=dabs(p_max-p_min)

      if(dev_max==0)Then
       TV(I)=0
      else
       TV(I)=dev/dev_max
      endif

      if(TV(I)>=0.1)THEN
       E(I,J)=p_mean(i)
       SP_EN=E(I,J)-EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)
       VAR(I,J,3)=1.0d0+B1*(SP_EN-1)+B2*(SP_EN-1)**2+B3*(SP_EN-1)**3+B4*(SP_EN-1)**4+B5*(SP_EN-1)**5
       ROW(I,J)=((P(I,J)/PK)+1)/VAR(I,J,3)
!       P(I,J)=PK*(ROW(I,J)*T(I,J)-1)
       U_VEC(I,J,1)=ROW(I,J)
       U_VEC(I,J,2)=ROW(I,J)*VAR(I,J,1)
       U_VEC(I,J,3)=ROW(I,J)*VAR(I,J,2)
       U_VEC(I,J,4)=ROW(I,J)*E(I,J)
      ENDIF

      END DO
       E(np(1),J)=E(1,J)
       ROW(np(1),J)=ROW(1,J)
       VAR(np(1),J,3)=VAR(1,J,3)

       U_VEC(np(1),J,1)=U_VEC(1,J,1)
       U_VEC(np(1),J,2)=U_VEC(1,J,2)
       U_VEC(np(1),J,3)=U_VEC(1,J,3)
       U_VEC(np(1),J,4)=U_VEC(1,J,4)
  
!     -------------------------------------------------------------------
!         correcting pressures for grid scale osccilations
!     -------------------------------------------------------------------
     
      J=1
      DO I=ITU4,ITL4
      ip=i+1
      im=i-1
      if(i.eq.1)im=np(1)-1
      
      VAR1=P(IP,J)
      VAR2=P(I,J)
      VAR3=P(IM,J)
      
      p_max=MAX(VAR1,VAR2,VAR3)
      p_min=MIN(VAR1,VAR2,VAR3)
      p_mean(i)=(VAR1+VAR2+VAR3)/3.0
      dev=dabs(P(I,J)-p_mean(i))
      dev_max=dabs(p_max-p_min)

      if(dev_max==0)Then
       TV(I)=0
      else
       TV(I)=dev/dev_max
      endif

      if(TV(I)>=0.1)THEN
       P(I,J)=p_mean(i)
!       T(I,J)=E(I,J)-EK*(U(I,J)**2+V(I,J)**2)
       ROW(I,J)=((P(I,J)/PK)+1)/VAR(I,J,3)
!       P(I,J)=PK*(ROW(I,J)*T(I,J)-1)
       U_VEC(I,J,1)=ROW(I,J)
       U_VEC(I,J,2)=ROW(I,J)*VAR(I,J,1)
       U_VEC(I,J,3)=ROW(I,J)*VAR(I,J,2)
       U_VEC(I,J,4)=ROW(I,J)*E(I,J)
      ENDIF

      END DO
       E(np(1),J)=E(1,J)
       ROW(np(1),J)=ROW(1,J)
       VAR(np(1),J,3)=VAR(1,J,3)
       P(np(1),J)=P(1,J)
       U_VEC(np(1),J,1)=U_VEC(1,J,1)
       U_VEC(np(1),J,2)=U_VEC(1,J,2)
       U_VEC(np(1),J,3)=U_VEC(1,J,3)
       U_VEC(np(1),J,4)=U_VEC(1,J,4)
      ENDIF      
      endif
  
! IF(MODEL.EQ.'VISCOUS') THEN

       

!---------Finding the maximum Mach no. in entire domain----------

!    MAXM=MACH
! DO 1-1
! DO J=bs,es
! DO I=1,np(1)-1
!    M_LOCAL=MACH*SQRT((VAR(I,J,1)**2+VAR(I,J,2)**2+**2)/VAR(I,J,3))
! IF(M_LOCAL.GE.MAXM)MAXM=M_LOCAL
! END DO
! END DO
! END DO 

!===========SETTING BOUNDARY CONDITIONS AT ARTIFICIAL SURFACE===========

       J=np(2)     
    DO I=1,np(1)-1
       IP=I+1
       IM=I-1
    IF (I.EQ.1)IM=np(1)-1
 
    V_N=VAR(I,J,1)*NA_X(I)+VAR(I,J,2)*NA_Y(I)
    V_T=VAR(I,J,1)*TA_X(I)+VAR(I,J,2)*TA_Y(I)
    
    VN(1)=VAR(I,J-1,1)*NA_X(I)+VAR(I,J-1,2)*NA_Y(I)
    VN(2)=VAR(I,J-2,1)*NA_X(I)+VAR(I,J-2,2)*NA_Y(I)
!   VN(3)=VAR(I,J-3,1)*NA_X(I)+VAR(I,J-3,2)*NA_Y(I)
!   VN(4)=VAR(I,J-4,1)*NA_X(I)+VAR(I,J-4,2)*NA_Y(I)

    VT(1)=VAR(I,J-1,1)*TA_X(I)+VAR(I,J-1,2)*TA_Y(I)
    VT(2)=VAR(I,J-2,1)*TA_X(I)+VAR(I,J-2,2)*TA_Y(I)
!   VT(3)=VAR(I,J-3,1)*TA_X(I)+VAR(I,J-3,2)*TA_Y(I)
!   VT(4)=VAR(I,J-4,1)*TA_X(I)+VAR(I,J-4,2)*TA_Y(I)

    SONIC=SQRT(VAR(I,J,3))/MACH
    
!======================INFLOW===========================================  
    IF(V_N.GT.0.d0) THEN
    
    MN=ABS(V_N)/SONIC

    IF ((V_N-SONIC).GE.0.d0) THEN
!-------------------ALL WAVES ARE ENTERING THE FLOW DOMAIN--------------
    VAR(I,J,1)=uinf
    VAR(I,J,2)=vinf
    ROW(I,J)=1.d0  !+A0*SIN((2.D0*PI/LAMBDA_Z)*Z()) ! PERTURBATION 
    VAR(I,J,3)=1.d0
    P(I,J)=(ROW(I,J)*VAR(I,J,3)-1)*PK
    E(I,J)=VAR(I,J,3)+APH1*(VAR(I,J,3)-1)**2+APH2*(VAR(I,J,3)-1)**3 &
               +APH3*(VAR(I,J,3)-1)**4+EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)

    ELSE
!      WRITE(*,*)I,'SUB-IN'

!-----------------ONE ACOUSTIC-WAVE LEAVES THE FLOW DOMAIN (K1)---------
    V_NINF=UINF*NA_X(I)+VINF*NA_Y(I)
    RNP=V_NINF+2.d0/(MACH*(GAMA-1))  ! FIXING THE VALUE EQUAL TO FREE-STREAM
       
    RN(1)=VN(1)-2*SQRT(VAR(I,J-1,3))/(MACH*(GAMA-1)) ! EXTRAPOLATING FROM INTERIOR
    RN(2)=VN(2)-2*SQRT(VAR(I,J-2,3))/(MACH*(GAMA-1))
 !  RN(3)=VN(3)-2*SQRT(VAR(I,J-3,3))/(MACH*(GAMA-1))
 !  RN(4)=VN(4)-2*SQRT(VAR(I,J-4,3))/(MACH*(GAMA-1))
     
    RNM=2*RN(1)-RN(2)  !RN(1)
    V_N=(RNP+RNM)/2
    VAR(I,J,3)=((RNP-RNM)*MACH*(GAMA-1)*0.25d0)**2
    ROW(I,J)=1.d0        !+A0*SIN((2.D0*PI/LAMBDA_Z)*Z()) ! PERTURBATION 
    V_T=UINF*TA_X(I)+VINF*TA_Y(I)   ! FIXING THE VALUE EQUAL TO FREE-STREAM

    TERM=NA_X(I)*TA_Y(I)-NA_Y(I)*TA_X(I)
    VAR(I,J,1)=(V_N*TA_Y(I)-V_T*NA_Y(I))/TERM
    VAR(I,J,2)=-(V_N*TA_X(I)-V_T*NA_X(I))/TERM
    P(I,J)=(ROW(I,J)*VAR(I,J,3)-1)*PK
    E(I,J)=VAR(I,J,3)+APH1*(VAR(I,J,3)-1)**2+APH2*(VAR(I,J,3)-1)**3 &
              +APH3*(VAR(I,J,3)-1)**4+EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)

    ENDIF
!--------------------------------------------
    ELSE    	
!=======================OUTFLOW=========================================
    IF ((V_N+SONIC).LE.0.d0) THEN
       
!----------------ALL WAVES ARE LEAVING THE FLOW DOMAIN-------------------------------
    VAR(I,J,1)=2*VAR(I,J-1,1)-VAR(I,J-2,1)
    VAR(I,J,2)=2*VAR(I,J-1,2)-VAR(I,J-2,2)
    VAR(I,J,3)=2*VAR(I,J-1,3)-VAR(I,J-2,3)
    ROW(I,J)  =2*ROW(I,J-1)-ROW(I,J-2)
    P(I,J)=(ROW(I,J)*VAR(I,J,3)-1)*PK
    E(I,J)=VAR(I,J,3)+APH1*(VAR(I,J,3)-1)**2+APH2*(VAR(I,J,3)-1)**3 &
               +APH3*(VAR(I,J,3)-1)**4+EK*(VAR(I,J,1)**2+VAR(I,J,2)**2) 
    ELSE
!C     WRITE(*,*)I,'SUB-OUT'
!-----------------ONE ACOUSTIC-WAVE ENTERS THE FLOW DOMAIN------------------------
    V_NOLD=V_N
    V_NINF=UINF*NA_X(I)+VINF*NA_Y(I)
!    RNM=-V_NINF-2.d0/(MACH*(GAMA-1))   ! FIXING THE VALUE EQUAL TO FREE-STREAM

    RNP=V_NINF+2.d0/(MACH*(GAMA-1))   ! FIXING THE VALUE EQUAL TO FREE-STREAM
    RN(1)=VN(1)-2*SQRT(VAR(I,J-1,3))/(MACH*(GAMA-1)) ! EXTRAPOLATING FROM INTERIOR
    RN(2)=VN(2)-2*SQRT(VAR(I,J-2,3))/(MACH*(GAMA-1))
!   RN(3)=-VN(3)+2*SQRT(VAR(I,J-3,3))/(MACH*(GAMA-1))
!   RN(4)=-VN(4)+2*SQRT(VAR(I,J-4,3))/(MACH*(GAMA-1))
       
!    RNP=2*RN(1)-RN(2) !RN(1)
    RNM=2*RN(1)-RN(2)
    V_N=(RNP+RNM)/2

    VAR(I,J,3)=VAR(I,J-1,3)
    V_T=2*VT(1)-VT(2)  !EXTRAPOLATING TANGENTIAL VELOCITY
    
!--------CHARACTERISTICS BC FOR PRESSURE AT OUTFLOW---------------------  
   ! DVN_DT=(V_N-V_NOLD)/DT
   ! RTDT=MACH*SQRT(VAR(I,J,3))/DT
   ! P(I,J)=(RTDT*P(I,J)-DVN_DT)/((DVN_DT/PK)+RTDT)

    P(I,J)=P(I,J)-(V_N-V_NOLD)/MACH
!-----------------------------------------------------------------------    
    ROW(I,J)=(P(I,J)/PK+1)/VAR(I,J,3)
    
    TERM=NA_X(I)*TA_Y(I)-NA_Y(I)*TA_X(I)
    VAR(I,J,1)=(V_N*TA_Y(I)-V_T*NA_Y(I))/TERM
    VAR(I,J,2)=-(V_N*TA_X(I)-V_T*NA_X(I))/TERM
    E(I,J)=VAR(I,J,3)+APH1*(VAR(I,J,3)-1)**2+APH2*(VAR(I,J,3)-1)**3 &
               +APH3*(VAR(I,J,3)-1)**4+EK*(VAR(I,J,1)**2+VAR(I,J,2)**2)            
    ENDIF 

    ENDIF !  INFLOW/OUTFLOW DECISION
!-----------------------------------------------------------------------     
    V_CONT(I,J,1)=VAR(I,J,1)*ZI_X(I,J)+VAR(I,J,2)*ZI_Y(I,J)
    V_CONT(I,J,2)=VAR(I,J,1)*ETA_X(I,J)+VAR(I,J,2)*ETA_Y(I,J)
 
!   IF(MODEL.EQ.'VISCOUS') THEN
    VIS(I,J)=(VAR(I,J,3)**1.5)*(Tr+(S_/T0))/(VAR(I,J,3)*Tr+(S_/T0))
    COND(I,J)=A1+A2*VAR(I,J,3)+A3*VAR(I,J,3)**2   
!   ENDIF

    END DO  !END OF BOUNDARY POINT LOOP (I LOOP)
     
!--------------------------PERIODICITY---------------------------------- 
    VAR(np(1),J,1:3)=VAR(1,J,1:3)
    ROW(np(1),J)=ROW(1,J)
    E(np(1),J)=E(1,J)
    P(np(1),J)=P(1,J)
    V_CONT(np(1),J,1:2)=V_CONT(1,J,1:2)

!   IF (MODEL.EQ.'VISCOUS') THEN
    VIS(np(1),J)=VIS(1,J)
    COND(np(1),J)=COND(1,J)
!   ENDIF 

!=======================================================================
      IF(MOD(nstep,100).EQ.0) THEN

!-----------------------------------------------------------------------
 OPEN(2,FILE="t_his5.DAT",STATUS='UNKNOWN',ACCESS='append')
  WRITE(2,'(6(e13.6,2x))')Time,VAR(246,143,1),VAR(246,143,2),VAR(246,143,3),row(246,143),p(246,143)
 CLOSE(2)
 !----------------------------------------------------------------------
 OPEN(2,FILE='t_his2.DAT',STATUS='UNKNOWN',ACCESS='append')
  WRITE(2,'(6(e13.6,2x))')Time,VAR(246,83,1),VAR(246,83,2),VAR(246,83,3),row(246,83),p(246,83)
 CLOSE(2)
!-----------------------------------------------------------------------

!===================Calculating the Force Coefficients==================
        
    PX=0.d0
    PY=0.d0
    PXY=0.d0
    VORTX=0.d0
    VORTY=0.d0
    VORTXY=0.d0
    DILATX=0.d0
    DILATY=0.d0
    DILATXY=0.d0
    HEATF=0.d0
               
       J=1
    DO I=1,np(1)
    
    IP=I+1
    IM=I-1
    IF(I.EQ.1)     IM=np(1)-1
    IF(I.EQ.np(1)) IP=2

    JACOB =ABS(ZI_X(I,J)*ETA_Y(I,J)-ZI_Y(I,J)*ETA_X(I,J))
    GRAETA=ETA_X(I,J)**2+ETA_Y(I,J)**2

    IF(I.GT.1) THEN
!--------------COEFFICIENTS DUE TO PRESSURE-----------------------------
    FUN2=P(I,J)*ETA_X(I,J)/JACOB
    FUN1=P(IM,J)*ETA_X(IM,J)/JACOBP
    PX=PX+0.5d0*ZI*(FUN2+FUN1)

    FUN2=P(I,J)*ETA_Y(I,J)/JACOB
    FUN1=P(IM,J)*ETA_Y(IM,J)/JACOBP
    PY=PY+0.5d0*ZI*(FUN2+FUN1)

    FUN2=P(I,J)*(X(I,J)*ETA_Y(I,J)-Y(I,J)*ETA_X(I,J))/JACOB
    FUN1=P(IM,J)*(X(IM,J)*ETA_Y(IM,J)-Y(IM,J)*ETA_X(IM,J))/JACOBP
    PXY=PXY+0.5d0*ZI*(FUN2+FUN1)
!-----------------------------------------------------------------------
    ENDIF
    
!   IF(MODEL.EQ.'VISCOUS') THEN
    DU_DZI=(VAR(IP,J,1)-VAR(IM,J,1))/(2.d0*ZI)
    DV_DZI=(VAR(IP,J,2)-VAR(IM,J,2))/(2.d0*ZI)
    DT_DZI=(VAR(IP,J,3)-VAR(IM,J,3))/(2.d0*ZI)

    DU_DETA=(VAR(I,J+1,1)-VAR(I,J,1))/ETA
    DV_DETA=(VAR(I,J+1,2)-VAR(I,J,2))/ETA
    DT_DETA=(VAR(I,J+1,3)-VAR(I,J,3))/ETA

    DU_DX=DU_DZI*ZI_X(I,J)+DU_DETA*ETA_X(I,J)
    DU_DY=DU_DZI*ZI_Y(I,J)+DU_DETA*ETA_Y(I,J)

    DV_DX=DV_DZI*ZI_X(I,J)+DV_DETA*ETA_X(I,J)
    DV_DY=DV_DZI*ZI_Y(I,J)+DV_DETA*ETA_Y(I,J)       

    DT_DX=DT_DZI*ZI_X(I,J)+DT_DETA*ETA_X(I,J)
    DT_DY=DT_DZI*ZI_Y(I,J)+DT_DETA*ETA_Y(I,J)

    VORT=DV_DX-DU_DY  !VORT_Z

    DILAT=DU_DX+DV_DY
    QFLUX=COND(I,J)*DT_DETA
  

    IF(I.GT.1) THEN
!--------------COEFFICIENTS DUE TO VORTICITY----------------------------
    FUN2=VIS(I,J)*VORT*ETA_Y(I,J)/JACOB
    FUN1=VIS(IM,J)*VORTP*ETA_Y(IM,J)/JACOBP
    VORTX=VORTX+(FUN2+FUN1)*0.5d0*ZI

    FUN2=VIS(I,J)*VORT*ETA_X(I,J)/JACOB
    FUN1=VIS(IM,J)*VORTP*ETA_X(IM,J)/JACOBP
    VORTY=VORTY+(FUN2+FUN1)*0.5d0*ZI

    FUN2=VIS(I,J)*VORT*(X(I,J)*ETA_X(I,J)+Y(I,J)*ETA_Y(I,J))/JACOB
    FUN1=VIS(IM,J)*VORTP*(X(IM,J)*ETA_X(IM,J)+Y(IM,J)*ETA_Y(IM,J))/JACOBP
    VORTXY=VORTXY+(FUN2+FUN1)*0.5d0*ZI
      
!--------------COEFFICIENTS DUE TO VOLUMETRIC STRAINING-----------------
    FUN2=VIS(I,J)*DILAT*ETA_X(I,J)/JACOB
    FUN1=VIS(IM,J)*DILATP*ETA_X(IM,J)/JACOBP
    DILATX=DILATX+(FUN2+FUN1)*0.5d0*ZI

    FUN2=VIS(I,J)*DILAT*ETA_Y(I,J)/JACOB
    FUN1=VIS(IM,J)*DILATP*ETA_Y(IM,J)/JACOBP
    DILATY=DILATY+(FUN2+FUN1)*0.5d0*ZI

    FUN2=VIS(I,J)*DILAT*(X(I,J)*ETA_Y(I,J)-Y(I,J)*ETA_X(I,J))/JACOB
    FUN1=VIS(IM,J)*DILATP*(X(IM,J)*ETA_Y(IM,J)-Y(IM,J)*ETA_X(IM,J))/JACOBP
    DILATXY=DILATXY+(FUN2+FUN1)*0.5d0*ZI
!-----------------------------------------------------------------------   

    FUN2=QFLUX*GRAETA/JACOB
    FUN1=QFLUXP*GRAETAP/JACOBP
    HEATF=HEATF+(FUN2+FUN1)*0.5d0*ZI 
    ENDIF
    
    VORTP=VORT
    DILATP=DILAT
    QFLUXP=QFLUX
!   ENDIF

   JACOBP=JACOB
   GRAETAP=GRAETA

   END DO

    CXP=-PX
    CYP=-PY
    CXYP=-PXY
   
    CXV=(-1.d0/RE)*VORTX
    CYV=(1.d0/RE)*VORTY
    CXYV=(1.d0/RE)*VORTXY

    CXD=(4.d0/(3*RE))*DILATX
    CYD=(4.d0/(3*RE))*DILATY
    CXYD=(4.d0/(3*RE))*DILATXY

    ANUSS=(-1.d0/(4.d0*EPS1))*HEATF
!-----------------------------------------------------------------------
   CX=CXP+CXV+CXD
   CY=CYP+CYV+CYD
     
    IF (EPS1.EQ.0)anuss=0.d0
    CD =2.d0*(CX*COS(THETAR)+CY*SIN(THETAR))
    CL =2.d0*(CY*COS(THETAR)-CX*SIN(THETAR))
    CDP=2.d0*(CXP*COS(THETAR)+CYP*SIN(THETAR))
    CDV=2.d0*(CXV*COS(THETAR)+CYV*SIN(THETAR))
    CDD=2.d0*(CXD*COS(THETAR)+CYD*SIN(THETAR))

!------------CALCULATING COEFFICIENT OF MOMENT--------------------------

    CMP=CXYP       !MOMENT CONTRIBUTION BY PRESSURE
    CMV=CXYV       !MOMENT CONTRIBUTION BY VORTICITY
    CMD=CXYD       !MOMENT CONTRIBUTION BY VOLUMETRIC STRAINING
    CM=CMP+CMV+CMD !TOTAL COEFFICIENT OF MOMENT
!-----------------------------------------------------------------------
 OPEN(2,FILE='COEFF_his.DAT',STATUS='UNKNOWN',ACCESS='append')
  WRITE(2,'(4(e13.6,2x))')Time,CD,CL,ANUSS
 CLOSE(2)

 OPEN(2,FILE='DRAG_COEFFS.DAT',STATUS='UNKNOWN',ACCESS='append')
   WRITE(2,'(4(e13.6,2x))')Time,CDP,CDV,CDD
 CLOSE(2)

 OPEN(2,FILE='MOMENT_COEFFS.DAT',STATUS='UNKNOWN',ACCESS='append')
   WRITE(2,'(5(e13.6,2x))')Time,CMP,CMV,CMD,CM
 CLOSE(2)
 
!======================SURFACE FILE====================================

 OPEN(7,FILE='SURFACE',STATUS='UNKNOWN')
  
       J=1
    DO I=1,np(1)
 
      IP=I+1
      IM=I-1
    IF (I.EQ.1)IM=np(1)-1
    IF (I.EQ.np(1))IP=2
 
    DU_DZI=(VAR(IP,J,1)-VAR(IM,J,1))/(2*ZI)
    DV_DZI=(VAR(IP,J,2)-VAR(IM,J,2))/(2*ZI)

    DU_DETA=(VAR(I,J+1,1)-VAR(I,J,1))/ETA
    DV_DETA=(VAR(I,J+1,2)-VAR(I,J,2))/ETA

    DU_DX=DU_DZI*ZI_X(I,J)+DU_DETA*ETA_X(I,J)
    DU_DY=DU_DZI*ZI_Y(I,J)+DU_DETA*ETA_Y(I,J)

    DV_DX=DV_DZI*ZI_X(I,J)+DV_DETA*ETA_X(I,J)
    DV_DY=DV_DZI*ZI_Y(I,J)+DV_DETA*ETA_Y(I,J)       

    DILAT=DU_DX+DV_DY
    VORT=DV_DX-DU_DY    
    SHEAR=DV_DX+DU_DY
       
!   Local Nusselt number profile on cylinder for Non-Boussinesq model
!-----------------------------------------------------------------------
    DT_DETA=(4*VAR(I,J+1,3)-3*VAR(I,J,3)-VAR(I,J+2,3))/ETA 
    Nu_local=(-1.d0/EPS1)*DT_DETA*SQRT(ETA_X(I,1)**2+ETA_Y(I,1)**2) 
    IF (EPS1.EQ.0)Nu_local=0.d0 
!-----------------------------------------------------------------------

    WRITE(7,20) (I-1)*ZI,P(I,1),VORT,DILAT,SHEAR,Nu_local
	    
    END DO
    CLOSE(7)
    
  20     FORMAT(1(F10.5,2X),5(E17.8,2X))
!-----------------------------------------------------------------------
    END IF   ! IF(MOD(nstep,1000).EQ.0) THEN
  
!==================write output files (starts)========================
 IF(mod(nstep,1000).EQ.0) THEN
      
    OPEN(2,FILE='VORT_DILA_SHEAR',STATUS='UNKNOWN')
          WRITE(2,*)'ZONE'
!write(2,*)'VARIABLES=,"X","Y","VORT-Z","DIV","SHEAR"'
          WRITE(2,*)'J=',NP(2),'I=',NP(1)

        DO J=1,NP(2)
        DO I=1,NP(1)

        IP=I+1
        IM=I-1

        IF(I.EQ.1)IM=NP(1)-1
        IF(I.EQ.NP(1))IP=2
        DU_DZI=(VAR(IP,J,1)-VAR(IM,J,1))/(2*ZI)
        DV_DZI=(VAR(IP,J,2)-VAR(IM,J,2))/(2*ZI)

        IF(J.EQ.1)THEN
        DU_DETA=(VAR(I,J+1,1)-VAR(I,J,1))/(ETA)
        DV_DETA=(VAR(I,J+1,2)-VAR(I,J,2))/(ETA)
        ENDIF

        IF(J.NE.1.AND.J.NE.NP(2))THEN
        DU_DETA=(VAR(I,J+1,1)-VAR(I,J-1,1))/(2*ETA)
        DV_DETA=(VAR(I,J+1,2)-VAR(I,J-1,2))/(2*ETA)
        ENDIF

        IF(J.EQ.NP(2))THEN
        DU_DETA=(VAR(I,J,1)-VAR(I,J-1,1))/(ETA)
        DV_DETA=(VAR(I,J,2)-VAR(I,J-1,2))/(ETA)
        ENDIF

        DU_DX=DU_DZI*ZI_X(I,J)+DU_DETA*ETA_X(I,J)
        DU_DY=DU_DZI*ZI_Y(I,J)+DU_DETA*ETA_Y(I,J)

        DV_DX=DV_DZI*ZI_X(I,J)+DV_DETA*ETA_X(I,J)
        DV_DY=DV_DZI*ZI_Y(I,J)+DV_DETA*ETA_Y(I,J)

        VORT=DV_DX-DU_DY
        DILAT=DU_DX+DV_DY
        SHEAR=DV_DX+DU_DY

        WRITE(2,24)X(I,J),Y(I,J),VORT,DILAT,SHEAR

        END DO
        WRITE(2,*)

        END DO
        CLOSE(2)

 OPEN(7,FILE='OUTPUT',STATUS='unknown')
 WRITE(7,103) 'TITLE ="',Time,nstep,'"'
 WRITE(7,*) 'VARIABLES=,"X","Y","ROW","P","T","U","V","E","ML"'
 WRITE(7,104) 'ZONE J=',np(2),',I=',np(1),',DATAPACKING="POINT"'
  
 DO J=1,np(2) !bs,es
 WRITE (7,*)
 DO I=1,np(1)
 
    FV=SQRT(VAR(I,J,1)**2+VAR(I,J,2)**2)
    ML=MACH*FV/SQRT(VAR(I,J,3))
  
!~    V_CONT(I,J,1)=VAR(I,J,1)*ZI_X(I,J)+VAR(I,J,2)*ZI_Y(I,J)        
!~    V_CONT(I,J,2)=VAR(I,J,1)*ETA_X(I,J)+VAR(I,J,2)*ETA_Y(I,J)
!~ 
!~    ML(1) =MACH*V_CONT(I,J,1)/SQRT(VAR(I,J,3))
!~    ML(2) =MACH*V_CONT(I,J,2)/SQRT(VAR(I,J,3))

    WRITE(7,22) X(I,J),Y(I,J),ROW(I,J),P(I,J),VAR(I,J,3),VAR(I,J,1),VAR(I,J,2),E(I,J),ML

    END DO	    
    END DO
    CLOSE(7)

 103	format(A,F20.10,I9,A)
 104	format(A,I3,A,I3,A,I3,A)
 22     FORMAT(2(F10.5,2X),7(E17.8,2X))
 24      FORMAT(2(F10.5,2X),3(E16.8,2X))
 
 END IF
!==================Write output files (ends)==========================
!--------------------------------------------------------
 !IF(mod(nstep,100)==0) THEN
 !CALL cpu_time(tEnd)
 !WRITE(6,114) nstep,time," | CPU Time:",tEnd-tStart,"sec for 100 iters | "
 !WRITE(6,'(A)')"-----------------------------------------------------"
 !ENDIF
!114   format(I10,F10.5,A,F6.3,2(1X,A))
	
 END DO
!========================Main time loop (ends)==========================
	     STOP
	     END
