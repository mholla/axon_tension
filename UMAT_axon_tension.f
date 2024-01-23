************************************************************************
!
! The umat for the axonal tension project
! 
************************************************************************
!     State Variables
!     --------------------------------------------------------------
!      statev(1) = growth parameters 
! 
!     Material Properties Vector
!     --------------------------------------------------------------
!      !!!!!Cortex!!!! 
!      mu_c      = props(3) shear modulus for cortex
!      lambda_c  = props(4) lame constant for cortex
!      Gctx      = props(5) rate constant for cortex growth 
!      !!!!Subcortex!!! 
!      mu_s        = props(1) shear modulus for subcortex
!      lambda_s    = props(2) lame constant for subcortex
***********************************************************************
***********************************************************************
      !--------------------------------------------------------
      ! initialize the growth variable 
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)

      include 'aba_param.inc'

      dimension statev(nstatv)
      
      ! growth_variable
      statev(1)=1.0d0    
          
      return
      end
      !--------------------------------------------------------
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

      include 'aba_param.inc'

      dimension stress(ntens),statev(nstatv),
     + field(1),direct(3,3),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      character*80 cmname
      character*256 CORTEX,SUBCORTEX


      !--------------------------------------------------------
      ! call particular user material to perform the analysis 
      IF (CMNAME(1:9) .EQ. 'SUBCORTEX') THEN

      call umat_subcortex(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

      ELSE IF(CMNAME(1:6) .EQ. 'CORTEX') THEN

      call umat_cortex(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

      Endif
      
      RETURN
      END subroutine umat
***********************************************************************
      subroutine umat_subcortex(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)

      include 'aba_param.inc'

      dimension stress(ntens),statev(nstatv),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      character*80 cmname
      integer i,j,k,l

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),T_tau(3,3)
      real*8 Finv(3,3),Fginv(3,3),C_tau(3,3)
      real*8 Be_tau(3,3),B_tau(3,3),Fg_tau(3,3),Fe_tau(3,3),Je,Jg,detF
      real*8 lamg_t,lamg_tau,jac(3,3,3,3),mu_s,lambda_s,Ce(3,3,3,3),Cs(3,3,3,3)


      ! Parameters
      real*8 zero,one,two,half,three,third,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0)

      ! Obtain old and new deformation gradients
      F_t = dfgrd0
      F_tau = dfgrd1

      ! Identity matrix
      call onem(Iden)

      ! Compute the relative volume change
      call mdet(F_tau,detF)

      ! Compute the inverse of the deformation gradient
      call m3inv(F_tau,Finv)

      ! Compute the right Cauchy Green tensor 
      C_tau = matmul(transpose(F_tau),F_tau)

      ! Obtain material properties 
      mu_s        = props(1)
      lambda_s    = props(2) 

      ! Subcortex does not grow
      lamg_tau = 1.0

      ! Update  kinematics 
      Fg_tau  = Iden
      B_tau = matmul(F_tau, transpose(F_tau))
      Be_tau = B_tau
 
      call mdet(Fg_tau,Jg)
      Je = detF/Jg

      call m3inv(Fg_tau,Fginv)     
      call m3inv(C_tau,Cinv)

      ! Compute Cauchy stress 
      T_tau = ((lambda_s*dlog(Je) - mu_s)*Iden  + mu_s*Be_tau)/Je

      ! Compute the analytical jacobian 
      Ce = zero 
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3            
                  Ce(i,j,k,l) =  
     +             +(1.0/Je)*lambda_s*Iden(i,j)*Iden(k,l)
     +             +(1.0/Je)*(mu_s - lambda_s*dlog(Je))*
     +             (Iden(i,k)*Iden(j,l) + Iden(j,k)*Iden(i,l))
             enddo
           enddo
         enddo
      enddo

      ! Compute the spatial tangent 
      Cs = zero  
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3            
                  Cs(i,j,k,l) =  
     +             + 0.5*(
     +             + Iden(i,k)*T_tau(j,l)
     +             + Iden(i,l)*T_tau(j,k)
     +             + Iden(j,k)*T_tau(i,l)
     +             + Iden(j,l)*T_tau(i,k)
     +                   )
             enddo
           enddo
         enddo
      enddo      

      jac  =  Ce + Cs

      ! Update state variables
      statev(1)  = lamg_tau

      ! Return Abaqus/Standard the Cauchy stress
      if(ntens.eq.6) then
         !
         ! 3D problem
         !
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
         stress(5) = T_tau(1,3)
         stress(6) = T_tau(2,3)
      elseif(ntens.eq.4) then
         !
         ! 2D problem
         !
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
      endif

      ! Return Abaqus/Standard the stress-deformation jacobian
      if(ntens.eq.6) then
         call jac3D(jac,ddsdde)
      elseif(ntens.eq.4) then
         call jac2D(jac,ddsdde)
      endif


      return
      end subroutine umat_subcortex

***********************************************************************
      subroutine umat_cortex(stress,statev,ddsdde,sse,spd,scd,
     + rpl,ddsddt,drplde,drpldt,
     + stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     + ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     + celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)


      include 'aba_param.inc'

      dimension stress(ntens),statev(nstatv),
     + ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     + stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     + props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      character*80 cmname

      integer i,j,k,l

      real*8 Iden(3,3),F_t(3,3),F_tau(3,3),T_tau(3,3)
      real*8 Finv(3,3),B_tau(3,3)
      real*8 lambda_c,mu_c,Be_tau(3,3),Fg_tau(3,3),Fe_tau(3,3),Je
      real*8 thetag_tau,thetag_t
      real*8 dtime,Jg,C_tau(3,3)
      real*8 Ce(3,3,3,3),Cg(3,3,3,3),Cs(3,3,3,3)
      real*8 jac(3,3,3,3),Gctx
      real*8 Fg_inv_tau(3,3),Ce_tau(3,3),TrCe,sse

      ! Parameters
      real*8 zero,one,two,half,three,third,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0)

      ! Obtain old and new deformation gradients
      F_t = dfgrd0
      F_tau = dfgrd1

      ! Identity matrix
      call onem(Iden)

      ! Compute the relative volume change
      call mdet(F_tau,detF)

      ! Compute the inverse of the deformation gradient
      call m3inv(F_tau,Finv)

      ! Compute the right Cauchy Green tensor 
      C_tau = matmul(transpose(F_tau),F_tau)

      ! Obtain material properties 
      mu_c      = props(3)
      lambda_c  = props(4)
      Gctx      = props(5)

      ! Time dependent rate 
      thetag_t = one
      thetag_tau = thetag_t + Gctx*(time(2)+dtime)     

      ! Update  kinematics 
      Fg_tau  = thetag_tau*Iden

      call m3inv(Fg_tau,Fg_inv_tau)
      
      Fe_tau = matmul(F_tau,Fg_inv_tau)

      Be_tau = matmul(Fe_tau,transpose(Fe_tau))

      Ce_tau = matmul(transpose(Fe_tau),Fe_tau)

      TrCe =  Ce_tau(1,1) + Ce_tau(2,2) + Ce_tau(3,3) 

      call mdet(Fg_tau,Jg)

      Je = detF/Jg
      
      ! Compute Cauchy stress 
      T_tau = ((lambda_c*dlog(Je) - mu_c)*Iden  + mu_c*Be_tau)/Je

      ! Compute the analytical jacobian 
      Ce = zero 
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3            
                  Ce(i,j,k,l) =  
     +             +(1.0/Je)*lambda_c*Iden(i,j)*Iden(k,l)
     +             +(1.0/Je)*(mu_c - lambda_c*dlog(Je))*
     +             (Iden(i,k)*Iden(j,l) + Iden(j,k)*Iden(i,l))
             enddo
           enddo
         enddo
      enddo
 
      ! Compute the spatial tangent 
      Cs = zero  
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3            
                  Cs(i,j,k,l) = 
     +             + 0.5*(
     +             + Iden(i,k)*T_tau(j,l)
     +             + Iden(i,l)*T_tau(j,k)
     +             + Iden(j,k)*T_tau(i,l)
     +             + Iden(j,l)*T_tau(i,k)
     +                   )
             enddo
           enddo
         enddo
      enddo      

      jac  =  Ce + Cs

      ! Compute the strain energy
      sse = (mu_c/2.0)*(TrCe-3.0-2.0*dlog(Je)) + (lambda_c)*(dlog(Je))**2.0/2.0

      ! Update state variables
      statev(1)  = thetag_tau

      ! Return Abaqus/Standard the Cauchy stress
      if(ntens.eq.6) then
         !
         ! 3D problem
         !
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
         stress(5) = T_tau(1,3)
         stress(6) = T_tau(2,3)
      elseif(ntens.eq.4) then
         !
         ! 2D problem
         !
         stress(1) = T_tau(1,1)
         stress(2) = T_tau(2,2)
         stress(3) = T_tau(3,3)
         stress(4) = T_tau(1,2)
      endif

      ! Return Abaqus/Standard the stress-deformation jacobian
      if(ntens.eq.6) then
         call jac3D(jac,ddsdde)
      elseif(ntens.eq.4) then
         call jac2D(jac,ddsdde)
      endif

      return
      end subroutine umat_cortex
****************************************************************************

      subroutine jac2D(SpTanMod,ddsdde)

      real*8 SpTanMod(3,3,3,3),ddsdde(4,4)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)

      end subroutine jac2D

***********************************************************************

      subroutine jac3D(SpTanMod,ddsdde)

      real*8 SpTanMod(3,3,3,3),ddsdde(6,6)

      ddsdde(1,1) = SpTanMod(1,1,1,1)
      ddsdde(1,2) = SpTanMod(1,1,2,2)
      ddsdde(1,3) = SpTanMod(1,1,3,3)
      ddsdde(1,4) = SpTanMod(1,1,1,2)
      ddsdde(1,5) = SpTanMod(1,1,1,3)
      ddsdde(1,6) = SpTanmod(1,1,2,3)

      ddsdde(2,1) = SpTanMod(2,2,1,1)
      ddsdde(2,2) = SpTanMod(2,2,2,2)
      ddsdde(2,3) = SpTanMod(2,2,3,3)
      ddsdde(2,4) = SpTanMod(2,2,1,2)
      ddsdde(2,5) = SpTanMod(2,2,1,3)
      ddsdde(2,6) = SpTanmod(2,2,2,3)

      ddsdde(3,1) = SpTanMod(3,3,1,1)
      ddsdde(3,2) = SpTanMod(3,3,2,2)
      ddsdde(3,3) = SpTanMod(3,3,3,3)
      ddsdde(3,4) = SpTanMod(3,3,1,2)
      ddsdde(3,5) = SpTanMod(3,3,1,3)
      ddsdde(3,6) = SpTanmod(3,3,2,3)

      ddsdde(4,1) = SpTanMod(1,2,1,1)
      ddsdde(4,2) = SpTanMod(1,2,2,2)
      ddsdde(4,3) = SpTanMod(1,2,3,3)
      ddsdde(4,4) = SpTanMod(1,2,1,2)
      ddsdde(4,5) = SpTanMod(1,2,1,3)
      ddsdde(4,6) = SpTanmod(1,2,2,3)

      ddsdde(5,1) = SpTanMod(1,3,1,1)
      ddsdde(5,2) = SpTanMod(1,3,2,2)
      ddsdde(5,3) = SpTanMod(1,3,3,3)
      ddsdde(5,4) = SpTanMod(1,3,1,2)
      ddsdde(5,5) = SpTanMod(1,3,1,3)
      ddsdde(5,6) = SpTanmod(1,3,2,3)

      ddsdde(6,1) = SpTanMod(2,3,1,1)
      ddsdde(6,2) = SpTanMod(2,3,2,2)
      ddsdde(6,3) = SpTanMod(2,3,3,3)
      ddsdde(6,4) = SpTanMod(2,3,1,2)
      ddsdde(6,5) = SpTanMod(2,3,1,3)
      ddsdde(6,6) = SpTanmod(2,3,2,3)

      end subroutine jac3D

C**********************************************************************
      SUBROUTINE ONEM(A)

C     THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C     3 BY 3 MATRIX [A]
C**********************************************************************

        REAL*8 A(3,3)
        DATA ZERO/0.D0/
        DATA ONE/1.D0/

      DO 1 I=1,3
        DO 1 J=1,3
          IF (I .EQ. J) THEN
              A(I,J) = 1.0
            ELSE
              A(I,J) = 0.0
            ENDIF
1       CONTINUE

      RETURN
      END
C**********************************************************************
      SUBROUTINE MDET(A,DET)
 
C     THIS SUBROUTINE CALCULATES THE DETERMINANT
C     OF A 3 BY 3 MATRIX [A].
C**********************************************************************

      REAL*8  A(3,3), DET

      DET =   A(1,1)*A(2,2)*A(3,3) 
     +              + A(1,2)*A(2,3)*A(3,1)
     +              + A(1,3)*A(2,1)*A(3,2)
     +            - A(3,1)*A(2,2)*A(1,3)
     +            - A(3,2)*A(2,3)*A(1,1)
     +            - A(3,3)*A(2,1)*A(1,2)

      RETURN
      END
C**********************************************************************
      SUBROUTINE M3INV(A,AINV)

C     THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C     [A] AND PLACES THE RESULT IN [AINV]. 
C     IF DET(A) IS ZERO, THE CALCULATION
C     IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

      REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C     A(3,3)              -- THE MATRIX WHOSE INVERSE IS DESIRED.
C     DET         -- THE COMPUTED DETERMINANT OF [A].
C     ACOFAC(3,3) -- THE MATRIX OF COFACTORS OF A(I,J).
C                    THE SIGNED MINOR (-1)**(I+J)*M_IJ
C                    IS CALLED THE COFACTOR OF A(I,J).
C     AADJ(3,3)   -- THE ADJOINT OF [A]. IT IS THE MATRIX
C                    OBTAINED BY REPLACING EACH ELEMENT OF
C                    [A] BY ITS COFACTOR, AND THEN TAKING
C                    TRANSPOSE OF THE RESULTING MATRIX.
C     AINV(3,3)   -- RETURNED AS INVERSE OF [A].
C                    [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------

      CALL MDET(A,DET)
      IF ( DET .EQ. 0.D0 ) THEN
        write(*,10)
        STOP
      ENDIF
      CALL MCOFAC(A,ACOFAC)
      CALL MTRANS(ACOFAC,AADJ)
      DO 1 I = 1,3
      DO 1 J = 1,3
           AINV(I,J) = AADJ(I,J)/DET
1     CONTINUE
10    FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +         10X,'PROGRAM TERMINATED')

      RETURN
      END

C**********************************************************************
      SUBROUTINE MCOFAC(A,ACOFAC)
 
C  THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C  AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************

      REAL*8  A(3,3), ACOFAC(3,3)

      ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
      ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
      ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
      ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
      ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
      ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
      ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
      ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
      ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

      RETURN
      END

C**********************************************************************
      SUBROUTINE MTRANS(A,ATRANS)
 
C  THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C  MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

      REAL*8 A(3,3),ATRANS(3,3)

      DO 1 I=1,3
         DO 1 J=1,3
            ATRANS(J,I) = A(I,J)
1     CONTINUE

      RETURN
      END

C**********************************************************************



