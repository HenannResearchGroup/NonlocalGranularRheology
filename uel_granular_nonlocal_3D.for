!************************************************************************
! User element subroutine (UEL) for large-deformation, steady-state, 
!  nonlocal granular flow in three-dimensions.
!************************************************************************
! Element details:
!************************************************************************
! 
! This subroutine is for a three-dimensional 8-node isoparametric
!  brick element as shown below with 8pt (full) integration.
!
! Solution variables (or nodal variables) are the displacements (DOFs 1-3)
!  and the granular fluidity (DOF 11).
!
! Material behavior is the nonlocal granular rheology or local inertial
!  rheology, specialized for steady-state flow.
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
! Mechanical, traction- and pressure-type boundary conditions 
!  may be applied to the dummy mesh using the Abaqus built-in 
!  commands *Dload or *Dsload.  Mechanical body forces, such
!  as gravity, may also be applied to the dummy mesh using the
!  Abaqus built-in command *Dload.
!
! Non-homogeneous flux-type boundary conditions for the granular
!  fluidity are not supported in this element. Energy outputs 
!  are not supported in this element.
!
! Element schematic:
!
!  8-node     8-----------7
!  brick     /|          /|       xi_3
!           / |         / |       
!          5-----------6  |       |     xi_2
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi_1
!          1-----------2        origin at cube center
!
! David L. Henann, April 2015
!
!************************************************************************
! Usage:
!************************************************************************
!
! User element statement in the input file:
!  *User Element,Nodes=8,Type=U1,Iproperties=1,Properties=7,Coordinates=3,Variables=72,Unsymm
!  1,2,3,11
!
! --------------------------------------------------------------
! State Variables
! --------------------------------------------------------------
! Global SDV's (for use in UVARM to make contour plots)
!   1) Stress ratio (mu)
!   2) Equivalent shear stress (tau)
!   3) Pressure (press)
!   4) Plastic shear strain rate (gammadotp)
!   5) Velocity - x-component - (v1)
!   6) Velocity - y-component - (v2)
!   7) Velocity - z-component - (v3)
!   8) Velocity - magnitude (vmag)
!
! Local SDV's (for use internal to the UEL)
!   j = 0
!   do k = 1,nIntPt
!      svars(1) = Fp(1,1) ---- Fp(1,1)
!      svars(2) = Fp(1,2) ---- Fp(1,2)
!      svars(3) = Fp(1,3) ---- Fp(1,3)
!      svars(4) = Fp(2,1) ---- Fp(2,1)
!      svars(5) = Fp(2,2) ---- Fp(2,2)
!      svars(6) = Fp(2,3) ---- Fp(2,3)
!      svars(7) = Fp(3,1) ---- Fp(3,1)
!      svars(8) = Fp(3,2) ---- Fp(3,2)
!      svars(9) = Fp(3,3) ---- Fp(3,3)
!      j = j + nlSdv
!   end loop over k
!
! In the input file, set '*User output variables' = 8
!
! --------------------------------------------------------------
! Material Properties Vector
! --------------------------------------------------------------
!  Gshear = props(1)  ! Shear modulus
!  Kbulk  = props(2)  ! Bulk modulus
!  rhos   = props(3)  ! Grain density
!  d      = props(4)  ! Particle diameter
!  mus    = props(5)  ! Static friction coeff
!  b      = props(6)  ! Rheology parameter
!  A      = props(7)  ! Nonlocal amplitude
!
!  matFlag = jprops(1) ! Nonlocal = 1, Local = 2
!
! In the input file, set the parameter matflag=1 for the steady-state 
!  nonlocal rheology or matflag=2 for local rheology.
!
! In the input file, the element type of the dummy mesh should be C3D8.
!
!************************************************************************

      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points. You must set that parameter value here.
      !
      !  offset
      !   Offset in the real mesh numbers and the dummy mesh 
      !   numbers. You must set that parameter value here.
      !
      !  elemPt
      !   Element pointer
      !

      integer numElem,offset,elemPt,elCount,err
      parameter(numElem=60,offset=60)

      real*8, allocatable :: globalSdv(:,:,:)

      
      end module global

!***********************************************************************

      subroutine UVARM(uvar,direct,t,time,dtime,cmname,orname,
     + nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord,
     + jmac,jmatyp,matlayo,laccfla)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.
     
      use global
     
      include 'ABA_PARAM.INC'

      character*80 cmname,orname
      character*3 flgray(15)
      dimension uvar(nuvarm),direct(3,3),t(3,3),time(2)
      dimension array(15),jarray(15),jmac(*),jmatyp(*),coord(*)

      ! The dimensions of the variables flgray, array and jarray
      !  must be set equal to or greater than 15.
      
      elCount = noel - offset


      uvar(1) = globalSdv(elCount,npt,1)
      uvar(2) = globalSdv(elCount,npt,2)
      uvar(3) = globalSdv(elCount,npt,3)
      uvar(4) = globalSdv(elCount,npt,4)
      uvar(5) = globalSdv(elCount,npt,5)
      uvar(6) = globalSdv(elCount,npt,6)
      uvar(7) = globalSdv(elCount,npt,7)
      uvar(8) = globalSdv(elCount,npt,8)


      return
      end subroutine uvarm

!****************************************************************************

      subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     +     props,nprops,coords,mcrd,nnode,uall,duall,vel,accn,jtype,
     +     time,dtime,kstep,kinc,jelem,params,ndload,jdltyp,adlmag,
     +     predef,npredf,lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,
     +     njprop,period)

      use global
      !
      implicit none
      !
      ! variables defined in uel, passed back to Abaqus
      !
      real*8 rhs(mlvarx,*),amatrx(ndofel,ndofel),svars(*),energy(8),
     +  pnewdt
      !
      ! variables passed into UEL
      !
      integer ndofel,nrhs,nsvars,nprops,mcrd,nnode,jtype,kstep,kinc,
     +  jelem,ndload,jdltyp(mdload,*),npredf,lflags(*),mlvarx,mdload,
     +  jprops(*),njprop
      !
      real*8 props(*),coords(mcrd,nnode),uall(ndofel),duall(mlvarx,*),
     +  vel(ndofel),accn(ndofel),time(2),dtime,params(*),
     +  adlmag(mdload,*),predef(2,npredf,nnode),ddlmag(mdload,*),period
      !
      ! variables defined and used in the UEL
      !
      integer i,j,k,l,jj,A11,A12,B11,B12,nDim,nlSdv,ngSdv,nInt,stat,
     +  matFlag,nIntPt,intpt
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      parameter(nDim=3)  ! number of spatial dimensions, do not change
      parameter(nlSdv=9) ! number of local state variables per integ pt, do not change
      parameter(ngSdv=8) ! number of global state variables per integ pt, do not change
      parameter(nInt=8)  ! number of integration points, do not change
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      real*8 Iden(3,3),Ru(3*nNode,1),Rg(nNode,1),Kuu(3*nNode,3*nNode),
     +  Kug(3*nNode,nNode),Kgu(nNode,3*nNode),Kgg(nNode,nNode),
     +  u(nNode,3),v(nNode,3),g(nNode),coordsC(mcrd,nNode),sh0(nNode),
     +  dshxi(nNode,3),dsh0(nNode,3),detMapJ0,dshC0(nNode,3),detMapJ0C,
     +  F0_tau(3,3),detF0_tau,xi(nInt,3),w(nInt),Fp_t(3,3),sh(nNode),
     +  dsh(nNode,3),detMapJ,dshC(nNode,3),detMapJC,F_tau(3,3),detF_tau,
     +  v_tau(3),v_mag,g_tau,dgdx(3,1),T_tau(3,3),gr(1,1),Fp_tau(3,3),
     +  mu_tau,taubar_tau,press,gammadotp_tau,Ctang(3,3,3,3),
     +  Lambda(3,3),Gamma(3,3),dgrdg(1,1),Smat(6,1),Bmat(6,3*nNode),
     +  Nmatg(1,nNode),Gmatg(3,nNode),Gmat(9,3*nNode),Gmat0(9,3*nNode),
     +  Cmat(9,9),Qmat(9,9),LambdaMat(9,1),Atang(3,3,3),Amat(3,9),
     +  GammaMat(1,9),Qmatg(1,9)
      !
      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,
     +     Pi=3.141592653d0,three=3.d0,third=1.d0/3.d0)


      ! Check the procedure type, this should be a coupled
      !  temperature displacement, which is either 72 or 73
      !
      if((lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! correct procedure specified
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and check the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear perturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear perturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear perturbation step'
         call xit         
      endif


      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.zero) return


      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF ELEMENTS -----------'
         write(*,*) '---------- numElem=',numElem
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF POINTS -------------'
         write(*,*) '---------- nInt=',nInt
         write(*,*) '-------------------------------------------------'
         write(*,*) '---------- YOU PUT NUMBER OF SDVs ---------------'
         write(*,*) '---------- ngSdv=',ngSdv
         write(*,*) '-------------------------------------------------'
      endif


      ! Counting for the SDV output
      !
      elemPt = jelem


      ! Get flag for material behavior
      !
      matFlag = jprops(1)


      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rg = zero
      Kuu = zero
      Kug = zero
      Kgu = zero
      Kgg = zero
      Energy = zero


      ! Obtain nodal displacements, velocities, and granular fluidities.
      !  (The velocity is approximated.)
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            v(i,j) = DUall(k,1)/dtime
         enddo
         k = k + 1
         g(i) = Uall(k)
      enddo


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures 33, 3277-3296.
      !
      ! Obtain shape functions and their local gradients at the element
      !  centroid, that means xi_1=xi_2=xi_3=0.0, and nIntPt=1
      !
      if(nNode.eq.8) then
         call calcShape3DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      call mapShape3D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Map shape functions from local to global current coordinate system
      !
      call mapShape3D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Calculate the deformation gradient at the element centroid
      !  at the end of the increment for use in the `F-bar' method.
      !  The subscript tau denotes the time at the end of the increment.
      !
      F0_tau = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               F0_tau(i,j) = F0_tau(i,j) + dsh0(k,j)*u(k,i)
            enddo
         enddo
      enddo
      !
      call mdet(F0_tau,detF0_tau)
      !
      ! With the deformation gradient known at the element centroid
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.8) then
         !
         ! Gauss integration for a brick element
         !
         if(nInt.eq.8) then
            call xint3D8pt(xi,w,nIntPt) ! 8-pt integration, nInt=8 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Loop over integration points
      !
      jj = 0 ! jj is used for tracking the state variables
      do intpt=1,nIntPt


         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! This is the first increment of the first step.
            !  Give initial conditions.
            !
            Fp_t = Iden
            !
         else
            !
            ! This is not the first increment; read old values.
            !
            Fp_t(1,1) = svars(1+jj)
            Fp_t(2,1) = svars(2+jj)
            Fp_t(3,1) = svars(3+jj)
            Fp_t(1,2) = svars(4+jj)
            Fp_t(2,2) = svars(5+jj)
            Fp_t(3,2) = svars(6+jj)
            Fp_t(1,3) = svars(7+jj)
            Fp_t(2,3) = svars(8+jj)
            Fp_t(3,3) = svars(9+jj)
            !
         endif


         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.8) then
            call calcShape3DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.8'
            call xit
         endif


         ! Map shape functions from local to global reference coordinate system
         !
         call mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Map shape functions from local to global current coordinate system
         !
         call mapShape3D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Obtain the deformation gradient at this integration point.
         !  The subscript tau denotes the time at the end of the increment.
         !
         F_tau = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
               enddo
            enddo
         enddo
         !
         ! Modify the deformation gradient for the `F-bar' method.
         !
         call mdet(F_tau,detF_tau)
         F_tau = ((detF0_tau/detF_tau)**third)*F_tau
            

         ! Obtain the (approximate) velocity at this integration point.
         !
         v_tau = zero
         do i=1,nDim
            do k=1,nNode
               v_tau(i) = v_tau(i) + v(k,i)*sh(k)
            end do
         end do
         !
         v_mag = dsqrt(sum(v_tau*v_tau))


         ! Obtain the granular fluidity and its spatial gradient at 
         !  this integration point.  
         !
         g_tau = zero
         dgdx = zero
         do k=1,nNode
            g_tau = g_tau + g(k)*sh(k)
            do i=1,nDim
               dgdx(i,1) = dgdx(i,1) + dshC(k,i)*g(k)
            end do
         end do


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive update at this integ. point
         !
         if (matflag.eq.1) then
            call nonlocal(props,nprops,dtime,F_tau,Fp_t,g_tau,
     +                 T_tau,gr,Fp_tau,mu_tau,taubar_tau,press,
     +                 gammadotp_tau,Ctang,Lambda,Gamma,dgrdg,stat)
         elseif (matflag.eq.2) then
            call local(props,nprops,dtime,F_tau,Fp_t,g_tau,
     +                 T_tau,gr,Fp_tau,mu_tau,taubar_tau,press,
     +                 gammadotp_tau,Ctang,Lambda,Gamma,dgrdg,stat)
         else
            write(*,*) 'Invalid matflag: matflag.ne.1 or 2'
            call xit
         endif
         !
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


         ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(elemPt,intPt,1) = mu_tau        ! stress ratio
         globalSdv(elemPt,intPt,2) = taubar_tau    ! equiv. shear stress
         globalSdv(elemPt,intPt,3) = press         ! pressure
         globalSdv(elemPt,intPt,4) = gammadotp_tau ! plastic shear strain rate
         globalSdv(elemPt,intPt,5) = v_tau(1)      ! velocity - 1-component
         globalSdv(elemPt,intPt,6) = v_tau(2)      ! velocity - 2-component
         globalSdv(elemPt,intPt,7) = v_tau(3)      ! velocity - 3-component
         globalSdv(elemPt,intPt,8) = v_mag         ! velocity - magnitude


         ! Save the state variables at this integ point
         !  at the end of the increment.
         !
         svars(1+jj)  = Fp_tau(1,1)
         svars(2+jj)  = Fp_tau(2,1)
         svars(3+jj)  = Fp_tau(3,1)
         svars(4+jj)  = Fp_tau(1,2)
         svars(5+jj)  = Fp_tau(2,2)
         svars(6+jj)  = Fp_tau(3,2)
         svars(7+jj)  = Fp_tau(1,3)
         svars(8+jj)  = Fp_tau(2,3)
         svars(9+jj)  = Fp_tau(3,3)
         jj = jj + nlSdv ! setup for the next intPt


         ! Compute/update the displacement residual vector.
         !
         ! First, compile the B-matrix, a matrix of shape function 
         !  derivatives in the current configuration which takes 
         !  advantage of the symmetry of the Cauchy stress.
         !
         Bmat = zero
         do k=1,nNode
            Bmat(1,1+nDim*(k-1)) = dshC(k,1)
            Bmat(2,2+nDim*(k-1)) = dshC(k,2)
            Bmat(3,3+nDim*(k-1)) = dshC(k,3)
            Bmat(4,1+nDim*(k-1)) = dshC(k,2)
            Bmat(4,2+nDim*(k-1)) = dshC(k,1)
            Bmat(5,2+nDim*(k-1)) = dshC(k,3)
            Bmat(5,3+nDim*(k-1)) = dshC(k,2)
            Bmat(6,1+nDim*(k-1)) = dshC(k,3)
            Bmat(6,3+nDim*(k-1)) = dshC(k,1)
         enddo
         !
         ! Next, compile the matrix of Cauchy stress components,
         !  called the S-matrix.
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(3,3)
         Smat(4,1) = T_tau(1,2)
         Smat(5,1) = T_tau(2,3)
         Smat(6,1) = T_tau(1,3)
         !
         ! Note: It would be straightforward to include body forces
         !       at this point; however, it is expected that body forces 
         !       will be applied to the dummy mesh in the input file.
         !
         ! Finally, add the contribution of this integration point
         !  to the displacement residual vector.
         !
         Ru = Ru - detMapJC*w(intpt)*matmul(transpose(Bmat),Smat)


         ! Compute/update the fluidity residual vector
         !
         ! Compile a G-matrix which is appropriately dimensioned for
         !  the fluidity residual.
         !
         Gmatg = zero
         do i=1,nDim
            do k=1,nNode
               Gmatg(i,k) = dshC(k,i)
            end do
         end do
         !
         ! Compile the N-matrix, a matrix of the shape functions 
         !  evaluated at this integration point
         !
         Nmatg = zero
         do i=1,nNode
            Nmatg(1,i) = sh(i)
         end do
         !
         ! Finally, add the contribution of this integration point
         !  to the fluidity residual vector.
         !
         Rg = Rg + detMapJC*w(intpt)*
     +     (matmul(transpose(Gmatg),dgdx)+matmul(transpose(Nmatg),gr))


         ! Compute/update the displacement tangent matrix
         !
         ! Compile the matrices of shape function derivatives 
         !  in the current configuration both at this integration point 
         !  as well as at the element centroid. Here, while we assume 
         !  that minor symmetries of the material tangents are present, 
         !  we do not assume that we can exploit major symmetries. We 
         !  refer to this matrix as the G-matrix.
         !
         Gmat = zero
         do k=1,nNode
            Gmat(1,1+nDim*(k-1)) = dshC(k,1)
            Gmat(2,2+nDim*(k-1)) = dshC(k,1)
            Gmat(3,3+nDim*(k-1)) = dshC(k,1)
            Gmat(4,1+nDim*(k-1)) = dshC(k,2)
            Gmat(5,2+nDim*(k-1)) = dshC(k,2)
            Gmat(6,3+nDim*(k-1)) = dshC(k,2)
            Gmat(7,1+nDim*(k-1)) = dshC(k,3)
            Gmat(8,2+nDim*(k-1)) = dshC(k,3)
            Gmat(9,3+nDim*(k-1)) = dshC(k,3)
         enddo
         !
         Gmat0 = zero
         do k=1,nNode
            Gmat0(1,1+nDim*(k-1)) = dshC0(k,1)
            Gmat0(2,2+nDim*(k-1)) = dshC0(k,1)
            Gmat0(3,3+nDim*(k-1)) = dshC0(k,1)
            Gmat0(4,1+nDim*(k-1)) = dshC0(k,2)
            Gmat0(5,2+nDim*(k-1)) = dshC0(k,2)
            Gmat0(6,3+nDim*(k-1)) = dshC0(k,2)
            Gmat0(7,1+nDim*(k-1)) = dshC0(k,3)
            Gmat0(8,2+nDim*(k-1)) = dshC0(k,3)
            Gmat0(9,3+nDim*(k-1)) = dshC0(k,3)
         enddo
         !
         ! Compile the matrix of spatial mechanical tangent components. 
         !  Call it the C-matrix.
         !
         Cmat = zero
         Cmat(1,1) = Ctang(1,1,1,1)
         Cmat(1,2) = Ctang(1,1,2,1)
         Cmat(1,3) = Ctang(1,1,3,1)
         Cmat(1,4) = Ctang(1,1,1,2)
         Cmat(1,5) = Ctang(1,1,2,2)
         Cmat(1,6) = Ctang(1,1,3,2)
         Cmat(1,7) = Ctang(1,1,1,3)
         Cmat(1,8) = Ctang(1,1,2,3)
         Cmat(1,9) = Ctang(1,1,3,3)
         Cmat(2,1) = Ctang(2,1,1,1)
         Cmat(2,2) = Ctang(2,1,2,1)
         Cmat(2,3) = Ctang(2,1,3,1)
         Cmat(2,4) = Ctang(2,1,1,2)
         Cmat(2,5) = Ctang(2,1,2,2)
         Cmat(2,6) = Ctang(2,1,3,2)
         Cmat(2,7) = Ctang(2,1,1,3)
         Cmat(2,8) = Ctang(2,1,2,3)
         Cmat(2,9) = Ctang(2,1,3,3)
         Cmat(3,1) = Ctang(3,1,1,1)
         Cmat(3,2) = Ctang(3,1,2,1)
         Cmat(3,3) = Ctang(3,1,3,1)
         Cmat(3,4) = Ctang(3,1,1,2)
         Cmat(3,5) = Ctang(3,1,2,2)
         Cmat(3,6) = Ctang(3,1,3,2)
         Cmat(3,7) = Ctang(3,1,1,3)
         Cmat(3,8) = Ctang(3,1,2,3)
         Cmat(3,9) = Ctang(3,1,3,3)
         Cmat(4,1) = Ctang(1,2,1,1)
         Cmat(4,2) = Ctang(1,2,2,1)
         Cmat(4,3) = Ctang(1,2,3,1)
         Cmat(4,4) = Ctang(1,2,1,2)
         Cmat(4,5) = Ctang(1,2,2,2)
         Cmat(4,6) = Ctang(1,2,3,2)
         Cmat(4,7) = Ctang(1,2,1,3)
         Cmat(4,8) = Ctang(1,2,2,3)
         Cmat(4,9) = Ctang(1,2,3,3)
         Cmat(5,1) = Ctang(2,2,1,1)
         Cmat(5,2) = Ctang(2,2,2,1)
         Cmat(5,3) = Ctang(2,2,3,1)
         Cmat(5,4) = Ctang(2,2,1,2)
         Cmat(5,5) = Ctang(2,2,2,2)
         Cmat(5,6) = Ctang(2,2,3,2)
         Cmat(5,7) = Ctang(2,2,1,3)
         Cmat(5,8) = Ctang(2,2,2,3)
         Cmat(5,9) = Ctang(2,2,3,3)
         Cmat(6,1) = Ctang(3,2,1,1)
         Cmat(6,2) = Ctang(3,2,2,1)
         Cmat(6,3) = Ctang(3,2,3,1)
         Cmat(6,4) = Ctang(3,2,1,2)
         Cmat(6,5) = Ctang(3,2,2,2)
         Cmat(6,6) = Ctang(3,2,3,2)
         Cmat(6,7) = Ctang(3,2,1,3)
         Cmat(6,8) = Ctang(3,2,2,3)
         Cmat(6,9) = Ctang(3,2,3,3)
         Cmat(7,1) = Ctang(1,3,1,1)
         Cmat(7,2) = Ctang(1,3,2,1)
         Cmat(7,3) = Ctang(1,3,3,1)
         Cmat(7,4) = Ctang(1,3,1,2)
         Cmat(7,5) = Ctang(1,3,2,2)
         Cmat(7,6) = Ctang(1,3,3,2)
         Cmat(7,7) = Ctang(1,3,1,3)
         Cmat(7,8) = Ctang(1,3,2,3)
         Cmat(7,9) = Ctang(1,3,3,3)
         Cmat(8,1) = Ctang(2,3,1,1)
         Cmat(8,2) = Ctang(2,3,2,1)
         Cmat(8,3) = Ctang(2,3,3,1)
         Cmat(8,4) = Ctang(2,3,1,2)
         Cmat(8,5) = Ctang(2,3,2,2)
         Cmat(8,6) = Ctang(2,3,3,2)
         Cmat(8,7) = Ctang(2,3,1,3)
         Cmat(8,8) = Ctang(2,3,2,3)
         Cmat(8,9) = Ctang(2,3,3,3)
         Cmat(9,1) = Ctang(3,3,1,1)
         Cmat(9,2) = Ctang(3,3,2,1)
         Cmat(9,3) = Ctang(3,3,3,1)
         Cmat(9,4) = Ctang(3,3,1,2)
         Cmat(9,5) = Ctang(3,3,2,2)
         Cmat(9,6) = Ctang(3,3,3,2)
         Cmat(9,7) = Ctang(3,3,1,3)
         Cmat(9,8) = Ctang(3,3,2,3)
         Cmat(9,9) = Ctang(3,3,3,3)
         !
         ! Compile the Q-matrix for use in the Fbar method.
         !  (See Eq. (15) of de Souza Neto et al., IJSS (1996))
         !
         Qmat = zero
         Qmat(1,1) = third*(Cmat(1,1)+Cmat(1,5)+Cmat(1,9)) 
     +        - (two/three)*T_tau(1,1)
         Qmat(2,1) = third*(Cmat(2,1)+Cmat(2,5)+Cmat(2,9))
     +        - (two/three)*T_tau(2,1)
         Qmat(3,1) = third*(Cmat(3,1)+Cmat(3,5)+Cmat(3,9))
     +        - (two/three)*T_tau(3,1)
         Qmat(4,1) = third*(Cmat(4,1)+Cmat(4,5)+Cmat(4,9))
     +        - (two/three)*T_tau(1,2)
         Qmat(5,1) = third*(Cmat(5,1)+Cmat(5,5)+Cmat(5,9))
     +        - (two/three)*T_tau(2,2)
         Qmat(6,1) = third*(Cmat(6,1)+Cmat(6,5)+Cmat(6,9))
     +        - (two/three)*T_tau(3,2)
         Qmat(7,1) = third*(Cmat(7,1)+Cmat(7,5)+Cmat(7,9))
     +        - (two/three)*T_tau(1,3)
         Qmat(8,1) = third*(Cmat(8,1)+Cmat(8,5)+Cmat(8,9))
     +        - (two/three)*T_tau(2,3)
         Qmat(9,1) = third*(Cmat(9,1)+Cmat(9,5)+Cmat(9,9))
     +        - (two/three)*T_tau(3,3)
         Qmat(1,5) = Qmat(1,1)
         Qmat(2,5) = Qmat(2,1)
         Qmat(3,5) = Qmat(3,1)
         Qmat(4,5) = Qmat(4,1)
         Qmat(5,5) = Qmat(5,1)
         Qmat(6,5) = Qmat(6,1)
         Qmat(7,5) = Qmat(7,1)
         Qmat(8,5) = Qmat(8,1)
         Qmat(9,5) = Qmat(9,1)
         Qmat(1,9) = Qmat(1,1)
         Qmat(2,9) = Qmat(2,1)
         Qmat(3,9) = Qmat(3,1)
         Qmat(4,9) = Qmat(4,1)
         Qmat(5,9) = Qmat(5,1)
         Qmat(6,9) = Qmat(6,1)
         Qmat(7,9) = Qmat(7,1)
         Qmat(8,9) = Qmat(8,1)
         Qmat(9,9) = Qmat(9,1)
         !
         ! Finally, add the contribution of this integration point
         !  to the tangent matrix.
         !
         Kuu = Kuu + detMapJC*w(intpt)*
     +              (
     +              matmul(matmul(transpose(Gmat),Cmat),Gmat)
     +              + matmul(transpose(Gmat),matmul(Qmat,(Gmat0-Gmat)))
     +              )


         ! Compute/update the displacement/fluidity tangent matrix
         !
         ! Compile the matrix of spatial displacement/fluidity tangent 
         !  components. Call it the Lambda-matrix.
         !
         LambdaMat = zero
         LambdaMat(1,1) = Lambda(1,1)
         LambdaMat(2,1) = Lambda(2,1)
         LambdaMat(3,1) = Lambda(3,1)
         LambdaMat(4,1) = Lambda(1,2)
         LambdaMat(5,1) = Lambda(2,2)
         LambdaMat(6,1) = Lambda(3,2)
         LambdaMat(7,1) = Lambda(1,3)
         LambdaMat(8,1) = Lambda(2,3)
         LambdaMat(9,1) = Lambda(3,3)
         !
         ! Add the contribution of this integration point to the 
         !  tangent matrix.
         !
         Kug = Kug + detMapJC*w(intPt)*
     +                (
     +                matmul(transpose(Gmat),matmul(LambdaMat,Nmatg))
     +                )
         
         
         ! Compute/update the fluidity/displacement tangent matrix
         !
         ! Compile the matrix of spatial fluidity/displacement tangent 
         !  components. Call it the Gamma-matrix.
         !
         GammaMat = zero
         GammaMat(1,1) = Gamma(1,1)
         GammaMat(1,2) = Gamma(2,1)
         GammaMat(1,3) = Gamma(3,1)
         GammaMat(1,4) = Gamma(1,2)
         GammaMat(1,5) = Gamma(2,2)
         GammaMat(1,6) = Gamma(3,2)
         GammaMat(1,7) = Gamma(1,3)
         GammaMat(1,8) = Gamma(2,3)
         GammaMat(1,9) = Gamma(3,3)
         !
         Qmatg = zero
         Qmatg(1,1) = third*(Gamma(1,1)+Gamma(2,2)+Gamma(3,3)) - gr(1,1)
         Qmatg(1,5) = third*(Gamma(1,1)+Gamma(2,2)+Gamma(3,3)) - gr(1,1)
         Qmatg(1,9) = third*(Gamma(1,1)+Gamma(2,2)+Gamma(3,3)) - gr(1,1)
         !
         ! Compile another matrix of spatial fluidity/displacement tangent 
         !  components. Call it the A-matrix.
         !
         Atang = zero
         do j=1,3
            do k=1,3
               do l=1,3
                  Atang(j,k,l) = Atang(j,k,l) - 
     +                  (dgdx(k,1)*iden(j,l) + 
     +                   dgdx(l,1)*iden(j,k) -
     +                   dgdx(j,1)*iden(k,l))
               end do
            end do
         end do
         !
         Amat = zero
         Amat(1,1) = Atang(1,1,1)
         Amat(1,2) = Atang(1,2,1)
         Amat(1,3) = Atang(1,3,1)
         Amat(1,4) = Atang(1,1,2)
         Amat(1,5) = Atang(1,2,2)
         Amat(1,6) = Atang(1,3,2)
         Amat(1,7) = Atang(1,1,3)
         Amat(1,8) = Atang(1,2,3)
         Amat(1,9) = Atang(1,3,3)
         Amat(2,1) = Atang(2,1,1)
         Amat(2,2) = Atang(2,2,1)
         Amat(2,3) = Atang(2,3,1)
         Amat(2,4) = Atang(2,1,2)
         Amat(2,5) = Atang(2,2,2)
         Amat(2,6) = Atang(2,3,2)
         Amat(2,7) = Atang(2,1,3)
         Amat(2,8) = Atang(2,2,3)
         Amat(2,9) = Atang(2,3,3)
         Amat(3,1) = Atang(3,1,1)
         Amat(3,2) = Atang(3,2,1)
         Amat(3,3) = Atang(3,3,1)
         Amat(3,4) = Atang(3,1,2)
         Amat(3,5) = Atang(3,2,2)
         Amat(3,6) = Atang(3,3,2)
         Amat(3,7) = Atang(3,1,3)
         Amat(3,8) = Atang(3,2,3)
         Amat(3,9) = Atang(3,3,3)
         !
         ! Add the contribution of this integration point
         !  to the tangent matrix.
         !
         Kgu = Kgu - detMapJC*w(intPt)*
     +              (
     +              matmul(matmul(transpose(Gmatg),Amat),Gmat) + 
     +              matmul(matmul(transpose(Nmatg),GammaMat),Gmat) + 
     +              matmul(matmul(transpose(Nmatg),Qmatg),(Gmat0-Gmat))
     +              )


         ! Compute/update the fluidity tangent matrix
         !
         Kgg = Kgg - detMapJC*w(intPt)*
     +               (
     +               matmul(transpose(Gmatg),Gmatg) + 
     +               matmul(transpose(Nmatg),matmul(dgrdg,Nmatg))
     +               )


      enddo
      !
      ! End the loop over integration points
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! Traction terms and surface fluidity terms are not implemented.
      ! Mechanical, traction- and pressure-type boundary conditions 
      !  may be applied to the dummy mesh using the Abaqus built-in 
      !  commands *Dload or *Dsload.
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.  This
      !  is essentially giving Abaqus the residual and the tangent matrix.
      !
      ! Return Abaqus the right hand side vector
      !
      do i=1,nNode
         A11 = (nDim+1)*(i-1)+1
         A12 = nDim*(i-1)+1
         !
         ! displacement
         !
         rhs(A11,1) = Ru(A12,1)
         rhs(A11+1,1) = Ru(A12+1,1)
         rhs(A11+2,1) = Ru(A12+2,1)
         !
         ! granular fluidity
         !
         rhs(A11+3,1) = Rg(i,1)
      enddo
      !
      ! Return Abaqus the tangent matrix
      !
      amatrx = zero
      do i=1,nNode
         do j=1,nNode
            A11 = (nDim+1)*(i-1)+1
            A12 = nDim*(i-1)+1
            B11 = (nDim+1)*(j-1)+1
            B12 = nDim*(j-1)+1
            !
            ! displacement
            !
            amatrx(A11,B11)     = Kuu(A12,B12)
            amatrx(A11,B11+1)   = Kuu(A12,B12+1)
            amatrx(A11,B11+2)   = Kuu(A12,B12+2)
            amatrx(A11+1,B11)   = Kuu(A12+1,B12)
            amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
            amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
            amatrx(A11+2,B11)   = Kuu(A12+2,B12)
            amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
            amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
            !
            ! displacement/granular fluidity
            !
            amatrx(A11,B11+3) = Kug(A12,j)
            amatrx(A11+1,B11+3) = Kug(A12+1,j)
            amatrx(A11+2,B11+3) = Kug(A12+2,j)
            !
            ! granular fluidity/displacement
            !
            amatrx(A11+3,B11) = Kgu(i,B12)
            amatrx(A11+3,B11+1) = Kgu(i,B12+1)
            amatrx(A11+3,B11+2) = Kgu(i,B12+2)
            !
            ! granular fluidity
            !
            amatrx(A11+3,B11+3) = Kgg(i,j)
            !
         enddo
      enddo
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine uel

!************************************************************************
!     Material subroutines
!************************************************************************

      subroutine nonlocal(props,nprops,dtime,F_tau,Fp_t,g_tau,
     +                 T_tau,gr,Fp_tau,mu_tau,taubar_tau,press,
     +                 gammadotp_tau,Ctang,Lambda,Gamma,dgrdg,stat)

      implicit none
      !
      integer i,j,k,l,m,n,nprops,stat
      !
      real*8 props(nprops),dtime,F_tau(3,3),Fp_t(3,3),g_tau,T_tau(3,3),
     +  gr(1,1),Fp_tau(3,3),mu_tau,press,gammadotp_tau,Ctang(3,3,3,3),
     +  Lambda(3,3),Gamma(3,3),dgrdg(1,1),Iden(3,3),Gshear,Kbulk,
     +  rhos,d,mus,b,A,pMin,Fp_t_inv(3,3),det_Fp_t,Fe_tr(3,3),
     +  Be_tr(3,3),Re_tr(3,3),Ue_tr(3,3),Ee_tr(3,3),trEe_tr,Ee0_tr(3,3),
     +  Me_tr(3,3),Me0_tr(3,3),taubar_tr,Np_tr(3,3),Npbar_tr(3,3),
     +  taubar_tau,Me_tau(3,3),det_F_tau,dtaudtau_tr,dtaudpress,
     +  dtaudE(3,3),dtaudg,dmudtau_tr,dmudpress,dmudE(3,3),dmudg,
     +  Dp_tau(3,3),Dp_eig(3),Dp_vec(3,3),expdtDp(3,3),det_Fp_tau,g_loc,
     +  dg_locdmu,dg_locdp,dg_locdE(3,3),dg_locdg,xicap,delmu,acoeff,xi,
     +  dxidmu,dxidE(3,3),dxidg,Dtang(3,3,3,3),Gtilde,fac,
     +  Ltang(3,3,3,3),Btang(3,3,3,3),tempm4(3,3,3,3),dgrdE(3,3)
      !
      real*8 zero,one,two,three,fourth,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0)


      ! Status
      !
      stat = 1
 

      ! Identity matrix
      !
      call onem(Iden)
 

      ! Obtain material properties
      !
      Gshear = props(1)    ! Shear modulus (Pa)
      Kbulk  = props(2)    ! Bulk modulus (Pa)
      rhos   = props(3)    ! Grain density (kg/m^3)
      d      = props(4)    ! Grain diameter (m)
      mus    = props(5)    ! Static friction coefficient
      b      = props(6)    ! Rheology parameter
      A      = props(7)    ! Nonlocal amplitude
      pMin   = Kbulk/1.d12 ! Minimum pressure (in terms of bulk modulus) (Pa)


      ! Step 1: Compute the trial elastic deformation gradient
      !
      call matInv3D(Fp_t,Fp_t_inv,det_Fp_t,stat)
      !
      if (stat==0) then
         !
         ! Fp_t is not physical, terminate job
         !
         write(*,*) 'Fp_t is not physical, terminating job'
         call xit
         !
      end if
      !
      Fe_tr = matmul(F_tau,Fp_t_inv)
      Be_tr = matmul(Fe_tr,transpose(Fe_tr)) ! Tr. elastic left C-G tensor


      ! Step 2: Compute the trial kinematics
      !
      call skinem(Fe_tr,Re_tr,Ue_tr,Ee_tr,stat)
      !
      if (stat==0) then
         !
         ! Problem in the kinematics, cut back the time increment
         !
         write(*,*) 'Problem in the kinematics, cutting back'
         return
         !
      end if
      !
      trEe_tr = Ee_tr(1,1) + Ee_tr(2,2) + Ee_tr(3,3) ! Tr. vol. strain
      Ee0_tr = Ee_tr - third*trEe_tr*Iden ! Tr. strain dev.


      ! Step 3: Calculate the trial stress and related quantities
      !
      ! Calculate the trial Mandel stress 
      !
      Me_tr = two*Gshear*Ee0_tr + Kbulk*trEe_tr*Iden
      !
      ! Calculate the trial pressure
      !
      press = -third*(Me_tr(1,1) + Me_tr(2,2) + Me_tr(3,3))
      !
      ! Calculate the trial stress deviator
      !
      Me0_tr = Me_tr + press*Iden
      !
      ! Calculate the trial equivalent shear stress
      !
      taubar_tr = dsqrt(half*sum(Me0_tr*Me0_tr))
      !
      ! Calculate the trial direction of plastic flow 
      !
      if(taubar_tr.le.zero) then
         Np_tr = zero
      else
         Np_tr = dsqrt(half)*(Me0_tr/taubar_tr)
      endif
      !
      ! Calculate the spatial tr. direction of plastic flow
      !
      Npbar_tr = matmul(Re_tr,matmul(Np_tr,transpose(Re_tr)))
      !
      ! Check that the pressure is positive
      !
      if ((press.ge.(-pMin)).and.(press<pMin)) then
         press = pMin
      elseif (press<(-pMin)) then
         !
         ! Pressure is too negative, cut back
         !
         stat = 0
         return
      end if
      
      
      ! Step 4: Update the stresses
      !
      ! Check that the fluidity isn't too negative
      !
      if ((press + Gshear*dtime*g_tau)<zero) then
         !
         ! The fluidity is too negative, cut back
         !
         stat = 0
         return
      end if
      !
      ! Calculate the equivalent shear stress and stress ratio
      !
      taubar_tau = taubar_tr*press/(press + Gshear*dtime*g_tau)
      mu_tau = taubar_tau/press
      !
      ! Calculate the Mandel stress
      !
      Me_tau = Me_tr - dsqrt(two)*(taubar_tr - taubar_tau)*Np_tr
      !
      ! Calculate the Cauchy stress
      !
      call mdet(F_tau,det_F_tau)
      T_tau = matmul(Re_tr,matmul(Me_tau,transpose(Re_tr)))/det_F_tau
      !
      ! Calculate the derivatives that will be necessary
      !  for the tangents
      !
      dtaudtau_tr = press/(press + Gshear*dtime*g_tau)
      dtaudpress = taubar_tr*Gshear*dtime*g_tau/
     +                  ((press + Gshear*dtime*g_tau)**two)
      dtaudE = dsqrt(two)*Gshear*dtaudtau_tr*Npbar_tr - 
     +                  Kbulk*dtaudpress*Iden
      dtaudg = -taubar_tr*press*Gshear*dtime/
     +                  ((press + Gshear*dtime*g_tau)**two)
      !
      dmudtau_tr = one/(press + Gshear*dtime*g_tau)
      dmudpress = -taubar_tr/((press + Gshear*dtime*g_tau)**two)
      dmudE = dsqrt(two)*Gshear*dmudtau_tr*Npbar_tr - 
     +                  Kbulk*dmudpress*Iden
      dmudg = -taubar_tr*Gshear*dtime/
     +                  ((press + Gshear*dtime*g_tau)**two)


      ! Step 5: Update the plastic distortion
      !
      ! Compute the plastic stretching
      !
      gammadotp_tau = g_tau*mu_tau
      Dp_tau = dsqrt(one/two)*gammadotp_tau*Np_tr
      !
      ! Compute the plastic deformation gradient at the
      !  end of the increment using the exponential map
      !
      if(gammadotp_tau.le.zero) then
         Fp_tau = Fp_t
      else
         call spectral(dtime*Dp_tau,Dp_eig,Dp_vec,stat)
         expdtDp = zero
         expdtDp(1,1) = dexp(Dp_eig(1))
         expdtDp(2,2) = dexp(Dp_eig(2))
         expdtDp(3,3) = dexp(Dp_eig(3))
         expdtDp = matmul(matmul(Dp_vec,expdtDp),transpose(Dp_vec))
         Fp_tau = matmul(expdtDp,Fp_t)
      endif
      !
      if (stat==0) then
         !
         write(*,*) 'Problem in updating Fp_tau'
         return
         !
      end if
      !
      ! Calculate det(Fp_tau), and check to make sure that det(Fp_tau)>0
      !
      call mdet(Fp_tau,det_Fp_tau)
      if(det_Fp_tau.le.zero) then
         stat = 0
         write(*,*) 'det(Fp_tau).le.zero in INTEG'
         return
      endif
      
      
      ! Step 6: Calculate the scalar contribution to the fluidity residual
      !
      ! Compute the local fluidity, g_loc, as well as the necessary
      !  derivatives
      !
      if (mu_tau > mus) then
         g_loc = (one/(b*d))*dsqrt(press/rhos)*(one - (mus/mu_tau))
         dg_locdmu = (one/(b*d))*dsqrt(press/rhos)*(mus/(mu_tau**two))
         dg_locdp = half*(one/(b*d))*dsqrt(one/(press*rhos))*
     +                               (one - (mus/mu_tau))
      else
         g_loc = zero
         dg_locdmu = zero
         dg_locdp = zero
      end if
      !
      dg_locdE = dg_locdmu*dmudE - Kbulk*dg_locdp*Iden
      dg_locdg = dg_locdmu*dmudg
      !
      ! Compute the cooperativity length, xi, as well as the necessary
      !  derivatives
      !
      xicap = 1000.d0
      delmu = ((5.d0/4.d0)*A/xicap)**two
      acoeff = 0.25d0*A/(((5.d0/4.d0)*A/xicap)**5.d0)
      if (mu_tau > (mus + delmu)) then
         xi = A*d/dsqrt(mu_tau - mus)
         dxidmu = -half*A*d*((mu_tau - mus)**(-3.d0/2.d0))
      elseif (mu_tau < (mus - delmu)) then
         xi = A*d/dsqrt(mus - mu_tau)
         dxidmu =  half*A*d*((mus - mu_tau)**(-3.d0/2.d0))
      else
         xi = xicap*d - acoeff*d*((mu_tau - mus)**two)
         dxidmu = - two*acoeff*d*(mu_tau - mus)
      end if
      !
      dxidE = dxidmu*dmudE
      dxidg = dxidmu*dmudg
      !
      ! Compute the fluidity residual
      !
      gr = (g_tau - g_loc)/(xi**two)


      ! Calculate the mechanical tangent modulus
      !
      ! First calculate the constitutive contribution.
      !
      Dtang = zero
      !
      if (taubar_tr.gt.zero) then
         Gtilde = Gshear*(taubar_tau/taubar_tr)
         fac = ((taubar_tau/taubar_tr) - dtaudtau_tr)
      else
         Gtilde = Gshear
         fac = zero
      end if
      !
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  Dtang(i,j,k,l) = Dtang(i,j,k,l)
     +                    + Kbulk*Iden(i,j)*Iden(k,l)
     +                    + Gtilde*
     +                       (Iden(i,k)*Iden(j,l) + 
     +                        Iden(i,l)*Iden(j,k) - 
     +                        (two/three)*Iden(i,j)*Iden(k,l))
     +                    - two*Gshear*fac*Npbar_tr(i,j)*Npbar_tr(k,l)
     +                    - dsqrt(two)*Kbulk*dtaudpress*
     +                       Npbar_tr(i,j)*Iden(k,l)
               enddo
            enddo
         enddo
      enddo
      !
      Dtang = Dtang/det_F_tau
      !
      ! Get the fourth order tensor L = d(ln(Be_tr))/d(Be_tr)
      !
      call dlnxdx(Be_tr,Ltang)
      !
      ! Get the fourth order tensor B
      !
      Btang = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  Btang(i,j,k,l) = Btang(i,j,k,l)
     +                + Iden(i,k)*Be_tr(j,l) + Iden(j,k)*Be_tr(i,l)
               enddo
            enddo
         enddo
      enddo
      !
      ! Finally, construct the mechanical tangent
      !
      call mprod4(Ltang,half*Btang,tempm4)
      call mprod4(Dtang,tempm4,Ctang)
      !
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  Ctang(i,j,k,l) = Ctang(i,j,k,l)
     +                - T_tau(i,l)*Iden(j,k)
               enddo
            enddo
         enddo
      enddo


      ! Calculate the stress/fluidity tangent
      !
      Lambda = dsqrt(two)*dtaudg*Npbar_tr/det_F_tau


      ! Calculate the fluidity residual/displacement tangent
      !
      ! First, calculate the constitutive contribution
      !
      dgrdE = -dg_locdE/(xi**two) - 
     +           two*((g_tau - g_loc)/(xi**three))*dxidE
      !
      ! Assemble the tangent
      !
      Gamma = zero
      do m=1,3
         do n=1,3
            do k=1,3
               do l=1,3
                  Gamma(k,l) = Gamma(k,l)
     +                + dgrdE(m,n)*tempm4(m,n,k,l)
               enddo
            enddo
         enddo
      enddo
      !
      Gamma = Gamma + gr(1,1)*Iden
      

      ! Calculate the fluidity residual/fluidity tangent
      !
      dgrdg = (one - dg_locdg)/(xi**two) - 
     +        two*((g_tau - g_loc)/(xi**three))*dxidg


      return
      end subroutine nonlocal
      
!****************************************************************************

      subroutine local(props,nprops,dtime,F_tau,Fp_t,g_tau,
     +                 T_tau,gr,Fp_tau,mu_tau,taubar_tau,press,
     +                 gammadotp_tau,Ctang,Lambda,Gamma,dgrdg,stat)

      implicit none
      !
      integer i,j,k,l,m,n,nprops,stat
      !
      real*8 props(nprops),dtime,F_tau(3,3),Fp_t(3,3),g_tau,T_tau(3,3),
     +  gr(1,1),Fp_tau(3,3),mu_tau,press,gammadotp_tau,Ctang(3,3,3,3),
     +  Lambda(3,3),Gamma(3,3),dgrdg(1,1),Iden(3,3),Gshear,Kbulk,
     +  rhos,d,mus,b,A,pMin,Fp_t_inv(3,3),det_Fp_t,Fe_tr(3,3),
     +  Be_tr(3,3),Re_tr(3,3),Ue_tr(3,3),Ee_tr(3,3),trEe_tr,Ee0_tr(3,3),
     +  Me_tr(3,3),Me0_tr(3,3),taubar_tr,Np_tr(3,3),Npbar_tr(3,3),
     +  taubar_tau,Me_tau(3,3),det_F_tau,dtaudtau_tr,dtaudpress,
     +  dtaudE(3,3),dmudtau_tr,dmudpress,dmudE(3,3),Dp_tau(3,3),
     +  Dp_eig(3),Dp_vec(3,3),expdtDp(3,3),det_Fp_tau,g_loc,dg_locdmu,
     +  dg_locdp,dg_locdE(3,3),xicap,delmu,acoeff,xi,dxidmu,dxidE(3,3),
     +  Dtang(3,3,3,3),Gtilde,fac,Ltang(3,3,3,3),Btang(3,3,3,3),
     +  tempm4(3,3,3,3),dgrdE(3,3)
      !
      real*8 zero,one,two,three,fourth,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,fourth=1.d0/4.d0,
     +     third=1.d0/3.d0,half=1.d0/2.d0)


      ! Status
      !
      stat = 1
 

      ! Identity matrix
      !
      call onem(Iden)
 

      ! Obtain material properties
      !
      Gshear = props(1)    ! Shear modulus (Pa)
      Kbulk  = props(2)    ! Bulk modulus (Pa)
      rhos   = props(3)    ! Grain density (kg/m^3)
      d      = props(4)    ! Grain diameter (m)
      mus    = props(5)    ! Static friction coefficient
      b      = props(6)    ! Rheology parameter
      A      = props(7)    ! Nonlocal amplitude
      pMin   = Kbulk/1.d12 ! Minimum pressure (in terms of bulk modulus) (Pa)


      ! Step 1: Compute the trial elastic deformation gradient
      !
      call matInv3D(Fp_t,Fp_t_inv,det_Fp_t,stat)
      !
      if (stat==0) then
         !
         ! Fp_t is not physical, terminate job
         !
         write(*,*) 'Fp_t is not physical, terminating job'
         call xit
         !
      end if
      !
      Fe_tr = matmul(F_tau,Fp_t_inv)
      Be_tr = matmul(Fe_tr,transpose(Fe_tr)) ! Tr. elastic left C-G tensor


      ! Step 2: Compute the trial kinematics
      !
      call skinem(Fe_tr,Re_tr,Ue_tr,Ee_tr,stat)
      !
      if (stat==0) then
         !
         ! Problem in the kinematics, cut back the time increment
         !
         write(*,*) 'Problem in the kinematics, cutting back'
         return
         !
      end if
      !
      trEe_tr = Ee_tr(1,1) + Ee_tr(2,2) + Ee_tr(3,3) ! Tr. vol. strain
      Ee0_tr = Ee_tr - third*trEe_tr*Iden ! Tr. strain dev.


      ! Step 3: Calculate the trial stress and related quantities
      !
      ! Calculate the trial Mandel stress 
      !
      Me_tr = two*Gshear*Ee0_tr + Kbulk*trEe_tr*Iden
      !
      ! Calculate the trial pressure
      !
      press = -third*(Me_tr(1,1) + Me_tr(2,2) + Me_tr(3,3))
      !
      ! Calculate the trial stress deviator
      !
      Me0_tr = Me_tr + press*Iden
      !
      ! Calculate the trial equivalent shear stress
      !
      taubar_tr = dsqrt(half*sum(Me0_tr*Me0_tr))
      !
      ! Calculate the trial direction of plastic flow 
      !
      if(taubar_tr.le.zero) then
         Np_tr = zero
      else
         Np_tr = dsqrt(half)*(Me0_tr/taubar_tr)
      endif
      !
      ! Calculate the spatial tr. direction of plastic flow
      !
      Npbar_tr = matmul(Re_tr,matmul(Np_tr,transpose(Re_tr)))
      !
      ! Check that the pressure is positive
      !
      if ((press.ge.(-pMin)).and.(press<pMin)) then
         press = pMin
      elseif (press<(-pMin)) then
         stat = 0
         write(*,*) 'Pressure is negative',press
         return
      end if
      
      
      ! Step 4: Update the stresses
      !
      ! Calculate the equivalent shear stress and necessary derivatives
      !
      if(taubar_tr.gt.(mus*press)) then
         taubar_tau = (b*d*dsqrt(rhos*press)*taubar_tr + 
     +              mus*press*Gshear*dtime)/
     +              (b*d*dsqrt(rhos*press) + Gshear*dtime)
         dtaudtau_tr = (b*d*dsqrt(rhos*press))/
     +              (b*d*dsqrt(rhos*press) + Gshear*dtime)
         dtaudpress = (half*b*d*dsqrt(rhos/press)*taubar_tr*
     +                 Gshear*dtime + 
     +                 half*b*d*dsqrt(rhos*press)*mus*Gshear*dtime + 
     +                 mus*((Gshear*dtime)**two))/
     +                ((b*d*dsqrt(rhos*press) + Gshear*dtime)**two)
      else
         taubar_tau = taubar_tr
         dtaudtau_tr = one
         dtaudpress = zero
      end if
      !
      dtaudE = dsqrt(two)*Gshear*dtaudtau_tr*Npbar_tr - 
     +                  Kbulk*dtaudpress*Iden
      !
      ! Calculate the stress ratio and necessary derivatives
      !
      mu_tau = taubar_tau/press
      !
      dmudtau_tr = dtaudtau_tr/press
      dmudpress = dtaudpress/press - taubar_tau/(press**two)
      dmudE = dsqrt(two)*Gshear*dmudtau_tr*Npbar_tr - 
     +                  Kbulk*dmudpress*Iden
      !
      ! Calculate the Mandel stress
      !
      Me_tau = Me_tr - dsqrt(two)*(taubar_tr - taubar_tau)*Np_tr
      !
      ! Calculate the Cauchy stress
      !
      call mdet(F_tau,det_F_tau)
      T_tau = matmul(Re_tr,matmul(Me_tau,transpose(Re_tr)))/det_F_tau


      ! Step 5: Update the plastic distortion
      !
      ! Compute the plastic stretching
      !
      if(taubar_tau.gt.(mus*press)) then
         gammadotp_tau = (taubar_tau-(mus*press))/
     +                    (b*d*dsqrt(rhos*press))
      else
         gammadotp_tau = zero
      end if
      !
      Dp_tau = dsqrt(one/two)*gammadotp_tau*Np_tr
      !
      ! Compute the plastic deformation gradient at the
      !  end of the increment using the exponential map
      !
      if(gammadotp_tau.le.zero) then
         Fp_tau = Fp_t
      else
         call spectral(dtime*Dp_tau,Dp_eig,Dp_vec,stat)
         expdtDp = zero
         expdtDp(1,1) = dexp(Dp_eig(1))
         expdtDp(2,2) = dexp(Dp_eig(2))
         expdtDp(3,3) = dexp(Dp_eig(3))
         expdtDp = matmul(matmul(Dp_vec,expdtDp),transpose(Dp_vec))
         Fp_tau = matmul(expdtDp,Fp_t)
      endif
      !
      if (stat==0) then
         !
         write(*,*) 'Problem in updating Fp_tau'
         return
         !
      end if
      !
      ! Calculate det(Fp_tau), and check to make sure that det(Fp_tau)>0
      !
      call mdet(Fp_tau,det_Fp_tau)
      if(det_Fp_tau.le.zero) then
         stat = 0
         write(*,*) 'det(Fp_tau).le.zero in INTEG'
         return
      endif
      
      
      ! Step 6: Calculate the scalar contribution to the fluidity residual
      !
      ! Compute the local fluidity, g_loc, as well as the necessary
      !  derivatives
      !
      if (mu_tau > mus) then
         g_loc = (one/(b*d))*dsqrt(press/rhos)*(one - (mus/mu_tau))
         dg_locdmu = (one/(b*d))*dsqrt(press/rhos)*(mus/(mu_tau**two))
         dg_locdp = half*(one/(b*d))*(one/dsqrt(rhos*press))*
     +                 (one - (mus/mu_tau))
      else
         g_loc = zero
         dg_locdmu = zero
         dg_locdp = zero
      end if
      !
      dg_locdE = dg_locdmu*dmudE - Kbulk*dg_locdp*Iden
      !
      ! Compute the cooperativity length, xi, as well as the necessary
      !  derivatives
      !
      xicap = 1000.d0
      delmu = ((5.d0/4.d0)*A/xicap)**two
      acoeff = 0.25d0*A/(((5.d0/4.d0)*A/xicap)**5.d0)
      if (mu_tau > (mus + delmu)) then
         xi = A*d/dsqrt(mu_tau - mus)
         dxidmu = -half*A*d*((mu_tau - mus)**(-3.d0/2.d0))
      else if (mu_tau < (mus - delmu)) then
         xi = A*d/dsqrt(mus - mu_tau)
         dxidmu =  half*A*d*((mus - mu_tau)**(-3.d0/2.d0))
      else
         xi = xicap*d - acoeff*d*((mu_tau - mus)**two)
         dxidmu = - two*acoeff*d*(mu_tau - mus)
      end if
      !
      dxidE = dxidmu*dmudE
      !
      ! Compute the fluidity residual
      !
      gr = (g_tau - g_loc)/(xi**two)


      ! Calculate the mechanical tangent modulus
      !
      ! First calculate the constitutive contribution.
      !  Since we only calculate the deviatoric part of the
      !  stress in this subroutine, the purely volumetric 
      !  tangent does not appear here.
      !
      Dtang = zero
      !
      if (taubar_tr.gt.zero) then
         Gtilde = Gshear*(taubar_tau/taubar_tr)
         fac = ((taubar_tau/taubar_tr) - dtaudtau_tr)
      else
         Gtilde = Gshear
         fac = zero
      end if
      !
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  Dtang(i,j,k,l) = Dtang(i,j,k,l)
     +                    + Kbulk*Iden(i,j)*Iden(k,l)
     +                    + Gtilde*
     +                       (Iden(i,k)*Iden(j,l) + 
     +                        Iden(i,l)*Iden(j,k) - 
     +                        (two/three)*Iden(i,j)*Iden(k,l))
     +                    - two*Gshear*fac*Npbar_tr(i,j)*Npbar_tr(k,l)
     +                    - dsqrt(two)*Kbulk*dtaudpress*
     +                       Npbar_tr(i,j)*Iden(k,l)
               enddo
            enddo
         enddo
      enddo
      !
      Dtang = Dtang/det_F_tau
      !
      ! Get the fourth order tensor L = d(ln(Be_tr))/d(Be_tr)
      !
      call dlnxdx(Be_tr,Ltang)
      !
      ! Get the fourth order tensor B
      !
      Btang = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  Btang(i,j,k,l) = Btang(i,j,k,l)
     +                + Iden(i,k)*Be_tr(j,l) + Iden(j,k)*Be_tr(i,l)
               enddo
            enddo
         enddo
      enddo
      !
      ! Finally, construct the mechanical tangent
      !
      call mprod4(Ltang,half*Btang,tempm4)
      call mprod4(Dtang,tempm4,Ctang)
      !
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  Ctang(i,j,k,l) = Ctang(i,j,k,l)
     +                - T_tau(i,l)*Iden(j,k)
               enddo
            enddo
         enddo
      enddo


      ! Calculate the stress/fluidity tangent
      !
      Lambda = zero


      ! Calculate the fluidity residual/displacement tangent
      !
      ! First, calculate the constitutive contribution
      !
      dgrdE = -dg_locdE/(xi**two) - 
     +           two*((g_tau - g_loc)/(xi**three))*dxidE
      !
      ! Assemble the tangent
      !
      Gamma = zero
      do m=1,3
         do n=1,3
            do k=1,3
               do l=1,3
                  Gamma(k,l) = Gamma(k,l)
     +                + dgrdE(m,n)*tempm4(m,n,k,l)
               enddo
            enddo
         enddo
      enddo
      !
      Gamma = Gamma + gr(1,1)*Iden
      

      ! Calculate the fluidity residual/fluidity tangent
      !
      dgrdg = one/(xi**two)


      return
      end subroutine local

!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(8,3),w(8)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

!************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt,i,j
      !
      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3),xi,eta,zeta
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      
      
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      
      
      return
      end subroutine calcShape3DLinear

!*************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode),mapJ(3,3),
     +  mapJ_inv(3,3),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      return
      end subroutine mapShape3D
      
!****************************************************************************
!     The next subroutine calculates various kinematical quantities
!      associated with the deformation gradient
!****************************************************************************

      subroutine skinem(F,R,U,E,istat)
      !
      ! This subroutine performs the right polar decomposition
      !  F = RU of the deformation gradient F into a rotation
      !  R and the right stretch tensor U.  The logarithmic 
      !  strain E = ln(U) is also returned.
      !
      !	F(3,3):       the deformation gradient; input
      !	detF:         the determinant of F; detF > 0
      !	R(3,3):       the rotation matrix; output
      !	U(3,3):       the right stretch tensor; output
      !	Uinv(3,3):    the inverse of U
      !	C(3,3):       the right Cauchy-Green tensor
      !	omega(3):     the squares of the principal stretches
      ! Ueigval(3):   the principal stretches
      !	eigvec(3,3):  matrix of eigenvectors of U
      !	E(3,3):       the logarithmic strain tensor; output
      ! istat:        success flag, istat=0 for a failed attempt; output
      !
      implicit none
      !
      integer istat
      !
      real*8 F(3,3),C(3,3),omega(3),Ueigval(3),eigvec(3,3),
     +  U(3,3),E(3,3),Uinv(3,3),R(3,3),detF
     

      !	Store the identity matrix in R, U, and Uinv
      !
      call onem(R)
      call onem(U)
      call onem(Uinv)
      

      ! Store the zero matrix in E
      !
      E = 0.d0
      

      ! Check if the determinant of F is greater than zero.
      !  If not, then print a diagnostic and cut back the 
      !  time increment.
      !
      call mdet(F,detF)
      if (detF.le.0.d0) then
        write(*,'(/5X,A/)') '--problem in kinematics-- the',
     +       ' determinant of F is not greater than 0'
        istat = 0
        return
      end if
      

      ! Calculate the right Cauchy-Green tensor C
      !
      C = matmul(transpose(F),F)
      
 
      ! Calculate the eigenvalues and eigenvectors of C
      !
      call spectral(C,omega,eigvec,istat)
      

      ! Calculate the principal values of U and E
      !
      Ueigval(1) = dsqrt(omega(1))
      Ueigval(2) = dsqrt(omega(2))
      Ueigval(3) = dsqrt(omega(3))
      !
      U(1,1) = Ueigval(1)
      U(2,2) = Ueigval(2)
      U(3,3) = Ueigval(3)
      !
      E(1,1) = dlog(Ueigval(1))
      E(2,2) = dlog(Ueigval(2))
      E(3,3) = dlog(Ueigval(3))
      

      ! Calculate the complete tensors U and E
      !
      U = matmul(matmul(eigvec,U),transpose(eigvec))
      E = matmul(matmul(eigvec,E),transpose(eigvec))
      

      ! Calculate Uinv
      !
      call matInv3D(U,Uinv,detF,istat)
      

      ! calculate R
      !
      R = matmul(F,Uinv)
      

      return
      end subroutine skinem

!****************************************************************************
!     The following subroutines calculate the spectral
!      decomposition of a symmetric 3 by 3 matrix
!****************************************************************************

      subroutine spectral(A,D,V,istat)
      !
      ! This subroutine calculates the eigenvalues and eigenvectors of
      !  a symmetric 3 by 3 matrix A.
      !
      ! The output consists of a vector D containing the three
      !  eigenvalues in ascending order, and a matrix V whose
      !  columns contain the corresponding eigenvectors.
      !
      implicit none
      !
      integer np,nrot,i,j,istat
      parameter(np=3)
      !
      real*8 D(3),V(3,3),A(3,3),E(3,3)


      E = A
      !
      call jacobi(E,3,np,D,V,nrot,istat)
      call eigsrt(D,V,3,np)
	

      return
      end subroutine spectral
	
!****************************************************************************

      subroutine jacobi(A,n,np,D,V,nrot,istat)
      !
      ! Computes all eigenvalues and eigenvectors of a real symmetric
      !  matrix A, which is of size n by n, stored in a physical
      !  np by np array.  On output, elements of A above the diagonal
      !  are destroyed, but the diagonal and sub-diagonal are unchanged
      !  and give full information about the original symmetric matrix.
      !  Vector D returns the eigenvalues of A in its first n elements.
      !  V is a matrix with the same logical and physical dimensions as
      !  A whose columns contain, upon output, the normalized
      !  eigenvectors of A.  nrot returns the number of Jacobi rotation
      !  which were required.
      !
      ! This subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer ip,iq,n,nmax,np,nrot,i,j,istat
      parameter (nmax=100)
      !
      real*8 A(np,np),D(np),V(np,np),B(nmax),Z(nmax),
     +  sm,tresh,G,T,H,theta,S,C,tau
     
      
      ! Initialize V to the identity matrix
      !
      call onem(V)
      
      
      ! Initialize B and D to the diagonal of A, and Z to zero.
      !  The vector Z will accumulate terms of the form T*A_PQ as
      !  in equation (11.1.14)
      !
      do ip = 1,n
	B(ip) = A(ip,ip)
	D(ip) = B(ip)
	Z(ip) = 0.d0
      end do
      
      
      ! Begin iteration
      !
      nrot = 0
      do i=1,50
          !
          ! Sum off-diagonal elements
          !
          sm = 0.d0
          do ip=1,n-1
            do iq=ip+1,n
	      sm = sm + dabs(A(ip,iq))
            end do
          end do
          !
          ! If sm = 0., then return.  This is the normal return,
          !  which relies on quadratic convergence to machine
          !  underflow.
          !
          if (sm.eq.0.d0) return
          !
          ! In the first three sweeps carry out the PQ rotation only if
          !  |A_PQ| > tresh, where tresh is some threshold value,
          !  see equation (11.1.25).  Thereafter tresh = 0.
          !
          if (i.lt.4) then
            tresh = 0.2d0*sm/n**2
          else
            tresh = 0.d0
          end if
          !
          do ip=1,n-1
            do iq=ip+1,n
              G = 100.d0*dabs(A(ip,iq))
              !
              ! After four sweeps, skip the rotation if the 
              !  off-diagonal element is small.
              !
	      if ((i.gt.4).and.(dabs(D(ip))+G.eq.dabs(D(ip)))
     +            .and.(dabs(D(iq))+G.eq.dabs(D(iq)))) then
                A(ip,iq) = 0.d0
              else if (dabs(A(ip,iq)).gt.tresh) then
                H = D(iq) - D(ip)
                if (dabs(H)+G.eq.dabs(H)) then
                  !
                  ! T = 1./(2.*theta), equation (11.1.10)
                  !
	          T =A(ip,iq)/H
	        else
	          theta = 0.5d0*H/A(ip,iq)
	          T =1.d0/(dabs(theta)+dsqrt(1.d0+theta**2.d0))
	          if (theta.lt.0.d0) T = -T
	        end if
	        C = 1.d0/dsqrt(1.d0 + T**2.d0)
	        S = T*C
	        tau = S/(1.d0 + C)
	        H = T*A(ip,iq)
	        Z(ip) = Z(ip) - H
	        Z(iq) = Z(iq) + H
	        D(ip) = D(ip) - H
	        D(iq) = D(iq) + H
	        A(ip,iq) = 0.d0
                !
                ! Case of rotations 1 <= J < P
		!		
	        do j=1,ip-1
	          G = A(j,ip)
	          H = A(j,iq)
	          A(j,ip) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations P < J < Q
                !
	        do j=ip+1,iq-1
	          G = A(ip,j)
	          H = A(j,iq)
	          A(ip,j) = G - S*(H + G*tau)
	          A(j,iq) = H + S*(G - H*tau)
	        end do
                !
                ! Case of rotations Q < J <= N
                !
	        do j=iq+1,n
                  G = A(ip,j)
	          H = A(iq,j)
	          A(ip,j) = G - S*(H + G*tau)
	          A(iq,j) = H + S*(G - H*tau)
	        end do
	        do j = 1,n
	          G = V(j,ip)
	          H = V(j,iq)
	          V(j,ip) = G - S*(H + G*tau)
	          V(j,iq) = H + S*(G - H*tau)
	        end do
	        nrot = nrot + 1
              end if
	    end do
	  end do
          !
          ! Update D with the sum of T*A_PQ, and reinitialize Z
          !
	  do ip=1,n
	    B(ip) = B(ip) + Z(ip)
	    D(ip) = B(ip)
	    Z(ip) = 0.d0
	  end do
	end do


      ! If the algorithm has reached this stage, then there
      !  are too many sweeps.  Print a diagnostic and cut the 
      !  time increment.
      !
      write (*,'(/1X,A/)') '50 iterations in jacobi should never happen'
      istat = 0
      

      return
      end subroutine jacobi
	
!****************************************************************************

      subroutine eigsrt(D,V,n,np)
      !
      ! Given the eigenvalues D and eigenvectors V as output from
      !  jacobi, this subroutine sorts the eigenvales into ascending
      !  order and rearranges the colmns of V accordingly.
      !
      ! The subroutine is taken from 'Numerical Recipes.'
      !
      implicit none
      !
      integer n,np,i,j,k
      !
      real*8 D(np),V(np,np),P
      

      do i=1,n-1
	k = i
	P = D(i)
	do j=i+1,n
	  if (D(j).ge.P) then
	    k = j
	    P = D(j)
	  end if
	end do
	if (k.ne.i) then
	  D(k) = D(i)
	  D(i) = P
	  do j=1,n
	    P = V(j,i)
	    V(j,i) = V(j,k)
	    V(j,k) = P
	  end do
  	end if
      end do
      

      return
      end subroutine eigsrt

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: SUBROUTINE matInv3:'
        write(*,*) 'WARNING: DET of MAT=',DET_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine matInv2D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse, and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(2,2),A_inv(2,2),det_A,det_A_inv

      
      istat = 1
      
      det_A = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: subroutine matInv2D:'
        write(*,*) 'WARNING: det of mat=',det_A
        istat = 0
        return
      end if
            
      det_A_inv = 1.d0/det_A
          
      A_inv(1,1) =  det_A_inv*A(2,2)
      A_inv(1,2) = -det_A_inv*A(1,2)
      A_inv(2,1) = -det_A_inv*A(2,1)
      A_inv(2,2) =  det_A_inv*A(1,1)


      return
      end subroutine matInv2D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

!****************************************************************************

      subroutine dlnxdx(X,DYDX)
      !
      ! This subroutine calculates the derivative of the logarithm
      ! of a symmetric tensor with respect to that tensor
      !
      implicit none
      !
      integer i,j,k,l
      !
      real*8 X(3,3),DYDX(3,3,3,3),Iden(3,3),Iden4(3,3,3,3),eigval(3),
     +  eigvec(3,3),ehat1(3),ehat2(3),ehat3(3),E1(3,3),E2(3,3),E3(3,3),
     +  y(3),DX2DX(3,3,3,3),s1,s2,s3,s4,s5,s6
      !
      real*8 zero,one,two,half,three,third,small
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,
     +     three=3.d0,third=1.d0/3.d0,small=1.d-12)
      
      
      ! Initialize
      !
      ! Second order identity tensor
      !
      call onem(Iden)
      !
      ! Fourth order symmetric identity tensor
      !
      Iden4 = zero
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              Iden4(i,j,k,l) = half*(Iden(i,k)*Iden(j,l) + 
     +                                Iden(i,l)*Iden(j,k))
            end do
          end do
        end do
      end do
      
      
      ! Calculate the eigenvalues and eigenvectors of X
      !
      call spectral(X,eigval,eigvec)
      !
      ! Extract the eigenvectors
      !
      do i=1,3
        ehat1(i) =  eigvec(i,1)
        ehat2(i) =  eigvec(i,2)
        ehat3(i) =  eigvec(i,3)
      end do
      !
      ! Assemble the eigenprojections
      !
      do i=1,3
        do j=1,3
	  E1(i,j) = ehat1(i)*ehat1(j)
	  E2(i,j) = ehat2(i)*ehat2(j)
	  E3(i,j) = ehat3(i)*ehat3(j)
	end do
      end do
      
      
      ! Calculate the eigenvalues of Y = ln(X)
      !
      y = dlog(eigval)
      
      
      ! Calculate the derivative of X^2 with respect to X
      !
      DX2DX = zero
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              DX2DX(i,j,k,l) = half*(Iden(i,k)*X(j,l) + 
     +                                Iden(i,l)*X(j,k) + 
     +                                X(i,k)*Iden(j,l) + 
     +                                X(i,l)*Iden(j,k))
            end do
          end do
        end do
      end do
         
            
      ! Calculate DYDX
      !
      DYDX = zero
      if (dabs(eigval(1)-eigval(3)).le.small) then
        !
        ! Three repeated eigenvalues
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = (one/eigval(1))*Iden4(i,j,k,l)
              end do
            end do
          end do
        end do
        !
      elseif (dabs(eigval(2)-eigval(3)).le.small) then
        !
        ! The eigenvalues 2 and 3 are repeated. Eigenvalue 1 is distinct.
        !
        s1 = (y(1) - y(2))/((eigval(1)-eigval(2))**two) - 
     +                  (one/eigval(2))/(eigval(1)-eigval(2))
        s2 = two*eigval(2)*(y(1)-y(2))/((eigval(1)-eigval(2))**two) - 
     +     (one/eigval(2))*(eigval(1)+eigval(2))/(eigval(1)-eigval(2))
        s3 = two*(y(1)-y(2))/((eigval(1)-eigval(2))**three) - 
     +   ((one/eigval(1))+(one/eigval(2)))/((eigval(1)-eigval(2))**two)
        s4 = eigval(2)*s3
        s5 = s4
        s6 = (eigval(2)**two)*s3
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = s1*DX2DX(i,j,k,l) - s2*Iden4(i,j,k,l) - 
     +                     s3*X(i,j)*X(k,l) + s4*X(i,j)*Iden(k,l) + 
     +                     s5*Iden(i,j)*X(k,l) - s6*Iden(i,j)*Iden(k,l)
              end do
            end do
          end do
        end do
        !
      elseif (dabs(eigval(1)-eigval(2)).le.small) then
        !
        ! The eigenvalues 1 and 2 are repeated. Eigenvalue 3 is distinct.
        !
        s1 = (y(3) - y(2))/((eigval(3)-eigval(2))**two) - 
     +                  (one/eigval(2))/(eigval(3)-eigval(2))
        s2 = two*eigval(2)*(y(3)-y(2))/((eigval(3)-eigval(2))**two) - 
     +     (one/eigval(2))*(eigval(3)+eigval(2))/(eigval(3)-eigval(2))
        s3 = two*(y(3)-y(2))/((eigval(3)-eigval(2))**three) - 
     +   ((one/eigval(3))+(one/eigval(2)))/((eigval(3)-eigval(2))**two)
        s4 = eigval(2)*s3
        s5 = s4
        s6 = (eigval(2)**two)*s3
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = s1*DX2DX(i,j,k,l) - s2*Iden4(i,j,k,l) - 
     +                     s3*X(i,j)*X(k,l) + s4*X(i,j)*Iden(k,l) + 
     +                     s5*Iden(i,j)*X(k,l) - s6*Iden(i,j)*Iden(k,l)
              end do
            end do
          end do
        end do
        !
      else
        !
        ! Eigenvalues are distinct.
        !
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                DYDX(i,j,k,l) = (y(1)/((eigval(1)-eigval(2))*
     +                                 (eigval(1)-eigval(3))))*
     +         (DX2DX(i,j,k,l) - (eigval(2)+eigval(3))*Iden4(i,j,k,l) - 
     +  ((eigval(1)-eigval(2))+(eigval(1)-eigval(3)))*E1(i,j)*E1(k,l) - 
     +       (eigval(2)-eigval(3))*(E2(i,j)*E2(k,l)-E3(i,j)*E3(k,l))) + 
     +                        (one/eigval(1))*E1(i,j)*E1(k,l) +
     +                          (y(2)/((eigval(2)-eigval(1))*
     +                                 (eigval(2)-eigval(3))))*
     +         (DX2DX(i,j,k,l) - (eigval(1)+eigval(3))*Iden4(i,j,k,l) - 
     +  ((eigval(2)-eigval(1))+(eigval(2)-eigval(3)))*E2(i,j)*E2(k,l) - 
     +       (eigval(1)-eigval(3))*(E1(i,j)*E1(k,l)-E3(i,j)*E3(k,l))) + 
     +                        (one/eigval(2))*E2(i,j)*E2(k,l) +
     +                          (y(3)/((eigval(3)-eigval(1))*
     +                                 (eigval(3)-eigval(2))))*
     +         (DX2DX(i,j,k,l) - (eigval(1)+eigval(2))*Iden4(i,j,k,l) - 
     +  ((eigval(3)-eigval(1))+(eigval(3)-eigval(2)))*E3(i,j)*E3(k,l) - 
     +       (eigval(1)-eigval(2))*(E1(i,j)*E1(k,l)-E2(i,j)*E2(k,l))) + 
     +                        (one/eigval(3))*E3(i,j)*E3(k,l)
              end do
            end do
          end do
        end do
        !
      end if
      
      return
      end subroutine dlnxdx

!****************************************************************************

      subroutine mprod4(A,B,C)
      !
      ! This subroutine calculates the product of two fourth order tensors,
      ! A and B, and places the product in C:
      !             C_{ijkl} = A_{ijmn}B_{mnkl}.
      !
      implicit none
      !
      integer i,j,k,l,m,n
      !
      real*8 A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
      
      
      C = 0.d0
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              do m=1,3
                do n=1,3
                  C(i,j,k,l) = C(i,j,k,l) + A(i,j,m,n)*B(m,n,k,l)
                end do
              end do
            end do
          end do
        end do
      end do
      
      
      return
      end subroutine mprod4

!****************************************************************************