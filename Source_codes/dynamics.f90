!****************************************************************************
  ! This file containing:
  ! Subroutine solve_output_for_dyn: solve the equation of dynamics problem step by step and output the results
  ! Subroutine ini_second_type_bc_dyn: add initial second type B.C. of dynamics
  ! Subroutine dyn_coef: form the matrix and initial vectors of dynamics problem
  ! Subroutine initial_vector: calculate the accelerate vector of initial step
  ! Subroutine pardiso_unsym_dyn: use PARDISO to solve the linear algebraic equation in dynamics problem (need MKL labrary in Intel FORTRAN or oneAPI)
  ! Subroutine dss_unsym_dyn: use DSS to solve the linear algebraic equation in dynamics problem (need MKL labrary in Intel FORTRAN or oneAPI)
  ! Subroutine initial_parameter: compute relative parameter of time marching scheme
!****************************************************************************

subroutine solve_output_for_dyn
  use fixed_values
  use vary_arrays
  use compute
  implicit real*8 (a-h,o-z)
  allocatable xval1(:),conve(:),er(:),x(:),xp(:),xpp(:),x1(:),xp1(:),xpp1(:)
  
  write(*,*)
  write(*,*)'         Dynamics response iteration'
  write(4,*)
  write(4,*)'         Dynamics response iteration'
  
  allocate (a(ka),b(ntf),isfirtyp(ntf),c(kc),rm(krm))
  a=0.d0
  b=0.d0
  isfirtyp=0
  c=0.d0
  rm=0.d0
  call coef_gen ! form original coeficient matrix
  call first_type_bc(0) ! add  first type of B.C.
  call third_type_bc(0) ! add third type of B.C.
  call source_item(0) ! add heat source item
    allocate(b0(ntf))
    b0=b
  call ini_second_type_bc_dyn ! add initial second type B.C. of dynamics
  allocate(x(ntf),xp(ntf),xpp(ntf))
  call dyn_coef(x,xp,xpp) ! prepare dynamics response iteration matrix and initial vector
  call solve_info ! information of matrix a
  
  open(99,file='results.qqq',status='unknown',form='unformatted') ! store data section of tecplot file
  rewind(99)
  
  if(iabs(msolver).eq.1) then ! pardiso
    call pardiso_unsym_dyn(istep,x,xp,xpp)
  elseif(iabs(msolver).eq.4) then ! dss
    call dss_unsym_dyn(istep,x,xp,xpp)
  else
    istep=0
    call output_trans(istep)
    call initial_parameter(deltat,mdyn,cdyn,ac0,ac1,ac2,ac3,ac4,ac5,ac6,ac7)
    allocate(x1(ntf),xp1(ntf),xpp1(ntf))
    allocate(xval1(ntf),conve(ntf),er(ndf))
    !-------------------------begin of transient iteration-------------------------
    do
      istep=istep+1
      ! record last step results
      xval1=xval
      x1=x
      xp1=xp
      xpp1=xpp
      ! get b of this step
      b=b0
      call second_type_bc(istep)
      call matmulmkl_mmk(ntf,ntf,irm,jrm,rm,krm,ac0*x+ac2*xp+ac3*xpp,b,1.d0,1.d0)
      call matmulmkl_mmk(ntf,ntf,ic,jc,c,kc,ac1*x+ac4*xp+ac5*xpp,b,1.d0,1.d0)
      x=0.d0
      ! solve
      if(iabs(msolver).eq.5) then ! ilu0 precondition fgmres without compress
        call fgmres_ilu0(ka,ia,ja,a,ntf,b,x)
      elseif(iabs(msolver).eq.6) then ! ilut precondition fgmres without compress
        call fgmres_ilut(ka,ia,ja,a,ntf,b,x,maxhb)
      endif
      ! get new xp, xpp, xval
      xpp=ac0*(x-x1)-ac2*xp1-ac3*xpp1
      xp=xp1+ac6*xpp1+ac7*xpp
      call matmulmkl_mmk(ntf,ntf,it,jtp,t,kt,x,xval,1.d0,0.d0)
      ! verify convergence and output the results of this step
      call output_trans(istep)
      call step_info(istep,1,ndf,ntf,conve,er,xval,xval1,1)
      if(tolt.lt.1)then ! judge by error
        if(maxval(er).lt.tolt)exit
      else ! judge by max step
        if(istep.ge.int(tolt))exit
      endif
    enddo
    !-------------------------end of transient iteration-------------------------
    deallocate(xval1,conve,er,x1,xp1,xpp1)
  endif
  deallocate(ia,ja,a,b,b0,ic,jc,c,irm,jrm,rm,x,xp,xpp)
  write(*,*)'Total transient iteration step=',istep
  write(4,*)'Total transient iteration step=',istep
  
  ! copy temporary file 99 to file 33 to complete data module of tecplot binary file
  call tecp_gen_trans2(33,99,ndim,ndf,nsig,mnode,ntp,nte,istep)
  
end subroutine solve_output_for_dyn

subroutine ini_second_type_bc_dyn ! add initial second type B.C. of dynamics
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  dimension ck(3,25),vtgrp(ndf)
  
  do ie=1,nte ! element loop
    ietype=netype(ie)
    node=nnode(ie)
    ndk=ndtok(node)
    do iface=1,nface(ndk) ! face loop
      igrp=nbcgrp(ie,iface)
      if(igrp.eq.0)cycle
      if(lret(igrp).ne.2.and.lret(igrp).ne.-2)cycle ! not second type of B.C.
      if(itef(igrp).eq.1)then
        acoef=tef(igrp)
        acoef=-acoef**2d0*dsin(acoef*0.d0)
      elseif(itef(igrp).eq.2)then
        acoef=tef(igrp)
        acoef=-acoef**2d0*dcos(acoef*0.d0)
      else
        cycle
      endif
      if(lret(igrp).eq.-2) then
        ! pression bc need the global coordinates of the nodes in the element
        do ipface=1,npface(ndk)
          inode=ndface(ipface,iface,ndk)
          ck(1:ndim,ipface)=cd(1:ndim,lnde(inode,ie))
        enddo
      endif
      do ipface=1,npface(ndk) ! point loop
        inode=ndface(ipface,iface,ndk)
        igp=lnde(inode,ie)
        ip0=(igp-1)*ndf
        ! get the value of second type bc and add the second type bc into system
        if(lret(igrp).eq.2) then ! normal b.c.
          vtgrp=pret(igrp,1:ndf)
        elseif(lret(igrp).eq.-2) then ! press is specified
          press=pret(igrp,1)
          call press_bc(ipface,ndim,ndf,npface(ndk),ck,press,vtgrp,info) ! treat press b.c.
        endif
        do idf=1,ndf
          if(kbcflag(idf,ie,iface).ne.2.or.isfirtyp(ip0+idf).eq.1) cycle ! jump the directions whose tractions are unknown or displacements are known
          b(ip0+idf)=b(ip0+idf)+vtgrp(idf)*acoef
        enddo
      enddo ! end of point loop
    enddo ! end of face loop
  enddo ! end of element loop
  
endsubroutine ini_second_type_bc_dyn

subroutine dyn_coef(x,xp,xpp) ! form the matrix and initial vectors of dynamics problem
  use fixed_values
  use vary_arrays
  use compute
  implicit real*8 (a-h,o-z)
  dimension x(ntf),xp(ntf),xpp(ntf)
  allocatable iap(:),jap(:),ap(:),bp(:)
  
  ! reorder a
  allocate(iap(ntf+1))
  call matmulmkl_mmm(ntf,ntf,ntf,it,jt,t,kt,ia,ja,a,ka,iap,jap,ap,kap)
  deallocate(ja,a)
  call matmulmkl_mmm(ntf,ntf,ntf,iap,jap,ap,kap,it,jtp,t,kt,ia,ja,a,ka)
  deallocate(iap,jap,ap)
  !call directoutput_csr_matrix(ntf,ntf,ka,ia,ja,a)
  ! reorder b
  allocate(bp(ntf))
  call matmulmkl_mmk(ntf,ntf,it,jt,t,kt,b,bp,1.d0,0.d0)
  b=bp
  deallocate(bp)
  ! get initial vector
  call initial_vector(x,xp,xpp)
  ! compute equivalent stiffness and stored into a
  ac0=1.d0/cdyn(2)/(deltat)**2.d0
  ac1=cdyn(1)/cdyn(2)/deltat
  call addmkl(ntf,ntf,ia,ja,a,ka,irm,jrm,ac0*rm,krm)
  call addmkl(ntf,ntf,ia,ja,a,ka,ic,jc,ac1*c,kc)
  
end subroutine dyn_coef

subroutine initial_vector(x,xp,xpp)
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  dimension x(ntf),xp(ntf),xpp(ntf)
  allocatable f0(:),irmp(:),jrmp(:),rmp(:)
  
  ! compute initial displacement and velocity
  call matmulmkl_mmk(ntf,ntf,it,jt,t,kt,xval,x,1.d0,0.d0)
  call matmulmkl_mmk(ntf,ntf,it,jt,t,kt,xvalp,xp,1.d0,0.d0)
  ! compute initial acceleration
  allocate(irmp(ntf+1))
  irmp(1:nobdf*ndf+1)=ia(1:nobdf*ndf+1)
  irmp(nobdf*ndf+2:ntf+1)=irm(nobdf*ndf+2:ntf+1)+ia(nobdf*ndf+1)-1
  krmp=irmp(ntf+1)-1
  allocate(jrmp(krmp),rmp(krmp))
  do i=1,ia(nobdf*ndf+1)-1
    jrmp(i)=ja(i)
    rmp(i)=a(i)
  enddo
  do i=1,krm
    jrmp(i+ia(nobdf*ndf+1)-1)=jrm(i)
    rmp(i+ia(nobdf*ndf+1)-1)=rm(i)
  enddo
  allocate(f0(ntf))
  f0=b
  call matmulmkl_mmk(ntf,ntf,ia,ja,a,ka,x,f0,-1.d0,1.d0)
  call matmulmkl_mmk(ntf,ntf,ic,jc,c,kc,xp,f0,-1.d0,1.d0)
  call pardiso_unsym(krmp,irmp,jrmp,rmp,ntf,f0,xpp)
  deallocate(f0,irmp,jrmp,rmp)
  
end subroutine initial_vector

subroutine pardiso_unsym_dyn(istep,x,xp,xpp)
  use fixed_values
  use vary_arrays
  implicit none
  
  integer i,istep,itf
  integer maxfct, mnum, mtype, phase, nrhs, n, error, msglvl
  integer iparm(64)
  integer idum(1)
  integer*8 pt(64)
  real*8 ac0,ac1,ac2,ac3,ac4,ac5,ac6,ac7
  real*8  ddum(1)
  real*8  x(ntf),xp(ntf),xpp(ntf)
  real*8,allocatable::xval1(:),conve(:),er(:),x1(:),xp1(:),xpp1(:)

  data nrhs /1/, maxfct /1/, mnum /1/
  
  ! set up pardiso control parameter
  mtype = 11 ! real unsymmetric
  error = 0 ! initialize error flag
  msglvl = 0 ! print statistical information
  call pardisoinit (pt, mtype, iparm) ! initiliaze pardiso solver
  iparm(1) = 1 ! no solver default
  iparm(2) = 3 ! fill-in reordering from metis
  iparm(8) = 9 ! numbers of iterative refinement steps
  n=ntf
  
  ! reordering and symbolic factorization, this step also allocates
  ! all memory that is necessary for the factorization
  phase = 11 ! only reordering and symbolic factorization
  call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
  if(error.ne.0)then
    call par_er(error,iparm(15),iparm(16),iparm(17),1)
  endif
  
  ! factorization.
  phase = 22 ! only factorization
  call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
  if(error.ne.0)then
    call par_er(error,iparm(15),iparm(16),iparm(17),2)
  endif
  
  ! back substitution and iterative refinement
  iparm(8) = 2 ! max numbers of iterative refinement steps
  phase = 33 ! only solve
  
  istep=0 ! transient step
  call output_trans(istep)
  call initial_parameter(deltat,mdyn,cdyn,ac0,ac1,ac2,ac3,ac4,ac5,ac6,ac7)
  allocate(x1(ntf),xp1(ntf),xpp1(ntf))
  allocate(xval1(ntf),conve(ntf),er(ndf))
  !-------------------------begin of transient iteration-------------------------
  do
    istep=istep+1
    ! record last step results
    xval1=xval
    x1=x
    xp1=xp
    xpp1=xpp
    ! get b of this step
    b=b0
    call second_type_bc(istep)
    call matmulmkl_mmk(ntf,ntf,irm,jrm,rm,krm,ac0*x+ac2*xp+ac3*xpp,b,1.d0,1.d0)
    call matmulmkl_mmk(ntf,ntf,ic,jc,c,kc,ac1*x+ac4*xp+ac5*xpp,b,1.d0,1.d0)
    x=0.d0
    ! solve
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, nrhs, iparm, msglvl, b, x, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),3)
    endif
    ! get new xp, xpp, xval
    xpp=ac0*(x-x1)-ac2*xp1-ac3*xpp1
    xp=xp1+ac6*xpp1+ac7*xpp
    call matmulmkl_mmk(ntf,ntf,it,jtp,t,kt,x,xval,1.d0,0.d0)
    ! verify convergence and output the results of this step
    call output_trans(istep)
    call step_info(istep,1,ndf,ntf,conve,er,xval,xval1,1)
    if(tolt.lt.1)then ! judge by error
      if(maxval(er).lt.tolt)exit
    else ! judge by max step
      if(istep.ge.int(tolt))exit
    endif
  enddo
  !-------------------------end of transient iteration-------------------------
  deallocate(xval1,conve,er,x1,xp1,xpp1)
  
  ! termination and release of memory
  phase = -1 ! release internal memory
  call pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)
  if(error.ne.0)then
    call par_er(error,iparm(15),iparm(16),iparm(17),4)
  endif

end subroutine pardiso_unsym_dyn

subroutine dss_unsym_dyn(istep,x,xp,xpp)
  use fixed_values
  use vary_arrays
  implicit none
  include 'mkl_dss.fi'
  
  integer nrhs,istep,itf,n
  parameter (nrhs=1)
  integer idum(1)
  integer*8 handle
  integer error
  character*100 statin
  double precision statout(6)
  integer buflen
  parameter(buflen=100)
  integer buff(buflen)
  double precision ac0,ac1,ac2,ac3,ac4,ac5,ac6,ac7
  double precision x(ntf),xp(ntf),xpp(ntf)
  double precision,allocatable::xval1(:),conve(:),er(:),x1(:),xp1(:),xpp1(:)
  n=ntf
  
  ! initialize the solver
  error=dss_create(handle,mkl_dss_defaults)
  if(error.ne.mkl_dss_success)then
    write(*,*)'Error code',error,'occurred in dss_create.'
    write(4,*)'Error code',error,'occurred in dss_create.'
    stop
  endif
  
  ! define the non-zero structure of the matrix
  error=dss_define_structure(handle,mkl_dss_non_symmetric,ia,n,n,ja,ka)
  if(error.ne.mkl_dss_success)then
    write(*,*)'Error code',error,'occurred in dss_define_structure.'
    write(4,*)'Error code',error,'occurred in dss_define_structure.'
    stop
  endif
  
  ! reorder the matrix
  error=dss_reorder(handle,mkl_dss_defaults,idum)
  !error=dss_reorder(handle,mkl_dss_metis_openmp_order,idum)
  !error=dss_reorder(handle,mkl_dss_get_order,idum)
  if(error.ne.mkl_dss_success)then
    write(*,*)'Error code',error,'occurred in dss_reorder.'
    write(4,*)'Error code',error,'occurred in dss_reorder.'
    stop
  endif

  ! factor the matrix
  error=dss_factor_real(handle,mkl_dss_defaults,a)
  if(error.ne.mkl_dss_success)then
    write(*,*)'Error code',error,'occurred in dss_factor_real.'
    write(4,*)'Error code',error,'occurred in dss_factor_real.'
    stop
  endif
  
  ! get the solution vector
  istep=0 ! transient step
  call output_trans(istep)
  call initial_parameter(deltat,mdyn,cdyn,ac0,ac1,ac2,ac3,ac4,ac5,ac6,ac7)
  allocate(x1(ntf),xp1(ntf),xpp1(ntf))
  allocate(xval1(ntf),conve(ntf),er(ndf))
  !-------------------------begin of transient iteration-------------------------
  do
    istep=istep+1
    ! record last step results
    xval1=xval
    x1=x
    xp1=xp
    xpp1=xpp
    ! get b of this step
    b=b0
    call second_type_bc(istep)
    call matmulmkl_mmk(ntf,ntf,irm,jrm,rm,krm,ac0*x+ac2*xp+ac3*xpp,b,1.d0,1.d0)
    call matmulmkl_mmk(ntf,ntf,ic,jc,c,kc,ac1*x+ac4*xp+ac5*xpp,b,1.d0,1.d0)
    x=0.d0
    ! solve
    error=dss_solve_real(handle,mkl_dss_defaults,b,nrhs,x)
    if(error.ne.mkl_dss_success)then
      write(*,*)'Error code',error,'occurred in dss_solve_real.'
      write(*,*)'Error step=',istep
      write(4,*)'Error code',error,'occurred in dss_solve_real.'
      write(4,*)'Error step=',istep
      stop
    endif
    ! get new xp, xpp, xval
    xpp=ac0*(x-x1)-ac2*xp1-ac3*xpp1
    xp=xp1+ac6*xpp1+ac7*xpp
    call matmulmkl_mmk(ntf,ntf,it,jtp,t,kt,x,xval,1.d0,0.d0)
    ! verify convergence and output the results of this step
    call output_trans(istep)
    call step_info(istep,1,ndf,ntf,conve,er,xval,xval1,1)
    if(tolt.lt.1)then ! judge by error
      if(maxval(er).lt.tolt)exit
    else ! judge by max step
      if(istep.ge.int(tolt))exit
    endif
  enddo
  !-------------------------end of transient iteration-------------------------
  deallocate(xval1,conve,er,x1,xp1,xpp1)
  
  ! print determinant of the matrix (no statistics for a diagonal matrix)
  statin='REORDERTIME,FACTORTIME,SOLVETIME,PEAKMEM,FACTORMEM,SOLVEMEM'
  call mkl_cvt_to_null_terminated_str(buff,buflen,statin)
  error=dss_statistics(handle,mkl_dss_defaults,buff,statout)
  write(4,"(' The reorder time is',f10.3,' second')") statout(1)
  write(4,"(' The factor time is',f10.3,' second')") statout(2)
  write(4,"(' The solve time is', f10.3,' second')") statout(3)
  write(4,"(' The peak memory is', f15.3,' kb')") statout(4)
  write(4,"(' The factor memory is',f15.3,' kb')") statout(5)
  write(4,"(' The solve memory is',f15.3,' kb')") statout(6)
  
  ! deallocate solver storage
  error = dss_delete( handle, mkl_dss_defaults )
  if(error.ne.mkl_dss_success)then
    write(*,*)'Error code',error,'occurred in dss_delete.'
    write(4,*)'Error code',error,'occurred in dss_delete.'
    stop
  endif
  
end subroutine dss_unsym_dyn

subroutine initial_parameter(deltat,mdyn,cdyn,ac0,ac1,ac2,ac3,ac4,ac5,ac6,ac7)
  implicit real*8 (a-h,o-z)
  dimension cdyn(5)
  
  ! compute relative parameter
  select case(mdyn)
  case(1) ! Newmark
    ac0=1.d0/cdyn(2)/(deltat)**2.d0
    ac1=cdyn(1)/cdyn(2)/deltat
    ac2=1.d0/cdyn(2)/deltat
    ac3=5.d-1/cdyn(2)-1.d0
    ac4=cdyn(1)/cdyn(2)-1.d0
    ac5=(cdyn(1)/cdyn(2)-2.d0)*deltat/2.d0
    ac6=(1-cdyn(1))*deltat
    ac7=cdyn(1)*deltat
  case(2) ! wilson-theta
    ! to be developed
  case(3) ! central difference
    ! to be developed
  end select
  
end subroutine initial_parameter