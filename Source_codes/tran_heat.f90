!****************************************************************************
  ! This file containing:
  ! Subroutine solve_output_for_trans_heat: solve the equation of transient heat problem step by step and output the results
  ! Subroutine trans_heat_coef: generate matrix of transient heat problem from steady heat problem's matrix
  ! Subroutine pardiso_unsym_trans_heat: use PARDISO to solve the linear algebraic equation in transient heat problem (need MKL labrary in Intel FORTRAN or oneAPI)
  ! Subroutine dss_unsym_trans_heat: use DSS to solve the linear algebraic equation in transient heat problem (need MKL labrary in Intel FORTRAN or oneAPI)
!****************************************************************************

subroutine solve_output_for_trans_heat
  use fixed_values
  use vary_arrays
  use compute
  implicit real*8 (a-h,o-z)
  real*8  xt(ntf)
  allocatable xval1(:),conve(:),er(:)
  
  write(*,*)
  write(*,*)'         Transient heat iteration'
  write(4,*)
  write(4,*)'         Transient heat iteration'
  
  allocate (a(ka),b(ntf),isfirtyp(ntf),c(kc),r(kr))
  a=0.d0
  b=0.d0
  isfirtyp=0
  c=0.d0
  r=0.d0
  call coef_gen ! form coeficient matrix of steady heat problem
  call first_type_bc(0) ! add  first type of B.C.
  call third_type_bc(0) ! add third type of B.C.
  call source_item(0) ! add heat source item
  call solve_info ! information of matrix a
  call second_type_bc(0) ! add  second type of B.C.
  call trans_heat_coef ! prepare transient heat iteration matrix
  allocate(b0(ntf))
  b0=b
  open(99,file='results.qqq',status='unknown',form='unformatted') ! store data section of tecplot file
  rewind(99)
  
  if(iabs(msolver).eq.1) then ! pardiso
    call pardiso_unsym_trans_heat(istep)
  elseif(iabs(msolver).eq.4) then ! dss
    call dss_unsym_trans_heat(istep)
  else
    istep=0
    call output_trans(istep)
    allocate(xval1(ntf),conve(ntf),er(ndf))
    !-------------------------begin of transient iteration-------------------------
    do
      istep=istep+1
      ! record last step results for comparision
      xval1=xval
      ! get b of this step
      b=b0
      call matmulmkl_mmk(ntf,ntf,ir,jr,r,kr,xval,b,-1.d0,1.d0)
      xval=0.d0
      ! solve
      if(iabs(msolver).eq.5) then ! ilu0 precondition fgmres without compress
        call fgmres_ilu0(ka,ia,ja,a,ntf,b,xval)
      elseif(iabs(msolver).eq.6) then ! ilut precondition fgmres without compress
        call fgmres_ilut(ka,ia,ja,a,ntf,b,xval,maxhb)
      endif
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
    deallocate(xval1,conve,er)
  endif
  deallocate(ia,ja,a,b,b0,ic,jc,c,ir,jr,r)
  write(*,*)'Total transient iteration step=',istep
  write(4,*)'Total transient iteration step=',istep
  
  ! write file 45 temperature.dat
  if(ndf.eq.1) then
    write(45,*) ntp
    do ip=1,ntp
      write(45,*)ip,xval(ip)
    enddo
  endif
  deallocate(xval)
  ! copy temporary file 99 to file 33 to complete data module of tecplot binary file
  call tecp_gen_trans2(33,99,ndim,ndf,nsig,mnode,ntp,nte,istep)
  
end subroutine solve_output_for_trans_heat
  
subroutine trans_heat_coef ! form the matrix of transient heat problem
  use fixed_values
  use vary_arrays
  use compute
  implicit real*8 (a-h,o-z)
  ! forming R matrix
  r=(1.d0-theta)*a
  call addmkl(ntf,ntf,ir,jr,r,kr,ic,jc,c/deltat,kc)
  ! forming A matrix
  a=theta*a
  call addmkl(ntf,ntf,ia,ja,a,ka,ic,jc,-c/deltat,kc)
  
end subroutine trans_heat_coef

subroutine pardiso_unsym_trans_heat(istep)
  use fixed_values
  use vary_arrays
  implicit none
  
  integer i,istep
  integer maxfct, mnum, mtype, phase, nrhs, n, error, msglvl
  integer iparm(64)
  integer idum(1)
  integer*8 pt(64)
  real*8  ddum(1)
  real*8,allocatable::xval1(:),conve(:),er(:)

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
  allocate(xval1(ntf),conve(ntf),er(ndf))
  !-------------------------begin of transient iteration-------------------------
  do
    istep=istep+1
    ! record last step results for comparision
    xval1=xval
    ! get b of this step
    b=b0
    call matmulmkl_mmk(ntf,ntf,ir,jr,r,kr,xval,b,-1.d0,1.d0)
    xval=0.d0
    ! solve
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, nrhs, iparm, msglvl, b, xval, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),3)
    endif
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
  deallocate(xval1,conve,er)
  
  ! termination and release of memory
  phase = -1 ! release internal memory
  call pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)
  if(error.ne.0)then
    call par_er(error,iparm(15),iparm(16),iparm(17),4)
  endif
  
end subroutine pardiso_unsym_trans_heat

subroutine dss_unsym_trans_heat(istep)
  use fixed_values
  use vary_arrays
  implicit none
  include 'mkl_dss.fi'
  
  integer nrhs,istep,n
  parameter (nrhs=1)
  integer idum(1)
  integer*8 handle
  integer error
  character*100 statin
  double precision statout(6)
  integer buflen
  parameter(buflen=100)
  integer buff(buflen)
  double precision,allocatable::xval1(:),conve(:),er(:)
  
  ! initialize the solver
  n=ntf
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
  allocate(xval1(ntf),conve(ntf),er(ndf))
   !-------------------------begin of transient iteration-------------------------
  do
    istep=istep+1
    ! record last step results for comparision
    xval1=xval
    ! get b of this step
    b=b0
    call matmulmkl_mmk(ntf,ntf,ir,jr,r,kr,xval,b,-1.d0,1.d0)
    xval=0.d0
    ! solve
    error=dss_solve_real(handle,mkl_dss_defaults,b,nrhs,xval)
    if(error.ne.mkl_dss_success)then
      write(*,*)'Error code',error,'occurred in dss_solve_real.'
      write(*,*)'Error step=',istep
      write(4,*)'Error code',error,'occurred in dss_solve_real.'
      write(4,*)'Error step=',istep
      stop
    endif
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
  deallocate(xval1,conve,er)
  
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
  
end subroutine dss_unsym_trans_heat