!****************************************************************************
  ! This file containing:
  ! Subroutine solve_system: chose the linear algebraic solver by msolver
  ! Subroutine solve_info: print the information of the final system
  ! Subroutine compress: further remove all non-zero elements in the matrix
  ! Subroutine output: output the results
  ! Subroutine output_trans: output the results for transient heat or dynamics problems
  ! Subroutine evaluate_stress: calculate the heat flux or stress
  ! Subroutine sigtitl: generate the title for output the heat flux or stress
  ! Subroutine strength: calculate the strength criteria by stress
!****************************************************************************

subroutine solve_system
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  allocatable ipiv(:),asys(:,:)
  
  select case(iabs(msolver))
  case(1) ! pardiso
    !call denseoutput_csr_matrix(ntf,ntf,ka,ia,ja,a)
    call pardiso_unsym(ka,ia,ja,a,ntf,b,xval)
  case(2) ! dgesv:full-populated matrix solver in mkl
    ! transform a sparse matrix into a dense matrix
    allocate(asys(ntf,ntf1))
    asys=0.d0
    do i=1,ntf
      do j=ia(i),ia(i+1)-1
        asys(i,ja(j))=a(j)
      enddo
      asys(i,ntf1)=b(i)
    enddo
    allocate (ipiv(ntf))
    call dgesv(ntf,1,asys,ntf,ipiv,asys(1,ntf1),ntf,info)
    if(info.ne.0)then
      write(*,*)' *** The solver has an error with info =',info
    endif
    xval=asys(:,ntf1)
    deallocate(asys,ipiv)
  case(3) ! an old full-populated matrix solver
    ! transform a sparse matrix into a dense matrix
    allocate(asys(ntf,ntf1))
    asys=0.d0
    do i=1,ntf
      do j=ia(i),ia(i+1)-1
        asys(i,ja(j))=a(j)
      enddo
      asys(i,ntf1)=b(i)
    enddo
    call cmlib_lud(ntf,ntf1,asys,ntf,1,info)
    if(info.ne.0)then
      write(*,*)' *** The solver has an error with info =',info
    endif
    xval=asys(:,ntf1)
    deallocate(asys)
  case(4) ! dss
    call dss_unsym(ka,ia,ja,a,ntf,b,xval)
  case(5) ! ilu0 precondition fgmres without compress
    xval=0.d0
    call fgmres_ilu0(ka,ia,ja,a,ntf,b,xval)
  case(6) ! ilut precondition fgmres without compress
    write(*,*)'The maximum halfband is',maxhb
    write(4,*)'The maximum halfband is',maxhb
    call fgmres_ilut(ka,ia,ja,a,ntf,b,xval,maxhb)
  end select
  
end subroutine solve_system

subroutine solve_info
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  integer*8 memory
  
  if(iabs(msolver).eq.2.or.iabs(msolver).eq.3) then ! full-populated matrix solver
    if(istrans.eq.1)then
      write(*,*)'Transient problems cannot use this solver'
      write(4,*)'Transient problems cannot use this solver'
      stop
    endif
  else ! sparse matrix solver
    if(msolver.lt.0)then ! sparse matrix solver with compress
      call compress
    endif
    if(isnonli.eq.0)then
      memory=sizeof(integer)*(ka+ntf1)+sizeof(real*8)*ka
      write(*,*)
      write(4,*)
      write(*,1)ka,memory
      write(4,1)ka,memory
    endif
  endif
  1 format('There are ',i10,' non-zero elements and need ',i12, ' byte memory to store the matrix.')
end subroutine solve_info

subroutine compress ! further remove all non-zero elements in the matrix
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  
  integer,dimension (:), allocatable ::iap,jap
  real*8,dimension (:), allocatable ::ap
  
  allocate(iap(ntf1),jap(ka),ap(ka))
  
  kap=0
  iap(1)=1
  do i=1,ntf1-1
    do j=ia(i),ia(i+1)-1
      if(dabs(a(j)).lt.1.0d-15)cycle
      kap=kap+1
      jap(kap)=ja(j)
      ap(kap)=a(j)
    enddo
    iap(i+1)=kap+1
  enddo
  ka=kap
  deallocate(ja,a)
  allocate(ja(ka),a(ka))
  ia=iap
  ja=jap(1:ka)
  a=ap(1:ka)
  deallocate(iap,jap,ap)
  
end subroutine compress

subroutine output
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  character pch*8,uch*19
  dimension uch(3)
  allocatable disp(:,:),sigma(:,:),streng(:,:)
  data pch/'   Node'/,uch/'         Ux       ',         '           Uy       ','              Uz       '/
  write(7,*)
  write(7,*)'                   Output results :'
  write(7,*)
  if(ndf.gt.1)then
    write(7,'(1x,a8,6a22)')pch,(uch(i),i=1,ndf)
  else ! ndf.eq.1
    write(7,'(1x,a8,a22)')pch,'         T        '
  endif
  allocate (disp(ndf,ntp))
  if(ndf.eq.1) write(45,*) ntp
  do ip=1,ntp
    ip0=ndf*(ip-1)
    do i=1,ndf
      disp(i,ip)=xval(ip0+i)
    enddo
    write(7,'(i8,1x,1p3e24.16)')ip,real(disp(1:ndf,ip))
    if(ndf.eq.1) write(45,*)ip,disp(1,ip)
  enddo
  if(ndim.eq.1) return
  ! compute stresses and output result to tecplot file
  ncomp=1
  if(ndf.gt.1) ncomp=4
  allocate (sigma(nsig,ntp),streng(ncomp,ntp))
  call evaluate_stress(ncomp,disp,sigma,streng) 
  ! generate plot file
  call tecp_gen(33,disp,sigma,streng,cd,lnde,ntp,nte,ndim,ndf,nsig,mnode,ncomp)
  deallocate (disp,sigma,streng)  
  
  return
end subroutine output

subroutine output_trans(istep)
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  character pch*8,uch*19
  dimension uch(3)
  allocatable disp(:,:),sigma(:,:),streng(:,:)
  data pch/'   Node'/,uch/'       Ux        ','         Uy        ','           Uz         '/
  write(7,*)
  write(7,*)'            Output the results of ',istep,'-th step:'
  write(7,*)
  if(ndf.gt.1)then
    write(7,'(1x,a8,6a22)')pch,(uch(i),i=1,ndf)
  else ! ndf.eq.1
    write(7,'(1x,a8,a22)')pch,'         T        '
  endif
  allocate (disp(ndf,ntp))
  do ip=1,ntp
    ip0=ndf*(ip-1)
    do i=1,ndf
      disp(i,ip)=xval(ip0+i)
    enddo
    write(7,'(i8,1x,1p3e24.16)')ip,real(disp(1:ndf,ip)) ! internal points
  enddo
  if(ndim.eq.1) return
  ! compute stresses and output result to tecplot file
  ncomp=1
  if(ndf.gt.1) ncomp=4
  allocate (sigma(nsig,ntp),streng(ncomp,ntp))
  call evaluate_stress(ncomp,disp,sigma,streng) 
  ! generate plot file
  call tecp_gen_trans(33,99,disp,sigma,streng,cd,lnde,ntp,nte,ndim,ndf,nsig,mnode,ncomp,istep,deltat)
  deallocate (disp,sigma,streng)
  
end subroutine output_trans

subroutine evaluate_stress(ncomp,disp,sigma,streng)
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  dimension disp(ndf,ntp),sigma(nsig,ntp),sigmap(nsigp),xi(3),dndx(3,mnode),valu(ndim),bval(nsigp,ndf),db(nsigp,ndf),dval(nsigp,nsigp,mnode)
  dimension dndxx(3,3,mnode),elcod(3,mnode),c2(9),ct0(ndf,mnode),rho(mnode),cp(mnode),pr(mnode),streng(ncomp,ntp)
  allocatable tval(:)
  
  if(ndf.ne.1.and.istm.eq.1) then
    rewind (45)
    read(45,*) ktp
    allocate (tval(max(ntp,ktp)))
    do ik=1,ktp
      read(45,*)ip,tval(ip) ! retriev nodal temperatures
    enddo
  endif

  write(7,*)
  call sigtitl(ndim,ndf,nsig)
  sigma=0.d0
  bval=0.d0
  do ie=1,nte ! element loop
    ietype=netype(ie)
    node=nnode(ie)
    ndk=ndtok(node)
    
    ! get material of each node in the element
    if(isnonli.eq.0)then ! linear problem
      if(nmatgrp.lt.0)then ! material properties are specified by global nodes
        do knode=1,node
          elcod(1:ndim,knode)=cd(1:ndim,lnde(knode,ie))
          ip=lnde(knode,ie)
          dval(:,:,knode)=pred(ip,1,:,:)
          if(nsig.eq.4) pr(knode)=0.5d0*dval(1,2,knode)/(dval(1,2,knode)+dval(3,3,knode))
          if(istrans.eq.1)then
            rho(knode)=prerho(ip,1)
            cp(knode)=precp(ip,1)
          endif
          if(istm.eq.1) ct0(:,knode)=prect0(ip,1,:)
        enddo
      else ! material properties are specified by elements
        do knode=1,node
          elcod(1:ndim,knode)=cd(1:ndim,lnde(knode,ie))
          dval(:,:,knode)=pred(nematg(ie),knode,:,:)
          if(nsig.eq.4) pr(knode)=0.5d0*dval(1,2,knode)/(dval(1,2,knode)+dval(3,3,knode))
          if(istrans.eq.1)then
            rho(knode)=prerho(nematg(ie),knode)
            cp(knode)=precp(nematg(ie),knode)
          endif
          if(istm.eq.1) ct0(:,knode)=prect0(nematg(ie),knode,:)
        enddo
      endif
    else ! non-linear problem
      imat=nematg(ie)
      do knode=1,node
        elcod(1:ndim,knode)=cd(1:ndim,lnde(knode,ie))
        ip=lnde(knode,ie)
        call val_mu_nonli(ndf,ntf,xval,ip,metype(imat),ne0(imat),ce(1:ne0(imat),imat),e0)
        if(ndf.eq.1)then
          call dmatrix(e0,0.d0,0.d0,ct0(:,knode),0,1,c2(1),dval(:,:,knode))
        else ! ndf.gt.1
          call val_mu_nonli(ndf,ntf,xval,ip,mptype(imat),npr(imat),cpr(1:npr(imat),imat),pr(knode))
          if(istm.eq.0)then
            call dmatrix(e0,pr(knode),0.d0,ct0(:,knode),0,1,c2(1),dval(:,:,knode))
          else ! istm.eq.0
            call val_mu_nonli(ndf,ntf,xval,ip,mcttype(imat),nct(imat),cct(1:nct(imat),imat),ct)
            call dmatrix(e0,pr(knode),ct,ct0(:,knode),0,1,c2(1),dval(:,:,knode))
          endif
        endif
      enddo
    endif
    
    do inode=1,node ! source point loop
      nn=nrn(inode,ndk)
      ip=lnde(inode,ie)
      dndx(:,1:node)=dndxa(:,1:node,inode,ie)
      do jn=1,nn ! field point loop
        jd=ndrn(jn,inode,ndk)
        call coefb_gen(bval,dndx(:,jd))
        call matmulmkl(nsigp,nsigp,ndf,dval(:,:,inode),bval,db,1.d0,0.d0)
        call matmulmkl(nsigp,ndf,1,db,disp(:,lnde(jd,ie)),sigmap,1.d0,0.d0)
        !sigmap=matmul(matmul(dval(:,:,inode),bval),disp(:,lnde(jd,ie)))
        sigma(1:nsigp,ip)=sigma(1:nsigp,ip)+sigmap
        if(nsig.eq.4) then
          sigma(4,ip)=sigma(4,ip)+pr(inode)*(sigmap(1)+sigmap(2))
        endif
      enddo
      if(ndf.ne.1.and.istm.eq.1)then
        do idim=1,ndim
          sigma(idim,ip)=sigma(idim,ip)-ct0(idim,inode)*tval(ip) ! thermoelasticity
        enddo
      endif
    enddo ! end of source point loop
  enddo ! end of element loop
  
  ! output stress or heat flux
  if(ndf.eq.1) sigma=-sigma
  do ip=1,ntp
    do isig=1,nsig
      sigma(isig,ip)=sigma(isig,ip)/dble(neop(ip))
    enddo
    write(7,'(i8,1x,1p6e19.10)')ip,(sigma(isig,ip),isig=1,nsig)
  enddo
  
  ! output strength
  if(ndf.eq.1) then
    do ip=1,ntp
     streng(1,ip)=dsqrt(dot_product(sigma(1:nsig,ip),sigma(1:nsig,ip))) ! size of flux
    enddo
  else
    write(7,*)
    write(7,*)'           Values of strength for 4 criteria :'
    write(7,"(4x,'node',8x,'Tresca',12x,'Von-Mises',9x,'Mohr-Coulomb',5x,'Drucker-Prager')")
    do ip=1,ntp
      call strength(ndim,nsig,dlt,iny,sigma(1,ip),streng(1,ip))
      write(7,'(i8,1p4e19.10)')ip,streng(:,ip)
    enddo
  endif
  
end subroutine evaluate_stress

subroutine sigtitl(ndim,ndf,nsig)
  implicit real*8 (a-h,o-z)
  character pch*7,ss*2
  dimension ss(6)
  pch=' Node '
  if(ndf.ne.1) then
    ss=(/'XX','YY','ZZ','XY','YZ','ZX'/)
    if(nsig.eq.6) write(7,'(1x,a7,6(10x,''stress-'',a2))')pch,(ss(i),i=1,nsig)
    if(nsig.eq.4) write(7,'(1x,a7,1x,4(10x,''stress-'',a2))')pch,(ss(i),i=1,2),ss(4),ss(3)
    if(nsig.eq.3) write(7,'(1x,a7,1x,3(10x,''stress-'',a2))')pch,(ss(i),i=1,2),ss(4)
  else
    ss=(/'X ','Y ','Z ','XY','YZ','ZX'/)
    if(ndim.eq.2) write(7,'(3x,a7,2(10x,'' flux-'',a2))')pch,(ss(i),i=1,nsig)
    if(ndim.eq.3) write(7,'(3x,a7,3(10x,'' flux-'',a2))')pch,(ss(i),i=1,nsig)
  endif
end subroutine sigtitl

subroutine strength(ndim,nsig,dlt,iny,sigma,yield)
  implicit real*8 (a-h,o-z)
  dimension dlt(3,3),iny(9),sigma(nsig),devia(3,3),yield(4)
  ! parameters for mohr-coulomb and drucker-prager criteria
  uniax=0.8d0 ! cohesion c
  frict=20.0d0 ! internal friction angle phi
  root3=dsqrt(3.d0)
  frict=frict*0.017453292d0
  sinf=dsin(frict)
  alfa=2.d0*sinf/(root3*(3.d0-sinf))
  ceqs3=0.5d0*(1.d0-sinf)
  ceqs4=1.d0/root3-alfa
  ! evaluate mean stress
  smean=0.d0
  do i=1,ndim
    smean=smean+sigma(i)/3.d0
  enddo
  if(nsig.eq.4) smean=smean+sigma(4)/3.d0 ! plane-strain case
  ! evaluate deviatoric stresses
  k=0
  do i=1,ndim
    do j=1,ndim
      k=k+1
      devia(i,j)=sigma(iny(k))-smean*dlt(i,j) ! eq.(7.8) of gao's book
    enddo
  enddo
  ! find second and third invariants of deviatoric stresses
  varj2=0.d0
  varj3=0.d0
  do i=1,ndim
    do j=1,ndim
      varj2=varj2+0.5d0*devia(i,j)*devia(i,j) ! eq.(7.9)2
      do k=1,ndim
        varj3=varj3+devia(i,j)*devia(j,k)*devia(k,i)/3.d0 ! eq.(7.9)3
      enddo
    enddo
  enddo
  if(nsig.eq.6) goto 35
  ! particular treatment for two-dimensional problems
  if(nsig.eq.3) devia(3,3)=-smean ! above eq.(10.1)
  if(nsig.eq.4) devia(3,3)=sigma(4)-smean ! for plane strain
  varj2=varj2+0.5d0*devia(3,3)*devia(3,3) ! above eq.(10.2)
  varj3=varj3+devia(3,3)**3/3.d0 ! above eq.(10.3)
  35  squj2=dsqrt(varj2)
  if(squj2.lt.1.d-6) goto 40
  sint3=-3.d0*root3*varj3/(2.d0*varj2*squj2) ! eq.(7.11)
  if(dabs(sint3).gt.1.d0) sint3=dsign(1.d0,sint3)
  goto 50
  40  sint3=0.0d0
  50  theta=dasin(sint3)/3.d0 ! eq. (7.11)
  ! find equivalent stress for four criteria
  yield(1)=2.*dcos(theta)*squj2 ! tresca
  yield(2)=root3*squj2 ! von mises
  snphi=dsin(frict) ! mohr-coulomb
  yield(3)=smean*snphi+squj2*(dcos(theta)-dsin(theta)*snphi/root3)
  yield(3)=yield(3)/ceqs3 ! consistent form
  yield(4)=(alfa*3.*smean+squj2)/ceqs4 ! drucker-prager
end subroutine strength