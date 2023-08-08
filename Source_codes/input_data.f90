!****************************************************************************
  ! This file containing:
  ! Subroutine input_data: read the input data
  ! Subroutine read_face_bc: read the faces that B.C.s are specified
  ! Subroutine read_bc_value: read the values of the specified B.C.
  ! Subroutine bc_mark: take part a mark of boundary conditions into a array
  ! Subroutine read_body: read source term information
  ! Subroutine read_material: read material properties, type and distributions
  ! Subroutine read_nonli_material: read non-linear material properties, type and distributions
  ! Subroutine dmatrix: generate the material matrix of one node
  ! Function ie_type: get type number of the element
  ! Function val_mu: calculate the material parameters
!****************************************************************************

subroutine input_data
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  character uch*2,tch*10,vch*2,cdch*13,ndch*5
  dimension uch(3),tch(3),vch(3),cdch(3)
  dimension npout(ndf),spu(2*ndf),vtgrp(ndf,mnode)
  dimension c1(ndim),c2(ndim)
  
  data uch/'ux','uy','uz'/,tch/'Component1','Component2','Component3'/,cdch/'X-coordinate ','Y-coordinate ','Z-coordinate '/,ndch/' Node'/,vch/'vx','vy','vz'/
  
  write(7,'(/4x,''node '',3a20)')(cdch(i),i=1,ndim)
  do n=1,ntp ! input nodal coordinates
    read(5,*)m,(cd(j,m),j=1,ndim)
    write(7,'(i8,1x,1p3e21.12)')m,(cd(j,m),j=1,ndim)
  enddo

  write(7,'(/3x,''elem   mat_group    nnode'',8(/,8x,8(1x,a5,i3)))') (ndch,i,i=1,mnode)
  avl=0.d0 ! max size of an element
  ! input element nodes and specified traction group and direction flag
  do ie=1,nte 
    read(5,*)m,iprop,ndin
    netype(m)=ie_type(ndin,ndim,node)
    nnode(m)=node
    nematg(m)=iprop
    backspace (5)
    read(5,*)m,iprop,ndin,(lnde(id,m),id=1,ndin)
    write(7,'(3i8,16(/,8x,8i8))')m,iprop,ndin,(lnde(id,m),id=1,node)
    c1=cd(1:ndim,lnde(1,m))
    c2=cd(1:ndim,lnde(node,m))
    avle=0.d0
    do i=1,ndim
      avle=avle+(c2(i)-c1(i))**2
    enddo
    avl=max(avl,dsqrt(avle))
  enddo
  timew='Read nodes and elements'
  call time_record(timew,datev,timev,4,2)
  
  ! read boundary conditions
  nbcgrp=0 ! B.C. group of each face
  kbcflag=0 ! B.C. flag of each face
  read(5,*)ntface
  write(7,'(/5x,''Total number of faces that are given B.C.s='',i5)') ntface
  call read_face_bc(ntface,ndf,nte,nfacem,nbcgrp,kbcflag)
  
  ! read prescribed B.C. values
  read(5,*) nbcv
  write(7,'(/3x,''Speicified B.C. groups ='',i5)') nbcv
  allocate(pret(nbcv,3),lret(nbcv),tef(nbcv),itef(nbcv))
  write(7,'(/2x,"Group   B.C. type ",3a16)')(tch(i),i=1,ndf)
  call read_bc_value(nbcv,ndf,pret,lret,itef,tef)
  
  ! read body forces
  read(5,*) nbody
  allocate(ibody(nbody),mbody(nbody),kbdelem(nte),vbody(10,nbody))
  call read_body(nbody,nte,ibody,mbody,kbdelem,vbody)
  
  ! read material properties
  read(5,*) nmatgrp,istrans,isnonli,istm ! whether the problem is transient, non-linerar or thermal-mechanics
  write(7,*)
  write(7,*) '   The number of material types =',nmatgrp
  write(7,*)
  if(isnonli.eq.0) then
    call read_material
  else
    call read_nonli_material
  endif
  
  ! read parameter of non-linear problem
  if(isnonli.eq.1) then
    read(5,*) omg,toln
    write(7,*)
    write(7,*)'Non-linear parameter:'
    if(toln.lt.1)then ! judge by error
      write(7,'("   omiga =",e15.6,"   toln =",e15.6)')omg,toln
    else ! judge by max step
      write(7,'("   omiga =",e15.6,"   toln =",i7)')omg,int(toln)
    endif
  endif
  
  ! read transient parameters
  if(istrans.eq.1) then
    if(ndf.eq.1)then ! read theta, deltat and tolerance of heat transient problem
      read(5,*) theta,deltat,tolt
      write(7,*)
      write(7,*)'Transient heat parameter:'
      if(tolt.lt.1)then ! judge by error
        write(7,'("   theta =",e15.6,"   deltat =",e15.6,"   tolt =",e15.6)')theta,deltat,tolt
      else ! judge by max step
        write(7,'("   theta =",e15.6,"   deltat =",e15.6,"   tolt =",i7)')theta,deltat,int(tolt)
      endif
    else ! read parameter of structure dynamic problem
      allocate(cdyn(5))
      read(5,*) mdyn,ndyn,cdyn(1:ndyn),deltat,tolt
      write(7,*)
      write(7,*)'Transient dynamic parameter:'
      select case(mdyn)
      case(1) ! newmark
        write(7,'("Newmark parameter delta =",e15.6)')cdyn(1)
        write(7,'("Newmark parameter alpha =",e15.6)')cdyn(2)
        write(7,'("deltat =",e15.6)')deltat
        if(toln.lt.1)then ! judge by error
          write(7,'("tolt =",e15.6)')tolt
        else ! judge by max step
          write(7,'("tolt =",i7)')int(tolt)
        endif
      case(2) ! wilson-theta, to be developed
        !write(7,'("wilson-theta parameter theta =",e15.6)')cdyn(1)
        !write(7,'("deltat =",e15.6)')deltat
        !if(toln.lt.1)then ! judge by error
        !  write(7,'("tolt =",e15.6)')tolt
        !else ! judge by max step
        !  write(7,'("tolt =",i7)')int(tolt)
        !endif
      case(3) ! central difference, to be developed
        !write(7,'("deltat =",e15.6)')deltat
        !if(toln.lt.1)then ! judge by error
        !  write(7,'("tolt =",e15.6)')tolt
        !else ! judge by max step
        !  write(7,'("tolt =",i7)')int(tolt)
        !endif
      end select
    endif
  endif
    
  if(istrans.eq.1.or.isnonli.eq.1) then
    ! read the initial condition or nonlinear initial vector
    if(istrans.eq.1.and.ndf.ne.1)then ! dynamics
      read(5,*)ninit
      write(7,'(/3x,''ip-bgn   ip-end  ip-step'',6a13)')(uch(i),i=1,ndf),(vch(i),i=1,ndf)
      allocate(xval(ntf),xvalp(ntf)) ! used to store initial conditions
      xval=0.d0
      xvalp=0.d0
      do iinit=1,ninit
        read(5,*)ipbgn,ipend,ipstp,(spu(j),j=1,2*ndf)
        write(7,'(1x,3i8,6e15.6)')ipbgn,ipend,ipstp,(spu(j),j=1,2*ndf)
        do ip=ipbgn,ipend,ipstp
          ip0=ndf*(ip-1)
          do idf=1,ndf
            xval(ip0+idf)=spu(idf)
            xvalp(ip0+idf)=spu(ndf+idf)
          enddo
        enddo
      enddo
    else
      read(5,*)ninit
      write(7,'(/3x,''ip-bgn   ip-end  ip-step  '',3a10)')(uch(i),i=1,ndf)
      allocate(xval(ntf)) ! used to store initial conditions
      xval=0.d0
      do iinit=1,ninit
        read(5,*)ipbgn,ipend,ipstp,(spu(j),j=1,ndf)
        write(7,'(1x,3i8,3e15.6)')ipbgn,ipend,ipstp,(spu(j),j=1,ndf)
        do ip=ipbgn,ipend,ipstp
          ip0=ndf*(ip-1)
          do idf=1,ndf
            xval(ip0+idf)=spu(idf)
          enddo
        enddo
      enddo
    endif
  endif
  
  timew='Get boundary conditions and material:'
  call time_record(timew,datev,timev,4,2)
  
  return
end subroutine input_data

integer function ie_type(node,nbdm,ndout)
  ndout=iabs(node)
  if(nbdm.eq.1) then ! line elements
    ie_type=node
  elseif(nbdm.eq.2) then ! surface elements
    ie_type=300+node
  else ! 3d volume cells
    ie_type=400+node
  endif
end function ie_type

subroutine read_face_bc(ntface,ndf,nte,nfacem,nbcgrp,kbcflag)
  ! read face that are given B.C.
  implicit real*8 (a-h,o-z)
  dimension nbcgrp(nte,nfacem),npout(ndf),kbcflag(ndf,nte,nfacem)
  
  write(7,'(/3x,''ie-bgn   ie-end  ie-step   face   group   mark'')')
  do iface=1,ntface
    read(5,*)iebgn,ieend,iestp,kface,ktgrp,kflag
    write(7,'(1x,6i8)')iebgn,ieend,iestp,kface,ktgrp,kflag
    do ie=iebgn,ieend,iestp
      nbcgrp(ie,kface)=ktgrp
      call bc_mark(ndf,kflag,npout)
      do idf=1,ndf
        if(kbcflag(idf,ie,kface).eq.0) kbcflag(idf,ie,kface)=npout(idf) ! only the first assignment is valid
      enddo
    enddo
  enddo
  
end subroutine read_face_bc

subroutine read_bc_value(nbcv,ndf,pret,lret,itef,tef)
  ! read prescribed traction information
  implicit real*8 (a-h,o-z)
  dimension pret(nbcv,3),lret(nbcv),tef(nbcv),itef(nbcv)
  character*3 :: func
  character*180 :: text
  
  itef=0
  do n=1,nbcv
    read(5,*)igrp,ibctype
    if(abs(igrp).gt.nbcv) then
      write(*,*)'The number of a B.C. group cannot exceed nbcv.'
      write(7,*)'The number of a B.C. group cannot exceed nbcv.'
      stop
    endif
    backspace (5)
    select case(ibctype)
    case(1) ! first type of B.C.
      read(5,*)igrp,lret(igrp),(pret(igrp,j),j=1,ndf)
      write(7,'(2(i5,5x),1p3e16.8)')igrp,lret(igrp),(pret(igrp,j),j=1,ndf)
    case(2) ! second type of B.C.
      read(5,*)igrp,lret(igrp),(pret(igrp,j),j=1,ndf)
      write(7,'(2(i5,5x),1p3e16.8)')igrp,lret(igrp),(pret(igrp,j),j=1,ndf)
    case(-2) ! press, only for mechanics problem
      read(5,*)igrp,lret(igrp),pret(igrp,1)
      write(7,'(2(i5,5x),1pe16.8)')igrp,lret(igrp),pret(igrp,1)
    case(3) ! third type of B.C., only for heat conduction problems now
      read(5,*)igrp,lret(igrp),pret(igrp,1),pret(igrp,2)
      write(7,'(2(i5,5x),1p2e16.8)')igrp,lret(igrp),pret(igrp,1),pret(igrp,2)
    end select
    read(5,'(a180)')text
    func=trim(adjustl(text))
    if(func.eq.'sin'.or.func.eq.'cos'.or.func.eq.'SIN'.or.func.eq.'COS')then
      write(7,'(5x,a3)')func
      if(func.eq.'sin'.or.func.eq.'SIN')then
        itef(iabs(igrp))=1
      elseif(func.eq.'cos'.or.func.eq.'COS')then
        itef(iabs(igrp))=2
      endif
      read(5,*)tef(iabs(igrp))
      write(7,'(5x,f11.6)')tef(iabs(igrp))
    else
      backspace(5)
    endif
  enddo
  
end subroutine read_bc_value

subroutine bc_mark(ndf,kinp,kout)
  ! take part a integer to an array
  implicit real*8 (a-h,o-z)
  dimension kout(ndf)
  if(kinp.gt.10**ndf)then
    write(*,'("The ndf should be",i4)')ndf
    write(*,'("kbcflag is wrong, please check input file!")')
    write(7,'("The ndf should be",i4)')ndf
    write(7,'("kbcflag is wrong, please check input file!")')
    stop
  endif
  k=kinp ! for keeping kinp's value not changing
  kout(1)=k/10**(ndf-1)
  do idf=2,ndf
    k=k-kout(idf-1)*10**(ndf-idf+1)
    kout(idf)=k/10**(ndf-idf)
  enddo
end subroutine bc_mark

subroutine read_body(nbody,nte,ibody,mbody,kbdelem,vbody)
  ! read body force or heat source
  implicit real*8 (a-h,o-z)
  dimension ibody(nbody),mbody(nbody),kbdelem(nte),vbody(10,nbody)
  
  if(nbody.eq.0) then
    write(7,'(/4x,"No body force information")')
  else
    write(7,'(/4x,"Body force groups =",i5)') nbody
    do ibo=1,nbody
      read(5,*)m,mbody(m),ibody(m),(vbody(j,m),j=1,ibody(m))
      write(7,'(3i9,1p10e15.6)')m,mbody(m),ibody(m),(vbody(j,m),j=1,ibody(m))
    enddo
    kbdelem=0
    read(5,*)nbdelem
    if(nbdelem.eq.0) then ! applied to all elements without more input
      write(7,'(/4x,"Body force is applied to all element")')
      do ibde=1,nte
        kbdelem(ibde)=1
      enddo
    else ! specify elements through multiple groups
      do ibdelem=1,nbdelem
        write(7,'(/4x,"Body force element groups =",i5)')nbdelem
        read(5,*) nebgn,neend,nestep,ibdgrp
        write(7,*) nebgn,neend,nestep,ibdgrp
        do ie=nebgn,neend,nestep
          kbdelem(ie)=ibdgrp
        enddo
      enddo
    endif
  endif
  
end subroutine read_body

subroutine read_material
  ! read material from input file if it is a linear problem
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  dimension dmat(nsigp,nsigp),c1(21),c2(21),alphab(ndf)
  
  if(nmatgrp.gt.0) then ! material properties are specified by elements' nodes
    allocate(pred(nmatgrp,mnode,nsigp,nsigp),lrem(nmatgrp))
    if(istrans.eq.1)then
      allocate(prerho(nmatgrp,mnode),precp(nmatgrp,mnode))
    endif
    if(istm.eq.1)allocate(prect0(nmatgrp,mnode,ndf))
    do imatgrp=1,nmatgrp
      read(5,*)imat
      if(abs(imat).gt.nmatgrp) then
        write(*,*)'The number of a material group exceed nmatgrp.'
        write(7,*)'The number of a material group exceed nmatgrp.'
        stop
      endif
      if(imat.gt.0) then
        backspace(5)
        if(istrans.eq.0.and.istm.eq.0)then ! steady state without thermal stress
          read(5,*)imat,e0,pr,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
          write(7,71)imat,e0,pr,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
        elseif(istrans.eq.0.and.istm.eq.1)then ! steady state with thermal stress
          read(5,*)imat,e0,pr,ct0,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
          write(7,72)imat,e0,pr,ct0,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
        elseif(istrans.eq.1.and.istm.eq.0)then ! transient or dynamic without thermal stress
          read(5,*)imat,e0,pr,rho,cp,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2) ! cp is specific heat in heat conduction and damp in mechanic
          write(7,73)imat,e0,pr,rho,cp,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
        elseif(istrans.eq.1.and.istm.eq.1)then ! dynamic with thermal stress
          read(5,*)imat,e0,pr,ct0,rho,cp,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
          write(7,74)imat,e0,pr,ct0,rho,cp,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
        endif
        if(mtype.ne.0) lrem(imat)=1 ! non-homogeneous material
        do id=1,mnode
          kp=1
          if(mtype.ne.0) then ! non-homogeneous material
            do ie=1,nte
              if(nematg(ie).eq.imat) then
                kp=lnde(id,ie)
                exit
              endif
            enddo
          endif
          e=val_mu(ndim,e0,n1,c1,cd(:,kp),mtype)
          call dmatrix(e,pr,ct0,alphab,misot,n2,c2,dmat)
          do i=1,nsigp
            do j=1,nsigp
              pred(imat,id,i,j)=dmat(i,j)
            enddo
          enddo
          if(istrans.eq.1)then
            prerho(imat,id)=rho
            precp(imat,id)=cp
          endif
          if(istm.eq.1)prect0(imat,id,:)=alphab(:)
        enddo
      else ! imat.lt.0, the material information of all nodes in the element are needed
        lrem(-imat)=1 ! must be non-homogeneous material
        do id=1,mnode
          if(istrans.eq.0.and.istm.eq.0)then ! steady state without thermal stress
            read(5,*)e0,pr,misot,n2,(c2(m),m=1,n2)
            write(7,75)e0,pr,misot,n2,(c2(m),m=1,n2)
          elseif(istrans.eq.0.and.istm.eq.1)then ! steady state with thermal stress
            read(5,*)e0,pr,ct0,misot,n2,(c2(m),m=1,n2)
            write(7,76)e0,pr,ct0,misot,n2,(c2(m),m=1,n2)
          elseif(istrans.eq.1.and.istm.eq.0)then ! transient or dynamic without thermal stress
            read(5,*)e0,pr,rho,cp,misot,n2,(c2(m),m=1,n2) ! cp is specific heat in heat conduction and damp in mechanic
            write(7,77)e0,pr,rho,cp,misot,n2,(c2(m),m=1,n2)
          elseif(istrans.eq.1.and.istm.eq.1)then ! dynamic with thermal stress
            read(5,*)e0,pr,ct0,rho,cp,misot,n2,(c2(m),m=1,n2)
            write(7,78)e0,pr,ct0,rho,cp,misot,n2,(c2(m),m=1,n2)
          endif
          call dmatrix(e0,pr,ct0,alphab,misot,n2,c2,dmat)
          do i=1,nsigp
            do j=1,nsigp
              pred(-imat,id,i,j)=dmat(i,j)
            enddo
          enddo
          if(istrans.eq.1)then
            prerho(-imat,id)=rho
            precp(-imat,id)=cp
          endif
          if(istm.eq.1)prect0(-imat,id,:)=alphab(:)
        enddo
      endif
    enddo
  else ! material properties are specified by global nodes
    allocate(pred(ntp,1,nsigp,nsigp))
    if(istrans.eq.1)then
      allocate(prerho(ntp,1),precp(ntp,1))
    endif
    if(istm.eq.1)allocate(prect0(ntp,1,ndf))
    do imatgrp=1,-nmatgrp
      if(istrans.eq.0.and.istm.eq.0)then ! steady state without thermal stress
        read(5,*)mpbgn,mpend,mpstp,e0,pr,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
        write(7,79)mpbgn,mpend,mpstp,e0,pr,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
      elseif(istrans.eq.0.and.istm.eq.1)then ! steady state with thermal stress
        read(5,*)mpbgn,mpend,mpstp,e0,pr,ct0,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
        write(7,80)mpbgn,mpend,mpstp,e0,pr,ct0,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
      elseif(istrans.eq.1.and.istm.eq.0)then ! transient or dynamic without thermal stress
        read(5,*)mpbgn,mpend,mpstp,e0,pr,rho,cp,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2) ! cp is specific heat in heat conduction and damp in mechanic
        write(7,81)mpbgn,mpend,mpstp,e0,pr,rho,cp,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
      elseif(istrans.eq.1.and.istm.eq.1)then ! dynamic with thermal stress
        read(5,*)mpbgn,mpend,mpstp,e0,pr,ct0,rho,cp,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
        write(7,82)mpbgn,mpend,mpstp,e0,pr,ct0,rho,cp,mtype,n1,misot,n2,(c1(m),m=1,n1),(c2(m),m=1,n2)
      endif
      if(mtype.eq.0)then
        write(*,*)'Homogeneous materials cannot be specified by ','global nodes'
        write(7,*)'Homogeneous materials cannot be specified by ','global nodes'
        stop
      endif
      do ip=mpbgn,mpend,mpstp
        e=val_mu(ndim,e0,n1,c1,cd(:,ip),mtype)
        call dmatrix(e,pr,ct0,alphab,misot,n2,c2,dmat)
        do i=1,nsigp
          do j=1,nsigp
            pred(ip,1,i,j)=dmat(i,j)
          enddo
        enddo
        if(istrans.eq.1)then
          prerho(ip,1)=rho
          precp(ip,1)=cp
        endif
        if(istm.eq.1)prect0(ip,1,:)=alphab(:)
      enddo
    enddo
  endif
  avd=maxval(pred(:,:,1,1)) ! max value of the material parameter
  
  71 format(i8,1p2e15.6,4i4,1p24e15.6)
  72 format(i8,1p3e15.6,4i4,1p24e15.6)
  73 format(i8,1p4e15.6,4i4,1p24e15.6)
  74 format(i8,1p5e15.6,4i4,1p24e15.6)
  75 format(1p2e15.6,2i4,1p21e15.6)
  76 format(1p3e15.6,2i4,1p21e15.6)
  77 format(1p4e15.6,2i4,1p21e15.6)
  78 format(1p5e15.6,2i4,1p21e15.6)
  79 format(3i9,1p2e15.6,4i4,1p24e15.6)
  80 format(3i9,1p3e15.6,4i4,1p24e15.6)
  81 format(3i9,1p4e15.6,4i4,1p24e15.6)
  82 format(3i9,1p5e15.6,4i4,1p24e15.6)
  
end subroutine read_material

subroutine read_nonli_material
  ! read material from input file if it is a non-linear problem
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  
  allocate(metype(nmatgrp),ne0(nmatgrp),ce(5,nmatgrp))
  if(ndf.gt.1)then
    allocate(mptype(nmatgrp),npr(nmatgrp),cpr(5,nmatgrp))
    if(istm.eq.1)allocate(mcttype(nmatgrp),nct(nmatgrp),cct(5,nmatgrp))
  endif
  if(istrans.eq.1)then
    allocate(mrtype(nmatgrp),nrh(nmatgrp),crh(5,nmatgrp))
    allocate(mcptype(nmatgrp),ncp(nmatgrp),ccp(5,nmatgrp))
  endif
  do imatgrp=1,nmatgrp
    read(5,*) imat
    write(7,'(a7,i6)')'imat: ',imat
    read(5,*) metype(imat),ne0(imat),(ce(m,imat),m=1,ne0(imat)) ! elasticity modulus or heat conductivity
    write(7,10)'   E:',metype(imat),ne0(imat),(ce(m,imat),m=1,ne0(imat))
    if(ndf.gt.1)then
      read(5,*) mptype(imat),npr(imat),(cpr(m,imat),m=1,npr(imat)) ! poisson ratio
      write(7,10)' PR:',mptype(imat),npr(imat),(cpr(m,imat),m=1,npr(imat))
      if(istm.eq.1)then
        read(5,*)mcttype(imat),nct(imat),(cct(m,imat),m=1,nct(imat)) ! thermal expansion coefficient
        write(7,10)'CT0:',mcttype(imat),nct(imat),(cct(m,imat),m=1,nct(imat))
      endif
    endif
    if(istrans.eq.1)then
      read(5,*) mrtype(imat),nrh(imat),(crh(m,imat),m=1,nrh(imat)) ! density
      write(7,10)'RHO:',mrtype(imat),nrh(imat),(crh(m,imat),m=1,nrh(imat))
      read(5,*) mcptype(imat),ncp(imat),(ccp(m,imat),m=1,ncp(imat)) ! cp is specific heat in heat conduction and damp in mechanic
      write(7,10)' Cp:',mcptype(imat),ncp(imat),(ccp(m,imat),m=1,ncp(imat))
    endif
  enddo
  avd=maxval(ce(1,:)) ! max value of the material parameter
      
  10 format(a6,2i4,1p5e15.6)
end subroutine read_nonli_material

subroutine dmatrix(e,pr,ct0,alphab,misot,n2,c2,dmat)
  ! this routine defines the elastic stiffness matrix d or heat conductivity kij
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  dimension c2(n2),dmat(nsigp,nsigp),alphab(ndf)
  
  dmat=0.d0
  if(misot.eq.0) then ! isotropic material
    if(ndf.eq.1) then ! heat conduction problems
      g=e ! heat conductivity
    else ! mechanics problems
      g=e*0.5d0/(1.d0+pr) ! shear modulus
      vlamda=2.d0*pr*g/(1.d0-2.d0*pr) ! lame coefficient
      alphab=ct0*e/(1.d0-2.d0*pr) ! eq.(3) of thermal stress paper in 2003  
      if(nsig.eq.3.and.ndf.eq.2) then ! plane stress problem, shear modulus donot be changed
        !e1=e*(2.d0*pr+1.d0)/(pr+1.d0)**2
        pr1=pr/(1.d0+pr)
        vlamda=2.d0*pr1*g/(1.d0-2.d0*pr1)
        alphab=ct0*e/(1.d0-pr) ! or alphab=ct0*2.d0*g/(1.d0-2.d0*pr1)
      endif
    endif
    if(ndf.eq.1) then ! potential problems
      do i=1,ndim
        dmat(i,i)=g
      enddo
    else ! elasticity problems
      do i=1,ndim
        do j=1,ndim
          if(i.eq.j)then
            dmat(i,j)=vlamda+2.d0*g
          else
            dmat(i,j)=vlamda
          endif
        enddo
      enddo
      do i=ndim+1,nsigp
        dmat(i,i)=g
      enddo
    endif
  else ! general anisotropic materials
    if(ndf.eq.1)then ! heat conduction problems
      write(*,*)'Anisotropic materials is not for heat conduction'
      write(7,*)'Anisotropic materials is not for heat conduction'
      stop
    elseif(ndf.eq.2)then ! plane stress anisotropic problems
      write(*,*)'Anisotropic materials is not for 2d problems'
      write(7,*)'Anisotropic materials is not for 2d problems'
      stop
    else ! 3d elasticity problems
      e1=c2(1)
      e2=c2(2)
      e3=c2(3)
      v12=c2(4)
      v13=c2(5)
      v23=c2(6)
      g12=c2(7)
      g13=c2(8)
      g23=c2(9)
      v21=v12*e2/e1
      v31=v13*e3/e1
      v32=v23*e3/e2
      delta=(1-v21*v32*v13-v31*v12*v23-v31*v13-v21*v12-v32*v23)/(e1*e2*e3)
      dmat(1,1)=(1-v23*v32)/(e2*e3*delta)
      dmat(2,2)=(1-v13*v31)/(e1*e3*delta)
      dmat(3,3)=(1-v12*v21)/(e1*e2*delta)
      dmat(1,2)=(v21+v31*v23)/(e2*e3*delta)
      dmat(2,1)=dmat(1,2)
      dmat(1,3)=(v13+v12*v23)/(e1*e2*delta)
      dmat(3,1)=dmat(1,3)
      dmat(2,3)=(v32+v12*v31)/(e1*e3*delta)
      dmat(3,2)=dmat(2,3)
      dmat(4,4)=g12
      dmat(5,5)=g23
      dmat(6,6)=g13
      alphab(1)=(dmat(1,1)+dmat(1,2)+dmat(1,3))*ct0
      alphab(2)=(dmat(2,1)+dmat(2,2)+dmat(2,3))*ct0
      alphab(3)=(dmat(3,1)+dmat(3,2)+dmat(3,3))*ct0
    endif
  endif
  
  return
end subroutine dmatrix

function val_mu(ndim,cmu0,n1,cmu,x,mtype)
  ! this routine determines mu or k for elasticity or heat conduction
  implicit real*8 (a-h,o-z)
  dimension cmu(n1),x(ndim)
  
  select case(mtype)
  case(0)
    val_mu=cmu0
  case(1)
    val_mu=cmu0*dexp(dot_product(cmu(1:ndim),x))
  case(2)
    val_mu=cmu0+dot_product(cmu(1:ndim),x)
  case(3)
    val_mu=cmu0+dot_product(cmu(1:ndim),x)
    n=ndim
    do i=1,ndim
      do j=i,ndim
          n=n+1
          val_mu=val_mu+cmu(n)*x(i)*x(j)
      enddo
    enddo
  case(4)
    r2=0.d0
    ix0=1+ndim
    do i=1,ndim
      r2=r2+cmu(1+i)*(x(i)-cmu(ix0+i))**2
    enddo
    val_mu=cmu0+cmu(1)*dsqrt(r2)
  case(5)
    r2=0.d0
    ix0=1+ndim
    do i=1,ndim
      r2=r2+cmu(1+i)*(x(i)-cmu(ix0+i))**2
    enddo
    val_mu=cmu0+cmu(1)*r2
  case(6)
    r2=0.d0
    ix0=1+ndim
    do i=1,ndim
      r2=r2+cmu(1+i)*(x(i)-cmu(ix0+i))**2
    enddo
    val_mu=cmu0*dexp(cmu(1)*dsqrt(r2))
  case(7)
    r2=0.d0
    ix0=1+ndim
    do i=1,ndim
      r2=r2+cmu(1+i)*(x(i)-cmu(ix0+i))**2
    enddo
    val_mu=cmu0*dexp(cmu(1)*r2)
  case(8)
    val_mu=(1.+x(1)/100.)**3 ! kassab for thermal
  case(9)
    val_mu=cmu0*dexp(-x(1)/6.d0)
  case(10)
    beta=1.d0/15.d0*cmu(1)
    val_mu=cmu0*dexp(beta*x(2))
  end select
  
end function val_mu