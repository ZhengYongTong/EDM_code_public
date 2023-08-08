!****************************************************************************
  ! This file containing:
  ! Subroutine coef_gen: generate the coefficient matrix
  ! Subroutine coefb_gen: transfer the derivate of the shape function from the vector to a matrix
  ! Subroutine treat_dispcond: set the coefficients in coefu into the whole coefficient matrix (find the right place of the vector a)
  ! Subroutine treat_dispcond_tran: set the coefficients in coefu into the whole coefficient matrix (find the right place of the vector a, c and rm) for transient problem
  ! Subroutine face_normal: calculate the normal of the face (specified by ie and iface) at source point (specified by inode)
  ! Subroutine first_type_bc: apply the first type of B.C.
  ! Subroutine second_type_bc: apply the second type of B.C.
  ! Subroutine third_type_bc: apply the third type of B.C. (only for heat convection B.C.)
  ! Subroutine press_bc: deal with press B.C. (related to outer face normal)
  ! Subroutine source_item: add the source term
  ! Subroutine val_body: define the value of source term
!****************************************************************************

subroutine coef_gen
  use cmd_progress
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  
  type( cls_cmd_progress ) ::progress ! for progress bar
  dimension xi(3),elcod(3,mnode),dndx(3,mnode),dndxx(3,3,mnode),lndb(25),ts(ndf),ck(3,25),coefu(ndf,ndf),bval(nsigp,ndf)
  dimension cosvalsum(nsigp,ndf),cosval(nsigp,ndf,nfacem),bvalf(nsigp,ndf),cosb(3),dval(nsigp,nsigp,mnode),db(nsigp,ndf)
  dimension coefut(ndf,ndf),ct0(ndf,mnode),rho(mnode),cp(mnode)
  allocatable tval(:)
  
  ! read nodal temperatures
  if(ndf.ne.1.and.istm.eq.1) then
    rewind (45)
    read(45,*) ktp
    allocate (tval(ntp))
    do ik=1,min(ktp,ntp)
     read(45,*)ip,tval(ip)
    enddo
  endif
  
  call progress % set( n = nte , l = 25 ) ! for progress bar
  do ie=1,nte ! element loop
    call progress % put( ie , cmd_progress_absolute ) ! for progress bar
    
    ietype=netype(ie)
    node=nnode(ie)
    ndk=ndtok(node)
    
    ! get material of each node in the element
    if(nmatgrp.lt.0)then ! material properties are specified by global nodes
      isnonho=1 ! must be non-homogeneous materials
      do knode=1,node
        elcod(1:ndim,knode)=cd(1:ndim,lnde(knode,ie))
        ip=lnde(knode,ie)
        dval(:,:,knode)=pred(ip,1,:,:)
        if(istm.eq.1) ct0(:,knode)=prect0(ip,1,:)
        if(istrans.eq.1)then
          rho(knode)=prerho(ip,1)
          cp(knode)=precp(ip,1)
        endif
      enddo
    else ! material properties are specified by elements
      if(lrem(nematg(ie)).eq.1) then ! non-homogeneous materials
        isnonho=1
      else ! homogeneous materials
        isnonho=0
      endif
      do knode=1,node
        elcod(1:ndim,knode)=cd(1:ndim,lnde(knode,ie))
        dval(:,:,knode)=pred(nematg(ie),knode,:,:)
        if(istm.eq.1) ct0(:,knode)=prect0(nematg(ie),knode,:)
        if(istrans.eq.1)then
          rho(knode)=prerho(nematg(ie),knode)
          cp(knode)=precp(nematg(ie),knode)
        endif
      enddo
    endif
    
    do inode=1,node ! source point loop
      igp=lnde(inode,ie)
      ip0=ndf*(igp-1)
      call node_xi(inode,ietype,xi) ! instrict coordinates at source point
      ncl=neoe(igp)*ndf
      idtype=ndtype(inode,ndk)
      
      select case(idtype)
      case(1) ! element internal nodes
       dndx(:,1:node)=dndxa(:,1:node,inode,ie)
       dndxx(:,:,1:node)=dndxxa(:,:,1:node,inode,ie)
       
       do jnode=1,node ! field point loop
        jgp=lnde(jnode,ie)
        jp0=ndf*(jgp-1)
        if(ndf.eq.1)then ! heat conduction problem
          coefu(1,1)=0.d0
          do i=1,ndim
            coefu(1,1)=coefu(1,1)+dval(i,i,inode)*dndxx(i,i,jnode)
          enddo
        else ! mechanics problem
          if(ndim.eq.1)then
            coefu(1,1)=dval(1,1,inode)*dndxx(1,1,jnode)
          elseif(ndim.eq.2)then
            coefu(1,1)=dval(1,1,inode)*dndxx(1,1,jnode)+dval(3,3,inode)*dndxx(2,2,jnode)
            coefu(1,2)=dval(1,2,inode)*dndxx(1,2,jnode)+dval(3,3,inode)*dndxx(1,2,jnode)
            coefu(2,1)=dval(2,1,inode)*dndxx(1,2,jnode)+dval(3,3,inode)*dndxx(1,2,jnode)
            coefu(2,2)=dval(2,2,inode)*dndxx(2,2,jnode)+dval(3,3,inode)*dndxx(1,1,jnode)
          else ! ndim.eq.3
            coefu(1,1)=dval(1,1,inode)*dndxx(1,1,jnode)+dval(4,4,inode)*dndxx(2,2,jnode)+dval(6,6,inode)*dndxx(3,3,jnode)
            coefu(1,2)=dval(1,2,inode)*dndxx(1,2,jnode)+dval(4,4,inode)*dndxx(1,2,jnode)
            coefu(1,3)=dval(1,3,inode)*dndxx(1,3,jnode)+dval(6,6,inode)*dndxx(1,3,jnode)
            coefu(2,1)=dval(2,1,inode)*dndxx(1,2,jnode)+dval(4,4,inode)*dndxx(1,2,jnode)
            coefu(2,2)=dval(2,2,inode)*dndxx(2,2,jnode)+dval(4,4,inode)*dndxx(1,1,jnode)+dval(5,5,inode)*dndxx(3,3,jnode)
            coefu(2,3)=dval(2,3,inode)*dndxx(2,3,jnode)+dval(5,5,inode)*dndxx(2,3,jnode)
            coefu(3,1)=dval(3,1,inode)*dndxx(1,3,jnode)+dval(6,6,inode)*dndxx(1,3,jnode)
            coefu(3,2)=dval(3,2,inode)*dndxx(2,3,jnode)+dval(5,5,inode)*dndxx(2,3,jnode)
            coefu(3,3)=dval(3,3,inode)*dndxx(3,3,jnode)+dval(6,6,inode)*dndxx(1,1,jnode)+dval(5,5,inode)*dndxx(2,2,jnode)
          endif
        endif
        
        if(isnonho.eq.1) then ! non-homogeneous or variable coefficient problems
          ! calculate the derivate matrix of the shape function with jnode's (field point) local coordinates
          bvalf=0.d0
          call coefb_gen(bvalf,dndx(:,jnode))
          do knode=1,node
            ! calculate the derivate matrix of the shape function with knode's local coordinates
            bval=0.d0
            call coefb_gen(bval,dndx(:,knode))
            call matmulmkl(nsigp,nsigp,ndf,dval(:,:,knode),bvalf,db,1.d0,0.d0)
            call matmulmkl(ndf,nsigp,ndf,transpose(bval),db,coefut,1.d0,0.d0)
            coefu=coefu+coefut
            !coefu=coefu+matmul(matmul(transpose(bval),dval(:,:,knode)),bvalf)
          enddo
        endif
        if(istrans.eq.1)then ! transient problem
          call treat_dispcond_tran(ncl,ip0,jp0,coefu,rho(inode),cp(inode),idtype)
        else
          call treat_dispcond(ncl,ip0,jp0,coefu)
        endif
       enddo ! end of field point loop
      
       ! add thermal stress term
       if(ndf.ne.1.and.istm.eq.1) then
         ts=0.d0
         do jd=1,node
           ts=ts+dndx(1:ndf,jd)*ct0(1:ndf,jd)*tval(lnde(jd,ie))
         enddo
         b(ip0+1:ip0+ndf)=b(ip0+1:ip0+ndf)+ts
       endif
       
      case(2) ! element boundary nodes of the problem
       
       cosval=0.d0
       cosvalsum=0.d0
       do iface=1,nface(ndk)
         call face_normal(node,ie,inode,iface,cosb,kpface,lndb,ietype,ck,ndk,info) ! kpface is source point in the face
         if(info.eq.1)then
           write(*,'(a18,i8,i3)')'Zero face jacobian, ie=',ie,' ,iface=',iface
           write(4,'(a18,i8,i3)')'Zero face jacobian, ie=',ie,' ,iface=',iface
           write(7,'(a18,i8,i3)')'Zero face jacobian, ie=',ie,' ,iface=',iface
           stop
         endif
         if(kpface.eq.0) cycle ! inode is not one of the nodes of face iface
         ! calculate face nomal matrix of each face --cosval, and add them to get cosvalsum
         call coefb_gen(cosval(:,:,iface),cosb)
         cosvalsum(:,:)=cosvalsum(:,:)+cosval(:,:,iface)
         ! treat thermal loads
         if(ndf.ne.1.and.istm.eq.1) then
           b(ip0+1:ip0+ndf)=b(ip0+1:ip0+ndf)+ct0(1:ndf,inode)*tval(igp)*cosb(1:ndf)
         endif
       enddo
       
       dndx(:,1:node)=dndxa(:,1:node,inode,ie)
       nn=nrn(inode,ndk)
       do jn=1,nn ! field point loop
         jnode=ndrn(jn,inode,ndk)
         jgp=lnde(jnode,ie)
         jp0=ndf*(jgp-1)
         ! calculate the derivate matrix of the shape function -- bval
         bval=0.d0
         call coefb_gen(bval,dndx(:,jnode))
         ! calculate coefficient coefu
         call matmulmkl(nsigp,nsigp,ndf,dval(:,:,inode),bval,db,1.d0,0.d0)
         call matmulmkl(ndf,nsigp,ndf,transpose(cosvalsum),db,coefu,1.d0,0.d0)
         !coefu=matmul(matmul(transpose(cosvalsum),dval(:,:,inode)),bval)
         if(ndf.eq.1) coefu=-coefu ! in heat conduction, q=-kij*
         if(istrans.eq.1)then ! transient problem
           call treat_dispcond_tran(ncl,ip0,jp0,coefu,0.d0,0.d0,idtype)
         else
           call treat_dispcond(ncl,ip0,jp0,coefu)
         endif
       enddo ! end of field point loop
       
      end select
    enddo ! end of source point loop
  enddo ! end of element loop
  if(ndf.ne.1.and.istm.eq.1) deallocate (tval)
  
  return
end subroutine coef_gen

subroutine coefb_gen(bval,dndx)
  ! transfer the derivate of the shape function from the vector to a matrix
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  dimension bval(nsigp,ndf),dndx(3)
  if(ndf.eq.1)then ! heat conduction problem
    bval(1:ndim,1)=dndx(1:ndim)
  else ! mechanics problem
    do idf=1,ndf
      bval(idf,idf)=dndx(idf)
    enddo
    if(ndim.eq.2.and.ndf.eq.2)then ! 2d mechanics problem
      bval(3,1)=dndx(2)
      bval(3,2)=dndx(1)
    elseif(ndim.eq.3.and.ndf.eq.3)then ! 3d mechanics problem
      bval(4,1)=dndx(2)
      bval(4,2)=dndx(1)
      bval(5,2)=dndx(3)
      bval(5,3)=dndx(2)
      bval(6,1)=dndx(3)
      bval(6,3)=dndx(1)
    endif
  endif
end subroutine coefb_gen

subroutine treat_dispcond(ncl,ip0,jp0,coefu)
  ! set the coefficients in coefu into the whole coefficient matrix (find the right place of the vector a)
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  dimension coefu(ndf,ndf)

  do j=ia(ip0+1),ia(ip0+2)-1,ndf
    if(ja(j).eq.jp0+1)exit
  enddo
  do idf=1,ndf
    do jdf=1,ndf
      a(j+(idf-1)*ncl+jdf-1)=a(j+(idf-1)*ncl+jdf-1)+coefu(idf,jdf)
    enddo
  enddo
  
end subroutine treat_dispcond

subroutine treat_dispcond_tran(ncl,ip0,jp0,coefu,rho,cp,idtype)
  ! set the coefficients in coefu into the whole coefficient matrix (find the right place of the vector a, c and rm) for transient problem
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  dimension coefu(ndf,ndf)
  
  do j=ia(ip0+1),ia(ip0+2)-1,ndf
    if(ja(j).eq.jp0+1)exit
  enddo
  do idf=1,ndf
    do jdf=1,ndf
      a(j+(idf-1)*ncl+jdf-1)=a(j+(idf-1)*ncl+jdf-1)+coefu(idf,jdf)
    enddo
  enddo
  if(ip0.eq.jp0.and.idtype.eq.1)then ! the diagnal of this line
    if(ndf.eq.1)then ! transient heat conduction
      c(ic(ip0+1))=c(ic(ip0+1))+rho*cp
    else ! dynamics
      ipf=lnof(ip0/ndf+1)
      ipi0=(ipf-nobdf-1)*ndf
      do idf=1,ndf
        rm(ipi0+idf)=rm(ipi0+idf)-rho
        c(ipi0+idf)=c(ipi0+idf)-cp
      enddo
    endif
  endif
  
end subroutine treat_dispcond_tran

subroutine face_normal(node,ie,inode,iface,cosb,kpface,lndb,ietype,ck,ndk,info)
  ! calculate the normal of the face (specified by ie and iface) at source point (specified by inode)
  use vary_arrays
  use fixed_values
  implicit real*8 (a-h,o-z)
  dimension dn(3,25),gd(3,3),xib(3),ck(3,npface(ndk)),cosb(3),lndb(npface(ndk))
  kpface=0
  do ipface=1,npface(ndk)
    if(ndim.eq.1) then
      if(iface.eq.1) knode=1
      if(iface.eq.2) knode=node
    else
      knode=ndface(ipface,iface,ndk)
    endif
    if(knode.eq.inode) kpface=ipface
    lndb(ipface)=lnde(knode,ie)
  enddo
  if(kpface.eq.0) return
  do ipface=1,npface(ndk)
    ck(1:ndim,ipface)=cd(1:ndim,lndb(ipface))
  enddo
  if(ndim.eq.1) then
    if(iface.eq.1) then
      term=cd(1,lnde(1,ie))-cd(1,lnde(node,ie))
    else
      term=cd(1,lnde(node,ie))-cd(1,lnde(1,ie))
    endif
    cosb(1)=term/dabs(term)
    return
  elseif(ndim.eq.2) then
    ketype=npface(ndk)
  else
    ketype=300+npface(ndk)
  endif
  call node_xi(kpface,ketype,xib)
  info=0
  call dshapef(ndim,ndim-1,ketype,npface(ndk),xib,ck,cosb,fjcb,info)
end subroutine face_normal

subroutine first_type_bc(istep) ! apply the first type of B.C.
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  
  do ie=1,nte ! element loop
    node=nnode(ie)
    ndk=ndtok(node)
    do iface=1,nface(ndk) ! face loop
      igrp=nbcgrp(ie,iface)
      if(igrp.eq.0)cycle
      if(lret(igrp).ne.1)cycle ! not first type of B.C.
      ! dynamics simple harmonic B.C.
      if(itef(igrp).eq.1)then
        acoef=tef(igrp)
        acoef=dsin(acoef*deltat*dble(istep))
      elseif(itef(igrp).eq.2)then
        acoef=tef(igrp)
        acoef=dcos(acoef*deltat*dble(istep))
      else
        acoef=1.d0
      endif
      do ipface=1,npface(ndk) ! point loop
        inode=ndface(ipface,iface,ndk)
        igp=lnde(inode,ie)
        ip0=(igp-1)*ndf
        ncl=neoe(igp)*ndf
        do j=ia(ip0+1),ia(ip0+2)-1,ndf
          if(ja(j).eq.ip0+1)exit
        enddo
        jdiag=j ! find the of diagonal position of (ip0+1)-th line
        ! get the value of first type B.C. and add the first type B.C. into system
        do idf=1,ndf
          if(kbcflag(idf,ie,iface).ne.1) cycle
          ipdf=jdiag+ncl*(idf-1)+idf-1
          a(ipdf)=a(ipdf)+1.d15*avd*avl
          b(ip0+idf)=pret(igrp,idf)*a(ipdf)*acoef
          isfirtyp(ip0+idf)=1 ! record the row
        enddo
      enddo ! end of point loop
    enddo ! end of face loop
  enddo ! end of element loop
  
end subroutine first_type_bc

subroutine second_type_bc(istep) ! apply the second type of B.C.
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  dimension ck(3,25),vtgrp(ndf)
  
  do ie=1,nte ! element loop
    node=nnode(ie)
    ndk=ndtok(node)
    do iface=1,nface(ndk) ! face loop
      igrp=nbcgrp(ie,iface)
      if(igrp.eq.0)cycle
      if(lret(igrp).ne.2.and.lret(igrp).ne.-2)cycle ! not second type of B.C.
      ! dynamics simple harmonic B.C.
      if(itef(igrp).eq.1)then
        acoef=tef(igrp)
        acoef=dsin(acoef*deltat*dble(istep))
      elseif(itef(igrp).eq.2)then
        acoef=tef(igrp)
        acoef=dcos(acoef*deltat*dble(istep))
      else
        acoef=1.d0
      endif
      ! press B.C. need the global coordinates of the nodes in the face
      if(ndf.gt.1.and.lret(igrp).eq.-2) then
        do ipface=1,npface(ndk)
          inode=ndface(ipface,iface,ndk)
          ck(1:ndim,ipface)=cd(1:ndim,lnde(inode,ie))
        enddo
      endif
      do ipface=1,npface(ndk) ! point loop
        inode=ndface(ipface,iface,ndk)
        igp=lnde(inode,ie)
        ip0=(igp-1)*ndf
        ! get the value of second type B.C. and add the second type B.C. into system
        if(ndf.gt.1)then ! mechanics problem
          if(lret(igrp).eq.2) then ! normal b.c.
            vtgrp=pret(igrp,1:ndf)
          elseif(lret(igrp).eq.-2) then ! press is specified
            press=pret(igrp,1)
            call press_bc(ipface,ndim,ndf,npface(ndk),ck,press,vtgrp,info) ! treat press b.c.
            if(info.eq.1)then
              write(*,'(a18,i8,i3)')'zero face jacobian',ie,iface
              write(4,'(a18,i8,i3)')'zero face jacobian',ie,iface
              write(7,'(a18,i8,i3)')'zero face jacobian',ie,iface
              stop
            endif
          endif
          do idf=1,ndf
            if(kbcflag(idf,ie,iface).ne.2.or.isfirtyp(ip0+idf).eq.1) cycle ! jump the directions whose tractions are unknown or displacements are known
            b(ip0+idf)=b(ip0+idf)+vtgrp(idf)*acoef
          enddo
        else ! heat conduction problem
          if(kbcflag(1,ie,iface).eq.2.and.isfirtyp(ip0+1).eq.0) b(ip0+1)=b(ip0+1)+pret(igrp,1)*acoef
        endif
      enddo ! end of point loop
    enddo ! end of face loop
  enddo ! end of element loop
  
end subroutine second_type_bc

subroutine third_type_bc(istep) ! apply the third type of B.C. (only for heat convection B.C.)
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  
  do ie=1,nte ! element loop
    node=nnode(ie)
    ndk=ndtok(node)
    do iface=1,nface(ndk) ! face loop
      igrp=nbcgrp(ie,iface)
      if(igrp.eq.0)cycle
      if(lret(igrp).ne.3)cycle ! not third type of B.C.
      do ipface=1,npface(ndk) ! point loop
        inode=ndface(ipface,iface,ndk)
        igp=lnde(inode,ie)
        ip0=(igp-1)*ndf
        do j=ia(ip0+1),ia(ip0+2)-1,ndf
          if(ja(j).eq.ip0+1)exit
        enddo
        jdiag=j ! find the of diagonal position of (ip0+1)-th line
        ! get the value of third type B.C. and add the third type B.C. into system
        if(ndf.gt.1)then ! mechanics problem
          ! to be developed
        else ! heat conduction problem
          coef_h=pret(igrp,1)
          tenv=pret(igrp,2)
          a(jdiag)=a(jdiag)-coef_h ! coef_h is the heat transfer coefficient h
          b(ip0+1)=b(ip0+1)-coef_h*tenv ! tenv is the enviroment temperature
        endif
      enddo ! end of point loop
    enddo ! end of face loop
  enddo ! end of element loop
  
end subroutine third_type_bc

subroutine press_bc(id,ndim,ndf,npface,ck,press,vtgrp,info)
  implicit real*8 (a-h,o-z)
  dimension cosn(3),ck(3,npface)
  dimension dn(3,25),gd(3,3),xib(3),vtgrp(ndf)
  if(ndim.eq.2) then
    ketype=npface
  else ! ndim.eq.3 
    ketype=300+npface
  endif
  call node_xi(id,ketype,xib)
  info=0
  call dshapef(ndim,ndim-1,ketype,npface,xib,ck,cosn,fjcb,info)
  vtgrp=-press*cosn(1:ndim) ! press is positive and pull is negative
end subroutine press_bc

subroutine source_item(istep) ! add the source term
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  dimension body(ndf)
  
  if(nbody.eq.0)return
  do ie=1,nte ! element loop
    node=nnode(ie)
    ndk=ndtok(node)
    do inode=1,node ! point loop
      ip=lnde(inode,ie)
      ip0=ndf*(ip-1)
      idtype=ndtype(inode,ndk)
      if(idtype.gt.1) cycle
      body=0.d0
      m=kbdelem(ie) ! body force group
      kbde=mbody(m) ! body force type
      call val_body(ndim,ndf,nbody,cd(1:ndim,ip),vbody,m,kbde,body)
      b(ip0+1:ip0+ndf)=b(ip0+1:ip0+ndf)-body
    enddo ! point loop
  enddo ! element loop
  
end subroutine source_item

subroutine val_body(ndim,ndf,nbody,x,vbody,m,kbde,bodyv) ! define the value of source term
  implicit real*8 (a-h,o-z)
  dimension x(ndim),vbody(10,nbody),bodyv(ndf)
  
  select case(kbde)
  case(1)
    do idf=1,ndf
      bodyv(idf)=bodyv(idf)+vbody(idf,m)
    enddo
  case(2)
    do idf=1,ndf
      bodyv(idf)=bodyv(idf)+vbody(idf,m)+vbody(idf+ndf,m)*x(1)
    enddo
  case(3) ! kassab for thermal
    xnorm=x(1)/100.d0
    bodyv(1)=xnorm*(1.-xnorm)
  case(4)
    rxy=dsqrt(x(1)*x(1)+x(2)*x(2)) ! cylinder 3d for thermal
    if(rxy.gt.0.1d0.or.x(3).gt.0.3d0) then 
     bodyv(1)=0.d0
    else
     bodyv(1)=(0.1d0-rxy)*(0.3d0-x(3))/0.03d0/0.01d0
    endif
  case(5)
    bodyv(1)=vbody(1,m)*x(1)*(6.d0-x(1))
  case(6)
    bodyv(1)=vbody(1,m)+vbody(2,m)*dsqrt(x(1)**2+x(2)**2)
  end select
  
end subroutine val_body