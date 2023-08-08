!****************************************************************************
  ! This file containing:
  ! Subroutine coef_nonli: generate the coefficient matrix of material nonlinear problem
  ! Subroutine val_mu_nonli: calculate the value of material parameters
!****************************************************************************

subroutine coef_nonli(istep,tval)
  use fixed_values
  use vary_arrays
  implicit real*8 (a-h,o-z)
  
  dimension xi(3),elcod(3,mnode),dndx(3,mnode),spu(ndf),dndxx(3,3,mnode),vtgrp(ndf,mnode),lndb(25),ts(ndf),ck(3,25)
  dimension coefu(ndf,ndf),bval(nsigp,ndf),cosvalsum(nsigp,ndf),cosval(nsigp,ndf,nfacem),bvalf(nsigp,ndf),ct0(ndf,mnode)
  dimension cosb(3),dval(nsigp,nsigp,mnode),db(nsigp,ndf),coefut(ndf,ndf),tval(ntp),c2(9),rho(mnode),cp(mnode)
  
  avd=0.d0 ! max value of the material parameter
  do ie=1,nte ! element loop
  
    ietype=netype(ie)
    node=nnode(ie)
    ndk=ndtok(node)
    imat=nematg(ie)
    
    ! get material of each node in the element
    isnonho=1 ! must be non-homogeneous materials
    do knode=1,node
      elcod(1:ndim,knode)=cd(1:ndim,lnde(knode,ie))
      ip=lnde(knode,ie)
      call val_mu_nonli(ndf,ntf,xval,ip,metype(imat),ne0(imat),ce(1:ne0(imat),imat),e0)
      if(ndf.eq.1)then
        call dmatrix(e0,0.d0,0.d0,ct0(:,knode),0,1,c2(1),dval(:,:,knode))
        if(istrans.eq.1)then
          call val_mu_nonli(ndf,ntf,xval,ip,mrtype(imat),nrh(imat),crh(1:nrh(imat),imat),rho(knode))
          call val_mu_nonli(ndf,ntf,xval,ip,mcptype(imat),ncp(imat),ccp(1:ncp(imat),imat),cp(knode))
        endif
      else ! ndf.gt.1
        call val_mu_nonli(ndf,ntf,xval,ip,mptype(imat),npr(imat),cpr(1:npr(imat),imat),pr)
        if(istm.eq.0)then
          call dmatrix(e0,pr,0.d0,ct0(:,knode),0,1,c2(1),dval(:,:,knode))
        else ! istm.ne.0
          call val_mu_nonli(ndf,ntf,xval,ip,mcttype(imat),nct(imat),cct(1:nct(imat),imat),ct)
          call dmatrix(e0,pr,ct,ct0(:,knode),0,1,c2(1),dval(:,:,knode))
        endif
      endif
    enddo
    avd=max(avd,maxval(dval(1,1,:))) ! max value of the material parameter
    
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
          rhocp=rho(inode)*cp(inode)/deltat
          call treat_dispcond_tran(ncl,ip0,jp0,coefu,rho(inode),cp(inode),idtype)
        else
          call treat_dispcond(ncl,ip0,jp0,coefu)
        endif
       enddo
       
       ts=0.d0
       if(ndf.ne.1.and.istm.eq.1) then ! thermal stresses
         do jd=1,node
           ts=ts+dndx(1:ndf,jd)*ct0(1:ndf,jd)*tval(lnde(jd,ie))
         enddo
       endif
       b(ip0+1:ip0+ndf)=b(ip0+1:ip0+ndf)+ts
       
      case(2)  ! element boundary nodes of the problem
       
       cosval=0.d0
       cosvalsum=0.d0
       do iface=1,nface(ndk)
         call face_normal(node,ie,inode,iface,cosb,kpface,lndb,ietype,ck,ndk,info) ! kpface源点在面内的编号
         if(info.eq.1)then
           write(*,'(a18,i8,i3)')'zero face jacobian',ie,iface
           write(4,'(a18,i8,i3)')'zero face jacobian',ie,iface
           write(7,'(a18,i8,i3)')'zero face jacobian',ie,iface
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
       enddo ! end of face loop
         
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
  
  return
end subroutine coef_nonli
   
subroutine val_mu_nonli(ndf,ntf,xval,ip,mtype,nr,cr,vr)
  implicit real*8 (a-h,o-z)
  
  dimension xval(ntf),cr(nr)
  
  ip0=(ip-1)*ndf
  select case(mtype)
  case(0)
    vr=cr(1)
  case(1)
    vr=cr(1)+cr(2)*xval(ip0+1)
  case(2)
    vr=cr(1)+cr(2)*xval(ip0+1)+cr(3)*xval(ip0+1)**2.d0
  end select
  
end subroutine val_mu_nonli