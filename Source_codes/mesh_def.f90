!****************************************************************************
  ! This file containing:
  ! Subroutine node_sort: generate the matrix pattern of the coefficient matrix
  ! Subroutine gen_dndx_dndxx: ! Get the first and second order derivates of the shape functions with global coordinates
  ! Function my_int_compar: comparing function in qsort
!****************************************************************************

subroutine node_sort(maxe,maxn)
  use vary_arrays
  use fixed_values
  use, intrinsic :: iso_c_binding, only: c_size_t ! for involking sort function in C language
  use ifport ! for invoking sort function in C language
  implicit real*8 (a-h,o-z)

  integer(2),external::my_int_compar
  integer(c_size_t)::icsize,jcsize
  
  allocatable nonode(:,:),nele(:,:),kn0sysb(:),nple(:,:),lbon(:),lion(:)
  allocate(nonode(maxe,ntp),nele(maxn,ntp),kn0sysb(ntp),nple(maxe,ntp))
  allocate(neop(ntp),neoe(ntp))
  
  ! find the element related to each node
  neop=0 ! neop(i) is the number of elements that are related to the i-th node
  nple=0 ! nple(j,i) is the j-th element related to the i-th node
  do ie=1,nte
    node=nnode(ie)
    do id=1,node
      ip=lnde(id,ie)
      neop(ip)=neop(ip)+1
      nple(neop(ip),ip)=ie
      nonode(neop(ip),ip)=id ! nonode(j,i) is the number of the i-th node in the j-th related element
    enddo
  enddo

  ! Find all nodes related to the ip-th node, sort them the nele(:,ip), and generate ia
  allocate(ia(ntf1))
  ia(1)=1
  kn0sysb=0
  maxhb=0 ! maximum half bandwidth
  do ip=1,ntp
    neoe(ip)=0 ! the number of nodes that are related to the ip-th node
    do je=1,neop(ip)
      ie=nple(je,ip)
      inode=nonode(je,ip)
      node=nnode(ie)
      ndk=ndtok(node)
      n=nrn(inode,ndk) ! nrn is the number of nodes that are related to the ip-th node in the element
      if(neoe(ip).eq.0)then
        do i=1,n
          id=ndrn(i,inode,ndk) ! the i-th node that are related to the node in the element
          jp=lnde(id,ie)
          kn0sysb(jp)=1
          neoe(ip)=neoe(ip)+1 ! the number of nodes that are related to the ip-th node + 1
          nele(neoe(ip),ip)=jp ! add the related node into nele
        enddo
      else
        do i=1,n
          id=ndrn(i,inode,ndk) ! the i-th node that are related to the node in the element
          jp=lnde(id,ie)
          if(kn0sysb(jp).eq.1) cycle ! judge if this node have been calculated in other element
          kn0sysb(jp)=1
          neoe(ip)=neoe(ip)+1 ! the number of nodes that are related to the ip-th node + 1
          nele(neoe(ip),ip)=jp ! add the related node into nele
        enddo
      endif
    enddo
    ! sort nele(:,ip)
    icsize=neoe(ip)
    jcsize=4
    call qsort(nele(1:neoe(ip),ip),icsize,jcsize,my_int_compar) ! invoke qsort function in C language
    maxhb=max(maxhb,(nele(neoe(ip),ip)-nele(1,ip)+1)*ndf)
    ! generate ia
    ip0=(ip-1)*ndf
    do i=1,ndf
      ia(ip0+i+1)=ia(ip0+i)+neoe(ip)*ndf
    enddo
    ! set kn0sysb be 0
    do i=1,neoe(ip)
      kn0sysb(nele(i,ip))=0
    enddo
  enddo
  maxhb=maxhb/2
  deallocate(kn0sysb)

  ! generate ja
  ka=ia(ntf1)-1
  allocate(ja(ka))
  do ip=1,ntp
    ip0=(ip-1)*ndf
    do idf=1,ndf
      k=ia(ip0+idf)-1
      do j=1,neoe(ip)
        jp0=(nele(j,ip)-1)*ndf
        do jdf=1,ndf
          k=k+1
          ja(k)=jp0+jdf
        enddo
      enddo
    enddo
  enddo
  deallocate(nele)
  
  if(istrans.eq.1.and.ndf.eq.1)then ! transient heat conduction
    ! generate ir,jr
    allocate(ir(ntf1))
    ir=ia
    kr=ka
    allocate(jr(kr))
    jr=ja
    ! generate ic
    allocate(ic(ntf1))
    ic(1)=1
    do ip=1,ntp
      ip0=(ip-1)*ndf
      if(neop(ip).ne.1)then ! internal nodes only related with one element
        do i=1,ndf
          ic(ip0+i+1)=ic(ip0+i)
        enddo
        cycle
      endif
      ie=nple(1,ip)
      inode=nonode(1,ip)
      node=nnode(ie)
      ndk=ndtok(node)
      idtype=ndtype(inode,ndk)
      if(idtype.ne.1)then ! not element internal node
        do i=1,ndf
          ic(ip0+i+1)=ic(ip0+i)
        enddo
        cycle
      endif
      do i=1,ndf
        ic(ip0+i+1)=ic(ip0+i)+1
      enddo
    enddo
    ! generate jc
    kc=ic(ntf1)-1
    allocate(jc(kc))
    do ip=1,ntp
      ip0=(ip-1)*ndf
      do idf=1,ndf
        k=ic(ip0+idf+1)-ic(ip0+idf)
        if(k.eq.1)jc(ic(ip0+idf))=ip0+idf
      enddo
    enddo
  elseif(istrans.eq.1.and.ndf.ne.1)then ! dynamics
    ! distinguish the nodes into element boundary nodes and element internal nodes
    allocate(lbon(ntp),lion(ntp))
    lbon=0
    lion=0
    nobdf=0 ! the number of element boundary nodes
    noidf=0 ! the number of element internal nodes
    do ip=1,ntp
      if(neop(ip).ne.1)then ! internal nodes only related with one element
        nobdf=nobdf+1
        lbon(nobdf)=ip ! transfer element boundary node number to global number
        cycle
      endif
      ie=nple(1,ip)
      inode=nonode(1,ip)
      node=nnode(ie)
      ndk=ndtok(node)
      idtype=ndtype(inode,ndk)
      if(idtype.ne.1)then ! not element internal node
        nobdf=nobdf+1
        lbon(nobdf)=ip ! transfer element boundary node number to global number
        cycle
      endif
      noidf=noidf+1
      lion(noidf)=ip ! transfer element internal node number to global number
    enddo
    ! re-number the node
    allocate(lnof(ntp))
    ! element boundary nodes first
    do iobdf=1,nobdf
      ip=lbon(iobdf)
      lnof(ip)=iobdf ! transfer global number to the equation number
    enddo
    ! element internal nodes next
    do ioidf=1,noidf
      ip=lion(ioidf)
      lnof(ip)=ioidf+nobdf ! transfer the equation number to global number
    enddo
    deallocate(lion,lbon)
    ! form transformation matrix it, jt, t
    kt=ntf
    allocate(it(ntf1),jt(kt),jtp(kt),t(kt))
    do i=1,ntf1
      it(i)=i
    enddo
    do ip=1,ntp
      ipf=lnof(ip)
      do idf=1,ndf
        jt((ipf-1)*ndf+idf)=(ip-1)*ndf+idf
        jtp((ip-1)*ndf+idf)=(ipf-1)*ndf+idf
      enddo
    enddo
    t=1.d0
    !call directoutput_csr_matrix(ntf,ntf,kt,it,jt,t)
    ! form matrix irm, ic, jrm and jc of the dynamic problem
    allocate(irm(ntf+1),ic(ntf+1))
    do i=1,nobdf*ndf+1
      irm(i)=1
      ic(i)=1
    enddo
    do i=nobdf*ndf+2,ntf+1
      irm(i)=i-nobdf*ndf
      ic(i)=i-nobdf*ndf
    enddo
    krm=noidf*ndf
    kc=noidf*ndf
    allocate(jrm(krm),jc(kc))
    do i=1,noidf*ndf
      jrm(i)=i+nobdf*ndf
      jc(i)=i+nobdf*ndf
    enddo
  endif
  deallocate(nonode,nple)
  
end subroutine node_sort

integer(2) function my_int_compar(a1,a2)
  integer*4::a1, a2,a3
  integer sign
  a3=a1-a2
  my_int_compar=sign(1,a3)
end function my_int_compar

!-------------------- gen_dndx_dndxx --------------------
subroutine gen_dndx_dndxx
  ! Get the first and second order derivates of the shape functions with global coordinates
  use vary_arrays
  use fixed_values
  implicit real*8 (a-h,o-z)
  
  dimension xi(3)
  allocatable elcod(:,:)
  
  allocate(dndxa(3,mnode,mnode,nte))
  allocate(dndxxa(3,3,mnode,mnode,nte))
  
  do ie=1,nte ! element loop
    ietype=netype(ie)
    node=nnode(ie)
    ndk=ndtok(node)
    
    ! get coordiates of each node in the element
    allocate(elcod(3,node))
    do inode=1,node
      ip=lnde(inode,ie)
      elcod(1:ndim,inode)=cd(1:ndim,ip)
    enddo
    
    ! get dndx and dndxx
    do id=1,node ! source point loop
      ! instrict coordinates at source point
      call node_xi(id,ietype,xi)
      idtype=ndtype(id,ndk)
      if(idtype.gt.1) idtype=2
      ! get dndx and gjacb at gauss points
      select case(idtype)
      case(1) ! element internal nodes
        call eval_dndx_dndxx(ie,ndim,node,xi,elcod,ietype,dndxa(:,1:node,id,ie),dndxxa(:,:,1:node,id,ie),2)
      case(2) ! element boundary nodes
        call eval_dndx_dndxx(ie,ndim,node,xi,elcod,ietype,dndxa(:,1:node,id,ie),dndxxa(:,:,1:node,id,ie),1)
      end select
    enddo ! end of gauss point loop
    deallocate(elcod)
  enddo ! end of element loop
  
end subroutine gen_dndx_dndxx