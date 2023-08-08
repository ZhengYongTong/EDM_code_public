!****************************************************************************
  ! This file containing:
  ! Subroutine nastran_mesh: plot the whole model into .nas file
  ! Subroutine gridout: write node into .nas file
  ! Subroutine ndecimal: set the .nas file with no comma scheme
  ! Subroutine netout: write mesh into .nas file
  ! Subroutine tecp_gen: plot the the results into tecplot file
  ! Subroutine tecp_gen_trans: plot the the transient results into tecplot file
  ! Subroutine tecp_gen_head: write the head section of the binary tecplot file
  ! Subroutine tecp_gen_zone: write the zone section of the binary tecplot file
  ! Subroutine tecp_gen_data: write the data section of the binary tecplot file
  ! Subroutine dumpstring: transfer characters in scheme needed in binary tecplot file
  ! Subroutine tecp_gen_trans2: copy temporary file 99 to file 33 to complete data module of tecplot binary file for transiet problem
!****************************************************************************

subroutine nastran_mesh(ifile,ndim,nbdm,ntp,nte,cd,lnde,nnode,mnode,netype)
  implicit real*8 (a-h,o-z)
  dimension cd(ndim,ntp),lnde(mnode,nte),nnode(nte),cdw(3),ndchg(3,20,20),netype(nte)
  data (ndchg(1,2,i),i=1,2)/1,2/
  data (ndchg(1,3,i),i=1,3)/1,3,2/
  data (ndchg(2,3,i),i=1,3)/1,2,3/
  data (ndchg(2,4,i),i=1,4)/1,2,3,4/
  data (ndchg(2,6,i),i=1,6)/1,3,5,2,4,6/
  data (ndchg(2,8,i),i=1,8)/1,3,5,7,2,4,6,8/
  data (ndchg(2,9,i),i=1,9)/1,3,5,7,2,4,6,8,9/
  data (ndchg(3,8,i),i=1,8)/1,2,3,4,5,6,7,8/
  data (ndchg(3,20,i),i=1,20)/1,3,5,7, 13,15,17,19, 2,4,6,8, 9,10,11,12, 14,16,18,20/
  
  ndigit=8 ! output grids using eight digits per number
  write(ifile,'(10a)')'BEGIN BULK'
  do ip=1,ntp
    do i=1,ndim
      cdw(i)=cd(i,ip)
    enddo
    call gridout(ip,0,cdw,ndim,ifile,ndigit)
  enddo
  ! output connectives
  write(ifile,'("$")')
  do ie=1,nte
    node=nnode(ie)
    ietype=netype(ie)
    ny=1
    nz=1
    if(nbdm.eq.1) then
      nx=node
    elseif(nbdm.eq.2) then
      nx=node/8+2
      ny=nx
    else
      nx=node/32+2
      if(node.le.27) nx=node/20+2
      ny=nx
      nz=nx
    endif
    call netout(ndim,nbdm,node,lnde(1,ie),ie,1,1,ifile,ndchg,cd,ntp,nte,ietype,nx,ny,nz)
  enddo
  write(ifile,'(7a)')'ENDDATA'
  return
end subroutine nastran_mesh

subroutine gridout(igrid,icor,coor,ndime,nfile,ndigit)
  implicit real*8 (a-h,o-z)
  dimension coor(3)
  do i=ndime+1,3
   coor(i)=0.d0
  enddo
  ndc1=1
  ndc2=1
  ndc3=1
  call ndecimal(coor(1),ndc1,ndigit)
  if(ndime.ge.2) call ndecimal(coor(2),ndc2,ndigit)
  if(ndime.ge.3) call ndecimal(coor(3),ndc3,ndigit)
  write(nfile,"('GRID    ',2i8,',',1pe20.12,',',1pe20.12,',',1pe20.12)")igrid,icor,(coor(i),i=1,3)
end subroutine gridout

subroutine ndecimal(a,ndecm,ndigit)
  implicit real*8 (a-h,o-z)
  numb=iabs(int(dabs(a)+0.00001e0)/10)
  npow=0
  do i=1,ndigit
    if(numb.eq.0) exit
    numb=numb/10
    npow=npow+1
  enddo
  ndecm=ndigit-3-npow+int((1.d0+sign(1.d0,a))/2.d0+0.00001d0)
  if(dabs(a)+0.1e0.lt.1.e0) ndecm=ndecm+1
  if(dabs(a)+1.e0.gt.999999.d0.and.ndecm.gt.0) ndecm=ndecm-1
  if(ndecm.lt.0) ndecm=0
  ndecm1=ndecm+1
  w1=dabs(a-dble(int(dabs(a))))*10.e0**ndecm1+1.e-6
  w=dble(int(w1))/10.e0**ndecm1
  do i=ndecm,0,-1
    ww=w*10.e0**i
    if(dble(int(ww)).ne.ww) exit
  enddo
  if(i.lt.ndecm) ndecm=i+1
  return
end subroutine ndecimal

subroutine netout(ndim,ndimb,node,netinp,ie,nele1,npoi1,nfile,ndchg,cd,npoint,nelement,ietype,nx,ny,nz)
  implicit real*8 (a-h,o-z)
  dimension netinp(node),net(node),ndchg(3,20,20),cd(ndim,npoint)
  
  net=iabs(netinp)
  if(ndimb.eq.1) then
    write(nfile,11)ie+nele1-1,1,net(1)+npoi1-1,net(node)+npoi1-1
  elseif(ndimb.eq.2) then
    if(node.eq.8) then ! 8 node element
      write(nfile,88)ie+nele1-1,1,(net(j)+npoi1-1,j=1,6)
      write(nfile,888)(net(j)+npoi1-1,j=7,8)
    elseif(node.eq.6.or.node.eq.7) then ! 6 or 7 node triangular element
      write(nfile,33)ie+nele1-1,1,(net(j)+npoi1-1,j=1,6)
    elseif(node.le.25) then ! only plot outer sides of the element
      net(1)=1
      net(2)=nx
      net(3)=node
      net(4)=net(3)-nx+1
      write(nfile,44)ie+nele1-1,1,(netinp(net(j)),j=1,4)
    elseif(node.gt.25) then ! general cases
      do jj=1,ny-1
        jj0=(jj-1)*nx
        do ii=1,nx-1
           ii0=ii-1
           nelement=nelement+1
           net(1)=jj0+ii0+1
           net(2)=net(1)+1
           net(3)=net(2)+nx
           net(4)=net(3)-1
           write(nfile,44)nelement,1,(netinp(net(j)),j=1,4)
        enddo
      enddo
    endif
  elseif(ndimb.eq.3) then
    if(node.eq.8) then
      write(nfile,20)ie+nele1-1,1,(net(j)+npoi1-1,j=1,6)
      write(nfile,888)(net(j)+npoi1-1,j=7,8)
    elseif(node.eq.11) then
      write(nfile,30)ie+nele1-1,1,(net(j)+npoi1-1,j=1,6)
      write(nfile,888)(net(j)+npoi1-1,j=7,10)
    elseif(node.eq.20.or.node.eq.21) then
      write(nfile,201)ie+nele1-1,1,(net(ndchg(ndimb,node,j))+npoi1-1,j=1,node)
    elseif(node.le.125) then
      net(1)=1
      net(2)=nx
      net(3)=nx*ny
      net(4)=net(3)-nx+1
      net(5)=net(1)+nx*ny*(nz-1)
      net(6)=net(5)+nx-1
      net(7)=node
      net(8)=net(7)-nx+1
      write(nfile,20)ie+nele1-1,1,(netinp(net(j))+npoi1-1,j=1,6)
      write(nfile,888)(netinp(net(j))+npoi1-1,j=7,8)
    endif
  endif
  
11 format('CBAR    ',6i8)
33 format('CTRIA6  ',8i8)
44 format('CQUAD4  ',6i8)
88 format('CQUAD8  ',8i8,'        ')
888 format('        ',8i8)
20 format('CHEXA   ',8i8,'        ')
201 format('CHEXA   ',8i8,'        ',/,'        ',8i8,'        ',/,'        ',8i8)
30 format('CTETRA  ',8i8,'        ')
  
end subroutine netout

! *******************************************************************
subroutine tecp_gen(ifile,disp,sigma,streng,cd,lnde,ntp,nte,ndim,ndf,nsig,mnode,ncomp)
  implicit none
  integer::ntp,nte,ndim,ndf,nsig,mnode,ncomp,ifile,zonetype,i
  integer::lnde(mnode,nte)
  real*8::disp(ndf,ntp),sigma(nsig,ntp),streng(ncomp,ntp),cd(ndim,ntp)
  real*8,allocatable::usum(:)
  real*4::zonemarker,eohmarker
  real*8:: solutiontimes
  
  zonemarker=299.0
  eohmarker=357.0
  ! zonetype=0-ordered,1-felineseg,2-fetriangle,3-fequadrilateral,
  ! 4-fetetrahedron,5-febrick,6-fepolygon,7-fepolyhedron
  if(ndim.eq.2)then
    if(mnode.eq.7)then
      zonetype=2 ! fetriangle
    else
      zonetype=3 ! fequadrilateral
    endif
  elseif(ndim.eq.3)then
    if(mnode.eq.11)then
      zonetype=4 ! fetetrahedron
    else
      zonetype=5 ! febrick
    endif
  endif
  ! compute total disp
  allocate(usum(ntp))
  if(ndf.gt.1)then
    if(ndim.eq.2)then
      do i=1,ntp
        usum(i)=dsqrt(disp(1,i)**2.0d0+disp(2,i)**2.0d0)
      enddo
    else
      do i=1,ntp
        usum(i)=dsqrt(disp(1,i)**2.0d0+disp(2,i)**2.0d0+disp(3,i)**2.0d0)
      enddo
    endif
  endif
  
  call tecp_gen_head(ifile,ndim,ndf,nsig,ncomp)
  
  solutiontimes=1.0
  call tecp_gen_zone(ifile,ntp,nte,zonetype,zonemarker,solutiontimes,1)
  
  write(ifile) eohmarker ! separate the header and the data with an eohmarker=357.0
  
  call tecp_gen_data(ifile,ntp,nte,ndim,ndf,nsig,mnode,ncomp,cd,lnde,disp,sigma,streng,usum,zonemarker)
  deallocate(usum)
  
end subroutine tecp_gen

subroutine tecp_gen_trans(ifile,ifile1,disp,sigma,streng,cd,lnde,ntp,nte,ndim,ndf,nsig,mnode,ncomp,istep,deltat)
  implicit none
  integer::ntp,nte,ndim,ndf,nsig,mnode,ncomp,ifile,ifile1,zonetype,i
  integer::istep
  integer::lnde(mnode,nte)
  real*8::disp(ndf,ntp),sigma(nsig,ntp),streng(ncomp,ntp),cd(ndim,ntp)
  real*8,allocatable::usum(:)
  real*4::zonemarker,eohmarker
  real*8:: solutiontimes,deltat
  
  zonemarker=299.0
  eohmarker=357.0
  ! zonetype=0-ordered,1-felineseg,2-fetriangle,3-fequadrilateral,
  ! 4-fetetrahedron,5-febrick,6-fepolygon,7-fepolyhedron
  if(ndim.eq.2)then
    if(mnode.eq.7)then
      zonetype=2 ! fetriangle
    else
      zonetype=3 ! fequadrilateral
    endif
  elseif(ndim.eq.3)then
    if(mnode.eq.11)then
      zonetype=4 ! fetetrahedron
    else
      zonetype=5 ! febrick
    endif
  endif
  ! compute total disp
  allocate(usum(ntp))
  if(ndf.gt.1)then
    if(ndim.eq.2)then
      do i=1,ntp
        usum(i)=dsqrt(disp(1,i)**2.0d0+disp(2,i)**2.0d0)
      enddo
    else
      do i=1,ntp
        usum(i)=dsqrt(disp(1,i)**2.0d0+disp(2,i)**2.0d0+disp(3,i)**2.0d0)
      enddo
    endif
  endif
  
  if(istep.eq.0)then
    call tecp_gen_head(ifile,ndim,ndf,nsig,ncomp)
  endif
  
  solutiontimes=dble(istep)*deltat
  call tecp_gen_zone(ifile,ntp,nte,zonetype,zonemarker,solutiontimes,istep)
  
  call tecp_gen_data(ifile1,ntp,nte,ndim,ndf,nsig,mnode,ncomp,cd,lnde,disp,sigma,streng,usum,zonemarker)
  deallocate(usum)
  
end subroutine tecp_gen_trans

subroutine tecp_gen_head(ifile,ndim,ndf,nsig,ncomp)
  implicit none
  integer::ndim,ndf,nsig,ncomp,ifile
  character(40)::title,varname
  
  ! 1. magic number, version number
  write(ifile) '#!TDV112'
  ! 2. integer value of 1.
  write(ifile) 1
  ! 3. title and variable names.
  write(ifile) 0
  title='Results' ! the title
  call dumpstring(title,ifile)
  if(ndf.gt.1)then
    write(ifile) ndim+ndf+1+nsig+ncomp ! number of variables in the datafile
  else
    write(ifile) ndim+ndf+nsig+ncomp ! number of variables in the datafile
  endif
  varname='X-coordinate'
  call dumpstring(varname,ifile)
  varname='Y-coordinate'
  call dumpstring(varname,ifile)
  if(ndim.eq.3)then
      varname='Z-coordinate'
      call dumpstring(varname,ifile)
  endif
  if(ndf.eq.1)then
    varname='Temperature'
    call dumpstring(varname,ifile)
    varname='X-heat-flux'
    call dumpstring(varname,ifile)
    varname='Y-heat-flux'
    call dumpstring(varname,ifile)
    if(ndim.eq.3)then
      varname='Z-heat-flux'
      call dumpstring(varname,ifile)
    endif
    varname='Heat-flux'
    call dumpstring(varname,ifile)
  else
    varname='X-displacement'
    call dumpstring(varname,ifile)
    varname='Y-displacement'
    call dumpstring(varname,ifile)
    if(ndim.eq.3)then
      varname='Z-displacement'
      call dumpstring(varname,ifile)
    endif
    varname='Displacement'
    call dumpstring(varname,ifile)
    if(nsig.eq.3.or.nsig.eq.4)then
      varname='Stress-XX'
      call dumpstring(varname,ifile)
      varname='Stress-YY'
      call dumpstring(varname,ifile)
      varname='Stress-XY'
      call dumpstring(varname,ifile)
      if(nsig.eq.4)then
        varname='Stress-ZZ'
        call dumpstring(varname,ifile)
      endif
    else ! nsig.eq.6
      varname='Stress-XX'
      call dumpstring(varname,ifile)
      varname='Stress-YY'
      call dumpstring(varname,ifile)
      varname='Stress-ZZ'
      call dumpstring(varname,ifile)
      varname='Stress-XY'
      call dumpstring(varname,ifile)
      varname='Stress-YZ'
      call dumpstring(varname,ifile)
      varname='Stress-XZ'
      call dumpstring(varname,ifile)
    endif
    varname='Tresca'
    call dumpstring(varname,ifile)
    varname='Mises'
    call dumpstring(varname,ifile)
    varname='Mohr-Coulomb'
    call dumpstring(varname,ifile)
    varname='Drucker-Prager'
    call dumpstring(varname,ifile)
  endif
  
end subroutine tecp_gen_head

subroutine tecp_gen_zone(ifile,ntp,nte,zonetype,zonemarker,solutiontimes,istep)
  implicit none
  integer::ntp,nte,ifile,zonetype
  integer::istep
  character(40)::zonename1
  real*4::zonemarker
  real*8:: solutiontimes
  
  ! 4. zones
  write(ifile) zonemarker
  write(zonename1,'(i6)') istep ! write istep into zone name
  zonename1="STEP 1 INCR "//zonename1 ! zone name.
  call dumpstring(zonename1,ifile)
  write(ifile) -1 ! parentzone
  write(ifile) 0  ! strandid
  write(ifile) solutiontimes ! solution times
  write(ifile) -1 ! not used. set to -1.
  write(ifile) zonetype ! zonetype
  write(ifile) 0 ! specify var location.0 = don¡¯t specify, all data is located at the nodes.
  write(ifile) 0 ! are raw local 1-to-1 face neighbors supplied
  write(ifile) 0 ! number of miscellaneous user-defined face neighbor connections
  write(ifile) ntp
  write(ifile) nte
  ! icelldim,jcelldim, kcelldim (for future use; set to zero)
  write(ifile) 0 
  write(ifile) 0
  write(ifile) 0
  write(ifile) 0 ! no more auxiliary name/value pairs.
  
end subroutine tecp_gen_zone

subroutine tecp_gen_data(ifile,ntp,nte,ndim,ndf,nsig,mnode,ncomp,cd,lnde,disp,sigma,streng,usum,zonemarker)
  implicit none
  integer::ntp,nte,ndim,ndf,nsig,mnode,ncomp,ifile,i
  integer::istep
  integer::lnde(mnode,nte)
  real*8::disp(ndf,ntp),sigma(nsig,ntp),streng(ncomp,ntp),cd(ndim,ntp),usum(ntp)
  real*4::zonemarker
  real*8:: maxvar,minvar
  
  ! 5. data section
  write(ifile) zonemarker
  do i=1,ndim+ndf+nsig+ncomp
    write(ifile)2 ! variable data format, 1=float, 2=double, 3=longint, 4=shortint, 5=byte, 6=bit
  enddo
  if(ndf.gt.1)then
    write(ifile)2 ! variable data format, 1=float, 2=double, 3=longint, 4=shortint, 5=byte, 6=bit
  endif
  write(ifile) 0 ! has passive variables: 0 = no, 1 = yes.
  write(ifile) 0 ! has variable sharing 0 = no, 1 = yes
  write(ifile) -1 ! zero based zone number to share connectivity list with (-1 = no sharing).
  do i=1,ndim
    minvar=minval(cd(i,:))
    maxvar=maxval(cd(i,:))
    write(ifile) minvar
    write(ifile) maxvar
  enddo
  do i=1,ndf
    minvar=minval(disp(i,:))
    maxvar=maxval(disp(i,:))
    write(ifile) minvar
    write(ifile) maxvar
  enddo
  if(ndf.gt.1)then
    minvar=minval(usum)
    maxvar=maxval(usum)
    write(ifile) minvar
    write(ifile) maxvar
  endif
  do i=1,nsig
    minvar=minval(sigma(i,:))
    maxvar=maxval(sigma(i,:))
    write(ifile) minvar
    write(ifile) maxvar
  enddo
  do i=1,ncomp
    minvar=minval(streng(i,:))
    maxvar=maxval(streng(i,:))
    write(ifile) minvar
    write(ifile) maxvar
  enddo
  do i=1,ndim
    write(ifile) cd(i,:)
  enddo
  do i=1,ndf
    write(ifile) disp(i,:)
  enddo
  if(ndf.gt.1)then
    write(ifile) usum(:)
  endif
  do i=1,nsig
    write(ifile) sigma(i,:)
  enddo
   do i=1,ncomp
    write(ifile) streng(i,:)
  enddo
  if(ndim.eq.2)then
    if(mnode.eq.7)then
      do i=1,nte
        write(ifile) lnde((/1,2,3/),i)-1
      enddo
    elseif(mnode.eq.8)then
      do i=1,nte
        write(ifile) lnde((/1,2,3,4/),i)-1
      enddo
    elseif(mnode.eq.9)then
      do i=1,nte
        write(ifile) lnde((/1,3,9,7/),i)-1
      enddo
    elseif(mnode.eq.16)then
      do i=1,nte
        write(ifile) lnde((/1,4,16,13/),i)-1
      enddo
    elseif(mnode.eq.25)then
      do i=1,nte
        write(ifile) lnde((/1,5,25,21/),i)-1
      enddo
    endif
  else ! ndim.eq.3
    if(mnode.eq.11)then
      do i=1,nte
        write(ifile) lnde((/1,2,3,4/),i)-1
      enddo
    elseif(mnode.eq.20.or.mnode.eq.21)then
      do i=1,nte
        write(ifile) lnde((/1,2,3,4,5,6,7,8/),i)-1
      enddo
    elseif(mnode.eq.27)then
      do i=1,nte
        write(ifile) lnde((/1,3,9,7,19,21,27,25/),i)-1
      enddo
    elseif(mnode.eq.64)then
      do i=1,nte
        write(ifile) lnde((/1,4,16,13,49,52,64,61/),i)-1
      enddo
    elseif(mnode.eq.125)then
      do i=1,nte
        write(ifile) lnde((/1,5,25,21,101,105,125,121/),i)-1
      enddo
    endif
  endif
  
end subroutine tecp_gen_data

subroutine dumpstring(instring,ifile)
  implicit none
  character(40)::instring
  integer::lenth,ii,i,ifile
  lenth=len_trim(instring)
  do ii=1,lenth
    i=ichar(instring(ii:ii))
    write(ifile) i
  end do
  write(ifile) 0
  return
end subroutine dumpstring

subroutine tecp_gen_trans2(ifile,ifile1,ndim,ndf,nsig,mnode,ntp,nte,istep)
  implicit none
  integer::ndim,ndf,nsig,ncomp,mnode,ifile,ifile1,ntp,nte,i,j,istep
  integer::n
  real*4::w
  real*8::z
  integer,allocatable::m(:)
  real*8,allocatable::y(:)
  
  write(ifile) 357.0 ! separate the header and the data with an eohmarker=357.0
  rewind(ifile1)
  ncomp=1
  if(ndf.gt.1) ncomp=4
  do j=0,istep
    read(ifile1) w
    write(ifile) w
    do i=1,ndim+ndf+nsig+ncomp
      read(ifile1) n
      write(ifile) n 
    enddo
    if(ndf.gt.1)then
      read(ifile1) n
      write(ifile) n
    endif
    read(ifile1) n
    write(ifile) n ! has passive variables: 0 = no, 1 = yes.
    read(ifile1) n
    write(ifile) n ! has variable sharing 0 = no, 1 = yes
    read(ifile1) n
    write(ifile) n ! zero based zone number to share connectivity list with (-1 = no sharing).
    do i=1,ndim
      read(ifile1) z
      write(ifile) z
      read(ifile1) z
      write(ifile) z
    enddo
    do i=1,ndf
      read(ifile1) z
      write(ifile) z
      read(ifile1) z
      write(ifile) z
    enddo
    if(ndf.gt.1)then
      read(ifile1) z
      write(ifile) z
      read(ifile1) z
      write(ifile) z
    endif
    do i=1,nsig
      read(ifile1) z
      write(ifile) z
      read(ifile1) z
      write(ifile) z
    enddo
    do i=1,ncomp
      read(ifile1) z
      write(ifile) z
      read(ifile1) z
      write(ifile) z
    enddo
    allocate(y(ntp))
    do i=1,ndim
      read(ifile1) y
      write(ifile) y
    enddo
    do i=1,ndf
      read(ifile1) y
      write(ifile) y
    enddo
    if(ndf.gt.1)then
      read(ifile1) y
      write(ifile) y
    endif
    do i=1,nsig
      read(ifile1) y
      write(ifile) y
    enddo
    do i=1,ncomp
      read(ifile1) y
      write(ifile) y
    enddo
    deallocate(y)
    if(ndim.eq.2)then
      if(mnode.eq.7)then
        allocate(m(3))
        do i=1,nte
          read(ifile1) m
          write(ifile) m
        enddo
        deallocate(m)
      elseif(mnode.eq.8.or.mnode.eq.9.or.mnode.eq.16.or.mnode.eq.25)then
        allocate(m(4))
        do i=1,nte
          read(ifile1) m
          write(ifile) m
        enddo
        deallocate(m)
      endif
    else ! ndim.eq.3
      if(mnode.eq.11)then
        allocate(m(4))
        do i=1,nte
          read(ifile1) m
          write(ifile) m
        enddo
        deallocate(m)
      elseif(mnode.eq.20.or.mnode.eq.21.or.mnode.eq.27.or.mnode.eq.64.or.mnode.eq.125)then
        allocate(m(8))
        do i=1,nte
          read(ifile1) m
          write(ifile) m
        enddo
        deallocate(m)
      endif
    endif
  enddo
  close(ifile1,status='delete')
end subroutine tecp_gen_trans2