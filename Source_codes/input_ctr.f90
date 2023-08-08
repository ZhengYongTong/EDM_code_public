!****************************************************************************
  ! This file containing:
  ! Subroutine input_ctr: read the control variables
  ! Subroutine block_data: set some values of contants
  ! Subroutine netypend: set some values about element
!****************************************************************************

subroutine input_ctr
  use fixed_values
  implicit real*8 (a-h,o-z)
  character*6 :: title(20)
  read(5,'(20a6)')(title(i),i=1,20) ! title of problem
  write(7,'(/1x,20a6)')(title(i),i=1,20)
  ! ndim is the dimension of the problem
  ! ndf is the number of the variables to be solved of one node in the problem
  ! nsig is the number of heat flux components or stress components
  ! mnode is the maximum number of nodes in the element
  ! ntp is the number of total nodes in the probelm
  ! nte is the number of total elements in the problem
  ! msolver is the type of linear algebraic solver
  write(7,'(//,'' ndim ndf nsig mnode    ntp     nte   msolver'')') ! control variables
  read(5,*)ndim,ndf,nsig,mnode,ntp,nte,msolver
  write(7,'(3i4,i7,2i8,i7,3i6)')ndim,ndf,nsig,mnode,ntp,nte,msolver
  ! ntf is the orders of the final system
  ntf=ntp*ndf
  ntf1=ntf+1
  if(ndf.eq.1) nsig=ndim ! In case the input file parameters are wrong
  ! nsigp:the order of material matrix of one node
  if(nsig.eq.4) then ! for plane strain problems
    nsigp=3
  else
    nsigp=nsig
  endif
end subroutine input_ctr

subroutine block_data
  use fixed_values
  implicit real*8 (a-h,o-z)
  ! set delta function and pi
  dlt=reshape((/1.,0.,0., 0.,1.,0., 0.,0.,1./),(/3,3/))
  pi=4.d0*datan(1.d0)
  ! set multipliers for repeated stress components
  if(ndim.eq.2.and.ndf.eq.2) then
    iny(1:4)=(/1,3,3,2/)
  elseif(ndim.eq.3.and.ndf.eq.3) then
    iny=(/1,4,6,  4,2,5,  6,5,3/) ! i,j to position of sigma(i,j)
  endif

  ! nfacem is the maximum number of the faces in an element
  if(mnode.eq.7)then
    nfacem=3
  elseif(mnode.eq.11)then
    nfacem=4
  else
    nfacem=2*ndim
  endif
  
  ! set the value about element
  call netypend(ndim,ndtok,ndface,nface,ndtype,npface,nrn,ndrn)
  
  return
end subroutine block_data

subroutine netypend(ndim,ndtok,ndface,nface,ndtype,npface,nrn,ndrn)
  implicit real*8 (a-h,o-z)
  dimension ndtok(125),ndface(25,6,5),nface(5),ndtype(125,5),npface(5),nrn(125,5),ndrn(125,125,5)
  
  if(ndim.eq.3) then ! 3d
    ndtok((/11,21,27,64,125/))=(/5,1,2,3,4/) ! for 3d, used to replace the number of nodes in corresponding element
   
    ! nface is the number of the faces in an element
    nface(1:4)=6
    nface(5)=4
    
    ! npface is the number of the nodes in an element face
    npface(1)=8
    npface(2)=9
    npface(3)=16
    npface(4)=25
    npface(5)=6
    
    ! ndface is the node number of each face in the element
    ! ndface(node in this face, face, element ndtok)
    ! node.eq.21
    ndface(1:8,1,1)=(/4,3,2,1, 11,10,9,12/)
    ndface(1:8,2,1)=(/2,3,7,6, 10,19,14,18/)
    ndface(1:8,3,1)=(/5,6,7,8, 13,14,15,16/)
    ndface(1:8,4,1)=(/5,8,4,1, 16,20,12,17/)
    ndface(1:8,5,1)=(/1,2,6,5, 9,18,13,17/)
    ndface(1:8,6,1)=(/8,7,3,4, 15,19,11,20/)
    ! node.eq.27
    ndface(1:9,1,2)=(/7,8,9, 4,5,6, 1,2,3/)
    ndface(1:9,2,2)=(/3,6,9, 12,15,18, 21,24,27/)
    ndface(1:9,3,2)=(/19,20,21, 22,23,24, 25,26,27/)
    ndface(1:9,4,2)=(/19,22,25, 10,13,16, 1,4,7/)
    ndface(1:9,5,2)=(/1,2,3, 10,11,12, 19,20,21/)
    ndface(1:9,6,2)=(/25,26,27, 16,17,18, 7,8,9/)
    ! node.eq.64
    ndface(1:16,1,3)=(/13,14,15,16, 9,10,11,12, 5,6,7,8, 1,2,3,4/)
    ndface(1:16,2,3)=(/4,8,12,16, 20,24,28,32, 36,40,44,48,52,56,60,64/)
    ndface(1:16,3,3)=(/49,50,51,52, 53,54,55,56, 57,58,59,60,61,62,63,64/)
    ndface(1:16,4,3)=(/49,53,57,61, 33,37,41,45, 17,21,25,29,1,5,9,13/)
    ndface(1:16,5,3)=(/1,2,3,4, 17,18,19,20, 33,34,35,36,49,50,51,52/)
    ndface(1:16,6,3)=(/61,62,63,64, 45,46,47,48, 29,30,31,32,13,14,15,16/)
    ! node.eq.125
    ndface(:,1,4)=(/21,22,23,24,25, 16,17,18,19,20, 11,12,13,14,15,6,7,8,9,10, 1,2,3,4,5/)
    ndface(:,2,4)=(/5,10,15,20,25, 30,35,40,45,50, 55,60,65,70,75,80,85,90,95,100, 105,110,115,120,125/)
    ndface(:,3,4)=(/101,102,103,104,105, 106,107,108,109,110,111,112,113,114,115, 116,117,118,119,120,121,122,123,124,125/)
    ndface(:,4,4)=(/101,106,111,116,121, 76,81,86,91,96,51,56,61,66,71, 26,31,36,41,46, 1,6,11,16,21/)
    ndface(:,5,4)=(/1,2,3,4,5, 26,27,28,29,30, 51,52,53,54,55,76,77,78,79,80, 101,102,103,104,105/)
    ndface(:,6,4)=(/121,122,123,124,125, 96,97,98,99,100,71,72,73,74,75, 46,47,48,49,50, 21,22,23,24,25/)
    ! node.eq.11
    ndface(1:6,1,5)=(/1,3,2,7,6,5/)
    ndface(1:6,2,5)=(/1,2,4,5,9,8/)
    ndface(1:6,3,5)=(/2,3,4,6,10,9/)
    ndface(1:6,4,5)=(/1,4,3,8,10,7/)
    
    ! ndtype is the type of nodes in the element
    ! ndtype(node, element ndtok)
    !        = 1, element internal nodes
    !        = 2, element surface nodes
    ! node.eq.21
    ndtype(21,1)=1
    ndtype(1:20,1)=2
    ! node.eq.27
    ndtype(1:27,2)=2
    ndtype(14,2)=1
    ! node.eq.64
    ndtype(1:64,3)=2
    ndtype((/22,23,26,27,38,39,42,43/),3)=1
    ! node.eq.125
    ndtype(1:125,4)=2
    ndtype((/32,33,34, 37,38,39, 42,43,44,57,58,59, 62,63,64, 67,68,69,82,83,84, 87,88,89, 92,93,94/),4)=1
    ! node.eq.11
    ndtype(1:10,5)=2 
    ndtype(11,5)=1
    
    ! nrn is the number of nodes that are related to the node in the element
    ! nrn(node, element ndtok)
    ! node=21
    nrn(1:8,1)=7
    nrn(9:21,1)=21
    ! node=27
    nrn(1:27,2)=7
    nrn(14,2)=27
    ! node=64
    nrn(1:64,3)=10
    nrn((/22,23,26,27,38,39,42,43/),3)=64
    ! node=125
    nrn(1:125,4)=13
    nrn((/32,33,34,37,38,39,42,43,44,57,58,59,62,63,64,67,68,69,82,83,84,87,88,89,92,93,94/),4)=125
    ! node=11
    nrn(1:11,5)=11
    
    ! ndrn is the node number that are related to the node in the element
    ! ndrn(node, related node, element ndtok)= related node number in the element
    ! node=21
    ndrn(1:7,1,1)=(/1,2,4,5,9,12,17/)
    ndrn(1:7,2,1)=(/1,2,3,6,9,10,18/)
    ndrn(1:7,3,1)=(/2,3,4,7,10,11,19/)
    ndrn(1:7,4,1)=(/1,3,4,8,11,12,20/)
    ndrn(1:7,5,1)=(/1,5,6,8,13,16,17/)
    ndrn(1:7,6,1)=(/2,5,6,7,13,14,18/)
    ndrn(1:7,7,1)=(/3,6,7,8,14,15,19/)
    ndrn(1:7,8,1)=(/4,5,7,8,15,16,20/)
    do i=9,21
      ndrn(1:21,i,1)=(/(j,j=1,21)/)
    enddo
    ! node=27
    ndrn(1:7,1,2)=(/1,2,3,4,7,10,19/)
    ndrn(1:7,2,2)=(/1,2,3,5,8,11,20/)
    ndrn(1:7,3,2)=(/1,2,3,6,9,12,21/)
    ndrn(1:7,4,2)=(/1,4,5,6,7,13,22/)
    ndrn(1:7,5,2)=(/2,4,5,6,8,14,23/)
    ndrn(1:7,6,2)=(/3,4,5,6,9,15,24/)
    ndrn(1:7,7,2)=(/1,4,7,8,9,16,25/)
    ndrn(1:7,8,2)=(/2,5,7,8,9,17,26/)
    ndrn(1:7,9,2)=(/3,6,7,8,9,18,27/)
    ndrn(1:7,10,2)=(/1,10,11,12,13,16,19/)
    ndrn(1:7,11,2)=(/2,10,11,12,14,17,20/)
    ndrn(1:7,12,2)=(/3,10,11,12,15,18,21/)
    ndrn(1:7,13,2)=(/4,10,13,14,15,16,22/)
    ndrn(1:27,14,2)=(/(i,i=1,27)/)
    ndrn(1:7,15,2)=(/6,12,13,14,15,18,24/)
    ndrn(1:7,16,2)=(/7,10,13,16,17,18,25/)
    ndrn(1:7,17,2)=(/8,11,14,16,17,18,26/)
    ndrn(1:7,18,2)=(/9,12,15,16,17,18,27/)
    ndrn(1:7,19,2)=(/1,10,19,20,21,22,25/)
    ndrn(1:7,20,2)=(/2,11,19,20,21,23,26/)
    ndrn(1:7,21,2)=(/3,12,19,20,21,24,27/)
    ndrn(1:7,22,2)=(/4,13,19,22,23,24,25/)
    ndrn(1:7,23,2)=(/5,14,20,22,23,24,26/)
    ndrn(1:7,24,2)=(/6,15,21,22,23,24,27/)
    ndrn(1:7,25,2)=(/7,16,19,22,25,26,27/)
    ndrn(1:7,26,2)=(/8,17,20,23,25,26,27/)
    ndrn(1:7,27,2)=(/9,18,21,24,25,26,27/)
    ! node=64
    ndrn(1:10,1,3)=(/1,2,3,4,5,9,13,17,33,49/)
    ndrn(1:10,2,3)=(/1,2,3,4,6,10,14,18,34,50/)
    ndrn(1:10,3,3)=(/1,2,3,4,7,11,15,19,35,51/)
    ndrn(1:10,4,3)=(/1,2,3,4,8,12,16,20,36,52/)
    ndrn(1:10,5,3)=(/1,5,6,7,8,9,13,21,37,53/)
    ndrn(1:10,6,3)=(/2,5,6,7,8,10,14,22,38,54/)
    ndrn(1:10,7,3)=(/3,5,6,7,8,11,15,23,39,55/)
    ndrn(1:10,8,3)=(/4,5,6,7,8,12,16,24,40,56/)
    ndrn(1:10,9,3)=(/1,5,9,10,11,12,13,25,41,57/)
    ndrn(1:10,10,3)=(/2,6,9,10,11,12,14,26,42,58/)
    ndrn(1:10,11,3)=(/3,7,9,10,11,12,15,27,43,59/)
    ndrn(1:10,12,3)=(/4,8,9,10,11,12,16,28,44,60/)
    ndrn(1:10,13,3)=(/1,5,9,13,14,15,16,29,45,61/)
    ndrn(1:10,14,3)=(/2,6,10,13,14,15,16,30,46,62/)
    ndrn(1:10,15,3)=(/3,7,11,13,14,15,16,31,47,63/)
    ndrn(1:10,16,3)=(/4,8,12,13,14,15,16,32,48,64/)
    ndrn(1:10,17,3)=(/1,17,18,19,20,21,25,29,33,49/)
    ndrn(1:10,18,3)=(/2,17,18,19,20,22,26,30,34,50/)
    ndrn(1:10,19,3)=(/3,17,18,19,20,23,27,31,35,51/)
    ndrn(1:10,20,3)=(/4,17,18,19,20,24,28,32,36,52/)
    ndrn(1:10,21,3)=(/5,17,21,22,23,24,25,29,37,53/)
    ndrn(1:64,22,3)=(/(j,j=1,64)/)
    ndrn(1:64,23,3)=(/(j,j=1,64)/)
    ndrn(1:10,24,3)=(/8,20,21,22,23,24,28,32,40,56/)
    ndrn(1:10,25,3)=(/9,17,21,25,26,27,28,29,41,57/)
    ndrn(1:64,26,3)=(/(j,j=1,64)/)
    ndrn(1:64,27,3)=(/(j,j=1,64)/)
    ndrn(1:10,28,3)=(/12,20,24,25,26,27,28,32,44,60/)
    ndrn(1:10,29,3)=(/13,17,21,25,29,30,31,32,45,61/)
    ndrn(1:10,30,3)=(/14,18,22,26,29,30,31,32,46,62/)
    ndrn(1:10,31,3)=(/15,19,23,27,29,30,31,32,47,63/)
    ndrn(1:10,32,3)=(/16,20,24,28,29,30,31,32,48,64/)
    ndrn(1:10,33,3)=(/1,17,33,34,35,36,37,41,45,49/)
    ndrn(1:10,34,3)=(/2,18,33,34,35,36,38,42,46,50/)
    ndrn(1:10,35,3)=(/3,19,33,34,35,36,39,43,47,51/)
    ndrn(1:10,36,3)=(/4,20,33,34,35,36,40,44,48,52/)
    ndrn(1:10,37,3)=(/5,21,33,37,38,39,40,41,45,53/)
    ndrn(1:64,38,3)=(/(j,j=1,64)/)
    ndrn(1:64,39,3)=(/(j,j=1,64)/)
    ndrn(1:10,40,3)=(/8,24,36,37,38,39,40,44,48,56/)
    ndrn(1:10,41,3)=(/9,25,33,37,41,42,43,44,45,57/)
    ndrn(1:64,42,3)=(/(j,j=1,64)/)
    ndrn(1:64,43,3)=(/(j,j=1,64)/)
    ndrn(1:10,44,3)=(/12,28,36,40,41,42,43,44,48,60/)
    ndrn(1:10,45,3)=(/13,29,33,37,41,45,46,47,48,61/)
    ndrn(1:10,46,3)=(/14,30,34,38,42,45,46,47,48,62/)
    ndrn(1:10,47,3)=(/15,31,35,39,43,45,46,47,48,63/)
    ndrn(1:10,48,3)=(/16,32,36,40,44,45,46,47,48,64/)
    ndrn(1:10,49,3)=(/1,17,33,49,50,51,52,53,57,61/)
    ndrn(1:10,50,3)=(/2,18,34,49,50,51,52,54,58,62/)
    ndrn(1:10,51,3)=(/3,19,35,49,50,51,52,55,59,63/)
    ndrn(1:10,52,3)=(/4,20,36,49,50,51,52,56,60,64/)
    ndrn(1:10,53,3)=(/5,21,37,49,53,54,55,56,57,61/)
    ndrn(1:10,54,3)=(/6,22,38,50,53,54,55,56,58,62/)
    ndrn(1:10,55,3)=(/7,23,39,51,53,54,55,56,59,63/)
    ndrn(1:10,56,3)=(/8,24,40,52,53,54,55,56,60,64/)
    ndrn(1:10,57,3)=(/9,25,41,49,53,57,58,59,60,61/)
    ndrn(1:10,58,3)=(/10,26,42,50,54,57,58,59,60,62/)
    ndrn(1:10,59,3)=(/11,27,43,51,55,57,58,59,60,63/)
    ndrn(1:10,60,3)=(/12,28,44,52,56,57,58,59,60,64/)
    ndrn(1:10,61,3)=(/13,29,45,49,53,57,61,62,63,64/)
    ndrn(1:10,62,3)=(/14,30,46,50,54,58,61,62,63,64/)
    ndrn(1:10,63,3)=(/15,31,47,51,55,59,61,62,63,64/)
    ndrn(1:10,64,3)=(/16,32,48,52,56,60,61,62,63,64/)
    ! node=125
    ndrn(1:13,1,4)=(/1,2,3,4,5,6,11,16,21,26,51,76,101/)
    ndrn(1:13,2,4)=(/1,2,3,4,5,7,12,17,22,27,52,77,102/)
    ndrn(1:13,3,4)=(/1,2,3,4,5,8,13,18,23,28,53,78,103/)
    ndrn(1:13,4,4)=(/1,2,3,4,5,9,14,19,24,29,54,79,104/)
    ndrn(1:13,5,4)=(/1,2,3,4,5,10,15,20,25,30,55,80,105/)
    ndrn(1:13,6,4)=(/1,6,7,8,9,10,11,16,21,31,56,81,106/)
    ndrn(1:13,7,4)=(/2,6,7,8,9,10,12,17,22,32,57,82,107/)
    ndrn(1:13,8,4)=(/3,6,7,8,9,10,13,18,23,33,58,83,108/)
    ndrn(1:13,9,4)=(/4,6,7,8,9,10,14,19,24,34,59,84,109/)
    ndrn(1:13,10,4)=(/5,6,7,8,9,10,15,20,25,35,60,85,110/)
    ndrn(1:13,11,4)=(/1,6,11,12,13,14,15,16,21,36,61,86,111/)
    ndrn(1:13,12,4)=(/2,7,11,12,13,14,15,17,22,37,62,87,112/)
    ndrn(1:13,13,4)=(/3,8,11,12,13,14,15,18,23,38,63,88,113/)
    ndrn(1:13,14,4)=(/4,9,11,12,13,14,15,19,24,39,64,89,114/)
    ndrn(1:13,15,4)=(/5,10,11,12,13,14,15,20,25,40,65,90,115/)
    ndrn(1:13,16,4)=(/1,6,11,16,17,18,19,20,21,41,66,91,116/)
    ndrn(1:13,17,4)=(/2,7,12,16,17,18,19,20,22,42,67,92,117/)
    ndrn(1:13,18,4)=(/3,8,13,16,17,18,19,20,23,43,68,93,118/)
    ndrn(1:13,19,4)=(/4,9,14,16,17,18,19,20,24,44,69,94,119/)
    ndrn(1:13,20,4)=(/5,10,15,16,17,18,19,20,25,45,70,95,120/)
    ndrn(1:13,21,4)=(/1,6,11,16,21,22,23,24,25,46,71,96,121/)
    ndrn(1:13,22,4)=(/2,7,12,17,21,22,23,24,25,47,72,97,122/)
    ndrn(1:13,23,4)=(/3,8,13,18,21,22,23,24,25,48,73,98,123/)
    ndrn(1:13,24,4)=(/4,9,14,19,21,22,23,24,25,49,74,99,124/)
    ndrn(1:13,25,4)=(/5,10,15,20,21,22,23,24,25,50,75,100,125/)
    ndrn(1:13,26,4)=(/1,26,27,28,29,30,31,36,41,46,51,76,101/)
    ndrn(1:13,27,4)=(/2,26,27,28,29,30,32,37,42,47,52,77,102/)
    ndrn(1:13,28,4)=(/3,26,27,28,29,30,33,38,43,48,53,78,103/)
    ndrn(1:13,29,4)=(/4,26,27,28,29,30,34,39,44,49,54,79,104/)
    ndrn(1:13,30,4)=(/5,26,27,28,29,30,35,40,45,50,55,80,105/)
    ndrn(1:13,31,4)=(/6,26,31,32,33,34,35,36,41,46,56,81,106/)
    ndrn(1:125,32,4)=(/(j,j=1,125)/)
    ndrn(1:125,33,4)=(/(j,j=1,125)/)
    ndrn(1:125,34,4)=(/(j,j=1,125)/)
    ndrn(1:13,35,4)=(/10,30,31,32,33,34,35,40,45,50,60,85,110/)
    ndrn(1:13,36,4)=(/11,26,31,36,37,38,39,40,41,46,61,86,111/)
    ndrn(1:125,37,4)=(/(j,j=1,125)/)
    ndrn(1:125,38,4)=(/(j,j=1,125)/)
    ndrn(1:125,39,4)=(/(j,j=1,125)/)
    ndrn(1:13,40,4)=(/15,30,35,36,37,38,39,40,45,50,65,90,115/)
    ndrn(1:13,41,4)=(/16,26,31,36,41,42,43,44,45,46,66,91,116/)
    ndrn(1:125,42,4)=(/(j,j=1,125)/)
    ndrn(1:125,43,4)=(/(j,j=1,125)/)
    ndrn(1:125,44,4)=(/(j,j=1,125)/)
    ndrn(1:13,45,4)=(/20,30,35,40,41,42,43,44,45,50,70,95,120/)
    ndrn(1:13,46,4)=(/21,26,31,36,41,46,47,48,49,50,71,96,121/)
    ndrn(1:13,47,4)=(/22,27,32,37,42,46,47,48,49,50,72,97,122/)
    ndrn(1:13,48,4)=(/23,28,33,38,43,46,47,48,49,50,73,98,123/)
    ndrn(1:13,49,4)=(/24,29,34,39,44,46,47,48,49,50,74,99,124/)
    ndrn(1:13,50,4)=(/25,30,35,40,45,46,47,48,49,50,75,100,125/)
    ndrn(1:13,51,4)=(/1,26,51,52,53,54,55,56,61,66,71,76,101/)
    ndrn(1:13,52,4)=(/2,27,51,52,53,54,55,57,62,67,72,77,102/)
    ndrn(1:13,53,4)=(/3,28,51,52,53,54,55,58,63,68,73,78,103/)
    ndrn(1:13,54,4)=(/4,29,51,52,53,54,55,59,64,69,74,79,104/)
    ndrn(1:13,55,4)=(/5,30,51,52,53,54,55,60,65,70,75,80,105/)
    ndrn(1:13,56,4)=(/6,31,51,56,57,58,59,60,61,66,71,81,106/)
    ndrn(1:125,57,4)=(/(j,j=1,125)/)
    ndrn(1:125,58,4)=(/(j,j=1,125)/)
    ndrn(1:125,59,4)=(/(j,j=1,125)/)
    ndrn(1:13,60,4)=(/10,35,55,56,57,58,59,60,65,70,75,85,110/)
    ndrn(1:13,61,4)=(/11,36,51,56,61,62,63,64,65,66,71,86,111/)
    ndrn(1:125,62,4)=(/(j,j=1,125)/)
    ndrn(1:125,63,4)=(/(j,j=1,125)/)
    ndrn(1:125,64,4)=(/(j,j=1,125)/)
    ndrn(1:13,65,4)=(/15,40,55,60,61,62,63,64,65,70,75,90,115/)
    ndrn(1:13,66,4)=(/16,41,51,56,61,66,67,68,69,70,71,91,116/)
    ndrn(1:125,67,4)=(/(j,j=1,125)/)
    ndrn(1:125,68,4)=(/(j,j=1,125)/)
    ndrn(1:125,69,4)=(/(j,j=1,125)/)
    ndrn(1:13,70,4)=(/20,45,55,60,65,66,67,68,69,70,75,95,120/)
    ndrn(1:13,71,4)=(/21,46,51,56,61,66,71,72,73,74,75,96,121/)
    ndrn(1:13,72,4)=(/22,47,52,57,62,67,71,72,73,74,75,97,122/)
    ndrn(1:13,73,4)=(/23,48,53,58,63,68,71,72,73,74,75,98,123/)
    ndrn(1:13,74,4)=(/24,49,54,59,64,69,71,72,73,74,75,99,124/)
    ndrn(1:13,75,4)=(/25,50,55,60,65,70,71,72,73,74,75,100,125/)
    ndrn(1:13,76,4)=(/01,26,51,76,77,78,79,80,81,86,91,96,101/)
    ndrn(1:13,77,4)=(/02,27,52,76,77,78,79,80,82,87,92,97,102/)
    ndrn(1:13,78,4)=(/03,28,53,76,77,78,79,80,83,88,93,98,103/)
    ndrn(1:13,79,4)=(/04,29,54,76,77,78,79,80,84,89,94,99,104/)
    ndrn(1:13,80,4)=(/05,30,55,76,77,78,79,80,85,90,95,100,105/)
    ndrn(1:13,81,4)=(/06,31,56,76,81,82,83,84,85,86,91,96,106/)
    ndrn(1:125,82,4)=(/(j,j=1,125)/)
    ndrn(1:125,83,4)=(/(j,j=1,125)/)
    ndrn(1:125,84,4)=(/(j,j=1,125)/)
    ndrn(1:13,85,4)=(/10,35,60,80,81,82,83,84,85,90,95,100,110/)
    ndrn(1:13,86,4)=(/11,36,61,76,81,86,87,88,89,90,91,96,111/)
    ndrn(1:125,87,4)=(/(j,j=1,125)/)
    ndrn(1:125,88,4)=(/(j,j=1,125)/)
    ndrn(1:125,89,4)=(/(j,j=1,125)/)
    ndrn(1:13,90,4)=(/15,40,65,80,85,86,87,88,89,90,95,100,115/)
    ndrn(1:13,91,4)=(/16,41,66,76,81,86,91,92,93,94,95,96,116/)
    ndrn(1:125,92,4)=(/(j,j=1,125)/)
    ndrn(1:125,93,4)=(/(j,j=1,125)/)
    ndrn(1:125,94,4)=(/(j,j=1,125)/)
    ndrn(1:13,95,4)=(/20,45,70,80,85,90,91,92,93,94,95,100,120/)
    ndrn(1:13,96,4)=(/21,46,71,76,81,86,91,96,97,98,99,100,121/)
    ndrn(1:13,97,4)=(/22,47,72,77,82,87,92,96,97,98,99,100,122/)
    ndrn(1:13,98,4)=(/23,48,73,78,83,88,93,96,97,98,99,100,123/)
    ndrn(1:13,99,4)=(/24,49,74,79,84,89,94,96,97,98,99,100,124/)
    ndrn(1:13,100,4)=(/25,50,75,80,85,90,95,96,97,98,99,100,125/)
    ndrn(1:13,101,4)=(/1,26,51,76,101,102,103,104,105,106,111,116,121/)
    ndrn(1:13,102,4)=(/2,27,52,77,101,102,103,104,105,107,112,117,122/)
    ndrn(1:13,103,4)=(/3,28,53,78,101,102,103,104,105,108,113,118,123/)
    ndrn(1:13,104,4)=(/4,29,54,79,101,102,103,104,105,109,114,119,124/)
    ndrn(1:13,105,4)=(/5,30,55,80,101,102,103,104,105,110,115,120,125/)
    ndrn(1:13,106,4)=(/6,31,56,81,101,106,107,108,109,110,111,116,121/)
    ndrn(1:13,107,4)=(/7,32,57,82,102,106,107,108,109,110,112,117,122/)
    ndrn(1:13,108,4)=(/8,33,58,83,103,106,107,108,109,110,113,118,123/)
    ndrn(1:13,109,4)=(/9,34,59,84,104,106,107,108,109,110,114,119,124/)
    ndrn(1:13,110,4)=(/10,35,60,85,105,106,107,108,109,110,115,120,125/)
    ndrn(1:13,111,4)=(/11,36,61,86,101,106,111,112,113,114,115,116,121/)
    ndrn(1:13,112,4)=(/12,37,62,87,102,107,111,112,113,114,115,117,122/)
    ndrn(1:13,113,4)=(/13,38,63,88,103,108,111,112,113,114,115,118,123/)
    ndrn(1:13,114,4)=(/14,39,64,89,104,109,111,112,113,114,115,119,124/)
    ndrn(1:13,115,4)=(/15,40,65,90,105,110,111,112,113,114,115,120,125/)
    ndrn(1:13,116,4)=(/16,41,66,91,101,106,111,116,117,118,119,120,121/)
    ndrn(1:13,117,4)=(/17,42,67,92,102,107,112,116,117,118,119,120,122/)
    ndrn(1:13,118,4)=(/18,43,68,93,103,108,113,116,117,118,119,120,123/)
    ndrn(1:13,119,4)=(/19,44,69,94,104,109,114,116,117,118,119,120,124/)
    ndrn(1:13,120,4)=(/20,45,70,95,105,110,115,116,117,118,119,120,125/)
    ndrn(1:13,121,4)=(/21,46,71,96,101,106,111,116,121,122,123,124,125/)
    ndrn(1:13,122,4)=(/22,47,72,97,102,107,112,117,121,122,123,124,125/)
    ndrn(1:13,123,4)=(/23,48,73,98,103,108,113,118,121,122,123,124,125/)
    ndrn(1:13,124,4)=(/24,49,74,99,104,109,114,119,121,122,123,124,125/)
    ndrn(1:13,125,4)=(/25,50,75,100,105,110,115,120,121,122,123,124,125/)
    ! node=11
    do i=1,11
      ndrn(1:11,i,5)=(/(j,j=1,11)/)
    enddo
    
  elseif(ndim.eq.2) then ! 2d
    ndtok((/7,8,9,16,25/))=(/5,1,2,3,4/) ! for 2d, used to replace the number of nodes in corresponding element
    
    ! nface is the number of the faces in an element
    nface(1:4)=4
    nface(5)=3
    
    ! npface is the number of the points in an element face
    npface(1)=3
    npface(2)=3
    npface(3)=4
    npface(4)=5
    npface(5)=3
    
    ! ndface is the node number of each face in the element
    ! ndface(node in this face, face, element ndtok)
    ! node.eq.8
    ndface(1:3,1:4,1)=reshape((/1,5,2, 2,6,3, 3,7,4, 4,8,1/),(/3,4/))
    ! node.eq.9
    ndface(1:3,1:4,2)=reshape((/1,2,3, 3,6,9, 9,8,7, 7,4,1/),(/3,4/))
    ! node.eq.16
    ndface(1:4,1:4,3)=reshape((/1,2,3,4, 4,8,12,16, 16,15,14,13, 13,9,5,1/),(/4,4/))
    ! node.eq.25
    ndface(1:5,1:4,4)=reshape((/1,2,3,4,5, 5,10,15,20,25,25,24,23,22,21, 21,16,11,6,1/),(/5,4/))
    ! node.eq.7
    ndface(1:3,1:3,5)=reshape((/1,4,2, 2,5,3, 3,6,1/),(/3,3/))
    
    ! ndtype is the type of nodes in the element
    ! ndtype(node, element ndtok)
    !        = 1, element internal nodes
    !        = 2, element surface nodes
    ! node.eq.8
    ndtype(1:8,1)=2
    ! node.eq.9
    ndtype(1:9,2)=2
    ndtype(5,2)=1
    ! node.eq.16
    ndtype(1:16,3)=(/2,2,2,2, 2,1,1,2, 2,1,1,2, 2,2,2,2/)
    ! node.eq.25
    ndtype(1:25,4)=(/2,2,2,2,2, 2,1,1,1,2, 2,1,1,1,2, 2,1,1,1,2,2,2,2,2,2/)
    ! node.eq.7
    ndtype(1:6,5)=2
    ndtype(7,5)=1
    
    ! nrn is the number of nodes that are related to the node in the element
    ! nrn(node, element ndtok)
    ! node=8
    nrn(1:8,1)=8
    ! node=9
    nrn(1:9,2)=5
    nrn(5,2)=9
    ! node=16
    nrn(1:16,3)=7
    nrn((/6,7,10,11/),3)=16
    ! node=25
    nrn(1:125,4)=9
    nrn((/7,8,9,12,13,14,17,18,19/),4)=25
    ! node=7
    nrn(1:7,5)=7
    
    ! ndrn is the node number that are related to the node in the element
    ! ndrn(node, related node, element ndtok)= related node number in the element
    ! node=8
    do i=1,8
      ndrn(1:8,i,1)=(/(j,j=1,8)/)
    enddo
    ! node=9
    ndrn(1:5,1,2)=(/1,2,3,4,7/)
    ndrn(1:5,2,2)=(/1,2,3,5,8/)
    ndrn(1:5,3,2)=(/1,2,3,6,9/)
    ndrn(1:5,4,2)=(/1,4,5,6,7/)
    ndrn(1:9,5,2)=(/(j,j=1,9)/)
    ndrn(1:5,6,2)=(/3,4,5,6,9/)
    ndrn(1:5,7,2)=(/1,4,7,8,9/)
    ndrn(1:5,8,2)=(/2,5,7,8,9/)
    ndrn(1:5,9,2)=(/3,6,7,8,9/)
    ! node=16
    ndrn(1:7,1,3)=(/1,2,3,4,5,9,13/)
    ndrn(1:7,2,3)=(/1,2,3,4,6,10,14/)
    ndrn(1:7,3,3)=(/1,2,3,4,7,11,15/)
    ndrn(1:7,4,3)=(/1,2,3,4,8,12,16/)
    ndrn(1:7,5,3)=(/1,5,6,7,8,9,13/)
    ndrn(1:16,6,3)=(/(j,j=1,16)/)
    ndrn(1:16,7,3)=(/(j,j=1,16)/)
    ndrn(1:7,8,3)=(/4,5,6,7,8,12,16/)
    ndrn(1:7,9,3)=(/1,5,9,10,11,12,13/)
    ndrn(1:16,10,3)=(/(j,j=1,16)/)
    ndrn(1:16,11,3)=(/(j,j=1,16)/)
    ndrn(1:7,12,3)=(/4,8,9,10,11,12,16/)
    ndrn(1:7,13,3)=(/1,5,9,13,14,15,16/)
    ndrn(1:7,14,3)=(/2,6,10,13,14,15,16/)
    ndrn(1:7,15,3)=(/3,7,11,13,14,15,16/)
    ndrn(1:7,16,3)=(/4,8,12,13,14,15,16/)
    ! node=25
    ndrn(1:9,1,4)=(/1,2,3,4,5,6,11,16,21/)
    ndrn(1:9,2,4)=(/1,2,3,4,5,7,12,17,22/)
    ndrn(1:9,3,4)=(/1,2,3,4,5,8,13,18,23/)
    ndrn(1:9,4,4)=(/1,2,3,4,5,9,14,19,24/)
    ndrn(1:9,5,4)=(/1,2,3,4,5,10,15,20,25/)
    ndrn(1:9,6,4)=(/1,6,7,8,9,10,11,16,21/)
    ndrn(1:25,7,4)=(/(j,j=1,25)/)
    ndrn(1:25,8,4)=(/(j,j=1,25)/)
    ndrn(1:25,9,4)=(/(j,j=1,25)/)
    ndrn(1:9,10,4)=(/5,6,7,8,9,10,15,20,25/)
    ndrn(1:9,11,4)=(/1,6,11,12,13,14,15,16,21/)
    ndrn(1:25,12,4)=(/(j,j=1,25)/)
    ndrn(1:25,13,4)=(/(j,j=1,25)/)
    ndrn(1:25,14,4)=(/(j,j=1,25)/)
    ndrn(1:9,15,4)=(/5,10,11,12,13,14,15,20,25/)
    ndrn(1:9,16,4)=(/1,6,11,16,17,18,19,20,21/)
    ndrn(1:25,17,4)=(/(j,j=1,25)/)
    ndrn(1:25,18,4)=(/(j,j=1,25)/)
    ndrn(1:25,19,4)=(/(j,j=1,25)/)
    ndrn(1:9,20,4)=(/5,10,15,16,17,18,19,20,25/)
    ndrn(1:9,21,4)=(/1,6,11,16,21,22,23,24,25/)
    ndrn(1:9,22,4)=(/2,7,12,17,21,22,23,24,25/)
    ndrn(1:9,23,4)=(/3,8,13,18,21,22,23,24,25/)
    ndrn(1:9,24,4)=(/4,9,14,19,21,22,23,24,25/)
    ndrn(1:9,25,4)=(/5,10,15,20,21,22,23,24,25/)
    ! node=7
    do i=1,7
      ndrn(1:7,i,5)=(/(j,j=1,7)/)
    enddo
    
  !else ! 1d
  !  ndtok((/3,4,5/))=(/1,2,3/) ! for 1d
  ! 
  !  ! nface is the number of the faces in an element
  !  nface(1:3)=2
  !  
  !  ! npface is the number of the points in an element face
  !  npface=1
  !
  !  ! node.eq.3
  !  ndtype(1:3,1)=(/2,1,2/)
  !  ! node.eq.4
  !  ndtype(1:4,2)=(/2,1,1,2/)
  !  ! node.eq.5
  !  ndtype(1:5,3)=(/2,1,1,1,2/)
  endif
  return
end subroutine netypend