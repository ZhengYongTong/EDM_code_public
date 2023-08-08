!---------------------------------------------------------------------------------------------------------------------------------------------
! vary_arrays: Global arrays in program EDM
!---------------------------------------------------------------------------------------------------------------------------------------------
module vary_arrays
  implicit real*8 (a-h,o-z)
  integer,dimension (:), allocatable :: netype,lret,nnode,nematg,kbdelem,neop,neoe,lrem,mbody,ibody,lnof,isfirtyp
  integer,dimension (:,:), allocatable ::nbcgrp,lnde
  integer,dimension (:,:,:), allocatable ::kbcflag
  integer,dimension (:), allocatable :: ne0,npr,nct,nrh,ncp,metype,mptype,mcttype,mrtype,mcptype
  integer,dimension (:), allocatable :: ia,ja,ic,jc,ir,jr,irm,jrm,it,jt,jtp,itef
  real*8,dimension (:), allocatable :: xval,b,b0,a,c,r,rm,t,cdyn,xvalp,tef
  real*8,dimension (:,:), allocatable :: ce,cpr,cct,crh,ccp
  real*8,dimension (:,:), allocatable :: cd,vbody,preu,prerho,precp,pret
  real*8,dimension (:,:,:), allocatable :: prect0
  real*8,dimension (:,:,:,:), allocatable :: pred,dndxa
  real*8,dimension (:,:,:,:,:), allocatable :: dndxxa
end module vary_arrays

!---------------------------------------------------------------------------------------------------------------------------------------------
! fixed_values: Global fixed values in program EDM
!---------------------------------------------------------------------------------------------------------------------------------------------
module fixed_values
  implicit real*8 (a-h,o-z)
  integer :: nsig,ndim,ndf,ntp,nte,mnode,ntf,ntf1,msolver,nfacem,nbody,nmatgrp,nsigp,istrans,isnonli,istm,mdyn,nobdf,noidf
  integer :: ka,kc,kr,krm,kt
  real*8 :: fjcb,pi,theta,deltat,tolt,toln,omg
  integer :: maxhb
  real*8 :: avd,avl
  dimension dlt(3,3),iny(9),ndtok(125),ndface(25,6,5),nface(5),ndtype(125,5),npface(5),nrn(125,5),ndrn(125,125,5)
  character*60 timew
  character*10 datev,timev
  character*1 yorn
end module fixed_values

!---------------------------------------------------------------------------------------------------------------------------------------------
module compute
  ! Module containing:
  ! Subroutine addmkl: add two CSR format matrixes.
  ! Subroutine matmulmkl_mmm: multiply two CSR format matrixes.
  ! Subroutine pardiso_inv: obtain the inverse of a CSR format matrix.
  ! Subroutine pardiso_inv_storless: obtain the inverse of a CSR format matrix (in case the insufficient memory).
  
  contains
  subroutine addmkl(m,n,ia,ja,a,ka,ib,jb,b,kb)
    implicit none
    integer::m,n,ka,kb,kc,info,nzmax
    integer::ia(m+1),ib(m+1),jb(kb),ic(m+1)
    real*8::beta
    real*8::b(kb)
    integer,allocatable::jc(:),ja(:)
    real*8,allocatable::c(:),a(:)
    
    beta=1.d0
    call mkl_dcsradd('N',1,3,m,n,a,ja,ia,beta,b,jb,ib,c,jc,ic,nzmax,info)
    kc=ic(m+1)-1
    allocate(jc(kc),c(kc))
    nzmax=kc
    call mkl_dcsradd('N',2,3,m,n,a,ja,ia,beta,b,jb,ib,c,jc,ic,nzmax,info)
    if(info.ne.0) then
      write(*,*) 'Error in calculating CSR matrix adding!'
      write(4,*) 'Error in calculating CSR matrix adding!'
      stop
    endif
    deallocate(ja,a)
    ka=kc
    allocate(ja(ka),a(ka))
    ia=ic
    ja=jc
    a=c
  end subroutine addmkl
  
  subroutine matmulmkl_mmm(m,n,k,ia,ja,a,ka,ib,jb,b,kb,ic,jc,c,kc)
    implicit none
    integer::m,n,k,ka,kb,kc,nzmax,info
    integer::ia(m+1),ja(ka),ib(n+1),jb(kb),ic(m+1)
    real*8::a(ka),b(kb)
    integer,allocatable::jc(:)
    real*8,allocatable::c(:)
    
    call mkl_dcsrmultcsr('N',1,3,m,n,k,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,info)
    kc=ic(m+1)-1
    allocate(jc(kc),c(kc))
    nzmax=kc
    call mkl_dcsrmultcsr('N',0,3,m,n,k,a,ja,ia,b,jb,ib,c,jc,ic,nzmax,info)
    if(info.ne.0) then
      write(*,*) 'Error in calculating CSR matrix multiplying!'
      write(4,*) 'Error in calculating CSR matrix multiplying!'
      stop
    endif
  end subroutine matmulmkl_mmm
  
  subroutine pardiso_inv(n,ia,ja,a,ka,m,ib,jb,b,kb,ix,jx,x,kx,c,y,error)
    implicit none
    integer,allocatable::ia(:),ja(:),ib(:),jb(:),ix(:),jx(:)
    real*8,allocatable::a(:),b(:),x(:),c(:),y(:),bp(:,:),xp(:,:)
    integer*8 pt(64)
    integer maxfct,mnum,mtype,phase,nrhs,error,msglvl
    integer n,m,ka,kb,kx
    integer iparm(64),job(8)
    integer idum(1)
    integer i,j
    real*8  ddum(1)
    
    data nrhs /1/, maxfct /1/, mnum /1/
    
    ! set up pardiso control parameter
    mtype = 11 ! real unsymmetric
    error = 0 ! initialize error flag
    msglvl = 0 ! print statistical information
    call pardisoinit (pt, mtype, iparm) ! initiliaze pardiso solver
    iparm(1) = 1 ! no solver default
    iparm(2) = 3 ! fill-in reordering from metis
    iparm(8) = 9 ! numbers of iterative refinement steps
    
    ! reordering and symbolic factorization, this step also allocates
    ! all memory that is necessary for the factorization
    phase = 11 ! only reordering and symbolic factorization
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),1)
      stop
    endif
    
    ! factorization.
    phase = 22 ! only factorization
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),2)
      stop
    endif
    
    ! back substitution and iterative refinement
    iparm(8) = 2 ! max numbers of iterative refinement steps
    phase = 33 ! only solve
    ! transform b into dense format
    allocate(bp(n,m))
    bp=0.d0
    do i=1,n
      do j=ib(i),ib(i+1)-1
        bp(i,jb(j))=b(j)
      enddo
    enddo
    ! solve ax=b
    allocate(xp(n,m))
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, m, iparm, msglvl, bp, xp, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),3)
      stop
    endif
    deallocate(bp)
    ! transform x into csr format
    allocate(ix(n+1))
    ix(1)=1
    kx=0
    do i=1,n
      do j=1,m
        if(dabs(xp(i,j)).gt.1.d-15)then
          kx=kx+1
        endif
      enddo
      ix(i+1)=kx+1
    enddo
    allocate(jx(kx),x(kx))
    kx=0
    do i=1,n
      do j=1,m
        if(dabs(xp(i,j)).gt.1.d-15)then
          kx=kx+1
          jx(kx)=j
          x(kx)=xp(i,j)
        endif
      enddo
    enddo
    deallocate(xp)
    ! solve ay=c
    allocate(y(n))
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, nrhs, iparm, msglvl, c, y, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),3)
      stop
    endif
    
    ! termination and release of memory
    phase = -1 ! release internal memory
    call pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, idum, nrhs, iparm, msglvl, ddum, ddum, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),4)
      stop
    endif
    
  end subroutine pardiso_inv
  
  subroutine pardiso_inv_storless(n,ia,ja,a,ka,m,ib,jb,b,kb,ix,jx,x,kx,c,y,error)
    ! in case the insufficient memory, notice that b matrix will be transposed after calculating
    implicit none
    integer,allocatable::ia(:),ja(:),ib(:),jb(:),ix(:),ix1(:),jx(:),ixp(:),jxp(:),ib1(:),jb1(:),ixp2(:),ibf(:),jbf(:),ibp(:),jbp(:),kxp(:)
    real*8,allocatable::a(:),b(:),x(:),xp(:),b1(:),c(:),y(:),d(:,:),e(:,:),bf(:),bp(:)
    integer*8 pt(64)
    integer maxfct,mnum,mtype,phase,nrhs,nrhsp,knrhs,error,msglvl
    integer n,m,ka,kb,kx
    integer iparm(64),job(8)
    integer idum(1)
    integer i,j,l,info
    real*8  ddum(1)
    character*1 yorn
    
    ! transpose b for easily reading b by each columns
    job(1)=0;job(2)=1;job(3)=1;job(6)=1
    if(n.ge.m) then
      allocate(ib1(n+1),jb1(kb),b1(kb))
      call mkl_dcsrcsc(job,n,b,jb,ib,b1,jb1,ib1,info)
      b=b1
      jb=jb1
      deallocate(ib)
      allocate(ib(m+1))
      ib(1:m+1)=ib1(1:m+1)
      deallocate(ib1,jb1,b1)
    else ! n.lt.m
      allocate(ib1(m+1),jb1(kb),b1(kb))
      ib1(1:n+1)=ib(1:n+1)
      ib1(n+2:m+1)=ib(n+1) ! fill row with 0
      jb1=jb
      b1=b
      deallocate(ib)
      allocate(ib(m+1))
      call mkl_dcsrcsc(job,m,b1,jb1,ib1,b,jb,ib,info)
      deallocate(ib1,jb1,b1)
    endif
    
    data maxfct /1/, mnum /1/
    nrhs=1
    ! set up pardiso control parameter
    mtype = 11 ! real unsymmetric
    error = 0 ! initialize error flag
    msglvl = 0 ! print statistical information
    call pardisoinit (pt, mtype, iparm) ! initiliaze pardiso solver
    iparm(1) = 1 ! no solver default
    iparm(2) = 3 ! fill-in reordering from metis
    iparm(8) = 9 ! numbers of iterative refinement steps
    
    phase = 11 ! only reordering and symbolic factorization
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),1)
      return
    endif
    
    phase = 22 ! only factorization
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, nrhs, iparm, msglvl, ddum, ddum, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),2)
      return
    endif
    
    iparm(8) = 2 ! max numbers of iterative refinement steps
    phase = 33 ! only factorization
    
    open(100,file='inv.qqq',status='unknown',form='unformatted')
    nrhsp=floor(10000000000./(8*n)) ! the max number of columns of d matrix, 10gb ram space for d
    if(nrhsp.eq.0)then
        write(*,*)'Matrix to be inversed is too large'
        write(*,*)'Program need to be changed'
        write(*,*)'Do you want to continue? (y or n)'
        read(*,*)yorn
        if(yorn.ne.'y'.and.yorn.ne.'y')then
           stop
        endif
    endif
    knrhs=floor(dble(m)/nrhsp)
    allocate(kxp(knrhs+1))
    kxp=0
    if(knrhs.eq.0) goto 10
    nrhs=nrhsp
    allocate(d(n,nrhs),e(n,nrhs))
    d=0.d0
    do i=1,knrhs
      ! get d from b
      do l=(i-1)*nrhsp+1,i*nrhsp
        do j=ib(l),ib(l+1)-1
          d(jb(j),l-(i-1)*nrhsp)=b(j)
        enddo
      enddo
      ! solve a*e=d
      call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,idum, nrhs, iparm, msglvl, d, e, error)
      if(error.ne.0)then
        call par_er(error,iparm(15),iparm(16),iparm(17),3)
        return
      endif
      ! store result e into file 100
      allocate(ix(nrhs+1))
      ix(1)=1
      do l=1,nrhs
        do j=1,n
          if(abs(e(j,l)).gt.1.0d-15) then
            kxp(i)=kxp(i)+1
          endif
        enddo
        ix(l+1)=1+kxp(i)
      enddo
      kx=0
      allocate(jx(kxp(i)),x(kxp(i)))
      do l=1,nrhs
        do j=1,n
          if(abs(e(j,l)).gt.1.0d-15) then
            kx=kx+1
            jx(kx)=j
            x(kx)=e(j,l)
          endif
        enddo
      enddo
      write(100)ix(1:nrhs+1)
      write(100)jx
      write(100)x
      deallocate(ix,jx,x)
      ! set d to be zero
       do l=(i-1)*nrhsp+1,i*nrhsp
         do j=ib(l),ib(l+1)-1
           d(jb(j),l-(i-1)*nrhsp)=0
         enddo
       enddo
    enddo
    deallocate(d,e)
     10 continue
    ! deal with rest m-knrhs*nrhsp columns
    nrhs=m-knrhs*nrhsp
    if(nrhs.ne.0) then
      allocate(d(n,nrhs),e(n,nrhs))
      ! get d from b
      do l=knrhs*nrhsp+1,m
        do j=ib(l),ib(l+1)-1
          d(jb(j),l-knrhs*nrhsp)=b(j)
        enddo
      enddo
      ! solve a*e=d
      call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja,idum, nrhs, iparm, msglvl, d, e, error)
      if(error.ne.0)then
        call par_er(error,iparm(15),iparm(16),iparm(17),3)
        return
      endif
      deallocate(d)
      allocate(ix1(nrhs+1))
      ix1(1)=1
      do l=1,nrhs
        do j=1,n
          if(abs(e(j,l)).gt.1.0d-15) then
            kxp(knrhs+1)=kxp(knrhs+1)+1
          endif
        enddo
        ix1(l+1)=1+kxp(knrhs+1)
      enddo
    endif
    ! solve a*y=c
    allocate(y(n))
    call pardiso (pt, maxfct, mnum, mtype, phase, n, a, ia, ja, idum, 1, iparm, msglvl, c, y, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),3)
      return
    endif
    ! release internal memory
    phase = -1
    call pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum,idum, idum, 1, iparm, msglvl, ddum, ddum, error)
    if(error.ne.0)then
      call par_er(error,iparm(15),iparm(16),iparm(17),4)
      return
    endif
    
    ! assembling each matrix
    kx=sum(kxp)
    if(n.ge.m)then
      allocate(ixp(n+1),jxp(kx),xp(kx))
    else
      allocate(ixp(m+1),jxp(kx),xp(kx))
    endif
    ixp(1)=1
    rewind(100)
    do i=1,knrhs
      allocate(ix(nrhsp+1),jx(kxp(i)),x(kxp(i)))
      read(100)ix
      read(100)jx
      read(100)x
      do l=(i-1)*nrhsp+1,i*nrhsp
        ixp(l+1)=ix(l-(i-1)*nrhsp+1)+ixp((i-1)*nrhsp+1)-1
      enddo
      do l=1,ix(nrhsp+1)-1
        jxp(ixp((i-1)*nrhsp+1)-1+l)=jx(l)
        xp(ixp((i-1)*nrhsp+1)-1+l)=x(l)
      enddo
      deallocate(ix,jx,x)
    enddo
    close(100,status='delete')
    
    ! deal with rest m-knrhs*nrhsp columns
    if(nrhs.ne.0)then
      kx=0
      allocate(jx(kxp(knrhs+1)),x(kxp(knrhs+1)))
      do l=1,nrhs
        do j=1,n
          if(abs(e(j,l)).gt.1.0d-15) then
            kx=kx+1
            jx(kx)=j
            x(kx)=e(j,l)
          endif
        enddo
      enddo
      deallocate(e)
      do l=knrhs*nrhsp+1,m
        ixp(l+1)=ix1(l-knrhs*nrhsp+1)+ixp(knrhs*nrhsp+1)-1
      enddo
      do l=1,ix1(m-knrhs*nrhsp+1)-1
        jxp(ixp(knrhs*nrhsp+1)-1+l)=jx(l)
        xp(ixp(knrhs*nrhsp+1)-1+l)=x(l)
      enddo
      deallocate(ix1,jx,x)
    endif
    deallocate(kxp)
    
    ! transpose x
    if(n.ge.m) then
      ixp(m+2:n+1)=ixp(m+1)
      allocate(ix(n+1),jx(kx),x(kx))
      job(1)=0;job(2)=1;job(3)=1;job(6)=1
      call mkl_dcsrcsc(job,n,xp,jxp,ixp,x,jx,ix,info)
      deallocate(ixp,jxp,xp)
    else ! n.lt.m
      allocate(ixp2(m+1),jx(kx),x(kx))
      job(1)=0;job(2)=1;job(3)=1;job(6)=1
      call mkl_dcsrcsc(job,m,xp,jxp,ixp,x,jx,ixp2,info)
      deallocate(ixp,jxp,xp)
      allocate(ix(n+1))
      ix(1:n+1)=ixp2(1:n+1)
      deallocate(ixp2)
    endif
    
  end subroutine pardiso_inv_storless
end module compute
!---------------------------------------------------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------------------------------------------------
! CMD_Progress: for printing the progress bar
!---------------------------------------------------------------------------------------------------------------------------------------------
Module CMD_Progress
  Implicit None
  private
  Logical , parameter , public :: CMD_PROGRESS_ABSOLUTE = .true.
  Type , public :: CLS_CMD_Progress
    Integer , private :: N , lens , i
    Character :: M = "*" , O = "."
    Character(len=64) :: Prefix
  Contains
    Procedure :: Set
    Procedure :: Put
  End Type CLS_CMD_Progress
   
contains
 
  Subroutine Set( this , N , L )
    Class( CLS_CMD_Progress ) :: this
    Integer , Intent( IN ) :: N , L
    this % N    = N
    this % lens = L
    this % i = 0
    this % Prefix = " Progress: "
  End Subroutine Set
   
  Subroutine Put( this , K , bAbsol )
    Class( CLS_CMD_Progress ) :: this
    Integer , Intent( IN ) :: K
    Logical , optional :: bAbsol
    Character(len=1) :: br
    integer :: jm
    this % i = this % i + K
    if ( present( bAbsol ) ) then
      if ( bAbsol ) this % i = K
    end if
    if ( this % i > this % n ) this % i = this % n
    jm = Nint( real( this%i * this%lens ) / real( this%N ) )
    if ( this%i < this%n ) then
      br = char(13)
    else
      br = char(10)
    end if
    write( * , '(5a,f6.2,2a\)') trim(this%Prefix) , '[' , & ! If your compiler do not support, replace this line by the next line
    !write( * , '(5a,f6.2,2a)',advance="no") trim(this%Prefix) , '[' , &
      repeat(this%M , jm ) , repeat( this%O , this%lens-jm ) , '] ' , this%i*100.0/this%N , "%" , br
  End Subroutine Put
   
End Module CMD_Progress