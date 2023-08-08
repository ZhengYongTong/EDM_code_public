!#########################################################################################
  ! This is a basic and general program for element differential method (EDM), a numerical algorithm proposed by Xiao-Wei Gao et al.
  ! The program is developed by Xiao-Wei Gao, Yong-Tong Zheng et al.
  
  ! The program can solve linear steady heat transfer problems (when ndf=1, istrans=0, isnonli=0),
  !                                         linear elasticity problems (when ndf=2 or 3, istrans=0, isnonli=0),
  !                                         linear transient heat transfer problems (when ndf=1, istrans=1, isnonli=0),
  !                                         linear elastodynamic response problems (when ndf=2 or 3, istrans=1, isnonli=0),
  !                                         non-linear steady heat transfer problems (when ndf=1, istrans=0, isnonli=1),
  !                                  and non-linear transient heat transfer problems (when ndf=1, istrans=1, isnonli=1).
  
  ! Before using the program:
  ! 1. Some basic knowledge (governing equation and boundary condition etc.) about heat conduction and mechanics should be known.
  ! 2. A FORTRAN compiler should be installed in the computer and it had better to be Intel Fortran or One API, otherwise the linear algebraic solver PARDISO, DSS,
  ! FGMRES_ILU0, FGMRES_ILUT and some other subroutines about matrix adding or multiplying may not work. It is highly recommended to use:
  ! Microsoft Visual Studio Community (MVS) (https://visualstudio.microsoft.com/downloads),
  ! Intel OneAPI Base Toolkit and HPC Toolkit (https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit) or Intel visual FORTRAN (IVF).
  ! If using these two softwares above, to compile the codes into a executable program, one can use the following compile command:
  ! IFORT /Qmkl:parallel /exe:EDM.exe *.F *.f90
  ! If compiling them in IDE, one still need to choose 'Use Intel Math Kernel Library' in 'Property' (not to be 'No').
  ! If one want to use some other compilers, It is recommended to find some other sparse matrix solvers.
  ! 3. It is highly recommended to read the EDM_User_Manual.doc to learn how to make a input file (EDM.inp) and what files will be generated after running the program.
  
  ! Before reading the code, it is highly recommended to refer to the following five papers:
  !     International Journal of Heat and Mass Transfer 115 (2017) 882每894 ----- linear steady heat transfer problems
  !     International Journal for Numerical Methods in Engineering vol.: 113, issue: 1,(2018) 82每108 ----- linear elasticity problems
  !     International Journal of Heat and Mass Transfer 127 (2018) 1189每1197 ----- linear transient heat transfer problems
  !     International Journal of Mechanical Sciences 151 (2019) 828每841 ----- linear elastodynamic response problems
  !     International Journal of Heat and Mass Transfer 126 (2018) 1111每1119 ----- non-linear steady and transient heat transfer problems
  
  ! If you want to get more details or give some suggestions, please contact prof. Gao via email: xwgao@dlut.edu.cn
!#########################################################################################

!****************************************************************************
  ! This file containing:
  ! FORTRAN Main program EDM
  ! Subroutine step_info: for generate the information (including residual, step number etc.) of each time step or non-linear iterative step
  ! Subroutine time_record: for output the time used in each part
  ! Function str2num: converts a string to a integer
!****************************************************************************

program EDM
  use vary_arrays
  use fixed_values
  implicit real*8 (a-h,o-z)
  allocatable xval1(:),tval(:),conve(:),xval0(:),er(:)
  character*80 title,t_name
  character*10 filenm
  character*10 date0,time0
  
  open(4,file='EDM.LOG',status='unknown') ! computational time recording
  open(5,file='EDM.INP',status='old') ! input file
  open(7,file='EDM.OUT',status='unknown') ! results' output file
  
  ! channels for plotting
  open(32,file='INPUT_GEO.NAS',status='unknown') ! plot input mesh to a NAS file
  open(33,file='RESULTS.PLT',status='unknown',form='binary') ! file for plotting contours in tecplot
  
  ! time record beggining
  call date_and_time (datev,timev) 
  time0=timev
  date0=datev
  
  ! read control variables
  call input_ctr
  
  ! assign values for commonly used arrays
  call block_data 
  
  allocate (cd(ndim,ntp),lnde(mnode,nte),nbcgrp(nte,nfacem),kbcflag(ndf,nte,nfacem),nnode(nte),netype(nte),nematg(nte))
  
  ! read co-ordinates, boundary condition et al.
  call input_data

  ! file for thermal stress analysis (temperatures of all nodes)
  if(ndf.eq.1) then ! for steady heat transfer problem, generate a temperature.dat file for possible next step thermal stress analysis
    open(45,file='TEMPERATURE.DAT',status='unknown')
  elseif(istm.eq.1) then ! thermal stress analysis needs a temperature.dat file
    open(45,file='TEMPERATURE.DAT',status='old')
  endif
  
  ! plot geometry
  ifile=32
  call nastran_mesh(ifile,ndim,ndim,ntp,nte,cd,lnde,nnode,mnode,netype)
  close(ifile,status='keep')
  timew='Generate femap plot of boundary :'
  call time_record(timew,datev,timev,4,2)
  
  ! get relationship of nodes and elements, generate ia and ja
  maxe=25 ! the maximum number of collocation elements for each node
  maxn=125 ! the maximum number of collocation nodes for each node
  call node_sort(maxe,maxn)
  timew='Get relationship of nodes and elements :'
  call time_record(timew,datev,timev,4,2)
  
  call gen_dndx_dndxx
  timew='Generate derivates of each element :'
  call time_record(timew,datev,timev,4,2)

  if(istrans.eq.0.and.isnonli.eq.0)then ! for steady heat transfer and elaticity problem
    write(*,*)
    write(*,*)'           Forming the coefficient matrix and the right-hand-side vector'
    allocate (a(ka),b(ntf),isfirtyp(ntf))
    a=0.d0
    b=0.d0
    isfirtyp=0
    call coef_gen
    call first_type_bc(-1)
    call second_type_bc(-1)
    call third_type_bc(-1)
    call source_item(-1)
    deallocate (kbdelem,vbody,neoe,pret,lret,tef,itef)
    timew='Calculating coefficients :'
    call time_record(timew,datev,timev,4,2)

    write(*,*)
    write(*,*)'           Solving the system of equations'
    allocate(xval(ntf))
    call solve_info
    call solve_system
    deallocate (ia,ja,a,b)
    timew='Solving the system of equations :'
    call time_record(timew,datev,timev,4,2)
    
    write(*,*)
    write(*,*)'          Outputing results'
    call output
    timew='Output results :'
    call time_record(timew,datev,timev,4,2)
    
  elseif(istrans.eq.1.and.isnonli.eq.0)then ! transient heat transfer problem
    
    if(ndf.eq.1)then ! transient heat
      call solve_output_for_trans_heat
    else ! dynamics
      call solve_output_for_dyn
    endif
    timew='Solve and output each step:'
    call time_record(timew,datev,timev,4,2)
    
  elseif(istrans.eq.0.and.isnonli.eq.1)then ! material non-linear heat transfer problem
    
    write(*,*)
    write(*,*)'          Non-linear iteration'
    write(4,*)
    write(4,*)'          Non-linear iteration'
    allocate (a(ka),b(ntf),isfirtyp(ntf),tval(ntp),xval1(ntf),conve(ntf),er(ndf),b0(ntf))
    if(ndf.ne.1.and.istm.eq.1) then
      ! read nodal temperatures
      rewind (45)
      read(45,*) ktp
      do ik=1,min(ktp,ntp)
        read(45,*)ip,tval(ip)
      enddo
    endif
    
    write(4,*)
    do istep=1,50
      xval1=xval ! record the results of last step for comparision
      a=0.d0
      b=0.d0
      isfirtyp=0
      call coef_nonli(istep,tval)
      call first_type_bc(istep)
      call second_type_bc(istep)
      call third_type_bc(istep)
      call source_item(istep)
      ! get b0=a(x)*x
      call matmulmkl_mmk(ntf,ntf,ia,ja,a,ka,xval,b0,1.d0,0.d0)
      b=b-b0 ! here b becomes residual
      call solve_info
      call solve_system
      xval=xval1+omg*xval
      ! judge convergence
      call step_info(istep,1,ndf,ntf,conve,er,xval,xval1,1)
      if(toln.lt.1)then ! judge by error
        if(maxval(er).lt.toln)exit
      else ! judge by max step
        if(istep.ge.int(toln))exit
      endif
      
    enddo
    timew='Time of non-linear iteration:'
    call time_record(timew,datev,timev,4,2)
    write(4,*)
    
    deallocate (ia,ja,a,b,pret,lret,tef,itef,tval)
    
    write(*,*)
    write(*,*)'          Outputing results'
    call output
    timew='Output results :'
    call time_record(timew,datev,timev,4,2)
    
  elseif(istrans.eq.1.and.isnonli.eq.1)then ! transient and non-linear problem
    
    write(*,*)
    write(*,*)'          Non-linear transient iteration'
    write(4,*)
    write(4,*)'          Non-linear transient iteration'
    allocate (a(ka),b(ntf),isfirtyp(ntf),c(kc),r(ka),tval(ntp),xval1(ntf),xval0(ntf),conve(ntf),er(ndf))
    allocate (b0(ntf))
    if(ndf.ne.1.and.istm.eq.1) then
      ! read nodal temperatures
      rewind (45)
      read(45,*) ktp
      do ik=1,min(ktp,ntp)
        read(45,*)ip,tval(ip)
      enddo
    endif
    open(99,file='RESULTS.QQQ',status='unknown',form='unformatted') ! store data section of tecplot file
    rewind(99)
    istep=0
    call output_trans(istep)
    do
      istep=istep+1
      xval0=xval ! record the result of last time step
      do jstep=1,50
        xval1=xval ! record the results of last step for comparision
        a=0.d0
        b=0.d0
        isfirtyp=0
        c=0.d0
        r=0.d0
        call coef_nonli(istep,tval)
        call first_type_bc(istep)
        call second_type_bc(istep)
        call third_type_bc(istep)
        call source_item(istep)
        if(ndf.eq.1)then ! transient heat
          call trans_heat_coef
        else ! dynamics
          write(*,*)'Non-linear dynamics is to be developed!'
          write(4,*)'Non-linear dynamics is to be developed!'
          ! to be developed
          stop
        endif
        ! forming matrix b
        call matmulmkl_mmk(ntf,ntf,ia,ja,r,ka,xval0,b,-1.d0,1.d0)
        ! get b0=a(x)*x
        call matmulmkl_mmk(ntf,ntf,ia,ja,a,ka,xval,b0,1.d0,0.d0)
        b=b-b0 ! here b becomes residual
        call solve_info
        call solve_system
        xval=xval1+omg*xval
        ! judge convergence
        call step_info(jstep,1,ndf,ntf,conve,er,xval,xval1,1)
        if(toln.lt.1)then ! judge by error
          if(maxval(er).lt.toln)exit
        else ! judge by max step
          if(jstep.ge.int(toln))exit
        endif
        
      enddo ! end of non-linear iteration
      
      !write(*,*)'  total non-linear iteration step=',jstep
      !write(4,*)'  total non-linear iteration step=',jstep
      call output_trans(istep)
      call step_info(istep,jstep,ndf,ntf,conve,er,xval,xval0,2)
      if(tolt.lt.1)then ! judge by error
        if(maxval(er).lt.tolt)exit
      else ! judge by max step
        if(istep.ge.int(tolt))exit
      endif
    enddo ! end of transien iteration
    
    write(*,*)'Total transient iteration steps=',istep
    write(4,*)'Total transient iteration steps=',istep
    timew='Time of nonlinear transient iteration:'
    call time_record(timew,datev,timev,4,2)
    
    ! write file 45 temperature.dat
    if(ndf.eq.1) then
      write(45,*) ntp
      do ip=1,ntp
        write(45,*)ip,xval(ip)
      enddo
    endif
    
    ! copy temporary file 99 to file 33 to complete data module of tecplot binary file
    call tecp_gen_trans2(33,99,ndim,ndf,nsig,mnode,ntp,nte,istep)
    
  endif
   
  write(*,*)
  write(*,*)'   Job terminated normally'
  timew='   Total computational time is'
  write(*,*)
  write(4,*)
  call time_record(timew,date0,time0,4,2)
  
end program EDM
      
subroutine step_info(istep,jstep,ndf,ntf,conve,er,xval,xval1,iflag)
  implicit none
  integer itf,idf,ip
  integer ndf,ntf,istep,jstep,iflag
  real*8 conve(ntf),er(ndf),xval(ntf),xval1(ntf)
  
  do itf=1,ntf
    if(abs(xval1(itf)).gt.1d-16.or.abs(xval(itf)).gt.1d-16)then
      conve(itf)=(xval(itf)-xval1(itf))/xval1(itf)
    else
      conve(itf)=0.d0
    endif
    !if(dabs(conve(itf)).gt.1.d-2)then
    !  write(4,*)itf/ndf,itf-itf/ndf*ndf,conve(itf)
    !endif
  enddo
  conve=abs(conve)
  do idf=1,ndf
    er(idf)=0.d0
    do ip=1,ntf/ndf
      itf=(ip-1)*ndf+idf
      if(conve(itf).gt.er(idf))er(idf)=conve(itf)
      !er(idf)=er(idf)+conve(itf)**2
    enddo
    !er(idf)=dsqrt(er(idf)/dble(ntf/ndf))
  enddo
  if(iflag.eq.2)then
    write(4,'(2i6,1p5e17.6)')istep,jstep,er(1:ndf)
    write(*,'(2i6,1p5e17.6)')istep,jstep,er(1:ndf)
  elseif(iflag.eq.1)then
    write(4,'(i5,1p5e15.4)')istep,er(1:ndf)
    write(*,'(i5,1p5e15.4)')istep,er(1:ndf)
  endif
  
end subroutine step_info

subroutine time_record(timew,datev,timev,ifile,iflag)
  ! ifile: file channal number of output time
  ! iflag=0,only output to ifile
  ! iflag=1,only output to screen
  ! iflag=2,output to both ifile and screen
  implicit none
  character*10 :: datev,timev,date0,time0
  character*60 :: timew
  integer :: day,hour,minute,other,iflag,ifile
  real*8 :: second,elapse
  integer,external :: str2num
  date0(:)=datev
  time0(:)=timev
  call date_and_time (datev,timev)
  day = str2num(datev)-str2num(date0) ! warning: may not be right when time cross two months
  hour = str2num(timev(1:2))-str2num(time0(1:2))
  minute = str2num(timev(3:4))-str2num(time0(3:4))
  second = dble(str2num(timev(5:6))-str2num(time0(5:6)))
  other = str2num(timev(8:10))-str2num(time0(8:10))
  elapse = day*24*3600.+hour*3600.+minute*60.+second+other/1000.
  hour=floor(elapse/3600.)
  minute=floor((elapse-hour*3600.)/60.)
  second=elapse-hour*3600.-minute*60.
  
  if(iflag.gt.0) write(*,360) timew,hour,minute,second
  if(iflag.ne.1) write(ifile,360) timew,hour,minute,second
360   format(1x,a,i6,' hours ',i3,' minites ',f7.3,' seconds')
end subroutine time_record

function str2num(s) result(f)
  implicit none
  character(len = *) :: s
  integer :: f
  read(s,*)f
end function str2num