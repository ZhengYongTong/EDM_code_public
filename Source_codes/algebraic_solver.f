      !****************************************************************************
      ! PARDISO_UNSYM: LU sparse solver from intel mkl library (recommended)
      ! DGESV: LU dense solver from intel mkl library
      ! CMLIB_LUD: linear equation solver and inverse using LU decomposition (dense solver)
      ! DSS_UNSYM: another LU sparse solver from intel mkl library
      ! FGMRES_ILU0: ILU0 precondition GMRES iteration sparse solver from intel mkl library
      ! FGMRES_ILUT: ILUT precondition GMRES iteration sparse solver from intel mkl library
      ! MATMULMKL: Dense matrix multiply dense matrix using DGEMM in MKL library
      ! MATMULMKL_MMK: Sparse matrix multiply sparse matrix using DGEMM in MKL library
      !****************************************************************************
      
! ##############################################################################
!     The begin of solver PARDISO_UNSYM
! ##############################################################################
      SUBROUTINE PARDISO_UNSYM(KA,IA,JA,A,N,B,X)
      IMPLICIT NONE
      
      INTEGER N,I,KA
      INTEGER MAXFCT, MNUM, MTYPE, PHASE, NRHS, ERROR, MSGLVL
      INTEGER IA(N+1),JA(KA)
      INTEGER IPARM(64)
      INTEGER IDUM(1)
      INTEGER*8 PT(64)
      REAL*8 X(N),A(KA),B(N)
      REAL*8  DDUM(1)

      DATA NRHS /1/, MAXFCT /1/, MNUM /1/
      
      ! SET UP PARDISO CONTROL PARAMETER
      MTYPE = 11 ! REAL UNSYMMETRIC
      ERROR = 0 ! INITIALIZE ERROR FLAG
      MSGLVL = 0 ! PRINT STATISTICAL INFORMATION
      CALL PARDISOINIT (PT, MTYPE, IPARM) ! INITILIAZE PARDISO SOLVER
      IPARM(1) = 1 ! NO SOLVER DEFAULT
      IPARM(2) = 3 ! FILL-IN REORDERING FROM METIS
      IPARM(8) = 9 ! NUMBERS OF ITERATIVE REFINEMENT STEPS
      
      ! REORDERING AND SYMBOLIC FACTORIZATION, THIS STEP ALSO ALLOCATES
      ! ALL MEMORY THAT IS NECESSARY FOR THE FACTORIZATION
      PHASE = 11 ! ONLY REORDERING AND SYMBOLIC FACTORIZATION
      CALL PARDISO (PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA,
     & IDUM, NRHS, IPARM, MSGLVL, DDUM, DDUM, ERROR)
      IF(ERROR.NE.0)THEN
        CALL PAR_ER(ERROR,IPARM(15),IPARM(16),IPARM(17),1)
        STOP
      ENDIF
      
      ! FACTORIZATION.
      PHASE = 22 ! ONLY FACTORIZATION
      CALL PARDISO (PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA,
     & IDUM, NRHS, IPARM, MSGLVL, DDUM, DDUM, ERROR)
      IF(ERROR.NE.0)THEN
        CALL PAR_ER(ERROR,IPARM(15),IPARM(16),IPARM(17),2)
        STOP
      ENDIF
      
      ! BACK SUBSTITUTION AND ITERATIVE REFINEMENT
      IPARM(8) = 2 ! MAX NUMBERS OF ITERATIVE REFINEMENT STEPS
      PHASE = 33 ! ONLY SOLVE
      CALL PARDISO (PT, MAXFCT, MNUM, MTYPE, PHASE, N, A, IA, JA,
     & IDUM, NRHS, IPARM, MSGLVL, B, X, ERROR)
      IF(ERROR.NE.0)THEN
        CALL PAR_ER(ERROR,IPARM(15),IPARM(16),IPARM(17),3)
        STOP
      ENDIF
      
      ! TERMINATION AND RELEASE OF MEMORY
      PHASE = -1 ! RELEASE INTERNAL MEMORY
      CALL PARDISO (PT, MAXFCT, MNUM, MTYPE, PHASE, N, DDUM, IDUM, IDUM,
     & IDUM, NRHS, IPARM, MSGLVL, DDUM, DDUM, ERROR)
      IF(ERROR.NE.0)THEN
        CALL PAR_ER(ERROR,IPARM(15),IPARM(16),IPARM(17),4)
        STOP
      ENDIF

      END
      
      SUBROUTINE PAR_ER(ERROR,IPARM15,IPARM16,IPARM17,NSTEP)
      IMPLICIT NONE
      
      INTEGER ERROR ! ERROR TYPE
      INTEGER NSTEP ! ERROR STEP
      INTEGER IPARM15,IPARM16,IPARM17 ! STORAGE USE
      
      WRITE(*,1)ERROR,NSTEP
      WRITE(4,1)ERROR,NSTEP
      IF(ERROR.EQ.-2)THEN
        WRITE(*,2)MAX(IPARM15,IPARM16+IPARM17),NSTEP
        WRITE(4,2)MAX(IPARM15,IPARM16+IPARM17),NSTEP
      ENDIF
      RETURN
      
1     FORMAT('ERROR NUMBER',I4,'OCCURRED IN STEP',I3,'OF PARDISO')
2     FORMAT(I10,'KB MEMORY IS NEEDED IN STEP',I3,'OF PARDISO')
      END SUBROUTINE PAR_ER
! ##############################################################################
!     The end of solver PARDISO_UNSYM
! ##############################################################################
      
! ##############################################################################
!     The begin of solver DSS_UNSYM
! ##############################################################################
      SUBROUTINE DSS_UNSYM(KA,IA,JA,A,NTF,B,X)
      IMPLICIT NONE
      INCLUDE 'MKL_DSS.FI'
      
      INTEGER NTF,KA,I,NRHS
      PARAMETER (NRHS=1)
      INTEGER IA(NTF+1),JA(KA),IDUM(1)
      DOUBLE PRECISION A(KA),B(NTF),X(NTF)
      INTEGER*8 HANDLE
      INTEGER ERROR
      CHARACTER*100 STATIN
      DOUBLE PRECISION STATOUT(6)
      INTEGER BUFLEN
      PARAMETER(BUFLEN=100)
      INTEGER BUFF(BUFLEN)
      ! INITIALIZE THE SOLVER
      ERROR=DSS_CREATE(HANDLE,MKL_DSS_DEFAULTS)
      IF(ERROR.NE.MKL_DSS_SUCCESS)THEN
        WRITE(*,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_CREATE.'
        WRITE(4,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_CREATE.'
        STOP
      ENDIF
      
      ! DEFINE THE NON-ZERO STRUCTURE OF THE MATRIX
      ERROR=DSS_DEFINE_STRUCTURE(HANDLE,MKL_DSS_NON_SYMMETRIC,
     &IA,NTF,NTF,JA,KA)
      IF(ERROR.NE.MKL_DSS_SUCCESS)THEN
        WRITE(*,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_DEFINE_STRUCTURE.'
        WRITE(4,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_DEFINE_STRUCTURE.'
        STOP
      ENDIF
      
      ! REORDER THE MATRIX
      ERROR=DSS_REORDER(HANDLE,MKL_DSS_DEFAULTS,IDUM)
      !ERROR=DSS_REORDER(HANDLE,MKL_DSS_METIS_OPENMP_ORDER,IDUM)
      !ERROR=DSS_REORDER(HANDLE,MKL_DSS_GET_ORDER,IDUM)
      IF(ERROR.NE.MKL_DSS_SUCCESS)THEN
        WRITE(*,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_REORDER.'
        WRITE(4,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_REORDER.'
        STOP
      ENDIF

      ! FACTOR THE MATRIX
      ERROR=DSS_FACTOR_REAL(HANDLE,MKL_DSS_DEFAULTS,A)
      IF(ERROR.NE.MKL_DSS_SUCCESS)THEN
        WRITE(*,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_FACTOR_REAL.'
        WRITE(4,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_FACTOR_REAL.'
        STOP
      ENDIF
      
      ! GET THE SOLUTION VECTOR
      ERROR=DSS_SOLVE_REAL(HANDLE,MKL_DSS_DEFAULTS,B,NRHS,X)
      IF(ERROR.NE.MKL_DSS_SUCCESS)THEN
        WRITE(*,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_SOLVE_REAL.'
        WRITE(4,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_SOLVE_REAL.'
        STOP
      ENDIF
      
      ! PRINT DETERMINANT OF THE MATRIX (NO STATISTICS FOR A DIAGONAL MATRIX)
      STATIN='REORDERTIME,FACTORTIME,SOLVETIME,PEAKMEM,FACTORMEM,
     &SOLVEMEM'
      CALL MKL_CVT_TO_NULL_TERMINATED_STR(BUFF,BUFLEN,STATIN)
      ERROR=DSS_STATISTICS(HANDLE,MKL_DSS_DEFAULTS,BUFF,STATOUT)
      WRITE(4,"(' THE REORDER TIME IS',F10.3,' SECOND')") STATOUT(1)
      WRITE(4,"(' THE FACTOR TIME IS',F10.3,' SECOND')") STATOUT(2)
      WRITE(4,"(' THE SOLVE TIME IS', F10.3,' SECOND')") STATOUT(3)
      WRITE(4,"(' THE PEAK MEMORY IS', F15.3,' KB')") STATOUT(4)
      WRITE(4,"(' THE FACTOR MEMORY IS',F15.3,' KB')") STATOUT(5)
      WRITE(4,"(' THE SOLVE MEMORY IS',F15.3,' KB')") STATOUT(6)
      
      ! DEALLOCATE SOLVER STORAGE
      ERROR = DSS_DELETE( HANDLE, MKL_DSS_DEFAULTS )
      IF(ERROR.NE.MKL_DSS_SUCCESS)THEN
        WRITE(*,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_DELETE.'
        WRITE(4,*)'ERROR CODE',ERROR,'OCCURRED IN DSS_DELETE.'
        STOP
      ENDIF
      
      END
! ##############################################################################
!     The end of solver DSS_UNSYM
! ##############################################################################
      
! ##############################################################################
!     The begin of solver FGMRES_ILU0
! ##############################################################################
      SUBROUTINE FGMRES_ILU0(KA,IA,JA,A,NTF,B,X)
      IMPLICIT NONE

      INTEGER NTF,KA,I
      INTEGER IA(NTF),JA(NTF)
      DOUBLE PRECISION A(KA),B(NTF),X(NTF)
      INTEGER IPAR(128),IERR
      DOUBLE PRECISION DPAR(128)
      DOUBLE PRECISION,ALLOCATABLE::RHS(:),RESIDUAL(:),TMP(:),TRVEC(:),
     &BILU0(:)
      INTEGER ITERCOUNT ! NUMBER OF ITERATIONS
      INTEGER RCI_REQUEST
      DOUBLE PRECISION DVAR ! THE EUCLIDEAN NORM OF VECTOR RESIDUAL
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2 ! COMPUTES THE EUCLIDEAN NORM OF A VECTOR
      
      ALLOCATE(RHS(NTF),RESIDUAL(NTF),TRVEC(NTF),BILU0(KA))
      
	CALL DCOPY(NTF,B,1,RHS,1) ! SAVE B IN VECTOR RHS FOR FUTURE USE
      ! INITIALIZE THE INITIAL GUESS SOLUTION
      DO I=1,NTF
         X(I)=0.D0
      ENDDO
      
      ! INITIALIZE THE SOLVER
      IPAR(15)=MIN(NTF,24000000000/NTF) ! DO THE RESTART AFTER MIN(NTF,24000000000/NTF) ITERATIONS
      ALLOCATE(TMP(NTF*(2*IPAR(15)+1)+IPAR(15)*(IPAR(15)+9)/2+1))
      CALL DFGMRES_INIT(NTF,X,B,RCI_REQUEST,IPAR,DPAR,TMP)
      IF(RCI_REQUEST.NE.0)THEN
        WRITE(*,*)'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES_INIT.'
        WRITE(4,*)'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES_INIT.'
        CALL MKL_FREE_BUFFERS
        STOP
      ENDIF

      ! SOME SPECIAL IPAR AND DPAR FOR DCSRILU0 ENTRIES ARE SET BELOW:
      ! IPAR(31)= 1 - CHANGE SMALL DIAGONAL VALUE TO THAT GIVEN BY DPAR(32)
	IPAR(31)=1
	DPAR(31)=1.D-16 ! DEFAUT VALUE
	DPAR(32)=1.D-10 ! DEFAUT VALUE
      
      ! CALCULATE ILU0 PRECONDITIONER.
      CALL DCSRILU0(NTF,A,IA,JA,BILU0,IPAR,DPAR,IERR)
	IF(IERR.NE.0) THEN
	  WRITE(*,*)'ERROR CODE',IERR,'OCCURRED IN PRECONDITION.'
        WRITE(4,*)'ERROR CODE',IERR,'OCCURRED IN PRECONDITION.'
        CALL MKL_FREE_BUFFERS
        STOP
      ENDIF

      IPAR(5)=MIN(NTF,24000000000/NTF) ! MAXIMAL NUMBER OF ITERATIONS
      IPAR(8)=1 ! DO THE STOPPING TEST FOR THE MAXIMAL NUMBER OF ITERATIONS
      IPAR(11)=1 ! SET PARAMETER IPAR(11) FOR PRECONDITIONER CALL
      
      ! CHECK THE CORRECTNESS AND CONSISTENCY OF THE NEWLY SET PARAMETERS
      CALL DFGMRES_CHECK(NTF,X,B,RCI_REQUEST,IPAR,DPAR,TMP)
      IF(RCI_REQUEST.NE.0)THEN
        WRITE(*,*)'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES_CHECK.'
        WRITE(4,*)'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES_CHECK.'
        CALL MKL_FREE_BUFFERS
        STOP
      ENDIF
      
      ! COMPUTE THE SOLUTION BY RCI (P)FGMRES SOLVER WITH PRECONDITIONING
      ! REVERSE COMMUNICATION STARTS HERE
1     CALL DFGMRES(NTF,X,B,RCI_REQUEST,IPAR,DPAR,TMP)
      
      IF (RCI_REQUEST.EQ.0) THEN
        ! IF RCI_REQUEST=0, THEN THE SOLUTION WAS FOUND WITH THE REQUIRED PRECISION
        GOTO 3
      ELSEIF(RCI_REQUEST.EQ.1)THEN
        ! IF RCI_REQUEST=1, THEN COMPUTE THE VECTOR A*TMP(IPAR(22))
        ! AND PUT THE RESULT IN VECTOR TMP(IPAR(23))
      	CALL MKL_DCSRGEMV('N',NTF,A,IA,JA,TMP(IPAR(22)),TMP(IPAR(23)))
      	GOTO 1
      ELSEIF(RCI_REQUEST.EQ.2) THEN
        ! IF RCI_REQUEST=2, THEN DO THE USER-DEFINED STOPPING TEST
        ! THE RESIDUAL STOPPING TEST FOR THE COMPUTED SOLUTION IS PERFORMED HERE
      	IPAR(13)=1 ! REQUEST TO THE DFGMRES_GET ROUTINE TO PUT THE SOLUTION INTO RHS(NTF) VIA IPAR(13)
        ! GET THE CURRENT FGMRES SOLUTION IN THE VECTOR RHS(NTF)
      	CALL DFGMRES_GET(NTF,X,RHS,RCI_REQUEST,IPAR,DPAR,TMP,ITERCOUNT)
        ! COMPUTE THE CURRENT TRUE RESIDUAL
      	CALL MKL_DCSRGEMV('N',NTF,A,IA,JA,RHS,RESIDUAL) ! RESIDUAL=A*RHS
      	CALL DAXPY(NTF,-1.0D0,B,1,RESIDUAL,1) ! RESIDUAL=RESIDUAL-B
      	DVAR=DNRM2(NTF,RESIDUAL,1)
      	IF(DVAR.LT.1.0E-6) THEN
      	   GOTO 3
      	ELSE
      	   GOTO 1
      	ENDIF
      ELSEIF(RCI_REQUEST.EQ.3)THEN
        ! IF RCI_REQUEST=3, THEN APPLY THE PRECONDITIONER ON THE VECTOR
        ! TMP(IPAR(22)) AND PUT THE RESULT IN VECTOR TMP(IPAR(23))
        CALL MKL_DCSRTRSV('L','N','U',NTF,BILU0,IA,JA,TMP(IPAR(22)),
     &  TRVEC)
        CALL MKL_DCSRTRSV('U','N','N',NTF,BILU0,IA,JA,TRVEC,
     & TMP(IPAR(23)))
        GOTO 1
      ELSEIF(RCI_REQUEST.EQ.4)THEN
        ! IF RCI_REQUEST=4, THEN CHECK IF THE NORM OF THE NEXT GENERATED VECTOR IS
        ! NOT ZERO UP TO ROUNDING AND COMPUTATIONAL ERRORS. THE NORM IS CONTAINED
        ! IN DPAR(7) PARAMETER.
      	IF(DPAR(7).LT.1.0D-12)THEN
      	   GOTO 3
      	ELSE
      	   GOTO 1
      	ENDIF
      ELSE ! IF RCI_REQUEST=ANYTHING ELSE, THEN DFGMRES SUBROUTINE FAILED
      	WRITE(*,*) 'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES.'
        WRITE(4,*) 'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES.'
        CALL MKL_FREE_BUFFERS
        STOP
      ENDIF
      ! GET THE CURRENT ITERATION NUMBER AND THE FGMRES SOLUTION
      ! REQUEST TO DFGMRES_GET TO PUT THE SOLUTION INTO VECTOR X(NTF) VIA IPAR(13)
3     IPAR(13)=0
      CALL DFGMRES_GET(NTF,X,B, RCI_REQUEST,IPAR,DPAR,TMP,ITERCOUNT)
      WRITE(*,*) 'TOTAL ITERATIONS IS ',ITERCOUNT
      WRITE(4,*) 'TOTAL ITERATIONS IS ',ITERCOUNT
      IF(ITERCOUNT.EQ.IPAR(5))THEN
        WRITE(*,*)'***WARNING: THE RESULTS MAY NOT BE RIGHT!***'
        WRITE(4,*)'***WARNING: THE RESULTS MAY NOT BE RIGHT!***'
      ENDIF
      
      ! RELEASE INTERNAL MKL MEMORY
      CALL MKL_FREE_BUFFERS

      RETURN
      END
! ##############################################################################
!     The end of solver FGMRES_ILU0
! ##############################################################################
      
! ##############################################################################
!     The begin of solver FGMRES_ILUT
! ##############################################################################
      SUBROUTINE FGMRES_ILUT(KA,IA,JA,A,NTF,B,X,MAXHB)
      IMPLICIT NONE

      INTEGER NTF,KA,I,MAXHB
      INTEGER IA(NTF),JA(NTF),IBILUT(NTF)
      DOUBLE PRECISION A(KA),B(NTF),X(NTF)
      INTEGER IPAR(128),IERR
      INTEGER*8 KBILUT
      DOUBLE PRECISION DPAR(128)
      INTEGER,ALLOCATABLE::JBILUT(:)
      DOUBLE PRECISION,ALLOCATABLE::RHS(:),RESIDUAL(:),TMP(:),TRVEC(:),
     &BILUT(:)
      INTEGER ITERCOUNT ! NUMBER OF ITERATIONS
      INTEGER RCI_REQUEST
      DOUBLE PRECISION DVAR ! THE EUCLIDEAN NORM OF VECTOR RESIDUAL
      DOUBLE PRECISION DNRM2
      EXTERNAL DNRM2 ! COMPUTES THE EUCLIDEAN NORM OF A VECTOR
      INTEGER MAXFIL
      DOUBLE PRECISION TOL
      
	MAXFIL=MAXHB
      KBILUT=(2*MAXFIL+1)*NTF-MAXFIL*(MAXFIL+1)+1
      ALLOCATE(RHS(NTF),RESIDUAL(NTF),TRVEC(NTF),BILUT(KBILUT),
     &JBILUT(KBILUT))
      
	CALL DCOPY(NTF,B,1,RHS,1) ! SAVE B IN VECTOR RHS FOR FUTURE USE
      ! INITIALIZE THE INITIAL GUESS SOLUTION
      DO I=1,NTF
         X(I)=0.D0
      ENDDO
      
      ! INITIALIZE THE SOLVER
      IPAR(15)=MIN(NTF,24000000000/NTF) ! DO THE RESTART AFTER MIN(NTF,24000000000/NTF) ITERATIONS
      ALLOCATE(TMP(NTF*(2*IPAR(15)+1)+IPAR(15)*(IPAR(15)+9)/2+1))
      CALL DFGMRES_INIT(NTF,X,B,RCI_REQUEST,IPAR,DPAR,TMP)
      IF(RCI_REQUEST.NE.0)THEN
        WRITE(*,*)'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES_INIT.'
        WRITE(4,*)'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES_INIT.'
        CALL MKL_FREE_BUFFERS
        STOP
      ENDIF

      ! SOME SPECIAL IPAR AND DPAR FOR DCSRILU0 ENTRIES ARE SET BELOW:
      ! IPAR(31)= 1 - CHANGE SMALL DIAGONAL VALUE TO THAT GIVEN BY DPAR(31),
      ! DPAR(31)= 1.D-5  INSTEAD OF THE DEFAULT VALUE SET BY DFGMRES_INIT.
      !                  IT IS THE TARGET VALUE OF THE DIAGONAL VALUE IF IT IS
      !                  SMALL AS COMPARED TO GIVEN TOLERANCE MULTIPLIED
      !                  BY THE MATRIX ROW NORM AND THE ROUTINE SHOULD
      !                  CHANGE IT RATHER THAN ABORT DCSRILUT CALCULATIONS.
	IPAR(31)=1
	DPAR(31)=1.D-5
	TOL=1.D-6
      
      ! CALCULATE ILUT PRECONDITIONER.
      CALL DCSRILUT(NTF,A,IA,JA,BILUT,IBILUT,JBILUT,TOL,MAXFIL,IPAR,DPAR
     &,IERR)
	IF(IERR.NE.0) THEN
	  WRITE(*,*)'ERROR CODE',IERR,'OCCURRED IN PRECONDITION.'
        WRITE(4,*)'ERROR CODE',IERR,'OCCURRED IN PRECONDITION.'
        CALL MKL_FREE_BUFFERS
        STOP
      ENDIF

      IPAR(5)=MIN(NTF,24000000000/NTF) ! MAXIMAL NUMBER OF ITERATIONS
      IPAR(8)=1 ! DO THE STOPPING TEST FOR THE MAXIMAL NUMBER OF ITERATIONS
      IPAR(11)=1 ! SET PARAMETER IPAR(11) FOR PRECONDITIONER CALL
      
      ! CHECK THE CORRECTNESS AND CONSISTENCY OF THE NEWLY SET PARAMETERS
      CALL DFGMRES_CHECK(NTF,X,B,RCI_REQUEST,IPAR,DPAR,TMP)
      IF(RCI_REQUEST.NE.0)THEN
        WRITE(*,*)'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES_CHECK.'
        WRITE(4,*)'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES_CHECK.'
        CALL MKL_FREE_BUFFERS
        STOP
      ENDIF
      
      ! COMPUTE THE SOLUTION BY RCI (P)FGMRES SOLVER WITH PRECONDITIONING
      ! REVERSE COMMUNICATION STARTS HERE
1     CALL DFGMRES(NTF,X,B,RCI_REQUEST,IPAR,DPAR,TMP)
      
      IF (RCI_REQUEST.EQ.0) THEN
        ! IF RCI_REQUEST=0, THEN THE SOLUTION WAS FOUND WITH THE REQUIRED PRECISION
        GOTO 3
      ELSEIF(RCI_REQUEST.EQ.1)THEN
        ! IF RCI_REQUEST=1, THEN COMPUTE THE VECTOR A*TMP(IPAR(22))
        ! AND PUT THE RESULT IN VECTOR TMP(IPAR(23))
      	CALL MKL_DCSRGEMV('N',NTF,A,IA,JA,TMP(IPAR(22)),TMP(IPAR(23)))
      	GOTO 1
      ELSEIF(RCI_REQUEST.EQ.2) THEN
        ! IF RCI_REQUEST=2, THEN DO THE USER-DEFINED STOPPING TEST
        ! THE RESIDUAL STOPPING TEST FOR THE COMPUTED SOLUTION IS PERFORMED HERE
      	IPAR(13)=1 ! REQUEST TO THE DFGMRES_GET ROUTINE TO PUT THE SOLUTION INTO RHS(NTF) VIA IPAR(13)
        ! GET THE CURRENT FGMRES SOLUTION IN THE VECTOR RHS(NTF)
      	CALL DFGMRES_GET(NTF,X,RHS,RCI_REQUEST,IPAR,DPAR,TMP,ITERCOUNT)
        ! COMPUTE THE CURRENT TRUE RESIDUAL
      	CALL MKL_DCSRGEMV('N',NTF,A,IA,JA,RHS,RESIDUAL) ! RESIDUAL=A*RHS
      	CALL DAXPY(NTF,-1.0D0,B,1,RESIDUAL,1) ! RESIDUAL=RESIDUAL-B
      	DVAR=DNRM2(NTF,RESIDUAL,1)
      	IF(DVAR.LT.1.0E-6) THEN
      	   GOTO 3
      	ELSE
      	   GOTO 1
      	ENDIF
      ELSEIF(RCI_REQUEST.EQ.3)THEN
        ! IF RCI_REQUEST=3, THEN APPLY THE PRECONDITIONER ON THE VECTOR
        ! TMP(IPAR(22)) AND PUT THE RESULT IN VECTOR TMP(IPAR(23))
        CALL MKL_DCSRTRSV('L','N','U',NTF,BILUT,IBILUT,JBILUT,
     &  TMP(IPAR(22)),TRVEC)
        CALL MKL_DCSRTRSV('U','N','N',NTF,BILUT,IBILUT,JBILUT,TRVEC,
     & TMP(IPAR(23)))
        GOTO 1
      ELSEIF(RCI_REQUEST.EQ.4)THEN
        ! IF RCI_REQUEST=4, THEN CHECK IF THE NORM OF THE NEXT GENERATED VECTOR IS
        ! NOT ZERO UP TO ROUNDING AND COMPUTATIONAL ERRORS. THE NORM IS CONTAINED
        ! IN DPAR(7) PARAMETER.
      	IF(DPAR(7).LT.1.0D-12)THEN
      	   GOTO 3
      	ELSE
      	   GOTO 1
      	ENDIF
      ELSE ! IF RCI_REQUEST=ANYTHING ELSE, THEN DFGMRES SUBROUTINE FAILED
      	WRITE(*,*) 'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES.'
        WRITE(4,*) 'ERROR CODE ',RCI_REQUEST,'OCCURRED IN DFGMRES.'
        CALL MKL_FREE_BUFFERS
        STOP
      ENDIF
      ! GET THE CURRENT ITERATION NUMBER AND THE FGMRES SOLUTION
      ! REQUEST TO DFGMRES_GET TO PUT THE SOLUTION INTO VECTOR X(NTF) VIA IPAR(13)
3     IPAR(13)=0
      CALL DFGMRES_GET(NTF,X,B, RCI_REQUEST,IPAR,DPAR,TMP,ITERCOUNT)
      WRITE(*,*) 'TOTAL ITERATIONS IS ',ITERCOUNT
      WRITE(4,*) 'TOTAL ITERATIONS IS ',ITERCOUNT
      IF(ITERCOUNT.EQ.IPAR(5))THEN
        WRITE(*,*)'***WARNING: THE RESULTS MAY NOT BE RIGHT!***'
        WRITE(4,*)'***WARNING: THE RESULTS MAY NOT BE RIGHT!***'
      ENDIF
      
      ! RELEASE INTERNAL MKL MEMORY
      CALL MKL_FREE_BUFFERS

      RETURN
      END
! ##############################################################################
!     The end of solver FGMRES_ILUT
! ##############################################################################
      
! ##############################################################################
!     The begin of solver CMLIB_LUD
! ##############################################################################
      SUBROUTINE CMLIB_LUD(NROW,NCOL,A,N,INDIC,INFO)
  ! -----------------------------------------------------------------
  !     INVERSE OF MATRIX, SOLVER OF EQUATIONS BY PARTIAL PIVOTING
  !     USING LU_DCOMPOSITION. INDIC=-1,ONLY INVERSE OF MATRIX;
  !     INIDC=0,INVERSE AND SOLVE EQUATION; INDIC=1,ONLY SOLVE EQUATION
  !    --------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NROW,NCOL),WORK(N),IPVT(N),DET(2)
	CALL DGEFA(A,NROW,N,IPVT,INFO)  ! LU-DECOMPOSITION
      IF(INFO.NE.0) WRITE(*,11)INFO
 11   FORMAT(/,'     INFO IN LU-DECOMPOSITION =',I8,/,                  &
     & '    (INFO = 0: NORMAL; OTHERWISE, IT IS NO. OF SINGULAR ROW)')
      IF(INDIC.GE.0) THEN
       DO IC=N+1,NCOL
        CALL DGESL(A,NROW,N,IPVT,A(1,IC),0)  ! SOLVING EQUATIONS
	 ENDDO
      ENDIF
      IF(INDIC.LE.0) CALL DGEDI(A,NROW,N,IPVT,DET,WORK,01) ! INVERT MATRIX
      END
  !    --------------------------------------------------------------
  !    LU-DECOMPOSITION PART
  !    --------------------------------------------------------------
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
C***BEGIN PROLOGUE  DGEFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2A1
C***KEYWORDS  DOUBLE PRECISION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Factors a double precision matrix by Gaussian elimination.
C***DESCRIPTION
C
C     DGEFA factors a double precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSCAL,IDAMAX
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSCAL,IDAMAX
C***END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER IDAMAXq,J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAXq(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCALq(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DAXPY
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A7
C***KEYWORDS  BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,TRIAD,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P computation y = a*x + y
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY), where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N
C       and LY is defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DAXPY
C
      DOUBLE PRECISION DX(*),DY(*),DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF(N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE
      RETURN
      END
      SUBROUTINE DSCALq(N,DA,DX,INCX)
C***BEGIN PROLOGUE  DSCAL
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A6
C***KEYWORDS  BLAS,LINEAR ALGEBRA,SCALE,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. vector scale x = a*x
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(1+I*INCX) with  DA * DX(1+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSCAL
C
      DOUBLE PRECISION DA,DX(*)
C***FIRST EXECUTABLE STATEMENT  DSCAL
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          DX(I) = DA*DX(I)
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I + 1) = DA*DX(I + 1)
        DX(I + 2) = DA*DX(I + 2)
        DX(I + 3) = DA*DX(I + 3)
        DX(I + 4) = DA*DX(I + 4)
   50 CONTINUE
      RETURN
      END
      INTEGER FUNCTION IDAMAXq(N,DX,INCX)
C***BEGIN PROLOGUE  IDAMAX
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A2
C***KEYWORDS  BLAS,DOUBLE PRECISION,LINEAR ALGEBRA,MAXIMUM COMPONENT,
C             VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Find largest component of d.p. vector
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX =  first I, I = 1 to N, to minimize  ABS(DX(1-INCX+I*INCX)
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  IDAMAX
C
      DOUBLE PRECISION DX(*),DMAX,XMAG
C***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAXq = 0
      IF(N.LE.0) RETURN
      IDAMAXq = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      DMAX = DABS(DX(1))
      NS = N*INCX
      II = 1
          DO 10 I = 1,NS,INCX
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 5
          IDAMAXq = II
          DMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
          XMAG = DABS(DX(I))
          IF(XMAG.LE.DMAX) GO TO 30
          IDAMAXq = I
          DMAX = XMAG
   30 CONTINUE
      RETURN
      END

      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
C***BEGIN PROLOGUE  DGESL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D2A1
C***KEYWORDS  DOUBLE PRECISION,LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Solves the double precision system  A*X=B or  TRANS(A)*X=B
C            using the factors computed by DGECO or DGEFA.
C***DESCRIPTION
C
C     DGESL solves the double precision system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B  where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGECO has set RCOND .GT. 0.0
C        or DGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DDOT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DDOT
C***END PROLOGUE  DGESL
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DDOT
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A4
C***KEYWORDS  BLAS,DOUBLE PRECISION,INNER PRODUCT,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  D.P. inner product of d.p. vectors
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY)
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DDOT
C
      DOUBLE PRECISION DX(*),DY(*)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     1   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
   70     CONTINUE
      RETURN
      END

      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
C***BEGIN PROLOGUE  DGEDI
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C***CATEGORY NO.  D3A1,D2A1
C***KEYWORDS  DETERMINANT,DOUBLE PRECISION,FACTOR,INVERSE,
C             LINEAR ALGEBRA,LINPACK,MATRIX
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  Computes the determinant and inverse of a matrix using
C            factors computed by DGECO or DGEFA.
C***DESCRIPTION
C
C     DGEDI computes the determinant and inverse of a matrix
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        WORK    DOUBLE PRECISION(N)
C                work vector.  Contents destroyed.
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        A       inverse of original matrix if requested.
C                Otherwise unchanged.
C
C        DET     DOUBLE PRECISION(2)
C                determinant of original matrix if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. DABS(DET(1)) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if DGECO has set RCOND .GT. 0.0 or DGEFA has set
C        INFO .EQ. 0 .
C
C     LINPACK.  This version dated 08/14/78 .
C     Cleve Moler, University of New Mexico, Argonne National Lab.
C
C     Subroutines and Functions
C
C     BLAS DAXPY,DSCAL,DSWAP
C     Fortran DABS,MOD
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  DAXPY,DSCAL,DSWAP
C***END PROLOGUE  DGEDI
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),DET(2),WORK(*)
C
      DOUBLE PRECISION T
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
C     COMPUTE DETERMINANT
C
C***FIRST EXECUTABLE STATEMENT  DGEDI
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCALq(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAPq(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
      SUBROUTINE DSWAPq(N,DX,INCX,DY,INCY)
C***BEGIN PROLOGUE  DSWAP
C***DATE WRITTEN   791001   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***REVISION HISTORY  (YYMMDD)
C   000330  Modified array declarations.  (JEC)
C
C***CATEGORY NO.  D1A5
C***KEYWORDS  BLAS,DOUBLE PRECISION,INTERCHANGE,LINEAR ALGEBRA,VECTOR
C***AUTHOR  LAWSON, C. L., (JPL)
C           HANSON, R. J., (SNLA)
C           KINCAID, D. R., (U. OF TEXAS)
C           KROGH, F. T., (JPL)
C***PURPOSE  Interchange d.p. vectors
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DX  input vector DY (unchanged if N .LE. 0)
C       DY  input vector DX (unchanged if N .LE. 0)
C
C     Interchange double precision DX and double precision DY.
C     For I = 0 to N-1, interchange  DX(LX+I*INCX) and DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = (-INCX)*N, and LY is
C     defined in a similar way using INCY.
C***REFERENCES  LAWSON C.L., HANSON R.J., KINCAID D.R., KROGH F.T.,
C                 *BASIC LINEAR ALGEBRA SUBPROGRAMS FOR FORTRAN USAGE*,
C                 ALGORITHM NO. 539, TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOLUME 5, NUMBER 3, SEPTEMBER 1979, 308-323
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSWAP
C
      DOUBLE PRECISION DX(*),DY(*),DTEMP1,DTEMP2,DTEMP3
C***FIRST EXECUTABLE STATEMENT  DSWAP
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C       CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP1 = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 3.
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP1 = DX(I)
        DTEMP2 = DX(I+1)
        DTEMP3 = DX(I+2)
        DX(I) = DY(I)
        DX(I+1) = DY(I+1)
        DX(I+2) = DY(I+2)
        DY(I) = DTEMP1
        DY(I+1) = DTEMP2
        DY(I+2) = DTEMP3
   50 CONTINUE
      RETURN
   60 CONTINUE
C
C     CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
      NS = N*INCX
        DO 70 I=1,NS,INCX
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   70   CONTINUE
      RETURN
      END
! ##############################################################################
!     The end of solver CMLIB_LUD
! ##############################################################################
      
! ##############################################################################
!     dense matrix multiply dense matrix using DGEMM in MKL library
! ##############################################################################
      SUBROUTINE MATMULMKL(M,K,N,A,B,C,ALPHA,BETA)
      ! TO COMPUTE (ALPHA*A*B+BETA*C) AND OUTPUT THE RESULTS INTO C
      IMPLICIT NONE
      INTEGER::M, N, K
      INTEGER::LDA, LDB, LDC
      REAL*8::ALPHA, BETA
      REAL*8::A(M,K), B(K,N), C(M,N)
      LDA=M
      LDB=K
      LDC=M
      CALL DGEMM('N','N',M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
      END SUBROUTINE MATMULMKL
      
! ##############################################################################
!     sparse matrix multiply sparse matrix using DGEMM in MKL library
! ##############################################################################
      SUBROUTINE MATMULMKL_MMK(M,K,IA,JA,A,KA,X,Y,ALPHA,BETA)
      ! TO COMPUTE (ALPHA*A*X+BETA*Y) AND OUTPUT THE RESULTS INTO Y
      IMPLICIT NONE
      INTEGER::M,K,KA
      INTEGER::JA(KA),IA(M+1),PNTRB(M),PNTRE(M)
      REAL*8::ALPHA,BETA
      REAL*8::A(KA),X(K),Y(M)
      CHARACTER::MATDESCRA(6)
      PNTRB(1:M)=IA(1:M)
      PNTRE(1:M)=IA(2:M+1)
      MATDESCRA(1)='G'
      MATDESCRA(4)='F'
      CALL MKL_DCSRMV('N',M,K,ALPHA,MATDESCRA,A,JA,PNTRB,PNTRE,X,BETA,Y)
      END SUBROUTINE MATMULMKL_MMK