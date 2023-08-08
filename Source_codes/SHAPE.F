      !****************************************************************************
      ! All the subroutines in this file are about element mapping and their shape functions
      !****************************************************************************
      
      SUBROUTINE EVAL_DNDX_DNDXX(IELEM,NDIM,NODE,XI,ELCOD,IETYPE,DNDX,
     &DNDXX,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(3),ELCOD(3,NODE),DNDX(3,NODE),DNDXI(3,NODE),
     &          DXIDX(3,3),DXDXI(3,3),DNDXX(3,3,NODE)
      CALL DN_DXI(IETYPE,NODE,XI,DNDXI)
      CALL VAL_JACOB(NDIM,NODE,ELCOD,IELEM,DNDXI,DXDXI,DJACB,DXIDX,DNDX)
      IF(IFLAG.EQ.2) CALL VAL_DNDXX(IETYPE,NDIM,NODE,XI,ELCOD,
     &                              DNDXI,DXDXI,DJACB,DXIDX,DNDXX)
      RETURN
      END

      SUBROUTINE VAL_JACOB(NDIM,NODE,ELCOD,IELEM,DNDXI,
     &                     DXDXI,DJACB,DXIDX,DNDX)
      IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION DNDXI(3,NODE),ELCOD(3,NODE),DNDX(3,NODE),DXIDX(3,3),
     &          DXDXI(3,3)
	DO IXI=1,NDIM
	 DO IDIM=1,NDIM
	  DXDXI(IXI,IDIM)=0.0D0
	  DO INODE=1,NODE
	   DXDXI(IXI,IDIM)=DXDXI(IXI,IDIM)+DNDXI(IXI,INODE)*
     &                   ELCOD(IDIM,INODE)
        ENDDO
       ENDDO
      ENDDO
      IF(NDIM.EQ.1)THEN
        DJACB=DXDXI(1,1)
      ELSEIF(NDIM.EQ.2)THEN
        DJACB=DXDXI(1,1)*DXDXI(2,2)-DXDXI(1,2)*DXDXI(2,1)
      ELSE
        DJACB=DXDXI(1,1)*(DXDXI(2,2)*DXDXI(3,3)-DXDXI(3,2)*
     & DXDXI(2,3))-DXDXI(1,2)*(DXDXI(2,1)*DXDXI(3,3)-DXDXI(3,1)*
     & DXDXI(2,3))+DXDXI(1,3)*(DXDXI(2,1)*DXDXI(3,2)-DXDXI(3,1)*
     & DXDXI(2,2))
      ENDIF
      IF(DABS(DJACB).LT.1.0D-12) THEN
        WRITE(*,600)REAL(DJACB),IELEM
        WRITE(4,600)REAL(DJACB),IELEM
        WRITE(7,600)REAL(DJACB),IELEM
        STOP
      ENDIF
      IF(NDIM.EQ.1) THEN
        DXIDX(1,1)=1.D0/DJACB
      ELSEIF(NDIM.EQ.2) THEN
	  DXIDX(1,1)=DXDXI(2,2)/DJACB
	  DXIDX(2,2)=DXDXI(1,1)/DJACB
	  DXIDX(1,2)=-DXDXI(1,2)/DJACB
	  DXIDX(2,1)=-DXDXI(2,1)/DJACB
      ELSE
        DXIDX(1,1)=(DXDXI(2,2)*DXDXI(3,3)-DXDXI(3,2)*DXDXI(2,3))/DJACB
        DXIDX(1,2)=(DXDXI(1,3)*DXDXI(3,2)-DXDXI(1,2)*DXDXI(3,3))/DJACB
        DXIDX(1,3)=(DXDXI(1,2)*DXDXI(2,3)-DXDXI(2,2)*DXDXI(1,3))/DJACB
        DXIDX(2,1)=(DXDXI(3,1)*DXDXI(2,3)-DXDXI(2,1)*DXDXI(3,3))/DJACB
        DXIDX(2,2)=(DXDXI(1,1)*DXDXI(3,3)-DXDXI(1,3)*DXDXI(3,1))/DJACB
        DXIDX(2,3)=(DXDXI(2,1)*DXDXI(1,3)-DXDXI(1,1)*DXDXI(2,3))/DJACB
        DXIDX(3,1)=(DXDXI(2,1)*DXDXI(3,2)-DXDXI(3,1)*DXDXI(2,2))/DJACB
        DXIDX(3,2)=(DXDXI(3,1)*DXDXI(1,2)-DXDXI(1,1)*DXDXI(3,2))/DJACB
        DXIDX(3,3)=(DXDXI(1,1)*DXDXI(2,2)-DXDXI(1,2)*DXDXI(2,1))/DJACB
      ENDIF
	DO 10 IDIM=1,NDIM
	DO 10 INODE=1,NODE
	DNDX(IDIM,INODE)=0.0D0
	DO 10 IXI=1,NDIM
	DNDX(IDIM,INODE)=DNDX(IDIM,INODE)+DXIDX(IDIM,IXI)*DNDXI(IXI,INODE)
10	CONTINUE
600   FORMAT(//,' Program halted in subroutine jacob :',/,8X,
     &          ' zero or negative area with jacob =',E14.6,/,8X,
     &          ' element number =',I7,/,5X,
     &          ' plotting file named "error_element.nas"')
	RETURN
      END

      SUBROUTINE VAL_DNDXX(IETYPE,NDIM,NODE,XI,ELCOD,DNDXI,
     &                     DXDXI,DJACB,DXIDX,DNDXX)
      IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION XI(3),DNDXI(3,NODE),ELCOD(3,NODE),DXIDX(3,3),DXDXI(3,3) ! INPUT
     &,DNDXX(3,3,NODE) ! OUTPUT
      DIMENSION DNDXIDXI(3,3,NODE),DXDXIDXI(3,3,NDIM),DXIDXDXI(3,3,3),
     &          DSH1D(NODE)  ! WORK ARRAY
      
      ! EVALUATE DN/(DXIKDXIL) --> DNDXIDXI
      CALL DN_DXI_DXI(IETYPE,NODE,XI,DNDXIDXI)
      ! EVALUATE DX/(DXIKDXIL) --> DXDXIDXI
      DO IDIM=1,NDIM
       DXDXIDXI(1:NDIM,1:NDIM,IDIM)=0.D0
       DO INODE=1,NODE
        DXDXIDXI(1:NDIM,1:NDIM,IDIM)=DXDXIDXI(1:NDIM,1:NDIM,IDIM)+
     &           DNDXIDXI(1:NDIM,1:NDIM,INODE)*ELCOD(IDIM,INODE)
       ENDDO
      ENDDO
      ! EVALUATE D(DXIK/DXI)/DXIL --> DXIDXDXI
      CALL EVAL_DXIDXDXI(NDIM,NODE,DXDXI,DXDXIDXI,DXIDXDXI,C_JACOBIAN)
      ! FORM DN/(DXIDXJ) --> DNDXX
      DO ID=1,NODE
       DO J=1,NDIM
        DO I=1,NDIM
         DNDXX(I,J,ID)=0.D0
         DO K=1,NDIM
          DO L=1,NDIM
           DNDXX(I,J,ID)=DNDXX(I,J,ID)+(DNDXIDXI(K,L,ID)*DXIDX(I,K)+
     &                   DNDXI(K,ID)*DXIDXDXI(L,I,K))*DXIDX(J,L)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE EVAL_DXIDXDXI(NDIM,NODE,DXDXI,DXDXIDXI,DXIDXDXI,
     &                         C_JACOBIAN)
      IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION DXDXI(3,3),DXDXIDXI(3,3,NDIM),DXIDXDXI(3,3,3),
     &	      A_MATRIX(3,3),DA_MATRIX(3,3,3),DC_JACOBIAN(3)
	! DXIDXDXI(L,I,K), DXDXI(IXI,IDIM),DXDXIDXI(IXI,IET,IDIM)
      IF(NDIM.EQ.1) THEN      ! 1D CASE
       DXIDXDXI(1,1,1)=-DXDXIDXI(1,1,1)/DXDXI(1,1)*DXDXI(1,1)
      ELSEIF(NDIM.EQ.2) THEN  ! 2D CASE
       C_JACOBIAN=DXDXI(1,1)*DXDXI(2,2)-DXDXI(1,2)*DXDXI(2,1)
       DXIDXDXI(1,1,1)=-DXDXIDXI(1,1,1)*DXDXI(2,2)*DXDXI(2,2)+
     &                  DXDXIDXI(1,2,1)*DXDXI(1,2)*DXDXI(2,2)-
     &                  DXDXIDXI(1,2,2)*DXDXI(2,1)*DXDXI(1,2)+
     &                  DXDXIDXI(1,1,2)*DXDXI(2,1)*DXDXI(2,2)
       DXIDXDXI(1,1,2)=DXDXIDXI(1,1,1)*DXDXI(1,2)*DXDXI(2,2)-
     &                 DXDXIDXI(1,2,1)*DXDXI(1,2)*DXDXI(1,2)+
     &                 DXDXIDXI(1,2,2)*DXDXI(1,1)*DXDXI(1,2)-
     &                 DXDXIDXI(1,1,2)*DXDXI(1,1)*DXDXI(2,2)
       DXIDXDXI(1,2,1)=-DXDXIDXI(1,2,1)*DXDXI(1,1)*DXDXI(2,2)+
     &                  DXDXIDXI(1,1,1)*DXDXI(2,1)*DXDXI(2,2)+
     &                  DXDXIDXI(1,2,2)*DXDXI(1,1)*DXDXI(2,1)-
     &                  DXDXIDXI(1,1,2)*DXDXI(2,1)*DXDXI(2,1)
       DXIDXDXI(1,2,2)=DXDXIDXI(1,2,1)*DXDXI(1,1)*DXDXI(1,2)-
     &                 DXDXIDXI(1,1,1)*DXDXI(2,1)*DXDXI(1,2)-
     &                 DXDXIDXI(1,2,2)*DXDXI(1,1)*DXDXI(1,1)+
     &                 DXDXIDXI(1,1,2)*DXDXI(1,1)*DXDXI(2,1)

       DXIDXDXI(2,1,1)=-DXDXIDXI(1,2,1)*DXDXI(2,2)*DXDXI(2,2)+
     &                  DXDXIDXI(2,2,1)*DXDXI(1,2)*DXDXI(2,2)-
     &                  DXDXIDXI(2,2,2)*DXDXI(2,1)*DXDXI(1,2)+
     &                  DXDXIDXI(1,2,2)*DXDXI(2,1)*DXDXI(2,2)
       DXIDXDXI(2,1,2)=DXDXIDXI(1,2,1)*DXDXI(1,2)*DXDXI(2,2)+
     &                 DXDXIDXI(2,2,2)*DXDXI(1,1)*DXDXI(1,2)-
     &                 DXDXIDXI(2,2,1)*DXDXI(1,2)*DXDXI(1,2)-
     &                 DXDXIDXI(1,2,2)*DXDXI(1,1)*DXDXI(2,2)
       DXIDXDXI(2,2,1)=-DXDXIDXI(2,2,1)*DXDXI(1,1)*DXDXI(2,2)+
     &                  DXDXIDXI(1,2,1)*DXDXI(2,1)*DXDXI(2,2)+
     &                  DXDXIDXI(2,2,2)*DXDXI(1,1)*DXDXI(2,1)-
     &                  DXDXIDXI(1,2,2)*DXDXI(2,1)*DXDXI(2,1)
       DXIDXDXI(2,2,2)=-DXDXIDXI(2,2,2)*DXDXI(1,1)*DXDXI(1,1)+
     &                  DXDXIDXI(2,2,1)*DXDXI(1,1)*DXDXI(1,2)-
     &                  DXDXIDXI(1,2,1)*DXDXI(2,1)*DXDXI(1,2)+
     &                  DXDXIDXI(1,2,2)*DXDXI(1,1)*DXDXI(2,1)
       DXIDXDXI(1:2,1:2,1:2)=DXIDXDXI(1:2,1:2,1:2)/C_JACOBIAN/C_JACOBIAN
      ELSE   ! 3D CASE      
       C_JACOBIAN = DXDXI(1,1)*DXDXI(2,2)*DXDXI(3,3)-
     &              DXDXI(2,1)*DXDXI(1,2)*DXDXI(3,3)-
     &              DXDXI(1,1)*DXDXI(3,2)*DXDXI(2,3)+
     &              DXDXI(3,1)*DXDXI(1,2)*DXDXI(2,3)+
     &              DXDXI(2,1)*DXDXI(3,2)*DXDXI(1,3)-
     &              DXDXI(3,1)*DXDXI(2,2)*DXDXI(1,3)

       A_MATRIX(1,1) = DXDXI(2,2)*DXDXI(3,3) - DXDXI(3,2)*DXDXI(2,3)
	 A_MATRIX(1,2) = DXDXI(3,2)*DXDXI(1,3) - DXDXI(1,2)*DXDXI(3,3)
	 A_MATRIX(1,3) = DXDXI(1,2)*DXDXI(2,3) - DXDXI(2,2)*DXDXI(1,3)
	 A_MATRIX(2,1) = DXDXI(3,1)*DXDXI(2,3) - DXDXI(2,1)*DXDXI(3,3)
	 A_MATRIX(2,2) = DXDXI(1,1)*DXDXI(3,3) - DXDXI(3,1)*DXDXI(1,3)
	 A_MATRIX(2,3) = DXDXI(2,1)*DXDXI(1,3) - DXDXI(1,1)*DXDXI(2,3)
	 A_MATRIX(3,1) = DXDXI(2,1)*DXDXI(3,2) - DXDXI(3,1)*DXDXI(2,2)
	 A_MATRIX(3,2) = DXDXI(3,1)*DXDXI(1,2) - DXDXI(1,1)*DXDXI(3,2)
	 A_MATRIX(3,3) = DXDXI(1,1)*DXDXI(2,2) - DXDXI(2,1)*DXDXI(1,2)

!      AIJ,K -> AIK,L
!  K=1 -> L=1

       DA_MATRIX(1,1,1) =  DXDXIDXI(2,1,2)*DXDXI(3,3) + 
     &                     DXDXIDXI(3,1,3)*DXDXI(2,2) -
     &                     DXDXIDXI(3,1,2)*DXDXI(2,3) -
     &                     DXDXIDXI(2,1,3)*DXDXI(3,2)
 
       DA_MATRIX(1,2,1) =  DXDXIDXI(3,1,2)*DXDXI(1,3) +
     &                     DXDXIDXI(1,1,3)*DXDXI(3,2) -
     &                     DXDXIDXI(1,1,2)*DXDXI(3,3) -
     &                     DXDXIDXI(3,1,3)*DXDXI(1,2)

       DA_MATRIX(1,3,1) =  DXDXIDXI(1,1,2)*DXDXI(2,3) +
     &                     DXDXIDXI(2,1,3)*DXDXI(1,2) -
     &                     DXDXIDXI(2,1,2)*DXDXI(1,3) -
     &                     DXDXIDXI(1,1,3)*DXDXI(2,2)

       DA_MATRIX(2,1,1) =  DXDXIDXI(3,1,1)*DXDXI(2,3) +
     &                     DXDXIDXI(2,1,3)*DXDXI(3,1) -
     &                     DXDXIDXI(2,1,1)*DXDXI(3,3) -
     &                     DXDXIDXI(3,1,3)*DXDXI(2,1)

       DA_MATRIX(2,2,1) =  DXDXIDXI(1,1,1)*DXDXI(3,3) +
     &                     DXDXIDXI(3,1,3)*DXDXI(1,1) -
     &                     DXDXIDXI(3,1,1)*DXDXI(1,3) -
     &                     DXDXIDXI(1,1,3)*DXDXI(3,1)

       DA_MATRIX(2,3,1) =  DXDXIDXI(2,1,1)*DXDXI(1,3) +
     &                     DXDXIDXI(1,1,3)*DXDXI(2,1) -
     &                     DXDXIDXI(1,1,1)*DXDXI(2,3) -
     &                     DXDXIDXI(2,1,3)*DXDXI(1,1)

       DA_MATRIX(3,1,1) =  DXDXIDXI(2,1,1)*DXDXI(3,2) +
     &                     DXDXIDXI(3,1,2)*DXDXI(2,1) -
     &                     DXDXIDXI(3,1,1)*DXDXI(2,2) -
     &                     DXDXIDXI(2,1,2)*DXDXI(3,1)

       DA_MATRIX(3,2,1) =  DXDXIDXI(3,1,1)*DXDXI(1,2) +
     &                     DXDXIDXI(1,1,2)*DXDXI(3,1) -
     &                     DXDXIDXI(1,1,1)*DXDXI(3,2) -
     &                     DXDXIDXI(3,1,2)*DXDXI(1,1)

       DA_MATRIX(3,3,1) =  DXDXIDXI(1,1,1)*DXDXI(2,2) +
     &                     DXDXIDXI(2,1,2)*DXDXI(1,1) -
     &                     DXDXIDXI(2,1,1)*DXDXI(1,2) -
     &                     DXDXIDXI(1,1,2)*DXDXI(2,1)

!  K=2 -> L=2

       DA_MATRIX(1,1,2) =  DXDXIDXI(2,2,2)*DXDXI(3,3) +
     &                     DXDXIDXI(3,2,3)*DXDXI(2,2) -
     &                     DXDXIDXI(3,2,2)*DXDXI(2,3) -
     &                     DXDXIDXI(2,2,3)*DXDXI(3,2)

       DA_MATRIX(1,2,2) =  DXDXIDXI(3,2,2)*DXDXI(1,3) +
     &                     DXDXIDXI(1,2,3)*DXDXI(3,2) -
     &                     DXDXIDXI(1,2,2)*DXDXI(3,3) -
     &                     DXDXIDXI(3,2,3)*DXDXI(1,2)

       DA_MATRIX(1,3,2) =  DXDXIDXI(1,2,2)*DXDXI(2,3) +
     &                     DXDXIDXI(2,2,3)*DXDXI(1,2) -
     &                     DXDXIDXI(2,2,2)*DXDXI(1,3) -
     &                     DXDXIDXI(1,2,3)*DXDXI(2,2)


       DA_MATRIX(2,1,2) =  DXDXIDXI(3,2,1)*DXDXI(2,3) +
     &                     DXDXIDXI(2,2,3)*DXDXI(3,1) -
     &                     DXDXIDXI(2,2,1)*DXDXI(3,3) -
     &                     DXDXIDXI(3,2,3)*DXDXI(2,1)

       DA_MATRIX(2,2,2) =  DXDXIDXI(1,2,1)*DXDXI(3,3) +
     &                     DXDXIDXI(3,2,3)*DXDXI(1,1) -
     &                     DXDXIDXI(3,2,1)*DXDXI(1,3) -
     &                     DXDXIDXI(1,2,3)*DXDXI(3,1)

       DA_MATRIX(2,3,2) =  DXDXIDXI(2,2,1)*DXDXI(1,3) +
     &                     DXDXIDXI(1,2,3)*DXDXI(2,1) -
     &                     DXDXIDXI(1,2,1)*DXDXI(2,3) -
     &                     DXDXIDXI(2,2,3)*DXDXI(1,1)

       DA_MATRIX(3,1,2) =  DXDXIDXI(2,2,1)*DXDXI(3,2) +
     &                     DXDXIDXI(3,2,2)*DXDXI(2,1) -
     &                     DXDXIDXI(3,2,1)*DXDXI(2,2) -
     &                     DXDXIDXI(2,2,2)*DXDXI(3,1)

       DA_MATRIX(3,2,2) =  DXDXIDXI(3,2,1)*DXDXI(1,2) +
     &                     DXDXIDXI(1,2,2)*DXDXI(3,1) -
     &                     DXDXIDXI(1,2,1)*DXDXI(3,2) -
     &                     DXDXIDXI(3,2,2)*DXDXI(1,1)

       DA_MATRIX(3,3,2) =  DXDXIDXI(1,2,1)*DXDXI(2,2) +
     &                     DXDXIDXI(2,2,2)*DXDXI(1,1) -
     &                     DXDXIDXI(2,2,1)*DXDXI(1,2) -
     &                     DXDXIDXI(1,2,2)*DXDXI(2,1)

!   K=3 -> L=3

       DA_MATRIX(1,1,3) =  DXDXIDXI(2,3,2)*DXDXI(3,3) + 
     &                     DXDXIDXI(3,3,3)*DXDXI(2,2) -
     &                     DXDXIDXI(3,3,2)*DXDXI(2,3) -
     &                     DXDXIDXI(2,3,3)*DXDXI(3,2)

       DA_MATRIX(1,2,3) =  DXDXIDXI(3,3,2)*DXDXI(1,3) +
     &                     DXDXIDXI(1,3,3)*DXDXI(3,2) -
     &                     DXDXIDXI(1,3,2)*DXDXI(3,3) -
     &                     DXDXIDXI(3,3,3)*DXDXI(1,2)


       DA_MATRIX(1,3,3) =  DXDXIDXI(1,3,2)*DXDXI(2,3) +
     &                     DXDXIDXI(2,3,3)*DXDXI(1,2) -
     &                     DXDXIDXI(2,3,2)*DXDXI(1,3) -
     &                     DXDXIDXI(1,3,3)*DXDXI(2,2)

       DA_MATRIX(2,1,3) =  DXDXIDXI(3,3,1)*DXDXI(2,3) +
     &                     DXDXIDXI(2,3,3)*DXDXI(3,1) -
     &                     DXDXIDXI(2,3,1)*DXDXI(3,3) -
     &                     DXDXIDXI(3,3,3)*DXDXI(2,1)

       DA_MATRIX(2,2,3) =  DXDXIDXI(1,3,1)*DXDXI(3,3) +
     &                     DXDXIDXI(3,3,3)*DXDXI(1,1) -
     &                     DXDXIDXI(3,3,1)*DXDXI(1,3) -
     &                     DXDXIDXI(1,3,3)*DXDXI(3,1)

       DA_MATRIX(2,3,3) =  DXDXIDXI(2,3,1)*DXDXI(1,3) +
     &                     DXDXIDXI(1,3,3)*DXDXI(2,1) -
     &                     DXDXIDXI(1,3,1)*DXDXI(2,3) -
     &                     DXDXIDXI(2,3,3)*DXDXI(1,1)

       DA_MATRIX(3,1,3) =  DXDXIDXI(2,3,1)*DXDXI(3,2) +
     &                     DXDXIDXI(3,3,2)*DXDXI(2,1) -
     &                     DXDXIDXI(3,3,1)*DXDXI(2,2) -
     &                     DXDXIDXI(2,3,2)*DXDXI(3,1)

       DA_MATRIX(3,2,3) =  DXDXIDXI(3,3,1)*DXDXI(1,2) +
     &                     DXDXIDXI(1,3,2)*DXDXI(3,1) -
     &                     DXDXIDXI(1,3,1)*DXDXI(3,2) -
     &                     DXDXIDXI(3,3,2)*DXDXI(1,1)

       DA_MATRIX(3,3,3) =  DXDXIDXI(1,3,1)*DXDXI(2,2) +
     &                     DXDXIDXI(2,3,2)*DXDXI(1,1) -
     &                     DXDXIDXI(2,3,1)*DXDXI(1,2) -
     &                     DXDXIDXI(1,3,2)*DXDXI(2,1)

!! FOR DC_JACOBIAN

       DC_JACOBIAN(1) = DXDXIDXI(1,1,1)*DXDXI(2,2)*DXDXI(3,3) + 
     &                  DXDXIDXI(2,1,2)*DXDXI(1,1)*DXDXI(3,3) +  
     &                  DXDXIDXI(3,1,3)*DXDXI(1,1)*DXDXI(2,2) -
     
     &                  DXDXIDXI(2,1,1)*DXDXI(1,2)*DXDXI(3,3) -  
     &                  DXDXIDXI(1,1,2)*DXDXI(2,1)*DXDXI(3,3) - 
     &                  DXDXIDXI(3,1,3)*DXDXI(2,1)*DXDXI(1,2) -
     
     &                  DXDXIDXI(1,1,1)*DXDXI(3,2)*DXDXI(2,3) -
     &                  DXDXIDXI(3,1,2)*DXDXI(1,1)*DXDXI(2,3) -
     &                  DXDXIDXI(2,1,3)*DXDXI(1,1)*DXDXI(3,2) +
     
     &                  DXDXIDXI(3,1,1)*DXDXI(1,2)*DXDXI(2,3) +
     &                  DXDXIDXI(1,1,2)*DXDXI(3,1)*DXDXI(2,3) +
     &                  DXDXIDXI(2,1,3)*DXDXI(3,1)*DXDXI(1,2) +
     
     &                  DXDXIDXI(2,1,1)*DXDXI(3,2)*DXDXI(1,3) +
     &                  DXDXIDXI(3,1,2)*DXDXI(2,1)*DXDXI(1,3) +
     &                  DXDXIDXI(1,1,3)*DXDXI(2,1)*DXDXI(3,2) -
     
     &                  DXDXIDXI(3,1,1)*DXDXI(2,2)*DXDXI(1,3) -
     &                  DXDXIDXI(2,1,2)*DXDXI(3,1)*DXDXI(1,3) -
     &                  DXDXIDXI(1,1,3)*DXDXI(3,1)*DXDXI(2,2)
     
       DC_JACOBIAN(2) = DXDXIDXI(1,2,1)*DXDXI(2,2)*DXDXI(3,3) + 
     &                  DXDXIDXI(2,2,2)*DXDXI(1,1)*DXDXI(3,3) +  
     &                  DXDXIDXI(3,2,3)*DXDXI(1,1)*DXDXI(2,2) -
     
     &                  DXDXIDXI(2,2,1)*DXDXI(1,2)*DXDXI(3,3) -  
     &                  DXDXIDXI(1,2,2)*DXDXI(2,1)*DXDXI(3,3) - 
     &                  DXDXIDXI(3,2,3)*DXDXI(2,1)*DXDXI(1,2) -
     
     &                  DXDXIDXI(1,2,1)*DXDXI(3,2)*DXDXI(2,3) -
     &                  DXDXIDXI(3,2,2)*DXDXI(1,1)*DXDXI(2,3) -
     &                  DXDXIDXI(2,2,3)*DXDXI(1,1)*DXDXI(3,2) +
     
     &                  DXDXIDXI(3,2,1)*DXDXI(1,2)*DXDXI(2,3) +
     &                  DXDXIDXI(1,2,2)*DXDXI(3,1)*DXDXI(2,3) +
     &                  DXDXIDXI(2,2,3)*DXDXI(3,1)*DXDXI(1,2) +
     
     &                  DXDXIDXI(2,2,1)*DXDXI(3,2)*DXDXI(1,3) +
     &                  DXDXIDXI(3,2,2)*DXDXI(2,1)*DXDXI(1,3) +
     &                  DXDXIDXI(1,2,3)*DXDXI(2,1)*DXDXI(3,2) -
     
     &                  DXDXIDXI(3,2,1)*DXDXI(2,2)*DXDXI(1,3) -
     &                  DXDXIDXI(2,2,2)*DXDXI(3,1)*DXDXI(1,3) -
     &                  DXDXIDXI(1,2,3)*DXDXI(3,1)*DXDXI(2,2)
     
     
       DC_JACOBIAN(3) = DXDXIDXI(1,3,1)*DXDXI(2,2)*DXDXI(3,3) + 
     &                  DXDXIDXI(2,3,2)*DXDXI(1,1)*DXDXI(3,3) +  
     &                  DXDXIDXI(3,3,3)*DXDXI(1,1)*DXDXI(2,2) -
     
     &                  DXDXIDXI(2,3,1)*DXDXI(1,2)*DXDXI(3,3) -  
     &                  DXDXIDXI(1,3,2)*DXDXI(2,1)*DXDXI(3,3) - 
     &                  DXDXIDXI(3,3,3)*DXDXI(2,1)*DXDXI(1,2) -
     
     &                  DXDXIDXI(1,3,1)*DXDXI(3,2)*DXDXI(2,3) -
     &                  DXDXIDXI(3,3,2)*DXDXI(1,1)*DXDXI(2,3) -
     &                  DXDXIDXI(2,3,3)*DXDXI(1,1)*DXDXI(3,2) +
     
     &                  DXDXIDXI(3,3,1)*DXDXI(1,2)*DXDXI(2,3) +
     &                  DXDXIDXI(1,3,2)*DXDXI(3,1)*DXDXI(2,3) +
     &                  DXDXIDXI(2,3,3)*DXDXI(3,1)*DXDXI(1,2) +
     
     &                  DXDXIDXI(2,3,1)*DXDXI(3,2)*DXDXI(1,3) +
     &                  DXDXIDXI(3,3,2)*DXDXI(2,1)*DXDXI(1,3) +
     &                  DXDXIDXI(1,3,3)*DXDXI(2,1)*DXDXI(3,2) -
     
     &                  DXDXIDXI(3,3,1)*DXDXI(2,2)*DXDXI(1,3) -
     &                  DXDXIDXI(2,3,2)*DXDXI(3,1)*DXDXI(1,3) -
     &                  DXDXIDXI(1,3,3)*DXDXI(3,1)*DXDXI(2,2)
 

!! CHECK THE MATRIX ELEMENT ORDER

       DXIDXDXI(1,1:3,1:3) = DA_MATRIX(1:3,1:3,1)/C_JACOBIAN - 
     &                       A_MATRIX(1:3,1:3)*DC_JACOBIAN(1) /
     &                       C_JACOBIAN**2

       DXIDXDXI(2,1:3,1:3) = DA_MATRIX(1:3,1:3,2)/C_JACOBIAN - 
     &                       A_MATRIX(1:3,1:3)*DC_JACOBIAN(2) /
     &                       C_JACOBIAN**2

       DXIDXDXI(3,1:3,1:3) = DA_MATRIX(1:3,1:3,3)/C_JACOBIAN - 
     &                       A_MATRIX(1:3,1:3)*DC_JACOBIAN(3) /
     &                       C_JACOBIAN**2

!       DO I=1,3
!        DXIDXDXI(1,I,1:3) = DA_MATRIX(I,1:3,1)/C_JACOBIAN - 
!     &                        A_MATRIX(I,1:3)*DC_JACOBIAN(1) /
!     &                        C_JACOBIAN**2
!
!        DXIDXDXI(2,I,1:3) = DA_MATRIX(I,1:3,2)/C_JACOBIAN - 
!     &                        A_MATRIX(I,1:3)*DC_JACOBIAN(2) /
!     &                        C_JACOBIAN**2
!
!        DXIDXDXI(3,I,1:3) = DA_MATRIX(I,1:3,3)/C_JACOBIAN - 
!     &                        A_MATRIX(I,1:3)*DC_JACOBIAN(3) /
!     &                        C_JACOBIAN**2
!       ENDDO

      ENDIF
      RETURN
      END
      
      SUBROUTINE DSHAPEF(NDIM,NDIMB,IETYPE,NODE,X,CK,COSB,FJCB,INFO)
      IMPLICIT REAL*8 (A-H,O-Z)
      ! IFCOS=2,ONLY FIND DX/DXI 
      ! IFCOS=0,FIND JACOB AND DX/DXI
      ! IFCOS=1,COMPUTE JACOB,DX/DXI AND OUTWARD NORMAL
      DIMENSION X(3),CK(3,NODE),DN(3,NODE),GD(3,3),COSB(NDIM),GR(3)
      ! COMPUTE DN_DXI
      CALL DN_DXI(IETYPE,NODE,X,DN) 
      ! COMPUTE JACOBIAN AND COSB
      DO I=1,NDIM
        DO J=1,NDIMB
          GD(I,J)=0.D0
          DO ID=1,NODE
            GD(I,J)=GD(I,J)+DN(J,ID)*CK(I,ID)
          ENDDO
        ENDDO
      ENDDO
      IF(NDIMB.EQ.1) THEN
        GR(1)=GD(2,1)
        GR(2)=-GD(1,1) ! FOR 2D BOUNDARY NORMALS
        FJCB=DSQRT(DOT_PRODUCT(GD(1:NDIM,1),GD(1:NDIM,1))) ! LINE JACOBIAN
      ELSE
        GR(1)=GD(2,1)*GD(3,2)-GD(3,1)*GD(2,2) ! FOR 3D SURFACE NORMALS
        GR(2)=GD(3,1)*GD(1,2)-GD(1,1)*GD(3,2)
        GR(3)=GD(1,1)*GD(2,2)-GD(2,1)*GD(1,2)
        FJCB=DSQRT(GR(1)*GR(1)+GR(2)*GR(2)+GR(3)*GR(3)) ! SURFACE JACOBIAN
      ENDIF
      COSB(1:NDIM)=GR(1:NDIM)/FJCB ! NORMAL COSINE
      IF(DABS(FJCB).LT.1.0D-10) INFO=1
      RETURN
      END
      
      SUBROUTINE NODE_XI(INODE,IETYPE,XI)
      ! DETERMINE XI FOR NODE INODE OF THE IETYPE ELEMENT
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XI(3)
      SELECT CASE(IETYPE)
      ! LINE ELEMENTS
      CASE(2) ! 2 NODED LINE ELEMENT
        XI(1)=-1.D0+DBLE(INODE-1)*2.D0
      CASE(3) ! 3 NODED LINE ELEMENT
        XI(1)=-1.D0+DBLE(INODE-1)*2.D0/2.D0
      CASE(4) ! 4 NODED LINE ELEMENT
        XI(1)=-1.D0+DBLE(INODE-1)*2.D0/3.D0
      CASE(5) ! 5 NODED LINE ELEMENT 
	  XI(1)=-1.D0+DBLE(INODE-1)*2.D0/4.D0
        
      ! PLATE ELEMENTS
      CASE(304) ! 4 NODE RECTANGULAR ELEMENT
	  SELECT CASE (INODE)
        CASE(1)
          XI(1)=-1.D0
          XI(2)=-1.D0
        CASE(2)
          XI(1)=1.D0
          XI(2)=-1.D0
        CASE(3)
          XI(1)=1.D0
          XI(2)=1.D0
        CASE(4)
          XI(1)=-1.D0
          XI(2)=1.D0
        END SELECT
      CASE(306) ! 6 NODE TRIANGULAR ELEMENT (ONLY FOR 11 NODE ELEMENT FACE)
        SELECT CASE (INODE)
        CASE(1)
          XI(1)=-1.D0
          XI(2)=-DSQRT(0.3D1) / 0.3D1
        CASE(2)
          XI(1)=1.D0
          XI(2)=-DSQRT(0.3D1) / 0.3D1
        CASE(3)
          XI(1)=0.D0
          XI(2)=0.2D1 / 0.3D1 * DSQRT(0.3D1)
        CASE(4)
          XI(1)=0.D0
          XI(2)=-DSQRT(0.3D1) / 0.3D1
        CASE(5)
          XI(1)=0.1D1 / 0.2D1
          XI(2)=DSQRT(0.3D1) / 0.6D1
        CASE(6)
          XI(1)=-0.1D1 / 0.2D1
          XI(2)=DSQRT(0.3D1) / 0.6D1
        END SELECT
      CASE(307) ! 7 NODE TRIANGULAR ELEMENT
        SELECT CASE (INODE)
        CASE(1)
          XI(1)=-1.D0
          XI(2)=-DSQRT(0.3D1) / 0.3D1
        CASE(2)
          XI(1)=1.D0
          XI(2)=-DSQRT(0.3D1) / 0.3D1
        CASE(3)
          XI(1)=0.D0
          XI(2)=0.2D1 / 0.3D1 * DSQRT(0.3D1)
        CASE(4)
          XI(1)=0.D0
          XI(2)=-DSQRT(0.3D1) / 0.3D1
        CASE(5)
          XI(1)=0.1D1 / 0.2D1
          XI(2)=DSQRT(0.3D1) / 0.6D1
        CASE(6)
          XI(1)=-0.1D1 / 0.2D1
          XI(2)=DSQRT(0.3D1) / 0.6D1
        CASE(7)
          XI(1)=0.D0
          XI(2)=0.D0
        END SELECT
      CASE(308) ! 8 NODE RECTANGULAR ELEMENT
	  SELECT CASE (INODE)
        CASE(1)
          XI(1)=-1.D0
          XI(2)=-1.D0
        CASE(2)
          XI(1)=1.D0
          XI(2)=-1.D0
        CASE(3)
          XI(1)=1.D0
          XI(2)=1.D0
        CASE(4)
          XI(1)=-1.D0
          XI(2)=1.D0
        CASE(5)
          XI(1)=0.D0
          XI(2)=-1.D0
        CASE(6)
          XI(1)=1.D0
          XI(2)=0.D0
        CASE(7)
          XI(1)=0.D0
          XI(2)=1.D0
        CASE(8)
          XI(1)=-1.D0
          XI(2)=0.D0
        END SELECT
      CASE(309) ! 9 NODE QUADRILATERAL ELEMENT
        ! FORM NX=3, NY=3
        ! NODAL VALUES OF KSI = {-1, 0, 1}
        ! NODAL VALUES OF ETA = {-1, 0, 1} 
        KY0=(INODE-1)/3
        KX=INODE-KY0*3
        XI(1)=-1.D0+2.D0/2.D0*(KX-1.D0)
        XI(2)=-1.D0+2.D0/2.D0*KY0
      CASE(316) ! 16 NODE RECTAGULAR SURFACE ELEMENT
        ! FORM NX=4, NY=4
        ! NODAL VALUES OF KSI = {-1, -1/3, 1/3, 1}
        ! NODAL VALUES OF ETA = {-1, -1/3, 1/3, 1} 
        KY0=(INODE-1)/4
        KX=INODE-KY0*4
        XI(1)=-1.D0+2.D0/3.D0*(KX-1.D0)
        XI(2)=-1.D0+2.D0/3.D0*KY0
      CASE(325) ! 25 NODE RECTAGULAR SURFACE ELEMENT
        ! FORM NX=5, NY=5
        ! NODAL VALUES OF KSI = {-1, -1/2, 0, 1/2, 1}
        ! NODAL VALUES OF ETA = {-1, -1/2, 0, 1/2, 1} 
        KY0=(INODE-1)/5
        KX=INODE-KY0*5
        XI(1)=-1.D0+2.D0/4.D0*(KX-1.D0)
        XI(2)=-1.D0+2.D0/4.D0*KY0

      ! BRICK ELEMENTS
      CASE(408) ! 8 NODED BRICK ELEMENT
        SELECT CASE (INODE)
        CASE(1)
          XI(1:3)=(/-1.D0,-1.D0,-1.D0/)
        CASE(2)
          XI(1:3)=(/1.D0,-1.D0,-1.D0/)
        CASE(3)
          XI(1:3)=(/1.D0,1.D0,-1.D0/)
        CASE(4)
          XI(1:3)=(/-1.D0, 1.D0,-1.D0/)
        CASE(5)
          XI(1:3)=(/-1.D0,-1.D0,1.D0/)
        CASE(6)
          XI(1:3)=(/1.D0,-1.D0,1.D0/)
        CASE(7)
          XI(1:3)=(/1.D0,1.D0,1.D0/)
        CASE(8)
          XI(1:3)=(/-1.D0,1.D0,1.D0/)
        END SELECT
      CASE(411) ! 11 NODED  TETRAHEDRAL ELEMENT
        SELECT CASE (INODE)
        CASE(1)
          XI(1)=-1.D0
          XI(2)=-DSQRT(3.D0)/3.D0
          XI(3)=-DSQRT(6.D0)/6.D0
        CASE(2)
          XI(1)=1.D0
          XI(2)=-DSQRT(3.D0)/3.D0
          XI(3)=-DSQRT(6.D0)/6.D0
        CASE(3)
          XI(1)=0.D0
          XI(2)=2.D0*DSQRT(3.D0)/3.D0
          XI(3)=-DSQRT(6.D0)/6.D0
        CASE(4)
          XI(1)=0.D0
          XI(2)=0.D0
          XI(3)=DSQRT(6.D0)/2.D0
        CASE(5)
          XI(1)=0.D0
          XI(2)=-DSQRT(3.D0)/3.D0
          XI(3)=-DSQRT(6.D0)/6.D0
        CASE(6)
          XI(1)=0.5D0
          XI(2)=DSQRT(3.D0)/6.D0
          XI(3)=-DSQRT(6.D0)/6.D0
        CASE(7)
          XI(1)=-0.5D0
          XI(2)=DSQRT(3.D0)/6.D0
          XI(3)=-DSQRT(6.D0)/6.D0
        CASE(8)
          XI(1)=-0.5D0
          XI(2)=-DSQRT(3.D0)/6.D0
          XI(3)=DSQRT(6.D0)/6.D0
        CASE(9)
          XI(1)=0.5D0
          XI(2)=-DSQRT(3.D0)/6.D0
          XI(3)=DSQRT(6.D0)/6.D0
        CASE(10)
          XI(1)=0.D0
          XI(2)=DSQRT(3.D0)/3.D0
          XI(3)=DSQRT(6.D0)/6.D0
        CASE(11)
          XI(1)=0.D0
          XI(2)=0.D0
          XI(3)=0.D0
        END SELECT
      CASE(420,421) ! 20 OR 21 NODED BRICK ELEMENT
	  SELECT CASE (INODE)
        CASE(1)
          XI(1:3)=(/-1.D0,-1.D0,-1.D0/)
        CASE(2)
          XI(1:3)=(/1.D0,-1.D0,-1.D0/)
        CASE(3)
          XI(1:3)=(/1.D0,1.D0,-1.D0/)
        CASE(4)
          XI(1:3)=(/-1.D0,1.D0,-1.D0/)
        CASE(5)
          XI(1:3)=(/-1.D0,-1.D0,1.D0/)
        CASE(6)
          XI(1:3)=(/1.D0,-1.D0,1.D0/)
        CASE(7)
          XI(1:3)=(/1.D0,1.D0,1.D0/)
        CASE(8)
          XI(1:3)=(/-1.D0,1.D0,1.D0/)
        CASE(9)
          XI(1:3)=(/0.D0,-1.D0,-1.D0/)
        CASE(10)
          XI(1:3)=(/1.D0,0.D0,-1.D0/)
        CASE(11)
          XI(1:3)=(/0.D0,1.D0,-1.D0/)
        CASE(12)
          XI(1:3)=(/-1.D0,0.D0,-1.D0/)
        CASE(13)
          XI(1:3)=(/0.D0,-1.D0,1.D0/)
        CASE(14)
          XI(1:3)=(/1.D0,0.D0,1.D0/)
        CASE(15)
          XI(1:3)=(/0.D0,1.D0,1.D0/)
        CASE(16)
          XI(1:3)=(/-1.D0,0.D0,1.D0/)
        CASE(17)
          XI(1:3)=(/-1.D0,-1.D0,0.D0/)
        CASE(18)
          XI(1:3)=(/1.D0,-1.D0,0.D0/)
        CASE(19)
          XI(1:3)=(/1.D0,1.D0,0.D0/)
        CASE(20)
          XI(1:3)=(/-1.D0,1.D0,0.D0/)
        CASE(21)
          XI(1:3)=(/0.D0,0.D0,0.D0/)
        END SELECT
      CASE(427) ! 27 NODED BRICK ELEMENT
        ! FORM NX=3, NY=3, NZ=3
        ! NODAL VALUES OF KSI = {-1, 0, 1}
        ! NODAL VALUES OF ETA = {-1, 0, 1} 
        ! NODAL VALUES OF ZTA = {-1, 0, 1} 
        KZ0=(INODE-1)/9
        KY0=(INODE-KZ0*9-1)/3
        KX=INODE-KZ0*9-KY0*3
        XI(1)=-1.D0+2.D0/2.D0*(KX-1.D0)
        XI(2)=-1.D0+2.D0/2.D0*KY0
        XI(3)=-1.D0+2.D0/2.D0*KZ0
        CASE(464) ! 64 NODED BRICK ELEMENT
        ! FORM NX=4, NY=4, NZ=4
        ! NODAL VALUES OF KSI = {-1, -1/3, 1/3, 1}
        ! NODAL VALUES OF ETA = {-1, -1/3, 1/3, 1} 
        ! NODAL VALUES OF ZTA = {-1, -1/3, 1/3, 1} 
        KZ0=(INODE-1)/16
        KY0=(INODE-KZ0*16-1)/4
        KX=INODE-KZ0*16-KY0*4
        XI(1)=-1.D0+2.D0/3.D0*(KX-1.D0)
        XI(2)=-1.D0+2.D0/3.D0*KY0
        XI(3)=-1.D0+2.D0/3.D0*KZ0
        CASE(525) ! 125 NODED BRICK ELEMENT
        ! FORM NX=5, NY=5, NZ=5
        ! NODAL VALUES OF KSI = {-1, -1/2, 0, 1/2, 1}
        ! NODAL VALUES OF ETA = {-1, -1/2, 0, 1/2, 1} 
        ! NODAL VALUES OF ZTA = {-1, -1/2, 0, 1/2, 1} 
        KZ0=(INODE-1)/25
        KY0=(INODE-KZ0*25-1)/5
        KX=INODE-KZ0*25-KY0*5
        XI(1)=-1.D0+2.D0/4.D0*(KX-1.D0)
        XI(2)=-1.D0+2.D0/4.D0*KY0
        XI(3)=-1.D0+2.D0/4.D0*KZ0
        
      END SELECT
      END

      SUBROUTINE VALUE_SHAPE(IETYPE,NODE,X,SP)
      ! CALCULATE VALUES OF SHAPE FUNCTIONS
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SP(NODE),X(3)
      SELECT CASE(IETYPE)
      ! LINE ELEMENTS
      CASE(2) ! 2 NODE RECTANGULAR ELEMENT
       CALL SHAPE_1D_2ND(NODE,X(1),SP)
      CASE(3) ! 3 NODE RECTANGULAR ELEMENT
       CALL SHAPE_1D_3ND(NODE,X(1),SP)
      CASE(4) ! 4 NODE QUADRILATERAL ELEMENT
       CALL SHAPE_1D_4ND(NODE,X(1),SP)
      CASE(5) ! 5 NODE RECTAGULAR SURFACE ELEMENT
       CALL SHAPE_1D_5ND(NODE,X(1),SP)
        
      ! PLATE ELEMENTS
      CASE(304) ! 4 NODE RECTANGULAR ELEMENT
       CALL SHAPE_2D(2,2,NODE,X(1),X(2),SP)
      CASE(307) ! 7 NODE TRIANGULAR ELEMENT
       CALL SHAPE_2D_7ND_TRI(NODE,X(1),X(2),SP)
      CASE(308) ! 8 NODE RECTANGULAR ELEMENT
       CALL SHAPE_2D_8ND(NODE,X(1),X(2),SP)
      CASE(309) ! 9 NODE QUADRILATERAL ELEMENT
       CALL SHAPE_2D(3,3,NODE,X(1),X(2),SP)
      CASE(316) ! 16 NODE RECTAGULAR SURFACE ELEMENT
       CALL SHAPE_2D(4,4,NODE,X(1),X(2),SP)
      CASE(325) ! 25 NODE RECTAGULAR SURFACE ELEMENT
       CALL SHAPE_2D(5,5,NODE,X(1),X(2),SP)

      ! BRICK ELEMENTS
      CASE(408)    ! 8 NODED BRICK ELEMENT
       CALL SHAPE_3D(2,2,2,NODE,X(1),X(2),X(3),SP)
      CASE(411) ! 11 NODED  TETRAHEDRAL ELEMENT
       CALL SHAPE_3D_11ND_TRI(NODE,X(1),X(2),X(3),SP)
      CASE(420)    ! 20 NODED BRICK ELEMENT
       CALL SHAPE_3D_20ND(NODE,X(1),X(2),X(3),SP) 
      CASE(421)    ! 21 NODED BRICK ELEMENT
       CALL SHAPE_3D_21ND(NODE,X(1),X(2),X(3),SP)
      CASE(427)    ! 27 NODED BRICK ELEMENT
       CALL SHAPE_3D(3,3,3,NODE,X(1),X(2),X(3),SP)
      CASE(464)    ! 64 NODED BRICK ELEMENT
       CALL SHAPE_3D(4,4,4,NODE,X(1),X(2),X(3),SP)
      CASE(525)    ! 125 NODED BRICK ELEMENT
       CALL SHAPE_3D(5,5,5,NODE,X(1),X(2),X(3),SP)
      END SELECT  
      END
      
      SUBROUTINE SHAPE_2D(NX,NY,NODE,X,Y,SHAPEF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SHAPEF(NODE),SI(NX),SJ(NY)
      
      SELECT CASE(NX)
      CASE(2)
        CALL SHAPE_1D_2ND(NX,X,SI)
      CASE(3)
        CALL SHAPE_1D_3ND(NX,X,SI)
      CASE(4)
        CALL SHAPE_1D_4ND(NX,X,SI)
      CASE(5)
        CALL SHAPE_1D_5ND(NX,X,SI)
      END SELECT
      
      SELECT CASE(NY)
      CASE(2)
        CALL SHAPE_1D_2ND(NY,Y,SJ)
      CASE(3)
        CALL SHAPE_1D_3ND(NY,Y,SJ)
      CASE(4)
        CALL SHAPE_1D_4ND(NY,Y,SJ)
      CASE(5)
        CALL SHAPE_1D_5ND(NY,Y,SJ)
      END SELECT
      
      ND=0
      DO IETA=1,NY
        DO IKXI=1,NX
          ND=ND+1
          SHAPEF(ND)=SI(IKXI)*SJ(IETA)
        ENDDO
      ENDDO
      END

      SUBROUTINE SHAPE_3D(NX,NY,NZ,NODE,X,Y,Z,SHAP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SHAP(NODE),SI(NX),SJ(NY),SK(NZ)
      
      SELECT CASE(NX)
      CASE(2)
        CALL SHAPE_1D_2ND(NX,X,SI)
      CASE(3)
        CALL SHAPE_1D_3ND(NX,X,SI)
      CASE(4)
        CALL SHAPE_1D_4ND(NX,X,SI)
      CASE(5)
        CALL SHAPE_1D_5ND(NX,X,SI)
      END SELECT
      
      SELECT CASE(NY)
      CASE(2)
        CALL SHAPE_1D_2ND(NY,Y,SJ)
      CASE(3)
        CALL SHAPE_1D_3ND(NY,Y,SJ)
      CASE(4)
        CALL SHAPE_1D_4ND(NY,Y,SJ)
      CASE(5)
        CALL SHAPE_1D_5ND(NY,Y,SJ)
      END SELECT
      
      SELECT CASE(NZ)
      CASE(2)
        CALL SHAPE_1D_2ND(NZ,Z,SK)
      CASE(3)
        CALL SHAPE_1D_3ND(NZ,Z,SK)
      CASE(4)
        CALL SHAPE_1D_4ND(NZ,Z,SK)
      CASE(5)
        CALL SHAPE_1D_5ND(NZ,Z,SK)
      END SELECT
      
      ND=0
      DO K=1,NZ
        DO J=1,NY
          DO I=1,NX
            ND=ND+1
            SHAP(ND)=SI(I)*SJ(J)*SK(K)
          ENDDO
        ENDDO
      ENDDO
      END
      
      SUBROUTINE DN_DXI(IETYPE,NODE,X,DN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(3),DN(3,NODE),DSP1D(5)
      SELECT CASE(IETYPE)
      ! LINE ELEMENTS
      CASE(2) ! 2 NODE LINE ELEMENT
       CALL DSHAPE_1D_2ND(NODE,X(1),DSP1D)
       DN(1,1:2)=DSP1D(1:2)
      CASE(3) ! 3 NODE LINE ELEMENT
       CALL DSHAPE_1D_3ND(NODE,X(1),DSP1D)
       DN(1,1:3)=DSP1D(1:3)
      CASE(4) ! 4 NODE LINE ELEMENT
       CALL DSHAPE_1D_4ND(NODE,X(1),DSP1D)
       DN(1,1:4)=DSP1D(1:4)
      CASE(5) ! 5 NODE LINE SURFACE ELEMENT
       CALL DSHAPE_1D_5ND(NODE,X(1),DSP1D)
       DN(1,1:5)=DSP1D(1:5)
       
      ! PLATE ELEMENTS
      CASE(304) ! 4 NODED RECTAGULAR ELEMENT
       CALL DSHAPE_2D(2,2,NODE,X(1),X(2),DN)
      CASE(306) ! 6 NODED TRIANGULAR ELEMENT (ONLY FOR 11-NODE ELEMENT FACE)
       CALL DSHAPE_2D_6ND_TRI(NODE,X(1),X(2),DN)
      CASE(307) ! 7 NODE TRIANGULAR ELEMENT
       CALL DSHAPE_2D_7ND_TRI(NODE,X(1),X(2),DN)
      CASE(308) ! 8 NODED RECTAGULAR ELEMENT
       CALL DSHAPE_2D_8ND(NODE,X(1),X(2),DN)   
      CASE(309) ! 9 NODED RECTAGULAR ELEMENT
       CALL DSHAPE_2D(3,3,NODE,X(1),X(2),DN)
      CASE(316) ! 16 NODED RECTAGULAR ELEMENT
       CALL DSHAPE_2D(4,4,NODE,X(1),X(2),DN)
      CASE(325) ! 25 NODED RECTAGULAR ELEMENT
       CALL DSHAPE_2D(5,5,NODE,X(1),X(2),DN)

      ! BRICK ELEMENTS
      CASE(408) ! 8 NODED BRICK ELEMENT
        CALL DSHAPE_3D(2,2,2,NODE,X(1),X(2),X(3),DN)
      CASE(411) ! 11 NODED  TETRAHEDRAL ELEMENT
       CALL DSHAPE_3D_11ND_TRI(NODE,X(1),X(2),X(3),DN)
      CASE(420) ! 20 NODED BRICK ELEMENT
       CALL DSHAPE_3D_20ND(NODE,X(1),X(2),X(3),DN)
      CASE(421) ! 21 NODED BRICK ELEMENT
       CALL DSHAPE_3D_21ND(NODE,X(1),X(2),X(3),DN)
      CASE(427) ! 27 NODED BRICK ELEMENT
       CALL DSHAPE_3D(3,3,3,NODE,X(1),X(2),X(3),DN)
      CASE(464) ! 27 NODED BRICK ELEMENT
       CALL DSHAPE_3D(4,4,4,NODE,X(1),X(2),X(3),DN)
      CASE(525) ! 27 NODED BRICK ELEMENT
       CALL DSHAPE_3D(5,5,5,NODE,X(1),X(2),X(3),DN)
      END SELECT
      RETURN
      END

      SUBROUTINE DSHAPE_2D(NX,NY,NODE,X,Y,DN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DN(3,NODE),SI(NX),SJ(NY)
      DN = 0.0D0
      ! FIND DN/DXI
      SELECT CASE(NX)
      CASE(2)
        CALL DSHAPE_1D_2ND(NX,X,SI)
      CASE(3)
        CALL DSHAPE_1D_3ND(NX,X,SI)
      CASE(4)
        CALL DSHAPE_1D_4ND(NX,X,SI)
      CASE(5)
        CALL DSHAPE_1D_5ND(NX,X,SI)
      END SELECT
      SELECT CASE(NY)
      CASE(2)
        CALL SHAPE_1D_2ND(NY,Y,SJ)
      CASE(3)
        CALL SHAPE_1D_3ND(NY,Y,SJ)
      CASE(4)
        CALL SHAPE_1D_4ND(NY,Y,SJ)
      CASE(5)
        CALL SHAPE_1D_5ND(NY,Y,SJ)
      END SELECT
      ND=0
      DO J=1,NY
        DO I=1,NX
          ND=ND+1
          DN(1,ND)=SI(I)*SJ(J)
        ENDDO
      ENDDO
      ! FIND DN/DETA
      SELECT CASE(NX)
      CASE(2)
        CALL SHAPE_1D_2ND(NX,X,SI)
      CASE(3)
        CALL SHAPE_1D_3ND(NX,X,SI)
      CASE(4)
        CALL SHAPE_1D_4ND(NX,X,SI)
      CASE(5)
        CALL SHAPE_1D_5ND(NX,X,SI)
      END SELECT
      SELECT CASE(NY)
      CASE(2)
        CALL DSHAPE_1D_2ND(NY,Y,SJ)
      CASE(3)
        CALL DSHAPE_1D_3ND(NY,Y,SJ)
      CASE(4)
        CALL DSHAPE_1D_4ND(NY,Y,SJ)
      CASE(5)
        CALL DSHAPE_1D_5ND(NY,Y,SJ)
      END SELECT
      ND=0
      DO J=1,NY
        DO I=1,NX
          ND=ND+1
          DN(2,ND)=SI(I)*SJ(J)
        ENDDO
      ENDDO
      END
      
      SUBROUTINE DSHAPE_3D(NX,NY,NZ,NODE,X,Y,Z,DN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DN(3,NODE),SI(NX),SJ(NY),SK(NZ)
      ! EVALUATE DN/DXI
      SELECT CASE(NX)
      CASE(2)
        CALL DSHAPE_1D_2ND(NX,X,SI)
      CASE(3)
        CALL DSHAPE_1D_3ND(NX,X,SI)
      CASE(4)
        CALL DSHAPE_1D_4ND(NX,X,SI)
      CASE(5)
        CALL DSHAPE_1D_5ND(NX,X,SI)
      END SELECT
      
      SELECT CASE(NY)
      CASE(2)
        CALL SHAPE_1D_2ND(NY,Y,SJ)
      CASE(3)
        CALL SHAPE_1D_3ND(NY,Y,SJ)
      CASE(4)
        CALL SHAPE_1D_4ND(NY,Y,SJ)
      CASE(5)
        CALL SHAPE_1D_5ND(NY,Y,SJ)
      END SELECT
      
      SELECT CASE(NZ)
      CASE(2)
        CALL SHAPE_1D_2ND(NZ,Z,SK)
      CASE(3)
        CALL SHAPE_1D_3ND(NZ,Z,SK)
      CASE(4)
        CALL SHAPE_1D_4ND(NZ,Z,SK)
      CASE(5)
        CALL SHAPE_1D_5ND(NZ,Z,SK)
      END SELECT
      ND=0
      DO K=1,NZ
        DO J=1,NY
          DO I=1,NX
           ND=ND+1
           DN(1,ND)=SI(I)*SJ(J)*SK(K)
          ENDDO
        ENDDO
      ENDDO      
      ! EVALUATE DN/DETA
      SELECT CASE(NX)
      CASE(2)
        CALL SHAPE_1D_2ND(NX,X,SI)
      CASE(3)
        CALL SHAPE_1D_3ND(NX,X,SI)
      CASE(4)
        CALL SHAPE_1D_4ND(NX,X,SI)
      CASE(5)
        CALL SHAPE_1D_5ND(NX,X,SI)
      END SELECT
      
      SELECT CASE(NY)
      CASE(2)
        CALL DSHAPE_1D_2ND(NY,Y,SJ)
      CASE(3)
        CALL DSHAPE_1D_3ND(NY,Y,SJ)
      CASE(4)
        CALL DSHAPE_1D_4ND(NY,Y,SJ)
      CASE(5)
        CALL DSHAPE_1D_5ND(NY,Y,SJ)
      END SELECT
      
      SELECT CASE(NZ)
      CASE(2)
        CALL SHAPE_1D_2ND(NZ,Z,SK)
      CASE(3)
        CALL SHAPE_1D_3ND(NZ,Z,SK)
      CASE(4)
        CALL SHAPE_1D_4ND(NZ,Z,SK)
      CASE(5)
        CALL SHAPE_1D_5ND(NZ,Z,SK)
      END SELECT
      ND=0
      DO K=1,NZ
        DO J=1,NY
          DO I=1,NX
            ND=ND+1
            DN(2,ND)=SI(I)*SJ(J)*SK(K)
          ENDDO
        ENDDO
      ENDDO      
      ! EVALUATE DN/DZTA
      SELECT CASE(NX)
      CASE(2)
        CALL SHAPE_1D_2ND(NX,X,SI)
      CASE(3)
        CALL SHAPE_1D_3ND(NX,X,SI)
      CASE(4)
        CALL SHAPE_1D_4ND(NX,X,SI)
      CASE(5)
        CALL SHAPE_1D_5ND(NX,X,SI)
      END SELECT
      
      SELECT CASE(NY)
      CASE(2)
        CALL SHAPE_1D_2ND(NY,Y,SJ)
      CASE(3)
        CALL SHAPE_1D_3ND(NY,Y,SJ)
      CASE(4)
        CALL SHAPE_1D_4ND(NY,Y,SJ)
      CASE(5)
        CALL SHAPE_1D_5ND(NY,Y,SJ)
      END SELECT
      
      SELECT CASE(NZ)
      CASE(2)
        CALL DSHAPE_1D_2ND(NZ,Z,SK)
      CASE(3)
        CALL DSHAPE_1D_3ND(NZ,Z,SK)
      CASE(4)
        CALL DSHAPE_1D_4ND(NZ,Z,SK)
      CASE(5)
        CALL DSHAPE_1D_5ND(NZ,Z,SK)
      END SELECT
      ND=0
      DO K=1,NZ
        DO J=1,NY
          DO I=1,NX
            ND=ND+1
            DN(3,ND)=SI(I)*SJ(J)*SK(K)
          ENDDO
        ENDDO
      ENDDO    
      RETURN
      END
      
      SUBROUTINE DN_DXI_DXI(IETYPE,NODE,X,DDN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(3),DDN(3,3,NODE),DDSP1D(5)
      SELECT CASE(IETYPE)
      ! LINE ELEMENTS
      CASE(3) ! 3 NODE LINE ELEMENT
       CALL DDSHAPE_1D_3ND(3,X(1),DDSP1D)
       DDN(1,1,1:3)=DDSP1D(1:3)
      CASE(4) ! 4 NODE LINE ELEMENT
       CALL DDSHAPE_1D_4ND(4,X(1),DDSP1D)
       DDN(1,1,1:4)=DDSP1D(1:4)
      CASE(5) ! 5 NODE LINE ELEMENT
       CALL DDSHAPE_1D_5ND(5,X(1),DDSP1D)
       DDN(1,1,1:5)=DDSP1D(1:5)
       
      ! PLATE ELEMENTS
      CASE(307) ! 7 NODE TRIANGULAR ELEMENT
       CALL DDSHAPE_2D_7ND_TRI(NODE,X(1),X(2),DDN)
      CASE(308) ! 8 NODED RECTAGULAR ELEMENT
       CALL DDSHAPE_2D_8ND(NODE,X(1),X(2),DDN)   
      CASE(309) ! 9 NODED RECTAGULAR ELEMENT
       CALL DDSHAPE_2D(3,3,NODE,X(1),X(2),DDN)
      CASE(316) ! 16 NODED RECTAGULAR ELEMENT
       CALL DDSHAPE_2D(4,4,NODE,X(1),X(2),DDN)
      CASE(325) ! 25 NODED RECTAGULAR ELEMENT
       CALL DDSHAPE_2D(5,5,NODE,X(1),X(2),DDN)

      ! BRICK ELEMENTS
      CASE(411) ! 11 NODED  TETRAHEDRAL ELEMENT
       CALL DDSHAPE_3D_11ND_TRI(NODE,X(1),X(2),X(3),DDN)
      !CASE(420) ! 20 NODED BRICK ELEMENT
      ! CALL DDSHAPE_3D_20ND(NODE,X(1),X(2),X(3),DDN)
      CASE(421) ! 21 NODED BRICK ELEMENT
       CALL DDSHAPE_3D_21ND(NODE,X(1),X(2),X(3),DDN)
      CASE(427) ! 27 NODED BRICK ELEMENT
       CALL DDSHAPE_3D(3,3,3,NODE,X(1),X(2),X(3),DDN)
      CASE(464) ! 27 NODED BRICK ELEMENT
       CALL DDSHAPE_3D(4,4,4,NODE,X(1),X(2),X(3),DDN)
      CASE(525) ! 27 NODED BRICK ELEMENT
       CALL DDSHAPE_3D(5,5,5,NODE,X(1),X(2),X(3),DDN)
      END SELECT
      RETURN
      
      END
      
      SUBROUTINE DDSHAPE_2D(NX,NY,NODE,X,Y,DDN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DDN(3,3,NODE),SI(NX),SJ(NY),
     &          DSI(NX),DSJ(NY),DDSI(NX),DDSJ(NY)
      ! EVALUATE DN/DXI
      IF(NX.EQ.3) THEN
        CALL SHAPE_1D_3ND(NX,X,SI)
        CALL DSHAPE_1D_3ND(NX,X,DSI)
        CALL DDSHAPE_1D_3ND(NX,X,DDSI)
      ELSEIF(NX.EQ.4) THEN
        CALL SHAPE_1D_4ND(NX,X,SI)
        CALL DSHAPE_1D_4ND(NX,X,DSI)
        CALL DDSHAPE_1D_4ND(NX,X,DDSI)
      ELSEIF(NX.EQ.5) THEN
        CALL SHAPE_1D_5ND(NX,X,SI)
        CALL DSHAPE_1D_5ND(NX,X,DSI)
        CALL DDSHAPE_1D_5ND(NX,X,DDSI)
      ENDIF
      
      IF(NY.EQ.3) THEN
        CALL SHAPE_1D_3ND(NY,Y,SJ)
        CALL DSHAPE_1D_3ND(NY,Y,DSJ)
        CALL DDSHAPE_1D_3ND(NY,Y,DDSJ)
      ELSEIF(NY.EQ.4) THEN
        CALL SHAPE_1D_4ND(NY,Y,SJ)
        CALL DSHAPE_1D_4ND(NY,Y,DSJ)
        CALL DDSHAPE_1D_4ND(NY,Y,DDSJ)
      ELSEIF(NY.EQ.5) THEN
        CALL SHAPE_1D_5ND(NY,Y,SJ)
        CALL DSHAPE_1D_5ND(NY,Y,DSJ)
        CALL DDSHAPE_1D_5ND(NY,Y,DDSJ)
      ENDIF
      ! EVALUATE SECOND DERIVATIVES      
      ND=0
      DO J=1,NY
        DO I=1,NX
          ND=ND+1
          DDN(1,1,ND)=DDSI(I)*SJ(J)
          DDN(1,2,ND)=DSI(I)*DSJ(J)
          DDN(2,2,ND)=SI(I)*DDSJ(J)
          DDN(2,1,ND)=DDN(1,2,ND)
        ENDDO
      ENDDO
      RETURN
      END

      SUBROUTINE DDSHAPE_3D(NX,NY,NZ,NODE,X,Y,Z,DDN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DDN(3,3,NODE),SI(NX),SJ(NY),SK(NZ),
     &          DSI(NX),DSJ(NY),DSK(NZ),DDSI(NX),DDSJ(NY),DDSK(NZ)
      ! EVALUATE DN/DXI
      IF(NX.EQ.3) THEN
       CALL SHAPE_1D_3ND(NX,X,SI)
       CALL DSHAPE_1D_3ND(NX,X,DSI)
       CALL DDSHAPE_1D_3ND(NX,X,DDSI)
      ELSEIF(NX.EQ.4) THEN
       CALL SHAPE_1D_4ND(NX,X,SI)
       CALL DSHAPE_1D_4ND(NX,X,DSI)
       CALL DDSHAPE_1D_4ND(NX,X,DDSI)
      ELSEIF(NX.EQ.5) THEN
       CALL SHAPE_1D_5ND(NX,X,SI)
       CALL DSHAPE_1D_5ND(NX,X,DSI)
       CALL DDSHAPE_1D_5ND(NX,X,DDSI)
      ENDIF
      
      IF(NY.EQ.3) THEN
       CALL SHAPE_1D_3ND(NY,Y,SJ)
       CALL DSHAPE_1D_3ND(NY,Y,DSJ)
       CALL DDSHAPE_1D_3ND(NY,Y,DDSJ)
      ELSEIF(NY.EQ.4) THEN
       CALL SHAPE_1D_4ND(NY,Y,SJ)
       CALL DSHAPE_1D_4ND(NY,Y,DSJ)
       CALL DDSHAPE_1D_4ND(NY,Y,DDSJ)
      ELSEIF(NY.EQ.5) THEN
       CALL SHAPE_1D_5ND(NY,Y,SJ)
       CALL DSHAPE_1D_5ND(NY,Y,DSJ)
       CALL DDSHAPE_1D_5ND(NY,Y,DDSJ)
      ENDIF
      
      IF(NZ.EQ.3) THEN
       CALL SHAPE_1D_3ND(NZ,Z,SK)
       CALL DSHAPE_1D_3ND(NZ,Z,DSK)
       CALL DDSHAPE_1D_3ND(NZ,Z,DDSK)
      ELSEIF(NZ.EQ.4) THEN
       CALL SHAPE_1D_4ND(NZ,Z,SK)
       CALL DSHAPE_1D_4ND(NZ,Z,DSK)
       CALL DDSHAPE_1D_4ND(NZ,Z,DDSK)
      ELSEIF(NZ.EQ.5) THEN
       CALL SHAPE_1D_5ND(NZ,Z,SK)
       CALL DSHAPE_1D_5ND(NZ,Z,DSK)
       CALL DDSHAPE_1D_5ND(NZ,Z,DDSK)
      ENDIF
      ! EVALUATE SECOND DERIVATIVES      
      ND=0
      DO K=1,NZ
       DO J=1,NY
        DO I=1,NX
         ND=ND+1
         DDN(1,1,ND)=DDSI(I)*SJ(J)*SK(K)
         DDN(1,2,ND)=DSI(I)*DSJ(J)*SK(K)
         DDN(1,3,ND)=DSI(I)*SJ(J)*DSK(K)
         DDN(2,2,ND)=SI(I)*DDSJ(J)*SK(K)
         DDN(2,3,ND)=SI(I)*DSJ(J)*DSK(K)
         DDN(3,3,ND)=SI(I)*SJ(J)*DDSK(K)
         DDN(2,1,ND)=DDN(1,2,ND)
         DDN(3,1,ND)=DDN(1,3,ND)
         DDN(3,2,ND)=DDN(2,3,ND)
        ENDDO
       ENDDO
      ENDDO                    
      RETURN
      END

      SUBROUTINE SHAPE_1D_2ND(NODE,X,SHAPEF) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SHAPEF(NODE)
      ! FORM FROM NX=2
      ! NODAL VALUES OF KSI = {-1, 1}
      SHAPEF(1)=0.5D0*(1.D0-X)
      SHAPEF(2)=0.5D0*(1.D0+X)
      END

      SUBROUTINE DSHAPE_1D_2ND(NODE,X,DSHAPE) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DSHAPE(NODE)
      ! FORM FROM NX=2
      ! NODAL VALUES OF KSI = {-1, 1}
      DSHAPE(1)=-0.5D0
      DSHAPE(2)=0.5D0
      END

      SUBROUTINE SHAPE_1D_3ND(NODE,X,SHAPEF) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SHAPEF(NODE)
      ! FORM FROM NX=3
      ! NODAL VALUES OF KSI = {-1, 0, 1}
      SHAPEF(1)=-0.5D0*X*(1.D0-X)
      SHAPEF(2)=1.D0-X*X 
      SHAPEF(3)=0.5D0*X*(1.D0+X)
      END

      SUBROUTINE DSHAPE_1D_3ND(NODE,X,DSHAPE) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DSHAPE(NODE)
      ! FORM FROM NX=3
      ! NODAL VALUES OF KSI = {-1, 0, 1}
      DSHAPE(1)=-0.5D0*(1.D0-2.D0*X)
      DSHAPE(2)=-2.D0*X
      DSHAPE(3)=0.5D0*(1.D0+2.D0*X)
      END

      SUBROUTINE DDSHAPE_1D_3ND(NODE,X,D2N) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D2N(NODE)
      ! FORM FROM NX=3
      ! NODAL VALUES OF KSI = {-1, 0, 1}
      D2N(1)=1.D0
      D2N(2)=-2.D0
      D2N(3)=1.D0
      END

      SUBROUTINE SHAPE_1D_4ND(NODE,X,SHAPEF) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SHAPEF(NODE)
      ! FORM FROM NX=4
      ! NODAL VALUES OF KSI = {-1, -1/3, 1/3, 1}
      SHAPEF(1)=-(1.D0-X)*(10.D0-9.D0*(X*X+1.D0))/16.D0
      SHAPEF(2)=9.D0*(1.D0-X*X)*(1.D0-3.D0*X)/16.D0
      SHAPEF(3)=9.D0*(1.D0-X*X)*(1.D0+3.D0*X)/16.D0
      SHAPEF(4)=-(1.D0+X)*(10.D0-9.D0*(X*X+1.D0))/16.D0
      END

      SUBROUTINE DSHAPE_1D_4ND(NODE,X,DSHAPE) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DSHAPE(NODE)
      ! FORM FROM NX=4
      ! NODAL VALUES OF KSI = {-1, -1/3, 1/3, 1}
      DSHAPE(1)=(1.D0+9.D0*X*(2.D0-3.D0*X))/16.D0
      DSHAPE(2)=-9.D0*(3.D0+X*(2.D0-9.D0*X))/16.D0
      DSHAPE(3)=9.D0*(3.D0-X*(2.D0+9.D0*X))/16.D0
      DSHAPE(4)=-(1.D0-9.D0*X*(2.D0+3.D0*X))/16.D0
      END

      SUBROUTINE DDSHAPE_1D_4ND(NODE,X,DSHAPE) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DSHAPE(NODE)
      ! FORM FROM NX=4
      ! NODAL VALUES OF KSI = {-1, -1/3, 1/3, 1}
      DSHAPE(1)=9.D0*(1.D0-3.D0*X)/8.D0
      DSHAPE(2)=-9.D0*(1.D0-9.D0*X)/8.D0
      DSHAPE(3)=-9.D0*(1.D0+9.D0*X)/8.D0
      DSHAPE(4)= 9.D0*(1.D0+3.D0*X)/8.D0
      END

      SUBROUTINE SHAPE_1D_5ND(NODE,X,SHAPEF) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SHAPEF(NODE)
      ! FORM FROM NX=5
      ! NODAL VALUES OF KSI = {-1, -1/2, 0, 1/2, 1}
      SHAPEF(1)= X*(1.D0-X)*(1.D0-4.D0*X*X)/6.D0
      SHAPEF(2)= -4.D0*X*(1.D0-X*X)*(1.D0-2.D0*X)/3.D0
      SHAPEF(3)= (1.D0-X*X)*(1.D0-4.D0*X*X)
      SHAPEF(4)= 4.D0*X*(1.D0-X*X)*(1.D0+2.D0*X)/3.D0
      SHAPEF(5)= -X*(1.D0+X)*(1.D0-4.D0*X*X)/6.D0
      END

      SUBROUTINE DSHAPE_1D_5ND(NODE,X,DSHAPE) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DSHAPE(NODE)
      ! FORM FROM NX=5
      ! NODAL VALUES OF KSI = {-1, -1/2, 0, 1/2, 1}
      DSHAPE(1)= (1.D0-4.D0*X)*(1.D0+2.D0*X-4.D0*X*X)/6.D0
      DSHAPE(2)= 4.D0*(-1.D0+X*(4.D0+(3.D0-8.D0*X)*X))/3.D0
      DSHAPE(3)=-2.D0*X*(5.D0-8.D0*X*X)
      DSHAPE(4)= 4.D0*(1.D0+X*(4.D0-X*(3.D0+8.D0*X)))/3.D0
      DSHAPE(5)= (1.D0+4.D0*X)*(-1.D0+2*X+4.D0*X*X)/6.D0
      END

      SUBROUTINE DDSHAPE_1D_5ND(NODE,X,DSHAPE) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DSHAPE(NODE)
      ! FORM FROM NX=5
      ! NODAL VALUES OF KSI = {-1, -1/2, 0, 1/2, 1}
      DSHAPE(1)=-1.D0/3.D0-4.D0*X+8.D0*X*X 
      DSHAPE(2)=16.D0/3.D0+8.D0*(1.D0-4.D0*X)*X 
      DSHAPE(3)=-10.D0+48.D0*X*X
      DSHAPE(4)= 16.D0/3.D0-8.D0*X*(1.D0+4.D0*X) 
      DSHAPE(5)=-1.D0/3.D0+4*X+8.D0*X*X 
      END

      SUBROUTINE DSHAPE_2D_6ND_TRI(NODE,X,Y,DSP)
      IMPLICIT NONE
      INTEGER::NODE
      REAL*8::X,Y
      REAL*8,DIMENSION(3,NODE)::DSP
      
      DSP(1,1)=X-0.1D1/0.6D1+DSQRT(0.3D1)*Y/0.3D1
      DSP(1,2)=0.1D1/0.6D1+X-DSQRT(0.3D1)*Y/0.3D1
      DSP(1,3)=0.D0
      DSP(1,4)=-0.2D1*X
      DSP(1,5)=0.2D1/0.3D1+0.2D1/0.3D1*DSQRT(0.3D1)*Y
      DSP(1,6)=-0.2D1/0.3D1-0.2D1/0.3D1*DSQRT(0.3D1)*Y
      
      DSP(2,1)=-DSQRT(0.3D1)/0.18D2+Y/0.3D1+X*DSQRT(0.3D1)/0.3D1
      DSP(2,2)=-DSQRT(0.3D1)*(0.1D1+0.6D1*X)/0.18D2+Y/0.3D1
      DSP(2,3)=DSQRT(0.3D1)/0.9D1+0.4D1/0.3D1*Y
      DSP(2,4)=-0.4D1/0.9D1*DSQRT(0.3D1)+0.2D1/0.3D1*Y
      DSP(2,5)=0.2D1/0.9D1*DSQRT(0.3D1)*(0.1D1+0.3D1*X)-0.4D1/0.3D1*Y
      DSP(2,6)=-0.2D1/0.9D1*DSQRT(0.3D1)*(-0.1D1+0.3D1*X)-0.4D1/0.3D1*Y

      RETURN
      END
      
      SUBROUTINE SHAPE_2D_7ND_TRI(NODE,X,Y,SP)
      IMPLICIT NONE
      INTEGER::NDIM,NODE
      REAL*8::X,Y
      REAL*8,DIMENSION(NODE)::SP
      SP(1) = X * (-0.1D1 + 0.2D1 * DSQRT(0.3D1) * Y + X * (0.12D2 * Y *
     &* 2D0 - 0.1D1)) / 0.6D1
      SP(2) = X * (0.1D1 - 0.2D1 * DSQRT(0.3D1) * Y + X * (0.12D2 * Y **
     & 2D0 - 0.1D1)) / 0.6D1
      SP(3) = Y * (DSQRT(0.3D1) + 0.3D1 * Y) / 0.6D1 + 0.2D1 / 0.3D1 * X
     & ** 2 * (0.3D1 * Y ** 2 - 0.1D1)
      SP(4) = -0.2D1 / 0.3D1 * DSQRT(0.3D1) * Y + Y ** 2D0 + X ** 2D0 * 
     &(0.5D1 / 0.3D1 - 0.8D1 * Y ** 2D0)
      SP(5) = 0.2D1 / 0.3D1 * X * (0.1D1 + DSQRT(0.3D1) * Y + X * (-0.12
     &D2 * Y ** 2D0 + 0.4D1))
      SP(6)=-2.D0/3.D0*X*(1.D0+ DSQRT(3.D0)*Y-X*(-12.D0*Y**2.D0+4.D0))
      SP(7) = 0.1D1 + (DSQRT(0.3D1) - 0.3D1 * Y) * Y / 0.2D1 + 0.6D1 * X
     & ** 2D0 * (0.3D1 * Y ** 2D0 - 0.1D1)
      END



      SUBROUTINE DSHAPE_2D_7ND_TRI(NODE,X,Y,DSP)
      IMPLICIT NONE
      INTEGER::NDIM,NODE
      REAL*8::X,Y
      REAL*8,DIMENSION(3,NODE)::DSP
    
      DSP=0.D0
      DSP(1,1) = -0.1D1 / 0.6D1 + DSQRT(0.3D1) * Y / 0.3D1 + X * (0.12D2
     &* Y ** 2D0 - 0.1D1) / 0.3D1
      DSP(1,2) = 0.1D1 / 0.6D1 - DSQRT(0.3D1) * Y / 0.3D1 + X * (0.12D2*
     & Y ** 2D0 - 0.1D1) / 0.3D1
      DSP(1,3) = 0.4D1 / 0.3D1 * X * (0.3D1 * Y ** 2D0 - 0.1D1)
      DSP(1,4) = 0.2D1 * X * (0.5D1 / 0.3D1 - 0.8D1 * Y ** 2D0)
      DSP(1,5) = 0.2D1 / 0.3D1 + 0.2D1 / 0.3D1 * DSQRT(0.3D1) * Y +0.4D1
     & / 0.3D1 * X * (-0.12D2 * Y ** 2D0 + 0.4D1)
      DSP(1,6) = -0.2D1 / 0.3D1 - 0.2D1 / 0.3D1 * DSQRT(0.3D1) * Y+0.2D1
     &/ 0.3D1 * X * (-0.12D2 * Y ** 2D0 + 0.4D1) - 0.2D1 / 0.3D1 * X *
     & (0.12D2 * Y ** 2D0 - 0.4D1)
      DSP(1,7) = 0.12D2 * X * (0.3D1 * Y ** 2D0 - 0.1D1)
      
      DSP(2,1) = X * (0.2D1 * DSQRT(0.3D1) + 0.24D2 * X * Y) / 0.6D1
      DSP(2,2) = X * (-0.2D1 * DSQRT(0.3D1) + 0.24D2 * X * Y) / 0.6D1
      DSP(2,3) = DSQRT(0.3D1) / 0.6D1 + Y + 0.4D1 * X ** 2D0 * Y
      DSP(2,4) = -0.2D1/0.3D1*DSQRT(0.3D1)+0.2D1*Y-0.16D2*X**2D0*Y
      DSP(2,5) = 0.2D1 / 0.3D1 * X * (DSQRT(0.3D1) - 0.24D2 * X * Y)
      DSP(2,6) = -0.2D1 / 0.3D1 * X * (DSQRT(0.3D1) + 0.24D2 * X * Y)
      DSP(2,7) = -0.3D1 * Y + DSQRT(0.3D1) / 0.2D1 + 0.36D2 * X ** 2D0*Y
      END


      SUBROUTINE DDSHAPE_2D_7ND_TRI(NODE,X,Y,DDSP)
      IMPLICIT NONE
      INTEGER::NDIM,NODE
      REAL*8::X,Y
      REAL*8,DIMENSION(3,3,NODE)::DDSP
      DDSP=0.D0
      DDSP(1,1,1) = 4.D0 * Y ** 2.D0 - 0.1D1 / 0.3D1
      DDSP(1,1,2) = 4.D0 * Y ** 2.D0 - 0.1D1 / 0.3D1
      DDSP(1,1,3) = 4.D0 * Y ** 2.D0 - 0.4D1 / 0.3D1
      DDSP(1,1,4) = 0.1D2 / 0.3D1 - 16.D0 * Y ** 2.D0
      DDSP(1,1,5) = -16D0 * Y ** 2D0 + 0.16D2 / 0.3D1
      DDSP(1,1,6) = -16D0 * Y ** 2D0 + 0.16D2 / 0.3D1
      DDSP(1,1,7) = 36D0 * Y ** 2D0 - 12D0
      
      DDSP(1,2,1) = DSQRT(0.3D1) / 0.3D1 + 8.D0 * X * Y
      DDSP(1,2,2) = -DSQRT(0.3D1) / 0.3D1 + 8.D0 * X * Y
      DDSP(1,2,3) = 8D0 * X * Y
      DDSP(1,2,4) = -32D0 * X * Y
      DDSP(1,2,5) = 0.2D1 / 0.3D1 * DSQRT(0.3D1) - 32.D0 * X * Y
      DDSP(1,2,6) = -0.2D1 / 0.3D1 * DSQRT(0.3D1) - 32.D0 * X * Y
      DDSP(1,2,7) = 72.D0 * X * Y
      
      DDSP(2,1,1) = DDSP(1,2,1)
      DDSP(2,1,2) = DDSP(1,2,2)
      DDSP(2,1,3) = DDSP(1,2,3)
      DDSP(2,1,4) = DDSP(1,2,4)
      DDSP(2,1,5) = DDSP(1,2,5)
      DDSP(2,1,6) = DDSP(1,2,6)
	DDSP(2,1,7) = DDSP(1,2,7)
      
      DDSP(2,2,1) = 4.D0 * X ** 2D0
      DDSP(2,2,2) = 4D0 * X ** 2D0
      DDSP(2,2,3) = 4D0 * X ** 2D0+ 1D0
      DDSP(2,2,4) = -16D0 * X ** 2D0 + 2D0
      DDSP(2,2,5) = -16D0 * X ** 2D0
      DDSP(2,2,6) = -16D0 * X ** 2D0
      DDSP(2,2,7) = 36D0 * X ** 2D0 - 3D0

      END
      
      SUBROUTINE SHAPE_3D_11ND_TRI(NODE,X,Y,Z,SP)
      IMPLICIT NONE
      INTEGER::NODE
      REAL*8::X,Y,Z
      REAL*8::SP(NODE)
      
      SP(1)=(-0.1D1+DSQRT(0.2D1))*((0.2D1-0.4D1*DSQRT(0.2D1))*
     &Z**2.D0+Y*Z*(0.8D1+DSQRT(0.2D1)-0.12D2*DSQRT(0.3D1)*Z)
     &+0.6D1*X**2.D0*(DSQRT(0.2D1)-0.2D1*DSQRT(0.3D1)*Y-0.2D1
     &*DSQRT(0.3D1)*Z)-0.2D1*Y**2.D0*(0.1D1+DSQRT(0.2D1)+0.6D1
     &*DSQRT(0.3D1)*Z)+DSQRT(0.3D1)*(0.3D1*Y-(-0.4D1+DSQRT(0
     &.2D1))*Z)+X*(0.3D1*DSQRT(0.3D1)+0.2D1*Y*(-0.3D1+0.4D1
     &*DSQRT(0.3D1)+0.4D1*DSQRT(0.6D1)-0.6D1*DSQRT(0.3D1)*Y)
     &+Z*(-0.3D1*DSQRT(0.2D1)+0.8D1*DSQRT(0.3D1)+0.4D1*DSQRT(0
     &.6D1)-0.12D2*DSQRT(0.3D1)*Z)))/0.24D2
      SP(2)=-(-0.1D1+DSQRT(0.2D1))*((-0.2D1+0.4D1*DSQRT(0.2D1))
     &*Z**2.D0-0.6D1*X**2.D0*(DSQRT(0.2D1)-0.2D1*DSQRT(0.3D1)*
     &Y-0.2D1*DSQRT(0.3D1)*Z)+0.2D1*Y**2.D0*(0.1D1+DSQRT(0.2D1)
     &+0.6D1*DSQRT(0.3D1)*Z)+Y*Z*(-0.8D1-DSQRT(0.2D1)+
     &0.12D2*DSQRT(0.3D1)*Z)+DSQRT(0.3D1)*(-0.3D1*Y+(-0.4D1+
     &DSQRT(0.2D1))*Z)+X*(-0.3D1*DSQRT(0.3D1)+0.2D1*Y*(0.3D1+
     &0.4D1*DSQRT(0.3D1)+0.4D1*DSQRT(0.6D1)+0.6D1*DSQRT(0.3D1)*
     &Y)+Z*(0.3D1*DSQRT(0.2D1)+0.8D1*DSQRT(0.3D1)+0.4D1*
     &DSQRT(0.6D1)+0.12D2*DSQRT(0.3D1)*Z)))/0.24D2
      SP(3)=-(-0.1D1+DSQRT(0.2D1))*(-0.10D2*(0.1D1+DSQRT(0.2D1)
     &)*Y**2.D0+Y*(0.16D2+0.11D2*DSQRT(0.2D1)+0.12D2*DSQRT(0.
     &3D1)*Y)*Z+0.2D1*(-0.1D1+0.2D1*DSQRT(0.2D1)+0.6D1*
     &DSQRT(0.3D1)*Y)*Z**2.D0+0.6D1*X**2.D0*(0.2D1+DSQRT(0.2D1)+
     &0.2D1*DSQRT(0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)+DSQRT(0.3D1)*
     &(-0.3D1*Y+(-0.4D1+DSQRT(0.2D1))*Z)+0.3D1*X*(-DSQRT(0.
     &3D1)+0.2D1*Y+0.4D1*DSQRT(0.3D1)*Y**2.D0+Z*(DSQRT(0.2D1)
     &+0.4D1*DSQRT(0.3D1)*Z)))/0.24D2
      SP(4)=-(-0.1D1+DSQRT(0.2D1))*(DSQRT(0.3D1)*(-0.3D1*Y+
     &(-0.4D1+DSQRT(0.2D1))*Z)-0.6D1*(0.3D1+0.2D1*DSQRT(0.2D1)) 
     &*Z**2.D0+0.6D1*Y**2.D0*(0.1D1+DSQRT(0.2D1)+0.2D1*DSQRT(
     &0.3D1)*Z)+0.6D1*X**2.D0*(0.2D1+DSQRT(0.2D1)+0.2D1*DSQRT(
     &0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)+0.3D1*Y*Z*(DSQRT(0.2D1
     &)+0.4D1*DSQRT(0.3D1)*Z)+0.3D1*X*(-DSQRT(0.3D1)+0.2D1 *
     &Y+0.4D1*DSQRT(0.3D1)*Y**2.D0+Z*(DSQRT(0.2D1)+0.4D1*
     &DSQRT(0.3D1)*Z)))/0.24D2
      SP(5)=(-0.1D1+DSQRT(0.2D1))*((0.2D1+0.8D1*DSQRT(0.2D1)) *
     &Z**2.D0-0.6D1*X**2.D0*(DSQRT(0.2D1)-0.2D1*DSQRT(0.3D1)*Y 
     &-0.2D1*DSQRT(0.3D1)*Z)+0.2D1*Y**2.D0*(0.5D1+0.5D1*
     &DSQRT(0.2D1)+0.6D1*DSQRT(0.3D1)*Z)+Y*Z*(0.8D1+0.7D1*
     &DSQRT(0.2D1)+0.12D2*DSQRT(0.3D1)*Z)-DSQRT(0.3D1)*((0.7D1+
     &0.4D1*DSQRT(0.2D1))*Y+(0.8D1+DSQRT(0.2D1))*Z)+0.3D1*X*
     &(-DSQRT(0.3D1)+0.2D1*Y+0.4D1*DSQRT(0.3D1)*Y**2.D0+Z*
     &(DSQRT(0.2D1)+0.4D1*DSQRT(0.3D1)*Z)))/0.12D2
      SP(6)=(-0.1D1+DSQRT(0.2D1))*((0.2D1+0.8D1*DSQRT(0.2D1)) *
     &Z**2.D0-0.2D1*Y**2.D0*(0.1D1+DSQRT(0.2D1)-0.6D1*DSQRT(
     &0.3D1)*Z)+0.6D1*X**2.D0*(0.2D1+DSQRT(0.2D1)+0.2D1*DSQRT(
     &0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)+Y*Z*(-0.4D1+DSQRT(0.2D1
     &)+0.12D2*DSQRT(0.3D1)*Z)+DSQRT(0.3D1)*((-0.1D1+0.2D1*
     &DSQRT(0.2D1))*Y-(0.8D1+DSQRT(0.2D1))*Z)+X*(0.6D1+0.6D1 
     &*DSQRT(0.2D1)-0.3D1*DSQRT(0.3D1)+0.2D1*Y*(0.3D1+0.4D1 *
     &DSQRT(0.3D1)+0.4D1*DSQRT(0.6D1)+0.6D1*DSQRT(0.3D1)*Y)+Z 
     &*(0.3D1*DSQRT(0.2D1)-0.2D1*DSQRT(0.3D1)*(0.2D1+DSQRT(0.2D1
     &))+0.12D2*DSQRT(0.3D1)*Z)))/0.12D2
      SP(7)=(-0.1D1+DSQRT(0.2D1))*((0.2D1+0.8D1*DSQRT(0.2D1)) *
     & Z**2.D0-0.2D1*Y**2.D0*(0.1D1+DSQRT(0.2D1)-0.6D1*DSQRT(
     &0.3D1)*Z)+0.6D1*X**2.D0*(0.2D1+DSQRT(0.2D1)+0.2D1*DSQRT(
     &0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)+Y*Z*(-0.4D1+DSQRT(0.2D1
     &)+0.12D2*DSQRT(0.3D1)*Z)+DSQRT(0.3D1)*((-0.1D1+0.2D1*
     &DSQRT(0.2D1))*Y-(0.8D1+DSQRT(0.2D1))*Z)+X*(-0.6D1-0.6D1
     &*DSQRT(0.2D1)-0.3D1*DSQRT(0.3D1)+0.2D1*Y*(0.3D1-0.4D1 
     &*DSQRT(0.3D1)-0.4D1*DSQRT(0.6D1)+0.6D1*DSQRT(0.3D1)*Y)+
     &DSQRT(0.3D1)*Z*(0.4D1+0.2D1*DSQRT(0.2D1)+DSQRT(0.6D1)+
     &0.12D2*Z)))/0.12D2
      SP(8)=(-0.1D1+DSQRT(0.2D1))*(0.3D1*(DSQRT(0.6D1)-0.2D1*
     &Z)*Z+0.6D1*Y**2.D0*(0.1D1+DSQRT(0.2D1)+0.2D1*DSQRT(
     &0.3D1)*Z)+0.6D1*X**2.D0*(0.2D1+DSQRT(0.2D1)+0.2D1*DSQRT(
     &0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)-Y*(0.5D1*DSQRT(0.3D1)+
     &0.2D1*DSQRT(0.6D1)+0.3D1*Z*(0.4D1+DSQRT(0.2D1)-0.4D1*
     &DSQRT(0.3D1)*Z))-0.3D1*X*(0.2D1+0.2D1*DSQRT(0.2D1)+DSQRT(
     &0.3D1)-0.2D1*Y*(0.1D1+0.2D1*DSQRT(0.3D1)*Y)-Z*(DSQRT(
     &0.2D1)-0.2D1*DSQRT(0.3D1)*(0.2D1+DSQRT(0.2D1))+0.4D1*
     &DSQRT(0.3D1)*Z)))/0.12D2
      SP(9)=(-0.1D1+DSQRT(0.2D1))*(0.3D1*(DSQRT(0.6D1)-0.2D1*
     &Z)*Z+0.6D1*Y**2.D0*(0.1D1+DSQRT(0.2D1)+0.2D1*DSQRT(
     &0.3D1)*Z)+0.6D1*X**2.D0*(0.2D1+DSQRT(0.2D1)+0.2D1*DSQRT(
     &0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)-Y*(0.5D1*DSQRT(0.3D1)+
     & 0.2D1*DSQRT(0.6D1)+0.3D1*Z*(0.4D1+DSQRT(0.2D1)-0.4D1*
     &DSQRT(0.3D1)*Z))+0.3D1*X*(0.2D1+0.2D1*DSQRT(0.2D1)-DSQRT(
     &0.3D1)+0.2D1*Y*(0.1D1+0.2D1*DSQRT(0.3D1)*Y)+Z*(DSQRT(
     &0.2D1)+0.2D1*DSQRT(0.3D1)*(0.2D1+DSQRT(0.2D1))+0.4D1*
     &DSQRT(0.3D1)*Z)))/0.12D2
      SP(10)=(-0.1D1+DSQRT(0.2D1))*(-0.3D1*DSQRT(0.3D1)*X+
     &0.3D1*(DSQRT(0.6D1)-0.2D1*Z)*Z+0.6D1*X**2.D0*(0.2D1+
     &DSQRT(0.2D1)+0.2D1*DSQRT(0.3D1)*Z)+0.6D1*Y**2.D0*(0.1D1+
     &DSQRT(0.2D1)+0.2D1*DSQRT(0.3D1)*X+0.2D1*DSQRT(0.3D1)*Z)+
     &0.3D1*X*Z*(DSQRT(0.2D1)+0.4D1*DSQRT(0.3D1)*Z)+Y*(
     &DSQRT(0.3D1)+0.4D1*DSQRT(0.6D1)+0.6D1*X*(0.1D1+0.2D1*
     &DSQRT(0.3D1)*X)+0.3D1*Z*(0.8D1+0.5D1*DSQRT(0.2D1)+0.4D1*
     &DSQRT(0.3D1)*Z)))/0.12D2
      SP(11)=-(-0.1D1+DSQRT(0.2D1))*(-0.3D1-0.3D1*DSQRT(0.2D1) 
     &-0.3D1*DSQRT(0.3D1)*Y+Z*(-0.4D1*DSQRT(0.3D1)+DSQRT(0.6D1
     &)+0.6D1*DSQRT(0.2D1)*Z)+0.6D1*Y**2.D0*(0.1D1+DSQRT(0.2D1
     &)+0.2D1*DSQRT(0.3D1)*Z)+0.6D1*X**2.D0*(0.2D1+DSQRT(0.2D1
     &)+0.2D1*DSQRT(0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)+0.3D1*
     &Y*Z*(DSQRT(0.2D1)+0.4D1*DSQRT(0.3D1)*Z)+0.3D1*X*(
     &-DSQRT(0.3D1)+0.2D1*Y+0.4D1*DSQRT(0.3D1)*Y**2.D0+Z*(DSQRT(
     &0.2D1)+0.4D1*DSQRT(0.3D1)*Z)))/0.3D1
      
      END
      
      SUBROUTINE DSHAPE_3D_11ND_TRI(NODE,X,Y,Z,DSP)  
      
      IMPLICIT NONE
      INTEGER::NODE
      REAL*8::X,Y,Z
      REAL*8::DSP(3,NODE)
      
      DSP(1,1)=(-0.1D1+DSQRT(0.2D1))*(0.12D2*X*(DSQRT(0.2D1)-0
     &.2D1*DSQRT(0.3D1)*Y-0.2D1*DSQRT(0.3D1)*Z)+0.3D1*DSQRT(0
     &.3D1)+0.2D1*Y*(-0.3D1+0.4D1*DSQRT(0.3D1)+0.4D1*DSQRT(0
     &.6D1)-0.6D1*DSQRT(0.3D1)*Y)+Z*(-0.3D1*DSQRT(0.2D1)+0.8
     &D1*DSQRT(0.3D1)+0.4D1*DSQRT(0.6D1)-0.12D2*DSQRT(0.3D1)*Z)
     &)/0.24D2
      DSP(1,2)=-(-0.1D1+DSQRT(0.2D1))*(-0.12D2*X*(DSQRT(0.2D1)-
     &0.2D1*DSQRT(0.3D1)*Y-0.2D1*DSQRT(0.3D1)*Z)-0.3D1*DSQRT
     &(0.3D1)+0.2D1*Y*(0.3D1+0.4D1*DSQRT(0.3D1)+0.4D1*DSQRT(
     &0.6D1)+0.6D1*DSQRT(0.3D1)*Y)+Z*(0.3D1*DSQRT(0.2D1)+0.8
     &D1*DSQRT(0.3D1)+0.4D1*DSQRT(0.6D1)+0.12D2*DSQRT(0.3D1)*Z)
     &)/0.24D2
      DSP(1,3)=-(-0.1D1+DSQRT(0.2D1))*(0.12D2*X*(0.2D1+DSQRT(0
     &.2D1)+0.2D1*DSQRT(0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)-0.3D1
     &*DSQRT(0.3D1)+0.6D1*Y+0.12D2*DSQRT(0.3D1)*Y**2.D0+0.3D1
     &*Z*(DSQRT(0.2D1)+0.4D1*DSQRT(0.3D1)*Z))/0.24D2
      DSP(1,4)=-(-0.1D1+DSQRT(0.2D1))*(0.12D2*X*(0.2D1+DSQRT(0
     &.2D1)+0.2D1*DSQRT(0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)-0.3D1
     &*DSQRT(0.3D1)+0.6D1*Y+0.12D2*DSQRT(0.3D1)*Y**2.D0+0.3D1
     &*Z*(DSQRT(0.2D1)+0.4D1*DSQRT(0.3D1)*Z))/0.24D2
      DSP(1,5)=(-0.1D1+DSQRT(0.2D1))*(-0.12D2*X*(DSQRT(0.2D1)-
     &0.2D1*DSQRT(0.3D1)*Y-0.2D1*DSQRT(0.3D1)*Z)-0.3D1*DSQRT(
     &0.3D1)+0.6D1*Y+0.12D2*DSQRT(0.3D1)*Y**2.D0+0.3D1*Z*(
     &DSQRT(0.2D1)+0.4D1*DSQRT(0.3D1)*Z))/0.12D2
      DSP(1,6)=(-0.1D1+DSQRT(0.2D1))*(0.12D2*X*(0.2D1+DSQRT(0.
     &2D1)+0.2D1*DSQRT(0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)+0.6D1
     &+0.6D1*DSQRT(0.2D1)-0.3D1*DSQRT(0.3D1)+0.2D1*Y*(0.3D1+
     &0.4D1*DSQRT(0.3D1)+0.4D1*DSQRT(0.6D1)+0.6D1*DSQRT(0.3D1)*
     &Y)+Z*(0.3D1*DSQRT(0.2D1)-0.2D1*DSQRT(0.3D1)*(0.2D1+DSQRT(
     &0.2D1))+0.12D2*DSQRT(0.3D1)*Z))/0.12D2
      DSP(1,7)=(-0.1D1+DSQRT(0.2D1))*(0.12D2*X*(0.2D1+DSQRT(0.
     &2D1)+0.2D1*DSQRT(0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)-0.6D1
     &-0.6D1*DSQRT(0.2D1)-0.3D1*DSQRT(0.3D1)+0.2D1*Y*(0.3D1-
     &0.4D1*DSQRT(0.3D1)-0.4D1*DSQRT(0.6D1)+0.6D1*DSQRT(0.3D1)*
     &Y)+DSQRT(0.3D1)*Z*(0.4D1+0.2D1*DSQRT(0.2D1)+DSQRT(0.6D1)
     &+0.12D2*Z))/0.12D2
      DSP(1,8)=(-0.1D1+DSQRT(0.2D1))*(0.12D2*X*(0.2D1+DSQRT(0.
     &2D1)+0.2D1*DSQRT(0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)-0.6D1
     &-0.6D1*DSQRT(0.2D1)-0.3D1*DSQRT(0.3D1)+0.6D1*Y*(0.1D1+
     &0.2D1*DSQRT(0.3D1)*Y)+0.3D1*Z*(DSQRT(0.2D1)-0.2D1*DSQRT(
     &0.3D1)*(0.2D1+DSQRT(0.2D1))+0.4D1*DSQRT(0.3D1)*Z))/0.1
     &2D2
      DSP(1,9)=(-0.1D1+DSQRT(0.2D1))*(0.12D2*X*(0.2D1+DSQRT(0.
     &2D1)+0.2D1*DSQRT(0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)+0.6D1
     &+ 0.6D1*DSQRT(0.2D1)-0.3D1*DSQRT(0.3D1)+0.6D1*Y*(0.1D1+
     &0.2D1*DSQRT(0.3D1)*Y)+0.3D1*Z*(DSQRT(0.2D1)+0.2D1*DSQRT(
     &0.3D1)*(0.2D1+DSQRT(0.2D1))+0.4D1*DSQRT(0.3D1)*Z))/0.1
     &2D2
      DSP(1,10)=(-0.1D1+DSQRT(0.2D1))*(-0.3D1*DSQRT(0.3D1)+0.12D
     &2*X*(0.2D1+DSQRT(0.2D1)+0.2D1*DSQRT(0.3D1)*Z)+0.12D2*
     &DSQRT(0.3D1)*Y**2.D0+0.3D1*Z*(DSQRT(0.2D1)+0.4D1*DSQRT(0.
     &3D1)*Z)+Y*(0.6D1+0.24D2*DSQRT(0.3D1)*X))/0.12D2
      DSP(1,11)=-(-0.1D1+DSQRT(0.2D1))*(0.12D2*X*(0.2D1+DSQRT(
     &0.2D1)+0.2D1*DSQRT(0.3D1)*Y+0.2D1*DSQRT(0.3D1)*Z)-0.3D
     &1*DSQRT(0.3D1)+0.6D1*Y+0.12D2*DSQRT(0.3D1)*Y**2.D0+0.3D
     &1*Z*(DSQRT(0.2D1)+0.4D1*DSQRT(0.3D1)*Z))/0.3D1
      DSP(2,1)=(-0.1D1+DSQRT(0.2D1))*(Z*(0.8D1+DSQRT(0.2D1)-0.
     &12D2*DSQRT(0.3D1)*Z)-0.12D2*X**2.D0*DSQRT(0.3D1)-0.4D1*
     &Y*(0.1D1+DSQRT(0.2D1)+0.6D1*DSQRT(0.3D1)*Z)+0.3D1*DSQRT(
     &0.3D1)+X*(-0.6D1+0.8D1*DSQRT(0.3D1)+0.8D1*DSQRT(0.6D1)
     &-0.24D2*DSQRT(0.3D1)*Y))/0.24D2
      DSP(2,2)=-(-0.1D1+DSQRT(0.2D1))*(0.12D2*X**2.D0*DSQRT(0.3D1
     &)+0.4D1*Y*(0.1D1+DSQRT(0.2D1)+0.6D1*DSQRT(0.3D1)*Z)+
     &Z*(-0.8D1-DSQRT(0.2D1)+0.12D2*DSQRT(0.3D1)*Z)-0.3D1*DSQRT(
     &0.3D1)+X*(0.6D1+0.8D1*DSQRT(0.3D1)+0.8D1*DSQRT(0.6D1)
     &+0.24D2*DSQRT(0.3D1)*Y))/0.24D2
      DSP(2,3)=-(-0.1D1+DSQRT(0.2D1))*(-0.20D2*(0.1D1+DSQRT(0.2D
     &1))*Y+(0.16D2+0.11D2*DSQRT(0.2D1)+0.12D2*DSQRT(0.3D1)*
     &Y)*Z+0.12D2*Y*Z*DSQRT(0.3D1)+0.12D2*DSQRT(0.3D1)*Z*
     &*2+0.12D2*X**2.D0*DSQRT(0.3D1)-0.3D1*DSQRT(0.3D1)+0.3D1
     &*X*(0.2D1+0.8D1*DSQRT(0.3D1)*Y))/0.24D2
      DSP(2,4)=-(-0.1D1+DSQRT(0.2D1))*(-0.3D1*DSQRT(0.3D1)+0.12D
     &2*Y*(0.1D1+DSQRT(0.2D1)+0.2D1*DSQRT(0.3D1)*Z)+0.12D2*
     &X**2.D0*DSQRT(0.3D1)+0.3D1*Z*(DSQRT(0.2D1)+0.4D1*DSQRT(0.
     &3D1)*Z)+0.3D1*X*(0.2D1+0.8D1*DSQRT(0.3D1)*Y))/0.24D2
      DSP(2,5)=(-0.1D1+DSQRT(0.2D1))*(0.12D2*X**2.D0*DSQRT(0.3D1)
     &+0.4D1*Y*(0.5D1+0.5D1*DSQRT(0.2D1)+0.6D1*DSQRT(0.3D1)
     &*Z)+Z*(0.8D1+0.7D1*DSQRT(0.2D1)+0.12D2*DSQRT(0.3D1)*Z
     &)-DSQRT(0.3D1)*(0.7D1+0.4D1*DSQRT(0.2D1))+0.3D1*X*(0.2
     &D1+0.8D1*DSQRT(0.3D1)*Y))/0.12D2
      DSP(2,6)=(-0.1D1+DSQRT(0.2D1))*(-0.4D1*Y*(0.1D1+DSQRT(0.
     &2D1)-0.6D1*DSQRT(0.3D1)*Z)+0.12D2*X**2.D0*DSQRT(0.3D1)+
     &Z*(-0.4D1+DSQRT(0.2D1)+0.12D2*DSQRT(0.3D1)*Z)+DSQRT(0.3D1
     &)*(-0.1D1+0.2D1*DSQRT(0.2D1))+X*(0.6D1+0.8D1*DSQRT(0.3
     &D1)+0.8D1*DSQRT(0.6D1)+0.24D2*DSQRT(0.3D1)*Y))/0.12D2
      DSP(2,7)=(-0.1D1+DSQRT(0.2D1))*(-0.4D1*Y*(0.1D1+DSQRT(0.
     &2D1)-0.6D1*DSQRT(0.3D1)*Z)+0.12D2*X**2.D0*DSQRT(0.3D1)+
     &Z*(-0.4D1+DSQRT(0.2D1)+0.12D2*DSQRT(0.3D1)*Z)+DSQRT(0.3D1
     &)*(-0.1D1+0.2D1*DSQRT(0.2D1))+X*(0.6D1-0.8D1*DSQRT(0.3
     &D1)-0.8D1*DSQRT(0.6D1)+0.24D2*DSQRT(0.3D1)*Y))/0.12D2
      DSP(2,8)=(-0.1D1+DSQRT(0.2D1))*(0.12D2*Y*(0.1D1+DSQRT(0.
     &2D1)+0.2D1*DSQRT(0.3D1)*Z)+0.12D2*X**2.D0*DSQRT(0.3D1)-
     &0.5D1*DSQRT(0.3D1)-0.2D1*DSQRT(0.6D1)-0.3D1*Z*(0.4D1+DSQRT(
     &0.2D1)-0.4D1*DSQRT(0.3D1)*Z)-0.3D1*X*(-0.2D1-0.8D
     &1*DSQRT(0.3D1)*Y))/0.12D2
      DSP(2,9)=(-0.1D1+DSQRT(0.2D1))*(0.12D2*Y*(0.1D1+DSQRT(0.
     &2D1)+0.2D1*DSQRT(0.3D1)*Z)+0.12D2*X**2.D0*DSQRT(0.3D1)-
     &0.5D1*DSQRT(0.3D1)-0.2D1*DSQRT(0.6D1)-0.3D1*Z*(0.4D1+DSQRT(0.2D1)-
     &0.4D1*DSQRT(0.3D1)*Z)+0.3D1*X*(0.2D1+0.8D1*DSQRT(0.3D1)*Y))/0.12D2
      DSP(2,10)=(-0.1D1+DSQRT(0.2D1))*(0.12D2*Y*(0.1D1+DSQRT(0
     &.2D1)+0.2D1*DSQRT(0.3D1)*X+0.2D1*DSQRT(0.3D1)*Z)+DSQRT(
     &0.3D1)+0.4D1*DSQRT(0.6D1)+0.6D1*X*(0.1D1+0.2D1*DSQRT(0.3D1)*X)
     &+0.3D1*Z*(0.8D1+0.5D1*DSQRT(0.2D1)+0.4D1*DSQRT(0.3D1)*Z))/0.12D2
      DSP(2,11)=-(-0.1D1+DSQRT(0.2D1))*(-0.3D1*DSQRT(0.3D1)+0.12
     &D2*Y*(0.1D1+DSQRT(0.2D1)+0.2D1*DSQRT(0.3D1)*Z)+0.12D2
     &*X**2.D0*DSQRT(0.3D1)+0.3D1*Z*(DSQRT(0.2D1)+0.4D1*DSQRT(0
     &.3D1)*Z)+0.3D1*X*(0.2D1+0.8D1*DSQRT(0.3D1)*Y))/0.3D1
      DSP(3,1)=(-0.1D1+DSQRT(0.2D1))*(0.2D1*(0.2D1-0.4D1*DSQRT(
     &0.2D1))*Z+Y*(0.8D1+DSQRT(0.2D1)-0.12D2*DSQRT(0.3D1)*Z
     &)-0.12D2*Y*Z*DSQRT(0.3D1)-0.12D2*X**2.D0*DSQRT(0.3D1)-
     &0.12D2*DSQRT(0.3D1)*Y**2.D0+DSQRT(0.3D1)*(0.4D1-DSQRT(0.2D1
     &))+X*(-0.3D1*DSQRT(0.2D1)+0.8D1*DSQRT(0.3D1)+0.4D1*DSQRT(
     &0.6D1)-0.24D2*DSQRT(0.3D1)*Z))/0.24D2
      DSP(3,2)=-(-0.1D1+DSQRT(0.2D1))*(0.2D1*(-0.2D1+0.4D1*DSQRT(
     &0.2D1))*Z+0.12D2*X**2.D0*DSQRT(0.3D1)+0.12D2*DSQRT(0.3
     &D1)*Y**2.D0+Y*(-0.8D1-DSQRT(0.2D1)+0.12D2*DSQRT(0.3D1)*
     &Z)+0.12D2*Y*Z*DSQRT(0.3D1)+DSQRT(0.3D1)*(-0.4D1+DSQRT(0
     &.2D1))+X*(0.3D1*DSQRT(0.2D1)+0.8D1*DSQRT(0.3D1)+0.4D1*
     &DSQRT(0.6D1)+0.24D2*DSQRT(0.3D1)*Z))/0.24D2
      DSP(3,3)=-(-0.1D1+DSQRT(0.2D1))*(Y*(0.16D2+0.11D2*DSQRT(
     &0.2D1)+0.12D2*DSQRT(0.3D1)*Y)+0.4D1*(-0.1D1+0.2D1*DSQRT(
     &0.2D1)+0.6D1*DSQRT(0.3D1)*Y)*Z+0.12D2*X**2.D0*DSQRT(0.3D1)
     &+DSQRT(0.3D1)*(-0.4D1+DSQRT(0.2D1))+0.3D1*X*(DSQRT(0.2D1)
     &+0.8D1*DSQRT(0.3D1)*Z))/0.24D2
      DSP(3,4)=-(-0.1D1+DSQRT(0.2D1))*(DSQRT(0.3D1)*(-0.4D1+DSQRT(
     &0.2D1))-0.12D2*(0.3D1+0.2D1*DSQRT(0.2D1))*Z+0.12D2*DSQRT(0.3D1)
     &*Y**2.D0+0.12D2*X**2.D0*DSQRT(0.3D1)+0.3D1*Y*(DSQRT(0.2D1)+
     &0.4D1*DSQRT(0.3D1)*Z)+0.12D2*Y*Z*DSQRT(0.3D1)+0.3D1*X*
     &(DSQRT(0.2D1)+0.8D1*DSQRT(0.3D1)*Z))/0.24D2
      DSP(3,5)=(-0.1D1+DSQRT(0.2D1))*(0.2D1*(0.2D1+0.8D1*DSQRT
     &(0.2D1))*Z+0.12D2*X**2.D0*DSQRT(0.3D1)+0.12D2*DSQRT(0.3D1
     &)*Y**2.D0+Y*(0.8D1+0.7D1*DSQRT(0.2D1)+0.12D2*DSQRT(0.3D
     &1)*Z)+0.12D2*Y*Z*DSQRT(0.3D1)-DSQRT(0.3D1)*(0.8D1+DSQRT(0.2D1))
     &+0.3D1*X*(DSQRT(0.2D1)+0.8D1*DSQRT(0.3D1)*Z))/0.12D2
      DSP(3,6)=(-0.1D1+DSQRT(0.2D1))*(0.2D1*(0.2D1+0.8D1*DSQRT(0.2D1))
     &*Z+0.12D2*DSQRT(0.3D1)*Y**2.D0+0.12D2*X**2.D0*DSQRT(0.3D1)+Y*(-0.4
     &D1+DSQRT(0.2D1)+0.12D2*DSQRT(0.3D1)*Z)+0.12D2*Y*Z*DSQRT(0.3D1)
     &+DSQRT(0.3D1)*(-0.8D1-DSQRT(0.2D1))+X*(0.3D1*DSQRT(0.2D1)-0.2D1*
     &DSQRT(0.3D1)*(0.2D1+DSQRT(0.2D1))+0.24D2*DSQRT(0.3D1)*Z))/0.12D2
      DSP(3,7)=(-0.1D1+DSQRT(0.2D1))*(0.2D1*(0.2D1+0.8D1*DSQRT(0.2D1))*Z
     &+0.12D2*DSQRT(0.3D1)*Y**2.D0+0.12D2*X**2.D0*DSQRT(0.3D1)+Y*(-0.4D1
     &+DSQRT(0.2D1)+0.12D2*DSQRT(0.3D1)*Z)+0.12D2*Y*Z*DSQRT(0.3D1)+
     &DSQRT(0.3D1)*(-0.8D1-DSQRT(0.2D1))+X*(DSQRT(0.3D1)*(0.4D1+0.2D1*
     &DSQRT(0.2D1)+DSQRT(0.6D1)+0.12D2*Z)+0.12D2*DSQRT(0.3D1)*Z))/0.12D2
      DSP(3,8)=(-0.1D1+DSQRT(0.2D1))*(-0.12D2*Z+0.3D1*DSQRT(0.
     &6D1)+0.12D2*DSQRT(0.3D1)*Y**2.D0+0.12D2*X**2.D0*DSQRT(0.3D
     &1)-Y*(0.12D2+0.3D1*DSQRT(0.2D1)-0.24D2*DSQRT(0.3D1)*Z)
     &-0.3D1*X*(-DSQRT(0.2D1)+0.2D1*DSQRT(0.3D1)*(0.2D1+DSQRT(0.2D1))
     &-0.8D1*DSQRT(0.3D1)*Z))/0.12D2
      DSP(3,9)=(-0.1D1+DSQRT(0.2D1))*(-0.12D2*Z+0.3D1*DSQRT(0.
     &6D1)+0.12D2*DSQRT(0.3D1)*Y**2.D0+0.12D2*X**2.D0*DSQRT(0.3D
     &1)-Y*(0.12D2+0.3D1*DSQRT(0.2D1)-0.24D2*DSQRT(0.3D1)*Z)
     &+0.3D1*X*(DSQRT(0.2D1)+0.2D1*DSQRT(0.3D1)*(0.2D1+DSQRT(
     &0.2D1))+0.8D1*DSQRT(0.3D1)*Z))/0.12D2
      DSP(3,10)=(-0.1D1+DSQRT(0.2D1))*(-0.12D2*Z+0.3D1*DSQRT(0
     &.6D1)+0.12D2*X**2.D0*DSQRT(0.3D1)+0.12D2*DSQRT(0.3D1)*Y*
     &*2+0.3D1*X*(DSQRT(0.2D1)+0.4D1*DSQRT(0.3D1)*Z)+0.12D2
     &*X*Z*DSQRT(0.3D1)+Y*(0.24D2+0.15D2*DSQRT(0.2D1)+0.24
     &D2*DSQRT(0.3D1)*Z))/0.12D2
      DSP(3,11)=-(-0.1D1+DSQRT(0.2D1))*(-0.4D1*DSQRT(0.3D1)+DSQRT(
     &0.6D1)+0.12D2*DSQRT(0.2D1)*Z+0.12D2*DSQRT(0.3D1)*Y**2.D0
     &+0.12D2*X**2.D0*DSQRT(0.3D1)+0.3D1*Y*(DSQRT(0.2D1)+0.4D
     &1*DSQRT(0.3D1)*Z)+0.12D2*Y*Z*DSQRT(0.3D1)+0.3D1*X*
     &(DSQRT(0.2D1)+0.8D1*DSQRT(0.3D1)*Z))/0.3D1

      END
      
      SUBROUTINE DDSHAPE_3D_11ND_TRI(NODE,X,Y,Z,DDSP)

      IMPLICIT NONE
      INTEGER::NODE
      REAL*8::X,Y,Z
      REAL*8::DDSP(3,3,NODE)
    
      DDSP(1,1,1)=-(-1.D0+DSQRT(2.D0))*(-DSQRT(2.D0)/2.D0+DSQRT(3.D0)*
     &(Y+Z))
      DDSP(1,1,2)=-(-1.D0+DSQRT(2.D0))*(-DSQRT(2.D0)/2.D0+DSQRT(3.D0)*
     &(Y+Z))
      DDSP(1,1,3)=-(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0+DSQRT(3.D0
     &)*(Y+Z))
      DDSP(1,1,4)=-(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0+DSQRT(3.D0
     &)*(Y+Z))
      DDSP(1,1,5)=2.D0*(-1.D0+DSQRT(2.D0))*(-DSQRT(2.D0)/2.D0+DSQRT(3.D0
     &)*(Y+Z))
      DDSP(1,1,6)=2.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0+DSQRT(
     &3.D0)*(Y+Z))
      DDSP(1,1,7)=2.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0+DSQRT(
     &3.D0)*(Y+Z))
      DDSP(1,1,8)=2.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0+DSQRT(
     &3.D0)*(Y+Z))
      DDSP(1,1,9)=2.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0+DSQRT(
     &3.D0)*(Y+Z))
      DDSP(1,1,10)=2.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+ 1.D0+
     &DSQRT(3.D0)*(Y+Z))
      DDSP(1,1,11)=-8.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+ 1.D0+
     &DSQRT(3.D0)*(Y+Z))
      DDSP(1,2,1)=-(-1.D0+DSQRT(2.D0))*(1.D0/4.D0+(X+Y-DSQRT(2.D0)/3.D0-
     &1.D0/3.D0)*DSQRT(3.D0))
      DDSP(1,2,2)=-(-1.D0+DSQRT(2.D0))*(1.D0/4.D0+(X+Y+DSQRT(2.D0)/3.D0+
     &1.D0/3.D0)*DSQRT(3.D0))
      DDSP(1,2,3)=-(1.D0/4.D0+DSQRT(3.D0)*(X+Y))*(-1.D0+DSQRT(2.D0))
      DDSP(1,2,4)=-(1.D0/4.D0+DSQRT(3.D0)*(X+Y))*(-1.D0+DSQRT(2.D0))
      DDSP(1,2,5)=2.D0*(1.D0/4.D0+DSQRT(3.D0)*(X+Y))*(-1.D0+DSQRT(2.D0))
      DDSP(1,2,6)=2.D0*(-1.D0+DSQRT(2.D0))*(1.D0/4.D0+(X+Y+DSQRT(2.D0)/
     &3.D0+1.D0/3.D0)*DSQRT(3.D0))
      DDSP(1,2,7)=2.D0*(-1.D0+DSQRT(2.D0))*(1.D0/4.D0+(X+Y-DSQRT(2.D0)/
     &3.D0-1.D0/3.D0)*DSQRT(3.D0))
      DDSP(1,2,8)=2.D0*(1.D0/4.D0+DSQRT(3.D0)*(X+Y))*(-1.D0+DSQRT(2.D0))
      DDSP(1,2,9)=2.D0*(1.D0/4.D0+DSQRT(3.D0)*(X+Y))*(-1.D0+DSQRT(2.D0))
      DDSP(1,2,10)=2.D0*(1.D0/4.D0+DSQRT(3.D0)*(X+Y))*(-1.D0+
     &DSQRT(2.D0))
      DDSP(1,2,11)=-8.D0*(1.D0/4.D0+DSQRT(3.D0)*(X+Y))*(-1.D0+
     &DSQRT(2.D0))
      DDSP(1,3,1)=-(-1.D0+DSQRT(2.D0))*((X+Z-DSQRT(2.D0)/6.D0-1.D0/3.D0)
     &*DSQRT(3.D0)+DSQRT(2.D0)/8.D0)
      DDSP(1,3,2)=-((X+Z+DSQRT(2.D0)/6.D0+1.D0/3.D0)*DSQRT(3.D0)+DSQRT(
     &2.D0)/8.D0)*(-1.D0+DSQRT(2.D0))
      DDSP(1,3,3)=-(DSQRT(2.D0)/8.D0+DSQRT(3.D0)*(X+Z))*(-1.D0+DSQRT(
     &2.D0))
      DDSP(1,3,4)=-(DSQRT(2.D0)/8.D0+DSQRT(3.D0)*(X+Z))*(-1.D0+DSQRT(
     &2.D0))
      DDSP(1,3,5)=2.D0*(DSQRT(2.D0)/8.D0+DSQRT(3.D0)*(X+Z))*(-1.D0+
     &DSQRT(2.D0))
      DDSP(1,3,6)=2.D0*(-1.D0+DSQRT(2.D0))*((X+Z-DSQRT(2.D0)/0.12D2-1.D0
     &/6.D0)*DSQRT(3.D0)+DSQRT(2.D0)/8.D0)
      DDSP(1,3,7)=(-1.D0+DSQRT(2.D0))*DSQRT(3.D0)*(DSQRT(2.D0)*DSQRT(
     &3.D0)+2.D0*DSQRT(2.D0)+0.24D2*X+0.24D2*Z+4.D0)/12.D0
      DDSP(1,3,8)=2.D0*(-1.D0+DSQRT(2.D0))*((X+Z-DSQRT(2.D0)/4.D0-1.D0/
     &2.D0)*DSQRT(3.D0)+DSQRT(2.D0)/8.D0)
      DDSP(1,3,9)=2.D0*(-1.D0+DSQRT(2.D0))*((X+Z+DSQRT(2.D0)/4.D0+1.D0/
     &2.D0)*DSQRT(3.D0)+DSQRT(2.D0)/8.D0)
      DDSP(1,3,10)=2.D0*(DSQRT(2.D0)/8.D0+DSQRT(3.D0)*(X+Z))*(-1.D0+
     &DSQRT(2.D0))
      DDSP(1,3,11)=-8.D0*(DSQRT(2.D0)/8.D0+DSQRT(3.D0)*(X+Z))*(-1.D0+
     &DSQRT(2.D0))
   
      DDSP(2,1,1)=DDSP(1,2,1)
      DDSP(2,1,2)=DDSP(1,2,2)
      DDSP(2,1,3)=DDSP(1,2,3)
      DDSP(2,1,4)=DDSP(1,2,4)
      DDSP(2,1,5)=DDSP(1,2,5)
      DDSP(2,1,6)=DDSP(1,2,6)
      DDSP(2,1,7)=DDSP(1,2,7)
      DDSP(2,1,8)=DDSP(1,2,8)
      DDSP(2,1,9)=DDSP(1,2,9)
      DDSP(2,1,10)=DDSP(1,2,10)
      DDSP(2,1,11)=DDSP(1,2,11)
      DDSP(2,2,1)=-(DSQRT(2.D0)/6.D0+1.D0/6.D0+DSQRT(3.D0)*(X+Z))*(-1.D0
     &+DSQRT(2.D0))
      DDSP(2,2,2)=-(DSQRT(2.D0)/6.D0+1.D0/6.D0+DSQRT(3.D0)*(X+Z))*(-1.D0
     &+DSQRT(2.D0))
      DDSP(2,2,3)=-(-1.D0+DSQRT(2.D0))*(-5.D0/6.D0*DSQRT(2.D0)-5.D0/6.D0
     &+DSQRT(3.D0)*(X+Z))
      DDSP(2,2,4)=-(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0/2.D0+
     &DSQRT(3.D0)*(X+Z))
      DDSP(2,2,5)=2.D0*(-1.D0+DSQRT(2.D0))*(5.D0/6.D0*DSQRT(2.D0)+5.D0/
     &6.D0+DSQRT(3.D0)*(X+Z))
      DDSP(2,2,6)=2.D0*(-1.D0+DSQRT(2.D0))*(-DSQRT(2.D0)/6.D0-1.D0/6.D0+
     &DSQRT(3.D0)*(X+Z))
      DDSP(2,2,7)=2.D0*(-1.D0+DSQRT(2.D0))*(-DSQRT(2.D0)/6.D0-1.D0/6.D0+
     &DSQRT(3.D0)*(X+Z))
      DDSP(2,2,8)=2.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0/2.D0+
     &DSQRT(3.D0)*(X+Z))
      DDSP(2,2,9)=2.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0/2.D0+
     &DSQRT(3.D0)*(X+Z))
      DDSP(2,2,10)=2.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0/2.D0+
     &DSQRT(3.D0)*(X+Z))
      DDSP(2,2,11)=-8.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/2.D0+1.D0/2.D0
     &+DSQRT(3.D0)*(X+Z))
      DDSP(2,3,1)=-(-1.D0+DSQRT(2.D0))*(-DSQRT(2.D0)/0.24D2-1.D0/3.D0+
     &DSQRT(3.D0)*(Y+Z))
      DDSP(2,3,2)=-(-1.D0+DSQRT(2.D0))*(-DSQRT(2.D0)/0.24D2-1.D0/3.D0+
     &DSQRT(3.D0)*(Y+Z))
      DDSP(2,3,3)=-(-1.D0+DSQRT(2.D0))*(0.11D2/0.24D2*DSQRT(2.D0)+2.D0/
     &3.D0+DSQRT(3.D0)*(Y+Z))
      DDSP(2,3,4)=-(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/8.D0+DSQRT(3.D0)*(Y+
     &Z))
      DDSP(2,3,5)=2.D0*(-1.D0+DSQRT(2.D0))*(7.D0/0.24D2*DSQRT(2.D0)+1.D0
     &/3.D0+DSQRT(3.D0)*(Y+Z))
      DDSP(2,3,6)=2.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/0.24D2-1.D0/6.D0
     &+DSQRT(3.D0)*(Y+Z))
      DDSP(2,3,7)=2.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/0.24D2-1.D0/6.D0
     &+DSQRT(3.D0)*(Y+Z))
      DDSP(2,3,8)=2.D0*(-1.D0+DSQRT(2.D0))*(-DSQRT(2.D0)/8.D0-1.D0/2.D0+
     &DSQRT(3.D0)*(Y+Z))
      DDSP(2,3,9)=2.D0*(-1.D0+DSQRT(2.D0))*(-DSQRT(2.D0)/8.D0-1.D0/2.D0+
     &DSQRT(3.D0)*(Y+Z))
      DDSP(2,3,10)=2.D0*(-1.D0+DSQRT(2.D0))*(5.D0/8.D0*DSQRT(2.D0)+1.D0+
     &DSQRT(3.D0)*(Y+Z))
      DDSP(2,3,11)=-8.D0*(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/8.D0+DSQRT(
     &3.D0)*(Y+Z))
      
      DDSP(3,1,1)=DDSP(1,3,1)
      DDSP(3,1,2)=DDSP(1,3,2)
      DDSP(3,1,3)=DDSP(1,3,3)
      DDSP(3,1,4)=DDSP(1,3,4)
      DDSP(3,1,5)=DDSP(1,3,5)
      DDSP(3,1,6)=DDSP(1,3,6)
      DDSP(3,1,7)=DDSP(1,3,7)
      DDSP(3,1,8)=DDSP(1,3,8)
      DDSP(3,1,9)=DDSP(1,3,9)
      DDSP(3,1,10)=DDSP(1,3,10)
      DDSP(3,1,11)=DDSP(1,3,11)
      DDSP(3,2,1)=DDSP(2,3,1)
      DDSP(3,2,2)=DDSP(2,3,2)
      DDSP(3,2,3)=DDSP(2,3,3)
      DDSP(3,2,4)=DDSP(2,3,4)
      DDSP(3,2,5)=DDSP(2,3,5)
      DDSP(3,2,6)=DDSP(2,3,6)
      DDSP(3,2,7)=DDSP(2,3,7)
      DDSP(3,2,8)=DDSP(2,3,8)
      DDSP(3,2,9)=DDSP(2,3,9)
      DDSP(3,2,10)=DDSP(2,3,10)
      DDSP(3,2,11)=DDSP(2,3,11)
      DDSP(3,3,1)=-(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/3.D0-1.D0/6.D0+
     &DSQRT(3.D0)*(X+Y))
      DDSP(3,3,2)=-(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/3.D0-1.D0/6.D0+
     &DSQRT(3.D0)*(X+Y))
      DDSP(3,3,3)=-(-1.D0+DSQRT(2.D0))*(DSQRT(2.D0)/3.D0-1.D0/6.D0+
     &DSQRT(3.D0)*(X+Y))
      DDSP(3,3,4)=-(-1.D0+DSQRT(2.D0))*(-DSQRT(2.D0)-3.D0/2.D0+DSQRT(
     &3.D0)*(X+Y))
      DDSP(3,3,5)=2.D0*(-1.D0+DSQRT(2.D0))*(2.D0/3.D0*DSQRT(2.D0)+1.D0/
     &6.D0+DSQRT(3.D0)*(X+Y))
      DDSP(3,3,6)=2.D0*(-1.D0+DSQRT(2.D0))*(2.D0/3.D0*DSQRT(2.D0)+1.D0/
     &6.D0+DSQRT(3.D0)*(X+Y))
      DDSP(3,3,7)=2.D0*(-1.D0+DSQRT(2.D0))*(2.D0/3.D0*DSQRT(2.D0)+1.D0/
     &6.D0+DSQRT(3.D0)*(X+Y))
      DDSP(3,3,8)=2.D0*(-1.D0+DSQRT(2.D0))*(-1.D0/2.D0+DSQRT(3.D0)*
     &(X+Y))
      DDSP(3,3,9)=2.D0*(-1.D0+DSQRT(2.D0))*(-1.D0/2.D0+DSQRT(3.D0)*
     &(X+Y))
      DDSP(3,3,10)=2.D0*(-1.D0+DSQRT(2.D0))*(-1.D0/2.D0+DSQRT(3.D0)*
     &(X+Y))
      DDSP(3,3,11)=-8.D0*(DSQRT(2.D0)/2.D0+DSQRT(3.D0)*(X+Y))*(-1.D0+
     &DSQRT(2.D0))

      END
      
      SUBROUTINE SHAPE_2D_8ND(NODE,X,Y,SP) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SP(NODE)
      SS=X*X; TT=Y*Y
      ST=X*Y; SST=X*X*Y
      STT=X*Y*Y
      SP(1)=(ST+SS+TT-SST-STT-1.0D0)*0.25D0
      SP(2)=(SS+TT-SST+STT-1.0D0-ST)*0.25D0
      SP(3)=(ST+SS+TT+SST+STT-1.0D0)*0.25D0
      SP(4)=(SS+TT+SST-STT-1.0D0-ST)*0.25D0
      SP(5)=(1.0D0-Y-SS+SST)*0.5D0
      SP(6)=(1.0D0+X-TT-STT)*0.5D0
      SP(7)=(1.0D0+Y-SS-SST)*0.5D0
      SP(8)=(1.0D0-X-TT+STT)*0.5D0
      END

      SUBROUTINE DSHAPE_2D_8ND(NODE,X,Y,DN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DN(3,NODE)
       X2=X*2.0D0
       T2=Y*2.0D0
       SS=X*X
       TT=Y*Y
       ST=X*Y
       ST2=X*Y*2.0D0
       DN(1,1)=(Y+X2-ST2-TT)*0.25D0
       DN(1,2)=(X2-ST2+TT-Y)*0.25D0
       DN(1,3)=(Y+X2+ST2+TT)*0.25D0
       DN(1,4)=(X2+ST2-TT-Y)*0.25D0
       DN(1,5)=ST-X
       DN(1,6)=(1.0-TT)*0.5D0
       DN(1,7)=-X-ST
       DN(1,8)=(TT-1.0)*0.5D0
       DN(2,1)=(X+T2-SS-ST2)*0.25D0
       DN(2,2)=(T2-SS+ST2-X)*0.25D0
       DN(2,3)=(X+T2+SS+ST2)*0.25D0
       DN(2,4)=(T2+SS-ST2-X)*0.25D0
       DN(2,5)=(SS-1.0D0)*0.5D0
       DN(2,6)=-Y-ST
       DN(2,7)=(1.0D0-SS)*0.5D0
       DN(2,8)=ST-Y
      END

      SUBROUTINE DDSHAPE_2D_8ND(NODE,X,Y,DDN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DDN(3,3,NODE)
      DDN(1,1,1)=0.5D0*(1.-Y)
      DDN(1,1,2)=0.5D0*(1.-Y)
      DDN(1,1,3)=0.5D0*(1.+Y)
      DDN(1,1,4)=0.5D0*(1.+Y)
      DDN(1,1,5)=-(1.D0-Y)
      DDN(1,1,6)=0.D0
      DDN(1,1,7)=-(1.D0+Y)
      DDN(1,1,8)=0.D0

      DDN(2,2,1)=0.5D0*(1.D0-X)
      DDN(2,2,2)=0.5D0*(1.D0+X)
      DDN(2,2,3)=0.5D0*(1.D0+X)
      DDN(2,2,4)=0.5D0*(1.D0-X)
      DDN(2,2,5)=0.D0
      DDN(2,2,6)=-(1.D0+X)
      DDN(2,2,7)=0.D0
      DDN(2,2,8)=-(1.D0-X)

      DDN(2,1,1)=0.25D0*(1.D0-2.D0*X-2.D0*Y)
      DDN(2,1,2)=0.25D0*(-1.D0-2.D0*X+2.D0*Y)
      DDN(2,1,3)=0.25D0*(1.D0+2.D0*X+2.D0*Y)
      DDN(2,1,4)=0.25D0*(-1.D0+2.D0*X-2.D0*Y)
      DDN(2,1,5)=X
      DDN(2,1,6)=-Y
      DDN(2,1,7)=-X
      DDN(2,1,8)=Y
            
      DDN(1,2,:)=DDN(2,1,:)
      END

      SUBROUTINE SHAPE_3D_20ND(NODE,X,Y,Z,SP) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SP(NODE)
      SP(1)=((1.D0-X)*(1.D0-Y)*(1.D0-Z)*(-2.D0-X-Y-Z))/8.D0
      SP(2)=((1.D0+X)*(1.D0-Y)*(1.D0-Z)*(-2.D0+X-Y-Z))/8.D0
      SP(3)=((1.D0+X)*(1.D0+Y)*(1.D0-Z)*(-2.D0+X+Y-Z))/8.D0
      SP(4)=((1.D0-X)*(1.D0+Y)*(1.D0-Z)*(-2.D0-X+Y-Z))/8.D0
      SP(5)=((1.D0-X)*(1.D0-Y)*(1.D0+Z)*(-2.D0-X-Y+Z))/8.D0
      SP(6)=((1.D0+X)*(1.D0-Y)*(1.D0+Z)*(-2.D0+X-Y+Z))/8.D0
      SP(7)=((1.D0+X)*(1.D0+Y)*(1.D0+Z)*(-2.D0+X+Y+Z))/8.D0
      SP(8)=((1.D0-X)*(1.D0+Y)*(1.D0+Z)*(-2.D0-X+Y+Z))/8.D0
      SP(9)=((1.D0-X**2.D0)*(1.D0-Y)*(1.D0-Z))/4.D0
      SP(10)=((1.D0+X)*(1.D0-Y**2.D0)*(1.D0-Z))/4.D0
      SP(11)=((1.D0-X**2.D0)*(1.D0+Y)*(1.D0-Z))/4.D0
      SP(12)=((1.D0-X)*(1.D0-Y**2.D0)*(1.D0-Z))/4.D0
      SP(13)=((1.D0-X**2.D0)*(1.D0-Y)*(1.D0+Z))/4.D0
      SP(14)=((1.D0+X)*(1.D0-Y**2.D0)*(1.D0+Z))/4.D0
      SP(15)=((1.D0-X**2.D0)*(1.D0+Y)*(1.D0+Z))/4.D0
      SP(16)=((1.D0-X)*(1.D0-Y**2.D0)*(1.D0+Z))/4.D0
      SP(17)=((1.D0-X)*(1.D0-Y)*(1.D0-Z**2.D0))/4.D0
      SP(18)=((1.D0+X)*(1.D0-Y)*(1.D0-Z**2.D0))/4.D0
      SP(19)=((1.D0+X)*(1.D0+Y)*(1.D0-Z**2.D0))/4.D0
      SP(20)=((1.D0-X)*(1.D0+Y)*(1.D0-Z**2.D0))/4.D0
      END

      SUBROUTINE DSHAPE_3D_20ND(NODE,X,Y,Z,DN) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DN(3,NODE)
      DN(1,1)=((-1.D0+Y)*(-1.D0+Z)*(1.D0+2.D0*X+Y+Z))/8.D0
      DN(1,2)=-((-1.D0+Y)*(-1.D0+Z)*(1.D0-2.D0*X+Y+Z))/8.D0
      DN(1,3)=((1.D0+Y)*(-1.D0+Z)*(1.D0-2.D0*X-Y+Z))/8.D0
      DN(1,4)=-((1.D0+Y)*(-1.D0+Z)*(1.D0+2.D0*X-Y+Z))/8.D0
      DN(1,5)=((-1.D0+Y)*(1.D0+Z)*(-1.D0-2.D0*X-Y+Z))/8.D0
      DN(1,6)=-((-1.D0+Y)*(1.D0+Z)*(-1.D0+2.D0*X-Y+Z))/8.D0
      DN(1,7)=((1.D0+Y)*(1.D0+Z)*(-1.D0+2.D0*X+Y+Z))/8.D0
      DN(1,8)=-((1.D0+Y)*(1.D0+Z)*(-1.D0-2.D0*X+Y+Z))/8.D0
      DN(1,9)=-(X*(-1.D0+Y)*(-1.D0+Z))/2.D0
      DN(1,10)=((-1.D0+Y**2.D0)*(-1.D0+Z))/4.D0
      DN(1,11)=(X*(1.D0+Y)*(-1.D0+Z))/2.D0
      DN(1,12)=-((-1.D0+Y**2.D0)*(-1.D0+Z))/4.D0
      DN(1,13)=(X*(-1.D0+Y)*(1.D0+Z))/2.D0
      DN(1,14)=-((-1.D0+Y**2.D0)*(1.D0+Z))/4.D0
      DN(1,15)=-(X*(1.D0+Y)*(1.D0+Z))/2.D0
      DN(1,16)=((-1.D0+Y**2.D0)*(1.D0+Z))/4.D0
      DN(1,17)=-((-1.D0+Y)*(-1.D0+Z**2.D0))/4.D0
      DN(1,18)=((-1.D0+Y)*(-1.D0+Z**2.D0))/4.D0
      DN(1,19)=-((1.D0+Y)*(-1.D0+Z**2.D0))/4.D0
      DN(1,20)=((1.D0+Y)*(-1.D0+Z**2.D0))/4.D0
      
      DN(2,1)=((-1.D0+X)*(-1.D0+Z)*(1.D0+X+2.D0*Y+Z))/8.D0
      DN(2,2)=-((1.D0+X)*(-1.D0+Z)*(1.D0-X+2.D0*Y+Z))/8.D0
      DN(2,3)=((1.D0+X)*(-1.D0+Z)*(1.D0-X-2.D0*Y+Z))/8.D0
      DN(2,4)=-((-1.D0+X)*(-1.D0+Z)*(1.D0+X-2.D0*Y+Z))/8.D0
      DN(2,5)=((-1.D0+X)*(1.D0+Z)*(-1.D0-X-2.D0*Y+Z))/8.D0
      DN(2,6)=-((1.D0+X)*(1.D0+Z)*(-1.D0+X-2.D0*Y+Z))/8.D0
      DN(2,7)=((1.D0+X)*(1.D0+Z)*(-1.D0+X+2.D0*Y+Z))/8.D0
      DN(2,8)=-((-1.D0+X)*(1.D0+Z)*(-1.D0-X+2.D0*Y+Z))/8.D0
      DN(2,9)=-((-1.D0+X**2.D0)*(-1.D0+Z))/4.D0
      DN(2,10)=((1.D0+X)*Y*(-1.D0+Z))/2.D0
      DN(2,11)=((-1.D0+X**2.D0)*(-1.D0+Z))/4.D0
      DN(2,12)=-((-1.D0+X)*Y*(-1.D0+Z))/2.D0
      DN(2,13)=((-1.D0+X**2.D0)*(1.D0+Z))/4.D0
      DN(2,14)=-((1.D0+X)*Y*(1.D0+Z))/2.D0
      DN(2,15)=-((-1.D0+X**2.D0)*(1.D0+Z))/4.D0
      DN(2,16)=((-1.D0+X)*Y*(1.D0+Z))/2.D0
      DN(2,17)=-((-1.D0+X)*(-1.D0+Z**2.D0))/4.D0
      DN(2,18)=((1.D0+X)*(-1.D0+Z**2.D0))/4.D0
      DN(2,19)=-((1.D0+X)*(-1.D0+Z**2.D0))/4.D0
      DN(2,20)=((-1.D0+X)*(-1.D0+Z**2.D0))/4.D0
      
      DN(3,1)=((-1.D0+X)*(-1.D0+Y)*(1.D0+X+Y+2.D0*Z))/8.D0
      DN(3,2)=-((1.D0+X)*(-1.D0+Y)*(1.D0-X+Y+2.D0*Z))/8.D0
      DN(3,3)=-((1.D0+X)*(1.D0+Y)*(-1.D0+X+Y-2.D0*Z))/8.D0
      DN(3,4)=((-1.D0+X)*(1.D0+Y)*(-1.D0-X+Y-2.D0*Z))/8.D0
      DN(3,5)=-((-1.D0+X)*(-1.D0+Y)*(1.D0+X+Y-2.D0*Z))/8.D0
      DN(3,6)=((1.D0+X)*(-1.D0+Y)*(1.D0-X+Y-2.D0*Z))/8.D0
      DN(3,7)=((1.D0+X)*(1.D0+Y)*(-1.D0+X+Y+2.D0*Z))/8.D0
      DN(3,8)=-((-1.D0+X)*(1.D0+Y)*(-1.D0-X+Y+2.D0*Z))/8.D0
      DN(3,9)=-((-1.D0+X**2.D0)*(-1.D0+Y))/4.D0
      DN(3,10)=((1.D0+X)*(-1.D0+Y**2.D0))/4.D0
      DN(3,11)=((-1.D0+X**2.D0)*(1.D0+Y))/4.D0
      DN(3,12)=-((-1.D0+X)*(-1.D0+Y**2.D0))/4.D0
      DN(3,13)=((-1.D0+X**2.D0)*(-1.D0+Y))/4.D0
      DN(3,14)=-((1.D0+X)*(-1.D0+Y**2.D0))/4.D0
      DN(3,15)=-((-1.D0+X**2.D0)*(1.D0+Y))/4.D0
      DN(3,16)=((-1.D0+X)*(-1.D0+Y**2.D0))/4.D0
      DN(3,17)=-((-1.D0+X)*(-1.D0+Y)*Z)/2.D0
      DN(3,18)=((1.D0+X)*(-1.D0+Y)*Z)/2.D0
      DN(3,19)=-((1.D0+X)*(1.D0+Y)*Z)/2.D0
      DN(3,20)=((-1.D0+X)*(1.D0+Y)*Z)/2.D0
      END
      
      SUBROUTINE SHAPE_3D_21ND(NNODE,X,Y,Z,SP)
      !ÐÎº¯ÊýN_A
      IMPLICIT NONE
      INTEGER::NNODE
      REAL*8::X,Y,Z
      REAL*8,DIMENSION(NNODE)::SP
      SP(1) =(-((-3 + Z)*Z) + 3*X*(-1 + Y)*(-1 + Z)*(1 + Y + Z) - 
     -     3*Y*(-1 + Z**2) + Y**2*(-1 + Z*(-3 + 2*Z)) + 
     -     X**2*(-1 + 2*Y**2 + 3*Y*(-1 + Z) + Z*(-3 + 2*Z)))/24.D0
      SP(2) =(-((-3 + Z)*Z) - 3*X*(-1 + Y)*(-1 + Z)*(1 + Y + Z) - 
     -     3*Y*(-1 + Z**2) + Y**2*(-1 + Z*(-3 + 2*Z)) + 
     -     X**2*(-1 + 2*Y**2 + 3*Y*(-1 + Z) + Z*(-3 + 2*Z)))/24.D0
      SP(3) =(-(Y*(3 + Y)) - 3*(-1 + Y**2)*Z+(-1 + Y*(3 + 2*Y))*Z**2+ 
     -     3*X*(1 + Y)*(-1 + Z)*(1 - Y + Z) + 
     -     X**2*(-1 + 3*Y + 2*Y**2 - 3*(1 + Y)*Z + 2*Z**2))/24.D0
      SP(4) =(-(Y*(3 + Y)) + 3*X*(1 + Y)*(-1 + Y - Z)*(-1 + Z) - 
     -     3*(-1 + Y**2)*Z + (-1 + Y*(3 + 2*Y))*Z**2 + 
     -     X**2*(-1 + 3*Y + 2*Y**2 - 3*(1 + Y)*Z + 2*Z**2))/24.D0
      SP(5) =(-3*X*(-1 + Y)*(1 + Y - Z)*(1 + Z) - Z*(3 + Z) - 
     -     3*Y*(-1 + Z**2) + Y**2*(-1 + Z*(3 + 2*Z)) + 
     -     X**2*(-1 + 2*Y**2 - 3*Y*(1 + Z) + Z*(3 + 2*Z)))/24.D0
      SP(6) =(3*X*(-1 + Y)*(1 + Y - Z)*(1 + Z) - Z*(3 + Z) - 
     -     3*Y*(-1 + Z**2) + Y**2*(-1 + Z*(3 + 2*Z)) + 
     -     X**2*(-1 + 2*Y**2 - 3*Y*(1 + Z) + Z*(3 + 2*Z)))/24.D0
      SP(7) =(-(Z*(3 + Z)) + 3*X*(1 + Y)*(1 + Z)*(-1 + Y + Z) + 
     -     3*Y*(-1 + Z**2) + Y**2*(-1 + Z*(3 + 2*Z)) + 
     -     X**2*(-1 + 2*Y**2 + 3*Y*(1 + Z) + Z*(3 + 2*Z)))/24.D0
      SP(8) =(-(Z*(3 + Z)) - 3*X*(1 + Y)*(1 + Z)*(-1 + Y + Z) + 
     -     3*Y*(-1 + Z**2) + Y**2*(-1 + Z*(3 + 2*Z)) + 
     -     X**2*(-1 + 2*Y**2 + 3*Y*(1 + Z) + Z*(3 + 2*Z)))/24.D0
      SP(9) =(3*Y*(-1 + Z) + Z*(-3 + 2*Z) - 
     -  X**2*(1 + Y**2 + 3*Y*(-1 + Z) + (-3 + Z)*Z) - 
     &  Y**2*(-2 + Z**2))/12.D0
      SP(10) =(3*X*(-1 + Y**2)*(-1 + Z) + Z*(-3 + 2*Z) - 
     -     Y**2*(1 + (-3 + Z)*Z) - X**2*(-2 + Y**2 + Z**2))/12.D0
      SP(11) =(Y*(3 + 2*Y) - 3*(1 + Y)*Z - (-2 + Y**2)*Z**2 - 
     -     X**2*(1 + 3*Y + Y**2 - 3*(1 + Y)*Z + Z**2))/12.D0
      SP(12) =(-3*X*(-1 + Y**2)*(-1 + Z) + Z*(-3 + 2*Z) - 
     -     Y**2*(1 + (-3 + Z)*Z) - X**2*(-2 + Y**2 + Z**2))/12.D0
      SP(13) =(-3*Y*(1 + Z) + Z*(3 + 2*Z) - Y**2*(-2 + Z**2) - 
     -     X**2*(1 + Y**2 - 3*Y*(1 + Z) + Z*(3 + Z)))/12.D0
      SP(14) =(-3*X*(-1 + Y**2)*(1 + Z) + Z*(3 + 2*Z) - 
     -     X**2*(-2 + Y**2 + Z**2) - Y**2*(1 + Z*(3 + Z)))/12.D0
      SP(15) =(3*Y*(1 + Z) + Z*(3 + 2*Z) - Y**2*(-2 + Z**2) - 
     -     X**2*(1 + Y**2 + 3*Y*(1 + Z) + Z*(3 + Z)))/12.D0
      SP(16) =(3*X*(-1 + Y**2)*(1 + Z) + Z*(3 + 2*Z) - 
     -     X**2*(-2 + Y**2 + Z**2) - Y**2*(1 + Z*(3 + Z)))/12.D0
      SP(17) =(Y*(-3 + 2*Y) - (1 + (-3 + Y)*Y)*Z**2 - 
     -     3*X*(-1 + Y)*(-1 + Z**2) - X**2*(-2 + Y**2 + Z**2))/12.D0
      SP(18) =(Y*(-3 + 2*Y) - (1 + (-3 + Y)*Y)*Z**2 + 
     -     3*X*(-1 + Y)*(-1 + Z**2) - X**2*(-2 + Y**2 + Z**2))/12.D0
      SP(19) =(Y*(3 + 2*Y) - (1 + Y*(3 + Y))*Z**2 - 
     -     3*X*(1 + Y)*(-1 + Z**2) - X**2*(-2 + Y**2 + Z**2))/12.D0
      SP(20) =(Y*(3 + 2*Y) - (1 + Y*(3 + Y))*Z**2 + 
     -     3*X*(1 + Y)*(-1 + Z**2) - X**2*(-2 + Y**2 + Z**2))/12.D0
      SP(21) = (3 - 2*Z**2 + Y**2*(-2 + Z**2) +
     &     X**2*(-2 + Y**2 + Z**2))/3.D0             
      END
      
      SUBROUTINE DSHAPE_3D_21ND(NNODE,X,Y,Z,DSP)
      !DNA/D(X,Y,Z)
      IMPLICIT NONE
      INTEGER::NNODE
      REAL*8::X,Y,Z
      REAL*8,DIMENSION(3,NNODE)::DSP!
      DSP(1,1)=(3*(-1 + Y)*(-1 + Z)*(1 + Y + Z) + 
     -     X*(-2 + 4*Y**2 + 6*Y*(-1 + Z) - 6*Z + 4*Z**2))/24.D0
      DSP(1,2)=(-3*(-1 + Y)*(-1 + Z)*(1 + Y + Z) + 
     -     X*(-2 + 4*Y**2 + 6*Y*(-1 + Z) - 6*Z + 4*Z**2))/24.D0
      DSP(1,3)=(3*(1 + Y)*(-1 + Z)*(1 - Y + Z) + 
     -     2*X*(-1 + 3*Y + 2*Y**2 - 3*(1 + Y)*Z + 2*Z**2))/24.D0
      DSP(1,4)=(3*(1 + Y)*(-1 + Y - Z)*(-1 + Z) + 
     -     2*X*(-1 + 3*Y + 2*Y**2 - 3*(1 + Y)*Z + 2*Z**2))/24.D0
      DSP(1,5)=(-3*(-1 + Y)*(1 + Y - Z)*(1 + Z) + 
     -     X*(-2 + 4*Y**2 + 6*Z + 4*Z**2 - 6*Y*(1 + Z)))/24.D0
      DSP(1,6)=(3*(-1 + Y)*(1 + Y - Z)*(1 + Z) + 
     -     X*(-2 + 4*Y**2 + 6*Z + 4*Z**2 - 6*Y*(1 + Z)))/24.D0
      DSP(1,7)=(3*(1 + Y)*(1 + Z)*(-1 + Y + Z) + 
     -     X*(-2 + 4*Y**2 + 6*Z + 4*Z**2 + 6*Y*(1 + Z)))/24.D0
      DSP(1,8)=(-3*(1 + Y)*(1 + Z)*(-1 + Y + Z) + 
     -     X*(-2 + 4*Y**2 + 6*Z + 4*Z**2 + 6*Y*(1 + Z)))/24.D0
      DSP(1,9)=-(X*(1 + Y**2 + 3*Y*(-1 + Z) + (-3 + Z)*Z))/6.D0
      DSP(1,10)=(3*(-1 + Y**2)*(-1 + Z) - 2*X*(-2 + Y**2 + Z**2))/12.D0
      DSP(1,11)=-(X*(1 + 3*Y + Y**2 - 3*(1 + Y)*Z + Z**2))/6.D0
      DSP(1,12)=(-3*(-1 + Y**2)*(-1 + Z) - 2*X*(-2 + Y**2 + Z**2))/12.D0
      DSP(1,13)=-(X*(1 + Y**2 - 3*Y*(1 + Z) + Z*(3 + Z)))/6.D0
      DSP(1,14)=(-3*(-1 + Y**2)*(1 + Z) - 2*X*(-2 + Y**2 + Z**2))/12.D0
      DSP(1,15)=-(X*(1 + Y**2 + 3*Y*(1 + Z) + Z*(3 + Z)))/6.D0
      DSP(1,16)=(3*(-1 + Y**2)*(1 + Z) - 2*X*(-2 + Y**2 + Z**2))/12.D0
      DSP(1,17)=(-3*(-1 + Y)*(-1 + Z**2) - 2*X*(-2 + Y**2 + Z**2))/12.D0
      DSP(1,18)=(3*(-1 + Y)*(-1 + Z**2) - 2*X*(-2 + Y**2 + Z**2))/12.D0
      DSP(1,19)=(-3*(1 + Y)*(-1 + Z**2) - 2*X*(-2 + Y**2 + Z**2))/12.D0
      DSP(1,20)=(3*(1 + Y)*(-1 + Z**2) - 2*X*(-2 + Y**2 + Z**2))/12.D0
      DSP(1,21)=(2*X*(-2 + Y**2 + Z**2))/3.D0

      DSP(2,1)=(3 - 2*Y - 6*Y*Z + (-3 + 4*Y)*Z**2+3*X*(-1+Z)*(2*Y+Z) + 
     -     X**2*(-3 + 4*Y + 3*Z))/24.D0
      DSP(2,2)=(3 - 2*Y + X*(6*Y + X*(-3 + 4*Y))+3*(1+X)*(X-2*Y)*Z + 
     -     (-3 - 3*X + 4*Y)*Z**2)/24.D0
      DSP(2,3)=(-3 - 2*Y + X*(6*Y + X*(3 + 4*Y)) - 3*(1 + X)*(X+2*Y)*Z+ 
     -     (3 + 3*X + 4*Y)*Z**2)/24.D0
      DSP(2,4)=(-3 - 2*Y + X*(-6*Y + X*(3 + 4*Y))-3*(-1+X)*(X-2*Y)*Z + 
     -     (3 - 3*X + 4*Y)*Z**2)/24.D0
      DSP(2,5)=(3 - 2*Y + X*(-6*Y + X*(-3 + 4*Y))-3*(-1+X)*(X+2*Y)*Z + 
     -     (-3 + 3*X + 4*Y)*Z**2)/24.D0
      DSP(2,6)=(3 - 2*Y + X*(6*Y + X*(-3 + 4*Y))-3*(1+X)*(X-2*Y)*Z + 
     -     (-3 - 3*X + 4*Y)*Z**2)/24.D0
      DSP(2,7)=(-3 - 2*Y + 6*Y*Z + (3 + 4*Y)*Z**2+3*X*(1+Z)*(2*Y+Z) + 
     -     X**2*(3 + 4*Y + 3*Z))/24.D0
      DSP(2,8)=(-3 - 2*Y + 6*Y*Z + (3 + 4*Y)*Z**2-3*X*(1+Z)*(2*Y+Z) + 
     -     X**2*(3 + 4*Y + 3*Z))/24.D0
      DSP(2,9)=(-3 + 4*Y + X**2*(3 - 2*Y - 3*Z) + 3*Z - 2*Y*Z**2)/12.D0
      DSP(2,10)=-(Y*(1 + 3*X + X**2 - 3*(1 + X)*Z + Z**2))/6.D0
      DSP(2,11)=(3 - 3*Z + X**2*(-3 - 2*Y + 3*Z)- 2*Y*(-2 + Z**2))/12.D0
      DSP(2,12)=-(Y*(1 + X**2 + 3*X*(-1 + Z) + (-3 + Z)*Z))/6.D0
      DSP(2,13)=(-3*(1 + Z) + X**2*(3-2*Y+3*Z)-2*Y*(-2+Z**2))/12.D0
      DSP(2,14)=-(Y*(1 + X**2 + 3*X*(1 + Z) + Z*(3 + Z)))/6.D0
      DSP(2,15)=(3 + 4*Y + 3*Z - 2*Y*Z**2 - X**2*(3 + 2*Y + 3*Z))/12.D0
      DSP(2,16)=-(Y*(1 + X**2 - 3*X*(1 + Z) + Z*(3 + Z)))/6.D0
      DSP(2,17)=(-3 + 3*X + 4*Y - 2*X**2*Y + (3 - 3*X - 2*Y)*Z**2)/12.D0
      DSP(2,18)=(-3 - 3*X + 4*Y - 2*X**2*Y + (3 + 3*X - 2*Y)*Z**2)/12.D0
      DSP(2,19)=(3 + 3*X + 4*Y - 2*X**2*Y - (3 + 3*X + 2*Y)*Z**2)/12.D0
      DSP(2,20)=(3 - 3*X + 4*Y - 2*X**2*Y + (-3 + 3*X - 2*Y)*Z**2)/12.D0
      DSP(2,21)=(2*Y*(-2 + X**2 + Z**2))/3.D0
      
      DSP(3,1)=(3*(-1 + X)*(-1 + Y)*(1 + X + Y) + 
     -     2*(-1 + 2*X**2 + 3*X*(-1 + Y) + Y*(-3 + 2*Y))*Z)/24.D0
      DSP(3,2)=(3*(1 + X)*(-1 + X - Y)*(-1 + Y) + 
     -     2*(-1 + 3*X + 2*X**2 - 3*(1 + X)*Y + 2*Y**2)*Z)/24.D0
      DSP(3,3)=(-3*(1 + X)*(1 + Y)*(-1 + X + Y) + 
     -     2*(-1 + 2*X**2 + 3*X*(1 + Y) + Y*(3 + 2*Y))*Z)/24.D0
      DSP(3,4)=(-3*(-1 + X)*(1 + X - Y)*(1 + Y) + 
     -     2*(-1 + 2*X**2 - 3*X*(1 + Y) + Y*(3 + 2*Y))*Z)/24.D0
      DSP(3,5)=(-3*(-1 + X)*(-1 + Y)*(1 + X + Y) + 
     -     2*(-1 + 2*X**2 + 3*X*(-1 + Y) + Y*(-3 + 2*Y))*Z)/24.D0
      DSP(3,6)=(3*(1 + X)*(-1 + Y)*(1 - X + Y) + 
     -     2*(-1 + 3*X + 2*X**2 - 3*(1 + X)*Y + 2*Y**2)*Z)/24.D0
      DSP(3,7)=(3*(1 + X)*(1 + Y)*(-1 + X + Y) + 
     -     2*(-1 + 2*X**2 + 3*X*(1 + Y) + Y*(3 + 2*Y))*Z)/24.D0
      DSP(3,8)=(3*(-1 + X)*(1 + X - Y)*(1 + Y) + 
     -     2*(-1 + 2*X**2 - 3*X*(1 + Y) + Y*(3 + 2*Y))*Z)/24.D0
      DSP(3,9)=-((-1 + X**2)*(-1 + Y))/4.D0 -((-2 + X**2 + Y**2)*Z)/6.D0
      DSP(3,10)=((1 + X)*(-1 + Y**2))/4.D0 - ((-2 + X**2 + Y**2)*Z)/6.D0
      DSP(3,11)=((-1 + X**2)*(1 + Y))/4.D0 - ((-2 + X**2 + Y**2)*Z)/6.D0
      DSP(3,12)=-((-1 + X)*(-1 + Y**2))/4.D0-((-2 + X**2 + Y**2)*Z)/6.D0
      DSP(3,13)=((-1 + X**2)*(-1 + Y))/4.D0-((-2 + X**2 + Y**2)*Z)/6.D0
      DSP(3,14)=-((1 + X)*(-1 + Y**2))/4.D0-((-2 + X**2 + Y**2)*Z)/6.D0
      DSP(3,15)=-((-1 + X**2)*(1 + Y))/4.D0-((-2 + X**2 + Y**2)*Z)/6.D0
      DSP(3,16)=((-1 + X)*(-1 + Y**2))/4.D0-((-2 + X**2 + Y**2)*Z)/6.D0
      DSP(3,17)=-((1 + X**2 + 3*X*(-1 + Y) + (-3 + Y)*Y)*Z)/6.D0
      DSP(3,18)=-((1 + 3*X + X**2 - 3*(1 + X)*Y + Y**2)*Z)/6.D0
      DSP(3,19)=-((1 + X**2 + 3*X*(1 + Y) + Y*(3 + Y))*Z)/6.D0
      DSP(3,20)=-((1 + X**2 - 3*X*(1 + Y) + Y*(3 + Y))*Z)/6.D0
      DSP(3,21)=(2*(-2 + X**2 + Y**2)*Z)/3.D0
      END
      
      SUBROUTINE DDSHAPE_3D_21ND(NNODE,X,Y,Z,DDSP)
      !D(DNA/D(X,Y,Z))/D(X,Y,Z)
      IMPLICIT NONE
      INTEGER::NNODE,ID
      REAL*8::X,Y,Z
      REAL*8,DIMENSION(3,3,NNODE)::DDSP
      !DXDXI
      DDSP(1,1,1) = (-1 + 2*Y**2 + 3*Y*(-1 + Z)+ Z*(-3 + 2*Z))/12.D0
      DDSP(1,1,2) = (-1 + 2*Y**2 + 3*Y*(-1 + Z)+ Z*(-3 + 2*Z))/12.D0
      DDSP(1,1,3) = (-1 + 3*Y + 2*Y**2 - 3*(1 + Y)*Z + 2*Z**2)/12.D0
      DDSP(1,1,4) = (-1 + 3*Y + 2*Y**2 - 3*(1 + Y)*Z + 2*Z**2)/12.D0
      DDSP(1,1,5) = (-1 + 2*Y**2 - 3*Y*(1 + Z) + Z*(3 + 2*Z))/12.D0
      DDSP(1,1,6) = (-1 + 2*Y**2 - 3*Y*(1 + Z) + Z*(3 + 2*Z))/12.D0
      DDSP(1,1,7) = (-1 + 2*Y**2 + 3*Y*(1 + Z) + Z*(3 + 2*Z))/12.D0
      DDSP(1,1,8) = (-1 + 2*Y**2 + 3*Y*(1 + Z) + Z*(3 + 2*Z))/12.D0
      DDSP(1,1,9) = (-1 - Y**2 - 3*Y*(-1 + Z) - (-3 + Z)*Z)/6.D0
      DDSP(1,1,10) = (2 - Y**2 - Z**2)/6.D0
      DDSP(1,1,11) =(-1 - Y*(3 + Y) + 3*Z + 3*Y*Z - Z**2)/6.D0 
      DDSP(1,1,12) = (2 - Y**2 - Z**2)/6.D0
      DDSP(1,1,13) = (-1 - Y**2 + 3*Y*(1 + Z) - Z*(3 + Z))/6.D0
      DDSP(1,1,14) = (2 - Y**2 - Z**2)/6.D0
      DDSP(1,1,15) = (-1 - Y**2 - 3*Y*(1 + Z) - Z*(3 + Z))/6.D0
      DDSP(1,1,16) = (2 - Y**2 - Z**2)/6.D0
      DDSP(1,1,17) = (2 - Y**2 - Z**2)/6.D0
      DDSP(1,1,18) = (2 - Y**2 - Z**2)/6.D0
      DDSP(1,1,19) = (2 - Y**2 - Z**2)/6.D0
      DDSP(1,1,20) = (2 - Y**2 - Z**2)/6.D0
      DDSP(1,1,21) = (2*(-2 + Y**2 + Z**2))/3.D0
      
      DDSP(1,2,1) = (3*(-1 + Z)*(2*Y + Z) + X*(-6 + 8*Y + 6*Z))/24.D0
      DDSP(1,2,2) = (-3*(-1 + Z)*(2*Y + Z) + X*(-6 + 8*Y + 6*Z))/24.D0
      DDSP(1,2,3) = (X*(6 + 8*Y - 6*Z) - 3*(2*Y - Z)*(-1 + Z))/24.D0
      DDSP(1,2,4) = (X*(6 + 8*Y - 6*Z) + 3*(2*Y - Z)*(-1 + Z))/24.D0
      DDSP(1,2,5) = (-3*(2*Y - Z)*(1 + Z) + X*(8*Y - 6*(1 + Z)))/24.D0
      DDSP(1,2,6) = (3*(2*Y - Z)*(1 + Z) + X*(8*Y - 6*(1 + Z)))/24.D0
      DDSP(1,2,7) = (3*(1 + Z)*(2*Y + Z) + X*(6 + 8*Y + 6*Z))/24.D0
      DDSP(1,2,8) = (-3*(1 + Z)*(2*Y + Z) + X*(6 + 8*Y + 6*Z))/24.D0
      DDSP(1,2,9) = -(X*(-3 + 2*Y + 3*Z))/6.D0
      DDSP(1,2,10) = -(Y*(3 + 2*X - 3*Z))/6.D0
      DDSP(1,2,11) = -(X*(3 + 2*Y - 3*Z))/6.D0
      DDSP(1,2,12) = -(Y*(-3 + 2*X + 3*Z))/6.D0
      DDSP(1,2,13) = (X*(3 - 2*Y + 3*Z))/6.D0
      DDSP(1,2,14) = -(Y*(3 + 2*X + 3*Z))/6.D0
      DDSP(1,2,15) = -(X*(3 + 2*Y + 3*Z))/6.D0
      DDSP(1,2,16) = (Y*(3 - 2*X + 3*Z))/6.D0
      DDSP(1,2,17) = (3 - 4*X*Y - 3*Z**2)/12.D0
      DDSP(1,2,18) = (-3 - 4*X*Y + 3*Z**2)/12.D0
      DDSP(1,2,19) = (3 - 4*X*Y - 3*Z**2)/12.D0
      DDSP(1,2,20) = (-3 - 4*X*Y + 3*Z**2)/12.D0
      DDSP(1,2,21) = (4*X*Y)/3.D0
      
      DDSP(1,3,1) = (3*(-1 + Y)*(Y + 2*Z) + X*(-6 + 6*Y + 8*Z))/24.D0
      DDSP(1,3,2) = (-3*(-1 + Y)*(Y + 2*Z) + X*(-6 + 6*Y + 8*Z))/24.D0
      DDSP(1,3,3) = (-3*(1 + Y)*(2*X + Y) + 2*(3 + 4*X + 3*Y)*Z)/24.D0
      DDSP(1,3,4) = (3*(1 + Y)*(Y - 2*Z) + X*(-6 - 6*Y + 8*Z))/24.D0
      DDSP(1,3,5) = (-3*(-1 + Y)*(Y - 2*Z) + X*(6 - 6*Y + 8*Z))/24.D0
      DDSP(1,3,6) = (3*(-1 + Y)*(Y - 2*Z) + X*(6 - 6*Y + 8*Z))/24.D0
      DDSP(1,3,7) = (3*(1 + Y)*(Y + 2*Z) + X*(6 + 6*Y + 8*Z))/24.D0
      DDSP(1,3,8) = (-3*(1 + Y)*(Y + 2*Z) + X*(6 + 6*Y + 8*Z))/24.D0
      DDSP(1,3,9) = -(X*(-3 + 3*Y + 2*Z))/6.D0
      DDSP(1,3,10) = (-3 + 3*Y**2 - 4*X*Z)/12.D0
      DDSP(1,3,11) = (X*(3 + 3*Y - 2*Z))/6.D0
      DDSP(1,3,12) = (3 - 3*Y**2 - 4*X*Z)/12.D0
      DDSP(1,3,13) = (X*(-3 + 3*Y - 2*Z))/6.D0
      DDSP(1,3,14) = (3 - 3*Y**2 - 4*X*Z)/12.D0
      DDSP(1,3,15) = -(X*(3 + 3*Y + 2*Z))/6.D0
      DDSP(1,3,16) = (-3 + 3*Y**2 - 4*X*Z)/12.D0
      DDSP(1,3,17) = -((-3 + 2*X + 3*Y)*Z)/6.D0
      DDSP(1,3,18) = -((3 + 2*X - 3*Y)*Z)/6.D0
      DDSP(1,3,19) = -((3 + 2*X + 3*Y)*Z)/6.D0
      DDSP(1,3,20) = ((3 - 2*X + 3*Y)*Z)/6.D0
      DDSP(1,3,21) = (4*X*Z)/3.D0
      
    !DYDXI      
      DDSP(2,2,1) = (-1 + 2*X**2 + 3*X*(-1 + Z) + Z*(-3 + 2*Z))/12.D0
      DDSP(2,2,2) = (-1 + 3*X + 2*X**2 - 3*(1 + X)*Z + 2*Z**2)/12.D0
      DDSP(2,2,3) = (-1 + 3*X + 2*X**2 - 3*(1 + X)*Z + 2*Z**2)/12.D0
      DDSP(2,2,4) = (-1 + 2*X**2 + 3*X*(-1 + Z) + Z*(-3 + 2*Z))/12.D0
      DDSP(2,2,5) = (-1 + 2*X**2 - 3*X*(1 + Z) + Z*(3 + 2*Z))/12.D0
      DDSP(2,2,6) = (-1 + 2*X**2 + 3*X*(1 + Z) + Z*(3 + 2*Z))/12.D0
      DDSP(2,2,7) = (-1 + 2*X**2 + 3*X*(1 + Z) + Z*(3 + 2*Z))/12.D0
      DDSP(2,2,8) = (-1 + 2*X**2 - 3*X*(1 + Z) + Z*(3 + 2*Z))/12.D0
      DDSP(2,2,9) = (2 - X**2 - Z**2)/6.D0
      DDSP(2,2,10) = (-1 - X*(3 + X) + 3*Z + 3*X*Z - Z**2)/6.D0
      DDSP(2,2,11) = (2 - X**2 - Z**2)/6.D0
      DDSP(2,2,12) = (-1 - X**2 - 3*X*(-1 + Z) - (-3 + Z)*Z)/6.D0
      DDSP(2,2,13) = (2 - X**2 - Z**2)/6.D0
      DDSP(2,2,14) = (-1 - X**2 - 3*X*(1 + Z) - Z*(3 + Z))/6.D0
      DDSP(2,2,15) = (2 - X**2 - Z**2)/6.D0
      DDSP(2,2,16) = (-1 - X**2 + 3*X*(1 + Z) - Z*(3 + Z))/6.D0
      DDSP(2,2,17) = (2 - X**2 - Z**2)/6.D0
      DDSP(2,2,18) = (2 - X**2 - Z**2)/6.D0
      DDSP(2,2,19) = (2 - X**2 - Z**2)/6.D0
      DDSP(2,2,20) = (2 - X**2 - Z**2)/6.D0
      DDSP(2,2,21) = (2*(-2 + X**2 + Z**2))/3.D0
      
      DDSP(2,3,1) = (3*(-1 + X)*(X + 2*Y) + 2*(-3 + 3*X + 4*Y)*Z)/24.D0
      DDSP(2,3,2) = (3*(1 + X)*(X - 2*Y) - 2*(3 + 3*X - 4*Y)*Z)/24.D0
      DDSP(2,3,3) = (-3*(1 + X)*(X + 2*Y) + 2*(3 + 3*X + 4*Y)*Z)/24.D0
      DDSP(2,3,4) = (-3*(-1 + X)*(X - 2*Y) + 2*(3 - 3*X + 4*Y)*Z)/24.D0
      DDSP(2,3,5) = (-3*(-1 + X)*(X + 2*Y) + 2*(-3 + 3*X + 4*Y)*Z)/24.D0
      DDSP(2,3,6) = (-3*(1 + X)*(X - 2*Y) - 2*(3 + 3*X - 4*Y)*Z)/24.D0
      DDSP(2,3,7) = (3*(1 + X)*(X + 2*Y) + 2*(3 + 3*X + 4*Y)*Z)/24.D0
      DDSP(2,3,8) = (3*(-1 + X)*(X - 2*Y) + 2*(3 - 3*X + 4*Y)*Z)/24.D0
      DDSP(2,3,9) = (3 - 3*X**2 - 4*Y*Z)/12.D0
      DDSP(2,3,10) = (Y*(3 + 3*X - 2*Z))/6.D0
      DDSP(2,3,11) = (-3 + 3*X**2 - 4*Y*Z)/12.D0
      DDSP(2,3,12) = -(Y*(-3 + 3*X + 2*Z))/6.D0
      DDSP(2,3,13) = (-3 + 3*X**2 - 4*Y*Z)/12.D0
      DDSP(2,3,14) = -(Y*(3 + 3*X + 2*Z))/6.D0
      DDSP(2,3,15) = (3 - 3*X**2 - 4*Y*Z)/12.D0
      DDSP(2,3,16) = (Y*(-3 + 3*X - 2*Z))/6.D0
      DDSP(2,3,17) = -((-3 + 3*X + 2*Y)*Z)/6.D0
      DDSP(2,3,18) = ((3 + 3*X - 2*Y)*Z)/6.D0
      DDSP(2,3,19) = -((3 + 3*X + 2*Y)*Z)/6.D0
      DDSP(2,3,20) = ((-3 + 3*X - 2*Y)*Z)/6.D0
      DDSP(2,3,21) = (4*Y*Z)/3.D0
      
    !DZDXI
      
      DDSP(3,3,1) = (-1 + 2*X**2 + 3*X*(-1 + Y) + Y*(-3 + 2*Y))/12.D0
      DDSP(3,3,2) = (-1 + 3*X + 2*X**2 - 3*(1 + X)*Y + 2*Y**2)/12.D0
      DDSP(3,3,3) = (-1 + 2*X**2 + 3*X*(1 + Y) + Y*(3 + 2*Y))/12.D0
      DDSP(3,3,4) = (-1 + 2*X**2 - 3*X*(1 + Y) + Y*(3 + 2*Y))/12.D0
      DDSP(3,3,5) = (-1 + 2*X**2 + 3*X*(-1 + Y) + Y*(-3 + 2*Y))/12.D0
      DDSP(3,3,6) = (-1 + 3*X + 2*X**2 - 3*(1 + X)*Y + 2*Y**2)/12.D0
      DDSP(3,3,7) = (-1 + 2*X**2 + 3*X*(1 + Y) + Y*(3 + 2*Y))/12.D0
      DDSP(3,3,8) = (-1 + 2*X**2 - 3*X*(1 + Y) + Y*(3 + 2*Y))/12.D0
      DDSP(3,3,9) = (2 - X**2 - Y**2)/6.D0
      DDSP(3,3,10) = (2 - X**2 - Y**2)/6.D0
      DDSP(3,3,11) = (2 - X**2 - Y**2)/6.D0
      DDSP(3,3,12) = (2 - X**2 - Y**2)/6.D0
      DDSP(3,3,13) = (2 - X**2 - Y**2)/6.D0
      DDSP(3,3,14) = (2 - X**2 - Y**2)/6.D0
      DDSP(3,3,15) = (2 - X**2 - Y**2)/6.D0
      DDSP(3,3,16) = (2 - X**2 - Y**2)/6.D0
      DDSP(3,3,17) = (-1 - X**2 - 3*X*(-1 + Y) - (-3 + Y)*Y)/6.D0
      DDSP(3,3,18) = (-1 - X*(3 + X) + 3*Y + 3*X*Y - Y**2)/6.D0
      DDSP(3,3,19) = (-1 - X**2 - 3*X*(1 + Y) - Y*(3 + Y))/6.D0
      DDSP(3,3,20) = (-1 - X**2 + 3*X*(1 + Y) - Y*(3 + Y))/6.D0
      DDSP(3,3,21) = (2*(-2 + X**2 + Y**2))/3.D0
      
      DO ID=1,21
       DDSP(3,1,ID) = DDSP(1,3,ID)
       DDSP(3,2,ID) = DDSP(2,3,ID)
       DDSP(2,1,ID) = DDSP(1,2,ID)
      ENDDO
      END