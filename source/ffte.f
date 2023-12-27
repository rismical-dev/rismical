C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@cs.tsukuba.ac.jp
C
C
C     3-D COMPLEX FFT ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     CALL ZFFT3D(A,NX,NY,NZ,IOPT)
C
C     A(NX,NY,NZ) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
C     NX IS THE LENGTH OF THE TRANSFORMS IN THE X-DIRECTION (INTEGER*4)
C     NY IS THE LENGTH OF THE TRANSFORMS IN THE Y-DIRECTION (INTEGER*4)
C     NZ IS THE LENGTH OF THE TRANSFORMS IN THE Z-DIRECTION (INTEGER*4)
C       ------------------------------------
C         NX = (2**IP) * (3**IQ) * (5**IR)
C         NY = (2**JP) * (3**JQ) * (5**JR)
C         NZ = (2**KP) * (3**KQ) * (5**KR)
C       ------------------------------------
C     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
C          = -1 FOR FORWARD TRANSFORM
C          = +1 FOR INVERSE TRANSFORM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE ZFFT3D(A,NX,NY,NZ,IOPT)
      IMPLICIT REAL*8 (A-H,O-Z)
c
C     HEADER FILE FOR PARAMETERS
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
C The maximum supported 2-D transform length is 65536.
      PARAMETER (NDA2=65536)
C The maximum supported 3-D transform length is 4096.
      PARAMETER (NDA3=4096)
      PARAMETER (NDA4=256)
C The parameter NBLK is a blocking parameter.
      PARAMETER (NBLK=16)
C      PARAMETER (NBLK=8)  (for PentiumIII and Athlon)
C      PARAMETER (NBLK=16) (for Pentium4, Athlon XP, Opteron, Itanium
C                           and Itanium2)
C The parameter NP is a padding parameter to avoid cache conflicts in
C the FFT routines.
      PARAMETER (NP=8)
C      PARAMETER (NP=2) (for PentiumIII)
C      PARAMETER (NP=4) (for Athlon, Athlon XP, Opteron and Itanium)
C      PARAMETER (NP=8) (for Pentium4 and Itanium2)
C Size of L2 cache
      PARAMETER (L2SIZE=1048576)
c
c
      COMPLEX*16 A(*)
      COMPLEX*16 B((NDA3+NP)*(NBLK+1)+NP)
      COMPLEX*16 WX(NDA3/2+NP),WY(NDA3/2+NP),WZ(NDA3/2+NP)
      DIMENSION LNX(3),LNY(3),LNZ(3)
      SAVE WX,WY,WZ
C
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
      CALL FACTOR(NZ,LNZ)
C
      IF (IOPT .EQ. 0) THEN
        CALL SETTBL(WX,NX)
        CALL SETTBL(WY,NY)
        CALL SETTBL(WZ,NZ)
        RETURN
      END IF
C
      IF (IOPT .EQ. 1) THEN
        DO 10 I=1,NX*NY*NZ
          A(I)=DCONJG(A(I))
   10   CONTINUE
      END IF
C
      NC=(MAX0(NY,NZ)+NP)*NBLK+NP
!$OMP PARALLEL PRIVATE(B)
      CALL ZFFT3D0(A,B,B,B(NC+1),WX,WY,WZ,NX,NY,NZ,LNX,LNY,LNZ)
!$OMP END PARALLEL
C
      IF (IOPT .EQ. 1) THEN
        DN=1.0D0/(DBLE(NX)*DBLE(NY)*DBLE(NZ))
        DO 20 I=1,NX*NY*NZ
          A(I)=DCONJG(A(I))*DN
   20   CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE ZFFT3D0(A,BY,BZ,C,WX,WY,WZ,NX,NY,NZ,LNX,LNY,LNZ)
      IMPLICIT REAL*8 (A-H,O-Z)
c
C     HEADER FILE FOR PARAMETERS
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
C The maximum supported 2-D transform length is 65536.
      PARAMETER (NDA2=65536)
C The maximum supported 3-D transform length is 4096.
      PARAMETER (NDA3=4096)
      PARAMETER (NDA4=256)
C The parameter NBLK is a blocking parameter.
      PARAMETER (NBLK=16)
C      PARAMETER (NBLK=8)  (for PentiumIII and Athlon)
C      PARAMETER (NBLK=16) (for Pentium4, Athlon XP, Opteron, Itanium
C                           and Itanium2)
C The parameter NP is a padding parameter to avoid cache conflicts in
C the FFT routines.
      PARAMETER (NP=8)
C      PARAMETER (NP=2) (for PentiumIII)
C      PARAMETER (NP=4) (for Athlon, Athlon XP, Opteron and Itanium)
C      PARAMETER (NP=8) (for Pentium4 and Itanium2)
C Size of L2 cache
      PARAMETER (L2SIZE=1048576)
c
c
      COMPLEX*16 A(NX,NY,*),BY(NY+NP,*),BZ(NZ+NP,*),C(*)
      COMPLEX*16 WX(*),WY(*),WZ(*)
      DIMENSION LNX(*),LNY(*),LNZ(*)
C
!$OMP DO
      DO 80 J=1,NY
        DO 70 II=1,NX,NBLK
          DO 30 KK=1,NZ,NBLK
            DO 20 I=II,MIN0(II+NBLK-1,NX)
              DO 10 K=KK,MIN0(KK+NBLK-1,NZ)
                BZ(K,I-II+1)=A(I,J,K)
   10         CONTINUE
   20       CONTINUE
   30     CONTINUE
          DO 40 I=II,MIN0(II+NBLK-1,NX)
            CALL FFT235(BZ(1,I-II+1),C,WZ,NZ,LNZ)
   40     CONTINUE
          DO 60 K=1,NZ
            DO 50 I=II,MIN0(II+NBLK-1,NX)
              A(I,J,K)=BZ(K,I-II+1)
   50       CONTINUE
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE
!$OMP DO
      DO 170 K=1,NZ
        DO 150 II=1,NX,NBLK
          DO 110 JJ=1,NY,NBLK
            DO 100 I=II,MIN0(II+NBLK-1,NX)
              DO 90 J=JJ,MIN0(JJ+NBLK-1,NY)
                BY(J,I-II+1)=A(I,J,K)
   90         CONTINUE
  100       CONTINUE
  110     CONTINUE
          DO 120 I=II,MIN0(II+NBLK-1,NX)
            CALL FFT235(BY(1,I-II+1),C,WY,NY,LNY)
  120     CONTINUE
          DO 140 J=1,NY
            DO 130 I=II,MIN0(II+NBLK-1,NX)
              A(I,J,K)=BY(J,I-II+1)
  130       CONTINUE
  140     CONTINUE
  150   CONTINUE
        DO 160 J=1,NY
          CALL FFT235(A(1,J,K),C,WX,NX,LNX)
  160   CONTINUE
  170 CONTINUE
      RETURN
      END
C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@cs.tsukuba.ac.jp
C
C
C     RADIX-2, 3, 4, 5 AND 8 FFT ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE FFT235(A,B,W,N,IP)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
      DIMENSION IP(*)
C
      IF (IP(1) .NE. 1) THEN
        KP4=2-MOD(IP(1)+2,3)
        KP8=(IP(1)-KP4)/3
      ELSE
        KP4=0
        KP8=0
      END IF
C
      KEY=1
      J=1
      L=N
      M=1
      DO 10 K=1,KP8
        L=L/8
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT8(A,B,W(J),M,L)
          ELSE
            CALL FFT8(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT8(A,A,W(J),M,L)
          ELSE
            CALL FFT8(B,A,W(J),M,L)
          END IF
        END IF
        M=M*8
        J=J+L
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT5(A,B,W(J),M,L)
          ELSE
            CALL FFT5(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT5(A,A,W(J),M,L)
          ELSE
            CALL FFT5(B,A,W(J),M,L)
          END IF
        END IF
        M=M*5
        J=J+L
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT4(A,B,W(J),M,L)
          ELSE
            CALL FFT4(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT4(A,A,W(J),M,L)
          ELSE
            CALL FFT4(B,A,W(J),M,L)
          END IF
        END IF
        M=M*4
        J=J+L
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        IF (L .GE. 2) THEN
          IF (KEY .GE. 0) THEN
            CALL FFT3(A,B,W(J),M,L)
          ELSE
            CALL FFT3(B,A,W(J),M,L)
          END IF
          KEY=-KEY
        ELSE
          IF (KEY .GE. 0) THEN
            CALL FFT3(A,A,W(J),M,L)
          ELSE
            CALL FFT3(B,A,W(J),M,L)
          END IF
        END IF
        M=M*3
        J=J+L
   40 CONTINUE
      IF (IP(1) .EQ. 1) THEN
        IF (KEY .GE. 0) THEN
          CALL FFT2(A,A,M)
        ELSE
          CALL FFT2(B,A,M)
        END IF
      END IF
      RETURN
      END
      SUBROUTINE FFT3(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT3A(A,B,W,L)
      ELSE
        CALL FFT3B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT4(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT4A(A,B,W,L)
      ELSE
        CALL FFT4B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT5(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT5A(A,B,W,L)
      ELSE
        CALL FFT5B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE FFT8(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 A(*),B(*),W(*)
C
      IF (M .EQ. 1) THEN
        CALL FFT8A(A,B,W,L)
      ELSE
        CALL FFT8B(A,B,W,M,L)
      END IF
      RETURN
      END
      SUBROUTINE SETTBL(W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 W(*)
      DIMENSION IP(3)
C
      CALL FACTOR(N,IP)
C
      IF (IP(1) .NE. 1) THEN
        KP4=2-MOD(IP(1)+2,3)
        KP8=(IP(1)-KP4)/3
      ELSE
        KP4=0
        KP8=0
      END IF
C
      J=1
      L=N
      DO 10 K=1,KP8
        L=L/8
        CALL SETTBL0(W(J),8,L)
        J=J+L
   10 CONTINUE
      DO 20 K=1,IP(3)
        L=L/5
        CALL SETTBL0(W(J),5,L)
        J=J+L
   20 CONTINUE
      DO 30 K=1,KP4
        L=L/4
        CALL SETTBL0(W(J),4,L)
        J=J+L
   30 CONTINUE
      DO 40 K=1,IP(2)
        L=L/3
        CALL SETTBL0(W(J),3,L)
        J=J+L
   40 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL0(W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(2,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(M)*DBLE(L))
!DIR$ VECTOR ALIGNED
      DO 10 I=1,L
        W(1,I)=DCOS(PX*DBLE(I-1))
        W(2,I)=DSIN(PX*DBLE(I-1))
   10 CONTINUE
      RETURN
      END
      SUBROUTINE SETTBL2(W,N1,N2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(2,N1,*)
C
      PI2=8.0D0*DATAN(1.0D0)
      PX=-PI2/(DBLE(N1)*DBLE(N2))
!$OMP PARALLEL DO
      DO 20 K=1,N2
!DIR$ VECTOR ALIGNED
        DO 10 J=1,N1
          W(1,J,K)=DCOS(PX*DBLE(J-1)*DBLE(K-1))
          W(2,J,K)=DSIN(PX*DBLE(J-1)*DBLE(K-1))
   10   CONTINUE
   20 CONTINUE
!$OMP END PARALLEL DO
      RETURN
      END
      SUBROUTINE FACTOR(N,IP)
      DIMENSION IP(*)
C
      IP(1)=0
      IP(2)=0
      IP(3)=0
      N2=N
      IF (MOD(N,2) .NE. 0 .AND. MOD(N,3) .NE. 0 .AND.
     1    MOD(N,5) .NE. 0) RETURN
   10 IF (N2 .LE. 1) RETURN
      IF (MOD(N2,2) .EQ. 0) THEN
        IP(1)=IP(1)+1
        N2=N2/2
        GO TO 10
      ELSE IF (MOD(N2,3) .EQ. 0) THEN
        IP(2)=IP(2)+1
        N2=N2/3
        GO TO 10
      ELSE IF (MOD(N2,5) .EQ. 0) THEN
        IP(3)=IP(3)+1
        N2=N2/5
        GO TO 10
      END IF
      RETURN
      END
      SUBROUTINE FACTOR8(N,IP)
      DIMENSION IP(*)
      INTEGER*8 N,N2
C
      IP(1)=0
      IP(2)=0
      IP(3)=0
      N2=N
      IF (MOD(N,2) .NE. 0 .AND. MOD(N,3) .NE. 0 .AND.
     1    MOD(N,5) .NE. 0) RETURN
   10 IF (N2 .LE. 1) RETURN
      IF (MOD(N2,2) .EQ. 0) THEN
        IP(1)=IP(1)+1
        N2=N2/2
        GO TO 10
      ELSE IF (MOD(N2,3) .EQ. 0) THEN
        IP(2)=IP(2)+1
        N2=N2/3
        GO TO 10
      ELSE IF (MOD(N2,5) .EQ. 0) THEN
        IP(3)=IP(3)+1
        N2=N2/5
        GO TO 10
      END IF
      RETURN
      END
C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         GRADUATE SCHOOL OF SYSTEMS AND INFORMATION ENGINEERING
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@cs.tsukuba.ac.jp
C
C
C     RADIX-2, 3, 4, 5 AND 8 FFT KERNEL ROUTINE
C
C     FORTRAN77 SOURCE PROGRAM
C
C     WRITTEN BY DAISUKE TAKAHASHI
C
      SUBROUTINE FFT2(A,B,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,*),B(2,M,*)
C
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        X0=A(1,I,1)
        Y0=A(2,I,1)
        X1=A(1,I,2)
        Y1=A(2,I,2)
        B(1,I,1)=X0+X1
        B(2,I,1)=Y0+Y1
        B(1,I,2)=X0-X1
        B(2,I,2)=Y0-Y1
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT3A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,3,*),W(2,*)
      DATA C31/0.86602540378443865D0/C32/0.5D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        X0=A(1,J,2)+A(1,J,3)
        Y0=A(2,J,2)+A(2,J,3)
        X1=A(1,J,1)-C32*X0
        Y1=A(2,J,1)-C32*Y0
        X2=C31*(A(2,J,2)-A(2,J,3))
        Y2=C31*(A(1,J,3)-A(1,J,2))
        B(1,1,J)=A(1,J,1)+X0
        B(2,1,J)=A(2,J,1)+Y0
        B(1,2,J)=WR1*(X1+X2)-WI1*(Y1+Y2)
        B(2,2,J)=WR1*(Y1+Y2)+WI1*(X1+X2)
        B(1,3,J)=WR2*(X1-X2)-WI2*(Y1-Y2)
        B(2,3,J)=WR2*(Y1-Y2)+WI2*(X1-X2)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT3B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,3,*),W(2,*)
      DATA C31/0.86602540378443865D0/C32/0.5D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        X0=A(1,I,1,2)+A(1,I,1,3)
        Y0=A(2,I,1,2)+A(2,I,1,3)
        X1=A(1,I,1,1)-C32*X0
        Y1=A(2,I,1,1)-C32*Y0
        X2=C31*(A(2,I,1,2)-A(2,I,1,3))
        Y2=C31*(A(1,I,1,3)-A(1,I,1,2))
        B(1,I,1,1)=A(1,I,1,1)+X0
        B(2,I,1,1)=A(2,I,1,1)+Y0
        B(1,I,2,1)=X1+X2
        B(2,I,2,1)=Y1+Y2
        B(1,I,3,1)=X1-X2
        B(2,I,3,1)=Y1-Y2
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          X0=A(1,I,J,2)+A(1,I,J,3)
          Y0=A(2,I,J,2)+A(2,I,J,3)
          X1=A(1,I,J,1)-C32*X0
          Y1=A(2,I,J,1)-C32*Y0
          X2=C31*(A(2,I,J,2)-A(2,I,J,3))
          Y2=C31*(A(1,I,J,3)-A(1,I,J,2))
          B(1,I,1,J)=A(1,I,J,1)+X0
          B(2,I,1,J)=A(2,I,J,1)+Y0
          B(1,I,2,J)=WR1*(X1+X2)-WI1*(Y1+Y2)
          B(2,I,2,J)=WR1*(Y1+Y2)+WI1*(X1+X2)
          B(1,I,3,J)=WR2*(X1-X2)-WI2*(Y1-Y2)
          B(2,I,3,J)=WR2*(Y1-Y2)+WI2*(X1-X2)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT4A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,4,*),W(2,*)
C
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        X0=A(1,J,1)+A(1,J,3)
        Y0=A(2,J,1)+A(2,J,3)
        X1=A(1,J,1)-A(1,J,3)
        Y1=A(2,J,1)-A(2,J,3)
        X2=A(1,J,2)+A(1,J,4)
        Y2=A(2,J,2)+A(2,J,4)
        X3=A(2,J,2)-A(2,J,4)
        Y3=A(1,J,4)-A(1,J,2)
        B(1,1,J)=X0+X2
        B(2,1,J)=Y0+Y2
        B(1,3,J)=WR2*(X0-X2)-WI2*(Y0-Y2)
        B(2,3,J)=WR2*(Y0-Y2)+WI2*(X0-X2)
        B(1,2,J)=WR1*(X1+X3)-WI1*(Y1+Y3)
        B(2,2,J)=WR1*(Y1+Y3)+WI1*(X1+X3)
        B(1,4,J)=WR3*(X1-X3)-WI3*(Y1-Y3)
        B(2,4,J)=WR3*(Y1-Y3)+WI3*(X1-X3)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT4B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,4,*),W(2,*)
C
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        X0=A(1,I,1,1)+A(1,I,1,3)
        Y0=A(2,I,1,1)+A(2,I,1,3)
        X1=A(1,I,1,1)-A(1,I,1,3)
        Y1=A(2,I,1,1)-A(2,I,1,3)
        X2=A(1,I,1,2)+A(1,I,1,4)
        Y2=A(2,I,1,2)+A(2,I,1,4)
        X3=A(2,I,1,2)-A(2,I,1,4)
        Y3=A(1,I,1,4)-A(1,I,1,2)
        B(1,I,1,1)=X0+X2
        B(2,I,1,1)=Y0+Y2
        B(1,I,3,1)=X0-X2
        B(2,I,3,1)=Y0-Y2
        B(1,I,2,1)=X1+X3
        B(2,I,2,1)=Y1+Y3
        B(1,I,4,1)=X1-X3
        B(2,I,4,1)=Y1-Y3
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          X0=A(1,I,J,1)+A(1,I,J,3)
          Y0=A(2,I,J,1)+A(2,I,J,3)
          X1=A(1,I,J,1)-A(1,I,J,3)
          Y1=A(2,I,J,1)-A(2,I,J,3)
          X2=A(1,I,J,2)+A(1,I,J,4)
          Y2=A(2,I,J,2)+A(2,I,J,4)
          X3=A(2,I,J,2)-A(2,I,J,4)
          Y3=A(1,I,J,4)-A(1,I,J,2)
          B(1,I,1,J)=X0+X2
          B(2,I,1,J)=Y0+Y2
          B(1,I,3,J)=WR2*(X0-X2)-WI2*(Y0-Y2)
          B(2,I,3,J)=WR2*(Y0-Y2)+WI2*(X0-X2)
          B(1,I,2,J)=WR1*(X1+X3)-WI1*(Y1+Y3)
          B(2,I,2,J)=WR1*(Y1+Y3)+WI1*(X1+X3)
          B(1,I,4,J)=WR3*(X1-X3)-WI3*(Y1-Y3)
          B(2,I,4,J)=WR3*(Y1-Y3)+WI3*(X1-X3)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT5A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,5,*),W(2,*)
      DATA C51/0.95105651629515357D0/C52/0.61803398874989485D0/
     1     C53/0.55901699437494742D0/C54/0.25D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
        X0=A(1,J,2)+A(1,J,5)
        Y0=A(2,J,2)+A(2,J,5)
        X1=A(1,J,3)+A(1,J,4)
        Y1=A(2,J,3)+A(2,J,4)
        X2=C51*(A(1,J,2)-A(1,J,5))
        Y2=C51*(A(2,J,2)-A(2,J,5))
        X3=C51*(A(1,J,3)-A(1,J,4))
        Y3=C51*(A(2,J,3)-A(2,J,4))
        X4=X0+X1
        Y4=Y0+Y1
        X5=C53*(X0-X1)
        Y5=C53*(Y0-Y1)
        X6=A(1,J,1)-C54*X4
        Y6=A(2,J,1)-C54*Y4
        X7=X6+X5
        Y7=Y6+Y5
        X8=X6-X5
        Y8=Y6-Y5
        X9=Y2+C52*Y3
        Y9=-X2-C52*X3
        X10=C52*Y2-Y3
        Y10=X3-C52*X2
        B(1,1,J)=A(1,J,1)+X4
        B(2,1,J)=A(2,J,1)+Y4
        B(1,2,J)=WR1*(X7+X9)-WI1*(Y7+Y9)
        B(2,2,J)=WR1*(Y7+Y9)+WI1*(X7+X9)
        B(1,3,J)=WR2*(X8+X10)-WI2*(Y8+Y10)
        B(2,3,J)=WR2*(Y8+Y10)+WI2*(X8+X10)
        B(1,4,J)=WR3*(X8-X10)-WI3*(Y8-Y10)
        B(2,4,J)=WR3*(Y8-Y10)+WI3*(X8-X10)
        B(1,5,J)=WR4*(X7-X9)-WI4*(Y7-Y9)
        B(2,5,J)=WR4*(Y7-Y9)+WI4*(X7-X9)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT5B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,5,*),W(2,*)
      DATA C51/0.95105651629515357D0/C52/0.61803398874989485D0/
     1     C53/0.55901699437494742D0/C54/0.25D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        X0=A(1,I,1,2)+A(1,I,1,5)
        Y0=A(2,I,1,2)+A(2,I,1,5)
        X1=A(1,I,1,3)+A(1,I,1,4)
        Y1=A(2,I,1,3)+A(2,I,1,4)
        X2=C51*(A(1,I,1,2)-A(1,I,1,5))
        Y2=C51*(A(2,I,1,2)-A(2,I,1,5))
        X3=C51*(A(1,I,1,3)-A(1,I,1,4))
        Y3=C51*(A(2,I,1,3)-A(2,I,1,4))
        X4=X0+X1
        Y4=Y0+Y1
        X5=C53*(X0-X1)
        Y5=C53*(Y0-Y1)
        X6=A(1,I,1,1)-C54*X4
        Y6=A(2,I,1,1)-C54*Y4
        X7=X6+X5
        Y7=Y6+Y5
        X8=X6-X5
        Y8=Y6-Y5
        X9=Y2+C52*Y3
        Y9=-X2-C52*X3
        X10=C52*Y2-Y3
        Y10=X3-C52*X2
        B(1,I,1,1)=A(1,I,1,1)+X4
        B(2,I,1,1)=A(2,I,1,1)+Y4
        B(1,I,2,1)=X7+X9
        B(2,I,2,1)=Y7+Y9
        B(1,I,3,1)=X8+X10
        B(2,I,3,1)=Y8+Y10
        B(1,I,4,1)=X8-X10
        B(2,I,4,1)=Y8-Y10
        B(1,I,5,1)=X7-X9
        B(2,I,5,1)=Y7-Y9
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          X0=A(1,I,J,2)+A(1,I,J,5)
          Y0=A(2,I,J,2)+A(2,I,J,5)
          X1=A(1,I,J,3)+A(1,I,J,4)
          Y1=A(2,I,J,3)+A(2,I,J,4)
          X2=C51*(A(1,I,J,2)-A(1,I,J,5))
          Y2=C51*(A(2,I,J,2)-A(2,I,J,5))
          X3=C51*(A(1,I,J,3)-A(1,I,J,4))
          Y3=C51*(A(2,I,J,3)-A(2,I,J,4))
          X4=X0+X1
          Y4=Y0+Y1
          X5=C53*(X0-X1)
          Y5=C53*(Y0-Y1)
          X6=A(1,I,J,1)-C54*X4
          Y6=A(2,I,J,1)-C54*Y4
          X7=X6+X5
          Y7=Y6+Y5
          X8=X6-X5
          Y8=Y6-Y5
          X9=Y2+C52*Y3
          Y9=-X2-C52*X3
          X10=C52*Y2-Y3
          Y10=X3-C52*X2
          B(1,I,1,J)=A(1,I,J,1)+X4
          B(2,I,1,J)=A(2,I,J,1)+Y4
          B(1,I,2,J)=WR1*(X7+X9)-WI1*(Y7+Y9)
          B(2,I,2,J)=WR1*(Y7+Y9)+WI1*(X7+X9)
          B(1,I,3,J)=WR2*(X8+X10)-WI2*(Y8+Y10)
          B(2,I,3,J)=WR2*(Y8+Y10)+WI2*(X8+X10)
          B(1,I,4,J)=WR3*(X8-X10)-WI3*(Y8-Y10)
          B(2,I,4,J)=WR3*(Y8-Y10)+WI3*(X8-X10)
          B(1,I,5,J)=WR4*(X7-X9)-WI4*(Y7-Y9)
          B(2,I,5,J)=WR4*(Y7-Y9)+WI4*(X7-X9)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
      SUBROUTINE FFT8A(A,B,W,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,L,*),B(2,8,*),W(2,*)
      DATA C81/0.70710678118654752D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 J=1,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
        WR5=WR2*WR3-WI2*WI3
        WI5=WR2*WI3+WI2*WR3
        WR6=WR3*WR3-WI3*WI3
        WI6=WR3*WI3+WR3*WI3
        WR7=WR3*WR4-WI3*WI4
        WI7=WR3*WI4+WI3*WR4
        X0=A(1,J,1)+A(1,J,5)
        Y0=A(2,J,1)+A(2,J,5)
        X1=A(1,J,1)-A(1,J,5)
        Y1=A(2,J,1)-A(2,J,5)
        X2=A(1,J,3)+A(1,J,7)
        Y2=A(2,J,3)+A(2,J,7)
        X3=A(2,J,3)-A(2,J,7)
        Y3=A(1,J,7)-A(1,J,3)
        U0=X0+X2
        V0=Y0+Y2
        U1=X0-X2
        V1=Y0-Y2
        X4=A(1,J,2)+A(1,J,6)
        Y4=A(2,J,2)+A(2,J,6)
        X5=A(1,J,2)-A(1,J,6)
        Y5=A(2,J,2)-A(2,J,6)
        X6=A(1,J,4)+A(1,J,8)
        Y6=A(2,J,4)+A(2,J,8)
        X7=A(1,J,4)-A(1,J,8)
        Y7=A(2,J,4)-A(2,J,8)
        U2=X4+X6
        V2=Y4+Y6
        U3=Y4-Y6
        V3=X6-X4
        B(1,1,J)=U0+U2
        B(2,1,J)=V0+V2
        B(1,5,J)=WR4*(U0-U2)-WI4*(V0-V2)
        B(2,5,J)=WR4*(V0-V2)+WI4*(U0-U2)
        B(1,3,J)=WR2*(U1+U3)-WI2*(V1+V3)
        B(2,3,J)=WR2*(V1+V3)+WI2*(U1+U3)
        B(1,7,J)=WR6*(U1-U3)-WI6*(V1-V3)
        B(2,7,J)=WR6*(V1-V3)+WI6*(U1-U3)
        U0=X1+C81*(X5-X7)
        V0=Y1+C81*(Y5-Y7)
        U1=X1-C81*(X5-X7)
        V1=Y1-C81*(Y5-Y7)
        U2=X3+C81*(Y5+Y7)
        V2=Y3-C81*(X5+X7)
        U3=X3-C81*(Y5+Y7)
        V3=Y3+C81*(X5+X7)
        B(1,2,J)=WR1*(U0+U2)-WI1*(V0+V2)
        B(2,2,J)=WR1*(V0+V2)+WI1*(U0+U2)
        B(1,6,J)=WR5*(U1+U3)-WI5*(V1+V3)
        B(2,6,J)=WR5*(V1+V3)+WI5*(U1+U3)
        B(1,4,J)=WR3*(U1-U3)-WI3*(V1-V3)
        B(2,4,J)=WR3*(V1-V3)+WI3*(U1-U3)
        B(1,8,J)=WR7*(U0-U2)-WI7*(V0-V2)
        B(2,8,J)=WR7*(V0-V2)+WI7*(U0-U2)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FFT8B(A,B,W,M,L)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(2,M,L,*),B(2,M,8,*),W(2,*)
      DATA C81/0.70710678118654752D0/
C
!DIR$ VECTOR ALIGNED
      DO 10 I=1,M
        X0=A(1,I,1,1)+A(1,I,1,5)
        Y0=A(2,I,1,1)+A(2,I,1,5)
        X1=A(1,I,1,1)-A(1,I,1,5)
        Y1=A(2,I,1,1)-A(2,I,1,5)
        X2=A(1,I,1,3)+A(1,I,1,7)
        Y2=A(2,I,1,3)+A(2,I,1,7)
        X3=A(2,I,1,3)-A(2,I,1,7)
        Y3=A(1,I,1,7)-A(1,I,1,3)
        U0=X0+X2
        V0=Y0+Y2
        U1=X0-X2
        V1=Y0-Y2
        X4=A(1,I,1,2)+A(1,I,1,6)
        Y4=A(2,I,1,2)+A(2,I,1,6)
        X5=A(1,I,1,2)-A(1,I,1,6)
        Y5=A(2,I,1,2)-A(2,I,1,6)
        X6=A(1,I,1,4)+A(1,I,1,8)
        Y6=A(2,I,1,4)+A(2,I,1,8)
        X7=A(1,I,1,4)-A(1,I,1,8)
        Y7=A(2,I,1,4)-A(2,I,1,8)
        U2=X4+X6
        V2=Y4+Y6
        U3=Y4-Y6
        V3=X6-X4
        B(1,I,1,1)=U0+U2
        B(2,I,1,1)=V0+V2
        B(1,I,5,1)=U0-U2
        B(2,I,5,1)=V0-V2
        B(1,I,3,1)=U1+U3
        B(2,I,3,1)=V1+V3
        B(1,I,7,1)=U1-U3
        B(2,I,7,1)=V1-V3
        U0=X1+C81*(X5-X7)
        V0=Y1+C81*(Y5-Y7)
        U1=X1-C81*(X5-X7)
        V1=Y1-C81*(Y5-Y7)
        U2=X3+C81*(Y5+Y7)
        V2=Y3-C81*(X5+X7)
        U3=X3-C81*(Y5+Y7)
        V3=Y3+C81*(X5+X7)
        B(1,I,2,1)=U0+U2
        B(2,I,2,1)=V0+V2
        B(1,I,6,1)=U1+U3
        B(2,I,6,1)=V1+V3
        B(1,I,4,1)=U1-U3
        B(2,I,4,1)=V1-V3
        B(1,I,8,1)=U0-U2
        B(2,I,8,1)=V0-V2
   10 CONTINUE
      DO 30 J=2,L
        WR1=W(1,J)
        WI1=W(2,J)
        WR2=WR1*WR1-WI1*WI1
        WI2=WR1*WI1+WR1*WI1
        WR3=WR1*WR2-WI1*WI2
        WI3=WR1*WI2+WI1*WR2
        WR4=WR2*WR2-WI2*WI2
        WI4=WR2*WI2+WR2*WI2
        WR5=WR2*WR3-WI2*WI3
        WI5=WR2*WI3+WI2*WR3
        WR6=WR3*WR3-WI3*WI3
        WI6=WR3*WI3+WR3*WI3
        WR7=WR3*WR4-WI3*WI4
        WI7=WR3*WI4+WI3*WR4
!DIR$ VECTOR ALIGNED
        DO 20 I=1,M
          X0=A(1,I,J,1)+A(1,I,J,5)
          Y0=A(2,I,J,1)+A(2,I,J,5)
          X1=A(1,I,J,1)-A(1,I,J,5)
          Y1=A(2,I,J,1)-A(2,I,J,5)
          X2=A(1,I,J,3)+A(1,I,J,7)
          Y2=A(2,I,J,3)+A(2,I,J,7)
          X3=A(2,I,J,3)-A(2,I,J,7)
          Y3=A(1,I,J,7)-A(1,I,J,3)
          U0=X0+X2
          V0=Y0+Y2
          U1=X0-X2
          V1=Y0-Y2
          X4=A(1,I,J,2)+A(1,I,J,6)
          Y4=A(2,I,J,2)+A(2,I,J,6)
          X5=A(1,I,J,2)-A(1,I,J,6)
          Y5=A(2,I,J,2)-A(2,I,J,6)
          X6=A(1,I,J,4)+A(1,I,J,8)
          Y6=A(2,I,J,4)+A(2,I,J,8)
          X7=A(1,I,J,4)-A(1,I,J,8)
          Y7=A(2,I,J,4)-A(2,I,J,8)
          U2=X4+X6
          V2=Y4+Y6
          U3=Y4-Y6
          V3=X6-X4
          B(1,I,1,J)=U0+U2
          B(2,I,1,J)=V0+V2
          B(1,I,5,J)=WR4*(U0-U2)-WI4*(V0-V2)
          B(2,I,5,J)=WR4*(V0-V2)+WI4*(U0-U2)
          B(1,I,3,J)=WR2*(U1+U3)-WI2*(V1+V3)
          B(2,I,3,J)=WR2*(V1+V3)+WI2*(U1+U3)
          B(1,I,7,J)=WR6*(U1-U3)-WI6*(V1-V3)
          B(2,I,7,J)=WR6*(V1-V3)+WI6*(U1-U3)
          U0=X1+C81*(X5-X7)
          V0=Y1+C81*(Y5-Y7)
          U1=X1-C81*(X5-X7)
          V1=Y1-C81*(Y5-Y7)
          U2=X3+C81*(Y5+Y7)
          V2=Y3-C81*(X5+X7)
          U3=X3-C81*(Y5+Y7)
          V3=Y3+C81*(X5+X7)
          B(1,I,2,J)=WR1*(U0+U2)-WI1*(V0+V2)
          B(2,I,2,J)=WR1*(V0+V2)+WI1*(U0+U2)
          B(1,I,6,J)=WR5*(U1+U3)-WI5*(V1+V3)
          B(2,I,6,J)=WR5*(V1+V3)+WI5*(U1+U3)
          B(1,I,4,J)=WR3*(U1-U3)-WI3*(V1-V3)
          B(2,I,4,J)=WR3*(V1-V3)+WI3*(U1-U3)
          B(1,I,8,J)=WR7*(U0-U2)-WI7*(V0-V2)
          B(2,I,8,J)=WR7*(V0-V2)+WI7*(U0-U2)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END

