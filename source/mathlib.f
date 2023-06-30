c---------------------------------------------------------
c     1-Dimensional Fourier Transform
c---------------------------------------------------------
      subroutine fft1d(ngrid,rdelta,func,ix)
c
c     ngrid    ... number of grid of r-space
c     rdelta   ... grid width of r-space [Angstrom]
c     func     ... Function are transformed
c     ix       ... flag for transform ix>0 Normal FFT
c                                     ix<0 Reverse FFT
c
c     for ix>0
c     F(k)=4*pi*integral(0 to infinity)[sin(kr)*r/k*f(r)]dr
c
c     for ix<0
c     f(r)=1/(2*pi**2)*integral(0 to infinity)[sin(kr)*k/r*F(k)]dk
c
      implicit real*8 (a-h,o-z)

      dimension func(ngrid),dum(ngrid)
c$$$      dimension wsave(4*ngrid+15)

      pi=dacos(-1.d0)
      deltak=pi/(dble(ngrid)*rdelta)
c---------------------------------------------------------
c$$$      call sinti(ngrid,wsave)
c
c     --- for ix>0  (Normal FFT)
c
      if (ix.gt.0) then
         dum(1)=0.d0
         do i=1,ngrid-1
            rr=rdelta*dble(i)
            dum(i+1)=rr*func(i)
         enddo

         call sinft(dum,ngrid)
c$$$         call sint(ngrid,dum,wsave)

         do j=1,ngrid-1
            rk=deltak*dble(j)
            func(j)=dum(j+1)*4.d0*pi*rdelta/rk
         enddo
         func(ngrid)=0.d0

         return
      endif
c     
c     --- for ix<0  (Reverse FFT)
c     
      if (ix.lt.0) then
         dum(1)=0.d0
         do j=1,ngrid-1
            rk=deltak*dble(j)
            dum(j+1)=rk*func(j)
         enddo
         
         call sinft(dum,ngrid)
c$$$         call sint(ngrid,dum,wsave)
         
         do i=1,ngrid-1
            rr=rdelta*dble(i)
            func(i)=dum(i+1)*deltak/(2.d0*pi*pi)/rr
         enddo
         func(ngrid)=0.d0

         return
      endif

c---------------------------------------------------------
      return
      end
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
c------+NOR>
C
C
C
c---------------------------------------------------------
c     Product of matrix
c---------------------------------------------------------
      subroutine matprd(array1,array2,result,n1,n2,n3,nx1,nx2,nx3,ill)
c
c     array1 x array2 = result 
c     (n1,n2)  (n2,n3)  (n1,n3)
c     
c     nx1,nx2,nx3 ... size of array
c
      implicit real*8 (a-h,o-z)

      dimension result(nx1,nx3)
      dimension array1(nx1,nx2)
      dimension array2(nx2,nx3)
c---------------------------------------------------------
      ill=0
      if (nx1.lt.n1) ill=30000
      if (nx2.lt.n2) ill=30000
      if (nx3.lt.n3) ill=30000
c---------------------------------------------------------
      do i=1,n1
         do j=1,n3
            sum=0.d0
            do k=1,n2
               sum=sum+array1(i,k)*array2(k,j)
            enddo
            result(i,j)=sum
         enddo
      enddo
c---------------------------------------------------------
      return
      end
c---------------------------------------------------------
c     Heaviside Function 
c---------------------------------------------------------
      real*8 function dheaviside(a)
      real*8 a
c---------------------------------------------------------
      dheaviside=0.d0
      if (a.gt.1.d0) dheaviside=1.d0
c---------------------------------------------------------
      return
      end
c**************************************************************
      subroutine matdet(a,n,nmax,det)
c     ------------------------------------------------------
c     Calculate Determinant
c     ------------------------------------------------------
c     a       ... Matrix
c     n,nmax  ... Size of Matrix, Size of Array Declaration
c     det     ... Determinant (OUTPUT)
c
      implicit real*8 (a-h,o-z)

      dimension a(nmax,nmax),indx(n)
c     ------------------------------------------------------
      call dgetrf(n,n,a,nmax,indx,info)

      det=1
      do i=1,n
         det=det*a(i,i)
      enddo
      detp=1
      do i=1,n
         if (indx(i).ne.i) then
            detp=-detp
         endif
      enddo
      det=det*detp
c     ------------------------------------------------------
      return
      end
c**************************************************************
      subroutine matinv(a,n,nmax)
c     ------------------------------------------------------
c     Inverse of a Matrix
c     ------------------------------------------------------
c     a       ... Matrix (IN/OUT)
c     n,nmax  ... Size of Matrix, Size of Array Declaration
c
      implicit real*8(a-h,o-z)
      
      dimension a(nmax,nmax),indx(nmax),work(nmax*64)
c     ------------------------------------------------------
      call dgetrf(n,n,a,nmax,indx,info)
      nwork=nmax*64
      call dgetri(n,a,nmax,indx,work,nwork,info)
c     ------------------------------------------------------
      return
      end
c**************************************************************
      subroutine matsolv(a,x,n,nmax,nx)
c     -------------------------------------------------------
c     Solve Linear Equation by using LU decomposition
c     -------------------------------------------------------
c     a       ... Matrix (IN)
c     x       ... Matrix (IN/OUT)
c     nmax    ... Size of Array Declaration
c
      implicit real*8(a-h,o-z)
      character*1 trans
      dimension a(nmax,nmax),x(nmax,nx),indx(nmax)
c     -------------------------------------------------------
      call dgetrf(n,n,a,nmax,indx,info)
      trans="N"
      call dgetrs(trans,n,nx,a,nmax,indx,x,nmax,info)
c     -------------------------------------------------------
      return
      end
c**************************************************************
c-------------------------------------------------------------
c     Rotate molecular frame to dipole princepal axis
c-------------------------------------------------------------
      subroutine rotate2zaxis(rv,q2,n2)
      
      implicit real*8 (a-h,o-z)
      LOGICAL debug

      include "solvent.i"

      dimension q2(n2),rv(3,n2)
      dimension qwork(100),rwork(3,100)
c-------------------------------------------------------------
c
c     count number of species
c
      ispc=0
      jspc=0
      do i=1,n2
         if (jspc.ne.nspc(i)) ispc=ispc+1
         jspc=nspc(i)
      enddo
c
      do i=1,ispc
         natom=0
         do j=1,n2
            if (nspc(j).eq.i) then
               natom=natom+1
               qwork(natom)=q2(j)
               rwork(1,natom)=rv(1,j)
               rwork(2,natom)=rv(2,j)
               rwork(3,natom)=rv(3,j)
            endif
         enddo

         call diprot(qwork,rwork,natom,100)

         natom=0
         do j=1,n2
            if (nspc(j).eq.i) then
               natom=natom+1
               rv(1,j)=rwork(1,natom)
               rv(2,j)=rwork(2,natom)
               rv(3,j)=rwork(3,natom)
            endif
         enddo

      enddo
c
c     debug write
c
      debug=.false.
      if (debug) then
         write(*,9998)
         do i=1,n2
            write(*,9999) i,nspc(i),rv(1,i),rv(2,i),rv(3,i)
         enddo
      endif
c-------------------------------------------------------------
      return
 9998 format (2x,"ATOM",2x,"SPC",8x,"X",15x,"Y",15x,"Z")
 9999 format (2x,i4,i5,3f16.8)
      end
c**************************************************************
      subroutine differential(wr,rdelta,nlong,nsize)
c     ------------------------------------------------------
c     Differentiation for functions 
c     ------------------------------------------------------
c     
c     INPUT ---------------
c     wr      -- numerical function
c     rdelta  -- width of grid
c     nlong   -- size of data of wr
c     nsize   -- size of array of wr
c     
c     TEMP ----------------
c     dr      -- dummy 
c     ax      -- numerical integration data
c
c     OUTPUT --------------
c     wr      -- derivative of input wr
c
c     REFERENCE -----------
c         HANDBOOK OF MATHEMATICAL FUNCTIONS  P-914
c
c         *** Coefficients For Differentiation
c          data ((ax(i,j),j=1,5),i=1,5)
c         $     /-50.d0, 96.d0,-72.d0, 32.d0, -6.d0,
c         $       -6.d0,-20.d0, 36.d0,-12.d0,  2.d0,
c         $        2.d0,-16.d0,  0.d0, 16.d0, -2.d0,
c         $       -2.d0, 12.d0,-36.d0, 20.d0,  6.d0,
c         $        6.d0,-32.d0, 72.d0,-96.d0, 50.d0/
c               
c     ------------------------------------------------------
      implicit real*8 (a-h,o-z)

      dimension wr(nsize),ax(5,5),dr(nlong)
c     ------------------------------------------------------
      ax(1,1)=-50.d0
      ax(1,2)=96.d0
      ax(1,3)=-72.d0
      ax(1,4)=32.d0
      ax(1,5)=-6.d0
      ax(2,1)=-6.d0
      ax(2,2)=-20.d0
      ax(2,3)=36.d0
      ax(2,4)=-12.d0
      ax(2,5)=2.d0
      ax(3,1)=2.d0
      ax(3,2)=-16.d0
      ax(3,3)=0.d0
      ax(3,4)=16.d0
      ax(3,5)=-2.d0
      ax(4,1)=-2.d0
      ax(4,2)=12.d0
      ax(4,3)=-36.d0
      ax(4,4)=20.d0
      ax(4,5)=6.d0
      ax(5,1)=6.d0
      ax(5,2)=-32.d0
      ax(5,3)=72.d0
      ax(5,4)=-96.d0
      ax(5,5)=50.d0
c     ------------------------------------------------------
      i=1
      dr(i)=(ax(1,1)*wr(i)
     &      +ax(1,2)*wr(i+1)
     &      +ax(1,3)*wr(i+2)
     &      +ax(1,4)*wr(i+3)
     &      +ax(1,5)*wr(i+4))
     &     /(24.d0*rdelta)

      i=2
      dr(i)=(ax(2,1)*wr(i-1)
     &      +ax(2,2)*wr(i)
     &      +ax(2,3)*wr(i+1)
     &      +ax(2,4)*wr(i+2)
     &      +ax(2,5)*wr(i+3))
     &     /(24.d0*rdelta)

      do i=3,nlong-2
         dr(i)=(ax(3,1)*wr(i-2)
     &         +ax(3,2)*wr(i-1)
     &         +ax(3,3)*wr(i)
     &         +ax(3,4)*wr(i+1)
     &         +ax(3,5)*wr(i+2))
     &      /(24.d0*rdelta)
      enddo

      i=nlong-1
      dr(i)=(ax(4,1)*wr(i-3)
     &      +ax(4,2)*wr(i-2)
     &      +ax(4,3)*wr(i-1)
     &      +ax(4,4)*wr(i)
     &      +ax(4,5)*wr(i+1))
     &     /(24.d0*rdelta)
      
      i=nlong
      dr(i)=(ax(5,1)*wr(i-4)
     &      +ax(5,2)*wr(i-3)
     &      +ax(5,3)*wr(i-2)
     &      +ax(5,4)*wr(i-1)
     &      +ax(5,5)*wr(i))
     &     /(24.d0*rdelta)

      do i=1,nlong
         wr(i)=dr(i)
      enddo
c     -----------------------------------------------
      return

      end
c**************************************************************
      subroutine zdifferential(zwr,rdelta,nlong,nsize)
c     ------------------------------------------------------
c     Differentiation for functions 
c     ------------------------------------------------------
c     
c     INPUT ---------------
c     wr      -- numerical function
c     rdelta  -- width of grid
c     nlong   -- size of data of wr
c     nsize   -- size of array of wr
c     
c     TEMP ----------------
c     dr      -- dummy 
c     ax      -- numerical integration data
c
c     OUTPUT --------------
c     wr      -- derivative of input wr
c
c     REFERENCE -----------
c         HANDBOOK OF MATHEMATICAL FUNCTIONS  P-914
c
c         *** Coefficients For Differentiation
c          data ((ax(i,j),j=1,5),i=1,5)
c         $     /-50.d0, 96.d0,-72.d0, 32.d0, -6.d0,
c         $       -6.d0,-20.d0, 36.d0,-12.d0,  2.d0,
c         $        2.d0,-16.d0,  0.d0, 16.d0, -2.d0,
c         $       -2.d0, 12.d0,-36.d0, 20.d0,  6.d0,
c         $        6.d0,-32.d0, 72.d0,-96.d0, 50.d0/
c               
c     ------------------------------------------------------
      implicit real*8 (a-h,o-y)
      implicit complex*16 (z)

      dimension zwr(nsize),ax(5,5),zdr(nlong)
c     ------------------------------------------------------
      ax(1,1)=-50.d0
      ax(1,2)=96.d0
      ax(1,3)=-72.d0
      ax(1,4)=32.d0
      ax(1,5)=-6.d0
      ax(2,1)=-6.d0
      ax(2,2)=-20.d0
      ax(2,3)=36.d0
      ax(2,4)=-12.d0
      ax(2,5)=2.d0
      ax(3,1)=2.d0
      ax(3,2)=-16.d0
      ax(3,3)=0.d0
      ax(3,4)=16.d0
      ax(3,5)=-2.d0
      ax(4,1)=-2.d0
      ax(4,2)=12.d0
      ax(4,3)=-36.d0
      ax(4,4)=20.d0
      ax(4,5)=6.d0
      ax(5,1)=6.d0
      ax(5,2)=-32.d0
      ax(5,3)=72.d0
      ax(5,4)=-96.d0
      ax(5,5)=50.d0
c     ------------------------------------------------------
      i=1
      zdr(i)=(ax(1,1)*zwr(i)
     &      +ax(1,2)*zwr(i+1)
     &      +ax(1,3)*zwr(i+2)
     &      +ax(1,4)*zwr(i+3)
     &      +ax(1,5)*zwr(i+4))
     &     /(24.d0*rdelta)

      i=2
      zdr(i)=(ax(2,1)*zwr(i-1)
     &      +ax(2,2)*zwr(i)
     &      +ax(2,3)*zwr(i+1)
     &      +ax(2,4)*zwr(i+2)
     &      +ax(2,5)*zwr(i+3))
     &     /(24.d0*rdelta)

      do i=3,nlong-2
         zdr(i)=(ax(3,1)*zwr(i-2)
     &          +ax(3,2)*zwr(i-1)
     &          +ax(3,3)*zwr(i)
     &          +ax(3,4)*zwr(i+1)
     &          +ax(3,5)*zwr(i+2))
     &      /(24.d0*rdelta)
      enddo

      i=nlong-1
      zdr(i)=(ax(4,1)*zwr(i-3)
     &      +ax(4,2)*zwr(i-2)
     &      +ax(4,3)*zwr(i-1)
     &      +ax(4,4)*zwr(i)
     &      +ax(4,5)*zwr(i+1))
     &     /(24.d0*rdelta)
      
      i=nlong
      zdr(i)=(ax(5,1)*zwr(i-4)
     &      +ax(5,2)*zwr(i-3)
     &      +ax(5,3)*zwr(i-2)
     &      +ax(5,4)*zwr(i-1)
     &      +ax(5,5)*zwr(i))
     &     /(24.d0*rdelta)

      do i=1,nlong
         zwr(i)=zdr(i)
      enddo
c     -----------------------------------------------
      return

      end
c**************************************************************
      subroutine trd(trace,dmat,n)
c     ------------------------------------------------------
c     Trace of matrix 
c     ------------------------------------------------------
c
c     trace   ... Trace of matrix (Output)
c     dmat    ... Matrix (Input)
c     n       ... Size of Matrix (Input)
c
      implicit real*8(a-h,o-z)

      dimension dmat(n,n)

c-----------------------------------------------------------
      trace=0.d0
      do i=1,n
         trace=trace+dmat(i,i)
      enddo
c-----------------------------------------------------------
      return
      end
c-------------------------------------------------------------
c     Rotate molecular frame to dipole princepal axis
c-------------------------------------------------------------
      subroutine diprot(q,r,natom,nmax)

      implicit real*8(a-h,o-z)

      dimension q(nmax),r(3,nmax),s(3)
      dimension rotmat(3,3),list(3),w(3)

      data zero/1.d-6/

      rx=0.d0
      ry=0.d0
      rz=0.d0
      do i=1,natom
         rx=rx+q(i)*r(1,i)
         ry=ry+q(i)*r(2,i)
         rz=rz+q(i)*r(3,i)
      enddo
      rr=dsqrt(rx*rx+ry*ry+rz*rz)

      if (dabs(rr).le.zero) return

      zdz=dsqrt(rr**2-rz**2)

      if (zdz.le.zero) then

         rotmat(1,1)= 1.d0
         rotmat(1,2)= 0.d0
         rotmat(1,3)= 0.d0

         rotmat(2,1)= 0.d0
         rotmat(2,2)= 1.d0
         rotmat(2,3)= 0.d0
     
         rotmat(3,1)= 0.d0
         rotmat(3,2)= 0.d0
         rotmat(3,3)= rz/rr

      else

         rotmat(1,1)= rx*rz/(rr*zdz)
         rotmat(1,2)= ry*rz/(rr*zdz)
         rotmat(1,3)=-zdz/rr

         rotmat(2,1)=-ry/zdz
         rotmat(2,2)= rx/zdz
         rotmat(2,3)= 0.d0
     
         rotmat(3,1)= rx/rr
         rotmat(3,2)= ry/rr
         rotmat(3,3)= rz/rr

      endif

      do i=1,natom
         do j=1,3
            t=0.d0
            do k=1,3
               t=t+rotmat(j,k)*r(k,i)
            enddo
            s(j)=t
         enddo
         do j=1,3
            r(j,i)=s(j)
         enddo
      enddo
      
      return
      end
c---------------------------------------------------------
c     Gauss-Hermite inter-polation
c---------------------------------------------------------
      real*8 function hrho(ngrid,dk1,hwv,rk,rk3,yd,ill)
      implicit real*8(a-h,o-z)
      dimension rk(0:ngrid),hwv(0:ngrid)
      dimension yd(0:ngrid)
      n = int(rk3/dk1)+1
      if (n.ge.ngrid+1) n=ngrid
      if (rk3.ge.(dk1*dble(ngrid))) then 
         hrho=hwv(ngrid)
         return
      endif
      m = 0
c-----<numpac+
c     n     ... previous grid point number of inter-polation point
c     rk3   ... coordinate of inter-polation point
c     y     ... inter-polated value (output)
c     m     ... 0
c     ngrid ... number of input data
c     rk    ... coordinate set of descreat point
c     hwv   ... descreat function 
c     yd    ... out put of 1st derivative
c     ill   ... not 0
c
c     this subroutine is modified from herm31 of numpac
c     you can see detail of this subroutine in
c     http://netnumpac.fuis.fukui-u.ac.jp/cgi-bin/numpac/htoh?083.html
c
      call herd31(n,rk3,y,m,ngrid+1,rk(0),hwv(0),yd(0),ngrid+1,ill)
c-----+numpac>
      hrho=y
      return
      end
c---------------------------------------------------------
c     3D-FFT
c
c     This is service routine to connect ZFFT3D and 3D-RISM
c     PREFACTOR, -1, denotes shift of coordinate to a half of box.
c 
C---------------------------------------------------------
      subroutine ffte3d(a,n,d,inv)
      implicit real*8(a-h,o-z)
      complex*16 a
      logical init/.true./
      dimension a(n,n,n)
      save init
c
      if (init) then
         init = .false.
         call zfft3d(a, n, n, n, 0)
      end if
c
!$omp parallel do
      do iz=1,n
      do iy=1,n
      do ix=1,n
        if(mod(ix+iy+iz,2).eq.0) a(ix,iy,iz)=-a(ix,iy,iz)
      enddo
      enddo
      enddo
!$omp end parallel do
c-----<ffte+
      if (inv.eq.0) inve=-1
      if (inv.eq.1) inve= 1

      call zfft3d(a,n,n,n,inve)
c-----+ffte>

c  inv=0  k -> r
c      1  r -> k

      if(inv.eq.0) then
        fac=1.0d0/(d*n)**3
      else
        fac=(n*d)**3
      end if
!$omp parallel do
      do iz=1,n
      do iy=1,n
      do  ix=1,n
        if(mod(ix+iy+iz,2).eq.0) then
          a(ix,iy,iz)=-fac*a(ix,iy,iz)
        else
          a(ix,iy,iz)= fac*a(ix,iy,iz)
        end if
      enddo
      enddo
      enddo
!$omp end parallel do
      return
      end
c----------------------------------------------------------------
c     Fourier transform of 3D function with grid shift
c----------------------------------------------------------------
      subroutine ft3dfunc(func,ngrid3d,rdelta3d,inv)
c
      implicit real*8(a-h,o-z)
      complex*16 func
c
      dimension func(ngrid3d**3)
c
      NGRID3D2=NGRID3D**2
      PI=DACOS(-1.d0)
      DK3D=2.D0*PI/(RDELTA3D*DBLE(NGRID3D))
      NGSHIFT=NGRID3D/2+1
c----------------------------------------------------------------
c
c     Transform r -> k  (inv=1)
c      
      if (inv.eq.1) then
         
!$OMP PARALLEL DO PRIVATE(rz,ry,rx,k,dkr)
      DO KZ=1,NGRID3D
      RZ=RDELTA3D*DBLE(KZ-NGSHIFT)
      DO KY=1,NGRID3D
      RY=RDELTA3D*DBLE(KY-NGSHIFT)
      DO KX=1,NGRID3D
      RX=RDELTA3D*DBLE(KX-NGSHIFT)

         K=KX+(KY-1)*NGRID3D+(KZ-1)*NGRID3D2
         DKR=DK3D/2.D0*(RX+RY+RZ) 
         func(k)=func(k)*cdexp(dcmplx(0.d0,dkr))
      
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      call ffte3d(func,ngrid3d,rdelta3d,inv)
c
c     Transform k -> r  (inv=0)
c      
      elseif (inv.eq.0) then

      call ffte3d(func,ngrid3d,rdelta3d,inv)

!$OMP PARALLEL DO PRIVATE(rz,ry,rx,k,dkr)
      DO KZ=1,NGRID3D
      RZ=RDELTA3D*DBLE(KZ-NGSHIFT)
      DO KY=1,NGRID3D
      RY=RDELTA3D*DBLE(KY-NGSHIFT)
      DO KX=1,NGRID3D
      RX=RDELTA3D*DBLE(KX-NGSHIFT)

         K=KX+(KY-1)*NGRID3D+(KZ-1)*NGRID3D2
         DKR=DK3D/2.D0*(RX+RY+RZ) 
         func(k)=func(k)*cdexp(dcmplx(0.d0,-dkr))
      
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      else

         write(*,*) "Invalid inv option in ft3dfunc."
         ierr=578
         call abrt(ierr)
         
      endif
c----------------------------------------------------------------
      return
      end
c----------------------------------------------------------------
