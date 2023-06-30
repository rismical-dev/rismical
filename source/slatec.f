*DECK SINT
      SUBROUTINE SINT (N, X, WSAVE)
C***BEGIN PROLOGUE  SINT
C***PURPOSE  Compute the sine transform of a real, odd sequence.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A3
C***TYPE      SINGLE PRECISION (SINT-S)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine SINT computes the discrete Fourier sine transform
C  of an odd sequence X(I).  The transform is defined below at
C  output parameter X.
C
C  SINT is the unnormalized inverse of itself since a call of SINT
C  followed by another call of SINT will multiply the input sequence
C  X by 2*(N+1).
C
C  The array WSAVE which is used by subroutine SINT must be
C  initialized by calling subroutine SINTI(N,WSAVE).
C
C  Input Parameters
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N+1 is the product of small primes.
C
C  X       an array which contains the sequence to be transformed
C
C
C  WSAVE   a work array with dimension at least INT(3.5*N+16)
C          in the program that calls SINT.  The WSAVE array must be
C          initialized by calling subroutine SINTI(N,WSAVE), and a
C          different WSAVE array must be used for each different
C          value of N.  This initialization does not have to be
C          repeated so long as N remains unchanged.  Thus subsequent
C          transforms can be obtained faster than the first.
C
C  Output Parameters
C
C  X       For I=1,...,N
C
C               X(I)= the sum from K=1 to K=N
C
C                    2*X(K)*SIN(K*I*PI/(N+1))
C
C               A call of SINT followed by another call of
C               SINT will multiply the sequence X by 2*(N+1).
C               Hence SINT is the unnormalized inverse
C               of itself.
C
C  WSAVE   contains initialization calculations which must not be
C          destroyed between calls of SINT.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTF
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing definition of variable SQRT3 by using
C               FORTRAN intrinsic function SQRT instead of a DATA
C               statement.
C   881128  Modified by Dick Valent to meet prologue standards.
C   891009  Removed unreferenced statement label.  (WRB)
C   891009  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SINT
      DIMENSION X(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINT
      SQRT3 = SQRT(3.)
      IF (N-2) 101,102,103
  101 X(1) = X(1)+X(1)
      RETURN
  102 XH = SQRT3*(X(1)+X(2))
      X(2) = SQRT3*(X(1)-X(2))
      X(1) = XH
      RETURN
  103 NP1 = N+1
      NS2 = N/2
      WSAVE(1) = 0.
      KW = NP1
      DO 104 K=1,NS2
         KW = KW+1
         KC = NP1-K
         T1 = X(K)-X(KC)
         T2 = WSAVE(KW)*(X(K)+X(KC))
         WSAVE(K+1) = T1+T2
         WSAVE(KC+1) = T2-T1
  104 CONTINUE
      MODN = MOD(N,2)
      IF (MODN .NE. 0) WSAVE(NS2+2) = 4.*X(NS2+1)
      NF = NP1+NS2+1
      CALL RFFTF (NP1,WSAVE,WSAVE(NF))
      X(1) = .5*WSAVE(1)
      DO 105 I=3,N,2
         X(I-1) = -WSAVE(I)
         X(I) = X(I-2)+WSAVE(I-1)
  105 CONTINUE
      IF (MODN .NE. 0) RETURN
      X(N) = -WSAVE(N+1)
      RETURN
      END
*DECK RFFTF
      SUBROUTINE RFFTF (N, R, WSAVE)
C***BEGIN PROLOGUE  RFFTF
C***SUBSIDIARY
C***PURPOSE  Compute the forward transform of a real, periodic sequence.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (RFFTF-S, CFFTF-C)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C   ********************************************************************
C   *   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   *
C   ********************************************************************
C   *                                                                  *
C   *   This routine uses non-standard Fortran 77 constructs and will  *
C   *   be removed from the library at a future date.  You are         *
C   *   requested to use RFFTF1.                                       *
C   *                                                                  *
C   ********************************************************************
C
C   Subroutine RFFTF computes the Fourier coefficients of a real
C   periodic sequence (Fourier analysis).  The transform is defined
C   below at output parameter R.
C
C   Input Arguments
C
C   N       the length of the array R to be transformed.  The method
C           is most efficient when N is a product of small primes.
C           N may change so long as different work arrays are provided.
C
C   R       a real array of length N which contains the sequence
C           to be transformed.
C
C   WSAVE   a work array which must be dimensioned at least 2*N+15
C           in the program that calls RFFTF.  The WSAVE array must be
C           initialized by calling subroutine RFFTI, and a different
C           WSAVE array must be used for each different value of N.
C           This initialization does not have to be repeated so long as
C           remains unchanged.  Thus subsequent transforms can be
C           obtained faster than the first.  Moreover, the same WSAVE
C           array can be used by RFFTF and RFFTB as long as N remains
C           unchanged.
C
C   Output Argument
C
C   R       R(1) = the sum from I=1 to I=N of R(I)
C
C           If N is even set L = N/2; if N is odd set L = (N+1)/2
C
C             then for K = 2,...,L
C
C                R(2*K-2) = the sum from I = 1 to I = N of
C
C                     R(I)*COS((K-1)*(I-1)*2*PI/N)
C
C                R(2*K-1) = the sum from I = 1 to I = N of
C
C                    -R(I)*SIN((K-1)*(I-1)*2*PI/N)
C
C           If N is even
C
C                R(N) = the sum from I = 1 to I = N of
C
C                     (-1)**(I-1)*R(I)
C
C   Note:  This transform is unnormalized since a call of RFFTF
C          followed by a call of RFFTB will multiply the input
C          sequence by N.
C
C   WSAVE  contains results which must not be destroyed between
C          calls of RFFTF or RFFTB.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTF1
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*).
C   861211  REVISION DATE from Version 3.2
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900131  Routine changed from user-callable to subsidiary
C           because of non-standard Fortran 77 arguments in the
C           call to CFFTB1.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RFFTF
      DIMENSION R(*), WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  RFFTF
      IF (N .EQ. 1) RETURN
      CALL RFFTF1 (N,R,WSAVE,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
*DECK RFFTF1
      SUBROUTINE RFFTF1 (N, C, CH, WA, IFAC)
C***BEGIN PROLOGUE  RFFTF1
C***PURPOSE  Compute the forward transform of a real, periodic sequence.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (RFFTF1-S, CFFTF1-C)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C   Subroutine RFFTF1 computes the Fourier coefficients of a real
C   periodic sequence (Fourier analysis).  The transform is defined
C   below at output parameter C.
C
C   The arrays WA and IFAC which are used by subroutine RFFTB1 must be
C   initialized by calling subroutine RFFTI1.
C
C   Input Arguments
C
C   N       the length of the array R to be transformed.  The method
C           is most efficient when N is a product of small primes.
C           N may change so long as different work arrays are provided.
C
C   C       a real array of length N which contains the sequence
C           to be transformed.
C
C   CH      a real work array of length at least N.
C
C   WA      a real work array which must be dimensioned at least N.
C
C   IFAC    an integer work array which must be dimensioned at least 15.
C
C           The WA and IFAC arrays must be initialized by calling
C           subroutine RFFTI1, and different WA and IFAC arrays must be
C           used for each different value of N.  This initialization
C           does not have to be repeated so long as N remains unchanged.
C           Thus subsequent transforms can be obtained faster than the
C           first.  The same WA and IFAC arrays can be used by RFFTF1
C           and RFFTB1.
C
C   Output Argument
C
C   C       C(1) = the sum from I=1 to I=N of R(I)
C
C           If N is even set L = N/2; if N is odd set L = (N+1)/2
C
C             then for K = 2,...,L
C
C                C(2*K-2) = the sum from I = 1 to I = N of
C
C                     C(I)*COS((K-1)*(I-1)*2*PI/N)
C
C                C(2*K-1) = the sum from I = 1 to I = N of
C
C                    -C(I)*SIN((K-1)*(I-1)*2*PI/N)
C
C           If N is even
C
C                C(N) = the sum from I = 1 to I = N of
C
C                     (-1)**(I-1)*C(I)
C
C   Notes:  This transform is unnormalized since a call of RFFTF1
C           followed by a call of RFFTB1 will multiply the input
C           sequence by N.
C
C           WA and IFAC contain initialization calculations which must
C           not be destroyed between calls of subroutine RFFTF1 or
C           RFFTB1.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RADF2, RADF3, RADF4, RADF5, RADFG
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*).
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900131  Routine changed from subsidiary to user-callable.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RFFTF1
      DIMENSION CH(*), C(*), WA(*), IFAC(*)
C***FIRST EXECUTABLE STATEMENT  RFFTF1
      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = IFAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL RADF2 (IDO,L1,C,CH,WA(IW))
         GO TO 110
  103    CALL RADF2 (IDO,L1,CH,C,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
         GO TO 110
  105    CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
         CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107    CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         NA = 1
         GO TO 110
  109    CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      IF (NA .EQ. 1) RETURN
      DO 112 I=1,N
         C(I) = CH(I)
  112 CONTINUE
      RETURN
      END
*DECK RADF2
      SUBROUTINE RADF2 (IDO, L1, CC, CH, WA1)
C***BEGIN PROLOGUE  RADF2
C***SUBSIDIARY
C***PURPOSE  Calculate the fast Fourier transform of subvectors of
C            length two.
C***LIBRARY   SLATEC (FFTPACK)
C***TYPE      SINGLE PRECISION (RADF2-S)
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*).
C   881128  Modified by Dick Valent to meet prologue standards.
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  RADF2
      DIMENSION CH(IDO,2,*), CC(IDO,L1,2), WA1(*)
C***FIRST EXECUTABLE STATEMENT  RADF2
      DO 101 K=1,L1
         CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 108
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=3,IDO,2
            IC = IDP2-I
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CH(I,1,K) = CC(I,K,1)+TI2
            CH(IC,2,K) = TI2-CC(I,K,1)
            CH(I-1,1,K) = CC(I-1,K,1)+TR2
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  103    CONTINUE
  104 CONTINUE
      GO TO 111
  108 DO 110 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 109 K=1,L1
            TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CH(I,1,K) = CC(I,K,1)+TI2
            CH(IC,2,K) = TI2-CC(I,K,1)
            CH(I-1,1,K) = CC(I-1,K,1)+TR2
            CH(IC-1,2,K) = CC(I-1,K,1)-TR2
  109    CONTINUE
  110 CONTINUE
  111 IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         CH(1,2,K) = -CC(IDO,K,2)
         CH(IDO,1,K) = CC(IDO,K,1)
  106 CONTINUE
  107 RETURN
      END
*DECK RADF3
      SUBROUTINE RADF3 (IDO, L1, CC, CH, WA1, WA2)
C***BEGIN PROLOGUE  RADF3
C***SUBSIDIARY
C***PURPOSE  Calculate the fast Fourier transform of subvectors of
C            length three.
C***LIBRARY   SLATEC (FFTPACK)
C***TYPE      SINGLE PRECISION (RADF3-S)
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing definition of variable TAUI by using
C               FORTRAN intrinsic function SQRT instead of a DATA
C               statement.
C   881128  Modified by Dick Valent to meet prologue standards.
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  RADF3
      DIMENSION CH(IDO,3,*), CC(IDO,L1,3), WA1(*), WA2(*)
C***FIRST EXECUTABLE STATEMENT  RADF3
      TAUR = -.5
      TAUI = .5*SQRT(3.)
      DO 101 K=1,L1
         CR2 = CC(1,K,2)+CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 104
      DO 103 K=1,L1
CDIR$ IVDEP
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2+DR3
            CI2 = DI2+DI3
            CH(I-1,1,K) = CC(I-1,K,1)+CR2
            CH(I,1,K) = CC(I,K,1)+CI2
            TR2 = CC(I-1,K,1)+TAUR*CR2
            TI2 = CC(I,K,1)+TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2+TR3
            CH(IC-1,2,K) = TR2-TR3
            CH(I,3,K) = TI2+TI3
            CH(IC,2,K) = TI3-TI2
  102    CONTINUE
  103 CONTINUE
      RETURN
  104 DO 106 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 105 K=1,L1
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2+DR3
            CI2 = DI2+DI3
            CH(I-1,1,K) = CC(I-1,K,1)+CR2
            CH(I,1,K) = CC(I,K,1)+CI2
            TR2 = CC(I-1,K,1)+TAUR*CR2
            TI2 = CC(I,K,1)+TAUR*CI2
            TR3 = TAUI*(DI2-DI3)
            TI3 = TAUI*(DR3-DR2)
            CH(I-1,3,K) = TR2+TR3
            CH(IC-1,2,K) = TR2-TR3
            CH(I,3,K) = TI2+TI3
            CH(IC,2,K) = TI3-TI2
  105    CONTINUE
  106 CONTINUE
      RETURN
      END
*DECK RADF4
      SUBROUTINE RADF4 (IDO, L1, CC, CH, WA1, WA2, WA3)
C***BEGIN PROLOGUE  RADF4
C***SUBSIDIARY
C***PURPOSE  Calculate the fast Fourier transform of subvectors of
C            length four.
C***LIBRARY   SLATEC (FFTPACK)
C***TYPE      SINGLE PRECISION (RADF4-S)
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*).
C           (b) changing definition of variable HSQT2 by using
C               FORTRAN intrinsic function SQRT instead of a DATA
C               statement.
C   881128  Modified by Dick Valent to meet prologue standards.
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  RADF4
      DIMENSION CC(IDO,L1,4), CH(IDO,4,*), WA1(*), WA2(*), WA3(*)
C***FIRST EXECUTABLE STATEMENT  RADF4
      HSQT2 = .5*SQRT(2.)
      DO 101 K=1,L1
         TR1 = CC(1,K,2)+CC(1,K,4)
         TR2 = CC(1,K,1)+CC(1,K,3)
         CH(1,1,K) = TR1+TR2
         CH(IDO,4,K) = TR2-TR1
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)
         CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 111
      DO 104 K=1,L1
CDIR$ IVDEP
         DO 103 I=3,IDO,2
            IC = IDP2-I
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            TR1 = CR2+CR4
            TR4 = CR4-CR2
            TI1 = CI2+CI4
            TI4 = CI2-CI4
            TI2 = CC(I,K,1)+CI3
            TI3 = CC(I,K,1)-CI3
            TR2 = CC(I-1,K,1)+CR3
            TR3 = CC(I-1,K,1)-CR3
            CH(I-1,1,K) = TR1+TR2
            CH(IC-1,4,K) = TR2-TR1
            CH(I,1,K) = TI1+TI2
            CH(IC,4,K) = TI1-TI2
            CH(I-1,3,K) = TI4+TR3
            CH(IC-1,2,K) = TR3-TI4
            CH(I,3,K) = TR4+TI3
            CH(IC,2,K) = TR4-TI3
  103    CONTINUE
  104 CONTINUE
      GO TO 110
  111 DO 109 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 108 K=1,L1
            CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            TR1 = CR2+CR4
            TR4 = CR4-CR2
            TI1 = CI2+CI4
            TI4 = CI2-CI4
            TI2 = CC(I,K,1)+CI3
            TI3 = CC(I,K,1)-CI3
            TR2 = CC(I-1,K,1)+CR3
            TR3 = CC(I-1,K,1)-CR3
            CH(I-1,1,K) = TR1+TR2
            CH(IC-1,4,K) = TR2-TR1
            CH(I,1,K) = TI1+TI2
            CH(IC,4,K) = TI1-TI2
            CH(I-1,3,K) = TI4+TR3
            CH(IC-1,2,K) = TR3-TI4
            CH(I,3,K) = TR4+TI3
            CH(IC,2,K) = TR4-TI3
  108    CONTINUE
  109 CONTINUE
  110 IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
         TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
         CH(IDO,1,K) = TR1+CC(IDO,K,1)
         CH(IDO,3,K) = CC(IDO,K,1)-TR1
         CH(1,2,K) = TI1-CC(IDO,K,3)
         CH(1,4,K) = TI1+CC(IDO,K,3)
  106 CONTINUE
  107 RETURN
      END
*DECK RADF5
      SUBROUTINE RADF5 (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
C***BEGIN PROLOGUE  RADF5
C***SUBSIDIARY
C***PURPOSE  Calculate the fast Fourier transform of subvectors of
C            length five.
C***LIBRARY   SLATEC (FFTPACK)
C***TYPE      SINGLE PRECISION (RADF5-S)
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing definition of variables PI, TI11, TI12,
C               TR11, TR12 by using FORTRAN intrinsic functions ATAN
C               and SIN instead of DATA statements.
C   881128  Modified by Dick Valent to meet prologue standards.
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  RADF5
      DIMENSION CC(IDO,L1,5), CH(IDO,5,*), WA1(*), WA2(*), WA3(*),
     +          WA4(*)
C***FIRST EXECUTABLE STATEMENT  RADF5
      PI = 4.*ATAN(1.)
      TR11 = SIN(.1*PI)
      TI11 = SIN(.4*PI)
      TR12 = -SIN(.3*PI)
      TI12 = SIN(.2*PI)
      DO 101 K=1,L1
         CR2 = CC(1,K,5)+CC(1,K,2)
         CI5 = CC(1,K,5)-CC(1,K,2)
         CR3 = CC(1,K,4)+CC(1,K,3)
         CI4 = CC(1,K,4)-CC(1,K,3)
         CH(1,1,K) = CC(1,K,1)+CR2+CR3
         CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
         CH(1,3,K) = TI11*CI5+TI12*CI4
         CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
         CH(1,5,K) = TI12*CI5-TI11*CI4
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      IF((IDO-1)/2.LT.L1) GO TO 104
      DO 103 K=1,L1
CDIR$ IVDEP
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2+DR5
            CI5 = DR5-DR2
            CR5 = DI2-DI5
            CI2 = DI2+DI5
            CR3 = DR3+DR4
            CI4 = DR4-DR3
            CR4 = DI3-DI4
            CI3 = DI3+DI4
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
            CH(I,1,K) = CC(I,K,1)+CI2+CI3
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
            TR5 = TI11*CR5+TI12*CR4
            TI5 = TI11*CI5+TI12*CI4
            TR4 = TI12*CR5-TI11*CR4
            TI4 = TI12*CI5-TI11*CI4
            CH(I-1,3,K) = TR2+TR5
            CH(IC-1,2,K) = TR2-TR5
            CH(I,3,K) = TI2+TI5
            CH(IC,2,K) = TI5-TI2
            CH(I-1,5,K) = TR3+TR4
            CH(IC-1,4,K) = TR3-TR4
            CH(I,5,K) = TI3+TI4
            CH(IC,4,K) = TI4-TI3
  102    CONTINUE
  103 CONTINUE
      RETURN
  104 DO 106 I=3,IDO,2
         IC = IDP2-I
CDIR$ IVDEP
         DO 105 K=1,L1
            DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2+DR5
            CI5 = DR5-DR2
            CR5 = DI2-DI5
            CI2 = DI2+DI5
            CR3 = DR3+DR4
            CI4 = DR4-DR3
            CR4 = DI3-DI4
            CI3 = DI3+DI4
            CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
            CH(I,1,K) = CC(I,K,1)+CI2+CI3
            TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
            TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
            TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
            TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
            TR5 = TI11*CR5+TI12*CR4
            TI5 = TI11*CI5+TI12*CI4
            TR4 = TI12*CR5-TI11*CR4
            TI4 = TI12*CI5-TI11*CI4
            CH(I-1,3,K) = TR2+TR5
            CH(IC-1,2,K) = TR2-TR5
            CH(I,3,K) = TI2+TI5
            CH(IC,2,K) = TI5-TI2
            CH(I-1,5,K) = TR3+TR4
            CH(IC-1,4,K) = TR3-TR4
            CH(I,5,K) = TI3+TI4
            CH(IC,4,K) = TI4-TI3
  105    CONTINUE
  106 CONTINUE
      RETURN
      END
*DECK RADFG
      SUBROUTINE RADFG (IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
C***BEGIN PROLOGUE  RADFG
C***SUBSIDIARY
C***PURPOSE  Calculate the fast Fourier transform of subvectors of
C            arbitrary length.
C***LIBRARY   SLATEC (FFTPACK)
C***TYPE      SINGLE PRECISION (RADFG-S)
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing references to intrinsic function FLOAT
C               to REAL, and
C           (c) changing definition of variable TPI by using
C               FORTRAN intrinsic function ATAN instead of a DATA
C               statement.
C   881128  Modified by Dick Valent to meet prologue standards.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900402  Added TYPE section.  (WRB)
C***END PROLOGUE  RADFG
      DIMENSION CH(IDO,L1,*), CC(IDO,IP,*), C1(IDO,L1,*),
     +          C2(IDL1,*), CH2(IDL1,*), WA(*)
C***FIRST EXECUTABLE STATEMENT  RADFG
      TPI = 8.*ATAN(1.)
      ARG = TPI/IP
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         CH2(IK,1) = C2(IK,1)
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            CH(1,K,J) = C1(1,K,J)
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
CDIR$ IVDEP
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
               CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
CDIR$ IVDEP
            DO 112 I=3,IDO,2
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
               C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
               C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
               C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)
            C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.
      AI1 = 0.
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)
            CH2(IK,LC) = AI1*C2(IK,IP)
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)
               CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+C2(IK,J)
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            CC(I,1,K) = CH(I,K,1)
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            CC(I,1,K) = CH(I,K,1)
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            CC(IDO,J2-2,K) = CH(1,K,J)
            CC(1,J2-1,K) = CH(1,K,JC)
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
CDIR$ IVDEP
            DO 138 I=3,IDO,2
               IC = IDP2-I
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
               CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
               CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
               CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END
*DECK SINTI
      SUBROUTINE SINTI (N, WSAVE)
C***BEGIN PROLOGUE  SINTI
C***PURPOSE  Initialize a work array for SINT.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A3
C***TYPE      SINGLE PRECISION (SINTI-S)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C  Subroutine SINTI initializes the array WSAVE which is used in
C  subroutine SINT.  The prime factorization of N together with
C  a tabulation of the trigonometric functions are computed and
C  stored in WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  The method
C          is most efficient when N+1 is a product of small primes.
C
C  Output Parameter
C
C  WSAVE   a work array with at least INT(3.5*N+16) locations.
C          Different WSAVE arrays are required for different values
C          of N.  The contents of WSAVE must not be changed between
C          calls of SINT.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTI
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing references to intrinsic function FLOAT
C               to REAL, and
C           (c) changing definition of variable PI by using
C               FORTRAN intrinsic function ATAN instead of a DATA
C               statement.
C   881128  Modified by Dick Valent to meet prologue standards.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  SINTI
      DIMENSION WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  SINTI
      IF (N .LE. 1) RETURN
      PI = 4.*ATAN(1.)
      NP1 = N+1
      NS2 = N/2
      DT = PI/NP1
      KS = N+2
      KF = KS+NS2-1
      FK = 0.
      DO 101 K=KS,KF
         FK = FK+1.
         WSAVE(K) = 2.*SIN(FK*DT)
  101 CONTINUE
      CALL RFFTI (NP1,WSAVE(KF+1))
      RETURN
      END
*DECK RFFTI
      SUBROUTINE RFFTI (N, WSAVE)
C***BEGIN PROLOGUE  RFFTI
C***SUBSIDIARY
C***PURPOSE  Initialize a work array for RFFTF and RFFTB.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (RFFTI-S, CFFTI-C)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C   ********************************************************************
C   *   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   NOTICE   *
C   ********************************************************************
C   *                                                                  *
C   *   This routine uses non-standard Fortran 77 constructs and will  *
C   *   be removed from the library at a future date.  You are         *
C   *   requested to use RFFTI1.                                       *
C   *                                                                  *
C   ********************************************************************
C
C   Subroutine RFFTI initializes the array WSAVE which is used in
C   both RFFTF and RFFTB.  The prime factorization of N together with
C   a tabulation of the trigonometric functions are computed and
C   stored in WSAVE.
C
C   Input Argument
C
C   N       the length of the sequence to be transformed.
C
C   Output Argument
C
C   WSAVE   a work array which must be dimensioned at least 2*N+15.
C           The same work array can be used for both RFFTF and RFFTB
C           as long as N remains unchanged.  Different WSAVE arrays
C           are required for different values of N.  The contents of
C           WSAVE must not be changed between calls of RFFTF or RFFTB.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  RFFTI1
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           changing dummy array size declarations (1) to (*).
C   861211  REVISION DATE from Version 3.2
C   881128  Modified by Dick Valent to meet prologue standards.
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900131  Routine changed from user-callable to subsidiary
C           because of non-standard Fortran 77 arguments in the
C           call to CFFTB1.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RFFTI
      DIMENSION WSAVE(*)
C***FIRST EXECUTABLE STATEMENT  RFFTI
      IF (N .EQ. 1) RETURN
      CALL RFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
      RETURN
      END
*DECK RFFTI1
      SUBROUTINE RFFTI1 (N, WA, IFAC)
C***BEGIN PROLOGUE  RFFTI1
C***PURPOSE  Initialize a real and an integer work array for RFFTF1 and
C            RFFTB1.
C***LIBRARY   SLATEC (FFTPACK)
C***CATEGORY  J1A1
C***TYPE      SINGLE PRECISION (RFFTI1-S, CFFTI1-C)
C***KEYWORDS  FFTPACK, FOURIER TRANSFORM
C***AUTHOR  Swarztrauber, P. N., (NCAR)
C***DESCRIPTION
C
C   Subroutine RFFTI1 initializes the work arrays WA and IFAC which are
C   used in both RFFTF1 and RFFTB1.  The prime factorization of N and a
C   tabulation of the trigonometric functions are computed and stored in
C   IFAC and WA, respectively.
C
C   Input Argument
C
C   N       the length of the sequence to be transformed.
C
C   Output Arguments
C
C   WA      a real work array which must be dimensioned at least N.
C
C   IFAC    an integer work array which must be dimensioned at least 15.
C
C   The same work arrays can be used for both RFFTF1 and RFFTB1 as long
C   as N remains unchanged.  Different WA and IFAC arrays are required
C   for different values of N.  The contents of WA and IFAC must not be
C   changed between calls of RFFTF1 or RFFTB1.
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C                 Computations (G. Rodrigue, ed.), Academic Press,
C                 1982, pp. 51-83.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790601  DATE WRITTEN
C   830401  Modified to use SLATEC library source file format.
C   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
C           (a) changing dummy array size declarations (1) to (*),
C           (b) changing references to intrinsic function FLOAT
C               to REAL, and
C           (c) changing definition of variable TPI by using
C               FORTRAN intrinsic functions instead of DATA
C               statements.
C   881128  Modified by Dick Valent to meet prologue standards.
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900131  Routine changed from subsidiary to user-callable.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RFFTI1
      DIMENSION WA(*), IFAC(*), NTRYH(4)
      SAVE NTRYH
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
C***FIRST EXECUTABLE STATEMENT  RFFTI1
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         IFAC(IB+2) = IFAC(IB+1)
  106 CONTINUE
      IFAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.*ATAN(1.)
      ARGH = TPI/N
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = LD*ARGH
            FI = 0.
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END
