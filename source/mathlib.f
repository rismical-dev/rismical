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
      
c-------< Use NUMPAC>
c      parameter (NFFTX=14)
c      parameter (NFFT=2**(NFFTX-1)-1)
c      common /FFTDAT/initw,nfftxx,trimat(nfft)
      
      dimension func(ngrid),dum(ngrid)
      
      pi=dacos(-1.d0)
      deltak=pi/(dble(ngrid)*rdelta)
c---------------------------------------------------------
c-------< Use NUMPAC>
c      if (initw.eq.0) then
c         icon=0
c         call trigqd(trimat,NFFTX-1,icon)
c         initw=1
c         ngtemp=ngrid
c         nfftxx=0
c         do while (ngtemp.ne.1)
c            ngtemp=ngtemp/2
c            nfftxx=nfftxx+1
c         enddo
c         if (2**nfftx.lt.ngrid) then
c            write(*,*) "Error, Number of FFT grid is not enough."
c            write(*,*) "You must recompile with larger parameter",
c     &                 " in fft1d."
c            stop
c         endif
c      endif

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
C--------------------------------------------------------------------
      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      REAL*8 data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL*8 tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=dble(wr)*data(j)-dble(wi)*data(j+1)
            tempi=dble(wr)*data(j+1)+dble(wi)*data(j)
c            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
c            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software.
C--------------------------------------------------------------------
      SUBROUTINE lubksbx(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software.
C--------------------------------------------------------------------
      SUBROUTINE ludcmpx(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) then
           write(*,*) 'singular matrix in ludcmp'
           call abrt
        endif
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software.
C--------------------------------------------------------------------
      SUBROUTINE realft(data,n,isign)
      INTEGER isign,n
      REAL*8 data(n)
CU    USES four1
      INTEGER i,i1,i2,i3,i4,n2p3
      REAL*8 c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5d0
      if (isign.eq.1) then
        c2=-0.5d0
        call four1(data,n/2,+1)
      else
        c2=0.5d0
        theta=-theta
      endif
      wpr=-2.0d0*dsin(0.5d0*theta)**2
      wpi=dsin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=dble(wr)
        wis=dble(wi)
c        wrs=sngl(wr)
c        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n/2,-1)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software.
C--------------------------------------------------------------------
      SUBROUTINE sinft(y,n)
      INTEGER n
      REAL*8 y(n)
CU    USES realft
      INTEGER j
      REAL*8 sum,y1,y2
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n)
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*dsin(0.5d0*theta)**2
      wpi=dsin(theta)
      y(1)=0.0d0
      do 11 j=1,n/2
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        y1=wi*(y(j+1)+y(n-j+1))
        y2=0.5d0*(y(j+1)-y(n-j+1))
        y(j+1)=y1+y2
        y(n-j+1)=y1-y2
11    continue
      call realft(y,n,+1)
      sum=0.0d0
      y(1)=0.5d0*y(1)
      y(2)=0.0d0
      do 12 j=1,n-1,2
        sum=sum+y(j)
        y(j)=y(j+1)
        y(j+1)=sum
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software.
C--------------------------------------------------------------------
      SUBROUTINE polin2d(x1a,x2a,ya,m,n,x1,x2,y,dy)
      INTEGER m,n,NMAX,MMAX
      REAL*8 dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      PARAMETER (NMAX=20,MMAX=20)
CU    USES polint
      INTEGER j,k
      REAL*8 ymtmp(MMAX),yntmp(NMAX)
      do 12 j=1,m
        do 11 k=1,n
          yntmp(k)=ya(j,k)
11      continue
        call polintd(x2a,yntmp,n,x2,ymtmp(j),dy)
12    continue
      call polintd(x1a,ymtmp,m,x1,y,dy)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
C--------------------------------------------------------------------
      SUBROUTINE polintd(xa,ya,n,x,y,dy)
      INTEGER n
c      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
c      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(N),d(N)
c      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0)then
             write(*,*) 'failure in polint'
             call abrt
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
C--------------------------------------------------------------------
      SUBROUTINE TRAPZD(FUNC,N,DEL,RESULT)
      INTEGER N
      REAL*8 DEL,RESULT,FUNC(N),SUM
      INTEGER J

      SUM=0.D0
      DO 11 J=2,N-1
         SUM=SUM+FUNC(J)
 11   CONTINUE
      RESULT=DEL*(SUM+(FUNC(1)+FUNC(N))*0.5D0)
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
C     MODIFIED By NY
C--------------------------------------------------------------------
      REAL*8 FUNCTION plgndr(l,m,x)
      INTEGER l,m
      REAL*8 x
      INTEGER i,ll
      REAL*8 fact,pll,pmm,pmmp1,somx2
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0) then
         write (*,*) 'bad arguments in plgndr'
         call abrt
      endif
      pmm=1.d0
      if(m.gt.0) then
        somx2=sqrt((1.d0-x)*(1.d0+x))
        fact=1.d0
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*dble(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm)/dble(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
C--------------------------------------------------------------------
      FUNCTION bessj(n,x)
      INTEGER n,IACC
      REAL*8 bessj,x,BIGNO,BIGNI
      PARAMETER (IACC=40,BIGNO=1.e10,BIGNI=1.e-10)
CU    USES bessj0,bessj1
      INTEGER j,jsum,m
      REAL*8 ax,bj,bjm,bjp,sum,tox,bessj0,bessj1
      if(n.lt.2) then
c     Mod by Norio
         if (n.eq.0) then
            bessj=bessj0(x)
         else
            bessj=bessj1(x)
         endif
         return
c
      endif
      ax=abs(x)
      if(ax.eq.0.d0)then
        bessj=0.d0
      else if(ax.gt.float(n))then
        tox=2.d0/ax
        bjm=bessj0(ax)
        bj=bessj1(ax)
        do 11 j=1,n-1
          bjp=j*tox*bj-bjm
          bjm=bj
          bj=bjp
11      continue
        bessj=bj
      else
        tox=2.d0/ax
        m=2*((n+int(sqrt(float(IACC*n))))/2)
        bessj=0.d0
        jsum=0
        sum=0.d0
        bjp=0.d0
        bj=1.d0
        do 12 j=m,1,-1
          bjm=j*tox*bj-bjp
          bjp=bj
          bj=bjm
          if(abs(bj).gt.BIGNO)then
            bj=bj*BIGNI
            bjp=bjp*BIGNI
            bessj=bessj*BIGNI
            sum=sum*BIGNI
          endif
          if(jsum.ne.0)sum=sum+bj
          jsum=1-jsum
          if(j.eq.n)bessj=bjp
12      continue
        sum=2.*sum-bj
        bessj=bessj/sum
      endif
      if(x.lt.0.d0.and.mod(n,2).eq.1)bessj=-bessj
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ?421.1.9.
C--------------------------------------------------------------------
      FUNCTION bessj0(x)
      REAL*8 bessj0,x
      REAL*8 ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     *-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     *.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     *651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,
     *s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,
     *59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.d0)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     *(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-.785398164d0
        bessj0=dsqrt(.636619772/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ?421.1.9.
C--------------------------------------------------------------------
      FUNCTION bessj1(x)
      REAL*8 bessj1,x
      REAL*8 ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     *s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     *99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     *-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.d0)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-2.356194491d0
        bessj1=dsqrt(.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.d0,x)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ?421.1.9.
C--------------------------------------------------------------------
C#NUMPAC#FFT3DB PARALLEL VERSION  REVISED ON 1987-08-07
      SUBROUTINE FFT3DB(A,KA,LA,M,INV,B,ILL)                            
      COMPLEX*16 A(KA,LA,1),B(1)                                        
      INTEGER M(3)                                                      
      COMPLEX*16, ALLOCATABLE::C(:,:,:)
C     
      IF(M(1).LE.1.OR.M(2).LE.1.OR.M(3).LE.1)GO TO 60                   
      NROW=2**M(1)                                                      
      NCOL=2**M(2)                                                      
      NLAY=2**M(3)                                                      
      IF(NROW.GT.KA.OR.NCOL.GT.LA) GO TO 60                             
      ALLOCATE (C(NLAY,KA,LA))
C
C     ---- FFT FOR X AXIS ----
C
      DO 15 K=1,NLAY
C
C
      DO 10 J=1,NCOL                                                    
 10   CALL FFTB(A(1,J,K),M(1),INV,ILL)                       
 15   CONTINUE
C
C     ---- FFT FOR Y AXIS ----
C
      DO 35 K=1,NLAY
C

      DO 30 I=1,NROW                             
      DO 20 J=1,NCOL                             
 20   B(J)=A(I,J,K)                              
      CALL FFTB(B,M(2),INV,ILL)                  
      DO 30 J=1,NCOL                             
 30      A(I,J,K)=B(J)                              
 35   CONTINUE                                  
C
      DO J=1,NCOL
         DO I=1,NROW
            DO K=1,NLAY
               C(K,I,J)=A(I,J,K)
            ENDDO
         ENDDO
      ENDDO
C
C     ---- FFT FOR Z AXIS ----
C
      DO 55 J=1,NCOL                             
      DO 50 I=1,NROW                             
      DO 40 K=1,NLAY                             
 40   B(K)=C(K,I,J)
      CALL FFTB(B,M(3),INV,ILL)
      DO 50 K=1,NLAY                             
 50   C(K,I,J)=B(K)
 55   CONTINUE
C
      DO J=1,NCOL
         DO I=1,NROW
            DO K=1,NLAY
               A(I,J,K)=C(K,I,J)
            ENDDO
         ENDDO
      ENDDO
C
      ILL=0                                      
C
      DEALLOCATE (C)
      RETURN                                     
C     
 60   ILL=30000                                  
      RETURN                                     
      END                                        
C#NUMPAC#FFT3DB ORIGINAL VERSION REVISED ON 1987-08-07
      SUBROUTINE FFT3DB_ORG(A,KA,LA,M,INV,B,ILL)                        
      COMPLEX*16 A(KA,LA,1),B(1)                                        
      INTEGER M(3)                                                      
      IF(M(1).LE.1.OR.M(2).LE.1.OR.M(3).LE.1)GO TO 60                   
      NROW=2**M(1)                                                      
      NCOL=2**M(2)                                                      
      NLAY=2**M(3)                                                      
      IF(NROW.GT.KA.OR.NCOL.GT.LA) GO TO 60                             
      DO 10 K=1,NLAY                                                    
      DO 10 J=1,NCOL                                                    
   10 CALL FFTB(A(1,J,K),M(1),INV,ILL)                                  
      DO 30 K=1,NLAY                                                    
      DO 30 I=1,NROW                                                    
      DO 20 J=1,NCOL                                                    
   20 B(J)=A(I,J,K)                                                     
      CALL FFTB(B,M(2),INV,ILL)                                         
      DO 30 J=1,NCOL                                                    
   30 A(I,J,K)=B(J)                                                     
      DO 50 J=1,NCOL                                                    
      DO 50 I=1,NROW                                                    
      DO 40 K=1,NLAY                                                    
   40 B(K)=A(I,J,K)                                                     
      CALL FFTB(B,M(3),INV,ILL)                                         
      DO 50 K=1,NLAY                                                    
   50 A(I,J,K)=B(K)                                                     
      ILL=0                                                             
      RETURN                                                            
   60 ILL=30000                                                         
      RETURN                                                            
      END                                                               
C#NUMPAC#HERD31              REVISED ON 1984-11-30                      
C      SUBROUTINE HERM31(I,X,Y,M,N,XI,YI,YD,ND,ILL)                     
      SUBROUTINE HERD31(I,X,Y,M,N,XI,YI,YD,ND,ILL)                      
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XI(N),YI(N),YD(ND),Z(2,2)                               
      IF(N.LE.1.OR.M.LT.0.OR.M.GE.3) GO TO 1010                         
      IF(ILL.EQ.0) GO TO 1000                                           
      MD=ILL                                                            
      IF(N.LE.MD) GO TO 1010                                            
C      CALL DERIV1(XI,YI,YD,N,1,1,MD)                             
       CALL DERID1(XI,YI,YD,N,1,1,MD)                                
 1000 CONTINUE                                                          
      Z(1,1)=YI(I)                                                      
      Z(1,2)=YD(I)                                                      
      Z(2,1)=YI(I+1)                                                    
      Z(2,2)=YD(I+1)                                                    
C      CALL PHER31(X,Y,M,XI(I),Z,ILL)                        
      CALL PHED31(X,Y,M,XI(I),Z,ILL)                                    
      RETURN                                                            
 1010 ILL=30000                                                         
      RETURN                                                            
      END                                                               
C#NUMPAC#DERIV1              REVISED ON 1984-11-30                      
C      SUBROUTINE DERIV1(X,Y,F,  N,M,L,MD)                      
      SUBROUTINE DERID1(X,Y,F,  N,M,L,MD)                               
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),Y(M),F(L)                                          
      N1=N-1                                                            
      IF(MD.GE.2) GO TO 1010                                            
      IE=N1*L+1                                                         
      I=N1*M+1                                                          
      IM=I-M                                                            
      F(IE)=(Y(I)-Y(IM))/(X(N)-X(N1))                                   
      DO 1000 I=1,N1                                                    
      IM=I-1                                                            
      II=IM*M+1                                                         
      IE=IM*L+1                                                         
      IP=II+M                                                           
 1000 F(IE)=(Y(IP)-Y(II))/(X(I+1)-X(I))                                 
      RETURN                                                            
 1010 CONTINUE                                                          
      H1=X(2)-X(1)                                                      
      H2=X(3)-X(2)                                                      
      H1P2=H1+H2                                                        
      T1=(Y(M+1)-Y(1))/H1                                               
      T2=(Y(2*M+1)-Y(M+1))/H2                                           
      F(1)=T1+H1/H1P2*(T1-T2)                                           
      IE=N1*L+1                                                         
      I=N1*M+1                                                          
      IM=I-M                                                            
      IMM=IM-M                                                          
      H1=X(N1)-X(N1-1)                                                  
      H2=X(N)-X(N1)                                                     
      H1P2=H1+H2                                                        
      T1=(Y(IM)-Y(IMM))/H1                                              
      T2=(Y(I)-Y(IM))/H2                                                
      F(IE)=(H1P2+H2)/H1P2*(T2-T1)+T1                                   
      M1=1                                                              
      IF(MD.GE.3) M1=N-3                                                
      DO 1020 I=2,N1,M1                                                 
      IM=I-1                                                            
      II=IM*M+1                                                         
      IE=IM*L+1                                                         
      H1=X(I)-X(IM)                                                     
      H2=X(I+1)-X(I)                                                    
      H1P2=H1+H2                                                        
      IPM=II+M                                                          
      IMM=II-M                                                          
      T1=(Y(II)-Y(IMM))/H1                                              
      T2=(Y(IPM)-Y(II))/H2                                              
      F(IE)=(T2*H1+T1*H2)/H1P2                                          
 1020 CONTINUE                                                          
      IF(MD.LE.2) RETURN                                                
 1030 CONTINUE                                                          
      M3=N-2                                                            
      DO 1040 I=3,M3                                                    
      IM=I-1                                                            
      II=IM*M+1                                                         
      IE=IM*L+1                                                         
      IPM1=II+M                                                         
      IPM2=IPM1+M                                                       
      IMM1=II-M                                                         
      IMM2=IMM1-M                                                       
      Y2=X(I-2)                                                         
      Y1=X(I-1)                                                         
      X0=X(I)                                                           
      X1=X(I+1)                                                         
      X2=X(I+2)                                                         
      G2=Y1-Y2                                                          
      G1=X0-Y1                                                          
      H0=X1-X0                                                          
      H1=X2-X1                                                          
      H21=G2+G1                                                         
      H10=G1+H0                                                         
      H01=H0+H1                                                         
      H210=H21+H0                                                       
      H101=H10+H1                                                       
      H2101=H210+H1                                                     
      AM2=Y(IMM2)/(G2*H21*H210*H2101)                                   
      AM1=-Y(IMM1)/(G2*G1*H10*H101)                                     
      A0=Y(II)/(H21*G1*H0*H01)                                          
      AP1=-Y(IPM1)/(H210*H10*H0*H1)                                     
      AP2=Y(IPM2)/(H2101*H101*H01*H1)                                   
      F(IE)=(AM2*G1*H0*H01+ AM1*H21*H0*H01                              
     *      +A0*(G1*H0*H01+H21*H0*H01-H21*G1*H01-H21*G1*H0)             
     *       +AP1*(-H21*G1*H01) +AP2*(-H21*G1*H0))                      
 1040 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C#NUMPAC#PHER31              REVISED ON 1984-11-30                      
C      SUBROUTINE PHER31(X,Y,M,XI,Z,ILL)             
      SUBROUTINE PHED31(X,Y,M,XI,Z,ILL)                                 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(M),XI(2),YI(2),Y1I(2),Z(2,2)                          
      YI(1)=Z(1,1)                                                      
      YI(2)=Z(2,1)                                                      
      Y1I(1)=Z(1,2)                                                     
      Y1I(2)=Z(2,2)                                                     
C      X           XCOORDINATE(INPUT)                                   
C      Y(M)        Y-VALUE(OUTPUT)                                      
C      M           ORDER OF DERIVATIVE(INPUT)                           
C      XI          X-RANGE(INPUT)                                       
C      YI           FUNCTION VALUE ON XI(INPUT)                         
C      Y1I         DERIVATIVE ON XI(INPUT)                              
C      ILL         OUTPUT)                                              
      ILL=0                                             
      EPS=1.D-6                                         
      X1=XI(1)                                          
      X2=XI(2)                                          
      H=X2-X1                                           
      IF(H) 1000 , 1060 , 1010                          
 1000 X1=XI(2)                                          
      X2=XI(1)                                          
      H=-H                                              
 1010 T=(X-X1)/H                                        
      IF(ABS(T).LE.EPS) T=0.0D0
      IF(ABS(T-1.D0).LE.EPS) T=1.0D0
      IF(T.LT.0.0d0.OR.T.GT.1.0d0) GO TO 1060      
      T2=T*T                                            
      T1=T-1.d0                                    
      M1=M+1                                            
      DO 1050 K=1,M1                                    
      GO TO ( 1020 , 1030 , 1040 ),K                    
 1020 P=1.d0+T2*(2.d0*T-3.d0)                                 
      Q=1.d0-P                                            
      R=T*T1**2                                         
      S=T2*T1                                           
      GO TO 1050                                        
 1030 P=6.d0*(T2-T)                                       
      Q=-P                                              
      R=T1*(3.d0*T-1.d0)                                    
      S=T*(3.d0*T-2.d0)                                     
      GO TO 1050                                        
 1040 P=-6.d0+12.d0*T                                       
      Q=-P                                              
      R=6.d0*T-4.d0                                         
      S=6.d0*T-2.d0                                         
 1050 Y(K)=(YI(1)*P+YI(2)*Q+(Y1I(1)*R+Y1I(2)*S)*H)/(H**(K-1)) 
      RETURN                                                            
 1060 ILL=30000                                                         
      RETURN                                                            
      END                                                               
C#NUMPAC#FFTB                REVISED ON 1984-11-30                      
      SUBROUTINE FFTB(A,M,INV,ILL)                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(1),DCC(32),DSS(32)                                    
      DATA DCC( 1)/ 0.707106781186547531D+00/                           
      DATA DCC( 2)/ 0.923879532511286752D+00/                           
      DATA DCC( 3)/ 0.980785280403230444D+00/                           
      DATA DCC( 4)/ 0.995184726672196887D+00/                           
      DATA DCC( 5)/ 0.998795456205172391D+00/                           
      DATA DCC( 6)/ 0.999698818696204222D+00/                           
      DATA DCC( 7)/ 0.999924701839144545D+00/                           
      DATA DCC( 8)/ 0.999981175282601137D+00/                           
      DATA DCC( 9)/ 0.999995293809576177D+00/                           
      DATA DCC(10)/ 0.999998823451701907D+00/                           
      DATA DCC(11)/ 0.999999705862882213D+00/                           
      DATA DCC(12)/ 0.999999926465717850D+00/                           
      DATA DCC(13)/ 0.999999981616429293D+00/                           
      DATA DCC(14)/ 0.999999995404107320D+00/                           
      DATA DCC(15)/ 0.999999998851026833D+00/                           
      DATA DCC(16)/ 0.999999999712756701D+00/                           
      DATA DCC(17)/ 0.999999999928189179D+00/                           
      DATA DCC(18)/ 0.999999999982047291D+00/                           
      DATA DCC(19)/ 0.999999999995511826D+00/                           
      DATA DCC(20)/ 0.999999999998877953D+00/                           
      DATA DCC(21)/ 0.999999999999719488D+00/                           
      DATA DCC(22)/ 0.999999999999929876D+00/                           
      DATA DCC(23)/ 0.999999999999982472D+00/                           
      DATA DCC(24)/ 0.999999999999995615D+00/                           
      DATA DCC(25)/ 0.999999999999998904D+00/                           
      DATA DCC(26)/ 0.999999999999999722D+00/                           
      DATA DCC(27)/ 0.999999999999999931D+00/                           
      DATA DCC(28)/ 0.999999999999999986D+00/                           
      DATA DCC(29)/ 0.100000000000000000D+01/                           
      DATA DCC(30)/ 0.100000000000000000D+01/                           
      DATA DCC(31)/ 0.100000000000000000D+01/                           
      DATA DCC(32)/ 0.100000000000000000D+01/                           
      DATA DSS( 1)/ 0.707106781186547531D+00/                           
      DATA DSS( 2)/ 0.382683432365089768D+00/                           
      DATA DSS( 3)/ 0.195090322016128262D+00/                           
      DATA DSS( 4)/ 0.980171403295606036D-01/                           
      DATA DSS( 5)/ 0.490676743274180141D-01/                           
      DATA DSS( 6)/ 0.245412285229122881D-01/                           
      DATA DSS( 7)/ 0.122715382857199263D-01/                           
      DATA DSS( 8)/ 0.613588464915447527D-02/                           
      DATA DSS( 9)/ 0.306795676296597625D-02/                           
      DATA DSS(10)/ 0.153398018628476561D-02/                           
      DATA DSS(11)/ 0.766990318742704540D-03/                           
      DATA DSS(12)/ 0.383495187571395563D-03/                           
      DATA DSS(13)/ 0.191747597310703308D-03/                           
      DATA DSS(14)/ 0.958737990959773447D-04/                           
      DATA DSS(15)/ 0.479368996030668847D-04/                           
      DATA DSS(16)/ 0.239684498084182193D-04/                           
      DATA DSS(17)/ 0.119842249050697064D-04/                           
      DATA DSS(18)/ 0.599211245264242774D-05/                           
      DATA DSS(19)/ 0.299605622633466084D-05/                           
      DATA DSS(20)/ 0.149802811316901114D-05/                           
      DATA DSS(21)/ 0.749014056584715715D-06/                           
      DATA DSS(22)/ 0.374507028292384129D-06/                           
      DATA DSS(23)/ 0.187253514146195347D-06/                           
      DATA DSS(24)/ 0.936267570730980836D-07/                           
      DATA DSS(25)/ 0.468133785365490931D-07/                           
      DATA DSS(26)/ 0.234066892682745532D-07/                           
      DATA DSS(27)/ 0.117033446341372770D-07/                           
      DATA DSS(28)/ 0.585167231706863850D-08/                           
      DATA DSS(29)/ 0.292583615853431935D-08/                           
      DATA DSS(30)/ 0.146291807926715968D-08/                           
      DATA DSS(31)/ 0.731459039633579864D-09/                           
      DATA DSS(32)/ 0.365729519816789906D-09/                           
      IF(M.LT.2) GO TO 70                                               
      NN=2**(M+1)                                                       
      L=NN                                                              
      G=INV+INV-1                                                       
      DO 40 I=2,M,2                                                     
      L4=L/4                                                            
      DO 10 K0=2,NN,L                                                   
      K1=K0+L4                                                          
      K2=K1+L4                                                          
      K3=K2+L4                                                          
      T1=A(K0-1)+A(K2-1)                                                
      U1=A(K0-1)-A(K2-1)                                                
      T2=A(K0)+A(K2)                                                    
      U2=A(K0)-A(K2)                                                    
      V1=A(K1-1)+A(K3-1)                                                
      W2=(A(K1-1)-A(K3-1))*G                                            
      V2=A(K1)+A(K3)                                                    
      W1=(A(K3)-A(K1))*G                                                
      A(K0-1)=T1+V1                                                     
      A(K0)=T2+V2                                                       
      A(K1-1)=T1-V1                                                     
      A(K1)=T2-V2                                                       
      A(K2-1)=U1+W1                                                     
      A(K2)=U2+W2                                                       
      A(K3-1)=U1-W1                                                     
   10 A(K3)=U2-W2                                                       
      IF(L4.LT.4) GO TO 40                                              
      DC=DCC(M-I)                                                       
      DS=DSS(M-I)*G                                                     
      CC=DC                                                             
      SS=DS                                                             
      PP=SS+SS                                                          
      CO=1.D0                                                           
      SO=0.D0                                                           
      DO 30 J=4,L4,2                                                    
      C1=CC                                                             
      S1=SS                                                             
      P=S1+S1                                                           
      C2=(C1-S1)*(C1+S1)                                                
      S2=C1*P                                                           
      C3=C1-P*S2                                                        
      S3=P*C2+S1                                                        
      DO 20 K0=J,NN,L                                                   
      K1=K0+L4                                                          
      K2=K1+L4                                                          
      K3=K2+L4                                                          
      T1=A(K0-1)+A(K2-1)                                                
      U1=A(K0-1)-A(K2-1)                                                
      T2=A(K0)+A(K2)                                                    
      U2=A(K0)-A(K2)                                                    
      V1=A(K1-1)+A(K3-1)                                                
      W2=(A(K1-1)-A(K3-1))*G                                            
      V2=A(K1)+A(K3)                                                    
      W1=(A(K3)-A(K1))*G                                                
      A(K0-1)=T1+V1                                                     
      A(K0)=T2+V2                                                       
      A(K1-1)=(T1-V1)*C2-(T2-V2)*S2                                     
      A(K1)=(T1-V1)*S2+(T2-V2)*C2                                       
      A(K2-1)=(U1+W1)*C1-(U2+W2)*S1                                     
      A(K2)=(U1+W1)*S1+(U2+W2)*C1                                       
      A(K3-1)=(U1-W1)*C3-(U2-W2)*S3                                     
   20 A(K3)=(U1-W1)*S3+(U2-W2)*C3                                       
      CN=CO-PP*SS                                                       
      SN=PP*CC+SO                                                       
      CO=CC                                                             
      CC=CN                                                             
      SO=SS                                                             
   30 SS=SN                                                             
   40 L=L4                                                              
      IF(L.EQ.2) GO TO 60                                               
      DO 50 K0=2,NN,4                                                   
      T1=A(K0-1)                                                        
      T2=A(K0)                                                          
      A(K0-1)=T1+A(K0+1)                                                
      A(K0)=T2+A(K0+2)                                                  
      A(K0+1)=T1-A(K0+1)                                                
   50 A(K0+2)=T2-A(K0+2)                                                
   60 CALL BITRVB(A,M,ILL)                                              
      RETURN                                                            
   70 ILL=30000                                                         
      RETURN                                                            
      END                                                               
C#NUMPAC#BITRVB              REVISED ON 1984-11-30                      
      SUBROUTINE BITRVB(A,L,ICON)                                       
      COMPLEX*16 A,W                                                    
      DIMENSION A(1),ITEST(20),INC(20)                                  
      IF(L.LE.0.OR.L.GT.23) GO TO 8000                                  
      NN=0                                                              
      NR=0                                                              
      I=L-1                                                             
      M=2**I                                                            
      K=2                                                               
      IF(I-2) 60,30,10                                                  
   10 I=I-1                                                             
      ITEST(I-1)=M-K                                                    
      K=K+K                                                             
      INC(I-1)=K-ITEST(I-1)                                             
      IF(I-3) 30,10,10                                                  
   20 NR=INC(I)+NR                                                      
   30 MR=M+NR                                                           
      IF(NR-NN) 50,50,40                                                
   40 W=A(NN+1)                                                         
      A(NN+1)=A(NR+1)                                                   
      A(NR+1)=W                                                         
      MN=M+NN                                                           
      W=A(MN+2)                                                         
      A(MN+2)=A(MR+2)                                                   
      A(MR+2)=W                                                         
   50 NN=NN+2                                                           
      W=A(NN)                                                           
      A(NN)=A(MR+1)                                                     
      A(MR+1)=W                                                         
      NR=K+NR                                                           
   60 MR=M+NR                                                           
      IF(NR-NN) 80,80,70                                                
   70 W=A(NN+1)                                                         
      A(NN+1)=A(NR+1)                                                   
      A(NR+1)=W                                                         
      MN=M+NN                                                           
      W=A(MN+2)                                                         
      A(MN+2)=A(MR+2)                                                   
      A(MR+2)=W                                                         
   80 NN=NN+2                                                           
      W=A(NN)                                                           
      A(NN)=A(MR+1)                                                     
      A(MR+1)=W                                                         
      I=1                                                               
      IF(NN-M) 100,110,110                                              
   90 I=I+1                                                             
  100 IF(NR-ITEST(I)) 20,90,90                                          
  110 ICON=0                                                            
      RETURN                                                            
 8000 ICON=30000                                                        
      RETURN                                                            
      END                                                               
C#NUMPAC#MINVD               REVISED ON 1984-11-30                      
      SUBROUTINE MINVD(A,KA,N,EPS,ILL)                                  
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(KA,N)                                                 
      INTEGER*2 MX(1000)                                                
      IF(N.LT.1.OR.N.GT.1000.OR.N.GT.KA.OR.EPS.LE.0.D0) GO TO 250         
C-----LU DECOMPOSITION--------------------------------------------------
      NM1=N-1                                                           
      DO 90 J=1,N                                                       
      IF(J.EQ.1) GO TO 30                                               
      JM1=J-1                                                           
      DO 20 I=1,JM1                                                     
      M=MX(I)                                                           
      S=A(M,J)                                                          
      A(M,J)=A(I,J)                                                     
      IF(I.EQ.1) GO TO 20                                               
      IM1=I-1                                                           
      DO 10 K=1,IM1                                                     
   10 S=A(I,K)*A(K,J)+S                                                 
   20 A(I,J)=S                                                          
   30 AM=0.D0                                                             
      DO 60 I=J,N                                                       
      S=A(I,J)                                                          
      IF(J.EQ.1) GO TO 50                                               
      DO 40 K=1,JM1                                                     
   40 S=A(I,K)*A(K,J)+S                                                 
      A(I,J)=S                                                          
   50 AA=DABS(S)                                                        
      IF(AA.LE.AM) GO TO 60                                             
      AM=AA                                                             
      M=I                                                               
   60 CONTINUE                                                          
      IF(AM.LT.EPS) GO TO 240                                           
      MX(J)=M                                                           
      IF(M.EQ.J) GO TO 80                                               
      DO 70 K=1,J                                                       
      W=A(M,K)                                                          
      A(M,K)=A(J,K)                                                     
   70 A(J,K)=W                                                          
   80 IF(J.EQ.N) GO TO 100                                              
      JP1=J+1                                                           
      W=-A(J,J)                                                         
      DO 90 I=JP1,N                                                     
   90 A(I,J)=A(I,J)/W                                                   
  100 IF(N.LE.2) GO TO 130                                              
C-----INPLACE INVERSION OF L-COMPONENT----------------------------------
      DO 120 I=3,N                                                      
      IM1=I-1                                                           
      IM2=I-2                                                           
      DO 120 J=1,IM2                                                    
      S=A(I,J)                                                          
      JP1=J+1                                                           
      DO 110 K=JP1,IM1                                                  
  110 S=A(I,K)*A(K,J)+S                                                 
  120 A(I,J)=S                                                          
C-----INPLACE INVERSION OF U-COMPONENT----------------------------------
  130 A(1,1)=1./A(1,1)                                                  
      IF(N.EQ.1) GO TO 230                                              
      DO 150 J=2,N                                                      
      A(J,J)=1./A(J,J)                                                  
      P=-A(J,J)                                                         
      JM1=J-1                                                           
      DO 150 I=1,JM1                                                    
      S=0.                                                              
      DO 140 K=I,JM1                                                    
  140 S=A(I,K)*A(K,J)+S                                                 
  150 A(I,J)=S*P                                                        
C-----INPLACE MULTIPLICATION OF L AND U COMPONENT-----------------------
      DO 190 J=1,NM1                                                    
      JP1=J+1                                                           
      DO 170 I=1,J                                                      
      S=A(I,J)                                                          
      DO 160 K=JP1,N                                                    
  160 S=A(I,K)*A(K,J)+S                                                 
  170 A(I,J)=S                                                          
      DO 190 I=JP1,N                                                    
      S=0.                                                              
      DO 180 K=I,N                                                      
  180 S=A(I,K)*A(K,J)+S                                                 
  190 A(I,J)=S                                                          
C------INTERCHANGE OF COLUMNS-------------------------------------------
      J=NM1                                                             
  200 M=MX(J)                                                           
      IF(M.EQ.J) GO TO 220                                              
      DO 210 I=1,N                                                      
      W=A(I,M)                                                          
      A(I,M)=A(I,J)                                                     
  210 A(I,J)=W                                                          
  220 J=J-1                                                             
      IF(J.GE.1) GO TO 200                                              
  230 ILL=0                                                             
      RETURN                                                            
  240 ILL=J                                                             
      RETURN                                                            
  250 ILL=30000                                                         
      RETURN                                                            
      END                                                               
C#NUMPAC#TRIGQD              REVISED ON 1984-11-30                      
      SUBROUTINE TRIGQD(C,M,ICON)                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C     TABLE OF TRIGONOMETRIC FUNCTION ARRANGED IN BIT REVERSE ORDER.    
C     AUTHOR...TATSUO TORII (DEPARTMENT OF INFORMATION SCIENCE ,        
C                            FACULTY OF ENGINEERING , NAGOYA UNIVERSITY)
C     DATE...1978.12.15                                                 
      DIMENSION C(*)                                                    
      ICON=30000                                                        
      IF(M.LE.0) RETURN                                                 
      P=0.5                                                             
      C(1)=DSQRT(P)                                                     
      C(2)=DSQRT((C(1)+1.0)*0.5)                                        
      C(3)=C(1)*0.5/C(2)                                                
      C(4)=DSQRT((C(2)+1.0)*0.5)                                        
      C(5)=C(3)*0.5/C(4)                                                
      C(6)=(C(4)+C(5))*C(1)                                             
      C(7)=(C(4)-C(5))*C(1)                                             
      ICON=0                                                            
      IF(M.LE.3) RETURN                                                 
      N0=4                                                              
      N1=8                                                              
      CA=C(4)                                                           
      SA=C(5)                                                           
      DO 100 K=4,M                                                      
      CA=DSQRT((1.0+CA)*0.5)                                            
      SA=SA*0.5/CA                                                      
      C(N1)=CA                                                          
      C(N1+1)=SA                                                        
      NM=N0+N1                                                          
      C(NM)=C(N0)*CA-C(N0+1)*SA                                         
      C(NM+1)=C(N0+1)*CA+C(N0)*SA                                       
      J1=N1+2                                                           
      J2=NM+2                                                           
      LP=N0+2                                                           
      MP=N1-2                                                           
      DO 200 J=LP,MP,2                                                  
      P=C(J)*CA                                                         
      Q=C(J)*SA                                                         
      R=C(J+1)*CA                                                       
      S=C(J+1)*SA                                                       
      C(J1)=P+S                                                         
      C(J1+1)=R-Q                                                       
      C(J2)=P-S                                                         
      C(J2+1)=R+Q                                                       
      J1=J1+2                                                           
      J2=J2+2                                                           
  200 CONTINUE                                                          
      N0=N1                                                             
      N1=N1+N1                                                          
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C#NUMPAC#FSINTD              REVISED ON 1984-11-30                      
      SUBROUTINE FSINTD(X,MX,W,MW,ICON)                                 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C     FAST FOURIER  SINE TRANSFORM BASED ON THE TRAPEZOIDAL RULE.       
C     AUTHOR...TATSUO TORII (DEPARTMENT OF INFORMATION SCIENCE ,        
C                            FACULTY OF ENGINEERING , NAGOYA UNIVERSITY)
C     DATE...1978.12.15                                                 
      DIMENSION X(*),W(*)                                   
C     PARAMETER CHECK                                                   
      ICON=30000                                                        
      MQ=MX-1                                                           
      IF(MQ.GT.MW.OR.MX.LE.0) RETURN                                    
      IF(X(1).NE.0.0) RETURN                                            
C     INITIALIZE                                                        
      CALL BITRVD(X,MX,ICON)                                            
      IF(ICON.NE.0) RETURN                                              
      IF(MX.EQ.1) RETURN                                                
      N1=2                                                              
      N2=4                                                              
      DO 10 M=1,MQ                                                      
      CALL FSINMD(X,M,N1,W, MW,ICON)                                    
      IF(ICON.NE.0) RETURN                                              
      IB=N2                                                             
      DO 20 IA=2,N1                                                     
      P=X(IB)-X(IA)                                                     
      X(IA)=X(IA)+X(IB)                                                 
      X(IB)=P                                                           
      IB=IB-1                                                           
   20 CONTINUE                                                          
      N1=N2                                                             
      N2=N2+N2                                                          
   10 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C#NUMPAC#FCOSTD              REVISED ON 1984-11-30                      
      SUBROUTINE FCOSTD(X,MX,W,MW,ICON)                                 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C     FAST FOURIER TRANSFORM BASED ON THE TRAPEZOIDAL RULE.             
C     AUTHOR...TATSUO TORII (DEPARTMENT OF INFORMATION SCIENCE ,        
C                            FACULTY OF ENGINEERING , NAGOYA UNIVERSITY)
C     DATE...1978.12.15                                                 
      DIMENSION X(*),W(*)                                               
C     PARAMETER CHECK                                                   
      ICON=30000                                                        
      MQ=MX-1                                                           
      IF(MQ.GT.MW.OR.MX.LE.0) RETURN                                    
C     INITIALIZE                                                        
      CALL BITRVD(X,MX,ICON)                                            
      IF(ICON.NE.0) RETURN                                              
      NP1=2**MX+1                                                       
      P=(X(1)+X(NP1))*0.5                                               
      Q=(X(1)-X(NP1))*0.5                                               
      X(1)=P+X(2)                                                       
      X(NP1)=P-X(2)                                                     
      X(2)=Q                                                            
      IF(MX.EQ.1) RETURN                                                
      CALL FCOSMD(X,1,2,W,MW,ICON)                                      
      IF(ICON.NE.0) RETURN                                              
      P=X(3)                                                            
      X(3)=X(NP1)                                                       
      X(NP1)=X(1)-P                                                     
      X(1)=X(1)+P                                                       
      Q=X(2)-X(4)                                                       
      X(2)=X(2)+X(4)                                                    
      X(4)=Q                                                            
      IF(MX.EQ.2) RETURN                                                
      N0=2                                                              
      N1=4                                                              
      N2=8                                                              
      DO 10 M=2,MQ                                                      
      CALL FCOSMD(X,M ,N1,W,MW,ICON)                                    
      IF(ICON.NE.0) RETURN                                              
      JA=N1+1                                                           
      P=X(JA)                                                           
      X(JA)=X(NP1)                                                      
      X(NP1)=X(1)-P                                                     
      X(1)=X(1)+P                                                       
      JA=JA+1                                                           
      IB=N1                                                             
      JB=N2                                                             
      DO 20 IA=2,N0                                                     
      P=X(IA)-X(JA)                                                     
      Q=X(IB)-X(JB)                                                     
      X(IA)=X(IA)+X(JA)                                                 
      X(IB)=X(IB)+X(JB)                                                 
      X(JA)=Q                                                           
      X(JB)=P                                                           
      IB=IB-1                                                           
      JA=JA+1                                                           
      JB=JB-1                                                           
   20 CONTINUE                                                          
      P=X(IB)-X(JB)                                                     
      X(IB)=X(IB)+X(JB)                                                 
      X(JB)=P                                                           
      N0=N1                                                             
      N1=N2                                                             
      N2=N2+N2                                                          
   10 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C#NUMPAC#BITRVD              REVISED ON 1984-11-30                      
      SUBROUTINE BITRVD(A,L,ICON)                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(*),ITEST(20),INC(20)                                  
      IF(L.LE.0.OR.L.GT.23) GO TO 8000                                  
      NN=0                                                              
      NR=0                                                              
      I=L-1                                                             
      M=2**I                                                            
      K=2                                                               
      IF(I-2) 60,30,10                                                  
   10 I=I-1                                                             
      ITEST(I-1)=M-K                                                    
      K=K+K                                                             
      INC(I-1)=K-ITEST(I-1)                                             
      IF(I-3) 30,10,10                                                  
   20 NR=INC(I)+NR                                                      
   30 MR=M+NR                                                           
      IF(NR-NN) 50,50,40                                                
   40 W=A(NN+1)                                                         
      A(NN+1)=A(NR+1)                                                   
      A(NR+1)=W                                                         
      MN=M+NN                                                           
      W=A(MN+2)                                                         
      A(MN+2)=A(MR+2)                                                   
      A(MR+2)=W                                                         
   50 NN=NN+2                                                           
      W=A(NN)                                                           
      A(NN)=A(MR+1)                                                     
      A(MR+1)=W                                                         
      NR=K+NR                                                           
   60 MR=M+NR                                                           
      IF(NR-NN) 80,80,70                                                
   70 W=A(NN+1)                                                         
      A(NN+1)=A(NR+1)                                                   
      A(NR+1)=W                                                         
      MN=M+NN                                                           
      W=A(MN+2)                                                         
      A(MN+2)=A(MR+2)                                                   
      A(MR+2)=W                                                         
   80 NN=NN+2                                                           
      W=A(NN)                                                           
      A(NN)=A(MR+1)                                                     
      A(MR+1)=W                                                         
      I=1                                                               
      IF(NN-M) 100,110,110                                              
   90 I=I+1                                                             
  100 IF(NR-ITEST(I)) 20,90,90                                          
  110 ICON=0                                                            
      RETURN                                                            
 8000 ICON=30000                                                        
      RETURN                                                            
      END                                                               
C#NUMPAC#MINVB               REVISED ON 1984-11-30                      
      SUBROUTINE MINVB(A,KA,N,EPS,ILL)                                  
      REAL*8 EPS,AM,AA                                                  
      COMPLEX*16 A(KA,N),S,W,P                                          
      INTEGER*2 MX(1000)                                                
      IF(N.LT.1.OR.N.GT.1000.OR.N.GT.KA.OR.EPS.LE.0.) GO TO 250         
C-----LU DECOMPOSITION--------------------------------------------------
      NM1=N-1                                                           
      DO 90 J=1,N                                                       
      IF(J.EQ.1) GO TO 30                                               
      JM1=J-1                                                           
      DO 20 I=1,JM1                                                     
      M=MX(I)                                                           
      S=A(M,J)                                                          
      A(M,J)=A(I,J)                                                     
      IF(I.EQ.1) GO TO 20                                               
      IM1=I-1                                                           
      DO 10 K=1,IM1                                                     
   10 S=A(I,K)*A(K,J)+S                                                 
   20 A(I,J)=S                                                          
   30 AM=0.                                                             
      DO 60 I=J,N                                                       
      S=A(I,J)                                                          
      IF(J.EQ.1) GO TO 50                                               
      DO 40 K=1,JM1                                                     
   40 S=A(I,K)*A(K,J)+S                                                 
      A(I,J)=S                                                          
   50 AA=CDABS(S)                                                       
      IF(AA.LE.AM) GO TO 60                                             
      AM=AA                                                             
      M=I                                                               
   60 CONTINUE                                                          
      IF(AM.LT.EPS) GO TO 240                                           
      MX(J)=M                                                           
      IF(M.EQ.J) GO TO 80                                               
      DO 70 K=1,J                                                       
      W=A(M,K)                                                          
      A(M,K)=A(J,K)                                                     
   70 A(J,K)=W                                                          
   80 IF(J.EQ.N) GO TO 100                                              
      JP1=J+1                                                           
      W=-A(J,J)                                                         
      DO 90 I=JP1,N                                                     
   90 A(I,J)=A(I,J)/W                                                   
  100 IF(N.LE.2) GO TO 130                                              
C-----INPLACE INVERSION OF L-COMPONENT----------------------------------
      DO 120 I=3,N                                                      
      IM1=I-1                                                           
      IM2=I-2                                                           
      DO 120 J=1,IM2                                                    
      S=A(I,J)                                                          
      JP1=J+1                                                           
      DO 110 K=JP1,IM1                                                  
  110 S=A(I,K)*A(K,J)+S                                                 
  120 A(I,J)=S                                                          
C-----INPLACE INVERSION OF U-COMPONENT----------------------------------
  130 A(1,1)=1./A(1,1)                                                  
      IF(N.EQ.1) GO TO 230                                              
      DO 150 J=2,N                                                      
      A(J,J)=1./A(J,J)                                                  
      P=-A(J,J)                                                         
      JM1=J-1                                                           
      DO 150 I=1,JM1                                                    
      S=0.                                                              
      DO 140 K=I,JM1                                                    
  140 S=A(I,K)*A(K,J)+S                                                 
  150 A(I,J)=S*P                                                        
C-----INPLACE MULTIPLICATION OF L AND U COMPONENT-----------------------
      DO 190 J=1,NM1                                                    
      JP1=J+1                                                           
      DO 170 I=1,J                                                      
      S=A(I,J)                                                          
      DO 160 K=JP1,N                                                    
  160 S=A(I,K)*A(K,J)+S                                                 
  170 A(I,J)=S                                                          
      DO 190 I=JP1,N                                                    
      S=0.                                                              
      DO 180 K=I,N                                                      
  180 S=A(I,K)*A(K,J)+S                                                 
  190 A(I,J)=S                                                          
C------INTERCHANGE OF COLUMNS-------------------------------------------
      J=NM1                                                             
  200 M=MX(J)                                                           
      IF(M.EQ.J) GO TO 220                                              
      DO 210 I=1,N                                                      
      W=A(I,M)                                                          
      A(I,M)=A(I,J)                                                     
  210 A(I,J)=W                                                          
  220 J=J-1                                                             
      IF(J.GE.1) GO TO 200                                              
  230 ILL=0                                                             
      RETURN                                                            
  240 ILL=J                                                             
      RETURN                                                            
  250 ILL=30000                                                         
      RETURN                                                            
      END                                                               
C#NUMPAC#MNORMB              REVISED ON 1984-11-30                      
      SUBROUTINE MNORMB(A,KA,N,M,S,ILL)                                 
C                                                80. 6.18               
      REAL*8 S(N),D,AM,F,CDABS1                                         
      COMPLEX*16 A(KA,M),Z                                     
      CDABS1(Z)=DABS(DREAL(Z))+DABS(DIMAG(Z))                           
      IF(N.LT.2.OR.N.GT.KA.OR.N.GT.M) GO TO 40                          
      F=1.0D0/DLOG(2.0D0)                                               
      DO 20 I=1,N                                                       
      AM=0.                                                             
      DO 10 J=1,N                                                       
   10 AM=DMAX1(CDABS1(A(I,J)),AM)                                       
      IF(AM.EQ.0.) GO TO 30                                             
      NPS=DLOG(AM)*F                                                    
      S(I)=2.0D0**NPS                                                   
      D=1.D0/S(I)                                                       
      DO 20 J=1,M                                                       
   20 A(I,J)=A(I,J)*D                                                   
      ILL=0                                                             
      RETURN                                                            
   30 ILL=I                                                             
      RETURN                                                            
   40 ILL=30000                                                         
      RETURN                                                            
      END                                                               
C#NUMPAC#FSINMD              REVISED ON 1984-11-30                      
      SUBROUTINE FSINMD(X,MX,LX,W,MW,ICON)                              
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C     FAST FOURIER  SINE TRANSFORM BASED ON THE MIDPOINT RULE.          
C     AUTHOR...TATSUO TORII (DEPARTMENT OF INFORMATION SCIENCE ,        
C                            FACULTY OF ENGINEERING , NAGOYA UNIVERSITY)
C     DATE...1978.12.15                                                 
      DIMENSION X(*),W(*) 
      DATA CHECK /0.0D0/                                                
      ICON=30000                                                        
      IF(MX.LE.0) RETURN                                                
      IF(MX.GT.MW) RETURN                                               
      IF(CHECK.EQ.0.0) CHECK=W(1)*0.5/W(2)                              
      IF(W(3).NE. CHECK) RETURN                                         
C     INITIALIZE                                                        
      ICON=0                                                            
      IF(MX.LE.2) GO TO 100                                             
      NQ=2**(MX-2)                                                      
      NH=NQ+NQ                                                          
      NHALF=LX+NH                                                       
      N=NH+NH                                                           
      LP1=LX+1                                                          
      NLX=N+LX                                                          
C     PRE-PROCESDSING                                                   
      J2=NLX-3                                                          
      JA=NH                                                             
      JB=N                                                              
      DO 10 J1=LP1,NHALF,4                                              
      P1=X(J1)-X(J2+3)                                                  
      Q1=X(J1)+X(J2+3)                                                  
      P2=X(J1+1)-X(J2+2)                                                
      Q2=-X(J1+1)-X(J2+2)                                               
      P3=X(J2)-X(J1+3)                                                  
      Q3=X(J2)+X(J1+3)                                                  
      P4=X(J2+1)-X(J1+2)                                                
      Q4=-X(J2+1)-X(J1+2)                                               
      X(J1)=P1+P2                                                       
      X(J1+2)=P1-P2                                                     
      X(J2)=P3+P4                                                       
      X(J2+2)=P3-P4                                                     
      X(J1+1)=Q1*W(JA)+Q2*W(JA+1)                                       
      X(J1+3)=Q2*W(JA)-Q1*W(JA+1)                                       
      JA=JA+2                                                           
      JB=JB-2                                                           
      X(J2+1)=Q3*W(JB)+Q4*W(JB+1)                                       
      X(J2+3)=Q4*W(JB)-Q3*W(JB+1)                                       
      J2=J2-4                                                           
   10 CONTINUE                                                          
      N0=8                                                              
      N1=4                                                              
      N2=2                                                              
      NR=NQ/2                                                           
   20 CONTINUE                                                          
      NN=NH                                                             
      NH=NQ                                                             
      NQ=NR                                                             
      NR=NR/2                                                           
      KP=N1-1                                                           
      KQ=N0-1                                                           
      NS=NQ+NR                                                          
      JQ=NR                                                             
      IS=-1                                                             
      JA=NH                                                             
      JB=NN                                                             
C     STAGE 1                                                           
      IF(N0.EQ.N) GO TO 70                                              
      I0=NLX-KQ                                                         
      DO 30 J0=LP1,NHALF,N0                                             
      J1=J0+N1                                                          
      I1=I0+N1                                                          
      P=X(J0)+X(I1)                                                     
      Q=X(I0)+X(J1)                                                     
      X(J0)=X(J0)-X(I1)                                                 
      X(I0)=X(I0)-X(J1)                                                 
      X(J1)=P                                                           
      X(I1)=Q                                                           
      LJQ=NS+JQ*IS                                                      
      CQ=W(LJQ)                                                         
      JQ=JQ-1                                                           
      IS=-IS                                                            
      LJQ1=LJQ+IS                                                       
      SQ=-W(LJQ1)                                                       
      KB=KP                                                             
      MP=N2-1                                                           
C     STAGE 2                                                           
      DO 40 KA=1,MP                                                     
      J0KA=J0+KA                                                        
      I0KA=I0+KA                                                        
      J1KA=J1+KA                                                        
      I1KA=I1+KA                                                        
      J0KB=J0+KB                                                        
      I0KB=I0+KB                                                        
      J1KB=J1+KB                                                        
      I1KB=I1+KB                                                        
      P2=X(J0KA)+X(I1KA)                                                
      Q2=X(J0KB)-X(I1KB)                                                
      P4=X(I0KA)+X(J1KA)                                                
      Q4=X(I0KB)-X(J1KB)                                                
      Q3=X(I0KB)+X(J1KB)                                                
      X(J0KA)=X(J0KA)-X(I1KA)                                           
      X(I0KA)=X(I0KA)-X(J1KA)                                           
      X(J1KB)=X(J0KB)+X(I1KB)                                           
      X(J0KB)=P2*CQ+Q2*SQ                                               
      X(J1KA)=P2*SQ-Q2*CQ                                               
      X(I0KB)=-P4*SQ-Q4*CQ                                              
      X(I1KA)=Q4*SQ-P4*CQ                                               
      X(I1KB)=Q3                                                        
      KB=KB-1                                                           
   40 CONTINUE                                                          
C     STAGE 3                                                           
      Q1=X(I1KA+1)*W(JA)+X(J0KA+1)*W(JA+1)                              
      X(J0KA+1)=X(J0KA+1)*W(JA)-X(I1KA+1)*W(JA+1)                       
      JA=JA+2                                                           
      JB=JB-2                                                           
      Q2=X(J1KA+1)*W(JB)+X(I0KA+1)*W(JB+1)                              
      X(I0KA+1)=X(I0KA+1)*W(JB)-X(J1KA+1)*W(JB+1)                       
      X(J1KA+1)=-Q1                                                     
      X(I1KA+1)=-Q2                                                     
      I0=I0-N0                                                          
   30 CONTINUE                                                          
      N2=N1                                                             
      N1=N0                                                             
      N0=N0+N0                                                          
      GO TO 20                                                          
   70 CONTINUE                                                          
C     POST-PROCESDSING                                                  
      J1=LP1+N1                                                         
      P=X(LP1)+X(J1)                                                    
      X(LP1)=X(LP1)-X(J1)                                               
      X(J1)=P*W(1)                                                      
      J0KB=LP1+KP                                                       
      J1KA=J1+1                                                         
      J1KB=J1+KP                                                        
      LP2=LP1+1                                                         
      MP=LP2+N2-2                                                       
      DO 80 J0KA=LP2,MP                                                 
      P2=X(J0KA)+X(J1KA)                                                
      Q2=X(J0KB)-X(J1KB)                                                
      X(J0KA)=X(J0KA)-X(J1KA)                                           
      X(J1KB)=-X(J1KB)-X(J0KB)                                          
      X(J0KB)=(P2-Q2)*W(1)                                              
      X(J1KA)=(P2+Q2)*W(1)                                              
      J0KB=J0KB-1                                                       
      J1KA=J1KA+1                                                       
      J1KB=J1KB-1                                                       
   80 CONTINUE                                                          
      Q1=X(J0KB)*W(3)+X(J1KA)*W(2)                                      
      X(J0KB)=X(J0KB)*W(2)-X(J1KA)*W(3)                                 
      X(J1KA)=Q1                                                        
      RETURN                                                            
  100 CONTINUE                                                          
      IF(MX.NE.1) GO TO 110                                             
      P=X(LX+1)+X(LX+2)                                                 
      X(LX+1)=X(LX+1)-X(LX+2)                                           
      X(LX+2)=P*W(1)                                                    
      RETURN                                                            
  110 CONTINUE                                                          
      P1=X(LX+1)-X(LX+4)                                                
      Q1=X(LX+1)+X(LX+4)                                                
      P2=X(LX+2)-X(LX+3)                                                
      Q2=X(LX+3)+X(LX+2)                                                
      X(LX+1)=P1+P2                                                     
      X(LX+3)=(P1-P2)*W(1)                                              
      X(LX+2)=Q1*W(2)-Q2*W(3)                                           
      X(LX+4)=Q1*W(3)+Q2*W(2)                                           
      RETURN                                                            
      END                                                               
C#NUMPAC#FCOSMD              REVISED ON 1984-11-30                      
      SUBROUTINE FCOSMD(X,MX,LX,W,MW,ICON)                              
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C     FAST FOURIER COSINE   TRANSFORM BASED ON THE MIDPOINT RULE.       
C     AUTHOR...TATSUO TORII (DEPARTMENT OF INFORMATION SCIENCE ,        
C                            FACULTY OF ENGINEERING , NAGOYA UNIVERSITY)
C     DATE...1978.12.15                                                 
      DIMENSION X(*),W(*)                                               
      DATA CHECK /0.0D0/                                                
      ICON=30000                                                        
      IF(MX.LE.0) RETURN                                                
      IF(MX.GT.MW) RETURN                                               
      IF(CHECK.EQ.0.0) CHECK=W(1)*0.5/W(2)                              
      IF(W(3).NE. CHECK) RETURN                                         
C     INITIALIZE                                                        
      ICON=0                                                            
      IF(MX.LE.2) GO TO 100                                             
      NQ=2**(MX-2)                                                      
      NH=NQ+NQ                                                          
      NHALF=LX+NH                                                       
      N=NH+NH                                                           
      LP1=LX+1                                                          
      NLX=N+LX                                                          
C     PRE-PROCESDSING                                                   
      J2=NLX-3                                                          
      JA=NH                                                             
      JB=N                                                              
      DO 10 J1=LP1,NHALF,4                                              
      P1=X(J1)+X(J2+3)                                                  
      Q1=X(J1)-X(J2+3)                                                  
      P2=X(J1+1)+X(J2+2)                                                
      Q2=X(J2+2)-X(J1+1)                                                
      P3=X(J2)+X(J1+3)                                                  
      Q3=X(J2)-X(J1+3)                                                  
      P4=X(J2+1)+X(J1+2)                                                
      Q4=X(J1+2)-X(J2+1)                                                
      X(J1)=P1+P2                                                       
      X(J1+2)=P1-P2                                                     
      X(J2)=P3+P4                                                       
      X(J2+2)=P3-P4                                                     
      X(J1+1)=Q1*W(JA)+Q2*W(JA+1)                                       
      X(J1+3)=Q2*W(JA)-Q1*W(JA+1)                                       
      JA=JA+2                                                           
      JB=JB-2                                                           
      X(J2+1)=Q3*W(JB)+Q4*W(JB+1)                                       
      X(J2+3)=Q4*W(JB)-Q3*W(JB+1)                                       
      J2=J2-4                                                           
   10 CONTINUE                                                          
      N0=8                                                              
      N1=4                                                              
      N2=2                                                              
      NR=NQ/2                                                           
   20 CONTINUE                                                          
      NN=NH                                                             
      NH=NQ                                                             
      NQ=NR                                                             
      NR=NR/2                                                           
      KP=N1-1                                                           
      KQ=N0-1                                                           
      NS=NQ+NR                                                          
      JQ=NR                                                             
      IS=-1                                                             
      JA=NH                                                             
      JB=NN                                                             
C     STAGE 1                                                           
      IF(N0.EQ.N) GO TO 70                                              
      I0=NLX-KQ                                                         
      DO 30 J0=LP1,NHALF,N0                                             
      J1=J0+N1                                                          
      I1=I0+N1                                                          
      P=X(J0)-X(I1)                                                     
      Q=X(I0)-X(J1)                                                     
      X(J0)=X(J0)+X(I1)                                                 
      X(I0)=X(I0)+X(J1)                                                 
      X(J1)=P                                                           
      X(I1)=Q                                                           
      LJQ=NS+JQ*IS                                                      
      CQ=W(LJQ)                                                         
      JQ=JQ-1                                                           
      IS=-IS                                                            
      LJQ1=LJQ+IS                                                       
      SQ=-W(LJQ1)                                                       
      KB=KP                                                             
      MP=N2-1                                                           
C     STAGE 2                                                           
      DO 40 KA=1,MP                                                     
      J0KA=J0+KA                                                        
      I0KA=I0+KA                                                        
      J1KA=J1+KA                                                        
      I1KA=I1+KA                                                        
      J0KB=J0+KB                                                        
      I0KB=I0+KB                                                        
      J1KB=J1+KB                                                        
      I1KB=I1+KB                                                        
      P2=X(J0KA)-X(I1KA)                                                
      Q2=X(J0KB)+X(I1KB)                                                
      P4=X(I0KA)-X(J1KA)                                                
      Q4=X(I0KB)+X(J1KB)                                                
      Q3=X(I0KB)-X(J1KB)                                                
      X(J0KA)=X(J0KA)+X(I1KA)                                           
      X(I0KA)=X(I0KA)+X(J1KA)                                           
      X(J1KB)=X(J0KB)-X(I1KB)                                           
      X(J0KB)=P2*CQ+Q2*SQ                                               
      X(J1KA)=P2*SQ-Q2*CQ                                               
      X(I0KB)=-P4*SQ-Q4*CQ                                              
      X(I1KA)=Q4*SQ-P4*CQ                                               
      X(I1KB)=Q3                                                        
      KB=KB-1                                                           
   40 CONTINUE                                                          
C     STAGE 3                                                           
      Q1=X(I1KA+1)*W(JA)-X(J0KA+1)*W(JA+1)                              
      X(J0KA+1)=X(J0KA+1)*W(JA)+X(I1KA+1)*W(JA+1)                       
      JA=JA+2                                                           
      JB=JB-2                                                           
      Q2=X(J1KA+1)*W(JB)-X(I0KA+1)*W(JB+1)                              
      X(I0KA+1)=X(I0KA+1)*W(JB)+X(J1KA+1)*W(JB+1)                       
      X(J1KA+1)=Q1                                                      
      X(I1KA+1)=Q2                                                      
      I0=I0-N0                                                          
   30 CONTINUE                                                          
      N2=N1                                                             
      N1=N0                                                             
      N0=N0+N0                                                          
      GO TO 20                                                          
   70 CONTINUE                                                          
C     POST-PROCESDSING                                                  
      J1=LP1+N1                                                         
      P=X(LP1)-X(J1)                                                    
      X(LP1)=X(LP1)+X(J1)                                               
      X(J1)=P*W(1)                                                      
      J0KB=LP1+KP                                                       
      J1KA=J1+1                                                         
      J1KB=J1+KP                                                        
      LP2=LP1+1                                                         
      MP=LP2+N2-2                                                       
      DO 80 J0KA=LP2,MP                                                 
      P2=X(J0KA)-X(J1KA)                                                
      Q2=X(J0KB)+X(J1KB)                                                
      X(J0KA)=X(J0KA)+X(J1KA)                                           
      X(J1KB)=X(J1KB)-X(J0KB)                                           
      X(J0KB)=(P2-Q2)*W(1)                                              
      X(J1KA)=(P2+Q2)*W(1)                                              
      J0KB=J0KB-1                                                       
      J1KA=J1KA+1                                                       
      J1KB=J1KB-1                                                       
   80 CONTINUE                                                          
      Q1=X(J0KB)*W(3)-X(J1KA)*W(2)                                      
      X(J0KB)=X(J0KB)*W(2)+X(J1KA)*W(3)                                 
      X(J1KA)=Q1                                                        
      RETURN                                                            
  100 CONTINUE                                                          
      IF(MX.NE.1) GO TO 110                                             
      P=X(LX+1)-X(LX+2)                                                 
      X(LX+1)=X(LX+1)+X(LX+2)                                           
      X(LX+2)=P*W(1)                                                    
      RETURN                                                            
  110 CONTINUE                                                          
      P1=X(LX+1)+X(LX+4)                                                
      Q1=X(LX+1)-X(LX+4)                                                
      P2=X(LX+2)+X(LX+3)                                                
      Q2=X(LX+3)-X(LX+2)                                                
      X(LX+1)=P1+P2                                                     
      X(LX+3)=(P1-P2)*W(1)                                              
      X(LX+2)=Q1*W(2)+Q2*W(3)                                           
      X(LX+4)=Q1*W(3)-Q2*W(2)                                           
      RETURN                                                            
      END                                                               
C#NUMPAC#WGLEGD              REVISED ON 1985-12-16                      
      SUBROUTINE WGLEGD(NP,PT,WT,ICON)                                  
C     LEGENDRE-GAUSS FORMULA                                            
c
c     PT..zero point
c     WT..Christoffel number
c     NP..n-point
c
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION X(25,50),W(25,50),PT(NP),WT(NP)                         
      DATA X( 1, 1)/  0.0                       /                       
      DATA W( 1, 1)/  0.20000000000000000000D+01/                       
      DATA X( 1, 2)/  0.57735026918962576451D+00/                       
      DATA W( 1, 2)/  0.10000000000000000000D+01/                       
      DATA X( 1, 3)/  0.77459666924148337704D+00/                       
      DATA X( 2, 3)/  0.0                       /                       
      DATA W( 1, 3)/  0.55555555555555555556D+00/                       
      DATA W( 2, 3)/  0.88888888888888888889D+00/                       
      DATA X( 1, 4)/  0.86113631159405257522D+00/                       
      DATA X( 2, 4)/  0.33998104358485626480D+00/                       
      DATA W( 1, 4)/  0.34785484513745385737D+00/                       
      DATA W( 2, 4)/  0.65214515486254614263D+00/                       
      DATA X( 1, 5)/  0.90617984593866399280D+00/                       
      DATA X( 2, 5)/  0.53846931010568309104D+00/                       
      DATA X( 3, 5)/  0.0                       /                       
      DATA W( 1, 5)/  0.23692688505618908751D+00/                       
      DATA W( 2, 5)/  0.47862867049936646804D+00/                       
      DATA W( 3, 5)/  0.56888888888888888889D+00/                       
      DATA X( 1, 6)/  0.93246951420315202781D+00/                       
      DATA X( 2, 6)/  0.66120938646626451366D+00/                       
      DATA X( 3, 6)/  0.23861918608319690863D+00/                       
      DATA W( 1, 6)/  0.17132449237917034504D+00/                       
      DATA W( 2, 6)/  0.36076157304813860757D+00/                       
      DATA W( 3, 6)/  0.46791393457269104739D+00/                       
      DATA X( 1, 7)/  0.94910791234275852453D+00/                       
      DATA X( 2, 7)/  0.74153118559939443986D+00/                       
      DATA X( 3, 7)/  0.40584515137739716691D+00/                       
      DATA X( 4, 7)/  0.0                       /                       
      DATA W( 1, 7)/  0.12948496616886969327D+00/                       
      DATA W( 2, 7)/  0.27970539148927666790D+00/                       
      DATA W( 3, 7)/  0.38183005050511894495D+00/                       
      DATA W( 4, 7)/  0.41795918367346938776D+00/                       
      DATA X( 1, 8)/  0.96028985649753623168D+00/                       
      DATA X( 2, 8)/  0.79666647741362673959D+00/                       
      DATA X( 3, 8)/  0.52553240991632898582D+00/                       
      DATA X( 4, 8)/  0.18343464249564980494D+00/                       
      DATA W( 1, 8)/  0.10122853629037625915D+00/                       
      DATA W( 2, 8)/  0.22238103445337447054D+00/                       
      DATA W( 3, 8)/  0.31370664587788728734D+00/                       
      DATA W( 4, 8)/  0.36268378337836198297D+00/                       
      DATA X( 1, 9)/  0.96816023950762608984D+00/                       
      DATA X( 2, 9)/  0.83603110732663579430D+00/                       
      DATA X( 3, 9)/  0.61337143270059039731D+00/                       
      DATA X( 4, 9)/  0.32425342340380892904D+00/                       
      DATA X( 5, 9)/  0.0                       /                       
      DATA W( 1, 9)/  0.81274388361574411972D-01/                       
      DATA W( 2, 9)/  0.18064816069485740406D+00/                       
      DATA W( 3, 9)/  0.26061069640293546232D+00/                       
      DATA W( 4, 9)/  0.31234707704000284007D+00/                       
      DATA W( 5, 9)/  0.33023935500125976316D+00/                       
      DATA X( 1,10)/  0.97390652851717172008D+00/                       
      DATA X( 2,10)/  0.86506336668898451073D+00/                       
      DATA X( 3,10)/  0.67940956829902440623D+00/                       
      DATA X( 4,10)/  0.43339539412924719080D+00/                       
      DATA X( 5,10)/  0.14887433898163121088D+00/                       
      DATA W( 1,10)/  0.66671344308688137594D-01/                       
      DATA W( 2,10)/  0.14945134915058059315D+00/                       
      DATA W( 3,10)/  0.21908636251598204400D+00/                       
      DATA W( 4,10)/  0.26926671930999635509D+00/                       
      DATA W( 5,10)/  0.29552422471475287017D+00/                       
      DATA X( 1,11)/  0.97822865814605699280D+00/                       
      DATA X( 2,11)/  0.88706259976809529908D+00/                       
      DATA X( 3,11)/  0.73015200557404932409D+00/                       
      DATA X( 4,11)/  0.51909612920681181593D+00/                       
      DATA X( 5,11)/  0.26954315595234497233D+00/                       
      DATA X( 6,11)/  0.0                       /                       
      DATA W( 1,11)/  0.55668567116173666483D-01/                       
      DATA W( 2,11)/  0.12558036946490462463D+00/                       
      DATA W( 3,11)/  0.18629021092773425143D+00/                       
      DATA W( 4,11)/  0.23319376459199047992D+00/                       
      DATA W( 5,11)/  0.26280454451024666218D+00/                       
      DATA W( 6,11)/  0.27292508677790063071D+00/                       
      DATA X( 1,12)/  0.98156063424671925069D+00/                       
      DATA X( 2,12)/  0.90411725637047485668D+00/                       
      DATA X( 3,12)/  0.76990267419430468704D+00/                       
      DATA X( 4,12)/  0.58731795428661744730D+00/                       
      DATA X( 5,12)/  0.36783149899818019375D+00/                       
      DATA X( 6,12)/  0.12523340851146891547D+00/                       
      DATA W( 1,12)/  0.47175336386511827195D-01/                       
      DATA W( 2,12)/  0.10693932599531843096D+00/                       
      DATA W( 3,12)/  0.16007832854334622633D+00/                       
      DATA W( 4,12)/  0.20316742672306592175D+00/                       
      DATA W( 5,12)/  0.23349253653835480876D+00/                       
      DATA W( 6,12)/  0.24914704581340278500D+00/                       
      DATA X( 1,13)/  0.98418305471858814947D+00/                       
      DATA X( 2,13)/  0.91759839922297796521D+00/                       
      DATA X( 3,13)/  0.80157809073330991279D+00/                       
      DATA X( 4,13)/  0.64234933944034022064D+00/                       
      DATA X( 5,13)/  0.44849275103644685288D+00/                       
      DATA X( 6,13)/  0.23045831595513479407D+00/                       
      DATA X( 7,13)/  0.0                       /                       
      DATA W( 1,13)/  0.40484004765315879520D-01/                       
      DATA W( 2,13)/  0.92121499837728447914D-01/                       
      DATA W( 3,13)/  0.13887351021978723846D+00/                       
      DATA W( 4,13)/  0.17814598076194573828D+00/                       
      DATA W( 5,13)/  0.20781604753688850231D+00/                       
      DATA W( 6,13)/  0.22628318026289723841D+00/                       
      DATA W( 7,13)/  0.23255155323087391019D+00/                       
      DATA X( 1,14)/  0.98628380869681233884D+00/                       
      DATA X( 2,14)/  0.92843488366357351734D+00/                       
      DATA X( 3,14)/  0.82720131506976499319D+00/                       
      DATA X( 4,14)/  0.68729290481168547015D+00/                       
      DATA X( 5,14)/  0.51524863635815409197D+00/                       
      DATA X( 6,14)/  0.31911236892788976044D+00/                       
      DATA X( 7,14)/  0.10805494870734366207D+00/                       
      DATA W( 1,14)/  0.35119460331751863032D-01/                       
      DATA W( 2,14)/  0.80158087159760209806D-01/                       
      DATA W( 3,14)/  0.12151857068790318469D+00/                       
      DATA W( 4,14)/  0.15720316715819353457D+00/                       
      DATA W( 5,14)/  0.18553839747793781374D+00/                       
      DATA W( 6,14)/  0.20519846372129560397D+00/                       
      DATA W( 7,14)/  0.21526385346315779020D+00/                       
      DATA X( 1,15)/  0.98799251802048542849D+00/                       
      DATA X( 2,15)/  0.93727339240070590431D+00/                       
      DATA X( 3,15)/  0.84820658341042721620D+00/                       
      DATA X( 4,15)/  0.72441773136017004742D+00/                       
      DATA X( 5,15)/  0.57097217260853884754D+00/                       
      DATA X( 6,15)/  0.39415134707756336990D+00/                       
      DATA X( 7,15)/  0.20119409399743452230D+00/                       
      DATA X( 8,15)/  0.0                       /                       
      DATA W( 1,15)/  0.30753241996117268355D-01/                       
      DATA W( 2,15)/  0.70366047488108124709D-01/                       
      DATA W( 3,15)/  0.10715922046717193501D+00/                       
      DATA W( 4,15)/  0.13957067792615431445D+00/                       
      DATA W( 5,15)/  0.16626920581699393355D+00/                       
      DATA W( 6,15)/  0.18616100001556221103D+00/                       
      DATA W( 7,15)/  0.19843148532711157646D+00/                       
      DATA W( 8,15)/  0.20257824192556127288D+00/                       
      DATA X( 1,16)/  0.98940093499164993260D+00/                       
      DATA X( 2,16)/  0.94457502307323257608D+00/                       
      DATA X( 3,16)/  0.86563120238783174388D+00/                       
      DATA X( 4,16)/  0.75540440835500303390D+00/                       
      DATA X( 5,16)/  0.61787624440264374845D+00/                       
      DATA X( 6,16)/  0.45801677765722738634D+00/                       
      DATA X( 7,16)/  0.28160355077925891323D+00/                       
      DATA X( 8,16)/  0.95012509837637440185D-01/                       
      DATA W( 1,16)/  0.27152459411754094852D-01/                       
      DATA W( 2,16)/  0.62253523938647892863D-01/                       
      DATA W( 3,16)/  0.95158511682492784810D-01/                       
      DATA W( 4,16)/  0.12462897125553387205D+00/                       
      DATA W( 5,16)/  0.14959598881657673208D+00/                       
      DATA W( 6,16)/  0.16915651939500253819D+00/                       
      DATA W( 7,16)/  0.18260341504492358887D+00/                       
      DATA W( 8,16)/  0.18945061045506849629D+00/                       
      DATA X( 1,17)/  0.99057547531441733568D+00/                       
      DATA X( 2,17)/  0.95067552176876776122D+00/                       
      DATA X( 3,17)/  0.88023915372698590212D+00/                       
      DATA X( 4,17)/  0.78151400389680140693D+00/                       
      DATA X( 5,17)/  0.65767115921669076585D+00/                       
      DATA X( 6,17)/  0.51269053708647696789D+00/                       
      DATA X( 7,17)/  0.35123176345387631530D+00/                       
      DATA X( 8,17)/  0.17848418149584785585D+00/                       
      DATA X( 9,17)/  0.0                       /                       
      DATA W( 1,17)/  0.24148302868547931960D-01/                       
      DATA W( 2,17)/  0.55459529373987201129D-01/                       
      DATA W( 3,17)/  0.85036148317179180884D-01/                       
      DATA W( 4,17)/  0.11188384719340397109D+00/                       
      DATA W( 5,17)/  0.13513636846852547329D+00/                       
      DATA W( 6,17)/  0.15404576107681028808D+00/                       
      DATA W( 7,17)/  0.16800410215645004451D+00/                       
      DATA W( 8,17)/  0.17656270536699264633D+00/                       
      DATA W( 9,17)/  0.17944647035620652546D+00/                       
      DATA X( 1,18)/  0.99156516842093094673D+00/                       
      DATA X( 2,18)/  0.95582394957139775518D+00/                       
      DATA X( 3,18)/  0.89260246649755573921D+00/                       
      DATA X( 4,18)/  0.80370495897252311568D+00/                       
      DATA X( 5,18)/  0.69168704306035320787D+00/                       
      DATA X( 6,18)/  0.55977083107394753461D+00/                       
      DATA X( 7,18)/  0.41175116146284264604D+00/                       
      DATA X( 8,18)/  0.25188622569150550959D+00/                       
      DATA X( 9,18)/  0.84775013041735301242D-01/                       
      DATA W( 1,18)/  0.21616013526483310313D-01/                       
      DATA W( 2,18)/  0.49714548894969796453D-01/                       
      DATA W( 3,18)/  0.76425730254889056529D-01/                       
      DATA W( 4,18)/  0.10094204410628716556D+00/                       
      DATA W( 5,18)/  0.12255520671147846018D+00/                       
      DATA W( 6,18)/  0.14064291467065065120D+00/                       
      DATA W( 7,18)/  0.15468467512626524493D+00/                       
      DATA W( 8,18)/  0.16427648374583272299D+00/                       
      DATA W( 9,18)/  0.16914238296314359184D+00/                       
      DATA X( 1,19)/  0.99240684384358440319D+00/                       
      DATA X( 2,19)/  0.96020815213483003085D+00/                       
      DATA X( 3,19)/  0.90315590361481790164D+00/                       
      DATA X( 4,19)/  0.82271465653714282498D+00/                       
      DATA X( 5,19)/  0.72096617733522937862D+00/                       
      DATA X( 6,19)/  0.60054530466168102347D+00/                       
      DATA X( 7,19)/  0.46457074137596094572D+00/                       
      DATA X( 8,19)/  0.31656409996362983199D+00/                       
      DATA X( 9,19)/  0.16035864564022537587D+00/                       
      DATA X(10,19)/  0.0                       /                       
      DATA W( 1,19)/  0.19461788229726477036D-01/                       
      DATA W( 2,19)/  0.44814226765699600333D-01/                       
      DATA W( 3,19)/  0.69044542737641226581D-01/                       
      DATA W( 4,19)/  0.91490021622449999464D-01/                       
      DATA W( 5,19)/  0.11156664554733399472D+00/                       
      DATA W( 6,19)/  0.12875396253933622768D+00/                       
      DATA W( 7,19)/  0.14260670217360661178D+00/                       
      DATA W( 8,19)/  0.15276604206585966678D+00/                       
      DATA W( 9,19)/  0.15896884339395434765D+00/                       
      DATA W(10,19)/  0.16105444984878369598D+00/                       
      DATA X( 1,20)/  0.99312859918509492479D+00/                       
      DATA X( 2,20)/  0.96397192727791379127D+00/                       
      DATA X( 3,20)/  0.91223442825132590587D+00/                       
      DATA X( 4,20)/  0.83911697182221882339D+00/                       
      DATA X( 5,20)/  0.74633190646015079261D+00/                       
      DATA X( 6,20)/  0.63605368072651502545D+00/                       
      DATA X( 7,20)/  0.51086700195082709800D+00/                       
      DATA X( 8,20)/  0.37370608871541956067D+00/                       
      DATA X( 9,20)/  0.22778585114164507808D+00/                       
      DATA X(10,20)/  0.76526521133497333755D-01/                       
      DATA W( 1,20)/  0.17614007139152118312D-01/                       
      DATA W( 2,20)/  0.40601429800386941331D-01/                       
      DATA W( 3,20)/  0.62672048334109063570D-01/                       
      DATA W( 4,20)/  0.83276741576704748725D-01/                       
      DATA W( 5,20)/  0.10193011981724043504D+00/                       
      DATA W( 6,20)/  0.11819453196151841731D+00/                       
      DATA W( 7,20)/  0.13168863844917662690D+00/                       
      DATA W( 8,20)/  0.14209610931838205133D+00/                       
      DATA W( 9,20)/  0.14917298647260374679D+00/                       
      DATA W(10,20)/  0.15275338713072585070D+00/                       
      DATA X( 1,21)/  0.99375217062038950026D+00/                       
      DATA X( 2,21)/  0.96722683856630629432D+00/                       
      DATA X( 3,21)/  0.92009933415040082879D+00/                       
      DATA X( 4,21)/  0.85336336458331728365D+00/                       
      DATA X( 5,21)/  0.76843996347567790862D+00/                       
      DATA X( 6,21)/  0.66713880419741231931D+00/                       
      DATA X( 7,21)/  0.55161883588721980706D+00/                       
      DATA X( 8,21)/  0.42434212020743878357D+00/                       
      DATA X( 9,21)/  0.28802131680240109660D+00/                       
      DATA X(10,21)/  0.14556185416089509094D+00/                       
      DATA X(11,21)/  0.0                       /                       
      DATA W( 1,21)/  0.16017228257774333324D-01/                       
      DATA W( 2,21)/  0.36953789770852493800D-01/                       
      DATA W( 3,21)/  0.57134425426857208284D-01/                       
      DATA W( 4,21)/  0.76100113628379302017D-01/                       
      DATA W( 5,21)/  0.93444423456033861553D-01/                       
      DATA W( 6,21)/  0.10879729916714837766D+00/                       
      DATA W( 7,21)/  0.12183141605372853420D+00/                       
      DATA W( 8,21)/  0.13226893863333746178D+00/                       
      DATA W( 9,21)/  0.13988739479107315472D+00/                       
      DATA W(10,21)/  0.14452440398997005906D+00/                       
      DATA W(11,21)/  0.14608113364969042719D+00/                       
      DATA X( 1,22)/  0.99429458548239929207D+00/                       
      DATA X( 2,22)/  0.97006049783542872712D+00/                       
      DATA X( 3,22)/  0.92695677218717400052D+00/                       
      DATA X( 4,22)/  0.86581257772030013654D+00/                       
      DATA X( 5,22)/  0.78781680597920816200D+00/                       
      DATA X( 6,22)/  0.69448726318668278005D+00/                       
      DATA X( 7,22)/  0.58764040350691159296D+00/                       
      DATA X( 8,22)/  0.46935583798675702641D+00/                       
      DATA X( 9,22)/  0.34193582089208422516D+00/                       
      DATA X(10,22)/  0.20786042668822128548D+00/                       
      DATA X(11,22)/  0.69739273319722221214D-01/                       
      DATA W( 1,22)/  0.14627995298272200685D-01/                       
      DATA W( 2,22)/  0.33774901584814154793D-01/                       
      DATA W( 3,22)/  0.52293335152683285940D-01/                       
      DATA W( 4,22)/  0.69796468424520488095D-01/                       
      DATA W( 5,22)/  0.85941606217067727414D-01/                       
      DATA W( 6,22)/  0.10041414444288096493D+00/                       
      DATA W( 7,22)/  0.11293229608053921839D+00/                       
      DATA W( 8,22)/  0.12325237681051242429D+00/                       
      DATA W( 9,22)/  0.13117350478706237073D+00/                       
      DATA W(10,22)/  0.13654149834601517135D+00/                       
      DATA W(11,22)/  0.13925187285563199338D+00/                       
      DATA X( 1,23)/  0.99476933499755212352D+00/                       
      DATA X( 2,23)/  0.97254247121811523196D+00/                       
      DATA X( 3,23)/  0.93297108682601610235D+00/                       
      DATA X( 4,23)/  0.87675235827044166738D+00/                       
      DATA X( 5,23)/  0.80488840161883989215D+00/                       
      DATA X( 6,23)/  0.71866136313195019446D+00/                       
      DATA X( 7,23)/  0.61960987576364615639D+00/                       
      DATA X( 8,23)/  0.50950147784600754969D+00/                       
      DATA X( 9,23)/  0.39030103803029083142D+00/                       
      DATA X(10,23)/  0.26413568097034493053D+00/                       
      DATA X(11,23)/  0.13325682429846611093D+00/                       
      DATA X(12,23)/  0.0                       /                       
      DATA W( 1,23)/  0.13411859487141772081D-01/                       
      DATA W( 2,23)/  0.30988005856979444311D-01/                       
      DATA W( 3,23)/  0.48037671731084668572D-01/                       
      DATA W( 4,23)/  0.64232421408525852127D-01/                       
      DATA W( 5,23)/  0.79281411776718954923D-01/                       
      DATA W( 6,23)/  0.92915766060035147477D-01/                       
      DATA W( 7,23)/  0.10489209146454141007D+00/                       
      DATA W( 8,23)/  0.11499664022241136494D+00/                       
      DATA W( 9,23)/  0.12304908430672953047D+00/                       
      DATA W(10,23)/  0.12890572218808214998D+00/                       
      DATA W(11,23)/  0.13246203940469661737D+00/                       
      DATA W(12,23)/  0.13365457218610617535D+00/                       
      DATA X( 1,24)/  0.99518721999702136018D+00/                       
      DATA X( 2,24)/  0.97472855597130949820D+00/                       
      DATA X( 3,24)/  0.93827455200273275852D+00/                       
      DATA X( 4,24)/  0.88641552700440103421D+00/                       
      DATA X( 5,24)/  0.82000198597390292195D+00/                       
      DATA X( 6,24)/  0.74012419157855436424D+00/                       
      DATA X( 7,24)/  0.64809365193697556925D+00/                       
      DATA X( 8,24)/  0.54542147138883953566D+00/                       
      DATA X( 9,24)/  0.43379350762604513849D+00/                       
      DATA X(10,24)/  0.31504267969616337439D+00/                       
      DATA X(11,24)/  0.19111886747361630916D+00/                       
      DATA X(12,24)/  0.64056892862605626085D-01/                       
      DATA W( 1,24)/  0.12341229799987199547D-01/                       
      DATA W( 2,24)/  0.28531388628933663181D-01/                       
      DATA W( 3,24)/  0.44277438817419806169D-01/                       
      DATA W( 4,24)/  0.59298584915436780746D-01/                       
      DATA W( 5,24)/  0.73346481411080305734D-01/                       
      DATA W( 6,24)/  0.86190161531953275917D-01/                       
      DATA W( 7,24)/  0.97618652104113888270D-01/                       
      DATA W( 8,24)/  0.10744427011596563478D+00/                       
      DATA W( 9,24)/  0.11550566805372560135D+00/                       
      DATA W(10,24)/  0.12167047292780339120D+00/                       
      DATA W(11,24)/  0.12583745634682829612D+00/                       
      DATA W(12,24)/  0.12793819534675215697D+00/                       
      DATA X( 1,25)/  0.99555696979049809791D+00/                       
      DATA X( 2,25)/  0.97666392145951751150D+00/                       
      DATA X( 3,25)/  0.94297457122897433941D+00/                       
      DATA X( 4,25)/  0.89499199787827536885D+00/                       
      DATA X( 5,25)/  0.83344262876083400142D+00/                       
      DATA X( 6,25)/  0.75925926303735763058D+00/                       
      DATA X( 7,25)/  0.67356636847346836449D+00/                       
      DATA X( 8,25)/  0.57766293024122296772D+00/                       
      DATA X( 9,25)/  0.47300273144571496052D+00/                       
      DATA X(10,25)/  0.36117230580938783774D+00/                       
      DATA X(11,25)/  0.24386688372098843205D+00/                       
      DATA X(12,25)/  0.12286469261071039639D+00/                       
      DATA X(13,25)/  0.0                       /                       
      DATA W( 1,25)/  0.11393798501026287948D-01/                       
      DATA W( 2,25)/  0.26354986615032137262D-01/                       
      DATA W( 3,25)/  0.40939156701306312656D-01/                       
      DATA W( 4,25)/  0.54904695975835191926D-01/                       
      DATA W( 5,25)/  0.68038333812356917207D-01/                       
      DATA W( 6,25)/  0.80140700335001018013D-01/                       
      DATA W( 7,25)/  0.91028261982963649811D-01/                       
      DATA W( 8,25)/  0.10053594906705064420D+00/                       
      DATA W( 9,25)/  0.10851962447426365312D+00/                       
      DATA W(10,25)/  0.11485825914571164834D+00/                       
      DATA W(11,25)/  0.11945576353578477223D+00/                       
      DATA W(12,25)/  0.12224244299031004169D+00/                       
      DATA W(13,25)/  0.12317605372671545120D+00/                       
      DATA X( 1,26)/  0.99588570114561692900D+00/                       
      DATA X( 2,26)/  0.97838544595647099110D+00/                       
      DATA X( 3,26)/  0.94715906666171425014D+00/                       
      DATA X( 4,26)/  0.90263786198430707422D+00/                       
      DATA X( 5,26)/  0.84544594278849801880D+00/                       
      DATA X( 6,26)/  0.77638594882067885619D+00/                       
      DATA X( 7,26)/  0.69642726041995726486D+00/                       
      DATA X( 8,26)/  0.60669229301761806323D+00/                       
      DATA X( 9,26)/  0.50844071482450571770D+00/                       
      DATA X(10,26)/  0.40305175512348630648D+00/                       
      DATA X(11,26)/  0.29200483948595689514D+00/                       
      DATA X(12,26)/  0.17685882035689018397D+00/                       
      DATA X(13,26)/  0.59230093429313207094D-01/                       
      DATA W( 1,26)/  0.10551372617343007156D-01/                       
      DATA W( 2,26)/  0.24417851092631908790D-01/                       
      DATA W( 3,26)/  0.37962383294362763950D-01/                       
      DATA W( 4,26)/  0.50975825297147811998D-01/                       
      DATA W( 5,26)/  0.63274046329574835539D-01/                       
      DATA W( 6,26)/  0.74684149765659745887D-01/                       
      DATA W( 7,26)/  0.85045894313485239210D-01/                       
      DATA W( 8,26)/  0.94213800355914148464D-01/                       
      DATA W( 9,26)/  0.10205916109442542324D+00/                       
      DATA W(10,26)/  0.10847184052857659066D+00/                       
      DATA W(11,26)/  0.11336181654631966655D+00/                       
      DATA W(12,26)/  0.11666044348529658204D+00/                       
      DATA W(13,26)/  0.11832141527926227652D+00/                       
      DATA X( 1,27)/  0.99617926288898856694D+00/                       
      DATA X( 2,27)/  0.97992347596150122286D+00/                       
      DATA X( 3,27)/  0.95090055781470500685D+00/                       
      DATA X( 4,27)/  0.90948232067749110430D+00/                       
      DATA X( 5,27)/  0.85620790801829449030D+00/                       
      DATA X( 6,27)/  0.79177163907050822714D+00/                       
      DATA X( 7,27)/  0.71701347373942369929D+00/                       
      DATA X( 8,27)/  0.63290797194649514093D+00/                       
      DATA X( 9,27)/  0.54055156457945689490D+00/                       
      DATA X(10,27)/  0.44114825175002688059D+00/                       
      DATA X(11,27)/  0.33599390363850889973D+00/                       
      DATA X(12,27)/  0.22645936543953685886D+00/                       
      DATA X(13,27)/  0.11397258560952996693D+00/                       
      DATA X(14,27)/  0.0                       /                       
      DATA W( 1,27)/  0.97989960512943602612D-02/                       
      DATA W( 2,27)/  0.22686231596180623196D-01/                       
      DATA W( 3,27)/  0.35297053757419711023D-01/                       
      DATA W( 4,27)/  0.47449412520615062704D-01/                       
      DATA W( 5,27)/  0.58983536859833599110D-01/                       
      DATA W( 6,27)/  0.69748823766245592984D-01/                       
      DATA W( 7,27)/  0.79604867773057771263D-01/                       
      DATA W( 8,27)/  0.88423158543756950194D-01/                       
      DATA W( 9,27)/  0.96088727370028507566D-01/                       
      DATA W(10,27)/  0.10250163781774579867D+00/                       
      DATA W(11,27)/  0.10757828578853318721D+00/                       
      DATA W(12,27)/  0.11125248835684519267D+00/                       
      DATA W(13,27)/  0.11347634610896514862D+00/                       
      DATA W(14,27)/  0.11422086737895698905D+00/                       
      DATA X( 1,28)/  0.99644249757395444995D+00/                       
      DATA X( 2,28)/  0.98130316537087275369D+00/                       
      DATA X( 3,28)/  0.95425928062893819725D+00/                       
      DATA X( 4,28)/  0.91563302639213207387D+00/                       
      DATA X( 5,28)/  0.86589252257439504894D+00/                       
      DATA X( 6,28)/  0.80564137091717917145D+00/                       
      DATA X( 7,28)/  0.73561087801363177203D+00/                       
      DATA X( 8,28)/  0.65665109403886496122D+00/                       
      DATA X( 9,28)/  0.56972047181140171931D+00/                       
      DATA X(10,28)/  0.47587422495511826103D+00/                       
      DATA X(11,28)/  0.37625151608907871022D+00/                       
      DATA X(12,28)/  0.27206162763517807768D+00/                       
      DATA X(13,28)/  0.16456928213338077128D+00/                       
      DATA X(14,28)/  0.55079289884034270427D-01/                       
      DATA W( 1,28)/  0.91242825930945177388D-02/                       
      DATA W( 2,28)/  0.21132112592771259752D-01/                       
      DATA W( 3,28)/  0.32901427782304379978D-01/                       
      DATA W( 4,28)/  0.44272934759004227840D-01/                       
      DATA W( 5,28)/  0.55107345675716745431D-01/                       
      DATA W( 6,28)/  0.65272923966999595793D-01/                       
      DATA W( 7,28)/  0.74646214234568779024D-01/                       
      DATA W( 8,28)/  0.83113417228901218390D-01/                       
      DATA W( 9,28)/  0.90571744393032840942D-01/                       
      DATA W(10,28)/  0.96930657997929915850D-01/                       
      DATA W(11,28)/  0.10211296757806076981D+00/                       
      DATA W(12,28)/  0.10605576592284641791D+00/                       
      DATA W(13,28)/  0.10871119225829413525D+00/                       
      DATA W(14,28)/  0.11004701301647519628D+00/                       
      DATA X( 1,29)/  0.99667944226059658616D+00/                       
      DATA X( 2,29)/  0.98254550526141317487D+00/                       
      DATA X( 3,29)/  0.95728559577808772580D+00/                       
      DATA X( 4,29)/  0.92118023295305878509D+00/                       
      DATA X( 5,29)/  0.87463780492010279042D+00/                       
      DATA X( 6,29)/  0.81818548761525244499D+00/                       
      DATA X( 7,29)/  0.75246285173447713391D+00/                       
      DATA X( 8,29)/  0.67821453760268651516D+00/                       
      DATA X( 9,29)/  0.59628179713822782038D+00/                       
      DATA X(10,29)/  0.50759295512422764210D+00/                       
      DATA X(11,29)/  0.41315288817400866389D+00/                       
      DATA X(12,29)/  0.31403163786763993495D+00/                       
      DATA X(13,29)/  0.21135228616600107451D+00/                       
      DATA X(14,29)/  0.10627823013267923017D+00/                       
      DATA X(15,29)/  0.0                       /                       
      DATA W( 1,29)/  0.85169038787464096543D-02/                       
      DATA W( 2,29)/  0.19732085056122705984D-01/                       
      DATA W( 3,29)/  0.30740492202093622644D-01/                       
      DATA W( 4,29)/  0.41402062518682836105D-01/                       
      DATA W( 5,29)/  0.51594826902497923913D-01/                       
      DATA W( 6,29)/  0.61203090657079138542D-01/                       
      DATA W( 7,29)/  0.70117933255051278570D-01/                       
      DATA W( 8,29)/  0.78238327135763783828D-01/                       
      DATA W( 9,29)/  0.85472257366172527545D-01/                       
      DATA W(10,29)/  0.91737757139258763348D-01/                       
      DATA W(11,29)/  0.96963834094408606302D-01/                       
      DATA W(12,29)/  0.10109127375991496612D+00/                       
      DATA W(13,29)/  0.10407331007772937391D+00/                       
      DATA W(14,29)/  0.10587615509732094141D+00/                       
      DATA W(15,29)/  0.10647938171831424425D+00/                       
      DATA X( 1,30)/  0.99689348407464954027D+00/                       
      DATA X( 2,30)/  0.98366812327974720997D+00/                       
      DATA X( 3,30)/  0.96002186496830751222D+00/                       
      DATA X( 4,30)/  0.92620004742927432588D+00/                       
      DATA X( 5,30)/  0.88256053579205268154D+00/                       
      DATA X( 6,30)/  0.82956576238276839744D+00/                       
      DATA X( 7,30)/  0.76777743210482619492D+00/                       
      DATA X( 8,30)/  0.69785049479331579693D+00/                       
      DATA X( 9,30)/  0.62052618298924286114D+00/                       
      DATA X(10,30)/  0.53662414814201989926D+00/                       
      DATA X(11,30)/  0.44703376953808917678D+00/                       
      DATA X(12,30)/  0.35270472553087811347D+00/                       
      DATA X(13,30)/  0.25463692616788984644D+00/                       
      DATA X(14,30)/  0.15386991360858354696D+00/                       
      DATA X(15,30)/  0.51471842555317695833D-01/                       
      DATA W( 1,30)/  0.79681924961666056155D-02/                       
      DATA W( 2,30)/  0.18466468311090959142D-01/                       
      DATA W( 3,30)/  0.28784707883323369350D-01/                       
      DATA W( 4,30)/  0.38799192569627049597D-01/                       
      DATA W( 5,30)/  0.48402672830594052903D-01/                       
      DATA W( 6,30)/  0.57493156217619066482D-01/                       
      DATA W( 7,30)/  0.65974229882180495128D-01/                       
      DATA W( 8,30)/  0.73755974737705206268D-01/                       
      DATA W( 9,30)/  0.80755895229420215355D-01/                       
      DATA W(10,30)/  0.86899787201082979802D-01/                       
      DATA W(11,30)/  0.92122522237786128718D-01/                       
      DATA W(12,30)/  0.96368737174644259639D-01/                       
      DATA W(13,30)/  0.99593420586795267063D-01/                       
      DATA W(14,30)/  0.10176238974840550460D+00/                       
      DATA W(15,30)/  0.10285265289355884034D+00/                       
      DATA X( 1,31)/  0.99708748181947707406D+00/                       
      DATA X( 2,31)/  0.98468590966515248400D+00/                       
      DATA X( 3,31)/  0.96250392509294966179D+00/                       
      DATA X( 4,31)/  0.93075699789664816496D+00/                       
      DATA X( 5,31)/  0.88976002994827104337D+00/                       
      DATA X( 6,31)/  0.83992032014626734009D+00/                       
      DATA X( 7,31)/  0.78173314841662494041D+00/                       
      DATA X( 8,31)/  0.71577678458685328391D+00/                       
      DATA X( 9,31)/  0.64270672292426034618D+00/                       
      DATA X(10,31)/  0.56324916140714926272D+00/                       
      DATA X(11,31)/  0.47819378204490248044D+00/                       
      DATA X(12,31)/  0.38838590160823294306D+00/                       
      DATA X(13,31)/  0.29471806998170161662D+00/                       
      DATA X(14,31)/  0.19812119933557062877D+00/                       
      DATA X(15,31)/  0.99555312152341520325D-01/                       
      DATA X(16,31)/  0.0                       /                       
      DATA W( 1,31)/  0.74708315792487758587D-02/                       
      DATA W( 2,31)/  0.17318620790310582463D-01/                       
      DATA W( 3,31)/  0.27009019184979421801D-01/                       
      DATA W( 4,31)/  0.36432273912385464024D-01/                       
      DATA W( 5,31)/  0.45493707527201102902D-01/                       
      DATA W( 6,31)/  0.54103082424916853712D-01/                       
      DATA W( 7,31)/  0.62174786561028426910D-01/                       
      DATA W( 8,31)/  0.69628583235410366168D-01/                       
      DATA W( 9,31)/  0.76390386598776616426D-01/                       
      DATA W(10,31)/  0.82392991761589263904D-01/                       
      DATA W(11,31)/  0.87576740608477876126D-01/                       
      DATA W(12,31)/  0.91890113893641478215D-01/                       
      DATA W(13,31)/  0.95290242912319512807D-01/                       
      DATA W(14,31)/  0.97743335386328725093D-01/                       
      DATA W(15,31)/  0.99225011226672307875D-01/                       
      DATA W(16,31)/  0.99720544793426451428D-01/                       
      DATA X( 1,32)/  0.99726386184948156354D+00/                       
      DATA X( 2,32)/  0.98561151154526833540D+00/                       
      DATA X( 3,32)/  0.96476225558750643077D+00/                       
      DATA X( 4,32)/  0.93490607593773968917D+00/                       
      DATA X( 5,32)/  0.89632115576605212397D+00/                       
      DATA X( 6,32)/  0.84936761373256997013D+00/                       
      DATA X( 7,32)/  0.79448379596794240696D+00/                       
      DATA X( 8,32)/  0.73218211874028968039D+00/                       
      DATA X( 9,32)/  0.66304426693021520098D+00/                       
      DATA X(10,32)/  0.58771575724076232904D+00/                       
      DATA X(11,32)/  0.50689990893222939002D+00/                       
      DATA X(12,32)/  0.42135127613063534536D+00/                       
      DATA X(13,32)/  0.33186860228212764978D+00/                       
      DATA X(14,32)/  0.23928736225213707454D+00/                       
      DATA X(15,32)/  0.14447196158279649349D+00/                       
      DATA X(16,32)/  0.48307665687738316235D-01/                       
      DATA W( 1,32)/  0.70186100094700966004D-02/                       
      DATA W( 2,32)/  0.16274394730905670605D-01/                       
      DATA W( 3,32)/  0.25392065309262059456D-01/                       
      DATA W( 4,32)/  0.34273862913021433103D-01/                       
      DATA W( 5,32)/  0.42835898022226680657D-01/                       
      DATA W( 6,32)/  0.50998059262376176196D-01/                       
      DATA W( 7,32)/  0.58684093478535547145D-01/                       
      DATA W( 8,32)/  0.65822222776361846838D-01/                       
      DATA W( 9,32)/  0.72345794108848506225D-01/                       
      DATA W(10,32)/  0.78193895787070306472D-01/                       
      DATA W(11,32)/  0.83311924226946755222D-01/                       
      DATA W(12,32)/  0.87652093004403811143D-01/                       
      DATA W(13,32)/  0.91173878695763884713D-01/                       
      DATA W(14,32)/  0.93844399080804565639D-01/                       
      DATA W(15,32)/  0.95638720079274859419D-01/                       
      DATA W(16,32)/  0.96540088514727800567D-01/                       
      DATA X( 1,33)/  0.99742469424645521727D+00/                       
      DATA X( 2,33)/  0.98645572623064248811D+00/                       
      DATA X( 3,33)/  0.96682290968999276893D+00/                       
      DATA X( 4,33)/  0.93869437261116835036D+00/                       
      DATA X( 5,33)/  0.90231676774343358304D+00/                       
      DATA X( 6,33)/  0.85800965267650406464D+00/                       
      DATA X( 7,33)/  0.80616235627416658980D+00/                       
      DATA X( 8,33)/  0.74723049644956215786D+00/                       
      DATA X( 9,33)/  0.68173195996974278627D+00/                       
      DATA X(10,33)/  0.61024234583637902731D+00/                       
      DATA X(11,33)/  0.53338990478634764355D+00/                       
      DATA X(12,33)/  0.45185001727245069573D+00/                       
      DATA X(13,33)/  0.36633925774807334107D+00/                       
      DATA X(14,33)/  0.27760909715249702940D+00/                       
      DATA X(15,33)/  0.18643929882799157234D+00/                       
      DATA X(16,33)/  0.93631065854733385671D-01/                       
      DATA X(17,33)/  0.0                       /                       
      DATA W( 1,33)/  0.66062278475873780586D-02/                       
      DATA W( 2,33)/  0.15321701512934676128D-01/                       
      DATA W( 3,33)/  0.23915548101749480351D-01/                       
      DATA W( 4,33)/  0.32300358632328953282D-01/                       
      DATA W( 5,33)/  0.40401541331669591563D-01/                       
      DATA W( 6,33)/  0.48147742818711695670D-01/                       
      DATA W( 7,33)/  0.55470846631663561285D-01/                       
      DATA W( 8,33)/  0.62306482530317480032D-01/                       
      DATA W( 9,33)/  0.68594572818656712806D-01/                       
      DATA W(10,33)/  0.74279854843954149342D-01/                       
      DATA W(11,33)/  0.79312364794886738364D-01/                       
      DATA W(12,33)/  0.83647876067038707614D-01/                       
      DATA W(13,33)/  0.87248287618844337607D-01/                       
      DATA W(14,33)/  0.90081958660638577240D-01/                       
      DATA W(15,33)/  0.92123986643316846213D-01/                       
      DATA W(16,33)/  0.93356426065596116161D-01/                       
      DATA W(17,33)/  0.93768446160209996567D-01/                       
      DATA X( 1,34)/  0.99757175379084191924D+00/                       
      DATA X( 2,34)/  0.98722781640630948505D+00/                       
      DATA X( 3,34)/  0.96870826253334428176D+00/                       
      DATA X( 4,34)/  0.94216239740510709163D+00/                       
      DATA X( 5,34)/  0.90780967771832446880D+00/                       
      DATA X( 6,34)/  0.86593463833456446926D+00/                       
      DATA X( 7,34)/  0.81688422790093366459D+00/                       
      DATA X( 8,34)/  0.76106487662987301419D+00/                       
      DATA X( 9,34)/  0.69893911321626290793D+00/                       
      DATA X(10,34)/  0.63102172708052854532D+00/                       
      DATA X(11,34)/  0.55787550066974664274D+00/                       
      DATA X(12,34)/  0.48010654519032703419D+00/                       
      DATA X(13,34)/  0.39835927775864594063D+00/                       
      DATA X(14,34)/  0.31331108133946324746D+00/                       
      DATA X(15,34)/  0.22566669161644948387D+00/                       
      DATA X(16,34)/  0.13615235725918297589D+00/                       
      DATA X(17,34)/  0.45509821953102542749D-01/                       
      DATA W( 1,34)/  0.62291405559086847186D-02/                       
      DATA W( 2,34)/  0.14450162748595035415D-01/                       
      DATA W( 3,34)/  0.22563721985494970084D-01/                       
      DATA W( 4,34)/  0.30491380638446131809D-01/                       
      DATA W( 5,34)/  0.38166593796387516322D-01/                       
      DATA W( 6,34)/  0.45525611523353272454D-01/                       
      DATA W( 7,34)/  0.52507414572678106168D-01/                       
      DATA W( 8,34)/  0.59054135827524493194D-01/                       
      DATA W( 9,34)/  0.65111521554076411379D-01/                       
      DATA W(10,34)/  0.70629375814255724999D-01/                       
      DATA W(11,34)/  0.75561974660031931271D-01/                       
      DATA W(12,34)/  0.79868444339771844739D-01/                       
      DATA W(13,34)/  0.83513099699845655187D-01/                       
      DATA W(14,34)/  0.86465739747035749784D-01/                       
      DATA W(15,34)/  0.88701897835693869287D-01/                       
      DATA W(16,34)/  0.90203044370640729574D-01/                       
      DATA W(17,34)/  0.90956740330259873615D-01/                       
      DATA X( 1,35)/  0.99770656909960029726D+00/                       
      DATA X( 2,35)/  0.98793576444385149804D+00/                       
      DATA X( 3,35)/  0.97043761603922983322D+00/                       
      DATA X( 4,35)/  0.94534514820782732954D+00/                       
      DATA X( 5,35)/  0.91285426135931761446D+00/                       
      DATA X( 6,35)/  0.87321912502522233152D+00/                       
      DATA X( 7,35)/  0.82674989909222540683D+00/                       
      DATA X( 8,35)/  0.77381025228691255527D+00/                       
      DATA X( 9,35)/  0.71481450155662878326D+00/                       
      DATA X(10,35)/  0.65022436466589038868D+00/                       
      DATA X(11,35)/  0.58054534474976450993D+00/                       
      DATA X(12,35)/  0.50632277324148861502D+00/                       
      DATA X(13,35)/  0.42813754151781425419D+00/                       
      DATA X(14,35)/  0.34660155443081394588D+00/                       
      DATA X(15,35)/  0.26235294120929605797D+00/                       
      DATA X(16,35)/  0.17605106116598956997D+00/                       
      DATA X(17,35)/  0.88371343275659263601D-01/                       
      DATA X(18,35)/  0.0                       /                       
      DATA W( 1,35)/  0.58834334204430849758D-02/                       
      DATA W( 2,35)/  0.13650828348361492266D-01/                       
      DATA W( 3,35)/  0.21322979911483580883D-01/                       
      DATA W( 4,35)/  0.28829260108894254049D-01/                       
      DATA W( 5,35)/  0.36110115863463380533D-01/                       
      DATA W( 6,35)/  0.43108422326170218782D-01/                       
      DATA W( 7,35)/  0.49769370401353529805D-01/                       
      DATA W( 8,35)/  0.56040816212370128578D-01/                       
      DATA W( 9,35)/  0.61873671966080188887D-01/                       
      DATA W(10,35)/  0.67222285269086903964D-01/                       
      DATA W(11,35)/  0.72044794772560064665D-01/                       
      DATA W(12,35)/  0.76303457155442053539D-01/                       
      DATA W(13,35)/  0.79964942242324262933D-01/                       
      DATA W(14,35)/  0.83000593728856588380D-01/                       
      DATA W(15,35)/  0.85386653392099125226D-01/                       
      DATA W(16,35)/  0.87104446997183534243D-01/                       
      DATA W(17,35)/  0.88140530430275462971D-01/                       
      DATA W(18,35)/  0.88486794907104290638D-01/                       
      DATA X( 1,36)/  0.99783046248408583620D+00/                       
      DATA X( 2,36)/  0.98858647890221223807D+00/                       
      DATA X( 3,36)/  0.97202769104969794934D+00/                       
      DATA X( 4,36)/  0.94827298439950754520D+00/                       
      DATA X( 5,36)/  0.91749777451565906608D+00/                       
      DATA X( 6,36)/  0.87992980089039713198D+00/                       
      DATA X( 7,36)/  0.83584716699247530642D+00/                       
      DATA X( 8,36)/  0.78557623013220651283D+00/                       
      DATA X( 9,36)/  0.72948917159355658209D+00/                       
      DATA X(10,36)/  0.66800123658552106210D+00/                       
      DATA X(11,36)/  0.60156765813598053508D+00/                       
      DATA X(12,36)/  0.53068028592624516164D+00/                       
      DATA X(13,36)/  0.45586394443342026721D+00/                       
      DATA X(14,36)/  0.37767254711968921632D+00/                       
      DATA X(15,36)/  0.29668499534402827050D+00/                       
      DATA X(16,36)/  0.21350089231686557894D+00/                       
      DATA X(17,36)/  0.12873610380938478865D+00/                       
      DATA X(18,36)/  0.43018198473708607227D-01/                       
      DATA W( 1,36)/  0.55657196642450453613D-02/                       
      DATA W( 2,36)/  0.12915947284065574405D-01/                       
      DATA W( 3,36)/  0.20181515297735471532D-01/                       
      DATA W( 4,36)/  0.27298621498568779094D-01/                       
      DATA W( 5,36)/  0.34213810770307229921D-01/                       
      DATA W( 6,36)/  0.40875750923644895474D-01/                       
      DATA W( 7,36)/  0.47235083490265978417D-01/                       
      DATA W( 8,36)/  0.53244713977759919092D-01/                       
      DATA W( 9,36)/  0.58860144245324817310D-01/                       
      DATA W(10,36)/  0.64039797355015489556D-01/                       
      DATA W(11,36)/  0.68745323835736442614D-01/                       
      DATA W(12,36)/  0.72941885005653061354D-01/                       
      DATA W(13,36)/  0.76598410645870674529D-01/                       
      DATA W(14,36)/  0.79687828912071601909D-01/                       
      DATA W(15,36)/  0.82187266704339709517D-01/                       
      DATA W(16,36)/  0.84078218979661934933D-01/                       
      DATA W(17,36)/  0.85346685739338627492D-01/                       
      DATA W(18,36)/  0.85983275670394747490D-01/                       
      DATA X( 1,37)/  0.99794458247791364894D+00/                       
      DATA X( 2,37)/  0.98918596321431918668D+00/                       
      DATA X( 3,37)/  0.97349303005648574433D+00/                       
      DATA X( 4,37)/  0.95097234326209482133D+00/                       
      DATA X( 5,37)/  0.92178143741246374267D+00/                       
      DATA X( 6,37)/  0.88612496215548607895D+00/                       
      DATA X( 7,37)/  0.84425298734055596799D+00/                       
      DATA X( 8,37)/  0.79645920050990229339D+00/                       
      DATA X( 9,37)/  0.74307883398196526255D+00/                       
      DATA X(10,37)/  0.68448630913095935745D+00/                       
      DATA X(11,37)/  0.62109260840892448315D+00/                       
      DATA X(12,37)/  0.55334239186158178124D+00/                       
      DATA X(13,37)/  0.48171087780320555415D+00/                       
      DATA X(14,37)/  0.40670050931832611010D+00/                       
      DATA X(15,37)/  0.32883742988370699950D+00/                       
      DATA X(16,37)/  0.24866779279136575881D+00/                       
      DATA X(17,37)/  0.16675393023985197697D+00/                       
      DATA X(18,37)/  0.83670408954769901943D-01/                       
      DATA X(19,37)/  0.0                       /                       
      DATA W( 1,37)/  0.52730572794979393517D-02/                       
      DATA W( 2,37)/  0.12238780100307556526D-01/                       
      DATA W( 3,37)/  0.19129044489083966044D-01/                       
      DATA W( 4,37)/  0.25886036990558933523D-01/                       
      DATA W( 5,37)/  0.32461639847521481067D-01/                       
      DATA W( 6,37)/  0.38809602501934544489D-01/                       
      DATA W( 7,37)/  0.44885364662437166657D-01/                       
      DATA W( 8,37)/  0.50646297654824601604D-01/                       
      DATA W( 9,37)/  0.56051987998274917809D-01/                       
      DATA W(10,37)/  0.61064516523225986131D-01/                       
      DATA W(11,37)/  0.65648722872751249484D-01/                       
      DATA W(12,37)/  0.69772451555700344885D-01/                       
      DATA W(13,37)/  0.73406777248488172725D-01/                       
      DATA W(14,37)/  0.76526207570529237886D-01/                       
      DATA W(15,37)/  0.79108861837529380767D-01/                       
      DATA W(16,37)/  0.81136624508465030510D-01/                       
      DATA W(17,37)/  0.82595272236437250891D-01/                       
      DATA W(18,37)/  0.83474573625862787252D-01/                       
      DATA W(19,37)/  0.83768360993138904797D-01/                       
      DATA X( 1,38)/  0.99804993053568761981D+00/                       
      DATA X( 2,38)/  0.98973945426638557194D+00/                       
      DATA X( 3,38)/  0.97484632859015350764D+00/                       
      DATA X( 4,38)/  0.95346633093352959567D+00/                       
      DATA X( 5,38)/  0.92574133204858439683D+00/                       
      DATA X( 6,38)/  0.89185573900463221679D+00/                       
      DATA X( 7,38)/  0.85203502193236218886D+00/                       
      DATA X( 8,38)/  0.80654416760531681555D+00/                       
      DATA X( 9,38)/  0.75568590375397068074D+00/                       
      DATA X(10,38)/  0.69979868037918435591D+00/                       
      DATA X(11,38)/  0.63925441582968170718D+00/                       
      DATA X(12,38)/  0.57445602104780708113D+00/                       
      DATA X(13,38)/  0.50583471792793110324D+00/                       
      DATA X(14,38)/  0.43384716943237648437D+00/                       
      DATA X(15,38)/  0.35897244047943501326D+00/                       
      DATA X(16,38)/  0.28170880979016526136D+00/                       
      DATA X(17,38)/  0.20257045389211670320D+00/                       
      DATA X(18,38)/  0.12208402533786741987D+00/                       
      DATA X(19,38)/  0.40785147904578239913D-01/                       
      DATA W( 1,38)/  0.50028807496393456759D-02/                       
      DATA W( 2,38)/  0.11613444716468674178D-01/                       
      DATA W( 3,38)/  0.18156577709613236899D-01/                       
      DATA W( 4,38)/  0.24579739738232375895D-01/                       
      DATA W( 5,38)/  0.30839500545175054659D-01/                       
      DATA W( 6,38)/  0.36894081594024738165D-01/                       
      DATA W( 7,38)/  0.42703158504674434236D-01/                       
      DATA W( 8,38)/  0.48228061860758683374D-01/                       
      DATA W( 9,38)/  0.53432019910332319974D-01/                       
      DATA W(10,38)/  0.58280399146997206022D-01/                       
      DATA W(11,38)/  0.62740933392133054053D-01/                       
      DATA W(12,38)/  0.66783937979140411935D-01/                       
      DATA W(13,38)/  0.70382507066898954739D-01/                       
      DATA W(14,38)/  0.73512692584743457145D-01/                       
      DATA W(15,38)/  0.76153663548446396066D-01/                       
      DATA W(16,38)/  0.78287844658210948075D-01/                       
      DATA W(17,38)/  0.79901033243527821586D-01/                       
      DATA W(18,38)/  0.80982493770597100623D-01/                       
      DATA W(19,38)/  0.81525029280385786699D-01/                       
      DATA X( 1,39)/  0.99814738306643290601D+00/                       
      DATA X( 2,39)/  0.99025153685468598364D+00/                       
      DATA X( 3,39)/  0.97609870933347105384D+00/                       
      DATA X( 4,39)/  0.95577521232465227711D+00/                       
      DATA X( 5,39)/  0.92940914848673822970D+00/                       
      DATA X( 6,39)/  0.89716711929299288785D+00/                       
      DATA X( 7,39)/  0.85925293799990615391D+00/                       
      DATA X( 8,39)/  0.81590629743014310435D+00/                       
      DATA X( 9,39)/  0.76740124293106349983D+00/                       
      DATA X(10,39)/  0.71404443589453467913D+00/                       
      DATA X(11,39)/  0.65617321343201091073D+00/                       
      DATA X(12,39)/  0.59415345495727798869D+00/                       
      DATA X(13,39)/  0.52837726866043747390D+00/                       
      DATA X(14,39)/  0.45926051230913604866D+00/                       
      DATA X(15,39)/  0.38724016397156145585D+00/                       
      DATA X(16,39)/  0.31277155924818592254D+00/                       
      DATA X(17,39)/  0.23632551246183576734D+00/                       
      DATA X(18,39)/  0.15838533999783779992D+00/                       
      DATA X(19,39)/  0.79443804608755477582D-01/                       
      DATA X(20,39)/  0.0                       /                       
      DATA W( 1,39)/  0.47529446916351013708D-02/                       
      DATA W( 2,39)/  0.11034788939164594243D-01/                       
      DATA W( 3,39)/  0.17256229093724919041D-01/                       
      DATA W( 4,39)/  0.23369384832178164595D-01/                       
      DATA W( 5,39)/  0.29334955983903378592D-01/                       
      DATA W( 6,39)/  0.35115111498131330761D-01/                       
      DATA W( 7,39)/  0.40673276847933843939D-01/                       
      DATA W( 8,39)/  0.45974301108916631884D-01/                       
      DATA W( 9,39)/  0.50984665292129405214D-01/                       
      DATA W(10,39)/  0.55672690340916299907D-01/                       
      DATA W(11,39)/  0.60008736088596149575D-01/                       
      DATA W(12,39)/  0.63965388138682388987D-01/                       
      DATA W(13,39)/  0.67517630966231265363D-01/                       
      DATA W(14,39)/  0.70643005970608760770D-01/                       
      DATA W(15,39)/  0.73321753414268617381D-01/                       
      DATA W(16,39)/  0.75536937322836057705D-01/                       
      DATA W(17,39)/  0.77274552544682016729D-01/                       
      DATA W(18,39)/  0.78523613287371176725D-01/                       
      DATA W(19,39)/  0.79276222568368471010D-01/                       
      DATA W(20,39)/  0.79527622139442852417D-01/                       
      DATA X( 1,40)/  0.99823770971055920035D+00/                       
      DATA X( 2,40)/  0.99072623869945700645D+00/                       
      DATA X( 3,40)/  0.97725994998377426266D+00/                       
      DATA X( 4,40)/  0.95791681921379165580D+00/                       
      DATA X( 5,40)/  0.93281280827867653336D+00/                       
      DATA X( 6,40)/  0.90209880696887429673D+00/                       
      DATA X( 7,40)/  0.86595950321225950382D+00/                       
      DATA X( 8,40)/  0.82461223083331166320D+00/                       
      DATA X( 9,40)/  0.77830565142651938769D+00/                       
      DATA X(10,40)/  0.72731825518992710328D+00/                       
      DATA X(11,40)/  0.67195668461417954838D+00/                       
      DATA X(12,40)/  0.61255388966798023795D+00/                       
      DATA X(13,40)/  0.54946712509512820208D+00/                       
      DATA X(14,40)/  0.48307580168617871291D+00/                       
      DATA X(15,40)/  0.41377920437160500152D+00/                       
      DATA X(16,40)/  0.34199409082575847301D+00/                       
      DATA X(17,40)/  0.26815218500725368114D+00/                       
      DATA X(18,40)/  0.19269758070137109972D+00/                       
      DATA X(19,40)/  0.11608407067525520848D+00/                       
      DATA X(20,40)/  0.38772417506050821933D-01/                       
      DATA W( 1,40)/  0.45212770985331912585D-02/                       
      DATA W( 2,40)/  0.10498284531152813615D-01/                       
      DATA W( 3,40)/  0.16421058381907888713D-01/                       
      DATA W( 4,40)/  0.22245849194166957262D-01/                       
      DATA W( 5,40)/  0.27937006980023401098D-01/                       
      DATA W( 6,40)/  0.33460195282547847393D-01/                       
      DATA W( 7,40)/  0.38782167974472017640D-01/                       
      DATA W( 8,40)/  0.43870908185673271992D-01/                       
      DATA W( 9,40)/  0.48695807635072232061D-01/                       
      DATA W(10,40)/  0.53227846983936824355D-01/                       
      DATA W(11,40)/  0.57439769099391551367D-01/                       
      DATA W(12,40)/  0.61306242492928939167D-01/                       
      DATA W(13,40)/  0.64804013456601038075D-01/                       
      DATA W(14,40)/  0.67912045815233903826D-01/                       
      DATA W(15,40)/  0.70611647391286779695D-01/                       
      DATA W(16,40)/  0.72886582395804059061D-01/                       
      DATA W(17,40)/  0.74723169057968264200D-01/                       
      DATA W(18,40)/  0.76110361900626242372D-01/                       
      DATA W(19,40)/  0.77039818164247965588D-01/                       
      DATA W(20,40)/  0.77505947978424811264D-01/                       
      DATA X( 1,41)/  0.99832158857477144152D+00/                       
      DATA X( 2,41)/  0.99116710969901630825D+00/                       
      DATA X( 3,41)/  0.97833867356108338447D+00/                       
      DATA X( 4,41)/  0.95990689173034622610D+00/                       
      DATA X( 5,41)/  0.93597698749785382568D+00/                       
      DATA X( 6,41)/  0.90668594475810117296D+00/                       
      DATA X( 7,41)/  0.87220151169244140883D+00/                       
      DATA X( 8,41)/  0.83272120040136133124D+00/                       
      DATA X( 9,41)/  0.78847114504740937274D+00/                       
      DATA X(10,41)/  0.73970480306992618106D+00/                       
      DATA X(11,41)/  0.68670150203495128958D+00/                       
      DATA X(12,41)/  0.62976483907219632049D+00/                       
      DATA X(13,41)/  0.56922094161021586965D+00/                       
      DATA X(14,41)/  0.50541659919940603271D+00/                       
      DATA X(15,41)/  0.43871727705140708852D+00/                       
      DATA X(16,41)/  0.36950502264048144143D+00/                       
      DATA X(17,41)/  0.29817627734182486592D+00/                       
      DATA X(18,41)/  0.22513960563342277561D+00/                       
      DATA X(19,41)/  0.15081335486399216357D+00/                       
      DATA X(20,41)/  0.75623258989162996924D-01/                       
      DATA X(21,41)/  0.0                       /                       
      DATA W( 1,41)/  0.43061403581648876840D-02/                       
      DATA W( 2,41)/  0.99999387739059453385D-02/                       
      DATA W( 3,41)/  0.15644938407818588531D-01/                       
      DATA W( 4,41)/  0.21201063368779553076D-01/                       
      DATA W( 5,41)/  0.26635899207110445468D-01/                       
      DATA W( 6,41)/  0.31918211731699281787D-01/                       
      DATA W( 7,41)/  0.37017716703507988435D-01/                       
      DATA W( 8,41)/  0.41905195195909689429D-01/                       
      DATA W( 9,41)/  0.46552648369014342061D-01/                       
      DATA W(10,41)/  0.50933454294617494781D-01/                       
      DATA W(11,41)/  0.55022519242578741880D-01/                       
      DATA W(12,41)/  0.58796420949871944991D-01/                       
      DATA W(13,41)/  0.62233542580966316472D-01/                       
      DATA W(14,41)/  0.65314196453527410436D-01/                       
      DATA W(15,41)/  0.68020736760876766736D-01/                       
      DATA W(16,41)/  0.70337660620817497482D-01/                       
      DATA W(17,41)/  0.72251696861023073396D-01/                       
      DATA W(18,41)/  0.73751882027223469939D-01/                       
      DATA W(19,41)/  0.74829623176221551891D-01/                       
      DATA W(20,41)/  0.75478747092715824027D-01/                       
      DATA W(21,41)/  0.75695535647298372319D-01/                       
      DATA X( 1,42)/  0.99839961899006241502D+00/                       
      DATA X( 2,42)/  0.99157728834086091979D+00/                       
      DATA X( 3,42)/  0.97934250806374819371D+00/                       
      DATA X( 4,42)/  0.96175936533820448875D+00/                       
      DATA X( 5,42)/  0.93892355735498817853D+00/                       
      DATA X( 6,42)/  0.91095972490412745258D+00/                       
      DATA X( 7,42)/  0.87802056981217274271D+00/                       
      DATA X( 8,42)/  0.84028598326181690093D+00/                       
      DATA X( 9,42)/  0.79796205325548741323D+00/                       
      DATA X(10,42)/  0.75127993568948048957D+00/                       
      DATA X(11,42)/  0.70049459055617121374D+00/                       
      DATA X(12,42)/  0.64588338886924783396D+00/                       
      DATA X(13,42)/  0.58774459748510932284D+00/                       
      DATA X(14,42)/  0.52639574993119228759D+00/                       
      DATA X(15,42)/  0.46217191207042192976D+00/                       
      DATA X(16,42)/  0.39542385204297505768D+00/                       
      DATA X(17,42)/  0.32651612446541151220D+00/                       
      DATA X(18,42)/  0.25582507934287908397D+00/                       
      DATA X(19,42)/  0.18373680656485455085D+00/                       
      DATA X(20,42)/  0.11064502720851986835D+00/                       
      DATA X(21,42)/  0.36948943165351775813D-01/                       
      DATA W( 1,42)/  0.41059986046490846106D-02/                       
      DATA W( 2,42)/  0.95362203017485024118D-02/                       
      DATA W( 3,42)/  0.14922443697357494145D-01/                       
      DATA W( 4,42)/  0.20227869569052644757D-01/                       
      DATA W( 5,42)/  0.25422959526113047887D-01/                       
      DATA W( 6,42)/  0.30479240699603468363D-01/                       
      DATA W( 7,42)/  0.35369071097592110833D-01/                       
      DATA W( 8,42)/  0.40065735180692261761D-01/                       
      DATA W( 9,42)/  0.44543577771965877874D-01/                       
      DATA W(10,42)/  0.48778140792803245027D-01/                       
      DATA W(11,42)/  0.52746295699174070344D-01/                       
      DATA W(12,42)/  0.56426369358018381646D-01/                       
      DATA W(13,42)/  0.59798262227586654313D-01/                       
      DATA W(14,42)/  0.62843558045002576409D-01/                       
      DATA W(15,42)/  0.65545624364908978927D-01/                       
      DATA W(16,42)/  0.67889703376521944855D-01/                       
      DATA W(17,42)/  0.69862992492594159766D-01/                       
      DATA W(18,42)/  0.71454714265170982922D-01/                       
      DATA W(19,42)/  0.72656175243804104888D-01/                       
      DATA W(20,42)/  0.73460813453467528264D-01/                       
      DATA W(21,42)/  0.73864234232172879996D-01/                       
      DATA X( 1,43)/  0.99847233224250771352D+00/                       
      DATA X( 2,43)/  0.99195955759324414642D+00/                       
      DATA X( 3,43)/  0.98027822098025533151D+00/                       
      DATA X( 4,43)/  0.96348661301407999341D+00/                       
      DATA X( 5,43)/  0.94167195684763786182D+00/                       
      DATA X( 6,43)/  0.91494790720613872946D+00/                       
      DATA X( 7,43)/  0.88345376521861686334D+00/                       
      DATA X( 8,43)/  0.84735371620931504900D+00/                       
      DATA X( 9,43)/  0.80683596413693863528D+00/                       
      DATA X(10,43)/  0.76211174719495512146D+00/                       
      DATA X(11,43)/  0.71341423526895705485D+00/                       
      DATA X(12,43)/  0.66099731375149813317D+00/                       
      DATA X(13,43)/  0.60513425963960093573D+00/                       
      DATA X(14,43)/  0.54611631666008471914D+00/                       
      DATA X(15,43)/  0.48425117678573472407D+00/                       
      DATA X(16,43)/  0.41986137602926925249D+00/                       
      DATA X(17,43)/  0.35328261286430380665D+00/                       
      DATA X(18,43)/  0.28486199803291362711D+00/                       
      DATA X(19,43)/  0.21495624486051820901D+00/                       
      DATA X(20,43)/  0.14392980951071331077D+00/                       
      DATA X(21,43)/  0.72152990874586235422D-01/                       
      DATA X(22,43)/  0.0                       /                       
      DATA W( 1,43)/  0.39194902538441272830D-02/                       
      DATA W( 2,43)/  0.91039966374014033189D-02/                       
      DATA W( 3,43)/  0.14248756431576486109D-01/                       
      DATA W( 4,43)/  0.19319901423683900396D-01/                       
      DATA W( 5,43)/  0.24290456613838815902D-01/                       
      DATA W( 6,43)/  0.29134413261498494916D-01/                       
      DATA W( 7,43)/  0.33826492086860292345D-01/                       
      DATA W( 8,43)/  0.38342222194132657572D-01/                       
      DATA W( 9,43)/  0.42658057197982083764D-01/                       
      DATA W(10,43)/  0.46751494754346580011D-01/                       
      DATA W(11,43)/  0.50601192784390156524D-01/                       
      DATA W(12,43)/  0.54187080318881786863D-01/                       
      DATA W(13,43)/  0.57490461956910519428D-01/                       
      DATA W(14,43)/  0.60494115249991294520D-01/                       
      DATA W(15,43)/  0.63182380449396112326D-01/                       
      DATA W(16,43)/  0.65541242126322797491D-01/                       
      DATA W(17,43)/  0.67558402229365169192D-01/                       
      DATA W(18,43)/  0.69223344193656684282D-01/                       
      DATA W(19,43)/  0.70527387765085028126D-01/                       
      DATA W(20,43)/  0.71463734252514141298D-01/                       
      DATA W(21,43)/  0.72027501971421974345D-01/                       
      DATA W(22,43)/  0.72215751693798987977D-01/                       
      DATA X( 1,44)/  0.99854020063677422494D+00/                       
      DATA X( 2,44)/  0.99231639213851580848D+00/                       
      DATA X( 3,44)/  0.98115183307791396666D+00/                       
      DATA X( 4,44)/  0.96509965042249313939D+00/                       
      DATA X( 5,44)/  0.94423950911819409920D+00/                       
      DATA X( 6,44)/  0.91867525998417577432D+00/                       
      DATA X( 7,44)/  0.88853423828604320234D+00/                       
      DATA X( 8,44)/  0.85396659500471037873D+00/                       
      DATA X( 9,44)/  0.81514453964513501049D+00/                       
      DATA X(10,44)/  0.77226147924875589902D+00/                       
      DATA X(11,44)/  0.72553105366071700261D+00/                       
      DATA X(12,44)/  0.67518607066612236533D+00/                       
      DATA X(13,44)/  0.62147734590357584780D+00/                       
      DATA X(14,44)/  0.56467245318547076842D+00/                       
      DATA X(15,44)/  0.50505439138820231798D+00/                       
      DATA X(16,44)/  0.44292017452541148383D+00/                       
      DATA X(17,44)/  0.37857935201470713251D+00/                       
      DATA X(18,44)/  0.31235246650278581224D+00/                       
      DATA X(19,44)/  0.24456945692820125151D+00/                       
      DATA X(20,44)/  0.17556801477551678575D+00/                       
      DATA X(21,44)/  0.10569190170865324712D+00/                       
      DATA X(22,44)/  0.35289236964135359058D-01/                       
      DATA W( 1,44)/  0.37454048031127775152D-02/                       
      DATA W( 2,44)/  0.87004813675248441226D-02/                       
      DATA W( 3,44)/  0.13619586755579985520D-01/                       
      DATA W( 4,44)/  0.18471481736814749172D-01/                       
      DATA W( 5,44)/  0.23231481902019210629D-01/                       
      DATA W( 6,44)/  0.27875782821281010081D-01/                       
      DATA W( 7,44)/  0.32381222812069820881D-01/                       
      DATA W( 8,44)/  0.36725347813808873643D-01/                       
      DATA W( 9,44)/  0.40886512310346218908D-01/                       
      DATA W(10,44)/  0.44843984081970031446D-01/                       
      DATA W(11,44)/  0.48578046448352037528D-01/                       
      DATA W(12,44)/  0.52070096091704461881D-01/                       
      DATA W(13,44)/  0.55302735563728052549D-01/                       
      DATA W(14,44)/  0.58259859877595495334D-01/                       
      DATA W(15,44)/  0.60926736701561968039D-01/                       
      DATA W(16,44)/  0.63290079733203854950D-01/                       
      DATA W(17,44)/  0.65338114879181434984D-01/                       
      DATA W(18,44)/  0.67060638906293652396D-01/                       
      DATA W(19,44)/  0.68449070269366660985D-01/                       
      DATA W(20,44)/  0.69496491861572578037D-01/                       
      DATA W(21,44)/  0.70197685473558212587D-01/                       
      DATA W(22,44)/  0.70549157789354068811D-01/                       
      DATA X( 1,45)/  0.99860364518193663816D+00/                       
      DATA X( 2,45)/  0.99264999844720374175D+00/                       
      DATA X( 3,45)/  0.98196871503454056824D+00/                       
      DATA X( 4,45)/  0.96660831039689460474D+00/                       
      DATA X( 5,45)/  0.94664169099562906178D+00/                       
      DATA X( 6,45)/  0.92216393671900038810D+00/                       
      DATA X( 7,45)/  0.89329167175324173846D+00/                       
      DATA X( 8,45)/  0.86016247596066422534D+00/                       
      DATA X( 9,45)/  0.82293422050208633704D+00/                       
      DATA X(10,45)/  0.78178431259390629131D+00/                       
      DATA X(11,45)/  0.73690884894549035262D+00/                       
      DATA X(12,45)/  0.68852168077120052523D+00/                       
      DATA X(13,45)/  0.63685339445322335927D+00/                       
      DATA X(14,45)/  0.58215021256935318668D+00/                       
      DATA X(15,45)/  0.52467282046291606709D+00/                       
      DATA X(16,45)/  0.46469512391963509858D+00/                       
      DATA X(17,45)/  0.40250294385854191408D+00/                       
      DATA X(18,45)/  0.33839265425060216164D+00/                       
      DATA X(19,45)/  0.27266976975237756061D+00/                       
      DATA X(20,45)/  0.20564748978326374572D+00/                       
      DATA X(21,45)/  0.13764520598325302876D+00/                       
      DATA X(22,45)/  0.68986980163144172490D-01/                       
      DATA X(23,45)/  0.0                       /                       
      DATA W( 1,45)/  0.35826631552835589311D-02/                       
      DATA W( 2,45)/  0.83231892962182416457D-02/                       
      DATA W( 3,45)/  0.13031104991582784321D-01/                       
      DATA W( 4,45)/  0.17677535257937590617D-01/                       
      DATA W( 5,45)/  0.22239847550578732394D-01/                       
      DATA W( 6,45)/  0.26696213967577664806D-01/                       
      DATA W( 7,45)/  0.31025374934515467163D-01/                       
      DATA W( 8,45)/  0.35206692201609016248D-01/                       
      DATA W( 9,45)/  0.39220236729302447564D-01/                       
      DATA W(10,45)/  0.43046880709164971152D-01/                       
      DATA W(11,45)/  0.46668387718373365268D-01/                       
      DATA W(12,45)/  0.50067499237952029799D-01/                       
      DATA W(13,45)/  0.53228016731268951946D-01/                       
      DATA W(14,45)/  0.56134878759786476644D-01/                       
      DATA W(15,45)/  0.58774232718841738574D-01/                       
      DATA W(16,45)/  0.61133500831066522502D-01/                       
      DATA W(17,45)/  0.63201440073819937750D-01/                       
      DATA W(18,45)/  0.64968195750723430854D-01/                       
      DATA W(19,45)/  0.66425348449842528083D-01/                       
      DATA W(20,45)/  0.67565954163607536271D-01/                       
      DATA W(21,45)/  0.68384577378669674532D-01/                       
      DATA W(22,45)/  0.68877316977661322882D-01/                       
      DATA W(23,45)/  0.69041824829232020111D-01/                       
      DATA X( 1,46)/  0.99866304213381798113D+00/                       
      DATA X( 2,46)/  0.99296234890617436407D+00/                       
      DATA X( 3,46)/  0.98273366980416686348D+00/                       
      DATA X( 4,46)/  0.96802139185399194274D+00/                       
      DATA X( 5,46)/  0.94889236344608979562D+00/                       
      DATA X( 6,46)/  0.92543379880675395098D+00/                       
      DATA X( 7,46)/  0.89775271153394196570D+00/                       
      DATA X( 8,46)/  0.86597539486685806292D+00/                       
      DATA X( 9,46)/  0.83024683706606605303D+00/                       
      DATA X(10,46)/  0.79073005707527425519D+00/                       
      DATA X(11,46)/  0.74760535961566605400D+00/                       
      DATA X(12,46)/  0.70106951202040569751D+00/                       
      DATA X(13,46)/  0.65133484620199771511D+00/                       
      DATA X(14,46)/  0.59862828971271515318D+00/                       
      DATA X(15,46)/  0.54319033026180263527D+00/                       
      DATA X(16,46)/  0.48527391838816466277D+00/                       
      DATA X(17,46)/  0.42514331328282839732D+00/                       
      DATA X(18,46)/  0.36307287702099571012D+00/                       
      DATA X(19,46)/  0.29934582270187001548D+00/                       
      DATA X(20,46)/  0.23425292220626976863D+00/                       
      DATA X(21,46)/  0.16809117946710352861D+00/                       
      DATA X(22,46)/  0.10116247530558423952D+00/                       
      DATA X(23,46)/  0.33772190016052041520D-01/                       
      DATA W( 1,46)/  0.34303008681070482860D-02/                       
      DATA W( 2,46)/  0.79698982297246224516D-02/                       
      DATA W( 3,46)/  0.12479883770988684207D-01/                       
      DATA W( 4,46)/  0.16933514007836238046D-01/                       
      DATA W( 5,46)/  0.21309998754136501054D-01/                       
      DATA W( 6,46)/  0.25589286397130010635D-01/                       
      DATA W( 7,46)/  0.29751829552202755799D-01/                       
      DATA W( 8,46)/  0.33778627999106896521D-01/                       
      DATA W( 9,46)/  0.37651305357386071328D-01/                       
      DATA W(10,46)/  0.41352190109678729704D-01/                       
      DATA W(11,46)/  0.44864395277318126767D-01/                       
      DATA W(12,46)/  0.48171895101712200530D-01/                       
      DATA W(13,46)/  0.51259598007143021335D-01/                       
      DATA W(14,46)/  0.54113415385856754492D-01/                       
      DATA W(15,46)/  0.56720325843991235817D-01/                       
      DATA W(16,46)/  0.59068434595546314808D-01/                       
      DATA W(17,46)/  0.61147027724650481015D-01/                       
      DATA W(18,46)/  0.62946621064394508179D-01/                       
      DATA W(19,46)/  0.64459003467139069588D-01/                       
      DATA W(20,46)/  0.65677274267781207379D-01/                       
      DATA W(21,46)/  0.66595874768454887376D-01/                       
      DATA W(22,46)/  0.67210613600678175862D-01/                       
      DATA W(23,46)/  0.67518685849036458820D-01/                       
      DATA X( 1,47)/  0.99871872858421210918D+00/                       
      DATA X( 2,47)/  0.99325521098776863469D+00/                       
      DATA X( 3,47)/  0.98345100307162370876D+00/                       
      DATA X( 4,47)/  0.96934678732656449715D+00/                       
      DATA X( 5,47)/  0.95100396925770844259D+00/                       
      DATA X( 6,47)/  0.92850269301236064820D+00/                       
      DATA X( 7,47)/  0.90194132943852535687D+00/                       
      DATA X( 8,47)/  0.87143601579689631694D+00/                       
      DATA X( 9,47)/  0.83712013989990212128D+00/                       
      DATA X(10,47)/  0.79914375416774194292D+00/                       
      DATA X(11,47)/  0.75767291844543863357D+00/                       
      DATA X(12,47)/  0.71288897340906430166D+00/                       
      DATA X(13,47)/  0.66498774739033272914D+00/                       
      DATA X(14,47)/  0.61417869995637360860D+00/                       
      DATA X(15,47)/  0.56068400593466419448D+00/                       
      DATA X(16,47)/  0.50473758386357791977D+00/                       
      DATA X(17,47)/  0.44658407310485570273D+00/                       
      DATA X(18,47)/  0.38647776408466713958D+00/                       
      DATA X(19,47)/  0.32468148633773590221D+00/                       
      DATA X(20,47)/  0.26146545921497457031D+00/                       
      DATA X(21,47)/  0.19710611027911180796D+00/                       
      DATA X(22,47)/  0.13188486655451489705D+00/                       
      DATA X(23,47)/  0.66086923916355675160D-01/                       
      DATA X(24,47)/  0.0                       /                       
      DATA W( 1,47)/  0.32874538425280148832D-02/                       
      DATA W( 2,47)/  0.76386162958488336141D-02/                       
      DATA W( 3,47)/  0.11962848464312320964D-01/                       
      DATA W( 4,47)/  0.16235333146433059671D-01/                       
      DATA W( 5,47)/  0.20436938147668427642D-01/                       
      DATA W( 6,47)/  0.24549211659658818538D-01/                       
      DATA W( 7,47)/  0.28554150700643386505D-01/                       
      DATA W( 8,47)/  0.32434235515184756768D-01/                       
      DATA W( 9,47)/  0.36172496584174951613D-01/                       
      DATA W(10,47)/  0.39752586122531003781D-01/                       
      DATA W(11,47)/  0.43158848648479538268D-01/                       
      DATA W(12,47)/  0.46376389086505911204D-01/                       
      DATA W(13,47)/  0.49391137747361169605D-01/                       
      DATA W(14,47)/  0.52189911780057144872D-01/                       
      DATA W(15,47)/  0.54760472781530225957D-01/                       
      DATA W(16,47)/  0.57091580293231540222D-01/                       
      DATA W(17,47)/  0.59173040942338875976D-01/                       
      DATA W(18,47)/  0.60995753008739645331D-01/                       
      DATA W(19,47)/  0.62551746220921662641D-01/                       
      DATA W(20,47)/  0.63834216605717030631D-01/                       
      DATA W(21,47)/  0.64837556238945726703D-01/                       
      DATA W(22,47)/  0.65557377766549740251D-01/                       
      DATA W(23,47)/  0.65990533588810474534D-01/                       
      DATA W(24,47)/  0.66135129623655479653D-01/                       
      DATA X( 1,48)/  0.99877100725242611860D+00/                       
      DATA X( 2,48)/  0.99353017226635075755D+00/                       
      DATA X( 3,48)/  0.98412458372282685774D+00/                       
      DATA X( 4,48)/  0.97059159254624725046D+00/                       
      DATA X( 5,48)/  0.95298770316043086072D+00/                       
      DATA X( 6,48)/  0.93138669070655433311D+00/                       
      DATA X( 7,48)/  0.90587913671556967282D+00/                       
      DATA X( 8,48)/  0.87657202027424788591D+00/                       
      DATA X( 9,48)/  0.84358826162439353071D+00/                       
      DATA X(10,48)/  0.80706620402944262708D+00/                       
      DATA X(11,48)/  0.76715903251574033925D+00/                       
      DATA X(12,48)/  0.72403413092381465467D+00/                       
      DATA X(13,48)/  0.67787237963266390521D+00/                       
      DATA X(14,48)/  0.62886739677651362400D+00/                       
      DATA X(15,48)/  0.57722472608397270382D+00/                       
      DATA X(16,48)/  0.52316097472223303368D+00/                       
      DATA X(17,48)/  0.46690290475095840454D+00/                       
      DATA X(18,48)/  0.40868648199071672992D+00/                       
      DATA X(19,48)/  0.34875588629216073816D+00/                       
      DATA X(20,48)/  0.28736248735545557674D+00/                       
      DATA X(21,48)/  0.22476379039468906122D+00/                       
      DATA X(22,48)/  0.16122235606889171806D+00/                       
      DATA X(23,48)/  0.97004699209462698930D-01/                       
      DATA X(24,48)/  0.32380170962869362033D-01/                       
      DATA W( 1,48)/  0.31533460523058386327D-02/                       
      DATA W( 2,48)/  0.73275539012762621024D-02/                       
      DATA W( 3,48)/  0.11477234579234539490D-01/                       
      DATA W( 4,48)/  0.15579315722943848728D-01/                       
      DATA W( 5,48)/  0.19616160457355527814D-01/                       
      DATA W( 6,48)/  0.23570760839324379141D-01/                       
      DATA W( 7,48)/  0.27426509708356948200D-01/                       
      DATA W( 8,48)/  0.31167227832798088902D-01/                       
      DATA W( 9,48)/  0.34777222564770438893D-01/                       
      DATA W(10,48)/  0.38241351065830706317D-01/                       
      DATA W(11,48)/  0.41545082943464749214D-01/                       
      DATA W(12,48)/  0.44674560856694280419D-01/                       
      DATA W(13,48)/  0.47616658492490474826D-01/                       
      DATA W(14,48)/  0.50359035553854474958D-01/                       
      DATA W(15,48)/  0.52890189485193667096D-01/                       
      DATA W(16,48)/  0.55199503699984162868D-01/                       
      DATA W(17,48)/  0.57277292100403215705D-01/                       
      DATA W(18,48)/  0.59114839698395635746D-01/                       
      DATA W(19,48)/  0.60704439165893880053D-01/                       
      DATA W(20,48)/  0.62039423159892663904D-01/                       
      DATA W(21,48)/  0.63114192286254025657D-01/                       
      DATA W(22,48)/  0.63924238584648186624D-01/                       
      DATA W(23,48)/  0.64466164435950082207D-01/                       
      DATA W(24,48)/  0.64737696812683922503D-01/                       
      DATA X( 1,49)/  0.99882015060663537936D+00/                       
      DATA X( 2,49)/  0.99378866194416779076D+00/                       
      DATA X( 3,49)/  0.98475789591421300436D+00/                       
      DATA X( 4,49)/  0.97176220090155538014D+00/                       
      DATA X( 5,49)/  0.95485365867413723356D+00/                       
      DATA X( 6,49)/  0.93410029475581014906D+00/                       
      DATA X( 7,49)/  0.90958565582807328521D+00/                       
      DATA X( 8,49)/  0.88140844557300891004D+00/                       
      DATA X( 9,49)/  0.84968211984416570103D+00/                       
      DATA X(10,49)/  0.81453442735985543154D+00/                       
      DATA X(11,49)/  0.77610689434544663502D+00/                       
      DATA X(12,49)/  0.73455425423740269621D+00/                       
      DATA X(13,49)/  0.69004382442513211350D+00/                       
      DATA X(14,49)/  0.64275483241923766406D+00/                       
      DATA X(15,49)/  0.59287769410890071246D+00/                       
      DATA X(16,49)/  0.54061324699172606656D+00/                       
      DATA X(17,49)/  0.48617194145249204218D+00/                       
      DATA X(18,49)/  0.42977299334157652466D+00/                       
      DATA X(19,49)/  0.37164350126228488886D+00/                       
      DATA X(20,49)/  0.31201753211974876221D+00/                       
      DATA X(21,49)/  0.25113517861257727351D+00/                       
      DATA X(22,49)/  0.18924159246181358649D+00/                       
      DATA X(23,49)/  0.12658599726967205107D+00/                       
      DATA X(24,49)/  0.63420684982686786029D-01/                       
      DATA X(25,49)/  0.0                       /                       
      DATA W( 1,49)/  0.30272789889229050775D-02/                       
      DATA W( 2,49)/  0.70350995900864514735D-02/                       
      DATA W( 3,49)/  0.11020551031593580498D-01/                       
      DATA W( 4,49)/  0.14962144935624651030D-01/                       
      DATA W( 5,49)/  0.18843595853089458444D-01/                       
      DATA W( 6,49)/  0.22649201587446676499D-01/                       
      DATA W( 7,49)/  0.26363618927066016961D-01/                       
      DATA W( 8,49)/  0.29971884620583825351D-01/                       
      DATA W( 9,49)/  0.33459466791622174342D-01/                       
      DATA W(10,49)/  0.36812320963000689819D-01/                       
      DATA W(11,49)/  0.40016945766373021369D-01/                       
      DATA W(12,49)/  0.43060436981259597988D-01/                       
      DATA W(13,49)/  0.45930539355595853542D-01/                       
      DATA W(14,49)/  0.48615695887828240278D-01/                       
      DATA W(15,49)/  0.51105094330144590675D-01/                       
      DATA W(16,49)/  0.53388710708258968528D-01/                       
      DATA W(17,49)/  0.55457349674803588690D-01/                       
      DATA W(18,49)/  0.57302681530187475485D-01/                       
      DATA W(19,49)/  0.58917275760027266025D-01/                       
      DATA W(20,49)/  0.60294630953152017303D-01/                       
      DATA W(21,49)/  0.61429200979192936297D-01/                       
      DATA W(22,49)/  0.62316417320057267401D-01/                       
      DATA W(23,49)/  0.62952707465195699474D-01/                       
      DATA W(24,49)/  0.63335509296491748591D-01/                       
      DATA W(25,49)/  0.63463281404790597718D-01/                       
      DATA X( 1,50)/  0.99886640442007105019D+00/                       
      DATA X( 2,50)/  0.99403196943209071259D+00/                       
      DATA X( 3,50)/  0.98535408404800588231D+00/                       
      DATA X( 4,50)/  0.97286438510669207371D+00/                       
      DATA X( 5,50)/  0.95661095524280794300D+00/                       
      DATA X( 6,50)/  0.93665661894487793378D+00/                       
      DATA X( 7,50)/  0.91307855665579189309D+00/                       
      DATA X( 8,50)/  0.88596797952361304864D+00/                       
      DATA X( 9,50)/  0.85542976942994608461D+00/                       
      DATA X(10,50)/  0.82158207085933594836D+00/                       
      DATA X(11,50)/  0.78455583290039926391D+00/                       
      DATA X(12,50)/  0.74449430222606853826D+00/                       
      DATA X(13,50)/  0.70155246870682225109D+00/                       
      DATA X(14,50)/  0.65589646568543936078D+00/                       
      DATA X(15,50)/  0.60770292718495023918D+00/                       
      DATA X(16,50)/  0.55715830451465005432D+00/                       
      DATA X(17,50)/  0.50445814490746420165D+00/                       
      DATA X(18,50)/  0.44980633497403878915D+00/                       
      DATA X(19,50)/  0.39341431189756512739D+00/                       
      DATA X(20,50)/  0.33550024541943735684D+00/                       
      DATA X(21,50)/  0.27628819377953199033D+00/                       
      DATA X(22,50)/  0.21600723687604175685D+00/                       
      DATA X(23,50)/  0.15489058999814590207D+00/                       
      DATA X(24,50)/  0.93174701560086140854D-01/                       
      DATA X(25,50)/  0.31098338327188876112D-01/                       
      DATA W( 1,50)/  0.29086225531551409584D-02/                       
      DATA W( 2,50)/  0.67597991957454015028D-02/                       
      DATA W( 3,50)/  0.10590548383650969264D-01/                       
      DATA W( 4,50)/  0.14380822761485574419D-01/                       
      DATA W( 5,50)/  0.18115560713489390351D-01/                       
      DATA W( 6,50)/  0.21780243170124792982D-01/                       
      DATA W( 7,50)/  0.25360673570012390440D-01/                       
      DATA W( 8,50)/  0.28842993580535198030D-01/                       
      DATA W( 9,50)/  0.32213728223578016648D-01/                       
      DATA W(10,50)/  0.35459835615146154161D-01/                       
      DATA W(11,50)/  0.38568756612587675245D-01/                       
      DATA W(12,50)/  0.41528463090147697422D-01/                       
      DATA W(13,50)/  0.44327504338803275492D-01/                       
      DATA W(14,50)/  0.46955051303948432966D-01/                       
      DATA W(15,50)/  0.49400938449466314921D-01/                       
      DATA W(16,50)/  0.51655703069581138490D-01/                       
      DATA W(17,50)/  0.53710621888996246523D-01/                       
      DATA W(18,50)/  0.55557744806212517624D-01/                       
      DATA W(19,50)/  0.57189925647728383723D-01/                       
      DATA W(20,50)/  0.58600849813222445835D-01/                       
      DATA W(21,50)/  0.59785058704265457510D-01/                       
      DATA W(22,50)/  0.60737970841770216032D-01/                       
      DATA W(23,50)/  0.61455899590316663756D-01/                       
      DATA W(24,50)/  0.61936067420683243384D-01/                       
      DATA W(25,50)/  0.62176616655347262321D-01/                       
      IF(NP.LT.1.OR.NP.GT.50) GO TO 20                                  
      ICON=0                                                            
      NH=(NP+1)/2                                                       
      DO 10 I=1,NH                                                      
      PT(I)=X(I,NP)                                                     
      WT(I)=W(I,NP)                                                     
      PT(NP+1-I)=-PT(I)                                                 
      WT(NP+1-I)=WT(I)                                                  
   10 CONTINUE                                                          
      RETURN                                                            
   20 ICON=30000                                                        
      RETURN                                                            
      END                                                               
C#NUMPAC#DJN                 REVISED ON 1989-01-27                      
      FUNCTION DJN(N,X)                                                 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DATA EPS,XMAX,XMIN/1.D-16,3.53D15,-180.218D0/                     
      DATA RP/0.56418958354775628695D0/                                 
      IF(N.LT.0.OR.X.LT.0.D0.OR.X.GE.XMAX) GO TO 210                    
      IF(N.GT.0.AND.X.EQ.0.D0) GO TO 200                                
      IF(N.EQ.0.AND.X.EQ.0.D0) GO TO 190                                
      IF(N.GE.5) GO TO 60                                               
      GO TO (10,20,30,40,50),N+1                                        
   10 DJN=DJ0(X)                                                        
      RETURN                                                            
   20 DJN=DJ1(X)                                                        
      RETURN                                                            
   30 DJN=DJ2(X)                                                        
      RETURN                                                            
   40 DJN=DJ3(X)                                                        
      RETURN                                                            
   50 DJN=DJ4(X)                                                        
      RETURN                                                            
   60 FN=N                                                              
      XH=X*0.5D0                                                        
      XX=-XH*XH                                                         
      IF(FN+1.D0+XX.LT.0.D0) GO TO 100                                  
C     POWER SERIES                                                      
      IF(FN.GT.56.D0) GO TO 70                                          
      ALT=XH**N/DFCTRL(N)                                               
      IF(ALT.EQ.0.D0) GO TO 200                                         
      GO TO 80                                                          
   70 ALLT=DLOG(XH)*FN-DLGAMA(FN+1.D0)                                  
      IF(ALLT.LT.XMIN) GO TO 200                                        
      ALT=DEXP(ALLT)                                                    
   80 T=1.D0                                                            
      S=1.D0                                                            
      D=1.D0                                                            
   90 T=T*XX/((D+FN)*D)                                                 
      S=T+S                                                             
      D=D+1.D0                                                          
      IF(DABS(T).GE.EPS) GO TO 90                                       
      DJN=S*ALT                                                         
      RETURN                                                            
  100 RX=DSQRT(X)                                                       
      IF(X.LT.18.D0.OR.X.LT.FN*FN*0.55D0) GO TO 140                     
C     ASYMPTOTIC EXPANSION                                              
      R=X*8.D0                                                          
      T=RP                                                              
      P=T                                                               
      Q=0.D0                                                            
      AMU=FN*FN*4.D0                                                    
      F=-1.D0                                                           
  110 F=F+2.D0                                                          
      T=(AMU-(F+F-1.D0)**2)*T/(F*R)                                     
      Q=Q+T                                                             
      T=((F+F+1.D0)**2-AMU)*T/((F+1.D0)*R)                              
      P=P+T                                                             
      IF(DABS(T).GE.EPS) GO TO 110                                      
      S=DSIN(X)                                                         
      C=DCOS(X)                                                         
      SG=1-MOD(N/2,2)*2                                                 
      IF(MOD(N,2).EQ.1) GO TO 120                                       
      CS=C                                                              
      CD=S                                                              
      GO TO 130                                                         
  120 CS=S                                                              
      CD=-C                                                             
  130 DJN=((P+Q)*CS+(P-Q)*CD)*SG/RX                                     
      RETURN                                                            
C     BACKWARD RECURRENCE                                               
  140 IF(X.GE.10.D0) GO TO 150                                          
      U=4.7D0*X+43.D0                                                   
      AM=2.1D0*X+17.D0                                                  
      GO TO 170                                                         
  150 IF(X.GT.100.D0) GO TO 160                                         
      U=1.83D0*X+71.D0                                                  
      AM=1.26D0*X+30.D0                                                 
      GO TO 170                                                         
  160 IF(X.GT.200.D0) GO TO 220                                         
      U=1.5D0*X+116.D0                                                  
      AM=1.14D0*X+42.D0                                                 
  170 IF(FN.GE.U) GO TO 200                                             
      D=DINT((DMAX1(FN-(AM+X)*0.5D0,0.D0)+AM+1.D0)*0.5D0)*2.D0          
      P=0.D0                                                            
      T=0.D0                                                            
      S=1.D-75                                                          
      F=2.D0/X                                                          
      L=1                                                               
  180 IF(L.GT.0) P=S+P                                                  
      R=D*F*S-T                                                         
      D=D-1.D0                                                          
      IF(D.EQ.FN) Y=R                                                   
      T=S                                                               
      S=R                                                               
      L=-L                                                              
      IF(D.GT.0.D0) GO TO 180                                           
      DJN=Y/(R+P+P)                                                     
      RETURN                                                            
  190 DJN=1.D0                                                          
      RETURN                                                            
C     UNDERFLOW                                                         
  200 DJN=0.D0                                                          
      RETURN                                                            
C     ARGUMENT ERROR                                                    
  210 DJN=0.D0                                                          
cc      CALL MGIDD(5HDJN  ,N,X,DJN                                        
cc     *,36HARG1 OR ARG2 LT 0 OR ARG2 GT 3.53D15,9)                       
      RETURN                                                            
  220 DJN=0.D0                                                          
cc      CALL MGIDD(5HDJN  ,N,X,DJN                                        
cc     *,40HARG1 GT SQRT(ARG2)*1.384 AND ARG2 GT 200,10)                  
      RETURN                                                            
      END                                                               
C#NUMPAC#DJ0                 REVISED ON 1984-11-30                      
      FUNCTION DJ0(X)                                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DATA BIG / 3.53D15/                                               
      DATA QP  / 7.85398163397448309616D-01/                            
      DATA A0  / 9.99999999999999999862D-01/                            
      DATA A1  /-2.49999999999999993074D-01/                            
      DATA A2  / 1.56249999999999428124D-02/                            
      DATA A3  /-4.34027777777594605992D-04/                            
      DATA A4  / 6.78168402747980046740D-06/                            
      DATA A5  /-6.78168399992931283849D-08/                            
      DATA A6  / 4.70950121196549501372D-10/                            
      DATA A7  /-2.40275166735379467327D-12/                            
      DATA A8  / 9.37404284847115370170D-15/                            
      DATA A9  /-2.75547709733212714672D-17/                            
      DATA B0  /-2.60051954901933437087D-01/                            
      DATA B1  /-3.39058958525936427299D-01/                            
      DATA B2  / 1.86535803871956067624D-01/                            
      DATA B3  / 2.95047563884387202935D-02/                            
      DATA B4  /-1.35025350162735993399D-02/                            
      DATA B5  /-9.83491879601822413012D-04/                            
      DATA B6  / 3.95446062777725161898D-04/                            
      DATA B7  / 1.75948602150204955635D-05/                            
      DATA B8  /-6.33925272136885327244D-06/                            
      DATA B9  /-1.96844606705338199035D-07/                            
      DATA B10 / 6.43233596548903618701D-08/                            
      DATA B11 / 1.50721201189975222410D-09/                            
      DATA B12 /-4.50747057270707242904D-10/                            
      DATA B13 /-8.26231997847312156351D-12/                            
      DATA B14 / 2.27616275675398483717D-12/                            
      DATA C0  /-1.77596771314338304428D-01/                            
      DATA C1  / 3.27579137591465086417D-01/                            
      DATA C2  / 5.60404718980226390698D-02/                            
      DATA C3  /-5.61486934744975453127D-02/                            
      DATA C4  /-1.70738759685108596284D-03/                            
      DATA C5  / 2.52021197018631416159D-03/                            
      DATA C6  / 1.12022146861180112407D-05/                            
      DATA C7  /-5.37950091820866235019D-05/                            
      DATA C8  / 2.13330063500786110354D-07/                            
      DATA C9  / 6.78109662379437531731D-07/                            
      DATA C10 /-4.88571389013903589446D-09/                            
      DATA C11 /-5.66336220352842176862D-09/                            
      DATA C12 / 4.79998732771847989924D-11/                            
      DATA C13 / 3.31502164620686196353D-11/                            
      DATA C14 /-2.94070926658281415257D-13/                            
      DATA D0  / 3.00079270519555596180D-01/                            
      DATA D1  / 4.68282348234589820892D-03/                            
      DATA D2  /-1.50374122651373875987D-01/                            
      DATA D3  / 6.39612989784324898769D-03/                            
      DATA D4  / 1.17901293571020962207D-02/                            
      DATA D5  /-5.93148973882069399566D-04/                            
      DATA D6  /-3.52849100230102466451D-04/                            
      DATA D7  / 1.72261259581817565520D-05/                            
      DATA D8  / 5.66074614061640975883D-06/                            
      DATA D9  /-2.57978931323497772224D-07/                            
      DATA D10 /-5.70714255001103155952D-08/                            
      DATA D11 / 2.40512263081095495454D-09/                            
      DATA D12 / 3.96491159793042848309D-10/                            
      DATA D13 /-1.51780840423917197272D-11/                            
      DATA D14 /-1.98457116272515869927D-12/                            
      DATA E0  / 9.10688034326506962475D+04/                            
      DATA E1  / 1.27539507434834602368D+05/                            
      DATA E2  / 4.68143748172517300040D+04/                            
      DATA E3  / 4.84077384573519787279D+03/                            
      DATA E4  / 9.85802550097924276429D+01/                            
      DATA F0  / 1.14440433534091350885D+06/                            
      DATA F1  / 1.60494388311640608551D+06/                            
      DATA F2  / 5.91362505047859871932D+05/                            
      DATA F3  / 6.19081071744797093901D+04/                            
      DATA F4  / 1.33517907943909263871D+03/                            
      DATA G0  /-2.28911141416091769260D+03/                            
      DATA G1  /-3.65129617996166549134D+03/                            
      DATA G2  /-1.57616402953763678108D+03/                            
      DATA G3  /-2.01139695408233120405D+02/                            
      DATA G4  /-5.45851732770385664951D+00/                            
      DATA H0  / 1.46503130506299285047D+05/                            
      DATA H1  / 2.34875201045700477225D+05/                            
      DATA H2  / 1.02725950947053047345D+05/                            
      DATA H3  / 1.36201116594189451125D+04/                            
      DATA H4  / 4.28239987601685748267D+02/                            
      IF(X.LT.0.D0) GO TO 60                                            
      IF(X.GE.8.D0) GO TO 50                                            
      ISWICH=X+1.0D0                                                    
      GO TO (10,10,20,20,30,30,40,40),ISWICH                            
   10 T=X*X                                                             
      DJ0=((((((((A9*T+A8)*T+A7)*T+A6)*T+A5)*T+A4)*T+A3)*               
     *T+A2)*T+A1)*T+A0                                                  
      RETURN                                                            
   20 T=X-3.D0                                                          
      DJ0=(((((((((((((B14*T+B13)*T+B12)*T+B11)*T+B10)*T+B9)*           
     *T+B8)*T+B7)*T+B6)*T+B5)*T+B4)*T+B3)*T+B2)*T+B1)*T+B0              
      RETURN                                                            
   30 T=X-5.D0                                                          
      DJ0=(((((((((((((C14*T+C13)*T+C12)*T+C11)*T+C10)*T+C9)*           
     *T+C8)*T+C7)*T+C6)*T+C5)*T+C4)*T+C3)*T+C2)*T+C1)*T+C0              
      RETURN                                                            
   40 T=X-7.D0                                                          
      DJ0=(((((((((((((D14*T+D13)*T+D12)*T+D11)*T+D10)*T+D9)*           
     *T+D8)*T+D7)*T+D6)*T+D5)*T+D4)*T+D3)*T+D2)*T+D1)*T+D0              
      RETURN                                                            
   50 IF(X.GT.BIG) GO TO 60                                             
      Y=8.D0/X                                                          
      T=Y*Y                                                             
      H0A=((((E4*T+E3)*T+E2)*T+E1)*T+E0)/                               
     *(((((T+F4)*T+F3)*T+F2)*T+F1)*T+F0)                                
      H0P=((((G4*T+G3)*T+G2)*T+G1)*T+G0)/                               
     *(((((T+H4)*T+H3)*T+H2)*T+H1)*T+H0)*Y                              
      DJ0=DCOS(H0P-QP+X)*DSQRT(Y*H0A)                                   
      RETURN                                                            
   60 DJ0=0.D0                                                          
cc      CALL MGDD(5HDJ0  ,X,DJ0 ,22HARG LT 0 OR GT 3.53D15,6)             
      RETURN                                                            
      END                                                               
C#NUMPAC#DJ1                 REVISED ON 1984-11-30                      
      FUNCTION DJ1(X)                                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DATA BIG / 3.53D15/                                               
      DATA QP  / 7.85398163397448309616D-01/                            
      DATA A0  / 4.99999999999999996383D-01/                            
      DATA A1  /-6.24999999999998613001D-02/                            
      DATA A2  / 2.60416666666578150139D-03/                            
      DATA A3  /-5.42534722200295867810D-05/                            
      DATA A4  / 6.78168400040012821731D-07/                            
      DATA A5  /-5.65140142652668401715D-09/                            
      DATA A6  / 3.36385029716467844879D-11/                            
      DATA A7  /-1.49980003131123197292D-13/                            
      DATA A8  / 4.95609337066131684756D-16/                            
      DATA B0  / 3.39058958525936458805D-01/                            
      DATA B1  /-3.73071607743912126170D-01/                            
      DATA B2  /-8.85142691653197029928D-02/                            
      DATA B3  / 5.40101400650941119770D-02/                            
      DATA B4  / 4.91745939807306046434D-03/                            
      DATA B5  /-2.37267637666381859883D-03/                            
      DATA B6  /-1.23164021932485703549D-04/                            
      DATA B7  / 5.07140217614141078414D-05/                            
      DATA B8  / 1.77160280605398176232D-06/                            
      DATA B9  /-6.43233579231124631716D-07/                            
      DATA B10 /-1.65814887767894328214D-08/                            
      DATA B11 / 5.40894969477222524312D-09/                            
      DATA B12 / 1.09111684639251778886D-10/                            
      DATA B13 /-3.18613153743493889088D-11/                            
      DATA B14 /-5.24166874680411486937D-13/                            
      DATA C0  /-3.27579137591465221526D-01/                            
      DATA C1  /-1.12080943796045278594D-01/                            
      DATA C2  / 1.68446080423507774121D-01/                            
      DATA C3  / 6.82955038740436092818D-03/                            
      DATA C4  /-1.26010598512041435713D-02/                            
      DATA C5  /-6.72132881168797693386D-05/                            
      DATA C6  / 3.76565066092217933665D-04/                            
      DATA C7  /-1.70664050729667364464D-06/                            
      DATA C8  /-6.10299267509887011424D-06/                            
      DATA C9  / 4.88571375164175788776D-08/                            
      DATA C10 / 6.23061277290522475749D-08/                            
      DATA C11 /-5.75997209084160501511D-10/                            
      DATA C12 /-4.38157844782505841595D-10/                            
      DATA C13 / 4.11655209936008770079D-12/                            
      DATA C14 / 2.21721548469535395227D-12/                            
      DATA D0  /-4.68282348234583296870D-03/                            
      DATA D1  / 3.00748245302747745139D-01/                            
      DATA D2  /-1.91883896935370617772D-02/                            
      DATA D3  /-4.71605174284081669737D-02/                            
      DATA D4  / 2.96574486954212449006D-03/                            
      DATA D5  / 2.11709460137865931674D-03/                            
      DATA D6  /-1.20582882586380446227D-04/                            
      DATA D7  /-4.52859691174992149085D-05/                            
      DATA D8  / 2.32181314634470517627D-06/                            
      DATA D9  / 5.70714241401770307424D-07/                            
      DATA D10 /-2.64607740621387775965D-08/                            
      DATA D11 /-4.75788206801760283791D-09/                            
      DATA D12 / 2.00802907113180031958D-10/                            
      DATA D13 / 2.77800518104936964502D-11/                            
      DATA D14 /-1.07353087329324334075D-12/                            
      DATA E0  / 4.28548556943758496680D+03/                            
      DATA E1  / 7.53167724919758395134D+03/                            
      DATA E2  / 3.69887485576240589863D+03/                            
      DATA E3  / 5.76643772967497887526D+02/                            
      DATA E4  / 2.42480433265881372760D+01/                            
      DATA E5  / 1.72185745103162331665D-01/                            
      DATA F0  / 5.38529999280407546520D+04/                            
      DATA F1  / 9.43303027397007158896D+04/                            
      DATA F2  / 4.59333379160879599250D+04/                            
      DATA F3  / 6.98495917314242939982D+03/                            
      DATA F4  / 2.67220025782266488678D+02/                            
      DATA G0  / 8.92150268629996415700D+03/                            
      DATA G1  / 1.36398688858717213460D+04/                            
      DATA G2  / 5.55180423676669770378D+03/                            
      DATA G3  / 6.48920811086677200345D+02/                            
      DATA G4  / 1.51548413247039541977D+01/                            
      DATA H0  / 1.90325390641066217999D+05/                            
      DATA H1  / 2.92284922040317697279D+05/                            
      DATA H2  / 1.20390573908313954794D+05/                            
      DATA H3  / 1.46006356254225211587D+04/                            
      DATA H4  / 4.00160109244578676504D+02/                            
      IF(X.LT.0.D0) GO TO 60                                            
      IF(X.GE.8.D0) GO TO 50                                            
      ISWICH=X+1.0D0                                                    
      GO TO (10,10,20,20,30,30,40,40),ISWICH                            
   10 T=X*X                                                             
      DJ1=((((((((A8*T+A7)*T+A6)*T+A5)*T+A4)*T+A3)*                     
     *T+A2)*T+A1)*T+A0)*X                                               
      RETURN                                                            
   20 T=X-3.D0                                                          
      DJ1=(((((((((((((B14*T+B13)*T+B12)*T+B11)*T+B10)*T+B9)*           
     *T+B8)*T+B7)*T+B6)*T+B5)*T+B4)*T+B3)*T+B2)*T+B1)*T+B0              
      RETURN                                                            
   30 T=X-5.D0                                                          
      DJ1=(((((((((((((C14*T+C13)*T+C12)*T+C11)*T+C10)*T+C9)*           
     *T+C8)*T+C7)*T+C6)*T+C5)*T+C4)*T+C3)*T+C2)*T+C1)*T+C0              
      RETURN                                                            
   40 T=X-7.D0                                                          
      DJ1=(((((((((((((D14*T+D13)*T+D12)*T+D11)*T+D10)*T+D9)*           
     *T+D8)*T+D7)*T+D6)*T+D5)*T+D4)*T+D3)*T+D2)*T+D1)*T+D0              
      RETURN                                                            
   50 IF(X.GT.BIG) GO TO 60                                             
      Y=8.D0/X                                                          
      T=Y*Y                                                             
      H1A=(((((E5*T+E4)*T+E3)*T+E2)*T+E1)*T+E0)/                        
     *(((((T+F4)*T+F3)*T+F2)*T+F1)*T+F0)                                
      H1P=((((G4*T+G3)*T+G2)*T+G1)*T+G0)/                               
     *(((((T+H4)*T+H3)*T+H2)*T+H1)*T+H0)*Y                              
      DJ1=DSIN(H1P-QP+X)*DSQRT(Y*H1A)                                   
      RETURN                                                            
   60 DJ1=0.D0                                                          
cc      CALL MGDD(5HDJ1  ,X,DJ1 ,22HARG LT 0 OR GT 3.53D15,6)             
      RETURN                                                            
      END                                                               
C#NUMPAC#DJ2                 REVISED ON 1988-12-17                      
      FUNCTION DJ2(X)                                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DATA BIG / 3.53D15/                                               
      DATA QP  / 7.85398163397448309616D-01/                            
      DATA A0/  1.24999999999999999850D-01/                             
      DATA A1/ -1.04166666666666608017D-02/                             
      DATA A2/  3.25520833333295308420D-04/                             
      DATA A3/ -5.42534722212680006089D-06/                             
      DATA A4/  5.65140334443668074661D-08/                             
      DATA A5/ -4.03671582625879396994D-10/                             
      DATA A6/  2.10242068188985477036D-12/                             
      DATA A7/ -8.33427531430941569015D-15/                             
      DATA A8/  2.49018592208268580443D-17/                             
      DATA B0/  4.86091260585891076428D-01/                             
      DATA B1/  1.49981181353423822280D-02/                             
      DATA B2/ -1.37525036518637868182D-01/                             
      DATA B3/ -9.83491879614568161300D-03/                             
      DATA B4/  1.02242287508915107066D-02/                             
      DATA B5/  4.94476383587244257391D-04/                             
      DATA B6/ -3.14550245395826188643D-04/                             
      DATA B7/ -1.07507846788630492121D-05/                             
      DATA B8/  5.23896275040840193719D-06/                             
      DATA B9/  1.34785163249016553819D-07/                             
      DATA B10/ -5.46912098128292907181D-08/                            
      DATA B11/ -1.11146346238454581615D-09/                            
      DATA B12/  3.91575966580219614513D-10/                            
      DATA B13/  6.41267583393998950065D-12/                            
      DATA B14/ -2.01019323420397766983D-12/                            
      DATA C0/  4.65651162777522156182D-02/                             
      DATA C1/ -3.46205184102565988159D-01/                             
      DATA C2/  1.50631695736006705374D-02/                             
      DATA C3/  4.46597853351349264516D-02/                             
      DATA C4/ -1.03525471575836246070D-03/                             
      DATA C5/ -1.99856882291438206604D-03/                             
      DATA C6/  3.50951822982765484565D-05/                             
      DATA C7/  4.38528735974776020174D-05/                             
      DATA C8/ -6.66100022041917063041D-07/                             
      DATA C9/ -5.68012852638644331202D-07/                             
      DATA C10/  7.78881100336825263363D-09/                            
      DATA C11/  4.85239213130467278990D-09/                            
      DATA C12/ -6.10748024931186410871D-11/                            
      DATA C13/ -2.89206704603501388844D-11/                            
      DATA C14/  3.36709046153362593459D-13/                            
      DATA D0/ -3.01417220085940119864D-01/                             
      DATA D1/  8.14363822564941382554D-02/                             
      DATA D2/  1.32588981919100422517D-01/                             
      DATA D3/ -1.73298290584935191164D-02/                             
      DATA D4/ -9.38081665713998946123D-03/                             
      DATA D5/  8.53845617152432731402D-04/                             
      DATA D6/  2.81154470452317801923D-04/                             
      DATA D7/ -1.99228843754342707901D-05/                             
      DATA D8/ -4.61211975292357408243D-06/                             
      DATA D9/  2.71236535378185466378D-07/                             
      DATA D10/  4.76172600846326956261D-08/                            
      DATA D11/ -2.41413440645675455117D-09/                            
      DATA D12/ -3.37830909030876090705D-10/                            
      DATA D13/  1.48765237790396156913D-11/                            
      DATA D14/  1.72075812141360072788D-12/                            
      DATA E0/ -1.43371195263276446603D+04/                             
      DATA E1/ -2.24429804067515444935D+04/                             
      DATA E2/ -9.47372881573607965583D+03/                             
      DATA E3/ -1.22732832216812921247D+03/                             
      DATA E4/ -4.82584462361262805636D+01/                             
      DATA E5/ -8.25341486258634771351D-01/                             
      DATA F0/ -1.80165557510198819471D+05/                             
      DATA F1/ -2.76748521664359452706D+05/                             
      DATA F2/ -1.10834274415964139032D+05/                             
      DATA F3/ -1.20128613388289117320D+04/                             
      DATA F4/ -1.92416981016247577003D+02/                             
      DATA G0/  5.60975928436673365537D+03/                             
      DATA G1/  1.00684419438033930486D+04/                             
      DATA G2/  4.98653643523871678567D+03/                             
      DATA G3/  7.54262779548808089577D+02/                             
      DATA G4/  2.65889279842112002528D+01/                             
      DATA G5/  6.72742178111518345382D-02/                             
      DATA H0/  2.39349729466313955879D+04/                             
      DATA H1/  4.30288076179492217865D+04/                             
      DATA H2/  2.14064693316191034542D+04/                             
      DATA H3/  3.28862261076042955090D+03/                             
      DATA H4/  1.26450191906754194015D+02/                             
      IF(X.LT.0.D0) GO TO 60                                            
      IF(X.GE.8.D0) GO TO 50                                            
      ISWICH=X+1.0D0                                                    
      GO TO (10,10,20,20,30,30,40,40),ISWICH                            
   10 T=X*X                                                             
      DJ2=((((((((A8*T+A7)*T+A6)*T+A5)*T+A4)*T+A3)*                     
     *T+A2)*T+A1)*T+A0)*T                                               
      RETURN                                                            
   20 T=X-3.D0                                                          
      DJ2=(((((((((((((B14*T+B13)*T+B12)*T+B11)*T+B10)*T+B9)*           
     *T+B8)*T+B7)*T+B6)*T+B5)*T+B4)*T+B3)*T+B2)*T+B1)*T+B0              
      RETURN                                                            
   30 T=X-5.D0                                                          
      DJ2=(((((((((((((C14*T+C13)*T+C12)*T+C11)*T+C10)*T+C9)*           
     *T+C8)*T+C7)*T+C6)*T+C5)*T+C4)*T+C3)*T+C2)*T+C1)*T+C0              
      RETURN                                                            
   40 T=X-7.D0                                                          
      DJ2=(((((((((((((D14*T+D13)*T+D12)*T+D11)*T+D10)*T+D9)*           
     *T+D8)*T+D7)*T+D6)*T+D5)*T+D4)*T+D3)*T+D2)*T+D1)*T+D0              
      RETURN                                                            
   50 IF(X.GT.BIG) GO TO 60                                             
      Y=8.D0/X                                                          
      T=Y*Y                                                             
      H2A=(((((E5*T+E4)*T+E3)*T+E2)*T+E1)*T+E0)/                        
     *(((((T+F4)*T+F3)*T+F2)*T+F1)*T+F0)                                
      H2P=(((((G5*T+G4)*T+G3)*T+G2)*T+G1)*T+G0)/                        
     *(((((T+H4)*T+H3)*T+H2)*T+H1)*T+H0)*Y                              
      DJ2=-DCOS(H2P-QP+X)*DSQRT(Y*H2A)                                  
      RETURN                                                            
   60 DJ2=0.D0                                                          
cc      CALL MGDD(5HDJ2  ,X,DJ2,22HARG LT 0 OR GT 3.53D15,6)              
      RETURN                                                            
      END                                                               
C#NUMPAC#DJ3                 REVISED ON 1988-12-17                      
      FUNCTION DJ3(X)                                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DATA BIG / 3.53D15/                                               
      DATA QP  / 7.85398163397448309616D-01/                            
      DATA A0/  2.08333333333333333273D-02/                             
      DATA A1/ -1.30208333333333309671D-03/                             
      DATA A2/  3.25520833333317878475D-05/                             
      DATA A3/ -4.52112268514616640789D-07/                             
      DATA A4/  4.03671667825103238805D-09/                             
      DATA A5/ -2.52294757333466516750D-11/                             
      DATA A6/  1.16801656880176194754D-13/                             
      DATA A7/ -4.16788266907527181938D-16/                             
      DATA A8/  1.13637355722007137148D-18/                             
      DATA B0/  3.09062722255251643697D-01/                             
      DATA B1/  1.77028538330639329726D-01/                             
      DATA B2/ -2.95047563884399145099D-02/                             
      DATA B3/ -2.77836899420374464211D-02/                             
      DATA B4/ -2.73044379023257999864D-05/                             
      DATA B5/  1.40192656808142988672D-03/                             
      DATA B6/  2.73469642598496181657D-05/                             
      DATA B7/ -3.31093822275883150665D-05/                             
      DATA B8/ -6.54532300577873528847D-07/                             
      DATA B9/  4.50590585247672012623D-07/                             
      DATA B10/  7.87418322258387786173D-09/                            
      DATA B11/ -3.98884603526475783967D-09/                            
      DATA B12/ -6.03609575457832249117D-11/                            
      DATA B13/  2.44150142485683034319D-11/                            
      DATA B14/  3.21055558571731373531D-13/                            
      DATA C0/  3.64831230613666994050D-01/                             
      DATA C1/ -1.72333622090447959139D-01/                             
      DATA C2/ -9.95126315873285930712D-02/                             
      DATA C3/  1.51115881134712075560D-02/                             
      DATA C4/  7.38462837842238765168D-03/                             
      DATA C5/ -4.88355475695703152094D-04/                             
      DATA C6/ -2.37375167491371397076D-04/                             
      DATA C7/  8.95095984342474734000D-06/                             
      DATA C8/  4.12124879108359320770D-06/                             
      DATA C9/ -1.06919078875331794237D-07/                             
      DATA C10/ -4.44626919066586516635D-08/                            
      DATA C11/  8.89794765079101473739D-10/                            
      DATA C12/  3.26539411112591727469D-10/                            
      DATA C13/ -5.31018327965503557215D-12/                            
      DATA C14/ -1.70938692121214854372D-12/                            
      DATA D0/ -1.67555587995334235783D-01/                             
      DATA D1/ -2.29607682373653933018D-01/                             
      DATA D2/  8.47905846574382125254D-02/                             
      DATA D3/  2.78860158287113686069D-02/                             
      DATA D4/ -5.57271130223726681037D-03/                             
      DATA D5/ -1.25675904404574015195D-03/                             
      DATA D6/  1.58337500371148022842D-04/                             
      DATA D7/  2.85079469162923786821D-05/                             
      DATA D8/ -2.56044984050468414679D-06/                             
      DATA D9/ -3.81630936514958642808D-07/                             
      DATA D10/  2.66587464593458891363D-08/                            
      DATA D11/  3.35003902135922773719D-09/                            
      DATA D12/ -1.92736116492975254598D-10/                            
      DATA D13/ -2.03942728287784599590D-11/                            
      DATA D14/  1.00382736709225010974D-12/                            
      DATA E0/  1.69364906519095270055D+04/                             
      DATA E1/  2.14846111816965821861D+04/                             
      DATA E2/  6.89281029747769803640D+03/                             
      DATA E3/  7.03346665574037790191D+02/                             
      DATA E4/  4.40023290425369719677D+01/                             
      DATA E5/  1.88577604483265604622D+00/                             
      DATA F0/  2.12830218438524709222D+05/                             
      DATA F1/  2.55434645901034435429D+05/                             
      DATA F2/  6.80054120822196897219D+04/                             
      DATA F3/  2.76727851405537195297D+03/                             
      DATA F4/ -5.15745113868748413656D+01/                             
      DATA G0/  1.31038735399658407632D+04/                             
      DATA G1/  1.71763145681739238313D+04/                             
      DATA G2/  5.41037758027535673915D+03/                             
      DATA G3/  4.50304091753467408014D+02/                             
      DATA G4/  1.83978344912001082916D+01/                             
      DATA G5/  3.04569588535784760451D-01/                             
      DATA H0/  2.39613687587946836985D+04/                             
      DATA H1/  3.13223188955297746519D+04/                             
      DATA H2/  9.79696395157882846807D+03/                             
      DATA H3/  8.09581658036249480692D+02/                             
      DATA H4/  3.77586303399260839521D+01/                             
      IF(X.LT.0.D0) GO TO 60                                            
      IF(X.GE.8.D0) GO TO 50                                            
      ISWICH=X+1.0D0                                                    
      GO TO (10,10,20,20,30,30,40,40),ISWICH                            
   10 T=X*X                                                             
      DJ3=((((((((A8*T+A7)*T+A6)*T+A5)*T+A4)*T+A3)*                     
     *T+A2)*T+A1)*T+A0)*T*X                                             
      RETURN                                                            
   20 T=X-3.D0                                                          
      DJ3=(((((((((((((B14*T+B13)*T+B12)*T+B11)*T+B10)*T+B9)*           
     *T+B8)*T+B7)*T+B6)*T+B5)*T+B4)*T+B3)*T+B2)*T+B1)*T+B0              
      RETURN                                                            
   30 T=X-5.D0                                                          
      DJ3=(((((((((((((C14*T+C13)*T+C12)*T+C11)*T+C10)*T+C9)*           
     *T+C8)*T+C7)*T+C6)*T+C5)*T+C4)*T+C3)*T+C2)*T+C1)*T+C0              
      RETURN                                                            
   40 T=X-7.D0                                                          
      DJ3=(((((((((((((D14*T+D13)*T+D12)*T+D11)*T+D10)*T+D9)*           
     *T+D8)*T+D7)*T+D6)*T+D5)*T+D4)*T+D3)*T+D2)*T+D1)*T+D0              
      RETURN                                                            
   50 IF(X.GT.BIG) GO TO 60                                             
      Y=8.D0/X                                                          
      T=Y*Y                                                             
      H3A=(((((E5*T+E4)*T+E3)*T+E2)*T+E1)*T+E0)/                        
     *(((((T+F4)*T+F3)*T+F2)*T+F1)*T+F0)                                
      H3P=(((((G5*T+G4)*T+G3)*T+G2)*T+G1)*T+G0)/                        
     *(((((T+H4)*T+H3)*T+H2)*T+H1)*T+H0)*Y                              
      DJ3=DCOS(H3P+QP+X)*DSQRT(Y*H3A)                                   
      RETURN                                                            
   60 DJ3=0.D0                                                          
cc      CALL MGDD(5HDJ3  ,X,DJ3,22HARG LT 0 OR GT 3.53D15,6)              
      RETURN                                                            
      END                                                               
C#NUMPAC#DJ4                 REVISED ON 1989-01-27                      
      FUNCTION DJ4(X)                                                   
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DATA BIG / 3.53D15/                                               
      DATA QP  / 7.85398163397448309616D-01/                            
      DATA A0/  2.60416666666666666644D-03/                             
      DATA A1/ -1.30208333333333324394D-04/                             
      DATA A2/  2.71267361111105247379D-06/                             
      DATA A3/ -3.22937334654599034251D-08/                             
      DATA A4/  2.52294792511047853242D-10/                             
      DATA A5/ -1.40163760184641997114D-12/                             
      DATA A6/  5.84010017006001074402D-15/                             
      DATA A7/ -1.89474953667861887900D-17/                             
      DATA A8/  4.75050304525352146350D-20/                             
      DATA B0/  1.32034183924612210664D-01/                             
      DATA B1/  1.33017143689102041788D-01/                             
      DATA B2/  2.91771031336099274643D-02/                             
      DATA B3/ -9.61648329292712520211D-03/                             
      DATA B4/ -3.79503693033903419368D-03/                             
      DATA B5/  1.66312812469513689092D-04/                             
      DATA B6/  1.48981108566116224939D-04/                             
      DATA B7/ -2.78267871441911678704D-07/                             
      DATA B8/ -2.87167650958189611698D-06/                             
      DATA B9/ -2.26984977689211333829D-08/                             
      DATA B10/  3.30773663331192357416D-08/                            
      DATA B11/  3.37196454506533603958D-10/                            
      DATA B12/ -2.54217530508694296813D-10/                            
      DATA B13/ -2.57583872004111347381D-12/                            
      DATA B14/  1.37582936041668812839D-12/                            
      DATA C0/  3.91232360458648177730D-01/                             
      DATA C1/  5.18453422467483695895D-02/                             
      DATA C2/ -7.56063591072315065340D-02/                             
      DATA C3/ -1.44172416922437194145D-02/                             
      DATA C4/  3.84830004128775642846D-03/                             
      DATA C5/  8.49933186978045035355D-04/                             
      DATA C6/ -9.02182561052840360704D-05/                             
      DATA C7/ -2.20871070447176924182D-05/                             
      DATA C8/  1.25844527403471841424D-06/                             
      DATA C9/  3.21240958047904476544D-07/                             
      DATA C10/ -1.17896818200417424025D-08/                            
      DATA C11/ -2.98453001094431210129D-09/                            
      DATA C12/  7.93638279527648407000D-11/                            
      DATA C13/  1.89343199071609859533D-11/                            
      DATA C14/ -3.94751016537232077821D-13/                            
      DATA D0/  1.57798144661367917688D-01/                             
      DATA D1/ -2.57725956373258705087D-01/                             
      DATA D2/ -3.47271130531870090315D-02/                             
      DATA D3/  2.72518613594043983193D-02/                             
      DATA D4/  3.18677378366351532443D-03/                             
      DATA D5/ -1.04620438729938481735D-03/                             
      DATA D6/ -1.17956788683907078060D-04/                             
      DATA D7/  2.10443130651624511999D-05/                             
      DATA D8/  2.25724436047429505740D-06/                             
      DATA D9/ -2.61938380073812647261D-07/                             
      DATA D10/ -2.60952108618968802432D-08/                            
      DATA D11/  2.21152038439633210401D-09/                            
      DATA D12/  2.01571190425871738818D-10/                            
      DATA D13/ -1.32266357175402394279D-11/                            
      DATA D14/ -1.09540620620791083727D-12/                            
      DATA E0/  3.60438505629335953856D+04/                             
      DATA E1/  3.82600583078652390218D+04/                             
      DATA E2/  8.96939186557248718185D+03/                             
      DATA E3/  1.00716298901914860942D+03/                             
      DATA E4/  1.63463922745439384465D+02/                             
      DATA E5/  1.03482342314919507273D+01/                             
      DATA F0/  4.52940384542402072608D+05/                             
      DATA F1/  4.25057173544383415249D+05/                             
      DATA F2/  5.14303467221580674206D+04/                             
      DATA F3/ -3.23961115701570965134D+03/                             
      DATA F4/  3.04711396414673849956D+02/                             
      DATA F5/ -2.40375626100426314713D+01/                             
      DATA G0/ -1.36987974486182787316D+04/                             
      DATA G1/ -1.65012863730579455202D+04/                             
      DATA G2/ -3.65792471072781434564D+03/                             
      DATA G3/ -1.47630098689480609176D+02/                             
      DATA G4/ -4.24886555853072408691D+01/                             
      DATA G5/  1.24503683257393075346D+00/                             
      DATA H0/ -1.39162386779614265338D+04/                             
      DATA H1/ -1.65865405547413437955D+04/                             
      DATA H2/ -3.51670713483149870412D+03/                             
      DATA H3/ -1.22213659982884073844D+02/                             
      DATA H4/ -4.85427921152916547231D+01/                             
      IF(X.LT.0.D0) GO TO 60                                            
      IF(X.GE.8.D0) GO TO 50                                            
      ISWICH=X+1.0D0                                                    
      GO TO (10,10,20,20,30,30,40,40),ISWICH                            
   10 T=X*X                                                             
      DJ4=((((((((A8*T+A7)*T+A6)*T+A5)*T+A4)*T+A3)*                     
     *T+A2)*T+A1)*T+A0)*T*T                                             
      RETURN                                                            
   20 T=X-3.D0                                                          
      DJ4=(((((((((((((B14*T+B13)*T+B12)*T+B11)*T+B10)*T+B9)*           
     *T+B8)*T+B7)*T+B6)*T+B5)*T+B4)*T+B3)*T+B2)*T+B1)*T+B0              
      RETURN                                                            
   30 T=X-5.D0                                                          
      DJ4=(((((((((((((C14*T+C13)*T+C12)*T+C11)*T+C10)*T+C9)*           
     *T+C8)*T+C7)*T+C6)*T+C5)*T+C4)*T+C3)*T+C2)*T+C1)*T+C0              
      RETURN                                                            
   40 T=X-7.D0                                                          
      DJ4=(((((((((((((D14*T+D13)*T+D12)*T+D11)*T+D10)*T+D9)*           
     *T+D8)*T+D7)*T+D6)*T+D5)*T+D4)*T+D3)*T+D2)*T+D1)*T+D0              
      RETURN                                                            
   50 IF(X.GT.BIG) GO TO 60                                             
      Y=8.D0/X                                                          
      T=Y*Y                                                             
      H4A=(((((E5*T+E4)*T+E3)*T+E2)*T+E1)*T+E0)/                        
     *((((((T+F5)*T+F4)*T+F3)*T+F2)*T+F1)*T+F0)                         
      H4P=(((((G5*T+G4)*T+G3)*T+G2)*T+G1)*T+G0)/                        
     *(((((T+H4)*T+H3)*T+H2)*T+H1)*T+H0)*Y                              
      DJ4=DCOS(H4P-QP+X)*DSQRT(Y*H4A)                                   
      RETURN                                                            
   60 DJ4=0.D0                                                          
cc      CALL MGDD(5HDJ4  ,X,DJ4,22HARG LT 0 OR GT 3.53D15,6)              
      RETURN                                                            
      END                                                               
C#NUMPAC#DFCTRL              REVISED ON 1984-11-30                      
      FUNCTION DFCTRL(N)                                                
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION FACT(57)                                                
      DATA (FACT(J),J=1,30)/                                            
     *  0.10000000000000000D+01,  0.10000000000000000D+01,              
     *  0.20000000000000000D+01,  0.60000000000000000D+01,              
     *  0.24000000000000000D+02,  0.12000000000000000D+03,              
     *  0.72000000000000000D+03,  0.50400000000000000D+04,              
     *  0.40320000000000000D+05,  0.36288000000000000D+06,              
     *  0.36288000000000000D+07,  0.39916800000000000D+08,              
     *  0.47900160000000000D+09,  0.62270208000000000D+10,              
     *  0.87178291200000000D+11,  0.13076743680000000D+13,              
     *  0.20922789888000000D+14,  0.35568742809600000D+15,              
     *  0.64023737057280000D+16,  0.12164510040883200D+18,              
     *  0.24329020081766400D+19,  0.51090942171709440D+20,              
     *  0.11240007277776077D+22,  0.25852016738884977D+23,              
     *  0.62044840173323944D+24,  0.15511210043330986D+26,              
     *  0.40329146112660564D+27,  0.10888869450418352D+29,              
     *  0.30488834461171386D+30,  0.88417619937397020D+31/              
      DATA (FACT(J),J=31,57)/                                           
     *  0.26525285981219106D+33,  0.82228386541779228D+34,              
     *  0.26313083693369353D+36,  0.86833176188118865D+37,              
     *  0.29523279903960414D+39,  0.10333147966386145D+41,              
     *  0.37199332678990122D+42,  0.13763753091226345D+44,              
     *  0.52302261746660111D+45,  0.20397882081197443D+47,              
     *  0.81591528324789773D+48,  0.33452526613163807D+50,              
     *  0.14050061177528799D+52,  0.60415263063373836D+53,              
     *  0.26582715747884488D+55,  0.11962222086548019D+57,              
     *  0.55026221598120890D+58,  0.25862324151116818D+60,              
     *  0.12413915592536073D+62,  0.60828186403426756D+63,              
     *  0.30414093201713378D+65,  0.15511187532873823D+67,              
     *  0.80658175170943879D+68,  0.42748832840600256D+70,              
     *  0.23084369733924138D+72,  0.12696403353658276D+74,              
     *  0.71099858780486345D+75/                                        
      IF(N.LT.0.OR.N.GT.56) GO TO 10                                    
      DFCTRL=FACT(N+1)                                                  
      RETURN                                                            
   10 DFCTRL=0.                                                         
cc      CALL MGID(6HDFCTRL,N,DFCTRL,17HARG LT 0 OR GT 56,5)               
      RETURN                                                            
      END                                                               
      FUNCTION DLGAMA(X)                                                
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DATA P  / 0.91893853320467274178D+00/                             
      DATA A0 / 0.67945554609656524353D+04/                             
      DATA A1 / 0.43009657105424900568D+04/                             
      DATA A2 / 0.20187899154002913727D+04/                             
      DATA A3 / 0.59965748999418715670D+03/                             
      DATA A4 / 0.15269750682610185725D+03/                             
      DATA A5 / 0.27474380854875363526D+02/                             
      DATA A6 / 0.45626974981667575657D+01/                             
      DATA A7 / 0.45637593637129094235D+00/                             
      DATA A8 / 0.50405460836766191827D-01/                             
      DATA B0 / 0.67945554609656524350D+04/                             
      DATA B1 / 0.14283340976884678912D+04/                             
      DATA B2 /-0.13833593325342669889D+04/                             
      DATA B3 / 0.41995656691668842732D+02/                             
      DATA B4 / 0.83657448811029595617D+02/                             
      DATA B5 /-0.16578359883065247032D+02/                             
      DATA C0 / 0.28606295245643387431D+06/                             
      DATA C1 / 0.49887307361909404998D+03/                             
      DATA C2 / 0.69389051213390013235D-01/                             
      DATA D0 / 0.34327554294772064918D+07/                             
      DATA D1 / 0.60982202242324210839D+04/                             
      IF(X.LE.0.D0) GO TO 40                                            
      IF(X.GT.32.D0) GO TO 30                                           
      N=X                                                               
      F=X-DFLOAT(N)                                                     
      Y=((((((((A8*F+A7)*F+A6)*F+A5)*F+A4)*F+A3)*F+A2)*F+A1)*F+A0)      
     */((((((F+B5)*F+B4)*F+B3)*F+B2)*F+B1)*F+B0)                        
      IF(N.GE.3) GO TO 10                                               
      IF(N.EQ.1) Y=Y/X                                                  
      IF(N.EQ.0) Y=Y/((X+1.D0)*X)                                       
      DLGAMA=DLOG(Y)                                                    
      RETURN                                                            
   10 V=X                                                               
      DO 20 J=3,N                                                       
      V=V-1.D0                                                          
   20 Y=Y*V                                                             
      DLGAMA=DLOG(Y)                                                    
      RETURN                                                            
   30 T=1024.D0/(X*X)                                                   
      DLGAMA=((C2*T+C1)*T+C0)/((T+D1)*T+D0)/X+P                         
     *-X+DLOG(X)*(X-0.5D0)                                              
      RETURN                                                            
   40 DLGAMA=0.D0                                                       
      RETURN                                                            
      END                                                               
C#NUMPAC#FLPOWD              REVISED ON 1984-11-30                      
      SUBROUTINE FLPOWD(X,N,B,KB,FUNC,GRAD,LF,NF,FLB,EPS,FM,ILL)        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION X(N),B(KB,N),G(1000),G1(1000),G2(1000)                  
     *,X1(1000),X2(1000),S(1000),Y(1000)                                
      ILL=30000                                                         
      IF(N.LT.2.OR.N.GT.MIN0(1000,KB,LF).OR.EPS.LE.0.D0) RETURN         
      ILL=0                                                             
      NF=1                                                              
      F=FUNC(X)                                                         
      CALL GRAD(X,G)                                                    
      DO 20 J=1,N                                                       
      DO 10 I=1,N                                                       
   10 B(I,J)=0.D0                                                       
   20 B(J,J)=1.D0                                                       
   30 IF(NF.GT.LF) GO TO 240                                            
      T=1.D0                                                            
      Q=0.D0                                                            
      DO 50 I=1,N                                                       
      W=0.D0                                                            
      DO 40 J=1,N                                                       
   40 W=B(J,I)*G(J)+W                                                   
      S(I)=-W                                                           
   50 Q=S(I)*G(I)+Q                                                     
      U=(FLB-F)*2.D0/Q                                                  
      FL=DSIGN(DMIN1(DABS(U),1.D0),U)                                   
   60 FL=FL*T                                                           
      DO 70 I=1,N                                                       
   70 X1(I)=S(I)*FL+X(I)                                                
      NF=NF+1                                                           
      F1=FUNC(X1)                                                       
      CALL GRAD(X1,G1)                                                  
      Q1=0.D0                                                           
      DO 80 I=1,N                                                       
   80 Q1=S(I)*G1(I)+Q1                                                  
      IF(FL*Q1.LT.0.D0.AND.F1.LT.F) GO TO 170                           
C     INTERPOLATION                                                     
      Z=(F-F1)*3.D0/FL+Q+Q1                                             
      W=Z*Z-Q*Q1                                                        
      IF(W.LT.0.D0) GO TO 230                                           
      W=DSQRT(W)                                                        
      A=(Z+W-Q)*FL/(Q1+W+W-Q)                                           
      DO 90 I=1,N                                                       
   90 X2(I)=S(I)*A+X(I)                                                 
      NF=NF+1                                                           
      F2=FUNC(X2)                                                       
      IF(F2.GT.DMIN1(F,F1)) GO TO 230                                   
      CALL GRAD(X2,G2)                                                  
      DO 100 I=1,N                                                      
  100 Y(I)=G2(I)-G(I)                                                   
C     UPDATE B                                                          
      SY=0.D0                                                           
      YBY=0.D0                                                          
      DO 120 I=1,N                                                      
      SY=S(I)*Y(I)+SY                                                   
      G1(I)=0.D0                                                        
      DO 110 J=1,N                                                      
  110 G1(I)=B(J,I)*Y(J)+G1(I)                                           
  120 YBY=Y(I)*G1(I)+YBY                                                
      IF(SY*YBY.EQ.0.D0) GO TO 140                                      
      DO 130 J=1,N                                                      
      W=S(J)*A/SY                                                       
      T=G1(J)/YBY                                                       
      DO 130 I=1,J                                                      
      B(I,J)=S(I)*W-G1(I)*T+B(I,J)                                      
  130 B(J,I)=B(I,J)                                                     
  140 DO 150 I=1,N                                                      
      X(I)=X2(I)                                                        
  150 G(I)=G2(I)                                                        
      F=F2                                                              
C     TEST OF CONVERGENCE                                               
      IF(DABS(F-F1).GT.DMAX1(DABS(F)*EPS,1.D-17)) GO TO 30              
      DO 160 I=1,N                                                      
      IF(DABS(A*S(I)).GT.DMAX1(DABS(X(I))*EPS,1.D-17)) GO TO 30         
  160 CONTINUE                                                          
      FM=F                                                              
      RETURN                                                            
C     EXTRAPOLATION                                                     
  170 DO 180 J=1,N                                                      
      W=-S(J)/Q                                                         
      DO 180 I=1,J                                                      
      B(I,J)=S(I)*W+B(I,J)                                              
  180 B(J,I)=B(I,J)                                                     
  190 DO 200 I=1,N                                                      
      X(I)=X1(I)                                                        
  200 G(I)=G1(I)                                                        
      F=F1                                                              
      IF(T.LT.0.97D0) GO TO 30                                          
      FL=FL+FL                                                          
      Q=0.D0                                                            
      DO 220 I=1,N                                                      
      W=0.D0                                                            
      DO 210 J=1,N                                                      
  210 W=B(J,I)*G(J)+W                                                   
      S(I)=-W                                                           
  220 Q=S(I)*G(I)+Q                                                     
      IF(FL*Q*F.GT.0.D0) FL=(FLB-F)*4.D0/Q                              
      GO TO 60                                                          
C     STEP REDUCTION                                                    
  230 IF(F1.LT.F) GO TO 190                                             
      IF(T.LT.1.D-17) GO TO 240                                         
      T=T*T*0.954D0                                                     
      GO TO 60                                                          
  240 ILL=1                                                             
      RETURN                                                            
      END                                                               
C# ORIGINAL#NUMPAC#FLPOWD              REVISED ON 1984-11-30
C# MODIFIED BY N. YOSHIDA
      SUBROUTINE FLPOWDX(X,N,B,KB,GRAD,LF,NF,FLB,EPS,FM,ILL)        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION X(N),B(KB,N),G(1000),G1(1000),G2(1000)                  
     *,X1(1000),X2(1000),S(1000),Y(1000)                                
      ILL=30000                                                         
      IF(N.LT.2.OR.N.GT.MIN0(1000,KB,LF).OR.EPS.LE.0.D0) RETURN         
      ILL=0                                                             
      NF=1                                                              
      CALL GRAD(X,F,G)                                                 
      DO 20 J=1,N                                                       
      DO 10 I=1,N                                                       
   10 B(I,J)=0.D0                                                       
   20 B(J,J)=1.D0                                                       
   30 IF(NF.GT.LF) GO TO 240                                            
      T=1.D0                                                            
      Q=0.D0                                                            
      DO 50 I=1,N                                                       
      W=0.D0                                                            
      DO 40 J=1,N                                                       
   40 W=B(J,I)*G(J)+W                                                   
      S(I)=-W                                        
   50 Q=S(I)*G(I)+Q                                                     
      U=(FLB-F)*2.D0/Q                                                  
      FL=DSIGN(DMIN1(DABS(U),1.D0),U)                      
   60 FL=FL*T                                              
      DO 70 I=1,N                                                       
   70 X1(I)=S(I)*FL+X(I)                                                
      NF=NF+1                                                           
      CALL GRAD(X1,F1,G1)                                            
      Q1=0.D0                                                           
      DO 80 I=1,N                                                       
   80 Q1=S(I)*G1(I)+Q1                                                  
      IF(FL*Q1.LT.0.D0.AND.F1.LT.F) GO TO 170                           
C     INTERPOLATION                                                     
      Z=(F-F1)*3.D0/FL+Q+Q1                                             
      W=Z*Z-Q*Q1                                                        
      IF(W.LT.0.D0) GO TO 230                                           
      W=DSQRT(W)                                                        
      A=(Z+W-Q)*FL/(Q1+W+W-Q)                                           
      DO 90 I=1,N                                                       
   90 X2(I)=S(I)*A+X(I)                                                 
      NF=NF+1                                                           
      CALL GRAD(X2,F2,G2)                                     
      IF(F2.GT.DMIN1(F,F1)) GO TO 230                                   
      DO 100 I=1,N                                                      
  100 Y(I)=G2(I)-G(I)                                                   
C     UPDATE B                                                          
      SY=0.D0                                                           
      YBY=0.D0                                                          
      DO 120 I=1,N                                                      
      SY=S(I)*Y(I)+SY                                                   
      G1(I)=0.D0                                                        
      DO 110 J=1,N                                                      
  110 G1(I)=B(J,I)*Y(J)+G1(I)                                           
  120 YBY=Y(I)*G1(I)+YBY                                                
      IF(SY*YBY.EQ.0.D0) GO TO 140                                      
      DO 130 J=1,N                                                      
      W=S(J)*A/SY                                                       
      T=G1(J)/YBY                                                       
      DO 130 I=1,J                                                      
      B(I,J)=S(I)*W-G1(I)*T+B(I,J)                                      
  130 B(J,I)=B(I,J)                                                     
  140 DO 150 I=1,N                                                      
      X(I)=X2(I)                                                        
  150 G(I)=G2(I)                                                        
      F=F2                                                              
C     TEST OF CONVERGENCE                                               
      IF(DABS(F-F1).GT.DMAX1(DABS(F)*EPS,1.D-17)) GO TO 30              
      DO 160 I=1,N                                                      
      IF(DABS(A*S(I)).GT.DMAX1(DABS(X(I))*EPS,1.D-17)) GO TO 30         
  160 CONTINUE                                                          
      FM=F                                                              
      RETURN                                                            
C     EXTRAPOLATION                                                     
  170 DO 180 J=1,N                                                      
      W=-S(J)/Q                                                         
      DO 180 I=1,J                                                      
      B(I,J)=S(I)*W+B(I,J)                                              
  180 B(J,I)=B(I,J)                                                     
  190 DO 200 I=1,N                                                      
      X(I)=X1(I)                                                        
  200 G(I)=G1(I)                                                        
      F=F1                                                              
      IF(T.LT.0.97D0) GO TO 30                                          
      FL=FL+FL                                                          
      Q=0.D0                                                            
      DO 220 I=1,N                                                      
      W=0.D0                                                            
      DO 210 J=1,N                                                      
  210 W=B(J,I)*G(J)+W                                                   
      S(I)=-W                                                           
  220 Q=S(I)*G(I)+Q                                                     
      IF(FL*Q*F.GT.0.D0) FL=(FLB-F)*4.D0/Q                              
      GO TO 60                                                          
C     STEP REDUCTION                                                    
  230 IF(F1.LT.F) GO TO 190                                             
      IF(T.LT.1.D-17) GO TO 240                                         
      T=T*T*0.954D0                                                     
      GO TO 60                                                          
  240 ILL=1                                                             
      RETURN                                                            
      END                                                               
C#NUMPAC#LEQLUD              REVISED ON 1988-06-06                      
      SUBROUTINE LEQLUD(A,KA,N,X,KX,M,D,MAX,EPS,IND)                    
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(KA,N),X(KX,M),MAX(N)                                  
      IF(N.LT.2.OR.N.GT.KA.OR.N.GT.KX.OR.M.LT.0.OR.EPS.LE.0.) GO TO 170 
      IF(D.NE.0.D0) D=1.D0                                              
      IF(IND.NE.0) GO TO 90                                             
      DO 80 J=1,N                                                       
      DO 20 I=1,J-1                                                     
      L=MAX(I)                                                          
      S=A(L,J)                                                          
      A(L,J)=A(I,J)                                                     
      A(I,J)=S                                                          
      DO 10 K=1,I-1                                                     
   10 A(I,J)=A(I,K)*A(K,J)+A(I,J)                                       
   20 CONTINUE                                                          
      DO 40 I=J,N                                                       
      DO 30 K=1,J-1                                                     
   30 A(I,J)=A(I,K)*A(K,J)+A(I,J)                                       
   40 CONTINUE                                                          
      AM=0.D0                                                           
      DO 50 I=J,N                                                       
      AA=DABS(A(I,J))                                                   
      IF(AA.LE.AM) GO TO 50                                             
      AM=AA                                                             
      L=I                                                               
   50 CONTINUE                                                          
      IF(AM.LT.EPS) GO TO 160                                           
      MAX(J)=L                                                          
      IF(L.EQ.J) GO TO 70                                               
      DO 60 K=1,J                                                       
      W=A(L,K)                                                          
      A(L,K)=A(J,K)                                                     
   60 A(J,K)=W                                                          
      D=-D                                                              
   70 D=A(J,J)*D                                                        
      W=-A(J,J)                                                         
      DO 80 I=J+1,N                                                     
   80 A(I,J)=A(I,J)/W                                                   
   90 IF(M.EQ.0) GO TO 150                                              
      DO 140 J=1,M                                                      
      DO 110 I=1,N                                                      
      L=MAX(I)                                                          
      S=X(L,J)                                                          
      X(L,J)=X(I,J)                                                     
      X(I,J)=S                                                          
      DO 100 K=1,I-1                                                    
  100 X(I,J)=A(I,K)*X(K,J)+X(I,J)                                       
  110 CONTINUE                                                          
      DO 130 I=N,1,-1                                                   
      DO 120 K=I+1,N                                                    
  120 X(I,J)=X(I,J)-A(I,K)*X(K,J)                                       
  130 X(I,J)=X(I,J)/A(I,I)                                              
  140 CONTINUE                                                          
  150 IND=0                                                             
      RETURN                                                            
  160 IND=J                                                             
      RETURN                                                            
  170 IND=30000                                                         
      RETURN                                                            
      END                                                               
C#NUMPAC#MULMMW              REVISED ON 1987-08-03                      
C*    *** MULMMW *******************************************************
C     *                                                                *
C     *   MULMMW......NUMPAC S810-VERSION.                             *
C     *                                                                *
C     *   AUTHOR....I.NINOMIYA (CHUBU UNIVERSITY)  1987.07             *
C     *   CODER.....I.NINOMIYA (CHUBU UNIVERSITY)  1987.07             *
C     *                                                                *
C     *   'MULMMW' COMPUTES THE PRODUCT C=A*B OF MATRICES A AND B      *
C     *                                                                *
C     *   USAGE                                                        *
C     *        CALL MULMMW(A,B,C,KA,KB,KC,L,M,N,ILL)                   *
C     *          A   ....GIVEN MULTIPLICAND MATRIX OF SIZE L X M.      *
C     *          B   ....GIVEN MULTIPLIER MATRIX OF SIZE M X N.        *
C     *          C   ....RESULTANT PRODUCT MATRIX OF SIZE L X N        *
C     *          KA  ....GIVEN ADJUSTABLE DIMENSION OF A.              *
C     *          KB  ....GIVEN ADJUSTABLE DIMENSION OF B.              *
C     *          KC  ....GIVEN ADJUSTABLE DIMENSION OF C.              *
C     *          L   ....GIVEN NUMBER OF ROWS OF A AND C.              *
C     *          M   ....GIVEN NUMBER OF COLUMNS(ROWS) OF A(B).        *
C     *          N   ....GIVEN NUMBER OF COLUMNS OF B AND C.           *
C     *          ILL ....RESULTANT ERROR CODE.                         *
C     *             0          ....NORMAL TERMINATION.                 *
C     *             30000      ....PARAMETER ERROR.                    *
C     *                                                                *
C     *   SLAVE SUBROUTINE                                             *
C     *        NONE.                                                   *
C     *                                                                *
C     ******************************************************************
C                                                                       
C     ******************************************************************
C                                                                       
      SUBROUTINE MULMMW(A,B,C,KA,KB,KC,L,M,N,ILL)                       
C                                                                       
C     ******************************************************************
C     ----------------------------------------------------------------- 
C     DECLARATION                                                       
C     ----------------------------------------------------------------- 
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(KA,M),B(KB,N),C(KC,N)                                 
C     ----------------------------------------------------------------- 
C     PARAMETER ERROR CHECK                                             
C     ----------------------------------------------------------------- 
      IF(L.LT.1.OR.L.GT.KA.OR.L.GT.KC) GO TO 70                         
      IF(M.LT.1.OR.M.GT.KB) GO TO 70                                    
      IF(N.LT.1) GO TO 70                                               
C     ----------------------------------------------------------------- 
C     OUTER PRODUCT SCHEME (LOOP UNROLLING OF MULTIPLICITY 8)           
C     ----------------------------------------------------------------- 
      MM=M/8*8                                                          
      DO 60 J=1,N                                                       
      DO 10 I=1,L                                                       
   10 C(I,J)=0.D0                                                       
      DO 30 K=1,MM,8                                                    
      DO 20 I=1,L                                                       
   20 C(I,J)=C(I,J)+A(I,K)*B(K,J)+A(I,K+1)*B(K+1,J)+A(I,K+2)*B(K+2,J)   
     *+A(I,K+3)*B(K+3,J)+A(I,K+4)*B(K+4,J)+A(I,K+5)*B(K+5,J)            
     *+A(I,K+6)*B(K+6,J)+A(I,K+7)*B(K+7,J)                              
   30 CONTINUE                                                          
      DO 50 K=MM+1,M                                                    
      DO 40 I=1,L                                                       
   40 C(I,J)=C(I,J)+A(I,K)*B(K,J)                                       
   50 CONTINUE                                                          
   60 CONTINUE                                                          
C     ----------------------------------------------------------------- 
C     NORMAL EXIT                                                       
C     ----------------------------------------------------------------- 
      ILL=0                                                             
      RETURN                                                            
C     ----------------------------------------------------------------- 
C     ERROR EXIT                                                        
C     ----------------------------------------------------------------- 
   70 ILL=30000                                                         
      RETURN                                                            
      END                                                               
C*
C********* CFS1A  ** SCFS1A*1 ******************************
C*                                                         *
C*       CURVE FITTING BY SPLINES WITH FIXED KNOTS         *
C*                                                         *
C*         THE KNOTS MUST PREVIOUSLY BE GIVEN              *
C*                                                         *
C*    DIMENSION XR(KOSU),FR(KOSU),SIGMAR(KOSU),XI(N+1)     *
C*             ,CJ(N+K-1),DRESP(N+K-1),STATI(3)            *
C*             ,IHIST(2,25),PERCT(10),PERLIM(10)           *
C*             ,WORKC((N+K-1)*(K+1)+K)                     *
C*             ,IWORKC(KOSU+N+K-1)                         *
C*                                                         *
C*             (N=100 , K=10 , KOSU=2001)                  *
C*                                                         *
C*********        **          ******************************
      SUBROUTINE CFS1A(XR,FR,SIGMAR,XI,CJ,DRESP,STATI,IHIST
     1          ,PERCT,WORKC,IWORKC,N,K,KOSU,IWR,ICON)
      DIMENSION XR(2001),FR(2001),SIGMAR(2001)
     1         ,XI(101),CJ(109),DRESP(109),STATI(3)
     2         ,IHIST(2,25),PERCT(10),WORKC(1209)
     3         ,IWORKC(2109)
      DATA HISINT,KOHIS / 0.2 , 25 /
C--------- PARAMETER CHECK ---------------------------------
      ICON=0
      IF(K.LT.1) ICON=-1
      IF(N.LT.1) ICON=-2
      IF(ICON.LT.0) RETURN
      XLEF=XI(1)
      XRIG=XI(1)
      NP1=N+1
      DO 1100 IP1=1,NP1
      IF(XI(IP1).LT.XLEF) XLEF=XI(IP1)
      IF(XI(IP1).GT.XRIG) XRIG=XI(IP1)
 1100 CONTINUE
      IF(XR(1   ).LT.XLEF) ICON=-3
      IF(XR(KOSU).GT.XRIG) ICON=-4
      IF(ICON.LT.0) RETURN
      KOSU1=KOSU-1
      DO 1110 K1=1,KOSU1
      IF(XR(K1).GT.XR(K1+1)) ICON=-5
 1110 CONTINUE
      IF(ICON.LT.0) RETURN
C--------- END OF PARAMETER CHECK --------------------------
      NA=N
      KA=K
      NAKM1=NA+KA-1
      IW1=1
      IW2=IW1+NAKM1*KA
      IW3=IW2+NAKM1
C     IW4=IW3+KA
      IIW1=1
      IIW2=IIW1+KOSU
C     IIW3=IIW2+NA+KA-1
      CALL REOK1(XI,XR,IWORKC(IIW1),NA,KOSU)
      CALL SPCA1(XR,FR,SIGMAR,XI,CJ,DRESP,NA,KA,KOSU,IWR
     1    ,IWORKC(IIW1),IWORKC(IIW2),WORKC(IW1),WORKC(IW2)
     2    ,WORKC(IW3),NAKM1)
      CALL SPOB1(XR,FR,SIGMAR,XI,CJ,STATI,IHIST,XINS
     1    ,HISINT,KOHIS,PERCT,NA,KA,KOSU,IWR,IWORKC(IIW1)
     2    ,WORKC(IW2),WORKC(IW3))
      RETURN
      END
C*
C********* SPCA1  ** SSPCA1*1 ******************************
C*                                                         *
C*            CURVE FITTING BY SPLINES                     *
C*                                                         *
C*       CALCULATION OF MULTIPLIERS OF B-SPLINES           *
C*                                                         *
C*    DIMENSION PMAT(N+K-1,K),UVEC(N+K-1),XI(N+1)          *
C*             ,XR(KOSU),FR(KOSU),SIGMAR(KOSU),IXR(KOSU)   *
C*             ,BN(K),CJ(N+K-1),DRESP(N+K-1)               *
C*             ,INDEP(N+K-1)                               *
C*                                                         *
C*             (N=100 , K=10 , KOSU=2001 )                 *
C*                                                         *
C*********        **          ******************************
      SUBROUTINE SPCA1(XR,FR,SIGMAR,XI,CJ,DRESP,N,K,KOSU,IWR
     1          ,IXR,INDEP,PMAT,UVEC,BN,NK1D)
      DIMENSION PMAT(NK1D,10),UVEC(NK1D),XI(101),XR(2001)
     1         ,FR(2001),SIGMAR(2001),IXR(2001),BN(10)
     2         ,CJ(NK1D),DRESP(NK1D),INDEP(NK1D)
      DATA ZEROP / 1.0E-20 /
      NM1=N-1
      NKM1=N+K-1
      DO 1110 IAPK=1,NKM1
      DO 1100 IB1=1,K
      PMAT(IAPK,IB1)=0.0
 1100 CONTINUE
      UVEC(IAPK)=0.0
 1110 CONTINUE
      DO 1150 K1=1,KOSU
      IN=IXR(K1)
      IF(IN.LT.0) GO TO 1140
      IF(IN.GT.NM1) GO TO 1140
      K1S=1
      IF(IWR.NE.0) K1S=K1
      WR1=(1.0/SIGMAR(K1S))**2
      FR1=FR(K1)
      CALL BAS0B(XR(K1),IN,N,K,XI,BN)
      IAKL=IN+1
      IAKU=IN+K
      DO 1130 IAPK=IAKL,IAKU
      IAO=IAPK-IN
      WB1=WR1*BN(IAO)
      DO 1120 IBPK=IAKL,IAPK
      IB1=IBPK-IAPK+K
      IBO=IBPK-IN
      PMAT(IAPK,IB1)=PMAT(IAPK,IB1)+WB1*BN(IBO)
 1120 CONTINUE
      UVEC(IAPK)=UVEC(IAPK)+WB1*FR1
 1130 CONTINUE
 1140 CONTINUE
 1150 CONTINUE
      AMAX=0.0
      DO 1160 IAPK=1,NKM1
      IF(PMAT(IAPK,K).GT.AMAX) AMAX=PMAT(IAPK,K)
      IF(PMAT(IAPK,K).EQ.0.0) PMAT(IAPK,K)=1.0
 1160 CONTINUE
      PZERO=AMAX*ZEROP
      DO 1170 IAPK=1,NKM1
      DRESP(IAPK)=UVEC(IAPK)
 1170 CONTINUE
      CALL CDB2A(PMAT,UVEC,PZERO,INDEP,NKM1,K,0,ICON,NK1D)
      DO 1180 IAPK=1,NKM1
      CJ(IAPK)=UVEC(IAPK)
      DRESP(IAPK)=DRESP(IAPK)*CJ(IAPK)
 1180 CONTINUE
      RETURN
      END
C*
C********* SPOB1  ** SSPOB1*1 ******************************
C*                                                         *
C*             CURVE FITTING BY SPLINES                    *
C*                                                         *
C*     CALCULATION OF STATISTICS , HISTGRAMS  AND          *
C*      PREPARATION FOR THE NEW KNOTS INSERTION            *
C*                                                         *
C*    DIMENSION STATI(3),IHIST(2,KOHIS),XI(N+1),CJ(N+K-1)  *
C*             ,XR(KOSU),FR(KOSU),SIGMAR(KOSU),IXR(KOSU)   *
C*             ,SIGW(N),BN(K),PERCT(10)                    *
C*                                                         *
C*             (N=100 , K=10 , KOSU=2001 , KOHIS=26)       *
C*                                                         *
C*********        **          ******************************
      SUBROUTINE SPOB1(XR,FR,SIGMAR,XI,CJ,STATI,IHIST,XINS
     1          ,HISINT,KOHIS,PERCT,N,K,KOSU,IWR,IXR,SIGW
     2          ,BN)
      DIMENSION STATI(3),IHIST(2,26),XI(101),CJ(109)
     1         ,XR(2001),FR(2001),SIGMAR(2001),IXR(2001)
     2         ,SIGW(100),BN(10),PERCT(10),IPERC(10)
      NM1=N-1
      RESID=0.0
      DO 1100 IP1=1,N
      SIGW(IP1)=0.0
 1100 CONTINUE
      DO 1110 IC1=1,KOHIS
      IHIST(1,IC1)=0
      IHIST(2,IC1)=0
 1110 CONTINUE
      DO 1120 IC1=1,10
      IPERC(IC1)=0
 1120 CONTINUE
      DO 1170 K1=1,KOSU
      IN=IXR(K1)
      IF(IN.LT.0) GO TO 1160
      IF(IN.GT.NM1) GO TO 1160
      CALL BAS0B(XR(K1),IN,N,K,XI,BN)
      SUM1=0.0
      DO 1130 IA=1,K
      IA1=IN+IA
      SUM1=SUM1+CJ(IA1)*BN(IA)
 1130 CONTINUE
      FAPP=SUM1
      K1S=1
      IF(IWR.NE.0) K1S=K1
      EPSI=(FR(K1)-FAPP)/SIGMAR(K1S)
      ERVAR=EPSI**2
      RESID=RESID+ERVAR
      IP1=IN+1
      SIGW(IP1)=SIGW(IP1)+ERVAR
      IF(EPSI.GT.0.0) GO TO 1140
      INTD=IFIX(-EPSI/HISINT)+1
      IF(INTD.GT.KOHIS) INTD=KOHIS
      IHIST(1,INTD)=IHIST(1,INTD)+1
      GO TO 1150
 1140 INTD=IFIX( EPSI/HISINT)+1
      IF(INTD.GT.KOHIS) INTD=KOHIS
      IHIST(2,INTD)=IHIST(2,INTD)+1
 1150 INTC=IFIX(ABS(EPSI))+1
      IF(INTC.GT.10) INTC=10
      IPERC(INTC)=IPERC(INTC)+1
 1160 CONTINUE
 1170 CONTINUE
      STATI(1)=RESID
      STATI(2)=RESID/FLOAT(KOSU-(N+K-1))
      STATI(3)=FLOAT(KOSU)*ALOG(RESID)+FLOAT(2*(N+K-1))
      AMAX=0.0
      DO 1190 IP1=1,N
      IF(AMAX.GT.SIGW(IP1)) GO TO 1180
      AMAX=SIGW(IP1)
      IMAX=IP1
 1180 CONTINUE
 1190 CONTINUE
      XINS=(XI(IMAX)+XI(IMAX+1))*0.5
      DO 1200 IC1=2,10
      IPERC(IC1)=IPERC(IC1)+IPERC(IC1-1)
 1200 CONTINUE
      DO 1210 IC1=1,10
      PERCT(IC1)=FLOAT(IPERC(IC1))/FLOAT(KOSU)
 1210 CONTINUE
      RETURN
      END
C*
C********* BAS0B  ** SBAS0B*1 ******************************
C*                                                         *
C*      COMPUTATION OF THE BASIS FUNCTION OF SPLINES       *
C*    (1) SIMPLE KNOTS BUT MULTIPLICITY-K ON BOTH ENDS     *
C*                                                         *
C*    @@@@@           N(J)(XP) ONLY              @@@@@     *
C*                                                         *
C*    DIMENSION XI(N+1),BN(K)                              *
C*                                                         *
C*    N.........THE NUMBER OF THE KNOTS OF B-SPLINES       *
C*    K.........ORDER OF THE B-SPLINES                     *
C*    XI........KNOTS OF B-SPLINES                         *
C*                                                         *
C*       MUST SATISFIES  XI(I+1)<=XP<XI(I+2)               *
C*                                                         *
C*    BN(1,2,...,K) CONTAINES N(I-K+1,I-K+2,...,I)         *
C*                                                         *
C*********        **          ******************************
      SUBROUTINE BAS0B(XP,I,N,K,XI,BN)
      DIMENSION XI(101),BN(10)
      BN(K)=1.0/(XI(I+2)-XI(I+1))
      DO 1120 IS=2,K
      JRS=K-IS+1
      DO 1110 JR=JRS,K
      IR  =I-K+JR
      IRPS=IR+IS
      IF(IR  .LT.0) IR  =0
      IF(IRPS.GT.N) IRPS=N
      TR  =XI(IR  +1)
      TRPS=XI(IRPS+1)
      AN1=0.0
      IF(JR.NE.JRS) AN1=BN(JR)
      AN2=0.0
      IF(JR.NE.K) AN2=BN(JR+1)
      BN(JR)=(XP-TR)*AN1+(TRPS-XP)*AN2
      IF(IS.NE.K) BN(JR)=BN(JR)/(TRPS-TR)
 1100 CONTINUE
 1110 CONTINUE
 1120 CONTINUE
      RETURN
      END
C*
C********* REOK1  ** SREOK1*1 ******************************
C*                                                         *
C*           CURVE FITTING BY SPLINES                      *
C*                                                         *
C*       REORDER THE KNOTS AND SEARCH INTERVAL             *
C*                                                         *
C*    DIMENSION XI(N+1),XR(KOSU),IXR(KOSU)                 *
C*                                                         *
C*              (N=100 , KOSU=2001)                        *
C*                                                         *
C*********        **          ******************************
      SUBROUTINE REOK1(XI,XR,IXR,N,KOSU)
      DIMENSION XI(101),XR(2001),IXR(2001)
      DATA ZERO / 1.0E-05 /
      NM1=N-1
      NP1=N+1
      EPSX=(XI(NP1)-XI(1))*ZERO
C *** REORDER THE KNOTS
      ICOUN=0
      I1=0
 1100 ISK=0
      I1=I1+1
      DO 1120 I2D=I1,N
      I2=N+I1-I2D
      IF(XI(I2).LE.XI(I2+1)) GO TO 1110
      ISK=ISK+1
      WK1=XI(I2)
      XI(I2)=XI(I2+1)
      XI(I2+1)=WK1
 1110 CONTINUE
 1120 CONTINUE
      IF(ISK.EQ.0) GO TO 1140
      ICOUN=ICOUN+1
      IF(ICOUN.LE.NM1) GO TO 1100
      DO 1130 IP1=1,N
      IF(XI(IP1+1)-XI(IP1).LT.EPSX) XI(IP1+1)=XI(IP1)
 1130 CONTINUE
C *** SEARCH INTERVAL
 1140 DO 1240 K1=1,KOSU
      XP=XR(K1)
      IN=IXR(K1)
      IF(IN.LT.  0) GO TO 1150
      IF(IN.GT.NM1) GO TO 1150
      IF(XP.LT.XI(IN+1)) GO TO 1150
      IF(XP.GE.XI(IN+2)) GO TO 1150
      GO TO 1200
 1150 IF(XP.LT.XI(  1)) GO TO 1210
      IF(XP.GT.XI(NP1)) GO TO 1220
      IL=  1
      IU=NP1
 1160 IF((IU-IL).LE.1) GO TO 1180
      IM=(IL+IU)/2
      IF(XP.GE.XI(IM)) GO TO 1170
      IU=IM
      GO TO 1160
 1170 IL=IM
      GO TO 1160
 1180 IF(XI(IL+1)-XI(IL).GT.EPSX) GO TO 1190
      IL=IL+1
      GO TO 1180
 1190 IN=IL-1
 1200 IXR(K1)=IN
      GO TO 1230
 1210 IXR(K1)= -1
      GO TO 1230
 1220 IXR(K1)=NP1
 1230 CONTINUE
 1240 CONTINUE
      RETURN
      END
C*
C********* CDB2A  ** SCDB2A*1 ******************************
C*                                                         *
C*      SOLUTION OF SIMULTANEOUS LINEAR EQUATION           *
C*    CHOLESKI DECOMPOSITION  (BAND MATRIX - LOWER)        *
C*    (SET THE VALUE 'ZERO' TO DEPENDENT VARIABLES)        *
C*                                                         *
C*    KIND=0...SIMULTANEOUS LINEAR EQUATION                *
C*    KIND=1...CHOLESKI DECOMPOSITION                      *
C*    KIND=2...FOWARD AND BACKWARD SUBSTITUTION            *
C*                                                         *
C*********        **          ******************************
      SUBROUTINE CDB2A(A,B,PZERO,INDEP,N,K,KIND,ICON,ND)
      DIMENSION A(ND,10),B(ND),INDEP(ND)
      DOUBLE PRECISION SUM1,THR
      ICON=0
      THR=PZERO
      IF(N.GE.2) GO TO 1110
      IF(N.LE.0) GO TO 1340
      IF(KIND.EQ.1) GO TO 1100
      B(1)=B(1)/A(1,K)
      GO TO 1340
 1100 INDEP(1)=1
      GO TO 1340
 1110 IF(KIND.EQ.2) GO TO 1260
      DO 1120 IP=1,N
      INDEP(IP)=1
 1120 CONTINUE
C     CHOLESKI DECOMPOSITION
      DO 1250 IP=1,N
      SUM1=A(IP,K)
      IF(IP.EQ.1) GO TO 1140
      ICL=1
      IF(IP.LT.K) ICL=K-IP+1
      ICU=K-1
      DO 1130 IC=ICL,ICU
      SUM1=SUM1-DBLE(A(IP,IC))**2
 1130 CONTINUE
 1140 IF(SUM1.GT.THR) GO TO 1190
      INDEP(IP)=0
      A(IP,K)=1.0
      ICL=IP-K+1
      IF(ICL.LT.1) ICL=1
      ICU=IP-1
      IF(ICU.LT.1) GO TO 1160
      DO 1150 IC=ICL,ICU
      ICS=IC-IP+K
      A(IP,ICS)=0.0
 1150 CONTINUE
 1160 IRU=IP+K-1
      IF(IRU.GT.N) IRU=N
      IRL=IP+1
      IF(IRL.GT.N) GO TO 1180
      DO 1170 IR=IRL,IRU
      ICS=IP-IR+K
      A(IR,ICS)=0.0
 1170 CONTINUE
 1180 ICON=-IP
      GO TO 1240
 1190 A(IP,K)=DSQRT(SUM1)
      IRL=IP+1
      IF(IRL.GT.N) GO TO 1240
      IRU=IP+K-1
      IF(IRU.GT.N) IRU=N
      DO 1230 IR=IRL,IRU
      IP1=IP+K-IR
      SUM1=A(IR,IP1)
      ICU=IP-1
      IF(ICU.LT.1) GO TO 1220
      ICL=IP-K+1
      IF(ICL.LT.1) ICL=1
      DO 1210 IC=ICL,ICU
      IC1=IC+K-IP
      IC2=IC+K-IR
      IF(IC2.LT.1) GO TO 1200
      SUM1=SUM1-DBLE(A(IP,IC1))*DBLE(A(IR,IC2))
 1200 CONTINUE
 1210 CONTINUE
 1220 A(IR,IP1)=SUM1/DBLE(A(IP,K))
 1230 CONTINUE
 1240 CONTINUE
 1250 CONTINUE
C     FOWARD AND BACKWARD SUBSTITUTION
 1260 IF(KIND.EQ.1) GO TO 1340
      DO 1270 IP=1,N
      IF(INDEP(IP).EQ.0) B(IP)=0.0
 1270 CONTINUE
      NP1=N+1
C     FOWARD SUBSTITUTION
      DO 1300 IP=1,N
      SUM1=B(IP)
      ICU=IP-1
      IF(ICU.LT.1) GO TO 1290
      ICL=IP-K+1
      IF(ICL.LT.1) ICL=1
      DO 1280 IC=ICL,ICU
      IC1=IC+K-IP
      SUM1=SUM1-DBLE(A(IP,IC1))*DBLE(B(IC))
 1280 CONTINUE
 1290 B(IP)=SUM1/DBLE(A(IP,K))
 1300 CONTINUE
C     BACKWARD SUBSTITUTION
      DO 1330 IPD=1,N
      IP=NP1-IPD
      SUM1=B(IP)
      ICL=IP+1
      IF(ICL.GT.N) GO TO 1320
      ICU=IP+K-1
      IF(ICU.GT.N) ICU=N
      DO 1310 IC=ICL,ICU
      IC1=IP+K-IC
      IP1=IC
      SUM1=SUM1-DBLE(A(IP1,IC1))*DBLE(B(IC))
 1310 CONTINUE
 1320 B(IP)=SUM1/DBLE(A(IP,K))
 1330 CONTINUE
 1340 RETURN
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
      call ludcmpx(a,n,nmax,indx,det)
      do i=1,n
         det=det*a(i,i)
      enddo
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
      
      dimension a(nmax,nmax),y(nmax,nmax),indx(nmax)
c     ------------------------------------------------------
c
c     --- Set up identity matrix
c
      do i=1,n
         do j=1,n
            y(i,j)=0.d0
         enddo
         y(i,i)=1.d0
      enddo
c
c     --- Decompose the matrix 
c
      call ludcmpx(a,n,nmax,indx,d)
c
c     --- Find inverse
c
      do i=1,n
         call lubksbx(a,n,nmax,indx,y(1,i))
      enddo
c
c     --- y -> a
c
      do i=1,n
         do j=1,n
            a(i,j)=y(i,j)
         enddo
      enddo
c     ------------------------------------------------------
      return
      end
c**************************************************************
      subroutine matsolv(a,x,n,nmax,nx)
c     -------------------------------------------------------
c     Solve Linear Equation by using LU decomposition
c     -------------------------------------------------------
c     a       ... Matrix (IN/OUT)
c     n,nmax  ... Size of Matrix, Size of Array Declaration
c
      implicit real*8(a-h,o-z)

      dimension a(nmax,nmax),x(nmax,nx),indx(nmax)
c     -------------------------------------------------------
      call ludcmpx(a,n,nmax,indx,d)
      do i=1,nx
         call lubksbx(a,n,nmax,indx,x(1,i))
      enddo
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
C----------------------------------------------------------------
      real*8 FUNCTION xERF(X)
C     ERROR FUNCTION 
C                                   WRITTEN BY I.NINOMIYA   1990.10.10  
C     ----------------------------------------------------------------
      implicit real*8(a-h,o-z)
      PARAMETER(BIG=4.1)
C     ------------------------------------------------------------- 
      DATA A0/ 2.4646979237E+01/
      DATA A1/-6.1755215006E+01/
      DATA A2/-7.1903866174E+00/
      DATA A3/-1.0539111539E+00/
      DATA A4/ 4.8348018815E-03/
      DATA B0/ 1.9198579814E+02/
      DATA B1/ 8.1444739854E+01/
      DATA B2/ 1.3858875126E+01/
      DATA C0/ 6.052726918657E+02/
      DATA D0/ 5.820282538495E+02/
      DATA D1/ 9.805022334795E+01/
      DATA D2/-2.799556174376E+01/
      DATA D3/ 3.221171518211E+01/
      DATA D4/-6.840230916011E+00/
C     --------ENTRY POINT--------------------------------------------
      T=ABS(X)
      IF (T.LE.2.0) THEN
C     --------RATIONAL APPROXIMATION TO ERF(X) IN (-2<=X<=2)---------
           U=T*T
           Y=((((A4*U+A3)*U+A2)*U+A1)*U+A0)*T/(((U+B2)*U+B1)*U+B0)+T 
      ELSE IF (T.LE.BIG) THEN
C     ---RATIONAL APPROXIMATION TO ERF(X) IN (X>2)(ABS.ERROR.CRIT.)--
           Y=C0/(((((T+D4)*T+D3)*T+D2)*T+D1)*T+D0) 
           Y=Y*Y 
           Y=Y*Y 
           Y=Y*Y 
           Y=1.0-Y*Y 
      ELSE
           Y=1.0
      ENDIF
      xERF=SIGN(Y,X) 
      RETURN
      END 
