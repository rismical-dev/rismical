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
      dimension ip(0:2+ngrid/2),work(0:ngrid/2-1),twork(0:ngrid*5/8-1)

      pi=dacos(-1.d0)
      deltak=pi/(dble(ngrid)*rdelta)
c---------------------------------------------------------
      ip(0)=0
c
c     --- for ix>0  (Normal FFT)
c
      if (ix.gt.0) then
         dum(1)=0.d0
         do i=1,ngrid-1
            rr=rdelta*dble(i)
            dum(i+1)=rr*func(i)
         enddo

c     Ooura FFT
         call dfst(ngrid,dum,work,ip,twork)

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
         
c     Ooura FFT
         call dfst(ngrid,dum,work,ip,twork)
         
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
