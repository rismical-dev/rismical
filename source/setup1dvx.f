c**************************************************************
c------------------------------------------------------------------
c     Read k-space 1D-VV Total Correlation Function for 3D-RISM
c------------------------------------------------------------------
      subroutine setup1dvx(ngr1d,n2uq,nx22,ngr3d,listxvv,xvv,delhv0)
c
c     ngrid3d   ... number of grid of 3d-rdf
c     nv        ... number of site of solvent
c     xvv       ... (wk+dens*hvk)
c     dens      ... number density of solvent 
c     listxvv   ... list vector of xvv
c     rdelta    ... grid width of 1d-rdf
c     rdelta3d  ... grid width of 3d-rdf
c     delhv0    ... q*chi/k^2 k=0 for PBC background charge
c
      implicit real*8 (a-h,o-z)
      character*1 char1
      character*2 char2
      logical skip

      include "rismio.i"
      include "rismrun.i"
      include "solvent.i"
      include "solute.i"
      include "phys_const.i"
c
      real*8 ,allocatable :: xvk2(:,:,:)
c      
      dimension rk(0:ngr1d)
      dimension xvv(nx22,n2uq,n2uq),xvk(0:ngr1d,n2uq,n2uq)
      dimension listxvv(ngr3d/2+1,ngr3d/2+1,ngr3d/2+1)
      dimension delhv0(n2uq)
c     
      dimension rwork(0:ngr1d)
c
      dimension ic(2),vc(2),xvkd(0:ngr1d)
      dimension wk(2*ngr1d+2)
      
C----------------------------------------------------------------
c     parametr for Hermite interpolation
      xvkd=0.d0
      skip=.false.
      ic(1)=3
      ic(2)=3
      nwk=2*ngrid+2
c
c     --- Make list for Xvv
c
      n=ngrid3d/2+1
      k=0
      do i=1,n
         k=k+1
         listxvv(i,i,i)=k
      enddo
      do i=1,n
         do j=1,n
            if (i.ne.j) then
               k=k+1
               listxvv(i,i,j)=k
               listxvv(i,j,i)=k
               listxvv(j,i,i)=k
            endif
         enddo
      enddo
      do i=1,n
         do j=1,i-1
            do l=1,j-1
               k=k+1
               listxvv(i,j,l)=k
               listxvv(i,l,j)=k
               listxvv(j,i,l)=k
               listxvv(j,l,i)=k
               listxvv(l,i,j)=k
               listxvv(l,j,i)=k
            enddo
         enddo
      enddo
      if (nxvv.ne.k) then
         write(*,*) "error(1492), in setup1dvx."
         call abrt
      endif
c     
c     --- setup xvv grid point
c     
      deltak=pi/(dble(ngrid)*rdelta)
      dk3d=2.d0*pi/(rdelta3d*dble(ngrid3d))

      do k=0,ngrid
         rk(k)=deltak*dble(k)
      enddo
c     
c     --- read 1D-xvv from file and make 3D-xvv
c     
      ift=45
      open (ift,file=solventxvv)
c
c     Skip header
c
      do i=1,5+nv
         read(ift,*) char2
      enddo
c
c     Read reduced xvv(k)
c
      do ig=1,ngrid
         read(ift,*) ((xvk(ig,i1,i2),i1=1,nvuq),i2=1,nvuq)
      enddo
c
c     Approximate k=0 value by k=1 value
c     Since xvk is an even function, 
c     this is to ensure that the slope is zero at k = 0
c
      do i2=1,nvuq
      do i1=1,nvuq
         xvk(0,i1,i2)=xvk(1,i1,i2)
      enddo
      enddo
c
      close(ift)
c     
c     --- make 3D-Xvv by interpolating 1D-Xvv
c     
      gshift=0.5d0                 ! Grid shift to avoid k=0
      if (pbc) gshift=1.d0

      do j=1,nvuq
         do i=1,nvuq
            
c
c           calculate derivative of xvk
c
c$$$c           Hermite interpolation     
c$$$            call dpchsp(ic,vc,ngrid+1,rk,xvk(0,i,j),xvkd,1,wk,nwk,ierr)
c           Spline interpolation
            call spline_init(ngrid,rk,xvk(0,i,j),xvkd,rwork)
c
            sum=0.d0
            
            ill=3
            
            do kz=1,ngrid3d/2+1
            do ky=1,ngrid3d/2+1
            do kx=1,ngrid3d/2+1

               rkx=dble(kx)-gshift
               rky=dble(ky)-gshift
               rkz=dble(kz)-gshift
               rk3=dsqrt(rkx*rkx+rky*rky+rkz*rkz)*dk3d
               
c$$$c           Hermite interpolation     
c$$$               call dpchfe(ngrid+1,rk,xvk(0,i,j),xvkd
c$$$     &                     ,1,skip,1,rk3,dum,ierr)
c           Spline interpolation
               call spline_eval(ngrid,rk,xvk(0,i,j),xvkd,rk3,dum)
               
               k=listxvv(kx,ky,kz)
               xvv(k,i,j)=dum

            enddo
            enddo
            enddo
            
         enddo                  ! of j
      enddo                     ! of i
c
c     --- delhv0 ---
c
      if (pbc) then
         nfit=10
         ndeg=2
         allocate (xvk2(ngrid,n2uq,n2uq))
         do iv=1,nvuq
            do jv=1,nvuq
               do ig=1,ngrid
                  xvk2(ig,iv,jv)=xvk(ig,iv,jv)
               enddo
            enddo
         enddo

         call chi_klim(ngrid,rdelta,nvuq
     &        ,q2uq,xvk2,nfit,ndeg,delhv0,ierr)
         deallocate (xvk2)

      endif
c----------------------------------------------------------------

 8000 continue
c
      return
      end
c----------------------------------------------------------------
c     k->0 extrapolation of  S_i = sum_j q_j chi_ji(k) / k^2
c----------------------------------------------------------------
      subroutine chi_klim(ngrid,rdelta,nv,qv,chi,nfit,ndeg,
     &                    sval,ierr)
      implicit real*8 (a-h,o-z)
c
      dimension qv(nv)
      dimension chi(ngrid,nv,nv)
      dimension sval(nv)
c
c     --- local workspace for the least-squares normal equations ---
c         (npmax = max polynomial terms = ndeg+1; ndeg<=4 supported)
      parameter (npmax=5)
      dimension amat(npmax,npmax),bvec(npmax),csol(npmax)
      dimension xpow(npmax)
c----------------------------------------------------------------
      pi = dacos(-1.d0)
      ierr = 0
c
c     --- sanity checks ---
      if (nfit.lt.3) then
         ierr = 1
         return
      endif
      if (nfit.gt.ngrid) then
         ierr = 2
         return
      endif
      if (ndeg.lt.1 .or. ndeg.gt.npmax-1) then
         ierr = 3
         return
      endif
      nterm = ndeg + 1
c
c     --- 1D k-grid spacing ---
      dk = pi/(rdelta*dble(ngrid))
c
c================================================================
c     loop over solvent sites i : build g_i(k) and fit vs x=k^2
c================================================================
      do i=1,nv
c
c        --- clear normal-equation accumulators ---
         do ia=1,nterm
            bvec(ia) = 0.d0
            do ja=1,nterm
               amat(ia,ja) = 0.d0
            enddo
         enddo
c
c        --- accumulate least-squares normal equations over the
c            nfit smallest-k points ---
         do ig=1,nfit
            rk  = dk*dble(ig)
            rk2 = rk*rk
c
c           --- charge-weighted susceptibility N_i(k)=sum_j q_j chi_ij ---
            an = 0.d0
            do j=1,nv
               an = an + qv(j)*chi(ig,j,i)
            enddo
c
c           --- ratio g = N_i(k)/k^2  (data point to be fitted) ---
            g = an/rk2
c
c           --- independent variable x = k^2 ---
            x = rk2
c
c           --- powers x^0, x^1, ..., x^ndeg ---
            xpow(1) = 1.d0
            do ip=2,nterm
               xpow(ip) = xpow(ip-1)*x
            enddo
c
c           --- accumulate A = X^T X  and  b = X^T g ---
            do ia=1,nterm
               bvec(ia) = bvec(ia) + xpow(ia)*g
               do ja=1,nterm
                  amat(ia,ja) = amat(ia,ja) + xpow(ia)*xpow(ja)
               enddo
            enddo
c
         enddo
c
c        --- solve the (nterm x nterm) normal equations by Gauss
c            elimination with partial pivoting ---
         call gauss_solve(amat,bvec,csol,nterm,npmax,info)
         if (info.ne.0) then
            ierr = 10 + i
            return
         endif
c
c        --- intercept (x=0 value) is the k->0 limit S_i ---
         sval(i) = csol(1)
c
      enddo
c----------------------------------------------------------------
      return
      end
c----------------------------------------------------------------
c     small dense linear solver: Gauss elimination, partial pivot
c     solves A c = b  (n x n), A,b destroyed; lda = leading dim
c----------------------------------------------------------------
      subroutine gauss_solve(a,b,c,n,lda,info)
      implicit real*8 (a-h,o-z)
      dimension a(lda,lda),b(lda),c(lda)
c----------------------------------------------------------------
      info = 0
c
c     --- forward elimination with partial pivoting ---
      do kk=1,n-1
c
c        --- find pivot ---
         ipiv = kk
         pmax = dabs(a(kk,kk))
         do ii=kk+1,n
            if (dabs(a(ii,kk)).gt.pmax) then
               pmax = dabs(a(ii,kk))
               ipiv = ii
            endif
         enddo
         if (pmax.lt.1.d-30) then
            info = 1
            return
         endif
c
c        --- swap rows kk and ipiv ---
         if (ipiv.ne.kk) then
            do jj=kk,n
               tmp        = a(kk,jj)
               a(kk,jj)   = a(ipiv,jj)
               a(ipiv,jj) = tmp
            enddo
            tmp      = b(kk)
            b(kk)    = b(ipiv)
            b(ipiv)  = tmp
         endif
c
c        --- eliminate below pivot ---
         do ii=kk+1,n
            fac = a(ii,kk)/a(kk,kk)
            do jj=kk,n
               a(ii,jj) = a(ii,jj) - fac*a(kk,jj)
            enddo
            b(ii) = b(ii) - fac*b(kk)
         enddo
c
      enddo
c
      if (dabs(a(n,n)).lt.1.d-30) then
         info = 1
         return
      endif
c
c     --- back substitution ---
      c(n) = b(n)/a(n,n)
      do ii=n-1,1,-1
         s = b(ii)
         do jj=ii+1,n
            s = s - a(ii,jj)*c(jj)
         enddo
         c(ii) = s/a(ii,ii)
      enddo
c----------------------------------------------------------------
      return
      end
