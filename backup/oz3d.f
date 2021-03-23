c**************************************************************
c-----------------------------------------------------------------
c     3-Dimensional Ornstein-Zernike equation
c-----------------------------------------------------------------
      subroutine oz3d(ngr3d,nx22,n2,n2uq
     &               ,listxvv,ck,cr,tr,xvv)
c     
c     ngrid3d   ... number of grid of 3d-rdf
c     nxvv      ... number of grid points of xvv
c     rdelta3d  ... grid width of r-space for 3d grid [angstrom]
c     nv        ... number of site of solvent
c     nvuq      ... number of symmetry uniq site of solvent
c     ck        ... (in)  k-space direct correlation function
c     cr        ... (out) r-space direct correlation function
c     tr        ... (in)  r-space direct correlation function
c               ... (out) tau bond =hr-cr
c     xvv       ... (wk2+dens*hvk)
c     iuniq     ... symmetry uniq flag
c     listxvv   ... list vector of xvv
c     
      implicit real*8 (a-h,o-z)
      complex*16 cdum,ck,cr,dsum
      complex*16 ,allocatable :: dumfft(:)

      include "solvent.i"
      include "rismrun.i"
      include "phys_const.i"

      dimension ck(ngr3d**3,n2uq)
      dimension cr(ngr3d**3,n2uq)
      dimension tr(ngr3d**3,n2uq)
      dimension xvv(nx22,n2,n2)
      dimension listxvv(ngr3d/2+1,ngr3d/2+1,ngr3d/2+1)
      dimension cdum(n2uq)

      ng3d=ngrid3d**3
      allocate (dumfft(ng3d))

      ngrid3d2=ngrid3d*ngrid3d
      dnshift=dble(ngrid3d+1)/2.d0 
      ngshift=ngrid3d/2+1
      dk3d=2.d0*pi/(rdelta3d*dble(ngrid3d))
c-----------------------------------------------------------------
c
c     --- parameter for fft
c
      m=nint(dlog(dble(ngrid3d))/dlog(2.d0))
c
      call vclr_mp(cr,1,nvuq*ngrid3d**3*2)
c
c     --- k-space 3d uv-oz [cr(k) --> hr(k)]
c            
!$omp parallel do private(nkz,nky,nkx,k,kxvv,cdum,jj,dsum,jj2)
      do kz=1,ngrid3d
      nkz=nint(abs(dble(kz)-dnshift)+0.5d0)
      do ky=1,ngrid3d
      nky=nint(abs(dble(ky)-dnshift)+0.5d0)
      do kx=1,ngrid3d
      nkx=nint(abs(dble(kx)-dnshift)+0.5d0)
         
         k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d2
         kxvv=listxvv(nkx,nky,nkz)

         do j=1,nvuq
            cdum(j)=ck(k,j)     ! here, c(k) is in "ck"
         enddo
         
         do j=1,nv
            
            jj=iuniq(j)
            if (jj.le.0) goto 8000

            dsum=(0.d0,0.d0)

            do j2=1,nv

               jj2=abs(iuniq(j2))
               dsum=dsum+cdum(jj2)*dcmplx(xvv(kxvv,j2,j),0.d0)

            enddo
            
            cr(k,jj)=dsum       ! here, hr(k) is in "cr"

 8000       continue

         enddo

      enddo
      enddo
      enddo
!$omp end parallel do
c
c     --- fourier transform [h(k) --> h(r) --> t(r)]
c     
      do j=1,nvuq
         
         call vclr_mp(dumfft,1,ng3d*2)

!$omp parallel do 
         do k=1,ng3d

            dumfft(k)=cr(k,j)

         enddo
!$omp end parallel do

         inv=0
         call ffte3d(dumfft,ngrid3d,rdelta3d,inv)
         
!$omp parallel do private(rz,ry,rx,k,dkr)

         do kz=1,ngrid3d
         rz=rdelta3d*dble(kz-ngshift)
         do ky=1,ngrid3d
         ry=rdelta3d*dble(ky-ngshift)
         do kx=1,ngrid3d
         rx=rdelta3d*dble(kx-ngshift)

            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d2
            dkr=dk3d/2.d0*(rx+ry+rz) 
            cr(k,j)=dcmplx(tr(k,j),0.d0)
            tr(k,j)=dble(dumfft(k)*cdexp(dcmplx(0.d0,-dkr)))
     *           -tr(k,j)
                                ! here, c(r) is in "cr"
                                !   and t(r) is in "tr"
         enddo
         enddo
         enddo     ! of kz
!$omp end parallel do

      enddo                     ! of j for nvuq
c-------------------------------------------------------------------
      deallocate (dumfft)
c
      return
      end
