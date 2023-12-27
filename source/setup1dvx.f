c**************************************************************
c------------------------------------------------------------------
c     Read k-space 1D-VV Total Correlation Function for 3D-RISM
c------------------------------------------------------------------
      subroutine setup1dvx(ngr1d,n2uq,nx22,ngr3d,listxvv,xvv)
c
c     ngrid3d   ... number of grid of 3d-rdf
c     nv        ... number of site of solvent
c     xvv       ... (wk+dens*hvk)
c     dens      ... number density of solvent 
c     listxvv   ... list vector of xvv
c     rdelta    ... grid width of 1d-rdf
c     rdelta3d  ... grid width of 3d-rdf
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
      
      dimension rk(0:ngr1d)
      dimension xvv(nx22,n2uq,n2uq),xvk(0:ngr1d,n2uq,n2uq)
      dimension listxvv(ngr3d/2+1,ngr3d/2+1,ngr3d/2+1)
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
      n=ngrid3d/2  
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
      do i2=1,nvuq
      do i1=1,nvuq
         xvk(0,i1,i2)=xvk(1,i1,i2)
      enddo
      enddo
      close(ift)
c     
c     --- make 3D-Xvv by interpolating 1D-Xvv
c     
      do j=1,nvuq
         do i=1,nvuq
            
c
c           calculate derivative of xvk
c
            call dpchsp(ic,vc,ngrid+1,rk,xvk(0,i,j),xvkd,1,wk,nwk,ierr)
c
            sum=0.d0
            
            ill=3
            
            do kz=1,ngrid3d/2
            do ky=1,ngrid3d/2
            do kx=1,ngrid3d/2

               rkx=dble(kx)-0.5d0
               rky=dble(ky)-0.5d0
               rkz=dble(kz)-0.5d0
               rk3=dsqrt(rkx*rkx+rky*rky+rkz*rkz)*dk3d
               
               call dpchfe(ngrid+1,rk,xvk(0,i,j),xvkd
     &                     ,1,skip,1,rk3,dum,ierr)
               
               k=listxvv(kx,ky,kz)
               xvv(k,i,j)=dum
               
            enddo
            enddo
            enddo
            
         enddo                  ! of j
      enddo                     ! of i
c----------------------------------------------------------------

 8000 continue
c
      return
      end
