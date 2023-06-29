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

      include "rismio.i"
      include "rismrun.i"
      include "solvent.i"
      include "solute.i"
      include "phys_const.i"
      
      dimension rk(0:ngr1d)
      dimension xvv(nx22,n2uq,n2uq),xvk(0:ngr1d,n2uq,n2uq)
      dimension listxvv(ngr3d/2+1,ngr3d/2+1,ngr3d/2+1)
      dimension yd(0:ngr1d)
      
C----------------------------------------------------------------
      yd=0.d0
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
c     --- setup hvv grid point
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
      read(ift,*) char2
      read(ift,*) char2
      read(ift,*) char2,ndum
      do i=1,ndum
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
c     
c     --- make 3D-Xvv by interpolating 1D-Xvv
c     
      do j=1,nvuq
         do i=1,nvuq
            
            sum=0.d0
            
            ill=3
            
            do kz=1,ngrid3d/2
            do ky=1,ngrid3d/2
            do kx=1,ngrid3d/2

               rkx=dble(kx)-0.5d0
               rky=dble(ky)-0.5d0
               rkz=dble(kz)-0.5d0
               rk3=dsqrt(rkx*rkx+rky*rky+rkz*rkz)*dk3d
               
               dum=hrho(ngrid,deltak,xvk(0,i,j),rk,rk3,yd,ill)
               
               k=listxvv(kx,ky,kz)
               xvv(k,i,j)=dum
               
            enddo
            enddo
            enddo
            
         enddo                  ! of j
      enddo                     ! of i
c$$$c
c$$$c     --- Non-reduced xvv version
c$$$c
c$$$      do i=1,nv
c$$$         do j=1,nv
c$$$            
c$$$            sum=0.d0
c$$$            
c$$$            do k=1,ngrid
c$$$               read(ift,*) dum
c$$$               xv1d(k)=dum
c$$$            enddo
c$$$            
c$$$            xv1d(0)=xv1d(1) ! to check later
c$$$            
c$$$            ill=3
c$$$            
c$$$            do kz=1,ngrid3d/2
c$$$            do ky=1,ngrid3d/2
c$$$            do kx=1,ngrid3d/2
c$$$
c$$$               rkx=dble(kx)-0.5d0
c$$$               rky=dble(ky)-0.5d0
c$$$               rkz=dble(kz)-0.5d0
c$$$               rk3=dsqrt(rkx*rkx+rky*rky+rkz*rkz)*dk3d
c$$$               
c$$$               dum=hrho(ngrid,deltak,xv1d,rk,rk3,yd,ill)
c$$$               
c$$$               k=listxvv(kx,ky,kz)
c$$$               xvv(k,i,j)=dum
c$$$               
c$$$            enddo
c$$$            enddo
c$$$            enddo
c$$$            
c$$$         enddo                  ! of j
c$$$      enddo                     ! of i
c$$$      close(ift)
c----------------------------------------------------------------

 8000 continue
c
      return
      end
