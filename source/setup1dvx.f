c**************************************************************
c------------------------------------------------------------------
c     Read k-space 1D-VV Total Correlation Function for 3D-RISM
c------------------------------------------------------------------
      subroutine setup1dvx(ngr1d,n2,nx22,ngr3d,listxvv
     &                    ,xvv,wkvv,wk2)
c
c     ngrid3d   ... number of grid of 3d-rdf
c     nv        ... number of site of solvent
c     xvv       ... (wk+dens*hvk)
c     wk2       ... intra-molecular correlation function of solvent
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
      
      dimension wk2(ngr1d,n2,n2),rk(0:ngr1d)
      dimension xvv(nx22,n2,n2),hv(0:ngr1d)
      dimension wkvv(nx22,n2,n2)
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
c     --- read solvent h from file and make xvv
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
      do i=1,nv
         do j=1,nv
            
            sum=0.d0
            
            do k=1,ngrid
               read(ift,*) dum
               hv(k)=wk2(k,i,j)+dens(nspc(i))*dum
            enddo
            
            hv(0)=hv(1) ! to check later
            
            ill=3
            
            do kz=1,ngrid3d/2
            do ky=1,ngrid3d/2
            do kx=1,ngrid3d/2

               rkx=dble(kx)-0.5d0
               rky=dble(ky)-0.5d0
               rkz=dble(kz)-0.5d0
               rk3=dsqrt(rkx*rkx+rky*rky+rkz*rkz)*dk3d
               
               dum=hrho(ngrid,deltak,hv,rk,rk3,yd,ill)
               
               k=listxvv(kx,ky,kz)
               xvv(k,i,j)=dum
               
            enddo
            enddo
            enddo
            
         enddo                  ! of j
      enddo                     ! of i
      close(ift)

c----------------------------------------------------------------

      do i=1,nv
         do j=1,nv
            
            sum=0.d0
            do k=1,ngrid
               hv(k)=wk2(k,i,j)
            enddo

            if (nspc(i).eq.nspc(j)) then
               hv(0)=1.d0
            else
               hv(0)=0.d0
            endif
            
            ill=3
            
            do kz=1,ngrid3d/2
            do ky=1,ngrid3d/2
            do kx=1,ngrid3d/2

               rkx=dble(kx)-0.5d0
               rky=dble(ky)-0.5d0
               rkz=dble(kz)-0.5d0
               rk3=dsqrt(rkx*rkx+rky*rky+rkz*rkz)*dk3d
               
               dum=hrho(ngrid,deltak,hv,rk,rk3,yd,ill)
               
               k=listxvv(kx,ky,kz)
               wkvv(k,i,j)=dum
               
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
