c----------------------------------------------------------------
c     Make 3D F-Bond for 3D-UV
c----------------------------------------------------------------
      subroutine fbond3duv(ng3d,fr,fk)
c     
c     ir            ... file number of output (STDOUT)
c     iw            ... file number of input (STDIN)
c     ng3d          ... ngrid3d**3
c     rdelta3d      ... grid width of r-space [Angstrom]
c     alp           ... parameter of f-bond
c     
      implicit real*8 (a-h,o-z)
      complex*16 fk
C      
      include "rismrun.i"
      include "solute.i"
      include "phys_const.i"
c
      dimension fr(ng3d),fk(ng3d)

      dk3d=2.d0*pi/(rdelta3d*dble(ngrid3d))
c-----------------------------------------------------------------
c     
c     --- Initialize
c     
      do k=1,ng3d
         fr(k)=0.d0
         fk(k)=(0.d0,0.d0)
      enddo
c
      if (pbc) return
c     
c     --- make 3d-phi bond
c     
!$omp parallel do default(none)
!$omp&private(kz,ky,kx,k,rx,ry,rz,rkx,rky,rkz,rk, 
!$omp&           i,rix,riy,riz,rri,rrik) 
!$omp&   shared(ngrid3d,rdelta3d,dk3d,nu,xyzu,qu,fr,fk)
      do kz=1,ngrid3d
      do ky=1,ngrid3d
      do kx=1,ngrid3d
         
         k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
         
         fr(k)=0.d0
         fk(k)=(0.d0,0.d0)
         
         rx=rdelta3d*dble(kx-1-ngrid3d/2)
         ry=rdelta3d*dble(ky-1-ngrid3d/2)
         rz=rdelta3d*dble(kz-1-ngrid3d/2)
         
         rkx=dk3d*(dble(kx-1-ngrid3d/2)+0.5d0)
         rky=dk3d*(dble(ky-1-ngrid3d/2)+0.5d0)
         rkz=dk3d*(dble(kz-1-ngrid3d/2)+0.5d0)
         
         rk=dsqrt(rkx**2+rky**2+rkz**2)
         
         do i=1,nu
            
            rix=xyzu(1,i)
            riy=xyzu(2,i)
            riz=xyzu(3,i)
            
            rri=dsqrt((rx-rix)**2+(ry-riy)**2+(rz-riz)**2)
c$$$c    --- erf version ---
c$$$            if (rri.lt.1.d-5) then
c$$$               fr(k)=fr(k)+qu(i)*alp3d*2.d0/dsqrt(pi)
c$$$            else
c$$$               fr(k)=fr(k)+qu(i)/rri*erf(alp3d*rri)
c$$$            endif
c     --- exp version ---
            if (rri.lt.1.d-5) then
               fr(k)=fr(k)+qu(i)
            else
               fr(k)=fr(k)+qu(i)/rri*(1.d0-dexp(-rri))
            endif
            
            rrik=rix*rkx+riy*rky+riz*rkz

c$$$c     --- erf version ---
c$$$            fk(k)=fk(k)+4.d0*pi*qu(i)/rk**2
c$$$     &                 *dexp(-(rk/2.d0/alp3d)**2)
c$$$     &                 *cdexp(dcmplx(0.d0,rrik))
c     --- exp version ---
            fk(k)=fk(k)+4.d0*pi*qu(i)/(rk**2*(rk**2+1.d0))
     &           *cdexp(dcmplx(0.d0,rrik))
            
         enddo
         
         fr(k)=fr(k)*fel     
         fk(k)=fk(k)*fel     
         
      enddo
      enddo
      enddo
!$omp end parallel do
c
c----------------------------------------------------------------
      return
      end
