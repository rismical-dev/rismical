c----------------------------------------------------------------
c     Calculate Site Contribution to Solute-Solvent Binding Energy
c----------------------------------------------------------------
      subroutine uvbind3d(ng3d,n2,n2uq,vres,urlj,cr,tr,listcore)
c     
      implicit real*8 (a-h,o-z)
      complex*16 cr

      include "rismrun.i"
      include "solvent.i"
      include "solute.i"
      include "phys_const.i"

      dimension vres(ng3d)
      dimension urlj(ng3d,n2uq)
      dimension listcore(ng3d)
      dimension cr(ng3d,n2uq)
      dimension tr(ng3d,n2uq)
      
      real*8 ,allocatable :: ebind(:,:)
      real*8 ,allocatable :: epsig6(:,:)
      real*8 ,allocatable :: epsig12(:,:)

      allocate (ebind(nv,nu))
      allocate (epsig6(nu,nv),epsig12(nu,nv))


      rd33=rdelta3d**3
c-----------------------------------------------------------------
c
c     --- LJ Parameter Settings
c
      call rsmljcomb(6,nu,nv
     &     ,siglju,epslju,sigljv,epsljv,epsig6)
      call rsmljcomb(12,nu,nv
     &     ,siglju,epslju,sigljv,epsljv,epsig12)
C
C
c     --- Initialize ebind
c
      do i=1,nu
         do j=1,nv
            ebind(j,i)=0.d0
         enddo
      enddo
c     
c     --- Calculate U-V Binding
c     
      k0=ngrid3d/2+1
c
      do j=1,nv
         jj=abs(iuniq(j))
         do i=1,nu
c     
c     --- LJ
c     
            do kz=1,ngrid3d
            do ky=1,ngrid3d
            do kx=1,ngrid3d

               k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2

               if (listcore(k).eq.0) goto 6000

               rx=rdelta3d*dble(kx-k0)-xyzu(1,i)
               ry=rdelta3d*dble(ky-k0)-xyzu(2,i)
               rz=rdelta3d*dble(kz-k0)-xyzu(3,i)
               rr2=rx**2+ry**2+rz**2

               rrinv2=1.d0/rr2

               rrinv6=rrinv2*rrinv2*rrinv2
               rrinv12=rrinv6*rrinv6
               rr6=epsig6(i,j)*rrinv6
               rr12=epsig12(i,j)*rrinv12
               gr=dble(cr(k,jj))+tr(k,jj)+1.d0
               if (abs(gr).lt.g0limit) gr=0.d0
               ebind(j,i)=ebind(j,i)+4.d0*(rr12-rr6) 
     &              *gr*rd33*dens(nspc(j))

 6000          continue

            enddo
            enddo
            enddo
c     
c     --- Electro Static
c     
            do kz=1,ngrid3d
            do ky=1,ngrid3d
            do kx=1,ngrid3d
C
               k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
C
               if (listcore(k).eq.0) goto 6100
C
               rx=rdelta3d*dble(kx-k0)-xyzu(1,i)
               ry=rdelta3d*dble(ky-k0)-xyzu(2,i)
               rz=rdelta3d*dble(kz-k0)-xyzu(3,i)
               rr=dsqrt(rx**2+ry**2+rz**2)
C
               gr=dble(cr(k,jj))+tr(k,jj)+1.d0
               if (abs(gr).lt.g0limit) gr=0.d0

               ebind(j,i)=ebind(j,i)
     &              +qu(i)*qv(j)/rr*fel  
     &              *gr*rd33*dens(nspc(j))
C
 6100          continue
C
            enddo
            enddo
            enddo

         enddo
      enddo
c
c     --- print out 
c
      write(*,9997)
      ebindtot=0.d0
      do i=1,nu
         esum=0.d0
         do j=1,nv
            esum=esum+ebind(j,i)
         enddo
         ebindtot=ebindtot+esum
         write(*,9998) i,esum,(ebind(j,i)*1.d-3,j=1,nv)
      enddo
      write(*,9999) ebindtot*1.d-3
c----------------------------------------------------------------
      return
 9997 format (/,4x,"======== UV BINDING ENERGY COMPONENT =========")
 9998 format (4x,I5,10F20.8)
 9999 format (/,4x,"----------------------------------------------"
     &       ,/,4x,"TOTAL UV BINDING ENERGY : ",F20.8,"[kJ/mol]")
      end
