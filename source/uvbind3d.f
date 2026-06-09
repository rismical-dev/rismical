c----------------------------------------------------------------
c     Calculate Site Contribution to Solute-Solvent Binding Energy
c----------------------------------------------------------------
      subroutine uvbind3d(namef,ng3d,n2uq,cr,tr)
c     
      implicit real*8 (a-h,o-z)
      complex*16 cr
      character*256 namef

      include "rismrun.i"
      include "solvent.i"
      include "solute.i"
      include "phys_const.i"

      dimension cr(ng3d,n2uq)
      dimension tr(ng3d,n2uq)
      
      real*8 ,allocatable :: ebindes(:,:),ebindlj(:,:)
      real*8 ,allocatable :: epsig6(:,:)
      real*8 ,allocatable :: epsig12(:,:)

      allocate (ebindlj(nvuq,nu))
      allocate (ebindes(nvuq,nu))
      allocate (epsig6(nu,nvuq),epsig12(nu,nvuq))

      rd33=rdelta3d**3
c-----------------------------------------------------------------
c
c     --- LJ Parameter Settings
c
      call rsmljcomb(6,nu,nvuq
     &     ,siglju,epslju,sigljvuq,epsljvuq,epsig6)
      call rsmljcomb(12,nu,nvuq
     &     ,siglju,epslju,sigljvuq,epsljvuq,epsig12)
C
C
c     --- Initialize ebind
c
      do i=1,nu
         do j=1,nvuq
            ebindes(j,i)=0.d0
            ebindlj(j,i)=0.d0
         enddo
      enddo
c     
c     --- Calculate U-V Binding
c     
      k0=ngrid3d/2+1
c
      do j=1,nvuq
         do i=1,nu

            ebindlj_loc=0.d0
            ebindes_loc=0.d0

!$OMP PARALLEL DO REDUCTION(+:ebindlj_loc,ebindes_loc)
!$OMP&   PRIVATE(kx,ky,k,rx,ry,rz,rr2,rr,grd,
!$OMP&           rrinv2,rrinv6,rrinv12,rr6,rr12)
!$OMP&   SCHEDULE(static)
            do kz=1,ngrid3d
            do ky=1,ngrid3d
            do kx=1,ngrid3d

               k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2

               rx=rdelta3d*dble(kx-k0)-xyzu(1,i)
               ry=rdelta3d*dble(ky-k0)-xyzu(2,i)
               rz=rdelta3d*dble(kz-k0)-xyzu(3,i)
               rr2=rx**2+ry**2+rz**2
               rr=dsqrt(rx**2+ry**2+rz**2)

               if (rr2.gt.1.d-10) then

                  grd=(dble(cr(k,j))+tr(k,j)+1.d0)
     &                 *rd33*densuq(j)
c     --- LJ
                  rrinv2=1.d0/rr2
                  rrinv6=rrinv2*rrinv2*rrinv2
                  rrinv12=rrinv6*rrinv6
                  rr6=epsig6(i,j)*rrinv6
                  rr12=epsig12(i,j)*rrinv12

                  ebindlj_loc=ebindlj_loc+4.d0*(rr12-rr6)*grd

c     --- Electro Static
                  ebindes_loc=ebindes_loc
     &                 +qu(i)*q2uq(j)/rr*fel*grd

               endif
C
            enddo
            enddo
            enddo
!$OMP END PARALLEL DO

            ebindlj(j,i)=ebindlj_loc
            ebindes(j,i)=ebindes_loc

         enddo
      enddo
c
c     --- print out 
c
      ift=45
      open (ift,file=namef)
      write(ift,9997)
      do i=1,nu
         write(ift,9998) (ebindlj(j,i)*1.d-3,j=1,nvuq)
     &                  ,(ebindes(j,i)*1.d-3,j=1,nvuq)
      enddo
      close(ift)
c----------------------------------------------------------------
      return
 9997 format ("## Solute-solvent interaction energy in [kJ/mol]"
     &     ,/,"## ((E_lj(j,i),j=1,nv),(E_es(j,i),j=1,nv),i=1,nu)")
 9998 format (10F20.8)
      end
