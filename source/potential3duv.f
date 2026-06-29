c----------------------------------------------------------------
c     Make intraction potential for 3D-UV
c----------------------------------------------------------------
      subroutine potential3duv(ng3d,n2uq,vres,urlj,tpbc,listcore)
c     
c     ng3d=ngrid3d**3 ... number of grid of 3D-RDF
c     n2uq=nvuq       ... number of symmetry uniq site of solvent
c     vres            ... electro static potential
c     urlj            ... LJ potential energy 
c     
      implicit real*8 (a-h,o-z)
c
      include "rismrun.i"
      include "rismio.i"
      include "solvent.i"
      include "solute.i"
      include "phys_const.i"
c
      real*8 ,allocatable :: epsig6(:,:)
      real*8 ,allocatable :: epsig12(:,:)
c      
      dimension vres(ng3d)
      dimension tpbc(ng3d)
      dimension urlj(ng3d,n2uq)
      dimension listcore(ng3d)
c     
      allocate (epsig6(nu,nv),epsig12(nu,nv))
c-----------------------------------------------------------------
c
c     --- LJ Parameter Settings
c
      call rsmljcomb(6,nu,nv
     &     ,siglju,epslju,sigljv,epsljv,epsig6)
      call rsmljcomb(12,nu,nv
     &     ,siglju,epslju,sigljv,epsljv,epsig12)
C
C     --- Initialize vres and urlj
C
      call vclr_mp(vres,1,ng3d)
      call vclr_mp(urlj,1,ng3d*nvuq)
C
c     --- Setup Core Region List
c
      do k=1,ng3d
         listcore(k)=1
      enddo

!$omp parallel do default(none) 
!$omp&private(kz,ky,kx,k,i,rx,ry,rz,rr) 
!$omp&shared(ngrid3d,rdelta3d,xyzu,nu,siglju,listcore,rcore3d)
      do kz=1,ngrid3d
      do ky=1,ngrid3d
      do kx=1,ngrid3d
         k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
         do i=1,nu
            rx=rdelta3d*dble(kx-1-ngrid3d/2)-xyzu(1,i)
            ry=rdelta3d*dble(ky-1-ngrid3d/2)-xyzu(2,i)
            rz=rdelta3d*dble(kz-1-ngrid3d/2)-xyzu(3,i)
            rr=dsqrt(rx**2+ry**2+rz**2)
            if (rr.lt.siglju(i)*rcore3d) then
               listcore(k)=0
            endif
         enddo
         enddo
         enddo
      enddo
!$omp end parallel do
c     
c     --- LJ
c     
      if (pbc) then
         call cal_ulj_pbc(urlj,epsig6,epsig12,ng3d,n2uq,listcore)
      else
         call cal_ulj(urlj,epsig6,epsig12,ng3d,n2uq,listcore)
      endif
c     
c     --- Electro Static 
c     
c
c     By partial charge
c
      if (ipot3d.eq.0) then

         if (pbc) then
            call cal_vres_pbc(vres,tpbc,listcore,ng3d)
         else
            call cal_vres(vres,listcore,ng3d)
         endif
c
c     Read from external file
c
      elseif (ipot3d.eq.1) then

         if (pbc) goto 9000
         call readespmap(espfile,vres,rdelta3d,ngrid3d)

      elseif (ipot3d.eq.2) then

         if (pbc) goto 9000
         call readespcube(espfile,vres,rdelta3d,ngrid3d)

      else
         goto 9000
      endif
c-------------------------------------------------------------------
      deallocate (epsig6,epsig12)
      RETURN
 9000 continue
      write(*,*) "Error. Invalid potential option."
      ierr=700
      call abrt(ierr)
c-------------------------------------------------------------------
 9999 format (/,4x,"======== Short Range Potential =========",/,
     &        /,4x,"Lennard-Jones (12-6) Potential")
      end
c-------------------------------------------------------------------
c     calculate LJ potential
c-------------------------------------------------------------------
      subroutine cal_ulj(urlj,epsig6,epsig12,ng3d,n2uq,listcore)
c
      implicit real*8 (a-h,o-z)
c
      include "rismrun.i"
      include "rismio.i"
      include "solvent.i"
      include "solute.i"
      include "phys_const.i"

      dimension urlj(ng3d,n2uq)
      dimension epsig6(nu,nv),epsig12(nu,nv)
      dimension listcore(ng3d)
c-------------------------------------------------------------------
      k0=ngrid3d/2+1

      do j=1,nv
         jj=iuniq(j)
         if (jj.lt.0) goto 6500
         do i=1,nu
c     
c     --- LJ
c     
!$omp parallel do default(none) 
!$omp&private(kz,ky,kx,k,rx,ry,rz,rr2,rrinv2,rrinv6,rrinv12,rr6,rr12) 
!$omp&shared(ngrid3d,rdelta3d,xyzu,listcore)
!$omp&shared(epsig6,epsig12,urlj,i,j,jj,k0)
            do kz=1,ngrid3d
            do ky=1,ngrid3d
            do kx=1,ngrid3d

               k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2

               if (listcore(k).eq.0) then
                  
                  rr2=1.d-6

               else

                  rx=rdelta3d*dble(kx-k0)-xyzu(1,i)
                  ry=rdelta3d*dble(ky-k0)-xyzu(2,i)
                  rz=rdelta3d*dble(kz-k0)-xyzu(3,i)
                  rr2=rx**2+ry**2+rz**2

               endif

               rrinv2=1.d0/rr2
               rrinv6=rrinv2*rrinv2*rrinv2
               rrinv12=rrinv6*rrinv6

               rr6=epsig6(i,j)*rrinv6
               rr12=epsig12(i,j)*rrinv12
               urlj(k,jj)=urlj(k,jj)+4.d0*(rr12-rr6) ![J/mol]

            enddo               ! of kx
            enddo               ! of ky
            enddo               ! of kz
!$omp end parallel do

         enddo                  ! of nu

 6500    continue
C
      enddo                     ! of nv
c-------------------------------------------------------------------
      return
      end
c-------------------------------------------------------------------
c     calculate LJ potential for PBC
c-------------------------------------------------------------------
      subroutine cal_ulj_pbc(urlj,epsig6,epsig12,ng3d,n2uq,listcore)
c
      implicit real*8 (a-h,o-z)
c
      include "rismrun.i"
      include "rismio.i"
      include "solvent.i"
      include "solute.i"
      include "phys_const.i"

      dimension urlj(ng3d,n2uq)
      dimension epsig6(nu,nv),epsig12(nu,nv)
      dimension listcore(ng3d)
c-------------------------------------------------------------------
      k0=ngrid3d/2+1
      box=dble(ngrid3d)*rdelta3d

      do j=1,nv
         jj=iuniq(j)
         if (jj.lt.0) goto 6500
         do i=1,nu
c     
c     --- LJ
c     
!$omp parallel do default(none) 
!$omp&private(kz,ky,kx,k,rx,ry,rz,rr2,rrinv2,rrinv6,rrinv12,rr6,rr12) 
!$omp&shared(ngrid3d,rdelta3d,box,xyzu,listcore)
!$omp&shared(epsig6,epsig12,urlj,i,j,jj,k0)
            do kz=1,ngrid3d
            do ky=1,ngrid3d
            do kx=1,ngrid3d

               k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2

               if (listcore(k).eq.0) then
                  
                  rr2=1.d-6

               else

                  rx=rdelta3d*dble(kx-k0)-xyzu(1,i)
                  ry=rdelta3d*dble(ky-k0)-xyzu(2,i)
                  rz=rdelta3d*dble(kz-k0)-xyzu(3,i)
c              -- apply minimal image convention --
                  rx = rx - anint( rx/box)*box
                  ry = ry - anint( ry/box)*box
                  rz = rz - anint( rz/box)*box

                  rr2=rx**2+ry**2+rz**2

               endif

               rrinv2=1.d0/rr2
               rrinv6=rrinv2*rrinv2*rrinv2
               rrinv12=rrinv6*rrinv6

               rr6=epsig6(i,j)*rrinv6
               rr12=epsig12(i,j)*rrinv12
               urlj(k,jj)=urlj(k,jj)+4.d0*(rr12-rr6) ![J/mol]

            enddo               ! of kx
            enddo               ! of ky
            enddo               ! of kz
!$omp end parallel do

         enddo                  ! of nu

 6500    continue
C
      enddo                     ! of nv
c-------------------------------------------------------------------
      return
      end
c-------------------------------------------------------------------
c     calculate ESP potential
c-------------------------------------------------------------------
      subroutine cal_vres(vres,listcore,ng3d)
      implicit real*8 (a-h,o-z)
c
      include "rismrun.i"
      include "rismio.i"
      include "solvent.i"
      include "solute.i"
      include "phys_const.i"

      dimension vres(ng3d)
      dimension listcore(ng3d)
c-------------------------------------------------------------------
      k0=ngrid3d/2+1
      do i=1,nu

!$omp parallel do default(none)
!$omp&private(ky,kx,k,rx,ry,rz,rr)
!$omp&shared(ngrid3d,rdelta3d,xyzu,listcore)
!$omp&shared(vres,qu,i,k0)
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d
C
            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
C
            if (listcore(k).ne.0) then
C
               rx=rdelta3d*dble(kx-k0)-xyzu(1,i)
               ry=rdelta3d*dble(ky-k0)-xyzu(2,i)
               rz=rdelta3d*dble(kz-k0)-xyzu(3,i)
               rr=dsqrt(rx**2+ry**2+rz**2)
C     
               vres(k)=vres(k) + qu(i)/rr*fel  
C
            endif
C
         enddo               ! of kx
         enddo               ! of ky
         enddo               ! of kz
!$omp end parallel do
C
      enddo                     ! of nu
c-------------------------------------------------------------------
      return
      end
c-------------------------------------------------------------------
c     calculate ESP potential for PBC
c-------------------------------------------------------------------
      subroutine cal_vres_pbc(vres,tpbc,listcore,ng3d)
      implicit real*8 (a-h,o-z)
      logical shift
c
      include "rismrun.i"
      include "rismio.i"
      include "solvent.i"
      include "solute.i"
      include "phys_const.i"

      complex*16,allocatable::fkwork(:)

      dimension vres(ng3d)
      dimension tpbc(ng3d)
      dimension listcore(ng3d)
c-------------------------------------------------------------------
      allocate (fkwork(ng3d))

      k0=ngrid3d/2+1
      box=dble(ngrid3d)*rdelta3d
      dk3d=2.d0*pi/(rdelta3d*dble(ngrid3d))
c
c     --- Short range (Real space) part
c
      do i=1,nu

!$omp parallel do default(none)
!$omp&private(ky,kx,k,rx,ry,rz,rr)
!$omp&shared(ngrid3d,rdelta3d,box,alp3d,xyzu)
!$omp&shared(vres,tpbc,qu,i,k0)
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d
C
            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
C
            rx=rdelta3d*dble(kx-k0)-xyzu(1,i)
            ry=rdelta3d*dble(ky-k0)-xyzu(2,i)
            rz=rdelta3d*dble(kz-k0)-xyzu(3,i)
c              -- apply minimal image convention --
            rx = rx - anint( rx/box)*box
            ry = ry - anint( ry/box)*box
            rz = rz - anint( rz/box)*box
c
            rr=dsqrt(rx**2+ry**2+rz**2)
            if (rr.lt.1.d-8) then
               vres(k)=vres(k) - qu(i)*2.d0*alp3d/dsqrt(pi) *fel
            else
               vres(k)=vres(k) + qu(i)/rr* erfc(alp3d*rr)      *fel
               tpbc(k)=tpbc(k) + qu(i)/rr*(erfc(alp3d*rr)-1.d0)*fel

            endif
C
         enddo               ! of kx
         enddo               ! of ky
         enddo               ! of kz
!$omp end parallel do
C
      enddo                     ! of nu
c
c     --- Long range (Reciprocal space) part
c
!$omp parallel do default(none)
!$omp&private(kz,ky,kx,k,rkx,rky,rkz,rk, rk2,
!$omp&           i,rix,riy,riz,rri,rrik) 
!$omp&   shared(ngrid3d,rdelta3d,dk3d,nu,xyzu,qu,fkwork,alp3d,box)
      do kz=1,ngrid3d
      do ky=1,ngrid3d
      do kx=1,ngrid3d
         
         k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
         
         fkwork(k)=(0.d0,0.d0)
         
c     --- non grid shift ver
         rkx=dk3d*(dble(kx-1-ngrid3d/2))
         rky=dk3d*(dble(ky-1-ngrid3d/2))
         rkz=dk3d*(dble(kz-1-ngrid3d/2))
         
         rk=dsqrt(rkx**2+rky**2+rkz**2)
         rk2=rkx**2+rky**2+rkz**2

         if (rk.gt.1.d-8) then
            do i=1,nu
            
               rix=xyzu(1,i)
               riy=xyzu(2,i)
               riz=xyzu(3,i)
               rrik=rix*rkx+riy*rky+riz*rkz

               fkwork(k)=fkwork(k)+4.d0*pi*qu(i)/rk**2
     &              *dexp(-(rk/2.d0/alp3d)**2)
     &              *cdexp(dcmplx(0.d0,rrik))
            enddo
         endif
c
         fkwork(k)=fkwork(k)*fel
c         
      enddo
      enddo
      enddo
!$omp end parallel do

      inv=0
      shift=.not.pbc
      call ft3dfunc(fkwork,ngrid3d,rdelta3d,inv,shift)

!$omp parallel do default(none)
!$omp&private(ky,kx,k)
!$omp&shared(ngrid3d)
!$omp&shared(vres,tpbc,fkwork)
      do kz=1,ngrid3d
      do ky=1,ngrid3d
      do kx=1,ngrid3d
C
         k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
         vres(k)=vres(k) + dble(fkwork(k))
         tpbc(k)=tpbc(k) + dble(fkwork(k))
C
      enddo                     ! of kx
      enddo                     ! of ky
      enddo                     ! of kz
!$omp end parallel do

c-------------------------------------------------------------------
      deallocate (fkwork)
      return
      end
