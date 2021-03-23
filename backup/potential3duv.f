c----------------------------------------------------------------
c     Make intraction potential for 3D-UV
c----------------------------------------------------------------
      subroutine potential3duv(ng3d,n2uq,vres,urlj,listcore)
c     
c     ng3d=ngrid3d**3 ... number of grid of 3D-RDF
c     n2uq=nvuq       ... number of symmetry uniq site of solvent
c     vres            ... electro static potential [erg/e]
c     urlj            ... LJ potential energy  [erg]
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
      
      dimension vres(ng3d)
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

      do kz=1,ngrid3d
      do ky=1,ngrid3d
      do kx=1,ngrid3d
         k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
         do i=1,nu
            rx=rdelta3d*dble(kx-1-ngrid3d/2)-xyzu(1,i)
            ry=rdelta3d*dble(ky-1-ngrid3d/2)-xyzu(2,i)
            rz=rdelta3d*dble(kz-1-ngrid3d/2)-xyzu(3,i)
            rr=dsqrt(rx**2+ry**2+rz**2)
            if (rr.lt.siglju(i)*0.1d0) then
               listcore(k)=0
            endif
         enddo
         enddo
         enddo
      enddo
C
c     --- Make Potential
c     
      k0=ngrid3d/2+1

      do j=1,nv
         jj=iuniq(j)
         if (jj.lt.0) goto 6500
         do i=1,nu
c     
c     --- LJ
c     
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

         enddo                  ! of nu

 6500    continue
C
      enddo                     ! of nv
c     
c     --- Electro Static 
c     
c
c     By partial charge
c
      if (ipot3d.eq.0) then

      do i=1,nu

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
            vres(k)=vres(k)
     &           +qu(i)/rr*fel  
C
 6100       continue
C
         enddo                  ! of kx
         enddo               ! of ky
         enddo               ! of kz
C
      enddo                     ! of nu
c
c     Read from external file
c
      elseif (ipot3d.eq.1) then

         ift=45
         open (ift,file=espfile,status='old')
         read(ift,*) n
         do i=1,n
            read(ift,*) dum
         enddo
         read(ift,*) (vres(k),k=1,ng3d)

         close(ift)
C
      else

         write(*,*) "Error. Invalid potential option."
         ierr=700
         call abrt(ierr)

      endif
c-------------------------------------------------------------------
      deallocate (epsig6,epsig12)
      RETURN
c-------------------------------------------------------------------
 9999 format (/,4x,"======== Short Range Potential =========",/,
     &        /,4x,"Lennard-Jones (12-6) Potential")
      end
