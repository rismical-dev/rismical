c----------------------------------------------------------------
c     Calculate and Output free energy gradient 
c----------------------------------------------------------------
      subroutine write3dgrad(namef,gr3d,rdelta3d
     &     ,ngrid3d,nvuqx,listcore)
      implicit real*8(a-h,o-z)
      character*256 namef

      include "phys_const.i"
      include "solute.i"
      include "solvent.i"
c      
      real*8 ,allocatable :: epsig6(:,:)
      real*8 ,allocatable :: epsig12(:,:)
c      
      dimension gr3d(ngrid3d**3,nvuqx)
      dimension delj(3),dees(3)
      dimension listcore(ngrid3d**3)

      allocate (epsig6(nu,nvuq),epsig12(nu,nvuq))
c----------------------------------------------------------------
c
c     --- LJ Parameter Settings
c
      call rsmljcomb(6,nu,nvuq
     &     ,siglju,epslju,sigljvuq,epsljvuq,epsig6)
      call rsmljcomb(12,nu,nvuq
     &     ,siglju,epslju,sigljvuq,epsljvuq,epsig12)
c
c     --- Calculate and output free energy derivative
c
      ift=45
      open (ift,file=namef)

      rd3=rdelta3d**3
      k0=ngrid3d/2+1
      
      do iu=1,nu

         do i=1,3
            delj(i)=0.d0
            dees(i)=0.d0
         enddo
            
!$omp parallel do collapse(3) default(shared)
!$omp& private(kz,ky,kx,k,rx,ry,rz,rr,rr2,rr3,rr6,rr12)
!$omp& private(jv,fac,dulj,dues)
!$omp& reduction(+:delj,dees)
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d

            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2

            rx=rdelta3d*dble(kx-k0)-xyzu(1,iu)
            ry=rdelta3d*dble(ky-k0)-xyzu(2,iu)
            rz=rdelta3d*dble(kz-k0)-xyzu(3,iu)
            rr2=rx**2+ry**2+rz**2
            rr=sqrt(rr2)
            rr3=rr2*rr
            rr6=rr2**3
            rr12=rr6*rr6

            if (listcore(k).eq.1) then

            do jv=1,nvuq

               fac=densuq(jv)*rd3*gr3d(k,jv) ! Unitless
c
c     LJ
c               
               dulj=24.d0*(2.d0*epsig12(iu,jv)/rr12
     &                         -epsig6(iu,jv) /rr6 )/rr2*fac
               
               delj(1)=delj(1)+dulj*rx ! J/mol/Angs
               delj(2)=delj(2)+dulj*ry
               delj(3)=delj(3)+dulj*rz
c
c     ES
c
               dues=qu(iu)*q2uq(jv)/rr3*hart2jmol/angtobohr*fac
               dees(1)=dees(1)+dues*rx  ! J/mol/Angs
               dees(2)=dees(2)+dues*ry
               dees(3)=dees(3)+dues*rz
               
            enddo   ! of nvuq

            endif
            
         enddo      ! of kx
         enddo      ! of ky
         enddo      ! of kz
!$omp end parallel do
c
c     Output derivative  [J/mol/Angs]
c     NOTE: this is not a force.
c
         write (ift,'(6e16.8)') (delj(i),i=1,3),(dees(i),i=1,3)

      enddo         ! of nu

      close(ift)
c----------------------------------------------------------------
      deallocate (epsig6,epsig12)
      return
      end
