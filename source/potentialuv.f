c**************************************************************
c----------------------------------------------------------------
c     Make intraction potential for uv system : reduced solvent
c----------------------------------------------------------------
      subroutine potentialuv(ngrid,n1,n2uq,rdelta,ures,urlj)
c     
c     ngrid         ... number of grid of RDF
c     n1            ... number of site of 1 (solute or solvent)
c     n2uq          ... number of site of 2 (solvent)
c     rdelta        ... grid width of r-space [Angstrom]
c     ures          ... electro static potential [erg]
c     urlj          ... LJ potential [erg]
c     siglj1,siglj2 ... sigma of LJ parameter [Angstrom]
c     epslj1,epslj2 ... epsilon of LJ parameter [erg]
c     q1,q2         ... partial charge of site [e]
c     epsig12,epsig6 ... combined parameter

      implicit real*8 (a-h,o-z)
      character*8 char8

      include "phys_const.i"
      include "solvent.i"
      include "solute.i"
      
      dimension ures(ngrid,n1,n2uq),urlj(ngrid,n1,n2uq)
      dimension epsig12(n1,n2uq),epsig6(n1,n2uq)
      dimension ftfunc(ngrid)
c
c-----------------------------------------------------------------
c
c     --- LJ Parameter Settings
c
      call rsmljcomb(6,n1,n2uq
     &     ,siglju,epslju,sigljvuq,epsljvuq,epsig6)
      call rsmljcomb(12,n1,n2uq
     &     ,siglju,epslju,sigljvuq,epsljvuq,epsig12)
c     
c     --- Initialize
c     
      do j=1,n2uq
         do i=1,n1
            do k=1,ngrid
               ures(k,i,j)=0.d0
               urlj(k,i,j)=0.d0
            enddo
         enddo
      enddo
c     
c     --- Make Potential
c     
      do j=1,n2uq
         do i=1,n1
            do k=1,ngrid
               
               rr=rdelta*dble(k)
c     
c     --- LJ
c     
               rr6=epsig6(i,j)/rr**6
               rr12=epsig12(i,j)/rr**12
               urlj(k,i,j)=4.d0*(rr12-rr6)    ![J/mol]
c     
c     --- Electro Static
c     
               ures(k,i,j)=qu(i)*q2uq(j)/rr*fel ![e**2/Ang --> J/mol]
               
            enddo
         enddo
      enddo
c----------------------------------------------------------------
      return
      end
c**************************************************************
