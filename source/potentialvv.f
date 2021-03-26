c**************************************************************
c----------------------------------------------------------------
c     Make intraction potential for vv system
c----------------------------------------------------------------
      subroutine potentialvv(ngrid,n2,rdelta,ures,urlj)
c     
c     ir            ... file number of output (STDOUT)
c     iw            ... file number of input (STDIN)
c     ilj           ... Short range potential type 0..LJ 1..X6 2..RF
c     iuv           ... OZ type 0...vv, 1...uv
c     ngrid         ... number of grid of RDF
c     n1            ... number of site of 1 (solute or solvent)
c     n2            ... number of site of 2 (solvent)
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
      
      dimension ures(ngrid,n2,n2),urlj(ngrid,n2,n2)
      dimension epsig12(n2,n2),epsig6(n2,n2)
      dimension dumn2(n2),ftfunc(ngrid)
c
c-----------------------------------------------------------------
c
c     --- LJ Parameter Settings
c
      call rsmljcomb(6,n2,n2,
     &     sigljv,epsljv,sigljv,epsljv,epsig6)
      call rsmljcomb(12,n2,n2,
     &     sigljv,epsljv,sigljv,epsljv,epsig12)
c     
c     --- Initialize
c     
      call vclr_mp(urlj,1,ngrid*n2*n2)
      call vclr_mp(ures,1,ngrid*n2*n2)
c     
c     --- Make Potential
c     
      do j=1,n2
         do i=1,n2
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
               ures(k,i,j)=qv(i)*qv(j)/rr*fel ![e**2/Ang --> J/mol]
               
            enddo
         enddo
      enddo
c----------------------------------------------------------------
      return
      end
c**************************************************************
