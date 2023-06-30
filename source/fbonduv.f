c------------------------------------------------------------
c     Make f-bond
c------------------------------------------------------------
      subroutine fbonduv(ngrid,n1,n2,rdelta,fr,fk,alp)
c
c     ngrid     ... number of grid of RDF
c     n1,n2     ... number of site
c     rdelta    ... grid width of r-space [Angstrom]
c     fr        ... f-bond [erg]
c     fk        ... k-space f-bond [erg*Ang^3]
c     q1,q2     ... partial charge of site [e]
c     alp       ... parameter for f-bond
c
      implicit real*8 (a-h,o-z)
      
      include "phys_const.i"
      include "solvent.i"
      include "solute.i"
      
      dimension fr(ngrid,n1,n2),fk(ngrid,n1,n2)
c     
c---------------------------------------------------------
c
c     r^{-1} Asymptoticity
c
      deltak=pi/(rdelta*dble(ngrid))

      do j=1,n2
         do i=1,n1
            do k=1,ngrid
               rr=rdelta*dble(k)
               rk=deltak*dble(k)
               
               fr(k,i,j)=qu(i)*qv(j)/rr*erf(alp*rr)*fel
               fk(k,i,j)=4.d0*pi*qu(i)*qv(j)/rk**2
     &              *dexp(-(rk/(2.d0*alp))**2)*fel
            enddo
         enddo
      enddo
c
c---------------------------------------------------------
      return
      end
