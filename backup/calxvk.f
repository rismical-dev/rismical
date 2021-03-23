c-----------------------------------------------------------------
c     Calculate Xvv(k)
c
c     h_ij(k)=w_ij(k)+rho_j*h_ij(k)
c
c-----------------------------------------------------------------
      subroutine calxvk(ngrid,n2,hvk,xvk,wk2)
c
c     ngrid     ... number of grid of RDF
c     nv        ... number of site of solvent
c     hvk       ... k-space V-V total correlation function
c     xvvk      ... k-space solvent susceptibility
c     wk2       ... k-space intramolecular correlation functions
c     dens      ... number density[molecule/angstrom3]
c
      implicit real*8 (a-h,o-z)

      include "phys_const.i"
      include "solvent.i"

      dimension hvk(ngrid,n2,n2)
      dimension xvk(ngrid,n2,n2)
      dimension wk2(ngrid,n2,n2)

c-----------------------------------------------------------------
c     
      do k=1,ngrid
c     
         do i1=1,n2
            do i2=1,n2
               xvk(k,i1,i2)=(wk2(k,i1,i2)+dens(nspc(i2))*hvk(k,i1,i2))
     &                     *dens(nspc(i1))
            enddo
         enddo
c     
      enddo
c-------------------------------------------------------------------
      return
      end
