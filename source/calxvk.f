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
               xvk(k,i1,i2)=(wk2(k,i1,i2)+dens(nspc(i1))*hvk(k,i1,i2))
            enddo
         enddo
c     
      enddo
c-------------------------------------------------------------------
      return
      end
c-----------------------------------------------------------------
c     Calculate reduced Xvv(k)
c
c     h_ij(k)=w_ij(k)+rho_j*h_ij(k)
c
c-----------------------------------------------------------------
      subroutine calredxvk(ngrid,n2,hvk,xvk,wk2)
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

      dimension redwk2(n2,n2)
      dimension redhvk(n2,n2)
c-----------------------------------------------------------------
c     
      do k=1,ngrid
c     
c     reducing wk2
c
         redwk2=0.d0
         do i1=1,n2
         do i2=1,n2
            i1uq=abs(iuniq(i1))
            i2uq=iuniq(i2)
            if (i2uq.gt.0) then
               redwk2(i1uq,i2uq)=redwk2(i1uq,i2uq)+wk2(k,i1,i2)
            endif
         enddo
         enddo
c     
c     reducing hvk
c
         redhvk=0.d0
         do i1=1,n2
         do i2=1,n2
            i1uq=iuniq(i1)
            i2uq=iuniq(i2)
            if ((i1uq.gt.0).and.(i2uq.gt.0)) then
               redhvk(i1uq,i2uq)=hvk(k,i1,i2)
            endif
         enddo
         enddo
c
c     gettin reduced xvv(k)
c
         do i1=1,nvuq
         do i2=1,nvuq
            xvk(k,i1,i2)=(redwk2(i1,i2)+densuq(i1)*redhvk(i1,i2))
         enddo
         enddo
c     
      enddo
c-------------------------------------------------------------------
      return
      end
