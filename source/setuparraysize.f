c---------------------------------------------------------
c     Setup Array Size
c---------------------------------------------------------
      subroutine setuparraysize

      implicit real*8(a-h,o-z)

      include "phys_const.i"
      include "solute.i"
      include "solvent.i"
      include "rismrun.i"
      include "rismio.i"

      dimension iuniqlabel(maxslv)

c---------------------------------------------------------
c
c     --- Count Number of Xvv
c
      n=ngrid3d/2 
      nxvv=n+n*(n-1)+n*(n-1)*(n-2)/6
c
c     --- Make Symmetry Uniq Solvent Point Charge
c
      do i=1,nv
         if (iuniq(i).gt.0) then
            q2uq(iuniq(i))=qv(i)
         endif
      enddo
c
c     --- Set symmetry multiplicity of solvent site
c
      do i=1,maxspc
         nmulsite(i)=0
      enddo
      do i=1,nv
         j=abs(iuniq(i))
         nmulsite(j)=nmulsite(j)+1
      enddo
      do i=1,nv
         j=abs(iuniq(i))
         if (iuniq(i).gt.0) then
            densuq(j)=dens(nspc(i))*nmulsite(j)
         endif
      enddo
c---------------------------------------------------------
C
      return
      end
