c----------------------------------------------------------------
c     Read reduced xvv(k) function from file
c----------------------------------------------------------------
      subroutine readxvk(ngrid,n2uq,xvk)
c
      implicit real*8(a-h,o-z)
      character*2 char2
      character*80 char80
c
      include "rismio.i"
      include "solvent.i"
      include "solute.i"
c
      dimension xvk(ngrid,n2uq,n2uq)
c
c----------------------------------------------------------------
      ift=45
      open (ift,file=solventxvv)
c
c     Skip header
c
      do i=1,5+nv
         read(ift,*) char2
      enddo
c
c     Read reduced xvv(k)
c
      do ig=1,ngrid
         read(ift,*) ((xvk(ig,i1,i2),i1=1,nvuq),i2=1,nvuq)
      enddo

      close(ift)
c----------------------------------------------------------------
      return
      end
