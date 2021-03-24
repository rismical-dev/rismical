c----------------------------------------------------------------
c     Read  xvv(k) function from file
c----------------------------------------------------------------
      subroutine readxvk(ngrid,n2,xvk)
c
      implicit real*8(a-h,o-z)
      character*2 char2
      character*80 char80
c
      include "rismio.i"
      include "solvent.i"
      include "solute.i"
c
      dimension xvk(ngrid,n2,n2)
c
c----------------------------------------------------------------
      ift=45
      open (ift,file=solventxvv)
c
c     Skip header
c
      read(ift,*) char2
      read(ift,*) char2
      read(ift,*) char2,ndum
      do i=1,ndum
         read(ift,*) char2
      enddo
c
c     Read xvv(k)
c
      do i2=1,n2
         do i1=1,n2
            do ig=1,ngrid
               read (ift,*) dum
               xvk(ig,i1,i2)=dum
            enddo
         enddo
      enddo

      close(ift)
c----------------------------------------------------------------
      return
      end
