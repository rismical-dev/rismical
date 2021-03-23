c----------------------------------------------------------------
c     Read  hv(k) function from file
c----------------------------------------------------------------
      subroutine readhvk(ngrid,n2,hvk)
c
      implicit real*8(a-h,o-z)
      character*2 char2
      character*80 char80
c
      include "rismio.i"
      include "solvent.i"
      include "solute.i"
c
      dimension hvk(ngrid,n2,n2)
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
c     Read hvk
c
      do i2=1,n2
         do i1=1,n2
            do ig=1,ngrid
               read (ift,*) dum
               hvk(ig,i1,i2)=dum
            enddo
         enddo
      enddo

      close(ift)
c----------------------------------------------------------------
      return
      end
