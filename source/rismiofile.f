c**************************************************************
c--------------------------------------------------------------
c     Set RISM IO files
c--------------------------------------------------------------
c**************************************************************
      subroutine rismiofile
c
      implicit real*8(a-h,o-z)
      character*256 char256
      include "rismio.i"
c--------------------------------------------------------------
c
c     Get job basename
c
      call getarg(2,inpfile)
      idot = index(inpfile,".",back=.true.)
      name_len=len_trim(inpfile)
      if (idot.eq.0) then
         basename=inpfile
      elseif (idot.gt.1) then
         basename=inpfile(1:idot-1)
      else
         write(*,*) "No input file given."
         ierr=1
         call abrt(ierr)
      endif
c
c     Get home path
c
      call getenv('RISMICALHOME',homepath)
c
      write(*,*)
      write(*,*) "RISMiCal home dir:",homepath
c--------------------------------------------------------------
      return
      end
