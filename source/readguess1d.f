c----------------------------------------------------------------
c     Read Guess For tr From UVDATA or VVDATA
c----------------------------------------------------------------
      subroutine readguess1d(ngrid,n1,n2,tr)
c
c     ngrid         ... number of grid of RDF
c     n1            ... number of site of 1 (solute or solvent)
c     n2            ... number of site of 2 (solvent)
c     tr            ... (tau bond =hr-cr)-fr
c
      implicit real*8 (a-h,o-z)
      character*1 char1
      character*6 char6
      character*256 scrjob


      include "rismio.i"

      dimension tr(ngrid,n1,n2)

C----------------------------------------------------------------
c
c     Set guess file name
c
      if (len_trim(guessfile).eq.0) then
         scrjob=trim(basename)//".tuv"
      else
         scrjob=guessfile
      endif
c
      write(*,*)
      write(*,*) "   Reading Guess from ",trim(scrjob)
c
c     Read guess t(r)
c
      call read1dfunc(scrjob,tr,n1,n2,ngrid)

C----------------------------------------------------------------
      return
 9998 format (e16.8,1x,7e20.12)
      end

      
