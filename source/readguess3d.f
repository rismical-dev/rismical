c----------------------------------------------------------------
c     Read Guess For tr From 3D-UVDATA 
c----------------------------------------------------------------
      subroutine readguess3d(ng3d,n2uq,tr)
c
c     ir            ... file number of output (STDOUT)
c     iw            ... file number of input (STDIN)
c     ngrid3d       ... number of grid of 3D-RDF
c     nvuq          ... number of symmetry uniq site of solvent
c     tr            ... tau bond =hr-cr
c     cr            ... direct correlation function
c     ck            ... k-space direct correlation function
c
      implicit real*8 (a-h,o-z)
      complex*16 cr
      character*1 char1
      character*256 scrjob

      include "rismio.i"
      include "rismrun.i"

      dimension tr(ng3d,n2uq)
      dimension cr(ng3d,n2uq)
c----------------------------------------------------------------
c
c     Set guess file name
c
      if (len_trim(guessfile).eq.0) then
         scrjob=trim(basename)//".tuv"
      else
         scrjob=guessfile
      endif
c
c     
c
      write(*,*)
      write(*,*) "   Reading Guess from ",trim(scrjob)
c
c     Read guess t(r)
c     
      call read3dfunc(scrjob,tr,n2uq,ng3d,outtype)

c----------------------------------------------------------------
      return
      end
