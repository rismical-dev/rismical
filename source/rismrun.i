c     --------------------------------------------
c     Runcontrol parameters
c     --------------------------------------------

      integer icl,iguess,ngrid,ngrid3d
      real*8 rdelta,rdelta3d,conv,rcore,alp1d,alp3d
      character*24 grid
     
      integer nsub
      real*8 dumpmax,dumpmin,dumpnume

      real*8 chgstep,chgratio,chgconv

      common /rismrun_i/icl,iguess
     &               ,ngrid,ngrid3d
     &               ,itrmax
      common /rismrun_r/rdelta,rdelta3d,conv,rcore
     &               ,alp1d,alp3d
c
      common /rismmdiis_i/nsub
      common /rismmdiis_r/dumpmax,dumpmin,dumpnume
c
      common /rismchgup/chgstep,chgratio,chgconv
      common /rismgrid/grid
c
