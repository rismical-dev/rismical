c**************************************************************
c--------------------------------------------------------------
c     READ RISM INPUT
c--------------------------------------------------------------
c**************************************************************
      subroutine readinput
c
c     icl        ... closure type 0...HNC, 1...MSA, 2...KH, 3...HNC+RBC
c     iguess     ... guess type of tr 0..f-bond 1..read file
c     ngrid      ... number of grid of RDF
c     nsub       ... number of mdiis sub-space
c     dumpmax    ... maximum value of dumping parameter eta
c     dumpmin    ... manimum value of dumping parameter eta
c     dumpnume   ... a numerator of dumping function
c
      implicit real*8 (a-h,o-z)
      character*2 CHAR2
      character*3 CHAR3,CLOSURE,CHARGEUP
      character*6 CHAR6
      character*9 CHAR9

      LOGICAL READERROR
c
      include "rismio.i"
      include "solvent.i"
      include "rismrun.i"

      namelist /RISM/CLOSURE
     &              ,ITRMAX,CONV,CHARGEUP,IGUESS
     &              ,ALP1D,ALP3D,iolist,guessfile,grid,outtype
      namelist /MDIIS/nsub,dumpmax,dumpmin,dumpnume
      namelist /CHARGEUPOPT/chgstep,chgconv
C
c---------------------------------------------------------------
c     Default Value
c---------------------------------------------------------------
      icl=2
      nu=0
      nv=0
      nsub=10
      dumpmax=0.8d0
      dumpmin=0.1d0
      dumpnume=0.1d0
      CLOSURE="KH"
      CHARGEUP="ON"
      IGUESS=0
      itrmax=1000
      conv=1.d-8
      alp1d=1.5d0
      alp3d=1.0d0
      iolist=""
      guessfile=""
      grid="standard"
      outtype="ASCII"
c---------------------------------------------------------------
c     Read
c---------------------------------------------------------------
      READERROR=.FALSE.
C
c     --- Read $RISM
c
      ir=45
      open (ir,file=inpfile)
c
      rewind ir
      read(ir,RISM,end=7000)

 
      if (CLOSURE.eq."HNC") then
         icl=0
      elseif (CLOSURE.eq."MSA") then
         icl=1
      elseif (CLOSURE.eq."KH") then
         icl=2
      elseif (CLOSURE.eq."RBC") then
         icl=3
         if (iruntyp.eq.0) then
            write(*,*) "Error. CLOSURE=RBC can be used for UV system."
            readerror=.true.
            goto 8000
         endif
      else
         write (*,*) "Error. Wrong CLOSURE was required"
         readerror=.true.
         goto 8000
      endif

      call upcasex(chargeup)
      if (chargeup.eq."ON") then
         chgstep=0.1d0
         chgconv=conv*1.d+3
         chgratio=0.d0
         rewind ir
         read(ir,chargeupopt,end=8700)
 8700    continue

      elseif (chargeup.eq."OFF") then
         chgstep=0.d0
         chgratio=1.d0
         chgconv=conv
      else
         write (*,*) "Error. Wrong CHARGEUP parameter was required"
         readerror=.true.
         goto 8000
      endif

      if (ALP1D.LE.0.D0) THEN
         write (*,*) "ALP1D in $RISM must be greater than 0.d0."
         readerror=.true.
         goto 8000
      endif
      if (ALP3D.LE.0.D0) THEN
         write (*,*) "ALP3D in $RISM must be greater than 0.d0."
         readerror=.true.
         goto 8000
      endif
c
      call upcasex(grid)
      call upcasex(outtype)
c
c     --- Read $MDIIS
c
      rewind ir
      read (ir,mdiis,end=7100)
 7100 continue
c
      close (ir)
c---------------------------------------------------------------
c     Output Summary of Input
c---------------------------------------------------------------
      write (*,9997)
      write(*,9995) CLOSURE,itrmax,conv

      write(*,9992) CHARGEUP

      if (iguess.eq.0) then
         char9="F-BOND   "
      else
         char9="READ FILE"
      endif
      write(*,9989) char9
c---------------------------------------------------------------
C
      return
c
c---------------------------------------------------------------
 7000 continue
      write(*,*) "NAMELIST $RISM was not found"
      READERROR=.TRUE.
      GOTO 8000
 8000 continue
      write(*,*) "INPUT ERROR"
      ierr=800
      call abrt(ierr)
c---------------------------------------------------------------
 9989 format (  4x,"Guess of tau-bond      :",1x,a9)
 9992 format (  4x,"Use CHARGEUP Procedure :",1x,a3)
 9995 format (  4x,"Closure Equation       :",1x,a3,
     &        /,4x,"Maximum Iteration      :",i10,
     &        /,4x,"Convergence Criterion  :",1pe12.4)
 9997 format (/,4x,"======= Computational Conditions =======")
 9998 format (i4)
 9999 format (a6)
      end
