c
c     Read solvent data from input file
c
      subroutine readvvdata
c
      implicit real*8(a-h,o-z)
      character*6 char6a,char6b
      character*256 char256
      character*80 solvent
      logical userdef

      include "phys_const.i"
      include "solvent.i"
      include "rismio.i"
      include "rismrun.i"
c
      dimension solvent(maxspc)
c
      namelist /rismsolvent/solvent,temp,dens,numspc
      namelist /GRID1D/ngrid,rdelta
c----------------------------------------------------------------------
c     
c     Read solvent grid data
c
c
c     --- Read $GRID1D 
c
      if (grid.eq."USER") then
         ngrid=4096
         rdelta=0.05d0

         ir=45
         open (ir,file=inpfile)
         rewind ir
         read (ir,grid1d,end=900)
 900     continue
         close(ir)

      elseif (grid.eq."FINE") then
         ngrid=8192
         rdelta=0.025d0

      elseif (grid.eq."LFINE") then

         ngrid=16384
         rdelta=0.025d0

      elseif (grid.eq."STANDARD") then
         ngrid=4096
         rdelta=0.05d0

      elseif (grid.eq."LSTANDARD") then
         ngrid=8192
         rdelta=0.05d0

      elseif (grid.eq."TEST") then
         ngrid=2048
         rdelta=0.1d0

      elseif (grid.eq."LTEST") then
         ngrid=4096
         rdelta=0.1d0

      endif

c
      write(*,*) "   --------------------------------------"
      write(*,'(4x,A19,A24)')   "Grid preset       :",grid
      write(*,'(4x,A19,i12)')   "Number of 1D-Grid :",ngrid
      write(*,'(4x,A19,f12.4,A5)') "Grid width       :",rdelta," [A]"
      write(*,*) "   --------------------------------------"
c
      ift=45
      open(ift,file=inpfile,status='old')
c
c     Set defaults
c
      numspc=1
      solvent="TIP3P"
      temp=298.15d0
      dens(1)=55.d0
c
c     Read namelist rismsolvent
c
      rewind ift
      read (ift,rismsolvent,end=1000)
 1000 continue
      userdef=.false.
      do i=1,numspc
         call upcasex(solvent(i))
         if (trim(solvent(i)).eq."USER") userdef=.true.
      enddo
c
c     Read user define solvent $VDATA
c     
      if (userdef) then
         char6a="$VDATA"
         rewind ift
 7000    read(ift,*,end=9901) char6b
         if (char6b.ne.char6a) goto 7000
c
         iv=0
         do j=1,numspc
            read(ift,*) n,solvent(j)
            do i=1,n
               iv=iv+1
               read(ift,*) nsitev(iv)
     &              ,sigljv(iv),epsljv(iv),qv(iv)
     &              ,xyzv(1,iv),xyzv(2,iv),xyzv(3,iv)
               nspc(iv)=j
            enddo
         enddo
         nv=iv
c
c     Read preset solvent parameter
c      
      else
         iv=0
         do j=1,numspc

            ift2=46
            open (ift2,file=trim(homepath)//"/params/"
     &           //trim(solvent(j))//".txt",err=9900)
            char6a="$VDATA"
 7100       read(ift2,*,end=9901) char6b
            if (char6b.ne.char6a) goto 7100
c
            read(ift2,*) n
            do i=1,n
               iv=iv+1
               read(ift2,*) nsitev(iv)
     &              ,sigljv(iv),epsljv(iv),qv(iv)
     &              ,xyzv(1,iv),xyzv(2,iv),xyzv(3,iv)
               nspc(iv)=j
            enddo

            close(ift2)

         enddo

         nv=iv

      endif
c
      close(ift)
c
      write(*,*) "   -----------------------------------------------"
      write(*,*) "                Solvent parameters"
      write(*,*) "   -----------------------------------------------"
      write(*,8000) numspc,nv
 8000 format (4x,"Solvent species:",i4,4x
     &       ,"Total number of solvent site:",i4)
      write(*,*) "     Density[M]  : Solvent"
      do i=1,numspc
         write(*,'(4x,f12.5,4x,A80)') dens(i),solvent(i)
      enddo
      write(*,9802) temp
c
c     Set inverse temperature
c      
      beta=1.d0/(gasconst*temp)  ![mol/J]
c
c     Print solvent parameters
c
      write(*,9801)
      do i=1,nv
         write(*,9800) nsitev(i),nspc(i),sigljv(i),epsljv(i),qv(i)
     &           ,xyzv(1,i),xyzv(2,i),xyzv(3,i),dens(nspc(i))
      enddo
c
c     Convert densty [M] -> [/Ang^3]
c
      do i=1,numspc
         dens(i)=dens(i)*avognum*1.D-27
      enddo
c----------------------------------------------------------------------
      return
c
c     errors
c
 9900 write(*,*) "Invalid solvent name.",solvent
      ierr=2
      call abrt(ierr)
 9901 write(*,*) "$VDATA is not given."
      write(*,*) "solvent=USER requires $VDATA."
      ierr=2
      call abrt(ierr)
 9800 format (4x,A4,1x,i3,7f12.5)
 9801 format (4x,"ATOM"," SPC"," sig[Angs]  "
     &     ," eps[J/mol] "," charge[e]  "
     &     ,"  ---X---   ","  ---Y---   ","  ---Z---   "," density[M] ")
 9802 format (4x,"Temperature :",f12.5," [K]")
      end
