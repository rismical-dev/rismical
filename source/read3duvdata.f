c**************************************************************
c--------------------------------------------------------------
c    Read solute and solvent parameters for 3D-RISM
c--------------------------------------------------------------
c**************************************************************
      subroutine read3duvdata
c
      implicit real*8(a-h,o-z)
      character*2 char2a,char2b,char2c,char2d
      character*6 char6a,char6b,chardum
      character*6 esptype
      character*256 solute,solutexyz,soluteesp,solutelj
      character*256 soluteepc
      character*256 solvent,solute_up,ljparam

      include "phys_const.i"
      include "solvent.i"
      include "solute.i"
      include "rismrun.i"
      include "rismio.i"
c
      namelist /rismsolution/solute,solutexyz,soluteesp,solutelj
     &     ,soluteepc,solvent,ljparam,esptype
      namelist /GRID3D/ngrid3d,rdelta3d
c--------------------------------------------------------------
c     
c     Read solvent grid data
c
c
c     --- Read $GRID3D 
c
      if (grid.eq."USER") then

         ngrid3d=128
         rdelta3d=0.5d0

         ir=45
         open (ir,file=inpfile)
         rewind ir
         read (ir,grid3d,end=900)
 900     continue
         close(ir)

      elseif (grid.eq."FINE") then

         ngrid3d=256
         rdelta3d=0.25d0

      elseif (grid.eq."LFINE") then

         ngrid3d=512
         rdelta3d=0.25d0

      elseif (grid.eq."STANDARD") then

         ngrid3d=128
         rdelta3d=0.5d0

      elseif (grid.eq."LSTANDARD") then

         ngrid3d=256
         rdelta3d=0.5d0

      elseif (grid.eq."TEST") then

         ngrid3d=64
         rdelta3d=1.d0

      elseif (grid.eq."LTEST") then

         ngrid3d=128
         rdelta3d=1.d0

      endif
c
      write(*,*) "     --------------------------------------"
      write(*,'(A19,A24)')   "Grid preset       :",grid
      write(*,'(A19,i12)')   "Number of 3D-Grid :",ngrid3d
      write(*,'(A19,f12.4,A5)') "3D-Grid width    :",rdelta3d," [A]"
      write(*,*) "     --------------------------------------"
c
c     --------------------------------------------------------
c
c     Read option for parameters
c
      ift=45
      open (ift,file=inpfile,status='old')
c
c     Set defaults
c
      solute="udata"
      ljparam="mm2.prm"
      esptype="PC"
      ipot3d=0
c
c     Read namelist rismsolvent
c
      rewind ift
      read (ift,rismsolution,end=1000)
 1000 continue
      solute_up=solute
      call upcasex(solute_up)
      call upcasex(esptype)
      if (esptype.eq."MAP") then
         ipot3d=1
         espfile=soluteesp
      endif
c
      if (len_trim(solvent).eq.0) then
         write(*,*) "Error. solvent in $rismsolution is empty."
         write(*,*) "solvent Xvv file must be given."
         ierr=30
         call abrt(ierr)
      endif
      solventxvv=solvent
c
c     -----------------------
c     Get Solute Params
c     -----------------------
c
c     Read solute parameter from $UDATA
c
      if (trim(solute_up).eq."UDATA") then

         char6a="$UDATA"
         rewind ift
 2000    read(ift,*,end=9999) char6b
         call upcasex(char6b)
         if (char6b.ne.char6a) goto 2000
c
         read(ift,*) nu
         do iu=1,nu
            read(ift,*) nsiteu(iu)
     &           ,siglju(iu),epslju(iu),qu(iu)
     &           ,xyzu(1,iu),xyzu(2,iu),xyzu(3,iu)
         enddo
c     
c     Read solute parameter from external files
c
      else
c
c     Set default if not given
c
         if (len_trim(soluteepc).eq.0) soluteepc=trim(solute)//".epc"
         if (len_trim(solutexyz).eq.0) solutexyz=trim(solute)//".xyz"
         if (len_trim(solutelj ).eq.0) solutelj =trim(solute)//".lj"
c
c     Read xyz file to get coordinate
c
         ift2=46
         open(ift2,file=solutexyz,status='old')
         read(ift2,*) nu
         read(ift2,*) chardum
         do i=1,nu
            read (ift2,*) nsiteu(i),xyzu(1,i),xyzu(2,i),xyzu(3,i)
         enddo
         close(ift2)
c
c     Read EPC file to get RESP point charge
c
         call readresp(soluteepc,qu,maxslu)
c
c     Get LJ parameter 
c
         solute_up=solutelj
         call upcasex(solute_up)
         if (solute_up.eq."BUILTIN") then

            call readbuiltinsoluteparam(ljparam)

         elseif (solute_up.eq."UDATA") then
            
            char6a="$UDATA"
            rewind ift
 3000       read(ift,*,end=9999) char6b
            call upcasex(char6b)
            if (char6b.ne.char6a) goto 3000
c
            read(ift,*) n
            do iu=1,nu
               read(ift,*) dum,siglju(iu),epslju(iu)
            enddo

         else

            ift2=46
            open(ift2,file=solutelj,status='old')
            read(ift2,*) n
            if (n.ne.nu) goto 9998
            do i=1,nu
               read (ift2,*) siglju(i),epslju(i)
            enddo
            close(ift2)
            
         endif

      endif

      close(ift)
c
c     Print Solute parameters
c
      write(*,*) "-----------------------------------------------"
      write(*,*) "            Solute parameters"
      write(*,*) "-----------------------------------------------"
      write(*,8000) nu
 8000 format ("Number of solute site:",i4)
c
      write(*,9801)
      do i=1,nu
         write(*,9800) nsiteu(i),siglju(i),epslju(i),qu(i)
     &           ,xyzu(1,i),xyzu(2,i),xyzu(3,i)
      enddo
c
c     -----------------------
c     Get Solvent Params
c     -----------------------
c
      ift=45
      open(ift,file=solvent,status='old')
      read(ift,*) char2a
      read(ift,*) char2a,nv,nvuq,ngrid,rdelta,deltak
      read(ift,*) char2a,ndum,temp,xt
      read(ift,*) char2a
      do i=1,nv
         read(ift,*) char2a,nsitev(i),nspc(i),iuniq(i)
     &        ,sigljv(i),epsljv(i)
     &        ,qv(i),xyzv(1,i),xyzv(2,i),xyzv(3,i),densv
         dens(nspc(i))=densv*avognum*1.D-27
      enddo
      close(ift)
c
c     Set reduced solvent parameters
c
      nmulsite=0
      densuq=0
      do i=1,nv
         iuq=iuniq(i)
         if (iuq.gt.0) then
            q2uq(iuq)=qv(i)
            epsljvuq(iuq)=epsljv(i)
            sigljvuq(iuq)=sigljv(i)
         endif
         nmulsite(abs(iuq))=nmulsite(abs(iuq))+1
         densuq(abs(iuq))=densuq(abs(iuq))+dens(nspc(i))
      enddo
c
c     Set inverse temperature
c      
      beta=1.d0/(gasconst*temp)  ![mol/J]
c
c     Print Solvent parameters
c
      write(*,*) "-----------------------------------------------"
      write(*,*) "            Solvent parameters"
      write(*,*) "-----------------------------------------------"
      write(*,8001) nv
 8001 format ("Number of solvnt site:",i4)
c
      write(*,9802)
      do i=1,nv
         densv=dens(nspc(i))/(avognum*1.D-27)
         write(*,9803) nsitev(i),nspc(i),iuniq(i)
     &           ,sigljv(i),epsljv(i),qv(i)
     &           ,xyzv(1,i),xyzv(2,i),xyzv(3,i),densv
      enddo
      write(*,9804) temp
      write(*,*) "Solvent xvv file:",trim(solventxvv)
c--------------------------------------------------------------
      return
 9998 write(*,*) "Error. Solute data missmatch"
      ierr=23
      call abrt(ierr)
 9999 write(*,*) "Error. $UDATA is not given."
      write(*,*) "solute=UDATA requires $UDATA."
      ierr=22
      call abrt(ierr)
 9800 format (A4,1x,6f12.5)
 9801 format ("ATOM"," sig[Angs]  "," eps[J/mol] "," charge[e]  "
     &     ,"  ---X---   ","  ---Y---   ","  ---Z---   ")
 9802 format ("ATOM"," SPC"," sig[Angs]  "," eps[J/mol] "," charge[e]  "
     &     ,"  ---X---   ","  ---Y---   ","  ---Z---   "," density[M] ")
 9803 format (A4,1x,i3,1x,i3,7f12.5)
 9804 format ("Temperature :",f12.5,"[K]")
      end
