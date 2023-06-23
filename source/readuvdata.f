c**************************************************************
c--------------------------------------------------------------
c    Read solute and solvent parameters
c--------------------------------------------------------------
c**************************************************************
      subroutine readuvdata
c
      implicit real*8(a-h,o-z)
      character*2 char2a,char2b,char2c,char2d
      character*6 char6a,char6b,chardum
      character*256 solute,solutexyz,soluteesp,solutelj
      character*256 solvent,solute_up,ljparam

      include "phys_const.i"
      include "solvent.i"
      include "solute.i"
      include "rismrun.i"
      include "rismio.i"
c
      namelist /rismsolution/solute,solutexyz,soluteesp,solutelj
     &     ,solvent,ljparam
c--------------------------------------------------------------
c
      ift=45
      open (ift,file=inpfile,status='old')
c
c     Set defaults
c
      solute="udata"
      ljparam="mm2.prm"
c
c     Read namelist rismsolvent
c
      rewind ift
      read (ift,rismsolution,end=1000)
 1000 continue
      solute_up=solute
      call upcasex(solute_up)
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
         if (len_trim(soluteesp).eq.0) soluteesp=trim(solute)//".esp"
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
c     Read EPS file to get point charge
c
         call readresp(soluteesp,qu,maxslu)
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
            read (ift2,*) siglju(iu),epslju(iu)
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
      read(ift,*) char2a,char2b,char2c,char2d,nv,ndum,ngrid,rdelta
      read(ift,*) char2a
      read(ift,*) char2a
      read(ift,*) char2a,ndum,ndum2,nvuq,temp,xt
      read(ift,*) char2a
      do i=1,nv
         read(ift,*) char2a,nsitev(i),nspc(i),iuniq(i)
     &            ,sigljv(i),epsljv(i)
     &            ,qv(i),xyzv(1,i),xyzv(2,i),xyzv(3,i),densv
         dens(nspc(i))=densv*avognum*1.D-27
      enddo
      close(ift)
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
     &        ,sigljv(i),epsljv(i),qv(i)
     &        ,xyzv(1,i),xyzv(2,i),xyzv(3,i),densv
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
c--------------------------------------------------------------
c     Read built in LJ parameter file
c--------------------------------------------------------------
      subroutine readbuiltinsoluteparam(ljparam)

      implicit real*8(a-h,o-z)
      character*4 natom,nameup
      character*10 char10,rtype,rsize
      character*256 ljparam
      parameter (maxparam=100)
c
      include "solute.i"
      include "rismio.i"
c
      dimension natom(maxparam)
      dimension sigp(maxparam),epsp(maxparam)
c--------------------------------------------------------------
c
c     Read parameter file
c
      ift=47
      open (ift,file=trim(homepath)//
     &     "/params/"//trim(ljparam),status='old')
      read(ift,*) char10
      read(ift,*) char10,rtype
      read(ift,*) char10,rsize
      call upcasex(rtype)
      call upcasex(rsize)
      read(ift,*) char10,nparam
      read(ift,*) char10
      do i=1,nparam
         read(ift,*) sigp(i),epsp(i),natom(i)
         call upcasex(natom(i))
      enddo
c
c     Mod params to DIAMETER-SIGMA
c
      if (trim(rtype).eq."R-MIN".and.trim(rsize).eq."RADIUS") then
         do i=1,nparam
            sigp(i)=sigp(i)*1.781797d0
         enddo
      endif
c
      if (trim(rtype).eq."SIGMA".and.trim(rsize).eq."RADIUS") then
         do i=1,nparam
            sigp(i)=sigp(i)*2.d0
         enddo
      endif
c
      if (trim(rtype).eq."SIGMA".and.trim(rsize).eq."DIAMETER") then
c        Nothing to do
      endif
c
      close(ift)
c
c     Assign parameter
c
      do i=1,nu
         nameup=nsiteu(i)
         call upcasex(nameup)
         do j=1,nparam
            if (trim(nameup).eq.trim(natom(j))) then
               epslju(i)=epsp(j)
               siglju(i)=sigp(j)
            endif
         enddo
         
      enddo
c--------------------------------------------------------------
      return 
      end
