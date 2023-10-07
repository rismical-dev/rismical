c**************************************************************
c--------------------------------------------------------------
c     Set RISM IO files
c--------------------------------------------------------------
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
      if (len(trim(homepath)).eq.0) then
         write(*,*) "Error. Please set environment variable "
     &        ,"RISMICALHOME to root directory of RISMical "
     &        ,"package."
         write(*,*) "ie: export RISMICALHOME=/foo/bar/rismical"
         call abrt
      endif
c
      write(*,*)
      write(*,*) "    RISMiCal home dir:",trim(homepath)
c--------------------------------------------------------------
      return
      end
c**************************************************************
c----------------------------------------------------------------
c     Read 1D function from file
c----------------------------------------------------------------
      subroutine read1dfunc(namef,func1d,n1,n2,ngrid)
c
      implicit real*8(a-h,o-z)
      character*256 namef
      character*2 char2
c
      include "phys_const.i"

      dimension func1d(ngrid,n1,n2)
c----------------------------------------------------------------
      ift=45
      nremark=0
      open (ift,file=namef)
c
c     read header
      read (ift,*) char2,n1x,n2x,ngridx,rdeltax
      read (ift,*) char2
      read (ift,*) char2,nremark
c
c     check consistency
      if (n1.ne.n1x .or. n2.ne.n2x
     &     .or. ngrid.ne.ngridx ) then
         write(*,*) "Error."
         write(*,*) "Reading 1d function inconsistency."
         ierr=567
         call abrt(ierr)
      endif
c
c     skip remark lines
      do i=1,nremark
         read (ift,*) char2
      enddo
c
c     read function
      do i2=1,n2
      do i1=1,n1
         do ig=1,ngrid
            read (ift,*) func1d(ig,i1,i2)
         enddo
      enddo
      enddo

      close(ift)
c----------------------------------------------------------------
      return
      end

c**************************************************************
c----------------------------------------------------------------
c     Write 1D function to file
c----------------------------------------------------------------
      subroutine write1dfunc(namef,func1d,rdelta,n1,n2,ngrid
     &                      ,nremark,char80,char802)
c
      implicit real*8(a-h,o-z)
      character*256 namef
      character*2 char2
      character*80 char80,char802
c
      include "phys_const.i"

      dimension func1d(ngrid,n1,n2)
      dimension char802(nremark)
c
c----------------------------------------------------------------
      ift=45
      open (ift,file=namef)
      write(ift,9990) n1,n2,ngrid,rdelta
      write(ift,9991) char80
      write(ift,9992) nremark
      do i=1,nremark
         write (ift,9993) adjustl(char802(i))
      enddo

      do i2=1,n2
      do i1=1,n1
         do ig=1,ngrid
            write (ift,9995) func1d(ig,i1,i2)
         enddo
         write(ift,*)
         write(ift,*)
      enddo
      enddo
c----------------------------------------------------------------
      return
 9990 format("## 1D Function :",3i8,f16.8)
 9991 format("##  ",a80)
 9992 format("##  ",i4)
 9993 format("##  ",a80)
 9994 format(e16.8e3)
 9995 format(e16.8e3,2x,e16.8e3)
 9996 format(e16.8e3,2x,16e16.8e3)
      end
c**************************************************************
c----------------------------------------------------------------
c     Write 1D function to file with solvent parameter
c----------------------------------------------------------------
      subroutine writehvvfunc(namef,hvk,rdelta,n2,ngrid,char80)
c
      implicit real*8(a-h,o-z)
      character*256 namef
      character*2 char2
      character*80 char80
c
      include "phys_const.i"
      include "solvent.i"
c
      dimension hvk(ngrid,n2,n2)
c
c----------------------------------------------------------------
      ift=45
      nremark=0
      open (ift,file=namef)
c
c     --- write hvk
c
      write(ift,9990) n2,n2,ngrid,rdelta
      write(ift,9991) char80
      write(ift,9992) 3+n2
c
c     --- write solvent parameters
c
      write(ift,'(A3,i4,1x,i4,1x,i4,1x,f16.8,e16.8)') 
     &             "## ",numspc,nv,nvuq,temp,xt
      write(ift,9801)
      do i=1,n2
         densm=dens(nspc(i))/(avognum*1.D-27)    ! to M
         write(ift,9800) nsitev(i),nspc(i),iuniq(i),sigljv(i),epsljv(i)
     &           ,qv(i),xyzv(1,i),xyzv(2,i),xyzv(3,i),densm
      enddo
      write(ift,'(A2)') "##"
 9801 format ("## ATOM"," SPC"," SYM"
     &     ," sig[Angs]  "," eps[J/mol] "," charge[e]  "
     &     ,"  ---X---   ","  ---Y---   ","  ---Z---   "," density[M] ")
 9800 format ("## ",A4,1x,i3,1x,i3,7f12.5)
c
      do i2=1,n2
      do i1=1,n2
         do ig=1,ngrid
            write (ift,9995) hvk(ig,i1,i2)
         enddo
      enddo
      enddo

      close(ift)
c----------------------------------------------------------------
      return
 9990 format("## 1D Function :",3i8,f16.8)
 9991 format("##  ",a80)
 9992 format("##  ",i4)
 9993 format("##  ",a80)
 9994 format(e16.8e3)
 9995 format(e16.8e3,2x,e16.8e3)
 9996 format("##  ",i5)
 9997 format("##  ",10i5)
      end
c----------------------------------------------------------------
c     Write Reduced Solvent Suseptibility Function
c----------------------------------------------------------------
      subroutine writexvvfunc(namef,xvk,rdelta,n2,ngrid)
c
      implicit real*8(a-h,o-z)
      character*256 namef
      character*2 char2
      character*80 char80
c
      include "phys_const.i"
      include "solvent.i"
c
      dimension xvk(ngrid,n2,n2)
c
c----------------------------------------------------------------
      char80="reduced xvv(k) data with solvent parameters"

      ift=45
      nremark=0
      open (ift,file=namef)
c
c     --- write xvk
c
      deltak=pi/(dble(ngrid)*rdelta)
      write(ift,9991) char80
c
c     --- array size and grid parameters
c
      write(ift,9990) nv,nvuq,ngrid,rdelta,deltak
c
c     --- solvent parameters
c
      write(ift,'(A3,i4,1x,f16.8,e16.8)') 
     &             "## ",numspc,temp,xt
      write(ift,9801)
      do i=1,n2
         densm=dens(nspc(i))/(avognum*1.D-27)    ! to M
         write(ift,9800) nsitev(i),nspc(i),iuniq(i),sigljv(i),epsljv(i)
     &           ,qv(i),xyzv(1,i),xyzv(2,i),xyzv(3,i),densm
      enddo
      write(ift,'(A2)') "##"
 9801 format ("## ATOM"," SPC"," SYM"
     &     ," sig[Angs]  "," eps[J/mol] "," charge[e]  "
     &     ,"  ---X---   ","  ---Y---   ","  ---Z---   "," density[M] ")
 9800 format ("## ",A4,1x,i3,1x,i3,7f12.5)
c
c     --- write xvv(k)
c
      do ig=1,ngrid
         write (ift,9995) ((xvk(ig,i1,i2),i1=1,nvuq),i2=1,nvuq)
      enddo

      close(ift)
c----------------------------------------------------------------
      return
 9990 format("##  ",3i8,2f16.8)
 9991 format("##  ",a80)
 9992 format("##  ",i4)
 9993 format("##  ",a80)
 9994 format(e16.8e3)
 9995 format(10(1x,g22.15e3))
 9996 format("##  ",i5)
 9997 format("##  ",10i5)
      end
c**************************************************************
c----------------------------------------------------------------
c     Read 3D function from file
c----------------------------------------------------------------
      subroutine read3dfunc(namef,func3d,nvuq,ng3d,outtype)
c
      implicit real*8(a-h,o-z)
      character*256 namef,outtype,char256
      character*2 char2
c
      dimension func3d(ng3d,nvuq)
c
c----------------------------------------------------------------
      if (trim(outtype).eq."ASCII") then
         ift=45
         open (ift,file=namef)
c
c     read header
 100     continue
         read(ift,*) char256
         if (char256(1:2).eq."##") goto 100
c
c     read solvent site, grid
         read(char256,*) nvuqx,ngrid3dx,ngrid3dy,ngrid3dz,
     &        rnx,rny,rnz,shiftx,shifty,shiftz
         ng3dx=ngrid3dx*ngrid3dy*ngrid3dz
c
c     check consistency
         if (nvuq.ne.nvuqx .or. ng3d.ne.ng3dx) then
            write(*,*) "Error."
            write(*,*) "Reading 3d function inconsistency."
            ierr=568
            call abrt(ierr)
         endif
c
c     read function
         do j=1,nvuq
            do k=1,ng3d
               read(ift,*)  func3d(k,j)
            enddo
         enddo
         
         close(ift)
c
c
      elseif (trim(outtype(1:3)).eq."BIN") then

         ift=45
         open (ift,file=namef,form="unformatted",access="direct",recl=4)
c
c     read solvent site, grid
         read(ift,rec=1) nvuqx
         read(ift,rec=2) ngrid3dx
         read(ift,rec=3) ngrid3dy
         read(ift,rec=4) ngrid3dz
         close(ift)
         ng3dx=ngrid3dx*ngrid3dy*ngrid3dz
c
c     check consistency
         if (nvuq.ne.nvuqx .or. ng3d.ne.ng3dx) then
            write(*,*) "Error."
            write(*,*) "Reading 3d function inconsistency."
            ierr=568
            call abrt(ierr)
         endif
c
         open (ift,file=namef,form="unformatted",access="direct",recl=8)
         read(ift,rec=3) rnx
         read(ift,rec=4) rny
         read(ift,rec=5) rnz
         read(ift,rec=6) shiftx
         read(ift,rec=7) shifty
         read(ift,rec=8) shiftz
c
c     read function
         do j=1,nvuq
            do k=1,ng3d
               irec=(iv-1)*ng3d+ig+8
               read(ift,rec=irec)  func3d(k,j)
            enddo
         enddo
         
         close(ift)

      else
         write(*,*) "Error. Invalid outtype in $RISM."
         ierr=987
         call abrt(ierr)
      endif
c----------------------------------------------------------------
      return
      end
c**************************************************************
c----------------------------------------------------------------
c     Write 3D function to file
c----------------------------------------------------------------
      subroutine write3dfunc(namef,func3d,rdelta3d,nvuq,ngrid3d
     &                      ,outtype,char80)
c
      implicit real*8(a-h,o-z)
      character*256 namef,outtype,char256
      character*2 char2
      character*80 char80
c
      dimension func3d(ngrid3d**3,nvuq)
c
      ng3d=ngrid3d**3
      rn=ngrid3d*rdelta3d
      xyzshift=0.d0
c----------------------------------------------------------------
      if (outtype(1:5).eq."ASCII") then
         ift=45
         open (ift,file=namef)
         write(ift,9991) char80
         write(ift,9990) nvuq,ngrid3d,ngrid3d,ngrid3d
     &        rn,rn,rn,xyzshift,xyzshift,xyzshift
      
         do iv=1,nvuq
            do ig=1,ng3d
               write (ift,9992) func3d(ig,iv)
            enddo
         enddo

         close(ift)

      elseif (outtype(1:3).eq."BIN") then

         ift=45
         open (ift,file=namef,form="unformatted",access="direct",recl=4)
         write(ift,rec=1) nvuq
         write(ift,rec=2) ngrid3d
         write(ift,rec=3) ngrid3d
         write(ift,rec=4) ngrid3d
         close(ift)

         open (ift,file=namef,form="unformatted",access="direct",recl=8)
         write(ift,rec=3) rn
         write(ift,rec=4) rn
         write(ift,rec=5) rn
         write(ift,rec=6) xyzshift
         write(ift,rec=7) xyzshift
         write(ift,rec=8) xyzshift

         do iv=1,nvuq
            do ig=1,ng3d
               irec=(iv-1)*ng3d+ig+8
               write (ift,rec=irec) func3d(ig,iv)
            enddo
         enddo
         close(ift)

      else
         write(*,*) "Error. Invalid outtype in $RISM.",outtype,namef
         ierr=987
         call abrt(ierr)
      endif

c----------------------------------------------------------------
      return
 9990 format(4i8,4f16.8)
 9991 format("##  ",a80)
 9992 format(e16.8e3)
      end
c**************************************************************
c----------------------------------------------------------------
c     Write 3D function to file (complex)
c----------------------------------------------------------------
      subroutine write3dfuncz(namef,zfunc3d,rdelta3d,nvuq,ngrid3d
     &                      ,outtype,char80)
c
      implicit real*8(a-h,o-z)
      character*256 namef,outtype
      character*2 char2
      character*80 char80
      complex*16 zfunc3d
c
      dimension zfunc3d(ngrid3d**3,nvuq)
c
      ng3d=ngrid3d**3
      rn=ngrid3d*rdelta3d
      xyzshift=0.d0
c----------------------------------------------------------------
      if (outtype(1:5).eq."ASCII") then
         ift=45
         open (ift,file=namef)
         write(ift,9991) char80
         write(ift,9990) nvuq,ngrid3d,ngrid3d,ngrid3d
     &        rn,rn,rn,0.d0,0.d0,0.d0

         do iv=1,nvuq
            do ig=1,ng3d
               write (ift,9992) zfunc3d(ig,iv)
            enddo
         enddo
         close(ift)

      elseif (outtype(1:3).eq."BIN") then

         ift=45
         open (ift,file=namef,form="unformatted",access="direct",recl=4)
         write(ift,rec=1) nvuq
         write(ift,rec=2) ngrid3d
         write(ift,rec=3) ngrid3d
         write(ift,rec=4) ngrid3d
         close(ift)

         open (ift,file=namef,form="unformatted",access="direct",recl=8)
         write(ift,rec=3) rn
         write(ift,rec=4) rn
         write(ift,rec=5) rn
         write(ift,rec=6) xyzshift
         write(ift,rec=7) xyzshift
         write(ift,rec=8) xyzshift
         close(ift)

         open (ift,file=namef,form="unformatted",access="direct"
     &        ,recl=16)
         do iv=1,nvuq
            do ig=1,ng3d
               irec=(iv-1)*ng3d+ig+4
               write (ift,rec=irec) zfunc3d(ig,iv)
            enddo
         enddo

         close(ift)

      else
         write(*,*) "Error. Invalid outtype in $RISM."
         ierr=987
         call abrt(ierr)
      endif

c----------------------------------------------------------------
      return
 9990 format(4i8,4f16.8)
 9991 format("##  ",a80)
 9992 format(e16.8e3,2x,e16.8e3)
      end
c**************************************************************
c----------------------------------------------------------------
c     Write 3D function to file with xyz coordinate
c----------------------------------------------------------------
      subroutine write3dfuncxyz(namef,func3d,rdelta3d,nvuq,ngrid3d
     &                      ,char80)
c
      implicit real*8(a-h,o-z)
      character*256 namef
      character*2 char2
      character*80 char80
c
      dimension func3d(ngrid3d**3,nvuq)
c
c----------------------------------------------------------------
      ift=45
      nremark=2
      open (ift,file=namef)
      write(ift,9990) nvuq,ngrid3d,1,rdelta3d
      write(ift,9991) char80
      write(ift,9992) nremark

      write (ift,9993) "Number of points:",ngrid3d**3*nvuq
      write (ift,9996) "   x[Ang]   ","   y[Ang]   "
     &                ,"   z[Ang]   ","     q[e]       "
      
      k0=ngrid3d/2+1
      do iv=1,nvuq
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d

            rx=rdelta3d*dble(kx-k0)
            ry=rdelta3d*dble(ky-k0)
            rz=rdelta3d*dble(kz-k0)
            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2

            write (ift,9995) rx,ry,rz,func3d(k,iv)

         enddo
         enddo
         enddo
      enddo
      
      close(ift)
c----------------------------------------------------------------
      return
 9990 format("## 3D Function :",3i8,f16.8)
 9991 format("##  ",a80)
 9992 format("##  REMARKS ",i4)
 9993 format("##  ",a20,i15)
 9994 format(e16.8e3)
 9995 format(4x,3f12.4, 2x,e16.8e3,2x,e16.8e3)
 9996 format("##  ",3A12,2x,2A16)
      end
c**************************************************************
c----------------------------------------------------------------
c     Read electrostatic potential map
c----------------------------------------------------------------
      subroutine readespmap(namef,vres,rdelta3d,ngrid3d)
c
      implicit real*8(a-h,o-z)
c
      include "phys_const.i"
c
      character*256 namef
      character*2 char2
      character*4 char4
      character*15 char15
      character*20 char20
c
      dimension vres(ngrid3d**3)
c
      
c----------------------------------------------------------------
      ift=45
      open (ift,file=namef,status='old',err=999)
c
      write(*,*) "Reading electrostatic potential map from :",namef
c
c     Read ESP map
c
c     Format notes:
c     x y z V(x,y,z)
c     
c     Unit: xyz=[Angstrom], V=[Hartree/e]
c
      k0=ngrid3d/2+1
      do iv=1,ngrid3d**3

         read (ift,'(3f20.12,A20)') rx,ry,rz,char20
         read(char20,*,err=100) val
         goto 200
 100     continue
         write(*,'("ESP map contains non-numeric value ",A20
     &        " at ",3f20.12)') char20,rx,ry,rz
         write(*,*) "The ESP value set to 0.d0."
         val=0.d0
 200     continue

         kx=nint(rx/rdelta3d)+k0
         ky=nint(ry/rdelta3d)+k0
         kz=nint(rz/rdelta3d)+k0
         k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
         vres(k)=val*hart2jmol
 
      enddo
c      
      close(ift)
c----------------------------------------------------------------
      return
 999  continue
      write(*,*) "Error during read esp file:",trim(namef)
      ierr=4728
      call abrt(ierr)
      end
c**************************************************************
c----------------------------------------------------------------
c     Read RESP point charge
c----------------------------------------------------------------
      subroutine readresp(namef,qu,maxslu)
c
      implicit real*8(a-h,o-z)
      character*256 namef
      character*2 char2
c
      dimension qu(maxslu)
c----------------------------------------------------------------
      ift=45
      open (ift,file=namef,status='old')
c
c     Read point charges
c
      read(ift,*) char2
      read(ift,*) nu
      do i=1,nu
         read(ift,*) qu(i)
      enddo
c      
      close(ift)
c----------------------------------------------------------------
      return
      end
