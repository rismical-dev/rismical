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
      subroutine writexvvfunc(namef,hvk,rdelta,n2,ngrid,char80)
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
      write(ift,'(A3,i4,1x,i4,1x,f16.8,e16.8)') "## ",numspc,nv,temp,xt
      write(ift,9801)
      do i=1,n2
         densm=dens(nspc(i))/(avognum*1.D-27)    ! to M
         write(ift,9800) nsitev(i),nspc(i),sigljv(i),epsljv(i),qv(i)
     &           ,xyzv(1,i),xyzv(2,i),xyzv(3,i),densm
      enddo
      write(ift,'(A2)') "##"
 9801 format ("## ATOM"," SPC"
     &     ," sig[Angs]  "," eps[J/mol] "," charge[e]  "
     &     ,"  ---X---   ","  ---Y---   ","  ---Z---   "," density[M] ")
 9800 format ("## ",A4,1x,i3,7f12.5)
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
      end
c**************************************************************
c----------------------------------------------------------------
c     Read 3D function from file
c----------------------------------------------------------------
      subroutine read3dfunc(namef,func3d,nvuq,ng3d,ncmp)
c
c     ncmp : 1... real function,  2...complex function
c      
      implicit real*8(a-h,o-z)
      character*256 namef
      character*2 char2
c
      dimension func3d(ncmp,ng3d,nvuq)
      dimension dum(ncmp)
c
c----------------------------------------------------------------
      ift=45
      open (ift,file=namef)
c
c     read header
      read(ift,*) char2 nvuqx,ng3dx,ncmpx,rdelta3dx
      read(ift,*) char2
      read(ift,*) char2,ndum
c
c     check consistency
      if (nvuq.ne.nvuqx .or. ng3d.ne.ng3dx
     &     .or. ncmp.ne.ncmpx ) then
         write(*,*) "Error."
         write(*,*) "Reading 3d function inconsistency."
         ierr=568
         call abrt(ierr)
      endif
c
c     skip remark
      do i=1,ndum
         read(ift,*) char2
      enddo
c
c     read function
      do j=1,nvuq
         do k=1,ng3d
            read(ift,*) (dum(i),i=1,ncmp)
            do i=1,ncmp
               func3d(i,k,j)=dum(i)
            enddo
         enddo
      enddo
      
      close(ift)
c----------------------------------------------------------------
      return
      end
c**************************************************************
c----------------------------------------------------------------
c     Write 3D function to file
c----------------------------------------------------------------
      subroutine write3dfunc(namef,func3d,rdelta3d,nvuq,ng3d,ncmp
     &                      ,char80)
c
c     ncmp : 1... real function,  2...complex function
c      
      implicit real*8(a-h,o-z)
      character*256 namef
      character*2 char2
      character*80 char80
c
      dimension func3d(ncmp,ng3d,nvuq)
c
c----------------------------------------------------------------
      ift=45
      nremark=0
      open (ift,file=namef)
      write(ift,9990) nvuq,ng3d,ncmp,rdelta3d
      write(ift,9991) char80
      write(ift,9992) nremark
      do i=1,nremark
         write (ift,9993) "remarks "
      enddo
      
      do iv=1,nvuq
         do ig=1,ng3d
            write (ift,9995) (func3d(icmp,ig,iv),icmp=1,ncmp)
         enddo
      enddo
      
      close(ift)
c----------------------------------------------------------------
      return
 9990 format("## 3D Function :",3i8,f16.8)
 9991 format("##  ",a80)
 9992 format("##  ",i4)
 9993 format("##  ",a80)
 9994 format(e16.8e3)
 9995 format(e16.8e3,2x,e16.8e3)
      end
c**************************************************************
c----------------------------------------------------------------
c     Write 3D function to file with xyz coordinate
c----------------------------------------------------------------
      subroutine write3dfuncxyz(namef,func3d,rdelta3d,nvuq,ngrid3d,ncmp
     &                      ,char80)
c
c     ncmp : 1... real function,  2...complex function
c      
      implicit real*8(a-h,o-z)
      character*256 namef
      character*2 char2
      character*80 char80
c
      dimension func3d(ncmp,ngrid3d**3,nvuq)
c
c----------------------------------------------------------------
      ift=45
      nremark=2
      open (ift,file=namef)
      write(ift,9990) nvuq,ngrid3d,ncmp,rdelta3d
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

            write (ift,9995) rx,ry,rz,(func3d(icmp,k,iv),icmp=1,ncmp)

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
