c----------------------------------------------------------------
c     Output solvent-solvent RISM result
c----------------------------------------------------------------
      subroutine outputvv(ngrid,rdelta,n2
     &                   ,cr,tr,hvk,xvk,ures,urlj,fr,ck,fk)
c
c     ngrid         ... number of grid of RDF
c     rdelta        ... grid width of r-space
c     n2            ... number of site of 2 (solvent)
c     nsite1,nsite2 ... name of site
c     cr            ... direct correlation function 
c     ck            ... k-space direct correlation function 
c     tr            ... tau bond =hr-cr
c     ures          ... electro static potential [erg]
c     urlj          ... LJ potential [erg]
c     hvk           ... k-space total correlation function
c     
      implicit real*8 (a-h,o-z)
      character*4 nsite1,nsite2
      character*6 char6
      character*256 scrjob
      character*80 char80
      real*8 ,allocatable :: gbuff(:,:,:)

      include "phys_const.i"
      include "rismio.i"
      include "solvent.i"

      dimension nsite2(n2)
      dimension cr(ngrid,n2,n2),tr(ngrid,n2,n2),fr(ngrid,n2,n2)
      dimension hvk(ngrid,n2,n2),ures(ngrid,n2,n2),urlj(ngrid,n2,n2)
      dimension xvk(ngrid,n2,n2)
      dimension ck(ngrid,n2,n2),fk(ngrid,n2,n2)

c----------------------------------------------------------------
c$$$c
c$$$c     --- Output For File 
c$$$c
c$$$      char6=".rsmvv"
c$$$
c$$$      ift=45
c$$$      open (ift,file=trim(basename)//char6)
c$$$
c$$$      write(ift,'(A2,1x,i10,f10.5)') "##",ngrid,rdelta
c$$$
c$$$      do i=1,n2
c$$$         do j=1,n2
c$$$
c$$$            write(ift,9997) nsitev(i),nsitev(j)
c$$$
c$$$            do k=1,ngrid
c$$$
c$$$               gr=tr(k,i,j)+cr(k,i,j)+1.d0
c$$$
c$$$               write(ift,'(4e16.8)') cr(k,i,j),tr(k,i,j),gr,hvk(k,i,j)
c$$$            enddo
c$$$            write(ift,*)
c$$$
c$$$         enddo
c$$$      enddo
c$$$      close(ift)
C
C     --- Separate style
C      
      allocate (gbuff(ngrid,n2,n2))
c
      koutg = index(iolist,'g') + index(iolist,'G')
      kouth = index(iolist,'h') + index(iolist,'H')
      koutu = index(iolist,'u') + index(iolist,'U')
      koutv = index(iolist,'v') + index(iolist,'V')
      koutc = index(iolist,'c') + index(iolist,'C')
      koutt = index(iolist,'t') + index(iolist,'T')
c
c     write gr
c         
      if (koutg.ne.0) then

         scrjob=trim(basename)//".gvv"

         do j=1,n2
            do i=1,n2
               do ig=1,ngrid
                  gbuff(ig,i,j)=tr(ig,i,j)+cr(ig,i,j)+1.d0
               enddo
            enddo
         enddo
         char80="g(r) data"
         call write1dfunc(scrjob,gbuff,rdelta
     &        ,n2,n2,ngrid,char80)
      endif
c
c     write ur
c         
      if (koutu.ne.0) then

         scrjob=trim(basename)//".uvv"

         do j=1,n2
            do i=1,n2
               do ig=1,ngrid
                  gbuff(ig,i,j)=beta*(ures(ig,i,j)+urlj(ig,i,j))
               enddo
            enddo
         enddo

         char80="beta x u(r) data"
         call write1dfunc(scrjob,gbuff,rdelta
     &        ,n2,n2,ngrid,char80)

      endif
c
c     write cr
c         
      if (koutc.ne.0) then
         
         scrjob=trim(basename)//".cvv"
         
         do j=1,n2
            do i=1,n2
               do ig=1,ngrid
                  gbuff(ig,i,j)=cr(ig,i,j)
               enddo
            enddo
         enddo
         char80="c(r) data"
         call write1dfunc(scrjob,gbuff,rdelta
     &           ,n2,n2,ngrid,char80)
      endif
c
c     write tr
c         
      if (koutt.ne.0) then
         
         scrjob=trim(basename)//".tvv"
         
         do j=1,n2
            do i=1,n2
               do ig=1,ngrid
                  gbuff(ig,i,j)=tr(ig,i,j)
               enddo
            enddo
         enddo
         char80="t(r) data"
         call write1dfunc(scrjob,gbuff,rdelta
     &           ,n2,n2,ngrid,char80)
      endif
c
c     write hvk
c         
      if (kouth.ne.0) then

         scrjob=trim(basename)//".hvk"

         char80="hv(k) data with solvent parameters"
         call writexvvfunc(scrjob,hvk,rdelta,n2,ngrid,char80)
      endif
c
      deallocate (gbuff)
      
c----------------------------------------------------------------
      return
 9997 format ("##<",4x,a4,">-<",a4,">")
 9998 format (E16.8,1x,10E20.12,1x,e16.8)
 9999 format ("##<",4x,a4,">-<",a4,">",/,
     &        "## r[Ang]",10x,"cr(r)",15x,"tr(r)",15x,"gr(r)",15x,
     &        "hv(k)",15x,"ur(r)",15x,"f-bond(r)",11x,"CN(r)",15x,
     &        "xvv(k)",14x,"c(k)",15x,"f(k)",15x,"k[/Ang]")
      end
c----------------------------------------------------------------
c     Write 1D function to file
c----------------------------------------------------------------
      subroutine write1dfunc(namef,func1d,rdelta,n1,n2,ngrid
     &                      ,char80)
c
      implicit real*8(a-h,o-z)
      character*256 namef
      character*2 char2
      character*80 char80
c
      include "phys_const.i"

      dimension func1d(ngrid,n1,n2)
c
c----------------------------------------------------------------
      ift=45
      nremark=0
      open (ift,file=namef)
      write(ift,9990) n1,n2,ngrid,rdelta
      write(ift,9991) char80
      write(ift,9992) nremark
      do i=1,nremark
         write (ift,9993) "remarks "
      enddo

      do i2=1,n2
      do i1=1,n1
         do ig=1,ngrid
            write (ift,9995) func1d(ig,i1,i2)
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
