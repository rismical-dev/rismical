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
c     xvk           ... k-space solvent susceptibility
c     
      implicit real*8 (a-h,o-z)
      character*6 char6
      character*256 scrjob
      character*80 char80,char802
      real*8 ,allocatable :: gbuff(:,:,:)

      include "phys_const.i"
      include "rismio.i"
      include "solvent.i"

      dimension cr(ngrid,n2,n2),tr(ngrid,n2,n2),fr(ngrid,n2,n2)
      dimension hvk(ngrid,n2,n2),ures(ngrid,n2,n2),urlj(ngrid,n2,n2)
      dimension xvk(ngrid,n2,n2)
      dimension ck(ngrid,n2,n2),fk(ngrid,n2,n2)
c----------------------------------------------------------------
      char802="REMARKS"
C
C     --- Individual format
C      
      allocate (gbuff(ngrid,n2,n2))
c
      koutg = index(iolist,'g') + index(iolist,'G')
      kouth = index(iolist,'h') + index(iolist,'H')
      koutu = index(iolist,'u') + index(iolist,'U')
      koutc = index(iolist,'c') + index(iolist,'C')
      koutt = index(iolist,'t') + index(iolist,'T')
      koutx = index(iolist,'x') + index(iolist,'X')
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
     &        ,n2,n2,ngrid,1,char80,char802)
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
     &        ,n2,n2,ngrid,1,char80,char802)

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
     &           ,n2,n2,ngrid,1,char80,char802)
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
     &           ,n2,n2,ngrid,1,char80,char802)
      endif
c
c     write hvk
c         
      if (kouth.ne.0) then

         scrjob=trim(basename)//".hvk"

         char80="hv(k) data with solvent parameters"
         call writehvvfunc(scrjob,hvk,rdelta,n2,ngrid,char80)
      endif
c
c
c     write xvvk
c         
      if (koutx.ne.0) then

         scrjob=trim(basename)//".xvk"

         char80="reduced xvv(k) data with solvent parameters"
c$$$         call writehvvfunc(scrjob,xvk,rdelta,n2,ngrid,char80)
         call writexvvfunc(scrjob,xvk,rdelta,n2,ngrid,char80)
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
