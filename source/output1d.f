c----------------------------------------------------------------
c     Output 1D-RISM result
c----------------------------------------------------------------
      subroutine output1d(ngrid,rdelta,n1,n2
     &                   ,cr,tr,ures,urlj,fr,ck,fk)
c
c     iuv           ... OZ type 0...vv, 1...uv
c     ngrid         ... number of grid of RDF
c     rdelta        ... grid width of r-space
c     n1            ... number of site of 1 (solute or solvent)
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
      character*6 char6
      character*256 scrjob
      character*80 char80
      real*8 ,allocatable :: gbuff(:,:,:)

      include "phys_const.i"
      include "rismio.i"
      include "solvent.i"
      include "solute.i"

      dimension cr(ngrid,n1,n2),tr(ngrid,n1,n2),fr(ngrid,n1,n2)
      dimension ures(ngrid,n1,n2),urlj(ngrid,n1,n2)
      dimension ck(ngrid,n1,n2),fk(ngrid,n1,n2)

c----------------------------------------------------------------
C
C     --- Individual format
C      
      allocate (gbuff(ngrid,n1,n2))
c
      koutg = index(iolist,'g') + index(iolist,'G')
      koutu = index(iolist,'u') + index(iolist,'U')
      koutv = index(iolist,'v') + index(iolist,'V')
      koutc = index(iolist,'c') + index(iolist,'C')
      koutt = index(iolist,'t') + index(iolist,'T')
      koutk = index(iolist,'k') + index(iolist,'K')
c
c     write gr
c         
      if (koutg.ne.0) then

         scrjob=trim(basename)//".guv"
         
         do j=1,n2
            do i=1,n1
               do ig=1,ngrid
                  gbuff(ig,i,j)=tr(ig,i,j)+cr(ig,i,j)+1.d0
               enddo
            enddo
         enddo
         char80="g(r) data"
         call write1dfunc(scrjob,gbuff,rdelta
     &           ,n1,n2,ngrid,char80)
      endif
c
c     write ur
c         
      if (koutu.ne.0) then

         scrjob=trim(basename)//".uuv"

         do j=1,n2
            do i=1,n1
               do ig=1,ngrid
                  gbuff(ig,i,j)=beta*(ures(ig,i,j)+urlj(ig,i,j))
               enddo
            enddo
         enddo

         char80="beta x u(r) data"
         call write1dfunc(scrjob,gbuff,rdelta
     &           ,n1,n2,ngrid,char80)

      endif
c
c     write cr
c         
      if (koutc.ne.0) then
            
         scrjob=trim(basename)//".cuv"

         do j=1,n2
            do i=1,n1
               do ig=1,ngrid
                  gbuff(ig,i,j)=cr(ig,i,j)
               enddo
            enddo
         enddo
         char80="c(r) data"
         call write1dfunc(scrjob,gbuff,rdelta
     &           ,n1,n2,ngrid,char80)
      endif
c
c     write ck
c         
      if (koutk.ne.0) then
            
         scrjob=trim(basename)//".cuvk"

         do j=1,n2
            do i=1,n1
               do ig=1,ngrid
                  gbuff(ig,i,j)=ck(ig,i,j)
               enddo
            enddo
         enddo
         char80="c(k) data"
         call write1dfunc(scrjob,gbuff,rdelta
     &           ,n1,n2,ngrid,char80)
      endif
c
c
c     write tr
c         
      if (koutt.ne.0) then
            
         scrjob=trim(basename)//".tuv"

         do j=1,n2
            do i=1,n1
               do ig=1,ngrid
                  gbuff(ig,i,j)=tr(ig,i,j)
               enddo
            enddo
         enddo
         char80="t(r) data"
         call write1dfunc(scrjob,gbuff,rdelta
     &           ,n1,n2,ngrid,char80)
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
