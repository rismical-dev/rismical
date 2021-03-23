c----------------------------------------------------------------
c     Output 3D-RISM result
c----------------------------------------------------------------
      subroutine output3d(ng3d,n2,n2uq
     &                   ,cr,tr,urlj,vres,fr,fk)
c
c     ngr3d=ngrid3d ... number of grid of RDF
c     n2=nv         ... number of site of solvent
c     ck            ... k-space direct correlation function 
c     cr            ... direct correlation function 
c     tr            ... tau bond =hr-cr
c     vres          ... electro static field [erg/e]
c     urlj          ... LJ potential energy [erg]
c     
      implicit real*8 (a-h,o-z)
      complex*16 cr,fk
      character*6 char6
      character*80 char80
      character*256 scrjob

      real*8 ,allocatable :: gbuff(:,:)
      complex*16 ,allocatable :: dumfft(:,:,:),hk(:,:)
c
      include "phys_const.i"
      include "rismio.i"
      include "rismrun.i"
      include "solvent.i"
      include "solute.i"
c
      dimension cr(ng3d,n2uq),tr(ng3d,n2uq)
      dimension urlj(ng3d,n2uq)
      dimension fr(ng3d)
      dimension fk(ng3d)
      dimension vres(ng3d)
c
      DK3D=2.D0*PI/(RDELTA3D*DBLE(NGRID3D))
      DNSHIFT=DBLE(NGRID3D+1)/2.D0 
c----------------------------------------------------------------
c
C     --- Separate style
C      
      allocate (gbuff(ng3d,nvuq))
c
      koutg = index(iolist,'g') + index(iolist,'G')
      kouth = index(iolist,'h') + index(iolist,'H')
      koutu = index(iolist,'u') + index(iolist,'U')
      koutv = index(iolist,'v') + index(iolist,'V')
      koutc = index(iolist,'c') + index(iolist,'C')
      koutf = index(iolist,'f') + index(iolist,'F')
      koutt = index(iolist,'t') + index(iolist,'T')
      koutq = index(iolist,'q') + index(iolist,'Q')
c
c     write guv
c         
      if (koutg.ne.0) then

         scrjob=trim(basename)//".guv"
            
         do iv=1,nvuq
            do ig=1,ng3d
               gbuff(ig,iv)=tr(ig,iv)+dble(cr(ig,iv))+1.d0
            enddo
         enddo
         char80="guv data"
         call write3dfunc(scrjob,gbuff,rdelta3d
     &        ,nvuq,ngrid3d,1,char80)
      endif
c
c     write ur
c         
      if (koutu.ne.0) then

         scrjob=trim(basename)//".uuv"

         do iv=1,nvuq
            do ig=1,ng3d
               gbuff(ig,iv)=(vres(ig)*q2uq(iv)+urlj(ig,iv))*beta
            enddo
         enddo
         char80="uuv data"
         call write3dfunc(scrjob,gbuff,rdelta3d
     &           ,nvuq,ngrid3d,1,char80)
      endif
c         
c     write vres
c         
      if (koutv.ne.0) then

         scrjob=trim(basename)//".vuv"

         do ig=1,ng3d
            gbuff(ig,1)=vres(ig)
         enddo
         char80="vres data (electrostatic field due to solute)"
         call write3dfunc(scrjob,gbuff,rdelta3d,1,ngrid3d,1,char80)
      endif
c
c     write cr
c         
      if (koutc.ne.0) then

         scrjob=trim(basename)//".cuv"
         
         char80="cuv data"
         call write3dfunc(scrjob,cr,rdelta3d,nvuq,ngrid3d,2,char80)
            
      endif
c
c     write tr
c         
      if (koutt.ne.0) then

         scrjob=trim(basename)//".tuv"

         char80="tuv data"
         call write3dfunc(scrjob,tr,rdelta3d,nvuq,ngrid3d,1,char80)
         
      endif
c     
c     write hk
c         
      if (kouth.ne.0) then
            
c     calc H(k)
c
         allocate (hk(ng3d,nvuq))

         inv=1                  ! r -> k
         do iv=1,nvuq
            do ig=1,ng3d
               hk(ig,iv)=dcmplx(tr(ig,iv),0.d0)+cr(ig,iv)
            enddo
         enddo

         do iv=1,nvuq
            call ft3dfunc(hk(1,iv),ngrid3d,rdelta3d,inv)
         enddo

         scrjob=trim(basename)//".huvk"
         char80="huv(k) data"
         call write3dfunc(scrjob,hk,rdelta3d,nvuq,ngrid3d,2,char80)
            
         deallocate (hk)

      endif
c
c     write fr and fk
c         
      if (koutf.ne.0) then
         scrjob=trim(basename)//".fuv"
         
         char80="f(r) data"
         call write3dfunc(scrjob,fr,rdelta3d,1,ngrid3d,1,char80)
         
         scrjob=trim(basename)//".fuvk"

         char80="f(k) data"
         call write3dfunc(scrjob,fk,rdelta3d,1,ngrid3d,2,char80)
            
      endif
c
c     write charge distribution of solvent
c         
      if (koutq.ne.0) then

         scrjob=trim(basename)//".qv"
            
         do ig=1,ng3d

            sum=0.d0
            do iv=1,nvuq
               sum=sum + densuq(iv)*q2uq(iv)
     &              *(tr(ig,iv)+dble(cr(ig,iv))+1.d0)
     &              *rdelta3d**3
            enddo

            gbuff(ig,1)=sum
         
         enddo

         char80="charge distribution of solvent data"
         call write3dfuncxyz(scrjob,gbuff,rdelta3d
     &                      ,1,ngrid3d,1,char80)
      endif
c
c
      deallocate (gbuff)
C      
c----------------------------------------------------------------
      return
      end
c**************************************************************
c----------------------------------------------------------------
c     Write 3D function to file
c----------------------------------------------------------------
      subroutine write3dfunc(namef,func3d,rdelta3d,nvuq,ngrid3d,ncmp
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
      nremark=0
      open (ift,file=namef)
      write(ift,9990) nvuq,ngrid3d,ncmp,rdelta3d
      write(ift,9991) char80
      write(ift,9992) nremark
      do i=1,nremark
         write (ift,9993) "remarks "
      enddo
      
      do iv=1,nvuq
         do ig=1,ngrid3d**3
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
