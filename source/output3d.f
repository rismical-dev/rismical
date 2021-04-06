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
     &        ,nvuq,ng3d,char80)
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
     &           ,nvuq,ng3d,char80)
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
         call write3dfunc(scrjob,gbuff,rdelta3d,1,ng3d,char80)
      endif
c
c     write cr
c         
      if (koutc.ne.0) then

         scrjob=trim(basename)//".cuv"
         
         char80="cuv data"
         call write3dfuncz(scrjob,cr,rdelta3d,nvuq,ng3d,char80)
            
      endif
c
c     write tr
c         
      if (koutt.ne.0) then

         scrjob=trim(basename)//".tuv"

         char80="tuv data"
         call write3dfunc(scrjob,tr,rdelta3d,nvuq,ng3d,char80)
         
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
         call write3dfuncz(scrjob,hk,rdelta3d,nvuq,ng3d,char80)
            
         deallocate (hk)

      endif
c
c     write fr and fk
c         
      if (koutf.ne.0) then
         scrjob=trim(basename)//".fuv"
         
         char80="f(r) data"
         call write3dfunc(scrjob,fr,rdelta3d,1,ng3d,char80)
         
         scrjob=trim(basename)//".fuvk"

         char80="f(k) data"
         call write3dfuncz(scrjob,fk,rdelta3d,1,ng3d,char80)
            
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
     &                      ,1,ngrid3d,char80)
      endif
c
c
      deallocate (gbuff)
C      
c----------------------------------------------------------------
      return
      end
