c**************************************************************
c--------------------------------------------------------------
c     3-Dimensional Closure Relations
c--------------------------------------------------------------
      subroutine closure3d(ng3d,n2uq,listcore
     &      ,ck,tr,vres,urlj, fr,fk)
c     
c     icl       ... closure type 0...hnc, 1...msa, 2...kh 3..hnc+rbc
c     ngrid3d   ... number of grid of 3d-rdf
c     rdelta3d  ... grid width of r-space [angstrom]
c     n2uq      ... number of symmetry uniq site of solvent
c     ck        ... (in)  toku ni nashi (nanndemo ok)
c               ... (out) k-space direct correlation function
c     tr        ... (in)  tau bond =h-c
c               ... (out) r-space direct correlation function
c     ures      ... electro static potential energy [erg]
c     urlj      ... lj potential energy [erg]
c     beta      ... inverse of kbt[1/erg]
c     dumfft    ... work space for 3d-fft
c     listcore  ... list of core region 0..core 1..not core
c     
      implicit real*8 (a-h,o-z)
      complex*16 ck,fk,fkk
      complex*16 ,allocatable :: dumfft(:)

      include "phys_const.i"
      include "rismrun.i"
      include "solvent.i"

      dimension ck(ng3d,n2uq),tr(ng3d,n2uq)
      dimension vres(ng3d)
      dimension urlj(ng3d,n2uq)
      dimension listcore(ng3d)
      dimension fr(ng3d),fk(ng3d)
c
      allocate (dumfft(ng3d))
c
      ngrid3d2=ngrid3d**2
      dk3d=2.d0*pi/(rdelta3d*dble(ngrid3d))
      ngshift=ngrid3d/2+1
c----------------------------------------------------------------
c
c     --- parameter for fft
c
      inv=1
      m=nint(dlog(dble(ngrid3d))/dlog(2.d0))
c
      call vclr_mp(ck,1,ng3d*2*nvuq)
c
c     --- solve closure relation ( t(r) --> c(r) )
c
      do j=1,nvuq
c
         call vclr_mp(dumfft,1,ng3d*2)
c
!$omp parallel do private(rz,ry,rx,k,bur,d,trs,dkr)
         do kz=1,ngrid3d
         rz=rdelta3d*dble(kz-ngshift)
         do ky=1,ngrid3d
         ry=rdelta3d*dble(ky-ngshift)
         do kx=1,ngrid3d
         rx=rdelta3d*dble(kx-ngshift)

            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d2
c     
c     --- inner core
c
            if (listcore(k).eq.0) then
               
               tr(k,j)=-tr(k,j)-1.d0
               
            else
c     
c     --- HNC
c     
               if (icl.eq.0) then
                  
                  bur=beta*(vres(k)*q2uq(j)*chgratio+urlj(k,j))
                  tr(k,j)=dexp(-bur+tr(k,j))-tr(k,j)-1.d0
c     
c     --- MSA
c     
               elseif (icl.eq.1) then
                  
                  bur=beta*(vres(k)*q2uq(j)*chgratio+urlj(k,j))
                  tr(k,j)=-bur
c     
c     --- KH
c     
               elseif (icl.eq.2) then
                  
                  bur=beta*(vres(k)*q2uq(j)*chgratio+urlj(k,j))
                  d=-bur+tr(k,j)
                  
                  if (d.gt.0.d0) then
                     tr(k,j)=d-tr(k,j)
                  else
                     tr(k,j)=dexp(d)-tr(k,j)-1.d0
                  endif
c     
c     --- HNC+RBC
c     
               elseif (icl.eq.3) then
                  
                  write(*,*) "RBC is not implemented yet."
                  ierr=44
                  call abrt(ierr)
c
c     --- ---
c
               endif
               
            endif
                                ! here, c(r) is in "tr"
c     
c     --- fourier transform ( c(r) --> c(k) )
c     
            trs=tr(k,j)+beta*fr(k)*q2uq(j)*chgratio
            dkr=dk3d*0.5d0*(rx+ry+rz)
            dumfft(k)=dcmplx(trs,0.d0)
     &           *cdexp(dcmplx(0.d0,dkr))

         enddo                  ! of kx
         enddo                  ! of ky
         enddo                  ! of kz
!$omp end parallel do

         call ffte3d(dumfft,ngrid3d,rdelta3d,inv)

!$omp parallel do
         do k=1,ng3d
            ck(k,j)=dumfft(k)-beta*fk(k)*q2uq(j)*chgratio ! c(k) is in "ck"
         enddo
!$omp end parallel do

      enddo                     ! of j to nvuq

c----------------------------------------------------------------
      deallocate (dumfft)
c
      return
      end
c**************************************************************
c--------------------------------------------------------------
c     3-Dimensional Closure Relations for postprocess
c--------------------------------------------------------------
      subroutine closure3d2(ng3d,n2uq,listcore
     &                    ,cr,tr,vres,urlj)
c     
c     icl       ... closure type 0...hnc, 1...msa, 2...kh 3..hnc+rbc
c     ngrid3d   ... number of grid of 3d-rdf
c     rdelta3d  ... grid width of r-space [angstrom]
c     nvuq      ... number of symmetry uniq site of solvent
c     cr        ... (out) r-space direct correlation function
c     tr        ... (in)  tau bond =h-c
c     ures      ... electro static potential energy [erg]
c     urlj      ... lj potential energy [erg]
c     listcore  ... list of core region 0..core 1..not core
c     
      implicit real*8 (a-h,o-z)
      complex*16 cr
      
      include "rismrun.i"
      include "solvent.i"
      include "phys_const.i"

      dimension cr(ng3d,n2uq),tr(ng3d,n2uq)
      dimension vres(ng3d)
      dimension urlj(ng3d,n2uq)
      dimension listcore(ng3d)
c
c----------------------------------------------------------------
c
c     --- solve closure relation ( t(r) --> c(r) )
c
      do j=1,nvuq
c
!$omp parallel do private(bur,d,dkr)
         do k=1,ng3d
c
c     --- inner core
c
            if (listcore(k).eq.0) then
               cr(k,j)=-tr(k,j)-1.d0
            else
c     
c     --- hnc
c     
               if (icl.eq.0) then
                  bur=beta*(vres(k)*q2uq(j)*chgratio+urlj(k,j))
                  cr(k,j)=dexp(-bur+tr(k,j))-tr(k,j)-1.d0
c     
c     --- msa
c     
               elseif (icl.eq.1) then
                  bur=beta*(vres(k)*q2uq(j)*chgratio+urlj(k,j))
                  cr(k,j)=-bur
c     
c     --- kh
c     
               elseif (icl.eq.2) then
                  
                  bur=beta*(vres(k)*q2uq(j)*chgratio+urlj(k,j))
                  d=-bur+tr(k,j)
                  
                  if (d.gt.0.d0) then
                     cr(k,j)=d-tr(k,j)
                  else
                     cr(k,j)=dexp(d)-tr(k,j)-1.d0
                  endif
c     
c     --- hnc+rbc
c     
               elseif (icl.eq.3) then
                  ierr=44
                  call abrt(ierr)
c
c     --- ---
c
               endif
               
            endif

         enddo   ! of kz
!$omp end parallel do
c
      enddo                     ! of j to nvuq

c----------------------------------------------------------------
      return
      end
