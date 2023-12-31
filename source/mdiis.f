c--------------------------------------------------------------
c     MDIIS
c--------------------------------------------------------------
      subroutine mdiis(ng,tr,rmsnew,cconv,nds)
c
c     ng        ... number of element in tr  [IN]
c     tr        ... tau bond =hr-cr          [IN/OUT]
c     nds       ... flag for initialize of MDIIS procedure [IN]
c     cconv     ... convergence criterion    [IN]
c     rmsnew    ... root mean squre error    [OUT]
c
c
      implicit real*8(a-h,o-z)
      character*1 char1

      logical init
      
      include "rismrun.i"
c
      real*8 ,allocatable,save::trmdiis(:,:)
      real*8 ,allocatable,save::rsmdiis(:,:)
      real*8 ,allocatable,save::rsdmat(:,:)
      integer ,allocatable,save::listrs(:)

      real*8 ,allocatable::dmat(:,:)
      real*8 ,allocatable::x(:)
c
      dimension tr(ng)

      data init/.true./
      save init,itr,idc,irmsmin
      save rmsmin

      facres=10.d0
c---------------------------------------------------------------
c
c     --- Setup MDIIS
c
      if (init) then
         
         allocate (trmdiis(ng,nsub))
         allocate (rsmdiis(ng,nsub))
         allocate (rsdmat(nsub+1,nsub+1))
         allocate (listrs(nsub))
         
         init=.false.

      endif
c     
c     --- Initialize MDIIS
c     
      if (nds.eq.0) then
         rsdmat(1,1)=0.d0
         do i=2,nsub+1
            do j=2,nsub+1
               rsdmat(i,j)=0.d0
            enddo
            rsdmat(i,1)=-1.d0
            rsdmat(1,i)=-1.d0
         enddo

         idc=0
         rmsmin=0.d0
         irmsmin=1
         itr=0

         listrs(1)=1
         do i=2,nsub
            listrs(i)=0
         enddo

         call vclr_mp(rsmdiis,1,ng*nsub)
         call vclr_mp(trmdiis,1,ng*nsub)
         do i=1,ng
            trmdiis(i,1)=tr(i)
         enddo

         return

      endif

      itr=itr+1
c
      allocate (dmat(nsub+1,nsub+1))
      allocate (x(nsub+1))
c      
      x(1)=-1.d0
      do i=2,nsub+1
         x(i)=0.d0
      enddo
c     
c     --- Calculate New RMS of Residual
c     
      rmsnew=0.d0
      irms=listrs(1)
!$OMP PARALLEL DO REDUCTION(+: RMSNEW)
      do jk=1,ng
         rmsnew=rmsnew+(tr(jk)-trmdiis(jk,irms))**2
      enddo
!$OMP END PARALLEL DO
      rmsnew=dsqrt(rmsnew/dble(ng))
      
      if ((rmsnew.le.cconv).and.(itr.ge.3)) return
c     
c     --- Dumping Parameter Setting
c     
      eta=dumpnume/rmsnew
      if (eta.gt.dumpmax) eta=dumpmax
      if (eta.lt.dumpmin) eta=dumpmin
c     
c     --- Check the lowest rms
c     
      if (idc.eq.0) rmsmin=rmsnew
      
      char1="F"
      if (rmsnew.le.rmsmin) then
         irmsmin=1
         rmsmin=rmsnew
         char1="T"
      endif
c     
c     --- If Convergence is Poor, Restart MDIIS
c     
      if ((rmsnew.gt.rmsmin*facres).and.(idc.gt.0)) then

         idc=1
         sum=0.0d0
         irms=listrs(irmsmin)

!$OMP PARALLEL DO REDUCTION(+: SUM)
         do ijk=1,ng
            sum=sum+dble(rsmdiis(ijk,irms)*rsmdiis(ijk,irms))
         enddo
!$OMP END PARALLEL DO

         rsdmat(3,3)=sum
         irmsmin=2
         rmsmin=sum
         listrs(2)=irms
      endif
c     
c     --- Increment Counter of MDIIS
c     
      irms=listrs(2)+1
      if (irms.gt.nsub) irms=1
      listrs(1)=irms
      idc=min(idc+1,nsub)
c     
c     --- Set New Residual
c     
      irms=listrs(1)
!$OMP PARALLEL DO
      do ijk=1,ng
         rsmdiis(ijk,irms)=dble(tr(ijk)-trmdiis(ijk,irms))
         trmdiis(ijk,irms)=tr(ijk)
      enddo
!$OMP END PARALLEL DO
c     
c     --- Make Residual Overlap Matrix
c     
      irms1=listrs(1)
      do i=1,idc
         sum=0.0d0
         irms=listrs(i)
!$OMP PARALLEL
!$OMP DO REDUCTION(+: SUM)
         do ijk=1,ng
            sum=sum+dble(rsmdiis(ijk,irms1)*rsmdiis(ijk,irms))
         enddo
!$OMP END DO
!$OMP END PARALLEL
         rsdmat(i+1,2)=sum
         rsdmat(2,i+1)=sum
      enddo
c     
c     --- Solve MDIIS Matrix Equation
c     
      do i=1,nsub+1
         do j=1,nsub+1
            dmat(i,j)=rsdmat(i,j)
         enddo
      enddo
      
      call matsolv(dmat,x,idc+1,nsub+1,1)
c     
c     --- Make New Estimate 
c     
!$OMP PARALLEL DO PRIVATE(SUM,IRMS)
      do ijk=1,ng

         sum=0.d0
         do m=1,idc
            irms=listrs(m)
            sum=sum+x(m+1)
     &           *(trmdiis(ijk,irms)+eta*dble(rsmdiis(ijk,irms)))
         enddo
         tr(ijk)=sum
      enddo
!$OMP END PARALLEL DO
c     
c     --- Shift Data rsdmat and Increment Minimum Flag
c     
      do i=nsub,2,-1
         do j=2,nsub
            rsdmat(i+1,j+1)=rsdmat(i,j)
         enddo
         listrs(i)=listrs(i-1)
      enddo
      irmsmin=irmsmin+1
c     
c     --- Save Newest Estimation
c     
      irmsnext=listrs(1)+1
      if (irmsnext.gt.nsub) irmsnext=1
      listrs(1)=irmsnext
!$OMP PARALLEL DO
      do ijk=1,ng
         trmdiis(ijk,irmsnext)=tr(ijk)
      enddo
!$OMP END PARALLEL DO
c     
c     --- Reset minimum point (if minimum point flow out)
c     
      if (irmsmin.gt.nsub) then 
         rmin=rsdmat(2,2)
         irmsmin=1
         do i=3,nsub+1
            if (rsdmat(i,i).lt.rmin) then
               rmin=rsdmat(i,i)
               irmsmin=i
            endif
         enddo
      endif
c     
c     --- Print Information
c     
      write(*,9999) itr,rmsnew,idc,char1,eta
c---------------------------------------------------------------
      deallocate (dmat)
      deallocate (x)
      return
 9999 format (4x,i6,f20.12,2x,i4,2x,a1,2x,f10.5)
      end
