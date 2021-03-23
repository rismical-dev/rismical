c--------------------------------------------------------------
c     1-Dimensional MDIIS
c--------------------------------------------------------------
      subroutine mdiis1d(itr,nsub,ngrid,n1,n2
     &                  ,nds,idc,irmsmin
     &                  ,tr,trnew,trmdiis,rsmdiis,rmsmin,rsdmat
     &                  ,rmsnew,dumpmax,dumpmin,dumpnume)
c
c     ir        ... file number of output (STDOUT)
c     iw        ... file number of input (STDIN)
c     ngrid     ... number of grid of RDF
c     n1        ... number of site of 1 (solute or solvent)
c     n2        ... number of site of 2 (solvent)
c     nsub      ... number of MDIIS sub-space
c     tr        ... tau bond =hr-cr
c     trnew     ... latest estimation of tau bond
c     trmdiis   ... stock of tr for MDIIS sub-space
c     rsmdiis   ... stock of residual for MDIIS sub-space
c     rsdmat    ... residual overlap matrix
c     nds       ... flag for initialize of MDIIS procedure
c     idc       ... counter for MDIIS sub-space
c     irmsmin   ... index of rmsmin
c     dumpmax   ... maximum value of dumping parameter eta
c     dumpmin   ... manimum value of dumping parameter eta
c     dumpnume  ... numerator of dumping function
c
      implicit real*8(a-h,o-z)
      character*1 char1

      dimension tr(ngrid,n1,n2),trnew(ngrid,n1,n2)
      dimension trmdiis(ngrid*n1*n2,nsub),rsmdiis(ngrid*n1*n2,nsub)
      dimension rsdmat(nsub+1,nsub+1),dmat(nsub+1,nsub+1),x(nsub+1)
      dimension list(nsub+1),w(nsub+1)
c
c     list      ... list vector for LU decomposition
c     w         ... work space for LU decomposition
c
c     --- parameter of check residual ---
      facres=10.d0
c---------------------------------------------------------------
c     
c     --- Initialize MDIIS
c     
      if (nds.eq.0) then
         nds=1
         idiiscount=0
         rsdmat(1,1)=0.d0
         do i=2,nsub+1
            do j=2,nsub+1
               rsdmat(i,j)=0.d0
            enddo
            rsdmat(i,1)=-1.d0
            rsdmat(1,i)=-1.d0
         enddo
      endif

      x(1)=-1.d0
      do i=2,nsub+1
         x(i)=0.d0
      enddo
c     
c     --- Calculate New RMS of Residual
c     
      rmsnew=0.d0
      do i=1,n1
         do j=1,n2
            do k=1,ngrid
              rmsnew=rmsnew+(tr(k,i,j)-trnew(k,i,j))**2
           enddo
        enddo
      enddo
      rmsnew=dsqrt(rmsnew/dble(n1*n2*ngrid))
c
c     --- Dumping Parameter Setting
c
cc      eta=0.1d0/rmsnew
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
         sum=0.d0
         do ijk=1,n1*n2*ngrid
            sum=sum+rsmdiis(ijk,irmsmin)*rsmdiis(ijk,irmsmin)
            trmdiis(ijk,2)=trmdiis(ijk,irmsmin)
            rsmdiis(ijk,2)=rsmdiis(ijk,irmsmin)
         enddo
         rsdmat(3,3)=sum
         irmsmin=1
         rmsmin=sum
      endif
c
c     --- Increment Counter of MDIIS
c
      idc=idc+1
      if (idc.gt.nsub) idc=nsub
c
c     --- Set New Residual and New tr 
c
      do i=1,n1
         do j=1,n2
            do k=1,ngrid
               ijk=(i-1)*n2*ngrid+(j-1)*ngrid+k
               trmdiis(ijk,1)=trnew(k,i,j)
               rsmdiis(ijk,1)=trnew(k,i,j)-tr(k,i,j)
            enddo
         enddo
      enddo
c
c     --- Make Residual Overlap Matrix
c
      do i=1,idc
         sum=0.d0
         do ijk=1,n1*n2*ngrid
            sum=sum+rsmdiis(ijk,1)*rsmdiis(ijk,i)
         enddo
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
      
      eps=1.d-16
      ind=0
      det=0.d0

      call matsolv(dmat,x,idc+1,nsub+1,1)
c     
c     --- Make New Estimate 
c     
      do i=1,n1
         do j=1,n2
            do k=1,ngrid
               ijk=(i-1)*n2*ngrid+(j-1)*ngrid+k
               sum=0.d0
               do m=1,idc
                  sum=sum+x(m+1)*(trmdiis(ijk,m)+eta*rsmdiis(ijk,m))
               enddo
               tr(k,i,j)=sum
            enddo
         enddo
      enddo
c
c     --- Shift Data trmdiis, rsmdiis and rsdmat
c
      do i=nsub,2,-1
         do j=2,nsub
            rsdmat(i+1,j+1)=rsdmat(i,j)
         enddo
      enddo
      
      do m=nsub,2,-1
         do ijk=1,n1*n2*ngrid
            trmdiis(ijk,m)=trmdiis(ijk,m-1)
            rsmdiis(ijk,m)=rsmdiis(ijk,m-1)
         enddo
      enddo
      irmsmin=irmsmin+1
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
c               irmsmin=irmsmin
            endif
         enddo
      endif
c
c     --- Print Information
c
      write(*,9999) itr,rmsnew,idc,char1,eta
c---------------------------------------------------------------
      return
 9999 format (i6,f20.12,2x,i4,2x,a1,2x,f10.5)
      end
c-------------------------------------------------------------
c     Calculate residual
c-------------------------------------------------------------
      subroutine residucalc(ngrid,n1,n2,tr,trnew,residu)
c
c     ngrid     ... number of grid of RDF
c     n1        ... number of site of 1 (solute or solvent)
c     n2        ... number of site of 2 (solvent)
c     tr        ... tau bond =hr-cr
c     trnew     ... latest estimation of tau bond
c     residu    ... residual of tau-bond
c
      implicit real*8 (a-h,o-z)

      dimension tr(ngrid,n1,n2),trnew(ngrid,n1,n2)
c-----------------------------------------------------------------
      residu=0.d0

      do i=1,n1
         do j=1,n2
            do k=1,ngrid
               residu=residu+(trnew(k,i,j)-tr(k,i,j))**2
            enddo
         enddo
      enddo

      residu=dsqrt(residu/dble(ngrid*n1*n2))
c-----------------------------------------------------------------
      return
      end
