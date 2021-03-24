c**************************************************************
c---------------------------------------------------------
c     3D-RISM for U-V system 
c---------------------------------------------------------
      subroutine rismical3d
c
      implicit real*8 (a-h,o-z)
      real*8 ,allocatable :: tr(:,:),fr(:)
     &                      ,xvv(:,:,:)
     &                      ,urlj(:,:),vres(:)
      complex*16 ,allocatable :: cr(:,:),ck(:,:),fk(:)
      integer ,allocatable :: listcore(:),listxvv(:,:,:)
c
      include "phys_const.i"
      include "solute.i"
      include "solvent.i"
      include "rismrun.i"
c
c----------------------------------------------------------------------
c
c     --- get RISM run parameters
c
      call readinput
c      
c     --- get solution  parameters
c
      call read3duvdata
c
c     --- allocate array
c
      call setuparraysize

      ng3d=ngrid3d**3
      allocate (ck(ng3d,nvuq))
      allocate (cr(ng3d,nvuq))
      allocate (tr(ng3d,nvuq))
      allocate (vres(ng3d))
      allocate (urlj(ng3d,nvuq))
      allocate (xvv(nxvv,nv,nv))
      allocate (listxvv(ngrid3d/2+1,ngrid3d/2+1,ngrid3d/2+1))
      allocate (listcore(ng3d))
      allocate (fr(ng3d),fk(ng3d))
c
c     --- Initialize
c     
      call vclr_mp(cr,1,ng3d*2*nvuq)
      call vclr_mp(ck,1,ng3d*2*nvuq)
      call vclr_mp(tr,1,ng3d*nvuq)
c     
c     --- Setup V-V Total Correlation Function 
c     
      call setup1dvx(ngrid,nv,nxvv,ngrid3d,listxvv,xvv)
c     
c     --- Make 3D f-Bond
c     
      call fbond3duv(ng3d,fr,fk)
c     
c     --- Make 3D-Potential 
c     
      call potential3duv(ng3d,nvuq,vres,urlj,listcore)
c     
c     --- Make initial guess for tr(r)
c     
      frfac=beta*chgratio
c     
      if (iguess.eq.0) then
         do j=1,nvuq
            do k=1,ng3d
               tr(k,j)=0.d0
               jk=(j-1)*ng3d+k
               trguess=fr(k)*q2uq(j)
               tr(k,j)=trguess*frfac
            enddo
         enddo
         
      elseif (iguess.eq.1) then

         call readguess3d(ng3d,nvuq,tr)

      else
         write(*,9990)
         ierr=199
         call abrt(ierr)
      endif
c---------------------------------------------------------
c     Charge Up Procedure
c---------------------------------------------------------
      chgratio=chgratio-chgstep
 1000 continue
c
c     --- Charge up
c
      chgratio=chgratio+chgstep
      if (chgratio.gt.1.d0) chgratio=1.d0
c
c     --- Make initial guess for tr(r) in charge up cycle
c
      prefac=frfac
      frfac =beta*chgratio
      do j=1,nvuq
!$OMP PARALLEL DO PRIVATE(trguess)
         do k=1,ng3d
            trguess=fr(k)*q2uq(j)
            tr(k,j)=tr(k,j)-trguess*prefac
     &                     +trguess*frfac
         enddo
!$OMP END PARALLEL DO
      enddo
c
c     --- Set Convergence Criterion
c
      if (chgratio.eq.1.d0) chgconv=conv
      cconv=chgconv
c
c     --- Setup mdiis
c     
      call mdiis(ng3d*nvuq,tr,residu,cconv,0)  
c---------------------------------------------------------
c     RISM Iteration Cycle
c---------------------------------------------------------
      write(*,9995)
      write(*,9999) temp,chgratio
      write(*,9991)
c     
      do itr=1,itrmax
c     
c     --- closure ( cr(r)=f[tr(r)] )
c     
c                 in      out
c          "ck"  ----     c(k)
c          "tr"  t(r)     c(r)
c     
         call closure3d(ng3d,nvuq,listcore
     &        ,ck,tr,vres,urlj,fr,fk)
c
c     --- oz ( tr(r)=f[cr(r)] )
c     
c                 in      out
c          "ck"  c(k)     c(k)
c          "cr"  ----     c(r)
c          "tr"  c(r)     t(r)
c     
         call oz3d(ngrid3d,nxvv,nv,nvuq
     &            ,listxvv,ck,cr,tr,xvv)
c    
c     --- check convergence and make guess for next loop
c     
         call mdiis(ng3d*nvuq,tr,residu,cconv,1)
c     
         if ((residu.le.cconv).and.(itr.ge.3)) goto 8000
c     
      enddo
c---------------------------------------------------------
c     not converged
c---------------------------------------------------------
      write(*,9997)
      ierr=799
      call abrt(ierr)
c---------------------------------------------------------
c     converged
c---------------------------------------------------------
 8000 continue
c
      write(*,9989) itr,residu,0,"x"
      write(*,9993)
c
c     go to next charge up cycle
c
      if (chgratio.ne.1.d0) goto 1000
c
c---------------------------------------------------------
c     output 3d-rism result
c---------------------------------------------------------
c
c     calc cr for final tr
c
      call closure3d2(ng3d,nvuq,listcore
     &               ,cr,tr,vres,urlj)
c
c     calc property
c
      call prop3duv(ng3d,nv,nvuq
     &             ,vres,urlj,listcore,cr,tr)
c
c     output 3D-DFs
c
      call output3d(ng3d,nv,nvuq
     &             ,cr,tr,urlj,vres,fr,fk)
c---------------------------------------------------------

      deallocate (ck)
      deallocate (cr)
      deallocate (tr)
      deallocate (xvv)
      deallocate (vres)
      deallocate (urlj)
      deallocate (listcore)
      deallocate (listxvv)
      deallocate (fr,fk)
c---------------------------------------------------------
      return
c---------------------------------------------------------
 9989 format (4x,i6,f20.12,2x,i4,2x,a1,4x,"converged")
 9990 format (/,4x,"error, wrong parameter for iguess in $rismguess")
 9991 format (/,3x,"itr",6x,"residual",7x,"#-sub",1x,"min",6x,"dump")
 9993 format (/,4x,"rism cycle is converged")
 9995 format (/,4x,"========== u-v 3drism ==========")
 9997 format (/,4x,"rism cycle is not converged.",
     &        /,4x,"----- e n d -----")
 9999 format (/,4x,"temperature      :",f10.5,"[k]",
     &        /,4x,"charge up factor :",f10.5)
      end
