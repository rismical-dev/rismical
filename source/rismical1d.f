c**************************************************************
c--------------------------------------------------------------
c     1D-RISM for solute-solvent systems
c--------------------------------------------------------------
c**************************************************************
      subroutine rismical1d
c
      implicit real*8 (a-h,o-z)
      real*8 ,allocatable :: cr(:,:,:),tr(:,:,:)
     &                      ,hvk(:,:,:),fr(:,:,:),fk(:,:,:),huvk(:,:,:)
     &                      ,ures(:,:,:),urlj(:,:,:)
     &                      ,wk1(:,:,:),wk2(:,:,:),ck(:,:,:)
     &                      ,rbc(:,:,:),rbc2(:,:,:)
     &                      ,zrk(:,:,:),wk2org(:,:,:)
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
      call readuvdata
c
c     --- allocate array
c
      allocate (cr(ngrid,nu,nv),tr(ngrid,nu,nv))
      allocate (hvk(ngrid,nv,nv),fr(ngrid,nu,nv),fk(ngrid,nu,nv))
      allocate (huvk(ngrid,nu,nv))
      allocate (ures(ngrid,nu,nv),urlj(ngrid,nu,nv))
      allocate (wk1(ngrid,nu,nu),wk2(ngrid,nv,nv),ck(ngrid,nu,nv))
      allocate (wk2org(ngrid,nv,nv))
      allocate (zrk(ngrid,nu,nv))
c     
c     --- Initialize
c     
      call vclr(hvk,1,nv*nv*ngrid)
c
      idrism=0
c     
c     --- Intramolecular correlation function
c     
      call makewxv(nv,ngrid,rdelta,wk2)
      call makewxu(nu,ngrid,rdelta,wk1)
c     
c     --- Make Potential and F-bond
c     
      call potentialuv(ngrid,nu,nv,rdelta,ures,urlj)
      
      call fbonduv(ngrid,nu,nv,rdelta,fr,fk,alp1d)
c
c     --- Read solvent suseptibility
c
      call readhvk(ngrid,nv,hvk)
c
c     --- Make initial guess for tr(r)
c
      frfac=beta*chgratio

      if (iguess.eq.0) then
         do j=1,nv
            do i=1,nu
               do k=1,ngrid
                  tr(k,i,j)=frfac*fr(k,i,j)
                  cr(k,i,j)=-tr(k,i,j)
               enddo
            enddo
         enddo
      elseif (iguess.eq.1) then
         call readguess1d(ngrid,nu,nv,tr)
      else
         write(*,9990)
         ierr=4
         CALL ABRT(ierr)
      endif
c---------------------------------------------------------
c     Charge Up Procedure
c---------------------------------------------------------
      chgratio=chgratio-chgstep
 1000 continue
      
      write(*,9995)
c     
c     --- Charge up
c     
      chgratio=chgratio+chgstep
      if (chgratio.gt.1.d0) chgratio=1.d0
      write(*,9999) temp,chgratio
c     
c     --- Make initial guess for tr(r) in charge up cycle
c     
      prefac=frfac
      frfac=beta*chgratio
      do j=1,nv
         do i=1,nu
            do k=1,ngrid
               tr(k,i,j)=tr(k,i,j)-prefac*fr(k,i,j)
     &                            +frfac*fr(k,i,j)
            enddo
         enddo
      enddo
c     
c     --- Set Convergence Criterion
c     
      if (chgratio.eq.1.d0) chgconv=conv
      cconv=chgconv
c     
c     --- Setup mdiis 
c     
      ng=ngrid*nu*nv
      call mdiis(ng,tr,residu,cconv,0)  
c---------------------------------------------------------
c     RISM Iteration Cycle
c---------------------------------------------------------
      write(*,9991)
      do itr=1,itrmax
c     
c     --- Closure - SSOZ [INPUT tr(r) -> OUTPUT trnew(r)]
c     
         iuv=1
         call cl_oz1d(icl,iuv,idrism,ngrid,rdelta,nu,nv
     &               ,chgratio,ck,hvk,fr,fk,wk1,wk2,zrk
     &               ,cr,tr,ures,urlj)
c     
c     --- Check Convergence and Make Guess For Next Loop
c     
         call mdiis(ng,tr,residu,cconv,1)  
         if (residu.le.cconv) goto 8000
         
      enddo
c---------------------------------------------------------
c     Not Converged
c---------------------------------------------------------
      write(*,9997)
      CALL ABRT
c---------------------------------------------------------
c     Converged
c---------------------------------------------------------
 8000 continue

      call  cl_1d(icl,ngrid,rdelta,nu,nv,chgratio
     &     ,cr,tr,ures,urlj)

      write(*,9989) itr,residu,0,"X"
      write(*,9993)
c
c     --- Go to next charge up cycle
c
      if (chgratio.ne.1.d0) goto 1000
c
 8100 continue
c---------------------------------------------------------
c     Output 1d-rism result
c---------------------------------------------------------
 9000 continue
c
      call prop1duv(icl,ngrid,rdelta,nu,nv
     &     ,cr,tr,ures,urlj)

      call output1d(ngrid,rdelta,nu,nv
     &     ,cr,tr,ures,urlj,fr,ck,fk)

c---------------------------------------------------------
      deallocate (cr,tr)
      deallocate (hvk,fr,fk,huvk)
      deallocate (ures,urlj)
      deallocate (wk1,wk2,ck)
      deallocate (zrk)
c---------------------------------------------------------
      return
c---------------------------------------------------------
 9989 format (4x,i6,f20.12,2x,i4,2x,a1,4x,"CONVERGED")
 9990 format (/,4x,"Error, Wrong Parameter For IGUESS IN $RISMGUESS")
 9991 format (/,3x,"ITR",6x,"RESIDUAL",7x,"#-SUB",1x,"MIN",6x,"DUMP")
 9993 format (/,4x,"RISM CYCLE IS CONVERGED")
 9995 format (/,4x,"========== U-V 1DRISM ==========")
 9997 format (/,4x,"RISM CYCLE IS NOT CONVERGED.",
     &        /,4x,"----- E N D -----")
 9999 format (/,4x,"Temperature      :",f10.5,"[K]",
     &        /,4x,"Charge Up Factor :",f10.5)
      end
