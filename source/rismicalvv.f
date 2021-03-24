c**************************************************************
c--------------------------------------------------------------
c     1D-RISM for solvent-solvent systems
c--------------------------------------------------------------
c**************************************************************
      subroutine rismicalvv
c
      implicit real*8 (a-h,o-z)
      real*8 ,allocatable :: cr(:,:,:),tr(:,:,:)
     &                      ,hvk(:,:,:),fr(:,:,:),fk(:,:,:),huvk(:,:,:)
     &                      ,ures(:,:,:),urlj(:,:,:)
     &                      ,wk1(:,:,:),wk2(:,:,:),ck(:,:,:)
     &                      ,zrk(:,:,:),wk2org(:,:,:),xvk(:,:,:)
c
      include "phys_const.i"
      include "solvent.i"
      include "rismrun.i"
c
c----------------------------------------------------------------------
c
c     --- get RISM run parameters
c
      call readinput
c      
c     --- get solvent parameters
c
      call readvvdata
c
c     --- allocate array
c
      allocate (cr(ngrid,nv,nv),tr(ngrid,nv,nv))
      allocate (hvk(ngrid,nv,nv),fr(ngrid,nv,nv),fk(ngrid,nv,nv))
      allocate (huvk(ngrid,nv,nv))
      allocate (xvk(ngrid,nv,nv))
      allocate (ures(ngrid,nv,nv),urlj(ngrid,nv,nv))
      allocate (wk1(ngrid,nv,nv),wk2(ngrid,nv,nv),ck(ngrid,nv,nv))
      allocate (wk2org(ngrid,nv,nv))
      allocate (zrk(ngrid,nv,nv))
c     
c     --- Initialize
c     
      call vclr(xvk,1,nv*nv*ngrid)
c
      idrism=0
c
c     --- Intramolecular correlation function of solvent
c
      call makewxv(nv,ngrid,rdelta,wk2)
c
c     --- Setup DRISM
c
      call drismzeta(idrism,nv,ngrid,rdelta,wk2,zrk)
c     
c     --- Make Potential and F-bond
c     
      call potentialvv(ngrid,nv,rdelta,ures,urlj)
      
      call fbondvv(ngrid,nv,rdelta,fr,fk,alp1d)
c
c     --- Make initial guess for tr(r)
c
      frfac=beta*chgratio

      if (iguess.eq.0) then
         do i=1,nv
            do j=1,nv
               do k=1,ngrid
                  tr(k,i,j)=frfac*fr(k,i,j)
                  cr(k,i,j)=-tr(k,i,j)
               enddo
            enddo
         enddo
      elseif (iguess.eq.1) then
         call readguess1d(ngrid,nv,nv,tr)
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
      
      write(*,9994)
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
      do i=1,nv
         do j=1,nv
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
      ng=ngrid*nv*nv
      call mdiis(ng,tr,residu,cconv,0)  

c---------------------------------------------------------
c     RISM Iteration Cycle
c---------------------------------------------------------
      write(*,9991)
      do itr=1,itrmax
c     
c     --- Closure - SSOZ [INPUT tr(r) -> OUTPUT tr(r)]
c     
         iuv=0
         call cl_oz1d(icl,iuv,idrism,ngrid,rdelta,nv,nv
     &               ,chgratio,ck,hvk,fr,fk,wk2,wk2,zrk
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

      call  cl_1d(icl,ngrid,rdelta,nv,nv,chgratio
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
      call calxvk(ngrid,nv,hvk,xvk,wk2)

      if (idrism.eq.1) then
         do iv2=1,nv
         do iv1=1,nv
            do k=1,ngrid
               wk2(k,iv1,iv2)=wk2(k,iv1,iv2)
     &              -dens(nspc(iv1))*zrk(k,iv1,iv2)
            enddo
         enddo
        enddo
      endif

      call prop1dvv(icl,ngrid,rdelta,rcore,nv
     &     ,wk2,ck,cr,tr,ures,urlj)

      call outputvv(ngrid,rdelta,nv
     &     ,cr,tr,hvk,xvk,ures,urlj,fr,ck,fk)

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
 9994 format (/,4x,"========== V-V 1DRISM ==========")
 9995 format (/,4x,"========== U-V 1DRISM ==========")
 9997 format (/,4x,"RISM CYCLE IS NOT CONVERGED.",
     &        /,4x,"----- E N D -----")
 9999 format (/,4x,"Temperature      :",f10.5,"[K]",
     &        /,4x,"Charge Up Factor :",f10.5)
      end
