c**************************************************************
c--------------------------------------------------------------
c     1D-RISM for solute-solvent systems
c--------------------------------------------------------------
c**************************************************************
      subroutine rismical1d
c
      implicit real*8 (a-h,o-z)
      real*8 ,allocatable :: cr(:,:,:),tr(:,:,:)
     &                      ,xvk(:,:,:),fr(:,:,:),fk(:,:,:),huvk(:,:,:)
     &                      ,ures(:,:,:),urlj(:,:,:)
     &                      ,wk1(:,:,:),ck(:,:,:)
     &                      ,zrk(:,:,:)
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
      allocate (cr(ngrid,nu,nvuq),tr(ngrid,nu,nvuq))
      allocate (fr(ngrid,nu,nvuq),fk(ngrid,nu,nvuq))
      allocate (xvk(ngrid,nvuq,nvuq))
      allocate (huvk(ngrid,nu,nvuq))
      allocate (ures(ngrid,nu,nvuq),urlj(ngrid,nu,nvuq))
      allocate (wk1(ngrid,nu,nu),ck(ngrid,nu,nvuq))
      allocate (zrk(ngrid,nu,nvuq))
c     
c     --- Initialize
c     
      call vclr(xvk,1,nvuq*nvuq*ngrid)
c
      idrism=0
c     
c     --- Intramolecular correlation function
c     
      call makewxu(nu,ngrid,rdelta,wk1)
c     
c     --- Make Potential and F-bond
c     
      call potentialuv(ngrid,nu,nvuq,rdelta,ures,urlj)
      
      call fbonduv(ngrid,nu,nvuq,rdelta,fr,fk,alp1d)
c
c     --- Read solvent suseptibility
c
      call readxvk(ngrid,nvuq,xvk)
c
c     --- Make initial guess for tr(r)
c
      frfac=beta*chgratio

      if (iguess.eq.0) then
         do j=1,nvuq
            do i=1,nu
               do k=1,ngrid
                  tr(k,i,j)=frfac*fr(k,i,j)
                  cr(k,i,j)=-tr(k,i,j)
               enddo
            enddo
         enddo
      elseif (iguess.eq.1) then
         call readguess1d(ngrid,nu,nvuq,tr)
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
      do j=1,nvuq
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
      ng=ngrid*nu*nvuq
      call mdiis(ng,tr,residu,cconv,0)  
c---------------------------------------------------------
c     RISM Iteration Cycle
c---------------------------------------------------------
      write(*,9991)
      do itr=1,itrmax
c     
c     --- Closure - SSOZ [INPUT tr(r) -> OUTPUT trnew(r)]
c     
         call cl_oz1duv(icl,ngrid,rdelta,nu,nvuq
     &               ,chgratio,ck,xvk,fr,fk,wk1,zrk
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
      ierr=425
      CALL ABRT(ierr)
c---------------------------------------------------------
c     Converged
c---------------------------------------------------------
 8000 continue

      call  cl_1d(icl,ngrid,rdelta,nu,nvuq,chgratio
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
      call prop1duv(icl,ngrid,rdelta,nu,nvuq
     &     ,cr,tr,ures,urlj)
c$$$      call prop1duv_old(icl,ngrid,rdelta,nu,nvuq
c$$$     &     ,cr,tr,ures,urlj)

      call output1d(ngrid,rdelta,nu,nvuq
     &     ,cr,tr,ures,urlj,fr,ck,fk)

c---------------------------------------------------------
      deallocate (cr,tr)
      deallocate (xvk,fr,fk,huvk)
      deallocate (ures,urlj)
      deallocate (wk1,ck)
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
