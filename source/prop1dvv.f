c------------------------------------------------------------
c     Physical Property of V-V System for 1D
c------------------------------------------------------------
      subroutine prop1dvv(icl,ngrid,rdelta,rcore,n
     &                    ,wk,ck,cr,tr,ures,urlj)
c
c     icl           ... closure type 0...HNC, 1...MSA, 2...KH 3...RBC
c     ngrid         ... number of grid of RDF
c     rdelta        ... grid width of r-space
c     rcore         ... ratio of core diameter (less than 1)
c     n             ... number of site of 2 (solvent)
c     dens          ... number density[molecule/angstrom3]
c     cr            ... direct correlation function 
c     ck            ... k-space direct correlation function 
c     tr            ... tau bond =hr-cr
c     wk            ... intra-molecular correlation function
c     ures          ... electro static potential [erg]
c     urlj          ... LJ potential [erg]
c     beta          ... inverse of kbT [/erg]
c     siglj         ... sigma of LJ parameter [Angstrom]
c     q             ... partial charge of site [e]
c     
      implicit real*8 (a-h,o-z)

      include "phys_const.i"
      include "solvent.i"

      dimension cr(ngrid,n,n),tr(ngrid,n,n),ck(ngrid,n,n)
      dimension ures(ngrid,n,n),urlj(ngrid,n,n)
      dimension esolv(n),cn(n,n),cnd(n,n),ebind(n),egf(n)
      dimension ud(ngrid)
      dimension wk(ngrid,n,n)
      dimension dum1(n,n),dum2(n,n),dum3(n,n)
      dimension esolvv(n,n)
      dimension esolspc(n)

      deltak=pi/(dble(ngrid)*rdelta)
c---------------------------------------------------------
      do i=1,n
         do j=1,n
            esolvv(i,j)=0.d0
         enddo
         esolspc(i)=0.d0
      enddo
c
c     --- Excess Helmholtz Free Energy per Particle
c
      esolvtot=0.d0
      do i=1,n
         esolv(i)=0.d0
      enddo
c
      if (numspc.eq.1) then
c
c
c     --- Excess Helmholtz Free Energy per Particle for pure solvent
c
c     --- 1st term (0.5*hr**2*(heviside)-cr)
c
c
c     --- HNC
c      
      if (icl.eq.0) then
         do i=1,n
            sum=0.d0
            do j=1,n
               do k=1,ngrid
                  rr=(rdelta*dble(k))**2*rdelta
                  hr=tr(k,i,j)+cr(k,i,j)
                  sum=sum+rr*(0.5d0*hr**2-cr(k,i,j))
               enddo
            enddo
            sum=sum*4.d0*pi*dens(1)/beta
            esolvtot=esolvtot+sum
         enddo
      endif
c
c     --- MSA
c      
      if (icl.eq.1) then
         do i=1,n
            sum=0.d0
            do j=1,n
               do k=1,ngrid
                  rr=(rdelta*dble(k))**2*rdelta
                  sum=sum+rr*(-cr(k,i,j))
               enddo
            enddo
            sum=sum*4.d0*pi*dens(1)/beta
            esolvtot=esolvtot+sum
         enddo
      endif
c
c     --- KH
c      
      if (icl.eq.2) then
         do i=1,n
            sum=0.d0
            do j=1,n
               do k=1,ngrid
                  rr=(rdelta*dble(k))**2*rdelta
                  hr=tr(k,i,j)+cr(k,i,j)
                  hevi=1.d0
                  if (hr.gt.0.d0) hevi=0.d0
                  sum=sum+rr*(0.5d0*hr**2*hevi-cr(k,i,j))
               enddo
            enddo
            sum=sum*4.d0*pi*dens(1)/beta
            esolvtot=esolvtot+sum
         enddo
      endif
c
c     --- GF
c      
      egftot=0.d0
      do i=1,n
         sum=0.d0
         do j=1,n
            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               sum=sum+rr*(-cr(k,i,j))
            enddo
         enddo
         sum=sum*4.d0*pi*dens(1)/beta
         egftot=egftot+sum
      enddo
c
c     --- 2nd term (Tr[wc]+lndet[I-wc*rho])
c
      ierr=0
      sumk=0.d0
      do k=1,ngrid

         do i=1,n
            do j=1,n
               dum1(i,j)=wk(k,i,j)
               dum2(i,j)=ck(k,i,j)
            enddo
         enddo

         call matprd(dum1,dum2,dum3,n,n,n,n,n,n,ill)

         do i=1,n
            do j=1,n
               dum2(i,j)=-dum3(i,j)*dens(1)
            enddo
            dum2(i,i)=1.d0+dum2(i,i)
         enddo
         
         call matdet(dum2,n,n,d)

         if (d.le.0.d0) then
            ierr=1
            goto 90
         endif
         
         call trd(t,dum3,n)
         
         rk2=(deltak*dble(k))**2
         sumk=sumk+rk2*(t+dlog(d)/dens(1))

      enddo
      
      esolvtot=esolvtot+sumk/(4.d0*pi**2*beta)*deltak
      egftot=egftot+sumk/(4.d0*pi**2*beta)*deltak
 90   continue
c
      else 
c
c     --- Solvation Free Energy (Excess Chemical Potential) for mixture
c
c
c     --- HNC
c
      if (icl.eq.0) then
         esolvtot=0.d0
         do i=1,n
            esolv(i)=0.d0
            sum=0.d0
            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               do j=1,n
                  hr=tr(k,i,j)+cr(k,i,j)
                  sum=sum+rr*(-cr(k,i,j)+0.5d0*hr*tr(k,i,j))
     &                 *dens(nspc(j))
                  esolvv(i,j)=esolvv(i,j)
     &                 +rr*(-cr(k,i,j)+0.5d0*hr*tr(k,i,j))
     &                 *dens(nspc(j))*4.d0*pi/beta
               enddo
            enddo
            esolv(i)=sum*4.d0*pi/beta
            esolvtot=esolvtot+esolv(i)
         enddo
      endif
c
c     --- MSA
c
      if (icl.eq.1) then
         esolvtot=0.d0
         do i=1,n
            esolv(i)=0.d0
            sum=0.d0
            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               do j=1,n
                  hr=tr(k,i,j)+cr(k,i,j)
                  sum=sum+rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j))
     &                 *dens(nspc(j))
                  esolvv(i,j)=esolvv(i,j)
     &                 +rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j))
     &                 *dens(nspc(j))*4.d0*pi/beta
               enddo
            enddo
            esolv(i)=sum*4.d0*pi/beta
            esolvtot=esolvtot+esolv(i)
         enddo
      endif
c
c     --- KH
c
      if (icl.eq.2) then
         esolvtot=0.d0
         do i=1,n
            esolv(i)=0.d0
            sum=0.d0
            do j=1,n
               do k=1,ngrid
                  rr=(rdelta*dble(k))**2*rdelta
                  hr=tr(k,i,j)+cr(k,i,j)
                  if (hr.ge.0.d0) then 
                     dhevi=0.d0
                  else
                     dhevi=1.d0
                  endif
                  sum=sum+rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j)
     &                 +0.5d0*hr**2*dhevi)
     &                 *dens(nspc(j))
                  esolvv(i,j)=esolvv(i,j)
     &                 +rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j)
     &                 +0.5d0*hr**2*dhevi)
     &                 *dens(nspc(j))*4.d0*pi/beta
               enddo
            enddo
            esolv(i)=sum*4.d0*pi/beta
            esolvtot=esolvtot+esolv(i)
         enddo
      endif
c
c     --- GF (Gaussian Fluctuation  Same as MSA expression)
c
      egftot=0.d0
      do i=1,n
         egf(i)=0.d0
         sum=0.d0
         do k=1,ngrid
            rr=(rdelta*dble(k))**2*rdelta
            do j=1,n
               hr=tr(k,i,j)+cr(k,i,j)
               sum=sum+rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j))
     &                 *dens(nspc(j))
            enddo
         enddo
         egf(i)=sum*4.d0*pi/beta
         egftot=egftot+egf(i)
      enddo
c
      endif
c
c     --- Coordination Number
c
      do i=1,n
         do j=1,n
            nmin=0
            iflag=0
            do 100 k=2,ngrid
               gr   =tr(k  ,i,j)+cr(k  ,i,j)+1.d0
               grpre=tr(k-1,i,j)+cr(k-1,i,j)+1.d0
               if  (gr.lt.grpre) iflag=1
               if ((gr.gt.grpre).and.(iflag.eq.1)) then
                  nmin=k-1
                  goto 110
               endif
 100        continue
 110        continue
            sum=0.d0
            do 120 k=1,nmin-1
               rr1=(rdelta*dble(k))**2
               rr2=(rdelta*dble(k+1))**2
               gr1=tr(k,i,j)+cr(k,i,j)+1.d0
               gr2=tr(k+1,i,j)+cr(k+1,i,j)+1.d0
               sum=sum+(gr1*rr1+gr2*rr2)*rdelta*0.5d0
 120        continue
            cnd(i,j)=dble(nmin)*rdelta
            cn(i,j)=sum*4.d0*pi*dens(nspc(j))
         enddo
      enddo
c
      pvir=0.d0
      pfree=0.d0
      comp=0.d0
      dfluc=0.d0
c
c     --- Isothermal Compressibility and Density Fluctuation
c
      ck0=0.d0
      do i=1,n
         do j=1,n
            sum=0.d0
            do k=1,ngrid
               rr=(rdelta*dble(k))**2
               sum=sum+rr*cr(k,i,j)*rdelta
     &              *dens(nspc(i))*dens(nspc(j))
            enddo
            ck0=ck0+sum
         enddo
      enddo
      ck0=ck0*4.d0*pi
      
      totaldens=0.d0
      ndum=0
      do i=1,n
         if (ndum.ne.nspc(i)) then
            totaldens=totaldens+dens(nspc(i))
            ndum=nspc(i)
         endif
      enddo
      dfluc=totaldens/(totaldens-ck0)
c      comp=beta*dfluc/totaldens*101325.d0*1.d-23*1.d+5
      comp=dfluc/totaldens*beta*1d-30*avognum ![/Pa]
      xt=comp   ! for record in common block
c
      if (numspc.eq.1) then
c     
c     --- Pressure (From Virial Equation : P-route)
c
c     NOTE: Pressure from Virial Equation in site-site
c           representation is not correct.
c     
      pvir=0.d0
      do i=1,n
         do j=1,n
            do k=1,ngrid
               ud(k)=ures(k,i,j)+urlj(k,i,j)
            enddo
            call differential(ud,rdelta,ngrid,ngrid)
            sum=0.d0
            do k=1,ngrid
               gr=tr(k,i,j)+cr(k,i,j)+1.d0
               rr3=(rdelta*dble(k))**3
               sum=sum+rr3*gr*ud(k)
            enddo
            pvir=pvir+sum
         enddo
      enddo
      pvir=(1.d0-4.d0*pi*dens(1)*beta/6.d0*pvir*rdelta)
     &     *1.d+23/101325.d0*dens(1)/beta
c
c     --- Pressure (From Free Energy Derivative : A-route)
c     
      pfree=dens(1)/beta /(1d-30*avognum) ![Pa]
c
c     --- 1st term (0.5*hr**2-cr)
c
      sumr=0.d0
      do i=1,n
         do j=1,n
            do k=1,ngrid
               rr2=(rdelta*dble(k))**2
               hr=tr(k,i,j)+cr(k,i,j)
               if ((icl.eq.0).or.(icl.eq.2.and.hr.le.0.d0)) then
                  sumr=sumr+rr2*0.5d0*hr**2
               endif
               sumr=sumr-rr2*cr(k,i,j)
            enddo
         enddo
      enddo
      pfree=pfree+2.d0*pi*dens(1)**2/beta*sumr*rdelta /(1d-30*avognum) ![Pa]
c
c     --- 2nd term (ln[det(1-wcp)]-Tr[wcp(1-wcp)^-1]
c
      ierr=0
      sumk=0.d0
      do k=1,ngrid
         
         do i=1,n
            do j=1,n
               dum1(i,j)=wk(k,i,j)
               dum2(i,j)=ck(k,i,j)*dens(1)
            enddo
         enddo
         
         call matprd(dum1,dum2,dum3,n,n,n,n,n,n,ill)

         do i=1,n
            do j=1,n
               dum1(i,j)=-dum3(i,j)
               dum2(i,j)=-dum3(i,j)
            enddo
            dum1(i,i)=1.d0+dum1(i,i)
            dum2(i,i)=1.d0+dum2(i,i)
         enddo

         call matdet(dum2,n,n,d)

         if (d.le.0.d0) then
            ierr=1
            pfree=0.d0
            goto 140
         endif
         
         call matinv(dum1,n,n)

         call matprd(dum3,dum1,dum2,n,n,n,n,n,n,ill)
         call trd(t,dum2,n)
         
         rk2=(deltak*dble(k))**2
         sumk=sumk+rk2*(dlog(d)+t)

      enddo
      
      pfree=pfree-sumk/(4.d0*pi**2*beta)*deltak/(1d-30*avognum) ![Pa]

 140  continue
c
      endif
c
c     --- Solvent-Solvent Binding Energy
c
      ebindtot=0.d0
      do i=1,n
         ebind(i)=0.d0
         sum=0.d0
         do j=1,n
            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               gr=tr(k,i,j)+cr(k,i,j)+1.d0
               sum=sum+rr*gr*(ures(k,i,j)+urlj(k,i,j))*dens(nspc(j))
            enddo
         enddo
         ebind(i)=sum*4.d0*pi
         ebindtot=ebindtot+ebind(i)
      enddo
c---------------------------------------------------------
c     Print Property of V-V System for 1D
c---------------------------------------------------------
C
      write(*,9999)
c
c     --- Excess Chemical Potential
c      
      write(*,9998) esolvtot
c
c     --- GF Formula Excess Chemical Potential
c      
      write(*,9990) egftot
c
c     --- V-V Binding Energy
c
      write(*,9992) ebindtot
      do i=1,n
         write(*,9991) nsitev(i),ebind(i)
      enddo
c
c     --- Coordination Number
c
      do i=1,n,10
         imax=min(i+9,n)
         write(*,9996) (nsitev(j),j=i,imax)
         do j=1,n
            write(*,9995) nsitev(j),(cn(j,k),cnd(j,k),k=i,imax)
         enddo
      enddo
c
c
c     --- Isothermal Compressibility
c
      write(*,9994) dfluc,comp*1.d+9,comp/9.86923d-6
c
c
      if (numspc.eq.1) then
c
c     --- Pressure
c
         write(*,9993) pfree,pfree*9.86923d-6
c
      endif
c
c     --- Excess Chemical Potential Components
c
      if (numspc.gt.1) then
         do i=1,n,10
            imax=min(i+9,n)
            write(*,9986) (nsitev(j),j=i,imax)
            do j=1,n
               write(*,9985) nsitev(j),(esolvv(j,k) ,k=i,imax)
                           
            enddo
         enddo

         do i=1,n
            do j=1,n
               esolspc(nspc(i))=esolspc(nspc(i))+esolvv(i,j)
            enddo
         enddo
         
         write(*,9984)
         do i=1,numspc
            write(*,9983) i,esolspc(i)
         enddo

      endif

c---------------------------------------------------------
      return
c---------------------------------------------------------
 9983 format (4x,i8,4x,F20.9)
 9984 format (/,4x,"Excess Chemical Potential on Each Species",
     &        /,4x,"# of Spc",14x,"[J/mol]")
 9985 format (4x,A8,1x,10F16.5)
 9986 format (/,4x,"Excess Chemical Potential Components [J/mol]",
     &        /,13x,10A16)
 9990 format (/,4x,"GF Form Ex.Chem.Pot.     :",f16.8,"[J/mol]")
 9991 format (  6x,"             ----->",a4,":",f16.8,"[J/mol]")
 9992 format (/,4x,"V-V Binding Energy       :",f16.8,"[J/mol]",
     &        /,5x,"Site Contribution (On Solvent Site)")
 9993 format (/,4x,"Pressure (A route)         :",e16.8,"[Pa]"
     &        /,4x,"                           :",e16.8,"[atm]")
 9994 format (/,4x,"Density Fluctuation        :",f16.8,
     &        /,4x,"Isothermal Compressibility :",f16.8,"[/GPa]",
     &        /,4x,"                           :",f16.8,"[/atm]")
 9995 format (5x,a4,1x,10(f8.5,"[",f7.3,"]"))
 9996 format (/,4x,"Coordination Number ",
     &        /,5x,"NAME",1x,10(A8,"[Radius ]"))
 9997 format (  6x,"             ----->",a4,":",f16.8,"[J/mol]")
 9998 format (/,4x,"Excess Chemical Potential:",f16.8,"[J/mol]")
 9999 format (/,4x,"======= Physical Property of V-V System =======")
      end
