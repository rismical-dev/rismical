c------------------------------------------------------------
c     Physical Property of U-V System for 1D
c------------------------------------------------------------
      subroutine prop1duv(icl,ngrid,rdelta,n1,n2,
     &                    cr,tr,ures,urlj)
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
      include "rismio.i"
      include "solvent.i"
      include "solute.i"

      dimension cr(ngrid,n1,n2),tr(ngrid,n1,n2)
      dimension esolv(n1),cn(n1,n2),cnd(n1,n2),ebind(n1),egf(n1)
      dimension elec(n1)
      dimension ures(ngrid,n1,n2),urlj(ngrid,n1,n2)
      dimension esoluv(n1,n2)

c----------------------------------------------------------------
      do i=1,n1
         do j=1,n2
            esoluv(i,j)=0.d0
         enddo
      enddo
c
c     --- Solvation Free Energy (Excess Chemical Potential)
c
c
c     --- HNC
c
      if (icl.eq.0) then
         esolvtot=0.d0
         do i=1,n1
            esolv(i)=0.d0
            sum=0.d0
            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               do j=1,n2
                  hr=tr(k,i,j)+cr(k,i,j)
                  sum=sum+rr*(-cr(k,i,j)+0.5d0*hr*tr(k,i,j))
     &                 *dens(nspc(j))
                  esoluv(i,j)=esoluv(i,j)
     &                 +rr*(-cr(k,i,j)+0.5d0*hr*tr(k,i,j))
     &                 *densuq(j)*4.d0*pi/beta
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
         do i=1,n1
            esolv(i)=0.d0
            sum=0.d0
            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               do j=1,n2
                  hr=tr(k,i,j)+cr(k,i,j)
                  sum=sum+rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j))
     &                 *dens(nspc(j))
                  esoluv(i,j)=esoluv(i,j)
     &                 +rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j))
     &                 *densuq(j)*4.d0*pi/beta
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
         do i=1,n1
            esolv(i)=0.d0
            sum=0.d0
            do j=1,n2
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
     &                 *densuq(j)
                  esoluv(i,j)=esoluv(i,j)
     &                 +rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j)
     &                 +0.5d0*hr**2*dhevi)
     &                 *densuq(j)*4.d0*pi/beta
               enddo
            enddo
            esolv(i)=sum*4.d0*pi/beta
            esolvtot=esolvtot+esolv(i)
         enddo
      endif
c     
c     --- HNC+RBC (Same as HNC case, this is not correct.)
c     
      if (icl.eq.3.or.icl.eq.4) then
         
         write(*,*) "RBC is not implemented yet."
         ierr=30
         call abrt(ierr)

      endif
c     
c     --- GF (Gaussian Fluctuation  Same as MSA expression)
c     
      egftot=0.d0
      do i=1,n1
         egf(i)=0.d0
         sum=0.d0
         do k=1,ngrid
            rr=(rdelta*dble(k))**2*rdelta
            do j=1,n2
               hr=tr(k,i,j)+cr(k,i,j)
               sum=sum+rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j))
     &                 *densuq(j)
            enddo
         enddo
         egf(i)=sum*4.d0*pi/beta
         egftot=egftot+egf(i)
      enddo
c
c     --- Coordination Number
c
      do i=1,n1
         do j=1,n2
            nmin=0
            iflag=0
            do 100 k=2,ngrid
               gr=tr(k,i,j)+cr(k,i,j)+1.d0
               grpre=tr(k-1,i,j)+cr(k-1,i,j)+1.d0
               if (gr.lt.grpre) iflag=1
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
            cn(i,j)=sum*4.d0*pi*densuq(j)
         enddo
      enddo
c
c     --- Solute-Solvent Binding Energy
c
      ebindtot=0.d0
      do i=1,n1
         ebind(i)=0.d0
         sum=0.d0
         do j=1,n2

            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               gr=tr(k,i,j)+cr(k,i,j)+1.d0
               sum=sum+rr*gr*(ures(k,i,j)+urlj(k,i,j))
     &              *densuq(j)
            enddo
         enddo
         ebind(i)=sum*4.d0*pi
         ebindtot=ebindtot+ebind(i)
      enddo
c
c     --- Electrostatic interaction between U and V
c
      electot=0.d0
      do i=1,n1
         elec(i)=0.d0
         sum=0.d0
         do j=1,n2
            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               gr=tr(k,i,j)+cr(k,i,j)+1.d0
               sum=sum+rr*gr*ures(k,i,j)
     &              *densuq(j)
            enddo
         enddo
         elec(i)=sum*4.d0*pi
         electot=electot+elec(i)
      enddo
c     
c     --- Partial molar volume (PMV)
c     
      ck0=0.d0
      do j=1,n2
         do i=1,n1
            sum=0.d0
            do k=1,ngrid
               rr=(rdelta*dble(k))**2
               sum=sum+rr*cr(k,i,j)*rdelta
     &              *densuq(j)
            enddo
            ck0=ck0+sum
         enddo
      enddo
      ck0=ck0*4.d0*pi
      pmv=xt*(1.d0-ck0)/beta
c     
c     --- Total Solvent Charge
c     
      do i=1,n1
         chgtot=0.d0
         do j=1,n2
            do k=1,ngrid
               rr=(rdelta*dble(k))**2
               gr=tr(k,i,j)+cr(k,i,j)+1.d0
               chgtot=chgtot+4.d0*pi*rr*gr*rdelta
     &              *qv(j)*densuq(j)
            enddo
         enddo
      enddo

c----------------------------------------------------------------
c     Print Property of U-V System for 1D
c----------------------------------------------------------------
      write(*,9999) 
c
c     --- Excess Chemical Potential
c      
      write(*,9998) esolvtot
      do i=1,n1
         write(*,9997) nsiteu(i),esolv(i)
      enddo
c
c     --- GF Formula Excess Chemical Potential
c      
      write(*,9992) egftot
      do i=1,n1
         write(*,9997) nsiteu(i),egf(i)
      enddo
c
c     --- U-V Binding Energy
c
      write(*,9994) ebindtot
      do i=1,n1
         write(*,9993) nsiteu(i),ebind(i)
      enddo
c
c     --- Electrostatic interaction
c
      write(*,9984) electot
      do i=1,n1
         write(*,9983) nsiteu(i),elec(i)*1.d-3
      enddo
c
c     --- Coordination Number
c
      do i=1,n2,10
         imax=min(i+9,n2)
         write(*,9996) (nsitev(k),k=i,imax)
         do j=1,n1
            write(*,9995) nsiteu(j),(cn(j,k),cnd(j,k),k=i,imax)
         enddo
      enddo
C
c
c     --- Excess Chemical Potential Components
c
      do i=1,n2,10
         imax=min(i+9,n2)
         write(*,9986) (nsitev(j),j=i,imax)
         do j=1,n1
            write(*,9985) nsiteu(j),(esoluv(j,k),k=i,imax)
         enddo
      enddo
c
c     --- Partial Molar Volume and Compressibility
c
      write(*,9989) pmv*1.d+3 ! from [m^3/mol] to [L/mol]
     &              ,xt*1.d+9 ! from [/Pa] to [/GPa]
c
c---------------------------------------------------------
c     Ourput .xmu file
c---------------------------------------------------------
      ift=45
      open (ift,file=trim(basename)//".xmu")
      write(ift,'(1X,A8)') "$RESULT"
      write(ift,'(A7,F16.5,A8)') "SFE_SC=",esolvtot,"!(J/mol)"
      do i=1,nvuq
         write(ift,'(A8,i3,A2,f16.5)') "SFEC_SC(",i-1,")=",esolv(i)
      enddo
      write(ift,'(A7,F16.5,A8)') "SFE_GF=",egftot,"!(J/mol)"
      do i=1,nvuq
         write(ift,'(A8,i3,A2,f16.5)') "SFEC_GF(",i-1,")=",egf(i)
      enddo

      write(ift,'(A4,f16.5,A8)') "PMV=",pmv*1.d+3,"!(L/mol)"

      write(ift,'(1X,A4)') "$END"
      close(ift)
c----------------------------------------------------------------
      return
c----------------------------------------------------------------
 9983 format (  6x,"              ----->",a4,":",e16.8,"[J/mol]")
 9984 format (/,4x,"Electrostatic Intraction  :",e16.8,"[J/mol]",
     &        /,5x,"Site Contribution (On Solute Site)")
 9985 format (4x,A8,1x,10F16.5)
 9986 format (/,4x,"Excess Chemical Potential Components [J/mol]",
     &     /,4x,"U\V",4x,10A16)
 9989 format (/,4x,"============= Partial Molar Volume ============="
     &       ,/,4x,"Partial Molar Volume       = ",G16.8," [L/mol]"
     &       ,/,4x,"Isothermal Compressibility = ",G16.8," [/GPa]")
 9992 format (/,4x,"GF Form Ex.Chem.Pot.      :",f16.8,"[J/mol]",
     &        /,5x,"Site Contribution (On Solute Site)")
 9993 format (  6x,"              ----->",a4,":",f16.8,"[J/mol]")
 9994 format (/,4x,"U-V Binding Energy        :",f16.8,"[J/mol]",
     &        /,5x,"Site Contribution (On Solute Site)")
 9995 format (5x,a4,1x,10(f8.5,"[",f7.3,"]"))
 9996 format (/,4x,"Coordination Number ",
     &        /,5x,"U\V ",1x,10(A8,"[Radius ]"))
 9997 format (  6x,"              ----->",a4,":",f16.8,"[J/mol]") 
 9998 format (/,4x,"Excess Chemical Potential :",f16.8,"[J/mol]",
     &        /,4x,"                  RBC/TPT :",f16.8,"[J/mol]",
     &        /,4x,"  Ex.Chem.Pot. +  RBC/TPT :",f16.8,"[J/mol]",
     &        /,5x,"Site Contribution (On Solute Site)")
 9999 format (/,4x,"======= Physical Property of U-V System =======")
      end
c------------------------------------------------------------
c     Physical Property of U-V System for 1D
c     -- Non-reduced solvent site version --
c------------------------------------------------------------
      subroutine prop1duv_old(icl,ngrid,rdelta,n1,n2,
     &                    cr,tr,ures,urlj)
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
      include "solute.i"

      dimension cr(ngrid,n1,n2),tr(ngrid,n1,n2)
      dimension esolv(n1),cn(n1,n2),cnd(n1,n2),ebind(n1),egf(n1)
      dimension elec(n1)
      dimension ures(ngrid,n1,n2),urlj(ngrid,n1,n2)
      dimension esoluv(n1,n2)

c----------------------------------------------------------------
      do i=1,n1
         do j=1,n2
            esoluv(i,j)=0.d0
         enddo
      enddo
c
c     --- Solvation Free Energy (Excess Chemical Potential)
c
c
c     --- HNC
c
      if (icl.eq.0) then
         esolvtot=0.d0
         do i=1,n1
            esolv(i)=0.d0
            sum=0.d0
            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               do j=1,n2
                  hr=tr(k,i,j)+cr(k,i,j)
                  sum=sum+rr*(-cr(k,i,j)+0.5d0*hr*tr(k,i,j))
     &                 *dens(nspc(j))
                  esoluv(i,j)=esoluv(i,j)
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
         do i=1,n1
            esolv(i)=0.d0
            sum=0.d0
            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               do j=1,n2
                  hr=tr(k,i,j)+cr(k,i,j)
                  sum=sum+rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j))
     &                 *dens(nspc(j))
                  esoluv(i,j)=esoluv(i,j)
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
         do i=1,n1
            esolv(i)=0.d0
            sum=0.d0
            do j=1,n2
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
                  esoluv(i,j)=esoluv(i,j)
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
c     --- HNC+RBC (Same as HNC case, this is not correct.)
c     
      if (icl.eq.3.or.icl.eq.4) then
         
         write(*,*) "RBC is not implemented yet."
         ierr=30
         call abrt(ierr)

      endif
c     
c     --- GF (Gaussian Fluctuation  Same as MSA expression)
c     
      egftot=0.d0
      do i=1,n1
         egf(i)=0.d0
         sum=0.d0
         do k=1,ngrid
            rr=(rdelta*dble(k))**2*rdelta
            do j=1,n2
               hr=tr(k,i,j)+cr(k,i,j)
               sum=sum+rr*(-cr(k,i,j)-0.5d0*hr*cr(k,i,j))
     &                 *dens(nspc(j))
            enddo
         enddo
         egf(i)=sum*4.d0*pi/beta
         egftot=egftot+egf(i)
      enddo
c
c     --- Coordination Number
c
      do i=1,n1
         do j=1,n2
            nmin=0
            iflag=0
            do 100 k=2,ngrid
               gr=tr(k,i,j)+cr(k,i,j)+1.d0
               grpre=tr(k-1,i,j)+cr(k-1,i,j)+1.d0
               if (gr.lt.grpre) iflag=1
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
c     --- Solute-Solvent Binding Energy
c
      ebindtot=0.d0
      do i=1,n1
         ebind(i)=0.d0
         sum=0.d0
         do j=1,n2

            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               gr=tr(k,i,j)+cr(k,i,j)+1.d0
               sum=sum+rr*gr*(ures(k,i,j)+urlj(k,i,j))
     &              *dens(nspc(j))
            enddo
         enddo
         ebind(i)=sum*4.d0*pi
         ebindtot=ebindtot+ebind(i)
      enddo
c
c     --- Electrostatic interaction between U and V
c
      electot=0.d0
      do i=1,n1
         elec(i)=0.d0
         sum=0.d0
         do j=1,n2
            do k=1,ngrid
               rr=(rdelta*dble(k))**2*rdelta
               gr=tr(k,i,j)+cr(k,i,j)+1.d0
               sum=sum+rr*gr*ures(k,i,j)
     &              *dens(nspc(j))
            enddo
         enddo
         elec(i)=sum*4.d0*pi
         electot=electot+elec(i)
      enddo
c     
c     --- Partial molar volume (PMV)
c     
      ck0=0.d0
      do j=1,n2
         do i=1,n1
            sum=0.d0
            do k=1,ngrid
               rr=(rdelta*dble(k))**2
               sum=sum+rr*cr(k,i,j)*rdelta
     &              *dens(nspc(j))
            enddo
            ck0=ck0+sum
         enddo
      enddo
      ck0=ck0*4.d0*pi
      pmv=xt*(1.d0-ck0)/beta
c     
c     --- Total Solvent Charge
c     
      do i=1,n1
         chgtot=0.d0
         do j=1,n2
            do k=1,ngrid
               rr=(rdelta*dble(k))**2
               gr=tr(k,i,j)+cr(k,i,j)+1.d0
               chgtot=chgtot+4.d0*pi*rr*gr*rdelta
     &              *qv(j)*dens(nspc(j))
            enddo
         enddo
      enddo

c----------------------------------------------------------------
c     Print Property of U-V System for 1D
c----------------------------------------------------------------
      write(*,9999) 
c
c     --- Excess Chemical Potential
c      
      write(*,9998) esolvtot*1.d-3
      do i=1,n1
         write(*,9997) nsiteu(i),esolv(i)*1.d-3
      enddo
c
c     --- GF Formula Excess Chemical Potential
c      
      write(*,9992) egftot*1.d-3
      do i=1,n1
         write(*,9997) nsiteu(i),egf(i)*1.d-3
      enddo
c
c     --- U-V Binding Energy
c
      write(*,9994) ebindtot*1.d-3
      do i=1,n1
         write(*,9993) nsiteu(i),ebind(i)*1.d-3
      enddo
c
c     --- Electrostatic interaction
c
      write(*,9984) electot*1.d-3
      do i=1,n1
         write(*,9983) nsiteu(i),elec(i)*1.d-3
      enddo
c
c     --- Coordination Number
c
      do i=1,n2,10
         imax=min(i+9,n2)
         write(*,9996) (nsitev(k),k=i,imax)
         do j=1,n1
            write(*,9995) nsiteu(j),(cn(j,k),cnd(j,k),k=i,imax)
         enddo
      enddo
C
c
c     --- Excess Chemical Potential Components
c
      do i=1,n2,10
         imax=min(i+9,n2)
         write(*,9986) (nsitev(j),j=i,imax)
         do j=1,n1
            write(*,9985) nsiteu(j),(esoluv(j,k)*1.d-3,k=i,imax)
         enddo
      enddo
c
c     --- Partial Molar Volume and Compressibility
c
      write(*,9989) pmv*1.d+3 ! from [m^3/mol] to [L/mol]
     &              ,xt*1.d+9 ! from [/Pa] to [/GPa]
c
c----------------------------------------------------------------
      return
c----------------------------------------------------------------
 9983 format (  6x,"              ----->",a4,":",e16.8,"[kJ/mol]")
 9984 format (/,4x,"Electrostatic Intraction  :",e16.8,"[kJ/mol]",
     &        /,5x,"Site Contribution (On Solute Site)")
 9985 format (4x,A8,1x,10F16.5)
 9986 format (/,4x,"Excess Chemical Potential Components [kJ/mol]",
     &     /,4x,"U\V",4x,10A16)
 9989 format (/,4x,"============= Partial Molar Volume ============="
     &       ,/,4x,"Partial Molar Volume       = ",G16.8," [L/mol]"
     &       ,/,4x,"Isothermal Compressibility = ",G16.8," [/GPa]")
 9992 format (/,4x,"GF Form Ex.Chem.Pot.      :",f16.8,"[kJ/mol]",
     &        /,5x,"Site Contribution (On Solute Site)")
 9993 format (  6x,"              ----->",a4,":",f16.8,"[kJ/mol]")
 9994 format (/,4x,"U-V Binding Energy        :",f16.8,"[kJ/mol]",
     &        /,5x,"Site Contribution (On Solute Site)")
 9995 format (5x,a4,1x,10(f8.5,"[",f7.3,"]"))
 9996 format (/,4x,"Coordination Number ",
     &        /,5x,"U\V ",1x,10(A8,"[Radius ]"))
 9997 format (  6x,"              ----->",a4,":",f16.8,"[kJ/mol]") 
 9998 format (/,4x,"Excess Chemical Potential :",f16.8,"[kJ/mol]",
     &        /,4x,"                  RBC/TPT :",f16.8,"[kJ/mol]",
     &        /,4x,"  Ex.Chem.Pot. +  RBC/TPT :",f16.8,"[kJ/mol]",
     &        /,5x,"Site Contribution (On Solute Site)")
 9999 format (/,4x,"======= Physical Property of U-V System =======")
      end
