c---------------------------------------------------------
c     Physical Property of U-V System for 3D
c---------------------------------------------------------
      subroutine prop3duv(ng3d,n2uq
     &                   ,vres,urlj,listcore,cr,tr)
c
c     icl           ... closure type 0...HNC, 1...MSA, 2...KH
c     ngrid3d       ... number of grid of 3D-RDF
c     nv            ... number of site of solvent
c     cr            ... direct correlation function 
c     tr            ... tau bond =hr-cr
c     
      implicit real*8 (a-h,o-z)
      complex*16 cr

      include "phys_const.i"
      include "rismio.i"
      include "rismrun.i"
      include "solvent.i"
      include "solute.i"

      dimension cr(ng3d,n2uq),tr(ng3d,n2uq)
      dimension vres(ng3d)
      dimension urlj(ng3d,n2uq)
      dimension listcore(ng3d)
      dimension esolvi(n2uq),egfi(n2uq)
c
      namelist /RISMUC/ahnc,akh,agf,bhnc,bkh,bgf
c
      data ahnc,akh,agf,bhnc,bkh,bgf
     &     /-4.58,-4.58,-3.16,0.340,0.340,0.894/
      ! Ref JPCB, 2015, 119, 5588

      rd33=rdelta3d**3
c---------------------------------------------------------
      call vclr(esolvi,1,nvuq)
      call vclr(egfi,1,nvuq)
c
c     --- Solvation Free Energy (Excess Chemical Potential)
c
c
c     --- HNC
c
      if (icl.eq.0.or.icl.eq.3) then

         sum=0.d0

         do j=1,nvuq

         sumi=0.d0
         
!$OMP PARALLEL DO PRIVATE(K,JJ,CRR,HR) REDUCTION(+: SUMI)
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d

            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d*ngrid3d

            crr=dble(cr(k,j))
            hr=tr(k,j)+crr
            sumi=sumi+(-crr+0.5d0*hr*tr(k,j))*densuq(j)

         enddo
         enddo
         enddo
!$OMP END PARALLEL DO

         sum=sum+sumi
         esolvi(j)=sumi/beta*rd33         

         enddo

         esolvtot=sum/beta*rd33
         
      endif
c
c     --- MSA
c
      if (icl.eq.1) then

         sum=0.d0
         do j=1,nvuq
         sumi=0.d0

!$OMP PARALLEL DO PRIVATE(K,JJ,CRR,HR) REDUCTION(+: SUMI)
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d

            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d*ngrid3d

            crr=dble(cr(k,j))
            hr=tr(k,j)+crr
            sumi=sumi+(-crr-0.5d0*hr*crr)*densuq(j)

         enddo
         enddo
         enddo
!$OMP END PARALLEL DO

         sum=sum+sumi
         esolvi(j)=sumi/beta*rd33

         enddo
         esolvtot=sum/beta*rd33
         
      endif
c
c     --- KH
c
      if (icl.eq.2) then

         sum=0.d0

         do j=1,nvuq

         sumi=0.d0

!$OMP PARALLEL DO PRIVATE(K,JJ,CRR,HR) REDUCTION(+: SUMI)
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d
            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d*ngrid3d

            crr=dble(cr(k,j))
            hr=tr(k,j)+crr
            if (hr.ge.0.d0) then 
               dhevi=0.d0
            else
               dhevi=1.d0
            endif

            sumi=sumi
     &           +(-crr-0.5d0*hr*crr+0.5d0*hr**2*dhevi)
     &           *densuq(j)
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
         
         sum=sum+sumi 
         esolvi(j)=sumi/beta*rd33

         enddo

         esolvtot=sum/beta*rd33
         
      endif
c
c     --- GF
c
      sum=0.d0

      do j=1,nvuq

         sumi=0.d0

!$OMP PARALLEL DO PRIVATE(K,JJ,CRR,HR) REDUCTION(+: SUMi)
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d
            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d*ngrid3d

            crr=dble(cr(k,j))
            hr=tr(k,j)+crr
            sumi=sumi+(-crr-0.5d0*hr*crr)
     &           *densuq(j)

         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
      
         sum=sum+sumi
         egfi(j)=sumi/beta*rd33
      enddo

      egftot=sum/beta*rd33
c     
c     ----- Solute-Solvent Binding Energy
c     
      ebtot=0.d0
      eblj=0.d0
!$OMP PARALLEL DO PRIVATE(JJ,CRR,GR) REDUCTION(+:ebtot,eblj)
      do k=1,ngrid3d**3
         if (listcore(k).eq.0) goto 8000

         do j=1,nvuq

            crr=dble(cr(k,j))
            gr=tr(k,j)+crr+1.d0
            if (gr.lt.0.d0) gr=0.d0
            ebtot=ebtot+densuq(j)
     &           *gr*(vres(k)*q2uq(j)+urlj(k,j))*rd33
            eblj=eblj+densuq(j)
     &           *gr*(urlj(k,j))*rd33
         enddo
 8000    continue
      enddo
!$OMP END PARALLEL DO
      ebindtot=ebtot
      ebindlj=eblj
c
c     
c     ----- Partial Molar Volume
c
      ck0=0.d0
      do j=1,nvuq
      sumi=0.d0
!$OMP PARALLEL DO PRIVATE(K,JJ) REDUCTION(+: SUMi)
      do kz=1,ngrid3d
      do ky=1,ngrid3d
      do kx=1,ngrid3d

         k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d*ngrid3d

         sumi=sumi+dble(cr(k,j))*densuq(j)
         
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
      ck0=ck0+sumi*rd33
      enddo

      pmv=xt*(1.d0-ck0)/beta    ![m^3/mol]
c
c     --- Total solvent density
c
      totaldens=0.d0
      do i=1,nvuq
         totaldens=totaldens+densuq(i)
      enddo
      totaldens=totaldens*1.D+30 ![molecule/m^3]
c
c     --- Universal correction for excess chemical potential
c
      ir=45
      open (ir,file=inpfile,status='old')
      
      rewind ir
      read (ir,RISMUC,end=9000)
 9000 continue

      close(ir)
      
      pmvx=pmv*totaldens/avognum        ! Dimensionless PMV 

      if (icl.eq.0.or.icl.eq.3) then ! for HNC
         auc=ahnc
         buc=bhnc
      elseif (icl.eq.2) then    ! for KH
         auc=akh
         buc=bkh
      endif
      ucterm=(auc*pmvx+buc)*1.d+3 ! for HNC or KH [J/mol]
      uctermgf=(agf*pmvx+bgf)*1.d+3 ! for GF      [J/mol]

      esolvuc=esolvtot/ekcal+ucterm ! for HNC or KH
      esolvgfuc=egftot/ekcal+uctermgf ! for GF  
c
c     --- Pressure correction for excess chemical potential
c
      pressure=0.5d0/beta*(beta/xt+totaldens/avognum) 
                                                 ![Pa]=[J/m^3]
      pcterm=-pressure*pmv ![J/mol]
c
c     --- Total Solvent Charge
c
      chgtot=0.d0
      do i=1,nvuq
         sum=0.d0
         do ig=1,ngrid3d**3
            sum=sum+tr(ig,i)+dble(cr(ig,i))+1.d0
         enddo
         chgtot=chgtot+sum*rd33*densuq(i)
     &        *q2uq(i)
      enddo 

      write(*,'(/,4x,A21,f16.8)') "Total Solvent Charge:",chgtot

c---------------------------------------------------------
c     Print Property of U-V System for 3D
c---------------------------------------------------------
      write(*,9999) 
      write(*,9998) esolvtot
      write(*,9993) egftot
      write(*,9997) ebindtot
      write(*,9996) ebindlj
      write(*,9995) (ebindtot-ebindlj)
      if (icl.eq.3) write(*,9994)
c
      write(*,9988) esolvuc,ucterm,auc,buc
     &             ,esolvgfuc,uctermgf*1.d-3,agf,bgf
c
      write(*,9987) (pcterm+esolvtot),pcterm
c
      write(*,9992)
      do i=1,nvuq
         write(*,9991) i,esolvi(i)
      enddo
      write(*,9990)
      do i=1,nvuq
         write(*,9991) i,egfi(i)
      enddo

      write(*,9989) pmv*1.d+3      ! from [m^3/mol] to [L/mol]
     &             ,xt*1.d+9       ! from [/Pa] to [/GPa]
     &             ,pressure/1.d+9 ! from [Pa] to [GPa]
c---------------------------------------------------------
c     Ourput .xmu file
c---------------------------------------------------------
      ift=45
      open (ift,file=trim(basename)//".xmu")
      write(ift,'(1X,A7)') "$RESULT"
      write(ift,'(A7,F16.5,2x,A8)') "SFE_SC=",esolvtot,"!(J/mol)"
      do i=1,nvuq
         write(ift,'(A8,i3,A2,f16.5)') "SFEC_SC(",i-1,")=",esolvi(i)
      enddo
      write(ift,'(A7,F16.5,2x,A8)') "SFE_GF=",egftot,"!(J/mol)"
      do i=1,nvuq
         write(ift,'(A8,i3,A2,f16.5)') "SFEC_GF(",i-1,")=",egfi(i)
      enddo

      write(ift,'(A4,f16.5,2x,A8)') "PMV=",pmv*1.d+3,"!(L/mol)"
      write(ift,'(A9,f16.5,2x,A11)') "Pressure=",pressure,"!(Pa=L/m^3)"
      write(ift,'(A16,f16.5,2x,A8)') 
     &     "Correction_term=",pcterm,"!(J/mol)"

      write(ift,'(1X,A4)') "$END"
      close(ift)
c---------------------------------------------------------
      return
c---------------------------------------------------------
 9987 format (/,4x,"Solvation Free Energy(PC)     :",g16.8,"[J/mol]"
     &         ,4x,"PC term :",g16.8,"[J/mol]")
 9988 format (/,4x,"Solvation Free Energy(UC)     :",g16.8,"[J/mol]"
     &         ,4x,"UC term :",g16.8,"[J/mol]"
     &         ,4x,"UC param: A=",g16.8
     &            ,"[J/mol], B=",g16.8,"[J/mol]"
     &        /,4x,"Solvation Free Energy(GF-UC)  :",g16.8,"[J/mol]"
     &         ,4x,"UC term :",g16.8,"[J/mol]"
     &         ,4x,"UC param: A=",g16.8
     &            ,"[J/mol], B=",g16.8,"[J/mol]")
 9989 format (/,4x,"============= Partial Molar Volume ============="
     &       ,/,4x,"Partial Molar Volume       = ",G12.4," [L/mol]"
     &       ,/,4x,"Isothermal Compressibility = ",G12.4," [/GPa]"
     &       ,/,4x,"Pressure                   = ",G12.4," [GPa]")
 9990 format (/,4x,"======= GF Solvation Free Energy Component =====")
 9991 format (  4x,i4,":",g16.8,"[kJ/mol]")
 9992 format (/,4x,"======= Solvation Free Energy Component  =======")
 9993 format (/,4x,"Solvation Free Energy(GF)     :",g16.8,"[J/mol]")
 9994 format (/,4x,"----------------------------------------------"
     &       ,/,4x,"NOTE:When HNC+RBC closure are selected,       "  
     &       ,/,4x,"     Solvation Free Energy is not correct.    "  
     &       ,/,4x,"     This value is evaluated based on HNC.    "  
     &       ,/,4x,"----------------------------------------------")
 9995 format (/,4x,"UV Binding Energy Component ES:",g16.8,"[J/mol]")
 9996 format (/,4x,"UV Binding Energy Component LJ:",g16.8,"[J/mol]")
 9997 format (/,4x,"Solute-Solvent Binding Energy :",g16.8,"[J/mol]")
 9998 format (/,4x,"Solvation Free Energy         :",g16.8,"[J/mol]")
 9999 format (/,4x,"======= Physical Property of U-V System =======")
      end
