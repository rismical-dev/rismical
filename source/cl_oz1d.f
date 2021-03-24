c--------------------------------------------------------------
c     Closure-OZ for VV
c--------------------------------------------------------------
      subroutine cl_oz1dvv(icl,idrism,ngrid,rdelta,n1,n2
     &     ,chgratio,ck,hvk,fr,fk,wk1,wk2,zrk
     &     ,cr,tr,ures,urlj)
c     
c     ir        ... file number of output (STDOUT)
c     iw        ... file number of input (STDIN)
c     icl       ... closure type 0...HNC, 1...MSA, 2...KH 3...RBC
c     ngrid     ... number of grid of RDF
c     rdelta    ... grid width of r-space [Angstrom]
c     n1        ... number of site of 1 (solute or solvent)
c     n2        ... number of site of 2 (solvent)
c     chgratio  ... charge up ratio
c     cr        ... direct correlation function 
c     tr        ... tau bond =hr-cr
c     ures      ... electro static potential [erg]
c     urlj      ... LJ potential [erg]
c     beta      ... inverse of kbT[1/erg]
c     rbc       ... bridge correction term exp(b(r))
c     
      implicit real*8 (a-h,o-z)

      include "phys_const.i"
      include "solvent.i"

      dimension ures(ngrid,n1,n2),urlj(ngrid,n1,n2)
      dimension cr(ngrid,n1,n2),ck(ngrid,n1,n2)
      dimension tr(ngrid,n1,n2)
      dimension hvk(ngrid,n2,n2)
      dimension fr(ngrid,n1,n2),fk(ngrid,n1,n2)
      dimension wk1(ngrid,n1,n1),wk2(ngrid,n2,n2)
      dimension tk(ngrid,n1,n2),ftfunc(ngrid)
      dimension dum1(n1,n2),dum2(n1,n2),dum3(n1,n2)
      dimension wktemp1(n1,n1),wktemp2(n2,n2),dumvv(n2,n2)
      dimension list(n2)
      dimension zrk(ngrid,n2,n2)
c----------------------------------------------------------------
      do j=1,n2
         do i=1,n1
c     
c     --- HNC
c     
            if (icl.eq.0) then
               do k=1,ngrid
                  bur=beta*(ures(k,i,j)*chgratio+urlj(k,i,j))
                  cr(k,i,j)=dexp(-bur+tr(k,i,j))
     &                 -tr(k,i,j)-1.d0
               enddo
c     
c     --- MSA
c     
            elseif (icl.eq.1) then
               do k=1,ngrid
                  bur=beta*(ures(k,i,j)*chgratio+urlj(k,i,j))
                  cr(k,i,j)=-bur
               enddo
c     
c     --- KH
c     
            elseif (icl.eq.2) then
               do k=1,ngrid
                  
                  bur=beta*(ures(k,i,j)*chgratio+urlj(k,i,j))
                  d=-bur+tr(k,i,j)
                  
                  if (d.gt.0.d0) then
                     cr(k,i,j)=d-tr(k,i,j)
                  else
                     cr(k,i,j)=dexp(d)-tr(k,i,j)-1.d0
                  endif
                  
               enddo
c     
c     --- HNC+RBC or CB-RBC
c     
            elseif (icl.eq.3) then

               write(*,*) "Error. RBC is not implemented yet."
               ierr=6
               call abrt(ierr)

            endif
            
         enddo
      enddo
c
c     ---------------------------------------
c
      do j=1,n2
         do i=1,n1
c
c     --- Create short range part
c
            do k=1,ngrid
               ftfunc(k)=cr(k,i,j)+beta*fr(k,i,j)*chgratio
            enddo
c     
c     --- Fourier transform 
c     
            call fft1d(ngrid,rdelta,ftfunc,1)
c     
c     --- Add long range part to short range part 
c     
            do k=1,ngrid
               ck(k,i,j)=ftfunc(k)-beta*fk(k,i,j)*chgratio
            enddo

         enddo
      enddo
c     
c     --- k-space VV-OZ
c     
      do k=1,ngrid
            
         if (n2.ge.2) then
            
            do i=1,n2
               do j=1,n2
                  dum1(i,j)=ck(k,i,j)
                  wktemp2(i,j)=wk2(k,i,j)
               enddo
            enddo
c
            do i=1,n2
               do j=1,n2
                  sum=0.d0
                  do kk=1,n2
                     sum=sum+wktemp2(kk,j)*dum1(kk,i)
                  enddo
                  dum2(j,i)=sum
               enddo
            enddo
c
            do i=1,n2
               do j=1,n2
                  dum3(i,j)=-dum2(i,j)*dens(nspc(j))
               enddo
               dum3(i,i)=1.d0+dum3(i,i)
            enddo
c      
            call matprd(dum2,wktemp2,dum1,n2,n2,n2,n2,n2,n2,ill)

            eps=1.d-16
            ind=0
            det=0.d0

            call matsolv(dum3,dum1,n2,n2,n2)

            do i=1,n2
               do j=1,n2
                  hvk(k,i,j)=dum1(i,j)
                  tk(k,i,j) =dum1(i,j)-ck(k,i,j)
               enddo
            enddo

         else
               
            hvk(k,1,1)=ck(k,1,1)/(1.d0-dens(nspc(1))*ck(k,1,1))
            tk(k,1,1)=hvk(k,1,1)-ck(k,1,1)

         endif
            
      enddo
c     
      if (idrism.eq.1) then
         do j=1,n2
            do i=1,n2
               do k=1,ngrid
                  hvk(k,i,j)=hvk(k,i,j)+zrk(k,i,j)
                  tk(k,i,j) =tk(k,i,j) +zrk(k,i,j)
               enddo
            enddo
         enddo
      endif
c     
      do j=1,n2
         do i=1,n1
c     
c     --- Create short range part
c     
            do k=1,ngrid
               ftfunc(k)=tk(k,i,j)-beta*fk(k,i,j)*chgratio
            enddo
c     
c     --- Reverse Fourier transform 
c     
            call fft1d(ngrid,rdelta,ftfunc,-1)
c     
c     --- Add long range part to short range part 
c     
            do k=1,ngrid
               tr(k,i,j)=ftfunc(k)+beta*fr(k,i,j)*chgratio
            enddo
            
         enddo
      enddo
c--------------------------------------------------------------
      return
      end
c**************************************************************
c--------------------------------------------------------------
c     Closure
c--------------------------------------------------------------
      subroutine cl_1d(icl,ngrid,rdelta,n1,n2,chgratio
     &     ,cr,tr,ures,urlj)
c     
c     ir        ... file number of output (STDOUT)
c     iw        ... file number of input (STDIN)
c     icl       ... closure type 0...HNC, 1...MSA, 2...KH 3...RBC
c     ngrid     ... number of grid of RDF
c     rdelta    ... grid width of r-space [Angstrom]
c     n1        ... number of site of 1 (solute or solvent)
c     n2        ... number of site of 2 (solvent)
c     chgratio  ... charge up ratio
c     cr        ... direct correlation function 
c     tr        ... tau bond =hr-cr
c     ures      ... electro static potential [erg]
c     urlj      ... LJ potential [erg]
c     beta      ... inverse of kbT[1/erg]
c     rbc       ... bridge correction term exp(b(r))
c     
      implicit real*8 (a-h,o-z)

      include "phys_const.i"
      include "solvent.i"

      dimension ures(ngrid,n1,n2),urlj(ngrid,n1,n2)
      dimension cr(ngrid,n1,n2)
      dimension tr(ngrid,n1,n2)
c----------------------------------------------------------------
      do j=1,n2
         do i=1,n1
c     
c     --- HNC
c     
            if (icl.eq.0) then
               do k=1,ngrid
                  bur=beta*(ures(k,i,j)*chgratio+urlj(k,i,j))
                  cr(k,i,j)=dexp(-bur+tr(k,i,j))
     &                 -tr(k,i,j)-1.d0
               enddo
c     
c     --- MSA
c     
            elseif (icl.eq.1) then
               do k=1,ngrid
                  bur=beta*(ures(k,i,j)*chgratio+urlj(k,i,j))
                  cr(k,i,j)=-bur
               enddo
c     
c     --- KH
c     
            elseif (icl.eq.2) then
               do k=1,ngrid
                  
                  bur=beta*(ures(k,i,j)*chgratio+urlj(k,i,j))
                  d=-bur+tr(k,i,j)
                  
                  if (d.gt.0.d0) then
                     cr(k,i,j)=d-tr(k,i,j)
                  else
                     cr(k,i,j)=dexp(d)-tr(k,i,j)-1.d0
                  endif
                  
               enddo
c     
c     --- HNC+RBC or CB-RBC
c     
            elseif (icl.eq.3) then

               write(*,*) "Error. RBC is not implemented yet."
               ierr=6
               call abrt(ierr)

            endif
            
         enddo
      enddo
c
c     ---------------------------------------
c
c--------------------------------------------------------------
      return
      end
c**************************************************************
c--------------------------------------------------------------
c     Closure-OZ for UV
c--------------------------------------------------------------
      subroutine cl_oz1duv(icl,ngrid,rdelta,n1,n2
     &     ,chgratio,ck,xvk,fr,fk,wk1,zrk
     &     ,cr,tr,ures,urlj)
c     
c     icl       ... closure type 0...HNC, 1...MSA, 2...KH 3...RBC
c     ngrid     ... number of grid of RDF
c     rdelta    ... grid width of r-space [Angstrom]
c     n1        ... number of site of 1 (solute or solvent)
c     n2        ... number of site of 2 (solvent)
c     chgratio  ... charge up ratio
c     cr        ... direct correlation function 
c     tr        ... tau bond =hr-cr
c     ures      ... electro static potential [erg]
c     urlj      ... LJ potential [erg]
c     beta      ... inverse of kbT[1/erg]
c     rbc       ... bridge correction term exp(b(r))
c     
      implicit real*8 (a-h,o-z)

      include "phys_const.i"
      include "solvent.i"

      dimension ures(ngrid,n1,n2),urlj(ngrid,n1,n2)
      dimension cr(ngrid,n1,n2),ck(ngrid,n1,n2)
      dimension tr(ngrid,n1,n2)
      dimension xvk(ngrid,n2,n2)
      dimension fr(ngrid,n1,n2),fk(ngrid,n1,n2)
      dimension wk1(ngrid,n1,n1)
      dimension tk(ngrid,n1,n2),ftfunc(ngrid)
      dimension dum1(n1,n2),dum2(n1,n2),dum3(n1,n2)
      dimension wktemp1(n1,n1),wktemp2(n2,n2),dumvv(n2,n2)
      dimension list(n2)
      dimension zrk(ngrid,n2,n2)
c----------------------------------------------------------------
      do j=1,n2
         do i=1,n1
c     
c     --- HNC
c     
            if (icl.eq.0) then
               do k=1,ngrid
                  bur=beta*(ures(k,i,j)*chgratio+urlj(k,i,j))
                  cr(k,i,j)=dexp(-bur+tr(k,i,j))
     &                 -tr(k,i,j)-1.d0
               enddo
c     
c     --- MSA
c     
            elseif (icl.eq.1) then
               do k=1,ngrid
                  bur=beta*(ures(k,i,j)*chgratio+urlj(k,i,j))
                  cr(k,i,j)=-bur
               enddo
c     
c     --- KH
c     
            elseif (icl.eq.2) then
               do k=1,ngrid
                  
                  bur=beta*(ures(k,i,j)*chgratio+urlj(k,i,j))
                  d=-bur+tr(k,i,j)
                  
                  if (d.gt.0.d0) then
                     cr(k,i,j)=d-tr(k,i,j)
                  else
                     cr(k,i,j)=dexp(d)-tr(k,i,j)-1.d0
                  endif
                  
               enddo
c     
c     --- HNC+RBC or CB-RBC
c     
            elseif (icl.eq.3) then

               write(*,*) "Error. RBC is not implemented yet."
               ierr=6
               call abrt(ierr)

            endif
            
         enddo
      enddo
c
c     ---------------------------------------
c
      do j=1,n2
         do i=1,n1
c
c     --- Create short range part
c
            do k=1,ngrid
               ftfunc(k)=cr(k,i,j)+beta*fr(k,i,j)*chgratio
            enddo
c     
c     --- Fourier transform 
c     
            call fft1d(ngrid,rdelta,ftfunc,1)
c     
c     --- Add long range part to short range part 
c     
            do k=1,ngrid
               ck(k,i,j)=ftfunc(k)-beta*fk(k,i,j)*chgratio
            enddo

         enddo
      enddo
c     
c     --- k-space UV-OZ
c     
      do k=1,ngrid
c     
         do i=1,n1
            do j=1,n2
               dum1(i,j)=ck(k,i,j)
            enddo
         enddo
         do i=1,n1
            do i2=1,n1
               wktemp1(i,i2)=wk1(k,i,i2)
            enddo
         enddo
c     
         call matprd(wktemp1,dum1,dum2,n1,n1,n2,n1,n1,n2,ill)
            
         do j=1,n2
            do j2=1,n2
               dumvv(j,j2)=xvk(k,j,j2)
            enddo
         enddo
            
         call matprd(dum2,dumvv,dum3,n1,n2,n2,n1,n2,n2,ill)
         
         do i=1,n1
            do j=1,n2
               tk(k,i,j)=dum3(i,j)-ck(k,i,j)
            enddo
         enddo
c     
      enddo
c     
      do j=1,n2
         do i=1,n1
c     
c     --- Create short range part
c     
            do k=1,ngrid
               ftfunc(k)=tk(k,i,j)-beta*fk(k,i,j)*chgratio
            enddo
c     
c     --- Reverse Fourier transform 
c     
            call fft1d(ngrid,rdelta,ftfunc,-1)
c     
c     --- Add long range part to short range part 
c     
            do k=1,ngrid
               tr(k,i,j)=ftfunc(k)+beta*fr(k,i,j)*chgratio
            enddo
            
         enddo
      enddo
c--------------------------------------------------------------
      return
      end
