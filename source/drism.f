c**************************************************************
c----------------------------------------------------------------
c     Make Zeta and Modify wk2 for VV-DRISM
c----------------------------------------------------------------
      subroutine drismzeta(idrism,n2,ngrid,rdelta
     &                    ,wk2,zrk)
c
c     beta [J/mol]
c
      implicit real*8 (a-h,o-z)
      character*4 char4
      character*6 char6
c
      include "phys_const.i"
      include "solvent.i"
c
      dimension wk2(ngrid,n2,n2)
      dimension zrk(ngrid,n2,n2)
      dimension d0x(n2),d0y(n2),d1z(n2)
      dimension dipspc(100)
      dimension rv(3,n2),q2(n2)
c
c---------------------------------------------------------------
c
c     ---- Read NAMELIST DRISM
c
      call drismio(idrism,delec,sparam)
c
      call vclr(zrk,1,ngrid*n2*n2)
c
c     Skip if not a drism run
c
      if (idrism.eq.0) goto 9000
c
c     set array of solvent charge and coordinate
c
      do i=1,n2
         q2(i)=qv(i)
         do ixyz=1,3
            rv(ixyz,i)=xyzv(ixyz,i)
         enddo
      enddo
c
c     --- Shift to center of absolute charge
c     (See Perkyns and Pettitt, JCP 97. 7656. 1992)
c
      cofchgx=0.d0
      cofchgy=0.d0
      cofchgz=0.d0
      qtotabs=0.d0
      do i=1,n2
         qtotabs=qtotabs+abs(q2(i))
         cofchgx=cofchgx+abs(q2(i))*rv(1,i)
         cofchgy=cofchgy+abs(q2(i))*rv(2,i)
         cofchgz=cofchgz+abs(q2(i))*rv(3,i)
      enddo
      if (qtotabs.ne.0.d0) then
         cofchgx=cofchgx/qtotabs
         cofchgy=cofchgy/qtotabs
         cofchgz=cofchgz/qtotabs
         r=cofchgx**2+cofchgy**2+cofchgz**2
         r=dsqrt(r)
         do i=1,n2
            rv(1,i)=rv(1,i)-cofchgx
            rv(2,i)=rv(2,i)-cofchgy
            rv(3,i)=rv(3,i)-cofchgz
         enddo
      endif
c
c     --- Rotate solvent 
c
      call rotate2zaxis(rv,q2,n2)
c
      do i=1,100
         dipspc(i)=0.d0
      enddo
      ispc=0
      jspc=0
      do i=1,n2
         if (jspc.ne.nspc(i)) ispc=ispc+1
         dipspc(ispc)=dipspc(ispc)+q2(i)*rv(3,i)
         jspc=nspc(i)
      enddo
c
      dipdens=0.d0
      densdip2=0.d0
      do jspc=1,ispc
         dipdens=dipdens+dens(jspc)
         densdip2=densdip2+dens(jspc)*dipspc(jspc)**2
      enddo
c
      y=4.d0*pi*beta*densdip2/9.d0
     &     /angtobohr*hart2jmol
c     &     /(ergtoau*angtobohr)
c
      if (dabs(dipdens).lt.1.d-10) then
         hcfac=0.d0
         write(*,9999) dipdens
      else
         hcfac=((delec-1.d0)/y-3.d0)/dipdens
      endif
c
c     ---- Make zeta and Modify wk2
c
      do k=1,ngrid
        rk = dble(k)*pi/(dble(ngrid)*rdelta)
        hck = hcfac * dexp(-(sparam*rk/2.d0)**2)   
c
c     ---- calculating spherical Bessel functions
c
        do iv=1,n2

          trx = rk*rv(1,iv)
          if (trx.eq.0.d0)  then
            d0x(iv) = 1.d0           
          else
            d0x(iv) = dsin(trx)/trx
          endif

          try = rk*rv(2,iv)
          if (try.eq.0.d0)  then
            d0y(iv) = 1.d0 
          else
            d0y(iv) = dsin(try)/try
          endif                  

          trz = rk*rv(3,iv)
          if (trz.eq.0.d0)  then
            d1z(iv) = 0.d0
          else 
            d1z(iv) = dsin(trz)/trz**2 - dcos(trz)/trz
          endif
        enddo 
c
c     ---- getting Zvv(k)
c
        do iv2=1,n2
        do iv1=1,n2
           zrk(k,iv1,iv2) = d0x(iv1)*d0y(iv1)*d1z(iv1) * hck
     &                    * d0x(iv2)*d0y(iv2)*d1z(iv2)
        enddo                           
        enddo
c
c     ---- getting Wvv(k)+Zvv(k)
c
        do iv2=1,n2
        do iv1=1,n2
          wk2(k,iv1,iv2) = wk2(k,iv1,iv2)
     &          + dens(nspc(iv1))*zrk(k,iv1,iv2)
        enddo
        enddo

      enddo
c----------------------------------------------------------------
 9000 continue
c      
      return
 9500 continue
      write(ir,*) "ERROR, $VDATA is missing."
      call abrt
 9998 format (/, "Error! Center of charge of solvent molecule is",
     &           " too far from molecular origin."
     &        /, "Center of charge of solvent molecule is :"
     &            3f15.5)
 9999 format (/,4x,"-----------------------------------------",
     &        /,4x,"WARNING:",
     &        /,4x,"Although you requested DRISM run,",
     &        /,4x,"Dipole moment is too small.(",f20.10,")",
     &        /,4x,"Please Check Solvent Orientation.....")
      end
c-------------------------------------------------------------
c     Read namelist drism
c-------------------------------------------------------------
      subroutine drismio(idrism,delec,sparam)
c
      implicit real*8(a-h,o-z)
c
      include "rismio.i"
c
      namelist /drism/idrism,delec,sparam
c
c     idrism      ... Flag to calculate DRISM
c                     0...Not Perform(default)
c                     1...Perform
c     delec       ... Dielecric constant(default=78.5)
c     sparam      ... parameter controlling the length
c                     (default=0.5)
c-------------------------------------------------------------
c
c     Read NAMELIST DRISM
c
      idrism=0
      delec=78.5d0
      sparam=0.5d0
c
      ir=45
      open (ir,file=inpfile)
      read (ir,drism,end=1000)
 1000 continue
      close (ir)
c
      if (idrism.eq.1) then
         write(*,9999) delec
      endif
 9999 format (/,4x,"Dielectric consistent rism option activated."
     &       ,/,4x,"Dielectric constant =",f10.4)
c-------------------------------------------------------------
      return
      end
c
