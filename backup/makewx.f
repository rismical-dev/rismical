c**************************************************************
c----------------------------------------------------------------
c     Make intramolecular correlation function of solvent
c----------------------------------------------------------------
c**************************************************************
      subroutine makewxv(n,ngrid,rdelta,wk)
c
c     n          ... number of site (=nv)
c     wk         ... k-space intramolecular correlation function
c     rl         ... intramolecular distances [Angstrom]
c     ngrid      ... number of grid
c     rdelta     ... grid width of r-space
c
      implicit real*8 (a-h,o-z)

      include "solvent.i"

      dimension wk(ngrid,n,n),rl(n,n)
c----------------------------------------------------------------
C
C     ---- Calculate Inter-site distance
C
      do i=1,n
         do j=1,n
c            
            if (nspc(i).ne.nspc(j)) then
               rl(i,j)=-1.d0
            else
c
               rl(i,j)=dsqrt((xyzv(1,i)-xyzv(1,j))**2
     &                      +(xyzv(2,i)-xyzv(2,j))**2
     &                      +(xyzv(3,i)-xyzv(3,j))**2)
c
            endif

         enddo
      enddo
C
      call makewx(n,ngrid,rdelta,wk,rl)
c----------------------------------------------------------------
      return
      end
c**************************************************************
c----------------------------------------------------------------
c     Make intramolecular correlation function of solute
c----------------------------------------------------------------
c**************************************************************
      subroutine makewxu(n,ngrid,rdelta,wk)
c
c     n          ... number of site (=nv)
c     wk         ... k-space intramolecular correlation function
c     rl         ... intramolecular distances [Angstrom]
c     ngrid      ... number of grid
c     rdelta     ... grid width of r-space
c
      implicit real*8 (a-h,o-z)

      include "solute.i"

      dimension wk(ngrid,n,n),rl(n,n)
c----------------------------------------------------------------
C
C     ---- Calculate Inter-site distance
C
      do i=1,n
         do j=1,n
c            
            rl(i,j)=dsqrt((xyzu(1,i)-xyzu(1,j))**2
     &                   +(xyzu(2,i)-xyzu(2,j))**2
     &                   +(xyzu(3,i)-xyzu(3,j))**2)

         enddo
      enddo
C
      call makewx(n,ngrid,rdelta,wk,rl)
c----------------------------------------------------------------
      return
      end
c**************************************************************
c----------------------------------------------------------------
c     Make intramolecular correlation function
c----------------------------------------------------------------
c**************************************************************
      subroutine makewx(n,ngrid,rdelta,wk,rl)
c
c     n          ... number of site
c     wk         ... k-space intramolecular correlation function
c     rl         ... intramolecular distances [Angstrom]
c     ngrid      ... number of grid
c     rdelta     ... grid width of r-space
c
      implicit real*8 (a-h,o-z)

      dimension wk(ngrid,n,n),rl(n,n)

      include "phys_const.i"

      deltak=pi/(rdelta*dble(ngrid))
c----------------------------------------------------------------
C
c     ---- Calc intramolecular correlation function
c
      do i=1,n
         do j=1,n
c
c     Identical site
c
            if (i.eq.j) then
               do k=1,ngrid
                  wk(k,i,j)=1.d0
               enddo
            else
c
c     Site on different species
c
               if (rl(i,j).lt.0.d0) then

                  do k=1,ngrid
                     wk(k,i,j)=0.d0
                  enddo

               else
c
c     Site on identical species
c
                  do k=1,ngrid
                     rk=deltak*dble(k)
                     rkl=rk*rl(i,j)
                     wk(k,i,j)=dsin(rkl)/rkl
                  enddo
c
               endif

            endif
         enddo
      enddo
c----------------------------------------------------------------
      return
      end
