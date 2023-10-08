c
c
c     Make .txt file to make gOpenmol style .plt file from UVDATA
c
c     Usage: % makeplt.x [ngrid3d] [rdelta3d] [nvuq] [nplt] [format] < [INPUT] > [OUTPUT]
c     Example: % makeplt.x 128 0.5 2 1 SIMPLE < UVDATA > gr_O.txt
c
c
      program main

      implicit real*8 (a-h,o-z)
      character*1 char1
      character*8 char8,fmat
      real*8, allocatable :: gr(:,:,:)
      real*8, allocatable :: tr(:,:,:)
c
cc      read (*,*) ngrid3d,rdelta3d,nvuq,nplt
      call getarg(1,char8)
      read (char8,*) ngrid3d
      call getarg(2,char8)
      read (char8,*) rdelta3d
      call getarg(3,char8)
      read (char8,*) nvuq
      call getarg(4,char8)
      read (char8,*) nplt
      call getarg(5,char8)
      read (char8,*) fmat
c
c
      allocate (gr(ngrid3d,ngrid3d,ngrid3d))
      allocate (tr(ngrid3d,ngrid3d,ngrid3d))
c
      if (fmat.eq."DETAIL") then
      do j=1,nplt
         
         read(*,*) char1
         read(*,*) char1
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d
            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
            read(*,*) rx,ry,rz,dum1,dum2,g,dum3,dum4
            if (j.eq.nplt) gr(kx,ky,kz)=g
         enddo
         read(*,*)
         enddo
         read(*,*)
         enddo
         
      enddo
c
      elseif (fmat.eq."SIMPLE") then
      read(*,*) char1
      do j=1,nplt
         
         read(*,*) char1
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d
            k=kx+(ky-1)*ngrid3d+(kz-1)*ngrid3d**2
            read(*,*) dum1,g
            if (j.eq.nplt) gr(kx,ky,kz)=g
         enddo
         read(*,*)
         enddo
         read(*,*)
         enddo
         
      enddo
c
      elseif (fmat.eq."UNFORMAT") then
         write(*,*) "UNFORMAT is not available. Please use SIMPLE."
         stop
      endif
c
      write(*,*) "3 200"
      write(*,*) ngrid3d,ngrid3d,ngrid3d

      xmin=rdelta3d*dble(-ngrid3d/2)

      xmax=rdelta3d*dble(ngrid3d/2-1)
c
      write(*,9999) xmin,xmax,xmin,xmax,xmin,xmax
c
      do iz=1,ngrid3d
      do iy=1,ngrid3d
      do ix=1,ngrid3d

         write(*,9998) gr(ix,iy,iz)
         
      enddo
      enddo
      enddo
c
      deallocate (gr)
c
      stop
 9998 format (e16.8)
 9999 format (6f9.3)
      end

      

