c
c
c     Make .dx file from UVDATA
c
c     Usage: % makedx.x [ngrid3d] [rdelta3d] [nvuq] [nplt] [format] < [INPUT] > [OUTPUT]
c     Example: % makedx.x 128 0.5 2 1 < UVDATA > gr_O.txt
c
c
      program main

      implicit real*8 (a-h,o-z)
      character*1 char1
      character*8 char8,fmat
      real*8, allocatable :: gr(:,:,:)
c
cc      read (*,*) ngrid3d,rdelta3d,nvuq,nplt
      call getarg(1,char8)
c
      if (char8.eq."-h") then
         write(*,*) "Usage: % makedx.x ",
     . "[ngrid3d] [rdelta3d] [nvuq] [nplt] [format]",
     . " < [INPUT] > [OUTPUT]"
         write(*,*) "Example: % makedx.x 128 0.5 2 1",
     . "< UVDATA > gr_O.txt"
         stop
      endif
c
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
      allocate (gr(ngrid3d,ngrid3d,ngrid3d))
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
      write(*,8000)
      write(*,8001) ngrid3d,ngrid3d,ngrid3d

      xmin=rdelta3d*dble(-ngrid3d/2)

      xmax=rdelta3d*dble(ngrid3d/2-1)
c
      write(*,8002) xmin,xmin,xmin
      write(*,8003) rdelta3d
      write(*,8004) rdelta3d
      write(*,8005) rdelta3d
      write(*,8006) ngrid3d,ngrid3d,ngrid3d
      write(*,8007) ngrid3d**3
c
      do ix=1,ngrid3d
      do iy=1,ngrid3d
      do iz=1,ngrid3d

         write(*,8008) gr(ix,iy,iz)
         
      enddo
      enddo
      enddo

      write(*,8009)
      write(*,8010)
      write(*,8011)
      write(*,8012)
      write(*,8013)
c
      deallocate (gr)
c
      stop
 9998 format (e16.8)
 9999 format (6f9.3)

 8000 format ("# Comments")
 8001 format ("object 1 class gridpositions counts",3i5)
 8002 format ("origin",3f16.8)
 8003 format ("delta ",f5.2," 0.0 0.0")
 8004 format ("delta 0.0 ",f5.2," 0.0 ")
 8005 format ("delta 0.0 0.0 ",f5.2)
 8006 format ("object 2 class gridconnections counts",3i5)
 8007 format ("object 3 class array type double rank 0 items "
     +        ,i10," data follows")
 8008 format (f16.8)
 8009 format ('attribute "dep" string "positions"')
 8010 format 
     + ('object "regular positions regular connections" class field')
 8011 format ('component "positions" value 1')
 8012 format ('component "connections" value 2')
 8013 format ('component "data" value 3')

      end

      

