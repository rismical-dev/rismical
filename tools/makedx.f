c
c
c     Make .dx file from UVDATA
c
c
      program main

      implicit real*8 (a-h,o-z)
      character*1 char1
      character*2 char2,char2a
      character*3 char3,charv
      character*4 char4,char4a,char4b
      character*8 char8,fmat
      character*80 char80
      real*8, allocatable :: gr(:,:,:,:)
c
      call getarg(1,char80)
c
      if (char8.eq."-h") then
         write(*,*) "Usage: % makedx.x [guv file]"
         stop
      endif
c
      ift=45
      open (ift, file=char80)
c
c     Read Header
c
      read (ift,*) char2,char2a,char8,char1,nv,ng3d,idum,rdelta3d
      read (ift,*) char2
      read (ift,*) char2,nremark
      do i=1,nremark
         read (ift,*) char2
      enddo
c
      ngrid3d=int(dble(ng3d)**(1.d0/3.d0))
c
      write(*,*) "ngrid3d:",ngrid3d
      write(*,*) "nv:",nv
      write(*,*) "rdelta3d:",rdelta3d
c
c     Read g(r)
c
      allocate (gr(ngrid3d,ngrid3d,ngrid3d,nv))
c
      do j=1,nv
         do kz=1,ngrid3d
         do ky=1,ngrid3d
         do kx=1,ngrid3d
            read(ift,*) g
            gr(kx,ky,kz,j)=g
         enddo
         enddo
         enddo
      enddo
c
      close(ift)
c
c     Output in dx format
c
      do j=1,nv
 
      write(charv,"(I3.3)") j

      jft=46
      open (jft,file=trim(char80)//"."//charv//".dx")
      write(jft,8000)
      write(jft,8001) ngrid3d,ngrid3d,ngrid3d

      xmin=rdelta3d*dble(-ngrid3d/2)

      xmax=rdelta3d*dble(ngrid3d/2-1)
c
      write(jft,8002) xmin,xmin,xmin
      write(jft,8003) rdelta3d
      write(jft,8004) rdelta3d
      write(jft,8005) rdelta3d
      write(jft,8006) ngrid3d,ngrid3d,ngrid3d
      write(jft,8007) ng3d
c
      do ix=1,ngrid3d
      do iy=1,ngrid3d
      do iz=1,ngrid3d

         write(jft,8008) gr(ix,iy,iz,j)
         
      enddo
      enddo
      enddo

      write(jft,8009)
      write(jft,8010)
      write(jft,8011)
      write(jft,8012)
      write(jft,8013)
c
      close (jft)
c
      enddo                     ! of j
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

      

