c**************************************************************
c--------------------------------------------------------------
c     Abrt
c--------------------------------------------------------------
c**************************************************************
      subroutine abrt(i)

      implicit real*8(a-h,o-z)
c--------------------------------------------------------------
C
C     ----- PRINT ACCOUNTING INFO AND GOODBYE MESSAGE -----
C
      write (*,900) i
  900 FORMAT(1X,'EXECUTION OF RISMiCal TERMINATED -ABNORMALLY-'
     &     /,1X,'ERROR CODE:',i10)
      stop "IN ABRT"
c--------------------------------------------------------------
      end
c--------------------------------------------------------------
c     Converts text string to capital letters
c
      subroutine upcasex(string)
      implicit none
      character*1 char
      character*(*) string
      integer leng_str,icode_a,icode_z,icode_aa,icode,i
c
c
      icode_a=ichar("a") 
      icode_z=ichar("z") 
      icode_aa=ichar("A")-icode_a

      leng_str = len(string)
      do i = 1, leng_str
         icode = ichar(string(i:i))
         if (icode.ge.icode_a.and. icode.le.icode_z) then
            string(i:i) = char(icode+icode_aa)
         endif
      end do
      return
c--------------------------------------------------------------
      end
c
c-------------------------------------------------------------
c     Combine the LJ parameter
c-------------------------------------------------------------
      subroutine rsmljcomb(norder,n1,n2,
     &     siglj1,epslj1,siglj2,epslj2,epsig)
      
      implicit real*8 (a-h,o-z)
      character*17 char17
      character*8 char8

      dimension siglj1(n1),siglj2(n2)
      dimension epslj1(n1),epslj2(n2)
      dimension epsig(n1,n2)
c$$$      dimension epsig12(n1,n2)
c$$$      dimension epsig6(n1,n2)

c-------------------------------------------------------------
      do i=1,n1
         do j=1,n2
            eps=dsqrt(epslj1(i)*epslj2(j))
            sig=(siglj1(i)+siglj2(j))*0.5d0
c$$$            epsig12(i,j)=eps*sig**12
c$$$            epsig6(i,j)=eps*sig**6
            epsig(i,j)=eps*sig**norder
         enddo
      enddo
c$$$
c$$$      if (norder.eq.6) then
c$$$         do i=1,n1
c$$$            do j=1,n2
c$$$               epsig(i,j)=epsig6(i,j)
c$$$            enddo
c$$$         enddo
c$$$      elseif (norder.eq.12) then
c$$$         do i=1,n1
c$$$            do j=1,n2
c$$$               epsig(i,j)=epsig12(i,j)
c$$$            enddo
c$$$         enddo
c$$$      else
c$$$         write(*,*) "Wrong order of LJ combination."
c$$$         call abrt
c$$$      endif
c-------------------------------------------------------------
      return
      end
c**************************************************************
c--------------------------------------------------------------
c     debug print
c--------------------------------------------------------------
c**************************************************************
      subroutine dbp(i)

      implicit real*8 (a-h,o-z)
      
      write(*,'(A5,i8)') "For Debug:",i
      
      return
      end
c---------------------------------------------------------------
      subroutine dbpf(i,n,f)

      implicit real*8 (a-h,o-z)
      dimension f(n)
      
      write(*,'(A5,i8,5f16.8)') "For Debug:",i,(f(j),j=1,n)
      
      return
      end
c---------------------------------------------------------------
      subroutine dbpi(i,n,integ)

      implicit real*8 (a-h,o-z)
      dimension integ(n)
      
      write(*,'(A5,i8,5i8)') "For Debug:",i,(integ(j),j=1,n)
      
      return
      end
c---------------------------------------------------------------
      subroutine vclr(a,inca,n)
c
      implicit double precision(a-h,o-z)
c
      dimension a(*)
c
      parameter (zero=0.0d+00)
c
c     ----- zero out vector -a-, using increment -inca- -----
c
      if (inca .ne. 1) go to 200
      do 110 l=1,n
         a(l) = zero
  110 continue
      return
c
  200 continue
      la=1-inca
      do 210 l=1,n
         la=la+inca
         a(la) = zero
  210 continue
      return
      end
c---------------------------------------------------------------
C*MODULE MTHLIB  *DECK VCLR for OpenMP
c---------------------------------------------------------------
      subroutine vclr_mp(a,inca,n)
c
      implicit double precision(a-h,o-z)
c
      dimension a(*)
c
      parameter (zero=0.0d+00)
c
c     ----- zero out vector -a-, using increment -inca- -----
c
      if (inca .ne. 1) go to 200
!$omp parallel do
      do l=1,n
         a(l) = zero
      enddo
!$omp end parallel do
      return
c
  200 continue
      la=1-inca
      do l=1,n
         la=la+inca
         a(la) = zero
      enddo
      continue
      return
      end
c---------------------------------------------------------------
C*MODULE MTHLIB  *DECK VCLR for OpenMP
c---------------------------------------------------------------
      subroutine vclrz_mp(a,inca,n)
c
      implicit double precision(a-h,o-z)
      complex*16 a
c
      dimension a(*)
c
      parameter (zero=0.0d+00)
c
c     ----- zero out vector -a-, using increment -inca- -----
c
      if (inca .ne. 1) go to 200
!$omp parallel do
      do l=1,n
         a(l) = zero
      enddo
!$omp end parallel do
      return
c
  200 continue
      la=1-inca
      do l=1,n
         la=la+inca
         a(la) = zero
      enddo
      continue
      return
      end

