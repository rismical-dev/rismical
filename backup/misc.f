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
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine upcase  --  convert string to all upper case  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "upcase" converts a text string to all upper case letters
c
c
      subroutine upcasex (string)
      implicit none
      integer i,length,code,ichar
      character*1 char
      character*(*) string
c
c
c     move through the string one character at a time,
c     converting lower case letters to upper case
c
      length = len(string)
      do i = 1, length
         code = ichar(string(i:i))
         if (code.ge.97 .and. code.le.122)
     &      string(i:i) = char(code-32)
      end do
      return
      end
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
      
      write(*,'(A5,i8)') "BUHI:",i
      
      return
      end
c---------------------------------------------------------------
      subroutine dbpf(i,n,f)

      implicit real*8 (a-h,o-z)
      dimension f(n)
      
      write(*,'(A5,i8,5f16.8)') "BUHI:",i,(f(j),j=1,n)
      
      return
      end
c---------------------------------------------------------------
      subroutine dbpi(i,n,integ)

      implicit real*8 (a-h,o-z)
      dimension integ(n)
      
      write(*,'(A5,i8,5i8)') "BUHI:",i,(integ(j),j=1,n)
      
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

