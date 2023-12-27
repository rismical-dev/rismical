c--------------------------------------------------------------
c     Converts text string to capital letters
c
      subroutine convcaps(string)
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
         write(*,*) icode
         if (icode.ge.icode_a.and. icode.le.icode_z) then
            string(i:i) = char(icode+icode_aa)
         endif
      end do
c
      return
      end
c
      program main

      implicit real*8(a-h,o-z)
      character*8 string

      string="abz$&Z"

      write(*,*) string      
      call convcaps(string)

      write(*,*) string

      stop
      end
