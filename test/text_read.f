      program main

      implicit real*8(a-h,o-z)

      character*256 char256a,char256b
      character*2 char2

c--------------------------------
      char256a="## 123  456 789"
      char256b="00 123  456 789"

      if (char256a(1:2).eq."##") then
         write(*,*) "This line is comment line."
         write(*,*) trim(char256a)
      endif

      read(char256b,*) a,b,c,d
      write(*,*) "YOMETA?",a,b,c,d

      stop
      end
      
