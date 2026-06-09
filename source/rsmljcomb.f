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

c-------------------------------------------------------------
      do i=1,n1
         do j=1,n2
            eps=dsqrt(epslj1(i)*epslj2(j))
            sig=(siglj1(i)+siglj2(j))*0.5d0
            epsig(i,j)=eps*sig**norder
         enddo
      enddo
c-------------------------------------------------------------
      return
      end
