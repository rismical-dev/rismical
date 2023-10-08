c**************************************************************
c---------------------------------------------------------
c     External program launcher service routine
c---------------------------------------------------------
      subroutine extern3d(iflag)

      implicit real*8 (a-h,o-z)
      
      include "rismio.i"
c      
      logical cuda
      character*256 cupath
c
      namelist /curism/cuda,cupath,ma,param1,param2
c---------------------------------------------------------
      ir=45
      open (ir,file=inpfile)
c
c     Check Requirement of External Program
c
      cuda=.false.
      ma=2
      param1=0.6
      param2=0.2
      rewind ir
      read(ir,curism,end=1000)
 1000 continue
c
      close(ir)
c
c     CUDA 3D-RISM
c
      if (cuda) then
         iflag=1
         call extern3d_cuda(cupath,iflag,ma,param1,param2)
      endif
c
c     Add some other program
c

c
c
c
c
c

c
c
c---------------------------------------------------------
      return
 9999 continue
      iflag=999
      return
      end
c---------------------------------------------------------
c     External CUDA program launcher
c---------------------------------------------------------
      subroutine extern3d_cuda(cupath,iflag,ma,param1,param2)

      implicit real*8 (a-h,o-z)
      character*256 cupath,iolistcu
      character*10 curisminp,curismxmu
      logical fexist

      include "phys_const.i"
      include "rismrun.i"
      include "rismio.i"
      include "solute.i"
      include "solvent.i"
c---------------------------------------------------------
c
c     Check CUDA RISM path
c
      inquire (file=cupath,exist=fexist)
      if (.not.fexist) goto 9000
c
c     IOLIST conversion
c
      iolistcu="m"
      if (index(iolist,'g') + index(iolist,'G')>0) then
         iolistcu=trim(iolistcu)//"g"
      endif
      if (index(iolist,'h') + index(iolist,'H')>0) then
         iolistcu=trim(iolistcu)//"h"
      endif
      if (index(iolist,'c') + index(iolist,'C')>0) then
         iolistcu=trim(iolistcu)//"c"
      endif
      if (index(iolist,'u') + index(iolist,'U')>0) then
         iolistcu=trim(iolistcu)//"u"
      endif
c
c     Generate CUDA RISM Input
c
      curisminp="curism.inp"
      ift=45
      open(ift,file=curisminp,status='replace')
      write(ift,*) trim(iolistcu)//"   0  20220301"
      if (icl.eq.0) then
         write(ift,'(A4)') "HNC"
      elseif (icl.eq.1) then
         write(ift,*) "Unavailable closure type."
         iflag=999
         return
      elseif (icl.eq.2) then
         write(ift,'(A4)') "KH"
      else
         write(ift,*) "Invalid closure type."
         iflag=999
         return
      endif

      write(ift,*) trim(solventxvv)
      write(ift,'(e12.4,i10)') conv,itrmax
      write(ift,'(i3,2f8.3)') ma,param1,param2
      rn=ngrid3d*rdelta3d
      write(ift,'(3f10.4)')  rn,rn,rn
      write(ift,'(3i10)')  ngrid3d,ngrid3d,ngrid3d
      write(ift,'(i10)') nu
      do i=1,nu
         write(ift,'(6f16.8)') siglju(i),epslju(i),qu(i)
     &        ,xyzu(1,i),xyzu(2,i),xyzu(3,i)
      enddo

      close(ift)
c
c     Run CUDA RISM
c
      call system(trim(cupath)//" curism.inp")
c
c     Output CUDA RISM results
c
      call prop3duv_cu

c---------------------------------------------------------
      return
 9000 continue
      write(*,*) "Error. no CUDA RISM binary found at ",cupath
      iflag=999
      stop
      end
c---------------------------------------------------------
c     Physical Property of U-V System for 3D (CUDA)
c---------------------------------------------------------
      subroutine prop3duv_cu
c
      implicit real*8 (a-h,o-z)

      include "phys_const.i"
      include "rismio.i"
      include "rismrun.i"
      include "solvent.i"
      include "solute.i"
c     
      dimension sfec_sc(0:maxslv-1), sfec_gf(0:maxslv-1)
c
      namelist /results/sfe_sc,sfec_sc,sfe_gf,sfec_gf
     &     ,pmv,pressure,correction_term

c---------------------------------------------------------
c
c     Read curism.xmu
c     
      ift=45
      open(ift,file="curism.xmu")
      rewind ift
      read (ift,results)
      close(ift)
c
c     Print Property of U-V System for 3D
c
      write(*,9999) 
      write(*,9998) sfe_sc
      write(*,9993) sfe_gf
c
      write(*,9987) sfe_sc+correction_term,correction_term
c
      write(*,9992)
      do i=0,nvuq-1
         write(*,9991) i+1,sfec_sc(i)
      enddo
      write(*,9990)
      do i=0,nvuq-1
         write(*,9991) i+1,sfec_gf(i)
      enddo
c
      write(*,9989) pmv            ! [L/mol]
     &             ,xt*1.d+9       ! from [/Pa] to [/GPa]
     &             ,pressure       ! [J/m^3]
c---------------------------------------------------------
      return
c---------------------------------------------------------
 9987 format (/,4x,"Solvation Free Energy(PC)     :",g16.8,"[J/mol]"
     &         ,4x,"PC term :",g16.8,"[J/mol]")
 9988 format (/,4x,"Solvation Free Energy(UC)     :",g16.8,"[J/mol]"
     &         ,4x,"UC term :",g16.8,"[J/mol]"
     &         ,4x,"UC param: A=",g16.8
     &            ,"[kJ/mol], B=",g16.8,"[J/mol]"
     &        /,4x,"Solvation Free Energy(GF-UC)  :",g16.8,"[J/mol]"
     &         ,4x,"UC term :",g16.8,"[J/mol]"
     &         ,4x,"UC param: A=",g16.8
     &            ,"[J/mol], B=",g16.8,"[J/mol]")
 9989 format (/,4x,"============= Partial Molar Volume ============="
     &       ,/,4x,"Partial Molar Volume       = ",G12.4," [L/mol]"
     &       ,/,4x,"Isothermal Compressibility = ",G12.4," [/GPa]"
     &       ,/,4x,"Pressure                   = ",G12.4," [GPa]")
 9990 format (/,4x,"======= GF Solvation Free Energy Component =====")
 9991 format (  4x,i4,":",g16.8,"[J/mol]")
 9992 format (/,4x,"======= Solvation Free Energy Component  =======")
 9993 format (/,4x,"Solvation Free Energy(GF)     :",g16.8,"[J/mol]")
 9994 format (/,4x,"----------------------------------------------"
     &       ,/,4x,"NOTE:When HNC+RBC closure are selected,       "  
     &       ,/,4x,"     Solvation Free Energy is not correct.    "  
     &       ,/,4x,"     This value is evaluated based on HNC.    "  
     &       ,/,4x,"----------------------------------------------")
 9995 format (/,4x,"UV Binding Energy Component ES:",g16.8,"[J/mol]")
 9996 format (/,4x,"UV Binding Energy Component LJ:",g16.8,"[J/mol]")
 9997 format (/,4x,"Solute-Solvent Binding Energy :",g16.8,"[J/mol]")
 9998 format (/,4x,"Solvation Free Energy         :",g16.8,"[J/mol]")
 9999 format (/,4x,"======= Physical Property of U-V System =======")
      end
