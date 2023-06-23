c------------------------------------------------------------------
c                                                                   
c        #####   #   ####   #    #      ####           #    
c        #    #  #  #    #  ##  ##  #  #    #          #   
c        #    #  #  #       # ## #     #         ##    #   
c        #####   #   ####   # ## #  #  #           #   #   
c        # #     #       #  #    #  #  #        ####   #   
c        #  #    #  #    #  #    #  #  #    #  #   #   #   
c        #    #  #   ####   #    #  #   ####    ### #  ##  
c
c     The Reference Interaction Site-Model integrated Calculator
c     Copyright(C) 2021 -- Norio Yoshida -- All Rights Reserved.
c     Copyright(C) 2023 -- Yutaka Maruyama, Norio Yoshida 
c     All Rights Reserved.
c
c------------------------------------------------------------------
c
      program rismical
c
      implicit real*8(a-h,o-z)
      character*8 carg8
c
      call printlogo(0)

c
      call getarg(1,carg8)
      call upcasex(carg8)
c
      call rismiofile
c
      if (trim(carg8).eq."VV") then

         write(*,9999) " 1D-RISM run for solvent-solvent system"
         call rismicalvv

      elseif (trim(carg8).eq."1D") then

         write(*,9999) " 1D-RISM run for solute-solvent system"
         call rismical1d

      elseif (trim(carg8).eq."3D") then

         write(*,9999) " 3D-RISM run for solute-solvent system"
         call rismical3d
         
c$$$      elseif (trim(carg8).eq."MOZVV") then
c$$$         call rismicalmozvv
c$$$      elseif (trim(carg8).eq."MOZUV") then
c$$$         call rismicalmozuv
c$$$      elseif (trim(carg8).eq."XMOZ") then
c$$$         call rismicalxmoz
         
      else

         write(*,*) "Invalid option."
         write(*,*) "1st argument option: VV  1D  3D"
         ierr=999
         call abrt(ierr)

      endif
c
      call printlogo(1)
c
      stop
 9999 format (/,2x,"-------------------------------------------------"
     &       ,/,2x,A40
     &       ,/,2x,"-------------------------------------------------")
      end


      
