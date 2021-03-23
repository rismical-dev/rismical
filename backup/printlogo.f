c**************************************************************
c--------------------------------------------------------------
c     Printing Software LOGO
c--------------------------------------------------------------
c**************************************************************
      subroutine printlogo(i)
c
c     i=0 : Opening logo
c     i=1 : Finalizing message
c    
c------------------------------------------------------------------
c
      if (i.eq.0) then

      write(*,*) 
     &     "-----------------------------------------------------------"
      write(*,*) "  #####   #   ####   #    #      ####           #    "
      write(*,*) "  #    #  #  #    #  ##  ##  #  #    #          #    " 
      write(*,*) "  #    #  #  #       # ## #     #         ##    #    " 
      write(*,*) "  #####   #   ####   # ## #  #  #           #   #    "
      write(*,*) "  # #     #       #  #    #  #  #        ####   #    "
      write(*,*) "  #  #    #  #    #  #    #  #  #    #  #   #   #    "
      write(*,*) "  #    #  #   ####   #    #  #   ####    ### #  ##   "
      write(*,*) "                                                     "
      write(*,*) "                                                     "
      write(*,*) 
     &     " The Reference Interaction Site-Model integrated Calculator"
      write(*,*) 
     &     " Copyright(C) 2021 -- Norio Yoshida -- All Rights Reserved."
      write(*,*) 
     &     "-----------------------------------------------------------"

      else
         
         write(*,*) 
         write(*,*) 
     &     "-----------------------------------------------------------"
         write(*,*) 
     &     "        RISMiCal computation is completed normally.        "   
         write(*,*) 
     &     "           --- Keep casting, tight lines. ---              "   
         write(*,*) 
     &     "-----------------------------------------------------------"
         write(*,*) 
     &     " Users must cite the following papers for any publications,"   
         write(*,*) 
     &     " which include the results obtained with RISMiCal code.    "   
         write(*,*) 
     &     " References:                                               "   
         write(*,*) 
     &        "N. Yoshida, IOP Conf. Series: Mat. Sci. Eng.773 (2020)",
     &        " 012062 (DOI:10.1088/1757-899X/773/1/012062)  "
         write(*,*) 
     &        "N. Yoshida, J. Chem. Info. Model 57 (2017) 2646-2656 ",
     &        "(DOI:10.1021/acs.jcim.7b00389)  "
         write(*,*) 
     &        "N. Yoshida, J. Chem. Phys. 140 (2014) 214118 ",
     &        "(DOI: 10.1063/1.4879795)  "
         write(*,*) 
     &     "-----------------------------------------------------------"
      endif

c------------------------------------------------------------------
c      return
      end
