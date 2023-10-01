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
     &        "Y. Maruyama, N. Yoshida, H. Tadano, D. Takahashi, ",
     &        "M. Sato, F. Hirata, J. Comput. Chem. 35 (2014) 1347."
         write(*,*) 
     &        "Y. Maruyama, F. Hirata, J. Chem. Theory Comput. ",
     &        "8 (2012) 3015. (DOI: 10.1063/1.4879795)  "
         write(*,*) 
     &     "-----------------------------------------------------------"
         write(*,*) 
     &     " Implimented Math Libraries"
         write(*,*) 
     &     " 3D Fast Fourier Transform: "
         write(*,*) 
     &     " Copyright(C) 2000-2004,2008-2014,2020 Daisuke Takahashi"
     &    ," (e-mail: daisuke[at]cs.tsukuba.ac.jp or ffte[at]ffte.jp)"
         write(*,*) 
     &     " Fast Fourier Transform: "
         write(*,*) 
     &        " Copyright Takuya OOURA, 1996-2001"
         write(*,*)
     &        " Mathmatical Libra: SLATEC Common Mathmatical Library"
         write(*,*)
     &     "Fong, Kirby W.; Jefferson, Thomas H.; Suyehiro, Tokihiko;"
     &    ,"Walton, Lee (July 1993). Guide to the SLATEC"
     &    ," Common Mathematical Library"
         write(*,*) 
     &     "-----------------------------------------------------------"
      endif

c------------------------------------------------------------------
c      return
      end
