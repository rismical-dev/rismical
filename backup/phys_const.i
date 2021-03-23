c     Physical constants
c
      real*8 pi
      parameter (pi=3.1415926535897932d0)
      complex*16 ii
      parameter (ii = (0.0d0, 1.0d0))

      real*8 angtobohr
      parameter (ANGTOBOHR= 1.889726342114D+00) ![bohr/Ang]

      real*8 ekcal,boltz,hart2j,avognum,hart2jmol
      parameter (ekcal=4.184d0)                ![J/cal]
      parameter (hart2j=4.359744722207185d-18) ![J/hartree]
      parameter (avognum=6.02214076d+23)       ![/mol]
      parameter (hart2jmol=hart2j*avognum)     ![J/mol/hartree]
      parameter (boltz=1.380649d-23)           ![J/K]
      parameter (gasconst=boltz*avognum)       ![J/mol/K]

      parameter (fel=1.d0/angtobohr*hart2jmol) ![e**2/Ang -> au -> J/mol]
      
