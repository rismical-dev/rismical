c
c     Solvent parameter
c
      character*4 nsitev

      parameter (maxspc=100)  ! Solvent species
      parameter (maxslv=100)  ! Solvent site
c
c     dens: solvent density [/Ang^3]
c     densuq : density of symmetry unique site [/Ang^3]
c     nspc(i): numbering of solvent species of solvent site i
c     numspc : number of species in solvent
c     nmulspc(i): symetry multipicity of species i
c     qv  : solvent site charge [e]
c     epsljv :solvent site LJ parameter [J/mol]
c     sigljv :solvent site LJ parameter [Ang]
c     temp   : Temperature [K]
c     beta   : 1/kBT [mol/J]
c     xt     : isothermal compressibility [/Pa]
c
      common /rismslv/
     &      nspc(maxslv),nmulsite(maxslv)
     &     ,dens(maxspc),densuq(maxslv)
     &     ,temp,beta,xt

      common /rismslvmol/nv,nsitev(maxslv),numspc
     &     ,qv(maxslv),xyzv(3,maxslv),epsljv(maxslv),sigljv(maxslv)
c
c     solvent parameters for 3D-RISM
c
      common /rismslv_3d/nxvv,nvuq,iuniq(maxslv),q2uq(maxslv)
