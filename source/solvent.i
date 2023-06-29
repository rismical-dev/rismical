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
c     qv  : solvent site charge [e]
c     epsljv :solvent site LJ parameter [J/mol]
c     sigljv :solvent site LJ parameter [Ang]
c     temp   : Temperature [K]
c     beta   : 1/kBT [mol/J]
c     xt     : isothermal compressibility [/Pa]
c
      common /rismslv/
     &      nspc(maxslv),temp,beta,xt
     &     ,dens(maxspc)

      common /rismslvmol/nv,nsitev(maxslv),numspc,nvsym(maxslv)
     &     ,qv(maxslv),xyzv(3,maxslv),epsljv(maxslv),sigljv(maxslv)
c
c     reduced solvent parameters 
c
      common /rismslvred/nxvv,nvuq,iuniq(maxslv),q2uq(maxslv)
     &     ,epsljvuq(maxslv),sigljvuq(maxslv),densuq(maxslv)
     &     ,nmulsite(maxslv)
