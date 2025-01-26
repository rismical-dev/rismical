c
c     Solute parameter
c
      character*4 nsiteu
      character*256 solventxvv

      parameter (maxslu=70000)  ! Solute site
c
c     qu     : solute site charge [e]
c     epslju : solute site LJ parameter [J/mol]
c     siglju : solute site LJ parameter [Ang]
c
c     ipot3d : ESP calculation option for 3D-RISM
c
c
      common /rismslu/qu(maxslu),xyzu(3,maxslu)
     &               ,epslju(maxslu),siglju(maxslu)
     &               ,nu,nsiteu(maxslu),ipot3d

      common /rismslu2/solventxvv



