c-------------------------------------------------------------     
c     namd2rsm
c
c     parameter converter from namd to rism
c
c     usage:
c     $ namd2rsm [type] [pdbfile] [psffile] [prmfile]
c
c     type: 1..RISMiCal 2..CUDA
c
c-------------------------------------------------------------

      program main

      implicit real*8 (a-h,o-z)
      character*2 char2
      character*4 char4,char4a,char4b,char4c
      character*4 atomtype,atomname
      character*5 char5
      character*6 char6,char6a,char6b,char6c
      character*9 char9
      character*80 pdbfil,psffil,prmfil
c
      parameter (mxatom=50000)
c
      dimension xyz(3,mxatom)
      dimension chg(mxatom)
      dimension sigma(mxatom)
      dimension epsi(mxatom)
      dimension atomtype(mxatom)
      dimension atomname(mxatom)

c-------------------------------------------------------------
      nargc = iargc()
      if (nargc.ne.4) then
         write(0,*) "Error. Insuffient arguments. "
         write(0,*) "namd2rsm requires pdbfile, psffile and prmfile."
         write(0,*) "Usage: namd2rsm [type] [pdb] [psf] [prm]"
         stop
      endif
      call  getarg (1,char2)                                       
      call  getarg (2,pdbfil)                                       
      call  getarg (3,psffil)                                       
      call  getarg (4,prmfil)
      read (char2,*) itype
c
c     -------- read pdb ---------
c
      open (3,file=pdbfil,err=9000)

c     skip non "ATOM" line
      iremark=0
 80   continue
      read (3,*) char6
      if (char6.ne."ATOM  ") then
         iremark=iremark+1
         goto 80
      endif
      rewind (3)
      do i=1,iremark
         read (3,*) char6
      enddo
c      
      iatom=0
 100  continue

      read(3,101,end=190) char6,ndum1,char4a,char5,ndum2,x,y,z
 101  format (A6,i5,1x,a4,1x,a5,1x,i3,4x,3f8.3)

      if (char6.eq."ATOM".or.char6.eq."HETATM") then
         iatom=iatom+1
         xyz(1,iatom)=x
         xyz(2,iatom)=y
         xyz(3,iatom)=z
         atomname(iatom)=char4a
      endif
 110  goto 100

 190  continue
      close(3)
c
      natom_pdb=iatom
c
c     -- For debug
c     
c      write(*,*) "NATOM_PDB:",natom_pdb
c      do i=1,natom_pdb
c         write(*,'(i6,2x,3g16.8)') i
c     +        ,xyz(1,i),xyz(2,i),xyz(3,i)
c      enddo
c
c     -------- read psf ---------
c
      open (3,file=psffil,err=9000)

      read (3,*) char4
      read (3,*) ndum
      do i=1,ndum
         read (3,*) char4
      enddo

      read (3,*) natom_psf
      if (natom_psf.ne.natom_pdb) then
         write(0,*) "Inconsistent atom number in psf and pdb."
         write(0,*) "psf:",natom_psf
         write(0,*) "pdb:",natom_pdb
         stop
      endif

      do i=1,natom_psf

         read (3,*) iatom,char2,ndum1,char4a,char4b,char4c,charge
         chg(i)=charge
         atomtype(i)=char4c

      enddo

      close(3)
cc
cc     -- For debug
cc     
c      write(*,*) "NATOM_PDB:",natom_pdb
c      do i=1,natom_pdb
c         write(*,'(i6,2x,a4,2x,4g16.8)') i,atomtype(i)
c     +        ,xyz(1,i),xyz(2,i),xyz(3,i),chg(i)
c      enddo
c
c
c     -------- read prm ---------
c
      open (3,file=prmfil,err=9000)

      do iatom=1,natom_pdb
         rewind (3)

c     --- Find NONBONDED record
 300     continue
         read (3,*) char9
         if (char9.ne."NONBONDED") goto 300

 310     continue
         read(3,*,err=310,end=9010) char4a
         if (char4a.eq.atomtype(iatom)) then
            backspace (3)
            read(3,*,end=9010) char4a, anum1,anum2,anum3
            epsi(iatom)=dabs(anum2)
            sigma(iatom)=anum3*1.781997d0
         else
            goto 310
         endif

      enddo

      close(3)
c
c     -------- write rsm --------
c
      if (itype.eq.1) then ! For NORISM"
         write(*,'(i8,4x,a14)') natom_pdb
         do i=1,natom_pdb
            write(*,'(i4,2x,6f16.8)') 
     +           atomname(i),sigma(i),epsi(i)*4184.d0,chg(i)
     +           ,xyz(1,i),xyz(2,i),xyz(3,i)
         enddo
c
      elseif (itype.eq.2)  then ! For CUDA"
         write(*,'(i8,4x,a14)') natom_pdb
         do i=1,natom_pdb
            write(*,'(6f16.8)') 
     +           sigma(i),epsi(i)*4184.d0,chg(i)
     +           ,xyz(1,i),xyz(2,i),xyz(3,i)
         enddo
c
      endif
c
c     ---- ende ----
c
      stop
c
c     --- error ---
c
 9000 write(0,*) "File does not exist."
      stop
 9010 write(0,*) "End of File during read prm file.",atomtype(iatom)
      stop
      end


