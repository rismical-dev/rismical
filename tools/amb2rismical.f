c
c     convert parm7 and pdb or crd to 3drism input
c
c     Usage: amb2rismcal [parm7] [pdb|crd] ([format])
c
c     Note: crd file should be ASCII format
c          (corresponding amber option is ntxo=1)
c
c     format=0  Kovarism/Marurism (default)
c     format=1  Norism
c
      program main

      implicit real*8(a-h,o-z)

      character*256 char256
      character*4 char4a,char4b
      character*5 char5
      character*8 char8
      character*78 char78
      character*256 parmfile,pdbfile

      integer, allocatable::nonbondedparamindex(:)
      integer, allocatable::atomtypeindex(:)
      real*8, allocatable::lja(:),ljb(:)

      character*4, allocatable:: atomname(:)
      character*4, allocatable:: atomnamepdb(:)
      real*8, allocatable::sigma(:),epsilon(:),charge(:)
      real*8, allocatable::xyz(:,:)

c-----------------
c
c     Get arguments
c
      nargc=iarg()
      if (nargc.lt.2) then
         write(0,*) "Error. Insufficient arguments."
         stop
      endif

      call getarg(1,parmfile)
      call getarg(2,pdbfile)

      iformat=0
      if (nargc.ge.3) then
         call getarg(3,char5)
         read(char5,*) iformat
         if (iformat.ne.0.and.iformat.ne.1) then
            write(0,*) "Error. Illegal format request."
         endif
      endif

c
c     check coordinate format
      p=scan(pdbfile,".",back=.true.)
      char8=pdbfile(p+1:p+9)
      if (trim(adjustl(char8)).eq."rst7".or.
     &     trim(adjustl(char8)).eq."rst".or.
     &     trim(adjustl(char8)).eq."crd".or.
     &     trim(adjustl(char8)).eq."inpcrd") then
         ifiletype=1
      elseif (trim(adjustl(char8)).eq."pdb".or.
     &     trim(adjustl(char8)).eq."PDB") then
         ifiletype=0
      else
         write(0,*) "Unexpected coordinate file type."
         stop
      endif
c
c     Read parmfile and make LJ and charge parameters
c
      open(3,file=parmfile)
c     read pointers
      rewind (3)
 100  continue
      read(3,'(A256)') char256
      if (char256(1:14).eq."%FLAG POINTERS") then
         read(3,*) char4a
         read(3,*) numatoms,numtypes
      else
         goto 100
      endif
c
      numnonbondedparam=numtypes**2
      numnonbondedparamtype=numtypes*(numtypes+1)/2
      allocate (nonbondedparamindex(numnonbondedparam))
      allocate (lja(numnonbondedparamtype))
      allocate (ljb(numnonbondedparamtype))
      allocate (atomtypeindex(numatoms))
      allocate (atomname(numatoms))
      allocate (sigma(numatoms))
      allocate (epsilon(numatoms))
      allocate (charge(numatoms))
c
c     read atom_name
      rewind(3)
 200  continue
      read(3,'(A256)') char256
      if (char256(1:15).eq."%FLAG ATOM_NAME") then
         read(3,*) char4a
         do i=1,numatoms,20
            jmax=min(numatoms,i+19)
            read(3,'(20a4)') (atomname(j),j=i,jmax)
         enddo
      else
         goto 200
      endif
c     read charge
      rewind(3)
 300  continue
      read(3,'(A256)') char256
      if (char256(1:12).eq."%FLAG CHARGE") then
         read(3,*) char4a
         do i=1,numatoms,5
            jmax=min(numatoms,i+4)
            read(3,'(5e16.8)') (charge(j),j=i,jmax)
         enddo
      else
         goto 300
      endif
c
c     read ATOM_TYPE_INDEX
      rewind(3)
 400  continue
      read(3,'(A256)') char256
      if (char256(1:21).eq."%FLAG ATOM_TYPE_INDEX") then
         read(3,*) char4a
         do i=1,numatoms,10
            jmax=min(numatoms,i+9)
            read(3,'(10i8)') (atomtypeindex(j),j=i,jmax)
         enddo
      else
         goto 400
      endif
c
c     read NONBONDED_PARM_INDEX
      rewind(3)
 500  continue
      read(3,'(A256)') char256
      if (char256(1:26).eq."%FLAG NONBONDED_PARM_INDEX") then
          read(3,*) char4a
         do i=1,numnonbondedparam,10
            jmax=min(numnonbondedparam,i+9)
            read(3,'(10i8)') (nonbondedparamindex(j),j=i,jmax)
         enddo
      else
         goto 500
      endif
c
c     read FLAG LENNARD_JONES_ACOEF                                                       
      rewind(3)
 600  continue
      read(3,'(A256)') char256
      if (char256(1:25).eq."%FLAG LENNARD_JONES_ACOEF") then
         read(3,*) char4a
         do i=1,numnonbondedparamtype,5
            jmax=min(numnonbondedparamtype,i+4)
            read(3,'(5e16.8)') (lja(j),j=i,jmax)
         enddo
      else
         goto 600
      endif
c
c     read FLAG LENNARD_JONES_BCOEF                       
      rewind(3)
 700  continue
      read(3,'(A256)') char256
      if (char256(1:25).eq."%FLAG LENNARD_JONES_BCOEF") then
         read(3,*) char4a
         do i=1,numnonbondedparamtype,5
            jmax=min(numnonbondedparamtype,i+4)
            read(3,'(5e16.8)') (ljb(j),j=i,jmax)
         enddo
      else
         goto 700
      endif
c
      close(3)
c
c     convert to rism format
c
      do iu=1,numatoms
         ij=numtypes*(atomtypeindex(iu)-1)+atomtypeindex(iu)
         id=nonbondedparamindex(ij)
         if (ljb(id).eq.0.d0) then
            sigma(iu)=0.7d0
         else
            sigma(iu)=(2.d0*lja(id)/ljb(id))**(1.d0/6.d0)/2.d0
         endif
         sigma(iu)=sigma(iu)*1.781797d0
         if (ljb(id).eq.0.d0) then
            epsilon(iu)=1.d-2
         else
            epsilon(iu)=ljb(id)**2/(4.d0*lja(id))
         endif

         charge(iu)=charge(iu)/18.2223d0

      enddo
c
c     Read coordinate
c
      allocate (atomnamepdb(numatoms))
      allocate (xyz(3,numatoms))

c
c     case: pdb format
c
      if (ifiletype.eq.0) then
c
         open(3,file=pdbfile)
c     
         maxiu=0
 1000    continue
         read(3,"(A78)",end=2000),char78
         if (char78(1:4).eq."ATOM".or.char78(1:6).eq."HETATM") then
            read(char78(7:11),*) iu
            maxiu=max(iu,maxiu)
            atomnamepdb(iu)=char78(13:16)
            read(char78(31:38),*) xyz(1,iu)
            read(char78(39:46),*) xyz(2,iu)
            read(char78(47:54),*) xyz(3,iu)
            if (trim(adjustl(atomnamepdb(iu)))
     &           .ne.trim(adjustl(atomname(iu)))) then
               write(0,*) "Error. ATOM name inconsistency."
               write(0,*) iu,atomnamepdb(iu),atomname(iu)
               stop
            endif
         endif
         goto 1000
 2000    continue
         close(3)
c
c     case: amber restart format
c
      elseif (ifiletype.eq.1) then
         open(3,file=pdbfile)

         read(3,*) char4a
         read(3,*) maxiu

         do iu=1,maxiu,2
            jmax=min(maxiu,iu+1)
            read(3,'(6f12.7)') ((xyz(i,j),i=1,3),j=iu,jmax)
         enddo

         close(3)
      endif
c
c     chack consistency of number of atoms in parm and coordinate files
c
      if (maxiu.ne.numatoms) then
         write(0,*) "Error. parm7 and pdb inconsistency."
         write(0,*) maxiu,numatoms
         stop
      endif
c
c     write rism inp
c
      write(*,'(i8)') numatoms
      if (iformat.eq.0) then
c     Kovarism/Marurism
         do iu=1,numatoms
            write(*,'(3f10.4,3f12.7)') charge(iu),sigma(iu),epsilon(iu)
     &           ,(xyz(i,iu),i=1,3)
         enddo
      else
c     Norism
         do iu=1,numatoms
            write(*,'(a4,3f10.4,3f12.7)') atomname(iu)
     &           ,sigma(iu),epsilon(iu),charge(iu)
     &           ,(xyz(i,iu),i=1,3)
         enddo
      endif
c
c-------------------------------------------------------------
      stop
      end
      
