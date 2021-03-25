************************************************************************
*           CALCULATION OF ANGULAR AVERAGES OF 3D PROFILES             *
*           ----------------------------------------------             *
*                                                                      *
*  Input:   3D site profiles   (i.e. UVDATA)                           *
*  Input:   Status sheet       (i.e. STATINP)                          *
*  Output:  radial distribution profiles around given centers          *
*                                                                      *
*  Usage: # g3d1d.x STATINP UVDATA                                     *
*                                                                      *
*                                                                      *
*  ------------------------------------------------------------------  *
*                                                                      *
*                                                                      *
*  Parameters in the Input file Status sheet:                          *
*                                                                      *
*   NR0  -  number of centers to average around and of output files    *
*   OUTEXT, X0,Y0,Z0  -  1-st output file suffix and                   *
*     ...                averaging center coordinates                  *
*   OUTEXT, X0,Y0,Z0  -  Nr0-th "-"                                    *
*   ngrid3d rdelta3d natv  - # of 3D-RISM grid, width, # of solvent    *
*   RMIN,RMAX,DELR  -  lower and upper limits, and step of the output  *
*                      radial grid                                     *
*                                                                      *
*                                                                      *
*  ------------------------------------------------------------------  *
*                Written by Andriy Kovalenko, April 1999               *
*                Modified by Norio Yoshida, Novenber 2006              *
*                Modified by Norio Yoshida, Febrary 2010               *
************************************************************************
                                                                        
      Program  ANGULAR_AVERAGE                                          

      implicit real*8(a-h,o-z)
      real*8, allocatable :: guv(:,:)
c      parameter (maxgx=256,maxgy=256,maxgz=256, maxatv=10,               
      parameter (maxgx=512,maxgy=512,maxgz=512, maxatv=10,               
c     parameter (maxgx=64,maxgy=64,maxgz=64, maxatv=3,                  
     .           maxr0=1000, maxnod=1000)                                
      parameter (maxgr=maxgx*maxgy*maxgz)                               
      character  filnam*40, quadra*40                                    
cIMAI      character*12  inpfil, angnod, guvext, guvdat,                     
      character*80  inpfil, angnod, guvext, guvdat,                     
     .              outext(maxr0), outfil(maxr0)                        
      character*2 char2
      character*10 char10
      character*16 char16
      integer igc(3,0:1),ndum
      real  box(3), origin(3), rv(3), delc(3,0:1),                      
     .      ravce0(3,maxr0), ravcen(3,maxr0),                           
     .      phi(maxnod), theta(maxnod), angwei(maxnod),                 
     .      cosphi(maxnod), sinphi(maxnod),                             
     .      costhe(maxnod), sinthe(maxnod),                             
     .      gvv(maxatv), gvv8(maxatv)                                  
cIMAI
      REAL XXX
c
                                                                        
c........................ reading input data ...........................
      call  getarg (1,inpfil)                                       
      call  getarg (2,filnam)                                       
                                                                        
      open (3,file=inpfil,status='old',err=90)                          
        read (3,*,err=90)  nr0                                          
        if (nr0.gt.maxr0)  then                                         
          print *,'O3DAV:  NR0 > MaxR0'                                 
          goto 99                                                       
        endif                                                           
        do ir0=1,nr0                                                    
          read (3,*,err=90)  outext(ir0), (ravce0(id,ir0),id=1,3)       
        enddo                                                           
        read (3,*,err=90)  rmin, rmax, delr                             
      close (3,err=90)                                                  
c
c
c
c........................ reading profile grid .........................
      print *,'reading Guv data:  ',guvdat,'...'                        

      open (10,file=guvdat,status='old',err=90)                         
      read (10,'(A16,3i8,f16.8)',err=90) 
     &     char16,natv,ngrid3d,idum,rdelta3d

      read (10,*) char2
      read (10,*) char2,iskip
      do i=1,iskip
         read(10,*) char2
      enddo

      ng33=ngrid3d**3
      allocate (guv(ng33,natv))

      do iv=1,natv
         do ig=1,ng33
            read (10,*,err=90)  dum1,guv(ig,iv)
         enddo                                                           
      enddo                                                           
      close (10,err=90)                                                 
                                                                        
c.......................................................................
c
      do id=1,3
         origin(id)=rdelta3d*ngrid3d/2
         box(id)=dble(ngrid3d)*rdelta3d
      enddo
c.......................... setting I/O files ..........................
      do il=1,len(filnam)                                               
        if (filnam(il:il).eq.' ')  then                                 
          lfil = il-1                                                   
          goto 10                                                       
        endif                                                           
      enddo                                                             
      lfil = len(filnam)                                                
   10 if (lfil.eq.0)  then                                              
        print *,'O3DAV:  input file not specified'                      
        goto 99                                                         
      endif                                                             
                                                                        
      guvdat = filnam(1:lfil)   !! // guvext
      do ir0=1,nr0                                                      
        outfil(ir0) = filnam(1:lfil) // outext(ir0)                     
      enddo                                                             
                                                                        
c..................... reading angular nodes data ......................
      call defaultread(nod,phi,theta,angwei)
c............... tabulating Sin and Cos of nodes angles ................
      do inod=1,nod                                                     
        cosphi(inod) = cos( phi(inod))                                 
        sinphi(inod) = sin( phi(inod))                                 
        costhe(inod) = cos( theta(inod))                               
        sinthe(inod) = sin( theta(inod))                               
      enddo                                                             
                                                                        
c.................... shifting origin to box corner ....................
      do ir0=1,nr0                                                      
      do id=1,3                                                         
        ravcen(id,ir0) = ravce0(id,ir0) + origin(id)                    
      enddo                                                             
      enddo                                                             
                                                                        
c....................enumerating averaging centers ....................
      print *,'output file       averaging center (x,y,z)'              
      do ir0=1,nr0                                                      
        print *, outfil(ir0),'   ',(ravce0(id,ir0),id=1,3)              
                                                                        
        open (11,file=outfil(ir0),err=90)                               
                                                                        
c................... enumerating radial grid points ....................
          do r=rmin,rmax+0.5*delr,delr                                  
c..................... summing over angular nodes ......................
            do iv=1,natv                                                
              gvv(iv) = 0.                                              
            enddo                                                       
                                                                        
            do inod=1,nod                                               
                                                                        
c................... getting vector at angular node ....................
              rv(1) = ravcen(1,ir0) + r*sinthe(inod)*cosphi(inod)       
              rv(2) = ravcen(2,ir0) + r*sinthe(inod)*sinphi(inod)       
              rv(3) = ravcen(3,ir0) + r*costhe(inod)                    
                                                                        
c................ applying periodic boundary conditions ................
              do id=1,3                                                 
                rv(id) = rv(id) - box(id)*anint( rv(id)/box(id)-0.5)    
              enddo                                                     
                                                                        
c......... specifying corners and sides of interpolation cell ..........
              do id=1,3                                                 
                rv0 = ngrid3d*rv(id)/box(id)                            
                igc(id,0) = mod( int(rv0), ngrid3d)                    
                igc(id,1) = mod( igc(id,0)+1, ngrid3d)                   
                delc(id,1) = rv0 - igc(id,0)                            
                delc(id,0) = 1. - delc(id,1)                            
              enddo                                                     
                                                                        
c............. summing Guv over interpolation cell corners .............
              do iv=1,natv                                              
                gvv8(iv) = 0.                                           
              enddo                                                     
                                                                        
              do id1=0,1                                                
              do id2=0,1                                                
              do id3=0,1                                                
                                                                        
c................. getting cell corner and its weight ..................
                ig = 1 + igc(1,id1) + igc(2,id2)*ngrid3d                 
     .                              + igc(3,id3)*ngrid3d**2          
                weight = delc(1,id1)*delc(2,id2)*delc(3,id3)            
                                                                        
c.................... maintaining interpolated Guv .....................
                do iv=1,natv                                            
                  gvv8(iv) = gvv8(iv) + weight*guv(ig,iv)               
                enddo                                                   
                                                                        
              enddo                                                     
              enddo                                                     
              enddo                                                     
                                                                        
c...................... maintaining averaged Gvv .......................
              do iv=1,natv                                              
                gvv(iv) = gvv(iv) + angwei(inod)*gvv8(iv)               
              enddo                                                     
                                                                        
            enddo                                                       
                                                                        
c.......................... outputting Gvv(r) ..........................
            write (11,111)  r, (gvv(iv),iv=1,natv)                      
  111       format( f8.4, 20(1x,f11.7) )                                
                                                                        
          enddo                                                         
                                                                        
        close (11,err=90)                                               
                                                                        
      enddo                                                             
c
      deallocate (guv)
c                                                                        
   99 stop                                                              
   90 print *,'O3DAV:  I/O error'                                       
      goto 99                                                           
                                 e n d                                  
      subroutine defaultread(nod,phi,theta,angwei)
      parameter (maxnod=1000)                                
      real   phi(maxnod), theta(maxnod), angwei(maxnod)
c
      nod=700

       phi(  1)   =  0.44280643E+01
       theta(  1) =  0.65374666E+00
       angwei(  1)=  0.14271057E-02

       phi(  2)   =  0.53645711E+01
       theta(  2) =  0.25645879E+00
       angwei(  2)=  0.14485514E-02

       phi(  3)   =  0.37072277E+01
       theta(  3) =  0.13156301E+01
       angwei(  3)=  0.14340515E-02

       phi(  4)   =  0.15577928E+01
       theta(  4) =  0.18487531E+01
       angwei(  4)=  0.14202882E-02

       phi(  5)   =  0.18717331E+01
       theta(  5) =  0.92745876E+00
       angwei(  5)=  0.14396844E-02

       phi(  6)   =  0.26844831E+01
       theta(  6) =  0.11569905E+01
       angwei(  6)=  0.14570996E-02

       phi(  7)   =  0.60148783E+01
       theta(  7) =  0.24038584E+01
       angwei(  7)=  0.14231028E-02

       phi(  8)   =  0.31313857E+00
       theta(  8) =  0.11268555E+01
       angwei(  8)=  0.14246602E-02

       phi(  9)   =  0.35869280E+00
       theta(  9) =  0.17075888E+01
       angwei(  9)=  0.14401994E-02

       phi( 10)   =  0.59223957E+01
       theta( 10) =  0.11804222E+01
       angwei( 10)=  0.14287453E-02

       phi( 11)   =  0.57193203E+01
       theta( 11) =  0.17844205E+01
       angwei( 11)=  0.14325900E-02

       phi( 12)   =  0.17817345E+01
       theta( 12) =  0.80102587E+00
       angwei( 12)=  0.14327910E-02

       phi( 13)   =  0.41235247E+01
       theta( 13) =  0.24884918E+01
       angwei( 13)=  0.14255943E-02

       phi( 14)   =  0.57310166E+01
       theta( 14) =  0.15027913E+01
       angwei( 14)=  0.14386674E-02

       phi( 15)   =  0.32443554E+01
       theta( 15) =  0.10519085E+01
       angwei( 15)=  0.14307159E-02

       phi( 16)   =  0.61475911E+01
       theta( 16) =  0.25207539E+01
       angwei( 16)=  0.14469452E-02

       phi( 17)   =  0.17387435E+00
       theta( 17) =  0.10599604E+01
       angwei( 17)=  0.14299089E-02

       phi( 18)   =  0.28774965E+01
       theta( 18) =  0.22469168E+01
       angwei( 18)=  0.14323026E-02

       phi( 19)   =  0.29376523E+01
       theta( 19) =  0.14705707E+01
       angwei( 19)=  0.14166960E-02

       phi( 20)   =  0.12336500E+01
       theta( 20) =  0.22694201E+01
       angwei( 20)=  0.15035272E-02

       phi( 21)   =  0.92751664E+00
       theta( 21) =  0.24198983E+01
       angwei( 21)=  0.14053179E-02

       phi( 22)   =  0.45858727E+01
       theta( 22) =  0.26266153E+01
       angwei( 22)=  0.15375474E-02

       phi( 23)   =  0.35033250E+01
       theta( 23) =  0.16438702E+01
       angwei( 23)=  0.14197548E-02

       phi( 24)   =  0.60484159E+00
       theta( 24) =  0.15543230E+01
       angwei( 24)=  0.14410883E-02

       phi( 25)   =  0.20665193E+00
       theta( 25) =  0.13541772E+01
       angwei( 25)=  0.14370964E-02

       phi( 26)   =  0.23809373E+01
       theta( 26) =  0.10243875E+01
       angwei( 26)=  0.14422248E-02

       phi( 27)   =  0.45255029E+00
       theta( 27) =  0.89273375E+00
       angwei( 27)=  0.14107986E-02

       phi( 28)   =  0.20663807E+01
       theta( 28) =  0.13323878E+01
       angwei( 28)=  0.14292873E-02

       phi( 29)   =  0.35634283E-01
       theta( 29) =  0.26334403E+01
       angwei( 29)=  0.14275265E-02

       phi( 30)   =  0.59795833E+01
       theta( 30) =  0.70057172E+00
       angwei( 30)=  0.14543344E-02

       phi( 31)   =  0.91808575E+00
       theta( 31) =  0.44961619E+00
       angwei( 31)=  0.14598066E-02

       phi( 32)   =  0.21721392E+01
       theta( 32) =  0.83381391E+00
       angwei( 32)=  0.14311024E-02

       phi( 33)   =  0.25419325E-01
       theta( 33) =  0.22624841E+01
       angwei( 33)=  0.14274773E-02

       phi( 34)   =  0.52620387E+01
       theta( 34) =  0.12487100E+01
       angwei( 34)=  0.14391389E-02

       phi( 35)   =  0.62644110E+01
       theta( 35) =  0.85835236E+00
       angwei( 35)=  0.14411777E-02

       phi( 36)   =  0.39943495E+01
       theta( 36) =  0.18410881E+01
       angwei( 36)=  0.14505351E-02

       phi( 37)   =  0.70746481E+00
       theta( 37) =  0.29024251E+01
       angwei( 37)=  0.14134037E-02

       phi( 38)   =  0.36517384E+01
       theta( 38) =  0.89180726E+00
       angwei( 38)=  0.14256298E-02

       phi( 39)   =  0.45573478E+01
       theta( 39) =  0.20075052E+01
       angwei( 39)=  0.14372853E-02

       phi( 40)   =  0.36063747E+01
       theta( 40) =  0.15427278E+01
       angwei( 40)=  0.14321239E-02

       phi( 41)   =  0.37505736E+01
       theta( 41) =  0.25469460E+01
       angwei( 41)=  0.12576147E-02

       phi( 42)   =  0.16028181E+01
       theta( 42) =  0.21387801E+01
       angwei( 42)=  0.14567050E-02

       phi( 43)   =  0.73707062E+00
       theta( 43) =  0.16183933E+01
       angwei( 43)=  0.14173348E-02

       phi( 44)   =  0.49461145E+01
       theta( 44) =  0.99037468E+00
       angwei( 44)=  0.14332122E-02

       phi( 45)   =  0.53017282E+01
       theta( 45) =  0.21647341E+01
       angwei( 45)=  0.14239332E-02

       phi( 46)   =  0.61818795E+01
       theta( 46) =  0.27594748E+01
       angwei( 46)=  0.14242479E-02

       phi( 47)   =  0.12513140E+01
       theta( 47) =  0.11167758E+01
       angwei( 47)=  0.15147439E-02

       phi( 48)   =  0.31657171E+01
       theta( 48) =  0.18542650E+00
       angwei( 48)=  0.14422931E-02

       phi( 49)   =  0.59641933E+01
       theta( 49) =  0.22469656E+01
       angwei( 49)=  0.15697015E-02

       phi( 50)   =  0.25092125E+01
       theta( 50) =  0.64665574E+00
       angwei( 50)=  0.14327602E-02

       phi( 51)   =  0.45696826E+01
       theta( 51) =  0.13806660E+01
       angwei( 51)=  0.14370638E-02

       phi( 52)   =  0.26264286E+01
       theta( 52) =  0.19555103E+01
       angwei( 52)=  0.14236935E-02

       phi( 53)   =  0.52031450E+01
       theta( 53) =  0.76522291E+00
       angwei( 53)=  0.14362399E-02

       phi( 54)   =  0.51018734E+01
       theta( 54) =  0.12343985E+01
       angwei( 54)=  0.14296726E-02

       phi( 55)   =  0.20554714E+01
       theta( 55) =  0.17238516E+01
       angwei( 55)=  0.14296651E-02

       phi( 56)   =  0.62033019E+01
       theta( 56) =  0.71266836E+00
       angwei( 56)=  0.14169238E-02

       phi( 57)   =  0.18311423E+00
       theta( 57) =  0.28738670E+01
       angwei( 57)=  0.14776009E-02

       phi( 58)   =  0.21597095E+01
       theta( 58) =  0.19957911E+01
       angwei( 58)=  0.14431784E-02

       phi( 59)   =  0.27834427E+01
       theta( 59) =  0.29573071E+01
       angwei( 59)=  0.14428508E-02

       phi( 60)   =  0.31943719E+01
       theta( 60) =  0.14234252E+01
       angwei( 60)=  0.14146473E-02

       phi( 61)   =  0.28931243E+01
       theta( 61) =  0.21069801E+01
       angwei( 61)=  0.14310009E-02

       phi( 62)   =  0.41406693E+01
       theta( 62) =  0.28804684E+01
       angwei( 62)=  0.14025802E-02

       phi( 63)   =  0.26726303E+01
       theta( 63) =  0.86437732E+00
       angwei( 63)=  0.14410461E-02

       phi( 64)   =  0.33640995E+01
       theta( 64) =  0.33215660E+00
       angwei( 64)=  0.14197140E-02

       phi( 65)   =  0.39609215E+01
       theta( 65) =  0.11763502E+01
       angwei( 65)=  0.14360076E-02

       phi( 66)   =  0.18243589E+01
       theta( 66) =  0.28322580E+01
       angwei( 66)=  0.14299372E-02

       phi( 67)   =  0.14907879E+00
       theta( 67) =  0.90938318E+00
       angwei( 67)=  0.14193122E-02

       phi( 68)   =  0.11904752E+01
       theta( 68) =  0.20658755E+01
       angwei( 68)=  0.13721366E-02

       phi( 69)   =  0.44778166E+01
       theta( 69) =  0.21259813E+01
       angwei( 69)=  0.14304656E-02

       phi( 70)   =  0.25418339E+01
       theta( 70) =  0.22230828E+01
       angwei( 70)=  0.14231077E-02

       phi( 71)   =  0.62946558E+00
       theta( 71) =  0.21207232E+01
       angwei( 71)=  0.14256808E-02

       phi( 72)   =  0.45754504E+01
       theta( 72) =  0.22476466E+01
       angwei( 72)=  0.14352525E-02

       phi( 73)   =  0.36260948E+01
       theta( 73) =  0.20453844E+01
       angwei( 73)=  0.14440342E-02

       phi( 74)   =  0.11109715E+01
       theta( 74) =  0.13819329E+01
       angwei( 74)=  0.14432662E-02

       phi( 75)   =  0.30008790E+01
       theta( 75) =  0.12397895E+01
       angwei( 75)=  0.14590552E-02

       phi( 76)   =  0.71676475E+00
       theta( 76) =  0.33072332E+00
       angwei( 76)=  0.14505428E-02

       phi( 77)   =  0.48635578E+00
       theta( 77) =  0.20604506E+01
       angwei( 77)=  0.14237369E-02

       phi( 78)   =  0.32879703E+01
       theta( 78) =  0.18419805E+01
       angwei( 78)=  0.14271941E-02

       phi( 79)   =  0.38305502E+01
       theta( 79) =  0.12549341E+01
       angwei( 79)=  0.14297583E-02

       phi( 80)   =  0.23220222E+01
       theta( 80) =  0.73297036E+00
       angwei( 80)=  0.14313650E-02

       phi( 81)   =  0.47784643E+01
       theta( 81) =  0.18757085E+01
       angwei( 81)=  0.14307266E-02

       phi( 82)   =  0.51645741E+01
       theta( 82) =  0.15929849E+01
       angwei( 82)=  0.13177810E-02

       phi( 83)   =  0.35661860E+01
       theta( 83) =  0.12638683E+01
       angwei( 83)=  0.14349418E-02

       phi( 84)   =  0.48235053E+00
       theta( 84) =  0.16307325E+01
       angwei( 84)=  0.14267901E-02

       phi( 85)   =  0.28908281E+01
       theta( 85) =  0.18182906E+01
       angwei( 85)=  0.14317775E-02

       phi( 86)   =  0.50831389E+00
       theta( 86) =  0.60975730E+00
       angwei( 86)=  0.14426104E-02

       phi( 87)   =  0.31579382E+01
       theta( 87) =  0.12908852E+01
       angwei( 87)=  0.14066888E-02

       phi( 88)   =  0.37251940E+01
       theta( 88) =  0.26778104E+01
       angwei( 88)=  0.14638340E-02

       phi( 89)   =  0.11027013E+01
       theta( 89) =  0.21730115E+01
       angwei( 89)=  0.13840999E-02

       phi( 90)   =  0.16346878E+01
       theta( 90) =  0.16264950E+01
       angwei( 90)=  0.14476192E-02

       phi( 91)   =  0.43819103E+01
       theta( 91) =  0.11584808E+01
       angwei( 91)=  0.14493105E-02

       phi( 92)   =  0.26989343E+01
       theta( 92) =  0.72324848E+00
       angwei( 92)=  0.14307892E-02

       phi( 93)   =  0.19134238E+01
       theta( 93) =  0.24604578E+01
       angwei( 93)=  0.14202428E-02

       phi( 94)   =  0.45672379E+01
       theta( 94) =  0.75747317E+00
       angwei( 94)=  0.14273319E-02

       phi( 95)   =  0.85957870E-01
       theta( 95) =  0.14325124E+01
       angwei( 95)=  0.14268853E-02

       phi( 96)   =  0.20696199E+01
       theta( 96) =  0.21083033E+01
       angwei( 96)=  0.14337987E-02

       phi( 97)   =  0.39133890E+01
       theta( 97) =  0.23173547E+01
       angwei( 97)=  0.13075450E-02

       phi( 98)   =  0.51580501E+01
       theta( 98) =  0.63012630E+00
       angwei( 98)=  0.14116742E-02

       phi( 99)   =  0.32712562E+01
       theta( 99) =  0.80233192E+00
       angwei( 99)=  0.14069019E-02

       phi(100)   =  0.14979795E+01
       theta(100) =  0.94786787E+00
       angwei(100)=  0.14324078E-02

       phi(101)   =  0.13098097E+01
       theta(101) =  0.97850865E+00
       angwei(101)=  0.14690572E-02

       phi(102)   =  0.40301967E+01
       theta(102) =  0.22091589E+01
       angwei(102)=  0.14388489E-02

       phi(103)   =  0.50823960E+01
       theta(103) =  0.28311570E+01
       angwei(103)=  0.13556631E-02

       phi(104)   =  0.51067715E+01
       theta(104) =  0.49409881E+00
       angwei(104)=  0.14213886E-02

       phi(105)   =  0.32584207E+01
       theta(105) =  0.17019588E+01
       angwei(105)=  0.14297735E-02

       phi(106)   =  0.74624306E+00
       theta(106) =  0.17583426E+01
       angwei(106)=  0.14158072E-02

       phi(107)   =  0.47846389E+01
       theta(107) =  0.12373986E+01
       angwei(107)=  0.14247082E-02

       phi(108)   =  0.72751421E+00
       theta(108) =  0.14781761E+01
       angwei(108)=  0.14240218E-02

       phi(109)   =  0.56208477E+01
       theta(109) =  0.15986755E+01
       angwei(109)=  0.14318154E-02

       phi(110)   =  0.17641463E+01
       theta(110) =  0.54817367E+00
       angwei(110)=  0.14199641E-02

       phi(111)   =  0.52329955E+01
       theta(111) =  0.89279658E+00
       angwei(111)=  0.12998633E-02

       phi(112)   =  0.34222439E+00
       theta(112) =  0.14197528E+01
       angwei(112)=  0.14281275E-02

       phi(113)   =  0.25959427E+01
       theta(113) =  0.16969677E+01
       angwei(113)=  0.14144569E-02

       phi(114)   =  0.28936689E+01
       theta(114) =  0.13416491E+01
       angwei(114)=  0.14078066E-02

       phi(115)   =  0.10354477E+01
       theta(115) =  0.27886360E+01
       angwei(115)=  0.14514817E-02

       phi(116)   =  0.37045467E+01
       theta(116) =  0.28244216E+01
       angwei(116)=  0.13824876E-02

       phi(117)   =  0.16689831E+01
       theta(117) =  0.17630450E+01
       angwei(117)=  0.14281705E-02

       phi(118)   =  0.61253023E+01
       theta(118) =  0.16586015E+01
       angwei(118)=  0.15154738E-02

       phi(119)   =  0.32534313E+01
       theta(119) =  0.24852107E+01
       angwei(119)=  0.13956878E-02

       phi(120)   =  0.31691210E+01
       theta(120) =  0.19342898E+01
       angwei(120)=  0.14365948E-02

       phi(121)   =  0.56146235E+01
       theta(121) =  0.18827296E+01
       angwei(121)=  0.14357158E-02

       phi(122)   =  0.29942140E+00
       theta(122) =  0.67491221E+00
       angwei(122)=  0.14496899E-02

       phi(123)   =  0.35058556E+01
       theta(123) =  0.21543491E+01
       angwei(123)=  0.14197314E-02

       phi(124)   =  0.62272224E+01
       theta(124) =  0.13806978E+01
       angwei(124)=  0.14057586E-02

       phi(125)   =  0.55416131E+01
       theta(125) =  0.76678538E+00
       angwei(125)=  0.14261004E-02

       phi(126)   =  0.28200569E+01
       theta(126) =  0.95051306E+00
       angwei(126)=  0.14048755E-02

       phi(127)   =  0.24522555E+01
       theta(127) =  0.13097230E+01
       angwei(127)=  0.14167981E-02

       phi(128)   =  0.55997906E+01
       theta(128) =  0.14480397E+01
       angwei(128)=  0.14285019E-02

       phi(129)   =  0.34334626E+01
       theta(129) =  0.27248065E+01
       angwei(129)=  0.14103857E-02

       phi(130)   =  0.40004549E+01
       theta(130) =  0.61510253E+00
       angwei(130)=  0.14453228E-02

       phi(131)   =  0.13685530E+01
       theta(131) =  0.15027646E+01
       angwei(131)=  0.14313877E-02

       phi(132)   =  0.59254956E+01
       theta(132) =  0.28740766E+01
       angwei(132)=  0.14167853E-02

       phi(133)   =  0.17004578E+01
       theta(133) =  0.19047906E+01
       angwei(133)=  0.14349866E-02

       phi(134)   =  0.27419529E+01
       theta(134) =  0.58533931E+00
       angwei(134)=  0.14388036E-02

       phi(135)   =  0.52136874E+01
       theta(135) =  0.25712898E+01
       angwei(135)=  0.14243542E-02

       phi(136)   =  0.46767664E+01
       theta(136) =  0.87075466E+00
       angwei(136)=  0.14219325E-02

       phi(137)   =  0.23190885E+01
       theta(137) =  0.20111089E+01
       angwei(137)=  0.14219383E-02

       phi(138)   =  0.30675325E+01
       theta(138) =  0.27251694E+01
       angwei(138)=  0.14251870E-02

       phi(139)   =  0.10335904E+01
       theta(139) =  0.58343190E+00
       angwei(139)=  0.14915690E-02

       phi(140)   =  0.40867505E+01
       theta(140) =  0.17329091E+01
       angwei(140)=  0.14454941E-02

       phi(141)   =  0.84098065E+00
       theta(141) =  0.81290728E+00
       angwei(141)=  0.13782646E-02

       phi(142)   =  0.34699659E+01
       theta(142) =  0.15031841E+01
       angwei(142)=  0.14289212E-02

       phi(143)   =  0.31985310E+00
       theta(143) =  0.25898757E+01
       angwei(143)=  0.14573450E-02

       phi(144)   =  0.25535755E+01
       theta(144) =  0.25770013E+01
       angwei(144)=  0.14354691E-02

       phi(145)   =  0.16697046E+01
       theta(145) =  0.68106574E+00
       angwei(145)=  0.14225863E-02

       phi(146)   =  0.10222669E+00
       theta(146) =  0.15702744E+01
       angwei(146)=  0.14297988E-02

       phi(147)   =  0.31609969E+01
       theta(147) =  0.67714059E+00
       angwei(147)=  0.14564090E-02

       phi(148)   =  0.27314189E+00
       theta(148) =  0.27317252E+01
       angwei(148)=  0.14684921E-02

       phi(149)   =  0.82105291E+00
       theta(149) =  0.94557506E+00
       angwei(149)=  0.14913034E-02

       phi(150)   =  0.33652050E+01
       theta(150) =  0.16044481E+01
       angwei(150)=  0.14334671E-02

       phi(151)   =  0.18819206E+01
       theta(151) =  0.15602293E+01
       angwei(151)=  0.14356235E-02

       phi(152)   =  0.43393507E+01
       theta(152) =  0.12972610E+01
       angwei(152)=  0.13095445E-02

       phi(153)   =  0.34687195E+01
       theta(153) =  0.20174837E+01
       angwei(153)=  0.14289760E-02

       phi(154)   =  0.17283856E+01
       theta(154) =  0.25545602E+01
       angwei(154)=  0.14339051E-02

       phi(155)   =  0.42786078E+01
       theta(155) =  0.26294200E+01
       angwei(155)=  0.14542799E-02

       phi(156)   =  0.37614119E+01
       theta(156) =  0.53204685E+00
       angwei(156)=  0.14006625E-02

       phi(157)   =  0.22079604E+01
       theta(157) =  0.13659425E+01
       angwei(157)=  0.14313748E-02

       phi(158)   =  0.20753270E+00
       theta(158) =  0.20669827E+01
       angwei(158)=  0.14412845E-02

       phi(159)   =  0.42364273E+01
       theta(159) =  0.17572056E+01
       angwei(159)=  0.14115442E-02

       phi(160)   =  0.11616672E+01
       theta(160) =  0.19322293E+01
       angwei(160)=  0.14408785E-02

       phi(161)   =  0.44696908E+01
       theta(161) =  0.18867072E+01
       angwei(161)=  0.14384263E-02

       phi(162)   =  0.58395219E+01
       theta(162) =  0.14098372E+01
       angwei(162)=  0.14469652E-02

       phi(163)   =  0.14232147E+01
       theta(163) =  0.10730859E+01
       angwei(163)=  0.14303358E-02

       phi(164)   =  0.57020321E+01
       theta(164) =  0.13554912E+01
       angwei(164)=  0.14175592E-02

       phi(165)   =  0.35571539E+01
       theta(165) =  0.60862666E+00
       angwei(165)=  0.14059814E-02

       phi(166)   =  0.15043054E+01
       theta(166) =  0.56564802E+00
       angwei(166)=  0.14468855E-02

       phi(167)   =  0.34394777E+00
       theta(167) =  0.22923889E+01
       angwei(167)=  0.14310076E-02

       phi(168)   =  0.57633104E+01
       theta(168) =  0.11253381E+01
       angwei(168)=  0.14325170E-02

       phi(169)   =  0.42858434E+01
       theta(169) =  0.41527990E+00
       angwei(169)=  0.13850789E-02

       phi(170)   =  0.60815778E+01
       theta(170) =  0.82482928E+00
       angwei(170)=  0.13938381E-02

       phi(171)   =  0.35836089E+01
       theta(171) =  0.19123602E+01
       angwei(171)=  0.14366875E-02

       phi(172)   =  0.24686007E+01
       theta(172) =  0.23477755E+01
       angwei(172)=  0.14174937E-02

       phi(173)   =  0.51304960E+01
       theta(173) =  0.22085128E+01
       angwei(173)=  0.14434033E-02

       phi(174)   =  0.30378181E+00
       theta(174) =  0.97210622E+00
       angwei(174)=  0.14507516E-02

       phi(175)   =  0.57481849E+00
       theta(175) =  0.26380055E+01
       angwei(175)=  0.14069000E-02

       phi(176)   =  0.47168365E+01
       theta(176) =  0.13671939E+01
       angwei(176)=  0.14559741E-02

       phi(177)   =  0.33982904E+01
       theta(177) =  0.10832126E+01
       angwei(177)=  0.14115018E-02

       phi(178)   =  0.48044485E+00
       theta(178) =  0.75864393E+00
       angwei(178)=  0.14313898E-02

       phi(179)   =  0.62450399E+01
       theta(179) =  0.15106615E+01
       angwei(179)=  0.14023052E-02

       phi(180)   =  0.28489516E+01
       theta(180) =  0.23836615E+01
       angwei(180)=  0.14415939E-02

       phi(181)   =  0.12633353E+01
       theta(181) =  0.17232842E+01
       angwei(181)=  0.14532934E-02

       phi(182)   =  0.24256501E+01
       theta(182) =  0.21159823E+01
       angwei(182)=  0.14575879E-02

       phi(183)   =  0.21098959E+01
       theta(183) =  0.10813547E+01
       angwei(183)=  0.14342616E-02

       phi(184)   =  0.27533975E+01
       theta(184) =  0.17470139E+01
       angwei(184)=  0.14558969E-02

       phi(185)   =  0.27777343E+01
       theta(185) =  0.14182973E+01
       angwei(185)=  0.14408763E-02

       phi(186)   =  0.40573454E+01
       theta(186) =  0.19716325E+01
       angwei(186)=  0.14323455E-02

       phi(187)   =  0.15227398E+01
       theta(187) =  0.17110815E+01
       angwei(187)=  0.14395473E-02

       phi(188)   =  0.57478199E+01
       theta(188) =  0.34198520E+00
       angwei(188)=  0.14410201E-02

       phi(189)   =  0.23059685E+01
       theta(189) =  0.12653667E+01
       angwei(189)=  0.14282782E-02

       phi(190)   =  0.57734246E+01
       theta(190) =  0.22253642E+01
       angwei(190)=  0.14058809E-02

       phi(191)   =  0.46547861E+01
       theta(191) =  0.21269419E+01
       angwei(191)=  0.14318404E-02

       phi(192)   =  0.37414358E+01
       theta(192) =  0.24811065E+00
       angwei(192)=  0.12729131E-02

       phi(193)   =  0.61113548E+01
       theta(193) =  0.15122023E+01
       angwei(193)=  0.14894537E-02

       phi(194)   =  0.25052643E+01
       theta(194) =  0.35842916E+00
       angwei(194)=  0.14351369E-02

       phi(195)   =  0.97026181E+00
       theta(195) =  0.10213143E+01
       angwei(195)=  0.14955119E-02

       phi(196)   =  0.60908031E+01
       theta(196) =  0.13678761E+01
       angwei(196)=  0.14708336E-02

       phi(197)   =  0.48810716E+01
       theta(197) =  0.19934826E+01
       angwei(197)=  0.14351697E-02

       phi(198)   =  0.15354016E+01
       theta(198) =  0.11675013E+01
       angwei(198)=  0.14300131E-02

       phi(199)   =  0.36939971E+01
       theta(199) =  0.11765542E+01
       angwei(199)=  0.14125622E-02

       phi(200)   =  0.47664183E+00
       theta(200) =  0.99758834E-01
       angwei(200)=  0.14306620E-02

       phi(201)   =  0.17442094E+01
       theta(201) =  0.15337566E+01
       angwei(201)=  0.14290902E-02

       phi(202)   =  0.57587218E+01
       theta(202) =  0.70751840E+00
       angwei(202)=  0.14312655E-02

       phi(203)   =  0.41295652E+01
       theta(203) =  0.97088820E+00
       angwei(203)=  0.14297101E-02

       phi(204)   =  0.24034948E+01
       theta(204) =  0.15422043E+01
       angwei(204)=  0.14646936E-02

       phi(205)   =  0.47762737E+01
       theta(205) =  0.27345183E+01
       angwei(205)=  0.14010384E-02

       phi(206)   =  0.11410867E+01
       theta(206) =  0.17972955E+01
       angwei(206)=  0.14007820E-02

       phi(207)   =  0.56125050E+01
       theta(207) =  0.21672392E+01
       angwei(207)=  0.14632058E-02

       phi(208)   =  0.29138241E+01
       theta(208) =  0.31471652E+00
       angwei(208)=  0.14418926E-02

       phi(209)   =  0.10502627E+00
       theta(209) =  0.17105343E+01
       angwei(209)=  0.14363743E-02

       phi(210)   =  0.55843544E+01
       theta(210) =  0.62224436E+00
       angwei(210)=  0.14373658E-02

       phi(211)   =  0.57565775E+01
       theta(211) =  0.96493524E+00
       angwei(211)=  0.15728464E-02

       phi(212)   =  0.58075986E+01
       theta(212) =  0.12657385E+01
       angwei(212)=  0.14388956E-02

       phi(213)   =  0.33439243E+00
       theta(213) =  0.23871675E+00
       angwei(213)=  0.14440432E-02

       phi(214)   =  0.50750289E+01
       theta(214) =  0.18451633E+01
       angwei(214)=  0.14072228E-02

       phi(215)   =  0.38183844E+01
       theta(215) =  0.95982474E+00
       angwei(215)=  0.14461310E-02

       phi(216)   =  0.41585355E+01
       theta(216) =  0.83462447E+00
       angwei(216)=  0.14259199E-02

       phi(217)   =  0.50913978E+01
       theta(217) =  0.23438694E+01
       angwei(217)=  0.14423057E-02

       phi(218)   =  0.28678496E+01
       theta(218) =  0.16763947E+01
       angwei(218)=  0.14092498E-02

       phi(219)   =  0.11159049E+01
       theta(219) =  0.15222754E+01
       angwei(219)=  0.14130006E-02

       phi(220)   =  0.20178132E+01
       theta(220) =  0.15863721E+01
       angwei(220)=  0.14286144E-02

       phi(221)   =  0.18403593E+01
       theta(221) =  0.14208010E+01
       angwei(221)=  0.14211218E-02

       phi(222)   =  0.56563110E+01
       theta(222) =  0.12144716E+01
       angwei(222)=  0.14014371E-02

       phi(223)   =  0.41554036E+01
       theta(223) =  0.13372555E+00
       angwei(223)=  0.14255798E-02

       phi(224)   =  0.16083627E+01
       theta(224) =  0.10481833E+01
       angwei(224)=  0.14037213E-02

       phi(225)   =  0.42634149E+01
       theta(225) =  0.10619690E+01
       angwei(225)=  0.14373332E-02

       phi(226)   =  0.16042082E+01
       theta(226) =  0.14868171E+01
       angwei(226)=  0.13508130E-02

       phi(227)   =  0.50207791E+01
       theta(227) =  0.15982584E+01
       angwei(227)=  0.14302200E-02

       phi(228)   =  0.37376773E+01
       theta(228) =  0.19337329E+01
       angwei(228)=  0.14290127E-02

       phi(229)   =  0.40292053E+01
       theta(229) =  0.16024301E+01
       angwei(229)=  0.14194859E-02

       phi(230)   =  0.53196459E+01
       theta(230) =  0.20293250E+01
       angwei(230)=  0.14239224E-02

       phi(231)   =  0.23152726E+01
       theta(231) =  0.28517361E+01
       angwei(231)=  0.14157010E-02

       phi(232)   =  0.46853895E+01
       theta(232) =  0.23701050E+01
       angwei(232)=  0.14279082E-02

       phi(233)   =  0.14560748E+01
       theta(233) =  0.70715654E+00
       angwei(233)=  0.14268936E-02

       phi(234)   =  0.60396318E+01
       theta(234) =  0.10817840E+01
       angwei(234)=  0.14476447E-02

       phi(235)   =  0.56983519E+00
       theta(235) =  0.46127391E+00
       angwei(235)=  0.14045514E-02

       phi(236)   =  0.37966385E+01
       theta(236) =  0.16988823E+01
       angwei(236)=  0.14156494E-02

       phi(237)   =  0.45382242E+01
       theta(237) =  0.11327947E+01
       angwei(237)=  0.14259525E-02

       phi(238)   =  0.57473426E+01
       theta(238) =  0.20753314E+01
       angwei(238)=  0.14174954E-02

       phi(239)   =  0.43884187E+01
       theta(239) =  0.17643502E+01
       angwei(239)=  0.14474795E-02

       phi(240)   =  0.62555394E+01
       theta(240) =  0.24026043E+01
       angwei(240)=  0.14336564E-02

       phi(241)   =  0.27742093E+01
       theta(241) =  0.26638558E+01
       angwei(241)=  0.14331640E-02

       phi(242)   =  0.20829082E+01
       theta(242) =  0.30363017E+00
       angwei(242)=  0.14247312E-02

       phi(243)   =  0.11401064E+01
       theta(243) =  0.98883206E+00
       angwei(243)=  0.13748177E-02

       phi(244)   =  0.22856941E+01
       theta(244) =  0.58918995E+00
       angwei(244)=  0.14306197E-02

       phi(245)   =  0.75797223E-01
       theta(245) =  0.19880508E+01
       angwei(245)=  0.14333487E-02

       phi(246)   =  0.21139371E+01
       theta(246) =  0.14739276E+01
       angwei(246)=  0.14310253E-02

       phi(247)   =  0.10372509E+01
       theta(247) =  0.20286078E+01
       angwei(247)=  0.14003088E-02

       phi(248)   =  0.89925760E+00
       theta(248) =  0.19774069E+01
       angwei(248)=  0.14394532E-02

       phi(249)   =  0.53784523E+01
       theta(249) =  0.29497895E+01
       angwei(249)=  0.14106609E-02

       phi(250)   =  0.22297725E+00
       theta(250) =  0.14957281E+01
       angwei(250)=  0.14299448E-02

       phi(251)   =  0.13800980E+01
       theta(251) =  0.16491755E+01
       angwei(251)=  0.14019426E-02

       phi(252)   =  0.46874080E+01
       theta(252) =  0.17562010E+01
       angwei(252)=  0.14350482E-02

       phi(253)   =  0.58356208E+00
       theta(253) =  0.97319323E+00
       angwei(253)=  0.13593758E-02

       phi(254)   =  0.22613497E-01
       theta(254) =  0.10060606E+01
       angwei(254)=  0.14419234E-02

       phi(255)   =  0.85467702E+00
       theta(255) =  0.22886877E+01
       angwei(255)=  0.14457622E-02

       phi(256)   =  0.18487558E-01
       theta(256) =  0.47159988E+00
       angwei(256)=  0.14408096E-02

       phi(257)   =  0.31467719E+01
       theta(257) =  0.43996504E+00
       angwei(257)=  0.14076688E-02

       phi(258)   =  0.44143066E+00
       theta(258) =  0.10475791E+01
       angwei(258)=  0.14400799E-02

       phi(259)   =  0.61555200E+01
       theta(259) =  0.21718516E+01
       angwei(259)=  0.14914206E-02

       phi(260)   =  0.22415335E+01
       theta(260) =  0.18778553E+01
       angwei(260)=  0.14271538E-02

       phi(261)   =  0.22543285E+01
       theta(261) =  0.15037255E+01
       angwei(261)=  0.14111444E-02

       phi(262)   =  0.42054124E+01
       theta(262) =  0.22336400E+01
       angwei(262)=  0.14435990E-02

       phi(263)   =  0.75702053E+00
       theta(263) =  0.18979386E+01
       angwei(263)=  0.14209149E-02

       phi(264)   =  0.55906758E+01
       theta(264) =  0.17298098E+01
       angwei(264)=  0.14444101E-02

       phi(265)   =  0.25078869E+01
       theta(265) =  0.79094523E+00
       angwei(265)=  0.14306750E-02

       phi(266)   =  0.84524643E+00
       theta(266) =  0.11058514E+01
       angwei(266)=  0.13901920E-02

       phi(267)   =  0.46269450E+01
       theta(267) =  0.18855658E+01
       angwei(267)=  0.14316423E-02

       phi(268)   =  0.45995355E+01
       theta(268) =  0.10009903E+01
       angwei(268)=  0.14324927E-02

       phi(269)   =  0.49473724E+00
       theta(269) =  0.22107072E+01
       angwei(269)=  0.14357894E-02

       phi(270)   =  0.28091307E+01
       theta(270) =  0.10807991E+01
       angwei(270)=  0.13400286E-02

       phi(271)   =  0.42189107E+01
       theta(271) =  0.13550603E+01
       angwei(271)=  0.15132560E-02

       phi(272)   =  0.26393111E+01
       theta(272) =  0.18253615E+01
       angwei(272)=  0.13278092E-02

       phi(273)   =  0.40919271E+01
       theta(273) =  0.23525960E+01
       angwei(273)=  0.14069632E-02

       phi(274)   =  0.27309299E+01
       theta(274) =  0.12894034E+01
       angwei(274)=  0.14424904E-02

       phi(275)   =  0.35019875E+00
       theta(275) =  0.19945741E+01
       angwei(275)=  0.14292259E-02

       phi(276)   =  0.24651895E+01
       theta(276) =  0.27162597E+01
       angwei(276)=  0.14291662E-02

       phi(277)   =  0.41314011E+01
       theta(277) =  0.20987251E+01
       angwei(277)=  0.14187782E-02

       phi(278)   =  0.20514374E+01
       theta(278) =  0.94237560E+00
       angwei(278)=  0.14313535E-02

       phi(279)   =  0.26847150E+01
       theta(279) =  0.23095977E+01
       angwei(279)=  0.14453395E-02

       phi(280)   =  0.46504812E+01
       theta(280) =  0.28797944E+01
       angwei(280)=  0.14029525E-02

       phi(281)   =  0.44118710E+01
       theta(281) =  0.27548540E+01
       angwei(281)=  0.14218985E-02

       phi(282)   =  0.61261230E+01
       theta(282) =  0.22990680E+01
       angwei(282)=  0.12818399E-02

       phi(283)   =  0.19803243E+01
       theta(283) =  0.80837518E+00
       angwei(283)=  0.14365111E-02

       phi(284)   =  0.87479901E+00
       theta(284) =  0.67998135E+00
       angwei(284)=  0.13051174E-02

       phi(285)   =  0.54989753E+01
       theta(285) =  0.15441504E+01
       angwei(285)=  0.12906935E-02

       phi(286)   =  0.59763904E+01
       theta(286) =  0.14584100E+01
       angwei(286)=  0.14035688E-02

       phi(287)   =  0.19165926E+01
       theta(287) =  0.16990209E+01
       angwei(287)=  0.14424355E-02

       phi(288)   =  0.65805566E+00
       theta(288) =  0.85479569E+00
       angwei(288)=  0.15454573E-02

       phi(289)   =  0.30857525E+01
       theta(289) =  0.10307834E+01
       angwei(289)=  0.13089490E-02

       phi(290)   =  0.42236595E+01
       theta(290) =  0.19887866E+01
       angwei(290)=  0.14303226E-02

       phi(291)   =  0.28236926E+01
       theta(291) =  0.15464159E+01
       angwei(291)=  0.14486561E-02

       phi(292)   =  0.25502403E+01
       theta(292) =  0.15689306E+01
       angwei(292)=  0.14030627E-02

       phi(293)   =  0.49270916E+01
       theta(293) =  0.18633194E+01
       angwei(293)=  0.14397879E-02

       phi(294)   =  0.99892098E+00
       theta(294) =  0.30250111E+01
       angwei(294)=  0.12978980E-02

       phi(295)   =  0.51160293E+01
       theta(295) =  0.17142231E+01
       angwei(295)=  0.14227778E-02

       phi(296)   =  0.11273929E+01
       theta(296) =  0.16599492E+01
       angwei(296)=  0.14062135E-02

       phi(297)   =  0.33750186E+01
       theta(297) =  0.22623489E+01
       angwei(297)=  0.14503959E-02

       phi(298)   =  0.53734622E+01
       theta(298) =  0.54856861E+00
       angwei(298)=  0.14491816E-02

       phi(299)   =  0.38973794E+01
       theta(299) =  0.19542016E+01
       angwei(299)=  0.14273108E-02

       phi(300)   =  0.17326663E+01
       theta(300) =  0.20459023E+01
       angwei(300)=  0.14611213E-02

       phi(301)   =  0.35732108E+00
       theta(301) =  0.18510175E+01
       angwei(301)=  0.14325534E-02

       phi(302)   =  0.39951966E+01
       theta(302) =  0.26004539E+01
       angwei(302)=  0.14469345E-02

       phi(303)   =  0.19090688E+00
       theta(303) =  0.12098907E+01
       angwei(303)=  0.14275865E-02

       phi(304)   =  0.36214397E+01
       theta(304) =  0.29722095E+01
       angwei(304)=  0.14313913E-02

       phi(305)   =  0.37953944E+01
       theta(305) =  0.37956074E+00
       angwei(305)=  0.15390803E-02

       phi(306)   =  0.47382757E+00
       theta(306) =  0.14876498E+01
       angwei(306)=  0.14271811E-02

       phi(307)   =  0.11071491E+01
       theta(307) =  0.12463511E+01
       angwei(307)=  0.13897103E-02

       phi(308)   =  0.42310638E+01
       theta(308) =  0.12022878E+01
       angwei(308)=  0.14460871E-02

       phi(309)   =  0.24128325E+01
       theta(309) =  0.11696250E+01
       angwei(309)=  0.14308210E-02

       phi(310)   =  0.44915870E+00
       theta(310) =  0.11984411E+01
       angwei(310)=  0.14328293E-02

       phi(311)   =  0.54591260E+01
       theta(311) =  0.25233803E+01
       angwei(311)=  0.14476196E-02

       phi(312)   =  0.13635889E+01
       theta(312) =  0.12018212E+01
       angwei(312)=  0.14266835E-02

       phi(313)   =  0.97591889E+00
       theta(313) =  0.13236827E+01
       angwei(313)=  0.14320284E-02

       phi(314)   =  0.71478343E+00
       theta(314) =  0.24003332E+01
       angwei(314)=  0.14480475E-02

       phi(315)   =  0.52237134E+01
       theta(315) =  0.17951726E+01
       angwei(315)=  0.14602470E-02

       phi(316)   =  0.48897023E+01
       theta(316) =  0.23681352E+01
       angwei(316)=  0.14363804E-02

       phi(317)   =  0.43882113E+01
       theta(317) =  0.15294087E+01
       angwei(317)=  0.12836265E-02

       phi(318)   =  0.31201782E+01
       theta(318) =  0.16516768E+01
       angwei(318)=  0.14285111E-02

       phi(319)   =  0.58699322E+01
       theta(319) =  0.16979020E+01
       angwei(319)=  0.14584924E-02

       phi(320)   =  0.51423144E-01
       theta(320) =  0.21273928E+01
       angwei(320)=  0.14336696E-02

       phi(321)   =  0.75040025E+00
       theta(321) =  0.56510425E+00
       angwei(321)=  0.14245670E-02

       phi(322)   =  0.61877861E+01
       theta(322) =  0.11092085E+01
       angwei(322)=  0.14007958E-02

       phi(323)   =  0.44291058E+01
       theta(323) =  0.10235684E+01
       angwei(323)=  0.14354188E-02

       phi(324)   =  0.31465867E+01
       theta(324) =  0.17936125E+01
       angwei(324)=  0.14314285E-02

       phi(325)   =  0.53799806E+01
       theta(325) =  0.40381297E+00
       angwei(325)=  0.14354588E-02

       phi(326)   =  0.24869716E+01
       theta(326) =  0.19830668E+01
       angwei(326)=  0.14491073E-02

       phi(327)   =  0.14422017E+01
       theta(327) =  0.19272357E+01
       angwei(327)=  0.14239576E-02

       phi(328)   =  0.56281085E+01
       theta(328) =  0.23117735E+01
       angwei(328)=  0.14116634E-02

       phi(329)   =  0.49433460E+01
       theta(329) =  0.14765370E+01
       angwei(329)=  0.14226101E-02

       phi(330)   =  0.38973196E+01
       theta(330) =  0.24492033E+01
       angwei(330)=  0.15650251E-02

       phi(331)   =  0.14413244E+01
       theta(331) =  0.27599738E+01
       angwei(331)=  0.14249188E-02

       phi(332)   =  0.19030166E+01
       theta(332) =  0.20782800E+01
       angwei(332)=  0.14229818E-02

       phi(333)   =  0.54614553E+01
       theta(333) =  0.21023030E+01
       angwei(333)=  0.14348985E-02

       phi(334)   =  0.22751191E+01
       theta(334) =  0.26047311E+01
       angwei(334)=  0.14336207E-02

       phi(335)   =  0.46243665E+00
       theta(335) =  0.13437709E+01
       angwei(335)=  0.14364955E-02

       phi(336)   =  0.60682993E+01
       theta(336) =  0.12248591E+01
       angwei(336)=  0.14582901E-02

       phi(337)   =  0.61920351E+00
       theta(337) =  0.18365624E+01
       angwei(337)=  0.14328395E-02

       phi(338)   =  0.41959758E+01
       theta(338) =  0.69848633E+00
       angwei(338)=  0.14246255E-02

       phi(339)   =  0.29799817E+01
       theta(339) =  0.93141145E+00
       angwei(339)=  0.14697644E-02

       phi(340)   =  0.50498176E+01
       theta(340) =  0.87131166E+00
       angwei(340)=  0.14277005E-02

       phi(341)   =  0.23475642E+01
       theta(341) =  0.22389965E+01
       angwei(341)=  0.14447100E-02

       phi(342)   =  0.43121624E+01
       theta(342) =  0.16386267E+01
       angwei(342)=  0.14452413E-02

       phi(343)   =  0.53612080E+01
       theta(343) =  0.11394585E+01
       angwei(343)=  0.13937283E-02

       phi(344)   =  0.30274241E+01
       theta(344) =  0.18824122E+01
       angwei(344)=  0.14364505E-02

       phi(345)   =  0.26619031E+01
       theta(345) =  0.10093025E+01
       angwei(345)=  0.14236324E-02

       phi(346)   =  0.55173864E+01
       theta(346) =  0.11593233E+01
       angwei(346)=  0.14810262E-02

       phi(347)   =  0.21980937E+01
       theta(347) =  0.17444750E+01
       angwei(347)=  0.14464136E-02

       phi(348)   =  0.30133197E+01
       theta(348) =  0.55899096E+00
       angwei(348)=  0.14106212E-02

       phi(349)   =  0.48073926E+01
       theta(349) =  0.48552734E+00
       angwei(349)=  0.14358648E-02

       phi(350)   =  0.40166478E+01
       theta(350) =  0.47473758E+00
       angwei(350)=  0.13341465E-02

       phi(351)   =  0.44847636E+01
       theta(351) =  0.12671937E+01
       angwei(351)=  0.14187500E-02

       phi(352)   =  0.50210757E+01
       theta(352) =  0.35788515E+00
       angwei(352)=  0.14065154E-02

       phi(353)   =  0.29569080E+01
       theta(353) =  0.10980997E+01
       angwei(353)=  0.15141660E-02

       phi(354)   =  0.49127688E+01
       theta(354) =  0.26119456E+01
       angwei(354)=  0.14779354E-02

       phi(355)   =  0.52576661E+01
       theta(355) =  0.24367771E+01
       angwei(355)=  0.14088656E-02

       phi(356)   =  0.31884420E+01
       theta(356) =  0.20741112E+01
       angwei(356)=  0.14262176E-02

       phi(357)   =  0.35482857E+01
       theta(357) =  0.11239818E+01
       angwei(357)=  0.14232772E-02

       phi(358)   =  0.53644891E+01
       theta(358) =  0.15405257E+01
       angwei(358)=  0.14761802E-02

       phi(359)   =  0.19908726E+01
       theta(359) =  0.25961878E+01
       angwei(359)=  0.14386257E-02

       phi(360)   =  0.17902092E+01
       theta(360) =  0.21754942E+01
       angwei(360)=  0.14025021E-02

       phi(361)   =  0.84718667E-01
       theta(361) =  0.61218423E+00
       angwei(361)=  0.14311698E-02

       phi(362)   =  0.34689646E+01
       theta(362) =  0.83910203E+00
       angwei(362)=  0.14414849E-02

       phi(363)   =  0.51703310E+01
       theta(363) =  0.13577095E+01
       angwei(363)=  0.14272069E-02

       phi(364)   =  0.61285353E+00
       theta(364) =  0.16953936E+01
       angwei(364)=  0.14418586E-02

       phi(365)   =  0.26125503E+01
       theta(365) =  0.13591708E+01
       angwei(365)=  0.14453956E-02

       phi(366)   =  0.94027817E+00
       theta(366) =  0.21375482E+01
       angwei(366)=  0.15602508E-02

       phi(367)   =  0.12387656E+01
       theta(367) =  0.12760528E+01
       angwei(367)=  0.14892892E-02

       phi(368)   =  0.54007015E+01
       theta(368) =  0.86173153E+00
       angwei(368)=  0.15786700E-02

       phi(369)   =  0.43863525E+01
       theta(369) =  0.22436972E+01
       angwei(369)=  0.14281876E-02

       phi(370)   =  0.10175459E+01
       theta(370) =  0.18846152E+01
       angwei(370)=  0.14309689E-02

       phi(371)   =  0.49722075E+01
       theta(371) =  0.17312418E+01
       angwei(371)=  0.14270146E-02

       phi(372)   =  0.18837428E+01
       theta(372) =  0.42090312E+00
       angwei(372)=  0.14248034E-02

       phi(373)   =  0.20326297E+01
       theta(373) =  0.55536509E+00
       angwei(373)=  0.14326109E-02

       phi(374)   =  0.58793831E+01
       theta(374) =  0.25065715E+01
       angwei(374)=  0.14065208E-02

       phi(375)   =  0.47684507E+01
       theta(375) =  0.22982448E+00
       angwei(375)=  0.14562866E-02

       phi(376)   =  0.14907203E+01
       theta(376) =  0.14210914E+01
       angwei(376)=  0.15206340E-02

       phi(377)   =  0.43111963E+01
       theta(377) =  0.18802044E+01
       angwei(377)=  0.14277794E-02

       phi(378)   =  0.25131414E+01
       theta(378) =  0.93699342E+00
       angwei(378)=  0.14224066E-02

       phi(379)   =  0.59759483E+01
       theta(379) =  0.94030881E+00
       angwei(379)=  0.14393063E-02

       phi(380)   =  0.56974645E+01
       theta(380) =  0.25962946E+01
       angwei(380)=  0.14274591E-02

       phi(381)   =  0.18104858E+01
       theta(381) =  0.18071423E+01
       angwei(381)=  0.14166650E-02

       phi(382)   =  0.52846355E+01
       theta(382) =  0.16512506E+01
       angwei(382)=  0.15473595E-02

       phi(383)   =  0.21073582E+01
       theta(383) =  0.27272644E+01
       angwei(383)=  0.14295267E-02

       phi(384)   =  0.22207234E+01
       theta(384) =  0.44546643E+00
       angwei(384)=  0.14363035E-02

       phi(385)   =  0.54542994E+01
       theta(385) =  0.26647117E+01
       angwei(385)=  0.14461564E-02

       phi(386)   =  0.32277770E+01
       theta(386) =  0.15613447E+01
       angwei(386)=  0.14332799E-02

       phi(387)   =  0.18946743E+01
       theta(387) =  0.67586964E+00
       angwei(387)=  0.14279035E-02

       phi(388)   =  0.52404599E+01
       theta(388) =  0.14848645E+01
       angwei(388)=  0.14144962E-02

       phi(389)   =  0.49501014E+01
       theta(389) =  0.22379889E+01
       angwei(389)=  0.14354091E-02

       phi(390)   =  0.26608891E+01
       theta(390) =  0.14873067E+01
       angwei(390)=  0.14457791E-02

       phi(391)   =  0.36770506E+01
       theta(391) =  0.10344514E+01
       angwei(391)=  0.14255462E-02

       phi(392)   =  0.29804711E+01
       theta(392) =  0.16010789E+01
       angwei(392)=  0.14176888E-02

       phi(393)   =  0.42672524E+01
       theta(393) =  0.14803891E+01
       angwei(393)=  0.15247550E-02

       phi(394)   =  0.45735703E+01
       theta(394) =  0.24959750E+01
       angwei(394)=  0.13428723E-02

       phi(395)   =  0.10009441E+01
       theta(395) =  0.88426000E+00
       angwei(395)=  0.14629283E-02

       phi(396)   =  0.16015631E+01
       theta(396) =  0.29344261E+00
       angwei(396)=  0.14437586E-02

       phi(397)   =  0.10244695E+01
       theta(397) =  0.21150395E+00
       angwei(397)=  0.14279189E-02

       phi(398)   =  0.16911321E+01
       theta(398) =  0.11720835E+01
       angwei(398)=  0.14601209E-02

       phi(399)   =  0.19486376E+01
       theta(399) =  0.10549579E+01
       angwei(399)=  0.14194081E-02

       phi(400)   =  0.40929127E+01
       theta(400) =  0.12461292E+01
       angwei(400)=  0.14329416E-02

       phi(401)   =  0.62581062E+01
       theta(401) =  0.16391308E+01
       angwei(401)=  0.13681385E-02

       phi(402)   =  0.54534559E+01
       theta(402) =  0.14137331E+01
       angwei(402)=  0.14569453E-02

       phi(403)   =  0.55644422E+01
       theta(403) =  0.92206025E+00
       angwei(403)=  0.13990415E-02

       phi(404)   =  0.11873386E+01
       theta(404) =  0.86333561E+00
       angwei(404)=  0.13871206E-02

       phi(405)   =  0.27552960E+01
       theta(405) =  0.20370193E+01
       angwei(405)=  0.14433416E-02

       phi(406)   =  0.47216315E+01
       theta(406) =  0.20036821E+01
       angwei(406)=  0.14280491E-02

       phi(407)   =  0.16898608E+01
       theta(407) =  0.92647356E+00
       angwei(407)=  0.14131286E-02

       phi(408)   =  0.51973052E+01
       theta(408) =  0.11225855E+01
       angwei(408)=  0.14468201E-02

       phi(409)   =  0.39769385E+01
       theta(409) =  0.14699109E+01
       angwei(409)=  0.14185011E-02

       phi(410)   =  0.46049390E+01
       theta(410) =  0.16347784E+01
       angwei(410)=  0.14385099E-02

       phi(411)   =  0.56615152E+01
       theta(411) =  0.24509807E+01
       angwei(411)=  0.14588651E-02

       phi(412)   =  0.43037643E+01
       theta(412) =  0.92491817E+00
       angwei(412)=  0.14315611E-02

       phi(413)   =  0.42277751E+01
       theta(413) =  0.28593963E+00
       angwei(413)=  0.14920069E-02

       phi(414)   =  0.32277424E+01
       theta(414) =  0.28563645E+01
       angwei(414)=  0.14464766E-02

       phi(415)   =  0.44793401E+01
       theta(415) =  0.23697152E+01
       angwei(415)=  0.14230440E-02

       phi(416)   =  0.20123355E+01
       theta(416) =  0.11913288E+01
       angwei(416)=  0.14367239E-02

       phi(417)   =  0.40828314E+01
       theta(417) =  0.13785076E+01
       angwei(417)=  0.13101019E-02

       phi(418)   =  0.33983223E+01
       theta(418) =  0.17442786E+01
       angwei(418)=  0.14363712E-02

       phi(419)   =  0.25329368E+01
       theta(419) =  0.10865189E+01
       angwei(419)=  0.14246020E-02

       phi(420)   =  0.29230769E+01
       theta(420) =  0.68412250E+00
       angwei(420)=  0.14124115E-02

       phi(421)   =  0.46335373E+01
       theta(421) =  0.12486376E+01
       angwei(421)=  0.14205059E-02

       phi(422)   =  0.38463230E+01
       theta(422) =  0.21885452E+01
       angwei(422)=  0.14284118E-02

       phi(423)   =  0.61582823E+01
       theta(423) =  0.30134978E+01
       angwei(423)=  0.15285515E-02

       phi(424)   =  0.66101241E+00
       theta(424) =  0.22659955E+01
       angwei(424)=  0.14327555E-02

       phi(425)   =  0.48021426E+01
       theta(425) =  0.14878632E+01
       angwei(425)=  0.14278215E-02

       phi(426)   =  0.36145308E+01
       theta(426) =  0.74978673E+00
       angwei(426)=  0.14381639E-02

       phi(427)   =  0.27703438E+01
       theta(427) =  0.18948412E+01
       angwei(427)=  0.14230757E-02

       phi(428)   =  0.11409125E+01
       theta(428) =  0.24055047E+01
       angwei(428)=  0.14413357E-02

       phi(429)   =  0.27599649E+01
       theta(429) =  0.28122356E+01
       angwei(429)=  0.14298749E-02

       phi(430)   =  0.19555250E+00
       theta(430) =  0.22107065E+01
       angwei(430)=  0.14259184E-02

       phi(431)   =  0.42479963E+01
       theta(431) =  0.56161964E+00
       angwei(431)=  0.14198254E-02

       phi(432)   =  0.23000684E+01
       theta(432) =  0.16345419E+01
       angwei(432)=  0.14041009E-02

       phi(433)   =  0.54572582E+01
       theta(433) =  0.23837395E+01
       angwei(433)=  0.14275002E-02

       phi(434)   =  0.13471869E+01
       theta(434) =  0.28966792E+01
       angwei(434)=  0.14322373E-02

       phi(435)   =  0.24531462E+01
       theta(435) =  0.16770039E+01
       angwei(435)=  0.14588326E-02

       phi(436)   =  0.43909488E+01
       theta(436) =  0.20034516E+01
       angwei(436)=  0.14284679E-02

       phi(437)   =  0.62051840E+01
       theta(437) =  0.20430145E+01
       angwei(437)=  0.13573597E-02

       phi(438)   =  0.47957602E+01
       theta(438) =  0.24984238E+01
       angwei(438)=  0.13969282E-02

       phi(439)   =  0.54850454E+01
       theta(439) =  0.18227590E+01
       angwei(439)=  0.14376859E-02

       phi(440)   =  0.13510191E+01
       theta(440) =  0.13605390E+01
       angwei(440)=  0.14062502E-02

       phi(441)   =  0.11909195E+00
       theta(441) =  0.75782728E+00
       angwei(441)=  0.14293730E-02

       phi(442)   =  0.43017063E+01
       theta(442) =  0.21164651E+01
       angwei(442)=  0.14359591E-02

       phi(443)   =  0.56484871E+01
       theta(443) =  0.48216826E+00
       angwei(443)=  0.14303393E-02

       phi(444)   =  0.23871017E+01
       theta(444) =  0.18893808E+01
       angwei(444)=  0.13168593E-02

       phi(445)   =  0.61136246E+01
       theta(445) =  0.57974929E+00
       angwei(445)=  0.14414493E-02

       phi(446)   =  0.23014672E+00
       theta(446) =  0.17821351E+01
       angwei(446)=  0.14219364E-02

       phi(447)   =  0.19762847E+01
       theta(447) =  0.29722757E+01
       angwei(447)=  0.14351010E-02

       phi(448)   =  0.29550397E+00
       theta(448) =  0.37943617E+00
       angwei(448)=  0.14225354E-02

       phi(449)   =  0.30430615E+01
       theta(449) =  0.23049495E+01
       angwei(449)=  0.14367941E-02

       phi(450)   =  0.61768198E+00
       theta(450) =  0.27719712E+01
       angwei(450)=  0.14039931E-02

       phi(451)   =  0.17759906E+01
       theta(451) =  0.12832919E+01
       angwei(451)=  0.13923714E-02

       phi(452)   =  0.11756015E+01
       theta(452) =  0.34078100E+00
       angwei(452)=  0.14153994E-02

       phi(453)   =  0.31218781E+01
       theta(453) =  0.11608856E+01
       angwei(453)=  0.14163393E-02

       phi(454)   =  0.46283035E+01
       theta(454) =  0.36652914E+00
       angwei(454)=  0.14040928E-02

       phi(455)   =  0.33969204E+01
       theta(455) =  0.69987893E+00
       angwei(455)=  0.14144792E-02

       phi(456)   =  0.38847029E+01
       theta(456) =  0.15694035E+01
       angwei(456)=  0.14742485E-02

       phi(457)   =  0.30377038E+01
       theta(457) =  0.24425309E+01
       angwei(457)=  0.14499451E-02

       phi(458)   =  0.37522118E+01
       theta(458) =  0.15717243E+01
       angwei(458)=  0.14011149E-02

       phi(459)   =  0.50885949E+01
       theta(459) =  0.14715016E+01
       angwei(459)=  0.14380616E-02

       phi(460)   =  0.57012906E+01
       theta(460) =  0.83953291E+00
       angwei(460)=  0.12918019E-02

       phi(461)   =  0.15443320E+01
       theta(461) =  0.42824858E+00
       angwei(461)=  0.14394724E-02

       phi(462)   =  0.53738594E+01
       theta(462) =  0.69702822E+00
       angwei(462)=  0.14331281E-02

       phi(463)   =  0.19669857E+01
       theta(463) =  0.22158589E+01
       angwei(463)=  0.14375967E-02

       phi(464)   =  0.29613683E+00
       theta(464) =  0.52607238E+00
       angwei(464)=  0.14335935E-02

       phi(465)   =  0.41066971E+01
       theta(465) =  0.11077297E+01
       angwei(465)=  0.14092609E-02

       phi(466)   =  0.37170246E+01
       theta(466) =  0.14469254E+01
       angwei(466)=  0.13383316E-02

       phi(467)   =  0.20445292E+01
       theta(467) =  0.23527415E+01
       angwei(467)=  0.14292718E-02

       phi(468)   =  0.34336288E+01
       theta(468) =  0.18811939E+01
       angwei(468)=  0.14377101E-02

       phi(469)   =  0.47653103E+01
       theta(469) =  0.22463646E+01
       angwei(469)=  0.14284020E-02

       phi(470)   =  0.12884139E+01
       theta(470) =  0.18616039E+01
       angwei(470)=  0.14493981E-02

       phi(471)   =  0.13131108E+01
       theta(471) =  0.21525939E+01
       angwei(471)=  0.13765151E-02

       phi(472)   =  0.30585737E+01
       theta(472) =  0.79886860E+00
       angwei(472)=  0.14576918E-02

       phi(473)   =  0.16110834E+01
       theta(473) =  0.13005909E+01
       angwei(473)=  0.15167029E-02

       phi(474)   =  0.35508785E+01
       theta(474) =  0.22928624E+01
       angwei(474)=  0.13838719E-02

       phi(475)   =  0.53220916E+01
       theta(475) =  0.13791350E+01
       angwei(475)=  0.14054454E-02

       phi(476)   =  0.50315685E+01
       theta(476) =  0.24800267E+01
       angwei(476)=  0.14304026E-02

       phi(477)   =  0.21456289E+01
       theta(477) =  0.24818335E+01
       angwei(477)=  0.14263762E-02

       phi(478)   =  0.54578090E+01
       theta(478) =  0.22426689E+01
       angwei(478)=  0.14406610E-02

       phi(479)   =  0.47802258E+01
       theta(479) =  0.73994112E+00
       angwei(479)=  0.14354571E-02

       phi(480)   =  0.31706223E+01
       theta(480) =  0.91309398E+00
       angwei(480)=  0.14158532E-02

       phi(481)   =  0.12468325E+00
       theta(481) =  0.25022867E+01
       angwei(481)=  0.14045945E-02

       phi(482)   =  0.46157331E+01
       theta(482) =  0.30191925E+01
       angwei(482)=  0.14656375E-02

       phi(483)   =  0.30087221E+01
       theta(483) =  0.17401700E+01
       angwei(483)=  0.14277359E-02

       phi(484)   =  0.87801874E+00
       theta(484) =  0.26651385E+01
       angwei(484)=  0.14267355E-02

       phi(485)   =  0.34440825E+01
       theta(485) =  0.13614568E+01
       angwei(485)=  0.14292999E-02

       phi(486)   =  0.28025663E+01
       theta(486) =  0.25186775E+01
       angwei(486)=  0.14114676E-02

       phi(487)   =  0.51267214E-01
       theta(487) =  0.11518378E+01
       angwei(487)=  0.14362916E-02

       phi(488)   =  0.23500414E+01
       theta(488) =  0.87768453E+00
       angwei(488)=  0.14248228E-02

       phi(489)   =  0.62455291E+00
       theta(489) =  0.19779546E+01
       angwei(489)=  0.14274061E-02

       phi(490)   =  0.21509836E+01
       theta(490) =  0.22409163E+01
       angwei(490)=  0.14215575E-02

       phi(491)   =  0.61453772E+01
       theta(491) =  0.96612263E+00
       angwei(491)=  0.14329579E-02

       phi(492)   =  0.39884269E+01
       theta(492) =  0.75378019E+00
       angwei(492)=  0.14313515E-02

       phi(493)   =  0.59501886E+00
       theta(493) =  0.14123455E+01
       angwei(493)=  0.14287485E-02

       phi(494)   =  0.78863138E+00
       theta(494) =  0.21663592E+01
       angwei(494)=  0.13081025E-02

       phi(495)   =  0.14079399E+01
       theta(495) =  0.17887923E+01
       angwei(495)=  0.14366105E-02

       phi(496)   =  0.10035330E+01
       theta(496) =  0.17432408E+01
       angwei(496)=  0.14377482E-02

       phi(497)   =  0.91848545E-01
       theta(497) =  0.18495408E+01
       angwei(497)=  0.14463478E-02

       phi(498)   =  0.84306091E+00
       theta(498) =  0.12585239E+01
       angwei(498)=  0.14098776E-02

       phi(499)   =  0.30439312E+01
       theta(499) =  0.13756356E+01
       angwei(499)=  0.14591672E-02

       phi(500)   =  0.29154406E+01
       theta(500) =  0.31052418E+01
       angwei(500)=  0.14531678E-02

       phi(501)   =  0.17051401E+01
       theta(501) =  0.24128301E+01
       angwei(501)=  0.14400225E-02

       phi(502)   =  0.37902975E+01
       theta(502) =  0.67496204E+00
       angwei(502)=  0.14349465E-02

       phi(503)   =  0.11020994E+01
       theta(503) =  0.11147425E+01
       angwei(503)=  0.13510743E-02

       phi(504)   =  0.88217020E+00
       theta(504) =  0.18282648E+01
       angwei(504)=  0.14481343E-02

       phi(505)   =  0.71380430E+00
       theta(505) =  0.11902037E+01
       angwei(505)=  0.14619579E-02

       phi(506)   =  0.41789637E+01
       theta(506) =  0.16235517E+01
       angwei(506)=  0.14240327E-02

       phi(507)   =  0.59792953E+01
       theta(507) =  0.20810947E+00
       angwei(507)=  0.14235998E-02

       phi(508)   =  0.25955100E+01
       theta(508) =  0.20932755E+01
       angwei(508)=  0.14080334E-02

       phi(509)   =  0.32995996E+01
       theta(509) =  0.56413239E+00
       angwei(509)=  0.14596641E-02

       phi(510)   =  0.17478095E+01
       theta(510) =  0.16144346E+00
       angwei(510)=  0.14407421E-02

       phi(511)   =  0.35428643E+01
       theta(511) =  0.17790889E+01
       angwei(511)=  0.14285496E-02

       phi(512)   =  0.58997231E+01
       theta(512) =  0.19865141E+01
       angwei(512)=  0.15650293E-02

       phi(513)   =  0.18581340E+01
       theta(513) =  0.11692768E+01
       angwei(513)=  0.14155287E-02

       phi(514)   =  0.15793765E+01
       theta(514) =  0.81738865E+00
       angwei(514)=  0.14344337E-02

       phi(515)   =  0.58906207E+01
       theta(515) =  0.21253586E+01
       angwei(515)=  0.13004532E-02

       phi(516)   =  0.58259025E+01
       theta(516) =  0.23578997E+01
       angwei(516)=  0.14578798E-02

       phi(517)   =  0.98387593E+00
       theta(517) =  0.14637889E+01
       angwei(517)=  0.14288406E-02

       phi(518)   =  0.27289562E+01
       theta(518) =  0.21759973E+01
       angwei(518)=  0.14277661E-02

       phi(519)   =  0.10515878E+01
       theta(519) =  0.74461287E+00
       angwei(519)=  0.15244258E-02

       phi(520)   =  0.13816199E+01
       theta(520) =  0.84646136E+00
       angwei(520)=  0.14627384E-02

       phi(521)   =  0.48790678E+00
       theta(521) =  0.19165231E+01
       angwei(521)=  0.14368552E-02

       phi(522)   =  0.17780305E+01
       theta(522) =  0.16690907E+01
       angwei(522)=  0.14165124E-02

       phi(523)   =  0.26777208E+01
       theta(523) =  0.53406414E-01
       angwei(523)=  0.14479023E-02

       phi(524)   =  0.39379692E+01
       theta(524) =  0.17083206E+01
       angwei(524)=  0.14520395E-02

       phi(525)   =  0.10372455E+01
       theta(525) =  0.25419853E+01
       angwei(525)=  0.14255500E-02

       phi(526)   =  0.32035980E+01
       theta(526) =  0.22130294E+01
       angwei(526)=  0.14122707E-02

       phi(527)   =  0.46632380E+01
       theta(527) =  0.14994409E+01
       angwei(527)=  0.14172882E-02

       phi(528)   =  0.53386173E+01
       theta(528) =  0.18940302E+01
       angwei(528)=  0.14356187E-02

       phi(529)   =  0.46681166E+01
       theta(529) =  0.62221366E+00
       angwei(529)=  0.14251893E-02

       phi(530)   =  0.56262293E+01
       theta(530) =  0.10568452E+01
       angwei(530)=  0.14448031E-02

       phi(531)   =  0.25021160E+01
       theta(531) =  0.18171146E+01
       angwei(531)=  0.15294061E-02

       phi(532)   =  0.58616943E+01
       theta(532) =  0.15541148E+01
       angwei(532)=  0.14559596E-02

       phi(533)   =  0.44879560E+01
       theta(533) =  0.89014512E+00
       angwei(533)=  0.14329241E-02

       phi(534)   =  0.19521799E+01
       theta(534) =  0.18370736E+01
       angwei(534)=  0.14333200E-02

       phi(535)   =  0.39587500E+01
       theta(535) =  0.20813687E+01
       angwei(535)=  0.14307173E-02

       phi(536)   =  0.22560933E+01
       theta(536) =  0.23610024E+01
       angwei(536)=  0.14389768E-02

       phi(537)   =  0.54596982E+01
       theta(537) =  0.10211383E+01
       angwei(537)=  0.14173810E-02

       phi(538)   =  0.33047929E+01
       theta(538) =  0.13242222E+01
       angwei(538)=  0.14280634E-02

       phi(539)   =  0.30424151E+01
       theta(539) =  0.21653633E+01
       angwei(539)=  0.14358205E-02

       phi(540)   =  0.21988007E+00
       theta(540) =  0.19235036E+01
       angwei(540)=  0.14207974E-02

       phi(541)   =  0.38079422E+01
       theta(541) =  0.81697512E+00
       angwei(541)=  0.14362298E-02

       phi(542)   =  0.12361592E+01
       theta(542) =  0.14358110E+01
       angwei(542)=  0.14212852E-02

       phi(543)   =  0.61937480E+01
       theta(543) =  0.33464894E+00
       angwei(543)=  0.14175157E-02

       phi(544)   =  0.25643942E+01
       theta(544) =  0.12286274E+01
       angwei(544)=  0.14105676E-02

       phi(545)   =  0.48646770E+01
       theta(545) =  0.86370116E+00
       angwei(545)=  0.14189477E-02

       phi(546)   =  0.34976404E+01
       theta(546) =  0.46262968E+00
       angwei(546)=  0.14613882E-02

       phi(547)   =  0.14554665E+01
       theta(547) =  0.20780504E+01
       angwei(547)=  0.14355081E-02

       phi(548)   =  0.44641414E+01
       theta(548) =  0.16478494E+01
       angwei(548)=  0.14271125E-02

       phi(549)   =  0.71247691E+00
       theta(549) =  0.10420808E+01
       angwei(549)=  0.13688814E-02

       phi(550)   =  0.62293768E+01
       theta(550) =  0.19075186E+01
       angwei(550)=  0.14181738E-02

       phi(551)   =  0.51457539E+01
       theta(551) =  0.27004724E+01
       angwei(551)=  0.14046977E-02

       phi(552)   =  0.86983621E+00
       theta(552) =  0.16849056E+01
       angwei(552)=  0.14476425E-02

       phi(553)   =  0.29007015E+01
       theta(553) =  0.19651550E+01
       angwei(553)=  0.14225882E-02

       phi(554)   =  0.50375776E+01
       theta(554) =  0.19767528E+01
       angwei(554)=  0.14215950E-02

       phi(555)   =  0.49143839E+01
       theta(555) =  0.61141342E+00
       angwei(555)=  0.14340070E-02

       phi(556)   =  0.23192015E+00
       theta(556) =  0.16388171E+01
       angwei(556)=  0.14224106E-02

       phi(557)   =  0.21605594E+01
       theta(557) =  0.12236048E+01
       angwei(557)=  0.14316599E-02

       phi(558)   =  0.36406956E+01
       theta(558) =  0.24305406E+01
       angwei(558)=  0.14518516E-02

       phi(559)   =  0.30401165E+01
       theta(559) =  0.20250418E+01
       angwei(559)=  0.14273189E-02

       phi(560)   =  0.19237430E+01
       theta(560) =  0.13045634E+01
       angwei(560)=  0.14478368E-02

       phi(561)   =  0.40439363E+01
       theta(561) =  0.27378726E+01
       angwei(561)=  0.14545859E-02

       phi(562)   =  0.33443382E+01
       theta(562) =  0.21206708E+01
       angwei(562)=  0.14445275E-02

       phi(563)   =  0.36918321E+01
       theta(563) =  0.18051356E+01
       angwei(563)=  0.14283616E-02

       phi(564)   =  0.54587259E+01
       theta(564) =  0.16688569E+01
       angwei(564)=  0.15356718E-02

       phi(565)   =  0.32702941E+00
       theta(565) =  0.12757547E+01
       angwei(565)=  0.14237009E-02

       phi(566)   =  0.18475275E+01
       theta(566) =  0.23174572E+01
       angwei(566)=  0.14376913E-02

       phi(567)   =  0.39678366E+01
       theta(567) =  0.10330670E+01
       angwei(567)=  0.14245578E-02

       phi(568)   =  0.48284864E+01
       theta(568) =  0.21211987E+01
       angwei(568)=  0.14307086E-02

       phi(569)   =  0.77047318E+00
       theta(569) =  0.20360379E+01
       angwei(569)=  0.14355909E-02

       phi(570)   =  0.57446833E+01
       theta(570) =  0.16493822E+01
       angwei(570)=  0.14066685E-02

       phi(571)   =  0.33343158E+01
       theta(571) =  0.14633075E+01
       angwei(571)=  0.14278633E-02

       phi(572)   =  0.33686113E+00
       theta(572) =  0.24428291E+01
       angwei(572)=  0.14393009E-02

       phi(573)   =  0.58474178E+01
       theta(573) =  0.57858521E+00
       angwei(573)=  0.14141405E-02

       phi(574)   =  0.12404848E+01
       theta(574) =  0.47349858E+00
       angwei(574)=  0.14059953E-02

       phi(575)   =  0.19772161E+01
       theta(575) =  0.14459172E+01
       angwei(575)=  0.14231557E-02

       phi(576)   =  0.38458059E+01
       theta(576) =  0.14156382E+01
       angwei(576)=  0.15396695E-02

       phi(577)   =  0.11942124E+01
       theta(577) =  0.26579921E+01
       angwei(577)=  0.14326890E-02

       phi(578)   =  0.35797517E+01
       theta(578) =  0.14026684E+01
       angwei(578)=  0.14307816E-02

       phi(579)   =  0.85077071E+00
       theta(579) =  0.14022880E+01
       angwei(579)=  0.14318903E-02

       phi(580)   =  0.42824626E+01
       theta(580) =  0.23692474E+01
       angwei(580)=  0.14128550E-02

       phi(581)   =  0.18479732E+01
       theta(581) =  0.19461192E+01
       angwei(581)=  0.14409972E-02

       phi(582)   =  0.58024325E+01
       theta(582) =  0.27299659E+01
       angwei(582)=  0.14140148E-02

       phi(583)   =  0.32770228E+01
       theta(583) =  0.26156900E+01
       angwei(583)=  0.13172831E-02

       phi(584)   =  0.28503489E+01
       theta(584) =  0.12121782E+01
       angwei(584)=  0.14139842E-02

       phi(585)   =  0.32266045E+01
       theta(585) =  0.23515017E+01
       angwei(585)=  0.14070267E-02

       phi(586)   =  0.14824526E+01
       theta(586) =  0.26218703E+01
       angwei(586)=  0.14179181E-02

       phi(587)   =  0.49397922E+01
       theta(587) =  0.12324299E+01
       angwei(587)=  0.14273359E-02

       phi(588)   =  0.59520144E+01
       theta(588) =  0.13190655E+01
       angwei(588)=  0.14090156E-02

       phi(589)   =  0.17635058E+01
       theta(589) =  0.26933696E+01
       angwei(589)=  0.14274383E-02

       phi(590)   =  0.51206636E+01
       theta(590) =  0.10023299E+01
       angwei(590)=  0.14362539E-02

       phi(591)   =  0.32758999E+01
       theta(591) =  0.11884954E+01
       angwei(591)=  0.14482429E-02

       phi(592)   =  0.41308904E+01
       theta(592) =  0.15011920E+01
       angwei(592)=  0.14140557E-02

       phi(593)   =  0.34310672E+01
       theta(593) =  0.24082410E+01
       angwei(593)=  0.14783198E-02

       phi(594)   =  0.48820591E+01
       theta(594) =  0.16102546E+01
       angwei(594)=  0.14244320E-02

       phi(595)   =  0.28053842E+01
       theta(595) =  0.44913766E+00
       angwei(595)=  0.14404260E-02

       phi(596)   =  0.59944577E+01
       theta(596) =  0.15967981E+01
       angwei(596)=  0.13944731E-02

       phi(597)   =  0.59762483E+01
       theta(597) =  0.45215464E+00
       angwei(597)=  0.14245824E-02

       phi(598)   =  0.21559427E+01
       theta(598) =  0.16114011E+01
       angwei(598)=  0.14336521E-02

       phi(599)   =  0.62503181E+01
       theta(599) =  0.17704509E+01
       angwei(599)=  0.13886028E-02

       phi(600)   =  0.17003530E+01
       theta(600) =  0.13975332E+01
       angwei(600)=  0.13842725E-02

       phi(601)   =  0.60340405E+01
       theta(601) =  0.20867679E+01
       angwei(601)=  0.14219985E-02

       phi(602)   =  0.47723708E+01
       theta(602) =  0.99008787E+00
       angwei(602)=  0.14366416E-02

       phi(603)   =  0.43518491E+01
       theta(603) =  0.25031712E+01
       angwei(603)=  0.14264456E-02

       phi(604)   =  0.15063552E+01
       theta(604) =  0.24839454E+01
       angwei(604)=  0.14257591E-02

       phi(605)   =  0.49931169E+01
       theta(605) =  0.74005717E+00
       angwei(605)=  0.14432019E-02

       phi(606)   =  0.39777308E+01
       theta(606) =  0.89269328E+00
       angwei(606)=  0.14259585E-02

       phi(607)   =  0.51624942E+01
       theta(607) =  0.20730674E+01
       angwei(607)=  0.14433162E-02

       phi(608)   =  0.51315856E+00
       theta(608) =  0.23586497E+01
       angwei(608)=  0.14157140E-02

       phi(609)   =  0.54166117E+01
       theta(609) =  0.12692925E+01
       angwei(609)=  0.14118016E-02

       phi(610)   =  0.54700246E+01
       theta(610) =  0.19616468E+01
       angwei(610)=  0.14237922E-02

       phi(611)   =  0.13205495E+01
       theta(611) =  0.19990258E+01
       angwei(611)=  0.14083047E-02

       phi(612)   =  0.71807390E+00
       theta(612) =  0.13361496E+01
       angwei(612)=  0.14347347E-02

       phi(613)   =  0.58472004E+01
       theta(613) =  0.18345515E+01
       angwei(613)=  0.14270891E-02

       phi(614)   =  0.35029843E+01
       theta(614) =  0.25542691E+01
       angwei(614)=  0.15660587E-02

       phi(615)   =  0.36477268E+01
       theta(615) =  0.16755241E+01
       angwei(615)=  0.14400824E-02

       phi(616)   =  0.53606768E+01
       theta(616) =  0.17586391E+01
       angwei(616)=  0.12670832E-02

       phi(617)   =  0.14940208E+01
       theta(617) =  0.15718086E+01
       angwei(617)=  0.14239190E-02

       phi(618)   =  0.41486855E+01
       theta(618) =  0.18641332E+01
       angwei(618)=  0.14271220E-02

       phi(619)   =  0.33167303E+01
       theta(619) =  0.19808934E+01
       angwei(619)=  0.14304267E-02

       phi(620)   =  0.35169251E+01
       theta(620) =  0.98131001E+00
       angwei(620)=  0.14392302E-02

       phi(621)   =  0.37351673E+01
       theta(621) =  0.23154078E+01
       angwei(621)=  0.14211626E-02

       phi(622)   =  0.12657359E+01
       theta(622) =  0.60877508E+00
       angwei(622)=  0.13809147E-02

       phi(623)   =  0.46973691E+01
       theta(623) =  0.11177912E+01
       angwei(623)=  0.14342704E-02

       phi(624)   =  0.57727319E+00
       theta(624) =  0.11233628E+01
       angwei(624)=  0.14213069E-02

       phi(625)   =  0.56056471E+01
       theta(625) =  0.20226030E+01
       angwei(625)=  0.14152161E-02

       phi(626)   =  0.15234364E+01
       theta(626) =  0.23446119E+01
       angwei(626)=  0.13113004E-02

       phi(627)   =  0.54186571E+00
       theta(627) =  0.24998009E+01
       angwei(627)=  0.14276823E-02

       phi(628)   =  0.55550632E+01
       theta(628) =  0.13103620E+01
       angwei(628)=  0.14375577E-02

       phi(629)   =  0.25017431E+01
       theta(629) =  0.14406101E+01
       angwei(629)=  0.14130498E-02

       phi(630)   =  0.26207154E+01
       theta(630) =  0.24397230E+01
       angwei(630)=  0.14308761E-02

       phi(631)   =  0.44181857E+01
       theta(631) =  0.14014628E+01
       angwei(631)=  0.15441464E-02

       phi(632)   =  0.59782114E+01
       theta(632) =  0.18697641E+01
       angwei(632)=  0.13173234E-02

       phi(633)   =  0.68278924E-01
       theta(633) =  0.12939001E+01
       angwei(633)=  0.14314738E-02

       phi(634)   =  0.24809768E+01
       theta(634) =  0.21162730E+00
       angwei(634)=  0.14104458E-02

       phi(635)   =  0.14723073E+01
       theta(635) =  0.12891775E+01
       angwei(635)=  0.12928731E-02

       phi(636)   =  0.50252228E+01
       theta(636) =  0.11152328E+01
       angwei(636)=  0.14230025E-02

       phi(637)   =  0.45255194E+01
       theta(637) =  0.15122312E+01
       angwei(637)=  0.14150772E-02

       phi(638)   =  0.20969374E+01
       theta(638) =  0.18590479E+01
       angwei(638)=  0.14149612E-02

       phi(639)   =  0.38257120E+01
       theta(639) =  0.11053438E+01
       angwei(639)=  0.14507603E-02

       phi(640)   =  0.14509438E+01
       theta(640) =  0.22254913E+01
       angwei(640)=  0.15521252E-02

       phi(641)   =  0.66881722E+00
       theta(641) =  0.71433651E+00
       angwei(641)=  0.14188326E-02

       phi(642)   =  0.61173124E+01
       theta(642) =  0.18100532E+01
       angwei(642)=  0.14619096E-02

       phi(643)   =  0.58826222E+01
       theta(643) =  0.10442461E+01
       angwei(643)=  0.12874934E-02

       phi(644)   =  0.99225485E+00
       theta(644) =  0.16038045E+01
       angwei(644)=  0.14300061E-02

       phi(645)   =  0.53033743E+01
       theta(645) =  0.10097818E+01
       angwei(645)=  0.14424101E-02

       phi(646)   =  0.22413275E+01
       theta(646) =  0.21270871E+01
       angwei(646)=  0.14120250E-02

       phi(647)   =  0.58811936E+01
       theta(647) =  0.83094931E+00
       angwei(647)=  0.14519484E-02

       phi(648)   =  0.52829838E+01
       theta(648) =  0.23004861E+01
       angwei(648)=  0.14195412E-02

       phi(649)   =  0.33483741E+01
       theta(649) =  0.94062102E+00
       angwei(649)=  0.14404807E-02

       phi(650)   =  0.78018945E+00
       theta(650) =  0.25352333E+01
       angwei(650)=  0.14269268E-02

       phi(651)   =  0.35350287E+00
       theta(651) =  0.15635340E+01
       angwei(651)=  0.14377113E-02

       phi(652)   =  0.60075622E+01
       theta(652) =  0.26322780E+01
       angwei(652)=  0.14286712E-02

       phi(653)   =  0.29399753E+00
       theta(653) =  0.82083327E+00
       angwei(653)=  0.14125253E-02

       phi(654)   =  0.10313065E+01
       theta(654) =  0.22897425E+01
       angwei(654)=  0.14158565E-02

       phi(655)   =  0.96903217E+00
       theta(655) =  0.11795936E+01
       angwei(655)=  0.14199740E-02

       phi(656)   =  0.59999876E+01
       theta(656) =  0.17337426E+01
       angwei(656)=  0.13898304E-02

       phi(657)   =  0.34233754E+01
       theta(657) =  0.12224236E+01
       angwei(657)=  0.14429486E-02

       phi(658)   =  0.15828084E+01
       theta(658) =  0.19944525E+01
       angwei(658)=  0.14222232E-02

       phi(659)   =  0.47437100E+01
       theta(659) =  0.16226993E+01
       angwei(659)=  0.14251589E-02

       phi(660)   =  0.39526696E+01
       theta(660) =  0.13215281E+01
       angwei(660)=  0.14255054E-02

       phi(661)   =  0.48642397E+01
       theta(661) =  0.13570019E+01
       angwei(661)=  0.14316170E-02

       phi(662)   =  0.21152999E+01
       theta(662) =  0.69418180E+00
       angwei(662)=  0.14293096E-02

       phi(663)   =  0.36739533E+01
       theta(663) =  0.21792083E+01
       angwei(663)=  0.14502127E-02

       phi(664)   =  0.38444841E+01
       theta(664) =  0.18266829E+01
       angwei(664)=  0.14161218E-02

       phi(665)   =  0.51920834E+01
       theta(665) =  0.19360954E+01
       angwei(665)=  0.14430539E-02

       phi(666)   =  0.50155826E+01
       theta(666) =  0.13505858E+01
       angwei(666)=  0.14339829E-02

       phi(667)   =  0.23526587E+01
       theta(667) =  0.14069036E+01
       angwei(667)=  0.14366622E-02

       phi(668)   =  0.54959559E+01
       theta(668) =  0.28023660E+01
       angwei(668)=  0.14203831E-02

       phi(669)   =  0.12806358E+01
       theta(669) =  0.25184035E+01
       angwei(669)=  0.14473184E-02

       phi(670)   =  0.45156908E+01
       theta(670) =  0.51196778E+00
       angwei(670)=  0.14466171E-02

       phi(671)   =  0.48303475E+01
       theta(671) =  0.17442203E+01
       angwei(671)=  0.14330748E-02

       phi(672)   =  0.57512469E+01
       theta(672) =  0.19383826E+01
       angwei(672)=  0.14455302E-02

       phi(673)   =  0.60807877E+01
       theta(673) =  0.19583306E+01
       angwei(673)=  0.14809545E-02

       phi(674)   =  0.48716891E+00
       theta(674) =  0.17735687E+01
       angwei(674)=  0.14287338E-02

       phi(675)   =  0.23449488E+01
       theta(675) =  0.17619830E+01
       angwei(675)=  0.14138463E-02

       phi(676)   =  0.86025435E+00
       theta(676) =  0.15434074E+01
       angwei(676)=  0.14468083E-02

       phi(677)   =  0.37879457E+01
       theta(677) =  0.20615716E+01
       angwei(677)=  0.14191865E-02

       phi(678)   =  0.16637945E+01
       theta(678) =  0.22640245E+01
       angwei(678)=  0.14161667E-02

       phi(679)   =  0.58573580E+00
       theta(679) =  0.12684261E+01
       angwei(679)=  0.14193770E-02

       phi(680)   =  0.12492483E+01
       theta(680) =  0.15801109E+01
       angwei(680)=  0.14739583E-02

       phi(681)   =  0.22210889E+01
       theta(681) =  0.97638190E+00
       angwei(681)=  0.14361239E-02

       phi(682)   =  0.17776165E+01
       theta(682) =  0.10491619E+01
       angwei(682)=  0.14526796E-02

       phi(683)   =  0.20040908E+01
       theta(683) =  0.19722826E+01
       angwei(683)=  0.14289016E-02

       phi(684)   =  0.27079394E+01
       theta(684) =  0.16157929E+01
       angwei(684)=  0.14446330E-02

       phi(685)   =  0.25109684E+01
       theta(685) =  0.50294071E+00
       angwei(685)=  0.14373082E-02

       phi(686)   =  0.22638640E+01
       theta(686) =  0.11206176E+01
       angwei(686)=  0.14284834E-02

       phi(687)   =  0.53211126E+01
       theta(687) =  0.11127435E+00
       angwei(687)=  0.13972082E-02

       phi(688)   =  0.30201983E+01
       theta(688) =  0.25796287E+01
       angwei(688)=  0.14472946E-02

       phi(689)   =  0.48593755E+01
       theta(689) =  0.11116065E+01
       angwei(689)=  0.14256048E-02

       phi(690)   =  0.43566246E+01
       theta(690) =  0.78907841E+00
       angwei(690)=  0.14377240E-02

       phi(691)   =  0.17090061E+00
       theta(691) =  0.23573208E+01
       angwei(691)=  0.14329718E-02

       phi(692)   =  0.30862660E+01
       theta(692) =  0.15109882E+01
       angwei(692)=  0.14339414E-02

       phi(693)   =  0.62093682E+01
       theta(693) =  0.12470018E+01
       angwei(693)=  0.14128769E-02

       phi(694)   =  0.49973130E+01
       theta(694) =  0.21076887E+01
       angwei(694)=  0.14237501E-02

       phi(695)   =  0.45419474E+01
       theta(695) =  0.17667595E+01
       angwei(695)=  0.14167115E-02

       phi(696)   =  0.12489823E+01
       theta(696) =  0.73740309E+00
       angwei(696)=  0.13408141E-02

       phi(697)   =  0.13417096E+01
       theta(697) =  0.23743668E+01
       angwei(697)=  0.14165894E-02

       phi(698)   =  0.28591774E+01
       theta(698) =  0.81415552E+00
       angwei(698)=  0.14072180E-02

       phi(699)   =  0.34754613E+00
       theta(699) =  0.21415970E+01
       angwei(699)=  0.14281906E-02

       phi(700)   =  0.23850479E+01
       theta(700) =  0.24759812E+01
       angwei(700)=  0.14265561E-02
c
       return
       end
