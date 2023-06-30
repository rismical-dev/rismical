C#NUMPAC#HERD31              REVISED ON 1984-11-30                      
C      SUBROUTINE HERM31(I,X,Y,M,N,XI,YI,YD,ND,ILL)                     
      SUBROUTINE HERD31(I,X,Y,M,N,XI,YI,YD,ND,ILL)                      
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION XI(N),YI(N),YD(ND),Z(2,2)                               
      IF(N.LE.1.OR.M.LT.0.OR.M.GE.3) GO TO 1010                         
      IF(ILL.EQ.0) GO TO 1000                                           
      MD=ILL                                                            
      IF(N.LE.MD) GO TO 1010                                            
C      CALL DERIV1(XI,YI,YD,N,1,1,MD)                             
       CALL DERID1(XI,YI,YD,N,1,1,MD)                                
 1000 CONTINUE                                                          
      Z(1,1)=YI(I)                                                      
      Z(1,2)=YD(I)                                                      
      Z(2,1)=YI(I+1)                                                    
      Z(2,2)=YD(I+1)                                                    
C      CALL PHER31(X,Y,M,XI(I),Z,ILL)                        
      CALL PHED31(X,Y,M,XI(I),Z,ILL)                                    
      RETURN                                                            
 1010 ILL=30000                                                         
      RETURN                                                            
      END                                                               
C#NUMPAC#DERIV1              REVISED ON 1984-11-30                      
C      SUBROUTINE DERIV1(X,Y,F,  N,M,L,MD)                      
      SUBROUTINE DERID1(X,Y,F,  N,M,L,MD)                               
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),Y(M),F(L)                                          
      N1=N-1                                                            
      IF(MD.GE.2) GO TO 1010                                            
      IE=N1*L+1                                                         
      I=N1*M+1                                                          
      IM=I-M                                                            
      F(IE)=(Y(I)-Y(IM))/(X(N)-X(N1))                                   
      DO 1000 I=1,N1                                                    
      IM=I-1                                                            
      II=IM*M+1                                                         
      IE=IM*L+1                                                         
      IP=II+M                                                           
 1000 F(IE)=(Y(IP)-Y(II))/(X(I+1)-X(I))                                 
      RETURN                                                            
 1010 CONTINUE                                                          
      H1=X(2)-X(1)                                                      
      H2=X(3)-X(2)                                                      
      H1P2=H1+H2                                                        
      T1=(Y(M+1)-Y(1))/H1                                               
      T2=(Y(2*M+1)-Y(M+1))/H2                                           
      F(1)=T1+H1/H1P2*(T1-T2)                                           
      IE=N1*L+1                                                         
      I=N1*M+1                                                          
      IM=I-M                                                            
      IMM=IM-M                                                          
      H1=X(N1)-X(N1-1)                                                  
      H2=X(N)-X(N1)                                                     
      H1P2=H1+H2                                                        
      T1=(Y(IM)-Y(IMM))/H1                                              
      T2=(Y(I)-Y(IM))/H2                                                
      F(IE)=(H1P2+H2)/H1P2*(T2-T1)+T1                                   
      M1=1                                                              
      IF(MD.GE.3) M1=N-3                                                
      DO 1020 I=2,N1,M1                                                 
      IM=I-1                                                            
      II=IM*M+1                                                         
      IE=IM*L+1                                                         
      H1=X(I)-X(IM)                                                     
      H2=X(I+1)-X(I)                                                    
      H1P2=H1+H2                                                        
      IPM=II+M                                                          
      IMM=II-M                                                          
      T1=(Y(II)-Y(IMM))/H1                                              
      T2=(Y(IPM)-Y(II))/H2                                              
      F(IE)=(T2*H1+T1*H2)/H1P2                                          
 1020 CONTINUE                                                          
      IF(MD.LE.2) RETURN                                                
 1030 CONTINUE                                                          
      M3=N-2                                                            
      DO 1040 I=3,M3                                                    
      IM=I-1                                                            
      II=IM*M+1                                                         
      IE=IM*L+1                                                         
      IPM1=II+M                                                         
      IPM2=IPM1+M                                                       
      IMM1=II-M                                                         
      IMM2=IMM1-M                                                       
      Y2=X(I-2)                                                         
      Y1=X(I-1)                                                         
      X0=X(I)                                                           
      X1=X(I+1)                                                         
      X2=X(I+2)                                                         
      G2=Y1-Y2                                                          
      G1=X0-Y1                                                          
      H0=X1-X0                                                          
      H1=X2-X1                                                          
      H21=G2+G1                                                         
      H10=G1+H0                                                         
      H01=H0+H1                                                         
      H210=H21+H0                                                       
      H101=H10+H1                                                       
      H2101=H210+H1                                                     
      AM2=Y(IMM2)/(G2*H21*H210*H2101)                                   
      AM1=-Y(IMM1)/(G2*G1*H10*H101)                                     
      A0=Y(II)/(H21*G1*H0*H01)                                          
      AP1=-Y(IPM1)/(H210*H10*H0*H1)                                     
      AP2=Y(IPM2)/(H2101*H101*H01*H1)                                   
      F(IE)=(AM2*G1*H0*H01+ AM1*H21*H0*H01                              
     *      +A0*(G1*H0*H01+H21*H0*H01-H21*G1*H01-H21*G1*H0)             
     *       +AP1*(-H21*G1*H01) +AP2*(-H21*G1*H0))                      
 1040 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C#NUMPAC#PHER31              REVISED ON 1984-11-30                      
C      SUBROUTINE PHER31(X,Y,M,XI,Z,ILL)             
      SUBROUTINE PHED31(X,Y,M,XI,Z,ILL)                                 
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(M),XI(2),YI(2),Y1I(2),Z(2,2)                          
      YI(1)=Z(1,1)                                                      
      YI(2)=Z(2,1)                                                      
      Y1I(1)=Z(1,2)                                                     
      Y1I(2)=Z(2,2)                                                     
C      X           XCOORDINATE(INPUT)                                   
C      Y(M)        Y-VALUE(OUTPUT)                                      
C      M           ORDER OF DERIVATIVE(INPUT)                           
C      XI          X-RANGE(INPUT)                                       
C      YI           FUNCTION VALUE ON XI(INPUT)                         
C      Y1I         DERIVATIVE ON XI(INPUT)                              
C      ILL         OUTPUT)                                              
      ILL=0                                             
      EPS=1.D-6                                         
      X1=XI(1)                                          
      X2=XI(2)                                          
      H=X2-X1                                           
      IF(H) 1000 , 1060 , 1010                          
 1000 X1=XI(2)                                          
      X2=XI(1)                                          
      H=-H                                              
 1010 T=(X-X1)/H                                        
      IF(ABS(T).LE.EPS) T=0.0D0
      IF(ABS(T-1.D0).LE.EPS) T=1.0D0
      IF(T.LT.0.0d0.OR.T.GT.1.0d0) GO TO 1060      
      T2=T*T                                            
      T1=T-1.d0                                    
      M1=M+1                                            
      DO 1050 K=1,M1                                    
      GO TO ( 1020 , 1030 , 1040 ),K                    
 1020 P=1.d0+T2*(2.d0*T-3.d0)                                 
      Q=1.d0-P                                            
      R=T*T1**2                                         
      S=T2*T1                                           
      GO TO 1050                                        
 1030 P=6.d0*(T2-T)                                       
      Q=-P                                              
      R=T1*(3.d0*T-1.d0)                                    
      S=T*(3.d0*T-2.d0)                                     
      GO TO 1050                                        
 1040 P=-6.d0+12.d0*T                                       
      Q=-P                                              
      R=6.d0*T-4.d0                                         
      S=6.d0*T-2.d0                                         
 1050 Y(K)=(YI(1)*P+YI(2)*Q+(Y1I(1)*R+Y1I(2)*S)*H)/(H**(K-1)) 
      RETURN                                                            
 1060 ILL=30000                                                         
      RETURN                                                            
      END                                                               
