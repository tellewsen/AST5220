
  �  =   k820309    P          16.0        �=W                                                                                                           
       ode_solver.f90 ODE_SOLVER                                                     
                                                           
                                                           
                                                                                                                                                                                                                                                                                                                         #         @                                                    	   #Y    #DYDX 	   #X 
   #HTRY    #EPS    #YSCAL    #HDID    #HNEXT    #DERIVS              
                                                  
               &                                                     
                                 	                   
              &                                                     
                                
     
                 
                                      
                
                                      
                
                                                    
              &                                                                                          
                                                      
       #         @                                        	               #X    #Y    #DYDX              
                                     
                
                                                   
 	             &                                                                                                       
 
              &                                                      @                                                       @                                                       @                                                                                                               @                                   
                 @P                                                 
                &                                                     @ P                                                 
                &                   &                                           #         @                                                    	   #YSTART    #X1    #X2    #EPS    #H1     #HMIN !   #DERIVS "   #RKQS &   #OUTPUT 3          0  
D@                                                 
               &                                                     
                                      
                
                                      
                
  @                                   
                
  @                                    
                
                                 !     
      #         @    @                             "     	               #X #   #Y $   #DYDX %                                   
                                #     
                
                                $                   
              &                                                                                    %                   
               &                                           #         @                                  &     	            	   #Y '   #DYDX (   #X )   #HTRY *   #EPS +   #YSCAL ,   #HDID -   #HNEXT .   #DERIVS /                                 
                               '                   
               &                                                     
                                (                   
 	             &                                                     
                               )     
                 
                                *     
                
                                +     
                
                                ,                   
 
             &                                                                                    -     
                                                .     
       #         @                                   /     	               #X 0   #Y 1   #DYDX 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
                                0     
                
                                1                   
              &                                                                                    2                   
               &                                           #         @                                   3     	               #X 4   #Y 5                                   
                                4     
                
                                5                   
              &                                           (         D                              6                                   
    #P 7   #N 8             &                                                    DP                              7                   
               &                                                     
  @                              8           (         D                              9                                   
    #P :   #N ;   #M <             &                   &                                                    DP                              :                   
               &                   &                                                     
  @                              ;                     
  @                              <              �   "      fn#fn    �   @   J   HEALPIX_TYPES      @   J   RK_MOD    B  @   J   BS_MOD "   �  p       I4B+HEALPIX_TYPES "   �  p       LGT+HEALPIX_TYPES !   b  p       DP+HEALPIX_TYPES    �  �       RKQS+RK_MOD    q  �   a   RKQS%Y+RK_MOD !   �  �   a   RKQS%DYDX+RK_MOD    �  @   a   RKQS%X+RK_MOD !   �  @   a   RKQS%HTRY+RK_MOD     	  @   a   RKQS%EPS+RK_MOD "   I  �   a   RKQS%YSCAL+RK_MOD !   �  @   a   RKQS%HDID+RK_MOD "     @   a   RKQS%HNEXT+RK_MOD #   U  `      RKQS%DERIVS+RK_MOD %   �  @   a   RKQS%DERIVS%X+RK_MOD %   �  �   a   RKQS%DERIVS%Y+RK_MOD (   �  �   a   RKQS%DERIVS%DYDX+RK_MOD      @       NOK    M  @       NBAD    �  @       KOUNT    �  @       SAVE_STEPS    	  @       DXSAV    M	  �       XP    �	  �       YP    }
  �       ODEINT      �   a   ODEINT%YSTART    �  @   a   ODEINT%X1    �  @   a   ODEINT%X2    *  @   a   ODEINT%EPS    j  @   a   ODEINT%H1    �  @   a   ODEINT%HMIN    �  v      ODEINT%DERIVS     `  @   a   ODEINT%DERIVS%X     �  �   a   ODEINT%DERIVS%Y #   ,  �   a   ODEINT%DERIVS%DYDX    �  �      ODEINT%RKQS    k  �   a   ODEINT%RKQS%Y !   �  �   a   ODEINT%RKQS%DYDX    �  @   a   ODEINT%RKQS%X !   �  @   a   ODEINT%RKQS%HTRY       @   a   ODEINT%RKQS%EPS "   C  �   a   ODEINT%RKQS%YSCAL !   �  @   a   ODEINT%RKQS%HDID "     @   a   ODEINT%RKQS%HNEXT #   O  K     ODEINT%RKQS%DERIVS %   �  @   a   ODEINT%RKQS%DERIVS%X %   �  �   a   ODEINT%RKQS%DERIVS%Y (   f  �   a   ODEINT%RKQS%DERIVS%DYDX    �  l      ODEINT%OUTPUT     ^  @   a   ODEINT%OUTPUT%X     �  �   a   ODEINT%OUTPUT%Y    *  �       REALLOCATE1    �  �   a   REALLOCATE1%P    `  @   a   REALLOCATE1%N    �  �       REALLOCATE2    i  �   a   REALLOCATE2%P      @   a   REALLOCATE2%N    M  @   a   REALLOCATE2%M 