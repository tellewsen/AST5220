	  ÝB  ¥   k820309    Ö          15.0        èÖV                                                                                                           
       rec_mod.f90 REC_MOD                                                    
                                                          
                                                          
                                                         
                                                         
                                                                                                                                                                                                                                                                                                                         %         @                               	                   
       #GET_H%SQRT 
   #GET_H%EXP    #X                                               
     SQRT                                                 EXP           
                                      
                                                       
                 
                 Zd;ßO§?        0.046D0                                                 
                   
                  ¬<kM$Ë:                                                         
                 
                 NUªé`:        1.673534D-27                                                 
                 
                 ÍÌÌÌÌÌ@        2.725D0                                                 
                 
                 07»Üy²9        9.10938188D-31                                                 
                 
                 n@
Ýç°0;        1.3806503D-23                                                 
                 
                 ë,]§¦á8        1.05457148D-34                                                 
                
          )       -DTû!	@        3.141592653589793238462643383279502884197                                                 
                   
                  zÒQD<        #         @                                                  	   #ODEINT%SIZE    #ODEINT%SIGN    #ODEINT%ABS    #YSTART    #X1    #X2    #EPS    #H1    #HMIN    #DERIVS     #RKQS $   #OUTPUT 1                                                   SIZE                                                 SIGN                                                 ABS           
                                                  
               &                                                     
                                      
                
                                      
                
                                      
                
                                      
                
                                      
      #         @                                         	               #X !   #Y "   #DYDX #             
                                !     
                
                                "                   
              &                                                                                    #                   
               &                                           #         @                                   $     	            	   #Y %   #DYDX &   #X '   #HTRY (   #EPS )   #YSCAL *   #HDID +   #HNEXT ,   #DERIVS -             
                               %                   
               &                                                     
                                &                   
 	             &                                                     
                               '     
                 
                                (     
                
                                )     
                
                                *                   
 
             &                                                                                    +     
                                                ,     
       #         @                                   -     	               #X .   #Y /   #DYDX 0             
                                .     
                
                                /                   
              &                                                                                    0                   
               &                                           #         @                                   1     	               #X 2   #Y 3             
                                2     
                
                                3                   
              &                                           #         @    @                             4                	   #BSSTEP%SIZE 5   #BSSTEP%MAXVAL 6   #BSSTEP%ABS 7   #BSSTEP%MAX 8   #BSSTEP%MIN 9   #Y :   #DYDX ;   #X <   #HTRY =   #EPS >   #YSCAL ?   #HDID @   #HNEXT A   #DERIVS B                                              5     SIZE                                            6     MAXVAL                                            7     ABS                                            8     MAX                                            9     MIN           
                               :                   
 
              &                                                     
                                 ;                   
              &                                                     
                                <     
                 
                                 =     
                
                                 >     
                
                                 ?                   
              &                                                                                     @     
                                                 A     
       #         @                                   B     	               #X C   #Y D   #DYDX E             
                                C     
                
                                D                   
              &                                                                                    E                   
 	              &                                           #         @                                  F                   #SPLINE%SIZE G   #X H   #Y I   #YP1 J   #YPN K   #Y2 L                                              G     SIZE           
                                H                   
              &                                                     
                                 I                   
              &                                                     
                                 J     
                
                                 K     
                                                L                   
               &                                                                                       M     
                 
                 Év@M:        6.652462D-29                                            N     
                 
                    JxÞ±A        2.99792458D8                                            O     
                 
                 ñ&7MÔã}?        7.29735308D-3                                            P     
                 
                    `9t @        8.227%         @                               Q                   
       #SPLINT%SIZE R   #SPLINT%MAX S   #SPLINT%MIN T   #XA U   #YA V   #Y2A W   #X X                                              R     SIZE                                            S     MAX                                            T     MIN           
                                U                   
              &                                                     
                                 V                   
              &                                                     
                                 W                   
              &                                                     
                                 X     
      %         @                               Y                   
       #SPLINT_DERIV%SIZE Z   #SPLINT_DERIV%MAX [   #SPLINT_DERIV%MIN \   #XA ]   #YA ^   #Y2A _   #X `                                              Z     SIZE                                            [     MAX                                            \     MIN           
                                ]                   
              &                                                     
                                 ^                   
              &                                                     
                                 _                   
              &                                                     
                                 `     
      #         @                                   a                    #X b   #Y c             
                                 b     
                
                                 c                   
              &                                                      @                                d                     @ @                              e                   
                &                                                    @                                f                   
                &                                                    @                                g                   
                &                                                    @ @                              h                   
                &                                                    @                                i                   
                &                                                    @ @                              j                   
                &                                                    @                                k                   
                &                                                    @                                l                   
                &                                                    @                                m                   
                &                                                    @ @                              n                   
                &                                                    @                                o                   
                &                                                    @ @                              p                   
                &                                                    @                                q                   
                &                                                    @ @                              r                   
                &                                                    @ @                              s                   
                &                                                    @                                t                   
                &                                                    @                                u                   
                &                                                    @                                v                   
                &                                                    @                                w                   
                &                                                    @ @                              x                   
                &                                                      @                                y                       @                                z                     @ @                              {                   
                &                                                    @                                |                   
                &                                                    @                                }                   
                &                                                    @                                ~                   
                &                                                    @                                                   
                &                                                    @                                                   
                &                                                                                                        
                &                                                                                            
       #         @                                                       #INITIALIZE_REC_MOD%SQRT    #INITIALIZE_REC_MOD%ABS    #INITIALIZE_REC_MOD%EXP    #INITIALIZE_REC_MOD%LOG                                                   SQRT                                                ABS                                                EXP                                                LOG #         @    @                                                #DX_EDX%SQRT    #DX_EDX%EXP    #DX_EDX%LOG    #X_REC    #X_E    #DYDX                                                   SQRT                                                EXP                                                LOG           
                                      
                
                                                    
              &                                                     D                                                   
               &                                           #         @    @                                                 #X    #Y              
                                      
                
                                                    
               &                                           %         @                                                  
       #GET_N_E%EXP    #X_IN                                                   EXP           
  @                                   
      #         @    @                                                 #X_REC    #TAU    #DYDX              
                                      
                
                                                    
              &                                                     D                                                   
               &                                           %         @                                                   
       #X_IN              
  @                                   
      %         @                                                   
       #X_IN              
  @                                   
      %         @                                                   
       #X              
                                      
      %         @                                                    
       #X               
                                       
      %         @                                ¡                    
       #X ¢             
                                 ¢     
      %         @                                £                    
       #X ¤             
                                 ¤     
                   fn#fn    ¼   @   J   HEALPIX_TYPES    ü   @   J   PARAMS    <  @   J   TIME_MOD    |  @   J   ODE_SOLVER    ¼  @   J   SPLINE_1D_MOD "   ü  p       I4B+HEALPIX_TYPES !   l  p       DP+HEALPIX_TYPES "   Ü  p       LGT+HEALPIX_TYPES    L  v       GET_H+TIME_MOD $   Â  =      GET_H%SQRT+TIME_MOD #   ÿ  <      GET_H%EXP+TIME_MOD !   ;  @   a   GET_H%X+TIME_MOD    {  w       OMEGA_B+PARAMS    ò  p       RHO_C+PARAMS    b  |       M_H+PARAMS    Þ  w       T_0+PARAMS    U  ~       M_E+PARAMS    Ó  }       K_B+PARAMS    P  ~       HBAR+PARAMS !   Î         PI+HEALPIX_TYPES !   g  p       EPSILON_0+PARAMS "   ×  Ó       ODEINT+ODE_SOLVER '   ª	  =      ODEINT%SIZE+ODE_SOLVER '   ç	  =      ODEINT%SIGN+ODE_SOLVER &   $
  <      ODEINT%ABS+ODE_SOLVER )   `
     a   ODEINT%YSTART+ODE_SOLVER %   ì
  @   a   ODEINT%X1+ODE_SOLVER %   ,  @   a   ODEINT%X2+ODE_SOLVER &   l  @   a   ODEINT%EPS+ODE_SOLVER %   ¬  @   a   ODEINT%H1+ODE_SOLVER '   ì  @   a   ODEINT%HMIN+ODE_SOLVER )   ,  `      ODEINT%DERIVS+ODE_SOLVER +     @   a   ODEINT%DERIVS%X+ODE_SOLVER +   Ì     a   ODEINT%DERIVS%Y+ODE_SOLVER .   X     a   ODEINT%DERIVS%DYDX+ODE_SOLVER '   ä        ODEINT%RKQS+ODE_SOLVER )        a   ODEINT%RKQS%Y+ODE_SOLVER ,        a   ODEINT%RKQS%DYDX+ODE_SOLVER )     @   a   ODEINT%RKQS%X+ODE_SOLVER ,   Û  @   a   ODEINT%RKQS%HTRY+ODE_SOLVER +     @   a   ODEINT%RKQS%EPS+ODE_SOLVER -   [     a   ODEINT%RKQS%YSCAL+ODE_SOLVER ,   ç  @   a   ODEINT%RKQS%HDID+ODE_SOLVER -   '  @   a   ODEINT%RKQS%HNEXT+ODE_SOLVER .   g  `      ODEINT%RKQS%DERIVS+ODE_SOLVER 0   Ç  @   a   ODEINT%RKQS%DERIVS%X+ODE_SOLVER 0        a   ODEINT%RKQS%DERIVS%Y+ODE_SOLVER 3        a   ODEINT%RKQS%DERIVS%DYDX+ODE_SOLVER )     V      ODEINT%OUTPUT+ODE_SOLVER +   u  @   a   ODEINT%OUTPUT%X+ODE_SOLVER +   µ     a   ODEINT%OUTPUT%Y+ODE_SOLVER    A  ó       BSSTEP+BS_MOD #   4  =      BSSTEP%SIZE+BS_MOD %   q  ?      BSSTEP%MAXVAL+BS_MOD "   °  <      BSSTEP%ABS+BS_MOD "   ì  <      BSSTEP%MAX+BS_MOD "   (  <      BSSTEP%MIN+BS_MOD     d     a   BSSTEP%Y+BS_MOD #   ð     a   BSSTEP%DYDX+BS_MOD     |  @   a   BSSTEP%X+BS_MOD #   ¼  @   a   BSSTEP%HTRY+BS_MOD "   ü  @   a   BSSTEP%EPS+BS_MOD $   <     a   BSSTEP%YSCAL+BS_MOD #   È  @   a   BSSTEP%HDID+BS_MOD $     @   a   BSSTEP%HNEXT+BS_MOD %   H  `      BSSTEP%DERIVS+BS_MOD '   ¨  @   a   BSSTEP%DERIVS%X+BS_MOD '   è     a   BSSTEP%DERIVS%Y+BS_MOD *   t     a   BSSTEP%DERIVS%DYDX+BS_MOD %             SPLINE+SPLINE_1D_MOD *     =      SPLINE%SIZE+SPLINE_1D_MOD '   ¾     a   SPLINE%X+SPLINE_1D_MOD '   J     a   SPLINE%Y+SPLINE_1D_MOD )   Ö  @   a   SPLINE%YP1+SPLINE_1D_MOD )     @   a   SPLINE%YPN+SPLINE_1D_MOD (   V     a   SPLINE%Y2+SPLINE_1D_MOD    â  |       SIGMA_T+PARAMS    ^  |       C+PARAMS    Ú  }       ALPHA+PARAMS #   W  u       LAMBDA_2S1S+PARAMS %   Ì  ¡       SPLINT+SPLINE_1D_MOD *   m   =      SPLINT%SIZE+SPLINE_1D_MOD )   ª   <      SPLINT%MAX+SPLINE_1D_MOD )   æ   <      SPLINT%MIN+SPLINE_1D_MOD (   "!     a   SPLINT%XA+SPLINE_1D_MOD (   ®!     a   SPLINT%YA+SPLINE_1D_MOD )   :"     a   SPLINT%Y2A+SPLINE_1D_MOD '   Æ"  @   a   SPLINT%X+SPLINE_1D_MOD +   #  ³       SPLINT_DERIV+SPLINE_1D_MOD 0   ¹#  =      SPLINT_DERIV%SIZE+SPLINE_1D_MOD /   ö#  <      SPLINT_DERIV%MAX+SPLINE_1D_MOD /   2$  <      SPLINT_DERIV%MIN+SPLINE_1D_MOD .   n$     a   SPLINT_DERIV%XA+SPLINE_1D_MOD .   ú$     a   SPLINT_DERIV%YA+SPLINE_1D_MOD /   %     a   SPLINT_DERIV%Y2A+SPLINE_1D_MOD -   &  @   a   SPLINT_DERIV%X+SPLINE_1D_MOD     R&  V       OUTPUT+TIME_MOD "   ¨&  @   a   OUTPUT%X+TIME_MOD "   è&     a   OUTPUT%Y+TIME_MOD    t'  @       N    ´'         X_REC    @(         A_REC    Ì(         Z_REC    X)         TAU    ä)         DTAU    p*         DDTAU    ü*         LOGTAU    +         LOGDTAU    ,         LOGDDTAU     ,         D4TAU    ,-         LOGD4TAU    ¸-         N_E    D.         N_E2    Ð.         LOGN_E    \/         LOGN_E2    è/         G    t0         G2     1         G22    1         H_REC    2         X_E    ¤2  @       J    ä2  @       K    $3         X_TEST    °3         N_ETEST    <4         Z_TEST    È4         A_TEST    T5         TAU_TEST    à5         DTAU_TEST    l6         DDTAU_TEST    ø6  @       X_0 #   87  ¹       INITIALIZE_REC_MOD (   ñ7  =      INITIALIZE_REC_MOD%SQRT '   .8  <      INITIALIZE_REC_MOD%ABS '   j8  <      INITIALIZE_REC_MOD%EXP '   ¦8  <      INITIALIZE_REC_MOD%LOG    â8         DX_EDX    y9  =      DX_EDX%SQRT    ¶9  <      DX_EDX%EXP    ò9  <      DX_EDX%LOG    .:  @   a   DX_EDX%X_REC    n:     a   DX_EDX%X_E    ú:     a   DX_EDX%DYDX    ;  V       OUTPUT1    Ü;  @   a   OUTPUT1%X    <     a   OUTPUT1%Y    ¨<  k       GET_N_E    =  <      GET_N_E%EXP    O=  @   a   GET_N_E%X_IN    =  f       DTAUDX    õ=  @   a   DTAUDX%X_REC    5>     a   DTAUDX%TAU    Á>     a   DTAUDX%DYDX    M?  Z       GET_TAU    §?  @   a   GET_TAU%X_IN    ç?  Z       GET_DTAU    A@  @   a   GET_DTAU%X_IN    @  W       GET_DDTAU    Ø@  @   a   GET_DDTAU%X    A  W       GET_G    oA  @   a   GET_G%X    ¯A  W       GET_DG    B  @   a   GET_DG%X    FB  W       GET_DDG    B  @   a   GET_DDG%X 