  /A  6   k820309    ĥ
          2021.4.0    wŜa                                                                                                          
       module_lightning_nox_driver.f90 MODULE_LIGHTNING_NOX_DRIVER #         @                                                   3   #LIGHTNING_NOX_DRIVER%NUM_MOIST    #CURR_SECS    #DT    #DX    #DY    #XLAT    #XLON    #XLAND    #HT    #T_PHY    #P_PHY    #RHO    #U    #V    #W    #Z    #MOIST    #IC_FLASHRATE    #CG_FLASHRATE    #REFL    #LIGHTNING_OPTION    #LIGHTNING_DT    #LIGHTNING_START_SECONDS    #N_IC    #N_CG     #LNOX_OPT !   #LNOX_PASSIVE "   #LTNG_TEMP_UPPER #   #LTNG_TEMP_LOWER $   #CELLCOUNT_METHOD %   #IDS &   #IDE '   #JDS (   #JDE )   #KDS *   #KDE +   #IMS 
   #IME 	   #JMS    #JME    #KMS    #KME    #ITS ,   #ITE -   #JTS .   #JTE /   #KTS 0   #KTE 1   #C_NO 2   #LNOX_TOTAL 3   #LNOX_IC 4   #LNOX_CG 5                                                                                                                                                      
                                      
                
  @                                    	                
  @                                    	                
  @                                    	               
  @                                                   	      5  p &       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p '       r    5  p &       r    p                                   
                                                      	      5  p &       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p '       r    5  p &       r    p                                   
  @                                                   	      5  p &       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p '       r    5  p &       r    p                                   
  @                                                   	      5  p &       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p '       r    5  p &       r    p                                   
  @                                                   	        5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                   
  @                                                   	        5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                   
  @                                                   	        5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                   
                                                      	        5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                   
                                                      	 	       5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                   
                                                      	 
       5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                   
  @                                                   	        5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                   
                                                      	          p          5  p '       r    5  p &       r    p        5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r      5 r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p          5 r                               
  @                                                   	      5  p &       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p '       r    5  p &       r    p                                   
  @                                                   	      5  p &       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p '       r    5  p &       r    p                                   
  @                                                   	        5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                    
                                                       
                                       	                
  @                                    	                
  @                                    	                
  @                                     	                
                                  !                     
                                  "                     
  @                               #     	                
  @                               $     	                
  @                               %                     
  @                               &                     
  @                               '                     
  @                               (                     
  @                               )                     
  @                               *                     
  @                               +                     
  @                     @         
                     
  @                     @         	                     
  @                     @                              
  @                     @                              
  @                     @                              
  @                     @                              
  @                               ,                     
  @                               -                     
  @                               .                     
  @                               /                     
  @                               0                     
  @                               1                    
D                                 2                    	         5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                   
F                                 3                    	         5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                   
F                                 4                    	         5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                   
F                                 5                    	         5  p &       r      5  p )       r    5  p (       r    p        5  p (       r      5  p %       r 	   5  p $       r 
   p        5  p $       r 
     & 5  p $       r 
   5  p %       r 	     & 5  p (       r    5  p )       r      & 5  p &       r    5  p '       r          5  p %       r 	   5  p $       r 
   p            5  p )       r    5  p (       r    p            5  p '       r    5  p &       r    p                                 D      fn#fn %   ä         LIGHTNING_NOX_DRIVER H   ó  @     LIGHTNING_NOX_DRIVER%NUM_MOIST+MODULE_STATE_DESCRIPTION /   3  @   a   LIGHTNING_NOX_DRIVER%CURR_SECS (   s  @   a   LIGHTNING_NOX_DRIVER%DT (   ³  @   a   LIGHTNING_NOX_DRIVER%DX (   ó  @   a   LIGHTNING_NOX_DRIVER%DY *   3    a   LIGHTNING_NOX_DRIVER%XLAT *   G    a   LIGHTNING_NOX_DRIVER%XLON +   [	    a   LIGHTNING_NOX_DRIVER%XLAND (   o    a   LIGHTNING_NOX_DRIVER%HT +       a   LIGHTNING_NOX_DRIVER%T_PHY +       a   LIGHTNING_NOX_DRIVER%P_PHY )   Ğ    a   LIGHTNING_NOX_DRIVER%RHO '   ż    a   LIGHTNING_NOX_DRIVER%U '   Ó    a   LIGHTNING_NOX_DRIVER%V '   ç    a   LIGHTNING_NOX_DRIVER%W '   û    a   LIGHTNING_NOX_DRIVER%Z +   #    a   LIGHTNING_NOX_DRIVER%MOIST 2   £&    a   LIGHTNING_NOX_DRIVER%IC_FLASHRATE 2   ·(    a   LIGHTNING_NOX_DRIVER%CG_FLASHRATE *   Ë*    a   LIGHTNING_NOX_DRIVER%REFL 6   ß-  @   a   LIGHTNING_NOX_DRIVER%LIGHTNING_OPTION 2   .  @   a   LIGHTNING_NOX_DRIVER%LIGHTNING_DT =   _.  @   a   LIGHTNING_NOX_DRIVER%LIGHTNING_START_SECONDS *   .  @   a   LIGHTNING_NOX_DRIVER%N_IC *   ß.  @   a   LIGHTNING_NOX_DRIVER%N_CG .   /  @   a   LIGHTNING_NOX_DRIVER%LNOX_OPT 2   _/  @   a   LIGHTNING_NOX_DRIVER%LNOX_PASSIVE 5   /  @   a   LIGHTNING_NOX_DRIVER%LTNG_TEMP_UPPER 5   ß/  @   a   LIGHTNING_NOX_DRIVER%LTNG_TEMP_LOWER 6   0  @   a   LIGHTNING_NOX_DRIVER%CELLCOUNT_METHOD )   _0  @   a   LIGHTNING_NOX_DRIVER%IDS )   0  @   a   LIGHTNING_NOX_DRIVER%IDE )   ß0  @   a   LIGHTNING_NOX_DRIVER%JDS )   1  @   a   LIGHTNING_NOX_DRIVER%JDE )   _1  @   a   LIGHTNING_NOX_DRIVER%KDS )   1  @   a   LIGHTNING_NOX_DRIVER%KDE )   ß1  @   a   LIGHTNING_NOX_DRIVER%IMS )   2  @   a   LIGHTNING_NOX_DRIVER%IME )   _2  @   a   LIGHTNING_NOX_DRIVER%JMS )   2  @   a   LIGHTNING_NOX_DRIVER%JME )   ß2  @   a   LIGHTNING_NOX_DRIVER%KMS )   3  @   a   LIGHTNING_NOX_DRIVER%KME )   _3  @   a   LIGHTNING_NOX_DRIVER%ITS )   3  @   a   LIGHTNING_NOX_DRIVER%ITE )   ß3  @   a   LIGHTNING_NOX_DRIVER%JTS )   4  @   a   LIGHTNING_NOX_DRIVER%JTE )   _4  @   a   LIGHTNING_NOX_DRIVER%KTS )   4  @   a   LIGHTNING_NOX_DRIVER%KTE *   ß4    a   LIGHTNING_NOX_DRIVER%C_NO 0   ó7    a   LIGHTNING_NOX_DRIVER%LNOX_TOTAL -   ;    a   LIGHTNING_NOX_DRIVER%LNOX_IC -   >    a   LIGHTNING_NOX_DRIVER%LNOX_CG 