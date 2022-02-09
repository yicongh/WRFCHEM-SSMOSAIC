MODULE MODULE_GOCART_SETTLING








CONTAINS

SUBROUTINE gocart_settling_driver(dt,config_flags,t_phy,moist,      &
           chem,rho_phy,dz8w,p8w,p_phy,                             &
           dustin,seasin,dx,g,                                      &
           dustgraset_1,dustgraset_2,dustgraset_3,                  & 
           dustgraset_4,dustgraset_5,                               & 
           setvel_1,setvel_2,setvel_3,setvel_4,setvel_5,imod,       & 
           ids,ide,jds,jde,kds,kde,                                 &
           ims,ime,jms,jme,kms,kme,                                 &
           its,ite,jts,jte,kts,kte)

  USE module_configure
  USE module_state_description
  USE module_data_gocart_dust
  USE module_data_gocart_seas
  USE module_model_constants, ONLY: mwdry
  IMPLICIT NONE

  TYPE(grid_config_rec_type), INTENT(IN) :: config_flags

  INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,                  &
                         ims,ime, jms,jme, kms,kme,                  &
                         its,ite, jts,jte, kts,kte

  REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_moist), INTENT(IN) :: moist
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_chem), INTENT(INOUT) :: chem
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: t_phy,p_phy,dz8w,p8w,rho_phy
  REAL,  DIMENSION( ims:ime , jms:jme, 5 ),                          &
          INTENT(IN   ) ::  dustin,seasin  

  INTEGER, INTENT(IN) :: imod          
  REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)::                  &
        dustgraset_1,dustgraset_2,dustgraset_3,                      &
        dustgraset_4,dustgraset_5,                                   &
        setvel_1,setvel_2,setvel_3,setvel_4,setvel_5 
                       
  REAL*8,  DIMENSION (1,1,5)          :: graset_dust,grasetvel_dust 
  REAL*8,  DIMENSION (1,1,4)          :: graset_ss,grasetvel_ss

  REAL, INTENT(IN) :: dt,dx,g
  INTEGER          :: kkk,nmx,i,j,k,kk,lmx,iseas,idust
  REAL*8, DIMENSION (1,1,kte-kts+1)   :: tmp,airden,p_mid,delz,rh
  REAL*8, DIMENSION (1,1,kte-kts+1,5) :: ddust
  REAL*8, DIMENSION (1,1,kte-kts+1,4) :: sea_salt
  
  INTEGER :: uoc_flag  






  REAL*8 conver, converi
  conver=1.e-9
  converi=1.e9
    
  uoc_flag = 0
  if (config_flags%dust_opt .eq. 4) uoc_flag = 1

  lmx=kte-kts+1

  do j=jts,jte
    do i=its,ite

       IF(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11) THEN
         do kkk=1,5
            ddust(1,1,kts,kkk)=dustin(i,j,kkk)*conver
         enddo
         kk=0
         do k=kts,kte
            kk=kk+1
            
            p_mid(1,1,kk) =.01*p_phy(i,kte-k+kts,j)
            delz(1,1,kk)  =dz8w(i,kte-k+kts,j)
            airden(1,1,kk)=rho_phy(i,k,j)
            tmp(1,1,kk)= t_phy(i,k,j)
            rh(1,1,kk) = .95
            rh(1,1,kk) = MIN( .95, moist(i,k,j,p_qv) / &
                 (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
                 (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
            rh(1,1,kk) = max(1.0D-1,rh(1,1,kk))
         enddo
       ELSE  
         kk=0       
         DO k=kts,kte
            kk=kk+1

            ddust(1,1,kk,1)=chem(i,k,j,p_dust_1)                         
            ddust(1,1,kk,2)=chem(i,k,j,p_dust_2)
            ddust(1,1,kk,3)=chem(i,k,j,p_dust_3)
            ddust(1,1,kk,4)=chem(i,k,j,p_dust_4)
            ddust(1,1,kk,5)=chem(i,k,j,p_dust_5)
            
            p_mid(1,1,kk)=.01*p_phy(i,kte-k+kts,j)
            delz(1,1,kk)=dz8w(i,kte-k+kts,j)                            
            airden(1,1,kk)=rho_phy(i,k,j)
            tmp(1,1,kk)   =t_phy(i,k,j)
            rh(1,1,kk)    = .95
            rh(1,1,kk)    = MIN( .95, moist(i,k,j,p_qv) /               &
                 (3.80*exp(17.27*(t_phy(i,k,j)-273.)/                   &
                 (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
            rh(1,1,kk)=max(1.0D-1,rh(1,1,kk))
         ENDDO
       ENDIF
       graset_dust(1,1,:)=0.
       graset_ss(1,1,:)=0.                                           

       iseas=0
       idust=1
       CALL settling(1,1,lmx,5,g,dyn_visc,ddust,tmp,p_mid,delz,       &
                     imod,graset_dust,grasetvel_dust, uoc_flag,       & 
                     den_dust,reff_dust,dt,rh,idust,iseas,airden)
          
       IF (config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) THEN
          kk=1
          do kkk=1,5
             if (kkk .le. 4) sea_salt(1,1,kts,kkk)=seasin(i,j,kkk)*conver
             if(ddust(1,1,kk,kkk) .ge. dustin(i,j,kkk)) ddust(1,1,kk,kkk)=dustin(i,j,kkk)
          enddo
          chem(i,kts,j,p_p25i)=chem(i,kts,j,p_p25i) &
                       +.25*(ddust(1,1,kk,1)+.286*ddust(1,1,kk,2))*converi
          chem(i,kts,j,p_p25i)=max(chem(i,kts,j,p_p25i),1.e-16)
          chem(i,kts,j,p_p25j)=chem(i,kts,j,p_p25j) &
                       +.75*(ddust(1,1,kk,1)+.286*ddust(1,1,kk,2))*converi
          chem(i,kts,j,p_p25j)=max(chem(i,kts,j,p_p25j),1.e-16)
          chem(i,kts,j,p_soila)=chem(i,kts,j,p_soila) &
                       +(.714*ddust(1,1,kk,2)+ddust(1,1,kk,3))*converi
          chem(i,kts,j,p_soila)=max(chem(i,kts,j,p_soila),1.e-16)
       ELSE                     
          kk = 0
          DO k = kts,kte
             kk = kk+1
             chem(i,k,j,p_dust_1)=ddust(1,1,kk,1)                         
             chem(i,k,j,p_dust_2)=ddust(1,1,kk,2)                         
             chem(i,k,j,p_dust_3)=ddust(1,1,kk,3)                         
             chem(i,k,j,p_dust_4)=ddust(1,1,kk,4)                         
             chem(i,k,j,p_dust_5)=ddust(1,1,kk,5)                         
             
             sea_salt(1,1,kk,1)=chem(i,k,j,p_seas_1)                    
             sea_salt(1,1,kk,2)=chem(i,k,j,p_seas_2)
             sea_salt(1,1,kk,3)=chem(i,k,j,p_seas_3)
             sea_salt(1,1,kk,4)=chem(i,k,j,p_seas_4)              
          ENDDO
       ENDIF



       dustgraset_1(i,j)=graset_dust(1,1,1)*airden(1,1,1)*(-1.d0)
       dustgraset_2(i,j)=graset_dust(1,1,2)*airden(1,1,1)*(-1.d0)
       dustgraset_3(i,j)=graset_dust(1,1,3)*airden(1,1,1)*(-1.d0)
       dustgraset_4(i,j)=graset_dust(1,1,4)*airden(1,1,1)*(-1.d0)
       dustgraset_5(i,j)=graset_dust(1,1,5)*airden(1,1,1)*(-1.d0)

       setvel_1(i,j)=grasetvel_dust(1,1,1)                            
       setvel_2(i,j)=grasetvel_dust(1,1,2)
       setvel_3(i,j)=grasetvel_dust(1,1,3)
       setvel_4(i,j)=grasetvel_dust(1,1,4)
       setvel_5(i,j)=grasetvel_dust(1,1,5)

       iseas=1
       idust=0

       CALL settling(1,1,lmx,4,g,dyn_visc,sea_salt,tmp,p_mid,delz,        &
                     imod,graset_ss,grasetvel_ss, uoc_flag,               &
                     den_seas,reff_seas,dt,rh,idust,iseas,airden)
       IF (config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) THEN
          kk=1
          do kkk=1,4
             if(sea_salt(1,1,kk,kkk) .ge. seasin(i,j,kkk))sea_salt(1,1,kk,kkk)=seasin(i,j,kkk)
          enddo
          chem(i,kts,j,p_naai)=chem(i,kts,j,p_naai) &
                       +.25*(sea_salt(1,1,kk,1)+.942*sea_salt(1,1,kk,2))*converi
          chem(i,kts,j,p_naai)=max(1.e-16,chem(i,kts,j,p_naai))
          chem(i,kts,j,p_naaj)=chem(i,kts,j,p_naaj) &
                       +.75*(sea_salt(1,1,kk,1)+.942*sea_salt(1,1,kk,2))*converi
          chem(i,kts,j,p_naaj)=max(1.e-16,chem(i,kts,j,p_naaj))
          chem(i,kts,j,p_seas)=chem(i,kts,j,p_seas) &
                       +(.058*sea_salt(1,1,kk,2)+sea_salt(1,1,kk,3))*converi
          chem(i,kts,j,p_seas)=max(1.e-16,chem(i,kts,j,p_seas))
       ELSE
          kk=0
          DO k=kts,kte
             kk=kk+1
             chem(i,k,j,p_seas_1)=sea_salt(1,1,kk,1)
             chem(i,k,j,p_seas_2)=sea_salt(1,1,kk,2)
             chem(i,k,j,p_seas_3)=sea_salt(1,1,kk,3)
             chem(i,k,j,p_seas_4)=sea_salt(1,1,kk,4)
          ENDDO
       ENDIF 
       
    enddo  
  enddo  

 END SUBROUTINE gocart_settling_driver


 subroutine settling(imx,jmx,lmx,nmx,g0,dyn_visc,tc,tmp,p_mid,delz, &
                    imod,graset, grasetvel, uoc,                    & 
                    den_in,reff_in,dt,rh,idust,iseas,airden)













  IMPLICIT  NONE

  INTEGER, INTENT(IN) :: imx, jmx, lmx, nmx,iseas,idust
  INTEGER             :: ntdt
  REAL,    INTENT(IN) :: dt,g0,dyn_visc
  REAL*8,  INTENT(IN) :: tmp(imx,jmx,lmx), delz(imx,jmx,lmx),  &
                         rh(imx,jmx,lmx), p_mid(imx,jmx,lmx),airden(imx,jmx,lmx)

  REAL*8, INTENT(IN)  :: den_in(nmx), reff_in(nmx)

  REAL*8, INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)

  INTEGER, INTENT(IN)   :: imod, uoc
  REAL*8, INTENT(INOUT) :: graset(imx,jmx,nmx)
  REAL*8, INTENT(OUT)   :: grasetvel(imx,jmx,nmx)
  
  REAL*8    :: den(nmx), reff(nmx)                    
  REAL*8    :: dt_settl(nmx), rcm(nmx), rho(nmx)
  INTEGER   :: ndt_settl(nmx)
  REAL*8    :: dzmin, vsettl, dtmax, pres, rhb, rwet(nmx), ratio_r(nmx)
  REAL*8    :: c_stokes, free_path, c_cun, viscosity, growth_fac
  REAL*8    :: vd_cor(lmx), vd_wk1
  INTEGER   :: k, n, i, j, l, l2

  REAL*8    :: transfer_to_below_level,temp_tc


  REAL*8, PARAMETER :: c1=0.7674, c2=3.079, c3=2.573E-11, c4=-1.424 


  REAL*8 :: rwet_priv(nmx), rho_priv(nmx)


  IF ( idust.ne.1 .and. iseas.ne.1 ) RETURN
  WHERE ( tc(:,:,:,:) < 0.0 ) tc(:,:,:,:) = 1.0D-32

  den = den_in  
  reff = reff_in  
  dzmin = MINVAL(delz(:,:,:))
  IF (idust == 1) growth_fac = 1.0
  IF (iseas == 1) growth_fac = 3.0

  IF (idust == 1 .and. uoc == 1) then
      den(1) = 2650.            

  ENDIF

  DO k = 1, nmx                                                   









     vsettl = 4.0/9.0 * g0 * den(k) * (growth_fac*reff(k))**2 / dyn_visc



     ntdt = INT(dt)
     dtmax = dzmin / vsettl
     ndt_settl(k) = MAX( 1,INT(ntdt/dtmax) )


     IF (ndt_settl(k) > 12) ndt_settl(k) = 12
     dt_settl(k) = REAL(ntdt) / REAL(ndt_settl(k))


     IF (iseas.eq.1)rcm(k) = reff(k)*100.0
     IF (idust.eq.1) THEN
      rwet(k) = reff(k)
      ratio_r(k) = 1.0
      rho(k) = den(k)
     ENDIF
  ENDDO



!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( i,   j,   l,   l2, n,   k,   rhb, rwet_priv, ratio_r, c_stokes)&
!$OMP PRIVATE( free_path, c_cun, viscosity, rho_priv, vd_cor )



  DO j = 1,jmx                      
    DO i = 1,imx                    
      DO k = 1,nmx                  
        graset(i,j,k)=0.
        grasetvel(i,j,k)=0.
    
        IF (idust.eq.1) THEN
           rwet_priv(k) = rwet(k)
           rho_priv(k)  = rho(k)
        END IF

        DO n = 1,ndt_settl(k)        

          transfer_to_below_level=0

          DO l = lmx,1,-1            
             l2 = lmx - l + 1

             IF (iseas.eq.1) THEN
                rhb = MIN(9.9D-1, rh(i,j,l))                             
                rwet_priv(k) = 0.01*(c1*rcm(k)**c2/(c3*rcm(k)**c4 -  &
                LOG10(rhb)) + rcm(k)**3)**0.33                           
                ratio_r(k) = (reff(k)/rwet_priv(k))**3.0
             END IF

             c_stokes = 1.458E-6*tmp(i,j,l)**1.5/(tmp(i,j,l) + 110.4)    
             free_path = 1.1E-3/p_mid(i,j,l2)/SQRT(tmp(i,j,l))           

             c_cun = 1.0+free_path/rwet_priv(k)*                     &
             (1.257 + 0.4*EXP(-1.1*rwet_priv(k)/free_path))              
             viscosity = c_stokes / c_cun                                

             IF (iseas.eq.1) THEN
                rho_priv(k) = ratio_r(k)*den(k) + (1.-ratio_r(k))*1000.
             END IF

             vd_cor(l) = 2./9.*g0*rho_priv(k)*rwet_priv(k)**2/viscosity  

            
            temp_tc=tc(i,j,l,k)                                          
            vd_wk1 = dt_settl(k)*vd_cor(l)/delz(i,j,l2)                  

            tc(i,j,l,k)   =  tc(i,j,l,k)*(1.- vd_wk1)+transfer_to_below_level          
            transfer_to_below_level = (temp_tc*vd_wk1)*((delz(i,j,l2)*airden(i,j,l))/(delz(i,j,l2+1)*airden(i,j,l-1))) 

             IF (l==1) THEN
                graset(i,j,k)=graset(i,j,k)+vd_cor(l)*(temp_tc*vd_wk1)/ndt_settl(k) 
                grasetvel(i,j,k)=vd_cor(l)                                     
             ENDIF
          ENDDO                 
        ENDDO                   
      ENDDO                     

    ENDDO                     
  END DO                      
!$OMP END PARALLEL DO

END SUBROUTINE settling

END MODULE MODULE_GOCART_SETTLING
