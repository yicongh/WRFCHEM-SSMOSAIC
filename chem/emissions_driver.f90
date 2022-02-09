






MODULE module_emissions_driver
   IMPLICIT NONE
CONTAINS

    subroutine emissions_driver(id,ktau,dtstep,DX,                         &
         adapt_step_flag,curr_secs,                                        &
         plumerisefire_frq,stepfirepl,                                     &
         bioemdt,stepbioe,                                                 &
         config_flags,gmt,julday,alt,t_phy,moist,p8w,t8w,u_phy,            &
         v_phy,vvel,e_bio,p_phy,chem,rho_phy,dz8w,ne_area,emis_ant,        &
         emis_vol,tsk,erod,erod_dri,lai_vegmask,                           &
         g,emis_seas,emis_dust,tracer,                                     &
         emis_seas2,                                                       &
         ebu, ebu_in,mean_fct_agtf,mean_fct_agef,                          &
         mean_fct_agsv,mean_fct_aggr,firesize_agtf,firesize_agef,          &
         firesize_agsv,firesize_aggr,                                      &
         u10,v10,ivgtyp,isltyp,gsw,vegfra,rmol,ust,znt,dms_0,              &
         erup_beg,erup_end,                                                &
         xland,xlat,xlong,z_at_w,z,smois,dustin,seasin,                    &
         sebio_iso,sebio_oli,sebio_api,sebio_lim,sebio_xyl,                &
         sebio_hc3,sebio_ete,sebio_olt,sebio_ket,sebio_ald,                &   
         sebio_hcho,sebio_eth,sebio_ora2,sebio_co,sebio_nr,                &
         sebio_sesq,sebio_mbo,                                             & 
         noag_grow,noag_nongrow,nononag,slai,                              &
         ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,                     &
         ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,                     &
         ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,             &
         ebio_sesq, ebio_mbo,ebio_bpi,ebio_myrc,                           &
         ebio_c10h16,ebio_tol,ebio_bigalk,ebio_ch3oh,ebio_acet,            &
         ebio_nh3,ebio_no2,ebio_c2h5oh,ebio_ch3cooh,ebio_mek,              &
         ebio_bigene,ebio_c2h6,ebio_c2h4,ebio_c3h6,ebio_c3h8,ebio_so2,     &
         ebio_dms,ebio_hcn,                                                &
         ebio_alk3, ebio_alk4, ebio_alk5, ebio_ole1, ebio_ole2,            &    
         ebio_aro1, ebio_aro2, ebio_ccho, ebio_meoh,                       &    
         ebio_ethene, ebio_hcooh, ebio_terp, ebio_bald,                    &    
         ebio_cco_oh, ebio_rco_oh,                                         &    
         clayfrac,sandfrac,dust_alpha,dust_gamma,dust_smtune,dust_ustune,  &
         clayfrac_nga,sandfrac_nga,                                        &
         snowh,zs,afwa_dustloft,tot_dust,tot_edust,vis_dust,               &
         soil_top_cat, ust_t, rough_cor, smois_cor,                        & 
         ebio_c5h8,ebio_apinene,ebio_bpinene,ebio_toluene,                 &
         ebio_ch3cho,ebio_ch3co2h,ebio_tbut2ene,ebio_c2h5cho,              &
         ebio_nc4h10,							   &
         
         T2,swdown,                                                        &
         nmegan,EFmegan,                                                   &
         msebio_isop,                                                      &
         mlai,                                                             &
         pftp_bt, pftp_nt, pftp_sb, pftp_hb,                               &
         mtsa,                                                             &
         mswdown,                                                          &
         mebio_isop, mebio_apin, mebio_bpin, mebio_bcar,                   &
         mebio_acet, mebio_mbo, mebio_no,                                  &
         current_month,                                                    &
         
         
         ht, refl_10cm,                                                    &
         ic_flashrate, cg_flashrate,                                       &
         
         
         emis_aircraft,                                                    &
         
         vprm_in,rad_vprm,lambda_vprm,alpha_vprm,resp_vprm,               &
         xtime,tslb,wet_in,rainc,rainnc,potevp,sfcevp,lu_index,            &
         biomt_par,emit_par,ebio_co2oce,eghg_bio,                          &
         dust_flux, seas_flux,                                             &

         plume_hgt_agtf,plume_hgt_agef,plume_hgt_agsv,plume_hgt_aggr,      &

         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )

  USE module_configure
  USE module_state_description
  USE module_data_radm2
  USE module_data_sorgam, only : mw_so4_aer,anthfac,factnumn,factnuma,factnumc
  USE module_model_constants, only : mwdry
  USE module_emissions_anthropogenics
  USE module_bioemi_simple
  USE module_bioemi_beis314
  USE module_bioemi_megan2
  USE module_aerosols_sorgam, only: sorgam_addemiss
  USE module_cbmz_addemiss
  USE module_cb05_addemiss
  USE gocart_dust
  USE gocart_dust_afwa
  USE gocart_seasalt
  USE uoc_dust  
  USE module_dms_emis
  USE module_mosaic_addemiss
  USE module_add_emis_cptec
  USE module_add_emiss_burn
  USE module_plumerise1
  USE module_aerosols_sorgam_vbs, only: sorgam_vbs_addemiss
  USE module_aerosols_soa_vbs, only: soa_vbs_addemiss
  USE module_ghg_fluxes
  USE module_lightning_nox_driver
  USE module_cam_mam_addemiss, only: cam_mam_addemiss

  USE shr_megan_mod, only : shr_megan_mechcomps_n, shr_megan_mechcomps

  
  IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER,      INTENT(IN   ) :: id,julday, ne_area,                      &
                                  ids,ide, jds,jde, kds,kde,               &
                                  ims,ime, jms,jme, kms,kme,               &
                                  its,ite, jts,jte, kts,kte
   INTEGER,      INTENT(IN   ) ::                                          &
                                  ktau,stepbioe,stepfirepl
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),                &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),                 &
         INTENT(INOUT ) ::                                   chem
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_tracer ),               &
         INTENT(INOUT ) ::                                   tracer
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_ebu ),                  &
         INTENT(INOUT ) ::                                   ebu
   REAL, DIMENSION( ims:ime, 1, jms:jme, num_ebu_in ),                     &
         INTENT(INOUT ) ::                                   ebu_in
   REAL, DIMENSION( ims:ime, jms:jme, ne_area ),                           &
         INTENT(INOUT ) ::                               e_bio

   REAL, DIMENSION( ims:ime, jms:jme),                     &
         INTENT(INOUT ) :: plume_hgt_agtf,plume_hgt_agef,plume_hgt_agsv,plume_hgt_aggr

   REAL, DIMENSION( ims:ime, 1:config_flags%kemit, jms:jme,num_emis_ant),&
         INTENT(IN ) ::                                                    &
         emis_ant
   REAL, DIMENSION( ims:ime,  kms:kme, jms:jme,num_emis_vol),              &
         INTENT(INOUT ) ::                                                 &
         emis_vol
   REAL, DIMENSION( ims:ime, jms:jme),&
         INTENT(IN ) ::                                                 &
         dms_0,tsk,erup_beg,erup_end
   REAL, DIMENSION( ims:ime, jms:jme,3),&
         INTENT(IN ) ::                                                 &
         erod, erod_dri
   REAL, DIMENSION( ims:ime, jms:jme), &
         INTENT(IN ) ::                                                    &
         lai_vegmask
   REAL, DIMENSION( ims:ime, jms:jme,5),&
         INTENT(INOUT ) ::                                                 &
         dustin,seasin
   REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_dust),   &
         OPTIONAL, INTENT(INOUT ) ::                                       &
         emis_dust
   REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_seas),   &
         OPTIONAL,                                                         &
         INTENT(INOUT ) ::                                                 &
         emis_seas
   REAL, DIMENSION( ims:ime, 1, jms:jme,num_emis_seas2),   &
         OPTIONAL,                                                         &
         INTENT(INOUT ) ::                                                 &
         emis_seas2
   REAL, DIMENSION( ims:ime,  jms:jme ),                                   &
         OPTIONAL,                                                         &
         INTENT(IN ) ::                                                    &
           mean_fct_agtf,mean_fct_agef,                                    &
           mean_fct_agsv,mean_fct_aggr,firesize_agtf,firesize_agef,        &
           firesize_agsv,firesize_aggr




   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,               &
          INTENT(IN   ) ::                                                 &
                                                        alt,               &
                                                      t_phy,               &
                                                      p_phy,               &
                                                      dz8w,                &
                                              t8w,p8w,z_at_w , z ,         &
                                              u_phy,v_phy,vvel,rho_phy
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,               &
          INTENT(IN   ) ::                                                 &
                                                     ivgtyp,               &
                                                     isltyp
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
          INTENT(IN   ) ::                                                 &
                                                     u10,                  &
                                                     v10,                  &
                                                     gsw,                  &
                                                  vegfra,                  &
                                                     rmol,                 &
                                                     ust,                  &
                                                     xland,                &
                                                     xlat,                 &
                                                     xlong,                &
                                                     znt,                  &


                                                     rainc,                &
                                                     rainnc,               &
                                                     potevp,               &
                                                     sfcevp,               &
                                                     lu_index

   REAL, DIMENSION( ims:ime , jms:jme )                  ,                 &
         OPTIONAL,                                                         &
         INTENT(IN   ) ::                                                  &
                                                     clayfrac,             &
                                                     sandfrac,             &
                                                     clayfrac_nga,         &
                                                     sandfrac_nga,         &
                                                     snowh
   REAL, INTENT(IN   ) ::                            dust_alpha,           &
                                                     dust_gamma,           &
                                                     dust_smtune,          &
                                                     dust_ustune

  REAL, DIMENSION( config_flags%num_soil_layers ) ,                        &
      INTENT(IN   ) ::                               zs
  REAL, DIMENSION( ims:ime , jms:jme )                   ,                 &
         OPTIONAL,                                                         &
         INTENT(  OUT) ::                                                  &
                                                     tot_edust,            &
                                                     afwa_dustloft
   REAL, DIMENSION( ims:ime , kms:kme , jms:jme )        ,                 &
         OPTIONAL,                                                         &
         INTENT(  OUT) ::                                                  &
                                                     tot_dust,             &
                                                     vis_dust
  REAL, DIMENSION( ims:ime, config_flags%num_soil_layers, jms:jme ) ,      &
      INTENT(INOUT ) ::                             smois, tslb    

   REAL,  DIMENSION( ims:ime , jms:jme )                   ,               &
         OPTIONAL,                                                         &
          INTENT(INOUT   ) ::                                                 &
               sebio_iso,sebio_oli,sebio_api,sebio_lim,sebio_xyl,      &
               sebio_hc3,sebio_ete,sebio_olt,sebio_ket,sebio_ald,      &
               sebio_hcho,sebio_eth,sebio_ora2,sebio_co,sebio_nr,      &
               sebio_sesq,sebio_mbo,                                   & 
               noag_grow,noag_nongrow,nononag,slai,                    &
               ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,           &
               ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,           &
               ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,   &
               ebio_sesq,ebio_mbo,ebio_bpi,ebio_myrc,                  &
               ebio_c10h16,ebio_tol,ebio_bigalk, ebio_ch3oh,ebio_acet, &
               ebio_nh3,ebio_no2,ebio_c2h5oh,ebio_ch3cooh,ebio_mek,    &
               ebio_bigene,ebio_c2h6,ebio_c2h4,ebio_c3h6,ebio_c3h8,    &
               ebio_so2,ebio_dms, ebio_co2oce , ebio_hcn,              &
               ebio_alk3, ebio_alk4, ebio_alk5, ebio_ole1, ebio_ole2,  &    
               ebio_aro1, ebio_aro2, ebio_ccho, ebio_meoh,             &    
               ebio_ethene, ebio_hcooh, ebio_terp, ebio_bald,          &    
               ebio_cco_oh, ebio_rco_oh,                               &
               ebio_c5h8,ebio_apinene,ebio_bpinene,ebio_toluene,       &
               ebio_ch3cho,ebio_ch3co2h,ebio_tbut2ene,ebio_c2h5cho,    &
               ebio_nc4h10
               
   REAL,  DIMENSION( ims:ime , jms:jme ) , OPTIONAL  ,                 &
          INTENT(INOUT) ::                                    ust_t,   & 
                                                          rough_cor,   &
                                                          smois_cor          

   REAL, DIMENSION(ims:ime,1:config_flags%num_soil_cat,jms:jme) ,      & 
          INTENT(IN)::                                 soil_top_cat 

   
   

   integer, intent(in   ) :: nmegan
   real, dimension (ims:ime, jms:jme , nmegan) ,                       &
         OPTIONAL,                                                         &
        intent(inout) ::                                               &
        EFmegan


   real, dimension (ims:ime, jms:jme ) ,                               &
         OPTIONAL,                                                         &
        intent(in) ::                                                  &
        msebio_isop,                                                   &
        pftp_bt, pftp_nt, pftp_sb, pftp_hb

   real, dimension (ims:ime, jms:jme, 12 ) ,                           &
         OPTIONAL,                                                         &
        intent(in) ::                                                  &
        mlai, mtsa, mswdown

   real, dimension (ims:ime, jms:jme ) ,                               &
         OPTIONAL,                                                         &
        intent(inout) ::                                               &
        mebio_isop, mebio_apin, mebio_bpin, mebio_bcar,                &
        mebio_acet, mebio_mbo, mebio_no

   real, dimension (ims:ime, jms:jme ) ,                               &
        intent(in) ::                                                  &
        T2, swdown

   integer, intent(in) :: current_month

   

      REAL(KIND=8), INTENT(IN   ) ::                                   &
           curr_secs

      REAL :: gmtp,gmtm
      integer :: curr_hours,ivolcano
      Integer :: endhr,endmin,beghr,begmin,ko,kk4,kl,k_initial,k_final
      real :: emiss_ash_mass,emiss_ash_height,so2_mass,vert_mass_dist(kts:kte)
      real :: eh
      real :: area,x1,percen_mass_umbrel,base_umbrel,ashz_above_vent

      REAL, INTENT(IN   ) ::                                           &
           bioemdt, dtstep, dx, gmt, g

      INTEGER, INTENT(IN   ) ::                                        &
           plumerisefire_frq

      LOGICAL, INTENT(IN   ) ::                                        &
           adapt_step_flag


      REAL, DIMENSION( ims:ime, 1:config_flags%kemit_aircraft, jms:jme,num_emis_aircraft), &
            OPTIONAL, INTENT(IN ) :: emis_aircraft 



      REAL, DIMENSION(ims:ime, 8, jms:jme, num_vprm_in), INTENT(IN)     ::  vprm_in
      REAL, DIMENSION(ims:ime, 1,jms:jme, num_eghg_bio), INTENT(INOUT ) ::  eghg_bio

      REAL, DIMENSION(ims:ime, jms:jme), INTENT(OUT)   :: dust_flux, seas_flux

      REAL, DIMENSION(8) :: rad_vprm,lambda_vprm,alpha_vprm,resp_vprm

      REAL, DIMENSION(14), INTENT(IN) :: biomt_par, emit_par
      REAL, DIMENSION(ims:ime,1,jms:jme,num_wet_in), INTENT(IN)  :: wet_in
      REAL, INTENT(IN) :: xtime


     REAL, DIMENSION( ims:ime,          jms:jme ),           INTENT(IN   ) :: ht, ic_flashrate, cg_flashrate
     REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: refl_10cm




      INTEGER :: begday,endday,i, j, k, m, p_in_chem, ksub, dust_emiss_active, seasalt_emiss_active,emiss_ash_hgt
      REAL :: conv,conv3,conv4,oconv3,oconv4

      REAL :: convert2(its:ite,jts:jte)

      CHARACTER (LEN=80) :: message
      LOGICAL :: do_bioemiss, do_plumerisefire,do_ex_volcanoe

      INTEGER :: imod  

            











       percen_mass_umbrel=.75
       base_umbrel=.25    

      ivolcano=0
      area=dx*dx
      dust_emiss_active    = 0
      seasalt_emiss_active = 0
      if(config_flags%dust_opt >= 2 )dust_emiss_active    = 1
      if(config_flags%seas_opt == 2 )seasalt_emiss_active = 1
      if(config_flags%seas_opt == 3 )seasalt_emiss_active = 3
      if(config_flags%seas_opt == 4 )seasalt_emiss_active = 4




      gmtp=curr_secs/3600.
      curr_hours=curr_secs/3600.
      gmtp=mod(gmt+gmtp,24.)
      gmtm=mod(gmtp,60.)


      if(config_flags%emiss_opt_vol == 1 .or. config_flags%emiss_opt_vol == 2)then
         do_ex_volcanoe = .false.

      emiss_ash_height = config_flags%emiss_ash_hgt
      if(emiss_ash_height.gt. 1.)then
      write(message,'(" ADJUSTED ASH HEIGHT: ",2f15.3)') emiss_ash_height, area
      CALL WRF_DEBUG (15,message)


      do j=jts,jte
      do i=its,ite
        if(erup_end(i,j).gt.0)then
        so2_mass=1.5e4*3600.*1.e9/64./area
        eh=2600.*(emiss_ash_height*.0005)**4.1494
        emiss_ash_mass=eh*1.e9/area
             


        ashz_above_vent=emiss_ash_height - z_at_w(i,kts,j)                              
      write(message,'("Found and adjusted active volcano at j,kts,kpe = ",3i8)') j,kts,kte 
      call wrf_message (message)

        do k=kte-1,kts,-1
           if(z_at_w(i,k,j) < emiss_ash_height)then                                    
             k_final=k+1
             exit
           endif 
        enddo
        do k=kte-1,kts,-1
          if(z_at_w(i,k,j) < ((1.-base_umbrel)*ashz_above_vent)+z_at_w(i,kts,j))then  
             k_initial=k
             exit
           endif
        enddo
        vert_mass_dist=0.

      
          kk4 = k_final-k_initial+2
          do ko=1,kk4-1
              kl=ko+k_initial-1
              vert_mass_dist(kl) = 6.*percen_mass_umbrel* float(ko)    &
                           /float(kk4)**2 * (1. - float(ko)/float(kk4))
          enddo
          if(sum(vert_mass_dist(kts:kte)) .ne. percen_mass_umbrel) then
            x1= ( percen_mass_umbrel- sum(vert_mass_dist(kts:kte)) )   &
                 /float(k_final-k_initial+1)
              do ko=k_initial,k_final

                vert_mass_dist(ko) = vert_mass_dist(ko)+ x1
              enddo
          endif 


          do ko=1,k_initial-1
             vert_mass_dist(ko)=float(ko)/float(k_initial-1)
          enddo
          x1=sum(vert_mass_dist(1:k_initial-1))
          do ko=1,k_initial-1
              vert_mass_dist(ko)=(1.-percen_mass_umbrel)*vert_mass_dist(ko)/x1
          enddo
          do ko=1,k_final
            emis_vol(i,ko,j,p_e_vash1)=.22*vert_mass_dist(ko)*emiss_ash_mass
            emis_vol(i,ko,j,p_e_vash2)=.05*vert_mass_dist(ko)*emiss_ash_mass
            emis_vol(i,ko,j,p_e_vash3)=.4*vert_mass_dist(ko)*emiss_ash_mass
            emis_vol(i,ko,j,p_e_vash4)=.05*vert_mass_dist(ko)*emiss_ash_mass
            emis_vol(i,ko,j,p_e_vash5)=.245*vert_mass_dist(ko)*emiss_ash_mass
            emis_vol(i,ko,j,p_e_vash6)=.12*vert_mass_dist(ko)*emiss_ash_mass
            emis_vol(i,ko,j,p_e_vash7)=.11*vert_mass_dist(ko)*emiss_ash_mass
            emis_vol(i,ko,j,p_e_vash8)=.08*vert_mass_dist(ko)*emiss_ash_mass
            emis_vol(i,ko,j,p_e_vash9)=.05*vert_mass_dist(ko)*emiss_ash_mass
            emis_vol(i,ko,j,p_e_vash10)=.035*vert_mass_dist(ko)*emiss_ash_mass
            if(config_flags%emiss_opt_vol == 2)emis_vol(i,ko,j,p_e_vso2)=vert_mass_dist(ko)*so2_mass
          enddo
          do ko=k_final+1,kte
           emis_vol(i,ko,j,p_e_vash1)=0.
           emis_vol(i,ko,j,p_e_vash2)=0.
           emis_vol(i,ko,j,p_e_vash3)=0.
           emis_vol(i,ko,j,p_e_vash4)=0.
           emis_vol(i,ko,j,p_e_vash5)=0.
           emis_vol(i,ko,j,p_e_vash6)=0.
           emis_vol(i,ko,j,p_e_vash7)=0.
           emis_vol(i,ko,j,p_e_vash8)=0.
           emis_vol(i,ko,j,p_e_vash9)=0.
           emis_vol(i,ko,j,p_e_vash10)=0.
           if(config_flags%emiss_opt_vol == 2)emis_vol(i,ko,j,p_e_vso2)=0.
         enddo
      endif  
      enddo 
      enddo 
     else


     endif 






      do j=jts,jte
      do i=its,ite
         ivolcano = 0
        if(erup_end(i,j).le.0)cycle




         begday=int(erup_beg(i,j)/1000.)-1
         beghr=int(erup_beg(i,j))-(begday+1)*1000
         begmin=00.
         endhr=beghr+int(erup_end(i,j)/60.)
         endday=int(begday+endhr/24)-1
         endmin=00.


         ivolcano = 1
         if(julday.le.begday .or. julday.ge.endday)then

            if( julday.lt.begday)then
                 write(message,'("before volcano stuff at julday = ",i8)') julday
                 call wrf_debug(15,message)
                 ivolcano=0
            elseif(julday.eq.begday)then
               if(beghr.gt.int(gmtp))then
                 write(message,'("before volcano stuff at gmtp = ",i8)') int(gmtp)
                 call wrf_debug(15,message)
                 ivolcano=0
               elseif(beghr.eq.int(gmtp))then
                  if(begmin.gt.gmtm)then
                     write(message,'("before volcano stuff at gmtp,begmin = ",2i8)') int(gmtp),int(begmin)
                     call wrf_debug(15,message)
                     ivolcano=0
                  endif
               endif
            endif
            if( julday.gt.endday)then
                 write(message,'("after volcano stuff at julday = ",i8)') julday
                 call wrf_debug(15,message)
                 ivolcano=0
            elseif(julday.eq.endday)then
               if(endhr.lt.int(gmtp))then
                 write(message,'("after volcano stuff at gmtp = ",i8)') int(gmtp)
                 call wrf_debug(15,message)
                 ivolcano=0
               elseif(endhr.eq.int(gmtp))then
                  if(endmin.lt.gmtm)then
                     write(message,'("after volcano stuff at gmtm,endmin = ",2i8)') int(gmtm),int(endmin)
                     call wrf_debug(15,message)
                     ivolcano=0
                  endif
               endif
            endif
         endif 

      volc_select:  SELECT CASE(config_flags%chem_opt)
      CASE (GOCART_SIMPLE,MOZCART_KPP,T1_MOZCART_KPP,GOCARTRADM2,GOCARTRACM_KPP)
        CALL wrf_debug(15,'Adding volcanic emissions')
                  do k=kts,kte
                    conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                    chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                   &
                     +emis_vol(i,k,j,p_e_vso2)*conv
                  enddo
               do k=kts,kte
                conv=float(ivolcano)*alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_p25)=chem(i,k,j,p_p25)+.5*emis_vol(i,k,j,p_e_vash10)*conv
                chem(i,k,j,p_p10)=chem(i,k,j,p_p10)     &
                                 +.5*emis_vol(i,k,j,p_e_vash10)*conv &
                                 +emis_vol(i,k,j,p_e_vash9)*conv    &
                                 +.5*emis_vol(i,k,j,p_e_vash8)*conv
               enddo
      CASE (RADM2SORG,RADM2SORG_AQ,RADM2SORG_KPP,RACMSORG_KPP,RACMSORG_AQ,RACM_ESRLSORG_KPP, &
            RACMSORG_AQCHEM_KPP,RACM_ESRLSORG_AQCHEM_KPP)

                  do k=kts,kte
                    conv = float(ivolcano)*4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                    chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                 &
                     +emis_vol(i,k,j,p_e_vso2)*conv

                    conv=alt(i,k,j)*dtstep/dz8w(i,k,j)
                    chem(i,k,j,p_p25i) = chem(i,k,j,p_p25i)                &
                     +.25*emis_vol(i,k,j,p_e_vash10)*conv
                    chem(i,k,j,p_nu0) = chem(i,k,j,p_nu0)                  &
                     +.25*anthfac*factnumn*emis_vol(i,k,j,p_e_vash10)*conv
                    chem(i,k,j,p_ac0) = chem(i,k,j,p_ac0)                  &
                     +.75*anthfac*factnuma*emis_vol(i,k,j,p_e_vash10)*conv
                    chem(i,k,j,p_p25j) = chem(i,k,j,p_p25j)                &
                     +.75*emis_vol(i,k,j,p_e_vash10)*conv
                    chem(i,k,j,p_antha) = chem(i,k,j,p_antha)              &
                     +emis_vol(i,k,j,p_e_vash9)*conv 

                    chem(i,k,j,p_corn) = chem(i,k,j,p_corn)                &
                     +anthfac*factnumc*emis_vol(i,k,j,p_e_vash9)*conv 

                  enddo
      CASE (CHEM_VOLC)
              CALL wrf_debug(15,'Adding volcanic emissions to case chem_volc')
               do k=kts,kte
                 conv = float(ivolcano)*4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                 chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                   &
                  +emis_vol(i,k,j,p_e_vso2)*conv
               enddo
               do k=kts,kte
                conv=float(ivolcano)*alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_vash_1)=chem(i,k,j,p_vash_1)+emis_vol(i,k,j,p_e_vash1)*conv
                chem(i,k,j,p_vash_2)=chem(i,k,j,p_vash_2)+emis_vol(i,k,j,p_e_vash2)*conv
                chem(i,k,j,p_vash_3)=chem(i,k,j,p_vash_3)+emis_vol(i,k,j,p_e_vash3)*conv
                chem(i,k,j,p_vash_4)=chem(i,k,j,p_vash_4)+emis_vol(i,k,j,p_e_vash4)*conv
                chem(i,k,j,p_vash_5)=chem(i,k,j,p_vash_5)+emis_vol(i,k,j,p_e_vash5)*conv
                chem(i,k,j,p_vash_6)=chem(i,k,j,p_vash_6)+emis_vol(i,k,j,p_e_vash6)*conv
                chem(i,k,j,p_vash_7)=chem(i,k,j,p_vash_7)+emis_vol(i,k,j,p_e_vash7)*conv
                chem(i,k,j,p_vash_8)=chem(i,k,j,p_vash_8)+emis_vol(i,k,j,p_e_vash8)*conv
                chem(i,k,j,p_vash_9)=chem(i,k,j,p_vash_9)+emis_vol(i,k,j,p_e_vash9)*conv
                chem(i,k,j,p_vash_10)=chem(i,k,j,p_vash_10)+emis_vol(i,k,j,p_e_vash10)*conv
               enddo
      CASE (CHEM_VOLC_4BIN)
               CALL wrf_debug(15,'Adding volcanic emissions to case chem_volc_4bin')
               do k=kts,kte
                conv=float(ivolcano)*alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_vash_7)=chem(i,k,j,p_vash_7)+emis_vol(i,k,j,p_e_vash7)*conv
                chem(i,k,j,p_vash_8)=chem(i,k,j,p_vash_8)+emis_vol(i,k,j,p_e_vash8)*conv
                chem(i,k,j,p_vash_9)=chem(i,k,j,p_vash_9)+emis_vol(i,k,j,p_e_vash9)*conv
                chem(i,k,j,p_vash_10)=chem(i,k,j,p_vash_10)+emis_vol(i,k,j,p_e_vash10)*conv
               enddo
      CASE (CHEM_VASH)
               CALL wrf_debug(15,'Adding volcanic emissions to case chem_volc')
               do k=kts,kte
                conv=float(ivolcano)*alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_vash_1)=chem(i,k,j,p_vash_1)+emis_vol(i,k,j,p_e_vash1)*conv
                chem(i,k,j,p_vash_2)=chem(i,k,j,p_vash_2)+emis_vol(i,k,j,p_e_vash2)*conv
                chem(i,k,j,p_vash_3)=chem(i,k,j,p_vash_3)+emis_vol(i,k,j,p_e_vash3)*conv
                chem(i,k,j,p_vash_4)=chem(i,k,j,p_vash_4)+emis_vol(i,k,j,p_e_vash4)*conv
                chem(i,k,j,p_vash_5)=chem(i,k,j,p_vash_5)+emis_vol(i,k,j,p_e_vash5)*conv
                chem(i,k,j,p_vash_6)=chem(i,k,j,p_vash_6)+emis_vol(i,k,j,p_e_vash6)*conv
                chem(i,k,j,p_vash_7)=chem(i,k,j,p_vash_7)+emis_vol(i,k,j,p_e_vash7)*conv
                chem(i,k,j,p_vash_8)=chem(i,k,j,p_vash_8)+emis_vol(i,k,j,p_e_vash8)*conv
                chem(i,k,j,p_vash_9)=chem(i,k,j,p_vash_9)+emis_vol(i,k,j,p_e_vash9)*conv
                chem(i,k,j,p_vash_10)=chem(i,k,j,p_vash_10)+emis_vol(i,k,j,p_e_vash10)*conv
               enddo
      CASE DEFAULT
      END SELECT volc_select

      enddo
      enddo
      ENDIF

      do_plumerisefire = .false.
      IF ( config_flags%biomass_burn_opt == BIOMASSB_MOZC .OR. &
           config_flags%biomass_burn_opt == BIOMASSB_T1_MOZCART  .OR. &
           config_flags%biomass_burn_opt == BIOMASSB_MOZ  .OR. &
           config_flags%biomass_burn_opt == BIOMASSB_SAPRC  .OR. &
           config_flags%biomass_burn_opt == BIOMASSB_GHG  .OR. &
           config_flags%biomass_burn_opt == BIOMASSB ) then
        IF ( ktau==1 ) then
           do_plumerisefire = .true.
        ELSE IF ( adapt_step_flag ) THEN
           IF ( (plumerisefire_frq<=0) .or. &
                ( curr_secs+real(dtstep,8)+0.01 >= &
                ( INT( curr_secs/real(plumerisefire_frq*60.,8)+1,8 )*real(plumerisefire_frq*60.,8) ) ) &
                ) then
              do_plumerisefire = .true.
           ENDIF
        ELSE IF ( (MOD(ktau,stepfirepl)==0) .or. (stepfirepl==1) ) THEN
           do_plumerisefire = .true.
        ENDIF
      ENDIF

      do_bioemiss = .false.
      IF ( ktau==1 ) then
         do_bioemiss = .true.
      ELSE IF ( adapt_step_flag ) THEN
         IF ( (bioemdt<=0) .or. &
              ( curr_secs+real(dtstep,8)+0.01 >= &
              ( INT( curr_secs/real(bioemdt*60.,8)+1,8 )*real(bioemdt*60.,8) ) ) &
              ) then
            do_bioemiss = .true.
         ENDIF
      ELSE IF ( (MOD(ktau,stepbioe)==0) .or. (stepbioe==1) ) THEN
         do_bioemiss = .true.
      ENDIF



       if( do_plumerisefire )then
          CALL wrf_debug(15,'fire emissions: calling biomassb')
          write(0,*)ktau,stepfirepl
          call plumerise_driver (id,ktau,dtstep,                           &
           ebu,ebu_in,                                                     &
           mean_fct_agtf,mean_fct_agef,mean_fct_agsv,mean_fct_aggr,        &
           firesize_agtf,firesize_agef,firesize_agsv,firesize_aggr,        &
           config_flags, t_phy,moist,                                      &
           rho_phy,vvel,u_phy,v_phy,p_phy,                                 &
           emis_ant,z_at_w,z,config_flags%scale_fire_emiss,                &

           plume_hgt_agtf,plume_hgt_agef,plume_hgt_agsv,plume_hgt_aggr,    &

           ids,ide, jds,jde, kds,kde,                                      &
           ims,ime, jms,jme, kms,kme,                                      &
           its,ite, jts,jte, kts,kte                                       )
        endif



      tracer_select:  SELECT CASE(config_flags%tracer_opt)
      CASE (TRACER_SMOKE,TRACER_TEST2)
          CALL wrf_debug(15,'tracer fire emissions: calling biomassb, only CO')



       call add_emis_burn(id,dtstep,ktau,dz8w,rho_phy,tracer,&
            julday,gmt,xlat,xlong,t_phy,p_phy,                           &
            ebu,0,config_flags%tracer_opt,config_flags%biomass_burn_opt,     &
            num_tracer,ids,ide, jds,jde, kds,kde,                                   &
            ims,ime, jms,jme, kms,kme,                                   &
            its,ite, jts,jte, kts,kte                                    )
      CASE DEFAULT
        CALL wrf_debug(15,'No tracer option selected')
      END SELECT tracer_select




      seasalt_select:  SELECT CASE(config_flags%seas_opt)
      CASE (SEASGOCART)
        CALL wrf_debug(15,'Gocart sea salt emissions')
         call gocart_seasalt_driver(ktau,dtstep,config_flags,julday,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,u10,v10,p8w,z_at_w,                  &
         xland,xlat,xlong,dx,g,emis_seas, seasin,&
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )

      CASE DEFAULT 
        if(seasalt_emiss_active.eq.1) then 
           CALL wrf_debug(15,'MOSAIC or SORGAM sea salt emissions')
        elseif(seasalt_emiss_active.eq.3) then
           CALL wrf_debug(15,'MOSAIC sea salt emissions (Fuentes et al) - low activity')
        elseif(seasalt_emiss_active.eq.4) then
           CALL wrf_debug(15,'MOSAIC sea salt emissions (Fuentes et al) - high activity')
        else
           CALL wrf_debug(15,'no sea salt emissions')
        end if
      END SELECT seasalt_select

      dust_select:  SELECT CASE(config_flags%dust_opt)
      CASE (DUSTGOCART)
        CALL wrf_debug(15,'Gocart dust emissions')
        call gocart_dust_driver(ktau,dtstep,config_flags,julday,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,dustin,           &
         ivgtyp,isltyp,vegfra,xland,xlat,xlong,gsw,dx,g,emis_dust,        &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
      CASE (DUSTGOCARTAFWA)
        CALL wrf_debug(15,'AFWA modified Gocart dust emissions')
        call gocart_dust_afwa_driver(ktau,dtstep,config_flags,julday,alt,t_phy,moist,u_phy,  &
         v_phy,chem,rho_phy,dz8w,smois,u10,v10,p8w,erod,erod_dri,dustin,snowh,zs,   &
         ivgtyp,isltyp,vegfra,lai_vegmask,xland,xlat,xlong,gsw,dx,g,emis_dust,      &
         ust,znt,clayfrac,sandfrac,clayfrac_nga,sandfrac_nga,afwa_dustloft,         &
         tot_dust,tot_edust,vis_dust,dust_alpha,dust_gamma,dust_smtune,dust_ustune, &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
      CASE (DUSTUOC)
       CALL wrf_debug(15,'UoC dust emission schemes')

       scheme_select:  SELECT CASE(config_flags%dust_schme)
       CASE (SHAO_2001)
        imod = 1
       CASE (SHAO_2004)
        imod = 2
       CASE (SHAO_2011)
        imod = 3
       CASE DEFAULT
        imod = 2
       END SELECT scheme_select
       call uoc_dust_driver (ktau,dtstep,config_flags,                     &
         chem,rho_phy,dz8w,smois,ust, isltyp,vegfra,g,emis_dust,           &
         ust_t, imod, rough_cor, smois_cor, soil_top_cat, erod,            &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )             
      CASE DEFAULT 
        if(dust_emiss_active.eq.1) then
            CALL wrf_debug(15,'MOSAIC or SORGAM dust emissions')
        else
             CALL wrf_debug(15,'no dust emissions')
        end if
      END SELECT dust_select

      dms_select:  SELECT CASE(config_flags%dmsemis_opt)
      CASE (DMSGOCART)
        CALL wrf_debug(15,'Gocart dms emissions')
        call gocart_dmsemis(dtstep,config_flags,alt,t_phy,u_phy,  &
         v_phy,chem,rho_phy,dz8w,u10,v10,p8w,dms_0,tsk,                  &
         ivgtyp,isltyp,xland,dx,g, &
         ids,ide, jds,jde, kds,kde,                                        &
         ims,ime, jms,jme, kms,kme,                                        &
         its,ite, jts,jte, kts,kte                                         )
      CASE DEFAULT 
        CALL wrf_debug(15,'no dms emissions')
      END SELECT dms_select

    ksub=0



    fire_select:  SELECT CASE(config_flags%biomass_burn_opt)
     CASE (BIOMASSB,BIOMASSB_MOZC,BIOMASSB_MOZ,BIOMASSB_T1_MOZCART,BIOMASSB_GHG,BIOMASSB_SAPRC)















       CALL wrf_debug(15,'fire emissions: adding biomassb emissions')
       call add_emis_burn(id,dtstep,ktau,dz8w,rho_phy,chem,&
            julday,gmt,xlat,xlong,t_phy,p_phy,                           &
            ebu,config_flags%chem_opt,0,config_flags%biomass_burn_opt,     &
            num_chem,ids,ide, jds,jde, kds,kde,                                   &
            ims,ime, jms,jme, kms,kme,                                   &
            its,ite, jts,jte, kts,kte                                    )
     CASE DEFAULT 
       CALL wrf_debug(15,'no biomass burning')
    END SELECT fire_select



    bioem_select: SELECT CASE(config_flags%bio_emiss_opt)
     CASE (GUNTHER1)
       if(ktau.eq.1.or.mod(ktau,stepbioe).eq.0)then
          CALL wrf_debug(15,'biogenic emissions: calling Gunther1')
          call bio_emissions(id,ktau,dtstep,DX,config_flags,               &
               gmt,julday,t_phy,moist,p8w,t8w,                             &
               e_bio,p_phy,chem,rho_phy,dz8w,ne_area,                      &
               ivgtyp,gsw,vegfra,rmol,ust,znt,xlat,xlong,z_at_w,           &
               ids,ide, jds,jde, kds,kde,                                  &
               ims,ime, jms,jme, kms,kme,                                  &
               its,ite, jts,jte, kts,kte                                   )
       endif
     CASE (BEIS314)
       if( do_bioemiss ) then
         beis314_check_mechanism_ok: SELECT CASE(config_flags%chem_opt) 
            CASE (RADM2, RADM2_KPP, RADM2SORG, RADM2SORG_AQ, RADM2SORG_AQCHEM, RADM2SORG_KPP, &
                  RACM_KPP, RACMPM_KPP, RACM_MIM_KPP, RACMSORG_AQ,RACMSORG_AQCHEM_KPP,        &
                  RACM_ESRLSORG_AQCHEM_KPP, RACMSORG_KPP,RACM_ESRLSORG_KPP, RACM_SOA_VBS_KPP, &
                  RACM_SOA_VBS_AQCHEM_KPP, RACM_SOA_VBS_HET_KPP, CBM4_KPP, NMHC9_KPP, GOCARTRACM_KPP,GOCARTRADM2)
            CASE DEFAULT 
               CALL wrf_error_fatal3("<stdin>",845,&
                  "emissions_driver: beis3.1.4 biogenic emis. implemented for RADM2 & RACM only")
         END SELECT beis314_check_mechanism_ok
         CALL wrf_debug(15,'biogenic emissions: calling beis3.1.4')
         call bio_emissions_beis314(id,config_flags,ktau,curr_secs,    &
               dtstep,julday,gmt,xlat,xlong,t_phy,p_phy,gsw,           &
               sebio_iso,sebio_oli,sebio_api,sebio_lim,sebio_xyl,      &
               sebio_hc3,sebio_ete,sebio_olt,sebio_ket,sebio_ald,      &
               sebio_hcho,sebio_eth,sebio_ora2,sebio_co,sebio_nr,      &
               sebio_sesq,sebio_mbo,                                   &
               noag_grow,noag_nongrow,nononag,slai,                    &
               ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,           &
               ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,           &
               ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,   &
               ebio_sesq,ebio_mbo,                                     &
               ids,ide, jds,jde, kds,kde,                              &
               ims,ime, jms,jme, kms,kme,                              &
               its,ite, jts,jte, kts,kte                               )
       endif

     CASE (MEGAN2)
       if(ktau.eq.1.or.mod(ktau,stepbioe).eq.0)then                        
         CALL wrf_debug(15,'biogenic emissions: calling megan v2.04')  
         call bio_emissions_megan2(id,config_flags,ktau,dtstep,        &
               curr_secs,julday,gmt,xlat,xlong,p_phy,rho_phy,dz8w,     &
               chem,ne_area,                                           &
               current_month,                                          &
               T2,swdown,                                              &
               nmegan, EFmegan, msebio_isop,                           &
               mlai,                                                   &
               pftp_bt, pftp_nt, pftp_sb, pftp_hb,                     &
               mtsa,                                                   &
               mswdown,                                                &
               mebio_isop, mebio_apin, mebio_bpin, mebio_bcar,         &
               mebio_acet, mebio_mbo, mebio_no,                        &
               ebio_iso,ebio_oli,ebio_api,ebio_lim,                    &
               ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,           &
               ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_no,           &
               ebio_c10h16,ebio_tol,ebio_bigalk, ebio_ch3oh,ebio_acet,         &
               ebio_nh3,ebio_no2,ebio_c2h5oh,ebio_ch3cooh,ebio_mek,            &
               ebio_bigene,ebio_c2h6,ebio_c2h4,ebio_c3h6,ebio_c3h8,ebio_so2,   &
               ebio_dms,ebio_hcn,                                              &
               ebio_c5h8,ebio_apinene,ebio_bpinene,ebio_toluene,       &
               ebio_ch3cho,ebio_ch3co2h,ebio_tbut2ene,ebio_c2h5cho,    &
               ebio_nc4h10, &
               ebio_sesq, ebio_mbo,ebio_bpi,ebio_myrc,                 &
               ebio_alk3, ebio_alk4, ebio_alk5, ebio_ole1, ebio_ole2,    &    
               ebio_aro1, ebio_aro2, ebio_ccho, ebio_meoh,               &    
               ebio_ethene, ebio_hcooh, ebio_terp, ebio_bald,            &    
               ebio_cco_oh, ebio_rco_oh,                                 &
               e_bio,                                                  &
               ids,ide, jds,jde, kds,kde,                              &
               ims,ime, jms,jme, kms,kme,                              &
               its,ite, jts,jte, kts,kte                               )
       endif

     CASE (MEGAN2_CLM)


       convert2(its:ite,jts:jte) = rho_phy(its:ite,kts,jts:jte)*60./.02897
       print*,'jdf shr_megan_mechcomps_n',shr_megan_mechcomps_n
       print*,'jdf index iso     ',p_isoprene
       print*,'jdf index no      ',p_no
       print*,'jdf index no2     ',p_no2
       print*,'jdf index co      ',p_co
       print*,'jdf index hcho    ',p_hcho
       print*,'jdf index ald     ',p_ald
       print*,'jdf index acet    ',p_acet
       print*,'jdf index tol     ',p_tol
       print*,'jdf index c10h16  ',p_c10h16
       print*,'jdf index so2     ',p_so2
       print*,'jdf index dms     ',p_dms
       print*,'jdf index bigalk  ',p_bigalk
       print*,'jdf index bigene  ',p_bigene
       print*,'jdf index nh3     ',p_nh3
       print*,'jdf index c3hoh   ',p_ch3oh
       print*,'jdf index c2h5oh  ',p_c2h5oh
       print*,'jdf index ch3co2h ',p_ch3co2h
       print*,'jdf index mek     ',p_mek
       print*,'jdf index c2h4    ',p_c2h4
       print*,'jdf index ceh6    ',p_c3h6
       print*,'jdf index c3h8    ',p_c3h8
       do m = 1,shr_megan_mechcomps_n
         p_in_chem = shr_megan_mechcomps(m)%index
         print*,'jdf index',m,p_in_chem

         IF ( p_in_chem+1 == p_isoprene ) THEN
           ebio_iso(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_no ) THEN
           ebio_no(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_no2 ) THEN
           ebio_no2(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_co ) THEN
           ebio_co(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_hcho ) THEN
           ebio_hcho(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_ald ) THEN
           ebio_ald(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_acet ) THEN
           ebio_acet(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_tol ) THEN
           ebio_tol(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_c10h16 ) THEN
           ebio_c10h16(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_so2 ) THEN
           ebio_so2(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_dms ) THEN
           ebio_dms(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_bigalk ) THEN
           ebio_bigalk(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_bigene ) THEN
           ebio_bigene(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_nh3 ) THEN
           ebio_nh3(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_ch3oh ) THEN
           ebio_ch3oh(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_c2h5oh ) THEN
           ebio_c2h5oh(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_ch3co2h ) THEN
           ebio_ch3cooh(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_mek ) THEN
           ebio_mek(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_c2h4 ) THEN
           ebio_c2h4(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_c2h6 ) THEN
           ebio_c2h6(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_c3h6 ) THEN
           ebio_c3h6(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ELSEIF ( p_in_chem+1 == p_c3h8 ) THEN
           ebio_c3h8(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
         ENDIF

         
         IF ( p_in_chem+1 == p_terp ) THEN
           e_bio(its:ite,jts:jte,p_in_chem) = e_bio(its:ite,jts:jte,p_in_chem)*0.5
           ebio_terp(its:ite,jts:jte) = e_bio(its:ite,jts:jte,p_in_chem)*convert2(its:ite,jts:jte)
           print*,"Reduce terp emissions by half!"
         ENDIF

       end do


     CASE DEFAULT 
       if( do_bioemiss ) &
            e_bio(its:ite,jts:jte,1:ne_area) = 0.


                                                     
    END SELECT bioem_select



    gas_addemiss_select: SELECT CASE(config_flags%chem_opt)
    CASE (RADM2, RADM2_KPP, RADM2SORG, RADM2SORG_AQ, RADM2SORG_AQCHEM, RADM2SORG_KPP, &
          RACM_KPP, RACMPM_KPP, RACM_MIM_KPP, RACMSORG_AQ, RACMSORG_AQCHEM_KPP, RACM_ESRLSORG_AQCHEM_KPP, RACMSORG_KPP, &
          RACM_SOA_VBS_KPP, RACM_SOA_VBS_AQCHEM_KPP, RACM_SOA_VBS_HET_KPP, RACM_ESRLSORG_KPP, MOZART_KPP, MOZCART_KPP,  &
          T1_MOZCART_KPP, MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP, &
          CRIMECH_KPP, CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP )
       IF(config_flags%emiss_inpt_opt /= 3 ) then
       IF(config_flags%kemit .GT. kte-ksub) THEN
         k=config_flags%kemit
         write(message,'(" WARNING: EMISSIONS_DRIVER: KEMIT > KTE ",3i6)') kme,kte-ksub,k
         CALL WRF_MESSAGE (message)
       ENDIF
       call wrf_debug(15,'emissions_driver calling add_anthropogenics')
       call add_anthropogenics(id,dtstep,dz8w,config_flags,rho_phy,alt, &
            chem, emis_ant,emis_aircraft,                               &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )
       call wrf_debug(15,'emissions_driver calling add_biogenics')
       call add_biogenics(id,dtstep,dz8w,config_flags, rho_phy,chem,    &
            e_bio,ne_area,                                              &
            ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,               &
            ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,               &
            ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,       &
            ebio_sesq,ebio_mbo,                                         & 
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )

       end if 




    CASE (CBMZ, CBMZ_BB, CBMZ_BB_KPP, CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_8BIN, &
          CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, & 
          CBMZSORG, CBMZSORG_AQ, CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, &
          CBMZ_MOSAIC_KPP, &
          CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &
          CBMZ_CAM_MAM3_NOAQ, CBMZ_CAM_MAM3_AQ, CBMZ_CAM_MAM7_NOAQ, CBMZ_CAM_MAM7_AQ)
       IF(config_flags%kemit .GT. kte-ksub) THEN
          message = ' EMISSIONS_DRIVER: KEMIT > KME '
          CALL wrf_error_fatal3("<stdin>",1039,&
message)
       ENDIF
       call wrf_debug(15,'emissions_driver calling cbmz_addemiss_anthro')
       call cbmz_addemiss_anthro( id, dtstep, dz8w, config_flags,        &
            rho_phy, chem,                                               &
            emis_ant,alt,ids,ide, jds,jde, kds,kde,                      &
            ims,ime, jms,jme, kms,kme,                                   &
            its,ite, jts,jte, kts,kte                                    )
       call wrf_debug(15,'emissions_driver calling cbmz_addemiss_bio')
       
       
       
       
       
       
       
       
       
       
       call add_biogenics(id,dtstep,dz8w,config_flags, rho_phy,chem,    &
            e_bio,ne_area,                                              &
            ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,               &
            ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,               &
            ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,       &
            ebio_sesq,ebio_mbo,                                         & 
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )

    CASE (CB05_SORG_AQ_KPP, CB05_SORG_VBS_AQ_KPP)
       IF(config_flags%kemit .GT. kte-ksub) THEN
         message = ' EMISSIONS_DRIVER: KEMIT > KME '
         CALL wrf_error_fatal3("<stdin>",1072,&
message)
       ENDIF
       call wrf_debug(15,'emissions_driver calling cb05_addemiss_anthro')
       call cb05_addemiss_anthro( id, dtstep, dz8w, config_flags,        &
            rho_phy, chem,                                               &
            emis_ant,ids,ide, jds,jde, kds,kde,                                   &
            ims,ime, jms,jme, kms,kme,                                   &
            its,ite, jts,jte, kts,kte                                    )
       call wrf_debug(15,'emissions_driver calling cb05_addemiss_bio')
       
       if (config_flags%bio_emiss_opt .ne. 0 .and.                      &
           config_flags%bio_emiss_opt .ne. GUNTHER1) then
         call add_biogenics(id,dtstep,dz8w,config_flags, rho_phy,chem,    &
            e_bio,ne_area,                                              &
            ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,               &
            ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,               &
            ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,       &
            ebio_sesq,ebio_mbo,                                         &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )
       endif

       if ( config_flags%bio_emiss_opt .eq. GUNTHER1 ) then
         call cb05_addemiss_bio( id, dtstep, dz8w, config_flags,         &
              rho_phy, chem, e_bio, ne_area, emis_ant(ims,kms,jms,p_e_iso),&
              ids,ide, jds,jde, kds,kde,                                 &
              ims,ime, jms,jme, kms,kme,                                 &
              its,ite, jts,jte, kts,kte                                  )
       endif

    CASE (CHEM_TRACER)
       do j=jts,jte  
          do i=its,ite  
             do k=kts,min(config_flags%kemit,kte-ksub)
                conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                   &
                     +emis_ant(i,k,j,p_e_so2)*conv
                chem(i,k,j,p_co)  = chem(i,k,j,p_co)                     &
                     +emis_ant(i,k,j,p_e_co)*conv
                chem(i,k,j,p_no)  = chem(i,k,j,p_no)                     &
                     +emis_ant(i,k,j,p_e_co)*conv
                chem(i,k,j,p_ald)  = chem(i,k,j,p_ald)                   &
                     +emis_ant(i,k,j,p_e_co)*conv
                chem(i,k,j,p_hcho)  = chem(i,k,j,p_hcho)                 &
                     +emis_ant(i,k,j,p_e_co)*conv
                chem(i,k,j,p_ora2)  = chem(i,k,j,p_ora2)                 &
                     +emis_ant(i,k,j,p_e_co)*conv
             end do
          end do
       end do

    CASE(CO2_TRACER,GHG_TRACER)  

      
   CALL VPRM(            ids,ide, jds,jde,                   &
                         ims,ime, jms,jme,                   &
                         its,ite, jts,jte,                   &

                         vprm_in,rad_vprm,lambda_vprm,       &
                         alpha_vprm,resp_vprm,               &
                         T2,swdown,                          &
                         eghg_bio                            )

  
   if (p_ch4_bio .GT. 1) then

   CALL KAPLAN(          ids,ide, jds,jde,                                        &
                         ims,ime, jms,jme,                                        &
                         its,ite, jts,jte,                                        &

                         xtime, tslb, smois, wet_in,                              &
                         isltyp,tsk,eghg_bio,                                     &
                         config_flags%num_soil_layers,config_flags%wpeat,         &
                         config_flags%wflood                                      )

   CALL SOILUPTAKE(      ids,ide, jds,jde,                                        &
                         ims,ime, jms,jme,                                        &
                         its,ite, jts,jte,                                        &

                         smois, isltyp, eghg_bio,                                 &
                         rainc, rainnc,                                           &
                         potevp, sfcevp, lu_index, T2, xtime,                     &
                         config_flags%num_soil_layers, wet_in                     )

   CALL termite(         ids,ide, jds,jde,                                        &
                         ims,ime, jms,jme,                                        &
                         its,ite, jts,jte,                                        &

                         xtime,eghg_bio,ivgtyp,                                   &
                         biomt_par,emit_par                                       )

   end if

   
   
   IF (config_flags%emiss_inpt_opt==16) THEN
      CALL add_ghg_fluxes(  ids,ide, jds,jde, kds,kde,          &
                            ims,ime, jms,jme, kms,kme,          &
                            its,ite, jts,jte, kts,kte,          &

                            dtstep,dz8w,config_flags,rho_phy,   &
                            chem,emis_ant,eghg_bio,ebio_co2oce  )
   END IF

    CASE (SAPRC99_KPP) 
     if(config_flags%emiss_opt == 13 ) then
       do j=jts,jte
          do i=its,ite
             do k=kts,min(config_flags%kemit,kte-ksub)
                conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                             &
                     +emis_ant(i,k,j,p_e_so2)*conv
                chem(i,k,j,p_c2h6)  = chem(i,k,j,p_c2h6)                           &
                     +emis_ant(i,k,j,p_e_c2h6)*conv
                chem(i,k,j,p_c3h8)  = chem(i,k,j,p_c3h8)                           &
                     +emis_ant(i,k,j,p_e_c3h8)*conv
                chem(i,k,j,p_c2h2)  = chem(i,k,j,p_c2h2)                           &
                     +emis_ant(i,k,j,p_e_c2h2)*conv
                chem(i,k,j,p_alk3)  = chem(i,k,j,p_alk3)                           &
                     +emis_ant(i,k,j,p_e_alk3)*conv
                chem(i,k,j,p_alk4)  = chem(i,k,j,p_alk4)                           &
                     +emis_ant(i,k,j,p_e_alk4)*conv
                chem(i,k,j,p_alk5)  = chem(i,k,j,p_alk5)                           &
                     +emis_ant(i,k,j,p_e_alk5)*conv
                chem(i,k,j,p_ethene)  = chem(i,k,j,p_ethene)                       &
                     +emis_ant(i,k,j,p_e_ethene)*conv
                chem(i,k,j,p_c3h6)  = chem(i,k,j,p_c3h6)                           &
                     +emis_ant(i,k,j,p_e_c3h6)*conv
                chem(i,k,j,p_ole1)  = chem(i,k,j,p_ole1)                           &
                     +emis_ant(i,k,j,p_e_ole1)*conv
                chem(i,k,j,p_ole2)  = chem(i,k,j,p_ole2)                           &
                     +emis_ant(i,k,j,p_e_ole2)*conv
                chem(i,k,j,p_aro1)  = chem(i,k,j,p_aro1)                           &
                     +emis_ant(i,k,j,p_e_aro1)*conv
                chem(i,k,j,p_aro2)  = chem(i,k,j,p_aro2)                           &
                     +emis_ant(i,k,j,p_e_aro2)*conv
                chem(i,k,j,p_hcho)  = chem(i,k,j,p_hcho)                           &
                     +emis_ant(i,k,j,p_e_hcho)*conv
                chem(i,k,j,p_ccho)  = chem(i,k,j,p_ccho)                           &
                     +emis_ant(i,k,j,p_e_ccho)*conv
                chem(i,k,j,p_rcho)  = chem(i,k,j,p_rcho)                           &
                     +emis_ant(i,k,j,p_e_rcho)*conv
                chem(i,k,j,p_acet)  = chem(i,k,j,p_acet)                           &
                     +emis_ant(i,k,j,p_e_acet)*conv
                chem(i,k,j,p_mek)  = chem(i,k,j,p_mek)                             &
                     +emis_ant(i,k,j,p_e_mek)*conv
                chem(i,k,j,p_isoprene)  = chem(i,k,j,p_isoprene)                   &
                     +emis_ant(i,k,j,p_e_isoprene)*conv
                chem(i,k,j,p_terp)  = chem(i,k,j,p_terp)                           &
                     +emis_ant(i,k,j,p_e_terp)*conv
                chem(i,k,j,p_sesq)  = chem(i,k,j,p_sesq)                           &
                     +emis_ant(i,k,j,p_e_sesq)*conv
                chem(i,k,j,p_co)  = chem(i,k,j,p_co)                               &
                     +emis_ant(i,k,j,p_e_co)*conv
                chem(i,k,j,p_no)  = chem(i,k,j,p_no)                               &
                     +emis_ant(i,k,j,p_e_no)*conv
                chem(i,k,j,p_no2)  = chem(i,k,j,p_no2)                             &
                     +emis_ant(i,k,j,p_e_no2)*conv
                chem(i,k,j,p_phen)  = chem(i,k,j,p_phen)                           &
                     +emis_ant(i,k,j,p_e_phen)*conv
                chem(i,k,j,p_cres)  = chem(i,k,j,p_cres)                           &
                     +emis_ant(i,k,j,p_e_cres)*conv
                chem(i,k,j,p_meoh)  = chem(i,k,j,p_meoh)                           &
                     +emis_ant(i,k,j,p_e_meoh)*conv
                chem(i,k,j,p_gly)  = chem(i,k,j,p_gly)                             &
                     +emis_ant(i,k,j,p_e_gly)*conv
                chem(i,k,j,p_mgly)  = chem(i,k,j,p_mgly)                           &
                     +emis_ant(i,k,j,p_e_mgly)*conv
                chem(i,k,j,p_bacl)  = chem(i,k,j,p_bacl)                           &
                     +emis_ant(i,k,j,p_e_bacl)*conv
                chem(i,k,j,p_isoprod)  = chem(i,k,j,p_isoprod)                     &
                     +emis_ant(i,k,j,p_e_isoprod)*conv
                chem(i,k,j,p_methacro)  = chem(i,k,j,p_methacro)                   &
                     +emis_ant(i,k,j,p_e_methacro)*conv
                chem(i,k,j,p_mvk)  = chem(i,k,j,p_mvk)                             &
                     +emis_ant(i,k,j,p_e_mvk)*conv
                chem(i,k,j,p_prod2)  = chem(i,k,j,p_prod2)                         &
                     +emis_ant(i,k,j,p_e_prod2)*conv
                chem(i,k,j,p_ch4)  = chem(i,k,j,p_ch4)                             &
                     +emis_ant(i,k,j,p_e_ch4)*conv
                chem(i,k,j,p_bald)  = chem(i,k,j,p_bald)                           &
                     +emis_ant(i,k,j,p_e_bald)*conv
                chem(i,k,j,p_hcooh)  = chem(i,k,j,p_hcooh)                         &
                     +emis_ant(i,k,j,p_e_hcooh)*conv
                chem(i,k,j,p_cco_oh)  = chem(i,k,j,p_cco_oh)                       &
                     +emis_ant(i,k,j,p_e_cco_oh)*conv
                chem(i,k,j,p_rco_oh)  = chem(i,k,j,p_rco_oh)                       &
                     +emis_ant(i,k,j,p_e_rco_oh)*conv

             end do
          end do
       end do
      else
       do j=jts,jte
          do i=its,ite
             do k=kts,min(config_flags%kemit,kte-ksub)
                conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                   &
                     +emis_ant(i,k,j,p_e_so2)*conv
                chem(i,k,j,p_co)  = chem(i,k,j,p_co)                     &
                     +emis_ant(i,k,j,p_e_co)*conv
                chem(i,k,j,p_no)  = chem(i,k,j,p_no)                     &
                     +emis_ant(i,k,j,p_e_no)*conv
                chem(i,k,j,p_hcho)  = chem(i,k,j,p_hcho)                 &
                     +emis_ant(i,k,j,p_e_hcho)*conv
             end do
          end do
       end do
      endif

      
      call add_biogenics(id,dtstep,dz8w,config_flags, rho_phy,chem,     &
            e_bio,ne_area,                                              &
            ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,               &
            ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,               &
            ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,       &
            ebio_sesq,ebio_mbo,                                         & 
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )
      
    CASE (SAPRC99_MOSAIC_4BIN_VBS2_KPP, SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, &
         SAPRC99_MOSAIC_8BIN_VBS2_KPP, & 
         SAPRC99_MOSAIC_20BIN_VBS2_AQ_KPP,SAPRC99_MOSAIC_20BIN_VBS2_KPP, & 
         SAPRC99_MOSAIC_12BIN_VBS2_AQ_KPP,SAPRC99_MOSAIC_12BIN_VBS2_KPP)
     if(config_flags%emiss_opt == 13 ) then
       do j=jts,jte
          do i=its,ite
             do k=kts,min(config_flags%kemit,kte-ksub)
                conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                conv3 = (dtstep/dz8w(i,k,j))*alt(i,k,j)*28/250*1e-3
                conv4= (dtstep/dz8w(i,k,j))*alt(i,k,j)*28/226*1e-3  
                oconv3= (dtstep/dz8w(i,k,j))*alt(i,k,j)*28/283*1e-3*4.5 
                oconv4=(dtstep/dz8w(i,k,j))*alt(i,k,j)*28/226*1e-3*0.9  
                chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                             &
                     +emis_ant(i,k,j,p_e_so2)*conv
                chem(i,k,j,p_c2h6)  = chem(i,k,j,p_c2h6)                           &
                     +emis_ant(i,k,j,p_e_c2h6)*conv
                chem(i,k,j,p_c3h8)  = chem(i,k,j,p_c3h8)                           &
                     +emis_ant(i,k,j,p_e_c3h8)*conv
                chem(i,k,j,p_c2h2)  = chem(i,k,j,p_c2h2)                           &
                     +emis_ant(i,k,j,p_e_c2h2)*conv
                chem(i,k,j,p_alk3)  = chem(i,k,j,p_alk3)                           &
                     +emis_ant(i,k,j,p_e_alk3)*conv
                chem(i,k,j,p_alk4)  = chem(i,k,j,p_alk4)                           &
                     +emis_ant(i,k,j,p_e_alk4)*conv
                chem(i,k,j,p_alk5)  = chem(i,k,j,p_alk5)                           &
                     +emis_ant(i,k,j,p_e_alk5)*conv
                chem(i,k,j,p_ethene)  = chem(i,k,j,p_ethene)                       &
                     +emis_ant(i,k,j,p_e_ethene)*conv
                chem(i,k,j,p_c3h6)  = chem(i,k,j,p_c3h6)                           &
                     +emis_ant(i,k,j,p_e_c3h6)*conv
                chem(i,k,j,p_ole1)  = chem(i,k,j,p_ole1)                           &
                     +emis_ant(i,k,j,p_e_ole1)*conv
                chem(i,k,j,p_ole2)  = chem(i,k,j,p_ole2)                           &
                     +emis_ant(i,k,j,p_e_ole2)*conv
                chem(i,k,j,p_aro1)  = chem(i,k,j,p_aro1)                           &
                     +emis_ant(i,k,j,p_e_aro1)*conv
                chem(i,k,j,p_aro2)  = chem(i,k,j,p_aro2)                           &
                     +emis_ant(i,k,j,p_e_aro2)*conv
                chem(i,k,j,p_hcho)  = chem(i,k,j,p_hcho)                           &
                     +emis_ant(i,k,j,p_e_hcho)*conv
                chem(i,k,j,p_ccho)  = chem(i,k,j,p_ccho)                           &
                     +emis_ant(i,k,j,p_e_ccho)*conv
                chem(i,k,j,p_rcho)  = chem(i,k,j,p_rcho)                           &
                     +emis_ant(i,k,j,p_e_rcho)*conv
                chem(i,k,j,p_acet)  = chem(i,k,j,p_acet)                           &
                     +emis_ant(i,k,j,p_e_acet)*conv
                chem(i,k,j,p_mek)  = chem(i,k,j,p_mek)                             &
                     +emis_ant(i,k,j,p_e_mek)*conv
                chem(i,k,j,p_isoprene)  = chem(i,k,j,p_isoprene)                   &
                     +emis_ant(i,k,j,p_e_isoprene)*conv
                chem(i,k,j,p_terp)  = chem(i,k,j,p_terp)                           &
                     +emis_ant(i,k,j,p_e_terp)*conv
                chem(i,k,j,p_sesq)  = chem(i,k,j,p_sesq)                           &
                     +emis_ant(i,k,j,p_e_sesq)*conv
                chem(i,k,j,p_co)  = chem(i,k,j,p_co)                               &
                     +emis_ant(i,k,j,p_e_co)*conv
                chem(i,k,j,p_no)  = chem(i,k,j,p_no)                               &
                     +emis_ant(i,k,j,p_e_no)*conv




              if((k.eq.kts).and.xland(i,j).eq.1.0)        & 
              chem(i,k,j,p_no)  = chem(i,k,j,p_no) +15.0/30.01*conv 

              if((k.eq.kts).and.xland(i,j).eq.1.0)        & 
              chem(i,k,j,p_nh3)  = chem(i,k,j,p_nh3) +15.0/17.0*conv

              if((k.eq.kts).and.xland(i,j).eq.1.0)        & 
              chem(i,k,j,p_dms)  = chem(i,k,j,p_dms) +0.24*conv  


                chem(i,k,j,p_no2)  = chem(i,k,j,p_no2)                             &
                     +emis_ant(i,k,j,p_e_no2)*conv
                chem(i,k,j,p_phen)  = chem(i,k,j,p_phen)                           &
                     +emis_ant(i,k,j,p_e_phen)*conv
                chem(i,k,j,p_cres)  = chem(i,k,j,p_cres)                           &
                     +emis_ant(i,k,j,p_e_cres)*conv
                chem(i,k,j,p_meoh)  = chem(i,k,j,p_meoh)                           &
                     +emis_ant(i,k,j,p_e_meoh)*conv
                chem(i,k,j,p_gly)  = chem(i,k,j,p_gly)                             &
                     +emis_ant(i,k,j,p_e_gly)*conv
                chem(i,k,j,p_mgly)  = chem(i,k,j,p_mgly)                           &
                     +emis_ant(i,k,j,p_e_mgly)*conv
                chem(i,k,j,p_bacl)  = chem(i,k,j,p_bacl)                           &
                     +emis_ant(i,k,j,p_e_bacl)*conv
                chem(i,k,j,p_isoprod)  = chem(i,k,j,p_isoprod)                     &
                     +emis_ant(i,k,j,p_e_isoprod)*conv
                chem(i,k,j,p_methacro)  = chem(i,k,j,p_methacro)                   &
                     +emis_ant(i,k,j,p_e_methacro)*conv
                chem(i,k,j,p_mvk)  = chem(i,k,j,p_mvk)                             &
                     +emis_ant(i,k,j,p_e_mvk)*conv
                chem(i,k,j,p_prod2)  = chem(i,k,j,p_prod2)                         &
                     +emis_ant(i,k,j,p_e_prod2)*conv
                chem(i,k,j,p_ch4)  = chem(i,k,j,p_ch4)                             &
                     +emis_ant(i,k,j,p_e_ch4)*conv
                chem(i,k,j,p_bald)  = chem(i,k,j,p_bald)                           &
                     +emis_ant(i,k,j,p_e_bald)*conv
                chem(i,k,j,p_hcooh)  = chem(i,k,j,p_hcooh)                         &
                     +emis_ant(i,k,j,p_e_hcooh)*conv
                chem(i,k,j,p_cco_oh)  = chem(i,k,j,p_cco_oh)                       &
                     +emis_ant(i,k,j,p_e_cco_oh)*conv
                chem(i,k,j,p_rco_oh)  = chem(i,k,j,p_rco_oh)                       &
                     +emis_ant(i,k,j,p_e_rco_oh)*conv
                chem(i,k,j,p_nh3)  = chem(i,k,j,p_nh3)                       &
                     +emis_ant(i,k,j,p_e_nh3)*conv





         chem(i,k,j,p_pcg1_b_c)  =  chem(i,k,j,p_pcg1_b_c)                        &
        +(emis_ant(i,k,j,p_e_orgi_bb)/1.57+emis_ant(i,k,j,p_e_orgj_bb)/1.57)*conv3*1.17
        chem(i,k,j,p_pcg2_b_c)  =  chem(i,k,j,p_pcg2_b_c)                        &
        +(emis_ant(i,k,j,p_e_orgi_bb)/1.57+emis_ant(i,k,j,p_e_orgj_bb)/1.57)*conv3*7.605
        chem(i,k,j,p_pcg1_f_c)  =  chem(i,k,j,p_pcg1_f_c)                        &
        +(emis_ant(i,k,j,p_e_orgi_a)/1.25+emis_ant(i,k,j,p_e_orgj_a)/1.25)*conv3*1.17
        chem(i,k,j,p_pcg2_f_c)  =  chem(i,k,j,p_pcg2_f_c)                        &
        +(emis_ant(i,k,j,p_e_orgi_a)/1.25+emis_ant(i,k,j,p_e_orgj_a)/1.25)*conv3*7.605
        chem(i,k,j,p_pcg1_b_o)  =  chem(i,k,j,p_pcg1_b_o)                        &
        +(emis_ant(i,k,j,p_e_orgi_bb)/1.57+emis_ant(i,k,j,p_e_orgj_bb)/1.57)*conv3*0.40
        chem(i,k,j,p_pcg2_b_o)  =  chem(i,k,j,p_pcg2_b_o)                        &
        +(emis_ant(i,k,j,p_e_orgi_bb)/1.57+emis_ant(i,k,j,p_e_orgj_bb)/1.57)*conv3*2.60
        chem(i,k,j,p_pcg1_f_o)  =  chem(i,k,j,p_pcg1_f_o)                        &
        +(emis_ant(i,k,j,p_e_orgi_a)/1.25+emis_ant(i,k,j,p_e_orgj_a)/1.25)*conv3*0.08
        chem(i,k,j,p_pcg2_f_o)  =  chem(i,k,j,p_pcg2_f_o)                        &
        +(emis_ant(i,k,j,p_e_orgi_a)/1.25+emis_ant(i,k,j,p_e_orgj_a)/1.25)*conv3*0.52


             end do
          end do
       end do
      else
       do j=jts,jte
          do i=its,ite
             do k=kts,min(config_flags%kemit,kte-ksub)
                conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                   &
                     +emis_ant(i,k,j,p_e_so2)*conv
                chem(i,k,j,p_co)  = chem(i,k,j,p_co)                     &
                     +emis_ant(i,k,j,p_e_co)*conv
                chem(i,k,j,p_no)  = chem(i,k,j,p_no)                     &
                     +emis_ant(i,k,j,p_e_no)*conv
                chem(i,k,j,p_hcho)  = chem(i,k,j,p_hcho)                 &
                     +emis_ant(i,k,j,p_e_hcho)*conv
             end do
          end do
       end do
      endif
      
      call add_biogenics(id,dtstep,dz8w,config_flags, rho_phy,chem,     &
            e_bio,ne_area,                                              &
            ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,               &
            ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,               &
            ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,       &
            ebio_sesq,ebio_mbo,                                         & 
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )


    CASE (SAPRC99_MOSAIC_4BIN_IEPOX_VBS_AQ_KPP,SAPRC99_SIMPLESOM_MOSAIC_4BIN_AQ_KPP) 
     if(config_flags%emiss_opt == 13 ) then


        



       do j=jts,jte
          do i=its,ite
             do k=kts,min(config_flags%kemit,kte-ksub)


                

                conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                conv3 = (dtstep/dz8w(i,k,j))*alt(i,k,j)*28/250*1e-3
                conv4= (dtstep/dz8w(i,k,j))*alt(i,k,j)*28/226*1e-3  
                oconv3= (dtstep/dz8w(i,k,j))*alt(i,k,j)*28/283*1e-3*4.5 
                oconv4=(dtstep/dz8w(i,k,j))*alt(i,k,j)*28/226*1e-3*0.9  


                chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                             &
                     +emis_ant(i,k,j,p_e_so2)*conv*config_flags%uq_so2
                chem(i,k,j,p_c2h6)  = chem(i,k,j,p_c2h6)                           &
                     +emis_ant(i,k,j,p_e_c2h6)*conv
                chem(i,k,j,p_c3h8)  = chem(i,k,j,p_c3h8)                           &
                     +emis_ant(i,k,j,p_e_c3h8)*conv
                chem(i,k,j,p_c2h2)  = chem(i,k,j,p_c2h2)                           &
                     +emis_ant(i,k,j,p_e_c2h2)*conv
                chem(i,k,j,p_alk3)  = chem(i,k,j,p_alk3)                           &
                     +emis_ant(i,k,j,p_e_alk3)*conv
                chem(i,k,j,p_alk4)  = chem(i,k,j,p_alk4)                           &
                     +emis_ant(i,k,j,p_e_alk4)*conv


                


                chem(i,k,j,p_alk5)  = chem(i,k,j,p_alk5)                           &
                     +emis_ant(i,k,j,p_e_alk5)*conv
                chem(i,k,j,p_ethene)  = chem(i,k,j,p_ethene)                       &
                     +emis_ant(i,k,j,p_e_ethene)*conv
                chem(i,k,j,p_c3h6)  = chem(i,k,j,p_c3h6)                           &
                     +emis_ant(i,k,j,p_e_c3h6)*conv
                chem(i,k,j,p_ole1)  = chem(i,k,j,p_ole1)                           &
                     +emis_ant(i,k,j,p_e_ole1)*conv
                chem(i,k,j,p_ole2)  = chem(i,k,j,p_ole2)                           &
                     +emis_ant(i,k,j,p_e_ole2)*conv
                chem(i,k,j,p_aro1)  = chem(i,k,j,p_aro1)                           &
                     +emis_ant(i,k,j,p_e_aro1)*conv
                chem(i,k,j,p_aro2)  = chem(i,k,j,p_aro2)                           &
                     +emis_ant(i,k,j,p_e_aro2)*conv
                chem(i,k,j,p_hcho)  = chem(i,k,j,p_hcho)                           &
                     +emis_ant(i,k,j,p_e_hcho)*conv
                chem(i,k,j,p_ccho)  = chem(i,k,j,p_ccho)                           &
                     +emis_ant(i,k,j,p_e_ccho)*conv
                chem(i,k,j,p_rcho)  = chem(i,k,j,p_rcho)                           &
                     +emis_ant(i,k,j,p_e_rcho)*conv
                chem(i,k,j,p_acet)  = chem(i,k,j,p_acet)                           &
                     +emis_ant(i,k,j,p_e_acet)*conv


                

                chem(i,k,j,p_mek)  = chem(i,k,j,p_mek)                             &
                     +emis_ant(i,k,j,p_e_mek)*conv

                chem(i,k,j,p_isoprene)  = chem(i,k,j,p_isoprene)                   &
                     +emis_ant(i,k,j,p_e_isoprene)*conv

                chem(i,k,j,p_terp)  = chem(i,k,j,p_terp)                           &
                     +emis_ant(i,k,j,p_e_terp)*conv

                chem(i,k,j,p_sesq)  = chem(i,k,j,p_sesq)                           &
                     +emis_ant(i,k,j,p_e_sesq)*conv

                chem(i,k,j,p_co)  = chem(i,k,j,p_co)                               &
                     +emis_ant(i,k,j,p_e_co)*conv

                chem(i,k,j,p_no)  = chem(i,k,j,p_no)                               &
                     +emis_ant(i,k,j,p_e_no)*conv*config_flags%uq_nox

                chem(i,k,j,p_etoh)  = chem(i,k,j,p_etoh)&
                     +emis_ant(i,k,j,p_e_etoh)*conv
                
                









                


                




              if((k.eq.kts).and.xland(i,j).eq.1.0)        & 
              chem(i,k,j,p_no)  = chem(i,k,j,p_no) +15.0/30.01*conv 
                                                             
              if((k.eq.kts).and.xland(i,j).eq.1.0)        & 
              chem(i,k,j,p_nh3)  = chem(i,k,j,p_nh3) +15.0/17.0*conv   

              if((k.eq.kts).and.xland(i,j).eq.1.0)        & 
              chem(i,k,j,p_dms)  = chem(i,k,j,p_dms) +0.24*conv  


                chem(i,k,j,p_no2)  = chem(i,k,j,p_no2)                             &
                     +emis_ant(i,k,j,p_e_no2)*conv*config_flags%uq_nox
                chem(i,k,j,p_phen)  = chem(i,k,j,p_phen)                           &
                     +emis_ant(i,k,j,p_e_phen)*conv
                chem(i,k,j,p_cres)  = chem(i,k,j,p_cres)                           &
                     +emis_ant(i,k,j,p_e_cres)*conv
                chem(i,k,j,p_meoh)  = chem(i,k,j,p_meoh)                           &
                     +emis_ant(i,k,j,p_e_meoh)*conv
                chem(i,k,j,p_gly)  = chem(i,k,j,p_gly)                             &
                     +emis_ant(i,k,j,p_e_gly)*conv
                chem(i,k,j,p_mgly)  = chem(i,k,j,p_mgly)                           &
                     +emis_ant(i,k,j,p_e_mgly)*conv
                chem(i,k,j,p_bacl)  = chem(i,k,j,p_bacl)                           &
                     +emis_ant(i,k,j,p_e_bacl)*conv
                chem(i,k,j,p_isoprod)  = chem(i,k,j,p_isoprod)                     &
                     +emis_ant(i,k,j,p_e_isoprod)*conv
                chem(i,k,j,p_methacro)  = chem(i,k,j,p_methacro)                   &
                     +emis_ant(i,k,j,p_e_methacro)*conv
                chem(i,k,j,p_mvk)  = chem(i,k,j,p_mvk)                             &
                     +emis_ant(i,k,j,p_e_mvk)*conv
                chem(i,k,j,p_prod2)  = chem(i,k,j,p_prod2)                         &
                     +emis_ant(i,k,j,p_e_prod2)*conv
                chem(i,k,j,p_ch4)  = chem(i,k,j,p_ch4)                             &
                     +emis_ant(i,k,j,p_e_ch4)*conv
                chem(i,k,j,p_bald)  = chem(i,k,j,p_bald)                           &
                     +emis_ant(i,k,j,p_e_bald)*conv
                chem(i,k,j,p_hcooh)  = chem(i,k,j,p_hcooh)                         &
                     +emis_ant(i,k,j,p_e_hcooh)*conv
                chem(i,k,j,p_cco_oh)  = chem(i,k,j,p_cco_oh)                       &
                     +emis_ant(i,k,j,p_e_cco_oh)*conv
                chem(i,k,j,p_rco_oh)  = chem(i,k,j,p_rco_oh)                       &
                     +emis_ant(i,k,j,p_e_rco_oh)*conv
                chem(i,k,j,p_nh3)  = chem(i,k,j,p_nh3)                       &
                     +emis_ant(i,k,j,p_e_nh3)*conv






         chem(i,k,j,p_pcg1_f_c)  =  chem(i,k,j,p_pcg1_f_c)                        &


        +(emis_ant(i,k,j,p_e_orgi_a)+emis_ant(i,k,j,p_e_orgj_a))*conv3*1.4          
                                                                                    
                                                                                    

         chem(i,k,j,p_pcg7_f_c)  =  chem(i,k,j,p_pcg7_f_c)                        &
                                   +(emis_ant(i,k,j,p_e_c2h6)*conv  &
                                   +emis_ant(i,k,j,p_e_c3h8)*conv  &
                                   +emis_ant(i,k,j,p_e_c2h2)*conv  &
                                   +emis_ant(i,k,j,p_e_alk3)*conv  &
                                   +emis_ant(i,k,j,p_e_alk4)*conv  &
                                   +emis_ant(i,k,j,p_e_alk5)*conv  &
                                   +emis_ant(i,k,j,p_e_ethene)*conv &
                                   +emis_ant(i,k,j,p_e_c3h6)*conv   &
                                   +emis_ant(i,k,j,p_e_ole1)*conv   &
                                   +emis_ant(i,k,j,p_e_ole2)*conv   &
                                   +emis_ant(i,k,j,p_e_aro1)*conv   &
                                   +emis_ant(i,k,j,p_e_aro2)*conv   &
                                   +emis_ant(i,k,j,p_e_hcho)*conv   &
                                   +emis_ant(i,k,j,p_e_ccho)*conv   &
                                   +emis_ant(i,k,j,p_e_rcho)*conv   &
                                   +emis_ant(i,k,j,p_e_acet)*conv   &
                                   +emis_ant(i,k,j,p_e_mek)*conv    &
                                   +emis_ant(i,k,j,p_e_isoprene)*conv &
                                   +emis_ant(i,k,j,p_e_terp)*conv     &
                                   +emis_ant(i,k,j,p_e_phen)*conv     &
                                   +emis_ant(i,k,j,p_e_cres)*conv     &
                                   +emis_ant(i,k,j,p_e_meoh)*conv     &
                                   +emis_ant(i,k,j,p_e_gly)*conv      &
                                   +emis_ant(i,k,j,p_e_mgly)*conv    &
                                   +emis_ant(i,k,j,p_e_bacl)*conv     &
                                   +emis_ant(i,k,j,p_e_isoprod)*conv  &
                                   +emis_ant(i,k,j,p_e_methacro)*conv &
                                   +emis_ant(i,k,j,p_e_mvk)*conv      &
                                   +emis_ant(i,k,j,p_e_prod2)*conv    &
                                   +emis_ant(i,k,j,p_e_bald)*conv     &
                                   +emis_ant(i,k,j,p_e_etoh)*conv     &
                                   +emis_ant(i,k,j,p_e_hcooh)*conv) *0.20 

             end do
          end do
       end do
       
456    CONTINUE

      else


       

       do j=jts,jte
          do i=its,ite
             do k=kts,min(config_flags%kemit,kte-ksub)
                conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                   &
                     +emis_ant(i,k,j,p_e_so2)*conv
                chem(i,k,j,p_co)  = chem(i,k,j,p_co)                     &
                     +emis_ant(i,k,j,p_e_co)*conv
                chem(i,k,j,p_no)  = chem(i,k,j,p_no)                     &
                     +emis_ant(i,k,j,p_e_no)*conv
                chem(i,k,j,p_hcho)  = chem(i,k,j,p_hcho)                 &
                     +emis_ant(i,k,j,p_e_hcho)*conv
             end do
          end do
       end do
      endif
      call add_biogenics(id,dtstep,dz8w,config_flags, rho_phy,chem,     &
            e_bio,ne_area,                                              &
            ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,               &
            ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,               &
            ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,       &
            ebio_sesq,ebio_mbo,                                         &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )



    CASE (GOCARTRACM_KPP,GOCARTRADM2)
       IF(config_flags%emiss_inpt_opt /= 3 ) then
       IF(config_flags%kemit .GT. kte-ksub) THEN
         k=config_flags%kemit
         write(message,'(" WARNING: EMISSIONS_DRIVER: KEMIT > KTE ",3i6)') kme,kte-ksub,k
         CALL WRF_MESSAGE (message)
       ENDIF
       call wrf_debug(15,'emissions_driver calling add_anthropogenics')
       call add_anthropogenics(id,dtstep,dz8w,config_flags,rho_phy,alt, &
            chem, emis_ant,emis_aircraft,                               &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )
       call wrf_debug(15,'emissions_driver calling add_biogenics')
       
       
       call add_biogenics(id,dtstep,dz8w,config_flags,rho_phy,chem,     &
            e_bio,ne_area,                                              &
            ebio_iso,ebio_oli,ebio_api,ebio_lim,ebio_xyl,               &
            ebio_hc3,ebio_ete,ebio_olt,ebio_ket,ebio_ald,               &
            ebio_hcho,ebio_eth,ebio_ora2,ebio_co,ebio_nr,ebio_no,       &
            ebio_sesq,ebio_mbo,                                         & 
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )

       end if 




       do j=jts,jte  
          do i=its,ite  
             do k=kts,min(config_flags%kemit,kte-ksub)
                conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                chem(i,k,j,p_bc1)  = chem(i,k,j,p_bc1)                     &
                     +(emis_ant(i,k,j,p_e_eci)+emis_ant(i,k,j,p_e_ecj))*alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_oc1)  = chem(i,k,j,p_oc1)                     &
                     +(emis_ant(i,k,j,p_e_orgj)+emis_ant(i,k,j,p_e_orgi))*alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_p25)  = chem(i,k,j,p_p25)                     &
                     +(emis_ant(i,k,j,p_e_pm25j)+emis_ant(i,k,j,p_e_pm25i) &
                     + emis_ant(i,k,j,p_e_no3j)+emis_ant(i,k,j,p_e_no3i))  &
                     *alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_sulf)  = chem(i,k,j,p_sulf)                   &
                     +(emis_ant(i,k,j,p_e_so4i)+emis_ant(i,k,j,p_e_so4j))*alt(i,k,j)*dtstep/dz8w(i,k,j)*mwdry/mw_so4_aer*1.e-3
             end do
          end do
       end do

    CASE (GOCART_SIMPLE)




       if(config_flags%emiss_opt <=  5  ) then
       do j=jts,jte  
          do i=its,ite  
             do k=kts,min(config_flags%kemit,kte-ksub)
                conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                   &
                     +emis_ant(i,k,j,p_e_so2)*conv
                chem(i,k,j,p_bc1)  = chem(i,k,j,p_bc1)                     &
                     +(emis_ant(i,k,j,p_e_eci)+emis_ant(i,k,j,p_e_ecj))*alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_oc1)  = chem(i,k,j,p_oc1)                     &
                     +(emis_ant(i,k,j,p_e_orgj)+emis_ant(i,k,j,p_e_orgi))*alt(i,k,j)*dtstep/dz8w(i,k,j)

                chem(i,k,j,p_p25)  = chem(i,k,j,p_p25)                     &
                     +(emis_ant(i,k,j,p_e_pm25j)+emis_ant(i,k,j,p_e_pm25i) &
                     +emis_ant(i,k,j,p_e_no3j)+emis_ant(i,k,j,p_e_no3i))   &
                     *alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_sulf)  = chem(i,k,j,p_sulf)                   &
                     +(emis_ant(i,k,j,p_e_so4i)+emis_ant(i,k,j,p_e_so4j))  &
                     *alt(i,k,j)*dtstep/dz8w(i,k,j)*mwdry/mw_so4_aer*1.e-3
             end do
          end do
       end do
       endif



       if(config_flags%emiss_opt == 6  ) then
       do j=jts,jte  
          do i=its,ite  
             do k=kts,min(config_flags%kemit,kte-ksub)
                conv = 4.828e-4/rho_phy(i,k,j)*dtstep/(dz8w(i,k,j)*60.)
                chem(i,k,j,p_so2)  = chem(i,k,j,p_so2)                   &
                     +emis_ant(i,k,j,p_e_so2)*conv
                chem(i,k,j,p_bc1)  = chem(i,k,j,p_bc1)                     &
                     +(emis_ant(i,k,j,p_e_bc))*alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_oc1)  = chem(i,k,j,p_oc1)                     &
                     +(emis_ant(i,k,j,p_e_oc))*alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_p25)  = chem(i,k,j,p_p25)                     &
                     +(emis_ant(i,k,j,p_e_pm_25))   &
                     *alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_p10)  = chem(i,k,j,p_p10)                     &
                     +(emis_ant(i,k,j,p_e_pm_10))   &
                     *alt(i,k,j)*dtstep/dz8w(i,k,j)
                chem(i,k,j,p_sulf)  = chem(i,k,j,p_sulf)                   &
                     +emis_ant(i,k,j,p_e_sulf)*conv
             end do
          end do
       end do
       endif


    CASE DEFAULT
       call wrf_debug(15,'emissions_driver NOT CALLING gas add_... routines')

    END SELECT gas_addemiss_select




    emiss_select: SELECT CASE(config_flags%emiss_inpt_opt)
    CASE (EMISS_INPT_CPTEC)
       call wrf_debug(15,'emissions_driver calling add_emiss_cptec')
       call add_emis_cptec(id,dtstep,ktau,dz8w,config_flags,curr_secs,   &
            rho_phy,chem,                                                &
            julday,gmt,xlat,xlong,t_phy,p_phy,                           &
            emis_ant,                                                    &
            ids,ide, jds,jde, kds,kde,                                   &
            ims,ime, jms,jme, kms,kme,                                   &
            its,ite, jts,jte, kts,kte                                    )
    CASE DEFAULT
       call wrf_debug(15,'emissions_driver not calling add_emiss_cptec')
    END SELECT emiss_select

    aer_addemiss_select: SELECT CASE(config_flags%chem_opt)

    CASE (RADM2SORG,RADM2SORG_AQ,RADM2SORG_AQCHEM,RADM2SORG_KPP, &
          RACMSORG_AQ,RACMSORG_AQCHEM_KPP,RACM_ESRLSORG_AQCHEM_KPP,RACMSORG_KPP,RACM_ESRLSORG_KPP,CBMZSORG,CBMZSORG_AQ, &
          CB05_SORG_AQ_KPP)
       call wrf_debug(15,'emissions_driver calling sorgam_addemiss')
       call sorgam_addemiss( id, dtstep, u10, v10, alt, dz8w, xland, chem,  &
            ebu,                                                        &
            slai,ust,smois,ivgtyp,isltyp,                               &
            emis_ant,dust_emiss_active,                                 &
            seasalt_emiss_active,config_flags%kemit,                    &
            config_flags%biomass_burn_opt,                              &
            config_flags%num_soil_layers,config_flags%emiss_opt,        &
            config_flags%dust_opt,                                      &
            ktau,p8w,u_phy,v_phy,rho_phy,g,dx,erod,                     &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )

    CASE (CB05_SORG_VBS_AQ_KPP)
       call wrf_debug(15,'emissions_driver calling sorgam_vbs_addemiss')
       call sorgam_vbs_addemiss( id, dtstep, u10, v10, alt, dz8w, xland, chem,  &
            ebu,                                                        &
            slai,ust,smois,ivgtyp,isltyp,                               &
            emis_ant,dust_emiss_active,                                 &
            seasalt_emiss_active,config_flags%kemit,                    &
            config_flags%biomass_burn_opt,                              &
            config_flags%num_soil_layers,config_flags%emiss_opt,        &
            config_flags%dust_opt,                                      &
            ktau,p8w,u_phy,v_phy,rho_phy,g,dx,erod,                     &
            emis_seas2,                                                 &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )


    CASE (RACM_SOA_VBS_KPP,RACM_SOA_VBS_AQCHEM_KPP,RACM_SOA_VBS_HET_KPP)
       call wrf_debug(15,'emissions_driver calling soa_vbs_addemiss')
       call soa_vbs_addemiss( id, dtstep, u10, v10, alt, dz8w, xland, chem,  &
            ebu,                                                        &
            slai,ust,smois,ivgtyp,isltyp,                               &
            emis_ant,dust_emiss_active,                                 &
            seasalt_emiss_active,config_flags%kemit,                    &
            config_flags%biomass_burn_opt,                              &
            config_flags%num_soil_layers,config_flags%emiss_opt,        &
            config_flags%dust_opt,                                      &
            ktau,p8w,u_phy,v_phy,rho_phy,g,dx,erod,                     &
            ids,ide, jds,jde, kds,kde,                                  &
            ims,ime, jms,jme, kms,kme,                                  &
            its,ite, jts,jte, kts,kte                                   )



    CASE (CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_KPP, CBMZ_MOSAIC_8BIN, & 
          CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
          CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, & 
          CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ,SAPRC99_MOSAIC_4BIN_VBS2_KPP,&
          MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP, &
          CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP,SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, &
          SAPRC99_MOSAIC_8BIN_VBS2_KPP, &
          SAPRC99_MOSAIC_20BIN_VBS2_AQ_KPP,SAPRC99_MOSAIC_20BIN_VBS2_KPP, & 
          SAPRC99_MOSAIC_12BIN_VBS2_AQ_KPP,SAPRC99_MOSAIC_12BIN_VBS2_KPP,SAPRC99_MOSAIC_4BIN_IEPOX_VBS_AQ_KPP, &
          SAPRC99_SIMPLESOM_MOSAIC_4BIN_AQ_KPP) 
       call wrf_debug(15,'emissions_driver calling mosaic_addemiss')
       call mosaic_addemiss( id, dtstep, u10, v10, alt, dz8w, xland,     &
            config_flags, chem, slai, ust, smois, ivgtyp, isltyp,        &
            emis_ant,ebu,config_flags%biomass_burn_opt,                  &
            config_flags%dust_opt,                                       &
            ktau,p8w,u_phy,v_phy,rho_phy,g,dx,erod,                      &
            dust_emiss_active, seasalt_emiss_active,                     &
            dust_flux, seas_flux,                                        &
            ids,ide, jds,jde, kds,kde,                                   &
            ims,ime, jms,jme, kms,kme,                                   &
            its,ite, jts,jte, kts,kte                                    )

    CASE ( CBMZ_CAM_MAM3_NOAQ, CBMZ_CAM_MAM3_AQ, CBMZ_CAM_MAM7_NOAQ, CBMZ_CAM_MAM7_AQ )
       call wrf_debug(15,'emissions_driver calling cam_mam_addemiss')
       call cam_mam_addemiss(                                            &
            id, dtstep, u10, v10, alt, dz8w, xland,                      &
            config_flags, chem, slai, ust, smois, ivgtyp, isltyp,        &
            emis_ant,ebio_iso,ebio_olt,ebio_oli,rho_phy,                 &
            dust_emiss_active, seasalt_emiss_active,                     &
            ids,ide, jds,jde, kds,kde,                                   &
            ims,ime, jms,jme, kms,kme,                                   &
            its,ite, jts,jte, kts,kte                                    )
       call wrf_debug(15,'emissions_driver backfrm cam_mam_addemiss')

    CASE DEFAULT
       call wrf_debug(15,'emissions_driver NOT CALLING aer add_... routines')

    END SELECT aer_addemiss_select


    CALL lightning_nox_driver ( &
          
            curr_secs=curr_secs, dt=dtstep, dx=dx, dy=dx,             &
            xlat=xlat, xlon=xlong, xland=xland, ht=ht,                &
            t_phy=t_phy, p_phy=p_phy, rho=rho_phy, u=u_phy, v=v_phy, w=vvel,        &
            z=z, moist=moist,                         &
            ic_flashrate=ic_flashrate, cg_flashrate=cg_flashrate,    &
          
            refl=refl_10cm,                                     &
          
            lightning_option=config_flags%lightning_opt,          &
            lightning_dt=config_flags%lightning_time_step,                  &
            lightning_start_seconds=config_flags%lightning_start_seconds,                    &
            N_IC=config_flags%N_IC, N_CG=config_flags%N_CG,          &
            lnox_opt=config_flags%lnox_opt, lnox_passive=config_flags%lnox_passive, &
          
            ltng_temp_upper=config_flags%ltng_temp_upper,                 &
            ltng_temp_lower=config_flags%ltng_temp_lower,                 &
            cellcount_method=config_flags%cellcount_method,               &
          
            ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,         &
            ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,         &
            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte,         &
          
            c_no=chem(:,:,:,p_no),                                 & 
            lnox_total=tracer(:,:,:,p_lnox_total),       &
            lnox_ic=tracer(:,:,:,p_lnox_ic),             &
            lnox_cg=tracer(:,:,:,p_lnox_cg)          &
          )

    END subroutine emissions_driver

END module module_emissions_driver
