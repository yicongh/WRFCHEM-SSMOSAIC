



      MODULE MODULE_MOSAIC_DRIVER
      
      USE MODULE_MOSAIC_SUMPM
      USE MODULE_STATE_DESCRIPTION, ONLY: num_moist, num_chem


      USE MODULE_CONFIGURE, ONLY: & 
      p_TERP_CN3_g, p_TERP_CN2_g, p_TERP_CN1_g, p_TERP_C00_g, p_TERP_C01_g, &
      p_TERP_C02_g, p_TERP_C03_g, p_TERP_C04_g, p_TERP_C05_g, p_TERP_C06_g, &
      p_TERP_C07_g, p_TERP_C08_g, p_TERP_C09_g, &
      p_TERP_CN3_a01, p_TERP_CN3_a02, p_TERP_CN3_a03, p_TERP_CN3_a04, &
      p_TERP_CN2_a01, p_TERP_CN2_a02, p_TERP_CN2_a03, p_TERP_CN2_a04, &
      p_TERP_CN1_a01, p_TERP_CN1_a02, p_TERP_CN1_a03, p_TERP_CN1_a04, &
      p_TERP_C00_a01, p_TERP_C00_a02, p_TERP_C00_a03, p_TERP_C00_a04, &
      p_TERP_C01_a01, p_TERP_C01_a02, p_TERP_C01_a03, p_TERP_C01_a04, &
      p_TERP_C02_a01, p_TERP_C02_a02, p_TERP_C02_a03, p_TERP_C02_a04, &
      p_TERP_C03_a01, p_TERP_C03_a02, p_TERP_C03_a03, p_TERP_C03_a04, &
      p_TERP_C04_a01, p_TERP_C04_a02, p_TERP_C04_a03, p_TERP_C04_a04, &
      p_TERP_C05_a01, p_TERP_C05_a02, p_TERP_C05_a03, p_TERP_C05_a04, &
      p_TERP_C06_a01, p_TERP_C06_a02, p_TERP_C06_a03, p_TERP_C06_a04, &
      p_TERP_C07_a01, p_TERP_C07_a02, p_TERP_C07_a03, p_TERP_C07_a04, &
      p_TERP_C08_a01, p_TERP_C08_a02, p_TERP_C08_a03, p_TERP_C08_a04, &
      p_TERP_C09_a01, p_TERP_C09_a02, p_TERP_C09_a03, p_TERP_C09_a04
      

      USE MODULE_DATA_MOSAIC_OTHER, ONLY: cAIRCLM

      CONTAINS
      
      SUBROUTINE drive_MOSAIC0(                         &
		id, curr_secs, ktau, dtstep, ktauc, dtstepc, config_flags, &
		t_phy, rho_phy, p_phy,                            &
		moist, chem, vbs_nbin,                            &
               swdown,                                            & 
               phw_bin1, phw_bin2,phw_bin3,phw_bin4,              &
               ph_no2,                                            &
               gammaaq,gammacld,                                  &
               clddropdia,                                        &
               betaaq,betacld,coatbin1,                           &
               coatbin2,coatbin3,coatbin4,                        &
                ion_prod_rate,                                    &
                JBNRATE, JBCHRATE, JTNRATE, JTCHRATE,             &
                JBNORGRATE, JBCHORGRATE, JTNRICRATE, JRATE,       &
                JBNRATEAC, JBCHRATEAC, JTNRATEAC, JTCHRATEAC,     &
                JBNORGRATEAC, JBCHORGRATEAC, JTNRICRATEAC, JRATEAC, & 
		ids,ide, jds,jde, kds,kde,                        &
		ims,ime, jms,jme, kms,kme,                        &
		its,ite, jts,jte, kts,kte                         )


      USE module_configure, only:  grid_config_rec_type,            &
			p_qv,                                         &
			p_so2, p_ho2, p_so4aj, p_corn, p_hcl, p_mtf,  &
			p_so4_a01, p_water_a01, p_num_a01,            &
			p_so4_a04, p_water_a04, p_num_a04

      USE module_state_description, only:  num_moist, num_chem
      USE module_data_mosaic_asect
      USE module_data_mosaic_other
      USE module_mosaic_therm, only:  aerchemistry, print_mosaic_stats, &
         iprint_mosaic_fe1, iprint_mosaic_perform_stats, &
         iprint_mosaic_diag1, iprint_mosaic_input_ok
      use module_mosaic2_driver, only:  mosaic2_aerchem_driver
      use module_mosaic_newnuc, only:  mosaic_newnuc_1clm
      use module_mosaic_coag, only:  mosaic_coag_1clm
      use module_peg_util, only:  peg_error_fatal, peg_message

      use module_data_mosaic_therm, only: glysoa_param,                 &
               glysoa_param_off, glysoa_param_simple, glysoa_param_complex
      use module_state_description, only: mozart_mosaic_4bin_kpp, &
         mozart_mosaic_4bin_aq_kpp

      implicit none




      INTEGER, INTENT(in) ::              &
		id, ktau, ktauc,                &
		ids, ide, jds, jde, kds, kde,   &
		ims, ime, jms, jme, kms, kme,   &
		its, ite, jts, jte, kts, kte,   &
                vbs_nbin(1)


    real(kind=8), intent(in) :: curr_secs
	real, intent(in) :: dtstep, dtstepc



	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme ) :: &
		t_phy, rho_phy, p_phy




    real, dimension (ims:ime, jms:jme ) ,                               &
        intent(in)             :: swdown


	real, intent(in),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
		moist



        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
           phw_bin1, phw_bin2, phw_bin3, phw_bin4,gammaaq,gammacld, &  
           clddropdia,                                              &
           betaaq,betacld,coatbin1,coatbin2,coatbin3,coatbin4

        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: ph_no2
 
	real, intent(inout),   &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem



	type(grid_config_rec_type), intent(in) :: config_flags


        real, intent(in),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                ion_prod_rate  


        real, intent(inout),   &
                dimension( ims:ime, kms:kme, jms:jme ) :: &
                JBNRATE, JBCHRATE, JTNRATE, JTCHRATE, &
                JBNORGRATE, JBCHORGRATE, JTNRICRATE, JRATE, &
                JBNRATEAC, JBCHRATEAC, JTNRATEAC, JTCHRATEAC, &
                JBNORGRATEAC, JBCHORGRATEAC, JTNRICRATEAC, JRATEAC 





	integer :: i, idum, istat, it, j, jt, k, l, n
	integer :: k_pegshift, kclm_calcbgn, kclm_calcend
	integer :: ktmaps, ktmape
	integer :: levdbg_err, levdbg_info
	integer :: i_force_dump, mode_force_dump
	integer :: idiagaa_dum, idiagbb_dum, ijcount_dum
	integer, parameter :: debug_level=0


	integer :: aercoag_onoff = 1
	integer :: aernewnuc_onoff = 1
	
	real :: dtchem, dtcoag, dtnuc
	real :: dum
	real :: rsub0(l2maxd,kmaxd,nsubareamaxd)

	character*100 msg

      INTEGER :: iNOW,jNOW


        if ( config_flags%aerchem_onoff <= 0 ) return   




         phvar1(:,:)=huge(1.0)
         phvar2(:,:)=huge(1.0)
         phvar3(:,:)=huge(1.0)
         phvar4(:,:)=huge(1.0)
         gammaiepoxaq(:,:)=huge(1.0)
         gammaiepoxcld(:,:)=huge(1.0)
         effdropdia(:,:)=huge(1.0)
         betaiepoxaq(:,:)=huge(1.0)
         betaiepoxcld(:,:)=huge(1.0)
         iepoxcoatbin1(:,:)=huge(1.0)
        iepoxcoatbin2(:,:)=huge(1.0)
        iepoxcoatbin3(:,:)=huge(1.0)
        iepoxcoatbin4(:,:)=huge(1.0)


 




        i = max( config_flags%mosaic_aerchem_optaa, 0 )
        if ( mod(i,10) == 2 ) then   
            call mosaic2_aerchem_driver(                          &
                id, curr_secs, ktau, dtstep, ktauc, dtstepc,      &
                config_flags,                                     &
                t_phy, rho_phy, p_phy,                            &
                moist, chem,                                      &
                ids,ide, jds,jde, kds,kde,                        &
                ims,ime, jms,jme, kms,kme,                        &
                its,ite, jts,jte, kts,kte                         )
            return
        end if




    if (debug_level .lt. 15) then 
       iprint_mosaic_fe1 = 1
       iprint_mosaic_perform_stats = 0
       iprint_mosaic_diag1 = 0
       iprint_mosaic_input_ok = 0
    end if

   glysoa_param = glysoa_param_off

   if (config_flags%chem_opt == mozart_mosaic_4bin_kpp) &
     glysoa_param = glysoa_param_simple
   if (config_flags%chem_opt == mozart_mosaic_4bin_aq_kpp) &
     glysoa_param = glysoa_param_complex


        aercoag_onoff = 1

        if ( config_flags%newnuc_method > 0 ) then  
           aernewnuc_onoff = 1
        else
           aernewnuc_onoff = 0
        end if

        i = max( config_flags%mosaic_aerchem_optaa, 0 )
        if (i >= 10000) then
            j = mod(i,100)/10  
            if (j == 0) msectional =  0  
            if (j == 1) msectional = 10  
            if (j == 2) msectional = 20  

            j = mod(i,1000)/100  
            if (j == 0) aernewnuc_onoff = 0

            j = mod(i,10000)/1000  
            if (j == 0) aercoag_onoff = 0
        end if 



	ktmaps = kts
	ktmape = kte




	k_pegshift = k_pegbegin - kts 
	kclm_calcbgn = kts + k_pegshift
	kclm_calcend = kte + k_pegshift


	mode_force_dump = 0
	levdbg_err = 0
        levdbg_info = 15



	iymdcur = 1 + int( curr_secs/86400._8 )
	ihmscur = nint( mod( curr_secs, 86400._8 ) )

	t = curr_secs
	ncorecnt = ktau - 1



	itot = ite
	jtot = jte
	nsubareas = 1


      dtchem = dtstepc

      i_force_dump = 0



      LOOP1: DO jt = jts, jte
      LOOP2: DO it = its, ite



      CALL mapaer_tofrom_host( 0,                       &
		ims,ime, jms,jme, kms,kme,                    &
		its,ite, jts,jte, kts,kte,                    &
		it,      jt,      ktmaps,ktmape,              &
		num_moist, num_chem, moist, chem,             &
		t_phy, p_phy, rho_phy,swdown,phw_bin1,        &
                phw_bin2, phw_bin3, phw_bin4,ph_no2,gammaaq,gammacld, &
                clddropdia,                                    &
                betaaq,betacld,coatbin1,coatbin2,coatbin3,coatbin4)                         


      rsub0(:,:,:) = rsub(:,:,:)
      
      idiagbb_dum = 110



      call aerchemistry( it, jt, kclm_calcbgn, kclm_calcend,   &
                         dtchem, vbs_nbin,config_flags  )
      


      if (aernewnuc_onoff > 0) then
         dtnuc = dtchem
	 call mosaic_newnuc_1clm( istat,  &
		it, jt, kclm_calcbgn, kclm_calcend,   &
		idiagbb_dum, dtchem, dtnuc, rsub0,   &
		id, ktau, ktauc, its, ite, jts, jte, kts, kte, &
                config_flags%newnuc_method, ion_prod_rate(it,kclm_calcbgn:kclm_calcend,jt), &
                JBNRATE(it,kclm_calcbgn:kclm_calcend,jt), JBCHRATE(it,kclm_calcbgn:kclm_calcend,jt), &
                JTNRATE(it,kclm_calcbgn:kclm_calcend,jt), JTCHRATE(it,kclm_calcbgn:kclm_calcend,jt), &
                JBNORGRATE(it,kclm_calcbgn:kclm_calcend,jt), JBCHORGRATE(it,kclm_calcbgn:kclm_calcend,jt), &
                JTNRICRATE(it,kclm_calcbgn:kclm_calcend,jt), JRATE(it,kclm_calcbgn:kclm_calcend,jt), &
                JBNRATEAC(it,kclm_calcbgn:kclm_calcend,jt), JBCHRATEAC(it,kclm_calcbgn:kclm_calcend,jt), &
                JTNRATEAC(it,kclm_calcbgn:kclm_calcend,jt), JTCHRATEAC(it,kclm_calcbgn:kclm_calcend,jt), &
                JBNORGRATEAC(it,kclm_calcbgn:kclm_calcend,jt), JBCHORGRATEAC(it,kclm_calcbgn:kclm_calcend,jt), &
                JTNRICRATEAC(it,kclm_calcbgn:kclm_calcend,jt), JRATEAC(it,kclm_calcbgn:kclm_calcend,jt) ) 
          
      end if
      


      if (aercoag_onoff > 0) then
	    dtcoag = dtchem
	    call mosaic_coag_1clm( istat,  &
		it, jt, kclm_calcbgn, kclm_calcend,   &
		idiagbb_dum, dtchem, dtcoag,   &
		id, ktau, ktauc, its, ite, jts, jte, kts, kte )
      end if



      call mapaer_tofrom_host( 1,                       &
		ims,ime, jms,jme, kms,kme,                    &
		its,ite, jts,jte, kts,kte,                    &
		it,      jt,      ktmaps,ktmape,              &
		num_moist, num_chem, moist, chem,             &
		t_phy, p_phy, rho_phy,swdown,phw_bin1,        &
                phw_bin2,phw_bin3,phw_bin4,ph_no2,gammaaq,gammacld,  &
                clddropdia,                                   &
                betaaq,betacld,coatbin1,coatbin2,coatbin3,coatbin4 )


      END DO LOOP2
      END DO LOOP1

      RETURN
      END SUBROUTINE drive_MOSAIC0








	subroutine mapaer_tofrom_host( imap,                  &
		ims,ime, jms,jme, kms,kme,                    &
		its,ite, jts,jte, kts,kte,                    &
		it,      jt,      ktmaps,ktmape,              &
		num_moist, num_chem, moist, chem,             &
		t_phy, p_phy, rho_phy,swdown,phw_bin1,        &
                phw_bin2,phw_bin3,phw_bin4,ph_no2,gammaaq,gammacld,  &
                clddropdia,                                   &
                betaaq,betacld,coatbin1,coatbin2,coatbin3,coatbin4 )

        use module_configure, only:   &
		p_qv, p_qc, p_h2so4,p_sulf, p_hno3, p_hcl, p_nh3, p_o3,   &
		p_so2, p_h2o2, p_hcho, p_ho, p_ho2, p_no3,   &
		p_no, p_no2, p_hono, p_pan,  &
       p_pcg1_b_c,p_pcg2_b_c,p_pcg3_b_c,p_pcg4_b_c,p_pcg5_b_c,p_pcg6_b_c, &
       p_pcg7_b_c,p_pcg8_b_c,p_pcg9_b_c,p_opcg1_b_c,p_opcg2_b_c,p_opcg3_b_c, &
       p_opcg4_b_c,p_opcg5_b_c,p_opcg6_b_c,p_opcg7_b_c,p_opcg8_b_c,p_pcg1_b_o,&
       p_pcg2_b_o,p_pcg3_b_o,p_pcg4_b_o,p_pcg5_b_o,p_pcg6_b_o,p_pcg7_b_o,&
       p_pcg8_b_o,p_pcg9_b_o,p_opcg1_b_o,p_opcg2_b_o,p_opcg3_b_o,p_opcg4_b_o, &
       p_opcg5_b_o,p_opcg6_b_o,p_opcg7_b_o,p_opcg8_b_o,p_pcg1_f_c,p_pcg2_f_c, &
       p_pcg3_f_c,p_pcg4_f_c,p_pcg5_f_c,p_pcg6_f_c,p_pcg7_f_c,p_pcg8_f_c, &
       p_pcg9_f_c, p_opcg1_f_c,p_opcg2_f_c,p_opcg3_f_c,p_opcg4_f_c, &
       p_opcg5_f_c,p_opcg6_f_c,p_opcg7_f_c,p_opcg8_f_c,p_pcg1_f_o, &
       p_pcg2_f_o,p_pcg3_f_o,p_pcg4_f_o,p_pcg5_f_o,p_pcg6_f_o,p_pcg7_f_o, &
       p_pcg8_f_o,p_pcg9_f_o,p_opcg1_f_o,p_opcg2_f_o,p_opcg3_f_o,p_opcg4_f_o,&
       p_opcg5_f_o,p_opcg6_f_o,p_opcg7_f_o,p_opcg8_f_o, &
       p_smpa,p_smpbb, &
       p_iepox,p_tetrol,p_cldtetrol, &
       p_gly, &
       p_ant1_c,p_ant2_c,p_ant3_c,p_ant4_c,p_ant1_o,p_ant2_o,p_ant3_o,p_ant4_o,&
       p_biog1_c,p_biog2_c,p_biog3_c,p_biog4_c,p_biog1_o, &
       p_biog2_o,p_biog3_o,p_biog4_o, &

       p_n2o5, p_clno2, &
       p_asoaX_a01, p_asoaX_a02, p_asoaX_a03, p_asoaX_a04, &
       p_asoa1_a01, p_asoa1_a02, p_asoa1_a03, p_asoa1_a04, &
       p_asoa2_a01, p_asoa2_a02, p_asoa2_a03, p_asoa2_a04, &
       p_asoa3_a01, p_asoa3_a02, p_asoa3_a03, p_asoa3_a04, &
       p_asoa4_a01, p_asoa4_a02, p_asoa4_a03, p_asoa4_a04, &
       p_bsoaX_a01, p_bsoaX_a02, p_bsoaX_a03, p_bsoaX_a04, &
       p_bsoa1_a01, p_bsoa1_a02, p_bsoa1_a03, p_bsoa1_a04, &
       p_bsoa2_a01, p_bsoa2_a02, p_bsoa2_a03, p_bsoa2_a04, &
       p_bsoa3_a01, p_bsoa3_a02, p_bsoa3_a03, p_bsoa3_a04, &
       p_bsoa4_a01, p_bsoa4_a02, p_bsoa4_a03, p_bsoa4_a04, &
       p_cvasoaX, p_cvasoa1, p_cvasoa2, p_cvasoa3, p_cvasoa4, &
       p_cvbsoaX, p_cvbsoa1, p_cvbsoa2, p_cvbsoa3, p_cvbsoa4, &
       p_TERP_CN3_g 

	use module_state_description, only:  param_first_scalar
	use module_data_mosaic_asect
	use module_data_mosaic_other
	use module_mosaic_csuesat, only:  esat_gchm
	use module_peg_util, only:  peg_error_fatal, peg_message

	implicit none




	integer, intent(in) :: imap

	integer, intent(in) :: num_moist, num_chem
	integer, intent(in) :: ims, ime, jms, jme, kms, kme
	integer, intent(in) :: its, ite, jts, jte, kts, kte

	integer, intent(in) :: it, jt, ktmaps, ktmape

	real, intent(in), dimension( ims:ime, kms:kme, jms:jme ) :: &
		t_phy, rho_phy, p_phy

    real, dimension (ims:ime, jms:jme ) ,                               &
        intent(in),optional             :: swdown

	real, intent(in), &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_moist ) :: &
		moist
 
	real, intent(inout), &
		dimension( ims:ime, kms:kme, jms:jme, 1:num_chem ) :: &
		chem

   real, dimension (ims:ime, kms:kme,jms:jme ) ,                               &
        intent(inout),optional             :: phw_bin1, phw_bin2,           &
         phw_bin3, phw_bin4,gammaaq,gammacld,clddropdia,betaaq,betacld,coatbin1, &
         coatbin2,coatbin3,coatbin4

   real, dimension (ims:ime, kms:kme,jms:jme ) ,                               &
        intent(in),optional             ::  ph_no2


	integer ido_l, idum, iphase, itype,   &
		k, k1, k2, kt, kt1, kt2, k_pegshift, l, n
	integer p1st
	real dum, dumesat, dumrsat, dumrelhum, onemeps

	real factdens, factpres, factmoist, factgas,   &
		factaerso4, factaerno3, factaercl, factaermsa,   &
		factaerco3, factaernh4, factaerna, factaerca,   &
		factaeroin, factaeroc, factaerbc,   &
       factaerpcg1_b_c,factaerpcg2_b_c,factaerpcg3_b_c,factaerpcg4_b_c,factaerpcg5_b_c,factaerpcg6_b_c, &
       factaerpcg7_b_c,factaerpcg8_b_c,factaerpcg9_b_c,factaeropcg1_b_c,factaeropcg2_b_c,factaeropcg3_b_c, &
       factaeropcg4_b_c,factaeropcg5_b_c,factaeropcg6_b_c,factaeropcg7_b_c,factaeropcg8_b_c, &
       factaerpcg1_b_o,factaerpcg2_b_o,factaerpcg3_b_o,factaerpcg4_b_o,factaerpcg5_b_o,factaerpcg6_b_o, &
       factaerpcg7_b_o,factaerpcg8_b_o,factaerpcg9_b_o,factaeropcg1_b_o,factaeropcg2_b_o,factaeropcg3_b_o,&
       factaeropcg4_b_o,factaeropcg5_b_o,factaeropcg6_b_o,factaeropcg7_b_o,factaeropcg8_b_o, &
       factaerpcg1_f_c,factaerpcg2_f_c,factaerpcg3_f_c,factaerpcg4_f_c,factaerpcg5_f_c,factaerpcg6_f_c,&
       factaerpcg7_f_c,factaerpcg8_f_c,factaerpcg9_f_c,factaeropcg1_f_c,factaeropcg2_f_c,factaeropcg3_f_c,&
       factaeropcg4_f_c,factaeropcg5_f_c,factaeropcg6_f_c,factaeropcg7_f_c,factaeropcg8_f_c,&
       factaerpcg1_f_o,factaerpcg2_f_o,factaerpcg3_f_o,factaerpcg4_f_o,factaerpcg5_f_o,factaerpcg6_f_o,&
       factaerpcg7_f_o,factaerpcg8_f_o,factaerpcg9_f_o,factaeropcg1_f_o,factaeropcg2_f_o,factaeropcg3_f_o,&
       factaeropcg4_f_o,factaeropcg5_f_o,factaeropcg6_f_o,factaeropcg7_f_o,factaeropcg8_f_o, &
       factaersmpa,factaersmpbb, &
       factaeriepox,factaeriepoxos,factaertetrol,factaergly, &
       factaercldiepoxos,factaercldtetrol, &
       factaerglyr1, factaerglyr2, factaerglysfc, factaerglynh4, factaerglyoh, &
       factaerant1_c,factaerant2_c,factaerant3_c,factaerant4_c, &
       factaerant1_o,factaerant2_o,factaerant3_o,factaerant4_o, &
       factaerbiog1_c,factaerbiog2_c,factaerbiog3_c,factaerbiog4_c, &
       factaerbiog1_o,factaerbiog2_o,factaerbiog3_o,factaerbiog4_o, &
       factaertanv, factaersopnv,factaerternv,factaerseqnv,factaersianv, &
       factaerasoaX,factaerasoa1,factaerasoa2,factaerasoa3,factaerasoa4, &
       factaerbsoaX,factaerbsoa1,factaerbsoa2,factaerbsoa3,factaerbsoa4, &
		factaerhysw, factaerwater, factaernum, &
       factaerterp_Cn3 

	real, parameter :: eps=0.622

	character*80 msg






	factdens = 28.966e3      
	factpres = 0.1           
	factmoist = eps          
	factgas = 1.0e6          





	factaernum = 1000./28.966 
	

	dum = factaernum*1.0e6   
	factaerso4   = dum*mw_so4_aer
	factaerno3   = dum*mw_no3_aer
	factaercl    = dum*mw_cl_aer
	factaermsa   = dum*mw_msa_aer
	factaerco3   = dum*mw_co3_aer
	factaernh4   = dum*mw_nh4_aer
	factaerna    = dum*mw_na_aer
	factaerca    = dum*mw_ca_aer
	factaeroin   = dum*mw_oin_aer
	factaeroc    = dum*mw_oc_aer
	factaerbc    = dum*mw_bc_aer
	factaerhysw  = dum*mw_water_aer
	factaerwater = dum*mw_water_aer
        factaertanv=   dum*mw_tanv_aer
        factaersopnv=   dum*mw_sopnv_aer
        factaerternv=   dum*mw_ternv_aer
        factaerseqnv=   dum*mw_seqnv_aer
        factaersianv=   dum*mw_sianv_aer
        factaerpcg1_b_c=dum*mw_pcg1_b_c_aer
        factaerpcg2_b_c=dum*mw_pcg2_b_c_aer
        factaerpcg3_b_c=dum*mw_pcg3_b_c_aer
        factaerpcg4_b_c=dum*mw_pcg4_b_c_aer
        factaerpcg5_b_c=dum*mw_pcg5_b_c_aer
        factaerpcg6_b_c=dum*mw_pcg6_b_c_aer
        factaerpcg7_b_c=dum*mw_pcg7_b_c_aer
        factaerpcg8_b_c=dum*mw_pcg8_b_c_aer
        factaerpcg9_b_c=dum*mw_pcg9_b_c_aer
        factaerpcg1_b_o=dum*mw_pcg1_b_o_aer
        factaerpcg2_b_o=dum*mw_pcg2_b_o_aer
        factaerpcg3_b_o=dum*mw_pcg3_b_o_aer
        factaerpcg4_b_o=dum*mw_pcg4_b_o_aer
        factaerpcg5_b_o=dum*mw_pcg5_b_o_aer
        factaerpcg6_b_o=dum*mw_pcg6_b_o_aer
        factaerpcg7_b_o=dum*mw_pcg7_b_o_aer
        factaerpcg8_b_o=dum*mw_pcg8_b_o_aer
        factaerpcg9_b_o=dum*mw_pcg9_b_o_aer
        factaeropcg1_b_c=dum*mw_opcg1_b_c_aer
        factaeropcg2_b_c=dum*mw_opcg2_b_c_aer
        factaeropcg3_b_c=dum*mw_opcg3_b_c_aer
        factaeropcg4_b_c=dum*mw_opcg4_b_c_aer
        factaeropcg5_b_c=dum*mw_opcg5_b_c_aer
        factaeropcg6_b_c=dum*mw_opcg6_b_c_aer
        factaeropcg7_b_c=dum*mw_opcg7_b_c_aer
        factaeropcg8_b_c=dum*mw_opcg8_b_c_aer
        factaeropcg1_b_o=dum*mw_opcg1_b_o_aer
        factaeropcg2_b_o=dum*mw_opcg2_b_o_aer
        factaeropcg3_b_o=dum*mw_opcg3_b_o_aer
        factaeropcg4_b_o=dum*mw_opcg4_b_o_aer
        factaeropcg5_b_o=dum*mw_opcg5_b_o_aer
        factaeropcg6_b_o=dum*mw_opcg6_b_o_aer
        factaeropcg7_b_o=dum*mw_opcg7_b_o_aer
        factaeropcg8_b_o=dum*mw_opcg8_b_o_aer
        factaerpcg1_f_c=dum*mw_pcg1_f_c_aer
        factaerpcg2_f_c=dum*mw_pcg2_f_c_aer
        factaerpcg3_f_c=dum*mw_pcg3_f_c_aer
        factaerpcg4_f_c=dum*mw_pcg4_f_c_aer
        factaerpcg5_f_c=dum*mw_pcg5_f_c_aer
        factaerpcg6_f_c=dum*mw_pcg6_f_c_aer
        factaerpcg7_f_c=dum*mw_pcg7_f_c_aer
        factaerpcg8_f_c=dum*mw_pcg8_f_c_aer
        factaerpcg9_f_c=dum*mw_pcg9_f_c_aer
        factaerpcg1_f_o=dum*mw_pcg1_f_o_aer
        factaerpcg2_f_o=dum*mw_pcg2_f_o_aer
        factaerpcg3_f_o=dum*mw_pcg3_f_o_aer
        factaerpcg4_f_o=dum*mw_pcg4_f_o_aer
        factaerpcg5_f_o=dum*mw_pcg5_f_o_aer
        factaerpcg6_f_o=dum*mw_pcg6_f_o_aer
        factaerpcg7_f_o=dum*mw_pcg7_f_o_aer
        factaerpcg8_f_o=dum*mw_pcg8_f_o_aer
        factaerpcg9_f_o=dum*mw_pcg9_f_o_aer
        factaeropcg1_f_c=dum*mw_opcg1_f_c_aer
        factaeropcg2_f_c=dum*mw_opcg2_f_c_aer
        factaeropcg3_f_c=dum*mw_opcg3_f_c_aer
        factaeropcg4_f_c=dum*mw_opcg4_f_c_aer
        factaeropcg5_f_c=dum*mw_opcg5_f_c_aer
        factaeropcg6_f_c=dum*mw_opcg6_f_c_aer
        factaeropcg7_f_c=dum*mw_opcg7_f_c_aer
        factaeropcg8_f_c=dum*mw_opcg8_f_c_aer
        factaeropcg1_f_o=dum*mw_opcg1_f_o_aer
        factaeropcg2_f_o=dum*mw_opcg2_f_o_aer
        factaeropcg3_f_o=dum*mw_opcg3_f_o_aer
        factaeropcg4_f_o=dum*mw_opcg4_f_o_aer
        factaeropcg5_f_o=dum*mw_opcg5_f_o_aer
        factaeropcg6_f_o=dum*mw_opcg6_f_o_aer
        factaeropcg7_f_o=dum*mw_opcg7_f_o_aer
        factaeropcg8_f_o=dum*mw_opcg8_f_o_aer
        factaersmpa=dum*mw_smpa_aer
        factaersmpbb=dum*mw_smpbb_aer
        factaeriepox=dum*mw_iepox_aer
        factaeriepoxos=dum*mw_iepoxos_aer
        factaertetrol=dum*mw_tetrol_aer
        factaercldiepoxos=dum*mw_iepoxos_aer
        factaercldtetrol=dum*mw_tetrol_aer
        factaergly=dum*mw_gly_aer
        factaerglyr1=dum*mw_glysoa_r1_aer
        factaerglyr2=dum*mw_glysoa_r2_aer
        factaerglysfc=dum*mw_glysoa_sfc_aer
        factaerglynh4=dum*mw_glysoa_nh4_aer
        factaerglyoh=dum*mw_glysoa_oh_aer
        factaerant1_c=dum*mw_ant1_c_aer
        factaerant2_c=dum*mw_ant2_c_aer
        factaerant3_c=dum*mw_ant3_c_aer
        factaerant4_c=dum*mw_ant4_c_aer
        factaerant1_o=dum*mw_ant1_o_aer
        factaerant2_o=dum*mw_ant2_o_aer
        factaerant3_o=dum*mw_ant3_o_aer
        factaerant4_o=dum*mw_ant4_o_aer
        factaerbiog1_c=dum*mw_biog1_c_aer
        factaerbiog2_c=dum*mw_biog2_c_aer
        factaerbiog3_c=dum*mw_biog3_c_aer
        factaerbiog4_c=dum*mw_biog4_c_aer
        factaerbiog1_o=dum*mw_biog1_o_aer
        factaerbiog2_o=dum*mw_biog2_o_aer
        factaerbiog3_o=dum*mw_biog3_o_aer
        factaerbiog4_o=dum*mw_biog4_o_aer
        factaerasoaX=dum*mw_asoaX_aer
        factaerasoa1=dum*mw_asoa1_aer
        factaerasoa2=dum*mw_asoa2_aer
        factaerasoa3=dum*mw_asoa3_aer
        factaerasoa4=dum*mw_asoa4_aer
        factaerbsoaX=dum*mw_bsoaX_aer
        factaerbsoa1=dum*mw_bsoa1_aer
        factaerbsoa2=dum*mw_bsoa2_aer
        factaerbsoa3=dum*mw_bsoa3_aer
        factaerbsoa4=dum*mw_bsoa4_aer

	
	factaerterp_Cn3=dum*mw_terp_Cn3_aer





	if (aboxtest_units_convert .eq. 10) then
	    factdens = 1.0
	    factpres = 1.0
	    factmoist = 1.0
	    factgas = 1.0
	    factaernum = 1.0
	    factaerso4   = 1.0
	    factaerno3   = 1.0
	    factaercl    = 1.0
	    factaermsa   = 1.0
	    factaerco3   = 1.0
	    factaernh4   = 1.0
	    factaerna    = 1.0
	    factaerca    = 1.0
	    factaeroin   = 1.0
	    factaeroc    = 1.0
	    factaerbc    = 1.0
	    factaerhysw  = 1.0
	    factaerwater = 1.0
            factaerpcg1_b_c=1.0
            factaerpcg2_b_c=1.0
            factaerpcg3_b_c=1.0
            factaerpcg4_b_c=1.0
            factaerpcg5_b_c=1.0
            factaerpcg6_b_c=1.0
            factaerpcg7_b_c=1.0
            factaerpcg8_b_c=1.0
            factaerpcg9_b_c=1.0
            factaerpcg1_b_o=1.0
            factaerpcg2_b_o=1.0
            factaerpcg3_b_o=1.0
            factaerpcg4_b_o=1.0
            factaerpcg5_b_o=1.0
            factaerpcg6_b_o=1.0
            factaerpcg7_b_o=1.0
            factaerpcg8_b_o=1.0
            factaerpcg9_b_o=1.0
            factaeropcg1_b_c=1.0
            factaeropcg2_b_c=1.0
            factaeropcg3_b_c=1.0
            factaeropcg4_b_c=1.0
            factaeropcg5_b_c=1.0
            factaeropcg6_b_c=1.0
            factaeropcg7_b_c=1.0
            factaeropcg8_b_c=1.0
            factaeropcg1_b_o=1.0
            factaeropcg2_b_o=1.0
            factaeropcg3_b_o=1.0
            factaeropcg4_b_o=1.0
            factaeropcg5_b_o=1.0
            factaeropcg6_b_o=1.0
            factaeropcg7_b_o=1.0
            factaeropcg8_b_o=1.0
            factaerpcg1_f_c=1.0
            factaerpcg2_f_c=1.0
            factaerpcg3_f_c=1.0
            factaerpcg4_f_c=1.0
            factaerpcg5_f_c=1.0
            factaerpcg6_f_c=1.0
            factaerpcg7_f_c=1.0
            factaerpcg8_f_c=1.0
            factaerpcg9_f_c=1.0
            factaerpcg1_f_o=1.0
            factaerpcg2_f_o=1.0
            factaerpcg3_f_o=1.0
            factaerpcg4_f_o=1.0
            factaerpcg5_f_o=1.0
            factaerpcg6_f_o=1.0
            factaerpcg7_f_o=1.0
            factaerpcg8_f_o=1.0
            factaerpcg9_f_o=1.0
            factaeropcg1_f_c=1.0
            factaeropcg2_f_c=1.0
            factaeropcg3_f_c=1.0
            factaeropcg4_f_c=1.0
            factaeropcg5_f_c=1.0
            factaeropcg6_f_c=1.0
            factaeropcg7_f_c=1.0
            factaeropcg8_f_c=1.0
            factaeropcg1_f_o=1.0
            factaeropcg2_f_o=1.0
            factaeropcg3_f_o=1.0
            factaeropcg4_f_o=1.0
            factaeropcg5_f_o=1.0
            factaeropcg6_f_o=1.0
            factaeropcg7_f_o=1.0
            factaeropcg8_f_o=1.0
            factaersmpa=1.0
            factaersmpbb=1.0
            factaeriepox=1.0
            factaeriepoxos=1.0
            factaertetrol=1.0
            factaercldiepoxos=1.0
            factaercldtetrol=1.0
            factaergly=1.0
            factaertanv=   1.0
            factaersopnv=   1.0
            factaerternv=   1.0
            factaerseqnv=   1.0
            factaersianv=   1.0
            factaerglyr1=1.0
            factaerglyr2=1.0
            factaerglysfc=1.0
            factaerglynh4=1.0
            factaerglyoh=1.0
            factaerant1_c=1.0
            factaerant2_c=1.0
            factaerant3_c=1.0
            factaerant4_c=1.0
            factaerant1_o=1.0
            factaerant2_o=1.0
            factaerant3_o=1.0
            factaerant4_o=1.0
            factaerbiog1_c=1.0
            factaerbiog2_c=1.0
            factaerbiog3_c=1.0
            factaerbiog4_c=1.0
            factaerbiog1_o=1.0
            factaerbiog2_o=1.0
            factaerbiog3_o=1.0
            factaerbiog4_o=1.0
            factaerasoaX=1.0
            factaerasoa1=1.0
            factaerasoa2=1.0
            factaerasoa3=1.0
            factaerasoa4=1.0
            factaerbsoaX=1.0
            factaerbsoa1=1.0
            factaerbsoa2=1.0
            factaerbsoa3=1.0
            factaerbsoa4=1.0
	    
	    factaerterp_Cn3=1.0



	end if





	k_pegshift = k_pegbegin - kts



	ktot = (kte-1) + k_pegshift

	if ((kte > kmaxd) .or. (ktot > kmaxd) .or. (ktot <= 0)) then
	    write( msg, '(a,4i5)' )   &
		'*** subr mapaer_tofrom_host -- ' //   &
		'ktot, kmaxd, kts, kte', ktot, kmaxd, kts, kte
	    call peg_message( lunerr, msg )
	    msg = '*** subr mosaic_aerchem_driver -- ' //   &
		'kte>kmaxd OR ktot>kmaxd OR ktot<=0'
	    call peg_error_fatal( lunerr, msg )
	end if



	kt1 = ktmaps
	kt2 = ktmape
	k1 = kt1 + k_pegshift
	k2 = kt2 + k_pegshift

	if (imap .gt. 0) goto 2000
 





	rsub(:,:,:) = 0.0
	cairclm(:) = 0.0
	ptotclm(:) = 0.0
	afracsubarea(:,:) = 0.0
	relhumclm(:) = aboxtest_min_relhum
	rcldwtr_sub(:,:) = 0.0
        swdown_cell=0.0
        phno2_box(:) = 0.0

	adrydens_sub( :,:,:,:) = 0.0
	aqmassdry_sub(:,:,:,:) = 0.0
	aqvoldry_sub( :,:,:,:) = 0.0







	if ((aboxtest_map_method .eq. 2) .or.   &
	    (aboxtest_map_method .eq. 3)) then
	    do l = 2, num_chem
		rsub(l,k1:k2,1) = chem(it,kt1:kt2,jt,l)/factgas
	    end do
	end if

	p1st = param_first_scalar
	if (aboxtest_map_method .ne. 2) then
            if (p_h2so4 .ge. p1st) then  
                rsub(kh2so4,k1:k2,1) = chem(it,kt1:kt2,jt,p_h2so4)/factgas
            elseif (p_sulf .ge. p1st)   then
                rsub(kh2so4,k1:k2,1) = chem(it,kt1:kt2,jt,p_sulf)/factgas
            endif

	    if (p_hno3 .ge. p1st)   &
		rsub(khno3,k1:k2,1)  = chem(it,kt1:kt2,jt,p_hno3)/factgas
	    if (p_hcl .ge. p1st)   &
		rsub(khcl,k1:k2,1)   = chem(it,kt1:kt2,jt,p_hcl)/factgas
	    if (p_nh3 .ge. p1st)   &
		rsub(knh3,k1:k2,1)   = chem(it,kt1:kt2,jt,p_nh3)/factgas
	    if (p_n2o5 .ge. p1st)   &
		rsub(kn2o5,k1:k2,1)   = chem(it,kt1:kt2,jt,p_n2o5)/factgas
	    if (p_clno2 .ge. p1st)   &
		rsub(kclno2,k1:k2,1)   = chem(it,kt1:kt2,jt,p_clno2)/factgas






	    if (p_o3 .ge. p1st)   &
		rsub(ko3,k1:k2,1)   = chem(it,kt1:kt2,jt,p_o3)/factgas
	    if (p_so2 .ge. p1st)   &
		rsub(kso2,k1:k2,1)   = chem(it,kt1:kt2,jt,p_so2)/factgas
	    if (p_h2o2 .ge. p1st)   &
		rsub(kh2o2,k1:k2,1)   = chem(it,kt1:kt2,jt,p_h2o2)/factgas
	    if (p_hcho .ge. p1st)   &
		rsub(khcho,k1:k2,1)   = chem(it,kt1:kt2,jt,p_hcho)/factgas
	    if (p_ho .ge. p1st)   &
		rsub(koh,k1:k2,1)   = chem(it,kt1:kt2,jt,p_ho)/factgas
	    if (p_ho2 .ge. p1st)   &
		rsub(kho2,k1:k2,1)   = chem(it,kt1:kt2,jt,p_ho2)/factgas
	    if (p_no3 .ge. p1st)   &
		rsub(kno3,k1:k2,1)   = chem(it,kt1:kt2,jt,p_no3)/factgas
	    if (p_no .ge. p1st)   &
		rsub(kno,k1:k2,1)   = chem(it,kt1:kt2,jt,p_no)/factgas
	    if (p_no2 .ge. p1st)   &
		rsub(kno2,k1:k2,1)   = chem(it,kt1:kt2,jt,p_no2)/factgas
	    if (p_hono .ge. p1st)   &
		rsub(khono,k1:k2,1)   = chem(it,kt1:kt2,jt,p_hono)/factgas
	    if (p_pan .ge. p1st)   &
		rsub(kpan,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pan)/factgas
            if (p_pcg1_b_c .ge. p1st)   &
                rsub(kpcg1_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg1_b_c)/factgas
            if (p_pcg2_b_c .ge. p1st)   &
                rsub(kpcg2_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg2_b_c)/factgas
            if (p_pcg3_b_c .ge. p1st)   &
                rsub(kpcg3_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg3_b_c)/factgas
            if (p_pcg4_b_c .ge. p1st)   &
                rsub(kpcg4_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg4_b_c)/factgas
            if (p_pcg5_b_c .ge. p1st)   &
                rsub(kpcg5_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg5_b_c)/factgas
            if (p_pcg6_b_c .ge. p1st)   &
                rsub(kpcg6_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg6_b_c)/factgas
            if (p_pcg7_b_c .ge. p1st)   &
                rsub(kpcg7_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg7_b_c)/factgas
            if (p_pcg8_b_c .ge. p1st)   &
                rsub(kpcg8_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg8_b_c)/factgas
            if (p_pcg9_b_c .ge. p1st)   &
                rsub(kpcg9_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg9_b_c)/factgas
            if (p_pcg1_b_o .ge. p1st)   &
                rsub(kpcg1_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg1_b_o)/factgas
            if (p_pcg2_b_o .ge. p1st)   &
                rsub(kpcg2_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg2_b_o)/factgas
            if (p_pcg3_b_o .ge. p1st)   &
                rsub(kpcg3_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg3_b_o)/factgas
            if (p_pcg4_b_o .ge. p1st)   &
                rsub(kpcg4_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg4_b_o)/factgas
            if (p_pcg5_b_o .ge. p1st)   &
                rsub(kpcg5_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg5_b_o)/factgas
            if (p_pcg6_b_o .ge. p1st)   &
                rsub(kpcg6_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg6_b_o)/factgas
            if (p_pcg7_b_o .ge. p1st)   &
                rsub(kpcg7_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg7_b_o)/factgas
            if (p_pcg8_b_o .ge. p1st)   &
                rsub(kpcg8_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg8_b_o)/factgas
            if (p_pcg9_b_o .ge. p1st)   &
                rsub(kpcg9_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg9_b_o)/factgas
            if (p_opcg1_b_c .ge. p1st)   &
                rsub(kopcg1_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg1_b_c)/factgas
            if (p_opcg2_b_c .ge. p1st)   &
                rsub(kopcg2_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg2_b_c)/factgas
            if (p_opcg3_b_c .ge. p1st)   &
                rsub(kopcg3_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg3_b_c)/factgas
            if (p_opcg4_b_c .ge. p1st)   &
                rsub(kopcg4_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg4_b_c)/factgas
            if (p_opcg5_b_c .ge. p1st)   &
                rsub(kopcg5_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg5_b_c)/factgas
            if (p_opcg6_b_c .ge. p1st)   &
                rsub(kopcg6_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg6_b_c)/factgas
            if (p_opcg7_b_c .ge. p1st)   &
                rsub(kopcg7_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg7_b_c)/factgas
            if (p_opcg8_b_c .ge. p1st)   &
                rsub(kopcg8_b_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg8_b_c)/factgas
            if (p_opcg1_b_o .ge. p1st)   &
                rsub(kopcg1_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg1_b_o)/factgas
            if (p_opcg2_b_o .ge. p1st)   &
                rsub(kopcg2_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg2_b_o)/factgas
            if (p_opcg3_b_o .ge. p1st)   &
                rsub(kopcg3_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg3_b_o)/factgas
            if (p_opcg4_b_o .ge. p1st)   &
                rsub(kopcg4_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg4_b_o)/factgas
            if (p_opcg5_b_o .ge. p1st)   &
                rsub(kopcg5_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg5_b_o)/factgas
            if (p_opcg6_b_o .ge. p1st)   &
                rsub(kopcg6_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg6_b_o)/factgas
            if (p_opcg7_b_o .ge. p1st)   &
                rsub(kopcg7_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg7_b_o)/factgas
            if (p_opcg8_b_o .ge. p1st)   &
                rsub(kopcg8_b_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg8_b_o)/factgas
            if (p_pcg1_f_c .ge. p1st)   &
                rsub(kpcg1_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg1_f_c)/factgas
            if (p_pcg2_f_c .ge. p1st)   &
                rsub(kpcg2_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg2_f_c)/factgas
            if (p_pcg3_f_c .ge. p1st)   &
                rsub(kpcg3_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg3_f_c)/factgas
            if (p_pcg4_f_c .ge. p1st)   &
                rsub(kpcg4_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg4_f_c)/factgas
            if (p_pcg5_f_c .ge. p1st)   &
                rsub(kpcg5_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg5_f_c)/factgas
            if (p_pcg6_f_c .ge. p1st)   &
                rsub(kpcg6_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg6_f_c)/factgas
            if (p_pcg7_f_c .ge. p1st)   &
                rsub(kpcg7_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg7_f_c)/factgas
            if (p_pcg8_f_c .ge. p1st)   &
                rsub(kpcg8_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg8_f_c)/factgas
            if (p_pcg9_f_c .ge. p1st)   &
                rsub(kpcg9_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg9_f_c)/factgas
            if (p_pcg1_f_o .ge. p1st)   &
                rsub(kpcg1_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg1_f_o)/factgas
            if (p_pcg2_f_o .ge. p1st)   &
                rsub(kpcg2_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg2_f_o)/factgas
            if (p_pcg3_f_o .ge. p1st)   &
                rsub(kpcg3_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg3_f_o)/factgas
            if (p_pcg4_f_o .ge. p1st)   &
                rsub(kpcg4_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg4_f_o)/factgas
            if (p_pcg5_f_o .ge. p1st)   &
                rsub(kpcg5_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg5_f_o)/factgas
            if (p_pcg6_f_o .ge. p1st)   &
                rsub(kpcg6_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg6_f_o)/factgas
            if (p_pcg7_f_o .ge. p1st)   &
                rsub(kpcg7_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg7_f_o)/factgas
            if (p_pcg8_f_o .ge. p1st)   &
                rsub(kpcg8_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg8_f_o)/factgas
            if (p_pcg9_f_o .ge. p1st)   &
                rsub(kpcg9_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_pcg9_f_o)/factgas
            if (p_opcg1_f_c .ge. p1st)   &
                rsub(kopcg1_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg1_f_c)/factgas
            if (p_opcg2_f_c .ge. p1st)   &
                rsub(kopcg2_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg2_f_c)/factgas
            if (p_opcg3_f_c .ge. p1st)   &
                rsub(kopcg3_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg3_f_c)/factgas
            if (p_opcg4_f_c .ge. p1st)   &
                rsub(kopcg4_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg4_f_c)/factgas
            if (p_opcg5_f_c .ge. p1st)   &
                rsub(kopcg5_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg5_f_c)/factgas
            if (p_opcg6_f_c .ge. p1st)   &
                rsub(kopcg6_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg6_f_c)/factgas
            if (p_opcg7_f_c .ge. p1st)   &
                rsub(kopcg7_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg7_f_c)/factgas
            if (p_opcg8_f_c .ge. p1st)   &
                rsub(kopcg8_f_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg8_f_c)/factgas
            if (p_opcg1_f_o .ge. p1st)   &
                rsub(kopcg1_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg1_f_o)/factgas
            if (p_opcg2_f_o .ge. p1st)   &
                rsub(kopcg2_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg2_f_o)/factgas
            if (p_opcg3_f_o .ge. p1st)   &
                rsub(kopcg3_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg3_f_o)/factgas
            if (p_opcg4_f_o .ge. p1st)   &
                rsub(kopcg4_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg4_f_o)/factgas
            if (p_opcg5_f_o .ge. p1st)   &
                rsub(kopcg5_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg5_f_o)/factgas
            if (p_opcg6_f_o .ge. p1st)   &
                rsub(kopcg6_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg6_f_o)/factgas
            if (p_opcg7_f_o .ge. p1st)   &
                rsub(kopcg7_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg7_f_o)/factgas
            if (p_opcg8_f_o .ge. p1st)   &
                rsub(kopcg8_f_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_opcg8_f_o)/factgas
            if (p_smpa .ge. p1st)   &
                rsub(ksmpa,k1:k2,1)   = chem(it,kt1:kt2,jt,p_smpa)/factgas
            if (p_smpbb .ge. p1st)   &
                rsub(ksmpbb,k1:k2,1)   = chem(it,kt1:kt2,jt,p_smpbb)/factgas
            if (p_iepox .ge. p1st)   &
                rsub(kiepox,k1:k2,1)   = chem(it,kt1:kt2,jt,p_iepox)/factgas
            if (p_tetrol .ge. p1st)   &
                rsub(ktetrol,k1:k2,1)   = chem(it,kt1:kt2,jt,p_tetrol)/factgas
            if (p_cldtetrol .ge. p1st)   &
                rsub(kcldtetrol,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cldtetrol)/factgas
            if (p_gly .ge. p1st)   &
                rsub(kgly,k1:k2,1)   = chem(it,kt1:kt2,jt,p_gly)/factgas
            if (p_ant1_c .ge. p1st)   &
                rsub(kant1_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_ant1_c)/factgas
            if (p_ant2_c .ge. p1st)   &
                rsub(kant2_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_ant2_c)/factgas
            if (p_ant3_c .ge. p1st)   &
                rsub(kant3_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_ant3_c)/factgas
            if (p_ant4_c .ge. p1st)   &
                rsub(kant4_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_ant4_c)/factgas
            if (p_ant1_o .ge. p1st)   &
                rsub(kant1_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_ant1_o)/factgas
            if (p_ant2_o .ge. p1st)   &
                rsub(kant2_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_ant2_o)/factgas
            if (p_ant3_o .ge. p1st)   &
                rsub(kant3_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_ant3_o)/factgas
            if (p_ant4_o .ge. p1st)   &
                rsub(kant4_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_ant4_o)/factgas
            if (p_biog1_c .ge. p1st)   &
                rsub(kbiog1_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_biog1_c)/factgas
            if (p_biog2_c .ge. p1st)   &
                rsub(kbiog2_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_biog2_c)/factgas
            if (p_biog3_c .ge. p1st)   &
                rsub(kbiog3_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_biog3_c)/factgas
            if (p_biog4_c .ge. p1st)   &
                rsub(kbiog4_c,k1:k2,1)   = chem(it,kt1:kt2,jt,p_biog4_c)/factgas
            if (p_biog1_o .ge. p1st)   &
                rsub(kbiog1_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_biog1_o)/factgas
            if (p_biog2_o .ge. p1st)   &
                rsub(kbiog2_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_biog2_o)/factgas
            if (p_biog3_o .ge. p1st)   &
                rsub(kbiog3_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_biog3_o)/factgas
            if (p_biog4_o .ge. p1st)   &
                rsub(kbiog4_o,k1:k2,1)   = chem(it,kt1:kt2,jt,p_biog4_o)/factgas
            if (p_cvasoaX .ge. p1st)   &
                rsub(kasoaX,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cvasoaX)/factgas
            if (p_cvasoa1 .ge. p1st)   &
                rsub(kasoa1,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cvasoa1)/factgas
            if (p_cvasoa2 .ge. p1st)   &
                rsub(kasoa2,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cvasoa2)/factgas
            if (p_cvasoa3 .ge. p1st)   &
                rsub(kasoa3,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cvasoa3)/factgas
            if (p_cvasoa4 .ge. p1st)   &
                rsub(kasoa4,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cvasoa4)/factgas
            if (p_cvbsoaX .ge. p1st)   &
                rsub(kbsoaX,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cvbsoaX)/factgas
            if (p_cvbsoa1 .ge. p1st)   &
                rsub(kbsoa1,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cvbsoa1)/factgas
            if (p_cvbsoa2 .ge. p1st)   &
                rsub(kbsoa2,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cvbsoa2)/factgas
            if (p_cvbsoa3 .ge. p1st)   &
                rsub(kbsoa3,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cvbsoa3)/factgas
            if (p_cvbsoa4 .ge. p1st)   &
                rsub(kbsoa4,k1:k2,1)   = chem(it,kt1:kt2,jt,p_cvbsoa4)/factgas

	    
	    if (p_TERP_CN3_g .ge. p1st)   &
	        rsub(kterp_g_Cn3,k1:k2,1)   = chem(it,kt1:kt2,jt,p_TERP_CN3_g)/factgas 


                    if(present (phw_bin1)) &
            phw_bin1(it,kt1:kt2,jt)=phvar1(k1:k2,1)
                    if(present (phw_bin2)) &
            phw_bin2(it,kt1:kt2,jt)=phvar2(k1:k2,1)
                    if(present (phw_bin3)) &
            phw_bin3(it,kt1:kt2,jt)=phvar3(k1:k2,1)
                    if(present (phw_bin4)) &
            phw_bin4(it,kt1:kt2,jt)=phvar4(k1:k2,1)

               if(present (gammaaq)) &
            gammaaq(it,kt1:kt2,jt)=gammaiepoxaq(k1:k2,1)

               if(present (gammacld)) &
            gammacld(it,kt1:kt2,jt)=gammaiepoxcld(k1:k2,1)

               if(present (clddropdia)) &
            clddropdia(it,kt1:kt2,jt)=effdropdia(k1:k2,1)


               if(present (betaaq)) &
            betaaq(it,kt1:kt2,jt)=betaiepoxaq(k1:k2,1)

               if(present (betacld)) &
            betacld(it,kt1:kt2,jt)=betaiepoxcld(k1:k2,1)

               if(present (coatbin1)) &
            coatbin1(it,kt1:kt2,jt)=iepoxcoatbin1(k1:k2,1)
               if(present (coatbin2)) &
            coatbin2(it,kt1:kt2,jt)=iepoxcoatbin2(k1:k2,1)
               if(present (coatbin3)) &
            coatbin3(it,kt1:kt2,jt)=iepoxcoatbin3(k1:k2,1)
               if(present (coatbin4)) &
            coatbin4(it,kt1:kt2,jt)=iepoxcoatbin4(k1:k2,1)


	    do iphase=1,nphase_aer
	    do itype=1,ntype_aer
	    do n = 1, nsize_aer(itype)
		rsub(lptr_so4_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_so4_aer(n,itype,iphase))/factaerso4
		rsub(numptr_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,numptr_aer(n,itype,iphase))/factaernum

		if (lptr_no3_aer(n,itype,iphase) .ge. p1st)   &
		    rsub(lptr_no3_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_no3_aer(n,itype,iphase))/factaerno3
		if (lptr_cl_aer(n,itype,iphase) .ge. p1st)   &
		    rsub(lptr_cl_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_cl_aer(n,itype,iphase))/factaercl
		if (lptr_msa_aer(n,itype,iphase) .ge. p1st)   &
		    rsub(lptr_msa_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_msa_aer(n,itype,iphase))/factaermsa
		if (lptr_co3_aer(n,itype,iphase) .ge. p1st)   &
		    rsub(lptr_co3_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_co3_aer(n,itype,iphase))/factaerco3
		if (lptr_nh4_aer(n,itype,iphase) .ge. p1st)   &
		    rsub(lptr_nh4_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_nh4_aer(n,itype,iphase))/factaernh4
		if (lptr_na_aer(n,itype,iphase) .ge. p1st)   &
		    rsub(lptr_na_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_na_aer(n,itype,iphase))/factaerna
		if (lptr_ca_aer(n,itype,iphase) .ge. p1st)   &
		    rsub(lptr_ca_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_ca_aer(n,itype,iphase))/factaerca
		if (lptr_oin_aer(n,itype,iphase) .ge. p1st)   &
		    rsub(lptr_oin_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_oin_aer(n,itype,iphase))/factaeroin
		if (lptr_oc_aer(n,itype,iphase) .ge. p1st)   &
		    rsub(lptr_oc_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_oc_aer(n,itype,iphase))/factaeroc
		if (lptr_bc_aer(n,itype,iphase) .ge. p1st)   &
		    rsub(lptr_bc_aer(n,itype,iphase),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,lptr_bc_aer(n,itype,iphase))/factaerbc
		if (hyswptr_aer(n,itype) .ge. p1st)   &
		    rsub(hyswptr_aer(n,itype),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,hyswptr_aer(n,itype))/factaerhysw
		if (waterptr_aer(n,itype) .ge. p1st)   &
		    rsub(waterptr_aer(n,itype),k1:k2,1) =   &
		    chem(it,kt1:kt2,jt,waterptr_aer(n,itype))/factaerwater
                if (lptr_tanv_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_tanv_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_tanv_aer(n,itype,iphase))/factaertanv
                if (lptr_sopnv_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_sopnv_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_sopnv_aer(n,itype,iphase))/factaersopnv
                if (lptr_ternv_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_ternv_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_ternv_aer(n,itype,iphase))/factaerternv
                if (lptr_seqnv_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_seqnv_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_seqnv_aer(n,itype,iphase))/factaerseqnv
                if (lptr_sianv_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_sianv_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_sianv_aer(n,itype,iphase))/factaersianv
                if (lptr_pcg1_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg1_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg1_b_c_aer(n,itype,iphase))/factaerpcg1_b_c
                if (lptr_pcg2_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg2_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg2_b_c_aer(n,itype,iphase))/factaerpcg2_b_c
                if (lptr_pcg3_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg3_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg3_b_c_aer(n,itype,iphase))/factaerpcg3_b_c
                if (lptr_pcg4_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg4_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg4_b_c_aer(n,itype,iphase))/factaerpcg4_b_c
                if (lptr_pcg5_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg5_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg5_b_c_aer(n,itype,iphase))/factaerpcg5_b_c
                if (lptr_pcg6_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg6_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg6_b_c_aer(n,itype,iphase))/factaerpcg6_b_c
                if (lptr_pcg7_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg7_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg7_b_c_aer(n,itype,iphase))/factaerpcg7_b_c
                if (lptr_pcg8_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg8_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg8_b_c_aer(n,itype,iphase))/factaerpcg8_b_c
                if (lptr_pcg9_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg9_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg9_b_c_aer(n,itype,iphase))/factaerpcg9_b_c
                if (lptr_pcg1_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg1_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg1_b_o_aer(n,itype,iphase))/factaerpcg1_b_o
                if (lptr_pcg2_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg2_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg2_b_o_aer(n,itype,iphase))/factaerpcg2_b_o
                if (lptr_pcg3_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg3_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg3_b_o_aer(n,itype,iphase))/factaerpcg3_b_o
                if (lptr_pcg4_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg4_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg4_b_o_aer(n,itype,iphase))/factaerpcg4_b_o
                if (lptr_pcg5_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg5_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg5_b_o_aer(n,itype,iphase))/factaerpcg5_b_o
                if (lptr_pcg6_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg6_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg6_b_o_aer(n,itype,iphase))/factaerpcg6_b_o
                if (lptr_pcg7_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg7_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg7_b_o_aer(n,itype,iphase))/factaerpcg7_b_o
                if (lptr_pcg8_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg8_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg8_b_o_aer(n,itype,iphase))/factaerpcg8_b_o
                if (lptr_pcg9_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg9_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg9_b_o_aer(n,itype,iphase))/factaerpcg9_b_o
                if (lptr_opcg1_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg1_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg1_b_c_aer(n,itype,iphase))/factaeropcg1_b_c
                if (lptr_opcg2_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg2_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg2_b_c_aer(n,itype,iphase))/factaeropcg2_b_c
                if (lptr_opcg3_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg3_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg3_b_c_aer(n,itype,iphase))/factaeropcg3_b_c
                if (lptr_opcg4_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg4_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg4_b_c_aer(n,itype,iphase))/factaeropcg4_b_c
                if (lptr_opcg5_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg5_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg5_b_c_aer(n,itype,iphase))/factaeropcg5_b_c
                if (lptr_opcg6_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg6_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg6_b_c_aer(n,itype,iphase))/factaeropcg6_b_c
                if (lptr_opcg7_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg7_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg7_b_c_aer(n,itype,iphase))/factaeropcg7_b_c
                if (lptr_opcg8_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg8_b_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg8_b_c_aer(n,itype,iphase))/factaeropcg8_b_c
                if (lptr_opcg1_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg1_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg1_b_o_aer(n,itype,iphase))/factaeropcg1_b_o
                if (lptr_opcg2_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg2_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg2_b_o_aer(n,itype,iphase))/factaeropcg2_b_o
                if (lptr_opcg3_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg3_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg3_b_o_aer(n,itype,iphase))/factaeropcg3_b_o
                if (lptr_opcg4_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg4_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg4_b_o_aer(n,itype,iphase))/factaeropcg4_b_o
                if (lptr_opcg5_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg5_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg5_b_o_aer(n,itype,iphase))/factaeropcg5_b_o
                if (lptr_opcg6_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg6_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg6_b_o_aer(n,itype,iphase))/factaeropcg6_b_o
                if (lptr_opcg7_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg7_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg7_b_o_aer(n,itype,iphase))/factaeropcg7_b_o
                if (lptr_opcg8_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg8_b_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg8_b_o_aer(n,itype,iphase))/factaeropcg8_b_o
                if (lptr_pcg1_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg1_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg1_f_c_aer(n,itype,iphase))/factaerpcg1_f_c
                if (lptr_pcg2_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg2_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg2_f_c_aer(n,itype,iphase))/factaerpcg2_f_c
                if (lptr_pcg3_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg3_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg3_f_c_aer(n,itype,iphase))/factaerpcg3_f_c
                if (lptr_pcg4_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg4_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg4_f_c_aer(n,itype,iphase))/factaerpcg4_f_c
                if (lptr_pcg5_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg5_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg5_f_c_aer(n,itype,iphase))/factaerpcg5_f_c
                if (lptr_pcg6_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg6_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg6_f_c_aer(n,itype,iphase))/factaerpcg6_f_c
                if (lptr_pcg7_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg7_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg7_f_c_aer(n,itype,iphase))/factaerpcg7_f_c
                if (lptr_pcg8_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg8_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg8_f_c_aer(n,itype,iphase))/factaerpcg8_f_c
                if (lptr_pcg9_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg9_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg9_f_c_aer(n,itype,iphase))/factaerpcg9_f_c
                if (lptr_pcg1_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg1_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg1_f_o_aer(n,itype,iphase))/factaerpcg1_f_o
                if (lptr_pcg2_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg2_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg2_f_o_aer(n,itype,iphase))/factaerpcg2_f_o
                if (lptr_pcg3_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg3_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg3_f_o_aer(n,itype,iphase))/factaerpcg3_f_o
                if (lptr_pcg4_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg4_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg4_f_o_aer(n,itype,iphase))/factaerpcg4_f_o
                if (lptr_pcg5_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg5_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg5_f_o_aer(n,itype,iphase))/factaerpcg5_f_o
                if (lptr_pcg6_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg6_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg6_f_o_aer(n,itype,iphase))/factaerpcg6_f_o
                if (lptr_pcg7_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg7_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg7_f_o_aer(n,itype,iphase))/factaerpcg7_f_o
                if (lptr_pcg8_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg8_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg8_f_o_aer(n,itype,iphase))/factaerpcg8_f_o
                if (lptr_pcg9_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_pcg9_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_pcg9_f_o_aer(n,itype,iphase))/factaerpcg9_f_o
                if (lptr_opcg1_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg1_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg1_f_c_aer(n,itype,iphase))/factaeropcg1_f_c
                if (lptr_opcg2_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg2_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg2_f_c_aer(n,itype,iphase))/factaeropcg2_f_c
                if (lptr_opcg3_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg3_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg3_f_c_aer(n,itype,iphase))/factaeropcg3_f_c
                if (lptr_opcg4_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg4_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg4_f_c_aer(n,itype,iphase))/factaeropcg4_f_c
                if (lptr_opcg5_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg5_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg5_f_c_aer(n,itype,iphase))/factaeropcg5_f_c
                if (lptr_opcg6_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg6_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg6_f_c_aer(n,itype,iphase))/factaeropcg6_f_c
                if (lptr_opcg7_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg7_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg7_f_c_aer(n,itype,iphase))/factaeropcg7_f_c
                if (lptr_opcg8_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg8_f_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg8_f_c_aer(n,itype,iphase))/factaeropcg8_f_c
                if (lptr_opcg1_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg1_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg1_f_o_aer(n,itype,iphase))/factaeropcg1_f_o
                if (lptr_opcg2_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg2_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg2_f_o_aer(n,itype,iphase))/factaeropcg2_f_o
                if (lptr_opcg3_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg3_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg3_f_o_aer(n,itype,iphase))/factaeropcg3_f_o
                if (lptr_opcg4_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg4_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg4_f_o_aer(n,itype,iphase))/factaeropcg4_f_o
                if (lptr_opcg5_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg5_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg5_f_o_aer(n,itype,iphase))/factaeropcg5_f_o
                if (lptr_opcg6_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg6_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg6_f_o_aer(n,itype,iphase))/factaeropcg6_f_o
                if (lptr_opcg7_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg7_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg7_f_o_aer(n,itype,iphase))/factaeropcg7_f_o
                if (lptr_opcg8_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_opcg8_f_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_opcg8_f_o_aer(n,itype,iphase))/factaeropcg8_f_o
                if (lptr_smpa_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_smpa_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_smpa_aer(n,itype,iphase))/factaersmpa
                if (lptr_smpbb_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_smpbb_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_smpbb_aer(n,itype,iphase))/factaersmpbb
                if (lptr_iepox_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_iepox_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_iepox_aer(n,itype,iphase))/factaeriepox
                if (lptr_iepoxos_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_iepoxos_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_iepoxos_aer(n,itype,iphase))/factaeriepoxos
                if (lptr_tetrol_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_tetrol_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_tetrol_aer(n,itype,iphase))/factaertetrol
                if (lptr_cldiepoxos_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_cldiepoxos_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_cldiepoxos_aer(n,itype,iphase))/factaercldiepoxos
                if (lptr_cldtetrol_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_cldtetrol_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_cldtetrol_aer(n,itype,iphase))/factaercldtetrol
                if (lptr_gly_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_gly_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_gly_aer(n,itype,iphase))/factaergly
                if (lptr_glysoa_r1_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_glysoa_r1_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_glysoa_r1_aer(n,itype,iphase))/factaerglyr1
                if (lptr_glysoa_r2_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_glysoa_r2_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_glysoa_r2_aer(n,itype,iphase))/factaerglyr2
                if (lptr_glysoa_sfc_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_glysoa_sfc_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_glysoa_sfc_aer(n,itype,iphase))/factaerglysfc
                if (lptr_glysoa_oh_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_glysoa_oh_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_glysoa_oh_aer(n,itype,iphase))/factaerglyoh
                if (lptr_glysoa_nh4_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_glysoa_nh4_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_glysoa_nh4_aer(n,itype,iphase))/factaerglynh4
                if (lptr_ant1_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_ant1_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_ant1_c_aer(n,itype,iphase))/factaerant1_c
                if (lptr_ant2_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_ant2_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_ant2_c_aer(n,itype,iphase))/factaerant2_c
                if (lptr_ant3_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_ant3_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_ant3_c_aer(n,itype,iphase))/factaerant3_c
                if (lptr_ant4_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_ant4_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_ant4_c_aer(n,itype,iphase))/factaerant4_c
                if (lptr_ant1_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_ant1_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_ant1_o_aer(n,itype,iphase))/factaerant1_o
                if (lptr_ant2_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_ant2_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_ant2_o_aer(n,itype,iphase))/factaerant2_o
                if (lptr_ant3_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_ant3_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_ant3_o_aer(n,itype,iphase))/factaerant3_o
                if (lptr_ant4_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_ant4_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_ant4_o_aer(n,itype,iphase))/factaerant4_o
                if (lptr_biog1_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_biog1_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_biog1_c_aer(n,itype,iphase))/factaerbiog1_c
                if (lptr_biog2_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_biog2_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_biog2_c_aer(n,itype,iphase))/factaerbiog2_c
                if (lptr_biog3_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_biog3_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_biog3_c_aer(n,itype,iphase))/factaerbiog3_c
                if (lptr_biog4_c_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_biog4_c_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_biog4_c_aer(n,itype,iphase))/factaerbiog4_c
                if (lptr_biog1_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_biog1_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_biog1_o_aer(n,itype,iphase))/factaerbiog1_o
                if (lptr_biog2_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_biog2_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_biog2_o_aer(n,itype,iphase))/factaerbiog2_o
                if (lptr_biog3_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_biog3_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_biog3_o_aer(n,itype,iphase))/factaerbiog3_o
                if (lptr_biog4_o_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_biog4_o_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_biog4_o_aer(n,itype,iphase))/factaerbiog4_o
                if (lptr_asoaX_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_asoaX_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_asoaX_aer(n,itype,iphase))/factaerasoaX
                if (lptr_asoa1_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_asoa1_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_asoa1_aer(n,itype,iphase))/factaerasoa1
                if (lptr_asoa2_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_asoa2_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_asoa2_aer(n,itype,iphase))/factaerasoa2
                if (lptr_asoa3_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_asoa3_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_asoa3_aer(n,itype,iphase))/factaerasoa3
                if (lptr_asoa4_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_asoa4_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_asoa4_aer(n,itype,iphase))/factaerasoa4
                if (lptr_bsoaX_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_bsoaX_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_bsoaX_aer(n,itype,iphase))/factaerbsoaX
                if (lptr_bsoa1_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_bsoa1_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_bsoa1_aer(n,itype,iphase))/factaerbsoa1
                if (lptr_bsoa2_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_bsoa2_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_bsoa2_aer(n,itype,iphase))/factaerbsoa2
                if (lptr_bsoa3_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_bsoa3_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_bsoa3_aer(n,itype,iphase))/factaerbsoa3
                if (lptr_bsoa4_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_bsoa4_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_bsoa4_aer(n,itype,iphase))/factaerbsoa4

                
		
		if (lptr_terp_Cn3_aer(n,itype,iphase) .ge. p1st)   &
                    rsub(lptr_terp_Cn3_aer(n,itype,iphase),k1:k2,1) =   &
                    chem(it,kt1:kt2,jt,lptr_terp_Cn3_aer(n,itype,iphase))/factaerterp_Cn3

	    end do 
	    end do 
	    end do 
	end if


	afracsubarea(k1:k2,1) = 1.0
	rsub(ktemp,k1:k2,1) = t_phy(it,kt1:kt2,jt)
	rsub(kh2o,k1:k2,1) = moist(it,kt1:kt2,jt,p_qv)/factmoist
	cairclm(k1:k2) = rho_phy(it,kt1:kt2,jt)/factdens
	ptotclm(k1:k2) = p_phy(it,kt1:kt2,jt)/factpres
        if(present (swdown)) &
        swdown_cell=swdown(it,jt)

	if (p_qc .ge. p1st)   &
	    rcldwtr_sub(k1:k2,1) = moist(it,kt1:kt2,jt,p_qc)/factmoist

        if(present (ph_no2)) &
        phno2_box(k1:k2) =ph_no2(it,kt1:kt2,jt)











	if ((aboxtest_rh_method .gt. 0) .and.   &
	    (aboxtest_rh_method .ne. 2)) then
	    do kt = ktmaps, ktmape
		k = kt + k_pegshift
		onemeps = 1.0 - 0.622
		dumesat = esat_gchm( rsub(ktemp,k,1) )
		dumrsat = dumesat / (ptotclm(k) - onemeps*dumesat)
		dumrelhum = rsub(kh2o,k,1) / max( dumrsat, 1.e-20 )

		dumrelhum = max( 0.0, min( 0.99, dumrelhum ) )

		if (aboxtest_rh_method .eq. 3) then




		    continue
		else
		    relhumclm(k) = dumrelhum
		end if
		relhumclm(k) = max( relhumclm(k), aboxtest_min_relhum )
		relhumclm(k) = min( relhumclm(k), aboxtest_max_relhum )
	    end do
	end if


	do kt = ktmaps, ktmape
	    k = kt + k_pegshift
	    rsub(ktemp,k,1) =   &
		max( rsub(ktemp,k,1), aboxtest_min_temp )
	end do

	return








2000	continue






	if ((aboxtest_map_method .eq. 2) .or.   &
	    (aboxtest_map_method .eq. 3)) then
	    do l = 2, num_chem
		ido_l = 1
		if (aboxtest_gases_fixed .eq. 10) then
		    if ((l .eq. kh2so4  ) .or. (l .eq. khno3  ) .or.   &
		        (l .eq. khcl    ) .or. (l .eq. knh3   ) .or.   &
		        (l .eq. kclno2  ) .or. (l .eq. kn2o5  ) .or.   &
		        (l .eq. ko3     ) .or.                         &
		        (l .eq. kso2    ) .or. (l .eq. kh2o2  ) .or.   &
		        (l .eq. khcho   )  .or.   &
		        (l .eq. koh     ) .or. (l .eq. kho2   ) .or.   &
		        (l .eq. kno3    ) .or. (l .eq. kno    ) .or.   &
		        (l .eq. kno2    ) .or. (l .eq. khono  ) .or.   &
		        (l .eq. kpan    )  .or.   &
                        (l .eq. kpcg1_b_c    ) .or. (l .eq. kpcg2_b_c   ) .or.   &
                        (l .eq. kpcg3_b_c    ) .or. (l .eq. kpcg4_b_c   ) .or.   &
                        (l .eq. kpcg5_b_c    ) .or. (l .eq. kpcg6_b_c   ) .or.   &
                        (l .eq. kpcg7_b_c    ) .or. (l .eq. kpcg8_b_c   ) .or.   &
                        (l .eq. kpcg9_b_c    ) .or. (l .eq. kpcg1_b_o   ) .or.   &
                        (l .eq. kpcg2_b_o    ) .or. (l .eq. kpcg3_b_o   ) .or.   &
                        (l .eq. kpcg4_b_o    ) .or. (l .eq. kpcg5_b_o   ) .or.   &
                        (l .eq. kpcg6_b_o    ) .or. (l .eq. kpcg7_b_o   ) .or.   &
                        (l .eq. kpcg8_b_o    ) .or. (l .eq. kpcg9_b_o   ) .or.   &
                        (l .eq. kopcg1_b_c    ) .or. (l .eq. kopcg2_b_c   ) .or.   &
                        (l .eq. kopcg3_b_c    ) .or. (l .eq. kopcg4_b_c   ) .or.   &
                        (l .eq. kopcg5_b_c    ) .or. (l .eq. kopcg6_b_c   ) .or.   &
                        (l .eq. kopcg7_b_c    ) .or. (l .eq. kopcg8_b_c   ) .or.   &
                        (l .eq. kopcg1_b_o    ) .or. (l .eq. kopcg2_b_o   ) .or.   &
                        (l .eq. kopcg3_b_o    ) .or. (l .eq. kopcg4_b_o   ) .or.   &
                        (l .eq. kopcg5_b_o    ) .or. (l .eq. kopcg6_b_o   ) .or.   &
                        (l .eq. kopcg7_b_o    ) .or. (l .eq. kopcg8_b_o   ) .or.   &
                        (l .eq. kpcg1_f_c    ) .or. (l .eq. kpcg2_f_c   ) .or.   &
                        (l .eq. kpcg3_f_c    ) .or. (l .eq. kpcg4_f_c   ) .or.   &
                        (l .eq. kpcg5_f_c    ) .or. (l .eq. kpcg6_f_c   ) .or.   &
                        (l .eq. kpcg7_f_c    ) .or. (l .eq. kpcg8_f_c   ) .or.   &
                        (l .eq. kpcg9_f_c    ) .or. (l .eq. kpcg1_f_o   ) .or.   &
                        (l .eq. kpcg2_f_o    ) .or. (l .eq. kpcg3_f_o   ) .or.   &
                        (l .eq. kpcg4_f_o    ) .or. (l .eq. kpcg5_f_o   ) .or.   &
                        (l .eq. kpcg6_f_o    ) .or. (l .eq. kpcg7_f_o   ) .or.   &
                        (l .eq. kpcg8_f_o    ) .or. (l .eq. kpcg9_f_o   ) .or.   &
                        (l .eq. kopcg1_f_c    ) .or. (l .eq. kopcg2_f_c   ) .or.   &
                        (l .eq. kopcg3_f_c    ) .or. (l .eq. kopcg4_f_c   ) .or.   &
                        (l .eq. kopcg5_f_c    ) .or. (l .eq. kopcg6_f_c   ) .or.   &
                        (l .eq. kopcg7_f_c    ) .or. (l .eq. kopcg8_f_c   ) .or.   &
                        (l .eq. kopcg1_f_o    ) .or. (l .eq. kopcg2_f_o   ) .or.   &
                        (l .eq. kopcg3_f_o    ) .or. (l .eq. kopcg4_f_o   ) .or.   &
                        (l .eq. kopcg5_f_o    ) .or. (l .eq. kopcg6_f_o   ) .or.   &
                        (l .eq. kopcg7_f_o    ) .or. (l .eq. kopcg8_f_o   ) .or.   &
                        (l .eq. ksmpa    ) .or. (l .eq. ksmpbb   ) .or.         &
                        (l .eq. kiepox    ) .or. (l .eq. kgly    )          .or. &
                        (l .eq. ktetrol    ).or. (l .eq. kcldtetrol    )    .or. &
                        (l .eq. kant1_c    ) .or. (l .eq. kant2_c   ) .or.         &
                        (l .eq. kant3_c    ) .or. (l .eq. kant4_c   ) .or.         &
			
			(l .eq. kterp_g_Cn3  ) .or. &
                        (l .eq. kant1_o    ) .or. (l .eq. kant2_o   ) .or.         &
                        (l .eq. kant3_o    ) .or. (l .eq. kant4_o   ) .or.         &
                        (l .eq. kbiog1_c    ) .or. (l .eq. kbiog2_c   ) .or.         &
                        (l .eq. kbiog3_c    ) .or. (l .eq. kbiog4_c   ) .or.         &
                        (l .eq. kbiog1_o    ) .or. (l .eq. kbiog2_o   ) .or.         &
                        (l .eq. kbiog3_o    ) .or. (l .eq. kbiog4_o   ) .or. &
                        (l .eq. kasoaX    ) .or. &
                        (l .eq. kasoa1    ) .or. (l .eq. kasoa2   ) .or.         &
                        (l .eq. kasoa3    ) .or. (l .eq. kasoa4   ) .or.         &
                        (l .eq. kbsoaX    ) .or. &
                        (l .eq. kbsoa1    ) .or. (l .eq. kbsoa2   ) .or.         &
                        (l .eq. kbsoa3    ) .or. (l .eq. kbsoa4   )) then
			ido_l = 0
		    end if
		end if
		if (ido_l .gt. 0) then
		    chem(it,kt1:kt2,jt,l) = rsub(l,k1:k2,1)*factgas
		end if
	    end do
	end if

	p1st = param_first_scalar
	if (aboxtest_map_method .ne. 2) then
	  if (aboxtest_gases_fixed .ne. 10) then
            if (p_h2so4 .ge. p1st)   then
                chem(it,kt1:kt2,jt,p_h2so4) = rsub(kh2so4,k1:k2,1)*factgas
           elseif (p_sulf .ge. p1st)   then
                chem(it,kt1:kt2,jt,p_sulf) = rsub(kh2so4,k1:k2,1)*factgas
            endif
	    if (p_hno3 .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_hno3)  = rsub(khno3,k1:k2,1)*factgas
	    if (p_hcl .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_hcl)   = rsub(khcl,k1:k2,1)*factgas
	    if (p_nh3 .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_nh3)  = rsub(knh3,k1:k2,1)*factgas
	    if (p_n2o5 .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_n2o5)  = rsub(kn2o5,k1:k2,1)*factgas
	    if (p_clno2 .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_clno2)  = rsub(kclno2,k1:k2,1)*factgas

	    if (p_o3 .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_o3)  = rsub(ko3,k1:k2,1)*factgas
	    if (p_so2 .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_so2)  = rsub(kso2,k1:k2,1)*factgas
	    if (p_h2o2 .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_h2o2)  = rsub(kh2o2,k1:k2,1)*factgas
	    if (p_hcho .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_hcho)  = rsub(khcho,k1:k2,1)*factgas
	    if (p_ho .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_ho)  = rsub(koh,k1:k2,1)*factgas
	    if (p_ho2 .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_ho2)  = rsub(kho2,k1:k2,1)*factgas
	    if (p_no3 .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_no3)  = rsub(kno3,k1:k2,1)*factgas
	    if (p_no .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_no)  = rsub(kno,k1:k2,1)*factgas
	    if (p_no2 .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_no2)  = rsub(kno2,k1:k2,1)*factgas
	    if (p_hono .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_hono)  = rsub(khono,k1:k2,1)*factgas
	    if (p_pan .ge. p1st)   &
		chem(it,kt1:kt2,jt,p_pan)  = rsub(kpan,k1:k2,1)*factgas
            if (p_pcg1_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg1_b_c)  = rsub(kpcg1_b_c,k1:k2,1)*factgas
            if (p_pcg2_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg2_b_c)  = rsub(kpcg2_b_c,k1:k2,1)*factgas
            if (p_pcg3_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg3_b_c)  = rsub(kpcg3_b_c,k1:k2,1)*factgas
            if (p_pcg4_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg4_b_c)  = rsub(kpcg4_b_c,k1:k2,1)*factgas
            if (p_pcg5_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg5_b_c)  = rsub(kpcg5_b_c,k1:k2,1)*factgas
            if (p_pcg6_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg6_b_c)  = rsub(kpcg6_b_c,k1:k2,1)*factgas
            if (p_pcg7_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg7_b_c)  = rsub(kpcg7_b_c,k1:k2,1)*factgas
            if (p_pcg8_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg8_b_c)  = rsub(kpcg8_b_c,k1:k2,1)*factgas
            if (p_pcg9_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg9_b_c)  = rsub(kpcg9_b_c,k1:k2,1)*factgas
            if (p_pcg1_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg1_b_o)  = rsub(kpcg1_b_o,k1:k2,1)*factgas
            if (p_pcg2_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg2_b_o)  = rsub(kpcg2_b_o,k1:k2,1)*factgas
            if (p_pcg3_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg3_b_o)  = rsub(kpcg3_b_o,k1:k2,1)*factgas
            if (p_pcg4_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg4_b_o)  = rsub(kpcg4_b_o,k1:k2,1)*factgas
            if (p_pcg5_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg5_b_o)  = rsub(kpcg5_b_o,k1:k2,1)*factgas
            if (p_pcg6_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg6_b_o)  = rsub(kpcg6_b_o,k1:k2,1)*factgas
            if (p_pcg7_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg7_b_o)  = rsub(kpcg7_b_o,k1:k2,1)*factgas
            if (p_pcg8_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg8_b_o)  = rsub(kpcg8_b_o,k1:k2,1)*factgas
            if (p_pcg9_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg9_b_o)  = rsub(kpcg9_b_o,k1:k2,1)*factgas
            if (p_opcg1_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg1_b_c)  = rsub(kopcg1_b_c,k1:k2,1)*factgas
            if (p_opcg2_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg2_b_c)  = rsub(kopcg2_b_c,k1:k2,1)*factgas
            if (p_opcg3_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg3_b_c)  = rsub(kopcg3_b_c,k1:k2,1)*factgas
            if (p_opcg4_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg4_b_c)  = rsub(kopcg4_b_c,k1:k2,1)*factgas
            if (p_opcg5_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg5_b_c)  = rsub(kopcg5_b_c,k1:k2,1)*factgas
            if (p_opcg6_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg6_b_c)  = rsub(kopcg6_b_c,k1:k2,1)*factgas
            if (p_opcg7_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg7_b_c)  = rsub(kopcg7_b_c,k1:k2,1)*factgas
            if (p_opcg8_b_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg8_b_c)  = rsub(kopcg8_b_c,k1:k2,1)*factgas
            if (p_opcg1_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg1_b_o)  = rsub(kopcg1_b_o,k1:k2,1)*factgas
            if (p_opcg2_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg2_b_o)  = rsub(kopcg2_b_o,k1:k2,1)*factgas
            if (p_opcg3_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg3_b_o)  = rsub(kopcg3_b_o,k1:k2,1)*factgas
            if (p_opcg4_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg4_b_o)  = rsub(kopcg4_b_o,k1:k2,1)*factgas
            if (p_opcg5_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg5_b_o)  = rsub(kopcg5_b_o,k1:k2,1)*factgas
            if (p_opcg6_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg6_b_o)  = rsub(kopcg6_b_o,k1:k2,1)*factgas
            if (p_opcg7_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg7_b_o)  = rsub(kopcg7_b_o,k1:k2,1)*factgas
            if (p_opcg8_b_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg8_b_o)  = rsub(kopcg8_b_o,k1:k2,1)*factgas
            if (p_pcg1_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg1_f_c)  = rsub(kpcg1_f_c,k1:k2,1)*factgas
            if (p_pcg2_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg2_f_c)  = rsub(kpcg2_f_c,k1:k2,1)*factgas
            if (p_pcg3_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg3_f_c)  = rsub(kpcg3_f_c,k1:k2,1)*factgas
            if (p_pcg4_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg4_f_c)  = rsub(kpcg4_f_c,k1:k2,1)*factgas
            if (p_pcg5_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg5_f_c)  = rsub(kpcg5_f_c,k1:k2,1)*factgas
            if (p_pcg6_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg6_f_c)  = rsub(kpcg6_f_c,k1:k2,1)*factgas
            if (p_pcg7_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg7_f_c)  = rsub(kpcg7_f_c,k1:k2,1)*factgas
            if (p_pcg8_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg8_f_c)  = rsub(kpcg8_f_c,k1:k2,1)*factgas
            if (p_pcg9_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg9_f_c)  = rsub(kpcg9_f_c,k1:k2,1)*factgas
            if (p_pcg1_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg1_f_o)  = rsub(kpcg1_f_o,k1:k2,1)*factgas
            if (p_pcg2_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg2_f_o)  = rsub(kpcg2_f_o,k1:k2,1)*factgas
            if (p_pcg3_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg3_f_o)  = rsub(kpcg3_f_o,k1:k2,1)*factgas
            if (p_pcg4_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg4_f_o)  = rsub(kpcg4_f_o,k1:k2,1)*factgas
            if (p_pcg5_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg5_f_o)  = rsub(kpcg5_f_o,k1:k2,1)*factgas
            if (p_pcg6_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg6_f_o)  = rsub(kpcg6_f_o,k1:k2,1)*factgas
            if (p_pcg7_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg7_f_o)  = rsub(kpcg7_f_o,k1:k2,1)*factgas
            if (p_pcg8_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg8_f_o)  = rsub(kpcg8_f_o,k1:k2,1)*factgas
            if (p_pcg9_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_pcg9_f_o)  = rsub(kpcg9_f_o,k1:k2,1)*factgas
            if (p_opcg1_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg1_f_c)  = rsub(kopcg1_f_c,k1:k2,1)*factgas
            if (p_opcg2_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg2_f_c)  = rsub(kopcg2_f_c,k1:k2,1)*factgas
            if (p_opcg3_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg3_f_c)  = rsub(kopcg3_f_c,k1:k2,1)*factgas
            if (p_opcg4_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg4_f_c)  = rsub(kopcg4_f_c,k1:k2,1)*factgas
            if (p_opcg5_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg5_f_c)  = rsub(kopcg5_f_c,k1:k2,1)*factgas
            if (p_opcg6_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg6_f_c)  = rsub(kopcg6_f_c,k1:k2,1)*factgas
            if (p_opcg7_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg7_f_c)  = rsub(kopcg7_f_c,k1:k2,1)*factgas
            if (p_opcg8_f_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg8_f_c)  = rsub(kopcg8_f_c,k1:k2,1)*factgas
            if (p_opcg1_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg1_f_o)  = rsub(kopcg1_f_o,k1:k2,1)*factgas
            if (p_opcg2_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg2_f_o)  = rsub(kopcg2_f_o,k1:k2,1)*factgas
            if (p_opcg3_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg3_f_o)  = rsub(kopcg3_f_o,k1:k2,1)*factgas
            if (p_opcg4_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg4_f_o)  = rsub(kopcg4_f_o,k1:k2,1)*factgas
            if (p_opcg5_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg5_f_o)  = rsub(kopcg5_f_o,k1:k2,1)*factgas
            if (p_opcg6_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg6_f_o)  = rsub(kopcg6_f_o,k1:k2,1)*factgas
            if (p_opcg7_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg7_f_o)  = rsub(kopcg7_f_o,k1:k2,1)*factgas
            if (p_opcg8_f_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_opcg8_f_o)  = rsub(kopcg8_f_o,k1:k2,1)*factgas
            if (p_smpa .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_smpa)  = rsub(ksmpa,k1:k2,1)*factgas
            if (p_smpbb .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_smpbb)  = rsub(ksmpbb,k1:k2,1)*factgas
            if (p_iepox .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_iepox)  = rsub(kiepox,k1:k2,1)*factgas
            if (p_tetrol .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_tetrol)  = rsub(ktetrol,k1:k2,1)*factgas
            if (p_cldtetrol .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cldtetrol)  = rsub(kcldtetrol,k1:k2,1)*factgas
            if (p_gly .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_gly)  = rsub(kgly,k1:k2,1)*factgas
            if (p_ant1_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_ant1_c)  = rsub(kant1_c,k1:k2,1)*factgas
            if (p_ant2_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_ant2_c)  = rsub(kant2_c,k1:k2,1)*factgas
            if (p_ant3_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_ant3_c)  = rsub(kant3_c,k1:k2,1)*factgas
            if (p_ant4_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_ant4_c)  = rsub(kant4_c,k1:k2,1)*factgas
            if (p_ant1_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_ant1_o)  = rsub(kant1_o,k1:k2,1)*factgas
            if (p_ant2_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_ant2_o)  = rsub(kant2_o,k1:k2,1)*factgas
            if (p_ant3_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_ant3_o)  = rsub(kant3_o,k1:k2,1)*factgas
            if (p_ant4_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_ant4_o)  = rsub(kant4_o,k1:k2,1)*factgas
            if (p_biog1_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_biog1_c)  = rsub(kbiog1_c,k1:k2,1)*factgas
            if (p_biog2_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_biog2_c)  = rsub(kbiog2_c,k1:k2,1)*factgas
            if (p_biog3_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_biog3_c)  = rsub(kbiog3_c,k1:k2,1)*factgas
            if (p_biog4_c .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_biog4_c)  = rsub(kbiog4_c,k1:k2,1)*factgas
            if (p_biog1_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_biog1_o)  = rsub(kbiog1_o,k1:k2,1)*factgas
            if (p_biog2_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_biog2_o)  = rsub(kbiog2_o,k1:k2,1)*factgas
            if (p_biog3_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_biog3_o)  = rsub(kbiog3_o,k1:k2,1)*factgas
            if (p_biog4_o .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_biog4_o)  = rsub(kbiog4_o,k1:k2,1)*factgas
            if (p_cvasoaX .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cvasoaX)  = rsub(kasoaX,k1:k2,1)*factgas
            if (p_cvasoa1 .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cvasoa1)  = rsub(kasoa1,k1:k2,1)*factgas
            if (p_cvasoa2 .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cvasoa2)  = rsub(kasoa2,k1:k2,1)*factgas
            if (p_cvasoa3 .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cvasoa3)  = rsub(kasoa3,k1:k2,1)*factgas
            if (p_cvasoa4 .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cvasoa4)  = rsub(kasoa4,k1:k2,1)*factgas
            if (p_cvbsoaX .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cvbsoaX)  = rsub(kbsoaX,k1:k2,1)*factgas
            if (p_cvbsoa1 .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cvbsoa1)  = rsub(kbsoa1,k1:k2,1)*factgas
            if (p_cvbsoa2 .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cvbsoa2)  = rsub(kbsoa2,k1:k2,1)*factgas
            if (p_cvbsoa3 .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cvbsoa3)  = rsub(kbsoa3,k1:k2,1)*factgas
            if (p_cvbsoa4 .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_cvbsoa4)  = rsub(kbsoa4,k1:k2,1)*factgas

	    
	    if (p_TERP_CN3_g .ge. p1st)   &
                chem(it,kt1:kt2,jt,p_TERP_CN3_g)  = rsub(kterp_g_Cn3,k1:k2,1)*factgas
	  end if

	    do iphase=1,nphase_aer
	    do itype=1,ntype_aer
	    do n = 1, nsize_aer(itype)
		chem(it,kt1:kt2,jt,lptr_so4_aer(n,itype,iphase)) =   &
		    rsub(lptr_so4_aer(n,itype,iphase),k1:k2,1)*factaerso4
		chem(it,kt1:kt2,jt,numptr_aer(n,itype,iphase)) =   &
		    rsub(numptr_aer(n,itype,iphase),k1:k2,1)*factaernum

		if (lptr_no3_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_no3_aer(n,itype,iphase)) =   &
		    rsub(lptr_no3_aer(n,itype,iphase),k1:k2,1)*factaerno3
		if (lptr_cl_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_cl_aer(n,itype,iphase)) =   &
		    rsub(lptr_cl_aer(n,itype,iphase),k1:k2,1)*factaercl
		if (lptr_msa_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_msa_aer(n,itype,iphase)) =   &
		    rsub(lptr_msa_aer(n,itype,iphase),k1:k2,1)*factaermsa
		if (lptr_co3_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_co3_aer(n,itype,iphase)) =   &
		    rsub(lptr_co3_aer(n,itype,iphase),k1:k2,1)*factaerco3
		if (lptr_nh4_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_nh4_aer(n,itype,iphase)) =   &
		    rsub(lptr_nh4_aer(n,itype,iphase),k1:k2,1)*factaernh4
		if (lptr_na_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_na_aer(n,itype,iphase)) =   &
		    rsub(lptr_na_aer(n,itype,iphase),k1:k2,1)*factaerna
		if (lptr_ca_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_ca_aer(n,itype,iphase)) =   &
		    rsub(lptr_ca_aer(n,itype,iphase),k1:k2,1)*factaerca
		if (lptr_oin_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_oin_aer(n,itype,iphase)) =   &
		    rsub(lptr_oin_aer(n,itype,iphase),k1:k2,1)*factaeroin
		if (lptr_oc_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_oc_aer(n,itype,iphase)) =   &
		    rsub(lptr_oc_aer(n,itype,iphase),k1:k2,1)*factaeroc
		if (lptr_bc_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_bc_aer(n,itype,iphase)) =   &
		    rsub(lptr_bc_aer(n,itype,iphase),k1:k2,1)*factaerbc
		if (hyswptr_aer(n,itype) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,hyswptr_aer(n,itype)) =   &
		    rsub(hyswptr_aer(n,itype),k1:k2,1)*factaerhysw
		if (waterptr_aer(n,itype) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,waterptr_aer(n,itype)) =   &
		    rsub(waterptr_aer(n,itype),k1:k2,1)*factaerwater
                if (lptr_tanv_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_tanv_aer(n,itype,iphase)) =   &
                    rsub(lptr_tanv_aer(n,itype,iphase),k1:k2,1)*factaertanv
                if (lptr_sopnv_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_sopnv_aer(n,itype,iphase)) =   &
                    rsub(lptr_sopnv_aer(n,itype,iphase),k1:k2,1)*factaersopnv
                if (lptr_ternv_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_ternv_aer(n,itype,iphase)) =   &
                    rsub(lptr_ternv_aer(n,itype,iphase),k1:k2,1)*factaerternv
                if (lptr_seqnv_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_seqnv_aer(n,itype,iphase)) =   &
                    rsub(lptr_seqnv_aer(n,itype,iphase),k1:k2,1)*factaerseqnv
                if (lptr_sianv_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_sianv_aer(n,itype,iphase)) =   &
                    rsub(lptr_sianv_aer(n,itype,iphase),k1:k2,1)*factaersianv
                if (lptr_pcg1_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg1_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg1_b_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg1_b_c
                if (lptr_pcg2_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg2_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg2_b_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg2_b_c
                if (lptr_pcg3_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg3_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg3_b_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg3_b_c
                if (lptr_pcg4_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg4_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg4_b_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg4_b_c
                if (lptr_pcg5_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg5_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg5_b_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg5_b_c
                if (lptr_pcg6_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg6_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg6_b_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg6_b_c
                if (lptr_pcg7_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg7_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg7_b_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg7_b_c
                if (lptr_pcg8_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg8_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg8_b_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg8_b_c
                if (lptr_pcg9_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg9_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg9_b_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg9_b_c
                if (lptr_pcg1_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg1_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg1_b_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg1_b_o
                if (lptr_pcg2_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg2_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg2_b_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg2_b_o
                if (lptr_pcg3_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg3_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg3_b_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg3_b_o
                if (lptr_pcg4_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg4_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg4_b_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg4_b_o
                if (lptr_pcg5_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg5_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg5_b_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg5_b_o
                if (lptr_pcg6_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg6_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg6_b_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg6_b_o
                if (lptr_pcg7_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg7_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg7_b_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg7_b_o
                if (lptr_pcg8_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg8_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg8_b_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg8_b_o
                if (lptr_pcg9_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg9_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg9_b_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg9_b_o
                if (lptr_opcg1_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg1_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg1_b_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg1_b_c
                if (lptr_opcg2_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg2_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg2_b_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg2_b_c
                if (lptr_opcg3_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg3_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg3_b_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg3_b_c
                if (lptr_opcg4_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg4_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg4_b_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg4_b_c
                if (lptr_opcg5_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg5_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg5_b_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg5_b_c
                if (lptr_opcg6_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg6_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg6_b_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg6_b_c
                if (lptr_opcg7_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg7_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg7_b_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg7_b_c
                if (lptr_opcg8_b_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg8_b_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg8_b_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg8_b_c
                if (lptr_opcg1_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg1_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg1_b_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg1_b_o
                if (lptr_opcg2_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg2_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg2_b_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg2_b_o
                if (lptr_opcg3_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg3_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg3_b_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg3_b_o
                if (lptr_opcg4_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg4_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg4_b_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg4_b_o
                if (lptr_opcg5_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg5_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg5_b_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg5_b_o
                if (lptr_opcg6_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg6_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg6_b_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg6_b_o
                if (lptr_opcg7_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg7_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg7_b_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg7_b_o
                if (lptr_opcg8_b_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg8_b_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg8_b_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg8_b_o
                if (lptr_pcg1_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg1_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg1_f_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg1_f_c
                if (lptr_pcg2_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg2_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg2_f_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg2_f_c
                if (lptr_pcg3_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg3_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg3_f_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg3_f_c
                if (lptr_pcg4_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg4_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg4_f_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg4_f_c
                if (lptr_pcg5_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg5_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg5_f_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg5_f_c
                if (lptr_pcg6_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg6_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg6_f_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg6_f_c
                if (lptr_pcg7_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg7_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg7_f_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg7_f_c
                if (lptr_pcg8_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg8_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg8_f_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg8_f_c
                if (lptr_pcg9_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg9_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg9_f_c_aer(n,itype,iphase),k1:k2,1)*factaerpcg9_f_c
                if (lptr_pcg1_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg1_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg1_f_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg1_f_o
                if (lptr_pcg2_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg2_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg2_f_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg2_f_o
                if (lptr_pcg3_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg3_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg3_f_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg3_f_o
                if (lptr_pcg4_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg4_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg4_f_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg4_f_o
                if (lptr_pcg5_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg5_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg5_f_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg5_f_o
                if (lptr_pcg6_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg6_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg6_f_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg6_f_o
                if (lptr_pcg7_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg7_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg7_f_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg7_f_o
                if (lptr_pcg8_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg8_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg8_f_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg8_f_o
                if (lptr_pcg9_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_pcg9_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_pcg9_f_o_aer(n,itype,iphase),k1:k2,1)*factaerpcg9_f_o
                if (lptr_opcg1_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg1_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg1_f_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg1_f_c
                if (lptr_opcg2_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg2_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg2_f_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg2_f_c
                if (lptr_opcg3_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg3_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg3_f_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg3_f_c
                if (lptr_opcg4_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg4_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg4_f_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg4_f_c
                if (lptr_opcg5_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg5_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg5_f_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg5_f_c
                if (lptr_opcg6_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg6_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg6_f_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg6_f_c
                if (lptr_opcg7_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg7_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg7_f_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg7_f_c
                if (lptr_opcg8_f_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg8_f_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg8_f_c_aer(n,itype,iphase),k1:k2,1)*factaeropcg8_f_c
                if (lptr_opcg1_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg1_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg1_f_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg1_f_o
                if (lptr_opcg2_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg2_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg2_f_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg2_f_o
                if (lptr_opcg3_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg3_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg3_f_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg3_f_o
                if (lptr_opcg4_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg4_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg4_f_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg4_f_o
                if (lptr_opcg5_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg5_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg5_f_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg5_f_o
                if (lptr_opcg6_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg6_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg6_f_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg6_f_o
                if (lptr_opcg7_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg7_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg7_f_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg7_f_o
                if (lptr_opcg8_f_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_opcg8_f_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_opcg8_f_o_aer(n,itype,iphase),k1:k2,1)*factaeropcg8_f_o
               if (lptr_smpa_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_smpa_aer(n,itype,iphase)) =   &
                    rsub(lptr_smpa_aer(n,itype,iphase),k1:k2,1)*factaersmpa
                if (lptr_smpbb_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_smpbb_aer(n,itype,iphase)) =   &
                    rsub(lptr_smpbb_aer(n,itype,iphase),k1:k2,1)*factaersmpbb
               if (lptr_iepox_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_iepox_aer(n,itype,iphase)) =   &
                    rsub(lptr_iepox_aer(n,itype,iphase),k1:k2,1)*factaeriepox
               if (lptr_iepoxos_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_iepoxos_aer(n,itype,iphase)) =   &
                    rsub(lptr_iepoxos_aer(n,itype,iphase),k1:k2,1)*factaeriepoxos
               if (lptr_tetrol_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_tetrol_aer(n,itype,iphase)) =   &
                    rsub(lptr_tetrol_aer(n,itype,iphase),k1:k2,1)*factaertetrol
               if (lptr_cldiepoxos_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_cldiepoxos_aer(n,itype,iphase)) =   &
                    rsub(lptr_cldiepoxos_aer(n,itype,iphase),k1:k2,1)*factaercldiepoxos
               if (lptr_cldtetrol_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_cldtetrol_aer(n,itype,iphase)) =   &
                    rsub(lptr_cldtetrol_aer(n,itype,iphase),k1:k2,1)*factaercldtetrol
               if (lptr_gly_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_gly_aer(n,itype,iphase)) =   &
                    rsub(lptr_gly_aer(n,itype,iphase),k1:k2,1)*factaergly
                if (lptr_glysoa_r1_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_glysoa_r1_aer(n,itype,iphase)) =   &
                    rsub(lptr_glysoa_r1_aer(n,itype,iphase),k1:k2,1)*factaerglyr1
                if (lptr_glysoa_r2_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_glysoa_r2_aer(n,itype,iphase)) =   &
                    rsub(lptr_glysoa_r2_aer(n,itype,iphase),k1:k2,1)*factaerglyr2
                if (lptr_glysoa_sfc_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_glysoa_sfc_aer(n,itype,iphase)) =   &
                    rsub(lptr_glysoa_sfc_aer(n,itype,iphase),k1:k2,1)*factaerglysfc
                if (lptr_glysoa_oh_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_glysoa_oh_aer(n,itype,iphase)) =   &
                    rsub(lptr_glysoa_oh_aer(n,itype,iphase),k1:k2,1)*factaerglyoh
                if (lptr_glysoa_nh4_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_glysoa_nh4_aer(n,itype,iphase)) =   &
                    rsub(lptr_glysoa_nh4_aer(n,itype,iphase),k1:k2,1)*factaerglynh4
                if (lptr_ant1_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_ant1_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_ant1_c_aer(n,itype,iphase),k1:k2,1)*factaerant1_c
                if (lptr_ant2_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_ant2_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_ant2_c_aer(n,itype,iphase),k1:k2,1)*factaerant2_c
                if (lptr_ant3_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_ant3_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_ant3_c_aer(n,itype,iphase),k1:k2,1)*factaerant3_c
                if (lptr_ant4_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_ant4_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_ant4_c_aer(n,itype,iphase),k1:k2,1)*factaerant4_c
                if (lptr_ant1_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_ant1_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_ant1_o_aer(n,itype,iphase),k1:k2,1)*factaerant1_o
                if (lptr_ant2_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_ant2_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_ant2_o_aer(n,itype,iphase),k1:k2,1)*factaerant2_o
                if (lptr_ant3_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_ant3_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_ant3_o_aer(n,itype,iphase),k1:k2,1)*factaerant3_o
                if (lptr_ant4_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_ant4_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_ant4_o_aer(n,itype,iphase),k1:k2,1)*factaerant4_o
                if (lptr_biog1_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_biog1_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_biog1_c_aer(n,itype,iphase),k1:k2,1)*factaerbiog1_c
                if (lptr_biog2_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_biog2_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_biog2_c_aer(n,itype,iphase),k1:k2,1)*factaerbiog2_c
                if (lptr_biog3_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_biog3_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_biog3_c_aer(n,itype,iphase),k1:k2,1)*factaerbiog3_c
                if (lptr_biog4_c_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_biog4_c_aer(n,itype,iphase)) =   &
                    rsub(lptr_biog4_c_aer(n,itype,iphase),k1:k2,1)*factaerbiog4_c
                if (lptr_biog1_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_biog1_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_biog1_o_aer(n,itype,iphase),k1:k2,1)*factaerbiog1_o
                if (lptr_biog2_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_biog2_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_biog2_o_aer(n,itype,iphase),k1:k2,1)*factaerbiog2_o
                if (lptr_biog3_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_biog3_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_biog3_o_aer(n,itype,iphase),k1:k2,1)*factaerbiog3_o
                if (lptr_biog4_o_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_biog4_o_aer(n,itype,iphase)) =   &
                    rsub(lptr_biog4_o_aer(n,itype,iphase),k1:k2,1)*factaerbiog4_o
                if (lptr_asoaX_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_asoaX_aer(n,itype,iphase)) =   &
                    rsub(lptr_asoaX_aer(n,itype,iphase),k1:k2,1)*factaerasoaX
                if (lptr_asoa1_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_asoa1_aer(n,itype,iphase)) =   &
                    rsub(lptr_asoa1_aer(n,itype,iphase),k1:k2,1)*factaerasoa1
                if (lptr_asoa2_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_asoa2_aer(n,itype,iphase)) =   &
                    rsub(lptr_asoa2_aer(n,itype,iphase),k1:k2,1)*factaerasoa2
                if (lptr_asoa3_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_asoa3_aer(n,itype,iphase)) =   &
                    rsub(lptr_asoa3_aer(n,itype,iphase),k1:k2,1)*factaerasoa3
                if (lptr_asoa4_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_asoa4_aer(n,itype,iphase)) =   &
                    rsub(lptr_asoa4_aer(n,itype,iphase),k1:k2,1)*factaerasoa4
                if (lptr_bsoaX_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_bsoaX_aer(n,itype,iphase)) =   &
                    rsub(lptr_bsoaX_aer(n,itype,iphase),k1:k2,1)*factaerbsoaX
                if (lptr_bsoa1_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_bsoa1_aer(n,itype,iphase)) =   &
                    rsub(lptr_bsoa1_aer(n,itype,iphase),k1:k2,1)*factaerbsoa1
                if (lptr_bsoa2_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_bsoa2_aer(n,itype,iphase)) =   &
                    rsub(lptr_bsoa2_aer(n,itype,iphase),k1:k2,1)*factaerbsoa2
                if (lptr_bsoa3_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_bsoa3_aer(n,itype,iphase)) =   &
                    rsub(lptr_bsoa3_aer(n,itype,iphase),k1:k2,1)*factaerbsoa3
                if (lptr_bsoa4_aer(n,itype,iphase) .ge. p1st)   &
                    chem(it,kt1:kt2,jt,lptr_bsoa4_aer(n,itype,iphase)) =   &
                    rsub(lptr_bsoa4_aer(n,itype,iphase),k1:k2,1)*factaerbsoa4

                      
		
		if (lptr_terp_Cn3_aer(n,itype,iphase) .ge. p1st)   &
		    chem(it,kt1:kt2,jt,lptr_terp_Cn3_aer(n,itype,iphase)) =   &
                    rsub(lptr_terp_Cn3_aer(n,itype,iphase),k1:k2,1)*factaerterp_Cn3

	    end do 
	    end do 
	    end do 
	end if


	return

	end subroutine mapaer_tofrom_host




	subroutine init_data_mosaic_asect( id, config_flags, vbs_nbin, is_aerosol )



	use module_data_mosaic_asect
	use module_data_mosaic_other, only:  lunerr, lunout,   &
		aboxtest_testmode, aboxtest_units_convert,   &
		aboxtest_rh_method, aboxtest_map_method,   &
		aboxtest_gases_fixed, aboxtest_min_temp,   &
		aboxtest_min_relhum, aboxtest_max_relhum
	use module_data_mosaic_therm, only:  nbin_a, nbin_a_maxd, n2o5_flag
	use module_mosaic2_driver, only:  init_data_mosaic2_asect
	use module_mosaic_csuesat, only:  init_csuesat
	use module_mosaic_movesect, only:  move_sections, test_move_sections
	use module_peg_util, only:  peg_error_fatal


	use module_configure, only:   &
		grid_config_rec_type, &
		p_so4_a01, p_so4_a02, p_so4_a03, p_so4_a04,   &
		p_so4_a05, p_so4_a06, p_so4_a07, p_so4_a08,   &
                p_so4_a09, p_so4_a10, p_so4_a11, p_so4_a12,   &  
                p_so4_a13, p_so4_a14, p_so4_a15, p_so4_a16,   &
                p_so4_a17, p_so4_a18, p_so4_a19, p_so4_a20

        use module_configure, only:   &
	   p_so4_cw01, p_no3_cw01, p_cl_cw01, p_nh4_cw01, p_na_cw01,   &
	   p_so4_cw02, p_no3_cw02, p_cl_cw02, p_nh4_cw02, p_na_cw02,   &
	   p_so4_cw03, p_no3_cw03, p_cl_cw03, p_nh4_cw03, p_na_cw03,   &
	   p_so4_cw04, p_no3_cw04, p_cl_cw04, p_nh4_cw04, p_na_cw04,   &
	   p_so4_cw05, p_no3_cw05, p_cl_cw05, p_nh4_cw05, p_na_cw05,   &
	   p_so4_cw06, p_no3_cw06, p_cl_cw06, p_nh4_cw06, p_na_cw06,   &
	   p_so4_cw07, p_no3_cw07, p_cl_cw07, p_nh4_cw07, p_na_cw07,   &
	   p_so4_cw08, p_no3_cw08, p_cl_cw08, p_nh4_cw08, p_na_cw08,   &
           p_so4_cw09, p_no3_cw09, p_cl_cw09, p_nh4_cw09, p_na_cw09,   & 
           p_so4_cw10, p_no3_cw10, p_cl_cw10, p_nh4_cw10, p_na_cw10,   &
           p_so4_cw11, p_no3_cw11, p_cl_cw11, p_nh4_cw11, p_na_cw11,   &
           p_so4_cw12, p_no3_cw12, p_cl_cw12, p_nh4_cw12, p_na_cw12,   &
           p_so4_cw13, p_no3_cw13, p_cl_cw13, p_nh4_cw13, p_na_cw13,   &
           p_so4_cw14, p_no3_cw14, p_cl_cw14, p_nh4_cw14, p_na_cw14,   &
           p_so4_cw15, p_no3_cw15, p_cl_cw15, p_nh4_cw15, p_na_cw15,   &
           p_so4_cw16, p_no3_cw16, p_cl_cw16, p_nh4_cw16, p_na_cw16,   &
           p_so4_cw17, p_no3_cw17, p_cl_cw17, p_nh4_cw17, p_na_cw17,   &
           p_so4_cw18, p_no3_cw18, p_cl_cw18, p_nh4_cw18, p_na_cw18,   &
           p_so4_cw19, p_no3_cw19, p_cl_cw19, p_nh4_cw19, p_na_cw19,   &
           p_so4_cw20, p_no3_cw20, p_cl_cw20, p_nh4_cw20, p_na_cw20,   &
	   p_oin_cw01, p_oc_cw01,  p_bc_cw01, p_num_cw01,              &
	   p_oin_cw02, p_oc_cw02,  p_bc_cw02, p_num_cw02,              &
	   p_oin_cw03, p_oc_cw03,  p_bc_cw03, p_num_cw03,              &
	   p_oin_cw04, p_oc_cw04,  p_bc_cw04, p_num_cw04,              &
	   p_oin_cw05, p_oc_cw05,  p_bc_cw05, p_num_cw05,              &
	   p_oin_cw06, p_oc_cw06,  p_bc_cw06, p_num_cw06,              &
	   p_oin_cw07, p_oc_cw07,  p_bc_cw07, p_num_cw07,              &
	   p_oin_cw08, p_oc_cw08,  p_bc_cw08, p_num_cw08,              &
           p_oin_cw09, p_oc_cw09,  p_bc_cw09, p_num_cw09,              &
           p_oin_cw10, p_oc_cw10,  p_bc_cw10, p_num_cw10,              &
           p_oin_cw11, p_oc_cw11,  p_bc_cw11, p_num_cw11,              &
           p_oin_cw12, p_oc_cw12,  p_bc_cw12, p_num_cw12,              &
           p_oin_cw13, p_oc_cw13,  p_bc_cw13, p_num_cw13,              &
           p_oin_cw14, p_oc_cw14,  p_bc_cw14, p_num_cw14,              &
           p_oin_cw15, p_oc_cw15,  p_bc_cw15, p_num_cw15,              &
           p_oin_cw16, p_oc_cw16,  p_bc_cw16, p_num_cw16,              &
           p_oin_cw17, p_oc_cw17,  p_bc_cw17, p_num_cw17,              &
           p_oin_cw18, p_oc_cw18,  p_bc_cw18, p_num_cw18,              &
           p_oin_cw19, p_oc_cw19,  p_bc_cw19, p_num_cw19,              &
           p_oin_cw20, p_oc_cw20,  p_bc_cw20, p_num_cw20,              &
     p_glysoa_r1_cw01, p_glysoa_r2_cw01, p_glysoa_sfc_cw01, p_glysoa_nh4_cw01, p_glysoa_oh_cw01,    &
     p_glysoa_r1_cw02, p_glysoa_r2_cw02, p_glysoa_sfc_cw02, p_glysoa_nh4_cw02, p_glysoa_oh_cw02,    &
     p_glysoa_r1_cw03, p_glysoa_r2_cw03, p_glysoa_sfc_cw03, p_glysoa_nh4_cw03, p_glysoa_oh_cw03,    &
     p_glysoa_r1_cw04, p_glysoa_r2_cw04, p_glysoa_sfc_cw04, p_glysoa_nh4_cw04, p_glysoa_oh_cw04,    &
     p_asoaX_cw01, p_asoa1_cw01, p_asoa2_cw01, p_asoa3_cw01, p_asoa4_cw01,     &
     p_bsoaX_cw01, p_bsoa1_cw01, p_bsoa2_cw01, p_bsoa3_cw01, p_bsoa4_cw01,     &
     p_asoaX_cw02, p_asoa1_cw02, p_asoa2_cw02, p_asoa3_cw02, p_asoa4_cw02,     &
     p_bsoaX_cw02, p_bsoa1_cw02, p_bsoa2_cw02, p_bsoa3_cw02, p_bsoa4_cw02,     &
     p_asoaX_cw03, p_asoa1_cw03, p_asoa2_cw03, p_asoa3_cw03, p_asoa4_cw03,     &
     p_bsoaX_cw03, p_bsoa1_cw03, p_bsoa2_cw03, p_bsoa3_cw03, p_bsoa4_cw03,     &
     p_asoaX_cw04, p_asoa1_cw04, p_asoa2_cw04, p_asoa3_cw04, p_asoa4_cw04,     &
     p_bsoaX_cw04, p_bsoa1_cw04, p_bsoa2_cw04, p_bsoa3_cw04, p_bsoa4_cw04,              &
       p_pcg1_b_c_cw01,p_pcg1_b_o_cw01,p_opcg1_b_c_cw01,p_opcg1_b_o_cw01, &
           p_pcg1_f_c_cw01,p_pcg1_f_o_cw01,p_opcg1_f_c_cw01,p_opcg1_f_o_cw01, &
           p_opcg2_f_o_cw01,p_opcg3_f_o_cw01,p_opcg4_f_o_cw01,p_opcg5_f_o_cw01,&
           p_opcg6_f_o_cw01, &
           p_opcg1_b_o_cw01,p_opcg2_b_o_cw01, p_opcg3_b_o_cw01, &
           p_opcg4_b_o_cw01,p_opcg5_b_o_cw01, p_opcg6_b_o_cw01, &
           p_opcg7_b_o_cw01,p_opcg8_b_o_cw01, &
           p_ant1_c_cw01,p_ant2_c_cw01,p_ant3_c_cw01,p_ant4_c_cw01,    &
           p_ant2_o_cw01,p_ant3_o_cw01,   &
           p_biog1_c_cw01, p_biog2_c_cw01, p_biog3_c_cw01, p_biog4_c_cw01,    &
           p_biog1_o_cw01, p_biog2_o_cw01, p_biog3_o_cw01, p_biog4_o_cw01,    &
           p_tanv_cw01,p_sopnv_cw01,p_seqnv_cw01,p_ternv_cw01,p_sianv_cw01, &
           p_iepox_cw01,     &
           p_iepoxos_cw01,     &
           p_tetrol_cw01,     &
           p_cldiepoxos_cw01,     &
           p_cldtetrol_cw01,     &
           p_gly_cw01,    &
           p_pcg1_b_c_cw02,p_pcg1_b_o_cw02,p_opcg1_b_c_cw02,p_opcg1_b_o_cw02, &
           p_pcg1_f_c_cw02,p_pcg1_f_o_cw02,p_opcg1_f_c_cw02,p_opcg1_f_o_cw02, &
           p_opcg2_f_o_cw02,p_opcg3_f_o_cw02,p_opcg4_f_o_cw02,p_opcg5_f_o_cw02,&
           p_opcg6_f_o_cw02, &
           p_opcg1_b_o_cw02,p_opcg2_b_o_cw02, p_opcg3_b_o_cw02, &
           p_opcg4_b_o_cw02,p_opcg5_b_o_cw02, p_opcg6_b_o_cw02, &
           p_opcg7_b_o_cw02,p_opcg8_b_o_cw02, &
           p_ant1_c_cw02,p_ant2_c_cw02,p_ant3_c_cw02,p_ant4_c_cw02,    &
           p_ant2_o_cw02,p_ant3_o_cw02, &
           p_biog1_c_cw02, p_biog2_c_cw02, p_biog3_c_cw02, p_biog4_c_cw02,    &
           p_biog1_o_cw02, p_biog2_o_cw02, p_biog3_o_cw02, p_biog4_o_cw02,    &
           p_iepox_cw02,     &
           p_iepoxos_cw02,     &
           p_tetrol_cw02,    &
           p_cldiepoxos_cw02,     &
           p_cldtetrol_cw02,    &
           p_gly_cw02,    &
           p_tanv_cw02,p_sopnv_cw02,p_seqnv_cw02,p_ternv_cw02,p_sianv_cw02, &
           p_pcg1_b_c_cw03,p_pcg1_b_o_cw03,p_opcg1_b_c_cw03,p_opcg1_b_o_cw03, &
           p_pcg1_f_c_cw03,p_pcg1_f_o_cw03,p_opcg1_f_c_cw03,p_opcg1_f_o_cw03, &
           p_opcg2_f_o_cw03,p_opcg3_f_o_cw03,p_opcg4_f_o_cw03,p_opcg5_f_o_cw03,&
           p_opcg6_f_o_cw03, &
           p_opcg1_b_o_cw03,p_opcg2_b_o_cw03, p_opcg3_b_o_cw03, &
           p_opcg4_b_o_cw03,p_opcg5_b_o_cw03, p_opcg6_b_o_cw03, &
           p_opcg7_b_o_cw03,p_opcg8_b_o_cw03, &
           p_ant1_c_cw03,p_ant2_c_cw03,p_ant3_c_cw03,p_ant4_c_cw03,    &
           p_ant2_o_cw03,p_ant3_o_cw03, &
           p_biog1_c_cw03, p_biog2_c_cw03, p_biog3_c_cw03, p_biog4_c_cw03,    &
           p_biog1_o_cw03, p_biog2_o_cw03, p_biog3_o_cw03, p_biog4_o_cw03,    &
           p_iepox_cw03,     &
           p_iepoxos_cw03,     &
           p_tetrol_cw03,     &
           p_cldiepoxos_cw03,     &
           p_cldtetrol_cw03,     &
           p_gly_cw03,     &
           p_tanv_cw03,p_sopnv_cw03,p_seqnv_cw03,p_ternv_cw03,p_sianv_cw03, &
           p_pcg1_b_c_cw04,p_pcg1_b_o_cw04,p_opcg1_b_c_cw04,p_opcg1_b_o_cw04, &
           p_pcg1_f_c_cw04,p_pcg1_f_o_cw04,p_opcg1_f_c_cw04,p_opcg1_f_o_cw04, &
           p_opcg2_f_o_cw04,p_opcg3_f_o_cw04,p_opcg4_f_o_cw04,p_opcg5_f_o_cw04,&
           p_opcg6_f_o_cw04, &
           p_opcg1_b_o_cw04,p_opcg2_b_o_cw04, p_opcg3_b_o_cw04, &
           p_opcg4_b_o_cw04,p_opcg5_b_o_cw04, p_opcg6_b_o_cw04, &
           p_opcg7_b_o_cw04,p_opcg8_b_o_cw04, &
           p_ant1_c_cw04,p_ant2_c_cw04,p_ant3_c_cw04,p_ant4_c_cw04,    &
           p_ant2_o_cw04,p_ant3_o_cw04, &
           p_biog1_c_cw04, p_biog2_c_cw04, p_biog3_c_cw04, p_biog4_c_cw04,    &
           p_biog1_o_cw04, p_biog2_o_cw04, p_biog3_o_cw04, p_biog4_o_cw04,    &
           p_iepox_cw04,     &
           p_iepoxos_cw04,     &
           p_tetrol_cw04,     &
           p_cldiepoxos_cw04,     &
           p_cldtetrol_cw04,     &
           p_gly_cw04,     &
           p_tanv_cw04,p_sopnv_cw04,p_seqnv_cw04,p_ternv_cw04,p_sianv_cw04, &
           p_pcg1_b_c_cw05,p_pcg1_b_o_cw05,p_opcg1_b_c_cw05,p_opcg1_b_o_cw05, &
           p_pcg1_f_c_cw05,p_pcg1_f_o_cw05,p_opcg1_f_c_cw05,p_opcg1_f_o_cw05, &
           p_ant1_c_cw05,p_ant2_c_cw05,p_ant3_c_cw05,p_ant4_c_cw05,    &
           p_biog1_c_cw05, p_biog2_c_cw05, p_biog3_c_cw05, p_biog4_c_cw05,    &
           p_biog1_o_cw05, p_biog2_o_cw05, p_biog3_o_cw05, p_biog4_o_cw05,    &
           p_iepox_cw05,     &
           p_iepoxos_cw05,     &
           p_tetrol_cw05,     &
           p_gly_cw05,     &
           p_tanv_cw05,p_sopnv_cw05,p_seqnv_cw05,p_ternv_cw05,p_sianv_cw05, &
           p_pcg1_b_c_cw06,p_pcg1_b_o_cw06,p_opcg1_b_c_cw06,p_opcg1_b_o_cw06, &
           p_pcg1_f_c_cw06,p_pcg1_f_o_cw06,p_opcg1_f_c_cw06,p_opcg1_f_o_cw06, &
           p_ant1_c_cw06,p_ant2_c_cw06,p_ant3_c_cw06,p_ant4_c_cw06,    &
           p_biog1_c_cw06, p_biog2_c_cw06, p_biog3_c_cw06, p_biog4_c_cw06,    &
           p_biog1_o_cw06, p_biog2_o_cw06, p_biog3_o_cw06, p_biog4_o_cw06,    &
           p_iepox_cw06,     &
           p_iepoxos_cw06,     &
           p_tetrol_cw06,     &
           p_gly_cw06,     &
           p_tanv_cw06,p_sopnv_cw06,p_seqnv_cw06,p_ternv_cw06,p_sianv_cw06, &
           p_pcg1_b_c_cw07,p_pcg1_b_o_cw07,p_opcg1_b_c_cw07,p_opcg1_b_o_cw07, &
           p_pcg1_f_c_cw07,p_pcg1_f_o_cw07,p_opcg1_f_c_cw07,p_opcg1_f_o_cw07, &
           p_ant1_c_cw07,p_ant2_c_cw07,p_ant3_c_cw07,p_ant4_c_cw07,    &
           p_biog1_c_cw07, p_biog2_c_cw07, p_biog3_c_cw07, p_biog4_c_cw07,    &
           p_biog1_o_cw07, p_biog2_o_cw07, p_biog3_o_cw07, p_biog4_o_cw07,    &
           p_iepox_cw07,     &
           p_iepoxos_cw07,     &
           p_tetrol_cw07,    &
           p_gly_cw07,    &
           p_tanv_cw07,p_sopnv_cw07,p_seqnv_cw07,p_ternv_cw07,p_sianv_cw07, &
           p_pcg1_b_c_cw08,p_pcg1_b_o_cw08,p_opcg1_b_c_cw08,p_opcg1_b_o_cw08, &
           p_pcg1_f_c_cw08,p_pcg1_f_o_cw08,p_opcg1_f_c_cw08,p_opcg1_f_o_cw08, &
           p_ant1_c_cw08,p_ant2_c_cw08,p_ant3_c_cw08,p_ant4_c_cw08,    &
           p_biog1_c_cw08, p_biog2_c_cw08, p_biog3_c_cw08, p_biog4_c_cw08,    &
           p_biog1_o_cw08, p_biog2_o_cw08, p_biog3_o_cw08, p_biog4_o_cw08,    &
           p_iepox_cw08,     &
           p_iepoxos_cw08,     &
           p_tetrol_cw08,    &
           p_gly_cw08,     &
           p_tanv_cw08,p_sopnv_cw08,p_seqnv_cw08,p_ternv_cw08,p_sianv_cw08,   &

           p_pcg1_b_c_cw09,p_pcg1_b_o_cw09,p_opcg1_b_c_cw09,p_opcg1_b_o_cw09, &
           p_pcg1_f_c_cw09,p_pcg1_f_o_cw09,p_opcg1_f_c_cw09,p_opcg1_f_o_cw09, &
           p_ant1_c_cw09,p_biog1_c_cw09,                                      &
           p_pcg1_b_c_cw10,p_pcg1_b_o_cw10,p_opcg1_b_c_cw10,p_opcg1_b_o_cw10, &
           p_pcg1_f_c_cw10,p_pcg1_f_o_cw10,p_opcg1_f_c_cw10,p_opcg1_f_o_cw10, &
           p_ant1_c_cw10,p_biog1_c_cw10,                                      &
           p_pcg1_b_c_cw11,p_pcg1_b_o_cw11,p_opcg1_b_c_cw11,p_opcg1_b_o_cw11, &
           p_pcg1_f_c_cw11,p_pcg1_f_o_cw11,p_opcg1_f_c_cw11,p_opcg1_f_o_cw11, &
           p_ant1_c_cw11,p_biog1_c_cw11,                                      &
           p_pcg1_b_c_cw12,p_pcg1_b_o_cw12,p_opcg1_b_c_cw12,p_opcg1_b_o_cw12, &
           p_pcg1_f_c_cw12,p_pcg1_f_o_cw12,p_opcg1_f_c_cw12,p_opcg1_f_o_cw12, &
           p_ant1_c_cw12,p_biog1_c_cw12,                                      &
           p_pcg1_b_c_cw13,p_pcg1_b_o_cw13,p_opcg1_b_c_cw13,p_opcg1_b_o_cw13, &
           p_pcg1_f_c_cw13,p_pcg1_f_o_cw13,p_opcg1_f_c_cw13,p_opcg1_f_o_cw13, &
           p_ant1_c_cw13,p_biog1_c_cw13,                                      &
           p_pcg1_b_c_cw14,p_pcg1_b_o_cw14,p_opcg1_b_c_cw14,p_opcg1_b_o_cw14, &
           p_pcg1_f_c_cw14,p_pcg1_f_o_cw14,p_opcg1_f_c_cw14,p_opcg1_f_o_cw14, &
           p_ant1_c_cw14,p_biog1_c_cw14,                                      &
           p_pcg1_b_c_cw15,p_pcg1_b_o_cw15,p_opcg1_b_c_cw15,p_opcg1_b_o_cw15, &
           p_pcg1_f_c_cw15,p_pcg1_f_o_cw15,p_opcg1_f_c_cw15,p_opcg1_f_o_cw15, &
           p_ant1_c_cw15,p_biog1_c_cw15,                                      &
           p_pcg1_b_c_cw16,p_pcg1_b_o_cw16,p_opcg1_b_c_cw16,p_opcg1_b_o_cw16, &
           p_pcg1_f_c_cw16,p_pcg1_f_o_cw16,p_opcg1_f_c_cw16,p_opcg1_f_o_cw16, &
           p_ant1_c_cw16,p_biog1_c_cw16,                                      &
           p_pcg1_b_c_cw17,p_pcg1_b_o_cw17,p_opcg1_b_c_cw17,p_opcg1_b_o_cw17, &
           p_pcg1_f_c_cw17,p_pcg1_f_o_cw17,p_opcg1_f_c_cw17,p_opcg1_f_o_cw17, &
           p_ant1_c_cw17,p_biog1_c_cw17,                                      &
           p_pcg1_b_c_cw18,p_pcg1_b_o_cw18,p_opcg1_b_c_cw18,p_opcg1_b_o_cw18, &
           p_pcg1_f_c_cw18,p_pcg1_f_o_cw18,p_opcg1_f_c_cw18,p_opcg1_f_o_cw18, &
           p_ant1_c_cw18,p_biog1_c_cw18,                                      &
           p_pcg1_b_c_cw19,p_pcg1_b_o_cw19,p_opcg1_b_c_cw19,p_opcg1_b_o_cw19, &
           p_pcg1_f_c_cw19,p_pcg1_f_o_cw19,p_opcg1_f_c_cw19,p_opcg1_f_o_cw19, &
           p_ant1_c_cw19,p_biog1_c_cw19,                                      &
           p_pcg1_b_c_cw20,p_pcg1_b_o_cw20,p_opcg1_b_c_cw20,p_opcg1_b_o_cw20, &
           p_pcg1_f_c_cw20,p_pcg1_f_o_cw20,p_opcg1_f_c_cw20,p_opcg1_f_o_cw20, &
           p_ant1_c_cw20,p_biog1_c_cw20, &
	   
	   p_terp_cw01_Cn3, p_terp_cw02_Cn3, p_terp_cw03_Cn3, p_terp_cw04_Cn3

	   

	use module_state_description, only:  param_first_scalar, num_chem

	implicit none


	logical, intent(out) :: is_aerosol(num_chem)

	integer, intent(in) :: id, &
	                       vbs_nbin(1)

	type(grid_config_rec_type), intent(in) :: config_flags


	integer idum, itype, jdum, l, ldum, n, nhi, nsize_aer_dum
	real dum
	real, parameter :: pi = 3.14159265




	msectional = 20
	maerocoag = -2
	maerchem = 1
	maeroptical = 1
	maerchem_boxtest_output = -1

        idum = max( config_flags%mosaic_aerchem_optaa, 0 )
        if (idum >= 10000) then
            jdum = mod(idum,100)/10  
            if (jdum == 0) msectional =  0  
            if (jdum == 1) msectional = 10  
            if (jdum == 2) msectional = 20  

            jdum = mod(idum,1000)/100  

            jdum = mod(idum,10000)/1000  
        end if 




	ntype_aer = 1




	nsize_aer(:) = 0
        itype=1
	if (p_so4_a01 .ge. param_first_scalar) nsize_aer(itype) = 1
	if (p_so4_a02 .ge. param_first_scalar) nsize_aer(itype) = 2
	if (p_so4_a03 .ge. param_first_scalar) nsize_aer(itype) = 3
	if (p_so4_a04 .ge. param_first_scalar) nsize_aer(itype) = 4
	if (p_so4_a05 .ge. param_first_scalar) nsize_aer(itype) = 5
	if (p_so4_a06 .ge. param_first_scalar) nsize_aer(itype) = 6
	if (p_so4_a07 .ge. param_first_scalar) nsize_aer(itype) = 7
	if (p_so4_a08 .ge. param_first_scalar) nsize_aer(itype) = 8
        if (p_so4_a09 .ge. param_first_scalar) nsize_aer(itype) = 9 
        if (p_so4_a10 .ge. param_first_scalar) nsize_aer(itype) = 10
        if (p_so4_a11 .ge. param_first_scalar) nsize_aer(itype) = 11
        if (p_so4_a12 .ge. param_first_scalar) nsize_aer(itype) = 12
        if (p_so4_a13 .ge. param_first_scalar) nsize_aer(itype) = 13
        if (p_so4_a14 .ge. param_first_scalar) nsize_aer(itype) = 14
        if (p_so4_a15 .ge. param_first_scalar) nsize_aer(itype) = 15
        if (p_so4_a16 .ge. param_first_scalar) nsize_aer(itype) = 16
        if (p_so4_a17 .ge. param_first_scalar) nsize_aer(itype) = 17
        if (p_so4_a18 .ge. param_first_scalar) nsize_aer(itype) = 18
        if (p_so4_a19 .ge. param_first_scalar) nsize_aer(itype) = 19
        if (p_so4_a20 .ge. param_first_scalar) nsize_aer(itype) = 20

	if (nsize_aer(itype) .le. 0) then
	    call peg_error_fatal( lunerr,   &
		'init_data_mosaic_asect - nsize_aer = 0' )
	else if (nsize_aer(itype) .gt. maxd_asize) then
	    call peg_error_fatal( lunerr,   &
		'init_data_mosaic_asect - nsize_aer > maxd_asize' )
	end if




	nbin_a = 0
	do itype = 1, ntype_aer
	    nbin_a = nbin_a + nsize_aer(itype)
	end do
	if (nbin_a .gt. nbin_a_maxd) then
	    call peg_error_fatal( lunerr,   &
		'init_data_mosaic_asect - nbin_a > nbin_a_maxd' )
	end if





	nphase_aer = 0
	maerosolincw = 0
	if (nsize_aer(1) .gt. 0) then
	    nphase_aer = 1
	    ai_phase = 1

	    if (p_so4_cw01 .ge. param_first_scalar) then
		nphase_aer = 2
		cw_phase = 2
		maerosolincw = 1
	    end if
	end if




n2o5_flag = config_flags%n2o5_hetchem





	
	ntot_mastercomp_aer = 123

	l = 1
	mastercompindx_so4_aer = l
	name_mastercomp_aer( l ) = 'sulfate'
	dens_mastercomp_aer( l ) =  dens_so4_aer
	mw_mastercomp_aer(   l ) =    mw_so4_aer
	hygro_mastercomp_aer(l ) = hygro_so4_aer

	l = 2
	mastercompindx_no3_aer = l
	name_mastercomp_aer( l ) = 'nitrate'
	dens_mastercomp_aer( l ) =  dens_no3_aer
	mw_mastercomp_aer(   l ) =    mw_no3_aer
	hygro_mastercomp_aer(l ) = hygro_no3_aer

	l = 3
	mastercompindx_cl_aer = l
	name_mastercomp_aer( l ) = 'chloride'
	dens_mastercomp_aer( l ) =  dens_cl_aer
	mw_mastercomp_aer(   l ) =    mw_cl_aer
	hygro_mastercomp_aer(l ) = hygro_cl_aer

	l = 4
	mastercompindx_co3_aer = l
	name_mastercomp_aer( l ) = 'carbonate'
	dens_mastercomp_aer( l ) =  dens_co3_aer
	mw_mastercomp_aer(   l ) =    mw_co3_aer
	hygro_mastercomp_aer(l ) = hygro_co3_aer

	l = 5
	mastercompindx_nh4_aer = l
	name_mastercomp_aer( l ) = 'ammonium'
	dens_mastercomp_aer( l ) =  dens_nh4_aer
	mw_mastercomp_aer(   l ) =    mw_nh4_aer
	hygro_mastercomp_aer(l ) = hygro_nh4_aer

	l = 6
	mastercompindx_na_aer = l
	name_mastercomp_aer( l ) = 'sodium'
	dens_mastercomp_aer( l ) =  dens_na_aer
	mw_mastercomp_aer(   l ) =    mw_na_aer
	hygro_mastercomp_aer(l ) = hygro_na_aer

	l = 7
	mastercompindx_ca_aer = l
	name_mastercomp_aer( l ) = 'calcium'
	dens_mastercomp_aer( l ) =  dens_ca_aer
	mw_mastercomp_aer(   l ) =    mw_ca_aer
	hygro_mastercomp_aer(l ) = hygro_ca_aer

	l = 8
	mastercompindx_oin_aer = l
	name_mastercomp_aer( l ) = 'otherinorg'
	dens_mastercomp_aer( l ) =  dens_oin_aer
	mw_mastercomp_aer(   l ) =    mw_oin_aer
	hygro_mastercomp_aer(l ) = hygro_oin_aer

	l = 9
	mastercompindx_oc_aer = l
	name_mastercomp_aer( l ) = 'organic-c'
	dens_mastercomp_aer( l ) =  dens_oc_aer
	mw_mastercomp_aer(   l ) =    mw_oc_aer
	hygro_mastercomp_aer(l ) = hygro_oc_aer

	l = 10
	mastercompindx_bc_aer = l
	name_mastercomp_aer( l ) = 'black-c'
	dens_mastercomp_aer( l ) =  dens_bc_aer
	mw_mastercomp_aer(   l ) =    mw_bc_aer
	hygro_mastercomp_aer(l ) = hygro_bc_aer

        l = 11
        mastercompindx_pcg1_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg1_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg1_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg1_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg1_b_c_aer

        l = 12
        mastercompindx_pcg2_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg2_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg2_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg2_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg2_b_c_aer

        l = 13
        mastercompindx_pcg3_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg3_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg3_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg3_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg3_b_c_aer

        l = 14
        mastercompindx_pcg4_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg4_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg4_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg4_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg4_b_c_aer

        l = 15
        mastercompindx_pcg5_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg5_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg5_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg5_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg5_b_c_aer

        l = 16
        mastercompindx_pcg6_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg6_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg6_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg6_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg6_b_c_aer

        l = 17
        mastercompindx_pcg7_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg7_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg7_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg7_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg7_b_c_aer

        l = 18
        mastercompindx_pcg8_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg8_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg8_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg8_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg8_b_c_aer

        l = 19
        mastercompindx_pcg9_b_c_aer = l
        name_mastercomp_aer( l ) = 'pcg9_b_c'
        dens_mastercomp_aer( l ) =  dens_pcg9_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg9_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg9_b_c_aer

        l = 20
        mastercompindx_pcg1_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg1_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg1_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg1_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg1_b_o_aer

        l = 21
        mastercompindx_pcg2_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg2_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg2_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg2_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg2_b_o_aer

        l = 22
        mastercompindx_pcg3_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg3_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg3_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg3_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg3_b_o_aer

        l = 23
        mastercompindx_pcg4_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg4_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg4_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg4_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg4_b_o_aer

        l = 24
        mastercompindx_pcg5_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg5_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg5_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg5_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg5_b_o_aer

        l = 25
        mastercompindx_pcg6_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg6_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg6_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg6_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg6_b_o_aer

        l = 26
        mastercompindx_pcg7_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg7_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg7_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg7_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg7_b_o_aer

        l = 27
        mastercompindx_pcg8_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg8_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg8_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg8_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg8_b_o_aer

        l = 28
        mastercompindx_pcg9_b_o_aer = l
        name_mastercomp_aer( l ) = 'pcg9_b_o'
        dens_mastercomp_aer( l ) =  dens_pcg9_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg9_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg9_b_o_aer

        l = 29
        mastercompindx_opcg1_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg1_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg1_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg1_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg1_b_c_aer

        l = 30
        mastercompindx_opcg2_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg2_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg2_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg2_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg2_b_c_aer

        l = 31
        mastercompindx_opcg3_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg3_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg3_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg3_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg3_b_c_aer

        l = 32
        mastercompindx_opcg4_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg4_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg4_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg4_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg4_b_c_aer

        l = 33
        mastercompindx_opcg5_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg5_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg5_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg5_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg5_b_c_aer

        l = 34
        mastercompindx_opcg6_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg6_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg6_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg6_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg6_b_c_aer

        l = 35
        mastercompindx_opcg7_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg7_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg7_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg7_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg7_b_c_aer

        l = 36
        mastercompindx_opcg8_b_c_aer = l
        name_mastercomp_aer( l ) = 'opcg8_b_c'
        dens_mastercomp_aer( l ) =  dens_opcg8_b_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg8_b_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg8_b_c_aer

        l = 37
        mastercompindx_opcg1_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg1_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg1_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg1_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg1_b_o_aer

        l = 38
        mastercompindx_opcg2_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg2_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg2_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg2_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg2_b_o_aer

        l = 39
        mastercompindx_opcg3_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg3_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg3_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg3_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg3_b_o_aer

        l = 40
        mastercompindx_opcg4_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg4_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg4_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg4_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg4_b_o_aer

        l = 41
        mastercompindx_opcg5_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg5_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg5_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg5_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg5_b_o_aer

        l = 42
        mastercompindx_opcg6_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg6_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg6_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg6_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg6_b_o_aer

        l = 43
        mastercompindx_opcg7_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg7_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg7_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg7_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg7_b_o_aer

        l = 44
        mastercompindx_opcg8_b_o_aer = l
        name_mastercomp_aer( l ) = 'opcg8_b_o'
        dens_mastercomp_aer( l ) =  dens_opcg8_b_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg8_b_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg8_b_o_aer

        l = 45
        mastercompindx_pcg1_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg1_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg1_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg1_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg1_f_c_aer

        l = 46
        mastercompindx_pcg2_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg2_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg2_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg2_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg2_f_c_aer

        l = 47
        mastercompindx_pcg3_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg3_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg3_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg3_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg3_f_c_aer

        l = 48
        mastercompindx_pcg4_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg4_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg4_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg4_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg4_f_c_aer

        l = 49
        mastercompindx_pcg5_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg5_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg5_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg5_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg5_f_c_aer

        l = 50
        mastercompindx_pcg6_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg6_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg6_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg6_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg6_f_c_aer

        l = 51
        mastercompindx_pcg7_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg7_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg7_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg7_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg7_f_c_aer

        l = 52
        mastercompindx_pcg8_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg8_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg8_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg8_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg8_f_c_aer

        l = 53
        mastercompindx_pcg9_f_c_aer = l
        name_mastercomp_aer( l ) = 'pcg9_f_c'
        dens_mastercomp_aer( l ) =  dens_pcg9_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_pcg9_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_pcg9_f_c_aer

        l = 54
        mastercompindx_pcg1_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg1_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg1_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg1_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg1_f_o_aer

        l = 55
        mastercompindx_pcg2_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg2_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg2_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg2_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg2_f_o_aer

        l = 56
        mastercompindx_pcg3_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg3_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg3_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg3_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg3_f_o_aer

        l = 57
        mastercompindx_pcg4_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg4_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg4_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg4_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg4_f_o_aer

        l = 58
        mastercompindx_pcg5_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg5_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg5_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg5_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg5_f_o_aer

        l = 59
        mastercompindx_pcg6_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg6_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg6_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg6_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg6_f_o_aer

        l = 60
        mastercompindx_pcg7_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg7_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg7_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg7_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg7_f_o_aer

        l = 61
        mastercompindx_pcg8_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg8_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg8_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg8_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg8_f_o_aer

        l = 62
        mastercompindx_pcg9_f_o_aer = l
        name_mastercomp_aer( l ) = 'pcg9_f_o'
        dens_mastercomp_aer( l ) =  dens_pcg9_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_pcg9_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_pcg9_f_o_aer

        l = 63
        mastercompindx_opcg1_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg1_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg1_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg1_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg1_f_c_aer

        l = 64
        mastercompindx_opcg2_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg2_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg2_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg2_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg2_f_c_aer

        l = 65
        mastercompindx_opcg3_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg3_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg3_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg3_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg3_f_c_aer

        l = 66
        mastercompindx_opcg4_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg4_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg4_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg4_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg4_f_c_aer

        l = 67
        mastercompindx_opcg5_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg5_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg5_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg5_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg5_f_c_aer

        l = 68
        mastercompindx_opcg6_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg6_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg6_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg6_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg6_f_c_aer

        l = 69
        mastercompindx_opcg7_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg7_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg7_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg7_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg7_f_c_aer

        l = 70
        mastercompindx_opcg8_f_c_aer = l
        name_mastercomp_aer( l ) = 'opcg8_f_c'
        dens_mastercomp_aer( l ) =  dens_opcg8_f_c_aer
        mw_mastercomp_aer(   l ) =    mw_opcg8_f_c_aer
        hygro_mastercomp_aer(l ) = hygro_opcg8_f_c_aer

        l = 71
        mastercompindx_opcg1_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg1_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg1_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg1_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg1_f_o_aer

        l = 72
        mastercompindx_opcg2_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg2_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg2_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg2_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg2_f_o_aer

        l = 73
        mastercompindx_opcg3_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg3_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg3_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg3_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg3_f_o_aer

        l = 74
        mastercompindx_opcg4_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg4_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg4_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg4_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg4_f_o_aer

        l = 75
        mastercompindx_opcg5_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg5_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg5_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg5_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg5_f_o_aer

        l = 76
        mastercompindx_opcg6_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg6_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg6_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg6_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg6_f_o_aer

        l = 77
        mastercompindx_opcg7_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg7_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg7_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg7_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg7_f_o_aer

        l = 78
        mastercompindx_opcg8_f_o_aer = l
        name_mastercomp_aer( l ) = 'opcg8_f_o'
        dens_mastercomp_aer( l ) =  dens_opcg8_f_o_aer
        mw_mastercomp_aer(   l ) =    mw_opcg8_f_o_aer
        hygro_mastercomp_aer(l ) = hygro_opcg8_f_o_aer

        l = 79
        mastercompindx_ant1_c_aer = l
        name_mastercomp_aer( l ) = 'ant1_c'
        dens_mastercomp_aer( l ) =  dens_ant1_c_aer
        mw_mastercomp_aer(   l ) =    mw_ant1_c_aer
        hygro_mastercomp_aer(l ) = hygro_ant1_c_aer

        l = 80
        mastercompindx_ant2_c_aer = l
        name_mastercomp_aer( l ) = 'ant2_c'
        dens_mastercomp_aer( l ) =  dens_ant2_c_aer
        mw_mastercomp_aer(   l ) =    mw_ant2_c_aer
        hygro_mastercomp_aer(l ) = hygro_ant2_c_aer

        l = 81
        mastercompindx_ant3_c_aer = l
        name_mastercomp_aer( l ) = 'ant3_c'
        dens_mastercomp_aer( l ) =  dens_ant3_c_aer
        mw_mastercomp_aer(   l ) =    mw_ant3_c_aer
        hygro_mastercomp_aer(l ) = hygro_ant3_c_aer

        l = 82
        mastercompindx_ant4_c_aer = l
        name_mastercomp_aer( l ) = 'ant4_c'
        dens_mastercomp_aer( l ) =  dens_ant4_c_aer
        mw_mastercomp_aer(   l ) =    mw_ant4_c_aer
        hygro_mastercomp_aer(l ) = hygro_ant4_c_aer

        l = 83
        mastercompindx_ant1_o_aer = l
        name_mastercomp_aer( l ) = 'ant1_o'
        dens_mastercomp_aer( l ) =  dens_ant1_o_aer
        mw_mastercomp_aer(   l ) =    mw_ant1_o_aer
        hygro_mastercomp_aer(l ) = hygro_ant1_o_aer

        l = 84
        mastercompindx_ant2_o_aer = l
        name_mastercomp_aer( l ) = 'ant2_o'
        dens_mastercomp_aer( l ) =  dens_ant2_o_aer
        mw_mastercomp_aer(   l ) =    mw_ant2_o_aer
        hygro_mastercomp_aer(l ) = hygro_ant2_o_aer

        l = 85
        mastercompindx_ant3_o_aer = l
        name_mastercomp_aer( l ) = 'ant3_o'
        dens_mastercomp_aer( l ) =  dens_ant3_o_aer
        mw_mastercomp_aer(   l ) =    mw_ant3_o_aer
        hygro_mastercomp_aer(l ) = hygro_ant3_o_aer

        l = 86
        mastercompindx_ant4_o_aer = l
        name_mastercomp_aer( l ) = 'ant4_o'
        dens_mastercomp_aer( l ) =  dens_ant4_o_aer
        mw_mastercomp_aer(   l ) =    mw_ant4_o_aer
        hygro_mastercomp_aer(l ) = hygro_ant4_o_aer


        l = 87
        mastercompindx_biog1_c_aer = l
        name_mastercomp_aer( l ) = 'biog1_c'
        dens_mastercomp_aer( l ) =  dens_biog1_c_aer
        mw_mastercomp_aer(   l ) =    mw_biog1_c_aer
        hygro_mastercomp_aer(l ) = hygro_biog1_c_aer

        l = 88
        mastercompindx_biog2_c_aer = l
        name_mastercomp_aer( l ) = 'biog2_c'
        dens_mastercomp_aer( l ) =  dens_biog2_c_aer
        mw_mastercomp_aer(   l ) =    mw_biog2_c_aer
        hygro_mastercomp_aer(l ) = hygro_biog2_c_aer

        l = 89
        mastercompindx_biog3_c_aer = l
        name_mastercomp_aer( l ) = 'biog3_c'
        dens_mastercomp_aer( l ) =  dens_biog3_c_aer
        mw_mastercomp_aer(   l ) =    mw_biog3_c_aer
        hygro_mastercomp_aer(l ) = hygro_biog3_c_aer

        l = 90
        mastercompindx_biog4_c_aer = l
        name_mastercomp_aer( l ) = 'biog4_c'
        dens_mastercomp_aer( l ) =  dens_biog4_c_aer
        mw_mastercomp_aer(   l ) =    mw_biog4_c_aer
        hygro_mastercomp_aer(l ) = hygro_biog4_c_aer

        l = 91
        mastercompindx_biog1_o_aer = l
        name_mastercomp_aer( l ) = 'biog1_o'
        dens_mastercomp_aer( l ) =  dens_biog1_o_aer
        mw_mastercomp_aer(   l ) =    mw_biog1_o_aer
        hygro_mastercomp_aer(l ) = hygro_biog1_o_aer

        l = 92
        mastercompindx_biog2_o_aer = l
        name_mastercomp_aer( l ) = 'biog2_o'
        dens_mastercomp_aer( l ) =  dens_biog2_o_aer
        mw_mastercomp_aer(   l ) =    mw_biog2_o_aer
        hygro_mastercomp_aer(l ) = hygro_biog2_o_aer

        l = 93
        mastercompindx_biog3_o_aer = l
        name_mastercomp_aer( l ) = 'biog3_o'
        dens_mastercomp_aer( l ) =  dens_biog3_o_aer
        mw_mastercomp_aer(   l ) =    mw_biog3_o_aer
        hygro_mastercomp_aer(l ) = hygro_biog3_o_aer

        l = 94
        mastercompindx_biog4_o_aer = l
        name_mastercomp_aer( l ) = 'biog4_o'
        dens_mastercomp_aer( l ) =  dens_biog4_o_aer
        mw_mastercomp_aer(   l ) =    mw_biog4_o_aer
        hygro_mastercomp_aer(l ) = hygro_biog4_o_aer

        l = 95
        mastercompindx_smpa_aer = l
        name_mastercomp_aer( l ) = 'smpa'
        dens_mastercomp_aer( l ) =  dens_smpa_aer
        mw_mastercomp_aer(   l ) =    mw_smpa_aer
        hygro_mastercomp_aer(l ) = hygro_smpa_aer

        l = 96
        mastercompindx_smpbb_aer = l
        name_mastercomp_aer( l ) = 'smpbb'
        dens_mastercomp_aer( l ) =  dens_smpbb_aer
        mw_mastercomp_aer(   l ) =    mw_smpbb_aer
        hygro_mastercomp_aer(l ) = hygro_smpbb_aer

        l = 97
        mastercompindx_glysoa_r1_aer = l
        name_mastercomp_aer( l ) = 'glysoa_r1'
        dens_mastercomp_aer( l ) =  dens_glysoa_r1_aer
        mw_mastercomp_aer(   l ) =    mw_glysoa_r1_aer
        hygro_mastercomp_aer(l ) = hygro_glysoa_r1_aer

        l = 98
        mastercompindx_glysoa_r2_aer = l
        name_mastercomp_aer( l ) = 'glysoa_r2'
        dens_mastercomp_aer( l ) =  dens_glysoa_r2_aer
        mw_mastercomp_aer(   l ) =    mw_glysoa_r2_aer
        hygro_mastercomp_aer(l ) = hygro_glysoa_r2_aer

        l = 99
        mastercompindx_glysoa_sfc_aer = l
        name_mastercomp_aer( l ) = 'glysoa_sfc'
        dens_mastercomp_aer( l ) =  dens_glysoa_sfc_aer
        mw_mastercomp_aer(   l ) =    mw_glysoa_sfc_aer
        hygro_mastercomp_aer(l ) = hygro_glysoa_sfc_aer

        l = 100
        mastercompindx_glysoa_nh4_aer = l
        name_mastercomp_aer( l ) = 'glysoa_nh4'
        dens_mastercomp_aer( l ) =  dens_glysoa_nh4_aer
        mw_mastercomp_aer(   l ) =    mw_glysoa_nh4_aer
        hygro_mastercomp_aer(l ) = hygro_glysoa_nh4_aer

        l = 101
        mastercompindx_glysoa_oh_aer = l
        name_mastercomp_aer( l ) = 'glysoa_oh'
        dens_mastercomp_aer( l ) =  dens_glysoa_oh_aer
        mw_mastercomp_aer(   l ) =    mw_glysoa_oh_aer
        hygro_mastercomp_aer(l ) = hygro_glysoa_oh_aer

        l = 102
        mastercompindx_asoaX_aer = l
        name_mastercomp_aer( l ) = 'asoaX'
        dens_mastercomp_aer( l ) =  dens_asoaX_aer
        mw_mastercomp_aer(   l ) =    mw_asoaX_aer
        hygro_mastercomp_aer(l ) = hygro_asoaX_aer

        l = 103
        mastercompindx_asoa1_aer = l
        name_mastercomp_aer( l ) = 'asoa1'
        dens_mastercomp_aer( l ) =  dens_asoa1_aer
        mw_mastercomp_aer(   l ) =    mw_asoa1_aer
        hygro_mastercomp_aer(l ) = hygro_asoa1_aer

        l = 104
        mastercompindx_asoa2_aer = l
        name_mastercomp_aer( l ) = 'asoa2'
        dens_mastercomp_aer( l ) =  dens_asoa2_aer
        mw_mastercomp_aer(   l ) =    mw_asoa2_aer
        hygro_mastercomp_aer(l ) = hygro_asoa2_aer

        l = 105
        mastercompindx_asoa3_aer = l
        name_mastercomp_aer( l ) = 'asoa3'
        dens_mastercomp_aer( l ) =  dens_asoa3_aer
        mw_mastercomp_aer(   l ) =    mw_asoa3_aer
        hygro_mastercomp_aer(l ) = hygro_asoa3_aer

        l = 106
        mastercompindx_asoa4_aer = l
        name_mastercomp_aer( l ) = 'asoa4'
        dens_mastercomp_aer( l ) =  dens_asoa4_aer
        mw_mastercomp_aer(   l ) =    mw_asoa4_aer
        hygro_mastercomp_aer(l ) = hygro_asoa4_aer

        l = 107
        mastercompindx_bsoaX_aer = l
        name_mastercomp_aer( l ) = 'bsoaX'
        dens_mastercomp_aer( l ) =  dens_bsoaX_aer
        mw_mastercomp_aer(   l ) =    mw_bsoaX_aer
        hygro_mastercomp_aer(l ) = hygro_bsoaX_aer

        l = 108
        mastercompindx_bsoa1_aer = l
        name_mastercomp_aer( l ) = 'bsoa1'
        dens_mastercomp_aer( l ) =  dens_bsoa1_aer
        mw_mastercomp_aer(   l ) =    mw_bsoa1_aer
        hygro_mastercomp_aer(l ) = hygro_bsoa1_aer

        l = 109
        mastercompindx_bsoa2_aer = l
        name_mastercomp_aer( l ) = 'bsoa2'
        dens_mastercomp_aer( l ) =  dens_bsoa2_aer
        mw_mastercomp_aer(   l ) =    mw_bsoa2_aer
        hygro_mastercomp_aer(l ) = hygro_bsoa2_aer

        l = 110
        mastercompindx_bsoa3_aer = l
        name_mastercomp_aer( l ) = 'bsoa3'
        dens_mastercomp_aer( l ) =  dens_bsoa3_aer
        mw_mastercomp_aer(   l ) =    mw_bsoa3_aer
        hygro_mastercomp_aer(l ) = hygro_bsoa3_aer

        l = 111
        mastercompindx_bsoa4_aer = l
        name_mastercomp_aer( l ) = 'bsoa4'
        dens_mastercomp_aer( l ) =  dens_bsoa4_aer
        mw_mastercomp_aer(   l ) =    mw_bsoa4_aer
        hygro_mastercomp_aer(l ) = hygro_bsoa4_aer

        l = 112
        mastercompindx_iepox_aer = l
        name_mastercomp_aer( l ) = 'iepox'
        dens_mastercomp_aer( l ) =  dens_iepox_aer
        mw_mastercomp_aer(   l ) =    mw_iepox_aer
        hygro_mastercomp_aer(l ) = hygro_iepox_aer

        l = 113
        mastercompindx_iepoxos_aer = l
        name_mastercomp_aer( l ) = 'iepoxos'
        dens_mastercomp_aer( l ) =  dens_iepoxos_aer
        mw_mastercomp_aer(   l ) =    mw_iepoxos_aer
        hygro_mastercomp_aer(l ) = hygro_iepoxos_aer

        l = 114
        mastercompindx_cldiepoxos_aer = l
        name_mastercomp_aer( l ) = 'cldiepoxos'
        dens_mastercomp_aer( l ) =  dens_cldiepoxos_aer
        mw_mastercomp_aer(   l ) =    mw_cldiepoxos_aer
        hygro_mastercomp_aer(l ) = hygro_cldiepoxos_aer

        l = 115
        mastercompindx_tetrol_aer = l
        name_mastercomp_aer( l ) = 'tetrol'
        dens_mastercomp_aer( l ) =  dens_tetrol_aer
        mw_mastercomp_aer(   l ) =    mw_tetrol_aer
        hygro_mastercomp_aer(l ) = hygro_tetrol_aer

        l = 116
        mastercompindx_cldtetrol_aer = l
        name_mastercomp_aer( l ) = 'cldtetrol'
        dens_mastercomp_aer( l ) =  dens_cldtetrol_aer
        mw_mastercomp_aer(   l ) =    mw_cldtetrol_aer
        hygro_mastercomp_aer(l ) = hygro_cldtetrol_aer

        l = 117
        mastercompindx_gly_aer = l
        name_mastercomp_aer( l ) = 'gly'
        dens_mastercomp_aer( l ) =  dens_gly_aer
        mw_mastercomp_aer(   l ) =    mw_gly_aer
        hygro_mastercomp_aer(l ) = hygro_gly_aer

        l = 118
        mastercompindx_tanv_aer = l
        name_mastercomp_aer( l ) = 'tanv'
        dens_mastercomp_aer( l ) =  dens_tanv_aer
        mw_mastercomp_aer(   l ) =    mw_tanv_aer
        hygro_mastercomp_aer(l ) = hygro_tanv_aer

        l = 119
        mastercompindx_sopnv_aer = l
        name_mastercomp_aer( l ) = 'sopnv'
        dens_mastercomp_aer( l ) =  dens_sopnv_aer
        mw_mastercomp_aer(   l ) =    mw_sopnv_aer
        hygro_mastercomp_aer(l ) = hygro_sopnv_aer

        l = 120
        mastercompindx_ternv_aer = l
        name_mastercomp_aer( l ) = 'ternv'
        dens_mastercomp_aer( l ) =  dens_ternv_aer
        mw_mastercomp_aer(   l ) =    mw_ternv_aer
        hygro_mastercomp_aer(l ) = hygro_ternv_aer

        l = 121
        mastercompindx_seqnv_aer = l
        name_mastercomp_aer( l ) = 'seqnv'
        dens_mastercomp_aer( l ) =  dens_seqnv_aer
        mw_mastercomp_aer(   l ) =    mw_seqnv_aer
        hygro_mastercomp_aer(l ) = hygro_seqnv_aer

        l = 122
        mastercompindx_sianv_aer = l
        name_mastercomp_aer( l ) = 'sianv'
        dens_mastercomp_aer( l ) =  dens_sianv_aer
        mw_mastercomp_aer(   l ) =    mw_sianv_aer
        hygro_mastercomp_aer(l ) = hygro_sianv_aer

	
	l = 123
        mastercompindx_terp_Cn3_aer = l
        name_mastercomp_aer( l ) = 'terp_Cn3'
        dens_mastercomp_aer( l ) =  dens_terp_Cn3_aer
        mw_mastercomp_aer(   l ) =    mw_terp_Cn3_aer
        hygro_mastercomp_aer(l ) = hygro_terp_Cn3_aer




        do itype = 1, ntype_aer
	    nhi = nsize_aer(itype)


            if (nhi == 20) then
               dlo_sect(1,itype) = 0.9688727e-7
            else if (nhi == 12) then
               dlo_sect(1,itype) = 2.441406e-7
            else if ((nhi == 8) .or. (nhi == 4)) then
               dlo_sect(1,itype) = 39.0625e-7
            else
               call peg_error_fatal( lunerr,   &
                   'size bins must be 4, 8, 12, or 20' )
            end if

	    dhi_sect(nhi,itype) = 10.0e-4

	    dum = log( dhi_sect(nhi,itype)/dlo_sect(1,itype) ) / nhi 
	    do n = 2, nhi
		dlo_sect(n,itype) = dlo_sect(1,itype) * exp( (n-1)*dum )
		dhi_sect(n-1,itype) = dlo_sect(n,itype)
	    end do
	    do n = 1, nhi
		dcen_sect(n,itype) = sqrt( dlo_sect(n,itype)*dhi_sect(n,itype) )
		volumlo_sect(n,itype) = (pi/6.) * (dlo_sect(n,itype)**3)
		volumhi_sect(n,itype) = (pi/6.) * (dhi_sect(n,itype)**3)
		volumcen_sect(n,itype) = (pi/6.) * (dcen_sect(n,itype)**3)
		sigmag_aer(n,itype) = (dhi_sect(n,itype)/dlo_sect(n,itype))**0.289
	    end do
	end do




	call init_data_mosaic_ptr( id,is_aerosol )




	call init_csuesat





	call move_sections(    -1,    1,    1, 1, 1 )

	call test_move_sections( 1,   1,    1, 1, 1 )
    

        idum = mod( max( config_flags%mosaic_aerchem_optaa, 0 ), 10 )   
        if ( idum == 2 ) then
            call init_data_mosaic2_asect( id, config_flags, is_aerosol, vbs_nbin )
        end if


	end subroutine init_data_mosaic_asect



	subroutine init_data_mosaic_ptr( id,is_aerosol )

	use module_configure
	use module_state_description, only:  param_first_scalar,num_chem
	use module_data_mosaic_asect
	use module_data_mosaic_other, only:   &
		kh2so4, khno3, khcl, knh3, ko3, kh2o, ktemp,   &
		kso2, kh2o2, khcho, koh, kho2,   &
		kno3, kno, kno2, khono, kpan,   &
                lmaxd, l2maxd, ltot, ltot2, lunout, lunerr, &
                name,kpcg1_b_c,kpcg2_b_c,kpcg3_b_c,kpcg4_b_c, &
                kpcg5_b_c,kpcg6_b_c,kpcg7_b_c,kpcg8_b_c, &
                kpcg9_b_c,kpcg1_b_o,kpcg2_b_o,kpcg3_b_o, &
                kpcg4_b_o,kpcg5_b_o,kpcg6_b_o,kpcg7_b_o, &
                kpcg8_b_o,kpcg9_b_o,kopcg1_b_c,kopcg2_b_c,&
                kopcg3_b_c, kopcg4_b_c,kopcg5_b_c,kopcg6_b_c,&
                kopcg7_b_c,kopcg8_b_c,kopcg1_b_o,kopcg2_b_o,&
                kopcg3_b_o,kopcg4_b_o,kopcg5_b_o,kopcg6_b_o,&
                kopcg7_b_o,kopcg8_b_o,&
                kpcg1_f_c,kpcg2_f_c,kpcg3_f_c,kpcg4_f_c, &
                kpcg5_f_c,kpcg6_f_c,kpcg7_f_c,kpcg8_f_c, &
                kpcg9_f_c,kpcg1_f_o,kpcg2_f_o,kpcg3_f_o, &
                kpcg4_f_o,kpcg5_f_o,kpcg6_f_o,kpcg7_f_o, &
                kpcg8_f_o,kpcg9_f_o,kopcg1_f_c,kopcg2_f_c,&
                kopcg3_f_c, kopcg4_f_c,kopcg5_f_c,kopcg6_f_c,&
                kopcg7_f_c,kopcg8_f_c,kopcg1_f_o,kopcg2_f_o,&
                kopcg3_f_o,kopcg4_f_o,kopcg5_f_o,kopcg6_f_o,&
                kopcg7_f_o,kopcg8_f_o,                      &
                ksmpa,ksmpbb,                         &
                kiepox,kgly,ktetrol,kcldtetrol,               &
                kant1_c,kant2_c,kant3_c,kant4_c,kant1_o,kant2_o, &
                kant3_o,kant4_o,                            &
                kbiog1_c,kbiog2_c,kbiog3_c,kbiog4_c,kbiog1_o,kbiog2_o, &
                kbiog3_o,kbiog4_o, &

                kn2o5, kclno2, &                            
                kasoaX,kasoa1,kasoa2,kasoa3,kasoa4,&
                kbsoaX,kbsoa1,kbsoa2,kbsoa3,kbsoa4,&
		kterp_g_Cn3 

	use module_peg_util, only:  peg_error_fatal, peg_message
	use module_mosaic_wetscav, only: initwet
    use module_scalar_tables, only:  chem_dname_table
	implicit none


        logical, intent(out) :: is_aerosol(num_chem)
        integer, intent(in) :: id

	integer l, ll, n, p1st
	integer iaddto_ncomp, iaddto_ncomp_plustracer
	integer l_mastercomp, lptr_dum
	integer mcindx_dum
	integer isize, itype, iphase
	integer nphasetxt, nsizetxt, nspectxt, ntypetxt
	integer ncomp_dum(maxd_asize,maxd_aphase)
	integer ncomp_plustracer_dum(maxd_asize,maxd_aphase)

	integer y_so4, y_no3, y_cl, y_msa, y_co3, y_nh4, y_na,   &
		y_ca, y_oin, y_oc, y_bc, y_hysw, y_water, &
                y_num, &
      y_pcg1_b_c,y_pcg2_b_c,y_pcg3_b_c,y_pcg4_b_c, &
      y_pcg5_b_c,y_pcg6_b_c,y_pcg7_b_c,y_pcg8_b_c, &
      y_pcg9_b_c,y_pcg1_b_o,y_pcg2_b_o,y_pcg3_b_o, &
      y_pcg4_b_o,y_pcg5_b_o,y_pcg6_b_o,y_pcg7_b_o, &
      y_pcg8_b_o,y_pcg9_b_o,y_opcg1_b_c,y_opcg2_b_c,&
      y_opcg3_b_c, y_opcg4_b_c,y_opcg5_b_c,y_opcg6_b_c,&
      y_opcg7_b_c,y_opcg8_b_c,y_opcg1_b_o,y_opcg2_b_o,&
      y_opcg3_b_o,y_opcg4_b_o,y_opcg5_b_o,y_opcg6_b_o,&
      y_opcg7_b_o,y_opcg8_b_o,&
      y_pcg1_f_c,y_pcg2_f_c,y_pcg3_f_c,y_pcg4_f_c, &
      y_pcg5_f_c,y_pcg6_f_c,y_pcg7_f_c,y_pcg8_f_c, &
      y_pcg9_f_c,y_pcg1_f_o,y_pcg2_f_o,y_pcg3_f_o, &
      y_pcg4_f_o,y_pcg5_f_o,y_pcg6_f_o,y_pcg7_f_o, &
      y_pcg8_f_o,y_pcg9_f_o,y_opcg1_f_c,y_opcg2_f_c,&
      y_opcg3_f_c, y_opcg4_f_c,y_opcg5_f_c,y_opcg6_f_c,&
      y_opcg7_f_c,y_opcg8_f_c,y_opcg1_f_o,y_opcg2_f_o,&
      y_opcg3_f_o,y_opcg4_f_o,y_opcg5_f_o,y_opcg6_f_o,&
      y_opcg7_f_o,y_opcg8_f_o, &
      y_smpa,y_smpbb, &
      y_iepox,y_iepoxos,y_tetrol,y_cldiepoxos,y_cldtetrol,y_gly, &
      y_glysoa_r1, y_glysoa_r2, y_glysoa_sfc, y_glysoa_nh4, y_glysoa_oh, &
      y_ant1_c,y_ant2_c,y_ant3_c,y_ant4_c, &
      y_ant1_o,y_ant2_o,y_ant3_o,y_ant4_o, &
      y_biog1_c,y_biog2_c,y_biog3_c,y_biog4_c, &
      y_biog1_o,y_biog2_o,y_biog3_o,y_biog4_o, &
      y_tanv,y_sopnv,y_ternv,y_seqnv,y_sianv, &
      y_asoaX, y_asoa1, y_asoa2, y_asoa3, y_asoa4, &
      y_bsoaX, y_bsoa1, y_bsoa2, y_bsoa3, y_bsoa4, &
      y_terp_Cn3 

	integer y_cw_so4, y_cw_no3, y_cw_cl, y_cw_msa, y_cw_co3,   &
		y_cw_nh4, y_cw_na,   &
		y_cw_ca, y_cw_oin, y_cw_oc, y_cw_bc, y_cw_num

	character*200 msg
	character*20 phasetxt, sizetxt, spectxt, typetxt


	p1st = param_first_scalar




	itype=1
	lptr_so4_aer(:,itype,:)      = 1
	lptr_no3_aer(:,itype,:)      = 1
	lptr_cl_aer(:,itype,:)       = 1
	lptr_msa_aer(:,itype,:)      = 1
	lptr_co3_aer(:,itype,:)      = 1
	lptr_nh4_aer(:,itype,:)      = 1
	lptr_na_aer(:,itype,:)       = 1
	lptr_ca_aer(:,itype,:)       = 1
	lptr_oin_aer(:,itype,:)      = 1
	lptr_oc_aer(:,itype,:)       = 1
	lptr_bc_aer(:,itype,:)       = 1
	hyswptr_aer(:,itype)    = 1
	waterptr_aer(:,itype)        = 1
	numptr_aer(:,itype,:)        = 1
        lptr_tanv_aer(:,itype,:)       = 1
        lptr_sopnv_aer(:,itype,:)       = 1
        lptr_ternv_aer(:,itype,:)       = 1
        lptr_seqnv_aer(:,itype,:)       = 1
        lptr_sianv_aer(:,itype,:)       = 1
        lptr_pcg1_b_c_aer(:,itype,:)      = 1
        lptr_pcg2_b_c_aer(:,itype,:)      = 1
        lptr_pcg3_b_c_aer(:,itype,:)      = 1
        lptr_pcg4_b_c_aer(:,itype,:)      = 1
        lptr_pcg5_b_c_aer(:,itype,:)      = 1
        lptr_pcg6_b_c_aer(:,itype,:)      = 1
        lptr_pcg7_b_c_aer(:,itype,:)      = 1
        lptr_pcg8_b_c_aer(:,itype,:)      = 1
        lptr_pcg9_b_c_aer(:,itype,:)      = 1
        lptr_pcg1_b_o_aer(:,itype,:)      = 1
        lptr_pcg2_b_o_aer(:,itype,:)      = 1
        lptr_pcg3_b_o_aer(:,itype,:)      = 1
        lptr_pcg4_b_o_aer(:,itype,:)      = 1
        lptr_pcg5_b_o_aer(:,itype,:)      = 1
        lptr_pcg6_b_o_aer(:,itype,:)      = 1
        lptr_pcg7_b_o_aer(:,itype,:)      = 1
        lptr_pcg8_b_o_aer(:,itype,:)      = 1
        lptr_pcg9_b_o_aer(:,itype,:)      = 1
        lptr_opcg1_b_c_aer(:,itype,:)      = 1
        lptr_opcg2_b_c_aer(:,itype,:)      = 1
        lptr_opcg3_b_c_aer(:,itype,:)      = 1
        lptr_opcg4_b_c_aer(:,itype,:)      = 1
        lptr_opcg5_b_c_aer(:,itype,:)      = 1
        lptr_opcg6_b_c_aer(:,itype,:)      = 1
        lptr_opcg7_b_c_aer(:,itype,:)      = 1
        lptr_opcg8_b_c_aer(:,itype,:)      = 1
        lptr_opcg1_b_o_aer(:,itype,:)      = 1
        lptr_opcg2_b_o_aer(:,itype,:)      = 1
        lptr_opcg3_b_o_aer(:,itype,:)      = 1
        lptr_opcg4_b_o_aer(:,itype,:)      = 1
        lptr_opcg5_b_o_aer(:,itype,:)      = 1
        lptr_opcg6_b_o_aer(:,itype,:)      = 1
        lptr_opcg7_b_o_aer(:,itype,:)      = 1
        lptr_opcg8_b_o_aer(:,itype,:)      = 1
        lptr_pcg1_f_c_aer(:,itype,:)      = 1
        lptr_pcg2_f_c_aer(:,itype,:)      = 1
        lptr_pcg3_f_c_aer(:,itype,:)      = 1
        lptr_pcg4_f_c_aer(:,itype,:)      = 1
        lptr_pcg5_f_c_aer(:,itype,:)      = 1
        lptr_pcg6_f_c_aer(:,itype,:)      = 1
        lptr_pcg7_f_c_aer(:,itype,:)      = 1
        lptr_pcg8_f_c_aer(:,itype,:)      = 1
        lptr_pcg9_f_c_aer(:,itype,:)      = 1
        lptr_pcg1_f_o_aer(:,itype,:)      = 1
        lptr_pcg2_f_o_aer(:,itype,:)      = 1
        lptr_pcg3_f_o_aer(:,itype,:)      = 1
        lptr_pcg4_f_o_aer(:,itype,:)      = 1
        lptr_pcg5_f_o_aer(:,itype,:)      = 1
        lptr_pcg6_f_o_aer(:,itype,:)      = 1
        lptr_pcg7_f_o_aer(:,itype,:)      = 1
        lptr_pcg8_f_o_aer(:,itype,:)      = 1
        lptr_pcg9_f_o_aer(:,itype,:)      = 1
        lptr_opcg1_f_c_aer(:,itype,:)      = 1
        lptr_opcg2_f_c_aer(:,itype,:)      = 1
        lptr_opcg3_f_c_aer(:,itype,:)      = 1
        lptr_opcg4_f_c_aer(:,itype,:)      = 1
        lptr_opcg5_f_c_aer(:,itype,:)      = 1
        lptr_opcg6_f_c_aer(:,itype,:)      = 1
        lptr_opcg7_f_c_aer(:,itype,:)      = 1
        lptr_opcg8_f_c_aer(:,itype,:)      = 1
        lptr_opcg1_f_o_aer(:,itype,:)      = 1
        lptr_opcg2_f_o_aer(:,itype,:)      = 1
        lptr_opcg3_f_o_aer(:,itype,:)      = 1
        lptr_opcg4_f_o_aer(:,itype,:)      = 1
        lptr_opcg5_f_o_aer(:,itype,:)      = 1
        lptr_opcg6_f_o_aer(:,itype,:)      = 1
        lptr_opcg7_f_o_aer(:,itype,:)      = 1
        lptr_opcg8_f_o_aer(:,itype,:)      = 1
        lptr_smpa_aer(:,itype,:)      = 1
        lptr_smpbb_aer(:,itype,:)      = 1
        lptr_iepox_aer(:,itype,:)          = 1
        lptr_iepoxos_aer(:,itype,:)        = 1
        lptr_tetrol_aer(:,itype,:)         = 1
        lptr_cldiepoxos_aer(:,itype,:)        = 1
        lptr_cldtetrol_aer(:,itype,:)         = 1
        lptr_gly_aer(:,itype,:)            = 1
        lptr_glysoa_r1_aer(:,itype,:)  = 1
        lptr_glysoa_r2_aer(:,itype,:)  = 1
        lptr_glysoa_sfc_aer(:,itype,:) = 1
        lptr_glysoa_nh4_aer(:,itype,:) = 1
        lptr_glysoa_oh_aer(:,itype,:)  = 1
        lptr_ant1_c_aer(:,itype,:)      = 1
        lptr_ant2_c_aer(:,itype,:)      = 1
        lptr_ant3_c_aer(:,itype,:)      = 1
        lptr_ant4_c_aer(:,itype,:)      = 1
        lptr_ant1_o_aer(:,itype,:)      = 1
        lptr_ant2_o_aer(:,itype,:)      = 1
        lptr_ant3_o_aer(:,itype,:)      = 1
        lptr_ant4_o_aer(:,itype,:)      = 1
        lptr_biog1_c_aer(:,itype,:)      = 1
        lptr_biog2_c_aer(:,itype,:)      = 1
        lptr_biog3_c_aer(:,itype,:)      = 1
        lptr_biog4_c_aer(:,itype,:)      = 1
        lptr_biog1_o_aer(:,itype,:)      = 1
        lptr_biog2_o_aer(:,itype,:)      = 1
        lptr_biog3_o_aer(:,itype,:)      = 1
        lptr_biog4_o_aer(:,itype,:)      = 1
        lptr_asoaX_aer(:,itype,:)        = 1
        lptr_asoa1_aer(:,itype,:)        = 1
        lptr_asoa2_aer(:,itype,:)        = 1
        lptr_asoa3_aer(:,itype,:)        = 1
        lptr_asoa4_aer(:,itype,:)        = 1
        lptr_bsoaX_aer(:,itype,:)        = 1
        lptr_bsoa1_aer(:,itype,:)        = 1
        lptr_bsoa2_aer(:,itype,:)        = 1
        lptr_bsoa3_aer(:,itype,:)        = 1
        lptr_bsoa4_aer(:,itype,:)        = 1
        
	lptr_terp_Cn3_aer(:,itype,:)      = 1


	if (nsize_aer(itype) .ge. 1) then
	if (p_so4_a01 .ge. p1st) lptr_so4_aer(01,itype,ai_phase) = p_so4_a01
	if (p_no3_a01 .ge. p1st) lptr_no3_aer(01,itype,ai_phase) = p_no3_a01
	if (p_cl_a01 .ge. p1st) lptr_cl_aer(01,itype,ai_phase)  = p_cl_a01
	if (p_msa_a01 .ge. p1st) lptr_msa_aer(01,itype,ai_phase) = p_msa_a01
	if (p_co3_a01 .ge. p1st) lptr_co3_aer(01,itype,ai_phase) = p_co3_a01
	if (p_nh4_a01 .ge. p1st) lptr_nh4_aer(01,itype,ai_phase) = p_nh4_a01
	if (p_na_a01 .ge. p1st) lptr_na_aer(01,itype,ai_phase)  = p_na_a01
	if (p_ca_a01 .ge. p1st) lptr_ca_aer(01,itype,ai_phase)  = p_ca_a01
	if (p_oin_a01 .ge. p1st) lptr_oin_aer(01,itype,ai_phase) = p_oin_a01
	if (p_oc_a01 .ge. p1st) lptr_oc_aer(01,itype,ai_phase)  = p_oc_a01
	if (p_bc_a01 .ge. p1st) lptr_bc_aer(01,itype,ai_phase)  = p_bc_a01
	if (p_hysw_a01 .ge. p1st) hyswptr_aer(01,itype)         = p_hysw_a01
	if (p_water_a01 .ge. p1st) waterptr_aer(01,itype)       = p_water_a01
        if (p_tanv_a01 .ge. p1st) lptr_tanv_aer(01,itype,ai_phase)  = p_tanv_a01
        if (p_sopnv_a01 .ge. p1st) lptr_sopnv_aer(01,itype,ai_phase)  = p_sopnv_a01
        if (p_ternv_a01 .ge. p1st) lptr_ternv_aer(01,itype,ai_phase)  = p_ternv_a01
        if (p_seqnv_a01 .ge. p1st) lptr_seqnv_aer(01,itype,ai_phase)  = p_seqnv_a01
        if (p_sianv_a01 .ge. p1st) lptr_sianv_aer(01,itype,ai_phase)  = p_sianv_a01
        if (p_pcg1_b_c_a01 .ge. p1st) lptr_pcg1_b_c_aer(01,itype,ai_phase) = p_pcg1_b_c_a01
        if (p_pcg2_b_c_a01 .ge. p1st) lptr_pcg2_b_c_aer(01,itype,ai_phase)  = p_pcg2_b_c_a01
        if (p_pcg3_b_c_a01 .ge. p1st) lptr_pcg3_b_c_aer(01,itype,ai_phase)  = p_pcg3_b_c_a01
        if (p_pcg4_b_c_a01 .ge. p1st) lptr_pcg4_b_c_aer(01,itype,ai_phase)  = p_pcg4_b_c_a01
        if (p_pcg5_b_c_a01 .ge. p1st) lptr_pcg5_b_c_aer(01,itype,ai_phase)  = p_pcg5_b_c_a01
        if (p_pcg6_b_c_a01 .ge. p1st) lptr_pcg6_b_c_aer(01,itype,ai_phase)  = p_pcg6_b_c_a01
        if (p_pcg7_b_c_a01 .ge. p1st) lptr_pcg7_b_c_aer(01,itype,ai_phase)  = p_pcg7_b_c_a01
        if (p_pcg8_b_c_a01 .ge. p1st) lptr_pcg8_b_c_aer(01,itype,ai_phase)  = p_pcg8_b_c_a01
        if (p_pcg9_b_c_a01 .ge. p1st) lptr_pcg9_b_c_aer(01,itype,ai_phase)  = p_pcg9_b_c_a01
        if (p_pcg1_b_o_a01 .ge. p1st) lptr_pcg1_b_o_aer(01,itype,ai_phase) = p_pcg1_b_o_a01
        if (p_pcg2_b_o_a01 .ge. p1st) lptr_pcg2_b_o_aer(01,itype,ai_phase)  = p_pcg2_b_o_a01
        if (p_pcg3_b_o_a01 .ge. p1st) lptr_pcg3_b_o_aer(01,itype,ai_phase)  = p_pcg3_b_o_a01
        if (p_pcg4_b_o_a01 .ge. p1st) lptr_pcg4_b_o_aer(01,itype,ai_phase)  = p_pcg4_b_o_a01
        if (p_pcg5_b_o_a01 .ge. p1st) lptr_pcg5_b_o_aer(01,itype,ai_phase)  = p_pcg5_b_o_a01
        if (p_pcg6_b_o_a01 .ge. p1st) lptr_pcg6_b_o_aer(01,itype,ai_phase)  = p_pcg6_b_o_a01
        if (p_pcg7_b_o_a01 .ge. p1st) lptr_pcg7_b_o_aer(01,itype,ai_phase)  = p_pcg7_b_o_a01
        if (p_pcg8_b_o_a01 .ge. p1st) lptr_pcg8_b_o_aer(01,itype,ai_phase)  = p_pcg8_b_o_a01
        if (p_pcg9_b_o_a01 .ge. p1st) lptr_pcg9_b_o_aer(01,itype,ai_phase)  = p_pcg9_b_o_a01
        if (p_opcg1_b_c_a01 .ge. p1st) lptr_opcg1_b_c_aer(01,itype,ai_phase) = p_opcg1_b_c_a01
        if (p_opcg2_b_c_a01 .ge. p1st) lptr_opcg2_b_c_aer(01,itype,ai_phase)  = p_opcg2_b_c_a01
        if (p_opcg3_b_c_a01 .ge. p1st) lptr_opcg3_b_c_aer(01,itype,ai_phase)  = p_opcg3_b_c_a01
        if (p_opcg4_b_c_a01 .ge. p1st) lptr_opcg4_b_c_aer(01,itype,ai_phase)  = p_opcg4_b_c_a01
        if (p_opcg5_b_c_a01 .ge. p1st) lptr_opcg5_b_c_aer(01,itype,ai_phase)  = p_opcg5_b_c_a01
        if (p_opcg6_b_c_a01 .ge. p1st) lptr_opcg6_b_c_aer(01,itype,ai_phase)  = p_opcg6_b_c_a01
        if (p_opcg7_b_c_a01 .ge. p1st) lptr_opcg7_b_c_aer(01,itype,ai_phase)  = p_opcg7_b_c_a01
        if (p_opcg8_b_c_a01 .ge. p1st) lptr_opcg8_b_c_aer(01,itype,ai_phase)  = p_opcg8_b_c_a01
        if (p_opcg1_b_o_a01 .ge. p1st) lptr_opcg1_b_o_aer(01,itype,ai_phase) = p_opcg1_b_o_a01
        if (p_opcg2_b_o_a01 .ge. p1st) lptr_opcg2_b_o_aer(01,itype,ai_phase)  = p_opcg2_b_o_a01
        if (p_opcg3_b_o_a01 .ge. p1st) lptr_opcg3_b_o_aer(01,itype,ai_phase)  = p_opcg3_b_o_a01
        if (p_opcg4_b_o_a01 .ge. p1st) lptr_opcg4_b_o_aer(01,itype,ai_phase)  = p_opcg4_b_o_a01
        if (p_opcg5_b_o_a01 .ge. p1st) lptr_opcg5_b_o_aer(01,itype,ai_phase)  = p_opcg5_b_o_a01
        if (p_opcg6_b_o_a01 .ge. p1st) lptr_opcg6_b_o_aer(01,itype,ai_phase)  = p_opcg6_b_o_a01
        if (p_opcg7_b_o_a01 .ge. p1st) lptr_opcg7_b_o_aer(01,itype,ai_phase)  = p_opcg7_b_o_a01
        if (p_opcg8_b_o_a01 .ge. p1st) lptr_opcg8_b_o_aer(01,itype,ai_phase)  = p_opcg8_b_o_a01
        if (p_pcg1_f_c_a01 .ge. p1st) lptr_pcg1_f_c_aer(01,itype,ai_phase) = p_pcg1_f_c_a01
        if (p_pcg2_f_c_a01 .ge. p1st) lptr_pcg2_f_c_aer(01,itype,ai_phase)  = p_pcg2_f_c_a01
        if (p_pcg3_f_c_a01 .ge. p1st) lptr_pcg3_f_c_aer(01,itype,ai_phase)  = p_pcg3_f_c_a01
        if (p_pcg4_f_c_a01 .ge. p1st) lptr_pcg4_f_c_aer(01,itype,ai_phase)  = p_pcg4_f_c_a01
        if (p_pcg5_f_c_a01 .ge. p1st) lptr_pcg5_f_c_aer(01,itype,ai_phase)  = p_pcg5_f_c_a01
        if (p_pcg6_f_c_a01 .ge. p1st) lptr_pcg6_f_c_aer(01,itype,ai_phase)  = p_pcg6_f_c_a01
        if (p_pcg7_f_c_a01 .ge. p1st) lptr_pcg7_f_c_aer(01,itype,ai_phase)  = p_pcg7_f_c_a01
        if (p_pcg8_f_c_a01 .ge. p1st) lptr_pcg8_f_c_aer(01,itype,ai_phase)  = p_pcg8_f_c_a01
        if (p_pcg9_f_c_a01 .ge. p1st) lptr_pcg9_f_c_aer(01,itype,ai_phase)  = p_pcg9_f_c_a01
        if (p_pcg1_f_o_a01 .ge. p1st) lptr_pcg1_f_o_aer(01,itype,ai_phase) = p_pcg1_f_o_a01
        if (p_pcg2_f_o_a01 .ge. p1st) lptr_pcg2_f_o_aer(01,itype,ai_phase)  = p_pcg2_f_o_a01
        if (p_pcg3_f_o_a01 .ge. p1st) lptr_pcg3_f_o_aer(01,itype,ai_phase)  = p_pcg3_f_o_a01
        if (p_pcg4_f_o_a01 .ge. p1st) lptr_pcg4_f_o_aer(01,itype,ai_phase)  = p_pcg4_f_o_a01
        if (p_pcg5_f_o_a01 .ge. p1st) lptr_pcg5_f_o_aer(01,itype,ai_phase)  = p_pcg5_f_o_a01
        if (p_pcg6_f_o_a01 .ge. p1st) lptr_pcg6_f_o_aer(01,itype,ai_phase)  = p_pcg6_f_o_a01
        if (p_pcg7_f_o_a01 .ge. p1st) lptr_pcg7_f_o_aer(01,itype,ai_phase)  = p_pcg7_f_o_a01
        if (p_pcg8_f_o_a01 .ge. p1st) lptr_pcg8_f_o_aer(01,itype,ai_phase)  = p_pcg8_f_o_a01
        if (p_pcg9_f_o_a01 .ge. p1st) lptr_pcg9_f_o_aer(01,itype,ai_phase)  = p_pcg9_f_o_a01
        if (p_opcg1_f_c_a01 .ge. p1st) lptr_opcg1_f_c_aer(01,itype,ai_phase) = p_opcg1_f_c_a01
        if (p_opcg2_f_c_a01 .ge. p1st) lptr_opcg2_f_c_aer(01,itype,ai_phase)  = p_opcg2_f_c_a01
        if (p_opcg3_f_c_a01 .ge. p1st) lptr_opcg3_f_c_aer(01,itype,ai_phase)  = p_opcg3_f_c_a01
        if (p_opcg4_f_c_a01 .ge. p1st) lptr_opcg4_f_c_aer(01,itype,ai_phase)  = p_opcg4_f_c_a01
        if (p_opcg5_f_c_a01 .ge. p1st) lptr_opcg5_f_c_aer(01,itype,ai_phase)  = p_opcg5_f_c_a01
        if (p_opcg6_f_c_a01 .ge. p1st) lptr_opcg6_f_c_aer(01,itype,ai_phase)  = p_opcg6_f_c_a01
        if (p_opcg7_f_c_a01 .ge. p1st) lptr_opcg7_f_c_aer(01,itype,ai_phase)  = p_opcg7_f_c_a01
        if (p_opcg8_f_c_a01 .ge. p1st) lptr_opcg8_f_c_aer(01,itype,ai_phase)  = p_opcg8_f_c_a01
        if (p_opcg1_f_o_a01 .ge. p1st) lptr_opcg1_f_o_aer(01,itype,ai_phase) = p_opcg1_f_o_a01
        if (p_opcg2_f_o_a01 .ge. p1st) lptr_opcg2_f_o_aer(01,itype,ai_phase)  = p_opcg2_f_o_a01
        if (p_opcg3_f_o_a01 .ge. p1st) lptr_opcg3_f_o_aer(01,itype,ai_phase)  = p_opcg3_f_o_a01
        if (p_opcg4_f_o_a01 .ge. p1st) lptr_opcg4_f_o_aer(01,itype,ai_phase)  = p_opcg4_f_o_a01
        if (p_opcg5_f_o_a01 .ge. p1st) lptr_opcg5_f_o_aer(01,itype,ai_phase)  = p_opcg5_f_o_a01
        if (p_opcg6_f_o_a01 .ge. p1st) lptr_opcg6_f_o_aer(01,itype,ai_phase)  = p_opcg6_f_o_a01
        if (p_opcg7_f_o_a01 .ge. p1st) lptr_opcg7_f_o_aer(01,itype,ai_phase)  = p_opcg7_f_o_a01
        if (p_opcg8_f_o_a01 .ge. p1st) lptr_opcg8_f_o_aer(01,itype,ai_phase)  = p_opcg8_f_o_a01
        if (p_smpa_a01 .ge. p1st) lptr_smpa_aer(01,itype,ai_phase)  = p_smpa_a01
        if (p_smpbb_a01 .ge. p1st) lptr_smpbb_aer(01,itype,ai_phase)  = p_smpbb_a01
        if (p_iepox_a01 .ge. p1st) lptr_iepox_aer(01,itype,ai_phase)    =p_iepox_a01
        if (p_iepoxos_a01 .ge. p1st) lptr_iepoxos_aer(01,itype,ai_phase)  = p_iepoxos_a01
        if (p_tetrol_a01 .ge. p1st) lptr_tetrol_aer(01,itype,ai_phase)  = p_tetrol_a01
        if (p_cldiepoxos_a01 .ge. p1st) lptr_cldiepoxos_aer(01,itype,ai_phase)  =p_cldiepoxos_a01
        if (p_cldtetrol_a01 .ge. p1st) lptr_cldtetrol_aer(01,itype,ai_phase)  =p_cldtetrol_a01

        if (p_gly_a01 .ge. p1st) lptr_gly_aer(01,itype,ai_phase)       = p_gly_a01
        if (p_glysoa_r1_a01 .ge. p1st) lptr_glysoa_r1_aer (01,itype,ai_phase) = p_glysoa_r1_a01
        if (p_glysoa_r2_a01 .ge. p1st) lptr_glysoa_r2_aer (01,itype,ai_phase) = p_glysoa_r2_a01
        if (p_glysoa_sfc_a01 .ge. p1st) lptr_glysoa_sfc_aer (01,itype,ai_phase) = p_glysoa_sfc_a01
        if (p_glysoa_nh4_a01 .ge. p1st) lptr_glysoa_nh4_aer (01,itype,ai_phase) = p_glysoa_nh4_a01
        if (p_glysoa_oh_a01 .ge. p1st) lptr_glysoa_oh_aer (01,itype,ai_phase) = p_glysoa_oh_a01
        if (p_ant1_c_a01 .ge. p1st) lptr_ant1_c_aer(01,itype,ai_phase)  = p_ant1_c_a01
        if (p_ant2_c_a01 .ge. p1st) lptr_ant2_c_aer(01,itype,ai_phase)  = p_ant2_c_a01
        if (p_ant3_c_a01 .ge. p1st) lptr_ant3_c_aer(01,itype,ai_phase)  = p_ant3_c_a01
        if (p_ant4_c_a01 .ge. p1st) lptr_ant4_c_aer(01,itype,ai_phase)  = p_ant4_c_a01
        if (p_ant1_o_a01 .ge. p1st) lptr_ant1_o_aer(01,itype,ai_phase)  = p_ant1_o_a01
        if (p_ant2_o_a01 .ge. p1st) lptr_ant2_o_aer(01,itype,ai_phase)  = p_ant2_o_a01
        if (p_ant3_o_a01 .ge. p1st) lptr_ant3_o_aer(01,itype,ai_phase)  = p_ant3_o_a01
        if (p_ant4_o_a01 .ge. p1st) lptr_ant4_o_aer(01,itype,ai_phase)  = p_ant4_o_a01
        if (p_biog1_c_a01 .ge. p1st) lptr_biog1_c_aer(01,itype,ai_phase)  = p_biog1_c_a01
        if (p_biog2_c_a01 .ge. p1st) lptr_biog2_c_aer(01,itype,ai_phase)  = p_biog2_c_a01
        if (p_biog3_c_a01 .ge. p1st) lptr_biog3_c_aer(01,itype,ai_phase)  = p_biog3_c_a01
        if (p_biog4_c_a01 .ge. p1st) lptr_biog4_c_aer(01,itype,ai_phase)  = p_biog4_c_a01
        if (p_biog1_o_a01 .ge. p1st) lptr_biog1_o_aer(01,itype,ai_phase)  = p_biog1_o_a01
        if (p_biog2_o_a01 .ge. p1st) lptr_biog2_o_aer(01,itype,ai_phase)  = p_biog2_o_a01
        if (p_biog3_o_a01 .ge. p1st) lptr_biog3_o_aer(01,itype,ai_phase)  = p_biog3_o_a01
        if (p_biog4_o_a01 .ge. p1st) lptr_biog4_o_aer(01,itype,ai_phase)  = p_biog4_o_a01
        if (p_asoaX_a01 .ge. p1st) lptr_asoaX_aer(01,itype,ai_phase)  = p_asoaX_a01
        if (p_asoa1_a01 .ge. p1st) lptr_asoa1_aer(01,itype,ai_phase)  = p_asoa1_a01
        if (p_asoa2_a01 .ge. p1st) lptr_asoa2_aer(01,itype,ai_phase)  = p_asoa2_a01
        if (p_asoa3_a01 .ge. p1st) lptr_asoa3_aer(01,itype,ai_phase)  = p_asoa3_a01
        if (p_asoa4_a01 .ge. p1st) lptr_asoa4_aer(01,itype,ai_phase)  = p_asoa4_a01
        if (p_bsoaX_a01 .ge. p1st) lptr_bsoaX_aer(01,itype,ai_phase)  = p_bsoaX_a01
        if (p_bsoa1_a01 .ge. p1st) lptr_bsoa1_aer(01,itype,ai_phase)  = p_bsoa1_a01
        if (p_bsoa2_a01 .ge. p1st) lptr_bsoa2_aer(01,itype,ai_phase)  = p_bsoa2_a01
        if (p_bsoa3_a01 .ge. p1st) lptr_bsoa3_aer(01,itype,ai_phase)  = p_bsoa3_a01
        if (p_bsoa4_a01 .ge. p1st) lptr_bsoa4_aer(01,itype,ai_phase)  = p_bsoa4_a01
	if (p_num_a01 .ge. p1st)  numptr_aer(01,itype,ai_phase)        = p_num_a01

	
	if (p_TERP_CN3_a01 .ge. p1st) lptr_terp_Cn3_aer(01,itype,ai_phase)  = p_TERP_CN3_a01

	end if

	if (nsize_aer(itype) .ge. 2) then
	if (p_so4_a02 .ge. p1st) lptr_so4_aer(02,itype,ai_phase) = p_so4_a02
        if (p_no3_a02 .ge. p1st) lptr_no3_aer(02,itype,ai_phase) = p_no3_a02
        if (p_cl_a02 .ge. p1st) lptr_cl_aer(02,itype,ai_phase)  = p_cl_a02
        if (p_msa_a02 .ge. p1st) lptr_msa_aer(02,itype,ai_phase) = p_msa_a02
        if (p_co3_a02 .ge. p1st) lptr_co3_aer(02,itype,ai_phase) = p_co3_a02
        if (p_nh4_a02 .ge. p1st) lptr_nh4_aer(02,itype,ai_phase) = p_nh4_a02
        if (p_na_a02 .ge. p1st) lptr_na_aer(02,itype,ai_phase)  = p_na_a02
        if (p_ca_a02 .ge. p1st) lptr_ca_aer(02,itype,ai_phase)  = p_ca_a02
        if (p_oin_a02 .ge. p1st) lptr_oin_aer(02,itype,ai_phase) = p_oin_a02
        if (p_oc_a02 .ge. p1st) lptr_oc_aer(02,itype,ai_phase)  = p_oc_a02
        if (p_bc_a02 .ge. p1st) lptr_bc_aer(02,itype,ai_phase)  = p_bc_a02
        if (p_hysw_a02 .ge. p1st) hyswptr_aer(02,itype)         = p_hysw_a02
        if (p_water_a02 .ge. p1st) waterptr_aer(02,itype)       = p_water_a02
        if (p_tanv_a02 .ge. p1st) lptr_tanv_aer(02,itype,ai_phase)  = p_tanv_a02
        if (p_sopnv_a02 .ge. p1st) lptr_sopnv_aer(02,itype,ai_phase)  = p_sopnv_a02
        if (p_ternv_a02 .ge. p1st) lptr_ternv_aer(02,itype,ai_phase)  = p_ternv_a02
        if (p_seqnv_a02 .ge. p1st) lptr_seqnv_aer(02,itype,ai_phase)  = p_seqnv_a02
        if (p_sianv_a02 .ge. p1st) lptr_sianv_aer(02,itype,ai_phase)  = p_sianv_a02
        if (p_pcg1_b_c_a02 .ge. p1st) lptr_pcg1_b_c_aer(02,itype,ai_phase) = p_pcg1_b_c_a02
        if (p_pcg2_b_c_a02 .ge. p1st) lptr_pcg2_b_c_aer(02,itype,ai_phase)  = p_pcg2_b_c_a02
        if (p_pcg3_b_c_a02 .ge. p1st) lptr_pcg3_b_c_aer(02,itype,ai_phase)  = p_pcg3_b_c_a02
        if (p_pcg4_b_c_a02 .ge. p1st) lptr_pcg4_b_c_aer(02,itype,ai_phase)  = p_pcg4_b_c_a02
        if (p_pcg5_b_c_a02 .ge. p1st) lptr_pcg5_b_c_aer(02,itype,ai_phase)  = p_pcg5_b_c_a02
        if (p_pcg6_b_c_a02 .ge. p1st) lptr_pcg6_b_c_aer(02,itype,ai_phase)  = p_pcg6_b_c_a02
        if (p_pcg7_b_c_a02 .ge. p1st) lptr_pcg7_b_c_aer(02,itype,ai_phase)  = p_pcg7_b_c_a02
        if (p_pcg8_b_c_a02 .ge. p1st) lptr_pcg8_b_c_aer(02,itype,ai_phase)  = p_pcg8_b_c_a02
        if (p_pcg9_b_c_a02 .ge. p1st) lptr_pcg9_b_c_aer(02,itype,ai_phase)  = p_pcg9_b_c_a02
        if (p_pcg1_b_o_a02 .ge. p1st) lptr_pcg1_b_o_aer(02,itype,ai_phase) = p_pcg1_b_o_a02
        if (p_pcg2_b_o_a02 .ge. p1st) lptr_pcg2_b_o_aer(02,itype,ai_phase)  = p_pcg2_b_o_a02
        if (p_pcg3_b_o_a02 .ge. p1st) lptr_pcg3_b_o_aer(02,itype,ai_phase)  = p_pcg3_b_o_a02
        if (p_pcg4_b_o_a02 .ge. p1st) lptr_pcg4_b_o_aer(02,itype,ai_phase)  = p_pcg4_b_o_a02
        if (p_pcg5_b_o_a02 .ge. p1st) lptr_pcg5_b_o_aer(02,itype,ai_phase)  = p_pcg5_b_o_a02
        if (p_pcg6_b_o_a02 .ge. p1st) lptr_pcg6_b_o_aer(02,itype,ai_phase)  = p_pcg6_b_o_a02
        if (p_pcg7_b_o_a02 .ge. p1st) lptr_pcg7_b_o_aer(02,itype,ai_phase)  = p_pcg7_b_o_a02
        if (p_pcg8_b_o_a02 .ge. p1st) lptr_pcg8_b_o_aer(02,itype,ai_phase)  = p_pcg8_b_o_a02
        if (p_pcg9_b_o_a02 .ge. p1st) lptr_pcg9_b_o_aer(02,itype,ai_phase)  = p_pcg9_b_o_a02
        if (p_opcg1_b_c_a02 .ge. p1st) lptr_opcg1_b_c_aer(02,itype,ai_phase) = p_opcg1_b_c_a02
        if (p_opcg2_b_c_a02 .ge. p1st) lptr_opcg2_b_c_aer(02,itype,ai_phase)  = p_opcg2_b_c_a02
        if (p_opcg3_b_c_a02 .ge. p1st) lptr_opcg3_b_c_aer(02,itype,ai_phase)  = p_opcg3_b_c_a02
        if (p_opcg4_b_c_a02 .ge. p1st) lptr_opcg4_b_c_aer(02,itype,ai_phase)  = p_opcg4_b_c_a02
        if (p_opcg5_b_c_a02 .ge. p1st) lptr_opcg5_b_c_aer(02,itype,ai_phase)  = p_opcg5_b_c_a02
        if (p_opcg6_b_c_a02 .ge. p1st) lptr_opcg6_b_c_aer(02,itype,ai_phase)  = p_opcg6_b_c_a02
        if (p_opcg7_b_c_a02 .ge. p1st) lptr_opcg7_b_c_aer(02,itype,ai_phase)  = p_opcg7_b_c_a02
        if (p_opcg8_b_c_a02 .ge. p1st) lptr_opcg8_b_c_aer(02,itype,ai_phase)  = p_opcg8_b_c_a02
        if (p_opcg1_b_o_a02 .ge. p1st) lptr_opcg1_b_o_aer(02,itype,ai_phase) = p_opcg1_b_o_a02
        if (p_opcg2_b_o_a02 .ge. p1st) lptr_opcg2_b_o_aer(02,itype,ai_phase)  = p_opcg2_b_o_a02
        if (p_opcg3_b_o_a02 .ge. p1st) lptr_opcg3_b_o_aer(02,itype,ai_phase)  = p_opcg3_b_o_a02
        if (p_opcg4_b_o_a02 .ge. p1st) lptr_opcg4_b_o_aer(02,itype,ai_phase)  = p_opcg4_b_o_a02
        if (p_opcg5_b_o_a02 .ge. p1st) lptr_opcg5_b_o_aer(02,itype,ai_phase)  = p_opcg5_b_o_a02
        if (p_opcg6_b_o_a02 .ge. p1st) lptr_opcg6_b_o_aer(02,itype,ai_phase)  = p_opcg6_b_o_a02
        if (p_opcg7_b_o_a02 .ge. p1st) lptr_opcg7_b_o_aer(02,itype,ai_phase)  = p_opcg7_b_o_a02
        if (p_opcg8_b_o_a02 .ge. p1st) lptr_opcg8_b_o_aer(02,itype,ai_phase)  = p_opcg8_b_o_a02
        if (p_pcg1_f_c_a02 .ge. p1st) lptr_pcg1_f_c_aer(02,itype,ai_phase) = p_pcg1_f_c_a02
        if (p_pcg2_f_c_a02 .ge. p1st) lptr_pcg2_f_c_aer(02,itype,ai_phase)  = p_pcg2_f_c_a02
        if (p_pcg3_f_c_a02 .ge. p1st) lptr_pcg3_f_c_aer(02,itype,ai_phase)  = p_pcg3_f_c_a02
        if (p_pcg4_f_c_a02 .ge. p1st) lptr_pcg4_f_c_aer(02,itype,ai_phase)  = p_pcg4_f_c_a02
        if (p_pcg5_f_c_a02 .ge. p1st) lptr_pcg5_f_c_aer(02,itype,ai_phase)  = p_pcg5_f_c_a02
        if (p_pcg6_f_c_a02 .ge. p1st) lptr_pcg6_f_c_aer(02,itype,ai_phase)  = p_pcg6_f_c_a02
        if (p_pcg7_f_c_a02 .ge. p1st) lptr_pcg7_f_c_aer(02,itype,ai_phase)  = p_pcg7_f_c_a02
        if (p_pcg8_f_c_a02 .ge. p1st) lptr_pcg8_f_c_aer(02,itype,ai_phase)  = p_pcg8_f_c_a02
        if (p_pcg9_f_c_a02 .ge. p1st) lptr_pcg9_f_c_aer(02,itype,ai_phase)  = p_pcg9_f_c_a02
        if (p_pcg1_f_o_a02 .ge. p1st) lptr_pcg1_f_o_aer(02,itype,ai_phase) = p_pcg1_f_o_a02
        if (p_pcg2_f_o_a02 .ge. p1st) lptr_pcg2_f_o_aer(02,itype,ai_phase)  = p_pcg2_f_o_a02
        if (p_pcg3_f_o_a02 .ge. p1st) lptr_pcg3_f_o_aer(02,itype,ai_phase)  = p_pcg3_f_o_a02
        if (p_pcg4_f_o_a02 .ge. p1st) lptr_pcg4_f_o_aer(02,itype,ai_phase)  = p_pcg4_f_o_a02
        if (p_pcg5_f_o_a02 .ge. p1st) lptr_pcg5_f_o_aer(02,itype,ai_phase)  = p_pcg5_f_o_a02
        if (p_pcg6_f_o_a02 .ge. p1st) lptr_pcg6_f_o_aer(02,itype,ai_phase)  = p_pcg6_f_o_a02
        if (p_pcg7_f_o_a02 .ge. p1st) lptr_pcg7_f_o_aer(02,itype,ai_phase)  = p_pcg7_f_o_a02
        if (p_pcg8_f_o_a02 .ge. p1st) lptr_pcg8_f_o_aer(02,itype,ai_phase)  = p_pcg8_f_o_a02
        if (p_pcg9_f_o_a02 .ge. p1st) lptr_pcg9_f_o_aer(02,itype,ai_phase)  = p_pcg9_f_o_a02
        if (p_opcg1_f_c_a02 .ge. p1st) lptr_opcg1_f_c_aer(02,itype,ai_phase) = p_opcg1_f_c_a02
        if (p_opcg2_f_c_a02 .ge. p1st) lptr_opcg2_f_c_aer(02,itype,ai_phase)  = p_opcg2_f_c_a02
        if (p_opcg3_f_c_a02 .ge. p1st) lptr_opcg3_f_c_aer(02,itype,ai_phase)  = p_opcg3_f_c_a02
        if (p_opcg4_f_c_a02 .ge. p1st) lptr_opcg4_f_c_aer(02,itype,ai_phase)  = p_opcg4_f_c_a02
        if (p_opcg5_f_c_a02 .ge. p1st) lptr_opcg5_f_c_aer(02,itype,ai_phase)  = p_opcg5_f_c_a02
        if (p_opcg6_f_c_a02 .ge. p1st) lptr_opcg6_f_c_aer(02,itype,ai_phase)  = p_opcg6_f_c_a02
        if (p_opcg7_f_c_a02 .ge. p1st) lptr_opcg7_f_c_aer(02,itype,ai_phase)  = p_opcg7_f_c_a02
        if (p_opcg8_f_c_a02 .ge. p1st) lptr_opcg8_f_c_aer(02,itype,ai_phase)  = p_opcg8_f_c_a02
        if (p_opcg1_f_o_a02 .ge. p1st) lptr_opcg1_f_o_aer(02,itype,ai_phase) = p_opcg1_f_o_a02
        if (p_opcg2_f_o_a02 .ge. p1st) lptr_opcg2_f_o_aer(02,itype,ai_phase)  = p_opcg2_f_o_a02
        if (p_opcg3_f_o_a02 .ge. p1st) lptr_opcg3_f_o_aer(02,itype,ai_phase)  = p_opcg3_f_o_a02
        if (p_opcg4_f_o_a02 .ge. p1st) lptr_opcg4_f_o_aer(02,itype,ai_phase)  = p_opcg4_f_o_a02
        if (p_opcg5_f_o_a02 .ge. p1st) lptr_opcg5_f_o_aer(02,itype,ai_phase)  = p_opcg5_f_o_a02
        if (p_opcg6_f_o_a02 .ge. p1st) lptr_opcg6_f_o_aer(02,itype,ai_phase)  = p_opcg6_f_o_a02
        if (p_opcg7_f_o_a02 .ge. p1st) lptr_opcg7_f_o_aer(02,itype,ai_phase)  = p_opcg7_f_o_a02
        if (p_opcg8_f_o_a02 .ge. p1st) lptr_opcg8_f_o_aer(02,itype,ai_phase)  = p_opcg8_f_o_a02
        if (p_smpa_a02 .ge. p1st) lptr_smpa_aer(02,itype,ai_phase)  = p_smpa_a02
        if (p_smpbb_a02 .ge. p1st) lptr_smpbb_aer(02,itype,ai_phase)  = p_smpbb_a02
        if (p_iepox_a02 .ge. p1st) lptr_iepox_aer(02,itype,ai_phase)          =p_iepox_a02
        if (p_iepoxos_a02 .ge. p1st) lptr_iepoxos_aer(02,itype,ai_phase)      =p_iepoxos_a02
        if (p_tetrol_a02 .ge. p1st) lptr_tetrol_aer(02,itype,ai_phase)        =p_tetrol_a02
        if (p_cldiepoxos_a02 .ge. p1st) lptr_cldiepoxos_aer(02,itype,ai_phase)  =p_cldiepoxos_a02
        if (p_cldtetrol_a02 .ge. p1st) lptr_cldtetrol_aer(02,itype,ai_phase)  =p_cldtetrol_a02
        if (p_gly_a02 .ge. p1st) lptr_gly_aer(02,itype,ai_phase)              =p_gly_a02
        if (p_glysoa_r1_a02 .ge. p1st) lptr_glysoa_r1_aer (02,itype,ai_phase) = p_glysoa_r1_a02
        if (p_glysoa_r2_a02 .ge. p1st) lptr_glysoa_r2_aer (02,itype,ai_phase) = p_glysoa_r2_a02
        if (p_glysoa_sfc_a02 .ge. p1st) lptr_glysoa_sfc_aer(02,itype,ai_phase)  = p_glysoa_sfc_a02
        if (p_glysoa_nh4_a02 .ge. p1st) lptr_glysoa_nh4_aer(02,itype,ai_phase)  = p_glysoa_nh4_a02
        if (p_glysoa_oh_a02 .ge. p1st) lptr_glysoa_oh_aer(02,itype,ai_phase)  = p_glysoa_oh_a02
        if (p_ant1_c_a02 .ge. p1st) lptr_ant1_c_aer(02,itype,ai_phase)  = p_ant1_c_a02
        if (p_ant2_c_a02 .ge. p1st) lptr_ant2_c_aer(02,itype,ai_phase)  = p_ant2_c_a02
        if (p_ant3_c_a02 .ge. p1st) lptr_ant3_c_aer(02,itype,ai_phase)  = p_ant3_c_a02
        if (p_ant4_c_a02 .ge. p1st) lptr_ant4_c_aer(02,itype,ai_phase)  = p_ant4_c_a02
        if (p_ant1_o_a02 .ge. p1st) lptr_ant1_o_aer(02,itype,ai_phase)  = p_ant1_o_a02
        if (p_ant2_o_a02 .ge. p1st) lptr_ant2_o_aer(02,itype,ai_phase)  = p_ant2_o_a02
        if (p_ant3_o_a02 .ge. p1st) lptr_ant3_o_aer(02,itype,ai_phase)  = p_ant3_o_a02
        if (p_ant4_o_a02 .ge. p1st) lptr_ant4_o_aer(02,itype,ai_phase)  = p_ant4_o_a02
        if (p_biog1_c_a02 .ge. p1st) lptr_biog1_c_aer(02,itype,ai_phase)  = p_biog1_c_a02
        if (p_biog2_c_a02 .ge. p1st) lptr_biog2_c_aer(02,itype,ai_phase)  = p_biog2_c_a02
        if (p_biog3_c_a02 .ge. p1st) lptr_biog3_c_aer(02,itype,ai_phase)  = p_biog3_c_a02
        if (p_biog4_c_a02 .ge. p1st) lptr_biog4_c_aer(02,itype,ai_phase)  = p_biog4_c_a02
        if (p_biog1_o_a02 .ge. p1st) lptr_biog1_o_aer(02,itype,ai_phase)  = p_biog1_o_a02
        if (p_biog2_o_a02 .ge. p1st) lptr_biog2_o_aer(02,itype,ai_phase)  = p_biog2_o_a02
        if (p_biog3_o_a02 .ge. p1st) lptr_biog3_o_aer(02,itype,ai_phase)  = p_biog3_o_a02
        if (p_biog4_o_a02 .ge. p1st) lptr_biog4_o_aer(02,itype,ai_phase)  = p_biog4_o_a02
        if (p_asoaX_a02 .ge. p1st) lptr_asoaX_aer(02,itype,ai_phase)  = p_asoaX_a02
        if (p_asoa1_a02 .ge. p1st) lptr_asoa1_aer(02,itype,ai_phase)  = p_asoa1_a02
        if (p_asoa2_a02 .ge. p1st) lptr_asoa2_aer(02,itype,ai_phase)  = p_asoa2_a02
        if (p_asoa3_a02 .ge. p1st) lptr_asoa3_aer(02,itype,ai_phase)  = p_asoa3_a02
        if (p_asoa4_a02 .ge. p1st) lptr_asoa4_aer(02,itype,ai_phase)  = p_asoa4_a02
        if (p_bsoaX_a02 .ge. p1st) lptr_bsoaX_aer(02,itype,ai_phase)  = p_bsoaX_a02
        if (p_bsoa1_a02 .ge. p1st) lptr_bsoa1_aer(02,itype,ai_phase)  = p_bsoa1_a02
        if (p_bsoa2_a02 .ge. p1st) lptr_bsoa2_aer(02,itype,ai_phase)  = p_bsoa2_a02
        if (p_bsoa3_a02 .ge. p1st) lptr_bsoa3_aer(02,itype,ai_phase)  = p_bsoa3_a02
        if (p_bsoa4_a02 .ge. p1st) lptr_bsoa4_aer(02,itype,ai_phase)  = p_bsoa4_a02
        if (p_num_a02 .ge. p1st)  numptr_aer(02,itype,ai_phase)        = p_num_a02
        
        if (p_TERP_CN3_a02 .ge. p1st) lptr_terp_Cn3_aer(02,itype,ai_phase)  = p_TERP_CN3_a02
        end if

	if (nsize_aer(itype) .ge. 3) then
        if (p_so4_a03 .ge. p1st) lptr_so4_aer(03,itype,ai_phase) = p_so4_a03
        if (p_no3_a03 .ge. p1st) lptr_no3_aer(03,itype,ai_phase) = p_no3_a03
        if (p_cl_a03 .ge. p1st) lptr_cl_aer(03,itype,ai_phase)  = p_cl_a03
        if (p_msa_a03 .ge. p1st) lptr_msa_aer(03,itype,ai_phase) = p_msa_a03
        if (p_co3_a03 .ge. p1st) lptr_co3_aer(03,itype,ai_phase) = p_co3_a03
        if (p_nh4_a03 .ge. p1st) lptr_nh4_aer(03,itype,ai_phase) = p_nh4_a03
        if (p_na_a03 .ge. p1st) lptr_na_aer(03,itype,ai_phase)  = p_na_a03
        if (p_ca_a03 .ge. p1st) lptr_ca_aer(03,itype,ai_phase)  = p_ca_a03
        if (p_oin_a03 .ge. p1st) lptr_oin_aer(03,itype,ai_phase) = p_oin_a03
        if (p_oc_a03 .ge. p1st) lptr_oc_aer(03,itype,ai_phase)  = p_oc_a03
        if (p_bc_a03 .ge. p1st) lptr_bc_aer(03,itype,ai_phase)  = p_bc_a03
        if (p_hysw_a03 .ge. p1st) hyswptr_aer(03,itype)         = p_hysw_a03
        if (p_water_a03 .ge. p1st) waterptr_aer(03,itype)       = p_water_a03
        if (p_tanv_a03 .ge. p1st) lptr_tanv_aer(03,itype,ai_phase)  = p_tanv_a03
        if (p_sopnv_a03 .ge. p1st) lptr_sopnv_aer(03,itype,ai_phase)  = p_sopnv_a03
        if (p_ternv_a03 .ge. p1st) lptr_ternv_aer(03,itype,ai_phase)  = p_ternv_a03
        if (p_seqnv_a03 .ge. p1st) lptr_seqnv_aer(03,itype,ai_phase)  = p_seqnv_a03
        if (p_sianv_a03 .ge. p1st) lptr_sianv_aer(03,itype,ai_phase)  = p_sianv_a03
        if (p_pcg1_b_c_a03 .ge. p1st) lptr_pcg1_b_c_aer(03,itype,ai_phase) = p_pcg1_b_c_a03
        if (p_pcg2_b_c_a03 .ge. p1st) lptr_pcg2_b_c_aer(03,itype,ai_phase)  = p_pcg2_b_c_a03
        if (p_pcg3_b_c_a03 .ge. p1st) lptr_pcg3_b_c_aer(03,itype,ai_phase)  = p_pcg3_b_c_a03
        if (p_pcg4_b_c_a03 .ge. p1st) lptr_pcg4_b_c_aer(03,itype,ai_phase)  = p_pcg4_b_c_a03
        if (p_pcg5_b_c_a03 .ge. p1st) lptr_pcg5_b_c_aer(03,itype,ai_phase)  = p_pcg5_b_c_a03
        if (p_pcg6_b_c_a03 .ge. p1st) lptr_pcg6_b_c_aer(03,itype,ai_phase)  = p_pcg6_b_c_a03
        if (p_pcg7_b_c_a03 .ge. p1st) lptr_pcg7_b_c_aer(03,itype,ai_phase)  = p_pcg7_b_c_a03
        if (p_pcg8_b_c_a03 .ge. p1st) lptr_pcg8_b_c_aer(03,itype,ai_phase)  = p_pcg8_b_c_a03
        if (p_pcg9_b_c_a03 .ge. p1st) lptr_pcg9_b_c_aer(03,itype,ai_phase)  = p_pcg9_b_c_a03
        if (p_pcg1_b_o_a03 .ge. p1st) lptr_pcg1_b_o_aer(03,itype,ai_phase) = p_pcg1_b_o_a03
        if (p_pcg2_b_o_a03 .ge. p1st) lptr_pcg2_b_o_aer(03,itype,ai_phase)  = p_pcg2_b_o_a03
        if (p_pcg3_b_o_a03 .ge. p1st) lptr_pcg3_b_o_aer(03,itype,ai_phase)  = p_pcg3_b_o_a03
        if (p_pcg4_b_o_a03 .ge. p1st) lptr_pcg4_b_o_aer(03,itype,ai_phase)  = p_pcg4_b_o_a03
        if (p_pcg5_b_o_a03 .ge. p1st) lptr_pcg5_b_o_aer(03,itype,ai_phase)  = p_pcg5_b_o_a03
        if (p_pcg6_b_o_a03 .ge. p1st) lptr_pcg6_b_o_aer(03,itype,ai_phase)  = p_pcg6_b_o_a03
        if (p_pcg7_b_o_a03 .ge. p1st) lptr_pcg7_b_o_aer(03,itype,ai_phase)  = p_pcg7_b_o_a03
        if (p_pcg8_b_o_a03 .ge. p1st) lptr_pcg8_b_o_aer(03,itype,ai_phase)  = p_pcg8_b_o_a03
        if (p_pcg9_b_o_a03 .ge. p1st) lptr_pcg9_b_o_aer(03,itype,ai_phase)  = p_pcg9_b_o_a03
        if (p_opcg1_b_c_a03 .ge. p1st) lptr_opcg1_b_c_aer(03,itype,ai_phase) = p_opcg1_b_c_a03
        if (p_opcg2_b_c_a03 .ge. p1st) lptr_opcg2_b_c_aer(03,itype,ai_phase)  = p_opcg2_b_c_a03
        if (p_opcg3_b_c_a03 .ge. p1st) lptr_opcg3_b_c_aer(03,itype,ai_phase)  = p_opcg3_b_c_a03
        if (p_opcg4_b_c_a03 .ge. p1st) lptr_opcg4_b_c_aer(03,itype,ai_phase)  = p_opcg4_b_c_a03
        if (p_opcg5_b_c_a03 .ge. p1st) lptr_opcg5_b_c_aer(03,itype,ai_phase)  = p_opcg5_b_c_a03
        if (p_opcg6_b_c_a03 .ge. p1st) lptr_opcg6_b_c_aer(03,itype,ai_phase)  = p_opcg6_b_c_a03
        if (p_opcg7_b_c_a03 .ge. p1st) lptr_opcg7_b_c_aer(03,itype,ai_phase)  = p_opcg7_b_c_a03
        if (p_opcg8_b_c_a03 .ge. p1st) lptr_opcg8_b_c_aer(03,itype,ai_phase)  = p_opcg8_b_c_a03
        if (p_opcg1_b_o_a03 .ge. p1st) lptr_opcg1_b_o_aer(03,itype,ai_phase) = p_opcg1_b_o_a03
        if (p_opcg2_b_o_a03 .ge. p1st) lptr_opcg2_b_o_aer(03,itype,ai_phase)  = p_opcg2_b_o_a03
        if (p_opcg3_b_o_a03 .ge. p1st) lptr_opcg3_b_o_aer(03,itype,ai_phase)  = p_opcg3_b_o_a03
        if (p_opcg4_b_o_a03 .ge. p1st) lptr_opcg4_b_o_aer(03,itype,ai_phase)  = p_opcg4_b_o_a03
        if (p_opcg5_b_o_a03 .ge. p1st) lptr_opcg5_b_o_aer(03,itype,ai_phase)  = p_opcg5_b_o_a03
        if (p_opcg6_b_o_a03 .ge. p1st) lptr_opcg6_b_o_aer(03,itype,ai_phase)  = p_opcg6_b_o_a03
        if (p_opcg7_b_o_a03 .ge. p1st) lptr_opcg7_b_o_aer(03,itype,ai_phase)  = p_opcg7_b_o_a03
        if (p_opcg8_b_o_a03 .ge. p1st) lptr_opcg8_b_o_aer(03,itype,ai_phase)  = p_opcg8_b_o_a03
        if (p_pcg1_f_c_a03 .ge. p1st) lptr_pcg1_f_c_aer(03,itype,ai_phase) = p_pcg1_f_c_a03
        if (p_pcg2_f_c_a03 .ge. p1st) lptr_pcg2_f_c_aer(03,itype,ai_phase)  = p_pcg2_f_c_a03
        if (p_pcg3_f_c_a03 .ge. p1st) lptr_pcg3_f_c_aer(03,itype,ai_phase)  = p_pcg3_f_c_a03
        if (p_pcg4_f_c_a03 .ge. p1st) lptr_pcg4_f_c_aer(03,itype,ai_phase)  = p_pcg4_f_c_a03
        if (p_pcg5_f_c_a03 .ge. p1st) lptr_pcg5_f_c_aer(03,itype,ai_phase)  = p_pcg5_f_c_a03
        if (p_pcg6_f_c_a03 .ge. p1st) lptr_pcg6_f_c_aer(03,itype,ai_phase)  = p_pcg6_f_c_a03
        if (p_pcg7_f_c_a03 .ge. p1st) lptr_pcg7_f_c_aer(03,itype,ai_phase)  = p_pcg7_f_c_a03
        if (p_pcg8_f_c_a03 .ge. p1st) lptr_pcg8_f_c_aer(03,itype,ai_phase)  = p_pcg8_f_c_a03
        if (p_pcg9_f_c_a03 .ge. p1st) lptr_pcg9_f_c_aer(03,itype,ai_phase)  = p_pcg9_f_c_a03
        if (p_pcg1_f_o_a03 .ge. p1st) lptr_pcg1_f_o_aer(03,itype,ai_phase) = p_pcg1_f_o_a03
        if (p_pcg2_f_o_a03 .ge. p1st) lptr_pcg2_f_o_aer(03,itype,ai_phase)  = p_pcg2_f_o_a03
        if (p_pcg3_f_o_a03 .ge. p1st) lptr_pcg3_f_o_aer(03,itype,ai_phase)  = p_pcg3_f_o_a03
        if (p_pcg4_f_o_a03 .ge. p1st) lptr_pcg4_f_o_aer(03,itype,ai_phase)  = p_pcg4_f_o_a03
        if (p_pcg5_f_o_a03 .ge. p1st) lptr_pcg5_f_o_aer(03,itype,ai_phase)  = p_pcg5_f_o_a03
        if (p_pcg6_f_o_a03 .ge. p1st) lptr_pcg6_f_o_aer(03,itype,ai_phase)  = p_pcg6_f_o_a03
        if (p_pcg7_f_o_a03 .ge. p1st) lptr_pcg7_f_o_aer(03,itype,ai_phase)  = p_pcg7_f_o_a03
        if (p_pcg8_f_o_a03 .ge. p1st) lptr_pcg8_f_o_aer(03,itype,ai_phase)  = p_pcg8_f_o_a03
        if (p_pcg9_f_o_a03 .ge. p1st) lptr_pcg9_f_o_aer(03,itype,ai_phase)  = p_pcg9_f_o_a03
        if (p_opcg1_f_c_a03 .ge. p1st) lptr_opcg1_f_c_aer(03,itype,ai_phase) = p_opcg1_f_c_a03
        if (p_opcg2_f_c_a03 .ge. p1st) lptr_opcg2_f_c_aer(03,itype,ai_phase)  = p_opcg2_f_c_a03
        if (p_opcg3_f_c_a03 .ge. p1st) lptr_opcg3_f_c_aer(03,itype,ai_phase)  = p_opcg3_f_c_a03
        if (p_opcg4_f_c_a03 .ge. p1st) lptr_opcg4_f_c_aer(03,itype,ai_phase)  = p_opcg4_f_c_a03
        if (p_opcg5_f_c_a03 .ge. p1st) lptr_opcg5_f_c_aer(03,itype,ai_phase)  = p_opcg5_f_c_a03
        if (p_opcg6_f_c_a03 .ge. p1st) lptr_opcg6_f_c_aer(03,itype,ai_phase)  = p_opcg6_f_c_a03
        if (p_opcg7_f_c_a03 .ge. p1st) lptr_opcg7_f_c_aer(03,itype,ai_phase)  = p_opcg7_f_c_a03
        if (p_opcg8_f_c_a03 .ge. p1st) lptr_opcg8_f_c_aer(03,itype,ai_phase)  = p_opcg8_f_c_a03
        if (p_opcg1_f_o_a03 .ge. p1st) lptr_opcg1_f_o_aer(03,itype,ai_phase) = p_opcg1_f_o_a03
        if (p_opcg2_f_o_a03 .ge. p1st) lptr_opcg2_f_o_aer(03,itype,ai_phase)  = p_opcg2_f_o_a03
        if (p_opcg3_f_o_a03 .ge. p1st) lptr_opcg3_f_o_aer(03,itype,ai_phase)  = p_opcg3_f_o_a03
        if (p_opcg4_f_o_a03 .ge. p1st) lptr_opcg4_f_o_aer(03,itype,ai_phase)  = p_opcg4_f_o_a03
        if (p_opcg5_f_o_a03 .ge. p1st) lptr_opcg5_f_o_aer(03,itype,ai_phase)  = p_opcg5_f_o_a03
        if (p_opcg6_f_o_a03 .ge. p1st) lptr_opcg6_f_o_aer(03,itype,ai_phase)  = p_opcg6_f_o_a03
        if (p_opcg7_f_o_a03 .ge. p1st) lptr_opcg7_f_o_aer(03,itype,ai_phase)  = p_opcg7_f_o_a03
        if (p_opcg8_f_o_a03 .ge. p1st) lptr_opcg8_f_o_aer(03,itype,ai_phase)  = p_opcg8_f_o_a03
        if (p_smpa_a03 .ge. p1st) lptr_smpa_aer(03,itype,ai_phase)  = p_smpa_a03
        if (p_smpbb_a03 .ge. p1st) lptr_smpbb_aer(03,itype,ai_phase)  = p_smpbb_a03
        if (p_iepox_a03 .ge. p1st) lptr_iepox_aer(03,itype,ai_phase)          =p_iepox_a03
        if (p_iepoxos_a03 .ge. p1st) lptr_iepoxos_aer(03,itype,ai_phase)      =p_iepoxos_a03
        if (p_tetrol_a03 .ge. p1st) lptr_tetrol_aer(03,itype,ai_phase)        =p_tetrol_a03
        if (p_cldiepoxos_a03 .ge. p1st) lptr_cldiepoxos_aer(03,itype,ai_phase)  =p_cldiepoxos_a03
        if (p_cldtetrol_a03 .ge. p1st) lptr_cldtetrol_aer(03,itype,ai_phase)  =p_cldtetrol_a03
        if (p_gly_a03 .ge. p1st) lptr_gly_aer(03,itype,ai_phase)              =p_gly_a03
        if (p_glysoa_r1_a03 .ge. p1st) lptr_glysoa_r1_aer (03,itype,ai_phase) = p_glysoa_r1_a03
        if (p_glysoa_r2_a03 .ge. p1st) lptr_glysoa_r2_aer (03,itype,ai_phase) = p_glysoa_r2_a03
        if (p_glysoa_sfc_a03 .ge. p1st) lptr_glysoa_sfc_aer(03,itype,ai_phase)  = p_glysoa_sfc_a03
        if (p_glysoa_nh4_a03 .ge. p1st) lptr_glysoa_nh4_aer(03,itype,ai_phase)  = p_glysoa_nh4_a03
        if (p_glysoa_oh_a03 .ge. p1st) lptr_glysoa_oh_aer(03,itype,ai_phase)  = p_glysoa_oh_a03
        if (p_ant1_c_a03 .ge. p1st) lptr_ant1_c_aer(03,itype,ai_phase)  = p_ant1_c_a03
        if (p_ant2_c_a03 .ge. p1st) lptr_ant2_c_aer(03,itype,ai_phase)  = p_ant2_c_a03
        if (p_ant3_c_a03 .ge. p1st) lptr_ant3_c_aer(03,itype,ai_phase)  = p_ant3_c_a03
        if (p_ant4_c_a03 .ge. p1st) lptr_ant4_c_aer(03,itype,ai_phase)  = p_ant4_c_a03
        if (p_ant1_o_a03 .ge. p1st) lptr_ant1_o_aer(03,itype,ai_phase)  = p_ant1_o_a03
        if (p_ant2_o_a03 .ge. p1st) lptr_ant2_o_aer(03,itype,ai_phase)  = p_ant2_o_a03
        if (p_ant3_o_a03 .ge. p1st) lptr_ant3_o_aer(03,itype,ai_phase)  = p_ant3_o_a03
        if (p_ant4_o_a03 .ge. p1st) lptr_ant4_o_aer(03,itype,ai_phase)  = p_ant4_o_a03
        if (p_biog1_c_a03 .ge. p1st) lptr_biog1_c_aer(03,itype,ai_phase)  = p_biog1_c_a03
        if (p_biog2_c_a03 .ge. p1st) lptr_biog2_c_aer(03,itype,ai_phase)  = p_biog2_c_a03
        if (p_biog3_c_a03 .ge. p1st) lptr_biog3_c_aer(03,itype,ai_phase)  = p_biog3_c_a03
        if (p_biog4_c_a03 .ge. p1st) lptr_biog4_c_aer(03,itype,ai_phase)  = p_biog4_c_a03
        if (p_biog1_o_a03 .ge. p1st) lptr_biog1_o_aer(03,itype,ai_phase)  = p_biog1_o_a03
        if (p_biog2_o_a03 .ge. p1st) lptr_biog2_o_aer(03,itype,ai_phase)  = p_biog2_o_a03
        if (p_biog3_o_a03 .ge. p1st) lptr_biog3_o_aer(03,itype,ai_phase)  = p_biog3_o_a03
        if (p_biog4_o_a03 .ge. p1st) lptr_biog4_o_aer(03,itype,ai_phase)  = p_biog4_o_a03
        if (p_asoaX_a03 .ge. p1st) lptr_asoaX_aer(03,itype,ai_phase)  = p_asoaX_a03
        if (p_asoa1_a03 .ge. p1st) lptr_asoa1_aer(03,itype,ai_phase)  = p_asoa1_a03
        if (p_asoa2_a03 .ge. p1st) lptr_asoa2_aer(03,itype,ai_phase)  = p_asoa2_a03
        if (p_asoa3_a03 .ge. p1st) lptr_asoa3_aer(03,itype,ai_phase)  = p_asoa3_a03
        if (p_asoa4_a03 .ge. p1st) lptr_asoa4_aer(03,itype,ai_phase)  = p_asoa4_a03
        if (p_bsoaX_a03 .ge. p1st) lptr_bsoaX_aer(03,itype,ai_phase)  = p_bsoaX_a03
        if (p_bsoa1_a03 .ge. p1st) lptr_bsoa1_aer(03,itype,ai_phase)  = p_bsoa1_a03
        if (p_bsoa2_a03 .ge. p1st) lptr_bsoa2_aer(03,itype,ai_phase)  = p_bsoa2_a03
        if (p_bsoa3_a03 .ge. p1st) lptr_bsoa3_aer(03,itype,ai_phase)  = p_bsoa3_a03
        if (p_bsoa4_a03 .ge. p1st) lptr_bsoa4_aer(03,itype,ai_phase)  = p_bsoa4_a03
        if (p_num_a03 .ge. p1st)  numptr_aer(03,itype,ai_phase)        = p_num_a03
        
        if (p_TERP_CN3_a03 .ge. p1st) lptr_terp_Cn3_aer(03,itype,ai_phase)  = p_TERP_CN3_a03
        end if


	if (nsize_aer(itype) .ge. 4) then
	if (p_so4_a04 .ge. p1st) lptr_so4_aer(04,itype,ai_phase) = p_so4_a04
        if (p_no3_a04 .ge. p1st) lptr_no3_aer(04,itype,ai_phase) = p_no3_a04
        if (p_cl_a04 .ge. p1st) lptr_cl_aer(04,itype,ai_phase)  = p_cl_a04
        if (p_msa_a04 .ge. p1st) lptr_msa_aer(04,itype,ai_phase) = p_msa_a04
        if (p_co3_a04 .ge. p1st) lptr_co3_aer(04,itype,ai_phase) = p_co3_a04
        if (p_nh4_a04 .ge. p1st) lptr_nh4_aer(04,itype,ai_phase) = p_nh4_a04
        if (p_na_a04 .ge. p1st) lptr_na_aer(04,itype,ai_phase)  = p_na_a04
        if (p_ca_a04 .ge. p1st) lptr_ca_aer(04,itype,ai_phase)  = p_ca_a04
        if (p_oin_a04 .ge. p1st) lptr_oin_aer(04,itype,ai_phase) = p_oin_a04
        if (p_oc_a04 .ge. p1st) lptr_oc_aer(04,itype,ai_phase)  = p_oc_a04
        if (p_bc_a04 .ge. p1st) lptr_bc_aer(04,itype,ai_phase)  = p_bc_a04
        if (p_hysw_a04 .ge. p1st) hyswptr_aer(04,itype)         = p_hysw_a04
        if (p_water_a04 .ge. p1st) waterptr_aer(04,itype)       = p_water_a04
        if (p_tanv_a04 .ge. p1st) lptr_tanv_aer(04,itype,ai_phase)  = p_tanv_a04
        if (p_sopnv_a04 .ge. p1st) lptr_sopnv_aer(04,itype,ai_phase)  = p_sopnv_a04
        if (p_ternv_a04 .ge. p1st) lptr_ternv_aer(04,itype,ai_phase)  = p_ternv_a04
        if (p_seqnv_a04 .ge. p1st) lptr_seqnv_aer(04,itype,ai_phase)  = p_seqnv_a04
        if (p_sianv_a04 .ge. p1st) lptr_sianv_aer(04,itype,ai_phase)  = p_sianv_a04
        if (p_pcg1_b_c_a04 .ge. p1st) lptr_pcg1_b_c_aer(04,itype,ai_phase) = p_pcg1_b_c_a04
        if (p_pcg2_b_c_a04 .ge. p1st) lptr_pcg2_b_c_aer(04,itype,ai_phase)  = p_pcg2_b_c_a04
        if (p_pcg3_b_c_a04 .ge. p1st) lptr_pcg3_b_c_aer(04,itype,ai_phase)  = p_pcg3_b_c_a04
        if (p_pcg4_b_c_a04 .ge. p1st) lptr_pcg4_b_c_aer(04,itype,ai_phase)  = p_pcg4_b_c_a04
        if (p_pcg5_b_c_a04 .ge. p1st) lptr_pcg5_b_c_aer(04,itype,ai_phase)  = p_pcg5_b_c_a04
        if (p_pcg6_b_c_a04 .ge. p1st) lptr_pcg6_b_c_aer(04,itype,ai_phase)  = p_pcg6_b_c_a04
        if (p_pcg7_b_c_a04 .ge. p1st) lptr_pcg7_b_c_aer(04,itype,ai_phase)  = p_pcg7_b_c_a04
        if (p_pcg8_b_c_a04 .ge. p1st) lptr_pcg8_b_c_aer(04,itype,ai_phase)  = p_pcg8_b_c_a04
        if (p_pcg9_b_c_a04 .ge. p1st) lptr_pcg9_b_c_aer(04,itype,ai_phase)  = p_pcg9_b_c_a04
        if (p_pcg1_b_o_a04 .ge. p1st) lptr_pcg1_b_o_aer(04,itype,ai_phase) = p_pcg1_b_o_a04
        if (p_pcg2_b_o_a04 .ge. p1st) lptr_pcg2_b_o_aer(04,itype,ai_phase)  = p_pcg2_b_o_a04
        if (p_pcg3_b_o_a04 .ge. p1st) lptr_pcg3_b_o_aer(04,itype,ai_phase)  = p_pcg3_b_o_a04
        if (p_pcg4_b_o_a04 .ge. p1st) lptr_pcg4_b_o_aer(04,itype,ai_phase)  = p_pcg4_b_o_a04
        if (p_pcg5_b_o_a04 .ge. p1st) lptr_pcg5_b_o_aer(04,itype,ai_phase)  = p_pcg5_b_o_a04
        if (p_pcg6_b_o_a04 .ge. p1st) lptr_pcg6_b_o_aer(04,itype,ai_phase)  = p_pcg6_b_o_a04
        if (p_pcg7_b_o_a04 .ge. p1st) lptr_pcg7_b_o_aer(04,itype,ai_phase)  = p_pcg7_b_o_a04
        if (p_pcg8_b_o_a04 .ge. p1st) lptr_pcg8_b_o_aer(04,itype,ai_phase)  = p_pcg8_b_o_a04
        if (p_pcg9_b_o_a04 .ge. p1st) lptr_pcg9_b_o_aer(04,itype,ai_phase)  = p_pcg9_b_o_a04
        if (p_opcg1_b_c_a04 .ge. p1st) lptr_opcg1_b_c_aer(04,itype,ai_phase) = p_opcg1_b_c_a04
        if (p_opcg2_b_c_a04 .ge. p1st) lptr_opcg2_b_c_aer(04,itype,ai_phase)  = p_opcg2_b_c_a04
        if (p_opcg3_b_c_a04 .ge. p1st) lptr_opcg3_b_c_aer(04,itype,ai_phase)  = p_opcg3_b_c_a04
        if (p_opcg4_b_c_a04 .ge. p1st) lptr_opcg4_b_c_aer(04,itype,ai_phase)  = p_opcg4_b_c_a04
        if (p_opcg5_b_c_a04 .ge. p1st) lptr_opcg5_b_c_aer(04,itype,ai_phase)  = p_opcg5_b_c_a04
        if (p_opcg6_b_c_a04 .ge. p1st) lptr_opcg6_b_c_aer(04,itype,ai_phase)  = p_opcg6_b_c_a04
        if (p_opcg7_b_c_a04 .ge. p1st) lptr_opcg7_b_c_aer(04,itype,ai_phase)  = p_opcg7_b_c_a04
        if (p_opcg8_b_c_a04 .ge. p1st) lptr_opcg8_b_c_aer(04,itype,ai_phase)  = p_opcg8_b_c_a04
        if (p_opcg1_b_o_a04 .ge. p1st) lptr_opcg1_b_o_aer(04,itype,ai_phase) = p_opcg1_b_o_a04
        if (p_opcg2_b_o_a04 .ge. p1st) lptr_opcg2_b_o_aer(04,itype,ai_phase)  = p_opcg2_b_o_a04
        if (p_opcg3_b_o_a04 .ge. p1st) lptr_opcg3_b_o_aer(04,itype,ai_phase)  = p_opcg3_b_o_a04
        if (p_opcg4_b_o_a04 .ge. p1st) lptr_opcg4_b_o_aer(04,itype,ai_phase)  = p_opcg4_b_o_a04
        if (p_opcg5_b_o_a04 .ge. p1st) lptr_opcg5_b_o_aer(04,itype,ai_phase)  = p_opcg5_b_o_a04
        if (p_opcg6_b_o_a04 .ge. p1st) lptr_opcg6_b_o_aer(04,itype,ai_phase)  = p_opcg6_b_o_a04
        if (p_opcg7_b_o_a04 .ge. p1st) lptr_opcg7_b_o_aer(04,itype,ai_phase)  = p_opcg7_b_o_a04
        if (p_opcg8_b_o_a04 .ge. p1st) lptr_opcg8_b_o_aer(04,itype,ai_phase)  = p_opcg8_b_o_a04
        if (p_pcg1_f_c_a04 .ge. p1st) lptr_pcg1_f_c_aer(04,itype,ai_phase) = p_pcg1_f_c_a04
        if (p_pcg2_f_c_a04 .ge. p1st) lptr_pcg2_f_c_aer(04,itype,ai_phase)  = p_pcg2_f_c_a04
        if (p_pcg3_f_c_a04 .ge. p1st) lptr_pcg3_f_c_aer(04,itype,ai_phase)  = p_pcg3_f_c_a04
        if (p_pcg4_f_c_a04 .ge. p1st) lptr_pcg4_f_c_aer(04,itype,ai_phase)  = p_pcg4_f_c_a04
        if (p_pcg5_f_c_a04 .ge. p1st) lptr_pcg5_f_c_aer(04,itype,ai_phase)  = p_pcg5_f_c_a04
        if (p_pcg6_f_c_a04 .ge. p1st) lptr_pcg6_f_c_aer(04,itype,ai_phase)  = p_pcg6_f_c_a04
        if (p_pcg7_f_c_a04 .ge. p1st) lptr_pcg7_f_c_aer(04,itype,ai_phase)  = p_pcg7_f_c_a04
        if (p_pcg8_f_c_a04 .ge. p1st) lptr_pcg8_f_c_aer(04,itype,ai_phase)  = p_pcg8_f_c_a04
        if (p_pcg9_f_c_a04 .ge. p1st) lptr_pcg9_f_c_aer(04,itype,ai_phase)  = p_pcg9_f_c_a04
        if (p_pcg1_f_o_a04 .ge. p1st) lptr_pcg1_f_o_aer(04,itype,ai_phase) = p_pcg1_f_o_a04
        if (p_pcg2_f_o_a04 .ge. p1st) lptr_pcg2_f_o_aer(04,itype,ai_phase)  = p_pcg2_f_o_a04
        if (p_pcg3_f_o_a04 .ge. p1st) lptr_pcg3_f_o_aer(04,itype,ai_phase)  = p_pcg3_f_o_a04
        if (p_pcg4_f_o_a04 .ge. p1st) lptr_pcg4_f_o_aer(04,itype,ai_phase)  = p_pcg4_f_o_a04
        if (p_pcg5_f_o_a04 .ge. p1st) lptr_pcg5_f_o_aer(04,itype,ai_phase)  = p_pcg5_f_o_a04
        if (p_pcg6_f_o_a04 .ge. p1st) lptr_pcg6_f_o_aer(04,itype,ai_phase)  = p_pcg6_f_o_a04
        if (p_pcg7_f_o_a04 .ge. p1st) lptr_pcg7_f_o_aer(04,itype,ai_phase)  = p_pcg7_f_o_a04
        if (p_pcg8_f_o_a04 .ge. p1st) lptr_pcg8_f_o_aer(04,itype,ai_phase)  = p_pcg8_f_o_a04
        if (p_pcg9_f_o_a04 .ge. p1st) lptr_pcg9_f_o_aer(04,itype,ai_phase)  = p_pcg9_f_o_a04
        if (p_opcg1_f_c_a04 .ge. p1st) lptr_opcg1_f_c_aer(04,itype,ai_phase) = p_opcg1_f_c_a04
        if (p_opcg2_f_c_a04 .ge. p1st) lptr_opcg2_f_c_aer(04,itype,ai_phase)  = p_opcg2_f_c_a04
        if (p_opcg3_f_c_a04 .ge. p1st) lptr_opcg3_f_c_aer(04,itype,ai_phase)  = p_opcg3_f_c_a04
        if (p_opcg4_f_c_a04 .ge. p1st) lptr_opcg4_f_c_aer(04,itype,ai_phase)  = p_opcg4_f_c_a04
        if (p_opcg5_f_c_a04 .ge. p1st) lptr_opcg5_f_c_aer(04,itype,ai_phase)  = p_opcg5_f_c_a04
        if (p_opcg6_f_c_a04 .ge. p1st) lptr_opcg6_f_c_aer(04,itype,ai_phase)  = p_opcg6_f_c_a04
        if (p_opcg7_f_c_a04 .ge. p1st) lptr_opcg7_f_c_aer(04,itype,ai_phase)  = p_opcg7_f_c_a04
        if (p_opcg8_f_c_a04 .ge. p1st) lptr_opcg8_f_c_aer(04,itype,ai_phase)  = p_opcg8_f_c_a04
        if (p_opcg1_f_o_a04 .ge. p1st) lptr_opcg1_f_o_aer(04,itype,ai_phase) = p_opcg1_f_o_a04
        if (p_opcg2_f_o_a04 .ge. p1st) lptr_opcg2_f_o_aer(04,itype,ai_phase)  = p_opcg2_f_o_a04
        if (p_opcg3_f_o_a04 .ge. p1st) lptr_opcg3_f_o_aer(04,itype,ai_phase)  = p_opcg3_f_o_a04
        if (p_opcg4_f_o_a04 .ge. p1st) lptr_opcg4_f_o_aer(04,itype,ai_phase)  = p_opcg4_f_o_a04
        if (p_opcg5_f_o_a04 .ge. p1st) lptr_opcg5_f_o_aer(04,itype,ai_phase)  = p_opcg5_f_o_a04
        if (p_opcg6_f_o_a04 .ge. p1st) lptr_opcg6_f_o_aer(04,itype,ai_phase)  = p_opcg6_f_o_a04
        if (p_opcg7_f_o_a04 .ge. p1st) lptr_opcg7_f_o_aer(04,itype,ai_phase)  = p_opcg7_f_o_a04
        if (p_opcg8_f_o_a04 .ge. p1st) lptr_opcg8_f_o_aer(04,itype,ai_phase)  = p_opcg8_f_o_a04
        if (p_smpa_a04 .ge. p1st) lptr_smpa_aer(04,itype,ai_phase)  = p_smpa_a04
        if (p_smpbb_a04 .ge. p1st) lptr_smpbb_aer(04,itype,ai_phase)  = p_smpbb_a04
        if (p_iepox_a04 .ge. p1st) lptr_iepox_aer(04,itype,ai_phase)          =p_iepox_a04
        if (p_iepoxos_a04 .ge. p1st) lptr_iepoxos_aer(04,itype,ai_phase)      =p_iepoxos_a04
        if (p_tetrol_a04 .ge. p1st) lptr_tetrol_aer(04,itype,ai_phase)        =p_tetrol_a04
        if (p_cldiepoxos_a04 .ge. p1st) lptr_cldiepoxos_aer(04,itype,ai_phase)  =p_cldiepoxos_a04
        if (p_cldtetrol_a04 .ge. p1st) lptr_cldtetrol_aer(04,itype,ai_phase)  =p_cldtetrol_a04
        if (p_gly_a04 .ge. p1st) lptr_gly_aer(04,itype,ai_phase)              =p_gly_a04
        if (p_glysoa_r1_a04 .ge. p1st) lptr_glysoa_r1_aer (04,itype,ai_phase) = p_glysoa_r1_a04
        if (p_glysoa_r2_a04 .ge. p1st) lptr_glysoa_r2_aer (04,itype,ai_phase) = p_glysoa_r2_a04
        if (p_glysoa_sfc_a04 .ge. p1st) lptr_glysoa_sfc_aer(04,itype,ai_phase)  = p_glysoa_sfc_a04
        if (p_glysoa_nh4_a04 .ge. p1st) lptr_glysoa_nh4_aer(04,itype,ai_phase)  = p_glysoa_nh4_a04
        if (p_glysoa_oh_a04 .ge. p1st) lptr_glysoa_oh_aer(04,itype,ai_phase)  = p_glysoa_oh_a04
        if (p_ant1_c_a04 .ge. p1st) lptr_ant1_c_aer(04,itype,ai_phase)  = p_ant1_c_a04
        if (p_ant2_c_a04 .ge. p1st) lptr_ant2_c_aer(04,itype,ai_phase)  = p_ant2_c_a04
        if (p_ant3_c_a04 .ge. p1st) lptr_ant3_c_aer(04,itype,ai_phase)  = p_ant3_c_a04
        if (p_ant4_c_a04 .ge. p1st) lptr_ant4_c_aer(04,itype,ai_phase)  = p_ant4_c_a04
        if (p_ant1_o_a04 .ge. p1st) lptr_ant1_o_aer(04,itype,ai_phase)  = p_ant1_o_a04
        if (p_ant2_o_a04 .ge. p1st) lptr_ant2_o_aer(04,itype,ai_phase)  = p_ant2_o_a04
        if (p_ant3_o_a04 .ge. p1st) lptr_ant3_o_aer(04,itype,ai_phase)  = p_ant3_o_a04
        if (p_ant4_o_a04 .ge. p1st) lptr_ant4_o_aer(04,itype,ai_phase)  = p_ant4_o_a04
        if (p_biog1_c_a04 .ge. p1st) lptr_biog1_c_aer(04,itype,ai_phase)  = p_biog1_c_a04
        if (p_biog2_c_a04 .ge. p1st) lptr_biog2_c_aer(04,itype,ai_phase)  = p_biog2_c_a04
        if (p_biog3_c_a04 .ge. p1st) lptr_biog3_c_aer(04,itype,ai_phase)  = p_biog3_c_a04
        if (p_biog4_c_a04 .ge. p1st) lptr_biog4_c_aer(04,itype,ai_phase)  = p_biog4_c_a04
        if (p_biog1_o_a04 .ge. p1st) lptr_biog1_o_aer(04,itype,ai_phase)  = p_biog1_o_a04
        if (p_biog2_o_a04 .ge. p1st) lptr_biog2_o_aer(04,itype,ai_phase)  = p_biog2_o_a04
        if (p_biog3_o_a04 .ge. p1st) lptr_biog3_o_aer(04,itype,ai_phase)  = p_biog3_o_a04
        if (p_biog4_o_a04 .ge. p1st) lptr_biog4_o_aer(04,itype,ai_phase)  = p_biog4_o_a04
        if (p_asoaX_a04 .ge. p1st) lptr_asoaX_aer(04,itype,ai_phase)  = p_asoaX_a04
        if (p_asoa1_a04 .ge. p1st) lptr_asoa1_aer(04,itype,ai_phase)  = p_asoa1_a04
        if (p_asoa2_a04 .ge. p1st) lptr_asoa2_aer(04,itype,ai_phase)  = p_asoa2_a04
        if (p_asoa3_a04 .ge. p1st) lptr_asoa3_aer(04,itype,ai_phase)  = p_asoa3_a04
        if (p_asoa4_a04 .ge. p1st) lptr_asoa4_aer(04,itype,ai_phase)  = p_asoa4_a04
        if (p_bsoaX_a04 .ge. p1st) lptr_bsoaX_aer(04,itype,ai_phase)  = p_bsoaX_a04
        if (p_bsoa1_a04 .ge. p1st) lptr_bsoa1_aer(04,itype,ai_phase)  = p_bsoa1_a04
        if (p_bsoa2_a04 .ge. p1st) lptr_bsoa2_aer(04,itype,ai_phase)  = p_bsoa2_a04
        if (p_bsoa3_a04 .ge. p1st) lptr_bsoa3_aer(04,itype,ai_phase)  = p_bsoa3_a04
        if (p_bsoa4_a04 .ge. p1st) lptr_bsoa4_aer(04,itype,ai_phase)  = p_bsoa4_a04
        if (p_num_a04 .ge. p1st)  numptr_aer(04,itype,ai_phase)        = p_num_a04
        
        if (p_TERP_CN3_a04 .ge. p1st) lptr_terp_Cn3_aer(04,itype,ai_phase)  = p_TERP_CN3_a04
        end if


	if (nsize_aer(itype) .ge. 5) then
	    lptr_so4_aer(05,itype,ai_phase)      = p_so4_a05
	    lptr_no3_aer(05,itype,ai_phase)      = p_no3_a05
	    lptr_cl_aer(05,itype,ai_phase)       = p_cl_a05
	    lptr_msa_aer(05,itype,ai_phase)      = p_msa_a05
	    lptr_co3_aer(05,itype,ai_phase)      = p_co3_a05
	    lptr_nh4_aer(05,itype,ai_phase)      = p_nh4_a05
	    lptr_na_aer(05,itype,ai_phase)       = p_na_a05
	    lptr_ca_aer(05,itype,ai_phase)       = p_ca_a05
	    lptr_oin_aer(05,itype,ai_phase)      = p_oin_a05
	    lptr_oc_aer(05,itype,ai_phase)       = p_oc_a05
	    lptr_bc_aer(05,itype,ai_phase)       = p_bc_a05
	    hyswptr_aer(05,itype)                = p_hysw_a05
	    waterptr_aer(05,itype)               = p_water_a05
        if (p_tanv_a05 .ge. p1st) lptr_tanv_aer(05,itype,ai_phase)  = p_tanv_a05
        if (p_sopnv_a05 .ge. p1st) lptr_sopnv_aer(05,itype,ai_phase)  = p_sopnv_a05
        if (p_ternv_a05 .ge. p1st) lptr_ternv_aer(05,itype,ai_phase)  = p_ternv_a05
        if (p_seqnv_a05 .ge. p1st) lptr_seqnv_aer(05,itype,ai_phase)  = p_seqnv_a05
        if (p_sianv_a05 .ge. p1st) lptr_sianv_aer(05,itype,ai_phase)  = p_sianv_a05


            if (p_pcg1_b_c_a05 .ge. p1st) lptr_pcg1_b_c_aer(05,itype,ai_phase) = p_pcg1_b_c_a05
            if (p_pcg2_b_c_a05 .ge. p1st) lptr_pcg2_b_c_aer(05,itype,ai_phase)  = p_pcg2_b_c_a05
            if (p_pcg3_b_c_a05 .ge. p1st) lptr_pcg3_b_c_aer(05,itype,ai_phase)  = p_pcg3_b_c_a05
            if (p_pcg4_b_c_a05 .ge. p1st) lptr_pcg4_b_c_aer(05,itype,ai_phase)  = p_pcg4_b_c_a05
            if (p_pcg5_b_c_a05 .ge. p1st) lptr_pcg5_b_c_aer(05,itype,ai_phase)  = p_pcg5_b_c_a05
            if (p_pcg6_b_c_a05 .ge. p1st) lptr_pcg6_b_c_aer(05,itype,ai_phase)  = p_pcg6_b_c_a05
            if (p_pcg7_b_c_a05 .ge. p1st) lptr_pcg7_b_c_aer(05,itype,ai_phase)  = p_pcg7_b_c_a05
            if (p_pcg8_b_c_a05 .ge. p1st) lptr_pcg8_b_c_aer(05,itype,ai_phase)  = p_pcg8_b_c_a05
            if (p_pcg9_b_c_a05 .ge. p1st) lptr_pcg9_b_c_aer(05,itype,ai_phase)  = p_pcg9_b_c_a05
            if (p_pcg1_b_o_a05 .ge. p1st) lptr_pcg1_b_o_aer(05,itype,ai_phase) = p_pcg1_b_o_a05
            if (p_pcg2_b_o_a05 .ge. p1st) lptr_pcg2_b_o_aer(05,itype,ai_phase)  = p_pcg2_b_o_a05
            if (p_pcg3_b_o_a05 .ge. p1st) lptr_pcg3_b_o_aer(05,itype,ai_phase)  = p_pcg3_b_o_a05
            if (p_pcg4_b_o_a05 .ge. p1st) lptr_pcg4_b_o_aer(05,itype,ai_phase)  = p_pcg4_b_o_a05
            if (p_pcg5_b_o_a05 .ge. p1st) lptr_pcg5_b_o_aer(05,itype,ai_phase)  = p_pcg5_b_o_a05
            if (p_pcg6_b_o_a05 .ge. p1st) lptr_pcg6_b_o_aer(05,itype,ai_phase)  = p_pcg6_b_o_a05
            if (p_pcg7_b_o_a05 .ge. p1st) lptr_pcg7_b_o_aer(05,itype,ai_phase)  = p_pcg7_b_o_a05
            if (p_pcg8_b_o_a05 .ge. p1st) lptr_pcg8_b_o_aer(05,itype,ai_phase)  = p_pcg8_b_o_a05
            if (p_pcg9_b_o_a05 .ge. p1st) lptr_pcg9_b_o_aer(05,itype,ai_phase)  = p_pcg9_b_o_a05
            if (p_opcg1_b_c_a05 .ge. p1st) lptr_opcg1_b_c_aer(05,itype,ai_phase) = p_opcg1_b_c_a05
            if (p_opcg2_b_c_a05 .ge. p1st) lptr_opcg2_b_c_aer(05,itype,ai_phase)  = p_opcg2_b_c_a05
            if (p_opcg3_b_c_a05 .ge. p1st) lptr_opcg3_b_c_aer(05,itype,ai_phase)  = p_opcg3_b_c_a05
            if (p_opcg4_b_c_a05 .ge. p1st) lptr_opcg4_b_c_aer(05,itype,ai_phase)  = p_opcg4_b_c_a05
            if (p_opcg5_b_c_a05 .ge. p1st) lptr_opcg5_b_c_aer(05,itype,ai_phase)  = p_opcg5_b_c_a05
            if (p_opcg6_b_c_a05 .ge. p1st) lptr_opcg6_b_c_aer(05,itype,ai_phase)  = p_opcg6_b_c_a05
            if (p_opcg7_b_c_a05 .ge. p1st) lptr_opcg7_b_c_aer(05,itype,ai_phase)  = p_opcg7_b_c_a05
            if (p_opcg8_b_c_a05 .ge. p1st) lptr_opcg8_b_c_aer(05,itype,ai_phase)  = p_opcg8_b_c_a05
            if (p_opcg1_b_o_a05 .ge. p1st) lptr_opcg1_b_o_aer(05,itype,ai_phase) = p_opcg1_b_o_a05
            if (p_opcg2_b_o_a05 .ge. p1st) lptr_opcg2_b_o_aer(05,itype,ai_phase)  = p_opcg2_b_o_a05
            if (p_opcg3_b_o_a05 .ge. p1st) lptr_opcg3_b_o_aer(05,itype,ai_phase)  = p_opcg3_b_o_a05
            if (p_opcg4_b_o_a05 .ge. p1st) lptr_opcg4_b_o_aer(05,itype,ai_phase)  = p_opcg4_b_o_a05
            if (p_opcg5_b_o_a05 .ge. p1st) lptr_opcg5_b_o_aer(05,itype,ai_phase)  = p_opcg5_b_o_a05
            if (p_opcg6_b_o_a05 .ge. p1st) lptr_opcg6_b_o_aer(05,itype,ai_phase)  = p_opcg6_b_o_a05
            if (p_opcg7_b_o_a05 .ge. p1st) lptr_opcg7_b_o_aer(05,itype,ai_phase)  = p_opcg7_b_o_a05
            if (p_opcg8_b_o_a05 .ge. p1st) lptr_opcg8_b_o_aer(05,itype,ai_phase)  = p_opcg8_b_o_a05
            if (p_pcg1_f_c_a05 .ge. p1st) lptr_pcg1_f_c_aer(05,itype,ai_phase) = p_pcg1_f_c_a05
            if (p_pcg2_f_c_a05 .ge. p1st) lptr_pcg2_f_c_aer(05,itype,ai_phase)  = p_pcg2_f_c_a05
            if (p_pcg3_f_c_a05 .ge. p1st) lptr_pcg3_f_c_aer(05,itype,ai_phase)  = p_pcg3_f_c_a05
            if (p_pcg4_f_c_a05 .ge. p1st) lptr_pcg4_f_c_aer(05,itype,ai_phase)  = p_pcg4_f_c_a05
            if (p_pcg5_f_c_a05 .ge. p1st) lptr_pcg5_f_c_aer(05,itype,ai_phase)  = p_pcg5_f_c_a05
            if (p_pcg6_f_c_a05 .ge. p1st) lptr_pcg6_f_c_aer(05,itype,ai_phase)  = p_pcg6_f_c_a05
            if (p_pcg7_f_c_a05 .ge. p1st) lptr_pcg7_f_c_aer(05,itype,ai_phase)  = p_pcg7_f_c_a05
            if (p_pcg8_f_c_a05 .ge. p1st) lptr_pcg8_f_c_aer(05,itype,ai_phase)  = p_pcg8_f_c_a05
            if (p_pcg9_f_c_a05 .ge. p1st) lptr_pcg9_f_c_aer(05,itype,ai_phase)  = p_pcg9_f_c_a05
            if (p_pcg1_f_o_a05 .ge. p1st) lptr_pcg1_f_o_aer(05,itype,ai_phase) = p_pcg1_f_o_a05
            if (p_pcg2_f_o_a05 .ge. p1st) lptr_pcg2_f_o_aer(05,itype,ai_phase)  = p_pcg2_f_o_a05
            if (p_pcg3_f_o_a05 .ge. p1st) lptr_pcg3_f_o_aer(05,itype,ai_phase)  = p_pcg3_f_o_a05
            if (p_pcg4_f_o_a05 .ge. p1st) lptr_pcg4_f_o_aer(05,itype,ai_phase)  = p_pcg4_f_o_a05
            if (p_pcg5_f_o_a05 .ge. p1st) lptr_pcg5_f_o_aer(05,itype,ai_phase)  = p_pcg5_f_o_a05
            if (p_pcg6_f_o_a05 .ge. p1st) lptr_pcg6_f_o_aer(05,itype,ai_phase)  = p_pcg6_f_o_a05
            if (p_pcg7_f_o_a05 .ge. p1st) lptr_pcg7_f_o_aer(05,itype,ai_phase)  = p_pcg7_f_o_a05
            if (p_pcg8_f_o_a05 .ge. p1st) lptr_pcg8_f_o_aer(05,itype,ai_phase)  = p_pcg8_f_o_a05
            if (p_pcg9_f_o_a05 .ge. p1st) lptr_pcg9_f_o_aer(05,itype,ai_phase)  = p_pcg9_f_o_a05
            if (p_opcg1_f_c_a05 .ge. p1st) lptr_opcg1_f_c_aer(05,itype,ai_phase) = p_opcg1_f_c_a05
            if (p_opcg2_f_c_a05 .ge. p1st) lptr_opcg2_f_c_aer(05,itype,ai_phase)  = p_opcg2_f_c_a05
            if (p_opcg3_f_c_a05 .ge. p1st) lptr_opcg3_f_c_aer(05,itype,ai_phase)  = p_opcg3_f_c_a05
            if (p_opcg4_f_c_a05 .ge. p1st) lptr_opcg4_f_c_aer(05,itype,ai_phase)  = p_opcg4_f_c_a05
            if (p_opcg5_f_c_a05 .ge. p1st) lptr_opcg5_f_c_aer(05,itype,ai_phase)  = p_opcg5_f_c_a05
            if (p_opcg6_f_c_a05 .ge. p1st) lptr_opcg6_f_c_aer(05,itype,ai_phase)  = p_opcg6_f_c_a05
            if (p_opcg7_f_c_a05 .ge. p1st) lptr_opcg7_f_c_aer(05,itype,ai_phase)  = p_opcg7_f_c_a05
            if (p_opcg8_f_c_a05 .ge. p1st) lptr_opcg8_f_c_aer(05,itype,ai_phase)  = p_opcg8_f_c_a05
            if (p_opcg1_f_o_a05 .ge. p1st) lptr_opcg1_f_o_aer(05,itype,ai_phase) = p_opcg1_f_o_a05
            if (p_opcg2_f_o_a05 .ge. p1st) lptr_opcg2_f_o_aer(05,itype,ai_phase)  = p_opcg2_f_o_a05
            if (p_opcg3_f_o_a05 .ge. p1st) lptr_opcg3_f_o_aer(05,itype,ai_phase)  = p_opcg3_f_o_a05
            if (p_opcg4_f_o_a05 .ge. p1st) lptr_opcg4_f_o_aer(05,itype,ai_phase)  = p_opcg4_f_o_a05
            if (p_opcg5_f_o_a05 .ge. p1st) lptr_opcg5_f_o_aer(05,itype,ai_phase)  = p_opcg5_f_o_a05
            if (p_opcg6_f_o_a05 .ge. p1st) lptr_opcg6_f_o_aer(05,itype,ai_phase)  = p_opcg6_f_o_a05
            if (p_opcg7_f_o_a05 .ge. p1st) lptr_opcg7_f_o_aer(05,itype,ai_phase)  = p_opcg7_f_o_a05
            if (p_opcg8_f_o_a05 .ge. p1st) lptr_opcg8_f_o_aer(05,itype,ai_phase)  = p_opcg8_f_o_a05
            if (p_ant1_c_a05 .ge. p1st) lptr_ant1_c_aer(05,itype,ai_phase)  = p_ant1_c_a05
            if (p_ant2_c_a05 .ge. p1st) lptr_ant2_c_aer(05,itype,ai_phase)  = p_ant2_c_a05
            if (p_ant3_c_a05 .ge. p1st) lptr_ant3_c_aer(05,itype,ai_phase)  = p_ant3_c_a05
            if (p_ant4_c_a05 .ge. p1st) lptr_ant4_c_aer(05,itype,ai_phase)  = p_ant4_c_a05
            if (p_biog1_c_a05 .ge. p1st) lptr_biog1_c_aer(05,itype,ai_phase)  = p_biog1_c_a05
            if (p_biog2_c_a05 .ge. p1st) lptr_biog2_c_aer(05,itype,ai_phase)  = p_biog2_c_a05
            if (p_biog3_c_a05 .ge. p1st) lptr_biog3_c_aer(05,itype,ai_phase)  = p_biog3_c_a05
            if (p_biog4_c_a05 .ge. p1st) lptr_biog4_c_aer(05,itype,ai_phase)  = p_biog4_c_a05
            if (p_ant1_o_a05 .ge. p1st) lptr_ant1_o_aer(05,itype,ai_phase)  = p_ant1_o_a05
            if (p_ant2_o_a05 .ge. p1st) lptr_ant2_o_aer(05,itype,ai_phase)  = p_ant2_o_a05
            if (p_ant3_o_a05 .ge. p1st) lptr_ant3_o_aer(05,itype,ai_phase)  = p_ant3_o_a05
            if (p_ant4_o_a05 .ge. p1st) lptr_ant4_o_aer(05,itype,ai_phase)  = p_ant4_o_a05
            if (p_biog1_o_a05 .ge. p1st) lptr_biog1_o_aer(05,itype,ai_phase)  = p_biog1_o_a05
            if (p_biog2_o_a05 .ge. p1st) lptr_biog2_o_aer(05,itype,ai_phase)  = p_biog2_o_a05
            if (p_biog3_o_a05 .ge. p1st) lptr_biog3_o_aer(05,itype,ai_phase)  = p_biog3_o_a05
            if (p_biog4_o_a05 .ge. p1st) lptr_biog4_o_aer(05,itype,ai_phase)  = p_biog4_o_a05
        if (p_iepox_a05 .ge. p1st) lptr_iepox_aer(05,itype,ai_phase)          =p_iepox_a05
        if (p_iepoxos_a05 .ge. p1st) lptr_iepoxos_aer(05,itype,ai_phase)      =p_iepoxos_a05
        if (p_tetrol_a05 .ge. p1st) lptr_tetrol_aer(05,itype,ai_phase)        =p_tetrol_a05
        if (p_gly_a05 .ge. p1st) lptr_gly_aer(05,itype,ai_phase)              =p_gly_a05

	    numptr_aer(05,itype,ai_phase)        = p_num_a05
	end if

	if (nsize_aer(itype) .ge. 6) then
	    lptr_so4_aer(06,itype,ai_phase)      = p_so4_a06
	    lptr_no3_aer(06,itype,ai_phase)      = p_no3_a06
	    lptr_cl_aer(06,itype,ai_phase)       = p_cl_a06
	    lptr_msa_aer(06,itype,ai_phase)      = p_msa_a06
	    lptr_co3_aer(06,itype,ai_phase)      = p_co3_a06
	    lptr_nh4_aer(06,itype,ai_phase)      = p_nh4_a06
	    lptr_na_aer(06,itype,ai_phase)       = p_na_a06
	    lptr_ca_aer(06,itype,ai_phase)       = p_ca_a06
	    lptr_oin_aer(06,itype,ai_phase)      = p_oin_a06
	    lptr_oc_aer(06,itype,ai_phase)       = p_oc_a06
	    lptr_bc_aer(06,itype,ai_phase)       = p_bc_a06
	    hyswptr_aer(06,itype)                = p_hysw_a06
	    waterptr_aer(06,itype)               = p_water_a06
        if (p_tanv_a06 .ge. p1st) lptr_tanv_aer(06,itype,ai_phase)  = p_tanv_a06
        if (p_sopnv_a06 .ge. p1st) lptr_sopnv_aer(06,itype,ai_phase)  = p_sopnv_a06
        if (p_ternv_a06 .ge. p1st) lptr_ternv_aer(06,itype,ai_phase)  = p_ternv_a06
        if (p_seqnv_a06 .ge. p1st) lptr_seqnv_aer(06,itype,ai_phase)  = p_seqnv_a06
        if (p_sianv_a06 .ge. p1st) lptr_sianv_aer(06,itype,ai_phase)  = p_sianv_a06

            if (p_pcg1_b_c_a06 .ge. p1st) lptr_pcg1_b_c_aer(06,itype,ai_phase) = p_pcg1_b_c_a06
            if (p_pcg2_b_c_a06 .ge. p1st) lptr_pcg2_b_c_aer(06,itype,ai_phase)  = p_pcg2_b_c_a06
            if (p_pcg3_b_c_a06 .ge. p1st) lptr_pcg3_b_c_aer(06,itype,ai_phase)  = p_pcg3_b_c_a06
            if (p_pcg4_b_c_a06 .ge. p1st) lptr_pcg4_b_c_aer(06,itype,ai_phase)  = p_pcg4_b_c_a06
            if (p_pcg5_b_c_a06 .ge. p1st) lptr_pcg5_b_c_aer(06,itype,ai_phase)  = p_pcg5_b_c_a06
            if (p_pcg6_b_c_a06 .ge. p1st) lptr_pcg6_b_c_aer(06,itype,ai_phase)  = p_pcg6_b_c_a06
            if (p_pcg7_b_c_a06 .ge. p1st) lptr_pcg7_b_c_aer(06,itype,ai_phase)  = p_pcg7_b_c_a06
            if (p_pcg8_b_c_a06 .ge. p1st) lptr_pcg8_b_c_aer(06,itype,ai_phase)  = p_pcg8_b_c_a06
            if (p_pcg9_b_c_a06 .ge. p1st) lptr_pcg9_b_c_aer(06,itype,ai_phase)  = p_pcg9_b_c_a06
            if (p_pcg1_b_o_a06 .ge. p1st) lptr_pcg1_b_o_aer(06,itype,ai_phase) = p_pcg1_b_o_a06
            if (p_pcg2_b_o_a06 .ge. p1st) lptr_pcg2_b_o_aer(06,itype,ai_phase)  = p_pcg2_b_o_a06
            if (p_pcg3_b_o_a06 .ge. p1st) lptr_pcg3_b_o_aer(06,itype,ai_phase)  = p_pcg3_b_o_a06
            if (p_pcg4_b_o_a06 .ge. p1st) lptr_pcg4_b_o_aer(06,itype,ai_phase)  = p_pcg4_b_o_a06
            if (p_pcg5_b_o_a06 .ge. p1st) lptr_pcg5_b_o_aer(06,itype,ai_phase)  = p_pcg5_b_o_a06
            if (p_pcg6_b_o_a06 .ge. p1st) lptr_pcg6_b_o_aer(06,itype,ai_phase)  = p_pcg6_b_o_a06
            if (p_pcg7_b_o_a06 .ge. p1st) lptr_pcg7_b_o_aer(06,itype,ai_phase)  = p_pcg7_b_o_a06
            if (p_pcg8_b_o_a06 .ge. p1st) lptr_pcg8_b_o_aer(06,itype,ai_phase)  = p_pcg8_b_o_a06
            if (p_pcg9_b_o_a06 .ge. p1st) lptr_pcg9_b_o_aer(06,itype,ai_phase)  = p_pcg9_b_o_a06
            if (p_opcg1_b_c_a06 .ge. p1st) lptr_opcg1_b_c_aer(06,itype,ai_phase) = p_opcg1_b_c_a06
            if (p_opcg2_b_c_a06 .ge. p1st) lptr_opcg2_b_c_aer(06,itype,ai_phase)  = p_opcg2_b_c_a06
            if (p_opcg3_b_c_a06 .ge. p1st) lptr_opcg3_b_c_aer(06,itype,ai_phase)  = p_opcg3_b_c_a06
            if (p_opcg4_b_c_a06 .ge. p1st) lptr_opcg4_b_c_aer(06,itype,ai_phase)  = p_opcg4_b_c_a06
            if (p_opcg5_b_c_a06 .ge. p1st) lptr_opcg5_b_c_aer(06,itype,ai_phase)  = p_opcg5_b_c_a06
            if (p_opcg6_b_c_a06 .ge. p1st) lptr_opcg6_b_c_aer(06,itype,ai_phase)  = p_opcg6_b_c_a06
            if (p_opcg7_b_c_a06 .ge. p1st) lptr_opcg7_b_c_aer(06,itype,ai_phase)  = p_opcg7_b_c_a06
            if (p_opcg8_b_c_a06 .ge. p1st) lptr_opcg8_b_c_aer(06,itype,ai_phase)  = p_opcg8_b_c_a06
            if (p_opcg1_b_o_a06 .ge. p1st) lptr_opcg1_b_o_aer(06,itype,ai_phase) = p_opcg1_b_o_a06
            if (p_opcg2_b_o_a06 .ge. p1st) lptr_opcg2_b_o_aer(06,itype,ai_phase)  = p_opcg2_b_o_a06
            if (p_opcg3_b_o_a06 .ge. p1st) lptr_opcg3_b_o_aer(06,itype,ai_phase)  = p_opcg3_b_o_a06
            if (p_opcg4_b_o_a06 .ge. p1st) lptr_opcg4_b_o_aer(06,itype,ai_phase)  = p_opcg4_b_o_a06
            if (p_opcg5_b_o_a06 .ge. p1st) lptr_opcg5_b_o_aer(06,itype,ai_phase)  = p_opcg5_b_o_a06
            if (p_opcg6_b_o_a06 .ge. p1st) lptr_opcg6_b_o_aer(06,itype,ai_phase)  = p_opcg6_b_o_a06
            if (p_opcg7_b_o_a06 .ge. p1st) lptr_opcg7_b_o_aer(06,itype,ai_phase)  = p_opcg7_b_o_a06
            if (p_opcg8_b_o_a06 .ge. p1st) lptr_opcg8_b_o_aer(06,itype,ai_phase)  = p_opcg8_b_o_a06
            if (p_pcg1_f_c_a06 .ge. p1st) lptr_pcg1_f_c_aer(06,itype,ai_phase) = p_pcg1_f_c_a06
            if (p_pcg2_f_c_a06 .ge. p1st) lptr_pcg2_f_c_aer(06,itype,ai_phase)  = p_pcg2_f_c_a06
            if (p_pcg3_f_c_a06 .ge. p1st) lptr_pcg3_f_c_aer(06,itype,ai_phase)  = p_pcg3_f_c_a06
            if (p_pcg4_f_c_a06 .ge. p1st) lptr_pcg4_f_c_aer(06,itype,ai_phase)  = p_pcg4_f_c_a06
            if (p_pcg5_f_c_a06 .ge. p1st) lptr_pcg5_f_c_aer(06,itype,ai_phase)  = p_pcg5_f_c_a06
            if (p_pcg6_f_c_a06 .ge. p1st) lptr_pcg6_f_c_aer(06,itype,ai_phase)  = p_pcg6_f_c_a06
            if (p_pcg7_f_c_a06 .ge. p1st) lptr_pcg7_f_c_aer(06,itype,ai_phase)  = p_pcg7_f_c_a06
            if (p_pcg8_f_c_a06 .ge. p1st) lptr_pcg8_f_c_aer(06,itype,ai_phase)  = p_pcg8_f_c_a06
            if (p_pcg9_f_c_a06 .ge. p1st) lptr_pcg9_f_c_aer(06,itype,ai_phase)  = p_pcg9_f_c_a06
            if (p_pcg1_f_o_a06 .ge. p1st) lptr_pcg1_f_o_aer(06,itype,ai_phase) = p_pcg1_f_o_a06
            if (p_pcg2_f_o_a06 .ge. p1st) lptr_pcg2_f_o_aer(06,itype,ai_phase)  = p_pcg2_f_o_a06
            if (p_pcg3_f_o_a06 .ge. p1st) lptr_pcg3_f_o_aer(06,itype,ai_phase)  = p_pcg3_f_o_a06
            if (p_pcg4_f_o_a06 .ge. p1st) lptr_pcg4_f_o_aer(06,itype,ai_phase)  = p_pcg4_f_o_a06
            if (p_pcg5_f_o_a06 .ge. p1st) lptr_pcg5_f_o_aer(06,itype,ai_phase)  = p_pcg5_f_o_a06
            if (p_pcg6_f_o_a06 .ge. p1st) lptr_pcg6_f_o_aer(06,itype,ai_phase)  = p_pcg6_f_o_a06
            if (p_pcg7_f_o_a06 .ge. p1st) lptr_pcg7_f_o_aer(06,itype,ai_phase)  = p_pcg7_f_o_a06
            if (p_pcg8_f_o_a06 .ge. p1st) lptr_pcg8_f_o_aer(06,itype,ai_phase)  = p_pcg8_f_o_a06
            if (p_pcg9_f_o_a06 .ge. p1st) lptr_pcg9_f_o_aer(06,itype,ai_phase)  = p_pcg9_f_o_a06
            if (p_opcg1_f_c_a06 .ge. p1st) lptr_opcg1_f_c_aer(06,itype,ai_phase) = p_opcg1_f_c_a06
            if (p_opcg2_f_c_a06 .ge. p1st) lptr_opcg2_f_c_aer(06,itype,ai_phase)  = p_opcg2_f_c_a06
            if (p_opcg3_f_c_a06 .ge. p1st) lptr_opcg3_f_c_aer(06,itype,ai_phase)  = p_opcg3_f_c_a06
            if (p_opcg4_f_c_a06 .ge. p1st) lptr_opcg4_f_c_aer(06,itype,ai_phase)  = p_opcg4_f_c_a06
            if (p_opcg5_f_c_a06 .ge. p1st) lptr_opcg5_f_c_aer(06,itype,ai_phase)  = p_opcg5_f_c_a06
            if (p_opcg6_f_c_a06 .ge. p1st) lptr_opcg6_f_c_aer(06,itype,ai_phase)  = p_opcg6_f_c_a06
            if (p_opcg7_f_c_a06 .ge. p1st) lptr_opcg7_f_c_aer(06,itype,ai_phase)  = p_opcg7_f_c_a06
            if (p_opcg8_f_c_a06 .ge. p1st) lptr_opcg8_f_c_aer(06,itype,ai_phase)  = p_opcg8_f_c_a06
            if (p_opcg1_f_o_a06 .ge. p1st) lptr_opcg1_f_o_aer(06,itype,ai_phase) = p_opcg1_f_o_a06
            if (p_opcg2_f_o_a06 .ge. p1st) lptr_opcg2_f_o_aer(06,itype,ai_phase)  = p_opcg2_f_o_a06
            if (p_opcg3_f_o_a06 .ge. p1st) lptr_opcg3_f_o_aer(06,itype,ai_phase)  = p_opcg3_f_o_a06
            if (p_opcg4_f_o_a06 .ge. p1st) lptr_opcg4_f_o_aer(06,itype,ai_phase)  = p_opcg4_f_o_a06
            if (p_opcg5_f_o_a06 .ge. p1st) lptr_opcg5_f_o_aer(06,itype,ai_phase)  = p_opcg5_f_o_a06
            if (p_opcg6_f_o_a06 .ge. p1st) lptr_opcg6_f_o_aer(06,itype,ai_phase)  = p_opcg6_f_o_a06
            if (p_opcg7_f_o_a06 .ge. p1st) lptr_opcg7_f_o_aer(06,itype,ai_phase)  = p_opcg7_f_o_a06
            if (p_opcg8_f_o_a06 .ge. p1st) lptr_opcg8_f_o_aer(06,itype,ai_phase)  = p_opcg8_f_o_a06
            if (p_ant1_c_a06 .ge. p1st) lptr_ant1_c_aer(06,itype,ai_phase)  = p_ant1_c_a06
            if (p_ant2_c_a06 .ge. p1st) lptr_ant2_c_aer(06,itype,ai_phase)  = p_ant2_c_a06
            if (p_ant3_c_a06 .ge. p1st) lptr_ant3_c_aer(06,itype,ai_phase)  = p_ant3_c_a06
            if (p_ant4_c_a06 .ge. p1st) lptr_ant4_c_aer(06,itype,ai_phase)  = p_ant4_c_a06
            if (p_biog1_c_a06 .ge. p1st) lptr_biog1_c_aer(06,itype,ai_phase)  = p_biog1_c_a06
            if (p_biog2_c_a06 .ge. p1st) lptr_biog2_c_aer(06,itype,ai_phase)  = p_biog2_c_a06
            if (p_biog3_c_a06 .ge. p1st) lptr_biog3_c_aer(06,itype,ai_phase)  = p_biog3_c_a06
            if (p_biog4_c_a06 .ge. p1st) lptr_biog4_c_aer(06,itype,ai_phase)  = p_biog4_c_a06
            if (p_ant1_o_a06 .ge. p1st) lptr_ant1_o_aer(06,itype,ai_phase)  = p_ant1_o_a06
            if (p_ant2_o_a06 .ge. p1st) lptr_ant2_o_aer(06,itype,ai_phase)  = p_ant2_o_a06
            if (p_ant3_o_a06 .ge. p1st) lptr_ant3_o_aer(06,itype,ai_phase)  = p_ant3_o_a06
            if (p_ant4_o_a06 .ge. p1st) lptr_ant4_o_aer(06,itype,ai_phase)  = p_ant4_o_a06
            if (p_biog1_o_a06 .ge. p1st) lptr_biog1_o_aer(06,itype,ai_phase)  = p_biog1_o_a06
            if (p_biog2_o_a06 .ge. p1st) lptr_biog2_o_aer(06,itype,ai_phase)  = p_biog2_o_a06
            if (p_biog3_o_a06 .ge. p1st) lptr_biog3_o_aer(06,itype,ai_phase)  = p_biog3_o_a06
            if (p_biog4_o_a06 .ge. p1st) lptr_biog4_o_aer(06,itype,ai_phase)  = p_biog4_o_a06
        if (p_iepox_a06 .ge. p1st) lptr_iepox_aer(06,itype,ai_phase)          =p_iepox_a06
        if (p_iepoxos_a06 .ge. p1st) lptr_iepoxos_aer(06,itype,ai_phase)      =p_iepoxos_a06
        if (p_tetrol_a06 .ge. p1st) lptr_tetrol_aer(06,itype,ai_phase)        =p_tetrol_a06
        if (p_gly_a06 .ge. p1st) lptr_gly_aer(06,itype,ai_phase)              =p_gly_a06

	    numptr_aer(06,itype,ai_phase)        = p_num_a06
	end if

	if (nsize_aer(itype) .ge. 7) then
	    lptr_so4_aer(07,itype,ai_phase)      = p_so4_a07
	    lptr_no3_aer(07,itype,ai_phase)      = p_no3_a07
	    lptr_cl_aer(07,itype,ai_phase)       = p_cl_a07
	    lptr_msa_aer(07,itype,ai_phase)      = p_msa_a07
	    lptr_co3_aer(07,itype,ai_phase)      = p_co3_a07
	    lptr_nh4_aer(07,itype,ai_phase)      = p_nh4_a07
	    lptr_na_aer(07,itype,ai_phase)       = p_na_a07
	    lptr_ca_aer(07,itype,ai_phase)       = p_ca_a07
	    lptr_oin_aer(07,itype,ai_phase)      = p_oin_a07
	    lptr_oc_aer(07,itype,ai_phase)       = p_oc_a07
	    lptr_bc_aer(07,itype,ai_phase)       = p_bc_a07
	    hyswptr_aer(07,itype)                = p_hysw_a07
	    waterptr_aer(07,itype)               = p_water_a07
        if (p_tanv_a07 .ge. p1st) lptr_tanv_aer(07,itype,ai_phase)  = p_tanv_a07
        if (p_sopnv_a07 .ge. p1st) lptr_sopnv_aer(07,itype,ai_phase)  = p_sopnv_a07
        if (p_ternv_a07 .ge. p1st) lptr_ternv_aer(07,itype,ai_phase)  = p_ternv_a07
        if (p_seqnv_a07 .ge. p1st) lptr_seqnv_aer(07,itype,ai_phase)  = p_seqnv_a07
        if (p_sianv_a07 .ge. p1st) lptr_sianv_aer(07,itype,ai_phase)  = p_sianv_a07

            if (p_pcg1_b_c_a07 .ge. p1st) lptr_pcg1_b_c_aer(07,itype,ai_phase) = p_pcg1_b_c_a07
            if (p_pcg2_b_c_a07 .ge. p1st) lptr_pcg2_b_c_aer(07,itype,ai_phase)  = p_pcg2_b_c_a07
            if (p_pcg3_b_c_a07 .ge. p1st) lptr_pcg3_b_c_aer(07,itype,ai_phase)  = p_pcg3_b_c_a07
            if (p_pcg4_b_c_a07 .ge. p1st) lptr_pcg4_b_c_aer(07,itype,ai_phase)  = p_pcg4_b_c_a07
            if (p_pcg5_b_c_a07 .ge. p1st) lptr_pcg5_b_c_aer(07,itype,ai_phase)  = p_pcg5_b_c_a07
            if (p_pcg6_b_c_a07 .ge. p1st) lptr_pcg6_b_c_aer(07,itype,ai_phase)  = p_pcg6_b_c_a07
            if (p_pcg7_b_c_a07 .ge. p1st) lptr_pcg7_b_c_aer(07,itype,ai_phase)  = p_pcg7_b_c_a07
            if (p_pcg8_b_c_a07 .ge. p1st) lptr_pcg8_b_c_aer(07,itype,ai_phase)  = p_pcg8_b_c_a07
            if (p_pcg9_b_c_a07 .ge. p1st) lptr_pcg9_b_c_aer(07,itype,ai_phase)  = p_pcg9_b_c_a07
            if (p_pcg1_b_o_a07 .ge. p1st) lptr_pcg1_b_o_aer(07,itype,ai_phase) = p_pcg1_b_o_a07
            if (p_pcg2_b_o_a07 .ge. p1st) lptr_pcg2_b_o_aer(07,itype,ai_phase)  = p_pcg2_b_o_a07
            if (p_pcg3_b_o_a07 .ge. p1st) lptr_pcg3_b_o_aer(07,itype,ai_phase)  = p_pcg3_b_o_a07
            if (p_pcg4_b_o_a07 .ge. p1st) lptr_pcg4_b_o_aer(07,itype,ai_phase)  = p_pcg4_b_o_a07
            if (p_pcg5_b_o_a07 .ge. p1st) lptr_pcg5_b_o_aer(07,itype,ai_phase)  = p_pcg5_b_o_a07
            if (p_pcg6_b_o_a07 .ge. p1st) lptr_pcg6_b_o_aer(07,itype,ai_phase)  = p_pcg6_b_o_a07
            if (p_pcg7_b_o_a07 .ge. p1st) lptr_pcg7_b_o_aer(07,itype,ai_phase)  = p_pcg7_b_o_a07
            if (p_pcg8_b_o_a07 .ge. p1st) lptr_pcg8_b_o_aer(07,itype,ai_phase)  = p_pcg8_b_o_a07
            if (p_pcg9_b_o_a07 .ge. p1st) lptr_pcg9_b_o_aer(07,itype,ai_phase)  = p_pcg9_b_o_a07
            if (p_opcg1_b_c_a07 .ge. p1st) lptr_opcg1_b_c_aer(07,itype,ai_phase) = p_opcg1_b_c_a07
            if (p_opcg2_b_c_a07 .ge. p1st) lptr_opcg2_b_c_aer(07,itype,ai_phase)  = p_opcg2_b_c_a07
            if (p_opcg3_b_c_a07 .ge. p1st) lptr_opcg3_b_c_aer(07,itype,ai_phase)  = p_opcg3_b_c_a07
            if (p_opcg4_b_c_a07 .ge. p1st) lptr_opcg4_b_c_aer(07,itype,ai_phase)  = p_opcg4_b_c_a07
            if (p_opcg5_b_c_a07 .ge. p1st) lptr_opcg5_b_c_aer(07,itype,ai_phase)  = p_opcg5_b_c_a07
            if (p_opcg6_b_c_a07 .ge. p1st) lptr_opcg6_b_c_aer(07,itype,ai_phase)  = p_opcg6_b_c_a07
            if (p_opcg7_b_c_a07 .ge. p1st) lptr_opcg7_b_c_aer(07,itype,ai_phase)  = p_opcg7_b_c_a07
            if (p_opcg8_b_c_a07 .ge. p1st) lptr_opcg8_b_c_aer(07,itype,ai_phase)  = p_opcg8_b_c_a07
            if (p_opcg1_b_o_a07 .ge. p1st) lptr_opcg1_b_o_aer(07,itype,ai_phase) = p_opcg1_b_o_a07
            if (p_opcg2_b_o_a07 .ge. p1st) lptr_opcg2_b_o_aer(07,itype,ai_phase)  = p_opcg2_b_o_a07
            if (p_opcg3_b_o_a07 .ge. p1st) lptr_opcg3_b_o_aer(07,itype,ai_phase)  = p_opcg3_b_o_a07
            if (p_opcg4_b_o_a07 .ge. p1st) lptr_opcg4_b_o_aer(07,itype,ai_phase)  = p_opcg4_b_o_a07
            if (p_opcg5_b_o_a07 .ge. p1st) lptr_opcg5_b_o_aer(07,itype,ai_phase)  = p_opcg5_b_o_a07
            if (p_opcg6_b_o_a07 .ge. p1st) lptr_opcg6_b_o_aer(07,itype,ai_phase)  = p_opcg6_b_o_a07
            if (p_opcg7_b_o_a07 .ge. p1st) lptr_opcg7_b_o_aer(07,itype,ai_phase)  = p_opcg7_b_o_a07
            if (p_opcg8_b_o_a07 .ge. p1st) lptr_opcg8_b_o_aer(07,itype,ai_phase)  = p_opcg8_b_o_a07
            if (p_pcg1_f_c_a07 .ge. p1st) lptr_pcg1_f_c_aer(07,itype,ai_phase) = p_pcg1_f_c_a07
            if (p_pcg2_f_c_a07 .ge. p1st) lptr_pcg2_f_c_aer(07,itype,ai_phase)  = p_pcg2_f_c_a07
            if (p_pcg3_f_c_a07 .ge. p1st) lptr_pcg3_f_c_aer(07,itype,ai_phase)  = p_pcg3_f_c_a07
            if (p_pcg4_f_c_a07 .ge. p1st) lptr_pcg4_f_c_aer(07,itype,ai_phase)  = p_pcg4_f_c_a07
            if (p_pcg5_f_c_a07 .ge. p1st) lptr_pcg5_f_c_aer(07,itype,ai_phase)  = p_pcg5_f_c_a07
            if (p_pcg6_f_c_a07 .ge. p1st) lptr_pcg6_f_c_aer(07,itype,ai_phase)  = p_pcg6_f_c_a07
            if (p_pcg7_f_c_a07 .ge. p1st) lptr_pcg7_f_c_aer(07,itype,ai_phase)  = p_pcg7_f_c_a07
            if (p_pcg8_f_c_a07 .ge. p1st) lptr_pcg8_f_c_aer(07,itype,ai_phase)  = p_pcg8_f_c_a07
            if (p_pcg9_f_c_a07 .ge. p1st) lptr_pcg9_f_c_aer(07,itype,ai_phase)  = p_pcg9_f_c_a07
            if (p_pcg1_f_o_a07 .ge. p1st) lptr_pcg1_f_o_aer(07,itype,ai_phase) = p_pcg1_f_o_a07
            if (p_pcg2_f_o_a07 .ge. p1st) lptr_pcg2_f_o_aer(07,itype,ai_phase)  = p_pcg2_f_o_a07
            if (p_pcg3_f_o_a07 .ge. p1st) lptr_pcg3_f_o_aer(07,itype,ai_phase)  = p_pcg3_f_o_a07
            if (p_pcg4_f_o_a07 .ge. p1st) lptr_pcg4_f_o_aer(07,itype,ai_phase)  = p_pcg4_f_o_a07
            if (p_pcg5_f_o_a07 .ge. p1st) lptr_pcg5_f_o_aer(07,itype,ai_phase)  = p_pcg5_f_o_a07
            if (p_pcg6_f_o_a07 .ge. p1st) lptr_pcg6_f_o_aer(07,itype,ai_phase)  = p_pcg6_f_o_a07
            if (p_pcg7_f_o_a07 .ge. p1st) lptr_pcg7_f_o_aer(07,itype,ai_phase)  = p_pcg7_f_o_a07
            if (p_pcg8_f_o_a07 .ge. p1st) lptr_pcg8_f_o_aer(07,itype,ai_phase)  = p_pcg8_f_o_a07
            if (p_pcg9_f_o_a07 .ge. p1st) lptr_pcg9_f_o_aer(07,itype,ai_phase)  = p_pcg9_f_o_a07
            if (p_opcg1_f_c_a07 .ge. p1st) lptr_opcg1_f_c_aer(07,itype,ai_phase) = p_opcg1_f_c_a07
            if (p_opcg2_f_c_a07 .ge. p1st) lptr_opcg2_f_c_aer(07,itype,ai_phase)  = p_opcg2_f_c_a07
            if (p_opcg3_f_c_a07 .ge. p1st) lptr_opcg3_f_c_aer(07,itype,ai_phase)  = p_opcg3_f_c_a07
            if (p_opcg4_f_c_a07 .ge. p1st) lptr_opcg4_f_c_aer(07,itype,ai_phase)  = p_opcg4_f_c_a07
            if (p_opcg5_f_c_a07 .ge. p1st) lptr_opcg5_f_c_aer(07,itype,ai_phase)  = p_opcg5_f_c_a07
            if (p_opcg6_f_c_a07 .ge. p1st) lptr_opcg6_f_c_aer(07,itype,ai_phase)  = p_opcg6_f_c_a07
            if (p_opcg7_f_c_a07 .ge. p1st) lptr_opcg7_f_c_aer(07,itype,ai_phase)  = p_opcg7_f_c_a07
            if (p_opcg8_f_c_a07 .ge. p1st) lptr_opcg8_f_c_aer(07,itype,ai_phase)  = p_opcg8_f_c_a07
            if (p_opcg1_f_o_a07 .ge. p1st) lptr_opcg1_f_o_aer(07,itype,ai_phase) = p_opcg1_f_o_a07
            if (p_opcg2_f_o_a07 .ge. p1st) lptr_opcg2_f_o_aer(07,itype,ai_phase)  = p_opcg2_f_o_a07
            if (p_opcg3_f_o_a07 .ge. p1st) lptr_opcg3_f_o_aer(07,itype,ai_phase)  = p_opcg3_f_o_a07
            if (p_opcg4_f_o_a07 .ge. p1st) lptr_opcg4_f_o_aer(07,itype,ai_phase)  = p_opcg4_f_o_a07
            if (p_opcg5_f_o_a07 .ge. p1st) lptr_opcg5_f_o_aer(07,itype,ai_phase)  = p_opcg5_f_o_a07
            if (p_opcg6_f_o_a07 .ge. p1st) lptr_opcg6_f_o_aer(07,itype,ai_phase)  = p_opcg6_f_o_a07
            if (p_opcg7_f_o_a07 .ge. p1st) lptr_opcg7_f_o_aer(07,itype,ai_phase)  = p_opcg7_f_o_a07
            if (p_opcg8_f_o_a07 .ge. p1st) lptr_opcg8_f_o_aer(07,itype,ai_phase)  = p_opcg8_f_o_a07
            if (p_ant1_c_a07 .ge. p1st) lptr_ant1_c_aer(07,itype,ai_phase)  = p_ant1_c_a07
            if (p_ant2_c_a07 .ge. p1st) lptr_ant2_c_aer(07,itype,ai_phase)  = p_ant2_c_a07
            if (p_ant3_c_a07 .ge. p1st) lptr_ant3_c_aer(07,itype,ai_phase)  = p_ant3_c_a07
            if (p_ant4_c_a07 .ge. p1st) lptr_ant4_c_aer(07,itype,ai_phase)  = p_ant4_c_a07
            if (p_biog1_c_a07 .ge. p1st) lptr_biog1_c_aer(07,itype,ai_phase)  = p_biog1_c_a07
            if (p_biog2_c_a07 .ge. p1st) lptr_biog2_c_aer(07,itype,ai_phase)  = p_biog2_c_a07
            if (p_biog3_c_a07 .ge. p1st) lptr_biog3_c_aer(07,itype,ai_phase)  = p_biog3_c_a07
            if (p_biog4_c_a07 .ge. p1st) lptr_biog4_c_aer(07,itype,ai_phase)  = p_biog4_c_a07
            if (p_ant1_o_a07 .ge. p1st) lptr_ant1_o_aer(07,itype,ai_phase)  = p_ant1_o_a07
            if (p_ant2_o_a07 .ge. p1st) lptr_ant2_o_aer(07,itype,ai_phase)  = p_ant2_o_a07
            if (p_ant3_o_a07 .ge. p1st) lptr_ant3_o_aer(07,itype,ai_phase)  = p_ant3_o_a07
            if (p_ant4_o_a07 .ge. p1st) lptr_ant4_o_aer(07,itype,ai_phase)  = p_ant4_o_a07
            if (p_biog1_o_a07 .ge. p1st) lptr_biog1_o_aer(07,itype,ai_phase)  = p_biog1_o_a07
            if (p_biog2_o_a07 .ge. p1st) lptr_biog2_o_aer(07,itype,ai_phase)  = p_biog2_o_a07
            if (p_biog3_o_a07 .ge. p1st) lptr_biog3_o_aer(07,itype,ai_phase)  = p_biog3_o_a07
            if (p_biog4_o_a07 .ge. p1st) lptr_biog4_o_aer(07,itype,ai_phase)  = p_biog4_o_a07
        if (p_iepox_a07 .ge. p1st) lptr_iepox_aer(07,itype,ai_phase)          =p_iepox_a07
        if (p_iepoxos_a07 .ge. p1st) lptr_iepoxos_aer(07,itype,ai_phase)      =p_iepoxos_a07
        if (p_tetrol_a07 .ge. p1st) lptr_tetrol_aer(07,itype,ai_phase)        =p_tetrol_a07
        if (p_gly_a07 .ge. p1st) lptr_gly_aer(07,itype,ai_phase)              =p_gly_a07

	    numptr_aer(07,itype,ai_phase)        = p_num_a07
         end if

	if (nsize_aer(itype) .ge. 8) then
	    lptr_so4_aer(08,itype,ai_phase)      = p_so4_a08
	    lptr_no3_aer(08,itype,ai_phase)      = p_no3_a08
	    lptr_cl_aer(08,itype,ai_phase)       = p_cl_a08
	    lptr_msa_aer(08,itype,ai_phase)      = p_msa_a08
	    lptr_co3_aer(08,itype,ai_phase)      = p_co3_a08
	    lptr_nh4_aer(08,itype,ai_phase)      = p_nh4_a08
	    lptr_na_aer(08,itype,ai_phase)       = p_na_a08
	    lptr_ca_aer(08,itype,ai_phase)       = p_ca_a08
	    lptr_oin_aer(08,itype,ai_phase)      = p_oin_a08
	    lptr_oc_aer(08,itype,ai_phase)       = p_oc_a08
	    lptr_bc_aer(08,itype,ai_phase)       = p_bc_a08
	    hyswptr_aer(08,itype)                = p_hysw_a08
	    waterptr_aer(08,itype)               = p_water_a08
        if (p_tanv_a08 .ge. p1st) lptr_tanv_aer(08,itype,ai_phase)  = p_tanv_a08
        if (p_sopnv_a08 .ge. p1st) lptr_sopnv_aer(08,itype,ai_phase)  = p_sopnv_a08
        if (p_ternv_a08 .ge. p1st) lptr_ternv_aer(08,itype,ai_phase)  = p_ternv_a08
        if (p_seqnv_a08 .ge. p1st) lptr_seqnv_aer(08,itype,ai_phase)  = p_seqnv_a08
        if (p_sianv_a08 .ge. p1st) lptr_sianv_aer(08,itype,ai_phase)  = p_sianv_a08

            if (p_pcg1_b_c_a08 .ge. p1st) lptr_pcg1_b_c_aer(08,itype,ai_phase) = p_pcg1_b_c_a08
            if (p_pcg2_b_c_a08 .ge. p1st) lptr_pcg2_b_c_aer(08,itype,ai_phase)  = p_pcg2_b_c_a08
            if (p_pcg3_b_c_a08 .ge. p1st) lptr_pcg3_b_c_aer(08,itype,ai_phase)  = p_pcg3_b_c_a08
            if (p_pcg4_b_c_a08 .ge. p1st) lptr_pcg4_b_c_aer(08,itype,ai_phase)  = p_pcg4_b_c_a08
            if (p_pcg5_b_c_a08 .ge. p1st) lptr_pcg5_b_c_aer(08,itype,ai_phase)  = p_pcg5_b_c_a08
            if (p_pcg6_b_c_a08 .ge. p1st) lptr_pcg6_b_c_aer(08,itype,ai_phase)  = p_pcg6_b_c_a08
            if (p_pcg7_b_c_a08 .ge. p1st) lptr_pcg7_b_c_aer(08,itype,ai_phase)  = p_pcg7_b_c_a08
            if (p_pcg8_b_c_a08 .ge. p1st) lptr_pcg8_b_c_aer(08,itype,ai_phase)  = p_pcg8_b_c_a08
            if (p_pcg9_b_c_a08 .ge. p1st) lptr_pcg9_b_c_aer(08,itype,ai_phase)  = p_pcg9_b_c_a08
            if (p_pcg1_b_o_a08 .ge. p1st) lptr_pcg1_b_o_aer(08,itype,ai_phase) = p_pcg1_b_o_a08
            if (p_pcg2_b_o_a08 .ge. p1st) lptr_pcg2_b_o_aer(08,itype,ai_phase)  = p_pcg2_b_o_a08
            if (p_pcg3_b_o_a08 .ge. p1st) lptr_pcg3_b_o_aer(08,itype,ai_phase)  = p_pcg3_b_o_a08
            if (p_pcg4_b_o_a08 .ge. p1st) lptr_pcg4_b_o_aer(08,itype,ai_phase)  = p_pcg4_b_o_a08
            if (p_pcg5_b_o_a08 .ge. p1st) lptr_pcg5_b_o_aer(08,itype,ai_phase)  = p_pcg5_b_o_a08
            if (p_pcg6_b_o_a08 .ge. p1st) lptr_pcg6_b_o_aer(08,itype,ai_phase)  = p_pcg6_b_o_a08
            if (p_pcg7_b_o_a08 .ge. p1st) lptr_pcg7_b_o_aer(08,itype,ai_phase)  = p_pcg7_b_o_a08
            if (p_pcg8_b_o_a08 .ge. p1st) lptr_pcg8_b_o_aer(08,itype,ai_phase)  = p_pcg8_b_o_a08
            if (p_pcg9_b_o_a08 .ge. p1st) lptr_pcg9_b_o_aer(08,itype,ai_phase)  = p_pcg9_b_o_a08
            if (p_opcg1_b_c_a08 .ge. p1st) lptr_opcg1_b_c_aer(08,itype,ai_phase) = p_opcg1_b_c_a08
            if (p_opcg2_b_c_a08 .ge. p1st) lptr_opcg2_b_c_aer(08,itype,ai_phase)  = p_opcg2_b_c_a08
            if (p_opcg3_b_c_a08 .ge. p1st) lptr_opcg3_b_c_aer(08,itype,ai_phase)  = p_opcg3_b_c_a08
            if (p_opcg4_b_c_a08 .ge. p1st) lptr_opcg4_b_c_aer(08,itype,ai_phase)  = p_opcg4_b_c_a08
            if (p_opcg5_b_c_a08 .ge. p1st) lptr_opcg5_b_c_aer(08,itype,ai_phase)  = p_opcg5_b_c_a08
            if (p_opcg6_b_c_a08 .ge. p1st) lptr_opcg6_b_c_aer(08,itype,ai_phase)  = p_opcg6_b_c_a08
            if (p_opcg7_b_c_a08 .ge. p1st) lptr_opcg7_b_c_aer(08,itype,ai_phase)  = p_opcg7_b_c_a08
            if (p_opcg8_b_c_a08 .ge. p1st) lptr_opcg8_b_c_aer(08,itype,ai_phase)  = p_opcg8_b_c_a08
            if (p_opcg1_b_o_a08 .ge. p1st) lptr_opcg1_b_o_aer(08,itype,ai_phase) = p_opcg1_b_o_a08
            if (p_opcg2_b_o_a08 .ge. p1st) lptr_opcg2_b_o_aer(08,itype,ai_phase)  = p_opcg2_b_o_a08
            if (p_opcg3_b_o_a08 .ge. p1st) lptr_opcg3_b_o_aer(08,itype,ai_phase)  = p_opcg3_b_o_a08
            if (p_opcg4_b_o_a08 .ge. p1st) lptr_opcg4_b_o_aer(08,itype,ai_phase)  = p_opcg4_b_o_a08
            if (p_opcg5_b_o_a08 .ge. p1st) lptr_opcg5_b_o_aer(08,itype,ai_phase)  = p_opcg5_b_o_a08
            if (p_opcg6_b_o_a08 .ge. p1st) lptr_opcg6_b_o_aer(08,itype,ai_phase)  = p_opcg6_b_o_a08
            if (p_opcg7_b_o_a08 .ge. p1st) lptr_opcg7_b_o_aer(08,itype,ai_phase)  = p_opcg7_b_o_a08
            if (p_opcg8_b_o_a08 .ge. p1st) lptr_opcg8_b_o_aer(08,itype,ai_phase)  = p_opcg8_b_o_a08
            if (p_pcg1_f_c_a08 .ge. p1st) lptr_pcg1_f_c_aer(08,itype,ai_phase) = p_pcg1_f_c_a08
            if (p_pcg2_f_c_a08 .ge. p1st) lptr_pcg2_f_c_aer(08,itype,ai_phase)  = p_pcg2_f_c_a08
            if (p_pcg3_f_c_a08 .ge. p1st) lptr_pcg3_f_c_aer(08,itype,ai_phase)  = p_pcg3_f_c_a08
            if (p_pcg4_f_c_a08 .ge. p1st) lptr_pcg4_f_c_aer(08,itype,ai_phase)  = p_pcg4_f_c_a08
            if (p_pcg5_f_c_a08 .ge. p1st) lptr_pcg5_f_c_aer(08,itype,ai_phase)  = p_pcg5_f_c_a08
            if (p_pcg6_f_c_a08 .ge. p1st) lptr_pcg6_f_c_aer(08,itype,ai_phase)  = p_pcg6_f_c_a08
            if (p_pcg7_f_c_a08 .ge. p1st) lptr_pcg7_f_c_aer(08,itype,ai_phase)  = p_pcg7_f_c_a08
            if (p_pcg8_f_c_a08 .ge. p1st) lptr_pcg8_f_c_aer(08,itype,ai_phase)  = p_pcg8_f_c_a08
            if (p_pcg9_f_c_a08 .ge. p1st) lptr_pcg9_f_c_aer(08,itype,ai_phase)  = p_pcg9_f_c_a08
            if (p_pcg1_f_o_a08 .ge. p1st) lptr_pcg1_f_o_aer(08,itype,ai_phase) = p_pcg1_f_o_a08
            if (p_pcg2_f_o_a08 .ge. p1st) lptr_pcg2_f_o_aer(08,itype,ai_phase)  = p_pcg2_f_o_a08
            if (p_pcg3_f_o_a08 .ge. p1st) lptr_pcg3_f_o_aer(08,itype,ai_phase)  = p_pcg3_f_o_a08
            if (p_pcg4_f_o_a08 .ge. p1st) lptr_pcg4_f_o_aer(08,itype,ai_phase)  = p_pcg4_f_o_a08
            if (p_pcg5_f_o_a08 .ge. p1st) lptr_pcg5_f_o_aer(08,itype,ai_phase)  = p_pcg5_f_o_a08
            if (p_pcg6_f_o_a08 .ge. p1st) lptr_pcg6_f_o_aer(08,itype,ai_phase)  = p_pcg6_f_o_a08
            if (p_pcg7_f_o_a08 .ge. p1st) lptr_pcg7_f_o_aer(08,itype,ai_phase)  = p_pcg7_f_o_a08
            if (p_pcg8_f_o_a08 .ge. p1st) lptr_pcg8_f_o_aer(08,itype,ai_phase)  = p_pcg8_f_o_a08
            if (p_pcg9_f_o_a08 .ge. p1st) lptr_pcg9_f_o_aer(08,itype,ai_phase)  = p_pcg9_f_o_a08
            if (p_opcg1_f_c_a08 .ge. p1st) lptr_opcg1_f_c_aer(08,itype,ai_phase) = p_opcg1_f_c_a08
            if (p_opcg2_f_c_a08 .ge. p1st) lptr_opcg2_f_c_aer(08,itype,ai_phase)  = p_opcg2_f_c_a08
            if (p_opcg3_f_c_a08 .ge. p1st) lptr_opcg3_f_c_aer(08,itype,ai_phase)  = p_opcg3_f_c_a08
            if (p_opcg4_f_c_a08 .ge. p1st) lptr_opcg4_f_c_aer(08,itype,ai_phase)  = p_opcg4_f_c_a08
            if (p_opcg5_f_c_a08 .ge. p1st) lptr_opcg5_f_c_aer(08,itype,ai_phase)  = p_opcg5_f_c_a08
            if (p_opcg6_f_c_a08 .ge. p1st) lptr_opcg6_f_c_aer(08,itype,ai_phase)  = p_opcg6_f_c_a08
            if (p_opcg7_f_c_a08 .ge. p1st) lptr_opcg7_f_c_aer(08,itype,ai_phase)  = p_opcg7_f_c_a08
            if (p_opcg8_f_c_a08 .ge. p1st) lptr_opcg8_f_c_aer(08,itype,ai_phase)  = p_opcg8_f_c_a08
            if (p_opcg1_f_o_a08 .ge. p1st) lptr_opcg1_f_o_aer(08,itype,ai_phase) = p_opcg1_f_o_a08
            if (p_opcg2_f_o_a08 .ge. p1st) lptr_opcg2_f_o_aer(08,itype,ai_phase)  = p_opcg2_f_o_a08
            if (p_opcg3_f_o_a08 .ge. p1st) lptr_opcg3_f_o_aer(08,itype,ai_phase)  = p_opcg3_f_o_a08
            if (p_opcg4_f_o_a08 .ge. p1st) lptr_opcg4_f_o_aer(08,itype,ai_phase)  = p_opcg4_f_o_a08
            if (p_opcg5_f_o_a08 .ge. p1st) lptr_opcg5_f_o_aer(08,itype,ai_phase)  = p_opcg5_f_o_a08
            if (p_opcg6_f_o_a08 .ge. p1st) lptr_opcg6_f_o_aer(08,itype,ai_phase)  = p_opcg6_f_o_a08
            if (p_opcg7_f_o_a08 .ge. p1st) lptr_opcg7_f_o_aer(08,itype,ai_phase)  = p_opcg7_f_o_a08
            if (p_opcg8_f_o_a08 .ge. p1st) lptr_opcg8_f_o_aer(08,itype,ai_phase)  = p_opcg8_f_o_a08
            if (p_ant1_c_a08 .ge. p1st) lptr_ant1_c_aer(08,itype,ai_phase)  = p_ant1_c_a08
            if (p_ant2_c_a08 .ge. p1st) lptr_ant2_c_aer(08,itype,ai_phase)  = p_ant2_c_a08
            if (p_ant3_c_a08 .ge. p1st) lptr_ant3_c_aer(08,itype,ai_phase)  = p_ant3_c_a08
            if (p_ant4_c_a08 .ge. p1st) lptr_ant4_c_aer(08,itype,ai_phase)  = p_ant4_c_a08
            if (p_biog1_c_a08 .ge. p1st) lptr_biog1_c_aer(08,itype,ai_phase)  = p_biog1_c_a08
            if (p_biog2_c_a08 .ge. p1st) lptr_biog2_c_aer(08,itype,ai_phase)  = p_biog2_c_a08
            if (p_biog3_c_a08 .ge. p1st) lptr_biog3_c_aer(08,itype,ai_phase)  = p_biog3_c_a08
            if (p_biog4_c_a08 .ge. p1st) lptr_biog4_c_aer(08,itype,ai_phase)  = p_biog4_c_a08
            if (p_ant1_o_a08 .ge. p1st) lptr_ant1_o_aer(08,itype,ai_phase)  = p_ant1_o_a08
            if (p_ant2_o_a08 .ge. p1st) lptr_ant2_o_aer(08,itype,ai_phase)  = p_ant2_o_a08
            if (p_ant3_o_a08 .ge. p1st) lptr_ant3_o_aer(08,itype,ai_phase)  = p_ant3_o_a08
            if (p_ant4_o_a08 .ge. p1st) lptr_ant4_o_aer(08,itype,ai_phase)  = p_ant4_o_a08
            if (p_biog1_o_a08 .ge. p1st) lptr_biog1_o_aer(08,itype,ai_phase)  = p_biog1_o_a08
            if (p_biog2_o_a08 .ge. p1st) lptr_biog2_o_aer(08,itype,ai_phase)  = p_biog2_o_a08
            if (p_biog3_o_a08 .ge. p1st) lptr_biog3_o_aer(08,itype,ai_phase)  = p_biog3_o_a08
            if (p_biog4_o_a08 .ge. p1st) lptr_biog4_o_aer(08,itype,ai_phase)  = p_biog4_o_a08
        if (p_iepox_a08 .ge. p1st) lptr_iepox_aer(08,itype,ai_phase)          =p_iepox_a08
        if (p_iepoxos_a08 .ge. p1st) lptr_iepoxos_aer(08,itype,ai_phase)      =p_iepoxos_a08
        if (p_tetrol_a08 .ge. p1st) lptr_tetrol_aer(08,itype,ai_phase)        =p_tetrol_a08
        if (p_gly_a08 .ge. p1st) lptr_gly_aer(08,itype,ai_phase)              =p_gly_a08

            
	    numptr_aer(08,itype,ai_phase)        = p_num_a08
	end if

	if (nsize_aer(itype) .ge. 9) then
	    lptr_so4_aer(09,itype,ai_phase)      = p_so4_a09
	    lptr_no3_aer(09,itype,ai_phase)      = p_no3_a09
	    lptr_cl_aer(09,itype,ai_phase)       = p_cl_a09
	    lptr_msa_aer(09,itype,ai_phase)      = p_msa_a09
	    lptr_co3_aer(09,itype,ai_phase)      = p_co3_a09
	    lptr_nh4_aer(09,itype,ai_phase)      = p_nh4_a09
	    lptr_na_aer(09,itype,ai_phase)       = p_na_a09
	    lptr_ca_aer(09,itype,ai_phase)       = p_ca_a09
	    lptr_oin_aer(09,itype,ai_phase)      = p_oin_a09
	    lptr_oc_aer(09,itype,ai_phase)       = p_oc_a09
	    lptr_bc_aer(09,itype,ai_phase)       = p_bc_a09
	    hyswptr_aer(09,itype)                = p_hysw_a09
	    waterptr_aer(09,itype)               = p_water_a09

            if (p_pcg1_b_c_a09 .ge. p1st) lptr_pcg1_b_c_aer(09,itype,ai_phase) = p_pcg1_b_c_a09
            if (p_pcg2_b_c_a09 .ge. p1st) lptr_pcg2_b_c_aer(09,itype,ai_phase)  = p_pcg2_b_c_a09
            if (p_pcg3_b_c_a09 .ge. p1st) lptr_pcg3_b_c_aer(09,itype,ai_phase)  = p_pcg3_b_c_a09
            if (p_pcg4_b_c_a09 .ge. p1st) lptr_pcg4_b_c_aer(09,itype,ai_phase)  = p_pcg4_b_c_a09
            if (p_pcg5_b_c_a09 .ge. p1st) lptr_pcg5_b_c_aer(09,itype,ai_phase)  = p_pcg5_b_c_a09
            if (p_pcg6_b_c_a09 .ge. p1st) lptr_pcg6_b_c_aer(09,itype,ai_phase)  = p_pcg6_b_c_a09
            if (p_pcg7_b_c_a09 .ge. p1st) lptr_pcg7_b_c_aer(09,itype,ai_phase)  = p_pcg7_b_c_a09
            if (p_pcg8_b_c_a09 .ge. p1st) lptr_pcg8_b_c_aer(09,itype,ai_phase)  = p_pcg8_b_c_a09
            if (p_pcg9_b_c_a09 .ge. p1st) lptr_pcg9_b_c_aer(09,itype,ai_phase)  = p_pcg9_b_c_a09
            if (p_pcg1_b_o_a09 .ge. p1st) lptr_pcg1_b_o_aer(09,itype,ai_phase) = p_pcg1_b_o_a09
            if (p_pcg2_b_o_a09 .ge. p1st) lptr_pcg2_b_o_aer(09,itype,ai_phase)  = p_pcg2_b_o_a09
            if (p_pcg3_b_o_a09 .ge. p1st) lptr_pcg3_b_o_aer(09,itype,ai_phase)  = p_pcg3_b_o_a09
            if (p_pcg4_b_o_a09 .ge. p1st) lptr_pcg4_b_o_aer(09,itype,ai_phase)  = p_pcg4_b_o_a09
            if (p_pcg5_b_o_a09 .ge. p1st) lptr_pcg5_b_o_aer(09,itype,ai_phase)  = p_pcg5_b_o_a09
            if (p_pcg6_b_o_a09 .ge. p1st) lptr_pcg6_b_o_aer(09,itype,ai_phase)  = p_pcg6_b_o_a09
            if (p_pcg7_b_o_a09 .ge. p1st) lptr_pcg7_b_o_aer(09,itype,ai_phase)  = p_pcg7_b_o_a09
            if (p_pcg8_b_o_a09 .ge. p1st) lptr_pcg8_b_o_aer(09,itype,ai_phase)  = p_pcg8_b_o_a09
            if (p_pcg9_b_o_a09 .ge. p1st) lptr_pcg9_b_o_aer(09,itype,ai_phase)  = p_pcg9_b_o_a09
            if (p_opcg1_b_c_a09 .ge. p1st) lptr_opcg1_b_c_aer(09,itype,ai_phase) = p_opcg1_b_c_a09
            if (p_opcg2_b_c_a09 .ge. p1st) lptr_opcg2_b_c_aer(09,itype,ai_phase)  = p_opcg2_b_c_a09
            if (p_opcg3_b_c_a09 .ge. p1st) lptr_opcg3_b_c_aer(09,itype,ai_phase)  = p_opcg3_b_c_a09
            if (p_opcg4_b_c_a09 .ge. p1st) lptr_opcg4_b_c_aer(09,itype,ai_phase)  = p_opcg4_b_c_a09
            if (p_opcg5_b_c_a09 .ge. p1st) lptr_opcg5_b_c_aer(09,itype,ai_phase)  = p_opcg5_b_c_a09
            if (p_opcg6_b_c_a09 .ge. p1st) lptr_opcg6_b_c_aer(09,itype,ai_phase)  = p_opcg6_b_c_a09
            if (p_opcg7_b_c_a09 .ge. p1st) lptr_opcg7_b_c_aer(09,itype,ai_phase)  = p_opcg7_b_c_a09
            if (p_opcg8_b_c_a09 .ge. p1st) lptr_opcg8_b_c_aer(09,itype,ai_phase)  = p_opcg8_b_c_a09
            if (p_opcg1_b_o_a09 .ge. p1st) lptr_opcg1_b_o_aer(09,itype,ai_phase) = p_opcg1_b_o_a09
            if (p_opcg2_b_o_a09 .ge. p1st) lptr_opcg2_b_o_aer(09,itype,ai_phase)  = p_opcg2_b_o_a09
            if (p_opcg3_b_o_a09 .ge. p1st) lptr_opcg3_b_o_aer(09,itype,ai_phase)  = p_opcg3_b_o_a09
            if (p_opcg4_b_o_a09 .ge. p1st) lptr_opcg4_b_o_aer(09,itype,ai_phase)  = p_opcg4_b_o_a09
            if (p_opcg5_b_o_a09 .ge. p1st) lptr_opcg5_b_o_aer(09,itype,ai_phase)  = p_opcg5_b_o_a09
            if (p_opcg6_b_o_a09 .ge. p1st) lptr_opcg6_b_o_aer(09,itype,ai_phase)  = p_opcg6_b_o_a09
            if (p_opcg7_b_o_a09 .ge. p1st) lptr_opcg7_b_o_aer(09,itype,ai_phase)  = p_opcg7_b_o_a09
            if (p_opcg8_b_o_a09 .ge. p1st) lptr_opcg8_b_o_aer(09,itype,ai_phase)  = p_opcg8_b_o_a09
            if (p_pcg1_f_c_a09 .ge. p1st) lptr_pcg1_f_c_aer(09,itype,ai_phase) = p_pcg1_f_c_a09
            if (p_pcg2_f_c_a09 .ge. p1st) lptr_pcg2_f_c_aer(09,itype,ai_phase)  = p_pcg2_f_c_a09
            if (p_pcg3_f_c_a09 .ge. p1st) lptr_pcg3_f_c_aer(09,itype,ai_phase)  = p_pcg3_f_c_a09
            if (p_pcg4_f_c_a09 .ge. p1st) lptr_pcg4_f_c_aer(09,itype,ai_phase)  = p_pcg4_f_c_a09
            if (p_pcg5_f_c_a09 .ge. p1st) lptr_pcg5_f_c_aer(09,itype,ai_phase)  = p_pcg5_f_c_a09
            if (p_pcg6_f_c_a09 .ge. p1st) lptr_pcg6_f_c_aer(09,itype,ai_phase)  = p_pcg6_f_c_a09
            if (p_pcg7_f_c_a09 .ge. p1st) lptr_pcg7_f_c_aer(09,itype,ai_phase)  = p_pcg7_f_c_a09
            if (p_pcg8_f_c_a09 .ge. p1st) lptr_pcg8_f_c_aer(09,itype,ai_phase)  = p_pcg8_f_c_a09
            if (p_pcg9_f_c_a09 .ge. p1st) lptr_pcg9_f_c_aer(09,itype,ai_phase)  = p_pcg9_f_c_a09
            if (p_pcg1_f_o_a09 .ge. p1st) lptr_pcg1_f_o_aer(09,itype,ai_phase) = p_pcg1_f_o_a09
            if (p_pcg2_f_o_a09 .ge. p1st) lptr_pcg2_f_o_aer(09,itype,ai_phase)  = p_pcg2_f_o_a09
            if (p_pcg3_f_o_a09 .ge. p1st) lptr_pcg3_f_o_aer(09,itype,ai_phase)  = p_pcg3_f_o_a09
            if (p_pcg4_f_o_a09 .ge. p1st) lptr_pcg4_f_o_aer(09,itype,ai_phase)  = p_pcg4_f_o_a09
            if (p_pcg5_f_o_a09 .ge. p1st) lptr_pcg5_f_o_aer(09,itype,ai_phase)  = p_pcg5_f_o_a09
            if (p_pcg6_f_o_a09 .ge. p1st) lptr_pcg6_f_o_aer(09,itype,ai_phase)  = p_pcg6_f_o_a09
            if (p_pcg7_f_o_a09 .ge. p1st) lptr_pcg7_f_o_aer(09,itype,ai_phase)  = p_pcg7_f_o_a09
            if (p_pcg8_f_o_a09 .ge. p1st) lptr_pcg8_f_o_aer(09,itype,ai_phase)  = p_pcg8_f_o_a09
            if (p_pcg9_f_o_a09 .ge. p1st) lptr_pcg9_f_o_aer(09,itype,ai_phase)  = p_pcg9_f_o_a09
            if (p_opcg1_f_c_a09 .ge. p1st) lptr_opcg1_f_c_aer(09,itype,ai_phase) = p_opcg1_f_c_a09
            if (p_opcg2_f_c_a09 .ge. p1st) lptr_opcg2_f_c_aer(09,itype,ai_phase)  = p_opcg2_f_c_a09
            if (p_opcg3_f_c_a09 .ge. p1st) lptr_opcg3_f_c_aer(09,itype,ai_phase)  = p_opcg3_f_c_a09
            if (p_opcg4_f_c_a09 .ge. p1st) lptr_opcg4_f_c_aer(09,itype,ai_phase)  = p_opcg4_f_c_a09
            if (p_opcg5_f_c_a09 .ge. p1st) lptr_opcg5_f_c_aer(09,itype,ai_phase)  = p_opcg5_f_c_a09
            if (p_opcg6_f_c_a09 .ge. p1st) lptr_opcg6_f_c_aer(09,itype,ai_phase)  = p_opcg6_f_c_a09
            if (p_opcg7_f_c_a09 .ge. p1st) lptr_opcg7_f_c_aer(09,itype,ai_phase)  = p_opcg7_f_c_a09
            if (p_opcg8_f_c_a09 .ge. p1st) lptr_opcg8_f_c_aer(09,itype,ai_phase)  = p_opcg8_f_c_a09
            if (p_opcg1_f_o_a09 .ge. p1st) lptr_opcg1_f_o_aer(09,itype,ai_phase) = p_opcg1_f_o_a09
            if (p_opcg2_f_o_a09 .ge. p1st) lptr_opcg2_f_o_aer(09,itype,ai_phase)  = p_opcg2_f_o_a09
            if (p_opcg3_f_o_a09 .ge. p1st) lptr_opcg3_f_o_aer(09,itype,ai_phase)  = p_opcg3_f_o_a09
            if (p_opcg4_f_o_a09 .ge. p1st) lptr_opcg4_f_o_aer(09,itype,ai_phase)  = p_opcg4_f_o_a09
            if (p_opcg5_f_o_a09 .ge. p1st) lptr_opcg5_f_o_aer(09,itype,ai_phase)  = p_opcg5_f_o_a09
            if (p_opcg6_f_o_a09 .ge. p1st) lptr_opcg6_f_o_aer(09,itype,ai_phase)  = p_opcg6_f_o_a09
            if (p_opcg7_f_o_a09 .ge. p1st) lptr_opcg7_f_o_aer(09,itype,ai_phase)  = p_opcg7_f_o_a09
            if (p_opcg8_f_o_a09 .ge. p1st) lptr_opcg8_f_o_aer(09,itype,ai_phase)  = p_opcg8_f_o_a09
            if (p_ant1_c_a09 .ge. p1st) lptr_ant1_c_aer(09,itype,ai_phase)  = p_ant1_c_a09
            if (p_ant2_c_a09 .ge. p1st) lptr_ant2_c_aer(09,itype,ai_phase)  = p_ant2_c_a09
            if (p_ant3_c_a09 .ge. p1st) lptr_ant3_c_aer(09,itype,ai_phase)  = p_ant3_c_a09
            if (p_ant4_c_a09 .ge. p1st) lptr_ant4_c_aer(09,itype,ai_phase)  = p_ant4_c_a09
            if (p_biog1_c_a09 .ge. p1st) lptr_biog1_c_aer(09,itype,ai_phase)  = p_biog1_c_a09
            if (p_biog2_c_a09 .ge. p1st) lptr_biog2_c_aer(09,itype,ai_phase)  = p_biog2_c_a09
            if (p_biog3_c_a09 .ge. p1st) lptr_biog3_c_aer(09,itype,ai_phase)  = p_biog3_c_a09
            if (p_biog4_c_a09 .ge. p1st) lptr_biog4_c_aer(09,itype,ai_phase)  = p_biog4_c_a09
            if (p_ant1_o_a09 .ge. p1st) lptr_ant1_o_aer(09,itype,ai_phase)  = p_ant1_o_a09
            if (p_ant2_o_a09 .ge. p1st) lptr_ant2_o_aer(09,itype,ai_phase)  = p_ant2_o_a09
            if (p_ant3_o_a09 .ge. p1st) lptr_ant3_o_aer(09,itype,ai_phase)  = p_ant3_o_a09
            if (p_ant4_o_a09 .ge. p1st) lptr_ant4_o_aer(09,itype,ai_phase)  = p_ant4_o_a09
            if (p_biog1_o_a09 .ge. p1st) lptr_biog1_o_aer(09,itype,ai_phase)  = p_biog1_o_a09
            if (p_biog2_o_a09 .ge. p1st) lptr_biog2_o_aer(09,itype,ai_phase)  = p_biog2_o_a09
            if (p_biog3_o_a09 .ge. p1st) lptr_biog3_o_aer(09,itype,ai_phase)  = p_biog3_o_a09
            if (p_biog4_o_a09 .ge. p1st) lptr_biog4_o_aer(09,itype,ai_phase)  = p_biog4_o_a09
            
	    numptr_aer(09,itype,ai_phase)        = p_num_a09
	end if

	if (nsize_aer(itype) .ge. 10) then
	    lptr_so4_aer(10,itype,ai_phase)      = p_so4_a10
	    lptr_no3_aer(10,itype,ai_phase)      = p_no3_a10
	    lptr_cl_aer(10,itype,ai_phase)       = p_cl_a10
	    lptr_msa_aer(10,itype,ai_phase)      = p_msa_a10
	    lptr_co3_aer(10,itype,ai_phase)      = p_co3_a10
	    lptr_nh4_aer(10,itype,ai_phase)      = p_nh4_a10
	    lptr_na_aer(10,itype,ai_phase)       = p_na_a10
	    lptr_ca_aer(10,itype,ai_phase)       = p_ca_a10
	    lptr_oin_aer(10,itype,ai_phase)      = p_oin_a10
	    lptr_oc_aer(10,itype,ai_phase)       = p_oc_a10
	    lptr_bc_aer(10,itype,ai_phase)       = p_bc_a10
	    hyswptr_aer(10,itype)                = p_hysw_a10
	    waterptr_aer(10,itype)               = p_water_a10

            if (p_pcg1_b_c_a10 .ge. p1st) lptr_pcg1_b_c_aer(10,itype,ai_phase) = p_pcg1_b_c_a10
            if (p_pcg2_b_c_a10 .ge. p1st) lptr_pcg2_b_c_aer(10,itype,ai_phase)  = p_pcg2_b_c_a10
            if (p_pcg3_b_c_a10 .ge. p1st) lptr_pcg3_b_c_aer(10,itype,ai_phase)  = p_pcg3_b_c_a10
            if (p_pcg4_b_c_a10 .ge. p1st) lptr_pcg4_b_c_aer(10,itype,ai_phase)  = p_pcg4_b_c_a10
            if (p_pcg5_b_c_a10 .ge. p1st) lptr_pcg5_b_c_aer(10,itype,ai_phase)  = p_pcg5_b_c_a10
            if (p_pcg6_b_c_a10 .ge. p1st) lptr_pcg6_b_c_aer(10,itype,ai_phase)  = p_pcg6_b_c_a10
            if (p_pcg7_b_c_a10 .ge. p1st) lptr_pcg7_b_c_aer(10,itype,ai_phase)  = p_pcg7_b_c_a10
            if (p_pcg8_b_c_a10 .ge. p1st) lptr_pcg8_b_c_aer(10,itype,ai_phase)  = p_pcg8_b_c_a10
            if (p_pcg9_b_c_a10 .ge. p1st) lptr_pcg9_b_c_aer(10,itype,ai_phase)  = p_pcg9_b_c_a10
            if (p_pcg1_b_o_a10 .ge. p1st) lptr_pcg1_b_o_aer(10,itype,ai_phase) = p_pcg1_b_o_a10
            if (p_pcg2_b_o_a10 .ge. p1st) lptr_pcg2_b_o_aer(10,itype,ai_phase)  = p_pcg2_b_o_a10
            if (p_pcg3_b_o_a10 .ge. p1st) lptr_pcg3_b_o_aer(10,itype,ai_phase)  = p_pcg3_b_o_a10
            if (p_pcg4_b_o_a10 .ge. p1st) lptr_pcg4_b_o_aer(10,itype,ai_phase)  = p_pcg4_b_o_a10
            if (p_pcg5_b_o_a10 .ge. p1st) lptr_pcg5_b_o_aer(10,itype,ai_phase)  = p_pcg5_b_o_a10
            if (p_pcg6_b_o_a10 .ge. p1st) lptr_pcg6_b_o_aer(10,itype,ai_phase)  = p_pcg6_b_o_a10
            if (p_pcg7_b_o_a10 .ge. p1st) lptr_pcg7_b_o_aer(10,itype,ai_phase)  = p_pcg7_b_o_a10
            if (p_pcg8_b_o_a10 .ge. p1st) lptr_pcg8_b_o_aer(10,itype,ai_phase)  = p_pcg8_b_o_a10
            if (p_pcg9_b_o_a10 .ge. p1st) lptr_pcg9_b_o_aer(10,itype,ai_phase)  = p_pcg9_b_o_a10
            if (p_opcg1_b_c_a10 .ge. p1st) lptr_opcg1_b_c_aer(10,itype,ai_phase) = p_opcg1_b_c_a10
            if (p_opcg2_b_c_a10 .ge. p1st) lptr_opcg2_b_c_aer(10,itype,ai_phase)  = p_opcg2_b_c_a10
            if (p_opcg3_b_c_a10 .ge. p1st) lptr_opcg3_b_c_aer(10,itype,ai_phase)  = p_opcg3_b_c_a10
            if (p_opcg4_b_c_a10 .ge. p1st) lptr_opcg4_b_c_aer(10,itype,ai_phase)  = p_opcg4_b_c_a10
            if (p_opcg5_b_c_a10 .ge. p1st) lptr_opcg5_b_c_aer(10,itype,ai_phase)  = p_opcg5_b_c_a10
            if (p_opcg6_b_c_a10 .ge. p1st) lptr_opcg6_b_c_aer(10,itype,ai_phase)  = p_opcg6_b_c_a10
            if (p_opcg7_b_c_a10 .ge. p1st) lptr_opcg7_b_c_aer(10,itype,ai_phase)  = p_opcg7_b_c_a10
            if (p_opcg8_b_c_a10 .ge. p1st) lptr_opcg8_b_c_aer(10,itype,ai_phase)  = p_opcg8_b_c_a10
            if (p_opcg1_b_o_a10 .ge. p1st) lptr_opcg1_b_o_aer(10,itype,ai_phase) = p_opcg1_b_o_a10
            if (p_opcg2_b_o_a10 .ge. p1st) lptr_opcg2_b_o_aer(10,itype,ai_phase)  = p_opcg2_b_o_a10
            if (p_opcg3_b_o_a10 .ge. p1st) lptr_opcg3_b_o_aer(10,itype,ai_phase)  = p_opcg3_b_o_a10
            if (p_opcg4_b_o_a10 .ge. p1st) lptr_opcg4_b_o_aer(10,itype,ai_phase)  = p_opcg4_b_o_a10
            if (p_opcg5_b_o_a10 .ge. p1st) lptr_opcg5_b_o_aer(10,itype,ai_phase)  = p_opcg5_b_o_a10
            if (p_opcg6_b_o_a10 .ge. p1st) lptr_opcg6_b_o_aer(10,itype,ai_phase)  = p_opcg6_b_o_a10
            if (p_opcg7_b_o_a10 .ge. p1st) lptr_opcg7_b_o_aer(10,itype,ai_phase)  = p_opcg7_b_o_a10
            if (p_opcg8_b_o_a10 .ge. p1st) lptr_opcg8_b_o_aer(10,itype,ai_phase)  = p_opcg8_b_o_a10
            if (p_pcg1_f_c_a10 .ge. p1st) lptr_pcg1_f_c_aer(10,itype,ai_phase) = p_pcg1_f_c_a10
            if (p_pcg2_f_c_a10 .ge. p1st) lptr_pcg2_f_c_aer(10,itype,ai_phase)  = p_pcg2_f_c_a10
            if (p_pcg3_f_c_a10 .ge. p1st) lptr_pcg3_f_c_aer(10,itype,ai_phase)  = p_pcg3_f_c_a10
            if (p_pcg4_f_c_a10 .ge. p1st) lptr_pcg4_f_c_aer(10,itype,ai_phase)  = p_pcg4_f_c_a10
            if (p_pcg5_f_c_a10 .ge. p1st) lptr_pcg5_f_c_aer(10,itype,ai_phase)  = p_pcg5_f_c_a10
            if (p_pcg6_f_c_a10 .ge. p1st) lptr_pcg6_f_c_aer(10,itype,ai_phase)  = p_pcg6_f_c_a10
            if (p_pcg7_f_c_a10 .ge. p1st) lptr_pcg7_f_c_aer(10,itype,ai_phase)  = p_pcg7_f_c_a10
            if (p_pcg8_f_c_a10 .ge. p1st) lptr_pcg8_f_c_aer(10,itype,ai_phase)  = p_pcg8_f_c_a10
            if (p_pcg9_f_c_a10 .ge. p1st) lptr_pcg9_f_c_aer(10,itype,ai_phase)  = p_pcg9_f_c_a10
            if (p_pcg1_f_o_a10 .ge. p1st) lptr_pcg1_f_o_aer(10,itype,ai_phase) = p_pcg1_f_o_a10
            if (p_pcg2_f_o_a10 .ge. p1st) lptr_pcg2_f_o_aer(10,itype,ai_phase)  = p_pcg2_f_o_a10
            if (p_pcg3_f_o_a10 .ge. p1st) lptr_pcg3_f_o_aer(10,itype,ai_phase)  = p_pcg3_f_o_a10
            if (p_pcg4_f_o_a10 .ge. p1st) lptr_pcg4_f_o_aer(10,itype,ai_phase)  = p_pcg4_f_o_a10
            if (p_pcg5_f_o_a10 .ge. p1st) lptr_pcg5_f_o_aer(10,itype,ai_phase)  = p_pcg5_f_o_a10
            if (p_pcg6_f_o_a10 .ge. p1st) lptr_pcg6_f_o_aer(10,itype,ai_phase)  = p_pcg6_f_o_a10
            if (p_pcg7_f_o_a10 .ge. p1st) lptr_pcg7_f_o_aer(10,itype,ai_phase)  = p_pcg7_f_o_a10
            if (p_pcg8_f_o_a10 .ge. p1st) lptr_pcg8_f_o_aer(10,itype,ai_phase)  = p_pcg8_f_o_a10
            if (p_pcg9_f_o_a10 .ge. p1st) lptr_pcg9_f_o_aer(10,itype,ai_phase)  = p_pcg9_f_o_a10
            if (p_opcg1_f_c_a10 .ge. p1st) lptr_opcg1_f_c_aer(10,itype,ai_phase) = p_opcg1_f_c_a10
            if (p_opcg2_f_c_a10 .ge. p1st) lptr_opcg2_f_c_aer(10,itype,ai_phase)  = p_opcg2_f_c_a10
            if (p_opcg3_f_c_a10 .ge. p1st) lptr_opcg3_f_c_aer(10,itype,ai_phase)  = p_opcg3_f_c_a10
            if (p_opcg4_f_c_a10 .ge. p1st) lptr_opcg4_f_c_aer(10,itype,ai_phase)  = p_opcg4_f_c_a10
            if (p_opcg5_f_c_a10 .ge. p1st) lptr_opcg5_f_c_aer(10,itype,ai_phase)  = p_opcg5_f_c_a10
            if (p_opcg6_f_c_a10 .ge. p1st) lptr_opcg6_f_c_aer(10,itype,ai_phase)  = p_opcg6_f_c_a10
            if (p_opcg7_f_c_a10 .ge. p1st) lptr_opcg7_f_c_aer(10,itype,ai_phase)  = p_opcg7_f_c_a10
            if (p_opcg8_f_c_a10 .ge. p1st) lptr_opcg8_f_c_aer(10,itype,ai_phase)  = p_opcg8_f_c_a10
            if (p_opcg1_f_o_a10 .ge. p1st) lptr_opcg1_f_o_aer(10,itype,ai_phase) = p_opcg1_f_o_a10
            if (p_opcg2_f_o_a10 .ge. p1st) lptr_opcg2_f_o_aer(10,itype,ai_phase)  = p_opcg2_f_o_a10
            if (p_opcg3_f_o_a10 .ge. p1st) lptr_opcg3_f_o_aer(10,itype,ai_phase)  = p_opcg3_f_o_a10
            if (p_opcg4_f_o_a10 .ge. p1st) lptr_opcg4_f_o_aer(10,itype,ai_phase)  = p_opcg4_f_o_a10
            if (p_opcg5_f_o_a10 .ge. p1st) lptr_opcg5_f_o_aer(10,itype,ai_phase)  = p_opcg5_f_o_a10
            if (p_opcg6_f_o_a10 .ge. p1st) lptr_opcg6_f_o_aer(10,itype,ai_phase)  = p_opcg6_f_o_a10
            if (p_opcg7_f_o_a10 .ge. p1st) lptr_opcg7_f_o_aer(10,itype,ai_phase)  = p_opcg7_f_o_a10
            if (p_opcg8_f_o_a10 .ge. p1st) lptr_opcg8_f_o_aer(10,itype,ai_phase)  = p_opcg8_f_o_a10
            if (p_ant1_c_a10 .ge. p1st) lptr_ant1_c_aer(10,itype,ai_phase)  = p_ant1_c_a10
            if (p_ant2_c_a10 .ge. p1st) lptr_ant2_c_aer(10,itype,ai_phase)  = p_ant2_c_a10
            if (p_ant3_c_a10 .ge. p1st) lptr_ant3_c_aer(10,itype,ai_phase)  = p_ant3_c_a10
            if (p_ant4_c_a10 .ge. p1st) lptr_ant4_c_aer(10,itype,ai_phase)  = p_ant4_c_a10
            if (p_biog1_c_a10 .ge. p1st) lptr_biog1_c_aer(10,itype,ai_phase)  = p_biog1_c_a10
            if (p_biog2_c_a10 .ge. p1st) lptr_biog2_c_aer(10,itype,ai_phase)  = p_biog2_c_a10
            if (p_biog3_c_a10 .ge. p1st) lptr_biog3_c_aer(10,itype,ai_phase)  = p_biog3_c_a10
            if (p_biog4_c_a10 .ge. p1st) lptr_biog4_c_aer(10,itype,ai_phase)  = p_biog4_c_a10
            if (p_ant1_o_a10 .ge. p1st) lptr_ant1_o_aer(10,itype,ai_phase)  = p_ant1_o_a10
            if (p_ant2_o_a10 .ge. p1st) lptr_ant2_o_aer(10,itype,ai_phase)  = p_ant2_o_a10
            if (p_ant3_o_a10 .ge. p1st) lptr_ant3_o_aer(10,itype,ai_phase)  = p_ant3_o_a10
            if (p_ant4_o_a10 .ge. p1st) lptr_ant4_o_aer(10,itype,ai_phase)  = p_ant4_o_a10
            if (p_biog1_o_a10 .ge. p1st) lptr_biog1_o_aer(10,itype,ai_phase)  = p_biog1_o_a10
            if (p_biog2_o_a10 .ge. p1st) lptr_biog2_o_aer(10,itype,ai_phase)  = p_biog2_o_a10
            if (p_biog3_o_a10 .ge. p1st) lptr_biog3_o_aer(10,itype,ai_phase)  = p_biog3_o_a10
            if (p_biog4_o_a10 .ge. p1st) lptr_biog4_o_aer(10,itype,ai_phase)  = p_biog4_o_a10
            
	    numptr_aer(10,itype,ai_phase)        = p_num_a10
	end if

	if (nsize_aer(itype) .ge. 11) then
	    lptr_so4_aer(11,itype,ai_phase)      = p_so4_a11
	    lptr_no3_aer(11,itype,ai_phase)      = p_no3_a11
	    lptr_cl_aer(11,itype,ai_phase)       = p_cl_a11
	    lptr_msa_aer(11,itype,ai_phase)      = p_msa_a11
	    lptr_co3_aer(11,itype,ai_phase)      = p_co3_a11
	    lptr_nh4_aer(11,itype,ai_phase)      = p_nh4_a11
	    lptr_na_aer(11,itype,ai_phase)       = p_na_a11
	    lptr_ca_aer(11,itype,ai_phase)       = p_ca_a11
	    lptr_oin_aer(11,itype,ai_phase)      = p_oin_a11
	    lptr_oc_aer(11,itype,ai_phase)       = p_oc_a11
	    lptr_bc_aer(11,itype,ai_phase)       = p_bc_a11
	    hyswptr_aer(11,itype)                = p_hysw_a11
	    waterptr_aer(11,itype)               = p_water_a11

            if (p_pcg1_b_c_a11 .ge. p1st) lptr_pcg1_b_c_aer(11,itype,ai_phase) = p_pcg1_b_c_a11
            if (p_pcg2_b_c_a11 .ge. p1st) lptr_pcg2_b_c_aer(11,itype,ai_phase)  = p_pcg2_b_c_a11
            if (p_pcg3_b_c_a11 .ge. p1st) lptr_pcg3_b_c_aer(11,itype,ai_phase)  = p_pcg3_b_c_a11
            if (p_pcg4_b_c_a11 .ge. p1st) lptr_pcg4_b_c_aer(11,itype,ai_phase)  = p_pcg4_b_c_a11
            if (p_pcg5_b_c_a11 .ge. p1st) lptr_pcg5_b_c_aer(11,itype,ai_phase)  = p_pcg5_b_c_a11
            if (p_pcg6_b_c_a11 .ge. p1st) lptr_pcg6_b_c_aer(11,itype,ai_phase)  = p_pcg6_b_c_a11
            if (p_pcg7_b_c_a11 .ge. p1st) lptr_pcg7_b_c_aer(11,itype,ai_phase)  = p_pcg7_b_c_a11
            if (p_pcg8_b_c_a11 .ge. p1st) lptr_pcg8_b_c_aer(11,itype,ai_phase)  = p_pcg8_b_c_a11
            if (p_pcg9_b_c_a11 .ge. p1st) lptr_pcg9_b_c_aer(11,itype,ai_phase)  = p_pcg9_b_c_a11
            if (p_pcg1_b_o_a11 .ge. p1st) lptr_pcg1_b_o_aer(11,itype,ai_phase) = p_pcg1_b_o_a11
            if (p_pcg2_b_o_a11 .ge. p1st) lptr_pcg2_b_o_aer(11,itype,ai_phase)  = p_pcg2_b_o_a11
            if (p_pcg3_b_o_a11 .ge. p1st) lptr_pcg3_b_o_aer(11,itype,ai_phase)  = p_pcg3_b_o_a11
            if (p_pcg4_b_o_a11 .ge. p1st) lptr_pcg4_b_o_aer(11,itype,ai_phase)  = p_pcg4_b_o_a11
            if (p_pcg5_b_o_a11 .ge. p1st) lptr_pcg5_b_o_aer(11,itype,ai_phase)  = p_pcg5_b_o_a11
            if (p_pcg6_b_o_a11 .ge. p1st) lptr_pcg6_b_o_aer(11,itype,ai_phase)  = p_pcg6_b_o_a11
            if (p_pcg7_b_o_a11 .ge. p1st) lptr_pcg7_b_o_aer(11,itype,ai_phase)  = p_pcg7_b_o_a11
            if (p_pcg8_b_o_a11 .ge. p1st) lptr_pcg8_b_o_aer(11,itype,ai_phase)  = p_pcg8_b_o_a11
            if (p_pcg9_b_o_a11 .ge. p1st) lptr_pcg9_b_o_aer(11,itype,ai_phase)  = p_pcg9_b_o_a11
            if (p_opcg1_b_c_a11 .ge. p1st) lptr_opcg1_b_c_aer(11,itype,ai_phase) = p_opcg1_b_c_a11
            if (p_opcg2_b_c_a11 .ge. p1st) lptr_opcg2_b_c_aer(11,itype,ai_phase)  = p_opcg2_b_c_a11
            if (p_opcg3_b_c_a11 .ge. p1st) lptr_opcg3_b_c_aer(11,itype,ai_phase)  = p_opcg3_b_c_a11
            if (p_opcg4_b_c_a11 .ge. p1st) lptr_opcg4_b_c_aer(11,itype,ai_phase)  = p_opcg4_b_c_a11
            if (p_opcg5_b_c_a11 .ge. p1st) lptr_opcg5_b_c_aer(11,itype,ai_phase)  = p_opcg5_b_c_a11
            if (p_opcg6_b_c_a11 .ge. p1st) lptr_opcg6_b_c_aer(11,itype,ai_phase)  = p_opcg6_b_c_a11
            if (p_opcg7_b_c_a11 .ge. p1st) lptr_opcg7_b_c_aer(11,itype,ai_phase)  = p_opcg7_b_c_a11
            if (p_opcg8_b_c_a11 .ge. p1st) lptr_opcg8_b_c_aer(11,itype,ai_phase)  = p_opcg8_b_c_a11
            if (p_opcg1_b_o_a11 .ge. p1st) lptr_opcg1_b_o_aer(11,itype,ai_phase) = p_opcg1_b_o_a11
            if (p_opcg2_b_o_a11 .ge. p1st) lptr_opcg2_b_o_aer(11,itype,ai_phase)  = p_opcg2_b_o_a11
            if (p_opcg3_b_o_a11 .ge. p1st) lptr_opcg3_b_o_aer(11,itype,ai_phase)  = p_opcg3_b_o_a11
            if (p_opcg4_b_o_a11 .ge. p1st) lptr_opcg4_b_o_aer(11,itype,ai_phase)  = p_opcg4_b_o_a11
            if (p_opcg5_b_o_a11 .ge. p1st) lptr_opcg5_b_o_aer(11,itype,ai_phase)  = p_opcg5_b_o_a11
            if (p_opcg6_b_o_a11 .ge. p1st) lptr_opcg6_b_o_aer(11,itype,ai_phase)  = p_opcg6_b_o_a11
            if (p_opcg7_b_o_a11 .ge. p1st) lptr_opcg7_b_o_aer(11,itype,ai_phase)  = p_opcg7_b_o_a11
            if (p_opcg8_b_o_a11 .ge. p1st) lptr_opcg8_b_o_aer(11,itype,ai_phase)  = p_opcg8_b_o_a11
            if (p_pcg1_f_c_a11 .ge. p1st) lptr_pcg1_f_c_aer(11,itype,ai_phase) = p_pcg1_f_c_a11
            if (p_pcg2_f_c_a11 .ge. p1st) lptr_pcg2_f_c_aer(11,itype,ai_phase)  = p_pcg2_f_c_a11
            if (p_pcg3_f_c_a11 .ge. p1st) lptr_pcg3_f_c_aer(11,itype,ai_phase)  = p_pcg3_f_c_a11
            if (p_pcg4_f_c_a11 .ge. p1st) lptr_pcg4_f_c_aer(11,itype,ai_phase)  = p_pcg4_f_c_a11
            if (p_pcg5_f_c_a11 .ge. p1st) lptr_pcg5_f_c_aer(11,itype,ai_phase)  = p_pcg5_f_c_a11
            if (p_pcg6_f_c_a11 .ge. p1st) lptr_pcg6_f_c_aer(11,itype,ai_phase)  = p_pcg6_f_c_a11
            if (p_pcg7_f_c_a11 .ge. p1st) lptr_pcg7_f_c_aer(11,itype,ai_phase)  = p_pcg7_f_c_a11
            if (p_pcg8_f_c_a11 .ge. p1st) lptr_pcg8_f_c_aer(11,itype,ai_phase)  = p_pcg8_f_c_a11
            if (p_pcg9_f_c_a11 .ge. p1st) lptr_pcg9_f_c_aer(11,itype,ai_phase)  = p_pcg9_f_c_a11
            if (p_pcg1_f_o_a11 .ge. p1st) lptr_pcg1_f_o_aer(11,itype,ai_phase) = p_pcg1_f_o_a11
            if (p_pcg2_f_o_a11 .ge. p1st) lptr_pcg2_f_o_aer(11,itype,ai_phase)  = p_pcg2_f_o_a11
            if (p_pcg3_f_o_a11 .ge. p1st) lptr_pcg3_f_o_aer(11,itype,ai_phase)  = p_pcg3_f_o_a11
            if (p_pcg4_f_o_a11 .ge. p1st) lptr_pcg4_f_o_aer(11,itype,ai_phase)  = p_pcg4_f_o_a11
            if (p_pcg5_f_o_a11 .ge. p1st) lptr_pcg5_f_o_aer(11,itype,ai_phase)  = p_pcg5_f_o_a11
            if (p_pcg6_f_o_a11 .ge. p1st) lptr_pcg6_f_o_aer(11,itype,ai_phase)  = p_pcg6_f_o_a11
            if (p_pcg7_f_o_a11 .ge. p1st) lptr_pcg7_f_o_aer(11,itype,ai_phase)  = p_pcg7_f_o_a11
            if (p_pcg8_f_o_a11 .ge. p1st) lptr_pcg8_f_o_aer(11,itype,ai_phase)  = p_pcg8_f_o_a11
            if (p_pcg9_f_o_a11 .ge. p1st) lptr_pcg9_f_o_aer(11,itype,ai_phase)  = p_pcg9_f_o_a11
            if (p_opcg1_f_c_a11 .ge. p1st) lptr_opcg1_f_c_aer(11,itype,ai_phase) = p_opcg1_f_c_a11
            if (p_opcg2_f_c_a11 .ge. p1st) lptr_opcg2_f_c_aer(11,itype,ai_phase)  = p_opcg2_f_c_a11
            if (p_opcg3_f_c_a11 .ge. p1st) lptr_opcg3_f_c_aer(11,itype,ai_phase)  = p_opcg3_f_c_a11
            if (p_opcg4_f_c_a11 .ge. p1st) lptr_opcg4_f_c_aer(11,itype,ai_phase)  = p_opcg4_f_c_a11
            if (p_opcg5_f_c_a11 .ge. p1st) lptr_opcg5_f_c_aer(11,itype,ai_phase)  = p_opcg5_f_c_a11
            if (p_opcg6_f_c_a11 .ge. p1st) lptr_opcg6_f_c_aer(11,itype,ai_phase)  = p_opcg6_f_c_a11
            if (p_opcg7_f_c_a11 .ge. p1st) lptr_opcg7_f_c_aer(11,itype,ai_phase)  = p_opcg7_f_c_a11
            if (p_opcg8_f_c_a11 .ge. p1st) lptr_opcg8_f_c_aer(11,itype,ai_phase)  = p_opcg8_f_c_a11
            if (p_opcg1_f_o_a11 .ge. p1st) lptr_opcg1_f_o_aer(11,itype,ai_phase) = p_opcg1_f_o_a11
            if (p_opcg2_f_o_a11 .ge. p1st) lptr_opcg2_f_o_aer(11,itype,ai_phase)  = p_opcg2_f_o_a11
            if (p_opcg3_f_o_a11 .ge. p1st) lptr_opcg3_f_o_aer(11,itype,ai_phase)  = p_opcg3_f_o_a11
            if (p_opcg4_f_o_a11 .ge. p1st) lptr_opcg4_f_o_aer(11,itype,ai_phase)  = p_opcg4_f_o_a11
            if (p_opcg5_f_o_a11 .ge. p1st) lptr_opcg5_f_o_aer(11,itype,ai_phase)  = p_opcg5_f_o_a11
            if (p_opcg6_f_o_a11 .ge. p1st) lptr_opcg6_f_o_aer(11,itype,ai_phase)  = p_opcg6_f_o_a11
            if (p_opcg7_f_o_a11 .ge. p1st) lptr_opcg7_f_o_aer(11,itype,ai_phase)  = p_opcg7_f_o_a11
            if (p_opcg8_f_o_a11 .ge. p1st) lptr_opcg8_f_o_aer(11,itype,ai_phase)  = p_opcg8_f_o_a11
            if (p_ant1_c_a11 .ge. p1st) lptr_ant1_c_aer(11,itype,ai_phase)  = p_ant1_c_a11
            if (p_ant2_c_a11 .ge. p1st) lptr_ant2_c_aer(11,itype,ai_phase)  = p_ant2_c_a11
            if (p_ant3_c_a11 .ge. p1st) lptr_ant3_c_aer(11,itype,ai_phase)  = p_ant3_c_a11
            if (p_ant4_c_a11 .ge. p1st) lptr_ant4_c_aer(11,itype,ai_phase)  = p_ant4_c_a11
            if (p_biog1_c_a11 .ge. p1st) lptr_biog1_c_aer(11,itype,ai_phase)  = p_biog1_c_a11
            if (p_biog2_c_a11 .ge. p1st) lptr_biog2_c_aer(11,itype,ai_phase)  = p_biog2_c_a11
            if (p_biog3_c_a11 .ge. p1st) lptr_biog3_c_aer(11,itype,ai_phase)  = p_biog3_c_a11
            if (p_biog4_c_a11 .ge. p1st) lptr_biog4_c_aer(11,itype,ai_phase)  = p_biog4_c_a11
            if (p_ant1_o_a11 .ge. p1st) lptr_ant1_o_aer(11,itype,ai_phase)  = p_ant1_o_a11
            if (p_ant2_o_a11 .ge. p1st) lptr_ant2_o_aer(11,itype,ai_phase)  = p_ant2_o_a11
            if (p_ant3_o_a11 .ge. p1st) lptr_ant3_o_aer(11,itype,ai_phase)  = p_ant3_o_a11
            if (p_ant4_o_a11 .ge. p1st) lptr_ant4_o_aer(11,itype,ai_phase)  = p_ant4_o_a11
            if (p_biog1_o_a11 .ge. p1st) lptr_biog1_o_aer(11,itype,ai_phase)  = p_biog1_o_a11
            if (p_biog2_o_a11 .ge. p1st) lptr_biog2_o_aer(11,itype,ai_phase)  = p_biog2_o_a11
            if (p_biog3_o_a11 .ge. p1st) lptr_biog3_o_aer(11,itype,ai_phase)  = p_biog3_o_a11
            if (p_biog4_o_a11 .ge. p1st) lptr_biog4_o_aer(11,itype,ai_phase)  = p_biog4_o_a11
            
	    numptr_aer(11,itype,ai_phase)        = p_num_a11
	end if

	if (nsize_aer(itype) .ge. 12) then
	    lptr_so4_aer(12,itype,ai_phase)      = p_so4_a12
	    lptr_no3_aer(12,itype,ai_phase)      = p_no3_a12
	    lptr_cl_aer(12,itype,ai_phase)       = p_cl_a12
	    lptr_msa_aer(12,itype,ai_phase)      = p_msa_a12
	    lptr_co3_aer(12,itype,ai_phase)      = p_co3_a12
	    lptr_nh4_aer(12,itype,ai_phase)      = p_nh4_a12
	    lptr_na_aer(12,itype,ai_phase)       = p_na_a12
	    lptr_ca_aer(12,itype,ai_phase)       = p_ca_a12
	    lptr_oin_aer(12,itype,ai_phase)      = p_oin_a12
	    lptr_oc_aer(12,itype,ai_phase)       = p_oc_a12
	    lptr_bc_aer(12,itype,ai_phase)       = p_bc_a12
	    hyswptr_aer(12,itype)                = p_hysw_a12
	    waterptr_aer(12,itype)               = p_water_a12

            if (p_pcg1_b_c_a12 .ge. p1st) lptr_pcg1_b_c_aer(12,itype,ai_phase) = p_pcg1_b_c_a12
            if (p_pcg2_b_c_a12 .ge. p1st) lptr_pcg2_b_c_aer(12,itype,ai_phase)  = p_pcg2_b_c_a12
            if (p_pcg3_b_c_a12 .ge. p1st) lptr_pcg3_b_c_aer(12,itype,ai_phase)  = p_pcg3_b_c_a12
            if (p_pcg4_b_c_a12 .ge. p1st) lptr_pcg4_b_c_aer(12,itype,ai_phase)  = p_pcg4_b_c_a12
            if (p_pcg5_b_c_a12 .ge. p1st) lptr_pcg5_b_c_aer(12,itype,ai_phase)  = p_pcg5_b_c_a12
            if (p_pcg6_b_c_a12 .ge. p1st) lptr_pcg6_b_c_aer(12,itype,ai_phase)  = p_pcg6_b_c_a12
            if (p_pcg7_b_c_a12 .ge. p1st) lptr_pcg7_b_c_aer(12,itype,ai_phase)  = p_pcg7_b_c_a12
            if (p_pcg8_b_c_a12 .ge. p1st) lptr_pcg8_b_c_aer(12,itype,ai_phase)  = p_pcg8_b_c_a12
            if (p_pcg9_b_c_a12 .ge. p1st) lptr_pcg9_b_c_aer(12,itype,ai_phase)  = p_pcg9_b_c_a12
            if (p_pcg1_b_o_a12 .ge. p1st) lptr_pcg1_b_o_aer(12,itype,ai_phase) = p_pcg1_b_o_a12
            if (p_pcg2_b_o_a12 .ge. p1st) lptr_pcg2_b_o_aer(12,itype,ai_phase)  = p_pcg2_b_o_a12
            if (p_pcg3_b_o_a12 .ge. p1st) lptr_pcg3_b_o_aer(12,itype,ai_phase)  = p_pcg3_b_o_a12
            if (p_pcg4_b_o_a12 .ge. p1st) lptr_pcg4_b_o_aer(12,itype,ai_phase)  = p_pcg4_b_o_a12
            if (p_pcg5_b_o_a12 .ge. p1st) lptr_pcg5_b_o_aer(12,itype,ai_phase)  = p_pcg5_b_o_a12
            if (p_pcg6_b_o_a12 .ge. p1st) lptr_pcg6_b_o_aer(12,itype,ai_phase)  = p_pcg6_b_o_a12
            if (p_pcg7_b_o_a12 .ge. p1st) lptr_pcg7_b_o_aer(12,itype,ai_phase)  = p_pcg7_b_o_a12
            if (p_pcg8_b_o_a12 .ge. p1st) lptr_pcg8_b_o_aer(12,itype,ai_phase)  = p_pcg8_b_o_a12
            if (p_pcg9_b_o_a12 .ge. p1st) lptr_pcg9_b_o_aer(12,itype,ai_phase)  = p_pcg9_b_o_a12
            if (p_opcg1_b_c_a12 .ge. p1st) lptr_opcg1_b_c_aer(12,itype,ai_phase) = p_opcg1_b_c_a12
            if (p_opcg2_b_c_a12 .ge. p1st) lptr_opcg2_b_c_aer(12,itype,ai_phase)  = p_opcg2_b_c_a12
            if (p_opcg3_b_c_a12 .ge. p1st) lptr_opcg3_b_c_aer(12,itype,ai_phase)  = p_opcg3_b_c_a12
            if (p_opcg4_b_c_a12 .ge. p1st) lptr_opcg4_b_c_aer(12,itype,ai_phase)  = p_opcg4_b_c_a12
            if (p_opcg5_b_c_a12 .ge. p1st) lptr_opcg5_b_c_aer(12,itype,ai_phase)  = p_opcg5_b_c_a12
            if (p_opcg6_b_c_a12 .ge. p1st) lptr_opcg6_b_c_aer(12,itype,ai_phase)  = p_opcg6_b_c_a12
            if (p_opcg7_b_c_a12 .ge. p1st) lptr_opcg7_b_c_aer(12,itype,ai_phase)  = p_opcg7_b_c_a12
            if (p_opcg8_b_c_a12 .ge. p1st) lptr_opcg8_b_c_aer(12,itype,ai_phase)  = p_opcg8_b_c_a12
            if (p_opcg1_b_o_a12 .ge. p1st) lptr_opcg1_b_o_aer(12,itype,ai_phase) = p_opcg1_b_o_a12
            if (p_opcg2_b_o_a12 .ge. p1st) lptr_opcg2_b_o_aer(12,itype,ai_phase)  = p_opcg2_b_o_a12
            if (p_opcg3_b_o_a12 .ge. p1st) lptr_opcg3_b_o_aer(12,itype,ai_phase)  = p_opcg3_b_o_a12
            if (p_opcg4_b_o_a12 .ge. p1st) lptr_opcg4_b_o_aer(12,itype,ai_phase)  = p_opcg4_b_o_a12
            if (p_opcg5_b_o_a12 .ge. p1st) lptr_opcg5_b_o_aer(12,itype,ai_phase)  = p_opcg5_b_o_a12
            if (p_opcg6_b_o_a12 .ge. p1st) lptr_opcg6_b_o_aer(12,itype,ai_phase)  = p_opcg6_b_o_a12
            if (p_opcg7_b_o_a12 .ge. p1st) lptr_opcg7_b_o_aer(12,itype,ai_phase)  = p_opcg7_b_o_a12
            if (p_opcg8_b_o_a12 .ge. p1st) lptr_opcg8_b_o_aer(12,itype,ai_phase)  = p_opcg8_b_o_a12
            if (p_pcg1_f_c_a12 .ge. p1st) lptr_pcg1_f_c_aer(12,itype,ai_phase) = p_pcg1_f_c_a12
            if (p_pcg2_f_c_a12 .ge. p1st) lptr_pcg2_f_c_aer(12,itype,ai_phase)  = p_pcg2_f_c_a12
            if (p_pcg3_f_c_a12 .ge. p1st) lptr_pcg3_f_c_aer(12,itype,ai_phase)  = p_pcg3_f_c_a12
            if (p_pcg4_f_c_a12 .ge. p1st) lptr_pcg4_f_c_aer(12,itype,ai_phase)  = p_pcg4_f_c_a12
            if (p_pcg5_f_c_a12 .ge. p1st) lptr_pcg5_f_c_aer(12,itype,ai_phase)  = p_pcg5_f_c_a12
            if (p_pcg6_f_c_a12 .ge. p1st) lptr_pcg6_f_c_aer(12,itype,ai_phase)  = p_pcg6_f_c_a12
            if (p_pcg7_f_c_a12 .ge. p1st) lptr_pcg7_f_c_aer(12,itype,ai_phase)  = p_pcg7_f_c_a12
            if (p_pcg8_f_c_a12 .ge. p1st) lptr_pcg8_f_c_aer(12,itype,ai_phase)  = p_pcg8_f_c_a12
            if (p_pcg9_f_c_a12 .ge. p1st) lptr_pcg9_f_c_aer(12,itype,ai_phase)  = p_pcg9_f_c_a12
            if (p_pcg1_f_o_a12 .ge. p1st) lptr_pcg1_f_o_aer(12,itype,ai_phase) = p_pcg1_f_o_a12
            if (p_pcg2_f_o_a12 .ge. p1st) lptr_pcg2_f_o_aer(12,itype,ai_phase)  = p_pcg2_f_o_a12
            if (p_pcg3_f_o_a12 .ge. p1st) lptr_pcg3_f_o_aer(12,itype,ai_phase)  = p_pcg3_f_o_a12
            if (p_pcg4_f_o_a12 .ge. p1st) lptr_pcg4_f_o_aer(12,itype,ai_phase)  = p_pcg4_f_o_a12
            if (p_pcg5_f_o_a12 .ge. p1st) lptr_pcg5_f_o_aer(12,itype,ai_phase)  = p_pcg5_f_o_a12
            if (p_pcg6_f_o_a12 .ge. p1st) lptr_pcg6_f_o_aer(12,itype,ai_phase)  = p_pcg6_f_o_a12
            if (p_pcg7_f_o_a12 .ge. p1st) lptr_pcg7_f_o_aer(12,itype,ai_phase)  = p_pcg7_f_o_a12
            if (p_pcg8_f_o_a12 .ge. p1st) lptr_pcg8_f_o_aer(12,itype,ai_phase)  = p_pcg8_f_o_a12
            if (p_pcg9_f_o_a12 .ge. p1st) lptr_pcg9_f_o_aer(12,itype,ai_phase)  = p_pcg9_f_o_a12
            if (p_opcg1_f_c_a12 .ge. p1st) lptr_opcg1_f_c_aer(12,itype,ai_phase) = p_opcg1_f_c_a12
            if (p_opcg2_f_c_a12 .ge. p1st) lptr_opcg2_f_c_aer(12,itype,ai_phase)  = p_opcg2_f_c_a12
            if (p_opcg3_f_c_a12 .ge. p1st) lptr_opcg3_f_c_aer(12,itype,ai_phase)  = p_opcg3_f_c_a12
            if (p_opcg4_f_c_a12 .ge. p1st) lptr_opcg4_f_c_aer(12,itype,ai_phase)  = p_opcg4_f_c_a12
            if (p_opcg5_f_c_a12 .ge. p1st) lptr_opcg5_f_c_aer(12,itype,ai_phase)  = p_opcg5_f_c_a12
            if (p_opcg6_f_c_a12 .ge. p1st) lptr_opcg6_f_c_aer(12,itype,ai_phase)  = p_opcg6_f_c_a12
            if (p_opcg7_f_c_a12 .ge. p1st) lptr_opcg7_f_c_aer(12,itype,ai_phase)  = p_opcg7_f_c_a12
            if (p_opcg8_f_c_a12 .ge. p1st) lptr_opcg8_f_c_aer(12,itype,ai_phase)  = p_opcg8_f_c_a12
            if (p_opcg1_f_o_a12 .ge. p1st) lptr_opcg1_f_o_aer(12,itype,ai_phase) = p_opcg1_f_o_a12
            if (p_opcg2_f_o_a12 .ge. p1st) lptr_opcg2_f_o_aer(12,itype,ai_phase)  = p_opcg2_f_o_a12
            if (p_opcg3_f_o_a12 .ge. p1st) lptr_opcg3_f_o_aer(12,itype,ai_phase)  = p_opcg3_f_o_a12
            if (p_opcg4_f_o_a12 .ge. p1st) lptr_opcg4_f_o_aer(12,itype,ai_phase)  = p_opcg4_f_o_a12
            if (p_opcg5_f_o_a12 .ge. p1st) lptr_opcg5_f_o_aer(12,itype,ai_phase)  = p_opcg5_f_o_a12
            if (p_opcg6_f_o_a12 .ge. p1st) lptr_opcg6_f_o_aer(12,itype,ai_phase)  = p_opcg6_f_o_a12
            if (p_opcg7_f_o_a12 .ge. p1st) lptr_opcg7_f_o_aer(12,itype,ai_phase)  = p_opcg7_f_o_a12
            if (p_opcg8_f_o_a12 .ge. p1st) lptr_opcg8_f_o_aer(12,itype,ai_phase)  = p_opcg8_f_o_a12
            if (p_ant1_c_a12 .ge. p1st) lptr_ant1_c_aer(12,itype,ai_phase)  = p_ant1_c_a12
            if (p_ant2_c_a12 .ge. p1st) lptr_ant2_c_aer(12,itype,ai_phase)  = p_ant2_c_a12
            if (p_ant3_c_a12 .ge. p1st) lptr_ant3_c_aer(12,itype,ai_phase)  = p_ant3_c_a12
            if (p_ant4_c_a12 .ge. p1st) lptr_ant4_c_aer(12,itype,ai_phase)  = p_ant4_c_a12
            if (p_biog1_c_a12 .ge. p1st) lptr_biog1_c_aer(12,itype,ai_phase)  = p_biog1_c_a12
            if (p_biog2_c_a12 .ge. p1st) lptr_biog2_c_aer(12,itype,ai_phase)  = p_biog2_c_a12
            if (p_biog3_c_a12 .ge. p1st) lptr_biog3_c_aer(12,itype,ai_phase)  = p_biog3_c_a12
            if (p_biog4_c_a12 .ge. p1st) lptr_biog4_c_aer(12,itype,ai_phase)  = p_biog4_c_a12
            if (p_ant1_o_a12 .ge. p1st) lptr_ant1_o_aer(12,itype,ai_phase)  = p_ant1_o_a12
            if (p_ant2_o_a12 .ge. p1st) lptr_ant2_o_aer(12,itype,ai_phase)  = p_ant2_o_a12
            if (p_ant3_o_a12 .ge. p1st) lptr_ant3_o_aer(12,itype,ai_phase)  = p_ant3_o_a12
            if (p_ant4_o_a12 .ge. p1st) lptr_ant4_o_aer(12,itype,ai_phase)  = p_ant4_o_a12
            if (p_biog1_o_a12 .ge. p1st) lptr_biog1_o_aer(12,itype,ai_phase)  = p_biog1_o_a12
            if (p_biog2_o_a12 .ge. p1st) lptr_biog2_o_aer(12,itype,ai_phase)  = p_biog2_o_a12
            if (p_biog3_o_a12 .ge. p1st) lptr_biog3_o_aer(12,itype,ai_phase)  = p_biog3_o_a12
            if (p_biog4_o_a12 .ge. p1st) lptr_biog4_o_aer(12,itype,ai_phase)  = p_biog4_o_a12
            
	    numptr_aer(12,itype,ai_phase)        = p_num_a12
	end if

	if (nsize_aer(itype) .ge. 13) then
	    lptr_so4_aer(13,itype,ai_phase)      = p_so4_a13
	    lptr_no3_aer(13,itype,ai_phase)      = p_no3_a13
	    lptr_cl_aer(13,itype,ai_phase)       = p_cl_a13
	    lptr_msa_aer(13,itype,ai_phase)      = p_msa_a13
	    lptr_co3_aer(13,itype,ai_phase)      = p_co3_a13
	    lptr_nh4_aer(13,itype,ai_phase)      = p_nh4_a13
	    lptr_na_aer(13,itype,ai_phase)       = p_na_a13
	    lptr_ca_aer(13,itype,ai_phase)       = p_ca_a13
	    lptr_oin_aer(13,itype,ai_phase)      = p_oin_a13
	    lptr_oc_aer(13,itype,ai_phase)       = p_oc_a13
	    lptr_bc_aer(13,itype,ai_phase)       = p_bc_a13
	    hyswptr_aer(13,itype)                = p_hysw_a13
	    waterptr_aer(13,itype)               = p_water_a13

            if (p_pcg1_b_c_a13 .ge. p1st) lptr_pcg1_b_c_aer(13,itype,ai_phase) = p_pcg1_b_c_a13
            if (p_pcg2_b_c_a13 .ge. p1st) lptr_pcg2_b_c_aer(13,itype,ai_phase)  = p_pcg2_b_c_a13
            if (p_pcg3_b_c_a13 .ge. p1st) lptr_pcg3_b_c_aer(13,itype,ai_phase)  = p_pcg3_b_c_a13
            if (p_pcg4_b_c_a13 .ge. p1st) lptr_pcg4_b_c_aer(13,itype,ai_phase)  = p_pcg4_b_c_a13
            if (p_pcg5_b_c_a13 .ge. p1st) lptr_pcg5_b_c_aer(13,itype,ai_phase)  = p_pcg5_b_c_a13
            if (p_pcg6_b_c_a13 .ge. p1st) lptr_pcg6_b_c_aer(13,itype,ai_phase)  = p_pcg6_b_c_a13
            if (p_pcg7_b_c_a13 .ge. p1st) lptr_pcg7_b_c_aer(13,itype,ai_phase)  = p_pcg7_b_c_a13
            if (p_pcg8_b_c_a13 .ge. p1st) lptr_pcg8_b_c_aer(13,itype,ai_phase)  = p_pcg8_b_c_a13
            if (p_pcg9_b_c_a13 .ge. p1st) lptr_pcg9_b_c_aer(13,itype,ai_phase)  = p_pcg9_b_c_a13
            if (p_pcg1_b_o_a13 .ge. p1st) lptr_pcg1_b_o_aer(13,itype,ai_phase) = p_pcg1_b_o_a13
            if (p_pcg2_b_o_a13 .ge. p1st) lptr_pcg2_b_o_aer(13,itype,ai_phase)  = p_pcg2_b_o_a13
            if (p_pcg3_b_o_a13 .ge. p1st) lptr_pcg3_b_o_aer(13,itype,ai_phase)  = p_pcg3_b_o_a13
            if (p_pcg4_b_o_a13 .ge. p1st) lptr_pcg4_b_o_aer(13,itype,ai_phase)  = p_pcg4_b_o_a13
            if (p_pcg5_b_o_a13 .ge. p1st) lptr_pcg5_b_o_aer(13,itype,ai_phase)  = p_pcg5_b_o_a13
            if (p_pcg6_b_o_a13 .ge. p1st) lptr_pcg6_b_o_aer(13,itype,ai_phase)  = p_pcg6_b_o_a13
            if (p_pcg7_b_o_a13 .ge. p1st) lptr_pcg7_b_o_aer(13,itype,ai_phase)  = p_pcg7_b_o_a13
            if (p_pcg8_b_o_a13 .ge. p1st) lptr_pcg8_b_o_aer(13,itype,ai_phase)  = p_pcg8_b_o_a13
            if (p_pcg9_b_o_a13 .ge. p1st) lptr_pcg9_b_o_aer(13,itype,ai_phase)  = p_pcg9_b_o_a13
            if (p_opcg1_b_c_a13 .ge. p1st) lptr_opcg1_b_c_aer(13,itype,ai_phase) = p_opcg1_b_c_a13
            if (p_opcg2_b_c_a13 .ge. p1st) lptr_opcg2_b_c_aer(13,itype,ai_phase)  = p_opcg2_b_c_a13
            if (p_opcg3_b_c_a13 .ge. p1st) lptr_opcg3_b_c_aer(13,itype,ai_phase)  = p_opcg3_b_c_a13
            if (p_opcg4_b_c_a13 .ge. p1st) lptr_opcg4_b_c_aer(13,itype,ai_phase)  = p_opcg4_b_c_a13
            if (p_opcg5_b_c_a13 .ge. p1st) lptr_opcg5_b_c_aer(13,itype,ai_phase)  = p_opcg5_b_c_a13
            if (p_opcg6_b_c_a13 .ge. p1st) lptr_opcg6_b_c_aer(13,itype,ai_phase)  = p_opcg6_b_c_a13
            if (p_opcg7_b_c_a13 .ge. p1st) lptr_opcg7_b_c_aer(13,itype,ai_phase)  = p_opcg7_b_c_a13
            if (p_opcg8_b_c_a13 .ge. p1st) lptr_opcg8_b_c_aer(13,itype,ai_phase)  = p_opcg8_b_c_a13
            if (p_opcg1_b_o_a13 .ge. p1st) lptr_opcg1_b_o_aer(13,itype,ai_phase) = p_opcg1_b_o_a13
            if (p_opcg2_b_o_a13 .ge. p1st) lptr_opcg2_b_o_aer(13,itype,ai_phase)  = p_opcg2_b_o_a13
            if (p_opcg3_b_o_a13 .ge. p1st) lptr_opcg3_b_o_aer(13,itype,ai_phase)  = p_opcg3_b_o_a13
            if (p_opcg4_b_o_a13 .ge. p1st) lptr_opcg4_b_o_aer(13,itype,ai_phase)  = p_opcg4_b_o_a13
            if (p_opcg5_b_o_a13 .ge. p1st) lptr_opcg5_b_o_aer(13,itype,ai_phase)  = p_opcg5_b_o_a13
            if (p_opcg6_b_o_a13 .ge. p1st) lptr_opcg6_b_o_aer(13,itype,ai_phase)  = p_opcg6_b_o_a13
            if (p_opcg7_b_o_a13 .ge. p1st) lptr_opcg7_b_o_aer(13,itype,ai_phase)  = p_opcg7_b_o_a13
            if (p_opcg8_b_o_a13 .ge. p1st) lptr_opcg8_b_o_aer(13,itype,ai_phase)  = p_opcg8_b_o_a13
            if (p_pcg1_f_c_a13 .ge. p1st) lptr_pcg1_f_c_aer(13,itype,ai_phase) = p_pcg1_f_c_a13
            if (p_pcg2_f_c_a13 .ge. p1st) lptr_pcg2_f_c_aer(13,itype,ai_phase)  = p_pcg2_f_c_a13
            if (p_pcg3_f_c_a13 .ge. p1st) lptr_pcg3_f_c_aer(13,itype,ai_phase)  = p_pcg3_f_c_a13
            if (p_pcg4_f_c_a13 .ge. p1st) lptr_pcg4_f_c_aer(13,itype,ai_phase)  = p_pcg4_f_c_a13
            if (p_pcg5_f_c_a13 .ge. p1st) lptr_pcg5_f_c_aer(13,itype,ai_phase)  = p_pcg5_f_c_a13
            if (p_pcg6_f_c_a13 .ge. p1st) lptr_pcg6_f_c_aer(13,itype,ai_phase)  = p_pcg6_f_c_a13
            if (p_pcg7_f_c_a13 .ge. p1st) lptr_pcg7_f_c_aer(13,itype,ai_phase)  = p_pcg7_f_c_a13
            if (p_pcg8_f_c_a13 .ge. p1st) lptr_pcg8_f_c_aer(13,itype,ai_phase)  = p_pcg8_f_c_a13
            if (p_pcg9_f_c_a13 .ge. p1st) lptr_pcg9_f_c_aer(13,itype,ai_phase)  = p_pcg9_f_c_a13
            if (p_pcg1_f_o_a13 .ge. p1st) lptr_pcg1_f_o_aer(13,itype,ai_phase) = p_pcg1_f_o_a13
            if (p_pcg2_f_o_a13 .ge. p1st) lptr_pcg2_f_o_aer(13,itype,ai_phase)  = p_pcg2_f_o_a13
            if (p_pcg3_f_o_a13 .ge. p1st) lptr_pcg3_f_o_aer(13,itype,ai_phase)  = p_pcg3_f_o_a13
            if (p_pcg4_f_o_a13 .ge. p1st) lptr_pcg4_f_o_aer(13,itype,ai_phase)  = p_pcg4_f_o_a13
            if (p_pcg5_f_o_a13 .ge. p1st) lptr_pcg5_f_o_aer(13,itype,ai_phase)  = p_pcg5_f_o_a13
            if (p_pcg6_f_o_a13 .ge. p1st) lptr_pcg6_f_o_aer(13,itype,ai_phase)  = p_pcg6_f_o_a13
            if (p_pcg7_f_o_a13 .ge. p1st) lptr_pcg7_f_o_aer(13,itype,ai_phase)  = p_pcg7_f_o_a13
            if (p_pcg8_f_o_a13 .ge. p1st) lptr_pcg8_f_o_aer(13,itype,ai_phase)  = p_pcg8_f_o_a13
            if (p_pcg9_f_o_a13 .ge. p1st) lptr_pcg9_f_o_aer(13,itype,ai_phase)  = p_pcg9_f_o_a13
            if (p_opcg1_f_c_a13 .ge. p1st) lptr_opcg1_f_c_aer(13,itype,ai_phase) = p_opcg1_f_c_a13
            if (p_opcg2_f_c_a13 .ge. p1st) lptr_opcg2_f_c_aer(13,itype,ai_phase)  = p_opcg2_f_c_a13
            if (p_opcg3_f_c_a13 .ge. p1st) lptr_opcg3_f_c_aer(13,itype,ai_phase)  = p_opcg3_f_c_a13
            if (p_opcg4_f_c_a13 .ge. p1st) lptr_opcg4_f_c_aer(13,itype,ai_phase)  = p_opcg4_f_c_a13
            if (p_opcg5_f_c_a13 .ge. p1st) lptr_opcg5_f_c_aer(13,itype,ai_phase)  = p_opcg5_f_c_a13
            if (p_opcg6_f_c_a13 .ge. p1st) lptr_opcg6_f_c_aer(13,itype,ai_phase)  = p_opcg6_f_c_a13
            if (p_opcg7_f_c_a13 .ge. p1st) lptr_opcg7_f_c_aer(13,itype,ai_phase)  = p_opcg7_f_c_a13
            if (p_opcg8_f_c_a13 .ge. p1st) lptr_opcg8_f_c_aer(13,itype,ai_phase)  = p_opcg8_f_c_a13
            if (p_opcg1_f_o_a13 .ge. p1st) lptr_opcg1_f_o_aer(13,itype,ai_phase) = p_opcg1_f_o_a13
            if (p_opcg2_f_o_a13 .ge. p1st) lptr_opcg2_f_o_aer(13,itype,ai_phase)  = p_opcg2_f_o_a13
            if (p_opcg3_f_o_a13 .ge. p1st) lptr_opcg3_f_o_aer(13,itype,ai_phase)  = p_opcg3_f_o_a13
            if (p_opcg4_f_o_a13 .ge. p1st) lptr_opcg4_f_o_aer(13,itype,ai_phase)  = p_opcg4_f_o_a13
            if (p_opcg5_f_o_a13 .ge. p1st) lptr_opcg5_f_o_aer(13,itype,ai_phase)  = p_opcg5_f_o_a13
            if (p_opcg6_f_o_a13 .ge. p1st) lptr_opcg6_f_o_aer(13,itype,ai_phase)  = p_opcg6_f_o_a13
            if (p_opcg7_f_o_a13 .ge. p1st) lptr_opcg7_f_o_aer(13,itype,ai_phase)  = p_opcg7_f_o_a13
            if (p_opcg8_f_o_a13 .ge. p1st) lptr_opcg8_f_o_aer(13,itype,ai_phase)  = p_opcg8_f_o_a13
            if (p_ant1_c_a13 .ge. p1st) lptr_ant1_c_aer(13,itype,ai_phase)  = p_ant1_c_a13
            if (p_ant2_c_a13 .ge. p1st) lptr_ant2_c_aer(13,itype,ai_phase)  = p_ant2_c_a13
            if (p_ant3_c_a13 .ge. p1st) lptr_ant3_c_aer(13,itype,ai_phase)  = p_ant3_c_a13
            if (p_ant4_c_a13 .ge. p1st) lptr_ant4_c_aer(13,itype,ai_phase)  = p_ant4_c_a13
            if (p_biog1_c_a13 .ge. p1st) lptr_biog1_c_aer(13,itype,ai_phase)  = p_biog1_c_a13
            if (p_biog2_c_a13 .ge. p1st) lptr_biog2_c_aer(13,itype,ai_phase)  = p_biog2_c_a13
            if (p_biog3_c_a13 .ge. p1st) lptr_biog3_c_aer(13,itype,ai_phase)  = p_biog3_c_a13
            if (p_biog4_c_a13 .ge. p1st) lptr_biog4_c_aer(13,itype,ai_phase)  = p_biog4_c_a13
            if (p_ant1_o_a13 .ge. p1st) lptr_ant1_o_aer(13,itype,ai_phase)  = p_ant1_o_a13
            if (p_ant2_o_a13 .ge. p1st) lptr_ant2_o_aer(13,itype,ai_phase)  = p_ant2_o_a13
            if (p_ant3_o_a13 .ge. p1st) lptr_ant3_o_aer(13,itype,ai_phase)  = p_ant3_o_a13
            if (p_ant4_o_a13 .ge. p1st) lptr_ant4_o_aer(13,itype,ai_phase)  = p_ant4_o_a13
            if (p_biog1_o_a13 .ge. p1st) lptr_biog1_o_aer(13,itype,ai_phase)  = p_biog1_o_a13
            if (p_biog2_o_a13 .ge. p1st) lptr_biog2_o_aer(13,itype,ai_phase)  = p_biog2_o_a13
            if (p_biog3_o_a13 .ge. p1st) lptr_biog3_o_aer(13,itype,ai_phase)  = p_biog3_o_a13
            if (p_biog4_o_a13 .ge. p1st) lptr_biog4_o_aer(13,itype,ai_phase)  = p_biog4_o_a13
            
	    numptr_aer(13,itype,ai_phase)        = p_num_a13
	end if

	if (nsize_aer(itype) .ge. 14) then
	    lptr_so4_aer(14,itype,ai_phase)      = p_so4_a14
	    lptr_no3_aer(14,itype,ai_phase)      = p_no3_a14
	    lptr_cl_aer(14,itype,ai_phase)       = p_cl_a14
	    lptr_msa_aer(14,itype,ai_phase)      = p_msa_a14
	    lptr_co3_aer(14,itype,ai_phase)      = p_co3_a14
	    lptr_nh4_aer(14,itype,ai_phase)      = p_nh4_a14
	    lptr_na_aer(14,itype,ai_phase)       = p_na_a14
	    lptr_ca_aer(14,itype,ai_phase)       = p_ca_a14
	    lptr_oin_aer(14,itype,ai_phase)      = p_oin_a14
	    lptr_oc_aer(14,itype,ai_phase)       = p_oc_a14
	    lptr_bc_aer(14,itype,ai_phase)       = p_bc_a14
	    hyswptr_aer(14,itype)                = p_hysw_a14
	    waterptr_aer(14,itype)               = p_water_a14

            if (p_pcg1_b_c_a14 .ge. p1st) lptr_pcg1_b_c_aer(14,itype,ai_phase) = p_pcg1_b_c_a14
            if (p_pcg2_b_c_a14 .ge. p1st) lptr_pcg2_b_c_aer(14,itype,ai_phase)  = p_pcg2_b_c_a14
            if (p_pcg3_b_c_a14 .ge. p1st) lptr_pcg3_b_c_aer(14,itype,ai_phase)  = p_pcg3_b_c_a14
            if (p_pcg4_b_c_a14 .ge. p1st) lptr_pcg4_b_c_aer(14,itype,ai_phase)  = p_pcg4_b_c_a14
            if (p_pcg5_b_c_a14 .ge. p1st) lptr_pcg5_b_c_aer(14,itype,ai_phase)  = p_pcg5_b_c_a14
            if (p_pcg6_b_c_a14 .ge. p1st) lptr_pcg6_b_c_aer(14,itype,ai_phase)  = p_pcg6_b_c_a14
            if (p_pcg7_b_c_a14 .ge. p1st) lptr_pcg7_b_c_aer(14,itype,ai_phase)  = p_pcg7_b_c_a14
            if (p_pcg8_b_c_a14 .ge. p1st) lptr_pcg8_b_c_aer(14,itype,ai_phase)  = p_pcg8_b_c_a14
            if (p_pcg9_b_c_a14 .ge. p1st) lptr_pcg9_b_c_aer(14,itype,ai_phase)  = p_pcg9_b_c_a14
            if (p_pcg1_b_o_a14 .ge. p1st) lptr_pcg1_b_o_aer(14,itype,ai_phase) = p_pcg1_b_o_a14
            if (p_pcg2_b_o_a14 .ge. p1st) lptr_pcg2_b_o_aer(14,itype,ai_phase)  = p_pcg2_b_o_a14
            if (p_pcg3_b_o_a14 .ge. p1st) lptr_pcg3_b_o_aer(14,itype,ai_phase)  = p_pcg3_b_o_a14
            if (p_pcg4_b_o_a14 .ge. p1st) lptr_pcg4_b_o_aer(14,itype,ai_phase)  = p_pcg4_b_o_a14
            if (p_pcg5_b_o_a14 .ge. p1st) lptr_pcg5_b_o_aer(14,itype,ai_phase)  = p_pcg5_b_o_a14
            if (p_pcg6_b_o_a14 .ge. p1st) lptr_pcg6_b_o_aer(14,itype,ai_phase)  = p_pcg6_b_o_a14
            if (p_pcg7_b_o_a14 .ge. p1st) lptr_pcg7_b_o_aer(14,itype,ai_phase)  = p_pcg7_b_o_a14
            if (p_pcg8_b_o_a14 .ge. p1st) lptr_pcg8_b_o_aer(14,itype,ai_phase)  = p_pcg8_b_o_a14
            if (p_pcg9_b_o_a14 .ge. p1st) lptr_pcg9_b_o_aer(14,itype,ai_phase)  = p_pcg9_b_o_a14
            if (p_opcg1_b_c_a14 .ge. p1st) lptr_opcg1_b_c_aer(14,itype,ai_phase) = p_opcg1_b_c_a14
            if (p_opcg2_b_c_a14 .ge. p1st) lptr_opcg2_b_c_aer(14,itype,ai_phase)  = p_opcg2_b_c_a14
            if (p_opcg3_b_c_a14 .ge. p1st) lptr_opcg3_b_c_aer(14,itype,ai_phase)  = p_opcg3_b_c_a14
            if (p_opcg4_b_c_a14 .ge. p1st) lptr_opcg4_b_c_aer(14,itype,ai_phase)  = p_opcg4_b_c_a14
            if (p_opcg5_b_c_a14 .ge. p1st) lptr_opcg5_b_c_aer(14,itype,ai_phase)  = p_opcg5_b_c_a14
            if (p_opcg6_b_c_a14 .ge. p1st) lptr_opcg6_b_c_aer(14,itype,ai_phase)  = p_opcg6_b_c_a14
            if (p_opcg7_b_c_a14 .ge. p1st) lptr_opcg7_b_c_aer(14,itype,ai_phase)  = p_opcg7_b_c_a14
            if (p_opcg8_b_c_a14 .ge. p1st) lptr_opcg8_b_c_aer(14,itype,ai_phase)  = p_opcg8_b_c_a14
            if (p_opcg1_b_o_a14 .ge. p1st) lptr_opcg1_b_o_aer(14,itype,ai_phase) = p_opcg1_b_o_a14
            if (p_opcg2_b_o_a14 .ge. p1st) lptr_opcg2_b_o_aer(14,itype,ai_phase)  = p_opcg2_b_o_a14
            if (p_opcg3_b_o_a14 .ge. p1st) lptr_opcg3_b_o_aer(14,itype,ai_phase)  = p_opcg3_b_o_a14
            if (p_opcg4_b_o_a14 .ge. p1st) lptr_opcg4_b_o_aer(14,itype,ai_phase)  = p_opcg4_b_o_a14
            if (p_opcg5_b_o_a14 .ge. p1st) lptr_opcg5_b_o_aer(14,itype,ai_phase)  = p_opcg5_b_o_a14
            if (p_opcg6_b_o_a14 .ge. p1st) lptr_opcg6_b_o_aer(14,itype,ai_phase)  = p_opcg6_b_o_a14
            if (p_opcg7_b_o_a14 .ge. p1st) lptr_opcg7_b_o_aer(14,itype,ai_phase)  = p_opcg7_b_o_a14
            if (p_opcg8_b_o_a14 .ge. p1st) lptr_opcg8_b_o_aer(14,itype,ai_phase)  = p_opcg8_b_o_a14
            if (p_pcg1_f_c_a14 .ge. p1st) lptr_pcg1_f_c_aer(14,itype,ai_phase) = p_pcg1_f_c_a14
            if (p_pcg2_f_c_a14 .ge. p1st) lptr_pcg2_f_c_aer(14,itype,ai_phase)  = p_pcg2_f_c_a14
            if (p_pcg3_f_c_a14 .ge. p1st) lptr_pcg3_f_c_aer(14,itype,ai_phase)  = p_pcg3_f_c_a14
            if (p_pcg4_f_c_a14 .ge. p1st) lptr_pcg4_f_c_aer(14,itype,ai_phase)  = p_pcg4_f_c_a14
            if (p_pcg5_f_c_a14 .ge. p1st) lptr_pcg5_f_c_aer(14,itype,ai_phase)  = p_pcg5_f_c_a14
            if (p_pcg6_f_c_a14 .ge. p1st) lptr_pcg6_f_c_aer(14,itype,ai_phase)  = p_pcg6_f_c_a14
            if (p_pcg7_f_c_a14 .ge. p1st) lptr_pcg7_f_c_aer(14,itype,ai_phase)  = p_pcg7_f_c_a14
            if (p_pcg8_f_c_a14 .ge. p1st) lptr_pcg8_f_c_aer(14,itype,ai_phase)  = p_pcg8_f_c_a14
            if (p_pcg9_f_c_a14 .ge. p1st) lptr_pcg9_f_c_aer(14,itype,ai_phase)  = p_pcg9_f_c_a14
            if (p_pcg1_f_o_a14 .ge. p1st) lptr_pcg1_f_o_aer(14,itype,ai_phase) = p_pcg1_f_o_a14
            if (p_pcg2_f_o_a14 .ge. p1st) lptr_pcg2_f_o_aer(14,itype,ai_phase)  = p_pcg2_f_o_a14
            if (p_pcg3_f_o_a14 .ge. p1st) lptr_pcg3_f_o_aer(14,itype,ai_phase)  = p_pcg3_f_o_a14
            if (p_pcg4_f_o_a14 .ge. p1st) lptr_pcg4_f_o_aer(14,itype,ai_phase)  = p_pcg4_f_o_a14
            if (p_pcg5_f_o_a14 .ge. p1st) lptr_pcg5_f_o_aer(14,itype,ai_phase)  = p_pcg5_f_o_a14
            if (p_pcg6_f_o_a14 .ge. p1st) lptr_pcg6_f_o_aer(14,itype,ai_phase)  = p_pcg6_f_o_a14
            if (p_pcg7_f_o_a14 .ge. p1st) lptr_pcg7_f_o_aer(14,itype,ai_phase)  = p_pcg7_f_o_a14
            if (p_pcg8_f_o_a14 .ge. p1st) lptr_pcg8_f_o_aer(14,itype,ai_phase)  = p_pcg8_f_o_a14
            if (p_pcg9_f_o_a14 .ge. p1st) lptr_pcg9_f_o_aer(14,itype,ai_phase)  = p_pcg9_f_o_a14
            if (p_opcg1_f_c_a14 .ge. p1st) lptr_opcg1_f_c_aer(14,itype,ai_phase) = p_opcg1_f_c_a14
            if (p_opcg2_f_c_a14 .ge. p1st) lptr_opcg2_f_c_aer(14,itype,ai_phase)  = p_opcg2_f_c_a14
            if (p_opcg3_f_c_a14 .ge. p1st) lptr_opcg3_f_c_aer(14,itype,ai_phase)  = p_opcg3_f_c_a14
            if (p_opcg4_f_c_a14 .ge. p1st) lptr_opcg4_f_c_aer(14,itype,ai_phase)  = p_opcg4_f_c_a14
            if (p_opcg5_f_c_a14 .ge. p1st) lptr_opcg5_f_c_aer(14,itype,ai_phase)  = p_opcg5_f_c_a14
            if (p_opcg6_f_c_a14 .ge. p1st) lptr_opcg6_f_c_aer(14,itype,ai_phase)  = p_opcg6_f_c_a14
            if (p_opcg7_f_c_a14 .ge. p1st) lptr_opcg7_f_c_aer(14,itype,ai_phase)  = p_opcg7_f_c_a14
            if (p_opcg8_f_c_a14 .ge. p1st) lptr_opcg8_f_c_aer(14,itype,ai_phase)  = p_opcg8_f_c_a14
            if (p_opcg1_f_o_a14 .ge. p1st) lptr_opcg1_f_o_aer(14,itype,ai_phase) = p_opcg1_f_o_a14
            if (p_opcg2_f_o_a14 .ge. p1st) lptr_opcg2_f_o_aer(14,itype,ai_phase)  = p_opcg2_f_o_a14
            if (p_opcg3_f_o_a14 .ge. p1st) lptr_opcg3_f_o_aer(14,itype,ai_phase)  = p_opcg3_f_o_a14
            if (p_opcg4_f_o_a14 .ge. p1st) lptr_opcg4_f_o_aer(14,itype,ai_phase)  = p_opcg4_f_o_a14
            if (p_opcg5_f_o_a14 .ge. p1st) lptr_opcg5_f_o_aer(14,itype,ai_phase)  = p_opcg5_f_o_a14
            if (p_opcg6_f_o_a14 .ge. p1st) lptr_opcg6_f_o_aer(14,itype,ai_phase)  = p_opcg6_f_o_a14
            if (p_opcg7_f_o_a14 .ge. p1st) lptr_opcg7_f_o_aer(14,itype,ai_phase)  = p_opcg7_f_o_a14
            if (p_opcg8_f_o_a14 .ge. p1st) lptr_opcg8_f_o_aer(14,itype,ai_phase)  = p_opcg8_f_o_a14
            if (p_ant1_c_a14 .ge. p1st) lptr_ant1_c_aer(14,itype,ai_phase)  = p_ant1_c_a14
            if (p_ant2_c_a14 .ge. p1st) lptr_ant2_c_aer(14,itype,ai_phase)  = p_ant2_c_a14
            if (p_ant3_c_a14 .ge. p1st) lptr_ant3_c_aer(14,itype,ai_phase)  = p_ant3_c_a14
            if (p_ant4_c_a14 .ge. p1st) lptr_ant4_c_aer(14,itype,ai_phase)  = p_ant4_c_a14
            if (p_biog1_c_a14 .ge. p1st) lptr_biog1_c_aer(14,itype,ai_phase)  = p_biog1_c_a14
            if (p_biog2_c_a14 .ge. p1st) lptr_biog2_c_aer(14,itype,ai_phase)  = p_biog2_c_a14
            if (p_biog3_c_a14 .ge. p1st) lptr_biog3_c_aer(14,itype,ai_phase)  = p_biog3_c_a14
            if (p_biog4_c_a14 .ge. p1st) lptr_biog4_c_aer(14,itype,ai_phase)  = p_biog4_c_a14
            if (p_ant1_o_a14 .ge. p1st) lptr_ant1_o_aer(14,itype,ai_phase)  = p_ant1_o_a14
            if (p_ant2_o_a14 .ge. p1st) lptr_ant2_o_aer(14,itype,ai_phase)  = p_ant2_o_a14
            if (p_ant3_o_a14 .ge. p1st) lptr_ant3_o_aer(14,itype,ai_phase)  = p_ant3_o_a14
            if (p_ant4_o_a14 .ge. p1st) lptr_ant4_o_aer(14,itype,ai_phase)  = p_ant4_o_a14
            if (p_biog1_o_a14 .ge. p1st) lptr_biog1_o_aer(14,itype,ai_phase)  = p_biog1_o_a14
            if (p_biog2_o_a14 .ge. p1st) lptr_biog2_o_aer(14,itype,ai_phase)  = p_biog2_o_a14
            if (p_biog3_o_a14 .ge. p1st) lptr_biog3_o_aer(14,itype,ai_phase)  = p_biog3_o_a14
            if (p_biog4_o_a14 .ge. p1st) lptr_biog4_o_aer(14,itype,ai_phase)  = p_biog4_o_a14
            
	    numptr_aer(14,itype,ai_phase)        = p_num_a14
	end if

	if (nsize_aer(itype) .ge. 15) then
	    lptr_so4_aer(15,itype,ai_phase)      = p_so4_a15
	    lptr_no3_aer(15,itype,ai_phase)      = p_no3_a15
	    lptr_cl_aer(15,itype,ai_phase)       = p_cl_a15
	    lptr_msa_aer(15,itype,ai_phase)      = p_msa_a15
	    lptr_co3_aer(15,itype,ai_phase)      = p_co3_a15
	    lptr_nh4_aer(15,itype,ai_phase)      = p_nh4_a15
	    lptr_na_aer(15,itype,ai_phase)       = p_na_a15
	    lptr_ca_aer(15,itype,ai_phase)       = p_ca_a15
	    lptr_oin_aer(15,itype,ai_phase)      = p_oin_a15
	    lptr_oc_aer(15,itype,ai_phase)       = p_oc_a15
	    lptr_bc_aer(15,itype,ai_phase)       = p_bc_a15
	    hyswptr_aer(15,itype)                = p_hysw_a15
	    waterptr_aer(15,itype)               = p_water_a15

            if (p_pcg1_b_c_a15 .ge. p1st) lptr_pcg1_b_c_aer(15,itype,ai_phase) = p_pcg1_b_c_a15
            if (p_pcg2_b_c_a15 .ge. p1st) lptr_pcg2_b_c_aer(15,itype,ai_phase)  = p_pcg2_b_c_a15
            if (p_pcg3_b_c_a15 .ge. p1st) lptr_pcg3_b_c_aer(15,itype,ai_phase)  = p_pcg3_b_c_a15
            if (p_pcg4_b_c_a15 .ge. p1st) lptr_pcg4_b_c_aer(15,itype,ai_phase)  = p_pcg4_b_c_a15
            if (p_pcg5_b_c_a15 .ge. p1st) lptr_pcg5_b_c_aer(15,itype,ai_phase)  = p_pcg5_b_c_a15
            if (p_pcg6_b_c_a15 .ge. p1st) lptr_pcg6_b_c_aer(15,itype,ai_phase)  = p_pcg6_b_c_a15
            if (p_pcg7_b_c_a15 .ge. p1st) lptr_pcg7_b_c_aer(15,itype,ai_phase)  = p_pcg7_b_c_a15
            if (p_pcg8_b_c_a15 .ge. p1st) lptr_pcg8_b_c_aer(15,itype,ai_phase)  = p_pcg8_b_c_a15
            if (p_pcg9_b_c_a15 .ge. p1st) lptr_pcg9_b_c_aer(15,itype,ai_phase)  = p_pcg9_b_c_a15
            if (p_pcg1_b_o_a15 .ge. p1st) lptr_pcg1_b_o_aer(15,itype,ai_phase) = p_pcg1_b_o_a15
            if (p_pcg2_b_o_a15 .ge. p1st) lptr_pcg2_b_o_aer(15,itype,ai_phase)  = p_pcg2_b_o_a15
            if (p_pcg3_b_o_a15 .ge. p1st) lptr_pcg3_b_o_aer(15,itype,ai_phase)  = p_pcg3_b_o_a15
            if (p_pcg4_b_o_a15 .ge. p1st) lptr_pcg4_b_o_aer(15,itype,ai_phase)  = p_pcg4_b_o_a15
            if (p_pcg5_b_o_a15 .ge. p1st) lptr_pcg5_b_o_aer(15,itype,ai_phase)  = p_pcg5_b_o_a15
            if (p_pcg6_b_o_a15 .ge. p1st) lptr_pcg6_b_o_aer(15,itype,ai_phase)  = p_pcg6_b_o_a15
            if (p_pcg7_b_o_a15 .ge. p1st) lptr_pcg7_b_o_aer(15,itype,ai_phase)  = p_pcg7_b_o_a15
            if (p_pcg8_b_o_a15 .ge. p1st) lptr_pcg8_b_o_aer(15,itype,ai_phase)  = p_pcg8_b_o_a15
            if (p_pcg9_b_o_a15 .ge. p1st) lptr_pcg9_b_o_aer(15,itype,ai_phase)  = p_pcg9_b_o_a15
            if (p_opcg1_b_c_a15 .ge. p1st) lptr_opcg1_b_c_aer(15,itype,ai_phase) = p_opcg1_b_c_a15
            if (p_opcg2_b_c_a15 .ge. p1st) lptr_opcg2_b_c_aer(15,itype,ai_phase)  = p_opcg2_b_c_a15
            if (p_opcg3_b_c_a15 .ge. p1st) lptr_opcg3_b_c_aer(15,itype,ai_phase)  = p_opcg3_b_c_a15
            if (p_opcg4_b_c_a15 .ge. p1st) lptr_opcg4_b_c_aer(15,itype,ai_phase)  = p_opcg4_b_c_a15
            if (p_opcg5_b_c_a15 .ge. p1st) lptr_opcg5_b_c_aer(15,itype,ai_phase)  = p_opcg5_b_c_a15
            if (p_opcg6_b_c_a15 .ge. p1st) lptr_opcg6_b_c_aer(15,itype,ai_phase)  = p_opcg6_b_c_a15
            if (p_opcg7_b_c_a15 .ge. p1st) lptr_opcg7_b_c_aer(15,itype,ai_phase)  = p_opcg7_b_c_a15
            if (p_opcg8_b_c_a15 .ge. p1st) lptr_opcg8_b_c_aer(15,itype,ai_phase)  = p_opcg8_b_c_a15
            if (p_opcg1_b_o_a15 .ge. p1st) lptr_opcg1_b_o_aer(15,itype,ai_phase) = p_opcg1_b_o_a15
            if (p_opcg2_b_o_a15 .ge. p1st) lptr_opcg2_b_o_aer(15,itype,ai_phase)  = p_opcg2_b_o_a15
            if (p_opcg3_b_o_a15 .ge. p1st) lptr_opcg3_b_o_aer(15,itype,ai_phase)  = p_opcg3_b_o_a15
            if (p_opcg4_b_o_a15 .ge. p1st) lptr_opcg4_b_o_aer(15,itype,ai_phase)  = p_opcg4_b_o_a15
            if (p_opcg5_b_o_a15 .ge. p1st) lptr_opcg5_b_o_aer(15,itype,ai_phase)  = p_opcg5_b_o_a15
            if (p_opcg6_b_o_a15 .ge. p1st) lptr_opcg6_b_o_aer(15,itype,ai_phase)  = p_opcg6_b_o_a15
            if (p_opcg7_b_o_a15 .ge. p1st) lptr_opcg7_b_o_aer(15,itype,ai_phase)  = p_opcg7_b_o_a15
            if (p_opcg8_b_o_a15 .ge. p1st) lptr_opcg8_b_o_aer(15,itype,ai_phase)  = p_opcg8_b_o_a15
            if (p_pcg1_f_c_a15 .ge. p1st) lptr_pcg1_f_c_aer(15,itype,ai_phase) = p_pcg1_f_c_a15
            if (p_pcg2_f_c_a15 .ge. p1st) lptr_pcg2_f_c_aer(15,itype,ai_phase)  = p_pcg2_f_c_a15
            if (p_pcg3_f_c_a15 .ge. p1st) lptr_pcg3_f_c_aer(15,itype,ai_phase)  = p_pcg3_f_c_a15
            if (p_pcg4_f_c_a15 .ge. p1st) lptr_pcg4_f_c_aer(15,itype,ai_phase)  = p_pcg4_f_c_a15
            if (p_pcg5_f_c_a15 .ge. p1st) lptr_pcg5_f_c_aer(15,itype,ai_phase)  = p_pcg5_f_c_a15
            if (p_pcg6_f_c_a15 .ge. p1st) lptr_pcg6_f_c_aer(15,itype,ai_phase)  = p_pcg6_f_c_a15
            if (p_pcg7_f_c_a15 .ge. p1st) lptr_pcg7_f_c_aer(15,itype,ai_phase)  = p_pcg7_f_c_a15
            if (p_pcg8_f_c_a15 .ge. p1st) lptr_pcg8_f_c_aer(15,itype,ai_phase)  = p_pcg8_f_c_a15
            if (p_pcg9_f_c_a15 .ge. p1st) lptr_pcg9_f_c_aer(15,itype,ai_phase)  = p_pcg9_f_c_a15
            if (p_pcg1_f_o_a15 .ge. p1st) lptr_pcg1_f_o_aer(15,itype,ai_phase) = p_pcg1_f_o_a15
            if (p_pcg2_f_o_a15 .ge. p1st) lptr_pcg2_f_o_aer(15,itype,ai_phase)  = p_pcg2_f_o_a15
            if (p_pcg3_f_o_a15 .ge. p1st) lptr_pcg3_f_o_aer(15,itype,ai_phase)  = p_pcg3_f_o_a15
            if (p_pcg4_f_o_a15 .ge. p1st) lptr_pcg4_f_o_aer(15,itype,ai_phase)  = p_pcg4_f_o_a15
            if (p_pcg5_f_o_a15 .ge. p1st) lptr_pcg5_f_o_aer(15,itype,ai_phase)  = p_pcg5_f_o_a15
            if (p_pcg6_f_o_a15 .ge. p1st) lptr_pcg6_f_o_aer(15,itype,ai_phase)  = p_pcg6_f_o_a15
            if (p_pcg7_f_o_a15 .ge. p1st) lptr_pcg7_f_o_aer(15,itype,ai_phase)  = p_pcg7_f_o_a15
            if (p_pcg8_f_o_a15 .ge. p1st) lptr_pcg8_f_o_aer(15,itype,ai_phase)  = p_pcg8_f_o_a15
            if (p_pcg9_f_o_a15 .ge. p1st) lptr_pcg9_f_o_aer(15,itype,ai_phase)  = p_pcg9_f_o_a15
            if (p_opcg1_f_c_a15 .ge. p1st) lptr_opcg1_f_c_aer(15,itype,ai_phase) = p_opcg1_f_c_a15
            if (p_opcg2_f_c_a15 .ge. p1st) lptr_opcg2_f_c_aer(15,itype,ai_phase)  = p_opcg2_f_c_a15
            if (p_opcg3_f_c_a15 .ge. p1st) lptr_opcg3_f_c_aer(15,itype,ai_phase)  = p_opcg3_f_c_a15
            if (p_opcg4_f_c_a15 .ge. p1st) lptr_opcg4_f_c_aer(15,itype,ai_phase)  = p_opcg4_f_c_a15
            if (p_opcg5_f_c_a15 .ge. p1st) lptr_opcg5_f_c_aer(15,itype,ai_phase)  = p_opcg5_f_c_a15
            if (p_opcg6_f_c_a15 .ge. p1st) lptr_opcg6_f_c_aer(15,itype,ai_phase)  = p_opcg6_f_c_a15
            if (p_opcg7_f_c_a15 .ge. p1st) lptr_opcg7_f_c_aer(15,itype,ai_phase)  = p_opcg7_f_c_a15
            if (p_opcg8_f_c_a15 .ge. p1st) lptr_opcg8_f_c_aer(15,itype,ai_phase)  = p_opcg8_f_c_a15
            if (p_opcg1_f_o_a15 .ge. p1st) lptr_opcg1_f_o_aer(15,itype,ai_phase) = p_opcg1_f_o_a15
            if (p_opcg2_f_o_a15 .ge. p1st) lptr_opcg2_f_o_aer(15,itype,ai_phase)  = p_opcg2_f_o_a15
            if (p_opcg3_f_o_a15 .ge. p1st) lptr_opcg3_f_o_aer(15,itype,ai_phase)  = p_opcg3_f_o_a15
            if (p_opcg4_f_o_a15 .ge. p1st) lptr_opcg4_f_o_aer(15,itype,ai_phase)  = p_opcg4_f_o_a15
            if (p_opcg5_f_o_a15 .ge. p1st) lptr_opcg5_f_o_aer(15,itype,ai_phase)  = p_opcg5_f_o_a15
            if (p_opcg6_f_o_a15 .ge. p1st) lptr_opcg6_f_o_aer(15,itype,ai_phase)  = p_opcg6_f_o_a15
            if (p_opcg7_f_o_a15 .ge. p1st) lptr_opcg7_f_o_aer(15,itype,ai_phase)  = p_opcg7_f_o_a15
            if (p_opcg8_f_o_a15 .ge. p1st) lptr_opcg8_f_o_aer(15,itype,ai_phase)  = p_opcg8_f_o_a15
            if (p_ant1_c_a15 .ge. p1st) lptr_ant1_c_aer(15,itype,ai_phase)  = p_ant1_c_a15
            if (p_ant2_c_a15 .ge. p1st) lptr_ant2_c_aer(15,itype,ai_phase)  = p_ant2_c_a15
            if (p_ant3_c_a15 .ge. p1st) lptr_ant3_c_aer(15,itype,ai_phase)  = p_ant3_c_a15
            if (p_ant4_c_a15 .ge. p1st) lptr_ant4_c_aer(15,itype,ai_phase)  = p_ant4_c_a15
            if (p_biog1_c_a15 .ge. p1st) lptr_biog1_c_aer(15,itype,ai_phase)  = p_biog1_c_a15
            if (p_biog2_c_a15 .ge. p1st) lptr_biog2_c_aer(15,itype,ai_phase)  = p_biog2_c_a15
            if (p_biog3_c_a15 .ge. p1st) lptr_biog3_c_aer(15,itype,ai_phase)  = p_biog3_c_a15
            if (p_biog4_c_a15 .ge. p1st) lptr_biog4_c_aer(15,itype,ai_phase)  = p_biog4_c_a15
            if (p_ant1_o_a15 .ge. p1st) lptr_ant1_o_aer(15,itype,ai_phase)  = p_ant1_o_a15
            if (p_ant2_o_a15 .ge. p1st) lptr_ant2_o_aer(15,itype,ai_phase)  = p_ant2_o_a15
            if (p_ant3_o_a15 .ge. p1st) lptr_ant3_o_aer(15,itype,ai_phase)  = p_ant3_o_a15
            if (p_ant4_o_a15 .ge. p1st) lptr_ant4_o_aer(15,itype,ai_phase)  = p_ant4_o_a15
            if (p_biog1_o_a15 .ge. p1st) lptr_biog1_o_aer(15,itype,ai_phase)  = p_biog1_o_a15
            if (p_biog2_o_a15 .ge. p1st) lptr_biog2_o_aer(15,itype,ai_phase)  = p_biog2_o_a15
            if (p_biog3_o_a15 .ge. p1st) lptr_biog3_o_aer(15,itype,ai_phase)  = p_biog3_o_a15
            if (p_biog4_o_a15 .ge. p1st) lptr_biog4_o_aer(15,itype,ai_phase)  = p_biog4_o_a15
            
	    numptr_aer(15,itype,ai_phase)        = p_num_a15
	end if

	if (nsize_aer(itype) .ge. 16) then
	    lptr_so4_aer(16,itype,ai_phase)      = p_so4_a16
	    lptr_no3_aer(16,itype,ai_phase)      = p_no3_a16
	    lptr_cl_aer(16,itype,ai_phase)       = p_cl_a16
	    lptr_msa_aer(16,itype,ai_phase)      = p_msa_a16
	    lptr_co3_aer(16,itype,ai_phase)      = p_co3_a16
	    lptr_nh4_aer(16,itype,ai_phase)      = p_nh4_a16
	    lptr_na_aer(16,itype,ai_phase)       = p_na_a16
	    lptr_ca_aer(16,itype,ai_phase)       = p_ca_a16
	    lptr_oin_aer(16,itype,ai_phase)      = p_oin_a16
	    lptr_oc_aer(16,itype,ai_phase)       = p_oc_a16
	    lptr_bc_aer(16,itype,ai_phase)       = p_bc_a16
	    hyswptr_aer(16,itype)                = p_hysw_a16
	    waterptr_aer(16,itype)               = p_water_a16

            if (p_pcg1_b_c_a16 .ge. p1st) lptr_pcg1_b_c_aer(16,itype,ai_phase) = p_pcg1_b_c_a16
            if (p_pcg2_b_c_a16 .ge. p1st) lptr_pcg2_b_c_aer(16,itype,ai_phase)  = p_pcg2_b_c_a16
            if (p_pcg3_b_c_a16 .ge. p1st) lptr_pcg3_b_c_aer(16,itype,ai_phase)  = p_pcg3_b_c_a16
            if (p_pcg4_b_c_a16 .ge. p1st) lptr_pcg4_b_c_aer(16,itype,ai_phase)  = p_pcg4_b_c_a16
            if (p_pcg5_b_c_a16 .ge. p1st) lptr_pcg5_b_c_aer(16,itype,ai_phase)  = p_pcg5_b_c_a16
            if (p_pcg6_b_c_a16 .ge. p1st) lptr_pcg6_b_c_aer(16,itype,ai_phase)  = p_pcg6_b_c_a16
            if (p_pcg7_b_c_a16 .ge. p1st) lptr_pcg7_b_c_aer(16,itype,ai_phase)  = p_pcg7_b_c_a16
            if (p_pcg8_b_c_a16 .ge. p1st) lptr_pcg8_b_c_aer(16,itype,ai_phase)  = p_pcg8_b_c_a16
            if (p_pcg9_b_c_a16 .ge. p1st) lptr_pcg9_b_c_aer(16,itype,ai_phase)  = p_pcg9_b_c_a16
            if (p_pcg1_b_o_a16 .ge. p1st) lptr_pcg1_b_o_aer(16,itype,ai_phase) = p_pcg1_b_o_a16
            if (p_pcg2_b_o_a16 .ge. p1st) lptr_pcg2_b_o_aer(16,itype,ai_phase)  = p_pcg2_b_o_a16
            if (p_pcg3_b_o_a16 .ge. p1st) lptr_pcg3_b_o_aer(16,itype,ai_phase)  = p_pcg3_b_o_a16
            if (p_pcg4_b_o_a16 .ge. p1st) lptr_pcg4_b_o_aer(16,itype,ai_phase)  = p_pcg4_b_o_a16
            if (p_pcg5_b_o_a16 .ge. p1st) lptr_pcg5_b_o_aer(16,itype,ai_phase)  = p_pcg5_b_o_a16
            if (p_pcg6_b_o_a16 .ge. p1st) lptr_pcg6_b_o_aer(16,itype,ai_phase)  = p_pcg6_b_o_a16
            if (p_pcg7_b_o_a16 .ge. p1st) lptr_pcg7_b_o_aer(16,itype,ai_phase)  = p_pcg7_b_o_a16
            if (p_pcg8_b_o_a16 .ge. p1st) lptr_pcg8_b_o_aer(16,itype,ai_phase)  = p_pcg8_b_o_a16
            if (p_pcg9_b_o_a16 .ge. p1st) lptr_pcg9_b_o_aer(16,itype,ai_phase)  = p_pcg9_b_o_a16
            if (p_opcg1_b_c_a16 .ge. p1st) lptr_opcg1_b_c_aer(16,itype,ai_phase) = p_opcg1_b_c_a16
            if (p_opcg2_b_c_a16 .ge. p1st) lptr_opcg2_b_c_aer(16,itype,ai_phase)  = p_opcg2_b_c_a16
            if (p_opcg3_b_c_a16 .ge. p1st) lptr_opcg3_b_c_aer(16,itype,ai_phase)  = p_opcg3_b_c_a16
            if (p_opcg4_b_c_a16 .ge. p1st) lptr_opcg4_b_c_aer(16,itype,ai_phase)  = p_opcg4_b_c_a16
            if (p_opcg5_b_c_a16 .ge. p1st) lptr_opcg5_b_c_aer(16,itype,ai_phase)  = p_opcg5_b_c_a16
            if (p_opcg6_b_c_a16 .ge. p1st) lptr_opcg6_b_c_aer(16,itype,ai_phase)  = p_opcg6_b_c_a16
            if (p_opcg7_b_c_a16 .ge. p1st) lptr_opcg7_b_c_aer(16,itype,ai_phase)  = p_opcg7_b_c_a16
            if (p_opcg8_b_c_a16 .ge. p1st) lptr_opcg8_b_c_aer(16,itype,ai_phase)  = p_opcg8_b_c_a16
            if (p_opcg1_b_o_a16 .ge. p1st) lptr_opcg1_b_o_aer(16,itype,ai_phase) = p_opcg1_b_o_a16
            if (p_opcg2_b_o_a16 .ge. p1st) lptr_opcg2_b_o_aer(16,itype,ai_phase)  = p_opcg2_b_o_a16
            if (p_opcg3_b_o_a16 .ge. p1st) lptr_opcg3_b_o_aer(16,itype,ai_phase)  = p_opcg3_b_o_a16
            if (p_opcg4_b_o_a16 .ge. p1st) lptr_opcg4_b_o_aer(16,itype,ai_phase)  = p_opcg4_b_o_a16
            if (p_opcg5_b_o_a16 .ge. p1st) lptr_opcg5_b_o_aer(16,itype,ai_phase)  = p_opcg5_b_o_a16
            if (p_opcg6_b_o_a16 .ge. p1st) lptr_opcg6_b_o_aer(16,itype,ai_phase)  = p_opcg6_b_o_a16
            if (p_opcg7_b_o_a16 .ge. p1st) lptr_opcg7_b_o_aer(16,itype,ai_phase)  = p_opcg7_b_o_a16
            if (p_opcg8_b_o_a16 .ge. p1st) lptr_opcg8_b_o_aer(16,itype,ai_phase)  = p_opcg8_b_o_a16
            if (p_pcg1_f_c_a16 .ge. p1st) lptr_pcg1_f_c_aer(16,itype,ai_phase) = p_pcg1_f_c_a16
            if (p_pcg2_f_c_a16 .ge. p1st) lptr_pcg2_f_c_aer(16,itype,ai_phase)  = p_pcg2_f_c_a16
            if (p_pcg3_f_c_a16 .ge. p1st) lptr_pcg3_f_c_aer(16,itype,ai_phase)  = p_pcg3_f_c_a16
            if (p_pcg4_f_c_a16 .ge. p1st) lptr_pcg4_f_c_aer(16,itype,ai_phase)  = p_pcg4_f_c_a16
            if (p_pcg5_f_c_a16 .ge. p1st) lptr_pcg5_f_c_aer(16,itype,ai_phase)  = p_pcg5_f_c_a16
            if (p_pcg6_f_c_a16 .ge. p1st) lptr_pcg6_f_c_aer(16,itype,ai_phase)  = p_pcg6_f_c_a16
            if (p_pcg7_f_c_a16 .ge. p1st) lptr_pcg7_f_c_aer(16,itype,ai_phase)  = p_pcg7_f_c_a16
            if (p_pcg8_f_c_a16 .ge. p1st) lptr_pcg8_f_c_aer(16,itype,ai_phase)  = p_pcg8_f_c_a16
            if (p_pcg9_f_c_a16 .ge. p1st) lptr_pcg9_f_c_aer(16,itype,ai_phase)  = p_pcg9_f_c_a16
            if (p_pcg1_f_o_a16 .ge. p1st) lptr_pcg1_f_o_aer(16,itype,ai_phase) = p_pcg1_f_o_a16
            if (p_pcg2_f_o_a16 .ge. p1st) lptr_pcg2_f_o_aer(16,itype,ai_phase)  = p_pcg2_f_o_a16
            if (p_pcg3_f_o_a16 .ge. p1st) lptr_pcg3_f_o_aer(16,itype,ai_phase)  = p_pcg3_f_o_a16
            if (p_pcg4_f_o_a16 .ge. p1st) lptr_pcg4_f_o_aer(16,itype,ai_phase)  = p_pcg4_f_o_a16
            if (p_pcg5_f_o_a16 .ge. p1st) lptr_pcg5_f_o_aer(16,itype,ai_phase)  = p_pcg5_f_o_a16
            if (p_pcg6_f_o_a16 .ge. p1st) lptr_pcg6_f_o_aer(16,itype,ai_phase)  = p_pcg6_f_o_a16
            if (p_pcg7_f_o_a16 .ge. p1st) lptr_pcg7_f_o_aer(16,itype,ai_phase)  = p_pcg7_f_o_a16
            if (p_pcg8_f_o_a16 .ge. p1st) lptr_pcg8_f_o_aer(16,itype,ai_phase)  = p_pcg8_f_o_a16
            if (p_pcg9_f_o_a16 .ge. p1st) lptr_pcg9_f_o_aer(16,itype,ai_phase)  = p_pcg9_f_o_a16
            if (p_opcg1_f_c_a16 .ge. p1st) lptr_opcg1_f_c_aer(16,itype,ai_phase) = p_opcg1_f_c_a16
            if (p_opcg2_f_c_a16 .ge. p1st) lptr_opcg2_f_c_aer(16,itype,ai_phase)  = p_opcg2_f_c_a16
            if (p_opcg3_f_c_a16 .ge. p1st) lptr_opcg3_f_c_aer(16,itype,ai_phase)  = p_opcg3_f_c_a16
            if (p_opcg4_f_c_a16 .ge. p1st) lptr_opcg4_f_c_aer(16,itype,ai_phase)  = p_opcg4_f_c_a16
            if (p_opcg5_f_c_a16 .ge. p1st) lptr_opcg5_f_c_aer(16,itype,ai_phase)  = p_opcg5_f_c_a16
            if (p_opcg6_f_c_a16 .ge. p1st) lptr_opcg6_f_c_aer(16,itype,ai_phase)  = p_opcg6_f_c_a16
            if (p_opcg7_f_c_a16 .ge. p1st) lptr_opcg7_f_c_aer(16,itype,ai_phase)  = p_opcg7_f_c_a16
            if (p_opcg8_f_c_a16 .ge. p1st) lptr_opcg8_f_c_aer(16,itype,ai_phase)  = p_opcg8_f_c_a16
            if (p_opcg1_f_o_a16 .ge. p1st) lptr_opcg1_f_o_aer(16,itype,ai_phase) = p_opcg1_f_o_a16
            if (p_opcg2_f_o_a16 .ge. p1st) lptr_opcg2_f_o_aer(16,itype,ai_phase)  = p_opcg2_f_o_a16
            if (p_opcg3_f_o_a16 .ge. p1st) lptr_opcg3_f_o_aer(16,itype,ai_phase)  = p_opcg3_f_o_a16
            if (p_opcg4_f_o_a16 .ge. p1st) lptr_opcg4_f_o_aer(16,itype,ai_phase)  = p_opcg4_f_o_a16
            if (p_opcg5_f_o_a16 .ge. p1st) lptr_opcg5_f_o_aer(16,itype,ai_phase)  = p_opcg5_f_o_a16
            if (p_opcg6_f_o_a16 .ge. p1st) lptr_opcg6_f_o_aer(16,itype,ai_phase)  = p_opcg6_f_o_a16
            if (p_opcg7_f_o_a16 .ge. p1st) lptr_opcg7_f_o_aer(16,itype,ai_phase)  = p_opcg7_f_o_a16
            if (p_opcg8_f_o_a16 .ge. p1st) lptr_opcg8_f_o_aer(16,itype,ai_phase)  = p_opcg8_f_o_a16
            if (p_ant1_c_a16 .ge. p1st) lptr_ant1_c_aer(16,itype,ai_phase)  = p_ant1_c_a16
            if (p_ant2_c_a16 .ge. p1st) lptr_ant2_c_aer(16,itype,ai_phase)  = p_ant2_c_a16
            if (p_ant3_c_a16 .ge. p1st) lptr_ant3_c_aer(16,itype,ai_phase)  = p_ant3_c_a16
            if (p_ant4_c_a16 .ge. p1st) lptr_ant4_c_aer(16,itype,ai_phase)  = p_ant4_c_a16
            if (p_biog1_c_a16 .ge. p1st) lptr_biog1_c_aer(16,itype,ai_phase)  = p_biog1_c_a16
            if (p_biog2_c_a16 .ge. p1st) lptr_biog2_c_aer(16,itype,ai_phase)  = p_biog2_c_a16
            if (p_biog3_c_a16 .ge. p1st) lptr_biog3_c_aer(16,itype,ai_phase)  = p_biog3_c_a16
            if (p_biog4_c_a16 .ge. p1st) lptr_biog4_c_aer(16,itype,ai_phase)  = p_biog4_c_a16
            if (p_ant1_o_a16 .ge. p1st) lptr_ant1_o_aer(16,itype,ai_phase)  = p_ant1_o_a16
            if (p_ant2_o_a16 .ge. p1st) lptr_ant2_o_aer(16,itype,ai_phase)  = p_ant2_o_a16
            if (p_ant3_o_a16 .ge. p1st) lptr_ant3_o_aer(16,itype,ai_phase)  = p_ant3_o_a16
            if (p_ant4_o_a16 .ge. p1st) lptr_ant4_o_aer(16,itype,ai_phase)  = p_ant4_o_a16
            if (p_biog1_o_a16 .ge. p1st) lptr_biog1_o_aer(16,itype,ai_phase)  = p_biog1_o_a16
            if (p_biog2_o_a16 .ge. p1st) lptr_biog2_o_aer(16,itype,ai_phase)  = p_biog2_o_a16
            if (p_biog3_o_a16 .ge. p1st) lptr_biog3_o_aer(16,itype,ai_phase)  = p_biog3_o_a16
            if (p_biog4_o_a16 .ge. p1st) lptr_biog4_o_aer(16,itype,ai_phase)  = p_biog4_o_a16
            
	    numptr_aer(16,itype,ai_phase)        = p_num_a16
	end if

	if (nsize_aer(itype) .ge. 17) then
	    lptr_so4_aer(17,itype,ai_phase)      = p_so4_a17
	    lptr_no3_aer(17,itype,ai_phase)      = p_no3_a17
	    lptr_cl_aer(17,itype,ai_phase)       = p_cl_a17
	    lptr_msa_aer(17,itype,ai_phase)      = p_msa_a17
	    lptr_co3_aer(17,itype,ai_phase)      = p_co3_a17
	    lptr_nh4_aer(17,itype,ai_phase)      = p_nh4_a17
	    lptr_na_aer(17,itype,ai_phase)       = p_na_a17
	    lptr_ca_aer(17,itype,ai_phase)       = p_ca_a17
	    lptr_oin_aer(17,itype,ai_phase)      = p_oin_a17
	    lptr_oc_aer(17,itype,ai_phase)       = p_oc_a17
	    lptr_bc_aer(17,itype,ai_phase)       = p_bc_a17
	    hyswptr_aer(17,itype)                = p_hysw_a17
	    waterptr_aer(17,itype)               = p_water_a17

            if (p_pcg1_b_c_a17 .ge. p1st) lptr_pcg1_b_c_aer(17,itype,ai_phase) = p_pcg1_b_c_a17
            if (p_pcg2_b_c_a17 .ge. p1st) lptr_pcg2_b_c_aer(17,itype,ai_phase)  = p_pcg2_b_c_a17
            if (p_pcg3_b_c_a17 .ge. p1st) lptr_pcg3_b_c_aer(17,itype,ai_phase)  = p_pcg3_b_c_a17
            if (p_pcg4_b_c_a17 .ge. p1st) lptr_pcg4_b_c_aer(17,itype,ai_phase)  = p_pcg4_b_c_a17
            if (p_pcg5_b_c_a17 .ge. p1st) lptr_pcg5_b_c_aer(17,itype,ai_phase)  = p_pcg5_b_c_a17
            if (p_pcg6_b_c_a17 .ge. p1st) lptr_pcg6_b_c_aer(17,itype,ai_phase)  = p_pcg6_b_c_a17
            if (p_pcg7_b_c_a17 .ge. p1st) lptr_pcg7_b_c_aer(17,itype,ai_phase)  = p_pcg7_b_c_a17
            if (p_pcg8_b_c_a17 .ge. p1st) lptr_pcg8_b_c_aer(17,itype,ai_phase)  = p_pcg8_b_c_a17
            if (p_pcg9_b_c_a17 .ge. p1st) lptr_pcg9_b_c_aer(17,itype,ai_phase)  = p_pcg9_b_c_a17
            if (p_pcg1_b_o_a17 .ge. p1st) lptr_pcg1_b_o_aer(17,itype,ai_phase) = p_pcg1_b_o_a17
            if (p_pcg2_b_o_a17 .ge. p1st) lptr_pcg2_b_o_aer(17,itype,ai_phase)  = p_pcg2_b_o_a17
            if (p_pcg3_b_o_a17 .ge. p1st) lptr_pcg3_b_o_aer(17,itype,ai_phase)  = p_pcg3_b_o_a17
            if (p_pcg4_b_o_a17 .ge. p1st) lptr_pcg4_b_o_aer(17,itype,ai_phase)  = p_pcg4_b_o_a17
            if (p_pcg5_b_o_a17 .ge. p1st) lptr_pcg5_b_o_aer(17,itype,ai_phase)  = p_pcg5_b_o_a17
            if (p_pcg6_b_o_a17 .ge. p1st) lptr_pcg6_b_o_aer(17,itype,ai_phase)  = p_pcg6_b_o_a17
            if (p_pcg7_b_o_a17 .ge. p1st) lptr_pcg7_b_o_aer(17,itype,ai_phase)  = p_pcg7_b_o_a17
            if (p_pcg8_b_o_a17 .ge. p1st) lptr_pcg8_b_o_aer(17,itype,ai_phase)  = p_pcg8_b_o_a17
            if (p_pcg9_b_o_a17 .ge. p1st) lptr_pcg9_b_o_aer(17,itype,ai_phase)  = p_pcg9_b_o_a17
            if (p_opcg1_b_c_a17 .ge. p1st) lptr_opcg1_b_c_aer(17,itype,ai_phase) = p_opcg1_b_c_a17
            if (p_opcg2_b_c_a17 .ge. p1st) lptr_opcg2_b_c_aer(17,itype,ai_phase)  = p_opcg2_b_c_a17
            if (p_opcg3_b_c_a17 .ge. p1st) lptr_opcg3_b_c_aer(17,itype,ai_phase)  = p_opcg3_b_c_a17
            if (p_opcg4_b_c_a17 .ge. p1st) lptr_opcg4_b_c_aer(17,itype,ai_phase)  = p_opcg4_b_c_a17
            if (p_opcg5_b_c_a17 .ge. p1st) lptr_opcg5_b_c_aer(17,itype,ai_phase)  = p_opcg5_b_c_a17
            if (p_opcg6_b_c_a17 .ge. p1st) lptr_opcg6_b_c_aer(17,itype,ai_phase)  = p_opcg6_b_c_a17
            if (p_opcg7_b_c_a17 .ge. p1st) lptr_opcg7_b_c_aer(17,itype,ai_phase)  = p_opcg7_b_c_a17
            if (p_opcg8_b_c_a17 .ge. p1st) lptr_opcg8_b_c_aer(17,itype,ai_phase)  = p_opcg8_b_c_a17
            if (p_opcg1_b_o_a17 .ge. p1st) lptr_opcg1_b_o_aer(17,itype,ai_phase) = p_opcg1_b_o_a17
            if (p_opcg2_b_o_a17 .ge. p1st) lptr_opcg2_b_o_aer(17,itype,ai_phase)  = p_opcg2_b_o_a17
            if (p_opcg3_b_o_a17 .ge. p1st) lptr_opcg3_b_o_aer(17,itype,ai_phase)  = p_opcg3_b_o_a17
            if (p_opcg4_b_o_a17 .ge. p1st) lptr_opcg4_b_o_aer(17,itype,ai_phase)  = p_opcg4_b_o_a17
            if (p_opcg5_b_o_a17 .ge. p1st) lptr_opcg5_b_o_aer(17,itype,ai_phase)  = p_opcg5_b_o_a17
            if (p_opcg6_b_o_a17 .ge. p1st) lptr_opcg6_b_o_aer(17,itype,ai_phase)  = p_opcg6_b_o_a17
            if (p_opcg7_b_o_a17 .ge. p1st) lptr_opcg7_b_o_aer(17,itype,ai_phase)  = p_opcg7_b_o_a17
            if (p_opcg8_b_o_a17 .ge. p1st) lptr_opcg8_b_o_aer(17,itype,ai_phase)  = p_opcg8_b_o_a17
            if (p_pcg1_f_c_a17 .ge. p1st) lptr_pcg1_f_c_aer(17,itype,ai_phase) = p_pcg1_f_c_a17
            if (p_pcg2_f_c_a17 .ge. p1st) lptr_pcg2_f_c_aer(17,itype,ai_phase)  = p_pcg2_f_c_a17
            if (p_pcg3_f_c_a17 .ge. p1st) lptr_pcg3_f_c_aer(17,itype,ai_phase)  = p_pcg3_f_c_a17
            if (p_pcg4_f_c_a17 .ge. p1st) lptr_pcg4_f_c_aer(17,itype,ai_phase)  = p_pcg4_f_c_a17
            if (p_pcg5_f_c_a17 .ge. p1st) lptr_pcg5_f_c_aer(17,itype,ai_phase)  = p_pcg5_f_c_a17
            if (p_pcg6_f_c_a17 .ge. p1st) lptr_pcg6_f_c_aer(17,itype,ai_phase)  = p_pcg6_f_c_a17
            if (p_pcg7_f_c_a17 .ge. p1st) lptr_pcg7_f_c_aer(17,itype,ai_phase)  = p_pcg7_f_c_a17
            if (p_pcg8_f_c_a17 .ge. p1st) lptr_pcg8_f_c_aer(17,itype,ai_phase)  = p_pcg8_f_c_a17
            if (p_pcg9_f_c_a17 .ge. p1st) lptr_pcg9_f_c_aer(17,itype,ai_phase)  = p_pcg9_f_c_a17
            if (p_pcg1_f_o_a17 .ge. p1st) lptr_pcg1_f_o_aer(17,itype,ai_phase) = p_pcg1_f_o_a17
            if (p_pcg2_f_o_a17 .ge. p1st) lptr_pcg2_f_o_aer(17,itype,ai_phase)  = p_pcg2_f_o_a17
            if (p_pcg3_f_o_a17 .ge. p1st) lptr_pcg3_f_o_aer(17,itype,ai_phase)  = p_pcg3_f_o_a17
            if (p_pcg4_f_o_a17 .ge. p1st) lptr_pcg4_f_o_aer(17,itype,ai_phase)  = p_pcg4_f_o_a17
            if (p_pcg5_f_o_a17 .ge. p1st) lptr_pcg5_f_o_aer(17,itype,ai_phase)  = p_pcg5_f_o_a17
            if (p_pcg6_f_o_a17 .ge. p1st) lptr_pcg6_f_o_aer(17,itype,ai_phase)  = p_pcg6_f_o_a17
            if (p_pcg7_f_o_a17 .ge. p1st) lptr_pcg7_f_o_aer(17,itype,ai_phase)  = p_pcg7_f_o_a17
            if (p_pcg8_f_o_a17 .ge. p1st) lptr_pcg8_f_o_aer(17,itype,ai_phase)  = p_pcg8_f_o_a17
            if (p_pcg9_f_o_a17 .ge. p1st) lptr_pcg9_f_o_aer(17,itype,ai_phase)  = p_pcg9_f_o_a17
            if (p_opcg1_f_c_a17 .ge. p1st) lptr_opcg1_f_c_aer(17,itype,ai_phase) = p_opcg1_f_c_a17
            if (p_opcg2_f_c_a17 .ge. p1st) lptr_opcg2_f_c_aer(17,itype,ai_phase)  = p_opcg2_f_c_a17
            if (p_opcg3_f_c_a17 .ge. p1st) lptr_opcg3_f_c_aer(17,itype,ai_phase)  = p_opcg3_f_c_a17
            if (p_opcg4_f_c_a17 .ge. p1st) lptr_opcg4_f_c_aer(17,itype,ai_phase)  = p_opcg4_f_c_a17
            if (p_opcg5_f_c_a17 .ge. p1st) lptr_opcg5_f_c_aer(17,itype,ai_phase)  = p_opcg5_f_c_a17
            if (p_opcg6_f_c_a17 .ge. p1st) lptr_opcg6_f_c_aer(17,itype,ai_phase)  = p_opcg6_f_c_a17
            if (p_opcg7_f_c_a17 .ge. p1st) lptr_opcg7_f_c_aer(17,itype,ai_phase)  = p_opcg7_f_c_a17
            if (p_opcg8_f_c_a17 .ge. p1st) lptr_opcg8_f_c_aer(17,itype,ai_phase)  = p_opcg8_f_c_a17
            if (p_opcg1_f_o_a17 .ge. p1st) lptr_opcg1_f_o_aer(17,itype,ai_phase) = p_opcg1_f_o_a17
            if (p_opcg2_f_o_a17 .ge. p1st) lptr_opcg2_f_o_aer(17,itype,ai_phase)  = p_opcg2_f_o_a17
            if (p_opcg3_f_o_a17 .ge. p1st) lptr_opcg3_f_o_aer(17,itype,ai_phase)  = p_opcg3_f_o_a17
            if (p_opcg4_f_o_a17 .ge. p1st) lptr_opcg4_f_o_aer(17,itype,ai_phase)  = p_opcg4_f_o_a17
            if (p_opcg5_f_o_a17 .ge. p1st) lptr_opcg5_f_o_aer(17,itype,ai_phase)  = p_opcg5_f_o_a17
            if (p_opcg6_f_o_a17 .ge. p1st) lptr_opcg6_f_o_aer(17,itype,ai_phase)  = p_opcg6_f_o_a17
            if (p_opcg7_f_o_a17 .ge. p1st) lptr_opcg7_f_o_aer(17,itype,ai_phase)  = p_opcg7_f_o_a17
            if (p_opcg8_f_o_a17 .ge. p1st) lptr_opcg8_f_o_aer(17,itype,ai_phase)  = p_opcg8_f_o_a17
            if (p_ant1_c_a17 .ge. p1st) lptr_ant1_c_aer(17,itype,ai_phase)  = p_ant1_c_a17
            if (p_ant2_c_a17 .ge. p1st) lptr_ant2_c_aer(17,itype,ai_phase)  = p_ant2_c_a17
            if (p_ant3_c_a17 .ge. p1st) lptr_ant3_c_aer(17,itype,ai_phase)  = p_ant3_c_a17
            if (p_ant4_c_a17 .ge. p1st) lptr_ant4_c_aer(17,itype,ai_phase)  = p_ant4_c_a17
            if (p_biog1_c_a17 .ge. p1st) lptr_biog1_c_aer(17,itype,ai_phase)  = p_biog1_c_a17
            if (p_biog2_c_a17 .ge. p1st) lptr_biog2_c_aer(17,itype,ai_phase)  = p_biog2_c_a17
            if (p_biog3_c_a17 .ge. p1st) lptr_biog3_c_aer(17,itype,ai_phase)  = p_biog3_c_a17
            if (p_biog4_c_a17 .ge. p1st) lptr_biog4_c_aer(17,itype,ai_phase)  = p_biog4_c_a17
            if (p_ant1_o_a17 .ge. p1st) lptr_ant1_o_aer(17,itype,ai_phase)  = p_ant1_o_a17
            if (p_ant2_o_a17 .ge. p1st) lptr_ant2_o_aer(17,itype,ai_phase)  = p_ant2_o_a17
            if (p_ant3_o_a17 .ge. p1st) lptr_ant3_o_aer(17,itype,ai_phase)  = p_ant3_o_a17
            if (p_ant4_o_a17 .ge. p1st) lptr_ant4_o_aer(17,itype,ai_phase)  = p_ant4_o_a17
            if (p_biog1_o_a17 .ge. p1st) lptr_biog1_o_aer(17,itype,ai_phase)  = p_biog1_o_a17
            if (p_biog2_o_a17 .ge. p1st) lptr_biog2_o_aer(17,itype,ai_phase)  = p_biog2_o_a17
            if (p_biog3_o_a17 .ge. p1st) lptr_biog3_o_aer(17,itype,ai_phase)  = p_biog3_o_a17
            if (p_biog4_o_a17 .ge. p1st) lptr_biog4_o_aer(17,itype,ai_phase)  = p_biog4_o_a17
            
	    numptr_aer(17,itype,ai_phase)        = p_num_a17
	end if

	if (nsize_aer(itype) .ge. 18) then
	    lptr_so4_aer(18,itype,ai_phase)      = p_so4_a18
	    lptr_no3_aer(18,itype,ai_phase)      = p_no3_a18
	    lptr_cl_aer(18,itype,ai_phase)       = p_cl_a18
	    lptr_msa_aer(18,itype,ai_phase)      = p_msa_a18
	    lptr_co3_aer(18,itype,ai_phase)      = p_co3_a18
	    lptr_nh4_aer(18,itype,ai_phase)      = p_nh4_a18
	    lptr_na_aer(18,itype,ai_phase)       = p_na_a18
	    lptr_ca_aer(18,itype,ai_phase)       = p_ca_a18
	    lptr_oin_aer(18,itype,ai_phase)      = p_oin_a18
	    lptr_oc_aer(18,itype,ai_phase)       = p_oc_a18
	    lptr_bc_aer(18,itype,ai_phase)       = p_bc_a18
	    hyswptr_aer(18,itype)                = p_hysw_a18
	    waterptr_aer(18,itype)               = p_water_a18

            if (p_pcg1_b_c_a18 .ge. p1st) lptr_pcg1_b_c_aer(18,itype,ai_phase) = p_pcg1_b_c_a18
            if (p_pcg2_b_c_a18 .ge. p1st) lptr_pcg2_b_c_aer(18,itype,ai_phase)  = p_pcg2_b_c_a18
            if (p_pcg3_b_c_a18 .ge. p1st) lptr_pcg3_b_c_aer(18,itype,ai_phase)  = p_pcg3_b_c_a18
            if (p_pcg4_b_c_a18 .ge. p1st) lptr_pcg4_b_c_aer(18,itype,ai_phase)  = p_pcg4_b_c_a18
            if (p_pcg5_b_c_a18 .ge. p1st) lptr_pcg5_b_c_aer(18,itype,ai_phase)  = p_pcg5_b_c_a18
            if (p_pcg6_b_c_a18 .ge. p1st) lptr_pcg6_b_c_aer(18,itype,ai_phase)  = p_pcg6_b_c_a18
            if (p_pcg7_b_c_a18 .ge. p1st) lptr_pcg7_b_c_aer(18,itype,ai_phase)  = p_pcg7_b_c_a18
            if (p_pcg8_b_c_a18 .ge. p1st) lptr_pcg8_b_c_aer(18,itype,ai_phase)  = p_pcg8_b_c_a18
            if (p_pcg9_b_c_a18 .ge. p1st) lptr_pcg9_b_c_aer(18,itype,ai_phase)  = p_pcg9_b_c_a18
            if (p_pcg1_b_o_a18 .ge. p1st) lptr_pcg1_b_o_aer(18,itype,ai_phase) = p_pcg1_b_o_a18
            if (p_pcg2_b_o_a18 .ge. p1st) lptr_pcg2_b_o_aer(18,itype,ai_phase)  = p_pcg2_b_o_a18
            if (p_pcg3_b_o_a18 .ge. p1st) lptr_pcg3_b_o_aer(18,itype,ai_phase)  = p_pcg3_b_o_a18
            if (p_pcg4_b_o_a18 .ge. p1st) lptr_pcg4_b_o_aer(18,itype,ai_phase)  = p_pcg4_b_o_a18
            if (p_pcg5_b_o_a18 .ge. p1st) lptr_pcg5_b_o_aer(18,itype,ai_phase)  = p_pcg5_b_o_a18
            if (p_pcg6_b_o_a18 .ge. p1st) lptr_pcg6_b_o_aer(18,itype,ai_phase)  = p_pcg6_b_o_a18
            if (p_pcg7_b_o_a18 .ge. p1st) lptr_pcg7_b_o_aer(18,itype,ai_phase)  = p_pcg7_b_o_a18
            if (p_pcg8_b_o_a18 .ge. p1st) lptr_pcg8_b_o_aer(18,itype,ai_phase)  = p_pcg8_b_o_a18
            if (p_pcg9_b_o_a18 .ge. p1st) lptr_pcg9_b_o_aer(18,itype,ai_phase)  = p_pcg9_b_o_a18
            if (p_opcg1_b_c_a18 .ge. p1st) lptr_opcg1_b_c_aer(18,itype,ai_phase) = p_opcg1_b_c_a18
            if (p_opcg2_b_c_a18 .ge. p1st) lptr_opcg2_b_c_aer(18,itype,ai_phase)  = p_opcg2_b_c_a18
            if (p_opcg3_b_c_a18 .ge. p1st) lptr_opcg3_b_c_aer(18,itype,ai_phase)  = p_opcg3_b_c_a18
            if (p_opcg4_b_c_a18 .ge. p1st) lptr_opcg4_b_c_aer(18,itype,ai_phase)  = p_opcg4_b_c_a18
            if (p_opcg5_b_c_a18 .ge. p1st) lptr_opcg5_b_c_aer(18,itype,ai_phase)  = p_opcg5_b_c_a18
            if (p_opcg6_b_c_a18 .ge. p1st) lptr_opcg6_b_c_aer(18,itype,ai_phase)  = p_opcg6_b_c_a18
            if (p_opcg7_b_c_a18 .ge. p1st) lptr_opcg7_b_c_aer(18,itype,ai_phase)  = p_opcg7_b_c_a18
            if (p_opcg8_b_c_a18 .ge. p1st) lptr_opcg8_b_c_aer(18,itype,ai_phase)  = p_opcg8_b_c_a18
            if (p_opcg1_b_o_a18 .ge. p1st) lptr_opcg1_b_o_aer(18,itype,ai_phase) = p_opcg1_b_o_a18
            if (p_opcg2_b_o_a18 .ge. p1st) lptr_opcg2_b_o_aer(18,itype,ai_phase)  = p_opcg2_b_o_a18
            if (p_opcg3_b_o_a18 .ge. p1st) lptr_opcg3_b_o_aer(18,itype,ai_phase)  = p_opcg3_b_o_a18
            if (p_opcg4_b_o_a18 .ge. p1st) lptr_opcg4_b_o_aer(18,itype,ai_phase)  = p_opcg4_b_o_a18
            if (p_opcg5_b_o_a18 .ge. p1st) lptr_opcg5_b_o_aer(18,itype,ai_phase)  = p_opcg5_b_o_a18
            if (p_opcg6_b_o_a18 .ge. p1st) lptr_opcg6_b_o_aer(18,itype,ai_phase)  = p_opcg6_b_o_a18
            if (p_opcg7_b_o_a18 .ge. p1st) lptr_opcg7_b_o_aer(18,itype,ai_phase)  = p_opcg7_b_o_a18
            if (p_opcg8_b_o_a18 .ge. p1st) lptr_opcg8_b_o_aer(18,itype,ai_phase)  = p_opcg8_b_o_a18
            if (p_pcg1_f_c_a18 .ge. p1st) lptr_pcg1_f_c_aer(18,itype,ai_phase) = p_pcg1_f_c_a18
            if (p_pcg2_f_c_a18 .ge. p1st) lptr_pcg2_f_c_aer(18,itype,ai_phase)  = p_pcg2_f_c_a18
            if (p_pcg3_f_c_a18 .ge. p1st) lptr_pcg3_f_c_aer(18,itype,ai_phase)  = p_pcg3_f_c_a18
            if (p_pcg4_f_c_a18 .ge. p1st) lptr_pcg4_f_c_aer(18,itype,ai_phase)  = p_pcg4_f_c_a18
            if (p_pcg5_f_c_a18 .ge. p1st) lptr_pcg5_f_c_aer(18,itype,ai_phase)  = p_pcg5_f_c_a18
            if (p_pcg6_f_c_a18 .ge. p1st) lptr_pcg6_f_c_aer(18,itype,ai_phase)  = p_pcg6_f_c_a18
            if (p_pcg7_f_c_a18 .ge. p1st) lptr_pcg7_f_c_aer(18,itype,ai_phase)  = p_pcg7_f_c_a18
            if (p_pcg8_f_c_a18 .ge. p1st) lptr_pcg8_f_c_aer(18,itype,ai_phase)  = p_pcg8_f_c_a18
            if (p_pcg9_f_c_a18 .ge. p1st) lptr_pcg9_f_c_aer(18,itype,ai_phase)  = p_pcg9_f_c_a18
            if (p_pcg1_f_o_a18 .ge. p1st) lptr_pcg1_f_o_aer(18,itype,ai_phase) = p_pcg1_f_o_a18
            if (p_pcg2_f_o_a18 .ge. p1st) lptr_pcg2_f_o_aer(18,itype,ai_phase)  = p_pcg2_f_o_a18
            if (p_pcg3_f_o_a18 .ge. p1st) lptr_pcg3_f_o_aer(18,itype,ai_phase)  = p_pcg3_f_o_a18
            if (p_pcg4_f_o_a18 .ge. p1st) lptr_pcg4_f_o_aer(18,itype,ai_phase)  = p_pcg4_f_o_a18
            if (p_pcg5_f_o_a18 .ge. p1st) lptr_pcg5_f_o_aer(18,itype,ai_phase)  = p_pcg5_f_o_a18
            if (p_pcg6_f_o_a18 .ge. p1st) lptr_pcg6_f_o_aer(18,itype,ai_phase)  = p_pcg6_f_o_a18
            if (p_pcg7_f_o_a18 .ge. p1st) lptr_pcg7_f_o_aer(18,itype,ai_phase)  = p_pcg7_f_o_a18
            if (p_pcg8_f_o_a18 .ge. p1st) lptr_pcg8_f_o_aer(18,itype,ai_phase)  = p_pcg8_f_o_a18
            if (p_pcg9_f_o_a18 .ge. p1st) lptr_pcg9_f_o_aer(18,itype,ai_phase)  = p_pcg9_f_o_a18
            if (p_opcg1_f_c_a18 .ge. p1st) lptr_opcg1_f_c_aer(18,itype,ai_phase) = p_opcg1_f_c_a18
            if (p_opcg2_f_c_a18 .ge. p1st) lptr_opcg2_f_c_aer(18,itype,ai_phase)  = p_opcg2_f_c_a18
            if (p_opcg3_f_c_a18 .ge. p1st) lptr_opcg3_f_c_aer(18,itype,ai_phase)  = p_opcg3_f_c_a18
            if (p_opcg4_f_c_a18 .ge. p1st) lptr_opcg4_f_c_aer(18,itype,ai_phase)  = p_opcg4_f_c_a18
            if (p_opcg5_f_c_a18 .ge. p1st) lptr_opcg5_f_c_aer(18,itype,ai_phase)  = p_opcg5_f_c_a18
            if (p_opcg6_f_c_a18 .ge. p1st) lptr_opcg6_f_c_aer(18,itype,ai_phase)  = p_opcg6_f_c_a18
            if (p_opcg7_f_c_a18 .ge. p1st) lptr_opcg7_f_c_aer(18,itype,ai_phase)  = p_opcg7_f_c_a18
            if (p_opcg8_f_c_a18 .ge. p1st) lptr_opcg8_f_c_aer(18,itype,ai_phase)  = p_opcg8_f_c_a18
            if (p_opcg1_f_o_a18 .ge. p1st) lptr_opcg1_f_o_aer(18,itype,ai_phase) = p_opcg1_f_o_a18
            if (p_opcg2_f_o_a18 .ge. p1st) lptr_opcg2_f_o_aer(18,itype,ai_phase)  = p_opcg2_f_o_a18
            if (p_opcg3_f_o_a18 .ge. p1st) lptr_opcg3_f_o_aer(18,itype,ai_phase)  = p_opcg3_f_o_a18
            if (p_opcg4_f_o_a18 .ge. p1st) lptr_opcg4_f_o_aer(18,itype,ai_phase)  = p_opcg4_f_o_a18
            if (p_opcg5_f_o_a18 .ge. p1st) lptr_opcg5_f_o_aer(18,itype,ai_phase)  = p_opcg5_f_o_a18
            if (p_opcg6_f_o_a18 .ge. p1st) lptr_opcg6_f_o_aer(18,itype,ai_phase)  = p_opcg6_f_o_a18
            if (p_opcg7_f_o_a18 .ge. p1st) lptr_opcg7_f_o_aer(18,itype,ai_phase)  = p_opcg7_f_o_a18
            if (p_opcg8_f_o_a18 .ge. p1st) lptr_opcg8_f_o_aer(18,itype,ai_phase)  = p_opcg8_f_o_a18
            if (p_ant1_c_a18 .ge. p1st) lptr_ant1_c_aer(18,itype,ai_phase)  = p_ant1_c_a18
            if (p_ant2_c_a18 .ge. p1st) lptr_ant2_c_aer(18,itype,ai_phase)  = p_ant2_c_a18
            if (p_ant3_c_a18 .ge. p1st) lptr_ant3_c_aer(18,itype,ai_phase)  = p_ant3_c_a18
            if (p_ant4_c_a18 .ge. p1st) lptr_ant4_c_aer(18,itype,ai_phase)  = p_ant4_c_a18
            if (p_biog1_c_a18 .ge. p1st) lptr_biog1_c_aer(18,itype,ai_phase)  = p_biog1_c_a18
            if (p_biog2_c_a18 .ge. p1st) lptr_biog2_c_aer(18,itype,ai_phase)  = p_biog2_c_a18
            if (p_biog3_c_a18 .ge. p1st) lptr_biog3_c_aer(18,itype,ai_phase)  = p_biog3_c_a18
            if (p_biog4_c_a18 .ge. p1st) lptr_biog4_c_aer(18,itype,ai_phase)  = p_biog4_c_a18
            if (p_ant1_o_a18 .ge. p1st) lptr_ant1_o_aer(18,itype,ai_phase)  = p_ant1_o_a18
            if (p_ant2_o_a18 .ge. p1st) lptr_ant2_o_aer(18,itype,ai_phase)  = p_ant2_o_a18
            if (p_ant3_o_a18 .ge. p1st) lptr_ant3_o_aer(18,itype,ai_phase)  = p_ant3_o_a18
            if (p_ant4_o_a18 .ge. p1st) lptr_ant4_o_aer(18,itype,ai_phase)  = p_ant4_o_a18
            if (p_biog1_o_a18 .ge. p1st) lptr_biog1_o_aer(18,itype,ai_phase)  = p_biog1_o_a18
            if (p_biog2_o_a18 .ge. p1st) lptr_biog2_o_aer(18,itype,ai_phase)  = p_biog2_o_a18
            if (p_biog3_o_a18 .ge. p1st) lptr_biog3_o_aer(18,itype,ai_phase)  = p_biog3_o_a18
            if (p_biog4_o_a18 .ge. p1st) lptr_biog4_o_aer(18,itype,ai_phase)  = p_biog4_o_a18
            
	    numptr_aer(18,itype,ai_phase)        = p_num_a18
	end if

	if (nsize_aer(itype) .ge. 19) then
	    lptr_so4_aer(19,itype,ai_phase)      = p_so4_a19
	    lptr_no3_aer(19,itype,ai_phase)      = p_no3_a19
	    lptr_cl_aer(19,itype,ai_phase)       = p_cl_a19
	    lptr_msa_aer(19,itype,ai_phase)      = p_msa_a19
	    lptr_co3_aer(19,itype,ai_phase)      = p_co3_a19
	    lptr_nh4_aer(19,itype,ai_phase)      = p_nh4_a19
	    lptr_na_aer(19,itype,ai_phase)       = p_na_a19
	    lptr_ca_aer(19,itype,ai_phase)       = p_ca_a19
	    lptr_oin_aer(19,itype,ai_phase)      = p_oin_a19
	    lptr_oc_aer(19,itype,ai_phase)       = p_oc_a19
	    lptr_bc_aer(19,itype,ai_phase)       = p_bc_a19
	    hyswptr_aer(19,itype)                = p_hysw_a19
	    waterptr_aer(19,itype)               = p_water_a19

            if (p_pcg1_b_c_a19 .ge. p1st) lptr_pcg1_b_c_aer(19,itype,ai_phase) = p_pcg1_b_c_a19
            if (p_pcg2_b_c_a19 .ge. p1st) lptr_pcg2_b_c_aer(19,itype,ai_phase)  = p_pcg2_b_c_a19
            if (p_pcg3_b_c_a19 .ge. p1st) lptr_pcg3_b_c_aer(19,itype,ai_phase)  = p_pcg3_b_c_a19
            if (p_pcg4_b_c_a19 .ge. p1st) lptr_pcg4_b_c_aer(19,itype,ai_phase)  = p_pcg4_b_c_a19
            if (p_pcg5_b_c_a19 .ge. p1st) lptr_pcg5_b_c_aer(19,itype,ai_phase)  = p_pcg5_b_c_a19
            if (p_pcg6_b_c_a19 .ge. p1st) lptr_pcg6_b_c_aer(19,itype,ai_phase)  = p_pcg6_b_c_a19
            if (p_pcg7_b_c_a19 .ge. p1st) lptr_pcg7_b_c_aer(19,itype,ai_phase)  = p_pcg7_b_c_a19
            if (p_pcg8_b_c_a19 .ge. p1st) lptr_pcg8_b_c_aer(19,itype,ai_phase)  = p_pcg8_b_c_a19
            if (p_pcg9_b_c_a19 .ge. p1st) lptr_pcg9_b_c_aer(19,itype,ai_phase)  = p_pcg9_b_c_a19
            if (p_pcg1_b_o_a19 .ge. p1st) lptr_pcg1_b_o_aer(19,itype,ai_phase) = p_pcg1_b_o_a19
            if (p_pcg2_b_o_a19 .ge. p1st) lptr_pcg2_b_o_aer(19,itype,ai_phase)  = p_pcg2_b_o_a19
            if (p_pcg3_b_o_a19 .ge. p1st) lptr_pcg3_b_o_aer(19,itype,ai_phase)  = p_pcg3_b_o_a19
            if (p_pcg4_b_o_a19 .ge. p1st) lptr_pcg4_b_o_aer(19,itype,ai_phase)  = p_pcg4_b_o_a19
            if (p_pcg5_b_o_a19 .ge. p1st) lptr_pcg5_b_o_aer(19,itype,ai_phase)  = p_pcg5_b_o_a19
            if (p_pcg6_b_o_a19 .ge. p1st) lptr_pcg6_b_o_aer(19,itype,ai_phase)  = p_pcg6_b_o_a19
            if (p_pcg7_b_o_a19 .ge. p1st) lptr_pcg7_b_o_aer(19,itype,ai_phase)  = p_pcg7_b_o_a19
            if (p_pcg8_b_o_a19 .ge. p1st) lptr_pcg8_b_o_aer(19,itype,ai_phase)  = p_pcg8_b_o_a19
            if (p_pcg9_b_o_a19 .ge. p1st) lptr_pcg9_b_o_aer(19,itype,ai_phase)  = p_pcg9_b_o_a19
            if (p_opcg1_b_c_a19 .ge. p1st) lptr_opcg1_b_c_aer(19,itype,ai_phase) = p_opcg1_b_c_a19
            if (p_opcg2_b_c_a19 .ge. p1st) lptr_opcg2_b_c_aer(19,itype,ai_phase)  = p_opcg2_b_c_a19
            if (p_opcg3_b_c_a19 .ge. p1st) lptr_opcg3_b_c_aer(19,itype,ai_phase)  = p_opcg3_b_c_a19
            if (p_opcg4_b_c_a19 .ge. p1st) lptr_opcg4_b_c_aer(19,itype,ai_phase)  = p_opcg4_b_c_a19
            if (p_opcg5_b_c_a19 .ge. p1st) lptr_opcg5_b_c_aer(19,itype,ai_phase)  = p_opcg5_b_c_a19
            if (p_opcg6_b_c_a19 .ge. p1st) lptr_opcg6_b_c_aer(19,itype,ai_phase)  = p_opcg6_b_c_a19
            if (p_opcg7_b_c_a19 .ge. p1st) lptr_opcg7_b_c_aer(19,itype,ai_phase)  = p_opcg7_b_c_a19
            if (p_opcg8_b_c_a19 .ge. p1st) lptr_opcg8_b_c_aer(19,itype,ai_phase)  = p_opcg8_b_c_a19
            if (p_opcg1_b_o_a19 .ge. p1st) lptr_opcg1_b_o_aer(19,itype,ai_phase) = p_opcg1_b_o_a19
            if (p_opcg2_b_o_a19 .ge. p1st) lptr_opcg2_b_o_aer(19,itype,ai_phase)  = p_opcg2_b_o_a19
            if (p_opcg3_b_o_a19 .ge. p1st) lptr_opcg3_b_o_aer(19,itype,ai_phase)  = p_opcg3_b_o_a19
            if (p_opcg4_b_o_a19 .ge. p1st) lptr_opcg4_b_o_aer(19,itype,ai_phase)  = p_opcg4_b_o_a19
            if (p_opcg5_b_o_a19 .ge. p1st) lptr_opcg5_b_o_aer(19,itype,ai_phase)  = p_opcg5_b_o_a19
            if (p_opcg6_b_o_a19 .ge. p1st) lptr_opcg6_b_o_aer(19,itype,ai_phase)  = p_opcg6_b_o_a19
            if (p_opcg7_b_o_a19 .ge. p1st) lptr_opcg7_b_o_aer(19,itype,ai_phase)  = p_opcg7_b_o_a19
            if (p_opcg8_b_o_a19 .ge. p1st) lptr_opcg8_b_o_aer(19,itype,ai_phase)  = p_opcg8_b_o_a19
            if (p_pcg1_f_c_a19 .ge. p1st) lptr_pcg1_f_c_aer(19,itype,ai_phase) = p_pcg1_f_c_a19
            if (p_pcg2_f_c_a19 .ge. p1st) lptr_pcg2_f_c_aer(19,itype,ai_phase)  = p_pcg2_f_c_a19
            if (p_pcg3_f_c_a19 .ge. p1st) lptr_pcg3_f_c_aer(19,itype,ai_phase)  = p_pcg3_f_c_a19
            if (p_pcg4_f_c_a19 .ge. p1st) lptr_pcg4_f_c_aer(19,itype,ai_phase)  = p_pcg4_f_c_a19
            if (p_pcg5_f_c_a19 .ge. p1st) lptr_pcg5_f_c_aer(19,itype,ai_phase)  = p_pcg5_f_c_a19
            if (p_pcg6_f_c_a19 .ge. p1st) lptr_pcg6_f_c_aer(19,itype,ai_phase)  = p_pcg6_f_c_a19
            if (p_pcg7_f_c_a19 .ge. p1st) lptr_pcg7_f_c_aer(19,itype,ai_phase)  = p_pcg7_f_c_a19
            if (p_pcg8_f_c_a19 .ge. p1st) lptr_pcg8_f_c_aer(19,itype,ai_phase)  = p_pcg8_f_c_a19
            if (p_pcg9_f_c_a19 .ge. p1st) lptr_pcg9_f_c_aer(19,itype,ai_phase)  = p_pcg9_f_c_a19
            if (p_pcg1_f_o_a19 .ge. p1st) lptr_pcg1_f_o_aer(19,itype,ai_phase) = p_pcg1_f_o_a19
            if (p_pcg2_f_o_a19 .ge. p1st) lptr_pcg2_f_o_aer(19,itype,ai_phase)  = p_pcg2_f_o_a19
            if (p_pcg3_f_o_a19 .ge. p1st) lptr_pcg3_f_o_aer(19,itype,ai_phase)  = p_pcg3_f_o_a19
            if (p_pcg4_f_o_a19 .ge. p1st) lptr_pcg4_f_o_aer(19,itype,ai_phase)  = p_pcg4_f_o_a19
            if (p_pcg5_f_o_a19 .ge. p1st) lptr_pcg5_f_o_aer(19,itype,ai_phase)  = p_pcg5_f_o_a19
            if (p_pcg6_f_o_a19 .ge. p1st) lptr_pcg6_f_o_aer(19,itype,ai_phase)  = p_pcg6_f_o_a19
            if (p_pcg7_f_o_a19 .ge. p1st) lptr_pcg7_f_o_aer(19,itype,ai_phase)  = p_pcg7_f_o_a19
            if (p_pcg8_f_o_a19 .ge. p1st) lptr_pcg8_f_o_aer(19,itype,ai_phase)  = p_pcg8_f_o_a19
            if (p_pcg9_f_o_a19 .ge. p1st) lptr_pcg9_f_o_aer(19,itype,ai_phase)  = p_pcg9_f_o_a19
            if (p_opcg1_f_c_a19 .ge. p1st) lptr_opcg1_f_c_aer(19,itype,ai_phase) = p_opcg1_f_c_a19
            if (p_opcg2_f_c_a19 .ge. p1st) lptr_opcg2_f_c_aer(19,itype,ai_phase)  = p_opcg2_f_c_a19
            if (p_opcg3_f_c_a19 .ge. p1st) lptr_opcg3_f_c_aer(19,itype,ai_phase)  = p_opcg3_f_c_a19
            if (p_opcg4_f_c_a19 .ge. p1st) lptr_opcg4_f_c_aer(19,itype,ai_phase)  = p_opcg4_f_c_a19
            if (p_opcg5_f_c_a19 .ge. p1st) lptr_opcg5_f_c_aer(19,itype,ai_phase)  = p_opcg5_f_c_a19
            if (p_opcg6_f_c_a19 .ge. p1st) lptr_opcg6_f_c_aer(19,itype,ai_phase)  = p_opcg6_f_c_a19
            if (p_opcg7_f_c_a19 .ge. p1st) lptr_opcg7_f_c_aer(19,itype,ai_phase)  = p_opcg7_f_c_a19
            if (p_opcg8_f_c_a19 .ge. p1st) lptr_opcg8_f_c_aer(19,itype,ai_phase)  = p_opcg8_f_c_a19
            if (p_opcg1_f_o_a19 .ge. p1st) lptr_opcg1_f_o_aer(19,itype,ai_phase) = p_opcg1_f_o_a19
            if (p_opcg2_f_o_a19 .ge. p1st) lptr_opcg2_f_o_aer(19,itype,ai_phase)  = p_opcg2_f_o_a19
            if (p_opcg3_f_o_a19 .ge. p1st) lptr_opcg3_f_o_aer(19,itype,ai_phase)  = p_opcg3_f_o_a19
            if (p_opcg4_f_o_a19 .ge. p1st) lptr_opcg4_f_o_aer(19,itype,ai_phase)  = p_opcg4_f_o_a19
            if (p_opcg5_f_o_a19 .ge. p1st) lptr_opcg5_f_o_aer(19,itype,ai_phase)  = p_opcg5_f_o_a19
            if (p_opcg6_f_o_a19 .ge. p1st) lptr_opcg6_f_o_aer(19,itype,ai_phase)  = p_opcg6_f_o_a19
            if (p_opcg7_f_o_a19 .ge. p1st) lptr_opcg7_f_o_aer(19,itype,ai_phase)  = p_opcg7_f_o_a19
            if (p_opcg8_f_o_a19 .ge. p1st) lptr_opcg8_f_o_aer(19,itype,ai_phase)  = p_opcg8_f_o_a19
            if (p_ant1_c_a19 .ge. p1st) lptr_ant1_c_aer(19,itype,ai_phase)  = p_ant1_c_a19
            if (p_ant2_c_a19 .ge. p1st) lptr_ant2_c_aer(19,itype,ai_phase)  = p_ant2_c_a19
            if (p_ant3_c_a19 .ge. p1st) lptr_ant3_c_aer(19,itype,ai_phase)  = p_ant3_c_a19
            if (p_ant4_c_a19 .ge. p1st) lptr_ant4_c_aer(19,itype,ai_phase)  = p_ant4_c_a19
            if (p_biog1_c_a19 .ge. p1st) lptr_biog1_c_aer(19,itype,ai_phase)  = p_biog1_c_a19
            if (p_biog2_c_a19 .ge. p1st) lptr_biog2_c_aer(19,itype,ai_phase)  = p_biog2_c_a19
            if (p_biog3_c_a19 .ge. p1st) lptr_biog3_c_aer(19,itype,ai_phase)  = p_biog3_c_a19
            if (p_biog4_c_a19 .ge. p1st) lptr_biog4_c_aer(19,itype,ai_phase)  = p_biog4_c_a19
            if (p_ant1_o_a19 .ge. p1st) lptr_ant1_o_aer(19,itype,ai_phase)  = p_ant1_o_a19
            if (p_ant2_o_a19 .ge. p1st) lptr_ant2_o_aer(19,itype,ai_phase)  = p_ant2_o_a19
            if (p_ant3_o_a19 .ge. p1st) lptr_ant3_o_aer(19,itype,ai_phase)  = p_ant3_o_a19
            if (p_ant4_o_a19 .ge. p1st) lptr_ant4_o_aer(19,itype,ai_phase)  = p_ant4_o_a19
            if (p_biog1_o_a19 .ge. p1st) lptr_biog1_o_aer(19,itype,ai_phase)  = p_biog1_o_a19
            if (p_biog2_o_a19 .ge. p1st) lptr_biog2_o_aer(19,itype,ai_phase)  = p_biog2_o_a19
            if (p_biog3_o_a19 .ge. p1st) lptr_biog3_o_aer(19,itype,ai_phase)  = p_biog3_o_a19
            if (p_biog4_o_a19 .ge. p1st) lptr_biog4_o_aer(19,itype,ai_phase)  = p_biog4_o_a19
            
	    numptr_aer(19,itype,ai_phase)        = p_num_a19
	end if

	if (nsize_aer(itype) .ge. 20) then
	    lptr_so4_aer(20,itype,ai_phase)      = p_so4_a20
	    lptr_no3_aer(20,itype,ai_phase)      = p_no3_a20
	    lptr_cl_aer(20,itype,ai_phase)       = p_cl_a20
	    lptr_msa_aer(20,itype,ai_phase)      = p_msa_a20
	    lptr_co3_aer(20,itype,ai_phase)      = p_co3_a20
	    lptr_nh4_aer(20,itype,ai_phase)      = p_nh4_a20
	    lptr_na_aer(20,itype,ai_phase)       = p_na_a20
	    lptr_ca_aer(20,itype,ai_phase)       = p_ca_a20
	    lptr_oin_aer(20,itype,ai_phase)      = p_oin_a20
	    lptr_oc_aer(20,itype,ai_phase)       = p_oc_a20
	    lptr_bc_aer(20,itype,ai_phase)       = p_bc_a20
	    hyswptr_aer(20,itype)                = p_hysw_a20
	    waterptr_aer(20,itype)               = p_water_a20

            if (p_pcg1_b_c_a20 .ge. p1st) lptr_pcg1_b_c_aer(20,itype,ai_phase) = p_pcg1_b_c_a20
            if (p_pcg2_b_c_a20 .ge. p1st) lptr_pcg2_b_c_aer(20,itype,ai_phase)  = p_pcg2_b_c_a20
            if (p_pcg3_b_c_a20 .ge. p1st) lptr_pcg3_b_c_aer(20,itype,ai_phase)  = p_pcg3_b_c_a20
            if (p_pcg4_b_c_a20 .ge. p1st) lptr_pcg4_b_c_aer(20,itype,ai_phase)  = p_pcg4_b_c_a20
            if (p_pcg5_b_c_a20 .ge. p1st) lptr_pcg5_b_c_aer(20,itype,ai_phase)  = p_pcg5_b_c_a20
            if (p_pcg6_b_c_a20 .ge. p1st) lptr_pcg6_b_c_aer(20,itype,ai_phase)  = p_pcg6_b_c_a20
            if (p_pcg7_b_c_a20 .ge. p1st) lptr_pcg7_b_c_aer(20,itype,ai_phase)  = p_pcg7_b_c_a20
            if (p_pcg8_b_c_a20 .ge. p1st) lptr_pcg8_b_c_aer(20,itype,ai_phase)  = p_pcg8_b_c_a20
            if (p_pcg9_b_c_a20 .ge. p1st) lptr_pcg9_b_c_aer(20,itype,ai_phase)  = p_pcg9_b_c_a20
            if (p_pcg1_b_o_a20 .ge. p1st) lptr_pcg1_b_o_aer(20,itype,ai_phase) = p_pcg1_b_o_a20
            if (p_pcg2_b_o_a20 .ge. p1st) lptr_pcg2_b_o_aer(20,itype,ai_phase)  = p_pcg2_b_o_a20
            if (p_pcg3_b_o_a20 .ge. p1st) lptr_pcg3_b_o_aer(20,itype,ai_phase)  = p_pcg3_b_o_a20
            if (p_pcg4_b_o_a20 .ge. p1st) lptr_pcg4_b_o_aer(20,itype,ai_phase)  = p_pcg4_b_o_a20
            if (p_pcg5_b_o_a20 .ge. p1st) lptr_pcg5_b_o_aer(20,itype,ai_phase)  = p_pcg5_b_o_a20
            if (p_pcg6_b_o_a20 .ge. p1st) lptr_pcg6_b_o_aer(20,itype,ai_phase)  = p_pcg6_b_o_a20
            if (p_pcg7_b_o_a20 .ge. p1st) lptr_pcg7_b_o_aer(20,itype,ai_phase)  = p_pcg7_b_o_a20
            if (p_pcg8_b_o_a20 .ge. p1st) lptr_pcg8_b_o_aer(20,itype,ai_phase)  = p_pcg8_b_o_a20
            if (p_pcg9_b_o_a20 .ge. p1st) lptr_pcg9_b_o_aer(20,itype,ai_phase)  = p_pcg9_b_o_a20
            if (p_opcg1_b_c_a20 .ge. p1st) lptr_opcg1_b_c_aer(20,itype,ai_phase) = p_opcg1_b_c_a20
            if (p_opcg2_b_c_a20 .ge. p1st) lptr_opcg2_b_c_aer(20,itype,ai_phase)  = p_opcg2_b_c_a20
            if (p_opcg3_b_c_a20 .ge. p1st) lptr_opcg3_b_c_aer(20,itype,ai_phase)  = p_opcg3_b_c_a20
            if (p_opcg4_b_c_a20 .ge. p1st) lptr_opcg4_b_c_aer(20,itype,ai_phase)  = p_opcg4_b_c_a20
            if (p_opcg5_b_c_a20 .ge. p1st) lptr_opcg5_b_c_aer(20,itype,ai_phase)  = p_opcg5_b_c_a20
            if (p_opcg6_b_c_a20 .ge. p1st) lptr_opcg6_b_c_aer(20,itype,ai_phase)  = p_opcg6_b_c_a20
            if (p_opcg7_b_c_a20 .ge. p1st) lptr_opcg7_b_c_aer(20,itype,ai_phase)  = p_opcg7_b_c_a20
            if (p_opcg8_b_c_a20 .ge. p1st) lptr_opcg8_b_c_aer(20,itype,ai_phase)  = p_opcg8_b_c_a20
            if (p_opcg1_b_o_a20 .ge. p1st) lptr_opcg1_b_o_aer(20,itype,ai_phase) = p_opcg1_b_o_a20
            if (p_opcg2_b_o_a20 .ge. p1st) lptr_opcg2_b_o_aer(20,itype,ai_phase)  = p_opcg2_b_o_a20
            if (p_opcg3_b_o_a20 .ge. p1st) lptr_opcg3_b_o_aer(20,itype,ai_phase)  = p_opcg3_b_o_a20
            if (p_opcg4_b_o_a20 .ge. p1st) lptr_opcg4_b_o_aer(20,itype,ai_phase)  = p_opcg4_b_o_a20
            if (p_opcg5_b_o_a20 .ge. p1st) lptr_opcg5_b_o_aer(20,itype,ai_phase)  = p_opcg5_b_o_a20
            if (p_opcg6_b_o_a20 .ge. p1st) lptr_opcg6_b_o_aer(20,itype,ai_phase)  = p_opcg6_b_o_a20
            if (p_opcg7_b_o_a20 .ge. p1st) lptr_opcg7_b_o_aer(20,itype,ai_phase)  = p_opcg7_b_o_a20
            if (p_opcg8_b_o_a20 .ge. p1st) lptr_opcg8_b_o_aer(20,itype,ai_phase)  = p_opcg8_b_o_a20
            if (p_pcg1_f_c_a20 .ge. p1st) lptr_pcg1_f_c_aer(20,itype,ai_phase) = p_pcg1_f_c_a20
            if (p_pcg2_f_c_a20 .ge. p1st) lptr_pcg2_f_c_aer(20,itype,ai_phase)  = p_pcg2_f_c_a20
            if (p_pcg3_f_c_a20 .ge. p1st) lptr_pcg3_f_c_aer(20,itype,ai_phase)  = p_pcg3_f_c_a20
            if (p_pcg4_f_c_a20 .ge. p1st) lptr_pcg4_f_c_aer(20,itype,ai_phase)  = p_pcg4_f_c_a20
            if (p_pcg5_f_c_a20 .ge. p1st) lptr_pcg5_f_c_aer(20,itype,ai_phase)  = p_pcg5_f_c_a20
            if (p_pcg6_f_c_a20 .ge. p1st) lptr_pcg6_f_c_aer(20,itype,ai_phase)  = p_pcg6_f_c_a20
            if (p_pcg7_f_c_a20 .ge. p1st) lptr_pcg7_f_c_aer(20,itype,ai_phase)  = p_pcg7_f_c_a20
            if (p_pcg8_f_c_a20 .ge. p1st) lptr_pcg8_f_c_aer(20,itype,ai_phase)  = p_pcg8_f_c_a20
            if (p_pcg9_f_c_a20 .ge. p1st) lptr_pcg9_f_c_aer(20,itype,ai_phase)  = p_pcg9_f_c_a20
            if (p_pcg1_f_o_a20 .ge. p1st) lptr_pcg1_f_o_aer(20,itype,ai_phase) = p_pcg1_f_o_a20
            if (p_pcg2_f_o_a20 .ge. p1st) lptr_pcg2_f_o_aer(20,itype,ai_phase)  = p_pcg2_f_o_a20
            if (p_pcg3_f_o_a20 .ge. p1st) lptr_pcg3_f_o_aer(20,itype,ai_phase)  = p_pcg3_f_o_a20
            if (p_pcg4_f_o_a20 .ge. p1st) lptr_pcg4_f_o_aer(20,itype,ai_phase)  = p_pcg4_f_o_a20
            if (p_pcg5_f_o_a20 .ge. p1st) lptr_pcg5_f_o_aer(20,itype,ai_phase)  = p_pcg5_f_o_a20
            if (p_pcg6_f_o_a20 .ge. p1st) lptr_pcg6_f_o_aer(20,itype,ai_phase)  = p_pcg6_f_o_a20
            if (p_pcg7_f_o_a20 .ge. p1st) lptr_pcg7_f_o_aer(20,itype,ai_phase)  = p_pcg7_f_o_a20
            if (p_pcg8_f_o_a20 .ge. p1st) lptr_pcg8_f_o_aer(20,itype,ai_phase)  = p_pcg8_f_o_a20
            if (p_pcg9_f_o_a20 .ge. p1st) lptr_pcg9_f_o_aer(20,itype,ai_phase)  = p_pcg9_f_o_a20
            if (p_opcg1_f_c_a20 .ge. p1st) lptr_opcg1_f_c_aer(20,itype,ai_phase) = p_opcg1_f_c_a20
            if (p_opcg2_f_c_a20 .ge. p1st) lptr_opcg2_f_c_aer(20,itype,ai_phase)  = p_opcg2_f_c_a20
            if (p_opcg3_f_c_a20 .ge. p1st) lptr_opcg3_f_c_aer(20,itype,ai_phase)  = p_opcg3_f_c_a20
            if (p_opcg4_f_c_a20 .ge. p1st) lptr_opcg4_f_c_aer(20,itype,ai_phase)  = p_opcg4_f_c_a20
            if (p_opcg5_f_c_a20 .ge. p1st) lptr_opcg5_f_c_aer(20,itype,ai_phase)  = p_opcg5_f_c_a20
            if (p_opcg6_f_c_a20 .ge. p1st) lptr_opcg6_f_c_aer(20,itype,ai_phase)  = p_opcg6_f_c_a20
            if (p_opcg7_f_c_a20 .ge. p1st) lptr_opcg7_f_c_aer(20,itype,ai_phase)  = p_opcg7_f_c_a20
            if (p_opcg8_f_c_a20 .ge. p1st) lptr_opcg8_f_c_aer(20,itype,ai_phase)  = p_opcg8_f_c_a20
            if (p_opcg1_f_o_a20 .ge. p1st) lptr_opcg1_f_o_aer(20,itype,ai_phase) = p_opcg1_f_o_a20
            if (p_opcg2_f_o_a20 .ge. p1st) lptr_opcg2_f_o_aer(20,itype,ai_phase)  = p_opcg2_f_o_a20
            if (p_opcg3_f_o_a20 .ge. p1st) lptr_opcg3_f_o_aer(20,itype,ai_phase)  = p_opcg3_f_o_a20
            if (p_opcg4_f_o_a20 .ge. p1st) lptr_opcg4_f_o_aer(20,itype,ai_phase)  = p_opcg4_f_o_a20
            if (p_opcg5_f_o_a20 .ge. p1st) lptr_opcg5_f_o_aer(20,itype,ai_phase)  = p_opcg5_f_o_a20
            if (p_opcg6_f_o_a20 .ge. p1st) lptr_opcg6_f_o_aer(20,itype,ai_phase)  = p_opcg6_f_o_a20
            if (p_opcg7_f_o_a20 .ge. p1st) lptr_opcg7_f_o_aer(20,itype,ai_phase)  = p_opcg7_f_o_a20
            if (p_opcg8_f_o_a20 .ge. p1st) lptr_opcg8_f_o_aer(20,itype,ai_phase)  = p_opcg8_f_o_a20
            if (p_ant1_c_a20 .ge. p1st) lptr_ant1_c_aer(20,itype,ai_phase)  = p_ant1_c_a20
            if (p_ant2_c_a20 .ge. p1st) lptr_ant2_c_aer(20,itype,ai_phase)  = p_ant2_c_a20
            if (p_ant3_c_a20 .ge. p1st) lptr_ant3_c_aer(20,itype,ai_phase)  = p_ant3_c_a20
            if (p_ant4_c_a20 .ge. p1st) lptr_ant4_c_aer(20,itype,ai_phase)  = p_ant4_c_a20
            if (p_biog1_c_a20 .ge. p1st) lptr_biog1_c_aer(20,itype,ai_phase)  = p_biog1_c_a20
            if (p_biog2_c_a20 .ge. p1st) lptr_biog2_c_aer(20,itype,ai_phase)  = p_biog2_c_a20
            if (p_biog3_c_a20 .ge. p1st) lptr_biog3_c_aer(20,itype,ai_phase)  = p_biog3_c_a20
            if (p_biog4_c_a20 .ge. p1st) lptr_biog4_c_aer(20,itype,ai_phase)  = p_biog4_c_a20
            if (p_ant1_o_a20 .ge. p1st) lptr_ant1_o_aer(20,itype,ai_phase)  = p_ant1_o_a20
            if (p_ant2_o_a20 .ge. p1st) lptr_ant2_o_aer(20,itype,ai_phase)  = p_ant2_o_a20
            if (p_ant3_o_a20 .ge. p1st) lptr_ant3_o_aer(20,itype,ai_phase)  = p_ant3_o_a20
            if (p_ant4_o_a20 .ge. p1st) lptr_ant4_o_aer(20,itype,ai_phase)  = p_ant4_o_a20
            if (p_biog1_o_a20 .ge. p1st) lptr_biog1_o_aer(20,itype,ai_phase)  = p_biog1_o_a20
            if (p_biog2_o_a20 .ge. p1st) lptr_biog2_o_aer(20,itype,ai_phase)  = p_biog2_o_a20
            if (p_biog3_o_a20 .ge. p1st) lptr_biog3_o_aer(20,itype,ai_phase)  = p_biog3_o_a20
            if (p_biog4_o_a20 .ge. p1st) lptr_biog4_o_aer(20,itype,ai_phase)  = p_biog4_o_a20
            
	    numptr_aer(20,itype,ai_phase)        = p_num_a20
	end if


	if (nsize_aer(itype) .ge. 1) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(01,itype,cw_phase)      = p_so4_cw01
	    lptr_no3_aer(01,itype,cw_phase)      = p_no3_cw01
	    lptr_cl_aer(01,itype,cw_phase)       = p_cl_cw01
	    lptr_msa_aer(01,itype,cw_phase)      = p_msa_cw01
	    lptr_co3_aer(01,itype,cw_phase)      = p_co3_cw01
	    lptr_nh4_aer(01,itype,cw_phase)      = p_nh4_cw01
	    lptr_na_aer(01,itype,cw_phase)       = p_na_cw01
	    lptr_ca_aer(01,itype,cw_phase)       = p_ca_cw01
	    lptr_oin_aer(01,itype,cw_phase)      = p_oin_cw01
	    lptr_oc_aer(01,itype,cw_phase)       = p_oc_cw01
	    lptr_bc_aer(01,itype,cw_phase)       = p_bc_cw01
           if(p_pcg1_b_c_cw01.ge.p1st) lptr_pcg1_b_c_aer(01,itype,cw_phase) = p_pcg1_b_c_cw01
            if(p_opcg1_b_c_cw01.ge.p1st) lptr_opcg1_b_c_aer(01,itype,cw_phase) = p_opcg1_b_c_cw01
            if(p_pcg1_b_o_cw01.ge.p1st) lptr_pcg1_b_o_aer(01,itype,cw_phase) = p_pcg1_b_o_cw01
            if(p_opcg1_b_o_cw01.ge.p1st) lptr_opcg1_b_o_aer(01,itype,cw_phase) = p_opcg1_b_o_cw01
            if(p_pcg1_f_c_cw01.ge.p1st) lptr_pcg1_f_c_aer(01,itype,cw_phase) = p_pcg1_f_c_cw01
            if(p_opcg1_f_c_cw01.ge.p1st) lptr_opcg1_f_c_aer(01,itype,cw_phase) = p_opcg1_f_c_cw01
            if(p_pcg1_f_o_cw01.ge.p1st) lptr_pcg1_f_o_aer(01,itype,cw_phase) = p_pcg1_f_o_cw01
            if(p_opcg1_f_o_cw01.ge.p1st) lptr_opcg1_f_o_aer(01,itype,cw_phase) = p_opcg1_f_o_cw01
            if(p_ant1_c_cw01.ge.p1st) lptr_ant1_c_aer(01,itype,cw_phase) = p_ant1_c_cw01
            if(p_ant2_c_cw01.ge.p1st) lptr_ant2_c_aer(01,itype,cw_phase) = p_ant2_c_cw01
            if(p_ant3_c_cw01.ge.p1st) lptr_ant3_c_aer(01,itype,cw_phase) = p_ant3_c_cw01
            if(p_ant4_c_cw01.ge.p1st) lptr_ant4_c_aer(01,itype,cw_phase) = p_ant4_c_cw01
            if(p_ant2_o_cw01.ge.p1st) lptr_ant2_o_aer(01,itype,cw_phase) =p_ant2_o_cw01
            if(p_ant3_o_cw01.ge.p1st) lptr_ant3_o_aer(01,itype,cw_phase) =p_ant3_o_cw01
            if(p_biog1_c_cw01.ge.p1st) lptr_biog1_c_aer(01,itype,cw_phase) = p_biog1_c_cw01
            if(p_biog2_c_cw01.ge.p1st) lptr_biog2_c_aer(01,itype,cw_phase) = p_biog2_c_cw01
            if(p_biog3_c_cw01.ge.p1st) lptr_biog3_c_aer(01,itype,cw_phase) = p_biog3_c_cw01
            if(p_biog4_c_cw01.ge.p1st) lptr_biog4_c_aer(01,itype,cw_phase) = p_biog4_c_cw01
            if(p_iepox_cw01.ge.p1st) lptr_iepox_aer(01,itype,cw_phase) = p_iepox_cw01
            if(p_iepoxos_cw01.ge.p1st) lptr_iepoxos_aer(01,itype,cw_phase) = p_iepoxos_cw01
            if(p_tetrol_cw01.ge.p1st) lptr_tetrol_aer(01,itype,cw_phase) = p_tetrol_cw01
            if(p_cldiepoxos_cw01.ge.p1st) lptr_cldiepoxos_aer(01,itype,cw_phase) = p_cldiepoxos_cw01
            if(p_cldtetrol_cw01.ge.p1st) lptr_cldtetrol_aer(01,itype,cw_phase) = p_cldtetrol_cw01
            if(p_gly_cw01.ge.p1st) lptr_gly_aer(01,itype,cw_phase) = p_gly_cw01
            if(p_biog1_o_cw01.ge.p1st) lptr_biog1_o_aer(01,itype,cw_phase) = p_biog1_o_cw01
            if(p_biog2_o_cw01.ge.p1st) lptr_biog2_o_aer(01,itype,cw_phase) = p_biog2_o_cw01
            if(p_biog3_o_cw01.ge.p1st) lptr_biog3_o_aer(01,itype,cw_phase) = p_biog3_o_cw01
            if(p_biog4_o_cw01.ge.p1st) lptr_biog4_o_aer(01,itype,cw_phase) = p_biog4_o_cw01
            if(p_opcg1_b_o_cw01.ge.p1st) lptr_opcg1_b_o_aer(01,itype,cw_phase) = p_opcg1_b_o_cw01
            if(p_opcg2_b_o_cw01.ge.p1st) lptr_opcg2_b_o_aer(01,itype,cw_phase) = p_opcg2_b_o_cw01
            if(p_opcg3_b_o_cw01.ge.p1st) lptr_opcg3_b_o_aer(01,itype,cw_phase) = p_opcg3_b_o_cw01
            if(p_opcg4_b_o_cw01.ge.p1st) lptr_opcg4_b_o_aer(01,itype,cw_phase) = p_opcg4_b_o_cw01
            if(p_opcg5_b_o_cw01.ge.p1st) lptr_opcg5_b_o_aer(01,itype,cw_phase) = p_opcg5_b_o_cw01
            if(p_opcg6_b_o_cw01.ge.p1st) lptr_opcg6_b_o_aer(01,itype,cw_phase) = p_opcg6_b_o_cw01
            if(p_opcg7_b_o_cw01.ge.p1st) lptr_opcg7_b_o_aer(01,itype,cw_phase) = p_opcg7_b_o_cw01
            if(p_opcg8_b_o_cw01.ge.p1st) lptr_opcg8_b_o_aer(01,itype,cw_phase) = p_opcg8_b_o_cw01
            if(p_opcg2_f_o_cw01.ge.p1st) lptr_opcg2_f_o_aer(01,itype,cw_phase) = p_opcg2_f_o_cw01
            if(p_opcg3_f_o_cw01.ge.p1st) lptr_opcg3_f_o_aer(01,itype,cw_phase) = p_opcg3_f_o_cw01
            if(p_opcg4_f_o_cw01.ge.p1st) lptr_opcg4_f_o_aer(01,itype,cw_phase) = p_opcg4_f_o_cw01
            if(p_opcg5_f_o_cw01.ge.p1st) lptr_opcg5_f_o_aer(01,itype,cw_phase) = p_opcg5_f_o_cw01
            if(p_opcg6_f_o_cw01.ge.p1st) lptr_opcg6_f_o_aer(01,itype,cw_phase) = p_opcg6_f_o_cw01
            if(p_tanv_cw01.ge.p1st) lptr_tanv_aer(01,itype,cw_phase) = p_tanv_cw01
            if(p_sopnv_cw01.ge.p1st) lptr_sopnv_aer(01,itype,cw_phase) = p_sopnv_cw01
            if(p_ternv_cw01.ge.p1st) lptr_ternv_aer(01,itype,cw_phase) = p_ternv_cw01
            if(p_seqnv_cw01.ge.p1st) lptr_seqnv_aer(01,itype,cw_phase) = p_seqnv_cw01
            if(p_sianv_cw01.ge.p1st) lptr_sianv_aer(01,itype,cw_phase) = p_sianv_cw01

      if (p_glysoa_r1_cw01 .ge. p1st) lptr_glysoa_r1_aer(01,itype,cw_phase) = p_glysoa_r1_cw01
      if (p_glysoa_r2_cw01 .ge. p1st) lptr_glysoa_r2_aer(01,itype,cw_phase) = p_glysoa_r2_cw01
      if (p_glysoa_sfc_cw01 .ge. p1st) lptr_glysoa_sfc_aer(01,itype,cw_phase) = p_glysoa_sfc_cw01
      if (p_glysoa_nh4_cw01 .ge. p1st) lptr_glysoa_nh4_aer(01,itype,cw_phase) = p_glysoa_nh4_cw01
      if (p_glysoa_oh_cw01 .ge. p1st) lptr_glysoa_oh_aer(01,itype,cw_phase) = p_glysoa_oh_cw01
      if (p_asoaX_cw01 .ge. p1st)    lptr_asoaX_aer(01,itype,cw_phase)    = p_asoaX_cw01
      if (p_asoa1_cw01 .ge. p1st)    lptr_asoa1_aer(01,itype,cw_phase)    = p_asoa1_cw01
      if (p_asoa2_cw01 .ge. p1st)    lptr_asoa2_aer(01,itype,cw_phase)    = p_asoa2_cw01
      if (p_asoa3_cw01 .ge. p1st)    lptr_asoa3_aer(01,itype,cw_phase)    = p_asoa3_cw01
      if (p_asoa4_cw01 .ge. p1st)    lptr_asoa4_aer(01,itype,cw_phase)    = p_asoa4_cw01
      if (p_bsoaX_cw01 .ge. p1st)    lptr_bsoaX_aer(01,itype,cw_phase)    = p_bsoaX_cw01
      if (p_bsoa1_cw01 .ge. p1st)    lptr_bsoa1_aer(01,itype,cw_phase)    = p_bsoa1_cw01
      if (p_bsoa2_cw01 .ge. p1st)    lptr_bsoa2_aer(01,itype,cw_phase)    = p_bsoa2_cw01
      if (p_bsoa3_cw01 .ge. p1st)    lptr_bsoa3_aer(01,itype,cw_phase)    = p_bsoa3_cw01
      if (p_bsoa4_cw01 .ge. p1st)    lptr_bsoa4_aer(01,itype,cw_phase)    = p_bsoa4_cw01

	    numptr_aer(01,itype,cw_phase)        = p_num_cw01
	  
	  if(p_terp_cw01_Cn3.ge.p1st) lptr_terp_Cn3_aer(01,itype,cw_phase) = p_terp_cw01_Cn3

	  end if
	end if

	if (nsize_aer(itype) .ge. 2) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(02,itype,cw_phase)      = p_so4_cw02
	    lptr_no3_aer(02,itype,cw_phase)      = p_no3_cw02
	    lptr_cl_aer(02,itype,cw_phase)       = p_cl_cw02
	    lptr_msa_aer(02,itype,cw_phase)      = p_msa_cw02
	    lptr_co3_aer(02,itype,cw_phase)      = p_co3_cw02
	    lptr_nh4_aer(02,itype,cw_phase)      = p_nh4_cw02
	    lptr_na_aer(02,itype,cw_phase)       = p_na_cw02
	    lptr_ca_aer(02,itype,cw_phase)       = p_ca_cw02
	    lptr_oin_aer(02,itype,cw_phase)      = p_oin_cw02
	    lptr_oc_aer(02,itype,cw_phase)       = p_oc_cw02
	    lptr_bc_aer(02,itype,cw_phase)       = p_bc_cw02
            if(p_pcg1_b_c_cw02.ge.p1st) lptr_pcg1_b_c_aer(02,itype,cw_phase) = p_pcg1_b_c_cw02
            if(p_opcg1_b_c_cw02.ge.p1st) lptr_opcg1_b_c_aer(02,itype,cw_phase) = p_opcg1_b_c_cw02
            if(p_pcg1_b_o_cw02.ge.p1st) lptr_pcg1_b_o_aer(02,itype,cw_phase) = p_pcg1_b_o_cw02
            if(p_opcg1_b_o_cw02.ge.p1st) lptr_opcg1_b_o_aer(02,itype,cw_phase) = p_opcg1_b_o_cw02
            if(p_pcg1_f_c_cw02.ge.p1st) lptr_pcg1_f_c_aer(02,itype,cw_phase) = p_pcg1_f_c_cw02
            if(p_opcg1_f_c_cw02.ge.p1st) lptr_opcg1_f_c_aer(02,itype,cw_phase) = p_opcg1_f_c_cw02
            if(p_pcg1_f_o_cw02.ge.p1st) lptr_pcg1_f_o_aer(02,itype,cw_phase) = p_pcg1_f_o_cw02
            if(p_opcg1_f_o_cw02.ge.p1st) lptr_opcg1_f_o_aer(02,itype,cw_phase) = p_opcg1_f_o_cw02
            if(p_ant1_c_cw02.ge.p1st) lptr_ant1_c_aer(02,itype,cw_phase) = p_ant1_c_cw02
            if(p_ant2_c_cw02.ge.p1st) lptr_ant2_c_aer(02,itype,cw_phase) = p_ant2_c_cw02
            if(p_ant3_c_cw02.ge.p1st) lptr_ant3_c_aer(02,itype,cw_phase) = p_ant3_c_cw02
            if(p_ant4_c_cw02.ge.p1st) lptr_ant4_c_aer(02,itype,cw_phase) = p_ant4_c_cw02
            if(p_ant2_o_cw02.ge.p1st) lptr_ant2_o_aer(02,itype,cw_phase)  =p_ant2_o_cw02
            if(p_ant3_o_cw02.ge.p1st) lptr_ant3_o_aer(02,itype,cw_phase)  =p_ant3_o_cw02
            if(p_biog1_c_cw02.ge.p1st) lptr_biog1_c_aer(02,itype,cw_phase) = p_biog1_c_cw02
            if(p_biog2_c_cw02.ge.p1st) lptr_biog2_c_aer(02,itype,cw_phase) = p_biog2_c_cw02
            if(p_biog3_c_cw02.ge.p1st) lptr_biog3_c_aer(02,itype,cw_phase) = p_biog3_c_cw02
            if(p_biog4_c_cw02.ge.p1st) lptr_biog4_c_aer(02,itype,cw_phase) = p_biog4_c_cw02
            if(p_iepox_cw02.ge.p1st) lptr_iepox_aer(02,itype,cw_phase) = p_iepox_cw02
            if(p_iepoxos_cw02.ge.p1st) lptr_iepoxos_aer(02,itype,cw_phase) = p_iepoxos_cw02
            if(p_tetrol_cw02.ge.p1st) lptr_tetrol_aer(02,itype,cw_phase) = p_tetrol_cw02
            if(p_cldiepoxos_cw02.ge.p1st) lptr_cldiepoxos_aer(02,itype,cw_phase) = p_cldiepoxos_cw02
            if(p_cldtetrol_cw02.ge.p1st) lptr_cldtetrol_aer(02,itype,cw_phase) = p_cldtetrol_cw02
            if(p_gly_cw02.ge.p1st) lptr_gly_aer(02,itype,cw_phase) = p_gly_cw02
            if(p_biog1_o_cw02.ge.p1st) lptr_biog1_o_aer(02,itype,cw_phase) = p_biog1_o_cw02
            if(p_biog2_o_cw02.ge.p1st) lptr_biog2_o_aer(02,itype,cw_phase) = p_biog2_o_cw02
            if(p_biog3_o_cw02.ge.p1st) lptr_biog3_o_aer(02,itype,cw_phase) = p_biog3_o_cw02
            if(p_biog4_o_cw02.ge.p1st) lptr_biog4_o_aer(02,itype,cw_phase) = p_biog4_o_cw02
            if(p_opcg1_b_o_cw02.ge.p1st) lptr_opcg1_b_o_aer(02,itype,cw_phase) = p_opcg1_b_o_cw02
            if(p_opcg2_b_o_cw02.ge.p1st) lptr_opcg2_b_o_aer(02,itype,cw_phase) = p_opcg2_b_o_cw02
            if(p_opcg3_b_o_cw02.ge.p1st) lptr_opcg3_b_o_aer(02,itype,cw_phase) = p_opcg3_b_o_cw02
            if(p_opcg4_b_o_cw02.ge.p1st) lptr_opcg4_b_o_aer(02,itype,cw_phase) = p_opcg4_b_o_cw02
            if(p_opcg5_b_o_cw02.ge.p1st) lptr_opcg5_b_o_aer(02,itype,cw_phase) = p_opcg5_b_o_cw02
            if(p_opcg6_b_o_cw02.ge.p1st) lptr_opcg6_b_o_aer(02,itype,cw_phase) = p_opcg6_b_o_cw02
            if(p_opcg7_b_o_cw02.ge.p1st) lptr_opcg7_b_o_aer(02,itype,cw_phase) = p_opcg7_b_o_cw02
            if(p_opcg8_b_o_cw02.ge.p1st) lptr_opcg8_b_o_aer(02,itype,cw_phase) = p_opcg8_b_o_cw02
            if(p_opcg2_f_o_cw02.ge.p1st) lptr_opcg2_f_o_aer(02,itype,cw_phase) = p_opcg2_f_o_cw02
            if(p_opcg3_f_o_cw02.ge.p1st) lptr_opcg3_f_o_aer(02,itype,cw_phase) = p_opcg3_f_o_cw02
            if(p_opcg4_f_o_cw02.ge.p1st) lptr_opcg4_f_o_aer(02,itype,cw_phase) = p_opcg4_f_o_cw02
            if(p_opcg5_f_o_cw02.ge.p1st) lptr_opcg5_f_o_aer(02,itype,cw_phase) = p_opcg5_f_o_cw02
            if(p_opcg6_f_o_cw02.ge.p1st) lptr_opcg6_f_o_aer(02,itype,cw_phase) = p_opcg6_f_o_cw02

            if(p_tanv_cw02.ge.p1st) lptr_tanv_aer(02,itype,cw_phase) = p_tanv_cw02
            if(p_sopnv_cw02.ge.p1st) lptr_sopnv_aer(02,itype,cw_phase) = p_sopnv_cw02
            if(p_ternv_cw02.ge.p1st) lptr_ternv_aer(02,itype,cw_phase) = p_ternv_cw02
            if(p_seqnv_cw02.ge.p1st) lptr_seqnv_aer(02,itype,cw_phase) = p_seqnv_cw02
            if(p_sianv_cw02.ge.p1st) lptr_sianv_aer(02,itype,cw_phase) = p_sianv_cw02


      if (p_glysoa_r1_cw02 .ge. p1st) lptr_glysoa_r1_aer(02,itype,cw_phase) = p_glysoa_r1_cw02
      if (p_glysoa_r2_cw02 .ge. p1st) lptr_glysoa_r2_aer(02,itype,cw_phase) = p_glysoa_r2_cw02
      if (p_glysoa_sfc_cw02 .ge. p1st) lptr_glysoa_sfc_aer(02,itype,cw_phase) = p_glysoa_sfc_cw02
      if (p_glysoa_nh4_cw02 .ge. p1st) lptr_glysoa_nh4_aer(02,itype,cw_phase) = p_glysoa_nh4_cw02
      if (p_glysoa_oh_cw02 .ge. p1st) lptr_glysoa_oh_aer(02,itype,cw_phase) = p_glysoa_oh_cw02
      if (p_asoaX_cw02 .ge. p1st)    lptr_asoaX_aer(02,itype,cw_phase)    = p_asoaX_cw02
      if (p_asoa1_cw02 .ge. p1st)    lptr_asoa1_aer(02,itype,cw_phase)    = p_asoa1_cw02
      if (p_asoa2_cw02 .ge. p1st)    lptr_asoa2_aer(02,itype,cw_phase)    = p_asoa2_cw02
      if (p_asoa3_cw02 .ge. p1st)    lptr_asoa3_aer(02,itype,cw_phase)    = p_asoa3_cw02
      if (p_asoa4_cw02 .ge. p1st)    lptr_asoa4_aer(02,itype,cw_phase)    = p_asoa4_cw02
      if (p_bsoaX_cw02 .ge. p1st)    lptr_bsoaX_aer(02,itype,cw_phase)    = p_bsoaX_cw02
      if (p_bsoa1_cw02 .ge. p1st)    lptr_bsoa1_aer(02,itype,cw_phase)    = p_bsoa1_cw02
      if (p_bsoa2_cw02 .ge. p1st)    lptr_bsoa2_aer(02,itype,cw_phase)    = p_bsoa2_cw02
      if (p_bsoa3_cw02 .ge. p1st)    lptr_bsoa3_aer(02,itype,cw_phase)    = p_bsoa3_cw02
      if (p_bsoa4_cw02 .ge. p1st)    lptr_bsoa4_aer(02,itype,cw_phase)    = p_bsoa4_cw02
	    numptr_aer(02,itype,cw_phase)        = p_num_cw02
          
          if(p_terp_cw02_Cn3.ge.p1st) lptr_terp_Cn3_aer(02,itype,cw_phase) = p_terp_cw02_Cn3

	  end if
	end if

	if (nsize_aer(itype) .ge. 3) then
	  if (cw_phase .gt. 0) then

	    lptr_so4_aer(03,itype,cw_phase)      = p_so4_cw03
	    lptr_no3_aer(03,itype,cw_phase)      = p_no3_cw03
	    lptr_cl_aer(03,itype,cw_phase)       = p_cl_cw03
	    lptr_msa_aer(03,itype,cw_phase)      = p_msa_cw03
	    lptr_co3_aer(03,itype,cw_phase)      = p_co3_cw03
	    lptr_nh4_aer(03,itype,cw_phase)      = p_nh4_cw03
	    lptr_na_aer(03,itype,cw_phase)       = p_na_cw03
	    lptr_ca_aer(03,itype,cw_phase)       = p_ca_cw03
	    lptr_oin_aer(03,itype,cw_phase)      = p_oin_cw03
	    lptr_oc_aer(03,itype,cw_phase)       = p_oc_cw03
	    lptr_bc_aer(03,itype,cw_phase)       = p_bc_cw03
            if(p_pcg1_b_c_cw03.ge.p1st) lptr_pcg1_b_c_aer(03,itype,cw_phase) = p_pcg1_b_c_cw03
            if(p_opcg1_b_c_cw03.ge.p1st) lptr_opcg1_b_c_aer(03,itype,cw_phase) = p_opcg1_b_c_cw03
            if(p_pcg1_b_o_cw03.ge.p1st) lptr_pcg1_b_o_aer(03,itype,cw_phase) = p_pcg1_b_o_cw03
            if(p_opcg1_b_o_cw03.ge.p1st) lptr_opcg1_b_o_aer(03,itype,cw_phase) = p_opcg1_b_o_cw03
            if(p_pcg1_f_c_cw03.ge.p1st) lptr_pcg1_f_c_aer(03,itype,cw_phase) = p_pcg1_f_c_cw03
            if(p_opcg1_f_c_cw03.ge.p1st) lptr_opcg1_f_c_aer(03,itype,cw_phase) = p_opcg1_f_c_cw03
            if(p_pcg1_f_o_cw03.ge.p1st) lptr_pcg1_f_o_aer(03,itype,cw_phase) = p_pcg1_f_o_cw03
            if(p_opcg1_f_o_cw03.ge.p1st) lptr_opcg1_f_o_aer(03,itype,cw_phase) = p_opcg1_f_o_cw03
            if(p_ant1_c_cw03.ge.p1st) lptr_ant1_c_aer(03,itype,cw_phase) = p_ant1_c_cw03
            if(p_ant2_c_cw03.ge.p1st) lptr_ant2_c_aer(03,itype,cw_phase) = p_ant2_c_cw03
            if(p_ant3_c_cw03.ge.p1st) lptr_ant3_c_aer(03,itype,cw_phase) = p_ant3_c_cw03
            if(p_ant4_c_cw03.ge.p1st) lptr_ant4_c_aer(03,itype,cw_phase) = p_ant4_c_cw03
            if(p_ant2_o_cw03.ge.p1st) lptr_ant2_o_aer(03,itype,cw_phase)  =p_ant2_o_cw03
            if(p_ant3_o_cw03.ge.p1st) lptr_ant3_o_aer(03,itype,cw_phase)  =p_ant3_o_cw03
            if(p_biog1_c_cw03.ge.p1st) lptr_biog1_c_aer(03,itype,cw_phase) = p_biog1_c_cw03
            if(p_biog2_c_cw03.ge.p1st) lptr_biog2_c_aer(03,itype,cw_phase) = p_biog2_c_cw03
            if(p_biog3_c_cw03.ge.p1st) lptr_biog3_c_aer(03,itype,cw_phase) = p_biog3_c_cw03
            if(p_biog4_c_cw03.ge.p1st) lptr_biog4_c_aer(03,itype,cw_phase) = p_biog4_c_cw03
            if(p_iepox_cw03.ge.p1st) lptr_iepox_aer(03,itype,cw_phase) = p_iepox_cw03
            if(p_iepoxos_cw03.ge.p1st) lptr_iepoxos_aer(03,itype,cw_phase) = p_iepoxos_cw03
            if(p_tetrol_cw03.ge.p1st) lptr_tetrol_aer(03,itype,cw_phase) = p_tetrol_cw03
            if(p_cldiepoxos_cw03.ge.p1st) lptr_cldiepoxos_aer(03,itype,cw_phase) = p_cldiepoxos_cw03
            if(p_cldtetrol_cw03.ge.p1st) lptr_cldtetrol_aer(03,itype,cw_phase) = p_cldtetrol_cw03
            if(p_gly_cw03.ge.p1st) lptr_gly_aer(03,itype,cw_phase) = p_gly_cw03
            if(p_biog1_o_cw03.ge.p1st) lptr_biog1_o_aer(03,itype,cw_phase) = p_biog1_o_cw03
            if(p_biog2_o_cw03.ge.p1st) lptr_biog2_o_aer(03,itype,cw_phase) = p_biog2_o_cw03
            if(p_biog3_o_cw03.ge.p1st) lptr_biog3_o_aer(03,itype,cw_phase) = p_biog3_o_cw03
            if(p_biog4_o_cw03.ge.p1st) lptr_biog4_o_aer(03,itype,cw_phase) = p_biog4_o_cw03
            if(p_opcg1_b_o_cw03.ge.p1st) lptr_opcg1_b_o_aer(03,itype,cw_phase) = p_opcg1_b_o_cw03
            if(p_opcg2_b_o_cw03.ge.p1st) lptr_opcg2_b_o_aer(03,itype,cw_phase) = p_opcg2_b_o_cw03
            if(p_opcg3_b_o_cw03.ge.p1st) lptr_opcg3_b_o_aer(03,itype,cw_phase) = p_opcg3_b_o_cw03
            if(p_opcg4_b_o_cw03.ge.p1st) lptr_opcg4_b_o_aer(03,itype,cw_phase) = p_opcg4_b_o_cw03
            if(p_opcg5_b_o_cw03.ge.p1st) lptr_opcg5_b_o_aer(03,itype,cw_phase) = p_opcg5_b_o_cw03
            if(p_opcg6_b_o_cw03.ge.p1st) lptr_opcg6_b_o_aer(03,itype,cw_phase) = p_opcg6_b_o_cw03
            if(p_opcg7_b_o_cw03.ge.p1st) lptr_opcg7_b_o_aer(03,itype,cw_phase) = p_opcg7_b_o_cw03
            if(p_opcg8_b_o_cw03.ge.p1st) lptr_opcg8_b_o_aer(03,itype,cw_phase) = p_opcg8_b_o_cw03
            if(p_opcg2_f_o_cw03.ge.p1st) lptr_opcg2_f_o_aer(03,itype,cw_phase) = p_opcg2_f_o_cw03
            if(p_opcg3_f_o_cw03.ge.p1st) lptr_opcg3_f_o_aer(03,itype,cw_phase) = p_opcg3_f_o_cw03
            if(p_opcg4_f_o_cw03.ge.p1st) lptr_opcg4_f_o_aer(03,itype,cw_phase) = p_opcg4_f_o_cw03
            if(p_opcg5_f_o_cw03.ge.p1st) lptr_opcg5_f_o_aer(03,itype,cw_phase) = p_opcg5_f_o_cw03
            if(p_opcg6_f_o_cw03.ge.p1st) lptr_opcg6_f_o_aer(03,itype,cw_phase) = p_opcg6_f_o_cw03
            if(p_tanv_cw03.ge.p1st) lptr_tanv_aer(03,itype,cw_phase) = p_tanv_cw03
            if(p_sopnv_cw03.ge.p1st) lptr_sopnv_aer(03,itype,cw_phase) = p_sopnv_cw03
            if(p_ternv_cw03.ge.p1st) lptr_ternv_aer(03,itype,cw_phase) = p_ternv_cw03
            if(p_seqnv_cw03.ge.p1st) lptr_seqnv_aer(03,itype,cw_phase) = p_seqnv_cw03
            if(p_sianv_cw03.ge.p1st) lptr_sianv_aer(03,itype,cw_phase) = p_sianv_cw03
            

      if (p_glysoa_r1_cw03 .ge. p1st) lptr_glysoa_r1_aer(03,itype,cw_phase) = p_glysoa_r1_cw03
      if (p_glysoa_r2_cw03 .ge. p1st) lptr_glysoa_r2_aer(03,itype,cw_phase) = p_glysoa_r2_cw03
      if (p_glysoa_sfc_cw03 .ge. p1st) lptr_glysoa_sfc_aer(03,itype,cw_phase) = p_glysoa_sfc_cw03
      if (p_glysoa_nh4_cw03 .ge. p1st) lptr_glysoa_nh4_aer(03,itype,cw_phase) = p_glysoa_nh4_cw03
      if (p_glysoa_oh_cw03 .ge. p1st) lptr_glysoa_oh_aer(03,itype,cw_phase) = p_glysoa_oh_cw03
      if (p_asoaX_cw03 .ge. p1st)    lptr_asoaX_aer(03,itype,cw_phase)    = p_asoaX_cw03
      if (p_asoa1_cw03 .ge. p1st)    lptr_asoa1_aer(03,itype,cw_phase)    = p_asoa1_cw03
      if (p_asoa2_cw03 .ge. p1st)    lptr_asoa2_aer(03,itype,cw_phase)    = p_asoa2_cw03
      if (p_asoa3_cw03 .ge. p1st)    lptr_asoa3_aer(03,itype,cw_phase)    = p_asoa3_cw03
      if (p_asoa4_cw03 .ge. p1st)    lptr_asoa4_aer(03,itype,cw_phase)    = p_asoa4_cw03
      if (p_bsoaX_cw03 .ge. p1st)    lptr_bsoaX_aer(03,itype,cw_phase)    = p_bsoaX_cw03
      if (p_bsoa1_cw03 .ge. p1st)    lptr_bsoa1_aer(03,itype,cw_phase)    = p_bsoa1_cw03
      if (p_bsoa2_cw03 .ge. p1st)    lptr_bsoa2_aer(03,itype,cw_phase)    = p_bsoa2_cw03
      if (p_bsoa3_cw03 .ge. p1st)    lptr_bsoa3_aer(03,itype,cw_phase)    = p_bsoa3_cw03
      if (p_bsoa4_cw03 .ge. p1st)    lptr_bsoa4_aer(03,itype,cw_phase)    = p_bsoa4_cw03
	    numptr_aer(03,itype,cw_phase)        = p_num_cw03
	      
	      if(p_terp_cw03_Cn3.ge.p1st) lptr_terp_Cn3_aer(03,itype,cw_phase) = p_terp_cw03_Cn3

	  end if
	end if

	if (nsize_aer(itype) .ge. 4) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(04,itype,cw_phase)      = p_so4_cw04
	    lptr_no3_aer(04,itype,cw_phase)      = p_no3_cw04
	    lptr_cl_aer(04,itype,cw_phase)       = p_cl_cw04
	    lptr_msa_aer(04,itype,cw_phase)      = p_msa_cw04
	    lptr_co3_aer(04,itype,cw_phase)      = p_co3_cw04
	    lptr_nh4_aer(04,itype,cw_phase)      = p_nh4_cw04
	    lptr_na_aer(04,itype,cw_phase)       = p_na_cw04
	    lptr_ca_aer(04,itype,cw_phase)       = p_ca_cw04
	    lptr_oin_aer(04,itype,cw_phase)      = p_oin_cw04
	    lptr_oc_aer(04,itype,cw_phase)       = p_oc_cw04
	    lptr_bc_aer(04,itype,cw_phase)       = p_bc_cw04

            if(p_pcg1_b_c_cw04.ge.p1st) lptr_pcg1_b_c_aer(04,itype,cw_phase) = p_pcg1_b_c_cw04
            if(p_opcg1_b_c_cw04.ge.p1st) lptr_opcg1_b_c_aer(04,itype,cw_phase) = p_opcg1_b_c_cw04
            if(p_pcg1_b_o_cw04.ge.p1st) lptr_pcg1_b_o_aer(04,itype,cw_phase) = p_pcg1_b_o_cw04
            if(p_opcg1_b_o_cw04.ge.p1st) lptr_opcg1_b_o_aer(04,itype,cw_phase) = p_opcg1_b_o_cw04
            if(p_pcg1_f_c_cw04.ge.p1st) lptr_pcg1_f_c_aer(04,itype,cw_phase) = p_pcg1_f_c_cw04
            if(p_opcg1_f_c_cw04.ge.p1st) lptr_opcg1_f_c_aer(04,itype,cw_phase) = p_opcg1_f_c_cw04
            if(p_pcg1_f_o_cw04.ge.p1st) lptr_pcg1_f_o_aer(04,itype,cw_phase) = p_pcg1_f_o_cw04
            if(p_opcg1_f_o_cw04.ge.p1st) lptr_opcg1_f_o_aer(04,itype,cw_phase) = p_opcg1_f_o_cw04
            if(p_ant1_c_cw04.ge.p1st) lptr_ant1_c_aer(04,itype,cw_phase) = p_ant1_c_cw04
            if(p_ant2_c_cw04.ge.p1st) lptr_ant2_c_aer(04,itype,cw_phase) = p_ant2_c_cw04
            if(p_ant3_c_cw04.ge.p1st) lptr_ant3_c_aer(04,itype,cw_phase) = p_ant3_c_cw04
            if(p_ant4_c_cw04.ge.p1st) lptr_ant4_c_aer(04,itype,cw_phase) = p_ant4_c_cw04
            if(p_ant2_o_cw04.ge.p1st) lptr_ant2_o_aer(04,itype,cw_phase)   =p_ant2_o_cw04
            if(p_ant3_o_cw04.ge.p1st) lptr_ant3_o_aer(04,itype,cw_phase)   =p_ant3_o_cw04
            if(p_biog1_c_cw04.ge.p1st) lptr_biog1_c_aer(04,itype,cw_phase) = p_biog1_c_cw04
            if(p_biog2_c_cw04.ge.p1st) lptr_biog2_c_aer(04,itype,cw_phase) = p_biog2_c_cw04
            if(p_biog3_c_cw04.ge.p1st) lptr_biog3_c_aer(04,itype,cw_phase) = p_biog3_c_cw04
            if(p_biog4_c_cw04.ge.p1st) lptr_biog4_c_aer(04,itype,cw_phase) = p_biog4_c_cw04
            if(p_iepox_cw04.ge.p1st) lptr_iepox_aer(04,itype,cw_phase) = p_iepox_cw04
            if(p_iepoxos_cw04.ge.p1st) lptr_iepoxos_aer(04,itype,cw_phase) = p_iepoxos_cw04
            if(p_tetrol_cw04.ge.p1st) lptr_tetrol_aer(04,itype,cw_phase) = p_tetrol_cw04
            if(p_cldiepoxos_cw04.ge.p1st) lptr_cldiepoxos_aer(04,itype,cw_phase) = p_cldiepoxos_cw04
            if(p_cldtetrol_cw04.ge.p1st) lptr_cldtetrol_aer(04,itype,cw_phase) = p_cldtetrol_cw04
            if(p_gly_cw04.ge.p1st) lptr_gly_aer(04,itype,cw_phase) = p_gly_cw04
            if(p_biog1_o_cw04.ge.p1st) lptr_biog1_o_aer(04,itype,cw_phase) = p_biog1_o_cw04
            if(p_biog2_o_cw04.ge.p1st) lptr_biog2_o_aer(04,itype,cw_phase) = p_biog2_o_cw04
            if(p_biog3_o_cw04.ge.p1st) lptr_biog3_o_aer(04,itype,cw_phase) = p_biog3_o_cw04
            if(p_biog4_o_cw04.ge.p1st) lptr_biog4_o_aer(04,itype,cw_phase) = p_biog4_o_cw04
            if(p_opcg1_b_o_cw04.ge.p1st) lptr_opcg1_b_o_aer(04,itype,cw_phase) = p_opcg1_b_o_cw04
            if(p_opcg2_b_o_cw04.ge.p1st) lptr_opcg2_b_o_aer(04,itype,cw_phase) = p_opcg2_b_o_cw04
            if(p_opcg3_b_o_cw04.ge.p1st) lptr_opcg3_b_o_aer(04,itype,cw_phase) = p_opcg3_b_o_cw04
            if(p_opcg4_b_o_cw04.ge.p1st) lptr_opcg4_b_o_aer(04,itype,cw_phase) = p_opcg4_b_o_cw04
            if(p_opcg5_b_o_cw04.ge.p1st) lptr_opcg5_b_o_aer(04,itype,cw_phase) = p_opcg5_b_o_cw04
            if(p_opcg6_b_o_cw04.ge.p1st) lptr_opcg6_b_o_aer(04,itype,cw_phase) = p_opcg6_b_o_cw04
            if(p_opcg7_b_o_cw04.ge.p1st) lptr_opcg7_b_o_aer(04,itype,cw_phase) = p_opcg7_b_o_cw04
            if(p_opcg8_b_o_cw04.ge.p1st) lptr_opcg8_b_o_aer(04,itype,cw_phase) = p_opcg8_b_o_cw04
            if(p_opcg2_f_o_cw04.ge.p1st) lptr_opcg2_f_o_aer(04,itype,cw_phase) = p_opcg2_f_o_cw04
            if(p_opcg3_f_o_cw04.ge.p1st) lptr_opcg3_f_o_aer(04,itype,cw_phase) = p_opcg3_f_o_cw04
            if(p_opcg4_f_o_cw04.ge.p1st) lptr_opcg4_f_o_aer(04,itype,cw_phase) = p_opcg4_f_o_cw04
            if(p_opcg5_f_o_cw04.ge.p1st) lptr_opcg5_f_o_aer(04,itype,cw_phase) = p_opcg5_f_o_cw04
            if(p_opcg6_f_o_cw04.ge.p1st) lptr_opcg6_f_o_aer(04,itype,cw_phase) = p_opcg6_f_o_cw04
            if(p_tanv_cw04.ge.p1st) lptr_tanv_aer(04,itype,cw_phase) = p_tanv_cw04
            if(p_sopnv_cw04.ge.p1st) lptr_sopnv_aer(04,itype,cw_phase) = p_sopnv_cw04
            if(p_ternv_cw04.ge.p1st) lptr_ternv_aer(04,itype,cw_phase) = p_ternv_cw04
            if(p_seqnv_cw04.ge.p1st) lptr_seqnv_aer(04,itype,cw_phase) = p_seqnv_cw04
            if(p_sianv_cw04.ge.p1st) lptr_sianv_aer(04,itype,cw_phase) = p_sianv_cw04


      if (p_glysoa_r1_cw04 .ge. p1st) lptr_glysoa_r1_aer(04,itype,cw_phase) = p_glysoa_r1_cw04
      if (p_glysoa_r2_cw04 .ge. p1st) lptr_glysoa_r2_aer(04,itype,cw_phase) = p_glysoa_r2_cw04
      if (p_glysoa_sfc_cw04 .ge. p1st) lptr_glysoa_sfc_aer(04,itype,cw_phase) = p_glysoa_sfc_cw04
      if (p_glysoa_nh4_cw04 .ge. p1st) lptr_glysoa_nh4_aer(04,itype,cw_phase) = p_glysoa_nh4_cw04
      if (p_glysoa_oh_cw04 .ge. p1st) lptr_glysoa_oh_aer(04,itype,cw_phase) = p_glysoa_oh_cw04
      if (p_asoaX_cw04 .ge. p1st)    lptr_asoaX_aer(04,itype,cw_phase)    = p_asoaX_cw04
      if (p_asoa1_cw04 .ge. p1st)    lptr_asoa1_aer(04,itype,cw_phase)    = p_asoa1_cw04
      if (p_asoa2_cw04 .ge. p1st)    lptr_asoa2_aer(04,itype,cw_phase)    = p_asoa2_cw04
      if (p_asoa3_cw04 .ge. p1st)    lptr_asoa3_aer(04,itype,cw_phase)    = p_asoa3_cw04
      if (p_asoa4_cw04 .ge. p1st)    lptr_asoa4_aer(04,itype,cw_phase)    = p_asoa4_cw04
      if (p_bsoaX_cw04 .ge. p1st)    lptr_bsoaX_aer(04,itype,cw_phase)    = p_bsoaX_cw04
      if (p_bsoa1_cw04 .ge. p1st)    lptr_bsoa1_aer(04,itype,cw_phase)    = p_bsoa1_cw04
      if (p_bsoa2_cw04 .ge. p1st)    lptr_bsoa2_aer(04,itype,cw_phase)    = p_bsoa2_cw04
      if (p_bsoa3_cw04 .ge. p1st)    lptr_bsoa3_aer(04,itype,cw_phase)    = p_bsoa3_cw04
      if (p_bsoa4_cw04 .ge. p1st)    lptr_bsoa4_aer(04,itype,cw_phase)    = p_bsoa4_cw04
	    numptr_aer(04,itype,cw_phase)        = p_num_cw04
	    if(p_terp_cw04_Cn3.ge.p1st) lptr_terp_Cn3_aer(04,itype,cw_phase) = p_terp_cw04_Cn3
	  end if
	end if

	if (nsize_aer(itype) .ge. 5) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(05,itype,cw_phase)      = p_so4_cw05
	    lptr_no3_aer(05,itype,cw_phase)      = p_no3_cw05
	    lptr_cl_aer(05,itype,cw_phase)       = p_cl_cw05
	    lptr_msa_aer(05,itype,cw_phase)      = p_msa_cw05
	    lptr_co3_aer(05,itype,cw_phase)      = p_co3_cw05
	    lptr_nh4_aer(05,itype,cw_phase)      = p_nh4_cw05
	    lptr_na_aer(05,itype,cw_phase)       = p_na_cw05
	    lptr_ca_aer(05,itype,cw_phase)       = p_ca_cw05
	    lptr_oin_aer(05,itype,cw_phase)      = p_oin_cw05
	    lptr_oc_aer(05,itype,cw_phase)       = p_oc_cw05
	    lptr_bc_aer(05,itype,cw_phase)       = p_bc_cw05
          if(p_pcg1_b_c_cw05.ge.p1st) lptr_pcg1_b_c_aer(05,itype,cw_phase) = p_pcg1_b_c_cw05
            if(p_opcg1_b_c_cw05.ge.p1st) lptr_opcg1_b_c_aer(05,itype,cw_phase) = p_opcg1_b_c_cw05
            if(p_pcg1_b_o_cw05.ge.p1st) lptr_pcg1_b_o_aer(05,itype,cw_phase) = p_pcg1_b_o_cw05
            if(p_opcg1_b_o_cw05.ge.p1st) lptr_opcg1_b_o_aer(05,itype,cw_phase) = p_opcg1_b_o_cw05
            if(p_pcg1_f_c_cw05.ge.p1st) lptr_pcg1_f_c_aer(05,itype,cw_phase) = p_pcg1_f_c_cw05
            if(p_opcg1_f_c_cw05.ge.p1st) lptr_opcg1_f_c_aer(05,itype,cw_phase) = p_opcg1_f_c_cw05
            if(p_pcg1_f_o_cw05.ge.p1st) lptr_pcg1_f_o_aer(05,itype,cw_phase) = p_pcg1_f_o_cw05
            if(p_opcg1_f_o_cw05.ge.p1st) lptr_opcg1_f_o_aer(05,itype,cw_phase) = p_opcg1_f_o_cw05
            if(p_ant1_c_cw05.ge.p1st) lptr_ant1_c_aer(05,itype,cw_phase) = p_ant1_c_cw05
            if(p_ant2_c_cw05.ge.p1st) lptr_ant2_c_aer(05,itype,cw_phase) = p_ant2_c_cw05
            if(p_ant3_c_cw05.ge.p1st) lptr_ant3_c_aer(05,itype,cw_phase) = p_ant3_c_cw05
            if(p_ant4_c_cw05.ge.p1st) lptr_ant4_c_aer(05,itype,cw_phase) = p_ant4_c_cw05
            if(p_biog1_c_cw05.ge.p1st) lptr_biog1_c_aer(05,itype,cw_phase) = p_biog1_c_cw05
            if(p_biog2_c_cw05.ge.p1st) lptr_biog2_c_aer(05,itype,cw_phase) = p_biog2_c_cw05
            if(p_biog3_c_cw05.ge.p1st) lptr_biog3_c_aer(05,itype,cw_phase) = p_biog3_c_cw05
            if(p_biog4_c_cw05.ge.p1st) lptr_biog4_c_aer(05,itype,cw_phase) = p_biog4_c_cw05
            if(p_iepox_cw05.ge.p1st) lptr_iepox_aer(05,itype,cw_phase) = p_iepox_cw05
            if(p_iepoxos_cw05.ge.p1st) lptr_iepoxos_aer(05,itype,cw_phase) = p_iepoxos_cw05
            if(p_tetrol_cw05.ge.p1st) lptr_tetrol_aer(05,itype,cw_phase) = p_tetrol_cw05
            if(p_gly_cw05.ge.p1st) lptr_gly_aer(05,itype,cw_phase) = p_gly_cw05
            if(p_biog1_o_cw05.ge.p1st) lptr_biog1_o_aer(05,itype,cw_phase) = p_biog1_o_cw05
            if(p_biog2_o_cw05.ge.p1st) lptr_biog2_o_aer(05,itype,cw_phase) = p_biog2_o_cw05
            if(p_biog3_o_cw05.ge.p1st) lptr_biog3_o_aer(05,itype,cw_phase) = p_biog3_o_cw05
            if(p_biog4_o_cw05.ge.p1st) lptr_biog4_o_aer(05,itype,cw_phase) = p_biog4_o_cw05
            if(p_opcg2_b_o_cw05.ge.p1st) lptr_opcg2_b_o_aer(05,itype,cw_phase) = p_opcg2_b_o_cw05
            if(p_opcg3_b_o_cw05.ge.p1st) lptr_opcg3_b_o_aer(05,itype,cw_phase) = p_opcg3_b_o_cw05
            if(p_opcg4_b_o_cw05.ge.p1st) lptr_opcg4_b_o_aer(05,itype,cw_phase) = p_opcg4_b_o_cw05
            if(p_opcg5_b_o_cw05.ge.p1st) lptr_opcg5_b_o_aer(05,itype,cw_phase) = p_opcg5_b_o_cw05
            if(p_opcg6_b_o_cw05.ge.p1st) lptr_opcg6_b_o_aer(05,itype,cw_phase) = p_opcg6_b_o_cw05
            if(p_tanv_cw05.ge.p1st) lptr_tanv_aer(05,itype,cw_phase) = p_tanv_cw05
            if(p_sopnv_cw05.ge.p1st) lptr_sopnv_aer(05,itype,cw_phase) = p_sopnv_cw05
            if(p_ternv_cw05.ge.p1st) lptr_ternv_aer(05,itype,cw_phase) = p_ternv_cw05
            if(p_seqnv_cw05.ge.p1st) lptr_seqnv_aer(05,itype,cw_phase) = p_seqnv_cw05
            if(p_sianv_cw05.ge.p1st) lptr_sianv_aer(05,itype,cw_phase) = p_sianv_cw05
            

	    numptr_aer(05,itype,cw_phase)        = p_num_cw05
	  end if
	end if

	if (nsize_aer(itype) .ge. 6) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(06,itype,cw_phase)      = p_so4_cw06
	    lptr_no3_aer(06,itype,cw_phase)      = p_no3_cw06
	    lptr_cl_aer(06,itype,cw_phase)       = p_cl_cw06
	    lptr_msa_aer(06,itype,cw_phase)      = p_msa_cw06
	    lptr_co3_aer(06,itype,cw_phase)      = p_co3_cw06
	    lptr_nh4_aer(06,itype,cw_phase)      = p_nh4_cw06
	    lptr_na_aer(06,itype,cw_phase)       = p_na_cw06
	    lptr_ca_aer(06,itype,cw_phase)       = p_ca_cw06
	    lptr_oin_aer(06,itype,cw_phase)      = p_oin_cw06
	    lptr_oc_aer(06,itype,cw_phase)       = p_oc_cw06
	    lptr_bc_aer(06,itype,cw_phase)       = p_bc_cw06
            if(p_pcg1_b_c_cw06.ge.p1st) lptr_pcg1_b_c_aer(06,itype,cw_phase) = p_pcg1_b_c_cw06
            if(p_opcg1_b_c_cw06.ge.p1st) lptr_opcg1_b_c_aer(06,itype,cw_phase) = p_opcg1_b_c_cw06
            if(p_pcg1_b_o_cw06.ge.p1st) lptr_pcg1_b_o_aer(06,itype,cw_phase) = p_pcg1_b_o_cw06
            if(p_opcg1_b_o_cw06.ge.p1st) lptr_opcg1_b_o_aer(06,itype,cw_phase) = p_opcg1_b_o_cw06
            if(p_pcg1_f_c_cw06.ge.p1st) lptr_pcg1_f_c_aer(06,itype,cw_phase) = p_pcg1_f_c_cw06
            if(p_opcg1_f_c_cw06.ge.p1st) lptr_opcg1_f_c_aer(06,itype,cw_phase) = p_opcg1_f_c_cw06
            if(p_pcg1_f_o_cw06.ge.p1st) lptr_pcg1_f_o_aer(06,itype,cw_phase) = p_pcg1_f_o_cw06
            if(p_opcg1_f_o_cw06.ge.p1st) lptr_opcg1_f_o_aer(06,itype,cw_phase) = p_opcg1_f_o_cw06
            if(p_ant1_c_cw06.ge.p1st) lptr_ant1_c_aer(06,itype,cw_phase) = p_ant1_c_cw06
            if(p_ant2_c_cw06.ge.p1st) lptr_ant2_c_aer(06,itype,cw_phase) = p_ant2_c_cw06
            if(p_ant3_c_cw06.ge.p1st) lptr_ant3_c_aer(06,itype,cw_phase) = p_ant3_c_cw06
            if(p_ant4_c_cw06.ge.p1st) lptr_ant4_c_aer(06,itype,cw_phase) = p_ant4_c_cw06
            if(p_biog1_c_cw06.ge.p1st) lptr_biog1_c_aer(06,itype,cw_phase) = p_biog1_c_cw06
            if(p_biog2_c_cw06.ge.p1st) lptr_biog2_c_aer(06,itype,cw_phase) = p_biog2_c_cw06
            if(p_biog3_c_cw06.ge.p1st) lptr_biog3_c_aer(06,itype,cw_phase) = p_biog3_c_cw06
            if(p_biog4_c_cw06.ge.p1st) lptr_biog4_c_aer(06,itype,cw_phase) = p_biog4_c_cw06
            if(p_iepox_cw06.ge.p1st) lptr_iepox_aer(06,itype,cw_phase) = p_iepox_cw06
            if(p_iepoxos_cw06.ge.p1st) lptr_iepoxos_aer(06,itype,cw_phase) = p_iepoxos_cw06
            if(p_tetrol_cw06.ge.p1st) lptr_tetrol_aer(06,itype,cw_phase) = p_tetrol_cw06
            if(p_gly_cw06.ge.p1st) lptr_gly_aer(06,itype,cw_phase) = p_gly_cw06
            if(p_biog1_o_cw06.ge.p1st) lptr_biog1_o_aer(06,itype,cw_phase) = p_biog1_o_cw06
            if(p_biog2_o_cw06.ge.p1st) lptr_biog2_o_aer(06,itype,cw_phase) = p_biog2_o_cw06
            if(p_biog3_o_cw06.ge.p1st) lptr_biog3_o_aer(06,itype,cw_phase) = p_biog3_o_cw06
            if(p_biog4_o_cw06.ge.p1st) lptr_biog4_o_aer(06,itype,cw_phase) = p_biog4_o_cw06
            if(p_opcg2_b_o_cw06.ge.p1st) lptr_opcg2_b_o_aer(06,itype,cw_phase) = p_opcg2_b_o_cw06
            if(p_opcg3_b_o_cw06.ge.p1st) lptr_opcg3_b_o_aer(06,itype,cw_phase) = p_opcg3_b_o_cw06
            if(p_opcg4_b_o_cw06.ge.p1st) lptr_opcg4_b_o_aer(06,itype,cw_phase) = p_opcg4_b_o_cw06
            if(p_opcg5_b_o_cw06.ge.p1st) lptr_opcg5_b_o_aer(06,itype,cw_phase) = p_opcg5_b_o_cw06
            if(p_opcg6_b_o_cw06.ge.p1st) lptr_opcg6_b_o_aer(06,itype,cw_phase) = p_opcg6_b_o_cw06
            if(p_tanv_cw06.ge.p1st) lptr_tanv_aer(06,itype,cw_phase) = p_tanv_cw06
            if(p_sopnv_cw06.ge.p1st) lptr_sopnv_aer(06,itype,cw_phase) = p_sopnv_cw06
            if(p_ternv_cw06.ge.p1st) lptr_ternv_aer(06,itype,cw_phase) = p_ternv_cw06
            if(p_seqnv_cw06.ge.p1st) lptr_seqnv_aer(06,itype,cw_phase) = p_seqnv_cw06
            if(p_sianv_cw06.ge.p1st) lptr_sianv_aer(06,itype,cw_phase) = p_sianv_cw06

            

	    numptr_aer(06,itype,cw_phase)        = p_num_cw06
	  end if
	end if

	if (nsize_aer(itype) .ge. 7) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(07,itype,cw_phase)      = p_so4_cw07
	    lptr_no3_aer(07,itype,cw_phase)      = p_no3_cw07
	    lptr_cl_aer(07,itype,cw_phase)       = p_cl_cw07
	    lptr_msa_aer(07,itype,cw_phase)      = p_msa_cw07
	    lptr_co3_aer(07,itype,cw_phase)      = p_co3_cw07
	    lptr_nh4_aer(07,itype,cw_phase)      = p_nh4_cw07
	    lptr_na_aer(07,itype,cw_phase)       = p_na_cw07
	    lptr_ca_aer(07,itype,cw_phase)       = p_ca_cw07
	    lptr_oin_aer(07,itype,cw_phase)      = p_oin_cw07
	    lptr_oc_aer(07,itype,cw_phase)       = p_oc_cw07
	    lptr_bc_aer(07,itype,cw_phase)       = p_bc_cw07
        if(p_pcg1_b_c_cw07.ge.p1st) lptr_pcg1_b_c_aer(07,itype,cw_phase) = p_pcg1_b_c_cw07
            if(p_opcg1_b_c_cw07.ge.p1st) lptr_opcg1_b_c_aer(07,itype,cw_phase) = p_opcg1_b_c_cw07
            if(p_pcg1_b_o_cw07.ge.p1st) lptr_pcg1_b_o_aer(07,itype,cw_phase) = p_pcg1_b_o_cw07
            if(p_opcg1_b_o_cw07.ge.p1st) lptr_opcg1_b_o_aer(07,itype,cw_phase) = p_opcg1_b_o_cw07
            if(p_pcg1_f_c_cw07.ge.p1st) lptr_pcg1_f_c_aer(07,itype,cw_phase) = p_pcg1_f_c_cw07
            if(p_opcg1_f_c_cw07.ge.p1st) lptr_opcg1_f_c_aer(07,itype,cw_phase) = p_opcg1_f_c_cw07
            if(p_pcg1_f_o_cw07.ge.p1st) lptr_pcg1_f_o_aer(07,itype,cw_phase) = p_pcg1_f_o_cw07
            if(p_opcg1_f_o_cw07.ge.p1st) lptr_opcg1_f_o_aer(07,itype,cw_phase) = p_opcg1_f_o_cw07
            if(p_ant1_c_cw07.ge.p1st) lptr_ant1_c_aer(07,itype,cw_phase) = p_ant1_c_cw07
            if(p_ant2_c_cw07.ge.p1st) lptr_ant2_c_aer(07,itype,cw_phase) = p_ant2_c_cw07
            if(p_ant3_c_cw07.ge.p1st) lptr_ant3_c_aer(07,itype,cw_phase) = p_ant3_c_cw07
            if(p_ant4_c_cw07.ge.p1st) lptr_ant4_c_aer(07,itype,cw_phase) = p_ant4_c_cw07
            if(p_biog1_c_cw07.ge.p1st) lptr_biog1_c_aer(07,itype,cw_phase) = p_biog1_c_cw07
            if(p_biog2_c_cw07.ge.p1st) lptr_biog2_c_aer(07,itype,cw_phase) = p_biog2_c_cw07
            if(p_biog3_c_cw07.ge.p1st) lptr_biog3_c_aer(07,itype,cw_phase) = p_biog3_c_cw07
            if(p_biog4_c_cw07.ge.p1st) lptr_biog4_c_aer(07,itype,cw_phase) = p_biog4_c_cw07
            if(p_iepox_cw07.ge.p1st) lptr_iepox_aer(07,itype,cw_phase) = p_iepox_cw07
            if(p_iepoxos_cw07.ge.p1st) lptr_iepoxos_aer(07,itype,cw_phase) = p_iepoxos_cw07
            if(p_tetrol_cw07.ge.p1st) lptr_tetrol_aer(07,itype,cw_phase) = p_tetrol_cw07
            if(p_gly_cw07.ge.p1st) lptr_gly_aer(07,itype,cw_phase) = p_gly_cw07
            if(p_biog1_o_cw07.ge.p1st) lptr_biog1_o_aer(07,itype,cw_phase) = p_biog1_o_cw07
            if(p_biog2_o_cw07.ge.p1st) lptr_biog2_o_aer(07,itype,cw_phase) = p_biog2_o_cw07
            if(p_biog3_o_cw07.ge.p1st) lptr_biog3_o_aer(07,itype,cw_phase) = p_biog3_o_cw07
            if(p_biog4_o_cw07.ge.p1st) lptr_biog4_o_aer(07,itype,cw_phase) = p_biog4_o_cw07
            if(p_opcg2_b_o_cw07.ge.p1st) lptr_opcg2_b_o_aer(07,itype,cw_phase) = p_opcg2_b_o_cw07
            if(p_opcg3_b_o_cw07.ge.p1st) lptr_opcg3_b_o_aer(07,itype,cw_phase) = p_opcg3_b_o_cw07
            if(p_opcg4_b_o_cw07.ge.p1st) lptr_opcg4_b_o_aer(07,itype,cw_phase) = p_opcg4_b_o_cw07
            if(p_opcg5_b_o_cw07.ge.p1st) lptr_opcg5_b_o_aer(07,itype,cw_phase) = p_opcg5_b_o_cw07
            if(p_opcg6_b_o_cw07.ge.p1st) lptr_opcg6_b_o_aer(07,itype,cw_phase) = p_opcg6_b_o_cw07
            if(p_tanv_cw07.ge.p1st) lptr_tanv_aer(07,itype,cw_phase) = p_tanv_cw07
            if(p_sopnv_cw07.ge.p1st) lptr_sopnv_aer(07,itype,cw_phase) = p_sopnv_cw07
            if(p_ternv_cw07.ge.p1st) lptr_ternv_aer(07,itype,cw_phase) = p_ternv_cw07
            if(p_seqnv_cw07.ge.p1st) lptr_seqnv_aer(07,itype,cw_phase) = p_seqnv_cw07
            if(p_sianv_cw07.ge.p1st) lptr_sianv_aer(07,itype,cw_phase) = p_sianv_cw07



	    numptr_aer(07,itype,cw_phase)        = p_num_cw07
	  end if
	end if

	if (nsize_aer(itype) .ge. 8) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(08,itype,cw_phase)      = p_so4_cw08
	    lptr_no3_aer(08,itype,cw_phase)      = p_no3_cw08
	    lptr_cl_aer(08,itype,cw_phase)       = p_cl_cw08
	    lptr_msa_aer(08,itype,cw_phase)      = p_msa_cw08
	    lptr_co3_aer(08,itype,cw_phase)      = p_co3_cw08
	    lptr_nh4_aer(08,itype,cw_phase)      = p_nh4_cw08
	    lptr_na_aer(08,itype,cw_phase)       = p_na_cw08
	    lptr_ca_aer(08,itype,cw_phase)       = p_ca_cw08
	    lptr_oin_aer(08,itype,cw_phase)      = p_oin_cw08
	    lptr_oc_aer(08,itype,cw_phase)       = p_oc_cw08
	    lptr_bc_aer(08,itype,cw_phase)       = p_bc_cw08
        if(p_pcg1_b_c_cw08.ge.p1st) lptr_pcg1_b_c_aer(08,itype,cw_phase) = p_pcg1_b_c_cw08
            if(p_opcg1_b_c_cw08.ge.p1st) lptr_opcg1_b_c_aer(08,itype,cw_phase) = p_opcg1_b_c_cw08
            if(p_pcg1_b_o_cw08.ge.p1st) lptr_pcg1_b_o_aer(08,itype,cw_phase) = p_pcg1_b_o_cw08
            if(p_opcg1_b_o_cw08.ge.p1st) lptr_opcg1_b_o_aer(08,itype,cw_phase) = p_opcg1_b_o_cw08
            if(p_pcg1_f_c_cw08.ge.p1st) lptr_pcg1_f_c_aer(08,itype,cw_phase) = p_pcg1_f_c_cw08
            if(p_opcg1_f_c_cw08.ge.p1st) lptr_opcg1_f_c_aer(08,itype,cw_phase) = p_opcg1_f_c_cw08
            if(p_pcg1_f_o_cw08.ge.p1st) lptr_pcg1_f_o_aer(08,itype,cw_phase) = p_pcg1_f_o_cw08
            if(p_opcg1_f_o_cw08.ge.p1st) lptr_opcg1_f_o_aer(08,itype,cw_phase) = p_opcg1_f_o_cw08
            if(p_ant1_c_cw08.ge.p1st) lptr_ant1_c_aer(08,itype,cw_phase) = p_ant1_c_cw08
            if(p_ant2_c_cw08.ge.p1st) lptr_ant2_c_aer(08,itype,cw_phase) = p_ant2_c_cw08
            if(p_ant3_c_cw08.ge.p1st) lptr_ant3_c_aer(08,itype,cw_phase) = p_ant3_c_cw08
            if(p_ant4_c_cw08.ge.p1st) lptr_ant4_c_aer(08,itype,cw_phase) = p_ant4_c_cw08
            if(p_biog1_c_cw08.ge.p1st) lptr_biog1_c_aer(08,itype,cw_phase) = p_biog1_c_cw08
            if(p_biog2_c_cw08.ge.p1st) lptr_biog2_c_aer(08,itype,cw_phase) = p_biog2_c_cw08
            if(p_biog3_c_cw08.ge.p1st) lptr_biog3_c_aer(08,itype,cw_phase) = p_biog3_c_cw08
            if(p_biog4_c_cw08.ge.p1st) lptr_biog4_c_aer(08,itype,cw_phase) = p_biog4_c_cw08
            if(p_iepox_cw08.ge.p1st) lptr_iepox_aer(08,itype,cw_phase) = p_iepox_cw08
            if(p_iepoxos_cw08.ge.p1st) lptr_iepoxos_aer(08,itype,cw_phase) = p_iepoxos_cw08
            if(p_tetrol_cw08.ge.p1st) lptr_tetrol_aer(08,itype,cw_phase) = p_tetrol_cw08
            if(p_gly_cw08.ge.p1st) lptr_gly_aer(08,itype,cw_phase) = p_gly_cw08
            if(p_biog1_o_cw08.ge.p1st) lptr_biog1_o_aer(08,itype,cw_phase) = p_biog1_o_cw08
            if(p_biog2_o_cw08.ge.p1st) lptr_biog2_o_aer(08,itype,cw_phase) = p_biog2_o_cw08
            if(p_biog3_o_cw08.ge.p1st) lptr_biog3_o_aer(08,itype,cw_phase) = p_biog3_o_cw08
            if(p_biog4_o_cw08.ge.p1st) lptr_biog4_o_aer(08,itype,cw_phase) = p_biog4_o_cw08
            if(p_opcg2_b_o_cw08.ge.p1st) lptr_opcg2_b_o_aer(08,itype,cw_phase) = p_opcg2_b_o_cw08
            if(p_opcg3_b_o_cw08.ge.p1st) lptr_opcg3_b_o_aer(08,itype,cw_phase) = p_opcg3_b_o_cw08
            if(p_opcg4_b_o_cw08.ge.p1st) lptr_opcg4_b_o_aer(08,itype,cw_phase) = p_opcg4_b_o_cw08
            if(p_opcg5_b_o_cw08.ge.p1st) lptr_opcg5_b_o_aer(08,itype,cw_phase) = p_opcg5_b_o_cw08
            if(p_opcg6_b_o_cw08.ge.p1st) lptr_opcg6_b_o_aer(08,itype,cw_phase) = p_opcg6_b_o_cw08
            if(p_tanv_cw08.ge.p1st) lptr_tanv_aer(08,itype,cw_phase) = p_tanv_cw08
            if(p_sopnv_cw08.ge.p1st) lptr_sopnv_aer(08,itype,cw_phase) = p_sopnv_cw08
            if(p_ternv_cw08.ge.p1st) lptr_ternv_aer(08,itype,cw_phase) = p_ternv_cw08
            if(p_seqnv_cw08.ge.p1st) lptr_seqnv_aer(08,itype,cw_phase) = p_seqnv_cw08
            if(p_sianv_cw08.ge.p1st) lptr_sianv_aer(08,itype,cw_phase) = p_sianv_cw08



	    numptr_aer(08,itype,cw_phase)        = p_num_cw08
	  end if
	end if

	if (nsize_aer(itype) .ge. 9) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(09,itype,cw_phase)      = p_so4_cw09
	    lptr_no3_aer(09,itype,cw_phase)      = p_no3_cw09
	    lptr_cl_aer(09,itype,cw_phase)       = p_cl_cw09
	    lptr_msa_aer(09,itype,cw_phase)      = p_msa_cw09
	    lptr_co3_aer(09,itype,cw_phase)      = p_co3_cw09
	    lptr_nh4_aer(09,itype,cw_phase)      = p_nh4_cw09
	    lptr_na_aer(09,itype,cw_phase)       = p_na_cw09
	    lptr_ca_aer(09,itype,cw_phase)       = p_ca_cw09
	    lptr_oin_aer(09,itype,cw_phase)      = p_oin_cw09
	    lptr_oc_aer(09,itype,cw_phase)       = p_oc_cw09
	    lptr_bc_aer(09,itype,cw_phase)       = p_bc_cw09

            lptr_pcg1_b_c_aer(09,itype,cw_phase) = p_pcg1_b_c_cw09
            lptr_opcg1_b_c_aer(09,itype,cw_phase) = p_opcg1_b_c_cw09
            lptr_pcg1_b_o_aer(09,itype,cw_phase) = p_pcg1_b_o_cw09
            lptr_opcg1_b_o_aer(09,itype,cw_phase) = p_opcg1_b_o_cw09
            lptr_pcg1_f_c_aer(09,itype,cw_phase) = p_pcg1_f_c_cw09
            lptr_opcg1_f_c_aer(09,itype,cw_phase) = p_opcg1_f_c_cw09
            lptr_pcg1_f_o_aer(09,itype,cw_phase) = p_pcg1_f_o_cw09
            lptr_opcg1_f_o_aer(09,itype,cw_phase) = p_opcg1_f_o_cw09
            lptr_ant1_c_aer(09,itype,cw_phase) = p_ant1_c_cw09
            lptr_biog1_c_aer(09,itype,cw_phase) = p_biog1_c_cw09

	    numptr_aer(09,itype,cw_phase)        = p_num_cw09
	  end if
	end if

	if (nsize_aer(itype) .ge. 10) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(10,itype,cw_phase)      = p_so4_cw10
	    lptr_no3_aer(10,itype,cw_phase)      = p_no3_cw10
	    lptr_cl_aer(10,itype,cw_phase)       = p_cl_cw10
	    lptr_msa_aer(10,itype,cw_phase)      = p_msa_cw10
	    lptr_co3_aer(10,itype,cw_phase)      = p_co3_cw10
	    lptr_nh4_aer(10,itype,cw_phase)      = p_nh4_cw10
	    lptr_na_aer(10,itype,cw_phase)       = p_na_cw10
	    lptr_ca_aer(10,itype,cw_phase)       = p_ca_cw10
	    lptr_oin_aer(10,itype,cw_phase)      = p_oin_cw10
	    lptr_oc_aer(10,itype,cw_phase)       = p_oc_cw10
	    lptr_bc_aer(10,itype,cw_phase)       = p_bc_cw10

            lptr_pcg1_b_c_aer(10,itype,cw_phase) = p_pcg1_b_c_cw10
            lptr_opcg1_b_c_aer(10,itype,cw_phase) = p_opcg1_b_c_cw10
            lptr_pcg1_b_o_aer(10,itype,cw_phase) = p_pcg1_b_o_cw10
            lptr_opcg1_b_o_aer(10,itype,cw_phase) = p_opcg1_b_o_cw10
            lptr_pcg1_f_c_aer(10,itype,cw_phase) = p_pcg1_f_c_cw10
            lptr_opcg1_f_c_aer(10,itype,cw_phase) = p_opcg1_f_c_cw10
            lptr_pcg1_f_o_aer(10,itype,cw_phase) = p_pcg1_f_o_cw10
            lptr_opcg1_f_o_aer(10,itype,cw_phase) = p_opcg1_f_o_cw10
            lptr_ant1_c_aer(10,itype,cw_phase) = p_ant1_c_cw10
            lptr_biog1_c_aer(10,itype,cw_phase) = p_biog1_c_cw10

	    numptr_aer(10,itype,cw_phase)        = p_num_cw10
	  end if
	end if

	if (nsize_aer(itype) .ge. 11) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(11,itype,cw_phase)      = p_so4_cw11
	    lptr_no3_aer(11,itype,cw_phase)      = p_no3_cw11
	    lptr_cl_aer(11,itype,cw_phase)       = p_cl_cw11
	    lptr_msa_aer(11,itype,cw_phase)      = p_msa_cw11
	    lptr_co3_aer(11,itype,cw_phase)      = p_co3_cw11
	    lptr_nh4_aer(11,itype,cw_phase)      = p_nh4_cw11
	    lptr_na_aer(11,itype,cw_phase)       = p_na_cw11
	    lptr_ca_aer(11,itype,cw_phase)       = p_ca_cw11
	    lptr_oin_aer(11,itype,cw_phase)      = p_oin_cw11
	    lptr_oc_aer(11,itype,cw_phase)       = p_oc_cw11
	    lptr_bc_aer(11,itype,cw_phase)       = p_bc_cw11

            lptr_pcg1_b_c_aer(11,itype,cw_phase) = p_pcg1_b_c_cw11
            lptr_opcg1_b_c_aer(11,itype,cw_phase) = p_opcg1_b_c_cw11
            lptr_pcg1_b_o_aer(11,itype,cw_phase) = p_pcg1_b_o_cw11
            lptr_opcg1_b_o_aer(11,itype,cw_phase) = p_opcg1_b_o_cw11
            lptr_pcg1_f_c_aer(11,itype,cw_phase) = p_pcg1_f_c_cw11
            lptr_opcg1_f_c_aer(11,itype,cw_phase) = p_opcg1_f_c_cw11
            lptr_pcg1_f_o_aer(11,itype,cw_phase) = p_pcg1_f_o_cw11
            lptr_opcg1_f_o_aer(11,itype,cw_phase) = p_opcg1_f_o_cw11
            lptr_ant1_c_aer(11,itype,cw_phase) = p_ant1_c_cw11
            lptr_biog1_c_aer(11,itype,cw_phase) = p_biog1_c_cw11

	    numptr_aer(11,itype,cw_phase)        = p_num_cw11
	  end if
	end if

	if (nsize_aer(itype) .ge. 12) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(12,itype,cw_phase)      = p_so4_cw12
	    lptr_no3_aer(12,itype,cw_phase)      = p_no3_cw12
	    lptr_cl_aer(12,itype,cw_phase)       = p_cl_cw12
	    lptr_msa_aer(12,itype,cw_phase)      = p_msa_cw12
	    lptr_co3_aer(12,itype,cw_phase)      = p_co3_cw12
	    lptr_nh4_aer(12,itype,cw_phase)      = p_nh4_cw12
	    lptr_na_aer(12,itype,cw_phase)       = p_na_cw12
	    lptr_ca_aer(12,itype,cw_phase)       = p_ca_cw12
	    lptr_oin_aer(12,itype,cw_phase)      = p_oin_cw12
	    lptr_oc_aer(12,itype,cw_phase)       = p_oc_cw12
	    lptr_bc_aer(12,itype,cw_phase)       = p_bc_cw12

            lptr_pcg1_b_c_aer(12,itype,cw_phase) = p_pcg1_b_c_cw12
            lptr_opcg1_b_c_aer(12,itype,cw_phase) = p_opcg1_b_c_cw12
            lptr_pcg1_b_o_aer(12,itype,cw_phase) = p_pcg1_b_o_cw12
            lptr_opcg1_b_o_aer(12,itype,cw_phase) = p_opcg1_b_o_cw12
            lptr_pcg1_f_c_aer(12,itype,cw_phase) = p_pcg1_f_c_cw12
            lptr_opcg1_f_c_aer(12,itype,cw_phase) = p_opcg1_f_c_cw12
            lptr_pcg1_f_o_aer(12,itype,cw_phase) = p_pcg1_f_o_cw12
            lptr_opcg1_f_o_aer(12,itype,cw_phase) = p_opcg1_f_o_cw12
            lptr_ant1_c_aer(12,itype,cw_phase) = p_ant1_c_cw12
            lptr_biog1_c_aer(12,itype,cw_phase) = p_biog1_c_cw12

	    numptr_aer(12,itype,cw_phase)        = p_num_cw12
	  end if
	end if

	if (nsize_aer(itype) .ge. 13) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(13,itype,cw_phase)      = p_so4_cw13
	    lptr_no3_aer(13,itype,cw_phase)      = p_no3_cw13
	    lptr_cl_aer(13,itype,cw_phase)       = p_cl_cw13
	    lptr_msa_aer(13,itype,cw_phase)      = p_msa_cw13
	    lptr_co3_aer(13,itype,cw_phase)      = p_co3_cw13
	    lptr_nh4_aer(13,itype,cw_phase)      = p_nh4_cw13
	    lptr_na_aer(13,itype,cw_phase)       = p_na_cw13
	    lptr_ca_aer(13,itype,cw_phase)       = p_ca_cw13
	    lptr_oin_aer(13,itype,cw_phase)      = p_oin_cw13
	    lptr_oc_aer(13,itype,cw_phase)       = p_oc_cw13
	    lptr_bc_aer(13,itype,cw_phase)       = p_bc_cw13

            lptr_pcg1_b_c_aer(13,itype,cw_phase) = p_pcg1_b_c_cw13
            lptr_opcg1_b_c_aer(13,itype,cw_phase) = p_opcg1_b_c_cw13
            lptr_pcg1_b_o_aer(13,itype,cw_phase) = p_pcg1_b_o_cw13
            lptr_opcg1_b_o_aer(13,itype,cw_phase) = p_opcg1_b_o_cw13
            lptr_pcg1_f_c_aer(13,itype,cw_phase) = p_pcg1_f_c_cw13
            lptr_opcg1_f_c_aer(13,itype,cw_phase) = p_opcg1_f_c_cw13
            lptr_pcg1_f_o_aer(13,itype,cw_phase) = p_pcg1_f_o_cw13
            lptr_opcg1_f_o_aer(13,itype,cw_phase) = p_opcg1_f_o_cw13
            lptr_ant1_c_aer(13,itype,cw_phase) = p_ant1_c_cw13
            lptr_biog1_c_aer(13,itype,cw_phase) = p_biog1_c_cw13

	    numptr_aer(13,itype,cw_phase)        = p_num_cw13
	  end if
	end if

	if (nsize_aer(itype) .ge. 14) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(14,itype,cw_phase)      = p_so4_cw14
	    lptr_no3_aer(14,itype,cw_phase)      = p_no3_cw14
	    lptr_cl_aer(14,itype,cw_phase)       = p_cl_cw14
	    lptr_msa_aer(14,itype,cw_phase)      = p_msa_cw14
	    lptr_co3_aer(14,itype,cw_phase)      = p_co3_cw14
	    lptr_nh4_aer(14,itype,cw_phase)      = p_nh4_cw14
	    lptr_na_aer(14,itype,cw_phase)       = p_na_cw14
	    lptr_ca_aer(14,itype,cw_phase)       = p_ca_cw14
	    lptr_oin_aer(14,itype,cw_phase)      = p_oin_cw14
	    lptr_oc_aer(14,itype,cw_phase)       = p_oc_cw14
	    lptr_bc_aer(14,itype,cw_phase)       = p_bc_cw14

            lptr_pcg1_b_c_aer(14,itype,cw_phase) = p_pcg1_b_c_cw14
            lptr_opcg1_b_c_aer(14,itype,cw_phase) = p_opcg1_b_c_cw14
            lptr_pcg1_b_o_aer(14,itype,cw_phase) = p_pcg1_b_o_cw14
            lptr_opcg1_b_o_aer(14,itype,cw_phase) = p_opcg1_b_o_cw14
            lptr_pcg1_f_c_aer(14,itype,cw_phase) = p_pcg1_f_c_cw14
            lptr_opcg1_f_c_aer(14,itype,cw_phase) = p_opcg1_f_c_cw14
            lptr_pcg1_f_o_aer(14,itype,cw_phase) = p_pcg1_f_o_cw14
            lptr_opcg1_f_o_aer(14,itype,cw_phase) = p_opcg1_f_o_cw14
            lptr_ant1_c_aer(14,itype,cw_phase) = p_ant1_c_cw14
            lptr_biog1_c_aer(14,itype,cw_phase) = p_biog1_c_cw14

	    numptr_aer(14,itype,cw_phase)        = p_num_cw14
	  end if
	end if

	if (nsize_aer(itype) .ge. 15) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(15,itype,cw_phase)      = p_so4_cw15
	    lptr_no3_aer(15,itype,cw_phase)      = p_no3_cw15
	    lptr_cl_aer(15,itype,cw_phase)       = p_cl_cw15
	    lptr_msa_aer(15,itype,cw_phase)      = p_msa_cw15
	    lptr_co3_aer(15,itype,cw_phase)      = p_co3_cw15
	    lptr_nh4_aer(15,itype,cw_phase)      = p_nh4_cw15
	    lptr_na_aer(15,itype,cw_phase)       = p_na_cw15
	    lptr_ca_aer(15,itype,cw_phase)       = p_ca_cw15
	    lptr_oin_aer(15,itype,cw_phase)      = p_oin_cw15
	    lptr_oc_aer(15,itype,cw_phase)       = p_oc_cw15
	    lptr_bc_aer(15,itype,cw_phase)       = p_bc_cw15

            lptr_pcg1_b_c_aer(15,itype,cw_phase) = p_pcg1_b_c_cw15
            lptr_opcg1_b_c_aer(15,itype,cw_phase) = p_opcg1_b_c_cw15
            lptr_pcg1_b_o_aer(15,itype,cw_phase) = p_pcg1_b_o_cw15
            lptr_opcg1_b_o_aer(15,itype,cw_phase) = p_opcg1_b_o_cw15
            lptr_pcg1_f_c_aer(15,itype,cw_phase) = p_pcg1_f_c_cw15
            lptr_opcg1_f_c_aer(15,itype,cw_phase) = p_opcg1_f_c_cw15
            lptr_pcg1_f_o_aer(15,itype,cw_phase) = p_pcg1_f_o_cw15
            lptr_opcg1_f_o_aer(15,itype,cw_phase) = p_opcg1_f_o_cw15
            lptr_ant1_c_aer(15,itype,cw_phase) = p_ant1_c_cw15
            lptr_biog1_c_aer(15,itype,cw_phase) = p_biog1_c_cw15

	    numptr_aer(15,itype,cw_phase)        = p_num_cw15
	  end if
	end if

	if (nsize_aer(itype) .ge. 16) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(16,itype,cw_phase)      = p_so4_cw16
	    lptr_no3_aer(16,itype,cw_phase)      = p_no3_cw16
	    lptr_cl_aer(16,itype,cw_phase)       = p_cl_cw16
	    lptr_msa_aer(16,itype,cw_phase)      = p_msa_cw16
	    lptr_co3_aer(16,itype,cw_phase)      = p_co3_cw16
	    lptr_nh4_aer(16,itype,cw_phase)      = p_nh4_cw16
	    lptr_na_aer(16,itype,cw_phase)       = p_na_cw16
	    lptr_ca_aer(16,itype,cw_phase)       = p_ca_cw16
	    lptr_oin_aer(16,itype,cw_phase)      = p_oin_cw16
	    lptr_oc_aer(16,itype,cw_phase)       = p_oc_cw16
	    lptr_bc_aer(16,itype,cw_phase)       = p_bc_cw16

            lptr_pcg1_b_c_aer(16,itype,cw_phase) = p_pcg1_b_c_cw16
            lptr_opcg1_b_c_aer(16,itype,cw_phase) = p_opcg1_b_c_cw16
            lptr_pcg1_b_o_aer(16,itype,cw_phase) = p_pcg1_b_o_cw16
            lptr_opcg1_b_o_aer(16,itype,cw_phase) = p_opcg1_b_o_cw16
            lptr_pcg1_f_c_aer(16,itype,cw_phase) = p_pcg1_f_c_cw16
            lptr_opcg1_f_c_aer(16,itype,cw_phase) = p_opcg1_f_c_cw16
            lptr_pcg1_f_o_aer(16,itype,cw_phase) = p_pcg1_f_o_cw16
            lptr_opcg1_f_o_aer(16,itype,cw_phase) = p_opcg1_f_o_cw16
            lptr_ant1_c_aer(16,itype,cw_phase) = p_ant1_c_cw16
            lptr_biog1_c_aer(16,itype,cw_phase) = p_biog1_c_cw16

	    numptr_aer(16,itype,cw_phase)        = p_num_cw16
	  end if
	end if

	if (nsize_aer(itype) .ge. 17) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(17,itype,cw_phase)      = p_so4_cw17
	    lptr_no3_aer(17,itype,cw_phase)      = p_no3_cw17
	    lptr_cl_aer(17,itype,cw_phase)       = p_cl_cw17
	    lptr_msa_aer(17,itype,cw_phase)      = p_msa_cw17
	    lptr_co3_aer(17,itype,cw_phase)      = p_co3_cw17
	    lptr_nh4_aer(17,itype,cw_phase)      = p_nh4_cw17
	    lptr_na_aer(17,itype,cw_phase)       = p_na_cw17
	    lptr_ca_aer(17,itype,cw_phase)       = p_ca_cw17
	    lptr_oin_aer(17,itype,cw_phase)      = p_oin_cw17
	    lptr_oc_aer(17,itype,cw_phase)       = p_oc_cw17
	    lptr_bc_aer(17,itype,cw_phase)       = p_bc_cw17

            lptr_pcg1_b_c_aer(17,itype,cw_phase) = p_pcg1_b_c_cw17
            lptr_opcg1_b_c_aer(17,itype,cw_phase) = p_opcg1_b_c_cw17
            lptr_pcg1_b_o_aer(17,itype,cw_phase) = p_pcg1_b_o_cw17
            lptr_opcg1_b_o_aer(17,itype,cw_phase) = p_opcg1_b_o_cw17
            lptr_pcg1_f_c_aer(17,itype,cw_phase) = p_pcg1_f_c_cw17
            lptr_opcg1_f_c_aer(17,itype,cw_phase) = p_opcg1_f_c_cw17
            lptr_pcg1_f_o_aer(17,itype,cw_phase) = p_pcg1_f_o_cw17
            lptr_opcg1_f_o_aer(17,itype,cw_phase) = p_opcg1_f_o_cw17
            lptr_ant1_c_aer(17,itype,cw_phase) = p_ant1_c_cw17
            lptr_biog1_c_aer(17,itype,cw_phase) = p_biog1_c_cw17

	    numptr_aer(17,itype,cw_phase)        = p_num_cw17
	  end if
	end if

	if (nsize_aer(itype) .ge. 18) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(18,itype,cw_phase)      = p_so4_cw18
	    lptr_no3_aer(18,itype,cw_phase)      = p_no3_cw18
	    lptr_cl_aer(18,itype,cw_phase)       = p_cl_cw18
	    lptr_msa_aer(18,itype,cw_phase)      = p_msa_cw18
	    lptr_co3_aer(18,itype,cw_phase)      = p_co3_cw18
	    lptr_nh4_aer(18,itype,cw_phase)      = p_nh4_cw18
	    lptr_na_aer(18,itype,cw_phase)       = p_na_cw18
	    lptr_ca_aer(18,itype,cw_phase)       = p_ca_cw18
	    lptr_oin_aer(18,itype,cw_phase)      = p_oin_cw18
	    lptr_oc_aer(18,itype,cw_phase)       = p_oc_cw18
	    lptr_bc_aer(18,itype,cw_phase)       = p_bc_cw18

            lptr_pcg1_b_c_aer(18,itype,cw_phase) = p_pcg1_b_c_cw18
            lptr_opcg1_b_c_aer(18,itype,cw_phase) = p_opcg1_b_c_cw18
            lptr_pcg1_b_o_aer(18,itype,cw_phase) = p_pcg1_b_o_cw18
            lptr_opcg1_b_o_aer(18,itype,cw_phase) = p_opcg1_b_o_cw18
            lptr_pcg1_f_c_aer(18,itype,cw_phase) = p_pcg1_f_c_cw18
            lptr_opcg1_f_c_aer(18,itype,cw_phase) = p_opcg1_f_c_cw18
            lptr_pcg1_f_o_aer(18,itype,cw_phase) = p_pcg1_f_o_cw18
            lptr_opcg1_f_o_aer(18,itype,cw_phase) = p_opcg1_f_o_cw18
            lptr_ant1_c_aer(18,itype,cw_phase) = p_ant1_c_cw18
            lptr_biog1_c_aer(18,itype,cw_phase) = p_biog1_c_cw18

	    numptr_aer(18,itype,cw_phase)        = p_num_cw18
	  end if
	end if

	if (nsize_aer(itype) .ge. 19) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(19,itype,cw_phase)      = p_so4_cw19
	    lptr_no3_aer(19,itype,cw_phase)      = p_no3_cw19
	    lptr_cl_aer(19,itype,cw_phase)       = p_cl_cw19
	    lptr_msa_aer(19,itype,cw_phase)      = p_msa_cw19
	    lptr_co3_aer(19,itype,cw_phase)      = p_co3_cw19
	    lptr_nh4_aer(19,itype,cw_phase)      = p_nh4_cw19
	    lptr_na_aer(19,itype,cw_phase)       = p_na_cw19
	    lptr_ca_aer(19,itype,cw_phase)       = p_ca_cw19
	    lptr_oin_aer(19,itype,cw_phase)      = p_oin_cw19
	    lptr_oc_aer(19,itype,cw_phase)       = p_oc_cw19
	    lptr_bc_aer(19,itype,cw_phase)       = p_bc_cw19

            lptr_pcg1_b_c_aer(19,itype,cw_phase) = p_pcg1_b_c_cw19
            lptr_opcg1_b_c_aer(19,itype,cw_phase) = p_opcg1_b_c_cw19
            lptr_pcg1_b_o_aer(19,itype,cw_phase) = p_pcg1_b_o_cw19
            lptr_opcg1_b_o_aer(19,itype,cw_phase) = p_opcg1_b_o_cw19
            lptr_pcg1_f_c_aer(19,itype,cw_phase) = p_pcg1_f_c_cw19
            lptr_opcg1_f_c_aer(19,itype,cw_phase) = p_opcg1_f_c_cw19
            lptr_pcg1_f_o_aer(19,itype,cw_phase) = p_pcg1_f_o_cw19
            lptr_opcg1_f_o_aer(19,itype,cw_phase) = p_opcg1_f_o_cw19
            lptr_ant1_c_aer(19,itype,cw_phase) = p_ant1_c_cw19
            lptr_biog1_c_aer(19,itype,cw_phase) = p_biog1_c_cw19

	    numptr_aer(19,itype,cw_phase)        = p_num_cw19
	  end if
	end if

	if (nsize_aer(itype) .ge. 20) then
	  if (cw_phase .gt. 0) then
	    lptr_so4_aer(20,itype,cw_phase)      = p_so4_cw20
	    lptr_no3_aer(20,itype,cw_phase)      = p_no3_cw20
	    lptr_cl_aer(20,itype,cw_phase)       = p_cl_cw20
	    lptr_msa_aer(20,itype,cw_phase)      = p_msa_cw20
	    lptr_co3_aer(20,itype,cw_phase)      = p_co3_cw20
	    lptr_nh4_aer(20,itype,cw_phase)      = p_nh4_cw20
	    lptr_na_aer(20,itype,cw_phase)       = p_na_cw20
	    lptr_ca_aer(20,itype,cw_phase)       = p_ca_cw20
	    lptr_oin_aer(20,itype,cw_phase)      = p_oin_cw20
	    lptr_oc_aer(20,itype,cw_phase)       = p_oc_cw20
	    lptr_bc_aer(20,itype,cw_phase)       = p_bc_cw20

            lptr_pcg1_b_c_aer(20,itype,cw_phase) = p_pcg1_b_c_cw20
            lptr_opcg1_b_c_aer(20,itype,cw_phase) = p_opcg1_b_c_cw20
            lptr_pcg1_b_o_aer(20,itype,cw_phase) = p_pcg1_b_o_cw20
            lptr_opcg1_b_o_aer(20,itype,cw_phase) = p_opcg1_b_o_cw20
            lptr_pcg1_f_c_aer(20,itype,cw_phase) = p_pcg1_f_c_cw20
            lptr_opcg1_f_c_aer(20,itype,cw_phase) = p_opcg1_f_c_cw20
            lptr_pcg1_f_o_aer(20,itype,cw_phase) = p_pcg1_f_o_cw20
            lptr_opcg1_f_o_aer(20,itype,cw_phase) = p_opcg1_f_o_cw20
            lptr_ant1_c_aer(20,itype,cw_phase) = p_ant1_c_cw20
            lptr_biog1_c_aer(20,itype,cw_phase) = p_biog1_c_cw20

	    numptr_aer(20,itype,cw_phase)        = p_num_cw20
	  end if
	end if








	do l = 1, l2maxd
	    write( name(l), '(a,i4.4,15x)' ) 'r', l
	end do
	massptr_aer(:,:,:,:) = -999888777
	mastercompptr_aer(:,:) = -999888777

	do 2800 itype = 1, ntype_aer

	if (itype .eq. 1) then
	    typetxt = ' '
	    ntypetxt = 1
	    if (ntype_aer .gt. 1) then
		typetxt = '_t1'
		ntypetxt = 3
	    end if
	else if (itype .le. 9) then
	    write(typetxt,'(a,i1)') '_t', itype
	    ntypetxt = 3
	else if (itype .le. 99) then
	    write(typetxt,'(a,i2)') '_t', itype
	    ntypetxt = 4
	else
	    typetxt = '_t?'
	    ntypetxt = 3
	end if

	ncomp_dum(:,:) = 0
	ncomp_plustracer_dum(:,:) = 0

	do 2700 isize = 1, nsize_aer(itype)
	n =isize

	if (isize .le. 9) then
	    write(sizetxt,'(i1)') isize
	    nsizetxt = 1
	else if (isize .le. 99) then
	    write(sizetxt,'(i2)') isize
	    nsizetxt = 2
	else if (isize .le. 999) then
	    write(sizetxt,'(i3)') isize
	    nsizetxt = 3
	else
	    sizetxt = 's?'
	    nsizetxt = 2
	end if


	do 2600 iphase = 1, nphase_aer

	if (iphase .eq. ai_phase) then
	    phasetxt = 'a'
	    nphasetxt = 1
	else if (iphase .eq. cw_phase) then
	    phasetxt = 'cw'
	    nphasetxt = 2
	else 
	    phasetxt = 'p?'
	    nphasetxt = 2
	end if


	do 2500 l_mastercomp = -2, ntot_mastercomp_aer

	iaddto_ncomp = 1
	iaddto_ncomp_plustracer = 1

	if (l_mastercomp .eq. -2) then
	    iaddto_ncomp = 0
	    iaddto_ncomp_plustracer = 0
	    lptr_dum = numptr_aer(n,itype,iphase)
	    mcindx_dum = -2
	    spectxt = 'numb_'
	    nspectxt = 5

	else if (l_mastercomp .eq. -1) then
	    if (iphase .ne. ai_phase) goto 2500
	    iaddto_ncomp = 0
	    iaddto_ncomp_plustracer = 0
	    lptr_dum = waterptr_aer(n,itype)
	    mcindx_dum = -1
	    spectxt = 'water_'
	    nspectxt = 6

	else if (l_mastercomp .eq. 0) then
	    if (iphase .ne. ai_phase) goto 2500
	    iaddto_ncomp = 0
	    iaddto_ncomp_plustracer = 0
	    lptr_dum = hyswptr_aer(n,itype)
	    mcindx_dum = 0
	    spectxt = 'hysw_'
	    nspectxt = 5

	else if (l_mastercomp .eq. mastercompindx_so4_aer) then
	    lptr_dum = lptr_so4_aer(n,itype,iphase)
	    mcindx_dum = mastercompindx_so4_aer
	    spectxt = 'so4_'
	    nspectxt = 4

	else if (l_mastercomp .eq. mastercompindx_no3_aer) then
	    lptr_dum = lptr_no3_aer(n,itype,iphase)
	    mcindx_dum = mastercompindx_no3_aer
	    spectxt = 'no3_'
	    nspectxt = 4

	else if (l_mastercomp .eq. mastercompindx_cl_aer) then
	    lptr_dum = lptr_cl_aer(n,itype,iphase)
	    mcindx_dum = mastercompindx_cl_aer
	    spectxt = 'cl_'
	    nspectxt = 3

	else if (l_mastercomp .eq. mastercompindx_co3_aer) then
	    lptr_dum = lptr_co3_aer(n,itype,iphase)
	    mcindx_dum = mastercompindx_co3_aer
	    spectxt = 'co3_'
	    nspectxt = 4

	else if (l_mastercomp .eq. mastercompindx_nh4_aer) then
	    lptr_dum = lptr_nh4_aer(n,itype,iphase)
	    mcindx_dum = mastercompindx_nh4_aer
	    spectxt = 'nh4_'
	    nspectxt = 4

	else if (l_mastercomp .eq. mastercompindx_na_aer) then
	    lptr_dum = lptr_na_aer(n,itype,iphase)
	    mcindx_dum = mastercompindx_na_aer
	    spectxt = 'na_'
	    nspectxt = 3

	else if (l_mastercomp .eq. mastercompindx_ca_aer) then
	    lptr_dum = lptr_ca_aer(n,itype,iphase)
	    mcindx_dum = mastercompindx_ca_aer
	    spectxt = 'ca_'
	    nspectxt = 3

	else if (l_mastercomp .eq. mastercompindx_oin_aer) then
	    lptr_dum = lptr_oin_aer(n,itype,iphase)
	    mcindx_dum = mastercompindx_oin_aer
	    spectxt = 'oin_'
	    nspectxt = 4

	else if (l_mastercomp .eq. mastercompindx_oc_aer) then
	    lptr_dum = lptr_oc_aer(n,itype,iphase)
	    mcindx_dum = mastercompindx_oc_aer
	    spectxt = 'oc_'
	    nspectxt = 3

	else if (l_mastercomp .eq. mastercompindx_bc_aer) then
	    lptr_dum = lptr_bc_aer(n,itype,iphase)
	    mcindx_dum = mastercompindx_bc_aer
	    spectxt = 'bc_'
	    nspectxt = 3


        else if (l_mastercomp .eq. mastercompindx_pcg1_b_c_aer) then
            lptr_dum = lptr_pcg1_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg1_b_c_aer
            spectxt = 'pcg1_b_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg2_b_c_aer) then
            lptr_dum = lptr_pcg2_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg2_b_c_aer
            spectxt = 'pcg2_b_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg3_b_c_aer) then
            lptr_dum = lptr_pcg3_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg3_b_c_aer
            spectxt = 'pcg3_b_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg4_b_c_aer) then
            lptr_dum = lptr_pcg4_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg4_b_c_aer
            spectxt = 'pcg4_b_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg5_b_c_aer) then
            lptr_dum = lptr_pcg5_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg5_b_c_aer
            spectxt = 'pcg5_b_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg6_b_c_aer) then
            lptr_dum = lptr_pcg6_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg6_b_c_aer
            spectxt = 'pcg6_b_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg7_b_c_aer) then
            lptr_dum = lptr_pcg7_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg7_b_c_aer
            spectxt = 'pcg7_b_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg8_b_c_aer) then
            lptr_dum = lptr_pcg8_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg8_b_c_aer
            spectxt = 'pcg8_b_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg9_b_c_aer) then
            lptr_dum = lptr_pcg9_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg9_b_c_aer
            spectxt = 'pcg9_b_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg1_b_o_aer) then
            lptr_dum = lptr_pcg1_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg1_b_o_aer
            spectxt = 'pcg1_b_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg2_b_o_aer) then
            lptr_dum = lptr_pcg2_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg2_b_o_aer
            spectxt = 'pcg2_b_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg3_b_o_aer) then
            lptr_dum = lptr_pcg3_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg3_b_o_aer
            spectxt = 'pcg3_b_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg4_b_o_aer) then
            lptr_dum = lptr_pcg4_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg4_b_o_aer
            spectxt = 'pcg4_b_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg5_b_o_aer) then
            lptr_dum = lptr_pcg5_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg5_b_o_aer
            spectxt = 'pcg5_b_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg6_b_o_aer) then
            lptr_dum = lptr_pcg6_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg6_b_o_aer
            spectxt = 'pcg6_b_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg7_b_o_aer) then
            lptr_dum = lptr_pcg7_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg7_b_o_aer
            spectxt = 'pcg7_b_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg8_b_o_aer) then
            lptr_dum = lptr_pcg8_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg8_b_o_aer
            spectxt = 'pcg8_b_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg9_b_o_aer) then
            lptr_dum = lptr_pcg9_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg9_b_o_aer
            spectxt = 'pcg9_b_o_'
            nspectxt = 9
        else if (l_mastercomp .eq. mastercompindx_opcg1_b_c_aer) then
            lptr_dum = lptr_opcg1_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg1_b_c_aer
            spectxt = 'opcg1_b_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg2_b_c_aer) then
            lptr_dum = lptr_opcg2_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg2_b_c_aer
            spectxt = 'opcg2_b_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg3_b_c_aer) then
            lptr_dum = lptr_opcg3_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg3_b_c_aer
            spectxt = 'opcg3_b_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg4_b_c_aer) then
            lptr_dum = lptr_opcg4_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg4_b_c_aer
            spectxt = 'opcg4_b_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg5_b_c_aer) then
            lptr_dum = lptr_opcg5_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg5_b_c_aer
            spectxt = 'opcg5_b_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg6_b_c_aer) then
            lptr_dum = lptr_opcg6_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg6_b_c_aer
            spectxt = 'opcg6_b_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg7_b_c_aer) then
            lptr_dum = lptr_opcg7_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg7_b_c_aer
            spectxt = 'opcg7_b_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg8_b_c_aer) then
            lptr_dum = lptr_opcg8_b_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg8_b_c_aer
            spectxt = 'opcg8_b_c_'
            nspectxt = 10
        else if (l_mastercomp .eq. mastercompindx_opcg1_b_o_aer) then
            lptr_dum = lptr_opcg1_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg1_b_o_aer
            spectxt = 'opcg1_b_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg2_b_o_aer) then
            lptr_dum = lptr_opcg2_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg2_b_o_aer
            spectxt = 'opcg2_b_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg3_b_o_aer) then
            lptr_dum = lptr_opcg3_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg3_b_o_aer
            spectxt = 'opcg3_b_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg4_b_o_aer) then
            lptr_dum = lptr_opcg4_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg4_b_o_aer
            spectxt = 'opcg4_b_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg5_b_o_aer) then
            lptr_dum = lptr_opcg5_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg5_b_o_aer
            spectxt = 'opcg5_b_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg6_b_o_aer) then
            lptr_dum = lptr_opcg6_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg6_b_o_aer
            spectxt = 'opcg6_b_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg7_b_o_aer) then
            lptr_dum = lptr_opcg7_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg7_b_o_aer
            spectxt = 'opcg7_b_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg8_b_o_aer) then
            lptr_dum = lptr_opcg8_b_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg8_b_o_aer
            spectxt = 'opcg8_b_o_'
            nspectxt = 10
        else if (l_mastercomp .eq. mastercompindx_pcg1_f_c_aer) then
            lptr_dum = lptr_pcg1_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg1_f_c_aer
            spectxt = 'pcg1_f_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg2_f_c_aer) then
            lptr_dum = lptr_pcg2_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg2_f_c_aer
            spectxt = 'pcg2_f_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg3_f_c_aer) then
            lptr_dum = lptr_pcg3_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg3_f_c_aer
            spectxt = 'pcg3_f_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg4_f_c_aer) then
            lptr_dum = lptr_pcg4_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg4_f_c_aer
            spectxt = 'pcg4_f_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg5_f_c_aer) then
            lptr_dum = lptr_pcg5_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg5_f_c_aer
            spectxt = 'pcg5_f_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg6_f_c_aer) then
            lptr_dum = lptr_pcg6_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg6_f_c_aer
            spectxt = 'pcg6_f_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg7_f_c_aer) then
            lptr_dum = lptr_pcg7_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg7_f_c_aer
            spectxt = 'pcg7_f_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg8_f_c_aer) then
            lptr_dum = lptr_pcg8_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg8_f_c_aer
            spectxt = 'pcg8_f_c_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg9_f_c_aer) then
            lptr_dum = lptr_pcg9_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg9_f_c_aer
            spectxt = 'pcg9_f_c_'
            nspectxt = 9
        else if (l_mastercomp .eq. mastercompindx_pcg1_f_o_aer) then
            lptr_dum = lptr_pcg1_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg1_f_o_aer
            spectxt = 'pcg1_f_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg2_f_o_aer) then
            lptr_dum = lptr_pcg2_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg2_f_o_aer
            spectxt = 'pcg2_f_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg3_f_o_aer) then
            lptr_dum = lptr_pcg3_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg3_f_o_aer
            spectxt = 'pcg3_f_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg4_f_o_aer) then
            lptr_dum = lptr_pcg4_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg4_f_o_aer
            spectxt = 'pcg4_f_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg5_f_o_aer) then
            lptr_dum = lptr_pcg5_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg5_f_o_aer
            spectxt = 'pcg5_f_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg6_f_o_aer) then
            lptr_dum = lptr_pcg6_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg6_f_o_aer
            spectxt = 'pcg6_f_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg7_f_o_aer) then
            lptr_dum = lptr_pcg7_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg7_f_o_aer
            spectxt = 'pcg7_f_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg8_f_o_aer) then
            lptr_dum = lptr_pcg8_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg8_f_o_aer
            spectxt = 'pcg8_f_o_'
            nspectxt = 9

        else if (l_mastercomp .eq. mastercompindx_pcg9_f_o_aer) then
            lptr_dum = lptr_pcg9_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_pcg9_f_o_aer
            spectxt = 'pcg9_f_o_'
            nspectxt = 9
        else if (l_mastercomp .eq. mastercompindx_opcg1_f_c_aer) then
            lptr_dum = lptr_opcg1_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg1_f_c_aer
            spectxt = 'opcg1_f_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg2_f_c_aer) then
            lptr_dum = lptr_opcg2_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg2_f_c_aer
            spectxt = 'opcg2_f_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg3_f_c_aer) then
            lptr_dum = lptr_opcg3_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg3_f_c_aer
            spectxt = 'opcg3_f_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg4_f_c_aer) then
            lptr_dum = lptr_opcg4_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg4_f_c_aer
            spectxt = 'opcg4_f_c_'
            nspectxt = 10
        else if (l_mastercomp .eq. mastercompindx_opcg5_f_c_aer) then
            lptr_dum = lptr_opcg5_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg5_f_c_aer
            spectxt = 'opcg5_f_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg6_f_c_aer) then
            lptr_dum = lptr_opcg6_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg6_f_c_aer
            spectxt = 'opcg6_f_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg7_f_c_aer) then
            lptr_dum = lptr_opcg7_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg7_f_c_aer
            spectxt = 'opcg7_f_c_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg8_f_c_aer) then
            lptr_dum = lptr_opcg8_f_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg8_f_c_aer
            spectxt = 'opcg8_f_c_'
            nspectxt = 10
        else if (l_mastercomp .eq. mastercompindx_opcg1_f_o_aer) then
            lptr_dum = lptr_opcg1_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg1_f_o_aer
            spectxt = 'opcg1_f_o_'
            nspectxt = 10
        else if (l_mastercomp .eq. mastercompindx_opcg2_f_o_aer) then
            lptr_dum = lptr_opcg2_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg2_f_o_aer
            spectxt = 'opcg2_f_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg3_f_o_aer) then
            lptr_dum = lptr_opcg3_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg3_f_o_aer
            spectxt = 'opcg3_f_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg4_f_o_aer) then
            lptr_dum = lptr_opcg4_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg4_f_o_aer
            spectxt = 'opcg4_f_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg5_f_o_aer) then
            lptr_dum = lptr_opcg5_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg5_f_o_aer
            spectxt = 'opcg5_f_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg6_f_o_aer) then
            lptr_dum = lptr_opcg6_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg6_f_o_aer
            spectxt = 'opcg6_f_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg7_f_o_aer) then
            lptr_dum = lptr_opcg7_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg7_f_o_aer
            spectxt = 'opcg7_f_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_opcg8_f_o_aer) then
            lptr_dum = lptr_opcg8_f_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_opcg8_f_o_aer
            spectxt = 'opcg8_f_o_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_ant1_c_aer) then
            lptr_dum = lptr_ant1_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_ant1_c_aer
            spectxt = 'ant1_c_'
            nspectxt = 7
        else if (l_mastercomp .eq. mastercompindx_ant2_c_aer) then
            lptr_dum = lptr_ant2_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_ant2_c_aer
            spectxt = 'ant2_c_'
            nspectxt = 7
        else if (l_mastercomp .eq. mastercompindx_ant3_c_aer) then
            lptr_dum = lptr_ant3_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_ant3_c_aer
            spectxt = 'ant3_c_'
            nspectxt = 7
        else if (l_mastercomp .eq. mastercompindx_ant4_c_aer) then
            lptr_dum = lptr_ant4_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_ant4_c_aer
            spectxt = 'ant4_c_'
            nspectxt = 7
        else if (l_mastercomp .eq. mastercompindx_ant1_o_aer) then
            lptr_dum = lptr_ant1_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_ant1_o_aer
            spectxt = 'ant1_o_'
            nspectxt = 7
        else if (l_mastercomp .eq. mastercompindx_ant2_o_aer) then
            lptr_dum = lptr_ant2_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_ant2_o_aer
            spectxt = 'ant2_o_'
            nspectxt = 7
        else if (l_mastercomp .eq. mastercompindx_ant3_o_aer) then
            lptr_dum = lptr_ant3_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_ant3_o_aer
            spectxt = 'ant3_o_'
            nspectxt = 7
        else if (l_mastercomp .eq. mastercompindx_ant4_o_aer) then
            lptr_dum = lptr_ant4_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_ant4_o_aer
            spectxt = 'ant4_o_'
            nspectxt = 7
        else if (l_mastercomp .eq. mastercompindx_biog1_c_aer) then
            lptr_dum = lptr_biog1_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_biog1_c_aer
            spectxt = 'biog1_c_'
            nspectxt = 8
        else if (l_mastercomp .eq. mastercompindx_biog2_c_aer) then
            lptr_dum = lptr_biog2_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_biog2_c_aer
            spectxt = 'biog2_c_'
            nspectxt = 8
        else if (l_mastercomp .eq. mastercompindx_biog3_c_aer) then
            lptr_dum = lptr_biog3_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_biog3_c_aer
            spectxt = 'biog3_c_'
            nspectxt = 8
        else if (l_mastercomp .eq. mastercompindx_biog4_c_aer) then
            lptr_dum = lptr_biog4_c_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_biog4_c_aer
            spectxt = 'biog4_c_'
            nspectxt = 8
        else if (l_mastercomp .eq. mastercompindx_biog1_o_aer) then
            lptr_dum = lptr_biog1_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_biog1_o_aer
            spectxt = 'biog1_o_'
            nspectxt = 8
        else if (l_mastercomp .eq. mastercompindx_biog2_o_aer) then
            lptr_dum = lptr_biog2_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_biog2_o_aer
            spectxt = 'biog2_o_'
            nspectxt = 8
        else if (l_mastercomp .eq. mastercompindx_biog3_o_aer) then
            lptr_dum = lptr_biog3_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_biog3_o_aer
            spectxt = 'biog3_o_'
            nspectxt = 8
        else if (l_mastercomp .eq. mastercompindx_biog4_o_aer) then
            lptr_dum = lptr_biog4_o_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_biog4_o_aer
            spectxt = 'biog4_o_'
            nspectxt = 8

        else if (l_mastercomp .eq. mastercompindx_iepox_aer) then
            lptr_dum = lptr_iepox_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_iepox_aer
            spectxt = 'iepox_'
            nspectxt = 6
        else if (l_mastercomp .eq. mastercompindx_iepoxos_aer) then
            lptr_dum = lptr_iepoxos_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_iepoxos_aer
            spectxt = 'iepoxos_'
            nspectxt = 8
        else if (l_mastercomp .eq. mastercompindx_tetrol_aer) then
            lptr_dum = lptr_tetrol_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_tetrol_aer
            spectxt = 'tetrol_'
            nspectxt = 7

        else if (l_mastercomp .eq. mastercompindx_cldiepoxos_aer) then
            lptr_dum = lptr_cldiepoxos_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_cldiepoxos_aer
            spectxt = 'cldiepoxos_'
            nspectxt = 11
        else if (l_mastercomp .eq. mastercompindx_cldtetrol_aer) then
            lptr_dum = lptr_cldtetrol_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_cldtetrol_aer
            spectxt = 'cldtetrol_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_gly_aer) then
            lptr_dum = lptr_gly_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_gly_aer
            spectxt = 'gly_'
            nspectxt = 4

       else if (l_mastercomp .eq. mastercompindx_tanv_aer) then
            lptr_dum = lptr_tanv_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_tanv_aer
            spectxt = 'tanv_'
            nspectxt = 5

        else if (l_mastercomp .eq. mastercompindx_sopnv_aer) then
            lptr_dum = lptr_sopnv_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_sopnv_aer
            spectxt = 'sopnv_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_ternv_aer) then
            lptr_dum = lptr_ternv_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_ternv_aer
            spectxt = 'ternv_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_seqnv_aer) then
            lptr_dum = lptr_seqnv_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_seqnv_aer
            spectxt = 'seqnv_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_sianv_aer) then
            lptr_dum = lptr_sianv_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_sianv_aer
            spectxt = 'sianv_'
            nspectxt = 6


        else if (l_mastercomp .eq. mastercompindx_smpa_aer) then
            lptr_dum = lptr_smpa_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_smpa_aer
            spectxt = 'smpa_'
            nspectxt = 5

        else if (l_mastercomp .eq. mastercompindx_smpbb_aer) then
            lptr_dum = lptr_smpbb_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_smpbb_aer
            spectxt = 'smpbb_'
            nspectxt = 5

         else if (l_mastercomp .eq. mastercompindx_glysoa_r1_aer) then
            lptr_dum = lptr_glysoa_r1_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_glysoa_r1_aer
            spectxt = 'glysoa_r1_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_glysoa_r2_aer) then
            lptr_dum = lptr_glysoa_r2_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_glysoa_r2_aer
            spectxt = 'glysoa_r2_'
            nspectxt = 10

        else if (l_mastercomp .eq. mastercompindx_glysoa_sfc_aer) then
            lptr_dum = lptr_glysoa_sfc_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_glysoa_sfc_aer
            spectxt = 'glysoa_sfc_'
            nspectxt = 11

        else if (l_mastercomp .eq. mastercompindx_glysoa_nh4_aer) then
            lptr_dum = lptr_glysoa_nh4_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_glysoa_nh4_aer
            spectxt = 'glysoa_nh4_'
            nspectxt = 11

        else if (l_mastercomp .eq. mastercompindx_glysoa_oh_aer) then
            lptr_dum = lptr_glysoa_oh_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_glysoa_oh_aer
            spectxt = 'glysoa_oh_'
            nspectxt = 10


        else if (l_mastercomp .eq. mastercompindx_asoaX_aer) then
            lptr_dum = lptr_asoaX_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_asoaX_aer
            spectxt = 'asoaX_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_asoa1_aer) then
            lptr_dum = lptr_asoa1_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_asoa1_aer
            spectxt = 'asoa1_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_asoa2_aer) then
            lptr_dum = lptr_asoa2_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_asoa2_aer
            spectxt = 'asoa2_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_asoa3_aer) then
            lptr_dum = lptr_asoa3_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_asoa3_aer
            spectxt = 'asoa3_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_asoa4_aer) then
            lptr_dum = lptr_asoa4_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_asoa4_aer
            spectxt = 'asoa4_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_bsoaX_aer) then
            lptr_dum = lptr_bsoaX_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_bsoaX_aer
            spectxt = 'bsoaX_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_bsoa1_aer) then
            lptr_dum = lptr_bsoa1_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_bsoa1_aer
            spectxt = 'bsoa1_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_bsoa2_aer) then
            lptr_dum = lptr_bsoa2_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_bsoa2_aer
            spectxt = 'bsoa2_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_bsoa3_aer) then
            lptr_dum = lptr_bsoa3_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_bsoa3_aer
            spectxt = 'bsoa3_'
            nspectxt = 6

        else if (l_mastercomp .eq. mastercompindx_bsoa4_aer) then
            lptr_dum = lptr_bsoa4_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_bsoa4_aer
            spectxt = 'bsoa4_'
            nspectxt = 6
	
        else if (l_mastercomp .eq. mastercompindx_terp_Cn3_aer) then
            lptr_dum = lptr_terp_Cn3_aer(n,itype,iphase)
            mcindx_dum = mastercompindx_terp_Cn3_aer
            spectxt = 'terp_Cn3'
            nspectxt = 8

	else
	    goto 2500
	end if
	
	if (lptr_dum .gt. lmaxd) then

	    write( msg, '(a,3(1x,i4))' ) 'itype, isize, iphase =',   &
		itype, isize, iphase
	    call peg_message( lunout, msg )
	    write( msg, '(a,3(1x,i4))' ) 'l_mastercomp, lptr_dum, lmaxd =',   &
		l_mastercomp, lptr_dum, lmaxd
	    call peg_message( lunout, msg )
	    msg = '*** subr init_data_mosaic_ptr error - lptr_dum > lmaxd'
	    call peg_error_fatal( lunerr, msg )

	else if (lptr_dum .ge. p1st) then

	    ncomp_dum(isize,iphase) = ncomp_dum(isize,iphase) + iaddto_ncomp
	    ncomp_plustracer_dum(isize,iphase) =   &
		ncomp_plustracer_dum(isize,iphase) + iaddto_ncomp_plustracer

	    name(lptr_dum) =   &
		spectxt(1:nspectxt) // phasetxt(1:nphasetxt) //   &
		sizetxt(1:nsizetxt) //  typetxt(1:ntypetxt)

	    if (l_mastercomp .eq. -2) then

		mprognum_aer(n,itype,iphase) = 1

	    else if (l_mastercomp .eq. -1) then

		continue

	    else if (l_mastercomp .eq. 0) then

		continue

	    else if (l_mastercomp .gt. 0) then
		ll = ncomp_plustracer_dum(isize,iphase)

		massptr_aer(ll,n,itype,iphase) = lptr_dum
		mastercompptr_aer(ll,itype) = mcindx_dum

		name_aer(ll,itype) = name_mastercomp_aer(mcindx_dum)
		dens_aer(ll,itype) = dens_mastercomp_aer(mcindx_dum)
		mw_aer(ll,itype) = mw_mastercomp_aer(mcindx_dum)
		hygro_aer(ll,itype) = hygro_mastercomp_aer(mcindx_dum)

	    end if

	end if

2500	continue	

2600	continue	

2700	continue	




	ncomp_aer(itype) = ncomp_dum(1,ai_phase)
	ncomp_plustracer_aer(itype) = ncomp_plustracer_dum(1,ai_phase)

	do iphase = 1, nphase_aer
	do isize = 1, nsize_aer(itype)
	   PRINT *, 'iphase = ', iphase
	   PRINT *, 'itype = ', itype
	   PRINT *, 'isize = ', isize
	    if (ncomp_aer(itype) .ne. ncomp_dum(isize,iphase)) then
	        msg =  '*** subr init_data_mosaic_ptr - ' //   &
		    'ncomp_aer .ne. ncomp_dum'
		call peg_message( lunerr, msg )
		write(msg,9350) 'isize, itype, iphase =', isize, itype, iphase
		call peg_message( lunerr, msg )
		write(msg,9350) 'ncomp_aer, ncomp_dum =',   &
		    ncomp_aer(itype), ncomp_dum(isize,iphase)
		call peg_error_fatal( lunerr, msg )
	    end if
	    if (ncomp_plustracer_aer(itype) .ne.   &
			ncomp_plustracer_dum(isize,iphase)) then
	        msg = '*** subr init_data_mosaic_ptr - ' //   &
		    'ncomp_plustracer_aer .ne. ncomp_plustracer_dum'
		call peg_message( lunerr, msg )
		write(msg,9350) 'isize, itype, iphase =', isize, itype, iphase
		call peg_message( lunerr, msg )
		write(msg,9350)   &
		    'ncomp_plustracer_aer, ncomp_plustracer_dum =',   &
		    ncomp_plustracer_aer(itype),   &
		    ncomp_plustracer_dum(isize,iphase)
		call peg_error_fatal( lunerr, msg )
	    end if
	end do
	end do


2800	continue	


9320	format( a, i1, i1, a, 8x )




9350	format( a, 32(1x,i4) )
	msg = ' '
	call peg_message( lunout, msg )
	msg = 'output from subr init_data_mosaic_ptr'
	call peg_message( lunout, msg )
	write(msg,9350) 'nphase_aer =     ', nphase_aer
	call peg_message( lunout, msg )

	do iphase=1,nphase_aer

	write(msg,9350) 'iphase =     ', iphase
	call peg_message( lunout, msg )
	write(msg,9350) 'ntype_aer =     ', ntype_aer
	call peg_message( lunout, msg )

        write(msg,9350) 'ncomp_aer =     ', ncomp_aer
        call peg_message( lunout, msg )

	do itype=1,ntype_aer

	write(msg,9350) 'itype =     ', itype
	call peg_message( lunout, msg )
	write(msg,9350) 'nsize_aer = ', nsize_aer(itype)
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_so4_aer ',   &
		(lptr_so4_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_no3_aer ',   &
		(lptr_no3_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_cl_aer  ',   &
		(lptr_cl_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_msa_aer ',   &
		(lptr_msa_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_co3_aer ',   &
		(lptr_co3_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_nh4_aer ',   &
		(lptr_nh4_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_na_aer  ',   &
		(lptr_na_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_ca_aer  ',   &
		(lptr_ca_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_oin_aer ',   &
		(lptr_oin_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_oc_aer  ',   &
		(lptr_oc_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'lptr_bc_aer  ',   &
		(lptr_bc_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'hyswptr_aer',   &
		(hyswptr_aer(n,itype), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'waterptr_aer  ',   &
		(waterptr_aer(n,itype), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )
	write(msg,9350) 'numptr_aer     ',   &
		(numptr_aer(n,itype,iphase), n=1,nsize_aer(itype))
	call peg_message( lunout, msg )

        write(msg,9350) 'lptr_tanv_aer  ',   &
                (lptr_tanv_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_sopnv_aer  ',   &
                (lptr_sopnv_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_ternv_aer  ',   &
                (lptr_ternv_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_seqnv_aer  ',   &
                (lptr_seqnv_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_sianv_aer  ',   &
                (lptr_sianv_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )

	

        write(msg,9350) 'lptr_pcg1_b_c_aer ',   &
                (lptr_pcg1_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg2_b_c_aer ',   &
                (lptr_pcg2_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg3_b_c_aer ',   &
                (lptr_pcg3_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg4_b_c_aer ',   &
                (lptr_pcg4_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg5_b_c_aer ',   &
                (lptr_pcg5_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg6_b_c_aer ',   &
                (lptr_pcg6_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg7_b_c_aer ',   &
                (lptr_pcg7_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg8_b_c_aer ',   &
                (lptr_pcg8_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg9_b_c_aer ',   &
                (lptr_pcg9_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg1_b_o_aer ',   &
                (lptr_pcg1_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg2_b_o_aer ',   &
                (lptr_pcg2_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg3_b_o_aer ',   &
                (lptr_pcg3_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg4_b_o_aer ',   &
                (lptr_pcg4_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg5_b_o_aer ',   &
                (lptr_pcg5_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg6_b_o_aer ',   &
                (lptr_pcg6_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg7_b_o_aer ',   &
                (lptr_pcg7_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg8_b_o_aer ',   &
                (lptr_pcg8_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg9_b_o_aer ',   &
                (lptr_pcg9_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg1_b_c_aer ',   &
                (lptr_opcg1_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg2_b_c_aer ',   &
                (lptr_opcg2_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg3_b_c_aer ',   &
                (lptr_opcg3_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg4_b_c_aer ',   &
                (lptr_opcg4_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg5_b_c_aer ',   &
                (lptr_opcg5_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg6_b_c_aer ',   &
                (lptr_opcg6_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg7_b_c_aer ',   &
                (lptr_opcg7_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg8_b_c_aer ',   &
                (lptr_opcg8_b_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg1_b_o_aer ',   &
                (lptr_opcg1_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg2_b_o_aer ',   &
                (lptr_opcg2_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg3_b_o_aer ',   &
                (lptr_opcg3_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg4_b_o_aer ',   &
                (lptr_opcg4_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg5_b_o_aer ',   &
                (lptr_opcg5_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg6_b_o_aer ',   &
                (lptr_opcg6_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg7_b_o_aer ',   &
                (lptr_opcg7_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg8_b_o_aer ',   &
                (lptr_opcg8_b_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg1_f_c_aer ',   &
                (lptr_pcg1_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg2_f_c_aer ',   &
                (lptr_pcg2_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg3_f_c_aer ',   &
                (lptr_pcg3_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg4_f_c_aer ',   &
                (lptr_pcg4_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg5_f_c_aer ',   &
                (lptr_pcg5_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg6_f_c_aer ',   &
                (lptr_pcg6_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg7_f_c_aer ',   &
                (lptr_pcg7_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg8_f_c_aer ',   &
                (lptr_pcg8_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg9_f_c_aer ',   &
                (lptr_pcg9_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg1_f_o_aer ',   &
                (lptr_pcg1_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg2_f_o_aer ',   &
                (lptr_pcg2_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg3_f_o_aer ',   &
                (lptr_pcg3_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg4_f_o_aer ',   &
                (lptr_pcg4_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg5_f_o_aer ',   &
                (lptr_pcg5_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg6_f_o_aer ',   &
                (lptr_pcg6_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg7_f_o_aer ',   &
                (lptr_pcg7_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg8_f_o_aer ',   &
                (lptr_pcg8_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_pcg9_f_o_aer ',   &
                (lptr_pcg9_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg1_f_c_aer ',   &
                (lptr_opcg1_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg2_f_c_aer ',   &
                (lptr_opcg2_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg3_f_c_aer ',   &
                (lptr_opcg3_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg4_f_c_aer ',   &
                (lptr_opcg4_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg5_f_c_aer ',   &
                (lptr_opcg5_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg6_f_c_aer ',   &
                (lptr_opcg6_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg7_f_c_aer ',   &
                (lptr_opcg7_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg8_f_c_aer ',   &
                (lptr_opcg8_f_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg1_f_o_aer ',   &
                (lptr_opcg1_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg2_f_o_aer ',   &
                (lptr_opcg2_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg3_f_o_aer ',   &
                (lptr_opcg3_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg4_f_o_aer ',   &
                (lptr_opcg4_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg5_f_o_aer ',   &
                (lptr_opcg5_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg6_f_o_aer ',   &
                (lptr_opcg6_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg7_f_o_aer ',   &
                (lptr_opcg7_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'lptr_opcg8_f_o_aer ',   &
                (lptr_opcg8_f_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )

        write(msg,9350) 'ant1_c_aer ',   &
                (lptr_ant1_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'ant2_c_aer ',   &
                (lptr_ant2_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'ant3_c_aer ',   &
                (lptr_ant3_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'ant4_c_aer ',   &
                (lptr_ant4_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'ant1_o_aer ',   &
                (lptr_ant1_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'ant2_o_aer ',   &
                (lptr_ant2_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'ant3_o_aer ',   &
                (lptr_ant3_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'ant4_o_aer ',   &
                (lptr_ant4_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'biog1_c_aer ',   &
                (lptr_biog1_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'biog2_c_aer ',   &
                (lptr_biog2_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'biog3_c_aer ',   &
                (lptr_biog3_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'biog4_c_aer ',   &
                (lptr_biog4_c_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'biog1_o_aer ',   &
                (lptr_biog1_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'biog2_o_aer ',   &
                (lptr_biog2_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'biog3_o_aer ',   &
                (lptr_biog3_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'biog4_o_aer ',   &
                (lptr_biog4_o_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )

        write(msg,9350) 'iepox_aer ',   &
                (lptr_iepox_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )

        write(msg,9350) 'iepoxos_aer ',   &
                (lptr_iepoxos_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )

        write(msg,9350) 'tetrol_aer ',   &
                (lptr_tetrol_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )

        write(msg,9350) 'cldiepoxos_aer ',   &
                (lptr_cldiepoxos_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )

        write(msg,9350) 'cldtetrol_aer ',   &
                (lptr_cldtetrol_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )

        write(msg,9350) 'gly_aer ',   &
                (lptr_gly_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )


        write(msg,9350) 'smpa_aer ',   &
                (lptr_smpa_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'smpbb_aer ',   &
                (lptr_smpbb_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'glysoa_r1_aer ',   &
                (lptr_glysoa_r1_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'glysoa_r2_aer ',   &
                (lptr_glysoa_r2_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'glysoa_sfc_aer ',   &
                (lptr_glysoa_sfc_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'glysoa_nh4_aer ',   &
                (lptr_glysoa_nh4_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'glysoa_oh_aer ',   &
                (lptr_glysoa_oh_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'asoaX_aer ',   &
                (lptr_asoaX_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'asoa1_aer ',   &
                (lptr_asoa1_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'asoa2_aer ',   &
                (lptr_asoa2_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'asoa3_aer ',   &
                (lptr_asoa3_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'asoa4_aer ',   &
                (lptr_asoa4_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )

        write(msg,9350) 'bsoaX_aer ',   &
                (lptr_bsoaX_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'bsoa1_aer ',   &
                (lptr_bsoa1_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'bsoa2_aer ',   &
                (lptr_bsoa2_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'bsoa3_aer ',   &
                (lptr_bsoa3_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
        write(msg,9350) 'bsoa4_aer ',   &
                (lptr_bsoa4_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )
	
	write(msg,9350) 'terp_Cn3_aer ',   &
	        (lptr_terp_Cn3_aer(n,itype,iphase), n=1,nsize_aer(itype))
        call peg_message( lunout, msg )


  write(*,*) " ---------------------- "




	do ll = 1, ncomp_plustracer_aer(itype)


	    write(msg,9350) name_aer(ll,itype),   &
		(massptr_aer(ll,n,itype,iphase), n=1,nsize_aer(itype)), ll
	    call peg_message( lunout, msg )
	end do
	end do 
	end do 




	do iphase=1,nphase_aer
	do itype=1,ntype_aer
	y_so4 = 0
	y_no3 = 0
	y_cl = 0
	y_msa = 0
	y_co3 = 0
	y_nh4 = 0
	y_na = 0
	y_ca = 0
	y_oin = 0
	y_oc = 0
	y_bc = 0
	y_hysw = 0
	y_water = 0
	y_num = 0
        y_tanv = 0
        y_sopnv = 0
        y_ternv = 0
        y_seqnv = 0
        y_sianv = 0

        y_pcg1_b_c=0
        y_pcg2_b_c=0
        y_pcg3_b_c=0
        y_pcg4_b_c=0
        y_pcg5_b_c=0
        y_pcg6_b_c=0
        y_pcg7_b_c=0
        y_pcg8_b_c=0
        y_pcg9_b_c=0
        y_pcg1_b_o=0
        y_pcg2_b_o=0
        y_pcg3_b_o=0
        y_pcg4_b_o=0
        y_pcg5_b_o=0
        y_pcg6_b_o=0
        y_pcg7_b_o=0
        y_pcg8_b_o=0
        y_pcg9_b_o=0
        y_opcg1_b_c=0
        y_opcg2_b_c=0
        y_opcg3_b_c=0
        y_opcg4_b_c=0
        y_opcg5_b_c=0
        y_opcg6_b_c=0
        y_opcg7_b_c=0
        y_opcg8_b_c=0
        y_opcg1_b_o=0
        y_opcg2_b_o=0
        y_opcg3_b_o=0
        y_opcg4_b_o=0
        y_opcg5_b_o=0
        y_opcg6_b_o=0
        y_opcg7_b_o=0
        y_opcg8_b_o=0
        y_pcg1_f_c=0
        y_pcg2_f_c=0
        y_pcg3_f_c=0
        y_pcg4_f_c=0
        y_pcg5_f_c=0
        y_pcg6_f_c=0
        y_pcg7_f_c=0
        y_pcg8_f_c=0
        y_pcg9_f_c=0
        y_pcg1_f_o=0
        y_pcg2_f_o=0
        y_pcg3_f_o=0
        y_pcg4_f_o=0
        y_pcg5_f_o=0
        y_pcg6_f_o=0
        y_pcg7_f_o=0
        y_pcg8_f_o=0
        y_pcg9_f_o=0
        y_opcg1_f_c=0
        y_opcg2_f_c=0
        y_opcg3_f_c=0
        y_opcg4_f_c=0
        y_opcg5_f_c=0
        y_opcg6_f_c=0
        y_opcg7_f_c=0
        y_opcg8_f_c=0
        y_opcg1_f_o=0
        y_opcg2_f_o=0
        y_opcg3_f_o=0
        y_opcg4_f_o=0
        y_opcg5_f_o=0
        y_opcg6_f_o=0
        y_opcg7_f_o=0
        y_opcg8_f_o=0
        y_ant1_c=0
        y_ant2_c=0
        y_ant3_c=0
        y_ant4_c=0
        y_ant1_o=0
        y_ant2_o=0
        y_ant3_o=0
        y_ant4_o=0
        y_biog1_c=0
        y_biog2_c=0
        y_biog3_c=0
        y_biog4_c=0
        y_biog1_o=0
        y_biog2_o=0
        y_biog3_o=0
        y_biog4_o=0
        y_iepox = 0
        y_iepoxos = 0
        y_tetrol = 0
        y_cldiepoxos = 0
        y_cldtetrol = 0
        y_gly = 0
        y_smpa=0
        y_smpbb=0
        y_glysoa_r1=0
        y_glysoa_r2=0
        y_glysoa_sfc=0
        y_glysoa_nh4=0
        y_glysoa_oh=0
        y_asoaX=0
        y_asoa1=0
        y_asoa2=0
        y_asoa3=0
        y_asoa4=0
        y_bsoaX=0
        y_bsoa1=0
        y_bsoa2=0
        y_bsoa3=0
        y_bsoa4=0
	
	y_terp_Cn3=0
    

	do n = 1, nsize_aer(itype)
	    if (lptr_so4_aer(n,itype,iphase) .ge. p1st) y_so4 = y_so4 + 1
	    if (lptr_no3_aer(n,itype,iphase) .ge. p1st) y_no3 = y_no3 + 1
	    if (lptr_cl_aer(n,itype,iphase)  .ge. p1st) y_cl  = y_cl + 1
	    if (lptr_msa_aer(n,itype,iphase) .ge. p1st) y_msa = y_msa + 1
	    if (lptr_co3_aer(n,itype,iphase) .ge. p1st) y_co3 = y_co3 + 1
	    if (lptr_nh4_aer(n,itype,iphase) .ge. p1st) y_nh4 = y_nh4 + 1
	    if (lptr_na_aer(n,itype,iphase)  .ge. p1st) y_na  = y_na + 1
	    if (lptr_ca_aer(n,itype,iphase)  .ge. p1st) y_ca  = y_ca + 1
	    if (lptr_oin_aer(n,itype,iphase) .ge. p1st) y_oin = y_oin + 1
	    if (lptr_oc_aer(n,itype,iphase)  .ge. p1st) y_oc  = y_oc + 1
	    if (lptr_bc_aer(n,itype,iphase)  .ge. p1st) y_bc  = y_bc + 1
	    if (hyswptr_aer(n,itype)    .ge. p1st) y_hysw = y_hysw + 1
	    if (waterptr_aer(n,itype)        .ge. p1st) y_water = y_water + 1
	    if (numptr_aer(n,itype,iphase)   .ge. p1st) y_num = y_num + 1
            if (lptr_tanv_aer(n,itype,iphase)  .ge. p1st) y_tanv  = y_tanv + 1
            if (lptr_sopnv_aer(n,itype,iphase)  .ge. p1st) y_sopnv  = y_sopnv + 1
            if (lptr_ternv_aer(n,itype,iphase)  .ge. p1st) y_ternv  = y_ternv + 1
            if (lptr_seqnv_aer(n,itype,iphase)  .ge. p1st) y_seqnv  = y_seqnv + 1
            if (lptr_sianv_aer(n,itype,iphase)  .ge. p1st) y_sianv  = y_sianv + 1
            if (lptr_pcg1_b_c_aer(n,itype,iphase) .ge. p1st) y_pcg1_b_c = y_pcg1_b_c + 1
            if (lptr_pcg2_b_c_aer(n,itype,iphase) .ge. p1st) y_pcg2_b_c = y_pcg2_b_c + 1
            if (lptr_pcg3_b_c_aer(n,itype,iphase) .ge. p1st) y_pcg3_b_c = y_pcg3_b_c + 1
            if (lptr_pcg4_b_c_aer(n,itype,iphase) .ge. p1st) y_pcg4_b_c = y_pcg4_b_c + 1
            if (lptr_pcg5_b_c_aer(n,itype,iphase) .ge. p1st) y_pcg5_b_c = y_pcg5_b_c + 1
            if (lptr_pcg6_b_c_aer(n,itype,iphase) .ge. p1st) y_pcg6_b_c = y_pcg6_b_c + 1
            if (lptr_pcg7_b_c_aer(n,itype,iphase) .ge. p1st) y_pcg7_b_c = y_pcg7_b_c + 1
            if (lptr_pcg8_b_c_aer(n,itype,iphase) .ge. p1st) y_pcg8_b_c = y_pcg8_b_c + 1
            if (lptr_pcg9_b_c_aer(n,itype,iphase) .ge. p1st) y_pcg9_b_c = y_pcg9_b_c + 1
            if (lptr_pcg1_b_o_aer(n,itype,iphase) .ge. p1st) y_pcg1_b_o = y_pcg1_b_o + 1
            if (lptr_pcg2_b_o_aer(n,itype,iphase) .ge. p1st) y_pcg2_b_o = y_pcg2_b_o + 1
            if (lptr_pcg3_b_o_aer(n,itype,iphase) .ge. p1st) y_pcg3_b_o = y_pcg3_b_o + 1
            if (lptr_pcg4_b_o_aer(n,itype,iphase) .ge. p1st) y_pcg4_b_o = y_pcg4_b_o + 1
            if (lptr_pcg5_b_o_aer(n,itype,iphase) .ge. p1st) y_pcg5_b_o = y_pcg5_b_o + 1
            if (lptr_pcg6_b_o_aer(n,itype,iphase) .ge. p1st) y_pcg6_b_o = y_pcg6_b_o + 1
            if (lptr_pcg7_b_o_aer(n,itype,iphase) .ge. p1st) y_pcg7_b_o = y_pcg7_b_o + 1
            if (lptr_pcg8_b_o_aer(n,itype,iphase) .ge. p1st) y_pcg8_b_o = y_pcg8_b_o + 1
            if (lptr_pcg9_b_o_aer(n,itype,iphase) .ge. p1st) y_pcg9_b_o = y_pcg9_b_o + 1
            if (lptr_opcg1_b_c_aer(n,itype,iphase) .ge. p1st) y_opcg1_b_c = y_opcg1_b_c + 1
            if (lptr_opcg2_b_c_aer(n,itype,iphase) .ge. p1st) y_opcg2_b_c = y_opcg2_b_c + 1
            if (lptr_opcg3_b_c_aer(n,itype,iphase) .ge. p1st) y_opcg3_b_c = y_opcg3_b_c + 1
            if (lptr_opcg4_b_c_aer(n,itype,iphase) .ge. p1st) y_opcg4_b_c = y_opcg4_b_c + 1
            if (lptr_opcg5_b_c_aer(n,itype,iphase) .ge. p1st) y_opcg5_b_c = y_opcg5_b_c + 1
            if (lptr_opcg6_b_c_aer(n,itype,iphase) .ge. p1st) y_opcg6_b_c = y_opcg6_b_c + 1
            if (lptr_opcg7_b_c_aer(n,itype,iphase) .ge. p1st) y_opcg7_b_c = y_opcg7_b_c + 1
            if (lptr_opcg8_b_c_aer(n,itype,iphase) .ge. p1st) y_opcg8_b_c = y_opcg8_b_c + 1
            if (lptr_opcg1_b_o_aer(n,itype,iphase) .ge. p1st) y_opcg1_b_o = y_opcg1_b_o + 1
            if (lptr_opcg2_b_o_aer(n,itype,iphase) .ge. p1st) y_opcg2_b_o = y_opcg2_b_o + 1
            if (lptr_opcg3_b_o_aer(n,itype,iphase) .ge. p1st) y_opcg3_b_o = y_opcg3_b_o + 1
            if (lptr_opcg4_b_o_aer(n,itype,iphase) .ge. p1st) y_opcg4_b_o = y_opcg4_b_o + 1
            if (lptr_opcg5_b_o_aer(n,itype,iphase) .ge. p1st) y_opcg5_b_o = y_opcg5_b_o + 1
            if (lptr_opcg6_b_o_aer(n,itype,iphase) .ge. p1st) y_opcg6_b_o = y_opcg6_b_o + 1
            if (lptr_opcg7_b_o_aer(n,itype,iphase) .ge. p1st) y_opcg7_b_o = y_opcg7_b_o + 1
            if (lptr_opcg8_b_o_aer(n,itype,iphase) .ge. p1st) y_opcg8_b_o = y_opcg8_b_o + 1
            if (lptr_pcg1_f_c_aer(n,itype,iphase) .ge. p1st) y_pcg1_f_c = y_pcg1_f_c + 1
            if (lptr_pcg2_f_c_aer(n,itype,iphase) .ge. p1st) y_pcg2_f_c = y_pcg2_f_c + 1
            if (lptr_pcg3_f_c_aer(n,itype,iphase) .ge. p1st) y_pcg3_f_c = y_pcg3_f_c + 1
            if (lptr_pcg4_f_c_aer(n,itype,iphase) .ge. p1st) y_pcg4_f_c = y_pcg4_f_c + 1
            if (lptr_pcg5_f_c_aer(n,itype,iphase) .ge. p1st) y_pcg5_f_c = y_pcg5_f_c + 1
            if (lptr_pcg6_f_c_aer(n,itype,iphase) .ge. p1st) y_pcg6_f_c = y_pcg6_f_c + 1
            if (lptr_pcg7_f_c_aer(n,itype,iphase) .ge. p1st) y_pcg7_f_c = y_pcg7_f_c + 1
            if (lptr_pcg8_f_c_aer(n,itype,iphase) .ge. p1st) y_pcg8_f_c = y_pcg8_f_c + 1
            if (lptr_pcg9_f_c_aer(n,itype,iphase) .ge. p1st) y_pcg9_f_c = y_pcg9_f_c + 1
            if (lptr_pcg1_f_o_aer(n,itype,iphase) .ge. p1st) y_pcg1_f_o = y_pcg1_f_o + 1
            if (lptr_pcg2_f_o_aer(n,itype,iphase) .ge. p1st) y_pcg2_f_o = y_pcg2_f_o + 1
            if (lptr_pcg3_f_o_aer(n,itype,iphase) .ge. p1st) y_pcg3_f_o = y_pcg3_f_o + 1
            if (lptr_pcg4_f_o_aer(n,itype,iphase) .ge. p1st) y_pcg4_f_o = y_pcg4_f_o + 1
            if (lptr_pcg5_f_o_aer(n,itype,iphase) .ge. p1st) y_pcg5_f_o = y_pcg5_f_o + 1
            if (lptr_pcg6_f_o_aer(n,itype,iphase) .ge. p1st) y_pcg6_f_o = y_pcg6_f_o + 1
            if (lptr_pcg7_f_o_aer(n,itype,iphase) .ge. p1st) y_pcg7_f_o = y_pcg7_f_o + 1
            if (lptr_pcg8_f_o_aer(n,itype,iphase) .ge. p1st) y_pcg8_f_o = y_pcg8_f_o + 1
            if (lptr_pcg9_f_o_aer(n,itype,iphase) .ge. p1st) y_pcg9_f_o = y_pcg9_f_o + 1
            if (lptr_opcg1_f_c_aer(n,itype,iphase) .ge. p1st) y_opcg1_f_c = y_opcg1_f_c + 1
            if (lptr_opcg2_f_c_aer(n,itype,iphase) .ge. p1st) y_opcg2_f_c = y_opcg2_f_c + 1
            if (lptr_opcg3_f_c_aer(n,itype,iphase) .ge. p1st) y_opcg3_f_c = y_opcg3_f_c + 1
            if (lptr_opcg4_f_c_aer(n,itype,iphase) .ge. p1st) y_opcg4_f_c = y_opcg4_f_c + 1
            if (lptr_opcg5_f_c_aer(n,itype,iphase) .ge. p1st) y_opcg5_f_c = y_opcg5_f_c + 1
            if (lptr_opcg6_f_c_aer(n,itype,iphase) .ge. p1st) y_opcg6_f_c = y_opcg6_f_c + 1
            if (lptr_opcg7_f_c_aer(n,itype,iphase) .ge. p1st) y_opcg7_f_c = y_opcg7_f_c + 1
            if (lptr_opcg8_f_c_aer(n,itype,iphase) .ge. p1st) y_opcg8_f_c = y_opcg8_f_c + 1
            if (lptr_opcg1_f_o_aer(n,itype,iphase) .ge. p1st) y_opcg1_f_o = y_opcg1_f_o + 1
            if (lptr_opcg2_f_o_aer(n,itype,iphase) .ge. p1st) y_opcg2_f_o = y_opcg2_f_o + 1
            if (lptr_opcg3_f_o_aer(n,itype,iphase) .ge. p1st) y_opcg3_f_o = y_opcg3_f_o + 1
            if (lptr_opcg4_f_o_aer(n,itype,iphase) .ge. p1st) y_opcg4_f_o = y_opcg4_f_o + 1
            if (lptr_opcg5_f_o_aer(n,itype,iphase) .ge. p1st) y_opcg5_f_o = y_opcg5_f_o + 1
            if (lptr_opcg6_f_o_aer(n,itype,iphase) .ge. p1st) y_opcg6_f_o = y_opcg6_f_o + 1
            if (lptr_opcg7_f_o_aer(n,itype,iphase) .ge. p1st) y_opcg7_f_o = y_opcg7_f_o + 1
            if (lptr_opcg8_f_o_aer(n,itype,iphase) .ge. p1st) y_opcg8_f_o = y_opcg8_f_o + 1
            if (lptr_ant1_c_aer(n,itype,iphase) .ge. p1st) y_ant1_c = y_ant1_c + 1
            if (lptr_ant2_c_aer(n,itype,iphase) .ge. p1st) y_ant2_c = y_ant2_c + 1
            if (lptr_ant3_c_aer(n,itype,iphase) .ge. p1st) y_ant3_c = y_ant3_c + 1
            if (lptr_ant4_c_aer(n,itype,iphase) .ge. p1st) y_ant4_c = y_ant4_c + 1
            if (lptr_ant1_o_aer(n,itype,iphase) .ge. p1st) y_ant1_o = y_ant1_o + 1
            if (lptr_ant2_o_aer(n,itype,iphase) .ge. p1st) y_ant2_o = y_ant2_o + 1
            if (lptr_ant3_o_aer(n,itype,iphase) .ge. p1st) y_ant3_o = y_ant3_o + 1
            if (lptr_ant4_o_aer(n,itype,iphase) .ge. p1st) y_ant4_o = y_ant4_o + 1
            if (lptr_biog1_c_aer(n,itype,iphase) .ge. p1st) y_biog1_c = y_biog1_c + 1
            if (lptr_biog2_c_aer(n,itype,iphase) .ge. p1st) y_biog2_c = y_biog2_c + 1
            if (lptr_biog3_c_aer(n,itype,iphase) .ge. p1st) y_biog3_c = y_biog3_c + 1
            if (lptr_biog4_c_aer(n,itype,iphase) .ge. p1st) y_biog4_c = y_biog4_c + 1
            if (lptr_biog1_o_aer(n,itype,iphase) .ge. p1st) y_biog1_o = y_biog1_o + 1
            if (lptr_biog2_o_aer(n,itype,iphase) .ge. p1st) y_biog2_o = y_biog2_o + 1
            if (lptr_biog3_o_aer(n,itype,iphase) .ge. p1st) y_biog3_o = y_biog3_o + 1
            if (lptr_biog4_o_aer(n,itype,iphase) .ge. p1st) y_biog4_o = y_biog4_o + 1
            if (lptr_iepox_aer(n,itype,iphase) .ge. p1st) y_iepox = y_iepox + 1
            if (lptr_iepoxos_aer(n,itype,iphase) .ge. p1st) y_iepoxos = y_iepoxos + 1
            if (lptr_tetrol_aer(n,itype,iphase) .ge. p1st) y_tetrol = y_tetrol + 1
            if (lptr_cldiepoxos_aer(n,itype,iphase) .ge. p1st) y_cldiepoxos = y_cldiepoxos + 1
            if (lptr_cldtetrol_aer(n,itype,iphase) .ge. p1st) y_cldtetrol = y_cldtetrol + 1
            if (lptr_gly_aer(n,itype,iphase) .ge. p1st) y_gly = y_gly + 1
            if (lptr_smpa_aer(n,itype,iphase) .ge. p1st) y_smpa = y_smpa + 1
            if (lptr_smpbb_aer(n,itype,iphase) .ge. p1st) y_smpbb = y_smpbb + 1
            if (lptr_glysoa_r1_aer(n,itype,iphase) .ge. p1st) y_glysoa_r1 = y_glysoa_r1 + 1
            if (lptr_glysoa_r2_aer(n,itype,iphase) .ge. p1st) y_glysoa_r2 = y_glysoa_r2 + 1
            if (lptr_glysoa_sfc_aer(n,itype,iphase) .ge. p1st) y_glysoa_sfc = y_glysoa_sfc + 1
            if (lptr_glysoa_nh4_aer(n,itype,iphase) .ge. p1st) y_glysoa_nh4 = y_glysoa_nh4 + 1
            if (lptr_glysoa_oh_aer(n,itype,iphase) .ge. p1st) y_glysoa_oh = y_glysoa_oh + 1
            if (lptr_asoaX_aer(n,itype,iphase) .ge. p1st) y_asoaX = y_asoaX + 1
            if (lptr_asoa1_aer(n,itype,iphase) .ge. p1st) y_asoa1 = y_asoa1 + 1
            if (lptr_asoa2_aer(n,itype,iphase) .ge. p1st) y_asoa2 = y_asoa2 + 1
            if (lptr_asoa3_aer(n,itype,iphase) .ge. p1st) y_asoa3 = y_asoa3 + 1
            if (lptr_asoa4_aer(n,itype,iphase) .ge. p1st) y_asoa4 = y_asoa4 + 1
            if (lptr_bsoaX_aer(n,itype,iphase) .ge. p1st) y_bsoaX = y_bsoaX + 1
            if (lptr_bsoa1_aer(n,itype,iphase) .ge. p1st) y_bsoa1 = y_bsoa1 + 1
            if (lptr_bsoa2_aer(n,itype,iphase) .ge. p1st) y_bsoa2 = y_bsoa2 + 1
            if (lptr_bsoa3_aer(n,itype,iphase) .ge. p1st) y_bsoa3 = y_bsoa3 + 1
            if (lptr_bsoa4_aer(n,itype,iphase) .ge. p1st) y_bsoa4 = y_bsoa4 + 1
	    
	    if (lptr_terp_Cn3_aer(n,itype,iphase) .ge. p1st) y_terp_Cn3 = y_terp_Cn3 + 1

	end do


	if (y_so4 .ne. nsize_aer(itype)) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for so4'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
	else if (y_water .ne. nsize_aer(itype)) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for water'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
	else if (y_num .ne. nsize_aer(itype)) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for num'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
	end if




	if      ((y_no3 .ne. 0) .and.   &
	         (y_no3 .ne. nsize_aer(itype))) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for no3'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
        else if ((y_cl .ne. 0) .and.   &
                 (y_cl .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for cl'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_pcg1_b_c .ne. 0) .and.   &
                 (y_pcg1_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg1_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg2_b_c .ne. 0) .and.   &
                 (y_pcg2_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg2_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg3_b_c .ne. 0) .and.   &
                 (y_pcg3_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg3_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg4_b_c .ne. 0) .and.   &
                 (y_pcg4_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg4_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg5_b_c .ne. 0) .and.   &
                 (y_pcg5_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg5_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg6_b_c .ne. 0) .and.   &
                 (y_pcg6_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg6_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg7_b_c .ne. 0) .and.   &
                 (y_pcg7_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg7_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg8_b_c .ne. 0) .and.   &
                 (y_pcg8_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg8_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg9_b_c .ne. 0) .and.   &
                 (y_pcg9_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg9_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg1_b_o .ne. 0) .and.   &
                 (y_pcg1_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg1_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg2_b_o .ne. 0) .and.   &
                 (y_pcg2_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg2_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg3_b_o .ne. 0) .and.   &
                 (y_pcg3_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg3_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg4_b_o .ne. 0) .and.   &
                 (y_pcg4_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg4_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg5_b_o .ne. 0) .and.   &
                 (y_pcg5_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg5_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg6_b_o .ne. 0) .and.   &
                 (y_pcg6_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg6_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg7_b_o .ne. 0) .and.   &
                 (y_pcg7_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg7_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg8_b_o .ne. 0) .and.   &
                 (y_pcg8_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg8_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg9_b_o .ne. 0) .and.   &
                 (y_pcg9_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg9_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg1_b_o .ne. 0) .and.   &
                 (y_opcg1_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg1_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg2_b_o .ne. 0) .and.   &
                 (y_opcg2_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg2_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg3_b_o .ne. 0) .and.   &
                 (y_opcg3_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg3_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg4_b_o .ne. 0) .and.   &
                 (y_opcg4_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg4_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg5_b_o .ne. 0) .and.   &
                 (y_opcg5_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg5_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg6_b_o .ne. 0) .and.   &
                 (y_opcg6_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg6_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg7_b_o .ne. 0) .and.   &
                 (y_opcg7_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg7_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg8_b_o .ne. 0) .and.   &
                 (y_opcg8_b_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg8_b_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg1_b_c .ne. 0) .and.   &
                 (y_opcg1_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg1_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg2_b_c .ne. 0) .and.   &
                 (y_opcg2_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg2_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg3_b_c .ne. 0) .and.   &
                 (y_opcg3_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg3_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg4_b_c .ne. 0) .and.   &
                 (y_opcg4_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg4_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg5_b_c .ne. 0) .and.   &
                 (y_opcg5_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg5_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg6_b_c .ne. 0) .and.   &
                 (y_opcg6_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg6_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg7_b_c .ne. 0) .and.   &
                 (y_opcg7_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg7_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg8_b_c .ne. 0) .and.   &
                 (y_opcg8_b_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg8_b_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg1_f_c .ne. 0) .and.   &
                 (y_pcg1_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg1_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg2_f_c .ne. 0) .and.   &
                 (y_pcg2_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg2_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg3_f_c .ne. 0) .and.   &
                 (y_pcg3_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg3_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg4_f_c .ne. 0) .and.   &
                 (y_pcg4_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg4_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg5_f_c .ne. 0) .and.   &
                 (y_pcg5_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg5_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg6_f_c .ne. 0) .and.   &
                 (y_pcg6_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg6_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg7_f_c .ne. 0) .and.   &
                 (y_pcg7_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg7_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg8_f_c .ne. 0) .and.   &
                 (y_pcg8_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg8_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg9_f_c .ne. 0) .and.   &
                 (y_pcg9_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg9_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg1_f_o .ne. 0) .and.   &
                 (y_pcg1_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg1_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg2_f_o .ne. 0) .and.   &
                 (y_pcg2_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg2_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg3_f_o .ne. 0) .and.   &
                 (y_pcg3_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg3_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg4_f_o .ne. 0) .and.   &
                 (y_pcg4_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg4_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg5_f_o .ne. 0) .and.   &
                 (y_pcg5_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg5_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg6_f_o .ne. 0) .and.   &
                 (y_pcg6_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg6_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg7_f_o .ne. 0) .and.   &
                 (y_pcg7_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg7_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg8_f_o .ne. 0) .and.   &
                 (y_pcg8_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg8_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_pcg9_f_o .ne. 0) .and.   &
                 (y_pcg9_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for pcg9_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg1_f_o .ne. 0) .and.   &
                 (y_opcg1_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg1_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg2_f_o .ne. 0) .and.   &
                 (y_opcg2_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg2_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg3_f_o .ne. 0) .and.   &
                 (y_opcg3_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg3_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg4_f_o .ne. 0) .and.   &
                 (y_opcg4_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg4_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg5_f_o .ne. 0) .and.   &
                 (y_opcg5_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg5_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg6_f_o .ne. 0) .and.   &
                 (y_opcg6_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg6_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg7_f_o .ne. 0) .and.   &
                 (y_opcg7_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg7_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg8_f_o .ne. 0) .and.   &
                 (y_opcg8_f_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg8_f_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg1_f_c .ne. 0) .and.   &
                 (y_opcg1_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg1_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg2_f_c .ne. 0) .and.   &
                 (y_opcg2_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg2_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg3_f_c .ne. 0) .and.   &
                 (y_opcg3_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg3_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg4_f_c .ne. 0) .and.   &
                 (y_opcg4_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg4_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg5_f_c .ne. 0) .and.   &
                 (y_opcg5_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg5_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg6_f_c .ne. 0) .and.   &
                 (y_opcg6_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg6_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg7_f_c .ne. 0) .and.   &
                 (y_opcg7_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg7_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_opcg8_f_c .ne. 0) .and.   &
                 (y_opcg8_f_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for opcg8_f_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_ant1_c .ne. 0) .and.   &
                 (y_ant1_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for ant1_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_ant2_c .ne. 0) .and.   &
                 (y_ant2_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for ant2_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_ant3_c .ne. 0) .and.   &
                 (y_ant3_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for ant3_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_ant4_c .ne. 0) .and.   &
                 (y_ant4_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for ant4_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

	
        else if ((y_terp_Cn3 .ne. 0) .and.   &
                 (y_terp_Cn3 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error terp_Cn3'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_ant1_o .ne. 0) .and.   &
                 (y_ant1_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for ant1_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_ant2_o .ne. 0) .and.   &
                 (y_ant2_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for ant2_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_ant3_o .ne. 0) .and.   &
                 (y_ant3_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for ant3_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_ant4_o .ne. 0) .and.   &
                 (y_ant4_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for ant4_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_biog1_c .ne. 0) .and.   &
                 (y_biog1_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for biog1_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_biog2_c .ne. 0) .and.   &
                 (y_biog2_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for biog2_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_biog3_c .ne. 0) .and.   &
                 (y_biog3_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for biog3_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_biog4_c .ne. 0) .and.   &
                 (y_biog4_c .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for biog4_c'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )


        else if ((y_biog1_o .ne. 0) .and.   &
                 (y_biog1_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for biog1_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_biog2_o .ne. 0) .and.   &
                 (y_biog2_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for biog2_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_biog3_o .ne. 0) .and.   &
                 (y_biog3_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for biog3_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_biog4_o .ne. 0) .and.   &
                 (y_biog4_o .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for biog4_o'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_smpa .ne. 0) .and.   &
                 (y_smpa .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for smpa'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_smpbb .ne. 0) .and.   &
                 (y_smpbb .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for smpbb'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_glysoa_r1 .ne. 0) .and.   &
                 (y_glysoa_r1 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for glysoa_r1'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_glysoa_r2 .ne. 0) .and.   &
                 (y_glysoa_r2 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for glysoa_r2'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_glysoa_sfc .ne. 0) .and.   &
                 (y_glysoa_sfc .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for glysoa_sfc'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
            
        else if ((y_glysoa_nh4 .ne. 0) .and.   &
                 (y_glysoa_nh4 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for glysoa_nh4'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_glysoa_oh .ne. 0) .and.   &
                 (y_glysoa_oh .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for glysoa_oh'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )


        else if ((y_asoaX .ne. 0) .and.   &
                 (y_asoaX .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for y_asoaX'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_asoa1 .ne. 0) .and.   &
                 (y_asoa1 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for y_asoa1'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_asoa2 .ne. 0) .and.   &
                 (y_asoa2 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for y_asoa2'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_asoa3 .ne. 0) .and.   &
                 (y_asoa3 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for y_asoa3'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_asoa4 .ne. 0) .and.   &
                 (y_asoa4 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for y_asoa4'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_bsoaX .ne. 0) .and.   &
                 (y_bsoaX .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for y_bsoaX'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_bsoa1 .ne. 0) .and.   &
                 (y_bsoa1 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for y_bsoa1'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_bsoa2 .ne. 0) .and.   &
                 (y_bsoa2 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for y_bsoa2'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_bsoa3 .ne. 0) .and.   &
                 (y_bsoa3 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for y_bsoa3'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_bsoa4 .ne. 0) .and.   &
                 (y_bsoa4 .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for y_bsoa4'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
	else if ((y_msa .ne. 0) .and.   &
	         (y_msa .ne. nsize_aer(itype))) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for msa'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
	else if ((y_co3 .ne. 0) .and.   &
	         (y_co3 .ne. nsize_aer(itype))) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for co3'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
	else if ((y_nh4 .ne. 0) .and.   &
	         (y_nh4 .ne. nsize_aer(itype))) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for nh4'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
	else if ((y_na .ne. 0) .and.   &
	         (y_na .ne. nsize_aer(itype))) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for na'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
	else if ((y_ca .ne. 0) .and.   &
	         (y_ca .ne. nsize_aer(itype))) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for ca'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
	else if ((y_oin .ne. 0) .and.   &
	         (y_oin .ne. nsize_aer(itype))) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for oin'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
       else if ((y_tanv .ne. 0) .and.   &
                 (y_tanv .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for tanv'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_sopnv .ne. 0) .and.   &
                 (y_sopnv .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for sopnv'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_ternv .ne. 0) .and.   &
                 (y_ternv .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for ternv'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_seqnv .ne. 0) .and.   &
                 (y_seqnv .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for seqnv'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_sianv .ne. 0) .and.   &
                 (y_sianv .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for sianv'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
	else if ((y_oc .ne. 0) .and.   &
	         (y_oc .ne. nsize_aer(itype))) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for oc'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
	else if ((y_bc .ne. 0) .and.   &
	         (y_bc .ne. nsize_aer(itype))) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for bc'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
	else if ((y_hysw .ne. 0) .and.   &
	         (y_hysw .ne. nsize_aer(itype))) then
	    msg = '*** subr init_data_mosaic_ptr - ptr error for hysw'
	    call peg_message( lunerr, msg )
	    write(msg,9350) 'phase, type=', iphase,itype
	    call peg_error_fatal( lunerr, msg )
        else if ((y_iepox .ne. 0) .and.   &
                 (y_iepox .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for iepox'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_iepoxos .ne. 0) .and.   &
                 (y_iepoxos .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for iepoxos'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_tetrol .ne. 0) .and.   &
                 (y_tetrol .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for tetrol'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

        else if ((y_cldiepoxos .ne. 0) .and.   &
                 (y_cldiepoxos .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for cldiepoxos'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )
        else if ((y_cldtetrol .ne. 0) .and.   &
                 (y_cldtetrol .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for cldtetrol'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )


        else if ((y_gly .ne. 0) .and.   &
                 (y_gly .ne. nsize_aer(itype))) then
            msg = '*** subr init_data_mosaic_ptr - ptr error for gly'
            call peg_message( lunerr, msg )
            write(msg,9350) 'phase, type=', iphase,itype
            call peg_error_fatal( lunerr, msg )

	end if

	enddo 
	enddo 




        if (p_h2so4 .ge. p1st) then 
            kh2so4 = p_h2so4
       
       elseif (p_sulf .ge. p1st) then
            kh2so4 = p_sulf
       
 



	end if
	if (p_hno3 .ge. p1st) then
	    khno3 = p_hno3
      endif 
       if (p_pcg1_b_c .ge. p1st) then
            kpcg1_b_c = p_pcg1_b_c
        endif
       if (p_pcg2_b_c .ge. p1st) then
            kpcg2_b_c = p_pcg2_b_c
        endif
       if (p_pcg3_b_c .ge. p1st) then
            kpcg3_b_c = p_pcg3_b_c
        endif
       if (p_pcg4_b_c .ge. p1st) then
            kpcg4_b_c = p_pcg4_b_c
        endif
       if (p_pcg5_b_c .ge. p1st) then
            kpcg5_b_c = p_pcg5_b_c
        endif
       if (p_pcg6_b_c .ge. p1st) then
            kpcg6_b_c = p_pcg6_b_c
        endif
       if (p_pcg7_b_c .ge. p1st) then
            kpcg7_b_c = p_pcg7_b_c
        endif
       if (p_pcg8_b_c .ge. p1st) then
            kpcg8_b_c = p_pcg8_b_c
        endif
       if (p_pcg9_b_c .ge. p1st) then
            kpcg9_b_c = p_pcg9_b_c
        endif
       if (p_pcg1_b_o .ge. p1st) then
            kpcg1_b_o = p_pcg1_b_o
        endif
       if (p_pcg2_b_o .ge. p1st) then
            kpcg2_b_o = p_pcg2_b_o
        endif
       if (p_pcg3_b_o .ge. p1st) then
            kpcg3_b_o = p_pcg3_b_o
        endif
       if (p_pcg4_b_o .ge. p1st) then
            kpcg4_b_o = p_pcg4_b_o
        endif
       if (p_pcg5_b_o .ge. p1st) then
            kpcg5_b_o = p_pcg5_b_o
        endif
       if (p_pcg6_b_o .ge. p1st) then
            kpcg6_b_o = p_pcg6_b_o
        endif
       if (p_pcg7_b_o .ge. p1st) then
            kpcg7_b_o = p_pcg7_b_o
        endif
       if (p_pcg8_b_o .ge. p1st) then
            kpcg8_b_o = p_pcg8_b_o
        endif
       if (p_pcg9_b_o .ge. p1st) then
            kpcg9_b_o = p_pcg9_b_o
        endif
       if (p_opcg1_b_o .ge. p1st) then
            kopcg1_b_o = p_opcg1_b_o
        endif
       if (p_opcg2_b_o .ge. p1st) then
            kopcg2_b_o = p_opcg2_b_o
        endif
       if (p_opcg3_b_o .ge. p1st) then
            kopcg3_b_o = p_opcg3_b_o
        endif
       if (p_opcg4_b_o .ge. p1st) then
            kopcg4_b_o = p_opcg4_b_o
        endif
       if (p_opcg5_b_o .ge. p1st) then
            kopcg5_b_o = p_opcg5_b_o
        endif
       if (p_opcg6_b_o .ge. p1st) then
            kopcg6_b_o = p_opcg6_b_o
        endif
       if (p_opcg7_b_o .ge. p1st) then
            kopcg7_b_o = p_opcg7_b_o
        endif
       if (p_opcg8_b_o .ge. p1st) then
            kopcg8_b_o = p_opcg8_b_o
        endif
       if (p_opcg1_b_c .ge. p1st) then
            kopcg1_b_c = p_opcg1_b_c
        endif
       if (p_opcg2_b_c .ge. p1st) then
            kopcg2_b_c = p_opcg2_b_c
        endif
       if (p_opcg3_b_c .ge. p1st) then
            kopcg3_b_c = p_opcg3_b_c
        endif
       if (p_opcg4_b_c .ge. p1st) then
            kopcg4_b_c = p_opcg4_b_c
        endif
       if (p_opcg5_b_c .ge. p1st) then
            kopcg5_b_c = p_opcg5_b_c
        endif
       if (p_opcg6_b_c .ge. p1st) then
            kopcg6_b_c = p_opcg6_b_c
        endif
       if (p_opcg7_b_c .ge. p1st) then
            kopcg7_b_c = p_opcg7_b_c
        endif
       if (p_opcg8_b_c .ge. p1st) then
            kopcg8_b_c = p_opcg8_b_c
        endif
       if (p_pcg1_f_c .ge. p1st) then
            kpcg1_f_c = p_pcg1_f_c
        endif
       if (p_pcg2_f_c .ge. p1st) then
            kpcg2_f_c = p_pcg2_f_c
        endif
       if (p_pcg3_f_c .ge. p1st) then
            kpcg3_f_c = p_pcg3_f_c
        endif
       if (p_pcg4_f_c .ge. p1st) then
            kpcg4_f_c = p_pcg4_f_c
        endif
       if (p_pcg5_f_c .ge. p1st) then
            kpcg5_f_c = p_pcg5_f_c
        endif
       if (p_pcg6_f_c .ge. p1st) then
            kpcg6_f_c = p_pcg6_f_c
        endif
       if (p_pcg7_f_c .ge. p1st) then
            kpcg7_f_c = p_pcg7_f_c
        endif
       if (p_pcg8_f_c .ge. p1st) then
            kpcg8_f_c = p_pcg8_f_c
        endif
       if (p_pcg9_f_c .ge. p1st) then
            kpcg9_f_c = p_pcg9_f_c
        endif
       if (p_pcg1_f_o .ge. p1st) then
            kpcg1_f_o = p_pcg1_f_o
        endif
       if (p_pcg2_f_o .ge. p1st) then
            kpcg2_f_o = p_pcg2_f_o
        endif
       if (p_pcg3_f_o .ge. p1st) then
            kpcg3_f_o = p_pcg3_f_o
        endif
       if (p_pcg4_f_o .ge. p1st) then
            kpcg4_f_o = p_pcg4_f_o
        endif
       if (p_pcg5_f_o .ge. p1st) then
            kpcg5_f_o = p_pcg5_f_o
        endif
       if (p_pcg6_f_o .ge. p1st) then
            kpcg6_f_o = p_pcg6_f_o
        endif
       if (p_pcg7_f_o .ge. p1st) then
            kpcg7_f_o = p_pcg7_f_o
        endif
       if (p_pcg8_f_o .ge. p1st) then
            kpcg8_f_o = p_pcg8_f_o
        endif
       if (p_pcg9_f_o .ge. p1st) then
            kpcg9_f_o = p_pcg9_f_o
        endif
       if (p_opcg1_f_o .ge. p1st) then
            kopcg1_f_o = p_opcg1_f_o
        endif
       if (p_opcg2_f_o .ge. p1st) then
            kopcg2_f_o = p_opcg2_f_o
        endif
       if (p_opcg3_f_o .ge. p1st) then
            kopcg3_f_o = p_opcg3_f_o
        endif
       if (p_opcg4_f_o .ge. p1st) then
            kopcg4_f_o = p_opcg4_f_o
        endif
       if (p_opcg5_f_o .ge. p1st) then
            kopcg5_f_o = p_opcg5_f_o
        endif
       if (p_opcg6_f_o .ge. p1st) then
            kopcg6_f_o = p_opcg6_f_o
        endif
       if (p_opcg7_f_o .ge. p1st) then
            kopcg7_f_o = p_opcg7_f_o
        endif
       if (p_opcg8_f_o .ge. p1st) then
            kopcg8_f_o = p_opcg8_f_o
        endif
       if (p_opcg1_f_c .ge. p1st) then
            kopcg1_f_c = p_opcg1_f_c
        endif
       if (p_opcg2_f_c .ge. p1st) then
            kopcg2_f_c = p_opcg2_f_c
        endif
       if (p_opcg3_f_c .ge. p1st) then
            kopcg3_f_c = p_opcg3_f_c
        endif
       if (p_opcg4_f_c .ge. p1st) then
            kopcg4_f_c = p_opcg4_f_c
        endif
       if (p_opcg5_f_c .ge. p1st) then
            kopcg5_f_c = p_opcg5_f_c
        endif
       if (p_opcg6_f_c .ge. p1st) then
            kopcg6_f_c = p_opcg6_f_c
        endif
       if (p_opcg7_f_c .ge. p1st) then
            kopcg7_f_c = p_opcg7_f_c
        endif
       if (p_opcg8_f_c .ge. p1st) then
            kopcg8_f_c = p_opcg8_f_c
        endif

       if (p_smpa .ge. p1st) then
            ksmpa = p_smpa
        endif
       if (p_smpbb .ge. p1st) then
            ksmpbb = p_smpbb
        endif

       if (p_gly .ge. p1st) then
           kgly = p_gly
       endif

       if (p_cvasoaX .ge. p1st) then
           kasoaX = p_cvasoaX
       endif

       if (p_cvasoa1 .ge. p1st) then
           kasoa1 = p_cvasoa1
       endif

       if (p_cvasoa2 .ge. p1st) then
           kasoa2 = p_cvasoa2
       endif

       if (p_cvasoa3 .ge. p1st) then
           kasoa3 = p_cvasoa3
       endif

       if (p_cvasoa4 .ge. p1st) then
           kasoa4 = p_cvasoa4
       endif

       if (p_cvbsoaX .ge. p1st) then
           kbsoaX = p_cvbsoaX
       endif

       if (p_cvbsoa1 .ge. p1st) then
           kbsoa1 = p_cvbsoa1
       endif

       if (p_cvbsoa2 .ge. p1st) then
           kbsoa2 = p_cvbsoa2
       endif

       if (p_cvbsoa3 .ge. p1st) then
           kbsoa3 = p_cvbsoa3
       endif

       if (p_cvbsoa4 .ge. p1st) then
           kbsoa4 = p_cvbsoa4
       endif
       if (p_ant1_c .ge. p1st) then
            kant1_c = p_ant1_c
        endif

       if (p_ant2_c .ge. p1st) then
            kant2_c = p_ant2_c
        endif
       if (p_ant3_c .ge. p1st) then
            kant3_c = p_ant3_c
        endif
       if (p_ant4_c .ge. p1st) then
            kant4_c = p_ant4_c
        endif

	
       if (p_TERP_CN3_g .ge. p1st) then
            kterp_g_Cn3 = p_TERP_CN3_g
        endif

       if (p_ant1_o .ge. p1st) then
            kant1_o = p_ant1_o
        endif

       if (p_ant2_o .ge. p1st) then
            kant2_o = p_ant2_o
        endif
       if (p_ant3_o .ge. p1st) then
            kant3_o = p_ant3_o
        endif
       if (p_ant4_o .ge. p1st) then
            kant4_o = p_ant4_o
        endif

       if (p_biog1_c .ge. p1st) then
            kbiog1_c = p_biog1_c
        endif

       if (p_biog2_c .ge. p1st) then
            kbiog2_c = p_biog2_c
        endif
       if (p_biog3_c .ge. p1st) then
            kbiog3_c = p_biog3_c
        endif
       if (p_biog4_c .ge. p1st) then
            kbiog4_c = p_biog4_c
        endif
       if (p_biog1_o .ge. p1st) then
            kbiog1_o = p_biog1_o
        endif

       if (p_biog2_o .ge. p1st) then
            kbiog2_o = p_biog2_o
        endif
       if (p_biog3_o .ge. p1st) then
            kbiog3_o = p_biog3_o
        endif
       if (p_biog4_o .ge. p1st) then
            kbiog4_o = p_biog4_o
        endif

       if (p_iepox .ge. p1st) then
            kiepox = p_iepox
        endif

       if (p_tetrol .ge. p1st) then
            ktetrol = p_tetrol
        endif

       if (p_cldtetrol .ge. p1st) then
            kcldtetrol = p_cldtetrol
        endif


        if (p_gly .ge. p1st) then
            kgly = p_gly
        endif





	if (p_hcl .ge. p1st) then
	    khcl = p_hcl



	end if
	if (p_nh3 .ge. p1st) then
	    knh3 = p_nh3



	end if
	if (p_n2o5 .ge. p1st) then
	    kn2o5 = p_n2o5



	end if
	if (p_clno2 .ge. p1st) then
	    kclno2 = p_clno2



	end if
	if (p_o3 .ge. p1st) then
	    ko3 = p_o3



	end if




	if (p_so2    .ge. p1st) kso2    = p_so2
	if (p_h2o2   .ge. p1st) kh2o2   = p_h2o2
	if (p_hcho   .ge. p1st) khcho   = p_hcho
	if (p_ho     .ge. p1st) koh     = p_ho
	if (p_ho2    .ge. p1st) kho2    = p_ho2
	if (p_no3    .ge. p1st) kno3    = p_no3
	if (p_no     .ge. p1st) kno     = p_no
	if (p_no2    .ge. p1st) kno2    = p_no2
	if (p_hono   .ge. p1st) khono   = p_hono
	if (p_pan    .ge. p1st) kpan    = p_pan




	is_aerosol(:) = .false.
	ltot = 0
	ltot = max( ltot, kh2so4 )
	ltot = max( ltot, khno3 )
	ltot = max( ltot, khcl )
	ltot = max( ltot, knh3 )
	ltot = max( ltot, kn2o5 )
	ltot = max( ltot, kclno2 )
	ltot = max( ltot, ko3 )
	ltot = max( ltot, kso2    )
	ltot = max( ltot, kh2o2   )
	ltot = max( ltot, khcho   )
	ltot = max( ltot, koh     )
	ltot = max( ltot, kho2    )
	ltot = max( ltot, kno3    )
	ltot = max( ltot, kno     )
	ltot = max( ltot, kno2    )
	ltot = max( ltot, khono   )
	ltot = max( ltot, kpan    )
        ltot = max( ltot, kpcg1_b_c )
        ltot = max( ltot, kpcg2_b_c )
        ltot = max( ltot, kpcg3_b_c )
        ltot = max( ltot, kpcg4_b_c )
        ltot = max( ltot, kpcg5_b_c )
        ltot = max( ltot, kpcg6_b_c )
        ltot = max( ltot, kpcg7_b_c )
        ltot = max( ltot, kpcg8_b_c )
        ltot = max( ltot, kpcg9_b_c )
        ltot = max( ltot, kpcg1_b_o )
        ltot = max( ltot, kpcg2_b_o )
        ltot = max( ltot, kpcg3_b_o )
        ltot = max( ltot, kpcg4_b_o )
        ltot = max( ltot, kpcg5_b_o )
        ltot = max( ltot, kpcg6_b_o )
        ltot = max( ltot, kpcg7_b_o )
        ltot = max( ltot, kpcg8_b_o )
        ltot = max( ltot, kpcg9_b_o )
        ltot = max( ltot, kopcg1_b_c )
        ltot = max( ltot, kopcg2_b_c )
        ltot = max( ltot, kopcg3_b_c )
        ltot = max( ltot, kopcg4_b_c )
        ltot = max( ltot, kopcg5_b_c )
        ltot = max( ltot, kopcg6_b_c )
        ltot = max( ltot, kopcg7_b_c )
        ltot = max( ltot, kopcg8_b_c )
        ltot = max( ltot, kopcg1_b_o )
        ltot = max( ltot, kopcg2_b_o )
        ltot = max( ltot, kopcg3_b_o )
        ltot = max( ltot, kopcg4_b_o )
        ltot = max( ltot, kopcg5_b_o )
        ltot = max( ltot, kopcg6_b_o )
        ltot = max( ltot, kopcg7_b_o )
        ltot = max( ltot, kopcg8_b_o )
        ltot = max( ltot, kpcg1_f_c )
        ltot = max( ltot, kpcg2_f_c )
        ltot = max( ltot, kpcg3_f_c )
        ltot = max( ltot, kpcg4_f_c )
        ltot = max( ltot, kpcg5_f_c )
        ltot = max( ltot, kpcg6_f_c )
        ltot = max( ltot, kpcg7_f_c )
        ltot = max( ltot, kpcg8_f_c )
        ltot = max( ltot, kpcg9_f_c )
        ltot = max( ltot, kpcg1_f_o )
        ltot = max( ltot, kpcg2_f_o )
        ltot = max( ltot, kpcg3_f_o )
        ltot = max( ltot, kpcg4_f_o )
        ltot = max( ltot, kpcg5_f_o )
        ltot = max( ltot, kpcg6_f_o )
        ltot = max( ltot, kpcg7_f_o )
        ltot = max( ltot, kpcg8_f_o )
        ltot = max( ltot, kpcg9_f_o )
        ltot = max( ltot, kopcg1_f_c )
        ltot = max( ltot, kopcg2_f_c )
        ltot = max( ltot, kopcg3_f_c )
        ltot = max( ltot, kopcg4_f_c )
        ltot = max( ltot, kopcg5_f_c )
        ltot = max( ltot, kopcg6_f_c )
        ltot = max( ltot, kopcg7_f_c )
        ltot = max( ltot, kopcg8_f_c )
        ltot = max( ltot, kopcg1_f_o )
        ltot = max( ltot, kopcg2_f_o )
        ltot = max( ltot, kopcg3_f_o )
        ltot = max( ltot, kopcg4_f_o )
        ltot = max( ltot, kopcg5_f_o )
        ltot = max( ltot, kopcg6_f_o )
        ltot = max( ltot, kopcg7_f_o )
        ltot = max( ltot, kopcg8_f_o )
        ltot = max( ltot, ksmpa )
        ltot = max( ltot, ksmpbb )
        ltot = max( ltot, kgly    )
        ltot = max( ltot, kant1_c )
        ltot = max( ltot, kant2_c )
        ltot = max( ltot, kant3_c )
	
	ltot = max( ltot, kterp_g_Cn3 )
        ltot = max( ltot, kant4_c )
        ltot = max( ltot, kant1_o )
        ltot = max( ltot, kant2_o )
        ltot = max( ltot, kant3_o )
        ltot = max( ltot, kant4_o )
        ltot = max( ltot, kbiog1_c )
        ltot = max( ltot, kbiog2_c )
        ltot = max( ltot, kbiog3_c )
        ltot = max( ltot, kbiog4_c )
        ltot = max( ltot, kbiog1_o )
        ltot = max( ltot, kbiog2_o )
        ltot = max( ltot, kbiog3_o )
        ltot = max( ltot, kbiog4_o )
        ltot = max( ltot, kiepox )
        ltot = max( ltot, ktetrol )
        ltot = max( ltot, kcldtetrol )
        ltot = max( ltot, kgly )
        ltot = max( ltot, kasoaX )
        ltot = max( ltot, kasoa1 )
        ltot = max( ltot, kasoa2 )
        ltot = max( ltot, kasoa3 )
        ltot = max( ltot, kasoa4 )
        ltot = max( ltot, kbsoaX )
        ltot = max( ltot, kbsoa1 )
        ltot = max( ltot, kbsoa2 )
        ltot = max( ltot, kbsoa3 )
        ltot = max( ltot, kbsoa4 )

	do iphase=1,nphase_aer
	    do itype=1,ntype_aer
		do n = 1, nsize_aer(itype)
		    do ll = 1, ncomp_plustracer_aer(itype)
		       ltot = max( ltot, massptr_aer(ll,n,itype,iphase) )
		       is_aerosol(massptr_aer(ll,n,itype,iphase))=.true.
		    end do
		    ltot = max( ltot, hyswptr_aer(n,itype) )
		    ltot = max( ltot, waterptr_aer(n,itype) )
		    ltot = max( ltot, numptr_aer(n,itype,iphase) )
		    l = hyswptr_aer(n,itype)
		    if (l .ge. p1st) is_aerosol(l)=.true.
		    l = waterptr_aer(n,itype)
		    if (l .ge. p1st) is_aerosol(l)=.true.
		    l = numptr_aer(n,itype,iphase)
		    if (l .ge. p1st) is_aerosol(l)=.true.
		end do
	    end do
	end do

	kh2o = ltot + 1
	ktemp = ltot + 2
	ltot2 = ktemp

	write( msg, '(a,4(1x,i4))' ) 'ltot, ltot2, lmaxd, l2maxd =',   &
		ltot, ltot2, lmaxd, l2maxd
	call peg_message( lunout, msg )
	if ((ltot .gt. lmaxd) .or. (ltot2 .gt. l2maxd)) then
	    msg = '*** subr init_data_mosaic_ptr - ltot/ltot2 too big'
	    call peg_error_fatal( lunerr, msg )
	end if

        if (p_h2so4   .ge. p1st)then 
           name(kh2so4 ) = 'h2so4'
        elseif (p_sulf   .ge. p1st) then
             name(kh2so4 ) = 'h2so4'
        endif
	if (p_hno3   .ge. p1st) name(khno3  ) = 'hno3'
        if (p_pcg1_b_c   .ge. p1st) name(kpcg1_b_c  ) = 'pcg1_b_c'
        if (p_pcg2_b_c   .ge. p1st) name(kpcg2_b_c  ) = 'pcg2_b_c'
        if (p_pcg3_b_c   .ge. p1st) name(kpcg3_b_c  ) = 'pcg3_b_c'
        if (p_pcg4_b_c   .ge. p1st) name(kpcg4_b_c  ) = 'pcg4_b_c'
        if (p_pcg5_b_c   .ge. p1st) name(kpcg5_b_c  ) = 'pcg5_b_c'
        if (p_pcg6_b_c   .ge. p1st) name(kpcg6_b_c  ) = 'pcg6_b_c'
        if (p_pcg7_b_c   .ge. p1st) name(kpcg7_b_c  ) = 'pcg7_b_c'
        if (p_pcg8_b_c   .ge. p1st) name(kpcg8_b_c  ) = 'pcg8_b_c'
        if (p_pcg9_b_c   .ge. p1st) name(kpcg9_b_c  ) = 'pcg9_b_c'
        if (p_pcg1_b_o   .ge. p1st) name(kpcg1_b_o  ) = 'pcg1_b_o'
        if (p_pcg2_b_o   .ge. p1st) name(kpcg2_b_o  ) = 'pcg2_b_o'
        if (p_pcg3_b_o   .ge. p1st) name(kpcg3_b_o  ) = 'pcg3_b_o'
        if (p_pcg4_b_o   .ge. p1st) name(kpcg4_b_o  ) = 'pcg4_b_o'
        if (p_pcg5_b_o   .ge. p1st) name(kpcg5_b_o  ) = 'pcg5_b_o'
        if (p_pcg6_b_o   .ge. p1st) name(kpcg6_b_o  ) = 'pcg6_b_o'
        if (p_pcg7_b_o   .ge. p1st) name(kpcg7_b_o  ) = 'pcg7_b_o'
        if (p_pcg8_b_o   .ge. p1st) name(kpcg8_b_o  ) = 'pcg8_b_o'
        if (p_pcg9_b_o   .ge. p1st) name(kpcg9_b_o  ) = 'pcg9_b_o'
        if (p_opcg1_b_c   .ge. p1st) name(kopcg1_b_c  ) = 'opcg1_b_c'
        if (p_opcg2_b_c   .ge. p1st) name(kopcg2_b_c  ) = 'opcg2_b_c'
        if (p_opcg3_b_c   .ge. p1st) name(kopcg3_b_c  ) = 'opcg3_b_c'
        if (p_opcg4_b_c   .ge. p1st) name(kopcg4_b_c  ) = 'opcg4_b_c'
        if (p_opcg5_b_c   .ge. p1st) name(kopcg5_b_c  ) = 'opcg5_b_c'
        if (p_opcg6_b_c   .ge. p1st) name(kopcg6_b_c  ) = 'opcg6_b_c'
        if (p_opcg7_b_c   .ge. p1st) name(kopcg7_b_c  ) = 'opcg7_b_c'
        if (p_opcg8_b_c   .ge. p1st) name(kopcg8_b_c  ) = 'opcg8_b_c'
        if (p_opcg1_b_o   .ge. p1st) name(kopcg1_b_o  ) = 'opcg1_b_o'
        if (p_opcg2_b_o   .ge. p1st) name(kopcg2_b_o  ) = 'opcg2_b_o'
        if (p_opcg3_b_o   .ge. p1st) name(kopcg3_b_o  ) = 'opcg3_b_o'
        if (p_opcg4_b_o   .ge. p1st) name(kopcg4_b_o  ) = 'opcg4_b_o'
        if (p_opcg5_b_o   .ge. p1st) name(kopcg5_b_o  ) = 'opcg5_b_o'
        if (p_opcg6_b_o   .ge. p1st) name(kopcg6_b_o  ) = 'opcg6_b_o'
        if (p_opcg7_b_o   .ge. p1st) name(kopcg7_b_o  ) = 'opcg7_b_o'
        if (p_opcg8_b_o   .ge. p1st) name(kopcg8_b_o  ) = 'opcg8_b_o'
        if (p_pcg1_f_c   .ge. p1st) name(kpcg1_f_c  ) = 'pcg1_f_c'
        if (p_pcg2_f_c   .ge. p1st) name(kpcg2_f_c  ) = 'pcg2_f_c'
        if (p_pcg3_f_c   .ge. p1st) name(kpcg3_f_c  ) = 'pcg3_f_c'
        if (p_pcg4_f_c   .ge. p1st) name(kpcg4_f_c  ) = 'pcg4_f_c'
        if (p_pcg5_f_c   .ge. p1st) name(kpcg5_f_c  ) = 'pcg5_f_c'
        if (p_pcg6_f_c   .ge. p1st) name(kpcg6_f_c  ) = 'pcg6_f_c'
        if (p_pcg7_f_c   .ge. p1st) name(kpcg7_f_c  ) = 'pcg7_f_c'
        if (p_pcg8_f_c   .ge. p1st) name(kpcg8_f_c  ) = 'pcg8_f_c'
        if (p_pcg9_f_c   .ge. p1st) name(kpcg9_f_c  ) = 'pcg9_f_c'
        if (p_pcg1_f_o   .ge. p1st) name(kpcg1_f_o  ) = 'pcg1_f_o'
        if (p_pcg2_f_o   .ge. p1st) name(kpcg2_f_o  ) = 'pcg2_f_o'
        if (p_pcg3_f_o   .ge. p1st) name(kpcg3_f_o  ) = 'pcg3_f_o'
        if (p_pcg4_f_o   .ge. p1st) name(kpcg4_f_o  ) = 'pcg4_f_o'
        if (p_pcg5_f_o   .ge. p1st) name(kpcg5_f_o  ) = 'pcg5_f_o'
        if (p_pcg6_f_o   .ge. p1st) name(kpcg6_f_o  ) = 'pcg6_f_o'
        if (p_pcg7_f_o   .ge. p1st) name(kpcg7_f_o  ) = 'pcg7_f_o'
        if (p_pcg8_f_o   .ge. p1st) name(kpcg8_f_o  ) = 'pcg8_f_o'
        if (p_pcg9_f_o   .ge. p1st) name(kpcg9_f_o  ) = 'pcg9_f_o'
        if (p_opcg1_f_c   .ge. p1st) name(kopcg1_f_c  ) = 'opcg1_f_c'
        if (p_opcg2_f_c   .ge. p1st) name(kopcg2_f_c  ) = 'opcg2_f_c'
        if (p_opcg3_f_c   .ge. p1st) name(kopcg3_f_c  ) = 'opcg3_f_c'
        if (p_opcg4_f_c   .ge. p1st) name(kopcg4_f_c  ) = 'opcg4_f_c'
        if (p_opcg5_f_c   .ge. p1st) name(kopcg5_f_c  ) = 'opcg5_f_c'
        if (p_opcg6_f_c   .ge. p1st) name(kopcg6_f_c  ) = 'opcg6_f_c'
        if (p_opcg7_f_c   .ge. p1st) name(kopcg7_f_c  ) = 'opcg7_f_c'
        if (p_opcg8_f_c   .ge. p1st) name(kopcg8_f_c  ) = 'opcg8_f_c'
        if (p_opcg1_f_o   .ge. p1st) name(kopcg1_f_o  ) = 'opcg1_f_o'
        if (p_opcg2_f_o   .ge. p1st) name(kopcg2_f_o  ) = 'opcg2_f_o'
        if (p_opcg3_f_o   .ge. p1st) name(kopcg3_f_o  ) = 'opcg3_f_o'
        if (p_opcg4_f_o   .ge. p1st) name(kopcg4_f_o  ) = 'opcg4_f_o'
        if (p_opcg5_f_o   .ge. p1st) name(kopcg5_f_o  ) = 'opcg5_f_o'
        if (p_opcg6_f_o   .ge. p1st) name(kopcg6_f_o  ) = 'opcg6_f_o'
        if (p_opcg7_f_o   .ge. p1st) name(kopcg7_f_o  ) = 'opcg7_f_o'
        if (p_opcg8_f_o   .ge. p1st) name(kopcg8_f_o  ) = 'opcg8_f_o'
        if (p_smpa   .ge. p1st) name(ksmpa  ) = 'smpa'
        if (p_smpbb   .ge. p1st) name(ksmpbb  ) = 'smpbb'
        if (p_iepox   .ge. p1st) name(kiepox  ) = 'iepox'
        if (p_tetrol   .ge. p1st) name(ktetrol  ) = 'tetrol'
       if (p_cldtetrol   .ge. p1st) name(kcldtetrol  ) = 'cldtetrol'
        if (p_gly    .ge. p1st) name(kgly   ) = 'gly'
        if (p_ant1_c   .ge. p1st) name(kant1_c  ) = 'ant1_c'
        if (p_ant2_c   .ge. p1st) name(kant2_c  ) = 'ant2_c'
        if (p_ant3_c   .ge. p1st) name(kant3_c  ) = 'ant3_c'
	
	if (p_TERP_CN3_g   .ge. p1st) name(kterp_g_Cn3  ) = 'terp_Cn3'
        if (p_ant4_c   .ge. p1st) name(kant4_c  ) = 'ant4_c'
        if (p_ant1_o   .ge. p1st) name(kant1_o  ) = 'ant1_o'
        if (p_ant2_o   .ge. p1st) name(kant2_o  ) = 'ant2_o'
        if (p_ant3_o   .ge. p1st) name(kant3_o  ) = 'ant3_o'
        if (p_ant4_o   .ge. p1st) name(kant4_o  ) = 'ant4_o'
        if (p_biog1_c   .ge. p1st) name(kbiog1_c  ) = 'biog1_c'
        if (p_biog2_c   .ge. p1st) name(kbiog2_c  ) = 'biog2_c'
        if (p_biog3_c   .ge. p1st) name(kbiog3_c  ) = 'biog3_c'
        if (p_biog4_c   .ge. p1st) name(kbiog4_c  ) = 'biog4_c'
        if (p_biog1_o   .ge. p1st) name(kbiog1_o  ) = 'biog1_o'
        if (p_biog2_o   .ge. p1st) name(kbiog2_o  ) = 'biog2_o'
        if (p_biog3_o   .ge. p1st) name(kbiog3_o  ) = 'biog3_o'
        if (p_biog4_o   .ge. p1st) name(kbiog4_o  ) = 'biog4_o'
        if (p_cvasoaX   .ge. p1st) name(kasoaX  ) = 'cvasoaX'
        if (p_cvasoa1   .ge. p1st) name(kasoa1  ) = 'cvasoa1'
        if (p_cvasoa2   .ge. p1st) name(kasoa2  ) = 'cvasoa2'
        if (p_cvasoa3   .ge. p1st) name(kasoa3  ) = 'cvasoa3'
        if (p_cvasoa4   .ge. p1st) name(kasoa4  ) = 'cvasoa4'
        if (p_cvbsoaX   .ge. p1st) name(kbsoaX  ) = 'cvbsoaX'
        if (p_cvbsoa1   .ge. p1st) name(kbsoa1  ) = 'cvbsoa1'
        if (p_cvbsoa2   .ge. p1st) name(kbsoa2  ) = 'cvbsoa2'
        if (p_cvbsoa3   .ge. p1st) name(kbsoa3  ) = 'cvbsoa3'
        if (p_cvbsoa4   .ge. p1st) name(kbsoa4  ) = 'cvbsoa4'
	if (p_hcl    .ge. p1st) name(khcl   ) = 'hcl'
	if (p_nh3    .ge. p1st) name(knh3   ) = 'nh3'
	if (p_n2o5   .ge. p1st) name(kn2o5  ) = 'n2o5'
	if (p_clno2  .ge. p1st) name(kclno2 ) = 'clno2'
	if (p_o3     .ge. p1st) name(ko3    ) = 'o3'
	if (p_so2    .ge. p1st) name(kso2   ) = 'so2'
	if (p_h2o2   .ge. p1st) name(kh2o2  ) = 'h2o2'
	if (p_hcho   .ge. p1st) name(khcho  ) = 'hcho'
	if (p_ho     .ge. p1st) name(koh    ) = 'oh'
	if (p_ho2    .ge. p1st) name(kho2   ) = 'ho2'
	if (p_no3    .ge. p1st) name(kno3   ) = 'no3'
	if (p_no     .ge. p1st) name(kno    ) = 'no'
	if (p_no2    .ge. p1st) name(kno2   ) = 'no2'
	if (p_hono   .ge. p1st) name(khono  ) = 'hono'
	if (p_pan    .ge. p1st) name(kpan   ) = 'pan'
	name(ktemp)  = 'temp'
	name(kh2o)   = 'h2o'

        call initwet(   &
	    ntype_aer, nsize_aer, ncomp_aer,   &
	    massptr_aer, dens_aer, numptr_aer,           &
	    maxd_acomp, maxd_asize,maxd_atype, maxd_aphase,   &
	    dcen_sect, sigmag_aer, &
	    waterptr_aer, dens_water_aer, &
	    scavimptblvol, scavimptblnum, nimptblgrow_mind,   &
	    nimptblgrow_maxd, dlndg_nimptblgrow)

	return
	end subroutine init_data_mosaic_ptr










	subroutine aerchem_debug_dump(   &
      		iflag, iclm, jclm, dtchem )

	use module_data_mosaic_asect
	use module_data_mosaic_other
	implicit none





	integer iflag, iclm, jclm
	real dtchem


	integer ientryno
	save ientryno
	integer iphase, itype, k, l, m, n

	real dtchem_sv1
	save dtchem_sv1
	real rsub_sv1(l2maxd,kmaxd,nsubareamaxd)

	data ientryno / -13579 /






	if (ientryno .ne. -13579) goto 1000

	ientryno = +1

95010	format( a )
95020	format( 8( 1x, i8 ) )
95030	format( 4( 1pe18.10 ) )

	print 95010, 'aerchem_debug_dump start -- first time'
	print 95020, ltot, ltot2, itot, jtot, ktot
	print 95010, (name(l), l=1,ltot2)

	print 95020, maerocoag, maerchem, maeroptical
	print 95020, msectional, maerosolincw
	do iphase = 1, nphase_aer
	do itype=1,ntype_aer
	print 95020, iphase, itype, nsize_aer(itype),   &
     		ncomp_plustracer_aer(itype)

	do n = 1, ncomp_plustracer_aer(itype)
	    print 95010,   &
      		name_aer(n,itype)
	    print 95030,   &
      		dens_aer(n,itype),     mw_aer(n,itype)
	end do

	do n = 1, nsize_aer(itype)
	    print 95020,   &
      		ncomp_plustracer_aer(n),       ncomp_aer(n),   &
      		waterptr_aer(n,itype),   numptr_aer(n,itype,iphase),    &
      		mprognum_aer(n,itype,iphase)
	    print 95020,   &
      		(mastercompptr_aer(l,itype), massptr_aer(l,n,itype,iphase),    &
      		l=1,ncomp_plustracer_aer(itype))
	    print 95030,   &
      		volumcen_sect(n,itype),     volumlo_sect(n,itype),   &
      		volumhi_sect(n,itype),      dcen_sect(n,itype),   &
      		dlo_sect(n,itype),          dhi_sect(n,itype)
	    print 95020,   &
      		lptr_so4_aer(n,itype,iphase),  lptr_msa_aer(n,itype,iphase),  &
      		lptr_no3_aer(n,itype,iphase),  lptr_cl_aer(n,itype,iphase),   &
      		lptr_co3_aer(n,itype,iphase),  lptr_nh4_aer(n,itype,iphase),  &
      		lptr_na_aer(n,itype,iphase),   lptr_ca_aer(n,itype,iphase),   &
      		lptr_oin_aer(n,itype,iphase),  lptr_oc_aer(n,itype,iphase),   &
      		lptr_bc_aer(n,itype,iphase),   hyswptr_aer(n,itype), &
                lptr_tanv_aer(n,itype,iphase),  lptr_sopnv_aer(n,itype,iphase),   &
                lptr_ternv_aer(n,itype,iphase),  lptr_seqnv_aer(n,itype,iphase),   &
                lptr_sianv_aer(n,itype,iphase),   &
                lptr_pcg1_b_c_aer(n,itype,iphase),  lptr_pcg2_b_c_aer(n,itype,iphase),&
                lptr_pcg3_b_c_aer(n,itype,iphase),  lptr_pcg4_b_c_aer(n,itype,iphase),&
                lptr_pcg5_b_c_aer(n,itype,iphase),  lptr_pcg6_b_c_aer(n,itype,iphase),&
                lptr_pcg7_b_c_aer(n,itype,iphase),  lptr_pcg8_b_c_aer(n,itype,iphase),&
                lptr_pcg9_b_c_aer(n,itype,iphase),  lptr_pcg1_b_o_aer(n,itype,iphase),&
                lptr_pcg2_b_o_aer(n,itype,iphase), lptr_pcg3_b_o_aer(n,itype,iphase), &
                lptr_pcg4_b_o_aer(n,itype,iphase), lptr_pcg5_b_o_aer(n,itype,iphase), &
                lptr_pcg6_b_o_aer(n,itype,iphase), lptr_pcg7_b_o_aer(n,itype,iphase), &
                lptr_pcg8_b_o_aer(n,itype,iphase), lptr_pcg9_b_o_aer(n,itype,iphase), &
                lptr_opcg1_b_c_aer(n,itype,iphase),&
                lptr_opcg2_b_c_aer(n,itype,iphase),  lptr_opcg3_b_c_aer(n,itype,iphase),&
                lptr_opcg4_b_c_aer(n,itype,iphase),  lptr_opcg5_b_c_aer(n,itype,iphase),&
                lptr_opcg6_b_c_aer(n,itype,iphase),  lptr_opcg7_b_c_aer(n,itype,iphase),&
                lptr_opcg8_b_c_aer(n,itype,iphase),  lptr_opcg1_b_o_aer(n,itype,iphase),&
                lptr_opcg2_b_o_aer(n,itype,iphase),  lptr_opcg3_b_o_aer(n,itype,iphase),&
                lptr_opcg4_b_o_aer(n,itype,iphase),  lptr_opcg5_b_o_aer(n,itype,iphase),&
                lptr_opcg6_b_o_aer(n,itype,iphase),  lptr_opcg7_b_o_aer(n,itype,iphase),&
                lptr_opcg8_b_o_aer(n,itype,iphase),  &
                lptr_pcg1_f_c_aer(n,itype,iphase),  lptr_pcg2_f_c_aer(n,itype,iphase),&
                lptr_pcg3_f_c_aer(n,itype,iphase),  lptr_pcg4_f_c_aer(n,itype,iphase),&
                lptr_pcg5_f_c_aer(n,itype,iphase),  lptr_pcg6_f_c_aer(n,itype,iphase),&
                lptr_pcg7_f_c_aer(n,itype,iphase),  lptr_pcg8_f_c_aer(n,itype,iphase),&
                lptr_pcg9_f_c_aer(n,itype,iphase),  lptr_pcg1_f_o_aer(n,itype,iphase),&
                lptr_pcg2_f_o_aer(n,itype,iphase), lptr_pcg3_f_o_aer(n,itype,iphase), &
                lptr_pcg4_f_o_aer(n,itype,iphase), lptr_pcg5_f_o_aer(n,itype,iphase), &
                lptr_pcg6_f_o_aer(n,itype,iphase), lptr_pcg7_f_o_aer(n,itype,iphase), &
                lptr_pcg8_f_o_aer(n,itype,iphase), lptr_pcg9_f_o_aer(n,itype,iphase), &
                lptr_opcg1_f_c_aer(n,itype,iphase),&
                lptr_opcg2_f_c_aer(n,itype,iphase),  lptr_opcg3_f_c_aer(n,itype,iphase),&
                lptr_opcg4_f_c_aer(n,itype,iphase),  lptr_opcg5_f_c_aer(n,itype,iphase),&
                lptr_opcg6_f_c_aer(n,itype,iphase),  lptr_opcg7_f_c_aer(n,itype,iphase),&
                lptr_opcg8_f_c_aer(n,itype,iphase),  lptr_opcg1_f_o_aer(n,itype,iphase),&
                lptr_opcg2_f_o_aer(n,itype,iphase),  lptr_opcg3_f_o_aer(n,itype,iphase),&
                lptr_opcg4_f_o_aer(n,itype,iphase),  lptr_opcg5_f_o_aer(n,itype,iphase),&
                lptr_opcg6_f_o_aer(n,itype,iphase),  lptr_opcg7_f_o_aer(n,itype,iphase),&
                lptr_opcg8_f_o_aer(n,itype,iphase),                                     &
                lptr_smpa_aer(n,itype,iphase),lptr_smpbb_aer(n,itype,iphase),        &
                lptr_iepox_aer(n,itype,iphase),lptr_iepoxos_aer(n,itype,iphase),        &
                lptr_tetrol_aer(n,itype,iphase),lptr_gly_aer(n,itype,iphase),        &
                lptr_cldiepoxos_aer(n,itype,iphase),lptr_cldtetrol_aer(n,itype,iphase), &
                lptr_glysoa_r1_aer(n,itype,iphase), &
                lptr_glysoa_r2_aer(n,itype,iphase), &
                lptr_glysoa_sfc_aer(n,itype,iphase), &
                lptr_glysoa_nh4_aer(n,itype,iphase), &
                lptr_glysoa_oh_aer(n,itype,iphase), &
                lptr_ant1_c_aer(n,itype,iphase),lptr_ant2_c_aer(n,itype,iphase),        &
                lptr_ant3_c_aer(n,itype,iphase),lptr_ant4_c_aer(n,itype,iphase),        &
                lptr_ant1_o_aer(n,itype,iphase),lptr_ant2_o_aer(n,itype,iphase),        &
                lptr_ant3_o_aer(n,itype,iphase),lptr_ant4_o_aer(n,itype,iphase),        &
                lptr_biog1_c_aer(n,itype,iphase),lptr_biog2_c_aer(n,itype,iphase),        &
                lptr_biog3_c_aer(n,itype,iphase),lptr_biog4_c_aer(n,itype,iphase),        &
                lptr_biog1_o_aer(n,itype,iphase),lptr_biog2_o_aer(n,itype,iphase),        &
                lptr_biog3_o_aer(n,itype,iphase),lptr_biog4_o_aer(n,itype,iphase),        &
                lptr_asoaX_aer(n,itype,iphase), &
                lptr_asoa1_aer(n,itype,iphase),lptr_asoa2_aer(n,itype,iphase),          &
                lptr_asoa3_aer(n,itype,iphase),lptr_asoa4_aer(n,itype,iphase),          &
                lptr_bsoaX_aer(n,itype,iphase), &
                lptr_bsoa1_aer(n,itype,iphase),lptr_bsoa2_aer(n,itype,iphase),          &
                lptr_bsoa3_aer(n,itype,iphase),lptr_bsoa4_aer(n,itype,iphase), &
		lptr_terp_Cn3_aer(n,itype,iphase) 
                         
	end do 
	end do 
	end do 
	print 95010, 'aerchem_debug_dump end -- first time'




1000	continue
	if (iflag .eq. 1) goto 1010
	if (iflag .eq. 2) goto 2000
	if (iflag .eq. 3) goto 3000
	return





1010	continue
	dtchem_sv1 = dtchem
	do m = 1, nsubareas
	do k = 1, ktot
	do l = 1, ltot2
	    rsub_sv1(l,k,m) = rsub(l,k,m)
	end do
	end do
	end do

	print 95010, 'aerchem_debug_dump start -- iflag=1'
	do m = 1, nsubareas
	do k = 1, ktot
	    print 95020, iymdcur, ihmscur,   &
		iclm, jclm, k, m, nsubareas, iflag
	    print 95030, t, dtchem_sv1, cairclm(k), relhumclm(k),   &
      		ptotclm(k), afracsubarea(k,m)
	    print 95030, (rsub_sv1(l,k,m), rsub(l,k,m), l=1,ltot2)
	end do
	end do
	print 95010, 'aerchem_debug_dump end -- iflag=1'

	return





2000	continue
	return





3000	continue
	print 95010, 'aerchem_debug_dump start -- iflag=3'
	do m = 1, nsubareas
	do k = 1, ktot
	    print 95020, iymdcur, ihmscur,   &
		iclm, jclm, k, m, nsubareas, iflag
	    print 95030, t, dtchem_sv1, cairclm(k), relhumclm(k),   &
      		ptotclm(k), afracsubarea(k,m)
	    print 95030, (rsub_sv1(l,k,m), rsub(l,k,m), l=1,ltot2)
	end do
	end do
	print 95010, 'aerchem_debug_dump end -- iflag=3'


	return
	end subroutine aerchem_debug_dump 




	end module module_mosaic_driver
