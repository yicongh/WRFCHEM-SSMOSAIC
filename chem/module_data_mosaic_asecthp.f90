      module module_data_mosaic_asecthp






      use module_data_mosaic_kind, only:  r8
      use module_data_mosaic_aero, only:  naer, nbin_a_max, ngas_aerchtot

      implicit none

























































































































































	integer, parameter :: maxd_atype = 1
        integer, parameter :: maxd_atype_md1 = 1
        integer, parameter :: maxd_atype_md2 = 1

	integer, parameter :: maxd_asize = 8
	integer, parameter :: maxd_acomp = 120 
	integer, parameter :: maxd_aphase = 2

	integer, save :: ai_phase = -999888777
	integer, save :: cw_phase = -999888777
	integer, save :: ci_phase = -999888777
	integer, save :: rn_phase = -999888777
	integer, save :: sn_phase = -999888777
	integer, save :: gr_phase = -999888777

	integer, save :: ntype_aer = 0 
	integer, save :: ntot_mastercomp_aer = 0 
	integer, save :: nphase_aer = 0 

                                                          
        integer, parameter :: identical_comps_optaa  = 1  
        integer, parameter :: identical_phases_optaa = 1  
        integer, parameter :: identical_sizes_optaa  = 1  


	integer, save ::   &
      	  nsize_aer( maxd_atype ),   & 
      	  ncomp_aer( maxd_atype ),   & 
      	  ncomp_plustracer_aer( maxd_atype ),   &
          mastercompptr_aer(maxd_acomp, maxd_atype), &         
      	  massptr_aer( maxd_acomp, maxd_asize, maxd_atype, maxd_aphase ), & 
		
      	  waterptr_aer( maxd_asize, maxd_atype ), &            
      	  hyswptr_aer( maxd_asize, maxd_atype ), &
      	  numptr_aer( maxd_asize, maxd_atype, maxd_aphase ), & 
          mprognum_aer(maxd_asize,maxd_atype,maxd_aphase)

      integer, save :: ntype_md1_aer = 0 
      integer, save :: ntype_md2_aer = 0 

      integer, save ::   &
         ibin_of_isize_itype(maxd_asize,maxd_atype), &
         isize_of_ibin(nbin_a_max), &
         itype_of_ibin(nbin_a_max), &
         itype_of_itype_md1md2(maxd_atype_md1,maxd_atype_md2), &
         itype_md1_of_itype(maxd_atype), &
         itype_md2_of_itype(maxd_atype)

      integer, save :: lunerr=-1, lunout=-1

      integer, save ::   &
         rbox_aer_ptr(-3:naer,nbin_a_max),  &   
         rbox_gas_ptr(ngas_aerchtot)            


      logical, save ::   &
          is_tracer_mastercomp_aer(  maxd_acomp )


      real(r8), save ::   &
         dens_aer( maxd_acomp, maxd_atype ),  &
         dens_mastercomp_aer( maxd_acomp ),   &
         mw_mastercomp_aer( maxd_acomp ),     &
         mw_aer( maxd_acomp, maxd_atype ),    &
         hygro_mastercomp_aer( maxd_acomp ),  &
         hygro_aer( maxd_acomp, maxd_atype )

      real(r8), save ::   &
         volumcen_sect( maxd_asize, maxd_atype ),  &
         volumlo_sect( maxd_asize, maxd_atype ),   &
         volumhi_sect( maxd_asize, maxd_atype ),   &
         dcen_sect( maxd_asize, maxd_atype ),      &
         dlo_sect( maxd_asize, maxd_atype ),       &
         dhi_sect( maxd_asize, maxd_atype ),       &
         sigmag_aer(maxd_asize, maxd_atype)

      real(r8), save ::   &
         volumcut_sect( 0:maxd_asize, maxd_atype ),  &
         dcut_sect( 0:maxd_asize, maxd_atype )





      real(r8), save ::    &
         xcut_atype_md1(0:maxd_atype_md1),   &
         xcut_atype_md2(0:maxd_atype_md2)


	character*20, save ::   &
      	  name_mastercomp_aer( maxd_acomp ),   &
      	  sname_mastercomp_aer( maxd_acomp ),  &
      	  name_aer( maxd_acomp, maxd_atype ),  &
      	  sname_aer( maxd_acomp, maxd_atype )

        character*20, save ::   &
          phase_chars_aer( maxd_acomp ),  &
          type_chars_aer( maxd_acomp )









        real,     parameter :: smallmassaa = 1.0e-22



        real,     parameter :: smallmassbb = 1.0e-32











	real, parameter :: dens_water_aer  = 1.0


        real, parameter :: dens_so4_aer = 1.80
        real, parameter :: dens_no3_aer = 1.80
        real, parameter :: dens_cl_aer  = 2.20
        real, parameter :: dens_msa_aer = 1.80
        real, parameter :: dens_co3_aer = 2.60
        real, parameter :: dens_nh4_aer = 1.80
        real, parameter :: dens_na_aer  = 2.20
        real, parameter :: dens_ca_aer  = 2.60
        real, parameter :: dens_oin_aer = 2.60





        real, parameter :: dens_oc_aer  = 1.00
        real, parameter :: dens_bc_aer  = 1.70
       real, parameter :: dens_tanv_aer  = 1.00
       real, parameter :: dens_sopnv_aer  = 1.00
       real, parameter :: dens_ternv_aer  = 1.00
       real, parameter :: dens_seqnv_aer  = 1.00
       real, parameter :: dens_sianv_aer  = 1.00
        real, parameter :: dens_pcg1_b_c_aer = 1.0
        real, parameter :: dens_pcg2_b_c_aer = 1.0
        real, parameter :: dens_pcg3_b_c_aer = 1.0
        real, parameter :: dens_pcg4_b_c_aer = 1.0
        real, parameter :: dens_pcg5_b_c_aer = 1.0
        real, parameter :: dens_pcg6_b_c_aer = 1.0
        real, parameter :: dens_pcg7_b_c_aer = 1.0
        real, parameter :: dens_pcg8_b_c_aer = 1.0
        real, parameter :: dens_pcg9_b_c_aer = 1.0
        real, parameter :: dens_pcg1_b_o_aer = 1.0
        real, parameter :: dens_pcg2_b_o_aer = 1.0
        real, parameter :: dens_pcg3_b_o_aer = 1.0
        real, parameter :: dens_pcg4_b_o_aer = 1.0
        real, parameter :: dens_pcg5_b_o_aer = 1.0
        real, parameter :: dens_pcg6_b_o_aer = 1.0
        real, parameter :: dens_pcg7_b_o_aer = 1.0
        real, parameter :: dens_pcg8_b_o_aer = 1.0
        real, parameter :: dens_pcg9_b_o_aer = 1.0
        real, parameter :: dens_opcg1_b_c_aer = 1.0
        real, parameter :: dens_opcg2_b_c_aer = 1.0
        real, parameter :: dens_opcg3_b_c_aer = 1.0
        real, parameter :: dens_opcg4_b_c_aer = 1.0
        real, parameter :: dens_opcg5_b_c_aer = 1.0
        real, parameter :: dens_opcg6_b_c_aer = 1.0
        real, parameter :: dens_opcg7_b_c_aer = 1.0
        real, parameter :: dens_opcg8_b_c_aer = 1.0
        real, parameter :: dens_opcg1_b_o_aer = 1.0
        real, parameter :: dens_opcg2_b_o_aer = 1.0
        real, parameter :: dens_opcg3_b_o_aer = 1.0
        real, parameter :: dens_opcg4_b_o_aer = 1.0
        real, parameter :: dens_opcg5_b_o_aer = 1.0
        real, parameter :: dens_opcg6_b_o_aer = 1.0
        real, parameter :: dens_opcg7_b_o_aer = 1.0
        real, parameter :: dens_opcg8_b_o_aer = 1.0
        real, parameter :: dens_pcg1_f_c_aer = 1.0
        real, parameter :: dens_pcg2_f_c_aer = 1.0
        real, parameter :: dens_pcg3_f_c_aer = 1.0
        real, parameter :: dens_pcg4_f_c_aer = 1.0
        real, parameter :: dens_pcg5_f_c_aer = 1.0
        real, parameter :: dens_pcg6_f_c_aer = 1.0
        real, parameter :: dens_pcg7_f_c_aer = 1.0
        real, parameter :: dens_pcg8_f_c_aer = 1.0
        real, parameter :: dens_pcg9_f_c_aer = 1.0
        real, parameter :: dens_pcg1_f_o_aer = 1.0
        real, parameter :: dens_pcg2_f_o_aer = 1.0
        real, parameter :: dens_pcg3_f_o_aer = 1.0
        real, parameter :: dens_pcg4_f_o_aer = 1.0
        real, parameter :: dens_pcg5_f_o_aer = 1.0
        real, parameter :: dens_pcg6_f_o_aer = 1.0
        real, parameter :: dens_pcg7_f_o_aer = 1.0
        real, parameter :: dens_pcg8_f_o_aer = 1.0
        real, parameter :: dens_pcg9_f_o_aer = 1.0
        real, parameter :: dens_opcg1_f_c_aer = 1.0
        real, parameter :: dens_opcg2_f_c_aer = 1.0
        real, parameter :: dens_opcg3_f_c_aer = 1.0
        real, parameter :: dens_opcg4_f_c_aer = 1.0
        real, parameter :: dens_opcg5_f_c_aer = 1.0
        real, parameter :: dens_opcg6_f_c_aer = 1.0
        real, parameter :: dens_opcg7_f_c_aer = 1.0
        real, parameter :: dens_opcg8_f_c_aer = 1.0
        real, parameter :: dens_opcg1_f_o_aer = 1.0
        real, parameter :: dens_opcg2_f_o_aer = 1.0
        real, parameter :: dens_opcg3_f_o_aer = 1.0
        real, parameter :: dens_opcg4_f_o_aer = 1.0
        real, parameter :: dens_opcg5_f_o_aer = 1.0
        real, parameter :: dens_opcg6_f_o_aer = 1.0
        real, parameter :: dens_opcg7_f_o_aer = 1.0
        real, parameter :: dens_opcg8_f_o_aer = 1.0
        real, parameter :: dens_smpa_aer = 1.0
        real, parameter :: dens_smpbb_aer = 1.0
       real, parameter :: dens_iepox_aer = 1.0
       real, parameter :: dens_iepoxos_aer = 1.0
       real, parameter :: dens_tetrol_aer = 1.0
       real, parameter :: dens_gly_aer = 1.0
        real, parameter :: dens_ant1_c_aer = 1.0
        real, parameter :: dens_ant2_c_aer = 1.0
        real, parameter :: dens_ant3_c_aer = 1.0
        real, parameter :: dens_ant4_c_aer = 1.0
        real, parameter :: dens_ant1_o_aer = 1.0
        real, parameter :: dens_ant2_o_aer = 1.0
        real, parameter :: dens_ant3_o_aer = 1.0
        real, parameter :: dens_ant4_o_aer = 1.0
        real, parameter :: dens_biog1_c_aer = 1.0
        real, parameter :: dens_biog2_c_aer = 1.0
        real, parameter :: dens_biog3_c_aer = 1.0
        real, parameter :: dens_biog4_c_aer = 1.0
        real, parameter :: dens_biog1_o_aer = 1.0
        real, parameter :: dens_biog2_o_aer = 1.0
        real, parameter :: dens_biog3_o_aer = 1.0
        real, parameter :: dens_biog4_o_aer = 1.0

        real, parameter :: dens_glysoa_r1_aer   = 1.0
        real, parameter :: dens_glysoa_r2_aer   = 1.0
        real, parameter :: dens_glysoa_oh_aer   = 1.0
        real, parameter :: dens_glysoa_nh4_aer   = 1.0
        real, parameter :: dens_glysoa_sfc_aer   = 1.0

        real, parameter :: dens_asoaX_aer = 1.5
        real, parameter :: dens_asoa1_aer = 1.5
        real, parameter :: dens_asoa2_aer = 1.5
        real, parameter :: dens_asoa3_aer = 1.5
        real, parameter :: dens_asoa4_aer = 1.5
        real, parameter :: dens_bsoaX_aer = 1.5
        real, parameter :: dens_bsoa1_aer = 1.5
        real, parameter :: dens_bsoa2_aer = 1.5
        real, parameter :: dens_bsoa3_aer = 1.5
        real, parameter :: dens_bsoa4_aer = 1.5



        real, parameter :: mw_oc_aer  = 250.0   

        real, parameter :: mw_so4_aer = 96.066
        real, parameter :: mw_no3_aer = 62.007
        real, parameter :: mw_cl_aer  = 35.450
        real, parameter :: mw_msa_aer = 96.109
        real, parameter :: mw_co3_aer = 60.007
        real, parameter :: mw_nh4_aer = 18.042
        real, parameter :: mw_na_aer  = 22.990
        real, parameter :: mw_ca_aer  = 40.080
        real, parameter :: mw_oin_aer = 1.0

        real, parameter :: mw_bc_aer  = 1.0
        real, parameter :: mw_water_aer  = 18.016





       real, parameter :: mw_tanv_aer  = 250.0
       real, parameter :: mw_sopnv_aer  = 250.0
       real, parameter :: mw_ternv_aer  = 250.0
       real, parameter :: mw_seqnv_aer  = 250.0
       real, parameter :: mw_sianv_aer  = 250.0
        real, parameter :: mw_pcg1_b_c_aer = 250.0
        real, parameter :: mw_pcg2_b_c_aer = 250.0
        real, parameter :: mw_pcg3_b_c_aer = 250.0
        real, parameter :: mw_pcg4_b_c_aer = 250.0
        real, parameter :: mw_pcg5_b_c_aer = 250.0
        real, parameter :: mw_pcg6_b_c_aer = 250.0
        real, parameter :: mw_pcg7_b_c_aer = 250.0
        real, parameter :: mw_pcg8_b_c_aer = 250.0
        real, parameter :: mw_pcg9_b_c_aer = 250.0
        real, parameter :: mw_pcg1_b_o_aer = 250.0
        real, parameter :: mw_pcg2_b_o_aer = 250.0
        real, parameter :: mw_pcg3_b_o_aer = 250.0
        real, parameter :: mw_pcg4_b_o_aer = 250.0
        real, parameter :: mw_pcg5_b_o_aer = 250.0
        real, parameter :: mw_pcg6_b_o_aer = 250.0
        real, parameter :: mw_pcg7_b_o_aer = 250.0
        real, parameter :: mw_pcg8_b_o_aer = 250.0
        real, parameter :: mw_pcg9_b_o_aer = 250.0
        real, parameter :: mw_opcg1_b_c_aer = 250.0
        real, parameter :: mw_opcg2_b_c_aer = 250.0
        real, parameter :: mw_opcg3_b_c_aer = 250.0
        real, parameter :: mw_opcg4_b_c_aer = 250.0
        real, parameter :: mw_opcg5_b_c_aer = 250.0
        real, parameter :: mw_opcg6_b_c_aer = 250.0
        real, parameter :: mw_opcg7_b_c_aer = 250.0
        real, parameter :: mw_opcg8_b_c_aer = 250.0
        real, parameter :: mw_opcg1_b_o_aer = 250.0
        real, parameter :: mw_opcg2_b_o_aer = 250.0
        real, parameter :: mw_opcg3_b_o_aer = 250.0
        real, parameter :: mw_opcg4_b_o_aer = 250.0
        real, parameter :: mw_opcg5_b_o_aer = 250.0
        real, parameter :: mw_opcg6_b_o_aer = 250.0
        real, parameter :: mw_opcg7_b_o_aer = 250.0
        real, parameter :: mw_opcg8_b_o_aer = 250.0
        real, parameter :: mw_pcg1_f_c_aer = 250.0
        real, parameter :: mw_pcg2_f_c_aer = 250.0
        real, parameter :: mw_pcg3_f_c_aer = 250.0
        real, parameter :: mw_pcg4_f_c_aer = 250.0
        real, parameter :: mw_pcg5_f_c_aer = 250.0
        real, parameter :: mw_pcg6_f_c_aer = 250.0
        real, parameter :: mw_pcg7_f_c_aer = 250.0
        real, parameter :: mw_pcg8_f_c_aer = 250.0
        real, parameter :: mw_pcg9_f_c_aer = 250.0
        real, parameter :: mw_pcg1_f_o_aer = 250.0
        real, parameter :: mw_pcg2_f_o_aer = 250.0
        real, parameter :: mw_pcg3_f_o_aer = 250.0
        real, parameter :: mw_pcg4_f_o_aer = 250.0
        real, parameter :: mw_pcg5_f_o_aer = 250.0
        real, parameter :: mw_pcg6_f_o_aer = 250.0
        real, parameter :: mw_pcg7_f_o_aer = 250.0
        real, parameter :: mw_pcg8_f_o_aer = 250.0
        real, parameter :: mw_pcg9_f_o_aer = 250.0
        real, parameter :: mw_opcg1_f_c_aer = 250.0
        real, parameter :: mw_opcg2_f_c_aer = 250.0
        real, parameter :: mw_opcg3_f_c_aer = 250.0
        real, parameter :: mw_opcg4_f_c_aer = 250.0
        real, parameter :: mw_opcg5_f_c_aer = 250.0
        real, parameter :: mw_opcg6_f_c_aer = 250.0
        real, parameter :: mw_opcg7_f_c_aer = 250.0
        real, parameter :: mw_opcg8_f_c_aer = 250.0
        real, parameter :: mw_opcg1_f_o_aer = 250.0
        real, parameter :: mw_opcg2_f_o_aer = 250.0
        real, parameter :: mw_opcg3_f_o_aer = 250.0
        real, parameter :: mw_opcg4_f_o_aer = 250.0
        real, parameter :: mw_opcg5_f_o_aer = 250.0
        real, parameter :: mw_opcg6_f_o_aer = 250.0
        real, parameter :: mw_opcg7_f_o_aer = 250.0
        real, parameter :: mw_opcg8_f_o_aer = 250.0
        real, parameter :: mw_smpa_aer = 250.0
        real, parameter :: mw_smpbb_aer = 250.0
       real, parameter :: mw_iepox_aer = 118.0
       real, parameter :: mw_iepoxos_aer = 216.0
       real, parameter :: mw_tetrol_aer = 136.0
       real, parameter :: mw_gly_aer = 58.0
        real, parameter :: mw_ant1_c_aer = 250.0
        real, parameter :: mw_ant2_c_aer = 250.0
        real, parameter :: mw_ant3_c_aer = 250.0
        real, parameter :: mw_ant4_c_aer = 250.0
        real, parameter :: mw_ant1_o_aer = 250.0
        real, parameter :: mw_ant2_o_aer = 250.0
        real, parameter :: mw_ant3_o_aer = 250.0
        real, parameter :: mw_ant4_o_aer = 250.0
        real, parameter :: mw_biog1_c_aer = 250.0
        real, parameter :: mw_biog2_c_aer = 250.0
        real, parameter :: mw_biog3_c_aer = 250.0
        real, parameter :: mw_biog4_c_aer = 250.0
        real, parameter :: mw_biog1_o_aer = 250.0
        real, parameter :: mw_biog2_o_aer = 250.0
        real, parameter :: mw_biog3_o_aer = 250.0
        real, parameter :: mw_biog4_o_aer = 250.0

        real, parameter :: mw_glysoa_r1_aer   = 250.0
        real, parameter :: mw_glysoa_r2_aer   = 250.0
        real, parameter :: mw_glysoa_oh_aer   = 250.0
        real, parameter :: mw_glysoa_nh4_aer   = 250.0
        real, parameter :: mw_glysoa_sfc_aer   = 250.0
        real, parameter :: mw_asoaX_aer = 250.0
        real, parameter :: mw_asoa1_aer = 250.0
        real, parameter :: mw_asoa2_aer = 250.0
        real, parameter :: mw_asoa3_aer = 250.0
        real, parameter :: mw_asoa4_aer = 250.0
        real, parameter :: mw_bsoaX_aer = 250.0
        real, parameter :: mw_bsoa1_aer = 250.0
        real, parameter :: mw_bsoa2_aer = 250.0
        real, parameter :: mw_bsoa3_aer = 250.0
        real, parameter :: mw_bsoa4_aer = 250.0


        real, parameter :: hygro_so4_aer = 0.5
        real, parameter :: hygro_no3_aer = 0.5
        real, parameter :: hygro_ca_aer  = 0.1
        real, parameter :: hygro_co3_aer = 0.1
        real, parameter :: hygro_nh4_aer = 0.5
        real, parameter :: hygro_msa_aer = 0.58
        real, parameter :: hygro_cl_aer  = 1.16
        real, parameter :: hygro_na_aer  = 1.16
        real, parameter :: hygro_oin_aer = 0.068 






        real, parameter :: hygro_oc_aer  = 1.0e-4 
        real, parameter :: hygro_bc_aer  = 1.0e-6
        real, parameter :: hygro_smpa_aer  = 0.10
        real, parameter :: hygro_smpbb_aer  = 0.140
       real, parameter :: hygro_iepox_aer  = 0.10
       real, parameter :: hygro_iepoxos_aer  = 0.10
       real, parameter :: hygro_tetrol_aer  = 0.10
       real, parameter :: hygro_gly_aer  = 0.10
       real, parameter :: hygro_tanv_aer  = 0.10
       real, parameter :: hygro_sopnv_aer  = 0.10
       real, parameter :: hygro_ternv_aer  = 0.10
       real, parameter :: hygro_seqnv_aer  = 0.10
       real, parameter :: hygro_sianv_aer  = 0.10
        real, parameter :: hygro_pcg1_b_c_aer = 0.04
        real, parameter :: hygro_pcg2_b_c_aer = 0.04
        real, parameter :: hygro_pcg3_b_c_aer = 0.04
        real, parameter :: hygro_pcg4_b_c_aer = 0.04
        real, parameter :: hygro_pcg5_b_c_aer = 0.04
        real, parameter :: hygro_pcg6_b_c_aer = 0.04
        real, parameter :: hygro_pcg7_b_c_aer = 0.04
        real, parameter :: hygro_pcg8_b_c_aer = 0.04
        real, parameter :: hygro_pcg9_b_c_aer = 0.04
        real, parameter :: hygro_pcg1_b_o_aer = 0.04
        real, parameter :: hygro_pcg2_b_o_aer = 0.04
        real, parameter :: hygro_pcg3_b_o_aer = 0.04
        real, parameter :: hygro_pcg4_b_o_aer = 0.04
        real, parameter :: hygro_pcg5_b_o_aer = 0.04
        real, parameter :: hygro_pcg6_b_o_aer = 0.04
        real, parameter :: hygro_pcg7_b_o_aer = 0.04
        real, parameter :: hygro_pcg8_b_o_aer = 0.04
        real, parameter :: hygro_pcg9_b_o_aer = 0.04
        real, parameter :: hygro_opcg1_b_c_aer = 0.10
        real, parameter :: hygro_opcg2_b_c_aer = 0.10
        real, parameter :: hygro_opcg3_b_c_aer = 0.10
        real, parameter :: hygro_opcg4_b_c_aer = 0.10
        real, parameter :: hygro_opcg5_b_c_aer = 0.10
        real, parameter :: hygro_opcg6_b_c_aer = 0.10
        real, parameter :: hygro_opcg7_b_c_aer = 0.10
        real, parameter :: hygro_opcg8_b_c_aer = 0.10
        real, parameter :: hygro_opcg1_b_o_aer = 0.10
        real, parameter :: hygro_opcg2_b_o_aer = 0.10
        real, parameter :: hygro_opcg3_b_o_aer = 0.10
        real, parameter :: hygro_opcg4_b_o_aer = 0.10
        real, parameter :: hygro_opcg5_b_o_aer = 0.10
        real, parameter :: hygro_opcg6_b_o_aer = 0.10
        real, parameter :: hygro_opcg7_b_o_aer = 0.10
        real, parameter :: hygro_opcg8_b_o_aer = 0.10
        real, parameter :: hygro_pcg1_f_c_aer = 1.0e-6
        real, parameter :: hygro_pcg2_f_c_aer = 1.0e-6
        real, parameter :: hygro_pcg3_f_c_aer = 1.0e-6
        real, parameter :: hygro_pcg4_f_c_aer = 1.0e-6
        real, parameter :: hygro_pcg5_f_c_aer = 1.0e-6
        real, parameter :: hygro_pcg6_f_c_aer = 1.0e-6
        real, parameter :: hygro_pcg7_f_c_aer = 1.0e-6
        real, parameter :: hygro_pcg8_f_c_aer = 1.0e-6
        real, parameter :: hygro_pcg9_f_c_aer = 1.0e-6
        real, parameter :: hygro_pcg1_f_o_aer = 1.0e-6
        real, parameter :: hygro_pcg2_f_o_aer = 1.0e-6
        real, parameter :: hygro_pcg3_f_o_aer = 1.0e-6
        real, parameter :: hygro_pcg4_f_o_aer = 1.0e-6
        real, parameter :: hygro_pcg5_f_o_aer = 1.0e-6
        real, parameter :: hygro_pcg6_f_o_aer = 1.0e-6
        real, parameter :: hygro_pcg7_f_o_aer = 1.0e-6
        real, parameter :: hygro_pcg8_f_o_aer = 1.0e-6
        real, parameter :: hygro_pcg9_f_o_aer = 1.0e-6
        real, parameter :: hygro_opcg1_f_c_aer = 0.10
        real, parameter :: hygro_opcg2_f_c_aer = 0.10
        real, parameter :: hygro_opcg3_f_c_aer = 0.10
        real, parameter :: hygro_opcg4_f_c_aer = 0.10
        real, parameter :: hygro_opcg5_f_c_aer = 0.10
        real, parameter :: hygro_opcg6_f_c_aer = 0.10
        real, parameter :: hygro_opcg7_f_c_aer = 0.10
        real, parameter :: hygro_opcg8_f_c_aer = 0.10
        real, parameter :: hygro_opcg1_f_o_aer = 0.10
        real, parameter :: hygro_opcg2_f_o_aer = 0.10
        real, parameter :: hygro_opcg3_f_o_aer = 0.10
        real, parameter :: hygro_opcg4_f_o_aer = 0.10
        real, parameter :: hygro_opcg5_f_o_aer = 0.10
        real, parameter :: hygro_opcg6_f_o_aer = 0.10
        real, parameter :: hygro_opcg7_f_o_aer = 0.10
        real, parameter :: hygro_opcg8_f_o_aer = 0.10
        real, parameter :: hygro_ant1_c_aer = 0.10
        real, parameter :: hygro_ant2_c_aer = 0.10
        real, parameter :: hygro_ant3_c_aer = 0.10
        real, parameter :: hygro_ant4_c_aer = 0.10
        real, parameter :: hygro_ant1_o_aer = 0.10
        real, parameter :: hygro_ant2_o_aer = 0.10
        real, parameter :: hygro_ant3_o_aer = 0.10
        real, parameter :: hygro_ant4_o_aer = 0.10
        real, parameter :: hygro_biog1_c_aer = 0.10
        real, parameter :: hygro_biog2_c_aer = 0.10
        real, parameter :: hygro_biog3_c_aer = 0.10
        real, parameter :: hygro_biog4_c_aer = 0.10
        real, parameter :: hygro_biog1_o_aer = 0.10
        real, parameter :: hygro_biog2_o_aer = 0.10
        real, parameter :: hygro_biog3_o_aer = 0.10
        real, parameter :: hygro_biog4_o_aer = 0.10

        real, parameter :: hygro_glysoa_r1_aer   = 0.14
        real, parameter :: hygro_glysoa_r2_aer   = 0.14
        real, parameter :: hygro_glysoa_oh_aer   = 0.14
        real, parameter :: hygro_glysoa_nh4_aer   = 0.14
        real, parameter :: hygro_glysoa_sfc_aer   = 0.14
        real, parameter :: hygro_asoaX_aer = 0.14
        real, parameter :: hygro_asoa1_aer = 0.14
        real, parameter :: hygro_asoa2_aer = 0.14
        real, parameter :: hygro_asoa3_aer = 0.14
        real, parameter :: hygro_asoa4_aer = 0.14
        real, parameter :: hygro_bsoaX_aer = 0.14
        real, parameter :: hygro_bsoa1_aer = 0.14
        real, parameter :: hygro_bsoa2_aer = 0.14
        real, parameter :: hygro_bsoa3_aer = 0.14
        real, parameter :: hygro_bsoa4_aer = 0.14



        integer, save :: mastercompindx_so4_aer = -999888777
        integer, save :: mastercompindx_no3_aer = -999888777
        integer, save :: mastercompindx_cl_aer  = -999888777
        integer, save :: mastercompindx_msa_aer = -999888777
        integer, save :: mastercompindx_co3_aer = -999888777
        integer, save :: mastercompindx_nh4_aer = -999888777
        integer, save :: mastercompindx_na_aer  = -999888777
        integer, save :: mastercompindx_ca_aer  = -999888777
        integer, save :: mastercompindx_oin_aer = -999888777





        integer, save :: mastercompindx_oc_aer  = -999888777
        integer, save :: mastercompindx_bc_aer  = -999888777
       integer, save :: mastercompindx_tanv_aer  = -999888777
       integer, save :: mastercompindx_sopnv_aer  = -999888777
       integer, save :: mastercompindx_ternv_aer  = -999888777
       integer, save :: mastercompindx_seqnv_aer  = -999888777
       integer, save :: mastercompindx_sianv_aer  = -999888777
        integer, save :: mastercompindx_pcg1_b_c_aer = -999888777
        integer, save :: mastercompindx_pcg2_b_c_aer = -999888777
        integer, save :: mastercompindx_pcg3_b_c_aer = -999888777
        integer, save :: mastercompindx_pcg4_b_c_aer = -999888777
        integer, save :: mastercompindx_pcg5_b_c_aer = -999888777
        integer, save :: mastercompindx_pcg6_b_c_aer = -999888777
        integer, save :: mastercompindx_pcg7_b_c_aer = -999888777
        integer, save :: mastercompindx_pcg8_b_c_aer = -999888777
        integer, save :: mastercompindx_pcg9_b_c_aer = -999888777
        integer, save :: mastercompindx_pcg1_b_o_aer = -999888777
        integer, save :: mastercompindx_pcg2_b_o_aer = -999888777
        integer, save :: mastercompindx_pcg3_b_o_aer = -999888777
        integer, save :: mastercompindx_pcg4_b_o_aer = -999888777
        integer, save :: mastercompindx_pcg5_b_o_aer = -999888777
        integer, save :: mastercompindx_pcg6_b_o_aer = -999888777
        integer, save :: mastercompindx_pcg7_b_o_aer = -999888777
        integer, save :: mastercompindx_pcg8_b_o_aer = -999888777
        integer, save :: mastercompindx_pcg9_b_o_aer = -999888777
        integer, save :: mastercompindx_opcg1_b_c_aer = -999888777
        integer, save :: mastercompindx_opcg2_b_c_aer = -999888777
        integer, save :: mastercompindx_opcg3_b_c_aer = -999888777
        integer, save :: mastercompindx_opcg4_b_c_aer = -999888777
        integer, save :: mastercompindx_opcg5_b_c_aer = -999888777
        integer, save :: mastercompindx_opcg6_b_c_aer = -999888777
        integer, save :: mastercompindx_opcg7_b_c_aer = -999888777
        integer, save :: mastercompindx_opcg8_b_c_aer = -999888777
        integer, save :: mastercompindx_opcg1_b_o_aer = -999888777
        integer, save :: mastercompindx_opcg2_b_o_aer = -999888777
        integer, save :: mastercompindx_opcg3_b_o_aer = -999888777
        integer, save :: mastercompindx_opcg4_b_o_aer = -999888777
        integer, save :: mastercompindx_opcg5_b_o_aer = -999888777
        integer, save :: mastercompindx_opcg6_b_o_aer = -999888777
        integer, save :: mastercompindx_opcg7_b_o_aer = -999888777
        integer, save :: mastercompindx_opcg8_b_o_aer = -999888777
        integer, save :: mastercompindx_pcg1_f_c_aer = -999888777
        integer, save :: mastercompindx_pcg2_f_c_aer = -999888777
        integer, save :: mastercompindx_pcg3_f_c_aer = -999888777
        integer, save :: mastercompindx_pcg4_f_c_aer = -999888777
        integer, save :: mastercompindx_pcg5_f_c_aer = -999888777
        integer, save :: mastercompindx_pcg6_f_c_aer = -999888777
        integer, save :: mastercompindx_pcg7_f_c_aer = -999888777
        integer, save :: mastercompindx_pcg8_f_c_aer = -999888777
        integer, save :: mastercompindx_pcg9_f_c_aer = -999888777
        integer, save :: mastercompindx_pcg1_f_o_aer = -999888777
        integer, save :: mastercompindx_pcg2_f_o_aer = -999888777
        integer, save :: mastercompindx_pcg3_f_o_aer = -999888777
        integer, save :: mastercompindx_pcg4_f_o_aer = -999888777
        integer, save :: mastercompindx_pcg5_f_o_aer = -999888777
        integer, save :: mastercompindx_pcg6_f_o_aer = -999888777
        integer, save :: mastercompindx_pcg7_f_o_aer = -999888777
        integer, save :: mastercompindx_pcg8_f_o_aer = -999888777
        integer, save :: mastercompindx_pcg9_f_o_aer = -999888777
        integer, save :: mastercompindx_opcg1_f_c_aer = -999888777
        integer, save :: mastercompindx_opcg2_f_c_aer = -999888777
        integer, save :: mastercompindx_opcg3_f_c_aer = -999888777
        integer, save :: mastercompindx_opcg4_f_c_aer = -999888777
        integer, save :: mastercompindx_opcg5_f_c_aer = -999888777
        integer, save :: mastercompindx_opcg6_f_c_aer = -999888777
        integer, save :: mastercompindx_opcg7_f_c_aer = -999888777
        integer, save :: mastercompindx_opcg8_f_c_aer = -999888777
        integer, save :: mastercompindx_opcg1_f_o_aer = -999888777
        integer, save :: mastercompindx_opcg2_f_o_aer = -999888777
        integer, save :: mastercompindx_opcg3_f_o_aer = -999888777
        integer, save :: mastercompindx_opcg4_f_o_aer = -999888777
        integer, save :: mastercompindx_opcg5_f_o_aer = -999888777
        integer, save :: mastercompindx_opcg6_f_o_aer = -999888777
        integer, save :: mastercompindx_opcg7_f_o_aer = -999888777
        integer, save :: mastercompindx_opcg8_f_o_aer = -999888777
        integer, save :: mastercompindx_smpa_aer = -999888777
        integer, save :: mastercompindx_smpbb_aer = -999888777
       integer, save :: mastercompindx_iepox_aer = -999888777
       integer, save :: mastercompindx_iepoxos_aer = -999888777
       integer, save :: mastercompindx_tetrol_aer = -999888777
       integer, save :: mastercompindx_gly_aer = -999888777
        integer, save :: mastercompindx_ant1_c_aer = -999888777
        integer, save :: mastercompindx_ant2_c_aer = -999888777
        integer, save :: mastercompindx_ant3_c_aer = -999888777
        integer, save :: mastercompindx_ant4_c_aer = -999888777
        integer, save :: mastercompindx_ant1_o_aer = -999888777
        integer, save :: mastercompindx_ant2_o_aer = -999888777
        integer, save :: mastercompindx_ant3_o_aer = -999888777
        integer, save :: mastercompindx_ant4_o_aer = -999888777
        integer, save :: mastercompindx_biog1_c_aer = -999888777
        integer, save :: mastercompindx_biog2_c_aer = -999888777
        integer, save :: mastercompindx_biog3_c_aer = -999888777
        integer, save :: mastercompindx_biog4_c_aer = -999888777
        integer, save :: mastercompindx_biog1_o_aer = -999888777
        integer, save :: mastercompindx_biog2_o_aer = -999888777
        integer, save :: mastercompindx_biog3_o_aer = -999888777
        integer, save :: mastercompindx_biog4_o_aer = -999888777

        integer, save :: mastercompindx_glysoa_r1_aer = -999888777
        integer, save :: mastercompindx_glysoa_r2_aer = -999888777
        integer, save :: mastercompindx_glysoa_oh_aer = -999888777
        integer, save :: mastercompindx_glysoa_nh4_aer = -999888777
        integer, save :: mastercompindx_glysoa_sfc_aer = -999888777
        integer, save :: mastercompindx_asoaX_aer = -999888777
        integer, save :: mastercompindx_asoa1_aer = -999888777
        integer, save :: mastercompindx_asoa2_aer = -999888777
        integer, save :: mastercompindx_asoa3_aer = -999888777
        integer, save :: mastercompindx_asoa4_aer = -999888777
        integer, save :: mastercompindx_bsoaX_aer = -999888777
        integer, save :: mastercompindx_bsoa1_aer = -999888777
        integer, save :: mastercompindx_bsoa2_aer = -999888777
        integer, save :: mastercompindx_bsoa3_aer = -999888777
        integer, save :: mastercompindx_bsoa4_aer = -999888777



        integer, save ::                     &
          lptr_so4_aer(maxd_asize, maxd_atype, maxd_aphase),      &
          lptr_msa_aer(maxd_asize, maxd_atype, maxd_aphase),      &
          lptr_no3_aer(maxd_asize, maxd_atype, maxd_aphase),      &
          lptr_cl_aer(maxd_asize, maxd_atype, maxd_aphase),       &
          lptr_co3_aer(maxd_asize, maxd_atype, maxd_aphase),      &
          lptr_nh4_aer(maxd_asize, maxd_atype, maxd_aphase),      &
          lptr_na_aer(maxd_asize, maxd_atype, maxd_aphase),       &
          lptr_ca_aer(maxd_asize, maxd_atype, maxd_aphase),       &
          lptr_oin_aer(maxd_asize, maxd_atype, maxd_aphase),      &
          lptr_aro1_aer(maxd_asize, maxd_atype, maxd_aphase),     &
          lptr_aro2_aer(maxd_asize, maxd_atype, maxd_aphase),     &
          lptr_alk1_aer(maxd_asize, maxd_atype, maxd_aphase),     &
          lptr_ole1_aer(maxd_asize, maxd_atype, maxd_aphase),     &
          lptr_api1_aer(maxd_asize, maxd_atype, maxd_aphase),     &
          lptr_api2_aer(maxd_asize, maxd_atype, maxd_aphase),     &
          lptr_lim1_aer(maxd_asize, maxd_atype, maxd_aphase),     &
          lptr_lim2_aer(maxd_asize, maxd_atype, maxd_aphase),     &





          lptr_oc_aer(maxd_asize, maxd_atype, maxd_aphase),       &
          lptr_bc_aer(maxd_asize, maxd_atype, maxd_aphase),       &
         lptr_tanv_aer(maxd_asize, maxd_atype, maxd_aphase),     &
         lptr_sopnv_aer(maxd_asize, maxd_atype, maxd_aphase),    &
         lptr_ternv_aer(maxd_asize, maxd_atype, maxd_aphase),    &
         lptr_seqnv_aer(maxd_asize, maxd_atype, maxd_aphase),    &
         lptr_sianv_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_pcg1_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg2_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg3_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg4_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg5_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg6_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg7_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg8_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg9_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg1_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg2_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg3_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg4_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg5_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg6_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg7_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg8_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg9_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg1_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg2_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg3_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg4_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg5_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg6_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg7_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg8_b_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg1_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg2_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg3_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg4_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg5_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg6_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg7_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg8_b_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg1_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg2_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg3_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg4_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg5_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg6_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg7_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg8_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg9_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg1_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg2_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg3_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg4_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg5_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg6_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg7_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg8_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_pcg9_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg1_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg2_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg3_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg4_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg5_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg6_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg7_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg8_f_c_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg1_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg2_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg3_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg4_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg5_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg6_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg7_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_opcg8_f_o_aer(maxd_asize, maxd_atype, maxd_aphase), &
          lptr_smpa_aer(maxd_asize, maxd_atype, maxd_aphase),      &
          lptr_smpbb_aer(maxd_asize, maxd_atype, maxd_aphase),     &
         lptr_iepox_aer(maxd_asize, maxd_atype, maxd_aphase),     &
         lptr_iepoxos_aer(maxd_asize, maxd_atype, maxd_aphase),   &
         lptr_tetrol_aer(maxd_asize, maxd_atype, maxd_aphase),    &
         lptr_gly_aer(maxd_asize, maxd_atype, maxd_aphase),       &
          lptr_ant1_c_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_ant2_c_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_ant3_c_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_ant4_c_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_ant1_o_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_ant2_o_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_ant3_o_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_ant4_o_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_biog1_c_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_biog2_c_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_biog3_c_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_biog4_c_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_biog1_o_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_biog2_o_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_biog3_o_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_biog4_o_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_glysoa_r1_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_glysoa_r2_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_glysoa_oh_aer(maxd_asize, maxd_atype, maxd_aphase),   &
          lptr_glysoa_nh4_aer(maxd_asize, maxd_atype, maxd_aphase),  &
          lptr_glysoa_sfc_aer(maxd_asize, maxd_atype, maxd_aphase),  &
          lptr_asoaX_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_asoa1_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_asoa2_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_asoa3_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_asoa4_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_bsoaX_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_bsoa1_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_bsoa2_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_bsoa3_aer(maxd_asize, maxd_atype, maxd_aphase),    &
          lptr_bsoa4_aer(maxd_asize, maxd_atype, maxd_aphase)


      end module module_data_mosaic_asecthp

