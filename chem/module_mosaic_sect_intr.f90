module module_mosaic_sect_intr


















  use module_data_mosaic_kind, only:  r8


  implicit none


  integer, parameter :: iunits_flagaa = 4
  integer, parameter :: iunits_flagbb = 1
  integer, parameter :: iunits_diagout_flagaa = 1

  
  real(r8), save, allocatable ::   &
	  dens_aer_tmp(:,:),  &
	  volumcut_sect_tmp(:,:),  &
	  volumcen_sect_tmp(:,:),  &
	  volumlo_sect_tmp(:,:),  &
	  volumhi_sect_tmp(:,:),  &
	  dcut_sect_tmp(:,:),  &
	  dcen_sect_tmp(:,:),  &
	  dlo_sect_tmp(:,:),  &
	  dhi_sect_tmp(:,:)
  

  contains


  

	subroutine sect_intr_allocate_memory

	use module_data_mosaic_asecthp, only:  maxd_acomp, maxd_asize, maxd_atype

	allocate( dens_aer_tmp( maxd_acomp, maxd_atype ) )
	allocate( volumcut_sect_tmp( 0:maxd_asize, maxd_atype ) )
	allocate( volumcen_sect_tmp(   maxd_asize, maxd_atype ) )
	allocate( volumlo_sect_tmp(    maxd_asize, maxd_atype ) )
	allocate( volumhi_sect_tmp(    maxd_asize, maxd_atype ) )
	allocate( dcut_sect_tmp( 0:maxd_asize, maxd_atype ) )
	allocate( dcen_sect_tmp(   maxd_asize, maxd_atype ) )
	allocate( dlo_sect_tmp(    maxd_asize, maxd_atype ) )
	allocate( dhi_sect_tmp(    maxd_asize, maxd_atype ) )

	return
	end subroutine sect_intr_allocate_memory
        

  
  subroutine sectional_interface_1( &
       idiagbb_host, &
       id, ii, jj, kk, ktau, ktauc, dtchem, &
       rbox_sv1, rbox, &
       it_mosaic, jaerosolstate_bgn, jaerosolstate, &
       jhyst_leg, jhyst_leg_sv1, &
       dp_dry_a, dp_dry_a_sv1, &
       sigmag_a, sigmag_a_sv1, &
       mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a, &
       fact_apmassmr, fact_apnumbmr, fact_apdens_in, fact_apdiam_in, fact_gasmr, &
       pr_atm, rh, te, cair_mol_m3, cair_mol_cc )

    use module_data_mosaic_boxmod, only: idiag_sect_coag, idiag_sect_movesect, idiag_sect_newnuc
    use module_data_mosaic_main, only: mw_air, ntot_max, ntot_used
    use module_data_mosaic_aero, only: &
       inh4_a, iso4_a, jh2o, &
       dens_aer_mac, mw_aer_mac, mw_comp_a, &
       mcoag_flag1, mmovesect_flag1, mnewnuc_flag1, &
       naer, naercomp, nbin_a, nbin_a_max
    use module_data_mosaic_asecthp, only:  ai_phase, &
       hygro_so4_aer, itype_of_itype_md1md2, &
       maxd_asize, maxd_atype, nsize_aer, ntype_aer, ntype_md2_aer, &
       xcut_atype_md2
    use module_mosaic_movesect1d, only:  move_sections_x3
    use module_mosaic_movesect3d, only:  move_sect_3d_x1
    use module_mosaic_coag1d,     only:  mosaic_coag_1box
    use module_mosaic_coag3d,     only:  mosaic_coag_3d_1box
    use module_mosaic_newnucb,    only:  mosaic_newnuc_1box

    implicit none

    
    integer, intent(in) :: idiagbb_host
    integer, intent(in) :: id, ii, jj, kk, ktau, ktauc, it_mosaic
    integer, intent(in),    dimension(nbin_a_max) :: jaerosolstate, jaerosolstate_bgn
    integer, intent(in),    dimension(nbin_a_max) :: jhyst_leg_sv1
    integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg

    real(r8), intent(in)    :: dtchem  
    real(r8), intent(in)    :: pr_atm  
    real(r8), intent(in)    :: rh      
    real(r8), intent(in)    :: te      
    real(r8), intent(in)    :: cair_mol_m3, cair_mol_cc  
    real(r8), intent(in)    :: fact_apmassmr, fact_apnumbmr, fact_apdens_in, fact_apdiam_in, fact_gasmr
    real(r8), intent(in)    :: rbox_sv1(ntot_used)
    real(r8), intent(in),    dimension(nbin_a_max) :: dp_dry_a_sv1, sigmag_a_sv1
    real(r8), intent(inout) :: rbox(ntot_used)  
    real(r8), intent(inout), dimension(nbin_a_max) :: dens_dry_a_bgn, dens_dry_a  
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_dry_a_bgn, mass_dry_a  
    real(r8), intent(inout), dimension(nbin_a_max) :: dp_dry_a, sigmag_a          

    
    integer  :: idiagbb_coag, istat_coag
    integer  :: idiagbb_newnuc, istat_newnuc, itype_newnuc, it2
    integer  :: idiagbb_movesect
    integer  :: jhyst_leg_sv2(nbin_a_max)
    integer  :: jsv
    integer  :: method_coag, method_movesect

    logical :: ldiag1, ldiag2
    logical :: ldiag_movesect, ldiag_newnuc, ldiag_coag

    real(r8) :: dens_nh4so4a_newnuc
    real(r8) :: fact_apdiam, fact_apdens
    real(r8) :: rh_box, rhoair_g_cc
    real(r8) :: tmpa

    real(r8) :: cnn(ntot_used), cnn_sv1(ntot_used), cnn_sv2(ntot_used)
    real(r8) :: rbox_sv2(ntot_used)
    real(r8) :: rbox_svx(ntot_used,4)

    real(r8) :: dp_dry_a_sv2(nbin_a_max)
    real(r8) :: sigmag_a_sv2(nbin_a_max)

    real(r8) :: drydens_pregrow(maxd_asize,maxd_atype)
    real(r8) :: drydens_aftgrow(maxd_asize,maxd_atype)
    real(r8) :: drymass_pregrow(maxd_asize,maxd_atype)
    real(r8) :: drymass_aftgrow(maxd_asize,maxd_atype)
    
    
    
    
    
    

    real(r8) :: adrydens_box(maxd_asize,maxd_atype)
    real(r8) :: adrydpav_box(maxd_asize,maxd_atype)
    real(r8) :: adryqmas_box(maxd_asize,maxd_atype)
    real(r8) :: awetdens_box(maxd_asize,maxd_atype)
    real(r8) :: awetdpav_box(maxd_asize,maxd_atype)
    
    
    
    
    
    
    real(r8) :: adrydens_box_sv(maxd_asize,maxd_atype,4)
    real(r8) :: adrydpav_box_sv(maxd_asize,maxd_atype,4)
    real(r8) :: adryqmas_box_sv(maxd_asize,maxd_atype,4)
    real(r8) :: awetdens_box_sv(maxd_asize,maxd_atype,4)
    real(r8) :: awetdpav_box_sv(maxd_asize,maxd_atype,4)


    if (idiagbb_host >= 100) then
       ldiag1 = .true. ; ldiag2 = .true.
       ldiag_movesect = .true. ; ldiag_newnuc = .true. ; ldiag_coag = .true. 
    else
       ldiag1 = .false. ; ldiag2 = .false.
       ldiag_movesect = .false. ; ldiag_newnuc = .false. ; ldiag_coag = .false. 
    end if

    jhyst_leg_sv2(1:nbin_a) = jhyst_leg(1:nbin_a)
    dp_dry_a_sv2( 1:nbin_a) = dp_dry_a( 1:nbin_a)
    sigmag_a_sv2( 1:nbin_a) = sigmag_a( 1:nbin_a)
    rbox_sv2(1:ntot_used) = rbox(1:ntot_used)
    fact_apdiam = fact_apdiam_in
    fact_apdens = fact_apdens_in


    
    
    
    if (mmovesect_flag1 <= 0) then
       write(*,*)   &
            '*** skipping sectional_interface_1 -- mmovesect_flag1 <= 0'
       return
    end if

    
    if (iunits_flagaa /= 4) then
       call wrf_error_fatal3("<stdin>",197,&
'*** sectional_interface_1 error -- iunits_flagaa /= 4')       
    end if
    if (iunits_flagbb /= 1) then
       call wrf_error_fatal3("<stdin>",201,&
'*** sectional_interface_1 error -- iunits_flagbb /= 1')       
    end if
    if (iunits_diagout_flagaa /= 1) then
       call wrf_error_fatal3("<stdin>",205,&
'*** sectional_interface_1 error -- iunits_diagout_flagaa /= 1')
    end if


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


    dens_nh4so4a_newnuc = dens_aer_mac(iso4_a)


    
    if ( ldiag2 ) call dump_secti( 'dym', &
       id, ii, jj, kk, ktau, ktauc, it_mosaic, dtchem, &
       rbox, cair_mol_m3 )
    if ( ldiag1 ) write(*,*) 'sectiface call sect_iface_bgn_end_calcs( 0 )'
    call sect_iface_bgn_end_calcs( 0, rbox, &
         dp_dry_a, sigmag_a, &
         drydens_aftgrow, drydens_pregrow,   &
         drymass_aftgrow, drymass_pregrow,   &
         adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
         adryqmas_box,   &
         fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
         it_mosaic, jaerosolstate_bgn, jaerosolstate, jhyst_leg, &
         mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a, &
         mw_aer_mac, cair_mol_cc, cair_mol_m3)

    jsv = 1
    rbox_svx(:,jsv) = rbox(:)
    adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
    awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
    adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
    awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
    adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


    
    
    
    method_movesect = mod( max(0,mmovesect_flag1), 100 )
    idiagbb_movesect = 0
    if ( ldiag_movesect .or. idiag_sect_movesect > 0) idiagbb_movesect = 70
    idiagbb_movesect = 70 
    if ( ii*jj*kk /= 1 ) idiagbb_movesect = 0
    if (method_movesect < 50) then
       if ( ldiag2 ) call dump_secti( 'dyn', &
          id, ii, jj, kk, ktau, ktauc, it_mosaic, dtchem, &
          rbox, cair_mol_m3 )
       if ( ldiag1 ) write(*,*) 'sectiface call move_sections_x3'







       idiagbb_movesect = 0
       if ( ldiag_movesect .and. ii*jj*kk == 1 ) idiagbb_movesect = 70
       call move_sections_x3( ai_phase, ii, jj, kk, 1, rbox,   &
            fact_apmassmr, fact_apnumbmr,   &
            fact_apdens, fact_apdiam,   &
            drydens_aftgrow, drydens_pregrow,   &
            drymass_aftgrow, drymass_pregrow,   &
            adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
            adryqmas_box, it_mosaic, mmovesect_flag1, idiagbb_movesect )
    else
       if ( ldiag2 ) call dump_secti( 'dyo', &
          id, ii, jj, kk, ktau, ktauc, it_mosaic, dtchem, &
          rbox, cair_mol_m3 )
       if ( ldiag1 ) write(*,*) 'sectiface call move_sect_3d_x1'
       call move_sect_3d_x1( 1, 1, 1, 1, 1, rbox,   &
            fact_apmassmr, fact_apnumbmr,   &
            fact_apdens, fact_apdiam,   &
            drydens_aftgrow, drydens_pregrow,   &
            drymass_aftgrow, drymass_pregrow,   &
            adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
            adryqmas_box, it_mosaic, mmovesect_flag1, idiagbb_movesect )
    end if

    jsv = 2
    rbox_svx(:,jsv) = rbox(:)
    adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
    awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
    adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
    awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
    adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


    
    
    if (mnewnuc_flag1 > 0) then
       idiagbb_newnuc = 0
       if ( ldiag_newnuc .or. idiag_sect_newnuc > 0) idiagbb_newnuc = 200
       if ( ii*jj*kk /= 1 ) idiagbb_newnuc = 0

       if (ntype_aer == 1) then
          itype_newnuc = 1
       else
          
          itype_newnuc = 0
          do it2 = 1, ntype_md2_aer
             if (hygro_so4_aer <= xcut_atype_md2(it2)) then
                itype_newnuc = itype_of_itype_md1md2(1,it2)
                exit
             end if
          end do
          if (itype_newnuc == 0) itype_newnuc = itype_of_itype_md1md2(1,ntype_md2_aer)
       end if

       rh_box = rh*0.01
       if ( ldiag2 ) &
       call dump_secti( 'dyp', id, ii, jj, kk, ktau, ktauc, it_mosaic, dtchem, &
          rbox, cair_mol_m3 )
       if ( ldiag1 ) write(*,*) 'sectiface call mosaic_newnuc_1box'
       call mosaic_newnuc_1box(   &
            istat_newnuc, idiagbb_newnuc,   &
            it_mosaic, itype_newnuc, dtchem,   &
            te, pr_atm, cair_mol_cc, rh_box,   &
            dens_nh4so4a_newnuc, mw_aer_mac(iso4_a), mw_aer_mac(inh4_a),   &
            mw_comp_a(jh2o), mw_air,  &
            fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr,   &
            rbox_sv1, rbox,   &
            adrydens_box, awetdens_box,   &
            adrydpav_box, awetdpav_box, adryqmas_box )
    end if

    jsv = 3
    rbox_svx(:,jsv) = rbox(:)
    adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
    awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
    adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
    awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
    adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


    
    
    if (mcoag_flag1 > 0) then
       idiagbb_coag = 0
       if ( ldiag_coag .or. idiag_sect_coag > 0) idiagbb_coag = 200
       if ( ii*jj*kk /= 1 ) idiagbb_coag = 0

       rhoair_g_cc = cair_mol_cc*mw_air
       method_coag = mod( max(0,mcoag_flag1), 100 )
       if (method_coag < 50) then
          if ( ldiag2 ) call dump_secti( 'dyq', &
             id, ii, jj, kk, ktau, ktauc, it_mosaic, dtchem, &
             rbox, cair_mol_m3 )
          if ( ldiag1 ) write(*,*) 'sectiface call mosaic_coag_1box'
          call mosaic_coag_1box( istat_coag,   &
               idiagbb_coag, it_mosaic,   &
               dtchem, te, pr_atm, rhoair_g_cc, rbox,   &
               fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam,   &
               adrydens_box, awetdens_box,   &
               adrydpav_box, awetdpav_box, adryqmas_box, mcoag_flag1 )
       else
          if ( ldiag2 ) call dump_secti( 'dyr', &
             id, ii, jj, kk, ktau, ktauc, it_mosaic, dtchem, &
             rbox, cair_mol_m3 )
          if ( ldiag1 ) write(*,*) 'sectiface call mosaic_coag_3d_1box'
          call mosaic_coag_3d_1box( istat_coag,   &
               idiagbb_coag, it_mosaic,   &
               dtchem, te, pr_atm, rhoair_g_cc, rbox,   &
               fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam,   &
               adrydens_box, awetdens_box,   &
               adrydpav_box, awetdpav_box, adryqmas_box, mcoag_flag1 )
       end if
    end if

    jsv = 4
    rbox_svx(:,jsv) = rbox(:)
    adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
    awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
    adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
    awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
    adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


    
    if ( ldiag2 ) call dump_secti( 'dys', &
       id, ii, jj, kk, ktau, ktauc, it_mosaic, dtchem, &
       rbox, cair_mol_m3 )
    if ( ldiag1 ) write(*,*) 'sectiface call sect_iface_bgn_end_calcs( 1 )'
    call sect_iface_bgn_end_calcs( 1, rbox, &
         dp_dry_a, sigmag_a, &
         drydens_aftgrow, drydens_pregrow,   &
         drymass_aftgrow, drymass_pregrow,   &
         adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
         adryqmas_box,   &
         fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
         it_mosaic, jaerosolstate_bgn, jaerosolstate, jhyst_leg, &
         mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a, &
         mw_aer_mac, cair_mol_cc, cair_mol_m3)


    if ( ldiag2 ) call dump_secti( 'dyt', &
       id, ii, jj, kk, ktau, ktauc, it_mosaic, dtchem, &
       rbox, cair_mol_m3 )
    if ( ldiag1 ) write(*,*) 'sectiface done'
    return
  end subroutine sectional_interface_1


  
  subroutine sect_iface_bgn_end_calcs( ibgn_end, rbox, &
       dp_dry_a, sigmag_a, &
       drydens_aftgrow, drydens_pregrow,   &
       drymass_aftgrow, drymass_pregrow,   &
       adrydens_box, awetdens_box, adrydpav_box,   &
       awetdpav_box, adryqmas_box,   &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       it_mosaic, jaerosolstate_bgn, jaerosolstate, jhyst_leg, &
       mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a, &
       mw_aer_mac, cair_mol_cc, cair_mol_m3)
    
    
    
    use module_data_mosaic_boxmod, only: idiag_sect_movesect
    use module_data_mosaic_main, only: &
       avogad, &
       mw_air, naer_tot, ngas_max, ntot_max, ntot_used
    use module_data_mosaic_aero, only: &
       jhyst_lo, jhyst_up, &
       mhyst_method, mhyst_force_lo, mhyst_force_up, &
       mhyst_uporlo_jhyst, mhyst_uporlo_waterhyst, &
       mmovesect_flag1, naer, naercomp, nbin_a, nbin_a_max, &
       all_solid, all_liquid, mixed, no_aerosol
    use module_data_mosaic_asecthp, only: &
       ai_phase, isize_of_ibin, itype_of_ibin, lunerr, &
       maxd_acomp, maxd_asize, maxd_atype, ncomp_aer, &
       hyswptr_aer, massptr_aer, &
       dens_aer, dens_water_aer, hygro_aer 

    use module_mosaic_movesect1d, only:  test_move_sections

    use module_peg_util, only:   peg_error_fatal

    implicit none

    
    integer, intent(in) :: ibgn_end, it_mosaic
    integer, intent(in), dimension(nbin_a_max) :: jaerosolstate, jaerosolstate_bgn
    integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg

    real(r8), intent(in) :: cair_mol_cc, cair_mol_m3
    real(r8), intent(inout), dimension(nbin_a_max)  :: mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a
    real(r8), intent(inout), dimension(nbin_a_max)  :: dp_dry_a, sigmag_a
    real(r8), intent(inout), dimension(naer) :: mw_aer_mac
    real(r8), intent(inout) :: rbox(ntot_used)

    real(r8), intent(inout) :: drydens_aftgrow(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: drydens_pregrow(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: drymass_aftgrow(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: drymass_pregrow(maxd_asize,maxd_atype)

    real(r8), intent(inout) :: adrydens_box(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: awetdens_box(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: adrydpav_box(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: awetdpav_box(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: adryqmas_box(maxd_asize,maxd_atype)

    real(r8), intent(in) :: fact_apmassmr, fact_apnumbmr, &
                            fact_apdens, fact_apdiam, fact_gasmr


    
    integer :: ibin, iphase, isize, itype, jhyst_tmp
    integer :: l, ll, mtmp
    real(r8) :: tmpa, tmph, tmpj, tmpr
    character(len=200) :: msg



    
    
    if ((it_mosaic == 1) .and. (ibgn_end == 0)) then
       call test_move_sections( 1, 1, 1, 1, 1, it_mosaic, mmovesect_flag1, idiag_sect_movesect )
    end if


    if ((ibgn_end < 0) .or. (ibgn_end > 1)) return

    if (ibgn_end == 1) goto 20000


    
    
    
    
    
    
    adrydens_box(:,:) = 0.0
    awetdens_box(:,:) = 0.0
    adrydpav_box(:,:) = 0.0
    awetdpav_box(:,:) = 0.0
    adryqmas_box(:,:) = 0.0

    do ibin = 1, nbin_a
       isize = isize_of_ibin(ibin)
       itype = itype_of_ibin(ibin)

       if (jaerosolstate_bgn(ibin) /= no_aerosol) then
          drydens_pregrow(isize,itype) = dens_dry_a_bgn(ibin)
          drymass_pregrow(isize,itype) = mass_dry_a_bgn(ibin)
       else
          drydens_pregrow(isize,itype) = -1.0
          drymass_pregrow(isize,itype) =  0.0
       end if
       if (jaerosolstate(ibin) /= no_aerosol) then
          drydens_aftgrow(isize,itype) = dens_dry_a(ibin)
          drymass_aftgrow(isize,itype) = mass_dry_a(ibin)
       else
          drydens_aftgrow(isize,itype) = -1.0
          drymass_aftgrow(isize,itype) =  0.0
       end if

       
       
       tmpa = 1.0e6_r8/(cair_mol_m3*mw_air*fact_apmassmr)
       drymass_pregrow(isize,itype) = drymass_pregrow(isize,itype) * tmpa
       drymass_aftgrow(isize,itype) = drymass_aftgrow(isize,itype) * tmpa

       
       
       tmpa = 1.0_r8/fact_apdens
       drydens_pregrow(isize,itype) = drydens_pregrow(isize,itype) * tmpa
       drydens_aftgrow(isize,itype) = drydens_aftgrow(isize,itype) * tmpa

       if (mhyst_method == mhyst_uporlo_jhyst) then
          
          msg = '*** sect_iface_bgn_end_calcs - mhyst_uporlo_jhyst not allowed in wrf-chem'
          call peg_error_fatal( lunerr, msg )
       end if

    end do

    return


20000 continue
    
    
    
    
    
    
    
    do ibin = 1, nbin_a
       isize = isize_of_ibin(ibin)
       itype = itype_of_ibin(ibin)

       
       dp_dry_a(ibin) = adrydpav_box(isize,itype)

       if (mhyst_method == mhyst_uporlo_waterhyst) then
          continue
       else if (mhyst_method == mhyst_force_up) then
          jhyst_leg(ibin) = jhyst_up
       else if (mhyst_method == mhyst_force_lo) then
          jhyst_leg(ibin) = jhyst_lo
       else if (mhyst_method == mhyst_uporlo_jhyst) then
          
          msg = '*** sect_iface_bgn_end_calcs - mhyst_uporlo_jhyst not allowed in wrf-chem'
          call peg_error_fatal( lunerr, msg )
       else
          write(msg,'(a,i10)') &
               '*** sect_iface_bgn_end_calcs - bad mhyst_method =', mhyst_method
          call peg_error_fatal( lunerr, msg )
       end if

    end do

    return
  end subroutine sect_iface_bgn_end_calcs


  
  subroutine dump_secti( dumpchaa, &
       id, ii, jj, kk, ktau, ktauc, it_mosaic, dtchem, &
       rbox, cair_mol_m3 )

    use module_data_mosaic_main, only: mw_air, ntot_max, ntot_used
    use module_data_mosaic_aero, only: &
       inh4_a, iso4_a, jh2o, &
       dens_aer_mac, mw_aer_mac, mw_comp_a, &
       mcoag_flag1, mmovesect_flag1, mnewnuc_flag1, &
       naer, naercomp, nbin_a, nbin_a_max
    use module_data_mosaic_asecthp, only:  maxd_asize, maxd_atype, &
       nsize_aer, ntype_aer, &
       numptr_aer, lptr_bc_aer, lptr_oc_aer, lptr_oin_aer

    implicit none

    
    character(len=*), intent(in) :: dumpchaa

    integer, intent(in) :: id, ii, jj, kk, ktau, ktauc, it_mosaic

    real(r8), intent(in)    :: dtchem  
    real(r8), intent(in)    :: cair_mol_m3  
    real(r8), intent(inout) :: rbox(ntot_used)  

    
    integer :: isize, itype
    real(r8) :: tmpa
    character(len=80) :: fmtaa

    if (ii*jj*kk /= 1) return

    fmtaa = '(2a,1p,e10.3,7e11.3)'
    if (nsize_aer(1)*ntype_aer > 8) fmtaa = '(2a,1p,e10.3,7e11.3/(7x,1p,e10.3,7e11.3))'


    tmpa = (cair_mol_m3*mw_air*1.0e-3_r8)*1.0e-6_r8
    write(131,fmtaa) dumpchaa, ' num', &
       ( ( rbox(numptr_aer(isize,itype,1))*tmpa, &
           isize=1,nsize_aer(itype) ), itype=1,ntype_aer )



    tmpa = (cair_mol_m3*mw_air*1.0e-3_r8)*1.0e3_r8
    write(131,fmtaa) dumpchaa, ' oin', &
       ( ( rbox(lptr_oin_aer(isize,itype,1))*tmpa, &
           isize=1,nsize_aer(itype) ), itype=1,ntype_aer )
    write(131,fmtaa) dumpchaa, ' oc ', &
       ( ( rbox(lptr_oc_aer(isize,itype,1))*tmpa, &
           isize=1,nsize_aer(itype) ), itype=1,ntype_aer )
    write(131,fmtaa) dumpchaa, ' bc ', &
       ( ( rbox(lptr_bc_aer(isize,itype,1))*tmpa, &
           isize=1,nsize_aer(itype) ), itype=1,ntype_aer )

    return
  end subroutine dump_secti


  


  end module module_mosaic_sect_intr

