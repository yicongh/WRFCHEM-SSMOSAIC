	module module_mosaic_movesect3d


	use module_data_mosaic_kind, only:  r8
	use module_data_mosaic_main, only:  mw_air, ntot_used, pi
	use module_data_mosaic_asecthp, only:  lunerr, lunout
	use module_peg_util


	implicit none


	contains










	subroutine move_sect_3d_x1( iphase_flag, iclm, jclm, k, m, rbox,   &
          fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam,   &
      	  drydens_aftgrow, drydens_pregrow,   &
      	  drymass_aftgrow, drymass_pregrow,   &
      	  adrydens_tmp, awetdens_tmp, adrydpav_tmp, awetdpav_tmp,   &
      	  adryqmas_tmp, it_mosaic, mmovesect_flag1, idiag_movesect )


























































	use module_data_mosaic_asecthp, only:  &
	    maxd_acomp, maxd_asize, maxd_atype,   &
	    ntype_aer, nsize_aer, nphase_aer, ai_phase, cw_phase,   &
	    dens_aer, dcen_sect,   &
	    smallmassaa_asecthp => smallmassaa, smallmassbb_asecthp => smallmassbb,   &
	    volumcen_sect, volumlo_sect, volumhi_sect

	implicit none


	integer, intent(in) :: iphase_flag, iclm, jclm, k, m, it_mosaic, &
	                       mmovesect_flag1, idiag_movesect

	real(r8), intent(inout) :: rbox(ntot_used)

	real(r8), intent(in) :: fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam

	real(r8), intent(in) :: drymass_pregrow(maxd_asize,maxd_atype)
        real(r8), intent(in) :: drydens_pregrow(maxd_asize,maxd_atype)
        real(r8), intent(in) :: drymass_aftgrow(maxd_asize,maxd_atype)
        real(r8), intent(in) :: drydens_aftgrow(maxd_asize,maxd_atype)

	real(r8), intent(inout) :: adrydens_tmp(maxd_asize,maxd_atype), awetdens_tmp(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: adrydpav_tmp(maxd_asize,maxd_atype), awetdpav_tmp(maxd_asize,maxd_atype)
	real(r8), intent(inout) :: adryqmas_tmp(maxd_asize,maxd_atype)







	integer iphase, isize, itype,   &
	  l, ll, llhysw, llwater, lnew, lold, l3,   &
      	  method_movesect, n, nnew, nold
	integer nnewsave(2,maxd_asize,maxd_atype)
	integer itypenewsave(4,maxd_asize,maxd_atype)

	real(r8) densh2o
	real(r8) delta_water_conform1, delta_numb_conform1
	real(r8) fact_apvolu
	real(r8) smallmassaa, smallmassbb

	real(r8) dcen_stmp(maxd_asize,maxd_atype)
	real(r8) densdefault(maxd_atype)
	real(r8) drydenspp(maxd_asize,maxd_atype), drydensxx0(maxd_asize,maxd_atype),   &
      	    drydensxx(maxd_asize,maxd_atype), drydensyy(maxd_asize,maxd_atype)
	real(r8) drymasspp(maxd_asize,maxd_atype), drymassxx0(maxd_asize,maxd_atype),   &
      	    drymassxx(maxd_asize,maxd_atype), drymassyy(maxd_asize,maxd_atype)
	real(r8) dryvolxx(maxd_asize,maxd_atype), dryvolyy(maxd_asize,maxd_atype)
	real(r8) rmassxx(maxd_acomp+2,maxd_asize,maxd_atype),   &
      	    rmassyy(maxd_acomp+2,maxd_asize,maxd_atype)
	real(r8) rnumbpp(maxd_asize,maxd_atype), rnumbxx0(maxd_asize,maxd_atype),   &
      	    rnumbxx(maxd_asize,maxd_atype), rnumbyy(maxd_asize,maxd_atype)
	real(r8) specdensxx(maxd_acomp)
	real(r8) xferfracvol(2,4,maxd_asize,maxd_atype), xferfracnum(2,4,maxd_asize,maxd_atype)
	real(r8) volumcen_stmp(maxd_asize,maxd_atype), &
	    volumlo_stmp(maxd_asize,maxd_atype), volumhi_stmp(maxd_asize,maxd_atype)
	real(r8) wetvolxx(maxd_asize,maxd_atype), wetvolyy(maxd_asize,maxd_atype)
	real(r8) wetmassxx(maxd_asize,maxd_atype), wetmassyy(maxd_asize,maxd_atype)

	character*160 msg





	if (mmovesect_flag1 <= 0) return
	if (ntype_aer       <= 0) return
	if (nphase_aer      <= 0) return



	method_movesect = mod( max(0,mmovesect_flag1), 100 )

	if ((method_movesect/10 .eq. 5) .or.   &
	    (method_movesect/10 .eq. 6)) then
	    continue
	else
	    msg = '*** subr move_sect_3d error - ' //   &
		'illegal value for mmovesect_flag1'
	    call peg_error_fatal( lunerr, msg )
	end if



	if (iabs(iphase_flag) .eq. 1) then
	    iphase = ai_phase
	else if (iabs(iphase_flag) .eq. 2) then










	    msg = '*** subr move_sect_3d error - ' //   &
		'iphase_flag=2 (after aqueous chemistry) is not implemented'
	    call peg_error_fatal( lunerr, msg )
	else
	    msg = '*** subr move_sect_3d error - ' //   &
		'iabs(iphase_flag) must be 1 or 2'
	    call peg_error_fatal( lunerr, msg )
	end if




	if (iphase_flag .le. 0) then
	    write(msg,9040) 'method, idiag', method_movesect, idiag_movesect
	    call peg_message( lunout, msg )
	    call move_sect_3d_checkptrs( iphase_flag, iclm, jclm, k, m )
	    return
	else if (it_mosaic .le. 1) then


	    call move_sect_3d_checkptrs( -iphase_flag, iclm, jclm, k, m )
	end if
9040	format( '*** subr move_sect_3d - ', a, ' =', 2(1x,i6) )



	if (idiag_movesect .ge. 70) then
	    msg = ' '
	    call peg_message( lunout, msg )
	    write(msg,9060) mmovesect_flag1, iclm, jclm, k, m, it_mosaic
	    call peg_message( lunout, msg )
	end if
9060	format( '*** move_sect_3d diags - ', &
	    'msflag, ijkm, itime =', i7, 4i4, i7 )


	densh2o = 1.0
	smallmassaa = smallmassaa_asecthp
	smallmassbb = smallmassbb_asecthp











	fact_apvolu = fact_apdiam*fact_apdiam*fact_apdiam


	do itype = 1, ntype_aer
	densdefault(itype) = dens_aer(1,itype)*fact_apdens   
	do isize = 1, nsize_aer(itype)
	    dcen_stmp(    isize,itype) = dcen_sect(    isize,itype)*fact_apdiam
	    volumcen_stmp(isize,itype) = volumcen_sect(isize,itype)*fact_apvolu
	    volumlo_stmp( isize,itype) = volumlo_sect( isize,itype)*fact_apvolu
	    volumhi_stmp( isize,itype) = volumhi_sect( isize,itype)*fact_apvolu
	end do
	end do



	do itype = 1, ntype_aer
	do isize = 1, nsize_aer(itype)
	    drydenspp(isize,itype) = drydens_pregrow(isize,itype)*fact_apdens
	    drydensxx(isize,itype) = drydens_aftgrow(isize,itype)*fact_apdens

	    drymasspp(isize,itype) = drymass_pregrow(isize,itype)*fact_apmassmr
	    drymassxx(isize,itype) = drymass_aftgrow(isize,itype)*fact_apmassmr
	end do
	end do


	call move_sect_3d_initial_conform(   &
	  iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
          fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_apvolu,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  rbox,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drydensxx0, drymassxx0, rnumbxx0,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx, specdensxx,   &
	  volumcen_stmp, volumlo_stmp, volumhi_stmp )

	if (method_movesect/10 .eq. 5) then

	call move_sect_3d_calc_movingcenter(   &
	  iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
	  nnewsave, itypenewsave,    &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx,   &
      	  xferfracvol, xferfracnum,   &
	  volumcen_stmp, volumlo_stmp, volumhi_stmp )
	else
	call move_sect_3d_calc_masnumadvect(   &
	  iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
	  nnewsave, itypenewsave,    &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drydensxx0, drymassxx0, rnumbxx0,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx,   &
      	  xferfracvol, xferfracnum,   &
	  volumcen_stmp, volumlo_stmp, volumhi_stmp )
	end if


	call move_sect_3d_apply_moves(   &
	  iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
	  nnewsave, itypenewsave,    &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
          fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_apvolu,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  rbox,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy,   &
      	  xferfracvol, xferfracnum,   &
	  dcen_stmp, volumcen_stmp, volumlo_stmp, volumhi_stmp,   &
      	  adrydens_tmp, awetdens_tmp, adrydpav_tmp, awetdpav_tmp,   &
      	  adryqmas_tmp )


	return
	end subroutine move_sect_3d_x1



	subroutine move_sect_3d_initial_conform(   &
	  iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
          fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_apvolu,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  rbox,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drydensxx0, drymassxx0, rnumbxx0,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx, specdensxx,   &
	  volumcen_stmp, volumlo_stmp, volumhi_stmp )











	use module_data_mosaic_asecthp, only:  &
	    maxd_acomp, maxd_asize, maxd_atype,   &
	    nsize_aer, ntype_aer, ncomp_aer, ncomp_plustracer_aer, ai_phase,   &
	    massptr_aer, numptr_aer, waterptr_aer, hyswptr_aer,   &
	    dens_aer

	implicit none


	integer iphase_flag, iclm, jclm, iphase, k,   &
      	  m, method_movesect, idiag_movesect, llhysw, llwater
	real(r8) densh2o, smallmassaa, smallmassbb
	real(r8) fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_apvolu
	real(r8) delta_water_conform1, delta_numb_conform1
	real(r8) densdefault(maxd_atype)
	real(r8) drydenspp(maxd_asize,maxd_atype), drydensxx0(maxd_asize,maxd_atype),   &
      	     drydensxx(maxd_asize,maxd_atype)
	real(r8) drymasspp(maxd_asize,maxd_atype), drymassxx0(maxd_asize,maxd_atype),   &
      	     drymassxx(maxd_asize,maxd_atype)
	real(r8) dryvolxx(maxd_asize,maxd_atype)
	real(r8) rbox(ntot_used)
	real(r8) rmassxx(maxd_acomp+2,maxd_asize,maxd_atype)
	real(r8) rnumbpp(maxd_asize,maxd_atype), rnumbxx0(maxd_asize,maxd_atype),   &
      	     rnumbxx(maxd_asize,maxd_atype)
	real(r8) specdensxx(maxd_acomp)
	real(r8) volumcen_stmp(maxd_asize,maxd_atype), &
	    volumlo_stmp(maxd_asize,maxd_atype), volumhi_stmp(maxd_asize,maxd_atype)
	real(r8) wetvolxx(maxd_asize,maxd_atype)
	real(r8) wetmassxx(maxd_asize,maxd_atype)



	integer itype, l, ll, lnew, lold, l3, n, nnew, nold

	real(r8) dummass, dumnum, dumnum_at_dhi, dumnum_at_dlo, dumr,   &
      	  dumvol, dumvol1p, dumwatrmass





	itype = 1
	do ll = 1, ncomp_plustracer_aer(itype)
	    specdensxx(ll) = dens_aer(ll,itype)*fact_apdens
	end do


	do l = 1, ntot_used
	    rbox(l) = max( 0.0_r8, rbox(l) )
	end do


	delta_water_conform1 = 0.0
	delta_numb_conform1 = 0.0




	do 1900 itype = 1, ntype_aer


	llhysw = ncomp_plustracer_aer(itype) + 1
	llwater = ncomp_plustracer_aer(itype) + 2
	do n = 1, nsize_aer(itype)
	    do ll = 1, ncomp_plustracer_aer(itype)
		l = massptr_aer(ll,n,itype,iphase)
		rmassxx(ll,n,itype) = rbox(l)*fact_apmassmr
	    end do
	    rmassxx(llhysw,n,itype) = 0.
	    l = 0
	    if (iphase .eq. ai_phase) l = hyswptr_aer(n,itype)
	    if (l .gt. 0) rmassxx(llhysw,n,itype) = rbox(l)*fact_apmassmr
	    rmassxx(llwater,n,itype) = 0.
	    l = 0
	    if (iphase .eq. ai_phase) l = waterptr_aer(n,itype)
	    if (l .gt. 0) rmassxx(llwater,n,itype) = rbox(l)*fact_apmassmr

	    rnumbxx(n,itype)  = rbox(numptr_aer(n,itype,iphase))*fact_apnumbmr
	    rnumbxx0(n,itype) = rnumbxx(n,itype)
	    rnumbpp(n,itype)  = rnumbxx(n,itype)

	    drydensxx0(n,itype) = drydensxx(n,itype)
	    drymassxx0(n,itype) = drymassxx(n,itype)
	end do

	do 1390 n = 1, nsize_aer(itype)








	if ( (drydensxx(n,itype) .lt.  0.1) .or.   &
	     (drydensxx(n,itype) .gt. 20.0) .or.   &
      	     (drymassxx(n,itype) .le. smallmassaa) ) then
	    dummass = 0.
	    dumvol = 0.
	    do ll = 1, ncomp_aer(itype)
		dumr = rmassxx(ll,n,itype)
		dummass = dummass + dumr
		dumvol  = dumvol  + dumr/specdensxx(ll)
	    end do
	    drymassxx(n,itype) = dummass
	    if (min(dummass,dumvol) .le. smallmassbb) then
		drydensxx(n,itype) = densdefault(itype)
		dumvol = dummass/densdefault(itype)
		dumnum = dummass/(volumcen_stmp(n,itype)*densdefault(itype))
	    else
		drydensxx(n,itype) = dummass/dumvol
		dumnum = rnumbxx(n,itype)
		dumnum_at_dhi = dumvol/volumhi_stmp(n,itype)
		dumnum_at_dlo = dumvol/volumlo_stmp(n,itype)
		dumnum = max( dumnum_at_dhi, min( dumnum_at_dlo, dumnum ) )
	    end if
	    delta_numb_conform1 = delta_numb_conform1 + dumnum - rnumbxx(n,itype)
	    rnumbxx(n,itype) = dumnum
	    rnumbpp(n,itype) = rnumbxx(n,itype)
	    delta_water_conform1 = delta_water_conform1 - rmassxx(llwater,n,itype) 
	    rmassxx(llwater,n,itype) = 0.
	end if



	dryvolxx(n,itype) = drymassxx(n,itype)/drydensxx(n,itype)
	dumwatrmass = rmassxx(llwater,n,itype)
	wetmassxx(n,itype) = drymassxx(n,itype) + dumwatrmass
	wetvolxx(n,itype) = dryvolxx(n,itype) + dumwatrmass/densh2o

1390	continue

1900	continue

	return
	end subroutine move_sect_3d_initial_conform                          



	subroutine move_sect_3d_calc_movingcenter(   &
	  iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
	  nnewsave, itypenewsave,    &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx,   &
      	  xferfracvol, xferfracnum,   &
	  volumcen_stmp, volumlo_stmp, volumhi_stmp )








	use module_data_mosaic_aero, only:  &
	    method_bcfrac, method_kappa, msectional_flag2
	use module_data_mosaic_asecthp, only:  &
	    dens_aer, hygro_aer,   &
	    itype_md1_of_itype, itype_md2_of_itype, itype_of_itype_md1md2,   &
	    lptr_bc_aer, massptr_aer, &
	    maxd_acomp, maxd_asize, maxd_atype,   &
	    ncomp_aer, nsize_aer, ntype_aer, ntype_md1_aer, ntype_md2_aer,   &
	    xcut_atype_md1, xcut_atype_md2

	implicit none


	integer iphase_flag, iclm, jclm, iphase, k,   &
      	  m, method_movesect, idiag_movesect, llhysw, llwater
	integer nnewsave(2,maxd_asize,maxd_atype)
	integer itypenewsave(4,maxd_asize,maxd_atype)
	real(r8) densh2o, smallmassaa, smallmassbb
	real(r8) densdefault(maxd_atype)
	real(r8) drydensxx(maxd_asize,maxd_atype)
	real(r8) drymassxx(maxd_asize,maxd_atype)
	real(r8) dryvolxx(maxd_asize,maxd_atype)
	real(r8) rmassxx(maxd_acomp+2,maxd_asize,maxd_atype)
	real(r8) rnumbxx(maxd_asize,maxd_atype)
	real(r8) xferfracvol(2,4,maxd_asize,maxd_atype), xferfracnum(2,4,maxd_asize,maxd_atype)
	real(r8) volumcen_stmp(maxd_asize,maxd_atype), &
	    volumlo_stmp(maxd_asize,maxd_atype), volumhi_stmp(maxd_asize,maxd_atype)
	real(r8) wetmassxx(maxd_asize,maxd_atype)
	real(r8) wetvolxx(maxd_asize,maxd_atype)


	integer isize, itype, itype1, itype2, itypenew, itype1new, itype2new
	integer itmpa, itmpb
	integer jtype1, jtype2
	integer ll, n, ndum, ndumb, nnew, nold
	real(r8) dumnum, dumvol, dumvol1p, sixoverpi, third
	real(r8) tmpa, tmpb, tmpc, tmpd, tmpe, tmpf
	real(r8) tmp_kappa, tmp_wbc
	character*160 fmtaa, msg
	character*11 txt11


	sixoverpi = 6.0/pi
	third = 1.0/3.0

	if (idiag_movesect .ge. 70) then
	    call peg_message( lunout, ' ' )
	    call peg_message( lunout, &
		'move_sect_3d_calc_movingcenter diagnostics' )
	    call peg_message( lunout, ' ' )
	end if




	do 1900 itype = 1, ntype_aer





	do 1390 n = 1, nsize_aer(itype)

	nnew = n
	itypenew = itype


	if (drymassxx(n,itype) .le. smallmassaa) goto 1290

	dumvol = dryvolxx(n,itype)
	dumnum = rnumbxx(n,itype)



	isize = nsize_aer(itype)
	if ( dumnum .le. dumvol/volumhi_stmp(isize,itype) ) then
	    nnew = nsize_aer(itype)
	    goto 1250

	else if ( dumnum .ge. dumvol/volumlo_stmp(1,itype) ) then
	    nnew = 1
	    goto 1250
	end if


	dumvol1p = dumvol/dumnum
	if ( dumvol1p .gt. volumhi_stmp(n,itype) ) then
	    do while ( ( nnew .lt. nsize_aer(itype) ) .and.   &
      		       ( dumvol1p .gt. volumhi_stmp(nnew,itype) ) )
		nnew = nnew + 1
	    end do

	else if ( dumvol1p .lt. volumlo_stmp(n,itype) ) then
	    do while ( ( nnew .gt. 1 ) .and.   &
      		       ( dumvol1p .lt. volumlo_stmp(nnew,itype) ) )
		nnew = nnew - 1
	    end do

	end if

1250	continue
	if ((ntype_aer <= 1) .or. (msectional_flag2 <=0)) then
	    goto 1290   
	end if


	tmpa = 0.0   
	tmpb = 0.0   
	tmpc = 0.0   
	tmpd = 0.0   
	do ll = 1, ncomp_aer(itype)
	    tmpe = rmassxx(ll,n,itype)
	    tmpb = tmpb + tmpe


	    if (massptr_aer(ll,n,itype,iphase) == lptr_bc_aer(n,itype,iphase)) then
		tmpa = tmpa + tmpe
		if (method_kappa == 12) cycle
		
	    end if
	    tmpf = tmpe/dens_aer(ll,itype)
	    tmpc = tmpc + tmpf*hygro_aer(ll,itype)
	    tmpd = tmpd + tmpf
	end do
	tmp_wbc   = tmpa/max(tmpb,1.0e-35_r8)
	tmp_kappa = tmpc/max(tmpd,1.0e-35_r8)

	itype1new = ntype_md1_aer
	do jtype1 = 1, ntype_md1_aer - 1
	    if (tmp_wbc <= xcut_atype_md1(jtype1)) then
		itype1new = jtype1
		exit
	    end if
	end do

	itype2new = ntype_md2_aer
	do jtype2 = 1, ntype_md2_aer - 1
	    if (tmp_kappa <= xcut_atype_md2(jtype2)) then
		itype2new = jtype2
		exit
	    end if
	end do

	itypenew = itype_of_itype_md1md2(itype1new,itype2new)
	itype1 = itype_md1_of_itype(itype)
	itype2 = itype_md2_of_itype(itype)

	if (idiag_movesect .ge. 70) then
	    write(msg,'(a,3i5,1p,2e10.2)') 'isz,ity,itynew,wbc,kap', n, itype, itypenew, tmp_wbc, tmp_kappa
	    call peg_message( lunout, msg )
	end if










1290	continue

	nnewsave(:,n,itype) = 0
	nnewsave(1,n,itype) = nnew
	itypenewsave(:,n,itype) = 0
	itypenewsave(1,n,itype) = itypenew
	if (method_movesect == 51) itypenewsave(1,n,itype) = itype

	xferfracvol(:,:,n,itype) = 0.0
	xferfracvol(1,1,n,itype) = 1.0
	xferfracnum(:,:,n,itype) = 0.0
	xferfracnum(1,1,n,itype) = 1.0

1390	continue



	if (idiag_movesect .ge. 70) then
	    ndum = 0 ; ndumb = 0
	    do n = 1, nsize_aer(itype)
		if (nnewsave(1,n,itype) .ne. n) ndum = ndum + 1
		if (itypenewsave(1,n,itype) .ne. itype) ndumb = ndumb + 1
	    end do
	    if (ndum .gt. 0 .and. ndumb .gt. 0) then
		txt11 = 'movesectYS3'
	    else if (ndum .gt. 0) then
		txt11 = 'movesectYS1'
	    else if (ndumb .gt. 0) then
		txt11 = 'movesectYS2'
	    else
		txt11 = 'movesectNO '
	    end if
	    do itmpa = 1, nsize_aer(itype), 24
	        itmpb = min( itmpa+23, nsize_aer(itype) )
		if (itmpa == 1) then
		    fmtaa = '( a, 4i3, i5, i6, 3x, 24i3 )'
		    write(msg,fmtaa) txt11, iclm, jclm, k, m, itype,   &
	      		ndum, (nnewsave(1,n,itype), n=itmpa,itmpb)
		else
		    fmtaa = '( 37x, 24i3 )'
		    write(msg,fmtaa) (nnewsave(1,n,itype), n=itmpa,itmpb)
		end if
		call peg_message( lunout, msg )
	    end do

	    txt11 = 'itypenew   '
	    do itmpa = 1, nsize_aer(itype), 24
	        itmpb = min( itmpa+23, nsize_aer(itype) )
		if (itmpa == 1) then
		    fmtaa = '( 26x, a, 24i3 )'
		    write(msg,fmtaa) txt11, (itypenewsave(1,n,itype), n=itmpa,itmpb)
		else
		    fmtaa = '( 37x, 24i3 )'
		    write(msg,fmtaa) (itypenewsave(1,n,itype), n=itmpa,itmpb)
		end if
		call peg_message( lunout, msg )
	    end do

	end if


1900	continue

	return
	end subroutine move_sect_3d_calc_movingcenter                          



	subroutine move_sect_3d_calc_masnumadvect(   &
	  iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
	  nnewsave, itypenewsave,    &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drydensxx0, drymassxx0, rnumbxx0,   &
      	  drydensxx, drymassxx, dryvolxx, rmassxx, rnumbxx,   &
      	  wetmassxx, wetvolxx,   &
      	  xferfracvol, xferfracnum,   &
	  volumcen_stmp, volumlo_stmp, volumhi_stmp )













	use module_data_mosaic_aero, only:  &
	    method_bcfrac, method_kappa, msectional_flag2
	use module_data_mosaic_asecthp, only:  &
	    dens_aer, hygro_aer,   &
	    itype_md1_of_itype, itype_md2_of_itype, itype_of_itype_md1md2,   &
	    lptr_bc_aer, massptr_aer, &
	    maxd_acomp, maxd_asize, maxd_atype,   &
	    ncomp_aer, nsize_aer, ntype_aer, ntype_md1_aer, ntype_md2_aer,   &
	    xcut_atype_md1, xcut_atype_md2

	implicit none


	integer iphase_flag, iclm, jclm, iphase, k,   &
      	  m, method_movesect, idiag_movesect, llhysw, llwater
	integer nnewsave(2,maxd_asize,maxd_atype)
	integer itypenewsave(4,maxd_asize,maxd_atype)

	real(r8) densh2o, smallmassaa, smallmassbb
	real(r8) densdefault(maxd_atype)
	real(r8) drydenspp(maxd_asize,maxd_atype), drydensxx0(maxd_asize,maxd_atype),   &
      	     drydensxx(maxd_asize,maxd_atype)
	real(r8) drymasspp(maxd_asize,maxd_atype), drymassxx0(maxd_asize,maxd_atype),   &
      	     drymassxx(maxd_asize,maxd_atype)
	real(r8) dryvolxx(maxd_asize,maxd_atype)
	real(r8) rmassxx(maxd_acomp+2,maxd_asize,maxd_atype)
	real(r8) rnumbpp(maxd_asize,maxd_atype), rnumbxx0(maxd_asize,maxd_atype),   &
      	     rnumbxx(maxd_asize,maxd_atype)
	real(r8) xferfracvol(2,4,maxd_asize,maxd_atype), xferfracnum(2,4,maxd_asize,maxd_atype)
	real(r8) volumcen_stmp(maxd_asize,maxd_atype), &
	    volumlo_stmp(maxd_asize,maxd_atype), volumhi_stmp(maxd_asize,maxd_atype)
	real(r8) wetvolxx(maxd_asize,maxd_atype)
	real(r8) wetmassxx(maxd_asize,maxd_atype)


	integer ierr, itype, itype1, itype2, itypenew, itype1new, itype2new
	integer iforce_movecenter(maxd_asize,maxd_atype)
	integer jtype1, jtype2
	integer kktypenew
	integer ll
	integer n, nnew, nnew2

	real(r8) dum1, dum2, dum3
	real(r8) dumaa, dumbb, dumgamma, dumratio
	real(r8) dumfracnum, dumfracvol
	real(r8) dumntot
	real(r8) dumv
	real(r8) dumvbar_aft, dumvbar_pre
	real(r8) dumvcutlo_nnew_pre, dumvcuthi_nnew_pre
	real(r8) dumvlo_pre, dumvhi_pre, dumvdel_pre
	real(r8) dumvtot_aft, dumvtot_pre
	real(r8) dumzlo, dumzhi
	real(r8) sixoverpi, third, tmpa
	real(r8) tmpb, tmpc, tmpd, tmpe, tmpf, tmp_kappa, tmp_wbc

	character*4 dumch4
	character*1 dumch1
	character*160 msg


	sixoverpi = 6.0/pi
	third = 1.0/3.0
	itypenewsave(:,:,:) = 0











	do 4900 itype = 1, ntype_aer

	if (idiag_movesect .ge. 70) then
	    call peg_message( lunout, ' ' )
	    write( msg, '(a,i5)' ) &
		'move_sect_3d_calc_massnumadvect diagnostics - itype ', &
		itype
	    call peg_message( lunout, msg )
	end if

	do 3900 n = 1, nsize_aer(itype)

	nnew = n
	iforce_movecenter(n,itype) = 0

	kktypenew = 1
	itypenewsave(1,n,itype) = itype




	xferfracvol(:,:,n,itype) = 0.0
	xferfracnum(:,:,n,itype) = 0.0
	xferfracvol(1,kktypenew,n,itype) = 1.0
	xferfracnum(1,kktypenew,n,itype) = 1.0

	dumvtot_aft = -1.0
	dumvtot_pre = -1.0
	dumvbar_aft = -1.0
	dumvbar_pre = -1.0
	dumvlo_pre = -1.0
	dumvhi_pre = -1.0
	dumgamma = -1.0
	dumratio = -1.0
	dumvcutlo_nnew_pre = volumlo_stmp(nnew,itype)*(dumvbar_pre/dumvbar_aft)
	dumvcuthi_nnew_pre = volumhi_stmp(nnew,itype)*(dumvbar_pre/dumvbar_aft)
	dumfracvol = -1.0
	dumfracnum = -1.0


	if (drymassxx(n,itype) .le. smallmassaa) then
	    iforce_movecenter(n,itype) = 1
	    goto 1290
	end if

	dumvtot_aft = dryvolxx(n,itype)
	dumntot = rnumbxx(n,itype)



	if (dumntot .le. dumvtot_aft/volumhi_stmp(nsize_aer(itype),itype)) then
	    nnew = nsize_aer(itype)
	    iforce_movecenter(n,itype) = 2
	    goto 1290
	else if (dumntot .ge. dumvtot_aft/volumlo_stmp(1,itype)) then
	    nnew = 1
	    iforce_movecenter(n,itype) = 3
	    goto 1290
	end if



	dumvbar_aft = dumvtot_aft/dumntot
	if (dumvbar_aft .gt. volumhi_stmp(n,itype)) then
	    do while ( (nnew .lt. nsize_aer(itype)) .and.   &
      		       (dumvbar_aft .gt. volumhi_stmp(nnew,itype)) )
		nnew = nnew + 1
	    end do

	else if (dumvbar_aft .lt. volumlo_stmp(n,itype)) then
	    do while ( (nnew .gt. 1) .and.   &
      		       (dumvbar_aft .lt. volumlo_stmp(nnew,itype)) )
		nnew = nnew - 1
	    end do

	end if

1290	nnewsave(:,n,itype) = 0
	nnewsave(1,n,itype) = nnew

	if (iforce_movecenter(n,itype) .gt. 0) goto 3700








	if ( (drydenspp(n,itype) .lt.  0.1) .or.   &
	     (drydenspp(n,itype) .gt. 20.0) .or.   &
      	     (drymasspp(n,itype) .le. smallmassaa) ) then
	    iforce_movecenter(n,itype) = 11
	    goto 3700
	end if

	dumvtot_pre = drymasspp(n,itype)/drydenspp(n,itype)

	dumvlo_pre = volumlo_stmp(n,itype)
	dumvhi_pre = volumhi_stmp(n,itype)
	dumvdel_pre = dumvhi_pre - dumvlo_pre



	dumv = dumvhi_pre - 0.01*dumvdel_pre
	if (dumntot .le. dumvtot_pre/dumv) then
	    iforce_movecenter(n,itype) = 12
	    goto 3700
	end if
	dumv = dumvlo_pre + 0.01*dumvdel_pre
	if (dumntot .ge. dumvtot_pre/dumv) then
	    iforce_movecenter(n,itype) = 13
	    goto 3700
	end if


	dumvbar_pre = dumvtot_pre/dumntot
	dumgamma = (dumvhi_pre/dumvlo_pre) - 1.0
	dumratio = dumvbar_pre/dumvlo_pre

	if (dumratio .le. (1.0001 + dumgamma/3.0)) then
	    dumv = dumvlo_pre + 3.0*(dumvbar_pre-dumvlo_pre)
	    dumvhi_pre = min( dumvhi_pre, dumv )
	    dumvdel_pre = dumvhi_pre - dumvlo_pre
	    dumgamma = (dumvhi_pre/dumvlo_pre) - 1.0
	    dumratio = dumvbar_pre/dumvlo_pre
	else if (dumratio .ge. (0.9999 + dumgamma*2.0/3.0)) then
	    dumv = dumvhi_pre + 3.0*(dumvbar_pre-dumvhi_pre)
	    dumvlo_pre = max( dumvlo_pre, dumv )
	    dumvdel_pre = dumvhi_pre - dumvlo_pre
	    dumgamma = (dumvhi_pre/dumvlo_pre) - 1.0
	    dumratio = dumvbar_pre/dumvlo_pre
	end if

	dumbb = (dumratio - 1.0 - 0.5*dumgamma)*12.0/dumgamma
	dumaa = 1.0 - 0.5*dumbb



	dumvcutlo_nnew_pre = volumlo_stmp(nnew,itype)*(dumvbar_pre/dumvbar_aft)
	dumvcuthi_nnew_pre = volumhi_stmp(nnew,itype)*(dumvbar_pre/dumvbar_aft)




	if (nnew .eq. 1) then
	    if (dumvhi_pre .le. dumvcuthi_nnew_pre) then
		iforce_movecenter(n,itype) = 21
	    else
		nnew2 = nnew + 1
	    end if
	else if (nnew .eq. nsize_aer(itype)) then
	    if (dumvlo_pre .ge. dumvcutlo_nnew_pre) then
		iforce_movecenter(n,itype) = 22
	    else
		nnew2 = nnew - 1
	    end if
	else
	    if ((dumvlo_pre .ge. dumvcutlo_nnew_pre) .and.   &
      		(dumvhi_pre .le. dumvcuthi_nnew_pre)) then
		iforce_movecenter(n,itype) = 23
	    else if (dumvlo_pre .lt. dumvcutlo_nnew_pre) then
		nnew2 = nnew - 1
	    else
		nnew2 = nnew + 1
	    end if
	end if
	if (iforce_movecenter(n,itype) .gt. 0) goto 3700



	dumzlo = (dumvcutlo_nnew_pre - dumvlo_pre)/dumvdel_pre
	dumzhi = (dumvcuthi_nnew_pre - dumvlo_pre)/dumvdel_pre
	dumzlo = max( dumzlo, 0.0_r8 )
	dumzhi = min( dumzhi, 1.0_r8 )
	dum1 =  dumzhi    - dumzlo
	dum2 = (dumzhi**2 - dumzlo**2)*0.5
	dum3 = (dumzhi**3 - dumzlo**3)/3.0
	dumfracnum = dumaa*dum1 + dumbb*dum2
	dumfracvol = (dumvlo_pre/dumvbar_pre) * (dumaa*dum1 +   &
      		(dumaa*dumgamma + dumbb)*dum2 + (dumbb*dumgamma)*dum3)

	if ((dumfracnum .le. 0.0) .or. (dumfracvol .le. 0.0)) then
	    iforce_movecenter(n,itype) = 31
	    nnewsave(1,n,itype) = nnew2
	else if ((dumfracnum .ge. 1.0) .or. (dumfracvol .ge. 1.0)) then
	    iforce_movecenter(n,itype) = 32
	end if
	if (iforce_movecenter(n,itype) .gt. 0) goto 3700

	nnewsave(2,n,itype) = nnew2












	xferfracvol(1,kktypenew,n,itype) = dumfracvol
	xferfracvol(2,kktypenew,n,itype) = 1.0 - dumfracvol
	xferfracnum(1,kktypenew,n,itype) = dumfracnum
	xferfracnum(2,kktypenew,n,itype) = 1.0 - dumfracnum

3700	continue


	tmpa = 0.0   
	tmpb = 0.0   
	tmpc = 0.0   
	tmpd = 0.0   
	do ll = 1, ncomp_aer(itype)
	    tmpe = rmassxx(ll,n,itype)
	    tmpb = tmpb + tmpe


	    if (massptr_aer(ll,n,itype,iphase) == lptr_bc_aer(n,itype,iphase)) then
		tmpa = tmpa + tmpe
		if (method_kappa == 12) cycle
		
	    end if
	    tmpf = tmpe/dens_aer(ll,itype)
	    tmpc = tmpc + tmpf*hygro_aer(ll,itype)
	    tmpd = tmpd + tmpf
	end do
	tmp_wbc   = tmpa/max(tmpb,1.0e-35_r8)
	tmp_kappa = tmpc/max(tmpd,1.0e-35_r8)

	itype1new = ntype_md1_aer
	do jtype1 = 1, ntype_md1_aer - 1
	    if (tmp_wbc <= xcut_atype_md1(jtype1)) then
		itype1new = jtype1
		exit
	    end if
	end do

	itype2new = ntype_md2_aer
	do jtype2 = 1, ntype_md2_aer - 1
	    if (tmp_kappa <= xcut_atype_md2(jtype2)) then
		itype2new = jtype2
		exit
	    end if
	end do

	itypenew = itype_of_itype_md1md2(itype1new,itype2new)
	itype1 = itype_md1_of_itype(itype)
	itype2 = itype_md2_of_itype(itype)

	itypenewsave(:,n,itype) = 0
	itypenewsave(1,n,itype) = itypenew
	if (method_movesect == 61) itypenewsave(1,n,itype) = itype



	if (idiag_movesect .lt. 70) goto 3800

	if (nnewsave(2,n,itype) .eq. 0) then
	    if (nnewsave(1,n,itype) .eq. 0) then
		dumch4 = 'NO X'
	    else if (nnewsave(1,n,itype) .eq. n) then
		dumch4 = 'NO A'
	    else
		dumch4 = 'YESA'
	    end if
	else if (nnewsave(1,n,itype) .eq. 0) then
	    if (nnewsave(2,n,itype) .eq. n) then
		dumch4 = 'NO B'
	    else
		dumch4 = 'YESB'
	    end if
	else if (nnewsave(2,n,itype) .eq. n) then
	    if (nnewsave(1,n,itype) .eq. n) then
		dumch4 = 'NO Y'
	    else
		dumch4 = 'YESC'
	    end if
	else if (nnewsave(1,n,itype) .eq. n) then
	    dumch4 = 'YESD'
	else
	    dumch4 = 'YESE'
	end if

	dumch1 = '+'
	if (drymasspp(n,itype) .gt. drymassxx(n,itype)) dumch1 = '-'
		
	msg = ' '
	call peg_message( lunout, msg )
	write(msg,97010) dumch1, dumch4, iclm, jclm, k, m,   &
      		n, nnewsave(1,n,itype), nnewsave(2,n,itype), &
		iforce_movecenter(n,itype), itypenewsave(1,n,itype)
	call peg_message( lunout, msg )
	write(msg,97020) 'pre mass, dens      ',   &
      		drymasspp(n,itype), drydenspp(n,itype)
	call peg_message( lunout, msg )
	write(msg,97020) 'aft mass, dens, numb',   &
      		drymassxx(n,itype), drydensxx(n,itype), rnumbxx(n,itype)
	call peg_message( lunout, msg )
	if ((drydensxx(n,itype) .ne. drydensxx0(n,itype)) .or.   &
      	    (drymassxx(n,itype) .ne. drymassxx0(n,itype)) .or.   &
      	    (rnumbxx(n,itype)   .ne. rnumbxx0(n,itype)  )) then
      	    write(msg,97020) 'aft0 mas, dens, numb',   &
      		drymassxx0(n,itype), drydensxx0(n,itype), rnumbxx0(n,itype)
	    call peg_message( lunout, msg )
	end if
	write(msg,97020) 'vlop0, vbarp,  vhip0',   &
      		volumlo_stmp(n,itype), dumvbar_pre, volumhi_stmp(n,itype)
	call peg_message( lunout, msg )
	write(msg,97020) 'vlop , vbarp,  vhip ',   &
      		dumvlo_pre, dumvbar_pre, dumvhi_pre
	call peg_message( lunout, msg )
	write(msg,97020) 'vloax, vbarax, vhiax',   &
      		dumvcutlo_nnew_pre, dumvbar_pre, dumvcuthi_nnew_pre
	call peg_message( lunout, msg )
	write(msg,97020) 'vloa0, vbara,  vhia0',   &
      		volumlo_stmp(nnew,itype), dumvbar_aft, volumhi_stmp(nnew,itype)
	call peg_message( lunout, msg )
	write(msg,97020) 'dumfrvol, num, ratio',   &
      		dumfracvol, dumfracnum, dumratio
	call peg_message( lunout, msg )
	write(msg,97020) 'frvol,num1; vol,num2',   &
      		xferfracvol(1,kktypenew,n,itype), xferfracnum(1,kktypenew,n,itype),   &
      		xferfracvol(2,kktypenew,n,itype), xferfracnum(2,kktypenew,n,itype)
	call peg_message( lunout, msg )

97010	format( 'movesect', 2a, 7x, 4i3, 4x,   &
      		'n,nnews', 3i3, 4x, 'iforce', i3.2, 4x, 'itypenew', i5 )
97020	format( a, 1p, 4e13.4 )

3800	continue






	ierr = 0
	if (nnewsave(1,n,itype) .eq. nnewsave(2,n,itype)) then
	    ierr = 1
	else if (nnewsave(1,n,itype)*nnewsave(2,n,itype) .ne. 0) then
	    if (iabs(nnewsave(1,n,itype)-nnewsave(2,n,itype)) .ne. 1) ierr = 1
	end if
	if (ierr .gt. 0) then
	    write(msg,97010) 'E', 'RROR', iclm, jclm, k, m,   &
      		n, nnewsave(1,n,itype), nnewsave(2,n,itype), iforce_movecenter(n,itype)
	    call peg_message( lunout, msg )
	end if




	if ((method_movesect .ge. 30) .and. (method_movesect .le. 39)) then
	    nnewsave(1,n,itype) = nnew
	    nnewsave(2,n,itype) = 0
	    xferfracvol(:,:,n,itype) = 0.0
	    xferfracnum(:,:,n,itype) = 0.0
	    xferfracvol(1,kktypenew,n,itype) = 1.0
	    xferfracnum(1,kktypenew,n,itype) = 1.0
	end if

3900	continue   

4900	continue   

	return
	end subroutine move_sect_3d_calc_masnumadvect                          



	subroutine move_sect_3d_apply_moves(   &
	  iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
	  nnewsave, itypenewsave,    &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
          fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_apvolu,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  rbox,   &
      	  drydenspp, drymasspp, rnumbpp,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy,   &
      	  xferfracvol, xferfracnum,   &
	  dcen_stmp, volumcen_stmp, volumlo_stmp, volumhi_stmp,   &
      	  adrydens_tmp, awetdens_tmp, adrydpav_tmp, awetdpav_tmp,   &
      	  adryqmas_tmp )




	use module_data_mosaic_asecthp, only:  &
	    maxd_acomp, maxd_asize, maxd_atype,   &
	    nsize_aer, ntype_aer, ncomp_plustracer_aer, ai_phase,   &
	    massptr_aer, numptr_aer, waterptr_aer, hyswptr_aer

	implicit none


	integer iphase_flag, iclm, jclm, iphase, k,   &
      	  m, method_movesect, idiag_movesect, llhysw, llwater
	integer nnewsave(2,maxd_asize,maxd_atype)
	integer itypenewsave(4,maxd_asize,maxd_atype)

	real(r8) densh2o, smallmassaa, smallmassbb
	real(r8) fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_apvolu
	real(r8) delta_water_conform1, delta_numb_conform1
	real(r8) densdefault(maxd_atype)
	real(r8) drydenspp(maxd_asize,maxd_atype)
	real(r8) drymasspp(maxd_asize,maxd_atype)
	real(r8) drymassxx(maxd_asize,maxd_atype), drymassyy(maxd_asize,maxd_atype)
	real(r8) dryvolxx(maxd_asize,maxd_atype), dryvolyy(maxd_asize,maxd_atype)
	real(r8) rbox(ntot_used)
	real(r8) rmassxx(maxd_acomp+2,maxd_asize,maxd_atype),   &
      	     rmassyy(maxd_acomp+2,maxd_asize,maxd_atype)
	real(r8) rnumbpp(maxd_asize,maxd_atype)
	real(r8) rnumbxx(maxd_asize,maxd_atype), rnumbyy(maxd_asize,maxd_atype)
	real(r8) xferfracvol(2,4,maxd_asize,maxd_atype), xferfracnum(2,4,maxd_asize,maxd_atype)
	real(r8) dcen_stmp(maxd_asize,maxd_atype), volumcen_stmp(maxd_asize,maxd_atype), &
	    volumlo_stmp(maxd_asize,maxd_atype), volumhi_stmp(maxd_asize,maxd_atype)
	real(r8) wetvolxx(maxd_asize,maxd_atype), wetvolyy(maxd_asize,maxd_atype)
	real(r8) wetmassxx(maxd_asize,maxd_atype), wetmassyy(maxd_asize,maxd_atype)
	real(r8) adrydens_tmp(maxd_asize,maxd_atype),  awetdens_tmp(maxd_asize,maxd_atype)
	real(r8) adrydpav_tmp(maxd_asize,maxd_atype), awetdpav_tmp(maxd_asize,maxd_atype)
	real(r8) adryqmas_tmp(maxd_asize,maxd_atype)


	integer itmpa, itype, itypenew, jj, kk, l, ll, n, ndum, nnew, nold
	integer jja, jjb, jjc

	real(r8) delta_numb_conform2,   &
	  dumbot, dumnum, dumnum_at_dhi, dumnum_at_dlo,   &
      	  dumvol, dumvol1p, dumxfvol, dumxfnum, sixoverpi, third
	real(r8) dumpp(maxd_asize), dumxx(maxd_asize), dumyy(maxd_asize),   &
	  dumout(maxd_asize)

	character*160 msg
	character*8 dumch8
	character*4 dumch4


	sixoverpi = 6.0/pi
	third = 1.0/3.0






	do 1950 itype = 1, ntype_aer
	do 1940 n = 1, nsize_aer(itype)

      	itmpa = sum( itypenewsave(2:4,n,itype) )
	if ( (nnewsave(1,n,itype) .eq. n) .and.   &
      	     (nnewsave(2,n,itype) .eq. 0) .and.   &
      	     (itypenewsave(1,n,itype) .eq. itype) .and.   &
      	     (itmpa                   .eq.     0) ) then




	    drymassyy(n,itype) = drymassxx(n,itype)
	    dryvolyy(n,itype) = dryvolxx(n,itype)
	    wetmassyy(n,itype) = wetmassxx(n,itype)
	    wetvolyy(n,itype) = wetvolxx(n,itype)
	    rnumbyy(n,itype) = rnumbxx(n,itype)
	    do ll = 1, ncomp_plustracer_aer(itype) + 2
		rmassyy(ll,n,itype) = rmassxx(ll,n,itype)
	    end do

	else




	    drymassyy(n,itype) = 0.0
	    dryvolyy(n,itype) = 0.0
	    wetmassyy(n,itype) = 0.0
	    wetvolyy(n,itype) = 0.0
	    rnumbyy(n,itype) = 0.0
	    do ll = 1, ncomp_plustracer_aer(itype) + 2
		rmassyy(ll,n,itype) = 0.0
	    end do

	end if


1940	continue   
1950	continue   




	do 2950 itype = 1, ntype_aer
	do 2900 n = 1, nsize_aer(itype)


      	itmpa = sum( itypenewsave(2:4,n,itype) )
	if ( (nnewsave(1,n,itype) .eq. n) .and.   &
      	     (nnewsave(2,n,itype) .eq. 0) .and.   &
      	     (itypenewsave(1,n,itype) .eq. itype) .and.   &
      	     (itmpa                   .eq.     0) ) goto 2900

	do 2820 jj = 1, 2
	nnew = nnewsave(jj,n,itype)
	if (nnew .le. 0) goto 2820

	do 2810 kk = 1, 4
	itypenew = itypenewsave(kk,n,itype)
	if (itypenew .le. 0) goto 2810

	dumxfvol = xferfracvol(jj,kk,n,itype)
	dumxfnum = xferfracnum(jj,kk,n,itype)
	if ((dumxfvol .le. 0.0) .and. (dumxfnum .le. 0.0)) goto 2810

	do ll = 1, ncomp_plustracer_aer(itype) + 2
	    rmassyy(ll,nnew,itypenew) = rmassyy(ll,nnew,itypenew) + rmassxx(ll,n,itype)*dumxfvol
	end do
	rnumbyy(nnew,itypenew) = rnumbyy(nnew,itypenew) + rnumbxx(n,itype)*dumxfnum

	drymassyy(nnew,itypenew) = drymassyy(nnew,itypenew) + drymassxx(n,itype)*dumxfvol
	dryvolyy( nnew,itypenew) = dryvolyy( nnew,itypenew) + dryvolxx( n,itype)*dumxfvol
	wetmassyy(nnew,itypenew) = wetmassyy(nnew,itypenew) + wetmassxx(n,itype)*dumxfvol
	wetvolyy( nnew,itypenew) = wetvolyy( nnew,itypenew) + wetvolxx( n,itype)*dumxfvol

2810	continue
2820	continue

2900	continue   
2950	continue   









	call move_sect_3d_conserve_check(   &
	  1, iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
	  nnewsave, itypenewsave,    &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
	  fact_apmassmr, fact_apnumbmr,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  rbox,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy )

	delta_numb_conform2 = 0.0

	do 3950 itype = 1, ntype_aer
	do 3900 n = 1, nsize_aer(itype)

	dumvol = dryvolyy(n,itype)
	if (min(drymassyy(n,itype),dumvol) .le. smallmassbb) then
	    dumvol = drymassyy(n,itype)/densdefault(itype)
	    dumnum = drymassyy(n,itype)/(volumcen_stmp(n,itype)*densdefault(itype))
	    delta_numb_conform2 = delta_numb_conform2 + dumnum - rnumbyy(n,itype)
	    rnumbyy(n,itype) = dumnum
	    adrydens_tmp(n,itype) = densdefault(itype)/fact_apdens
	    awetdens_tmp(n,itype) = densdefault(itype)/fact_apdens
	    adrydpav_tmp(n,itype) = dcen_stmp(n,itype)/fact_apdiam
	    awetdpav_tmp(n,itype) = dcen_stmp(n,itype)/fact_apdiam
	else
	    dumnum = rnumbyy(n,itype)
	    dumnum_at_dhi = dumvol/volumhi_stmp(n,itype)
	    dumnum_at_dlo = dumvol/volumlo_stmp(n,itype)
	    dumnum = max( dumnum_at_dhi, min( dumnum_at_dlo, dumnum ) )
	    delta_numb_conform2 = delta_numb_conform2 + dumnum - rnumbyy(n,itype)
	    rnumbyy(n,itype) = dumnum
	    adrydens_tmp(n,itype) = drymassyy(n,itype)/dumvol/fact_apdens
	    dumvol1p = dumvol/dumnum
	    adrydpav_tmp(n,itype) = ((dumvol1p*sixoverpi)**third) / fact_apdiam
	    awetdens_tmp(n,itype) = wetmassyy(n,itype)/wetvolyy(n,itype)/fact_apdens
	    dumvol1p = wetvolyy(n,itype)/dumnum
	    awetdpav_tmp(n,itype) = min( 100.0_r8*dcen_stmp(n,itype),   &
      			(dumvol1p*sixoverpi)**third ) / fact_apdiam
	end if
	adryqmas_tmp(n,itype) = drymassyy(n,itype)/fact_apmassmr

	do ll = 1, ncomp_plustracer_aer(itype)
	    l = massptr_aer(ll,n,itype,iphase)
	    rbox(l) = rmassyy(ll,n,itype)/fact_apmassmr
	end do
	l = 0
	if (iphase .eq. ai_phase) then
	    l = waterptr_aer(n,itype)
	    if (l .gt. 0) rbox(l) = rmassyy(llwater,n,itype)/fact_apmassmr
	    l = hyswptr_aer(n,itype)
	    if (l .gt. 0) rbox(l) = rmassyy(llhysw,n,itype)/fact_apmassmr
	end if
	rbox(numptr_aer(n,itype,iphase)) = rnumbyy(n,itype)/fact_apnumbmr

3900	continue   
3950	continue   

	delta_numb_conform1 = delta_numb_conform1 + delta_numb_conform2

	call move_sect_3d_conserve_check(   &
	  2, iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
	  nnewsave, itypenewsave,    &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
	  fact_apmassmr, fact_apnumbmr,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  rbox,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy )



	if (idiag_movesect .lt. 70) goto 4900

	do 4800 itype = 1, ntype_aer

	msg = ' '
	call peg_message( lunout, msg )
	write(msg,97005) itype
	call peg_message( lunout, msg )

	ndum = 0
	do n = 1, nsize_aer(itype)
	    if (nnewsave(1,n,itype)+nnewsave(2,n,itype) .ne. n) ndum = ndum + 1
	end do

	dumch4 = 'NONE'
	if (ndum .gt. 0) dumch4 = 'SOME'
	msg = ' '
	call peg_message( lunout, msg )
	write(msg,97010) dumch4, iclm, jclm, k, m, ndum
	call peg_message( lunout, msg )
	do jjb = 1, nsize_aer(itype), 10
	    jjc = min( jjb+9, nsize_aer(itype) )
	    write(msg,97011) (nnewsave(1,n,itype), nnewsave(2,n,itype), n=jjb,jjc)
	    call peg_message( lunout, msg )
	end do











	dumbot = log( volumhi_stmp(1,itype)/volumlo_stmp(1,itype) )
	do n = 1, nsize_aer(itype)
	    dumpp(n) = -9.99
	    dumxx(n) = -9.99
	    dumyy(n) = -9.99
	    if ( (drydenspp(n,itype) .gt. 0.5) .and.   &
      	         (drymasspp(n,itype) .gt. smallmassaa) ) then
      		dumvol = drymasspp(n,itype)/drydenspp(n,itype)
		if ((rnumbpp(n,itype) .ge. 1.0e-35) .and.   &
      		    (dumvol .ge. 1.0e-35)) then
		    dumvol1p = dumvol/rnumbpp(n,itype)
		    dumpp(n) = 1.0 + log(dumvol1p/volumlo_stmp(1,itype))/dumbot
		end if
	    end if
	    if ((rnumbxx(n,itype) .ge. 1.0e-35) .and.   &
      		(dryvolxx(n,itype) .ge. 1.0e-35)) then
		dumvol1p = dryvolxx(n,itype)/rnumbxx(n,itype)
		dumxx(n) = 1.0 + log(dumvol1p/volumlo_stmp(1,itype))/dumbot
	    end if
	    if ((rnumbyy(n,itype) .ge. 1.0e-35) .and.   &
      		(dryvolyy(n,itype) .ge. 1.0e-35)) then
		dumvol1p = dryvolyy(n,itype)/rnumbyy(n,itype)
		dumyy(n) = 1.0 + log(dumvol1p/volumlo_stmp(1,itype))/dumbot
	    end if
	end do






	do jja = 1, 7

	    if      (jja .eq. 1) then
		dumch8 = 'rnumbold'
		dumout(:) = rnumbxx(:,itype)
	    else if (jja .eq. 2) then
		dumch8 = 'rnumbnew'
		dumout(:) = rnumbyy(:,itype)
	    else if (jja .eq. 3) then
		dumch8 = 'drvolold'
		dumout(:) = dryvolxx(:,itype)
	    else if (jja .eq. 4) then
		dumch8 = 'drvolnew'
		dumout(:) = dryvolyy(:,itype)
	    else if (jja .eq. 5) then
		dumch8 = 'lnvolold'
		dumout(:) = dumxx(:)
	    else if (jja .eq. 6) then
		dumch8 = 'lnvolnew'
		dumout(:) = dumyy(:)
	    else if (jja .eq. 7) then
		dumch8 = 'lnvolpre'
		dumout(:) = dumpp(:)
	    else if (jja .eq. 8) then
		dumch8 = 'wtvolold'
		dumout(:) = wetvolxx(:,itype)
	    else if (jja .eq. 9) then
		dumch8 = 'wtvolnew'
		dumout(:) = wetvolyy(:,itype)
	    else if (jja .eq. 10) then
		dumch8 = 'wtmasold'
		dumout(:) = wetmassxx(:,itype)
	    else if (jja .eq. 11) then
		dumch8 = 'wtmasnew'
		dumout(:) = wetmassyy(:,itype)
	    else if (jja .eq. 12) then
		dumch8 = 'drmasold'
		dumout(:) = drymassxx(:,itype)
	    else if (jja .eq. 13) then
		dumch8 = 'drmasnew'
		dumout(:) = drymassyy(:,itype)
	    end if
	    do jjb = 1, nsize_aer(itype), 10
		jjc = min( jjb+9, nsize_aer(itype) )
		if ((jja .le. 4) .or. (jja .ge. 8)) then
		    if (jjc == nsize_aer(itype)) then
			write(msg,97020) dumch8, (dumout(n), n=jjb,jjc), &
			    sum( dumout(1:jjc) )
		    else
			write(msg,97020) dumch8, (dumout(n), n=jjb,jjc)
		    end if
		else
		    write(msg,97030) dumch8, (dumout(n), n=jjb,jjc)
		end if
		call peg_message( lunout, msg )
		dumch8 = ' '
	    end do
	end do

97005	format( 'movesectapply - itype =', i6 )
97010	format( 'movesectapply', a, 4i3, i6 )



97011	format( 5x, 10(5x,2i3) )
97020	format( a, 1p, 11e11.3 )
97030	format( a,     10f11.5 )

4800	continue   

4900	continue
	return
	end subroutine move_sect_3d_apply_moves                          



	subroutine move_sect_3d_conserve_check( ipass,   &
	  iphase_flag, iclm, jclm, k, m, iphase,   &
      	  method_movesect, idiag_movesect, llhysw, llwater,   &
	  nnewsave, itypenewsave,    &
      	  densdefault, densh2o, smallmassaa, smallmassbb,   &
	  fact_apmassmr, fact_apnumbmr,   &
      	  delta_water_conform1, delta_numb_conform1,   &
      	  rbox,   &
      	  drymassxx, dryvolxx, rmassxx, rnumbxx, wetmassxx, wetvolxx,   &
      	  drymassyy, dryvolyy, rmassyy, rnumbyy, wetmassyy, wetvolyy )



















	use module_data_mosaic_asecthp, only:  &
	    maxd_acomp, maxd_asize, maxd_atype,   &
	    nsize_aer, ntype_aer, ncomp_plustracer_aer, ai_phase,   &
	    massptr_aer, numptr_aer, waterptr_aer, hyswptr_aer,   &
	    name_aer

	implicit none


	integer ipass, iphase_flag, iclm, jclm, iphase, k,   &
	  m, method_movesect, idiag_movesect, llhysw, llwater
	integer nnewsave(2,maxd_asize,maxd_atype)
	integer itypenewsave(2,maxd_asize,maxd_atype)
	real(r8) densh2o, smallmassaa, smallmassbb
	real(r8) fact_apmassmr, fact_apnumbmr
	real(r8) delta_water_conform1, delta_numb_conform1
	real(r8) densdefault(maxd_atype)
	real(r8) drymassxx(maxd_asize,maxd_atype), drymassyy(maxd_asize,maxd_atype)
	real(r8) dryvolxx(maxd_asize,maxd_atype), dryvolyy(maxd_asize,maxd_atype)
	real(r8) rbox(ntot_used)
	real(r8) rmassxx(maxd_acomp+2,maxd_asize,maxd_atype),   &
      	     rmassyy(maxd_acomp+2,maxd_asize,maxd_atype)
	real(r8) rnumbxx(maxd_asize,maxd_atype), rnumbyy(maxd_asize,maxd_atype)
	real(r8) wetvolxx(maxd_asize,maxd_atype), wetvolyy(maxd_asize,maxd_atype)
	real(r8) wetmassxx(maxd_asize,maxd_atype), wetmassyy(maxd_asize,maxd_atype)


	integer ii, itype, itmpa, itmpb, jj, jtype, l, ll, llworst, llworstb, n
	integer nerr, nerrmax
	save nerr, nerrmax
	data nerr, nerrmax / 0, 999 /

	real(r8) dumbot, dumtop, dumtoler, dumerr, dumworst, dumworstb
	real(r8) duma, dumb, dumc, dume
	real(r8) dumerrsv(maxd_acomp+7)
	real(r8) thesum(4,maxd_acomp+7)

	character*8 dumname(maxd_acomp+7)
	character*160 fmtaa, msg


	if (ipass .eq. 2) goto 2000

	itype = 1
	do ll = 1, ncomp_plustracer_aer(itype)+7
	do jj = 1, 4
	    thesum(jj,ll) = 0.0
	end do
	end do

	do itype = 1, ntype_aer
	do n = 1, nsize_aer(itype)
	    do ll = 1, ncomp_plustracer_aer(itype)+2
		thesum(1,ll) = thesum(1,ll) + rmassxx(ll,n,itype)
		thesum(2,ll) = thesum(2,ll) + rmassyy(ll,n,itype)
	    end do
	    ll = ncomp_plustracer_aer(itype)+3
	    thesum(1,ll) = thesum(1,ll) + rnumbxx(n,itype)
	    thesum(2,ll) = thesum(2,ll) + rnumbyy(n,itype)
	    ll = ncomp_plustracer_aer(itype)+4
	    thesum(1,ll) = thesum(1,ll) + drymassxx(n,itype)
	    thesum(2,ll) = thesum(2,ll) + drymassyy(n,itype)
	    ll = ncomp_plustracer_aer(itype)+5
	    thesum(1,ll) = thesum(1,ll) + dryvolxx(n,itype)
	    thesum(2,ll) = thesum(2,ll) + dryvolyy(n,itype)
	    ll = ncomp_plustracer_aer(itype)+6
	    thesum(1,ll) = thesum(1,ll) + wetmassxx(n,itype)
	    thesum(2,ll) = thesum(2,ll) + wetmassyy(n,itype)
	    ll = ncomp_plustracer_aer(itype)+7
	    thesum(1,ll) = thesum(1,ll) + wetvolxx(n,itype)
	    thesum(2,ll) = thesum(2,ll) + wetvolyy(n,itype)
	end do
	end do






2000	continue
	do itype = 1, ntype_aer
	do n = 1, nsize_aer(itype)
	    do ll = 1, ncomp_plustracer_aer(itype)+3
		duma = fact_apmassmr
		if (ll .le. ncomp_plustracer_aer(itype)) then
		    l = massptr_aer(ll,n,itype,iphase)
		else if (ll .eq. ncomp_plustracer_aer(itype)+1) then
		    l = 0
		    if (iphase .eq. ai_phase) l = hyswptr_aer(n,itype)
		else if (ll .eq. ncomp_plustracer_aer(itype)+2) then
		    l = 0
		    if (iphase .eq. ai_phase) l = waterptr_aer(n,itype)
		else
		    l = numptr_aer(n,itype,iphase)
		    duma = fact_apnumbmr
		end if
		if (l .gt. 0)   &
		    thesum(ipass+2,ll) = thesum(ipass+2,ll) + rbox(l)*duma
	    end do
	end do
	end do

	itype = 1
	if (ipass .eq. 2) then
	    ll = ncomp_plustracer_aer(itype)+2
	    thesum(3,ll) = thesum(3,ll) + delta_water_conform1
	    ll = ncomp_plustracer_aer(itype)+3
	    thesum(3,ll) = thesum(3,ll) + delta_numb_conform1
	end if






	itype = 1
	do ll = 1, ncomp_plustracer_aer(itype)+7
	    dumname(ll) = ' '
	    write(dumname(ll),'(i4.4)') ll


	    if (ll .le. ncomp_plustracer_aer(itype))   dumname(ll) = name_aer(ll,1)(1:4)
	    if (ll .eq. ncomp_plustracer_aer(itype)+1) dumname(ll) = 'hysw'
	    if (ll .eq. ncomp_plustracer_aer(itype)+2) dumname(ll) = 'watr'
	    if (ll .eq. ncomp_plustracer_aer(itype)+3) dumname(ll) = 'numb'
	    if (ll .eq. ncomp_plustracer_aer(itype)+4) dumname(ll) = 'drymass'
	    if (ll .eq. ncomp_plustracer_aer(itype)+5) dumname(ll) = 'dryvol'
	    if (ll .eq. ncomp_plustracer_aer(itype)+6) dumname(ll) = 'wetmass'
	    if (ll .eq. ncomp_plustracer_aer(itype)+7) dumname(ll) = 'wetvol'
	end do

	jj = 2*ipass - 1
	dumworst = 0.0
	dumworstb = 0.0
	dumerrsv(:) = 0.0
	llworst = 0
	llworstb = 0
	do ll = 1, ncomp_plustracer_aer(itype)+7
	    dumtop = thesum(jj+1,ll) - thesum(jj,ll)
	    dumbot = max( abs(thesum(jj,ll)), abs(thesum(jj+1,ll)), 1.0e-35_r8 )
	    dumerr = dumtop/dumbot







	    if (ipass .eq. 2) then
		dumc = 1.0
		if (ll .eq. ncomp_plustracer_aer(itype)+2) then
		    dumc = delta_water_conform1
		else if (ll .eq. ncomp_plustracer_aer(itype)+3) then
		    dumc = delta_numb_conform1
		end if
		if (dumc .lt. 0.0) then
		    duma = thesum(3,ll) - dumc
		    dumb = thesum(4,ll) - dumc
		    dumtop = dumb - duma
		    dumbot = max( abs(duma), abs(dumb), 1.0e-35_r8 )
		    dume = dumtop/dumbot
		    if (abs(dume) .lt. abs(dumerr)) dumerr = dume
		end if
	    end if
	    dumerrsv(ll) = dumerr

	    if (abs(dumerr) .gt. abs(dumworst)) then
		llworstb = llworst
		dumworstb = dumworst
		llworst = ll
		dumworst = dumerr
	    end if
	end do

	dumtoler = 1.0e-6
	if (abs(dumworst) .gt. dumtoler) then
	    nerr = nerr + 1
	    if (nerr .le. nerrmax) then
		msg = ' '
		call peg_message( lunout, msg )
		write(msg,97110) iclm, jclm, k, m, ipass, llworst
		call peg_message( lunout, msg )
		write(msg,'(a)') '    jtype, ii, nnew(ii,1:n,jtype)'
		call peg_message( lunout, msg )

		do jtype = 1, ntype_aer
		do ii = 1, 2
		do itmpa = 1, nsize_aer(jtype), 24
		    itmpb = min( itmpa+23, nsize_aer(jtype) )  
		    if (itmpa == 1) then
			fmtaa = '(i4,i2,2x,24i4)'
			write(msg,fmtaa) jtype, ii, &
			    (nnewsave(ii,n,jtype), n=itmpa,itmpb)
		    else
			fmtaa = '(8x,24i4)'
			write(msg,fmtaa) (nnewsave(ii,n,jtype), n=itmpa,itmpb)
		    end if
	    	    call peg_message( lunout, msg )
		end do
		end do
		end do

		ll = llworst
		if (ll .eq. 0) ll = ncomp_plustracer_aer(itype)+2
		write(msg,97130) 'name/relerrA/thesum', jj, '/thesum', jj+1,   &
      			dumname(ll), dumworst, thesum(jj,ll), thesum(jj+1,ll)
		call peg_message( lunout, msg )

		if ( (ll .eq. ncomp_plustracer_aer(itype)+3) .and.   &
		     (abs(dumworstb) .gt. dumtoler) ) then
		    ll = max( 1, llworstb )
		    dumtop = thesum(jj+1,ll) - thesum(jj,ll)
		    dumbot = max( abs(thesum(jj,ll)), abs(thesum(jj+1,ll)), 1.0e-35_r8 )
		    dumerr = dumtop/dumbot
		    write(msg,97130) 'name/relerrB/thesum', jj, '/thesum', jj+1,   &
      			dumname(ll), dumerr, thesum(jj,ll), thesum(jj+1,ll)
		    call peg_message( lunout, msg )
		end if
	    end if
	end if


	do ll = 1, ncomp_plustracer_aer(itype)+7

	    if (abs(dumerrsv(ll)) .le. dumtoler) cycle
	    nerr = nerr + 1
	    if (nerr .le. nerrmax) then
		write(msg,97130) 'name/relerrC/thesum', jj, '/thesum', jj+1,   &
      			dumname(ll), dumworst, thesum(jj,ll), thesum(jj+1,ll)
		call peg_message( lunout, msg )
	    end if
	end do 

97110	format( 'movesect conserve ERROR - i/j/k/m/pass/llworst',   &
		4i3, 2x, 2i3 )
97120	format( 64i3 )
97130	format( a, i1, a, i1, 2x, a, 1p, 3e16.7 )


	return
	end subroutine move_sect_3d_conserve_check        



	subroutine move_sect_3d_checkptrs( iphase_flag, iclm, jclm, k, m )



	use module_data_mosaic_asecthp, only:  &
	    ai_phase, ntype_aer, nsize_aer, nphase_aer,   &
	    ncomp_aer, ncomp_plustracer_aer,   &
	    numptr_aer, waterptr_aer, mastercompptr_aer

	implicit none


	integer iphase_flag, iclm, jclm, k, ll, m


	integer l, itype, iphase, n, ndum
	character*160 msg

	do 1900 itype = 1, ntype_aer
	do 1800 iphase = 1, nphase_aer

	ndum = 0
	do n = 1, nsize_aer(itype)
	    l = numptr_aer(n,itype,iphase)
	    if ((l .lt. 1) .or. (l .gt. ntot_used)) then
		msg = '*** subr move_sect_3d error - ' //   &
			'numptr_amode not defined'
		call peg_message( lunerr, msg )
		write(msg,9030) 'mode, numptr =', n, l
		call peg_message( lunerr, msg )
		write(msg,9030) 'iphase, itype =', iphase, itype
		call peg_message( lunerr, msg )
		call peg_error_fatal( lunerr, msg )
	    end if


	    l = 0
	    if (iphase .eq. ai_phase) l = waterptr_aer(n,itype)
	    if ((l .ge. 1) .and. (l .le. ntot_used)) ndum = ndum + 1
	end do
	if ((ndum .ne. 0) .and. (ndum .ne. nsize_aer(itype))) then
	    msg = '*** subr move_sect_3d error - ' //   &
      		'waterptr_aer must be on/off for all modes'
	    call peg_message( lunerr, msg )
	    write(msg,9030) 'iphase, itype =', iphase, itype
	    call peg_message( lunerr, msg )
	    call peg_error_fatal( lunerr, msg )
	end if
9030	format( a, 6(1x,i6) )

	msg = ' '
	if (nsize_aer(itype) /= nsize_aer(1)) then
	    write(msg,9030) 'nsize diffs', &
		1, itype, nsize_aer(1), nsize_aer(itype)
	else if (ncomp_aer(itype) /= ncomp_aer(1)) then
	    write(msg,9030) 'ncomp diffs', &
		1, itype, ncomp_aer(1), ncomp_aer(itype)
	else if (ncomp_plustracer_aer(itype) /= ncomp_plustracer_aer(1)) then
	    write(msg,9030) 'ncomp_plustracer diffs', &
		1, itype, ncomp_plustracer_aer(1), ncomp_plustracer_aer(itype)
	end if
	do ll = 1, ncomp_plustracer_aer(1)
	    if (mastercompptr_aer(ll,itype) /= mastercompptr_aer(ll,1)) then
		write(msg,9030) 'mastercompptrcomp diffs', &
		    1, itype, ll, mastercompptr_aer(ll,1), mastercompptr_aer(ll,itype)
	    end if
	end do
	if (msg /= ' ') then
	    call peg_message( lunerr, msg )
	    call peg_error_fatal( lunerr, msg )
	end if

1800	continue
1900	continue

	return
	end subroutine move_sect_3d_checkptrs                           



	subroutine test_move_sect_3d( iphase_flag_inp, iclm, jclm, k, m, &
             it_mosaic, mmovesect_flag1, idiag_movesect )




	use module_data_mosaic_asecthp, only:  &
	    maxd_atype, maxd_asize,   &
	    ntype_aer, nsize_aer, ncomp_plustracer_aer,   &
	    ai_phase, cw_phase,   &
	    massptr_aer, numptr_aer, waterptr_aer,   &
	    dens_aer, volumlo_sect

	implicit none


	integer, intent(in) :: iphase_flag_inp, iclm, jclm, k, m, it_mosaic, &
	                       mmovesect_flag1, idiag_movesect


	integer iphase_flag, ii, iphase, itype, jj,   &
	  l, ll, n, nn
	integer, save :: ientryno = 0

	real(r8) dumnumb, dumvolpre, dumvolaft
	real(r8) adrydens_tmp(maxd_asize,maxd_atype),  awetdens_tmp(maxd_asize,maxd_atype)
	real(r8) adrydpav_tmp(maxd_asize,maxd_atype), awetdpav_tmp(maxd_asize,maxd_atype)
	real(r8) adryqmas_tmp(maxd_asize,maxd_atype)
	real(r8) fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam
	real(r8) tmp_drymass_pregrow(maxd_asize,maxd_atype)
	real(r8) tmp_drymass_aftgrow(maxd_asize,maxd_atype)
	real(r8) tmp_drydens_pregrow(maxd_asize,maxd_atype)
	real(r8) tmp_drydens_aftgrow(maxd_asize,maxd_atype)
	real(r8) rbox(ntot_used)

	character*160 msg

	integer maxvolfactpre, maxvolfactaft
	parameter (maxvolfactpre=15, maxvolfactaft=23)

	real(r8) dumvolfactpre(maxvolfactpre)
	data dumvolfactpre /   &
	2.0, 0.0, 1.0e-20, 0.5, 0.9,   &
	1.0, 1.01, 1.1, 2.0, 4.0, 7.9, 7.99, 8.0,   &
	8.1, 16.0 /

	real(r8) dumvolfactaft(maxvolfactaft)
	data dumvolfactaft /   &
	4.0, 0.0, 1.0e-20, 0.01, 0.02, 0.05, 0.1, 0.5, 0.9,   &
	1.0, 1.01, 1.1, 2.0, 4.0, 7.9, 7.99, 8.0,   &
	8.1, 16.0, 32., 64., 128., 256. /









	if (mmovesect_flag1 .le. 0) return
	if (nsize_aer(1) .le. 0) return

	ientryno = ientryno + 1
	if (ientryno .gt. 1) return

	write(*,*) '*** doing move_sect_3d_checkptrs ***'
	call move_sect_3d_checkptrs( -1, iclm, jclm, k, m )

	if (idiag_movesect .ne. 80) return
	
	write(*,*) '*** executing test_move_sect_3d ***'





	do 3900 iphase_flag = 1, 1

	iphase = ai_phase
	if (iabs(iphase_flag) .eq. 2) iphase = cw_phase

	do 3800 itype = 1, ntype_aer

	do 2900 nn = 1, nsize_aer(itype)

	do 2800 ii = 1, maxvolfactpre

	do 2700 jj = 1, maxvolfactaft














	do n = 1, nsize_aer(itype)
	    do ll = 1, ncomp_plustracer_aer(itype)
		rbox(massptr_aer(ll,n,itype,iphase)) = 0.0
	    end do
	    l = 0
	    if (iphase .eq. ai_phase) l = waterptr_aer(n,itype)
	    if (l .gt. 0) rbox(l) = 0.0
	    l = numptr_aer(n,itype,iphase)
	    if (l .gt. 0) rbox(l) = 0.0
	    tmp_drymass_pregrow(n,itype) = 0.0
	    tmp_drymass_aftgrow(n,itype) = 0.0
	    tmp_drydens_pregrow(n,itype) = -1.0
	    tmp_drydens_aftgrow(n,itype) = -1.0
	end do


	n = nn
	dumnumb = 1.0e7/mw_air   
	rbox(numptr_aer(n,itype,iphase)) = dumnumb
	ll = 1
	l = massptr_aer(ll,n,itype,iphase)

	dumvolpre = volumlo_sect(n,itype)*dumvolfactpre(ii)*dumnumb
	tmp_drydens_pregrow(n,itype) = dens_aer(ll,itype)
	tmp_drymass_pregrow(n,itype) = dumvolpre*tmp_drydens_pregrow(n,itype)
	if (ii .eq. 1) tmp_drydens_pregrow(n,itype) = -1.0
	
	dumvolaft = volumlo_sect(n,itype)*dumvolfactaft(jj)*dumnumb
	tmp_drydens_aftgrow(n,itype) = dens_aer(ll,itype)
	tmp_drymass_aftgrow(n,itype) = dumvolaft*tmp_drydens_aftgrow(n,itype)
	if (jj .eq. 1) tmp_drydens_aftgrow(n,itype) = -1.0
	
	rbox(l) = tmp_drymass_aftgrow(n,itype)   
	
	msg = ' '
	call peg_message( lunout, msg )
	write(msg,98010) nn, ii, jj
	call peg_message( lunout, msg )
	write(msg,98011) dumvolfactpre(ii), dumvolfactaft(jj)
	call peg_message( lunout, msg )

	fact_apmassmr = 1.0   
	fact_apnumbmr = 1.0   
	fact_apdens = 1.0     
	fact_apdiam = 1.0     


	call move_sect_3d_x1( iphase_flag, iclm, jclm, k, m, rbox,   &
          fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam,   &
      	  tmp_drydens_aftgrow, tmp_drydens_pregrow,   &
      	  tmp_drymass_aftgrow, tmp_drymass_pregrow,   &
      	  adrydens_tmp, awetdens_tmp, adrydpav_tmp, awetdpav_tmp,   &
      	  adryqmas_tmp, it_mosaic, mmovesect_flag1, idiag_movesect )

	msg = ' '
	call peg_message( lunout, msg )
	write(msg,98010) nn, ii, jj
	call peg_message( lunout, msg )
	write(msg,98011) dumvolfactpre(ii), dumvolfactaft(jj)
	call peg_message( lunout, msg )

2700	continue
2800	continue
2900	continue

3800	continue
3900	continue

98010	format( 'test_move_sect_3d output - nn, ii, jj =', 3i3 )
98011	format( 'volfactpre, volfactaft =', 1p, 2e12.4 )


	write(*,*) '*** leaving   test_move_sect_3d at end ***'
	return
	end subroutine test_move_sect_3d                           





	end module module_mosaic_movesect3d

