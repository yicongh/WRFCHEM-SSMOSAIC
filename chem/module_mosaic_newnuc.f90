







	module module_mosaic_newnuc



	use module_peg_util



	implicit none



	contains




	subroutine mosaic_newnuc_1clm( istat_newnuc,   &
		it, jt, kclm_calcbgn, kclm_calcend,   &
		idiagbb_in,   &
		dtchem, dtnuc_in, rsub0,   &
		id, ktau, ktauc, its, ite, jts, jte, kts, kte, &
                newnuc_method, ion_prod_rate_1d, &
                JBNRATE_1D, JBCHRATE_1D, JTNRATE_1D, JTCHRATE_1D, &
                JBNORGRATE_1D, JBCHORGRATE_1D, JTNRICRATE_1D, JRATE_1D, &
                JBNRATEAC_1D, JBCHRATEAC_1D, JTNRATEAC_1D, JTCHRATEAC_1D, &
                JBNORGRATEAC_1D, JBCHORGRATEAC_1D, JTNRICRATEAC_1D, JRATEAC_1D ) 
















	use module_data_mosaic_asect
	use module_data_mosaic_other
	use module_state_description, only:  param_first_scalar


	integer, intent(inout) :: istat_newnuc    
	integer, intent(in) ::   &
		it, jt, kclm_calcbgn, kclm_calcend,   &
		idiagbb_in,   &
		id, ktau, ktauc, its, ite, jts, jte, kts, kte
	real, intent(in) :: dtchem, dtnuc_in
	real, intent(in) :: rsub0(l2maxd,kmaxd,nsubareamaxd)



        real, intent(in) :: ion_prod_rate_1d(kclm_calcbgn:kclm_calcend)
        integer, intent(in) :: newnuc_method
        real, intent(inout) :: JBNRATE_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JBCHRATE_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JTNRATE_1D(kclm_calcbgn:kclm_calcend) 
        real, intent(inout) :: JTCHRATE_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JBNORGRATE_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JBCHORGRATE_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JTNRICRATE_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JRATE_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JBNRATEAC_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JBCHRATEAC_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JTNRATEAC_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JTCHRATEAC_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JBNORGRATEAC_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JBCHORGRATEAC_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JTNRICRATEAC_1D(kclm_calcbgn:kclm_calcend)
        real, intent(inout) :: JRATEAC_1D(kclm_calcbgn:kclm_calcend)
















	integer :: k, l, ll, m, n   
	integer :: isize, itype, iphase
	integer :: iconform_numb
	integer :: idiagbb
	integer :: nsize, ntau_nuc
	integer, save :: ncount(10)
	integer :: p1st

	real, parameter :: densdefault = 2.0
	real, parameter :: smallmassbb = 1.0e-30


        real, parameter :: a_zsr_xx1 =  1.15510
        real, parameter :: a_zsr_xx2 = -3.20815
        real, parameter :: a_zsr_xx3 =  2.71141
        real, parameter :: a_zsr_xx4 =  2.01155
        real, parameter :: a_zsr_xx5 = -4.71014
        real, parameter :: a_zsr_xx6 =  2.04616
        real, parameter :: b_zsr_xx  = 29.4779


        real, parameter :: MU = 1.2E-4 
        real, parameter :: E_ELEC = 1.6022E-19 
        real, parameter :: PPI = 3.1415926536 
        real, parameter :: ZBOLTZ = 1.38e-23 


	real :: aw
	real :: cair_box
	real :: dens_nh4so4a, dtnuc
	real :: duma, dumb, dumc
	real :: rh_box
	real :: qh2so4_avg, qh2so4_cur, qh2so4_del 
	real :: qnh3_avg, qnh3_cur, qnh3_del
        real :: qnuma_del, qso4a_del, qnh4a_del
 
        real :: qbiog1_c_avg, qbiog1_c_cur, qbiog1_c_del
        real :: qbiog2_c_avg, qbiog2_c_cur, qbiog2_c_del
	real :: qbiog1_ca_cur, qbiog2_ca_cur, qbiog1_ca_del, qbiog2_ca_del

	real :: temp_box
	real :: xxdens, xxmass, xxnumb, xxvolu

	real,save :: dumveca(10), dumvecb(10), dumvecc(10), dumvecd(10), dumvece(10)
	real :: volumlo_nuc(maxd_asize), volumhi_nuc(maxd_asize)


        real :: CRII_G  
        real :: CSINK_ION  
        real :: TMP_CS 
        real :: ND 
        real :: DP 


	character(len=100) :: msg

    p1st = PARAM_FIRST_SCALAR


	istat_newnuc = 0
	if ((newnuc_method .ne. 2) .and. (newnuc_method .ne. 3)) then  
	    if ((it .eq. its) .and. (jt .eq. jts))   &
		call peg_message( lunerr,   &
		'*** mosaic_newnuc_1clm -- illegal newnuc_method' )
	    istat_newnuc = -1
	    return
	end if

        if ((newnuc_method == 3) .and. (nsize_aer(1) .ne. 12) .and. (nsize_aer(1) .ne. 20)) then
            if ((it .eq. its) .and. (jt .eq. jts))   &
                call peg_error_fatal( lunerr,   &
                'The glomap nucleation option only works for 20-bin and 12-bin schemes' )
        end if



	ntau_nuc = nint( dtnuc_in/dtchem )
	ntau_nuc = max( 1, ntau_nuc )
	if (mod(ktauc,ntau_nuc) .ne. 0) return  
	dtnuc = dtchem*ntau_nuc



	idiagbb = idiagbb_in

	itype = 1
	iphase = ai_phase
	nsize = nsize_aer(itype)
	volumlo_nuc(1:nsize) = volumlo_sect(1:nsize,itype)
	volumhi_nuc(1:nsize) = volumhi_sect(1:nsize,itype)



	do 2900 m = 1, nsubareas

	do 2800 k = kclm_calcbgn, kclm_calcend



	if ((it .eq. its) .and.   &
	    (jt .eq. jts) .and. (k .eq. kclm_calcbgn)) then
	    dumveca(:) = 0.0         
	    dumvecb(:) = +1.0e35     
	    dumvecc(:) = -1.0e35     
	    dumvecd(:) = 0.0         
	    dumvece(:) = 0.0         
	    ncount(:) = 0
	end if


	ncount(1) = ncount(1) + 1
	if (afracsubarea(k,m) .lt. 1.e-4) goto 2700

	cair_box = cairclm(k)
	temp_box = rsub(ktemp,k,m)
	rh_box = relhumclm(k)

	qh2so4_cur = max(0.0,rsub(kh2so4,k,m))
	qnh3_cur   = max(0.0,rsub(knh3,k,m))
	qh2so4_avg = 0.5*( qh2so4_cur + max(0.0,rsub0(kh2so4,k,m)) )
	qnh3_avg   = 0.5*( qnh3_cur   + max(0.0,rsub0(knh3,k,m)) )

        if (newnuc_method == 3) then
        qbiog1_c_cur   = max(0.0,rsub(kbiog1_c,k,m))
        qbiog1_c_avg   = 0.5*( qbiog1_c_cur   + max(0.0,rsub0(kbiog1_c,k,m)) )
        qbiog2_c_cur   = max(0.0,rsub(kbiog2_c,k,m))
        qbiog2_c_avg   = 0.5*( qbiog2_c_cur   + max(0.0,rsub0(kbiog2_c,k,m)) )

        qbiog1_ca_cur = 0.0
        qbiog2_ca_cur = 0.0
        do n = 1,nsize
           l = lptr_biog1_c_aer(n,itype,iphase)
           if (l .ge. p1st) then
           qbiog1_ca_cur = qbiog1_ca_cur + rsub(l,k,m)
           end if
           l = lptr_biog2_c_aer(n,itype,iphase)
           if (l .ge. p1st) then
           qbiog2_ca_cur = qbiog2_ca_cur + rsub(l,k,m)
           end if
        end do

        CRII_G = ion_prod_rate_1d(k)
        CSINK_ION = 0.0
        TMP_CS = 4 * PPI * ZBOLTZ * temp_box * MU / E_ELEC 
        do n = 1, nsize
          ND = rsub(numptr_aer(n,itype,iphase),k,m)* cair_box  
          DP = dcen_sect(n,itype) *0.01 
          CSINK_ION = CSINK_ION +TMP_CS *0.5* DP *ND * 1.E6
        end do
        end if


	qh2so4_del = 0.0
	qnh3_del = 0.0
	qnuma_del = 0.0
	qso4a_del = 0.0
	qnh4a_del = 0.0

        qbiog1_c_del = 0.0
        qbiog1_ca_del = 0.0
        qbiog2_c_del = 0.0
        qbiog2_ca_del = 0.0


	dens_nh4so4a = dens_so4_aer

	isize = 0


	if (newnuc_method .eq. 1) then
            call ternary_nuc_mosaic_1box(   &
               dtnuc, temp_box, rh_box, cair_box,   &
               qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
               nsize, maxd_asize, volumlo_nuc, volumhi_nuc,   &
               isize, qnuma_del, qso4a_del, qnh4a_del,   &
               qh2so4_del, qnh3_del, dens_nh4so4a )
	else if (newnuc_method .eq. 2) then
            call wexler_nuc_mosaic_1box(   &
               dtnuc, temp_box, rh_box, cair_box,   &
               qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
               nsize, maxd_asize, volumlo_nuc, volumhi_nuc,   &
               isize, qnuma_del, qso4a_del, qnh4a_del,   &
               qh2so4_del, qnh3_del, dens_nh4so4a )

            print*, "k=",k
            if ((it == ite) .and. (jt == jte) .and. (mod(k,6)==1)) then
            print*, "qh2so4_cur=",qh2so4_cur," qnh3_cur=",qnh3_cur
            print*, "qh2so4_del=",qh2so4_del," qnh3_del=",qnh3_del
            print*, "qso4a_del=",qso4a_del," qnh4a_del=",qnh4a_del
            print*, "qnuma_del=",qnuma_del
            end if
        else if (newnuc_method .eq. 3) then  
            call glomap_nuc_mosaic_1box(   &
               dtnuc, temp_box, rh_box, cair_box,   &
               qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur, &
               qbiog1_c_avg, qbiog1_c_cur, qbiog2_c_avg, qbiog2_c_cur, &
               nsize, maxd_asize, volumlo_nuc, volumhi_nuc,   &
               isize, qnuma_del, qso4a_del, qnh4a_del, &
               qbiog1_ca_del, qbiog2_ca_del, &
               qh2so4_del, qnh3_del, &
               qbiog1_c_del, qbiog2_c_del, &
               dens_so4_aer, dens_nh4_aer, dens_biog1_c_aer, dens_biog2_c_aer, &
               mw_so4_aer, mw_nh4_aer, mw_biog1_c_aer, mw_biog2_c_aer, &
               CRII_G, CSINK_ION, &
               JBNRATE_1D(k), JBCHRATE_1D(k), JTNRATE_1D(k), JTCHRATE_1D(k), &
               JBNORGRATE_1D(k), JBCHORGRATE_1D(k), JTNRICRATE_1D(k), JRATE_1D(k))

            if (ktauc == 1) then
               JBNRATEAC_1D(k) = 0.0
               JBCHRATEAC_1D(k) = 0.0
               JTNRATEAC_1D(k) = 0.0
               JTCHRATEAC_1D(k) = 0.0
               JBNORGRATEAC_1D(k) = 0.0
               JBCHORGRATEAC_1D(k) = 0.0
               JTNRICRATEAC_1D(k) = 0.0
               JRATEAC_1D(k) = 0.0
            end if
            JBNRATEAC_1D(k) = JBNRATEAC_1D(k) + JBNRATE_1D(k)*dtnuc
            JBCHRATEAC_1D(k) = JBCHRATEAC_1D(k) + JBCHRATE_1D(k)*dtnuc
            JTNRATEAC_1D(k) = JTNRATEAC_1D(k) + JTNRATE_1D(k)*dtnuc
            JTCHRATEAC_1D(k) = JTCHRATEAC_1D(k) + JTCHRATE_1D(k)*dtnuc
            JBNORGRATEAC_1D(k) = JBNORGRATEAC_1D(k) + JBNORGRATE_1D(k)*dtnuc
            JBCHORGRATEAC_1D(k) = JBCHORGRATEAC_1D(k) + JBCHORGRATE_1D(k)*dtnuc
            JTNRICRATEAC_1D(k) = JTNRICRATEAC_1D(k) + JTNRICRATE_1D(k)*dtnuc
            JRATEAC_1D(k) = JRATEAC_1D(k) + JRATE_1D(k)*dtnuc

            print*, "ktauc=",ktauc
            print*, "k=",k
            if ((it == ite) .and. (jt == jte) .and. (mod(k,6)==1)) then
            print*, "qh2so4_cur=",qh2so4_cur," qnh3_cur=",qnh3_cur
            print*, "qbiog1_c_cur=",qbiog1_c_cur," qbiog2_c_cur=",qbiog2_c_cur
            print*, "qbiog1_ca_cur=",qbiog1_ca_cur," qbiog2_ca_cur=",qbiog2_ca_cur
            print*, "gasfracbiog1=",qbiog1_c_cur/(qbiog1_c_cur+qbiog1_ca_cur), &
                   "gasfracbiog2=",qbiog2_c_cur/(qbiog2_c_cur+qbiog2_ca_cur)
            print*, "qh2so4_del=",qh2so4_del," qnh3_del=",qnh3_del
            print*, "qbiog1_c_del=",qbiog1_c_del," qbiog2_c_del=",qbiog2_c_del
            print*, "qso4a_del=",qso4a_del," qnh4a_del=",qnh4a_del
            print*, "qbiog1_ca_del=",qbiog1_ca_del,"qbiog2_ca_del=",qbiog2_ca_del
            print*, "qnuma_del=",qnuma_del
            end if
	else
	    istat_newnuc = -1
	    return
	end if



	dumveca(1) = temp_box
	dumveca(2) = rh_box
	dumveca(3) = rsub(kso2,k,m)
	dumveca(4) = qh2so4_avg
	dumveca(5) = qnh3_avg
	dumveca(6) = qnuma_del
	do l = 1, 6
	    dumvecb(l) = min( dumvecb(l), dumveca(l) )
	    dumvecc(l) = max( dumvecc(l), dumveca(l) )
	    dumvecd(l) = dumvecd(l) + dumveca(l)
	    if (qnuma_del .gt. dumvece(6)) dumvece(l) = dumveca(l)
	end do



	if (qnuma_del .le. 0.0) goto 2700


	if (isize .ne. 1) ncount(3) = ncount(3) + 1
	if ((isize .lt. 1) .or. (isize .gt. nsize)) then
	    write(msg,93010) 'newnucxx bad isize_nuc' , it, jt, k,   &
		isize, nsize
	    call peg_message( lunerr, msg )
	    goto 2700
	end if
93010	format( a, 3i3, 1p, 9e10.2 )


	ncount(2) = ncount(2) + 1


	rsub(kh2so4,k,m) = max( 0.0, rsub(kh2so4,k,m) + qh2so4_del )
	rsub(knh3,  k,m) = max( 0.0, rsub(knh3,  k,m) + qnh3_del )

        if (newnuc_method == 3) then
        rsub(kbiog1_c,k,m) = max( 0.0, rsub(kbiog1_c,k,m) + qbiog1_c_del )
        rsub(kbiog2_c,k,m) = max( 0.0, rsub(kbiog2_c,k,m) + qbiog2_c_del )
        end if


	l = lptr_so4_aer(isize,itype,iphase)
	if (l .ge. p1st) then
	    rsub(l,k,m) = rsub(l,k,m) + qso4a_del
	end if
	l = lptr_nh4_aer(isize,itype,iphase)
	if (l .ge. p1st) then
	    rsub(l,k,m) = rsub(l,k,m) + qnh4a_del
	end if

        if (newnuc_method == 3) then
        l = lptr_biog1_c_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            rsub(l,k,m) = rsub(l,k,m) + qbiog1_ca_del
        end if
        l = lptr_biog2_c_aer(isize,itype,iphase)
        if (l .ge. p1st) then
            rsub(l,k,m) = rsub(l,k,m) + qbiog2_ca_del
        end if
        end if

	l = numptr_aer(isize,itype,iphase)
	rsub(l,k,m) = rsub(l,k,m) + qnuma_del
	xxnumb = rsub(l,k,m)





	l = waterptr_aer(isize,itype)
	if ((rh_box .gt. 0.10) .and. (l .ge. p1st)) then
	    aw = min( rh_box, 0.98 )
	    if (aw .lt. 0.97) then
		duma =       a_zsr_xx1 +   &
		        aw*( a_zsr_xx2 +   &
		        aw*( a_zsr_xx3 +   &
		        aw*( a_zsr_xx4 +   &
		        aw*( a_zsr_xx5 +   &
		        aw*  a_zsr_xx6 ))))
	    else
		dumb = -b_zsr_xx*log(aw)
	        dumb = max( dumb, 0.5 )
		duma = 1.0/(1.0 + 55.509/dumb)
	    end if
	    duma = max( duma, 0.01 )
	    dumc = (1.0 - duma)/duma
	    rsub(l,k,m) = rsub(l,k,m) + qso4a_del*dumc
	end if






	xxmass = aqmassdry_sub(isize,itype,k,m)
	xxdens = adrydens_sub( isize,itype,k,m)
	iconform_numb = 1

	if ((xxdens .lt. 0.1) .or. (xxdens .gt. 20.0)) then

	    continue
	else


	    xxvolu = xxmass/xxdens
	    duma = qso4a_del*mw_so4_aer + qnh4a_del*mw_nh4_aer + qbiog1_ca_del*mw_biog1_c_aer + qbiog2_ca_del*mw_biog2_c_aer 
	    xxmass = xxmass + duma
	    xxvolu  = xxvolu  + qso4a_del*mw_so4_aer/dens_so4_aer + qnh4a_del*mw_nh4_aer/dens_nh4_aer &
                    + qbiog1_ca_del*mw_biog1_c_aer/dens_biog1_c_aer + qbiog2_ca_del*mw_biog2_c_aer/dens_biog2_c_aer 
	    if (xxmass .le. smallmassbb) then

		xxdens = 0.001
	    else if (xxmass .gt. 1000.0*xxvolu) then


		xxdens = 1000.0
	    else 
		xxdens = xxmass/xxvolu
	    end if
	end if

	if ((xxdens .lt. 0.1) .or. (xxdens .gt. 20.0)) then


	    ncount(4) = ncount(4) + 1
	    xxmass = 0.0
	    xxvolu  = 0.0
	    do ll = 1, ncomp_aer(itype)
		l = massptr_aer(ll,isize,itype,iphase)
		if (l .ge. p1st) then
		    duma = max( 0.0, rsub(l,k,m) )*mw_aer(ll,itype)
		    xxmass = xxmass + duma
		    xxvolu = xxvolu + duma/dens_aer(ll,itype)
		end if
	    end do
	end if

	if (xxmass .le. smallmassbb) then


	    ncount(5) = ncount(5) + 1
	    xxdens = densdefault
	    xxvolu = xxmass/xxdens
	    xxnumb = xxmass/(volumcen_sect(isize,itype)*xxdens)
	    iconform_numb = 0
	    l = waterptr_aer(isize,itype)
	    if (l .ge. p1st) rsub(l,k,m) = 0.0
	    l = hyswptr_aer(isize,itype)
	    if (l .ge. p1st) rsub(l,k,m) = 0.0
	else
	    xxdens = xxmass/xxvolu
	end if

	if (iconform_numb .gt. 0) then

	    if (xxnumb .gt. xxvolu/volumlo_sect(isize,itype)) then
		ncount(6) = ncount(6) + 1
		xxnumb = xxvolu/volumlo_sect(isize,itype)
	    else if (xxnumb .lt. xxvolu/volumhi_sect(isize,itype)) then
		ncount(7) = ncount(7) + 1
		xxnumb = xxvolu/volumhi_sect(isize,itype)
	    end if
	end if


	l = numptr_aer(isize,itype,iphase)
	rsub(l,k,m) = xxnumb
	adrydens_sub( isize,itype,k,m) = xxdens
	aqmassdry_sub(isize,itype,k,m) = xxmass
	aqvoldry_sub( isize,itype,k,m) = xxvolu


2700	continue


	if ((idiagbb .ge. 100) .and.   &
	    (it .eq. ite) .and.    &
	    (jt .eq. jte) .and. (k .eq. kclm_calcend)) then
	    if (idiagbb .ge. 110) then
	      write(msg,93020) 'newnucbb mins ', dumvecb(1:6)
	      call peg_message( lunerr, msg )
	      write(msg,93020) 'newnucbb maxs ', dumvecc(1:6)
	      call peg_message( lunerr, msg )
	      duma = max( 1, ncount(1) ) 
	      write(msg,93020) 'newnucbb avgs ', dumvecd(1:6)/duma
	      call peg_message( lunerr, msg )
	      write(msg,93020) 'newnucbb hinuc', dumvece(1:6)
	      call peg_message( lunerr, msg )
	      write(msg,93020) 'newnucbb dtnuc', dtnuc
	      call peg_message( lunerr, msg )
	    end if
	    write(msg,93030) 'newnucbb ncnt ', ncount(1:7)
	    call peg_message( lunerr, msg )
	end if
93020	format( a, 1p, 10e10.2 )
93030	format( a, 1p, 10i10 )


2800	continue	

2900	continue	


	return
	end subroutine mosaic_newnuc_1clm






        subroutine ternary_nuc_mosaic_1box(   &
           dtnuc, temp_in, rh_in, cair,   &
           qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
           nsize, maxd_asize, volumlo_sect, volumhi_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a )
























        implicit none


        real, intent(in) :: dtnuc             
        real, intent(in) :: temp_in           
        real, intent(in) :: rh_in             
        real, intent(in) :: cair              

        real, intent(in) :: qh2so4_avg, qh2so4_cur   
        real, intent(in) :: qnh3_avg, qnh3_cur       
             
             

        integer, intent(in) :: nsize                    
        integer, intent(in) :: maxd_asize               
        real, intent(in) :: volumlo_sect(maxd_asize)    
        real, intent(in) :: volumhi_sect(maxd_asize)    


        integer, intent(out) :: isize_nuc     
        real, intent(out) :: qnuma_del        
        real, intent(out) :: qso4a_del        
        real, intent(out) :: qnh4a_del        
        real, intent(out) :: qh2so4_del       
        real, intent(out) :: qnh3_del         
                                              


        real, intent(inout) :: dens_nh4so4a   





        real :: ratenuclt        
        real :: rateloge         
        real :: cnum_h2so4       
        real :: cnum_nh3         
        real :: cnum_tot         
        real :: radius_cluster   







        integer i
        integer, save :: icase = 0, icase_reldiffmax = 0

        real, parameter :: pi = 3.1415926536
        real, parameter :: avogad = 6.022e23   
        real, parameter :: mw_air = 28.966     



        real, parameter :: dens_ammsulf   = 1.769
        real, parameter :: dens_ammbisulf = 1.78
        real, parameter :: dens_sulfacid  = 1.841




        real, parameter :: mw_ammsulf   = 132.0
        real, parameter :: mw_ammbisulf = 114.0
        real, parameter :: mw_sulfacid  =  96.0

        real, parameter :: mw_so4a      =  96.0
        real, parameter :: mw_nh4a      =  18.0

        real, save :: reldiffmax = 0.0

        real dens_part                
        real duma, dumb, dumc, dume
        real dum_m1, dum_m2, dum_m3, dum_n1, dum_n2, dum_n3
        real fogas, foso4a, fonh4a, fonuma
        real freduce                  
                                      
        real freducea, freduceb
        real gramaero_per_moleso4a    
        real mass_part                
        real molenh4a_per_moleso4a    
        real nh3conc_in               
        real so4vol_in                
        real qmolnh4a_del_max         
        real qmolso4a_del_max         
        real vol_cluster              
        real vol_part                 






        isize_nuc = 1
        qnuma_del = 0.0
        qso4a_del = 0.0
        qnh4a_del = 0.0
        qh2so4_del = 0.0
        qnh3_del = 0.0
        if (qh2so4_avg .le. 4.0e-16) return
        if (qh2so4_cur .le. 4.0e-16) return








        nh3conc_in = qnh3_avg/1.0e-12



        so4vol_in  = (qh2so4_avg) * cair * avogad

        call ternary_nuc_napari(   &
            temp_in, rh_in, nh3conc_in, so4vol_in,   &
            ratenuclt, rateloge,   &
            cnum_h2so4, cnum_nh3, cnum_tot, radius_cluster )



        if (ratenuclt .le. 1.0e-6) return




        vol_cluster = (pi*4.0/3.0)* (radius_cluster**3) * 1.0e-21
        isize_nuc = 1
        vol_part = volumlo_sect(1)
        if (vol_cluster .le. volumlo_sect(1)) then
           continue
        else if (vol_cluster .ge. volumhi_sect(nsize)) then
           isize_nuc = nsize
           vol_part = volumhi_sect(nsize)
        else
           do i = 1, nsize
              if (vol_cluster .lt. volumhi_sect(i)) then
                 isize_nuc = i
                 vol_part = vol_cluster
                 vol_part = min( vol_part, volumhi_sect(i) )
                 vol_part = max( vol_part, volumlo_sect(i) )
                 exit
              end if
           end do
        end if










        if (qnh3_cur .ge. qh2so4_cur) then


           dum_n1 = (qnh3_cur/qh2so4_cur) - 1.0
           dum_n1 = max( 0.0, min( 1.0, dum_n1 ) )
           dum_n2 = 1.0 - dum_n1
           dum_n3 = 0.0
        else


           dum_n1 = 0.0
           dum_n2 = (qnh3_cur/qh2so4_cur)
           dum_n2 = max( 0.0, min( 1.0, dum_n2 ) )
           dum_n3 = 1.0 - dum_n2
	end if

        dum_m1 = dum_n1*mw_ammsulf
        dum_m2 = dum_n2*mw_ammbisulf
        dum_m3 = dum_n3*mw_sulfacid
        dens_part = (dum_m1 + dum_m2 + dum_m3)/   &
           ((dum_m1/dens_ammsulf) + (dum_m2/dens_ammbisulf)   &
                                  + (dum_m3/dens_sulfacid))

	if (abs(dens_nh4so4a-1.8) .le. 0.2) then
	    dens_part = dens_nh4so4a
	else
            dens_nh4so4a = dens_part
	end if
        mass_part  = vol_part*dens_part 
        molenh4a_per_moleso4a = 2.0*dum_n1 + dum_n2  
        gramaero_per_moleso4a = dum_m1 + dum_m2 + dum_m3  



        duma = max( 0.0, (ratenuclt*dtnuc*mass_part) )

        dumc = duma/gramaero_per_moleso4a

        dume = dumc/cair


        qmolso4a_del_max = dume


        freducea = 1.0
        if (qmolso4a_del_max .gt. qh2so4_cur) then
           freducea = qh2so4_cur/qmolso4a_del_max
        end if


        freduceb = 1.0
        if (molenh4a_per_moleso4a .ge. 1.0e-10) then

           qmolnh4a_del_max = qmolso4a_del_max*molenh4a_per_moleso4a
           if (qmolnh4a_del_max .gt. qnh3_cur) then
              freduceb = qnh3_cur/qmolnh4a_del_max
           end if
        end if
        freduce = min( freducea, freduceb )



        if (freduce*ratenuclt .le. 1.0e-6) return














        duma = 0.9999
        qh2so4_del = min( duma*qh2so4_cur, freduce*qmolso4a_del_max )
        qnh3_del   = min( duma*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
        qh2so4_del = -qh2so4_del
        qnh3_del   = -qnh3_del


        qso4a_del = -qh2so4_del
        qnh4a_del =   -qnh3_del

        qnuma_del = (qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part


        return
        end subroutine ternary_nuc_mosaic_1box




        subroutine ternary_nuc_napari(   &
           temp_in, rh_in, nh3conc_in, so4vol_in,   &
           ratenuclt, rateloge,   &
           cnum_h2so4, cnum_nh3, cnum_tot, radius_cluster )












































      implicit none


        real, intent(in) :: temp_in           
        real, intent(in) :: rh_in             
        real, intent(in) :: nh3conc_in        
        real, intent(in) :: so4vol_in         


        real, intent(out) :: ratenuclt        
        real, intent(out) :: rateloge         

        real, intent(out) :: cnum_h2so4       
                                              
        real, intent(out) :: cnum_nh3         
                                              
        real, intent(out) :: cnum_tot         
                                              
        real, intent(out) :: radius_cluster   


        integer ncoeff
        parameter ( ncoeff = 4 )      
                                      
        integer npoly
        parameter ( npoly = 20 )      

        integer n                     


        real f ( npoly )


        real a  ( ncoeff, npoly )

        real temp, rh, nh3conc, so4vol     
        real log_rh, log_nh3conc, log_so4vol


        data a  / -0.355297,  -33.8449,      0.34536,    -8.24007e-4,   &
                   3.13735,    -0.772861,    5.61204e-3, -9.74576e-6,   &
                  19.0359,     -0.170957,    4.79808e-4, -4.14699e-7,   &
                   1.07605,     1.48932,    -7.96052e-3,  7.61229e-6,   &
                   6.0916,     -1.25378,     9.39836e-3, -1.74927e-5,   &
                   0.31176,     1.64009,    -3.43852e-3, -1.09753e-5,   &
                  -2.00738e-2, -0.752115,    5.25813e-3, -8.98038e-6,   &
                   0.165536,    3.26623,    -4.89703e-2,  1.46967e-4,   &
                   6.52645,    -0.258002,    1.43456e-3, -2.02036e-6,   &
                   3.68024,    -0.204098,    1.06259e-3, -1.2656e-6 ,   &
                  -6.6514e-2,  -7.82382,     1.22938e-2,  6.18554e-5,   &
                   0.65874,     0.190542,   -1.65718e-3,  3.41744e-6,   &
                   5.99321e-2,  5.96475,    -3.62432e-2,  4.93337e-5,   &
                  -0.732731,   -1.84179e-2,  1.47186e-4, -2.37711e-7,   &
                   0.728429,    3.64736,    -2.7422e-2,   4.93478e-5,   &
                  41.3016,     -0.35752,     9.04383e-4, -5.73788e-7,   &
                  -0.160336,    8.89881e-3, -5.39514e-5,  8.39522e-8,   &
                   8.57868,    -0.112358,    4.72626e-4, -6.48365e-7,   &
                   5.30167e-2, -1.98815,     1.57827e-2, -2.93564e-5,   &
                  -2.32736,     2.34646e-2, -7.6519e-5,   8.0459e-8 /





        temp    = max( 240.15, min (300.15, temp_in ) )
        rh      = max( 0.05,   min (0.95,   rh_in ) )
        so4vol  = max( 1.0e4,  min (1.0e9,  so4vol_in ) )
        nh3conc = max( 0.1,    min (100.0,  nh3conc_in ) )





        do n = 1, npoly
            f ( n )   = a ( 1, n ) + a ( 2, n ) * temp   &
                      + a ( 3, n ) * ( temp ) ** 2.0   &
                      + a ( 4, n ) * ( temp ) ** 3.0
        end do




        log_rh = log ( rh )
        log_nh3conc = log ( nh3conc )
        log_so4vol = log ( so4vol )
        rateloge = -84.7551   &
                 + f ( 1 ) / log_so4vol   &
                 + f ( 2 )  * ( log_so4vol )   &
                 + f ( 3 )  * ( log_so4vol ) **2.0   &
                 + f ( 4 )  * ( log_nh3conc )   &
                 + f ( 5 )  * ( log_nh3conc ) **2.0   &
                 + f ( 6 )  * rh   &
                 + f ( 7 )  * ( log_rh )   &
                 + f ( 8 )  * ( log_nh3conc /   &
                              log_so4vol )   &
                 + f ( 9 )  * ( log_nh3conc *   &
                              log_so4vol )   &
                 + f ( 10 ) * rh  *   &
                              ( log_so4vol )   &
                 + f ( 11 ) * ( rh /   &
                              log_so4vol )   &
                 + f ( 12 ) * ( rh *   &
                              log_nh3conc )   &
                 + f ( 13 ) * ( log_rh /   &
                              log_so4vol )   &
                 + f ( 14 ) * ( log_rh *   &
                              log_nh3conc )   &
                 + f ( 15 ) * (( log_nh3conc ) ** 2.0   &
                              / log_so4vol )   &
                 + f ( 16 ) * ( log_so4vol *   &
                              ( log_nh3conc ) ** 2.0 )   &
                 + f ( 17 ) * (( log_so4vol ) ** 2.0 *   &
                              log_nh3conc )   &
                 + f ( 18 ) * ( rh  *   &
                              ( log_nh3conc ) ** 2.0 )   &
                 + f ( 19 ) * ( rh  *  log_nh3conc   &
                              / log_so4vol )   &
                 + f ( 20 ) * (( log_so4vol ) ** 2.0 *   &
                              ( log_nh3conc ) ** 2.0 )

        ratenuclt = exp ( rateloge )



        cnum_h2so4 = 38.1645 + 0.774106 * rateloge   &
                   + 0.00298879 * ( rateloge ) ** 2.0   &
                   - 0.357605 * temp   &
                   - 0.00366358 * temp * rateloge   &
                   + 0.0008553 * ( temp ) ** 2.0

        cnum_nh3   = 26.8982 + 0.682905 * rateloge   &
                   + 0.00357521 * ( rateloge ) ** 2.0   &
                   - 0.265748 * temp   &
                   - 0.00341895 * temp * rateloge   &
                   + 0.000673454 * ( temp ) ** 2.0

        cnum_tot   = 79.3484 + 1.7384 * rateloge   &
                   + 0.00711403 * ( rateloge ) ** 2.0   &
                   - 0.744993 * temp   &
                   - 0.00820608 * temp * rateloge   &
                   + 0.0017855 * ( temp ) ** 2.0

        radius_cluster = 0.141027 - 0.00122625 * rateloge   &
                   - 7.82211e-6 * ( rateloge ) ** 2.0   &
                   - 0.00156727 * temp   &
                   - 0.00003076 * temp * rateloge   &
                   + 0.0000108375 * ( temp ) ** 2.0

        return
        end subroutine ternary_nuc_napari




        subroutine wexler_nuc_mosaic_1box(   &
           dtnuc, temp_in, rh_in, cair,   &
           qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur,   &
           nsize, maxd_asize, volumlo_sect, volumhi_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del,   &
           qh2so4_del, qnh3_del, dens_nh4so4a )






















        implicit none


        real, intent(in) :: dtnuc             
        real, intent(in) :: temp_in           
        real, intent(in) :: rh_in             
        real, intent(in) :: cair              

        real, intent(in) :: qh2so4_avg, qh2so4_cur   
        real, intent(in) :: qnh3_avg, qnh3_cur       
             
             

        integer, intent(in) :: nsize                    
        integer, intent(in) :: maxd_asize               
        real, intent(in) :: volumlo_sect(maxd_asize)    
        real, intent(in) :: volumhi_sect(maxd_asize)    


        integer, intent(out) :: isize_nuc     
        real, intent(out) :: qnuma_del        
        real, intent(out) :: qso4a_del        
        real, intent(out) :: qnh4a_del        
        real, intent(out) :: qh2so4_del       
        real, intent(out) :: qnh3_del         
                                              


        real, intent(inout) :: dens_nh4so4a   



        integer i
        integer, save :: icase = 0, icase_reldiffmax = 0

        real, parameter :: pi = 3.1415926536
        real, parameter :: avogad = 6.022e23   
        real, parameter :: mw_air = 28.966     



        real, parameter :: dens_ammsulf   = 1.769
        real, parameter :: dens_ammbisulf = 1.78
        real, parameter :: dens_sulfacid  = 1.841




        real, parameter :: mw_ammsulf   = 132.0
        real, parameter :: mw_ammbisulf = 114.0
        real, parameter :: mw_sulfacid  =  96.0

        real, parameter :: mw_so4a      =  96.0
        real, parameter :: mw_nh4a      =  18.0

        real, save :: reldiffmax = 0.0

        real ch2so4_crit              
        real dens_part                
        real duma, dumb, dumc, dume
        real dum_m1, dum_m2, dum_m3, dum_n1, dum_n2, dum_n3
        real fogas, foso4a, fonh4a, fonuma
        real mass_part                
        real molenh4a_per_moleso4a    
        real qh2so4_crit              
        real qh2so4_avail             
        real vol_part                 
        real diameter_part            




        isize_nuc = 1
        qnuma_del = 0.0
        qso4a_del = 0.0
        qnh4a_del = 0.0
        qh2so4_del = 0.0
        qnh3_del = 0.0




        ch2so4_crit = 0.16 * exp( 0.1*temp_in - 3.5*rh_in - 27.7 )


        qh2so4_crit = (ch2so4_crit*1.0e-12/98.0)/cair
        qh2so4_avail = (qh2so4_cur - qh2so4_crit)*dtnuc



        if (qh2so4_avail .le. 4.0e-18) then
           return
        end if



        if (nsize == 20) then
           isize_nuc = 2
           diameter_part = 1.7 
           vol_part = (pi/6.0)* (diameter_part**3) * 1.0e-21  
        else if (nsize == 12) then
           isize_nuc = 1
           diameter_part = 3.0 
           vol_part = (pi/6.0)* (diameter_part**3) * 1.0e-21  
        else
           isize_nuc = 1
           vol_part = volumlo_sect(1)
        end if









        if (qnh3_cur .ge. qh2so4_avail) then


           dum_n1 = (qnh3_cur/qh2so4_avail) - 1.0
           dum_n1 = max( 0.0, min( 1.0, dum_n1 ) )
           dum_n2 = 1.0 - dum_n1
           dum_n3 = 0.0
        else


           dum_n1 = 0.0
           dum_n2 = (qnh3_cur/qh2so4_avail)
           dum_n2 = max( 0.0, min( 1.0, dum_n2 ) )
           dum_n3 = 1.0 - dum_n2
	end if

        dum_m1 = dum_n1*mw_ammsulf
        dum_m2 = dum_n2*mw_ammbisulf
        dum_m3 = dum_n3*mw_sulfacid
        dens_part = (dum_m1 + dum_m2 + dum_m3)/   &
           ((dum_m1/dens_ammsulf) + (dum_m2/dens_ammbisulf)   &
                                  + (dum_m3/dens_sulfacid))

	if (abs(dens_nh4so4a-1.8) .le. 0.2) then
	    dens_part = dens_nh4so4a
	else
            dens_nh4so4a = dens_part
	end if
        mass_part  = vol_part*dens_part 
        molenh4a_per_moleso4a = 2.0*dum_n1 + dum_n2



        duma = 0.9999
        qh2so4_del = min( duma*qh2so4_cur, qh2so4_avail )
        qnh3_del   = min( duma*qnh3_cur, qh2so4_del*molenh4a_per_moleso4a )
        qh2so4_del = -qh2so4_del
        qnh3_del   = -qnh3_del


        qso4a_del = -qh2so4_del
        qnh4a_del =   -qnh3_del

        qnuma_del = (qso4a_del*mw_so4a + qnh4a_del*mw_nh4a)/mass_part


        return
        end subroutine wexler_nuc_mosaic_1box



        subroutine glomap_nuc_mosaic_1box(   &
           dtnuc, temp_in, rh_in, cair,   &
           qh2so4_avg, qh2so4_cur, qnh3_avg, qnh3_cur, &
           qbiog1_c_avg, qbiog1_c_cur, qbiog2_c_avg, qbiog2_c_cur, &
           nsize, maxd_asize, volumlo_sect, volumhi_sect,   &
           isize_nuc, qnuma_del, qso4a_del, qnh4a_del, &
           qbiog1_ca_del, qbiog2_ca_del, &
           qh2so4_del, qnh3_del, & 
           qbiog1_c_del, qbiog2_c_del, &
           dens_so4_aer, dens_nh4_aer, dens_biog1_c_aer, dens_biog2_c_aer, &
           mw_so4_aer, mw_nh4_aer, mw_biog1_c_aer, mw_biog2_c_aer, &
           CRII_G, CSINK_ION, &
           J1POINT7BNRATE, J1POINT7BCHRATE, J1POINT7TNRATE, J1POINT7TCHRATE, &
           J1POINT7BNORGRATE, J1POINT7BCHORGRATE, J1POINT7TNRICRATE, J1POINT7RATE)
























        implicit none


        real, intent(in) :: dtnuc             
        real, intent(in) :: temp_in           
        real, intent(in) :: rh_in             
        real, intent(in) :: cair              

        real, intent(in) :: qh2so4_avg, qh2so4_cur   
        real, intent(in) :: qnh3_avg, qnh3_cur       
        real, intent(in) :: qbiog1_c_avg, qbiog1_c_cur 
        real, intent(in) :: qbiog2_c_avg, qbiog2_c_cur 
             
             

        integer, intent(in) :: nsize                    
        integer, intent(in) :: maxd_asize               
        real, intent(in) :: volumlo_sect(maxd_asize)    
        real, intent(in) :: volumhi_sect(maxd_asize)    

        real, intent(in) :: dens_so4_aer   
        real, intent(in) :: dens_nh4_aer   
        real, intent(in) :: dens_biog1_c_aer 
        real, intent(in) :: dens_biog2_c_aer 

        real, intent(in) :: mw_so4_aer   
        real, intent(in) :: mw_nh4_aer   
        real, intent(in) :: mw_biog1_c_aer   
        real, intent(in) :: mw_biog2_c_aer   

        real, intent(in) :: CRII_G  
        real, intent(in) :: CSINK_ION 


        integer, intent(out) :: isize_nuc     
        real, intent(out) :: qnuma_del        
        real, intent(out) :: qso4a_del        
        real, intent(out) :: qnh4a_del        
        real, intent(out) :: qh2so4_del       
        real, intent(out) :: qnh3_del         
        real, intent(out) :: qbiog1_ca_del    
        real, intent(out) :: qbiog2_ca_del    
        real, intent(out) :: qbiog1_c_del     
        real, intent(out) :: qbiog2_c_del     
                                              

        real, intent(inout) :: J1POINT7RATE         
        real, intent(inout) :: J1POINT7BNRATE       
        real, intent(inout) :: J1POINT7BCHRATE      
        real, intent(inout) :: J1POINT7TNRATE       
        real, intent(inout) :: J1POINT7TCHRATE      
        real, intent(inout) :: J1POINT7BNORGRATE    
        real, intent(inout) :: J1POINT7BCHORGRATE   
        real, intent(inout) :: J1POINT7TNRICRATE    


        real, parameter :: pi = 3.1415926536
        real, parameter :: avogad = 6.022e23   
        real, parameter :: mw_air = 28.966     

        real, parameter :: SO_YIELD_NUC =  0.0145 
        real, parameter :: SO_YIELD_NUCOH = 0.006 
        real, parameter :: CONC_EPS = 1.0e2       

        real dens_part                
        real mass_part                
        real vol_part                 
        real duma, dumb, dumc, dume
        real dum_n1, dum_n2, dum_n3
        real molenh4a_per_moleso4a    
        real gramaero_per_moleso4a    
        real qh2so4_avail             
        real diameter_part            
        real mass_part_so4            
        real mass_part_so4org         
        real mass_part_org            
        real qmolso4a_del_max         
        real qmolorga_del_max         
        real freduce1,freduce2,freduce3,freduce 


        REAL :: H2SO4                
        REAL :: NH3                  
        REAL :: BIOG1_C              
        REAL :: BIOG2_C              
        REAL :: SECORGFROMOHO3       
        REAL :: SECORGFROMOH         
        REAL :: TEMPT                
        REAL :: RH                   
        REAL :: AIRD                 

        REAL :: TEMP                 
        REAL :: RHF                  
        double precision NH3_SCALE            
        double precision SA_SCALE             
        REAL :: KINLIM               
        REAL :: ION_CONC             
        REAL :: ALFA                 
        REAL :: org_bi_ch
        double precision X, nh3_coeff_n,nh3_coeff_ch
        double precision k_bn,k_tn,k_bch,k_tch,sa_bn,sa_bi,sa_tn,sa_ti
        REAL :: ftnuc                
        REAL :: fthom                
        double precision A(20)


        qnuma_del = 0.0
        qso4a_del = 0.0
        qnh4a_del = 0.0
        qh2so4_del = 0.0
        qnh3_del = 0.0
        qbiog1_ca_del = 0.0
        qbiog1_c_del = 0.0
        qbiog2_ca_del = 0.0
        qbiog2_c_del = 0.0


        J1POINT7RATE=0
        J1POINT7BNRATE=0
        J1POINT7BNORGRATE=0
        J1POINT7BCHRATE=0
        J1POINT7BCHORGRATE=0
        J1POINT7TNRATE=0
        J1POINT7TNRICRATE=0
        J1POINT7TCHRATE=0


        TEMPT = temp_in
        RH = rh_in
        AIRD = cair * avogad 
        H2SO4 = qh2so4_cur * cair * avogad  
        NH3 = qnh3_cur * cair * avogad  
        BIOG1_C = qbiog1_c_cur * cair * avogad  
        BIOG2_C = qbiog2_c_cur * cair * avogad  
        SECORGFROMOH = BIOG1_C/0.2065 
        SECORGFROMOHO3 = BIOG1_C/0.2065*SO_YIELD_NUCOH + BIOG2_C/0.0545*SO_YIELD_NUC
                 
                 



        IF (H2SO4.LE.CONC_EPS .AND. SECORGFROMOHO3 .LE. CONC_EPS) RETURN


        ALFA = 6e-8*sqrt(300/TEMPT) + 6e-26*AIRD*(300/TEMPT)**4.0 
        IF (TEMPT .LT. 208) TEMPT = 208  
        RHF = (RH/0.38)*(1+0.02*(RH-0.38)*(TEMPT-208)**1.2)
        IF (RHF .LT. 0) RHF=0.0
        NH3_SCALE = NH3*1E-6   
        SA_SCALE = 1.e-6 * H2SO4
        TEMP = TEMPT / 1000.



        A = (/3.95451,9.702973,12.62259,-0.007066146,182.4495,1.203451, &
           -4.188065,9.003471,636801.6,2.891024,-11.48166,25.49469,0.1810722,3.373738,4.071246, &
           -23.8002,37.03029,0.227413,206.9792,3.138719 /)

        nh3_coeff_n = A(9) * (NH3_SCALE**(1-A(8)))                         &
        / (1 + A(9) * (NH3_SCALE**(1-A(8))) * (SA_SCALE**A(10)) )

        nh3_coeff_ch =A(19) * (NH3_SCALE**(1-A(15)))                       &
        / (1 + A(19) * (NH3_SCALE**(1-A(15))) * (SA_SCALE**A(20)) )

        k_bn = exp(A(2) - exp(A(3)*(TEMP-A(4))))

        k_bch = exp(A(11) - exp(A(12)*(TEMP-A(13))))

        k_tn =(exp(A(5) - exp(A(6)*(TEMP-A(7))))) * nh3_coeff_n *         &
        (NH3_SCALE / 100.)**A(8) * 100**A(8)

        k_tch =(exp(A(16)-exp(A(17)*(TEMP-A(18)))))*nh3_coeff_ch*(NH3_SCALE**A(15))

        sa_bn = SA_SCALE**A(1)
        sa_bi = SA_SCALE**A(14)
        sa_tn = SA_SCALE**A(10)
        sa_ti = SA_SCALE**A(20) 

        org_bi_ch = 0.00136641*(SECORGFROMOHO3/1.0e7)**(1.56588+0.186303/(SECORGFROMOHO3/1.0e7))

        ftnuc = EXP(-(TEMPT-278.0)/10.0)        

        fthom = 1.0

        X = CSINK_ION + sa_bi * k_bch + k_tch * sa_ti + org_bi_ch

        ION_CONC = (SQRT(X**2 + 4. * ALFA * CRII_G) - X) / (2 * ALFA)

        IF (ION_CONC .EQ. 0.) THEN
           IF (X*X .GT. CRII_G*ALFA) THEN
              ION_CONC = CRII_G / x
           ELSE
              ION_CONC = SQRT(CRII_G / ALFA)
           ENDIF
        ENDIF

        J1POINT7BNRATE  = k_bn * sa_bn*rhf
        J1POINT7BCHRATE = ION_CONC * sa_bi * k_bch*rhf
        J1POINT7TNRATE  = k_tn * sa_tn*rhf
        J1POINT7TCHRATE = k_tch * ION_CONC * sa_ti*rhf
        
        J1POINT7BNORGRATE  = 0.0400097*(SECORGFROMOHO3/1.0e7)**(1.84826+0.186303/(SECORGFROMOHO3/1.0e7)) *ftnuc*fthom  
        J1POINT7BCHORGRATE = ION_CONC*org_bi_ch *ftnuc*fthom  
        
        J1POINT7TNRICRATE  = 0.5*3.27E-21*H2SO4*H2SO4*SECORGFROMOH *ftnuc*fthom 
        
        KINLIM = (49.19*SA_SCALE**2.016)/(1+3.126*SA_SCALE**-3.594)
        IF(J1POINT7TNRATE .GT. KINLIM) THEN
           J1POINT7TNRATE = KINLIM
        ENDIF
        IF(J1POINT7TCHRATE .GT. KINLIM) THEN
           J1POINT7TCHRATE = KINLIM
        ENDIF


        if (nsize == 20) then    
           isize_nuc = 2
           diameter_part = 1.7 
           vol_part = (pi/6.0)* (diameter_part**3) * 1.0e-21  
        else if (nsize == 12) then
           isize_nuc = 1
           diameter_part = 3.0 
           vol_part = (pi/6.0)* (diameter_part**3) * 1.0e-21  
           
           J1POINT7BNRATE = J1POINT7BNRATE/5.5
           J1POINT7BCHRATE = J1POINT7BCHRATE/5.5
           J1POINT7TNRATE = J1POINT7TNRATE/5.5
           J1POINT7TCHRATE = J1POINT7TCHRATE/5.5
           J1POINT7BNORGRATE = J1POINT7BNORGRATE/5.5
           J1POINT7BCHORGRATE = J1POINT7BCHORGRATE/5.5
           J1POINT7TNRICRATE = J1POINT7TNRICRATE/5.5
        else
           call peg_error_fatal( -1,   &
               'The glomap nucleation option only works for 20-bin and 12-bin schemes' )
        end if


        mass_part_so4 = vol_part*(dens_so4_aer+dens_nh4_aer)/2
        mass_part_so4org = vol_part*(dens_so4_aer/2+dens_biog1_c_aer/4+dens_biog2_c_aer/4)
        mass_part_org = vol_part*(dens_biog1_c_aer+dens_biog2_c_aer)/2


        duma = max( 0.0, (J1POINT7BNORGRATE + J1POINT7BCHORGRATE)*dtnuc*mass_part_org/(mw_biog1_c_aer/2+mw_biog2_c_aer/2)/cair )
        dumb = max( 0.0, J1POINT7TNRICRATE*dtnuc*mass_part_so4org/(mw_so4_aer+mw_biog1_c_aer/2+mw_biog2_c_aer/2)/cair )
        qmolorga_del_max = duma+dumb
        qmolso4a_del_max = dumb


        dume = 0.9999
        freduce1 = 1.0
        if (qmolso4a_del_max .gt. qh2so4_cur*dume) then
           freduce1 = qh2so4_cur*dume/qmolso4a_del_max
        end if
        freduce2 = 1.0
        if (qmolorga_del_max .gt. (qbiog1_c_cur+qbiog2_c_cur)*dume) then
           freduce2 = (qbiog1_c_cur+qbiog2_c_cur)*dume/qmolorga_del_max
        end if
        freduce=min(freduce1,freduce2)


        qh2so4_del = freduce*qmolso4a_del_max
        qbiog1_c_del = (freduce2*duma+freduce*dumb) *qbiog1_c_cur/(qbiog1_c_cur+qbiog2_c_cur)
        qbiog2_c_del = (freduce2*duma+freduce*dumb) *qbiog2_c_cur/(qbiog1_c_cur+qbiog2_c_cur)


        qh2so4_avail = qh2so4_cur - qh2so4_del 
        if (qh2so4_avail .lt. 0.0) then
           call peg_error_fatal( -1,   &
               'fatal error: qh2so4_avail less than zero' )
        end if







        if (qnh3_cur .ge. qh2so4_avail) then


           dum_n1 = (qnh3_cur/qh2so4_avail) - 1.0
           dum_n1 = max( 0.0, min( 1.0, dum_n1 ) )
           dum_n2 = 1.0 - dum_n1
           dum_n3 = 0.0
        else


           dum_n1 = 0.0
           dum_n2 = (qnh3_cur/qh2so4_avail)
           dum_n2 = max( 0.0, min( 1.0, dum_n2 ) )
           dum_n3 = 1.0 - dum_n2
	end if

        molenh4a_per_moleso4a = 2.0*dum_n1 + dum_n2  
        gramaero_per_moleso4a = mw_so4_aer + molenh4a_per_moleso4a*mw_nh4_aer  


        dumc = max( 0.0, (J1POINT7BNRATE +J1POINT7BCHRATE+ J1POINT7TNRATE +J1POINT7TCHRATE) &
              *dtnuc*mass_part_so4/gramaero_per_moleso4a/cair )
        qmolso4a_del_max = dumc


        freduce3 = 1.0
        if (qmolso4a_del_max .gt. qh2so4_avail*dume) then
           freduce3 = qh2so4_avail*dume/qmolso4a_del_max
        end if


        qh2so4_del = qh2so4_del + freduce3*qmolso4a_del_max
        qnh3_del   = min( dume*qnh3_cur, freduce3*qmolso4a_del_max*molenh4a_per_moleso4a )


        qh2so4_del = -qh2so4_del
        qnh3_del   = -qnh3_del
        qbiog1_c_del = -qbiog1_c_del
        qbiog2_c_del = -qbiog2_c_del


        qso4a_del = -qh2so4_del
        qnh4a_del =   -qnh3_del
        qbiog1_ca_del = -qbiog1_c_del
        qbiog2_ca_del = -qbiog2_c_del


        qnuma_del = (J1POINT7BNORGRATE + J1POINT7BCHORGRATE)*dtnuc/cair*freduce2 &
                  + J1POINT7TNRICRATE*dtnuc/cair*freduce &
                  + (J1POINT7BNRATE +J1POINT7BCHRATE+ J1POINT7TNRATE +J1POINT7TCHRATE)*dtnuc/cair*freduce3




        J1POINT7BNRATE = J1POINT7BNRATE*freduce3
        J1POINT7BCHRATE = J1POINT7BCHRATE*freduce3
        J1POINT7TNRATE = J1POINT7TNRATE*freduce3
        J1POINT7TCHRATE = J1POINT7TCHRATE*freduce3
        J1POINT7BNORGRATE = J1POINT7BNORGRATE*freduce2
        J1POINT7BCHORGRATE = J1POINT7BCHORGRATE*freduce2
        J1POINT7TNRICRATE = J1POINT7TNRICRATE*freduce

        J1POINT7RATE= J1POINT7BNRATE +J1POINT7BCHRATE+    &
        J1POINT7TNRATE +J1POINT7TCHRATE + J1POINT7BNORGRATE+ &
        J1POINT7BCHORGRATE +  J1POINT7TNRICRATE

        return
        end subroutine glomap_nuc_mosaic_1box






	end module module_mosaic_newnuc
