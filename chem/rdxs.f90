


























      module module_xsections

      use module_params

      IMPLICIT NONE

      public :: o3xs, rdo2xs, rdso2xs, no2xs_jpl06a
      public :: rdxs_init
      private

      REAL, allocatable :: rei218(:), rei228(:), rei243(:), rei295(:)
      REAL :: v195, v345, v830
      REAL, allocatable :: wmo203(:), wmo273(:)
      REAL :: v176, v850

      REAL, allocatable :: jpl295(:), jpl218(:)
      REAL :: v186, v825

      REAL, allocatable :: mol226(:), mol263(:), mol298(:)
      REAL :: v185, v240, v350

      REAL, allocatable :: c0(:), c1(:), c2(:)
      REAL vb245, vb342

      REAL, allocatable :: no2xs_a(:), no2xs_b(:)

      CONTAINS

      SUBROUTINE rdxs_init( nw, wl )

      integer, intent(in) :: nw
      real, intent(in)    :: wl(nw)

      integer :: istat, astat

      istat = 0
      if( .not. allocated( rei218 ) ) then
        allocate( rei218(nw-1),rei228(nw-1),rei243(nw-1),rei295(nw-1),stat=astat )
        istat = istat + astat
      endif
      if( .not. allocated( wmo203 ) ) then
        allocate( wmo203(nw-1),wmo273(nw-1),stat=astat )
        istat = istat + astat
      endif
      if( .not. allocated( jpl218 ) ) then
        allocate( jpl218(nw-1),jpl295(nw-1),stat=astat )
        istat = istat + astat
      endif
      if( .not. allocated( mol226 ) ) then
        allocate( mol226(nw-1),mol263(nw-1),mol298(nw-1),stat=astat )
        istat = istat + astat
      endif
      if( .not. allocated( c0 ) ) then
        allocate( c0(nw-1),c1(nw-1),c2(nw-1),stat=astat )
        istat = istat + astat
      endif
      if( .not. allocated( no2xs_a ) ) then
        allocate( no2xs_a(nw-1),no2xs_b(nw-1),stat=astat )
        istat = istat + astat
      endif
      








      CALL o3_rei(nw,wl)
      CALL o3_jpl(nw,wl)
      CALL o3_wmo(nw,wl)
      CALL o3_mol(nw,wl)
      CALL o3_bas(nw,wl)

      CALL rdno2xs(nw,wl)

      END SUBROUTINE rdxs_init

      SUBROUTINE o3xs(nz,t,nw,wl, xs)




















      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(nw)

      INTEGER, intent(in) :: nz
      REAL, intent(in) :: t(nz)







      REAL, intent(inout) :: xs(:,:)



      INTEGER :: iz, iw
      REAL    :: factor, factor1
      REAL    :: tc(nz)









      DO iw = 1, nw-1
        IF(wl(iw) < v185) THEN
          factor = (wmo273(iw) - wmo203(iw))/(273. - 203.)
          xs(1:nz,iw) = wmo203(iw) + (t(1:nz) - 203.)*factor
          WHERE (t(1:nz) <= 203.) 
            xs(1:nz,iw) = wmo203(iw)
          ELSEWHERE (t(1:nz) >= 273.) 
            xs(1:nz,iw) = wmo273(iw)
          ENDWHERE
        ELSEIF(wl(iw) >= v185 .AND. wl(iw) < v195) THEN
          factor = (jpl295(iw) - jpl218(iw))/(295. - 218.)
          xs(1:nz,iw) = jpl218(iw) + (t(1:nz) - 218.)*factor
          WHERE (t(1:nz) <= 218.) 
            xs(1:nz,iw) = jpl218(iw)
          ELSEWHERE (t(1:nz) >= 295.) 
            xs(1:nz,iw) = jpl295(iw)
          ENDWHERE
        ELSEIF(wl(iw) >= v195 .AND. wl(iw) < v345) THEN
          factor = .1*(rei228(iw) - rei218(iw))
          WHERE( t(1:nz) < 218. ) 
            xs(1:nz,iw) = rei218(iw)
          ELSEWHERE( t(1:nz) >= 218. .AND. t(1:nz) < 228. )
            xs(1:nz,iw) = rei218(iw) + (t(1:nz) - 218.)*factor
          ENDWHERE
          factor = (rei243(iw) - rei228(iw))/15.
          WHERE( t(1:nz) >= 228. .AND. t(1:nz) < 243. )
            xs(1:nz,iw) = rei228(iw) + (t(1:nz) - 228.)*factor
          ENDWHERE
          factor = (rei295(iw) - rei243(iw))/(295. - 243.)
          WHERE( t(1:nz) >= 243. .AND. t(1:nz) < 295.)
            xs(1:nz,iw) = rei243(iw) + (t(1:nz) - 243.)*factor
          ELSEWHERE( t(1:nz) >= 295. )
            xs(1:nz,iw) = rei295(iw)
          ENDWHERE
        ELSEIF(wl(iw) >= v345) THEN
          xs(1:nz,iw) = rei295(iw)
        ENDIF
      END DO

      END SUBROUTINE o3xs



      SUBROUTINE o3_rei(nw,wl)





















      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(nw)



      INTEGER, PARAMETER :: kdata = 70000

      INTEGER n1, n2, n3, n4, iw
      REAL x1(kdata), x2(kdata), x3(kdata), x4(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata), y4(kdata)

      INTEGER i
      INTEGER ierr



      REAL ri(kdata)
      CHARACTER(len=256) :: emsg






      OPEN(UNIT=kin,FILE='DATAE1/O3/1995Malicet_O3.txt',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
        call wrf_error_fatal3("<stdin>",247,&
'o3_rei: Failed to open DATAE1/O3/1985Malicet_O3.txt' )
      endif
      DO i = 1, 2
         READ(kin,*,iostat=ierr)
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",253,&
'o3_rei: Failed to read DATAE1/O3/1985Malicet_O3.txt' )
         endif
      ENDDO
      n1 = 15001
      n2 = 15001
      n3 = 15001
      n4 = 15001
      DO i = 1, n1
         READ(kin,*,iostat=ierr) x1(i), y1(i), y2(i), y3(i), y4(i)
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",264,&
'o3_rei: Failed to read DATAE1/O3/1985Malicet_O3.txt' )
         endif
         x2(i) = x1(i)
         x3(i) = x1(i)
         x4(i) = x1(i)
      ENDDO
      CLOSE (kin)




      OPEN(UNIT=kin,FILE='DATAE1/O3/1998Brion_295.txt',STATUS='old')
      DO i = 1, 15
         READ(kin,*)
      ENDDO
      DO i = 1, 48515-15
         n1 = n1 + 1
         READ(kin,*) x1(n1), y1(n1)
      ENDDO
      CLOSE (kin)

      DO i = 1, n1
         ri(i) = refrac(x1(i), 2.45E19)
      ENDDO
      DO i = 1, n1
         x1(i) = x1(i) * ri(i)
      ENDDO

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,            1.e+38,0.)
      CALL inter2(nw,wl,rei295,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_rei: interp err = '',i5,'' in O3 xsect - Reims 295K'')') ierr
         call wrf_error_fatal3("<stdin>",300,&
trim(emsg) )
      ENDIF

      DO i = 1, n2
         ri(i) = refrac(x2(i), 2.45E19)
      ENDDO
      DO i = 1, n2
         x2(i) = x2(i) * ri(i)
         x3(i) = x2(i)
         x4(i) = x2(i)
      ENDDO

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,0.)
      CALL inter2(nw,wl,rei243,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_wmo: interp err = '',i5,'' in O3 xsect - Reims 243K'')') ierr
         call wrf_error_fatal3("<stdin>",320,&
trim(emsg) )
      ENDIF

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,               0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,            1.e+38,0.)
      CALL inter2(nw,wl,rei228,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_wmo: interp err = '',i5,'' in O3 xswct - Reims 228K'')') ierr
         call wrf_error_fatal3("<stdin>",331,&
trim(emsg) )
      ENDIF

      CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,               0.,0.)
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,            1.e+38,0.)
      CALL inter2(nw,wl,rei218,n4,x4,y4,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_wmo: interp err = '',i5,'' in O3 xswct - Reims 218K'')') ierr
         call wrf_error_fatal3("<stdin>",342,&
trim(emsg) )
      ENDIF



      v195 = 195.00 * refrac(195.00, 2.45E19)
      v345 = 345.00 * refrac(345.00, 2.45E19)
      v830 = 830.00 * refrac(830.00, 2.45E19)

      END SUBROUTINE o3_rei



      SUBROUTINE o3_wmo(nw,wl)



















      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(nw)



      INTEGER, parameter :: kdata = 200

      INTEGER n1, n2
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

      INTEGER i, idum, iw
      REAL a1, a2, dum
      INTEGER ierr



      REAL ri(kdata)
      CHARACTER(len=256) :: emsg







      OPEN(UNIT=kin,FILE='DATAE1/wmo85',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
        call wrf_error_fatal3("<stdin>",404,&
'o3_wmo: Failed to open DATAE1/wmo85' )
      endif
      DO i = 1, 3
         read(kin,*,iostat=ierr)
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",410,&
'o3_wmo: Failed to read DATAE1/wmo85' )
         endif
      ENDDO
      n1 = 158
      n2 = 158
      DO i = 1, n1
         READ(kin,*,iostat=ierr) idum, a1, a2, dum, dum, dum, y1(i), y2(i)
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",419,&
'o3_wmo: Failed to read DATAE1/wmo85' )
         endif
         x1(i) = (a1+a2)/2.
         x2(i) = (a1+a2)/2.
      ENDDO
      CLOSE (kin)



      DO i = 1, n1
         ri(i) = refrac(x1(i), 2.45E19)
      ENDDO
      DO i = 1, n1
         x1(i) = x1(i) * ri(i)
         x2(i) = x2(i) * ri(i)
      ENDDO

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
      CALL inter2(nw,wl,wmo203,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_wmo: interp err = '',i5,'' in O3 cross section - WMO - 203K'')') ierr
         call wrf_error_fatal3("<stdin>",444,&
trim(emsg) )
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
      CALL inter2(nw,wl,wmo273,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_wmo: interp err = '',i5,'' in O3 cross section - WMO - 273K'')') ierr
         call wrf_error_fatal3("<stdin>",455,&
trim(emsg) )
      ENDIF


      
      a1 = (175.438 + 176.991) / 2.
      v176 = a1 * refrac(a1,2.45E19)

      a1 = (847.5 + 852.5) / 2.
      v850 = a1 * refrac(a1, 2.45E19)

      END SUBROUTINE o3_wmo



      SUBROUTINE o3_jpl(nw,wl)


















      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(nw)



      INTEGER, parameter :: kdata = 200

      INTEGER n1, n2, iw
      REAL x1(kdata), x2(kdata)
      REAL y1(kdata), y2(kdata)

      INTEGER i
      REAL dum
      INTEGER ierr
      CHARACTER(len=256) :: emsg



      REAL ri(kdata)



      OPEN(UNIT=kin,FILE='DATAE1/O3/2006JPL_O3.txt',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
        call wrf_error_fatal3("<stdin>",514,&
'o3_jpl: Failed to open DATAE1/O3/2006JPL_O3.txt' )
      endif
      DO i = 1, 2
         read(kin,*,iostat=ierr)
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",520,&
'o3_jpl: Failed to read DATAE1/O3/2006JPL_O3.txt' )
         endif
      ENDDO
      n1 = 167
      n2 = 167
      DO i = 1, n1
         READ(kin,*,iostat=ierr) dum, dum, x1(i), y1(i), y2(i)
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",529,&
'o3_jpl: Failed to read DATAE1/O3/2006JPL_O3.txt' )
         endif
         y1(i) = y1(i) * 1.e-20
         y2(i) = y2(i) * 1.e-20
      ENDDO
      CLOSE (kin)



      DO i = 1, n1
         ri(i) = refrac(x1(i), 2.45E19)
      ENDDO
      DO i = 1, n1
         x1(i) = x1(i) * ri(i)
         x2(i) = x1(i)
      ENDDO

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
      CALL inter2(nw,wl,jpl295,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_jpl: interp err = '',i5,'' in file O3 cross section - WMO - 295K'')') ierr
         call wrf_error_fatal3("<stdin>",554,&
trim(emsg) )
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
      CALL inter2(nw,wl,jpl218,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_jpl: interp err = '',i5,'' in file O3 cross section - WMO - 218K'')') ierr
         call wrf_error_fatal3("<stdin>",565,&
trim(emsg) )
      ENDIF



      v186 = 186.051 * refrac(186.051, 2.45E19)
      v825 = 825.    * refrac(825.   , 2.45E19)


      END SUBROUTINE o3_jpl



      SUBROUTINE o3_mol(nw,wl)




















      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(nw)



      INTEGER i
      INTEGER ierr

      INTEGER, parameter :: kdata = 335
      INTEGER n1, n2, n3, iw
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)



      REAL ri(kdata)
      CHARACTER(len=256) :: emsg



      OPEN(UNIT=kin,FILE='DATAE1/O3/1986Molina.txt',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
        call wrf_error_fatal3("<stdin>",622,&
'o3_mol: Failed to open DATAE1/O3/1986Molina.txt' )
      endif
      DO i = 1, 10
        READ(kin,*,iostat=ierr)
        if( ierr /= 0 ) then
          call wrf_error_fatal3("<stdin>",628,&
'o3_mol: Failed to read DATAE1/O3/1986Molina.txt' )
        endif
      ENDDO
      n1 = 0
      n2 = 0
      n3 = 0
      DO i = 1, 121-10
         n1 = n1 + 1
         n3 = n3 + 1
         READ(kin,*,iostat=ierr) x1(n1), y1(n1),  y3(n3)
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",640,&
'o3_mol: Failed to read DATAE1/O3/1986Molina.txt' )
         endif
         x3(n3) = x1(n1)
      ENDDO
      DO i = 1, 341-122
         n1 = n1 + 1
         n2 = n2 + 1
         n3 = n3 + 1
         READ(kin,*,iostat=ierr) x1(n1), y1(n1), y2(n2), y3(n3)
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",651,&
'o3_mol: Failed to read DATAE1/O3/1986Molina.txt' )
         endif
         x2(n2) = x1(n1)
         x3(n3) = x1(n1)
      ENDDO
      CLOSE (kin)



      DO i = 1, n1
         ri(i) = refrac(x1(i), 2.45E19)
      ENDDO
      DO i = 1, n1
         x1(i) = x1(i) * ri(i)
      ENDDO

      DO i = 1, n2
         ri(i) = refrac(x2(i), 2.45E19)
      ENDDO
      DO i = 1, n2
         x2(i) = x2(i) * ri(i)
      ENDDO

      DO i = 1, n3
         ri(i) = refrac(x3(i), 2.45E19)
      ENDDO
      DO i = 1, n3
         x3(i) = x3(i) * ri(i)
      ENDDO



      v185 = 185.  * refrac(185. , 2.45E19)
      v240 = 240.5 * refrac(240.5, 2.45E19)
      v350 = 350.  * refrac(350. , 2.45E19)



      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,            1.e+38,0.)
      CALL inter2(nw,wl,mol226,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_mol: interp err = '',i5,'' in O3 xsect - 226K Molina'')') ierr
         call wrf_error_fatal3("<stdin>",697,&
trim(emsg) )
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,0.)
      CALL inter2(nw,wl,mol263,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_mol: interp err = '',i5,'' in O3 xsect - 263K Molina'')') ierr
         call wrf_error_fatal3("<stdin>",708,&
trim(emsg) )
      ENDIF

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,               0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,            1.e+38,0.)
      CALL inter2(nw,wl,mol298,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_mol: interp err = '',i5,'' in O3 xsect - 298K Molina'')') ierr
         call wrf_error_fatal3("<stdin>",719,&
trim(emsg) )
      ENDIF

      END SUBROUTINE o3_mol



      SUBROUTINE o3_bas(nw,wl)



















      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(nw)



      INTEGER, parameter :: kdata = 2000

      INTEGER i, iw
      INTEGER ierr

      INTEGER n1, n2, n3
      REAL x1(kdata), x2(kdata), x3(kdata)
      REAL y1(kdata), y2(kdata), y3(kdata)



      REAL ri(kdata)
      CHARACTER(len=256) :: emsg

      OPEN(UNIT=kin,FILE='DATAE1/O3/1985Bass_O3.txt',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
        call wrf_error_fatal3("<stdin>",768,&
'o3_bas: Failed to open DATAE1/O3/1985Bass_O3.txt' )
      endif
      DO i = 1, 8
         READ(kin,*,iostat=ierr)
      ENDDO
      n1 = 1915
      n2 = 1915
      n3 = 1915
      DO i = 1, n1
        READ(kin,*) x1(i), y1(i), y2(i), y3(i)
        if( ierr /= 0 ) then
          call wrf_error_fatal3("<stdin>",780,&
'o3_bas: Failed to read DATAE1/O3/1985Bass_O3.txt' )
        endif
      ENDDO
      CLOSE (kin)
      y1(1:n1) = 1.e-20 * y1(1:n1)
      y2(1:n1) = 1.e-20 * y2(1:n1)
      y3(1:n1) = 1.e-20 * y3(1:n1)



      DO i = 1, n1
         ri(i) = refrac(x1(i), 2.45E19)
      ENDDO
      x1(1:n1) = x1(1:n1) * ri(1:n1)
      x2(1:n1) = x1(1:n1)
      x3(1:n1) = x1(1:n1)



      vb245 = 245.018 * refrac(245.018, 2.45E19)
      vb342 = 341.981 * refrac(341.981, 2.45E19)



      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,            1.e+38,0.)
      CALL inter2(nw,wl,c0,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_bas: interp err = '',i5,'' in O3 xsect - c0 Bass'')') ierr
         call wrf_error_fatal3("<stdin>",812,&
trim(emsg) )
      ENDIF

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,0.)
      CALL inter2(nw,wl,c1,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_bas: interp err = '',i5,'' in O3 xsect - c1 Bass'')') ierr
         call wrf_error_fatal3("<stdin>",823,&
trim(emsg) )
      ENDIF

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,               0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,            1.e+38,0.)
      CALL inter2(nw,wl,c2,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''o3_bas: interp err = '',i5,'' in O3 xsect - c2 Bass'')') ierr
         call wrf_error_fatal3("<stdin>",834,&
trim(emsg) )
      ENDIF

      END SUBROUTINE o3_bas



      SUBROUTINE rdo2xs(nw,wl,o2xs1)

















      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(nw)




      REAL, intent(inout) :: o2xs1(:)



      INTEGER, parameter :: kdata = 200
      INTEGER :: i, n
      INTEGER :: ierr
      REAL    :: x, y
      REAL    :: x1(kdata), y1(kdata)
      CHARACTER(len=256) :: emsg










      n = 0

      OPEN(UNIT=kin,FILE='DATAE1/O2/O2_brasseur.abs',iostat=ierr)
      if( ierr /= 0 ) then
        call wrf_error_fatal3("<stdin>",890,&
'rdso2xs: Failed to open DATAE1/O2/O2_brasseur.abs' )
      endif
      DO i = 1, 7
         READ(kin,*,iostat=ierr)
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",896,&
'rdso2xs: Failed to read DATAE1/O2/O2_brasseur.abs' )
         endif
      ENDDO
      DO i = 1, 78
         READ(kin,*,iostat=ierr) x, y
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",903,&
'rdso2xs: Failed to read DATAE1/O2/O2_brasseur.abs' )
         endif
         IF (x .LE. 204.) THEN
            n = n + 1
            x1(n) = x
            y1(n) = y
         ENDIF
      ENDDO
      CLOSE(kin)

      OPEN(UNIT=kin,FILE='DATAE1/O2/O2_yoshino.abs',STATUS='old',iostat=ierr)
      if( ierr /= 0 ) then
         call wrf_error_fatal3("<stdin>",916,&
'rdso2xs: Failed to open DATAE1/O2/O2_yoshino.abs' )
      endif

      DO i = 1, 8
         READ(kin,*,iostat=ierr)
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",923,&
'rdso2xs: Failed to read DATAE1/O2/O2_yoshino.abs' )
         endif
      ENDDO
      DO i = 1, 36
         n = n + 1
         READ(kin,*,iostat=ierr) x, y
         if( ierr /= 0 ) then
           call wrf_error_fatal3("<stdin>",931,&
'rdso2xs: Failed to read DATAE1/O2/O2_yoshino.abs' )
         endif
         y1(n) = y*1.E-24
         x1(n) = x
      END DO
      CLOSE (kin)




      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,0.               ,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,              1.E+38,0.)
      CALL inter2(nw,wl,o2xs1, n,x1,y1, ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''rdo2xs: interp err = '',i5,'' in O2 -> O + O'')') ierr
         call wrf_error_fatal3("<stdin>",949,&
trim(emsg) )
      ENDIF

      END SUBROUTINE rdo2xs



      SUBROUTINE rdno2xs(nw,wl)
















      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(nw)


      INTEGER, parameter :: kdata = 100
      INTEGER :: iw
      INTEGER :: i, n1, n2, ierr
      REAL    :: dum1, dum2
      REAL    :: x1(kdata), x2(kdata), y1(kdata), y2(kdata)




      OPEN(UNIT=kin,FILE='DATAE1/NO2/NO2_jpl2006.abs',STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO 
      n1 = 81
      DO i = 1, n1
         READ(kin,*) dum1, dum2, y1(i), y2(i)
         x1(i) = 0.5 * (dum1 + dum2)
         x2(i) = x1(i) 
         y1(i) = y1(i)*1.E-20
         y2(i) = y2(i)*1.E-20
      ENDDO
      CLOSE(kin)
      n2 = n1

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),   0.)
      CALL addpnt(x1,y1,kdata,n1,            1.e+38,   0.)
      CALL inter2(nw,wl,no2xs_a,n1,x1,y1,ierr)
      
      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),   0.)
      CALL addpnt(x2,y2,kdata,n2,            1.e+38,   0.)
      CALL inter2(nw,wl,no2xs_b,n2,x2,y2,ierr)

      END SUBROUTINE rdno2xs



      SUBROUTINE no2xs_jpl06a(nz,t,nw,wl, no2xs)





      INTEGER, intent(in) :: nz
      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: t(nz)
      REAL, intent(in)    :: wl(nw)



      REAL, intent(inout) :: no2xs(:,:)



      INTEGER :: iw
      REAL    :: tfac(nz)
      
      tfac(1:nz) = (t(1:nz) - 220.)/74.
      DO iw = 1, nw-1
        no2xs(1:nz,iw) = no2xs_a(iw) + (no2xs_b(iw)-no2xs_a(iw))*tfac(1:nz)
      ENDDO 

      END SUBROUTINE no2xs_jpl06a



      SUBROUTINE rdso2xs(nw,wl,so2xs)






      INTEGER, parameter :: kdata = 1000


      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(nw)



      REAL, intent(inout) :: so2xs(nw)


      REAL x1(kdata)
      REAL y1(kdata)
      REAL yg(kw)
      REAL dum
      INTEGER ierr
      INTEGER i, l, n, idum
      CHARACTER(len=40)  :: fil
      CHARACTER(len=256) :: emsg





      fil = 'DATA/McGee87'
      OPEN(UNIT=kin,FILE='DATAE1/SO2/SO2xs.all',STATUS='old')
      DO i = 1,3 
        read(kin,*)
      ENDDO
      n = 704 
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
        x1(i) = .1*x1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
      CALL inter2(nw,wl,so2xs,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(emsg,'(''rdso2xs: interp err = '',i5,'' in file '',a)') ierr, fil
         call wrf_error_fatal3("<stdin>",1097,&
trim(emsg) )
      ENDIF

      END SUBROUTINE rdso2xs

      real FUNCTION refrac(w,airden)

      IMPLICIT NONE



      REAL, intent(in) :: w, airden






      REAL :: sig,  sigsq, dum





      IF (w < 200.) then
        dum = 5.e-3
      ELSEIF (w > 2000.) then
        dum = 5.e-4
      ELSE
        dum = 1./w
      ENDIF
      sig = 1.E3*dum
      sigsq = sig * sig

      dum = 8342.13 + 2406030./(130. - sigsq) + 15997./(38.9 - sigsq)


      dum = dum * airden/(2.69e19 * 273.15/288.15)


      refrac = 1. + 1.E-8 * dum

      END function refrac

      end module module_xsections
