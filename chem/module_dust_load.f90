MODULE module_dust_load


 CONTAINS

   SUBROUTINE dust_load_driver ( config_flags,                         &
            alt, chem, dz8w, dustload_1, dustload_2, dustload_3,       &
            dustload_4, dustload_5,                                    &
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )

   USE module_configure
   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::                                   &
                                      ids,ide, jds,jde, kds,kde,       &
                                      ims,ime, jms,jme, kms,kme,       &
                                      its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(IN ) :: chem

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(IN ) :: alt, dz8w

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: dustload_1,   &
                dustload_2, dustload_3, dustload_4, dustload_5

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER :: i, j, k



       dustload_1(its:ite,jts:jte) = 0.
       dustload_2(its:ite,jts:jte) = 0.
       dustload_3(its:ite,jts:jte) = 0.
       dustload_4(its:ite,jts:jte) = 0.
       dustload_5(its:ite,jts:jte) = 0.


     do j=jts,jte
      do i=its,ite
       do k=kts,kte

        dustload_1(i,j)= dustload_1(i,j) + chem(i,k,j,p_dust_1)/alt(i,k,j) * dz8w(i,k,j)
        dustload_2(i,j)= dustload_2(i,j) + chem(i,k,j,p_dust_2)/alt(i,k,j) * dz8w(i,k,j)
        dustload_3(i,j)= dustload_3(i,j) + chem(i,k,j,p_dust_3)/alt(i,k,j) * dz8w(i,k,j)
        dustload_4(i,j)= dustload_4(i,j) + chem(i,k,j,p_dust_4)/alt(i,k,j) * dz8w(i,k,j)
        dustload_5(i,j)= dustload_5(i,j) + chem(i,k,j,p_dust_5)/alt(i,k,j) * dz8w(i,k,j)

       enddo
      enddo
     enddo

   END SUBROUTINE dust_load_driver

END MODULE module_dust_load

