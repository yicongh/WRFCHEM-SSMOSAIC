MODULE module_vertmx_wrf

CONTAINS


    SUBROUTINE vertmx( dt, phi, kt_turb, dryrho,  &
       zsigma, zsigma_half, vd, kts, ktem1 )























      IMPLICIT NONE 


      INTEGER, INTENT(IN) :: kts,ktem1
      REAL, INTENT(IN) :: dt, vd


      REAL, INTENT(IN), DIMENSION (kts:ktem1+1) :: kt_turb, zsigma
      REAL, INTENT(IN), DIMENSION (kts:ktem1) :: dryrho, zsigma_half
      REAL, INTENT(INOUT), DIMENSION (kts:ktem1) :: phi


      INTEGER :: k


      REAL, DIMENSION (kts+1:ktem1) :: a_coeff
      REAL, DIMENSION (kts:ktem1) :: b_coeff, lhs1, lhs2, lhs3, rhs




      CALL coeffs( kts, ktem1, dryrho, zsigma, zsigma_half, a_coeff, b_coeff )

      CALL rlhside( kts, ktem1, kt_turb, dryrho, a_coeff, b_coeff, &
        phi, dt, vd, rhs, lhs1, lhs2, lhs3 )

      CALL tridiag( kts, ktem1, lhs1, lhs2, lhs3, rhs )

      DO k = kts,ktem1
        phi(k) = rhs(k)
      END DO

    END SUBROUTINE vertmx



    SUBROUTINE rlhside( kts, ktem1, k_turb, dryrho, a_coeff, b_coeff, &
      phi, dt, vd, rhs, lhs1, lhs2, lhs3 )
  
  
  
  
  
  
      IMPLICIT NONE 


      INTEGER, INTENT(IN) :: kts,ktem1
      REAL, INTENT(IN) :: dt, vd


      REAL, INTENT(IN), DIMENSION (kts:ktem1+1) :: k_turb
      REAL, INTENT(IN), DIMENSION (kts+1:ktem1) :: a_coeff
      REAL, INTENT(IN), DIMENSION (kts:ktem1) :: b_coeff, dryrho, phi      
      REAL, INTENT(OUT), DIMENSION (kts:ktem1) :: lhs1, lhs2, lhs3, rhs


      REAL :: a1, a2, alfa_explicit = .25, beta_implicit = .75
      INTEGER :: i


      i = kts
      a2 = a_coeff(i+1)*k_turb(i+1)
      rhs(i) = (1./(dt*b_coeff(i)) - alfa_explicit*(vd*dryrho(i)+a2))*phi(i) + &
        alfa_explicit*(a2*phi(i+1))
      lhs1(i) = 0.
      lhs2(i) = 1./(dt*b_coeff(i)) + beta_implicit*(vd*dryrho(i)+a2)
      lhs3(i) = -beta_implicit*a2

      DO i = kts+1, ktem1-1
        a1 = a_coeff(i)*k_turb(i)
        a2 = a_coeff(i+1)*k_turb(i+1)

        rhs(i) = (1./(dt*b_coeff(i)) - alfa_explicit*(a1+a2))*phi(i) + &
          alfa_explicit*(a1*phi(i-1) + a2*phi(i+1))

        lhs1(i) = -beta_implicit*a1
        lhs2(i) = 1./(dt*b_coeff(i)) + beta_implicit*(a1+a2)
        lhs3(i) = -beta_implicit*a2
      END DO

      i = ktem1
      a1 = a_coeff(i)*k_turb(i)
      rhs(i) = (1./(dt*b_coeff(i)) - alfa_explicit*(a1   ))*phi(i) + &
        alfa_explicit*(a1*phi(i-1))
      lhs1(i) = -beta_implicit*a1
      lhs2(i) = 1./(dt*b_coeff(i)) + beta_implicit*(a1   )
      lhs3(i) = 0.

    END SUBROUTINE rlhside



    SUBROUTINE tridiag( kts, ktem1, a, b, c, f )
  
  
  
  
  
  
  
  
  
      IMPLICIT NONE 


      INTEGER, INTENT(IN) :: kts,ktem1


      REAL, INTENT(IN), DIMENSION (kts:ktem1) :: a, b, c
      REAL, INTENT(INOUT), DIMENSION (kts:ktem1) :: f


      REAL :: p
      INTEGER :: i


      REAL, DIMENSION (kts:ktem1) :: q

      q(kts) = -c(kts)/b(kts)
      f(kts) = f(kts)/b(kts)

      DO i = kts+1, ktem1
        p = 1./(b(i)+a(i)*q(i-1))
        q(i) = -c(i)*p
        f(i) = (f(i)-a(i)*f(i-1))*p
      END DO

      DO i = ktem1 - 1, kts, -1
        f(i) = f(i) + q(i)*f(i+1)
      END DO

    END SUBROUTINE tridiag



    SUBROUTINE coeffs( kts, ktem1, dryrho, &
      z_sigma, z_sigma_half, a_coeff, b_coeff )






      IMPLICIT NONE 

      INTEGER, INTENT(IN) :: kts,ktem1


      REAL, INTENT(IN), DIMENSION (kts:ktem1+1) :: z_sigma
      REAL, INTENT(IN), DIMENSION (kts:ktem1) :: z_sigma_half, dryrho
      REAL, INTENT(OUT), DIMENSION (kts+1:ktem1) :: a_coeff
      REAL, INTENT(OUT), DIMENSION (kts:ktem1) :: b_coeff


      INTEGER :: i
      REAL :: dryrho_at_w

      DO i = kts, ktem1
        b_coeff(i) = 1./(dryrho(i)*(z_sigma(i+1)-z_sigma(i)))
      END DO

      DO i = kts+1, ktem1
        dryrho_at_w = 0.5*(dryrho(i)+dryrho(i-1))
        a_coeff(i) = dryrho_at_w/(z_sigma_half(i)-z_sigma_half(i-1))
      END DO

    END SUBROUTINE coeffs


END MODULE module_vertmx_wrf
