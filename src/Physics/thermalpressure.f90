!Pseudo code for TP

#ifdef BG 
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif


MODULE Thermal_pressure_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  INTERFACE Thermal_pressure_3D
     MODULE PROCEDURE Thermal_pressure_3D
  END INTERFACE
  INTERFACE heat_source
     MODULE PROCEDURE heat_source
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC  :: Thermal_pressure_3D, heat_source

  !---------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE Thermal_pressure_3D(DISC,dt,nz,alpha_th, alpha_hy,rho_c, Lambda, theta, sigma, Sh, SR, Dwn, DFinv,temperature, pressure)
    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tDiscretization)         :: DISC
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: i
    INTEGER     :: nz
    REAL        :: omega(nz)                                                  ! shear heating source
    REAL        :: Sh, SR, dt, tauV                                           ! shear stress, slip rate
    REAL        :: Dwn(nz), DFinv(nz)                                         !
    REAL        :: alpha_th, alpha_hy, rho_c                                  ! thermal and hydraulic diffusivities
    REAL        :: theta(nz), sigma(nz)                                       ! stored diffusion from previous timestep
    REAL        :: theta_current(nz), sigma_current(nz)                       ! diffusion next timestep
    REAL        :: tmp(nz)
    REAL        :: Lambda, Lambda_prime
    REAL        :: T,p, hwid
    REAL        :: temperature, pressure                                      ! temperatur, pressure in space domain
    !-------------------------------------------------------------------------!
    INTENT(IN)  :: DISC, nz, alpha_th, alpha_hy, rho_c, Lambda, Sh, SR, Dwn, DFinv
    INTENT(INOUT):: theta, sigma, temperature, pressure
    !-------------------------------------------------------------------------!


    tauV = 0D0 ! Sh*SR !fault strenght*slip rate
    Lambda_prime = Lambda*alpha_th/(alpha_hy-alpha_th)
    hwid = DISC%DynRup%hwid
    tmp = (Dwn/hwid)**2

    !1. Calculate diffusion of the field at previous timestep
    logInfo0(*) 'theta before', theta(1)
    logInfo0(*) 'sigma before', sigma(1)

    !temperature
    theta_current = theta*exp(-alpha_th*dt*tmp)
    !pore pressure + lambda'*temp
    sigma_current = sigma*exp(-alpha_hy*dt*tmp)

    !2. Add current contribution and get new temperature
    CALL heat_source(hwid,alpha_th,dt,Dwn,nz,omega)
    theta = theta_current + (tauV/rho_c)*omega
    CALL heat_source(DISC%DynRup%hwid,alpha_hy,dt,Dwn,nz,omega)
    sigma = sigma_current + ((Lambda+Lambda_prime)*tauV)/(rho_c)*omega

    !3. Recover temperature and pressure using inverse Fourier
    ! transformation with the calculated fourier coefficients

    T = 0.0
    p = 0.0

    !new contribution
    DO i=1,nz
       T = T + (DFinv(i)/hwid)*theta(i)
       p = p + (DFinv(i)/hwid)*sigma(i)
    ENDDO

    !Update pore pressure change (sigma = pore pressure + lambda'*temp)
    !In the BIEM code (Lapusta) they use T without initial value
    p = p - Lambda_prime*T

    !Temp and pore pressure change at single GP on the fault + initial values
    temperature = T + DISC%DynRup%IniTemp
    pressure = p + DISC%DynRup%IniPP_xx
    logInfo0(*) 'Theta', theta(1)
    logInfo0(*) 'Sigma', sigma(1)
    logInfo0(*) 'Theta_current', theta_current(1)
    logInfo0(*) 'Sigma_current', sigma_current(1)
    logInfo0(*) 'pore pressure in TP', pressure


  END SUBROUTINE Thermal_pressure_3D


   SUBROUTINE  heat_source(hwid, alpha, dt, Dwn, nz, omega)

    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    REAL        :: hwid
    INTEGER     :: nz                                                         ! number of points at the fault
    REAL        :: Dwn(nz)                                                    ! insert DISC%DynRup%TP_grid
    REAL        :: Sh, SR                                                     ! current shear stress and slip rate
    REAL        :: alpha, dt                                                      ! difussion parameter
    REAL        :: omega(nz)                                                  !
    REAL        :: tmp(nz)
    REAL, PARAMETER :: pi=3.141592653589793     ! CONSTANT pi
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: hwid, alpha, dt, Dwn, nz
    INTENT(OUT)   :: omega
    !-------------------------------------------------------------------------!
    !Gaussian shear zone in spectral domain, normalized by w
    tmp = (Dwn/hwid)**2
    !original function in spatial domain
    !omega = 1/(w*sqrt(2*pi))*exp(-0.5*(z/hwid).^2);
    !function in the wavenumber domain
    !omega = 1/(sqrt(2.0*pi))*exp(-0.5*(Dwn*hwid)**2)*(1-exp(-alpha_th*tmp))
    !inserting Dwn/hwid (scaled) for Dwn cancels out hwid
    omega = 1.0/(alpha*tmp*(sqrt(2.0*pi)))*exp(-0.5*(Dwn)**2)*(1.0 - exp(-alpha*dt*tmp))

    RETURN

   END SUBROUTINE

END MODULE Thermal_pressure_mod
