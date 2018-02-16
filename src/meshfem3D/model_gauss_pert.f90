!
!                     add gaussain perturbation oin 1D reference model  
!



module model_gauss_pert_par

  use constants

  !! parameters of gaussian peeturbation 
  double precision              ::  x_center_gauss, y_center_gauss, z_center_gauss
  double precision              ::  Ampl_rho, Ampl_vp, Ampl_vs
  double precision              ::  sigma_gauss 

end module model_gauss_pert_par


!
!--------------------------------------------------------------------------------------------------
!
  subroutine model_gauss_pert_broadcast()
    
    use model_gauss_pert_par
    
    double precision              :: lon_center_gauss_in_degree, lat_center_gauss_in_degree
    double precision              :: depth_center_gauss, radius_gauss
    double precision              :: theta,phi
    double precision              :: USED_RADIUS_EARTH = 6371000.d0

    !!  INPUTS FOR GAUSSAIN PERTURBATION 
    !! define gaussain perturbation  : need to read a file but for now just 
    !! hardcoded here ---------------------------------------------------------
    
    !! center of gaussian perturbation (geocentric values) 
    lon_center_gauss_in_degree = 1.5   ! (degree)
    lat_center_gauss_in_degree = 42.5  ! (degree)
    !! depth (in m)
    depth_center_gauss = 1500000.
    !! spatial standard deviation of gaussian perturbation (m)
    sigma_gauss = 500000. 
   
    !! relative perturbation 
    Ampl_vp = -0.2
    Ampl_vs = -0.1
    Ampl_rho = -0.1
    !!! -------------------------------------------------------------------
    
    !! switch to relative values  
    sigma_gauss = sigma_gauss /  USED_RADIUS_EARTH

    !! relative radius of gaussian perturbation 
    radius_gauss = (USED_RADIUS_EARTH - depth_center_gauss) / USED_RADIUS_EARTH
   
    !! convert to (theta, phi) (co-latitude and longitude in radian 
    phi = lon_center_gauss_in_degree * DEGREES_TO_RADIANS
    call lat_2_geocentric_colat_dble(lat_center_gauss_in_degree, theta)

    if (myrank == 0 ) then 
       write(IMAIN,*)
       write(IMAIN,*) '       incorporating gaussian perturbation '
       write(IMAIN,*)
       write(IMAIN,*) '             center in degree  : ',  lat_center_gauss_in_degree, lon_center_gauss_in_degree
       write(IMAIN,*) '             center in radian  : ',  theta, phi
       write(IMAIN,*) '             normalized radius : ',  radius_gauss
    end if

    !! switch to specfem3D_globe coordiante system ----
    call rthetaphi_2_xyz_dble(x_center_gauss, y_center_gauss, z_center_gauss, &
         radius_gauss,  theta, phi)

    if (myrank == 0) then 
       write(IMAIN,*) ' gaussian pert center in specfem3D_globe coordinate system ', &
            x_center_gauss, y_center_gauss, z_center_gauss
       write(IMAIN,*)
    end if

  end subroutine model_gauss_pert_broadcast

!
!--------------------------------------------------------------------------------------------------
!

  ! compute gaussian perturbation 
  subroutine model_gaussian_pert(x, y, z, drho, dvp, dvs)
    
    use model_gauss_pert_par

    double precision, intent(in)    :: x, y, z
    double precision, intent(inout) :: drho, dvp, dvs
    double precision                :: gauss_value
    
    gauss_value = ( ( x - x_center_gauss ) / sigma_gauss )**2 + &
                  ( ( y - y_center_gauss ) / sigma_gauss )**2 + &
                  ( ( z - z_center_gauss ) / sigma_gauss )**2

    gauss_value = exp ( -0.5d0 * gauss_value )

    
    drho =  gauss_value * Ampl_rho
    dvp =   gauss_value * Ampl_vp
    dvs =   gauss_value * Ampl_vs

  end subroutine model_gaussian_pert
    
  
!
!--------------------------------------------------------------------------------------------------
!

