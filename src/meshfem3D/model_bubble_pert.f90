!
!       add a sharp bubble perturbation on a 1D reference model
!       the bubble is tapered by a gaussian (half width half max should be 5-10% of the bubble radius)
!



module model_bubble_pert_par

  use constants

  !! MPC gaussain and bubble perturbation definitions
  double precision :: Ampl_pert_vp_read, Ampl_pert_vs_read, Ampl_pert_rho_read
  double precision :: Ampl_pert_vp, Ampl_pert_vs, Ampl_pert_rho
  double precision :: x_center_gauss, y_center_gauss, z_center_gauss
  double precision :: bubble_radius, HWHM

end module model_bubble_pert_par


!
!--------------------------------------------------------------------------------------------------
!
  subroutine model_bubble_pert_broadcast()

    use model_bubble_pert_par

    double precision              :: lon_center_gauss_in_degree, lat_center_gauss_in_degree
    double precision              :: depth_center_gauss_in_m, radius_center_gauss
    double precision              :: theta,phi
    double precision              :: USED_RADIUS_EARTH = 6371000.d0
    character(len=10)             :: line


    !! MPC  INPUTS FOR BUBBLE PERTURBATION
    open(27,file='DATA/ParFileBubblePertRelative',action='read')
    read(27, '(a)') line
    read(27, '(a)') line
    read(27, *)     lat_center_gauss_in_degree, lon_center_gauss_in_degree, depth_center_gauss_in_m
    read(27, '(a)') line
    read(27, *)     HWHM
    read(27, '(a)') line
    read(27, *)     bubble_radius
    read(27, '(a)') line
    read(27, *)     Ampl_pert_vp_read
    read(27, '(a)') line
    read(27, *)     Ampl_pert_vs_read
    read(27, '(a)') line
    read(27, *)     Ampl_pert_rho_read
    close(27)

    Ampl_pert_vp  = dble(Ampl_pert_vp_read)  / 100.d0
    Ampl_pert_vs  = dble(Ampl_pert_vs_read)  / 100.d0
    Ampl_pert_rho = dble(Ampl_pert_rho_read) / 100.d0

    !! switch to relative values
     HWHM = HWHM /  USED_RADIUS_EARTH
    bubble_radius = bubble_radius /  USED_RADIUS_EARTH
    !! relative radius of gaussian perturbation
    radius_center_gauss = (USED_RADIUS_EARTH - depth_center_gauss_in_m) / USED_RADIUS_EARTH

    !! convert to (theta, phi) (co-latitude and longitude in radian
    phi = lon_center_gauss_in_degree * DEGREES_TO_RADIANS
    call lat_2_geocentric_colat_dble(lat_center_gauss_in_degree, theta)

    if (myrank == 0 ) then
       write(IMAIN,*)
       write(IMAIN,*) '       incorporating bubble perturbation '
       write(IMAIN,*)
       write(IMAIN,*) '             center in degree  : ',  lat_center_gauss_in_degree, lon_center_gauss_in_degree
       write(IMAIN,*) '             center in radian  : ',  theta, phi
       write(IMAIN,*) '             normalized radius : ',  radius_center_gauss
    end if

    !! switch to specfem3D_globe coordiante system ----
    call rthetaphi_2_xyz_dble(x_center_gauss, y_center_gauss, z_center_gauss, &
         radius_center_gauss,  theta, phi)

    if (myrank == 0) then
       write(IMAIN,*) ' bubble pert center in specfem3D_globe coordinate system ', &
            x_center_gauss, y_center_gauss, z_center_gauss
       write(IMAIN,*)
    end if

  end subroutine model_bubble_pert_broadcast

!
!--------------------------------------------------------------------------------------------------
!

  ! compute the bubble perturbation with a gaussian taper
  subroutine model_bubble_pert(x, y, z, drho, dvp, dvs)

    use model_bubble_pert_par

    double precision,       intent(in)    :: x, y, z
    double precision, intent(inout)       :: drho, dvp, dvs
    double precision                      :: stddev
    double precision                      :: distance
    double precision                      :: gauss_value

    stddev = HWHM / sqrt(2. * log(2.))

    distance = sqrt((x - x_center_gauss) ** 2 + &
                    (y - y_center_gauss) ** 2 + &
                    (z - z_center_gauss) ** 2) - bubble_radius

    if (distance < 0.) then
      distance = 0.
    endif

    if (distance > 4. * HWHM) then
      dvp = 0.
      dvs = 0.
      drho = 0.
      return
    endif

    gauss_value = exp( -0.5d0 * ((distance/stddev) ** 2))

    dvp = Ampl_pert_vp * gauss_value
    dvs = Ampl_pert_vs * gauss_value
    drho = Ampl_pert_rho * gauss_value


  end subroutine model_bubble_pert


!
!--------------------------------------------------------------------------------------------------
!
