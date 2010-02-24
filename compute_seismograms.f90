!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine compute_seismograms(nrec_local,nrec,displ_crust_mantle, &
                                nu,hxir_store,hetar_store,hgammar_store, &
                                scale_displ,ibool_crust_mantle, &
                                ispec_selected_rec,number_receiver_global, &
                                seismo_current,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
                                seismograms)

  implicit none
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer nrec_local,nrec
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
    displ_crust_mantle

  double precision, dimension(NDIM,NDIM,nrec) :: nu

  double precision, dimension(nrec_local,NGLLX) :: hxir_store
  double precision, dimension(nrec_local,NGLLY) :: hetar_store
  double precision, dimension(nrec_local,NGLLZ) :: hgammar_store
  
  double precision scale_displ

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  integer, dimension(nrec) :: ispec_selected_rec
  integer, dimension(nrec_local) :: number_receiver_global
  
  integer :: seismo_current  
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS
  
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: &
    seismograms
  
  ! local parameters
  double precision :: uxd,uyd,uzd,hlagrange
  integer :: i,j,k,iglob,irec_local,irec

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    ! perform the general interpolation using Lagrange polynomials
    uxd = ZERO
    uyd = ZERO
    uzd = ZERO

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          iglob = ibool_crust_mantle(i,j,k,ispec_selected_rec(irec))

          hlagrange = hxir_store(irec_local,i)*hetar_store(irec_local,j)*hgammar_store(irec_local,k)

          uxd = uxd + dble(displ_crust_mantle(1,iglob))*hlagrange
          uyd = uyd + dble(displ_crust_mantle(2,iglob))*hlagrange
          uzd = uzd + dble(displ_crust_mantle(3,iglob))*hlagrange

        enddo
      enddo
    enddo
    ! store North, East and Vertical components

    ! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      seismograms(:,irec_local,seismo_current) = sngl(scale_displ*(nu(:,1,irec)*uxd + &
                 nu(:,2,irec)*uyd + nu(:,3,irec)*uzd))
    else
      seismograms(:,irec_local,seismo_current) = scale_displ*(nu(:,1,irec)*uxd + &
                 nu(:,2,irec)*uyd + nu(:,3,irec)*uzd)
    endif  

  enddo
  
  end subroutine compute_seismograms

!
!-------------------------------------------------------------------------------------------------
!
  
  subroutine compute_seismograms_backward(nrec_local,nrec,b_displ_crust_mantle, &
                                nu,hxir_store,hetar_store,hgammar_store, &
                                scale_displ,ibool_crust_mantle, &
                                ispec_selected_rec,number_receiver_global, &
                                seismo_current,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
                                seismograms)

  implicit none
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer nrec_local,nrec
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_displ_crust_mantle

  double precision, dimension(NDIM,NDIM,nrec) :: nu

  double precision, dimension(nrec_local,NGLLX) :: hxir_store
  double precision, dimension(nrec_local,NGLLY) :: hetar_store
  double precision, dimension(nrec_local,NGLLZ) :: hgammar_store
  
  double precision scale_displ

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  integer, dimension(nrec) :: ispec_selected_rec
  integer, dimension(nrec_local) :: number_receiver_global
  
  integer :: seismo_current
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS
  
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: &
    seismograms
  
  ! local parameters
  double precision :: uxd,uyd,uzd,hlagrange
  integer :: i,j,k,iglob,irec_local,irec

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    ! perform the general interpolation using Lagrange polynomials
    uxd = ZERO
    uyd = ZERO
    uzd = ZERO

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          iglob = ibool_crust_mantle(i,j,k,ispec_selected_rec(irec))

          hlagrange = hxir_store(irec_local,i)*hetar_store(irec_local,j)*hgammar_store(irec_local,k)

          uxd = uxd + dble(b_displ_crust_mantle(1,iglob))*hlagrange
          uyd = uyd + dble(b_displ_crust_mantle(2,iglob))*hlagrange
          uzd = uzd + dble(b_displ_crust_mantle(3,iglob))*hlagrange

        enddo
      enddo
    enddo
    ! store North, East and Vertical components

    ! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      seismograms(:,irec_local,seismo_current) = sngl(scale_displ*(nu(:,1,irec)*uxd + &
           nu(:,2,irec)*uyd + nu(:,3,irec)*uzd))
    else
      seismograms(:,irec_local,seismo_current) = scale_displ*(nu(:,1,irec)*uxd + &
           nu(:,2,irec)*uyd + nu(:,3,irec)*uzd)
    endif


  enddo
  
  end subroutine compute_seismograms_backward

!
!-------------------------------------------------------------------------------------------------
!
    
  subroutine compute_seismograms_adjoint(NSOURCES,nrec_local,displ_crust_mantle, &
                    eps_trace_over_3_crust_mantle,epsilondev_crust_mantle, &
                    nu_source,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                    hxir_store,hetar_store,hgammar_store, &
                    hpxir_store,hpetar_store,hpgammar_store, &
                    t_cmt,hdur_gaussian,DT,t0,scale_displ, &
                    hprime_xx,hprime_yy,hprime_zz, &
                    xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                    etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                    gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                    moment_der,sloc_der, &
                    NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismograms,deltat, &
                    ibool_crust_mantle,ispec_selected_source,number_receiver_global, &
                    NSTEP,it,nit_written)

  implicit none
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer NSOURCES,nrec_local
  
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
    displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY) :: &
    eps_trace_over_3_crust_mantle
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT) :: &
    epsilondev_crust_mantle

  double precision, dimension(NDIM,NDIM,NSOURCES) :: nu_source
  double precision, dimension(NSOURCES) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  double precision, dimension(nrec_local,NGLLX) :: hxir_store,hpxir_store
  double precision, dimension(nrec_local,NGLLY) :: hetar_store,hpetar_store
  double precision, dimension(nrec_local,NGLLZ) :: hgammar_store,hpgammar_store

  double precision, dimension(NSOURCES) :: t_cmt,hdur_gaussian
  double precision :: DT,t0
  double precision :: scale_displ

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
        etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
        gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nrec_local) :: moment_der
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local) :: sloc_der
  
  integer NTSTEP_BETWEEN_OUTPUT_SEISMOS
  
  real(kind=CUSTOM_REAL), dimension(NDIM*NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: &
    seismograms
  real(kind=CUSTOM_REAL) :: deltat
    
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  integer,dimension(NSOURCES) :: ispec_selected_source
  integer, dimension(nrec_local) :: number_receiver_global  
  integer :: NSTEP,it,nit_written  
  
  ! local parameters
  double precision :: uxd,uyd,uzd,hlagrange
  double precision :: eps_trace,dxx,dyy,dxy,dxz,dyz
  double precision :: eps_loc(NDIM,NDIM), eps_loc_new(NDIM,NDIM)
  double precision :: stf
  real(kind=CUSTOM_REAL) :: displ_s(NDIM,NGLLX,NGLLY,NGLLZ)
  real(kind=CUSTOM_REAL) :: eps_s(NDIM,NDIM), eps_m_s(NDIM), stf_deltat
  integer :: i,j,k,iglob,irec_local,irec,ispec

  double precision, external :: comp_source_time_function

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    ! perform the general interpolation using Lagrange polynomials
    uxd = ZERO
    uyd = ZERO
    uzd = ZERO


    eps_trace = ZERO
    dxx = ZERO
    dyy = ZERO
    dxy = ZERO
    dxz = ZERO
    dyz = ZERO

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          iglob = ibool_crust_mantle(i,j,k,ispec_selected_source(irec))

          hlagrange = hxir_store(irec_local,i)*hetar_store(irec_local,j)*hgammar_store(irec_local,k)

          uxd = uxd + dble(displ_crust_mantle(1,iglob))*hlagrange
          uyd = uyd + dble(displ_crust_mantle(2,iglob))*hlagrange
          uzd = uzd + dble(displ_crust_mantle(3,iglob))*hlagrange

          eps_trace = eps_trace + dble(eps_trace_over_3_crust_mantle(i,j,k,ispec_selected_source(irec)))*hlagrange
          dxx = dxx + dble(epsilondev_crust_mantle(1,i,j,k,ispec_selected_source(irec)))*hlagrange
          dyy = dyy + dble(epsilondev_crust_mantle(2,i,j,k,ispec_selected_source(irec)))*hlagrange
          dxy = dxy + dble(epsilondev_crust_mantle(3,i,j,k,ispec_selected_source(irec)))*hlagrange
          dxz = dxz + dble(epsilondev_crust_mantle(4,i,j,k,ispec_selected_source(irec)))*hlagrange
          dyz = dyz + dble(epsilondev_crust_mantle(5,i,j,k,ispec_selected_source(irec)))*hlagrange

          displ_s(:,i,j,k) = displ_crust_mantle(:,iglob)

        enddo
      enddo
    enddo

    eps_loc(1,1) = eps_trace + dxx
    eps_loc(2,2) = eps_trace + dyy
    eps_loc(3,3) = eps_trace - dxx - dyy
    eps_loc(1,2) = dxy
    eps_loc(1,3) = dxz
    eps_loc(2,3) = dyz
    eps_loc(2,1) = dxy
    eps_loc(3,1) = dxz
    eps_loc(3,2) = dyz

    eps_loc_new(:,:) = eps_loc(:,:)
    ! rotate to the local cartesian coordinates (e-n-z):  eps_new=P*eps*P'
    eps_loc_new(:,:) = matmul(matmul(nu_source(:,:,irec),eps_loc(:,:)), transpose(nu_source(:,:,irec)))

    ! distinguish between single and double precision for reals
    if (CUSTOM_REAL == SIZE_REAL) then
      seismograms(1,irec_local,it-nit_written) = sngl(eps_loc_new(1,1))
      seismograms(2,irec_local,it-nit_written) = sngl(eps_loc_new(2,2))
      seismograms(3,irec_local,it-nit_written) = sngl(eps_loc_new(3,3))
      seismograms(4,irec_local,it-nit_written) = sngl(eps_loc_new(1,2))
      seismograms(5,irec_local,it-nit_written) = sngl(eps_loc_new(1,3))
      seismograms(6,irec_local,it-nit_written) = sngl(eps_loc_new(2,3))
      seismograms(7:9,irec_local,it-nit_written) = sngl(scale_displ*(nu_source(:,1,irec)*uxd + &
                  nu_source(:,2,irec)*uyd + nu_source(:,3,irec)*uzd))
    else
      seismograms(1,irec_local,it-nit_written) = eps_loc_new(1,1)
      seismograms(2,irec_local,it-nit_written) = eps_loc_new(2,2)
      seismograms(3,irec_local,it-nit_written) = eps_loc_new(3,3)
      seismograms(4,irec_local,it-nit_written) = eps_loc_new(1,2)
      seismograms(5,irec_local,it-nit_written) = eps_loc_new(1,3)
      seismograms(6,irec_local,it-nit_written) = eps_loc_new(2,3)
      seismograms(7:9,irec_local,it-nit_written) = scale_displ*(nu_source(:,1,irec)*uxd + &
                  nu_source(:,2,irec)*uyd + nu_source(:,3,irec)*uzd)
    endif

    ! frechet derviatives of the source
    ispec = ispec_selected_source(irec)

    call compute_adj_source_frechet(displ_s,Mxx(irec),Myy(irec),Mzz(irec), &
                Mxy(irec),Mxz(irec),Myz(irec),eps_s,eps_m_s, &
                hxir_store(irec_local,:),hetar_store(irec_local,:),hgammar_store(irec_local,:), &
                hpxir_store(irec_local,:),hpetar_store(irec_local,:),hpgammar_store(irec_local,:), &
                hprime_xx,hprime_yy,hprime_zz, &
                xix_crust_mantle(:,:,:,ispec),xiy_crust_mantle(:,:,:,ispec),xiz_crust_mantle(:,:,:,ispec), &
                etax_crust_mantle(:,:,:,ispec),etay_crust_mantle(:,:,:,ispec),etaz_crust_mantle(:,:,:,ispec), &
                gammax_crust_mantle(:,:,:,ispec),gammay_crust_mantle(:,:,:,ispec),gammaz_crust_mantle(:,:,:,ispec))

    stf = comp_source_time_function(dble(NSTEP-it)*DT-t0-t_cmt(irec),hdur_gaussian(irec))
    stf_deltat = stf * deltat
    
    moment_der(:,:,irec_local) = moment_der(:,:,irec_local) + eps_s(:,:) * stf_deltat
    sloc_der(:,irec_local) = sloc_der(:,irec_local) + eps_m_s(:) * stf_deltat

  enddo
  
  end subroutine compute_seismograms_adjoint  