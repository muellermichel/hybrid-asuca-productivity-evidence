module lbc
  use nrtype, only: rp
  real(rp) :: tratio_bnd
  real(rp) :: mtratio_bnd
contains
  subroutine lbc_run_rk_short_dens_rldamp
    use prm, only: nx_mn, nx_mx, ny_mn, ny_mx, nx, ny, nz_mn, nz_mx
    use mpi_control, only: mpi_calc, mpi_my_rank
    use timeset_vars, only : timestep_counter_s, dt_s, dt_rk_s
    use ref, only : dens_ref_f
    implicit none

    integer(4):: i
    integer(4):: j
    integer(4):: k

  ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    if( mpi_calc(mpi_my_rank) /= 1 ) then
      return
    end if

    tratio_bnd = ( &
      &   ( current_time + (timestep_counter_s - 1) * dt_s + dt_rk_s ) / 60._rp &
      & - ( iseq_bnd(1) - iseq_ini ) ) &
      & / ( iseq_bnd(2) - iseq_bnd(1) )
    mtratio_bnd = 1._rp - tratio_bnd

  !   ====================================================================
  !   >>>   update gpv for lateral and upper damipng                   <<<
  !   ====================================================================
  !$OMP PARALLEL DO
    do j = ny_mn, ny_mx
    do i = nx_mn, nx_mx
    do k = nz_mn, nz_mx
      dens_ptb_damp(k,i,j) = &
        &  mtratio_bnd * ( dens_ref_f(k,i,j)   &
        &                + dens_ptb_bnd(k,i,j,1)   ) &
        & + tratio_bnd * ( dens_ref_f(k,i,j)   &
        &                + dens_ptb_bnd(k,i,j,2)   ) &
        & - dens_ref_f(k,i,j)
    end do
    end do
    end do
  !$OMP END PARALLEL DO

    return
  end subroutine lbc_run_rk_short_dens_rldamp
end module