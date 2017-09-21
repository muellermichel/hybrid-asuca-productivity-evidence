module sf_slab_flx
  contains
  subroutine sf_slab_flx_tile_run( &
    & month, &
    & u1, v1, pt1, qv1, r_exner_surf, z_f1, &
    & skind, tg, beta, qsatg, solar, tlcvr, &
    & z_0m, z_0h, z_0q, r_mol, u_f, &
    & taux_surf_ex, tauy_surf_ex, ftl_surf_ex, fqw_surf_ex, &
    & r_mol_surf, u_f_surf, &
    & taux_tile_ex, tauy_tile_ex, ftl_tile_ex, fqw_tile_ex, &
    & dfm_tile, dfh_tile, dfq_tile)

    use pp_vardef
    use sf_slab_grid, only: ngm
    use sf_slab_skind, only: kind_sea
    use sf_slab_stile, only: ntlm, tile_land, tile_sea
    implicit none

    integer(4),   intent(in) :: month
    real(r_size), intent(in) :: u1
    real(r_size), intent(in) :: v1
    real(r_size), intent(in) :: pt1
    real(r_size), intent(in) :: qv1
    real(r_size), intent(in) :: r_exner_surf
    real(r_size), intent(in) :: z_f1

    integer(4),   intent(in) :: skind(ntlm)
    real(r_size), intent(in) :: tg(0:ngm, ntlm)
    real(r_size), intent(in) :: beta(ntlm)
    real(r_size), intent(in) :: qsatg(ntlm)
    real(r_size), intent(in) :: solar(ntlm)
    real(r_size), intent(in) :: tlcvr(ntlm)

    real(r_size), intent(inout) :: z_0m(ntlm)
    real(r_size), intent(inout) :: z_0h(ntlm)
    real(r_size), intent(inout) :: z_0q(ntlm)
    real(r_size), intent(inout) :: r_mol(ntlm)
    real(r_size), intent(inout) :: u_f(ntlm)

    real(r_size), intent(out) :: taux_surf_ex
    real(r_size), intent(out) :: tauy_surf_ex
    real(r_size), intent(out) :: ftl_surf_ex
    real(r_size), intent(out) :: fqw_surf_ex
    real(r_size), intent(out) :: r_mol_surf
    real(r_size), intent(out) :: u_f_surf

    real(r_size), intent(out) :: taux_tile_ex(ntlm)
    real(r_size), intent(out) :: tauy_tile_ex(ntlm)
    real(r_size), intent(out) :: ftl_tile_ex(ntlm)
    real(r_size), intent(out) :: fqw_tile_ex(ntlm)
    real(r_size), intent(out) :: dfm_tile(ntlm)
    real(r_size), intent(out) :: dfh_tile(ntlm)
    real(r_size), intent(out) :: dfq_tile(ntlm)

    integer(4) :: lt


    lt = tile_land
    if (tlcvr(lt) > 0.0_r_size) then
      call sf_slab_flx_land_run(&
        & month, &
        & u1, v1, pt1, qv1, &
        & tg(0, lt), beta(lt), qsatg(lt), r_exner_surf, solar(lt), &
        & z_f1, z_0m(lt), z_0h(lt), z_0q(lt), r_mol(lt), &
        & dfm_tile(lt), dfh_tile(lt), dfq_tile(lt), &
        & taux_tile_ex(lt), tauy_tile_ex(lt), &
        & ftl_tile_ex(lt), fqw_tile_ex(lt))

      u_f(lt) = sqrt(sqrt(taux_tile_ex(lt) ** 2 + tauy_tile_ex(lt) ** 2))
    else
      dfm_tile(lt) = 0.0_r_size
      dfh_tile(lt) = 0.0_r_size
      dfq_tile(lt) = 0.0_r_size
      taux_tile_ex(lt) = 0.0_r_size
      tauy_tile_ex(lt) = 0.0_r_size
      ftl_tile_ex(lt)  = 0.0_r_size
      fqw_tile_ex(lt)  = 0.0_r_size
    end if

    lt = tile_sea
    if (tlcvr(lt) > 0.0_r_size) then
      if (skind(lt) == kind_sea) then
        ! sea
        call sf_slab_flx_sea_run(&
          & u1, v1, pt1, qv1, &
          & tg(0, lt), beta(lt), qsatg(lt), r_exner_surf, u_f(lt), &
          & z_f1, r_mol(lt), z_0m(lt), z_0h(lt), z_0q(lt), &
          & dfm_tile(lt), dfh_tile(lt), dfq_tile(lt), &
          & taux_tile_ex(lt), tauy_tile_ex(lt), &
          & ftl_tile_ex(lt), fqw_tile_ex(lt))
      else
        ! seaice
        call sf_slab_flx_land_run(&
          & month, &
          & u1, v1, pt1, qv1, &
          & tg(0, lt), beta(lt), qsatg(lt), r_exner_surf, solar(lt), &
          & z_f1, z_0m(lt), z_0h(lt), z_0q(lt), r_mol(lt), &
          & dfm_tile(lt), dfh_tile(lt), dfq_tile(lt), &
          & taux_tile_ex(lt), tauy_tile_ex(lt), &
          & ftl_tile_ex(lt), fqw_tile_ex(lt))
      end if

      u_f(lt) = sqrt(sqrt(taux_tile_ex(lt) ** 2 + tauy_tile_ex(lt) ** 2))
    else
      dfm_tile(lt) = 0.0_r_size
      dfh_tile(lt) = 0.0_r_size
      dfq_tile(lt) = 0.0_r_size
      taux_tile_ex(lt) = 0.0_r_size
      tauy_tile_ex(lt) = 0.0_r_size
      ftl_tile_ex(lt)  = 0.0_r_size
      fqw_tile_ex(lt)  = 0.0_r_size
    end if

    taux_surf_ex = 0.0_r_size
    tauy_surf_ex = 0.0_r_size
    ftl_surf_ex  = 0.0_r_size
    fqw_surf_ex  = 0.0_r_size
    r_mol_surf = 0.0_r_size
    u_f_surf = 0.0_r_size
    do lt = 1, ntlm
      taux_surf_ex = taux_surf_ex + taux_tile_ex(lt) * tlcvr(lt)
      tauy_surf_ex = tauy_surf_ex + tauy_tile_ex(lt) * tlcvr(lt)
      ftl_surf_ex  = ftl_surf_ex  + ftl_tile_ex(lt)  * tlcvr(lt)
      fqw_surf_ex  = fqw_surf_ex  + fqw_tile_ex(lt)  * tlcvr(lt)

      r_mol_surf = r_mol_surf + r_mol(lt) * tlcvr(lt)
      u_f_surf = u_f_surf + u_f(lt) * tlcvr(lt)
    end do

    return
  end subroutine sf_slab_flx_tile_run
end module