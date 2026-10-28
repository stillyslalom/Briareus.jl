"""
    Briareus

View-factor radiation and laser-pointing design for high-energy-density experiments:
indirect-drive hohlraums, backlighter foils, and other laser-irradiated targets.

Units: lengths mm, time ns, power TW, flux TW/mm² (`to_Wcm2` converts), radiation
temperature eV, angles in degrees at the user interface.
"""
module Briareus

using LinearAlgebra
using Printf
using Random
using StaticArrays
using ForwardDiff
using SpecialFunctions: gamma, gamma_inc, gamma_inc_inv, dawson, erf
import Optim
using PrecompileTools: @setup_workload, @compile_workload

export CurvedSurface, SphereSurface, CylinderSurface, ConeSurface, TorusSurface, DiskSurface, surface_point, element_sector
export Vec3, Mesh, npatches, add_polygon!, merge_meshes, split_mesh, transform, translate, rotation_to,
       orient_mesh, coarsen_mesh, remesh_mesh, SurfaceRegion, fit_surfaces, curved_components,
       cylinder_mesh, annulus_mesh, disk_mesh, sphere_mesh, cubed_sphere_mesh, icosphere_mesh,
       icosphere_frequency, rectangle_mesh, frustum_mesh
export AlbedoModel, ConstantAlbedo, TabulatedAlbedo, WallLossAlbedo, BurnthroughAlbedo, HeatFront,
       RosselandTransmission, foil_transmission, GapClosure, LipClosure, lip_depth,
       heat_front, burnthrough_time, areal_mass, optical_depth, eddington_transmission, albedo, Material, gold, ablator, BLACK,
       TemperatureAlbedo, CapsuleAlbedo, RadiusSchedule, TabulatedRadius,
       ScaledRadius, TrajectoryRadius, radius_fraction, du_wall_loss, wall_loss_factor, uranium,
       nif_gold, nif_uranium,
       Conversion, conversion_efficiency, nif_conversion,
       Component, aperture, Target, nelements, elements, component_index, parts, part_pieces,
       element_quadrature
export hohlraum, hohlraum_counts, capsule_counts, foil
export WallWindow, Sombrero, conforming_grid, equal_area_radius
export SuperGaussianSpot, Beam, port_vector, port_angles, direction, spot_axes, remake, intensity, wavelength,
       beams_from_table, nif_quads, nif_beams, nif_spot, NIF_QUAD_NAMES, nif_port, nif_bundle, nif_bundles, hybride_mock_pointing, omega_ports, omega_beams, omega_spot
export DepositionOperator, Deposition, deposit, CoronaProfile, CoronaOperator, raytrace_deposition, raytrace_component_power, laser_flux,
       component_power
export Enclosure, enclosure, cad_target
export view_factors, Float32Matrix, HMatrix, closure_error, enforce_closure!, open_fraction, RadiationModel, RadiationSolution, solve,
       flux_response, energy_balance, radiation_temperature, to_Wcm2, SIGMA_SB, albedos, conversions,
       montecarlo_incident
export WallMotion, wall_radius, wall_displacement, moving_hohlraum, Keyframes, model_at, deposition_at
export PulseShape, pulse_energy, with_pulses, drive_history, energy_budget, cone_fraction_pulses,
       write_pulse_table, read_pulse_table, retune_drive
export enclosed_radius, envelope_radius, focal_profile, enclosed_fraction, with_axes,
       BeamEnvelope, envelope_distance, ClearanceRow, ClearanceReport, clearance,
       ApertureKeepOut, aperture_clearance, SphereKeepOut, sphere_clearance,
       beam_component_power, violations, min_clearance, write_clearance_csv,
       SpotWing, PowerLawWing, UniformScatter, nif_wings, DEFAULT_CLEARANCE_THRESHOLDS
export LaserPlasmaModel, apply_lpi, Backscatter, CalibratedTransfer, linear_transfer,
       tabulated_transfer, lpi_summary, BubbleAbsorption, bubble_speed, bubble_depth,
       bubble_transmission, bubble_fractions, iris_fraction, ridge_fraction, lpi_deposit,
       izumi2021_bubble, izumi2021_absorption, young2020_radii, young2020_albedos, pre109_fig1a_wall,
       kritcher2024_capsule, kritcher2024_radius, chen2020_leh, kemp2026_sbs
export CBET, TwoZoneCBET, n221204_cbet, hybride_cbet, hybride_cbet_iraum, pulse_rise, cbet_clamp_schedule, cbet_null, n221204_fig9, n220919_fig9, n210808_fig1, pre109_fig1, hybride_shots, n221204_fig10, n221204_pulse, zylstra2022_fig2, n221204_fig9_beams, CH_IONS, ion_wave_amplitude, receiver_ratio, calibrate_cbet, cbet_gain, cbet_coupling, cbet_plasma, plasma_dispersion
export UnconvertedColor, nif_unconverted, UnconvertedRow, UnconvertedReport, unconverted_light_check,
       unconverted_deposition
export TaLISRow, TaLISReport, talis_check, write_talis_csv, TALIS_SOURCE, TALIS_OUT_OF_SCOPE, talis_requirement
export BEAM_COLUMNS, write_beam_table, read_beam_table, read_csv, write_element_table, write_moments, read_stl, read_stl_solids, write_stl, write_results
export VISRAD_SST_COLUMNS, VISRAD_PHASE_PLATES, visrad_spot, write_visrad_sst, read_visrad_sst, write_visrad_summary, read_visrad_summary, write_visrad_powers, read_visrad_powers, read_visrad_pulse_shape, read_visrad_table, visrad_legendre_fit, visrad_ylm_fit, write_visrad_kit, compare_visrad
export VisRadObject, VisRadWorkspace, read_visrad_workspace, visrad_component, visrad_target,
       write_visrad_workspace, VisRadProfile, visrad_profile, visrad_envelope,
       visrad_power_sources, visrad_deposition
export legendre, legendre_moments, legendre_weights, solid_angle, ylm_moments, ylm_weights, mode_center,
       ylm_matrix, ylm_frame, mode_index, relative_mode_power, mode_amplitude, mode_power, mode_rms, rms_nonuniformity
export OrthoCamera, render, radiance, line_of_sight_power
export backlighter_foil, line_of_sight, LineOfSightReport
export Frame, CHAMBER, to_chamber, to_frame, compose, euler_rotation, euler_frame, port_frame,
       foil_frame, part_centroid, component_frame, coords_of, point_of, FrameDef, resolve_frames,
       FramePointing, placement_sweep, SweepRow, SweepResult, write_sweep_csv
export kshell_line_energy, kshell_ce_data, kshell_ce_model, KShellCEModel, kshell_ce,
       mass_attenuation, filter_transmission, backlighter_signal, BacklighterSignal,
       henke_mass_attenuation, mirror_reflectivity
export planck_fraction, planck_fractions, spectral_edges, MBandEmission, mband_conversion,
       nif_mband, mband_power, mband_albedo, MBandSplit, SpectralSolution, spectral_solve, ngroups, band_flux,
       mband_recirculation, DanteChannel, dante1_channels, DanteView, dante1, dante2, dante, visible_area,
       AbstractDanteChannel, DanteResponse, dante1_responses, xrd_sensitivity, dante_unfold, sentinel
export GroupPointing, nparams, parameter_names, initial_parameters, parameters, apply, bounds,
       Scenario, Evaluation, evaluate, subset, flux_ylm, asymmetry_objective, influence_matrix,
       clipped_fraction, missed_fraction, power_fraction,
       flux_moments, mean_flux, symmetry_objective, irradiance_objective, optimize_pointing, PointingResult, report,
       mode_sensitivity, mode_uncertainty, sample_modes, infer_errors
export HOTSPOT_VELOCITY_PER_PERCENT, mode1_vector, hotspot_velocity, ResponseKernel, flat_kernel,
       yield_factor, alpha_factor, rke_fraction, reactivity_exponent, burn_amplification,
       shape_amplification, young2020_green, young2020_kernel, young2020_pulses, young2020_dropped_quads,
       implosion_metrics, mode1_tolerance, swing_metrics,
       velocity_uncertainty, ntof_covariance, infer_velocity_errors, VelocityPenalty
export CapsuleDrive, capsule_drive

include("geometry.jl")
include("meshops.jl")
include("bvh.jl")
include("target.jl")
include("targets.jl")
include("cadfit.jl")
include("beams.jl")
include("facilities.jl")
include("frames.jl")
include("deposition.jl")
include("storage.jl")
include("viewfactors.jl")
include("hmatrix.jl")
include("enclosure.jl")
include("radiation.jl")
include("montecarlo.jl")
include("analysis.jl")
include("corona.jl")
include("imaging.jl")
include("optimize.jl")
include("frame_pointing.jl")
include("history.jl")
include("spectral.jl")
include("xray_optics.jl")
include("dante.jl")
include("dante_response.jl")
include("dynamics.jl")
include("retune.jl")
include("clearance.jl")
include("backlighter.jl")
include("sweep.jl")
include("backlighter_signal.jl")
include("keepout.jl")
include("uncertainty.jl")
include("implosion.jl")
include("capsule_drive.jl")
include("talis.jl")
include("lpi.jl")
include("bubble.jl")
include("io.jl")
include("visrad.jl")
include("visrad_workspace.jl")
include("cbet.jl")
include("unconverted.jl")

# Time to first design: a coarse NIF hohlraum through setup, evaluation, a short
# optimization with keep-outs, clearance, the radiation solve and a drive history: what
# most scripted sessions run. Niche paths (CBET, spectra and Dante, TALIS, CAD import,
# the VisRad kit) compile at first use instead; each costs a rebuild for every user
# otherwise ("Compile latency", research notes). Set the environment variable
# BRIAREUS_NO_WORKLOAD=true to skip it (faster rebuilds while developing).
@setup_workload begin
    get(ENV, "BRIAREUS_NO_WORKLOAD", "") == "true" || @compile_workload begin
        t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                     h = 1.2, capsule_h = 0.45)
        quads = nif_quads(power = 8.0)
        p = GroupPointing(quads; fields = (:z, :r, :power))
        s = Scenario(t, p; observe = :capsule, viewfactor_kw = (; closed = true))
        x0 = initial_parameters(p)
        ev = evaluate(s, x0)
        a = flux_moments(ev)
        flux_ylm(ev)
        ev.flux
        keepout = (ApertureKeepOut(t), SphereKeepOut(t))
        for obj in (symmetry_objective(; modes = (2 => 1.0, 4 => 1.0, 6 => 0.25), drive = a[1],
                                       keepout),
                    asymmetry_objective(; keepout))
            optimize_pointing(s, obj; x0, iterations = 1, evaluations = 2)
        end
        report(devnull, ev)
        sol = solve(s.model, ev.deposition)
        energy_balance(sol)
        legendre_moments(sol.H, t, s.observed; lmax = 4)
        # the large-model path (Float32 storage, iterative solves) on the same elements
        ml = RadiationModel(t; F = s.model.F, storage = Float32, closed = false)
        energy_balance(solve(ml, ev.deposition))
        flux_response(ml, s.observed; weights = s.moment_weights)
        # hierarchical view factors (the default above Briareus.HMATRIX_MODEL[] elements)
        mh = RadiationModel(t; storage = :hmatrix, leaf = 16, closed = true)
        energy_balance(solve(mh, ev.deposition))
        flux_response(mh, s.observed; weights = s.moment_weights)
        clearance(t, quads[1:2]; targets = (:wall,), nspot = 5, nangles = 4)
        raytrace_deposition(t, quads[1:2]; nspot = 5)
        # a drive history on a wall-loss wall (its march is compiled once for every
        # `beams_at` closure)
        tw = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                      h = 1.2, capsule_h = 0.45, wall = Material(:Au; albedo = WallLossAlbedo(), conversion = 0.8))
        mw = RadiationModel(tw; F = s.model.F, closed = true)
        pulse = PulseShape([0.0, 1.0, 2.0], [0.0, 4.0, 4.0])
        drive_history(mw, s.deposition, τ -> with_pulses(quads, pulse, τ), [0.5, 1.0];
                      record = sol -> sol)
    end
end

end
