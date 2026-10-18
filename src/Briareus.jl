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
using SpecialFunctions: gamma, gamma_inc, gamma_inc_inv, dawson
import Optim
using PrecompileTools: @setup_workload, @compile_workload

export CurvedSurface, SphereSurface, CylinderSurface, DiskSurface, surface_point, element_sector
export Vec3, Mesh, npatches, add_polygon!, merge_meshes, split_mesh, transform, translate, rotation_to,
       orient_mesh, coarsen_mesh, remesh_mesh,
       cylinder_mesh, annulus_mesh, disk_mesh, sphere_mesh, cubed_sphere_mesh, rectangle_mesh, frustum_mesh
export AlbedoModel, ConstantAlbedo, TabulatedAlbedo, WallLossAlbedo, BurnthroughAlbedo, HeatFront,
       RosselandTransmission, foil_transmission, GapClosure,
       heat_front, burnthrough_time, areal_mass, optical_depth, eddington_transmission, albedo, Material, gold, ablator, BLACK,
       Component, aperture, Target, nelements, elements, component_index, element_quadrature
export hohlraum, hohlraum_counts, capsule_counts, foil, disk_foil
export WallWindow, conforming_grid, equal_area_radius
export SuperGaussianSpot, Beam, port_vector, port_angles, direction, spot_axes, remake, intensity, wavelength,
       beams_from_table, nif_quads, nif_beams, nif_spot, NIF_QUAD_NAMES, nif_port, nif_bundle, nif_bundles, omega_ports, omega_beams, omega_spot
export DepositionOperator, Deposition, deposit, raytrace_deposition, raytrace_component_power, laser_flux,
       component_power
export view_factors, closure_error, enforce_closure!, RadiationModel, RadiationSolution, solve,
       flux_response, energy_balance, radiation_temperature, to_Wcm2, SIGMA_SB, albedos
export PulseShape, pulse_energy, with_pulses, drive_history, cone_fraction_pulses,
       write_pulse_table, read_pulse_table, retune_drive
export enclosed_radius, envelope_radius, focal_profile, enclosed_fraction, with_axes,
       BeamEnvelope, envelope_distance, ClearanceRow, ClearanceReport, clearance,
       ApertureKeepOut, aperture_clearance, SphereKeepOut, sphere_clearance,
       beam_component_power, violations, min_clearance, write_clearance_csv,
       SpotWing, PowerLawWing, UniformScatter, nif_wings, DEFAULT_CLEARANCE_THRESHOLDS
export LaserPlasmaModel, apply_lpi, Backscatter, CalibratedTransfer, linear_transfer,
       tabulated_transfer, lpi_summary
export CBET, CH_IONS, ion_wave_amplitude, receiver_ratio, calibrate_cbet, cbet_gain, cbet_coupling, cbet_plasma, plasma_dispersion
export UnconvertedColor, nif_unconverted, UnconvertedRow, UnconvertedReport, unconverted_light_check
export TaLISRow, TaLISReport, talis_check, write_talis_csv, TALIS_SOURCE, TALIS_OUT_OF_SCOPE, talis_requirement
export BEAM_COLUMNS, write_beam_table, read_beam_table, read_csv, write_element_table, write_moments, read_stl, read_stl_solids, write_stl, write_results
export VISRAD_SST_COLUMNS, VISRAD_PHASE_PLATES, visrad_spot, write_visrad_sst, read_visrad_sst, write_visrad_summary, read_visrad_summary, write_visrad_powers, read_visrad_powers, read_visrad_pulse_shape, read_visrad_table, visrad_legendre_fit, visrad_ylm_fit, write_visrad_kit, compare_visrad
export VisRadObject, VisRadWorkspace, read_visrad_workspace, visrad_component, visrad_target,
       write_visrad_workspace, VisRadProfile, visrad_profile
export legendre, legendre_moments, legendre_weights, solid_angle, ylm_moments, ylm_weights,
       ylm_matrix, ylm_frame, mode_index, relative_mode_power, mode_amplitude, mode_power, rms_nonuniformity
export OrthoCamera, camera_from_port, render, radiance, line_of_sight_power
export GroupPointing, nparams, parameter_names, initial_parameters, apply, bounds,
       Scenario, Evaluation, evaluate, subset, flux_ylm, asymmetry_objective, influence_matrix,
       clipped_fraction, missed_fraction, power_fraction,
       flux_moments, symmetry_objective, irradiance_objective, optimize_pointing, report,
       mode_sensitivity, mode_uncertainty, sample_modes

include("geometry.jl")
include("meshops.jl")
include("bvh.jl")
include("target.jl")
include("targets.jl")
include("beams.jl")
include("facilities.jl")
include("deposition.jl")
include("viewfactors.jl")
include("radiation.jl")
include("analysis.jl")
include("imaging.jl")
include("optimize.jl")
include("history.jl")
include("retune.jl")
include("clearance.jl")
include("keepout.jl")
include("uncertainty.jl")
include("talis.jl")
include("lpi.jl")
include("io.jl")
include("visrad.jl")
include("visrad_workspace.jl")
include("cbet.jl")
include("unconverted.jl")

# Time to first design: a coarse NIF hohlraum through setup, evaluation, a short
# optimization with keep-outs, clearance and the radiation solve. Set the environment
# variable BRIAREUS_NO_WORKLOAD=true to skip it (faster rebuilds while developing).
@setup_workload begin
    get(ENV, "BRIAREUS_NO_WORKLOAD", "") == "true" || @compile_workload begin
        t = hohlraum(; length = 10.1, radius = 2.875, leh_radius = 1.55, capsule_radius = 0.9,
                     h = 1.2, capsule_h = 0.45)
        quads = nif_quads(power = 8.0)
        p = GroupPointing(quads; fields = (:z, :r, :power))
        s = Scenario(t, quads; observe = :capsule, parameterization = p,
                     viewfactor_kw = (; closed = true))
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
        clearance(t, quads[1:2]; targets = (:wall,), nspot = 5, nangles = 4)
        raytrace_deposition(t, quads[1:2]; nspot = 5)
    end
end

end
