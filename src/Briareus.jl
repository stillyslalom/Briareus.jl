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
using SpecialFunctions: gamma
import Optim

export Vec3, Mesh, add_polygon!, merge_meshes, transform, translate, rotation_to,
       cylinder_mesh, annulus_mesh, disk_mesh, sphere_mesh, rectangle_mesh, frustum_mesh
export AlbedoModel, ConstantAlbedo, TabulatedAlbedo, WallLossAlbedo, albedo, Material, gold, ablator, BLACK,
       Component, aperture, Target, nelements, elements, component_index, element_quadrature
export hohlraum, foil, disk_foil
export SuperGaussianSpot, Beam, port_vector, port_angles, direction, remake, intensity,
       beams_from_table, nif_quads, nif_beams, nif_spot, omega_ports, omega_beams, omega_spot
export DepositionOperator, Deposition, deposit, raytrace_deposition, laser_flux,
       component_power
export view_factors, closure_error, enforce_closure!, RadiationModel, RadiationSolution, solve,
       flux_response, energy_balance, radiation_temperature, to_Wcm2, SIGMA_SB, albedos
export PulseShape, pulse_energy, with_pulses, drive_history
export legendre, legendre_moments, legendre_weights, solid_angle, ylm_moments, mode_power, rms_nonuniformity
export OrthoCamera, camera_from_port, render, radiance, line_of_sight_power
export GroupPointing, nparams, parameter_names, initial_parameters, apply, bounds,
       Scenario, Evaluation, evaluate, clipped_fraction, missed_fraction, power_fraction,
       flux_moments, symmetry_objective, irradiance_objective, optimize_pointing, report

include("geometry.jl")
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

end
