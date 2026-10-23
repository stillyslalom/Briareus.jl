# Kritcher et al., Phys. Plasmas 31, 070502 (2024) — what bears on Briareus

Source: full text read (H). doi:10.1063/5.0210904. Review of the HYBRID-E N221204 design
(first target gain > 1) and follow-ons.

## Numbers (H)

| Item | Value | Where |
|---|---|---|
| Hohlraum | Au-lined DU, 11.24 mm long × 6.4 mm diameter | Sec. II |
| Cones | inner 23°, 30°; outer 44°, 50° | Sec. II |
| Capsule | HDC, 1050 µm inner radius, ~85 µm thick; W dopant 0.62% (N221204), 0.4% (N230729) | Sec. II, Fig. 3 |
| Ablator thickness vs N210808 | +5.75 µm (N221204), +6.9 µm (N230729), +~5 µm more (N231029/N240210) | Sec. II, III |
| DT ice | 65 µm | Sec. II |
| Laser | N210808 1.92 MJ; N221204 2.05 MJ; N231029/N240210 2.2 MJ; all ~440 TW peak | Secs. I, II |
| Δλ (1ω, inner − outer) | 1.8 Å N210808; 2.5 Å N220919; 2.75 Å N221204; 2.9 Å N230729 | Secs. II, Fig. 1 |
| Pulse changes | "Peak" extended ~300 ps, trough lengthened (N221204 vs N210808); early pulse +~30 ps (N230729) | Fig. 1, Sec. II |
| Trad | > 300 eV, ~3×10⁶ K | Sec. II, Fig. 1 |
| v_imp | 380–400 km/s; 386/383/374 km/s for +0/+4/+10 µm HDC | Sec. II, Table III |
| Capsule absorbed energy | 245 / 257 / 274 kJ for +0/+4/+10 µm HDC (N230729 base, 2.05 MJ) | Table III |
| Hot-spot velocity (km/s) | N220919 51±9 (~40 no-α); N221204 36±6 (~27); N230618 107±7 (~84); N230729 104±9 (~75); N230904 82±17 (~63) | Table II |
| Hot-spot P2 (PNI, µm) | −10.5, ~0, ~0, +6, +6 (same order) | Table II |
| Mode-1 drive → velocity | applied mode one of ~0.7% reproduces burn-on >100 km/s | Fig. 5 caption |
| Mode-2 sensitivity | reducing mode two by 5 µm → ~90 % more yield (N230729 sim) | Fig. 6 |
| CBET clamp (HYDRA) | δn/n = 0.002 until start of rise, 0.008 thereafter | Appendix |
| Flux limiter | 0.15 | Appendix |
| Backscatter | energetically low on N210808/N221204; mostly inner-beam SBS early in peak | Appendix |
| N220919 | transfer limited by SBS optics-damage risk → waist drive deficit, oblate | Sec. II |
| Burn-on amplification | hot-spot velocity ratio α-on/α-off grows with yield (Fig. 5 right); burn-off ~70–75 km/s ↔ burn-on >100 | Sec. II |

## Statements about the modeling chain (H)

- Symmetry must be controlled to ~1 % time-varying flux asymmetry.
- The wall ingresses into the beam paths; this is "a main challenge" for symmetry.
- Late-time symmetry comes from a semi-empirical model (Callahan 2018, PoP 25) calibrated
  on years of data, not from first-principles LPI.
- Simulations cannot predict precisely enough; two capsule scales were planned.
- The P2 swing (time-varying P2) is the largest unknown; it is calibrated against
  in-flight radiography at ~200 µm, and for N221204 that calibration had not been done
  at the time of writing. Post-shot HYDRA under-predicted N230729 until the swing was
  reduced by hand.
- Thicker ablators implode later and sample more of the late pole-hot drive → more
  oblate, compensated by more Δλ (more waist-hot early).
- Mode-1 sources named: laser delivery imbalance, capsule thickness non-uniformity, DT
  thickness non-uniformity, non-axially symmetric diagnostic windows.
- HYDRA tuning: multipliers on input laser powers to match in-flight symmetry, in-flight
  velocity, shock velocities on two lines of sight; as-delivered per-beam powers.

## Gaps this exposes for a view-factor code

1. Static geometry: no wall motion, LEH closure or spot motion → no P2 swing. (see
   `hohlraum_dynamics.md`)
2. CBET: constant plasma state, clamp and multiplier are free knobs; the paper's own
   designers could not extrapolate 2.75 → 2.9 Å. Need measured P2-vs-Δλ data to
   calibrate against. (see `cbet_calibration.md`)
3. Absolute energetics: gold Marshak albedo on a U-Au wall, fixed η = 0.75, capsule
   albedo 0.3, open LEH; 245 kJ capsule absorbed energy is a checkable number. (see
   `energetics.md`)
4. Outputs stop at flux moments; experiments report hot-spot velocity, P2 in µm, RKE,
   yield. Need transfer functions (0.7 % → 100 km/s; 2.7 µm per 1 % cone fraction).
   (see `asymmetry_transfer.md`)
5. No benchmark case against a documented shot. (see `n221204_benchmark.md`)
