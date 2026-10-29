// The window's explanations, one entry per topic, grouped by page: the Help screen
// (HelpDialog.qml) shows them, and an ⓘ (Info.qml) shows its topic's entry on hover and
// opens the Help screen there. Pages carry short status lines only.
.pragma library

var pages = [
    { key: "start", title: "Getting around", items: [
        { key: "start.select", title: "Selection",
          text: "Select a component, beam or analysis in the model tree or the 3D view to see and " +
                "edit its properties. Ctrl + click adds beams to the selection. Each analysis has its " +
                "settings in Properties and its results in the tab of the same name below the 3D view." },
        { key: "start.docs", title: "Documentation",
          text: "Help ▸ Documentation opens the Desktop application page of the Briareus " +
                "documentation; Help ▸ Mouse and keys lists the bindings of the 3D view." }
    ] },
    { key: "target", title: "Target", items: [
        { key: "target.hohlraum", title: "Hohlraum with capsule",
          text: "Cylindrical hohlraum along z with entrance holes and a capsule at TCC (capsule radius 0: none)." },
        { key: "target.foil", title: "Backlighter foil",
          text: "The foil sits on the +x axis facing TCC; the sphere at TCC stands in for the main " +
                "target: it absorbs the laser light and x-rays that reach it and emits nothing (albedo " +
                "and conversion 0). Sphere radius 0: none." },
        { key: "target.custom", title: "Custom (parts only)",
          text: "A target made of free parts only: add them from the target's Properties or from Target ▸ Add part." },
        { key: "target.mesh", title: "Mesh presets",
          text: "Presets come from a work-precision study on the NIF hohlraum " +
                "(examples/mesh_convergence.jl): errors in the capsule P2/P0 and P4/P0 of about " +
                "0.5, 0.1 and 0.02 percentage points. A part can set its own element size." }
    ] },
    { key: "component", title: "Template components", items: [
        { key: "component.intro", title: "Wall, caps, entrance holes, capsule, foil",
          text: "Template components are sized from the target dimensions: select Target to edit " +
                "them. However a component is drawn (invisible included), it takes part in every calculation." }
    ] },
    { key: "part", title: "Parts", items: [
        { key: "part.placement", title: "Placement and frames",
          text: "The axis of a part is its local z: the part is turned by the rotation about " +
                "that axis, tilted to polar angle θ and azimuth φ (as VisRad orders them), then " +
                "moved to (x, y, z) in the frame it is placed in. A disk, washer or rectangle " +
                "faces along its axis. An entrance hole's frame has its origin at the hole and z " +
                "pointing out of the hohlraum, so a part there moves with the hole and mirrors " +
                "between the two; the backlighter foil's frame has z along its front normal, " +
                "toward TCC. Every part has its own frame (origin at its center, z its axis): " +
                "parts placed in it and beams aimed in it move and turn with the part. Changing " +
                "the frame keeps the part where it is." },
        { key: "part.steps", title: "Placement steps",
          text: "Steps apply in order after the placement fields, each in closed form. Relations " +
                "fix degrees of freedom: Coincident (the origin at a point of a frame), On axis, " +
                "Concentric (z along a frame's axis, the origin on it), Distance (from a point, " +
                "along the line through the origin), Facing (z toward a point, by the smallest " +
                "turn), Parallel (z along an axis) and Roll (x at an angle from a frame's axis). " +
                "Offset and Rotate are relative moves and fix nothing. What no step fixes keeps " +
                "the values of the placement fields. A relation that undoes an earlier one is " +
                "refused: reorder the steps or remove one. A new step leaves the part where it " +
                "is where it can (a coincident point, a distance and a roll start at the current " +
                "values). The part follows the frames its steps refer to, as it follows the " +
                "frame it is placed in." },
        { key: "part.role", title: "Role",
          text: "A role names the part as the analyses expect it: the capsule (a sphere, or a CAD " +
                "part a sphere is fitted to) whose x-ray drive is observed, or the backlighter foil " +
                "(a flat part or a mesh) whose spot the foil objective makes uniform. Renaming a " +
                "part, or giving it a role, keeps the parts placed in its frame and the beams aimed " +
                "in it." },
        { key: "part.cad", title: "Parts from CAD files",
          text: "Windings repaired (closed solids face out of their material, open sheets the " +
                "target center) and the surface remeshed into elements of the part's size; a part " +
                "left out of the calculations is drawn from its file's facets. A session saved with " +
                "its CAD files reads the saved copy; Read from the file reads the file on disk instead." },
        { key: "part.visrad", title: "VisRad objects",
          text: "A VisRad CAD part's element size (0: the target's) sets its remeshing; a " +
                "VisRad rectangle, sphere or hohlraum keeps VisRad's grid unless given a size. " +
                "Mounting hardware far from the target is drawn only after an import; put it back " +
                "in the calculations for the laser and clearance to see it (the drive's cavity " +
                "leaves it out of the x-ray view factors either way)." }
    ] },
    { key: "beams", title: "Beams", items: [
        { key: "beams.preset", title: "Presets and beam files",
          text: "NIF presets are pointed at illustrative example cone aims; OMEGA sets at the " +
                "LEH centers or the foil. Beam tables and VisRad pointing or beam-summary files are in File." },
        { key: "beams.display", title: "Rays and volumes",
          text: "Rays run along each beam axis to the surface where it lands (diamond); the dot " +
                "is the aim point. Volumes show the 5 % contour of each spot carried along its " +
                "f-cone through best focus, each edge cut where it meets the target, so clipping on " +
                "an LEH lip or a capsule shows directly." },
        { key: "beams.pointing", title: "Pointing a selection",
          text: "Values are the mean over the selection; an edit sets that field on every beam " +
                "and leaves the others. The beam direction is fixed by its port." },
        { key: "beams.aim", title: "Aim point and aim frame",
          text: "Coordinates in the aim frame. Point every selected beam at one point, as VisRad " +
                "points a cluster at a reference point (e.g. a backlighter foil). Shown: the mean " +
                "aim. Beams aimed in a part's or a template's frame keep these coordinates in it " +
                "when it moves or turns, or when the part is renamed." },
        { key: "beams.group", title: "Groups",
          text: "Groups gather beams for pointing, power and on/off together, and are the " +
                "optimizer's per-group parameters. Type a new name to move the selection to a new group." },
        { key: "beams.quad", title: "Splitting and merging quads",
          text: "A NIF quad drawn as one beam splits into its four beams (own ports, a quarter of " +
                "the power each) to point them apart; the drive is then recomputed. Merging brings " +
                "them back into one beam (mean aim and defocus, total power)." }
    ] },
    { key: "drive", title: "Drive", items: [
        { key: "drive.intro", title: "Computing the drive",
          text: "View factors, the laser deposition operator and the capsule flux response, " +
                "once per target and beam ports (seconds on the standard mesh). Pointing, power " +
                "and on/off edits then re-evaluate the drive at once. Coloring the target by a " +
                "drive field computes it when needed. The modes and the power accounting are in " +
                "the Drive tab below the 3D view." },
        { key: "drive.model", title: "Beam model",
          text: "Briareus (f-cone): the focal spot summed over the f-cone's directions, each with " +
                "its own incidence. VisRad profile: VisRad's beam model (a super-Gaussian cut at " +
                "each beam's envelope), so the design deposits in VisRad as here; an imported " +
                "workspace sets its envelopes and focal length." },
        { key: "drive.cavity", title: "Cavity",
          text: "Only the elements that bound the free space around the point take part in the " +
                "drive: a CAD assembly's outer faces and mounting hardware drop out of the view " +
                "factors. Mark thin films modeled as solids (a capsule tent, an entrance-hole " +
                "window) in their Properties, or they seal the point into a pocket. The clearance " +
                "and TaLIS checks use the whole target." },
        { key: "drive.unconverted", title: "Unconverted light",
          text: "Beams whose 1ω or 2ω dot is on in the calculations row (model tree, beam " +
                "properties) add that light's absorbed power to the drive, converted to x-rays as the " +
                "3ω light is. Only NIF beams and quads have footprints." },
        { key: "drive.modes", title: "Capsule modes",
          text: "εₗ: rotation-invariant amplitude over all m (= |Pₗ/P₀| when axisymmetric)." },
        { key: "drive.retune", title: "Retune to a reference drive",
          text: "The reference shot fired every beam of the table with the pulses of a pulse table " +
                "(power per beam in TW against time_ns, one column per cone group, quad or beam). " +
                "The disabled beams are the ones this shot loses: the enabled beams' powers are " +
                "retuned, at each time of the table, to give the capsule the reference's drive " +
                "level and its modes through l = 4 (Compute ▸ Retune to a reference drive)." }
    ] },
    { key: "optimize", title: "Pointing optimization", items: [
        { key: "optimize.intro", title: "What is optimized",
          text: "Free pointing fields per cone, quad or beam within ± bounds around the beam " +
                "table (power bounds: a fraction of the start) and minimize an objective with " +
                "L-BFGS-B and exact gradients, in the background. Canceling keeps the best design. " +
                "The aims are free either in the chamber (z, r, tangential offset) or in each " +
                "beam's aim frame (u, v, w: on the backlighter foil's frame, u and v across its " +
                "face and w along its normal); a group shares its frame coordinates." }
    ] },
    { key: "clearance", title: "Beam clearance", items: [
        { key: "clearance.intro", title: "The check",
          text: "Worst-case clearance of each enabled beam from the hardware over pointing, " +
                "positioning, hardware and focus errors that add linearly. The envelope is each " +
                "spot's contour at the given fraction of its peak grown by the f-cone; with wings, " +
                "the far field also counts. Only hardware before a beam lands is checked. Status " +
                "thresholds are engineering defaults (fail: 10¹² W/cm² or 10⁻³ of the beam power; " +
                "review: 10⁹ W/cm² or 10⁻⁵)." },
        { key: "clearance.errors", title: "Envelope and errors",
          text: "The errors and the sampling also apply to the TaLIS checks." }
    ] },
    { key: "sweep", title: "Placement sweep", items: [
        { key: "sweep.intro", title: "The sweep",
          text: "Step a part (with the parts placed in its frame and the beams aimed in those " +
                "frames) along an axis of any frame, and check each step in the background: beam " +
                "clearance with the clearance check's settings, the backlighter's line of sight, " +
                "the power the foil captures, and optionally the pointing re-optimized with the " +
                "optimization's settings. Along an axis, the sweep translates; it does not turn." },
        { key: "sweep.steps", title: "Sweeping a placement step",
          text: "Instead of a shift along an axis, the sweep can step a number of one of the " +
                "part's placement steps (a distance, an angle, a point's or an offset's " +
                "coordinate) through a range of values. At each value the frames are resolved " +
                "again, and every part whose frame changes moves rigidly with it, turns " +
                "included: the parts placed in it, and those whose steps refer to it (a part " +
                "facing it). Beams aimed in those frames move with them. A sweep along an axis " +
                "does the same when other parts' steps refer to the moving part." },
        { key: "sweep.results", title: "Results",
          text: "A row per step: the least clearance and where, the failing and review pairs, the " +
                "line of sight, the captured power and the objective. Click the chart or a row to " +
                "select a step; double-click a row to apply it (one edit, undone as one)." }
    ] },
    { key: "talis", title: "TaLIS rules", items: [
        { key: "talis.intro", title: "The checks",
          text: "The NIF target-configuration rules for each enabled beam: 3ω capture, incidence " +
                "and clearance, 1ω clearance and unconverted light, and the TAS stay-in zone, " +
                "each judged on the worst case over the clearance error set plus the rule tolerance " +
                "(not on incidence). The report is in the TaLIS tab below the 3D view." },
        { key: "talis.tolerance", title: "Rule tolerance",
          text: "The angular tolerance becomes a lateral allowance r·tan(tolerance) at each hardware " +
                "point, with r measured from the chosen point. The error set (pointing, positioning, " +
                "hardware, focus) and the sampling are the clearance check's." },
        { key: "talis.report", title: "The report",
          text: "Left: one line per note (the tolerance, the error set, the checks not evaluated, " +
                "the requirements out of scope), the whole note on hover, and the per-rule summary " +
                "with each rule's worst check. Right: every check, colored by status with its margin." }
    ] },
    { key: "view", title: "3D view", items: [
        { key: "view.cut", title: "Cutaway angles",
          text: "Azimuth about z runs from +x toward +y, about x from +y toward +z, about y from +z toward +x." }
    ] }
]

var _entries = {}
var _page = {}
for (var i = 0; i < pages.length; ++i)
    for (var j = 0; j < pages[i].items.length; ++j) {
        _entries[pages[i].items[j].key] = pages[i].items[j]
        _page[pages[i].items[j].key] = i
    }

// The explanation of a topic ("" for an unknown one).
function text(key) { return _entries[key] ? _entries[key].text : "" }
// The index into `pages` of a topic's page (0 for an unknown one).
function pageOf(key) { return _page[key] !== undefined ? _page[key] : 0 }
