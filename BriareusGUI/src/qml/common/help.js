// The window's explanations, one entry per topic, grouped by page: the Help screen
// (HelpDialog.qml) shows them, and an ⓘ (Info.qml) shows its topic's entry on hover and
// opens the Help screen there. Pages carry short status lines only.
.pragma library

var pages = [
    { key: "start", title: "Getting around", items: [
        { key: "start.select", title: "Selection",
          text: "Select a part, component, frame, beam or analysis in the model tree or the 3D view " +
                "to see and edit its properties. Ctrl + click adds a node of the same kind to the " +
                "selection, or takes it out (beams with beams, parts and components with parts and " +
                "components, frames with frames); in the model tree Shift + click selects the range " +
                "from the selected node. An edit then applies to everything selected. Each analysis " +
                "has its settings in Properties and its results in the tab of the same name below " +
                "the 3D view." },
        { key: "start.docs", title: "Documentation",
          text: "Help ▸ Documentation opens the Desktop application page of the Briareus " +
                "documentation; Help ▸ Mouse and keys lists the bindings of the 3D view." }
    ] },
    { key: "target", title: "Target", items: [
        { key: "target.insert", title: "Starters",
          text: "A target is made of parts. Insert ▸ adds a starter assembly in the selected part's " +
                "or piece's frame (TCC with nothing selected), its names made unique (hohlraum_2): " +
                "Hohlraum with capsule, a hohlraum part along z with a capsule sphere at its center " +
                "that holds the observed role; Backlighter foil, a foil 6 mm along x facing the frame's " +
                "origin and a sphere there standing in for the main target (it absorbs the laser " +
                "light and x-rays that reach it and emits nothing); Halfraum pair with shock tube, " +
                "two halfraums (hohlraum parts with one end open) facing each other across a " +
                "cylinder package whose end faces close them, the top face observed; Opacity " +
                "halfraum, a halfraum with a sample foil, observed, in its open end; Frustraum with " +
                "capsule, a hohlraum part with the frustraum profile (the stated size of a published " +
                "NIF frustraum) and an observed capsule sphere at its center. File ▸ New " +
                "target starts a target with one starter. The dimensions are illustrative " +
                "mock-ups." },
        { key: "target.mesh", title: "Mesh presets",
          text: "Presets come from a work-precision study on the NIF hohlraum " +
                "(examples/mesh_convergence.jl): errors in the capsule P2/P0 and P4/P0 of about " +
                "0.5, 0.1 and 0.02 percentage points. A part can set its own element size." }
    ] },
    { key: "assembly", title: "Assemblies", items: [
        { key: "assembly.intro", title: "The model tree",
          text: "An assembly groups parts and other assemblies, at any depth; parts and " +
                "assemblies can also sit at the top level. Each assembly has a frame, placed like a " +
                "part's: its members are placed in it by default, so moving or sweeping the " +
                "assembly moves them and the beams aimed in their frames. Membership and " +
                "placement are separate: a member may be placed in any frame (a capsule in its " +
                "hohlraum's). An assembly's box shows or hides everything in it (a dash: partly " +
                "shown); Include in calculations puts all of it in or leaves all of it out, and a " +
                "part left out gives up its role (observed, foil, wall) until it is put back. The " +
                "analyses see the parts alone, never the tree." },
        { key: "assembly.commands", title: "Tree commands",
          text: "Right-click in the model tree: New assembly (of the selected parts and " +
                "assemblies, or an empty one from Target), Move into assembly (the members stay " +
                "where they are), Ungroup (the assembly goes, its members stay, in its parent), " +
                "Rename (whatever is placed in it, aimed in it or a member of it follows) and " +
                "Remove with its contents. An imported VisRad workspace or multi-body STL file " +
                "arrives as one assembly named after the file; Insert with an assembly selected " +
                "puts the starter in it." }
    ] },
    { key: "component", title: "Pieces", items: [
        { key: "component.intro", title: "Wall, caps, entrance holes, windows",
          text: "A hohlraum part's pieces are sized from the part: select the part to edit them. " +
                "However a piece is drawn (invisible included), it takes part in every calculation." },
        { key: "component.role", title: "A piece's role",
          text: "A hohlraum part's wall and caps hold the Wall role (the drive history's and their " +
                "cavity's wall temperature) and its other pieces none. A piece's own role replaces " +
                "that default: a cap made Observed is the observed surface, a flat face with its face " +
                "map and flux statistics, and no longer one of its cavity's walls; a wall section or a " +
                "cap made None leaves the walls; a closed window made Wall joins them. One observed " +
                "surface at a time: a part or piece that held the role gives it up. Entrance holes, " +
                "their lips and open windows are apertures and take no role. Set it in the piece's " +
                "properties or from its row's right-click menu in the model tree." }
    ] },
    { key: "frames", title: "Frames", items: [
        { key: "frames.intro", title: "User frames",
          text: "A user frame is a named frame with no geometry, placed like a part: an origin " +
                "at (x, y, z) in the frame it is placed in, its z axis tilted to polar angle θ and " +
                "azimuth φ, turned by the rotation about z, then moved by its placement steps. " +
                "Parts can be placed in it and beams aimed in it, and they follow it; a sweep can " +
                "step it along an axis. A frame in use (parts or frames placed in it, or steps " +
                "that refer to it) cannot be removed; beams aimed in a removed frame keep their " +
                "aims in the chamber frame. Add frame places the new frame in the selected " +
                "frame, TCC with none selected." },
        { key: "frames.port", title: "Port frames",
          text: "A port frame has its origin at TCC and z toward the port (VisRad's diagnostic " +
                "port axis), so the port's beams travel along −z; x is the port's polar " +
                "direction θ̂ turned by the roll toward φ̂ (an assumption: VisRad's x for port " +
                "frames is not documented). The ports are those of the beam table, NIF quads at " +
                "the angles of the NIF port table. Only the roll can be edited; a frame placed in " +
                "a port frame carries offsets along the port axis." },
        { key: "frames.transform", title: "Coordinate transform",
          text: "Tools ▸ Coordinate transform gives the coordinates in one frame of a point or " +
                "a direction given in another, each in Cartesian, cylindrical (r, φ, z about the " +
                "frame's z) or spherical (r, θ, φ) coordinates, mm and degrees. A direction turns " +
                "with the frames and does not move with their origins. Create frame adds a user " +
                "frame placed in the To frame: at the point with the From frame's axes, or, for a " +
                "direction, at the To frame's origin with z along it." },
        { key: "frames.axes", title: "Frame axes in the view",
          text: "The 3D view draws the selected frame's axes, and those of the frames whose " +
                "axes are switched on (the frame's page, or its menu in the model tree): x red, y " +
                "green, z blue, with the frame's name at the origin. Drawing only." }
    ] },
    { key: "part", title: "Parts", items: [
        { key: "part.placement", title: "Placement and frames",
          text: "The axis of a part is its local z: the part is turned by the rotation about " +
                "that axis, tilted to polar angle θ and azimuth φ (as VisRad orders them), then " +
                "moved to (x, y, z) in the frame it is placed in. A disk, washer or rectangle " +
                "faces along its axis. A hohlraum part's entrance-hole frames have their origin at the " +
                "hole and z pointing out of the hohlraum, so a part there moves with the hole and " +
                "mirrors between the two; a foil's frame has z along its front normal. Every part " +
                "has its own frame (origin at its center, z its axis): " +
                "parts placed in it and beams aimed in it move and turn with the part. Changing " +
                "the frame keeps the part where it is." },
        { key: "part.triad", title: "Placement triad",
          text: "The 3D view draws a triad on the selected part, piece or assembly, in the axes of " +
                "the frame it is placed in (x red, y green, z blue). Dragging an arrow moves the part " +
                "along that axis in whole steps of the Move-along row's \"by\" field; dragging a ring " +
                "turns it about that axis, through the part's origin, in 15° steps. Ctrl held at the " +
                "press snaps to a tenth of either. The moved parts and the beams aimed in their frames " +
                "are previewed while dragging; the release makes one edit, undone as one, the edit " +
                "the Move-along row's buttons make. A press released without moving is a click, so " +
                "the beams behind a handle stay selectable." },
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
          text: "A role marks the part for the analyses; the part keeps its name. The observed " +
                "surface is the one whose x-ray drive the drive, the optimization and the drive " +
                "history analyze, one at a time: a capsule (a sphere, or a CAD part a sphere is " +
                "fitted to: modes, θ–φ map, keep-out) or a flat face (a one-sided rectangle or disk, " +
                "a VisRad rectangle, or a hohlraum part's cap given the role as a piece: the face's " +
                "map and flux statistics). A " +
                "backlighter foil is a flat part or a mesh. A wall is a re-emitting cavity surface: " +
                "the drive history's wall temperature, and a surface the clearance and TaLIS checks " +
                "land the beams on (a hohlraum part's wall and caps hold the role for the drive history). Any number of " +
                "parts can be foils or walls; " +
                "the active foil (Drive ▸ Active foil, or the part's own check) is the one whose " +
                "spot statistics, map and objective the drive and the optimization take, and whose " +
                "line of sight and captured power the sweep checks. The clearance and TaLIS checks " +
                "take every foil as a surface the beams land on. Renaming a part keeps the parts " +
                "placed in its frame, the beams aimed in it, and what its components were given " +
                "(left out, display modes, the cut, the checks' roles, thin films, the sweep)." },
        { key: "part.hohlraum", title: "Hohlraum parts",
          text: "A hohlraum along its part's z axis (a cylinder, a frustraum or sections: see " +
                "Profile), in pieces: the wall (one piece per section), a cap at each end that is " +
                "not open, each end's entrance hole, and the windows. An end is an entrance hole in " +
                "its cap, a closed cap, or open (no cap: a halfraum's end, closed by another part " +
                "or letting the radiation out). The wall and caps hold the wall role by default; a " +
                "piece's Role changes it (a cap observed as a flat face). The clearance and TaLIS " +
                "checks keep the caps as hardware, the entrance holes' edges. Each end has a " +
                "frame (origin at the end's center, z pointing out of the hohlraum, so one aim or " +
                "placement mirrors between the two), and each window one at its center on the " +
                "wall: beams aim and parts are placed in them, and they follow the hohlraum." },
        { key: "part.profile", title: "Hohlraum profiles",
          text: "Cylinder: one radius. Frustraum: two cones meeting at the waist (z = 0), the " +
                "waist radius and each end's radius (the bottom's 0: the top's). Sections: a " +
                "radius at each break, a cone between unequal radii, a cylinder between equal " +
                "ones and a flat step where two breaks share a height; the end breaks sit at " +
                "±length/2. Each section is its own wall piece with its exact surface. Windows " +
                "are cut into cylindrical sections only. The frustraum starter has the stated " +
                "size of a published NIF frustraum (Baker et al. 2023), its end radius derived " +
                "from the 23° wall angle: a mock-up of the shape, not that platform. The NIF " +
                "presets aim in its holes with a pointing derived for it (see Presets and beam files)." },
        { key: "part.windows", title: "Diagnostic windows",
          text: "A diagnostic window is cut into a hohlraum part's wall, centered at azimuth φ from the " +
                "part's +x and at z along its axis: round (a disk of the given diameter on the unrolled " +
                "wall, the default) or a rectangle (a width in mm of arc by a height along z). Closed, it " +
                "is an opaque wall patch with its own x-ray albedo and laser → x-ray conversion; open, it " +
                "is a hole that laser light and x-rays leave through. Windows are cut into cylindrical " +
                "sections only. Add window (the part's properties) places a closed 1 mm window at the " +
                "midplane, at the first of φ = 90°, 270°, 0°, 180°, … where it overlaps no other. Each " +
                "window has a frame at its center on the wall for aims and placements, and holds no role " +
                "unless given one (a closed window made Wall joins its cavity's walls)." },
        { key: "part.leh_closure", title: "Closing entrance holes",
          text: "With Closing entrance holes on, each entrance hole stays open to x-rays out to the open " +
                "radius (a fraction of its radius), and the lip beyond closes: it lets its incident flux " +
                "escape until the closure time, then across the closure duration turns into plasma that " +
                "re-emits with the wall's albedo. Laser light and lines of sight pass the lip throughout. " +
                "The drive history follows the closure through the pulse. Defaults: 0.84 of the LEH " +
                "radius, the measured radius at peak (Schneider et al. 2015; research notes, \"Closing " +
                "entrance holes\"); closure from 5 ns over 0.5 ns, an engineering example to set per design." },
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
        { key: "beams.intro", title: "The beam set",
          text: "The Beams node holds the beam set: the facility's ports, grouped by cone (or by user " +
                "group), with NIF quads or single beams. Its properties load a preset, switch every beam " +
                "on or off and set how the beams are drawn and colored. The Beam table tab below the 3D " +
                "view edits pointing, defocus and power per cone, quad or beam. File imports and exports " +
                "beam tables, VisRad SST pointing and VisRad beam summaries." },
        { key: "beams.nodes", title: "Groups, quads and beams",
          text: "Under Beams the model tree has a node per group, and under each its quads (NIF) or " +
                "beams. A node's box switches its beams on or off (a dash: some are off); its 2 × 3 dots " +
                "are the 1ω, 2ω and 3ω light in the calculations (top row) and drawn (bottom row), a pale " +
                "dot on for some of its beams. Selecting a node shows its properties: the pointing (the " +
                "mean over its beams), the aim point and its frame, the group and, for a NIF quad, " +
                "splitting and merging; an edit sets that field on every selected beam and leaves the " +
                "others. Right-click: Look along beam, Move to a new group, Split or Merge a quad." },
        { key: "beams.table", title: "Beam table",
          text: "The Beam table tab below the 3D view: a row per cone, quad or beam (Rows), with the " +
                "number of beams n, the port angles θ / φ, the frame the aims are held in (mixed: the " +
                "row's beams aim in different frames), and the aim z, aim r, tangential offset and " +
                "defocus (mm) and the power (TW per beam), edited in place: an edit sets the field on " +
                "every beam of the row. The box switches the row's beams on or off (pale: some are off). " +
                "A click on a heading sorts the rows by that column, again descending, a third time back " +
                "to the table's order; sorting changes the display only, not the beams' order. " +
                "Selecting a row selects its node in the model tree and the 3D view. Problems found in an " +
                "imported VisRad file show above the table." },
        { key: "beams.preset", title: "Presets and beam files",
          text: "The hohlraum presets aim in the entrance holes' frames: the upper hemisphere's " +
                "beams in the hole facing up, the lower hemisphere's in the one facing down, so a " +
                "moved hohlraum carries its beams. NIF beams keep the offsets of an illustrative " +
                "example cone pointing; in a frustraum's holes, a pointing derived from it that " +
                "puts each cone's spot at the same fraction of the wall height where the hole " +
                "allows (the inner cones' spots land farther past the waist than on the " +
                "cylinder). " +
                "OMEGA beams aim at the holes' centers, and a hemisphere " +
                "with no hole facing it is switched off (the opacity halfraum). The OMEGA " +
                "backlighter set aims at the active foil's center, in its frame. An engineering " +
                "default, not a facility's pointing. Beam tables and VisRad pointing or " +
                "beam-summary files are in File." },
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
                "aim. Beams aimed in a part's or a piece's frame keep these coordinates in it " +
                "when it moves or turns, or when the part is renamed. Choosing another aim frame keeps " +
                "the aims where they are; Aim at the frame's center points them at its center." },
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
        { key: "drive.walls", title: "Hohlraum wall",
          text: "The radiation temperature of the area-weighted mean x-ray flux incident on one " +
                "hohlraum part's pieces that hold the wall role (by default its wall and caps, but " +
                "not a cap observed as a flat face; never its entrance holes or windows). A " +
                "target with several cavities (a halfraum pair) has one per hohlraum part: the " +
                "panel opens on the one that holds the observed surface, and the drive history " +
                "charts each one. Parts given the wall role outside a hohlraum part (a shield) " +
                "are not a cavity of their own." },
        { key: "drive.modes", title: "Capsule modes",
          text: "εₗ: rotation-invariant amplitude over all m (= |Pₗ/P₀| when axisymmetric)." },
        { key: "drive.retune", title: "Retune to a reference drive",
          text: "The reference shot fired every beam of the table with the pulses of a pulse table " +
                "(power per beam in TW against time_ns, one column per cone group, quad or beam). " +
                "The disabled beams are the ones this shot loses: the enabled beams' powers are " +
                "retuned, at each time of the table, to restore the reference's drive: a capsule's " +
                "level and modes through l = 4, a flat observed face's flux over the face (mean, " +
                "RMS, peak-to-mean), or, laser only, the active foil's 3ω irradiance within the spot " +
                "radius (Compute ▸ Retune to a reference drive)." }
    ] },
    { key: "time", title: "Time", items: [
        { key: "time.pulses", title: "Pulses",
          text: "The Time tab's first chart: the laser power of each beam group against time (the sum " +
                "over its enabled beams, TW), colored as in the model tree, and the total. Pulse table… " +
                "fires the beams with a pulse table (power per beam in TW against time_ns, one column " +
                "per cone group, quad or beam; disabled beams stay off); N221204 shape fires them with " +
                "N221204's delivered pulse shape (Pak et al. 2024), each beam reaching its table power at " +
                "the peak. Run history solves the drive at each time of the pulses with the drive's view " +
                "factors, for the other charts." },
        { key: "time.drive", title: "Drive history",
          text: "From the latest history: the capsule's radiation temperature Tr and each recorded " +
                "diagnostic's (eV, dashed), and a₂/a₀ and a₄/a₀ (%). For a flat observed face, its Tr " +
                "(from its area-mean flux) and, on the right axis, its RMS nonuniformity and its peak " +
                "and lowest flux relative to its mean (%, dashed). A target with several cavities " +
                "charts each hohlraum part's wall temperature." },
        { key: "time.diagnostics", title: "Diagnostics chart",
          text: "From the latest history: each recorded diagnostic's Tr (eV) and the intensity toward " +
                "it (TW/sr, dashed). A diagnostic is recorded when it is ticked under Record before the " +
                "history runs." },
        { key: "time.energy", title: "Energy",
          text: "From the latest history: the laser energy delivered by each time and where it went " +
                "(kJ), and the unconverted 1ω and 2ω light's part of the delivered energy (dotted) when " +
                "that light takes part." },
        { key: "time.retune", title: "Retune chart",
          text: "From the latest retune (Compute ▸ Retune to a reference drive): each group's power on " +
                "the reference pulses (all beams, dashed) and on the retuned ones (enabled beams), and " +
                "the capsule's mode ε₁ before and after the retune (%). For a flat face or a foil spot: " +
                "its RMS nonuniformity before, after and in the reference, and its mean's departure from " +
                "the reference's before and after (%)." },
        { key: "time.cursor", title: "Time cursor",
          text: "A click on a chart, a drag along it, or the slider puts the time cursor at a step of " +
                "the history: the 3D view is colored by the drive at that step, and the images at the " +
                "cursor are taken there. Now returns to the drive at the beam table's powers. The " +
                "cursor is off while the history is out of date." },
        { key: "time.images", title: "Images list",
          text: "Images (n) lists the images the latest history holds of its recorded imagers: for " +
                "each, the image at the time cursor's step (detected radiance, W/cm²/sr; the last step " +
                "without a cursor) and each frame of its recording (detected fluence per steradian, " +
                "J/cm²/sr: the integral over the shot or against a gate). An image opens in the image " +
                "view on its region with signal and stays until another is chosen: a drive edit does " +
                "not replace it." },
        { key: "time.irradiance", title: "Irradiance checks",
          text: "Run checks computes the irradiance of the unconverted 1ω and 2ω light on each part " +
                "of the drive's target through the pulses: the drive's power ratios on the beams " +
                "whose 1ω or 2ω light takes part, plus the prepulse (Drive properties). The chart " +
                "shows each part's peak against the plasma threshold, 10⁸ W/cm² by default: the NIF " +
                "User Guide's guideline for the onset of weakly ionized plasma on a metal (§4.4), " +
                "taken as the end of specular reflection. It is an engineering default for metals, " +
                "not a facility criterion. Specular passes while every part stays below it; " +
                "pre-plasma, on a target with a foil, fails when the unconverted light reaches it on " +
                "the active foil before the 3ω light does. A dot is gray when the checks are out of date." }
    ] },
    { key: "diagnostics", title: "Diagnostics", items: [
        { key: "diagnostics.record", title: "Recording diagnostics",
          text: "Record ▸ Diagnostics (the Time tab) lists every diagnostic definition found, by kind: " +
                "the built-ins, those in the search path's folders and those saved with the session. " +
                "Each ticked one is recorded in the next drive history and shows as a chip beside the " +
                "list; ✕ stops recording it. A definition saved with the session that differs from the " +
                "search path's of the same id is marked * and is the one used; ↺ uses the search " +
                "path's instead. A definition that does not build is dimmed, the reason in its tooltip." },
        { key: "diagnostics.folders", title: "Folders and the search path",
          text: "Diagnostics ▸ Folders… shows the search path: the built-in definitions are read " +
                "first, then each folder's TOML files in order, and a definition with the id of an " +
                "earlier one replaces it. Folders added there are saved for later sessions; the " +
                "BRIAREUS_DIAGNOSTICS environment variable adds others (marked environment). Read again " +
                "reads every folder again, and Problems lists each definition's errors and warnings with " +
                "its file. A session saves the definitions it records, so it opens without the folders." }
    ] },
    { key: "optimize", title: "Pointing optimization", items: [
        { key: "optimize.intro", title: "What is optimized",
          text: "Free pointing fields per cone, quad or beam within ± bounds around the beam " +
                "table (power bounds: a fraction of the start) and minimize an objective with " +
                "L-BFGS-B and exact gradients, in the background. Canceling keeps the best design. " +
                "The aims are free either in the chamber (z, r, tangential offset) or in each " +
                "beam's aim frame (u, v, w: on the backlighter foil's frame, u and v across its " +
                "face and w along its normal); a group shares its frame coordinates." },
        { key: "optimize.weights", title: "Objective weights",
          text: "The capsule objectives weigh each mode's (Pₗ/P₀)², about 10⁻² on an untuned capsule. " +
                "The observed face's uniformity weight multiplies its squared RMS nonuniformity over the " +
                "start's, so the term starts at the weight on any face: an opacity sample uniform to 0.01 % " +
                "and a package face at 5 % weigh alike. The clipped, missed and direct-light penalties " +
                "weigh 10 on squared power fractions (10 % missed: 0.1) and the held drive 1 on the squared " +
                "relative departure. A weight of 1 puts the shape first; a smaller one gives the light lost " +
                "more say." },
        { key: "optimize.results", title: "The Optimization tab",
          text: "The tab below the 3D view: the objective's values, the modes εₗ before and after, the " +
                "clipping and the capsule clearance of the result, the foil spot or the observed face, and each free " +
                "parameter's start, optimized value and bounds (a value at a bound is marked). Apply to " +
                "the beam table writes the result into the table. \"Inputs changed since\" marks a result " +
                "computed for other beams, a drive or settings. Set it up and run it from Pointing " +
                "optimization in the model tree's Analyses, or the toolbar's Optimize…; it needs a " +
                "current drive." }
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
          text: "The errors and the sampling also apply to the TaLIS checks." },
        { key: "clearance.results", title: "The Clearance tab",
          text: "The tab below the 3D view: a summary line, the closest beam to each component (its " +
                "worst-case and nominal clearance, the beam's power and the status), and a row per beam " +
                "and component: the status, the worst-case and nominal clearance (mm), the power (TW) " +
                "and the peak intensity on the component (W/cm²). Export CSV… writes the rows. " +
                "\"Inputs changed since\" marks a report computed for other beams, target or settings. " +
                "Set the check up and run it from Beam clearance in the model tree's Analyses, or the " +
                "toolbar's Clearance…." }
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
        { key: "view.cutaway", title: "Cutaway",
          text: "Cutaway (toolbar) removes part of the target from the drawing to show inside. Drawing " +
                "only: every calculation sees the whole target, and a click in the removed part picks " +
                "nothing there. A sector removes a slice of the given width about the x, y or z axis, " +
                "centered on an azimuth, optionally between bounds along the axis; a plane removes the " +
                "side its normal points to, through a point and moved by the offset along the normal. A " +
                "plane can contain a line of sight (from TCC or the selected component, toward θ, φ or a " +
                "beam port) and turn about it. Right-click the Cutaway button for the common cuts, Cut " +
                "settings… and Components cut: by default the hohlraum parts' shells (the wall and the " +
                "caps of the ends that are not open), not the entrance holes, lips or windows; that " +
                "default follows the target as it changes." },
        { key: "view.display", title: "Display modes",
          text: "A component is drawn opaque, semi-transparent, as a wireframe or not at all (its " +
                "properties, or its row's right-click menu in the model tree). Drawing only: an " +
                "invisible component still takes part in every calculation it is included in; Include " +
                "in calculations is separate. The model tree's boxes show or hide a part or an assembly " +
                "with everything in it." },
        { key: "view.cut", title: "Cutaway angles",
          text: "Azimuth about z runs from +x toward +y, about x from +y toward +z, about y from +z toward +x." },
        { key: "view.image", title: "Maps and images",
          text: "The capsule map, the foil face and the x-ray and diagnostic images zoom with the wheel " +
                "and to a box dragged with the left button, and pan with the right. A diagnostic image " +
                "opens zoomed to its region with signal: Ctrl + Shift + click shows the whole field, " +
                "Ctrl + click returns to the opening view." }
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
