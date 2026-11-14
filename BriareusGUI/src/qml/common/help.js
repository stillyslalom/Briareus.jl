// The window's Help: one entry per topic, grouped by page. The Help window
// (HelpWindow.qml) shows the pages; an ⓘ (Info.qml) shows its topic's `tip` on hover and
// opens the Help window at the topic on a click.
//
// Writing an entry: `tip` is one short line for the hover. `text` is Markdown, written for
// someone who is stuck: short sentences, one idea each, lists for steps and choices, the
// control's name in bold where it helps to find it. Say what a control does and when to
// use it; leave derivations and provenance to the documentation.
.pragma library

// A text block written with the code's indentation: each line trimmed, so that Markdown
// sees paragraphs and list items at the left margin.
function md(s) { return s.split("\n").map(l => l.trim()).join("\n").trim() }

var pages = [
    { key: "start", title: "Getting around", items: [
        { key: "start.layout", title: "The window",
          tip: "Model tree, 3D view, Properties, result tabs",
          text: md(`
            - **Model** (left): the target's parts, the frames, the beams and the analyses.
            - **3D view** (center): the target and the beams.
            - **Properties** (right): the settings of whatever is selected.
            - **Result tabs** (bottom): the beam table and the results of each analysis.

            Nothing has to be done in a fixed order. An analysis that no longer matches the design is marked out of date; run it again when needed.`) },
        { key: "start.select", title: "Selecting",
          tip: "Click to select; Ctrl or Shift + click for several; Esc clears",
          text: md(`
            Click a row of the model tree, or an object in the 3D view, to select it. Its settings appear in Properties.

            - **Ctrl + click** adds a row to the selection or takes it out.
            - **Shift + click** in the tree selects the range from the selected row.
            - Rows of one kind select together: parts with parts, beams with beams, frames with frames.
            - **Esc**, or a click on empty space in the 3D view, clears the selection.

            The 3D view highlights everything selected. An edit in Properties, a right-click command or Remove applies to all of it.`) },
        { key: "start.keys", title: "Mouse and keys",
          tip: "Navigation and keyboard shortcuts",
          text: md(`
            **3D view**

            - Wheel: zoom toward the pointer.
            - Left or middle drag: orbit about the selection, or about the point under the pointer.
            - Right drag, Shift + drag or Ctrl + middle drag: pan.
            - Click: select. Ctrl + click: add or remove. Esc: clear the selection.
            - Double click: fit the view to what was clicked.
            - F: fit everything. Arrow keys: orbit by 15°.

            **Maps and images**

            - Wheel: zoom. Left drag: zoom to a box. Right drag: pan.
            - Ctrl + click: back to the opening view.
            - Ctrl + Shift + click: the whole field.

            **Keyboard**

            - Ctrl + Z, Ctrl + Y: undo, redo.
            - Delete: remove the selected parts, assemblies, frames and windows.
            - Ctrl + S, Ctrl + O: save, open a session.
            - Ctrl + R: compute the drive.
            - F1: Help.`) },
        { key: "start.turntable", title: "Turntable rotation",
          tip: "Keeps z up while orbiting",
          text: md(`
            **View ▸ Turntable rotation** keeps the chamber's z axis up while orbiting. It is on by default.

            Switched off, the view turns freely about the screen's axes and can roll. Switching it back on levels the view: z points up again.`) },
        { key: "start.files", title: "Sessions and exports",
          tip: "Saving the design and writing results",
          text: md(`
            - **File ▸ Save session** writes the whole design: target, beams, settings and the CAD files it reads (.briar).
            - **Save session without CAD files** keeps the CAD files by path instead.
            - **Export results (CSV)** writes the drive to a folder: fields per element, power per part and per beam, the capsule modes and a copy of the session.

            With a current drive history, Export results also writes the history: the drive, the energy and each recorded diagnostic against time (history.csv), the fields per element at every step, and the imagers' images. An image open in the image view is written too.`) },
        { key: "start.docs", title: "Documentation",
          tip: "Opens the documentation's page on this window",
          text: md(`
            **Documentation** opens the Desktop application page of the Briareus documentation: a first session, every panel, and the methods behind the numbers.`) }
    ] },
    { key: "target", title: "Target", items: [
        { key: "target.new", title: "New target",
          tip: "A ready-made target with its facility's beams",
          text: md(`
            **File ▸ New target** replaces the target with a ready-made one and loads the beams of the facility it belongs to.

            **NIF**

            - *Hohlraum with capsule*: a cylindrical hohlraum with two entrance holes and a capsule.
            - *Frustraum with capsule*: the same with a frustraum wall.
            - *Halfraum pair with shock tube*: two halfraums facing each other across a package.
            - *Opacity halfraum*: one halfraum with a sample in its open end.

            **OMEGA**

            - *Backlighter foil*: a foil 6 mm from a sphere that stands in for the main target.

            **Empty** starts with no parts; build the target from the Insert menu.

            The dimensions are illustrative mock-ups, not a facility's designs. Imported beams are kept.`) },
        { key: "target.insert", title: "Inserting parts",
          tip: "Insert adds a part, window, assembly or frame",
          text: md(`
            The **Insert** menu adds to the target:

            - a **part**: sphere, disk, washer, rectangle, cylinder, cone, a bare hohlraum, or a mesh from an STL file;
            - a **diagnostic window** in the selected hohlraum's wall;
            - an **assembly** of the selected parts;
            - a **frame**.

            A new part appears at the chamber center (TCC), selected. Set its size and place in Properties.`) },
        { key: "target.remove", title: "Removing",
          tip: "Remove takes out everything selected",
          text: md(`
            Select parts, assemblies, frames or windows, then right-click ▸ **Remove** (or press Delete). Everything selected is removed in one step, and Undo brings it back.

            - An assembly is removed with its contents.
            - Parts that were placed in a removed part's frame stay where they are. They are then placed in the next frame up (a capsule left behind by its hohlraum is placed at TCC).
            - Beams that aimed in a removed frame keep their aim points.
            - A hohlraum's pieces (wall, caps, entrance holes) cannot be removed; hide them or change the hohlraum's ends instead.`) },
        { key: "target.mesh", title: "Mesh",
          tip: "Element sizes: accuracy against time",
          text: md(`
            The mesh preset sets the element sizes of the whole target.

            - **Coarse**: about 0.5 % error in the capsule P2 and P4, under a second.
            - **Standard**: about 0.1 %, a few seconds.
            - **Fine**: converged, tens of seconds.

            The capsule element size sets the accuracy of P4. A part can set its own element size in its properties.`) }
    ] },
    { key: "assembly", title: "Assemblies", items: [
        { key: "assembly.intro", title: "Assemblies",
          tip: "A group of parts that moves together",
          text: md(`
            An assembly groups parts and other assemblies. It has its own frame, and its members are placed in that frame, so moving the assembly moves them and the beams aimed at them.

            - The assembly's box shows or hides everything in it. A dash means partly shown.
            - **Include in calculations** puts all of it in or leaves all of it out.
            - A member can still be placed in another frame, such as a capsule in its hohlraum's.

            The analyses see the parts only; the tree is for organizing.`) },
        { key: "assembly.commands", title: "Tree commands",
          tip: "Right-click a row for its commands",
          text: md(`
            Right-click a row of the model tree:

            - **New assembly of the selection** groups the selected parts.
            - **Move into assembly…** moves them to another assembly or to the top level. They stay where they are in space.
            - **Ungroup** removes the assembly and keeps its members.
            - **Rename…** renames a part, assembly or frame. Everything that refers to it follows.
            - **Remove** removes everything selected.

            An imported VisRad workspace or a multi-body STL file arrives as one assembly named after the file.`) }
    ] },
    { key: "component", title: "Hohlraum pieces", items: [
        { key: "component.intro", title: "Pieces",
          tip: "Wall, caps, entrance holes and windows",
          text: md(`
            A hohlraum part is drawn as pieces: wall sections, caps, entrance holes and windows. Their sizes come from the part, so select the part to change them.

            A piece can be shown, hidden or drawn semi-transparent on its own. However it is drawn, it stays in the calculations.`) },
        { key: "component.role", title: "A piece's role",
          tip: "Wall, observed face, or none",
          text: md(`
            By default the wall and the caps are **Wall**: they count toward the hohlraum's wall temperature. The other pieces have no role.

            - **Observed** on a cap makes it the observed surface: a flat face with its own flux map. Only one surface is observed at a time.
            - **None** on a wall section or cap takes it out of the wall temperature.
            - **Wall** on a closed window adds it to the walls.

            Entrance holes and open windows are openings and take no role.`) }
    ] },
    { key: "frames", title: "Frames", items: [
        { key: "frames.intro", title: "Frames",
          tip: "A named origin and axes to place parts and aim beams in",
          text: md(`
            A frame is a named origin with axes and no geometry. Parts can be placed in it and beams aimed in it; they follow when it moves.

            - **Insert ▸ Frame** adds one in the selected frame, or at TCC.
            - Place it like a part: position, tilt of its z axis (θ, φ), rotation about z, then placement steps.
            - Every part also has a frame of its own, with z along the part's axis.

            Removing a frame leaves what was placed in it where it is.`) },
        { key: "frames.port", title: "Port frames",
          tip: "Origin at TCC, z toward a port",
          text: md(`
            A port frame has its origin at TCC and its z axis toward a beam port, so the port's beams travel along −z.

            Only its roll about z can be edited. To work at a distance along the port axis, place another frame in it.`) },
        { key: "frames.transform", title: "Coordinate transform",
          tip: "Convert a point or direction between frames",
          text: md(`
            **Tools ▸ Coordinate transform** converts a point or a direction from one frame to another, in Cartesian, cylindrical or spherical coordinates (mm and degrees).

            A direction turns with the frames but ignores their origins.

            **Create frame** adds a frame at the point (or with z along the direction).`) },
        { key: "frames.axes", title: "Frame axes",
          tip: "x red, y green, z blue",
          text: md(`
            The 3D view draws the axes of the selected frame: x red, y green, z blue, with its name at the origin.

            **Draw axes in the 3D view** keeps a frame's axes on when it is not selected.`) }
    ] },
    { key: "part", title: "Parts", items: [
        { key: "part.placement", title: "Placing a part",
          tip: "Position, tilt and rotation in a frame",
          text: md(`
            A part is placed in a frame (**Placed in**): TCC, a user frame, or another part's frame.

            - **x, y, z**: its center in that frame. **Position in** switches to cylindrical or spherical coordinates.
            - **θ, φ**: the direction its axis (its own z) is tilted to.
            - **Rotation**: its turn about that axis.

            A disk, washer or rectangle faces along its axis. Choosing another frame keeps the part where it is.

            To move a part by a step, use **Move along** or drag the triad in the 3D view.`) },
        { key: "part.triad", title: "Moving with the triad",
          tip: "Drag an arrow to move, a ring to turn",
          text: md(`
            The 3D view draws a triad on the selected part, along the axes of the frame it is placed in.

            - Drag an **arrow** to move the part along that axis.
            - Drag a **square** to move it in that plane.
            - Drag a **ring** to turn it about that axis.

            Moves snap to the **Move along** step and turns to 15°. Hold Ctrl when pressing for a tenth of that. Each drag is one undo step.`) },
        { key: "part.steps", title: "Placement steps",
          tip: "Relations and moves applied after the placement fields",
          text: md(`
            Placement steps position a part relative to other frames. They apply in order, after the placement fields.

            **Relations** hold the part to something:

            - *Coincident*: its origin at a point of a frame.
            - *On axis*: its origin on a frame's axis.
            - *Concentric*: on the axis, with its z along it.
            - *Distance*: at a distance from a point.
            - *Facing*: its z toward a point.
            - *Parallel*: its z along an axis.
            - *Roll*: its x at an angle from an axis.

            **Offset** and **Rotate** are plain moves on top.

            Whatever the steps leave free is set by the placement fields. The note under the steps says what is still free. A step that contradicts an earlier one is refused: reorder or remove one.`) },
        { key: "part.role", title: "Role",
          tip: "What the analyses take the part as",
          text: md(`
            A role tells the analyses what a part is.

            - **Observed**: the surface whose x-ray drive is analyzed. A sphere is treated as a capsule (modes, θ–φ map); a flat one-sided part as a face (flux map). One at a time.
            - **Foil**: a backlighter foil, analyzed for its laser spot. With several foils, the **active foil** is the one the drive and the optimization report on.
            - **Wall**: a re-emitting surface that counts toward the wall temperature.
            - **None**: plain hardware.

            A part left out of the calculations gives up its role until it is put back.`) },
        { key: "part.hohlraum", title: "Hohlraum parts",
          tip: "A hohlraum along its z axis, in pieces",
          text: md(`
            A hohlraum part lies along its own z axis. Each end is one of:

            - **Entrance hole**: a cap with a hole of the LEH radius.
            - **Closed**: a full cap.
            - **Open**: no cap at all, as on a halfraum closed by a package.

            Each end has a frame with z pointing out of the hohlraum, and each window has one at its center. Beams aim in these frames and follow the hohlraum when it moves.`) },
        { key: "part.profile", title: "Wall profile",
          tip: "Cylinder, frustraum or sections",
          text: md(`
            - **Cylinder**: one radius.
            - **Frustraum**: two cones that meet at the waist. Set the waist radius and the end radius.
            - **Sections**: a list of breaks, each with a height and a radius. Unequal radii give a cone between them; two breaks at one height give a step.

            Windows can only be cut into cylindrical sections.`) },
        { key: "part.windows", title: "Diagnostic windows",
          tip: "A patch or hole in the hohlraum wall",
          text: md(`
            A diagnostic window sits in the hohlraum wall at azimuth φ and height z.

            - **Round** or rectangular.
            - **Closed**: a wall patch with its own albedo and conversion.
            - **Open**: a hole that light and x-rays leave through.

            **Insert ▸ Diagnostic window** adds a closed 1 mm window at the midplane, at the first free azimuth.`) },
        { key: "part.leh_closure", title: "Closing entrance holes",
          tip: "Plasma closes the holes to x-rays during the pulse",
          text: md(`
            With **Closing entrance holes** on, each hole keeps an open radius for x-rays, and the lip outside it closes during the pulse:

            - before the closure time, the lip lets x-rays out;
            - over the closure duration, it turns into plasma that re-emits like the wall.

            Laser light passes throughout. The drive history follows the closure in time.

            The defaults are 0.84 of the LEH radius, closing from 5 ns over 0.5 ns. Set them per design.`) },
        { key: "part.cad", title: "Parts from CAD files",
          tip: "STL meshes, remeshed to the element size",
          text: md(`
            **Insert ▸ Mesh from STL file** adds a part from a CAD export. Its surface is remeshed to the part's element size. Set **Scale** to the file's unit in mm (25.4 for inches).

            A session saved with its CAD files reads the saved copies. **Read from the file** reads the file on disk again, and **Relink…** points the part at another file.`) },
        { key: "part.visrad", title: "VisRad objects",
          tip: "Objects imported from a VisRad workspace",
          text: md(`
            Objects from an imported VisRad workspace keep VisRad's grid unless given an element size.

            Mounting hardware far from the target is only drawn after an import. Tick **Include in calculations** for the laser and the clearance check to see it.`) }
    ] },
    { key: "beams", title: "Beams", items: [
        { key: "beams.intro", title: "The beam set",
          tip: "The facility's beams, by cone, quad and beam",
          text: md(`
            The **Beams** row holds every beam, grouped by cone, then by quad (NIF), then by beam.

            - Its properties load a preset, switch all beams on or off and set how beams are drawn.
            - The **Beam table** tab edits pointing and power row by row.
            - **File** imports and exports beam tables and VisRad pointing files.`) },
        { key: "beams.facility", title: "Facility",
          tip: "Sets the beams, the diagnostics and the checks",
          text: md(`
            **Beams ▸ Facility** sets which facility the design is for.

            - **NIF**, **OMEGA**: loads the facility's beams. The Diagnostics tab lists its instruments, and the TaLIS rules apply to NIF only.
            - **OMEGA EP**, **Other**: no built-in beams. The beam table is emptied; import one with **File ▸ Import beam table**.

            Imported beams are kept when the facility changes.`) },
        { key: "beams.nodes", title: "Groups, quads and beams",
          tip: "A row edits all the beams under it",
          text: md(`
            Selecting a group, quad or beam shows its pointing and power. An edit sets that value on every beam under the row.

            - The **box** switches the beams on or off. A dash means some are off.
            - The **dots** are the 1ω, 2ω and 3ω light: top row in the calculations, bottom row drawn.

            Clicking a beam in the 3D view selects its quad while the quad's row is closed in the tree. Open the quad's row to pick single beams.`) },
        { key: "beams.table", title: "Beam table",
          tip: "Pointing and power per cone, quad or beam",
          text: md(`
            One row per cone, quad or beam (**Rows**). Edit a cell to set that value on every beam of the row.

            - **aim z, aim r, tang.**: the aim point relative to the quad's port azimuth (mm).
            - **defocus**: how far best focus lies beyond the aim point (mm).
            - **TW/quad** or **TW/beam**: the power. NIF rows give the power per quad.
            - **frame**: the frame the aims are held in.

            Click a heading to sort. Sorting changes the display only.`) },
        { key: "beams.preset", title: "Presets",
          tip: "A facility's beams, aimed for this target",
          text: md(`
            A preset replaces the beams with a facility's set, aimed for the current target.

            - Hohlraum presets aim into the entrance holes and follow the hohlraum when it moves.
            - The OMEGA backlighter set aims at the active foil.

            The pointing is an illustrative default, not a shot's.`) },
        { key: "beams.display", title: "Rays and volumes",
          tip: "How beams are drawn",
          text: md(`
            - **Rays**: each beam's axis to where it lands (diamond). The dot is the aim point.
            - **Volumes**: each beam's envelope, cut where it meets the target. Clipping on a lip or capsule shows directly.
            - **Hidden**: no beams drawn.

            Drawing never changes a calculation.`) },
        { key: "beams.pointing", title: "Pointing",
          tip: "Aim offsets, defocus and power of the selection",
          text: md(`
            The fields show the mean over the selected beams. An edit sets that field on all of them and leaves the others alone.

            A beam's direction is fixed by its port; pointing moves where it aims.`) },
        { key: "beams.aim", title: "Aim point and aim frame",
          tip: "Point the selection at one spot, in a frame",
          text: md(`
            **x, y, z** point every selected beam at one spot, in the coordinates of the **aim frame**.

            Beams aimed in a part's frame follow the part when it moves. Choosing another aim frame leaves the aim where it is. **Aim at the frame's center** moves it to the frame's origin.`) },
        { key: "beams.group", title: "Groups",
          tip: "Named sets of beams edited together",
          text: md(`
            A group is a named set of beams: a cone by default. Groups are pointed, powered and switched together, and the optimizer can treat each as one unit.

            To make a new group, type a name and press **Move**.`) },
        { key: "beams.quad", title: "Quads and beams",
          tip: "Split quads for per-beam deposition",
          text: md(`
            A NIF quad can be held as one beam or as its four beams.

            - **Split into beams** gives each beam its own port, so the laser deposition is computed per beam.
            - **Merge into quads** goes back to one beam per quad.

            Both work on a quad, a group, or all beams (right-click the row).

            The quad stays the unit either way:

            - Its row sets one aim point for its beams. Aim r and the tangential offset are measured in the quad's port azimuth.
            - Its power is given per quad and shared by its beams.
            - Pulse tables with cone or quad columns give the power per quad.
            - The optimizer per cone or per quad keeps each quad's beams on one aim point.

            To point a quad's beams apart, open its row and edit single beams. The drive must be computed again after a split or a merge.`) }
    ] },
    { key: "drive", title: "Drive", items: [
        { key: "drive.intro", title: "Computing the drive",
          tip: "Laser deposition and x-ray drive on the target",
          text: md(`
            **Compute ▸ Drive** (Ctrl + R) builds the view factors and the laser deposition for the target and the beam ports. This takes seconds on the standard mesh.

            After that, pointing and power edits update the drive at once. Changing the target or the ports marks it out of date.

            The **Drive** tab shows the capsule modes, the wall temperature and where the laser power lands.`) },
        { key: "drive.model", title: "Beam model",
          tip: "How a beam deposits its light",
          text: md(`
            - **Briareus (f-cone)**: sums the focal spot over the directions of the beam's cone.
            - **VisRad profile**: VisRad's beam model, for comparing a design with VisRad.`) },
        { key: "drive.cavity", title: "Cavity",
          tip: "Limit the drive to the surfaces around a point",
          text: md(`
            With **Compute the drive on the cavity around a point**, only the surfaces that bound the free space around that point take part in the x-ray drive. Outer faces and mounting hardware of a CAD assembly drop out.

            Thin films drawn as solids (a capsule tent, an LEH window) would seal the point in: tick **Thin film** on those parts.

            The clearance and TaLIS checks always use the whole target.`) },
        { key: "drive.unconverted", title: "Unconverted light",
          tip: "1ω and 2ω light in the drive",
          text: md(`
            Set a **1ω** or **2ω power ratio** (a fraction of each beam's 3ω power) to include that light. Then switch it on per beam with the top-row dots in the model tree.

            The **prepulse** fields add light that does not follow the 3ω pulse; it enters the drive history and the irradiance checks.

            Only NIF beams have unconverted-light footprints.`) },
        { key: "drive.walls", title: "Hohlraum wall",
          tip: "The wall's radiation temperature",
          text: md(`
            The wall temperature is the radiation temperature of the mean x-ray flux on a hohlraum's wall pieces.

            A target with several hohlraums has one per hohlraum: choose which to show. The drive history charts each.`) },
        { key: "drive.modes", title: "Capsule modes",
          tip: "εₗ: the drive asymmetry at mode l",
          text: md(`
            **εₗ** is the amplitude of mode l of the x-ray flux on the capsule, relative to the mean, over all azimuthal orders.

            **Pₗ/P₀** is the axisymmetric (Legendre) part. The two agree for an axisymmetric drive.`) },
        { key: "drive.retune", title: "Retune to a reference drive",
          tip: "Recover a reference shot's drive with fewer beams",
          text: md(`
            **Compute ▸ Retune to a reference drive** adjusts the pulses of the enabled beams so that they reproduce a reference shot's drive.

            1. Disable the beams this shot loses.
            2. Give the reference shot's pulse table (all beams firing).
            3. Choose one power per quad, cone or beam, and a power limit.
            4. Press **Retune**. The result is written as a new pulse table.

            **Target** chooses what to match: the reference's drive, or a uniform drive at the reference's level or a set one. **Fire the retuned pulses** loads the result as the session's pulses.

            Cone and quad columns give the power per quad, in the reference table and in the new one. If the Time tab is set to **TW per beam**, both are per beam.`) }
    ] },
    { key: "time", title: "Time", items: [
        { key: "time.pulses", title: "Pulses",
          tip: "The laser power of each group against time",
          text: md(`
            **Pulses ▾** sets the pulses the beams fire:

            - **Pulse table…** reads a CSV file: power in TW against time_ns, one column per cone, quad or beam.
            - **Published shapes** scales a published NIF pulse to the beams, each reaching its own power at the peak.

            A cone's or a quad's column gives the power **per quad**, shared by the quad's beams, so the pulses fire the same whether quads are split or not. Switch the box beside Pulses ▾ to **TW per beam** if the table gives each beam's power.

            The chart shows the power of each group and the total. **Run history** then computes the drive at each time.`) },
        { key: "time.drive", title: "Drive history",
          tip: "Tr and symmetry through the pulse",
          text: md(`
            After **Run history**: the radiation temperature of the observed surface, and the drive symmetry (a₂/a₀ and a₄/a₀ for a capsule; the nonuniformity for a flat face).

            A target with several hohlraums also shows each wall temperature.`) },
        { key: "time.energy", title: "Energy",
          tip: "Where the laser energy has gone by each time",
          text: md(`
            The laser energy delivered by each time, and where it went: the walls, the observed surface, the entrance holes and the windows (kJ).`) },
        { key: "time.retune", title: "Retune chart",
          tip: "Powers and symmetry before and after a retune",
          text: md(`
            From the latest retune: each group's power on the reference pulses (dashed) and on the retuned ones, and the asymmetry before and after.`) },
        { key: "time.cursor", title: "What the 3D view shows",
          tip: "Beam powers, or a time of the drive history",
          text: md(`
            The 3D view colors the target by one of two drives:

            - **Beam powers**: the drive at the powers set on the beams (Beams tree, Beam table). Pointing and power edits show at once.
            - **Drive history**: the drive at one time of the latest history. Move the slider, or click or drag on the chart, to choose the time.

            The history view is unavailable while the history is out of date.`) },
        { key: "time.irradiance", title: "Irradiance checks",
          tip: "Unconverted light against the plasma threshold",
          text: md(`
            **Run checks** computes the irradiance of the unconverted 1ω and 2ω light on each part through the pulse and compares it with the plasma threshold (Drive properties, 10⁸ W/cm² by default).

            - **Specular** passes while every part stays below the threshold.
            - **Pre-plasma** fails when unconverted light reaches the threshold on the active foil before the 3ω light does.

            The threshold is an engineering default for metals. A gray dot means the checks are out of date.`) }
    ] },
    { key: "diagnostics", title: "Diagnostics", items: [
        { key: "diagnostics.record", title: "Recording diagnostics",
          tip: "Add instruments, then run the drive history",
          text: md(`
            The **Diagnostics** tab records synthetic diagnostics in the drive history.

            1. Choose a category: x-ray power, x-ray imagers or radiography.
            2. Press **Add ▾** and pick an instrument. It appears as a chip.
            3. Run the drive history (Time tab).

            The list shows the instruments of the session's facility. **✕** on a chip stops recording it. A dimmed instrument cannot be built; its tooltip says why.`) },
        { key: "diagnostics.config", title: "Configurations",
          tip: "Right-click a chip to change its configuration",
          text: md(`
            An instrument can have several configurations, such as Dante's channel sets from different years and shots. It is listed once and recorded in one configuration.

            Right-click or double-click its chip to choose another configuration.

            A configuration marked * was saved with the session and differs from the one in the diagnostic directories. **Use the directories' definition** switches to the latter.`) },
        { key: "diagnostics.chart", title: "X-ray power",
          tip: "Each recorded diagnostic's Tr and intensity",
          text: md(`
            For each recorded power diagnostic: its radiation temperature (eV) and the intensity toward it (TW/sr, dashed), against time.

            Click or drag on the chart to move the time cursor.`) },
        { key: "diagnostics.images", title: "X-ray imagers",
          tip: "The images a history holds of its imagers",
          text: md(`
            Each recorded imager gives:

            - an image **at the cursor**: what it detects at the time the cursor shows;
            - its **recording**: one time-integrated image, or one frame per gate.

            Click an entry to show it. **Open in the main view** shows it large, where it zooms and pans.`) },
        { key: "diagnostics.folders", title: "Diagnostic directories",
          tip: "Where diagnostic definitions are read from",
          text: md(`
            Diagnostics are defined in TOML files. **Diagnostics ▸ Diagnostic directories** lists the directories they are read from, after the built-in ones.

            - **Add directory…** adds one and keeps it for later sessions.
            - **Refresh** reads every definition again.
            - **Problems** lists the errors found, with their files.

            A definition with the same id as an earlier one replaces it. A session saves the definitions it records, so it opens without the directories.`) }
    ] },
    { key: "optimize", title: "Pointing optimization", items: [
        { key: "optimize.intro", title: "What is optimized",
          tip: "Pointing and power, within bounds, for an objective",
          text: md(`
            The optimization varies the pointing fields that are ticked as free, within ± bounds around the current values, to minimize an objective. It runs in the background; canceling keeps the best design so far.

            - Parameters are per cone, quad or beam.
            - Aims are free either in the chamber (z, r, tangential offset) or in each beam's aim frame.

            It needs a current drive.`) },
        { key: "optimize.weights", title: "Objective weights",
          tip: "How much each term counts",
          text: md(`
            Each mode's weight multiplies its squared amplitude. A weight of 0 ignores the mode.

            The penalties for clipped, missed and direct light use a weight of 10 on the squared power fraction. A weight of 1 on the drive puts its shape first; smaller weights give lost light more say.`) },
        { key: "optimize.results", title: "Results",
          tip: "Before and after, and each parameter",
          text: md(`
            The **Optimization** tab shows the objective, the modes before and after, the clipping, and each free parameter with its bounds. A value at a bound is marked.

            **Apply to the beam table** writes the result into the beams.`) }
    ] },
    { key: "clearance", title: "Beam clearance", items: [
        { key: "clearance.intro", title: "The check",
          tip: "Worst-case clearance of each beam from the hardware",
          text: md(`
            The check finds how close each enabled beam comes to the hardware in the worst case over the pointing, positioning, hardware and focus errors.

            Give each component a role: a **target** the beams land on, **hardware** to stay clear of, or none. Only hardware a beam passes before it lands is checked.

            The fail and review thresholds are engineering defaults.`) },
        { key: "clearance.errors", title: "Envelope and errors",
          tip: "The beam envelope and the error set",
          text: md(`
            The envelope is each spot's contour at the given fraction of its peak, grown along the beam's cone. With **wings**, the far field counts too.

            The errors add linearly. The TaLIS checks use the same errors and sampling.`) },
        { key: "clearance.results", title: "Results",
          tip: "The closest beam per component, and every pair",
          text: md(`
            The **Clearance** tab lists the closest beam to each component, then every beam–component pair: status, worst-case and nominal clearance, power and peak intensity.

            **Export CSV…** writes the rows.`) }
    ] },
    { key: "sweep", title: "Placement sweep", items: [
        { key: "sweep.intro", title: "The sweep",
          tip: "Step a part along an axis and check each step",
          text: md(`
            A sweep moves a part along an axis in steps and checks each position in the background: beam clearance, the backlighter's line of sight, the power the foil captures, and optionally a re-optimized pointing.

            Parts placed in the moving part's frame, and beams aimed in it, move with it.`) },
        { key: "sweep.steps", title: "Sweeping a placement step",
          tip: "Vary one number of a placement step instead",
          text: md(`
            Instead of moving along an axis, a sweep can vary one number of a placement step: a distance, an angle or a coordinate. Everything that depends on the part follows at each value, turns included.`) },
        { key: "sweep.results", title: "Results",
          tip: "One row per step",
          text: md(`
            Each row is one step: the least clearance and where, failing pairs, the line of sight, the captured power and the objective.

            Click a row or the chart to select a step. Double-click a row, or press **Apply step**, to move the part there.`) }
    ] },
    { key: "talis", title: "TaLIS rules", items: [
        { key: "talis.intro", title: "The checks",
          tip: "NIF's target-configuration rules, per beam",
          text: md(`
            The TaLIS checks apply NIF's target-configuration rules to each enabled beam: 3ω capture, incidence and clearance, 1ω clearance and unconverted light, and the TAS stay-in zone.

            Each is judged in the worst case over the clearance check's errors plus the rule tolerance.`) },
        { key: "talis.tolerance", title: "Rule tolerance",
          tip: "A rotation of the component, as a lateral allowance",
          text: md(`
            The rule tolerance is an angle: the component rotated off nominal in the least favorable direction, as when a target is inserted rotated or assembled off-square.

            It becomes a lateral allowance r·tan(tolerance) at each hardware point, with r measured from the chosen pivot. For a component assembled off-square, set the pivot to its mount.

            It applies to the clearances and the capture margin, not to the 70° incidence limit.`) },
        { key: "talis.report", title: "The report",
          tip: "Rules on the left, every check on the right",
          text: md(`
            - **Left**: each rule with its pass, review and fail counts and its worst check, then the notes.
            - **Right**: every check with its status and margin.

            **Export CSV…** writes the checks.`) },
        { key: "talis.scope", title: "Out of scope",
          tip: "Requirements these checks do not cover",
          text: md(`
            These requirements of the TCC Configuration Requirements are not checked, because the hardware involved is not modeled:

            - positioner clearance inside and outside the CIVS field of view;
            - targets to the TAS, and the TAS withdrawal envelope;
            - diagnostics to the target shroud's swept volume;
            - ARC beams crossing a NIF beam footprint on a target feature.`) }
    ] },
    { key: "view", title: "3D view", items: [
        { key: "view.cutaway", title: "Cutaway",
          tip: "Hide part of the target to see inside",
          text: md(`
            **Cutaway** removes part of the target from the drawing. Calculations still see the whole target.

            - A **sector** removes a slice about an axis: set its center azimuth and width.
            - A **plane** removes one side: set its normal, a point and an offset.

            Right-click the Cutaway button for common cuts, **Cut settings…** and the list of components that are cut.`) },
        { key: "view.cut", title: "Cutaway angles",
          tip: "Azimuth runs from +x toward +y about z",
          text: md(`
            The azimuth about z runs from +x toward +y; about x from +y toward +z; about y from +z toward +x.`) },
        { key: "view.display", title: "Display modes",
          tip: "Opaque, semi-transparent, wireframe or hidden",
          text: md(`
            A component is drawn opaque, semi-transparent, as a wireframe, or not at all. Set it in the component's properties or its right-click menu.

            This only changes the drawing. Whether a component takes part in the calculations is set by **Include in calculations**.`) },
        { key: "view.image", title: "Maps and images",
          tip: "Wheel zooms, drag zooms to a box, right-drag pans",
          text: md(`
            The capsule map, the foil face and the x-ray images zoom with the wheel or a dragged box and pan with the right button.

            A diagnostic image opens on its region with signal. Ctrl + Shift + click shows the whole field; Ctrl + click returns.`) }
    ] }
]

var _entries = {}
var _page = {}
for (var i = 0; i < pages.length; ++i)
    for (var j = 0; j < pages[i].items.length; ++j) {
        _entries[pages[i].items[j].key] = pages[i].items[j]
        _page[pages[i].items[j].key] = i
    }

// Whether `key` is a topic.
function has(key) { return _entries[key] !== undefined }
// A topic's one-line hint for a hover: its `tip`, else its title ("" for an unknown one).
function hint(key) { return _entries[key] ? (_entries[key].tip || _entries[key].title) : "" }
// The explanation of a topic, Markdown ("" for an unknown one).
function text(key) { return _entries[key] ? _entries[key].text : "" }
// The index into `pages` of a topic's page (0 for an unknown one).
function pageOf(key) { return _page[key] !== undefined ? _page[key] : 0 }
