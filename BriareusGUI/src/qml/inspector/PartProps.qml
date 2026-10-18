// Properties of a free part: name, size, placement (center, axis angles, rotation),
// material, sidedness, element size.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note { text: app.partKind + " · " + app.compInfo; opacity: 0.9 }
    Note {
        visible: app.selectionCount > 1
        text: app.selectionCount + " components selected: edits apply to every selected part " +
              "that has the field (shown: the first)."
    }
    Note {
        visible: app.partFile !== ""
        text: "File: " + app.partFile + "\nWindings repaired (closed solids face out of their " +
              "material, open sheets the target center) and its surface remeshed into elements of " +
              "the part's size; a part left out of the calculations is drawn from its file's facets."
    }
    PropRow {
        label: "Name"
        value: app.partName
        fieldWidth: 140
        onEdited: (text) => Julia.br_edit_part("name", text)
    }
    RowLayout {
        Label { text: "Role"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: app.roleLabels
            currentIndex: app.partRole
            onActivated: (index) => Julia.br_part_role(index)
        }
    }
    Note {
        visible: app.partKind === "VisRad object"
        text: "A VisRad CAD part's element size (below; 0: the target's) sets its remeshing; a " +
              "VisRad rectangle, sphere or hohlraum keeps VisRad's grid unless given a size. " +
              "Mounting hardware far from the target is drawn only after an import; put it back " +
              "in the calculations for the laser and clearance to see it (the drive's cavity " +
              "leaves it out of the x-ray view factors either way)."
    }
    CheckBox {
        text: "Thin film: left out of the drive's cavity (a capsule tent, an LEH window)"
        checked: app.partFilm
        onToggled: Julia.br_part_film(checked)
    }
    CheckBox {
        text: "Include in calculations (laser, x-ray drive, clearance)"
        checked: app.compIncluded
        onToggled: Julia.br_node_include(app.selection, checked)
    }
    RowLayout {
        Label { text: "Placed in"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: app.anchorLabels
            // -1 after a rejected change: show the part's frame again
            currentIndex: app.partAnchor < 0 ? 0 : app.partAnchor
            onActivated: (index) => Julia.br_part_anchor(index)
        }
    }
    RowLayout {
        Label { text: "Position in"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: ["Cartesian (x, y, z)", "Cylindrical (r, φ, z)", "Spherical (r, θ, φ)"]
            currentIndex: app.partCoords
            onActivated: (index) => Julia.br_part_coords(index)
        }
    }
    Repeater {
        model: partFieldModel
        delegate: PropRow {
            label: model.label
            value: model.value
            unit: model.unit === "" ? "fraction" : model.unit
            onEdited: (text) => Julia.br_edit_part(model.key, text)
        }
    }
    CheckBox {
        text: "Two-sided (both sides radiate and absorb)"
        checked: app.partTwoSided
        enabled: !app.partAperture
        onToggled: Julia.br_part_flag("two_sided", checked)
    }
    CheckBox {
        text: "Inward (radiates from its inside, or a flat part from its back)"
        checked: app.partInward
        onToggled: Julia.br_part_flag("inward", checked)
    }
    CheckBox {
        text: "Aperture (transparent: rays pass, its flux is what escapes)"
        checked: app.partAperture
        onToggled: Julia.br_part_flag("aperture", checked)
    }
    Note {
        text: "The axis of a part is its local z: the part is turned by the rotation about " +
              "that axis, tilted to polar angle θ and azimuth φ (as VisRad orders them), then " +
              "moved to (x, y, z) in the frame it is placed in. A disk, washer or rectangle " +
              "faces along its axis. An entrance hole's frame has its origin at the hole and z " +
              "pointing out of the hohlraum, so a part there moves with the hole and mirrors " +
              "between the two; changing the frame keeps the part where it is."
    }
    CheckBox {
        text: "Show in the 3D view"
        checked: !app.compHidden
        onToggled: Julia.br_tree_on(app.selection, checked)
    }
    ErrorText { text: app.targetError }
    Button { text: "Remove part"; onClicked: Julia.br_remove_selected() }
}
