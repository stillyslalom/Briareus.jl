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
    PropRow {
        label: "Name"
        value: app.partName
        fieldWidth: 140
        onEdited: (text) => Julia.br_edit_part("name", text)
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
              "moved to (x, y, z). A disk, washer or rectangle faces along its axis."
    }
    CheckBox {
        text: "Show in the 3D view"
        checked: !app.compHidden
        onToggled: Julia.br_tree_on(app.selection, checked)
    }
    ErrorText { text: app.targetError }
    Button { text: "Remove part"; onClicked: Julia.br_remove_selected() }
}
