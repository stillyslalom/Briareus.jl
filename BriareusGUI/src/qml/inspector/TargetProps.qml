// Properties of the target: template, dimensions and materials, mesh element sizes,
// and adding windows and parts.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8

    RowLayout {
        spacing: 8
        Label { text: "Template"; Layout.preferredWidth: 150 }
        ComboBox {
            Layout.fillWidth: true
            textRole: "text"; valueRole: "value"
            model: [{ text: "Hohlraum with capsule", value: "hohlraum" },
                    { text: "Backlighter foil", value: "foil" },
                    { text: "Custom (parts only)", value: "custom" }]
            currentIndex: Math.max(0, ["hohlraum", "foil", "custom"].indexOf(app.targetKind))
            onActivated: Julia.br_set_target_kind(currentValue)
        }
    }
    Note {
        text: app.targetKind === "foil"
              ? "The foil sits on the +x axis facing TCC; the black sphere at TCC stands in for " +
                "the main target (it absorbs and emits nothing). Sphere radius 0: none."
              : app.targetKind === "custom" ? "A target made of free parts only: add them below or from Target ▸ Add part."
              : "Cylindrical hohlraum along z with entrance holes and a capsule at TCC (capsule radius 0: none)."
    }

    SectionLabel { text: "Dimensions and wall"; visible: app.targetKind !== "custom" }
    Repeater {
        model: targetFieldModel
        delegate: PropRow {
            label: model.label
            value: model.value
            unit: model.unit === "" ? "fraction" : model.unit
            onEdited: (text) => Julia.br_edit_target(model.key, text)
        }
    }

    SectionLabel { text: "Mesh" }
    RowLayout {
        spacing: 8
        Label { text: "Preset"; Layout.preferredWidth: 150 }
        ComboBox {
            Layout.fillWidth: true
            textRole: "text"; valueRole: "value"
            model: [{ text: "Coarse · ≈0.5 % in P2, P4 · < 1 s", value: "coarse" },
                    { text: "Standard · ≈0.1 % · seconds", value: "standard" },
                    { text: "Fine · converged · tens of s", value: "fine" },
                    { text: "Custom element sizes", value: "custom" }]
            currentIndex: Math.max(0, ["coarse", "standard", "fine", "custom"].indexOf(app.targetResolution))
            onActivated: if (currentValue !== "custom") Julia.br_set_resolution(currentValue)
        }
    }
    PropRow {
        label: "Element size"
        value: app.mesh_h
        unit: app.targetKind === "foil" ? "mm (sphere, parts)" : "mm (wall, caps, parts)"
        tip: "Largest element side on the hohlraum wall, end caps and LEHs, the main-target sphere and the parts"
        onEdited: (text) => Julia.br_edit_mesh("h", text)
    }
    PropRow {
        label: app.targetKind === "foil" ? "Foil element size" : "Capsule element size"
        value: app.mesh_capsule_h
        unit: "mm"
        tip: "The surface whose drive is studied. The capsule mesh sets the accuracy of P4: its error falls as the square of this size"
        onEdited: (text) => Julia.br_edit_mesh("capsule_h", text)
    }
    Note {
        text: "Presets come from a work-precision study on the NIF hohlraum " +
              "(examples/mesh_convergence.jl): errors in the capsule P2/P0 and P4/P0 of about " +
              "0.5, 0.1 and 0.02 percentage points. A part can set its own element size."
    }

    SectionLabel { text: "Add to the target" }
    Flow {
        Layout.fillWidth: true
        spacing: 6
        Button {
            text: "Window"
            visible: app.targetKind === "hohlraum"
            onClicked: Julia.br_add_window()
            ToolTip.visible: hovered
            ToolTip.text: "A closed 1 × 1 mm diagnostic window in the wall at the midplane, where it fits"
        }
        Repeater {
            model: [{ k: "sphere", t: "Sphere" }, { k: "disk", t: "Disk" }, { k: "annulus", t: "Washer" },
                    { k: "rectangle", t: "Rectangle" }, { k: "cylinder", t: "Cylinder" }, { k: "cone", t: "Cone" }]
            Button { text: modelData.t; onClicked: Julia.br_add_part(modelData.k) }
        }
    }
    ErrorText { text: app.targetError }
    Note { text: app.targetSummary }
}
