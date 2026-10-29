// Properties of the target: dimensions and materials, mesh element sizes, and adding
// windows and parts. The template (hohlraum, foil, parts only) starts a target: File ▸
// New target.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 5

    SectionLabel { text: "Dimensions and wall"; visible: app.targetKind !== "custom"; Layout.topMargin: 0 }
    Repeater {
        model: targetFieldModel
        delegate: PropRow {
            label: model.label
            value: model.value
            unit: model.unit === "" ? "fraction" : model.unit
            onEdited: (text) => Julia.br_edit_target(model.key, text)
        }
    }
    CheckLine {
        visible: app.targetKind === "hohlraum"
        text: "Closing entrance holes"
        checked: app.lehClosure
        onToggled: Julia.br_leh_closure(checked)
        tip: "Plasma closes each LEH to x-rays down to the open radius during the pulse; " +
                      "laser light still passes. Closed, the lip re-emits like the wall"
    }
    PropRow {
        visible: app.targetKind === "hohlraum" && app.lehClosure
        label: "Closure time"; value: app.leh_closure_time; unit: "ns"
        tip: "Middle of the closure; set it per design (e.g. from soft x-ray images)"
        onEdited: (text) => Julia.br_edit_target("leh_closure_time", text)
    }
    PropRow {
        visible: app.targetKind === "hohlraum" && app.lehClosure
        label: "Closure duration"; value: app.leh_closure_width; unit: "ns"
        onEdited: (text) => Julia.br_edit_target("leh_closure_width", text)
    }
    PropRow {
        visible: app.targetKind === "hohlraum" && app.lehClosure
        label: "Open radius"; value: app.leh_open_fraction; unit: "× LEH radius"
        tip: "Radius that stays open at peak: about 0.84 of the LEH radius on NIF (Schneider et al. 2015)"
        onEdited: (text) => Julia.br_edit_target("leh_open_fraction", text)
    }

    SectionLabel { text: "Mesh"; info: "target.mesh" }
    RowLayout {
        spacing: 8
        Label { text: "Preset"; Layout.preferredWidth: 140 }
        ComboBox {
            Layout.fillWidth: true
            // the name in the box, the accuracy and time in its list
            displayText: ["Coarse", "Standard", "Fine", "Custom element sizes"][currentIndex] ?? ""
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
