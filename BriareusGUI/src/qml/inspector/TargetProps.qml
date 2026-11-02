// Properties of the target: the mesh element sizes, inserting starter assemblies, and
// adding windows and parts. A hohlraum's dimensions and wall are its part's properties.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 5

    SectionLabel { text: "Mesh"; info: "target.mesh"; Layout.topMargin: 0 }
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
        unit: "mm (walls, caps, parts)"
        tip: "Largest element side on hohlraum walls, end caps and LEHs and on the parts, unless a part sets its own"
        onEdited: (text) => Julia.br_edit_mesh("h", text)
    }
    PropRow {
        label: "Capsule element size"
        value: app.mesh_capsule_h
        unit: "mm"
        tip: "Spheres meshed as a capsule (the starters' capsule), and a starter foil when inserted. " +
             "The capsule mesh sets the accuracy of P4: its error falls as the square of this size"
        onEdited: (text) => Julia.br_edit_mesh("capsule_h", text)
    }

    SectionLabel { text: "Insert"; info: "target.insert" }
    Flow {
        Layout.fillWidth: true
        spacing: 6
        Repeater {
            model: [{ k: "hohlraum", t: "Hohlraum with capsule" }, { k: "foil", t: "Backlighter foil" },
                    { k: "shock_tube", t: "Halfraum pair with shock tube" }, { k: "opacity", t: "Opacity halfraum" }]
            Button {
                text: modelData.t
                onClicked: Julia.br_insert(modelData.k)
                ToolTip.visible: hovered
                ToolTip.text: "Insert at TCC (Insert ▸ places it in the selected part's or piece's frame)"
            }
        }
    }
    SectionLabel { text: "Add to the target" }
    Flow {
        Layout.fillWidth: true
        spacing: 6
        Button {
            text: "Window"
            enabled: app.hasHohlraum
            onClicked: Julia.br_add_window()
            ToolTip.visible: hovered
            ToolTip.text: "A closed 1 mm diagnostic window in the first hohlraum's wall at its midplane, where it fits"
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
