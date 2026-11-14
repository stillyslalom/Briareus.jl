// Properties of the target: the mesh element sizes, and inserting parts, windows,
// assemblies and frames (the Insert menu's commands). A hohlraum's dimensions and wall
// are its part's properties; ready-made targets are in File ▸ New target.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 5

    Note { text: app.targetSummary; opacity: 0.9 }
    Note { visible: app.partCount === 0; text: "The target is empty: insert a part, or start from File ▸ New target." }

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
        tip: "Largest element side, unless a part sets its own"
        onEdited: (text) => Julia.br_edit_mesh("h", text)
    }
    PropRow {
        label: "Capsule element size"
        value: app.mesh_capsule_h
        unit: "mm"
        tip: "Element size on capsule spheres; sets the accuracy of P4"
        onEdited: (text) => Julia.br_edit_mesh("capsule_h", text)
    }

    SectionLabel { text: "Insert"; info: "target.insert" }
    Flow {
        Layout.fillWidth: true
        spacing: 6
        Repeater {
            model: [{ k: "sphere", t: "Sphere" }, { k: "disk", t: "Disk" }, { k: "annulus", t: "Washer" },
                    { k: "rectangle", t: "Rectangle" }, { k: "cylinder", t: "Cylinder" }, { k: "cone", t: "Cone" },
                    { k: "hohlraum", t: "Hohlraum" }]
            Button { text: modelData.t; onClicked: Julia.br_add_part(modelData.k); tip: "Insert at TCC" }
        }
        Button { text: "STL mesh…"; onClicked: win.openStl(); tip: "Insert a part from a CAD export" }
    }
    Flow {
        Layout.fillWidth: true
        spacing: 6
        Button {
            text: "Diagnostic window"
            enabled: app.hasHohlraum
            onClicked: Julia.br_add_window()
            tip: app.hasHohlraum ? "A closed 1 mm window in the hohlraum's wall" : "Needs a hohlraum"
        }
        Button { text: "Assembly"; onClicked: Julia.br_new_assembly(); tip: "Group the selected parts" }
        Button { text: "Frame"; onClicked: Julia.br_add_frame(); tip: "A named origin and axes" }
    }
    ErrorText { text: app.targetError }
}
