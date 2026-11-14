// Properties of a free part: name, size, placement (the frame it is placed in, center,
// axis angles, rotation, placement steps), a move along any frame's axis, material,
// sidedness, element size.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    readonly property bool hohlraum: app.partKind === "Hohlraum"
    RowLayout {
        Layout.fillWidth: true
        Note { text: app.partKind + " · " + app.compInfo; opacity: 0.9 }
        Info { topic: "part.visrad"; visible: app.partKind === "VisRad object"; Layout.alignment: Qt.AlignTop }
        Info { topic: "part.hohlraum"; visible: hohlraum; Layout.alignment: Qt.AlignTop }
    }
    Note {
        visible: app.selectionCount > 1
        text: app.selectionCount + " components selected: edits apply to every selected part " +
              "that has the field (shown: the first)."
    }
    RowLayout {
        visible: app.partFile !== ""
        Layout.fillWidth: true
        Cell { text: "File: " + app.partFile; opacity: 0.72; Layout.fillWidth: true }
        Info { topic: "part.cad" }
    }
    Note {
        visible: app.partFile !== "" && app.partCopy
        text: "Read from the copy saved in the session" +
              (app.partSourceExists ? "." : " (the file is not on this computer).")
    }
    RowLayout {
        visible: app.partFile !== ""
        Button {
            text: "Relink…"
            Tip { visible: parent.hovered; text: "Read this file's parts from another file" }
            onClicked: relinkDialog.open()
        }
        Button {
            text: "Read from the file"
            visible: app.partCopy
            enabled: app.partSourceExists
            Tip { visible: parent.hovered; text: "Read the file on disk instead of the saved copy" }
            onClicked: Julia.br_part_from_source()
        }
    }
    FileDialog {
        id: relinkDialog
        title: "Read the part's CAD from another file"
        fileMode: FileDialog.OpenFile
        nameFilters: ["CAD files (*.stl *.STL *.vrw)", "All files (*)"]
        onAccepted: Julia.br_relink_part(selectedFile.toString())
    }
    PropRow {
        label: "Name"
        value: app.partName
        fieldWidth: 140
        tip: "Letters, digits and _"
        onEdited: (text) => Julia.br_edit_part("name", text)
    }
    RowLayout {
        visible: !hohlraum
        Label { text: "Role"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: app.roleLabels
            currentIndex: app.partRole
            onActivated: (index) => Julia.br_part_role(index)
            Tip {
                visible: parent.hovered
                text: "What the analyses take the part as"
            }
        }
        Info { topic: "part.role" }
    }
    CheckLine {
        // with two or more foils: the one the foil analyses take (Drive ▸ Active foil too)
        visible: app.partRole === 2 && app.foilNames.length > 1
        text: "Active foil (the foil analyses take it)"
        checked: app.foilNames[app.activeFoil] === app.partName
        tip: "The foil the drive and the optimization report on"
        onToggled: {
            if (checked) Julia.br_active_foil(app.foilNames.indexOf(app.partName))
            checked = Qt.binding(() => app.foilNames[app.activeFoil] === app.partName)
        }
    }
    CheckLine {
        visible: !hohlraum
        text: "Thin film: left out of the drive's cavity (a capsule tent, an LEH window)"
        checked: app.partFilm
        onToggled: Julia.br_part_film(checked)
        tip: "Keeps a thin film from sealing the drive's cavity"
    }
    CheckLine {
        text: "Include in calculations (laser, x-ray drive, clearance)"
        checked: app.compIncluded
        onToggled: Julia.br_node_include(app.selection, checked)
        tip: "Left out: drawn, but in no calculation"
    }
    RowLayout {
        Label { text: "Placed in"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: app.frameLabels
            // -1 after a rejected change: show the part's frame again
            currentIndex: app.partFrame < 0 ? 0 : app.partFrame
            onActivated: (index) => Julia.br_part_frame(index)
            tip: "The frame the part is placed in; it stays where it is"
        }
        Info { topic: "part.placement" }
    }
    MoveAlongRow { what: "the part" }
    Button {
        text: "Sweep along this axis…"
        Layout.leftMargin: 114
        onClicked: Julia.br_sweep_part()
        tip: "Step the part along this axis and check each step"
    }
    RowLayout {
        Label { text: "Position in"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: ["Cartesian (x, y, z)", "Cylindrical (r, φ, z)", "Spherical (r, θ, φ)"]
            currentIndex: app.partCoords
            onActivated: (index) => Julia.br_part_coords(index)
            tip: "Coordinates of the position fields"
        }
    }
    // a hohlraum part's wall profile: one radius, a frustraum, or a radius per break
    RowLayout {
        visible: hohlraum
        Label { text: "Profile"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: ["Cylinder", "Frustraum", "Sections"]
            currentIndex: app.partProfile
            onActivated: (index) => Julia.br_edit_part("profile", String(index))
            Tip {
                visible: parent.hovered
                text: "The shape of the wall"
            }
        }
        Info { topic: "part.profile" }
    }
    ColumnLayout {
        id: sectionTable
        visible: hohlraum && app.partProfile === 2
        spacing: 2
        readonly property int count: app.partSectionZ.length
        TableHeader {
            id: sectionHead
            columns: [{ t: "Break", w: 40 }, { t: "z (mm)", w: 72 }, { t: "r (mm)", w: 72 }, { t: "", w: 30 }]
        }
        Repeater {
            model: sectionTable.count
            delegate: RowLayout {
                spacing: sectionHead.spacing
                readonly property bool inner: index > 0 && index < sectionTable.count - 1
                Label { text: String(index + 1); Layout.preferredWidth: sectionHead.widths[0] }
                TextField {
                    Layout.preferredWidth: sectionHead.widths[1]
                    text: app.partSectionZ[index] || ""
                    enabled: inner
                    selectByMouse: true
                    onEditingFinished: if (text !== app.partSectionZ[index]) Julia.br_section_edit(index, "z", text)
                    Tip { visible: parent.hovered && !inner; text: "An end break: at ±length/2 (edit the length)" }
                }
                TextField {
                    Layout.preferredWidth: sectionHead.widths[2]
                    text: app.partSectionR[index] || ""
                    selectByMouse: true
                    onEditingFinished: if (text !== app.partSectionR[index]) Julia.br_section_edit(index, "r", text)
                }
                Button {
                    text: "−"
                    Layout.preferredWidth: sectionHead.widths[3]
                    enabled: inner
                    onClicked: Julia.br_section_remove(index)
                    Tip { visible: parent.hovered; text: "Remove this break, joining the sections beside it" }
                }
            }
        }
        Button {
            text: "Add break"
            onClicked: Julia.br_section_add()
            Tip { visible: parent.hovered; text: "Split the longest section at its middle (the shape is unchanged)" }
        }
    }
    Repeater {
        model: partFieldModel
        delegate: PropRow {
            label: model.label
            value: model.value
            unit: model.unit
            tip: model.tip
            onEdited: (text) => Julia.br_edit_part(model.key, text)
        }
    }
    // a hohlraum part: each end an entrance hole or closed, and the closing holes
    Repeater {
        model: hohlraum ? [{ k: "top_end", t: "Top end", v: app.partTopEnd },
                           { k: "bottom_end", t: "Bottom end", v: app.partBottomEnd }] : []
        delegate: RowLayout {
            Label { text: modelData.t; Layout.preferredWidth: 110 }
            ComboBox {
                Layout.fillWidth: true
                model: ["Closed (the cap covers it)", "Entrance hole", "Open (to a package)"]
                currentIndex: modelData.v
                onActivated: (index) => Julia.br_edit_part(modelData.k, String(index))
                Tip {
                    visible: parent.hovered
                    text: "Entrance hole, closed cap, or fully open"
                }
            }
        }
    }
    Button {
        visible: hohlraum
        text: "Add window"
        onClicked: Julia.br_add_window()
        Tip {
            visible: parent.hovered
            text: "A closed 1 mm window in the wall, at the midplane"
        }
    }
    RowLayout {
        visible: hohlraum
        Layout.fillWidth: true
        CheckLine {
            text: "Closing entrance holes"
            checked: app.partLehClosure
            onToggled: Julia.br_edit_part("leh_closure", checked ? "1" : "0")
            tip: "Plasma closes the holes to x-rays during the pulse"
        }
        Info { topic: "part.leh_closure" }
    }
    CheckLine {
        visible: app.partKind === "Sphere"
        text: "Capsule mesh (geodesic, the capsule element size)"
        checked: app.partCapsuleMesh
        onToggled: Julia.br_edit_part("capsule_mesh", checked ? "1" : "0")
        tip: "Mesh the sphere as a capsule (geodesic)"
    }
    PlacementSteps { visible: app.selectionCount <= 1 }
    CheckLine {
        visible: !hohlraum
        text: "Two-sided (both sides radiate and absorb)"
        checked: app.partTwoSided
        enabled: !app.partAperture
        onToggled: Julia.br_part_flag("two_sided", checked)
    }
    CheckLine {
        visible: !hohlraum
        text: "Inward (radiates from its inside, or a flat part from its back)"
        checked: app.partInward
        onToggled: Julia.br_part_flag("inward", checked)
    }
    CheckLine {
        visible: !hohlraum
        text: "Aperture (transparent: rays pass, its flux is what escapes)"
        checked: app.partAperture
        onToggled: Julia.br_part_flag("aperture", checked)
    }
    DisplayRow { visible: !hohlraum }      // (a hohlraum's pieces have their own, in the tree)
    ErrorText { text: app.targetError }
    RowLayout {
        spacing: 8
        Button { text: app.selectionCount > 1 ? "Remove " + app.selectionCount + " selected" : "Remove part"
                 onClicked: Julia.br_remove_selected(); tip: "Undo brings it back" }
        Info { topic: "target.remove" }
    }
}
