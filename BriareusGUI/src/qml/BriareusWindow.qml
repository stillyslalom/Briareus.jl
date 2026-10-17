// The Briareus window, laid out like a CAD tool: menus and toolbar; the model tree
// (left), the 3D viewer over the result panels (center), the properties of the
// selection (right); a status bar. Nothing is a fixed sequence: any operation is
// available at any time, and analyses say when they are out of date.
//
// All state lives in Julia: `app` (property map) and the item models (`treeModel`,
// `targetFieldModel`, `partFieldModel`, `windowModel`, `beamModel`, `modeModel`,
// `powerModel`, `optModeModel`, `paramModel`, `compModel`, `clearModel`, `worstModel`,
// `talisCompModel`, `ruleModel`, `talisModel`, `noteModel`); actions call Julia.br_*.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Window
import QtQuick.Dialogs
import Makie
import jlqml

ApplicationWindow {
    id: win
    title: app.title
    width: 1640
    height: 1000
    minimumWidth: 1100
    minimumHeight: 680
    visible: true

    onClosing: Qt.quit()

    // The canvas takes the keyboard (and so the mouse) when the pointer enters it.
    function focusCanvas(hovered) {
        if (hovered && !canvas.activeFocus) canvas.forceActiveFocus(Qt.MouseFocusReason)
    }
    // Select an analysis: its settings in Properties, its results below.
    function showTask(key) {
        Julia.br_select("task:" + key)
        results.show(key)
    }
    function openKitDialog() { kitDialog.open() }
    function openClearanceExport() { fileDialog.run("export_clearance") }
    function openTalisExport() { fileDialog.run("export_talis") }

    // follow the selection with the result tabs
    readonly property string selectionWatch: app.selKind + ":" + app.selKey
    onSelectionWatchChanged: followSelection()
    function followSelection() {
        if (app.selKind === "task") results.show(app.selKey)
        else if (["beams", "group", "quad", "beam"].indexOf(app.selKind) >= 0) results.show("beams")
    }

    menuBar: MenuBar {
        Menu {
            title: "&File"
            Menu {
                title: "New target"
                Action { text: "Hohlraum"; onTriggered: { Julia.br_set_target_kind("hohlraum"); Julia.br_select("target") } }
                Action { text: "Backlighter foil"; onTriggered: { Julia.br_set_target_kind("foil"); Julia.br_select("target") } }
                Action { text: "Custom (parts only)"; onTriggered: { Julia.br_set_target_kind("custom"); Julia.br_select("target") } }
            }
            MenuSeparator {}
            Action { text: "Import beam table…"; onTriggered: fileDialog.run("import_beams") }
            Action { text: "Apply VisRad SST pointing…"; onTriggered: fileDialog.run("import_sst") }
            Action { text: "Apply VisRad beam summary…"; onTriggered: fileDialog.run("import_summary") }
            MenuSeparator {}
            Action { text: "Export beam table…"; onTriggered: fileDialog.run("export_beams") }
            Action { text: "Export VisRad SST pointing…"; onTriggered: fileDialog.run("export_sst") }
            Action { text: "Export VisRad beam summary…"; onTriggered: fileDialog.run("export_summary") }
            Action { text: "Write VisRad benchmark kit…"; enabled: app.kitReady && !app.kitRunning; onTriggered: kitDialog.open() }
            MenuSeparator {}
            Action { text: "Save view as image…"; onTriggered: saveViewDialog.open() }
            MenuSeparator {}
            Action { text: "&Quit"; shortcut: StandardKey.Quit; onTriggered: win.close() }
        }
        Menu {
            title: "&Target"
            Action { text: "Properties"; onTriggered: Julia.br_select("target") }
            Menu {
                title: "Add part"
                Action { text: "Sphere"; onTriggered: Julia.br_add_part("sphere") }
                Action { text: "Disk"; onTriggered: Julia.br_add_part("disk") }
                Action { text: "Washer (annulus)"; onTriggered: Julia.br_add_part("annulus") }
                Action { text: "Rectangle"; onTriggered: Julia.br_add_part("rectangle") }
                Action { text: "Cylinder (tube)"; onTriggered: Julia.br_add_part("cylinder") }
                Action { text: "Cone (frustum)"; onTriggered: Julia.br_add_part("cone") }
            }
            Action { text: "Add diagnostic window"; enabled: app.targetKind === "hohlraum"; onTriggered: Julia.br_add_window() }
            Action {
                text: "Remove selected part or window"
                enabled: app.selKind === "part" || app.selKind === "window"
                onTriggered: Julia.br_remove_selected()
            }
            MenuSeparator {}
            Menu {
                title: "Mesh"
                Action { text: "Coarse"; checkable: true; checked: app.targetResolution === "coarse"; onTriggered: Julia.br_set_resolution("coarse") }
                Action { text: "Standard"; checkable: true; checked: app.targetResolution === "standard"; onTriggered: Julia.br_set_resolution("standard") }
                Action { text: "Fine (converged)"; checkable: true; checked: app.targetResolution === "fine"; onTriggered: Julia.br_set_resolution("fine") }
            }
        }
        Menu {
            title: "&Beams"
            Action { text: "Properties"; onTriggered: Julia.br_select("beams") }
            Menu {
                title: "Load preset"
                Action { text: "NIF quads (48)"; onTriggered: Julia.br_load_preset("nif_quads") }
                Action { text: "NIF beams (192)"; onTriggered: Julia.br_load_preset("nif_beams") }
                Action { text: "OMEGA hohlraum (40)"; onTriggered: Julia.br_load_preset("omega_hohlraum") }
                Action { text: "OMEGA backlighter"; onTriggered: Julia.br_load_preset("omega_backlighter") }
            }
            Action { text: "All on"; onTriggered: Julia.br_tree_on("beams", true) }
            Action { text: "All off"; onTriggered: Julia.br_tree_on("beams", false) }
        }
        Menu {
            title: "&View"
            Action { text: "Fit all"; onTriggered: Julia.br_fit_all() }
            Action { text: "Fit selection"; onTriggered: Julia.br_fit_selection() }
            Action {
                text: "Look along selected beam"
                enabled: ["group", "quad", "beam"].indexOf(app.selKind) >= 0
                onTriggered: Julia.br_look_along_beam()
            }
            Menu {
                title: "Standard view"
                Action { text: "Isometric"; onTriggered: Julia.br_standard_view("iso") }
                Action { text: "Front (along +y)"; onTriggered: Julia.br_standard_view("front") }
                Action { text: "Back (along −y)"; onTriggered: Julia.br_standard_view("back") }
                Action { text: "Right (along −x)"; onTriggered: Julia.br_standard_view("right") }
                Action { text: "Left (along +x)"; onTriggered: Julia.br_standard_view("left") }
                Action { text: "Top (along −z)"; onTriggered: Julia.br_standard_view("top") }
                Action { text: "Bottom (along +z)"; onTriggered: Julia.br_standard_view("bottom") }
            }
            MenuSeparator {}
            Action { text: "Orthographic projection"; checkable: true; checked: app.ortho; onTriggered: Julia.br_projection(checked) }
            Action { text: "Turntable rotation (z stays up)"; checkable: true; checked: app.turntable; onTriggered: Julia.br_turntable(checked) }
            Action { text: "Cutaway"; checkable: true; checked: app.cutaway; onTriggered: Julia.br_set_cutaway(checked) }
            MenuSeparator {}
            Menu {
                title: "Beams"
                Action { text: "Rays"; checkable: true; checked: app.beamStyle === "rays"; onTriggered: Julia.br_beam_style("rays") }
                Action { text: "Volumes"; checkable: true; checked: app.beamStyle === "volumes"; onTriggered: Julia.br_beam_style("volumes") }
                Action { text: "Hidden"; checkable: true; checked: app.beamStyle === "off"; onTriggered: Julia.br_beam_style("off") }
                MenuSeparator {}
                Action { text: "Color by cone"; checkable: true; checked: app.beamColor === "group"; onTriggered: Julia.br_beam_color("group") }
                Action { text: "Color by clearance"; checkable: true; checked: app.beamColor === "clearance"; onTriggered: Julia.br_beam_color("clearance") }
                Action { text: "Color by TaLIS"; checkable: true; checked: app.beamColor === "talis"; onTriggered: Julia.br_beam_color("talis") }
            }
            Menu {
                title: "Color target by"
                Action { text: "Components"; checkable: true; checked: app.colorBy === "component"; onTriggered: Julia.br_set_color("component") }
                Action { text: "Absorbed laser"; checkable: true; checked: app.colorBy === "laser"; onTriggered: Julia.br_set_color("laser") }
                Action { text: "X-ray radiosity"; checkable: true; checked: app.colorBy === "radiosity"; onTriggered: Julia.br_set_color("radiosity") }
                Action { text: "Incident x-ray flux"; checkable: true; checked: app.colorBy === "incident"; onTriggered: Julia.br_set_color("incident") }
            }
            Action {
                text: app.hasFoil ? "Foil map" : "Capsule θ–φ map"
                checkable: true
                checked: app.driveView === "map"
                enabled: app.driveHasResult && (app.hasCapsule || app.hasFoil)
                onTriggered: Julia.br_set_view(checked ? "map" : "scene")
            }
        }
        Menu {
            title: "&Compute"
            Action { text: "Drive (view factors)"; shortcut: "Ctrl+R"; enabled: !app.driveBuilding; onTriggered: { Julia.br_build(); results.show("drive") } }
            Action { text: "Optimize pointing…"; onTriggered: win.showTask("optimize") }
            Action { text: "Check clearance…"; onTriggered: win.showTask("clearance") }
            Action { text: "TaLIS rules…"; onTriggered: win.showTask("talis") }
            MenuSeparator {}
            Action { text: "Cancel running analyses"; onTriggered: Julia.br_cancel_all() }
        }
        Menu {
            title: "&Help"
            Action { text: "Mouse and keys"; onTriggered: helpDialog.open() }
        }
    }

    header: MainToolbar {}

    footer: ToolBar {
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: 10
            anchors.rightMargin: 10
            spacing: 10
            BusyIndicator { running: app.busy !== ""; visible: running; Layout.preferredHeight: 20; Layout.preferredWidth: 20 }
            Label { text: app.busy; visible: text !== ""; color: "#2f7fd8" }
            Label { text: app.status; elide: Text.ElideRight; Layout.fillWidth: true; opacity: 0.8 }
            Label {
                text: "wheel zoom · drag orbit · right-drag pan · double-click fit · F fit all"
                opacity: 0.55
                font.pixelSize: 11
            }
        }
    }

    SplitView {
        id: body
        anchors.fill: parent
        orientation: Qt.Horizontal

        ModelTree {
            SplitView.preferredWidth: 290
            SplitView.minimumWidth: 200
        }

        SplitView {
            SplitView.fillWidth: true
            SplitView.minimumWidth: 420
            orientation: Qt.Vertical

            // The canvas item lives as long as the window: QMLMakie canvases must not
            // be destroyed while the application runs.
            Item {
                SplitView.fillHeight: true
                SplitView.minimumHeight: 260
                MakieArea {
                    id: canvas
                    anchors.fill: parent
                    scene: app.canvas
                    HoverHandler { onHoveredChanged: win.focusCanvas(hovered) }
                }
                Label {
                    anchors.top: parent.top
                    anchors.left: parent.left
                    anchors.margins: 8
                    text: app.driveView === "map" ? "" :
                          app.colorBy !== "component" && !app.driveCurrent
                          ? (app.driveBuilding ? "computing the drive… " + app.driveElapsed : "drive out of date: Compute ▸ Drive")
                          : ""
                    visible: text !== ""
                    color: "#d08a00"
                    background: Rectangle { color: palette.window; opacity: 0.85; radius: 3 }
                    padding: 4
                }
            }
            ResultsPanel {
                id: results
                SplitView.preferredHeight: 300
                SplitView.minimumHeight: 120
            }
        }

        Inspector {
            SplitView.preferredWidth: 400
            SplitView.minimumWidth: 300
        }
    }

    // one dialog for the beam files and report exports
    FileDialog {
        id: fileDialog
        property string mode: "import_beams"
        function run(m) { mode = m; open() }
        readonly property bool saving: mode.startsWith("export")
        title: {
            if (mode === "import_beams") return "Import a beam table"
            if (mode === "import_sst") return "Apply a NIF SST pointing file"
            if (mode === "import_summary") return "Apply a VisRad Laser Beam Summary"
            if (mode === "export_beams") return "Export the beam table"
            if (mode === "export_sst") return "Write the pointing as a NIF SST file"
            if (mode === "export_summary") return "Write a VisRad Laser Beam Summary"
            if (mode === "export_clearance") return "Export the clearance report"
            return "Export the TaLIS report"
        }
        fileMode: saving ? FileDialog.SaveFile : FileDialog.OpenFile
        defaultSuffix: "csv"
        nameFilters: ["CSV (*.csv *.tsv *.txt)", "All files (*)"]
        onAccepted: {
            const url = selectedFile.toString()
            if (mode === "import_beams") Julia.br_import_beams(url)
            else if (mode === "import_sst") Julia.br_import_sst(url)
            else if (mode === "import_summary") Julia.br_import_summary(url)
            else if (mode === "export_beams") Julia.br_export_beams(url)
            else if (mode === "export_sst") Julia.br_export_sst(url)
            else if (mode === "export_summary") Julia.br_export_summary(url)
            else if (mode === "export_clearance") Julia.br_rv_export(url)
            else Julia.br_tl_export(url)
        }
    }
    FolderDialog {
        id: kitDialog
        title: "Folder for the VisRad benchmark kit"
        onAccepted: Julia.br_write_kit(selectedFolder.toString())
    }
    FileDialog {
        id: saveViewDialog
        title: "Save the viewer as an image"
        fileMode: FileDialog.SaveFile
        defaultSuffix: "png"
        nameFilters: ["PNG image (*.png)"]
        onAccepted: {
            const target = selectedFile
            canvas.grabToImage(function (img) { img.saveToFile(target) })
        }
    }
    Dialog {
        id: helpDialog
        title: "Mouse and keys"
        anchors.centerIn: parent
        standardButtons: Dialog.Close
        Label {
            text: "In the 3D view\n\n" +
                  "Wheel: zoom toward the point under the pointer\n" +
                  "Left or middle drag: orbit about the selection, or the point under the pointer\n" +
                  "Right drag, Shift + drag, or Ctrl + middle drag: pan\n" +
                  "Click: select a component or beam (empty space clears); Ctrl + click adds beams\n" +
                  "Double click: fit the view to what was clicked (empty space: everything)\n" +
                  "F: fit everything · arrow keys: orbit 15°\n\n" +
                  "The model tree selects the same things; its boxes switch beams on and off\n" +
                  "and show or hide components. View ▸ Turntable keeps z up while orbiting."
        }
    }

    Timer {
        interval: 16
        running: true
        repeat: true
        onTriggered: {
            const r = Julia.br_tick()
            if (app.grabPath !== "") {
                const path = app.grabPath
                Julia.br_grab_started()
                body.grabToImage(function (img) { img.saveToFile(path) })
            }
            if (r === 2) {
                win.close()
            } else if (r === 1) {
                canvas.update()
            }
        }
    }
}
