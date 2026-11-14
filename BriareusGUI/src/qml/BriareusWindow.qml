// The Briareus window, laid out like a CAD tool: menus and toolbar; the model tree
// (left), the 3D viewer over the result panels (center), the properties of the
// selection (right); a status bar. Nothing is a fixed sequence: any operation is
// available at any time, and analyses say when they are out of date.
//
// All state lives in Julia: `app` (property map) and the item models (`treeModel`,
// `partFieldModel`, `windowModel`, `beamModel`, `modeModel`,
// `powerModel`, `optModeModel`, `paramModel`, `compModel`, `clearModel`, `worstModel`,
// `talisCompModel`, `ruleModel`, `talisModel`, `noteModel`, `timeLegendModel`); actions
// call Julia.br_*.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Window
import QtQuick.Dialogs
import Makie
import jlqml
import "common"

ApplicationWindow {
    id: win
    title: app.title
    // 1640 × 1000, centered on the screen; on a screen whose free area (less the taskbar)
    // is smaller than that with a margin (a laptop, a scaled display), maximized. Set once
    // at start: the bindings are cut in onCompleted, so moving the window to another
    // screen keeps its size.
    readonly property int availW: Math.min(Screen.width, Screen.desktopAvailableWidth)
    readonly property int availH: Math.min(Screen.height, Screen.desktopAvailableHeight)
    width: Math.min(1640, Math.round(0.9 * availW))
    height: Math.min(1000, Math.round(0.9 * availH))
    x: Screen.virtualX + Math.max(0, Math.round((availW - width) / 2))
    y: Screen.virtualY + Math.max(0, Math.round((availH - height) / 2))
    minimumWidth: Math.min(1100, availW)
    minimumHeight: Math.min(680, availH)
    // a denser text than the style's 14 px, as in an Office ribbon: the controls in "common"
    // follow it with their own tighter padding
    font.pixelSize: 13
    visible: true
    readonly property bool smallScreen: availW < 1800 || availH < 1100
    Component.onCompleted: {
        width = width; height = height; x = x; y = y
        if (smallScreen) showMaximized()
    }

    onClosing: Qt.quit()
    // keys released while another window had focus never reach the canvas
    onActiveChanged: if (!active) Julia.br_canvas_focus(false)

    // The canvas takes the keyboard (and so the mouse) when the pointer enters it.
    function focusCanvas(hovered) {
        if (hovered && !canvas.activeFocus) canvas.forceActiveFocus(Qt.MouseFocusReason)
    }
    // Select an analysis: its settings in Properties, its results below.
    function showTask(key) {
        Julia.br_select("task:" + key)
        results.show(key)
    }
    function openClearanceExport() { fileDialog.run("export_clearance") }
    function openTalisExport() { fileDialog.run("export_talis") }
    function openDiagnostics() { diagnosticsDialog.open() }
    function openStl() { stlDialog.open() }

    // a results tab asked for from Julia (`request_tab`)
    readonly property string resultsRequest: app.resultsRequest
    onResultsRequestChanged: if (resultsRequest !== "") results.show(resultsRequest)
    // a window size asked for from Julia (`request_window_size`, "width x height")
    readonly property string sizeRequest: app.sizeRequest
    onSizeRequestChanged: if (sizeRequest !== "") {
        const wh = sizeRequest.split("x")
        width = parseInt(wh[0])
        height = parseInt(wh[1])
    }
    // a dialog or menu asked for from Julia (`request_popup`), and an image of it when
    // `app.popupGrabPath` is set (scripted screenshots)
    property var openPopup: null
    readonly property string popupRequest: app.popupRequest
    onPopupRequestChanged: {
        if (openPopup) openPopup.close()
        openPopup = null
        const key = popupRequest
        if (key === "") return
        // "window": no popup (an image of the toolbar)
        const dialogs = { retune: retuneDialog, port: portDialog, transform: transformDialog }
        if (key.startsWith("topics")) { helpTopics.openAt(key.substring(7)); openPopup = helpTopics }
        else if (key === "help") { helpTopics.openAt("start.keys"); openPopup = helpTopics }
        else if (key === "group") { tree.openGroupDialog(); openPopup = tree.groupDialog }
        else if (key.startsWith("rowmenu:")) {
            // the model tree's row menu for a node id (a beam row, a frame, an analysis)
            const id = key.substring(8)
            tree.openMenu(id, id.split(":")[0], "")
            openPopup = tree.rowMenu
        }
        else if (key === "diagnostics") { results.show("diagnostics"); openPopup = results.openDiagnosticPicker() }
        else if (key === "images") { results.show("diagnostics"); Julia.br_diag_category("xray_imager") }
        else if (dialogs[key]) { openPopup = dialogs[key]; openPopup.open() }
        else if (key !== "window") {
            for (let i = 0; i < menuBar.count; ++i) {
                const m = menuBar.menuAt(i)
                if (m && m.title.replace("&", "") === key) {
                    const it = menuBar.itemAt(i)
                    m.popup(it, 0, it.height)
                    openPopup = m
                }
            }
        }
        if (app.popupGrabPath !== "") popupGrabTimer.start()
    }
    Timer {
        id: popupGrabTimer
        interval: 500
        onTriggered: {
            const path = app.popupGrabPath
            // the popup alone (the window's root item cannot be grabbed), the Help window's
            // content, or the toolbar
            const item = openPopup === helpTopics ? helpTopics.contentItem
                       : openPopup ? openPopup.contentItem.parent : win.header
            item.grabToImage(function (img) { img.saveToFile(path) })
        }
    }
    // the Help screen at a topic (an ⓘ clicked, `br_help`)
    readonly property string helpRequest: app.helpRequest
    onHelpRequestChanged: if (helpRequest.startsWith("topic:")) helpTopics.openAt(helpRequest.substring(6))
    HelpWindow { id: helpTopics }
    // follow the selection with the result tabs
    readonly property string selectionWatch: app.selKind + ":" + app.selKey
    onSelectionWatchChanged: followSelection()
    function followSelection() {
        if (app.selKind === "task") results.show(app.selKey)
        else if (["beams", "group", "quad", "beam"].indexOf(app.selKind) >= 0) results.show("beams")
    }

    // as dense as the toolbar: the style gives the bar a toolbar's height and each title
    // a button's
    menuBar: MenuBar {
        topPadding: 0
        bottomPadding: 0
        implicitHeight: contentItem.implicitHeight
        delegate: MenuBarItem {
            topPadding: 3
            bottomPadding: 3
            leftPadding: 8
            rightPadding: 8
            implicitHeight: 24
        }
        AutoMenu {
            title: "&File"
            Action { text: "&Open session…"; shortcut: StandardKey.Open; onTriggered: sessionDialog.run("open") }
            Action {
                text: "&Save session"
                shortcut: StandardKey.Save
                onTriggered: app.sessionPath === "" ? sessionDialog.run("save") : Julia.br_save_session("", false)
            }
            Action { text: "Save session &as…"; shortcut: StandardKey.SaveAs; onTriggered: sessionDialog.run("save") }
            Action { text: "Save session without CAD files…"; onTriggered: sessionDialog.run("save_refs") }
            MenuSeparator {}
            // ready-made targets, by the facility they are mock-ups for: each comes with
            // that facility's beams
            AutoMenu {
                title: "New target"
                AutoMenu {
                    title: "NIF"
                    Action { text: "Hohlraum with capsule"; onTriggered: Julia.br_new_target("hohlraum") }
                    Action { text: "Frustraum with capsule"; onTriggered: Julia.br_new_target("frustraum") }
                    Action { text: "Halfraum pair with shock tube"; onTriggered: Julia.br_new_target("shock_tube") }
                    Action { text: "Opacity halfraum"; onTriggered: Julia.br_new_target("opacity") }
                }
                AutoMenu {
                    title: "OMEGA"
                    Action { text: "Backlighter foil"; onTriggered: Julia.br_new_target("foil") }
                }
                MenuSeparator {}
                Action { text: "Empty"; onTriggered: Julia.br_new_target("empty") }
            }
            MenuSeparator {}
            Action { text: "Import VisRad workspace…"; onTriggered: vrwDialog.run("import") }
            Action {
                text: app.visradSource === "" ? "Export VisRad workspace… (import one first)"
                                              : "Export VisRad workspace (beams into " + app.visradSource + ")…"
                enabled: app.visradSource !== ""
                onTriggered: vrwDialog.run("export")
            }
            MenuSeparator {}
            Action { text: "Import beam table…"; onTriggered: fileDialog.run("import_beams") }
            Action { text: "Apply VisRad SST pointing…"; onTriggered: fileDialog.run("import_sst") }
            Action { text: "Apply VisRad beam summary…"; onTriggered: fileDialog.run("import_summary") }
            MenuSeparator {}
            Action { text: "Export beam table…"; onTriggered: fileDialog.run("export_beams") }
            Action { text: "Export VisRad SST pointing…"; onTriggered: fileDialog.run("export_sst") }
            Action { text: "Export VisRad beam summary…"; onTriggered: fileDialog.run("export_summary") }
            Action {
                // the drive, and the drive history with its diagnostics when it is current
                text: app.timeHasHistory && !app.timeStale ? "Export results with the drive history (CSV)…"
                                                           : "Export results (CSV)…"
                enabled: app.driveCurrent
                onTriggered: resultsDialog.open()
            }
            Action { text: "Export target as STL…"; onTriggered: exportStlDialog.open() }
            MenuSeparator {}
            Action { text: "Save view as image…"; onTriggered: saveViewDialog.open() }
            MenuSeparator {}
            Action { text: "&Quit"; shortcut: StandardKey.Quit; onTriggered: win.close() }
        }
        AutoMenu {
            title: "&Edit"
            Action { text: "&Undo"; shortcut: StandardKey.Undo; enabled: app.canUndo; onTriggered: Julia.br_undo() }
            Action { text: "&Redo"; shortcut: StandardKey.Redo; enabled: app.canRedo; onTriggered: Julia.br_redo() }
            MenuSeparator {}
            Action {
                text: app.selectionCount > 1 ? "Remove " + app.selectionCount + " selected" : "Remove selected"
                shortcut: StandardKey.Delete
                enabled: ["part", "assembly", "window", "frame"].indexOf(app.selKind) >= 0
                onTriggered: Julia.br_remove_selected()
            }
            Action { text: "Deselect"; enabled: app.selKind !== ""; onTriggered: Julia.br_select("") }
        }
        AutoMenu {
            // what a target is built from: parts, windows, assemblies and frames
            title: "&Insert"
            Action { text: "Sphere"; onTriggered: Julia.br_add_part("sphere") }
            Action { text: "Disk"; onTriggered: Julia.br_add_part("disk") }
            Action { text: "Washer (annulus)"; onTriggered: Julia.br_add_part("annulus") }
            Action { text: "Rectangle"; onTriggered: Julia.br_add_part("rectangle") }
            Action { text: "Cylinder (tube)"; onTriggered: Julia.br_add_part("cylinder") }
            Action { text: "Cone (frustum)"; onTriggered: Julia.br_add_part("cone") }
            Action { text: "Hohlraum"; onTriggered: Julia.br_add_part("hohlraum") }
            Action { text: "Mesh from STL file…"; onTriggered: stlDialog.open() }
            MenuSeparator {}
            Action { text: "Diagnostic window"; enabled: app.hasHohlraum; onTriggered: Julia.br_add_window() }
            MenuSeparator {}
            Action { text: "Assembly (of the selected parts)"; onTriggered: Julia.br_new_assembly() }
            Action { text: "Frame (in the selected frame)"; onTriggered: Julia.br_add_frame() }
        }
        AutoMenu {
            title: "&Target"
            Action { text: "Properties"; onTriggered: Julia.br_select("target") }
            AutoMenu {
                title: "Mesh"
                Action { text: "Coarse"; checkable: true; checked: app.targetResolution === "coarse"; onTriggered: Julia.br_set_resolution("coarse") }
                Action { text: "Standard"; checkable: true; checked: app.targetResolution === "standard"; onTriggered: Julia.br_set_resolution("standard") }
                Action { text: "Fine (converged)"; checkable: true; checked: app.targetResolution === "fine"; onTriggered: Julia.br_set_resolution("fine") }
            }
        }
        AutoMenu {
            title: "&Beams"
            Action { text: "Properties"; onTriggered: Julia.br_select("beams") }
            AutoMenu {
                title: "Facility"
                Action { text: "NIF"; checkable: true; checked: app.facility === "nif"; onTriggered: Julia.br_set_facility("nif") }
                Action { text: "OMEGA"; checkable: true; checked: app.facility === "omega"; onTriggered: Julia.br_set_facility("omega") }
                Action { text: "OMEGA EP"; checkable: true; checked: app.facility === "omega_ep"; onTriggered: Julia.br_set_facility("omega_ep") }
                Action { text: "Other"; checkable: true; checked: app.facility === "custom"; onTriggered: Julia.br_set_facility("custom") }
            }
            AutoMenu {
                title: "Load preset"
                enabled: app.facility === "nif" || app.facility === "omega"
                Action { text: "NIF quads (48)"; enabled: app.facility === "nif"; onTriggered: Julia.br_load_preset("nif_quads") }
                Action { text: "NIF beams (192)"; enabled: app.facility === "nif"; onTriggered: Julia.br_load_preset("nif_beams") }
                Action { text: "OMEGA hohlraum (40)"; enabled: app.facility === "omega"; onTriggered: Julia.br_load_preset("omega_hohlraum") }
                Action { text: "OMEGA backlighter"; enabled: app.facility === "omega"; onTriggered: Julia.br_load_preset("omega_backlighter") }
            }
            MenuSeparator {}
            Action { text: "Split all quads into beams"; enabled: app.hasQuads; onTriggered: Julia.br_split_quads("beams") }
            Action { text: "Merge all beams into quads"; enabled: app.hasQuads; onTriggered: Julia.br_merge_quads("beams") }
            MenuSeparator {}
            Action { text: "All on"; onTriggered: Julia.br_tree_on("beams", true) }
            Action { text: "All off"; onTriggered: Julia.br_tree_on("beams", false) }
        }
        AutoMenu {
            title: "&View"
            Action { text: "Fit all"; onTriggered: Julia.br_fit_all() }
            Action { text: "Fit selection"; onTriggered: Julia.br_fit_selection() }
            Action { text: "Look from port…"; onTriggered: portDialog.open() }
            Action {
                text: "Look along selected beam"
                enabled: ["group", "quad", "beam"].indexOf(app.selKind) >= 0
                onTriggered: Julia.br_look_along_beam()
            }
            AutoMenu {
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
            Action { text: "Turntable rotation (keeps z up)"; checkable: true; checked: app.turntable; onTriggered: Julia.br_turntable(checked) }
            Action { text: "Cutaway"; checkable: true; checked: app.cutaway; onTriggered: Julia.br_set_cutaway(checked) }
            MenuSeparator {}
            AutoMenu {
                title: "Beams"
                Action { text: "Rays"; checkable: true; checked: app.beamStyle === "rays"; onTriggered: Julia.br_beam_style("rays") }
                Action { text: "Volumes"; checkable: true; checked: app.beamStyle === "volumes"; onTriggered: Julia.br_beam_style("volumes") }
                Action { text: "Hidden"; checkable: true; checked: app.beamStyle === "off"; onTriggered: Julia.br_beam_style("off") }
                MenuSeparator {}
                Action { text: "Show beams that are off"; checkable: true; checked: app.showOff; onTriggered: Julia.br_show_off(checked) }
                MenuSeparator {}
                Action { text: "Color by cone"; checkable: true; checked: app.beamColor === "group"; onTriggered: Julia.br_beam_color("group") }
                Action { text: "Color by clearance"; checkable: true; checked: app.beamColor === "clearance"; onTriggered: Julia.br_beam_color("clearance") }
                Action { text: "Color by TaLIS"; checkable: true; checked: app.beamColor === "talis"; onTriggered: Julia.br_beam_color("talis") }
            }
            AutoMenu {
                title: "Color target by"
                Action { text: "Components"; checkable: true; checked: app.colorBy === "component"; onTriggered: Julia.br_set_color("component") }
                Action { text: "Absorbed laser"; checkable: true; checked: app.colorBy === "laser"; onTriggered: Julia.br_set_color("laser") }
                Action { text: "X-ray radiosity"; checkable: true; checked: app.colorBy === "radiosity"; onTriggered: Julia.br_set_color("radiosity") }
                Action { text: "Incident x-ray flux"; checkable: true; checked: app.colorBy === "incident"; onTriggered: Julia.br_set_color("incident") }
            }
            MenuSeparator {}
            Action {
                text: "3D view"
                checkable: true
                checked: app.driveView === "scene"
                onTriggered: Julia.br_set_view("scene")
            }
            Action {
                text: app.observedKind === "face" ? "Observed face map" :
                      app.observedKind === "capsule" ? "Capsule θ–φ map" : "Foil map"
                checkable: true
                checked: app.driveView === "map"
                enabled: app.driveHasResult && (app.hasCapsule || app.hasFoil)
                onTriggered: Julia.br_set_view(checked ? "map" : "scene")
            }
            Action {
                text: "Map or image beside the 3D view"
                checkable: true
                checked: app.sideView
                onTriggered: Julia.br_side_view(checked)
            }
        }
        AutoMenu {
            title: "&Compute"
            Action { text: "Drive"; shortcut: "Ctrl+R"; enabled: !app.driveBuilding; onTriggered: { Julia.br_build(); results.show("drive") } }
            Action { text: "Optimization…"; onTriggered: win.showTask("optimize") }
            Action { text: "Clearance…"; onTriggered: win.showTask("clearance") }
            Action { text: "Placement sweep…"; onTriggered: win.showTask("sweep") }
            Action {
                text: "X-ray image along the view"
                checkable: true
                checked: app.driveView === "image"
                // unchecked: back to the view shown before
                onTriggered: checked ? Julia.br_xray_image() : Julia.br_close_image()
            }
            Action { text: "TaLIS (NIF)…"; enabled: app.facility === "nif"; onTriggered: win.showTask("talis") }
            Action { text: "Retune to a reference drive…"; enabled: app.driveCurrent; onTriggered: retuneDialog.open() }
            Action {
                text: "Drive history"
                enabled: app.timeHasPulses && !app.timeRunning
                onTriggered: { Julia.br_time_run(); results.show("time") }
            }
            MenuSeparator {}
            Action { text: "Cancel running analyses"; onTriggered: Julia.br_cancel_all() }
        }
        AutoMenu {
            title: "&Tools"
            Action { text: "Coordinate transform…"; onTriggered: transformDialog.open() }
            Action { text: "Frames…"; onTriggered: Julia.br_select("frames") }
        }
        AutoMenu {
            title: "&Diagnostics"
            Action { text: "Show the Diagnostics tab"; onTriggered: results.show("diagnostics") }
            MenuSeparator {}
            Action { text: "Diagnostic directories…"; onTriggered: diagnosticsDialog.open() }
            Action { text: "Refresh"; onTriggered: Julia.br_diag_reload() }
        }
        AutoMenu {
            title: "&Help"
            Action { text: "Help"; shortcut: StandardKey.HelpContents; onTriggered: helpTopics.openAt("start.layout") }
            Action { text: "Mouse and keys"; onTriggered: helpTopics.openAt("start.keys") }
            Action { text: "Documentation"; onTriggered: Qt.openUrlExternally(Julia.br_docs_url()) }
        }
    }

    header: MainToolbar {}

    footer: ToolBar {
        topPadding: 0
        bottomPadding: 0
        contentHeight: 22
        RowLayout {
            anchors.fill: parent
            anchors.leftMargin: 10
            anchors.rightMargin: 10
            spacing: 10
            BusyIndicator { running: app.busy !== ""; visible: running; Layout.preferredHeight: 18; Layout.preferredWidth: 18 }
            Label { text: app.busy; visible: text !== ""; color: "#2f7fd8" }
            // the latest note, cut short to the bar: the whole of it on hover
            Label {
                id: statusLine
                text: app.status
                elide: Text.ElideRight
                Layout.fillWidth: true
                opacity: 0.8
                HoverHandler { id: statusHover }
                Tip { visible: statusLine.truncated && statusHover.hovered; delay: 300; text: statusLine.text; maxWidth: 640 }
            }
            Label {
                // the gestures of the view on the canvas: the 3D view's, or a map's or image's
                text: app.driveView === "scene" || app.sideView
                      ? "wheel zoom · drag orbit · right-drag pan · double-click fit · F fit all · Esc deselect"
                      : "wheel zoom · drag box zoom · right-drag pan · Ctrl+click opening view · Ctrl+Shift+click whole field"
                opacity: 0.55
                font.pixelSize: 11
            }
            Info { topic: "view.image"; visible: !(app.driveView === "scene" || app.sideView); font.pixelSize: 11 }
        }
    }

    SplitView {
        id: body
        anchors.fill: parent
        orientation: Qt.Horizontal
        // Esc clears the selection. The key reaches here from whatever had the keyboard
        // (a field, the tree) unless that used it; a popup closes on Esc before this, and
        // the 3D view, which takes its own keys, clears the selection itself.
        focus: true
        Keys.onEscapePressed: (event) => { Julia.br_select(""); event.accepted = true }

        ModelTree {
            id: tree
            // side panels and the results take a share of a small window, the 3D view the rest
            SplitView.preferredWidth: Math.min(290, Math.round(0.17 * win.availW))
            SplitView.minimumWidth: 170
        }

        SplitView {
            id: centerSplit
            SplitView.fillWidth: true
            SplitView.minimumWidth: 420
            orientation: Qt.Vertical

            // The canvas items live as long as the window: QMLMakie canvases must not
            // be destroyed while the application runs. The side canvas holds the map or
            // the x-ray image beside the 3D view when that is chosen (View menu).
            SplitView {
                SplitView.fillHeight: true
                SplitView.minimumHeight: 260
                orientation: Qt.Horizontal
                Item {
                    SplitView.fillWidth: true
                    SplitView.minimumWidth: 260
                    MakieArea {
                        id: canvas
                        anchors.fill: parent
                        scene: app.canvas
                        HoverHandler { onHoveredChanged: win.focusCanvas(hovered) }
                        onActiveFocusChanged: Julia.br_canvas_focus(activeFocus)
                    }
                    Label {
                        anchors.top: parent.top
                        anchors.left: parent.left
                        anchors.margins: 8
                        text: app.driveView !== "scene" && !app.sideView ? "" :
                              app.colorBy !== "component" && !app.driveCurrent
                              ? (app.driveBuilding ? "computing the drive… " + app.driveElapsed : "drive out of date: Compute ▸ Drive")
                              : ""
                        visible: text !== ""
                        color: "#d08a00"
                        background: Rectangle { color: palette.window; opacity: 0.85; radius: 3 }
                        padding: 4
                    }
                    // the drive history's step the 3D view shows (Time tab)
                    Label {
                        anchors.top: parent.top
                        anchors.horizontalCenter: parent.horizontalCenter
                        anchors.margins: 8
                        text: app.timeLabel !== "" && app.colorBy !== "component" &&
                              (app.driveView === "scene" || app.sideView) ? "drive history at " + app.timeLabel : ""
                        visible: text !== ""
                        font.weight: Font.DemiBold
                        color: "#2f7fd8"
                        background: Rectangle { color: palette.window; opacity: 0.85; radius: 3 }
                        padding: 4
                    }
                }
                Item {
                    visible: app.sideView
                    SplitView.preferredWidth: 420
                    SplitView.minimumWidth: app.sideView ? 220 : 0
                    MakieArea {
                        id: sideCanvas
                        anchors.fill: parent
                        scene: app.sideCanvas
                    }
                }
            }
            ResultsPanel {
                id: results
                SplitView.preferredHeight: Math.min(300, Math.round(0.27 * win.availH))
                SplitView.minimumHeight: 100
                // a tab that needs more height (Time's chart) grows the panel when it opens,
                // the 3D view keeping its minimum; the next tab returns the panel to its
                // height unless the splitter moved meanwhile (a shorter Time tab scrolls)
                property real heightBefore: -1
                property real grownTo: -1
                function fitTab() {
                    const need = Math.min(fitHeight, centerSplit.height - 260 - 8)
                    if (height < need - 1) {
                        if (heightBefore < 0) heightBefore = height
                        results.SplitView.preferredHeight = need
                        grownTo = need
                    } else if (fitHeight === 0 && heightBefore >= 0) {
                        if (Math.abs(height - grownTo) < 1) results.SplitView.preferredHeight = heightBefore
                        heightBefore = -1
                    }
                }
                onCurrentIndexChanged: Qt.callLater(fitTab)
            }
        }

        Inspector {
            SplitView.preferredWidth: Math.min(400, Math.round(0.24 * win.availW))
            SplitView.minimumWidth: 260
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
    // session files: the whole design (target, beams, settings) and the CAD files it reads
    // (.briar, a zip), or the design alone as TOML with the CAD files by path
    FileDialog {
        id: sessionDialog
        property string mode: "open"
        function run(m) { mode = m; open() }
        title: mode === "open" ? "Open a Briareus session" :
               mode === "save" ? "Save the Briareus session with its CAD files" :
               "Save the Briareus session without its CAD files (by path)"
        fileMode: mode === "open" ? FileDialog.OpenFile : FileDialog.SaveFile
        defaultSuffix: mode === "save_refs" ? "toml" : "briar"
        nameFilters: mode === "open" ? ["Briareus session (*.briar *.briar.toml *.toml)", "All files (*)"] :
                     mode === "save" ? ["Briareus session with its CAD files (*.briar)"] :
                     ["Briareus session, CAD files by path (*.briar.toml)"]
        onAccepted: {
            const url = selectedFile.toString()
            if (mode === "open") Julia.br_open_session(url)
            else Julia.br_save_session(url, mode === "save_refs")
        }
    }
    // after opening a session: CAD files changed on disk since it was saved
    Dialog {
        id: cadDialog
        title: "CAD files changed"
        anchors.centerIn: parent
        modal: true
        closePolicy: Popup.NoAutoClose
        footer: DialogButtonBox {
            Button { text: "Read the changed files"; DialogButtonBox.buttonRole: DialogButtonBox.AcceptRole }
            Button { text: "Keep the saved copies"; DialogButtonBox.buttonRole: DialogButtonBox.RejectRole }
        }
        Label { text: app.cadPrompt }
        onAccepted: Julia.br_reload_cad_sources()
        onRejected: Julia.br_keep_cad_copies()
    }
    FileDialog {
        id: stlDialog
        title: "Add a part from an STL mesh (CAD export)"
        fileMode: FileDialog.OpenFile
        nameFilters: ["STL mesh (*.stl)", "All files (*)"]
        onAccepted: Julia.br_add_stl(selectedFile.toString())
    }
    FileDialog {
        id: exportStlDialog
        title: "Export the target as STL (mm; the components in the calculations)"
        fileMode: FileDialog.SaveFile
        defaultSuffix: "stl"
        nameFilters: ["STL mesh (*.stl)", "All files (*)"]
        onAccepted: Julia.br_export_stl(selectedFile.toString())
    }
    // VisRad workspaces: import a design; export its beams into a copy of the imported one
    FileDialog {
        id: vrwDialog
        property string mode: "import"
        function run(m) { mode = m; open() }
        title: mode === "import" ? "Import a VisRad workspace" : "Write the beams into a copy of the VisRad workspace"
        fileMode: mode === "import" ? FileDialog.OpenFile : FileDialog.SaveFile
        defaultSuffix: "vrw"
        nameFilters: ["VisRad workspace (*.vrw)", "All files (*)"]
        onAccepted: {
            const url = selectedFile.toString()
            if (mode === "import") Julia.br_import_visrad_workspace(url)
            else Julia.br_export_visrad_workspace(url)
        }
    }
    FolderDialog {
        id: resultsDialog
        title: "Folder for the results (the drive, and the drive history when there is one)"
        onAccepted: Julia.br_export_results(selectedFolder.toString())
    }
    // retune the enabled beams' pulses to a reference shot's drive (capsule, face or foil)
    Dialog {
        id: retuneDialog
        title: "Retune to a reference drive"
        // what the retune holds: capsule, face or foil (the level's units), the first the
        // design offers unless several are chosen between
        readonly property int choice: app.retuneOn.length > 1 ? Math.max(retuneOn.currentIndex, 0) : 0
        readonly property string kind: app.retuneKinds.length > choice ? app.retuneKinds[choice] : ""
        anchors.centerIn: parent
        modal: true
        standardButtons: Dialog.Close
        FileDialog {
            id: retunePick
            property bool saving: false
            title: saving ? "Write the retuned pulses" : "Reference pulse table"
            fileMode: saving ? FileDialog.SaveFile : FileDialog.OpenFile
            defaultSuffix: "csv"
            nameFilters: ["Pulse table (*.csv)", "All files (*)"]
            onAccepted: {
                const url = selectedFile.toString()
                if (saving) retuneOut.text = url
                else {
                    retuneTable.text = url
                    if (retuneOut.text === "") retuneOut.text = url.replace(/\.csv$/i, "") + "_retuned.csv"
                }
            }
        }
        ColumnLayout {
            spacing: 8
            RowLayout {
                Label { text: "The enabled beams' pulses, retuned to a reference shot's drive" }
                Info { topic: "drive.retune" }
            }
            GridLayout {
                columns: 3
                Label { text: "Reference pulses" }
                TextField { id: retuneTable; Layout.preferredWidth: 340; placeholderText: "pulse table (CSV)"
                            tip: text !== "" ? text : "The reference shot's pulse table, all beams firing" }
                Button { text: "Browse…"; onClicked: { retunePick.saving = false; retunePick.open() }
                         tip: "Choose the reference pulse table" }
                Label { text: "Retuned pulses" }
                TextField { id: retuneOut; Layout.preferredWidth: 340; placeholderText: "file to write"
                            tip: text !== "" ? text : "Where the retuned pulse table is written" }
                Button { text: "Browse…"; onClicked: { retunePick.saving = true; retunePick.open() }
                         tip: "Choose the file to write" }
                Label { text: "One power per" }
                ComboBox {
                    id: retuneBy
                    Layout.preferredWidth: 160
                    textRole: "label"
                    valueRole: "key"
                    model: [{ label: "quad", key: "quad" }, { label: "cone group", key: "group" },
                            { label: "beam", key: "beam" }]
                    tip: "The beams that share one retuned power"
                }
                Item { width: 1 }
                Label { text: "Power limit (× reference)" }
                TextField { id: retuneCap; text: "1.3"; Layout.preferredWidth: 100
                            tip: "Upper limit of each retuned power, as a multiple of its reference" }
                Item { width: 1 }
                Label { text: "Stay near the reference" }
                TextField { id: retuneStay; text: "0"; Layout.preferredWidth: 100
                            tip: "Weight on staying near the reference powers (0: none; try 1e-4)" }
                Item { width: 1 }
                Label { text: "Holds"; visible: retuneOn.visible }
                ComboBox {
                    id: retuneOn
                    visible: app.retuneOn.length > 1
                    Layout.preferredWidth: 220
                    model: app.retuneLabels
                    tip: "Match the observed surface's x-rays, or the foil's laser spot"
                }
                Item { width: 1; visible: retuneOn.visible }
                Label { text: "Target" }
                ComboBox {
                    id: retuneTarget
                    Layout.preferredWidth: 160
                    textRole: "label"
                    valueRole: "key"
                    model: [{ label: "reference shape", key: "reference" }, { label: "uniform", key: "uniform" }]
                    tip: "Match the reference's drive as it is, or aim for a uniform drive"
                }
                Item { width: 1 }
                Label { text: "Level"; visible: retuneTarget.currentValue === "uniform" }
                TextField {
                    id: retuneLevel
                    visible: retuneTarget.currentValue === "uniform"
                    Layout.preferredWidth: 100
                    placeholderText: "reference's"
                    tip: "The level to hold at every time; empty: the reference's"
                }
                Label { text: retuneDialog.kind === "capsule" ? "Tᵣ (eV)" : "TW/mm²"; visible: retuneLevel.visible }
                Item { width: 1 }
                CheckBox {
                    id: retuneLoad
                    text: "Fire the retuned pulses"
                    checked: true
                    Layout.columnSpan: 2
                    tip: "Load the retuned pulses as the session's pulses"
                }
            }
            RowLayout {
                Button {
                    text: "Retune"
                    enabled: app.driveCurrent && !app.retuneRunning && retuneTable.text !== "" && retuneOut.text !== ""
                    onClicked: Julia.br_retune(retuneTable.text, retuneOut.text, retuneBy.currentValue,
                                               retuneCap.text, retuneStay.text,
                                               retuneOn.visible ? app.retuneOn[retuneDialog.choice] : "auto",
                                               retuneTarget.currentValue, retuneLevel.text, retuneLoad.checked)
                    tip: app.driveCurrent ? "Runs in the background" : "Compute the drive first"
                }
                BusyIndicator { running: app.retuneRunning; visible: running; Layout.preferredHeight: 24; Layout.preferredWidth: 24 }
                Label { text: app.retuneStatus; Layout.fillWidth: true; Layout.maximumWidth: 400; wrapMode: Text.WordWrap }
            }
        }
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
    TransformDialog { id: transformDialog; x: (win.width - width) / 2; y: 60 }
    // the view from a beam or diagnostic port: angles typed, or a table beam's port
    Dialog {
        id: portDialog
        title: "Look from port"
        anchors.centerIn: parent
        standardButtons: Dialog.Ok | Dialog.Cancel
        GridLayout {
            columns: 2
            Label { text: "Beam port" }
            ComboBox {
                Layout.preferredWidth: 180
                model: app.beamNames
                onActivated: (index) => {
                    const a = Julia.br_port_angles(index)
                    portTheta.text = a[0]
                    portPhi.text = a[1]
                }
            }
            Label { text: "θ (°)" }
            TextField { id: portTheta; text: "90"; Layout.preferredWidth: 100 }
            Label { text: "φ (°)" }
            TextField { id: portPhi; text: "0"; Layout.preferredWidth: 100 }
        }
        onAccepted: Julia.br_look_from_port(portTheta.text, portPhi.text)
    }
    // the directories diagnostic definitions are read from (saved for later sessions)
    // and the problems found in their definitions (check_diagnostics)
    Dialog {
        id: diagnosticsDialog
        title: "Diagnostic directories: " + app.diagCount + " definitions found"
        anchors.centerIn: parent
        width: Math.min(win.width - 40, 760)
        height: Math.min(win.height - 40, 520)
        standardButtons: Dialog.Close
        FolderDialog {
            id: diagFolderDialog
            title: "Directory of diagnostic definitions (.toml)"
            onAccepted: Julia.br_diag_add_folder(selectedFolder.toString())
        }
        ColumnLayout {
            anchors.fill: parent
            spacing: 4
            RowLayout {
                SectionLabel { text: "Directories"; Layout.fillWidth: true }
                Info { topic: "diagnostics.folders" }
                Button { text: "Add directory…"; onClicked: diagFolderDialog.open()
                         tip: "Read definitions from another directory, kept for later sessions" }
                Button { text: "Refresh"; onClicked: Julia.br_diag_reload(); tip: "Read every definition again" }
            }
            Label { visible: app.diagFolders === 0; text: "Built-in definitions only"; opacity: 0.6 }
            Repeater {
                model: folderModel
                delegate: RowLayout {
                    required property string path
                    required property bool saved
                    Layout.fillWidth: true
                    Cell { text: path; elide: Text.ElideMiddle; Layout.fillWidth: true }
                    Label { visible: !saved; text: "environment"; opacity: 0.6 }
                    Button { visible: saved; text: "Remove"; flat: true; onClicked: Julia.br_diag_remove_folder(path)
                             tip: "Stop reading this directory" }
                }
            }
            SectionLabel { text: "Problems (" + app.diagIssues + ")" }
            ListView {
                Layout.fillWidth: true
                Layout.fillHeight: true
                clip: true
                model: issueModel
                ScrollBar.vertical: ScrollBar {}
                delegate: RowLayout {
                    required property string level
                    required property string where_
                    required property string message
                    required property string file
                    width: ListView.view.width
                    Label { text: level; Layout.preferredWidth: 60; opacity: 0.75
                            color: level === "error" ? "#D55E00" : palette.windowText }
                    Cell { text: where_ + ": " + message; tip: file; Layout.fillWidth: true }
                }
            }
        }
    }

    Timer {
        interval: 16
        running: true
        repeat: true
        onTriggered: {
            const r = Julia.br_tick()
            if (app.cadPrompt !== "" && !cadDialog.opened) cadDialog.open()
            if (app.grabPath !== "") {
                const path = app.grabPath
                Julia.br_grab_started()
                body.grabToImage(function (img) { img.saveToFile(path) })
            }
            if (r === 2) {
                win.close()
            } else {
                // 1: the 3D view or the map changed, 4: the time chart, 8: the sweep chart,
                // 16: the Diagnostics tab's chart or image
                if (r & 1) {
                    canvas.update()
                    if (app.sideView) sideCanvas.update()
                }
                if (r & 4) results.redrawTime()
                if (r & 8) results.redrawSweep()
                if (r & 16) results.redrawDiagnostics()
            }
        }
    }
}
