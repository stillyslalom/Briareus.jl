// The toolbar: views, display, and the analyses. Everything here is also in the menus.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQml.Models
import jlqml
import "common"

ToolBar {
    id: bar
    // `app.cutSettingsRequest` opens the Cut settings popup and, with `app.cutGrabPath`
    // set, saves an image of it (scripted screenshots)
    readonly property bool cutSettingsRequested: app.cutSettingsRequest === true
    onCutSettingsRequestedChanged: if (cutSettingsRequested) {
        cutPopup.open()
        if (app.cutGrabPath !== "") cutGrabTimer.start()
    }
    Timer {
        id: cutGrabTimer
        interval: 500
        onTriggered: {
            const path = app.cutGrabPath
            cutPopup.contentItem.parent.grabToImage(function (img) { img.saveToFile(path) })
        }
    }
    // the tools wrap to a second row when the window is too narrow for one
    topPadding: 2
    bottomPadding: 2
    contentHeight: tools.implicitHeight
    // the tools' height, not the style's toolbar background (48 px)
    implicitHeight: contentHeight + topPadding + bottomPadding
    Flow {
        id: tools
        x: 4
        width: parent.width - 8
        spacing: 1

        component Tool: ToolButton {
            property string tip: ""
            Tip { visible: parent.hovered && parent.tip !== ""; text: parent.tip }
            focusPolicy: Qt.NoFocus
            topPadding: 3
            bottomPadding: 3
            leftPadding: 8
            rightPadding: 8
            implicitHeight: 28
        }

        Tool { text: "Fit"; tip: "Fit the whole target (F in the view)"; onClicked: Julia.br_fit_all() }
        Tool { text: "Fit sel."; tip: "Fit the selection (or double-click it)"; onClicked: Julia.br_fit_selection() }
        Tool { text: "Iso"; tip: "Isometric view"; onClicked: Julia.br_standard_view("iso") }
        Tool { text: "Front"; tip: "Look along +y"; onClicked: Julia.br_standard_view("front") }
        Tool { text: "Side"; tip: "Look along −x"; onClicked: Julia.br_standard_view("right") }
        Tool { text: "Top"; tip: "Look down −z"; onClicked: Julia.br_standard_view("top") }
        Tool {
            text: "Along beam"
            tip: "Look down the selected beam from its port"
            enabled: ["group", "quad", "beam"].indexOf(app.selKind) >= 0
            onClicked: Julia.br_look_along_beam()
        }
        ToolSeparator { height: 28 }
        Tool {
            text: "Ortho"; checkable: true; checked: app.ortho
            tip: "Orthographic projection (otherwise perspective)"
            onClicked: Julia.br_projection(checked)
        }
        Tool {
            id: cutTool
            text: "Cutaway"; checkable: true; checked: app.cutaway
            tip: "Cut the target open to see inside (right-click: options)"
            onClicked: Julia.br_set_cutaway(checked)
            TapHandler {
                acceptedButtons: Qt.RightButton
                onTapped: cutMenu.popup()
            }
            AutoMenu {
                id: cutMenu
                readonly property bool sector: app.cutKind === "sector"
                MenuItem { text: "Quadrant (90°)"; checkable: true; checked: cutMenu.sector && app.cutWidth === 90; onTriggered: Julia.br_cut_width(90) }
                MenuItem { text: "Half (180°)"; checkable: true; checked: cutMenu.sector && app.cutWidth === 180; onTriggered: Julia.br_cut_width(180) }
                MenuItem { text: "Three quarters (270°)"; checkable: true; checked: cutMenu.sector && app.cutWidth === 270; onTriggered: Julia.br_cut_width(270) }
                MenuSeparator {}
                MenuItem { text: "Facing the current view"; onTriggered: Julia.br_cut_toward_view() }
                MenuItem { text: "Facing −x −y"; checkable: true; checked: cutMenu.sector && app.cutAxis === "z" && app.cutAzimuth === 225; onTriggered: Julia.br_cut_azimuth(225) }
                MenuItem { text: "Facing −y"; checkable: true; checked: cutMenu.sector && app.cutAxis === "z" && app.cutAzimuth === 270; onTriggered: Julia.br_cut_azimuth(270) }
                MenuItem { text: "Facing +x"; checkable: true; checked: cutMenu.sector && app.cutAxis === "z" && app.cutAzimuth === 0; onTriggered: Julia.br_cut_azimuth(0) }
                MenuSeparator {}
                MenuItem { text: "Section plane x (removes +x)"; checkable: true; checked: !cutMenu.sector && app.cutNormal === "1, 0, 0"; onTriggered: Julia.br_cut_plane("x") }
                MenuItem { text: "Section plane y (removes +y)"; checkable: true; checked: !cutMenu.sector && app.cutNormal === "0, 1, 0"; onTriggered: Julia.br_cut_plane("y") }
                MenuItem { text: "Section plane z (removes +z)"; checkable: true; checked: !cutMenu.sector && app.cutNormal === "0, 0, 1"; onTriggered: Julia.br_cut_plane("z") }
                MenuItem { text: "Section plane facing the view"; onTriggered: Julia.br_cut_plane("view") }
                MenuSeparator {}
                MenuItem { text: "Cut settings…"; onTriggered: cutPopup.open() }
                MenuSeparator {}
                AutoMenu {
                    id: cutCompMenu
                    title: "Components cut"
                    Instantiator {
                        model: app.cutChoices
                        delegate: MenuItem {
                            required property string modelData
                            required property int index
                            text: app.cutLabels[index] || modelData
                            checkable: true
                            checked: app.cutSet.indexOf(modelData) >= 0
                            onTriggered: Julia.br_cut_component(modelData, checked)
                        }
                        onObjectAdded: (index, object) => cutCompMenu.insertItem(index, object)
                        onObjectRemoved: (index, object) => cutCompMenu.removeItem(object)
                    }
                }
            }
            // every setting of the cut, kept open while the view is turned
            Popup {
                id: cutPopup
                y: cutTool.height + 4
                padding: 12
                closePolicy: Popup.CloseOnEscape
                ColumnLayout {
                    spacing: 6
                    RowLayout {
                        Label { text: "Cut"; font.bold: true }
                        Info { topic: "view.cutaway" }
                        Item { Layout.fillWidth: true }
                        RadioButton { text: "Sector"; checked: app.cutKind === "sector"; onToggled: if (checked) Julia.br_cut_kind("sector")
                                      Tip { visible: parent.hovered; text: "Remove a slice about a coordinate axis" } }
                        RadioButton { text: "Plane"; checked: app.cutKind === "plane"; onToggled: if (checked) Julia.br_cut_kind("plane")
                                      Tip { visible: parent.hovered; text: "Remove the side a plane's normal points to" } }
                    }
                    ColumnLayout {
                        visible: app.cutKind === "sector"
                        spacing: 4
                        RowLayout {
                            Label { text: "About axis"; Layout.preferredWidth: 140 }
                            Repeater {
                                model: ["x", "y", "z"]
                                RadioButton {
                                    required property string modelData
                                    text: modelData
                                    checked: app.cutAxis === modelData
                                    onToggled: if (checked) Julia.br_cut_axis(modelData)
                                    Tip {
                                        visible: parent.hovered
                                        text: "The slice turns about the chamber's " + modelData + " axis"
                                    }
                                }
                            }
                        }
                        RowLayout {
                            Label { text: "Center azimuth"; Layout.preferredWidth: 140 }
                            Slider {
                                from: 0; to: 360; stepSize: 5
                                value: app.cutAzimuth
                                Layout.preferredWidth: 160
                                onMoved: Julia.br_cut_azimuth(value)
                                Tip {
                                    visible: parent.hovered && !parent.pressed
                                    text: "The azimuth the removed slice is centered on (°)"
                                }
                            }
                            Label { text: Math.round(app.cutAzimuth) + "°"; Layout.preferredWidth: 40 }
                            Info { topic: "view.cut" }
                        }
                        RowLayout {
                            Label { text: "Width"; Layout.preferredWidth: 140 }
                            Slider {
                                from: 5; to: 360; stepSize: 5
                                value: app.cutWidth
                                Layout.preferredWidth: 160
                                onMoved: Julia.br_cut_width(value)
                                Tip {
                                    visible: parent.hovered && !parent.pressed
                                    text: "The azimuthal width of the removed slice (°)"
                                }
                            }
                            Label { text: Math.round(app.cutWidth) + "°"; Layout.preferredWidth: 40 }
                        }
                        PropRow {
                            label: "From (along the axis)"; value: app.cutLo; unit: "mm (empty: no bound)"
                            tip: "The slice's lower bound along the axis"
                            onEdited: (text) => Julia.br_cut_range(text, app.cutHi)
                        }
                        PropRow {
                            label: "To"; value: app.cutHi; unit: "mm (empty: no bound)"
                            tip: "The slice's upper bound along the axis"
                            onEdited: (text) => Julia.br_cut_range(app.cutLo, text)
                        }
                    }
                    ColumnLayout {
                        visible: app.cutKind === "plane"
                        spacing: 4
                        RowLayout {
                            Cell { text: "Normal"; Layout.preferredWidth: 140
                                   tip: "The cut removes the side the normal (" + app.cutNormal + ") points to" }
                            Button { text: "x"; implicitWidth: 36; onClicked: Julia.br_cut_plane("x"); tip: "Normal along +x: removes +x" }
                            Button { text: "y"; implicitWidth: 36; onClicked: Julia.br_cut_plane("y"); tip: "Normal along +y: removes +y" }
                            Button { text: "z"; implicitWidth: 36; onClicked: Julia.br_cut_plane("z"); tip: "Normal along +z: removes +z" }
                            Button { text: "View"; onClicked: Julia.br_cut_plane("view"); tip: "Normal toward the camera: removes the near side" }
                            Button { text: "Flip"; onClicked: Julia.br_cut_plane("flip"); tip: "Reverse the normal: removes the other side" }
                        }
                        RowLayout {
                            Label { text: "Normal θ, φ"; Layout.preferredWidth: 140 }
                            TextField {
                                id: cutThetaField; text: app.cutTheta; Layout.preferredWidth: 70; selectByMouse: true
                                onEditingFinished: if (text !== app.cutTheta) Julia.br_cut_angles(text, cutPhiField.text)
                                tip: "The normal's polar angle (°)"
                            }
                            TextField {
                                id: cutPhiField; text: app.cutPhi; Layout.preferredWidth: 70; selectByMouse: true
                                onEditingFinished: if (text !== app.cutPhi) Julia.br_cut_angles(cutThetaField.text, text)
                                tip: "The normal's azimuth (°)"
                            }
                            Label { text: "°"; opacity: 0.7 }
                        }
                        PropRow {
                            label: "Normal vector"; value: app.cutNormal; unit: "x, y, z"; fieldWidth: 150
                            tip: "The normal's direction (any length)"
                            onEdited: (text) => Julia.br_cut_normal(text)
                        }
                        RowLayout {
                            Label { text: "Through point"; Layout.preferredWidth: 140 }
                            Repeater {
                                id: cutPointFields
                                model: 3
                                TextField {
                                    required property int index
                                    text: app.cutPoint ? app.cutPoint[index] : "0"
                                    Layout.preferredWidth: 72
                                    selectByMouse: true
                                    onEditingFinished: {
                                        const p = [0, 1, 2].map(i => cutPointFields.itemAt(i).text)
                                        if (text !== app.cutPoint[index]) Julia.br_cut_point(p[0], p[1], p[2])
                                    }
                                    tip: "The plane's point, " + ["x", "y", "z"][index] + " (mm)"
                                }
                            }
                            Label { text: "mm"; opacity: 0.7 }
                            Button { text: "TCC"; onClicked: Julia.br_cut_point("0", "0", "0"); tip: "The plane through the chamber center" }
                        }
                        Button {
                            text: app.cutSelected !== "" ? "Through " + app.cutSelected : "Through the selected component"
                            enabled: app.cutSelected !== ""
                            tip: "A plane through the selected component"
                            onClicked: Julia.br_cut_through_selected()
                        }
                        RowLayout {
                            Label { text: "Offset"; Layout.preferredWidth: 140 }
                            Slider {
                                from: -app.cutScale; to: app.cutScale
                                value: app.cutOffset
                                Layout.preferredWidth: 160
                                onMoved: Julia.br_cut_offset(value)
                                Tip {
                                    visible: parent.hovered && !parent.pressed
                                    text: "Moves the plane along its normal from the point (mm)"
                                }
                            }
                            Label { text: app.cutOffset.toFixed(2) + " mm from the point" }
                        }
                        MenuSeparator { Layout.fillWidth: true }
                        RowLayout {
                            Label { text: "Direction"; font.bold: true; Layout.preferredWidth: 140 }
                            ComboBox {
                                Layout.preferredWidth: 150
                                model: app.beamNames
                                displayText: "Beam port…"
                                onActivated: (index) => {
                                    const a = Julia.br_port_angles(index)
                                    dirTheta.text = a[0]
                                    dirPhi.text = a[1]
                                }
                                tip: "Take θ, φ from a beam's port"
                            }
                        }
                        RowLayout {
                            Label { text: "θ, φ"; Layout.preferredWidth: 140 }
                            TextField { id: dirTheta; text: app.cutLineTheta; Layout.preferredWidth: 70; selectByMouse: true
                                        tip: "The direction's polar angle (°)" }
                            TextField { id: dirPhi; text: app.cutLinePhi; Layout.preferredWidth: 70; selectByMouse: true
                                        tip: "The direction's azimuth (°)" }
                            Label { text: "°"; opacity: 0.7 }
                        }
                        RowLayout {
                            Label { text: ""; Layout.preferredWidth: 140 }
                            Button {
                                text: "Facing it"
                                tip: "The normal along the direction: removes the side toward it"
                                onClicked: Julia.br_cut_angles(dirTheta.text, dirPhi.text)
                            }
                            Button {
                                text: "Containing the line of sight"
                                tip: "The plane through the line from the origin toward the direction"
                                onClicked: Julia.br_cut_line(dirTheta.text, dirPhi.text, losFromSel.checked, app.cutRotation)
                            }
                        }
                        RowLayout {
                            Label { text: "Line from"; Layout.preferredWidth: 140 }
                            RadioButton { text: "TCC"; checked: app.cutAtTcc || app.cutSelected === ""
                                          Tip { visible: parent.hovered; text: "The line starts at the chamber center" } }
                            RadioButton {
                                id: losFromSel; text: "Selected component"
                                enabled: app.cutSelected !== ""; checked: !app.cutAtTcc && app.cutSelected !== ""
                                Tip { visible: parent.hovered; text: "The line of sight starts at the selected component" }
                            }
                        }
                        RowLayout {
                            enabled: app.cutLine
                            Label { text: "Rotation about it"; Layout.preferredWidth: 140 }
                            Slider {
                                from: 0; to: 360; stepSize: 5
                                value: app.cutRotation
                                Layout.preferredWidth: 160
                                onMoved: Julia.br_cut_rotation(value)
                                Tip {
                                    visible: parent.hovered && !parent.pressed
                                    text: "Turns the plane about the line of sight (°)"
                                }
                            }
                            Label { text: Math.round(app.cutRotation) + "°"; Layout.preferredWidth: 40 }
                        }
                        Cell {
                            visible: app.cutLine
                            text: "Line of sight: " + app.cutLineText
                            opacity: 0.7
                            Layout.maximumWidth: 520
                        }
                    }
                    Button { text: "Close"; Layout.alignment: Qt.AlignRight; onClicked: cutPopup.close() }
                }
            }
        }
        ToolSeparator { height: 28 }
        Label { text: "Beams"; opacity: 0.7; leftPadding: 2; rightPadding: 4; height: beamStyleBox.height; verticalAlignment: Text.AlignVCenter }
        ComboBox {
            id: beamStyleBox
            width: 96
            height: 28                                  // the tool buttons' height
            focusPolicy: Qt.NoFocus
            textRole: "text"; valueRole: "value"
            model: [{ text: "Rays", value: "rays" }, { text: "Volumes", value: "volumes" },
                    { text: "Hidden", value: "off" }]
            currentIndex: Math.max(0, ["rays", "volumes", "off"].indexOf(app.beamStyle))
            onActivated: Julia.br_beam_style(currentValue)
            tip: "Draw beams as rays, as volumes, or not at all"
        }
        Label { text: "Color"; opacity: 0.7; leftPadding: 8; rightPadding: 4; height: beamStyleBox.height; verticalAlignment: Text.AlignVCenter }
        ComboBox {
            width: 160
            height: 28
            focusPolicy: Qt.NoFocus
            textRole: "text"; valueRole: "value"
            model: [{ text: "Components", value: "component" },
                    { text: "Absorbed laser", value: "laser" },
                    { text: "X-ray radiosity", value: "radiosity" },
                    { text: "Incident x-ray flux", value: "incident" }]
            currentIndex: Math.max(0, ["component", "laser", "radiosity", "incident"].indexOf(app.colorBy))
            onActivated: Julia.br_set_color(currentValue)
            tip: "Color by component, or by a drive field (W/cm²)"
        }
        Tool {
            text: app.observedKind === "face" ? "Face map" : app.observedKind === "capsule" ? "Capsule map" : "Foil map"
            checkable: true
            checked: app.driveView === "map"
            enabled: app.driveHasResult && (app.hasCapsule || app.hasFoil)
            tip: "Show the flux map in place of the 3D view"
            onClicked: Julia.br_set_view(checked ? "map" : "scene")
        }
        ToolSeparator { height: 28 }
        Tool {
            text: app.driveCurrent ? "Drive ✓" : "Compute drive"
            enabled: !app.driveBuilding
            tip: "View factors and the laser and x-ray drive (Ctrl+R)"
            onClicked: Julia.br_build()
        }
        Tool { text: "Optimize…"; tip: "Pointing optimization settings"; onClicked: win.showTask("optimize") }
        Tool { text: "Clearance…"; tip: "Beam clearance check settings"; onClicked: win.showTask("clearance") }
        Tool { text: "Sweep…"; tip: "Placement sweep settings"; onClicked: win.showTask("sweep") }
        Tool { text: "TaLIS…"; tip: "TaLIS rule check settings"; onClicked: win.showTask("talis") }
    }
}
