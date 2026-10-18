// The toolbar: views, display, and the analyses. Everything here is also in the menus.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQml.Models
import jlqml
import "common"

ToolBar {
    id: bar
    RowLayout {
        anchors.fill: parent
        anchors.leftMargin: 8
        anchors.rightMargin: 8
        spacing: 2

        component Tool: ToolButton {
            property string tip: ""
            ToolTip.visible: hovered && tip !== ""
            ToolTip.delay: 500
            ToolTip.text: tip
            focusPolicy: Qt.NoFocus
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
        ToolSeparator {}
        Tool {
            text: "Ortho"; checkable: true; checked: app.ortho
            tip: "Orthographic projection (otherwise perspective)"
            onClicked: Julia.br_projection(checked)
        }
        Tool {
            id: cutTool
            text: "Cutaway"; checkable: true; checked: app.cutaway
            tip: "Cut away part of the target to see inside (right-click: where and what)"
            onClicked: Julia.br_set_cutaway(checked)
            TapHandler {
                acceptedButtons: Qt.RightButton
                onTapped: cutMenu.popup()
            }
            AutoMenu {
                id: cutMenu
                MenuItem { text: "Quadrant (90°)"; checkable: true; checked: app.cutWidth === 90; onTriggered: Julia.br_cut_width(90) }
                MenuItem { text: "Half (180°)"; checkable: true; checked: app.cutWidth === 180; onTriggered: Julia.br_cut_width(180) }
                MenuItem { text: "Three quarters (270°)"; checkable: true; checked: app.cutWidth === 270; onTriggered: Julia.br_cut_width(270) }
                MenuSeparator {}
                MenuItem { text: "Facing the current view"; onTriggered: Julia.br_cut_toward_view() }
                MenuItem { text: "Facing −x −y"; checkable: true; checked: app.cutAzimuth === 225; onTriggered: Julia.br_cut_azimuth(225) }
                MenuItem { text: "Facing −y"; checkable: true; checked: app.cutAzimuth === 270; onTriggered: Julia.br_cut_azimuth(270) }
                MenuItem { text: "Facing +x"; checkable: true; checked: app.cutAzimuth === 0; onTriggered: Julia.br_cut_azimuth(0) }
                MenuSeparator {}
                AutoMenu {
                    id: cutCompMenu
                    title: "Components cut"
                    Instantiator {
                        model: app.cutChoices
                        delegate: MenuItem {
                            required property string modelData
                            text: modelData
                            checkable: true
                            checked: app.cutSet.indexOf(modelData) >= 0
                            onTriggered: Julia.br_cut_component(modelData, checked)
                        }
                        onObjectAdded: (index, object) => cutCompMenu.insertItem(index, object)
                        onObjectRemoved: (index, object) => cutCompMenu.removeItem(object)
                    }
                }
            }
        }
        ToolSeparator {}
        Label { text: "Beams"; opacity: 0.7; leftPadding: 4; rightPadding: 4 }
        ComboBox {
            Layout.preferredWidth: 110
            focusPolicy: Qt.NoFocus
            textRole: "text"; valueRole: "value"
            model: [{ text: "Rays", value: "rays" }, { text: "Volumes", value: "volumes" },
                    { text: "Hidden", value: "off" }]
            currentIndex: Math.max(0, ["rays", "volumes", "off"].indexOf(app.beamStyle))
            onActivated: Julia.br_beam_style(currentValue)
        }
        Label { text: "Color"; opacity: 0.7; leftPadding: 8; rightPadding: 4 }
        ComboBox {
            Layout.preferredWidth: 180
            focusPolicy: Qt.NoFocus
            textRole: "text"; valueRole: "value"
            model: [{ text: "Components", value: "component" },
                    { text: "Absorbed laser", value: "laser" },
                    { text: "X-ray radiosity", value: "radiosity" },
                    { text: "Incident x-ray flux", value: "incident" }]
            currentIndex: Math.max(0, ["component", "laser", "radiosity", "incident"].indexOf(app.colorBy))
            onActivated: Julia.br_set_color(currentValue)
        }
        Tool {
            text: app.hasFoil ? "Foil map" : "Capsule map"; checkable: true
            checked: app.driveView === "map"
            enabled: app.driveHasResult && (app.hasCapsule || app.hasFoil)
            tip: "Show the capsule θ–φ flux map (or the foil face) instead of the 3D view"
            onClicked: Julia.br_set_view(checked ? "map" : "scene")
        }
        ToolSeparator {}
        Tool {
            text: app.driveCurrent ? "Drive ✓" : "Compute drive"
            enabled: !app.driveBuilding
            tip: "View factors and the laser and x-ray drive (Ctrl+R)"
            onClicked: Julia.br_build()
        }
        Tool { text: "Optimize…"; tip: "Pointing optimization settings"; onClicked: win.showTask("optimize") }
        Tool { text: "Clearance…"; tip: "Beam clearance check settings"; onClicked: win.showTask("clearance") }
        Tool { text: "TaLIS…"; tip: "TaLIS rule check settings"; onClicked: win.showTask("talis") }
        Item { Layout.fillWidth: true }
    }
}
