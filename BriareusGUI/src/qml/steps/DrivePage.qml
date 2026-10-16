import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

StepPage {
    title: "Drive"
    guidance: "Compute the view factors for this target and beam set (seconds to tens of " +
              "seconds; once per target and beam ports). Pointing, power and on/off edits then " +
              "re-evaluate the laser deposition and x-ray drive at once. ε_l is the " +
              "rotation-invariant mode amplitude of the capsule flux over all m (equal to " +
              "|P_l/P0| for an axisymmetric drive)."

    RowLayout {
        spacing: 8
        Button {
            text: app.driveBuilt && !app.driveStale ? "Compute again" : "Compute drive"
            highlighted: !app.driveBuilt || app.driveStale
            enabled: !app.driveBuilding
            onClicked: Julia.br_build()
        }
        BusyIndicator {
            running: app.driveBuilding || app.driveEvaluating
            visible: running
            Layout.preferredHeight: 32
            Layout.preferredWidth: 32
        }
        Label {
            text: app.driveBuilding ? "computing view factors… " + app.driveElapsed :
                  app.driveStale ? "The target or the beam ports changed since the view factors were computed." :
                  app.driveStatus
            opacity: 0.75
            wrapMode: Text.WordWrap
            Layout.fillWidth: true
        }
    }

    GridLayout {
        columns: 2
        columnSpacing: 16
        rowSpacing: 8
        enabled: app.driveHasResult
        Label { text: "Viewer" }
        RowLayout {
            ButtonGroup { id: viewGroup }
            Button {
                text: "3D target"; checkable: true; ButtonGroup.group: viewGroup
                checked: app.driveView === "scene"
                onClicked: Julia.br_set_view("scene")
            }
            Button {
                text: app.hasFoil ? "Foil face" : "Capsule θ–φ map"; checkable: true
                ButtonGroup.group: viewGroup
                checked: app.driveView === "map"
                enabled: app.hasCapsule || app.hasFoil
                onClicked: Julia.br_set_view("map")
            }
        }
        Label { text: "Color the target by" }
        ComboBox {
            Layout.preferredWidth: 280
            textRole: "text"; valueRole: "value"
            model: [{ text: "Absorbed laser intensity", value: "laser" },
                    { text: "X-ray radiosity", value: "radiosity" },
                    { text: "Incident x-ray flux", value: "incident" }]
            currentIndex: Math.max(0, ["laser", "radiosity", "incident"].indexOf(app.driveField))
            onActivated: Julia.br_set_field(currentValue)
        }
        Label { text: "Spot radius"; visible: app.hasFoil }
        RowLayout {
            visible: app.hasFoil
            TextField {
                Layout.preferredWidth: 90
                text: app.spotRadius
                onEditingFinished: if (text !== app.spotRadius) Julia.br_set_spot_radius(text)
            }
            Label {
                text: "mm about the foil center"
                opacity: 0.75
                ToolTip.visible: radiusHover.hovered
                ToolTip.text: "The spot over which the irradiance statistics and the foil objective are taken"
                HoverHandler { id: radiusHover }
            }
        }
    }

    Label { text: app.driveSummary; font.weight: Font.DemiBold; visible: app.driveHasResult; Layout.topMargin: 8 }

    // capsule drive and modes
    Label { text: "Capsule drive"; font.weight: Font.DemiBold; visible: app.driveMean !== ""; Layout.topMargin: 8 }
    Label { text: "Mean incident flux " + app.driveMean; visible: app.driveMean !== "" }
    GridLayout {
        visible: app.driveMean !== ""
        columns: 3
        columnSpacing: 24
        rowSpacing: 4
        Label { text: "l"; font.weight: Font.DemiBold }
        Label { text: "ε_l"; font.weight: Font.DemiBold }
        Label { text: "P_l / P0"; font.weight: Font.DemiBold }
        Repeater {
            model: modeModel
            delegate: Label { Layout.row: model.l; Layout.column: 0; text: model.l }
        }
        Repeater {
            model: modeModel
            delegate: Label { Layout.row: model.l; Layout.column: 1; text: model.epsilon }
        }
        Repeater {
            model: modeModel
            delegate: Label { Layout.row: model.l; Layout.column: 2; text: model.legendre }
        }
    }

    Label { text: "Foil spot"; font.weight: Font.DemiBold; visible: app.spotText !== ""; Layout.topMargin: 8 }
    Label { text: app.spotText; visible: text !== ""; wrapMode: Text.WordWrap; Layout.fillWidth: true }

    // laser power accounting
    Label { text: "Laser power"; font.weight: Font.DemiBold; visible: app.driveHasResult; Layout.topMargin: 8 }
    Label { text: app.delivered; visible: app.driveHasResult; opacity: 0.75 }
    GridLayout {
        visible: app.driveHasResult
        columns: 3
        columnSpacing: 24
        rowSpacing: 4
        Repeater {
            model: powerModel
            delegate: Label { Layout.row: index; Layout.column: 0; text: model.name }
        }
        Repeater {
            model: powerModel
            delegate: Label { Layout.row: index; Layout.column: 1; text: model.power + " TW" }
        }
        Repeater {
            model: powerModel
            delegate: Label { Layout.row: index; Layout.column: 2; text: model.percent }
        }
    }
}
