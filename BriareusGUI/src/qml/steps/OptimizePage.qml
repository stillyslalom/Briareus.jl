import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml

StepPage {
    id: page
    title: "Optimize"
    guidance: "Free the pointing fields to optimize per cone, quad or beam within ± bounds " +
              "around the current table (power bounds are a fraction of the start), choose " +
              "an objective, and run. L-BFGS-B with exact (forward-mode) gradients runs in " +
              "the background; canceling keeps the best design found so far. Apply the result " +
              "to put it into the beam table."

    GridLayout {
        columns: 2
        columnSpacing: 16
        rowSpacing: 8
        enabled: !app.optRunning
        Label { text: "Parameters per" }
        ComboBox {
            Layout.preferredWidth: 200
            textRole: "text"; valueRole: "value"
            model: [{ text: "cone", value: "group" }, { text: "quad", value: "quad" },
                    { text: "beam", value: "beam" }]
            currentIndex: Math.max(0, ["group", "quad", "beam"].indexOf(app.optBy))
            onActivated: Julia.br_opt_set_by(currentValue)
        }
    }

    Label { text: "Free fields and bounds"; font.weight: Font.DemiBold; Layout.topMargin: 8 }
    GridLayout {
        columns: 3
        columnSpacing: 12
        rowSpacing: 4
        enabled: !app.optRunning
        Repeater {
            model: [{ key: "z", label: "Aim z", unit: "± mm" },
                    { key: "r", label: "Aim r", unit: "± mm" },
                    { key: "phi", label: "Tangential offset", unit: "± mm" },
                    { key: "defocus", label: "Defocus", unit: "± mm" },
                    { key: "power", label: "Power", unit: "± fraction" }]
            delegate: RowLayout {
                Layout.columnSpan: 3
                spacing: 12
                CheckBox {
                    text: modelData.label
                    Layout.preferredWidth: 180
                    checked: app["optFree_" + modelData.key]
                    onToggled: Julia.br_opt_free(modelData.key, checked)
                }
                TextField {
                    Layout.preferredWidth: 80
                    text: app["optDelta_" + modelData.key]
                    onEditingFinished: if (text !== app["optDelta_" + modelData.key])
                                           Julia.br_opt_edit(modelData.key, text)
                }
                Label { text: modelData.unit; opacity: 0.75 }
            }
        }
    }

    Label { text: "Objective"; font.weight: Font.DemiBold; Layout.topMargin: 8 }
    GridLayout {
        columns: 2
        columnSpacing: 16
        rowSpacing: 6
        enabled: !app.optRunning
        Label { text: "Minimize" }
        ComboBox {
            Layout.preferredWidth: 320
            textRole: "text"; valueRole: "value"
            model: [{ text: "Capsule symmetry, axisymmetric (P_l/P0)", value: "legendre" },
                    { text: "Capsule symmetry, 3D (all m)", value: "ylm" },
                    { text: "Foil spot uniformity", value: "irradiance" }]
            currentIndex: Math.max(0, ["legendre", "ylm", "irradiance"].indexOf(app.optObjective))
            onActivated: Julia.br_opt_objective(currentValue)
        }
        Label { text: "Mode weights"; visible: app.optObjective !== "irradiance" }
        RowLayout {
            visible: app.optObjective !== "irradiance"
            Repeater {
                model: [1, 2, 3, 4]
                delegate: RowLayout {
                    Label { text: "l = " + modelData }
                    TextField {
                        Layout.preferredWidth: 56
                        text: app["optW" + modelData]
                        onEditingFinished: if (text !== app["optW" + modelData])
                                               Julia.br_opt_edit("w" + modelData, text)
                    }
                }
            }
        }
        Label { text: "Hold the drive"; visible: app.optObjective !== "irradiance" }
        Switch {
            visible: app.optObjective !== "irradiance"
            checked: app.optHold
            onToggled: Julia.br_opt_hold(checked)
            ToolTip.visible: hovered
            ToolTip.text: "Penalize a change of the mean capsule flux from its start value"
        }
        Label { text: "LEH keep-out"; visible: app.optObjective !== "irradiance" }
        RowLayout {
            visible: app.optObjective !== "irradiance"
            Switch {
                checked: app.optKeepout
                onToggled: Julia.br_opt_keepout(checked)
                ToolTip.visible: hovered
                ToolTip.text: "Penalize beams whose worst-case envelope (1 % of peak, with pointing, " +
                              "positioning and focus errors) comes closer to an LEH rim than this"
            }
            TextField {
                Layout.preferredWidth: 70
                enabled: app.optKeepout
                text: app.optKeepoutMin
                onEditingFinished: if (text !== app.optKeepoutMin) Julia.br_opt_edit("keepout_min", text)
            }
            Label { text: "mm minimum clearance"; opacity: 0.75 }
        }
        Label { text: "Spill weight"; visible: app.optObjective === "irradiance" }
        RowLayout {
            visible: app.optObjective === "irradiance"
            TextField {
                Layout.preferredWidth: 70
                text: app.optSpill
                onEditingFinished: if (text !== app.optSpill) Julia.br_opt_edit("spill", text)
            }
            Label { text: "on power outside the spot (radius set on Drive)"; opacity: 0.75 }
        }
        Label { text: "Iterations" }
        SpinBox {
            from: 1; to: 2000; editable: true
            value: app.optIterations
            onValueModified: Julia.br_opt_edit("iterations", value)
        }
    }
    Label {
        text: app.optError
        visible: text !== ""
        color: "#c42b1c"
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }

    RowLayout {
        spacing: 8
        Layout.topMargin: 8
        Button {
            text: "Optimize"
            highlighted: true
            enabled: !app.optRunning && app.driveBuilt && !app.driveStale
            onClicked: Julia.br_opt_start()
        }
        Button { text: "Cancel"; enabled: app.optRunning; onClicked: Julia.br_opt_cancel() }
        BusyIndicator {
            running: app.optRunning
            visible: running
            Layout.preferredHeight: 32
            Layout.preferredWidth: 32
        }
    }
    Label {
        text: !app.driveBuilt || app.driveStale ? "Compute the drive first (Drive step)." : app.optProgress
        visible: text !== ""
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }

    // result: before/after and the parameters
    Label { text: "Result"; font.weight: Font.DemiBold; visible: app.optHasResult; Layout.topMargin: 8 }
    Label { text: app.optValues; visible: app.optHasResult; wrapMode: Text.WordWrap; Layout.fillWidth: true }
    Label { text: app.optClipped; visible: app.optHasResult; opacity: 0.75 }
    Label { text: app.optSpot; visible: app.optSpot !== ""; opacity: 0.75 }
    GridLayout {
        visible: app.optHasResult
        columns: 3
        columnSpacing: 24
        rowSpacing: 4
        Label { text: "l"; font.weight: Font.DemiBold; visible: optModeRepeater.count > 0 }
        Label { text: "ε_l before"; font.weight: Font.DemiBold; visible: optModeRepeater.count > 0 }
        Label { text: "ε_l after"; font.weight: Font.DemiBold; visible: optModeRepeater.count > 0 }
        Repeater {
            id: optModeRepeater
            model: optModeModel
            delegate: Label { Layout.row: model.l; Layout.column: 0; text: model.l }
        }
        Repeater {
            model: optModeModel
            delegate: Label { Layout.row: model.l; Layout.column: 1; text: model.before }
        }
        Repeater {
            model: optModeModel
            delegate: Label { Layout.row: model.l; Layout.column: 2; text: model.after }
        }
    }
    RowLayout {
        visible: app.optHasResult
        spacing: 8
        Button {
            text: app.optApplied ? "Applied" : "Apply to the beam table"
            enabled: !app.optApplied && !app.optRunning
            onClicked: Julia.br_opt_accept()
        }
    }
    GridLayout {
        visible: app.optHasResult
        columns: 5
        columnSpacing: 16
        rowSpacing: 2
        Label { text: "Parameter"; font.weight: Font.DemiBold }
        Label { text: "Start"; font.weight: Font.DemiBold }
        Label { text: "Optimized"; font.weight: Font.DemiBold }
        Label { text: "Lower"; font.weight: Font.DemiBold }
        Label { text: "Upper"; font.weight: Font.DemiBold }
        Repeater {
            model: paramModel
            delegate: Label { Layout.row: index + 1; Layout.column: 0; text: model.name }
        }
        Repeater {
            model: paramModel
            delegate: Label { Layout.row: index + 1; Layout.column: 1; text: model.start }
        }
        Repeater {
            model: paramModel
            delegate: Label {
                Layout.row: index + 1; Layout.column: 2
                text: model.optimized + (model.atBound ? "  (at bound)" : "")
                color: model.atBound ? "#d08a00" : palette.windowText
            }
        }
        Repeater {
            model: paramModel
            delegate: Label { Layout.row: index + 1; Layout.column: 3; text: model.lower; opacity: 0.75 }
        }
        Repeater {
            model: paramModel
            delegate: Label { Layout.row: index + 1; Layout.column: 4; text: model.upper; opacity: 0.75 }
        }
    }
}
