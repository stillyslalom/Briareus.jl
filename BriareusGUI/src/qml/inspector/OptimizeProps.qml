// The pointing optimization: free fields and bounds, objective, keep-outs, limits, run.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note {
        text: "Free pointing fields per cone, quad or beam within ± bounds around the beam " +
              "table (power bounds: a fraction of the start) and minimize an objective with " +
              "L-BFGS-B and exact gradients, in the background. Canceling keeps the best design."
    }
    ColumnLayout {
        spacing: 6
        enabled: !app.optRunning
        Layout.fillWidth: true
        RowLayout {
            spacing: 8
            Label { text: "Parameters per"; Layout.preferredWidth: 150 }
            ComboBox {
                Layout.fillWidth: true
                textRole: "text"; valueRole: "value"
                model: [{ text: "cone", value: "group" }, { text: "quad", value: "quad" },
                        { text: "beam", value: "beam" }]
                currentIndex: Math.max(0, ["group", "quad", "beam"].indexOf(app.optBy))
                onActivated: Julia.br_opt_set_by(currentValue)
            }
        }

        SectionLabel { text: "Free fields and bounds" }
        Repeater {
            model: [{ key: "z", label: "Aim z", unit: "± mm" },
                    { key: "r", label: "Aim r", unit: "± mm" },
                    { key: "phi", label: "Tangential offset", unit: "± mm" },
                    { key: "defocus", label: "Defocus", unit: "± mm" },
                    { key: "power", label: "Power", unit: "± fraction" }]
            delegate: RowLayout {
                spacing: 8
                CheckBox {
                    text: modelData.label
                    Layout.preferredWidth: 170
                    checked: app["optFree_" + modelData.key]
                    onToggled: Julia.br_opt_free(modelData.key, checked)
                }
                TextField {
                    Layout.preferredWidth: 90
                    text: app["optDelta_" + modelData.key]
                    onEditingFinished: if (text !== app["optDelta_" + modelData.key])
                                           Julia.br_opt_edit(modelData.key, text)
                }
                Label { text: modelData.unit; opacity: 0.7 }
            }
        }

        SectionLabel { text: "Objective" }
        ComboBox {
            Layout.fillWidth: true
            textRole: "text"; valueRole: "value"
            model: [{ text: "Capsule symmetry, axisymmetric (Pₗ/P₀)", value: "legendre" },
                    { text: "Capsule symmetry, 3D (all m)", value: "ylm" },
                    { text: "Foil spot uniformity", value: "irradiance" }]
            currentIndex: Math.max(0, ["legendre", "ylm", "irradiance"].indexOf(app.optObjective))
            onActivated: Julia.br_opt_objective(currentValue)
        }
        RowLayout {
            visible: app.optObjective !== "irradiance"
            spacing: 4
            Label { text: "Mode weights"; Layout.preferredWidth: 100 }
            Repeater {
                model: [1, 2, 3, 4]
                delegate: RowLayout {
                    spacing: 2
                    Label { text: "ℓ" + modelData; opacity: 0.7 }
                    TextField {
                        Layout.preferredWidth: 38
                        text: app["optW" + modelData]
                        onEditingFinished: if (text !== app["optW" + modelData])
                                               Julia.br_opt_edit("w" + modelData, text)
                    }
                }
            }
        }
        CheckBox {
            visible: app.optObjective !== "irradiance"
            text: "Hold the mean capsule drive"
            checked: app.optHold
            onToggled: Julia.br_opt_hold(checked)
        }
        RowLayout {
            visible: app.optObjective !== "irradiance" && app.hasCapsule
            spacing: 8
            CheckBox {
                text: "Capsule keep-out"
                Layout.preferredWidth: 170
                checked: app.optCapsule
                onToggled: Julia.br_opt_capsule(checked)
                ToolTip.visible: hovered
                ToolTip.delay: 400
                ToolTip.text: "Keep every beam's worst-case envelope (1 % of peak, with pointing, " +
                              "positioning and focus errors) off the capsule by at least this much"
            }
            TextField {
                Layout.preferredWidth: 90
                enabled: app.optCapsule
                text: app.optCapsuleMin
                onEditingFinished: if (text !== app.optCapsuleMin) Julia.br_opt_edit("capsule_min", text)
            }
            Label { text: "mm clearance"; opacity: 0.7 }
        }
        RowLayout {
            visible: app.optObjective !== "irradiance"
            spacing: 8
            CheckBox {
                text: "LEH keep-out"
                Layout.preferredWidth: 170
                checked: app.optKeepout
                onToggled: Julia.br_opt_keepout(checked)
                ToolTip.visible: hovered
                ToolTip.delay: 400
                ToolTip.text: "Keep every beam's worst-case envelope off the LEH rims by at least this much"
            }
            TextField {
                Layout.preferredWidth: 90
                enabled: app.optKeepout
                text: app.optKeepoutMin
                onEditingFinished: if (text !== app.optKeepoutMin) Julia.br_opt_edit("keepout_min", text)
            }
            Label { text: "mm clearance"; opacity: 0.7 }
        }
        PropRow {
            visible: app.optObjective === "irradiance"
            label: "Spill weight"
            value: app.optSpill
            unit: "power outside the spot"
            onEdited: (text) => Julia.br_opt_edit("spill", text)
        }

        SectionLabel { text: "Limits" }
        PropRow {
            label: "Iterations"
            value: "" + app.optIterations
            unit: "L-BFGS-B iterations"
            onEdited: (text) => Julia.br_opt_edit("iterations", text)
        }
        PropRow {
            label: "Evaluations"
            value: "" + app.optMaxEvaluations
            unit: "of the objective"
            tip: "A line search can take many evaluations in one iteration; this bounds the run time"
            onEdited: (text) => Julia.br_opt_edit("max_evaluations", text)
        }
    }
    ErrorText { text: app.optError }
    RowLayout {
        spacing: 8
        Button {
            text: "Optimize"
            highlighted: true
            enabled: !app.optRunning && app.driveCurrent
            onClicked: Julia.br_opt_start()
        }
        Button { text: "Cancel"; enabled: app.optRunning; onClicked: Julia.br_opt_cancel() }
        Button {
            text: app.optApplied ? "Applied" : "Apply result"
            enabled: app.optHasResult && !app.optApplied && !app.optRunning
            onClicked: Julia.br_opt_accept()
        }
    }
    Note {
        text: !app.driveCurrent ? "Compute the drive first (Compute ▸ Drive)." : app.optProgress
        color: !app.driveCurrent ? "#d08a00" : palette.windowText
        opacity: 1
    }
}
