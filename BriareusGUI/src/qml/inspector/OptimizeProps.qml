// The pointing optimization: free fields and bounds, objective, keep-outs, limits, run.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
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
                model: app.hasQuads ? [{ text: "cone", value: "group" }, { text: "quad", value: "quad" },
                                       { text: "beam", value: "beam" }]
                                    : [{ text: "cone", value: "group" }, { text: "beam", value: "beam" }]
                currentIndex: Math.max(0, (app.hasQuads ? ["group", "quad", "beam"] : ["group", "beam"]).indexOf(app.optBy))
                onActivated: Julia.br_opt_set_by(currentValue)
                tip: "One set of free fields per cone, per quad or per beam: the beams of one share their values"
            }
        }

        SectionLabel { text: "Free fields and bounds" }
        Repeater {
            model: [{ key: "z", label: "Aim z (chamber)", unit: "± mm", tip: "The aim's z, free within ± this many mm of the beam table's value" },
                    { key: "r", label: "Aim r (chamber)", unit: "± mm", tip: "The aim's radial offset, free within ± this many mm of the beam table's value" },
                    { key: "phi", label: "Tangential offset", unit: "± mm", tip: "The aim's tangential offset, free within ± this many mm of the beam table's value" },
                    { key: "u", label: "Aim u (frame x)", unit: "± mm in the aim frame", tip: "The aim along its frame's x (on a foil's frame, across its face)" },
                    { key: "v", label: "Aim v (frame y)", unit: "± mm in the aim frame", tip: "The aim along its frame's y (on a foil's frame, across its face)" },
                    { key: "w", label: "Aim w (frame z)", unit: "± mm in the aim frame", tip: "The aim along its frame's z (on a foil's frame, along its normal)" },
                    { key: "defocus", label: "Defocus", unit: "± mm", tip: "The defocus, free within ± this many mm of the beam table's value" },
                    { key: "power", label: "Power", unit: "± fraction of the start", tip: "The power, free within ± this fraction of the beam table's" }]
            delegate: PropRow {
                check: true
                label: modelData.label
                labelWidth: 170
                checked: app["optFree_" + modelData.key]
                onToggled: (on) => Julia.br_opt_free(modelData.key, on)
                value: app["optDelta_" + modelData.key]
                unit: modelData.unit
                tip: modelData.tip
                onEdited: (text) => Julia.br_opt_edit(modelData.key, text)
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
            tip: "Axisymmetric: the weighted (Pₗ/P₀)² of the capsule flux; 3D: the weighted rotation-invariant mode " +
                 "power over all m (a dropped quad); foil: the spot's RMS nonuniformity plus the spill"
        }
        RowLayout {
            visible: app.optObjective !== "irradiance"
            spacing: 4
            Cell { text: "Mode weights"; Layout.preferredWidth: 100; tip: "The weight of each degree ℓ in the objective (0: left out)" }
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
                        tip: "The weight of mode ℓ = " + modelData + " (0: left out)"
                    }
                }
            }
        }
        CheckLine {
            visible: app.optObjective !== "irradiance"
            text: "Hold the mean capsule drive"
            checked: app.optHold
            onToggled: Julia.br_opt_hold(checked)
            tip: "Keep the capsule's mean drive at its start value"
        }
        PropRow {
            visible: app.optObjective !== "irradiance" && app.observedKind === "capsule"
            check: true
            label: "Capsule keep-out"
            labelWidth: 170
            checked: app.optCapsule
            onToggled: (on) => Julia.br_opt_capsule(on)
            fieldEnabled: app.optCapsule
            value: app.optCapsuleMin
            unit: "mm clearance"
            tip: "Keep every beam's worst-case envelope (1 % of peak, with pointing, " +
                 "positioning and focus errors) off the capsule by at least this much"
            onEdited: (text) => Julia.br_opt_edit("capsule_min", text)
        }
        PropRow {
            visible: app.optObjective !== "irradiance"
            check: true
            label: "LEH keep-out"
            labelWidth: 170
            checked: app.optKeepout
            onToggled: (on) => Julia.br_opt_keepout(on)
            fieldEnabled: app.optKeepout
            value: app.optKeepoutMin
            unit: "mm clearance"
            tip: "Keep every beam's worst-case envelope off the LEH rims by at least this much"
            onEdited: (text) => Julia.br_opt_edit("keepout_min", text)
        }
        PropRow {
            visible: app.optObjective === "irradiance"
            label: "Spill weight"
            value: app.optSpill
            unit: "power outside the spot"
            tip: "The weight on the squared fraction of the power that falls outside the spot radius"
            onEdited: (text) => Julia.br_opt_edit("spill", text)
        }

        SectionLabel { text: "Limits" }
        PropRow {
            label: "Iterations"
            value: "" + app.optIterations
            unit: "L-BFGS-B iterations"
            tip: "The most iterations of the optimizer"
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
            tip: "Run the optimization in the background; the result is in the Optimization tab"
        }
        Button { text: "Cancel"; enabled: app.optRunning; onClicked: Julia.br_opt_cancel()
                 tip: "Stop the optimization and keep the best design so far" }
        Button {
            text: app.optApplied ? "Applied" : "Apply result"
            enabled: app.optHasResult && !app.optApplied && !app.optRunning
            onClicked: Julia.br_opt_accept()
            tip: "Write the optimized pointing into the beam table"
        }
    }
    Note {
        text: !app.driveCurrent ? "Compute the drive first (Compute ▸ Drive)." : app.optProgress
        color: !app.driveCurrent ? "#d08a00" : palette.windowText
        opacity: 1
    }
}
