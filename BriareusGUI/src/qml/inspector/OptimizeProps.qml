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
                tip: "The beams that share one set of free values"
            }
        }

        SectionLabel { text: "Free fields and bounds" }
        Repeater {
            model: [{ key: "z", label: "Aim z (chamber)", unit: "± mm", tip: "The aim's z, free within ± this many mm of the beam table's value" },
                    { key: "r", label: "Aim r (chamber)", unit: "± mm", tip: "The aim's radial offset, free within ± this many mm" },
                    { key: "phi", label: "Tangential offset", unit: "± mm", tip: "The aim's tangential offset, free within ± this many mm" },
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

        SectionLabel { text: "Objective"; info: "optimize.weights" }
        ComboBox {
            Layout.fillWidth: true
            textRole: "text"; valueRole: "value"
            model: [{ text: "Capsule symmetry, axisymmetric (Pₗ/P₀)", value: "legendre" },
                    { text: "Capsule symmetry, 3D (all m)", value: "ylm" },
                    { text: "Observed face uniformity", value: "face" },
                    { text: "Foil spot uniformity", value: "irradiance" }]
            currentIndex: Math.max(0, ["legendre", "ylm", "face", "irradiance"].indexOf(app.optObjective))
            onActivated: Julia.br_opt_objective(currentValue)
            tip: "What to minimize: capsule modes, or the face's or foil's nonuniformity"
        }
        RowLayout {
            visible: app.optObjective === "legendre" || app.optObjective === "ylm"
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
            text: app.optObjective === "face" ? "Hold the mean face drive" : "Hold the mean capsule drive"
            checked: app.optHold
            onToggled: Julia.br_opt_hold(checked)
            tip: app.optObjective === "face" ? "Keep the face's area-mean x-ray flux at its start value"
                                             : "Keep the capsule's mean drive at its start value"
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
            tip: "Keep every beam's worst-case envelope this far off the capsule"
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
            tip: "Keep every beam's worst-case envelope this far off the LEH rims"
            onEdited: (text) => Julia.br_opt_edit("keepout_min", text)
        }
        PropRow {
            visible: app.optObjective === "face"
            label: "Uniformity weight"
            value: app.optUniformity
            unit: "× (RMS/start RMS)²"
            tip: "Weight on the face's nonuniformity (1: the shape comes first)"
            onEdited: (text) => Julia.br_opt_edit("uniformity", text)
        }
        PropRow {
            visible: app.optObjective === "irradiance"
            label: "Spill weight"
            value: app.optSpill
            unit: "power outside the spot"
            tip: "Weight on the power that misses the spot"
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
            tip: "Bounds the run time"
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
            tip: "Runs in the background; results in the Optimization tab"
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
