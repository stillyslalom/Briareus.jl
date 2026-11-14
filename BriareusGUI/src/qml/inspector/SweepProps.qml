// The placement sweep: what moves, along which frame's axis or through which number of
// its placement steps, the values, the checks at each step, the sampling, run.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    ColumnLayout {
        spacing: 6
        enabled: !app.swRunning
        Layout.fillWidth: true
        SectionLabel { text: "Move" }
        RowLayout {
            Label { text: "Part"; Layout.preferredWidth: 110 }
            ComboBox {
                Layout.fillWidth: true
                model: app.swSubjects
                currentIndex: app.swSubject
                displayText: app.swSubject < 0 ? "(add a part)" : currentText
                onActivated: (index) => Julia.br_sw_subject(index)
                tip: "The part to move; what is placed in it follows"
            }
        }
        RowLayout {
            Label { text: "Vary"; Layout.preferredWidth: 110 }
            ComboBox {
                Layout.fillWidth: true
                model: ["Its place along a frame's axis", "A number of its placement steps"]
                currentIndex: app.swMode
                onActivated: (index) => Julia.br_sw_mode(index)
                tip: "Move along an axis, or vary a placement step"
            }
            Info { topic: "sweep.steps" }
        }
        RowLayout {
            visible: app.swMode === 1
            Label { text: "Step"; Layout.preferredWidth: 110 }
            ComboBox {
                Layout.fillWidth: true
                model: app.swStepChoices
                currentIndex: app.swStepIndex
                displayText: app.swStepChoices.length === 0 ? "(the part has no placement steps)" :
                             app.swStepIndex < 0 ? "(choose a step)" : currentText
                onActivated: (index) => Julia.br_sw_step(index)
                tip: "The placement step whose number the sweep varies"
            }
        }
        RowLayout {
            visible: app.swMode === 1
            Label { text: "Number"; Layout.preferredWidth: 110 }
            ComboBox {
                Layout.fillWidth: true
                model: app.swNumbers
                currentIndex: app.swNumber
                displayText: app.swNumber < 0 ? "(choose a step)" : currentText
                onActivated: (index) => Julia.br_sw_number(index)
                tip: "The number of the step to vary"
            }
        }
        RowLayout {
            visible: app.swMode === 0
            spacing: 4
            Label { text: "Along"; Layout.preferredWidth: 110 }
            ComboBox {
                Layout.preferredWidth: 52
                model: ["x", "y", "z"]
                currentIndex: app.swAxis
                onActivated: (index) => Julia.br_sw_axis(index)
                tip: "The axis the part moves along"
            }
            Label { text: "of" }
            ComboBox {
                Layout.fillWidth: true
                model: app.swFrames
                currentIndex: app.swFrame
                onActivated: (index) => Julia.br_sw_frame(index)
                tip: "The frame whose axis the part moves along"
            }
        }
        PropRow { label: "From"; value: app.sw_start
                  unit: app.swMode === 1 ? app.swUnit : "mm (shift from where it is)"
                  tip: "The first value"
                  onEdited: (text) => Julia.br_sw_edit("start", text) }
        PropRow { label: "To"; value: app.sw_stop; unit: app.swMode === 1 ? app.swUnit : "mm"
                  tip: "The last value"
                  onEdited: (text) => Julia.br_sw_edit("stop", text) }
        PropRow { label: "Increment"; value: app.sw_step
                  unit: (app.swMode === 1 ? app.swUnit : "mm") + (app.swSteps > 0 ? " · " + app.swSteps + " steps" : "")
                  tip: "The change from one step to the next"
                  onEdited: (text) => Julia.br_sw_edit("step", text) }

        SectionLabel { text: "Checks at each step" }
        CheckLine {
            text: "Beam clearance (the clearance check's components, envelope and errors)"
            checked: app.sw_clearance
            onToggled: Julia.br_sw_check("clearance", checked)
            tip: "Check the beam clearance at each step"
        }
        CheckLine {
            text: "Line of sight from the backlighter foil"
            checked: app.sw_los
            enabled: app.hasFoil
            onToggled: Julia.br_sw_check("los", checked)
            tip: app.hasFoil ? "Check whether the target blocks the line of sight"
                             : "Needs a backlighter foil (a part with the foil role)"
        }
        ColumnLayout {
            visible: app.sw_los
            spacing: 6
            Layout.leftMargin: 24
            Layout.fillWidth: true
            RowLayout {
                Label { text: "Through"; Layout.preferredWidth: 86 }
                ComboBox {
                    Layout.fillWidth: true
                    model: app.swObjects
                    currentIndex: app.swObject
                    onActivated: (index) => Julia.br_sw_object(index)
                    tip: "The part the line of sight passes through"
                }
            }
            PropRow { label: "Toward θ"; labelWidth: 86; value: app.sw_los_theta; unit: "° (empty: foil → object)"
                      tip: "Polar angle of the line of sight; empty: from the foil to the object"
                      onEdited: (text) => Julia.br_sw_edit("los_theta", text) }
            PropRow { label: "φ"; labelWidth: 86; value: app.sw_los_phi; unit: "°"
                      tip: "The line of sight's azimuth"
                      onEdited: (text) => Julia.br_sw_edit("los_phi", text) }
        }
        CheckLine {
            text: "Power the foil captures (least fraction per beam)"
            checked: app.sw_capture
            enabled: app.hasFoil
            onToggled: Julia.br_sw_check("capture", checked)
            tip: app.hasFoil ? "The least fraction of a beam's power the active foil absorbs"
                             : "Needs a backlighter foil (a part with the foil role)"
        }
        CheckLine {
            text: "Re-optimize the pointing at each step (Optimization's settings)"
            checked: app.sw_optimize
            onToggled: Julia.br_sw_check("optimize", checked)
            tip: "Re-optimize the pointing at each step (slow)"
        }
        RowLayout {
            spacing: 8
            Label { text: "Sampling"; Layout.preferredWidth: 110 }
            ComboBox {
                Layout.fillWidth: true
                textRole: "text"; valueRole: "value"
                model: [{ text: "Quick", value: "quick" }, { text: "Standard", value: "standard" }]
                currentIndex: app.swQuality === "standard" ? 1 : 0
                onActivated: Julia.br_sw_quality(currentValue)
                tip: "Quick: coarser sampling, faster. Standard: finer"
            }
        }
    }
    ErrorText { text: app.swError }
    RowLayout {
        spacing: 8
        Button { text: "Sweep"; highlighted: true; enabled: !app.swRunning; onClicked: Julia.br_sw_start()
                 tip: "Check each step in the background; the steps are in the Sweep tab" }
        Button { text: "Cancel"; enabled: app.swRunning; onClicked: Julia.br_sw_cancel()
                 tip: "Stop the sweep and keep the steps checked so far" }
    }
    ProgressBar {
        Layout.fillWidth: true
        from: 0
        to: Math.max(1, app.swTotal)
        value: app.swDone
        visible: app.swRunning
    }
    Note {
        text: app.swRunning ? "step " + app.swDone + " of " + app.swTotal :
              app.swStale ? "Out of date. " + app.swBrief : app.swStatus
        tip: app.swStale && !app.swRunning ? "The design or the settings changed since this sweep" : ""
        color: app.swStale && !app.swRunning ? "#d08a00" : palette.windowText
        opacity: 1
    }
}
