// The placement sweep: what moves and along which frame's axis, the shifts, the checks
// at each step, the sampling, run.
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
            }
        }
        RowLayout {
            spacing: 4
            Label { text: "Along"; Layout.preferredWidth: 110 }
            ComboBox {
                Layout.preferredWidth: 52
                model: ["x", "y", "z"]
                currentIndex: app.swAxis
                onActivated: (index) => Julia.br_sw_axis(index)
            }
            Label { text: "of" }
            ComboBox {
                Layout.fillWidth: true
                model: app.swFrames
                currentIndex: app.swFrame
                onActivated: (index) => Julia.br_sw_frame(index)
            }
        }
        PropRow { label: "From"; value: app.sw_start; unit: "mm (shift from where it is)"
                  onEdited: (text) => Julia.br_sw_edit("start", text) }
        PropRow { label: "To"; value: app.sw_stop; unit: "mm"
                  onEdited: (text) => Julia.br_sw_edit("stop", text) }
        PropRow { label: "Step"; value: app.sw_step; unit: app.swSteps > 0 ? "mm · " + app.swSteps + " steps" : "mm"
                  onEdited: (text) => Julia.br_sw_edit("step", text) }

        SectionLabel { text: "Checks at each step" }
        CheckLine {
            text: "Beam clearance (the clearance check's components, envelope and errors)"
            checked: app.sw_clearance
            onToggled: Julia.br_sw_check("clearance", checked)
        }
        CheckLine {
            text: "Line of sight from the backlighter foil"
            checked: app.sw_los
            enabled: app.hasFoil
            onToggled: Julia.br_sw_check("los", checked)
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
                }
            }
            PropRow { label: "Toward θ"; labelWidth: 86; value: app.sw_los_theta; unit: "° (empty: foil → object)"
                      onEdited: (text) => Julia.br_sw_edit("los_theta", text) }
            PropRow { label: "φ"; labelWidth: 86; value: app.sw_los_phi; unit: "°"
                      onEdited: (text) => Julia.br_sw_edit("los_phi", text) }
        }
        CheckLine {
            text: "Power the foil captures (least fraction per beam)"
            checked: app.sw_capture
            enabled: app.hasFoil
            onToggled: Julia.br_sw_check("capture", checked)
        }
        CheckLine {
            text: "Re-optimize the pointing at each step (Optimization's settings)"
            checked: app.sw_optimize
            onToggled: Julia.br_sw_check("optimize", checked)
            tip: "Slow: an optimization per step. Aims free in their frame (u, v, w) warm-start " +
                 "from the previous step, as the frame moves with the part."
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
            }
        }
    }
    ErrorText { text: app.swError }
    RowLayout {
        spacing: 8
        Button { text: "Sweep"; highlighted: true; enabled: !app.swRunning; onClicked: Julia.br_sw_start() }
        Button { text: "Cancel"; enabled: app.swRunning; onClicked: Julia.br_sw_cancel() }
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
              app.swStale ? "The design or the settings changed since this sweep. " + app.swBrief : app.swStatus
        color: app.swStale && !app.swRunning ? "#d08a00" : palette.windowText
        opacity: 1
    }
}
