// The clearance check: component roles, envelope and error set, sampling, run.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    SectionLabel { text: "Components" }
    Repeater {
        model: compModel
        delegate: RoleRow {
            enabled: !app.rvRunning
            label: model.label
            role: model.role
            readonly property string compName: model.name
            roles: [{ text: "beams land on it", value: "target",
                      tip: "A surface the beams are meant to land on: each beam is checked only up to where it lands" },
                    { text: "check clearance", value: "hardware",
                      tip: "Hardware the beams must clear: each beam's worst-case envelope is measured from it" },
                    { text: "ignore", value: "none", tip: "Left out of the check" }]
            onChosen: (value) => Julia.br_rv_role(compName, value)
        }
    }
    SectionLabel { text: "Envelope and errors"; info: "clearance.errors" }
    ColumnLayout {
        spacing: 6
        enabled: !app.rvRunning && !app.tlRunning
        Layout.fillWidth: true
        PropRow { label: "Contour"; value: app.rvPeak; unit: "of each spot's peak"
                  tip: "The envelope is each spot's contour at this fraction of its peak, grown by the f-cone"
                  onEdited: (text) => Julia.br_rv_edit("peak_fraction", text) }
        PropRow { label: "Pointing error"; value: app.rv_pointing; unit: "mm"
                  tip: "Displacement of each aim point in any direction across the beam"
                  onEdited: (text) => Julia.br_rv_edit("pointing", text) }
        PropRow { label: "Target positioning"; value: app.rv_target_offset; unit: "mm"
                  tip: "Displacement of the whole target from its nominal position"
                  onEdited: (text) => Julia.br_rv_edit("target_offset", text) }
        PropRow { label: "Hardware tolerance"; value: app.rv_hardware_tolerance; unit: "mm"
                  tip: "Fabrication and assembly tolerance of the hardware, added to the other errors"
                  onEdited: (text) => Julia.br_rv_edit("hardware_tolerance", text) }
        PropRow { label: "Focus error"; value: app.rv_focus; unit: "± mm"
                  tip: "Shift of best focus along each beam, both ways"
                  onEdited: (text) => Julia.br_rv_edit("focus", text) }
        CheckLine { text: "Far-field wings (1 % halo, 10⁻³ scatter)"; checked: app.rvWings; onToggled: Julia.br_rv_wings(checked)
                    tip: "The NIF far-field wings (nif_wings) also count in the envelope, for the clearance check only" }
        RowLayout {
            spacing: 8
            Cell { text: "Sampling"; Layout.preferredWidth: 150 }
            ComboBox {
                Layout.fillWidth: true
                textRole: "text"; valueRole: "value"
                model: [{ text: "Quick", value: "quick" }, { text: "Standard (4× the points)", value: "standard" }]
                currentIndex: app.rvQuality === "standard" ? 1 : 0
                onActivated: Julia.br_rv_quality(currentValue)
                tip: "Quick: hardware at triangle vertices and 4 sub-triangle points, 8 pointing-error angles, a 31-point spot grid; " +
                     "Standard: 16 sub-points, 16 angles, a 61-point grid"
            }
        }
    }
    ErrorText { text: app.rvError }
    RowLayout {
        spacing: 8
        Button { text: "Check clearance"; highlighted: true; enabled: !app.rvRunning; onClicked: Julia.br_rv_start()
                 tip: "Check every enabled beam in the background; the report is in the Clearance tab" }
        Button { text: "Cancel"; enabled: app.rvRunning; onClicked: Julia.br_rv_cancel()
                 tip: "Stop the check and keep the beams checked so far" }
    }
    ProgressBar {
        Layout.fillWidth: true
        from: 0
        to: Math.max(1, app.rvTotal)
        value: app.rvDone
        visible: app.rvRunning
    }
    Note {
        text: app.rvRunning ? "checked " + app.rvDone + " of " + app.rvTotal + " beams" :
              app.rvStale ? "Out of date. " + app.rvBrief : app.rvStatus
        tip: app.rvStale && !app.rvRunning ? "The beams, the target or the settings changed since this check" : ""
        color: app.rvStale && !app.rvRunning ? "#d08a00" : palette.windowText
        opacity: 1
    }
}
