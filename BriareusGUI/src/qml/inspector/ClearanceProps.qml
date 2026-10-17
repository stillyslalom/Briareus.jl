// The clearance check: component roles, envelope and error set, sampling, run.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note {
        text: "Worst-case clearance of each enabled beam from the hardware over pointing, " +
              "positioning, hardware and focus errors that add linearly. The envelope is each " +
              "spot's contour at the given fraction of its peak grown by the f-cone; with wings, " +
              "the far field also counts. Only hardware before a beam lands is checked. Status " +
              "thresholds are engineering defaults (fail: 10¹² W/cm² or 10⁻³ of the beam power; " +
              "review: 10⁹ W/cm² or 10⁻⁵)."
    }
    SectionLabel { text: "Components" }
    Repeater {
        model: compModel
        delegate: RowLayout {
            id: compRow
            readonly property string compName: model.name
            readonly property string compRole: model.role
            enabled: !app.rvRunning
            spacing: 8
            Label { text: compRow.compName; Layout.preferredWidth: 120; elide: Text.ElideRight }
            ComboBox {
                Layout.fillWidth: true
                textRole: "text"; valueRole: "value"
                model: [{ text: "beams land on it", value: "target" },
                        { text: "check clearance", value: "hardware" },
                        { text: "ignore", value: "none" }]
                currentIndex: Math.max(0, ["target", "hardware", "none"].indexOf(compRow.compRole))
                onActivated: Julia.br_rv_role(compRow.compName, currentValue)
            }
        }
    }
    SectionLabel { text: "Envelope and errors" }
    ColumnLayout {
        spacing: 6
        enabled: !app.rvRunning && !app.tlRunning
        PropRow { label: "Contour"; value: app.rvPeak; unit: "of each spot's peak"
                  onEdited: (text) => Julia.br_rv_edit("peak_fraction", text) }
        PropRow { label: "Pointing error"; value: app.rv_pointing; unit: "mm"
                  onEdited: (text) => Julia.br_rv_edit("pointing", text) }
        PropRow { label: "Target positioning"; value: app.rv_target_offset; unit: "mm"
                  onEdited: (text) => Julia.br_rv_edit("target_offset", text) }
        PropRow { label: "Hardware tolerance"; value: app.rv_hardware_tolerance; unit: "mm"
                  onEdited: (text) => Julia.br_rv_edit("hardware_tolerance", text) }
        PropRow { label: "Focus error"; value: app.rv_focus; unit: "± mm"
                  onEdited: (text) => Julia.br_rv_edit("focus", text) }
        CheckBox { text: "Far-field wings (1 % halo, 10⁻³ scatter)"; checked: app.rvWings; onToggled: Julia.br_rv_wings(checked) }
        RowLayout {
            spacing: 8
            Label { text: "Sampling"; Layout.preferredWidth: 150 }
            ComboBox {
                Layout.fillWidth: true
                textRole: "text"; valueRole: "value"
                model: [{ text: "Quick", value: "quick" }, { text: "Standard (4× the points)", value: "standard" }]
                currentIndex: app.rvQuality === "standard" ? 1 : 0
                onActivated: Julia.br_rv_quality(currentValue)
            }
        }
    }
    Note { text: "The errors and the sampling also apply to the TaLIS checks." }
    ErrorText { text: app.rvError }
    RowLayout {
        spacing: 8
        Button { text: "Check clearance"; highlighted: true; enabled: !app.rvRunning; onClicked: Julia.br_rv_start() }
        Button { text: "Cancel"; enabled: app.rvRunning; onClicked: Julia.br_rv_cancel() }
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
              app.rvStale ? "The beams or the target changed since this check. " + app.rvBrief : app.rvStatus
        color: app.rvStale && !app.rvRunning ? "#d08a00" : palette.windowText
        opacity: 1
    }
}
