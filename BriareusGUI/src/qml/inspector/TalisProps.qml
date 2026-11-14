// The TaLIS rule checks: capture and checked components, the rule tolerance and its
// interpretation, run. The report is in the TaLIS tab.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    SectionLabel { text: "Components" }
    Repeater {
        model: talisCompModel
        delegate: RoleRow {
            enabled: !app.tlRunning
            label: model.label
            role: model.role
            readonly property string compName: model.name
            roles: [{ text: "capture: beams land on it", value: "capture",
                      tip: "The beams land on it: capture and incidence rules" },
                    { text: "check 3ω and 1ω clearance", value: "hardware",
                      tip: "Hardware checked for 3ω and 1ω clearance" },
                    { text: "diagnostic: 3ω, and 10–25 mm from 1ω", value: "diagnostic",
                      tip: "Diagnostic hardware: 3ω clearance, 10–25 mm from 1ω light" },
                    { text: "ignore", value: "none", tip: "Left out of the checks" }]
            onChosen: (value) => Julia.br_tl_role(compName, value)
        }
    }
    SectionLabel {
        text: "Rule tolerance"
        info: "talis.tolerance"
        infoTip: "Uses the clearance check's errors and " + app.rvQuality + " sampling"
    }
    ColumnLayout {
        spacing: 6
        enabled: !app.tlRunning
        Layout.fillWidth: true
        PropRow { label: "Tolerance"; value: app.tlTolerance; unit: "degrees (0: none)"
                  tip: "Rotation of the component off nominal, in the worst direction"
                  onEdited: (text) => Julia.br_tl_edit("tolerance", text) }
        RowLayout {
            spacing: 8
            Cell { text: "Distance r from"; Layout.preferredWidth: 150 }
            ComboBox {
                Layout.fillWidth: true
                textRole: "text"; valueRole: "value"
                model: [{ text: "TCC", value: "tcc" },
                        { text: "the aim point", value: "aim" },
                        { text: "the farther of the two", value: "max" },
                        { text: "a point (x, y, z mm)", value: "point" }]
                currentIndex: Math.max(0, ["tcc", "aim", "max", "point"].indexOf(app.tlPivot))
                onActivated: Julia.br_tl_pivot(currentValue)
                tip: "The point the tolerance rotates the component about"
            }
        }
        PropRow { visible: app.tlPivot === "point"; label: "Pivot point"; value: app.tlPivotPoint
                  fieldWidth: 140; unit: "x, y, z (mm)"
                  tip: "In the chamber frame (TCC)"
                  onEdited: (text) => Julia.br_tl_edit("pivot_point", text) }
        CheckLine { text: "1ω rules (NIF beams only)"; checked: app.tlOmega1; onToggled: Julia.br_tl_omega1(checked)
                    tip: "Also check the 1ω rules" }
    }
    ErrorText { text: app.tlError }
    RowLayout {
        spacing: 8
        Button { text: "Run TaLIS checks"; highlighted: true; enabled: !app.tlRunning; onClicked: Julia.br_tl_start()
                 tip: "Runs in the background; results in the TaLIS tab" }
        Button { text: "Cancel"; enabled: app.tlRunning; onClicked: Julia.br_tl_cancel(); tip: "Stop the checks" }
        BusyIndicator { running: app.tlRunning; visible: running; Layout.preferredHeight: 28; Layout.preferredWidth: 28 }
    }
    Note {
        text: app.tlRunning ? "checking… " + app.tlElapsed :
              app.tlStale ? "Out of date. " + app.tlBrief : app.tlStatus
        tip: app.tlStale && !app.tlRunning ? "The beams, the target or the settings changed since this check" : ""
        color: app.tlStale && !app.tlRunning ? "#d08a00" : palette.windowText
        opacity: 1
    }
}
