// The TaLIS rule checks: capture and checked components, the rule tolerance and its
// interpretation, run. The report is in the TaLIS tab.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note {
        text: "The NIF target-configuration rules for each enabled beam: 3ω capture, incidence " +
              "and clearance, 1ω clearance and unconverted light, and the TAS stay-in zone, " +
              "each judged on the worst case over the clearance error set plus the rule tolerance."
    }
    SectionLabel { text: "Components" }
    Repeater {
        model: talisCompModel
        delegate: RowLayout {
            id: compRow
            readonly property string compName: model.name
            readonly property string compRole: model.role
            readonly property string compLabel: model.label
            enabled: !app.tlRunning
            spacing: 8
            Label { text: compRow.compLabel; Layout.preferredWidth: 120; elide: Text.ElideRight }
            ComboBox {
                Layout.fillWidth: true
                textRole: "text"; valueRole: "value"
                model: [{ text: "capture: beams land on it", value: "capture" },
                        { text: "check 3ω and 1ω clearance", value: "hardware" },
                        { text: "ignore", value: "none" }]
                currentIndex: Math.max(0, ["capture", "hardware", "none"].indexOf(compRow.compRole))
                onActivated: Julia.br_tl_role(compRow.compName, currentValue)
            }
        }
    }
    SectionLabel { text: "Rule tolerance" }
    ColumnLayout {
        spacing: 6
        enabled: !app.tlRunning
        PropRow { label: "Tolerance"; value: app.tlTolerance; unit: "degrees (0: none)"
                  onEdited: (text) => Julia.br_tl_edit("tolerance", text) }
        RowLayout {
            spacing: 8
            Label { text: "Distance r from"; Layout.preferredWidth: 150 }
            ComboBox {
                Layout.fillWidth: true
                textRole: "text"; valueRole: "value"
                model: [{ text: "TCC", value: "tcc" },
                        { text: "the aim point", value: "aim" },
                        { text: "the farther of the two", value: "max" },
                        { text: "a point (x, y, z mm)", value: "point" }]
                currentIndex: Math.max(0, ["tcc", "aim", "max", "point"].indexOf(app.tlPivot))
                onActivated: Julia.br_tl_pivot(currentValue)
            }
        }
        PropRow { visible: app.tlPivot === "point"; label: "Pivot point"; value: app.tlPivotPoint
                  fieldWidth: 140; unit: "x, y, z (mm)"
                  onEdited: (text) => Julia.br_tl_edit("pivot_point", text) }
        CheckBox { text: "1ω rules (NIF beams only)"; checked: app.tlOmega1; onToggled: Julia.br_tl_omega1(checked) }
    }
    Note {
        text: "The angular tolerance becomes a lateral allowance r·tan(tolerance) at each hardware " +
              "point. Error set (from the clearance check): pointing " + app.rv_pointing +
              " mm, positioning " + app.rv_target_offset + " mm, hardware " + app.rv_hardware_tolerance +
              " mm, focus ± " + app.rv_focus + " mm; " + app.rvQuality + " sampling."
    }
    ErrorText { text: app.tlError }
    RowLayout {
        spacing: 8
        Button { text: "Run TaLIS checks"; highlighted: true; enabled: !app.tlRunning; onClicked: Julia.br_tl_start() }
        Button { text: "Cancel"; enabled: app.tlRunning; onClicked: Julia.br_tl_cancel() }
        BusyIndicator { running: app.tlRunning; visible: running; Layout.preferredHeight: 28; Layout.preferredWidth: 28 }
    }
    Note {
        text: app.tlRunning ? "checking… " + app.tlElapsed :
              app.tlStale ? "The beams or the target changed since this check. " + app.tlBrief : app.tlStatus
        color: app.tlStale && !app.tlRunning ? "#d08a00" : palette.windowText
        opacity: 1
    }
}
