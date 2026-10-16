import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Dialogs
import jlqml

StepPage {
    title: "Review"
    guidance: "Worst-case clearance of each enabled beam from the hardware, over pointing, " +
              "target-positioning, hardware and focus errors that add linearly. The envelope is " +
              "each spot's contour at the given fraction of its peak, grown by the f-cone; with " +
              "wings, the far field (a 1 % halo and a 10⁻³ scatter floor) also counts toward the " +
              "power and intensity on hardware. Only hardware before a beam's landing on the " +
              "target surfaces is checked. Status thresholds are engineering defaults " +
              "(fail: 10¹² W/cm² or 10⁻³ of the beam power; review: 10⁹ W/cm² or 10⁻⁵)."

    Label { text: "Components"; font.weight: Font.DemiBold }
    ColumnLayout {
        spacing: 2
        enabled: !app.rvRunning
        Repeater {
            model: compModel
            delegate: RowLayout {
                id: compRow
                // the ComboBox's own `model` hides the row's
                readonly property string compName: model.name
                readonly property string compRole: model.role
                spacing: 16
                Label { text: compRow.compName; Layout.preferredWidth: 120 }
                ComboBox {
                    Layout.preferredWidth: 220
                    textRole: "text"; valueRole: "value"
                    model: [{ text: "beams land on it", value: "target" },
                            { text: "check clearance", value: "hardware" },
                            { text: "ignore", value: "none" }]
                    currentIndex: Math.max(0, ["target", "hardware", "none"].indexOf(compRow.compRole))
                    onActivated: Julia.br_rv_role(compRow.compName, currentValue)
                }
            }
        }
    }

    Label { text: "Envelope and errors"; font.weight: Font.DemiBold; Layout.topMargin: 8 }
    GridLayout {
        columns: 3
        columnSpacing: 12
        rowSpacing: 4
        enabled: !app.rvRunning
        Repeater {
            model: [{ key: "peak_fraction", app: "rvPeak", label: "Contour", unit: "of each spot's peak" },
                    { key: "pointing", app: "rv_pointing", label: "Pointing error", unit: "mm" },
                    { key: "target_offset", app: "rv_target_offset", label: "Target positioning", unit: "mm" },
                    { key: "hardware_tolerance", app: "rv_hardware_tolerance", label: "Hardware tolerance", unit: "mm" },
                    { key: "focus", app: "rv_focus", label: "Focus error", unit: "± mm" }]
            delegate: RowLayout {
                Layout.columnSpan: 3
                spacing: 12
                Label { text: modelData.label; Layout.preferredWidth: 160 }
                TextField {
                    Layout.preferredWidth: 90
                    text: app[modelData.app]
                    onEditingFinished: if (text !== app[modelData.app]) Julia.br_rv_edit(modelData.key, text)
                }
                Label { text: modelData.unit; opacity: 0.75 }
            }
        }
        Label { text: "Far-field wings" }
        Switch { checked: app.rvWings; onToggled: Julia.br_rv_wings(checked) }
        Label { text: "" }
        Label { text: "Sampling" }
        ComboBox {
            Layout.preferredWidth: 160
            textRole: "text"; valueRole: "value"
            model: [{ text: "Quick", value: "quick" }, { text: "Standard", value: "standard" }]
            currentIndex: app.rvQuality === "standard" ? 1 : 0
            onActivated: Julia.br_rv_quality(currentValue)
        }
        Label { text: "standard: 4× the points"; opacity: 0.75 }
    }
    Label {
        text: app.rvError
        visible: text !== ""
        color: "#c42b1c"
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }

    RowLayout {
        spacing: 8
        Layout.topMargin: 8
        Button {
            text: "Check clearance"
            highlighted: true
            enabled: !app.rvRunning
            onClicked: Julia.br_rv_start()
        }
        Button { text: "Cancel"; enabled: app.rvRunning; onClicked: Julia.br_rv_cancel() }
        Button {
            text: "Export CSV…"
            enabled: app.rvHasReport && !app.rvRunning
            onClicked: exportDialog.open()
        }
    }
    ProgressBar {
        Layout.fillWidth: true
        from: 0
        to: Math.max(1, app.rvTotal)
        value: app.rvDone
        visible: app.rvRunning
    }
    Label {
        text: app.rvRunning ? "checked " + app.rvDone + " of " + app.rvTotal + " beams" :
              app.rvStale ? "The beams or the target changed since this check. " + app.rvBrief : app.rvStatus
        visible: text !== ""
        color: app.rvStale && !app.rvRunning ? "#d08a00" : palette.windowText
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }

    function statusColor(s) {
        if (s === "fail") return "#c42b1c"
        if (s === "review") return "#d08a00"
        return "#2e9e44"
    }

    Label { text: "Closest beam per component"; font.weight: Font.DemiBold; visible: app.rvHasReport; Layout.topMargin: 8 }
    Repeater {
        model: worstModel
        delegate: Label {
            text: model.component + ": " + model.beam + " · " +
                  (model.worst === "–" ? "no hardware point before the landing"
                                       : "clearance " + model.worst + " mm (nominal " + model.nominal + " mm)") +
                  " · " + model.power + " TW · " + model.status
            color: statusColor(model.status)
            wrapMode: Text.WordWrap
            Layout.fillWidth: true
        }
    }

    Label { text: "Every beam and component"; font.weight: Font.DemiBold; visible: app.rvHasReport; Layout.topMargin: 8 }
    Label {
        text: "Worst-case and nominal clearance in mm (– : no hardware point before the landing), " +
              "worst-case power on the component in TW, peak intensity in W/cm²."
        visible: app.rvHasReport
        opacity: 0.75
        wrapMode: Text.WordWrap
        Layout.fillWidth: true
    }
    RowLayout {
        visible: app.rvHasReport
        spacing: 4
        Repeater {
            model: [{ t: "Beam", w: 56 }, { t: "Component", w: 88 }, { t: "Status", w: 52 },
                    { t: "Worst", w: 60 }, { t: "Nominal", w: 60 },
                    { t: "Power", w: 64 }, { t: "Peak", w: 72 }]
            Label { text: modelData.t; Layout.preferredWidth: modelData.w; font.weight: Font.DemiBold }
        }
    }
    ListView {
        visible: app.rvHasReport
        Layout.fillWidth: true
        Layout.preferredHeight: Math.min(contentHeight, 420)
        clip: true
        model: clearModel
        ScrollBar.vertical: ScrollBar {}
        delegate: RowLayout {
            spacing: 4
            Label { text: model.beam; Layout.preferredWidth: 56 }
            Label { text: model.component; Layout.preferredWidth: 88; elide: Text.ElideRight }
            Label { text: model.status; Layout.preferredWidth: 52; color: statusColor(model.status) }
            Label { text: model.worst; Layout.preferredWidth: 60 }
            Label { text: model.nominal; Layout.preferredWidth: 60; opacity: 0.75 }
            Label { text: model.power; Layout.preferredWidth: 64 }
            Label { text: model.intensity; Layout.preferredWidth: 72 }
        }
    }

    FileDialog {
        id: exportDialog
        title: "Export the clearance report"
        fileMode: FileDialog.SaveFile
        defaultSuffix: "csv"
        nameFilters: ["Clearance report (*.csv)"]
        onAccepted: Julia.br_rv_export(selectedFile.toString())
    }
}
