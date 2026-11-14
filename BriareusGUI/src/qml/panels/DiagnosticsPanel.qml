// Diagnostics: the synthetic diagnostics a drive history records, by category. A toggle
// chooses the category (x-ray power, x-ray imagers, radiography); the row beside it adds
// instruments of that category and shows the recorded ones as chips (a chip's menu
// chooses its configuration); below, the category's view on one canvas: the power
// diagnostics' chart, or one of the imagers' images.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQml.Models
import Makie
import jlqml
import "../common"

Item {
    id: panel
    readonly property real fitHeight: column.implicitHeight + 2 * margin
    readonly property int margin: 6
    readonly property int viewMinimum: 175
    readonly property string category: app.diagCategory
    function redraw() { diagCanvas.update() }
    function openPicker() { picker.open(); return picker }

    ColumnLayout {
        id: column
        anchors.fill: parent
        anchors.margins: panel.margin
        spacing: 4

        // the categories, then what is recorded of the chosen one (wraps on a narrow panel
        // or with many instruments); the definitions' directories at the right
        RowLayout {
            Layout.fillWidth: true
            spacing: 6
            Flow {
                Layout.fillWidth: true
                Layout.alignment: Qt.AlignTop
                spacing: 6
                Repeater {
                    model: app.diagKinds
                    delegate: Button {
                        required property string modelData
                        required property int index
                        text: app.diagKindLabels[index] + (app.diagCounts[index] > 0 ? " (" + app.diagCounts[index] + ")" : "")
                        implicitWidth: implicitContentWidth + leftPadding + rightPadding
                        flat: panel.category !== modelData
                        highlighted: panel.category === modelData
                        onClicked: Julia.br_diag_category(modelData)
                        tip: app.diagCounts[index] + " recorded"
                    }
                }
                ToolSeparator { height: 26 }
                Button {
                    id: pickerButton
                    text: "Add ▾"
                    onClicked: picker.opened ? picker.close() : picker.open()
                    tip: "Record an instrument in the next drive history"
                }
                Repeater {
                    model: diagModel
                    delegate: Rectangle {
                        id: chip
                        required property string id
                        required property string name
                        required property string kind
                        required property bool ticked
                        required property bool usable
                        required property bool changed
                        required property bool foreign
                        required property string tip
                        required property string config
                        required property string configIds
                        required property string configLabels
                        readonly property var ids: configIds.split("|")
                        readonly property var labels: configLabels.split("|")
                        visible: ticked && kind === panel.category
                        height: 26
                        width: chipRow.implicitWidth + 10
                        radius: 13
                        color: "transparent"
                        border.color: palette.mid
                        function openMenu() { chipMenu.popup() }
                        TapHandler { acceptedButtons: Qt.RightButton; onTapped: chip.openMenu() }
                        TapHandler { acceptedButtons: Qt.LeftButton; onDoubleTapped: chip.openMenu() }
                        Row {
                            id: chipRow
                            anchors.verticalCenter: parent.verticalCenter
                            x: 9
                            spacing: 2
                            Label {
                                text: chip.name + (chip.config !== "" ? " · " + chip.config : "") + (chip.changed ? " *" : "")
                                anchors.verticalCenter: parent.verticalCenter
                                opacity: chip.usable && !chip.foreign ? 1 : 0.5
                                HoverHandler { id: chipHover }
                                Tip {
                                    visible: chipHover.hovered
                                    text: (chip.foreign ? "Another facility's instrument\n" : "") + chip.tip +
                                          (chip.ids.length > 1 ? "\nRight-click for its configurations" : "")
                                }
                            }
                            ToolButton {
                                text: "✕"
                                implicitHeight: 24
                                implicitWidth: 24
                                font.pixelSize: 10
                                onClicked: Julia.br_time_diagnostic(chip.id, false)
                                Tip { visible: parent.hovered; text: "Stop recording" }
                            }
                        }
                        // the instrument's configurations (when it has several), then its commands
                        AutoMenu {
                            id: chipMenu
                            Instantiator {
                                model: chip.ids.length > 1 ? chip.ids.length : 0
                                delegate: MenuItem {
                                    required property int index
                                    text: chip.labels[index]
                                    checkable: true
                                    checked: chip.ids[index] === chip.id
                                    onTriggered: Julia.br_diag_config(chip.id, chip.ids[index])
                                }
                                onObjectAdded: (index, object) => chipMenu.insertItem(index, object)
                                onObjectRemoved: (index, object) => chipMenu.removeItem(object)
                            }
                            MenuItem {
                                text: "Use the directories' definition"
                                enabled: chip.changed
                                onTriggered: Julia.br_diag_use_path(chip.id)
                            }
                            MenuItem { text: "Stop recording"; onTriggered: Julia.br_time_diagnostic(chip.id, false) }
                        MenuSeparator {}
                        MenuItem { text: "About configurations…"; onTriggered: Julia.br_help("diagnostics.config") }
                        }
                    }
                }
            }
            Button { text: "Diagnostic directories…"; flat: true; onClicked: win.openDiagnostics()
                     Layout.alignment: Qt.AlignTop
                     tip: "Where definitions are read from, and the problems found" }
            Button { text: "Refresh"; flat: true; onClicked: Julia.br_diag_reload(); tip: "Read every definition again"
                     Layout.alignment: Qt.AlignTop }
            Info { topic: "diagnostics.record"; Layout.alignment: Qt.AlignTop; Layout.topMargin: 4 }
        }

        // the category's instruments of this facility, a check box each (opens under its button)
        Popup {
            id: picker
            parent: pickerButton
            x: 0
            y: pickerButton.height + 2
            margins: 6
            padding: 6
            closePolicy: Popup.CloseOnEscape | Popup.CloseOnPressOutsideParent
            width: pickerList.implicitWidth + leftPadding + rightPadding + 14
            height: Math.min(pickerList.implicitHeight + topPadding + bottomPadding,
                             (Overlay.overlay ? Overlay.overlay.height : 800) - 12)
            contentItem: ScrollView {
                clip: true
                contentWidth: pickerList.implicitWidth
                Column {
                    id: pickerList
                    spacing: 0
                    Label {
                        visible: pickerRows.shown === 0
                        text: "No instrument of this kind for this facility"
                        opacity: 0.7
                        padding: 4
                    }
                    Repeater {
                        id: pickerRows
                        model: diagModel
                        // how many rows the category shows
                        readonly property int shown: {
                            let n = 0
                            for (let i = 0; i < count; ++i) {
                                const it = itemAt(i)
                                if (it && it.visible) n += 1
                            }
                            return n
                        }
                        delegate: CheckBox {
                            visible: model.first && model.kind === panel.category
                            height: visible ? implicitHeight : 0
                            text: model.name
                            checked: model.ticked
                            enabled: model.usable || model.ticked
                            onToggled: Julia.br_diag_record(model.name, checked)
                            // (why an instrument cannot be recorded)
                            tip: model.usable ? "" : model.tip
                        }
                    }
                }
            }
        }

        // the chart's or the image's title
        RowLayout {
            Layout.fillWidth: true
            spacing: 8
            Label {
                text: panel.category === "xray_power" ? app.diagTitle
                    : panel.category === "xray_imager" ? "Images of the recorded imagers" : "Radiography"
                font.weight: Font.DemiBold
                elide: Text.ElideRight
                Layout.maximumWidth: implicitWidth
                Layout.fillWidth: true
            }
            Info { topic: panel.category === "xray_power" ? "diagnostics.chart" : "diagnostics.images"
                   visible: panel.category !== "radiography" }
            Item { Layout.fillWidth: true }
        }

        // the category's view on one canvas (it lives as long as the window, and shows the
        // chart's figure or the image's): x-ray power, the chart and its key; x-ray
        // imagers, the history's images in a list and the one chosen; radiography, a note
        RowLayout {
            Layout.fillWidth: true
            Layout.fillHeight: true
            Layout.minimumHeight: panel.viewMinimum
            spacing: 10
            ColumnLayout {
                visible: panel.category === "xray_imager"
                Layout.fillHeight: true
                Layout.preferredWidth: 300
                Layout.maximumWidth: 300
                spacing: 4
                Note {
                    visible: imageList.count === 0
                    text: app.diagCounts[1] > 0 ? "No images yet: run the drive history (Time tab)"
                                                : "No imager recorded: Add ▾, then run the drive history (Time tab)"
                }
                ListView {
                    id: imageList
                    Layout.fillWidth: true
                    Layout.fillHeight: true
                    clip: true
                    model: diagImageModel
                    boundsBehavior: Flickable.StopAtBounds
                    ScrollBar.vertical: ScrollBar {}
                    delegate: ItemDelegate {
                        required property string id
                        required property string label
                        required property int frame
                        required property int index
                        width: imageList.width
                        text: label
                        topPadding: 3
                        bottomPadding: 3
                        highlighted: index === app.diagImageIndex
                        onClicked: Julia.br_diag_pick_image(id, frame)
                        onDoubleClicked: Julia.br_diag_image(id, frame)
                        Tip {
                            visible: parent.hovered
                            text: frame === 0 ? "Detected radiance at the time cursor (W/cm²/sr)"
                                              : "Detected fluence of the recording (J/cm²/sr)"
                        }
                    }
                }
                Button {
                    text: "Open in the main view"
                    visible: imageList.count > 0
                    enabled: app.diagImageIndex >= 0
                    onClicked: {
                        const it = imageList.itemAtIndex(app.diagImageIndex)
                        if (it) Julia.br_diag_image(it.id, it.frame)
                    }
                    tip: "Show the image in place of the 3D view, where it zooms and pans"
                }
            }
            Item {
                Layout.fillWidth: true
                Layout.fillHeight: true
                MakieArea {
                    id: diagCanvas
                    objectName: "diagCanvas"
                    anchors.fill: parent
                    scene: app.diagCanvas
                    visible: panel.category === "xray_power" ? app.diagEmpty === ""
                           : panel.category === "xray_imager" && imageList.count > 0
                }
                Note {
                    anchors.top: parent.top
                    anchors.left: parent.left
                    width: Math.min(parent.width - 20, 520)
                    text: panel.category === "xray_power" ? app.diagEmpty
                        : panel.category === "radiography"
                          ? "Radiography definitions can be added, but a drive history records nothing for them yet." : ""
                }
            }
            ColumnLayout {
                visible: panel.category === "xray_power" && app.diagEmpty === ""
                Layout.alignment: Qt.AlignTop
                Layout.maximumWidth: 260
                spacing: 3
                Repeater {
                    model: diagLegendModel
                    delegate: RowLayout {
                        id: entry
                        Layout.fillWidth: true
                        spacing: 6
                        readonly property string lineStyle: model.style
                        readonly property string lineColor: model.color
                        Row {
                            spacing: 2
                            Layout.preferredWidth: 22
                            Repeater {
                                model: entry.lineStyle === "solid" ? 1 : entry.lineStyle === "dash" ? 3 : 5
                                delegate: Rectangle {
                                    width: entry.lineStyle === "solid" ? 22 : entry.lineStyle === "dash" ? 6 : 2.4
                                    height: entry.lineStyle === "solid" ? 3 : 2
                                    color: entry.lineColor
                                }
                            }
                        }
                        Cell { text: model.label; font.pixelSize: 12; Layout.fillWidth: true }
                    }
                }
            }
        }
    }
}
