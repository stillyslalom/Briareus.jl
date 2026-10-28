// Properties of a diagnostic window in the hohlraum wall.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note { text: app.compInfo; opacity: 0.9 }
    Repeater {
        model: windowModel
        delegate: ColumnLayout {
            id: w
            visible: model.number === app.windowNumber
            spacing: 6
            Layout.fillWidth: true
            readonly property int rowNumber: model.number
            readonly property bool isOpen: model.open
            readonly property bool isRound: model.round
            PropRow { label: "Azimuth φ"; value: model.phi; unit: "° from +x"
                      onEdited: (text) => Julia.br_edit_window(w.rowNumber, "phi", text) }
            PropRow { label: "Center z"; value: model.z; unit: "mm"
                      onEdited: (text) => Julia.br_edit_window(w.rowNumber, "z", text) }
            CheckBox {
                text: "Round"
                checked: w.isRound
                onToggled: Julia.br_window_round(w.rowNumber, checked)
                ToolTip.visible: hovered
                ToolTip.text: "A disk on the unrolled wall (otherwise a rectangle)"
            }
            PropRow { label: w.isRound ? "Diameter" : "Width"; value: model.width
                      unit: w.isRound ? "mm" : "mm of arc"
                      onEdited: (text) => Julia.br_edit_window(w.rowNumber, "width", text) }
            PropRow { visible: !w.isRound; label: "Height"; value: model.height; unit: "mm along z"
                      onEdited: (text) => Julia.br_edit_window(w.rowNumber, "height", text) }
            CheckBox {
                text: "Open (a hole: laser light and x-rays leave through it)"
                checked: w.isOpen
                onToggled: Julia.br_window_open(w.rowNumber, checked)
            }
            PropRow { label: "Albedo"; value: model.albedo; unit: "closed window"; enabled: !w.isOpen
                      onEdited: (text) => Julia.br_edit_window(w.rowNumber, "albedo", text) }
            PropRow { label: "Conversion"; value: model.conversion; unit: "laser → x-ray"; enabled: !w.isOpen
                      onEdited: (text) => Julia.br_edit_window(w.rowNumber, "conversion", text) }
        }
    }
    DisplayRow {}
    ErrorText { text: app.targetError }
    Button { text: "Remove window"; onClicked: Julia.br_remove_selected() }
}
