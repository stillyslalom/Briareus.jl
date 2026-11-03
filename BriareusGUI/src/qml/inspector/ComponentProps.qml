// Properties of a component of the target (a hohlraum part's piece, a CAD part's): what it is,
// its role (a hohlraum part's piece), and how the 3D view draws it.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note { text: app.compInfo; opacity: 0.9 }
    RowLayout {
        visible: app.pieceRoleable
        Label { text: "Role"; Layout.preferredWidth: 110 }
        ComboBox {
            Layout.fillWidth: true
            model: app.pieceRoleLabels
            currentIndex: app.pieceRole
            onActivated: (index) => Julia.br_piece_role(index)
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "What the analyses take this piece as, in place of its part's default (a cap can be the observed face)"
        }
        Info { topic: "component.role" }
    }
    DisplayRow {}
    Button { text: "Target properties"; onClicked: Julia.br_select("target") }
}
