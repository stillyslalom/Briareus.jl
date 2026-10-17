// Properties of a template component (wall, caps, LEHs, capsule, foil): what it is,
// and whether it is shown.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note { text: app.compInfo; opacity: 0.9 }
    CheckBox {
        text: "Show in the 3D view"
        checked: !app.compHidden
        onToggled: Julia.br_tree_on(app.selection, checked)
    }
    Note {
        text: "Template components are sized from the target dimensions: select Target to edit " +
              "them. Hidden components still take part in every calculation."
    }
    Button { text: "Target properties"; onClicked: Julia.br_select("target") }
}
