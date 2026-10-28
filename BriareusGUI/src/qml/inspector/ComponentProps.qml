// Properties of a template component (wall, caps, LEHs, capsule, foil): what it is,
// and how the 3D view draws it.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note { text: app.compInfo; opacity: 0.9 }
    DisplayRow {}
    Note {
        text: "Template components are sized from the target dimensions: select Target to edit " +
              "them. However a component is drawn (invisible included), it takes part in every calculation."
    }
    Button { text: "Target properties"; onClicked: Julia.br_select("target") }
}
