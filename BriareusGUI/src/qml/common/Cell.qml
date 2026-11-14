// A table cell or a label of fixed width: one line, cut short to its column with the
// whole text in a tooltip; `tip`, when given, shows on hover (after the text when it is
// cut short).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import "."

Label {
    property string tip: ""
    elide: Text.ElideRight
    Tip {
        text: parent.tip === "" ? parent.text : parent.truncated ? parent.text + ": " + parent.tip : parent.tip
        visible: (parent.truncated || parent.tip !== "") && cellHover.hovered
    }
    HoverHandler { id: cellHover }
}
