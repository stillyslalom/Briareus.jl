// A table cell or a label of fixed width: one line, cut short to its column with the
// whole text in a tooltip; `tip`, when given, shows on hover (after the text when it is
// cut short).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Label {
    property string tip: ""
    elide: Text.ElideRight
    ToolTip.visible: (truncated || tip !== "") && cellHover.hovered
    ToolTip.delay: 400
    ToolTip.text: tip === "" ? text : truncated ? text + ": " + tip : tip
    HoverHandler { id: cellHover }
}
