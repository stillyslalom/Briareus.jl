// A table cell: one line, cut short to its column with the whole text in a tooltip.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Label {
    elide: Text.ElideRight
    ToolTip.visible: truncated && cellHover.hovered
    ToolTip.delay: 400
    ToolTip.text: text
    HoverHandler { id: cellHover }
}
