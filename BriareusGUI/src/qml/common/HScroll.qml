// A column of content (a table's header and rows) that scrolls sideways when the panel
// is narrower than `minimumContentWidth`, instead of cutting the last columns off. The
// column takes the full height and at least the full width.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Flickable {
    id: hs
    default property alias content: column.data
    property real minimumContentWidth: 0
    property alias spacing: column.spacing
    readonly property bool scrolls: minimumContentWidth > width
    contentWidth: Math.max(width, minimumContentWidth)
    contentHeight: height
    flickableDirection: Flickable.HorizontalFlick
    boundsBehavior: Flickable.StopAtBounds
    clip: true
    ScrollBar.horizontal: ScrollBar { id: bar; policy: hs.scrolls ? ScrollBar.AlwaysOn : ScrollBar.AlwaysOff }
    ColumnLayout {
        id: column
        width: hs.contentWidth
        height: hs.height - (hs.scrolls ? bar.height : 0)
        spacing: 4
    }
}
