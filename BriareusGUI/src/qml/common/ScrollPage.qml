// A vertically scrolling column for a panel or inspector page; children go into the
// column, which is as wide as the view.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

ScrollView {
    id: scroll
    default property alias content: column.data
    property int margin: 14
    contentWidth: availableWidth
    contentHeight: column.implicitHeight + 2 * margin
    clip: true
    ColumnLayout {
        id: column
        x: scroll.margin
        y: scroll.margin
        width: scroll.availableWidth - 2 * scroll.margin
        spacing: 8
    }
}
