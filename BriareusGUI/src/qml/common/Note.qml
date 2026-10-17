// Secondary explanatory text that wraps to the panel width.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Label {
    opacity: 0.72
    wrapMode: Text.WordWrap
    Layout.fillWidth: true
    visible: text !== ""
}
