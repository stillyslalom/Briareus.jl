// An error message (hidden when empty).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Label {
    color: "#c42b1c"
    wrapMode: Text.WordWrap
    Layout.fillWidth: true
    visible: text !== ""
}
