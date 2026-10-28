// A short status line (a summary, a state, an outcome) that wraps to the panel width;
// explanations go in help.js, shown by an ⓘ (Info.qml).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

Label {
    opacity: 0.72
    wrapMode: Text.WordWrap
    Layout.fillWidth: true
    visible: text !== ""
}
