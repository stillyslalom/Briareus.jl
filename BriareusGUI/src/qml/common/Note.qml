// A short status line, such as a summary, a state or an outcome, that wraps to the panel
// width. Explanations go in help.js, shown by an ⓘ (Info.qml), or in `tip`, shown on
// hover.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import "."

Label {
    property string tip: ""
    opacity: 0.72
    wrapMode: Text.WordWrap
    Layout.fillWidth: true
    visible: text !== ""
    Tip { text: parent.tip; visible: parent.tip !== "" && noteHover.hovered }
    HoverHandler { id: noteHover }
}
