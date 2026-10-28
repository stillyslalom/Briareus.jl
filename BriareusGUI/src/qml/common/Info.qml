// An ⓘ beside a heading or a control: its Help entry (`topic`, see help.js), or `tip`
// when given (an entry with current values), on hover; a click opens the Help screen
// at the topic.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "help.js" as Help

Label {
    id: info
    property string topic: ""
    property string tip: ""
    readonly property string explanation: tip !== "" ? tip : Help.text(topic)
    visible: explanation !== ""
    text: "ⓘ"
    color: palette.link
    opacity: infoHover.hovered ? 1 : 0.7
    font.weight: Font.Normal
    HoverHandler { id: infoHover; cursorShape: Qt.PointingHandCursor }
    TapHandler { onTapped: Julia.br_help(info.topic) }
    ToolTip {
        visible: infoHover.hovered
        delay: 300
        text: info.explanation
        // a paragraph wraps at a readable width
        width: Math.min(implicitWidth, 440)
    }
}
