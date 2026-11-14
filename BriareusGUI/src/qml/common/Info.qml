// An ⓘ beside a heading or a control. Hovering shows a one-line hint (`tip` when given,
// else the Help topic's own short `tip`, else its title); a click opens the Help window
// at `topic` (see help.js), where the explanation is.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "."
import "help.js" as Help

Label {
    id: info
    property string topic: ""
    property string tip: ""
    readonly property bool known: Help.has(topic)
    readonly property string hint: tip !== "" ? tip : Help.hint(topic)
    visible: known || tip !== ""
    text: "ⓘ"
    color: palette.link
    opacity: infoHover.hovered ? 1 : 0.7
    font.weight: Font.Normal
    HoverHandler { id: infoHover; cursorShape: Qt.PointingHandCursor }
    TapHandler { onTapped: if (info.known) Julia.br_help(info.topic) }
    Tip {
        visible: infoHover.hovered
        delay: 300
        text: info.hint + (info.known ? "\nClick for Help" : "")
    }
}
