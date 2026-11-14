// The style's Button at the window's density (see TextField.qml); `tip` shows on hover.
import QtQuick
import QtQuick.Controls as T
import "."

T.Button {
    property string tip: ""
    topPadding: 3
    bottomPadding: 3
    leftPadding: 10
    rightPadding: 10
    implicitHeight: Math.max(26, implicitContentHeight + topPadding + bottomPadding)
    Tip { text: parent.tip; visible: parent.tip !== "" && parent.hovered }
}
