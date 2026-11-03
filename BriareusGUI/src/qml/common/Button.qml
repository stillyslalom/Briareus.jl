// The style's Button at the window's density (see TextField.qml); `tip` shows on hover.
import QtQuick
import QtQuick.Controls as T

T.Button {
    property string tip: ""
    topPadding: 3
    bottomPadding: 3
    leftPadding: 10
    rightPadding: 10
    implicitHeight: Math.max(26, implicitContentHeight + topPadding + bottomPadding)
    T.ToolTip.visible: tip !== "" && hovered
    T.ToolTip.delay: 400
    T.ToolTip.text: tip
}
