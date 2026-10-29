// The style's Button at the window's density (see TextField.qml).
import QtQuick
import QtQuick.Controls as T

T.Button {
    topPadding: 3
    bottomPadding: 3
    leftPadding: 10
    rightPadding: 10
    implicitHeight: Math.max(26, implicitContentHeight + topPadding + bottomPadding)
}
