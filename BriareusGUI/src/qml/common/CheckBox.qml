// The style's CheckBox at the window's density (see TextField.qml).
import QtQuick
import QtQuick.Controls as T

T.CheckBox {
    topPadding: 2
    bottomPadding: 2
    implicitHeight: Math.max(implicitIndicatorHeight, implicitContentHeight) + topPadding + bottomPadding
}
