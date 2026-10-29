// The style's RadioButton at the window's density (see TextField.qml).
import QtQuick
import QtQuick.Controls as T

T.RadioButton {
    topPadding: 2
    bottomPadding: 2
    implicitHeight: Math.max(implicitIndicatorHeight, implicitContentHeight) + topPadding + bottomPadding
}
