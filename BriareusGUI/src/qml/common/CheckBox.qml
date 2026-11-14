// The style's CheckBox at the window's density (see TextField.qml); `tip` shows on hover.
import QtQuick
import QtQuick.Controls as T
import "."

T.CheckBox {
    property string tip: ""
    topPadding: 2
    bottomPadding: 2
    implicitHeight: Math.max(implicitIndicatorHeight, implicitContentHeight) + topPadding + bottomPadding
    Tip { text: parent.tip; visible: parent.tip !== "" && parent.hovered }
}
