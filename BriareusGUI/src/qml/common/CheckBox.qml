// The style's CheckBox at the window's density (see TextField.qml); `tip` shows on hover.
import QtQuick
import QtQuick.Controls as T

T.CheckBox {
    property string tip: ""
    topPadding: 2
    bottomPadding: 2
    implicitHeight: Math.max(implicitIndicatorHeight, implicitContentHeight) + topPadding + bottomPadding
    T.ToolTip.visible: tip !== "" && hovered
    T.ToolTip.delay: 400
    T.ToolTip.text: tip
}
