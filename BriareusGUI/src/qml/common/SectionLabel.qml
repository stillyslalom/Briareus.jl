// A section heading inside a panel, with an ⓘ for its Help entry when `info` names one
// (`infoTip`: the explanation with current values instead).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

RowLayout {
    property alias text: heading.text
    property string info: ""
    property string infoTip: ""
    spacing: 6
    Layout.topMargin: 6
    Label { id: heading; font.weight: Font.DemiBold }
    Info { visible: parent.info !== "" && explanation !== ""; topic: parent.info; tip: parent.infoTip }
    Item { Layout.fillWidth: true }
}
