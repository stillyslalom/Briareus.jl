// A section heading inside a panel, with an ⓘ for its Help entry when `info` names one
// (`infoTip`: a one-line hint with current values, in place of the entry's own).
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
    Info { topic: parent.info; tip: parent.infoTip; visible: parent.info !== "" && (known || tip !== "") }
    Item { Layout.fillWidth: true }
}
