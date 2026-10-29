// The display mode of the selected components in the 3D view: opaque, semi-transparent,
// wireframe or invisible. Drawing only: the calculations see every component included in
// them, however it is drawn.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import "."
import jlqml

RowLayout {
    Label { text: "Display"; Layout.preferredWidth: 110 }
    ComboBox {
        Layout.fillWidth: true
        model: app.displayLabels
        currentIndex: app.compDisplay
        onActivated: (index) => Julia.br_set_display(app.selection, index)
        ToolTip.visible: hovered
        ToolTip.delay: 500
        ToolTip.text: "How the 3D view draws the selected components (drawing only: " +
                      "whether they take part in the calculations is separate)"
    }
}
