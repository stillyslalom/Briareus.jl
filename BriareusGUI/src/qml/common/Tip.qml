// A tooltip for a quick hint: a line or two, wrapped at `maxWidth` and centered over the
// item it is declared in. The style's own ToolTip centers itself from its unwrapped
// width, which puts a long text far from its control; this one centers from the width
// it is drawn at. Details belong in Help (an ⓘ opens it), not here.
import QtQuick
import QtQuick.Controls

ToolTip {
    id: tip
    property int maxWidth: 320
    delay: 500
    margins: 6
    width: Math.min(implicitWidth, maxWidth)
    x: parent ? Math.round((parent.width - width) / 2) : 0
    y: -height - 4
}
