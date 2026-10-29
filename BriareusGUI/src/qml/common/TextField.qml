// The style's TextField at the window's density: every page that imports "common" gets
// this one instead of QtQuick.Controls' (a later import shadows an earlier one).
import QtQuick
import QtQuick.Controls as T

T.TextField {
    topPadding: 3
    bottomPadding: 3
    leftPadding: 6
    rightPadding: 6
    implicitHeight: Math.max(24, contentHeight + topPadding + bottomPadding)
}
