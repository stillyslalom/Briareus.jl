// Move along [axis] of [frame] by [Δ] mm [−] [+]: the selected part, assembly or frame
// moved along an axis of any frame, each click one edit (undone as one); the placement
// triad's arrows step by the same Δ. `what` names the moved thing in the tooltips ("the
// part", "the assembly, with its members"); `triad` shows the triad's ⓘ.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "."

ColumnLayout {
    id: move
    property string what: "the part"
    property int labelWidth: 110
    property bool triad: true
    spacing: 6
    Layout.fillWidth: true
    RowLayout {
        spacing: 4
        Layout.fillWidth: true
        Label { text: "Move along"; Layout.preferredWidth: move.labelWidth }
        ComboBox {
            Layout.preferredWidth: 52
            model: ["x", "y", "z"]
            currentIndex: app.moveAxis
            onActivated: (index) => Julia.br_move_axis(index)
            tip: "The axis to move " + move.what + " along"
        }
        Label { text: "of" }
        ComboBox {
            Layout.fillWidth: true
            Layout.minimumWidth: 90
            model: app.moveFrames
            currentIndex: app.moveFrame
            onActivated: (index) => Julia.br_move_frame(index)
            tip: "The frame whose axis " + move.what + " moves along: TCC, a part's own (its axis is z), a piece's or a user frame"
        }
    }
    RowLayout {
        spacing: 4
        Layout.fillWidth: true
        Label { text: "by"; Layout.preferredWidth: move.labelWidth; horizontalAlignment: Text.AlignRight }
        TextField {
            Layout.preferredWidth: 64
            text: app.moveStep
            selectByMouse: true
            onEditingFinished: if (text !== app.moveStep) Julia.br_move_step(text)
            tip: "The step of each click, and of the placement triad's arrows (mm)"
        }
        Label { text: "mm"; opacity: 0.7 }
        Button { text: "−"; Layout.preferredWidth: 34; onClicked: Julia.br_nudge_part(-1)
                 tip: "Move " + move.what + " back along the axis (one edit)" }
        Button { text: "+"; Layout.preferredWidth: 34; onClicked: Julia.br_nudge_part(1)
                 tip: "Move " + move.what + " forward along the axis (one edit)" }
        Info { topic: "part.triad"; visible: move.triad }
    }
}
