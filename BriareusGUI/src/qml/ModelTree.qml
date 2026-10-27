// The model tree: target components (with their colors in the 3D view), beams by
// group and quad, analyses. Click selects (in the 3D view too), Ctrl + click adds beams
// to the selection, double click fits the view to it, the arrow expands, the box
// switches a beam set on or off or shows and hides a component, right click opens the
// node's actions.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "common"

Pane {
    id: tree
    padding: 0

    ColumnLayout {
        anchors.fill: parent
        spacing: 0
        Pane {
            Layout.fillWidth: true
            padding: 12
            bottomPadding: 6
            Label { text: "MODEL"; font.pixelSize: 11; font.letterSpacing: 1; opacity: 0.6 }
        }
        Note {
            visible: app.selectionCount > 1
            text: app.selectionCount + " selected (Ctrl + click adds or removes)"
            Layout.leftMargin: 12
        }
        ListView {
            id: list
            Layout.fillWidth: true
            Layout.fillHeight: true
            clip: true
            model: treeModel
            boundsBehavior: Flickable.StopAtBounds
            ScrollBar.vertical: ScrollBar {}
            delegate: Rectangle {
                id: node
                width: list.width
                height: 28
                readonly property string nodeId: model.id
                readonly property string nodeKind: model.kind
                readonly property color hl: palette.highlight
                color: model.selected ? Qt.rgba(hl.r, hl.g, hl.b, 0.28)
                     : hover.hovered ? Qt.rgba(hl.r, hl.g, hl.b, 0.08) : "transparent"
                HoverHandler { id: hover }
                MouseArea {
                    anchors.fill: parent
                    acceptedButtons: Qt.LeftButton | Qt.RightButton
                    onClicked: (mouse) => {
                        if (mouse.button === Qt.RightButton) {
                            if (!model.selected) Julia.br_select(node.nodeId)
                            const f = Julia.br_node_flags(node.nodeId)
                            showItem.checked = f[0]
                            includeItem.checked = f[1]
                            menu.popup()
                        } else if (mouse.modifiers & Qt.ShiftModifier) {
                            Julia.br_select_range(node.nodeId)
                        } else if (mouse.modifiers & Qt.ControlModifier) {
                            Julia.br_select_add(node.nodeId)
                        } else {
                            Julia.br_select(node.nodeId)
                        }
                    }
                    onDoubleClicked: { Julia.br_select(node.nodeId); Julia.br_fit_selection() }
                }
                RowLayout {
                    anchors.fill: parent
                    anchors.leftMargin: 8 + 16 * model.depth
                    anchors.rightMargin: 8
                    spacing: 4
                    Label {
                        text: model.expandable ? (model.expanded ? "▾" : "▸") : ""
                        Layout.preferredWidth: 14
                        horizontalAlignment: Text.AlignHCenter
                        opacity: 0.7
                        TapHandler { onTapped: Julia.br_toggle_expand(node.nodeId) }
                    }
                    CheckBox {
                        visible: model.on >= 0
                        checked: model.on !== 0
                        opacity: model.on === 2 ? 0.55 : 1
                        padding: 0
                        Layout.preferredHeight: 22
                        Layout.preferredWidth: 22
                        onToggled: Julia.br_tree_on(node.nodeId, checked)
                        ToolTip.visible: hovered
                        ToolTip.delay: 500
                        ToolTip.text: node.nodeKind === "comp" || node.nodeKind === "part" || node.nodeKind === "window"
                                      ? "Show or hide in the 3D view" : "Switch these beams on or off"
                    }
                    StatusDot { visible: model.status !== ""; status: model.status }
                    Rectangle {
                        visible: model.color !== ""
                        implicitWidth: 12
                        implicitHeight: 12
                        radius: 2
                        color: model.color === "" ? "transparent" : model.color
                        border.color: Qt.rgba(0, 0, 0, 0.25)
                    }
                    Label {
                        text: model.label
                        font.weight: model.depth === 0 ? Font.DemiBold : Font.Normal
                        elide: Text.ElideRight
                        Layout.maximumWidth: list.width * (model.depth === 0 ? 0.7 : 0.55)   // a heading's detail is a count
                        ToolTip.visible: labelHover.hovered && truncated
                        ToolTip.delay: 400
                        ToolTip.text: model.label
                        HoverHandler { id: labelHover }
                    }
                    Label {
                        text: model.detail
                        opacity: 0.6
                        font.pixelSize: 11
                        elide: Text.ElideRight
                        Layout.fillWidth: true
                        ToolTip.visible: detailHover.hovered && truncated
                        ToolTip.delay: 400
                        ToolTip.text: model.detail
                        HoverHandler { id: detailHover }
                    }
                }
                AutoMenu {
                    id: menu
                    readonly property bool showable: ["comp", "part", "window", "group", "quad", "beam", "beams"].indexOf(node.nodeKind) >= 0
                    MenuItem {
                        id: showItem
                        text: "Show in the 3D view"
                        checkable: true
                        visible: menu.showable
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_node_show(node.nodeId, checked)
                    }
                    MenuItem {
                        id: includeItem
                        text: "Include in calculations"
                        checkable: true
                        visible: menu.showable
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_node_include(node.nodeId, checked)
                    }
                    MenuSeparator { visible: menu.showable; height: visible ? implicitHeight : 0 }
                    MenuItem { text: "Fit view"; onTriggered: Julia.br_fit_selection() }
                    MenuItem {
                        text: "Look along beam"
                        visible: ["group", "quad", "beam"].indexOf(node.nodeKind) >= 0
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_look_along_beam()
                    }
                    MenuItem {
                        text: "Move to a new group…"
                        visible: ["group", "quad", "beam"].indexOf(node.nodeKind) >= 0
                        height: visible ? implicitHeight : 0
                        onTriggered: groupDialog.open()
                    }
                    MenuItem {
                        text: "Split quad into beams"
                        visible: app.bnCanSplit && ["quad", "beam"].indexOf(node.nodeKind) >= 0
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_split_quad()
                    }
                    MenuItem {
                        text: "Merge quad"
                        visible: app.bnCanMerge && ["quad", "beam"].indexOf(node.nodeKind) >= 0
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_merge_quad()
                    }
                    MenuItem {
                        text: "Remove"
                        visible: node.nodeKind === "part" || node.nodeKind === "window"
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_remove_selected()
                    }
                    MenuItem {
                        text: "Run"
                        visible: node.nodeKind === "task"
                        height: visible ? implicitHeight : 0
                        onTriggered: {
                            const k = node.nodeId.substring(5)
                            if (k === "drive") Julia.br_build()
                            else if (k === "optimize") Julia.br_opt_start()
                            else if (k === "clearance") Julia.br_rv_start()
                            else Julia.br_tl_start()
                        }
                    }
                }
            }
        }
    }

    Dialog {
        id: groupDialog
        title: "Move the selected beams to a group"
        anchors.centerIn: Overlay.overlay
        modal: true
        standardButtons: Dialog.Ok | Dialog.Cancel
        onOpened: { groupName.text = ""; groupName.forceActiveFocus() }
        onAccepted: if (groupName.text.trim() !== "") Julia.br_set_cluster(groupName.text)
        ColumnLayout {
            spacing: 8
            Label { text: "Group name (an existing group, or a new one such as \"backlighter\")" }
            TextField {
                id: groupName
                Layout.fillWidth: true
                placeholderText: "backlighter"
                onAccepted: groupDialog.accept()
            }
        }
    }
}
