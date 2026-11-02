// The model tree: target components (with their colors in the 3D view) under their
// parts and assemblies, beams by group and quad, analyses. Click selects (in the 3D view too), Ctrl + click adds beams
// to the selection, double click fits the view to it, the arrow expands, a beam row's
// dots switch its 1ω, 2ω and 3ω light in the calculations (top) and in the view
// (bottom), a component's box (its color, filled by its display mode) shows or hides
// it, right click opens the node's actions (a component's display mode among them).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "common"

Pane {
    id: tree
    padding: 0
    property alias groupDialog: groupDialog
    function openGroupDialog() { groupDialog.open() }

    ColumnLayout {
        anchors.fill: parent
        spacing: 0
        // a title row as high as the Properties panel's, with the size of a multiple selection
        RowLayout {
            Layout.fillWidth: true
            Layout.leftMargin: 8
            Layout.rightMargin: 8
            Layout.topMargin: 6
            Label { text: "Model"; font.pixelSize: 15; font.weight: Font.DemiBold }
            Label {
                visible: app.selectionCount > 1
                text: app.selectionCount + " selected"
                opacity: 0.6
                Layout.fillWidth: true
                horizontalAlignment: Text.AlignRight
                HoverHandler { id: countHover }
                ToolTip.visible: countHover.hovered
                ToolTip.delay: 400
                ToolTip.text: "Ctrl + click adds or removes; Shift + click selects a range"
            }
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
                height: 22
                readonly property string nodeId: model.id
                readonly property string nodeKind: model.kind
                readonly property string nodeDisplay: model.drawmode
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
                            axesItem.checked = f[0]
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
                    anchors.leftMargin: 6 + 14 * model.depth
                    anchors.rightMargin: 6
                    spacing: 4
                    Label {
                        text: model.expandable ? (model.expanded ? "▾" : "▸") : ""
                        Layout.preferredWidth: 12
                        horizontalAlignment: Text.AlignHCenter
                        opacity: 0.7
                        TapHandler { onTapped: Julia.br_toggle_expand(node.nodeId) }
                    }
                    BeamDots {
                        visible: model.matrix !== ""
                        matrix: model.matrix
                        nodeId: node.nodeId
                        Layout.preferredWidth: implicitWidth
                    }
                    // a component: one box in its color, checked while drawn; the fill is its
                    // display mode (solid opaque, pale semi-transparent, an outline as
                    // wireframe, empty hidden). A click shows or hides it; its menu sets the mode.
                    // A part in pieces or an assembly: its components' box together ("mixed": some
                    // shown, some hidden; a click shows them all).
                    Rectangle {
                        id: compBox
                        visible: model.drawmode !== ""
                        readonly property color base: model.color === "" ? palette.mid : model.color
                        readonly property bool mixed: model.drawmode === "mixed"
                        readonly property bool shown: model.drawmode !== "invisible" && !mixed
                        readonly property bool dark: 0.299 * base.r + 0.587 * base.g + 0.114 * base.b < 0.6
                        implicitWidth: 14
                        implicitHeight: 14
                        radius: 3
                        color: model.drawmode === "opaque" ? base
                             : model.drawmode === "transparent" || mixed ? Qt.rgba(base.r, base.g, base.b, 0.35) : "transparent"
                        border.color: shown || mixed ? base : Qt.rgba(base.r, base.g, base.b, 0.6)
                        border.width: model.drawmode === "wireframe" ? 2 : 1.2
                        Label {
                            anchors.centerIn: parent
                            visible: compBox.shown || compBox.mixed
                            text: compBox.mixed ? "–" : "✓"
                            font.pixelSize: 11
                            font.bold: true
                            color: model.drawmode === "opaque" ? (compBox.dark ? "white" : "black") : palette.text
                        }
                        HoverHandler { id: boxHover; cursorShape: Qt.PointingHandCursor }
                        ToolTip.visible: boxHover.hovered
                        ToolTip.delay: 400
                        ToolTip.text: (compBox.mixed ? "Partly drawn: click to show all" :
                                       (compBox.shown ? "Drawn " : "Hidden") +
                                       ({ opaque: "opaque", transparent: "semi-transparent", wireframe: "as wireframe",
                                          invisible: "" })[model.drawmode] +
                                       ": click to " + (compBox.shown ? "hide" : "show")) +
                                      ", right-click for the display mode (drawing only: it stays in the calculations)"
                        TapHandler { onTapped: Julia.br_tree_on(node.nodeId, !compBox.shown) }
                    }
                    StatusDot { visible: model.status !== ""; status: model.status }
                    // a beam set's color (a component's is its box)
                    Rectangle {
                        visible: model.color !== "" && model.drawmode === ""
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
                    readonly property bool showable: ["comp", "part", "window", "assembly", "group", "quad", "beam", "beams"].indexOf(node.nodeKind) >= 0
                    // what the model tree's commands act on: parts and assemblies
                    readonly property bool treeNode: node.nodeKind === "part" || node.nodeKind === "assembly"
                    readonly property bool isComp: node.nodeDisplay !== ""
                    MenuItem {
                        id: showItem
                        text: "Show 3ω light in the 3D view"
                        checkable: true
                        visible: menu.showable && !menu.isComp
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_node_show(node.nodeId, checked)
                    }
                    MenuItem {
                        text: app.displayLabels[0]
                        checkable: true
                        checked: node.nodeDisplay === "opaque"
                        visible: menu.isComp; height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_set_display(node.nodeId, "opaque")
                    }
                    MenuItem {
                        text: app.displayLabels[1]
                        checkable: true
                        checked: node.nodeDisplay === "transparent"
                        visible: menu.isComp; height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_set_display(node.nodeId, "transparent")
                    }
                    MenuItem {
                        text: app.displayLabels[2]
                        checkable: true
                        checked: node.nodeDisplay === "wireframe"
                        visible: menu.isComp; height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_set_display(node.nodeId, "wireframe")
                    }
                    MenuItem {
                        text: app.displayLabels[3]
                        checkable: true
                        checked: node.nodeDisplay === "invisible"
                        visible: menu.isComp; height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_set_display(node.nodeId, "invisible")
                    }
                    MenuItem {
                        id: includeItem
                        text: "Include in calculations"
                        checkable: true
                        visible: menu.showable
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_node_include(node.nodeId, checked)
                    }
                    MenuItem {
                        text: "Draw axes in the 3D view"
                        checkable: true
                        id: axesItem
                        visible: node.nodeKind === "frame"
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_set_triad(node.nodeId, checked)
                    }
                    MenuItem {
                        text: node.nodeKind === "frame" ? "Add frame in it" : "Add frame"
                        visible: node.nodeKind === "frame" || node.nodeKind === "frames"
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_add_frame()
                    }
                    MenuItem {
                        text: "Coordinate transform…"
                        visible: node.nodeKind === "frame" || node.nodeKind === "frames"
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_open_transform()
                    }
                    MenuSeparator { visible: menu.showable; height: visible ? implicitHeight : 0 }
                    MenuItem {
                        text: "Fit view"
                        visible: node.nodeKind !== "frame" && node.nodeKind !== "frames"
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_fit_selection()
                    }
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
                    MenuSeparator { visible: menu.treeNode || node.nodeKind === "target"; height: visible ? implicitHeight : 0 }
                    MenuItem {
                        text: node.nodeKind === "target" ? "New assembly" : "New assembly of the selection"
                        visible: menu.treeNode || node.nodeKind === "target"
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_new_assembly()
                    }
                    MenuItem {
                        text: "Move into assembly…"
                        visible: menu.treeNode
                        height: visible ? implicitHeight : 0
                        onTriggered: assemblyDialog.open()
                    }
                    MenuItem {
                        text: "Ungroup"
                        visible: node.nodeKind === "assembly"
                        height: visible ? implicitHeight : 0
                        onTriggered: Julia.br_ungroup_assembly()
                    }
                    MenuItem {
                        text: "Rename…"
                        visible: menu.treeNode || node.nodeKind === "frame"
                        height: visible ? implicitHeight : 0
                        onTriggered: renameDialog.openFor(node.nodeId.substring(node.nodeId.indexOf(":") + 1))
                    }
                    MenuItem {
                        text: node.nodeKind === "assembly" ? "Remove with its contents" : "Remove"
                        visible: menu.treeNode || node.nodeKind === "window" || node.nodeKind === "frame"
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

    // Move into assembly: the selected parts and assemblies, kept where they are
    Dialog {
        id: assemblyDialog
        title: "Move into assembly"
        anchors.centerIn: Overlay.overlay
        modal: true
        standardButtons: Dialog.Ok | Dialog.Cancel
        onOpened: assemblyChoice.currentIndex = 0
        onAccepted: Julia.br_move_into_assembly(assemblyChoice.currentIndex)
        RowLayout {
            spacing: 8
            Label { text: "Into" }
            ComboBox {
                id: assemblyChoice
                model: app.assemblyChoices
                Layout.preferredWidth: Math.max(200, implicitWidth)
                ToolTip.visible: hovered
                ToolTip.delay: 400
                ToolTip.text: "The members stay where they are; the top level is in no assembly"
            }
        }
    }

    // Rename a part, assembly or frame: what is placed in it, aims in it or is a member follows
    Dialog {
        id: renameDialog
        title: "Rename"
        anchors.centerIn: Overlay.overlay
        modal: true
        standardButtons: Dialog.Ok | Dialog.Cancel
        function openFor(name) { newName.text = name; open() }
        onOpened: { newName.selectAll(); newName.forceActiveFocus() }
        onAccepted: if (newName.text.trim() !== "") Julia.br_rename_node(newName.text)
        ColumnLayout {
            spacing: 8
            Label { text: "New name (letters, digits and _)" }
            TextField {
                id: newName
                Layout.fillWidth: true
                Layout.preferredWidth: 220
                onAccepted: renameDialog.accept()
            }
        }
    }
}
