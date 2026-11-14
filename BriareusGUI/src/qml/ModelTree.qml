// The model tree: target components (with their colors in the 3D view) under their
// parts and assemblies, beams by group and quad, analyses. Click selects (in the 3D view
// too), Ctrl + click adds to the selection, Shift + click a range, double click fits the
// view to it, the arrow expands, a beam row's dots switch its 1ω, 2ω and 3ω light in the
// calculations (top) and in the view (bottom), a component's box (its color, filled by
// its display mode) shows or hides it, right click opens the row's commands. The menu is
// built for the row clicked, so it holds only the commands that apply to it.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "common"

Pane {
    id: tree
    padding: 0
    property alias groupDialog: groupDialog
    property alias rowMenu: menu
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
            Info { topic: "assembly.commands" }
            Item { Layout.fillWidth: app.selectionCount <= 1 }
            Label {
                visible: app.selectionCount > 1
                text: app.selectionCount + " selected"
                opacity: 0.6
                Layout.fillWidth: true
                horizontalAlignment: Text.AlignRight
                HoverHandler { id: countHover }
                Tip { visible: countHover.hovered; text: "Esc clears the selection" }
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
                            tree.openMenu(node.nodeId, node.nodeKind, node.nodeDisplay)
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
                        Tip {
                            visible: boxHover.hovered
                            text: compBox.mixed ? "Partly shown: click to show all"
                                                : compBox.shown ? "Click to hide" : "Click to show"
                        }
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
                        Tip { visible: labelHover.hovered && parent.truncated; text: model.label }
                        HoverHandler { id: labelHover }
                    }
                    Label {
                        text: model.detail
                        opacity: 0.6
                        font.pixelSize: 11
                        elide: Text.ElideRight
                        Layout.fillWidth: true
                        Tip { visible: detailHover.hovered && parent.truncated; text: model.detail }
                        HoverHandler { id: detailHover }
                    }
                }
            }
        }
    }

    // The row menu, built for the row clicked: only the commands that apply to it (a
    // hidden item of a fixed menu still takes its padding in the style's menus).
    // Each entry is { text, run, checkable, checked, enabled }; null is a separator.
    Component {
        id: menuItem
        MenuItem {
            property var run: null
            onTriggered: if (run) run(checked)
        }
    }
    Component { id: menuSeparator; MenuSeparator {} }
    AutoMenu { id: menu }
    function openMenu(id, kind, display) {
        const entries = menuEntries(id, kind, display)
        while (menu.count > 0) {
            const it = menu.takeItem(0)
            if (it) it.destroy()
        }
        let last = null
        for (const e of entries) {
            if (e === null && last === null) continue             // no leading or doubled separator
            // (created under the tree, an item: the menu takes them over)
            menu.addItem(e === null ? menuSeparator.createObject(tree)
                                    : menuItem.createObject(tree, { text: e.text, run: e.run, checkable: !!e.checkable,
                                                                    checked: !!e.checked, enabled: e.enabled !== false }))
            last = e
        }
        if (last === null && menu.count > 0) menu.takeItem(menu.count - 1).destroy()
        if (menu.count > 0) menu.popup()
    }
    function menuEntries(id, kind, display) {
        // [shown, included, quads to split, quads to merge, piece roleable, piece role,
        //  role 0 allowed, role 1 allowed, role 2 allowed] for the rows the menu acts on
        const f = Julia.br_node_menu(id)
        const several = app.selectionCount > 1
        const beamRow = ["beams", "group", "quad", "beam"].indexOf(kind) >= 0
        const treeNode = kind === "part" || kind === "assembly"
        const isComp = display !== ""
        const name = id.substring(id.indexOf(":") + 1)
        const e = []
        if (beamRow) {
            e.push({ text: "Show 3ω light in the 3D view", checkable: true, checked: f[0] === 1,
                     run: (on) => Julia.br_node_show(id, on) })
            e.push({ text: "Include in calculations", checkable: true, checked: f[1] === 1,
                     run: (on) => Julia.br_node_include(id, on) })
            e.push(null)
            e.push({ text: "Fit view", run: () => Julia.br_fit_selection() })
            if (kind !== "beams") {
                e.push({ text: "Look along beam", run: () => Julia.br_look_along_beam() })
                e.push({ text: "Move to a new group…", run: () => groupDialog.open() })
            }
            if (f[2] > 0 || f[3] > 0) e.push(null)
            if (f[2] > 0)
                e.push({ text: f[2] === 1 ? "Split quad into beams" : "Split " + f[2] + " quads into beams",
                         run: () => Julia.br_split_quads(id) })
            if (f[3] > 0)
                e.push({ text: f[3] === 1 ? "Merge beams into quad" : "Merge beams into " + f[3] + " quads",
                         run: () => Julia.br_merge_quads(id) })
        } else if (isComp) {
            const modes = ["opaque", "transparent", "wireframe", "invisible"]
            for (let k = 0; k < 4; ++k)
                e.push({ text: app.displayLabels[k], checkable: true, checked: display === modes[k],
                         run: () => Julia.br_set_display(id, modes[k]) })
            e.push(null)
            e.push({ text: "Include in calculations", checkable: true, checked: f[1] === 1,
                     run: (on) => Julia.br_node_include(id, on) })
            if (f[4] === 1) {
                e.push(null)
                for (let k = 0; k < 3; ++k)
                    e.push({ text: "Role: " + app.pieceRoleLabels[k], checkable: true, checked: f[5] === k,
                             enabled: f[6 + k] === 1, run: () => Julia.br_piece_role(k) })
            }
            e.push(null)
            e.push({ text: "Fit view", run: () => Julia.br_fit_selection() })
        } else if (kind === "frame") {
            e.push({ text: "Draw axes in the 3D view", checkable: true, checked: f[0] === 1,
                     run: (on) => Julia.br_set_triad(id, on) })
            e.push({ text: "Insert frame in it", run: () => Julia.br_add_frame() })
            e.push({ text: "Coordinate transform…", run: () => Julia.br_open_transform() })
        } else if (kind === "frames") {
            e.push({ text: "Insert frame", run: () => Julia.br_add_frame() })
            e.push({ text: "Coordinate transform…", run: () => Julia.br_open_transform() })
        } else if (kind === "target") {
            e.push({ text: "Fit view", run: () => Julia.br_fit_all() })
            e.push({ text: "New assembly", run: () => Julia.br_new_assembly() })
        } else if (kind === "task") {
            e.push({ text: "Run", run: () => Julia.br_run_task(name) })
        }
        if (treeNode) {
            e.push(null)
            e.push({ text: "New assembly of the selection", run: () => Julia.br_new_assembly() })
            e.push({ text: "Move into assembly…", run: () => assemblyDialog.open() })
            if (kind === "assembly" && !several) e.push({ text: "Ungroup", run: () => Julia.br_ungroup_assembly() })
        }
        if ((treeNode || kind === "frame") && !several) {
            if (kind === "frame") e.push(null)
            e.push({ text: "Rename…", run: () => renameDialog.openFor(name) })
        }
        if (treeNode || kind === "window" || kind === "frame") {
            if (kind === "window") e.push(null)
            e.push({ text: several ? "Remove " + app.selectionCount + " selected"
                                   : kind === "assembly" ? "Remove with its contents" : "Remove",
                     run: () => Julia.br_remove_selected() })
        }
        return e
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
                tip: "The members stay where they are in space"
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
