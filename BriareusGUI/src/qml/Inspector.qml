// The Properties panel: a title for the selection and its page (inspector/*.qml).
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "common"

Pane {
    id: inspector
    padding: 0

    readonly property string kind: app.selKind
    readonly property string title: {
        if (kind === "target") return "Target"
        if (kind === "comp" || kind === "part" || kind === "window") return app.selTitle
        if (kind === "beams") return "Beams"
        if (kind === "group" || kind === "quad" || kind === "beam") return app.selTitle
        if (kind === "task") {
            if (app.selKey === "drive") return "Drive"
            if (app.selKey === "optimize") return "Pointing optimization"
            if (app.selKey === "clearance") return "Beam clearance"
            if (app.selKey === "sweep") return "Placement sweep"
            if (app.selKey === "talis") return "TaLIS rules"
        }
        return "Properties"
    }
    readonly property string page: {
        if (kind === "target") return "inspector/TargetProps.qml"
        if (kind === "comp") return "inspector/ComponentProps.qml"
        if (kind === "window") return "inspector/WindowProps.qml"
        if (kind === "part") return "inspector/PartProps.qml"
        if (kind === "beams") return "inspector/BeamsProps.qml"
        if (kind === "group" || kind === "quad" || kind === "beam") return "inspector/BeamProps.qml"
        if (kind === "task") {
            if (app.selKey === "drive") return "inspector/DriveProps.qml"
            if (app.selKey === "optimize") return "inspector/OptimizeProps.qml"
            if (app.selKey === "clearance") return "inspector/ClearanceProps.qml"
            if (app.selKey === "sweep") return "inspector/SweepProps.qml"
            if (app.selKey === "talis") return "inspector/TalisProps.qml"
        }
        return ""
    }

    // the Help topic that introduces the page (its ⓘ beside the title)
    readonly property string intro: {
        if (kind === "target") return "target." + app.targetKind
        if (kind === "comp") return "component.intro"
        if (kind === "part") return "part.placement"
        if (kind === "task") return app.selKey + ".intro"
        if (kind === "") return "start.select"
        return ""
    }

    ColumnLayout {
        anchors.fill: parent
        spacing: 0
        // the selection's name heads the page (the panel needs no caption of its own)
        RowLayout {
            spacing: 6
            Layout.fillWidth: true
            Layout.leftMargin: 8
            Layout.rightMargin: 8
            Layout.topMargin: 6
            Label {
                text: inspector.title
                font.pixelSize: 15
                font.weight: Font.DemiBold
                elide: Text.ElideRight
                Layout.fillWidth: true
                Layout.maximumWidth: implicitWidth
                ToolTip.visible: truncated && titleHover.hovered
                ToolTip.delay: 400
                ToolTip.text: inspector.title
                HoverHandler { id: titleHover }
            }
            Info { topic: inspector.intro }
            Item { Layout.fillWidth: true }
        }
        ScrollPage {
            Layout.fillWidth: true
            Layout.fillHeight: true
            Loader {
                id: loader
                Layout.fillWidth: true
                // the page takes the panel's width, however wide its rows would like to be
                Layout.preferredWidth: 0
                Layout.preferredHeight: item ? item.implicitHeight : 0
                source: inspector.page
                visible: inspector.page !== ""
                onLoaded: item.width = Qt.binding(function () { return loader.width })
            }
            ColumnLayout {
                visible: inspector.page === ""
                Layout.fillWidth: true
                spacing: 8
                Note {
                    text: "Select a component, beam or analysis in the model tree or the 3D view " +
                          "to see and edit its properties."
                }
                Note { text: app.targetSummary }
                Note { text: app.beamsSummary }
                Button { text: "Target"; onClicked: Julia.br_select("target") }
                Button { text: "Beams"; onClicked: Julia.br_select("beams") }
            }
        }
    }
}
