// The Help screen: the explanations of every page (common/help.js), a page list on the
// left and the page's topics on the right; a filter above the list shows the topics of
// every page that contain its words instead. `openAt(topic)` shows the topic's page with
// the topic marked, as an ⓘ does when clicked; Documentation opens the documentation's
// Desktop application page.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "common"
import "common/help.js" as Help

Dialog {
    id: help
    title: "Help"
    modal: false
    anchors.centerIn: parent
    width: Math.min(parent ? parent.width - 80 : 900, 900)
    height: Math.min(parent ? parent.height - 80 : 640, 640)
    standardButtons: Dialog.Close
    property int page: 0
    property string marked: ""
    property string query: ""
    // the topics of every page whose title or text holds each word of the filter, titled
    // with their page
    readonly property var found: {
        const words = query.toLowerCase().split(/\s+/).filter(w => w !== "")
        if (words.length === 0) return []
        const out = []
        for (const p of Help.pages)
            for (const it of p.items) {
                const hay = (it.title + " " + it.text).toLowerCase()
                if (words.every(w => hay.indexOf(w) >= 0))
                    out.push({ key: it.key, title: p.title + " ▸ " + it.title, text: it.text })
            }
        return out
    }
    readonly property bool filtering: query.trim() !== ""
    function openAt(topic) {
        query = ""
        filter.text = ""
        marked = topic
        page = Help.pageOf(topic)
        open()
        topics.positionAt(topic)
    }

    RowLayout {
        anchors.fill: parent
        spacing: 12
        ColumnLayout {
            Layout.preferredWidth: 200
            Layout.fillWidth: false
            Layout.fillHeight: true
            spacing: 6
            TextField {
                id: filter
                Layout.fillWidth: true
                placeholderText: "Filter topics"
                selectByMouse: true
                onTextChanged: help.query = text
                tip: "Show the topics of every page that contain these words"
            }
            ListView {
                id: pageList
                Layout.fillWidth: true
                Layout.fillHeight: true
                clip: true
                model: Help.pages
                currentIndex: help.page
                delegate: ItemDelegate {
                    width: pageList.width
                    text: modelData.title
                    highlighted: !help.filtering && index === help.page
                    onClicked: { filter.text = ""; help.marked = ""; help.page = index }
                }
            }
            RowLayout {
                Layout.fillWidth: true
                Button {
                    text: "Documentation"
                    Layout.fillWidth: true
                    onClicked: Qt.openUrlExternally(Julia.br_docs_url())
                    tip: "The Desktop application page of the Briareus documentation (a local build, or online)"
                }
                Info { topic: "start.docs" }
            }
        }
        ListView {
            id: topics
            Layout.fillWidth: true
            Layout.fillHeight: true
            clip: true
            spacing: 14
            model: help.filtering ? help.found : Help.pages[help.page].items
            ScrollBar.vertical: ScrollBar {}
            function positionAt(topic) {
                const items = Help.pages[help.page].items
                for (let i = 0; i < items.length; ++i)
                    if (items[i].key === topic) positionViewAtIndex(i, ListView.Beginning)
            }
            header: Label {
                text: help.filtering ? (help.found.length === 0 ? "No topic matches" :
                                        help.found.length + (help.found.length === 1 ? " topic" : " topics"))
                                     : Help.pages[help.page].title
                font.pixelSize: 18
                font.weight: Font.DemiBold
                bottomPadding: 10
            }
            delegate: ColumnLayout {
                width: topics.width - 16
                spacing: 4
                Label {
                    text: modelData.title
                    font.weight: Font.DemiBold
                    color: modelData.key === help.marked ? palette.link : palette.windowText
                }
                Label {
                    text: modelData.text
                    wrapMode: Text.WordWrap
                    Layout.fillWidth: true
                }
            }
        }
    }
}
