// The Help window: a window of its own (it resizes, moves and stays open beside the
// design) with the pages of common/help.js on the left and the page's topics on the
// right; a filter above the list shows the topics of every page that contain its words.
// `openAt(topic)` shows the topic's page with the topic marked, as an ⓘ does when
// clicked.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import QtQuick.Window
import jlqml
import "common"
import "common/help.js" as Help

Window {
    id: help
    title: "Briareus Help"
    width: 860
    height: 620
    minimumWidth: 480
    minimumHeight: 320
    color: palette.window
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
                const hay = (it.title + " " + (it.tip || "") + " " + it.text).toLowerCase()
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
        show()
        raise()
        requestActivate()
        Qt.callLater(topics.positionAt, topic)
    }
    Shortcut { sequence: "Escape"; onActivated: help.close() }

    // the window's background as an item, so that an image grabbed of the content has it
    // too (the window's own color is not part of its content item)
    Rectangle { anchors.fill: parent; color: palette.window }

    SplitView {
        anchors.fill: parent
        anchors.margins: 10
        orientation: Qt.Horizontal
        ColumnLayout {
            SplitView.preferredWidth: 210
            SplitView.minimumWidth: 140
            spacing: 6
            TextField {
                id: filter
                Layout.fillWidth: true
                Layout.rightMargin: 8
                placeholderText: "Search Help"
                selectByMouse: true
                onTextChanged: help.query = text
            }
            ListView {
                id: pageList
                Layout.fillWidth: true
                Layout.fillHeight: true
                Layout.rightMargin: 8
                clip: true
                model: Help.pages
                currentIndex: help.page
                boundsBehavior: Flickable.StopAtBounds
                delegate: ItemDelegate {
                    width: pageList.width
                    text: modelData.title
                    topPadding: 5
                    bottomPadding: 5
                    highlighted: !help.filtering && index === help.page
                    onClicked: { filter.text = ""; help.marked = ""; help.page = index; topics.positionViewAtBeginning() }
                }
            }
            Button {
                text: "Documentation"
                Layout.fillWidth: true
                Layout.rightMargin: 8
                onClicked: Qt.openUrlExternally(Julia.br_docs_url())
                tip: "Open the documentation in the browser"
            }
        }
        ListView {
            id: topics
            SplitView.fillWidth: true
            SplitView.minimumWidth: 260
            leftMargin: 14
            rightMargin: 14
            clip: true
            spacing: 18
            boundsBehavior: Flickable.StopAtBounds
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
                font.pixelSize: 20
                font.weight: Font.DemiBold
                bottomPadding: 12
            }
            delegate: ColumnLayout {
                // a readable line length, however wide the window
                width: Math.min(topics.width - topics.leftMargin - topics.rightMargin - 12, 640)
                spacing: 6
                Label {
                    text: modelData.title
                    font.pixelSize: 15
                    font.weight: Font.DemiBold
                    color: modelData.key === help.marked ? palette.link : palette.windowText
                    wrapMode: Text.WordWrap
                    Layout.fillWidth: true
                }
                Label {
                    text: modelData.text
                    textFormat: Text.MarkdownText
                    wrapMode: Text.WordWrap
                    lineHeight: 1.15
                    Layout.fillWidth: true
                }
            }
            footer: Item { height: 16 }
        }
    }
}
