// The launcher's splash (BriareusLauncher.start_splash): shown by Qt's `qml` tool within
// about a second, while the launcher loads BriareusGUI (15-20 s). Arguments after `--`:
// the icon's URL and the status file's URL. The file holds the stage shown, "quit", or
// "error: <message>" (the splash then stays, with a Close link, until clicked). Empty or
// missing for a second, it quits; after 30 min without "quit" (a lost launcher), too.
import QtQuick
import QtQuick.Window

Window {
    id: win
    width: 460
    height: 170
    visible: true
    title: "Briareus (starting)"
    flags: Qt.SplashScreen | Qt.FramelessWindowHint
    color: "#fbfbfd"

    readonly property var args: Qt.application.arguments
    property string stage: "Loading Briareus"
    property bool failed: false
    property int quarters: 0          // quarter seconds since the splash opened
    property int empty: 0             // consecutive empty reads of the status file

    function poll() {
        var x = new XMLHttpRequest()
        x.onreadystatechange = function () {
            if (x.readyState !== XMLHttpRequest.DONE)
                return
            var t = x.responseText.trim()
            if (t === "quit" || (t === "" && ++empty >= 4)) {
                Qt.quit()
            } else if (t.startsWith("error:")) {
                failed = true
                stage = t.substring(6).trim()
                ticker.stop()
            } else if (t !== "") {
                empty = 0
                stage = t
            }
        }
        x.open("GET", args[args.length - 1])
        x.send()
    }

    Rectangle { anchors.fill: parent; color: "transparent"; border.color: "#d0d0d8" }
    Image {
        id: icon
        source: args[args.length - 2]
        sourceSize: Qt.size(96, 96)
        anchors { left: parent.left; leftMargin: 28; verticalCenter: parent.verticalCenter }
    }
    Column {
        anchors { left: icon.right; leftMargin: 24; right: parent.right; rightMargin: 24
                  verticalCenter: parent.verticalCenter }
        spacing: 8
        Text { text: "Briareus"; font.pixelSize: 26; font.bold: true; color: "#202028" }
        Text {
            width: parent.width
            wrapMode: Text.Wrap
            text: failed ? stage : stage + " … " + Math.floor(quarters / 4) + " s"
            font.pixelSize: 14
            color: failed ? "#b3261e" : "#404048"
        }
        Text {
            // a start after an update compiles first (checking for that costs 2 s)
            visible: !failed && quarters > 160
            width: parent.width
            wrapMode: Text.Wrap
            text: "The first start after an update compiles for a few minutes."
            font.pixelSize: 12
            color: "#606068"
        }
        Row {
            visible: !failed
            spacing: 6
            Repeater {
                model: 5
                Rectangle {
                    width: 8; height: 8; radius: 4; color: "#4063D8"
                    opacity: quarters % 5 === index ? 1.0 : 0.25
                }
            }
        }
        Text {
            visible: failed
            text: "<u>Close</u>"
            font.pixelSize: 14
            color: "#4063D8"
            MouseArea { anchors.fill: parent; onClicked: Qt.quit() }
        }
    }
    Timer {
        id: ticker
        interval: 250; repeat: true; running: true
        onTriggered: { quarters += 1; poll() }
    }
    Timer { interval: 1800000; running: !failed; onTriggered: Qt.quit() }
    Component.onCompleted: poll()
}
