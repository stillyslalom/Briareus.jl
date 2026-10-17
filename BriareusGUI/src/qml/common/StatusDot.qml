// A small dot colored by an analysis state: ok, attention, busy, todo.
import QtQuick

Rectangle {
    property string status: ""
    width: 9
    height: 9
    radius: 4.5
    color: status === "ok" ? "#2e9e44" : status === "attention" ? "#d08a00"
         : status === "busy" ? "#2f7fd8" : "#a0a0a0"
}
