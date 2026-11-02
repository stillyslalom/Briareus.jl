// The style's ComboBox at the window's density (see TextField.qml), its list as wide as
// its widest entry so that no choice is clipped there (the box itself may show a short
// `displayText`).
import QtQuick
import QtQuick.Controls as T

T.ComboBox {
    id: box
    topPadding: 2
    bottomPadding: 2
    implicitHeight: Math.max(26, implicitContentHeight + topPadding + bottomPadding)
    FontMetrics { id: metrics; font: box.font }
    // the widest entry of a list model, a JavaScript array or a list from Julia (an item
    // model: the box's own width)
    readonly property real widestEntry: {
        if (model === null || model === undefined || typeof model.length !== "number") return 0
        let w = 0
        for (let i = 0; i < model.length; i++) {
            const m = model[i]
            const t = typeof m === "object" && m !== null ? String(m[textRole] ?? "") : String(m)
            w = Math.max(w, metrics.advanceWidth(t))
        }
        return w
    }
    popup.width: Math.max(width, Math.ceil(widestEntry) + 48)
    // The FluentWinUI3 background is an ItemGroup around a nine-patch image. Placed by a
    // positioner (the toolbar's Flow) rather than a layout, the box keeps the image at
    // height 0 (measured: group 26 px, image 0 px), drawn as its top edge alone, so the
    // background's children follow its size here.
    function fitBackground() {
        if (!background) return
        for (const c of background.children) {
            c.width = background.width
            c.height = background.height
        }
    }
    Connections {
        target: box.background
        function onWidthChanged() { box.fitBackground() }
        function onHeightChanged() { box.fitBackground() }
    }
    Component.onCompleted: fitBackground()
}
