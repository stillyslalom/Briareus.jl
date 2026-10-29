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
    // the widest entry of a list model (an item model: the box's own width)
    readonly property real widestEntry: {
        if (!Array.isArray(model)) return 0
        let w = 0
        for (const m of model) {
            const t = typeof m === "object" && m !== null ? String(m[textRole] ?? "") : String(m)
            w = Math.max(w, metrics.advanceWidth(t))
        }
        return w
    }
    popup.width: Math.max(width, Math.ceil(widestEntry) + 48)
}
