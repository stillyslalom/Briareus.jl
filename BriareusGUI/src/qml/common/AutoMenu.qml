// A menu as wide as its widest entry. Qt Quick's Menu has a fixed default width, which
// clips longer entries ("Apply VisRad SST point…").
import QtQuick
import QtQuick.Controls

Menu {
    width: {
        let w = 0
        for (let i = 0; i < count; ++i) {
            const it = itemAt(i)
            if (it) w = Math.max(w, it.implicitWidth)
        }
        return Math.max(160, w + leftPadding + rightPadding)
    }
}
