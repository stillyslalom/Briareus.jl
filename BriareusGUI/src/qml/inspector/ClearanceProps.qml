// The clearance check page sets component roles, the envelope and error set, and the
// sampling, and runs the check.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    SectionLabel { text: "Components" }
    Repeater {
        model: compModel
        delegate: RoleRow {
            enabled: !app.rvRunning
            label: model.label
            role: model.role
            readonly property string compName: model.name
            roles: [{ text: "beams land on it", value: "target",
                      tip: "The beams land on it" },
                    { text: "check clearance", value: "hardware",
                      tip: "The beams must stay clear of it" },
                    { text: "ignore", value: "none", tip: "Left out of the check" }]
            onChosen: (value) => Julia.br_rv_role(compName, value)
        }
    }
    SectionLabel { text: "Envelope and errors"; info: "clearance.errors" }
    ColumnLayout {
        spacing: 6
        enabled: !app.rvRunning && !app.tlRunning
        Layout.fillWidth: true
        PropRow { label: "Contour"; value: app.rvPeak; unit: "of each spot's peak"
                  tip: "The spot contour that bounds the beam, as a fraction of its peak"
                  onEdited: (text) => Julia.br_rv_edit("peak_fraction", text) }
        PropRow { label: "Pointing error"; value: app.rv_pointing; unit: "mm"
                  tip: "Displacement of each aim point in any direction across the beam"
                  onEdited: (text) => Julia.br_rv_edit("pointing", text) }
        PropRow { label: "Target positioning"; value: app.rv_target_offset; unit: "mm"
                  tip: "Displacement of the whole target from its nominal position"
                  onEdited: (text) => Julia.br_rv_edit("target_offset", text) }
        PropRow { label: "Hardware tolerance"; value: app.rv_hardware_tolerance; unit: "mm"
                  tip: "Fabrication and assembly tolerance of the hardware"
                  onEdited: (text) => Julia.br_rv_edit("hardware_tolerance", text) }
        PropRow { label: "Focus error"; value: app.rv_focus; unit: "± mm"
                  tip: "Shift of best focus along each beam, both ways"
                  onEdited: (text) => Julia.br_rv_edit("focus", text) }
        CheckLine { text: "Far-field wings (1% halo, 10⁻³ scatter)"; checked: app.rvWings; onToggled: Julia.br_rv_wings(checked)
                    tip: "Count the NIF far-field wings in the envelope" }
        RowLayout {
            spacing: 8
            Cell { text: "Sampling"; Layout.preferredWidth: 150 }
            ComboBox {
                Layout.fillWidth: true
                textRole: "text"; valueRole: "value"
                model: [{ text: "Quick", value: "quick" }, { text: "Standard (4× the points)", value: "standard" }]
                currentIndex: app.rvQuality === "standard" ? 1 : 0
                onActivated: Julia.br_rv_quality(currentValue)
                tip: "Quick: coarser sampling, faster. Standard: finer"
            }
        }
    }
    ErrorText { text: app.rvError }
    RowLayout {
        spacing: 8
        Button { text: "Check clearance"; highlighted: true; enabled: !app.rvRunning; onClicked: Julia.br_rv_start()
                 tip: "Runs in the background; results in the Clearance tab" }
        Button { text: "Cancel"; enabled: app.rvRunning; onClicked: Julia.br_rv_cancel()
                 tip: "Stop the check and keep the beams checked so far" }
    }
    ProgressBar {
        Layout.fillWidth: true
        from: 0
        to: Math.max(1, app.rvTotal)
        value: app.rvDone
        visible: app.rvRunning
    }
    Note {
        text: app.rvRunning ? "checked " + app.rvDone + " of " + app.rvTotal + " beams" :
              app.rvStale ? "Out of date. " + app.rvBrief : app.rvStatus
        tip: app.rvStale && !app.rvRunning ? "The beams, the target or the settings changed since this check" : ""
        color: app.rvStale && !app.rvRunning ? "#d08a00" : palette.windowText
        opacity: 1
    }
    // the keep-out map: where new hardware fits between the beams and the lines of sight
    SectionLabel { text: "Keep-out map" }
    CheckLine { text: "Show the keep-out map"; checked: app.cmShown; onToggled: Julia.br_cm_show(checked)
                tip: "Colors a sphere about TCC or a plane by what crosses it: 3ω light (blue), 1ω light (vermilion), the lines of sight of the ticked diagnostics (purple). The open region is green, and yellow where the margin is below the tight threshold. Uses the errors above and the TaLIS rule tolerance. The Clearance tab draws it flat" }
    RowLayout {
        spacing: 8
        Cell { text: "Surface"; Layout.preferredWidth: 150 }
        ComboBox {
            Layout.fillWidth: true
            textRole: "text"; valueRole: "value"
            model: [{ text: "Sphere about TCC", value: "sphere" }, { text: "Plane", value: "plane" }]
            currentIndex: app.cmSurface === "plane" ? 1 : 0
            onActivated: Julia.br_cm_surface(currentValue)
            tip: "Sphere: the directions from TCC at one distance. Plane: a slice through the mounting zone"
        }
    }
    PropRow { visible: app.cmSurface !== "plane"; label: "Radius"; value: app.cmRadius; unit: "mm from TCC"
              tip: "The distance from TCC at which the hardware is to be placed"
              onEdited: (text) => Julia.br_cm_radius(text) }
    PropRow { visible: app.cmSurface === "plane"; label: "Point"; value: app.cmPoint; unit: "x, y, z mm"
              tip: "A point of the plane, at the center of the map"
              onEdited: (text) => Julia.br_cm_edit("point", text) }
    PropRow { visible: app.cmSurface === "plane"; label: "Normal"; value: app.cmNormal; unit: "x, y, z"
              tip: "The direction across the plane"
              onEdited: (text) => Julia.br_cm_edit("normal", text) }
    PropRow { visible: app.cmSurface === "plane"; label: "Half-width"; value: app.cmExtent; unit: "mm"
              tip: "The map spans this distance on each side of the point"
              onEdited: (text) => Julia.br_cm_edit("extent", text) }
    ErrorText { text: app.cmError }
    Note {
        visible: app.cmShown
        text: app.cmStale ? "Out of date. " + app.cmText : app.cmText
        tip: app.cmStale ? "The beams, the diagnostics or the settings changed since this map: switch it off and on to update" : ""
        color: app.cmStale ? "#d08a00" : palette.windowText
        opacity: 1
    }
}
