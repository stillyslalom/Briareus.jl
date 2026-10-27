// The drive analysis: compute the view factors, the beam model of the laser deposition,
// what the target is colored by, the capsule or foil map, and the foil spot radius.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    Note {
        text: "View factors, the laser deposition operator and the capsule flux response, " +
              "once per target and beam ports (seconds on the standard mesh). Pointing, power " +
              "and on/off edits then re-evaluate the drive at once. Coloring the target by a " +
              "drive field computes it when needed."
    }
    RowLayout {
        spacing: 8
        Button {
            text: app.driveCurrent ? "Compute again" : "Compute drive"
            highlighted: !app.driveCurrent
            enabled: !app.driveBuilding
            onClicked: Julia.br_build()
        }
        BusyIndicator {
            running: app.driveBuilding || app.driveEvaluating
            visible: running
            Layout.preferredHeight: 28
            Layout.preferredWidth: 28
        }
    }
    Note {
        text: app.driveBuilding ? app.driveStatus + " " + app.driveElapsed :
              app.driveStale ? "The target, the beam ports, the cavity or the beam model changed since the drive was computed." :
              app.driveStatus
        color: app.driveStale ? "#d08a00" : palette.windowText
    }

    SectionLabel { text: "Beam model" }
    ComboBox {
        Layout.fillWidth: true
        textRole: "text"; valueRole: "value"
        model: [{ text: "Briareus (f-cone)", value: "fcone" },
                { text: "VisRad profile", value: "visrad" }]
        currentIndex: Math.max(0, ["fcone", "visrad"].indexOf(app.beamModel))
        onActivated: Julia.br_beam_model(currentValue)
    }
    Note {
        text: app.beamModel === "visrad" ?
              "VisRad's beam model (a super-Gaussian cut at each beam's envelope), so the design " +
              "deposits in VisRad as here; an imported workspace sets its envelopes and focal length." :
              "The focal spot summed over the f-cone's directions, each with its own incidence."
    }

    SectionLabel { text: "Cavity" }
    CheckBox {
        text: "Compute the drive on the cavity around a point"
        checked: app.driveCavity
        onToggled: Julia.br_cavity(checked)
    }
    PropRow {
        visible: app.driveCavity
        label: "Point"
        value: app.cavityPoint
        unit: "x, y, z mm"
        fieldWidth: 140
        tip: "A point in the free space of the cavity (inside the capsule is fine)"
        onEdited: (text) => Julia.br_cavity_point(text)
    }
    Note {
        visible: app.driveCavity
        text: "Only the elements that bound the free space around the point take part in the " +
              "drive: a CAD assembly's outer faces and mounting hardware drop out of the view " +
              "factors. Mark thin films modeled as solids (a capsule tent, an entrance-hole " +
              "window) in their Properties, or they seal the point into a pocket. The clearance " +
              "and TaLIS checks use the whole target."
    }
    Note { visible: app.driveCavity && app.cavitySummary !== ""; text: app.cavitySummary }

    SectionLabel { text: "Show" }
    Label { text: "Color the target by" }
    RowLayout {
        spacing: 8
        ComboBox {
            Layout.fillWidth: true
            textRole: "text"; valueRole: "value"
            model: [{ text: "Components", value: "component" },
                    { text: "Absorbed laser intensity", value: "laser" },
                    { text: "X-ray radiosity", value: "radiosity" },
                    { text: "Incident x-ray flux", value: "incident" }]
            currentIndex: Math.max(0, ["component", "laser", "radiosity", "incident"].indexOf(app.colorBy))
            onActivated: Julia.br_set_color(currentValue)
        }
    }
    Flow {
        spacing: 8
        Layout.fillWidth: true
        Button {
            text: "3D target"; checkable: true; checked: app.driveView === "scene"
            onClicked: Julia.br_set_view("scene")
        }
        Button {
            text: app.hasFoil ? "Foil face" : "Capsule θ–φ map"; checkable: true
            checked: app.driveView === "map"
            enabled: app.driveHasResult && (app.hasCapsule || app.hasFoil)
            onClicked: Julia.br_set_view("map")
        }
    }
    PropRow {
        visible: app.hasFoil
        label: "Spot radius"
        value: app.spotRadius
        unit: "mm about the foil center"
        tip: "The spot over which the irradiance statistics and the foil objective are taken"
        onEdited: (text) => Julia.br_set_spot_radius(text)
    }
    Label { text: app.driveSummary; visible: app.driveHasResult; wrapMode: Text.WordWrap; Layout.fillWidth: true; font.weight: Font.DemiBold }
    Note { text: "The modes and the power accounting are in the Drive tab below." }
}
