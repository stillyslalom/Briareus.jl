// The drive analysis: compute the view factors, the beam model of the laser deposition,
// what the target is colored by, the capsule or foil map, and the foil spot radius.
import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
import jlqml
import "../common"

ColumnLayout {
    spacing: 8
    RowLayout {
        spacing: 8
        Button {
            text: app.driveCurrent ? "Compute again" : "Compute drive"
            highlighted: !app.driveCurrent
            enabled: !app.driveBuilding
            onClicked: Julia.br_build()
            tip: "View factors, the laser deposition operator and the flux response, then the drive (Ctrl+R)"
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
              app.driveStale ? "Out of date: compute again" :
              app.driveStatus
        tip: app.driveStale && !app.driveBuilding ? "The target, the beam ports, the cavity or the beam model changed since the drive was computed" : ""
        color: app.driveStale ? "#d08a00" : palette.windowText
    }

    SectionLabel { text: "Beam model"; info: "drive.model" }
    ComboBox {
        Layout.fillWidth: true
        textRole: "text"; valueRole: "value"
        model: [{ text: "Briareus (f-cone)", value: "fcone" },
                { text: "VisRad profile", value: "visrad" }]
        currentIndex: Math.max(0, ["fcone", "visrad"].indexOf(app.beamModel))
        onActivated: Julia.br_beam_model(currentValue)
        tip: "How each beam deposits its light: the f-cone sum, or VisRad's profile cut at each beam's envelope"
    }

    SectionLabel { text: "Cavity"; info: "drive.cavity" }
    CheckLine {
        text: "Compute the drive on the cavity around a point"
        checked: app.driveCavity
        onToggled: Julia.br_cavity(checked)
        tip: "Only the elements that bound the free space around the point take part in the view factors"
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
            tip: "Components: each its own color; or a drive field per element (W/cm²), computing the drive when needed"
        }
    }
    Flow {
        spacing: 8
        Layout.fillWidth: true
        Button {
            text: "3D target"; checkable: true; checked: app.driveView === "scene"
            onClicked: Julia.br_set_view("scene")
            tip: "The 3D view of the target"
        }
        Button {
            text: app.observedKind === "face" ? "Observed face" :
                  app.observedKind === "capsule" ? "Capsule θ–φ map" : "Foil face"
            checkable: true
            checked: app.driveView === "map"
            enabled: app.driveHasResult && (app.hasCapsule || app.hasFoil)
            onClicked: Julia.br_set_view("map")
            tip: "The observed surface's x-ray flux map (the capsule's θ–φ map or the flat face's), or the foil's laser irradiance"
        }
    }
    RowLayout {
        visible: app.foilNames.length > 1
        spacing: 8
        Layout.fillWidth: true
        Label { text: "Active foil"; Layout.preferredWidth: 140 }
        ComboBox {
            Layout.fillWidth: true
            model: app.foilNames
            currentIndex: app.activeFoil
            onActivated: (index) => Julia.br_active_foil(index)
            ToolTip.visible: hovered
            ToolTip.delay: 400
            ToolTip.text: "The foil whose spot statistics, map and objective the drive and the optimization take, and the sweep's line of sight and captured power"
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
    SectionLabel { text: "Unconverted light"; info: "drive.unconverted" }
    PropRow {
        label: "1ω power ratio"; value: app.ucRatio1; unit: "of each beam's 3ω power (empty: off)"
        tip: "No default: the NIF User Guide gives no fraction (most unconverted light is 1ω, more for shaped pulses)"
        onEdited: (text) => Julia.br_unconverted_ratio(1, text)
    }
    PropRow {
        label: "2ω power ratio"; value: app.ucRatio2; unit: "of each beam's 3ω power (empty: off)"
        tip: "The 2ω power as a fraction of the 3ω power, on the beams whose 2ω light takes part; no default"
        onEdited: (text) => Julia.br_unconverted_ratio(2, text)
    }
    PropRow {
        label: "1ω prepulse"; value: app.ppOmega1; unit: "TW per beam (empty: off)"
        tip: "Light not proportional to the 3ω pulse, flat over the prepulse window, on the beams whose 1ω light takes part; the drive history, the retune and the irradiance checks add it (the drive at the table's powers does not). No default level (derived scale: 4.8e-4 TW of 1ω on one beamlet puts its 22 mm footprint at 10⁸ W/cm²)"
        onEdited: (text) => Julia.br_prepulse(1, text)
    }
    PropRow {
        label: "2ω prepulse"; value: app.ppOmega2; unit: "TW per beam (empty: off)"
        tip: "As the 1ω prepulse, on the beams whose 2ω light takes part"
        onEdited: (text) => Julia.br_prepulse(2, text)
    }
    PropRow {
        label: "Prepulse from"; value: app.ppStart; unit: "ns"
        tip: "The NIF User Guide (§4.4) puts the 1ω prepulse above 10⁸ W/cm² 3–5 ns before the nominal start of the pulse; -5 to 0 ns is an engineering default"
        onEdited: (text) => Julia.br_prepulse(3, text)
    }
    PropRow {
        label: "Prepulse to"; value: app.ppStop; unit: "ns"
        tip: "The end of the prepulse window"
        onEdited: (text) => Julia.br_prepulse(4, text)
    }
    PropRow {
        label: "Plasma threshold"; value: app.irrThreshold; unit: "W/cm² (irradiance checks)"
        tip: "The irradiance the Time tab's checks compare with: the NIF User Guide's guideline for the onset of weakly ionized plasma on a metal (10⁸ W/cm², §4.4), taken as the end of specular reflection; an engineering default for metals only"
        onEdited: (text) => Julia.br_irr_threshold(text)
    }
    Label { text: app.unconvertedText; visible: app.unconvertedText !== ""; wrapMode: Text.WordWrap; Layout.fillWidth: true }
    Label { text: app.driveSummary; visible: app.driveHasResult; wrapMode: Text.WordWrap; Layout.fillWidth: true; font.weight: Font.DemiBold }
}
