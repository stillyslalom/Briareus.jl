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
            tip: "Build the view factors and compute the drive (Ctrl+R)"
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
        tip: app.driveStale && !app.driveBuilding ? "The target, the ports or the settings changed" : ""
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
        tip: "How each beam deposits its light"
    }

    SectionLabel { text: "Cavity"; info: "drive.cavity" }
    CheckLine {
        text: "Compute the drive on the cavity around a point"
        checked: app.driveCavity
        onToggled: Julia.br_cavity(checked)
        tip: "Use only the surfaces around the point for the x-ray drive"
    }
    PropRow {
        visible: app.driveCavity
        label: "Point"
        value: app.cavityPoint
        unit: "x, y, z mm"
        fieldWidth: 140
        tip: "A point in the cavity's free space"
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
            tip: "Color by component, or by a drive field (W/cm²)"
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
            tip: "The x-ray flux map, or the foil's laser irradiance"
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
            Tip {
                visible: parent.hovered
                text: "The foil the drive and the optimization report on"
            }
        }
    }
    PropRow {
        visible: app.hasFoil
        label: "Spot radius"
        value: app.spotRadius
        unit: "mm about the foil center"
        tip: "Radius of the spot the foil statistics cover"
        onEdited: (text) => Julia.br_set_spot_radius(text)
    }
    SectionLabel { text: "Unconverted light"; info: "drive.unconverted" }
    PropRow {
        label: "1ω power ratio"; value: app.ucRatio1; unit: "of each beam's 3ω power (empty: off)"
        tip: "1ω power as a fraction of the 3ω power; no default"
        onEdited: (text) => Julia.br_unconverted_ratio(1, text)
    }
    PropRow {
        label: "2ω power ratio"; value: app.ucRatio2; unit: "of each beam's 3ω power (empty: off)"
        tip: "2ω power as a fraction of the 3ω power; no default"
        onEdited: (text) => Julia.br_unconverted_ratio(2, text)
    }
    PropRow {
        label: "1ω prepulse"; value: app.ppOmega1; unit: "TW per beam (empty: off)"
        tip: "Flat 1ω power per beam over the prepulse window; no default"
        onEdited: (text) => Julia.br_prepulse(1, text)
    }
    PropRow {
        label: "2ω prepulse"; value: app.ppOmega2; unit: "TW per beam (empty: off)"
        tip: "As the 1ω prepulse, on the beams whose 2ω light takes part"
        onEdited: (text) => Julia.br_prepulse(2, text)
    }
    PropRow {
        label: "Prepulse from"; value: app.ppStart; unit: "ns"
        tip: "Start of the prepulse window"
        onEdited: (text) => Julia.br_prepulse(3, text)
    }
    PropRow {
        label: "Prepulse to"; value: app.ppStop; unit: "ns"
        tip: "The end of the prepulse window"
        onEdited: (text) => Julia.br_prepulse(4, text)
    }
    PropRow {
        label: "Plasma threshold"; value: app.irrThreshold; unit: "W/cm² (irradiance checks)"
        tip: "The irradiance the checks compare with (default 10⁸ W/cm²)"
        onEdited: (text) => Julia.br_irr_threshold(text)
    }
    Label { text: app.unconvertedText; visible: app.unconvertedText !== ""; wrapMode: Text.WordWrap; Layout.fillWidth: true }
    Label { text: app.driveSummary; visible: app.driveHasResult; wrapMode: Text.WordWrap; Layout.fillWidth: true; font.weight: Font.DemiBold }
}
