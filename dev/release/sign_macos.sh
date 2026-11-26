#!/usr/bin/env bash
# Signs a macOS build of `briareus` (build.jl) with a Developer ID certificate, submits it
# for notarization and wraps it in a disk image. Run by the Release workflow when its
# secrets are set; it needs, in the environment, MACOS_CERTIFICATE (a .p12, base64),
# MACOS_CERTIFICATE_PASSWORD, MACOS_IDENTITY ("Developer ID Application: …"), APPLE_ID,
# APPLE_TEAM_ID and APPLE_APP_PASSWORD (an app-specific password). Julia compiles code
# at run time, so the hardened runtime allows JIT and unsigned executable memory
# (entitlements.plist).
#
#     dev/release/sign_macos.sh build/Briareus
set -euo pipefail
app="$(cd "$1" && pwd)"
here="$(cd "$(dirname "$0")" && pwd)"
keychain="$RUNNER_TEMP/briareus.keychain-db"
echo "$MACOS_CERTIFICATE" | base64 --decode > "$RUNNER_TEMP/cert.p12"
security create-keychain -p briareus "$keychain"
security set-keychain-settings -lut 21600 "$keychain"
security unlock-keychain -p briareus "$keychain"
security import "$RUNNER_TEMP/cert.p12" -P "$MACOS_CERTIFICATE_PASSWORD" -A -t cert -f pkcs12 -k "$keychain"
security list-keychain -d user -s "$keychain"
rm "$RUNNER_TEMP/cert.p12"
sign() { codesign --force --timestamp --options runtime --entitlements "$here/entitlements.plist" \
                   --sign "$MACOS_IDENTITY" "$@"; }
# every library and executable, the innermost first
while IFS= read -r -d '' f; do sign "$f"; done < <(find "$app" -type f \( -name '*.dylib' -o -name '*.so' -o -perm -u+x \) -print0)
sign "$app/bin/briareus"
dmg="$(dirname "$app")/Briareus.dmg"
hdiutil create -volname Briareus -srcfolder "$app" -ov -format UDZO "$dmg"
sign "$dmg"
xcrun notarytool submit "$dmg" --apple-id "$APPLE_ID" --team-id "$APPLE_TEAM_ID" --password "$APPLE_APP_PASSWORD" --wait
xcrun stapler staple "$dmg"
