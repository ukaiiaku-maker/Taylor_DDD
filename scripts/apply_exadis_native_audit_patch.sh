#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXADIS_ROOT="${EXADIS_ROOT:-$ROOT_DIR/core/exadis}"
PATCH_FILE="${PATCH_FILE:-}"
REQUIRED_SHA="20ea2e82cdb919581c0611c338a6e46f6ad3f008"

if [[ ! -d "$EXADIS_ROOT/.git" && ! -f "$EXADIS_ROOT/.git" ]]; then
  echo "error: ExaDiS checkout not found at $EXADIS_ROOT" >&2
  exit 2
fi
if [[ -n "$PATCH_FILE" ]]; then
  patches=("$PATCH_FILE")
else
  patches=(
    "$ROOT_DIR/exadis_native_patches/0001-native-event-audit.patch"
    "$ROOT_DIR/exadis_native_patches/0002-native-arrhenius-exp-floor.patch"
  )
fi
for patch in "${patches[@]}"; do
  if [[ ! -f "$patch" ]]; then
    echo "error: native ExaDiS patch not found at $patch" >&2
    exit 2
  fi
done

actual_sha="$(git -C "$EXADIS_ROOT" rev-parse HEAD)"
if [[ "$actual_sha" != "$REQUIRED_SHA" && "${ALLOW_EXADIS_SHA_MISMATCH:-0}" != "1" ]]; then
  echo "error: patch is pinned to ExaDiS $REQUIRED_SHA; found $actual_sha" >&2
  echo "set ALLOW_EXADIS_SHA_MISMATCH=1 only after reviewing the patch against that source" >&2
  exit 2
fi

for patch in "${patches[@]}"; do
  patch_name="$(basename "$patch")"
  patch_present=0
  case "$patch_name" in
    0001-native-event-audit.patch)
      if [[ -f "$EXADIS_ROOT/src/audit/event_audit.h" ]] &&
         grep -q "NATIVE_EVENT_AUDIT_COMPILED" "$EXADIS_ROOT/python/exadis_pybind.cpp"; then
        patch_present=1
      fi
      ;;
    0002-native-arrhenius-exp-floor.patch)
      if [[ -f "$EXADIS_ROOT/src/arrhenius/arrhenius_exp_floor.h" ]] &&
         [[ -f "$EXADIS_ROOT/src/mobility_types/mobility_fcc0_arrhenius.h" ]] &&
         grep -q "arrhenius/arrhenius_exp_floor.cpp" "$EXADIS_ROOT/src/CMakeLists.txt"; then
        patch_present=1
      fi
      ;;
  esac
  if [[ $patch_present -eq 1 ]]; then
    echo "native patch is already applied: $patch_name"
    continue
  fi
  if git -C "$EXADIS_ROOT" apply --reverse --check "$patch" 2>/dev/null; then
    echo "native patch is already applied: $patch_name"
    continue
  fi
  git -C "$EXADIS_ROOT" apply --check "$patch"
  git -C "$EXADIS_ROOT" apply "$patch"
  echo "applied native ExaDiS patch: $patch_name"
done
