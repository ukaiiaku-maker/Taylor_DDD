#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXADIS_ROOT="${EXADIS_ROOT:-$ROOT_DIR/core/exadis}"
PATCH_FILE="${PATCH_FILE:-$ROOT_DIR/exadis_native_patches/0001-native-event-audit.patch}"
REQUIRED_SHA="20ea2e82cdb919581c0611c338a6e46f6ad3f008"

if [[ ! -d "$EXADIS_ROOT/.git" && ! -f "$EXADIS_ROOT/.git" ]]; then
  echo "error: ExaDiS checkout not found at $EXADIS_ROOT" >&2
  exit 2
fi
if [[ ! -f "$PATCH_FILE" ]]; then
  echo "error: native audit patch not found at $PATCH_FILE" >&2
  exit 2
fi

actual_sha="$(git -C "$EXADIS_ROOT" rev-parse HEAD)"
if [[ "$actual_sha" != "$REQUIRED_SHA" && "${ALLOW_EXADIS_SHA_MISMATCH:-0}" != "1" ]]; then
  echo "error: patch is pinned to ExaDiS $REQUIRED_SHA; found $actual_sha" >&2
  echo "set ALLOW_EXADIS_SHA_MISMATCH=1 only after reviewing the patch against that source" >&2
  exit 2
fi

if git -C "$EXADIS_ROOT" apply --reverse --check "$PATCH_FILE" 2>/dev/null; then
  echo "native audit patch is already applied"
  exit 0
fi

git -C "$EXADIS_ROOT" apply --check "$PATCH_FILE"
git -C "$EXADIS_ROOT" apply "$PATCH_FILE"
echo "applied native ExaDiS event-audit patch to $EXADIS_ROOT"
