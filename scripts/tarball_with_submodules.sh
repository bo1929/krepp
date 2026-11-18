#!/usr/bin/env bash

# Package a tagged source tree and its pinned submodules.
set -euo pipefail

fail() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

[[ $# -le 1 ]] || fail "usage: $0 [destination]"

case "$(uname -s)" in
  Darwin) tar_command=gtar ;;
  Linux)  tar_command=tar ;;
  *)      fail "unsupported platform" ;;
esac
for command_name in git gzip "$tar_command"; do
  command -v "$command_name" >/dev/null 2>&1 || fail "$command_name is required"
done

repo_root=$(git rev-parse --show-toplevel)
tag=$(git -C "$repo_root" describe --tags --exact-match --match 'v[0-9]*' HEAD 2>/dev/null) ||
  fail "check out a release tag before creating an archive"
[[ "$tag" != */* ]] || fail "release tag must not contain a slash: $tag"
git -C "$repo_root" diff --cached --quiet --ignore-submodules=none HEAD -- &&
  git -C "$repo_root" diff --quiet --ignore-submodules=all -- ||
  fail "tracked files must match the release tag"

archive_prefix="krepp-${tag#v}"
destination=${1:-"$archive_prefix-including-submodules.tar.gz"}
[[ ! -d "$destination" ]] || fail "destination must be a file path: $destination"

scratch=$(mktemp -d "${TMPDIR:-/tmp}/krepp-archive.XXXXXX")
trap 'rm -rf -- "$scratch"' EXIT
archive_file="$scratch/source.tar"
submodule_archive="$scratch/submodule.tar"
export archive_prefix archive_file submodule_archive tar_command

git -C "$repo_root" archive --prefix="$archive_prefix/" --output="$archive_file" "$tag"
git -C "$repo_root" submodule update --init --recursive --checkout
git -C "$repo_root" submodule foreach --recursive '
  git diff --cached --quiet --ignore-submodules=none "$sha1" -- &&
    git archive --prefix="$archive_prefix/$displaypath/" --output="$submodule_archive" "$sha1" &&
    "$tar_command" --concatenate --file="$archive_file" "$submodule_archive"
'

gzip -n -9 "$archive_file"
cp -- "$archive_file.gz" "$destination"
printf 'created %s\n' "$destination"
