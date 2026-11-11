#!/usr/bin/env bash
# Cut a spearfish release: bump the workspace version, date its CHANGELOG section, commit,
# publish the crates, tag `vX.Y.Z`, and push. Run from a clean, up-to-date `main`.
set -euo pipefail

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./bump_and_publish.sh <version> [--publish] [--dry-run] [--allow-same-version]
                                  [--github-release] [--python]

Bumps the workspace version (the spearfish crates and the Python package share it) and the
workspace's own dependency versions, dates the CHANGELOG section `## <version> (unreleased)`,
commits, checks that every crate builds from its own package (`cargo publish --workspace
--dry-run`), tags `v<version>`, and pushes the commit and the tag to origin.

A crate that sets a version of its own (crates/nudgy) keeps it: it is published with the
release when that version is not on crates.io yet, and the release stops if the crate has
changed since the version it names was published.

Options:
  --publish             Publish the crates to crates.io (`cargo publish --workspace`, in
                        dependency order) before tagging. Crates already published at
                        their version are skipped, so an interrupted release can be re-run.
  --allow-same-version  Proceed when the workspace is already at <version> (the version edits
                        are skipped; the CHANGELOG is still dated, and check/publish/tag run).
  --github-release      Also make a GitHub release from the tag, with the CHANGELOG section
                        as its notes (needs `gh`).
  --python              Also tag and push `py-spearfish-v<version>`, which publishes the
                        Python package to PyPI (.github/workflows/publish-python.yml).
  --dry-run             Show what would be done without changing files, committing,
                        publishing, tagging or pushing.
  -h, --help            Show this help message.

Examples:
  ./bump_and_publish.sh 0.1.0 --allow-same-version --dry-run
  ./bump_and_publish.sh 0.1.0 --allow-same-version --publish --github-release
  ./bump_and_publish.sh 0.2.0 --publish --github-release --python
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

VERSION=""
PUBLISH=false
DRY_RUN=false
ALLOW_SAME_VERSION=false
GITHUB_RELEASE=false
PYTHON=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --publish) PUBLISH=true ;;
        --dry-run) DRY_RUN=true ;;
        --allow-same-version) ALLOW_SAME_VERSION=true ;;
        --github-release) GITHUB_RELEASE=true ;;
        --python) PYTHON=true ;;
        -h | --help)
            usage
            exit 0
            ;;
        -*) die "unknown option: $1" ;;
        *)
            [[ -z "$VERSION" ]] || die "version specified more than once"
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || {
    usage
    exit 1
}

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ROOT_CARGO="Cargo.toml"
LOCKFILE="Cargo.lock"
CHANGELOG="CHANGELOG.md"
TAG="v${VERSION}"
PY_TAG="py-spearfish-v${VERSION}"
REPO_URL="https://github.com/COMBINE-lab/spearfish"
TODAY="$(date +%Y-%m-%d)"

for f in "$ROOT_CARGO" "$LOCKFILE" "$CHANGELOG"; do
    [[ -f "$f" ]] || die "not found: $f"
done

# ── Pre-flight: a clean main, in sync with origin, without the release's tags ─────────────────

[[ -z "$(git status --porcelain)" ]] || die "working tree is not clean; commit or stash existing changes first"

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
[[ "$BRANCH" == "main" ]] || die "releases are cut from main (on: $BRANCH)"
git fetch -q origin main
RESUMING=false
if [[ "$(git rev-parse HEAD)" != "$(git rev-parse origin/main)" ]]; then
    # An interrupted release leaves its release commit unpushed: main is then exactly that
    # commit ahead of origin/main, and the release resumes from it.
    if git merge-base --is-ancestor origin/main HEAD &&
        [[ "$(git rev-list --count origin/main..HEAD)" == 1 ]] &&
        [[ "$(git log -1 --format=%s)" == "chore(release): spearfish v${VERSION}" ]]; then
        RESUMING=true
        echo "Resuming the release of v${VERSION} from its unpushed release commit"
    else
        die "main is not at origin/main; pull (or push) first"
    fi
fi

for t in "$TAG" $([[ "$PYTHON" == true ]] && echo "$PY_TAG"); do
    if git rev-parse -q --verify "refs/tags/$t" >/dev/null ||
        [[ -n "$(git ls-remote --tags origin "refs/tags/$t")" ]]; then
        die "tag $t already exists"
    fi
done

# ── Versions ────────────────────────────────────────────────────────────────────────────────

# `version` in [workspace.package]: every crate inherits it.
workspace_version() {
    awk -F'"' '/^\[workspace\.package\]/{p=1; next} /^\[/{p=0} p && /^version[[:space:]]*=/{print $2; exit}' "$ROOT_CARGO"
}

# The versions of the workspace's own crates in [workspace.dependencies] (entries with both a
# `path = "crates/..."` and a `version`), one per line as `name version`.
internal_dep_versions() {
    awk '/^\[workspace\.dependencies\]/{p=1; next} /^\[/{p=0}
         p && /path = "crates\// && /version = "/ {
             v = $0; sub(/.*version = "/, "", v); sub(/".*/, "", v); print $1, v
         }' "$ROOT_CARGO"
}

CURRENT_VERSION="$(workspace_version)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine the [workspace.package] version from $ROOT_CARGO"
[[ -n "$(internal_dep_versions)" ]] || die "found no versioned workspace crates in [workspace.dependencies]"

SAME_VERSION=false
if [[ "$CURRENT_VERSION" == "$VERSION" ]]; then
    if [[ "$RESUMING" == true ]]; then
        SAME_VERSION=true
    elif [[ "$ALLOW_SAME_VERSION" == true ]]; then
        SAME_VERSION=true
        echo "Workspace already at $VERSION; skipping the version edits (--allow-same-version)"
    else
        die "workspace version is already $VERSION (pass --allow-same-version to release it as-is)"
    fi
fi

# ── CHANGELOG: the release needs its own section ──────────────────────────────────────────────

VERSION_RE="${VERSION//./\\.}"
UNRELEASED_HEADING="$(grep -m1 -E "^## ${VERSION_RE} \(unreleased\)[[:space:]]*$" "$CHANGELOG" || true)"
DATED_HEADING="$(grep -m1 -E "^## ${VERSION_RE} \([0-9]{4}-[0-9]{2}-[0-9]{2}\)[[:space:]]*$" "$CHANGELOG" || true)"
if [[ -n "$UNRELEASED_HEADING" ]]; then
    NEW_HEADING="## ${VERSION} (${TODAY})"
elif [[ -n "$DATED_HEADING" ]]; then
    NEW_HEADING="$DATED_HEADING"
else
    die "$CHANGELOG has no section for $VERSION: add '## ${VERSION} (unreleased)' with the release notes"
fi

# ── Tools and crates ────────────────────────────────────────────────────────────────────────

cargo publish --help | grep -q -- '--workspace' ||
    die "this cargo cannot publish a workspace in one command (needs cargo 1.90 or newer)"
if [[ "$GITHUB_RELEASE" == true ]]; then
    command -v gh >/dev/null || die "--github-release needs the GitHub CLI (gh)"
fi

# The crates to publish: every workspace member that does not say `publish = false`, one per
# line as `name version directory`. The version is `-` for a crate at the workspace version; a
# crate that sets a version of its own is released at that one.
PUBLISHABLE="$(cargo metadata --no-deps --format-version 1 --locked | python3 -c '
import json, os, re, sys
shared = re.compile(r"^version(\.workspace\s*=\s*true|\s*=\s*\{\s*workspace\s*=\s*true\s*\})", re.M)
for p in json.load(sys.stdin)["packages"]:
    if p.get("publish") == []:
        continue
    with open(p["manifest_path"]) as f:
        own = not shared.search(f.read())
    print(p["name"], p["version"] if own else "-", os.path.relpath(os.path.dirname(p["manifest_path"])))
')"
[[ -n "$PUBLISHABLE" ]] || die "found no publishable crates"

# The version a crate is released at, from the second field of its PUBLISHABLE line.
crate_version() {
    if [[ "$1" == "-" ]]; then echo "$VERSION"; else echo "$1"; fi
}

# Whether crate $1 is on crates.io at version $2.
on_crates_io() {
    local code
    code="$(curl -s -o /dev/null -w '%{http_code}' -A "spearfish-release (${REPO_URL})" \
        "https://crates.io/api/v1/crates/$1/$2")"
    [[ "$code" == "200" ]]
}

# Whether directory $3 holds what crate $1 was published with at version $2: the package's
# files are the directory's committed ones. (The package keeps the Cargo.toml as it was
# written as Cargo.toml.orig, and adds a lockfile and a note of the commit it was made from.)
unchanged_since_published() {
    local tmp same=true
    tmp="$(mktemp -d)"
    mkdir "$tmp/published" "$tmp/here"
    curl -sfL -A "spearfish-release (${REPO_URL})" \
        "https://crates.io/api/v1/crates/$1/$2/download" |
        tar -xzf - -C "$tmp/published" --strip-components 1 ||
        die "could not download $1 v$2 from crates.io"
    mv "$tmp/published/Cargo.toml.orig" "$tmp/published/Cargo.toml" ||
        die "the package of $1 v$2 has no Cargo.toml.orig"
    rm -f "$tmp/published/Cargo.lock" "$tmp/published/.cargo_vcs_info.json"
    git archive "HEAD:$3" | tar -xf - -C "$tmp/here" || die "could not read $3 at HEAD"
    diff -r -q "$tmp/published" "$tmp/here" >&2 || same=false
    rm -rf "${tmp:?}"
    [[ "$same" == true ]]
}

# A crate with a version of its own is uploaded only when that version is new. If the version
# is already on crates.io and the crate has changed since, the release would ship the other
# crates against the old code: the crate needs a new version first.
while read -r name ver dir; do
    [[ "$ver" != "-" ]] || continue
    if on_crates_io "$name" "$ver" && ! unchanged_since_published "$name" "$ver" "$dir"; then
        die "$dir has changed since $name v$ver was published: give it a new version (in $dir/Cargo.toml, and where another crate requires it)"
    fi
done <<<"$PUBLISHABLE"

echo "Current workspace version : $CURRENT_VERSION"
echo "New workspace version     : $VERSION"
echo "CHANGELOG section         : $NEW_HEADING"
echo "Crates                    : $(awk '{printf "%s%s ", $1, ($2 == "-" ? "" : " (v" $2 ")")}' <<<"$PUBLISHABLE")"
echo "Tag                       : $TAG$([[ "$PYTHON" == true ]] && echo ", $PY_TAG")"
echo "Publish crates            : $([[ "$PUBLISH" == true ]] && echo yes || echo no)"
echo "GitHub release            : $([[ "$GITHUB_RELEASE" == true ]] && echo yes || echo no)"
echo "Dry-run                   : $([[ "$DRY_RUN" == true ]] && echo yes || echo no)"
echo

# ── Edits ───────────────────────────────────────────────────────────────────────────────────

echo "Updating $ROOT_CARGO"
echo "  [workspace.package] version: $CURRENT_VERSION -> $VERSION"
while read -r name ver; do
    echo "  [workspace.dependencies] $name version: $ver -> $VERSION"
done <<<"$(internal_dep_versions)"
echo "Updating $CHANGELOG"
echo "  ${UNRELEASED_HEADING:-$DATED_HEADING} -> $NEW_HEADING"

if [[ "$DRY_RUN" == false ]]; then
    if [[ "$SAME_VERSION" == false ]]; then
        awk -v ver="$VERSION" '
            /^\[workspace\.package\]/ { pkg = 1; deps = 0; print; next }
            /^\[workspace\.dependencies\]/ { deps = 1; pkg = 0; print; next }
            /^\[/ { pkg = 0; deps = 0 }
            pkg && /^version[[:space:]]*=/ { sub(/"[^"]*"/, "\"" ver "\"") }
            deps && /path = "crates\// && /version = "/ { sub(/version = "[^"]*"/, "version = \"" ver "\"") }
            { print }
        ' "$ROOT_CARGO" >"${ROOT_CARGO}.tmp" && mv "${ROOT_CARGO}.tmp" "$ROOT_CARGO"
        [[ "$(workspace_version)" == "$VERSION" ]] || die "workspace version update failed"
        while read -r name ver; do
            [[ "$ver" == "$VERSION" ]] || die "workspace dependency $name version update failed"
        done <<<"$(internal_dep_versions)"
    fi
    if [[ -n "$UNRELEASED_HEADING" ]]; then
        awk -v old="$UNRELEASED_HEADING" -v new="$NEW_HEADING" '
            !done && $0 == old { print new; done = 1; next } { print }
        ' "$CHANGELOG" >"${CHANGELOG}.tmp" && mv "${CHANGELOG}.tmp" "$CHANGELOG"
    fi
else
    echo "Dry-run: would rewrite $ROOT_CARGO and $CHANGELOG and refresh $LOCKFILE"
fi

# Refreshes Cargo.lock for the new version.
run cargo check --workspace -q
if [[ "$DRY_RUN" == true || -n "$(git status --porcelain)" ]]; then
    run git add "$ROOT_CARGO" "$LOCKFILE" "$CHANGELOG"
    run git commit -m "chore(release): spearfish v${VERSION}"
else
    echo "Nothing to commit (workspace and CHANGELOG already at v${VERSION})"
fi

# ── Publish ─────────────────────────────────────────────────────────────────────────────────

# Crates already on crates.io at their version (published by an interrupted release or by an
# earlier attempt below, or a crate with a version of its own that has not changed) are left
# out of the check and the upload.
EXCLUDE=()
TO_PUBLISH=0
pending_crates() {
    local name ver
    EXCLUDE=()
    TO_PUBLISH=0
    while read -r name ver _; do
        if on_crates_io "$name" "$(crate_version "$ver")"; then
            EXCLUDE+=(--exclude "$name")
        else
            TO_PUBLISH=$((TO_PUBLISH + 1))
        fi
    done <<<"$PUBLISHABLE"
}

# How long to wait before trying again after a failed `cargo publish` whose output is in $1, in
# seconds; fails for errors that waiting does not fix. crates.io answers 429 when an account
# creates many new crates in a short time ("Please try again after <date> GMT"), and cargo gives
# up when a just-published dependency is slow to appear in the index.
retry_wait() {
    local after now until
    if grep -q "429 Too Many Requests" "$1"; then
        after="$(grep -o 'try again after [^.]* GMT' "$1" | head -1 | sed 's/^try again after //')"
        now="$(date +%s)"
        if [[ -n "$after" ]] && until="$(date -d "$after" +%s 2>/dev/null ||
            TZ=UTC date -j -f '%a, %d %b %Y %H:%M:%S GMT' "$after" +%s 2>/dev/null)"; then
            echo $((until > now ? until - now + 15 : 15))
        else
            echo 600
        fi
    elif grep -q -E "timed out waiting|timeout while waiting" "$1"; then
        echo 60
    else
        return 1
    fi
}

pending_crates
if [[ "$TO_PUBLISH" -eq 0 ]]; then
    echo "Every crate of v${VERSION} is already on crates.io"
else
    while read -r name ver _; do
        if [[ " ${EXCLUDE[*]-} " == *" $name "* ]]; then
            echo "$name v$(crate_version "$ver") is already on crates.io; leaving it out"
        fi
    done <<<"$PUBLISHABLE"
    # Every crate must build from its own package, as crates.io will build it.
    run cargo publish --workspace ${EXCLUDE[@]+"${EXCLUDE[@]}"} --dry-run
    if [[ "$PUBLISH" == true ]]; then
        LOG="$(mktemp)"
        for attempt in 1 2 3 4 5 6; do
            # Cargo publishes in dependency order and waits for each crate to be indexed.
            if run cargo publish --workspace ${EXCLUDE[@]+"${EXCLUDE[@]}"} 2>&1 | tee "$LOG"; then
                break
            fi
            WAIT="$(retry_wait "$LOG")" || die "cargo publish failed (see above)"
            [[ "$attempt" -lt 6 ]] || die "cargo publish still failing after $attempt attempts"
            echo "crates.io is not ready (attempt $attempt); trying again in ${WAIT}s"
            sleep "$WAIT"
            pending_crates
            [[ "$TO_PUBLISH" -gt 0 ]] || break
        done
        rm -f "$LOG"
    fi
fi

# ── Tag and push ────────────────────────────────────────────────────────────────────────────

run git tag -a "$TAG" -m "spearfish v${VERSION}"
run git push origin HEAD:main
run git push origin "$TAG"

if [[ "$GITHUB_RELEASE" == true ]]; then
    NOTES="$(mktemp)"
    trap 'rm -f "$NOTES"' EXIT
    # The section's body: from after its heading to the next `## ` heading.
    awk -v h="$NEW_HEADING" '
        $0 == h { p = 1; next } p && /^## / { exit } p { print }
    ' "$CHANGELOG" >"$NOTES"
    if [[ "$DRY_RUN" == true ]]; then
        # Not dated yet in a dry run: the section under its current heading.
        awk -v h="${UNRELEASED_HEADING:-$DATED_HEADING}" '
            $0 == h { p = 1; next } p && /^## / { exit } p { print }
        ' "$CHANGELOG" >"$NOTES"
    fi
    [[ -s "$NOTES" ]] || die "the CHANGELOG section for $VERSION is empty"
    run gh release create "$TAG" --verify-tag --title "spearfish ${VERSION}" --notes-file "$NOTES"
fi

if [[ "$PYTHON" == true ]]; then
    run git tag -a "$PY_TAG" -m "spearfish Python package v${VERSION}"
    run git push origin "$PY_TAG"
fi

# ── Next steps ──────────────────────────────────────────────────────────────────────────────

echo
if [[ "$DRY_RUN" == true ]]; then
    echo "Dry-run complete"
else
    echo "Released spearfish v${VERSION}"
    TARBALL="${REPO_URL}/archive/refs/tags/${TAG}.tar.gz"
    if command -v sha256sum >/dev/null; then SUM=(sha256sum); else SUM=(shasum -a 256); fi
    if SHA="$(curl -sfL "$TARBALL" | "${SUM[@]}" | cut -d' ' -f1)" && [[ -n "$SHA" ]]; then
        echo
        echo "For the bioconda recipe:"
        echo "  url: $TARBALL"
        echo "  sha256: $SHA"
    fi
fi
