"""Put the sdist's license files where PyPI looks for them, and check they are there.

maturin builds the sdist of a crate inside a Cargo workspace with the workspace layout: the
package lives at `<root>/crates/py-spearfish/`, LICENSE with it. PKG-INFO declares
`License-File: LICENSE`, and PyPI looks for that path at the archive root, so it rejects the
upload ("License-File LICENSE does not exist in distribution file"). This copies each declared
license file from the package directory to the root, then checks that every declared file is
present.

Usage: python3 sdist_license.py DIST_DIR   (rewrites each *.tar.gz in place)
"""

import email.parser
import gzip
import io
import sys
import tarfile
from pathlib import Path

PACKAGE_DIR = "crates/py-spearfish"


def fix(path: Path) -> None:
    with tarfile.open(path, "r:gz") as tf:
        members = [(m, tf.extractfile(m).read() if m.isfile() else None) for m in tf.getmembers()]
    names = {m.name for m, _ in members}
    root = members[0][0].name.split("/")[0]
    pkg_info = next(data for m, data in members if m.name == f"{root}/PKG-INFO")
    declared = email.parser.Parser().parsestr(pkg_info.decode()).get_all("License-File") or []
    added = []
    for lf in declared:
        dest = f"{root}/{lf}"
        if dest in names:
            continue
        src = f"{root}/{PACKAGE_DIR}/{lf}"
        data = next((d for m, d in members if m.name == src and m.isfile()), None)
        if data is None:
            sys.exit(f"{path.name}: License-File {lf} is neither at the root nor in {PACKAGE_DIR}")
        info = tarfile.TarInfo(dest)
        info.size, info.mode = len(data), 0o644
        info.mtime = next(m.mtime for m, _ in members if m.name == src)
        members.append((info, data))
        names.add(dest)
        added.append(lf)
    missing = [lf for lf in declared if f"{root}/{lf}" not in names]
    if missing:
        sys.exit(f"{path.name}: declared license files missing: {missing}")
    if added:
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
            with tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as out:
                for m, data in members:
                    out.addfile(m, io.BytesIO(data) if data is not None else None)
        path.write_bytes(buf.getvalue())
    print(f"{path.name}: license files {declared} present" + (f" (added {added} at the root)" if added else ""))


if __name__ == "__main__":
    sdists = sorted(Path(sys.argv[1]).glob("*.tar.gz"))
    if not sdists:
        sys.exit(f"no sdist in {sys.argv[1]}")
    for p in sdists:
        fix(p)
