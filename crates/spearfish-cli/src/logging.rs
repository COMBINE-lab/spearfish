//! A minimal stderr logger for the `log` facade: `-q` shows errors only, the default adds
//! warnings and progress summaries (info), `-v` debug and `-vv` trace messages.
//! [`warn_always`] writes the few warnings about the inputs that change results, which `-q`
//! does not hide.

use std::io::IsTerminal;
use std::sync::OnceLock;

use log::{Level, LevelFilter, Log, Metadata, Record};

struct Stderr {
    color: bool,
}

impl Log for Stderr {
    fn enabled(&self, m: &Metadata) -> bool {
        m.level() <= log::max_level()
    }

    fn log(&self, r: &Record) {
        if !self.enabled(r.metadata()) {
            return;
        }
        let (tag, code) = match r.level() {
            Level::Error => ("error", "31"),
            Level::Warn => ("warning", "33"),
            Level::Info => return eprintln!("{}", r.args()),
            Level::Debug => ("debug", "2"),
            Level::Trace => ("trace", "2"),
        };
        emit(self.color, tag, code, r.args());
    }

    fn flush(&self) {}
}

fn emit(color: bool, tag: &str, code: &str, msg: &std::fmt::Arguments) {
    if color {
        eprintln!("\x1b[{code}m{tag}\x1b[0m: {msg}");
    } else {
        eprintln!("{tag}: {msg}");
    }
}

static LOGGER: OnceLock<Stderr> = OnceLock::new();

/// A warning shown at every level, `-q` included: for a problem with the inputs that changes
/// the results and is fixed by rebuilding or rerunning (not for advice or ignored options,
/// which go through `log::warn!`).
pub fn warn_always(msg: std::fmt::Arguments) {
    let color = LOGGER.get().is_some_and(|l| l.color);
    emit(color, "warning", "33", &msg);
}

/// Install the logger at the level `-q`/`-v` select.
pub fn init(quiet: bool, verbose: u8) {
    let level = match (quiet, verbose) {
        (true, _) => LevelFilter::Error,
        (false, 0) => LevelFilter::Info,
        (false, 1) => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };
    let logger = LOGGER.get_or_init(|| Stderr {
        color: std::io::stderr().is_terminal() && std::env::var_os("NO_COLOR").is_none(),
    });
    if log::set_logger(logger).is_ok() {
        log::set_max_level(level);
    }
}

/// Whether progress lines should be drawn: info level, and stderr is a terminal.
pub fn progress_enabled() -> bool {
    log::max_level() >= LevelFilter::Info && std::io::stderr().is_terminal()
}
