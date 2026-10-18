//! A minimal stderr logger for the `log` facade: `-q` shows errors only, the default adds
//! warnings and progress summaries (info), `-v` debug and `-vv` trace messages.

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
        if self.color {
            eprintln!("\x1b[{code}m{tag}\x1b[0m: {}", r.args());
        } else {
            eprintln!("{tag}: {}", r.args());
        }
    }

    fn flush(&self) {}
}

static LOGGER: OnceLock<Stderr> = OnceLock::new();

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
