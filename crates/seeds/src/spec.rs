//! Seed-scheme spec strings: `name:key=value,...`.
//!
//! - `lcp:K=7`
//! - `minimizer:k=15,w=10`
//! - `closed-syncmer:k=15,s=5[,keep=0.5]`
//! - `open-syncmer:k=15,s=9,t=3[,keep=0.5]`
//! - `randstrobe:k=15,s=9[,t=3][,keep=1],wmin=1,wmax=4` (open syncmers if `t` is given,
//!   closed otherwise)

use std::collections::HashMap;
use std::fmt;

use lcp_core::LcpParser;

use crate::{LcpSeeds, Minimizer, Randstrobe, SeedScheme, Syncmer, SyncmerKind};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpecError(pub String);

impl fmt::Display for SpecError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "bad seed spec: {}", self.0)
    }
}

impl std::error::Error for SpecError {}

struct Params<'a> {
    spec: &'a str,
    kv: HashMap<&'a str, &'a str>,
}

impl Params<'_> {
    fn get<T: std::str::FromStr>(&self, key: &str) -> Result<Option<T>, SpecError> {
        self.kv
            .get(key)
            .map(|v| {
                v.parse()
                    .map_err(|_| SpecError(format!("{}: invalid {key}={v}", self.spec)))
            })
            .transpose()
    }

    fn req<T: std::str::FromStr>(&self, key: &str) -> Result<T, SpecError> {
        self.get(key)?
            .ok_or_else(|| SpecError(format!("{}: missing {key}", self.spec)))
    }
}

pub fn from_spec(spec: &str) -> Result<Box<dyn SeedScheme>, SpecError> {
    let (name, rest) = spec.split_once(':').unwrap_or((spec, ""));
    let mut kv = HashMap::new();
    for part in rest.split(',').filter(|p| !p.is_empty()) {
        let (k, v) = part
            .split_once('=')
            .ok_or_else(|| SpecError(format!("{spec}: expected key=value, got {part}")))?;
        kv.insert(k, v);
    }
    let p = Params { spec, kv };
    let check_k = |k: usize| {
        if (1..=32).contains(&k) {
            Ok(k)
        } else {
            Err(SpecError(format!("{spec}: k must be in 1..=32")))
        }
    };
    let syncmer = |kind: SyncmerKind| -> Result<Syncmer, SpecError> {
        let k = check_k(p.req("k")?)?;
        let s: usize = p.req("s")?;
        if s == 0 || s > k {
            return Err(SpecError(format!("{spec}: need 1 <= s <= k")));
        }
        if let SyncmerKind::Open { t } = kind
            && t > k - s
        {
            return Err(SpecError(format!("{spec}: need t <= k - s")));
        }
        let keep = p.get("keep")?.unwrap_or(1.0);
        if !(keep > 0.0 && keep <= 1.0) {
            return Err(SpecError(format!("{spec}: keep must be in (0, 1]")));
        }
        Ok(Syncmer { k, s, kind, keep })
    };
    Ok(match name {
        "lcp" => Box::new(LcpSeeds {
            parser: LcpParser::default(),
            max_level: p.req("K")?,
        }),
        "minimizer" => Box::new(Minimizer {
            k: check_k(p.req("k")?)?,
            w: p.req("w")?,
        }),
        "closed-syncmer" => Box::new(syncmer(SyncmerKind::Closed)?),
        "open-syncmer" => Box::new(syncmer(SyncmerKind::Open { t: p.req("t")? })?),
        "randstrobe" => {
            let kind = match p.get("t")? {
                Some(t) => SyncmerKind::Open { t },
                None => SyncmerKind::Closed,
            };
            let (w_min, w_max): (usize, usize) = (p.req("wmin")?, p.req("wmax")?);
            if w_min == 0 || w_max < w_min {
                return Err(SpecError(format!("{spec}: need 1 <= wmin <= wmax")));
            }
            Box::new(Randstrobe {
                base: syncmer(kind)?,
                w_min,
                w_max,
            })
        }
        _ => return Err(SpecError(format!("{spec}: unknown scheme {name}"))),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_names() {
        for s in [
            "lcp:K=7",
            "minimizer:k=15,w=10",
            "closed-syncmer:k=15,s=5",
            "closed-syncmer:k=15,s=5,keep=0.25",
            "open-syncmer:t=3,k=15,s=9",
            "randstrobe:t=3,k=15,s=9,wmin=1,wmax=4",
            "randstrobe:k=15,s=5,wmin=2,wmax=6",
        ] {
            let scheme = from_spec(s).unwrap();
            assert_eq!(from_spec(&scheme.name()).unwrap().name(), scheme.name(), "{s}");
        }
    }

    #[test]
    fn rejects_bad_specs() {
        for s in [
            "minimizer:k=40,w=10",
            "closed-syncmer:k=15,s=20",
            "foo:k=1",
            "minimizer:k=15",
            "randstrobe:k=15,s=5,wmin=3,wmax=2",
        ] {
            assert!(from_spec(s).is_err(), "{s}");
        }
    }
}
