//! Numerical core, free of extendr types so that it can be tested with
//! `cargo test` and shared by the panel cointegration tests of the package.

pub mod adf;
pub mod banerjee;
pub mod bootstrap;
pub mod breaks;
pub mod design;
pub mod ecm;
pub mod eigen;
pub mod factors;
pub mod linalg;
pub mod lrv;
pub mod mq;
pub mod nulldist;
pub mod panel;
pub mod parallel;
pub mod rng;
pub mod westerlund;
