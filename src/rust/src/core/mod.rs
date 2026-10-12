//! Numerical core, free of extendr types so that it can be tested with
//! `cargo test` and shared by the panel cointegration tests of the package.

pub mod bootstrap;
pub mod design;
pub mod ecm;
pub mod linalg;
pub mod lrv;
pub mod panel;
pub mod westerlund;
