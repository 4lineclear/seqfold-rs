#[derive(Debug, Default)]
pub struct ByteStr<B>(pub B);

impl<B: AsRef<[u8]>> std::fmt::Display for ByteStr<B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s: String = self.0.as_ref().iter().map(|&b| b as char).collect();
        write!(f, "{s}")
    }
}

pub trait ToIsize {
    fn to_isize(self) -> isize;
}

impl ToIsize for usize {
    fn to_isize(self) -> isize {
        self as isize
    }
}

impl ToIsize for isize {
    fn to_isize(self) -> isize {
        self
    }
}

impl ToIsize for i32 {
    fn to_isize(self) -> isize {
        self as isize
    }
}

/// Round to one decimal
pub fn round1(n: f64) -> f64 {
    (n * 10.0).round() / 10.0
}

/// Round to two decimal
pub fn round2(n: f64) -> f64 {
    (n * 100.0).round() / 100.0
}
