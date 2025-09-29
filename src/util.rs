use std::ops::Deref;
use std::ops::DerefMut;
use std::ops::Index;
use std::ops::IndexMut;

use crate::fold::{Value, Values};

#[derive(Debug, Default)]
pub struct ByteStr<B>(pub B);

impl<B: AsRef<[u8]>> std::fmt::Display for ByteStr<B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.0
                .as_ref()
                .iter()
                .map(|&b| b as char)
                .collect::<String>()
        )
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

#[derive(Debug)]
pub struct Matrix<V>(pub V);

impl<V> From<V> for Matrix<V> {
    fn from(value: V) -> Self {
        Self(value)
    }
}

impl<V> Matrix<V> {
    fn get(&self, i: usize, j: usize) -> Option<&Value>
    where
        V: Deref<Target = Values>,
    {
        self.0.get(i)?.get(j)
    }

    fn get_mut(&mut self, i: usize, j: usize) -> Option<&mut Value>
    where
        V: DerefMut<Target = Values>,
    {
        self.0.get_mut(i)?.get_mut(j)
    }
}

// impl<V> Index<usize> for Matrix<V>
// where
//     V: Deref<Target = Values>,
// {
//     type Output = Vec<Value>;
//
//     fn index(&self, i: usize) -> &Self::Output {
//         self.0.get(i).expect("index out of bounds")
//     }
// }
//
// impl<V> IndexMut<usize> for Matrix<V>
// where
//     V: DerefMut<Target = Values>,
// {
//     fn index_mut(&mut self, i: usize) -> &mut Self::Output {
//         self.0.get_mut(i).expect("index out of bounds")
//     }
// }

impl<V> Index<(usize, usize)> for Matrix<V>
where
    V: Deref<Target = Values>,
{
    type Output = Value;

    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        // &self[i][j]
        self.get(i, j).expect("index out of bounds")
    }
}

impl<V> IndexMut<(usize, usize)> for Matrix<V>
where
    V: DerefMut<Target = Values>,
{
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Self::Output {
        // &mut self[i][j]
        self.get_mut(i, j).expect("index out of bounds")
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
