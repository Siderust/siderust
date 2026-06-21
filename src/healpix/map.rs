use crate::coordinates::frames::ReferenceFrame;
use crate::healpix::{HealpixError, HealpixGrid, Result};
use std::marker::PhantomData;

#[derive(Debug, Clone, PartialEq)]
pub struct HealpixMap<F, T>
where F: ReferenceFrame,
{
    grid: HealpixGrid,
    values: Vec<T>,
    pub(crate) marker: PhantomData<F>,
}

impl<F, T> HealpixMap<F, T>
where F: ReferenceFrame,
{
    pub fn new(grid: HealpixGrid, values: Vec<T>) -> Result<Self> {
        let len = values.len();
        let npix = grid.npix();
        if u64::try_from(len).expect("usize length fits u64") != npix {
            return Err(HealpixError::MapLengthMismatch { len, npix });
        }
        Ok(Self { grid, values, marker: PhantomData })
    }
    #[must_use]
    pub fn grid(&self) -> HealpixGrid { self.grid }
    #[must_use]
    pub fn values(&self) -> &[T] { &self.values }
    #[must_use]
    pub fn values_mut(&mut self) -> &mut [T] { &mut self.values }
    #[must_use]
    pub fn into_values(self) -> Vec<T> { self.values }
}
