pub trait Kind {}
#[derive(Debug, Copy, Clone)]
pub struct PositionKind;
#[derive(Debug, Copy, Clone)]
pub struct DirectionKind; // ||v|| == 1
impl Kind for PositionKind {}
impl Kind for DirectionKind {}
