pub trait Kind {}
pub struct PositionKind;
pub struct DirectionKind; // ||v|| == 1
impl Kind for PositionKind {}
impl Kind for DirectionKind {}
