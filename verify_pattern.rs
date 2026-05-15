// Quick pattern verification - not a real file, just checking the API
// jd.to_chrono()  → jd.to::<UTC>().to_chrono()
// JulianDate::from_chrono(dt) → Time::<UTC>::from_chrono(dt).to::<TT>().to::<JD>()
// ModifiedJulianDate::from(jd) → jd.to::<MJD>()
