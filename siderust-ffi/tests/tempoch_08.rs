//! Validate the upstream time ABI linked by Siderust.
use std::ptr;
use tempoch_ffi::{
    tempoch_ffi_version, tempoch_time_from_civil, tempoch_time_from_format, tempoch_time_to_civil,
    tempoch_time_to_format, TempochFormatTag, TempochScaleTag, TempochStatus, TempochTime,
    TempochUtc,
};

#[test]
fn civil_utc_day_formats_use_the_published_abi() {
    assert_eq!(tempoch_ffi_version(), 1102);
    for (format, raw) in [
        (TempochFormatTag::JD, 2457754.5),
        (TempochFormatTag::MJD, 57754.0),
    ] {
        let mut time = TempochTime {
            hi_seconds: 0.0,
            lo_seconds: 0.0,
        };
        let mut civil = TempochUtc {
            year: 0,
            month: 0,
            day: 0,
            hour: 0,
            minute: 0,
            second: 0,
            nanosecond: 0,
        };
        let mut encoded = f64::NAN;
        // All pointers address live writable values; null context selects the upstream default.
        unsafe {
            assert_eq!(
                tempoch_time_from_format(
                    raw,
                    TempochScaleTag::UTC as i32,
                    format as i32,
                    ptr::null(),
                    &mut time
                ),
                TempochStatus::Ok
            );
            assert_eq!(
                tempoch_time_to_civil(time, ptr::null(), &mut civil),
                TempochStatus::Ok
            );
            assert_eq!(
                (
                    civil.year,
                    civil.month,
                    civil.day,
                    civil.hour,
                    civil.minute,
                    civil.second
                ),
                (2017, 1, 1, 0, 0, 0)
            );
            assert_eq!(
                tempoch_time_to_format(
                    time,
                    TempochScaleTag::UTC as i32,
                    format as i32,
                    ptr::null(),
                    &mut encoded
                ),
                TempochStatus::Ok
            );
            assert!((encoded - raw).abs() < 1e-9);
            civil = TempochUtc {
                year: 2016,
                month: 12,
                day: 31,
                hour: 23,
                minute: 59,
                second: 60,
                nanosecond: 0,
            };
            assert_eq!(
                tempoch_time_from_civil(civil, ptr::null(), &mut time),
                TempochStatus::Ok
            );
            encoded = -1.0;
            assert_ne!(
                tempoch_time_to_format(
                    time,
                    TempochScaleTag::UTC as i32,
                    format as i32,
                    ptr::null(),
                    &mut encoded
                ),
                TempochStatus::Ok
            );
            assert_eq!(encoded, -1.0);
            assert_ne!(
                tempoch_time_from_format(
                    f64::INFINITY,
                    TempochScaleTag::UTC as i32,
                    format as i32,
                    ptr::null(),
                    &mut time
                ),
                TempochStatus::Ok
            );
        }
    }
}

#[test]
fn siderust_tle_ffi_reports_the_published_civil_epoch() {
    use siderust_ffi::*;
    use std::ffi::CString;
    let l1 = CString::new("1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753")
        .unwrap();
    let l2 = CString::new("2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667")
        .unwrap();
    let mut tle = ptr::null_mut();
    let mut prop = ptr::null_mut();
    assert_eq!(
        siderust_tle_parse(l1.as_ptr(), l2.as_ptr(), &mut tle),
        SiderustStatus::Ok
    );
    assert_eq!(siderust_sgp4_new(tle, 0, &mut prop), SiderustStatus::Ok);
    let mut epoch = f64::NAN;
    assert_eq!(
        siderust_sgp4_epoch_jd_utc(prop, &mut epoch),
        SiderustStatus::Ok
    );
    assert!((epoch - 2451723.28495062).abs() < 1e-9);
    let mut pos = [0.0; 3];
    let mut vel = [0.0; 3];
    assert_eq!(
        siderust_sgp4_propagate_at(prop, epoch, pos.as_mut_ptr(), vel.as_mut_ptr()),
        SiderustStatus::Ok
    );
    assert!(pos.iter().chain(&vel).all(|v| v.is_finite()));
    pos = [-123.0; 3];
    assert_eq!(
        siderust_sgp4_propagate_at(prop, f64::NAN, pos.as_mut_ptr(), vel.as_mut_ptr()),
        SiderustStatus::InvalidArgument
    );
    assert_eq!(pos, [-123.0; 3]);
    siderust_sgp4_free(prop);
    siderust_tle_free(tle);
}

#[test]
fn precise_horizontal_rejects_invalid_scale_coordinates_without_writing_output() {
    use siderust_ffi::coordinates::siderust_spherical_dir_to_horizontal_precise;
    use siderust_ffi::error::SiderustStatus;
    use siderust_ffi::types::{SiderustFrame, SiderustGeodetict, SiderustSphericalDir};
    let observer = SiderustGeodetict {
        lon_deg: 0.0,
        lat_deg: 45.0,
        height_m: 0.0,
    };
    for (tt, ut1) in [
        (f64::NAN, 2451545.0),
        (2451545.0, f64::NAN),
        (2451545.0, f64::INFINITY),
    ] {
        let mut out = SiderustSphericalDir {
            polar_deg: -123.0,
            azimuth_deg: -456.0,
            frame: SiderustFrame::ICRS,
        };
        assert_eq!(
            siderust_spherical_dir_to_horizontal_precise(
                20.0,
                30.0,
                SiderustFrame::ICRS,
                tt,
                ut1,
                observer,
                &mut out
            ),
            SiderustStatus::InvalidArgument
        );
        assert_eq!(out.polar_deg, -123.0);
        assert_eq!(out.azimuth_deg, -456.0);
        assert_eq!(out.frame, SiderustFrame::ICRS);
    }
}
