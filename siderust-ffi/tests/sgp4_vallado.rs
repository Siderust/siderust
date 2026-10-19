//! Verify SGP4 model discriminants and the shared numerical backend at the C boundary.
use siderust::astro::sgp4::{GravityModel, Sgp4Propagator};
use siderust::formats::tle::parse_tle;
use siderust_ffi::*;
use std::{ffi::CString, ptr};

const L1: &str = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
const L2: &str = "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

#[test]
fn all_ffi_models_use_the_canonical_typed_propagator() {
    let line1 = CString::new(L1).unwrap();
    let line2 = CString::new(L2).unwrap();
    let mut tle_handle = ptr::null_mut();
    assert_eq!(
        siderust_tle_parse(line1.as_ptr(), line2.as_ptr(), &mut tle_handle),
        SiderustStatus::Ok
    );
    let tle = parse_tle(L1, L2).unwrap();
    for (code, model) in [
        (0, GravityModel::Wgs72),
        (1, GravityModel::Wgs72Iau),
        (2, GravityModel::Wgs84),
    ] {
        let typed = Sgp4Propagator::from_tle_with_model(&tle, model).unwrap();
        let mut handle = ptr::null_mut();
        assert_eq!(
            siderust_sgp4_new(tle_handle, code, &mut handle),
            SiderustStatus::Ok
        );
        let mut actual_model = -1;
        assert_eq!(
            siderust_sgp4_gravity_model(handle, &mut actual_model),
            SiderustStatus::Ok
        );
        assert_eq!(actual_model, code);
        let mut epoch = 0.0;
        assert_eq!(
            siderust_sgp4_epoch_jd_utc(handle, &mut epoch),
            SiderustStatus::Ok
        );
        for days in [-1.0, 0.0, 1.0] {
            let raw = epoch + days;
            let target = tempoch::JulianDate::<tempoch::UTC>::try_new_with(
                qtty::Day::new(raw),
                typed.time_context(),
            )
            .unwrap();
            let expected = typed.propagate_at(target).unwrap();
            let mut position = [0.0; 3];
            let mut velocity = [0.0; 3];
            assert_eq!(
                siderust_sgp4_propagate_at(
                    handle,
                    raw,
                    position.as_mut_ptr(),
                    velocity.as_mut_ptr()
                ),
                SiderustStatus::Ok
            );
            for i in 0..3 {
                assert_eq!(position[i], expected.position().as_array()[i].value());
                assert_eq!(velocity[i], expected.velocity().as_array()[i].value());
            }
        }
        siderust_sgp4_free(handle);
    }
    let mut handle = ptr::null_mut();
    // Preserve the established C ABI fallback for unrecognized model codes.
    assert_eq!(
        siderust_sgp4_new(tle_handle, 3, &mut handle),
        SiderustStatus::Ok
    );
    let mut model = -1;
    assert_eq!(
        siderust_sgp4_gravity_model(handle, &mut model),
        SiderustStatus::Ok
    );
    assert_eq!(model, 0);
    siderust_sgp4_free(handle);
    siderust_tle_free(tle_handle);
}

#[test]
fn divergent_trajectory_maps_to_invalid_argument_without_writing_outputs() {
    let l1 = CString::new("1 33333U 05037B   05333.02012661  .25992681  00000-0  24476-3 0  1532")
        .unwrap();
    let l2 = CString::new("2 33333  96.4736 157.9986 9950000 244.0492 110.6523  4.00004038 10700")
        .unwrap();
    let mut tle = ptr::null_mut();
    let mut prop = ptr::null_mut();
    assert_eq!(
        siderust_tle_parse(l1.as_ptr(), l2.as_ptr(), &mut tle),
        SiderustStatus::Ok
    );
    assert_eq!(siderust_sgp4_new(tle, 0, &mut prop), SiderustStatus::Ok);
    let mut epoch = 0.0;
    assert_eq!(
        siderust_sgp4_epoch_jd_utc(prop, &mut epoch),
        SiderustStatus::Ok
    );
    let mut position = [-123.0; 3];
    let mut velocity = [-456.0; 3];
    assert_eq!(
        siderust_sgp4_propagate_at(
            prop,
            epoch + 25.0 / 1440.0,
            position.as_mut_ptr(),
            velocity.as_mut_ptr()
        ),
        SiderustStatus::InvalidArgument
    );
    assert_eq!(position, [-123.0; 3]);
    assert_eq!(velocity, [-456.0; 3]);
    siderust_sgp4_free(prop);
    siderust_tle_free(tle);
}
