#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Vallés Puig, Ramon
"""Development-only independent CSPICE fixture generator; never used by Rust builds.

Requires SpiceyPy 8.3.0 and NumPy. Run at repository root with KERNEL_DIR set.
"""
# Independent reference generation: Python 3 + spiceypy 8.3.0 (CSPICE N0067).
# Run from repository root with KERNEL_DIR set to the five downloaded kernels.
import os, pathlib, struct, csv, hashlib, json, shutil
import numpy as np
import spiceypy as spice
source = pathlib.Path(os.environ['KERNEL_DIR'])
out = pathlib.Path('tests/fixtures/lunar_de440')
out.mkdir(parents=True, exist_ok=True)
files = ['de440s.bsp','moon_pa_de440_200625.bpc','moon_de440_220930.tf','pck00011.tpc','naif0012.tls']
for f in files: spice.furnsh(str(source/f))
# Five original BPC records. Preserve exact DE440 coefficients with PCKW02.
b = (source/files[1]).read_bytes()
r = struct.unpack_from('<i', b, 76)[0]
first,last,cls,base,typ,begin,finish = struct.unpack_from('<2d5i', b, (r-1)*1024+24)
init,interval,size,count = struct.unpack_from('<4d', b, (finish-4)*8)
k = int((820000000-init)//interval)
start = init+k*interval
end = start+4*interval
for name in ['excerpt.bpc','excerpt.bsp']:
    (out/name).unlink(missing_ok=True)
h = spice.pckopn(str(out/'excerpt.bpc'),'DE440 original coefficients excerpt',0)
records = np.array(struct.unpack_from('<'+str(5*int(size))+'d', b, (begin-1+k*int(size))*8)).reshape(5,int(size))
cdata = records[:,2:].flatten()
spice.pckw02(h,cls,'J2000',start,end,'DE440 five original records',interval,5,int((size-2)/3)-1,cdata,start)
spice.pckcls(h)
h = spice.spkopn(str(out/'excerpt.bsp'),'DE440s original segment excerpt',0)
for target in [10,3,301,399]:
    handle,desc,ident = spice.spksfs(target,(start+end)/2,80)
    spice.spksub(handle,desc,ident,start,end,h)
spice.spkcls(h)
for f in [files[2], files[4]]: shutil.copyfile(source/f,out/f)
# A deliberately large three-axis fixed rotation exercises composition order.
rotated = """KPL/FK
\\begindata
FRAME_TEST_LUNAR_ROTATED = 139999
FRAME_139999_NAME = 'TEST_LUNAR_ROTATED'
FRAME_139999_CLASS = 4
FRAME_139999_CLASS_ID = 139999
FRAME_139999_CENTER = 301
TKFRAME_139999_SPEC = 'ANGLES'
TKFRAME_139999_RELATIVE = 'MOON_ME'
TKFRAME_139999_ANGLES = ( 30 50 20 )
TKFRAME_139999_AXES = ( 1 2 3 )
TKFRAME_139999_UNITS = 'DEGREES'
"""
(out/'test_rotation.tf').write_text(rotated)
spice.furnsh(str(out/'test_rotation.tf'))
# References always query the ORIGINAL kernels, not the excerpts or Siderust.
sites = [(0,0,1737.4),(45,30,1738.0),(0,90,1737.4),(70,-90,1737.4),
         (32.3,-89,1738.0),(25.2,-80.13,1736.0),(179.999,-45,1737.4),(-179.999,45,1737.4)]
header = ['et','lon','lat','radius_km','target','frame']+[f'm{i}' for i in range(9)]+[f'p{i}' for i in range(3)]+[f'v{i}' for i in range(3)]+['alt','az','range','above']
def row(et,lon,lat,radius,target,frame):
    rot=spice.pxform('J2000',frame,et)
    pos=spice.spkpos(str(target),et,frame,'NONE','301')[0]
    site=spice.latrec(radius,np.deg2rad(lon),np.deg2rad(lat))
    rel=pos-site
    l,p=np.deg2rad([lon,lat]); east=np.array([-np.sin(l),np.cos(l),0]); up=site/radius; north=np.cross(up,east)
    e,n,u=np.dot(rel,east),np.dot(rel,north),np.dot(rel,up)
    alt=np.rad2deg(np.arctan2(u,np.hypot(n,e))); az=np.rad2deg(np.arctan2(e,n))%360
    return [et,lon,lat,radius,target,spice.namfrm(frame),*rot.flatten(),*pos,*rel,alt,az,np.linalg.norm(rel),int(alt>0)]
with (out/'reference.csv').open('w') as f:
    w=csv.writer(f);w.writerow(header)
    for et in np.linspace(start,end,17):
        for lon,lat,radius in sites:
            for target in [10,399]: w.writerow(row(et,lon,lat,radius,target,'MOON_ME'))
    # Explicit Sun crossings of the horizon by bisection, just above/below.
    for lon,lat,radius in [sites[0],sites[4]]:
        ts=np.linspace(start,end,257)
        for a,z in zip(ts[:-1],ts[1:]):
            if row(a,lon,lat,radius,10,'MOON_ME')[-4]*row(z,lon,lat,radius,10,'MOON_ME')[-4]>=0:continue
            for _ in range(45):
                mid=(a+z)/2
                if row(a,lon,lat,radius,10,'MOON_ME')[-4]*row(mid,lon,lat,radius,10,'MOON_ME')[-4]<=0:z=mid
                else:a=mid
            for et in [(a+z)/2-1,(a+z)/2+1]:w.writerow(row(et,lon,lat,radius,10,'MOON_ME'))
    for et in [start,end]:
        for frame in ['MOON_PA_DE440','MOON_ME_DE440_ME421','TEST_LUNAR_ROTATED']:
            w.writerow(row(et,45,30,1737.4,399,frame))
# Broad external-kernel validation, plus actual BPC coverage edges for rotations.
with (out/'external_reference.csv').open('w') as f:
    w=csv.writer(f);w.writerow(header)
    for et in [-4733596800.0,0,717768000.0,1700000000.0,4735368000.0]:
        for target in [10,399]: w.writerow(row(et,32.3,-89,1738.0,target,'MOON_ME'))
with (out/'coverage_rotations.csv').open('w') as f:
    w=csv.writer(f)
    for et in [first,first+1,last-1,last]:w.writerow([et,*spice.pxform('J2000','MOON_PA_DE440',et).flatten()])
manifest={'generator':'SpiceyPy '+spice.__version__,'toolkit':spice.tkvrsn('TOOLKIT'),'abcorr':'NONE','excerpt_et':[start,end],
          'source_sha256':{f:hashlib.sha256((source/f).read_bytes()).hexdigest() for f in files}}
(out/'provenance.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(manifest)
# Time references separately quantify tempoch vs LSK DELTET approximations.
with (out/'time_reference.csv').open('w') as f:
    w=csv.writer(f)
    for utc in ['1972-07-01T00:00:00Z','1999-01-01T00:00:00Z','2016-12-31T23:59:59Z','2016-12-31T23:59:60Z','2017-01-01T00:00:00Z','2026-10-29T12:00:00Z']:
        et=spice.str2et(utc)
        w.writerow([utc,spice.unitim(et,'TDB','TDT'),et])
