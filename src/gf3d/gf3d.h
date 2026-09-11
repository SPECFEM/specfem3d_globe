/*
 *=====================================================================
 *
 *                       S p e c f e m 3 D  G l o b e
 *                       ----------------------------
 *
 *     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 *                        Princeton University, USA
 *                and CNRS / University of Marseille, France
 *                 (there are currently many more authors!)
 * (c) Princeton University and CNRS / University of Marseille, April 2014
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *=====================================================================
 */

/*
 * libgf3d -- Green function extraction from a specfem3d_globe reciprocal
 * database, as a C-callable library.
 *
 * The Fortran behind every function here is the same code bin/xgf3d runs;
 * this is a facade over it (src/gf3d/gf3d_capi.F90), not a second
 * implementation. The Python package in utils/green_function/gf3d/ is a
 * ctypes wrapper over exactly these calls.
 *
 * Linking
 * -------
 *     cc myprog.c -I<repo>/include -L<repo>/lib -lgf3d <HDF5 link line>
 *
 * or dlopen()/ctypes.CDLL() on <repo>/lib/libgf3d.so, which has the HDF5
 * library directories baked in as an rpath. Link C code with the Fortran
 * driver (or add -lgfortran / -lifcore) so that the Fortran runtime comes
 * along.
 *
 * Conventions
 * -----------
 * Return value    every function returns a GF_* code; GF_OK (0) is success.
 *                 Nothing here can abort, exit or stop the calling process:
 *                 that is the point of the facade, since a Fortran stop
 *                 inside a shared object would kill a Python interpreter
 *                 with no traceback.
 * Handles         gf3d_open() returns a small positive integer, not a
 *                 pointer. A stale or forged handle is rejected with
 *                 GF_ERR_ARG rather than dereferenced.
 * Output arrays   the caller allocates every one of them, sized from
 *                 gf3d_get_plan() and gf3d_get_info(). Nothing is
 *                 allocated on the library side and handed back.
 * Array order     C order, i.e. the last index varies fastest:
 *                     seis[ista][icomp][it]
 *                     dp[ista][ip][icomp][it]
 *                 with icomp 0,1,2 = N,E,Z.
 * Indices         0-based (ista, icomp, ip), as C expects. The Fortran
 *                 underneath is 1-based; the facade converts.
 * Units           displacement in metres, time in seconds, latitude and
 *                 longitude in degrees, depth in km (source) or metres
 *                 (station burial), moment tensor in dyne-cm.
 * Time origin     t = 0 is the *centroid* time, i.e. the source file's
 *                 origin time plus its time shift. The time shift is
 *                 carried in gf3d_source.time_shift as metadata and is not
 *                 in the trace. (This is get_cmt()'s convention, and
 *                 getting it backwards puts a trace 29 s out against real
 *                 data for the shipped example.)
 * Threads         not thread-safe. The library keeps process-wide state:
 *                 the last error message, the kd-tree that serves whichever
 *                 database located last, and specfem's shared_parameters.
 *                 Serialise calls, as the Python wrapper does.
 * Two databases   supported: each handle owns its own metadata. But the
 *                 search tree is rebuilt whenever a locate switches
 *                 database, so alternating between two of them is slow, and
 *                 two databases of *different planets or topography grids*
 *                 must not be used alternately at all.
 */

#ifndef GF3D_H
#define GF3D_H

#ifdef __cplusplus
extern "C" {
#endif

/* bumped when anything below changes incompatibly */
#define GF3D_API_VERSION 1

/* fits 'NET.STA': MAX_LENGTH_NETWORK_NAME + 1 + MAX_LENGTH_STATION_NAME */
#define GF3D_STRLEN 64

/* the number of databases that may be open at once */
#define GF3D_MAX_HANDLES 32

/* a Morton code as 16 hex digits, plus room for a terminator */
#define GF3D_MORTON_STRLEN 24

/* an open database. Valid values are 1..GF3D_MAX_HANDLES. */
typedef int gf3d_handle;

/* error codes; these are gf_par.F90's GF_* values, unchanged */
enum gf3d_status {
  GF_OK             = 0,
  GF_ERR_NO_HDF5    = 1,   /* built without HDF5 */
  GF_ERR_NO_PATH    = 2,   /* database directory not found */
  GF_ERR_NO_FILE    = 3,
  GF_ERR_HDF5       = 4,
  GF_ERR_IO         = 5,
  GF_ERR_FORMAT     = 6,   /* malformed database */
  GF_ERR_MISMATCH   = 7,   /* incompatible database */
  GF_ERR_INCOMPLETE = 8,   /* database missing element/station files */
  GF_ERR_ALLOC      = 9,
  GF_ERR_ARG        = 10,  /* invalid argument, including a bad handle */
  GF_ERR_NO_ELEMENT = 11,  /* the point is not inside the database */
  GF_ERR_GEOMETRY   = 12   /* degenerate element geometry */
};

/* source kinds */
enum gf3d_source_type {
  GF_SRC_FORCE = 1,
  GF_SRC_CMT   = 2
};

/* fixed sizes */
enum {
  GF_NCOMP   = 3,   /* N, E, Z */
  GF_NDP_MT  = 6,   /* itypsokern = 1: the six moment-tensor partials */
  GF_NDP_LOC = 10   /* itypsokern = 2: those plus lat, lon, depth, time */
};

/*
 * A seismic source, of either kind.
 *
 * For GF_SRC_CMT the fields used are latitude, longitude, depth_km, hdur,
 * time_shift and moment; for GF_SRC_FORCE they are latitude, longitude,
 * depth_km, hdur (= the FORCESOLUTION's f0), time_shift, force_stf,
 * force_factor and force_dir. Set the rest to zero.
 *
 * hdur is the file's own field, unconverted: a CMTSOLUTION's *triangle*
 * half duration (the Gaussian width specfem uses is hdur/1.628, applied
 * inside), or a FORCESOLUTION's f0. The library applies specfem's own
 * clamps, so a zero hdur becomes 5*dt exactly as get_cmt() would make it.
 */
typedef struct {
  int    source_type;      /* GF_SRC_CMT or GF_SRC_FORCE */
  int    force_stf;        /* force only: 0 Gaussian, 1 Ricker, 2 step,
                              3 monochromatic, 4 Gaussian (Meschede) */
  double latitude;         /* degrees */
  double longitude;        /* degrees */
  double depth_km;         /* km below the surface */
  double hdur;             /* s (CMT: triangle half duration; force: f0) */
  double time_shift;       /* s; metadata, not part of the trace */
  double moment[6];        /* Mrr,Mtt,Mpp,Mrt,Mrp,Mtp in dyne-cm; CMT only */
  double force_factor;     /* N; force only */
  double force_dir[3];     /* E, N, Z-up; any length, force only */
} gf3d_source;

/* what a database says about itself */
typedef struct {
  int    nelem;            /* elements in the index */
  int    nstations;
  int    nstep;            /* solver time steps */
  int    nt_subsampled;    /* stored samples = nstep / subsample_step */
  int    subsample_step;
  int    ngllx, nglly, ngllz;
  int    topography;       /* 0 or 1 */
  int    ellipticity;
  int    rotation;
  int    attenuation;
  int    gravity;
  int    pad_;             /* keeps the doubles 8-byte aligned explicitly */
  double dt;               /* solver time step, s */
  double t0;               /* first stored sample, s before the origin */
  double r_planet;         /* m; the effective value, see below */
  double rhoav;            /* kg/m^3; likewise */
  double scale_displ;
} gf3d_info;
/*
 * r_planet and rhoav are the values the library is *using*, which come from
 * the database when it records them and from specfem's Earth defaults when
 * it does not. The current writer stores neither RHOAV nor the flattening,
 * so on the shipped example databases rhoav is the Earth default.
 */

typedef struct {
  char   id[GF3D_STRLEN];        /* 'NET.STA' */
  char   network[GF3D_STRLEN];
  char   station[GF3D_STRLEN];
  double latitude;               /* degrees */
  double longitude;              /* degrees */
  double depth_m;                /* burial, metres below the surface */
  double hdur;                   /* Gaussian width of the reciprocal run, s */
  double f_cutoff;               /* Hz */
  double factor_force_source;    /* non-dimensional */
  double time_shift;             /* s */
} gf3d_station;

/* where a source sits in the mesh */
typedef struct {
  int    ielem;                          /* 1..nelem, 0 if unset */
  char   morton_hex[GF3D_MORTON_STRLEN];
  double xi, eta, gamma;                 /* element coordinates */
  double xyz[3];                         /* mapped position, non-dimensional */
  double xyz_target[3];                  /* requested position */
  double distance_km;                    /* |mapped - target| */
  double anchor_err;                     /* 27-anchor reconstruction residual */
  double theta;                          /* geocentric colatitude, radians */
  double phi;                            /* longitude, radians */
  double r_surface;                      /* surface radius above the source */
} gf3d_location;

/*
 * The output time axis and the source time function conversion, decided
 * before any element is read. nt is the length every output array needs.
 *
 * The axis is the stored one extended to the left by whole samples:
 *     t[i] = t_first + i*dt_sub,   i = 0..nt-1
 * so t_first is at or just before the requested -t0, never after it.
 */
typedef struct {
  int    nt;               /* samples on the output axis */
  int    nt_db;            /* samples the database stores */
  int    npad;             /* samples prepended: nt = nt_db + npad */
  int    subsample_step;
  int    khalf;            /* conversion kernel half length, in samples */
  int    guard;            /* 1 when the requested source is narrower than
                              the database's own, so no width correction is
                              possible; the trace is still usable, but its
                              duration is the database's */
  int    kind_stf;         /* 0 none, 1 Gaussian, 2 Heaviside */
  int    pad_;
  double dt;               /* solver step, s */
  double dt_sub;           /* stored sample spacing = dt*subsample_step */
  double t0_db;
  double t0_req;           /* the start time asked for */
  double t0;
  double t_first;          /* t[0] */
  double hdur_src;         /* the source's own field */
  double hdur_target;      /* the Gaussian width specfem would use */
  double hdur_db;          /* the database's Gaussian width */
  double hdur_corr;        /* sqrt(hdur_target^2 - hdur_db^2) */
  double trunc;            /* kernel truncation, in units of hdur_corr */
} gf3d_plan;

/* ------------------------------------------------------------------ */
/* version, sizes and errors                                          */
/* ------------------------------------------------------------------ */

/* library version string, e.g. "0.1.0". Truncated to buflen, always
   null-terminated. */
int gf3d_version(char *buf, int buflen);

/* c_sizeof of each struct above, so that a binding written against this
   header can check at load time that it agrees with the library it found.
   Any pointer may be NULL. */
int gf3d_sizeof(int *source, int *info, int *station, int *location, int *plan);

/* the message from the most recent failure, process-wide */
int gf3d_last_error(char *buf, int buflen);

/* a short name for a status code, e.g. "invalid argument" */
int gf3d_error_string(int code, char *buf, int buflen);

/* ------------------------------------------------------------------ */
/* opening and interrogating a database                               */
/* ------------------------------------------------------------------ */

/*
 * Open the database rooted at `path` (the directory holding mesh_info.h5,
 * elements/ and stations/) and write its handle to *h.
 *
 * check_completion != 0 verifies that every element/station file the index
 * promises is present, which costs one stat() per file; xgf3d --info does
 * this, extraction does not.
 */
int gf3d_open(const char *path, int check_completion, gf3d_handle *h);

/* Close a database and release the search tree. Idempotent only in the
   sense that a second call returns GF_ERR_ARG; it is never an abort. */
int gf3d_close(gf3d_handle h);

int gf3d_get_info(gf3d_handle h, gf3d_info *info);

/* ista is 0-based, 0 .. info.nstations-1 */
int gf3d_get_station(gf3d_handle h, int ista, gf3d_station *sta);

/* ------------------------------------------------------------------ */
/* locating and planning                                              */
/* ------------------------------------------------------------------ */

/*
 * Find the element containing a geographic point, and the position in it.
 * Returns GF_ERR_NO_ELEMENT when the point is outside the database.
 *
 * This also loads the topography grid on first use (about 58 MB), which is
 * why it is worth doing once rather than per extraction if the source does
 * not move.
 */
int gf3d_locate(gf3d_handle h, double lat, double lon, double depth_km,
                gf3d_location *loc);

/*
 * The output axis and source time function conversion for a source.
 *
 * t0_req is the requested start time in seconds before the centroid time;
 * pass a negative value for specfem's own rule for a forward run (1.5*hdur
 * for a moment tensor). Call this first: plan.nt is how long the output
 * arrays must be.
 */
int gf3d_get_plan(gf3d_handle h, const gf3d_source *src, double t0_req,
                  gf3d_plan *plan);

/* number of partials for itypsokern = 0, 1 or 2: 0, GF_NDP_MT, GF_NDP_LOC */
int gf3d_ndp(int itypsokern, int *ndp);

/* name ("Mrr" .. "tim") and unit ("m/dyne-cm", "m/deg", "m/km", "m/s") of
   partial ip, 0-based. Either buffer may be NULL. */
int gf3d_partial_name(int ip, char *name, int namelen, char *unit, int unitlen);

/* ------------------------------------------------------------------ */
/* extraction                                                         */
/* ------------------------------------------------------------------ */

/*
 * Seismograms at every station.
 *
 *   nt      must equal plan.nt from gf3d_get_plan() with the same source
 *           and t0_req
 *   seis    [nstations][3][nt], metres, N/E/Z          (caller-allocated)
 *   t       [nt], seconds relative to the centroid time
 *   onset   [nstations]; the amplitude just before the record starts,
 *           relative to the trace peak. Anything above ~1e-3 means the
 *           conversion kernel is reaching back past the beginning of the
 *           stored record and the first arrivals may be contaminated.
 *   loc     may be NULL; otherwise receives where the source was located
 */
int gf3d_seismograms(gf3d_handle h, const gf3d_source *src, double t0_req,
                     int nt, double *seis, double *t, double *onset,
                     gf3d_location *loc);

/*
 * Seismograms and their partial derivatives, for a moment-tensor source.
 *
 *   itypsokern  1 for the six moment-tensor partials, 2 for those plus
 *               d/d(latitude), d/d(longitude), d/d(depth), d/d(time)
 *   ndp         must equal gf3d_ndp(itypsokern)
 *   dp          [nstations][ndp][3][nt]                (caller-allocated)
 *
 * The partials are analytic throughout -- no finite differences -- and are
 * per dyne-cm for slots 0-5, per degree for 6 and 7, per km for 8 and per
 * second for 9. So sum(moment[v]*dp[.][v][.][.]) over v = 0..5 reproduces
 * the seismogram, with the CMTSOLUTION's own numbers.
 *
 * A force source returns GF_ERR_ARG: there is nothing to differentiate
 * against a moment tensor.
 */
int gf3d_partials(gf3d_handle h, const gf3d_source *src, double t0_req,
                  int itypsokern, int nt, int ndp,
                  double *seis, double *dp, double *t, double *onset,
                  gf3d_location *loc);

#ifdef __cplusplus
}
#endif

#endif /* GF3D_H */
