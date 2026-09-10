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
 * test_gf_capi -- the C round trip, and that no path can stop
 *
 * Tier 2: needs the HDF5 build of the library and one of the example
 * databases; skips cleanly without them (see 9d.test_gf_capi.sh).
 *
 * This is the first test in the tree written in C, because it is the first
 * thing in the tree with a C ABI. It exists to check two different claims.
 *
 * The first is that the ABI is right: that a caller who has only
 * include/gf3d.h and lib/libgf3d.so can open a database, ask what is in it,
 * locate a source, and get seismograms and partials whose shape, order and
 * units are what the header says. The oracle for the numbers is an identity
 * rather than a reference file: summing the moment-tensor partials against
 * the CMTSOLUTION's own components must reproduce the seismogram, exactly
 * as it does in the Fortran tests, and that cannot come out right by
 * accident if an index is transposed.
 *
 * The second is that nothing here can end the process. Every deliberate
 * mistake below -- a missing database, a bad handle, a closed handle, a
 * wrong array length, a NaN, partials of a force source -- has to come back
 * as an error code with the program still running. That is the property
 * that lets this library be loaded into a Python interpreter, and the only
 * way to test it is to try each one and still reach the end.
 *
 * Usage: test_gf_capi <GFDB directory> <CMTSOLUTION>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gf3d.h"

static int nfail = 0;

static void ok(const char *name, int cond)
{
  if (cond) {
    printf("  ok   %s\n", name);
  } else {
    printf("  FAIL %s\n", name);
    nfail++;
  }
}

static void ok_err(const char *name, double err, double tol)
{
  if (err <= tol && !isnan(err)) {
    printf("  ok   %-38s error = %12.5e  (tol %12.5e)\n", name, err, tol);
  } else {
    printf("  FAIL %-38s error = %12.5e  (tol %12.5e)\n", name, err, tol);
    nfail++;
  }
}

/* the status code a call must return, reported with its message */
static void ok_status(const char *name, int got, int want)
{
  char msg[512], nm[64];
  gf3d_error_string(got, nm, (int)sizeof(nm));
  gf3d_last_error(msg, (int)sizeof(msg));
  if (got == want) {
    printf("  ok   %-38s -> %s\n", name, nm);
  } else {
    printf("  FAIL %-38s -> %s (wanted %d), %s\n", name, nm, want, msg);
    nfail++;
  }
}

/*
 * The CMTSOLUTION, parsed here in C.
 *
 * Deliberate: the C API takes numbers, never a path, precisely so that no
 * caller reaches get_cmt()'s stop statements. Parsing is the caller's job,
 * and this is what that costs -- about twenty lines. Every line of the
 * format is "name: value", so the value is what follows the last colon.
 */
static int read_cmtsolution(const char *path, gf3d_source *src)
{
  FILE *f = fopen(path, "r");
  char line[512];
  int iline = 0;

  if (f == NULL) return 1;

  memset(src, 0, sizeof(*src));
  src->source_type = GF_SRC_CMT;

  while (fgets(line, (int)sizeof(line), f) != NULL) {
    char *colon = strrchr(line, ':');
    double v = 0.0;
    if (iline > 0) {
      if (colon == NULL) { fclose(f); return 2; }
      v = atof(colon + 1);
    }
    switch (iline) {
      case 2:  src->time_shift = v; break;
      case 3:  src->hdur       = v; break;
      case 4:  src->latitude   = v; break;
      case 5:  src->longitude  = v; break;
      case 6:  src->depth_km   = v; break;
      case 7:  src->moment[0]  = v; break;
      case 8:  src->moment[1]  = v; break;
      case 9:  src->moment[2]  = v; break;
      case 10: src->moment[3]  = v; break;
      case 11: src->moment[4]  = v; break;
      case 12: src->moment[5]  = v; break;
      default: break;
    }
    iline++;
  }
  fclose(f);

  return (iline >= 13) ? 0 : 3;
}

int main(int argc, char **argv)
{
  const char *dbpath, *cmtpath;
  gf3d_handle h = 0, h2 = 0;
  gf3d_info info;
  gf3d_station sta;
  gf3d_location loc;
  gf3d_plan plan;
  gf3d_source src, force;
  char buf[128], name[16], unit[16];
  double *seis = NULL, *dp = NULL, *t = NULL, *onset = NULL;
  int ierr, ndp, i, j, k, ip, nsta, nt;
  int sz_source, sz_info, sz_station, sz_location, sz_plan;

  if (argc < 3) {
    fprintf(stderr, "usage: %s <GFDB> <CMTSOLUTION>\n", argv[0]);
    return 1;
  }
  dbpath = argv[1];
  cmtpath = argv[2];

  printf("\n ******************************\n");
  printf(" test_gf_capi\n");
  printf(" ******************************\n\n");

  /* ---------------------------------------------------------------- */
  printf(" 1. version, sizes and error strings\n");

  ierr = gf3d_version(buf, (int)sizeof(buf));
  ok("gf3d_version returns GF_OK", ierr == GF_OK);
  ok("version string is not empty", strlen(buf) > 0);
  printf("       library version: %s\n", buf);

  /*
   * The struct sizes: this header is transcribed by hand into Fortran
   * bind(C) types and again into Python ctypes.Structures, so the one thing
   * that can silently go wrong is a layout drift. Ask the library.
   */
  ierr = gf3d_sizeof(&sz_source, &sz_info, &sz_station, &sz_location, &sz_plan);
  ok("gf3d_sizeof returns GF_OK", ierr == GF_OK);
  ok("sizeof(gf3d_source) agrees", sz_source == (int)sizeof(gf3d_source));
  ok("sizeof(gf3d_info) agrees", sz_info == (int)sizeof(gf3d_info));
  ok("sizeof(gf3d_station) agrees", sz_station == (int)sizeof(gf3d_station));
  ok("sizeof(gf3d_location) agrees", sz_location == (int)sizeof(gf3d_location));
  ok("sizeof(gf3d_plan) agrees", sz_plan == (int)sizeof(gf3d_plan));

  /* NULL is allowed for every one of them */
  ierr = gf3d_sizeof(NULL, NULL, NULL, NULL, NULL);
  ok("gf3d_sizeof accepts NULL throughout", ierr == GF_OK);

  gf3d_error_string(GF_ERR_NO_ELEMENT, buf, (int)sizeof(buf));
  ok("error string for GF_ERR_NO_ELEMENT", strlen(buf) > 0);

  /* ---------------------------------------------------------------- */
  printf("\n 2. a database that is not there\n");

  ierr = gf3d_open("/nonexistent/gf3d/database", 0, &h);
  ok_status("opening a missing database", ierr, GF_ERR_NO_PATH);
  ok("no handle was handed out", h == 0);
  gf3d_last_error(buf, (int)sizeof(buf));
  ok("a message was left behind", strlen(buf) > 0);

  /* ---------------------------------------------------------------- */
  printf("\n 3. opening the example database\n");

  ierr = gf3d_open(dbpath, 0, &h);
  ok_status("gf3d_open", ierr, GF_OK);
  if (ierr != GF_OK) return 1;
  ok("handle is in range", h >= 1 && h <= GF3D_MAX_HANDLES);

  ierr = gf3d_get_info(h, &info);
  ok_status("gf3d_get_info", ierr, GF_OK);
  printf("       %d elements, %d stations, %d samples at dt = %g x %d\n",
         info.nelem, info.nstations, info.nt_subsampled, info.dt,
         info.subsample_step);
  ok("at least one element", info.nelem > 0);
  ok("at least one station", info.nstations > 0);
  ok("dt is positive", info.dt > 0.0);
  ok("subsample_step is positive", info.subsample_step > 0);
  ok("r_planet is a planet", info.r_planet > 1.0e6);
  ok("rhoav is a density", info.rhoav > 100.0);

  nsta = info.nstations;

  /* station names must be present and distinct */
  {
    int distinct = 1, nonempty = 1;
    char ids[64][GF3D_STRLEN];
    for (i = 0; i < nsta && i < 64; i++) {
      ierr = gf3d_get_station(h, i, &sta);
      if (ierr != GF_OK) { nonempty = 0; break; }
      strncpy(ids[i], sta.id, GF3D_STRLEN - 1);
      ids[i][GF3D_STRLEN - 1] = '\0';
      if (strlen(ids[i]) == 0) nonempty = 0;
      for (j = 0; j < i; j++) if (strcmp(ids[i], ids[j]) == 0) distinct = 0;
    }
    ok("every station has a name", nonempty);
    ok("station names are distinct", distinct);
    printf("       first station: %s at %.4f, %.4f\n", ids[0],
           sta.latitude, sta.longitude);
  }

  ok_status("station index -1 refused", gf3d_get_station(h, -1, &sta), GF_ERR_ARG);
  ok_status("station index nsta refused", gf3d_get_station(h, nsta, &sta), GF_ERR_ARG);

  /* ---------------------------------------------------------------- */
  printf("\n 4. the partial derivative table\n");

  ierr = gf3d_ndp(0, &ndp);
  ok("gf3d_ndp(0) = 0", ierr == GF_OK && ndp == 0);
  ierr = gf3d_ndp(1, &ndp);
  ok("gf3d_ndp(1) = 6", ierr == GF_OK && ndp == GF_NDP_MT);
  ierr = gf3d_ndp(2, &ndp);
  ok("gf3d_ndp(2) = 10", ierr == GF_OK && ndp == GF_NDP_LOC);
  ok_status("gf3d_ndp(3) refused", gf3d_ndp(3, &ndp), GF_ERR_ARG);

  ierr = gf3d_partial_name(0, name, (int)sizeof(name), unit, (int)sizeof(unit));
  ok("partial 0 is Mrr per dyne-cm",
     ierr == GF_OK && strcmp(name, "Mrr") == 0 && strncmp(unit, "m/dyne-cm", 9) == 0);
  ierr = gf3d_partial_name(9, name, (int)sizeof(name), unit, (int)sizeof(unit));
  ok("partial 9 is tim per second",
     ierr == GF_OK && strcmp(name, "tim") == 0 && strncmp(unit, "m/s", 3) == 0);
  ierr = gf3d_partial_name(6, name, (int)sizeof(name), NULL, 0);
  ok("partial 6 is lat, unit buffer NULL", ierr == GF_OK && strcmp(name, "lat") == 0);
  ok_status("partial index 10 refused",
            gf3d_partial_name(10, name, (int)sizeof(name), unit, (int)sizeof(unit)),
            GF_ERR_ARG);

  /* ---------------------------------------------------------------- */
  printf("\n 5. locating the validation source\n");

  ierr = read_cmtsolution(cmtpath, &src);
  ok("CMTSOLUTION parsed", ierr == 0);
  if (ierr != 0) return 1;
  printf("       %.4f, %.4f at %.2f km, hdur %.1f s, shift %.1f s\n",
         src.latitude, src.longitude, src.depth_km, src.hdur, src.time_shift);

  ierr = gf3d_locate(h, src.latitude, src.longitude, src.depth_km, &loc);
  ok_status("gf3d_locate", ierr, GF_OK);
  if (ierr != GF_OK) return 1;
  printf("       element %d (%s) at xi,eta,gamma = %.6f %.6f %.6f\n",
         loc.ielem, loc.morton_hex, loc.xi, loc.eta, loc.gamma);
  ok("element index is in range", loc.ielem >= 1 && loc.ielem <= info.nelem);
  ok("Morton code came back", strlen(loc.morton_hex) > 0);
  ok("the mapped point is where it was asked for", loc.distance_km < 1.0e-6);
  ok("xi is inside the accepted range", fabs(loc.xi) <= 1.1);
  ok("the 27-anchor residual is small", loc.anchor_err < 1.0e-5);

  /* the NaN that would otherwise reach the kd-tree's own stop */
  ok_status("NaN latitude refused",
            gf3d_locate(h, NAN, src.longitude, src.depth_km, &loc), GF_ERR_ARG);
  ok_status("a point on the far side of the planet refused",
            gf3d_locate(h, -src.latitude, src.longitude + 180.0, 10.0, &loc),
            GF_ERR_NO_ELEMENT);

  /* ---------------------------------------------------------------- */
  printf("\n 6. the plan\n");

  ierr = gf3d_get_plan(h, &src, -1.0, &plan);
  ok_status("gf3d_get_plan with specfem's own t0", ierr, GF_OK);
  if (ierr != GF_OK) return 1;
  nt = plan.nt;
  printf("       nt = %d (%d stored + %d prepended), dt_sub = %g, b = %g\n",
         plan.nt, plan.nt_db, plan.npad, plan.dt_sub, plan.t_first);
  ok("nt is positive", nt > 0);
  ok("nt is the stored length plus the padding", plan.nt == plan.nt_db + plan.npad);
  ok("dt_sub is dt times the subsampling",
     fabs(plan.dt_sub - plan.dt * plan.subsample_step) < 1.0e-12 * plan.dt_sub);
  ok("the axis starts at or before -1.5*hdur", plan.t_first <= -1.5 * src.hdur);
  ok("a Heaviside conversion for a moment tensor", plan.kind_stf == 2);

  /* ---------------------------------------------------------------- */
  printf("\n 7. seismograms\n");

  seis = (double *)malloc((size_t)nsta * GF_NCOMP * nt * sizeof(double));
  t = (double *)malloc((size_t)nt * sizeof(double));
  onset = (double *)malloc((size_t)nsta * sizeof(double));
  if (seis == NULL || t == NULL || onset == NULL) {
    fprintf(stderr, "out of memory\n");
    return 1;
  }

  ierr = gf3d_seismograms(h, &src, -1.0, nt, seis, t, onset, &loc);
  ok_status("gf3d_seismograms", ierr, GF_OK);
  if (ierr != GF_OK) return 1;

  {
    int allfinite = 1, nonzero = 0;
    double peak = 0.0, dtmax = 0.0;
    for (i = 0; i < nsta * GF_NCOMP * nt; i++) {
      if (!isfinite(seis[i])) allfinite = 0;
      if (seis[i] != 0.0) nonzero = 1;
      if (fabs(seis[i]) > peak) peak = fabs(seis[i]);
    }
    ok("every sample is finite", allfinite);
    ok("the traces are not all zero", nonzero);
    printf("       peak displacement %.4e m, onset ratio %.2e\n", peak, onset[0]);

    for (i = 1; i < nt; i++) {
      double d = fabs((t[i] - t[i-1]) - plan.dt_sub);
      if (d > dtmax) dtmax = d;
    }
    ok_err("the time axis is uniform at dt_sub", dtmax, 1.0e-12 * plan.dt_sub);
    ok_err("t[0] is the planned first sample", fabs(t[0] - plan.t_first),
           1.0e-12 * fabs(plan.t_first));
  }

  ok_status("the wrong nt refused",
            gf3d_seismograms(h, &src, -1.0, nt + 1, seis, t, onset, NULL), GF_ERR_ARG);

  /* ---------------------------------------------------------------- */
  printf("\n 8. partial derivatives, and the identity they satisfy\n");

  ierr = gf3d_ndp(2, &ndp);
  dp = (double *)malloc((size_t)nsta * ndp * GF_NCOMP * nt * sizeof(double));
  if (dp == NULL) { fprintf(stderr, "out of memory\n"); return 1; }

  ierr = gf3d_partials(h, &src, -1.0, 2, nt, ndp, seis, dp, t, onset, NULL);
  ok_status("gf3d_partials, itypsokern = 2", ierr, GF_OK);
  if (ierr != GF_OK) return 1;

  {
    int allfinite = 1;
    double worst = 0.0, peak = 0.0;

    for (i = 0; i < nsta * ndp * GF_NCOMP * nt; i++)
      if (!isfinite(dp[i])) allfinite = 0;
    ok("every partial is finite", allfinite);

    /*
     * The identity of Stage 6: the six moment-tensor partials are per
     * dyne-cm and linear, so contracting them with the CMTSOLUTION's own
     * components must give the seismogram back. It needs no reference file,
     * and no transposition of the four-dimensional index can survive it.
     */
    for (i = 0; i < nsta; i++) {
      for (j = 0; j < GF_NCOMP; j++) {
        for (k = 0; k < nt; k++) {
          double s = seis[(i * GF_NCOMP + j) * nt + k];
          double sum = 0.0;
          for (ip = 0; ip < GF_NDP_MT; ip++)
            sum += src.moment[ip] * dp[((i * ndp + ip) * GF_NCOMP + j) * nt + k];
          if (fabs(s) > peak) peak = fabs(s);
          if (fabs(sum - s) > worst) worst = fabs(sum - s);
        }
      }
    }
    ok_err("sum(M_v dp_v) reproduces the seismogram", worst / peak, 1.0e-12);
  }

  ok_status("the wrong ndp refused",
            gf3d_partials(h, &src, -1.0, 2, nt, GF_NDP_MT, seis, dp, t, onset, NULL),
            GF_ERR_ARG);
  ok_status("itypsokern 0 refused (use gf3d_seismograms)",
            gf3d_partials(h, &src, -1.0, 0, nt, 0, seis, dp, t, onset, NULL),
            GF_ERR_ARG);

  /* partials of a force source have nothing to differentiate against */
  memset(&force, 0, sizeof(force));
  force.source_type = GF_SRC_FORCE;
  force.latitude = src.latitude;
  force.longitude = src.longitude;
  force.depth_km = src.depth_km;
  force.hdur = 45.0;
  force.force_stf = 0;
  force.force_factor = 1.0e15;
  force.force_dir[0] = 1.0;
  ok_status("partials of a force source refused",
            gf3d_partials(h, &force, -1.0, 1, nt, GF_NDP_MT, seis, dp, t, onset, NULL),
            GF_ERR_ARG);

  /* but its seismograms are fine */
  {
    gf3d_plan fplan;
    ierr = gf3d_get_plan(h, &force, -1.0, &fplan);
    ok_status("a force source plans", ierr, GF_OK);
    if (ierr == GF_OK && fplan.nt <= nt) {
      ierr = gf3d_seismograms(h, &force, -1.0, fplan.nt, seis, t, onset, NULL);
      ok_status("a force source extracts", ierr, GF_OK);
    }
  }

  /* a source type that is neither */
  {
    gf3d_source bad = src;
    bad.source_type = 7;
    ok_status("an unknown source type refused",
              gf3d_get_plan(h, &bad, -1.0, &plan), GF_ERR_ARG);
  }

  /* ---------------------------------------------------------------- */
  printf("\n 9. handles\n");

  ok_status("handle 0 refused", gf3d_get_info(0, &info), GF_ERR_ARG);
  ok_status("handle -1 refused", gf3d_get_info(-1, &info), GF_ERR_ARG);
  ok_status("handle 99 refused", gf3d_get_info(99, &info), GF_ERR_ARG);

  /*
   * Two handles on the same database at once, then closing one and using
   * the other. Closing releases the process-wide search tree even though
   * the second handle is still open, so the extraction below is what
   * proves the tree is rebuilt on demand rather than silently missing.
   */
  ierr = gf3d_open(dbpath, 0, &h2);
  ok_status("a second handle on the same database", ierr, GF_OK);
  ok("the two handles differ", h2 != h);

  if (ierr == GF_OK) {
    ierr = gf3d_seismograms(h2, &src, -1.0, nt, seis, t, onset, NULL);
    ok_status("the second handle extracts", ierr, GF_OK);

    ierr = gf3d_close(h);
    ok_status("closing the first handle", ierr, GF_OK);

    ierr = gf3d_seismograms(h2, &src, -1.0, nt, seis, t, onset, NULL);
    ok_status("the second handle still extracts", ierr, GF_OK);

    ok_status("the closed handle is refused", gf3d_get_info(h, &info), GF_ERR_ARG);

    ierr = gf3d_close(h2);
    ok_status("closing the second handle", ierr, GF_OK);
    ok_status("closing it twice is refused", gf3d_close(h2), GF_ERR_ARG);
  }

  free(seis);
  free(dp);
  free(t);
  free(onset);

  printf("\n");
  if (nfail == 0) {
    printf(" test_gf_capi: all assertions passed\n\n");
    return 0;
  }
  printf(" test_gf_capi: %d assertion(s) FAILED\n\n", nfail);
  return 1;
}
