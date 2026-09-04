/*
!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
*/

/*
 * Directory enumeration for the Green function database reader.
 *
 * Why this is in C: the station list is not recorded anywhere in the
 * database. {GFDB}/mesh_info.h5 holds simulation-wide parameters and
 * {GFDB}/manifest.csv holds the element index, but the only record of
 * which stations exist is the set of files in {GFDB}/stations/. Fortran
 * has no portable way to read a directory, and shelling out to `ls`
 * would need a writable temporary file and a shell.
 *
 * Follows the FC_FUNC_ convention used by the rest of the tree
 * (src/shared/param_reader.c, src/shared/binary_c_io.c); there is no
 * bind(C) anywhere in specfem3d_globe.
 *
 * The interface is a single stateless call. Call it once with nmax = 0
 * to obtain the number of matching entries, then again with a buffer.
 * Names come back sorted with strcmp() so that station and element
 * ordering is reproducible across runs and filesystems.
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

/* upper bound on a single directory entry name we are willing to keep */
#define GF_DIRLIST_NAME_MAX 256

static int gf_dirlist_strcmp(const void *a, const void *b) {
  return strcmp(*(const char **)a, *(const char **)b);
}

/*
 * gf_dir_scan(path, path_len, suffix, suffix_len, want_dirs,
 *             nmax, namelen, names, count, ierr)
 *
 *   path       directory to scan (not NUL-terminated; path_len chars)
 *   suffix     if suffix_len > 0, keep only entries ending in this string
 *   want_dirs  1 -> keep only subdirectories, 0 -> keep only regular files
 *   nmax       number of slots in `names`; pass 0 to count only
 *   namelen    width of one slot in `names`, in characters
 *   names      output buffer, nmax * namelen chars, blank padded
 *              (Fortran character(len=namelen), dimension(nmax))
 *   count      number of matching entries found (may exceed nmax)
 *   ierr       0 on success
 *              1 directory could not be opened
 *              2 an entry name does not fit in namelen
 *              3 out of memory
 *
 * Entries "." and ".." are always skipped.
 */
void
FC_FUNC_(gf_dir_scan,GF_DIR_SCAN)(char *path, int *path_len,
                                  char *suffix, int *suffix_len,
                                  int *want_dirs,
                                  int *nmax, int *namelen, char *names,
                                  int *count, int *ierr) {

  DIR *dir;
  struct dirent *ent;
  char dirname[4096];
  char fullpath[4096 + GF_DIRLIST_NAME_MAX + 2];
  struct stat st;
  char **found = NULL;
  int nfound = 0, ncap = 0;
  int i, j, n, len, slen;

  *count = 0;
  *ierr = 0;

  n = *path_len;
  if (n < 0) n = 0;
  if (n >= (int) sizeof(dirname)) n = (int) sizeof(dirname) - 1;
  memcpy(dirname, path, (size_t) n);
  dirname[n] = '\0';

  slen = *suffix_len;
  if (slen < 0) slen = 0;

  dir = opendir(dirname);
  if (dir == NULL) { *ierr = 1; return; }

  while ((ent = readdir(dir)) != NULL) {

    if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) continue;

    len = (int) strlen(ent->d_name);
    if (len >= GF_DIRLIST_NAME_MAX) continue;

    /* suffix filter */
    if (slen > 0) {
      if (len < slen) continue;
      if (strncmp(ent->d_name + (len - slen), suffix, (size_t) slen) != 0) continue;
    }

    /* type filter.
       note: d_type is not portable (and is DT_UNKNOWN on some filesystems),
       so stat() the entry rather than trusting it */
    snprintf(fullpath, sizeof(fullpath), "%s/%s", dirname, ent->d_name);
    if (stat(fullpath, &st) != 0) continue;

    if (*want_dirs) {
      if (!S_ISDIR(st.st_mode)) continue;
    } else {
      if (!S_ISREG(st.st_mode)) continue;
    }

    if (nfound == ncap) {
      char **tmp;
      ncap = (ncap == 0) ? 64 : 2 * ncap;
      tmp = (char **) realloc(found, (size_t) ncap * sizeof(char *));
      if (tmp == NULL) { *ierr = 3; goto cleanup; }
      found = tmp;
    }
    found[nfound] = strdup(ent->d_name);
    if (found[nfound] == NULL) { *ierr = 3; goto cleanup; }
    nfound++;
  }

  /* readdir() order is filesystem-dependent; sort so that station and
     element indices are stable from run to run */
  if (nfound > 1) qsort(found, (size_t) nfound, sizeof(char *), gf_dirlist_strcmp);

  *count = nfound;

  for (i = 0; i < nfound && i < *nmax; i++) {
    len = (int) strlen(found[i]);
    if (len > *namelen) { *ierr = 2; break; }
    memcpy(names + (size_t) i * (size_t) (*namelen), found[i], (size_t) len);
    for (j = len; j < *namelen; j++) names[(size_t) i * (size_t) (*namelen) + j] = ' ';
  }

cleanup:
  for (i = 0; i < nfound; i++) free(found[i]);
  free(found);
  closedir(dir);
}

/*
 * gf_path_exists(path, path_len, want_dirs, exists)
 *
 *   exists = 1 if `path` exists and is of the requested kind
 *            (want_dirs = 1 -> directory, 0 -> regular file), else 0
 *
 * Fortran's inquire(file=...,exist=) is not required by the standard to
 * report on directories, and compilers disagree, so directory existence
 * is checked here.
 */
void
FC_FUNC_(gf_path_exists,GF_PATH_EXISTS)(char *path, int *path_len,
                                        int *want_dirs, int *exists) {

  char name[4096];
  struct stat st;
  int n;

  *exists = 0;

  n = *path_len;
  if (n < 0) n = 0;
  if (n >= (int) sizeof(name)) n = (int) sizeof(name) - 1;
  memcpy(name, path, (size_t) n);
  name[n] = '\0';

  if (stat(name, &st) != 0) return;

  if (*want_dirs) {
    if (S_ISDIR(st.st_mode)) *exists = 1;
  } else {
    if (S_ISREG(st.st_mode)) *exists = 1;
  }
}
