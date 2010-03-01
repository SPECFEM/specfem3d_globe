/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            March 2010
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

// after Brian's function

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static int fd;

void
FC_FUNC_(open_file_create,OPEN_FILE)(char *file) {
  /*    fprintf(stderr, "Opening file: %s\n", file); */
  fd = open(file, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if(fd == -1) {
    fprintf(stderr, "Error opening file: %s exiting\n", file);
    exit(-1);
  }
}

void
FC_FUNC_(open_file_append,OPEN_FILE)(char *file) {
  /*    fprintf(stderr, "Opening file: %s\n", file); */
  fd = open(file, O_WRONLY | O_CREAT | O_APPEND, 0644);
  if(fd == -1) {
    fprintf(stderr, "Error opening file: %s exiting\n", file);
    exit(-1);
  }
}

void
FC_FUNC_(close_file,CLOSE_FILE)() {
  /*    fprintf(stderr, "Closing file\n"); */
  close(fd);
}

void
FC_FUNC_(write_integer,WRITE_INTEGER)(int *z) {
  write(fd, z, sizeof(int));
}

void
FC_FUNC_(write_real,WRITE_REAL)(float *z) {
  write(fd, z, sizeof(float));
}

/* BS BS begin. Added section for writing SAC binary data*/
void
FC_FUNC_(write_n_real,WRITE_N_REAL)(float *z,int *n) {
  write(fd, z, *n*sizeof(float));
}

void
FC_FUNC_(write_character,WRITE_CHARACTER)(char *z, int *lchar) {
  write(fd, z, *lchar*sizeof(char));
}

// LQY -- added for combine_vol/surf_data to write multiple binary files simultaneously --

void
FC_FUNC_(open_file_fd,OPEN_FILE_FD)(char *file, int *pfd) {
  /*    fprintf(stderr, "Opening file: %s\n", file); */
  *pfd = open(file, O_WRONLY | O_CREAT, 0644);
  if(*pfd == -1) {
    fprintf(stderr, "Error opening file: %s exiting\n", file);
    exit(-1);
  }
}

void
FC_FUNC_(close_file_fd,CLOSE_FILE_FD)(int *pfd) {
  /*    fprintf(stderr, "Closing file\n"); */
  close(*pfd);
}

void
FC_FUNC_(write_integer_fd,WRITE_INTEGER_FD)(int *pfd, int *z) {
  write(*pfd, z, sizeof(int));
}

void
FC_FUNC_(write_real_fd,WRITE_REAL_FD)(int *pfd, float *z) {
  write(*pfd, z, sizeof(float));
}

/* BS BS begin. Added section for writing SAC binary data*/
void
FC_FUNC_(write_n_real_fd,WRITE_N_REAL_FD)(int *pfd, float *z,int *n) {
  write(*pfd, z, *n*sizeof(float));
}

void
FC_FUNC_(write_character_fd,WRITE_CHARACTER_FD)(int *pfd, char *z, int *lchar) {
  write(*pfd, z, *lchar*sizeof(char));
}
