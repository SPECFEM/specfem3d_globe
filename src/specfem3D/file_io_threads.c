/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

/* ---------------------------------------

// asynchronous file i/o

 --------------------------------------- */

// --------------------------------------------------------------------------------
// ------------------------- Non-blocking IO routines -----------------------------
// --------------------------------------------------------------------------------

// In order to overlap disk I/O with computation, we have defined
// several routines which can do disk I/O in a non-blocking way. The
// first launches a thread (using pthreads), which writes to
// disk using fwrite. The second version, attempts to use the linux
// system call aio_write, which produces a similar result, but is not
// able to use buffered I/O resulting in slower performance. Both
// these routines must copy the memory region they are saving, because
// it will change as the output is uptdated during the next time
// step. This could be solved by rotating pointers, but this is not
// very easy in Fortran (especially without changing the fortran code)


#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

// debug: outputs traces
#define DEBUG 0
#if DEBUG == 1
#define TRACE(x) printf("%s\n",x);
#else
#define TRACE(x) // printf("%s\n",x);
#endif

#define exit_error(msg) { fputs (msg,stderr); abort(); }

// --------------------------------------------------------------------------------

// Pthreads version -- faster than aio version

// --------------------------------------------------------------------------------

// threading version
#include <pthread.h>

struct threadInfo{
  int  pid;
  char* buffer;
  size_t bytes_to_rw;
  int it_sub;
  bool finished;
  bool started;
};

// --------------------------------------------------------------------------------

// adjoint sources threads

// --------------------------------------------------------------------------------

struct threadInfo ptDataAdj;

pthread_t adj_io_thread;

// fortran function
extern void FC_FUNC_(read_adjoint_sources_local,READ_ADJOINT_SOURCES_LOCAL)(char*, const int* , const int*);

// Thread that does actual read in adjoint sources. dummy argument is needed to avoid compiler warning.
void *fread_adj_thread(void* dummy) {

  // note: having it_sub_adj as function argument and using int it_sub = (* (int*) it_sub_adj); did not work...
  TRACE("fread_adj_thread");
  // debug info
  //printf("adjoint thread: nadj_rec_local = %i - it_sub_adj = %i \n",ptDataAdj.pid,ptDataAdj.it_sub);

  // calls fortran function
  // see file: src/specfem3D/read_adjoint_sources.f90
  FC_FUNC_(read_adjoint_sources_local,READ_ADJOINT_SOURCES_LOCAL)(ptDataAdj.buffer,&ptDataAdj.pid,&ptDataAdj.it_sub);

  // reading done
  ptDataAdj.finished = true;

  // exit thread
  pthread_exit(NULL);
  return NULL;  // Never used, but remove warnings.
}

// Waits until thread is finished with I/O
void wait_adj_io_thread() {

  TRACE("wait_adj_io_thread");

  int rc;
  // checks if thread still runs
  if (ptDataAdj.started ) {
    void* status;
    rc = pthread_join(adj_io_thread, &status);
    if (rc ) {
      printf("Error; return code from pthread_join() is %d\n", rc);
      exit_error("Error in wait_adj_io_thread: thread_join failed");
    }
    // checks finished flag
    assert(ptDataAdj.finished == true); // Adjoint thread has completed, but somehow it isn't finished?

    // reset
    ptDataAdj.started = false;
  }
}

// initializes adjoint thread
void
FC_FUNC_(prepare_adj_io_thread,CREATE_IO_ADJ_THREAD)(char *buffer, long* length, int* nadj_rec_local) {

  TRACE("prepare_adj_io_thread");

  size_t bytes_to_read = *length;

  // checks if buffer valid
  assert(buffer != NULL); // "Adjoint thread: associated buffer is invalid"
  if (bytes_to_read <= 0 ) exit_error("Adjoint thread: associated buffer length is invalid");

  // initializes thread info
  ptDataAdj.started = false;
  ptDataAdj.finished = false;
  ptDataAdj.pid = *nadj_rec_local;
  ptDataAdj.buffer = buffer;
  ptDataAdj.bytes_to_rw = bytes_to_read;
  ptDataAdj.it_sub = 0;
}

// creates thread for reading adjoint sources
void
FC_FUNC_(read_adj_io_thread,CREATE_IO_ADJ_THREAD)(int* it_sub_adj) {

  TRACE("read_adj_io_thread");
  // debug
  //printf("read_adj_io_thread: it_sub_adj = %i \n",*it_sub_adj);

  int rc;

  // checks if buffer valid
  assert(ptDataAdj.buffer != NULL); // Adjoint thread: associated buffer is invalid

  // waits until previous thread finishes
  wait_adj_io_thread();

  // prepares the thread
  pthread_attr_t attr;

  rc = pthread_attr_init(&attr);
  if (rc != 0 ) exit_error("Adjoint thread: initialization of thread failed");

  rc = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  if (rc != 0 ) exit_error("Adjoint thread: setting thread state failed");

  // sets new thread info
  ptDataAdj.started = true;
  ptDataAdj.finished = false;
  ptDataAdj.it_sub = *it_sub_adj;

  // create and launch the thread.
  // note: using it_sub_adj as argument (void*) it_sub_adj did not work...
  rc = pthread_create(&adj_io_thread,&attr,fread_adj_thread,NULL);
  if (rc != 0 ) exit_error("Adjoint thread: creating thread failed");

}

// synchronizes with read thread to make sure buffer is read in
void
FC_FUNC_(sync_adj_io_thread,SYNC_ADJ_IO_THREAD)(char *adj_sourcearrays) {

  TRACE("sync_adj_io_thread");

  // Wait until previous thread finishes reading
  wait_adj_io_thread();

  // the buffer is updated
  // copies the thread buffer to the output
  memcpy(adj_sourcearrays,ptDataAdj.buffer,ptDataAdj.bytes_to_rw);
}



// --------------------------------------------------------------------------------

// absorbing boundaries threads

// --------------------------------------------------------------------------------
// not used yet...
/*
struct threadInfo ptData[ABS_FILEID];

pthread_t io_thread[ABS_FILEID];

void *fwrite_thread(void *fileID);
void *fread_thread(void *fileID);
void prepare_thread_io(int *fid);
void write_abs_ptio(int *fid, char *buffer, int *length, int *index);
void wait_io_thread(int *fid);

// Thread that does actual fwrite.
void *fwrite_thread(void *fileID)
{
  int fid;
  fid = (int)fileID;

  fwrite(ptData[fid].buffer, 1, ptData[fid].bytes_to_rw,fp_abs[fid]);

  ptData[fid].finished = true;
  pthread_exit(NULL);
}

// Thread that does actual fread.
void *fread_thread(void *fileID)
{
  int fid;
  fid = (int)fileID;

  fread(ptData[fid].buffer, 1, ptData[fid].bytes_to_rw,fp_abs[fid]);

  ptData[fid].finished = true;
  pthread_exit(NULL);
}


// Setup thread values including the memory region
void prepare_thread_io(int *fid) {
  ptData[*fid].started = false;
  ptData[*fid].finished = false;
  ptData[*fid].pid = *fid;
  ptData[*fid].buffer = NULL;
}

// Starts a thread to write the desired buffer to disk
void write_abs_ptio(int *fid, char *buffer, int *length, int *index) {

  // allocates buffer if not done so yet
  size_t bytes_to_write = *length;
  if (ptData[*fid].buffer == NULL) {
    ptData[*fid].buffer = (char*)malloc(bytes_to_write);
  }

  // Wait until previous thread finishes writing
  wait_io_thread(fid);

  // prepare the thread for writing.
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  ptData[*fid].started = true;
  ptData[*fid].finished = false;
  ptData[*fid].bytes_to_rw = bytes_to_write;

  // the input buffer will be updated while we write the data to disk
  // -- Thus we have to make a copy so that the data doesn't change
  // while we are writing it to disk.

  memcpy(ptData[*fid].buffer,buffer,ptData[*fid].bytes_to_rw);

  // create and launch the thread.
  int rc = pthread_create(&io_thread[*fid],&attr,fwrite_thread,(void*) *fid);

}

// Starts a thread to read the desired buffer from disk
void read_abs_ptio(int *fid, char *buffer, int *length, int *index) {

  // allocates buffer if not done so yet
  size_t bytes_to_read = *length;
  if (ptData[*fid].buffer == NULL) {
    ptData[*fid].buffer = (char*)malloc(bytes_to_read);
  }

  // Wait until previous thread finishes reading
  wait_io_thread(fid);

  // prepare the thread for reading.
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  ptData[*fid].started = true;
  ptData[*fid].finished = false;
  ptData[*fid].bytes_to_rw = bytes_to_read;

  // create and launch the thread.
  int rc = pthread_create(&io_thread[*fid],&attr,fread_thread,(void*) *fid);

}


// Waits until thread is finished with I/O
void wait_io_thread(int *fid) {
  int rc;
  // checks if thread still runs
  if (ptData[*fid].started ) {
    void* status;
    rc = pthread_join(io_thread[*fid], &status);
    if (rc) {
      printf("Error; return code from pthread_join() is %d\n", rc);
      exit(-1);
    }
    // checks finished flag
    assert(ptData[*fid].finished == true && "Thread has completed, but somehow it isn't finished?");

    // reset
    ptData[*fid].started = false;
  }
}

// synchronizes with read thread to make sure buffer is read in
void sync_buffer_io_thread(int *fid, char *buffer, int *length) {

  // checks buffer length
  size_t bytes_to_read = *length;
  if (ptData[*fid].buffer == NULL) {
    assert("Associated thread has no buffer");
  }
  assert(ptData[*fid].bytes_to_rw != bytes_to_read && "Associated thread has different buffer length");

  // Wait until previous thread finishes reading
  wait_io_thread(fid);

  // the input buffer will be updated while we read the data to disk
  // -- Thus we have to make a copy so that the data doesn't change
  // while we are reading into it from disk.
  //
  // copies the thread buffer to the output
  memcpy(buffer,ptData[*fid].buffer,ptData[*fid].bytes_to_rw);
}
*/



// --------------------------------------------------------------------------------

// async IO

// --------------------------------------------------------------------------------
// unused yet...
/*
#include <aio.h>
#include <signal.h>
#include <errno.h>

// Async_IO
static int fp_aio[ABS_FILEID];
struct aiocb fIOw[ABS_FILEID];

// --------------------------------------------------------------------------------

// write routines using AIO

void write_abs_aio(int *fid, char *buffer, int *length, int *index) {

  printf("aio_error return code=%d\n",aio_error(&fIOw[*fid]));

  // wait for previous IO to finish -- aio_error returns 0 when the
  // job specified by fIOw is finished successfully
  while(aio_error(&fIOw[*fid])) {}

  // setup async IO for this transfer
  fIOw[*fid].aio_fildes = fp_aio[*fid];
  // fIOw[*fid].aio_offset = 0;
  fIOw[*fid].aio_buf = buffer;
  size_t write_length = *length;
  fIOw[*fid].aio_nbytes = write_length;
  fIOw[*fid].aio_reqprio = 0;
  fIOw[*fid].aio_sigevent.sigev_notify = SIGEV_NONE;
  printf("aio_writing %ld bytes\n",write_length);
  aio_write(&fIOw[*fid]);

}

void close_file_abs_aio(int * fid) {

  // wait for previous IO to finish -- aio_error returns 0 when the
  // job specified by fIOw is finished successfully
  while(aio_error(&fIOw[*fid])) {}

  // closes file

  close(fp_aio[*fid]);

}
*/
