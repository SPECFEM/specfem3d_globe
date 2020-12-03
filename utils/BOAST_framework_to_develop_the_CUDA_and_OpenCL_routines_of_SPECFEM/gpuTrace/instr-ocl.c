#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
#include "ldChecker.h"
#include "ocl_helper.h"

#include <inttypes.h>
#include <signal.h>
#include <pthread.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define NS_TO_MS(_ns_) ((_ns_)/1000000)
const char* clewErrorString (cl_int error);

static int mocl_errcode;
static inline cl_int _clCheck(cl_int errcode, const char *file, int line, const char *func) {
  mocl_errcode = errcode;
  if (mocl_errcode != CL_SUCCESS) {
    error ("%d/%s at %s:%d %s\n", mocl_errcode,
           clewErrorString(mocl_errcode),
           file, line, func);
  }
  return errcode;
}
#define clCheck(to_check) _clCheck(to_check,__FILE__, __LINE__,  __func__)
#define clck_(var) var); clCheck(*var

struct ld_ocl_s {
  cl_context context;
  cl_command_queue command_queue;
} ldOclEnv;

struct ld_program_s {
  unsigned int uid;
  cl_program handle;
  char *source;
  unsigned int released;
};

struct ld_mem_offset_s {
  unsigned int uid;
  cl_mem handle;
  size_t size;
  ld_flags flags;
  size_t offset;
  unsigned int released;
  struct ld_mem_s *parent;
};

struct ld_queue_s {
  unsigned int uid;
  cl_command_queue handle;
};


cl_command_queue (*real_clCreateCommandQueue)(cl_context context,
                                              cl_device_id device,
                                              cl_command_queue_properties properties,
                                              cl_int *errcode_ret);
cl_int (*real_clFinish) (cl_command_queue command_queue);
cl_int (*real_clReleaseCommandQueue) (cl_command_queue command_queue);
cl_int (*real_clSetKernelArg) (cl_kernel kernel, cl_uint arg_index,
                               size_t arg_size, const void *arg_value);
cl_mem (*real_clCreateBuffer) (cl_context context, cl_mem_flags flags,
                               size_t size, void *host_ptr,
                               cl_int *errcode_ret);
cl_mem (*real_clCreateSubBuffer) (cl_mem buffer,
                             cl_mem_flags flags,
                             cl_buffer_create_type buffer_create_type,
                             const void *buffer_create_info,
                             cl_int *errcode_ret);
cl_int (*real_clReleaseMemObject) (cl_mem memobj);
cl_kernel (*real_clCreateKernel) (cl_program  program, const char *kernel_name,
                                  cl_int *errcode_ret);
cl_program (*real_clCreateProgramWithSource) (cl_context context,
                                              cl_uint count,
                                              const char **strings,
                                              const size_t *lengths,
                                              cl_int *errcode_ret);
cl_int (*real_clEnqueueNDRangeKernel) (cl_command_queue command_queue,
                                       cl_kernel kernel, cl_uint work_dim,
                                       const size_t *global_work_offset,
                                       const size_t *global_work_size,
                                       const size_t *local_work_size,
                                       cl_uint num_events_in_wait_list,
                                       const cl_event *event_wait_list,
                                       cl_event *event);
cl_int (*real_clEnqueueWriteBuffer) (cl_command_queue command_queue,
                                     cl_mem buffer,
                                     cl_bool blocking_write,
                                     size_t offset,
                                     size_t size,
                                     const void *ptr,
                                     cl_uint num_events_in_wait_list,
                                     const cl_event *event_wait_list,
                                     cl_event *event);
cl_int (*real_clEnqueueReadBuffer) (cl_command_queue command_queue,
                                    cl_mem buffer,
                                    cl_bool blocking_read,
                                    size_t offset,
                                    size_t size,
                                    void *ptr,
                                    cl_uint num_events_in_wait_list,
                                    const cl_event *event_wait_list,
                                    cl_event *event);
cl_int (*real_clReleaseKernel) (cl_kernel kernel);
cl_int (*real_clReleaseProgram) (cl_program program);

struct ld_bindings_s ocl_bindings[] = {
  {"clSetKernelArg", (void **) &real_clSetKernelArg},
  {"clCreateBuffer", (void **) &real_clCreateBuffer},
  {"clCreateSubBuffer", (void **) &real_clCreateSubBuffer},
  {"clCreateKernel", (void **) &real_clCreateKernel},
  {"clCreateProgramWithSource", (void **) &real_clCreateProgramWithSource},
  {"clEnqueueNDRangeKernel", (void **) &real_clEnqueueNDRangeKernel},
  {"clEnqueueWriteBuffer", (void **) &real_clEnqueueWriteBuffer},
  {"clEnqueueReadBuffer", (void **) &real_clEnqueueReadBuffer},
  {"clCreateCommandQueue", (void **) &real_clCreateCommandQueue},
  {"clFinish", (void **) &real_clFinish},
  {"clReleaseCommandQueue", (void **) &real_clReleaseCommandQueue},
  {"clReleaseMemObject", (void **) &real_clReleaseMemObject},
  {"clReleaseKernel", (void **) &real_clReleaseKernel},
  {"clReleaseProgram", (void **) &real_clReleaseProgram},
  {NULL, NULL}
};

CREATE_HASHMAP(program, cl_program, 100)
CREATE_HASHMAP(kernel, cl_kernel, 100)
CREATE_HASHMAP(mem, cl_mem, 400)
CREATE_HASHMAP(mem_offset, cl_mem, 200)
CREATE_HASHMAP(queue, cl_command_queue, 10)

/* ************************************************************************* */

static void sig_handler(int signum);

#if ENABLE_KERNEL_PROFILING == 1
pthread_mutex_t stats_lock;
static void print_statistics(int force);
#endif

int ocl_getBufferContent (struct ld_mem_s *ldBuffer, void *buffer,
                          size_t offset, size_t size);
int ocl_triggerKernelExecution (struct ld_kernel_s *ldKernel,
                                const struct work_size_s *work_sizes,
                                unsigned int work_dim);
void *ocl_setParameterValue (struct ld_kernel_s *ldKernel,
                             struct ld_kern_param_s *ldParam,
                             void *buffer,  size_t size);
int ocl_getAndReleaseParameterValue (struct ld_kernel_s *ldKernel,
                                     struct ld_kern_param_s *ldParam,
                                     void *buffer_handle,
                                     void *buffer, size_t size);

struct ld_mem_s *create_ocl_buffer (cl_mem handle);

void init_ocl_ldchecker(void) {
  struct callback_s callbacks = {ocl_getBufferContent, ocl_setParameterValue, ocl_getAndReleaseParameterValue, ocl_triggerKernelExecution};

#if ENABLE_KERNEL_PROFILING == 1
  pthread_mutex_init(&stats_lock, NULL);

  signal(SIGUSR1, sig_handler);
#endif

  init_ldchecker(callbacks, ocl_bindings);
  create_ocl_buffer(NULL);
  ocl_init_helper();
}

void UNUSED sig_handler(int signum) {
  switch(signum) {
#if ENABLE_KERNEL_PROFILING == 1
  case SIGUSR1:
    print_statistics(1);
    break;
#endif
  }
}

cl_command_queue clCreateCommandQueue(cl_context context,
                                      cl_device_id device,
                                      cl_command_queue_properties properties,
                                      cl_int *errcode_ret)
{
  init_ocl_ldchecker();

#if ENABLE_KERNEL_PROFILING == 1
  properties |=  CL_QUEUE_PROFILING_ENABLE;
#endif

  ldOclEnv.command_queue = real_clCreateCommandQueue(context, device,
                                                     properties, errcode_ret);
  ldOclEnv.context = context;

  return ldOclEnv.command_queue;
}

cl_int clFinish (cl_command_queue command_queue) {
  return real_clFinish(command_queue);
}

#if ENABLE_KERNEL_PROFILING == 1
static void print_statistics(int force) {
  if (force || IS_MPI_MAIN()) {
    int i;

    pthread_mutex_lock(&stats_lock);

    for (i = 0; i < kernel_elt_count; i++) {
      struct ld_kernel_s *ldKernel = &kernel_map[i];

      if (ldKernel->exec_counter == 0) {
        continue;
      }

      info("Kernel %s:\n", ldKernel->name);
      info("\t- was executed %d times\n", ldKernel->exec_counter);

      unsigned long avg_ns = ldKernel->exec_span_ns / ldKernel->exec_counter;
      info("\t- it took %"PRId64"ms (%"PRId64"ns).\n", NS_TO_MS(ldKernel->exec_span_ns), ldKernel->exec_span_ns);
      info("\t- average time was %"PRId64"ms (%"PRId64"ns).\n", NS_TO_MS(avg_ns), avg_ns);
    }
    info("---------------\n");

    pthread_mutex_unlock(&stats_lock);
  }
}
#endif

cl_int clReleaseCommandQueue (cl_command_queue command_queue) {

#if ENABLE_KERNEL_PROFILING == 1
  print_statistics(0);
#endif
#if ENABLE_LEAK_DETECTION == 1
  if (IS_MPI_MAIN()) {
    int i;
    struct ld_program_s *ldprogram;
    struct ld_kernel_s *ldkernel;
    struct ld_mem_s *ldmem;
    //struct ld_mem_offset_s *ldmem_offset;

#define CHECK_ALL_RELEASED(_NAME)                               \
    FOR_ALL_MAP_ELTS(i, ld##_NAME, _NAME) {                     \
      if (!ld##_NAME->released) {                               \
        warning("Memory leak: " # _NAME                         \
                " %p not released. (#%d) %d\n",                 \
                ld##_NAME->handle, ld##_NAME->uid, ld##_NAME->released);     \
      }                                                         \
    }                                                           \

    CHECK_ALL_RELEASED(program);
    CHECK_ALL_RELEASED(kernel);
    CHECK_ALL_RELEASED(mem);
    //CHECK_ALL_RELEASED(mem_offset);
  }
#endif
  return real_clReleaseCommandQueue (command_queue);
}
/* ************************************************************************* */

cl_program clCreateProgramWithSource (cl_context context,
                                      cl_uint count,
                                      const char **strings,
                                      const size_t *lengths,
                                      cl_int *errcode_ret)
{
  struct ld_program_s *program = get_next_program_spot();

  program->handle = real_clCreateProgramWithSource(context, count, strings,
                                                   lengths, errcode_ret);

  if (!IS_MPI_MAIN()) {
    program_elt_count--;
    return program->handle;
  }

  ocl_handle_program(program->handle, count, strings, lengths);

  return program->handle;
}

cl_int clReleaseProgram (cl_program program) {
  struct ld_program_s *ldProgram = find_program_entry(program);

  assert(ldProgram);

  ldProgram->released = 1;

  return real_clReleaseProgram (program);
}

/* ************************************************************************* */

static inline ld_flags ocl_buffer_flags_to_ld (cl_mem_flags flags) {
  ld_flags ldFlags = 0;

  if (flags & CL_MEM_WRITE_ONLY) ldFlags |= LD_FLAG_WRITE_ONLY;
  if (flags & CL_MEM_READ_ONLY) ldFlags |= LD_FLAG_READ_ONLY;
  if (flags & CL_MEM_READ_WRITE) ldFlags |= LD_FLAG_READ_WRITE;

  return ldFlags;
}

struct ld_mem_s *create_ocl_buffer(cl_mem handle) {
  static unsigned int buffer_uid = -1;
  struct ld_mem_s *ldBuffer = get_next_mem_spot();

  ldBuffer->handle = handle;
  ldBuffer->uid = buffer_uid++;

  return ldBuffer;
}

cl_mem clCreateBuffer (cl_context context, cl_mem_flags flags, size_t size,
                       void *host_ptr, cl_int *errcode_ret)
{
  struct ld_mem_s *buffer;

  if (flags & CL_MEM_ALLOC_HOST_PTR) {
    goto unhandled;
  }

  buffer = create_ocl_buffer(real_clCreateBuffer(context, flags, size, host_ptr, errcode_ret));

  buffer->size = size;
  buffer->flags = ocl_buffer_flags_to_ld(flags);

  buffer_created_event(buffer);

  return buffer->handle;

unhandled:
  return real_clCreateBuffer(context, flags, size, host_ptr, errcode_ret);
}

/* ************************************************************************* */

cl_mem clCreateSubBuffer (cl_mem buffer,
                          cl_mem_flags flags,
                          cl_buffer_create_type buffer_create_type,
                          const void *buffer_create_info,
                          cl_int *errcode_ret)
{
  struct ld_mem_s *ldBuffer = find_mem_entry(buffer);
  struct ld_mem_offset_s *ldSubBuffer = NULL;
  cl_mem subbuffer;
  int i;

  assert(ldBuffer);

  subbuffer = real_clCreateSubBuffer(buffer, flags, buffer_create_type,
                                     buffer_create_info, errcode_ret);

  if (subbuffer == buffer) {
    return subbuffer;
  }
  for (i = 0; i < mem_offset_elt_count; i++) {
    if (mem_offset_map[i].released) {
      ldSubBuffer = &mem_offset_map[i];
    }
  }

  if (!ldSubBuffer) {
    ldSubBuffer = get_next_mem_offset_spot();
  }

  assert(ldSubBuffer);

  ldSubBuffer->handle = subbuffer;
  ldSubBuffer->flags = ocl_buffer_flags_to_ld(flags);
  ldSubBuffer->offset = ((cl_buffer_region *) buffer_create_info)->origin;
  ldSubBuffer->size = ((cl_buffer_region*) buffer_create_info)->size;
  ldSubBuffer->parent = ldBuffer;
  ldSubBuffer->released = 0;

  subbuffer_created_event(ldBuffer, ldSubBuffer->offset);

  return subbuffer;
}

/* ************************************************************************* */

cl_int clReleaseMemObject (cl_mem memobj) {
  struct ld_mem_offset_s *ldSubBuffer = find_mem_offset_entry(memobj);
  struct ld_mem_s *ldBuffer = find_mem_entry(memobj);

  if (ldSubBuffer) {
    ldSubBuffer->released = 1;
  }

  if (ldBuffer) {
    buffer_released(ldBuffer);
  }

  //ignore other buffer types

  return real_clReleaseMemObject (memobj);
}

/* ************************************************************************* */

cl_kernel clCreateKernel (cl_program  program,
                          const char *kernel_name,
                          cl_int *errcode_ret)
{
  int nb_params, i;
  char **types_and_names;

  struct ld_kernel_s *ldKernel = get_next_kernel_spot();

  ldKernel->handle = real_clCreateKernel(program, kernel_name, errcode_ret);
  ldKernel->name = kernel_name;

  if (!IS_MPI_MAIN()) {
    kernel_elt_count--;
    return ldKernel->handle;
  }

  types_and_names = ocl_handle_create_kernel(program, ldKernel->handle, kernel_name);
  for (nb_params = 0; types_and_names[nb_params]; nb_params++);
  nb_params /= 2;

  ldKernel->nb_params = nb_params;
  ldKernel->params = malloc(sizeof(struct ld_kern_param_s) * nb_params);

  for (i = 0; i < nb_params; i++) {
    ldKernel->params[i].name = types_and_names[i*2 + 1];
    ldKernel->params[i].type = types_and_names[i*2];
    ldKernel->params[i].index = i;
  }

  kernel_created_event(ldKernel);

  return ldKernel->handle;
}

cl_int clReleaseKernel (cl_kernel kernel) {
  struct ld_kernel_s *ldKernel = find_kernel_entry(kernel);

  assert(ldKernel);

  ldKernel->released = 1;

  return real_clReleaseKernel (kernel);
}

/* ************************************************************************* */

#define DEFAULT_SIZE(size) if (!size) size = 4;

int ocl_getBufferContent (struct ld_mem_s *ldBuffer, void *buffer,
                          size_t offset, size_t size)
{
  if (size == 0) {
    return 1;
  }

  //debug("*** Read %zub from buffer #%d at +%zub *** \n", size, ldBuffer->uid, offset);
  cl_int err = real_clEnqueueReadBuffer(ldOclEnv.command_queue,
                                        ldBuffer->handle, CL_TRUE,
                                        offset, size, buffer,
                                        0, NULL, NULL);
  assert(err == CL_SUCCESS);

  return err == CL_SUCCESS;
}

void *ocl_setParameterValue (struct ld_kernel_s *ldKernel,
                             struct ld_kern_param_s *ldParam,
                             void *buffer,  size_t size)
{
  cl_mem mem_obj = (void *) -1;
  cl_int errcode_ret;

  if (ldParam->is_pointer) {
    DEFAULT_SIZE(size)

    mem_obj = real_clCreateBuffer(ldOclEnv.context, CL_MEM_READ_WRITE, size, NULL, clck_(&errcode_ret));

    clCheck(real_clEnqueueWriteBuffer(ldOclEnv.command_queue,
                                      (cl_mem) mem_obj,
                                      CL_TRUE,
                                      0, size, buffer,
                                      0, NULL, NULL));
    buffer = &mem_obj;
    size = sizeof(cl_mem);
  }

  if (size == 0 && strstr(ldParam->name, "_tex")) {
    cl_image_format format = {CL_R, CL_UNSIGNED_INT32};
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    mem_obj = clCreateImage2D (ldOclEnv.context, CL_MEM_READ_ONLY, &format, 100, 1, 0, &format, clck_(&errcode_ret));
#pragma GCC diagnostic pop
    buffer = &mem_obj;
    size = sizeof(cl_mem);
  }

  clCheck(real_clSetKernelArg ((cl_kernel) ldKernel->handle, ldParam->index,
                               size, buffer));

  return mem_obj;

}

int ocl_triggerKernelExecution (struct ld_kernel_s *ldKernel,
                                const struct work_size_s *work_sizes,
                                unsigned int work_dim)
{
  static size_t global_work_size[MAX_WORK_DIM], local_work_size[MAX_WORK_DIM];
  unsigned int i;

  for (i = 0; i < work_dim; i++) {
    local_work_size[i] = work_sizes->local[i];
    global_work_size[i] = work_sizes->global[i] * work_sizes->local[i];
  }

  clCheck(real_clEnqueueNDRangeKernel(ldOclEnv.command_queue, (cl_kernel) ldKernel->handle,
                                      work_dim, NULL,
                                      global_work_size,
                                      local_work_size,
                                      0, NULL, NULL));
  return 1;
}

int ocl_getAndReleaseParameterValue (struct ld_kernel_s *ldKernel,
                                       struct ld_kern_param_s *ldParam,
                                       void *buffer_handle,
                                       void *buffer,  size_t size)
{
  if (buffer_handle == (void *) -1) {
    return 1;
  }

  if (size == 0) {
    goto do_release;
  }

  clCheck(real_clEnqueueReadBuffer(ldOclEnv.command_queue,
                                   (cl_mem) buffer_handle,
                                   CL_TRUE,
                                   0, size, buffer,
                                   0, NULL, NULL));
do_release:
  clCheck(real_clReleaseMemObject((cl_mem) buffer_handle));

  return 1;
}

cl_int clSetKernelArg (cl_kernel kernel,
                       cl_uint arg_index,
                       size_t arg_size,
                       const void *arg_value) {
  struct ld_kernel_s *ldKernel = find_kernel_entry(kernel);
  struct ld_kern_param_s *ldParam;

  assert(ldKernel);
  ldParam = &ldKernel->params[arg_index];

  if (ldParam->is_pointer) {
    struct ld_mem_s *ldBuffer = find_mem_entry(arg_value == NULL ? NULL : *(cl_mem *) arg_value);
    size_t offset = 0;

    if (!ldBuffer) {
      struct ld_mem_offset_s *ldSubBuffer = find_mem_offset_entry(*(cl_mem *) arg_value);

      assert(ldSubBuffer);
      ldBuffer = ldSubBuffer->parent;
      offset = ldSubBuffer->offset;
    }

    assert(ldBuffer);
    kernel_set_buffer_arg_event (ldKernel, ldParam, arg_index, ldBuffer, offset);
  } else {
    kernel_set_scalar_arg_event (ldKernel, ldParam, arg_index, (const void **) arg_value);
  }

  return real_clSetKernelArg(kernel, arg_index, arg_size, arg_value);
}

/* ************************************************************************* */

#if ENABLE_KERNEL_PROFILING == 1
static void CL_CALLBACK  kernel_profiler_cb (cl_event event,
                                             cl_int event_command_exec_status,
                                             void *user_data)
{
  static cl_ulong tstart, tstop, len;
  cl_int refcnt;
  struct ld_kernel_s *ldKernel = (struct ld_kernel_s *) user_data;

  pthread_mutex_lock(&stats_lock);
  clReleaseEvent(event);
  clCheck(clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(tstop), &tstop, NULL));
  clCheck(clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(tstart), &tstart, NULL));

  clCheck(clGetEventInfo(event,  CL_EVENT_REFERENCE_COUNT, sizeof(refcnt), &refcnt, NULL));

  len = tstop - tstart;
  if (tstart > tstop) {
    len = tstart - tstop;
  }

  if (tstart == 0ul || tstop == 0ul) {
    // invalid timestamps
    len = 0;
  }

  ldKernel->exec_span_ns += len;
  pthread_mutex_unlock(&stats_lock);
}
#endif

cl_int clEnqueueNDRangeKernel (cl_command_queue command_queue,
                               cl_kernel kernel,
                               cl_uint work_dim,
                               const size_t *global_work_offset,
                               const size_t *global_work_size,
                               const size_t *local_work_size,
                               cl_uint num_events_in_wait_list,
                               const cl_event *event_wait_list,
                               cl_event *event)
{
  static struct work_size_s work_sizes;

  struct ld_kernel_s *ldKernel = find_kernel_entry(kernel);
  int i;
  cl_int errcode;

  if (num_events_in_wait_list) {
    clCheck(clWaitForEvents(num_events_in_wait_list, event_wait_list));
  }

  assert(ldKernel);
  for (i = 0; i < work_dim; i++) {
    work_sizes.local[i] = local_work_size[i];
    work_sizes.global[i] = global_work_size[i]/work_sizes.local[i];
  }

#if ENABLE_KERNEL_PROFILING == 1
  static cl_event kern_event;

  if (!event) {
    event = &kern_event; // scope of the event is limited to this function.
  }
#endif

  kernel_executed_event(ldKernel, &work_sizes, work_dim);

  errcode = real_clEnqueueNDRangeKernel(command_queue, kernel, work_dim,
                                        global_work_offset, global_work_size,
                                        local_work_size, num_events_in_wait_list,
                                        event_wait_list, event);
#if ENABLE_KERNEL_PROFILING == 1
  clCheck(errcode);

  clRetainEvent(*event);
  clSetEventCallback(*event, CL_COMPLETE, kernel_profiler_cb, ldKernel);
#endif

#if FORCE_FINISH_KERNEL
  real_clFinish(command_queue);
#endif

  kernel_finished_event(ldKernel, &work_sizes, work_dim);

  return errcode;
}

/* ************************************************************************* */

struct ld_mem_s *readWriteMemory(cl_mem buffer, void **ptr, int direction,
                                 size_t size, size_t offset)
{
  struct ld_mem_s *ldBuffer = find_mem_entry(buffer);
  size_t real_offset = offset;

  if (!ldBuffer) {
    struct ld_mem_offset_s *ldSubBuffer = find_mem_offset_entry(buffer);

    assert(ldSubBuffer);
    ldBuffer = ldSubBuffer->parent;
    real_offset += ldSubBuffer->offset;
  }

  assert(ldBuffer);

  buffer_copy_event(ldBuffer, direction, (void **) ptr, size, real_offset);

  return ldBuffer;
}

cl_int clEnqueueWriteBuffer (cl_command_queue command_queue,
                             cl_mem buffer,
                             cl_bool blocking_write,
                             size_t offset,
                             size_t size,
                             const void *ptr,
                             cl_uint num_events_in_wait_list,
                             const cl_event *event_wait_list,
                             cl_event *event)
{
  readWriteMemory(buffer, (void **) ptr, LD_WRITE, size, offset);

  return real_clEnqueueWriteBuffer(command_queue, buffer, blocking_write,
                                   offset, size, ptr, num_events_in_wait_list,
                                   event_wait_list, event);
}

/* ************************************************************************* */

cl_int clEnqueueReadBuffer (cl_command_queue command_queue,
                            cl_mem buffer,
                            cl_bool blocking_read,
                            size_t offset,
                            size_t size,
                            void *ptr,
                            cl_uint num_events_in_wait_list,
                            const cl_event *event_wait_list,
                            cl_event *event)
{
  cl_int errcode;

  errcode = real_clEnqueueReadBuffer(command_queue, buffer, CL_TRUE /* blocking_read */,
                                     offset, size, ptr, num_events_in_wait_list,
                                     event_wait_list, event);

  readWriteMemory(buffer, ptr, LD_READ, size, offset);


  return errcode;
}


/* The OpenCL Extension Wrangler Library
 * https://code.google.com/p/clew/
 * MIT License
 * */

const char* clewErrorString (cl_int error) {
  static const char* strings[] = {
    // Error Codes
    "CL_SUCCESS"                                  //   0
    , "CL_DEVICE_NOT_FOUND"                         //  -1
    , "CL_DEVICE_NOT_AVAILABLE"                     //  -2
    , "CL_COMPILER_NOT_AVAILABLE"                   //  -3
    , "CL_MEM_OBJECT_ALLOCATION_FAILURE"            //  -4
    , "CL_OUT_OF_RESOURCES"                         //  -5
    , "CL_OUT_OF_HOST_MEMORY"                       //  -6
    , "CL_PROFILING_INFO_NOT_AVAILABLE"             //  -7
    , "CL_MEM_COPY_OVERLAP"                         //  -8
    , "CL_IMAGE_FORMAT_MISMATCH"                    //  -9
    , "CL_IMAGE_FORMAT_NOT_SUPPORTED"               //  -10
    , "CL_BUILD_PROGRAM_FAILURE"                    //  -11
    , "CL_MAP_FAILURE"                              //  -12

    , ""    //  -13
    , ""    //  -14
    , ""    //  -15
    , ""    //  -16
    , ""    //  -17
    , ""    //  -18
    , ""    //  -19

    , ""    //  -20
    , ""    //  -21
    , ""    //  -22
    , ""    //  -23
    , ""    //  -24
    , ""    //  -25
    , ""    //  -26
    , ""    //  -27
    , ""    //  -28
    , ""    //  -29

    , "CL_INVALID_VALUE"                            //  -30
    , "CL_INVALID_DEVICE_TYPE"                      //  -31
    , "CL_INVALID_PLATFORM"                         //  -32
    , "CL_INVALID_DEVICE"                           //  -33
    , "CL_INVALID_CONTEXT"                          //  -34
    , "CL_INVALID_QUEUE_PROPERTIES"                 //  -35
    , "CL_INVALID_COMMAND_QUEUE"                    //  -36
    , "CL_INVALID_HOST_PTR"                         //  -37
    , "CL_INVALID_MEM_OBJECT"                       //  -38
    , "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR"          //  -39
    , "CL_INVALID_IMAGE_SIZE"                       //  -40
    , "CL_INVALID_SAMPLER"                          //  -41
    , "CL_INVALID_BINARY"                           //  -42
    , "CL_INVALID_BUILD_OPTIONS"                    //  -43
    , "CL_INVALID_PROGRAM"                          //  -44
    , "CL_INVALID_PROGRAM_EXECUTABLE"               //  -45
    , "CL_INVALID_KERNEL_NAME"                      //  -46
    , "CL_INVALID_KERNEL_DEFINITION"                //  -47
    , "CL_INVALID_KERNEL"                           //  -48
    , "CL_INVALID_ARG_INDEX"                        //  -49
    , "CL_INVALID_ARG_VALUE"                        //  -50
    , "CL_INVALID_ARG_SIZE"                         //  -51
    , "CL_INVALID_KERNEL_ARGS"                      //  -52
    , "CL_INVALID_WORK_DIMENSION"                   //  -53
    , "CL_INVALID_WORK_GROUP_SIZE"                  //  -54
    , "CL_INVALID_WORK_ITEM_SIZE"                   //  -55
    , "CL_INVALID_GLOBAL_OFFSET"                    //  -56
    , "CL_INVALID_EVENT_WAIT_LIST"                  //  -57
    , "CL_INVALID_EVENT"                            //  -58
    , "CL_INVALID_OPERATION"                        //  -59
    , "CL_INVALID_GL_OBJECT"                        //  -60
    , "CL_INVALID_BUFFER_SIZE"                      //  -61
    , "CL_INVALID_MIP_LEVEL"                        //  -62
    , "CL_INVALID_GLOBAL_WORK_SIZE"                 //  -63
    , "CL_UNKNOWN_ERROR_CODE"
  };

  if (error >= -63 && error <= 0) {
    return strings[-error];
  } else {
    return strings[64];
  }
}
