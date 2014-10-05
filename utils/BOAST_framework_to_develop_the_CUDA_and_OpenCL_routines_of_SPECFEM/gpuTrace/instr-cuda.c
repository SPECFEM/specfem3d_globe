#include <stdio.h>
#include <stdlib.h>

#include <cuda_runtime.h>

#include "ldChecker.h"

cudaError_t (*real_cudaGetDeviceCount)(int *count);
cudaError_t (*real_cudaMalloc)  (void **devPtr, size_t size);
cudaError_t (*real_cudaMemcpy)  (void *dst, const void *src,
                                 size_t count, enum cudaMemcpyKind kind);
cudaError_t (*real_cudaMemcpyToSymbol)(const void *symbol, const void *src,
                                       size_t count, size_t offset,
                                       enum cudaMemcpyKind kind );
cudaError_t (*real_cudaMemcpyAsync)(void *dst, const void * src, size_t count,
                                    enum cudaMemcpyKind kind, cudaStream_t stream);
cudaError_t (*real_cudaMemset) (void *devPtr, int  value, size_t  count);
cudaError_t (*real_cudaConfigureCall) (dim3 gridDim, dim3 blockDim,
                                       size_t sharedMem, cudaStream_t stream);
cudaError_t (*real_cudaSetupArgument) (const void *arg, size_t size, size_t offset);
cudaError_t (*real_cudaLaunch) (const void *entry);
cudaError_t (*real_cudaGetSymbolAddress) (void **devPtr, const void *symbol);
cudaError_t (*real_cudaFree) (void *devPtr);
cudaError_t (*real_cudaStreamCreate) (cudaStream_t *stream);
cudaError_t (*real_cudaStreamSynchronize) (cudaStream_t stream);
cudaError_t (*real_cudaStreamDestroy) (cudaStream_t stream);

struct ld_bindings_s cuda_bindings[] = {
  {"cudaGetDeviceCount", (void **) &real_cudaGetDeviceCount},
  {"cudaMalloc", (void **) &real_cudaMalloc},
  {"cudaMemcpy", (void **) &real_cudaMemcpy},
  {"cudaMemcpyToSymbol", (void **) &real_cudaMemcpyToSymbol},
  {"cudaMemcpyAsync", (void **) &real_cudaMemcpyAsync},
  {"cudaMemset", (void **) &real_cudaMemset},
  {"cudaConfigureCall", (void **) &real_cudaConfigureCall},
  {"cudaSetupArgument", (void **) &real_cudaSetupArgument},
  {"cudaLaunch", (void **) &real_cudaLaunch},
  {"cudaGetSymbolAddress", (void **) &real_cudaGetSymbolAddress},
  {"cudaFree", (void **) &real_cudaFree},
  {"cudaStreamCreate", (void **) &real_cudaStreamCreate},
  {"cudaStreamSynchronize", (void **) &real_cudaStreamSynchronize},
  {"cudaStreamDestroy", (void **) &real_cudaStreamDestroy},
  {NULL, NULL}
};

#include "cuda_helper.h"

static struct kernel_lookup_s *cuda_lookup_table;

int cuda_getBufferContent (struct ld_mem_s *ldBuffer, void *buffer,
                           size_t offset, size_t size);
static struct work_size_s *configure_get_worksizes(dim3 *gridDim, dim3 *blockDim);

CREATE_HASHMAP(mem, void*, 200)
CREATE_HASHMAP(kernel, void*, 100)

struct ld_mem_s *create_cuda_buffer(const void *devPtr, size_t size);

static void init_cuda_ldchecker(void) {
  static int inited = 0;
  struct callback_s callbacks = {cuda_getBufferContent};
  int i, j;

  if (inited) {
    return;
  }

  init_ldchecker(callbacks, cuda_bindings);
  cuda_init_helper();

  inited = 1;

  cuda_lookup_table = cuda_get_lookup_table();
  for (i = 0; cuda_lookup_table[i].address != NULL; i++) {
    struct kernel_lookup_s *kernel = &cuda_lookup_table[i];
    struct ld_kernel_s *ldKernel = get_next_kernel_spot();

    ldKernel->handle = kernel;
    ldKernel->name = kernel->name;
    ldKernel->nb_params = kernel->nb_params;
    ldKernel->params = malloc(sizeof(struct ld_kern_param_s) * ldKernel->nb_params);

    for (j = 0; j < ldKernel->nb_params; j++) {
      ldKernel->params[j].name = kernel->params[j].name;
      ldKernel->params[j].type = kernel->params[j].type;
      ldKernel->params[i].index = j;
    }

    kernel_created_event(ldKernel);
  }

  create_cuda_buffer(NULL, 0);
}

#define ELT_SIZE 4
struct ld_mem_s *find_off_mem_entry(const void *handle, size_t *offset) {
  struct ld_mem_s *ldBuffer = find_mem_entry((void *) handle);
  int i;

  if (ldBuffer) {
    return ldBuffer;
  }

  *offset = -1;
  for (i = 0; i < mem_elt_count; i++) {
    struct ld_mem_s *current = &mem_map[i];
    size_t distance = ((float *) handle - (float *) current->handle) * ELT_SIZE;

    if (distance > 0 &&  distance < current->size) {
      if (*offset == -1 || distance < *offset) {
        if (*offset != -1) {
          warning("find_off_mem_entry found multiple matching buffers "
                  " ptr=%p, found=%p+%zu\n", handle, ldBuffer->handle,
                  *offset);
        }
        *offset = distance;
        ldBuffer = current;
      }
    }
  }

  return ldBuffer;
}

cudaError_t cudaGetDeviceCount (int *count) {
  init_cuda_ldchecker();

  return real_cudaGetDeviceCount (count);
}

int cuda_getBufferContent (struct ld_mem_s *ldBuffer, void *buffer,
                           size_t offset, size_t size)
{
  cudaError_t err = cudaSuccess;

  if (offset != 0) {
    /* need to use the type to multiply correctly the buffer handle.  */
    warning("offset == %zu != 0, not handled yet in Cuda's ldChecker\n",
            offset);
  }

  err = real_cudaMemcpy (buffer,
                         ldBuffer->handle,
                         size,
                         cudaMemcpyDeviceToHost);

  return err == cudaSuccess;
}

static struct {
  const void *arg;
  size_t size;
  size_t offset;
} temp_arguments[90];
static int temp_argument_cnt = 0;

cudaError_t cudaLaunch (const void *entry) {
  struct ld_kernel_s *ldKernel;
  struct kernel_lookup_s *kernel;
  cudaError_t err;
  int i;

  for (i = 0; cuda_lookup_table[i].address != NULL; i++) {
    if (cuda_lookup_table[i].address == entry) {
      break;
    }
  }
  if (!cuda_lookup_table[i].address) {
    warning("Couldn't find kernel @%p\n", entry);
    goto error;
  }

  kernel = &cuda_lookup_table[i];
  ldKernel = find_kernel_entry(kernel);
  assert(ldKernel);

  dbg_notify_event();

  for (i = 0; i < ldKernel->nb_params; i++) {
    struct ld_kern_param_s *ldParam = &ldKernel->params[i];
    const void **current_arg = (const void **) temp_arguments[i].arg;

    if (ldParam->is_pointer) {
      size_t offset = 0;
      struct ld_mem_s *ldBuffer = find_off_mem_entry(*current_arg, &offset);
      /* better handling of invalid buffers ?*/
      if (!ldBuffer) {
        error("in kernel %s, arg #%d, invalid bufer %s@%p", ldKernel->name, i, ldParam->name, *current_arg);
      }
      assert(ldBuffer);
      kernel_set_buffer_arg_event (ldKernel, ldParam, i, ldBuffer, offset);
  } else {
      kernel_set_scalar_arg_event (ldKernel, ldParam, i, current_arg);
    }
  }

  kernel_executed_event(ldKernel, configure_get_worksizes(NULL, NULL), 3);
  err = real_cudaLaunch(entry);
  kernel_finished_event(ldKernel, configure_get_worksizes(NULL, NULL), 3);

  return err;

error:
  return real_cudaLaunch(entry);
}

cudaError_t cudaSetupArgument (const void *arg, size_t size, size_t offset) {
  temp_arguments[temp_argument_cnt].arg = arg;
  temp_arguments[temp_argument_cnt].size = size;
  temp_arguments[temp_argument_cnt].offset = offset;
  temp_argument_cnt++;

  return real_cudaSetupArgument(arg, size, offset);
}

struct ld_mem_s *create_cuda_buffer(const void *devPtr, size_t size) {
  static int uid = 0;
  struct ld_mem_s *ldBuffer = get_next_mem_spot();

  assert(ldBuffer);
  ldBuffer->handle = (void *) devPtr;
  ldBuffer->uid = uid++;
  ldBuffer->size = size;
  ldBuffer->flags = 0;

  ldBuffer->has_values = 0;
  ldBuffer->values_outdated = 0;

  buffer_created_event(ldBuffer);

  return ldBuffer;
}

cudaError_t cudaStreamCreate (cudaStream_t *stream) {
  return real_cudaStreamCreate(stream);
}

cudaError_t cudaStreamSynchronize (cudaStream_t stream) {
  return real_cudaStreamSynchronize(stream);
}

cudaError_t cudaStreamDestroy (cudaStream_t stream) {
  return real_cudaStreamDestroy(stream);
}


cudaError_t cudaMalloc (void **devPtr, size_t size) {
  cudaError_t retcode;

  retcode = real_cudaMalloc(devPtr, size);

  create_cuda_buffer(*devPtr, size);

  return retcode;
}

cudaError_t cudaMemcpyToSymbol(const void *symbol, const void *src, size_t count, size_t offset , enum cudaMemcpyKind kind ) {

  if (kind == cudaMemcpyHostToDevice) {
    struct ld_mem_s *ldBuffer;
    void *devPtr;
    cudaError_t err = real_cudaGetSymbolAddress(&devPtr, symbol);

    if (err != cudaSuccess) {
      warning("cudaMemcpyToSymbol failed ...\n");
      goto finish;
    }

    ldBuffer = find_mem_entry(devPtr);

    if (!ldBuffer) {
      ldBuffer = create_cuda_buffer(devPtr, count + offset);
    }
    buffer_copy_event(ldBuffer, LD_WRITE, (void **) src, count, offset);

  } else {
    warning("cudaMemcpyToSymbol kind=%d not handled\n", kind);
  }

 finish:
  return real_cudaMemcpyToSymbol(symbol, src, count, offset, kind);
}

cudaError_t cudaGetSymbolAddress (void **devPtr, const void *symbol) {
  cudaError_t err;

  err = real_cudaGetSymbolAddress(devPtr, symbol);

  return err;
}

cudaError_t cudaMemset ( void *devPtr, int  value, size_t  count) {
  cudaError_t retcode = real_cudaMemset(devPtr, value, count);
  size_t offset;
  struct ld_mem_s *ldBuffer = find_off_mem_entry(devPtr, &offset);

  warning("cudaMemset %d (unhandled)\n", ldBuffer->uid);

  return retcode;
}

void handle_cudaMemcpy(void *dst, const void * src, size_t count,
                       enum cudaMemcpyKind kind, cudaStream_t stream)
{
  struct ld_mem_s *ldBuffer;
  void *pointer;
  int direction;
  size_t offset = 0;

  if (kind == cudaMemcpyHostToDevice) /* write */ {
    ldBuffer = find_off_mem_entry(dst, &offset);

    assert(ldBuffer);
    pointer = (void *) src;
    direction = LD_WRITE;
  } else if (kind == cudaMemcpyDeviceToHost) /* read */ {
    ldBuffer = find_off_mem_entry((void *) src, &offset);

    assert(ldBuffer);
    pointer = dst;
    direction = LD_READ;
  } else {
    warning("Unhandled memcopy\n");
    return;
  }

  if (!ldBuffer) {
    warning("Coundln't find buffer structure, maybe it's offsetted, "\
            "and it's time to handle it ?\n");
    return;
  }

  buffer_copy_event(ldBuffer, direction, (void **) pointer, count, offset);
}

cudaError_t cudaMemcpyAsync(void *dst, const void * src, size_t count,
                            enum cudaMemcpyKind kind, cudaStream_t stream)
{
  cudaError_t retcode = real_cudaMemcpyAsync(dst, src, count, kind, stream);

  handle_cudaMemcpy(dst, src, count, kind, stream);

  return retcode;
}

cudaError_t cudaMemcpy (void *dst, const void *src,
                        size_t count, enum cudaMemcpyKind kind)
{
  cudaError_t retcode = real_cudaMemcpy(dst, src, count, kind);

  handle_cudaMemcpy(dst, src, count, kind, (void *) -1);

  return retcode;
}

cudaError_t cudaFree (void *devPtr) {
  struct ld_mem_s *ldBuffer = find_mem_entry(devPtr);

  buffer_released(ldBuffer);

  return real_cudaFree(devPtr);
}

static struct work_size_s *configure_get_worksizes(dim3 *gridDim, dim3 *blockDim) {
  static struct work_size_s work_sizes;
  if (gridDim != NULL) {
    work_sizes.global[0] = gridDim->x;
    work_sizes.global[1] = gridDim->y;
    work_sizes.global[2] = gridDim->z;
  }

  if (blockDim != NULL) {
    work_sizes.local[0] = blockDim->x;
    work_sizes.local[1] = blockDim->y;
    work_sizes.local[2] = blockDim->z;
  }

  return &work_sizes;
}


cudaError_t cudaConfigureCall (dim3 gridDim, dim3 blockDim,
                               size_t sharedMem, cudaStream_t stream)
{
  cudaError_t errcode;

  configure_get_worksizes(&gridDim, &blockDim);

  temp_argument_cnt = 0;

  errcode = real_cudaConfigureCall(gridDim, blockDim, sharedMem, stream);

  return errcode;
  //return 1 ? cudaSuccess : errcode;
}
