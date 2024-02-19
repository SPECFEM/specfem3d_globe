#ifndef LD_CHECKER_H
#define LD_CHECKER_H

#include <stdio.h>
#include <string.h>
#include <ctype.h>

/******************************
 * Configuration
 ******************************/

#define PRINT_BUFFER_TRANSFER_FIRST_BYTES_AS_FLOAT 1

#define BUFFER_ZERO_IS_NULL 0

#define FORCE_FINISH_KERNEL 1


/******************************
 * Select what and when info
 * are traced.
 ******************************/

#define PRINT_KERNEL_BUFFER_CREATION 0
#define PRINT_BUFFER_DIRECTION 0
#define PRINT_BUFFER_TRANSFER 0
#define PRINT_BUFFER_RELEASE 0

#define PRINT_KERNEL_BEFORE_EXEC 1
#define PRINT_KERNEL_AFTER_EXEC 1
#define PRINT_KERNEL_AFTER_EXEC_IGNORE_CONST 0
#define PRINT_KERNEL_NAME_ONLY 0

/******************************
 * Filters for kernel execution
 ******************************/

#define USE_ADVANCED_KERNEL_FILTER 1
#define ENV__KERNEL_FILTER "LD_KERNEL_FILTER"

#define FILTER_BY_KERNEL_EXEC_CPT 1
#define KERNEL_EXEC_CPT_LOWER_BOUND 999
#define KERNEL_EXEC_CPT_UPPER_BOUND 999

#define FILTER_BY_KERNEL_NAME 1
#define KERNEL_NAME_FILTER "write_seismograms_transfer_from_device_kernel"

/******************************
 * Print full buffers to screen
 * or into a file.
 ******************************/

#define PRINT_KERNEL_PARAMS_TO_FILE 0

#define PRINT_KERNEL_ARG_FULL_BUFFER 0
#define FULL_BUFFER_SIZE_LIMIT 0

#define PARAM_FILE_DIRECTORY "kernel_test_data"
#define ENV__PARAM_FILE_SUBDIR_PREFIX "LD_KERNEL_SUBDIR_PREFIX"
#define PARAM_FILE_DIRECTORY_PERM 0777

#define ENABLE_KERNEL_TESTING 1

/******************************
 * MPI support
 ******************************/

#define WITH_MPI 0
#define ONLY_MPI_ROOT_OUTPUT 1

/******************************
 * Execution monitoring
 ******************************/

#define ENABLE_KERNEL_PROFILING 0
#define ENABLE_LEAK_DETECTION 0

/******************************
 ******************************
 ******************************/

typedef int ld_flags;

#define LD_FLAG_READ_ONLY 1
#define LD_FLAG_WRITE_ONLY 2
#define LD_FLAG_READ_WRITE 4

struct ld_bindings_s {
  const char *name;
  void **real_address;
};

enum type_info_type_e {
  TYPE_INFO_DEFAULT,
  TYPE_INFO_FLOAT,
  TYPE_INFO_DOUBLE,
  TYPE_INFO_INT,
  TYPE_INFO_UINT
};

extern struct type_info_s {
  const char *type_name;
  const char *format;
  const unsigned int size;
  const enum type_info_type_e type;
} TYPE_FORMATTERS[] ;

/** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **/

#define FIRST_BYTES_TO_READ 32
struct ld_mem_s {
  void *handle;
  unsigned int uid;
  size_t size;
  ld_flags flags;
  int has_values;
  int values_outdated;
  unsigned int released;
  char first_values[FIRST_BYTES_TO_READ];
};

#define CURRENT_VALUE_BUFFER_SZ 80
struct ld_kern_param_s {
  const char *name;
  const char *type;
  const struct type_info_s *type_info;
  int is_pointer;
  int has_current_value;
  struct ld_mem_s *current_buffer;
  size_t offset;
  size_t index;
  char current_value[CURRENT_VALUE_BUFFER_SZ];
  char current_binary_value[10];
};

struct ld_kernel_s {
  unsigned int uid;
  void *handle;
  int nb_params;
  const char *name;
  struct ld_kern_param_s *params;
  unsigned int exec_counter;
  unsigned int released;
#if ENABLE_KERNEL_PROFILING == 1
  unsigned long long int exec_span_ns;
#endif
};

#define MAX_WORK_DIM 3
struct work_size_s {
  size_t local[MAX_WORK_DIM];
  size_t global[MAX_WORK_DIM];
};

/** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **/

#if WITH_MPI == 1 && ONLY_MPI_ROOT_OUTPUT == 1
extern int myrank;
#define ONLY_MPI_MAIN()      \
  if (myrank != -1 && myrank != 0) { \
    return;                          \
  }
#define IS_MPI_MAIN() (myrank == -1 || myrank == 0)
#else
#define ONLY_MPI_MAIN()
#define IS_MPI_MAIN() 1
#endif

/** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **/

void dbg_crash_event(void);
void dbg_notify_event(void);

struct callback_s {
  int (*getBufferContent) (struct ld_mem_s *ldBuffer, void *buffer,
                           size_t offset, size_t size);
  void *(*setParameterValue) (struct ld_kernel_s *ldKernel,
                              struct ld_kern_param_s *ldParam,
                              void *buffer,  size_t size);
  int (*getAndReleaseParameterValue) (struct ld_kernel_s *ldKernel,
                                      struct ld_kern_param_s *ldParam,
                                      void *buffer_handle,
                                      void *buffer, size_t size);
  int (*triggerKernelExecution) (struct ld_kernel_s *ldKernel,
                                 const struct work_size_s *work_sizes,
                                 unsigned int work_dim);
};

void debug(const char *format, ...);
void info(const char *format, ...);
void warning(const char *format, ...);
void error(const char *format, ...);
void gpu_trace(const char *format, ...);

/*
 * struct ld_name_s *find_name_entry(cl_key_type key);
 * struct ld_name_s *get_next_name_spot();
 * int add_to_name_map (cl_key_type key);
 * struct ld_name_s *get_name_entry(int id));
 *
 * */

#define assert(val) {if (!val) {error ("pointer " #val " is null at %s:%d %s\n", \
                                       __FILE__, __LINE__,  __func__);}}

#define CREATE_HASHMAP(_NAME, _KEY_TYPE,_MAX_KEYS)                      \
  static struct ld_##_NAME##_s _NAME##_map[_MAX_KEYS];                  \
  static unsigned int _NAME##_elt_count = 0;                            \
                                                                        \
  static int find_##_NAME##_entry_id (_KEY_TYPE key);                   \
  static struct ld_##_NAME##_s *get_##_NAME##_entry (int id);           \
  static struct ld_##_NAME##_s *find_##_NAME##_entry(_KEY_TYPE key);    \
  static struct ld_##_NAME##_s *get_next_##_NAME##_spot(void);          \
  static int add_to_##_NAME##_map (_KEY_TYPE key);                      \
                                                                        \
  int find_##_NAME##_entry_id (_KEY_TYPE key) {                         \
    int i;                                                              \
    for (i = 0; i < _NAME##_elt_count; i++) {                           \
      if (_NAME##_map[i].handle == key) {                               \
        return i ;                                                      \
      }                                                                 \
    }                                                                   \
    return -1;                                                          \
  }                                                                     \
                                                                        \
  static inline struct ld_##_NAME##_s *get_##_NAME##_entry (int id) {   \
    return (id < _NAME##_elt_count && id >= 0 ?                         \
            &_NAME##_map[id] : NULL);                                   \
  }                                                                     \
                                                                        \
  static inline struct ld_##_NAME##_s *find_##_NAME##_entry(_KEY_TYPE key) { \
    return get_##_NAME##_entry(find_##_NAME##_entry_id (key));          \
  }                                                                     \
                                                                        \
  static inline struct ld_##_NAME##_s *get_next_##_NAME##_spot(void) {  \
    if (_NAME##_elt_count >= _MAX_KEYS) {                               \
      return NULL;                                                      \
    }                                                                   \
    _NAME##_map[_NAME##_elt_count].uid =_NAME##_elt_count;              \
    return &_NAME##_map[_NAME##_elt_count++];                           \
  }                                                                     \
                                                                        \
  static inline int add_to_##_NAME##_map (_KEY_TYPE key) {              \
    struct ld_##_NAME##_s *next;                                        \
      next = get_next_##_NAME##_spot();                                 \
        if (!next) {                                                    \
          return -1;                                                    \
        }                                                               \
        next->handle = key;                                             \
        return _NAME##_elt_count - 1;   /* next - start / size ? */     \
  }
#define FOR_ALL_MAP_ELTS(_CPT, _VAR, _NAME)                             \
  for (_CPT = 0, _VAR = &_NAME##_map[_CPT];                             \
       _CPT < _NAME##_elt_count;                                        \
       _CPT++, _VAR = &_NAME##_map[_CPT])

/** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **/

static inline char *buffer_flag_to_string(ld_flags flags) {
  switch(flags) {
  case LD_FLAG_READ_ONLY: return "READ_WRITE";
  case LD_FLAG_WRITE_ONLY: return "WRITE_ONLY";
  case LD_FLAG_READ_WRITE: return "READ_WRITE";
  default: return "UNKNWON";
  }
}

static inline char *buffer_flag_to_direction(ld_flags flags) {
#if  PRINT_BUFFER_DIRECTION
  switch(flags) {
  case LD_FLAG_READ_ONLY: return "--->";
  case LD_FLAG_WRITE_ONLY: return "<---";
  case LD_FLAG_READ_WRITE: return "<-->";
  default: return "?--?";
  }
#else
  return "";
#endif
}

void init_ldchecker(struct callback_s callbacks, struct ld_bindings_s *lib_bindings);
void buffer_created_event(struct ld_mem_s *buffer);
void subbuffer_created_event(struct ld_mem_s *buffer, size_t offset);
#define LD_WRITE 0
#define LD_READ 1
void buffer_copy_event(struct ld_mem_s *buffer, int is_read, void **ptr,
                       size_t size, size_t offset);
void buffer_released (struct ld_mem_s *ldBuffer);
void kernel_created_event(struct ld_kernel_s *kernel);
void kernel_set_scalar_arg_event (struct ld_kernel_s *ldKernel,
                                  struct ld_kern_param_s *ldParam,
                                  int arg_index,
                                  const void **arg_value);
void kernel_set_buffer_arg_event (struct ld_kernel_s *ldKernel,
                                  struct ld_kern_param_s *ldParam,
                                  int arg_index,
                                  struct ld_mem_s *ldBuffer,
                                  size_t offset);
void kernel_executed_event(struct ld_kernel_s *kernel,
                           const struct work_size_s *work_sizes,
                           unsigned int work_dim);
void kernel_finished_event(struct ld_kernel_s *ldKernel,
                           const struct work_size_s *work_sizes,
                           int work_dim);
/** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** ** **/

static inline const struct type_info_s *get_type_info (const char *type_name) {
  int i;

  for (i = 0; TYPE_FORMATTERS[i].type_name
         && (!type_name || !strstr(type_name, TYPE_FORMATTERS[i].type_name)); i++);

  return &TYPE_FORMATTERS[i];
}

static inline int is_pointer_type (const char *type_name) {
  int i;

  if (!type_name) {
    return 0;
  }

  for(i = 0; type_name[i]; i++) {
    if (type_name[i] == '*') {
      return 1;
    }
  }
  return 0;
}

#endif

