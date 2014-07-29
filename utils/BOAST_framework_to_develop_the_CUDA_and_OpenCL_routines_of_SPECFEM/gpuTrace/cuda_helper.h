struct kernel_lookup_s {
  void *address;
  const char *name;
  size_t nb_params;
  struct param_info_s *params;
};

struct param_info_s {
  const char *name;
  const char *type;
};

void cuda_init_helper(void);

int cuda_get_nb_kernels(void);

struct kernel_lookup_s *cuda_get_lookup_table(void);
