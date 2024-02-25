

void ocl_init_helper(void);

char **ocl_handle_create_kernel(void *program, void *kernel, const char *name);

void ocl_handle_program(void *program,
                    unsigned int count,
                    const char **strings,
                    const size_t *lengths);
