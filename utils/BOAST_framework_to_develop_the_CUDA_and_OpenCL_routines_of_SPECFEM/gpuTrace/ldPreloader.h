struct ld_bindings_s {
  const char *name;
  void **real_address;
} ;

extern struct ld_bindings_s bindings[];
