      type_f = BOAST::Real::new.decl
      if BOAST::get_default_real_size == 8 then
        type_i = "unsigned long int"
        cmpx_name = "atom_cmpxchg"
      else
        type_i = "unsigned int"
        cmpx_name = "atomic_cmpxchg"
      end
      BOAST::get_output.print <<EOF
inline void atomicAdd(volatile __global float *source, const float val) {
  union {
    #{type_i} iVal;
    #{type_f} fVal;
  } res, orig;
  do {
    orig.fVal = *source;
    res.fVal = orig.fVal + val;
  } while (#{cmpx_name}((volatile __global #{type_i} *)source, orig.iVal, res.iVal) != orig.iVal);
}
EOF

