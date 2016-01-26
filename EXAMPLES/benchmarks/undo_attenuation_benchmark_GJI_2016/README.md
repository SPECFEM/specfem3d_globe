Undo Attenuation Benchmark
==========================

This benchmark reproduces the results of section 3 of the paper "Exact
time-reversal viscoelastic sensitivity kernels with parsimonious storage for
full waveform inversion" by D. Komatitsch, Z. Xie, E. Bozdag, D. B. Peter, E.
Sales de Andrade, Q. Liu and J. Tromp, published in GJI 2016.

There are three stages to the benchmark:
  * Exact computation of adjoint kernel using a full dump of the attenuated
    forward field. As this method requires the full field at all time steps,
    the storage requirements are quite large (~16 TiB) and it is only
    implemented for the alpha kernel. On the other hand, since attenuation does
    not need to be reversed, the result is exact.
  * Calculation of the adjoint kernel using a forward field with partial
    physical dispersion only. Since true attenuation is not used, the result is
    not exact, but it avoids the instability with solving the backward
    equations.
  * Calculation of the adjoint kernel using the newly proposed undo attenuation
    method. Via checkpointing, this method achieves a balance between excessive
    storage space and excessive computation, and is mathematically exact.
