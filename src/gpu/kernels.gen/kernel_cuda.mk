cuda_kernels_OBJS := \
	$O/outer_core_impl_kernel_forward.cuda-kernel.o \
	$O/outer_core_impl_kernel_adjoint.cuda-kernel.o \
	$O/inner_core_impl_kernel_forward.cuda-kernel.o \
	$O/inner_core_impl_kernel_adjoint.cuda-kernel.o \
	$O/compute_rho_kernel.cuda-kernel.o \
	$O/compute_iso_kernel.cuda-kernel.o \
	$O/compute_ani_kernel.cuda-kernel.o \
	$O/compute_hess_kernel.cuda-kernel.o \
	$O/compute_acoustic_kernel.cuda-kernel.o \
	$O/compute_strength_noise_kernel.cuda-kernel.o \
	$O/crust_mantle_impl_kernel_forward.cuda-kernel.o \
	$O/crust_mantle_impl_kernel_adjoint.cuda-kernel.o \
	$O/compute_ani_undo_att_kernel.cuda-kernel.o \
	$O/compute_iso_undo_att_kernel.cuda-kernel.o \
	$(EMPTY_MACRO)
