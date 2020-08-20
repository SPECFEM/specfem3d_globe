require 'rubygems'
require 'BOAST'
require 'narray'
require './HEADER.rb'
require './FUNCTIONS.rb'

def rndup( val, div)
  return (val%div) == 0 ? val : val + div - (val%div)
end

$options = {:output_dir => "./output", :elem_per_thread => 1, :langs => [:CUDA, :CL] }

$parser = OptionParser::new do |opts|
  opts.on("-c","--check","Check kernels by building them") {
    $options[:check] = true
  }
  opts.on("-d","--display","Display kernels") {
    $options[:display] = true
  }
  opts.on("-v","--verbose","Verbosity") {
    $options[:verbose] = true
  }
  opts.on("-o","--output-dir DIR","Output directory") { |dir|
    $options[:output_dir] = dir
  }
  opts.on("-p","--platform PLATFORM","Selected Platform") { |platform|
    $options[:platform] = platform
  }
  opts.on("-k","--kernel KERNEL","Selected kernel") { |kernel|
    $options[:kernel] = kernel
  }
  opts.on("-e","--elem ELEM_PER_THREAD","Treat several elements in big kernels") { |elem_per_thread|
    $options[:elem_per_thread] = elem_per_thread.to_i
  }
  opts.on("-l","--lang LANG","Select language to use (CUDA or CL)") { |lang|
     $options[:langs] = [ lang.to_sym ]
  }
  opts.parse!
end

small_kernels = [
:assemble_boundary_accel_on_device,
:assemble_boundary_potential_on_device,
:prepare_boundary_potential_on_device,
:prepare_boundary_accel_on_device,
:get_maximum_scalar_kernel,
:get_maximum_vector_kernel,
:compute_add_sources_adjoint_kernel,
:compute_add_sources_kernel,
:compute_coupling_fluid_CMB_kernel,
:compute_coupling_fluid_ICB_kernel,
:compute_coupling_CMB_fluid_kernel,
:compute_coupling_ICB_fluid_kernel,
:compute_coupling_ocean_kernel,
:write_seismograms_transfer_from_device_kernel,
:write_seismograms_transfer_strain_from_device_kernel,
:noise_transfer_surface_to_host_kernel,
:noise_add_source_main_rec_kernel,
:noise_add_surface_movie_kernel,
:compute_seismograms_kernel,
:compute_stacey_acoustic_kernel,
:compute_stacey_acoustic_backward_kernel,
:compute_stacey_elastic_kernel,
:compute_stacey_elastic_backward_kernel,
:update_disp_veloc_kernel,
:update_potential_kernel,
:update_accel_elastic_kernel,
:update_veloc_elastic_kernel,
:update_accel_acoustic_kernel,
:update_veloc_acoustic_kernel,
:update_acoustic_lddrk_kernel,
:update_elastic_lddrk_kernel,
:compute_rho_kernel,
:compute_iso_kernel,
:compute_ani_kernel,
:compute_hess_kernel,
:compute_acoustic_kernel,
:compute_strength_noise_kernel,
:compute_ani_undoatt_kernel,
:compute_iso_undoatt_kernel,
:compute_strain_kernel,
:resort_array
]

big_kernels = [
:outer_core_impl_kernel_forward,
:outer_core_impl_kernel_adjoint,
:inner_core_impl_kernel_forward,
:inner_core_impl_kernel_adjoint,
:crust_mantle_impl_kernel_forward,
:crust_mantle_impl_kernel_adjoint
]

kernels = small_kernels + big_kernels

langs = $options[:langs]

BOAST::set_default_real_size(4)
BOAST::set_replace_constants(false)

kerns = kernels
kerns = kerns.select { |k,v| k.to_s.match($options[:kernel]) } if $options[:kernel]

# debug
#puts "kernels:"
#puts kerns

# output info
v = BOAST::get_boast_version()
puts ""
puts "BOAST version #{v}"
puts "-------------------------------"
puts "building kernel files:"
puts "-------------------------------"

# checks output directory
if File.exists? "#{$options[:output_dir]}/" then
  puts "directory exists: #{$options[:output_dir]}/"
else
  puts "directory does not exist: #{$options[:output_dir]}/"
  puts "please check your output-directory (--output-dir)"
  puts ""
  puts "exiting now..."
  abort "Error: output directory does not exist"
end

# loops over all kernels
kerns.each { |kern|
  puts kern.to_s
  # imports kernel ruby file
  require "./#{kern.to_s}.rb"
  langs.each { |lang|
    puts "  " + lang.to_s
    BOAST::set_lang( BOAST::const_get(lang))
    # outputs reference cuda kernel
    if $options[:display] && lang == :CUDA
      puts "  REF"
      k = BOAST::method(kern).call
      k.print
    end
    # generates kernels
    if lang == :CUDA then
      k = BOAST::method(kern).call(false)
      puts "  Generated"
      k.print if $options[:display]
      filename = "#{kern}.cu"
    elsif lang == :CL
      if big_kernels.include?(kern) then
        k = BOAST::method(kern).call(false, $options[:elem_per_thread])
      else
        k = BOAST::method(kern).call
      end
      puts "  Generated"
      k.print if $options[:display]
      filename = "#{kern}_cl.c"
    end

    # file name
    f = File::new("#{$options[:output_dir]}/#{filename}", "w+")

    # writes out specfem3d_globe info text at beginning of file
    v = BOAST::specfem3d_globe_header_info()

    # writes out generate kernel
    if lang == :CUDA then
      k_s = "#{v}" + k.to_s
      f.puts k_s
      if $options[:check] then
        puts "  building kernel"
        puts "  !! for checking with `make test_boast_kernels`, it might need BOAST version 2.0.2 !!"
        k.build(:LDFLAGS => " -L/usr/local/cuda-5.5.22/lib64", :NVCCFLAGS => "-arch sm_20 -O2 --compiler-options -Wall", :verbose => $options[:verbose] )
      end
    elsif lang == :CL then
      s = k.to_s
      res = "const char * #{kern}_program = \"\\\n"
      s.each_line { |line|
        res += line.sub("\n","\\n\\\n")
      }
      res += "\";\n"
      res = "#{v}\n" + res
      f.print res
      if $options[:check] then
        puts "  building kernel"
        k.build(:verbose => $options[:verbose], :platform_vendor => $options[:platform] )
      end
    end

    # regression testing
    if $options[:check] then
      puts "  testing kernel with ../kernels.test/ cases"
      inputs = k.load_ref_inputs("../kernels.test/")
      outputs_ref = k.load_ref_outputs("../kernels.test/")
      inputs.each_key { |key|
        puts key
        if big_kernels.include?(kern) then
          if lang == :CL then
            inputs[key].last[:local_work_size][0] /= $options[:elem_per_thread]
            inputs[key].last[:global_work_size][0] /= $options[:elem_per_thread]
          elsif lang == :CUDA then
            inputs[key].last[:block_size][0] /= $options[:elem_per_thread]
          end
        end
        puts k.run(*(inputs[key])).inspect
        puts k.compare_ref( outputs_ref[key], inputs[key] ).inspect
      }
    end
    f.close
  }
}

# output info
puts ""
puts "-------------------------------"
puts "building header & make files"
puts "-------------------------------"

# sorts kernels for file lists
kernels.sort!

langs.each { |lang|
  puts "  " + lang.to_s

  # opens output files
  if lang == :CUDA then
    suffix = ".cu"
    kern_proto_f = File::new("#{$options[:output_dir]}/kernel_proto.cu.h", "w+")
    kern_mk_f = File::new("#{$options[:output_dir]}/kernel_cuda.mk", "w+")
    kern_mk_f.puts "cuda_kernels_OBJS := \\"
  elsif lang == :CL
    suffix = "_cl.c"
    kern_inc_f = File::new("#{$options[:output_dir]}/kernel_inc"+suffix, "w+")
  end

  kern_list_f = File::new("#{$options[:output_dir]}/kernel_list.h", "w+")

  kernels.each { |kern|
    kern_list_f.puts "BOAST_KERNEL(#{kern.to_s});"

    if lang == :CUDA then
      require "./#{kern.to_s}.rb"
      BOAST::set_lang( BOAST::const_get(lang))
      k = BOAST::method(kern).call(false)
      BOAST::set_output( kern_proto_f )
      k.procedure.decl
      kern_mk_f.puts "\t$O/#{kern.to_s}.cuda-kernel.o \\"
    elsif lang == :CL
      kern_inc_f.puts "#include \"#{kern.to_s}#{suffix}\""
    end
  }

  elem_thread_check = "\n#if GPU_ELEM_PER_THREAD != #{$options[:elem_per_thread]}
#error \"Preprocessor macro mesh_constants_gpu.h::GPU-ELEM_PER_THREAD (GPU_ELEM_PER_THREAD) is different from BOAST's value (#{$options[:elem_per_thread]}).\"
#endif"

  if lang == :CUDA then
    kern_proto_f.puts elem_thread_check
  elsif lang == :CL
    kern_inc_f.puts elem_thread_check
  end


  # adds typedef-definitions for CUDA function pointers
  if lang == :CUDA then
    # array holding typedef-definitions for CUDA function pointers to **_forward and **_adjoint kernels
    proto_defs_cuda = []

    # loops over big kernels only
    big_kernels.each { |kern|
      require "./#{kern.to_s}.rb"
      BOAST::set_lang( BOAST::const_get(lang))
      k = BOAST::method(kern).call(false)

      # get procedure declaration as string
      my_typedef = ""
      pointer_name = ""
      proto = k.procedure.decl.to_s

      # typedefs for big kernels
      # kernel name string
      kernel_name = "#{kern.to_s}"

      # gets only function definition string, stripping leading __global__ .. declarations
      istart = proto.index(kernel_name)
      iend = proto.length
      function = proto[istart..iend]

      # gets arguments starting at "("
      istart = function.index("(")
      iend = function.length
      arguments = function[istart..iend]

      # makes pointer name string, e.g. outer_core_impl_kernel_forward -> (*outer_core_impl_kernel)
      istart = kernel_name.index("_impl_kernel")
      if istart then
        pointer_name = "(*" + kernel_name[0..istart] + "impl_kernel" + ")"
      else
        puts "pointer name not recognized in kernel: ",kernel_name
        abort "Error: kernel name invalid"
      end

      # declares typedef for function pointer, e.g. typedef void (*outer_core_impl_kernel) (...) ;
      my_typedef = "typedef void " + pointer_name + " " + arguments + " ;"

      # adds new typedef entries to list
      if my_typedef.length > 1 then
        #puts ""
        #puts "typedef: ",my_typedef
        #puts ""
        #puts "listed: ",proto_defs_cuda.any? { |s| s.include?(pointer_name)}
        #puts "index:",proto_defs.index(my_typedef)
        # checks if pointer already listed
        index = proto_defs_cuda.index { |s| s.include?(pointer_name) }
        if index then
          # already contained in list
          #puts "  #{kern.to_s} already listed: " + index.to_s
          # checks that arguments match
          if proto_defs_cuda[index] != my_typedef then
            puts "typedef-definition does not match for " + kernel_name
            puts "Please make sure that **_impl_kernel_forward and **_impl_kernel_adjoint have the same arguments"
            puts ""
            abort "Error: invalid arguments for forward and adjoint kernel call of " + kernel_name
          end
        else
          # adds new definition to list
          puts "  created typedef-definition for function pointer: " + pointer_name
          proto_defs_cuda.push(my_typedef)
        end
      end
    }

    # adds typedef-definitions
    if proto_defs_cuda.length > 0 then
      puts "  adding #{proto_defs_cuda.length.to_s} typedef-definitions for CUDA function pointers to file kernel_proto.cu.h"
      #puts proto_defs_cuda
      kern_proto_f.puts ""
      kern_proto_f.puts "// typedef-definitions for function pointers"
      kern_proto_f.puts ""
      # lists each definition
      proto_defs_cuda.each{ |mydef|
        kern_proto_f.puts mydef
        kern_proto_f.puts ""
      }
    end
  end

  # closing files
  kern_list_f.close
  if lang == :CUDA then
    kern_proto_f.close
    kern_mk_f.puts "\t$(EMPTY_MACRO)"
    kern_mk_f.close
  elsif lang == :CL
    kern_inc_f.close
  end
  puts "  Generated"
}

puts ""
puts "done"
puts ""
outdir = File.expand_path("#{$options[:output_dir]}/")
puts "see all files in directory: #{outdir}"
puts ""
