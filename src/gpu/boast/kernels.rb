require 'rubygems'
require 'BOAST'
require 'narray'
require './HEADER.rb'

def rndup( val, div)
  return (val%div) == 0 ? val : val + div - (val%div)
end

$options = {:output_dir => "./output"}

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
  opts.parse!
end

kernels = [
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
:noise_add_source_master_rec_kernel,
:noise_add_surface_movie_kernel,
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
:outer_core_impl_kernel_forward,
:outer_core_impl_kernel_adjoint,
:inner_core_impl_kernel_forward,
:inner_core_impl_kernel_adjoint,
:compute_rho_kernel,
:compute_iso_kernel,
:compute_ani_kernel,
:compute_hess_kernel,
:compute_acoustic_kernel,
:compute_strength_noise_kernel,
:crust_mantle_impl_kernel_forward,
:crust_mantle_impl_kernel_adjoint,
:compute_ani_undo_att_kernel,
:compute_iso_undo_att_kernel,
:compute_strain_kernel
]

langs = [ :CUDA, :CL]
BOAST::set_default_real_size(4)
BOAST::set_replace_constants(false)
class Float
  alias_method :old_to_s, :to_s
  def to_s
    s = ""+old_to_s
    s += "f" if BOAST::get_default_real_size == 4
    return s
  end
end

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
      k = BOAST::method(kern).call
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
      k = "#{v}" + k
      f.puts k
      if $options[:check] then
        puts "  building kernel"
        k.build( :LDFLAGS => " -L/usr/local/cuda-5.5.22/lib64", :NVCCFLAGS => "-arch sm_20 -O2 --compiler-options -Wall", :verbose => $options[:verbose] )
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
      proto = k.procedure.decl(false)[0..-3]+";"
      kern_proto_f.puts proto
      kern_mk_f.puts "\t$O/#{kern.to_s}.cuda-kernel.o \\"
    elsif lang == :CL
      kern_inc_f.puts "#include \"#{kern.to_s}#{suffix}\""
    end
  }

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
