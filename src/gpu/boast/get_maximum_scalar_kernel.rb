module BOAST
  def BOAST::get_maximum_scalar_kernel(ref = true, block_size_transfer = 256)
    BOAST::get_maximum_kernel(:scalar, ref, block_size_transfer)
  end
  def BOAST::get_maximum_kernel(type, ref = true, block_size_transfer = 256)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    size = Int("size", :dir => :in)

    if type == :scalar then
      function_name = "get_maximum_scalar_kernel"
      array = Real("array", :dir => :in, :dim => [ Dim(size)])
    elsif type == :vector then
      function_name = "get_maximum_vector_kernel"
      array = Real("array", :dir => :in, :dim => [ Dim(size*3)] )
    else
      raise "Unsupported maximum type: #{type}!"
    end
    
    d_max = Real("d_max", :dir => :out, :dim => [ Dim()])
    blocksize_transfer = Int("BLOCKSIZE_TRANSFER", :const => block_size_transfer)
    p = Procedure(function_name, [array, size, d_max])
    if(get_lang == CUDA and ref) then
      @@output.print File::read("references/#{function_name}.cu")
    elsif(get_lang == CUDA or get_lang == CL) then
      make_specfem3d_header( :blocksize_transfer => block_size_transfer )
      decl p
      sdata = Real("sdata",  :local => true, :dim => [Dim(blocksize_transfer)] )
      tid =   Int("tid")
      bx =    Int("bx")
      i =     Int("i")
      s =     Int("s")
      decl sdata
      decl tid
      decl bx
      decl i
      decl s
      print tid === get_local_id(0)
      print bx === get_group_id(1)*get_num_groups(0) + get_group_id(0)
      print i === tid + bx*get_local_size(0)
      if type == :scalar then
        print sdata[tid] === Ternary( i < size, fabs(array[i]), 0.0)
      else
        print sdata[tid] === Ternary( i < size, sqrt(array[i*3+0]*array[i*3+0]+array[i*3+1]*array[i*3+1]+array[i*3+2]*array[i*3+2]), 0.0)
      end
      print barrier(:local)
      print s === get_local_size(0)/2
      print While(s > 0) {
        print If(tid < s) {
          print If( sdata[tid] < sdata[tid + s] ) {
            print sdata[tid] === sdata[tid + s]
          }
        }
        print s === Expression(">>",s,1)
        print barrier(:local)
      }
      print If(tid == 0) {
        print d_max[bx] === sdata[0]
      }
      close p
    else
      raise "Unsupported language!"
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end
end
