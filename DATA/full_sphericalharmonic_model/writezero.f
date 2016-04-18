        dimension a(21,441),b(441)
        open(20,file='zero_sh')
        write(20,*) 20
        write(20,*) (b(i),i=1,441)
        open(21,file='zero_3d')
        write(21,*) 20,21
        do i=1,21
          write(21,*) (a(i,j),j=1,441)
        enddo
        end