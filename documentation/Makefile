LIBAIMS=${aims}/lib/libaims.190522.scalapack.mpi.so
include_dir=${anaconda}/anaconda3/envs/genarris_env/include

aims_w.so:
	aims_w.f90 f2py --f90exec=mpiifort --fcompiler=intelem -m aims_w \ -c aims_w.f90 ${LIBAIMS} -I${include_dir}
clean:
	rm aims_w.*.so
