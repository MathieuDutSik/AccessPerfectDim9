PROGRAM		= NC_AccessForms

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++17
#CC		= clang++ -g -Weverything -std=c++17 -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null
#CC		= clang++ -g -Weverything -std=c++17
#CC		= clang++ -g -Weverything -std=c++17 -fsanitize=bounds -fsanitize-memory-use-after-dtor -fsanitize-memory-track-origins
#CC		= clang++ -O2 -Weverything -std=c++17 -DDEBUG
#CC              = /usr/share/clang/scan-build-3.8/libexec/c++-analyzer -std=c++17 -O2
#CC              = /home/mathieu/opt/gcc-6.3.0/bin/g++ -std=c++17 -O3
#CC		= g++ -std=c++17 -Wall -Wextra -O2 -g -ftemplate-backtrace-limit=0 -fsanitize=float-divide-by-zero -fsignaling-nans
CC		= g++ -std=c++17 -Wall -Wextra -O3 -g
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -ggdb
#CC		= g++ -std=c++17 -Wall -Wextra -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment
#CC		= g++ -std=c++17 -Wall -Wextra -g -O2 -static-libgcc -static-libstdc++

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../polyhedral_common/src_basic -I../polyhedral_common/basic_common_cpp   -I../polyhedral_common/src_number -I../polyhedral_common/src_matrix -I../polyhedral_common/src_comb -I../polyhedral_common/src_graph -I../polyhedral_common/src_group -I../polyhedral_common/src_latt -I../polyhedral_common/src_dualdesc -I../polyhedral_common/src_poly -I../polyhedral_common/src_perfect -I../polyhedral_common/permutalib/src -I../polyhedral_common/src_lorentzian -I../polyhedral_common/src_isotropy -I.
CFLAGS 		= -I$(GMP_INCDIR) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(MPI_INCLUDE_CPP) $(CDDLIB_INCLUDE) $(NAUTY_INCLUDE) $(NETCDF_CXX_ALLINC)
CHASHMAP        = -I../polyhedral_common/basic_common_cpp/sparse-map/include/tsl/ -I../polyhedral_common/basic_common_cpp/robin-map/include/tsl/ -I../polyhedral_common/basic_common_cpp/hopscotch-map/include/tsl/
CDEFS           = -DUSE_CDDLIB -DERR_LOG
CTOTAL          = $(CLOCAL) $(CFLAGS) $(CHASHMAP) $(CDEFS)

LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(MPI_LINK_CPP) $(CDDLIB_DOUBLE_LINK) $(NAUTY_LINK) $(NETCDF_CXX_LINK)


%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CTOTAL) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
