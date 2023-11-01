#CC=clang++
CC=g++
CXXFLAGS=-lstdc++fs -std=c++17 -O3 -W
INCLUDES=-Iinclude
LINKERS=-L/usr/local/lib -L/usr/lib -lstdc++fs -lcrypto -lssl
WARNINGS=-Wno-unused-result -Wno-sign-compare -Wno-maybe-uninitialized -Wno-narrowing -Wno-return-type -Wno-format

SUPPORT_SRCFILES = src/Util.cpp \
                   src/EncryptionEngine.cpp \
                   src/FileSimulator.cpp \
                   src/PathORAM.cpp \
                   src/ObliviousHeap.cpp \
                   #TODO:
                   #src/MongoConnector.cpp

HEADER_ONLY_SOURCES = include/PrivIndexJoin.h \
                      include/ObliIndexJoin.h \
                      include/OneORAMIndexJoin.h \
                      include/ObliMultiwayJoin.h \
                      include/OneORAMMultiwayJoin.h \
                      include/ObliviousBTree.h \
                      include/OneORAMBTree.h \
                      include/ObliMultiwayBTree.h \
                      include/OneORAMMultiwayBTree.h \
                      include/ObliDatabaseJoin.h \
                      include/ObliviousSort.h \
                      include/App.h \
                      include/Basic.h \
                      include/Schema.h \
                      include/Config.h

SUPPORT_OBJFILES = $(patsubst src/%.cpp, obj/%.o, $(SUPPORT_SRCFILES))

EXECUTABLE = odbj_load \
             odbj_equi_join \
             odbj_band_join \
             oij_type1_load \
             oij_type1_inlj_equi_join \
             oij_type1_smj_equi_join \
             oij_type1_inlj_band_join \
             oij_type1_smj_band_join \
             oij_type1_multi_load \
             oij_type1_multi_join \
             one_type1_load \
             one_type1_inlj_equi_join \
             one_type1_smj_equi_join \
             one_type1_inlj_band_join \
             one_type1_smj_band_join \
             one_type1_multi_load \
             one_type1_multi_join \

all: $(EXECUTABLE)

%_load: %_load.o $(SUPPORT_OBJFILES)
	$(CC) $(CXXFLAGS) $(INCLUDES) $^ -o $@ $(LINKERS) $(WARNINGS)

%_join: %_join.o $(SUPPORT_OBJFILES)
	$(CC) $(CXXFLAGS) $(INCLUDES) $^ -o $@ $(LINKERS) $(WARNINGS)

%_padding: %_padding.o $(SUPPORT_OBJFILES)
	$(CC) $(CXXFLAGS) $(INCLUDES) $^ -o $@ $(LINKERS) $(WARNINGS)

%.o: %.cpp $(HEADER_ONLY_SOURCES)
	$(CC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ $(WARNINGS)

obj/%.o: src/%.cpp include/%.h
	@mkdir -p $(@D)
	$(CC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ $(WARNINGS)

clean:
	rm -rf obj $(EXECUTABLE) *.o
