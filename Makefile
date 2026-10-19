CC ?= cc
CXX ?= c++
AR ?= ar
RM ?= rm -f
MKDIR_P ?= mkdir -p
CURL ?= curl
TAR ?= tar

BUILD_DIR ?= build
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

SRC_DIR := src
SUBPROJECTS_DIR := $(SRC_DIR)/subprojects
TREELLH_DIR := $(SUBPROJECTS_DIR)/treellh
TSKIT_VERSION := 1.0.0
TSKIT_DIR := $(SUBPROJECTS_DIR)/tskit-$(TSKIT_VERSION)
TSKIT_ARCHIVE := $(SUBPROJECTS_DIR)/packagecache/tskit-$(TSKIT_VERSION).tar.xz
TSKIT_URL := https://github.com/tskit-dev/tskit/releases/download/C_$(TSKIT_VERSION)/tskit-$(TSKIT_VERSION).tar.xz
TSKIT_SHA256 := f8a085fc6ef170706b76a307e18be4ab425f31a1f50fa9765ab2c2d81eef65c2
KASTORE_DIR := $(TSKIT_DIR)/subprojects/kastore

OBJ_DIR := $(BUILD_DIR)/obj
LIB_DIR := $(BUILD_DIR)/lib
TSKIT_STAMP := $(BUILD_DIR)/.tskit-unpacked

CPPFLAGS ?=
CFLAGS ?= -O3 -DNDEBUG
CXXFLAGS ?= -O3 -DNDEBUG
CXXSTD ?= -std=c++20
LDFLAGS ?=
LDLIBS ?= -lm
NLOPT_PREFIX ?= $(shell brew --prefix nlopt 2>/dev/null)
NLOPT_CFLAGS ?= $(if $(NLOPT_PREFIX),-I$(NLOPT_PREFIX)/include,$(shell pkg-config --cflags nlopt 2>/dev/null))
NLOPT_LIBS ?= $(if $(NLOPT_PREFIX),-L$(NLOPT_PREFIX)/lib -lnlopt,$(shell pkg-config --libs nlopt 2>/dev/null || printf '%s' '-lnlopt'))

COMMON_CPPFLAGS := \
	-I$(TREELLH_DIR) \
	-I$(TSKIT_DIR) \
	-I$(KASTORE_DIR)

TSKIT_C_SOURCES := \
	$(TSKIT_DIR)/tskit/core.c \
	$(TSKIT_DIR)/tskit/tables.c \
	$(TSKIT_DIR)/tskit/trees.c \
	$(TSKIT_DIR)/tskit/genotypes.c \
	$(TSKIT_DIR)/tskit/stats.c \
	$(TSKIT_DIR)/tskit/convert.c \
	$(TSKIT_DIR)/tskit/haplotype_matching.c

KASTORE_C_SOURCE := $(KASTORE_DIR)/kastore.c
TREELLH_CPP_SOURCE := $(TREELLH_DIR)/treellh.cpp

TSKIT_OBJECTS := $(patsubst $(TSKIT_DIR)/tskit/%.c,$(OBJ_DIR)/tskit/%.o,$(TSKIT_C_SOURCES))
KASTORE_OBJECT := $(OBJ_DIR)/kastore/kastore.o
TREELLH_OBJECT := $(OBJ_DIR)/treellh/treellh.o $(OBJ_DIR)/treellh/likelihood_cache.o

TSKIT_LIB := $(LIB_DIR)/libtskit.a
KASTORE_LIB := $(LIB_DIR)/libkastore.a
TREELLH_LIB := $(LIB_DIR)/libtreellh.a

BINARIES := $(BUILD_DIR)/argly $(BUILD_DIR)/main $(BUILD_DIR)/argly-epoch-fit

.PHONY: all clean distclean install tskit

all: $(BINARIES)

$(BUILD_DIR)/argly-epoch-fit: $(SRC_DIR)/epoch_cli.cpp $(SRC_DIR)/epoch_model.cpp $(SRC_DIR)/epoch_fit.cpp $(SRC_DIR)/epoch_model.h $(SRC_DIR)/epoch_fit.h $(TREELLH_LIB) $(TSKIT_LIB) $(KASTORE_LIB)
	$(CXX) $(CPPFLAGS) $(COMMON_CPPFLAGS) $(NLOPT_CFLAGS) $(CXXFLAGS) $(CXXSTD) $(LDFLAGS) -o $@ $(SRC_DIR)/epoch_cli.cpp $(SRC_DIR)/epoch_model.cpp $(SRC_DIR)/epoch_fit.cpp $(TREELLH_LIB) $(TSKIT_LIB) $(KASTORE_LIB) $(NLOPT_LIBS) $(LDLIBS)

$(BUILD_DIR)/epoch_test: $(SRC_DIR)/epoch_test.cpp $(SRC_DIR)/epoch_model.cpp $(SRC_DIR)/epoch_fit.cpp $(SRC_DIR)/epoch_model.h $(SRC_DIR)/epoch_fit.h $(TREELLH_LIB) $(TSKIT_LIB) $(KASTORE_LIB)
	$(CXX) $(CPPFLAGS) $(COMMON_CPPFLAGS) $(NLOPT_CFLAGS) $(CXXFLAGS) $(CXXSTD) $(LDFLAGS) -o $@ $(SRC_DIR)/epoch_test.cpp $(SRC_DIR)/epoch_model.cpp $(SRC_DIR)/epoch_fit.cpp $(TREELLH_LIB) $(TSKIT_LIB) $(KASTORE_LIB) $(NLOPT_LIBS) $(LDLIBS)

.PHONY: test-epoch
test-epoch: $(BUILD_DIR)/epoch_test
	$(BUILD_DIR)/epoch_test

.PHONY: test-cache
test-cache: $(BUILD_DIR)/cache_test
	$(BUILD_DIR)/cache_test

$(BUILD_DIR)/cache_test: $(SRC_DIR)/cache_test.cpp $(TREELLH_DIR)/likelihood_cache.h $(TREELLH_LIB) $(TSKIT_LIB) $(KASTORE_LIB)
	$(CXX) $(CPPFLAGS) $(COMMON_CPPFLAGS) $(CXXFLAGS) $(CXXSTD) $(LDFLAGS) -o $@ $(SRC_DIR)/cache_test.cpp $(TREELLH_LIB) $(TSKIT_LIB) $(KASTORE_LIB) $(LDLIBS)

tskit: $(TSKIT_STAMP)

$(BUILD_DIR)/argly: $(OBJ_DIR)/argly_cli.o $(OBJ_DIR)/fit_search.o $(TREELLH_LIB) $(TSKIT_LIB) $(KASTORE_LIB)
	$(CXX) $(LDFLAGS) -o $@ $^ $(NLOPT_LIBS) $(LDLIBS)

$(BUILD_DIR)/main: $(OBJ_DIR)/main.o $(TREELLH_LIB) $(TSKIT_LIB) $(KASTORE_LIB)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

$(TREELLH_LIB): $(TREELLH_OBJECT) | $(LIB_DIR)
	$(AR) rcs $@ $^

$(TSKIT_LIB): $(TSKIT_OBJECTS) | $(LIB_DIR)
	$(AR) rcs $@ $^

$(KASTORE_LIB): $(KASTORE_OBJECT) | $(LIB_DIR)
	$(AR) rcs $@ $^

$(OBJ_DIR)/argly_cli.o: $(TSKIT_STAMP) $(SRC_DIR)/argly_cli.cpp $(SRC_DIR)/fit_search.h $(TREELLH_DIR)/treellh.h $(TREELLH_DIR)/likelihood_cache.h | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(COMMON_CPPFLAGS) $(CXXFLAGS) $(CXXSTD) -c $(SRC_DIR)/argly_cli.cpp -o $@

$(OBJ_DIR)/fit_search.o: $(SRC_DIR)/fit_search.cpp $(SRC_DIR)/fit_search.h $(TREELLH_DIR)/likelihood_cache.h | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(COMMON_CPPFLAGS) $(NLOPT_CFLAGS) $(CXXFLAGS) $(CXXSTD) -c $(SRC_DIR)/fit_search.cpp -o $@

$(OBJ_DIR)/main.o: $(TSKIT_STAMP) $(SRC_DIR)/main.cpp $(TREELLH_DIR)/treellh.h | $(OBJ_DIR)
	$(CXX) $(CPPFLAGS) $(COMMON_CPPFLAGS) $(CXXFLAGS) $(CXXSTD) -c $(SRC_DIR)/main.cpp -o $@

$(OBJ_DIR)/treellh/treellh.o: $(TSKIT_STAMP) $(TREELLH_CPP_SOURCE) $(TREELLH_DIR)/treellh.h | $(OBJ_DIR)/treellh
	$(CXX) $(CPPFLAGS) $(COMMON_CPPFLAGS) $(CXXFLAGS) $(CXXSTD) -c $(TREELLH_CPP_SOURCE) -o $@

$(OBJ_DIR)/treellh/likelihood_cache.o: $(TSKIT_STAMP) $(TREELLH_DIR)/likelihood_cache.cpp $(TREELLH_DIR)/likelihood_cache.h $(TREELLH_DIR)/treellh.h | $(OBJ_DIR)/treellh
	$(CXX) $(CPPFLAGS) $(COMMON_CPPFLAGS) $(CXXFLAGS) $(CXXSTD) -c $(TREELLH_DIR)/likelihood_cache.cpp -o $@

$(OBJ_DIR)/tskit/%.o: $(TSKIT_STAMP) $(TSKIT_DIR)/tskit/%.c | $(OBJ_DIR)/tskit
	$(CC) $(CPPFLAGS) $(COMMON_CPPFLAGS) $(CFLAGS) -std=c99 -c $(TSKIT_DIR)/tskit/$*.c -o $@

$(OBJ_DIR)/kastore/kastore.o: $(TSKIT_STAMP) $(KASTORE_C_SOURCE) | $(OBJ_DIR)/kastore
	$(CC) $(CPPFLAGS) $(COMMON_CPPFLAGS) $(CFLAGS) -std=c99 -c $(KASTORE_C_SOURCE) -o $@

$(TSKIT_STAMP): $(TSKIT_ARCHIVE) | $(BUILD_DIR)
	if [ ! -f "$(TSKIT_DIR)/tskit.h" ]; then $(TAR) -C $(SUBPROJECTS_DIR) -xf $(TSKIT_ARCHIVE); fi
	touch $@

$(TSKIT_ARCHIVE):
	$(MKDIR_P) $(dir $@)
	$(CURL) -L -o $@ $(TSKIT_URL)
	if command -v shasum >/dev/null 2>&1; then \
		printf '%s  %s\n' '$(TSKIT_SHA256)' '$@' | shasum -a 256 -c -; \
	else \
		printf '%s  %s\n' '$(TSKIT_SHA256)' '$@' | sha256sum -c -; \
	fi

$(BUILD_DIR) $(OBJ_DIR) $(OBJ_DIR)/treellh $(OBJ_DIR)/tskit $(OBJ_DIR)/kastore $(LIB_DIR):
	$(MKDIR_P) $@

install: all
	$(MKDIR_P) $(DESTDIR)$(BINDIR)
	cp $(BINARIES) $(DESTDIR)$(BINDIR)/

clean:
	$(RM) -r $(BUILD_DIR)

distclean: clean
	$(RM) -r $(TSKIT_DIR)
