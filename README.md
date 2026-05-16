# ARGLy
Maximum Likelihood method for detecting ghost population with ARG.

## Build

Preferred simple build: GNU Make and C/C++ compilers.

```sh
make
```

This builds `build/argly` and `build/main`. The Makefile also builds bundled
`tskit` and `kastore` from `src/subprojects/tskit-1.0.0`; if that directory is
missing, it unpacks `src/subprojects/packagecache/tskit-1.0.0.tar.xz`, or
downloads the same archive as a fallback.

Useful options:

```sh
make CXX=g++ CC=gcc
make install PREFIX="$HOME/.local"
make clean
```

The old Meson build is still available if needed. Requirements: Meson, Ninja,
and a C/C++ compiler.

```sh
meson setup builddir src
meson compile -C builddir
```

The main binaries are created in `builddir/`:

```sh
./builddir/argly --help
```
