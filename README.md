# ARGLy
Maximum Likelihood method for detecting ghost population with ARG.

## Build

Requirements: Meson, Ninja, and a C/C++ compiler.

```sh
meson setup builddir src
meson compile -C builddir
```

The main binaries are created in `builddir/`:

```sh
./builddir/argly --help
```
