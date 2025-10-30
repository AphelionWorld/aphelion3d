# Aphelion3D – CMake + MinGW + VS Code starter

## Prereqs (Windows/MinGW)
- Install MSYS2: https://www.msys2.org/
- From MSYS2 MinGW64 shell: `pacman -S --needed mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake mingw-w64-x86_64-ninja`
- Add `C:\msys64\mingw64\bin` to your PATH (or use the MinGW64 terminal inside VS Code).

## VS Code
- Install extensions: *C/C++* and *CMake Tools*.
- Open this folder, select preset **MinGW Debug (Ninja)** in the status bar if not auto-selected.
- Click **Build**, then **Debug** the `hello` or `viewer` targets.

## CLI
```bash
cmake --preset mingw-debug
cmake --build --preset mingw-debug -j
ctest --preset mingw-debug
cmake --install build/mingw-release --prefix out/install
```

## Structure
- `src/libutil` – header-only interface lib
- `src/libmath` – compiled library (static by default; toggle BUILD_SHARED_LIBS=ON for DLL)
- `src/app/hello` – console app using the libs
- `src/app/viewer` – second app using both libs
- `tests` – GoogleTest unit tests

## Packaging
- Installs headers to `include/aphelion3d/`, libs to `lib/`, apps to `bin/`
- Exports targets as `Aphelion3D::libmath` and `Aphelion3D::libutil`
- CMake package config installed to `lib/cmake/Aphelion3D`
