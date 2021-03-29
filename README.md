# CGM

CGM (Computer Graphics Maths) is a C99 maths library for 2D and 3D graphics applications. The idea is to provide a nice minimal implementation for mathematical functions so users can build their own systems for physics simulation and game logic.

# Features
- Linear Algebra types and functions for 2D & 3D
	- Vectors (CgmVec2/CgmVec3/CgmVec4)
	- Matrices (CgmMat3x2/CgmMat4x4)
- 2D & 3D Collision Checking that return both contact normal & distance.
- Float16 encode and decode for packing data into cache lines.
- Misc floating point helper functions such as lerp, sign, remap

# Usage

Everywhere you need to use this library, put this at the top of the file:
```
#include "cgm.h"
```
In a single compilation unit, include the source file.
```
#include "cgm.c"
```

