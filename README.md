pySurf
======
[![Build Status](https://dev.azure.com/mdolab/Private/_apis/build/status/mdolab.pysurf?repoName=mdolab%2Fpysurf&branchName=main)](https://dev.azure.com/mdolab/Private/_build/latest?definitionId=36&repoName=mdolab%2Fpysurf&branchName=main)
[![Documentation Status](https://readthedocs.com/projects/mdolab-pysurf/badge/?version=latest&token=067843d7a8abdc6f145c3207abe46d4d73bf44ad406656e48b06a15a4cfa37a7)](https://mdolab-pysurf.readthedocs-hosted.com/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/mdolab/pysurf/branch/main/graph/badge.svg?token=4WJA4N7CDQ)](https://codecov.io/gh/mdolab/pysurf)

pySurf provides a wide range of geometric operations for triangulated surfaces and piecewise linear curves, including computing intersections between components.
These operations are differentiated using AD, which enables their use in gradient-based aerodynamic shape optimization.

**Note: The hyperbolic collar mesh generation functionality has been deprecated since v1.2.0.**

Citation
--------

Please cite the [pySurf journal article](https://arc.aiaa.org/doi/abs/10.2514/1.J056550) if you use any part of this code.

Ney R. Secco, John P. Jasa, Gaetan K. W. Kenway, and Joaquim R. R. A. Martins.  "Component-Based Geometry Manipulation for Aerodynamic Shape Optimization with Overset Meshes", AIAA Journal, Vol. 56, No. 9 (2018), pp. 3667-3679.

```
@article{Secco2018b,
	Author = {Ney R. Secco and John P. Jasa and Gaetan K. W. Kenway and Joaquim R. R. A. Martins},
	Doi = {10.2514/1.J056550},
	Journal = {AIAA Journal},
	Month = {September},
	Number = {9},
	Pages = {3667--3679},
	Title = {Component-based Geometry Manipulation for Aerodynamic Shape Optimization with Overset Meshes},
	Volume = {56},
	Year = {2018}}
```

License
-------
Copyright 2016-2019 MDO Lab

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
