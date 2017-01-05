# strong-skel
A research implementation of RS-S and RS-WS from Minden et al., "A recursive skeletonization factorization based on strong admissibility".

Built on the [FLAM library](https://github.com/klho/FLAM/) by Ken L. Ho, this research code implements the factorizations described in "A recursive skeletonization factorization based on strong admissibility" by Minden et al., available at [arXiv:1609.08130](https://arxiv.org/abs/1609.08130).

The main functions provided are as follows:
- srskelf.m: the strong recursive skeletonization factorization (aka "RS-S") for symmetric matrices.
- srskelf_asym.m: RS-S for asymmetric matrices.
- srskelf_hybrid.m: the hybrid recursive skeletonization factorization (aka "RS-WS") for symmetric matrices.
- srskelf_asymhybrid.m: RS-WS for asymmetric matrices.

Tests are available as MATLAB functions in the directory "test" and may be run with no arguments, assuming that FLAM and the strong-skel directories are in your MATLAB path.  The tests provided are as follows:
- ie_square1_s.m: The setup of Example 1 from Minden et al.
- ie_cube1_s.m: The setup of Example 2 from Minden et al.
- ie_sphere_s.m: The setup of Example 3 from Minden et al.


This code is provided for reference under the GNU GPL.
