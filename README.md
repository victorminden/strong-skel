# strong-skel
A research implementation of RS-S and RS-WS from Minden et al., "A recursive skeletonization factorization based on strong admissibility" ([Journal version](http://dx.doi.org/10.1137/16M1095949), [arXiv](https://arxiv.org/abs/1609.08130))

We emphasize that this code is provided under the GNU GPL as an as-is research reference implementation.  In particular, it is not an actively maintained library.  

Contributors:
Victor Minden, Ken L. Ho, Anil Damle, Lexing Ying

## Details
Built on the [FLAM library](https://github.com/klho/FLAM/) by Ken L. Ho, this research code implements the factorizations described in "A recursive skeletonization factorization based on strong admissibility" by Minden et al.

The main functions provided are as follows:
- srskelf.m: the strong recursive skeletonization factorization (aka "RS-S") for symmetric matrices.  This function uses Cholesky factorizations of diagonal subblocks that can be assumed positive-definite in exact arithmetic.
- srskelf_asym.m: RS-S for asymmetric matrices.  This function uses LU factorizations where srskelf.m uses Cholesky.
- srskelf_hybrid.m: the hybrid recursive skeletonization factorization (aka "RS-WS") for symmetric matrices.  This function is like srskelf.m, but additionally interlaces levels of standard ("weak") skeletonization.
- srskelf_asymhybrid.m: RS-WS for asymmetric matrices.  This function is like srskelf_asym.m, but additionally interlaces levels of standard ("weak") skeletonization.

Tests are available as MATLAB functions in the directory "test" and may be run with no arguments, assuming that FLAM and the strong-skel directories are in your MATLAB path.  The tests provided are as follows:
- ie_square.m: The setup of Example 1 as described in Minden et al.
- ie_cube.m: The setup of Example 2 as described in Minden et al.
- ie_sphere.m: The setup of Example 3 as described in Minden et al.