# Velvet-Noise Convolution

Efficient implementation of a sparse pseudo-random filter.

Includes five sequential versions that optimize the code to take advantage of the fast caches of the processors and reduce both, computational and spatial cost.

Uses OpenMP to parallelize the processing of different channels.

More information in:

Jose A. Belloch, Jose M. Badia, German Leon and Vesa Välimäki, "Efficient Velvet-Noise Convolution in Multicore Processors".

submitted to the Journal of the Audio Engineering Society on January 2024.
