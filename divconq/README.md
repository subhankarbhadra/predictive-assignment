# Divide and conquer algorithms for community detection

We provide MATLAB implementations of the algorithms PACE and GALE [1]. The script `local_patching.m` contains a number of example simulations using spectral clustering.

The script `patching_large_sim.m` is used to run our large scale simulations, where we compare global spectral clustering against GALE (which also uses spectral clustering on sampled patches).

#### How to run?
Clone the repository:

```shell
git clone https://gitlab.com/soumendu041/divide-and-conquer-community-detection.git
```

To run the example (small scale) simulations, run the following in MATLAB:

```matlab
> cd <path-to-cloned-repository>
> local_patching
```
To run a large scale simulation, run the following in MATLAB (**warning:** the example graph in `patching_large_sim.m` has `n = 10^7` nodes, and `npar = 100` MATLAB workers are used by default for parallel computation; change these parameters according to your machine specs):

```matlab
> cd <path-to-cloned-repository>
> patching_large_sim
```

### References
1. Soumendu Sundar Mukherjee, Purnamrita Sarkar, and Peter J. Bickel. [Two provably consistent divide and conquer clustering algorithms for large networks](https://arxiv.org/abs/1708.05573). arXiv:1708.05573, 2017.
