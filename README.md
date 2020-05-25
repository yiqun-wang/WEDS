# WEDS

This code implements Wavelet Energy Decomposition Signature (WEDS) in our SIGGRAPH 2020 paper:

["MGCN: Descriptor Learning using Multiscale GCNs"](https://arxiv.org/abs/2001.10472) 

by Yiqun Wang, Jing Ren, Dong-Ming Yan, Jianwei Guo, Xiaopeng Zhang, Peter Wonka.

Please consider citing the above paper if this code/program (or part of it) benefits your project. 


## Usage

Please run "demo.m" using MATLAB for generating Wavelet Energy Decomposition Signature (WEDS) on Linux, Mac or Windows system.

1. Input: "off", "obj", or "ply" format of mesh model in "data_mesh" folder.

2. Output: WEDS on every vertex and wavelets are saved in "data_wavelet" folder. 
	
## Cite

    @article{wang2020mgcn,
      title=      {MGCN: Descriptor Learning using Multiscale GCNs},
      author=     {Wang, Yiqun and Ren, Jing and Yan, Dong-Ming and Guo, Jianwei and Zhang, Xiaopeng and Wonka, Peter},
      journal=    {ACM  Trans. on Graphics (Proc. {SIGGRAPH})},
      year=       {2020},
    }

## License

This program is free software; you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation; either version 2 of 
the License, or (at your option) any later version. 
