# SiteAmplification
This work is to calculate site amplification given a 1D model. 
Vertically incident body wave site amplification is computed using propagation matrix method (Aki & Richards, 2002).
Horizontally propagating surface wave site amplification is including body and surface waves (Bowden & Tsai, 2017).
The surface wave integrals are calcualated using mat_disperse function from the MATLAB package forked from [Yiran Ma ](https://github.com/yiran06/mat_disperse)

benchmark.m should produce an exact same figure as Figure 1C in Bowden and Tsai (2017):
<img width="390" alt="Screenshot 2023-02-02 at 14 39 00" src="https://user-images.githubusercontent.com/20884599/216465771-0428c777-f590-441c-aea3-0eb56e3c73ce.png">
