# CIT
Category to Instance Tracker

Copyright 2014 D.Hall and P.Perona.
[dhall-at-caltech-dot-edu]
Please email me if you find bugs, or have suggestions or questions!
Licensed under the Simplified BSD License [see bsd.txt]
********************************************************************************
External Dependencies:
MATLAB 2014b or greater
Piotr's Image & Video Matlab Toolbox(https://github.com/pdollar/toolbox)

The script example.m  gives a simple example of how CIT can be used to track
faces in a video. 

To run this script:
Download the FPOQ dataset found here:
http://www.vision.caltech.edu/Image_Datasets/FiftyPeopleOneQuestion/FPOQ.zip
Extract into the data/Datasets directory. 

Run compile.m to build the detector/acf_detect_all.cpp

A pre-trained face detector for the FPOQ dataset is provided along with 
corresponding spread values computed on the Berlinale video. Refer to [2] for details.
You can recompute these values by setting ESTIMATE_SPREAD = true in example.m

********************************************************************************
Please cite our publications if this code is used in any work.

[1] "Online, Real-Time Tracking Using a Category-to-Individual Detector". 
 D.Hall and P.Perona. ECCV 2014, Zurich, Switzerland.

[2] "From Categories to Individual in Real Time - A Unified Boosting Approach". 
 D.Hall and P.Perona. CVPR 2014, Columbus, USA.
********************************************************************************
