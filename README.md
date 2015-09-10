# CIT

The CIT (Category to Instance Tracker) is for tracking objects identified as targets by a pre-trained category detector.  

Created by David Hall and Pietro Perona at Caltech

### Introduction
CIT is a method for online, real-time tracking of objects. Tracking is treated as a repeated detection problem where potential target objects are identified with a pre-trained category detector
and object identity across frames is established by individual-specific detectors.
The individual detectors are (re-)trained online from a single positive example whenever there is a coincident category detection. This ensures that the tracker is robust to drift.

For all of technical details please refer to the following publications:
* [Online, Real-Time Tracking Using a Category-to-Individual Detector](http://www.vision.caltech.edu/~dhall/projects/CIT/Data/ECCV2014_HALL.pdf)
* [From Categories to Indviduals in Real Time---A Unified Boosting Approach](http://www.vision.caltech.edu/~dhall/projects/CategoriesToIndividuals/Data/CVPR2014_HALL.pdf)

### Citation
Please cite our publications if this code is used in any work.
```
Online, Real-Time Tracking Using a Category-to-Individual Detector. 
D.Hall and P.Perona. 
ECCV 2014, Zurich, Switzerland.
```
```
From Categories to Individual in Real Time - A Unified Boosting Approach.
D.Hall and P.Perona. 
CVPR 2014, Columbus, USA.
```

### License
CIT is licensed under the Simplified BSD License [see LICENSE for details]

### Prerequisites:
* MATLAB 2014b or greater
* [Piotr's Image & Video Matlab Toolbox](https://github.com/pdollar/toolbox)

### Installation
* Download the source code
* Run compile.m
* Download the [FPOQ dataset](http://www.vision.caltech.edu/Image_Datasets/FiftyPeopleOneQuestion/FPOQ.zip)
* Extract into the CIT/data/Datasets directory. 

### Usage
    
Run the script example.m. It gives a simple example of how CIT can be used to track
faces in a video. 

A pre-trained face detector for the FPOQ dataset is provided along with 
corresponding spread values computed on the Berlinale video.
You can recompute these values by setting ESTIMATE_SPREAD = true in example.m

To train your own category detector refer to detector/train_category_detector.m
