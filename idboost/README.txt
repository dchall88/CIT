IDBoost Code
Copyright 2014 D.Hall and P.Perona.
[dhall-at-caltech-dot-com]
Please email me if you find bugs, or have suggestions or questions!
Licensed under the Simplified BSD License [see bsd.txt]
********************************************************************************
External Dependencies:
Piotr's Image & Video Matlab Toolbox (http://vision.ucsd.edu/~pdollar/toolbox/doc/)
VideoTools (http://www.vision.caltech.edu/~dhall/projects/CategoriesToIndividuals/)

The functions create_individual_clf.m and compute_feature_spread.m  contains the IDBoost code as 
outlined in our work.

The function acfDetect is a modified version of the function found in Dollar's toolbox. This version 
supports video and returns the detected feature vectors.

The script example.m  gives a simple example of how a category-to-individual detector can be used to 
find the frames in a video containing a particular individual. To run this script
download our datasets and place them in Data/Datasets. A pre-trained face detector for the FPOQ 
dataset is provided. Also ensure the Code/ directory is added to the path with subfolders 

********************************************************************************
Please cite our publication if this code is used in any work.
"From Categories to Individual in Real Time - A Unified Boosting Approach". 
 D.Hall and P.Perona. CVPR 2014, Columbus, USA.
********************************************************************************