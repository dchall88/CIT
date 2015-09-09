%% A simple example script demonstrating the use of a 
% category-to-instance tracker (CIT). 
% 
% CIT is an application of the IDBoost Algorithm 
%
% Piotr's Image & Video Matlab Toolbox is required! It can be found here
% https://github.com/pdollar/toolbox
%
% A category detector is trained using ACF 
%   (for details refer to acfDemoCal in Piotr's toolbox.)   
% The IDBoost parameters are estimated using a validation set
%   (for details refer to create_individual_clf.m in the idboost directory) 
% Tracking is performed on 500 frames of a video

%% Parameters and Settings
addpath(genpath('utils'), genpath('idboost'), genpath('detector'), genpath('cit'))

TRAIN_CATEGORY_DETECTOR = false;
categoryDetectorName = 'data/Detector/acfFPOQDetector.mat'; % Pretrained face detector on the FPOQ dataset

ESTIMATE_SPREAD = false;
idBoostSpreadName = 'data/Detector/acfFPOQidBoostSpread.mat'; % Pretrained face detector on the FPOQ dataset

valVideoName = 'data/Datasets/FPOQ/Videos/Berlinale-uv5l-uu6SIw.mp4';
valGtName = 'data/Datasets/FPOQ/Annotations/Berlinale.vbb';
testVideoName = 'data/Datasets/FPOQ/Videos/Alabama-csHddXn91YE.mp4';

% Video Parameters
VideoParameters = struct;
VideoParameters.frames = 1001:1500;

%% Tuning parameter for IDBoost
% This is the most interesting parameter - the smaller it is the faster the
% instance detector is evaluated. The disadvantage is that the
% detector now has a narrower model of the individual so it is more likely to fail. 
TrackingParameters.idBoostBeta = 1;

% Other Parameters
TrackingParameters.destroyAfterNFrames = 5;
TrackingParameters.minCategoryDetectorScore = 30;
TrackingParameters.minDetectionOverlap = 0.4;
TrackingParameters.VIEW_RESULTS = false;

%% Load/Train Category Detector
if TRAIN_CATEGORY_DETECTOR
    % Refer to train_category_detector and acfDemoCal.m for details on how
    % to train your own category_detector.
    disp('Train Category Detector')
    CategoryDetector = train_category_detector(imgDir, gtDir, resDir, modelName);
    % In practice - after this step - you would also calibrate the
    % category detector on a validation set to operate at some desired 
    % value of precision/recall 
else
    disp('Load Category Detector')
    CategoryDetector = load(categoryDetectorName); 
    CategoryDetector = CategoryDetector.detector;
    CategoryDetector.opts.returnPatch = 1; 
end

%% Estimate sigma (the spread) from a validation set. For the IDBoost Algorithm.
if ESTIMATE_SPREAD
     disp('Estimate Spread')
    [features, thumbnails] = ...
        extract_feature_tracks(valVideoName, valGtName, CategoryDetector);

    TrackingParameters.idBoostSpread = compute_feature_spread(features);
else
    disp('Load idBoost Spread')
    idBoostSpread = load(idBoostSpreadName); 
    idBoostSpread = idBoostSpread.stdev;
    TrackingParameters.idBoostSpread = idBoostSpread;
end

%% Perform Tracking

[~, boundingBoxes, ~, ~] = track_instances_in_video(testVideoName,...
                                                    VideoParameters,...
                                                    CategoryDetector,...
                                                    TrackingParameters);



