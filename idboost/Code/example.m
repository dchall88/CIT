%% An simple example script demonstrating the use of a category-to-individual detector. 
% An individual is identified by a category detector. 
% This single example is used to generate an individual detector. 
% The individual detector is then evaluated on the entire video.  
% Results are then shown.
%
% Change 'individual_id' to experiment with different individuals in the dataset.
%
% The IDBoost algorithm is implemented in 'create_individual_clf'

%% Parameters and Settings

category_detector_name = 'Data/Detector/acfFPOQDetector.mat';

train_vid_name = 'Data/Datasets/FPOQ/Videos/Berlinale-uv5l-uu6SIw.mp4';
train_gt_name = 'Data/Datasets/FPOQ/Annotations/Berlinale-uv5l-uu6SIw.vbb';

test_vid_name = 'Data/Datasets/FPOQ/Videos/Alabama-csHddXn91YE.mp4';
test_gt_name = 'Data/Datasets/FPOQ/Annotations/Alabama-csHddXn91YE.vbb';

% The individual's id number, the first frame they occur in and their
% ground truth bounding box (to ensure we chose the correct exemplar if
% their are multiple people in the frame.)
individual_id = 1;
[frame_number, gt_bb] = get_first_frame(test_gt_name, individual_id);

% Tuning parameter for IDBoost
beta = 0.6;

%% Load Category Detector
disp('Load Category Detector')
category_detector = load(category_detector_name); 
category_detector = category_detector.detector;
category_detector.opts.returnPatch = 1; 

%% Estimate sigma (the spread) from a training set. Refer to Eqn (4).
disp('Estimate Spread')

[features, thumbnails] = ...
    extract_feature_tracks(train_vid_name, train_gt_name, category_detector);

[spread] = compute_feature_spread(features);

% View Tracks and Histograms
if(1)
   track_number = 10;
   feature_number = category_detector.clf.fids(1);
   figure(1); montage(thumbnails{track_number});
   title('Track of an Individual from the Training Set used to Estimate Sigma');
   figure(2); ksdensity(features{track_number}(:, feature_number));
   title('Density Estimate of First Feature used by the Category Detector for the Individual Track in Figure 1');
end
clear thumbnails

%% Example

% Test Video
disp('Load Test Video')
test_video = VideoReader(test_vid_name);
total_frames = test_video.NumberOfFrames;

% Take first exemplar of an individual (I am just taking the first example
% of the very first individual in the video - she starts at frame 496)
disp('Detect First Instance of Individual')
I = read(test_video, frame_number); 

[bbs_cat, patches] = acfDetect(I, category_detector);

% Choose example that corresponds to the correct individual
[~, matches] = bbGt('evalRes', gt_bb, bbs_cat(:,1:5));
individual_number = find(matches(:,6));

feature_vector = patches(individual_number, :); 

% View example of the individual
if(1)
    figure(3)
    Iexample = imcrop(I, bbs_cat(individual_number, 1:4));
    im(Iexample);
    title('Examplar used to train Individual Detector');
end

% Create Individual Detector for this example
disp('Create Individual Detector')
individual_detector = category_detector;
individual_detector.opts.returnPatch = 0; 
[individual_detector.clf] = ...
    create_individual_clf(feature_vector, category_detector.clf, spread, beta);

% Evaluate Individual Detector on the entire video (This may take some
% time.) If you wish to evaluate the detector on a smaller number of frames
% you can specify the frame indices in test_frames.
disp('Evaluate Individual Detector')
test_video_struct.vid = test_video;
test_frames = 1:total_frames;
test_video_struct.frmInds = test_frames;

bbs = acfDetect(test_video_struct, individual_detector);

% Makes the results compatible with bbGt->evalRes
for k = 1:length(bbs)
    bbs{k} = bbs{k}(:,1:5);
end

%% View results
disp('View Results')
test_gt = vbb('vbbLoad',test_gt_name);
stats = vbb('getStats', test_gt);

% Extract ground-truth information about the individual
frames = stats.frame(stats.ids == individual_id);
gt_bb = stats.pos(stats.ids == individual_id, :);
gt_occl = stats.occl(stats.ids == individual_id);
gt_bb = [gt_bb gt_occl'];

gt = cell(1, total_frames);
gt(:) = {zeros(0,5)};
for k = 1:length(frames)
    gt{frames(k)} = gt_bb(k, :);
end
gt = gt(test_video_struct.frmInds);

% Compute the ROC
ref = 10.^(-5:.25:0);
[gt1, dt1] = bbGt('evalRes', gt, bbs');
[fp, tp, score, miss] = bbGt('compRoc', gt1, dt1, 1, ref);
miss=exp(mean(log(max(1e-10, 1-miss))));
roc=[score, fp, tp];

% Plot ROC (If no curve is visible then change 'logx' to 0 and the first 
% value of 'lims' to 0) 
figure(4)
plotRoc([fp tp],'logx',1,'logy',1,'xLbl','fppi',...
    'lims',[10^-6 1e1 .005 1],'color','b','smooth',1,'fpTarget',ref);
title(sprintf('log-average miss rate = %.2f%%',miss*100));

% Show video with individual detections overlayed (Only show high scoring
% detections)
threshold_index = find(roc(:,3) == 1, 1); %Find where TPR = 1
if(roc(threshold_index, 2) < 10^-4) %If FPR is small for perfect TPR
    threshold = roc(threshold_index, 1);
else
    threshold_index = find(roc(:,3) > 0.8, 1); %Find where TPR = 0.8
    threshold = roc(threshold_index,1);
end

bbs_show = cell(1, total_frames);
bbs_show(:) = {zeros(0,6)};
for k = 1:length(bbs)
    valid = bbs{k}(:,5) >= threshold;
    bbs_show{test_frames(k)} = bbs{k}(valid,:);
    bbs_show{test_frames(k)}(:,6) = 1;
end
results_player(test_vid_name, bbs_show);

