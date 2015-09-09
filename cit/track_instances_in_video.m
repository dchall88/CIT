function [boundingBoxesPerInstance, boundingBoxes, InstanceDetectors, time] = ...
track_instances_in_video(videoName, VideoParameters, CategoryDetector,  TrackingParameters)
%
% INPUTS
%  videoName        - ['REQ'] Directory of Video (any format VideoReader supports)
%  VideoParameters
%     .frames     - []  Frames of video to run tracking system on
%     .start      - [1] Frame to start processing from (if frames=[])
%     .stride     - [1]  Step between frames (if frames=[])
%  CategoryDetector  - category detector trained via acfTrain
%  TrackingParameters
%     .idBoostBeta - ['REQ'] tuning parameter for idBoost (see create_individual_clf.m) 
%     .idBoostSpread - ['REQ'] parameter for idBoost (see create_individual_clf.m) 
%     .destroyAfterNFrames - [5] Destroy individual detector after this number of frames.
%     .minCategoryDetectorScore - [100] Minimum category detector score to initialise an individual detector
%     .minDetectionOverlap - [0.4] Minimum overlap between detections
%     .VIEW_RESULTS - [false] whether to display results while running
%
% OUTPUTS
%  boundingBoxesPerInstance - {nFrames}[nBoxesx6] cell array of bounding
%                               boxes for all instance and category
%                               detectors.
%  boundingBoxes - cell array of combined bounding boxes
%  InstanceDetectors - structure with details about each instance detector
%                       (consumes a bit of memory)
%  time - the time it took to run tracking on the video

if ~isfield(CategoryDetector.opts, 'returnPatch')
    CategoryDetector.opts.returnPatch=1;
end

if ~isfield(CategoryDetector.clf, 'cascThr')
    CategoryDetector.clf.cascThr = CategoryDetector.opts.cascThr;
end

%% Initialise Tracking Parameters
defaults = {'idBoostBeta', 'REQ', 'idBoostSpread', 'REQ',...
			'destroyAfterNFrames', 5, 'minCategoryDetectorScore', 0,...
			'minDetectionOverlap',0.5, 'VIEW_RESULTS', false};
[idBoostBeta,... 
 idBoostSpread,...
 destroyAfterNFrames,...
 minCategoryDetectorScore,...
 minDetectionOverlap,...
 VIEW_RESULTS] = getPrmDflt(TrackingParameters, defaults, 0);

TrackingParameters.destroyAfterNFrames = destroyAfterNFrames; 
TrackingParameters.idBoostBeta = idBoostBeta;
TrackingParameters.idBoostSpread = idBoostSpread;
TrackingParameters.clf = CategoryDetector.clf;

%% Initialise Video
defaults = {'frames',[], 'start',1, 'stride',1};
[frames, start, stride] = getPrmDflt(VideoParameters, defaults, 1);
vid = vid_io(videoName, 'reader');
info = vid.getinfo();
if(isempty(frames)), 
    frames = start:stride:info.('NumberOfFrames'); 
end
nFrames = numel(frames);

%% Initialise Variables
[InstanceDetectors, active_detectors] = deal([]);
[boundingBoxesPerInstance, boundingBoxes] = deal(cell(0, 1));
CategoryDetector.opts.pNms.separate = 1;

ticId=ticStatus('Begin...'); tic;
for iFrame = 1:nFrames
    
	%% Get Next Frame
    vid.seek(frames(iFrame));
    I = vid.getframe();
    
    %% Combine Detectors
    detector = combine_detectors(CategoryDetector, InstanceDetectors, active_detectors);
    
    %% Evaluate Category and Individual Detectors on frame I
    [bboxes1, patch1] = acf_detect_all(I, detector);
    
    %% Give Detections the Correct ID Number
    bboxes1 = update_identities(bboxes1, active_detectors); 
	boundingBoxesPerInstance{iFrame} = bboxes1;
    
    %% Combine Category and Individual Detections
    [Bboxes1, Patch1] = combine_category_and_instance_detections(bboxes1, patch1,  minDetectionOverlap, minCategoryDetectorScore);
    if ~isempty(Bboxes1.both) || ~isempty(Bboxes1.instance)
        Bboxes1.cat = []; 
		Patch1.cat = [];
    else
        if size(Bboxes1.cat, 1) > 1
            [~, maxind] = max(Bboxes1.cat(:, 5)); 
			Bboxes1.cat = Bboxes1.cat(maxind, :); 
            Patch1.cat = Patch1.cat(maxind, :);
        end
    end
    
    
    %% Track Individuals
    [Bboxes1, InstanceDetectors, active_detectors] = ...
		track_cit(frames(iFrame), Bboxes1, Patch1, InstanceDetectors, active_detectors, TrackingParameters);
    
    %% Prepare Results for Output
    boundingBoxes{iFrame} = [Bboxes1.cat; Bboxes1.both; Bboxes1.instance];
       
    %% View Results
    if VIEW_RESULTS == 1
        view_tracker(I, boundingBoxesPerInstance{iFrame}, boundingBoxes{iFrame},...
					InstanceDetectors, frames(iFrame));
        pause(1 / info.fps);
    end
    tocStatus(ticId, iFrame/nFrames);
end
time=toc;

end

	function Detector = combine_detectors(CategoryDetector, InstanceDetectors, active_detectors)
		Detector = CategoryDetector;
		if ~isempty(InstanceDetectors)
			Detector.clf = {CategoryDetector.clf, InstanceDetectors(active_detectors-1).clf};
		end
	end