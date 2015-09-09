function [Bboxes, InstanceDetectors, active_detectors] = track_cit(frame, Bboxes, Patches, InstanceDetectors, active_detectors, TrackingParameters)

	%% Create New Individual Detectors
	new_patches = Patches.cat; 
	bboxes_cat = Bboxes.cat;
	nNew = size(new_patches, 1); 
	nCurrent = size(InstanceDetectors, 1);
	newInstanceDetectors = create_detector(frame, nNew, nCurrent, new_patches, bboxes_cat, TrackingParameters);
	if ~isempty(newInstanceDetectors)
		Bboxes.cat(:,6) = [newInstanceDetectors.id]; 
	end

	%% Update Existing Individual Detectors with New Appearance Model
	renew_patches = Patches.both; 
	bboxes_both = Bboxes.both;
	InstanceDetectors = renew_detector(frame, InstanceDetectors, renew_patches, bboxes_both, TrackingParameters);
	Bboxes.both = Bboxes.both(:, 1:6);

	%% Update Existing Individual Detectors.
	% The appearance model is not updated since there is no corresponding category detection
	bboxes_instance = Bboxes.instance;
	InstanceDetectors = update_detector(frame, InstanceDetectors, bboxes_instance);

	%% Destroy Individual Detectors that have had no corresponding category
	% detection after a certain period of time
	ids = bboxes_both(:,6)';
	[active_detectors, InstanceDetectors] = destroy_detector(frame, InstanceDetectors, ids, active_detectors, TrackingParameters);

	%% Add New Detectors to Existing Ones
	InstanceDetectors = cat(1, InstanceDetectors, newInstanceDetectors);
	active_detectors = cat(1, active_detectors, [newInstanceDetectors.id]');
end

function detectors = create_detector(frame, nNew, nCurrent, patches, bboxes_cat, TrackingParameters)

	detectors(1:nNew,1)=deal(struct('clf',struct(),...
		'id',[],'frames',[],'F',[],...
		'score_cat',[],'bboxes_cat',[],'last_frame_cat',[],...
		'score_instance',[],'bboxes_instance',[],'last_frame_instance',[],...
		'bboxes_both',[]));

	%Assign ID numbers
	ids = num2cell((2:nNew+1) + nCurrent); %Add Total Number of Individuals
	[detectors(1:nNew).id] = deal(ids{:});

	%Create Classifiers
	for k=1:nNew
		F = patches(k,:); assert(isa(F,'single'));
		detectors(k).clf = create_individual_clf(F, TrackingParameters.clf, TrackingParameters.idBoostSpread, TrackingParameters.idBoostBeta);
		detectors(k).frames = frame;
		detectors(k).F = F;
		detectors(k).score_cat = bboxes_cat(k,5);
		detectors(k).bboxes_cat = bboxes_cat(k,1:5);
		detectors(k).last_frame_cat = frame;
		detectors(k).score_instance = -5;
		detectors(k).bboxes_both = bboxes_cat(k,1:5);
	end

end

function detectors = renew_detector(frame, detectors, patches, bboxes_both, TrackingParameters)

	%Determine which detectors to renew
	%Could do something fancy based on scores of category and individual
	%detection but for now will just renew all individual detectors that
	%overlap with a category detection regardless of score.
	ind=bboxes_both(:,6)-1; 
	N=length(ind);

	%Create Classifiers
	for k=1:N
		id = ind(k);
		F = patches(k,:); assert(isa(F,'single'));
		detectors(id).clf = create_individual_clf(F, TrackingParameters.clf, TrackingParameters.idBoostSpread, TrackingParameters.idBoostBeta);
		detectors(id).frames(end+1,1) = frame;
		detectors(id).F = [detectors(id).F;F];
		detectors(id).score_cat(end+1,1) = bboxes_both(k,11);
		detectors(id).bboxes_cat(end+1,:) = bboxes_both(k,7:11);
		detectors(id).last_frame_cat = frame;
		detectors(id).score_instance(end+1,1) = bboxes_both(k,16);
		detectors(id).bboxes_instance(end+1,:) = bboxes_both(k,12:16);
		detectors(id).last_frame_instance = frame;
		detectors(id).bboxes_both(end+1,:) = bboxes_both(k,1:5);
	end

end

function detectors =update_detector(frame, detectors, bboxes_instance)

	%Determine which detectors to update
	ind = bboxes_instance(:,6) - 1; 
	N=length(ind);

	%Update Detectors with Detection results
	for k = 1:N
		id = ind(k);
		detectors(id).frames(end+1,1) = frame;
		detectors(id).score_cat(end+1,1) = -5;
		detectors(id).score_instance(end+1,1) = bboxes_instance(k,5);
		detectors(id).bboxes_instance(end+1,:) = bboxes_instance(k,1:5);
		detectors(id).last_frame_instance = frame;
		detectors(id).bboxes_both(end+1,:) = bboxes_instance(k,1:5);
	end

end

function [active_detectors, detectors] = destroy_detector(frame, detectors, ids, active_detectors, TrackingParameters)

	destroyAfterNFrames = TrackingParameters.destroyAfterNFrames;
	[ids, ind] = setdiff(active_detectors, ids);
	if ~isempty(ids)
		v = (frame - [detectors(ids - 1).last_frame_cat]) > destroyAfterNFrames;
		if(any(v)), 
			[detectors(ids(v)-1).clf] = deal([]);
		end
		active_detectors( ind(v) ) = [];
	end
	if numel(active_detectors)==0
		active_detectors=[]; 
	end
end