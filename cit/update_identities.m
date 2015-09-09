function bboxes = update_identities(bboxes, active_detectors)
%Updates bounding box id number with correct id

bboxes_temp = bboxes;
for k=2:(numel(active_detectors)+1)
	bboxes_temp(bboxes(:,6)==k, 6) = active_detectors(k-1);
end
bboxes = bboxes_temp;


end