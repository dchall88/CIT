function [frame_number, gt_bb] = get_first_frame(gt_name, individual_id)

gt = vbb('vbbLoad', gt_name);

if(individual_id > max(gt.objId))
    error(['There are only ' num2str(max(gt.objId)) 'individuals in this video. Change individual_id accordingly']);
end

frame_number = gt.objStr(find(gt.objId == individual_id, 1));
gt_bb = [gt.objLists{frame_number}([gt.objLists{frame_number}.id] == individual_id).pos 0];


