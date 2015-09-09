function [bbs, patches] = acf_detect_all(I, detector, fileName )
% Run aggregate channel features object detector on given image(s).
%
% The input 'I' can either be a single image (or filename) or a cell array
% of images (or filenames). In the first case, the return is a set of bbs
% where each row has the format [x y w h score] and score is the confidence
% of detection. If the input is a cell array, the output is a cell array
% where each element is a set of bbs in the form above (in this case a
% parfor loop is used to speed execution). If 'fileName' is specified, the
% bbs are saved to a comma separated text file and the output is set to
% bbs=1. If saving detections for multiple images the output is stored in
% the format [imgId x y w h score] and imgId is a one-indexed image id.
%
% A cell of detectors trained with the same channels can be specified,
% detected bbs from each detector are concatenated. If using multiple
% detectors and opts.pNms.separate=1 then each bb has a sixth element
% bbType=j, where j is the j-th detector, see bbNms2.m for details.
%
% USAGE
%  bbs = acfDetect( I, detector, [fileName] )
%
% INPUTS
%  I          - input image(s) of filename(s) of input image(s)
%               or video structure
%  detector   - detector(s) trained via acfTrain
%  fileName   - [] target filename (if specified return is 1)
%  pPatch     - [0] return detected feature vectors if pPatch==1.
%
% OUTPUTS
%  bbs        - [nx5] array of bounding boxes or cell array of bbs
%
% EXAMPLE
%
% See also acfTrain, acfModify, bbGt>loadAll, bbNms2
%
% Copyright 2013 Piotr Dollar and Ron Appel. (Original Source Code)
%
% Copyright 2013 David Hall.  [dhall-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License
%
% Modified Dec-2013 David Hall
%   *changed to support video inputs
%   *returns feature vectors as well as bbs

% run detector on every image
if(nargin<3), fileName=''; end;
if(~isfield(detector.opts, 'returnPatch')), detector.opts.returnPatch=0; end;
multiple=iscell(I); video=isstruct(I);
if(~isempty(fileName) && exist(fileName,'file')), bbs=1; return; end
if(video)
    n=length(I.frmInds);
    [bbs,patches]=acfDetectVid(I,detector,fileName);
elseif(multiple),
    n=length(I); bbs=cell(n,1);
    parfor i=1:n, [bbs{i},patches{i}]=acfDetectImg(I{i},detector); end
else
    [bbs,patches]=acfDetectImg(I,detector);
end

% write results to disk if fileName specified
if(isempty(fileName)), return; end
d=fileparts(fileName); if(~isempty(d)&&~exist(d,'dir')), mkdir(d); end
if( multiple || video ) % add image index to each bb and flatten result
    for i=1:n, bbs{i}=[ones(size(bbs{i},1),1)*i bbs{i}]; end
    bbs=cell2mat(bbs);
end
dlmwrite(fileName,bbs); bbs=1;

end

function [bbs,patches] = acfDetectImg( I, detector)
% Run trained sliding-window object detector on given image.
Ds=detector; Clfs=detector.clf; if(~iscell(Clfs)), Clfs={Clfs}; end
% ROI=detector.ROI;
opts=Ds.opts; pPyramid=opts.pPyramid; pNms=opts.pNms;
pNms.separate=1; patches=[]; nClfs=length(Clfs);
modelDsPad=opts.modelDsPad; modelDs=opts.modelDs;
imreadf=opts.imreadf; imreadp=opts.imreadp;
shrink=pPyramid.pChns.shrink; pad=pPyramid.pad;
shift=(modelDsPad-modelDs)/2-pad;
% perform actual computations
if(all(ischar(I))), I=feval(imreadf,I,imreadp{:}); end
P = chnsPyramid(I,pPyramid); [bbs,scale,bbso] = deal(cell(P.nScales,nClfs));

for i=1:P.nScales
    bb = acf_detect_all1(P.data{i},Clfs,shrink,...
        modelDsPad(1),modelDsPad(2),opts.stride);
    bbso{i}=bb;
    bb(:,1)=(bb(:,1)+shift(2))/P.scaleshw(i,2);
    bb(:,2)=(bb(:,2)+shift(1))/P.scaleshw(i,1);
    bb(:,3)=modelDs(2)/P.scales(i);
    bb(:,4)=modelDs(1)/P.scales(i);
    bbs{i}=bb; scale{i}=i.*ones(size(bb,1),1);
end; bbs=cat(1,bbs{:}); bbso=cat(1,bbso{:}); scale=cat(1,scale{:});
if(~isempty(pNms)),
    [bbs,ind] = bb_nms_all(bbs,pNms);
    if(opts.returnPatch),
        patches=getPatches(P.data,bbso,scale,ind,shrink);
    end;
end
end

function [bbs,patches] = acfDetectVid( I, detector, fileName)
frmInd = I.frmInds; vid = I.vid; split=500;
n=numel(frmInd); [bbs,patches]=deal(cell(n,1));
frmInds=chunkify(1:n,split); N=length(frmInds);
ticId=ticStatus('Begin..');
for k=1:N,
    I=read(vid,[frmInds{k}(1),frmInds{k}(end)]);
    I=squeeze(mat2cell(I,size(I,1),size(I,2),size(I,3),ones(1,size(I,4))));
    [bbs1,patch] = acf_detect_all(I, detector, fileName);
    bbs(frmInds{k})=bbs1; patches(frmInds{k})=patch;
    tocStatus(ticId,k/N);
end
end

function [patches]=getPatches(data,bbs,scale,ind,shrink)
bb=bbs(ind,:); scale=scale(ind);
bb=bb./shrink;
patches=cell(1,size(bb,1));
for i=1:size(bb,1),
    patches{i}=data{scale(i)}((bb(i,2)+1):(bb(i,2)+bb(i,4)),(bb(i,1)+1):(bb(i,1)+bb(i,3)),:);
end
patches=cellfun(@(x) reshape(x,1,[]),patches,'UniformOutput',false);
patches=cat(1,patches{:});
end

