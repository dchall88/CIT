function detector = train_category_detector(imDir,gtDir,resDir,modelName)
% imDir - directory of images
% gtDir - directory of text files (for gt format refer to acfTrain())
% resDir - directory to save detector
% modelName - name to give detector
% Refer to acfDemoCal.m for details on how to train your own category_detector.
%% set up opts for training detector (see acfTrain)
opts=acfTrain();
% opts.modelDs=[132 99]; opts.modelDsPad=[168 124];
opts.modelDs=[36 36]; opts.modelDsPad=[45 45];
opts.pPyramid.smooth=1; opts.pPyramid.pChns.pColor.smooth=1;
opts.pPyramid.pChns.shrink=2;
opts.nPerNeg=3;
opts.posGtDir=[resDir '/' gtDir]; opts.nWeak=[32 128 512 2048];
opts.posImgDir=[resDir '/' imDir]; opts.pJitter=struct('flip',1);
opts.pBoost.pTree.fracFtrs=1/16; 
opts.pBoost.pTree.maxDepth=1;
opts.name=[resDir '/' modelName];

%% train detector (see acfTrain)
detector = acfTrain( opts );

end
