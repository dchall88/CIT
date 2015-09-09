function out = vid_io( fName, action, varargin )
% Utilities for reading and writing video files.
%
% vid_io contains a number of utility functions for working with video files.
% The format for accessing the various utility functions is:
% out = vid_io( fName, 'action', inputs );
% The list of functions and help for each is given below. Also, help on
% individual subfunctions can be accessed by: "help vid_io>action".
%
% The external packages mmread and mmwrite by Micah Richert, available at 
% Matlab Central, are required.
%
% Create interface sr for reading video files.
%   sr = vid_io( fName, 'reader', [cache] )
% Create interface sw for writing video files.
%   sw = vid_io( fName, 'writer', info )
% Get info about video file.
%   info = vid_io( fName, 'getInfo' )
% Crop sub-sequence from video file.
%   vid_io( fName, 'crop', tName, frames )
% Extract images from video file to target directory or array.
%   Is = vid_io( fName, 'toImgs', [tDir], [skip], [f0], [f1], [ext] )
% Create video file from an array or directory of images
%   vid_io( fName, 'frImgs', info, varargin )
% Convert video file by applying imgFun(I) to each frame I.
%   vid_io( fName, 'convert', tName, imgFun, varargin )
% Create interface sr for reading dual video files.
%   sr = seqIo( fNames, 'readerDual', [cache] )
%
% USAGE
%  out = vid_io( fName, action, varargin )
%
% INPUTS
%  fName      - video file to open (anything mmread supports)
%  action     - controls action (see above)
%  varargin   - additional inputs (see above)
%
% OUTPUTS
%  out       - depends on action (see above)
%
% EXAMPLE
%
% See also vid_io>reader, vid_io>writer, vid_io>getInfo, vid_io>crop,
% vid_io>toImgs, vid_io>frImgs, vid_io>convert, vid_io>newHeader,
% vid_io>readerDual, vid_player, vid_reader_plugin, vid_writer_plugin
%
% Copyright 2013 Piotr Dollar. (Original Source Code)
%
% Copyright 2013 David Hall.  [dhall-at-caltech.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Simplified BSD License
%
% Modified Dec-2013 David Hall
%   *changed to support video formats avi, mp4, etc

switch lower(action)
  case {'reader','r'}, out = reader( fName, varargin{:} );
  case {'writer','w'}, out = writer( fName, varargin{:} );
  case 'getinfo', out = getInfo( fName );
  case 'crop', crop( fName, varargin{:} ); out=1;
  case 'toimgs', out = toImgs( fName, varargin{:} );
  case 'frimgs', frImgs( fName, varargin{:} ); out=1;
  case 'convert', convert( fName, varargin{:} ); out=1;
  case {'readerdual','rdual'}, out=readerDual(fName,varargin{:});
  otherwise, error('vid_io unknown action: ''%s''',action);
end
end

function sr = reader( fName, cache )
% Create interface sr for reading video files.
%
% Create interface sr to video file with the following commands:
%  sr.close();            % Close video file (sr is useless after).
%  I = sr.getframe();  % Get current frame (returns [] if invalid).
%  ts = sr.getts();       % Return timestamps for all frames.
%  info = sr.getinfo();   % Return struct with info about video.
%  I = sr.getnext();   % Shortcut for next() followed by getframe().
%  out = sr.next();       % Go to next frame (out=0 on fail).
%  out = sr.seek(frame);  % Go to specified frame (out=0 on fail).
%  out = sr.step(delta);  % Go to current frame+delta (out=0 on fail).
%
% If cache>0, reader() will cache frames in memory, so that calls to
% getframe() can avoid disk IO for cached frames (note that only frames
% returned by getframe() are cached). This is useful if the same frames are
% accessed repeatedly. When the cache is full, the frame in the cache
% accessed least recently is discarded. Memory requirements are
% proportional to cache size.
%
% USAGE
%  sr = vid_io( fName, 'reader', [cache] )
%
% INPUTS
%  fName  - video file name
%  cache  - [10] size of cache
%
% OUTPUTS
%  sr     - interface for reading video file
%
% EXAMPLE
%
% See also vid_io, vid_reader_plugin
if(nargin<2 || isempty(cache)), cache=10; end
if( cache>0 ), [fs, Is, inds]=deal([]); end
r=@vid_reader_plugin; s=r('open',int32(-1),fName);
sr = struct( 'close',@() r('close',s), 'getframe',@getframe, ...
  'getts',@() r('getts',s), ...
  'getinfo',@() r('getinfo',s), 'getnext',@() r('getnext',s), ...
  'next',@() r('next',s), 'seek',@(f) r('seek',s,f), ...
  'step',@(d) r('step',s,d));

  function I=getframe()
    % if not using cache simply call 'getframe' and done
    if(cache<=0), I=r('getframe',s); return; end
    % if cache initialized and frame in cache perform lookup
    f=r('getinfo',s); numFrms=f.NumberOfFrames; 
    f=f.curFrame; i=find(f==fs,1);
    if(i), I=Is(inds{:},i); return; end
    % if image not in cache add
    fs=f:min(numFrms,f+cache-1); I=r('getframes',s,[fs(1):fs(end)]); 
    Is=I; I=Is(:,:,:,1); 
    inds=repmat({':'},1,ndims(I));
  end
end

function sw = writer( fName, info )
% Create interface sw for writing video files.
%
% Create interface sw to video file with the following commands:
%  sw.close();              % Close video file (sw is useless after).
%  sw.addframe(I,[ts]);     % Writes video frame (and timestamp)
%  info = sw.getinfo();     % Return struct with info about video.
%
% The following params must be specified in struct 'info' upon opening:
%  width          - frame width
%  height         - frame height
%  fps            - frames per second
%
% USAGE
%  sw = vid_io( fName, 'writer', info )
%
% INPUTS
%  fName  - video file name
%  info   - see above
%
% OUTPUTS
%  sw     - interface for writing video file
%
% EXAMPLE
%
% See also seqIo, seqWriterPlugin
if(nargin<2), info=[]; end
w=@vid_writer_plugin; s=w('open',int32(-1),fName,info);
sw = struct( 'close',@() w('close',s), 'getinfo',@() w('getinfo',s), ...
  'addframe',@(varargin) w('addframe',s,varargin{:}));
end

function info = getInfo( fName )
% Get info about video file.
%
% USAGE
%  info = vid_io( fName, 'getInfo' )
%
% INPUTS
%  fName  - video file name
%
% OUTPUTS
%  info   - information struct
%
% EXAMPLE
%
% See also seqIo
sr=reader(fName); info=sr.getinfo(); sr.close();
end

function crop( fName, tName, frames )
% Crop sub-sequence from video file.
%
% Frame indices are 1 indexed. frames need not be consecutive and can
% contain duplicates. An index of 0 indicates a blank (all 0) frame.
%
% USAGE
%  vid_io( fName, 'crop', tName, frames )
%
% INPUTS
%  fName      - video file name
%  tName      - cropped video file name
%  frames     - frame indices (1 indexed)
%
% OUTPUTS
%
% EXAMPLE
%
% See also vid_io
sr=reader(fName); info=sr.getinfo(); sw=writer(tName,info);
frames=frames(:)'; pad=sr.getnext(); pad(:)=0;
kp=frames>=0 & frames<=abs(info.NumberOfFrames); if(~all(kp)), frames=frames(kp);
  warning('piotr:vid_io:crop','%i out of bounds frames',sum(~kp)); end
n=length(frames); k=0; tid=ticStatus;
for f=frames
  if(f==0), sw.addframe(pad); continue; end
  sr.seek(f); I=sr.getframe(); k=k+1; tocStatus(tid,k/n);
  sw.addframe(I); 
end; sw.close(); sr.close();
end

function Is = toImgs( fName, tDir, skip, f0, f1, ext, frameArray )
% Extract images from video file to target directory or array.
%
% USAGE
%  Is = vid_io( fName, 'toImgs', [tDir], [skip], [f0], [f1], [ext] )
%
% INPUTS
%  fName      - vid file name
%  tDir       - [] target directory (if empty extract images to array)
%  skip       - [1] skip between written frames
%  f0         - [1] first frame to write
%  f1         - [numFrames] last frame to write
%  ext        - ['.png'] optionally save as given type (slow, reconverts)
%  frameArray - [] optionally specify exact frames to extract
%
% OUTPUTS
%  Is         - if isempty(tDir) outputs image array (else Is=[])
%
% EXAMPLE
%
% See also vid_io
if(nargin<2 || isempty(tDir)), tDir=[]; end
if(nargin<3 || isempty(skip)), skip=1; end
if(nargin<4 || isempty(f0)), f0=1; end
if(nargin<5 || isempty(f1)), f1=inf; end
if(nargin<6 || isempty(ext)), ext='png'; end
if(nargin<7 || isempty(frameArray)), frameArray=[]; end

sr=reader(fName); info=sr.getinfo(); f1=min(f1,info.NumberOfFrames);
if(isempty(frameArray)), frames=f0:skip:f1; else frames=frameArray; end 
n=length(frames); tid=ticStatus; k=0;
% output images to array
if(isempty(tDir))
  I=sr.getnext(); d=ndims(I); assert(d==2 || d==3);
  try Is=zeros([size(I,1)/2 size(I,2)/2 size(I,3) n],class(I)); catch e; sr.close(); throw(e); end
  for k=1:n, sr.seek(frames(k)); I=sr.getframe(); tocStatus(tid,k/n);
    if(d==2), Is(:,:,k)=I; else Is(:,:,:,k)=imresize(I,0.5); end; end
  sr.close(); return;
end
% output images to directory
if(~exist(tDir,'dir')), mkdir(tDir); end; Is=[];
for frame=frames
  f=[tDir '/I' int2str2(frame,5) '.']; sr.seek(frame);
  I=sr.getframe(); imwrite(I,[f ext]); 
  k=k+1; tocStatus(tid,k/n);
end; sr.close();
end

function frImgs( fName, info, varargin )
% Create video file from an array or directory of images.
%
% For info, if converting from array, only codec (e.g., 'jpg') and fps must
% be specified while width and height and determined automatically.
%
% USAGE
%  vid_io( fName, 'frImgs', info, varargin )
%
% INPUTS
%  fName      - video file name
%  info       - defines codec, etc, see seqIo>writer
%  varargin   - additional params (struct or name/value pairs)
%   .Is         - [] if specified create video from image array
%   .sDir       - [] source directory
%   .skip       - [1] skip between frames
%   .name       - ['I'] base name of images
%   .ext        - ['png'] extension of images
%   .nDigits    - [5] number of digits for filename index
%   .f0         - [0] first frame to read
%   .f1         - [10^6] last frame to read
%
% OUTPUTS
%
% EXAMPLE
%
% See also vid_io, vid_io>writer
dfs={'Is',[],'sDir',[],'skip',1,'name','I',...
  'ext','png','nDigits',5,'f0',0,'f1',10^6};
[Is,sDir,skip,name,ext,nDigits,f0,f1] ...
  = getPrmDflt(varargin,dfs,1);
if( isempty(Is) )
  assert(exist(sDir,'dir')==7); sw=writer(fName,info);
  f1=length(dir([sDir '/*.' ext]))+f0-1;
  frmStr=sprintf('%s/%s%%0%ii.%s',sDir,name,nDigits,ext);
  for frame = f0:skip:f1
    f=sprintf(frmStr,frame);
    if(~exist(f,'file')), sw.close(); assert(false); end
    I=imread(f); sw.addframe(I);
  end; sw.close();
  if(frame==f0), warning('No images found.'); end %#ok<WNTAG>
else
  nd=ndims(Is); if(nd==2), nd=3; end; assert(nd<=4); nFrm=size(Is,nd);
  info.height=size(Is,1); info.width=size(Is,2); sw=writer(fName,info);
  if(nd==3), for f=1:nFrm, sw.addframe(Is(:,:,f)); end; end
  if(nd==4), for f=1:nFrm, sw.addframe(Is(:,:,:,f)); end; end
  sw.close();
end
end

function convert( fName, tName, imgFun, varargin )
% Convert video file by applying imgFun(I) to each frame I.
%
% USAGE
%  vid_io( fName, 'convert', tName, imgFun, varargin )
%
% INPUTS
%  fName      - video file name
%  tName      - converted video file name
%  imgFun     - function to apply to each image
%  varargin   - additional params (struct or name/value pairs)
%   .info       - [] info for target video file
%   .skip       - [1] skip between frames
%   .f0         - [1] first frame to read
%   .f1         - [inf] last frame to read
%   .prm        - {} parameters for ImgFun
%
% OUTPUTS
%
% EXAMPLE
%
% See also seqIo
dfs={'info',[],'skip',1,'f0',1,'f1',inf,'prm',{}};
[info,skip,f0,f1,prm]=getPrmDflt(varargin,dfs,1);
assert(~strcmp(tName,fName)); sr=reader(fName); infor=sr.getinfo();
if(isempty(info)), info=infor; end; n=infor.NumberOfFrames; f1=min(f1,n);
I=sr.getnext(); I=imgFun(I,prm.dt{1},prm.p{1}); info.width=size(I,2); info.height=size(I,1);
sw=writer(tName,info); tid=ticStatus('converting video');
frames=f0:skip:f1; n=length(frames); k=0;
for f=frames, 
  sr.seek(f); I=sr.getframe(); I=imgFun(I,prm.dt{f},prm.p{f});
  sw.addframe(I);
  k=k+1; tocStatus(tid,k/n);
end; sw.close(); sr.close();
end

function sr = readerDual( fNames, cache )
% Create interface sr for reading dual video files.
%
% Wrapper for two video files of the same image dims and roughly the same
% frame counts that are treated as a single reader object. getframe()
% returns the concatentation of the two frames. For videos of different
% frame counts, the first video serves as the "dominant" video and the
% frame count of the second video is adjusted accordingly. Same general
% usage as in reader, but the only supported operations are: close(),
% getframe(), getinfo(), and seek().
%
% USAGE
%  sr = vid_io( fNames, 'readerDual', [cache] )
%
% INPUTS
%  fNames - two seq file names
%  cache  - [10] size of cache (see seqIo>reader)
%
% OUTPUTS
%  sr     - interface for reading video file
%
% EXAMPLE
%
% See also seqIo, seqIo>reader
if(nargin<2 || isempty(cache)), cache=10; end
s1=reader(fNames{1}, cache); i1=s1.getinfo();
s2=reader(fNames{2}, cache); i2=s2.getinfo();
info=i1; info.Width=i1.Width+i2.Width;
if( i1.Width~=i2.Width || i1.Height~=i2.Height )
  s1.close(); s2.close(); error('Mismatched videos'); end
if( i1.NumberOfFrames~=i2.NumberOfFrames )
  warning('video files of different lengths'); end %#ok<WNTAG>
frame2=@(f) round(f/(i1.NumberOfFrames-1)*(i2.NumberOfFrames-1));

sr=struct('close',@() min(s1.close(),s2.close()), ...
  'getframe',@getframe, 'getinfo',@() info, ...
  'seek',@(f) s1.seek(f) & s2.seek(frame2(f)) );

  function I=getframe()
    I1=s1.getframe(); I2=s2.getframe(); I=[I1 I2]; end
end