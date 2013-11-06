function inds=cellfind(cellarray,str,lengths,hashtable)
%
%
% CELLFIND - Search for strings in a cell array
%
% INDS=CELLFIND(CELLARRAY,STR)
%   or
% INDS=CELLFIND(CELLARRAY,STR,LENGTHS,HASHTABLE)
%
% Inputs: 
%   CELLARRAY - a cell array to search in
%   STR       - either a string or a cell array of strings
%               when a cell array, all strings are searched for.
%               assuming each string would be found only once.
%   LENGTHS   - vector of lengths of words in c
%   HASHTABLE - vector of hash values of words in c
%
% Note: 
%   If strings appear more than once, only the first occurence will be returned
% 
% See also 
%   CELLFIND_PREP, for preparing hashtable & lengths in advance
% 
% Gal Chechik (C) 2004
% 
if(min(size(cellarray))>1)
  disp('min(size(c))>1');  
  disp('cellfind searches only inside vectors');
  return  
end
if(min(size(str))>1)
  disp('min(size(str))>1');    
  disp('cellfind cannot search for matrices in the cell array');
  return    
end

% For debug enzyme names
%for(i=1:length(cellarray))
%  tmp = cellarray{i};
%  if(length(find(tmp=='_'))>0)
%    if(length(tmp)<=8), keyboard, end
%    if(tmp(8)=='_' & tmp(10)=='_'), keyboard; end
%  end
%end


% Cell array str
if(iscell(str))

  %for(i=1:length(str))
  %  tmp = str{i};
  %  if(length(find(tmp=='_'))>0)
  %    if(length(tmp)<=8), keyboard, end
  %    if(tmp(8)=='_' & tmp(10)=='_'), keyboard; end
  %  end
  %end
  
  if(nargin<3)  
    [lengths,hashtable]=cellfind_prep(cellarray);
  end  
  
  inds = zeros(1,length(str));
  hashstrs = hash(str);
  for(i=1:length(str))
    ind = cellfind1(cellarray,str{i},lengths,hashtable,hashstrs(i));
    if(~isempty(ind))
      inds(i)=ind(1);
    end
  end
  
  return;
end

% Simple string str
if(nargin<3)
  inds=cellfind1(cellarray,str);
elseif(nargin<4)
  inds=cellfind1(cellarray,str,lengths);
else
  inds=cellfind1(cellarray,str,lengths,hashtable);
end

if(isempty(inds)), inds=0; end

return

% =============================================
function inds=cellfind1(c,str,lengths,hashtable,hashstr)

lenc=length(c);
  

if(nargin==5)
  % ============

  inds = find(hashtable==hashstr);
  ninds = length(inds);
  
  b = zeros(lenc,1);
  for(ii=1:ninds)
    i = inds(ii);
    ci = c{i};
    if(length(ci)==length(str))
      if(ci==str)
	b(i)=1;
      end
    end
  end
  inds = find(b);
     
elseif(nargin==2)
  % ============
  lstr=length(str);
  
  b = zeros(lenc,1);
  for(i=1:lenc)
    ci = c{i};
    if(length(ci)==lstr)
      if(ci==str)
	b(i)=1;
      end
    end    
  end
  inds = find(b);
  
elseif(nargin==3)
  % ============

  inds = find(lengths==length(str));
  ninds = length(inds);
  
  b = zeros(lenc,1);
  for(ii=1:ninds)
    i = inds(ii);
    ci = c{i};
    if(ci==str)
      b(i)=1;
    end
  end
  inds = find(b);
  
elseif(nargin==4)
  % ============
  hashstr = hash(str);
  inds = find(hashtable==hashstr);
  ninds = length(inds);
  
  b = zeros(lenc,1);
  for(ii=1:ninds)
    i = inds(ii);
    ci = c{i};
    if(length(ci)==length(str))
      if(ci==str)
	b(i)=1;
      end
    end
  end
  inds = find(b);

end

return
