function citeme(fname)
%CITEME Display citation requests
% Processes CITATION_REQUEST blocks in the first comment block of a
% function. Must be called at citeme(mfilename). 
% Citation requests are displayed only once per session, and look like
% this:
%
% CITATION_REQUEST:
% AUTHOR: Author Name, Author Name, Author Name, Author Name, Author Name, 
% Author Name and Author Name
% TITLE: A very long and complicated title that necessarily goes over 
% several lines
% JOURNAL: A renowned international journal that is so full of itself that 
% its name is very long
% VOL: V
% NUM: N
% YEAR: 2013
% PAGE: 1-10
% LINK: http://www.journal-website.com/dir/dir/some_strange_numbers_and_
%       letters_that_goes_over_several_lines
% PDF: ftp://myserver.com/paper1.pdf
% DOI: 1234/qwre/123
% PUBMED: 12321345
% ARXIV:  0706.0001
%
% A CITATION_REQUEST block is terminated by an empty line. 
% The order of fields is irrelevant. All fields are optional.
% At the moment, only one CITATION_REQUEST block per file is supported.
%
% See also STARTUP
if ischar(fname) && strcmp(fname, 'UNIT_TEST')
    citeme(mfilename);
    return;
end
cm = warning('query','EIDORS:CITEME'); % one would think this is supported
cf = warning('query',sprintf('EIDORS:CITEME:%s',fname));
if strcmp(cm.state, 'off') || strcmp(cf.state, 'off');
    return
end
h = help(fname);
s = strfind(h,'CITATION_REQUEST:');
if isempty(s), return, end;
T = textscan(h(s:end),'%s','delimiter','\n');
T = T{1};
line = 2;
lastfld = [];
cite = empty_cite;
while line <= numel(T) && ~isempty(strtrim(T{line}))
    [tok rest] = strtok(T{line});
    if strfind(tok,':');
            fld = lower(tok(1:end-1));
            cite.(fld) = strtrim(rest);
            lastfld = fld;
    elseif ~isempty(lastfld)
        cite.(lastfld) = sprintf('%s %s%s', cite.(lastfld), tok, rest);
    else
        error('Error processing citation request of file %s',fname);
    end
    line = line + 1;
end
flds = {'link','pdf','pubmed','arxiv','doi'};
for i = 1:numel(flds)
    cite.(flds{i})(isspace(cite.(flds{i}))) = [];
end
pretty_print(cite,fname);
 
function str = pretty_print(cite, fname)
str = sprintf('%s (%s) "%s"',...
    cite.author,cite.year,cite.title);

if ~isempty(cite.link)
    str = sprintf('%s <a href="matlab: web %s -browser">%s</a>%s', str, cite.link, cite.journal);
else
    str = sprintf('%s %s',str, cite.journal);
end

if ~isempty(cite.vol)
    str = sprintf('%s %s',str,cite.vol);
end
if ~isempty(cite.num)
    if ~isempty(cite.vol)
        str = sprintf('%s(%s)',str, cite.num);
    else
        str = sptrinf('%s %s',str, cite.num);
    end
end
if ~isempty(cite.page)
    str = sprintf('%s: %s',str, cite.page);
end
if ~isempty(cite.pdf)
    str = sprintf('%s <a href="matlab: web %s -browser">PDF</a>',str, cite.pdf);
end
if ~isempty(cite.doi)
    str = sprintf('%s DOI:<a href="matlab: web http://dx.doi.org/%s -browser">%s</a>',...
                    str, cite.doi,cite.doi);
end
if ~isempty(cite.pubmed)
    str = sprintf('%s PubMed:<a href="matlab: web http://www.ncbi.nlm.nih.gov/pubmed/%s -browser">%s</a>',...
                    str, cite.doi,cite.doi);
end

if ~isempty(cite.arxiv)
    str = sprintf('%s arXiv:<a href="matlab: web http://arxiv.org/abs/%s -browser">%s</a>',...
                    str, cite.doi,cite.doi);
end

str = strtrim(str);                
sp = strfind(str,' ');
tp1 = strfind(str,'<');
tp2 = strfind(str,'>');
sp = [sp length(str)];

rsp = sp; % real spaces
idx = false(size(sp));
for i = 1:length(tp1)
    nidx =  rsp>tp1(i) & rsp<tp2(i);
    idx = idx | nidx;
    sp(rsp>tp1(i)) = sp(rsp>tp1(i)) - (tp2(i) + 2- tp1(i));
end
sp(idx) = [];
rsp(idx) = [];

[jnk, nl] = unique(floor(sp/82),'last');
nl = rsp(nl);
nl(nl==length(str)) = [];
str(nl) = sprintf('\n');

ws = warning('off', 'backtrace');
idstr = sprintf('EIDORS:CITEME:%s',fname);
stars(1:80) = '=';
msg = sprintf('\n%s\nIf you use %s in a publication, please cite:\n', stars,upper(fname));
warning(idstr,[msg str '\n' stars]);
if 0
   disp('Press any key to continue...');
   pause
else
   fprintf('Continuing in  ');
   for i = 2:-1:1
      fprintf('\b%d',i);
      pause(1);
   end
   fprintf('\b0\n');
end
warning(ws.state, 'backtrace');
if ~strcmp(fname,'citeme')
    warning('off',idstr); % only show once per session
end

function cite = empty_cite
cite.author = [];
cite.year   = [];
cite.journal= [];
cite.title  = [];
cite.vol    = [];
cite.num    = [];
cite.page   = [];
cite.link   = [];
cite.pdf    = [];
cite.doi    = [];
cite.pubmed = [];
cite.arxiv  = [];
