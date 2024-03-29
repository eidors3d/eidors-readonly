function obj= matlab_xml_test( filename)
%  if nargin==0; filename = 'sample.xml'; end
   if nargin==0; filename = '../docs/samples/sampleGoeFile2.oeit/info/frame_types.xml'; end
      try
         tree = xmlread(filename);
      catch
         error('Failed to read XML file %s.',filename);
      end

switch 2
   case 1;
      obj = print_recurse( tree, 0);

   case 2;
      obj = xpath_try(tree);

   case 2;
%   http://www.mathworks.com/matlabcentral/answers/99993


   otherwise;
      error('boo')
end

function obj = xpath_try(tree)
% From:
% http://blogs.mathworks.com/community/2010/11/01/xml-and-matlab-navigating-a-tree/
   import javax.xml.xpath.*
   factory = XPathFactory.newInstance;
   xpath = factory.newXPath;

   i = 1;
   % compile and evaluate the XPath Expression
   while 1;
      expression = xpath.compile(sprintf( ...
            'frame_type_list/frame_type/acquisition[%d]/stim_list/stim',i));
%           'frame_type_list/frame_type/acquisition/stim_list/stim[%d]',i));
      stimnode = expression.evaluate(tree, XPathConstants.NODE);
      if isempty(stimnode); break; end
      node= parseChildNodes( stimnode);
      jj= 1; for j= 1: length(node)
         if strcmp( node(j).Name, '#text');  continue ; end
if DEBUG
         fprintf('%d: %d: %s [',i, j, node(j).Name);
         for k=1:length(node(j).Attributes);
            fprintf('%s=%s ', ...
              node(j).Attributes(k).Name, node(j).Attributes(k).Value);
         end
         fprintf(']\n');
end
         obj.acquisition(1).stim(jj) = node(j);
         jj= jj+1;
      end
     
      i = i+1;
   end



function d = DEBUG; d=0;

function obj = print_recurse( node, level);
%  if level>5; obj.name ='empty'; return;end
   
   obj = struct();

   [name, atts,  kids] = non_text_nodes( node);
   obj.name = name;
   % deftault fields
   obj.comment = '';
   switch name % this is a hack
      case 'electrode'; obj.contact_impedance = '???';
   end

if DEBUG
   for i=1:level; fprintf('>'); end
   fprintf('[ %s: ',name);
end
   for i=1:length(atts);
if DEBUG
      fprintf('%s="%s" ',atts(i).Name, atts(i).Value);
end
      obj.(atts(i).Name) = atts(i).Value;
   end
if DEBUG
   fprintf(']\n');
end
   for i=1:length(kids);
     obj_i = print_recurse( kids{i}, level+1);
     on = obj_i.name;
     if strcmp(on, 'empty'); continue; end
     if strcmp(on, 'user_data'); continue; end

%    try   ll = length(obj.( on ) );
%    catch ll = 0;
%    end
     ll= 0;
     if isfield(obj, on)
        ll = length(obj.( on ) );
     end

     obj.( obj_i.name )( ll+1 ) = orderfields( obj_i);
   end

   switch name % this is a hack
      case 'acquisition';
         if ~isfield( obj,'stim_list');
            obj.stim_list.stim.elec.id = '??';
            obj.stim_list.stim.elec.gain = '??';
         end
   end
   obj = orderfields( obj);


function [name, atts,  kids] = non_text_nodes( node)

   name = char( node.getNodeName );
   atts = parseAttributes(node);

   cnodes =  node.getChildNodes;
   k = 1; kids = {};
   for i=1:cnodes.getLength;
      cnode_i = cnodes.getChildNodes.item(i-1);
      if strcmp(cnode_i.getNodeName, '#text'); continue; end 
      kids{k} = cnode_i;
      k=k+1;
   end
 



function old_stuff
   s = parseXML( filename );

s.Children(2).Children(4).Children(2).Attributes(1)

s.Children(2).Children(10)


% This function is from http://www.mathworks.de/de/help/matlab/ref/xmlread.html
function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
   tree = xmlread(filename);
catch
   error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems 
% with very deeply nested trees.
try
   theStruct = parseChildNodes(tree);
catch
   error('Unable to parse XML file %s.',filename);
end


% ----- Local function PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);

    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end

% ----- Local function MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
   'Name', char(theNode.getNodeName),       ...
   'Attributes', parseAttributes(theNode),  ...
   'Data', '',                              ...
   'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
   nodeStruct.Data = '';
end

% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end
