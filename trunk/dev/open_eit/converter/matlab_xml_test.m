function obj= matlab_xml_test( filename)
   if nargin==0; filename = 'sample.xml'; end
   try
      tree = xmlread(filename);
   catch
      error('Failed to read XML file %s.',filename);
   end

   obj = print_recurse( tree, 0);

function obj = print_recurse( node, level);
   if level>5; obj.name ='empty'; return;end
   
   obj = struct();

   [name, atts,  kids] = non_text_nodes( node);
   obj.name = name;
   % deftault fields
   obj.comment = '';
   switch name % this is a hack
      case 'electrode'; obj.contact_impedance = '???';
   end

   for i=1:level; fprintf('>'); end
   fprintf('[ %s: ',name);
   for i=1:length(atts);
      fprintf('%s="%s" ',atts(i).Name, atts(i).Value);
      obj.(atts(i).Name) = atts(i).Value;
   end
   fprintf(']\n');
   for i=1:length(kids);
     obj_i = print_recurse( kids{i}, level+1);
     on = obj_i.name;
     if strcmp(on, 'empty'); continue; end
     if strcmp(on, 'user_data'); continue; end

     try   ll = length(obj.( on ) );
     catch ll = 0;
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
