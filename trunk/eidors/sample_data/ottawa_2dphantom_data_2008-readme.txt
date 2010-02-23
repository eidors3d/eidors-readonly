Data collected Feb 26, 2008.

Data was averaged since the images are static (no movement during measurements).

matlab script:
  % Load measured data: Goe-MF II
  for i=20:28;
    vv=eval(sprintf('eidors_readdata(''phantom_data/data/999999100%d.get'');',i));
    eval(sprintf('phantomdata_%d=mean(vv,2);',i));
  end

  % electrode positions measured from photos
  % order is [ x y ] pair of columns per phantomdata_%d, in order
  % first three rows give a normalization basis as three points from the graph paper grid
  %  taped to the floor of the phantom, when used to form an x- and y-axis basis, they will correct
  %  for any camera angle skew
  % a single graph paper unit is 2.5mm
  % the remaining rows are the electrode locations 1 --> 16, then gnd node location
  %   and centre and outer location of target objects where appropriate
   d = [698,1496,1043,1193,1448,672,1396,738,1512,531,1492,553,1440,592,1426,507,1117,663;
  736,1498,1065,1192,1479,671,1426,736,1533,531,1516,553,1470,592,1455,506,1149,663;
  695,1531,1043,1212,1450,700,1398,768,1513,551,1492,577,1440,622,1427,535,1117,695;
  794,2094,935,1872,1465,1392,1407,1505,1448,1260,1402,1425,1424,1538,1533,1497,1320,1537;
  543,1985,783,1817,1251,1328,1232,1405,1328,1175,1248,1340,1236,1441,1352 ,1425,1119,1440;
  373,1776,681,1710,1115,1181,1126,1238,1251,1054,1165,1192,1105,1282,1219,1283,979,1290;
  302,1503,629,1550,1027,968,1076,1002,1216,884,1127,990,1025,1053,1139,1071,902,1033;
  348,1236,647,1397,1024,758,1077,765,1220,707,1130,790,1017,815,1138,842,902,780;
  485,990,715,1251,1102,541,1127,536,1258,538,1175,596,1066,577,1168,616,949,525;
  664,818,815,1155,1231,399,1236,389,1336,427,1269,470,1184,439,1274,479,1076,376;
  955,752,970,1107,1446,326,1427,308,1487,376,1436,408,1378,360,1457,395,1288,292;
  1194,771,1106,1116,1633,319,1592,313,1609,384,1577,408,1562,357,1633,381,1477,300;
  1439,893,1247,1172,1828,384,1765,397,1730,456,1720,490,1755,438,1819,450,1676,410;
  1622,1088,1358,1282,1992,529,1878,580,1808,599,1813,648,1903,616,1963,611,1835,578;
  1715,1353,1412,1433,2074,739,1930,814,1848,766,1861,838,1984,826,2054,808,1921,799;
  1680,1647,1403,1598,2088,965,1939,1063,1850,945,1863,1058,1995,1054,2063,1043,1931,1055;
  1574,1850,1352,1718,2021,1134,1901,1245,1820,1074,1819,1210,1960,1246,2024,1218,1869,1250;
  1335,2047,1233,1840,1871,1310,1771,1437,1707,1206,1708,1364,1819,1426,1901,1407,1728,1447;
  1071,2146,1086,1902,1677,1404,1599,1530,1586,1266,1563,1447,1639,1531,1732,1515,1534,1541;
  943,1461,1151,1356,1515,877,1581,861,1585,825,1564,875,1556,878,1835,894,1473,845;
  0,0,927,1560,1731,1058,0,0,1511,1036,1333,854,0,0,1547,1060,1284,1005;
  0,0,927,1417,1710,1084,0,0,1676,1037,1318,878,0,0,1535,867,1321,1005];
  
  % locations converted to complex coordinates: z = x + i*y
  for j = 1:size(d,2)/2;
    normalize(:,j) = d(1:3,2*j-1) + i*d(1:3,2*j);
    elec(:,j) = d(4:4+16-1,2*j-1) + i*d(4:4+16-1,2*j);
  end
  % normalize has [ origin ; unit x-axis; unit y-axis ] and
  %is to be used to normalize the electrode locations
  % in elec
  unitx = normalize(2,:) - normalize(1,:);
  unity = normalize(3,:) - normalize(1,:);
  elec2 = zeros(size(elec));
  for e = 1:length(unitx)
    elec2(:,e) = elec(:,e) ./ unitx(e) ./ unity(e) * 2.5; % 2.5mm per square on graph paper
    xm = mean(real(elec2(:,e)));
    ym = mean(imag(elec2(:,e)));
    elec2(:,e) = elec2(:,e) - xm - i*ym; % remove mean value to center electrode locations
  end
  
  phantomdata_true_electrode_positions = elec2;
  save ottawa200802phantoma_data phantomdata_*;

The averaged data is stored in variables phantomdata_* where * is:
  % 20 no compression, no target
  % 21      ''       , glass target 
  % 22      ''       , metal target

  % 23 2 pt compression, no target
  % 24      ''       , glass target 
  % 25      ''       , metal target

  % 26 3 pt compression, no target
  % 27      ''       , glass target 
  % 28      ''       , metal target

True electrode locations are stored in variable phantomdata_true_electrode_positions

