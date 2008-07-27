% Demo algorithms $Id$

algs = {'GREIT_NOSER_ndiff', ...
        'GREIT_Sheffield_backproj' };

alg = algs{1};

k=1;
[img,map] = feval(alg, v(k).vh, v(k).vi );
