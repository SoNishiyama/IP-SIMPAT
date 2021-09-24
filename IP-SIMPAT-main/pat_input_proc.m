function [indlist, markerlist, genmat, nummar, numind, chr, bp, len] = pat_input_proc(gen, pos, chrlen)
% process all the input data into the form for optimization

indlist = gen.textdata(1,3:end);
markerlist = gen.textdata(2:end,1);
genmat = gen.data(:,2:end);
[nummar, numind] = size(genmat);

chr = pos.data(:,1);
bp = pos.data(:,2);
len = chrlen.data(:,2);