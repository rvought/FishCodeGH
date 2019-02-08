function out = overture(in)
% Calculates the overlap in movements for all pairs of fish


diffhist = in(1).realhist - in(2).realhist;

overlapfishone = 1 - (sum(diffhist(diffhist > 0)) / sum(sum(in(1).realhist)));
overlapfishtwo = 1 - (sum(diffhist(diffhist < 0)) / sum(sum(in(2).realhist)));

