function M = histogram_match(I, J)
% A deterministic histogram matching algorithm

% Cumulative histogram of the first image
Ix = linspace(0,max(I(:)),64);
Jx = linspace(0,max(J(:)),64);
I_ch = cumsum(histc(I(I>0),Ix));
J_ch = cumsum(histc(J(J>0),Jx));

I_ch = I_ch / max(I_ch);
J_ch = J_ch / max(J_ch);

% We need to drop repeated values in the cumulative histogram
[I_ch_uniq I_ch_idx] = unique(I_ch);

% These are the intensities in I corresponding to bins in J
Q = interp1(I_ch_uniq, Ix(I_ch_idx), J_ch, 'cubic', 0);
%plot(Jx, Q);

% Now interpolate the intensities in J
J_remap = interp1(Jx, Q, J(J > 0), 'linear', 0);
M = J;
M(J > 0) = J_remap;