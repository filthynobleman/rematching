function N = rmt_resample(M, OutputRes)
%N = RMT_RESAMPLE(M, OutputRes) Preprocess the triangular mesh M and
%resamples it to make it suitable for a remeshing to OutputRes vertices.
%
%   This function applies a resampling strategy to the input triangular
%   mesh M that improves the uniformity of the remeshing when applying the
%   ReMathcing remesh algorithm for an output resolution of OutputRes
%   vertices. The input mesh M must have the following fields:
%   - n the number of vertices of M;
%   - m the number of triangles of M;
%   - VERT a n-by-3 matrix containing in each row the coordinates of the
%   vertices of M;
%   - TRIV a m-by-3 matrix containing in each row the indices of the
%   vertices forming the triangles of M.
%
%
%
%Author:        Filippo Maggioli 
%               'La Sapienza' Department of Computer Science
%EMail:         maggioli.filippo@gmail.com maggioli@di.uniroma1.it
%Last Revision: 27 October 2023

    [N.VERT, N.TRIV] = mex_resample(double(M.VERT), int32(M.TRIV), int32(OutputRes));
    N.n = size(N.VERT, 1);
    N.m = size(N.TRIV, 1);

end

