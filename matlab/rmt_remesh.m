function N = rmt_remesh(M, NSamples)
%N = RMT_REMESH(M, NSamples) Remeshes the input triangular mesh M with the
%ReMatching algorithm to NSamples vertices.
%
%   This function applies the ReMatching remeshing algorithm to the input
%   triangular mesh M, outputting a new triangular mesh N that has NSamples
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

    
    [N.VERT, N.TRIV] = mex_remesh(double(M.VERT), int32(M.TRIV), int32(NSamples));
    N.n = size(N.VERT, 1);
    N.m = size(N.TRIV, 1);


end

