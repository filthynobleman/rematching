function U = rmt_wmap(Src, Rem)
%U = RMT_WMAP(Src, Rem) Computes the weightmap that transfer functions
%from the remesh to the source full shape.
%
%   This function computes the sparse mapping that brings scalar
%   function defined on the remeshed shape Rem to the original full
%   resolution triangular mesh Src. Each input mesh M must have the 
%   following fields:
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
%Last Revision: 6 November 2023

    [I, J, V] = mex_wmap(double(Rem.VERT), int32(Rem.TRIV), double(Src.VERT));
    U = sparse(I, J, V, Src.n, Rem.n);
    
end