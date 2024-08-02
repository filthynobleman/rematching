#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <rmt/rmt.hpp>

namespace py = pybind11;

class PyRMTMesh
{
public:
    Eigen::MatrixXd Vertices;
    Eigen::MatrixXi Triangles;

    PyRMTMesh(Eigen::Ref<Eigen::MatrixXd> Verts,
              Eigen::Ref<Eigen::MatrixXi> Tris)
        : Vertices(Verts), Triangles(Tris) { }

    PyRMTMesh Remesh(int NumSamples)
    {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        rmt::Remesh(Vertices, Triangles, NumSamples, V, F);

        return PyRMTMesh(V, F);
    }

    void MakeManifold()
    {
        rmt::MakeManifold(Vertices, Triangles);
    }

    void CleanUp()
    {
        rmt::CleanUp(Vertices, Triangles);
    }

    Eigen::SparseMatrix<double> BarycMap(Eigen::Ref<Eigen::MatrixXd> Points)
    {
        return rmt::WeightMap(Points, Vertices, Triangles);
    }
};


PYBIND11_MODULE(PyRMT, m)
{
    m.doc() = "Python binding to the remeshing operations in ReMatching.";

    py::class_<PyRMTMesh>(m, "RMTMesh")
        .def(py::init<Eigen::Ref<Eigen::MatrixXd>, Eigen::Ref<Eigen::MatrixXi>>())
        .def_readwrite("vertices", &PyRMTMesh::Vertices)
        .def_readwrite("triangles", &PyRMTMesh::Triangles)
        .def("remesh", &PyRMTMesh::Remesh)
        .def("make_manifold", &PyRMTMesh::MakeManifold)
        .def("clean_up", &PyRMTMesh::CleanUp)
        .def("baryc_map", &PyRMTMesh::BarycMap);
}