/**
 * @file        batch.cpp
 * 
 * @brief       Application for batched remeshing.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-07-31
 */
#include <rmt/rmt.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <chrono>
#include <filesystem>



void StartTimer();
double StopTimer();

struct rmtArgs
{
    std::string InDir;
    std::string OutDir;
    std::string OutCSV;
    bool AllMeshes;
    std::vector<std::string> InMeshes;
    std::vector<std::string> OutMeshes;
    std::vector<std::string> WMaps;

    // int NumSamples;
    std::vector<int> NumSamples;
    // double Resolution;
    std::vector<double> Resolution;
    bool FixedSize;
    bool Resampling;
    bool Evaluate;
};

struct RMTime
{
    double Resampling;
    double Boundary;
    double VoronoiFPS;
    double FlatUnion;
    double Reconstruction;
    double WMap;
    double Total;
};

struct MeshStats
{
    int NVerts;
    int NEdges;
    int NTris;
    int RMSize;
};


rmtArgs ParseArgs(int argc, const char* const argv[]);
std::string OutputName(const rmtArgs& Args, int Run, const std::string& Name);
void SanityCheck(rmtArgs& Args);
void Usage(const std::string& Prog, bool IsError = false);




int main(int argc, const char* const argv[])
{
    // Parse arguments
    rmtArgs Args;
    try
    {
        Args = ParseArgs(argc, argv);
        SanityCheck(Args);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return -1;
    }

    int nMeshes = Args.InMeshes.size();
    int nReps = Args.Resolution.size();
    if (Args.FixedSize)
        nReps = Args.NumSamples.size();
    std::cout << "Remeshing " << nMeshes << " meshes from directory " << Args.InDir << std::endl;
    std::cout << "Output meshes will be saved to " << Args.OutDir << std::endl;
    if (Args.Resampling)
        std::cout << "Mehses will be resampled." << std::endl;
    else
        std::cout << "Meshes will not be resampled." << std::endl;
    if (Args.FixedSize)
    {
        std::cout << "The output meshes will have sizes: ";
        for (auto it = Args.NumSamples.begin(); it != Args.NumSamples.end(); it++)
        {
            std::cout << *it;
            if (it != Args.NumSamples.end() - 1)
                std::cout << ", ";
        }
    }
    else
    {
        std::cout << "The output meshes will have resolutions: ";
        for (auto it = Args.Resolution.begin(); it != Args.Resolution.end(); it++)
        {
            std::cout << ((*it) * 100) << "%";
            if (it != Args.Resolution.end() - 1)
                std::cout << ", ";
        }
    }
    std::cout << std::endl;


    // Execution time and metrics
    std::vector<RMTime> Times;
    Times.resize(nMeshes * nReps);
    std::vector<MeshStats> Stats;
    Stats.resize(nMeshes * nReps);
    std::vector<rmt::EvaluationMetrics> Metrics;
    Metrics.resize(nMeshes * nReps);
    std::vector<bool> Failures;
    Failures.resize(nMeshes * nReps, false);

    // For each mesh, apply the remeshing
    for (int i = 0; i < nMeshes; ++i)
    {
        Eigen::MatrixXd Vp;
        Eigen::MatrixXi Fp;

        if (!rmt::LoadMesh(Args.InMeshes[i], Vp, Fp))
        {
            std::cerr << "Cannot load mesh " << Args.InMeshes[i] << std::endl;
            for (int j = 0; j < nReps; ++j)
                Failures[i * nMeshes + j] = true;
            continue;
        }
        for (int j = 0; j < nReps; ++j)
        {
            int RunIdx = i * nReps + j;
            rmt::Mesh Mesh(Vp, Fp);

            int NSamples = Args.NumSamples[j];
            if (!Args.FixedSize)
                NSamples = Args.Resolution[j] * Mesh.NumVertices();

            std::cout << "Remeshing " << Args.InMeshes[i] << " to " << NSamples << " vertices... " << std::endl;

            int nVertsOrig = Mesh.NumVertices();
            if (Args.Resampling)
            {
                StartTimer();
                Mesh.Resample(NSamples);
                Times[RunIdx].Resampling = StopTimer();
            }

            StartTimer();
            Mesh.ComputeEdgesAndBoundaries();
            Times[RunIdx].Boundary = StopTimer();

            Stats[RunIdx] = { Mesh.NumVertices(), Mesh.NumEdges(), Mesh.NumTriangles(), NSamples };

            StartTimer();
            rmt::VoronoiPartitioning VPart(Mesh);
            while (VPart.NumSamples() < NSamples)
                VPart.AddSample(VPart.FarthestVertex());
            Times[RunIdx].VoronoiFPS = StopTimer();

            StartTimer();
            rmt::FlatUnion FU(Mesh, VPart);
            do
            {
                FU.DetermineRegions();
                FU.ComputeTopologies();
            } while (!FU.FixIssues());
            Times[RunIdx].FlatUnion = StopTimer();
            

            StartTimer();
            Eigen::MatrixXd VV;
            Eigen::MatrixXi FF;
            rmt::MeshFromVoronoi(Mesh.GetVertices(), Mesh.GetTriangles(), VPart, VV, FF);
            Times[RunIdx].Reconstruction = StopTimer();

            Times[RunIdx].Total = Times[RunIdx].Boundary + 
                                  Times[RunIdx].VoronoiFPS +
                                  Times[RunIdx].FlatUnion + 
                                  Times[RunIdx].Reconstruction + 
                                  Times[RunIdx].Resampling;
            std::cout << "\tElapsed time is " << Times[RunIdx].Total << " seconds." << std::endl;

            if (FF.rows() == 0)
            {
                std::cerr << "\tUnable to generate any face." << std::endl;
                Failures[RunIdx] = true;
                continue;
            }

            std::string OutMesh = OutputName(Args, j, Args.OutMeshes[i]);
            if (!rmt::ExportMesh(OutMesh, VV, FF))
                std::cerr << "\tCannot output the mesh to " << OutMesh << '.' << std::endl;
            else
                std::cout << "\tOutput mesh written to " << OutMesh << '.' << std::endl;

            std::string OutWMap = OutputName(Args, j, Args.WMaps[i]);
            StartTimer();
            auto W = rmt::WeightMap(Mesh.GetVertices(), VV, FF, nVertsOrig);
            Times[RunIdx].WMap = StopTimer();
            if (!rmt::ExportWeightmap(OutWMap, W))
                std::cerr << "\tCannot save the weightmap to " << OutWMap;
            else
                std::cout << "\tWeightmap written to " << OutWMap << '.' << std::endl;

            if (Args.Evaluate)
            {
                Mesh.RescaleInsideUnitBox();
                rmt::RescaleInsideUnitBox(VV);
                Metrics[RunIdx] = rmt::Evaluate(Mesh.GetVertices(), Fp, VV, FF, nVertsOrig);
            }
        }
    }


    // Output data
    std::ofstream Out;
    Out.open(Args.OutCSV, std::ios::out);
    if (!Out.is_open())
    {
        std::cerr << "Cannot open file " << Args.OutCSV << " for writing." << std::endl;
        return -1;
    }

    // Write header
    Out << "Mesh,";
    Out << "NVerts,NEdges,NTris,OutResolution,Success,";
    Out << "Time,Resample,Boundary,VoronoiFPS,FlatUnion,Reconstruct,WMap";
    if (Args.Evaluate)
    {
        Out << ",Hausdorff,Chamfer,";
        Out << "MinArea,MaxArea,AvgArea,StdArea,";
        Out << "MinQuality,MaxQuality,AvgQuality,StdQuality";
    }
    Out << '\n';

    for (int i = 0; i < nMeshes; ++i)
    {
        std::string MName = std::filesystem::path(Args.InMeshes[i]).filename().string();
        MName = MName.substr(0, MName.rfind('.'));
        for (int j = 0; j < nReps; ++j)
        {
            int RunIdx = i * nReps + j;

            Out << MName << ',';
            Out << Stats[RunIdx].NVerts << ',' << Stats[RunIdx].NEdges << ',' << Stats[RunIdx].NTris << ',' << Stats[RunIdx].RMSize << ',';
            Out << !Failures[RunIdx] << ',';
            Out << Times[RunIdx].Total << ',';
            Out << Times[RunIdx].Resampling << ',';
            Out << Times[RunIdx].Boundary << ',';
            Out << Times[RunIdx].VoronoiFPS << ',';
            Out << Times[RunIdx].FlatUnion << ',';
            Out << Times[RunIdx].Reconstruction << ',';
            Out << Times[RunIdx].WMap;
            if (Args.Evaluate)
            {
                Out << ',';
                Out << Metrics[RunIdx].Hausdorff << ',' << Metrics[RunIdx].Chamfer << ',';
                Out << Metrics[RunIdx].MinArea << ',';
                Out << Metrics[RunIdx].MaxArea << ',';
                Out << Metrics[RunIdx].AvgArea << ',';
                Out << Metrics[RunIdx].StdArea << ',';
                Out << Metrics[RunIdx].MinQuality << ',';
                Out << Metrics[RunIdx].MaxQuality << ',';
                Out << Metrics[RunIdx].AvgQuality << ',';
                Out << Metrics[RunIdx].StdQuality;
            }
            Out << '\n';
        }
    }
    

    Out.close();

    
    return 0;
}


































std::chrono::system_clock::time_point Start;
void StartTimer()
{
    Start = std::chrono::system_clock::now();
}

double StopTimer()
{
    std::chrono::system_clock::time_point End;
    End = std::chrono::system_clock::now();
    std::chrono::system_clock::duration ETA;
    ETA = End - Start;
    size_t ms;
    ms = std::chrono::duration_cast<std::chrono::milliseconds>(ETA).count();
    return ms * 1.0e-3;
}

void SanityCheck(rmtArgs& Args)
{
    static const std::vector<std::string> Exts = {
        ".obj",
        ".off",
        ".ply"
    };

    // Check existence of input directory
    if (!std::filesystem::exists(Args.InDir))
    {
        std::cerr << "Input directory " << Args.InDir << " does not exist." << std::endl;
        exit(-1);
    }
    if (!std::filesystem::is_directory(Args.InDir))
    {
        std::cerr << Args.InDir << " is not a directory." << std::endl;
        exit(-1);
    }

    // Check output directory
    if (!std::filesystem::exists(Args.OutDir))
        std::filesystem::create_directories(Args.OutDir);
    else
    {
        if (!std::filesystem::is_directory(Args.OutDir))
        {
            std::cerr << Args.OutDir << " is not a directory." << std::endl;
            exit(-1);
        }
    }

    // Must use all the meshes in the input directory?
    if (Args.AllMeshes)
    {
        for (const auto & entry : std::filesystem::directory_iterator(Args.InDir))
        {
            int ELen = entry.path().string().length();
            bool Enlist = false;
            for (auto Ext : Exts)
            {
                int EExt = entry.path().string().rfind(Ext);
                if (EExt == std::string::npos || EExt != ELen - 4)
                    continue;
                Enlist = true;
            }
            if (!Enlist)
                continue;
            Args.InMeshes.emplace_back(entry.path().filename().string());
        }
    }

    // Sort num samples/resolution
    if (Args.FixedSize)
        std::sort(Args.NumSamples.begin(), Args.NumSamples.end());
    else
        std::sort(Args.Resolution.begin(), Args.Resolution.end());

    // Join filenames and create output paths
    std::filesystem::path IMD, OMD, OWD;
    IMD = std::filesystem::path(Args.InDir);
    OMD = std::filesystem::path(Args.OutDir) / std::filesystem::path("meshes");
    OWD = std::filesystem::path(Args.OutDir) / std::filesystem::path("wmaps");
    Args.OutMeshes.resize(Args.InMeshes.size());
    Args.WMaps.resize(Args.InMeshes.size());
    for (int i = 0; i < Args.InMeshes.size(); ++i)
    {
        std::string Name = Args.InMeshes[i];
        Args.InMeshes[i] = (IMD / std::filesystem::path(Name)).string();
        Args.OutMeshes[i] = (OMD / std::filesystem::path(Name)).string();
        Name = Name.substr(0, Name.rfind('.')) + ".mat";
        Args.WMaps[i] = (OWD / std::filesystem::path(Name)).string();
    }

    // Create output directories
    if (!std::filesystem::exists(OMD.string()))
        std::filesystem::create_directories(OMD.string());
    if (!std::filesystem::exists(OWD.string()))
        std::filesystem::create_directories(OWD.string());

    // Name of the CSV data file
    Args.OutCSV = (std::filesystem::path(Args.OutDir) / std::filesystem::path("batch.csv")).string();
}

rmtArgs ParseFromFile(const std::string& Filename)
{
    std::ifstream Stream;
    Stream.open(Filename, std::ios::in);
    if (!Stream.is_open())
    {
        std::cerr << "Cannot open file " << Filename << " for reading." << std::endl;
        exit(-1);
    }

    nlohmann::json j;
    try
    {
        Stream >> j;
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        exit(-1);
    }
    
    
    Stream.close();

    rmtArgs Args;
    // Args.NumSamples = 0;
    // Args.Resolution = 0.0;
    Args.FixedSize = true;
    Args.Resampling = false;
    Args.Evaluate = false;

    std::vector<std::string> Attrs = {
        "input_dir",
        "output_dir",
        "meshes"
    };
    for (const auto& attr : Attrs)
    {
        if (!j.contains(attr))
        {
            std::cerr << Filename << " does not contain attribute \"" << attr << "\"." << std::endl;
            exit(-1);
        }
    }

    if (!j["input_dir"].is_string())
    {
        std::cerr << "Attribute \"input_dir\" in " << Filename << " is not a string." << std::endl;
        exit(-1);
    }
    if (!j["output_dir"].is_string())
    {
        std::cerr << "Attribute \"input_dir\" in " << Filename << " is not a string." << std::endl;
        exit(-1);
    }

    Args.InDir = j["input_dir"];
    Args.OutDir = j["output_dir"];
    Args.AllMeshes = false;
    if (j["meshes"].is_array())
        Args.InMeshes = j["meshes"];
    else if (j["meshes"].is_string())
    {
        if (std::string(j["meshes"]) == "*")
            Args.AllMeshes = true;
        else
            Args.InMeshes.emplace_back(j["meshes"]);
    }
    else
    {
        std::cerr << "Attribute \"meshes\" in " << Filename << " is neither a string nor an array." << std::endl;
        exit(-1);
    }
    
    if (j.contains("resampling"))
    {
        if (!j["resampling"].is_boolean())
        {
            std::cerr << Filename << " contains attribute \"resampling\", but it is not a boolean." << std::endl;
            exit(-1);
        }
        Args.Resampling = j["resampling"];
    }
    if (j.contains("evaluate"))
    {
        if (!j["evaluate"].is_boolean())
        {
            std::cerr << Filename << " contains attribute \"evaluate\", but it is not a boolean." << std::endl;
            exit(-1);
        }
        Args.Evaluate = j["evaluate"];
    }

    if (j.contains("fixed_size"))
    {
        if (!j["fixed_size"].is_boolean())
        {
            std::cerr << Filename << " contains attribute \"fixed_size\", but it is not a boolean." << std::endl;
            exit(-1);
        }
        Args.FixedSize = j["fixed_size"];
    }
    else
    {
        if (j.contains("num_samples") && !j.contains("resolution"))
            Args.FixedSize = true;
        else if (j.contains("resolution") && !j.contains("num_samples"))
            Args.FixedSize = false;
        else if (j.contains("resolution") && j.contains("num_samples"))
        {
            std::cerr << Filename << "contains both attributes \"resolution\" and \"num_samples\", but the ";
            std::cerr << "attribute \"fixed_size\" is not set to determine which mode must be used." << std::endl;
            exit(-1);
        }
        else
        {
            std::cerr << Filename << " does not define neither attribute \"resolution\" nor \"num_samples." << std::endl;
            exit(-1);
        }
    }
    if (Args.FixedSize)
    {
        if (!j.contains("num_samples"))
        {
            std::cerr << Filename << " has attribute \"fixed_size\" set to true, but does not contain attribute \"num_samples\"." << std::endl;
            exit(-1);
        }
        if (!j["num_samples"].is_number_integer() && !j["num_samples"].is_array())
        {
            std::cerr << Filename << " contains attribute \"num_samples\", but it is not an integer number neither an array." << std::endl;
            exit(-1);
        }
        if (j["num_samples"].is_number_integer())
            Args.NumSamples.emplace_back(j["num_samples"]);
        else
            Args.NumSamples = j["num_samples"].template get<std::vector<int>>();
    }
    else
    {
        if (!j.contains("resolution"))
        {
            std::cerr << Filename << " has attribute \"fixed_size\" set to false, but does not contain attribute \"resolution\"." << std::endl;
            exit(-1);
        }
        if (!j["resolution"].is_number_float() && !j["resolution"].is_array())
        {
            std::cerr << Filename << " contains attribute \"resolution\", but it is not a floating point number neither an array." << std::endl;
            exit(-1);
        }
        if (j["resolution"].is_number_float())
            Args.Resolution.emplace_back(j["resolution"]);
        else
            Args.Resolution = j["resolution"].template get<std::vector<double>>();
        for (int i = 0; i < Args.Resolution.size(); ++i)
        {
            if (Args.Resolution[i] <= 0.0 && Args.Resolution[i] > 1.0)
            {
                std::cerr << Filename << " contains a \"resolution\" attribute with value (" << Args.Resolution[i] << ") , which is outside the interval (0, 1]." << std::endl;
                exit(-1);
            }
        }
    }

    return Args;
}

rmtArgs ParseArgs(int argc, const char* const argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        std::string argvi(argv[i]);
        if (argvi == "-h" || argvi == "--help")
        {
            Usage(argv[0]);
            exit(0);
        }
    }

    if (argc < 2)
    {
        std::cerr << "No input file given." << std::endl;
        Usage(argv[0], true);
    }

    return ParseFromFile(argv[1]);
}


void Usage(const std::string& Prog, bool IsError)
{
    std::ostream* _out = &std::cout;
    if (IsError)
        _out = &std::cerr;
    std::ostream& out = *_out;

    out << std::endl;
    out << Prog << " usage:" << std::endl;
    out << std::endl;
    out << "\t" << Prog << " config_file" << std::endl;
    out << "\t" << Prog << " -h|--help" << std::endl;
    out << std::endl;
    out << "Arguments details:" << std::endl;
    out << "\t- config_file is a JSON file with the settings for the run;" << std::endl;
    out << "\t- -h|--help prints this message." << std::endl;

    if (IsError)
        exit(-1);
}


std::string OutputName(const rmtArgs& Args, int Run, const std::string& Name)
{
    // If only a single run must be performed, output the identity
    if (Args.FixedSize && Args.NumSamples.size() == 1)
        return Name;
    else if (!Args.FixedSize && Args.Resolution.size() == 1)
        return Name;

    // Otherwise, get a proper directory name
    std::stringstream ss;
    if (Args.FixedSize)
    {
        int MaxPrecision = 1;
        for (int i = 0; i < Args.NumSamples.size(); ++i)
            MaxPrecision = std::max(MaxPrecision, (int)std::floor(std::log10(Args.NumSamples[i])));
        MaxPrecision += 1;
        ss << "ns-" << std::setfill('0') << std::setw(MaxPrecision) << Args.NumSamples[Run];
    }
    else
    {
        int MaxPrecision = -1;
        for (int i = 0; i < Args.Resolution.size(); ++i)
            MaxPrecision = std::min(MaxPrecision, (int)std::ceil(std::log10(Args.Resolution[i])));
        MaxPrecision = -(MaxPrecision - 1);
        ss.precision(MaxPrecision);
        ss << "res-" << Args.Resolution[Run];
    }
    std::string Dir = ss.str();

    // Decompose the name into directory and basename
    std::filesystem::path NPath(Name);
    std::filesystem::path Parent = NPath.parent_path();
    std::filesystem::path BName = NPath.filename();
    std::filesystem::path DirPath(Dir);

    // Endure output path exists
    if (!std::filesystem::exists(Parent / DirPath))
        std::filesystem::create_directories(Parent / DirPath);

    return (Parent / DirPath / BName).string();
}