/**
 * @file        mex_remesh.cpp
 * 
 * @brief       MEX interface for the remeshing function.
 * 
 * @author      Filippo Maggioli\n
 *              (maggioli@di.uniroma1.it, maggioli.filippo@gmail.com)\n
 *              Sapienza, University of Rome - Department of Computer Science
 * 
 * @date        2023-10-26
 */
#include <mex.hpp>
#include <mexAdapter.hpp>
#include <rmt/mex/rmt.hpp>


class MexFunction : public matlab::mex::Function
{
public:
    void operator()(matlab::mex::ArgumentList outputs,
                    matlab::mex::ArgumentList inputs)
    {
        ValidateArgs(outputs, inputs);

        matlab::data::TypedArray<double> MVin = std::move(inputs[0]);
        matlab::data::TypedArray<int> MFin = std::move(inputs[1]);

        double* Vout;
        int* Fout;
        size_t NVout = 0;
        size_t NFout = 0;
        rmt::Remesh(MVin.release().get(), MVin.getDimensions()[0], 
                    MFin.release().get(), MFin.getDimensions()[0], 
                    inputs[2][0], 
                    &Vout,
                    &NVout, 
                    &Fout,
                    &NFout);

        matlab::data::ArrayFactory factory;
        matlab::data::ArrayDimensions VDims = { NVout, 3 };
        matlab::data::ArrayDimensions FDims = { NFout, 3 };
        outputs[0] = factory.createArray(VDims, Vout, Vout + 3 * NVout);
        outputs[1] = factory.createArray(FDims, Fout, Fout + 3 * NFout);
    }


    void ValidateArgs(matlab::mex::ArgumentList outputs,
                      matlab::mex::ArgumentList inputs)
    {
        // Get pointer to engine
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

        // Get array factory
        matlab::data::ArrayFactory factory;

        // First input must be a 3 columns double matrix
        if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[0].getDimensions()[1] != 3)
        {
            matlabPtr->feval(u"error", 
                             0, 
                             std::vector<matlab::data::Array>({ factory.createScalar("First input must be a N-by-3 double matrix.") }));
        }

        // Second input must be a 3 columns int32 matrix
        if (inputs[1].getType() != matlab::data::ArrayType::INT32 ||
            inputs[1].getDimensions()[1] != 3)
        {
            matlabPtr->feval(u"error", 
                             0, 
                             std::vector<matlab::data::Array>({ factory.createScalar("Second input must be a N-by-3 int32 matrix.") }));
        }

        // Third input must be a scalar int32
        if (inputs[2].getType() != matlab::data::ArrayType::INT32 ||
            inputs[2].getNumberOfElements() != 1)
        {
            matlabPtr->feval(u"error", 
                             0, 
                             std::vector<matlab::data::Array>({ factory.createScalar("This input must be a scalar int32.") }));
        }

        // Exactly two outputs
        if (outputs.size() != 2)
        {
            matlabPtr->feval(u"error", 
                             0, 
                             std::vector<matlab::data::Array>({ factory.createScalar("Function returns exactly two outputs.") }));
        }
    }
};