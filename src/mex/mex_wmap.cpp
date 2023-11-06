/**
 * @file        mex_wmap.cpp
 * 
 * @brief       MEX interface for the weightmap function.
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

        matlab::data::TypedArray<double> MVSrc = std::move(inputs[0]);
        matlab::data::TypedArray<int> MFSrc = std::move(inputs[1]);
        matlab::data::TypedArray<double> MVTrg = std::move(inputs[2]);

        int* UI;
        int* UJ;
        double* UV;
        size_t NNZ;
        rmt::WeightMap(MVSrc.release().get(), MVSrc.getDimensions()[0], 
                       MFSrc.release().get(), MFSrc.getDimensions()[0],
                       MVTrg.release().get(), MVTrg.getDimensions()[0],
                       &UI, &UJ, &UV, &NNZ);

        matlab::data::ArrayFactory factory;
        matlab::data::ArrayDimensions Dims = { NNZ, 1 };
        outputs[0] = factory.createArray(Dims, UI, UI + NNZ);
        outputs[1] = factory.createArray(Dims, UJ, UJ + NNZ);
        outputs[2] = factory.createArray(Dims, UV, UV + NNZ);
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

        // Third input must be a 3 columns double matrix
        if (inputs[2].getType() != matlab::data::ArrayType::DOUBLE ||
            inputs[2].getDimensions()[1] != 3)
        {
            matlabPtr->feval(u"error", 
                             0, 
                             std::vector<matlab::data::Array>({ factory.createScalar("Third input must be a N-by-3 double matrix.") }));
        }

        // Exactly three outputs
        if (outputs.size() != 3)
        {
            matlabPtr->feval(u"error", 
                             0, 
                             std::vector<matlab::data::Array>({ factory.createScalar("Function returns exactly three outputs.") }));
        }
    }
};