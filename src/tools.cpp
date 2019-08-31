#include "tools.h"

using Eigen::VectorXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
    /** Calculate the RMSE */
    VectorXd RMSE = VectorXd(4);
    RMSE << 0,0,0,0;

    VectorXd err = VectorXd(4);
    err << 0,0,0,0;

    // Error notification for different lengths of estimations and ground truths
    if(estimations.size() != ground_truth.size() || estimations.size() < 1)
    {
        cout<<"Cannot compute RMSE. Invalid input size" <<endl;
		return RMSE;
    }

    for(unsigned int i=0;i<estimations.size(); i++)
    {
        err = estimations[i] - ground_truth[i];
        err = err.array() * err.array();
        RMSE = RMSE + err;
    }
    RMSE = RMSE / estimations.size();
    RMSE = RMSE.array().sqrt();
    return RMSE;
}
