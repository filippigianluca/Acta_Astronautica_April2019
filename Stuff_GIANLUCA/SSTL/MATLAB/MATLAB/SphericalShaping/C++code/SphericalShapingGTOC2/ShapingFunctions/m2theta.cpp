#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <iostream>

#include "loaddata.h"

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

// Implementation of M2theta
//-------------------------------------------------------------------------//

double M2theta(double M, double e)

{
    using namespace std;

    double err = 1.0;
    double toll = 1.e-7;
    double E = 0.0;
    double ratio = 0.0;
    double tan_theta2 = 0.0;
    double theta = 0.0;

    if (M < tudat::basic_mathematics::mathematical_constants::PI)
       { E = M + e/2.0;  }
    else
       { E = M - e/2.0;  }


    while (err > toll)
    {
        ratio = (M - E + e*sin(E) ) / (1.0 - e*cos(E));
        E = E + ratio;
        err = abs(ratio);
    }

    tan_theta2 = sqrt((1.0+e) / (1.0-e))*tan(E/2.0);
    theta = 2.0*atan(tan_theta2);

    if (theta < 0.0)
      {  theta = theta + 2.0*tudat::basic_mathematics::mathematical_constants::PI;  }

    return theta;

}
