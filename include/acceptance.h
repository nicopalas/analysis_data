
#include "constants.h"

static bool loadAcceptanceCSV(const std::string& fname, Vec2D& dOmega)
{
    dOmega.assign(nbins_beam_fine,
        std::vector<double>(nbins_det_fine, 0.0));

    std::ifstream fin(fname);
    if(!fin.is_open()){
        std::cerr << "No se pudo abrir " << fname << "\n";
        return false;
    }

    for(int j = 0; j < nbins_beam_fine; ++j)
        for(int i = 0; i < nbins_det_fine; ++i)
            dOmega[j][i] = 0.0;

    std::string line;
    std::getline(fin, line);

    while(std::getline(fin, line)){
        if(line.empty()) continue;

        std::stringstream ss(line);
        double cb = 0.0, cd = 0.0, omega = 0.0;
        long long counts = 0;
        char comma1, comma2, comma3;

        ss >> cb >> comma1 >> cd >> comma2 >> counts >> comma3 >> omega;
        if(!ss) continue;

        int j = static_cast<int>(cb / dcos_beam_fine);
        int i = static_cast<int>(cd / dcos_det_fine);

        if(j == nbins_beam_fine) j = nbins_beam_fine - 1;
        if(i == nbins_det_fine)  i = nbins_det_fine  - 1;

        if(j < 0 || j >= nbins_beam_fine) continue;
        if(i < 0 || i >= nbins_det_fine)  continue;

        dOmega[j][i] = omega;
    }

    return true;
}

static Vec2D rebin(const Vec2D& fine)
{
    Vec2D coarse(nbins_beam, std::vector<double>(nbins_det, 0.0));
    const int ratio_beam = nbins_beam_fine / nbins_beam;
    const int ratio_det  = nbins_det_fine  / nbins_det;

    for(int ib = 0; ib < nbins_beam; ++ib){
        int jb0 = ib * ratio_beam;
        int jb1 = (ib + 1) * ratio_beam;

        for(int id = 0; id < nbins_det; ++id){
            int jd0 = id * ratio_det;
            int jd1 = (id + 1) * ratio_det;

            double sum = 0.0;
            for(int jf = jb0; jf < jb1; ++jf){
                for(int kf = jd0; kf < jd1; ++kf){
                    sum += fine[jf][kf];
                }
            }
            coarse[ib][id] = sum;
        }
    }
    return coarse;
}