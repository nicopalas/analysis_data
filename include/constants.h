#pragma once

constexpr int    nbins_beam      = 10;
constexpr int    nbins_det      = 10;
constexpr int    nbins_beam_fine = 100;
constexpr int    nbins_det_fine  = 20;
constexpr double dcos_beam       = 1.0 / nbins_beam;
constexpr double dcos_det        = 1.0 / nbins_det;
constexpr double dcos_beam_fine = 1.0 / nbins_beam_fine;
constexpr double dcos_det_fine   = 1.0 / nbins_det_fine;
