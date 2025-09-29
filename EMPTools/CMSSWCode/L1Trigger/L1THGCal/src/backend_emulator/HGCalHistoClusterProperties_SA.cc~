#include "L1Trigger/L1THGCal/interface/backend_emulator/HGCalHistoClusterProperties_SA.h"

#include <cmath>
#include <algorithm>

using namespace std;
using namespace l1thgcfirmware;

HGCalHistoClusterProperties::HGCalHistoClusterProperties(const ClusterAlgoConfig& config) : config_(config) {}

void HGCalHistoClusterProperties::runClusterProperties(
    const l1thgcfirmware::HGCalClusterSAPtrCollection& protoClustersIn,
    const CentroidHelperPtrCollection& readoutFlags,
    HGCalClusterSAPtrCollection& clustersOut) const {
  // Cluster properties
  HGCalClusterSAPtrCollection clusterAccumulation;
  clusterSum(protoClustersIn, readoutFlags, clusterAccumulation, clustersOut);
  calcProperties(clustersOut);
}

// Accumulates/combines inputs cluster objects (each corresponding to one TC belonging to a cluster) into clusters  (one per cluster made up of TCs)
void HGCalHistoClusterProperties::clusterSum(const HGCalClusterSAPtrCollection& protoClusters,
                                             const CentroidHelperPtrCollection& readoutFlags,
                                             HGCalClusterSAPtrCollection& clusterAccumulation,
                                             HGCalClusterSAPtrCollection& clusterSums) const {
  HGCalClusterSAPtrCollections protoClustersPerColumn(config_.cColumns());
  vector<unsigned int> clock(config_.cColumns(), 0);
  for (const auto& protoCluster : protoClusters) {
    auto index = protoCluster->index();
    // Do we have to make a copy of protoCluster here?
    protoClustersPerColumn.at(index).push_back(make_unique<HGCalCluster>(*protoCluster));
  }

  map<unsigned int, HGCalClusterSAPtr> sums;

  for (const auto& flag : readoutFlags) {
    auto accumulator = make_unique<HGCalCluster>(0, 0, true, true);
    const unsigned stepLatency = 23;
    flag->setClock(flag->clock() + stepLatency);

    for (const auto& protoCluster : protoClustersPerColumn.at(flag->index())) {
      if (protoCluster->clock() <= clock.at(flag->index()))
        continue;
      if (protoCluster->clock() > flag->clock())
        continue;
      *accumulator += *protoCluster;
    }

    clock.at(flag->index()) = flag->clock();
    accumulator->setClock(flag->clock());
    accumulator->setIndex(flag->index());
    accumulator->setDataValid(true);

    if (sums.find(flag->clock()) == sums.end()) {
      const unsigned stepLatency = 7;
      auto sum = make_unique<HGCalCluster>(flag->clock() + stepLatency, 0, true, true);
      sums[flag->clock()] = move(sum);
    }

    *(sums.at(flag->clock())) += *accumulator;

    clusterAccumulation.push_back(move(accumulator));
  }

  for (auto& sum : sums) {
    clusterSums.push_back(move(sum.second));
  }
}

// Calculates properties of clusters from accumulated quantities
void HGCalHistoClusterProperties::calcProperties(HGCalClusterSAPtrCollection& clusterSums) const {

  for (auto& c : clusterSums) {

    HGCalCluster_HW& hwCluster = c->hwCluster();
    hwCluster.e = Scales::HGCaltoL1_et(c->e());
    hwCluster.e_em = Scales::HGCaltoL1_et(c->e_em());
    hwCluster.fractionInCE_E = Scales::makeL1EFraction(c->e_em(), c->e());
    hwCluster.fractionInCoreCE_E = Scales::makeL1EFraction(c->e_em_core(), c->e_em());
    hwCluster.fractionInEarlyCE_E = Scales::makeL1EFraction(c->e_h_early(), c->e());
    hwCluster.setGCTBits();
    std::vector<int> layeroutput = showerLengthProperties(c->layerbits());
    c->set_firstLayer(layeroutput[0]);
    c->set_lastLayer(layeroutput[1]);
    c->set_showerLen(layeroutput[2]);
    c->set_coreShowerLen(layeroutput[3]);
    hwCluster.firstLayer = c->firstLayer();
    hwCluster.lastLayer = c->lastLayer();
    hwCluster.showerLength = c->showerLen();
    hwCluster.coreShowerLength = c->coreShowerLen();
    hwCluster.nTC = c->n_tc();

    if (c->n_tc_w() == 0)
      continue;

    hwCluster.w_eta = convertRozToEta( c );
    bool saturatedPhi = false;
    bool nominalPhi = false;
    hwCluster.w_phi = Scales::HGCaltoL1_phi(float(c->wphi())/c->w(), saturatedPhi, nominalPhi);
    hwCluster.w_z = Scales::HGCaltoL1_z( float(c->wz()) / c->w() );

    // Quality flags are placeholders at the moment
    hwCluster.setQualityFlags(Scales::HGCaltoL1_et(c->e_em_core()), Scales::HGCaltoL1_et(c->e_h_early()), c->sat_tc(), c->shapeq(), saturatedPhi, nominalPhi);

    const double sigma_E_scale = 0.008982944302260876 * 4;  // Et scale changed, multiplying existing sigma E scale by 4 to compensate
    hwCluster.sigma_E = sigma_coordinate(c->n_tc_w(), c->w2(), c->w(), sigma_E_scale);

    const double sigma_z_scale = 0.08225179463624954;
    hwCluster.sigma_z = sigma_coordinate(c->w(), c->wz2(), c->wz(), sigma_z_scale);

    const double sigma_phi_scale = 0.907465934753418;
    hwCluster.sigma_phi = sigma_coordinate(c->w(), c->wphi2(), c->wphi(), sigma_phi_scale);

    hwCluster.sigma_eta = convertSigmaRozRozToSigmaEtaEta(c);

    const double sigma_roz_scale = 0.5073223114013672;
    unsigned int sigma_roz = sigma_coordinate(c->w(), c->wroz2(), c->wroz(), sigma_roz_scale);
    // Emulation of a bug in firmware
    // if ( sigma_roz >=256 ) sigma_roz -= 256;
    while (sigma_roz >= 256) sigma_roz -= 256;
    if ( sigma_roz > 127 ) sigma_roz = 127;
    hwCluster.sigma_roz = sigma_roz;
  }
}

unsigned int HGCalHistoClusterProperties::sigma_coordinate(unsigned int w,
                                                            unsigned long int wc2,
                                                            unsigned int wc,
                                                            double scale ) const {
  if ( w == 0 ) return 0;
  unsigned int sigma = round(sqrt( (float(w)*float(wc2) - float(wc) * float(wc))  / ( float(w) * float(w) ) ) * scale);
  return sigma;
}

std::vector<int> HGCalHistoClusterProperties::showerLengthProperties(unsigned long int layerBits) const {
  int counter = 0;
  int firstLayer = 0;
  bool firstLayerFound = false;
  int lastLayer = 0;
  std::vector<int> layerBits_array;

  bitset<34> layerBitsBitset(layerBits);
  for (size_t i = 0; i < layerBitsBitset.size(); ++i) {
      bool bit = layerBitsBitset[34-1-i];
      if ( bit ) {
        if ( !firstLayerFound ) {
          firstLayer = i + 1;
          firstLayerFound = true;
        }
        lastLayer = i+1;
        counter += 1;
      } else {
        layerBits_array.push_back(counter);
        counter = 0;
      }
  }

  int showerLen = lastLayer - firstLayer + 1;
  int coreShowerLen = config_.nTriggerLayers();
  if (!layerBits_array.empty()) {
    coreShowerLen = *std::max_element(layerBits_array.begin(), layerBits_array.end());
  }
  return {firstLayer, lastLayer, showerLen, coreShowerLen};
}

double HGCalHistoClusterProperties::convertRozToEta( HGCalClusterSAPtr& cluster ) const {
  // TODO : named constants for magic numbers
  double roz = double(cluster->wroz())/cluster->w();
  if ( roz < 1026.9376220703125 ) roz = 1026.9376220703125;
  else if ( roz > 5412.17138671875 ) roz = 5412.17138671875;
  roz -= 1026.9376220703125;
  roz *= 0.233510936;
  roz = int(round(roz));
  if ( roz > 1023 ) roz = 1023;
  return config_.rozToEtaLUT(roz);
}

double HGCalHistoClusterProperties::convertSigmaRozRozToSigmaEtaEta( HGCalClusterSAPtr& cluster ) const {
  // TODO : named constants for magic numbers
  // Sigma eta eta calculation
  double roz = cluster->wroz()/cluster->w();
  const double min_roz = 809.9324340820312;
  const double max_roz = 4996.79833984375;
  if ( roz < min_roz ) roz = min_roz;
  else if ( roz > max_roz ) roz = max_roz;
  roz -= min_roz;
  const double scale = 0.015286154113709927;
  roz *= scale;
  roz = int(round(roz));
  if ( roz > 63 ) roz = 63;

  const double sigma_roz_scale = 0.220451220870018;
  double sigmaRoz = sigma_coordinate(cluster->w(), cluster->wroz2(), cluster->wroz(), sigma_roz_scale);

  sigmaRoz = int(round(sigmaRoz));
  if ( sigmaRoz > 63 ) sigmaRoz = 63;
  unsigned int lutAddress = roz * 64 + sigmaRoz;
  if ( lutAddress >= 4096 ) lutAddress = 4095;
  return config_.sigmaRozToSigmaEtaLUT(lutAddress);
}
