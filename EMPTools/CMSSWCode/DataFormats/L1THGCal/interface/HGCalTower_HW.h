#ifndef DataFormats_L1THGCal_HGCalTower_HW_h
#define DataFormats_L1THGCal_HGCalTower_HW_h

#include "DataFormats/L1TParticleFlow/interface/bit_encoding.h"

#include <array>
// #include <cstdint>

namespace l1thgcfirmware {

  typedef ap_uint<10> e_tower_t;
  typedef ap_uint<3> hoe_tower_t;
  typedef ap_uint<3> fb_tower_t;

  typedef ap_uint<16> TowerWord;
  
  namespace Scales {
    constexpr float ET_TOWER_L1_LSB = 0.25;

    inline float floatEt(e_tower_t et) { return et.to_float() * ET_TOWER_L1_LSB; }

    inline e_tower_t HGCaltoL1_towerEt(double et) {
      return ( ap_fixed<10,8>(et) << 2 );
    }
  }


  struct HGCalTower_HW {

    e_tower_t e;
    hoe_tower_t hoe;
    fb_tower_t featureBits;

    inline bool operator==(const HGCalTower_HW &other) const {
      return e == other.e && hoe == other.hoe && featureBits == other.featureBits;
    }

    inline bool operator!=(const HGCalTower_HW &other) const {
      return !(*this == other);
    }

    inline void clear() {
      e = 0;
      hoe = 0;
      featureBits = 0;
    }

    inline TowerWord pack() const {
      TowerWord ret;
      unsigned int start = 0;
      pack_into_bits(ret, start, e);
      pack_into_bits(ret, start, hoe);
      pack_into_bits(ret, start, featureBits);
      return ret;
    }

    inline static void unpackWord(const TowerWord &src, HGCalTower_HW& tower) {
      tower.initFromBits(src);
      return;
    }

    inline void initFromBits(const TowerWord &src ) {
      unsigned int start = 0;
      unpack_from_bits(src, start, e);
      unpack_from_bits(src, start, hoe);
      unpack_from_bits(src, start, featureBits);
    }

    inline static HGCalTower_HW unpack(const TowerWord &src) {
        HGCalTower_HW tower;
        unpackWord( src, tower );
        return tower;
    }

    inline void setEncodedHOE( double h, double e ) {
      double hoe_float = 999;
      if ( e > 0 ) {
        hoe_float = h / e;
      }

      if ( e == 0 ) hoe = 7;
      else if ( hoe_float < 0.4 ) hoe = 0;
      else if ( hoe_float < 0.5) hoe = 1;
      else if ( hoe_float < 0.55) hoe = 2;
      else if ( hoe_float < 0.6) hoe = 3;
      else if ( hoe_float < 0.65) hoe = 4;
      else if ( hoe_float < 0.7) hoe = 5;
      else hoe = 6;
    }

    inline void setFeatureBits( double e ) {
      ap_uint<1> fb0 = e > 255 ;
      ap_uint<1> fb1 = 0;
      ap_uint<1> fb2 = 0;
      featureBits = (fb2, fb1, fb0);
    }
  };

  inline void clear(HGCalTower_HW &c) { c.clear(); }

}  // namespace l1thgcfirmware

#endif