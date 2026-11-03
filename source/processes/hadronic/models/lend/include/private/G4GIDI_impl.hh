// include/private/G4GIDI_impl.hh
#ifndef G4GIDI_impl_hh
#define G4GIDI_impl_hh

#include "MCGIDI.hpp"   // brings in MCGIDI::, GIDI::, PoPI::

struct G4GIDI_target::Impl {
    MCGIDI::Protare*                                m_MCGIDI_protare;
    MCGIDI::DomainHash                              m_domainHash;
    MCGIDI::URR_protareInfos                        m_URR_protareInfos;
    MCGIDI::Probabilities::ProbabilityBase2d const* m_elasticAngular;

    std::string          m_target, m_fileName, m_evaluation;
    int                  m_targetZ, m_targetA, m_targetM;
    double               m_targetMass;
    std::vector<int>     m_elasticIndices, m_captureIndices,
                         m_fissionIndices, m_othersIndices;

    Impl( PoPI::Database const&, MCGIDI::DomainHash const&,
          GIDI::Protare const&, MCGIDI::Protare* );
    ~Impl();
};

struct G4GIDI::Impl {
    G4int                        m_projectileIP;
    std::string                  m_projectile;
    std::vector<GIDI::Map::Map*> m_maps;
    std::vector<G4GIDI_target*>  m_protares;
};

#endif
