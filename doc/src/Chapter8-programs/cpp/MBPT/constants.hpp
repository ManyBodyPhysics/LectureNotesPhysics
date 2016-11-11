#ifndef CONSTANTS_H
#define CONSTANTS_H

#define hbarc_MeVfm 197.3269788 // MeV fm
#define m_neutronc2 939.5654133 // MeV
#define m_protonc2 938.2720813 // MeV

const double neutron_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_neutronc2);
const double proton_prefac = hbarc_MeVfm * hbarc_MeVfm / (2.0 * m_protonc2);

#endif
