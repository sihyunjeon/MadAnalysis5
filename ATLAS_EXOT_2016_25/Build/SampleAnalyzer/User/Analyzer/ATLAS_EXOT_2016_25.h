#ifndef analysis_ATLAS_EXOT_2016_25_h
#define analysis_ATLAS_EXOT_2016_25_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"

namespace MA5
{
class ATLAS_EXOT_2016_25 : public AnalyzerBase
{
  INIT_ANALYSIS(ATLAS_EXOT_2016_25,"ATLAS_EXOT_2016_25")

 public:
  virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual bool Execute(SampleFormat& sample, const EventFormat& event);

 private:
  double math_pi;
  double GetdPhi(double phi1, double phi2);

};
}

#endif
